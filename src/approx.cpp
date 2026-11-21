// Approximate DNA matching (k-differences / edit distance): exact algorithms and heuristics in one driver.
//
//   sellers    : Sellers (1980) DP, O(mn)                                  baseline / oracle
//   banded     : Ukkonen (1985) cutoff DP, O(kn)
//   myers      : Myers (1999) bit-vector (blocks), O(ceil(m/w) n)
//   myers_pruned: Myers with Ukkonen cutoff on blocks (only blocks that can hold a cell <= k are computed)
//   lv         : Landau-Vishkin (1989) diagonal algorithm, O(nk) expected (2-bit packed text, XOR + ctz extension, temporal blocking)
//   lv_stream  : packed extension, layers streamed over the whole text (intermediate variant)
//   lv_direct  : the same with direct character comparison (the earlier implementation, kept for comparison)
//   filter     : pigeonhole filter, k+1 EQUAL-length pieces  + Myers verification, cost-model fallback
//   filter_fa  : pigeonhole filter, k+1 pieces chosen by DP to minimise total seed occurrences
//   filter_s1/s2/s3, filter_fs2/fs3: generalised pigeonhole: k+s seeds (equal split / frequency-aware), s distinct seed hits within a diagonal band of width k
//   filter_auto: per query, the cost model picks equal split, frequency-aware split (fast selection) or a full Myers scan
//   filter_fa_full: filter_fa with the original (slow) seed-selection DP, kept for comparison
//   filter_fa_wg: filter_fa_w only for repeat-heavy reads (gate on equal-split seed hit counts), fast filter_fa otherwise
//   filter_fa_w: as filter_fa but the DP minimises hits + lambda * windows (windows made by each seed's own hits)
//   blast1     : BLAST-style one-hit seed (11-mer) + ungapped X-drop extension + Myers on window (heuristic)
//   blast2     : same with BLAST two-hit trigger (two non-overlapping seeds on one diagonal, window A=40)
//   myers_best : Myers scan reporting only the ends with the minimum distance (<= k); used to check Edlib
//
// Semantics: report every END index e (0-based) in text T such that some substring of T ending at e has
// edit distance <= k from the pattern P (empty substrings are not reported).
//
//   map        : exact ends from `filter`, then best alignment per read with CIGAR (lines "A <read> <dist> <start> <end> <cigar> <nends>")
// Usage: approx <algo> <text-file> <patterns-file> <k> [--dump] [--times]
//   --times prints one line per pattern: T <i> <ns> <hits> <verified_bp> <pred_filter_ns> <pred_scan_ns> <fallbacks>
// Env:   APPROX_NOFALLBACK=1  never fall back to a full Myers scan (measure the pure filter)
//        APPROX_CM="c_hit,c_win,c_bp,a_scan,c_sel" cost-model constants in ns (see scripts/calibrate.py)
#include <bits/stdc++.h>
#include <malloc.h>
#if defined(__AVX512F__) && defined(__AVX512CD__) && defined(__AVX512VL__)
#include <immintrin.h>
#define LV_SIMD 1
#endif
using namespace std;
typedef unsigned long long u64;
typedef vector<uint8_t> Seq;
typedef long long ll;

static inline uint8_t code(char c) {
    switch (c | 32) { case 'a': return 0; case 'c': return 1; case 'g': return 2; case 't': return 3; }
    return 4;  // N / anything else: matches nothing
}
static inline double now_ns() {
    return chrono::duration<double, nano>(chrono::steady_clock::now().time_since_epoch()).count();
}

// ---------------------------------------------------------------- statistics
struct Stats {
    ll hits = 0, windows = 0, verified_bp = 0, fallbacks = 0, filtered = 0;
    ll cells_banded = 0;         // DP cells actually computed by the banded algorithm
    ll ungapped = 0, hsps = 0;   // heuristic: extensions tried / HSPs above cutoff
    ll auto_s2 = 0, auto_s3 = 0;   // filter_auto picked the 2-hit / 3-hit generalised filter
    ll auto_skip = 0, auto_eq = 0, auto_fa = 0, auto_scan = 0;   // filter_auto decisions (skip = DP not run)
    double gen_ns = 0, ver_ns = 0, sel_ns = 0, pred_filter_ns = 0, pred_scan_ns = 0;
} S;
static bool g_nofallback = false;
static int SELMODE = 2;    // seed selection: 0 = original DP (all lengths, per-seed lookups), 1 = prefix tables + reused buffers (all lengths), 2 = tables + capped lengths (env APPROX_SEL)
static int SELCAP = 2;     // mode 2: seed lengths in [s0 - cap, min(11, s0 + cap)], s0 = min(floor(m/(k+1)), 11) (env APPROX_CAP)
static double AUTO_SKIP = 1.0;   // filter_auto: skip the DP when predicted equal-split cost <= AUTO_SKIP * estimated selection cost (1.0: the DP cannot save more than the whole equal-split cost) (env APPROX_AUTOSKIP)
static int g_debug = 0;
static double g_sel_nsop = 4.0;  // online estimate of seed-selection ns per DP unit (updated after every selection)
static int AUTO_S = 1;         // filter_auto may choose the s-hit filters (env APPROX_AUTOS=0 disables)
static double AUTO_S_TRIG = 2;  // ... only when the best fixed cost exceeds AUTO_S_TRIG * estimated selection cost (env APPROX_AUTOSTRIG)
static int WBIN = 1;       // filter_fa_w: 1 = binned window estimate (no sort), 0 = sort the hit diagonals (env APPROX_WBIN)
static int GATE_RATIO = 20;  // filter_fa_wg: repeat-heavy if max equal-split seed count >= GATE_RATIO * expected random count and >= GATE_MIN (env APPROX_GATE)
static int GATE_MIN = 100;
static int WCAP = 20000;   // filter_fa_w: seeds with more hits than this are charged one window per hit (env APPROX_WCAP)
static double CM_HIT = 40, CM_WIN = 150, CM_BP = 6, CM_SCAN = 6, CM_SEL = 15;  // ns; overridden by APPROX_CM

// ---------------------------------------------------------------- 1. Sellers
static void sellers(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) {
    int m = P.size();
    vector<int> C(m + 1);
    for (int i = 0; i <= m; i++) C[i] = i;
    for (size_t j = 1; j <= n; j++) {
        uint8_t t = T[j - 1];
        int diag = 0, prev = 0;  // D[0][j-1]=0, D[0][j]=0 (free start in text)
        for (int i = 1; i <= m; i++) {
            int tmp = C[i];
            int cost = (P[i - 1] == t && t < 4) ? 0 : 1;
            int v = min(diag + cost, min(C[i] + 1, prev + 1));
            diag = tmp; C[i] = v; prev = v;
        }
        if (C[m] <= k) out.push_back(j - 1);
    }
}

// ---------------------------------------------------------------- 2. Ukkonen cutoff
static void banded(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) {
    int m = P.size();
    vector<int> C(m + 1, k + 1);
    for (int i = 0; i <= min(m, k); i++) C[i] = i;
    int lact = min(m, k);  // last active row: largest i with C[i] <= k
    ll cells = 0;
    for (size_t j = 1; j <= n; j++) {
        uint8_t t = T[j - 1];
        int diag = 0, prev = 0;
        int top = min(lact + 1, m);
        cells += top;
        for (int i = 1; i <= top; i++) {
            int tmp = C[i];
            int cost = (P[i - 1] == t && t < 4) ? 0 : 1;
            int v = min(diag + cost, min(C[i] + 1, prev + 1));
            diag = tmp; C[i] = v; prev = v;
        }
        lact = top;
        while (lact > 0 && C[lact] > k) lact--;
        if (lact == m) out.push_back(j - 1);
    }
    S.cells_banded += cells;
}

// ---------------------------------------------------------------- 3. Myers bit-vector (blocks)
struct MyersPat {
    int m = 0, B = 0;
    u64 lastHigh = 0;
    vector<u64> peq;  // peq[c*B + b], c in 0..4
    explicit MyersPat(const Seq& P) {
        m = P.size(); B = (m + 63) / 64;
        peq.assign(5 * (size_t)B, 0);
        for (int i = 0; i < m; i++) if (P[i] < 4) peq[(size_t)P[i] * B + i / 64] |= 1ULL << (i % 64);
        lastHigh = 1ULL << ((m - 1) % 64);
    }
};

// Advance one 64-row block by one text char; returns horizontal delta out of the block.
// Pv/Mv: vertical deltas (+1 / -1) of the current column; Eq: match mask of the text char;
// hin: horizontal delta entering the block's top row; high: the bit whose row's horizontal delta is returned.
static inline int advance_block(u64& Pv, u64& Mv, u64 Eq, int hin, u64 high) {
    u64 Xv = Eq | Mv;
    if (hin < 0) Eq |= 1ULL;
    u64 Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
    u64 Ph = Mv | ~(Xh | Pv);
    u64 Mh = Pv & Xh;
    int hout = 0;
    if (Ph & high) hout = 1;
    if (Mh & high) hout = -1;
    Ph <<= 1; Mh <<= 1;
    if (hin < 0) Mh |= 1ULL; else if (hin > 0) Ph |= 1ULL;
    Pv = Mh | ~(Xv | Ph);
    Mv = Ph & Xv;
    return hout;
}

// Scan T[0..len); call emit(index, score) for each column; score = D[m][j].
template <class F>
static inline void myers_columns(const MyersPat& mp, const uint8_t* T, size_t len, F emit) {
    int B = mp.B;
    static thread_local vector<u64> PvV, MvV;
    PvV.assign(B, ~0ULL); MvV.assign(B, 0ULL);
    u64* Pv = PvV.data(); u64* Mv = MvV.data();
    const u64* peq = mp.peq.data();
    int score = mp.m;
    for (size_t j = 0; j < len; j++) {
        const u64* eq = peq + (size_t)T[j] * B;
        int hin = 0;
        for (int b = 0; b < B - 1; b++) hin = advance_block(Pv[b], Mv[b], eq[b], hin, 1ULL << 63);
        hin = advance_block(Pv[B - 1], Mv[B - 1], eq[B - 1], hin, mp.lastHigh);
        score += hin;
        emit(j, score);
    }
}

// Single-word fast path (m <= 64, one block): all state in registers, no vectors, no carry bookkeeping.  Same recurrence as advance_block with
// hin = 0.  Used by plain Myers and by the pruned scan (with one block there is nothing to prune).  Peq is built once per read (MyersPat).
static void myers_scan_1(const MyersPat& mp, const uint8_t* T, size_t len, int k, size_t base, vector<size_t>& out) {
    const u64* peq = mp.peq.data();
    u64 Pv = ~0ULL, Mv = 0; int score = mp.m; const int sh = mp.m - 1;
    for (size_t j = 0; j < len; j++) {
        u64 Eq = peq[T[j]];
        u64 Xv = Eq | Mv;
        u64 Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
        u64 Ph = Mv | ~(Xh | Pv);
        u64 Mh = Pv & Xh;
        score += (int)((Ph >> sh) & 1) - (int)((Mh >> sh) & 1);
        Ph <<= 1; Mh <<= 1;
        Pv = Mh | ~(Xv | Ph);
        Mv = Ph & Xv;
        if (score <= k) out.push_back(base + j);
    }
}
static int MYERS_FAST1 = 1;   // env APPROX_FAST1=0 disables the single-word fast path (for before/after timing)

static void myers_scan(const MyersPat& mp, const uint8_t* T, size_t len, int k, size_t base, vector<size_t>& out) {
    if (MYERS_FAST1 && mp.B == 1) { myers_scan_1(mp, T, len, k, base, out); return; }
    myers_columns(mp, T, len, [&](size_t j, int sc) { if (sc <= k) out.push_back(base + j); });
}
static void myers(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) {
    MyersPat mp(P);
    myers_scan(mp, T, n, k, 0, out);
}

// ---------------------------------------------------------------- 3c. Myers with block pruning (Ukkonen cutoff on blocks)
// Only blocks 0..y are computed.  S[b] = D[last row of block b][current column].  Invariant: every cell of the current column in
// blocks above y has value > k.  Block y+1 is activated when S[y] <= k (a cell <= k can only enter it through its first row);
// its state is initialised as "vertical deltas all +1" (stand-in values that are upper bounds of true values, which are > k, so they
// can never create a wrong value <= k).  Block y is dropped when S[y] - (rows(y) - 1) > k, because vertical deltas are >= -1 so every
// cell of the block is then > k.  Reports D[m][j] <= k only when the last block is active.  Expected active blocks: about ceil(k/64) + 1.
static ll pruned_block_steps = 0;
static void myers_pruned_scan(const MyersPat& mp, const uint8_t* T, size_t len, int k, size_t base, vector<size_t>& out) {
    if (MYERS_FAST1 && mp.B == 1) { myers_scan_1(mp, T, len, k, base, out); return; }   // one block: nothing to prune
    int B = mp.B, m = mp.m;
    static thread_local vector<u64> PvV, MvV; static thread_local vector<int> S;
    PvV.assign(B, ~0ULL); MvV.assign(B, 0ULL); S.assign(B, 0);
    auto rows = [&](int b) { return b == B - 1 ? m - 64 * (B - 1) : 64; };
    int y = 0;
    if (k > 0) y = max(0, (min(m, k) - 1) / 64);              // blocks that contain a row i <= k at column 0 (D[i][0] = i)
    for (int b = 0; b <= y; b++) S[b] = min(m, 64 * (b + 1));
    u64* Pv = PvV.data(); u64* Mv = MvV.data();
    const u64* peq = mp.peq.data();
    for (size_t j = 0; j < len; j++) {
        if (y < B - 1 && S[y] <= k) {                           // activate the next block
            y++; Pv[y] = ~0ULL; Mv[y] = 0ULL; S[y] = S[y - 1] + rows(y);
        }
        const u64* eq = peq + (size_t)T[j] * B;
        int hin = 0;
        for (int b = 0; b <= y; b++) {
            u64 high = (b == B - 1) ? mp.lastHigh : (1ULL << 63);
            hin = advance_block(Pv[b], Mv[b], eq[b], hin, high);
            S[b] += hin;
        }
        pruned_block_steps += y + 1;
        while (y > 0 && S[y] >= k + rows(y)) y--;                // block y has no cell <= k: drop it
        if (y == B - 1 && S[y] <= k) out.push_back(base + j);
    }
}
static void myers_pruned(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) {
    MyersPat mp(P);
    myers_pruned_scan(mp, T, n, k, 0, out);
}

// Only the ends achieving the minimum distance, if that minimum is <= k (this is what Edlib HW reports).
static void myers_best(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) {
    MyersPat mp(P);
    int best = INT_MAX; vector<size_t> ends;
    myers_columns(mp, T, n, [&](size_t j, int sc) {
        if (sc < best) { best = sc; ends.clear(); }
        if (sc == best) ends.push_back(j);
    });
    if (best <= k) out.insert(out.end(), ends.begin(), ends.end());
}


// ---------------------------------------------------------------- 3b. Landau-Vishkin diagonal algorithm, O(nk) expected
// L[d][e] = furthest row i (pattern characters consumed) reached on diagonal d = j - i (j = text column) with at most e
// errors.  Row 0 is free on every diagonal d >= 0 (a match may start anywhere), D[i][0] = i makes d = -i reachable
// from d = 0 by deletions.  Recurrence (layers e = 0..k, two arrays over diagonals):
//     L[d][e] = slide( max( L[d][e-1] + 1,  L[d-1][e-1],  L[d+1][e-1] + 1 ) )          (substitution, insertion, deletion)
// where slide extends along the diagonal while P[i] == T[i + d].  Diagonal d reaches the last row m with e errors
// iff D[m][m + d] <= e, i.e. end index m + d - 1.  Extension uses direct character comparison: expected O(1) per cell on
// random text, so O(nk) expected; worst case O(nm) (long repeats).  An LCE structure (suffix array + LCP + RMQ) would
// give worst-case O(nk) but is NOT implemented.
static ll lv_cells = 0;
static void lv_direct(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) {   // original: direct character comparison, kept for comparison
    int m = P.size();
    ll lo = -k, hi = (ll)n - m + k;                       // diagonals that can matter
    if (hi < lo) return;
    size_t W = (size_t)(hi - lo + 1);
    const int NEG = -(1 << 29);
    vector<int> cur(W + 2, NEG), prv(W + 2, NEG);         // index = d - lo + 1 (guard cells at both ends)
    vector<char> done(W, 0);                              // diagonal already reported
    auto slide = [&](ll d, int i) {
        ll j = i + d;
        while (i < m && j < (ll)n && j >= 0 && P[i] == T[j] && T[j] < 4) { i++; j++; }
        return i;
    };
    for (int e = 0; e <= k; e++) {
        for (size_t x = 0; x < W; x++) {
            ll d = lo + (ll)x;
            int best = NEG;
            if (e == 0) { if (d >= 0) best = 0; }
            else {
                int a = prv[x + 1] <= NEG / 2 ? NEG : prv[x + 1] + 1;      // substitution (same diagonal)
                int b = prv[x];                                             // insertion: from diagonal d-1, row unchanged
                int c = prv[x + 2] <= NEG / 2 ? NEG : prv[x + 2] + 1;      // deletion: from diagonal d+1, row + 1
                best = max(a, max(b, c));
            }
            if (best <= NEG / 2) { cur[x + 1] = NEG; continue; }
            if (best > m) best = m;
            ll j = best + d;
            if (j > (ll)n) { cur[x + 1] = NEG; continue; }              // would run past the text end
            if (j < 0) { cur[x + 1] = NEG; continue; }
            best = slide(d, best);
            cur[x + 1] = best; lv_cells++;
            if (best >= m && !done[x]) {
                ll jend = (ll)m + d;                                        // exclusive end column
                if (jend >= 1) { out.push_back((size_t)(jend - 1)); done[x] = 1; }
            }
        }
        swap(cur, prv);
    }
}


// ---------------------------------------------------------------- 3b'. Landau-Vishkin with 2-bit packed text and pattern, XOR + ctz extension
// Same recurrence as lv_direct, restructured for speed:
//  * text packed 2 bits per base (32 bases per u64, built once per text); every pattern window P[i..i+32) precomputed once per read;
//  * extension = XOR of the two 32-base windows and count-trailing-zeros (a mismatch is the lowest set 2-bit group), no per-character loop
//    and no data-dependent branch in the common case;
//  * the layer step runs in tiles of 2048 diagonals: pass 1 (auto-vectorised) computes the three-way maximum into a small buffer,
//    pass 2 does the branch-free extension;
//  * N (code 4) is handled exactly: packed N is code 0, so the extension is truncated at the first N of the text (sorted position list)
//    and of the pattern (next-N table); both fix-ups are skipped when there is no N.
static vector<u64> g_T2; static vector<uint32_t> g_Npos; static const uint8_t* g_T2src = nullptr;
static inline u64 get32(const u64* a, size_t pos) {   // 32 bases starting at base index pos
    size_t w = pos >> 5; unsigned sft = (pos & 31) * 2;
    return (a[w] >> sft) | ((a[w + 1] << 1) << (63 - sft));
}
static void build_T2(const uint8_t* T, size_t n) {
    if (g_T2src == T) return;
    g_T2.assign(n / 32 + 3, 0); g_Npos.clear();
    for (size_t i = 0; i < n; i++) { u64 c = T[i]; if (c > 3) { g_Npos.push_back((uint32_t)i); c = 0; } g_T2[i >> 5] |= c << ((i & 31) * 2); }
    g_T2src = T;
}
static void lv_stream(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) {   // packed extension, layers streamed over the whole text
    int m = P.size();
    ll lo = -k, hi = (ll)n - m + k;
    if (hi < lo) return;
    size_t W = (size_t)(hi - lo + 1);
    build_T2(T, n);
    const u64* T2 = g_T2.data(); const bool hasN = !g_Npos.empty();
    vector<u64> P2((m >> 5) + 3, 0); bool pHasN = false;
    for (int i = 0; i < m; i++) { u64 c = P[i]; if (c > 3) { pHasN = true; c = 0; } P2[i >> 5] |= c << ((i & 31) * 2); }
    vector<u64> PW(m + 1);
    for (int i = 0; i <= m; i++) PW[i] = get32(P2.data(), i);
    vector<int> pnext;
    if (pHasN) { pnext.assign(m + 1, m); for (int i = m - 1; i >= 0; i--) pnext[i] = P[i] > 3 ? i : pnext[i + 1]; }
    const int NEG = -(1 << 29);
    vector<int> cur(W + 2, NEG), prv(W + 2, NEG);
    const size_t TILE = 2048;
    int cand[TILE];
    for (int e = 0; e <= k; e++) {
        for (size_t x0 = 0; x0 < W; x0 += TILE) {
            size_t tn = min(TILE, W - x0);
            if (e == 0) { for (size_t t = 0; t < tn; t++) cand[t] = (lo + (ll)(x0 + t) >= 0) ? 0 : NEG; }
            else {
                const int* pv = prv.data() + x0;
                for (size_t t = 0; t < tn; t++) { int a = pv[t + 1] + 1, b = pv[t], c = pv[t + 2] + 1; int v = max(a, max(b, c)); cand[t] = min(v, m); }
            }
            for (size_t t = 0; t < tn; t++) {
                ll d = lo + (ll)(x0 + t); int i = cand[t]; ll j = i + d;
                bool valid = (i >= 0) & (j <= (ll)n) & (j >= 0);
                int ii = valid ? i : 0; size_t jj = valid ? (size_t)j : 0;
                u64 xr = PW[ii] ^ get32(T2, jj);
                int ext = xr ? (__builtin_ctzll(xr) >> 1) : 32;
                if (ext == 32) {                                      // all 32 bases equal: continue block by block (rare on random text)
                    while (true) {
                        size_t ip = ii + ext, jp = jj + ext;
                        if (ip >= (size_t)m || jp >= n) break;
                        u64 x2 = get32(P2.data(), ip) ^ get32(T2, jp);
                        if (x2) { ext += __builtin_ctzll(x2) >> 1; break; }
                        ext += 32;
                    }
                }
                ext = min(ext, min(m - ii, (int)min<size_t>(n - jj, (size_t)m)));
                if (hasN) { auto it = lower_bound(g_Npos.begin(), g_Npos.end(), (uint32_t)jj); if (it != g_Npos.end() && (ll)*it - (ll)jj < ext) ext = (int)((ll)*it - (ll)jj); }
                if (pHasN) ext = min(ext, pnext[ii] - ii);
                int row = valid ? ii + ext : NEG;
                cur[x0 + t + 1] = row;
                if (row >= m) {                                       // reached the last row: report once (first layer that gets there)
                    bool already = (e > 0) && prv[x0 + t + 1] >= m;
                    ll jend = (ll)m + d;
                    if (!already && jend >= 1) out.push_back((size_t)(jend - 1));
                }
            }
        }
        lv_cells += (ll)W;
        swap(cur, prv);
        cur[0] = NEG; cur[W + 1] = NEG;
    }
}


// Final version: the same packed extension, but with TEMPORAL BLOCKING.  A chunk of C = 2048 diagonals is carried through all k + 1 layers
// in two small buffers that stay in L1; at layer e the chunk needs the diagonals [x0 - (k - e), x1 + (k - e)) because layer e reads the
// neighbours of layer e - 1 (a trapezoid, recomputing about 2k cells per chunk, 3% at k = 30).  Only cells inside [x0, x1) are reported.
// This removes the 40 MB per layer that lv_stream moves through memory.
static int LV_USE_SIMD = 1;   // env APPROX_LVSIMD=0 forces the scalar cell loop
static void lv(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) {
    int m = P.size();
    ll lo = -k, hi = (ll)n - m + k;
    if (hi < lo) return;
    ll W = hi - lo + 1;
    build_T2(T, n);
    const u64* T2 = g_T2.data(); const bool hasN = !g_Npos.empty();
    vector<u64> P2((m >> 5) + 3, 0); bool pHasN = false;
    for (int i = 0; i < m; i++) { u64 c = P[i]; if (c > 3) { pHasN = true; c = 0; } P2[i >> 5] |= c << ((i & 31) * 2); }
    vector<u64> PW(m + 1);
    for (int i = 0; i <= m; i++) PW[i] = get32(P2.data(), i);
    vector<int> pnext;
    if (pHasN) { pnext.assign(m + 1, m); for (int i = m - 1; i >= 0; i--) pnext[i] = P[i] > 3 ? i : pnext[i + 1]; }
    const int NEG = -(1 << 29);
    const ll C = 2048;
    const size_t LS = (size_t)(C + 2 * k + 4);
    static thread_local vector<int> bufA, bufB, candv;
    bufA.assign(LS, NEG); bufB.assign(LS, NEG); candv.assign(LS, 0);
    for (ll x0 = 0; x0 < W; x0 += C) {
        ll x1 = min(W, x0 + C), b0 = x0 - k;              // local index of diagonal x is x - b0 + 1
        fill(bufA.begin(), bufA.end(), NEG); fill(bufB.begin(), bufB.end(), NEG);
        int* prv = bufA.data(); int* cur = bufB.data(); int* cand = candv.data();
        for (int e = 0; e <= k; e++) {
            ll xa = max(0LL, x0 - (k - e)), xb = min(W, x1 + (k - e));
            size_t tn = (size_t)(xb - xa);
            if (e == 0) { for (size_t t = 0; t < tn; t++) cand[t] = (lo + xa + (ll)t >= 0) ? 0 : NEG; }
            else {
                const int* pv = prv + (xa - b0);            // pv[t + 1] = diagonal xa + t of layer e - 1
                for (size_t t = 0; t < tn; t++) { int a = pv[t + 1] + 1, b = pv[t], c = pv[t + 2] + 1; int v = max(a, max(b, c)); cand[t] = min(v, m); }
            }
            int* cu = cur + (xa - b0 + 1);
            const int* pr = prv + (xa - b0 + 1);
            auto scalar_cell = [&](size_t t) {
                ll d = lo + xa + (ll)t; int i = cand[t]; ll j = i + d;
                bool valid = (i >= 0) & (j <= (ll)n) & (j >= 0);
                int ii = valid ? i : 0; size_t jj = valid ? (size_t)j : 0;
                u64 xr = PW[ii] ^ get32(T2, jj);
                int ext = xr ? (__builtin_ctzll(xr) >> 1) : 32;
                if (ext == 32) {
                    while (true) {
                        size_t ip = ii + ext, jp = jj + ext;
                        if (ip >= (size_t)m || jp >= n) break;
                        u64 x2 = get32(P2.data(), ip) ^ get32(T2, jp);
                        if (x2) { ext += __builtin_ctzll(x2) >> 1; break; }
                        ext += 32;
                    }
                }
                ext = min(ext, min(m - ii, (int)min<size_t>(n - jj, (size_t)m)));
                if (hasN) { auto it = lower_bound(g_Npos.begin(), g_Npos.end(), (uint32_t)jj); if (it != g_Npos.end() && (ll)*it - (ll)jj < ext) ext = (int)((ll)*it - (ll)jj); }
                if (pHasN) ext = min(ext, pnext[ii] - ii);
                int row = valid ? ii + ext : NEG;
                cu[t] = row;
                if (row >= m) {
                    ll x = xa + (ll)t;
                    if (x >= x0 && x < x1) {
                        bool already = (e > 0) && pr[t] >= m;
                        ll jend = (ll)m + d;
                        if (!already && jend >= 1) out.push_back((size_t)(jend - 1));
                    }
                }
            };
            size_t t = 0;
#ifdef LV_SIMD
            if (LV_USE_SIMD && !hasN && !pHasN) {
                const __m512i zero = _mm512_setzero_si512(), nvec = _mm512_set1_epi64((ll)n), mvec = _mm512_set1_epi64(m);
                const __m512i three = _mm512_set1_epi64(3), mask58 = _mm512_set1_epi64((1LL << 58) - 1), c63 = _mm512_set1_epi64(63);
                const __m512i lanes = _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0), negv = _mm512_set1_epi64(NEG), c29 = _mm512_set1_epi64(29);
                const unsigned char* T2b = (const unsigned char*)T2;
                for (; t + 8 <= tn; t += 8) {
                    __m512i vd = _mm512_add_epi64(_mm512_set1_epi64(lo + xa + (ll)t), lanes);
                    __m512i vi = _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i*)(cand + t)));
                    __m512i vj = _mm512_add_epi64(vi, vd);
                    __mmask8 valid = _mm512_cmpge_epi64_mask(vi, zero) & _mm512_cmple_epi64_mask(vj, nvec) & _mm512_cmpge_epi64_mask(vj, zero);
                    vi = _mm512_maskz_mov_epi64(valid, vi); vj = _mm512_maskz_mov_epi64(valid, vj);
                    __m512i pw = _mm512_i64gather_epi64(vi, (const long long*)PW.data(), 8);
                    __m512i boff = _mm512_srli_epi64(vj, 2);
                    __m512i sh = _mm512_slli_epi64(_mm512_and_si512(vj, three), 1);
                    __m512i tw = _mm512_srlv_epi64(_mm512_i64gather_epi64(boff, (const long long*)T2b, 1), sh);
                    __m512i x = _mm512_and_si512(_mm512_xor_si512(pw, tw), mask58);
                    __m512i low = _mm512_and_si512(x, _mm512_sub_epi64(zero, x));
                    __m512i ext = _mm512_srli_epi64(_mm512_sub_epi64(c63, _mm512_lzcnt_epi64(low)), 1);
                    __mmask8 allEq = _mm512_cmpeq_epi64_mask(x, zero);
                    ext = _mm512_mask_mov_epi64(ext, allEq, c29);
                    ext = _mm512_min_epi64(ext, _mm512_min_epi64(_mm512_sub_epi64(mvec, vi), _mm512_sub_epi64(nvec, vj)));
                    alignas(64) ll ea[8], ia[8], ja[8];
                    __mmask8 cont = allEq & _mm512_cmpeq_epi64_mask(ext, c29);   // all 29 bases equal and not limited by the ends: continue per lane
                    if (cont) {
                        _mm512_store_si512(ea, ext); _mm512_store_si512(ia, vi); _mm512_store_si512(ja, vj);
                        for (int l = 0; l < 8; l++) if (cont >> l & 1) {
                            size_t ii = (size_t)ia[l], jj = (size_t)ja[l]; int ex = 29;
                            while (true) {
                                size_t ip = ii + ex, jp = jj + ex;
                                if (ip >= (size_t)m || jp >= n) break;
                                u64 x2 = get32(P2.data(), ip) ^ get32(T2, jp);
                                if (x2) { ex += __builtin_ctzll(x2) >> 1; break; }
                                ex += 32;
                            }
                            ex = min(ex, (int)min((ll)(m - (ll)ii), (ll)(n - jj)));
                            ea[l] = ex;
                        }
                        ext = _mm512_load_si512(ea);
                    }
                    __m512i row = _mm512_mask_blend_epi64(valid, negv, _mm512_add_epi64(vi, ext));
                    _mm256_storeu_si256((__m256i*)(cu + t), _mm512_cvtepi64_epi32(row));
                    __mmask8 rep = valid & _mm512_cmpge_epi64_mask(row, mvec);
                    if (rep) {
                        for (int l = 0; l < 8; l++) if (rep >> l & 1) {
                            ll x = xa + (ll)t + l;
                            if (x >= x0 && x < x1) {
                                bool already = (e > 0) && pr[t + l] >= m;
                                ll jend = (ll)m + lo + xa + (ll)t + l;
                                if (!already && jend >= 1) out.push_back((size_t)(jend - 1));
                            }
                        }
                    }
                }
            }
#endif
            for (; t < tn; t++) scalar_cell(t);
            lv_cells += (ll)tn;
            swap(cur, prv);
        }
    }
}

// Minimum edit distance of P to any substring of T (reported as a single "value"); verified against Sellers by scripts/real_divergence.py.
static void mindist(const uint8_t* T, size_t n, const Seq& P, int, vector<size_t>& out) {
    MyersPat mp(P);
    int best = INT_MAX;
    myers_columns(mp, T, n, [&](size_t, int sc) { if (sc < best) best = sc; });
    out.push_back((size_t)best);
}

// ---------------------------------------------------------------- seed index (direct-address 11-mer table)
static const int SMAX = 11;
struct Index {
    int s = 0;
    vector<uint32_t> start, pos;  // CSR: occurrences of 11-mer c are pos[start[c] .. start[c+1])
    vector<uint32_t> bad;         // text positions p whose s-mer is NOT indexed (tail, or contains N)
    vector<uint32_t> cntL[11];    // cntL[L][c] = occurrences of the L-mer c (6 <= L <= 10): precomputed prefix counts for seed selection
    void build(const uint8_t* T, size_t n, int s_) {
        s = s_;
        size_t SZ = (size_t)1 << (2 * s);
        start.assign(SZ + 1, 0);
        vector<char> indexed(n, 0);
        auto scan = [&](auto fn) {
            u64 h = 0, mask = SZ - 1; int valid = 0;
            for (size_t i = 0; i < n; i++) {
                if (T[i] > 3) { valid = 0; h = 0; continue; }
                h = ((h << 2) | T[i]) & mask;
                if (++valid >= s) fn(h, i + 1 - s);
            }
        };
        scan([&](u64 h, size_t p) { start[h + 1]++; indexed[p] = 1; });
        for (size_t i = 0; i < SZ; i++) start[i + 1] += start[i];
        pos.resize(start[SZ]);
        vector<uint32_t> fill(start.begin(), start.end() - 1);
        scan([&](u64 h, size_t p) { pos[fill[h]++] = p; });
        for (size_t p = 0; p < n; p++) if (!indexed[p]) bad.push_back(p);
        for (int L = 6; L <= 10 && L < s; L++) {                       // counts of every L-mer from the CSR offsets (4^L entries each)
            int sh = 2 * (s - L); size_t cnt_n = (size_t)1 << (2 * L);
            cntL[L].resize(cnt_n);
            for (size_t c = 0; c < cnt_n; c++) cntL[L][c] = start[(c + 1) << sh] - start[c << sh];
        }
    }
};
static const uint8_t* gT; static size_t gN;
static Index gIdx; static bool gIdxBuilt = false;
static const Index& get_index() {
    if (!gIdxBuilt) { gIdx.build(gT, gN, SMAX); gIdxBuilt = true; }
    return gIdx;
}

// ---------------------------------------------------------------- shared window machinery
typedef pair<ll, ll> Win;
static void merge_windows(vector<Win>& W) {
    sort(W.begin(), W.end());
    size_t o = 0;
    for (size_t i = 0; i < W.size(); i++) {
        if (o > 0 && W[i].first <= W[o - 1].second) W[o - 1].second = max(W[o - 1].second, W[i].second);
        else W[o++] = W[i];
    }
    W.resize(o);
}
static void verify_merged(const vector<Win>& M, const MyersPat& mp, const uint8_t* T, int k, vector<size_t>& out) {
    double t0 = now_ns();
    ll total = 0;
    for (auto& w : M) { total += w.second - w.first; myers_scan(mp, T + w.first, w.second - w.first, k, w.first, out); }
    S.windows += M.size(); S.verified_bp += total;
    S.ver_ns += now_ns() - t0;
}
static void scan_fallback(const MyersPat& mp, const uint8_t* T, size_t n, int k, vector<size_t>& out) {
    S.fallbacks++;
    myers_scan(mp, T, n, k, 0, out);
}

// ---------------------------------------------------------------- 4. Pigeonhole filter (equal split or frequency-aware split)
struct Seed { int a, len; };  // pattern offset, seed length (<= SMAX); seeds are pairwise disjoint

// number of index entries (11-mer occurrences) whose first `len` characters equal P[a..a+len)
static bool seed_range(const Index& ix, const Seq& P, const Seed& sd, uint32_t& lo, uint32_t& hi, u64& c) {
    c = 0;
    for (int i = 0; i < sd.len; i++) { if (P[sd.a + i] > 3) return false; c = (c << 2) | P[sd.a + i]; }
    u64 l = c << (2 * (ix.s - sd.len)), h = (c + 1) << (2 * (ix.s - sd.len));
    lo = ix.start[l]; hi = ix.start[h];
    return true;
}

// Frequency-aware split: choose k+1 disjoint seeds (length 1..SMAX) minimising total occurrences by DP.
// f[j][i] = min total count using j seeds inside P[0..i).  O(k * m * SMAX) after an O(m * SMAX) count table.
// wmode = true: objective is hits + lambda * windows, with lambda = (c_win + c_bp B (m + 2k)) / c_hit from the fitted cost model,
// and windows(seed) = number of merged windows made by that seed's own hits (clusters of hit diagonals no more than m + 2k
// apart).  Seeds shorter than 9 characters, or with more than 20000 hits, are charged one window per hit.
static bool select_fa(const Index& ix, const Seq& P, int k, vector<Seed>& seeds, bool wmode = false) {
    int m = P.size(), J = k + 1;
    if (m < J) return false;
    int B = (m + 63) / 64; ll Lw = m + 2LL * k;
    double lam = wmode ? (CM_WIN + CM_BP * B * (double)Lw) / CM_HIT : 0.0;
    vector<double> cnt((size_t)m * (SMAX + 1), 0.0);
    vector<ll> diag;
    for (int a = 0; a < m; a++)
        for (int L = 1; L <= SMAX && a + L <= m; L++) {
            uint32_t lo, hi; u64 c;
            Seed sd{a, L};
            double v = 0;
            if (seed_range(ix, P, sd, lo, hi, c)) {
                double hits = (double)(hi - lo), w = hits;
                if (wmode && L >= 9 && hits > 1 && hits <= WCAP) {
                    if (WBIN) {          // binned estimate: distinct bins of width Lw touched by this seed's hit diagonals (O(hits), no sort)
                        static thread_local vector<u64> wbits;
                        size_t nb = gN / Lw + 2;
                        if (wbits.size() < (nb + 63) / 64 + 1) wbits.assign((nb + 63) / 64 + 1, 0);
                        const u64 winv = ((1ULL << 32) + Lw - 1) / Lw;
                        ll occ = 0;
                        for (uint32_t x = lo; x < hi; x++) {
                            ll d = (ll)ix.pos[x] - a; if (d < 0) d = 0;
                            size_t b = (size_t)(((u64)d * winv) >> 32); u64 bit = 1ULL << (b & 63);
                            if (!(wbits[b >> 6] & bit)) { wbits[b >> 6] |= bit; occ++; }
                        }
                        for (uint32_t x = lo; x < hi; x++) {
                            ll d = (ll)ix.pos[x] - a; if (d < 0) d = 0;
                            size_t b = (size_t)(((u64)d * winv) >> 32); wbits[b >> 6] = 0;
                        }
                        w = (double)occ;
                    } else {
                        diag.clear();
                        for (uint32_t x = lo; x < hi; x++) diag.push_back((ll)ix.pos[x] - a);
                        sort(diag.begin(), diag.end());
                        w = 1;
                        for (size_t x = 1; x < diag.size(); x++) if (diag[x] - diag[x - 1] > Lw) w++;
                    }
                }
                v = hits + lam * w;
            }
            cnt[(size_t)a * (SMAX + 1) + L] = v;
        }
    const double INF = 4e18;
    vector<vector<double>> f(J + 1, vector<double>(m + 1, INF));
    vector<vector<int8_t>> bp(J + 1, vector<int8_t>(m + 1, 0));  // 0 = skip char i-1, else seed length
    for (int i = 0; i <= m; i++) f[0][i] = 0;
    for (int j = 1; j <= J; j++)
        for (int i = j; i <= m; i++) {
            double best = f[j][i - 1]; int8_t ch = 0;
            for (int L = 1; L <= SMAX && L <= i; L++) {
                double prev = f[j - 1][i - L];
                if (prev >= INF) continue;
                double v = prev + cnt[(size_t)(i - L) * (SMAX + 1) + L];
                if (v < best) { best = v; ch = L; }
            }
            f[j][i] = best; bp[j][i] = ch;
        }
    if (f[J][m] >= INF) return false;
    seeds.clear();
    int i = m, j = J;
    while (j > 0) {
        if (bp[j][i] == 0) i--;
        else { int L = bp[j][i]; seeds.push_back({i - L, L}); i -= L; j--; }
    }
    reverse(seeds.begin(), seeds.end());
    return true;
}

// ---- fast seed selection: precomputed prefix counts, capped piece lengths, reused buffers (hits objective)
static inline uint32_t seed_count_fast(const Index& ix, int L, u64 c) {
    if (L == 11) return ix.start[c + 1] - ix.start[c];
    if (L >= 6) return ix.cntL[L][c];
    int sh = 2 * (ix.s - L);
    return ix.start[(c + 1) << sh] - ix.start[c << sh];
}
static void sel_range(int m, int k, int cap, int& Lmin, int& Lmax) {
    int s0 = min(m / (k + 1), SMAX);
    Lmin = cap >= SMAX ? 1 : max(1, s0 - cap);
    Lmax = cap >= SMAX ? SMAX : min(SMAX, s0 + cap);
}
static double sel_ops(int m, int k, int cap) {   // DP units: count table m*nL lookups + (k+1)*m*nL transitions
    int Lmin, Lmax; sel_range(m, k, cap, Lmin, Lmax);
    double nL = Lmax - Lmin + 1;
    return (double)m * nL + (double)(k + 1) * m * nL;
}
static bool select_fa_fast(const Index& ix, const Seq& P, int k, vector<Seed>& seeds, int cap) {
    int m = P.size(), J = k + 1;
    if (m < J) return false;
    int Lmin, Lmax; sel_range(m, k, cap, Lmin, Lmax);
    int nL = Lmax - Lmin + 1;
    static thread_local vector<uint32_t> cnt; static thread_local vector<int64_t> f; static thread_local vector<int8_t> bp;
    static thread_local vector<int> nPre;
    cnt.assign((size_t)m * nL, 0);
    nPre.assign(m + 1, 0);
    for (int i = 0; i < m; i++) nPre[i + 1] = nPre[i] + (P[i] > 3);
    for (int a = 0; a + Lmin <= m; a++) {
        u64 c = 0;
        for (int t = 0; t < Lmin; t++) c = (c << 2) | (P[a + t] & 3);
        int top = min(Lmax, m - a);
        for (int L = Lmin; L <= top; L++) {
            if (L > Lmin) c = (c << 2) | (P[a + L - 1] & 3);
            if (nPre[a + L] - nPre[a] == 0) cnt[(size_t)a * nL + (L - Lmin)] = seed_count_fast(ix, L, c);   // a seed containing N has no occurrence: cost 0
        }
    }
    const int64_t INF = (int64_t)4e18;
    size_t W1 = m + 1;
    f.assign((size_t)(J + 1) * W1, INF); bp.assign((size_t)(J + 1) * W1, 0);
    for (int i = 0; i <= m; i++) f[i] = 0;
    for (int j = 1; j <= J; j++) {
        int64_t* fj = &f[(size_t)j * W1]; const int64_t* fp = &f[(size_t)(j - 1) * W1]; int8_t* bj = &bp[(size_t)j * W1];
        for (int i = j * Lmin; i <= m; i++) {
            int64_t best = fj[i - 1]; int8_t ch = 0;
            for (int L = Lmin; L <= Lmax && L <= i; L++) {
                int64_t prev = fp[i - L];
                if (prev >= INF) continue;
                int64_t v = prev + cnt[(size_t)(i - L) * nL + (L - Lmin)];
                if (v < best) { best = v; ch = (int8_t)L; }
            }
            fj[i] = best; bj[i] = ch;
        }
    }
    if (f[(size_t)J * W1 + m] >= INF) return false;
    seeds.clear();
    int i = m, j = J;
    while (j > 0) {
        int8_t ch = bp[(size_t)j * W1 + i];
        if (ch == 0) i--;
        else { seeds.push_back({i - ch, (int)ch}); i -= ch; j--; }
    }
    reverse(seeds.begin(), seeds.end());
    return true;
}

// ---- seed set: index ranges (plus brute-force hits at un-indexed positions) and the hit total H
struct SeedSet { vector<Seed> seeds; vector<uint32_t> lo, hi; vector<char> ok; vector<vector<uint32_t>> bad; ll H = 0; bool bad_done = false; };
static void collect_ranges(const Index& ix, const Seq& P, SeedSet& ss) {   // index ranges only (cheap: 2 lookups per seed)
    size_t ns = ss.seeds.size();
    ss.lo.assign(ns, 0); ss.hi.assign(ns, 0); ss.ok.assign(ns, 0); ss.H = 0; ss.bad_done = false;
    for (size_t i = 0; i < ns; i++) {
        u64 c; ss.ok[i] = seed_range(ix, P, ss.seeds[i], ss.lo[i], ss.hi[i], c);
        if (ss.ok[i]) ss.H += ss.hi[i] - ss.lo[i];
    }
}
static void add_bad(const Index& ix, const uint8_t* T, size_t n, const Seq& P, SeedSet& ss) {   // hits at the few un-indexed text positions (tail, N)
    if (ss.bad_done) return;
    ss.bad_done = true; ss.bad.assign(ss.seeds.size(), vector<uint32_t>());
    if (ix.bad.empty()) return;
    for (size_t i = 0; i < ss.seeds.size(); i++) {
        if (!ss.ok[i]) continue;
        for (uint32_t p : ix.bad) {
            if (p + ss.seeds[i].len > n) continue;
            bool eq = true;
            for (int t = 0; t < ss.seeds[i].len; t++) if (T[p + t] != P[ss.seeds[i].a + t] || T[p + t] > 3) { eq = false; break; }
            if (eq) ss.bad[i].push_back(p);
        }
        ss.H += ss.bad[i].size();
    }
}
static void collect(const Index& ix, const uint8_t* T, size_t n, const Seq& P, SeedSet& ss) { collect_ranges(ix, P, ss); add_bad(ix, T, n, P, ss); }
static double model_cost(ll H, int m, int k, int B, size_t n) {   // cost model, stage 1: expected merged windows / bases under uniform hit placement
    double Lw = m + 2.0 * k, lam = (double)H * Lw / (double)n;
    return CM_HIT * H + CM_WIN * H * exp(-lam) + CM_BP * B * n * (1.0 - exp(-lam));
}
// cost model with MEASURED occupancy: hits are dropped into bins of width m + 2k by diagonal; occupied bins estimate the verified length
// and runs of occupied bins the merged windows.  For uniformly placed hits this equals the Poisson estimate of model_cost; for clustered
// hits (repeats, low-complexity) it sees that they merge into few windows.  O(H).
static int OCC_SMALL = 64;   // below this many hits every hit is charged its own window (no bins): O(1) (env APPROX_OCCSMALL)
static double occupancy_cost(const Index& ix, const SeedSet& ss, int m, int k, int B, size_t n) {
    ll Lw = m + 2LL * k; size_t nb = n / Lw + 2;
    if (ss.H <= OCC_SMALL) return ss.H * (CM_HIT + CM_WIN + CM_BP * B * (double)Lw);
    static thread_local vector<u64> bits;
    bits.assign((nb + 63) / 64 + 1, 0);
    const u64 inv = ((1ULL << 32) + Lw - 1) / Lw;      // bin = d / Lw via multiply (exact enough for an estimate)
    for (size_t i = 0; i < ss.seeds.size(); i++) {
        if (!ss.ok[i]) continue;
        ll a = ss.seeds[i].a;
        for (uint32_t x = ss.lo[i]; x < ss.hi[i]; x++) {
            ll d = (ll)ix.pos[x] - a; if (d < 0) d = 0;
            size_t b = (size_t)(((u64)d * inv) >> 32); bits[b >> 6] |= 1ULL << (b & 63);
        }
    }
    ll occ = 0, runs = 0; u64 carry = 0;
    for (u64 w : bits) { occ += __builtin_popcountll(w); runs += __builtin_popcountll(w & ~((w << 1) | carry)); carry = w >> 63; }
    double V = min((double)n, (double)occ * (double)Lw);
    return CM_HIT * ss.H + CM_WIN * runs + CM_BP * B * V;
}
static void equal_seeds(int m, int k, vector<Seed>& seeds) {
    int q = m / (k + 1); seeds.clear();
    for (int i = 0; i <= k; i++) seeds.push_back({i * q, min(q, SMAX)});
}
// windows -> merge -> stage-2 fallback check -> verify.  pred_f is the stage-1 prediction (recorded in the statistics).
static void execute_seeds(const uint8_t* T, size_t n, const Seq& P, int k, const MyersPat& mp, SeedSet& ss, double pred_f, double pred_s, vector<size_t>& out) {
    int m = P.size(), B = mp.B; double tg = now_ns();
    add_bad(get_index(), T, n, P, ss);
    S.pred_filter_ns += pred_f; S.pred_scan_ns += pred_s; S.filtered++; S.hits += ss.H;
    if (!g_nofallback && pred_f > pred_s) { S.gen_ns += now_ns() - tg; scan_fallback(mp, T, n, k, out); return; }
    const Index& ix = get_index();
    vector<Win> W; W.reserve(ss.H);
    for (size_t i = 0; i < ss.seeds.size(); i++) {
        if (!ss.ok[i]) continue;
        ll a = ss.seeds[i].a;
        for (uint32_t x = ss.lo[i]; x < ss.hi[i]; x++) { ll p = ix.pos[x]; W.push_back({max(0LL, p - a - k), min((ll)n, p - a + m + k)}); }
        for (uint32_t p : ss.bad[i]) W.push_back({max(0LL, (ll)p - a - k), min((ll)n, (ll)p - a + m + k)});
    }
    merge_windows(W);
    ll total = 0; for (auto& w : W) total += w.second - w.first;
    S.gen_ns += now_ns() - tg;
    if (!g_nofallback && CM_HIT * ss.H + CM_WIN * W.size() + CM_BP * B * (double)total > pred_s) { scan_fallback(mp, T, n, k, out); return; }
    verify_merged(W, mp, T, k, out);
}
// expected candidate windows W and verified bases V of an s-hit filter under uniformly placed hits: a band of k + 1 diagonals contains hits of
// seed i with probability 1 - exp(-h_i (k + 1) / n); the number of bands with at least s distinct seeds is a Poisson-binomial tail, times n / (k + 1),
// plus the read's own locus.  Repeats make real hits cluster, which only lowers the true count (the estimate is conservative there).
static int S_FACTOR_X100 = 200;   // empirical multiplier on the predicted window count (env APPROX_SFACTOR, percent)
static void predict_s(const SeedSet& ss, int m, int k, int s, size_t n, double& W, double& V) {
    double b = k + 1.0;
    double dp[8] = {1, 0, 0, 0, 0, 0, 0, 0};
    for (size_t i = 0; i < ss.seeds.size(); i++) {
        if (!ss.ok[i]) continue;
        double mu = (double)(ss.hi[i] - ss.lo[i]) * b / (double)n, pr = 1.0 - exp(-mu);
        for (int j = s; j >= 1; j--) dp[j] = (j == s ? dp[j] : dp[j] * (1 - pr)) + dp[j - 1] * pr;
        dp[0] *= (1 - pr);
    }
    W = (double)n / b * dp[s] * (S_FACTOR_X100 / 100.0) + 1.0;
    V = min((double)n, W * (m + 3.0 * k));
}
// candidates from a prepared seed set (k + s seeds): sort hits, two-pointer band pass, merged windows, optional fallback, verification
static void run_s_candidates(const uint8_t* T, size_t n, const Seq& P, int k, int s, const MyersPat& mp, SeedSet& ss, double sel_cost, double pred_s, vector<size_t>& out) {
    int m = P.size(), B = mp.B, J = (int)ss.seeds.size();
    const Index& ix = get_index();
    add_bad(ix, T, n, P, ss);
    S.hits += ss.H; S.filtered++; S.pred_scan_ns += pred_s;
    if (!g_nofallback && CM_HIT * ss.H > pred_s) { scan_fallback(mp, T, n, k, out); return; }
    double tg = now_ns();
    static thread_local vector<u64> keys; static thread_local vector<int> cnt;
    keys.clear(); keys.reserve((size_t)ss.H);
    for (size_t i = 0; i < ss.seeds.size(); i++) {
        if (!ss.ok[i]) continue;
        ll a = ss.seeds[i].a;
        for (uint32_t x = ss.lo[i]; x < ss.hi[i]; x++) keys.push_back((u64)((ll)ix.pos[x] - a + m) << 8 | (u64)i);
        for (uint32_t p : ss.bad[i]) keys.push_back((u64)((ll)p - a + m) << 8 | (u64)i);
    }
    sort(keys.begin(), keys.end());
    cnt.assign(J, 0);
    vector<Win> W; int distinct = 0; size_t j = 0;
    for (size_t i = 0; i < keys.size(); i++) {
        ll di = (ll)(keys[i] >> 8); int pi = (int)(keys[i] & 255);
        if (cnt[pi]++ == 0) distinct++;
        while (di - (ll)(keys[j] >> 8) > k) { if (--cnt[(int)(keys[j] & 255)] == 0) distinct--; j++; }
        if (distinct >= s) {
            ll dlo = (ll)(keys[j] >> 8) - m, dhi = di - m;
            ll wlo = max(0LL, dlo - k), whi = min((ll)n, dhi + m + k);
            if (!W.empty() && wlo <= W.back().second) W.back().second = max(W.back().second, whi);
            else W.push_back({wlo, whi});
        }
    }
    ll total = 0; for (auto& w : W) total += w.second - w.first;
    S.gen_ns += now_ns() - tg;
    double cost = CM_HIT * ss.H + CM_WIN * W.size() + CM_BP * B * (double)total + sel_cost;
    S.pred_filter_ns += cost;
    if (g_debug > 2) fprintf(stderr, "SMEAS s=%d H=%lld windows=%zu bp=%lld\n", s, ss.H, W.size(), total);
    if (!g_nofallback && cost > pred_s) { scan_fallback(mp, T, n, k, out); return; }
    verify_merged(W, mp, T, k, out);
}
// mode 0 = equal split, 1 = frequency-aware (SELMODE), 2 = hits + lambda * windows, 3 = auto (cost-model choice per query), 4 = original DP for comparison
static void filter_core(const uint8_t* T, size_t n, const Seq& P, int k, int mode, vector<size_t>& out) {
    int m = P.size();
    MyersPat mp(P);
    int B = mp.B;
    int q = m / (k + 1);
    if (q < 1) { scan_fallback(mp, T, n, k, out); return; }
    const Index& ix = get_index();
    double pred_s = CM_SCAN * (double)n * B;
    SeedSet ss;
    auto select = [&](int how, SeedSet& dst) -> bool {   // how: 0 original DP, 1/2 fast (cap ignored for 1)
        double ts = now_ns(); bool ok;
        if (how == 2) ok = select_fa(ix, P, k, dst.seeds, true);
        else if (how == 0) ok = select_fa(ix, P, k, dst.seeds, false);
        else {
            int cap = (how == 1) ? SMAX : SELCAP;
            ok = select_fa_fast(ix, P, k, dst.seeds, cap);
            double dt = now_ns() - ts; g_sel_nsop = 0.9 * g_sel_nsop + 0.1 * (dt / sel_ops(m, k, cap));
        }
        S.sel_ns += now_ns() - ts;
        return ok;
    };
    if (mode == 0) {
        equal_seeds(m, k, ss.seeds); collect(ix, T, n, P, ss);
        execute_seeds(T, n, P, k, mp, ss, model_cost(ss.H, m, k, B, n), pred_s, out);
    } else if (mode == 1 || mode == 2 || mode == 4 || mode == 5) {
        if (mode == 5) {   // gated window objective: only repeat-heavy reads (max equal-split seed count far above the random expectation) pay for it
            SeedSet gq; equal_seeds(m, k, gq.seeds); collect_ranges(ix, P, gq);
            ll maxc = 0; for (size_t i = 0; i < gq.seeds.size(); i++) if (gq.ok[i]) maxc = max(maxc, (ll)(gq.hi[i] - gq.lo[i]));
            int s0 = min(q, SMAX); double lam = (double)n / (double)((u64)1 << (2 * s0));
            bool heavy = maxc >= GATE_MIN && (double)maxc >= GATE_RATIO * max(1.0, lam);
            if (g_debug > 1) fprintf(stderr, "GATE m=%d k=%d maxc=%lld lam=%.1f heavy=%d\n", m, k, maxc, lam, (int)heavy);
            if (heavy) S.auto_fa++; else S.auto_skip++;
            mode = heavy ? 2 : 1;
        }
        int how = mode == 2 ? 2 : (mode == 4 ? 0 : (SELMODE == 0 ? 0 : (SELMODE == 1 ? 1 : 3)));
        if (!select(how == 3 ? 3 : how, ss)) { scan_fallback(mp, T, n, k, out); return; }
        collect(ix, T, n, P, ss);
        double sel_term = (how == 0 || how == 2) ? CM_SEL * m * SMAX : g_sel_nsop * sel_ops(m, k, how == 1 ? SMAX : SELCAP);
        execute_seeds(T, n, P, k, mp, ss, model_cost(ss.H, m, k, B, n) + sel_term, pred_s, out);
    } else {   // mode 3: filter_auto
        static thread_local SeedSet eq, fa, sset[2];   // reused buffers: no allocation per query
        double cost_s = 0;
        equal_seeds(m, k, eq.seeds); collect_ranges(ix, P, eq);
        int pick = 0;   // 0 equal split, 1 frequency-aware, 2 scan
        const double INFC = 1e300;
        double cost_eq = INFC, cost_fa = INFC;
        const double LB_eq = CM_HIT * eq.H;   // lower bound of the equal-split cost: its hit cost alone
        if (!g_nofallback && LB_eq > pred_s) pick = 2;                     // hit cost alone already exceeds a scan
        else {
            double sel_est = g_sel_nsop * sel_ops(m, k, SELCAP);
            bool run_dp = true;
            if (LB_eq <= AUTO_SKIP * sel_est) {                              // equal split might be cheap: look at its real cost
                cost_eq = occupancy_cost(ix, eq, m, k, B, n);
                run_dp = cost_eq > AUTO_SKIP * sel_est;
            }
            if (run_dp) {
                if (select(3, fa)) { collect_ranges(ix, P, fa); if (g_nofallback || CM_HIT * fa.H <= pred_s) cost_fa = occupancy_cost(ix, fa, m, k, B, n); }
                if (cost_eq >= INFC && (cost_fa >= INFC || LB_eq < cost_fa)) cost_eq = occupancy_cost(ix, eq, m, k, B, n);   // equal split can only win if its bound is below cost_fa
                double best = cost_eq;
                if (cost_fa < best) { best = cost_fa; pick = 1; }
                if (!g_nofallback && pred_s < best) { pick = 2; best = pred_s; }
                if (AUTO_S && best > AUTO_S_TRIG * sel_est) {          // generalised pigeonhole: only worth looking at when the fixed filters are expensive
                    for (int sv = 2; sv <= 3; sv++) {
                        SeedSet& sx = sset[sv - 2];
                        if (m / (k + sv) < 1 || k + sv > 250) continue;
                        double ts = now_ns(); bool okx = select_fa_fast(ix, P, k + sv - 1, sx.seeds, SELCAP); S.sel_ns += now_ns() - ts;
                        if (!okx) continue;
                        collect_ranges(ix, P, sx);
                        if (CM_HIT * sx.H > best) continue;                   // hit cost alone already exceeds the best option
                        double Wp, Vp; predict_s(sx, m, k, sv, n, Wp, Vp);
                        double cs = CM_HIT * sx.H + CM_WIN * Wp + CM_BP * B * Vp;
                        if (g_debug) fprintf(stderr, "AUTOS s=%d H=%lld Wpred=%.1f cost=%.0f best=%.0f\n", sv, sx.H, Wp, cs, best);
                        if (cs < best) { best = cs; pick = sv + 1; cost_s = cs; }
                    }
                }
            } else {
                S.auto_skip++;
                if (!g_nofallback && pred_s < cost_eq) pick = 2;
            }
        }
        if (g_debug) fprintf(stderr, "AUTO m=%d k=%d H_eq=%lld cost_eq=%.0f H_fa=%lld cost_fa=%.0f scan=%.0f sel_est=%.0f pick=%d\n", m, k, eq.H, cost_eq, fa.H, cost_fa, pred_s, g_sel_nsop * sel_ops(m, k, SELCAP), pick);
        if (pick == 2) { S.auto_scan++; S.pred_filter_ns += min(min(cost_eq, cost_fa), 1e15); S.pred_scan_ns += pred_s; S.filtered++; S.hits += eq.H; scan_fallback(mp, T, n, k, out); }
        else if (pick == 3 || pick == 4) { int sv = pick - 1; (sv == 2 ? S.auto_s2 : S.auto_s3)++; run_s_candidates(T, n, P, k, sv, mp, sset[sv - 2], 0.0, pred_s, out); }
        else if (pick == 1) { S.auto_fa++; execute_seeds(T, n, P, k, mp, fa, cost_fa, pred_s, out); }
        else { S.auto_eq++; execute_seeds(T, n, P, k, mp, eq, cost_eq, pred_s, out); }
    }
}

// ---- generalised pigeonhole, s required hits (filter_s2, filter_s3, filter_fs2, filter_fs3)
// P is cut into J = k + s disjoint seeds.  An alignment with at most k edits leaves at least s seeds untouched (each edit touches at most one
// seed), so at least s DISTINCT seeds occur exactly, on diagonals (p - a) that differ by at most k (only indels between two untouched seeds shift
// the diagonal, and there are at most k edits in all).  So a candidate needs s distinct seeds with hits inside a diagonal band of width k;
// its window is [dmin - k, dmax + m + k).  s = 1 is the classic filter.  Hits are sorted by (diagonal, seed) and a two-pointer pass over the band
// counts distinct seeds; windows are produced in non-decreasing order and merged on the fly.
static void filter_s_core(const uint8_t* T, size_t n, const Seq& P, int k, int s, bool usefa, vector<size_t>& out) {
    int m = P.size(); MyersPat mp(P); int B = mp.B;
    int J = k + s;
    if (J > 250 || m / J < 1) { scan_fallback(mp, T, n, k, out); return; }
    const Index& ix = get_index();
    double pred_s = CM_SCAN * (double)n * B;
    static thread_local SeedSet ss;
    double sel_cost = 0;
    if (usefa) {
        double ts = now_ns(); bool ok = select_fa_fast(ix, P, J - 1, ss.seeds, SELCAP); S.sel_ns += now_ns() - ts;
        if (!ok) { scan_fallback(mp, T, n, k, out); return; }
        sel_cost = g_sel_nsop * sel_ops(m, J - 1, SELCAP);
    } else equal_seeds(m, J - 1, ss.seeds);
    collect_ranges(ix, P, ss);
    if (g_debug > 2) { double Wp, Vp; predict_s(ss, m, k, s, n, Wp, Vp); fprintf(stderr, "SPRED s=%d H=%lld W=%.1f V=%.0f\n", s, ss.H, Wp, Vp); }
    run_s_candidates(T, n, P, k, s, mp, ss, sel_cost, pred_s, out);
}
static void filter_s1(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) { filter_s_core(T, n, P, k, 1, false, out); }
static void filter_s2(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) { filter_s_core(T, n, P, k, 2, false, out); }
static void filter_s3(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) { filter_s_core(T, n, P, k, 3, false, out); }
static void filter_fs2(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) { filter_s_core(T, n, P, k, 2, true, out); }
static void filter_fs3(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) { filter_s_core(T, n, P, k, 3, true, out); }

static void filter_alg(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) { filter_core(T, n, P, k, 0, out); }
static void filter_fa(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) { filter_core(T, n, P, k, 1, out); }
static void filter_auto(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) { filter_core(T, n, P, k, 3, out); }
static void filter_fa_full(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) { filter_core(T, n, P, k, 4, out); }
static void filter_fa_wg(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) { filter_core(T, n, P, k, 5, out); }
static void filter_fa_w(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) { filter_core(T, n, P, k, 2, out); }

// ---------------------------------------------------------------- 5. BLAST-style heuristic
// Differences from real BLAST (blastn): exact 11-mer seeds only (no neighbourhood words / substitution matrix),
// ungapped X-drop extension with +2/-3 scoring and a fixed raw-score cutoff (no E-value / Karlin-Altschul
// statistics), no gapped X-drop DP: instead an HSP above the cutoff triggers an exact Myers scan of a
// window [diag-k, diag+m+k); no low-complexity masking; two-hit mode requires two NON-overlapping seeds on the
// SAME diagonal within A=40, as in BLAST.
static const int MATCH = 2, MISMATCH = -3, XDROP = 20, HSP_CUT = 34, TWOHIT_A = 40;

static int ungapped_extend(const uint8_t* T, ll n, const Seq& P, ll diag, int o, int w, int& end_o) {
    int m = P.size(); ll t0 = diag + o;
    int best = 0, cur = 0, bestI = o + w;
    for (int i = o + w; i < m && t0 + (i - o) < n; i++) {
        ll j = t0 + (i - o);
        cur += (P[i] == T[j] && T[j] < 4) ? MATCH : MISMATCH;
        if (cur > best) { best = cur; bestI = i + 1; } else if (best - cur >= XDROP) break;
    }
    end_o = bestI;
    int bestL = 0; cur = 0;
    for (int i = o - 1; i >= 0 && t0 - (o - i) >= 0; i--) {
        ll j = t0 - (o - i);
        cur += (P[i] == T[j] && T[j] < 4) ? MATCH : MISMATCH;
        if (cur > bestL) bestL = cur; else if (bestL - cur >= XDROP) break;
    }
    return w * MATCH + best + bestL;
}

static void blast_core(const uint8_t* T, size_t n, const Seq& P, int k, bool twohit, vector<size_t>& out) {
    int m = P.size();
    MyersPat mp(P);
    int w = SMAX;
    if (m < w + 4) { scan_fallback(mp, T, n, k, out); return; }  // too short to seed
    const Index& ix = get_index();
    vector<pair<ll, int>> hits;  // (diagonal, pattern offset)
    for (int o = 0; o + w <= m; o++) {
        Seed sd{o, w}; uint32_t lo, hi; u64 c;
        if (!seed_range(ix, P, sd, lo, hi, c)) continue;
        for (uint32_t x = lo; x < hi; x++) hits.push_back({(ll)ix.pos[x] - o, o});
    }
    sort(hits.begin(), hits.end());
    S.hits += hits.size();
    vector<Win> W;
    size_t i = 0;
    while (i < hits.size()) {
        size_t j = i; ll d = hits[i].first;
        int last_o = -1000000, ext_end = -1;
        for (; j < hits.size() && hits[j].first == d; j++) {
            int o = hits[j].second;
            if (o < ext_end) continue;                      // already covered by an extension on this diagonal
            bool trigger = true;
            if (twohit) {
                if (last_o >= 0 && o >= last_o + w && o - last_o <= TWOHIT_A) trigger = true;
                else { if (o >= last_o + w || last_o < 0) last_o = o; trigger = false; }
            }
            if (!trigger) continue;
            S.ungapped++;
            int e; int sc = ungapped_extend(T, (ll)n, P, d, o, w, e);
            ext_end = e;
            if (twohit) last_o = o;
            if (sc >= HSP_CUT) { S.hsps++; W.push_back({max(0LL, d - k), min((ll)n, d + m + k)}); }
        }
        i = j;
    }
    if (W.empty()) return;
    merge_windows(W);
    verify_merged(W, mp, T, k, out);
}
static void blast1(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) { blast_core(T, n, P, k, false, out); }
static void blast2(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) { blast_core(T, n, P, k, true, out); }


// ---------------------------------------------------------------- alignment output (traceback on verified windows)
// Given the exact set of end positions (distance <= k) from the filter, find the best alignment (minimum distance, ties:
// leftmost end) and a CIGAR.  Ends are grouped into clusters whose windows overlap; for each cluster a full semi-global DP
// (free start in the text) over the window [first_end + 1 - (m + k), last_end] is filled, the best end of the cluster is
// chosen, and a traceback gives the alignment.  CIGAR uses SAM semantics with extended operators:
//   =  match   X  mismatch   I  extra read (pattern) character   D  extra text (reference) character
struct Aln { int dist = -1; ll start = 0, end = 0; string cigar; size_t nends = 0; };

static int dp_window(const uint8_t* T, ll ws, ll we, const Seq& P, vector<int>& D) {
    int m = P.size(); ll W = we - ws;
    D.assign((size_t)(m + 1) * (W + 1), 0);
    for (int i = 0; i <= m; i++) D[(size_t)i * (W + 1)] = i;
    for (ll j = 1; j <= W; j++) {
        uint8_t t = T[ws + j - 1];
        for (int i = 1; i <= m; i++) {
            int cost = (P[i - 1] == t && t < 4) ? 0 : 1;
            int v = D[(size_t)(i - 1) * (W + 1) + j - 1] + cost;
            v = min(v, D[(size_t)(i - 1) * (W + 1) + j] + 1);
            v = min(v, D[(size_t)i * (W + 1) + j - 1] + 1);
            D[(size_t)i * (W + 1) + j] = v;
        }
    }
    return 0;
}

static void traceback(const vector<int>& D, const uint8_t* T, ll ws, ll W, const Seq& P, ll jend, Aln& a) {
    int m = P.size(); int i = m; ll j = jend; string ops;
    while (i > 0) {
        int cur = D[(size_t)i * (W + 1) + j];
        if (j > 0) {
            uint8_t t = T[ws + j - 1];
            int cost = (P[i - 1] == t && t < 4) ? 0 : 1;
            if (cur == D[(size_t)(i - 1) * (W + 1) + j - 1] + cost) { ops.push_back(cost ? 'X' : '='); i--; j--; continue; }
        }
        if (cur == D[(size_t)(i - 1) * (W + 1) + j] + 1) { ops.push_back('I'); i--; continue; }
        ops.push_back('D'); j--;   // only reachable when j > 0
    }
    reverse(ops.begin(), ops.end());
    a.start = ws + j;
    string cig; size_t x = 0;
    while (x < ops.size()) { size_t y = x; while (y < ops.size() && ops[y] == ops[x]) y++; cig += to_string(y - x); cig.push_back(ops[x]); x = y; }
    a.cigar = cig;
}

static bool align_best(const uint8_t* T, size_t n, const Seq& P, int k, const vector<size_t>& ends, Aln& best) {
    if (ends.empty()) return false;
    int m = P.size();
    best = Aln(); best.nends = ends.size();
    size_t a = 0;
    vector<int> D;
    while (a < ends.size()) {
        size_t b = a;
        while (b + 1 < ends.size() && (ll)(ends[b + 1] - ends[b]) <= (ll)(m + k)) b++;
        ll ws = max(0LL, (ll)ends[a] + 1 - (m + k)), we = (ll)ends[b] + 1;
        ll W = we - ws;
        if ((ll)(m + 1) * (W + 1) > 40000000LL) {   // very long cluster: align each end on its own window
            for (size_t c = a; c <= b; c++) {
                ll ws1 = max(0LL, (ll)ends[c] + 1 - (m + k)), we1 = (ll)ends[c] + 1;
                dp_window(T, ws1, we1, P, D);
                ll W1 = we1 - ws1; int d = D[(size_t)m * (W1 + 1) + W1];
                if (best.dist < 0 || d < best.dist) { best.dist = d; best.end = we1; Aln tmp; traceback(D, T, ws1, W1, P, W1, tmp); best.start = tmp.start; best.cigar = tmp.cigar; }
            }
        } else {
            dp_window(T, ws, we, P, D);
            for (size_t c = a; c <= b; c++) {
                ll jc = (ll)ends[c] + 1 - ws; int d = D[(size_t)m * (W + 1) + jc];
                if (best.dist < 0 || d < best.dist) {
                    best.dist = d; best.end = (ll)ends[c] + 1;
                    Aln tmp; traceback(D, T, ws, W, P, jc, tmp); best.start = tmp.start; best.cigar = tmp.cigar;
                }
            }
        }
        a = b + 1;
    }
    return best.dist >= 0;
}

// ---------------------------------------------------------------- driver
static long cur_rss_kb() {   // VmRSS now
    ifstream f("/proc/self/status"); string line;
    while (getline(f, line)) if (line.rfind("VmRSS:", 0) == 0) return atol(line.c_str() + 6);
    return -1;
}
static long peak_rss_kb() {   // VmHWM from /proc/self/status (peak resident set size)
    ifstream f("/proc/self/status"); string line;
    while (getline(f, line)) if (line.rfind("VmHWM:", 0) == 0) return atol(line.c_str() + 6);
    return -1;
}
static void noop_alg(const uint8_t*, size_t, const Seq&, int, vector<size_t>&) {}

static Seq read_text(const char* path) {
    ifstream f(path);
    if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
    Seq T; string line;
    while (getline(f, line)) {
        if (!line.empty() && line[0] == '>') continue;
        for (char c : line) if (!isspace((unsigned char)c)) T.push_back(code(c));
    }
    return T;
}

int main(int argc, char** argv) {
    if (argc < 5) { fprintf(stderr, "usage: %s algo text patterns k [--dump]\n", argv[0]); return 1; }
    string algo = argv[1]; int k = atoi(argv[4]);
    bool dump = false, times = false;
    for (int a = 5; a < argc; a++) { if (string(argv[a]) == "--dump") dump = true; if (string(argv[a]) == "--times") times = true; }
    if (getenv("APPROX_NOFALLBACK")) g_nofallback = atoi(getenv("APPROX_NOFALLBACK")) != 0;
    if (const char* wc = getenv("APPROX_WCAP")) WCAP = atoi(wc);
    if (const char* e = getenv("APPROX_LVSIMD")) LV_USE_SIMD = atoi(e);
    if (const char* e = getenv("APPROX_FAST1")) MYERS_FAST1 = atoi(e);
    if (const char* e = getenv("APPROX_WBIN")) WBIN = atoi(e);
    if (const char* e = getenv("APPROX_AUTOS")) AUTO_S = atoi(e);
    if (const char* e = getenv("APPROX_AUTOSTRIG")) AUTO_S_TRIG = atof(e);
    if (const char* e = getenv("APPROX_SFACTOR")) S_FACTOR_X100 = atoi(e);
    if (const char* e = getenv("APPROX_GATE")) GATE_RATIO = atoi(e);
    if (const char* e = getenv("APPROX_GATEMIN")) GATE_MIN = atoi(e);
    if (const char* e = getenv("APPROX_SEL")) SELMODE = atoi(e);
    if (const char* e = getenv("APPROX_CAP")) SELCAP = atoi(e);
    if (const char* e = getenv("APPROX_OCCSMALL")) OCC_SMALL = atoi(e);
    if (const char* e = getenv("APPROX_DEBUG")) g_debug = atoi(e);
    if (const char* e = getenv("APPROX_AUTOSKIP")) AUTO_SKIP = atof(e);
    if (const char* cm = getenv("APPROX_CM")) sscanf(cm, "%lf,%lf,%lf,%lf,%lf", &CM_HIT, &CM_WIN, &CM_BP, &CM_SCAN, &CM_SEL);
    Seq T = read_text(argv[2]);
    vector<Seq> pats;
    { ifstream f(argv[3]); string line;
      while (getline(f, line)) { istringstream is(line); string tok; if (!(is >> tok)) continue;
          Seq p; for (char c : tok) p.push_back(code(c)); pats.push_back(p); } }
    gT = T.data(); gN = T.size();
    malloc_trim(0);
    long loaded_kb = cur_rss_kb();   // steady RSS after loading text and patterns (transient read buffers returned)

    typedef void (*Fn)(const uint8_t*, size_t, const Seq&, int, vector<size_t>&);
    Fn fn = algo == "sellers" ? sellers : algo == "banded" ? banded : algo == "myers" ? myers : algo == "myers_pruned" ? myers_pruned : algo == "lv" ? lv : algo == "lv_stream" ? lv_stream : algo == "lv_direct" ? lv_direct
          : algo == "filter" ? filter_alg : algo == "filter_fa" ? filter_fa : algo == "filter_fa_w" ? filter_fa_w : algo == "filter_fa_wg" ? filter_fa_wg : algo == "filter_auto" ? filter_auto : algo == "filter_s1" ? filter_s1 : algo == "filter_s2" ? filter_s2 : algo == "filter_s3" ? filter_s3 : algo == "filter_fs2" ? filter_fs2 : algo == "filter_fs3" ? filter_fs3 : algo == "filter_fa_full" ? filter_fa_full
          : algo == "blast1" ? blast1 : algo == "blast2" ? blast2 : algo == "map" ? filter_alg : algo == "myers_best" ? myers_best : algo == "mindist" ? mindist : algo == "noop" ? noop_alg : nullptr;
    if (!fn) { fprintf(stderr, "unknown algo %s\n", algo.c_str()); return 1; }

    // index build is timed separately (one-time cost, amortised over queries)
    double t0 = now_ns();
    if (algo == "filter" || algo == "filter_fa" || algo == "filter_fa_w" || algo == "filter_fa_wg" || algo == "filter_auto" || algo == "filter_fa_full" || algo == "filter_s1" || algo == "filter_s2" || algo == "filter_s3" || algo == "filter_fs2" || algo == "filter_fs3" || algo == "blast1" || algo == "blast2") get_index();
    double t1 = now_ns();

    vector<vector<size_t>> res(pats.size());
    struct PerQ { double ns, hits, bp, pf, ps, fb; };
    vector<PerQ> perq(pats.size());
    for (size_t i = 0; i < pats.size(); i++) {
        Stats b = S; double ts = now_ns();
        fn(T.data(), T.size(), pats[i], k, res[i]);
        sort(res[i].begin(), res[i].end());
        res[i].erase(unique(res[i].begin(), res[i].end()), res[i].end());
        double te = now_ns();
        perq[i] = {te - ts, (double)(S.hits - b.hits), (double)(S.verified_bp - b.verified_bp),
                   S.pred_filter_ns - b.pred_filter_ns, S.pred_scan_ns - b.pred_scan_ns, (double)(S.fallbacks - b.fallbacks)};
    }
    double t2 = now_ns();
    size_t total = 0; for (auto& r : res) total += r.size();
    printf("RESULT algo=%s n=%zu patterns=%zu k=%d matches=%zu index_ms=%.3f time_ms=%.3f hits=%lld windows=%lld "
           "verified_bp=%lld fallbacks=%lld filtered=%lld gen_ms=%.3f ver_ms=%.3f sel_ms=%.3f pred_filter_ms=%.3f pred_scan_ms=%.3f "
           "cells_banded=%lld ungapped=%lld hsps=%lld lv_cells=%lld pruned_blocks=%lld peak_kb=%ld loaded_kb=%ld auto_skip=%lld auto_eq=%lld auto_fa=%lld auto_scan=%lld auto_s2=%lld auto_s3=%lld\n",
           algo.c_str(), T.size(), pats.size(), k, total, (t1 - t0) / 1e6, (t2 - t1) / 1e6, S.hits, S.windows,
           S.verified_bp, S.fallbacks, S.filtered, S.gen_ns / 1e6, S.ver_ns / 1e6, S.sel_ns / 1e6, S.pred_filter_ns / 1e6,
           S.pred_scan_ns / 1e6, S.cells_banded, S.ungapped, S.hsps, lv_cells, pruned_block_steps, peak_rss_kb(), loaded_kb, S.auto_skip, S.auto_eq, S.auto_fa, S.auto_scan, S.auto_s2, S.auto_s3);
    if (dump) for (size_t i = 0; i < res.size(); i++) for (size_t e : res[i]) printf(algo == "mindist" ? "D %zu %zu\n" : "M %zu %zu\n", i, e);
    if (algo == "map") for (size_t i = 0; i < pats.size(); i++) {
        Aln al;
        if (align_best(T.data(), T.size(), pats[i], k, res[i], al))
            printf("A %zu %d %lld %lld %s %zu\n", i, al.dist, al.start, al.end, al.cigar.c_str(), al.nends);
        else printf("A %zu * * * * 0\n", i);
    }
    if (times) for (size_t i = 0; i < perq.size(); i++)
        printf("T %zu %.0f %.0f %.0f %.1f %.1f %.0f\n", i, perq[i].ns, perq[i].hits, perq[i].bp, perq[i].pf, perq[i].ps, perq[i].fb);
    return 0;
}
