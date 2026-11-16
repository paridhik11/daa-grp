// Approximate DNA matching (k-differences / edit distance): exact algorithms and heuristics in one driver.
//
//   sellers    : Sellers (1980) DP, O(mn)                                  baseline / oracle
//   banded     : Ukkonen (1985) cutoff DP, O(kn)
//   myers      : Myers (1999) bit-vector (blocks), O(ceil(m/w) n)
//   myers_pruned: Myers with Ukkonen cutoff on blocks (only blocks that can hold a cell <= k are computed)
//   lv         : Landau-Vishkin (1989) diagonal algorithm, O(nk) expected (direct character extension)
//   filter     : pigeonhole filter, k+1 EQUAL-length pieces  + Myers verification, cost-model fallback
//   filter_fa  : pigeonhole filter, k+1 pieces chosen by DP to minimise total seed occurrences
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
    double gen_ns = 0, ver_ns = 0, sel_ns = 0, pred_filter_ns = 0, pred_scan_ns = 0;
} S;
static bool g_nofallback = false;
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

static void myers_scan(const MyersPat& mp, const uint8_t* T, size_t len, int k, size_t base, vector<size_t>& out) {
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
static void lv(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) {
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
                    diag.clear();
                    for (uint32_t x = lo; x < hi; x++) diag.push_back((ll)ix.pos[x] - a);
                    sort(diag.begin(), diag.end());
                    w = 1;
                    for (size_t x = 1; x < diag.size(); x++) if (diag[x] - diag[x - 1] > Lw) w++;
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

static void filter_core(const uint8_t* T, size_t n, const Seq& P, int k, int mode, vector<size_t>& out) {
    bool freq_aware = mode >= 1;
    int m = P.size();
    MyersPat mp(P);
    int B = mp.B;
    int q = m / (k + 1);
    if (q < 1) { scan_fallback(mp, T, n, k, out); return; }
    const Index& ix = get_index();
    vector<Seed> seeds;
    if (freq_aware) {
        double ts = now_ns();
        bool okfa = select_fa(ix, P, k, seeds, mode == 2);
        S.sel_ns += now_ns() - ts;
        if (!okfa) { scan_fallback(mp, T, n, k, out); return; }
    } else for (int i = 0; i <= k; i++) seeds.push_back({i * q, min(q, SMAX)});

    double tg = now_ns();
    // hits: positions from the index range, plus brute-force check of the (few) un-indexed positions
    ll H = 0;
    vector<uint32_t> lo(seeds.size()), hi(seeds.size()); vector<char> ok(seeds.size());
    vector<vector<uint32_t>> badhits(seeds.size());
    for (size_t i = 0; i < seeds.size(); i++) {
        u64 c; ok[i] = seed_range(ix, P, seeds[i], lo[i], hi[i], c);
        if (!ok[i]) continue;
        H += hi[i] - lo[i];
        for (uint32_t p : ix.bad) {
            if (p + seeds[i].len > n) continue;
            bool eq = true;
            for (int t = 0; t < seeds[i].len; t++) if (T[p + t] != P[seeds[i].a + t] || T[p + t] > 3) { eq = false; break; }
            if (eq) badhits[i].push_back(p);
        }
        H += badhits[i].size();
    }
    // cost model, stage 1 (before generating windows): expected merged windows/bp under uniform hit placement
    double Lw = m + 2.0 * k, lam = (double)H * Lw / (double)n;
    double V_est = n * (1.0 - exp(-lam)), W_est = H * exp(-lam);
    double pred_f = CM_HIT * H + CM_WIN * W_est + CM_BP * B * V_est + (freq_aware ? CM_SEL * m * SMAX : 0.0);
    double pred_s = CM_SCAN * (double)n * B;
    S.pred_filter_ns += pred_f; S.pred_scan_ns += pred_s; S.filtered++; S.hits += H;
    if (!g_nofallback && pred_f > pred_s) { S.gen_ns += now_ns() - tg; scan_fallback(mp, T, n, k, out); return; }

    vector<Win> W; W.reserve(H);
    for (size_t i = 0; i < seeds.size(); i++) {
        if (!ok[i]) continue;
        ll a = seeds[i].a;
        for (uint32_t x = lo[i]; x < hi[i]; x++) { ll p = ix.pos[x]; W.push_back({max(0LL, p - a - k), min((ll)n, p - a + m + k)}); }
        for (uint32_t p : badhits[i]) W.push_back({max(0LL, (ll)p - a - k), min((ll)n, (ll)p - a + m + k)});
    }
    merge_windows(W);
    ll total = 0; for (auto& w : W) total += w.second - w.first;
    S.gen_ns += now_ns() - tg;
    // stage 2 with exact merged size
    if (!g_nofallback && CM_HIT * H + CM_WIN * W.size() + CM_BP * B * (double)total > pred_s) { scan_fallback(mp, T, n, k, out); return; }
    verify_merged(W, mp, T, k, out);
}
static void filter_alg(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) { filter_core(T, n, P, k, 0, out); }
static void filter_fa(const uint8_t* T, size_t n, const Seq& P, int k, vector<size_t>& out) { filter_core(T, n, P, k, 1, out); }
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
    if (const char* cm = getenv("APPROX_CM")) sscanf(cm, "%lf,%lf,%lf,%lf,%lf", &CM_HIT, &CM_WIN, &CM_BP, &CM_SCAN, &CM_SEL);
    Seq T = read_text(argv[2]);
    vector<Seq> pats;
    { ifstream f(argv[3]); string line;
      while (getline(f, line)) { istringstream is(line); string tok; if (!(is >> tok)) continue;
          Seq p; for (char c : tok) p.push_back(code(c)); pats.push_back(p); } }
    gT = T.data(); gN = T.size();

    typedef void (*Fn)(const uint8_t*, size_t, const Seq&, int, vector<size_t>&);
    Fn fn = algo == "sellers" ? sellers : algo == "banded" ? banded : algo == "myers" ? myers : algo == "myers_pruned" ? myers_pruned : algo == "lv" ? lv
          : algo == "filter" ? filter_alg : algo == "filter_fa" ? filter_fa : algo == "filter_fa_w" ? filter_fa_w
          : algo == "blast1" ? blast1 : algo == "blast2" ? blast2 : algo == "map" ? filter_alg : algo == "myers_best" ? myers_best : algo == "mindist" ? mindist : nullptr;
    if (!fn) { fprintf(stderr, "unknown algo %s\n", algo.c_str()); return 1; }

    // index build is timed separately (one-time cost, amortised over queries)
    double t0 = now_ns();
    if (algo == "filter" || algo == "filter_fa" || algo == "filter_fa_w" || algo == "blast1" || algo == "blast2") get_index();
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
           "cells_banded=%lld ungapped=%lld hsps=%lld lv_cells=%lld pruned_blocks=%lld\n",
           algo.c_str(), T.size(), pats.size(), k, total, (t1 - t0) / 1e6, (t2 - t1) / 1e6, S.hits, S.windows,
           S.verified_bp, S.fallbacks, S.filtered, S.gen_ns / 1e6, S.ver_ns / 1e6, S.sel_ns / 1e6, S.pred_filter_ns / 1e6,
           S.pred_scan_ns / 1e6, S.cells_banded, S.ungapped, S.hsps, lv_cells, pruned_block_steps);
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
