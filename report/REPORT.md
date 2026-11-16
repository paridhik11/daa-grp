# Approximate DNA Sequence Matching: Dynamic Programming, Bit-Parallelism and Filtering Compared

*Design and Analysis of Algorithms, end-term project, Problem 8. All measurements: `results/summary.md`, `results/bench.csv`, plots in `results/`. Reproduce everything with `make bench`.*

## 1. Problem Statement

**Input.** A text T[0..n) (a genome), a pattern P[0..m) (a read), both over {A, C, G, T}, and an integer k >= 0. Edit distance ed(X, Y) is the minimum number of substitutions, insertions and deletions turning X into Y.

**k-differences problem.** Report every end position e of T such that some substring of T ending at e has ed(P, substring) <= k. This "all ends within k" definition is the one every exact algorithm here implements and is checked against.

Real sequences differ by point mutations and indels, so exact matching is not enough, and the textbook dynamic program (Needleman-Wunsch 1970 [18], Smith-Waterman 1981 [19], and its approximate-matching form by Sellers 1980 [1]) costs Theta(mn), too slow for many reads against a genome. We study how far exact methods can be pushed (banding [2], bit-parallelism [3], pigeonhole filtering [4, 5]), and what accuracy a BLAST-style heuristic [6] gives up for speed. Worst-case speedups beyond the DP are unlikely: edit distance has no O(n^(2-epsilon)) algorithm unless the Strong Exponential Time Hypothesis fails [7]. Gains come from small k / m, which is the regime of read mapping.

**Data.** Real genomes fetched from GitHub raw files (NCBI, EBI ENA, UCSC and Ensembl were unreachable from the sandbox, so *E. coli* K-12 MG1655 could not be obtained; no dataset here is called K-12):

| Dataset | Length (bp) | GC | Source |
|---|---|---|---|
| *E. coli* 536 (NC_008253.1), main text | 4,938,920 | 0.505 | BenLangmead/bowtie `genomes/NC_008253.fna` |
| SARS-CoV-2 isolate MT192765.1 | 29,829 | 0.380 | nf-core/test-datasets |
| phage lambda (NC_001416.1) | 48,502 | 0.499 | BenLangmead/bowtie2 example |
| human mitochondrial DNA | 16,569 | 0.444 | lh3/minimap2 test |

Headers, lengths, GC and SHA-256 are in `data/SOURCES.txt`; non-ACGT characters were stripped. Reads are simulated from the genome with exactly k random edits, so a true match within k always exists; the oracle for correctness is Sellers. A synthetic 4.6 Mbp genome (random, with injected repeats) is used only to fit the cost model.

**Scope and novelty.** This is a comparative study plus an improved variant of a known technique. The frequency-aware splitting (Section 2.5) follows the idea of the Optimal Seed Solver [13] (pick pigeonhole seeds with the fewest occurrences by dynamic programming) and relates to optimal search schemes [14]. Our differences are engineering choices: seeds capped at 11 characters so that one direct-address table answers every lookup, gaps allowed between seeds, verification by Myers windows, and a calibrated cost model deciding when to fall back to a full scan. We do not claim a new algorithm.

## 2. Algorithm

All algorithms are in `src/approx.cpp` (C++17). Notation: text T[0..n), pattern P[0..m), threshold k, alphabet {A,C,G,T}; any other character (N) matches nothing. w = 64 is the machine word size, B = ceil(m/w).

### 2.1 Sellers DP (baseline and oracle)
```
D[0][j] = 0 for all j ; D[i][0] = i
D[i][j] = min( D[i-1][j-1] + [P_i != T_j],  D[i-1][j] + 1,  D[i][j-1] + 1 )
report end j when D[m][j] <= k          (one column kept, O(m) space)
```

### 2.2 Ukkonen cutoff (banded) DP
```
lact = min(m, k)                                   // last row whose value is <= k
for each text char t:
    for i = 1 .. min(lact+1, m): C[i] = min(diag + [P_i != t], C[i] + 1, C[i-1]_new + 1)
    lact = min(lact+1, m); while lact > 0 and C[lact] > k: lact--
    if lact == m: report
```
Rows below `lact + 1` are never computed.

### 2.3 Myers bit-vector (blocks of 64 rows)
State per block: `Pv` (rows whose vertical delta is +1), `Mv` (rows whose vertical delta is -1). `Eq` is the match mask of the current text character, `hin` the horizontal delta entering the block from above (0 for block 0, because D[0][j] = 0).
```
advance_block(Pv, Mv, Eq, hin):                    // returns hout
    Xv = Eq | Mv
    if hin < 0: Eq |= 1
    Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq
    Ph = Mv | ~(Xh | Pv)          Mh = Pv & Xh
    hout = +1 if Ph has the block's top bit, -1 if Mh has it, else 0
    Ph <<= 1 ; Mh <<= 1 ; if hin < 0: Mh |= 1 elif hin > 0: Ph |= 1
    Pv = Mh | ~(Xv | Ph)          Mv = Ph & Xv
score = m ; for each text char: hin = 0; for each block: hin = advance_block(...); score += hin
report when score <= k
```

### 2.4 Pigeonhole filter, equal split (`filter`)
```
q = floor(m / (k+1)); seeds = { P[i*q .. i*q + min(q,11)) : i = 0..k }       // prefixes of k+1 equal pieces
H = total index hits of the seeds                                            // index: 11-mer table of T
if cost_model(H) predicts filter slower than scan: scan whole text with Myers (fallback)
windows = { [p - a - k, p - a + m + k) : seed at pattern offset a occurs at text position p }
merge overlapping windows; run Myers on each merged window
```
The index is a direct-address table over all 4^11 11-mers (CSR layout, 16 MB plus 4n bytes). A seed shorter than 11 characters is looked up as a contiguous range of the table; the few text positions that have no indexed 11-mer (text tail, or an N inside the next 11 characters) are checked directly.

### 2.5 Frequency-aware filter (`filter_fa`, our variant)
Instead of equal pieces, choose k+1 pairwise disjoint substrings of P (lengths 1 to 11, gaps allowed) that minimise the total number of occurrences in T.
```
cnt[a][L] = occurrences in T of P[a .. a+L)                 // one index lookup each, a < m, L <= 11
f[0][i] = 0
f[j][i] = min( f[j][i-1],                                   // character i-1 not used
               min_{1<=L<=min(11,i)}  f[j-1][i-L] + cnt[i-L][L] )   // last seed is P[i-L .. i)
choose the seeds by backtracking from f[k+1][m]; then proceed as in 2.4
```

### 2.6 Cost model and fallback rule
With H seed hits, window length L_w = m + 2k, and lambda = H L_w / n, the model assumes hits are spread uniformly: merged windows W ~ H e^(-lambda), verified length V ~ n (1 - e^(-lambda)).
```
t_filter = c_hit H + c_win W + c_bp B V   (+ c_sel m 11 for filter_fa)
t_scan   = a_scan n B
fall back to a full Myers scan when t_filter > t_scan
```
The constants (nanoseconds) are fitted by non-negative least squares (all constants >= 0) on the synthetic genome, with calibration points that cover the crossover (k up to 25), and evaluated on E. coli 536. The rule is applied once before windows are built (using the estimate) and once after (with the exact merged length).

### 2.7 BLAST-style heuristic (`blast1`, `blast2`)
```
for every 11-mer of P found in the index: record (diagonal = p - a, offset a)
one-hit mode: every hit triggers an ungapped extension
two-hit mode: trigger only when a second NON-overlapping hit lies on the SAME diagonal within A = 40 of the previous one
ungapped extension: score +2 per match, -3 per mismatch, extend both ways, stop when the score falls X = 20 below the best
if the extension score >= 34 (a hit-segment pair, HSP): run Myers on the window [diag - k, diag + m + k)
```
Differences from real BLAST (blastn): exact 11-mer seeds only (no neighbourhood words, no substitution matrix); a fixed raw-score cutoff instead of E-values and Karlin-Altschul statistics; no gapped X-drop alignment, an HSP instead triggers an exact Myers scan of a window (so every reported match is genuine, but recall depends only on the seeding and the cutoff); no low-complexity masking; no database-level statistics; DNA only. It is a model of BLAST's seed-and-extend logic, not a re-implementation.

### 2.8 Baselines from libraries
- **Edlib** (`pip install edlib`, HW mode, `task="locations"`): returns the minimum distance and only the end positions achieving it.
- **WFA2** (`pip install pywfa`, ends-free, Levenshtein): returns one best alignment per read.

Both therefore answer "best match per read", a different question from "all ends within k". They are timed (including about 0.1 ms of Python call overhead) and checked for consistency, never compared as all-match sets. A small case shows the difference (`tests/edlib_semantics.py`): for T = TTACGTTTTACGATT, P = ACGT, k = 1, Sellers reports ends {4, 5, 6, 11, 12, 13}, while Edlib reports distance 0 and only end {5}.

## 3. Proof of Correctness

Throughout, D[i][j] is the Sellers matrix: D[i][j] = min over starts s of ed(P[1..i], T[s+1..j]), with the empty substring allowed.

**Theorem 1 (Sellers).** The recurrence of 2.1 computes D.
*Proof.* Induction on i + j. For i = 0 the empty prefix matches the empty substring at cost 0; for j = 0 the only option is i deletions. For i, j >= 1 take an optimal alignment of P[1..i] with some T[s+1..j] and look at its last column. It is a match or mismatch of P_i with T_j (cost [P_i != T_j] plus an optimal alignment for (i-1, j-1)), or P_i deleted (cost 1 plus an optimal alignment for (i-1, j)), or T_j inserted (cost 1 plus an optimal alignment for (i, j-1)). Each remainder is again an alignment of a pattern prefix against a text substring with a free start, so it is counted by D. The minimum of the three is therefore D[i][j]. For the case s = j the alignment is only deletions, covered by the (i-1, j) term. Position j is reported iff D[m][j] <= k. QED

**Lemma 2 (neighbouring cells).** For all valid i, j: (a) |D[i][j] - D[i-1][j]| <= 1, (b) |D[i][j] - D[i][j-1]| <= 1, (c) 0 <= D[i][j] - D[i-1][j-1] <= 1.
*Proof.* "<= +1" in each direction is the recurrence (delete P_i, insert T_j, or add one match/mismatch column). For (a), D[i-1][j] <= D[i][j] + 1: in an optimal alignment for (i, j) delete the column of P_i. If P_i was matched or mismatched with a text character, that character becomes an insertion (cost +1 at most); if P_i was deleted the cost drops by 1. For (b), D[i][j-1] <= D[i][j] + 1: delete the column of T_j; if it was aligned to a pattern character that character becomes a deletion (+1 at most), if it was an insertion the cost drops by 1, and if T_j was the first character of the substring the substring becomes empty, which is allowed. For the lower bound of (c) take an optimal alignment of cost d for (i, j) and split on its last column. If it matches or mismatches P_i with T_j, dropping that column gives an alignment for (i-1, j-1) of cost at most d. If P_i is deleted, the rest is an alignment for (i-1, j) of cost d - 1, so D[i-1][j] <= d - 1 and by (b) D[i-1][j-1] <= D[i-1][j] + 1 <= d. If T_j is inserted, the rest is an alignment for (i, j-1) of cost d - 1, so D[i][j-1] <= d - 1 and by (a) D[i-1][j-1] <= D[i][j-1] + 1 <= d. QED

**Theorem 3 (Ukkonen cutoff).** A reports j iff D[m][j] <= k.
*Proof.* Let lact_j be the largest row with D[i][j] <= k (row 0 always qualifies). By Lemma 2(c), if D[i][j+1] <= k then D[i-1][j] <= k, so i - 1 <= lact_j; hence all cells of column j+1 with value <= k lie in rows [0, lact_j + 1], which are exactly the rows computed. Invariant: after column j, every stored C[i] is either the exact D[i][j] or a value greater than k. Inputs that are > k contribute at least k + 2 to the minimum, so they never create a wrong value <= k, and every output <= k is exact. A stored value <= k belongs to a row <= lact_j, which is recomputed in column j+1 (rows up to lact_j + 1 are recomputed), so stale values are always > k. Finally lact_j = m iff D[m][j] <= k. QED

**Theorem 4 (pigeonhole for ANY k+1 disjoint seeds, with window offset).** Let s_0, ..., s_k be k+1 pairwise disjoint substrings of P (a seed may be any substring, in particular a prefix of a piece), s_t = P[a_t .. a_t + len_t). If ed(P, T[i..j)) <= k then there is a t and a text position p with T[p..p + len_t) = s_t exactly, and the whole aligned substring lies inside W = [p - a_t - k, p - a_t + m + k).
*Proof.* Take an optimal alignment with at most k edit operations. An operation touches a seed if it substitutes or deletes a pattern character of that seed, or inserts a text character between two characters of that seed. Each operation touches at most one seed, since seeds are disjoint. So at most k seeds are touched and some seed s_t is untouched: its characters are aligned one-to-one to consecutive text characters, giving an exact occurrence at some p. The text aligned to P[0..a_t) has length a_t + delta1 and the text aligned to P[a_t + len_t .. m) has length (m - a_t - len_t) + delta2, where delta1 and delta2 are (insertions - deletions) in the two parts, so |delta1| + |delta2| <= k. Hence the start i = p - a_t - delta1 >= p - a_t - k and the end j = p + len_t + (m - a_t - len_t) + delta2 = p - a_t + m + delta2 <= p - a_t + m + k. QED

**Corollary 5 (the filters are exact).** `filter` and `filter_fa` report exactly the end positions reported by Sellers.
*Proof.* Completeness: by Theorem 4 every occurrence has a seed hit, and its window contains the whole alignment, so Myers on the (merged) window, which tries every start inside it, reports the same end. Soundness: windows are substrings of T, so any end reported by Myers on a window is the end of a substring of T within k edits. Merging windows only enlarges them. The seed lookup is exact: positions not covered by the 11-mer table are checked directly, and a seed containing N has no occurrence, because N matches nothing. The fallback is a full Myers scan, exact by Theorem 7. QED
The equal-split seeds (`min(q,11)`-character prefixes of k+1 pieces of length q) are disjoint, so Theorem 4 applies; the same holds for the DP-chosen seeds.

**Lemma 6 (the DP picks minimum-occurrence cuts).** Let cnt(a, L) be the number of occurrences of P[a..a+L) in T. Let f[j][i] be the minimum of the sum of cnt over all sets of j pairwise disjoint seeds of length at most 11 inside P[0..i). Then f satisfies f[0][i] = 0 and f[j][i] = min(f[j][i-1], min_{L} f[j-1][i-L] + cnt(i-L, L)), and f[k+1][m] is the minimum total number of seed hits over all valid (k+1)-seed choices. In particular the DP never has more hits than the equal split.
*Proof.* Induction on (j, i). An optimal set of j seeds inside P[0..i) either leaves character i-1 uncovered (then it lies inside P[0..i-1): value f[j][i-1]) or covers it. In the second case the seed that covers it ends exactly at i (the prefix ends there) and has some length L <= min(11, i); the other j-1 seeds are disjoint from it and lie inside P[0..i-L), and any optimal choice for them can be combined with this seed, so the value is f[j-1][i-L] + cnt(i-L, L). Taking the minimum over both cases gives the recurrence. The equal-split seeds are k+1 disjoint seeds of length at most 11, hence a feasible point, so the optimum is at most their total. QED
This optimises hit count, not running time (hits differ in how much they merge), and seeds longer than 11 are not considered because the index covers 11-mers.

**Theorem 7 (Myers bit-vector).** `advance_block` and the surrounding loop maintain: after processing T[1..j], (Pv, Mv) encode the vertical deltas of column j, and `score` = D[m][j].
*Setup.* For column j let dv_i = D[i][j] - D[i-1][j] and, for the step from column j-1, dh_i = D[i][j] - D[i][j-1]. By Lemma 2 all of these lie in {-1, 0, +1}. Row 0 has dh_0 = 0 because D[0][j] = 0 for all j. Bit i-1 of the word stands for row i. `Pv` has the bit set where dv = +1 in the old column (before the step) and `Mv` where dv = -1; `Eq` has the bit set where P_i = T_j (and no bit for N). Write pv_i, mv_i, e_i for those bits, and ph_i, mh_i for the indicators dh_i = +1, dh_i = -1 of the new step.
*Step 1: the cell rule.* Fix a cell (i, j) and write x = D[i-1][j-1]. Put a = D[i][j-1] - x = old dv_i and b = D[i-1][j] - x = dh_{i-1} (the horizontal delta entering the cell from above); both lie in {-1, 0, +1}. The recurrence is translation invariant, so D[i][j] - x = min([P_i != T_j], b + 1, a + 1). Because a, b >= -1 the last two terms are >= 0. Hence D[i][j] - x = 0 if e_i = 1 or a = -1 or b = -1, and 1 otherwise. The outputs are the new vertical delta D[i][j] - D[i-1][j] = (D[i][j] - x) - b and the horizontal delta dh_i = (D[i][j] - x) - a.
*Step 2: Boolean form of the outputs.* Let Xh_i = e_i OR mh_{i-1} (mh_0 = 0, or the injected hin, see Step 4), the condition "e_i or b = -1". Then D[i][j] - x = 0 exactly when Xh_i OR mv_i. Checking the cases: dh_i = +1 iff (D[i][j] - x = 0 and a = -1) or (D[i][j] - x = 1 and a = 0), i.e. ph_i = mv_i OR NOT(Xh_i OR pv_i) (when mv_i = 0 the second case is "not Xh and not pv"); dh_i = -1 iff (D[i][j] - x = 0 and a = +1), i.e. mh_i = pv_i AND Xh_i (pv_i and mv_i exclude each other). This is `Ph = Mv | ~(Xh | Pv)` and `Mh = Pv & Xh`. For the new vertical delta, with Xv_i = e_i OR mv_i: it is +1 iff (b = -1) or (D[i][j] - x = 1 and b = 0), i.e. new pv_i = mh_{i-1} OR NOT(Xv_i OR ph_{i-1}); it is -1 iff (D[i][j] - x = 0 and b = +1), i.e. new mv_i = ph_{i-1} AND Xv_i. Row i uses the horizontal deltas of row i-1, which is the shift by one bit: `Pv = (Mh<<1) | ~(Xv | (Ph<<1))` and `Mv = (Ph<<1) & Xv`.
*Step 3: the carry chain computes Xh.* From Step 2, mh_i = pv_i AND Xh_i, so Xh_{i+1} = e_{i+1} OR (pv_i AND Xh_i), a recurrence going down the rows. Let S = (Eq AND Pv) + Pv with carries c_i into bit i (c_1 = 0 for the lowest bit). The carry out of bit i is the majority of (e_i pv_i, pv_i, c_i), which simplifies to c_{i+1} = pv_i AND (e_i OR c_i). By induction c_i = mh_{i-1}: it holds for the lowest bit, and mh_i = pv_i AND (e_i OR mh_{i-1}) = pv_i AND (e_i OR c_i) = c_{i+1}. The sum bit is (e_i pv_i) XOR pv_i XOR c_i, so ((S XOR Pv) OR Eq)_i = ((e_i pv_i) XOR c_i) OR e_i, which equals 1 if e_i = 1 and equals c_i = mh_{i-1} if e_i = 0. That is exactly e_i OR mh_{i-1} = Xh_i. This is `Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq`.
*Step 4: blocks.* Rows are split into blocks of 64. The first row of a block receives hin = dh of the row above (0 for block 0 because dh_0 = 0). In Step 2 this replaces mh_{i-1} by [hin < 0] and ph_{i-1} by [hin > 0] for the lowest row of the block. In Step 3 it means a carry-in of [hin < 0] into the chain, which the code injects as `Eq |= 1` for the Xh computation only (Xv is computed from the original Eq first, because Xv_i = e_i OR mv_i does not involve mh), and in the new Pv/Mv through `Mh |= 1` or `Ph |= 1` after the shift. `hout` is read before the shift from the block's top row (bit 63, or bit (m-1) mod 64 in the last block), and is the hin of the next block. Rows above m in the last block are padding; information only flows to higher bit positions (carries and shifts go upward), so they cannot change any row <= m.
*Step 5: induction on columns and score.* Column 0 has D[i][0] = i, so every dv = +1: `Pv = all ones, Mv = 0`; the score is D[m][0] = m. The step j-1 to j, applied block by block (each block's hin is correct by the previous block's hout), maps the encoding of column j-1 to that of column j (Steps 1 to 4). The horizontal delta of row m is dh_m, the hout of the last block, so D[m][j] = D[m][j-1] + dh_m, which is how `score` is updated. Position j is reported iff score <= k, which by Theorem 1 is the same as Sellers. QED

**Proposition 8 (the BLAST-style heuristic is sound but incomplete).** Every match it reports is genuine (verification is an exact Myers scan of a window, a substring of T). It can miss a true match: take m = 100 and an occurrence that differs from P by a substitution at every 10th position (10 edits, k = 10). Every exact run shared with P has length 9, shorter than an 11-mer, so no seed hits and nothing is reported although Sellers reports the end. The measured recall is in Section 5 and Section 6.

## 4. Complexity Analysis

n = text length, m = pattern length, k = threshold, w = 64, B = ceil(m/w), L = 11 (seed cap), H = seed hits, V = total verified (merged window) length. Measured values are from *E. coli* 536 with m = 100 unless stated. Where a value is a median with a range, the range is the 10th to 90th percentile over 100 reads (the headline experiment); other experiments use means over 20 reads.

| Algorithm | Time (worst case) | Space | Measured behaviour |
|---|---|---|---|
| Sellers | Theta(mn) | O(m) | about 0.6 s per pattern (0.58 to 0.68 s), independent of k |
| Ukkonen banded | O(kn) expected, O(mn) worst | O(m) | about 2k + 1 cells per text character (3.1 at k = 1, 20 at k = 10, 59 at k = 30) |
| Myers | Theta(B n) | O(B sigma) | median 52 ms per pattern at m = 100 (B = 2), flat in k; 5.2 ns per block-character, 10 ns per text character |
| `filter` (equal split) | build O(n + 4^L); query O(k) lookups + O(H log H) + O(V B) | 4(4^L + 1) + 4n bytes (about 37 MB here) | median 0.007 to 0.036 ms for k <= 8, 0.46 ms at k = 10, then rises to Myers cost |
| `filter_fa` | as `filter` plus seed selection O(m L^2 + k m L) | same | selection costs about 0.02 to 0.05 ms per query (m = 100) |
| BLAST-style | seed lookups O(m) + sort of hits + extensions + O(V B) | same index | 0.02 to 0.08 ms per query, no accuracy guarantee |

(The O(kn) diagonal-transition algorithm of Landau and Vishkin [8] is related to banding but is not implemented here.)

**Seed-selection DP.** The count table needs m * L lookups, each of cost O(L) to form the seed code, so O(m L^2); the DP has (k+1) * m states with L transitions each, so O(k m L). For m = 100, k = 10, L = 11 this is about 11,000 transitions, measured at about 0.02 to 0.05 ms.

**Expected hits.** In a uniformly random text a seed of length s = min(floor(m/(k+1)), L) has expected (k+1) n 4^(-s) occurrences. For k = 10 (s = 9) that is about 207 hits for *E. coli* 536; we measure a median of 326 for the equal split and 212 for `filter_fa`. The measured value is above the random-text figure for the equal split because the genome is not random (repeats, composition bias) and because true occurrences add hits. The filter is selective while the hit count times the window length stays well below n, and the cost model of Section 2.6 predicts the crossover with the scan.

**Index build versus query time and break-even.** The index (4^11 table plus positions) takes 90 to 150 ms to build for *E. coli* 536 (about 20 ms for the 17 to 49 kbp genomes, where clearing the 16.8 MB table dominates). A query costs a median 0.007 to 0.46 ms (k <= 10) against 52 ms for a Myers scan. The break-even number of queries Q* = index time / (Myers median time - filter median time) is therefore small (medians over 100 reads; the index time is one measurement per row):

| k | index (ms) | Myers (ms/query) | filter (ms/query) | Q* filter | filter_fa (ms/query) | Q* filter_fa |
|---|---|---|---|---|---|---|
| 1 | 112 | 51.9 | 0.007 | 2.2 | 0.022 | 1.8 |
| 6 | 118 | 52.4 | 0.025 | 2.3 | 0.033 | 1.7 |
| 10 | 110 | 52.1 | 0.463 | 2.1 | 0.334 | 2.1 |
| 14 | 120 | 53.1 | 23.6 | 4.1 | 9.76 | 2.2 |
| 15 | 141 | 53.4 | 27.7 | 5.5 | 17.0 | 3.5 |
| 17 | 114 | 52.8 | 53.3 | never | 38.3 | 6.4 |
| 20 | 89 | 52.5 | 53.3 | never | 52.8 | never |

Two to three queries pay for the index for k <= 13, and at most 5.5 for k <= 15; at k >= 17 the equal-split filter never beats repeated scans (its fallback makes it equal to Myers plus small overhead), and `filter_fa` needs 6.4 queries at k = 17. Against Edlib (24 to 32 ms per query in the main sweep, 38 ms in one run, best-match semantics) the break-even is about 3.5 queries at k = 5 and 4.3 at k = 10 (index 121 to 134 ms). See `results/breakeven.png`. The numbers are for one fixed genome; a single query against a genome you will not query again is better served by Edlib or Myers directly.

**Cost model and crossover.** Section 2.6 gives the model. Its constants are fitted by non-negative least squares (c >= 0) on the synthetic genome, with calibration points covering the crossover (k up to 25 for m = 100 and 200), and evaluated on *E. coli* 536, a different text:

| constant | value |
|---|---|
| c_hit (per seed hit) | 82.1 ns |
| c_win (per merged window) | 404 ns |
| c_bp (per block-base verified) | 5.10 ns |
| a_scan (per block-character of a full scan) | 5.21 ns |
| c_sel (seed selection, per m * 11) | 45.2 ns |

NNLS did not set c_win to 0: the unconstrained least-squares estimate (404 ns) was already non-negative, so the sign constraint was inactive. Predicted against measured median time of the pure filter (fallback off), over 100 reads per k: median factor error 1.21x (worst 1.24x) for `filter` over 14 values of k, 1.24x (worst 2.38x) for `filter_fa` over 15 values, and 1.02x for the Myers scan. The worst `filter_fa` case is at small k, where the model charges about 0.05 ms for seed selection although the whole measured query takes 0.02 to 0.03 ms. The predicted crossover with a Myers scan is k/m = 0.166 (k = 16.6) for `filter` against 0.169 measured, and 0.184 against 0.191 for `filter_fa` (`results/costmodel.png`).

The fit is sensitive to the calibration range, and this matters. A first NNLS refit that used only k <= 15 gave c_win = 1,123 ns and c_bp = 2.56 ns (the number of windows and the number of bases verified are collinear, so the split between them is poorly determined). That model underestimated verification cost near the crossover: on the 100-read headline run at k = 20 only 13 of 100 reads fell back to a scan, and the filter's median was 84 ms against 53 ms for Myers, so the fallback did not do its job (`results/evidence/`). Earlier unconstrained fits on the same data had given c_win = 1,073, 339 and (after clipping a negative estimate) 0, which is why an earlier version of this report said c_win was zero. Widening the calibration range to k <= 25 gives the constants above, and with them the fallback fires for 99 to 100 of 100 reads at k >= 17 and the filter stays at Myers speed (53 ms). The range was widened after seeing the first fit fail on *E. coli* 536, so the evaluation is not fully blind to the choice of calibration range. The constants are machine specific; on other hardware re-fit with `scripts/calibrate.py`.

**Frequency-aware versus equal split.** Pure filter (fallback off), median over 100 reads. `filter_fa` has 1.4 to 9.0 times fewer seed hits than the equal split at every k tested (at least 1.5 times up to k = 20), and verifies 1.0 to 9.7 times fewer bases (the gap closes as the filter approaches a full scan). In time it loses when queries are cheap and wins when they are expensive:

| k | hits equal | hits fa | time equal (ms) | time fa (ms) | fa speedup |
|---|---|---|---|---|---|
| 1 | 5 | 1 | 0.006 | 0.022 | 0.29x |
| 4 | 14 | 3 | 0.018 | 0.036 | 0.50x |
| 8 | 24 | 15 | 0.033 | 0.061 | 0.55x |
| 10 | 326 | 212 | 0.458 | 0.328 | 1.40x |
| 12 | 4,828 | 1,656 | 6.47 | 2.41 | 2.69x |
| 15 | 23,647 | 14,000 | 27.2 | 18.5 | 1.46x |
| 20 | 435,776 | 144,428 | 82.9 | 59.3 | 1.40x |

`filter_fa` is slower for k <= 8 at m = 100 (0.29x to 0.76x, that is 1.3 to 3.4 times slower), because the selection step (about 0.02 to 0.05 ms) costs more than the verification it saves, and faster for k from 10 to 30 (1.35x to 2.7x). Across read lengths at 5% error it loses for m <= 100 (m = 50: 0.008 against 0.024 ms; m = 100: 0.032 against 0.054 ms) and wins for m >= 200 (m = 200, k = 10: 1.7x; m = 400, k = 20: 3.4x; m = 800, k = 40: 8.3x, 7.5 ms against 0.90 ms; means over 20 reads). Lemma 6 minimises hits, not verified length, and that is a real limitation: in the `lowcomplex_edge` stress case at k = 10, `filter_fa` has 4.5 times fewer hits (2,839 against 12,693) but they fall in 2.9 times more distinct merged windows (1,536 against 531), so it verifies 2.6 times more bases (193,687 against 75,968 per read) and is slower (2.42 ms against 1.19 ms). The equal split's hits cluster in the same low-complexity region and merge into few windows, while the DP's rarer seeds are scattered.

**Banded versus Myers: cells updated per second.** Counting cells computed by the banded algorithm and DP cells represented by Myers (m * n per pattern; one word operation updates 64 of them):

| k | banded cells computed | banded cells/s | Myers cells represented | Myers cells/s | ratio | banded ms | Myers ms |
|---|---|---|---|---|---|---|---|
| 1 | 3.1e8 | 2.8e8 | 9.8e9 | 8.7e9 | 31x | 54.6 | 56.0 |
| 4 | 8.6e8 | 4.5e8 | 9.8e9 | 9.2e9 | 20x | 95.6 | 53.4 |
| 10 | 2.0e9 | 6.0e8 | 9.7e9 | 9.2e9 | 15x | 167.4 | 52.5 |
| 20 | 4.0e9 | 6.9e8 | 9.5e9 | 8.9e9 | 13x | 287.3 | 53.1 |
| 30 | 5.8e9 | 7.1e8 | 9.9e9 | 8.9e9 | 12x | 410.7 | 55.8 |

(Cell counts are totals over 20 patterns; these rows are means over 20 reads.) Banded does 32 times fewer cell updates than Myers represents at k = 1, yet takes about the same time, and from k = 2 it is slower; at k = 30 it is 7.4 times slower. Myers processes 12 to 31 times more cells per second.

*What was measured.* Hardware counters (`perf`) are not available in this VM: perf is not installed and the virtual machine exposes no CPU performance-monitoring device, so instructions per cell and branch misses could not be read from hardware. Instead `scripts/instr_profile.py` ran both algorithms under `valgrind --tool=cachegrind --branch-sim=yes` on a generic x86-64 build (m = 100, 150 kbp of the text, 5 patterns; instruction counts are exact counts of executed instructions, branch mispredictions come from cachegrind's simulated predictor and are indicative only; the timings above use the `-march=native` build):

| k | algorithm | instructions per text character | instructions per cell | cells per text character | simulated mispredictions per text character |
|---|---|---|---|---|---|
| 5 | banded | 229.8 | 21.4 | 10.7 | 2.02 |
| 5 | Myers | 98.7 | 0.99 | 100 | 0.62 |
| 20 | banded | 754.3 | 18.9 | 39.9 | 1.95 |
| 20 | Myers | 98.8 | 0.99 | 100 | 0.64 |

Banded executes 2.3 times (k = 5) and 7.6 times (k = 20) the instructions of Myers per text character, while its measured time is 2.1 times (k = 5, m = 100 length sweep: 114 against 53 ms) and 5.4 times (k = 20: 287 against 53 ms) that of Myers. The time gap therefore goes in the same direction and is of similar size to the instruction gap, although the ratios are not equal (different builds, and instructions per cycle may differ). Per DP cell, banded costs about 19 to 21 instructions and Myers about 1. The simulated predictor reports more mispredictions per character for banded (about 2.0 against 0.6), but real misprediction penalties and cache behaviour were not measured, so branch misses are not established as a cause.

*Hypothesis (not measured).* Why a banded cell costs about 20 instructions while 64 Myers cells share about 50 (a block update, about 49 per block-character in the generic build) is explained here from the code, not from counters: (1) a block update packs 64 cells into a few dozen bitwise operations, while a banded cell needs three additions, two comparisons, a load and a store; (2) the banded column loop has a serial dependency (each cell needs the cell above it) that blocks vectorisation, whereas Myers carries that dependency inside the add of the carry chain; (3) the loop bound `lact` is data dependent and recomputed for every column, while Myers does the same work for every character; (4) for m = 100 the whole Myers state is two words, whereas banded reloads and stores an integer array. These are plausible but unverified. Fewer cells updated does not mean less time unless the per-cell cost is comparable. Edlib (which adds a banded cutoff to Myers' algorithm [9]) is about twice as fast as our unpruned Myers (24 to 32 ms against 52 ms at m = 100).

## 5. Conclusion

All exact methods return identical match sets: Sellers, banded, Myers, `filter` and `filter_fa` agree on every benchmark row (the benchmark aborts on any mismatch, with the filter run both with and without its scan fallback; the 100-read headline experiment is checked read by read) and on 300 randomised cross-checks (uniform, low-complexity, indel-heavy, clustered, N-containing; run with the fallback on and off). Edlib is checked differently, because it answers a different question: HW mode returns the minimum distance and only the ends at that minimum. It equals our `myers_best` (minimum-distance ends) in all 1,398 pattern checks, and for T = TTACGTTTTACGATT, P = ACGT, k = 1 it reports the single end {5} while the all-within-k set is {4, 5, 6, 11, 12, 13}. Edlib and WFA2 timings are therefore "best match per read", not all matches within k.

**What won where (*E. coli* 536, m = 100, time per query in ms).** Myers, `filter` and `filter_fa` are medians [10th, 90th percentile] over 100 reads; the other columns are means over 20 reads (Sellers 1, WFA2 5).

| k (k/m) | Sellers | Banded | Myers | Edlib (best/read) | WFA2 (best/read) | `filter` | `filter_fa` | one-hit / two-hit heuristic (read recall) |
|---|---|---|---|---|---|---|---|---|
| 2 (0.02) | 635 | 67 | 52.4 [49.9, 57.3] | 25 | 87 | **0.010** [0.006, 0.021] | 0.027 [0.024, 0.032] | 0.056 / 0.031 (1.00 / 1.00) |
| 8 (0.08) | 628 | 144 | 52.3 [50.1, 56.0] | 26 | 216 | **0.036** [0.019, 0.106] | 0.051 [0.041, 0.068] | 0.051 / 0.026 (1.00 / 0.95) |
| 10 (0.10) | 646 | 167 | 52.1 [50.2, 54.0] | 31 | 267 | 0.463 [0.292, 1.631] | **0.334** [0.209, 0.479] | 0.056 / 0.022 (1.00 / 0.95) |
| 14 (0.14) | 635 | 227 | 53.1 [50.7, 61.7] | 28 | 417 | 23.6 [16.9, 27.4] | **9.76** [7.19, 13.7] | 0.045 / 0.025 (1.00 / 0.70) |
| 20 (0.20) | 581 | 287 | 52.5 [50.4, 57.0] | 29 | 432 | 53.3 [51.3, 60.0] | 52.8 [51.1, 56.2] | 0.041 / 0.018 (0.95 / 0.40) |
| 30 (0.30) | 636 | 411 | 52.2 [50.1, 54.2] | 32 | 592 | 52.7 [50.5, 54.0] | 53.1 [51.5, 57.6] | 0.037 / 0.018 (0.30 / 0.00) |

- Bit-parallelism is the best exact scan: Myers is about 12 times faster than Sellers at m = 100 and flat in k; banded is slower than Myers from k = 2 despite updating far fewer cells (Section 4).
- With an index, the exact filters are 1,400 to 7,400 times faster than a Myers scan up to k/m = 0.08 (`filter` medians 0.007 to 0.036 ms against about 52 ms) and still about 110 times faster (`filter`) and 160 times (`filter_fa`) at k/m = 0.10. The **crossover with a full scan is at k/m of about 0.17** (measured 0.169 for `filter`, 0.191 for `filter_fa`; the cost model predicted 0.166 and 0.184). Beyond it the cost-model fallback keeps the filters at scan speed (median 53 ms at k = 20 and 30). Search time is nearly flat in n (0.007 ms at 250 kbp to 0.030 ms at 4.9 Mbp for `filter`, k = 5), against linear growth for the scans.
- **Where `filter_fa` wins and loses** (pure-filter medians, Section 4). It **loses for k <= 8 at m = 100** (k/m <= 0.08; 0.29x to 0.76x of the equal split's speed, that is 1.3 to 3.4 times slower) and for reads of 100 bp or less at 5% error, because seed selection (about 0.02 to 0.05 ms) costs more than the verification it saves. It **wins for k from 10 to 30 at m = 100** (k/m from 0.10 to 0.30; 1.35x to 2.7x faster; 1.4x at k = 10, 2.7x at k = 12, 1.5x at k = 15) and for **reads of 200 bp or more at 5% error** (k >= 10; 1.7x at m = 200, 3.4x at m = 400, 8.3x at m = 800). It has fewer hits than the equal split at every k tested. The crossover with Myers moves from k/m = 0.169 to 0.191, but above both the fallback decides and the two filters are equal to the scan.
- The index pays for itself after about 2 to 3 queries for k <= 13 and at most 5.5 for k <= 15; for the equal split it never pays at k >= 17 (Section 4).
- Edlib is 1.6 to 2.1 times as fast as our unpruned Myers (24 to 32 ms against about 52 ms), and WFA2 through Python is slow for scanning a genome (66 to 592 ms per query, growing with k, consistent with its O(ns) cost); both answer the best-match question only.
- **Results by dataset** (m = 100, k = 5 and 10): the same ordering holds on every genome. Per-query times for Myers / Edlib / `filter` / `filter_fa` at k = 5 are 52.6 / 37.9 / 0.026 / 0.040 ms on *E. coli* 536, 0.32 / 0.15 / 0.004 / 0.043 ms on SARS-CoV-2, 0.51 / 0.50 / 0.004 / 0.025 ms on phage lambda and 0.17 / 0.09 / 0.004 / 0.032 ms on human mtDNA (means over 20 reads). On the small genomes `filter_fa` is slower than `filter` at both k, because selection costs more than the tiny verification it saves, and the index build (about 20 ms) is dominated by clearing the 16.8 MB table. Two-hit read recall is 1.00 at k = 5 and 0.95 at k = 10 on all four genomes; one-hit recall is 1.00.

**Heuristic recall under stress** (*E. coli* 536 with about 1% injected microsatellite and homopolymer arrays, 20 reads per cell, so recall moves in steps of 0.05; read recall is the fraction of reads for which at least one true match is reported):

| scenario | k = 6 one-hit / two-hit | k = 10 one-hit / two-hit | k = 15 one-hit / two-hit |
|---|---|---|---|
| uniform | 1.00 / 1.00 | 1.00 / 0.95 | 1.00 / 0.45 |
| clustered | 1.00 / 1.00 | 1.00 / 1.00 | 1.00 / 1.00 |
| indel-heavy (90% indels) | 1.00 / 1.00 | 1.00 / 0.90 | 0.95 / 0.60 |
| low-complexity edge | 1.00 / 1.00 | 1.00 / 0.95 | 1.00 / 0.65 |
| low-complexity inside | 1.00 / 1.00 | 1.00 / 0.95 | 1.00 / 0.75 |

Clustered edits are the easiest case for a seed-based method, because they leave long clean stretches. The one-hit heuristic loses reads only under indel-heavy edits at k = 15 (0.95), and loses end positions inside low-complexity arrays (end recall 0.91 at k = 15, where one read can have thousands of matches). The two-hit trigger is much more fragile: recall falls to 0.45 to 0.75 at k = 15 in every scenario except clustered. Exact methods lose nothing in any scenario, and the filters get slower where hits pile up in repeats (`lowcomplex_inside`: 17,000 to 18,000 equal-split hits per read at k = 6 and 10, 0.6 to 1.6 ms). At the highest error rates the heuristic also collapses on uniform reads (k = 30: one-hit 0.30, two-hit 0.00). The heuristic's speed (0.02 to 0.08 ms) is similar to the filters' where the filters work, and it does not degrade at high error, which is exactly where its recall fails.

**Which method when.** Few queries or error above k/m of about 0.17: Myers, or Edlib if the best match per read is enough. Many queries against one genome with error up to about 0.15 and exact answers required: the index with `filter` (short reads, k <= 8) or `filter_fa` (k >= 10 or reads of 200 bp or more). A heuristic only when missing some matches is acceptable and the error is below about 0.15 (one-hit, recall 1.00 up to there on this data) or about 0.06 (two-hit, which starts losing reads at 0.08).

### 5.1 Limitations
- **Single machine, single thread.** The headline experiment (k sweep at m = 100) uses 100 reads per setting and reports medians with 10th to 90th percentile ranges; every other experiment averages 20 reads (Sellers 1, WFA2 5) in a single run, with no repetitions or confidence intervals, so differences under about 10% there should not be read as real. Timings are query time after the index is built; Edlib and WFA2 include Python call overhead (about 0.1 ms).
- **Simulated reads.** Reads carry exactly k random edits (uniform, 90% indels, or clustered), never real sequencer errors. Real long reads have context-dependent errors; no real reads were used.
- **Genome scale.** The index is a fixed 11-mer table. For *E. coli* (4.9 Mbp) an 11-mer occurs about once by chance; for a human-size genome it would occur hundreds of times, so the filters would lose selectivity much earlier. Longer seeds need a suffix array or FM-index (not implemented). Behaviour on large or highly repetitive genomes is untested.
- **Dataset naming.** The main genome is *E. coli* 536, not K-12 MG1655, which no reachable source provided. The three other genomes are small.
- **Cost model.** Constants are machine specific, fitted by NNLS on a synthetic genome, and the model assumes uniformly placed hits. It was validated on one other genome only (*E. coli* 536) and not on other CPUs. The fit is sensitive to the calibration range (Section 4): the first fit, restricted to k <= 15, made the fallback misfire near the crossover, and the range was widened afterwards, so the evaluation is not fully blind. The window and base-count terms are collinear, so their individual values are poorly determined.
- **Minimising hits is not minimising verified length.** `filter_fa` minimises the number of seed hits (Lemma 6), which is only a proxy for verification cost. In the `lowcomplex_edge` stress case at k = 10 it has 4.5 times fewer hits than the equal split (2,839 against 12,693) but 2.9 times more merged windows (1,536 against 531) and 2.6 times more verified bases (193,687 against 75,968), and is twice as slow (2.42 against 1.19 ms), because the equal split's hits cluster and merge while the rarer seeds are scattered. The DP does not model window merging.
- **Heuristic.** It models BLAST's seed-and-extend logic with hand-picked constants (match +2, mismatch -3, X-drop 20, HSP cutoff 34, two-hit window 40); it is not blastn. Recall numbers depend on those constants.
- **Edlib and WFA2** answer the best-match question and are not comparable in output to the all-within-k methods; WFA2 is used through `pywfa`, scanning a whole genome, which is not its design target.
- **Banded versus Myers.** Hardware counters were unavailable (no perf, no CPU performance-monitoring device in the VM). Instruction counts are exact but come from cachegrind on a generic build, and branch mispredictions come from a simulated predictor. The explanation of why a banded cell costs about 20 instructions and a Myers cell about 1 is a hypothesis reasoned from the code (Section 4), not a measurement.
- **Not implemented:** Landau-Vishkin diagonal transitions, block pruning for Myers, multithreading, FM-index search.

### 5.2 Novelty statement
This work is a comparative study of known algorithms plus a variant (`filter_fa`) that applies dynamic-programming seed selection, the idea of the Optimal Seed Solver [13], to an index of 11-mers with Myers verification and a calibrated fallback. The combination, the exactness and optimality-of-hit-count proofs written out in Section 3, the cost model, and the measurements are the contribution. It is not a new algorithm.

## 6. Gaps and Future Directions

1. **Optimal search schemes.** Equal or frequency-aware partitioning uses one split of the pattern into k+1 pieces. Search schemes over a bidirectional index split the pattern into more parts and allow errors in different places in different searches, which cuts the verified candidates further [14].
2. **FM-index and bidirectional search.** A 4^11 table limits seeds to 11 characters and scales poorly beyond bacterial genomes. A compressed full-text index (FM-index [16]) supports arbitrary-length seeds and inexact search by backtracking, as in BWA [17]; combining our DP seed selection with FM-index counts is the natural next step.
3. **Wavefront alignment.** WFA [10] runs in O(ns) time for score s and is exact for gap-affine alignment; here it was slow when scanning a genome (63 to 704 ms), but it is the stronger tool for verifying a candidate window, replacing Myers for long, similar sequences.
4. **Minimizers and chaining.** Sampling seeds as minimizers [15] and chaining them, as in minimap2 [11], reduces the index and the candidate count for long noisy reads, trading away our completeness guarantee.
5. **Pruned and SIMD/GPU bit-parallel verification.** Edlib [9], which adds a banded cutoff to Myers' algorithm, is about twice as fast as our unpruned Myers; adding such pruning, plus SIMD or GPU batching of many windows, would speed the verification step that dominates for k >= 10.
6. **Long noisy reads.** At 10 to 15% error the pigeonhole pieces are 6 to 9 characters long and every filter degenerates to a scan (our crossover is k/m of about 0.17 even on a small genome). Exact methods need either fewer-but-longer seeds with error-tolerant matching, or heuristics with a measured recall.
7. **Graphs and pangenomes.** Matching against labelled graphs has a quadratic conditional lower bound under SETH [12], so the same trade-off between exactness and speed returns in harder form.
8. **Minimise verified length or window count, not hit count.** Lemma 6 minimises seed hits, but verification cost depends on the merged window length and count, and the low-complexity stress case shows that fewer hits can mean more verified bases. A seed-selection DP whose objective is the (estimated) merged window length, the number of merged windows, or the cost-model prediction itself, for example by weighting each seed's occurrences by how clustered they are, is the natural fix.
9. **Better cost models.** Our model assumes uniform hits, was fitted on one machine, and its window and base terms are collinear (Section 4); modelling clustered hits (repeats), fitting with calibration over the crossover region and adapting seed length per query are open improvements.
10. **Limits of this setup:** see Section 5.1 (simulated reads, one run per point outside the headline experiment, bacterial genome scale, single thread, no hardware counters).

### References
Each citation was checked by web search (authors, year, venue, pages).

1. P. H. Sellers. The theory and computation of evolutionary distances: pattern recognition. *Journal of Algorithms* 1(4):359-373, 1980.
2. E. Ukkonen. Algorithms for approximate string matching. *Information and Control* 64(1-3):100-118, 1985.
3. G. Myers. A fast bit-vector algorithm for approximate string matching based on dynamic programming. *Journal of the ACM* 46(3):395-415, 1999.
4. S. Wu and U. Manber. Fast text searching allowing errors. *Communications of the ACM* 35(10):83-91, 1992.
5. G. Navarro. A guided tour to approximate string matching. *ACM Computing Surveys* 33(1):31-88, 2001.
6. S. F. Altschul, W. Gish, W. Miller, E. W. Myers and D. J. Lipman. Basic local alignment search tool. *Journal of Molecular Biology* 215(3):403-410, 1990.
7. A. Backurs and P. Indyk. Edit distance cannot be computed in strongly subquadratic time (unless SETH is false). *Proceedings of STOC 2015*, pp. 51-58.
8. G. M. Landau and U. Vishkin. Fast parallel and serial approximate string matching. *Journal of Algorithms* 10(2):157-169, 1989.
9. M. Sosic and M. Sikic. Edlib: a C/C++ library for fast, exact sequence alignment using edit distance. *Bioinformatics* 33(9):1394-1395, 2017.
10. S. Marco-Sola, J. C. Moure, M. Moreto and A. Espinosa. Fast gap-affine pairwise alignment using the wavefront algorithm. *Bioinformatics* 37(4):456-463, 2021.
11. H. Li. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34(18):3094-3100, 2018.
12. M. Equi, R. Grossi, V. Makinen and A. I. Tomescu. On the complexity of string matching for graphs. *Proceedings of ICALP 2019*, LIPIcs 132, pp. 55:1-55:15.
13. H. Xin, S. Nahar, R. Zhu, J. Emmons, G. Pekhimenko, C. Kingsford, C. Alkan and O. Mutlu. Optimal seed solver: optimizing seed selection in read mapping. *Bioinformatics* 32(11):1632-1642, 2016.
14. G. Kucherov, K. Salikhov and D. Tsur. Approximate string matching using a bidirectional index. *Theoretical Computer Science* 638:145-158, 2016.
15. M. Roberts, W. Hayes, B. R. Hunt, S. M. Mount and J. A. Yorke. Reducing storage requirements for biological sequence comparison. *Bioinformatics* 20(18):3363-3369, 2004.
16. P. Ferragina and G. Manzini. Opportunistic data structures with applications. *Proceedings of FOCS 2000*, pp. 390-398.
17. H. Li and R. Durbin. Fast and accurate short read alignment with Burrows-Wheeler transform. *Bioinformatics* 25(14):1754-1760, 2009.
18. S. B. Needleman and C. D. Wunsch. A general method applicable to the search for similarities in the amino acid sequence of two proteins. *Journal of Molecular Biology* 48(3):443-453, 1970.
19. T. F. Smith and M. S. Waterman. Identification of common molecular subsequences. *Journal of Molecular Biology* 147(1):195-197, 1981.
