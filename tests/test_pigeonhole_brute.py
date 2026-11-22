#!/usr/bin/env python3
"""Exhaustive check of the combinatorial core of Theorems 4 and 12 (no algorithm code is involved).
For small m, k, s: enumerate EVERY edit script with at most k operations (match / substitution / deletion per pattern character, any number of
insertions in each of the m + 1 gaps, including the gaps between pieces and at both ends) and EVERY split of the pattern into J = k + s
contiguous pieces.  A piece is touched when an operation substitutes or deletes one of its characters or inserts inside it (between two of its
characters); an insertion on a piece boundary touches neither neighbour.  Checks:
  (a) at least s pieces are untouched                      (Theorem 12; s = 1 is Theorem 4)
  (b) the text diagonals of the untouched pieces lie within a band of width k
  (c) the aligned substring [start, end) lies inside [dmin - k, dmax + m + k)   (the verification window of Theorem 12)
  (d) the bound s is tight: some script leaves exactly s untouched, so requiring s + 1 hits would lose real matches.
Prints the number of (split, script) pairs checked; exits non-zero on any violation."""
import itertools
import sys


def compositions(m, parts):
    for cuts in itertools.combinations(range(1, m), parts - 1):
        b = (0,) + cuts + (m,)
        yield [(b[i], b[i + 1]) for i in range(parts)]


def scripts(m, k):
    """All (op[i] in 0 match,1 sub,2 del; ins[g] for g in 0..m) with at most k edits."""
    slots = [("o", i) for i in range(m)] + [("g", g) for g in range(m + 1)]

    def rec(idx, left, op, ins):
        if idx == len(slots):
            yield tuple(op), tuple(ins)
            return
        kind, i = slots[idx]
        if kind == "o":
            for v in (0, 1, 2):
                c = 0 if v == 0 else 1
                if c <= left:
                    op[i] = v
                    yield from rec(idx + 1, left - c, op, ins)
            op[i] = 0
        else:
            for c in range(left + 1):
                ins[i] = c
                yield from rec(idx + 1, left - c, op, ins)
            ins[i] = 0
    yield from rec(0, k, [0] * m, [0] * (m + 1))


def check(m, k, s):
    J = k + s
    L = k + 3   # left flank: no clamping at the text start
    pairs, tight = 0, False
    splits = list(compositions(m, J))
    for op, ins in scripts(m, k):
        # diag[a] = L + (insertions in gaps <= a) - (deleted characters before a)
        diag = []
        ci, cd = 0, 0
        for a in range(m + 1):
            ci += ins[a]
            diag.append(L + ci - cd)
            if a < m and op[a] == 2:
                cd += 1
        start, end = L, L + m + sum(ins) - sum(1 for o in op if o == 2)
        for pcs in splits:
            un = [a for a, b in pcs
                  if all(op[i] == 0 for i in range(a, b)) and all(ins[g] == 0 for g in range(a + 1, b))]
            pairs += 1
            if len(un) < s:
                return f"(a) violated m={m} k={k} s={s} op={op} ins={ins} split={pcs}"
            ds = [diag[a] for a in un]
            if max(ds) - min(ds) > k:
                return f"(b) violated m={m} k={k} s={s} op={op} ins={ins} split={pcs}"
            if not (start >= min(ds) - k and end <= max(ds) + m + k):
                return f"(c) violated m={m} k={k} s={s} op={op} ins={ins} split={pcs}"
            if len(un) == s:
                tight = True
    if not tight:
        return f"(d) not tight m={m} k={k} s={s}"
    return pairs


def main():
    total, bad = 0, 0
    for k in (1, 2, 3):
        for s in (1, 2, 3):
            for m in range(k + s, k + s + 4):
                if m > 9:
                    continue
                r = check(m, k, s)
                if isinstance(r, str):
                    print("FAIL", r)
                    bad += 1
                else:
                    total += r
    print(f"{total} (split, edit script) pairs checked (m <= 9, k <= 3, s <= 3), {bad} failures; the bound of s untouched pieces is attained in every setting")
    return bad


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
