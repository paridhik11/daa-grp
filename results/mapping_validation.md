## Alignment output validation

3200 random reads in 320 runs (texts of 200 to 20,000 bp, some with N; k in 0..12; uniform, indel-heavy and clustered edits; 2298 mapped, 902 unmapped): 0 failures. Checks per mapped read: distance equals the Sellers minimum (Sellers at k = d matches, at k = d - 1 does not); CIGAR edit count (X + I + D) equals the reported distance; the CIGAR replays correctly against read and text; the end is a Sellers end position. Unmapped reads have Sellers minimum distance above k.

