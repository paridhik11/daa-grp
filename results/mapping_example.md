### Worked example: `scripts/map_reads` on E. coli 536

A 60 bp stretch of E. coli 536 (positions 2000000..2000060) was copied and given 2 substitutions, 1 insertion and 1 deletion (4 edits), then mapped with k = 8:

```
scripts/map_reads data/raw/ecoli536.fa reads.fa -k 8
```
TSV output:
```
read	length	status	distance	start	end	strand	cigar	best_ends
example_read	60	mapped	4	2000000	2000060	+	10=1X15=1I14=1D12=1X6=	10
unmappable_read	60	unmapped	*	*	*	*	*	0
```
Alignment reconstructed from the CIGAR (reference on top, read below; `x` mismatch, `-` gap):
```
ATATGGCAAAAGCGCTCAGGGCGGGA-TCATCAACATCGTCACCCAGCAGCCGGACAGCAC
||||||||||x||||||||||||||| |||||||||||||| ||||||||||||x||||||
ATATGGCAAACGCGCTCAGGGCGGGAGTCATCAACATCGTC-CCCAGCAGCCGGTCAGCAC
```
The reported distance is 4, equal to the number of X, I and D operations in `10=1X15=1I14=1D12=1X6=`, and the start (2000000) matches the planted position (2000000). The second read has no alignment within k = 8 and is reported as unmapped. SAM-like output (`--format sam`):
```
@HD	VN:1.6	SO:unsorted
@SQ	SN:gi|110640213|ref|NC_008253.1|	LN:4938920
@PG	ID:map_reads	PN:map_reads
example_read	0	gi|110640213|ref|NC_008253.1|	2000001	255	10=1X15=1I14=1D12=1X6=	*	0	0	ATATGGCAAACGCGCTCAGGGCGGGAGTCATCAACATCGTCCCCAGCAGCCGGTCAGCAC	*	NM:i:4
unmappable_read	4	*	0	0	*	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	*
```
