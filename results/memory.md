## Peak memory per method (measured)

E. coli 536 (n = 4,938,920), m = 100, k = 10, 5 reads; peak resident set size (VmHWM). working memory = peak RSS minus the steady RSS after the text and reads are loaded and the heap trimmed (8.4 MB); methods whose peak does not exceed the transient peak of loading the file (about 20 MB, the file is read line by line into a string) are shown as negligible. Edlib and WFA2 run in a Python process, so their baseline differs.

| method | peak RSS (MB) | working memory of the method (MB) | expected |
|---|---|---|---|
| noop | 20.2 | negligible (below the loading peak) | text only (n bytes) plus read buffers |
| sellers | 20.2 | negligible (below the loading peak) | O(m) |
| banded | 20.1 | negligible (below the loading peak) | O(m) |
| myers | 20.2 | negligible (below the loading peak) | O(B sigma) |
| myers_pruned | 20.1 | negligible (below the loading peak) | O(B sigma) |
| lv | 50.8 | 42.4 | 2 arrays of n + 2k ints (8n bytes) |
| filter | 63.9 | 55.5 | 4(4^11+1) + 4n bytes index (about 37 MB) |
| filter_fa | 63.9 | 55.5 | same index |
| filter_fa_w | 63.9 | 55.5 | same index |
| blast1 | 63.9 | 55.5 | same index |
| blast2 | 63.9 | 55.5 | same index |
| map | 63.9 | 55.5 | same index + O(m (m + k)) per cluster |
| edlib (Python process) | 39.3 | negligible (below the loading peak) | library; measured as ru_maxrss increase over loading the text |
