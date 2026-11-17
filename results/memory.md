## Peak memory per method (measured)

E. coli 536 (n = 4,938,920), m = 100, k = 10, 5 reads; peak resident set size (VmHWM). working memory = peak RSS minus the steady RSS after the text and reads are loaded and the heap trimmed (8.4 MB); methods whose peak does not exceed the transient peak of loading the file (about 20 MB, the file is read line by line into a string) are shown as negligible. Edlib and WFA2 run in a Python process, so their baseline differs.

| method | peak RSS (MB) | working memory of the method (MB) | expected |
|---|---|---|---|
| noop | 20.2 | negligible (below the loading peak) | text only (n bytes) plus read buffers |
| sellers | 20.2 | negligible (below the loading peak) | O(m) |
| banded | 20.2 | negligible (below the loading peak) | O(m) |
| myers | 20.2 | negligible (below the loading peak) | O(B sigma) |
| myers_pruned | 20.2 | negligible (below the loading peak) | O(B sigma) |
| lv | 20.2 | negligible (below the loading peak) | packed text n/4 bytes (1.2 MB) + L1-sized chunk buffers |
| lv_direct | 50.8 | 42.4 | 2 layers of n + 2k ints plus n flags (about 44 MB) |
| filter | 69.3 | 60.9 | index 4(4^11+1) + 4n bytes (37 MB) plus build temporaries (counter copy 17 MB, flags 5 MB) |
| filter_auto | 69.3 | 60.9 | same index |
| filter_fa | 69.3 | 60.9 | same index |
| filter_fa_w | 69.3 | 60.9 | same index |
| blast1 | 69.3 | 60.9 | same index |
| blast2 | 69.3 | 60.9 | same index |
| map | 69.3 | 60.9 | same index + O(m (m + k)) per cluster |
| edlib (Python process) | 38.9 | negligible (below the loading peak) | library; measured as ru_maxrss increase over loading the text |
