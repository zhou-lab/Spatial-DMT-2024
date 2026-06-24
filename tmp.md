| Metric | 16D, our pipeline | E65, our pipeline | E65, their reproduced pipeline |
|---|---:|---:|---:|
| Raw read pairs | 30.18M | 333.28M | 333.28M |
| Structure / detected | 86.7% | 90.4% detected, 72.3% passed | 76.0% structure-pass |
| Mapped | 48.65M primary / 50.15M primary reads | 158.76M primary / 241.02M reads | 52.05M before filter, 48.11M after filter |
| Mapping rate | **97.0% primary** | **65.9% primary** | nBismark unique-only |
| MAPQ > 0 | high; 36.77M of 49.57M mapped records | 118.53M of 158.76M primary mapped | unique Bismark only |
| Duplication | ~8.7% by flagstat primary duplicate | 18.2% | not comparable |
| CpG methylation | **0.681** | 0.553 | 0.657 after filter |
| Non-CpG / CH | **~1.0%** CA/CC/CT | **5.4–6.9%** CA/CC/CT | 5.89% before filter, **1.73% after filter** |
| Lambda retention | **0.037** | 0.965 (high -> better) | / |
