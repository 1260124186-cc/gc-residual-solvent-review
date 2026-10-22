# GC sequence review

The `gc_review` package turns a GC-FID run sequence into a QA review report:

```bash
python -m gc_review.cli data/sequence.csv data/chromatograms.csv report.json
```

The files under `data/` are a reproducible simulated sequence for checking the
command and the expected report shape. They are not laboratory records.
