# cuteSV memory failure and safe restart

`Cannot allocate memory` during `Rebuilding signatures of structural variants`
means cuteSV could not allocate memory while reconstructing candidate/read lists.
cuteSV 2.1.3 loads whole signature types in worker processes; each worker can hold
a large list. The ONTSeq default is therefore one cuteSV worker, while other tools
keep their independent `--threads` setting. The installed Desktop service inherits
this default after restart. `--cutesv-threads N` explicitly overrides it for
`run`, `analyze`, `serve`, and `watch`.

One worker reduces concurrency, but does not bound the largest list. More available
RAM/swap may still be required. Do not lower analytical thresholds, omit reads or
disable the second caller merely to make a failed analysis appear completed.

## WSL scratch storage

Reducing the worker count alone may not fix `pickle.load` ENOMEM. An isolated
one-worker diagnostic still failed while swap remained available. A controlled
synthetic file test then failed on a Windows-mounted directory after 19.9 million
records, while the same file loaded all 28 million records on native Linux storage.
This supports moving scratch I/O off the Windows mount; it does not establish the
precise kernel cause or guarantee arbitrary input sizes fit in memory.

cuteSV signature/pickle intermediates now use Python's system temporary directory
(`/tmp` in the qualified WSL environment). Keep `TMPDIR` on native Linux storage
with sufficient free disk space. The final VCF is still staged beside its destination,
validated, then atomically promoted. Both temporary directories are removed after
success, tool failure, invalid output or a runner exception. The effective scratch
root is recorded in tool provenance. No analytical threshold is changed.

## Known float-boundary defect

Changing worker count changes the original cuteSV partition boundaries. A synthetic
comparison exposed dropped reads at fractional boundaries. ONTSeq corrects this
for one exact pinned 2.1.3 script body by executing a temporary copy with integer
chunk sizes (`max(1, floor(...))`). Original installation files remain unchanged.
The source and execution hashes and the compatibility build ID are recorded in
provenance; the source identity also participates in the SV resume signature.
Unknown script bodies are passed through without transformation and are explicitly
identified as unqualified. This is a local compatibility correction, not an official
upstream release. No support threshold or read-selection rule is relaxed.

## Existing runs

- Inspect `provenance/run.json` for the whole stage verdict. A caller file alone
  does not prove that the configured SV workflow or fusion assessment completed.
- Keep the original report, normalized data and release checksums together. This
  correction changes future assembly; it does not rewrite archived reports.
- On exactly the same software version and Git commit, completed stages can resume
  only after their signatures and output checksums verify. Failed stages execute
  again. Changing only `--cutesv-threads` affects the SV plan, not methylation/QC.
- Sniffles2 and cuteSV share one SV stage. Restarting that stage reruns both callers;
  an isolated Sniffles2 JSON is not a separately validated resume checkpoint.
- Installing a new code revision invalidates all stage signatures by design.
  Never claim an old Git commit to bypass that check. Use a new run directory for
  a full rerun after an update, or perform an isolated cuteSV diagnostic in a
  separate directory. Such a diagnostic is not a completed ONTSeq result and must
  not be spliced into an old checksummed release.

The synthetic fixture in `tests/test_cutesv_real_tool.py` can be enabled with
`ONTSEQ_CUTESV_REAL_TOOL=1`. It requires cuteSV 2.1.3 and pysam and compares the
normalized deletion from one and four workers. It makes no assertion about memory
requirements on biological samples.
