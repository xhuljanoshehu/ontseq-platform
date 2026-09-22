# Qualified independent-MM-group modkit build

ONTSeq can recognize an explicitly qualified executable for the independent-group
scanner correction in [modkit PR #709](https://github.com/nanoporetech/modkit/pull/709).
The reviewed source is frozen at
[`9fb9aea763aa1b78ac1172fa6d6734e7955ccc70`](https://github.com/SuhasSrinivasan/modkit/tree/9fb9aea763aa1b78ac1172fa6d6734e7955ccc70).
Its version banner remains `modkit 0.6.4`. A banner alone never enables the exception.

## Build and qualification

`build-support/modkit-pr709/source.json` records source and build constraints;
`build-support/modkit-pr709/Cargo.lock` freezes the dependency graph and crate checksums.
The lock retains `hts-sys 2.2.0`, as specified by
[upstream PR #666](https://github.com/nanoporetech/modkit/pull/666): 2.2.1 causes ten
binding/type compilation errors with rust-htslib 0.46.0. The analytical Rust source is
unchanged.

Use Rust 1.93.1, a Linux C/C++ toolchain and the native development dependencies required
by rust-htslib/OpenSSL. ONTSeq's Python test dependencies, pysam and samtools must be
available. The local qualification uses the installed conda Linux compiler/sysroot so
that the resulting binary can run in the existing WSL runtime.

```sh
scripts/build_modkit_pr709.sh /absolute/build-directory
```

The script builds a candidate from the exact upstream source, runs the real executable
against synthetic raw-pileup and full-adapter tests, preserves the upstream licence and
writes checksums. It does not install the candidate or add it to the runtime allowlist.
A new compiler, dependency graph or executable checksum requires review and qualification;
only explicitly registered checksums in `modkit_build.py` can bypass the stock 0.6.4 guard.
The build directory, Cargo.lock, source revision, compiler information and test results
must accompany an installation. Retain the upstream
[ONT Public License 1.0](https://github.com/SuhasSrinivasan/modkit/blob/9fb9aea763aa1b78ac1172fa6d6734e7955ccc70/LICENCE.txt)
and notices with distributed source/binaries. This does not assign a licence to ONTSeq.

## Runtime selection and rollback

Install a qualified executable separately from conda's stock `modkit`. Point Desktop's
existing `modkitExecutableWsl` setting to its absolute WSL executable path, or use the
canonical CLI's `--modkit` argument. Keep the old setting, wheel and executable available
for rollback. The normal conda environment still supplies stock 0.6.4; installing that
environment alone does not install or qualify PR #709.

The backend records SHA-256, build ID and upstream commit in both stage planning and
normalized tool provenance. It checks identity before and after processing and rejects
replacements. Resume therefore cannot reuse results produced by a different executable
with the same version banner. Existing failed outputs must be recomputed, never relabelled.

## Scope of evidence

Literal synthetic oracles cover independent/combined groups, forward/reverse strands,
clipping, confidence filtering, measured zero, low coverage and absent modifications.
Separate adapter tests confirm normalization and deliberate refusal of mixed MM/noMM
inputs when modkit reports failed records. Such failures remain fatal even after PR #709;
no untagged reads are silently converted to canonical bases or removed from the BAM.
No clinical sensitivity, specificity, diagnostic accuracy or patient result is established
by this technical qualification. See `CLINICAL_VALIDATION.md` for validation impact.

## Registered local executable

The qualified Linux executable SHA-256 is
`f71542ee41ec964202a751919c3f47e0d54398ba0e5c5c1ccb0990c8a65769aa`. Its 76-case raw/adapter qualification passes; stock 0.6.4 fails 19 of the 73 correction-matrix cases. The other three cases are the existing real-tool integration checks. The qualified candidate does not establish a successful rerun of an existing user sample.

The registered local build requires glibc 2.38 or later and the existing conda runtime
`libgcc_s.so.1`. Its embedded runtime search path is specific to that installation; it
is not advertised as a generic portable Linux binary. The Desktop executable already
supports selecting this separate tool.

The exact upstream workspace suite passes 185 tests with its 14 declared ignored tests.
Run it serially because upstream tests share fixed temporary output filenames. In the
release profile, enable `profile.release.package.mod_kit.debug-assertions=true` for the
upstream `should_panic` test that checks a debug assertion. These test execution settings
avoid harness-only failures without editing upstream source or skipping active tests.
The installed candidate remains the separately checksummed ordinary release build.
