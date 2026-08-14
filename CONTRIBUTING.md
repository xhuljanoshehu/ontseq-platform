# Contributing

1. Open or link an issue describing the scientific or engineering change.
2. Work on a focused branch and keep data outside the repository.
3. Add or update tests for every contract, rule, filter or report change.
4. Record tool versions, parameter changes and reference-build assumptions.
5. Run `make safety`, `make test` and `make lint` before requesting review.
6. Require both code review and domain review for changes that alter biological output.

Do not import third-party or institutional source code until its license and ownership are
documented. A pull request is not a clinical release. Clinical release criteria are defined
separately in `docs/CLINICAL_VALIDATION.md`.
