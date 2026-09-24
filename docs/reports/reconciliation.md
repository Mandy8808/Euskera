# Repository reconciliation

## Source inventory

The initial local main and fetched origin/main both pointed to `5c5243e`.
The `6dd4` worktree contained uncommitted recoverable-consolidation work based
on `8c676fb`; it was not yet in GitHub main. The `6d5b` worktree contained no
additional tracked code, test, documentation, or CI changes. Both worktrees
were preserved rather than reset or deleted.

## Integrated changes

- Imported `euskera/io/consolidation.py` and merged the existing consolidation
  writer, configuration, parameter-validation, and failure-recovery tests.
- Connected `cleanup_policy` and `consolidation_batch_size` to legacy writers
  and every selective-output writer, including grid and diagnostics.
- Allowed serialized `rules` through the legacy update interface; previously
  even `rules=None` from `OutputConfig.to_dict()` was rejected.
- Added end-to-end selected-output recovery tests for both formats and policies.
- Integrated English consolidation documentation and added cleanup settings to
  the selective-output notebook. Historical implementation reports describe
  earlier snapshots; this report records the reconciled implementation.
- Fixed two relative documentation links in `profiles/README.md`. These broken
  links reproduce a failing test locally and affect every Python CI job.
- Added an explicit Git LFS download for workflow notebooks before execution in
  CI. Other large notebooks and simulation results are not downloaded by this step.

## Validation

- Full suite in `ciencia` with FFTW: 185 passed before adding the final four
  integration tests; the expanded output module then passed all 12 tests
  (189 distinct tests covered in total).
- The four workflow notebooks were executed after integration.
- Incremental-cleanup warnings are expected and describe its reduced protection
  against corruption after original files have been deleted.
- Local validation does not substitute for the subsequent GitHub Actions run.

## Remaining work

Single-sample reading and numerical evolution checkpoints remain separate future
features. Recoverable consolidation is now integrated; it should not be
reimplemented. No simulation outputs were regenerated in the project directory.
