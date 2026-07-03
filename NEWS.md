# sigrap 0.3.0 (2026-05-22)

- :star: Add MBS (Multi-Base Substitution) support: `sig_count_mbs`, `sig_plot_mbs`, `sig_mbs_table` ([pr14](https://github.com/umccr/sigrap/pull/14)).
- :star: Add `predefined_dbs_mbs` parameter to handle pre-merged DBS/MBS variants in the VCF ([pr14](https://github.com/umccr/sigrap/pull/14)).
- :bug: Fix crash when sample has zero DBS variants -- return placeholder plots and skip signature fitting ([pr14](https://github.com/umccr/sigrap/pull/14)).
- :bug: Fix `RelFreq` producing `NaN` in JSON output for zero-variant samples ([pr14](https://github.com/umccr/sigrap/pull/14)).

# sigrap 0.1.0 (2022-06-10)

- :star: Add Docker support ([pr6](https://github.com/umccr/sigrap/pull/6)).
- :star: Add conda-lock support ([pr6](https://github.com/umccr/sigrap/pull/6)).
- :star: Modularise CLI

# sigrap 0.0.4 (2022-06-07)

- :star: Update to `{signature.tools.lib}` v2.1.2
  - also refactored problematic code feeding into `tidyr::pivot_`.

# sigrap 0.0.3 (2022-05-19)

- :bug: Fix bug when running HRDetect
  ([issue3](https://github.com/umccr/sigrap/issues/3),
  [pr4](https://github.com/umccr/sigrap/pull/4)).

# sigrap 0.0.1 (2021-12-22)

- Initial release. Moving signature tool wrappers (CHORD, HRDetect,
  and MutationalPatterns) from gpgr (https://github.com/umccr/gpgr).
