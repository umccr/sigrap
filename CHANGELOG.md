# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- New `sig_plot_rainfall()` function for generating rainfall plots showing intermutational distances
- New `sig_plot_strand_bias()` function for transcriptional and replicative strand bias analysis
- Enhanced CLI support for MutationalPatterns with `--rainfall` and `--strand-bias` options
- Improved plot dimensions (4x larger) for better quality publication-ready figures
- Development Docker workflow for testing and deployment
- Enhanced reference genome handling with better error messages

### Changed

- Updated `sig_plot_snv()` function to include optional rainfall plot generation
- Enhanced `sig_workflow_run()` function with new `rainfall` and `strand_bias` parameters
- Improved plot saving with dimension-specific sizing for different plot types
- Better error handling for missing TxDb packages in strand bias analysis

### Technical Details

- Added dependency on `TxDb.Hsapiens.UCSC.hg38.knownGene` for strand bias analysis (optional)
- Enhanced CLI argument parsing with proper defaults for new options
- Improved test coverage with comprehensive Docker integration tests

### Dependencies

- Added optional dependency: `TxDb.Hsapiens.UCSC.hg38.knownGene` (for strand bias analysis)

## [0.2.0] - 2022-06-10

### Added

- Docker support
- conda-lock support
- Modularised CLI

### Changed

- Updated to `{signature.tools.lib}` v2.1.2

## [0.1.0] - 2022-06-07

### Changed

- Updated to `{signature.tools.lib}` v2.1.2
- Refactored problematic code feeding into `tidyr::pivot_`

## [0.0.4] - 2022-06-07

### Changed

- Updated to `{signature.tools.lib}` v2.1.2
- Refactored problematic code feeding into `tidyr::pivot_`

## [0.0.3] - 2022-05-19

### Fixed

- Fixed bug when running HRDetect ([issue3](https://github.com/umccr/sigrap/issues/3), [pr4](https://github.com/umccr/sigrap/pull/4))

## [0.0.1] - 2021-12-22

### Added

- Initial release
- Signature tool wrappers (CHORD, HRDetect, and MutationalPatterns) moved from gpgr
