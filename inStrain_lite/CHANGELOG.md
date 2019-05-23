# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project (attempts to) adhere to [Semantic Versioning](http://semver.org/).

## [0.3.0] - 2019-05-15
- Now calculates D', and normalized versions of D' and R2
- min_snp is now 
- Has the new gene_statistics.py and combine_samples.py scripts

## [0.2.8] - 2019-05-09
- fixed varbase = refbase error when there's a tie for most abundant variant
- made default output prefix be the fasta file prefix
- now generate scafold counts numpy array which show total variant counts for all positions by scaffold

## [0.2.7] - 2019-04-10
### Changed
- Change to how scaff2pair2info is handled 

## [0.2.6] - 2019-04-09
### Changed
- Minor change to the pickle protocol and stuff

## [0.2.5] - 2019-04-07
### Changed
- Allow changing filtering criteria, and some basic testing to make sure it works
- Will convert .sams to .bams and index for you

## [0.2.4] - 2019-03-29
### Changed
- Add a little bit of logging and efficiency into read report

## [0.2.3] - 2019-03-26
### Changed
- multiprocessing of scaffolds re-factored to be better
- now produces readable dataframe outputs
- make a nice read report
- add a --read_report option for filter_reads

## [0.2.2] - 2019-03-26
### Changed
- paired read filtering is now multi-processed

## [0.2.1] - 2019-03-26
### Changed
- argparse now displays the version number and defaults
- runs Matts way by default

## [0.2.0] - 2019-03-26
### Fixed
- implemented logging module
- made read filtering into a saved table
- implemented and tested "filter_reads.get_paired_reads"
- multiprocessing of SNP calling implemented
- setup.py is made; can now install itself
- changed the quoting on CC's table outputs to default
- changed "level" to "filter_cutoff"
- changed "min_coverage" to "min_cov"
- account for read overlap when calculating read length (half code is there, but reverted)
- changed to "no_filter" stepper
- clonality is now called the same with CC and Matt's versions
- min_cov is called with >= now

## [0.1.2] - 2019-03-12
- Stepper is now all

## [0.1.1] - 2019-03-12
- Tests now work

## [0.0.4] - 2018-06-25
- Adding some of the graph stuff and functions. The Python
is really slow, in the future it should call node2vec C++ library.

## [0.0.2] - 2018-06-21
### Fixed
- Test suite now runs on it's own and makes sure it succeeds

## [0.0.2] - 2018-06-21
### Fixed
- Test suite produces a command
- Program doesn't crash (lotsa small bug fixes)

## [0.0.1] - 2018-06-21
### Added
- Changelog and versioning is born
- Test suite is born

### Fixed
- Dumb print error is fixed
- Import statement is pretty :)
