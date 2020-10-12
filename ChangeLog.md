# Change log for the PCRamp program

- Version 0.3 (August 25, 2020)
	- Fixed a bug in random traget sequence selection step in PCR::random_assay().
	- Fixed formatting error in JSON output file ( missed a '\n' in write_json() ).
	- Added additional, synonymous command line arguments `--T.prefix` (a synonym for `--target.prefix`) and `--B.prefix` (a synonym for `--background.prefix`) to try to make the command line arguments internally consistent.
	- Modified confusing command line arguments to improve consistency
	- Removed the option of JSON-formatted input configuration file. This option can be added back, but this method of specifying the input configuration is currently not well-tested (and would need more detailed documentation to be useful).

- Version 0.2 (April 23, 2020)
	- Improved the multiplex compatibility algorithm to better exclude overlapping amplicons.
