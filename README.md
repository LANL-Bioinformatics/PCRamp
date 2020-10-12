# PCRamp
Software for designing multiplex-compatible, PCR-based enrichment assays

## Overview

The PCRamp program designs mulitplex PCR assays for amplicon sequencing-based target enrichment. 

Given one or more target nucleic acid sequences, PCRamp iteratively designs the requested number of multiplex-compatible PCR primers to amplify non-overlapping regions of the target sequences. At each iteration, the design algorithm selects a PCR primer pair that (a) satisfies all specified design constraints (based on melting temperature, hairpin formation, G+C content, length, etc.) and (b) provides the best enrichment by amplifing the largest number of target sequences. This strategy preferentially selects the most conserved regions of the target sequences for amplification.

Optionally, users can also specify one or more *background* nucleic acid sequences, which must *not* be amplified by a pair of PCR primers.

The output of the PCRamp program is a list of PCR primer pairs. Associated with each primer pair is a list of target sequences that are predicted to be amplified by this primer pair.

## Design strategy - how does PCRamp work?

The PCRamp program attempts to find highly conserved, target-specific, multiplex compatible PCR primers. There are many target enrichment scenarios that the PCRamp program seeks to address, including:
- **Gene target enrichment**: Challenges include potentially large numbers (tens of thousands) of moderately diverse target sequences. For many genes, the target sequences are short (less than 5 kb) and there may be closely related nearest neighbors.
- **Viral target enrichment**: Challenges include potentially large numbers (tens of thousands) of highly diverse target sequences. However, for many viruses, the target genomes sequences are fairly short (less than 50 kb) and the nearest neighbors are usually *not* closely related (which often means we do *not* need to provide background sequences to ensure PCR specificity). Degenerate nucleotides may be needed to cope with high sequence diversity.
- **Bacterial target enrichment**: Challenges include fairly large genome sequence lengths (around 5 mb) and often closely related near-neighbors. However, most bacterial targets have relatively few available genome sequences (typically less than five thousand).

The following algorithmic features are intended to address these challenges within a single assay design program:
- A random-sampling approach to designing assays that avoids the need to load all of the target and background genomes into main memory at the same time. The target coverage (i.e. inclusivity; number of targets amplified) and specificity (i.e. avoiding amplification of background sequences) of the randomly selected PCR assays are then improved with optional local optimization steps. The downside to this approach is that assay designs are no longer deterministic - it is expected that a different set of PCR primers for the same set of targets and background sequences will be generated every time PCRamp is run.
- A combination of heuristic and physics-based rules for predicting successful PCR.
- Where possible, previously designed PCR primers are "reused" to increase the number of diverse target sequences that can be enriched while minimizing the total number of primer oligos required for the multiplex PCR.
- Efficient, highly parallel implementation that uses multiple processors (via MPI), multiple cores (via OpenMP) and data parallel (SIMD) instructions for x86 CPUs.

## Building and installing PCRamp

PCRamp is written in C++ and requires a C++ compiler and a local installation of MPI ("Message Passing Interface") to compile. PCRamp was developed and tested using OpenMPI with both the gnu (on Linux) and clang (on OS X) C++ compilers. The included `Makefile` is intentionally very simple (and hopefully easy to read). With the exception of zlib (included on most Unix-like systems; needed for the on-the-fly reading of compressed fasta files), there are *no* other software dependencies. After downloading the PCRamp source code files, running the `make` command should build the `pcramp` program. There is no formal install process, but the `pcramp` program can be manually copied to any desired location. Please keep in mind that when running on a cluster computer with MPI, the version of MPI used to compile PCRamp must be the same as the `mpirun` program.

## Running PCRamp

PCRamp uses MPI to run on a cluster computer and will automatically use all of the available CPU cores on each compute node. As a result, MPI and any cluster scheduling software (i.e. slurm, UGE, SGE, torque, ...) should be configured to run a *single* instance of PCRamp on each computer in the cluster and let each instance use all of the available CPU cores.

When running on multiple computers, PCRamp is typically envoked via a batch script (submitted to a cluster scheduler) or the command line (when running in interactive mode) using some variant of `mpirun <MPI options> pcramp <PCRamp options>`. The particular MPI options will depend on the configuration of each particular cluster computer.

PCRamp can also run on a single computer (e.g. workstation, laptop). In this case, simply omit the `mpirun` and directly invoke `pcramp` from a script or the command line. By default, `pcramp` will use all available threads (unless the `--thread` option is provided to limit the number of threads).

## Input sequence formats

Target sequences (for which to be enriched) and optional background sequences (that are not to be amplified) are provided to PCRamp in the fasta file format. These fasta files may optionally be compressed using the gzip program (and they will automatically be decompressed "on-the-fly" to keep filesystem space requirements to a minimum). All fasta files must have one of the following file extensions: `.fna`, `.fasta` or `.fa` (with an optional `.gz` suffix to indication compression).

PCRamp accommodates two of the common ways that biological sequences are stored in filesystems, each with its own command line flag:
1. For relatively short sequences (like single-segment viruses and genes), sets of distinct target sequences are often stored in a *single* fasta file. A toy example of this arrangement is:
```
>virus1
ACTAGCGATGCGACGTAGCTAGCAGCGATGCAGCTAGCAGTCGTA
>virus2
ACTAGCGCATGCGACGTAGCTAGCAGCGAGCAGCTAGACAGTCGTAGCTA
>virus3
ATGCGATCATAGCGATGCGACGTAGTAGCAGCGATGAGCTAGCAGTCGTA
```
To read multiple, separate targets from a single fasta file, specify the fasta file on the command line with the *lower case* `-t` flag. For example, `pcramp -t <fasta file of targets> ...`. If the target sequences are contained in multiple fasta files, then the `-t` flag may be repeated on the command line, `pcramp -t <fasta file of targets 1> -t <fasta file of targets 2>`. Similarly, if distinct background sequences are stored in a single fasta file, they would be specified using `pcramp -b <fasta file of backgrounds>`. As with the `-t` flag, the `-b` flag may also be repeated to specify additional background sequences.

2. For longer sequences, as well as multi-segment viruses and multi-chromosomal bacteria, it is common to store one or more sequences/segments/chromosomes from an individual target in *different* files. All of the sequences found in a same directory as assumed to belong to the same target. For example, a directory listing of different *Bacillus* genomes might look like:
```
Bacillus/
	anthracis/
		Bacillus_anthracis_str._Ames_7845/
			NC_003997.3.fna.gz
		Bacillus_anthracis_str.__Ames_Ancestor__8445/
			NC_007322.2.fna.gz
			NC_007323.3.fna.gz
			NC_007530.2.fna.gz
		Bacillus_anthracis_str._Australia_94_167335/
			wgs.AAES.1.fna.gz
	cereus/
		Bacillus_cereus_G9241_167215/
			wgs.AAEK.1.fna.gz
		Bacillus_cereus_MC118_399245/
			wgs.AHEM.1.fna.gz
	thuringiensis/
		Bacillus_thuringiensis_serovar_tochigiensis_BGSC_4Y1_161475/
		 	wgs.ACMY.1.fna.gz
		Bacillus_thuringiensis_serovar_wuhanensis_2147375/
		 	wgs.NFEE.1.fna.gz
		Bacillus_thuringiensis_str._Al_Hakam_15065/
			NC_008598.1.fna.gz
			NC_008600.1.fna.gz
```

To read multiple, separate targets, where each target is a collection of one or more fasta files in a single directory, specify the directories (or parent directories) on the command line with the *capital* `-T` flag. For example, `pcramp -T Bacillus/anthracis/Bacillus_anthracis_str.__Ames_Ancestor__8445` will load all of the data in fasta files contained in the `Bacillus/anthracis/Bacillus_anthracis_str.__Ames_Ancestor__8445/` directory (i.e. the chromosome and plasmid sequences stored in `NC_007322.2.fna.gz`, `NC_007323.3.fna.gz` and `NC_007530.2.fna.gz`).

As another example, the command `pcramp -T Bacillus/anthracis` will *recursively* search the `Bacillus/anthracis` directory to load the three *B. anthracis* genomes `Bacillus_anthracis_str._Ames_7845`, `Bacillus_anthracis_str.__Ames_Ancestor__8445` and `Bacillus_anthracis_str._Australia_94_167335`. The sequences in each subdirectory will be associated with the three respective *B. anthracis* genomes. The ability to load all subdirectories is useful when the directory structure mirrors genome taxonomy. However, directories can still be specified individually. For example, `pcramp -T Bacillus/Bacillus_anthracis_str._Ames_7845 -T Bacillus/Bacillus_anthracis_str._Australia_94_167335` will load the `Bacillus_anthracis_str._Ames_7845` and `Bacillus_anthracis_str._Australia_94_167335` but would not include the `Bacillus_anthracis_str.__Ames_Ancestor__8445` genome.

Similar to the loading of target sequences that are stored in separate directories, background sequences that are stored in one or more files in a directory can be specified to PCRamp with the  *capital* `-B` command line flag.

When loading sequences that are grouped by directory (as in the above example of *Bacillus* genomes), PCRamp allows a shared target or background directory prefix to be specified using the `--T.prefix <path>` (for targets) and `--B.prefix <path>` (for backgrounds). These directory prefixes are useful for de-cluttering the command line invocations of the PCRamp file. For example, the following command line invocation:
```
./pcramp \
	-T Bacillus/anthracis/Bacillus_anthracis_str._Ames_7845 \
	-T Bacillus/anthracis/Bacillus_anthracis_str.__Ames_Ancestor__8445 \
	...
```
is equivalent to:
```
./pcramp \
	--T.prefix  Bacillus/anthracis \
	-T Bacillus_anthracis_str._Ames_7845 \
	-T Bacillus_anthracis_str.__Ames_Ancestor__8445 \
	...
```
In the second example, the use of `--T.prefix` allows the path `Bacillus/anthracis` to be specified just *once*, and all of the `-T` command arguments are treated as subdirectories of the specified target prefix path. Similarly, `--B.prefix` specifies a shared-background path that will be appended to all background genome directories. Finally, `--input.prefix` can be used to specify a common path for *both* target and background genomes.

### Sequence weights

As an extension to the traditional fasta file format, PCRamp allows for fasta files that include optional sequence weights. The default weight of a sequence is 1, but can be changed by including the string `[w=xxxx]` in the fasta defline, where `xxxx` is a floating point number greater than or equal to zero. The per-sequence weights are used to compute the fractional target and background coverage. Setting the per-sequence weights enables the design of PCR enrichment assays that preferentially amplify the target sequences that have been assigned larger weights.

## Output assay format

The PCRamp program outputs assay design results in either a human-readable text file (the default; or when specified by the `--o.text` command line flag) or a JSON-formated file (when specified by the `--o.json` commad line flag).

## Experimental constraints

PCRamp constrains PCR primer designs based on several experimental conditions:
- Buffer salt concentration, specified by the `--salt` parameter. This parameter influences the calculation of primer melting temperatures.
- Primer concentration, specified by the `--primer.strand` parameter. This parameter influences the calculation of primer/template and primer/dimer melting temperatures.
- The range of allowed melting temperatures for primer/template, primer/dimer and primer-hairpin can all be specified (see the list of command line options below). These temperatures should be adjusted by the operator to be consistent with the experimental PCR annealing temperature.

## Ensuring multiplex compatibility

PCRamp employs several strategies for ensuring multiplex compatibility. PCR primer pairs are iteratively added to the final multiplex pool of PCR primers, subject to the following constraints:
- PCR primer pairs are not allowed to form heterodimers with previously added primer pairs. The identification of heterodimers is controlled by the `--primer.dimer` melting temperature threshold.
- Individual PCR primers are not allowed to bind to the interior sequence of any amplicons produced by a previously added PCR primer pair. This effectively forbids the selection of primer pairs that might generate overlapping amplicons.
- Where possible, existing PCR primers are "reused" -- that is, if the target coverage can be increased by combining a previously selected primer oligo with a new primer oligo (rather than adding two new primers), then only a single new primer will be added to multiplex pool (assuming that both options increase the target coverage by the same amount). The design algorithm is "greedy", and selects the primer combination (either two new primers or a new plus an existing primer) that provides the largest increase in target coverage.

## Primer design and optimization

The PCRamp primer design algorithm starts by randomly selecting PCR primers. To generate a random PCR primer pair, the algorithm first randomly selects a single target sequence and then randomly selects a primer that satisfies the required length, melting temperature, hairpin and homodimer constraints. Then, the algorithm searches for a second primer that is in the correct orientation and within the allowed amplicon distance to the first primer and satisfies the length, melting temperature, hairpin, homodimer and heterodimer constraints.

After the requested number of trial PCR primer pairs have been selected (this number is determined by the `--trial` parameter), PCRamp attempts to increase the target coverage and decrease the background coverage of each primer by (a) adding and/or removing degenerate nucleotides (if allowed by the degeneracy level, `-d`) and (b) varying the primer sequence by adding/removing 5' bases (if `--optimize.5` is set) and adding/removing 3' bases (if `--optimize.3` is set). Finally each trial assay pair is assigned a score equal to the fraction of target sequences amplified minus the fraction of background sequences amplified. The primer pair with the highest score is added to the pool of multiplex primers. This process stops when (a) the requested number of PCR primer pairs has been selected (determined by the `--count` parameter), or (b) a valid primer pair cannot be found.

## Command line arguments - the gory details

The PCRamp program accepts a large number of command line arguments. However, there are *only two required arguments*: the target sequence(s) in fasta format and the name of the output file for writing results to. All other program options have (hopefully) sensible default values.

Here is a list of all of the PCRamp command arguments:

- **Sequence inputs and assay design output:** The allowed format and storage conventions for specifying the input DNA sequences for PCRamp are explained in more detail in the above section on "Input sequence formats".
	- **Required:** `-t <target fasta file>` or `-T <root directory of target subdirectories>` These arguments specify the target sequences that are to be enriched for/amplified. Multiple fasta files and target directories can be specified. 
	- Optional: `--T.prefix <directory prefix for target genomes>` Specifies a common filesystem path that will be appended as a prefix to all of the specified target directories.
	- Optional: `-b <background fasta file>` A file of background sequences, where each sequence is a distinct background.
	- Optional: `-B <root directory of background subdirectories>` A directory containing one or more fasta files or one or more subdirectories that contain fasta files. The lowest-level directory that contains fasta files will be treated as a single background and all of the sequences found within this directory will be treated as molecules within a single background.
	- Optional: `--B.prefix <directory prefix for background genomes>` Specifies a common filesystem path that will be appended as a prefix to all of the specified background directories.
	- **Required:** `-o <output file>` The name of the output file of assay designs that will be created by the PCRamp program.
	- Optional: `--o.text (output results as [poorly] structured text; default)` Specifies that the output file should be written in a human-readable, text-based format.
	- Optional: `--o.json (output results in JSON format)` Specifies that the output file should be written in the JSON file format.
	- Optional: `--input.prefix <directory prefix for both target and background input genomes>` A common filesystem path that will be appended as a prefix to *both* the target *and* the background directories that are specified with the `-T` and `-B` commands.
	- Optional: `-v <verbosity level: silent, verbose, everything>` Controls the amount of informative text written to standard error I/O stream during program execution.
- **Search and sampling strategy:**
	- Optional: `--trial <number of trials (default is 1000)>` The number of randomly-selected primer pairs that will be evaluated for each multiplex level. For example, if a 100-plex assay is requested (using the `--count` option), then 1000 trial primer pairs will be evaluated for each of the 100 multiplex levels. Increasing the number of trials typically improves the target coverage and specificity of the final assay design at the expense of a longer running time. When running on multiple computers/CPU cores, the requested number of trial assays are divided among the available parallel workers.
	- Optional: `--seed <random number seed> (default is time-based)` Specifies the random number seed that controls the sampling of random primer pairs. Please note that running in parallel on multiple computers/CPU cores can confound reproducibility, even when using the exact same seed value.
	- Optional: `--thread <maximum number of OpenMP threads> (default is all)` Limits the number of threads that PCRamp is allowed use on each computer.
	- Optional: `--optimize.5 (enable 5' oligo search to optimize assay coverage)` Performs a local search that attempts to maximize primer pair target coverage and minimize the background coverage by adding new, or removing existing, 5' bases. This option can improve the coverage and specificity of assay designs at the expense of a longer runtime.
	- Optional: `--no-optimize.5 (disable 5' oligo search to optimize assay coverage; default)` Disables 5' primer optimization.
	- Optional: `--optimize.3 (enable 3' oligo search to optimize assay coverage)` Performs a local search that attempts to maximize primer pair target coverage and minimize the background coverage by adding new, or removing existing, 3' bases. This option can improve the coverage and specificity of assay designs at the expense of a longer runtime.
	- Optional: `--no-optimize.3 (disable 3' oligo search to optimize assay coverage; default)` Disables 3' primer optimization.
	- Optional: `--optimize.top-down Search using maximally degenerate inital assay oligos (default is bottom up)` Starts the primer optimzation process by selecting maximally degenerate primers (and then reducing the degeneracy). Only applies when the allowed primer degneracy (determined the the `-d` parameter) is greater than 1.
	- Optional: `--pack.degen.max <max degen when packing words> (default is 256)` Determines the maximum allowed target and background sequence degeneracy when collecting all of the potential primer binding sites prior to primer optimization. This is intended to balance the inclusion of "legitimate" uncertain bases in input sequences (i.e. a lone iUPAC 'N' or partially degenerate base, like 'Y') and the exclusion of runs of completely uncertain bases (i.e. "NNNNNNNNNN") that would lead to spurous primer-template matches.
	- Optional: `--pack.gc.max <max fractional GC content when packing words> (default is 1, disabled)` The maximum allowed fractional G+C content (from 0 to 1) of potential primer binding sites. Used to exclude high G+C regions as primer binding sites.
	- Optional: `--pack.gc.min <min fractional GC content when packing words> (default is 0, disabled)` The minimum allowed fractional G+C content (from 0 to 1) of potential primer binding sites. Used to exclude low G+C regions as primer binding sites.
- **High-level assay and primer design constraints:**
	- Optional: `-d <max degeneracy> (default is 1)` The maximum allowed degeneracy of a primer pair. The degeneracy of a primer pair is the product of individual primer degeneracies. For example, if a primer has the sequence "YCAGCTGN", its degeneracy is 2x1x1x1x1x1x1x4 = 8. Combining this primer with different primer of degeneracy 3 (i.e. "ACGATVCGA") would yeild a primer pair with total degeneracy equal to 8x3 = 24.
	- Optional: `--salt <salt concentration> (default is 0.05` The molar salt concentration of the multiplex PCR reaction.
	- Optional: `--primer.hairpin <max oligo hairpin Tm> (default is 40)` The *maximum* allowed primer hairpin melting temperature in Celsius.
	- Optional: `--primer.dimer <max oligo hetero/homo-dimer Tm> (default is 40)` The *maximum* allowed primer homo and hetero dimer melting temperature in Celsius.
	- Optional: `--count <total number of assays to produce> (default is 100)` The desired multiplex level of the enrichment assay.
	- Optional: `--primer.size.min <minimum primer length> (default is 18)` The *minimum* allowed length of a PCR primer.
	- Optional: `--primer.size.max <maximum primer length> (default is 25)` The *maximum* allowed length of a PCR primer.
	- Optional: `--primer.tm.min <minimum primer melting temperature> (default is 50)` The *minimum* allowed primer-template melting temperature in Celsius.
	- Optional: `--primer.tm.max <maximum primer melting temperature> (default is 75)` The *maximum* allowed primer-template melting temperature in Celsius.
	- Optional: `--primer.strand <primer strand concentration> (default is 9e-07)` The molar primer strand concentration. This is the assumed concentration of *each* primer oligo in the multiplex reaction.
	- Optional: `--primer.taq-mama (use Taq MAMA rules for terminal primer mismatches; default is false)` When attempting to enrich targets with closely related backgrounds, it can be advantageous to exploit the dependancy of PCR on 3' terminal and penultimate primer mismatches. The identification of 3' primer-template mismatches that inhibit PCR is from Li et al "Genotyping with TaqMAMA", Genomics 83 (2004) 311-320.
- **Target-specific assay design constraints:**
	- Optional: `--target.amplicon.min <minimum amplicon length> (default is 80)` The *minimum* allowed target amplicon length.
	- Optional: `--target.amplicon.max <maximum amplicon length> (default is 200)` The *maximum* allowed target amplicon length.
	- Optional: `--target.threshold <target detection threshold> (default is 1)` The target detection threshold controls the number of allowed mismatches between a primer pair and each target sequence. The threshold controls the fraction of primer bases that must be a perfect match to a target sequence. *Gaps between primers and targets are not allowed*, and a primer pair is defined as detecting a target sequence if the sqrt(f x r) >= target_threshold, where "f" is the fraction of forward primer bases that are a perfect match to the template and "r" is the fraction of reverse primer bases that are perfect match to the template. 
	- Optional: `--target.search <target search multiplier> (default is 0.9)` PCRamp maximizes the target coverage and simultaneously minimizes the background coverage of each primer pair by adding degeneracies and adding and/or removing bases to the 5' and/or 3' ends of a primer sequence. To make this optimization process computationally efficient, only a relatively small set of *potential* primer binding sites are loaded into memory. The `--target.search` parameter controls the degree of complementarity a potential primer binding site must possess to be included in the assay design process. A value of 1 only includes perfect match, target/primer binding sites. A value of 0.5 includes target binding sites that match as little as 50% of the primer bases.
	- Optional: `--target.cover <*minimum* primer pair coverage> (default is 0)` Primer pairs are only accepted if they are predicted to detect the specified minimum fraction of target sequences. While useful for restricting the assay design search to high coverage primer pairs, setting the parameter to value greater than zero runs the risk of halting assay design while poorly conserved target sequences are still uncovered/unamplified.
	- Optional: `--target.ignore <defline key word to exclude a sequence>` Ignore target sequences with a fasta defline that match the supplied, case-insensitive, keyword. This command can be repeated to specify multiple keywords that will be used to exclude target sequences.
	- Optional: `--target.normalize (normalize target weights per fasta file)` Normalize the weight scores of all target sequences to sum to 1. This is only relevant when per-sequence target weights are being used.
	- Optional: `--target.size.min <length> (minimum input target length in bp)` Exclude target sequences that are *shorter* than the specified length from the assay design process.
	- Optional: `--target.size.max <length> (maximum input target length in bp)` Exclude target sequences that are *longer* than the specified length from the assay design process.
- **Background-specific assay design constraints:**
	- Optional: `--background.amplicon.min <minimum amplicon length>` Only search for background sequence amplicons that are *greater* the specified length. By default, the minimum background amplicon length is zero.
	- Optional: `--background.amplicon.max <maximum amplicon length>` Only search for background sequence amplicons that are *less* the specified length. By default, there is *no* maximum background amplicon length.
	- Optional: `--background.threshold <background detection threshold> (default is 0.8)` When counting background matches, only include background matches that have a primer pair score sqrt(f x r) >= background_threshold, where "f" is the fraction of forward primer bases that match the background sequence and "r" is the fraction of reverse primer bases that match the background sequence.
	- Optional: `--background.search <background search multiplier> (default is 0.9)` PCRamp maximizes the target coverage and simultaneously minimizes the background coverage of each primer pair by adding degeneracies and adding and/or removing bases to the 5' and/or 3' ends of a primer sequence. To make this optimization process computationally efficient, only a relatively small set of potential primer binding sites are loaded into memory. The `--background.search` parameter controls the degree of complementarity a potential primer binding site must possess to be included in the assay design process. A value of 1 only includes perfect match, background/primer binding sites. A value of 0.5 includes background binding sites that match as little as 50% of the primer bases.
	- Optional: `--background.cover <*maximum* primer pair coverage> (default is 0)` The maximum fraction of background sequences that PCR primer pair is allowed to amplify. By default, primer pairs are required to be 100% specific and not detect any background sequences.
	- Optional: `--background.ignore <defline key word to exclude a sequence>` Ignore background sequences with a fasta defline that match the supplied, case-insensitive, keyword. This command can be repeated to specify multiple keywords that will be used to exclude background sequences.
	- Optional: `--background.normalize (normalize background weights per fasta file)` Normalize the weight scores of all background sequences to sum to 1. This is only relevant when per-sequence background weights are being used and when the `--background.cover` value is greater than zero.
	- Optional: `--background.size.min <length> (minimum input background length in bp)` Exclude background sequences that are *shorter* than the specified length from the assay design process.
	- Optional: `--background.size.max <length> (maximum input background length in bp)` Exclude background sequences that are *longer* than the specified length from the assay design process.
