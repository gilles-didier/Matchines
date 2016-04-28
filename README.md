# Matching machines
Asymptotic behavior of pattern matching algorithms and optimal algorithms

Directory "data" contains the texts used in the articles.

Directory "src" contains the C sources of the software.

-----------------------------
| MATCHING MACHINES PACKAGE |
-----------------------------

--------------------------
REQUIREMENT

	The package needs the superlu library.

--------------------------
COMPILING

	Just type
	> make all
	to build all the binaries.

--------------------------
DESCRIPTION

	'theoretical' computes the asymptotic speeds wrt a given Bernoulli model of 9 standard algorithms (Naive, Morris-Pratt, Knuth-Morris-Pratt, Quicksearch, Horspool, FJSS, TVSBS, EBOM and Hashq), heuristic matching machines up to a given order, the fastest machine (wrt to the Bernoulli model) for patterns of length <= 4. If option '-u' is set, it provides the speeds of the heuristic and the fastest matching machines computed from the uniform Bernoulli model. The asymptotic speeds are computed for a single pattern if it is given (option '-p'), for N samples of a given length (set with option '-l') if option '-s N' is provided and for all the patterns of a given length. Results are provided as a table saved in CSV or LaTeX tabular format (option '-f').

	'empirical' works in a very similar way but deals with the average speeds over a text file given as parameter.

	'control' is for verification purpose. It computes 4 kinds of speeds:
		1) the asymptotic speeds wrt a given Bernoulli model B,
		2) the average speeds over a random text drawn from B of the corresponding matching machines (the length of the random text can be set with option '-r'),
		3) the average speeds over a random text drawn from B of the corresponding expanded matching machines (the length of the random text can be set with option '-r'),
		4) the average speeds over a random text drawn from B of the initial algorithms as implemented in SMART (S. Faro and T. Lecroq, http://www.dmi.unict.it/~faro/smart/).
	
		Each kind of speed is display on a line.
		Speeds 4) are computed only for standard algorithms (other approaches exist only as matching machines).
		Some small variations between speeds 2) and 3-4) are allowed (matching machines corresponding to standard algorithms do not always deal exactly in the same way with the ends of texts).

	'draw' provides lattices and matching machines as .gv format which can be display wit xdot and graphviz utilities (http://www.graphviz.org/) for illustration purposes.

	The format of the Bernoulli model files is basic. A line for each symbol containing the symbol, any separator and its probability, for instance:
		a 0.1
		b 0.9
		

--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	theoretical - compute the asymptotic speed of various pattern matching algorithms under a Bernoulli model
	
SYNOPSIS
	theoretical [OPTIONS]

DESCRIPTION
	return a table containing the asymptotic speeds under a Bernoulli model of a given pattern or of all the patterns of a given size.
	options -a and -b have to be filled
	
	-a [ALPHABET]
		set the alphabet
	-p [PATTERN]
		set the pattern
	-l [LENGTH]
		set the length of the patterns to be tested (not used if a pattern is given via option -p)
	-s [NUMBER]
		set the number of pattern samples to be tested (not used if a pattern is given via option -p)
		If this option is not used, all the pattern of the selected length are tested.
	-x [FILE]
		set the file name of the output (default is 'table_result.xxx')
	-n [ORDER]
		set the max order of the heuristic
	-u
		evaluate heuristic and the fastest strategy computed from the uniform Bernoulli model
	-f [format]
		set the format of the output (l: LaTex tabular format, c/default: CSV format)
	-r [TYPE]
		type of random used for sampling patterns (active with option 's'):
			m	use Bernoulli model from option b
			u	use uniform model
	-b [FILE]
		set the file containing the Bernoulli model
	-h
		display help

--------------------------

NAME
	empirical - compute the average speed over a text or a binary file of various pattern matching algorithms
	
SYNOPSIS
	empirical [OPTIONS] [FILE]

DESCRIPTION
	return a table containing the average speed over a text or a binary file of various pattern matching algorithms wrt a given pattern or of all the patterns of a given size.
	
	-a [ALPHABET]
		set the alphabet (if the option is not used then the alphabet is determined from the text file and the pattern)
	-p [PATTERN]
		set the pattern
	-l [LENGTH]
		set the length of the patterns to be tested (not used if a pattern is given via option -p)
	-s [NUMBER]
		set the number of pattern samples to be tested (not used if a pattern is given via option -p)
		If this option is not used, all the pattern of the selected length are tested.
	-d
		indicate that the file is binary (considered as a text file otherwise)
	-x [FILE]
		set the file name of the output (default is 'table_result.xxx')
	-n [ORDER]
		set the max order of the heuristic
	-u
		evaluate heuristic and the fastest strategy computed from the uniform Bernoulli model
	-f [FORMAT]
		set the format of the output (l: LaTex tabular format, c/default: CSV format)
	-r [TYPE]
		type of random used for sampling patterns (active with option 's'):
			t	sample a position uniformly in the text and read a pattern from it
			m	use Bernoulli model from option b
			u	use uniform model
	-L [LENGTH]
		read only the prefix of specified length from the file (the whole file if not used)
	-b [FILE]
		set the file containing the Bernoulli model (if the option is not used then the model is determined from the text file)
	-h
		display help

--------------------------

NAME
	control - compute the asymptotic speed and average speed over random texts of various pattern matching algorithms under a Bernoulli model
	
SYNOPSIS
	control [OPTIONS]

DESCRIPTION
	return a table containing
		- the asymptotic speeds under a Bernoulli model of a given pattern or of all the patterns of a given size,
		- the average speeds of the corresponding matching machines over a random text under the same model,
		- the average speeds of the initial algorithms over a random text under the same model.
	options -a and -b have to be filled
	
	-a [ALPHABET]
		set the alphabet
	-p [PATTERN]
		set the pattern
	-l [LENGTH]
		set the length of the patterns to be tested (not used if a pattern is given via option -p)
	-s [NUMBER]
		set the number of pattern samples to be tested (not used if a pattern is given via option -p)
		If this option is not used, all the pattern of the selected length are tested.
	-x [FILE]
		set the file name of the output (default is 'table_result.xxx')
	-n [ORDER]
		set the max order of the heuristic
	-f [format]
		set the format of the output (l: LaTex tabular format, c/default: CSV format)
	-b [FILE]
		set the file containing the Bernoulli model
	-L [LENGTH]
		set the length of the random text used for testing
	-h
		display help

--------------------------

NAME
	draw - return the full lattice and the matching machines of various pattern matching algorithms in dot '.gv' format
	
SYNOPSIS
	draw [OPTIONS] [PATTERN]

DESCRIPTION
	return the full lattice and the matching machines of various pattern matching algorithms corresponding to a given pattern in dot '.gv' format
	
	-a [ALPHABET]
		set the alphabet
	-b [FILE]
		set the file containing the Bernoulli model
	-n [ORDER]
		set the max order of the heuristic
	-h
		display help

--------------------------
