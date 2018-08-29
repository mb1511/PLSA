# MMLSA Masive-Multiple Locus Sequence Analysis

## Description

A Python script to perform MLSA on a much larger scale, using as many conserved genes as possible and repeated
random sampling to avoid bias

## Installation

This package utilises a variety of bionformatics tools and they must be built prior running setup

Download and build the BLAST package from [NCBI](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)

If possible - use [gpu-blast](http://archimedes.cheme.cmu.edu/?q=gpublast) for increased search speed for the f_con module

Download and build [MUSCLE](http://www.drive5.com/muscle/downloads.htm)

Download and build [MEGACC](http://www.megasoftware.net/) (commandline version of MEGA, not GUI; tested with v7)

If not already installed, run:
	
	pip install numpy
	pip install dill
	pip install ete3

Run:

	./setup.py install PATH_TO/blast-x.x.xx+/bin PATH_TO/muscle_executable PATH_TO/megacc_executable

## Usage

Use the MEGA prototyper (included with MEGACC) to generate a make_tree.mao file; one is included in the examples folder for 
expected output and if suitable can be used directly. Make sure number of threads is set to 1 - allows to run 1 instance of
MEGACC on each cpu core.

	import mmlsa
	
	if __name__ == '__main__':
		mmlsa.run(proteome_directory, output_directory, blast_data_storage_directory)

Proteome directory 	 = where the raw proteome files (.faa files usually) are kept

Output directory 	 = where to output alignments, trees, and the final comparison (usually named dists.txt)

Blast data directory = where the blast databases and conserved gene lists are built and stored

Other args:

	ext 	= file extension to search for in proteome dir, default = faa
	
	size	= sample size for repeats, default = 10
	
	repeats	= number of repeats to perform, default = 1000
	
	cov		= cut-off value (% similarity) for adding conserved sequences, default = 65

	
See examples for more information.

## Other Notes

This package originates from a larger project that had HPC integration for performing the 1000 odd repeats on
nodes containing 32 cores - would advise to do if possible so not to jam up your PC for a considerable length 
of time. The f_con module runs in an acceptable amount of time, usually, and is faster if there are fewer 
consevered genes
