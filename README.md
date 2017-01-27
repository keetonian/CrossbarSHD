# CrossbarSHD
This repository contains the code used for the CrossbarSHD research project.

## crossbar_shd.cpp
This file contains code that takes in a reference genome (fasta format) and a file of reads (fastq format) and finds potential mapping locations.
To use: download this repository and run "make" to compile the program.

  -g path to reference genome

  -r path to read file

  -t error threshold (default = 0)

  -s shift distance (default = 0)

  -4 use 4-bit encodings

  -16 use 16-bit encodings (more accurate, default)

By default, it will only look for perfect hamming distance matches. Change the error checking with the -s and -t flags.


