r-index: the run-length BWT index
===============

Modified by: Taher Mun
Fork of `r-index` (see below) that uses pfBWT as the BWT-construction algorithm
(https://gitlab.com/manzai/Big-BWT/tree/master) and to index+locate DNA
sequences. 

### To compile:

1) Clone repository

```
git clone --recursive https://github.com/alshai/r-index
```

2) Install required packages:

```
apt-get update -qq && \
apt-get install -y zlib1g-dev git cmake build-essential python3
```

3) Compile BigBWT and add to path

```
cd Big-BWT
make
export PATH=$(pwd):$PATH
cd ..
```

4) compile and install

```
mkdir build
cd build
cmake ..
make
make install
```

### To run:

To build index from a fasta file (outputs to `input.fa.ri`):

```
ri-buildfasta -b <bigbwt|sais|from_bwt> input.fa
```

To count queries in a fast[a|q] file:

```
ri-align count index_prefix reads.fq
```

To locate queries in a fast[a|q] file:

```
ri-align --max-hits (default:-1) --max-range (default:-1) count index_prefix reads.fq
```



ORIGINAL README
==============

Author: Nicola Prezza (nicola.prezza@gmail.com)
Joint work with Travis Gagie and Gonzalo Navarro

cite as:

Gagie T, Navarro G, Prezza N. Optimal-time text indexing in BWT-runs bounded space. In Proceedings of the Twenty-Ninth Annual ACM-SIAM Symposium on Discrete Algorithms, SODA 2018, New Orleans, NA, USA, January 7-10 2017.

### Brief description

The r-index is the first full-text index of size O(r), r being the number of BWT runs of the input text (of size n), supporting fast (almost optimal) locate of pattern occurrences. The r-index employs a novel suffix array sampling of size 2r; in classical FM-indexes, this sampling would result in a locate time of Omega(n/r) per occurrence. The r-index, on the other hand, reduces this time to O(log(n/r)).

Let s be the alphabet size and fix a constant eps>0. The r-index offers the following tradeoffs:

- Space: r * ( log s + (1+eps)log(n/r) + 2log n ) bits
- Count time: O( (m/eps) * (log (n/r) + log s) )
- Locate time: After count, O( log(n/r) ) time per occurrence 

On very repetitive datasets, the r-index locates orders of magnitude faster than the RLCSA (with a sampling rate resulting in the same size for the two indexes).

NEWS: refactored locate strategy. Let (l,r) be the SA range. Now, the index first finds SA[r] and then applies function Phi to locate SA[r-1], SA[r-2], ..., SA[l]. This is both faster and more space efficient than the strategy originally implemented and described in the paper.

### Download

To clone the repository, run:

> git clone http://github.com/nicolaprezza/r-index

### Compile

The library has been tested under linux using gcc 6.2.0. You need the SDSL library installed on your system (https://github.com/simongog/sdsl-lite).

We use cmake to generate the Makefile. Create a build folder in the main r-index folder:

> mkdir build

run cmake:

> cd build; cmake ..

and compile:

> make

### Run

After compiling, run 

>  ri-build input

This command will create the r-index of the text file "input" and will store it as "input.ri". Use option -o to specify a different basename for the index file. 

Run

> ri-count index.ri patterns

to count number of occurrences of the patterns, where <patterns> is a file containing the patterns in pizza&chili format (http://pizzachili.dcc.uchile.cl/experiments.html). To generate pattern files, use the tool http://pizzachili.dcc.uchile.cl/utils/genpatterns.c 

Run

> ri-locate index.ri patterns

to locate all occurrences of the patterns.

Be aware that the above executables are just benchmarking tools: no output is generated (pattern occurrences are deleted after being extracted and not printed to output).

