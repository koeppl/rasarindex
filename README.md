# 🗂️e  Accessing the Suffix Array via Φ−1 Forest


## ✔️ Prerequisites

 - libdivsufsort
 - g++

## 🚀 Complete Test Run

```shell
git clone https://github.com/koeppl/randSAbench.git
cd randSAbench
submodule update --init --recursive
```

In the file `bench.sh`,
you have to adjust the variable `kFullFasta` for the path to the FASTA file you want to index,
and `sequences` for the number of sequences you want to extract from this file.
Then you can run `bench.sh` measuring the query time for suffix array access with our proposed method and the standard method of the r-index.
Note that the default is to also build the plain suffix array to check whether the reported entry is correct.
Building the plain suffix array will not work with large inputs and a modest amount of memory.

