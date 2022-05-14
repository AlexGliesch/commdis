# Supplementary material for *A multistart alternating tabu search for commercial districting*

This repository holds the source code, detailed result tables and instances for the following paper:

> **A multistart alternating tabu search for commercial districting**<br>
> Alex Gliesch, Marcus Ritt, Mayron C. O. Moreira (2018)<br>
> European Conference on Evolutionary Computation in Combinatorial Optimization, pp. 158--173 <br>
> https://doi.org/10.1007/978-3-319-77449-7_11

**Bibtex**

```bibtex
@inproceedings{Gliesch.etal/2018,
  title        = {A Multistart Alternating Tabu Search for Commercial Districting},
  author       = {Gliesch, Alex and Ritt, Marcus and Moreira, Mayron CO},
  booktitle    = {European Conference on Evolutionary Computation in Combinatorial Optimization},
  pages        = {158--173},
  year         = {2018},
  organization = {Springer},
  doi          = {10.1007/978-3-319-77449-7_11}
}
```

Please use the reference above if you use this material in your research.

## Dependencies 

- [Boost](https://www.boost.org/)
- [fmtlib](https://github.com/fmtlib/fmt)

## Running the code 

1. Unpack the instances in `instances.tar.gz`.
2. Compile the code under `src` using `make release`. 
3. Run using `./commdis --in {instance} --time {timeLimit} --out {outFile}`. For more options, see `--help`. 

Our reimplementation of Ríos-Mercado & Escalante (2016)'s heuristic will be uploaded soon.

## Instance generator 

1. Compile the code under `instance-generator` using `make`. 
1. Run `./generate --n {numVertices} --p {numDistricts} --tau {balancingTolerance} --seed {randomSeed}`.
1. By default it generates a Delaunay graph. For more options, see `--help`.
