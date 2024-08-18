# PARALLEL AND SEQUENTIAL VERSIONS OF VF2PP ALGORITHM

## DATA GENERATOR
### Requirements 
To generate all graphs, `python` must be installed along with the `networkx library`. You can install `networkx` using the command:

```bash
pip3 install networkx
```
### Run
To execute the `generator`, use the following command:

```bash 
python3 data/graph_generator.py
```

## VF2PP ALGORITHM
### Requirements 
To execute the parallel version, `CUDA` must be installed. It is recommended to upload the repository to google drive and run notebooks using `google colab`.
### Compile and Run
To compile the whole project use the following command:
```bash 
make
```
You can compile and execute the parallel version by running the `vf2pp_parallel.ipynb` notebook, and perform all performance analyses between the sequential and parallel versions by running the `comparison.ipynb` notebook.
### TEST
To compile all tests use the following command:
```bash 
make test
```
You can compile and execute all tests by running the `test.ipynb` notebook. 
### Cleaning
All objects, executables and temporary files can be removed by running the following command:
```bash 
make clean
```