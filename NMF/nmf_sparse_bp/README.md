# NMF Block-Pivot Implementation

This program performs Non-negative Matrix Factorization (NMF) using block-pivot (BP) updates with support for different 
regularization and stopping rules. The code in this directory is based on the work of Jingu Kim. If used, please include
the citations below:

Reference:
Jingu Kim and Haesun Park, Toward Faster Nonnegative Matrix Factorization: A New Algorithm and Comparisons,
In Proceedings of the 2008 Eighth IEEE International Conference on Data Mining (ICDM'08), 353-362, 2008

---
## Dependencies

- C++11 or higher
- [GSL 2.6](https://www.gnu.org/software/gsl/)

Make sure to set your `LD_LIBRARY_PATH` to include the GSL library:

```bash
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/path/to/gsl/lib
```

---
## Compilation
```bash
# Create a build directory
mkdir -p build
cd build

# Configure the project
cmake ..

# Build the executable
make
mv run_nmf_bp ..

# After successful compilation, the executable will be at:
# ./run_nmf_bp
```

---
## Command-Line Arguments

### Required arguments

| Flag | Long option   | Description                        |
|------|---------------|------------------------------------|
| -x   | --data        | Path to input matrix file          |
| -n   | --nRows       | Number of samples/cells (rows)     |
| -N   | --nCols       | Number of features/genes (columns) |
| -k   | --nFactors    | Number of factors                  |

### Optional arguments

| Flag | Long option   | Default          | Description                                            |
|------|---------------|------------------|--------------------------------------------------------|
| -o   | --output      | `.`              | Output directory                                       |
| -p   | --inputPrefix | `""`             | Prefix for reading previous U/V matrices               |
| -r   | --seed        | 1010             | RNG seed                                               |
| -s   | --verbose     | false            | Print verbose logs                                     |
| -M   | --maxIter     | 100              | Maximum number of iterations                           |
| -m   | --minIter     | 20               | Minimum number of iterations                           |
| -t   | --tolerance   | 0.001            | Convergence tolerance                                  |
| -a   | --alpha       | 0.1              | Regularization alpha                                   |
| -b   | --beta        | 0.1              | Regularization beta                                    |
| -R   | --reg         | sparse           | Regularization type: `normal`, `sparse`, `regularized` |
| -S   | --stop        | normalized_pgrad | Stop rule: `normalized_pgrad`, `pgrad`                 |
| -h   | --help        |                  | Print usage information                                |

---

## Usage Examples

### Minimal run (sparse BP, default alpha/beta)

```bash
./run_nmf_bp --data input/toy/A.txt --nRows 95 --nCols 120 --nFactors 3 --output output/sparse
```

### NMF with normal regularization

```bash
./run_nmf_bp --data input/toy/A.txt --nRows 95 --nCols 120 --nFactors 3 --reg normal --output output/norm
```

### NMF with regularized BP (alpha=0.1, beta=0.1)

```bash
./run_nmf_bp --data input/toy/A.txt --nRows 95 --nCols 120 --nFactors 3 --reg regularized --output output/reg
```

### NMF with stop method `pgrad`

```bash
./run_nmf_bp --data input/toy/A.txt --nRows 95 --nCols 120 --nFactors 3 --stop pgrad --output output/stop_meth_grad
```

### NMF sparse BP with custom alpha/beta

```bash
./run_nmf_bp --data input/toy/A.txt --nRows 95 --nCols 120 --nFactors 3 --alpha 1 --beta 1 --output output/sparse_1_1
```

### Verbose run

```bash
./run_nmf_bp --data input/toy/A.txt --nRows 95 --nCols 120 --nFactors 3 --output output/sparse --verbose
```

---

## Output

- `usage.txt` : Memory and time usage
- `U.txt` : Final row embedding
- `V.txt` : Final column embedding

---

## Notes

1. All required arguments **must be provided**.
2. Optional arguments have **default values** if not specified.
3. Flags can appear **in any order**.
4. Nested output directories will be created automatically.

