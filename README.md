# GDL Tagging Decomposition

Python script for decomposing multi-copy trees in newick format using various methods based of tagging.

## Dependencies

- Python 3
- [treeswift](https://github.com/niemasd/TreeSwift)

## Usage

**Input**: File containing list of multi-copy trees in newick format

**Output**: File containing list resulting list of trees after decomposition in newick format

```
python tag_decomp.py -i <input_file> -o <ouput_file> -d <delimiter>
```

### Arguments

#### Required

- `-i`: Input newick tree file

#### Optional

- `-o`: Output newick tree file
- `-d`: Delimiter separating species name from rest of leaf label. Default None.
- `-m`: Output only single tree (discarding smallest duplicate clades).
- `-s`: Discards duplicate clades that have leafsets which are subsets of their counterpart.
- `-t`: Trim duplicate leaves under each duplication event from smallest clade.
- `-tb`: Trim duplicate leaves under each duplication event. Gives two single copy trees (trimming from smallest/largest).
- `-r`: Randomly samples single-copy trees from gene family trees
- `-rm`: Choose number of single copy trees per gene tree (used with `-r`). Options `linear`, `exp`, or number. (default: 5)
- `-v`: Enable verbose output
- `-rp`: Remove in-paralogs before rooting/scoring
- `--trivial`: Includes trivial trees in decomposition output
- `--outgroups`: Write outgroups (including ties) to txt file. (Might make program slower).

### Example

```
python tag_decomp.py -i example/gtrees-mult.trees
```
