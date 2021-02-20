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

- `-i`: Input newick tree file
- `-o`: (optional) Output newick tree file
- `-d`: (optional) Delimiter separating species name from rest of leaf label. Default None.
- `-m`: (optional) Output only single tree (discarding smallest duplicate clades).
- `-s`: (optional) Discards duplicate clades that have leafsets which are subsets of their counterpart.
- `-t`: (optional) Trim duplicate leaves under each duplication event from smallest clade.
- `-tb`: (optional) Trim duplicate leaves under each duplication event. Gives two single copy trees (trimming from smallest/largest).
- `-r`: (optional) Randomly samples single-copy trees from gene family trees
- `-rm`: (optional) Choose number of single copy trees per gene tree (used with `-r`). Options `linear`, `exp`, or number. (default: 5)
- `-v`: (optional) Enable verbose output

### Example

```
python tag_decomp.py -i example/gtrees-mult.trees
```
