# Synthetic CRE Generator (`generate_synCREs.py`)

This script generates synthetic cis-regulatory elements (CREs) for ENGRAM experiments, allowing you to specify motifs, repeat numbers, spacing, and targeting system (DTT or HEK3). Output is a tab-separated table ready for oligonucleotide ordering.

## Requirements

- Python 3.6+
- No external dependencies (uses only standard library)

## Usage

Run the script from the command line:

```bash
python3 generate_synCREs.py [OPTIONS]
```

### Options

- `-m`, `--motif`  
  Specify a single motif sequence (DNA).

- `-n`, `--motif-name`  
  Name for the motif (used with `--motif`). Default: `motif`.

- `-f`, `--motifs-file`  
  Tab-separated file with motif names and sequences. Columns: `name` and `sequence`.

- `-r`, `--repeats`  
  Number of repeats (single number, range like `1-7`, or comma-separated list like `1,3,5`). Default: `6`.

- `-s`, `--spacing`  
  Number of bases between motifs. Default: `15`.

- `-d`, `--distance`  
  Distance from last motif to minimal promoter. Default: `0`.

- `--system`  
  Targeting system: `DTT` or `HEK3`. Default: `HEK3`.

- `--background-seq`  
  Background DNA sequence for CRE generation. Default: built-in sequence.

- `-o`, `--output`  
  Output file path (TSV). If not specified, prints to stdout.

- `-v`, `--verbose`  
  Enable verbose output.

- `-h`, `--help`  
  Show help message.

### Examples

**Single motif:**
```bash
python3 Nprotocol_scripts/generate_synCREs.py \
  --motif AGGTCA \
  --motif-name DR1 \
  --repeats 3-5 \
  --spacing 10 \
  --distance 5 \
  --system HEK3 \
  --output synCREs.tsv \
  --verbose
```

**Motif file:**
Create a file `motifs.tsv`:
```
name	sequence
DR1	AGGTCA
DR4	AGGTCAAGGTCA
```
Run:
```bash
python3 Nprotocol_scripts/generate_synCREs.py \
  --motifs-file motifs.tsv \
  --repeats 2-4 \
  --spacing 15 \
  --system DTT \
  --output synCREs.tsv
```

### Output

A TSV file with columns:
- `CRE_Name`
- `Motif_Sequence`
- `Repeats`
- `System`
- `ready_to_order_oligo`

## Help

For full options, run:
```bash
python3 Nprotocol_scripts/generate_synCREs.py --help
```

## Author

Wei Chen