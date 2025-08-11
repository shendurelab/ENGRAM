# ENGRAM CRE Barcode Generator (`generate_ENGRAM_CRE_bc.py`)

This script generates synthetic cis-regulatory elements (CREs) for ENGRAM experiments or processes custom CREs provided by users. It assigns unique barcodes from predefined barcode libraries and outputs oligonucleotides ready for ordering. You can specify motifs, repeat numbers, spacing, targeting system (DTT or HEK3), or provide your own custom CRE sequences.

## Features

- **synCRE Generation**: Create synthetic CREs from motifs with configurable repeats and spacing
- **Custom CRE Processing**: Process user-provided CRE sequences (200-300bp recommended)  
- **Barcode Assignment**: Automatically assigns unique barcodes from 5N or 8N barcode libraries
- **System Support**: Compatible with both DTT and HEK3 targeting systems
- **Reproducible**: Optional random seed for consistent barcode assignment

## Requirements

- Python 3.6+
- No external dependencies (uses only standard library)
- Barcode files: `ENGRAM_5N_good_barcodes.csv` and `ENGRAM_8N_good_barcodes.csv`

## Usage

Run the script from the command line:

```bash
python3 generate_ENGRAM_CRE_bc.py [OPTIONS]
```

### Input Modes

**1. Single motif synCRE generation:**

```bash
python3 generate_ENGRAM_CRE_bc.py --motif AGGTCA --motif-name DR1 [OPTIONS]
```

**2. Multiple motifs synCRE generation:**

```bash
python3 generate_ENGRAM_CRE_bc.py --motifs-file motifs.tsv [OPTIONS]
```

**3. Custom CRE processing:**

```bash
python3 generate_ENGRAM_CRE_bc.py --custom-cres-file custom_cres.tsv [OPTIONS]
```

### Options

#### Input Options (mutually exclusive)
- `-m`, `--motif`  
  Specify a single motif sequence (DNA) for synCRE generation.

- `-f`, `--motifs-file`  
  Tab-separated file with motif names and sequences for synCRE generation.

- `-c`, `--custom-cres-file`  
  Tab-separated file with custom CRE names and sequences (200-300bp recommended).

#### Motif Options
- `-n`, `--motif-name`  
  Name for the motif (used with `--motif`). Default: `motif`.

#### synCRE Generation Parameters
- `-r`, `--repeats`  
  Number of repeats (single number, range like `1-7`, or comma-separated list like `1,3,5`). Default: `6`.

- `-s`, `--spacing`  
  Number of bases between motifs. Default: `15`.

- `-d`, `--distance`  
  Distance from last motif to minimal promoter. Default: `0`.

- `--background-seq`  
  Background DNA sequence for CRE generation. Default: built-in sequence.

#### System and Barcode Options
- `--system`  
  Targeting system: `DTT` or `HEK3`. Default: `HEK3`.

- `--barcode-file`  
  CSV file with barcodes. If not specified, auto-selects based on system:
  - HEK3: `ENGRAM_5N_good_barcodes.csv`
  - DTT: `ENGRAM_8N_good_barcodes.csv`

- `--seed`  
  Random seed for barcode assignment (for reproducibility).

#### Output Options
- `-o`, `--output`  
  Output file path (TSV). If not specified, prints to stdout.

- `-v`, `--verbose`  
  Enable verbose output.

- `-h`, `--help`  
  Show help message.

### File Formats

#### Motifs File (`motifs.tsv`)
Tab-separated file with columns:
```
name	sequence
DR1	AGGTCA
DR4	AGGTCAAGGTCA
AP1	TGAGTCA
```

#### Custom CREs File (`custom_cres.tsv`)  
Tab-separated file with columns:
```
name	sequence
enhancer1	GCAAACAGTCCCCATGTTCACATTA...
enhancer2	ATGCGTCTCAGTCTGCCATCAAAT...
control	TTTTTTTTTTTTTTTTTTTTTTTT...
```

#### Barcode Files (`ENGRAM_5N_good_barcodes.csv` / `ENGRAM_8N_good_barcodes.csv`)
CSV file with header `pBC_Seq`:
```
pBC_Seq
AAAAA
AAAAT
AAAAC
...
```

### Examples

**1. Generate synCREs from a single motif:**
```bash
python3 generate_ENGRAM_CRE_bc.py \
  --motif AGGTCA \
  --motif-name DR1 \
  --repeats 3-5 \
  --spacing 10 \
  --distance 5 \
  --system HEK3 \
  --output synCREs.tsv \
  --verbose
```

**2. Generate synCREs from multiple motifs:**
```bash
python3 generate_ENGRAM_CRE_bc.py \
  --motifs-file motifs.tsv \
  --repeats 2,4,6 \
  --spacing 15 \
  --system DTT \
  --output synCREs.tsv \
  --verbose
```

**3. Process custom CREs:**
```bash
python3 generate_ENGRAM_CRE_bc.py \
  --custom-cres-file my_enhancers.tsv \
  --system HEK3 \
  --output custom_oligos.tsv \
  --barcode-file my_barcodes.csv \
  --verbose
```

### Output

A TSV file with columns:
- `CRE_Name`: Name of the CRE or generated name
- `Motif_Sequence` / `CRE_Type`: Motif sequence (synCRE) or "N/A" (custom)
- `Repeats`: Number of repeats (synCRE) or "N/A" (custom)  
- `System`: DTT or HEK3
- `Barcode`: Assigned unique barcode
- `ready_to_order_oligo`: Complete oligonucleotide sequence

The script also outputs a summary of CRE-barcode associations when verbose mode is enabled.

## Help

For full options, run:
```bash
python3 generate_ENGRAM_CRE_bc.py --help
```

## Author

Wei Chen