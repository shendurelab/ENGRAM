#!/usr/bin/env python3
"""
ENGRAM PE Analysis - 5XTAPE WNT/DOX (Simplified)
Author: Will Chen
this code is first written by Will Chen and modified by Claude

This script processes sequencing data from ENGRAM experiments to analyze 
prime editing events in 5XTAPE constructs.
"""

import argparse
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import Align


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Process ENGRAM 5XTAPE WNT/DOX sequencing data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('sample', help='Sample name (without .txt extension)')
    parser.add_argument('workdir', type=Path, help='Working directory containing input files')
    # Removed verbose argument
    parser.add_argument('-o', '--output', type=Path, help='Output file path (default: workdir/sample_bc_count.csv)')
    parser.add_argument('-b', '--barcode-file', type=Path, 
                       help='Tab-separated file with barcode sequences and names (barcode<tab>name)')
    
    return parser.parse_args()
# Constants - predefined parameters as you liked
TAPE_5X = 'GGATGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCGCTTTAAGGCCGGTCCTAGCAA'
TAPE_6X = 'GGATGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCGCTTTAAGGCCGGTCCTAGCAA'
TAPE_4X = 'GGATGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCGCTTTAAGGCCGGTCCTAGCAA'

# Sequence patterns
UMI_PATTERN = re.compile(r".ACACC([ATCG]{8})(GGAT.*)")
EDIT_PATTERN = re.compile(r'CG------TGA')

# Barcode mapping (default - can be overridden by file)
DEFAULT_BARCODE_MAP = {'GTT': 'Dox', 'ACA': 'Wnt'}

# Alignment parameters
ALIGNMENT_PARAMS = {
    'match': 5,
    'mismatch': -4,
    'gap_open': -10,
    'gap_extend': -0.5
}

# Position mapping for insertions (from original script)
INSERT_POSITIONS = {
    17: 1,
    31: 2, 37: 2,
    45: 3, 51: 3, 57: 3,
    59: 4, 65: 4, 71: 4, 77: 4,
    73: 5, 79: 5, 85: 5, 91: 5, 97: 5
}


# Analysis thresholds
MIN_READ_LENGTH = 55
MIN_SCORE_THRESHOLD = 325

# Define aligner once at module level
aligner = Align.PairwiseAligner()
aligner.match_score = ALIGNMENT_PARAMS['match']
aligner.mismatch_score = ALIGNMENT_PARAMS['mismatch']
aligner.open_gap_score = ALIGNMENT_PARAMS['gap_open']
aligner.extend_gap_score = ALIGNMENT_PARAMS['gap_extend']
aligner.end_gap_score = 0.0



def align_tape(seq, ref):
    """
    Align a sequence to a reference tape using Bio.pairwise2.
    Returns aligned sequence, aligned reference, and score.
    """
    alignments = aligner.align(seq, ref)
    aln = alignments[-1]
    # Extract aligned sequences with gaps, as in original code
    aligned_seq = aln[0][:-20]
    aligned_ref = aln[1][:-20]
    score = aln.score
    return aligned_seq, aligned_ref, score
def load_barcode_map(file_path):
    """
    Load barcode mapping from a tab-separated file.
    Expected format: barcode_sequence<tab>barcode_name
    
    Args:
        file_path: Path to the barcode mapping file
        
    Returns:
        Dictionary mapping barcode sequences to names
    """
    barcode_map = {}
    
    try:
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):  # Skip empty lines and comments
                    continue
                
                parts = line.split('\t')
                if len(parts) != 2:
                    print(f"Warning: Invalid format in barcode file line {line_num}: {line}")
                    continue
                
                barcode_seq, barcode_name = parts[0].strip().upper(), parts[1].strip()
                barcode_map[barcode_seq] = barcode_name
                
    except Exception as e:
        print(f"Error reading barcode file {file_path}: {e}")
        sys.exit(1)
    
    return barcode_map

def simple_hamming_distance(seq1: str, seq2: str) -> int:
    """
    Calculate Hamming distance between two sequences of equal length.
    Much faster than edit distance for our specific use case.
    """
    if len(seq1) != len(seq2):
        return 999  # Return high number for length mismatch
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def min_edit(barcode_dict, seq):
    """Find closest matching barcode with minimum edit distance (simplified one-liner)."""
    return min(((s, simple_hamming_distance(s, seq)) for s in barcode_dict.keys()), key=lambda x: x[1])


def main():
    """Main execution function - simplified, no logging overhead."""
    args = parse_arguments()
    
    # Load barcode mapping
    if args.barcode_file:
        if not args.barcode_file.exists():
            print(f"Error: Barcode file does not exist: {args.barcode_file}")
            sys.exit(1)
        barcode_map = load_barcode_map(args.barcode_file)
    else:
        barcode_map = DEFAULT_BARCODE_MAP
    # Validate inputs
    if not args.workdir.exists():
        print(f"Error: Working directory does not exist: {args.workdir}")
        sys.exit(1)
    input_file = args.workdir / f"{args.sample}.txt"
    if not input_file.exists():
        print(f"Error: Input file does not exist: {input_file}")
        sys.exit(1)
    
    # Parse input file    
    sequence_data = {}
    with input_file.open('r') as f:
        for line_num, line in enumerate(f, 1):
            try:
                parts = line.rstrip().split('\t')
                read_name, _, read_seq = parts[0], parts[1], parts[2]
                if len(read_seq) > MIN_READ_LENGTH:
                    match = UMI_PATTERN.search(read_seq)
                    if match:
                        umi, sequence = match.groups()
                        sequence_data[read_name] = [None] * 4  # Match original size
                        sequence_data[read_name][0] = sequence
                        sequence_data[read_name][1] = TAPE_5X
                        sequence_data[read_name][2] = umi
            except Exception as e:
                
                continue
    
    if not sequence_data:
        print("Error: No valid sequences found in input file")
        sys.exit(1)
    
    
    
    # Align sequences
    
    
    # No need to check for pairwise2 anymore; Align module is imported above
    
    temp_matrix = np.array(list(sequence_data.values()), dtype=object)
    
    for i in range(len(temp_matrix)):
        seq, ref = temp_matrix[i, 0], temp_matrix[i, 1]
        # Initial alignment with 5X TAPE
        aligned_seq, aligned_ref, score = align_tape(seq, ref)
        # Check for gaps and try different TAPE lengths
        if '-----' in aligned_seq:
            aligned_seq, aligned_ref, score = align_tape(seq, TAPE_4X)
        elif '----------' in aligned_ref:
            aligned_seq, aligned_ref, score = align_tape(seq, TAPE_6X)
        temp_matrix[i, 0:2] = aligned_seq, aligned_ref
        temp_matrix[i, -1] = score
    
    # Use barcode map keys for position counting
    barcode_keys = list(barcode_map.keys())
    position_counts = {i: {barcode: 0 for barcode in barcode_keys} for i in range(1, 6)}
    # Initialize combination counts
    combination_counts = {f"{a}+{b} {i}to{i+1}": 0 for a in barcode_keys for b in barcode_keys for i in range(1, 5)}    
    combination_counts.update({f'{a}+{b} 1to3': 0 for a in barcode_keys for b in barcode_keys if a != b})
    for sequence_data in temp_matrix:
        # Find edit positions
        edit_matches = list(EDIT_PATTERN.finditer(sequence_data[1]))
        edit_indices = [(m.start() + 2, m.end() - 3) for m in edit_matches]
        # Extract inserted sequences using min_edit function, dynamically using barcode length
        barcode_length = len(barcode_keys[0])
        inserted_sequences = [
            min_edit(barcode_map, sequence_data[0][k[0]:k[0]+barcode_length])[0]
            for k in edit_indices
        ]
        # Map to positions
        position_indices = []
        for index_pair in edit_indices:
            if index_pair[0] in INSERT_POSITIONS:
                position_indices.append(INSERT_POSITIONS[index_pair[0]])
        # Count valid edits
        if edit_indices and sequence_data[-1] > MIN_SCORE_THRESHOLD and position_indices:
            for i, pos in enumerate(position_indices):
                if i < len(inserted_sequences):
                    barcode = inserted_sequences[i]
                    # Original script behavior: no explicit validation, just try to count
                    # This will only work if barcode is 'GTT' or 'ACA' (existing keys)
                    try:
                        position_counts[pos][barcode] += 1
                        # Count combinations
                        try:
                            if i + 1 < len(inserted_sequences):
                                next_barcode = inserted_sequences[i + 1]
                                next_pos = position_indices[i + 1]
                                combo_key = f"{barcode}+{next_barcode} {pos}to{next_pos}"
                                if combo_key in combination_counts:
                                    combination_counts[combo_key] += 1
                            
                            if i + 2 < len(inserted_sequences):
                                next2_barcode = inserted_sequences[i + 2]
                                next2_pos = position_indices[i + 2]
                                combo_key = f"{barcode}+{next2_barcode} {pos}to{next2_pos}"
                                if combo_key in combination_counts:
                                    combination_counts[combo_key] += 1
                        except (IndexError, KeyError):
                            continue
                    except KeyError:
                        # Barcode not found in position_counts (not GTT or ACA), skip like original
                        continue
    
    # Create output DataFrame
    
    
    # Position-based data
    position_data = []
    for pos in position_counts:
        for barcode in position_counts[pos]:
            count = position_counts[pos][barcode]
            position_data.append([barcode, pos, count])
    
    # Combination data
    combo_data = []
    for combo_key, count in combination_counts.items():
        parts = combo_key.split(' ')
        barcode = parts[0]
        position = parts[1]
        combo_data.append([barcode, position, count])
    
    # Combine all data
    all_data = position_data + combo_data
    result_df = pd.DataFrame(all_data, columns=['barcode', 'Position', 'count'])
    result_df['ratio'] = (result_df['count'] / len(temp_matrix)) * 100
    
    # Save results
    output_file = args.output or (args.workdir / f"{args.sample}_bc_count.csv")
    result_df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()
