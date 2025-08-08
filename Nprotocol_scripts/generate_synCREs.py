#!/usr/bin/env python3
"""
Generate Synthetic CREs (synCREs) for ENGRAM
Author: Will Chen

This script generates synthetic cis-regulatory elements (CREs) with specified motifs,
repeats, spacing, and targeting system (DTT or HEK3). Outputs a tab-separated table
ready for oligonucleotide ordering.
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import List, Tuple


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(complement)[::-1]


def gen_synCRE(sequence: str, motif: str, repeat: int, spacing: int, dist: int) -> str:
    """
    Generate a synthetic enhancer sequence by repeating a motif with specified spacing 
    and distance from the minimal promoter.

    Args:
        sequence (str): The base DNA sequence.
        motif (str): The motif to repeat.
        repeat (int): Number of motif repeats.
        spacing (int): Number of bases between motifs.
        dist (int): Distance from the end of the last motif to the minimal promoter.

    Returns:
        str: The constructed synthetic enhancer sequence.
    """
    motif_len = len(motif)
    total_motif_len = motif_len * repeat
    total_spacing_len = spacing * (repeat - 1)
    enhancer_len = total_motif_len + total_spacing_len

    start = len(sequence) - enhancer_len - dist
    parts = [sequence[:start]]
    
    # Add motif and spacing in order
    for idx in range(repeat):
        parts.append(motif)
        if idx < repeat - 1:
            # Calculate the start and end for the spacing region
            spacing_start = start + (motif_len + spacing) * idx + motif_len
            spacing_end = spacing_start + spacing
            parts.append(sequence[spacing_start:spacing_end])

    # Add the minimal promoter region if dist > 0
    if dist == 0:
        return ''.join(parts)
    else:
        return ''.join(parts) + sequence[-dist:]


def generate_oligo(cre_sequence: str, system: str) -> str:
    """
    Generate the full oligonucleotide sequence with the appropriate cloning components.

    Args:
        cre_sequence (str): The CRE sequence to insert
        system (str): Either 'DTT' or 'HEK3' for the targeting system

    Returns:
        str: Complete oligonucleotide sequence ready for ordering
    """
    # Sequence components for ENGRAM cloning
    BsaI_fwd = 'atactacGGTCTCagaac'
    BsaI_rev = 'actgcGAGACCgtaatgc'
    BsmBI_filler = 'CACTAgagacgattaaATGGACAGCAtgcgtctcTGTCC'
    
    # pegRNA sequences
    DTT_peg = 'CACCATCATCCNNNNNCGTGCTCACCATC'  # HA + typewriter_key + symbol (5N) + PBS
    HEK3_peg = 'TCTGCCATCANNNNNCGTGCTCAGTCTG'  # HA + symbol (5N) + PBS
    
    if system.upper() == 'DTT':
        peg = DTT_peg
    elif system.upper() == 'HEK3':
        peg = HEK3_peg
    else:
        raise ValueError(f"Unknown system: {system}. Must be 'DTT' or 'HEK3'")
    
    return BsaI_fwd + cre_sequence + BsmBI_filler + peg + BsaI_rev


def parse_motifs_file(file_path: Path) -> List[Tuple[str, str]]:
    """
    Parse a motifs file containing motif names and sequences.
    Expected format: tab-separated with columns 'name' and 'sequence'
    
    Args:
        file_path (Path): Path to the motifs file
        
    Returns:
        List[Tuple[str, str]]: List of (motif_name, motif_sequence) tuples
    """
    motifs = []
    try:
        with open(file_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if 'name' in row and 'sequence' in row:
                    motifs.append((row['name'].strip(), row['sequence'].strip().upper()))
                else:
                    # Try first two columns if headers are different
                    cols = list(row.values())
                    if len(cols) >= 2:
                        motifs.append((cols[0].strip(), cols[1].strip().upper()))
    except Exception as e:
        print(f"Error reading motifs file: {e}")
        sys.exit(1)
    
    return motifs


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate synthetic CREs (synCREs) for ENGRAM experiments",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Motif specification
    motif_group = parser.add_mutually_exclusive_group(required=True)
    motif_group.add_argument(
        '-m', '--motif',
        help='Single motif sequence (DNA)'
    )
    motif_group.add_argument(
        '-f', '--motifs-file',
        type=Path,
        help='Tab-separated file with motif names and sequences'
    )
    
    # Motif name (only used with single motif)
    parser.add_argument(
        '-n', '--motif-name',
        default='motif',
        help='Name for the motif (used with -m/--motif)'
    )
    
    # CRE parameters
    parser.add_argument(
        '-r', '--repeats',
        type=str,
        default='6',
        help='Number of repeats (single number or range like "1-7")'
    )
    
    parser.add_argument(
        '-s', '--spacing',
        type=int,
        default=15,
        help='Number of bases between motifs'
    )
    
    parser.add_argument(
        '-d', '--distance',
        type=int,
        default=0,
        help='Distance from last motif to minimal promoter'
    )
    
    # System selection
    parser.add_argument(
        '--system',
        choices=['DTT', 'HEK3'],
        default='HEK3',
        help='Targeting system to use'
    )
    
    # Background sequence
    parser.add_argument(
        '--background-seq',
        default='GCAAACAGTCCCCATGTTCACATTAGGTTCCCAATCAGTATTTCCCATCAAGGGAGAAGCTTAGGCTGGGGAATCCTGCATACATAATTCATACCCATTATCTCAGGCTTGCTTTCTTCAACTTAAGTAGATAGGTGTATATAAAGAGACACAGCATGAAGAAGAAAGATAATAAAGATGGGCCCTGGAGAGTGACCTCTGGTGACCACAAACCTGACCTC',
        help='Background sequence for CRE generation'
    )
    
    # Output options
    parser.add_argument(
        '-o', '--output',
        type=Path,
        help='Output file path (default: stdout)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )
    
    return parser.parse_args()


def parse_repeat_range(repeat_str: str) -> List[int]:
    """
    Parse repeat specification string.
    Args:
        repeat_str: String like "3" or "1-7" or "1,3,5"
    Returns:
        List of repeat numbers
    """
    if ',' in repeat_str:
        # Comma-separated list
        return [int(x.strip()) for x in repeat_str.split(',')]
    elif '-' in repeat_str:
        # Range
        start, end = repeat_str.split('-', 1)
        return list(range(int(start.strip()), int(end.strip()) + 1))
    else:
        # Single number
        return [int(repeat_str.strip())]


def main():
    """Main execution function."""
    args = parse_arguments()
    if args.verbose:
        print("Generating synthetic CREs...", file=sys.stderr)
    # Get motifs
    if args.motif:
        motifs = [(args.motif_name, args.motif.upper())]
    else:
        motifs = parse_motifs_file(args.motifs_file)
    # Parse repeat range
    try:
        repeat_list = parse_repeat_range(args.repeats)
    except ValueError as e:
        print(f"Error parsing repeat range '{args.repeats}': {e}", file=sys.stderr)
        sys.exit(1)
    # Generate synCREs
    results = []
    for motif_name, motif_seq in motifs:
        if args.verbose:
            print(f"Processing motif: {motif_name}", file=sys.stderr)
        for repeat_count in repeat_list:
            # Generate the synCRE
            syn_cre = gen_synCRE(
                args.background_seq,
                motif_seq,
                repeat_count,
                args.spacing,
                args.distance
            )
            # Generate full oligo
            oligo = generate_oligo(syn_cre, args.system)
            # Create name
            cre_name = f"{motif_name}_{repeat_count}x_{args.system}"
            # Store result (always include sequence)
            results.append([cre_name, motif_seq, str(repeat_count), args.system, oligo])
    
    # Write output
    if args.output:
        output_file = open(args.output, 'w', newline='')
        if args.verbose:
            print(f"Writing output to: {args.output}", file=sys.stderr)
    else:
        output_file = sys.stdout
    
    try:
        writer = csv.writer(output_file, delimiter='\t')
        # Write header (always include sequence)
        writer.writerow(['CRE_Name', 'Motif_Sequence', 'Repeats', 'System', 'ready_to_order_oligo'])
        # Write data
        writer.writerows(results)
        if args.verbose:
            print(f"Generated {len(results)} synCRE designs", file=sys.stderr)
    finally:
        if args.output:
            output_file.close()

if __name__ == "__main__":
    main()
