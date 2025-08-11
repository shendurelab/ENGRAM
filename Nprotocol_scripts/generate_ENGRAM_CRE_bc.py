#!/usr/bin/env python3
"""
Generate Synthetic CREs (synCREs) for ENGRAM
Author: Will Chen

This script generates synthetic cis-regulatory elements (CREs) with specified motifs,
repeats, spacing, and targeting system (DTT or HEK3). It can also process custom CREs
provided by the user. All oligonucleotides are assigned unique barcodes from predefined
barcode libraries. Outputs a tab-separated table ready for oligonucleotide ordering.
"""

import argparse
import csv
import sys
import random
from pathlib import Path
from typing import List, Tuple, Dict, Set, Optional


def load_barcodes(barcode_file: Path) -> List[str]:
    """
    Load barcodes from a CSV file.
    
    Args:
        barcode_file (Path): Path to the barcode CSV file
        
    Returns:
        List[str]: List of barcode sequences
    """
    barcodes = []
    try:
        with open(barcode_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if 'pBC_Seq' in row:
                    barcodes.append(row['pBC_Seq'].strip())
                else:
                    # Try first column if header is different
                    first_col = list(row.values())[0]
                    if first_col.strip():
                        barcodes.append(first_col.strip())
    except Exception as e:
        print(f"Error reading barcode file: {e}", file=sys.stderr)
        sys.exit(1)
    
    return barcodes


def assign_unique_barcodes(num_needed: int, barcodes: List[str], 
                          used_barcodes: Optional[Set[str]] = None) -> List[str]:
    """
    Assign unique barcodes from the available pool.
    
    Args:
        num_needed (int): Number of unique barcodes needed
        barcodes (List[str]): Available barcode pool
        used_barcodes (Set[str], optional): Previously used barcodes to avoid
        
    Returns:
        List[str]: List of unique assigned barcodes
        
    Raises:
        ValueError: If not enough unique barcodes are available
    """
    if used_barcodes is None:
        used_barcodes = set()
    
    available_barcodes = [bc for bc in barcodes if bc not in used_barcodes]
    
    if len(available_barcodes) < num_needed:
        raise ValueError(f"Not enough unique barcodes available. Need {num_needed}, "
                        f"but only {len(available_barcodes)} available after excluding used ones.")
    
    return random.sample(available_barcodes, num_needed)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(complement)[::-1]


def parse_custom_cres_file(file_path: Path) -> List[Tuple[str, str]]:
    """
    Parse a custom CREs file containing CRE names and sequences.
    Expected format: tab-separated with column 1 = name, column 2 = sequence
    
    Args:
        file_path (Path): Path to the custom CREs file
        
    Returns:
        List[Tuple[str, str]]: List of (cre_name, cre_sequence) tuples
    """
    custom_cres = []
    try:
        with open(file_path, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader)  # Skip header row
            for row in reader:
                if len(row) >= 2:
                    name = row[0].strip()
                    sequence = row[1].strip().upper()
                    custom_cres.append((name, sequence))
    except Exception as e:
        print(f"Error reading custom CREs file: {e}", file=sys.stderr)
        sys.exit(1)
    
    return custom_cres


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


def generate_oligo(cre_sequence: str, system: str, barcode: Optional[str] = None) -> str:
    """
    Generate the full oligonucleotide sequence with the appropriate cloning components.

    Args:
        cre_sequence (str): The CRE sequence to insert
        system (str): Either 'DTT' or 'HEK3' for the targeting system
        barcode (str, optional): Barcode sequence to replace NNNNN placeholders

    Returns:
        str: Complete oligonucleotide sequence ready for ordering
    """
    # Sequence components for ENGRAM cloning
    BsaI_fwd = 'atactacGGTCTCagaac'
    BsaI_rev = 'actgcGAGACCgtaatgc'
    BsmBI_filler = 'CACTAgagacgattaaATGGACAGCAtgcgtctcTGTCC'
    
    # pegRNA sequences with barcode placeholders
    # Use placeholder length based on barcode length if provided, otherwise default to 5N
    if barcode:
        barcode_length = len(barcode)
        n_placeholder = 'N' * barcode_length
    else:
        n_placeholder = 'NNNNN'  # Default 5N
    
    DTT_peg = f'CACCATCATCC{n_placeholder}CGTGCTCACCATC'  # HA + typewriter_key + symbol + PBS
    HEK3_peg = f'TCTGCCATCA{n_placeholder}CGTGCTCAGTCTG'  # HA + symbol + PBS
    
    if system.upper() == 'DTT':
        peg = DTT_peg
    elif system.upper() == 'HEK3':
        peg = HEK3_peg
    else:
        raise ValueError(f"Unknown system: {system}. Must be 'DTT' or 'HEK3'")
    
    # Replace N placeholders with barcode
    if barcode:
        peg = peg.replace(n_placeholder, barcode)
    
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
        description="Generate synthetic CREs (synCREs) for ENGRAM experiments or process custom CREs with barcode assignment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input mode selection
    input_group = parser.add_mutually_exclusive_group(required=True)
    
    # Motif specification (for synCRE generation)
    input_group.add_argument(
        '-m', '--motif',
        help='Single motif sequence (DNA) for synCRE generation'
    )
    input_group.add_argument(
        '-f', '--motifs-file',
        type=Path,
        help='Tab-separated file with motif names and sequences for synCRE generation'
    )
    
    # Custom CREs specification
    input_group.add_argument(
        '-c', '--custom-cres-file',
        type=Path,
        help='Tab-separated file with custom CRE names and sequences (200-300bp recommended)'
    )
    
    # Motif name (only used with single motif)
    parser.add_argument(
        '-n', '--motif-name',
        default='motif',
        help='Name for the motif (used with -m/--motif)'
    )
    
    # CRE parameters (only for synCRE generation)
    parser.add_argument(
        '-r', '--repeats',
        type=str,
        default='6',
        help='Number of repeats for synCRE generation (single number, range like "1-7", or comma-separated like "1,3,5")'
    )
    
    parser.add_argument(
        '-s', '--spacing',
        type=int,
        default=15,
        help='Number of bases between motifs in synCRE generation'
    )
    
    parser.add_argument(
        '-d', '--distance',
        type=int,
        default=0,
        help='Distance from last motif to minimal promoter in synCRE generation'
    )
    
    # System selection
    parser.add_argument(
        '--system',
        choices=['DTT', 'HEK3'],
        default='HEK3',
        help='Targeting system to use'
    )
    
    # Barcode options
    parser.add_argument(
        '--barcode-file',
        type=Path,
        help='CSV file with barcodes (auto-selects 5N or 8N file based on number of CREs)'
    )
    
    # Background sequence (only for synCRE generation)
    parser.add_argument(
        '--background-seq',
        default='GCAAACAGTCCCCATGTTCACATTAGGTTCCCAATCAGTATTTCCCATCAAGGGAGAAGCTTAGGCTGGGGAATCCTGCATACATAATTCATACCCATTATCTCAGGCTTGCTTTCTTCAACTTAAGTAGATAGGTGTATATAAAGAGACACAGCATGAAGAAGAAAGATAATAAAGATGGGCCCTGGAGAGTGACCTCTGGTGACCACAAACCTGACCTC',
        help='Background sequence for synCRE generation'
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
    
    # Always set random seed to 42 for reproducibility
    random.seed(42)
    
    if args.verbose:
        print("Processing ENGRAM CREs with barcode assignment...", file=sys.stderr)
    
    # First, collect all CREs to determine how many we need
    cres = []  # List of (name, sequence) tuples
    motifs = []  # Initialize for later use in results processing
    
    if args.custom_cres_file:
        # Custom CRE mode
        if args.verbose:
            print(f"Processing custom CREs from: {args.custom_cres_file}", file=sys.stderr)
        cres = parse_custom_cres_file(args.custom_cres_file)
        
    elif args.motif or args.motifs_file:
        # synCRE generation mode
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
                # Create name
                cre_name = f"{motif_name}_{repeat_count}x_{args.system}"
                cres.append((cre_name, syn_cre))
    
    if not cres:
        print("Error: No CREs to process", file=sys.stderr)
        sys.exit(1)
    
    # Determine barcode file to use based on number of CREs and system
    if args.barcode_file:
        barcode_file = args.barcode_file
    else:
        script_dir = Path(__file__).parent
        if args.system.upper() == 'DTT':
            # DTT always uses 5N barcodes
            barcode_file = script_dir / 'ENGRAM_5N_good_barcodes.csv'
        else:  # HEK3
            # HEK3 uses 5N by default, but switches to 8N for large numbers of CREs
            if len(cres) > 600:  # Switch to 8N if we need more than 600 barcodes
                barcode_file = script_dir / 'ENGRAM_8N_good_barcodes.csv'
                if args.verbose:
                    print(f"Switching to 8N barcodes due to large number of CREs ({len(cres)})", file=sys.stderr)
            else:
                barcode_file = script_dir / 'ENGRAM_5N_good_barcodes.csv'
    
    # Load barcodes
    if not barcode_file.exists():
        print(f"Error: Barcode file not found: {barcode_file}", file=sys.stderr)
        print("Please specify a barcode file with --barcode-file or ensure the default files exist", file=sys.stderr)
        sys.exit(1)
    
    barcodes = load_barcodes(barcode_file)
    if args.verbose:
        print(f"Loaded {len(barcodes)} barcodes", file=sys.stderr)

    try:
        assigned_barcodes = assign_unique_barcodes(len(cres), barcodes)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Generate oligonucleotides with barcodes
    results = []
    cre_barcode_map = {}  # Track which CRE gets which barcode
    
    for i, ((cre_name, cre_sequence), barcode) in enumerate(zip(cres, assigned_barcodes)):
        # Generate full oligo with barcode
        oligo = generate_oligo(cre_sequence, args.system, barcode)
        
        # Store the CRE-barcode association
        cre_barcode_map[cre_name] = barcode
        
        # For synCREs, extract repeat count from the name
        if args.motif or args.motifs_file:
            repeats = cre_name.split('_')[1]  # Get the repeat part (e.g., "6x")            
            results.append([cre_name, repeats, args.system, barcode, oligo])
        else:
            # Custom CRE mode - simplified output
            results.append([cre_name, args.system, barcode, oligo])
    
    # Write output
    if args.output:
        output_file = open(args.output, 'w', newline='')
        if args.verbose:
            print(f"Writing output to: {args.output}", file=sys.stderr)
    else:
        output_file = sys.stdout
    
    try:
        writer = csv.writer(output_file, delimiter='\t')
        # Write header
        if args.custom_cres_file:
            writer.writerow(['CRE_Name', 'System', 'Barcode', 'ready_to_order_oligo'])
        else:
            writer.writerow(['CRE_Name', 'Repeats', 'System', 'Barcode', 'ready_to_order_oligo'])
        
        # Write data
        writer.writerows(results)
        
        if args.verbose:
            print(f"Generated {len(results)} oligonucleotides with unique barcodes", file=sys.stderr)
            print(f"Used {len(set(assigned_barcodes))} unique barcodes", file=sys.stderr)
            
            # Print CRE-barcode mapping summary
            print("\nCRE-Barcode associations:", file=sys.stderr)
            for cre_name, barcode in sorted(cre_barcode_map.items()):
                print(f"  {cre_name} -> {barcode}", file=sys.stderr)
                
    finally:
        if args.output:
            output_file.close()

if __name__ == "__main__":
    main()
