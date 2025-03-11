# -*- coding: utf-8 -*-
# FASTA File Analysis Toolkit v2024.12.1
# Author: WJJ

"""
Features:
1. Multi-threaded processing of FASTA files/directories
2. Comprehensive statistical analysis
3. Sequence filtering and file management
4. Robust error handling and input validation
"""

import os
import argparse
import shutil
from typing import List, Tuple, Union
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Constants
VALID_EXTENSIONS = ('.fas', '.fasta', '.fa')
DEFAULT_THREADS = 4

def is_valid_fasta(path: str) -> bool:
    """Check if a file has valid FASTA extension.
    
    Args:
        path: File path to verify
        
    Returns:
        True if valid FASTA extension, False otherwise
    """
    return path.lower().endswith(VALID_EXTENSIONS)

def calculate_fasta_stats(fasta_file: str) -> Tuple[int, float, int, int, int]:
    """Calculate comprehensive statistics for a FASTA file.
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        Tuple containing:
        - sequence_count: Number of sequences
        - average_length: Average sequence length
        - total_length: Total base pairs
        - min_length: Shortest sequence length
        - max_length: Longest sequence length
        
    Raises:
        FileNotFoundError: If input file doesn't exist
        ValueError: If empty file is provided
    """
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"File not found: {fasta_file}")
    
    lengths = []
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            lengths.append(len(record.seq))
    except Exception as e:
        print(f"Error parsing {fasta_file}: {str(e)}")
        return (0, 0.0, 0, 0, 0)

    if not lengths:
        print(f"Warning: Empty file detected - {fasta_file}")
        return (0, 0.0, 0, 0, 0)

    sequence_count = len(lengths)
    total_length = sum(lengths)
    return (
        sequence_count,
        total_length / sequence_count,
        total_length,
        min(lengths),
        max(lengths)
    )

def parallel_process_files(directory: str) -> List[Tuple[str, tuple]]:
    """Process FASTA files in parallel using ThreadPoolExecutor.
    
    Args:
        directory: Path to directory containing FASTA files
        
    Returns:
        List of tuples (file_path, statistics)
    """
    fasta_files = [entry.path for entry in os.scandir(directory) 
                  if entry.is_file() and is_valid_fasta(entry.name)]
    
    with ThreadPoolExecutor(max_workers=DEFAULT_THREADS) as executor:
        futures = {executor.submit(calculate_fasta_stats, f): f for f in fasta_files}
        return [(f.result(), futures[f]) for f in as_completed(futures)]

def filter_sequences(fasta_file: str, threshold: int) -> List[SeqRecord]:
    """Retrieve sequences exceeding length threshold.
    
    Args:
        fasta_file: Input FASTA file path
        threshold: Minimum sequence length to include
        
    Returns:
        List of filtered SeqRecord objects
    """
    return [record for record in SeqIO.parse(fasta_file, "fasta")
           if len(record.seq) >= threshold]

def safe_write_output(records: List[SeqRecord], output_path: str) -> None:
    """Safely write sequences to output file with directory validation.
    
    Args:
        records: List of SeqRecord objects to write
        output_path: Target output file path
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as handle:
        SeqIO.write(records, handle, "fasta")

def directory_analysis(directory: str) -> None:
    """Analyze all FASTA files in a directory and display summary statistics.
    
    Args:
        directory: Path to directory containing FASTA files
    """
    results = parallel_process_files(directory)
    valid_results = [r for r in results if r[0][0] > 0]
    
    loci_count = len(valid_results)
    total_sites = sum(r[0][2] for r in valid_results)
    lengths = [r[0][2] for r in valid_results]
    
    print(f"\nDirectory Analysis: {directory}")
    print(f"Valid FASTA Files: {loci_count}")
    print(f"Total Base Pairs: {total_sites}")
    print(f"Average Locus Length: {total_sites/loci_count:.2f}" if loci_count else "No valid files")
    print(f"Size Range: [{min(lengths)} - {max(lengths)}]" if lengths else "No data")

def handle_user_choices(input_path: str, is_directory: bool) -> None:
    """Manage interactive user prompts for file operations.
    
    Args:
        input_path: Path being processed
        is_directory: Flag for directory vs file processing
    """
    if is_directory:
        choice = input("Copy files with average length threshold? [y/n]: ").lower()
        if choice == 'y':
            threshold = int(input("Minimum average length: "))
            output_dir = input("Output directory: ")
            copied = copy_qualified_files(input_path, threshold, output_dir)
            print(f"Copied {len(copied)} files to {output_dir}")
    else:
        choice = input("Extract long sequences? [y/n]: ").lower()
        if choice == 'y':
            threshold = int(input("Minimum sequence length: "))
            output_file = input("Output file path: ")
            records = filter_sequences(input_path, threshold)
            safe_write_output(records, output_file)
            print(f"Saved {len(records)} sequences to {output_file}")

def main():
    """Main execution flow with argument parsing."""
    parser = argparse.ArgumentParser(
        description="FASTA File Analysis System",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("inputs", nargs="+", help="FASTA files/directories to process")
    args = parser.parse_args()

    for path in args.inputs:
        if not os.path.exists(path):
            print(f"Invalid path: {path}")
            continue
            
        if os.path.isdir(path):
            directory_analysis(path)
            handle_user_choices(path, is_directory=True)
        else:
            try:
                stats = calculate_fasta_stats(path)
                print(f"\nFile Analysis: {path}")
                print(f"Sequences: {stats[0]}\nAvg Length: {stats[1]:.1f}")
                print(f"Total Length: {stats[2]}\nSize Range: [{stats[3]} - {stats[4]}]")
                handle_user_choices(path, is_directory=False)
            except Exception as e:
                print(f"Error processing {path}: {str(e)}")

if __name__ == "__main__":
    main()
