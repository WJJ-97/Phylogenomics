#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# FASTA File Filter v2024.12.1
# Author: WJJ

import os
import argparse
import glob
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
from time import perf_counter
import sys

DEFAULT_THREADS = min(os.cpu_count() or 4, 8)  # Cap at 8 threads by default

def count_protein_sequences(file_path):
    """
    Count protein sequences in a FASTA file (lines starting with '>').
    
    Args:
        file_path (str): Path to FASTA file
        
    Returns:
        tuple: (file_path: str, count: int, error: Exception)
    """
    try:
        with open(file_path, 'r') as f:
            return (file_path, sum(1 for line in f if line.startswith('>')), None)
    except Exception as e:
        return (file_path, 0, e)

def copy_validated_files(file_list, dest_dir, max_workers):
    """
    Parallel file copy operation with error handling.
    
    Args:
        file_list (list): List of source file paths
        dest_dir (str): Destination directory path
        max_workers (int): Maximum parallel threads
        
    Returns:
        tuple: (success_count: int, error_list: list)
    """
    os.makedirs(dest_dir, exist_ok=True)
    success = 0
    errors = []
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        copy_task = partial(shutil.copy, dst=dest_dir)
        future_map = {executor.submit(copy_task, src): src for src in file_list}
        
        for future in as_completed(future_map):
            src_path = future_map[future]
            try:
                future.result()
                success += 1
            except Exception as e:
                errors.append((src_path, str(e)))
    
    return success, errors

def validate_arguments(source_dir, target_ext, threshold):
    """Perform parameter validation checks"""
    if not os.path.isdir(source_dir):
        raise ValueError(f"Source directory not found: {source_dir}")
    if not target_ext.isalnum():
        raise ValueError(f"Invalid file extension: {target_ext}")
    if threshold < 1:
        raise ValueError("Threshold must be ≥ 1")

def main():
    """Main workflow execution"""
    parser = argparse.ArgumentParser(
        description="Filter FASTA files by protein sequence count",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("source_dir", help="Directory containing FASTA files")
    parser.add_argument("target_ext", help="File extension (without dot)")
    parser.add_argument("threshold", type=int, help="Minimum sequence count")
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS,
                       help="Maximum parallel processing threads")
    
    args = parser.parse_args()
    
    try:
        validate_arguments(args.source_dir, args.target_ext, args.threshold)
    except ValueError as e:
        sys.exit(f"Parameter error: {e}")

    file_pattern = os.path.join(args.source_dir, f"*.{args.target_ext}")
    fasta_files = glob.glob(file_pattern)
    
    if not fasta_files:
        sys.exit(f"No *.{args.target_ext} files found in {args.source_dir}")
    
    print(f"Analyzing {len(fasta_files)} FASTA files...")
    start_time = perf_counter()
    
    valid_files = []
    error_list = []
    sequence_counts = []

    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        results = executor.map(count_protein_sequences, fasta_files)
        
        for result in results:
            file_path, count, error = result
            if error:
                error_list.append((file_path, error))
            elif count >= args.threshold:
                valid_files.append(file_path)
                sequence_counts.append(count)
    
    # Generate statistics
    analysis_time = perf_counter() - start_time
    valid_count = len(valid_files)
    validation_ratio = valid_count / len(fasta_files)
    
    print(f"\nAnalysis completed in {analysis_time:.2f}s")
    print(f"Files meeting threshold ({args.threshold}+ sequences): {valid_count}")
    print(f"Validation rate: {validation_ratio:.1%}")
    
    if error_list:
        print("\nProcessing errors encountered:")
        for path, err in error_list[:3]:  # Show first 3 errors
            print(f" - {os.path.basename(path)}: {err}")
        if len(error_list) > 3:
            print(f" (...{len(error_list)-3} additional errors)")

    if valid_count > 0:
        if input("\nCopy validated files? [y/N]: ").lower() == 'y':
            dest_dir = input("Destination directory path: ").strip()
            if dest_dir:
                copy_start = perf_counter()
                copied, copy_errors = copy_validated_files(
                    valid_files, dest_dir, args.threads
                )
                copy_time = perf_counter() - copy_start
                
                print(f"\nCopied {copied}/{valid_count} files in {copy_time:.2f}s")
                if copy_errors:
                    print(f"Copy errors: {len(copy_errors)}")
            else:
                print("No destination directory specified")
    
    total_time = perf_counter() - start_time
    print(f"\nTotal execution time: {total_time:.2f} seconds")

if __name__ == "__main__":
    main()
