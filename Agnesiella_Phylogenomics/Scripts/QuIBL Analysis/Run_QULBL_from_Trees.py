#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Author: WJJ v2024.12.1

import os
import logging
import argparse
from itertools import combinations
from ete3 import Tree
import multiprocessing
import subprocess
import shutil
import time

# ==============================================
# Performance Tracking Infrastructure
# ==============================================
class Timer:
    """Context manager for code block timing"""
    def __init__(self, name):
        self.name = name
        self.start = None
        self.end = None

    def __enter__(self):
        self.start = time.time()
        logging.info(f"START: {self.name}")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end = time.time()
        elapsed = self.end - self.start
        logging.info(f"COMPLETE: {self.name} ({elapsed:.2f}s)")

def track_time(step_name):
    """Decorator for function-level timing"""
    def decorator(func):
        def wrapper(*args, **kwargs):
            start = time.time()
            result = func(*args, **kwargs)
            elapsed = time.time() - start
            logging.info(f"[TIMED] {step_name}: {elapsed:.2f}s")
            return result
        return wrapper
    return decorator

# ==============================================
# Core Processing Functions
# ==============================================
@track_time("Species Combination Generation")
def generate_species_combinations(species_list, output_file, outgroup):
    """Generate all 3-species combinations with outgroup
    
    Args:
        species_list: List of base species
        output_file: Output file path
        outgroup: Outgroup species name
    
    Returns:
        int: Number of generated combinations
    """
    combs_3 = list(combinations(species_list, 3))
    combs_4 = [comb + (outgroup,) for comb in combs_3]
    
    # Batch write for I/O efficiency
    with open(output_file, 'w') as f:
        f.write('\n'.join(' '.join(comb) for comb in combs_4))
    
    return len(combs_4)

@track_time("Tree Pruning Operation")
def process_line(args):
    """Process individual tree line with pruning
    
    Args:
        args: Tuple containing (line, line_num, tree_path, output_dir)
    
    Returns:
        int: Number of successfully processed trees
    """
    line, line_number, tree_file_path, output_dir = args
    try:
        subtree_taxa = line.strip().split()
        pruned_trees = []
        
        # Bulk tree loading
        with open(tree_file_path, 'r') as tree_file:
            trees = [t.strip() for t in tree_file if t.strip()]
            
        # Tree processing loop
        for tree_line in trees:
            try:
                tree = Tree(tree_line)
                tree.prune(subtree_taxa, preserve_branch_length=True)
                pruned_trees.append(tree.write())
            except Exception as e:
                logging.debug(f"Pruning error line {line_number}: {str(e)}")
                continue
        
        # Batch output writing
        if pruned_trees:
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, f"out_subtree_{line_number}.txt")
            with open(output_path, 'w') as f:
                f.write('\n'.join(pruned_trees))
                
        return len(pruned_trees)
        
    except Exception as e:
        logging.error(f"Critical error line {line_number}: {str(e)}")
        return 0

@track_time("Configuration Generation")
def generate_input_config_files(pruned_tree_dir, output_config_dir, outgroup):
    """Generate QuIBL configuration files
    
    Args:
        pruned_tree_dir: Directory with pruned trees
        output_config_dir: Output directory for configs
        outgroup: Outgroup species name
    """
    os.makedirs(output_config_dir, exist_ok=True)
    
    for file_name in os.listdir(pruned_tree_dir):
        if not file_name.startswith('out_subtree'):
            continue
            
        pruned_tree_file = os.path.join(pruned_tree_dir, file_name)
        config_file_path = os.path.join(output_config_dir, f"run_{file_name}")
        
        config_content = f"""[Input]
treefile: {pruned_tree_file}
numdistributions: 2
likelihoodthresh: 0.01
numsteps: 50
gradascentscalar: 0.5
totaloutgroup: {outgroup}
multiproc: True
maxcores: 1

[Output]
OutputPath: {os.path.join(output_config_dir, f"{file_name}.csv")}
"""
        with open(config_file_path, 'w') as config_file:
            config_file.write(config_content)

# ==============================================
# Execution Management
# ==============================================
@track_time("QuIBL Execution")
def run_quibl(filename, quibl_script, finished_dir):
    """Execute QuIBL analysis for a config file
    
    Args:
        filename: Config file path
        quibl_script: Path to QuIBL.py
        finished_dir: Directory for processed files
    """
    try:
        subprocess.call(['python', quibl_script, filename])
        logging.info(f"Processed: {filename}")
        
        os.makedirs(finished_dir, exist_ok=True)
        shutil.move(filename, finished_dir)
        logging.info(f"Moved: {filename} -> {finished_dir}")

    except Exception as e:
        logging.error(f"Processing failed: {filename} - {str(e)}")

# ==============================================
# Main Control Flow
# ==============================================
def main():
    """Main pipeline controller"""
    global_start = time.time()
    parser = argparse.ArgumentParser(description='Optimized phylogenetic pipeline')
    
    # Argument definitions
    parser.add_argument('--species_list_file', help='Species list file path')
    parser.add_argument('--tree_file_path', help='Input tree file path')
    parser.add_argument('--outgroup', required=True, help='Outgroup species')
    parser.add_argument('--pruned_tree_dir', help='Pruned trees directory')
    parser.add_argument('--output_path_base', help='Base output directory')
    parser.add_argument('--quibl_script', help='QuIBL script path')
    parser.add_argument('--finished_dir', help='Processed files directory')
    parser.add_argument('--pool_size', type=int, default=4, help='Process pool size')
    parser.add_argument('--steps', nargs='+', required=True, 
                       choices=['generate_combinations', 'prune_trees', 
                               'generate_config', 'run_quibl'],
                       help='Pipeline steps to execute')

    args = parser.parse_args()

    # Logging configuration
    logging.basicConfig(
        filename='pipeline_analysis.log',
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.INFO,
        filemode='a'
    )

    # Execute requested steps
    for step in args.steps:
        with Timer(f"STEP: {step.upper()}"):
            if step == 'generate_combinations':
                if not args.species_list_file:
                    raise ValueError("Missing species list file")
                
                with open(args.species_list_file) as f:
                    species_list = [line.strip() for line in f if line.strip()]
                
                comb_count = generate_species_combinations(
                    species_list, 
                    'temp_combinations.txt', 
                    args.outgroup
                )
                logging.info(f"Generated {comb_count} combinations")

            elif step == 'prune_trees':
                if not all([args.tree_file_path, args.pruned_tree_dir]):
                    raise ValueError("Missing required arguments for tree pruning")
                
                cpu_count = multiprocessing.cpu_count()
                pool_size = min(args.pool_size, cpu_count)
                
                with multiprocessing.Pool(pool_size) as pool:
                    tasks = []
                    with open('temp_combinations.txt') as taxa_file:
                        for line_num, line in enumerate(taxa_file, 1):
                            tasks.append((
                                line, 
                                line_num, 
                                os.path.abspath(args.tree_file_path), 
                                os.path.abspath(args.pruned_tree_dir)
                            ))
                    
                    # Process in chunks for memory efficiency
                    chunk_size = 100
                    results = []
                    for i in range(0, len(tasks), chunk_size):
                        chunk = tasks[i:i+chunk_size]
                        results.extend(pool.map(process_line, chunk))
                        logging.info(f"Processed chunk {i//chunk_size} ({sum(results)} trees)")
                
                logging.info(f"Total pruned trees: {sum(results)}")

            elif step == 'generate_config':
                generate_input_config_files(
                    os.path.abspath(args.pruned_tree_dir),
                    os.path.abspath(args.output_path_base),
                    args.outgroup
                )

            elif step == 'run_quibl':
                input_files = []
                for root, _, files in os.walk(args.output_path_base):
                    input_files.extend(
                        os.path.join(root, f) 
                        for f in files 
                        if f.startswith('run_out_subtree')
                    )

                with multiprocessing.Pool(args.pool_size) as pool:
                    pool.starmap(run_quibl, [
                        (f, args.quibl_script, args.finished_dir) 
                        for f in input_files
                    ])

    # Final reporting
    total_time = time.time() - global_start
    logging.info(f"TOTAL PIPELINE TIME: {total_time:.2f}s")
    print(f"\nPipeline completed in {total_time//3600:.0f}h "
          f"{(total_time%3600)//60:.0f}m {total_time%60:.2f}s")

if __name__ == '__main__':
    main()
