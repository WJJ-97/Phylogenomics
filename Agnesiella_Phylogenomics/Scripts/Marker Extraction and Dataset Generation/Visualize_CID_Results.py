#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# CID Distance Analyzer v2024.12.1
# Author: WJJ

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List

# Configure visual settings
plt.style.use('seaborn')
sns.set_palette('husl')
COLOR_SCATTER = '#4B8BBE'
COLOR_IQR = '#FF6347'
COLOR_MAD = '#FFA500'

def parse_cid_matrix(file_path: str) -> pd.DataFrame:
    """Parse CID distance matrix from CSV file.
    
    Args:
        file_path: Path to CSV file containing distance matrix
        
    Returns:
        pd.DataFrame: Distance matrix with string-type indices
        
    Raises:
        FileNotFoundError: If input file is not found
    """
    try:
        df = pd.read_csv(file_path, index_col=0)
        return df.astype(float)
    except FileNotFoundError:
        raise FileNotFoundError(f"Distance matrix file not found: {file_path}")

def validate_tree_names(matrix: pd.DataFrame, tree_names: List[str]) -> None:
    """Validate tree names against matrix dimensions.
    
    Args:
        matrix: CID distance matrix DataFrame
        tree_names: List of tree names to validate
        
    Raises:
        ValueError: If count mismatch occurs
    """
    if len(matrix) != len(tree_names):
        raise ValueError(
            f"Dimension mismatch: {len(matrix)} trees in matrix vs "
            f"{len(tree_names)} names provided"
        )

def compute_statistics(distances: pd.Series) -> Dict[str, float]:
    """Calculate distribution statistics for CID distances.
    
    Args:
        distances: Series of average distances
        
    Returns:
        Dictionary containing statistical measures
    """
    stats = {
        'mean': distances.mean(),
        'median': distances.median(),
        'q1': np.percentile(distances, 25),
        'q3': np.percentile(distances, 75),
        'mad': np.median(np.abs(distances - distances.median()))
    }
    stats['iqr'] = stats['q3'] - stats['q1']
    return stats

def generate_report(output_file: str, stats: Dict, thresholds: Dict, 
                    distances: pd.Series) -> None:
    """Generate statistical report with thresholds.
    
    Args:
        output_file: Path for report output
        stats: Statistical measures dictionary
        thresholds: Threshold values dictionary
        distances: Sorted distance values
    """
    with Path(output_file).open('w') as f:
        f.write("CID Distance Analysis Report\n")
        f.write("============================\n\n")
        f.write(f"Average Distance: {stats['mean']:.4f}\n")
        f.write(f"Median Distance:  {stats['median']:.4f}\n")
        f.write(f"IQR Threshold:    {thresholds['iqr']:.4f} (Median + 1.5*IQR)\n")
        f.write(f"MAD Threshold:    {thresholds['mad']:.4f} (Median + 3.5*MAD)\n\n")
        f.write("Per-Tree Averages:\n")
        f.write(distances.to_string(float_format="%.4f"))

def visualize_distances(distances: pd.Series, stats: Dict, 
                        thresholds: Dict) -> plt.Figure:
    """Create visualization of CID distance distribution.
    
    Args:
        distances: Sorted average distances
        stats: Statistical measures
        thresholds: Threshold values
        
    Returns:
        matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Scatter plot configuration
    sns.scatterplot(
        x=range(len(distances)),
        y=distances,
        color=COLOR_SCATTER,
        alpha=0.6,
        s=40,
        ax=ax
    )
    
    # Threshold lines
    ax.axhline(thresholds['iqr'], color=COLOR_IQR, 
               linestyle='--', label='IQR Threshold')
    ax.axhline(thresholds['mad'], color=COLOR_MAD,
               linestyle='--', label='MAD Threshold')
    
    # Style configuration
    ax.set_title('CID Distance Distribution Analysis', pad=20)
    ax.set_xlabel('Tree Index', labelpad=15)
    ax.set_ylabel('Average CID Distance', labelpad=15)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', frameon=True)
    
    # Dynamic y-axis limits
    buffer = max(distances.std(), 0.1) * 1.5
    ax.set_ylim(
        max(distances.min() - buffer, 0),
        distances.max() + buffer
    )
    
    return fig

def main():
    """Main execution routine."""
    parser = argparse.ArgumentParser(
        description='CID Distance Analysis and Visualization',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('matrix_file', type=str, 
                       help='Path to CID distance matrix CSV')
    parser.add_argument('tree_file', type=str,
                       help='Path to tree names list')
    
    args = parser.parse_args()
    
    # Data processing
    matrix = parse_cid_matrix(args.matrix_file)
    tree_names = Path(args.tree_file).read_text().splitlines()
    validate_tree_names(matrix, tree_names)
    
    # Add tree names and calculate averages
    matrix.index = tree_names
    avg_distances = matrix.mean(axis=1).sort_values()
    
    # Statistical analysis
    stats = compute_statistics(avg_distances)
    thresholds = {
        'iqr': stats['median'] + 1.5 * stats['iqr'],
        'mad': stats['median'] + 3.5 * stats['mad']
    }
    
    # Generate outputs
    generate_report('4-CID_distance_report.txt', stats, thresholds, avg_distances)
    
    # Visualization
    fig = visualize_distances(avg_distances, stats, thresholds)
    fig.savefig('4-CID_distance_trend.png', dpi=300, bbox_inches='tight')
    
    # Threshold outputs
    for threshold_type in ['iqr', 'mad']:
        outliers = avg_distances[avg_distances > thresholds[threshold_type]]
        outliers.to_csv(f'5-Exceed_{threshold_type.upper()}_CID_trees.txt',
                       header=False, float_format="%.4f")

if __name__ == '__main__':
    main()
