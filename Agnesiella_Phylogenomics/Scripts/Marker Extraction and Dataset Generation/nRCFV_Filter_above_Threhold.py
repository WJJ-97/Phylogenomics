# nRCFV Filter Script - Identifies outliers using Median Absolute Deviation (MAD) v2024.12.1
# Author: WJJ

import sys
import argparse
import statistics
from typing import List, Tuple

def read_data(filename: str) -> Tuple[List[str], List[float]]:
    """
    Read input file and extract identifiers with corresponding values
    
    Args:
        filename: Path to input file with two columns:
                  Column 1: Identifier (string)
                  Column 2: Numeric value (float)
    
    Returns:
        Tuple containing:
        - List of identifiers
        - List of numeric values
    
    Raises:
        FileNotFoundError: If input file doesn't exist
        ValueError: If invalid data format encountered
    """
    identifiers = []
    values = []
    
    with open(filename, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            
            parts = line.split()
            if len(parts) != 2:
                raise ValueError(f"Invalid column count at line {line_num}: {line}")
            
            identifier, value = parts
            try:
                values.append(float(value))
                identifiers.append(identifier)
            except ValueError:
                raise ValueError(f"Non-numeric value at line {line_num}: {value}")
    
    if len(values) < 4:
        raise ValueError("Insufficient data points (minimum 4 required)")
    
    return identifiers, values

def calculate_threshold(values: List[float]) -> Tuple[float, float, float]:
    """
    Calculate statistical thresholds using Median and MAD
    
    Args:
        values: List of numeric values
    
    Returns:
        Tuple containing:
        - Median
        - MAD
        - Modified Z-score threshold value
    """
    median = statistics.median(values)
    mad = statistics.median([abs(x - median) for x in values])
    threshold_value = median + 3.5 * mad  # Standard modified Z-score cutoff
    return median, mad, threshold_value

def write_outliers(identifiers: List[str], values: List[float], 
                  threshold: float, output_file: str) -> int:
    """
    Write outlier identifiers to output file
    
    Args:
        identifiers: List of original identifiers
        values: Corresponding numeric values
        threshold: Threshold value for outlier detection
        output_file: Output filename
    
    Returns:
        Number of outliers detected
    """
    outlier_count = 0
    with open(output_file, 'w') as f:
        for ident, val in zip(identifiers, values):
            if val >= threshold:
                f.write(f"{ident}\n")
                outlier_count += 1
    return outlier_count

def main():
    """Main execution routine"""
    parser = argparse.ArgumentParser(
        description="Identify outliers using MAD-based modified Z-scores")
    parser.add_argument("input_file", help="Path to input data file")
    parser.add_argument("-o", "--output", default="nRCFV_rm.txt",
                       help="Output filename (default: nRCFV_rm.txt)")
    args = parser.parse_args()

    try:
        # Data ingestion
        identifiers, values = read_data(args.input_file)
        
        # Statistical calculations
        median, mad, threshold = calculate_threshold(values)
        
        # Result output
        outlier_count = write_outliers(identifiers, values, threshold, args.output)
        
        # Summary report
        print(f"Statistical Report:\n"
              f"- Median: {median:.4f}\n"
              f"- MAD: {mad:.4f}\n"
              f"- Threshold Value: {threshold:.4f}\n"
              f"- Outliers Identified: {outlier_count}\n"
              f"Results written to: {args.output}")

    except Exception as e:
        print(f"\nError: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
