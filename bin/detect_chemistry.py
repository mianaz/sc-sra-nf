#!/usr/bin/env python3
"""
Detect single-cell chemistry/platform from FASTQ files
Identifies 10x, Seqwell, GeXscope, and other formats
"""
import sys
import gzip
import os
import re
import json
from pathlib import Path

def read_fastq_headers(fastq_path: str, n_reads: int = 1000) -> list:
    """Read first N headers from FASTQ file"""
    headers = []
    opener = gzip.open if fastq_path.endswith('.gz') else open

    try:
        with opener(fastq_path, 'rt') as f:
            count = 0
            for line in f:
                if line.startswith('@') and count < n_reads:
                    headers.append(line.strip())
                    count += 1
                    # Skip next 3 lines (seq, +, qual)
                    for _ in range(3):
                        next(f, None)
                if count >= n_reads:
                    break
    except Exception as e:
        print(f"Error reading {fastq_path}: {e}", file=sys.stderr)

    return headers

def detect_read_structure(fastq_dir: str) -> dict:
    """Analyze FASTQ files to determine read structure"""
    fastq_files = []
    for ext in ['*.fastq.gz', '*.fq.gz', '*.fastq', '*.fq']:
        fastq_files.extend(Path(fastq_dir).glob(f'**/{ext}'))

    if not fastq_files:
        return {'error': 'No FASTQ files found'}

    # Group files by sample/run
    file_groups = {}
    for f in fastq_files:
        # Extract sample/run identifier
        fname = f.name
        # Remove common suffixes to get base name
        base = re.sub(r'(_S\d+)?(_L\d+)?(_[IR][12])(_\d+)?\.(fastq|fq)(\.gz)?$', '', fname)
        if base not in file_groups:
            file_groups[base] = []
        file_groups[base].append(str(f))

    # Analyze read lengths
    read_info = {}
    for base, files in file_groups.items():
        read_info[base] = {}
        for f in sorted(files):
            fname = os.path.basename(f)
            # Estimate read length
            length = estimate_read_length(f, 100)
            read_info[base][fname] = length

    return read_info

def estimate_read_length(fastq_path: str, n_reads: int = 100) -> int:
    """Estimate average read length from FASTQ"""
    opener = gzip.open if fastq_path.endswith('.gz') else open
    lengths = []

    try:
        with opener(fastq_path, 'rt') as f:
            count = 0
            while count < n_reads:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip()
                plus = f.readline()
                qual = f.readline()
                if seq:
                    lengths.append(len(seq))
                    count += 1
    except Exception:
        return 0

    return int(sum(lengths) / len(lengths)) if lengths else 0

def detect_chemistry(fastq_dir: str, sample_name: str) -> dict:
    """
    Detect chemistry/platform from FASTQ structure
    Returns: {
        'chemistry': '10x3v3'|'10x3v2'|'10x5v2'|'seqwell_s3'|'gexscope'|'unknown',
        'platform': '10x'|'seqwell'|'gexscope'|'unknown',
        'quantifier': 'cellranger'|'starsolo',
        'read_structure': {...},
        'confidence': 'high'|'medium'|'low'
    }
    """
    result = {
        'sample': sample_name,
        'chemistry': 'unknown',
        'platform': 'unknown',
        'quantifier': 'starsolo',  # Default fallback
        'confidence': 'low',
        'details': {}
    }

    # Get read structure
    read_info = detect_read_structure(fastq_dir)
    if 'error' in read_info:
        result['error'] = read_info['error']
        return result

    result['read_structure'] = read_info

    # Analyze file patterns and read lengths
    all_files = []
    for base, files in read_info.items():
        all_files.extend(files.keys())

    # Check for 10x naming pattern
    has_10x_pattern = any(re.search(r'_S\d+_L\d+_[IR][12]_001', f) for f in all_files)

    # Analyze read lengths for chemistry detection
    r1_lengths = []
    r2_lengths = []
    i1_lengths = []
    i2_lengths = []

    for base, files in read_info.items():
        for fname, length in files.items():
            if '_R1' in fname or 'R1_' in fname or fname.endswith('_1.fastq.gz'):
                r1_lengths.append(length)
            elif '_R2' in fname or 'R2_' in fname or fname.endswith('_2.fastq.gz'):
                r2_lengths.append(length)
            elif '_I1' in fname or 'I1_' in fname:
                i1_lengths.append(length)
            elif '_I2' in fname or 'I2_' in fname:
                i2_lengths.append(length)

    # Get representative lengths
    r1_len = max(r1_lengths) if r1_lengths else 0
    r2_len = max(r2_lengths) if r2_lengths else 0
    i1_len = max(i1_lengths) if i1_lengths else 0
    i2_len = max(i2_lengths) if i2_lengths else 0

    # Detect chemistry based on read structure
    if has_10x_pattern:
        # Standard 10x format detected
        if r1_len >= 26 and r1_len <= 28:
            if i1_len == 8 or i1_len == 10:
                if i2_len == 0:
                    # 10x 3' v2 (28+8+91)
                    result['chemistry'] = '10x3v2'
                    result['platform'] = '10x'
                    result['quantifier'] = 'cellranger'
                    result['confidence'] = 'high'
                elif i2_len >= 8:
                    # 10x 3' v3 (28+8+8+91)
                    result['chemistry'] = '10x3v3'
                    result['platform'] = '10x'
                    result['quantifier'] = 'cellranger'
                    result['confidence'] = 'high'
        elif r1_len >= 150:
            # 10x 5' v2 (150bp R1, 150bp R2)
            result['chemistry'] = '10x5v2'
            result['platform'] = '10x'
            result['quantifier'] = 'cellranger'
            result['confidence'] = 'high'
    else:
        # Non-10x format - check for known platforms

        # Seqwell S^3 detection
        # Typically: R1=20bp (cell barcode), R2=50-150bp (transcript)
        if r1_len >= 16 and r1_len <= 24 and r2_len >= 50:
            result['chemistry'] = 'seqwell_s3'
            result['platform'] = 'seqwell'
            result['quantifier'] = 'starsolo'
            result['confidence'] = 'medium'
            result['starsolo_params'] = {
                'soloType': 'CB_UMI_Simple',
                'soloCBstart': 1,
                'soloCBlen': 12,
                'soloUMIstart': 13,
                'soloUMIlen': 8,
                'soloBarcodeMate': 1,
                'soloCBwhitelist': 'seqwell_whitelist.txt'
            }

        # GeXscope detection
        # Similar to 10x but different barcode structure
        elif r1_len >= 26 and r1_len <= 30 and r2_len >= 90:
            result['chemistry'] = 'gexscope'
            result['platform'] = 'gexscope'
            result['quantifier'] = 'starsolo'
            result['confidence'] = 'medium'
            result['starsolo_params'] = {
                'soloType': 'CB_UMI_Simple',
                'soloCBstart': 1,
                'soloCBlen': 16,
                'soloUMIstart': 17,
                'soloUMIlen': 12,
                'soloBarcodeMate': 1
            }

        # Generic non-10x single-cell
        elif r1_len > 0 and r2_len > 0:
            result['chemistry'] = 'custom'
            result['platform'] = 'unknown'
            result['quantifier'] = 'starsolo'
            result['confidence'] = 'low'
            # Use generic STARsolo parameters
            result['starsolo_params'] = {
                'soloType': 'CB_UMI_Simple',
                'soloCBstart': 1,
                'soloCBlen': 16,
                'soloUMIstart': 17,
                'soloUMIlen': 10,
                'soloBarcodeMate': 1
            }

    # Add detailed read structure
    result['details'] = {
        'r1_length': r1_len,
        'r2_length': r2_len,
        'i1_length': i1_len,
        'i2_length': i2_len,
        'has_10x_naming': has_10x_pattern,
        'file_count': len(all_files)
    }

    print(f"Chemistry detection for {sample_name}:", file=sys.stderr)
    print(f"  Platform: {result['platform']}", file=sys.stderr)
    print(f"  Chemistry: {result['chemistry']}", file=sys.stderr)
    print(f"  Quantifier: {result['quantifier']}", file=sys.stderr)
    print(f"  Confidence: {result['confidence']}", file=sys.stderr)
    print(f"  Read structure: R1={r1_len}bp, R2={r2_len}bp, I1={i1_len}bp, I2={i2_len}bp", file=sys.stderr)

    return result

def main():
    if len(sys.argv) != 4:
        print("Usage: detect_chemistry.py <fastq_dir> <output.json> <sample_name>", file=sys.stderr)
        sys.exit(1)

    fastq_dir = sys.argv[1]
    output_file = sys.argv[2]
    sample_name = sys.argv[3]

    # Detect chemistry
    result = detect_chemistry(fastq_dir, sample_name)

    # Write JSON output
    with open(output_file, 'w') as f:
        json.dump(result, f, indent=2)

    # Also write TSV for easy parsing
    tsv_file = output_file.replace('.json', '.tsv')
    with open(tsv_file, 'w') as f:
        f.write("Sample\tPlatform\tChemistry\tQuantifier\tConfidence\n")
        f.write(f"{sample_name}\t{result['platform']}\t{result['chemistry']}\t{result['quantifier']}\t{result['confidence']}\n")

    # Always exit 0 - routing is handled via JSON output
    # The quantifier field in JSON determines which tool to use
    sys.exit(0)

if __name__ == "__main__":
    main()