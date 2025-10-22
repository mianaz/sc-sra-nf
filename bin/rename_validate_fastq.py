#!/usr/bin/env python3
import gzip
import os
import re
import sys
from glob import glob

"""
Fixed version of rename_validate_fastq.py that ensures unique lane assignments per run.
This prevents lane collisions when multiple runs have the same lane numbers.

Key changes:
- Always uses position-based lane assignment from run group
- Removes header/filename lane detection to prevent collisions
- Maintains read type detection (I1/I2/R1/R2) via read lengths
"""

def detect_lane_from_run_position(run_id: str, run_list: list) -> str:
    """
    Assign unique lane number based on position in run group.
    This ensures no collisions between runs.
    """
    try:
        if run_id in run_list:
            lane_num = run_list.index(run_id) + 1
            return f"L{lane_num:03d}"
    except (ValueError, IndexError):
        pass
    return "L001"


def estimate_read_length(fastq_path: str, n_reads: int = 1000) -> int:
    """Estimate average read length from a FASTQ file"""
    opener = gzip.open if fastq_path.endswith('.gz') else open
    length_sum = 0
    count = 0
    try:
        with opener(fastq_path, 'rt') as fh:
            while count < n_reads:
                hdr = fh.readline()
                if not hdr:
                    break
                seq = fh.readline().strip()
                plus = fh.readline()
                qual = fh.readline()
                if not qual:
                    break
                length_sum += len(seq)
                count += 1
    except Exception:
        return 0
    return int(length_sum / count) if count else 0


def analyze_fastq_files(fastqs: list, run_id: str) -> dict:
    """Analyze FASTQ files and map to read types based on read lengths"""
    count = len(fastqs)
    print(f"Analyzing {count} FASTQ files for run {run_id}", file=sys.stderr)

    # Get read lengths for intelligent mapping
    read_info = []
    for f in fastqs:
        length = estimate_read_length(f, 200)
        read_info.append((f, length))
        print(f"  {os.path.basename(f)}: ~{length}bp", file=sys.stderr)

    # Sort by read length to help with mapping
    read_info.sort(key=lambda x: x[1])

    mapping = {}

    if count == 1:
        # Single-end sequencing
        mapping['R1'] = read_info[0][0]
        print("  Detected: Single-end sequencing", file=sys.stderr)

    elif count == 2:
        # Could be paired-end or index + read
        len1, len2 = read_info[0][1], read_info[1][1]

        if len1 <= 35 and len2 >= 50:
            # Short read is likely index, long read is R1
            mapping['I1'] = read_info[0][0]
            mapping['R1'] = read_info[1][0]
            print("  Detected: Index + R1 sequencing", file=sys.stderr)
        else:
            # Assume paired-end
            if len1 <= len2:
                mapping['R1'] = read_info[0][0]
                mapping['R2'] = read_info[1][0]
            else:
                mapping['R1'] = read_info[1][0]
                mapping['R2'] = read_info[0][0]
            print("  Detected: Paired-end sequencing", file=sys.stderr)

    elif count == 3:
        # Typical 10x: I1, R1, R2
        mapping['I1'] = read_info[0][0]  # shortest
        mapping['R1'] = read_info[1][0]  # middle
        mapping['R2'] = read_info[2][0]  # longest
        print("  Detected: Index + paired-end (typical 10x)", file=sys.stderr)

    elif count == 4:
        # Full 4-file: I1, I2, R1, R2
        mapping['I1'] = read_info[0][0]  # shortest
        mapping['I2'] = read_info[1][0]  # second shortest
        mapping['R1'] = read_info[2][0]  # third
        mapping['R2'] = read_info[3][0]  # longest
        print("  Detected: Dual index + paired-end", file=sys.stderr)

    else:
        print(f"WARNING: Unusual number of FASTQ files ({count}). Using best-effort mapping.", file=sys.stderr)
        # Best effort for unusual cases
        if count >= 2:
            mapping['R1'] = read_info[-2][0]  # second longest
            mapping['R2'] = read_info[-1][0]  # longest
        if count >= 3:
            mapping['I1'] = read_info[0][0]   # shortest
        if count >= 4:
            mapping['I2'] = read_info[1][0]   # second shortest

    # Validation warnings
    if 'R1' in mapping and 'R2' in mapping:
        r1_len = estimate_read_length(mapping['R1'], 100)
        r2_len = estimate_read_length(mapping['R2'], 100)
        if r1_len > r2_len + 10:  # Allow some tolerance
            print(f"WARNING: R1 ({r1_len}bp) longer than R2 ({r2_len}bp). Unusual for 10x.", file=sys.stderr)

    return mapping


def main():
    if len(sys.argv) < 4:
        print("Usage: rename_validate_fastq.py <run_id> <sample_id> <indir> <outdir> [run_group]", file=sys.stderr)
        sys.exit(2)

    run_id = sys.argv[1]
    sample_id = sys.argv[2]
    indir = sys.argv[3]
    outdir = sys.argv[4] if len(sys.argv) > 4 else indir
    run_group = sys.argv[5].split(',') if len(sys.argv) > 5 else None

    os.makedirs(outdir, exist_ok=True)

    # Find FASTQ files (search recursively for bamtofastq nested structure)
    fastqs = sorted(glob(os.path.join(indir, '*.fastq'))) + sorted(glob(os.path.join(indir, '*.fastq.gz')))

    # Check if this is bamtofastq output (nested structure with subdirectories)
    is_bamtofastq = False
    if not fastqs:
        # Try recursive search for nested bamtofastq outputs (e.g., SRR_ID/flowcell_name/*.fastq.gz)
        fastqs = sorted(glob(os.path.join(indir, '**', '*.fastq'), recursive=True)) + \
                 sorted(glob(os.path.join(indir, '**', '*.fastq.gz'), recursive=True))
        if fastqs:
            is_bamtofastq = True

    if not fastqs:
        print(f"ERROR: No FASTQ files found in {indir}", file=sys.stderr)
        sys.exit(1)

    print(f"Processing run {run_id} -> sample {sample_id}", file=sys.stderr)

    # Handle bamtofastq output: copy all files preserving structure for CellRanger
    if is_bamtofastq and len(fastqs) > 10:  # bamtofastq typically produces many files
        print(f"Detected bamtofastq output with {len(fastqs)} files across multiple lanes/flowcells", file=sys.stderr)
        print(f"Copying all files for CellRanger compatibility", file=sys.stderr)

        import shutil
        written = []
        for src in fastqs:
            # Extract just the filename (e.g., bamtofastq_S1_L001_R1_001.fastq.gz)
            basename = os.path.basename(src)
            # Rename to sample format: GSM_S1_L001_R1_001.fastq.gz
            parts = basename.split('_')
            if len(parts) >= 5 and parts[0] == 'bamtofastq':
                # Replace 'bamtofastq' with sample name
                new_basename = '_'.join([sample_id] + parts[1:])
                dest = os.path.join(outdir, new_basename)

                try:
                    shutil.copy2(src, dest)
                    written.append(dest)
                except Exception as e:
                    print(f"ERROR: Failed to copy {src} to {dest}: {e}", file=sys.stderr)
                    sys.exit(1)

        print(f"Successfully copied {len(written)} bamtofastq files", file=sys.stderr)

        # Print resulting files for Nextflow
        for f in written:
            print(f)

        sys.exit(0)

    # FIXED: Always use position-based lane assignment to prevent collisions
    if run_group:
        lane = detect_lane_from_run_position(run_id, run_group)
        print(f"Assigned lane {lane} based on run position in group (run {run_group.index(run_id)+1} of {len(run_group)})", file=sys.stderr)
    else:
        lane = "L001"
        print(f"No run group provided, using default lane {lane}", file=sys.stderr)

    # Analyze files and map to read types
    mapping = analyze_fastq_files(fastqs, run_id)

    # Write outputs with validation
    written = []
    required_files = ['R1']  # At minimum need R1

    for role in ['I1', 'I2', 'R1', 'R2']:
        src = mapping.get(role)
        if not src:
            continue

        dest = os.path.join(outdir, f"{sample_id}_S1_{lane}_{role}_001.fastq.gz")

        try:
            # Copy and compress if needed
            if src.endswith('.gz'):
                with open(src, 'rb') as s, open(dest, 'wb') as d:
                    d.write(s.read())
            else:
                with open(src, 'rb') as s, gzip.open(dest, 'wb') as d:
                    d.write(s.read())

            written.append(dest)
            print(f"Created: {os.path.basename(dest)}", file=sys.stderr)

        except Exception as e:
            print(f"ERROR: Failed to process {role} file {src}: {e}", file=sys.stderr)
            sys.exit(1)

    # Validate we have minimum required files
    created_roles = []
    for f in written:
        basename = os.path.basename(f)
        # Parse filename: GSM_TestSample_S1_L001_R1_001.fastq.gz
        parts = basename.split('_')
        if len(parts) >= 5:
            role = parts[3]  # R1, R2, I1, I2 (at index 3)
            created_roles.append(role)

    if 'R1' not in created_roles:
        print("ERROR: Failed to create required R1 file", file=sys.stderr)
        sys.exit(1)

    print(f"Successfully processed {len(written)} files for {run_id} with lane {lane}", file=sys.stderr)

    # Print resulting files one per line for Nextflow collection
    for f in written:
        print(f)


if __name__ == '__main__':
    main()