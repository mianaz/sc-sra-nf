#!/usr/bin/env python3
import csv
import sys
from collections import defaultdict

"""
Reads SraRunTable.csv robustly (handles embedded commas/quotes) and emits two TSVs:
- sample_run_map.tsv: columns SampleAccession (GSM), RunAccession (SRR)
- run_metadata.tsv: the cleaned full metadata with tab delimiter

Detects column names case-insensitively for 'Sample' (GSM*) and 'Run' (SRR*),
falling back to common SRA headers: 'Sample Name' or 'Sample Accession', 'Run' or 'Run accession'.
"""

def normalize_header(name: str) -> str:
    return name.strip().lower().replace(' ', '_')


def find_column(headers, candidates):
    norm = {normalize_header(h): h for h in headers}
    for cand in candidates:
        key = normalize_header(cand)
        if key in norm:
            return norm[key]
    # try partial contains
    for h in headers:
        nh = normalize_header(h)
        for cand in candidates:
            if normalize_header(cand).replace('_', '') in nh.replace('_', ''):
                return h
    return None


def main():
    if len(sys.argv) < 2:
        print("Usage: parse_sra_table.py SraRunTable.csv", file=sys.stderr)
        sys.exit(2)
    infile = sys.argv[1]

    with open(infile, newline='') as fh:
        # Dialect sniffing for robustness
        sample = fh.read(4096)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample)
        except Exception:
            dialect = csv.excel
        reader = csv.DictReader(fh, dialect=dialect)
        headers = reader.fieldnames or []
        if not headers:
            print("ERROR: No headers detected in SraRunTable.csv", file=sys.stderr)
            sys.exit(1)

        sample_col = find_column(headers, [
            'Sample Accession', 'Sample', 'Sample Name', 'GSM', 'sample_accession'
        ])
        run_col = find_column(headers, [
            'Run', 'Run accession', 'Run Accession', 'SRR', 'run'
        ])
        if not sample_col or not run_col:
            print(f"ERROR: Could not find sample/run columns. headers={headers}", file=sys.stderr)
            sys.exit(1)

        # Write full cleaned metadata as TSV
        with open('run_metadata.tsv', 'w', newline='') as meta_out:
            tsv_writer = csv.DictWriter(meta_out, fieldnames=headers, delimiter='\t')
            tsv_writer.writeheader()

            # Build mapping
            pairs = []
            for row in reader:
                tsv_writer.writerow(row)
                s = row.get(sample_col, '').strip()
                r = row.get(run_col, '').strip()
                if s and r:
                    pairs.append((s, r))

        # Deduplicate while preserving order
        seen = set()
        uniq_pairs = []
        for s, r in pairs:
            key = (s, r)
            if key not in seen:
                seen.add(key)
                uniq_pairs.append(key)

        with open('sample_run_map.tsv', 'w', newline='') as map_out:
            w = csv.writer(map_out, delimiter='\t')
            w.writerow(['SampleAccession', 'RunAccession'])
            w.writerows(uniq_pairs)


if __name__ == '__main__':
    main()
