#!/usr/bin/env python3
"""
Detect available data formats (FASTQ, BAM, or SRA) for SRA accessions
and extract download URLs when available.
"""
import sys
import re
import requests
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional
import json

def fetch_sra_metadata(sra_id: str) -> Optional[ET.Element]:
    """Fetch metadata from ENA for an SRA accession"""
    base_url = "https://www.ebi.ac.uk/ena/browser/api/xml"
    try:
        response = requests.get(f"{base_url}/{sra_id}", timeout=10)
        if response.status_code == 200:
            return ET.fromstring(response.text)
        else:
            # Try NCBI as fallback
            ncbi_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={sra_id}&retmode=xml"
            response = requests.get(ncbi_url, timeout=10)
            if response.status_code == 200:
                return ET.fromstring(response.text)
    except Exception as e:
        print(f"Error fetching metadata for {sra_id}: {e}", file=sys.stderr)
    return None

def extract_s3_bam_urls_from_page(html_text: str, run_id: str) -> List[str]:
    """Extract S3 BAM URLs directly from SRA page HTML (sra-pub-src-* buckets)"""
    if not html_text:
        return []
    urls = set()
    # Pattern for S3 public source bucket links that include BAMs
    # Example: https://sra-pub-src-1.s3.amazonaws.com/SRR7545330/D17PrPzF_BE.bam.1
    s3_pattern = re.compile(r"https?://sra-pub-src-\d+\.s3\.amazonaws\.com/" + re.escape(run_id) + r"/[^\s\"'<>]+?\.bam(?:\.1)?", re.IGNORECASE)
    for m in s3_pattern.finditer(html_text):
        urls.add(m.group(0))

    # Also catch protocol-relative //sra-pub-src-*.s3.amazonaws.com/... links
    s3_proto_rel = re.compile(r"//sra-pub-src-\d+\.s3\.amazonaws\.com/" + re.escape(run_id) + r"/[^\s\"'<>]+?\.bam(?:\.1)?", re.IGNORECASE)
    for m in s3_proto_rel.finditer(html_text):
        urls.add("https:" + m.group(0))

    return sorted(urls)

def detect_available_formats(sra_id: str) -> Dict[str, List[str]]:
    """
    Detect available formats and return URLs
    Returns: {
        'format': 'FASTQ'|'BAM'|'SRA'|'UNKNOWN',
        'urls': [...],
        'method': 'direct'|'sra-tools'
    }
    """
    result = {
        'format': 'UNKNOWN',
        'urls': [],
        'method': 'sra-tools',
        'paired': False
    }

    # Check SRA SDL (Sequence Data Locator) API for original format files with free egress
    # This API shows file types and S3 locations, allowing us to detect:
    # 1. If original format files (BAM, FASTQ) are available
    # 2. If they're in sra-pub-src-* buckets (free egress)
    try:
        sdl_url = f"https://www.ncbi.nlm.nih.gov/Traces/sdl/2/retrieve?acc={sra_id}"
        response = requests.get(sdl_url, timeout=15)
        if response.status_code == 200:
            data = response.json()

            # Check if any files have original format with free egress
            for result_item in data.get('result', []):
                for file_info in result_item.get('files', []):
                    file_type = file_info.get('type', '')
                    locations = file_info.get('locations', [])

                    # Check for original format types (TenX, bam, fastq, etc.)
                    # and verify they're in free egress buckets (sra-pub-src-*)
                    is_original_format = file_type.lower() in ['tenx', 'bam', 'fastq']
                    has_free_egress = any(
                        'sra-pub-src-' in loc.get('link', '')
                        for loc in locations
                    )

                    if is_original_format and has_free_egress:
                        print(f"  Detected original format file: type={file_type}", file=sys.stderr)
                        print(f"  Free egress available in sra-pub-src-* bucket", file=sys.stderr)

                        # Determine format based on type
                        if file_type.lower() in ['tenx', 'bam']:
                            result['format'] = 'BAM'
                        elif file_type.lower() == 'fastq':
                            result['format'] = 'FASTQ'

                        result['urls'] = []  # Will be found by extract_bam_urls.py
                        result['method'] = 'direct'
                        return result

            print(f"  No original format files with free egress found", file=sys.stderr)
    except Exception as e:
        print(f"  Warning: Could not check SDL API for {sra_id}: {e}", file=sys.stderr)

    # Check ENA for direct downloads
    root = fetch_sra_metadata(sra_id)
    if root is None:
        # Default to SRA if can't fetch metadata
        result['format'] = 'SRA'
        return result

    # Look for file URLs in ENA metadata
    fastq_urls = []
    bam_urls = []

    # Parse ENA XML structure
    for run in root.findall('.//RUN'):
        run_acc = run.get('accession', '')
        if run_acc != sra_id:
            continue

        # Check for file links
        for file_elem in run.findall('.//FILE'):
            filename = file_elem.get('filename', '')
            file_url = file_elem.get('url', '')

            if not file_url:
                # Try alternative URL locations
                for url_elem in file_elem.findall('.//URL'):
                    file_url = url_elem.text
                    if file_url:
                        break

            if file_url:
                if filename.endswith(('.fastq.gz', '.fq.gz')):
                    fastq_urls.append(file_url)
                elif filename.endswith('.bam'):
                    bam_urls.append(file_url)

    # Also check for direct FTP links
    for links in root.findall('.//XREF_LINK'):
        for url in links.findall('.//URL'):
            url_text = url.text
            if url_text:
                if 'fastq' in url_text.lower():
                    fastq_urls.append(url_text)
                elif '.bam' in url_text:
                    bam_urls.append(url_text)

    # Check AWS Open Data (processed files)
    aws_base = f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{sra_id}/{sra_id}"

    # Try to detect if original format files exist on AWS
    try:
        # Check for BAM
        response = requests.head(f"{aws_base}.bam", timeout=5)
        if response.status_code == 200:
            bam_urls.append(f"{aws_base}.bam")
    except:
        pass

    try:
        # Check for FASTQ (common patterns)
        for suffix in ['_1.fastq.gz', '_2.fastq.gz', '.fastq.gz']:
            response = requests.head(f"{aws_base}{suffix}", timeout=5)
            if response.status_code == 200:
                fastq_urls.append(f"{aws_base}{suffix}")
    except:
        pass

    # Determine best format
    if fastq_urls:
        result['format'] = 'FASTQ'
        result['urls'] = fastq_urls
        result['method'] = 'direct'
        result['paired'] = len(fastq_urls) >= 2
    elif bam_urls:
        result['format'] = 'BAM'
        result['urls'] = bam_urls
        result['method'] = 'direct'
    else:
        # Default to SRA download via sra-tools
        result['format'] = 'SRA'
        result['method'] = 'sra-tools'

    return result

def process_run_list(runs_file: str, output_file: str, sample_id: str):
    """Process a list of SRA runs and output format detection results"""

    results = {}
    with open(runs_file, 'r') as f:
        for line in f:
            sra_id = line.strip()
            if not sra_id:
                continue

            print(f"Detecting format for {sra_id}...", file=sys.stderr)
            format_info = detect_available_formats(sra_id)
            results[sra_id] = format_info

            print(f"  {sra_id}: {format_info['format']} via {format_info['method']}", file=sys.stderr)

    # Determine overall strategy for the sample
    formats = [r['format'] for r in results.values()]
    if all(f == 'FASTQ' for f in formats):
        overall_strategy = 'FASTQ_DIRECT'
    elif all(f == 'BAM' for f in formats):
        overall_strategy = 'BAM_DIRECT'
    elif 'FASTQ' in formats or 'BAM' in formats:
        overall_strategy = 'MIXED'
    else:
        overall_strategy = 'SRA_TOOLS'

    # Output results
    output = {
        'sample': sample_id,
        'strategy': overall_strategy,
        'runs': results
    }

    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2)

    # Also create TSV for easier parsing
    tsv_file = output_file.replace('.json', '.tsv')
    with open(tsv_file, 'w') as f:
        f.write("RunID\tFormat\tMethod\tURLs\n")
        for run_id, info in results.items():
            urls = ','.join(info['urls']) if info['urls'] else 'N/A'
            f.write(f"{run_id}\t{info['format']}\t{info['method']}\t{urls}\n")

    print(f"Sample {sample_id}: Strategy={overall_strategy}", file=sys.stderr)
    return overall_strategy

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: detect_data_format.py <runs.tsv> <output.json> <sample_id>", file=sys.stderr)
        sys.exit(1)

    process_run_list(sys.argv[1], sys.argv[2], sys.argv[3])