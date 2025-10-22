#!/usr/bin/env python3
import csv
import sys
import re
import os
import requests
from urllib.parse import urljoin
import time

"""
Extract BAM file URLs from SRA metadata for given run accessions.
This script attempts to find BAM files in the SRA S3 bucket using common naming patterns.
"""

def get_sra_metadata(run_id):
    """Fetch SRA metadata for a run ID"""
    try:
        url = f"https://www.ncbi.nlm.nih.gov/sra/{run_id}"
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            return response.text
    except Exception as e:
        print(f"Warning: Could not fetch metadata for {run_id}: {e}", file=sys.stderr)
    return None
def get_geo_metadata(gsm_id: str):
    """Fetch GEO metadata page for a GSM accession"""
    if not gsm_id:
        return None
    try:
        # Use the GEO accession page (HTML)
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm_id}"
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            return response.text
    except Exception as e:
        print(f"Warning: Could not fetch GEO page for {gsm_id}: {e}", file=sys.stderr)
    return None

def extract_candidate_names_from_geo(geo_html: str):
    """Extract candidate BAM/FASTQ-related names from GEO page (titles, sample names)"""
    if not geo_html:
        return []
    candidates = set()
    # Try to capture potential file basenames from GEO supplementary file links and titles
    # Example: D17PrPzF_BE in 'D17PrPzF_BE.bam.1'
    title_match = re.findall(r"<title>([^<]+)</title>", geo_html, re.IGNORECASE)
    for t in title_match:
        t = t.strip()
        t = re.sub(r"[^\w\-]+", "_", t)
        if t:
            candidates.add(t)
    # Supplementary file names
    supp_links = re.findall(r"href=\"([^\"]+?(?:\.bam(?:\.1)?|\.fastq(?:\.gz)?))\"", geo_html, re.IGNORECASE)
    for link in supp_links:
        base = os.path.basename(link)
        stem = re.sub(r"\.(bam|fastq(?:\.gz)?)$", "", base, flags=re.IGNORECASE)
        if stem:
            candidates.add(stem)
    # Also parse nearby text for likely basenames
    name_spans = re.findall(r">([A-Za-z0-9][A-Za-z0-9_\-]{3,})<", geo_html)
    for n in name_spans:
        candidates.add(n)
    return sorted(candidates)


def extract_sample_name_from_metadata(metadata_text, run_id):
    """Extract sample name from SRA metadata page"""
    if not metadata_text:
        return None
    
    # Look for sample name patterns in the metadata
    patterns = [
        r'<td[^>]*>Sample name</td>\s*<td[^>]*>([^<]+)</td>',
        r'<td[^>]*>Sample</td>\s*<td[^>]*>([^<]+)</td>',
        r'<span[^>]*class="[^"]*sample[^"]*"[^>]*>([^<]+)</span>',
        r'<a[^>]*href="[^"]*sample[^"]*"[^>]*>([^<]+)</a>'
    ]
    
    for pattern in patterns:
        match = re.search(pattern, metadata_text, re.IGNORECASE)
        if match:
            sample_name = match.group(1).strip()
            # Clean up the sample name
            sample_name = re.sub(r'[^\w\-_]', '_', sample_name)
            return sample_name
    
    return None

def extract_s3_bam_urls_from_metadata(metadata_text: str, run_id: str):
    """Parse SRA run page HTML for direct S3 BAM links (sra-pub-src-1)"""
    if not metadata_text:
        return []
    urls = set()
    # Common patterns for S3 public source bucket links that include BAMs
    # Examples:
    # https://sra-pub-src-1.s3.amazonaws.com/SRR7545330/D17PrPzF_BE.bam.1
    # https://sra-pub-src-1.s3.amazonaws.com/SRRXXXXXXX/<anything>.bam
    s3_pattern = re.compile(r"https?://sra-pub-src-1\.s3\.amazonaws\.com/" + re.escape(run_id) + r"/[^\s\"'<>]+?\.bam(?:\.1)?", re.IGNORECASE)
    for m in s3_pattern.finditer(metadata_text):
        urls.add(m.group(0))

    # Also catch protocol-relative //sra-pub-src-1.s3.amazonaws.com/... links
    s3_proto_rel = re.compile(r"//sra-pub-src-1\.s3\.amazonaws\.com/" + re.escape(run_id) + r"/[^\s\"'<>]+?\.bam(?:\.1)?", re.IGNORECASE)
    for m in s3_proto_rel.finditer(metadata_text):
        urls.add("https:" + m.group(0))

    return sorted(urls)

def generate_bam_urls(run_id, sample_name=None, extra_candidates=None):
    """Generate potential BAM URLs for a run ID"""
    base_url = "https://sra-pub-src-1.s3.amazonaws.com"
    
    # Common BAM file naming patterns
    potential_names = []
    
    if sample_name:
        # Try with sample name
        potential_names.extend([
            f"{sample_name}.bam.1",
            f"{sample_name}.bam",
            f"{run_id}_{sample_name}.bam.1",
            f"{run_id}_{sample_name}.bam"
        ])
    # Add any extra candidate basenames from GEO
    if extra_candidates:
        for cand in extra_candidates:
            potential_names.extend([
                f"{cand}.bam.1",
                f"{cand}.bam",
                f"{run_id}_{cand}.bam.1",
                f"{run_id}_{cand}.bam"
            ])
    
    # Try with run ID variations
    potential_names.extend([
        f"{run_id}.bam.1",
        f"{run_id}.bam",
        f"{run_id}_1.bam.1",
        f"{run_id}_1.bam"
    ])
    
    # Try common 10x Genomics patterns
    potential_names.extend([
        f"{run_id}_possorted_genome_bam.bam.1",
        f"{run_id}_possorted_genome_bam.bam",
        f"{run_id}_filtered_feature_bc_matrix.bam.1",
        f"{run_id}_filtered_feature_bc_matrix.bam"
    ])
    
    urls = []
    for name in potential_names:
        urls.append(f"{base_url}/{run_id}/{name}")
    
    return urls

def check_url_exists(url):
    """Check if a URL exists (HEAD request)"""
    try:
        response = requests.head(url, timeout=10, allow_redirects=True)
        return response.status_code == 200
    except:
        return False

def find_bam_url(run_id, gsm_id=None):
    """Find the correct BAM URL for a run ID"""
    print(f"Searching for BAM file for {run_id}...", file=sys.stderr)
    
    # First try to get metadata to extract sample name
    metadata = get_sra_metadata(run_id)
    sample_name = None
    if metadata:
        sample_name = extract_sample_name_from_metadata(metadata, run_id)
        if sample_name:
            print(f"Found sample name: {sample_name}", file=sys.stderr)
        # Try to parse direct S3 BAM URLs from the metadata page
        direct_urls = extract_s3_bam_urls_from_metadata(metadata, run_id)
        for url in direct_urls:
            if check_url_exists(url):
                print(f"Found BAM file via page link: {url}", file=sys.stderr)
                return url

    # Try GEO for additional candidate names
    geo_candidates = []
    if gsm_id:
        geo_html = get_geo_metadata(gsm_id)
        if geo_html:
            geo_candidates = extract_candidate_names_from_geo(geo_html)
            if geo_candidates:
                print(f"GEO candidates for {gsm_id}: {', '.join(geo_candidates[:5])}...", file=sys.stderr)
    
    # Generate potential URLs
    potential_urls = generate_bam_urls(run_id, sample_name, geo_candidates)
    
    # Check each URL
    for url in potential_urls:
        if check_url_exists(url):
            print(f"Found BAM file: {url}", file=sys.stderr)
            return url
        time.sleep(0.1)  # Be nice to the server
    
    print(f"No BAM file found for {run_id}", file=sys.stderr)
    return "N/A"

def main():
    if len(sys.argv) < 3:
        print("Usage: extract_bam_urls.py <runs_tsv> <output_tsv> [gsm_id] [hints_tsv]", file=sys.stderr)
        sys.exit(2)
    
    runs_file = sys.argv[1]
    output_file = sys.argv[2]
    gsm_id = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] not in ("", "-", "None", "null") else None
    hints_file = sys.argv[4] if len(sys.argv) > 4 and sys.argv[4] not in ("", "-", "None", "null") else None
    
    run_ids = []
    with open(runs_file, 'r') as f:
        for line in f:
            run_id = line.strip()
            if run_id:
                run_ids.append(run_id)
    
    # Optional hints mapping SRR -> URL
    hints = {}
    if hints_file and os.path.exists(hints_file):
        try:
            with open(hints_file, 'r') as hf:
                for line in hf:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = re.split(r"\t|,|\s+", line, maxsplit=1)
                    if len(parts) == 2:
                        hints[parts[0]] = parts[1]
        except Exception as e:
            print(f"Warning: Failed to read hints file {hints_file}: {e}", file=sys.stderr)

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['RunAccession', 'BamUrl'])
        
        for run_id in run_ids:
            if run_id in hints:
                bam_url = hints[run_id]
                print(f"Using hinted BAM URL for {run_id}: {bam_url}", file=sys.stderr)
            else:
                bam_url = find_bam_url(run_id, gsm_id)
            writer.writerow([run_id, bam_url])
            f.flush()  # Write immediately to see progress

if __name__ == '__main__':
    main()
