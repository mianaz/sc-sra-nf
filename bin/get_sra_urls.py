#!/usr/bin/env python3
"""
Get direct SRA file URLs from NCBI SDL API
"""
import sys
import requests
import json

def get_sra_url(sra_id: str) -> str:
    """Get direct S3 URL for SRA file"""
    try:
        sdl_url = f"https://www.ncbi.nlm.nih.gov/Traces/sdl/2/retrieve?acc={sra_id}"
        response = requests.get(sdl_url, timeout=15)
        if response.status_code == 200:
            data = response.json()

            # Find the SRA file
            for result_item in data.get('result', []):
                for file_info in result_item.get('files', []):
                    if file_info.get('type', '').lower() == 'sra':
                        locations = file_info.get('locations', [])
                        if locations:
                            return locations[0].get('link', '')
    except Exception as e:
        print(f"Error fetching URL for {sra_id}: {e}", file=sys.stderr)

    return ''

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: get_sra_urls.py <runs.tsv> <output.tsv>", file=sys.stderr)
        sys.exit(1)

    runs_file = sys.argv[1]
    output_file = sys.argv[2]

    with open(runs_file, 'r') as f_in, open(output_file, 'w') as f_out:
        f_out.write("RunID\tURL\n")
        for line in f_in:
            sra_id = line.strip()
            if not sra_id:
                continue

            print(f"Getting URL for {sra_id}...", file=sys.stderr)
            url = get_sra_url(sra_id)
            if url:
                print(f"  Found: {url}", file=sys.stderr)
                f_out.write(f"{sra_id}\t{url}\n")
            else:
                print(f"  No URL found, will use prefetch", file=sys.stderr)
                f_out.write(f"{sra_id}\tN/A\n")
