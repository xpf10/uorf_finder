# uORF Finder

A Python tool for extracting upstream open reading frames (uORFs) from gene transcripts. This tool searches for uORFs in the 5' UTR region of transcripts and outputs the sequences in FASTA format, with optional BED file output for genomic coordinates.

## Overview

uORF Finder is a bioinformatics tool that:
- Fetches transcript sequences for a given gene from NCBI
- Identifies upstream open reading frames (uORFs) in the 5' UTR region
- Supports both canonical (ATG) and non-canonical start codons
- Outputs results in FASTA format with detailed headers
- Optionally generates BED files with genomic coordinates

## Installation

This project uses `uv` for dependency management. To install:

1. Ensure you have Python 3.13+ installed
2. Install `uv` if not already installed:
   ```bash
   pip install uv
   ```
3. Clone the repository and install dependencies:
   ```bash
   git clone https://github.com/xpf10/uorf_finder
   cd uorf_finder
   uv sync
   ```

## Usage

Basic usage:
```bash
python main.py <gene_name> -o output.fasta
```

### Command Line Options

- `gene`: Gene name (e.g., TP53, BRCA1)
- `-o, --output`: Output FASTA file name (default: uorfs.fasta)
- `-b, --bed`: Output BED file name (default: same as FASTA with .bed extension)
- `-s, --species`: Species (default: human, also supports mouse)
- `-c, --canonical-only`: Only search for ATG start codons (default: includes non-canonical)
- `-m, --min-length`: Minimum uORF length in bp (default: 9)
- `-e, --email`: NCBI Entrez email (required by NCBI)

### Examples

```bash
# Extract uORFs for TP53 gene
python main.py TP53 -o tp53_uorfs.fasta

# Extract only canonical uORFs for BRCA1 in human
python main.py BRCA1 --species human --canonical-only

# Extract uORFs for Actb in mouse with minimum length of 15 bp
python main.py Actb --species mouse --min-length 15

# Extract uORFs with custom email for NCBI
python main.py TP53 -e your.email@example.com -o tp53_uorfs.fasta
```

## Features

- **Multiple start codons**: Supports both canonical (ATG) and non-canonical start codons (CTG, GTG, TTG, ACG, ATT, ATA)
- **Genomic coordinates**: Provides exact positions of uORFs in the transcript
- **Detailed output**: FASTA headers include gene name, transcript ID, uORF number, start codon, position, length, reading frame, and type
- **BED format support**: Outputs genomic coordinates in BED format for downstream analysis
- **Multiple transcript support**: Finds uORFs across all available transcripts for a gene

## Output Format

### FASTA Output
The FASTA header contains detailed information about each uORF:
```
>{gene_name}|{transcript_id}|uORF{number}|{start_codon}|pos:{start}-{end}|len:{length}bp|frame:{frame}|type:{type}
{uORF_sequence}
```

### BED Output
The BED file includes additional columns with start/stop codons, type, and frame information:
```
# chrom	start	end	name	score	strand	start_codon	stop_codon	type	frame
```

## Dependencies

- Python >= 3.13
- Biopython >= 1.86

## Requirements

- NCBI Entrez requires a valid email address for API access
- Internet connection for fetching transcript data from NCBI

## License

MIT License

Copyright (c) 2025

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.