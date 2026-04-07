#!/usr/bin/env python3
"""
01_extract_CYCD.py

D-type Cyclin Extraction Pipeline for Phylogenetic Analysis
GmJAG1 Soybean RNA-Seq Manuscript - Task 02

Strategy:
  1. Extract 10 known Arabidopsis CYCDs by gene ID
  2. Extract 40 known soybean D-type cyclins by gene ID
  3. BLAST Arabidopsis CYCDs against Medicago, Rice, Poplar proteomes
  4. Check LxCxE motif (RBR-binding) in all sequences
  5. Combine all sequences with clean headers + metadata CSV

Dependencies: Python 3.6+, BLAST+ (loaded via module)
No BioPython required.

Usage: python scripts/01_extract_CYCD.py
  (run from cyclin_phylogeny/ directory)
"""

import os
import sys
import re
import csv
import subprocess
from pathlib import Path
from collections import defaultdict

# ============================================================
# CONFIGURATION
# ============================================================

BASE_DIR = Path(__file__).resolve().parent.parent
SEQ_DIR = BASE_DIR / "01_sequences"
ALN_DIR = BASE_DIR / "02_alignment"
REPORT_DIR = BASE_DIR / "05_reports"

# Output files
COMBINED_FASTA = SEQ_DIR / "CYCD_all_species.fa"
ATH_QUERY_FASTA = SEQ_DIR / "Ath_CYCD_queries.fa"
METADATA_CSV = REPORT_DIR / "CYCD_sequence_metadata.csv"

# BLAST parameters
EVALUE_CUTOFF = 1e-20
MIN_SEQ_LENGTH = 200    # amino acids
MIN_ALN_LENGTH = 150    # amino acids in BLAST alignment

# ============================================================
# KNOWN GENE LISTS
# ============================================================

# Arabidopsis D-type cyclins (10 genes)
ATH_CYCD = {
    "AT1G70210": "CYCD1;1",
    "AT2G22490": "CYCD2;1",
    "AT4G34160": "CYCD3;1",
    "AT5G67260": "CYCD3;2",
    "AT3G50070": "CYCD3;3",   # Known JAG target (Schiessl et al. 2014)
    "AT5G65420": "CYCD4;1",
    "AT5G10440": "CYCD4;2",
    "AT4G37630": "CYCD5;1",
    "AT4G03270": "CYCD6;1",
    "AT5G02110": "CYCD7;1",
}

# Soybean D-type cyclins (40 genes from Cyclin_comprehensive_list.csv)
# Value = putative subclade from automated annotation (to be validated by tree)
GMAX_CYCD = {
    "Glyma.01G024700": "CYCD4",  "Glyma.01G035600": "CYCD1",
    "Glyma.01G100000": "CYCD7",  "Glyma.01G189700": "CYCD3",
    "Glyma.01G193400": "CYCD5",  "Glyma.02G030200": "CYCD1",
    "Glyma.02G040500": "CYCD4",  "Glyma.02G212400": "CYCD2",
    "Glyma.03G069200": "CYCD7",  "Glyma.04G042000": "CYCD3",
    "Glyma.04G045500": "CYCD5",  "Glyma.04G092800": "CYCD2",
    "Glyma.05G097300": "CYCD3",  "Glyma.05G106800": "CYCD1",
    "Glyma.06G042900": "CYCD3",  "Glyma.06G046200": "CYCD5",
    "Glyma.06G094600": "CYCD2",  "Glyma.06G198800": "CYCD7",
    "Glyma.06G299700": "CYCD6",  "Glyma.08G280600": "CYCD2",
    "Glyma.08G291000": "CYCD1",  # JAG1 target (CYCLIN-D1-1)
    "Glyma.10G256200": "CYCD1",  "Glyma.10G263500": "CYCD3",
    "Glyma.11G048500": "CYCD5",  "Glyma.11G052500": "CYCD3",
    "Glyma.12G105600": "CYCD6",  "Glyma.12G199300": "CYCD6",
    "Glyma.12G218000": "CYCD6",  "Glyma.13G247600": "CYCD6",  # JAG1 target (CYCD6)
    "Glyma.13G303000": "CYCD6",  "Glyma.14G086200": "CYCD3",
    "Glyma.14G180300": "CYCD2",  "Glyma.15G066200": "CYCD6",
    "Glyma.17G160000": "CYCD1",  "Glyma.17G167700": "CYCD3",
    "Glyma.17G238400": "CYCD3",  "Glyma.18G133100": "CYCD1",
    "Glyma.18G145600": "CYCD2",  "Glyma.20G126700": "CYCD3",
    "Glyma.20G135400": "CYCD1",
}

JAG1_TARGETS = {"Glyma.08G291000", "Glyma.13G247600"}

# Species configuration
SPECIES_CONFIG = {
    "soybean":     {"prefix": "Gmax", "fasta": "Gmax_protein.fa"},
    "arabidopsis": {"prefix": "Ath",  "fasta": "Araport11_protein.fa"},
    "medicago":    {"prefix": "Mtr",  "fasta": "Mtruncatula_protein.fa"},
    "rice":        {"prefix": "Osa",  "fasta": "Osativa_protein.fa"},
    "poplar":      {"prefix": "Ptr",  "fasta": "Ptrichocarpa_protein.fa"},
}

# Expected D-type cyclin counts (from literature, for sanity check)
EXPECTED_COUNTS = {
    "medicago": (12, 25),
    "rice":     (10, 20),
    "poplar":   (18, 30),
}

# ============================================================
# FASTA UTILITIES (pure Python, no BioPython)
# ============================================================

def parse_fasta(filepath):
    """Parse FASTA file yielding (header, sequence) tuples."""
    header = None
    seq_parts = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.rstrip('\n\r')
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, ''.join(seq_parts)


def write_fasta(filepath, sequences, mode='w'):
    """Write sequences to FASTA. sequences = list of (header, seq)."""
    with open(filepath, mode) as f:
        for header, seq in sequences:
            f.write(f'>{header}\n')
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')


def get_gene_id_phytozome(header):
    """Extract gene-level ID from Phytozome FASTA header.
    Example: 'Glyma.01G000100.2.p pacid=... locus=Glyma.01G000100'
    Returns: 'Glyma.01G000100'
    """
    match = re.search(r'locus=(\S+)', header)
    if match:
        return match.group(1)
    # Fallback: first field, strip transcript suffix
    first = header.split()[0]
    if first.endswith('.p'):
        first = first[:-2]
    parts = first.rsplit('.', 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return first


def get_gene_id_araport(header):
    """Extract gene-level ID from Araport11 FASTA header.
    Example: 'AT1G01010.1 | Symbols:NAC001 | ...'
    Returns: 'AT1G01010'
    """
    first = header.split()[0]
    return re.sub(r'\.\d+$', '', first)


def get_first_field(header):
    """Get the first field (BLAST sequence ID) from header."""
    return header.split()[0]


def has_lxcxe_motif(sequence):
    """Check for LxCxE RBR-binding motif in protein sequence."""
    return bool(re.search(r'L.C.E', sequence))


# ============================================================
# BLAST FUNCTIONS
# ============================================================

def make_blast_db(fasta_path, db_name):
    """Create BLAST protein database."""
    cmd = ['makeblastdb', '-in', str(fasta_path), '-dbtype', 'prot',
           '-out', str(db_name)]
    print(f"  Creating BLAST db: {Path(db_name).name}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr}")
        sys.exit(1)


def run_blastp(query, db, output, evalue=1e-20):
    """Run BLASTP search."""
    fmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
    cmd = ['blastp', '-query', str(query), '-db', str(db),
           '-outfmt', fmt, '-evalue', str(evalue),
           '-max_target_seqs', '100', '-num_threads', '4',
           '-out', str(output)]
    print(f"  Running BLASTP...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr}")
        sys.exit(1)


def parse_blast_hits(blast_file):
    """Parse BLAST tabular output. Returns {sseqid: best_hit_info}."""
    best = {}
    with open(blast_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 14:
                continue
            sseqid = parts[1]
            info = {
                'qseqid':   parts[0],
                'pident':   float(parts[2]),
                'length':   int(parts[3]),
                'evalue':   float(parts[10]),
                'bitscore': float(parts[11]),
                'slen':     int(parts[13]),
            }
            # Filter
            if info['evalue'] > EVALUE_CUTOFF:
                continue
            if info['length'] < MIN_ALN_LENGTH:
                continue
            if info['slen'] < MIN_SEQ_LENGTH:
                continue
            # Keep best per subject
            if sseqid not in best or info['bitscore'] > best[sseqid]['bitscore']:
                best[sseqid] = info
    return best


def get_subclade_from_query(qseqid):
    """Extract subclade from BLAST query ID like 'Ath_AT1G70210_CYCD1;1'."""
    parts = qseqid.split('_')
    if len(parts) >= 3:
        sc = parts[-1]
        return re.sub(r';.*', '', sc)  # CYCD1;1 -> CYCD1
    return "CYCD"


# ============================================================
# MAIN PIPELINE
# ============================================================

def main():
    print("=" * 60)
    print("D-type Cyclin Extraction Pipeline")
    print("GmJAG1 Manuscript - Task 02")
    print("=" * 60)

    ALN_DIR.mkdir(exist_ok=True)
    REPORT_DIR.mkdir(exist_ok=True)

    all_seqs = []       # (clean_header, sequence)
    metadata = []       # list of dicts

    # ----------------------------------------------------------
    # STEP 1: Arabidopsis CYCDs (10 known genes)
    # ----------------------------------------------------------
    print("\n[Step 1] Extracting 10 Arabidopsis D-type cyclins...")
    ath_fasta = SEQ_DIR / "arabidopsis" / SPECIES_CONFIG["arabidopsis"]["fasta"]
    ath_found = {}
    for header, seq in parse_fasta(ath_fasta):
        gene_id = get_gene_id_araport(header)
        if gene_id in ATH_CYCD:
            ath_found[gene_id] = (header, seq)

    print(f"  Found {len(ath_found)}/{len(ATH_CYCD)}")
    missing = set(ATH_CYCD) - set(ath_found)
    if missing:
        print(f"  WARNING: Missing: {missing}")

    # Write query FASTA for BLAST + add to combined
    query_seqs = []
    for gid in sorted(ath_found):
        header, seq = ath_found[gid]
        name = ATH_CYCD[gid]
        label = f"Ath_{gid}_{name}"
        query_seqs.append((label, seq))
        all_seqs.append((label, seq))
        metadata.append({
            'species': 'Arabidopsis', 'gene_id': gid, 'label': label,
            'subclade': name, 'length_aa': len(seq),
            'LxCxE': 'YES' if has_lxcxe_motif(seq) else 'NO',
            'JAG1_target': 'NO',
            'JAG_reference': 'YES' if gid == 'AT3G50070' else 'NO',
            'method': 'known',
        })

    write_fasta(ATH_QUERY_FASTA, query_seqs)
    print(f"  Query file: {ATH_QUERY_FASTA.name}")

    # ----------------------------------------------------------
    # STEP 2: Soybean CYCDs (40 known genes)
    # ----------------------------------------------------------
    print("\n[Step 2] Extracting 40 Soybean D-type cyclins...")
    gmax_fasta = SEQ_DIR / "soybean" / SPECIES_CONFIG["soybean"]["fasta"]
    gmax_found = {}
    for header, seq in parse_fasta(gmax_fasta):
        gene_id = get_gene_id_phytozome(header)
        if gene_id in GMAX_CYCD:
            gmax_found[gene_id] = (header, seq)

    print(f"  Found {len(gmax_found)}/{len(GMAX_CYCD)}")
    missing = set(GMAX_CYCD) - set(gmax_found)
    if missing:
        print(f"  WARNING: Missing: {missing}")

    skipped = 0
    for gid in sorted(gmax_found):
        header, seq = gmax_found[gid]
        subclade = GMAX_CYCD[gid]
        is_target = gid in JAG1_TARGETS
        label = f"Gmax_{gid}_{subclade}"
        if is_target:
            label += "_JAG1target"

        if len(seq) < MIN_SEQ_LENGTH:
            print(f"  SKIPPED {gid}: {len(seq)} aa (< {MIN_SEQ_LENGTH})")
            skipped += 1
            continue

        all_seqs.append((label, seq))
        metadata.append({
            'species': 'Soybean', 'gene_id': gid, 'label': label,
            'subclade': subclade, 'length_aa': len(seq),
            'LxCxE': 'YES' if has_lxcxe_motif(seq) else 'NO',
            'JAG1_target': 'YES' if is_target else 'NO',
            'JAG_reference': 'NO',
            'method': 'known',
        })

    if skipped:
        print(f"  Skipped {skipped} sequences (too short)")

    # ----------------------------------------------------------
    # STEP 3: BLAST for Medicago, Rice, Poplar
    # ----------------------------------------------------------
    for sp in ['medicago', 'rice', 'poplar']:
        cfg = SPECIES_CONFIG[sp]
        prefix = cfg['prefix']
        sp_fasta = SEQ_DIR / sp / cfg['fasta']

        print(f"\n[Step 3] Finding D-type cyclins in {sp.capitalize()} via BLAST...")

        # Make BLAST db
        db_dir = SEQ_DIR / sp / "blastdb"
        db_dir.mkdir(exist_ok=True)
        db_path = str(db_dir / f"{prefix}_db")
        make_blast_db(sp_fasta, db_path)

        # Run BLAST
        blast_out = SEQ_DIR / sp / f"{prefix}_blast_results.txt"
        run_blastp(ATH_QUERY_FASTA, db_path, blast_out)

        # Parse hits
        hits = parse_blast_hits(blast_out)
        print(f"  Candidate hits: {len(hits)}")

        # Sanity check
        exp = EXPECTED_COUNTS.get(sp, (0, 999))
        if len(hits) < exp[0]:
            print(f"  WARNING: Fewer than expected ({exp[0]}-{exp[1]})")
        elif len(hits) > exp[1]:
            print(f"  WARNING: More than expected ({exp[0]}-{exp[1]})")
        else:
            print(f"  OK: Within expected range ({exp[0]}-{exp[1]})")

        # Extract sequences
        hit_ids = set(hits.keys())
        seen_loci = set()
        count = 0

        for header, seq in parse_fasta(sp_fasta):
            first = get_first_field(header)
            if first not in hit_ids:
                continue

            locus = get_gene_id_phytozome(header)
            if locus in seen_loci:
                continue
            seen_loci.add(locus)

            hit_info = hits[first]
            subclade = get_subclade_from_query(hit_info['qseqid'])
            label = f"{prefix}_{locus}_{subclade}"

            if len(seq) < MIN_SEQ_LENGTH:
                continue

            all_seqs.append((label, seq))
            metadata.append({
                'species': sp.capitalize(),
                'gene_id': locus,
                'label': label,
                'subclade': subclade,
                'length_aa': len(seq),
                'LxCxE': 'YES' if has_lxcxe_motif(seq) else 'NO',
                'JAG1_target': 'NO',
                'JAG_reference': 'NO',
                'method': f'BLAST (pident={hit_info["pident"]:.1f}%, '
                          f'evalue={hit_info["evalue"]:.1e})',
            })
            count += 1

        print(f"  Extracted: {count} sequences")

    # ----------------------------------------------------------
    # STEP 4: Write output
    # ----------------------------------------------------------
    print(f"\n[Step 4] Writing output files...")

    write_fasta(COMBINED_FASTA, all_seqs)
    print(f"  Combined FASTA: {COMBINED_FASTA}")
    print(f"  Total sequences: {len(all_seqs)}")

    with open(METADATA_CSV, 'w', newline='') as f:
        fields = ['species', 'gene_id', 'label', 'subclade', 'length_aa',
                  'LxCxE', 'JAG1_target', 'JAG_reference', 'method']
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(metadata)
    print(f"  Metadata CSV: {METADATA_CSV}")

    # ----------------------------------------------------------
    # SUMMARY
    # ----------------------------------------------------------
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    sp_counts = defaultdict(int)
    for row in metadata:
        sp_counts[row['species']] += 1
    for sp, n in sorted(sp_counts.items()):
        print(f"  {sp:15s}: {n:3d} D-type cyclins")
    print(f"  {'TOTAL':15s}: {len(metadata):3d} sequences")

    # LxCxE motif check
    print("\n  LxCxE motif (RBR-binding domain):")
    for row in metadata:
        if row['JAG1_target'] == 'YES' or row['JAG_reference'] == 'YES':
            status = "HAS LxCxE" if row['LxCxE'] == 'YES' else "LACKS LxCxE"
            tag = "JAG1-target" if row['JAG1_target'] == 'YES' else "JAG-ref"
            print(f"    {row['label']:50s} {status:12s} [{tag}]")

    # CYCD6 motif check
    cycd6_seqs = [r for r in metadata if 'CYCD6' in r.get('subclade', '')]
    cycd6_no_lxcxe = [r for r in cycd6_seqs if r['LxCxE'] == 'NO']
    print(f"\n  CYCD6 subclade: {len(cycd6_seqs)} sequences, "
          f"{len(cycd6_no_lxcxe)} lack LxCxE (expected for true CYCD6)")

    print("\n  Key sequences for manuscript:")
    for row in metadata:
        if row['JAG1_target'] == 'YES':
            print(f"    * {row['label']}")
    for row in metadata:
        if row['JAG_reference'] == 'YES':
            print(f"    * {row['label']}")

    print("=" * 60)
    print("Done. Next: run 02_align_and_tree.sh")


if __name__ == '__main__':
    main()
