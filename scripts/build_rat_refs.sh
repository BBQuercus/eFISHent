#!/usr/bin/env bash
# Download and build Rat (GRCr8) reference indices for eFISHent.
#
# Requirements: bowtie2, jellyfish, gffread, makeblastdb, samtools
# Estimated time: ~1-2 hours (~2.9 GB genome)
# Estimated disk: ~15 GB
#
# Usage:
#   bash scripts/build_rat_refs.sh [output_dir]

set -euo pipefail

OUT="${1:-refs_rat}"
mkdir -p "$OUT"
cd "$OUT"

ENSEMBL_RELEASE=115
SPECIES="Rattus_norvegicus"
ASSEMBLY="GRCr8"

echo "=== Building Rat (${ASSEMBLY}) reference indices ==="
echo "Output directory: $(pwd)"
echo ""

# ── 1. Download genome FASTA ──
if [ -f "${SPECIES}.${ASSEMBLY}.dna.toplevel.fa" ]; then
    echo "[1/7] Genome FASTA — already downloaded"
else
    echo "[1/7] Downloading genome FASTA..."
    curl -LO "https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/${SPECIES,,}/dna/${SPECIES}.${ASSEMBLY}.dna.toplevel.fa.gz"
    echo "  Decompressing..."
    gunzip "${SPECIES}.${ASSEMBLY}.dna.toplevel.fa.gz"
    echo "  Done"
fi

# ── 2. Download GTF annotation ──
if [ -f "${SPECIES}.${ASSEMBLY}.${ENSEMBL_RELEASE}.gtf" ]; then
    echo "[2/7] GTF annotation — already downloaded"
else
    echo "[2/7] Downloading GTF annotation..."
    curl -LO "https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/${SPECIES,,}/${SPECIES}.${ASSEMBLY}.${ENSEMBL_RELEASE}.gtf.gz"
    gunzip "${SPECIES}.${ASSEMBLY}.${ENSEMBL_RELEASE}.gtf.gz"
    echo "  Done"
fi

# ── 3. Build FASTA index ──
if [ -f "${SPECIES}.${ASSEMBLY}.dna.toplevel.fa.fai" ]; then
    echo "[3/7] FASTA index — already built"
else
    echo "[3/7] Building FASTA index..."
    samtools faidx "${SPECIES}.${ASSEMBLY}.dna.toplevel.fa"
    echo "  Done"
fi

# ── 4. Build Bowtie2 index ──
if [ -f "${SPECIES}.${ASSEMBLY}.dna.toplevel.1.bt2" ]; then
    echo "[4/7] Bowtie2 index — already built"
else
    echo "[4/7] Building Bowtie2 index..."
    NPROC=1
    if command -v nproc >/dev/null 2>&1; then
        NPROC=$(nproc)
    elif command -v sysctl >/dev/null 2>&1; then
        NPROC=$(sysctl -n hw.ncpu 2>/dev/null || echo 1)
    fi
    bowtie2-build --threads "$NPROC" \
        "${SPECIES}.${ASSEMBLY}.dna.toplevel.fa" \
        "${SPECIES}.${ASSEMBLY}.dna.toplevel"
    echo "  Done"
fi

# ── 5. Build Jellyfish index (k=15) ──
if [ -f "${SPECIES}.${ASSEMBLY}.dna.toplevel_15.jf" ]; then
    echo "[5/7] Jellyfish index (k=15) — already built"
else
    echo "[5/7] Building Jellyfish index (k=15)..."
    jellyfish count -m 15 -s 3G -t 4 -C \
        "${SPECIES}.${ASSEMBLY}.dna.toplevel.fa" \
        -o "${SPECIES}.${ASSEMBLY}.dna.toplevel_15.jf"
    echo "  Done"
fi

# ── 6. Build transcriptome ──
if [ -f "transcriptome.fa" ]; then
    echo "[6/7] Transcriptome — already built"
else
    echo "[6/7] Building transcriptome with gffread..."
    if ! command -v gffread >/dev/null 2>&1; then
        echo "  ERROR: gffread not found. Install it first."
        exit 1
    fi
    gffread "${SPECIES}.${ASSEMBLY}.${ENSEMBL_RELEASE}.gtf" \
        -g "${SPECIES}.${ASSEMBLY}.dna.toplevel.fa" \
        -w transcriptome.fa
    echo "  Done"
fi

# ── 7. Build BLAST DB ──
if [ -f "transcriptome.fa.nsq" ]; then
    echo "[7/7] BLAST DB — already built"
else
    echo "[7/7] Building BLAST database..."
    if ! command -v makeblastdb >/dev/null 2>&1; then
        echo "  ERROR: makeblastdb not found. Install BLAST+ first."
        exit 1
    fi
    makeblastdb -in transcriptome.fa -dbtype nucl -parse_seqids
    echo "  Done"
fi

echo ""
echo "=== All Rat indices built ==="
echo "Total size: $(du -sh . | cut -f1)"
echo ""
echo "Files:"
ls -lh *.fa *.gtf *.bt2 *.jf transcriptome.fa.* 2>/dev/null
echo ""
echo "Next steps:"
echo "  1. Run: bash scripts/prepare_rat_hf_upload.sh $(pwd) /tmp/efishent-hf"
echo "  2. Upload via Python script"
