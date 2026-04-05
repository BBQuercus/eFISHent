#!/usr/bin/env bash
# Prepare refs/ directory contents for HuggingFace upload.
# Creates the expected directory structure with standardized filenames.
#
# Usage:
#   1. Clone the HF repo:  git clone https://huggingface.co/datasets/bbquercus/efishent
#   2. Run this script:     bash scripts/prepare_hf_upload.sh ./refs ./efishent
#   3. Push:                cd efishent && git add . && git commit -m "Add hg38 indices" && git push

set -euo pipefail

REFS_DIR="${1:?Usage: $0 <refs_dir> <hf_repo_dir>}"
HF_DIR="${2:?Usage: $0 <refs_dir> <hf_repo_dir>}"

GENOME_DIR="$HF_DIR/homo_sapiens/GRCh38"
mkdir -p "$GENOME_DIR"

echo "Copying and renaming files from $REFS_DIR to $GENOME_DIR..."

# Genome FASTA + index
cp "$REFS_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa" "$GENOME_DIR/genome.fa"
cp "$REFS_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai" "$GENOME_DIR/genome.fa.fai"

# Bowtie2 indices
for i in 1 2 3 4; do
    cp "$REFS_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.$i.bt2" "$GENOME_DIR/genome.$i.bt2"
done
for i in 1 2; do
    cp "$REFS_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.rev.$i.bt2" "$GENOME_DIR/genome.rev.$i.bt2"
done

# GTF annotation
cp "$REFS_DIR/Homo_sapiens.GRCh38.115.gtf" "$GENOME_DIR/annotation.gtf"

# Jellyfish index (k=15)
cp "$REFS_DIR/Homo_sapiens.GRCh38.dna.primary_assembly_15.jf" "$GENOME_DIR/genome_15.jf"

# Transcriptome + BLAST DB
cp "$REFS_DIR/transcriptome.fa" "$GENOME_DIR/transcriptome.fa"
for ext in ndb nhr nin njs not nsq ntf nto; do
    cp "$REFS_DIR/transcriptome.fa.$ext" "$GENOME_DIR/transcriptome.fa.$ext"
done

# Also copy the esther/stellaris comparison data
mkdir -p "$HF_DIR/analysis"
cp "$REFS_DIR/esther_probes.fasta" "$HF_DIR/analysis/"
cp "$REFS_DIR/stellaris_probes.fasta" "$HF_DIR/analysis/"
cp "$REFS_DIR/esther_vs_stellaris.csv" "$HF_DIR/analysis/"
cp "$REFS_DIR/probe_comparison_detailed.csv" "$HF_DIR/analysis/"
cp "$REFS_DIR/esther_vs_stellaris_Esther (21nt)_blast.tsv" "$HF_DIR/analysis/"
cp "$REFS_DIR/esther_vs_stellaris_Stellaris (20nt)_blast.tsv" "$HF_DIR/analysis/"

# Generate metadata.json with SHA256 checksums
echo "Computing checksums..."
cd "$GENOME_DIR"
python3 -c "
import hashlib, json, os

files = {}
for f in sorted(os.listdir('.')):
    if f == 'metadata.json' or os.path.isdir(f):
        continue
    sha = hashlib.sha256()
    with open(f, 'rb') as fh:
        for chunk in iter(lambda: fh.read(8192), b''):
            sha.update(chunk)
    files[f] = sha.hexdigest()
    print(f'  {f}: {sha.hexdigest()[:16]}...')

metadata = {
    'genome': 'Homo sapiens (GRCh38)',
    'source': 'Ensembl release 115',
    'bowtie2_version': '2.5.x',
    'jellyfish_kmer': 15,
    'blast_version': 'BLAST+ 2.16.x',
    'files': files,
}
with open('metadata.json', 'w') as f:
    json.dump(metadata, f, indent=2)
print('  metadata.json written')
"

echo ""
echo "Done! Files prepared at: $GENOME_DIR"
echo "Total size: $(du -sh "$GENOME_DIR" | cut -f1)"
echo ""
echo "Next steps:"
echo "  cd $HF_DIR"
echo "  git add ."
echo "  git commit -m 'Add hg38/GRCh38 genome indices'"
echo "  git push"
