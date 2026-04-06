#!/usr/bin/env bash
# Prepare yeast (S. cerevisiae R64-1-1) refs for HuggingFace upload.
# Creates the expected directory structure with standardized filenames.
#
# Usage:
#   bash scripts/prepare_yeast_hf_upload.sh ./refs_yeast /tmp/efishent-hf

set -euo pipefail

REFS_DIR="${1:?Usage: $0 <refs_dir> <hf_repo_dir>}"
HF_DIR="${2:?Usage: $0 <refs_dir> <hf_repo_dir>}"

GENOME_DIR="$HF_DIR/saccharomyces_cerevisiae/R64"
mkdir -p "$GENOME_DIR"

echo "Copying and renaming files from $REFS_DIR to $GENOME_DIR..."

# Genome FASTA + index
cp "$REFS_DIR/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa" "$GENOME_DIR/genome.fa"
cp "$REFS_DIR/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.fai" "$GENOME_DIR/genome.fa.fai"

# Bowtie2 indices
for i in 1 2 3 4; do
    cp "$REFS_DIR/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.$i.bt2" "$GENOME_DIR/genome.$i.bt2"
done
for i in 1 2; do
    cp "$REFS_DIR/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.rev.$i.bt2" "$GENOME_DIR/genome.rev.$i.bt2"
done

# GTF annotation
cp "$REFS_DIR/Saccharomyces_cerevisiae.R64-1-1.115.gtf" "$GENOME_DIR/annotation.gtf"

# Jellyfish index (k=15)
cp "$REFS_DIR/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel_15.jf" "$GENOME_DIR/genome_15.jf"

# Transcriptome + BLAST DB
cp "$REFS_DIR/transcriptome.fa" "$GENOME_DIR/transcriptome.fa"
for ext in ndb nhr nin njs nog nos not nsq ntf nto; do
    cp "$REFS_DIR/transcriptome.fa.$ext" "$GENOME_DIR/transcriptome.fa.$ext"
done

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
    'genome': 'Saccharomyces cerevisiae (R64-1-1)',
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
echo "  git commit -m 'Add R64 S. cerevisiae (yeast) genome indices'"
echo "  git push"
