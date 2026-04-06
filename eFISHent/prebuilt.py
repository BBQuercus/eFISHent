"""Pre-built genome index management via Hugging Face Hub.

Downloads and caches pre-built genome indices (bowtie2, jellyfish, BLAST,
transcriptome, GTF) so users can skip the index-building step with
``--genome hg38``.
"""

from typing import Dict, List, Optional, Tuple
import hashlib
import json
import logging
import os

logger = logging.getLogger("custom-logger")

# Mapping of user-friendly aliases to canonical genome identifiers.
# Only genomes with uploaded indices belong here.
GENOME_ALIASES: Dict[str, str] = {
    # Human
    "hg38": "homo_sapiens/GRCh38",
    "GRCh38": "homo_sapiens/GRCh38",
    "human": "homo_sapiens/GRCh38",
    # Mouse
    "mm39": "mus_musculus/GRCm39",
    "GRCm39": "mus_musculus/GRCm39",
    "mouse": "mus_musculus/GRCm39",
    # Yeast
    "R64": "saccharomyces_cerevisiae/R64",
    "sacCer3": "saccharomyces_cerevisiae/R64",
    "yeast": "saccharomyces_cerevisiae/R64",
    # C. elegans
    "WBcel235": "caenorhabditis_elegans/WBcel235",
    "ce11": "caenorhabditis_elegans/WBcel235",
    "worm": "caenorhabditis_elegans/WBcel235",
    "elegans": "caenorhabditis_elegans/WBcel235",
}

# Genomes we plan to support but have not uploaded indices for yet.
# Maps alias -> (display name, canonical id).
PLANNED_GENOMES: Dict[str, Tuple[str, str]] = {
    # Zebrafish
    "danRer11": ("Zebrafish (GRCz11)", "danio_rerio/GRCz11"),
    "GRCz11": ("Zebrafish (GRCz11)", "danio_rerio/GRCz11"),
    "zebrafish": ("Zebrafish (GRCz11)", "danio_rerio/GRCz11"),
    # Rat
    "mRatBN7.2": ("Rat (mRatBN7.2)", "rattus_norvegicus/mRatBN7.2"),
    "rn7": ("Rat (mRatBN7.2)", "rattus_norvegicus/mRatBN7.2"),
    "rat": ("Rat (mRatBN7.2)", "rattus_norvegicus/mRatBN7.2"),
    # Drosophila
    "dm6": ("Drosophila (BDGP6)", "drosophila_melanogaster/BDGP6"),
    "BDGP6": ("Drosophila (BDGP6)", "drosophila_melanogaster/BDGP6"),
    "fly": ("Drosophila (BDGP6)", "drosophila_melanogaster/BDGP6"),
}

WIKI_URL = "https://github.com/beichenberger/efishent/wiki/Building-Custom-Indices"

# HF Hub repository
HF_REPO_ID = "bbquercus/efishent"

# Default cache directory
DEFAULT_CACHE_DIR = os.path.join(os.path.expanduser("~"), ".local", "efishent", "indices")

# Files expected per genome
GENOME_FILES = [
    "genome.fa",
    "genome.fa.fai",
    "annotation.gtf",
    "transcriptome.fa",
    # Bowtie2 index files
    "genome.1.bt2", "genome.2.bt2", "genome.3.bt2", "genome.4.bt2",
    "genome.rev.1.bt2", "genome.rev.2.bt2",
    # Jellyfish index
    "genome_15.jf",
    # BLAST transcriptome DB (BLAST+ v5 format)
    "transcriptome.fa.nsq", "transcriptome.fa.nin", "transcriptome.fa.nhr",
    "transcriptome.fa.ndb", "transcriptome.fa.not", "transcriptome.fa.ntf",
    "transcriptome.fa.nto", "transcriptome.fa.njs",
    # Metadata
    "metadata.json",
]


def resolve_genome(alias: str) -> str:
    """Resolve a genome alias to a canonical identifier.

    Args:
        alias: User-provided genome name (e.g., "hg38", "human", "GRCh38").

    Returns:
        Canonical genome path (e.g., "homo_sapiens/GRCh38").

    Raises:
        ValueError: If the alias is not recognized.
    """
    canonical = GENOME_ALIASES.get(alias)
    if canonical is not None:
        return canonical

    # Check if this is a planned but not-yet-available genome
    planned = PLANNED_GENOMES.get(alias)
    if planned is not None:
        display_name, _canonical_id = planned
        raise ValueError(
            f"{display_name} indices are planned but not yet available. "
            f"Build your own indices using the instructions at {WIKI_URL}"
        )

    available = ", ".join(sorted(set(GENOME_ALIASES.values())))
    planned_names = ", ".join(
        sorted(set(name for name, _ in PLANNED_GENOMES.values()))
    )
    raise ValueError(
        f"Unknown genome '{alias}'. "
        f"Available genomes: {available}. "
        f"Planned (not yet available): {planned_names}. "
        f"Aliases: {', '.join(sorted(GENOME_ALIASES.keys()))}"
    )


def get_cache_dir(cache_dir: Optional[str] = None) -> str:
    """Get the genome cache directory, creating it if needed."""
    cache = cache_dir or DEFAULT_CACHE_DIR
    os.makedirs(cache, exist_ok=True)
    return cache


def get_genome_dir(genome_id: str, cache_dir: Optional[str] = None) -> str:
    """Get the directory for a specific genome's cached files."""
    return os.path.join(get_cache_dir(cache_dir), genome_id)


def is_genome_cached(genome_id: str, cache_dir: Optional[str] = None) -> bool:
    """Check if a genome's indices are already cached and valid."""
    genome_dir = get_genome_dir(genome_id, cache_dir)
    metadata_path = os.path.join(genome_dir, "metadata.json")

    if not os.path.isfile(metadata_path):
        return False

    # Check that all expected files exist
    try:
        with open(metadata_path) as f:
            metadata = json.load(f)
        for filename in metadata.get("files", {}).keys():
            if not os.path.isfile(os.path.join(genome_dir, filename)):
                return False
        return True
    except (json.JSONDecodeError, KeyError):
        return False


def verify_checksums(
    genome_id: str, cache_dir: Optional[str] = None
) -> Tuple[bool, List[str]]:
    """Verify SHA256 checksums for all cached genome files.

    Returns:
        Tuple of (all_valid, list_of_failed_files).
    """
    genome_dir = get_genome_dir(genome_id, cache_dir)
    metadata_path = os.path.join(genome_dir, "metadata.json")

    if not os.path.isfile(metadata_path):
        return False, ["metadata.json missing"]

    try:
        with open(metadata_path) as f:
            metadata = json.load(f)
    except (json.JSONDecodeError, KeyError):
        return False, ["metadata.json corrupt"]

    failed = []
    for filename, expected_hash in metadata.get("files", {}).items():
        filepath = os.path.join(genome_dir, filename)
        if not os.path.isfile(filepath):
            failed.append(filename)
            continue

        sha256 = hashlib.sha256()
        with open(filepath, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                sha256.update(chunk)
        if sha256.hexdigest() != expected_hash:
            failed.append(filename)

    return len(failed) == 0, failed


def download_genome(
    genome_id: str,
    cache_dir: Optional[str] = None,
    force: bool = False,
) -> str:
    """Download pre-built genome indices from Hugging Face Hub.

    Args:
        genome_id: Canonical genome identifier (e.g., "homo_sapiens/GRCh38").
        cache_dir: Override default cache directory.
        force: Re-download even if already cached.

    Returns:
        Path to the genome directory containing all index files.

    Raises:
        RuntimeError: If download fails.
    """
    from huggingface_hub import hf_hub_download, list_repo_tree

    genome_dir = get_genome_dir(genome_id, cache_dir)

    if not force and is_genome_cached(genome_id, cache_dir):
        logger.info(f"Genome {genome_id} already cached at {genome_dir}")
        return genome_dir

    os.makedirs(genome_dir, exist_ok=True)

    from .console import console

    console.print(f"  Downloading genome indices for [bold]{genome_id}[/bold]...")

    # List files in the genome subdirectory
    try:
        files = [
            item.rfilename
            for item in list_repo_tree(HF_REPO_ID, path_in_repo=genome_id)
            if hasattr(item, "rfilename")
        ]
    except Exception as e:
        raise RuntimeError(
            f"Could not list files for genome {genome_id} in {HF_REPO_ID}: {e}"
        )

    if not files:
        raise RuntimeError(
            f"No files found for genome {genome_id} in {HF_REPO_ID}. "
            f"This genome may not have pre-built indices yet."
        )

    # Download each file
    for filepath in files:
        filename = os.path.basename(filepath)
        console.print(f"    Downloading {filename}...", end="")
        try:
            downloaded = hf_hub_download(
                repo_id=HF_REPO_ID,
                filename=filepath,
                local_dir=get_cache_dir(cache_dir),
                local_dir_use_symlinks=False,
            )
            console.print(" [green]done[/green]")
        except Exception as e:
            console.print(f" [red]failed[/red]")
            raise RuntimeError(f"Failed to download {filepath}: {e}")

    logger.info(f"Genome {genome_id} downloaded to {genome_dir}")
    return genome_dir


def get_reference_paths(
    genome_id: str, cache_dir: Optional[str] = None
) -> Dict[str, str]:
    """Get paths to all reference files for a cached genome.

    Returns a dict with keys matching CLI parameter names:
    - reference_genome: path to genome.fa
    - reference_annotation: path to annotation.gtf
    - reference_transcriptome: path to transcriptome.fa

    Raises:
        FileNotFoundError: If the genome is not cached.
    """
    genome_dir = get_genome_dir(genome_id, cache_dir)

    if not is_genome_cached(genome_id, cache_dir):
        raise FileNotFoundError(
            f"Genome {genome_id} not cached. "
            f"Run: efishent --download-genome {genome_id}"
        )

    return {
        "reference_genome": os.path.join(genome_dir, "genome.fa"),
        "reference_annotation": os.path.join(genome_dir, "annotation.gtf"),
        "reference_transcriptome": os.path.join(genome_dir, "transcriptome.fa"),
    }


def list_available_genomes() -> List[Dict[str, str]]:
    """List all available genome aliases and their canonical IDs."""
    # Group aliases by canonical ID
    canonical_to_aliases: Dict[str, List[str]] = {}
    for alias, canonical in GENOME_ALIASES.items():
        if canonical not in canonical_to_aliases:
            canonical_to_aliases[canonical] = []
        canonical_to_aliases[canonical].append(alias)

    result = []
    for canonical, aliases in sorted(canonical_to_aliases.items()):
        result.append(
            {
                "id": canonical,
                "aliases": sorted(aliases),
                "display": f"{canonical} ({', '.join(sorted(aliases))})",
            }
        )
    return result
