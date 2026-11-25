"""Pre-built genome index management via Hugging Face Hub.

Downloads and caches pre-built genome indices (bowtie2, jellyfish, BLAST,
transcriptome, GTF) so users can skip the index-building step with
``--genome hg38``.
"""

import hashlib
import json
import logging
import os

logger = logging.getLogger("custom-logger")

# Canonical genome identifiers (species/assembly) of the hosted indices.
HUMAN = "homo_sapiens/GRCh38"
MOUSE = "mus_musculus/GRCm39"
YEAST = "saccharomyces_cerevisiae/R64"
WORM = "caenorhabditis_elegans/WBcel235"
FLY = "drosophila_melanogaster/BDGP6"
ZEBRAFISH = "danio_rerio/GRCz11"
RAT = "rattus_norvegicus/GRCr8"

# Mapping of user-friendly aliases to canonical genome identifiers.
# Only genomes with uploaded indices belong here.
GENOME_ALIASES: dict[str, str] = {
    "hg38": HUMAN,
    "GRCh38": HUMAN,
    "human": HUMAN,
    "mm39": MOUSE,
    "GRCm39": MOUSE,
    "mouse": MOUSE,
    "R64": YEAST,
    "sacCer3": YEAST,
    "yeast": YEAST,
    "WBcel235": WORM,
    "ce11": WORM,
    "worm": WORM,
    "elegans": WORM,
    "dm6": FLY,
    "BDGP6": FLY,
    "fly": FLY,
    "danRer11": ZEBRAFISH,
    "GRCz11": ZEBRAFISH,
    "zebrafish": ZEBRAFISH,
    "GRCr8": RAT,
    "rn7": RAT,
    "rat": RAT,
}

# Genomes we plan to support but have not uploaded indices for yet.
# Maps alias -> (display name, canonical id).
PLANNED_GENOMES: dict[str, tuple[str, str]] = {}

WIKI_URL = "https://github.com/beichenberger/efishent/wiki/Building-Custom-Indices"

# HF Hub repository
HF_REPO_ID = "bbquercus/efishent"

# Default cache directory
DEFAULT_CACHE_DIR = os.path.join(os.path.expanduser("~"), ".local", "efishent", "indices")

METADATA_FILENAME = "metadata.json"


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
        sorted({name for name, _ in PLANNED_GENOMES.values()})
    )
    raise ValueError(
        f"Unknown genome '{alias}'. "
        f"Available genomes: {available}. "
        f"Planned (not yet available): {planned_names}. "
        f"Aliases: {', '.join(sorted(GENOME_ALIASES.keys()))}"
    )


def get_cache_dir(cache_dir: str | None = None) -> str:
    """Get the genome cache directory, creating it if needed."""
    cache = cache_dir or os.environ.get("EFISHENT_INDEX_DIR") or DEFAULT_CACHE_DIR
    os.makedirs(cache, exist_ok=True)
    return cache


def get_genome_dir(genome_id: str, cache_dir: str | None = None) -> str:
    """Get the directory for a specific genome's cached files."""
    return os.path.join(get_cache_dir(cache_dir), genome_id)


def is_genome_cached(genome_id: str, cache_dir: str | None = None) -> bool:
    """Check if a genome's indices are already cached and valid."""
    genome_dir = get_genome_dir(genome_id, cache_dir)
    metadata_path = os.path.join(genome_dir, METADATA_FILENAME)

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


def _sha256(path: str) -> str:
    """Return the SHA256 hex digest of a file, read in 1 MiB chunks."""
    sha256 = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def _listed_checksums(genome_dir: str) -> dict[str, str]:
    """Return the filename -> SHA256 map in a genome's metadata.json.

    Raises FileNotFoundError if metadata.json is missing, ValueError or
    AttributeError if it is corrupt.
    """
    with open(os.path.join(genome_dir, METADATA_FILENAME)) as f:
        return json.load(f).get("files", {})


def verify_checksums(
    genome_id: str, cache_dir: str | None = None
) -> tuple[bool, list[str]]:
    """Verify SHA256 checksums for all cached genome files.

    Checks every file listed under "files" in metadata.json (written by
    scripts/prepare_*_hf_upload.sh). Files on disk without an entry there are
    not checked; download_genome warns about them instead of failing.

    Returns:
        Tuple of (all_valid, list_of_failed_files). A missing or corrupt
        metadata.json is reported as a single "metadata.json missing/corrupt"
        entry.
    """
    genome_dir = get_genome_dir(genome_id, cache_dir)

    try:
        expected = _listed_checksums(genome_dir)
    except FileNotFoundError:
        return False, [f"{METADATA_FILENAME} missing"]
    except (OSError, ValueError, AttributeError):
        return False, [f"{METADATA_FILENAME} corrupt"]

    failed = [
        filename
        for filename, expected_hash in expected.items()
        if not os.path.isfile(os.path.join(genome_dir, filename))
        or _sha256(os.path.join(genome_dir, filename)) != expected_hash
    ]
    return len(failed) == 0, failed


def download_genome(
    genome_id: str,
    cache_dir: str | None = None,
    force: bool = False,
) -> str:
    """Download pre-built genome indices from Hugging Face Hub.

    Args:
        genome_id: Canonical genome identifier (e.g., "homo_sapiens/GRCh38").
        cache_dir: Override default cache directory.
        force: Re-download even if already cached.

    Returns:
        Path to the genome directory containing all index files.

    After the download, every file is checked against the SHA256 checksums in
    metadata.json. Files that fail are deleted, so the genome does not count
    as cached and the next run downloads them again. Downloaded files without
    a checksum entry are kept with a warning.

    Raises:
        RuntimeError: If the download or the checksum check fails.
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
            for item in list_repo_tree(HF_REPO_ID, path_in_repo=genome_id, repo_type="dataset")
            if hasattr(item, "rfilename")
        ]
    except Exception as e:
        raise RuntimeError(
            f"Could not list files for genome {genome_id} in {HF_REPO_ID}: {e}"
        ) from e

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
            hf_hub_download(
                repo_id=HF_REPO_ID,
                filename=filepath,
                repo_type="dataset",
                revision="main",
                local_dir=get_cache_dir(cache_dir),
                local_dir_use_symlinks=False,
            )
            console.print(" [green]done[/green]")
        except Exception as e:
            console.print(" [red]failed[/red]")
            raise RuntimeError(f"Failed to download {filepath}: {e}") from e

    console.print("    Verifying checksums...", end="")
    valid, failed = verify_checksums(genome_id, cache_dir)
    if not valid:
        console.print(" [red]failed[/red]")
        for filename in failed:
            path = os.path.join(genome_dir, filename)
            if os.path.isfile(path):
                os.remove(path)
        raise RuntimeError(
            f"Checksum check failed for genome {genome_id}: {', '.join(failed)}. "
            "The bad files were deleted. Run the command again to re-download them."
        )
    console.print(" [green]done[/green]")

    unchecked = sorted(
        {os.path.basename(p) for p in files}
        - _listed_checksums(genome_dir).keys()
        - {METADATA_FILENAME}
    )
    if unchecked:
        console.print(
            f"    [yellow]No checksum in {METADATA_FILENAME} for "
            f"{', '.join(unchecked)}; kept unchecked.[/yellow]"
        )

    logger.info(f"Genome {genome_id} downloaded to {genome_dir}")
    return genome_dir


def get_reference_paths(
    genome_id: str, cache_dir: str | None = None
) -> dict[str, str]:
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


def try_download_supplemental_index(
    genome_id: str,
    kmer_length: int,
    cache_dir: str | None = None,
) -> bool:
    """Try to download a pre-built jellyfish index for a non-default kmer_length.

    The standard download includes only genome_15.jf. For other kmer lengths
    (e.g. 20), the index can be hosted separately on HF and fetched on demand
    here so users never have to build from the full genome locally.

    Returns True if the file is now present (either already there or freshly
    downloaded), False if it is unavailable (not yet hosted or network error)
    or failed its checksum. A file that fails is deleted, so the caller builds
    the index locally. A file without a checksum in metadata.json is kept
    with a warning.
    """
    filename = f"genome_{kmer_length}.jf"
    genome_dir = get_genome_dir(genome_id, cache_dir)
    out_path = os.path.join(genome_dir, filename)

    if os.path.isfile(out_path):
        return True

    hf_path = f"{genome_id}/{filename}"
    try:
        from huggingface_hub import hf_hub_download
        hf_hub_download(
            repo_id=HF_REPO_ID,
            filename=hf_path,
            repo_type="dataset",
            local_dir=get_cache_dir(cache_dir),
            local_dir_use_symlinks=False,
        )
    except Exception:
        return False
    if not os.path.isfile(out_path):
        return False

    from .console import console

    try:
        expected = _listed_checksums(genome_dir).get(filename)
    except (OSError, ValueError, AttributeError):
        expected = None
    if expected is None:
        console.print(
            f"\n  [yellow]No checksum in {METADATA_FILENAME} for {filename}; "
            "using it unchecked.[/yellow]"
        )
        return True
    if _sha256(out_path) != expected:
        os.remove(out_path)
        console.print(
            f"\n  [red]Checksum mismatch for {filename}; deleted the download.[/red]"
        )
        return False
    return True


def list_available_genomes() -> list[dict[str, str]]:
    """List all available genome aliases and their canonical IDs."""
    # Group aliases by canonical ID
    canonical_to_aliases: dict[str, list[str]] = {}
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
