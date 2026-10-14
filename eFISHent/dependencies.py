"""Discovery and checks of external tool dependencies."""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


def _efishent_home() -> Path:
    """Return eFISHent's per-user install directory (~/.local/efishent)."""
    return Path.home() / ".local" / "efishent"


def _ensure_deps_on_path() -> None:
    """Auto-discover eFISHent deps installed by install.sh.

    If deps were installed to ~/.local/efishent/deps/bin (the default),
    add them to PATH so users don't need to source activate.sh manually.
    """
    default_deps = _efishent_home() / "deps"
    deps_bin = default_deps / "bin"
    deps_edirect = default_deps / "edirect"

    if deps_bin.is_dir() and str(deps_bin) not in os.environ.get("PATH", ""):
        os.environ["PATH"] = f"{deps_bin}:{os.environ.get('PATH', '')}"

    if deps_edirect.is_dir() and str(deps_edirect) not in os.environ.get("PATH", ""):
        os.environ["PATH"] = f"{deps_edirect}:{os.environ.get('PATH', '')}"

    # Also set library path for jellyfish/glpk shared libs
    deps_lib = default_deps / "lib"
    if deps_lib.is_dir():
        if sys.platform == "darwin":
            lib_var = "DYLD_LIBRARY_PATH"
        else:
            lib_var = "LD_LIBRARY_PATH"
        if str(deps_lib) not in os.environ.get(lib_var, ""):
            os.environ[lib_var] = f"{deps_lib}:{os.environ.get(lib_var, '')}"


def _get_tool_version(_name: str, args: list, pattern: str = r"[0-9]+\.[0-9]+\.?[0-9]*") -> str:
    """Try to get a tool's version string."""
    import re

    try:
        result = subprocess.run(args, capture_output=True, text=True, timeout=5)
        output = result.stdout + result.stderr
        match = re.search(pattern, output)
        return match.group(0) if match else "installed"
    except Exception:
        return ""


# Dependency definitions: {name: (check_args, version_args, needed_for)}
DEPENDENCIES = {
    "bowtie": {
        "version_args": ["bowtie", "--version"],
        "needed_for": "probe alignment (legacy aligner)",
    },
    "bowtie2": {
        "version_args": ["bowtie2", "--version"],
        "needed_for": "probe alignment (default aligner)",
    },
    "jellyfish": {
        "version_args": ["jellyfish", "--version"],
        "needed_for": "k-mer counting",
    },
    "blastn": {
        "version_args": ["blastn", "-version"],
        "needed_for": "transcriptome off-target filtering (optional)",
    },
    "makeblastdb": {
        "version_args": ["makeblastdb", "-version"],
        "needed_for": "building BLAST databases for transcriptome filtering (optional)",
    },
    "dustmasker": {
        "version_args": ["dustmasker", "-version-full"],
        "needed_for": "repeat masking for off-target filtering (optional)",
    },
    "glpsol": {
        "version_args": ["glpsol", "--version"],
        "needed_for": "optimal optimization (optional if using greedy)",
    },
    "esearch": {
        "version_args": ["esearch", "-version"],
        "needed_for": "NCBI gene download (optional if providing sequence file)",
    },
    "gffread": {
        "version_args": ["gffread", "--version"],
        "needed_for": "building transcriptome from genome + GTF (optional)",
    },
}


_FOLD_PURPOSE = "secondary structure prediction"


def _dependency_status(name: str, info: dict) -> dict:
    """Look up a single external tool and report where it is found."""
    path = shutil.which(name)
    if not path:
        return {"found": False, "needed_for": info["needed_for"]}
    version = _get_tool_version(name, info["version_args"])
    return {
        "found": True,
        "path": path,
        "version": version or path,
        "needed_for": info["needed_for"],
    }


def _system_fold() -> str | None:
    """Return the RNAstructure Fold binary on PATH, if any."""
    from .secondary_structure import _is_rnastructure_fold

    system_fold = shutil.which("Fold")
    if system_fold and _is_rnastructure_fold(system_fold):
        return system_fold
    return None


def _bundled_fold() -> Path | None:
    """Return the executable Fold binary bundled for this platform, if any."""
    fold_path = Path(__file__).resolve().parent
    if sys.platform.startswith("linux"):
        fold_bin = fold_path / "Fold_linux"
    elif sys.platform == "darwin":
        fold_bin = fold_path / "Fold_osx"
    else:
        return None

    if fold_bin.exists() and os.access(fold_bin, os.X_OK):
        return fold_bin
    return None


def _fold_status() -> dict:
    """Check for Fold binary: prefer system install, fall back to bundled."""
    system_fold = _system_fold()
    if system_fold:
        return {
            "found": True,
            "version": "system",
            "path": system_fold,
            "needed_for": _FOLD_PURPOSE,
        }

    fold_bin = _bundled_fold()
    if fold_bin:
        return {
            "found": True,
            "version": "bundled",
            "path": str(fold_bin),
            "needed_for": _FOLD_PURPOSE,
        }

    return {"found": False, "needed_for": _FOLD_PURPOSE}


def check_all_dependencies() -> dict[str, dict]:
    """Check all external dependencies and return status dict."""
    results = {
        name: _dependency_status(name, info) for name, info in DEPENDENCIES.items()
    }
    results["Fold"] = _fold_status()
    return results


def _needs_ncbi_download(args: argparse.Namespace) -> bool:
    """Whether the gene sequence has to be downloaded from NCBI."""
    if args.build_indices or getattr(args, "sequence_file", None):
        return False
    return bool(getattr(args, "gene_name", None) or getattr(args, "ensembl_id", None))


def _required_tools(args: argparse.Namespace) -> list[str]:
    """List the external tools needed for the current run mode (excluding Fold)."""
    # Aligner depends on --aligner setting
    aligner = getattr(args, "aligner", "bowtie2")
    required = ["bowtie2" if aligner == "bowtie2" else "bowtie", "jellyfish"]

    # BLAST required if transcriptome filtering, repeat masking, or rDNA filter is enabled
    if getattr(args, "reference_transcriptome", "") or getattr(args, "filter_rdna_45s", True):
        required += ["blastn", "makeblastdb"]
    if getattr(args, "mask_repeats", False):
        required.append("dustmasker")

    # glpsol only needed for optimal optimization
    if getattr(args, "optimization_method", "greedy") == "optimal":
        required.append("glpsol")

    # esearch/efetch only needed for NCBI download
    if _needs_ncbi_download(args):
        required.append("esearch")

    return required


def check_required_dependencies(args: argparse.Namespace) -> list[str]:
    """Check only the dependencies required for the current run mode.

    Returns list of missing dependency names (empty = all good).
    """
    # Fold is required for the full pipeline (not for index building)
    if not args.build_indices and not (_system_fold() or _bundled_fold()):
        return ["Fold (not found in PATH and bundled binary missing)"]

    return [name for name in _required_tools(args) if not shutil.which(name)]
