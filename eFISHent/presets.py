"""Parameter presets for common FISH protocols."""

from typing import Dict

PRESETS: Dict[str, Dict] = {
    "smfish": {
        "description": "Standard smFISH (18-22nt probes, adaptive length, 10% formamide)",
        "params": {
            "min_length": 18,
            "max_length": 22,
            "min_tm": 40.0,
            "max_tm": 60.0,
            "formamide_concentration": 10.0,
            "spacing": 2,
            "adaptive_length": True,
        },
    },
    "merfish": {
        "description": "MERFISH encoding probes (tight Tm, 30% formamide)",
        "params": {
            "min_length": 20,
            "max_length": 22,
            "min_tm": 45.0,
            "max_tm": 55.0,
            "formamide_concentration": 30.0,
            "max_gc": 65.0,
        },
    },
    "dna-fish": {
        "description": "DNA FISH (longer probes, relaxed specificity)",
        "params": {
            "min_length": 25,
            "max_length": 40,
            "min_tm": 50.0,
            "max_tm": 70.0,
            "is_endogenous": True,
        },
    },
    "strict": {
        "description": "Maximum specificity (low k-mer tolerance, low-complexity filter, CpG depletion)",
        "params": {
            "max_off_targets": 0,
            "max_kmers": 3,
            "max_deltag": -5.0,
            "filter_low_complexity": True,
            "sequence_similarity": 80,
            "max_cpg_fraction": 0.10,
        },
    },
    "relaxed": {
        "description": "Maximum probe yield (permissive thresholds + rescue filters)",
        "params": {
            "max_off_targets": 2,
            "max_kmers": 10,
            "max_deltag": -15.0,
            "mask_repeats": True,
            "off_target_min_tm": 37.0,
        },
    },
    "exogenous": {
        "description": "Exogenous genes (GFP, Renilla, reporters — no k-mer filter, strict BLAST)",
        "params": {
            "min_length": 19,
            "max_length": 22,
            "adaptive_length": True,
            "is_endogenous": False,
            "max_transcriptome_off_targets": 0,
            "min_tm": 38.0,
            "max_tm": 62.0,
            "formamide_concentration": 10.0,
            "spacing": 2,
        },
    },
}


def get_preset_names():
    """Return list of available preset names."""
    return list(PRESETS.keys())


def get_preset(name: str) -> Dict:
    """Get a preset by name. Raises KeyError if not found."""
    return PRESETS[name]


def format_preset_list() -> str:
    """Format a human-readable list of available presets."""
    lines = ["Available presets:"]
    for name, info in PRESETS.items():
        lines.append(f"  {name:<10s} {info['description']}")
    lines.append("")
    lines.append("Use: efishent --preset <name> [--other-options...]")
    lines.append("Preset values can be overridden by explicit arguments.")
    return "\n".join(lines)
