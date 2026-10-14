"""Score how accessible (unpaired) each probe's target site is."""

import logging
import os
import subprocess
import sys

import pandas as pd


def fold_sequence(sequence: str) -> str | None:
    """Fold a sequence and return the dot-bracket notation.

    Uses the bundled Fold binary (RNAstructure). Returns None if
    the binary is not available or folding fails.
    """
    from pathlib import Path

    file_path = Path(__file__).resolve().parent.as_posix()
    data_table = os.path.join(file_path, "data_tables/")
    if sys.platform in ("linux", "linux2"):
        fold_path = os.path.join(file_path, "Fold_linux")
    elif sys.platform == "darwin":
        fold_path = os.path.join(file_path, "Fold_osx")
    else:
        return None

    if not os.path.isfile(fold_path):
        return None

    os.environ["DATAPATH"] = data_table
    fasta_input = f">target\n{sequence}\n"
    try:
        result = subprocess.run(
            [fold_path, "-", "-", "--bracket", "--MFE"],
            input=fasta_input,
            capture_output=True,
            text=True,
            check=True,
            timeout=30,
        )
        lines = result.stdout.strip().split("\n")
        if len(lines) >= 3:
            return lines[-1]
    except Exception as e:
        logging.getLogger("custom-logger").debug(f"Fold failed: {e}")
    return None


def compute_accessibility_scores(
    df: pd.DataFrame, target_sequence: str, window: int = 200
) -> pd.Series:
    """Compute target RNA accessibility for each probe's binding site.

    Uses local RNA folding (200nt window) around each probe binding site
    to determine what fraction of the target is unpaired (accessible).
    Returns a Series of accessibility scores in [0, 1].
    """
    scores = []
    for _, row in df.iterrows():
        start = int(row["start"])
        end = int(row["end"])

        win_start = max(0, start - window // 2)
        win_end = min(len(target_sequence), end + window // 2)
        local_seq = target_sequence[win_start:win_end]

        if len(local_seq) < 10:
            scores.append(1.0)
            continue

        try:
            dot_bracket = fold_sequence(local_seq)
            if dot_bracket is None:
                scores.append(1.0)
                continue

            probe_start_in_window = start - win_start
            probe_end_in_window = end - win_start
            binding_region = dot_bracket[probe_start_in_window:probe_end_in_window]

            if len(binding_region) == 0:
                scores.append(1.0)
                continue

            unpaired = binding_region.count(".")
            accessibility = unpaired / len(binding_region)
            scores.append(round(accessibility, 4))
        except Exception:
            scores.append(1.0)

    return pd.Series(scores, index=df.index)
