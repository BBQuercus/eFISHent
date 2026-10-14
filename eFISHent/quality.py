"""Score the quality of final probes and recommend whether to use them."""

import Bio.Seq
import pandas as pd

from .config import ProbeConfig

# Quality score weights per (has_accessibility, has_binding) combination
_QUALITY_WEIGHTS: dict[tuple[bool, bool], dict[str, float]] = {
    (True, True): {
        "tm": 0.16, "gc": 0.12, "cpg": 0.05, "dg": 0.11, "ot": 0.21,
        "kmer": 0.08, "lc": 0.08, "acc": 0.09, "binding": 0.10,
    },
    (True, False): {
        "tm": 0.18, "gc": 0.14, "cpg": 0.05, "dg": 0.13, "ot": 0.23,
        "kmer": 0.09, "lc": 0.09, "acc": 0.09,
    },
    (False, True): {
        "tm": 0.18, "gc": 0.13, "cpg": 0.05, "dg": 0.13, "ot": 0.23,
        "kmer": 0.08, "lc": 0.10, "binding": 0.10,
    },
    (False, False): {
        "tm": 0.20, "gc": 0.15, "cpg": 0.05, "dg": 0.15, "ot": 0.25,
        "kmer": 0.10, "lc": 0.10,
    },
}


def _gc_score(gc: float) -> float:
    """Score GC content: full score in 45-52%, linear penalty outside, steep above 55%."""
    if 45.0 <= gc <= 52.0:
        return 1.0
    elif gc < 45.0:
        # Linear decay from 45% down to 20% (score 0)
        return max(0.0, (gc - 20.0) / 25.0)
    elif gc <= 55.0:
        # Gentle penalty 52-55%
        return 1.0 - (gc - 52.0) / 10.0
    else:
        # Steep penalty above 55% (p=0.002 for unspecific binding)
        return max(0.0, 0.7 - (gc - 55.0) / 15.0)


def _score_column(df: pd.DataFrame, column: str, score_fn) -> pd.Series | float:
    """Apply ``score_fn`` to ``df[column]``, or return a perfect score if absent."""
    if column in df.columns:
        return score_fn(df[column])
    return 1.0


class QualityScoringMixin:
    """Per-probe quality scores and PASS/FLAG/FAIL recommendations."""

    @staticmethod
    def _compute_recommendation(row: pd.Series) -> str:
        """Compute a PASS/FLAG/FAIL recommendation for a probe."""
        quality = row.get("quality", 0)
        txome_ot = row.get("txome_off_targets", 0)
        expr_risk = row.get("expression_risk", "")

        # FAIL conditions
        if quality < 30:
            return "FAIL"

        # FLAG conditions
        flags = []
        if quality < 70:
            flags.append("low_quality")
        if txome_ot > 3:
            flags.append("many_off_targets")
        if "HIGH" in str(expr_risk):
            flags.append("high_expression_risk")

        if flags:
            return f"FLAG({','.join(flags)})"

        return "PASS"

    @staticmethod
    def _compute_quality_scores(df: pd.DataFrame) -> pd.Series:
        """Compute a composite quality score (0-100) for each probe.

        Data-driven weights (n=84 probe sets):
        - Tm: deviation from set median (not config midpoint) (20%)
        - GC: optimal 45-52%, steep penalty above 55% (15%)
        - CpG: CpG dinucleotide frequency penalty (5%)
        - deltaG: less negative is better, 0 is ideal (15%)
        - off-target count: genome + transcriptome (25%)
        - kmers: lower is better (10%)
        - low-complexity: penalize low-complexity probes (10%)
        """
        scores = QualityScoringMixin._compute_component_scores(df)

        # Data-driven weighted composite
        # Accessibility gets weight when enabled (redistributed from others)
        has_accessibility = "accessibility" in df.columns and (df["accessibility"] < 1.0).any()
        has_binding = "on_target_dg" in df.columns and (df["on_target_dg"] < 0.0).any()
        weights = _QUALITY_WEIGHTS[(bool(has_accessibility), bool(has_binding))]

        quality = sum(scores[col] * weight for col, weight in weights.items())
        return (quality * 100).round(1)

    @staticmethod
    def _compute_component_scores(df: pd.DataFrame) -> pd.DataFrame:
        """Compute the individual quality score components, each in [0, 1]."""
        cfg = ProbeConfig()
        scores = pd.DataFrame(index=df.index)

        # Tm score: distance from the *actual median Tm* of the selected set
        tm_median = df["TM"].median()
        tm_range = max((cfg.max_tm - cfg.min_tm) / 2, 1)
        scores["tm"] = 1.0 - (df["TM"] - tm_median).abs().clip(upper=tm_range) / tm_range

        # GC score: data-driven optimal range 45-52%
        scores["gc"] = df["GC"].apply(_gc_score)

        # CpG score: CpG frequency penalty (p=0.004 for nuclear background)
        # 0.0 CpG = perfect, 0.10+ = zero score
        scores["cpg"] = _score_column(
            df, "cpg_fraction", lambda cpg: 1.0 - cpg.clip(upper=0.10) / 0.10
        )

        # deltaG score: 0 is ideal, more negative is worse
        dg_max = abs(cfg.max_deltag) if cfg.max_deltag < 0 else 10.0
        scores["dg"] = 1.0 - df["deltaG"].abs().clip(upper=dg_max) / dg_max

        # K-mer score: 0 is ideal, higher is worse
        kmer_max = max(cfg.max_kmers, 1)
        scores["kmer"] = 1.0 - df["kmers"].clip(upper=kmer_max) / kmer_max

        # Off-target score: combines genome off-targets and transcriptome off-targets
        ot_max = max(cfg.max_off_targets + 1, 1)
        genome_ot_score = 1.0 - df["count"].clip(lower=0, upper=ot_max) / ot_max
        txome_ot_score = pd.Series(1.0, index=df.index)
        if "txome_off_targets" in df.columns:
            txome_ot_score = 1.0 - df["txome_off_targets"].clip(lower=0, upper=5) / 5
        scores["ot"] = (genome_ot_score * 0.5 + txome_ot_score * 0.5)

        # Low-complexity score: penalize probes with low Shannon entropy
        scores["lc"] = _score_column(
            df, "low_complexity", lambda lc: 1.0 - lc.clip(upper=0.3) / 0.3
        )

        # Accessibility score: higher = more unpaired target bases = better binding
        scores["acc"] = _score_column(
            df, "accessibility", lambda acc: acc.clip(lower=0.0, upper=1.0)
        )

        # On-target duplex ΔG score: more negative = more stable binding = better
        # Typical range: -20 to -60 kcal/mol for 20-25nt probes
        # Normalize: -60 or below → 1.0, 0 → 0.0
        dg_floor = -60.0
        scores["binding"] = _score_column(
            df,
            "on_target_dg",
            lambda dg: (dg.clip(upper=0.0) / dg_floor).clip(lower=0.0, upper=1.0),
        )

        # Clip all component scores to [0, 1] for robustness
        for col in scores.columns:
            scores[col] = scores[col].clip(lower=0.0, upper=1.0)
        return scores

    @staticmethod
    def _compute_low_complexity_score(sequence: Bio.Seq.Seq) -> float:
        """Compute a low-complexity fraction for a probe sequence.

        Returns the fraction of 10bp windows with Shannon entropy < 1.5 bits.
        Higher values indicate more repetitive/low-complexity sequence.
        """
        import math
        seq_str = str(sequence)
        window = 10
        if len(seq_str) < window:
            return 0.0
        n_windows = len(seq_str) - window + 1
        low_count = 0
        for i in range(n_windows):
            w = seq_str[i:i + window]
            freqs = [w.count(b) / window for b in "ACGT"]
            entropy = -sum(f * math.log2(f) for f in freqs if f > 0)
            if entropy < 1.5:
                low_count += 1
        return low_count / n_windows
