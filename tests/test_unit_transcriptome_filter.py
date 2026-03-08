import shutil

import pandas as pd

BLAST_AVAILABLE = shutil.which("blastn") is not None


def test_blast_output_parsing():
    """Parse sample BLAST -outfmt 6 output."""
    from io import StringIO

    columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch",
        "gapopen", "qstart", "qend", "sstart", "send",
        "evalue", "bitscore",
    ]
    sample = (
        "candidate-0-100\tGENE1\t90.0\t20\t2\t0\t1\t20\t100\t119\t0.001\t30\n"
        "candidate-0-100\tGENE2\t80.0\t18\t3\t1\t1\t18\t200\t217\t0.01\t25\n"
        "candidate-1-200\tGENE1\t95.0\t22\t1\t0\t1\t22\t300\t321\t0.0001\t35\n"
    )
    df = pd.read_csv(StringIO(sample), sep="\t", header=None, names=columns)
    assert len(df) == 3
    assert df.iloc[0]["qseqid"] == "candidate-0-100"
    assert df.iloc[0]["pident"] == 90.0


def test_identity_threshold_filtering():
    """Hits below threshold are ignored."""
    from io import StringIO

    columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch",
        "gapopen", "qstart", "qend", "sstart", "send",
        "evalue", "bitscore",
    ]
    sample = (
        "probe1\tGENE1\t90.0\t20\t2\t0\t1\t20\t100\t119\t0.001\t30\n"
        "probe1\tGENE2\t60.0\t18\t7\t1\t1\t18\t200\t217\t0.5\t10\n"
    )
    df = pd.read_csv(StringIO(sample), sep="\t", header=None, names=columns)

    # Apply threshold at 75%
    threshold = 75.0
    min_match_len = 15
    df["effective_len"] = df["length"] - df["gapopen"]
    df_hits = df[
        (df["effective_len"] >= min_match_len) & (df["pident"] >= threshold)
    ]
    assert len(df_hits) == 1
    assert df_hits.iloc[0]["sseqid"] == "GENE1"


def test_target_gene_excluded():
    """Self-hits to target gene don't count as off-targets."""
    from io import StringIO

    columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch",
        "gapopen", "qstart", "qend", "sstart", "send",
        "evalue", "bitscore",
    ]
    sample = (
        "probe1\trenilla_transcript\t100.0\t25\t0\t0\t1\t25\t1\t25\t0.0\t50\n"
        "probe1\tother_gene\t85.0\t20\t3\t0\t1\t20\t100\t119\t0.01\t30\n"
    )
    df = pd.read_csv(StringIO(sample), sep="\t", header=None, names=columns)

    gene_name = "renilla"
    off_targets = df[~df["sseqid"].str.lower().str.contains(gene_name, na=False)]
    assert len(off_targets) == 1
    assert off_targets.iloc[0]["sseqid"] == "other_gene"
