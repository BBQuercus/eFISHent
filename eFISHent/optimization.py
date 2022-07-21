"""Find the optimal probeset with the highest coverage (fast or optimal)."""

from typing import List, Tuple
import functools
import logging
import multiprocessing
import os

import Bio.Seq
import Bio.SeqIO
import Bio.pairwise2
import luigi
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyomo.environ as pe

from . import util
from .config import GeneralConfig
from .config import ProbeConfig
from .config import RunConfig
from .secondary_structure import SecondaryStructureFiltering


class OptimizeProbeCoverage(luigi.Task):
    """Find the optimal probeset with the highest coverage (fast or optimal)."""

    logger = logging.getLogger("custom-logger")

    def requires(self):
        return SecondaryStructureFiltering()

    def output(self):
        return {
            name: luigi.LocalTarget(os.path.join(util.get_output_dir(), filename))
            for name, filename in [
                ("probes", f"{util.get_gene_name()}_optimal.fasta"),
                ("coverage", f"{util.get_gene_name()}_optimization.png"),
            ]
        }

    def run_optimal(self, threads: int, time_limit: int):
        """Run the optimal model by parallelizing non-overlapping segments."""
        index_block = 0
        prev_probe = self.df.iloc[0]
        blocks = []

        for _, probe in self.df.iterrows():
            if not is_overlapping(
                prev_probe[["start", "end"]], probe[["start", "end"]]
            ):
                index_block += 1
            prev_probe = probe
            blocks.append(index_block)
        self.df["block"] = blocks

        self.time_limit = time_limit
        with multiprocessing.Pool(threads) as pool:
            _assigned = pool.map(self.run_optimal_block, self.df["block"].unique())

        # Flatten the list of lists
        assigned = functools.reduce(lambda x, y: x + y, _assigned)
        return assigned

    def run_optimal_block(self, block: int) -> List[str]:
        """Run a single non-overlapping block."""
        df = self.df[self.df["block"] == block]
        return optimal_model(df, self.time_limit)

    def run_greedy(self):
        """Run sequential greedy model."""
        return greedy_model(self.df)

    def filter_binding_probes(
        self, assigned: List[str], match_percentage: float
    ) -> List[str]:
        """Remove probes of too high rev/comp sequence similarity.

        It would be cleaner to have it evaluated at run time / not have them included
        to begin with but that'd be O(n2) with potentially 1k+ probes leading to a
        runtime of multiple hours (with 20ms per `is_binding` call).
        """
        assigned_sequences = self.df.loc[
            self.df["name"].isin(assigned), "sequence"
        ].values

        # Create probe vs probe matrix that could bind
        binding_matrix = np.array(
            [
                [is_binding(x, y, match_percentage) for x in assigned_sequences]
                for y in assigned_sequences
            ]
        )  # type: ignore

        trouble_makers = np.unique(np.where(binding_matrix)[1])  # type: ignore
        self.logger.debug(f"Found {len(trouble_makers)} trouble makers to be removed.")
        assigned = [
            name for idx, name in enumerate(assigned) if idx not in trouble_makers
        ]
        return assigned

    def run(self):
        sequences = list(Bio.SeqIO.parse(self.input().path, "fasta"))
        self.df = util.create_data_table(sequences)
        self.df["end"] += ProbeConfig().spacing

        if RunConfig().optimization_method == "greedy":
            assigned = self.run_greedy()
        elif RunConfig().optimization_method == "optimal":
            assigned = self.run_optimal(
                threads=GeneralConfig().threads,
                time_limit=RunConfig().optimization_time_limit,
            )
        else:
            raise ValueError(
                f"Invalid optimization method: {RunConfig().optimization_method}"
                "Must be greedy or optimal."
            )
        if ProbeConfig().sequence_similarity > 0:
            match_percentage = ProbeConfig().sequence_similarity / 100
            assigned = self.filter_binding_probes(assigned, match_percentage)

        visualize_assignment(self.df, assigned, self.output()["coverage"].path)

        candidates = [sequence for sequence in sequences if sequence.id in assigned]
        util.log_and_check_candidates(
            self.logger, "OptimizeProbeCoverage", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output()["probes"].path, format="fasta")


def is_overlapping(x: Tuple[int, int], y: Tuple[int, int]) -> bool:
    """Check if two ranges overlap."""
    x_set = set(range(x[0], x[1] + 1))
    return len(x_set.intersection(range(y[0], y[1] + 1))) != 0


def is_binding(seq1: str, seq2: str, match_percentage: float = 0.75) -> bool:
    """Check if seq1 is similar to (rev) complement of seq2 / if would bind."""
    if 0 > match_percentage > 1:
        raise ValueError(
            "Matching percentage must be between 0 and 100. "
            f"Found `{match_percentage * 100}`"
        )
    alignments = [
        *Bio.pairwise2.align.globalxx(seq1, Bio.Seq.Seq(seq2).complement()),
        *Bio.pairwise2.align.globalxx(seq1, Bio.Seq.Seq(seq2).reverse_complement()),
    ]
    if not alignments:
        return False

    max_score = max(map(lambda x: x.score, alignments))  # type: ignore
    if (max_score / len(seq1)) <= match_percentage:
        return False

    return True


def greedy_model(df: pd.DataFrame) -> List[str]:
    """Run the greedy/fast model."""
    rec_array = df[["start", "end"]].to_records(index=False)
    probes = df["name"].values

    # Create overlap matrix of all vs all
    vect = np.vectorize(is_overlapping)
    matrix = vect(rec_array[:, None], rec_array)
    overlap = {probe: row for probe, row in zip(probes, matrix)}

    # Assign non-ovelapping probes to the first probe
    assign = {}
    assigned = [probes[0]]
    for idx, probe in enumerate(probes):
        if not overlap[assigned[-1]][idx]:
            assigned.append(probe)
        assign[probe] = probe in assigned
    return assigned


# TODO add proper mathematical description
def optimal_model(df: pd.DataFrame, time_limit: int) -> List[str]:
    """Run the optimal mathematical model."""
    # Convert dataframe to usable model inputs
    sequence = list(range(df["start"].min(), df["end"].max()))
    probes = df["name"].values
    probe_starts = {k: v for _, (k, v) in df[["name", "start"]].iterrows()}
    probe_ends = {k: v for _, (k, v) in df[["name", "end"]].iterrows()}

    # Model to make non-contiguous connections across a sequence
    # with objective to "cover" as many points in sequence as possible
    coverages = {
        p: [t for t in sequence if t >= probe_starts[p] and t <= probe_ends[p]]
        for p in probes
    }

    # Model definition
    model = pe.ConcreteModel()
    model.sequence = pe.Set(initialize=sequence)
    model.probes = pe.Set(initialize=probes)
    model.covers = pe.Set(model.probes, within=model.sequence, initialize=coverages)
    model.covers_flat_set = pe.Set(
        initialize=[(p, s) for p in probes for s in model.covers[p]]
    )
    model.assign = pe.Var(model.probes, domain=pe.Binary)
    model.covered = pe.Var(model.covers_flat_set, domain=pe.Binary)

    # Objectiv
    obj = sum(model.covered[p, s] for (p, s) in model.covers_flat_set)
    model.objective = pe.Objective(expr=obj, sense=pe.maximize)

    # Constraints
    # Selected probe must cover the associated points between start and end, if assigned
    def cover(model, p):
        return (
            sum(model.covered[p, s] for s in model.covers[p])
            == len(model.covers[p]) * model.assign[p]
        )

    # Cannot cover any point by more than 1 probe
    def over_cover(model, s):
        cov_options = [(p, s) for p in model.probes if (p, s) in model.covers_flat_set]
        if not cov_options:
            return pe.Constraint.Skip  # no possible coverages
        return sum(model.covered[p, s] for (p, s) in cov_options) <= 1

    model.C1 = pe.Constraint(model.probes, rule=cover)
    model.C2 = pe.Constraint(model.sequence, rule=over_cover)
    solver = pe.SolverFactory("glpk")
    solver.options["tmlim"] = time_limit
    solver.solve(model, tee=False)

    assigned = [name for name, assign in model.assign.get_values().items() if assign]
    return assigned


def visualize_assignment(df: pd.DataFrame, assigned: List[str], filename: str) -> None:
    """Visualize the assignment as overlap matrix."""
    start = df["start"].min()
    end = df["end"].max()
    matrix = np.zeros((len(df) + 1, end - start))
    for idx, (_, row) in enumerate(df.iterrows()):
        matrix[idx + 1, row["start"] - start : row["end"] - start + 1] = (
            2 if row["name"] in assigned else 1
        )

    plt.figure()
    plt.title(f"Probe Coverage from {start} to {end}")
    plt.imshow(matrix)
    plt.tight_layout()
    plt.savefig(filename, dpi=600, bbox_inches="tight")
