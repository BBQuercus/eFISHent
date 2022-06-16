"""
Find the optimal probeset with the highest coverage (fast or optimal).
"""

from typing import List, Tuple
import functools
import logging
import multiprocessing
import os

import Bio.SeqIO
import luigi
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyomo.environ as pe

from config import GeneralConfig, ProbeConfig
from config import RunConfig
from secondary_structure import SecondaryStructureFiltering
import util


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
                ("coverage", f"{util.get_gene_name()}.png"),
            ]
        }

    def run_optimal(self):
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

        with multiprocessing.Pool(GeneralConfig().threads) as pool:
            _assigned = pool.map(self.run_optimal_block, blocks)

        # Flatten the list of lists
        assigned = functools.reduce(lambda x, y: x + y, _assigned)
        return assigned

    def run_optimal_block(self, block: int) -> List[str]:
        df = self.df[self.df["block"] == block]
        sequence = list(range(df["start"].min(), df["end"].max()))
        probes = df["name"].values
        probe_starts = {k: v for _, (k, v) in df[["name", "start"]].iterrows()}
        probe_ends = {k: v for _, (k, v) in df[["name", "end"]].iterrows()}
        assigned = optimal_model(sequence, probes, probe_starts, probe_ends)
        return assigned

    def run_greedy(self):
        assigned = greedy_model(self.df)
        return assigned

    def run(self):
        sequences = list(Bio.SeqIO.parse(self.input().path, "fasta"))
        self.df = util.create_data_table(sequences)
        self.df["end"] += ProbeConfig().spacing

        if RunConfig().optimization_method == "greedy":
            assigned = self.run_greedy()
        elif RunConfig().optimization_method == "optimal":
            assigned = self.run_optimal()
        else:
            raise ValueError(
                f"Invalid optimization method: {RunConfig().optimization_method}"
                "Must be greedy or optimal."
            )

        visualize_assignment(self.df, assigned, self.output()["coverage"].path)

        candidates = [sequence for sequence in sequences if sequence.id in assigned]
        util.log_and_check_candidates(
            self.logger, "OptimizeProbeCoverage", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output()["probes"].path, format="fasta")


def is_overlapping(x: Tuple[int, int], y: Tuple[int, int]) -> bool:
    """Check if two ranges overlap."""
    x_start, x_end = x
    y_start, y_end = y
    return (
        (x_start >= y_start and x_start <= y_end)
        or (x_end >= y_start and x_end <= y_end)
        or (y_start >= x_start and y_start <= x_end)
        or (y_end >= x_start and y_end <= x_end)
    )


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
def optimal_model(
    sequence: list, probes: list, probe_starts: dict, probe_ends: dict
) -> List[str]:
    """Run the optimal mathematical model."""
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
    solver.options["tmlim"] = RunConfig().optimization_time_limit
    solver.solve(model, tee=False)

    assigned = [name for name, assign in model.assign.get_values().items() if assign]
    return assigned


def visualize_assignment(
    df: pd.DataFrame, assigned: List[str], filename: str
) -> np.ndarray:
    """Visualize the assignment as overlap matrix."""
    start = df["start"].min()
    end = df["end"].max()
    matrix = np.zeros((len(df), end - start))
    for idx, (_, row) in enumerate(df.iterrows()):
        matrix[idx, row["start"] - start : row["end"] - start + 1] = (
            2 if row["name"] in assigned else 1
        )

    plt.figure(figsize=(10, 10))
    plt.title(f"Probe Coverage from {start} to {end}")
    plt.imshow(matrix)
    plt.colorbar()
    plt.savefig(filename, dpi=600)
