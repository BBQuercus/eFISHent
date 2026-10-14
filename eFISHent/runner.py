"""Selection and execution of the luigi tasks for a run."""

import argparse

import luigi


def _prefill_cached_steps(add_cached_step_fn) -> None:
    """Check which pipeline steps are already cached and add them to the display."""
    from .alignment import AlignProbeCandidates
    from .basic_filtering import BasicFiltering
    from .cleanup import CleanUpOutput
    from .generate_probes import GenerateAllProbes
    from .kmers import KMerFiltering
    from .optimization import OptimizeProbeCoverage

    # Map stage names to Luigi task classes
    from .prepare_sequence import PrepareSequence
    from .rdna_filter import FilterRibosomalRNA
    from .secondary_structure import SecondaryStructureFiltering
    from .transcriptome_filter import TranscriptomeFiltering
    from .util import _get_pipeline_stages

    task_classes = {
        "PrepareSequence": PrepareSequence,
        "GenerateAllProbes": GenerateAllProbes,
        "BasicFiltering": BasicFiltering,
        "FilterRibosomalRNA": FilterRibosomalRNA,
        "AlignProbeCandidates": AlignProbeCandidates,
        "TranscriptomeFiltering": TranscriptomeFiltering,
        "KMerFiltering": KMerFiltering,
        "SecondaryStructureFiltering": SecondaryStructureFiltering,
        "OptimizeProbeCoverage": OptimizeProbeCoverage,
        "CleanUpOutput": CleanUpOutput,
    }

    # Get config-aware stage definitions and sort by step number
    stages = _get_pipeline_stages()
    numbered_steps = sorted(
        [(v["order"], v["total"], v["desc"], task_classes.get(k))
         for k, v in stages.items() if v["order"] > 0],
        key=lambda x: x[0],
    )

    for order, total, desc, task_cls in numbered_steps:
        try:
            if task_cls().complete():
                add_cached_step_fn(order, total, desc)
            else:
                break  # Stop at first incomplete step
        except Exception:
            break


def _task_classes(args: argparse.Namespace) -> list[type[luigi.Task]]:
    """Select the luigi tasks to run for the current mode."""
    # Lazy imports - only load heavy modules after arg parsing
    from .analyze import AnalyzeProbeset
    from .cleanup import CleanUpOutput
    from .indexing import BuildBowtie2Index, BuildBowtieIndex
    from .kmers import BuildJellyfishIndex

    if args.build_indices:
        if args.aligner == "bowtie2":
            index_task = BuildBowtie2Index
        else:
            index_task = BuildBowtieIndex
        return [BuildJellyfishIndex, index_task]
    if args.analyze_probeset:
        return [AnalyzeProbeset]
    return [CleanUpOutput]


def _run_tasks(tasks: list[luigi.Task], args: argparse.Namespace) -> None:
    """Run the tasks, with progress tracking unless silent."""
    if args.silent:
        luigi.build(tasks, local_scheduler=True)
        return

    from .console import add_cached_step, pipeline_progress

    if args.analyze_probeset:
        # Analysis has 7 steps
        with pipeline_progress(total_stages=7, mode="analyze"):
            luigi.build(tasks, local_scheduler=True)
    elif args.build_indices:
        with pipeline_progress(total_stages=2, mode="pipeline"):
            luigi.build(tasks, local_scheduler=True)
    else:
        from .util import _get_pipeline_total
        with pipeline_progress(total_stages=_get_pipeline_total()):
            _prefill_cached_steps(add_cached_step)
            luigi.build(tasks, local_scheduler=True)
