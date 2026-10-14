"""Index creating step for the alignment."""

import logging
import os
import subprocess

import luigi
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from . import constants, util
from .config import GeneralConfig

# Columns written to the prepared annotation parquet. Repetitive strings are
# dictionary-encoded so they load as pandas categoricals.
_GTF_STRING = pa.string()
_GTF_CATEGORY = pa.dictionary(pa.int32(), pa.string())
_GTF_SCHEMA = pa.schema(
    [
        ("seqname", _GTF_CATEGORY),
        ("feature", _GTF_CATEGORY),
        ("start", pa.int64()),
        ("end", pa.int64()),
        ("frame", _GTF_CATEGORY),
        ("gene_id", _GTF_CATEGORY),
        ("gene_name", _GTF_CATEGORY),
        ("gene_biotype", _GTF_CATEGORY),
        ("gene_type", _GTF_CATEGORY),
        ("transcript_id", _GTF_STRING),
        ("exon_number", _GTF_CATEGORY),
    ]
)
_GTF_ATTRIBUTE_COLUMNS = set(constants.GTF_ATTRIBUTES)


def read_annotation(path: str, columns: list[str] | None = None) -> pd.DataFrame:
    """Read a prepared parquet annotation, loading only the requested columns.

    Columns missing from the file are skipped, so callers can ask for e.g.
    both gene_biotype and gene_type and check which one came back.
    """
    if columns is None:
        return pd.read_parquet(path)
    available = set(pq.read_schema(path).names)
    return pd.read_parquet(path, columns=[c for c in columns if c in available])


class BuildBowtieIndex(luigi.Task):
    """Create bowtie index for a given reference genome."""

    logger = logging.getLogger("custom-logger")

    def output(self):
        return [
            luigi.LocalTarget(util.get_genome_name() + ext)
            for ext in [
                ".1.ebwt",
                ".2.ebwt",
                ".3.ebwt",
                ".4.ebwt",
                ".rev.1.ebwt",
                ".rev.2.ebwt",
            ]
        ]

    def build_bowtie_index(self, fname: str, genome: str) -> None:
        """Build bowtie index for file fname titled genome."""
        from .console import spinner

        args_bowtie = ["bowtie-build", fname, genome]
        self.logger.debug(f"Running bowtie with - {''.join(args_bowtie)}")
        with spinner("Building Bowtie index..."):
            subprocess.check_call(
                args_bowtie, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
            )

    def run(self):
        util.log_stage_start(self.logger, "BuildBowtieIndex")
        self.build_bowtie_index(
            fname=os.path.abspath(GeneralConfig().reference_genome),
            genome=util.get_genome_name(),
        )
        self.logger.info("Finished building bowtie index.")


class BuildBowtie2Index(luigi.Task):
    """Create bowtie2 index for a given reference genome."""

    logger = logging.getLogger("custom-logger")

    def output(self):
        return luigi.LocalTarget(util.get_genome_name() + ".1.bt2")

    def run(self):
        util.log_stage_start(self.logger, "BuildBowtie2Index")
        from .console import spinner

        args = [
            "bowtie2-build",
            os.path.abspath(GeneralConfig().reference_genome),
            util.get_genome_name(),
        ]
        self.logger.debug(f"Running bowtie2-build with - {''.join(args)}")
        with spinner("Building Bowtie2 index..."):
            subprocess.check_call(
                args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
            )
        self.logger.info("Finished building Bowtie2 index.")


class PrepareAnnotationFile(luigi.Task):
    """Convert gtf annotation file to parquet.

    Doesn't change any data but greatly speeds up subsequent I/O (~10x increase).
    """

    logger = logging.getLogger("custom-logger")

    def output(self):
        return luigi.LocalTarget(GeneralConfig().reference_annotation + ".parq")

    @staticmethod
    def _parse_gtf_attributes(attr_str: str) -> dict:
        """Parse GTF attribute string into a dict."""
        attrs = {}
        for item in attr_str.strip().split(";"):
            item = item.strip()
            if not item:
                continue
            parts = item.split(" ", 1)
            if len(parts) == 2:
                attrs[parts[0]] = parts[1].strip('"')
        return attrs

    @staticmethod
    def _gtf_chunk_to_table(columns: dict[str, list]) -> pa.Table:
        """Convert a chunk of parsed GTF column lists into an arrow table."""
        arrays = []
        for field in _GTF_SCHEMA:
            values = columns[field.name]
            if field.type == _GTF_CATEGORY:
                arrays.append(pa.array(values, type=_GTF_STRING).dictionary_encode())
            else:
                arrays.append(pa.array(values, type=field.type))
        return pa.Table.from_arrays(arrays, schema=_GTF_SCHEMA)

    def prepare_gtf_file(
        self, fname_input: str, fname_output: str, chunk_size: int = 250_000
    ) -> None:
        """Save the GTF rows and attributes used downstream as parquet.

        Streams the GTF in chunks of column lists, so peak memory scales with
        the kept rows in compact arrow form rather than with one Python dict
        per GTF line (which exhausted 20 GB on the human GTF). Only
        constants.GTF_FEATURES rows and constants.GTF_ATTRIBUTES are kept.
        Attribute columns absent from the whole file are dropped, matching
        the columns a full parse would have produced.
        """
        features = set(constants.GTF_FEATURES)
        names = _GTF_SCHEMA.names
        columns: dict[str, list] = {name: [] for name in names}
        seen_attributes = set()
        tables = []

        with open(fname_input) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t", 8)
                if len(parts) < 9 or parts[2] not in features:
                    continue
                attrs = self._parse_gtf_attributes(parts[8])
                columns["seqname"].append(parts[0])
                columns["feature"].append(parts[2])
                columns["start"].append(int(parts[3]))
                columns["end"].append(int(parts[4]))
                columns["frame"].append(parts[7])
                for key in _GTF_ATTRIBUTE_COLUMNS:
                    value = attrs.get(key)
                    columns[key].append(value)
                    if value is not None:
                        seen_attributes.add(key)
                if len(columns["seqname"]) >= chunk_size:
                    tables.append(self._gtf_chunk_to_table(columns))
                    columns = {name: [] for name in names}

        if columns["seqname"] or not tables:
            tables.append(self._gtf_chunk_to_table(columns))
        table = pa.concat_tables(tables)
        unused = [c for c in _GTF_ATTRIBUTE_COLUMNS if c not in seen_attributes]
        table = table.drop_columns(unused)
        pq.write_table(table, fname_output)
        self.logger.debug(
            f"Wrote {table.num_rows} GTF rows ({', '.join(table.column_names)})"
        )

    def run(self):
        self.prepare_gtf_file(GeneralConfig().reference_annotation, self.output().path)
        self.logger.info("Finished parsing GTF annotation file.")
