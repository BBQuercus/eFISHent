"""Tests for parameter presets."""

import pytest

from eFISHent.presets import (
    PRESETS,
    format_preset_list,
    get_preset,
    get_preset_names,
)


class TestPresets:
    def test_all_presets_have_description(self):
        for name, preset in PRESETS.items():
            assert "description" in preset, f"Preset '{name}' missing description"
            assert "params" in preset, f"Preset '{name}' missing params"

    def test_all_presets_have_valid_params(self):
        """All preset params should be valid CLI parameter names."""
        valid_params = {
            "min_length", "max_length", "spacing", "min_tm", "max_tm",
            "min_gc", "max_gc", "formamide_concentration", "na_concentration",
            "max_off_targets", "max_kmers", "max_deltag", "kmer_length",
            "sequence_similarity", "is_endogenous", "optimization_method",
        }
        for name, preset in PRESETS.items():
            for param in preset["params"]:
                assert param in valid_params, (
                    f"Preset '{name}' has unknown param '{param}'"
                )

    def test_get_preset(self):
        preset = get_preset("smfish")
        assert preset["params"]["min_length"] == 20

    def test_get_preset_unknown(self):
        with pytest.raises(KeyError):
            get_preset("nonexistent")

    def test_get_preset_names(self):
        names = get_preset_names()
        assert "smfish" in names
        assert "merfish" in names
        assert "strict" in names

    def test_format_preset_list(self):
        output = format_preset_list()
        assert "smfish" in output
        assert "merfish" in output
        assert "dna-fish" in output
        assert "--preset" in output

    def test_smfish_preset_values(self):
        p = get_preset("smfish")["params"]
        assert p["min_length"] == 20
        assert p["max_length"] == 20
        assert p["formamide_concentration"] == pytest.approx(10.0)

    def test_merfish_preset_values(self):
        p = get_preset("merfish")["params"]
        assert p["min_tm"] == pytest.approx(45.0)
        assert p["max_tm"] == pytest.approx(55.0)
        assert p["formamide_concentration"] == pytest.approx(30.0)

    def test_strict_vs_relaxed(self):
        strict = get_preset("strict")["params"]
        relaxed = get_preset("relaxed")["params"]
        assert strict["max_kmers"] < relaxed["max_kmers"]
        assert strict["max_off_targets"] < relaxed["max_off_targets"]
