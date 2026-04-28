import numpy as np
import spatialperturb as sp


def test_build_signature_matrix_creates_binary_membership_matrix():
    signature_matrix = sp.build_signature_matrix(
        {
            "ifn": ["STAT1", "IRF1", "STAT1"],
            "cell_cycle": ["TOP2A"],
        }
    )

    assert list(signature_matrix.index) == ["ifn", "cell_cycle"]
    assert list(signature_matrix.columns) == ["IRF1", "STAT1", "TOP2A"]
    assert signature_matrix.loc["ifn", "STAT1"] == 1
    assert signature_matrix.loc["ifn", "TOP2A"] == 0
    assert signature_matrix.loc["cell_cycle", "TOP2A"] == 1


def test_build_signature_matrix_handles_empty_input():
    signature_matrix = sp.build_signature_matrix({})

    assert signature_matrix.empty
    assert signature_matrix.shape == (0, 0)


def test_derive_score_and_compare_programs(demo_adata):
    de_results = sp.intrinsic_de(demo_adata, perturbation="Lrrk2", control="control", cell_type="neuron", roi="hippocampus")
    programs = sp.derive_perturbation_programs(de_results, top_n=3)
    scores = sp.score_programs(demo_adata, programs)
    concordance = sp.compare_program_concordance(scores, scores)

    assert "Lrrk2" in programs
    assert len(programs["Lrrk2"]) == 3
    assert list(scores.columns) == ["Lrrk2"]
    assert np.isclose(concordance.loc[0, "score"], 1.0)


def test_effect_size_only_de_preserves_log2fc(demo_adata):
    full = sp.intrinsic_de(demo_adata, perturbation="Lrrk2", control="control", cell_type="neuron", roi="hippocampus")
    fast = sp.intrinsic_de(
        demo_adata,
        perturbation="Lrrk2",
        control="control",
        cell_type="neuron",
        roi="hippocampus",
        effect_size_only=True,
    )

    merged = full[["gene", "log2fc"]].merge(fast[["gene", "log2fc", "pvalue", "fdr"]], on="gene", suffixes=("_full", "_fast"))
    assert np.allclose(merged["log2fc_full"], merged["log2fc_fast"])
    assert np.allclose(merged["pvalue"], 1.0)
    assert np.allclose(merged["fdr"], 1.0)


def test_build_reference_programs_and_neighbor_aggregation(demo_pair):
    spatial, reference = demo_pair
    programs, de_results = sp.build_reference_programs(
        reference,
        control="control",
        groupby="cell_line",
        top_n=3,
        direction="both",
        return_de_results=True,
    )

    scores = sp.score_programs(spatial, programs)
    spatial.obsm["program_scores"] = scores
    neighbor_scores = sp.neighbor_program_scores(spatial)
    aggregated = sp.aggregate_program_scores(spatial, scores, groupby=["cell_type", "roi"])

    assert programs
    assert not de_results.empty
    assert not neighbor_scores.empty
    assert {"program", "mean_score", "n_cells"}.issubset(aggregated.columns)
