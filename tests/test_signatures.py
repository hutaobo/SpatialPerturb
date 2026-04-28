import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
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


def test_calibrate_program_scores_returns_null_calibrated_table(demo_pair):
    spatial, reference = demo_pair
    programs = sp.build_reference_programs(reference, control="control", top_n=3, effect_size_only=True)
    scores = sp.score_programs(spatial, programs)
    spatial.obsm["program_scores"] = scores

    calibrated, nulls = sp.calibrate_program_scores(
        spatial,
        programs,
        n_random=3,
        n_label_shuffles=2,
        min_cells=1,
        seed=7,
        return_nulls=True,
    )

    assert not calibrated.empty
    assert not nulls.empty
    assert {"z_score", "empirical_pvalue", "fdr", "is_claim_level"}.issubset(calibrated.columns)
    assert set(nulls["null_type"]) == {"expression_matched", "group_label_shuffle"}


def test_validate_reference_programs_recovers_planted_perturbation_signal():
    genes = ["GENE_A", "GENE_B", "NOISE1", "NOISE2"]
    obs_rows = []
    matrix = []
    for perturbation, signal_gene in [("control", None), ("PERT_A", "GENE_A"), ("PERT_B", "GENE_B")]:
        for guide in ["g1", "g2"]:
            for idx in range(8):
                values = np.ones(len(genes), dtype=float)
                if signal_gene is not None:
                    values[genes.index(signal_gene)] = 12.0
                matrix.append(values)
                obs_rows.append(
                    {
                        "perturbation": perturbation,
                        "perturbation_status": "single",
                        "guide_id": f"{perturbation}_{guide}",
                        "sample": guide,
                    }
                )
    reference = ad.AnnData(np.asarray(matrix), obs=pd.DataFrame(obs_rows), var=pd.DataFrame(index=genes))
    validation = sp.validate_reference_programs(reference, control="control", split_strategy="guide", top_n=1, seed=2)

    assert {"program", "auroc", "auprc", "delta_auroc"}.issubset(validation.columns)
    assert validation["auroc"].max() > 0.8


def test_spatial_autocorrelation_detects_clustered_scores():
    adata = ad.AnnData(np.ones((6, 2)), obs=pd.DataFrame(index=[f"cell{i}" for i in range(6)]), var=pd.DataFrame(index=["G1", "G2"]))
    graph = sparse.csr_matrix(
        np.array(
            [
                [0, 1, 1, 0, 0, 0],
                [1, 0, 1, 0, 0, 0],
                [1, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 1],
                [0, 0, 0, 1, 0, 1],
                [0, 0, 0, 1, 1, 0],
            ],
            dtype=float,
        )
    )
    adata.obsp["sp_knn"] = graph
    adata.obsm["program_scores"] = pd.DataFrame({"clustered": [10, 9, 8, 0, 1, 2]}, index=adata.obs_names)

    stats = sp.spatial_autocorrelation_scores(adata, n_permutations=10, seed=3)

    assert stats.loc[0, "program"] == "clustered"
    assert stats.loc[0, "moran_i"] > 0


def test_program_redundancy_table_marks_overlapping_programs():
    table = sp.program_redundancy_table({"a": ["G1", "G2"], "b": ["G1", "G2", "G3"], "c": ["X"]}, threshold=0.5)

    assert not table.empty
    assert table.loc[(table["program_a"] == "a") & (table["program_b"] == "b"), "redundant"].iloc[0]
