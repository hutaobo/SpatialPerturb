import spatialperturb as sp


def test_run_benchmarks_return_expected_tables(demo_pair):
    spatial, reference = demo_pair
    core = sp.run_core_benchmark(
        spatial,
        perturbations=["Lrrk2", "Srf"],
        control="control",
        target_map={"Lrrk2": "LRRK2", "Srf": "SRF"},
        sample_sizes=[2, 4],
        cell_type="neuron",
        roi="hippocampus",
    )
    concordance = sp.run_cross_platform_benchmark(
        sp.intrinsic_de(spatial, perturbation="Lrrk2", control="control", cell_type="neuron", roi="hippocampus"),
        sp.intrinsic_de(reference, perturbation="Lrrk2", control="control", cell_type="neuron", roi="hippocampus"),
    )

    assert set(core) == {"intrinsic_de", "neighbor_de", "differential_lr", "power_curve", "dataset_catalog"}
    assert not core["intrinsic_de"].empty
    assert not concordance.empty


def test_benchmark_report_and_figures(tmp_path, demo_pair):
    spatial, reference = demo_pair
    output_dir = tmp_path / "benchmark-report"
    results = sp.run_core_benchmark(
        spatial,
        perturbations=["Lrrk2"],
        control="control",
        config={"reference_adata": reference, "concordance_level": "both"},
        output_dir=output_dir,
    )

    assert (output_dir / "manifest.json").exists()
    assert (output_dir / "tables" / "intrinsic_de.tsv").exists()
    assert (output_dir / "figures" / "workflow_schema.png").exists()
    assert (output_dir / "figures" / "power_curve.png").exists()
    assert "platform_concordance" in results


def test_reference_projection_benchmark_writes_expected_outputs(tmp_path, demo_pair, demo_unannotated):
    _, reference = demo_pair
    spatial = demo_unannotated.copy()
    spatial.obs["sample"] = "xenium_demo"
    spatial.obs["cell_type"] = ["tumor"] * spatial.n_obs
    spatial.obs["roi"] = ["core"] * spatial.n_obs

    breast_reference = reference.copy()
    breast_reference.obs["cell_line"] = "HCC38"
    breast_reference.uns.setdefault("spatialperturb", {})
    breast_reference.uns["spatialperturb"]["dataset_name"] = "gse241115_breast_cropseq"

    pathway_reference = reference.copy()
    pathway_reference.obs["cell_line"] = "MCF7"
    pathway_reference.obs["stimulus"] = "TNFA"
    pathway_reference.uns.setdefault("spatialperturb", {})
    pathway_reference.uns["spatialperturb"]["dataset_name"] = "gse281048_pathway_atlas"

    output_dir = tmp_path / "reference-report"
    results = sp.run_reference_projection_benchmark(
        spatial,
        reference_datasets=["gse241115_breast_cropseq", "gse281048_pathway_atlas"],
        config={
            "reference_adatas": {
                "gse241115_breast_cropseq": breast_reference,
                "gse281048_pathway_atlas": pathway_reference,
            }
        },
        output_dir=output_dir,
    )

    assert not results["program_scores"].empty
    assert not results["neighbor_program_scores"].empty
    assert (output_dir / "manifest.json").exists()
    assert (output_dir / "input_spatial.h5ad").exists()
    assert (output_dir / "references" / "gse241115_breast_cropseq.h5ad").exists()
    assert (output_dir / "tables" / "program_scores_by_group.tsv").exists()
    assert (output_dir / "figures" / "program_scores_heatmap.png").exists()


def test_nature_methods_breast_analysis_writes_submission_outputs(tmp_path, demo_pair, demo_unannotated):
    _, reference = demo_pair
    spatial = demo_unannotated.copy()
    spatial.obs["sample"] = "xenium_demo"
    spatial.obs["cell_type"] = ["tumor"] * spatial.n_obs
    spatial.obs["roi"] = ["core"] * spatial.n_obs

    breast_reference = reference.copy()
    breast_reference.obs["cell_line"] = "HCC38"
    breast_reference.obs["guide_id"] = [
        f"{perturbation}_g{idx % 2}"
        for idx, perturbation in enumerate(breast_reference.obs["perturbation"].astype(str))
    ]
    breast_reference.uns.setdefault("spatialperturb", {})
    breast_reference.uns["spatialperturb"]["dataset_name"] = "gse241115_breast_cropseq"

    pathway_reference = reference.copy()
    pathway_reference.obs["cell_line"] = "MCF7"
    pathway_reference.obs["stimulus"] = "TNFA"
    pathway_reference.obs["guide_id"] = [
        f"{perturbation}_g{idx % 2}"
        for idx, perturbation in enumerate(pathway_reference.obs["perturbation"].astype(str))
    ]
    pathway_reference.uns.setdefault("spatialperturb", {})
    pathway_reference.uns["spatialperturb"]["dataset_name"] = "gse281048_pathway_atlas"

    output_dir = tmp_path / "nature-methods-report"
    results = sp.run_nature_methods_breast_analysis(
        spatial,
        reference_datasets=["gse241115_breast_cropseq", "gse281048_pathway_atlas"],
        config={
            "reference_adatas": {
                "gse241115_breast_cropseq": breast_reference,
                "gse281048_pathway_atlas": pathway_reference,
            },
            "n_random": 2,
            "n_label_shuffles": 1,
            "n_spatial_permutations": 2,
            "n_bootstrap": 3,
            "min_claim_cells": 1,
            "top_n": 3,
            "top_n_values": [2],
            "graph_k_values": [2],
            "seed": 5,
        },
        output_dir=output_dir,
    )

    assert results["manifest"]["benchmark"] == "nature_methods_breast_shortcomm"
    assert (output_dir / "nature_methods_summary.md").exists()
    assert (output_dir / "biological_interpretation.md").exists()
    assert (output_dir / "figures" / "main_figure_1.png").exists()
    assert (output_dir / "figures" / "main_figure_2.png").exists()
    assert (output_dir / "tables" / "reference_validation.tsv").exists()
    assert (output_dir / "tables" / "calibrated_program_scores_by_group.tsv").exists()
    assert (output_dir / "tables" / "spatial_autocorrelation.tsv").exists()
    assert (output_dir / "tables" / "ablation_summary.tsv").exists()
