import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import spatialperturb as sp
from spatialperturb.figurekit import FigureKitError

SHORTCOMM_CANDIDATES = [
    ("Mast Cells", "global", "gse281048_pathway_atlas:FOS"),
    ("Basal-like Structured DCIS Cells", "Structure 1", "gse281048_pathway_atlas:CEBPB"),
    ("Dendritic Cells", "global", "gse281048_pathway_atlas:SP1"),
    ("Dendritic Cells", "Structure 2", "gse281048_pathway_atlas:MTOR"),
    ("Dendritic Cells", "Structure 3", "gse281048_pathway_atlas:RPS6KB1"),
    ("Dendritic Cells", "Structure 4", "gse281048_pathway_atlas:MAPK3"),
    ("Luminal-like Amorphous DCIS Cells", "Structure 5", "gse281048_pathway_atlas:PTGS2"),
    ("CAFs, Invasive Associated", "Structure 6", "gse281048_pathway_atlas:MAPK8"),
    ("11q13 Invasive Tumor Cells (Mitotic)", "Structure 7", "gse281048_pathway_atlas:IFNAR1"),
    ("11q13 Invasive Tumor Cells (Mitotic)", "Structure 8", "gse281048_pathway_atlas:TYK2"),
]


def test_publication_rcparams_sets_editable_font_defaults():
    sp.set_publication_rcparams()

    assert matplotlib.rcParams["pdf.fonttype"] == 42
    assert matplotlib.rcParams["ps.fonttype"] == 42
    assert matplotlib.rcParams["svg.fonttype"] == "none"
    assert matplotlib.rcParams["axes.labelsize"] == 6
    assert matplotlib.rcParams["xtick.labelsize"] == 5


def test_save_panel_writes_outputs_source_data_and_manifest(tmp_path):
    spec = sp.PanelSpec("fig1a_test_panel", "fig1", "a", 50, 35)
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    source = pd.DataFrame({"x": [0, 1], "y": [0, 1]})

    result = sp.save_panel(fig, spec, source, tmp_path)

    assert (tmp_path / "panels" / "fig1a_test_panel.pdf").exists()
    assert (tmp_path / "panels" / "fig1a_test_panel.png").exists()
    assert (tmp_path / "panels" / "fig1a_test_panel.svg").exists()
    assert (tmp_path / "source_data" / "fig1a_test_panel.tsv").exists()
    assert (tmp_path / "panel_manifest.tsv").exists()
    assert {"pdf", "png", "svg", "source_data", "manifest"}.issubset(result)
    manifest = pd.read_csv(tmp_path / "panel_manifest.tsv", sep="\t")
    assert manifest.loc[0, "panel_id"] == "fig1a_test_panel"


def test_save_panel_requires_source_data_in_strict_mode(tmp_path):
    spec = sp.PanelSpec("fig1a_missing_source", "fig1", "a", 50, 35)
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])

    with pytest.raises(FigureKitError, match="source_data is required"):
        sp.save_panel(fig, spec, None, tmp_path, strict=True)
    plt.close(fig)


def test_dense_spatial_scatter_rasterizes_large_point_cloud():
    fig, ax = plt.subplots()
    coords = np.column_stack([np.arange(50_001), np.zeros(50_001)])

    collection = sp.dense_spatial_scatter(ax, coords)

    assert collection.get_rasterized() is True
    plt.close(fig)


def _write_synthetic_nature_methods_report(report_dir):
    tables = report_dir / "tables"
    tables.mkdir(parents=True)

    obs = pd.DataFrame(
        {
            "cell_type": [
                "Mast Cells",
                "Basal-like Structured DCIS Cells",
                "Dendritic Cells",
                "Dendritic Cells",
                "CAFs, Invasive Associated",
                "11q13 Invasive Tumor Cells (Mitotic)",
            ],
            "roi": ["global", "Structure 1", "Structure 2", "Structure 3", "Structure 6", "Structure 7"],
        },
        index=[f"cell{i}" for i in range(6)],
    )
    adata = ad.AnnData(
        X=np.ones((6, 3)),
        obs=obs,
        var=pd.DataFrame(index=["GeneA", "GeneB", "GeneC"]),
        obsm={"spatial": np.array([[0, 0], [1, 0], [0, 1], [1, 1], [2, 0], [2, 1]], dtype=float)},
    )
    report_dir.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(report_dir / "input_spatial.h5ad")

    pd.DataFrame(
        [
            {
                "program": "PRPF6",
                "reference_dataset": "gse241115_breast_cropseq",
                "auroc": 0.79,
                "auprc": 0.25,
                "shuffled_auroc": 0.42,
                "delta_auroc": 0.37,
            },
            *[
                {
                    "program": program.split(":", 1)[-1],
                    "reference_dataset": "gse281048_pathway_atlas",
                    "auroc": 0.68 + idx * 0.01,
                    "auprc": 0.18,
                    "shuffled_auroc": 0.51,
                    "delta_auroc": 0.17 + idx * 0.01,
                }
                for idx, (_, _, program) in enumerate(SHORTCOMM_CANDIDATES)
            ],
        ]
    ).to_csv(tables / "reference_validation.tsv", sep="\t", index=False)

    candidate_rows = [
        {
            "grouping": "cell_type | roi",
            "group": f"cell_type={cell_type} | roi={roi}",
            "cell_type": cell_type,
            "roi": roi,
            "program": program,
            "mean_score": 0.8 - idx * 0.02,
            "z_score": 5.0 - idx * 0.1,
            "fdr": 0.249,
            "n_cells": 60,
            "is_claim_level": True,
            "claim_status": "claim_ready",
        }
        for idx, (cell_type, roi, program) in enumerate(SHORTCOMM_CANDIDATES)
    ]
    pd.DataFrame(
        [
            *candidate_rows,
            {
                "grouping": "cell_type | roi",
                "group": "cell_type=Dendritic Cells | roi=edge",
                "cell_type": "Dendritic Cells",
                "roi": "edge",
                "program": "gse241115_breast_cropseq:PRPF6",
                "mean_score": 0.5,
                "z_score": 3.0,
                "fdr": 0.2,
                "n_cells": 55,
                "is_claim_level": True,
                "claim_status": "claim_ready",
            },
        ]
    ).to_csv(tables / "calibrated_program_scores_by_group.tsv", sep="\t", index=False)

    pd.DataFrame(
        [
            {"ablation_type": "top_n", "value": 25, "score_spearman_vs_primary": 0.91},
            {"ablation_type": "top_n", "value": 50, "score_spearman_vs_primary": 1.0},
            {"ablation_type": "graph_k", "value": 15, "score_spearman_vs_primary": 0.88},
        ]
    ).to_csv(tables / "ablation_summary.tsv", sep="\t", index=False)

    pd.DataFrame(
        [
            *[
                {
                    "program": program,
                    "moran_i": 0.32 - idx * 0.01,
                    "z_score": 8.0 - idx * 0.1,
                    "fdr": 0.04,
                    "n_permutations": 25,
                }
                for idx, (_, _, program) in enumerate(SHORTCOMM_CANDIDATES)
            ],
            {
                "program": "gse241115_breast_cropseq:PRPF6",
                "moran_i": 0.22,
                "z_score": 5.0,
                "fdr": 0.05,
                "n_permutations": 25,
            },
        ]
    ).to_csv(tables / "spatial_autocorrelation.tsv", sep="\t", index=False)

    score_columns = {
        program: np.linspace(0.8 - idx * 0.03, 0.2 - idx * 0.01, obs.shape[0])
        for idx, (_, _, program) in enumerate(SHORTCOMM_CANDIDATES)
    }
    score_columns["gse241115_breast_cropseq:PRPF6"] = [0.1, 0.2, 0.6, 0.5, 0.3, 0.4]
    pd.DataFrame({"cell": obs.index, **score_columns}).to_csv(
        tables / "program_scores_cell_level.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )


def test_render_nature_methods_panels_from_synthetic_report(tmp_path):
    report_dir = tmp_path / "report"
    output_dir = tmp_path / "panels"
    _write_synthetic_nature_methods_report(report_dir)

    outputs = sp.render_nature_methods_panels(report_dir, output_dir)

    assert len(outputs) == 8
    manifest = pd.read_csv(output_dir / "panel_manifest.tsv", sep="\t")
    assert len(manifest) == 8
    assert len(list((output_dir / "panels").glob("*.pdf"))) == 8
    assert len(list((output_dir / "panels").glob("*.png"))) == 8
    assert len(list((output_dir / "panels").glob("*.svg"))) == 8
    assert len(list((output_dir / "source_data").glob("*.tsv"))) == 8
    heatmap_source = pd.read_csv(output_dir / "source_data" / "fig2c_roi_celltype_heatmap.tsv", sep="\t")
    selected_genes = set(heatmap_source["program"].astype(str).str.split(":", n=1).str[-1])
    assert {"FOS", "CEBPB", "SP1", "MTOR", "RPS6KB1", "MAPK3", "PTGS2", "MAPK8", "IFNAR1", "TYK2"}.issubset(
        selected_genes
    )
    validation_source = pd.read_csv(output_dir / "source_data" / "fig1b_reference_validation.tsv", sep="\t")
    assert {"gse241115_breast_cropseq", "gse281048_pathway_atlas"}.issubset(
        set(validation_source["reference_dataset"].astype(str))
    )
