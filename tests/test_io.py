import json

import pandas as pd
import spatialperturb as sp


def test_from_tables_builds_schema_compliant_anndata():
    expression = pd.DataFrame(
        [[1, 2], [3, 4]],
        index=["cell_a", "cell_b"],
        columns=["GeneA", "GeneB"],
    )
    obs = pd.DataFrame({"cell_type": ["neuron", "astrocyte"], "roi": ["hippocampus", "cortex"], "x": [0, 1], "y": [1, 0]}, index=expression.index)

    adata = sp.from_tables(expression, obs=obs, metadata={"platform": "demo"})

    sp.validate_spatialperturb_schema(adata)
    assert adata.uns["spatialperturb"]["platform"] == "demo"
    assert list(adata.var_names) == ["GeneA", "GeneB"]


def test_read_directory_backends_support_simple_csv_exports(tmp_path):
    expression = pd.DataFrame(
        [[1, 0], [0, 2]],
        index=["cell_a", "cell_b"],
        columns=["GeneA", "GeneB"],
    )
    cells = pd.DataFrame({"cell_type": ["neuron", "astrocyte"], "roi": ["hippocampus", "cortex"], "x": [0, 1], "y": [1, 0]}, index=expression.index)

    xenium_dir = tmp_path / "xenium"
    stereoseq_dir = tmp_path / "stereoseq"
    xenium_dir.mkdir()
    stereoseq_dir.mkdir()

    expression.to_csv(xenium_dir / "counts.csv")
    cells.to_csv(xenium_dir / "cells.csv")
    expression.to_csv(stereoseq_dir / "counts.csv")
    cells.to_csv(stereoseq_dir / "cells.csv")

    xenium = sp.read_xenium(xenium_dir)
    stereoseq = sp.read_stereoseq(stereoseq_dir)

    assert xenium.shape == (2, 2)
    assert stereoseq.shape == (2, 2)


def test_read_xenium_supports_cell_groups_and_roi_annotations(tmp_path):
    expression = pd.DataFrame(
        [[3, 1], [0, 4]],
        index=["cell_a", "cell_b"],
        columns=["GeneA", "GeneB"],
    )
    cells = pd.DataFrame(
        {
            "cell_type": ["unknown", "unknown"],
            "roi": ["global", "global"],
            "x_centroid": [0.5, 5.0],
            "y_centroid": [0.5, 5.0],
        },
        index=expression.index,
    )
    cell_groups = pd.DataFrame(
        {
            "cell_id": ["cell_a", "cell_b"],
            "group": ["Tumor", "Immune"],
            "color": ["#ff0000", "#00ff00"],
        }
    )
    roi_geojson = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {"assigned_structure": "Tumor core"},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[[0.0, 0.0], [1.5, 0.0], [1.5, 1.5], [0.0, 1.5], [0.0, 0.0]]],
                },
            }
        ],
    }

    xenium_dir = tmp_path / "xenium"
    xenium_dir.mkdir()
    expression.to_csv(xenium_dir / "counts.csv")
    cells.to_csv(xenium_dir / "cells.csv")
    cell_groups_path = tmp_path / "cell_groups.csv"
    cell_groups.to_csv(cell_groups_path, index=False)
    roi_path = tmp_path / "roi.geojson"
    roi_path.write_text(json.dumps(roi_geojson), encoding="utf-8")

    adata = sp.read_xenium(
        xenium_dir,
        cell_group_path=cell_groups_path,
        roi_geojson_path=roi_path,
        sample_name="breast1",
    )

    assert list(adata.obs["cell_type"]) == ["Tumor", "Immune"]
    assert list(adata.obs["roi"]) == ["Tumor core", "global"]
    assert list(adata.obs["sample"]) == ["breast1", "breast1"]
    assert adata.uns["spatialperturb"]["cell_group_annotation"]["matched_cells"] == 2
    assert adata.uns["spatialperturb"]["roi_annotation"]["matched_cells"] == 1
