import spatialperturb as sp


def test_assign_perturbations_sets_single_multiple_and_unassigned(demo_unannotated):
    adata = sp.assign_perturbations(
        demo_unannotated,
        barcode_columns=["CTRL_BARCODE", "LRRK2_BARCODE", "SRF_BARCODE"],
        barcode_to_perturbation={
            "CTRL_BARCODE": "control",
            "LRRK2_BARCODE": "Lrrk2",
            "SRF_BARCODE": "Srf",
        },
    )

    counts = adata.obs["perturbation"].value_counts().to_dict()
    assert counts["control"] == 4
    assert counts["Lrrk2"] == 4
    assert counts["Srf"] == 4
    assert counts["multiple"] == 2
    assert counts["unassigned"] == 6


def test_qc_perturbations_reports_target_sanity_checks(demo_adata):
    qc = sp.qc_perturbations(
        demo_adata,
        control="control",
        target_map={"Lrrk2": "LRRK2", "Srf": "SRF"},
        min_cells=3,
    )

    lrrk2_row = qc[qc["perturbation"] == "Lrrk2"].iloc[0]
    assert lrrk2_row["target_gene"] == "LRRK2"
    assert bool(lrrk2_row["valid_for_inference"])
