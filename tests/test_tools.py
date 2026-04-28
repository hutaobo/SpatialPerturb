import matplotlib.pyplot as plt
import numpy as np
import spatialperturb as sp


def test_intrinsic_de_returns_tidy_results(demo_adata):
    result = sp.intrinsic_de(demo_adata, perturbation="Lrrk2", control="control", cell_type="neuron", roi="hippocampus")

    assert {"gene", "log2fc", "fdr", "effect_type"}.issubset(result.columns)
    assert result.iloc[0]["effect_type"] == "intrinsic"
    assert "STAT1" in set(result["gene"])


def test_neighbor_de_returns_tidy_results(demo_adata):
    result = sp.neighbor_de(demo_adata, perturbation="Lrrk2", control="control", cell_type="neuron", roi="hippocampus")

    assert {"gene", "log2fc", "shared_neighbors_removed"}.issubset(result.columns)
    assert result.iloc[0]["effect_type"] == "neighbor"


def test_differential_lr_platform_concordance_and_power_curve(demo_pair):
    spatial, reference = demo_pair
    spatial_de = sp.intrinsic_de(spatial, perturbation="Lrrk2", control="control", cell_type="neuron", roi="hippocampus")
    reference_de = sp.intrinsic_de(reference, perturbation="Lrrk2", control="control", cell_type="neuron", roi="hippocampus")
    lr = sp.differential_lr(spatial, perturbation="Lrrk2", control="control")
    concordance = sp.platform_concordance(spatial_de, reference_de)
    power = sp.power_curve(spatial, perturbation="Lrrk2", control="control", sample_sizes=[2, 4], n_boot=5)

    assert not lr.empty
    assert set(lr.columns) >= {"ligand", "receptor", "diff_score"}
    assert concordance.iloc[0]["perturbation"] == "Lrrk2"
    assert np.all((power["power"].dropna() >= 0) & (power["power"].dropna() <= 1))


def test_pseudobulk_and_custom_network_modes(demo_pair):
    spatial, reference = demo_pair
    intrinsic = sp.intrinsic_de(
        spatial,
        perturbation="Lrrk2",
        control="control",
        method="pseudobulk",
        sample_col="sample",
        cell_type="neuron",
        roi="hippocampus",
    )
    neighbor = sp.neighbor_de(
        spatial,
        perturbation="Lrrk2",
        control="control",
        method="pseudobulk",
        sample_col="sample",
        aggregate="pseudobulk",
        cell_type="neuron",
        roi="hippocampus",
    )
    network = sp.differential_lr(
        spatial,
        perturbation="Lrrk2",
        control="control",
        lr_network=sp.tl._FALLBACK_LR_NETWORK.assign(ligand="NOT_A_GENE"),
    )
    custom_lr = sp.differential_lr(
        spatial,
        perturbation="Lrrk2",
        control="control",
        lr_network=sp.tl._FALLBACK_LR_NETWORK.tail(1).copy(),
    )
    concordance = sp.platform_concordance(
        intrinsic,
        sp.intrinsic_de(
            reference,
            perturbation="Lrrk2",
            control="control",
            method="pseudobulk",
            sample_col="sample",
            cell_type="neuron",
            roi="hippocampus",
        ),
        level="both",
    )

    assert set(intrinsic["method"]) == {"pseudobulk"}
    assert set(neighbor["method"]) == {"pseudobulk"}
    assert set(neighbor["aggregate"]) == {"pseudobulk"}
    assert network.empty
    assert not custom_lr.empty
    assert {"program_jaccard", "top_gene_jaccard"}.issubset(concordance.columns)


def test_plotting_helpers_return_axes(demo_pair):
    spatial, reference = demo_pair
    intrinsic = sp.intrinsic_de(spatial, perturbation="Lrrk2", control="control", cell_type="neuron", roi="hippocampus")
    neighbor = sp.neighbor_de(spatial, perturbation="Lrrk2", control="control", cell_type="neuron", roi="hippocampus")
    lr = sp.differential_lr(spatial, perturbation="Lrrk2", control="control")
    concordance = sp.platform_concordance(
        intrinsic,
        sp.intrinsic_de(reference, perturbation="Lrrk2", control="control", cell_type="neuron", roi="hippocampus"),
    )
    power = sp.power_curve(spatial, perturbation="Lrrk2", control="control", sample_sizes=[2, 4], n_boot=5)

    axes = [
        sp.pl.barcode_spread(spatial),
        sp.pl.own_vs_neighbor(intrinsic, neighbor, perturbation="Lrrk2"),
        sp.pl.lr_pairs(lr),
        sp.pl.lr_map(spatial, lr, perturbation="Lrrk2"),
        sp.pl.platform_concordance(concordance),
        sp.pl.power_curve(power),
    ]

    assert all(hasattr(ax, "figure") for ax in axes)
    plt.close("all")
