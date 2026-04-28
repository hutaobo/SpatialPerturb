import spatialperturb as sp


def test_build_spatial_graph_creates_knn_graph(demo_unannotated):
    adata = sp.build_spatial_graph(demo_unannotated, mode="knn", k=2)

    assert "sp_knn" in adata.obsp
    assert adata.obsp["sp_knn"].nnz > 0


def test_collect_neighbors_excludes_perturbed_by_default(demo_adata):
    lrrk2_cells = demo_adata.obs_names[demo_adata.obs["perturbation"] == "Lrrk2"]
    neighbors = sp.collect_neighbors(demo_adata, cells=lrrk2_cells[:1])

    assert neighbors
    first_neighbors = next(iter(neighbors.values()))
    assert first_neighbors
    assert set(demo_adata.obs.loc[first_neighbors, "perturbation_status"]) == {"unassigned"}
