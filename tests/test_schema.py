import numpy as np
import pandas as pd
from anndata import AnnData

import spatialperturb as sp


def test_ensure_spatialperturb_schema_fills_defaults():
    adata = AnnData(
        X=np.ones((3, 2)),
        obs=pd.DataFrame(index=["c1", "c2", "c3"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )

    sp.ensure_spatialperturb_schema(adata)
    sp.validate_spatialperturb_schema(adata)

    assert list(adata.obs.columns[:4]) == ["perturbation", "perturbation_status", "cell_type", "roi"]
    assert adata.obsm["spatial"].shape == (3, 2)
    assert "spatialperturb" in adata.uns
