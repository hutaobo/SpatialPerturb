from typing import Iterable, Dict
import pandas as pd

def build_signature_matrix(gene_sets: Dict[str, Iterable[str]]) -> pd.DataFrame:
    """
    Build a binary signature-by-gene matrix from dict of gene sets.
    Parameters
    ----------
    gene_sets : dict
        Mapping signature name -> iterable of gene symbols.
    Returns
    -------
    pd.DataFrame
        DataFrame with signatures as rows and genes as columns (1/0).
    """
    genes = sorted({g for s in gene_sets.values() for g in s})
    df = pd.DataFrame(0, index=gene_sets.keys(), columns=genes, dtype=int)
    for name, gs in gene_sets.items():
        for g in gs:
            if g in df.columns:
                df.loc[name, g] = 1
    return df
