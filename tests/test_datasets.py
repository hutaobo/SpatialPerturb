import gzip
import io
import tarfile

import pandas as pd
from scipy import io as spio
from scipy import sparse

import spatialperturb as sp


def _gzip_bytes(payload: bytes) -> bytes:
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb") as handle:
        handle.write(payload)
    return buffer.getvalue()


def _add_tar_bytes(archive: tarfile.TarFile, name: str, payload: bytes) -> None:
    info = tarfile.TarInfo(name=name)
    info.size = len(payload)
    archive.addfile(info, io.BytesIO(payload))


def test_prepare_gse241115_breast_cropseq_parses_sparse_and_control_inputs(tmp_path):
    raw_dir = tmp_path / "gse241115_breast_cropseq" / "raw"
    raw_dir.mkdir(parents=True)
    archive_path = raw_dir / "GSE241115_RAW.tar"

    control_counts = pd.DataFrame(
        {
            "CTRL_A.1": [5, 0],
            "CTRL_B.1": [0, 4],
        },
        index=["GeneA", "GeneB"],
    )
    sparse_counts = sparse.csr_matrix(
        [
            [5, 0, 1],
            [0, 4, 1],
            [1, 1, 0],
        ]
    )
    matrix_buffer = io.BytesIO()
    spio.mmwrite(matrix_buffer, sparse_counts)
    protospacer_calls = pd.DataFrame(
        [
            {"cell_barcode": "PERT_A-1", "num_features": 1, "feature_call": "ARID1A_sgRNA1", "num_umis": 10},
            {"cell_barcode": "PERT_B-1", "num_features": 1, "feature_call": "NT_sgRNA1", "num_umis": 8},
        ]
    )

    with tarfile.open(archive_path, "w") as archive:
        _add_tar_bytes(
            archive,
            "GSM7716949_counts_HCC38.csv.gz",
            _gzip_bytes(control_counts.to_csv().encode("utf-8")),
        )
        _add_tar_bytes(
            archive,
            "GSM7716951_matrix_HCC38_aggrMH001-3.mtx.gz",
            _gzip_bytes(matrix_buffer.getvalue()),
        )
        _add_tar_bytes(
            archive,
            "GSM7716951_barcodes_HCC38_aggrMH001-3.tsv.gz",
            _gzip_bytes(b"PERT_A-1\nPERT_B-1\nPERT_C-1\n"),
        )
        _add_tar_bytes(
            archive,
            "GSM7716951_features_HCC38_aggrMH001-3.tsv.gz",
            _gzip_bytes(
                b"ENSG000001\tGeneA\tGene Expression\n"
                b"ENSG000002\tGeneB\tGene Expression\n"
                b"ARID1A_sgRNA1\tARID1A_sgRNA1\tCRISPR Guide Capture\n"
            ),
        )
        _add_tar_bytes(
            archive,
            "GSM7716951_protospacer_calls_per_cell_HCC38_aggrMH001-3.csv.gz",
            _gzip_bytes(protospacer_calls.to_csv(index=False).encode("utf-8")),
        )

    prepare = sp.prepare_dataset("gse241115_breast_cropseq", cache_dir=tmp_path)
    adata = sp.load_public_dataset("gse241115_breast_cropseq", cache_dir=tmp_path)

    assert prepare["status"] == "ready"
    assert "guide_id" in adata.obs.columns
    assert "cell_line" in adata.obs.columns
    assert "HCC38" in set(adata.obs["cell_line"].astype(str))
    assert "ARID1A" in set(adata.obs["perturbation"].astype(str))
    assert "control" in set(adata.obs["perturbation"].astype(str))
    assert "unassigned" in set(adata.obs["perturbation"].astype(str))
    assert "ARID1A_sgRNA1" in adata.uns["spatialperturb"]["barcode_columns"]

    perturbed = adata.obs[adata.obs["guide_id"].astype(str) == "ARID1A_sgRNA1"].iloc[0]
    assert perturbed["perturbation"] == "ARID1A"
    assert perturbed["perturbation_status"] == "single"
