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
