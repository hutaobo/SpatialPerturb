def test_import():
    import spatialperturb as sp
    assert hasattr(sp, "__version__")

def test_cli_import():
    from spatialperturb.cli import app
    assert app is not None
