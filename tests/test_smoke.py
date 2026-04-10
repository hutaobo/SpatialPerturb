import spatialperturb as sp
from typer.testing import CliRunner

from spatialperturb.cli import app


def test_import():
    assert hasattr(sp, "__version__")
    assert hasattr(sp, "build_signature_matrix")


def test_cli_version():
    runner = CliRunner()
    result = runner.invoke(app, ["version"])

    assert result.exit_code == 0
    assert sp.__version__ in result.stdout
