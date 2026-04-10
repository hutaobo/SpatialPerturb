import typer

app = typer.Typer(
    help="SpatialPerturb CLI",
    add_completion=False,
    no_args_is_help=True,
)


@app.callback()
def main() -> None:
    """Command group for SpatialPerturb utilities."""


@app.command()
def version() -> None:
    """Print package version."""
    from . import __version__

    typer.echo(__version__)


if __name__ == "__main__":
    app()
