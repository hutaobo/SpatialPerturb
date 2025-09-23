import typer

app = typer.Typer(help="SpatialPerturb CLI")


@app.command()
def version():
    """Print package version."""
    from . import __version__
    typer.echo(__version__)


if __name__ == "__main__":
    app()
