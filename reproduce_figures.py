"""Draw the 18 article and supplement PDFs from the archived research data."""

import argparse
from pathlib import Path
from subprocess import run
import sys


ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="Directory for the rendered figures")
    output = parser.parse_args().output.resolve()
    if output == ROOT/"figures":
        parser.error("Choose a new directory to preserve the archived figure PDFs.")
    output.mkdir(parents=True, exist_ok=True)
    commands = [
        ("figures/plot_geometry.py", output),
        ("figures/finite_element_band_validation.py", "--pdf",
         output/"finite_element_band_validation.pdf"),
        ("figures/plot_band_figures.py", output),
        ("figures/plot_control_figures.py", output),
        ("figures/plot_summary_figures.py", output),
        ("figures/plot_field_figures.py", output),
        ("bands/plot_coupling.py", output/"bandwidth_and_channel_attenuation.pdf"),
    ]
    for source, *arguments in commands:
        run([sys.executable, "-B", str(ROOT/source), *map(str, arguments)], check=True)


if __name__ == "__main__":
    main()
