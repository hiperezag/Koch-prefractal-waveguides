"""Export the eight Figure 2 point tables from the saved TE spectra."""

import argparse
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPECTRA = ROOT / "validation" / "dispersion" / "source_spectra"
sys.path.insert(0, str(ROOT / "gaps"))
from iem_bandas import cargar_bandas


CASES = (
    ("externa/triangulo/FK0/bandas_TE.dat",
     "bandas_TE_fem_externa_triangulo_FK0_L1_h0.04.dat"),
    ("igual_llenado/ref_triangulo/bandas_TE_interna_triangulo_FK3.dat",
     "bandas_TE_fem_interna_triangulo_FK3_L1_h0.04.dat"),
    ("igual_llenado/ref_triangulo/bandas_TE_externa_rectangulo_vertical_FK3.dat",
     "bandas_TE_fem_externa_rectangulo_vertical_FK3_2.19_h0.04.dat"),
    ("externa/hexagono/FK0/bandas_TE.dat",
     "bandas_TE_fem_externa_hexagono_FK0_L1_h0.04.dat"),
)


def write_points(path, k, frequencies):
    values = np.column_stack((np.repeat(k, frequencies.shape[1]), frequencies.ravel()))
    np.savetxt(path, values[np.isfinite(values[:, 1])], fmt="%.8f")


def export_tables(output):
    tables = []
    for panel, (iem_path, fem_name) in enumerate(CASES, start=1):
        iem_k, iem_frequencies = cargar_bandas(SPECTRA / iem_path)
        fem_k, fem_frequencies = cargar_bandas(SPECTRA / "fem" / fem_name)
        if len(iem_k) != 100 or len(fem_k) != 100:
            raise ValueError(f"Panel {panel} requires 100 wavenumbers in each spectrum.")
        np.testing.assert_allclose(iem_k, fem_k, atol=1e-10)

        if panel == 4:
            # The field test identifies this nearly flat branch as an interior resonance.
            interior_frequency = 2.682
            frequency_tolerance = 0.003
            interior = np.all(
                np.abs(iem_frequencies - interior_frequency) < frequency_tolerance,
                axis=0,
            )
            if np.count_nonzero(interior) != 1:
                raise ValueError("Expected one interior resonance in the FK0 hexagon IEM spectrum.")
            if np.ptp(iem_frequencies[:, interior]) >= 0.0002:
                raise ValueError("The FK0 hexagon interior resonance varies by 0.0002 or more.")
            if np.any(np.abs(fem_frequencies - interior_frequency) < frequency_tolerance):
                raise ValueError("The FK0 hexagon FEM spectrum contains an unexpected root near 2.682.")
            iem_frequencies = iem_frequencies[:, ~interior]

        tables.append((f"iem_{panel}.dat", iem_k, iem_frequencies))
        tables.append((f"fem_{panel}.dat", fem_k, fem_frequencies))

    output.mkdir(parents=True, exist_ok=True)
    for name, k, frequencies in tables:
        write_points(output / name, k, frequencies)
    print(f"Exported {len(tables)} Figure 2 point tables to {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="Directory for the eight point tables")
    export_tables(parser.parse_args().output)
