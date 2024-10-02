"""Submodule to test that the CLI command works as expected."""

import os
import subprocess
from typing import List
import numpy as np
from rdkit.Chem import MolFromSmiles  # pylint: disable=import-error, no-name-in-module
import pandas as pd
from map4 import MAP4


def test_cli():
    """Test that the CLI works as expected."""

    input_file = "tests/molecules.smi"
    output_file = "tests/fingerprints.csv"

    # Run the CLI command

    subprocess.run(
        [
            "map4",
            "-i",
            input_file,
            "-o",
            output_file,
            "-d",
            "1024",
            "-r",
            "2",
            "--clean-mols",
            "False",
            "--batch-size",
            "16",
        ],
        check=True,
    )

    # Check that the output file exists
    assert os.path.exists(output_file)

    df: pd.DataFrame = pd.read_csv(
        output_file, header=None, index_col=0
    )

    # The resulting dataframe will have the normalized SMILES as the
    # index and the fingerprint as the only 1024 columns
    assert df.shape == (64, 1024)

    # We compare the fingerprints obtained using the CLI with the ones
    # obtained using the Python API
    map4: MAP4 = MAP4(radius=2, dimensions=1024)

    # Load the molecules
    molecules = [
        MolFromSmiles(smile) for smile in pd.read_csv(input_file, header=None)[0]
    ]

    # Compute the fingerprints
    fingerprints: List[np.ndarray] = map4.calculate_many(molecules)

    for i, (smile, fingerprint) in enumerate(zip(df.index, df.values)):
        for left_bit, right_bit in zip(fingerprint, fingerprints[i]):
            assert left_bit == right_bit
