"""Test suite for the map4 package."""

from typing import List
from rdkit.Chem import MolFromSmiles  # pylint: disable=no-name-in-module
from rdkit.Chem import Mol  # pylint: disable=no-name-in-module
import pandas as pd
from tqdm.auto import tqdm
from map4 import MAP4


def test_compute_sequential_fingerprints():
    """Test to evaluate the computation of MAP4 fingerprints"""
    dimensions: int = 1024

    smiles_df: pd.DataFrame = pd.read_csv(
        "tests/smiles_np_classifier.csv.gz", nrows=1000
    )

    for include_duplicated_shingles in tqdm((True, False), desc="Counted", leave=False):
        map_calculator: MAP4 = MAP4(
            dimensions=dimensions,
            include_duplicated_shingles=include_duplicated_shingles,
        )

        for smiles in tqdm(
            smiles_df.smiles,
            desc="Computing fingerprints",
            leave=False,
        ):
            molecule: Mol = MolFromSmiles(smiles)
            _map_fingerprint = map_calculator.calculate(molecule)


def test_compute_many_fingerprints():
    """Test to evaluate the computation of MAP4 fingerprints"""
    dimensions: int = 1024

    smiles_df: pd.DataFrame = pd.read_csv(
        "tests/smiles_np_classifier.csv.gz", nrows=1000
    )

    for include_duplicated_shingles in tqdm((True, False), desc="Counted", leave=False):
        map_calculator: MAP4 = MAP4(
            dimensions=dimensions,
            include_duplicated_shingles=include_duplicated_shingles,
        )

        molecules: List[Mol] = [
            MolFromSmiles(smiles)
            for smiles in tqdm(
                smiles_df.smiles,
                desc="Computing fingerprints",
                leave=False,
            )
        ]
        _map_fingerprints = map_calculator.calculate_many(molecules)
