"""Submodule to test visualization utilities."""

from typing import List
import pandas as pd
import matplotlib.pyplot as plt
from rdkit.Chem import Mol  # pylint: disable=no-name-in-module
from rdkit.Chem import MolFromSmiles  # pylint: disable=no-name-in-module
from map4 import MAP4


def test_visualization():
    """Test the visualization utilities."""
    dimensions: int = 2048

    smiles_df: pd.DataFrame = pd.read_csv("tests/smiles_np_classifier.csv.gz")

    map_calculator: MAP4 = MAP4(
        dimensions=dimensions,
    )

    molecules: List[Mol] = [MolFromSmiles(smiles) for smiles in smiles_df["smiles"]]
    labels: List[str] = list(smiles_df["pathway_label"])

    fig, ax = map_calculator.visualize(
        molecules=molecules,
        labels=labels,
        metric="cosine",
        preliminarly_reduce_with_pca=True,
    )
    ax.set_title("t-SNE of MAP4 (Molecular pathways)")
    fig.savefig("tests/tsne.png")
    plt.close(fig)
