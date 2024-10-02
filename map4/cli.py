"""Command line interface for map4."""

import argparse
from typing import List, Optional
import numpy as np
from rdkit.Chem import MolFromSmiles, GetMolFrags  # pylint: disable=no-name-in-module
from rdkit.Chem import MolToSmiles  # pylint: disable=no-name-in-module
from rdkit.Chem import Mol  # pylint: disable=no-name-in-module
from map4.map4 import MAP4


def to_smiles(mol: Mol) -> str:
    """Converts a molecule to the canonical SMILES string."""
    return MolToSmiles(mol, canonical=True, isomericSmiles=False)


def main():
    """Main function for the command line interface."""
    args: argparse.Namespace = parse_args()

    def _parse_line(line: str) -> Optional[Mol]:
        line: str = line.strip()
        fields: List[str] = line.split(args.delimiter)
        molecule: Mol = MolFromSmiles(fields[0])
        if molecule is None:
            return None

        if args.clean_mols:
            molecule: Mol = sorted(
                GetMolFrags(molecule, asMols=True),
                key=lambda mol: molecule.GetNumHeavyAtoms(),
                reverse=True,
            )[0]
            molecule: Mol = MolFromSmiles(to_smiles(molecule))
        return molecule

    calculator = MAP4(args.dimensions, args.radius, args.include_duplicated_shingles)

    def process(batch: List[str], output_file):
        maybe_molecules: List[Optional[Mol]] = [
            _parse_line(line) for line in batch
        ]
        # Remove None values
        molecules: List[Mol] = [
            molecule
            for molecule in maybe_molecules
            if molecule is not None
        ]
        fingerprints: List[np.ndarray] = calculator.calculate_many(molecules)
        for molecule, fingerprint in zip(molecules, fingerprints):
            if len(fingerprint) > 0:
                fp_str: str = args.fp_delimiter.join(str(v) for v in fingerprint)
                output_file.write(
                    args.delimiter.join([to_smiles(molecule), fp_str]) + "\n"
                )

    with open(args.input_path, "r", encoding="utf8") as input_file:
        with open(args.output_path, "w+", encoding="utf8") as output_file:
            batch: List[str] = []
            for line in input_file:
                batch.append(line)
                if len(batch) >= args.batch_size:
                    process(batch, output_file)
                    batch: List[str] = []
            process(batch, output_file)


def parse_args() -> argparse.Namespace:
    """Parses the command line arguments using argparse."""
    parser = argparse.ArgumentParser(description="MAP4 calculator")
    parser.add_argument("--input-path", "-i", help="", type=str, required=True)
    parser.add_argument("--output-path", "-o", help="", type=str, required=True)
    parser.add_argument(
        "--dimensions",
        "-d",
        help="Number of dimensions of the MinHashed fingerprint [DEFAULT: 2048]",
        type=int,
        default=2048,
        choices=[128, 512, 1024, 2048],
    )
    parser.add_argument(
        "--radius",
        "-r",
        help="Radius of the fingerprint [DEFAULT: 2]",
        type=int,
        default=2,
    )
    parser.add_argument(
        "--include-duplicated-shingles",
        help="Whether to include duplicated shingles in the fingerprint [DEFAULT: False]",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--clean-mols",
        help="Molecules will be canonicalized, cleaned, and chirality information will be removed, \
    NECESSARY FOR FINGERPRINT CONSISTENCY ACROSS DIFFERENT SMILES INPUT [DEFAULT: True].",
        type=lambda x: (str(x).lower() == "true"),
        default="True",
        metavar="True/False",
    )
    parser.add_argument(
        "--delimiter",
        help="Delimiter used for both the input and output files [DEFAULT: ',']",
        type=str,
        default=",",
    )
    parser.add_argument(
        "--fp-delimiter",
        help="Delimiter used between the numbers in the fingerprint output [DEFAULT: ',']",
        type=str,
        default=",",
    )
    parser.add_argument(
        "--batch-size",
        "-b",
        help="Numbers of molecules to process in a batch [DEFAULT: 500]",
        type=int,
        default=500,
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
