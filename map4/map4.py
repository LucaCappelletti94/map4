"""Submodule for MAP4 fingerprint calculation."""

import itertools
from typing import List, Set, Dict, Iterable, Optional, Tuple
from collections import defaultdict, Counter
from multiprocessing import Pool

import numpy as np
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib import cm
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from mhfp.encoder import MHFPEncoder  # pylint: disable=no-name-in-module, import-error
from rdkit.Chem import MolToSmiles  # pylint: disable=no-name-in-module
from rdkit.Chem import PathToSubmol  # pylint: disable=no-name-in-module
from rdkit.Chem import Mol  # pylint: disable=no-name-in-module
from rdkit.Chem.rdmolops import GetDistanceMatrix  # pylint: disable=no-name-in-module
from rdkit.Chem.rdmolops import (  # pylint: disable=no-name-in-module
    FindAtomEnvironmentOfRadiusN,  # pylint: disable=no-name-in-module
)


class MAP4:
    """Class to calculate the MAP4 fingerprint of a molecule."""

    def __init__(
        self,
        dimensions: int = 2048,
        radius: int = 2,
        include_duplicated_shingles: bool = False,
        seed: int = 75434278,
    ):
        """Initializes the MAP4 calculator.

        Parameters
        ----------
        dimensions : int = 2048
            The number of dimensions for the fingerprint
        radius : int = 2
            The radius of atom environments to consider
        include_duplicated_shingles : bool = False
            If True, a counter is employed to make each shingle unique
            and therefore also consider duplicates in the fingerprint
        seed : int = 75434278
            The seed for the MinHash algorithm

        """
        self.dimensions: int = dimensions
        self.radius: int = radius
        self.include_duplicated_shingles: bool = include_duplicated_shingles
        self.encoder: MHFPEncoder = MHFPEncoder(dimensions, seed=seed)

    def calculate(self, mol: Mol) -> np.ndarray:
        """Calculates the atom pair minhashed fingerprint

        Parameters
        ----------
        mol: Mol
            The molecule to calculate the fingerprint for

        Returns
        -------
        np.ndarray
            The fingerprint for the molecule
        """
        atom_env_pairs: Set[str] = self._calculate(mol)
        return self._fold(atom_env_pairs)

    def calculate_many(
        self,
        mols: Iterable[Mol],
        number_of_threads: Optional[int] = None,
        verbose: bool = False,
    ) -> List[np.ndarray]:
        """Calculates the atom pair minhashed fingerprint.

        Parameters
        ----------
        mols: Iterable[Mol]
            The molecules to calculate the fingerprints for
        number_of_threads: Optional[int] = None
            The number of threads to use for the calculation.
            If None, the number of threads is set to the number of CPUs.

        Returns
        -------
        List[np.ndarray]
            The fingerprints for each molecule
        """
        with Pool(number_of_threads) as pool:
            fingerprints = list(
                tqdm(
                    pool.imap(self.calculate, mols),
                    total=len(mols),
                    desc="Calculating fingerprints",
                    disable=not verbose,
                )
            )
            pool.close()
            pool.join()

        return fingerprints

    def visualize(
        self,
        molecules: List[Mol],
        labels: List[str],
        number_of_threads: Optional[int] = None,
        metric: str = "approximated_jaccard",
        preliminarly_reduce_with_pca: bool = False,
        verbose: bool = True,
    ) -> Tuple[Figure, Axes]:
        """Visualizes the molecules using t-SNE.

        Parameters
        ----------
        molecules : List[Mol]
            The molecules to visualize
        labels : List[str]
            The labels for each molecule.
        number_of_threads : int = -1
            The number of threads to use for the calculation.
            If None, the number of threads is set to the number of CPUs.
        metric: str = "approximated_jaccard"
            The metric to use for the t-SNE calculation.
            If 'approximated_jaccard', the approximated Jaccard for MinHashed fingerprints is used.
            Other viable metrics include 'cosine'.
        preliminarly_reduce_with_pca : bool = False
            If True, the data is preliminarly reduced with PCA.
            This flag should NOT be employed if approximated Jaccard is requested, as it
            would no longer make sense.
        verbose : bool = True
            If True, the progress is shown.

        Returns
        -------
        Tuple[Figure, Axes]
            The figure and axes of the plot
        """
        if preliminarly_reduce_with_pca and metric == "approximated_jaccard":
            raise ValueError(
                "Preliminary reduction with PCA is not compatible with MinHashed fingerprints."
            )

        fingerprints: List[np.ndarray] = self.calculate_many(
            molecules,
            number_of_threads=number_of_threads,
            verbose=verbose,
        )

        fingerprints: np.ndarray = np.array(fingerprints)

        if preliminarly_reduce_with_pca:
            pca = PCA(n_components=50)
            fingerprints = pca.fit_transform(fingerprints)

        tsne = TSNE(
            n_components=2,
            metric=(
                self.encoder.distance if metric == "approximated_jaccard" else metric
            ),
            verbose=verbose,
            random_state=65467554,
            n_jobs=-1 if number_of_threads is None else number_of_threads,
        )

        fingerprints_embedded = tsne.fit_transform(fingerprints)

        # We determine the most common labels, as we may not be able to plot all of them
        # in different colors, and to avoid confusion we only plot the 19 most common ones,
        # with the 20-th color being used for all other labels.
        label_counter = Counter(labels)
        if len(label_counter) > 20:
            most_common_labels: List[str] = [
                label for label, _count in label_counter.most_common(19)
            ]
            most_common_labels.append("Other")
        else:
            most_common_labels: List[str] = [
                label for label, _count in label_counter.most_common()
            ]

        # We determine the color for each molecule using tab20.
        colormap = list(cm.get_cmap("tab20").colors)
        colors = [
            (
                colormap[most_common_labels.index(label)]
                if label in most_common_labels
                else colormap[19]
            )
            for label in labels
        ]

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10), dpi=200)

        ax.scatter(
            fingerprints_embedded[:, 0],
            fingerprints_embedded[:, 1],
            c=colors,
            marker=".",
            alpha=0.7,
        )

        ax.set_xlabel("t-SNE 1")
        ax.set_ylabel("t-SNE 2")

        # We prepare the legend
        legend_handles = [
            plt.Line2D(
                [0], [0], marker="o", color="w", markerfacecolor=color, markersize=10
            )
            for color in colormap
        ]

        ax.legend(
            legend_handles,
            most_common_labels + ["Other"],
            title="Labels",
            title_fontsize="large",
            fontsize="small",
            loc="upper right",
        )

        return fig, ax

    def _calculate(self, mol: Mol) -> Set[str]:
        """Returns the set of atom pairs in the molecule.

        Parameters
        ----------
        mol : Mol
            The molecule to calculate the fingerprint for.
        """
        return self._all_pairs(mol, self._get_atom_envs(mol))

    def _fold(self, pairs: Set[str]) -> np.ndarray:
        """Folds the fingerprint using the MinHash algorithm.

        Parameters
        ----------
        pairs : Set[str]
            The set of atom pairs in the molecule

        Returns
        -------
        np.ndarray
            The folded fingerprint
        """
        fp_hash = self.encoder.hash(pairs)
        return self.encoder.fold(fp_hash, self.dimensions)

    def _get_atom_envs(self, mol: Mol) -> Dict[int, List[Optional[str]]]:
        """Returns the atom environments for each atom in the molecule.

        Parameters
        ----------
        mol : Mol
            The molecule to calculate the fingerprint for

        Returns
        -------
        Dict[int, List[Optional[str]]]
            The atom environments for each atom in the molecule,
            containing
        """
        atoms_env: Dict[int, List[Optional[str]]] = {}
        for atom in mol.GetAtoms():
            atom_identifier: int = atom.GetIdx()
            for radius in range(1, self.radius + 1):
                if atom_identifier not in atoms_env:
                    atoms_env[atom_identifier] = []
                atoms_env[atom_identifier].append(
                    MAP4._find_env(mol, atom_identifier, radius)
                )
        return atoms_env

    @classmethod
    def _find_env(cls, mol: Mol, atom_identifier: int, radius: int) -> Optional[str]:
        """Returns a smile representation of the atom environment of a given radius.

        Parameters
        ----------
        mol : Mol
            The molecule to calculate the fingerprint for
        atom_identifier : int
            The index of the atom to calculate the environment for
        radius : int
            The radius of the environment

        Returns
        -------
        str
            The SMILES representation of the atom environment, or
            None if the atom is not found
        """
        atom_identifiers_within_radius: List[int] = FindAtomEnvironmentOfRadiusN(
            mol=mol, radius=radius, rootedAtAtom=atom_identifier
        )
        atom_map = {}

        sub_molecule: Mol = PathToSubmol(
            mol, atom_identifiers_within_radius, atomMap=atom_map
        )
        if atom_identifier in atom_map:
            smiles = MolToSmiles(
                sub_molecule,
                rootedAtAtom=atom_map[atom_identifier],
                canonical=True,
                isomericSmiles=False,
            )
            return smiles

        return None

    def _all_pairs(
        self, mol: Mol, atoms_env: Dict[int, List[Optional[str]]]
    ) -> Set[str]:
        """Returns the set of atom pairs in the molecule.

        Parameters
        ----------
        mol : Mol
            The molecule to calculate the fingerprint for
        atoms_env : Dict[int, List[Optional[str]]]
            The atom environments for each atom in the molecule

        Returns
        -------
        Set[str]
            The set of atom pairs in the molecule
        """
        atom_pairs: Set[str] = set()
        distance_matrix = GetDistanceMatrix(mol)
        num_atoms = mol.GetNumAtoms()
        shingle_dict = defaultdict(int)
        for idx1, idx2 in itertools.combinations(range(num_atoms), 2):
            dist = str(int(distance_matrix[idx1][idx2]))

            for i in range(self.radius):
                env_a: Optional[str] = atoms_env[idx1][i]
                env_b: Optional[str] = atoms_env[idx2][i]

                # None strings are treated as empty strings
                if env_a is None:
                    env_a = ""
                if env_b is None:
                    env_b = ""

                if len(env_a) > len(env_b):
                    larger_env: str = env_a
                    smaller_env: str = env_b
                else:
                    larger_env: str = env_b
                    smaller_env: str = env_a

                shingle: str = f"{smaller_env}|{dist}|{larger_env}"

                if self.include_duplicated_shingles:
                    shingle_dict[shingle] += 1
                    shingle += f"|{shingle_dict[shingle]}"

                atom_pairs.add(shingle.encode("utf-8"))
        return atom_pairs
