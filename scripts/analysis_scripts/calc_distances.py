import argparse
import os
import pickle

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd
import seaborn as sns
from MDAnalysis.analysis import distances

parser = argparse.ArgumentParser(
    description="Calculate distances between a pair of residues."
)
parser.add_argument(
    "-pdb_file",
    dest="pdb_file",
    type=str,
    help="the PDB file",
)
parser.add_argument(
    "-traj_path",
    dest="traj_path",
    type=str,
    help="the location of the trajectories",
)
parser.add_argument(
    "-residue_pair",
    dest="residue_pair",
    type=str,
    nargs="+",
    help="pair to analyse e.g. [1,2]",
)
parser.add_argument(
    "-plot_hist",
    dest="plot_hist",
    type=bool,
    default=False,
    help="specify whether to plot a histogram of distances",
)
parser.add_argument(
    "-analyse_lipid",
    dest="analyse_lipid",
    type=bool,
    default=False,
    help="specify whether to analyse lipid contacts",
)
parser.add_argument(
    "-cog",
    dest="cog",
    type=bool,
    default=True,
    help="specify whether to use the centre of geometry of terminal atoms in the distance calculation",
)
parser.add_argument(
    "-save_name",
    dest="save_name",
    type=str,
    default="distances",
    help="the name for the returned plot and .pkl file",
)
parser.add_argument(
    "-system_name",
    dest="system_name",
    type=str,
    default="system",
    help="the name of the system, e.g. WT",
)
parser.add_argument(
    "-discard_time",
    dest="discard_time",
    type=float,
    default=100,
)
parser.add_argument(
    "-alt_res0",
    dest="alt_res0",
    type=str,
)
parser.add_argument(
    "-alt_res1",
    dest="alt_res1",
    type=str,
)
parser.add_argument(
    "-alt_res_swap",
    dest="alt_res_swap",
    type=bool,
    default=False,
)

def get_trajs(traj_path: str) -> list:

    traj_list = []
    traj_list += [
        os.path.abspath(os.path.dirname(traj_path)) + "/" + f
        for f in os.listdir(traj_path)
        if f.endswith(".xtc")
    ]

    return traj_list


def load_universe(pdb_file: str, traj: str = None) -> mda.Universe:

    if traj is not None:
        u = mda.Universe(pdb_file, traj)
    else:
        u = mda.Universe(pdb_file)

    return u


def calc_min_dist(
    universe: mda.Universe,
    res1: int,
    res2: int,
    discard_time: float,
    analyse_lipid: bool = False,
    cog: bool = False,
) -> list:

    pair_min_distances = []

    # specify the atoms to use for centre of geometry calculation
    # TODO handle edge cases e.g. HID etc
    # TODO expand to more lipid groups
    # TODO support user defined atoms 
    cog_atoms_dict = {
        "ALA": "CB",
        "ARG": "NH1 NH2",
        "ASN": "ND2 OD1",
        "ASP": "OD1 OD2",
        "CYS": "SG",
        "GLN": "OE1 NE2",
        "GLU": "OE1 OE2",
        "GLY": "CA",
        "HIS": "CG CD1 CE1 NE2 CD2",
        "ILE": "CD CG1",
        "LEU": "CD CG1",
        "LYS": "NZ",
        "MET": "SD CE",
        "PHE": "CG CD1 CE1 CZ CE2 CD2",
        "PRO": "CB CD CG",
        "SER": "OG",
        # "SER": "HG",
        "THR": "CG2 OG1",
        # "THR": "HG1",
        "TRP": "CD1 NE1 CE2 CZ2 CH2 CZ3 CE3 CD2 CG",
        "TYR": "CG CD2 CE2 CZ OH CE1 CD1",
        "VAL": "CG1 CG2",
        "POPC": "P",  # lipid
    }

    res1_sel = universe.select_atoms(f"resid {res1}")

    if analyse_lipid:
        res2_sel = universe.select_atoms(f"resname {res2}")
    else:
        res2_sel = universe.select_atoms(f"resid {res2}")

    # Check if we want to use centre of geometry, see cog_atoms_dict
    if cog:

        res1_sel_atoms = res1_sel.select_atoms(
            f"resid {res1} and name {cog_atoms_dict[res1_sel.atoms.resnames[0]]}"
        )

        if analyse_lipid:
            res2_sel_atoms = res2_sel.select_atoms(
                f"resname {res2} and name {cog_atoms_dict[res2_sel.atoms.resnames[0]]}"
            )
        else:
            res2_sel_atoms = res2_sel.select_atoms(
                f"resid {res2} and name {cog_atoms_dict[res2_sel.atoms.resnames[0]]}"
            )

    else:
        res1_sel_atoms = res1_sel.atoms
        res2_sel_atoms = res2_sel.atoms

    for frame in universe.trajectory:

        # discard first N ns
        if frame.time / 1000 < discard_time:
            continue

        if cog:
            d = np.min(
                distances.distance_array(
                    res1_sel_atoms.atoms.center_of_geometry(),
                    res2_sel_atoms.atoms.positions,
                )
            )

        else:

            d = np.min(
                distances.distance_array(
                    res1_sel_atoms.positions, res2_sel_atoms.positions
                )
            )

        pair_min_distances.append(d)

    return pair_min_distances


def plot_pair_hist(df, residue_pair, save_name) -> None:

    # TODO take residue pair names in the form X000

    fig, ax = plt.subplots()

    sns.displot(
        df,
        x=r"Minimum distance ($\AA$)",
        alpha=0.75,
        palette="colorblind",
        stat="density",
    )

    plt.title(f"{residue_pair[0]} - {residue_pair[1]}")
    plt.xlabel(r"Minimum distance ($\AA$)")

    plt.tight_layout()

    plt.savefig(f"{save_name}.png")


def check_lipid_input(pdb_file: str, pair: list) -> None:
    assert (
        pair[0].isnumeric() and pair[1].isalpha()
    ), "for lipid analysis, check pair selection format: (resid, lipid resname)"

    u_ref = load_universe(pdb_file)
    lipid_resids = u_ref.select_atoms(f"resname {pair[1]}")

    assert len(lipid_resids) != 0, f"lipid resname {pair[1]} not recognised!"


def get_ol_res(pdb_file: str, pair: list, analyse_lipid: bool = False, alt_res0: str = None, alt_res1: str = None, alt_res_swap: bool = True) -> list:

    ref_u = mda.Universe(pdb_file)

    p1_tl = f"{ref_u.select_atoms(f'resid {pair[0]}').resnames[0]}"
    p1_ol = mda.lib.util.convert_aa_code(p1_tl)
    
    if alt_res0 is not None:
        if alt_res_swap:
            p1_str = f"{alt_res0}/{p1_ol}{pair[0]}"
        else:
            p1_str = f"{p1_ol}/{alt_res0}{pair[0]}"
    else:
        p1_str = f"{p1_ol}{pair[0]}"

    if analyse_lipid:
        p2_str = "POPC"
    else:
        p2_tl = f"{ref_u.select_atoms(f'resid {pair[1]}').resnames[0]}"
        p2_ol = mda.lib.util.convert_aa_code(p2_tl)
    
    if alt_res1 is not None:
        if alt_res_swap:
            p2_str = f"{alt_res1}/{p2_ol}{pair[1]}"
        else:
            p2_str = f"{p2_ol}/{alt_res1}{pair[1]}"
    else:
        p2_str = f"{p2_ol}{pair[1]}"

    ol_num_pair = [p1_str, p2_str]

    return ol_num_pair


if __name__ == "__main__":

    args = parser.parse_args()

    # Get trajectories from supplied directory
    traj_path = os.path.join(args.traj_path, "") # correct for if trailing "/" missed
    trajs = get_trajs(traj_path)

    pair = args.residue_pair

    discard_time = args.discard_time

    print(f"Will discard {discard_time} ns from trajectory...")

    if args.analyse_lipid:
        check_lipid_input(args.pdb_file, pair)

    dist_store = []

    for i, traj in enumerate(trajs):

        print(f"Analysing trajectory {i + 1} of {len(trajs)}...")

        u = load_universe(args.pdb_file, traj)

        # TODO get pair one letter AA codes
        dist_store.append(
            calc_min_dist(u, pair[0], pair[1], discard_time=discard_time, analyse_lipid=args.analyse_lipid, cog=args.cog)
        )

    dist_store_collapsed = [i for l in dist_store for i in l]

    # Get one-letter residue names, if lipid then this is just "POPC"
    print(args.alt_res0, args.alt_res1)
    ol_residues = get_ol_res(args.pdb_file, pair, args.analyse_lipid, alt_res0=args.alt_res0, alt_res1=args.alt_res1, alt_res_swap=args.alt_res_swap)

    df = pd.DataFrame(
        {
            r"Minimum distance ($\AA$)": dist_store_collapsed,
            "System": f"{args.system_name}",
        }
    )

    df.attrs["bond"] = f"{ol_residues[0]} - {ol_residues[1]}"

    print(df.attrs["bond"])

    if args.plot_hist:
        plot_pair_hist(df, ol_residues, save_name=args.save_name)

    # Save the results
    pickle_name = f"{args.save_name}.pkl"
    print(f"Saving results to {pickle_name}...")

    with open(pickle_name, "wb") as f:
        pickle.dump(df, f)
