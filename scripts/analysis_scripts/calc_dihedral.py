import argparse
import os

import matplotlib.pyplot as plt
import MDAnalysis as mda
import seaborn as sns
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.dihedrals import Dihedral, Janin

import pandas as pd

import pickle

from typing import Union

parser = argparse.ArgumentParser(description="Calculate the dihedral of a protein.")
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
    "-residue",
    dest="residue",
    type=str,
    help="the residue to use in the dihedral calculation",
)
parser.add_argument(
    "-angle_atoms",
    dest="angle_atoms",
    nargs="+",
    type=str,
    help="the atoms to use in the dihedral analysis",
)
parser.add_argument(
    "-plot_hist",
    dest="plot_hist",
    type=bool,
    default=False,
    help="specify whether to plot a histogram of dihedrals",
)
parser.add_argument(
    "-save_name",
    dest="save_name",
    type=str,
    default="Dihedral",
    help="the name for returned plots and .pkl files",
)
parser.add_argument(
    "-system_name",
    dest="system_name",
    type=str,
    default="WT",
    help="the name of the system: WT, D211K, T207S, or Double",
)
parser.add_argument(
    "-discard_time",
    dest="discard_time",
    type=float,
    default=100,
    help="the time, in ns, to discard",
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


def make_sel_string(sel_list):
    sel_string = " ".join(sel_list)

    return sel_string


def calc_dihedral(
    universe: mda.Universe,
    run_number: int,
    system_name: str,
    angle_atoms: str,
    discard_time: float,
    residue: Union[str, None] = None,
) -> list:

    sel = universe.select_atoms(f"resid {residue}")

    angle_dict = {'Time (ns)': [], 'Dihedral': []}

    sel_atoms = sel.select_atoms(f"name {angle_atoms}")

    for ts in u.trajectory:

        if ts.time / 1000 < discard_time:
            continue

        dihedral = sel_atoms.dihedral.value()

        angle_dict['Time (ns)'].append(ts.time / 1000)
        angle_dict['Dihedral'].append(dihedral)

    df = pd.DataFrame(angle_dict)

    df["Run"] = run_number
    df["System"] = system_name

    df["Time (ns)"] -= discard_time

    return df


def plot_hist(dihedrals: list, residue: str, save_name: str) -> None:

    df_concat = pd.concat(dihedrals)

    fig, ax = plt.subplots()

    sns.displot(
        df_concat["Dihedral"],
        legend=False,
        palette="colorblind",
        stat="density",
    )

    plt.xlabel(f"Dihedral angle")

    plt.tight_layout()

    plt.savefig(f"{save_name}_hist.png")


def plot_line(dihedrals: list, residue: str, save_name: str) -> None:

    df_concat = pd.concat(dihedrals)

    fig, ax = plt.subplots()

    sns.lineplot(x="Time (ns)", y=f"Dihedral", data=df_concat)

    plt.tight_layout()

    plt.savefig(f"{save_name}_line.png")


if __name__ == "__main__":

    args = parser.parse_args()

    trajs = get_trajs(args.traj_path)

    dihedral_store = []

    print(f"System name: {args.system_name}...")

    # if args.residues is not None:
    #     residue = make_sel_string(args.residues)

    residue = args.residue
    angle_atoms = " ".join(args.angle_atoms)
    discard_time = args.discard_time

    print(residue)
    print(angle_atoms)
    print(f"Discarding {discard_time} ns from start of trajectory...")

    for i, traj in enumerate(trajs):

        print(f"Analysing trajectory {i + 1} of {len(trajs)}...")

        u = load_universe(args.pdb_file, traj)

        dihedral_store.append(calc_dihedral(universe=u, run_number=i, system_name=args.system_name, residue=residue, discard_time=discard_time, angle_atoms=angle_atoms))

    if args.plot_hist:
        print("Plotting histogram...")
        plot_hist(dihedral_store, residue, args.save_name)

    print("Plotting timeseries...")
    plot_line(dihedral_store, residue, args.save_name)

    # Save the results
    pickle_name = f"{args.save_name}.pkl"
    print(f"Saving results to {pickle_name}...")

    with open(pickle_name, "wb") as f:
        pickle.dump(pd.concat(dihedral_store), f)
