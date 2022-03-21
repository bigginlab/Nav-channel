import argparse
import os

import matplotlib.pyplot as plt
import MDAnalysis as mda
import seaborn as sns
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.dihedrals import Dihedral, Janin, Ramachandran

import pandas as pd

import pickle

import numpy as np

from typing import Union

parser = argparse.ArgumentParser(description="Calculate the Ramachandran angles of a protein region.")
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
    "-residues",
    dest="residues",
    type=str,
    help="the residue to use in the Ramachandran calculation",
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


def calc_rama(
    universe: mda.Universe,
    run_number: int,
    system_name: str,
    # discard_time: float,
    residues: str = None,
) -> list:

    r = universe.select_atoms(f"resid {residues}")

    # TODO check discard time?
    R = Ramachandran(r).run() # can use start = start FRAME

    a = R.angles.reshape(
                np.prod(R.angles.shape[:2]), 2)

    df = pd.DataFrame(a)
    df.columns = ["Phi", "Psi"]
    df["System"] = system_name
    df["Run"] = run_number

    return df

def plot_joint(df: list, residues: str, save_name: str) -> None:

    fig, ax = plt.subplots()

    print(df)

    g = sns.jointplot(
        data=df,
        x="Phi",
        y="Psi",
        hue="System",
        palette="colorblind",
        xlim=(-180, 180),
        ylim=(-180, 180),
        ax=ax
    )

    # for ax in g.ax_joint:
    #     ax.axvline(0, color='black', lw=0.75)
    #     ax.axhline(0, color='black', lw=0.75)

    print(f"Saving plot as {save_name}.png")
    plt.savefig(f"{save_name}_jp.png")

    plt.clf()


if __name__ == "__main__":

    args = parser.parse_args()

    trajs = get_trajs(args.traj_path)

    rama_store = []

    print(f"System name: {args.system_name}...")

    residues = args.residues

    # print(f"Discarding {discard_time} ns from start of trajectory...")

    for i, traj in enumerate(trajs):

        print(f"Analysing trajectory {i + 1} of {len(trajs)}...")

        u = load_universe(args.pdb_file, traj)

        rama_store.append(calc_rama(universe=u, run_number=i, system_name=args.system_name, residues=residues))


    print("Plotting...")
    print(rama_store)
    df_concat = pd.concat(rama_store)
    plot_joint(df_concat, residues, args.save_name)

    # Save the results
    pickle_name = f"{args.save_name}.pkl"
    print(f"Saving results to {pickle_name}...")

    with open(pickle_name, "wb") as f:
        pickle.dump(pd.concat(rama_store), f)