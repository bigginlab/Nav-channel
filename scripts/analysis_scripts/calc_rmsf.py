import argparse
import os

import matplotlib.pyplot as plt
import MDAnalysis as mda
import seaborn as sns
from MDAnalysis.analysis import rms

import pandas as pd

import pickle

parser = argparse.ArgumentParser(description="Calculate the rmsf of a protein.")
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
    help="the residues to use in the rmsf calculation e.g. [1,2,3] or [1:3]",
)
parser.add_argument(
    "-plot_hist",
    dest="plot_hist",
    type=bool,
    default=False,
    help="specify whether to plot a histogram of rmsf",
)
parser.add_argument(
    "-save_name",
    dest="save_name",
    type=str,
    default="rmsf",
    help="the name for returned plots and .pkl files",
)
parser.add_argument(
    "-system_name",
    dest="system_name",
    type=str,
    default="WT",
    help="the name of the system: WT, D211K, T207S, or Double",
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


def calc_rmsf(
    universe: mda.Universe,
    residues: str,
    run_number: int,
    system_name: str,
) -> list:

    if residues is not None:
        rmsf_group = f"name CA and resid {residues}"

    else:
        rmsf_group = "name CA"

    sel = universe.select_atoms(rmsf_group)

    R = rms.RMSF(sel)

    R.run()

    df = pd.DataFrame(
        {"Residue" : sel.resids, "RMSF": R.rmsf}
        )


    df["Run"] = run_number
    df["System"] = system_name

    return df


def plot_line(rmsfs: list, save_name: str) -> None:

    for i in rmsfs:
        print(i)
    
    print("here")

    df_concat = pd.concat(rmsfs)

    fig, ax = plt.subplots()

    sns.lineplot(
        data=df_concat,
        x="Residue",
        y="RMSF",
        alpha=0.75,
        palette="colorblind",
    )

    plt.ylabel(f"Calpha RMSF ($\AA$)")
    plt.xlabel(f"Residue number")

    plt.tight_layout()

    plt.savefig(f"{save_name}_line.png")


if __name__ == "__main__":

    args = parser.parse_args()

    trajs = get_trajs(args.traj_path)

    rmsf_store = []

    print(f"System name: {args.system_name}...")

    for i, traj in enumerate(trajs):

        print(f"Analysing trajectory {i + 1} of {len(trajs)}...")

        u = load_universe(args.pdb_file, traj)

        rmsf_store.append(calc_rmsf(u, args.residues, i, args.system_name))

    print("Plotting...")
    plot_line(rmsf_store, args.save_name)

    # Save the results
    pickle_name = f"{args.save_name}.pkl"
    print(f"Saving results to {pickle_name}...")

    with open(pickle_name, "wb") as f:
        pickle.dump(pd.concat(rmsf_store), f)
