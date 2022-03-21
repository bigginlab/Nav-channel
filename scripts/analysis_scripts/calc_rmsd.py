import argparse
import os

import matplotlib.pyplot as plt
import MDAnalysis as mda
import seaborn as sns
from MDAnalysis.analysis import rms

import pandas as pd

import pickle

from typing import Union

parser = argparse.ArgumentParser(description="Calculate the rmsd of a protein.")
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
    nargs="+",
    type=str,
    help="the residues to use in the RMSD calculation e.g. [1,2,3] or [1:3]",
)
parser.add_argument(
    "-plot_hist",
    dest="plot_hist",
    type=bool,
    default=False,
    help="specify whether to plot a histogram of rmsd",
)
parser.add_argument(
    "-save_name",
    dest="save_name",
    type=str,
    default="RMSD",
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


def make_sel_string(sel_list):
    sel_string = " ".join(sel_list)

    return sel_string


def calc_rmsd(
    universe: mda.Universe,
    run_number: int,
    system_name: str,
    residues: Union[str, None] = None,
) -> list:

    if residues is not None:

        rmsd_group = f"backbone and resid {residues}"

        R = rms.RMSD(
            universe,  # universe to align
            universe,  # reference universe or atomgroup
            select="backbone",  # group to superimpose and calculate RMSD
            groupselections=[rmsd_group],  # group for RMSD e.g. a protein region
            ref_frame=0,
        )  # frame index of the reference

        col_names = ["Frame", "Time (ns)", "Protein backbone", f"Residues {residues}"]

    else:
        R = rms.RMSD(
            universe,  # universe to align
            universe,  # reference universe or atomgroup
            select="backbone",  # group to superimpose and calculate RMSD
            ref_frame=0,
        )  # frame index of the reference

        col_names = ["Frame", "Time (ns)", f"Protein backbone RMSD ($\AA$)"]

    R.run()

    df = pd.DataFrame(
        R.rmsd,
        columns=col_names,
    )
    df["Run"] = run_number
    df["System"] = system_name

    # convert returned time to ns
    df["Time (ns)"] /= 1000

    return df


def plot_hist(rmsds: list, residues: str, save_name: str) -> None:

    df_concat = pd.concat(rmsds)

    fig, ax = plt.subplots()

    if residues is not None:

        df_melt = df_concat.melt(
            value_vars=[f"Protein backbone RMSD ($\AA$)", f"Residues {residues}"],
            var_name="Region",
            value_name="RMSD",
        )

        sns.displot(
            df_melt,
            x="RMSD",
            hue="Region",
            alpha=0.75,
            palette="colorblind",
            stat="density",
        )

    else:
        sns.displot(
            df_concat[f"Protein backbone RMSD ($\AA$)"],
            legend=False,
            palette="colorblind",
            stat="density",
        )

    plt.xlabel(f"Protein backbone RMSD ($\AA$)")

    plt.tight_layout()

    plt.savefig(f"{save_name}_hist.png")


def plot_line(rmsds: list, residues: str, save_name: str) -> None:

    df_concat = pd.concat(rmsds)

    fig, ax = plt.subplots()

    sns.lineplot(x="Time (ns)", y=f"Residues {residues}", data=df_concat)

    plt.tight_layout()

    plt.savefig(f"{save_name}_line.png")


if __name__ == "__main__":

    args = parser.parse_args()

    trajs = get_trajs(args.traj_path)

    rmsd_store = []

    print(f"System name: {args.system_name}...")

    if args.residues is not None:
        residues = make_sel_string(args.residues)

    for i, traj in enumerate(trajs):

        print(f"Analysing trajectory {i + 1} of {len(trajs)}...")

        u = load_universe(args.pdb_file, traj)

        if args.residues is not None:
            rmsd_store.append(
                calc_rmsd(
                    universe=u,
                    run_number=i,
                    system_name=args.system_name,
                    residues=residues,
                )
            )
        else:
            rmsd_store.append(
                calc_rmsd(universe=u, run_number=i, system_name=args.system_name)
            )

    if args.plot_hist:
        print("Plotting histogram...")
        plot_hist(rmsd_store, residues, args.save_name)

    print("Plotting timeseries...")
    plot_line(rmsd_store, residues, args.save_name)

    # Save the results
    pickle_name = f"{args.save_name}.pkl"
    print(f"Saving results to {pickle_name}...")

    with open(pickle_name, "wb") as f:
        pickle.dump(pd.concat(rmsd_store), f)
