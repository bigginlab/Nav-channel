import argparse
import pickle

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

parser = argparse.ArgumentParser(description="Plot all Dihedral results.")
parser.add_argument(
    "-data_files",
    dest="data_files",
    nargs="+",
    help="data files (.pkl) to analyse",
)
parser.add_argument(
    "-residue",
    dest="residue",
    nargs="+",
    help="the residue used for dihedral calculation",
)
parser.add_argument(
    "-xlabel",
    dest="xlabel",
    nargs="+",
    default="Dihedral angle",
    help="the xlabel",
)
parser.add_argument(
    "-save_name",
    dest="save_name",
    type=str,
    default="all_dihedrals",
    help="the name for returned plots and .pkl file",
)


def load_pickle_files(pickle_files):

    df_list = []

    for pf in pickle_files:
        try:
            with open(pf, "rb") as f:
                df_list.append(pickle.load(f))
        except:
            print(f"Warning, couldn't find {pf}!")
            continue

    return df_list


def plot_hist(df, residues: str, save_name: str, xlabel: str) -> None:

    g = sns.kdeplot(
        data=df,
        x="Dihedral",
        hue="System",
        palette="colorblind",
        shade=True,
        # stat="density",
    )

    plt.ylabel("Density", size=14)

    if xlabel is not None:
        # plt.xlabel(f"{xlabel}" + r" ($\theta$)", size=14)
        plt.xlabel("T/S207 dihedral angle" + r" ($\theta$)", size=14)

    g.set_xticklabels(g.get_xticks(), size = 12)
    g.set_yticklabels(g.get_yticks(), size = 12)


    print(f"Saving plot as {save_name}.png")
    plt.savefig(f"{save_name}_hist.png")

    plt.clf()


def plot_line(df, residues: str, save_name: str) -> None:

    print("hi")
    print(residues)

    if residues:
        print("here")
        print(residues)
        sns.lineplot(
            data=df,
            x="Time (ns)",
            y=residues,
            hue="System",
            palette="colorblind",
        )
    else:
        sns.lineplot(
            data=df,
            x="Time (ns)",
            y="Protein backbone",
            hue="System",
            palette="colorblind",
        )

    plt.savefig(f"{save_name}_line.png")

    plt.clf()


if __name__ == "__main__":

    args = parser.parse_args()

    df_list = load_pickle_files(args.data_files)

    df_concat = pd.concat(df_list)

    print(df_concat)

    print("Creating plots...")
    if args.residue:
        residues = "Residue " + " ".join(args.residue)
    plot_hist(df_concat, residues, args.save_name, args.xlabel)
    # plot_line(df_concat, residues, args.save_name)

    # Save the results
    pickle_name = f"{args.save_name}.pkl"
    print(f"Saving results to {pickle_name}...")

    with open(pickle_name, "wb") as f:
        pickle.dump(df_concat, f)
