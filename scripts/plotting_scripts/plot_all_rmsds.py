import argparse
import pickle

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

parser = argparse.ArgumentParser(description="Plot all RMSD results.")
parser.add_argument(
    "-data_files",
    dest="data_files",
    nargs="+",
    help="data files (.pkl) to analyse",
)
parser.add_argument(
    "-residues",
    dest="residues",
    nargs="+",
    default=None,
    help="the residues used for RMSD calculation in the format 1:10 or 1:10 12:20",
)
parser.add_argument(
    "-save_name",
    dest="save_name",
    type=str,
    default="all_RMSDs",
    help="the name for returned plots and .pkl file",
)
parser.add_argument(
    "-xlabel_name",
    dest="xlabel_name",
    type=str,
    default="xlabel",
    help="the custom xlabel if wanted",
)
parser.add_argument(
    "-kde",
    dest="kde",
    type=bool,
    default=False
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


def plot_hist(df, residues: str, save_name: str, xlabel_name, kde:bool = False) -> None:

    xmin, xmax = 0, 10

    if residues:

        xmax = df[residues].max() + 1

        if kde:
            g = sns.kdeplot(
                data=df,
                x=residues,
                hue="System",
                palette="colorblind",
                shade=True,
            )
        else:
            g = sns.displot(
                df,
                x=residues,
                hue="System",
                palette="colorblind",
                stat="density",
            )
    

    else:

        xmax = df["Protein backbone"].max() + 1

        if kde:
            g = sns.kdeplot(
                data=df,
                x="Protein backbone",
                hue="System",
                palette="colorblind",
                shade=True,
            )
        else:
            g = sns.displot(
                df,
                x="Protein backbone",
                hue="System",
                palette="colorblind",
                stat="density",
            )

    plt.ylabel("Density", size=14)

    if xlabel_name is not None:
        plt.xlabel(f"{xlabel_name} " + r"($\AA$)", size=14)

    g.set_xticklabels(g.get_xticks(), size = 12)
    g.set_yticklabels(g.get_yticks(), size = 12)


    if kde:
        plt.xlim(xmin, xmax) # hides -ve KDE section
        print(f"Saving plot as {save_name}_kde.png")
        plt.savefig(f"{save_name}_kde.png")

    else:
        
        print(f"Saving plot as {save_name}_hist.png")
        plt.savefig(f"{save_name}_hist.png")

    plt.clf()


def plot_line(df, residues: str, save_name: str, xlabel_name) -> None:

    if residues:

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

    print("Creating plots...")
    if args.residues:
        residues = "Residues " + " ".join(args.residues)
    else:
        residues = None
    plot_hist(df_concat, residues, args.save_name, args.xlabel_name, kde=args.kde)
    plot_line(df_concat, residues, args.save_name, args.xlabel_name)

    # Save the results
    pickle_name = f"{args.save_name}.pkl"
    print(f"Saving results to {pickle_name}...")

    with open(pickle_name, "wb") as f:
        pickle.dump(df_concat, f)
