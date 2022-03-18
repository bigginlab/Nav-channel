import argparse
import pickle

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

parser = argparse.ArgumentParser(description="Plot results from distance data")
parser.add_argument(
    "-data_files",
    dest="data_files",
    nargs="+",
    help="data files (.pkl) to analyse",
)
parser.add_argument(
    "-save_name",
    dest="save_name",
    type=str,
    default="all_distances",
    help="the name for the returned plot and .pkl file",
)
parser.add_argument(
    "-kde",
    dest="kde",
    type=bool,
    default=False,
)
parser.add_argument(
    "-xlabel",
    dest="xlabel",
    nargs="+",
    default="Minimum distance",
    help="the xlabel",
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


def plot_hist(df, save_name: str, kde: bool=False, xlabel=None) -> None:

    if xlabel is None:
        xlabel = "Minimum distance"

    if kde:

        g = sns.kdeplot(
            data=df,
            x=r"Minimum distance ($\AA$)",
            hue="System",
            palette="colorblind",
            shade=True,
        )

    else:

        g = sns.displot(
            df,
            x=r"Minimum distance ($\AA$)",
            hue="System",
            palette="colorblind",
            stat="density",
        )

    plt.ylabel("Density", size=14)

    if xlabel is not None:
        plt.xlabel(f"{xlabel} " + r"($\AA$)", size=14)

    g.set_xticklabels(g.get_xticks(), size = 12)
    g.set_yticklabels(g.get_yticks(), size = 12)

    if kde:
        print(f"Saving plot as {save_name}_kde.png")
        plt.savefig(f"{save_name}_kde.png")

    else:
        print(f"Saving plot as {save_name}_hist.png")
        plt.savefig(f"{save_name}_hist.png")

    plt.clf()


if __name__ == "__main__":

    args = parser.parse_args()

    df_list = load_pickle_files(args.data_files)

    df_concat = pd.concat(df_list)

    xlabel = " ".join(args.xlabel)


    print("Creating plots...")
    plot_hist(df_concat, args.save_name, kde=args.kde, xlabel=xlabel)
