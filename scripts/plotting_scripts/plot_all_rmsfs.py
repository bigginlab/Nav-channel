import argparse
import pickle

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

parser = argparse.ArgumentParser(description="Plot all RMSF results.")
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
    default="all_RMSFs",
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


def plot_box(df, save_name: str) -> None:

    fig, ax = plt.subplots()

    from matplotlib.ticker import MultipleLocator

    ax.xaxis.set_minor_locator(MultipleLocator(5))

    g = sns.boxplot(
        data=df,
        x="Residue",
        y="RMSF",
        hue="System",
        palette="colorblind",
        ax=ax,
        width=.6,
        linewidth=0.5,
    )

    plt.legend(loc="upper right", fontsize=8)

    plt.xlabel("Residue number")
    plt.ylabel(f"Calpha RMSF ($\AA$)")

    print(f"Saving plot as {save_name}.png")
    plt.savefig(f"{save_name}_box.png")

    plt.clf()

if __name__ == "__main__":

    args = parser.parse_args()

    df_list = load_pickle_files(args.data_files)

    df_concat = pd.concat(df_list)

    print(df_concat)

    print("Creating plots...")
    plot_box(df_concat, args.save_name)

    # Save the results
    pickle_name = f"{args.save_name}.pkl"
    print(f"Saving results to {pickle_name}...")

    with open(pickle_name, "wb") as f:
        pickle.dump(df_concat, f)
