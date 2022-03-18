import argparse
import pickle

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

parser = argparse.ArgumentParser(description="Plot all Ramachandran results.")
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
    help="the residues used for Ramachandran calculation",
)
parser.add_argument(
    "-save_name",
    dest="save_name",
    type=str,
    default="all_Rama",
    help="the name for returned plots and .pkl file",
)
parser.add_argument(
    "-wt_col",
    dest="wt_col",
    type=str,
    default="blue",
    help="the colour of the WT plot",
)
parser.add_argument(
    "-mut_col",
    dest="mut_col",
    type=str,
    default="orange",
    help="the colour of the mutant plot",
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


def plot(df, residues: str, save_name: str, wt_col: str = "blue", mut_col: str = "orange") -> None:

    pal = sns.color_palette('colorblind')
    pal_hex = pal.as_hex()

    cb_pal = {
        "blue": pal_hex[0],
        "orange": pal_hex[1],
        "green": pal_hex[2],
        "red": pal_hex[3],
        "purple": pal_hex[4],
        "brown": pal_hex[5],
        "pink": pal_hex[6],
        "grey": pal_hex[7],
        "gray": pal_hex[7],
        "yellow": pal_hex[8],
        "cyan": pal_hex[9],
    }

    custom_pal = [cb_pal[wt_col], cb_pal[mut_col]]

    # plot the 2D density, Ramachandran style
    g = sns.jointplot(
        data=df,
        x="Phi",
        y="Psi",
        hue="System",
        palette=custom_pal,
        xlim=(-180, 180),
        ylim=(-180, 180),
        kind="kde",
    )

    # add a scatter plot over the top
    g.plot_joint(sns.scatterplot, alpha=0.1)

    print(f"Saving plot as {save_name}.png")
    plt.savefig(f"{save_name}_jp.png")

    plt.clf()




if __name__ == "__main__":

    args = parser.parse_args()

    df_list = load_pickle_files(args.data_files)

    df_concat = pd.concat(df_list)

    print(df_concat)

    print("Creating plots...")
    if args.residues:
        residues = "Residue " + " ".join(args.residues)
    plot(df_concat, residues, args.save_name, wt_col=args.wt_col, mut_col=args.mut_col)

    # Save the results
    pickle_name = f"{args.save_name}.pkl"
    print(f"Saving results to {pickle_name}...")

    with open(pickle_name, "wb") as f:
        pickle.dump(df_concat, f)