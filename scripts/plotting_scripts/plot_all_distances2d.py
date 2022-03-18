import argparse
from os import system
import pickle

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from typing import Union

parser = argparse.ArgumentParser(description="Plot results from distance data in 2D")
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
    "-xmax",
    dest="xmax",
    type=float,
    help="the maximum x value for the plot",
)
parser.add_argument(
    "-ymax",
    dest="ymax",
    type=float,
    help="the maximum y value for the plot",
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
    nargs="+",
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


def plot_joint(
    df,
    x: str,
    y: str,
    xmax: Union[float, None] = None,
    ymax: Union[float, None] = None,
    wt_col: str = "blue",
    mut_col: str = "orange",
    save_name: str = "jointplot",
) -> None:

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

    if xmax: 
        g = sns.jointplot(
            data=df,
            x=x,
            y=y,
            hue="System",
            palette=custom_pal,
            kind="kde",
            xlim=(0, xmax), 
        )
    if ymax: 
        g = sns.jointplot(
            data=df,
            x=x,
            y=y,
            hue="System",
            palette=custom_pal,
            kind="kde",
            ylim=(0, ymax), 
        )
    if xmax and ymax:
        g = sns.jointplot(
            data=df,
            x=x,
            y=y,
            hue="System",
            palette=custom_pal,
            kind="kde",
            xlim=(0, xmax),
            ylim=(0, ymax), 
        )
    else:
        g = sns.jointplot(
            data=df,
            x=x,
            y=y,
            hue="System",
            palette=custom_pal,
            kind="kde",
        )

    g.ax_joint.set_xlabel(fr"{x} minimum distance ($\AA$)", size=14)
    g.ax_joint.set_ylabel(fr"{y} minimum distance ($\AA$)", size=14)

    # g.set_xticklabels(g.get_xticks(), size = 12)
    # g.set_yticklabels(g.get_yticks(), size = 12)

    print(f"Saving plot as {save_name}_jp_v2.png")
    plt.savefig(f"{save_name}_jp_v2.png")

    plt.clf()


if __name__ == "__main__":

    args = parser.parse_args()

    df_list = load_pickle_files(args.data_files)

    new_df_dict = {}

    systems = [df_list[i]["System"].to_list()[0] for i in range(len(df_list))]

    for df in df_list:
        print(df.attrs["bond"])

        df.columns = [df.attrs["bond"], "System"]

        s = df["System"].to_list()[0]

        if s not in new_df_dict:
            new_df_dict[s] = []

        new_df_dict[s].append(df.drop(["System"], axis=1))

    new_df_concat_dict = {}

    for system_name in new_df_dict:

        new_df_concat_dict[system_name] = pd.concat(new_df_dict[system_name], axis=1)
        new_df_concat_dict[system_name]["System"] = system_name

    # Plotting
    print("Creating plot(s)")
    xmax = args.xmax
    ymax = args.ymax
    wt_col = args.wt_col
    mut_col_args = args.mut_col

    multiple_mut = False
    if len(mut_col_args) > 1:
        multiple_mut = True

    for i, system_name in enumerate(new_df_concat_dict):

        if system_name == "WT":
            continue

        df_to_plot = pd.concat(
            [new_df_concat_dict["WT"], new_df_concat_dict[system_name]]
        )

        # Get the residue distances to plot from df columns
        residue_distances = df_to_plot.columns.to_list()
        residue_distances.remove("System")

        print(residue_distances)

        

        # Get the distance label, form: 'X123 - Y456'
        x = residue_distances[0]
        y = residue_distances[1]

        # # Check if one residue in the pair has changed between axes, if so label appropriately
        # def _get_resid(str):
        #     alpha = ""
        #     num = ""
        #     for i in range(len(str)):
        #         if (str[i].isdigit()):
        #             num = num+ str[i]
        #         elif((str[i] >= 'A' and str[i] <= 'Z') or (str[i] >= 'a' and str[i] <= 'z')):
        #             alpha += str[i]
        #         else:
        #             continue

        #     return alpha, num

        # res0_x_resname, res0_x_resnum  = _get_resid(x.split(" - ")[0])
        # res1_x_resname, res1_x_resnum  = _get_resid(x.split(" - ")[1])

        # res0_y_resname, res0_y_resnum  = _get_resid(y.split(" - ")[0])
        # res1_y_resname, res1_y_resnum  = _get_resid(y.split(" - ")[1])

        # print(x, y)
        # print(res0_x_resname, res0_x_resnum, res1_x_resname, res1_x_resnum)

        # if res0_x_resname != res0_y_resname: # mutation

        #     new_res0_x = f"{res0_x_resname}/{res1_x_resname}{res0_x_resnum}"
        
        # if res0_y_resname != res1_y_resname: # mutation

        #     new_res0_y = f"{res0_y_resname}/{res1_y_resname}{res0_y_resnum}"






        # Check if the user has specified more than one colour
        # This implies multiple mutant systems have been provided
        if multiple_mut:
            mut_col = mut_col_args[i-1] # i-1 since we skip WT
        else:
            mut_col = mut_col_args

        # Check if alt names used and replace for saving later
        try:
            x_p = x.replace("/", "_")
            y_p = y.replace("/","_")
        except:
            continue
        
        x_pair = x_p.replace(" - ", "-")
        y_pair = y_p.replace(" - ", "-")


        plot_joint(
            df=df_to_plot,
            x=x,
            y=y,
            xmax=xmax,
            ymax=ymax,
            wt_col=wt_col,
            mut_col=mut_col,
            save_name=f"WT_{system_name}_{x_pair}_{y_pair}",
        )

