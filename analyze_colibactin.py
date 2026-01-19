import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="ticks", context="paper", palette="muted")

sns.set_context(
    "paper",
    rc={
        "font.size": 15,  # base font size
        "axes.titlesize": 17,  # title
        "axes.labelsize": 15,  # axis labels
        "xtick.labelsize": 14,  # x tick labels
        "ytick.labelsize": 14,  # y tick labels
        "legend.fontsize": 12,  # legend
        "legend.titlesize": 15,  # legend
        "lines.linewidth": 1.5,  # line width
        "lines.markersize": 10,  # marker size
    },
)

plt.rcParams["figure.figsize"] = (6.5, 4.5)  # width x height in inches


def compare_fg_to_gc():
    fg_data = pd.read_csv("./data/comparisons/fg_to_gc/fg.tsv", sep="\t", header=None, names=["gene", "amount"])
    fg_data["Type"] = "Default"

    fg_reg_data = pd.read_csv("./data/comparisons/fg_to_gc/fg_reg.tsv", sep="\t", header=None, names=["gene", "amount"])
    fg_reg_data["Type"] = "Regex"

    df = pd.concat([fg_data, fg_reg_data])
    df = df.sort_values(by="gene")

    plt.figure(figsize=(12, 6))
    ax = sns.barplot(data=df, x="gene", y="amount", hue="Type", edgecolor="black", linewidth=1.5)
    ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
    ax.xaxis.grid(False)
    plt.xticks(rotation=90)
    ax.set_title("Gain in Identified GC IDs using Regex Seach")
    ax.set_ylabel("Number of identified GC IDs")
    ax.set_xlabel("Query")
    ax.set_yscale("log")
    plt.tight_layout()
    plt.savefig("./plots/comp_search.png", dpi=300)



if __name__ == "__main__":
    compare_fg_to_gc()
