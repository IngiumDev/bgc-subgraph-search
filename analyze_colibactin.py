import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from pathlib import Path
import json
import numpy as np


sns.set_theme(style="ticks", context="paper", palette="muted")
sns.set_context(
    "paper",
    rc={
        "font.size": 18,  # base font size
        "axes.titlesize": 20,  # title
        "axes.labelsize": 18,  # axis labels
        "xtick.labelsize": 16,  # x tick labels
        "ytick.labelsize": 16,  # y tick labels
        "legend.fontsize": 12,  # legend
        "lines.linewidth": 1.5,  # line width
        "lines.markersize": 10,  # marker size
    },
)

plt.rcParams["figure.figsize"] = (6.5, 4.5)  # width x height in inches
mpl.rcParams["legend.title_fontsize"] = 13


def compare_fg_to_gc():
    fg_data = pd.read_csv("./data/comparisons/fg_to_gc/fg.tsv", sep="\t", header=None, names=["gene", "amount", "number_domains"])
    fg_data["Type"] = "Exact Search"

    fg_reg_data = pd.read_csv("./data/comparisons/fg_to_gc/fg_reg.tsv", sep="\t", header=None, names=["gene", "amount", "number_domains"])
    fg_reg_data["Type"] = "Regex Search"
    
    fg_perm_data = pd.read_csv("./data/comparisons/fg_to_gc/fg_perm.tsv", sep="\t", header=None, names=["gene", "amount", "number_domains"])
    fg_perm_data["Type"] = "Permutation Search"

    df = pd.concat([fg_data, fg_reg_data, fg_perm_data])
    df = df.sort_values(by="gene")

    plt.figure(figsize=(12, 6))
    ax = sns.barplot(data=df, x="gene", y="amount", hue="Type", edgecolor="black", linewidth=1.5)
    ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
    ax.xaxis.grid(False)
    plt.xticks(rotation=90)
    ax.set_title("Number of Identified GC IDs per Search Strategy")
    ax.set_ylabel("Number of identified GC IDs")
    ax.set_xlabel("Query")
    ax.set_yscale("log")
    plt.tight_layout()
    plt.xticks(rotation=45)
    plt.savefig("./colibactin_plots/comp_search.png", dpi=300)


def plot_gc_gain_vs_domain_complexity():

    fg_exact = pd.read_csv(
        "./data/comparisons/fg_to_gc/fg.tsv",
        sep="\t",
        header=None,
        names=["gene", "amount", "number_domains"]
    )

    fg_regex = pd.read_csv(
        "./data/comparisons/fg_to_gc/fg_reg.tsv",
        sep="\t",
        header=None,
        names=["gene", "amount", "number_domains"]
    )

    fg_perm = pd.read_csv(
        "./data/comparisons/fg_to_gc/fg_perm.tsv",
        sep="\t",
        header=None,
        names=["gene", "amount", "number_domains"]
    )

    fg_exact = fg_exact.rename(columns={"amount": "exact_amount"})
    fg_regex = fg_regex.rename(columns={"amount": "regex_amount"})
    fg_perm = fg_perm.rename(columns={"amount": "perm_amount"})

    merged = (
        fg_exact[["gene", "number_domains", "exact_amount"]]
        .merge(
            fg_regex[["gene", "number_domains", "regex_amount"]],
            on=["gene", "number_domains"],
            how="inner"
        )
        .merge(
            fg_perm[["gene", "number_domains", "perm_amount"]],
            on=["gene", "number_domains"],
            how="inner"
        )
    )

    merged["Regex Gain"] = merged["regex_amount"] - merged["exact_amount"]
    merged["Permutation Gain"] = merged["perm_amount"] - merged["exact_amount"]

    gain_df = merged.melt(
        id_vars=["gene", "number_domains"],
        value_vars=["Regex Gain", "Permutation Gain"],
        var_name="Search Strategy",
        value_name="GC Gain"
    )

    agg = (
        gain_df
        .groupby(["number_domains", "Search Strategy"], as_index=False)
        .agg(
            sum_gain=("GC Gain", "sum"),
            mean_gain=("GC Gain", "mean"),
            n_queries=("GC Gain", "count")
        )
    )

    plt.figure(figsize=(8, 6))
    ax = sns.lineplot(
        data=agg,
        x="number_domains",
        y="sum_gain",
        hue="Search Strategy",
        marker="o"
    )

    ax.set_title("Gain in Identified GCs vs Domain Complexity")
    ax.set_xlabel("Number of Domains per Query")
    ax.set_ylabel("GC Gain vs Exact Search")
    ax.yaxis.grid(True, linestyle="--", linewidth=0.5)


    plt.tight_layout()
    plt.savefig("./colibactin_plots/sum_gain.png", dpi=300)
    plt.close()


# OLD FUNC
def analyze_graphs():

    exact = pd.read_csv("./data/comparisons/graphs/default.tsv", sep="\t")
    exact["SearchType"] = "Exact"

    regex = pd.read_csv("./data/comparisons/graphs/regex.tsv", sep="\t")
    regex["SearchType"] = "Regex"

    perm = pd.read_csv("./data/comparisons/graphs/perm.tsv", sep="\t")
    perm["SearchType"] = "Permutation"

    df = pd.concat([exact, regex, perm], ignore_index=True)

    for search_type, sub in df.groupby("SearchType"):
        plt.figure(figsize=(10, 5))
        ax = sns.barplot(
            data=sub,
            x="Graph",
            y="Numdomains",
            hue="Taxonomy",
            edgecolor="black"
        )

        ax.set_title(f"Number of Domains per Graph ({search_type} Search)")
        ax.set_xlabel("Graph")
        ax.set_ylabel("Number of Domains")
        ax.legend(
            title="Taxonomy",
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False
        )
        plt.tight_layout()
        plt.savefig(f"./colibactin_plots/{search_type}_numdomains_per_graph.png", dpi=300)
        plt.close()

    for search_type, sub in df.groupby("SearchType"):
        plt.figure(figsize=(10, 5))
        ax = sns.barplot(
            data=sub,
            x="Graph",
            y="Numuniqdomains",
            hue="Taxonomy",
            edgecolor="black"
        )

        ax.set_title(f"Unique Domains per Graph ({search_type} Search)")
        ax.set_xlabel("Graph")
        ax.set_ylabel("Number of Unique Domains")
        ax.legend(
            title="Taxonomy",
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False
        )
        plt.tight_layout()
        plt.savefig(f"./colibactin_plots/{search_type}_numuniqdomains_per_graph.png", dpi=300)
        plt.close()

    graph_counts = (
        df.groupby("SearchType")["Graph"]
        .nunique()
        .reset_index(name="NumGraphs")
    )

    plt.figure(figsize=(6, 5))
    ax = sns.barplot(
        data=graph_counts,
        x="SearchType",
        y="NumGraphs",
        edgecolor="black"
    )

    ax.set_title("Number of Resulting Graphs \n per Search Strategy")
    ax.set_xlabel("Search Strategy")
    ax.set_ylabel("Number of Graphs")
    ax.yaxis.grid(True, linestyle="--", linewidth=0.5)

    plt.tight_layout()
    plt.savefig("./colibactin_plots/graphs_per_search_type.png", dpi=300)
    plt.close()

def analyze_graphs_from_json_names(base_dir="./data/comparisons/graphs_json"):

    search_types = {
        "colibactin_ecoli_default": "Default",
        "colibactin_ecoli_regex": "Regex",
        "colibactin_ecoli_perm": "Permutation"
    }

    rows = []

    # -------- Collect per-graph FG hits --------
    for search_key, search_label in search_types.items():
        search_dir = Path(base_dir) / search_key / "meta"

        for json_file in search_dir.glob("*.json"):
            graph_name = json_file.stem
            species = graph_name.split("_", 1)[1]

            with open(json_file) as f:
                node_map = json.load(f)["label_info"]

            # Count each FG occurrence
            for fg in node_map.values():
                rows.append({
                    "SearchType": search_label,
                    "Species": species,
                    "FG": fg
                })

    df = pd.DataFrame(rows)

    # -------- Aggregate per species + search type + FG --------
    agg = (
        df.groupby(["SearchType", "Species", "FG"])
        .size()
        .reset_index(name="Count")
    )

    # Order species by total matches
    species_order = (
        agg.groupby("Species")["Count"]
        .sum()
        .sort_values(ascending=False)
        .index
        .tolist()
    )


    import numpy as np
    for search_type in agg["SearchType"].unique():
        subset = agg[agg["SearchType"] == search_type]

        pivot = (
            subset.pivot(index="Species", columns="FG", values="Count")
            .fillna(0)
            .reindex(species_order, fill_value=0)
        )

        fig, ax = plt.subplots(figsize=(12, 6))

        left = np.zeros(len(pivot))  # track stacking position

        for fg in pivot.columns:
            values = pivot[fg].values
            bars = ax.barh(pivot.index, values, left=left, edgecolor="black", label=fg)

            # ---- Add FG name text inside each segment ----
            for i, (bar, val) in enumerate(zip(bars, values)):
                if val > 0:  # only label non-zero segments
                    ax.text(
                        left[i] + val / 2,
                        bar.get_y() + bar.get_height() / 2,
                        fg,
                        ha='center',
                        va='center',
                        fontsize=7,
                        color='black'
                    )

            left += values  # update stacking position

        ax.set_title(f"Matched Query Composition per Species ({search_type})")
        ax.set_xlabel("Number of Matches")
        ax.set_ylabel("Species")
        ax.xaxis.grid(True, linestyle="--", linewidth=0.5)
        plt.tight_layout()
        plt.savefig(f"./colibactin_plots/stacked_queries_labeled_{search_type.lower()}.png", dpi=300)
        plt.close()


        df_unique = df.drop_duplicates(["SearchType", "Species", "FG"])

        agg_unique = (
            df_unique.groupby(["SearchType", "Species", "FG"])
            .size()
            .reset_index(name="Count")
        )
    for search_type in agg_unique["SearchType"].unique():
        subset = agg_unique[agg_unique["SearchType"] == search_type]

        pivot = (
            subset.pivot(index="Species", columns="FG", values="Count")
            .fillna(0)
            .reindex(species_order, fill_value=0)
        )

        fig, ax = plt.subplots(figsize=(12, 6))
        left = np.zeros(len(pivot))

        for fg in pivot.columns:
            values = pivot[fg].values
            bars = ax.barh(pivot.index, values, left=left, edgecolor="black", label=fg)

            for i, (bar, val) in enumerate(zip(bars, values)):
                if val > 0:
                    ax.text(
                        left[i] + val / 2,
                        bar.get_y() + bar.get_height() / 2,
                        fg,
                        ha='center',
                        va='center',
                        fontsize=7
                    )

            left += values

        ax.set_title(f"Unique Query Presence per Species ({search_type})")
        ax.set_xlabel("Number of Unique Queries")
        ax.set_ylabel("Species")
        ax.xaxis.grid(True, linestyle="--", linewidth=0.5)
        plt.tight_layout()
        plt.savefig(f"./colibactin_plots/stacked_queries_unique_{search_type.lower()}.png", dpi=300)
        plt.close()

def analyze_graphs_from_json(base_dir="./data/comparisons/graphs_json"):
    """
    Aggregate JSONs per species and search type,
    compute unique matched queries per species, and plot.
    """
    search_types = {"colibactin_ecoli_default": "Default", "colibactin_ecoli_regex": "Regex", "colibactin_ecoli_perm":"Permutation"}

    # store sets per graph first
    graph_fg_sets = []

    for search_type in search_types.keys():
        search_dir = Path(base_dir) / search_type / "meta"
        for json_file in search_dir.glob("*.json"):
            graph_name = json_file.stem  # e.g. 0_Megasphaera_sp000417505
            species = graph_name.split("_", 1)[1]  # Megasphaera_sp000417505

            with open(json_file) as f:
                node_map = json.load(f)

            node_map = node_map["label_info"]
            fg_set = set(node_map.values())
            graph_fg_sets.append({
                "SearchType": search_types[search_type],
                "Graph": graph_name,
                "Species": species,
                "FGSet": fg_set
            })

    # ---------- Aggregate per species ----------
    species_agg = []

    for search_type in search_types.values():
        # filter graph sets for this search type
        df_search = [g for g in graph_fg_sets if g["SearchType"] == search_type]

        # group by species
        species_to_union = {}
        for g in df_search:
            species_to_union.setdefault(g["Species"], set()).update(g["FGSet"])

        # store results
        for species, fg_union in species_to_union.items():
            species_agg.append({
                "SearchType": search_type,
                "Species": species,
                "NumMatchedQueries": len(fg_union)
            })

    agg_df = pd.DataFrame(species_agg)

    species_order = (
        agg_df
        .groupby("Species")["NumMatchedQueries"]
        .max()
        .sort_values(ascending=False)
        .index
        .tolist()
    )

    # ---------- Plot ----------
    plt.figure(figsize=(12, 6))
    ax = sns.barplot(
        data=agg_df,
        order=species_order,
        y="Species",
        x="NumMatchedQueries",
        hue="SearchType",
        edgecolor="black"
    )

    ax.set_title("Number of Matched Queries per Species by Search Strategy")
    ax.set_ylabel("Species")
    ax.set_xlabel("Number of Matched Queries")
    ax.legend(title="Search Type", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    ax.xaxis.grid(True, linestyle="--", linewidth=0.5)
    plt.tight_layout()
    plt.savefig("./colibactin_plots/uniq_matched_queries_per_species.png", dpi=300)
    plt.close()


def plot_abundance_boxplots(
    base_dir="./data/runs",
):

    search_types = {
        "Exact": "colibactin_ecoli_default",
        "Regex": "colibactin_ecoli_regex",
        "Permutation": "colibactin_ecoli_perm",
    }

    for search_label, search_dirname in search_types.items():
        records = []

        # load graph name and id
        meta_dir = Path(base_dir) / search_dirname / "meta"
        for json_file in meta_dir.glob("*.json"):
            graph_name = json_file.stem  # e.g. 0_Megasphaera_sp000417505
            species = graph_name.split("_", 1)[1]  # Megasphaera_sp000417505
            graph_id = int(graph_name.split("_", 1)[0])  # 0

            node_info = {}
            with open(json_file) as f:
                node_info = json.load(f)["node_info"]

            for node_id, sample_vals in node_info.items():
                records.append({
                    "Graph": graph_name,
                    "Node": node_id,
                    "Abundance": sample_vals["abundance_mean"],
                })

        df = pd.DataFrame(records)

        graph_order = sorted(
            df["Graph"].unique(),
            key=lambda x: int(x.split("_", 1)[0])
        )

        plt.figure(figsize=(12, 6))
        ax = sns.boxplot(
            data=df,
            y="Graph",
            x="Abundance",
            order=graph_order,
            linewidth=1.5,
            showfliers=False
        )

        ax.set_title(f"Abundance per Graph ({search_label} Search)")
        ax.set_ylabel("Graph")
        ax.set_xlabel("Mean RPKM per Node")
        ax.xaxis.grid(True, linestyle="--", linewidth=0.5)

        plt.tight_layout()
        plt.savefig(f"./colibactin_plots/{search_label}_abundance_boxplots.png", dpi=300)
        plt.close()


def plot_prevalence_boxplots(
    base_dir="./data/runs",
):

    search_types = {
        "Exact": "colibactin_ecoli_default",
        "Regex": "colibactin_ecoli_regex",
        "Permutation": "colibactin_ecoli_perm",
    }

    for search_label, search_dirname in search_types.items():
        records = []

        # load graph name and id
        meta_dir = Path(base_dir) / search_dirname / "meta"
        for json_file in meta_dir.glob("*.json"):
            graph_name = json_file.stem  # e.g. 0_Megasphaera_sp000417505
            species = graph_name.split("_", 1)[1]  # Megasphaera_sp000417505
            graph_id = int(graph_name.split("_", 1)[0])  # 0

            node_info = {}
            with open(json_file) as f:
                data = json.load(f)
                node_info = data["node_info"]

            for node_id, sample_vals in node_info.items():
                records.append({
                    "Graph": graph_name,
                    "Node": node_id,
                    "Prevalence": float(sample_vals["gene_family_info"]["prevalence"].split(" ")[1].replace("%", "").replace("(", "").replace(")", "")),
                })

        df = pd.DataFrame(records)

        graph_order = sorted(
            df["Graph"].unique(),
            key=lambda x: int(x.split("_", 1)[0])
        )

        plt.figure(figsize=(12, 6))
        ax = sns.boxplot(
            data=df,
            y="Graph",
            x="Prevalence",
            order=graph_order,
            linewidth=1.5,
            showfliers=False
        )

        ax.set_title(f"Prevalence per Graph ({search_label} Search)")
        ax.set_ylabel("Graph")
        ax.set_xlabel("Node Prevalence")
        ax.xaxis.grid(True, linestyle="--", linewidth=0.5)

        plt.tight_layout()
        plt.savefig(f"./colibactin_plots/{search_label}_prevalence_boxplots.png", dpi=300)
        plt.close()


def plot_abundance_boxplots_per_d_group(
    base_dir="./data/runs",
    meta_data="./data/CRC_YachidaS_2019.metadata.csv"
):

    # disease group colors
    disease_palette = {
        "Healthy": "#2166ac",
        "CRC": "#762a83",
        "MP": "#9970ab"
    }

    search_types = {
        "Exact": "colibactin_ecoli_default",
        "Regex": "colibactin_ecoli_regex",
        "Permutation": "colibactin_ecoli_perm",
    }

    sample_meta_data = pd.read_csv(meta_data, sep=";")
    sample_to_group = dict(
        zip(sample_meta_data["sample_id"], sample_meta_data["disease_group"])
    )

    for search_label, search_dirname in search_types.items():
        records = []

        meta_dir = Path(base_dir) / search_dirname / "meta"

        for json_file in meta_dir.glob("*.json"):
            graph_name = json_file.stem
            graph_id = int(graph_name.split("_", 1)[0])

            data = ""
            with open(json_file) as f:
                data = json.load(f)

            node_info = data["node_info"]
            node_labels = data["label_info"]

            for node_id, node_data in node_info.items():
                if node_id not in node_labels.keys():
                    continue

                abundance_annot = node_data.get("abundance_annotation", {})

                # group abundances by disease group
                group_vals = {}
                for sample, abundance in abundance_annot.items():
                    group = sample_to_group.get(sample)
                    if group is None:
                        continue
                    group_vals.setdefault(group, []).append(abundance)

                # mean per disease group
                for group, values in group_vals.items():
                    if not values:
                        continue
                    records.append({
                        "Graph": graph_name,
                        "Node": node_id,
                        "Group": group,
                        "Abundance": sum(values) / len(values),
                    })

        if not records:
            continue

        df = pd.DataFrame(records)

        # sort graphs numerically
        graph_order = sorted(
            df["Graph"].unique(),
            key=lambda x: int(x.split("_", 1)[0])
        )

        hue_order = ["Healthy", "MP", "CRC"]


        df["Abundance_sqrt"] = np.sqrt(df["Abundance"])
        # ---------- Plot ----------



        plt.figure(figsize=(12, 9))
        ax = sns.boxplot(
            data=df,
            y="Graph",
            x="Abundance_sqrt",
            hue="Group",
            order=graph_order,
            hue_order=hue_order,
            palette=disease_palette,
            linewidth=1.5,
            showfliers=False
        )

        ax.set_title(f"Node Abundance per Graph by Disease Group \n ({search_label} Search)")
        ax.set_xlabel("Square Root Mean RPKM per Node")
        ax.set_ylabel("Graph")



        # shading
        for i in range(len(graph_order)):
            if i % 2 == 0:
                ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.4, zorder=0)

        plt.xticks(rotation=0)

        ax.legend(
            title="Disease group",
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False
        )

        ax.xaxis.grid(True, linestyle="--", linewidth=0.5)

        plt.tight_layout()
        plt.savefig(
            f"./colibactin_plots/{search_label}_abundance_boxplots_by_disease_group.png",
            dpi=300
        )
        plt.close()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.collections import PatchCollection

def analyze_contig_break():
    from matplotlib.colors import ListedColormap, BoundaryNorm
# Read the BLAST data

    
    columns = ['gene', 'contig', 'pident', 'length', 'mismatch', 'gapopen', 
               'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    df = pd.read_csv("./data/clbGenes/hits.tsv", sep='\t', names=columns)

    gene_order = ['clbB', 'clbC', 'clbD', 'clbE', 'clbF', 'clbG', 'clbH', 'clbI', 
                  'clbJ', 'clbK', 'clbL', 'clbM', 'clbN', 'clbO', 'clbP', 'clbQ', 'clbS']

    contig_order = ['k141_0', 'k141_1', 'k141_2', 'k141_3']

    plt.figure(figsize=(16, 6))

    hit_matrix = df.groupby(['gene', 'contig']).size().unstack(fill_value=0)
    hit_matrix = hit_matrix.reindex(index=gene_order, columns=contig_order, fill_value=0)

    presence_matrix = (hit_matrix > 0).astype(int)

    ## ORI
    # ax = sns.heatmap(hit_matrix, annot=True, fmt='d', cmap='YlOrRd', 
    #             cbar_kws={'label': 'Number of BLAST hits'})
    # ax.set_title('Number of BLAST Hits per Gene-Contig Pair', fontsize=14, fontweight='bold')
    # ax.set_xlabel('Contig', fontsize=12)
    # ax.set_ylabel('Gene (operon order)', fontsize=12)
    # ORI



    grey = sns.color_palette("muted")[7]  # or use "light:#808080" for custom grey
    green = sns.color_palette("muted")[2]  # seaborn's green
    red = sns.color_palette("muted")[3]    # seaborn's red
    cmap = ListedColormap([grey, green, red])
    # cmap = ListedColormap(['lightgrey', 'green', 'red'])

# Define boundaries so each integer gets its own color
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)

    ax = sns.heatmap(
        hit_matrix,
        annot=True,
        fmt='d',
        cmap=cmap,
        norm=norm,
        cbar_kws={'label': 'Number of BLAST hits', 'ticks': [0, 1, 2]}
    )

    ax.set_title('Number of BLAST Hits per Gene-Contig Pair', fontsize=14, fontweight='bold')
    ax.set_xlabel('Contig', fontsize=12)
    ax.set_ylabel('Gene order', fontsize=12)

    plt.tight_layout()
    plt.savefig('./colibactin_plots/contig_break/gene_contig_heatmap.png', dpi=300, bbox_inches='tight')
    plt.show()

# ==============================================================================
# 2. COVERAGE MAP - GENE TILING ACROSS CONTIGS
# ==============================================================================
    fig, axes = plt.subplots(len(contig_order), 1, figsize=(20, 18), sharex=False)

    colors = plt.cm.tab20(np.linspace(0, 1, len(gene_order)))
    gene_colors = dict(zip(gene_order, colors))

    for idx, contig in enumerate(contig_order):
        ax = axes[idx]
        contig_data = df[df['contig'] == contig].copy()
        
        if len(contig_data) == 0:
            ax.text(0.5, 0.5, 'No alignments', ha='center', va='center', 
                    fontsize=12, transform=ax.transAxes)
            ax.set_xlim(0, 1)
        else:
            # Get contig length
            contig_length = max(contig_data['send'].max(), contig_data['sstart'].max())
            
            for _, row in contig_data.iterrows():
                start = min(row['sstart'], row['send'])
                end = max(row['sstart'], row['send'])
                gene = row['gene']
                pident = row['pident']
                
                # Color intensity based on identity
                alpha = 0.4 + 0.6 * (pident / 100.0)
                
                rect = Rectangle((start, gene_order.index(gene)), 
                               end - start, 0.8,
                               facecolor=gene_colors[gene], 
                               edgecolor='black', 
                               linewidth=1.5,
                               alpha=alpha)
                ax.add_patch(rect)
                
                # Add gene label
                ax.text((start + end) / 2, gene_order.index(gene) + 0.4, 
                       gene, ha='center', va='center', fontsize=8, fontweight='bold')
        
            ax.set_xlim(0, contig_length)
        
        ax.set_ylim(-0.5, len(gene_order) - 0.5)
        ax.set_yticks(range(len(gene_order)))
        ax.set_yticklabels(gene_order)
        ax.set_ylabel('Genes', fontsize=10)
        ax.set_title(f'{contig} (length: {contig_length:,} bp)', fontsize=12, fontweight='bold')
        ax.grid(axis='x', alpha=0.3)
        ax.set_xlabel('Contig position (bp)', fontsize=10)

    plt.tight_layout()
    plt.savefig('./colibactin_plots/contig_break/gene_coverage_tiling.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    # compare_fg_to_gc()

    # analyze_graphs_from_json(base_dir="./data/comparisons/graphs/")
    # analyze_graphs_from_json_names(base_dir="./data/comparisons/graphs/")
    # plot_abundance_boxplots(base_dir="./data/comparisons/graphs/")
    # plot_prevalence_boxplots(base_dir="./data/comparisons/graphs/")

    # plot_abundance_boxplots_per_d_group(base_dir="./data/comparisons/graphs/")
    analyze_contig_break()

