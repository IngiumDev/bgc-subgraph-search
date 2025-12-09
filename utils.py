from tqdm import tqdm
import pandas as pd
from pymongo import MongoClient


from .single_cohort_API import get_FG2gc_id
from .parser import load_domain_search_result, parse_fmt6_blast_results

from typing import List, Dict, Tuple, FrozenSet, Optional, Union
from collections import defaultdict
import igraph as ig

import random
import re
import matplotlib.colors as mcolors

import numpy as np
import os
import sys

sys.path.append(os.path.abspath(os.path.dirname(__file__)))


# Run blast
def run_blast(
    query_fasta: str,
    blast_DB: str,
    threads: int,
    blast_result: str = None,
    blast_type: str = "blastn",
    e_value=0.01,
):
    if blast_type not in ["blastn", "blastp", "tblastn"]:
        print("not supported blast_type:", blast_type)
        return {}
    if blast_result == None:  # use the defult value if no output path is passed in
        blast_result = "{0}.{1}.out".format(query_fasta, blast_type)

    ## regarding evalue: default is 10
    os.system(
        f"{blast_type} -num_threads {threads} -db {blast_DB} -evalue {e_value} "
        f'-outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore" '
        f"-query {query_fasta} -out {blast_result}"
    )

    blast_restable = pd.read_csv(
        blast_result,
        sep="\t",
        names=[
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "qlen",
            "slen",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
        ],
    )
    blast_restable["cov"] = 100.0 * (
        blast_restable["length"] / blast_restable["qlen"]
    )  # make it 0~100%
    blast_restable.to_csv(blast_result, sep="\t", index=False, header=False)

    return


# moved to MAGraph_utils
def select_N_colors(N):
    random.seed(666)
    COLOR_LIST = list(mcolors.CSS4_COLORS)
    if N > len(COLOR_LIST):
        print("too many colors to select:", N, " max_N:", len(COLOR_LIST))
        print("warning: allowing duplicated colors")
    return random.sample(COLOR_LIST, N)


# count the elements in the list, return the counts
def count_lst(lst):
    ret = {}
    for k in lst:
        if k not in ret:
            ret[k] = 0
        ret[k] += 1
    return ret


# return {k:count}, k_lst, sorted by count and selected for the topN
def _topN_in_dic(count_dic, topN=1):
    if len(count_dic) == 0:
        return {}, []
    k_lst = list(count_dic.keys())
    count_lst = [count_dic[k] for k in k_lst]

    idx_lst = list(np.argsort(count_lst))
    idx_lst.reverse()
    topN = min(len(idx_lst), topN)

    sel_k = [k_lst[ii] for ii in idx_lst[:topN]]
    ret = {}
    for k in sel_k:
        ret[k] = count_dic[k]
    return ret, sel_k


class cohort_search:
    def __init__(self):
        self.prj_name = ""
        self.query_lst = []
        self.matched_gc_lst = []
        self.query2gc_id_dic = {}
        self.gc_id2query_dic = {}
        self.gc_id2color_dic = {}
        return


class cohort_domain_search(cohort_search):
    def __init__(self, domain_search_result_fp, prj_name, client):
        # the output file from api_cli.py helper convert_to_domain function
        self.prj_name = prj_name
        self.query2FG_dic = load_domain_search_result(
            domain_search_result_fp
        )  # {query: domain}
        self.FG_lst = list(set(self.query2FG_dic.values()))

        # color the node by matching to querys
        self.query_lst = list(self.query2FG_dic.keys())
        self.color_lst = select_N_colors(len(self.query_lst))
        self.query2color_dic = dict(zip(self.query_lst, self.color_lst))

        print(
            "Extracting matched gc_id for [{0}]\nsearch for {1} functional groups: \n\t{2}".format(
                prj_name, len(self.FG_lst), self.FG_lst
            )
        )
        self.FG2gc_id_dic = get_FG2gc_id(self.FG_lst, DB_name=prj_name, client=client)
        # now converting different id mappings
        self.matched_gc_lst = []

        self.query2gc_id_dic = {}  # generate a map {query: gc_id}
        self.gc_id2query_dic = {}  # a reverse map from gc_id to query
        for query in self.query2FG_dic:
            cur_FG = self.query2FG_dic[query]
            # handle when FG not present in the MGX
            cur_gc_lst = []
            if cur_FG in self.FG2gc_id_dic:
                cur_gc_lst = self.FG2gc_id_dic[cur_FG]
            self.matched_gc_lst = self.matched_gc_lst + cur_gc_lst
            N_gc = len(cur_gc_lst)
            print("{0}: matched to {1} gc_id".format(query, N_gc))

            if query not in self.query2gc_id_dic:
                self.query2gc_id_dic[query] = []
            self.query2gc_id_dic[query] = self.query2gc_id_dic[query] + cur_gc_lst
            for cur_gc in cur_gc_lst:
                self.gc_id2query_dic[cur_gc] = query

        # match gc_id to color
        self.gc_id2color_dic = {}
        for cur_gc in self.gc_id2query_dic:
            cur_q = self.gc_id2query_dic[cur_gc]
            self.gc_id2color_dic[cur_gc] = self.query2color_dic[cur_q]

        return

    def summary_overview(self):
        print(
            "{0} queries have assigned functional groups".format(len(self.query2FG_dic))
        )
        print(
            "{0} gc_id are selected from {1}:".format(
                len(self.gc_id2query_dic), self.prj_name
            )
        )

        for q in self.query2gc_id_dic:
            print(
                "\t {0}: matched to {1} gc_id".format(q, len(self.query2gc_id_dic[q]))
            )

        print("Most commonly used mapping dictionaries:")
        print("\t query2gc_id_dic: {query: gc_id_lst}")
        print("\t gc_id2query_dic: {gc_id: query}")
        print("\t gc_id2query_dic: {gc_id: query}")

        print("Other useful fields:")
        print("matched_gc_lst: a list of matched gc_id")

        #


class cohort_blast_search(cohort_search):
    def __init__(
        self,
        blast_search_result_fp: str,
        prj_name: str,
        client: MongoClient,
        eval_threshold: float = None,
        score_threshold: float = None,
        ident_threshold: float = None,
        cov_threshold: float = None,
    ):

        self.prj_name = prj_name
        self.blast_search_result_fp = blast_search_result_fp
        # the output file from api_cli.py helper run_blast function

        self.query2gc_id_dic = parse_fmt6_blast_results(
            blast_filepath=blast_search_result_fp,
            ident_threshold=ident_threshold,
            cov_threshold=cov_threshold,
            eval_threshold=eval_threshold,
            score_threshold=score_threshold,
        )  # Dict[str, List[str]] {query:[gc_id_lst]}

        # color the node by matching to querys
        self.query_lst = list(self.query2gc_id_dic.keys())
        self.query_lst.sort()
        self.color_lst = select_N_colors(len(self.query_lst))
        self.query2color_dic = dict(zip(self.query_lst, self.color_lst))

        print(
            "Extracting matched gc_id for [{0}]\nsearch for {1} queries: \n\t{2}".format(
                prj_name, len(self.query_lst), self.query_lst
            )
        )

        # now converting different id mappings
        self.matched_gc_lst = []
        self.gc_id2query_dic = {}  # a reverse map from gc_id to query
        for query in self.query2gc_id_dic:
            cur_gc_lst = self.query2gc_id_dic[query]
            self.matched_gc_lst = self.matched_gc_lst + cur_gc_lst
            N_gc = len(cur_gc_lst)
            print("{0}: matched to {1} gc_id".format(query, N_gc))
            for cur_gc in cur_gc_lst:
                self.gc_id2query_dic[cur_gc] = query

        # match gc_id to color
        self.gc_id2color_dic = {}
        for cur_gc in self.gc_id2query_dic:
            cur_q = self.gc_id2query_dic[cur_gc]
            self.gc_id2color_dic[cur_gc] = self.query2color_dic[cur_q]

        return

    def summary_overview(self):  # to fix here!
        print(
            "{0} queries have assigned matched to the catalog".format(
                len(self.query2gc_id_dic)
            )
        )
        print(
            "{0} gc_id are selected from {1}:".format(
                len(self.gc_id2query_dic), self.prj_name
            )
        )

        for q in self.query2gc_id_dic:
            print(
                "\t {0}: matched to {1} gc_id".format(q, len(self.query2gc_id_dic[q]))
            )

        print("Most commonly used mapping dictionaries:")
        print("\t query2gc_id_dic: {query: gc_id_lst}")
        print("\t gc_id2query_dic: {gc_id: query}")
        print("\t gc_id2query_dic: {gc_id: query}")

        print("Other useful fields:")
        print("matched_gc_lst: a list of matched gc_id")

        #


## SysBioMed 2025
# a. Given a list of Pfam domains, return how many unique proteins match each
# Pfam. Output should be a dictionary {pfam: N_matched_protein}
# ["PF06826", "PF05658"]
def count_prots_per_pfam(
    pfam_domains: list[str], db_host, db_port, user_name, user_pd
) -> dict[str, int]:
    """
    :arg1: pfam_domains → a list[str] of (mandatory) pfam domains
    """
    pfam_dict = defaultdict(set)

    # NOTE: we will probably pass client directly to this func, for now we hardcode it for simplicity
    client = MongoClient(
        host=db_host,
        port=db_port,
        username=user_name,
        password=user_pd,
        authSource="admin",
    )
    collection = client["CRC_YangY_2021"]["protein_pfam_domain"]

    # INFO: we need to watch out for proteins that have the same pfam more than once:
    # see: here one protein has 4 domains (but only two uniq)
    # {'_id': ObjectId('65524e0f04b7fa348fcfed15'), 'protein_id': 'cohort_merged__SRR16124168__k141_192152::2::261::1925::-', 'Length': 554, 'N_domains': 4, 'position': 1, 'start': 21, 'end': 183, 'Pfam_id': 'PF06826', 'Pfam_description': 'Predicted Permease Membrane Region', 'eval': 5e-43}
    # {'_id': ObjectId('65524e0f04b7fa348fcfed16'), 'protein_id': 'cohort_merged__SRR16124168__k141_192152::2::261::1925::-', 'Length': 554, 'N_domains': 4, 'position': 2, 'start': 212, 'end': 275, 'Pfam_id': 'PF02080', 'Pfam_description': 'TrkA-C domain', 'eval': 1.6e-07}
    # {'_id': ObjectId('65524e0f04b7fa348fcfed17'), 'protein_id': 'cohort_merged__SRR16124168__k141_192152::2::261::1925::-', 'Length': 554, 'N_domains': 4, 'position': 3, 'start': 303, 'end': 363, 'Pfam_id': 'PF02080', 'Pfam_description': 'TrkA-C domain', 'eval': 9.4e-12}
    # {'_id': ObjectId('65524e0f04b7fa348fcfed18'), 'protein_id': 'cohort_merged__SRR16124168__k141_192152::2::261::1925::-', 'Length': 554, 'N_domains': 4, 'position': 4, 'start': 379, 'end': 550, 'Pfam_id': 'PF06826', 'Pfam_description': 'Predicted Permease Membrane Region', 'eval': 2.3e-47}

    for p_id in pfam_domains:
        results = collection.find({"Pfam_id": p_id})
        print(results)

        for doc in results:
            pfam_dict[p_id].add(doc.protein_id)
            # "PF05658" → {prot1, prot2}

    res = {k: len(v) for k, v in pfam_dict.items()}  # convert to pfam: count
    return res


# b. Given an FG query (formatted as described in the slides), extract a list of mandatory domains.
#    If no mandatory domains are found, return an empty list and issue a warning. We will handle edge cases later.
def get_mandatory_domains(FG_query: str) -> list[str]:
    # we either make this relatively complex and work with the implied regex tree
    # or we use the knowledge, that any input regex of our format is less complex
    # - if the format is <token>, <next token>, ... and in our case
    #   one token == pfam + regex modifier, we can relatively easily check all cases
    # - it depends on how complex the input reges is allowed to get
    # - also, another problem are cases like (P1 | P2) since we are not allowed to use such cases as reducing pfams
    pass


# c. Given an FG query and a candidate FG, return True if the candidate FG matches the query; otherwise return False.
def matches_fg_query(FG_query: str, FG_candidate: str) -> bool:
    FG_query = FG_query.replace(", ", "").replace(" ", "")
    print(FG_query)
    pattern = re.compile(FG_query)
    match = pattern.search(FG_candidate)
    return match != None


if __name__ == "__main__":
    # examples
    FG_query = "PF03895, .*, (PF05658 | PF05662){1,20}, .*, PF03895, PF03895+, PF03895, (PF0023 | PF0055)"

    # a)
    db_host = sys.argv[1]
    db_port = sys.argv[2]
    user_name = sys.argv[3]
    user_pd = sys.argv[4]
    pfam_domains = ["PF03895", "PF05662"]
    count_prots_per_pfam(pfam_domains, db_host, db_port, user_name, user_pd)

    # b) TODO
    # get_mandatory_domains(
    #    FG_query
    # )

    # c)
    # FG_candidate = (
    #     "PF03895PF03895PF03895PF05658PF05662PF05662PF03895PF03895PF03895PF0023"
    # )
    # s = matches_fg_query(FG_query, FG_candidate)
    # print(s)
