import os, glob, sys
from collections import defaultdict

from pymongo import MongoClient
from pymongo.collection import Collection as pymongoCollection
from typing import List, Tuple, Dict, Optional
import warnings


MAX_BUFFER_SIZE = 5000  # split the input search list to avoid

### helper functions

"""
Search the values in search_lst against the search_key field in the target_collection
extract the values in the ret_key field
return {search_value:[target_val]}
"""


def _extract_mapping_from_DB(
    search_lst: list,
    search_key: str,
    target_collection: pymongoCollection,
    ret_key: str = "",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    """
    TODO: to check if search_key and ret_key valid
    """
    ret = {}

    cur_search_lst = []
    for ii, cur_ss in enumerate(search_lst):
        cur_search_lst.append(cur_ss)
        # when reaching the max_buffer_size, search the cur_search_lst and empty it afterwards
        if ii % max_buffer_size == 0 and len(cur_search_lst) > 0:
            DB_query_res = target_collection.find({search_key: {"$in": cur_search_lst}})
            for item in DB_query_res:
                if ret_key == "":
                    cur_res_val = item
                else:
                    cur_res_val = item[ret_key]
                cur_search_val = item[search_key]
                if cur_search_val not in ret:
                    ret[cur_search_val] = []
                ret[cur_search_val].append(cur_res_val)
            cur_search_lst = []
            continue

    # search the rest of search_lst when not empty
    if len(cur_search_lst) > 0:
        DB_query_res = target_collection.find({search_key: {"$in": cur_search_lst}})
        for item in DB_query_res:
            if ret_key == "":
                cur_res_val = item
            else:
                cur_res_val = item[ret_key]
            cur_search_val = item[search_key]
            if cur_search_val not in ret:
                ret[cur_search_val] = []
            ret[cur_search_val].append(cur_res_val)

    return ret


def gene_list_to_neighbours(gene_ids, client, prj_name, window_size=3):
    db = client[prj_name]
    # adjust for excluding
    window_size = window_size - 1

    pipeline = [
        # Match input gene documents
        {"$match": {"gene_id": {"$in": gene_ids}}},
        {
            "$lookup": {
                "from": "gene_location_info",
                "let": {
                    "t_contig": "$contig_id",
                    "t_loc": "$loc",
                    "t_gene": "$gene_id",
                },
                "pipeline": [
                    {
                        "$match": {
                            "$expr": {
                                "$and": [
                                    {"$eq": ["$contig_id", "$$t_contig"]},
                                    {
                                        "$gte": [
                                            "$loc",
                                            {"$subtract": ["$$t_loc", window_size]},
                                        ]
                                    },
                                    {
                                        "$lte": [
                                            "$loc",
                                            {"$add": ["$$t_loc", window_size]},
                                        ]
                                    },
                                ]
                            }
                        }
                    },
                    {
                        "$project": {
                            "_id": 1,
                            "gene_id": 1,
                            "loc": 1,
                            "start_pos": 1,
                            "end_pos": 1,
                            "strand": 1,
                        }
                    },
                ],
                "as": "neighbors",
            }
        },
        {"$project": {"gene_id": 1, "contig_id": 1, "loc": 1, "neighbors": 1}},
    ]
    return list(db.gene_location_info.aggregate(pipeline))


def _query_neighbours(gc_id_dict, window_size, client, prj_name):
    all_genes = []
    gene_to_gc_id = {}
    for gc_id in gc_id_dict:
        for gene in gc_id_dict[gc_id]:
            all_genes.append(gene)
            gene_to_gc_id[gene] = gc_id

    db_res = gene_list_to_neighbours(all_genes, client, prj_name, window_size)

    ret = defaultdict(set)
    for gene_info in db_res:
        gene_id = gene_info["gene_id"]
        neighbors = [n["gene_id"] for n in gene_info["neighbors"]]
        gc_id = gene_to_gc_id[gene_id]
        ret[gc_id].update(neighbors)

    return ret


### main functions
# return {taxa_id: linage}
def NCBI_taxa2lineage(
    taxa_id_lst: list,
    client: MongoClient,
    DB_name: str = "NCBI_database",
    col_name: str = "taxdump_info",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    taxa_id_lst = [int(it) for it in taxa_id_lst if str(it).isdigit()]
    ret = {}
    mycol = client[DB_name][col_name]
    pre_ret = _extract_mapping_from_DB(
        search_lst=taxa_id_lst,
        search_key="tax_id",
        ret_key="lineage",
        target_collection=mycol,
        max_buffer_size=max_buffer_size,
    )

    for k in pre_ret:
        ret[str(k)] = pre_ret[k][0]

    return ret


# return {acc_id: lineage}
# NCBI_refgenome_acc2lineage(acc_lst,client)
def NCBI_refgenome_acc2lineage(
    acc_lst: list,
    client: MongoClient,
    DB_name: str = "NCBI_refseq_bacteria_genomes",
    col_name: str = "access_id_taxonomy",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):

    ret = {}
    mycol = client[DB_name][col_name]
    pre_ret = _extract_mapping_from_DB(
        search_lst=acc_lst,
        search_key="access_id",
        ret_key="lineage",
        target_collection=mycol,
        max_buffer_size=max_buffer_size,
    )

    for k in pre_ret:
        ret[str(k)] = pre_ret[k][0]

    return ret


# return {gene_id: acc_id}
# NCBI_gene_id2acc(gene_id_lst,client)
def NCBI_gene_id2acc(
    NCBI_gene_id_lst: list,
    client: MongoClient,
    DB_name: str = "NCBI_refseq_bacteria_genomes",
    col_name: str = "gene_location_info",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    ret = {}
    mycol = client[DB_name][col_name]
    pre_ret = _extract_mapping_from_DB(
        search_lst=NCBI_gene_id_lst,
        search_key="gene_id",
        ret_key="genome_id",
        target_collection=mycol,
        max_buffer_size=max_buffer_size,
    )

    for k in pre_ret:
        ret[str(k)] = pre_ret[k][0]
    return ret


# return {contig_id: acc_id}
# NCBI_contig_id2acc(contig_id_lst,client)
def NCBI_contig_id2acc(
    NCBI_contig_id_lst: list,
    client: MongoClient,
    DB_name: str = "NCBI_refseq_bacteria_genomes",
    col_name: str = "gene_location_info",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    ret = {}
    mycol = client[DB_name][col_name]
    pre_ret = _extract_mapping_from_DB(
        search_lst=NCBI_contig_id_lst,
        search_key="contig_id",
        ret_key="genome_id",
        target_collection=mycol,
        max_buffer_size=max_buffer_size,
    )

    for k in pre_ret:
        ret[str(k)] = pre_ret[k][
            0
        ]  # only keep one accession id # assuming it is consistant
    return ret


# return {NCBI_gc_id: {gene_id:[acc_id, linage]}}
# NCBI_gc_id2taxa(gc_id_lst,client)
def NCBI_gc_id2taxa(
    NCBI_gc_id_lst: list,
    client: MongoClient,
    DB_name: str = "NCBI_refseq_bacteria_genomes",
    col_name: str = "access_id_taxonomy",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    # return {gc_id:[gene_id_lst]}
    gc_id2gene_id_dic = _gc_id2gene_id(NCBI_gc_id_lst, client, DB_name)
    gene_id_lst = []
    for gc_id in gc_id2gene_id_dic:
        gene_id_lst = gene_id_lst + gc_id2gene_id_dic[gc_id]
    gene_id_lst = list(set(gene_id_lst))

    # get gene2acc
    gene2acc_dic = NCBI_gene_id2acc(gene_id_lst, client)
    acc_lst = list(set(gene2acc_dic.values()))

    # get acc2lineage
    acc2lineage_dic = NCBI_refgenome_acc2lineage(acc_lst, client)

    # put everything up together
    ret = {}
    for gc_id in gc_id2gene_id_dic:
        if gc_id not in ret:
            ret[gc_id] = {}
        for cur_g in gc_id2gene_id_dic[gc_id]:
            cur_acc = gene2acc_dic[cur_g]
            cur_lineage = acc2lineage_dic[cur_acc]
            ret[gc_id][cur_g] = [cur_acc, cur_lineage]

    return ret


# return {NCBI_contig_id: [acc_id, linage]}
# NCBI_contig_id2taxa(contig_id_lst,client)
def NCBI_contig_id2taxa(
    NCBI_contig_id_lst: list,
    client: MongoClient,
    DB_name: str = "NCBI_refseq_bacteria_genomes",
    col_name: str = "access_id_taxonomy",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    # return {contig_id:a cc_id}
    contig_id2acc_dic = NCBI_contig_id2acc(NCBI_contig_id_lst, client)
    acc_lst = list(set(contig_id2acc_dic.values()))

    # get acc2lineage
    acc2lineage_dic = NCBI_refgenome_acc2lineage(acc_lst, client)

    # put everything up together
    ret = {}
    for contig_id in contig_id2acc_dic:
        cur_acc = contig_id2acc_dic[contig_id]
        cur_lineage = acc2lineage_dic[cur_acc]
        ret[contig_id] = [cur_acc, cur_lineage]

    return ret


# A.1.a
# return {gene_id:seq}
def _gene_id2seq(
    gene_id_lst: list,
    client: MongoClient,
    DB_name: str,
    col_name: str = "gene_sequences",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    ret = {}
    mycol = client[DB_name][col_name]
    pre_ret = _extract_mapping_from_DB(
        search_lst=gene_id_lst,
        search_key="gene_id",
        ret_key="seq",
        target_collection=mycol,
        max_buffer_size=max_buffer_size,
    )

    for k in pre_ret:
        ret[k] = pre_ret[k][0]

    return ret


# return {gene_id:gc_id}
def _gene_id2gc_id(
    gene_id_lst: list,
    client: MongoClient,
    DB_name: str,
    col_name: str = "nucleotide_clstr",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    mycol = client[DB_name][col_name]
    ret = {}

    pre_ret = _extract_mapping_from_DB(
        search_lst=gene_id_lst,
        search_key="gene_id",
        ret_key="representative",
        target_collection=mycol,
        max_buffer_size=max_buffer_size,
    )

    for k in pre_ret:
        ret[k] = pre_ret[k][0]

    return ret


# return {gc_id:[gene_id_lst]}
def _gc_id2gene_id(
    gc_id_lst: list,
    client: MongoClient,
    DB_name: str,
    col_name: str = "nucleotide_clstr",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    mycol = client[DB_name][col_name]
    ret = _extract_mapping_from_DB(
        search_lst=gc_id_lst,
        search_key="representative",
        ret_key="gene_id",
        target_collection=mycol,
        max_buffer_size=max_buffer_size,
    )

    return ret


# return {gc_id:{sample_id:RPKM}}
def _gc_id2RPKM(
    gc_id_lst: list,
    DB_name: str,
    client: MongoClient,
    col_name: str = "gc_mgx_RPKM",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    mycol = client[DB_name][col_name]
    ret_pre = _extract_mapping_from_DB(
        search_lst=gc_id_lst,
        search_key="gc_id",
        ret_key="RPKM",
        target_collection=mycol,
        max_buffer_size=max_buffer_size,
    )

    ret = {}
    for gc_id in ret_pre:
        ret[gc_id] = {}
        # confirming that only one set of RPKM records can be extracted for the abundance
        if len(ret_pre[gc_id]) != 1:
            if len(ret_pre[gc_id]) > 1:
                message = "Multiple abundance records! gc_id: {0} | DB: {1} | Collection: {2}".format(
                    gc_id, DB_name, col_name
                )
            else:
                message = "Missing abundance record! gc_id: {0} | DB: {1} | Collection: {2}".format(
                    gc_id, DB_name, col_name
                )
            warnings.warn(message)
        # extract the abundance record for current gc_id
        cur_abd_record = ret_pre[gc_id][0]
        for sample_id in cur_abd_record:
            ret[gc_id][sample_id] = float(cur_abd_record[sample_id])
    return ret


# return {gc_id:prot_id}
def _gc_id2prot_id(
    gc_id_lst: list,
    DB_name: str,
    client: MongoClient,
    col_name: str = "prot90_clstr",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    mycol = client[DB_name][col_name]

    ret = {}
    pre_ret = _extract_mapping_from_DB(
        search_lst=gc_id_lst,
        search_key="gc_id",
        ret_key="representative",
        target_collection=mycol,
        max_buffer_size=max_buffer_size,
    )

    for k in pre_ret:
        ret[k] = pre_ret[k][0]

    return ret


# return {gc_id:[msp_id,category]}
def _gc_id2msp_id(
    gc_id_lst: list,
    client: MongoClient,
    DB_name: str,
    col_name: str = "gc_to_msp",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    mycol = client[DB_name][col_name]

    pre_ret = _extract_mapping_from_DB(
        search_lst=gc_id_lst,
        search_key="gc_id",
        ret_key="",
        target_collection=mycol,
        max_buffer_size=max_buffer_size,
    )

    ret = {}
    for k in pre_ret:
        item = pre_ret[k][0]
        ret[k] = [item["msp_id"], item["gene_category"]]

    return ret


# return {prot_id:[gc_id_lst]}
def _prot_id2gc_id(
    prot_id_lst: list,
    DB_name: str,
    client: MongoClient,
    col_name: str = "prot90_clstr",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    mycol = client[DB_name][col_name]

    ret = _extract_mapping_from_DB(
        search_lst=prot_id_lst,
        search_key="representative",
        ret_key="gc_id",
        target_collection=mycol,
        max_buffer_size=max_buffer_size,
    )

    return ret


# return {gene_id:{contig_id,loc}}
def _gene_id2contig_loc(
    gene_id_lst: list,
    client: MongoClient,
    DB_name: str,
    col_name: str = "gene_location_info",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    mycol = client[DB_name][col_name]
    pre_ret = _extract_mapping_from_DB(
        search_lst=gene_id_lst,
        search_key="gene_id",
        ret_key="",
        target_collection=mycol,
        max_buffer_size=max_buffer_size,
    )

    ret = {}
    for gene_id in pre_ret:
        ret[gene_id] = {
            "contig_id": pre_ret[gene_id][0]["contig_id"],
            "loc": pre_ret[gene_id][0]["loc"],
        }
    return ret


# input: center_loc_dic {contig_id: [center_loc_lst]}


# return {contig_id:
# 			{
# 				gene_id:{"loc","start_pos","end_pos","strand"}
# 			}
# 	}
def _contig_id2gene_ids(
    contig_id_lst: list,
    client: MongoClient,
    DB_name: str,
    center_loc_dic: Dict,
    loc_window: int,
    col_name: str = "gene_location_info",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    # center_loc_dic {contig_id: [center_loc_lst]}

    # check if cur_loc belongs to the center_loc_lst within the loc_window
    def _check_loc(center_loc_lst, loc_window, cur_loc):
        for center_loc in center_loc_lst:
            if abs(center_loc - cur_loc) < loc_window:
                return True
        return False

    mycol = client[DB_name][col_name]
    # DB_query_res = mycol.find({"contig_id": {"$in": contig_id_lst}})
    pre_ret = _extract_mapping_from_DB(
        search_lst=contig_id_lst,
        search_key="contig_id",
        ret_key="",
        target_collection=mycol,
        max_buffer_size=max_buffer_size,
    )

    ret = {}
    for cur_contig_id in pre_ret:
        for item in pre_ret[cur_contig_id]:

            # check if cur_loc fit in the loc_window
            cur_loc = item["loc"]
            center_loc_lst = center_loc_dic[cur_contig_id]
            if not _check_loc(center_loc_lst, loc_window, cur_loc):
                continue
            cur_gene_info = {
                item["gene_id"]: {
                    "loc": item["loc"],
                    "start_pos": item["start_pos"],
                    "end_pos": item["end_pos"],
                    "strand": item["strand"],
                }
            }

            if cur_contig_id not in ret:
                ret[cur_contig_id] = {}

            ret[cur_contig_id][item["gene_id"]] = cur_gene_info

    return ret


# return a list of gene_id which within the window size of any gene_id in the input gene_id_lst
def _get_gene_neighbours_directly_connected(
    prj_name: str, gene_id_lst: list, neighbourhood_size: int, client: MongoClient
):
    # get all contiges
    # return {gene_id:{contig_id,loc}}
    gene_id2contig_loc_dic = _gene_id2contig_loc(gene_id_lst, client, prj_name)
    contig_id_lst = list(
        set(
            [
                gene_id2contig_loc_dic[gene_id]["contig_id"]
                for gene_id in gene_id2contig_loc_dic
            ]
        )
    )

    center_loc_dic = {}
    for gene_id in gene_id2contig_loc_dic:
        cur_contig_id = gene_id2contig_loc_dic[gene_id]["contig_id"]
        if cur_contig_id not in center_loc_dic:
            center_loc_dic[cur_contig_id] = []
        center_loc_dic[cur_contig_id].append(gene_id2contig_loc_dic[gene_id]["loc"])

    # return {contig_id:{gene_id:{"loc","start_pos","end_pos","strand"}}}
    contig_id2gene_ids_dic = _contig_id2gene_ids(
        contig_id_lst, client, prj_name, center_loc_dic, loc_window=neighbourhood_size
    )

    # return all uniq gene_ids
    ret = []
    for contig_id in contig_id2gene_ids_dic:
        ret = ret + list(contig_id2gene_ids_dic[contig_id].keys())
    return list(set(ret))


#######################
# user_functions
#######################


# return {FG:[protein_id]}
def get_FG2prot_id(
    FG_lst: list,
    DB_name: str,
    client: MongoClient,
    col_name: str = "FG_interproscan_Pfam",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    mycol = client[DB_name][col_name]

    ret = _extract_mapping_from_DB(
        search_lst=FG_lst,
        search_key="FG",
        ret_key="protein_id",
        target_collection=mycol,
        max_buffer_size=max_buffer_size,
    )
    return ret


# return {FG:[gc_id]}
## TODO:
def _search_protein_by_FG(
    queryFG: str,
    DB_name: str,
    client: MongoClient,
    col_name: str = "protein_pfam_domain",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):

    # return a list of protein_id
    pass


def get_FG2gc_id(
    FG_lst: list,
    DB_name: str,
    client: MongoClient,
    col_name: str = "FG_interproscan_Pfam",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    FG2prot_id_dic = get_FG2prot_id(
        FG_lst, DB_name, client, col_name, max_buffer_size=MAX_BUFFER_SIZE
    )
    ret = {}
    for FG in FG2prot_id_dic:
        cur_prot_id_lst = FG2prot_id_dic[FG]
        ret[FG] = []
        cur_prot_id2gc_id_dic = _prot_id2gc_id(
            cur_prot_id_lst, DB_name, client, max_buffer_size=MAX_BUFFER_SIZE
        )  # {prot_id:[gc_id_lst]}
        for cur_prot_id in cur_prot_id2gc_id_dic:
            ret[FG] = ret[FG] + cur_prot_id2gc_id_dic[cur_prot_id]

    return ret


# for each gc_id, get all the gene_id within the range of neighbourhood difinition
# return {gc_id:[gene_id]}
def get_gc_neighbourhood_directly_connected(
    prj_name: str, gc_id_lst: list, neighbourhood_size: int, client: MongoClient
):
    # 1 get all gene_id for each gc_id
    # return {gc_id:[gene_id_lst]}
    gc_id2gene_id_dic = _gc_id2gene_id(gc_id_lst, client, prj_name)

    ret = _query_neighbours(gc_id2gene_id_dic, neighbourhood_size, client, prj_name)
    return ret


# return {msp_id:{
# 	"gtdbtk_ann",
# 	"metaphlan_LM_spp",
# 	"metaphlan_LM_spp_r2"
# }}
def _msp_id2taxaAnn(
    msp_id_lst: list,
    DB_name: str,
    client: MongoClient,
    col_name: str = "msp_data",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    mycol = client[DB_name][col_name]
    if len(msp_id_lst) == 0:
        DB_query_res = mycol.find({})
    else:
        DB_query_res = mycol.find({"msp_id": {"$in": msp_id_lst}})
    ret = {}
    for item in DB_query_res:
        cur_msp_id = item["msp_id"]
        ret[cur_msp_id] = {}

        gtdbtk_ann = "NA"
        if len(item["gtdbtk_ann"]) > 1:
            gtdbtk_ann = item["gtdbtk_ann"][0]
        ret[cur_msp_id]["gtdbtk_ann"] = gtdbtk_ann

        metaphlan_LM_spp = "NA"
        metaphlan_LM_spp_r2 = "NA"
        if "metaphlan_LM_spp" in item["metaphlan_LM_spp"]:
            metaphlan_LM_spp = item["metaphlan_LM_spp"]["metaphlan_LM_spp"]
            metaphlan_LM_spp_r2 = round(float(item["metaphlan_LM_spp"]["r2"]), 4)
        ret[cur_msp_id]["metaphlan_LM_spp"] = metaphlan_LM_spp
        ret[cur_msp_id]["metaphlan_LM_spp_r2"] = metaphlan_LM_spp_r2

    return ret


#######################
# unsorted below
#######################


# A.2.e
# return {gc_id:FG}
def _gc_id2FG(gc_id_lst: list, DB_name: str, client: MongoClient):
    gc_id2prot_id_dic = _gc_id2prot_id(gc_id_lst, DB_name, client)
    # get all prot_ids and extract FGs
    prot_id_lst = list(set([gc_id2prot_id_dic[k] for k in gc_id2prot_id_dic]))
    prot_id2FG_dic = _prot_id2FG(prot_id_lst, DB_name, client)
    ret = {}
    for gc_id in gc_id_lst:
        prot_id = gc_id2prot_id_dic[gc_id]
        cur_FG = ""
        if prot_id in prot_id2FG_dic:
            cur_FG = prot_id2FG_dic[prot_id]
        ret[gc_id] = cur_FG
    return ret


# A.2.f
# return {gc_id:{gene_id:seq}}
def _gc_id2seq(
    gc_id_lst: list,
    prj_name: str,
    client: MongoClient,
):
    ret = {}
    # get all related gene_id for each gc_id
    # print("AUAS: get all related gene_id for each gc_id",prj_name)
    gc_id2gene_id_dic = _gc_id2gene_id(gc_id_lst, client, prj_name)
    for gc_id in gc_id2gene_id_dic:
        # print("AUAS: get seq for",gc_id)
        gene_id_lst = gc_id2gene_id_dic[gc_id]
        seq_dic = _gene_id2seq(gene_id_lst, client, prj_name)
        if gc_id not in ret:
            ret[gc_id] = {}
        for gene_id in gene_id_lst:
            if gene_id not in seq_dic:
                print(
                    "missing gene sequence for {0} in {1}, skipped".format(
                        gene_id, prj_name
                    )
                )
                continue
            ret[gc_id][gene_id] = seq_dic[gene_id]
    return ret


# A.2.f1
# return {gc_id:N_seq}
def _gc_id2seq_count(
    gc_id_lst: list,
    prj_name: str,
    client: MongoClient,
):
    ret = {}
    # get all related gene_id for each gc_id
    # print("AUAS: get all related gene_id for each gc_id",prj_name)
    gc_id2gene_id_dic = _gc_id2gene_id(gc_id_lst, client, prj_name)
    for gc_id in gc_id2gene_id_dic:
        ret[gc_id] = len(gc_id2gene_id_dic[gc_id])
    return ret


# A.2.g
# return {gc_id:{gene_id:sample_id}}
def _gc_id2sample_id(
    gc_id_lst: list,
    prj_name: str,
    client: MongoClient,
):
    ret = {}
    # get all related gene_id for each gc_id
    gc_id2gene_id_dic = _gc_id2gene_id(gc_id_lst, client, prj_name)
    for gc_id in gc_id2gene_id_dic:
        gene_id_lst = gc_id2gene_id_dic[gc_id]
        sample_id_dic = _gene_id2sample_id(
            gene_id_lst, client, prj_name
        )  # to improve effiency by merging all genes to search together!
        if gc_id not in ret:
            ret[gc_id] = {}
        for gene_id in gene_id_lst:
            if gene_id not in sample_id_dic:
                print(
                    "missing sample_id information for {0} in {1}, skipped".format(
                        gene_id, prj_name
                    )
                )
                continue
            ret[gc_id][gene_id] = sample_id_dic[gene_id]
    return ret

    # return {gc_id:{sample_id: if_assembled [bool]}}


def _gc_id2prev(
    gc_id_lst: list,
    prj_name: str,
    client: MongoClient,
):
    pre_ret = _gc_id2sample_id(
        gc_id_lst, prj_name, client
    )  # {gc_id:{gene_id:sample_id}}
    gc2sample_dic = {}
    for gc_id in pre_ret:
        gc2sample_dic[gc_id] = {}
        for sample_id in pre_ret.values():
            gc2sample_dic[sample_id] = (
                1  # {gc_id: {sample:1}} mark all samples with the gc_id assembled
            )

    # load metadata to get all sample ids
    metadata_dic = load_metadata(prj_name=prj_name, client=client)  # key is sample_id

    ret = {}
    for gc_id in gc2sample_dic:
        ret[gc_id] = {}
        for sample_id in sample_id_lst:
            ret[gc_id][sample_id] = sample_id in metadata_dic
    return ret  #


# A.2.h
# return {gc_id:{msp_id, gene_category, gtdbtk_ann, metaphlan_LM_ann, metaphlan_LM_ann_R2}}
def _gc_id2taxa_ann(gc_id_lst: list, prj_name: str, client: MongoClient):
    ret = {}
    # load the matching msp_id to each gc_id:
    gc_id2msp_id_dic = _gc_id2msp_id(gc_id_lst, client, prj_name)
    # get all relatd msp_id
    msp_id_lst = []
    for gc_id in gc_id2msp_id_dic:
        msp_id_lst.append(gc_id2msp_id_dic[gc_id][0])
    msp_id_lst = list(set(msp_id_lst))

    # load msp annotation info
    msp_ann_dic = _msp_id2taxaAnn(msp_id_lst, prj_name, client)

    for gc_id in gc_id_lst:
        msp_id, gene_category, gtdbtk_ann, metaphlan_LM_spp, metaphlan_LM_ann_R2 = [
            "NA",
            "NA",
            "NA",
            "NA",
            "NA",
        ]
        ret[gc_id] = {}
        if gc_id not in gc_id2msp_id_dic:
            ret[gc_id]["msp_id"] = msp_id
            ret[gc_id]["gene_category"] = gene_category
            ret[gc_id]["gtdbtk_ann"] = gtdbtk_ann
            ret[gc_id]["metaphlan_LM_ann"] = metaphlan_LM_spp
            ret[gc_id]["metaphlan_LM_ann_R2"] = metaphlan_LM_ann_R2
            continue
        msp_id, gene_category = gc_id2msp_id_dic[gc_id]
        if msp_id not in msp_ann_dic:
            ret[gc_id]["msp_id"] = msp_id
            ret[gc_id]["gene_category"] = gene_category
            ret[gc_id]["gtdbtk_ann"] = gtdbtk_ann
            ret[gc_id]["metaphlan_LM_ann"] = metaphlan_LM_spp
            ret[gc_id]["metaphlan_LM_ann_R2"] = metaphlan_LM_ann_R2
            continue
        ret[gc_id]["msp_id"] = msp_id
        ret[gc_id]["gene_category"] = gene_category
        ret[gc_id]["gtdbtk_ann"] = msp_ann_dic[msp_id]["gtdbtk_ann"]
        ret[gc_id]["metaphlan_LM_ann"] = msp_ann_dic[msp_id]["metaphlan_LM_spp"]
        ret[gc_id]["metaphlan_LM_ann_R2"] = msp_ann_dic[msp_id]["metaphlan_LM_spp_r2"]

    return ret


# A.3.a


# A.3.b
# return {prot_id:FG}
def _prot_id2FG(
    prot_id_lst: list,
    DB_name: str,
    client: MongoClient,
    col_name: str = "FG_interproscan_Pfam",
):
    mycol = client[DB_name][col_name]
    DB_query_res = mycol.find({"protein_id": {"$in": prot_id_lst}})
    ret = {}
    for item in DB_query_res:
        ret[item["protein_id"]] = item["FG"]
    return ret


# A.4.a


# A.4.b
# return {msp_id:{sample_id: RPKM}}
def _msp_id2RPKM(
    msp_id_lst: list, DB_name: str, client: MongoClient, col_name: str = "msp_data"
):
    ret = {}
    mycol = client[DB_name][col_name]
    DB_query_res = mycol.find({"msp_id": {"$in": msp_id_lst}})
    for item in DB_query_res:
        cur_msp_id = DB_query_res["msp_id"]
        ret[cur_msp_id] = DB_query_res["mgx_RPKM"]
        del ret[cur_msp_id]["N"]
    return ret


# A.4.c
# return {gc_id:gene_category}
def _msp_id2gc_id(
    msp_id: str, client: MongoClient, DB_name: str, col_name: str = "gc_to_msp"
):
    mycol = client[DB_name][col_name]
    DB_query_res = mycol.find({"msp_id": msp_id})
    ret = {}
    for item in DB_query_res:
        ret[item["gc_id"]] = item["gene_category"]
    return ret


# return {sample_id:{prj_name,disease_state,disease_grp}}
# {
#  prj_name,
#  metadata:{
# 		sample_id:{
# 			disease_group,
# 			disease_state,
# 			ethnic_group,
# 			gender,
# 			age,
# 			use_antibiotic,
# 			QC_status:{
# 				raw_pairs,
# 				trimmed_pairs,
# 				decontaminated_hg38_pairs,
# 				final_pairs
# 			}
# 		}
# 	}
# }


# B.1.a
def load_metadata(
    prj_name: str,
    client: MongoClient,
    DB_name: str = "multi_cohorts_db",
    col_name: str = "mcmd",
):

    mycol = client[DB_name][col_name]
    lst = list(mycol.find({"prj_name": prj_name}))

    if not len(lst) == 1:  # expected only one match for each project
        print("failed fetch metadata for: {0}\n".format(prj_name), lst)
        return {}

    return lst[0]["metadata"]


# B.1.b
# return {N_genes, N_gf}
def count_records(
    prj_name: str,
    client: MongoClient,
    nucleotide_clstr_name: str = "nucleotide_clstr",
    protine_clstr_name: str = "prot90_clstr",
):
    # prj_name = "Healthy_500FG"
    # client = MongoClient("localhost", 27017, maxPoolSize=50)
    ret = {}
    nucleotide_clstr_col = client[prj_name][nucleotide_clstr_name]
    ret["N_genes"] = nucleotide_clstr_col.count_documents({})
    protine_clstr_col = client[prj_name][protine_clstr_name]
    ret["N_gfs"] = protine_clstr_col.count_documents({})
    ret["prj_name"] = prj_name
    return ret


# B.1.c
# save return as _msp_id2taxaAnn, load all msp's information!
def load_all_MSP(prj_name: str, client: MongoClient, col_name: str = "msp_data"):
    return _msp_id2taxaAnn(
        msp_id_lst=[], DB_name=prj_name, client=client, col_name=col_name
    )


# B.2.a
def save_gc_RPKM_tsv(prj_name: str, gc_id_lst: list, client: MongoClient, sfp: str):

    # load RPKM table
    rpkm_dic = _gc_id2RPKM(gc_id_lst, prj_name, client)
    # save to tsv file
    _save_RPKM_tsv(rpkm_dic, prj_name, client, sfp)

    return


# helper function
# RPKM_dic {feature_id:{sample_id:RPKM}}
def _save_RPKM_tsv(rpkm_dic, prj_name, client, sfp):
    # load metadata to get all sample ids
    metadata_dic = load_metadata(prj_name=prj_name, client=client)
    sample_id_lst = list(metadata_dic.keys())
    sample_id_lst.sort()

    feature_id_lst = list(rpkm_dic.keys())
    feature_id_lst.sort()
    # file to save the result
    sf = open(sfp, "w")
    sf.write("\t".join(["sample_id"] + feature_id_lst) + "\n")
    for sample_id in sample_id_lst:
        cur_rec = [sample_id]
        for f_id in feature_id_lst:
            if f_id not in rpkm_dic:
                print("[WARNNING] missing RPKM record for:", f_id)
                continue
            cur_val = "0."
            if sample_id in rpkm_dic[f_id]:
                cur_val = str(rpkm_dic[f_id][sample_id])
            cur_rec.append(cur_val)
        sf.write("\t".join(cur_rec) + "\n")
    sf.close()
    return


# B.2.b
# save matched MSP info to a tsv file: header=[gc_id, msp_id, gene_category, gtdbtk_ann, metaphlan_LM_ann, metaphlan_LM_ann_R2]
def save_gc_taxa_ann_tsv(prj_name: str, gc_id_lst: list, client: MongoClient, sfp: str):
    ann_header_lst = [
        "msp_id",
        "gene_category",
        "gtdbtk_ann",
        "metaphlan_LM_ann",
        "metaphlan_LM_ann_R2",
    ]
    gc_id2taxa_ann_dic = _gc_id2taxa_ann(gc_id_lst, prj_name, client)
    # write header and save records
    sf = open(sfp, "w")
    sf.write("\t".join(["gc_id"] + ann_header_lst) + "\n")

    for gc_id in gc_id_lst:
        cur_rec = [gc_id] + [str(gc_id2taxa_ann_dic[gc_id][h]) for h in ann_header_lst]
        sf.write("\t".join(cur_rec) + "\n")
    sf.close()

    return


# B.2.c
# save the corresponding FG to a tsv file: header=[gc_id,prot_id,FG]
def save_gc_FG_tsv(
    prj_name: str,
    gc_id_lst: list,
    client: MongoClient,
    sfp: str,
    write_mode="w",
    save_header=True,
):
    # get matched prot_id first:
    gc_id2prot_id_dic = _gc_id2prot_id(gc_id_lst, prj_name, client)

    # get all prot_ids and extract FGs
    gc_id2FG_dic = _gc_id2FG(gc_id_lst, prj_name, client)  # {gc_id:FG}

    # write header and save to tsv files
    sf = open(sfp, write_mode)
    if save_header:
        sf.write("\t".join(["gc_id", "prot_id", "FG"]) + "\n")
    for gc_id in gc_id_lst:
        prot_id, FG = ["NA", "NA"]
        if gc_id not in gc_id2prot_id_dic:
            print("[WARNNING] no matched prot_id for gc_id:", gc_id)
            sf.write("\t".join([gc_id, prot_id, FG]) + "\n")
            continue
        prot_id = gc_id2prot_id_dic[gc_id]
        if gc_id not in gc_id2FG_dic:
            sf.write("\t".join([gc_id, prot_id, FG]) + "\n")
            continue
        FG = gc_id2FG_dic[gc_id]
        sf.write("\t".join([gc_id, prot_id, FG]) + "\n")
    sf.close()
    return


# B.2.d
# save the DNA sequences for all the genes matched to the gc_id_lst into a fasta file, header format: >{gene_id} {gc_id} {sample_id}
def save_gc_seq_fna(
    prj_name: str,
    gc_id_lst: list,
    client: MongoClient,
    sfp: str,
):
    sf = open(sfp, "w")
    # get seq for each genes belonging to any gc_id in the gc_id_lst
    gc_id2seq_dic = _gc_id2seq(
        gc_id_lst, prj_name, client
    )  # return {gc_id:{gene_id:seq}}
    gc_id2sample_id_dic = _gc_id2sample_id(
        gc_id_lst, prj_name, client
    )  # return {gc_id:{gene_id:sample_id}}
    for gc_id in gc_id2seq_dic:
        gene_id_lst = list(gc_id2seq_dic[gc_id].keys())
        for gene_id in gene_id_lst:
            sample_id = gc_id2sample_id_dic[gc_id][gene_id]
            seq = gc_id2seq_dic[gc_id][gene_id]
            sf.write(">{0} {1} {2}\n".format(gene_id, gc_id, sample_id))
            sf.write(seq + "\n")
    sf.close()
    return


# B.3.a
# return {domain:[prot_id]}
def get_prot_by_domain(
    prj_name: str,
    domain_lst: list,
    client: MongoClient,
    col_name: str = "protein_pfam_domain",
):
    ret = {}

    mycol = client[prj_name][col_name]
    DB_query_res = mycol.find({"Pfam_id": {"$in": domain_lst}})
    ret = {}
    for item in DB_query_res:
        cur_domain = item["Pfam_id"]
        if cur_domain not in ret:
            ret[cur_domain] = []
        ret[cur_domain].append(item["protein_id"])
    return ret


# B.3.b
# return {domain:[gc_id]}
def get_gc_by_domain(
    prj_name: str,
    domain_lst: list,
    client: MongoClient,
    col_name: str = "protein_pfam_domain",
):
    ret = {}
    prot_res = get_prot_by_domain(prj_name, domain_lst, client, col_name)

    for k in prot_res:
        prot_id_lst = prot_res[k]
        # return {prot_id:[gc_id_lst]}
        prot_id2gc_id_dic = _prot_id2gc_id(prot_id_lst, prj_name, client)
        ret[k] = []
        for prot_id in prot_id2gc_id_dic:
            ret[k] = ret[k] + prot_id2gc_id_dic[prot_id]
    return ret


# return {domain:domain description}
def _get_domain_by_pattern(
    prj_name: str,
    search_patterns: list,
    client: MongoClient,
    col_name: str = "protein_pfam_domain",
):
    ret = {}
    search_patterns = [it.lower() for it in search_patterns]
    mycol = client[prj_name][col_name]
    DB_query_res = mycol.find()
    ret = {}
    for ii, item in enumerate(DB_query_res):
        # if ii%1000000==0:
        # 	print(ii)
        for p in search_patterns:
            if p in item["Pfam_description"].lower():
                ret[item["Pfam_id"]] = item["Pfam_description"]
                break
    return ret


# return {domain:domain description}
def _get_domain_by_id(
    prj_name: str,
    domain_lst: list,
    client: MongoClient,
    col_name: str = "protein_pfam_domain",
):
    ret = {}

    mycol = client[prj_name][col_name]
    DB_query_res = mycol.find({"Pfam_id": {"$in": domain_lst}})
    ret = {}
    for item in DB_query_res:
        cur_domain = item["Pfam_id"]
        ret[cur_domain] = item["Pfam_description"]
    return ret


def get_mobile_element_domains(
    prj_name: str,
    client: MongoClient,
    search_patterns: list = ["integrase", "integrative", "transposon"],
):
    return _get_domain_by_pattern(prj_name, search_patterns, client)


def get_bacterial_toxin_domains(
    prj_name: str, client: MongoClient, search_patterns: list = ["toxin"]
):
    return _get_domain_by_pattern(prj_name, search_patterns, client)


# prj_name = "IBD_Elinav_2022" # IBD_Elinav_2022  IBD_HMP2
# get_prot_by_domain(prj_name,search_patterns,client)
# get_mobile_element_domains(prj_name,client)
if __name__ == "__main__":
    import os, glob, sys
    from pymongo import MongoClient
    from typing import List, Tuple, Dict, Optional
    import pandas as pd
    from scripts.single_cohort_API import get_prot_by_domain, get_gc_neighbourhood

    client = MongoClient("localhost", 27017, maxPoolSize=50)
    merged_taxa_fp = "/nfs/data/projects/multi_cohorts_db/genes_and_pathways/t6ss/t6ss.faa.FG.txt.taxa/merged.taxa.txt"
    gc_id_lst = [
        "Germany__3208__k119_68986::3::672::1016::+",
        "Israel__10051__k141_201970::2::1027::2475::+",
        "Israel__10051__k141_67351::14::17847::18332::+",
        "Israel__10051__k141_67351::15::18344::19813::+",
        "Israel__10051__k141_67351::5::5772::6191::+",
        "Israel__10051__k141_67351::6::6194::7912::+",
        "Israel__10051__k141_67351::8::9034::10698::+",
        "US__11515__k141_304139::3::2817::4385::+",
        "US__11515__k141_304139::4::4404::4982::+",
        "US__11515__k141_304139::5::5084::6559::+",
        "US__11515__k141_304139::6::6563::6952::+",
        "US__11515__k141_304139::7::6949::8724::+",
        "US__11518__k141_36348::2::215::733::-",
        "US__11524__k141_126756::26::29871::31928::-",
    ]  # gc_id from IBD_Elinav_2022 in merged_taxa_fp
    prj_name = "IBD_Elinav_2022"
    gc_neighbourhood = get_gc_neighbourhood(
        prj_name, gc_id_lst, neighbourhood_size=3, client=client
    )

    # {FG:{cohort:gene}}
    # input_fp_lst = glob.glob(os.path.join("/nfs/data/projects/mongoDB_data/MCPC/cohort_list","*.out"))
    # for input_fp in input_fp_lst:
    # 	print("working with:",input_fp)
    # 	_extract_FG_of_MCPC(input_fp=input_fp,client=client)
    # ann_fp_lst = glob.glob(os.path.join("/nfs/data/projects/mongoDB_data/MCPC/cohort_list","*.out.FG"))
    # _Pfam_MCGC(client,ann_fp_lst,DB_name="multi_cohorts_db")
    # metadata_dic = load_metadata(prj_name="IBD_HMP2",client=client)

    # gc_id_lst= ['cohort_merged__HV001__k119_79506::11::14828::18115::-','cohort_merged__HV001__k119_79506::7::7365::9515::+']
    # DB_name = "Healthy_500FG"
    # rpkm_dic = _gc_id2RPKM(gc_id_lst,DB_name,client)

    # prj_name = "Healthy_500FG"

    # msp_id_lst = ["msp_001","msp_002"]
    # msp_ann = _msp_id2taxaAnn(msp_id_lst,DB_name,client)

    # prj_name = "IBD_HMP2"
    # gc_id_lst = ["cohort_merged__CSM5MCVN_k105_13509::5::3124::14730::-","cohort_merged__MSM6J2ML_k105_21214::3::1703::13423::+"]
    # gc_id2prot_id_dic = _gc_id2prot_id(gc_id_lst,DB_name,client)
    #
    # prj_name = "IBD_RobertHM_2022"
    # fp = "/nfs/data/projects/{0}/analysis/assembly/genes_and_pathways/uniprotkb_GMPS_AZA6MP_related.faa.{0}.70cov_30id.blastp.out".format(prj_name)
    # fp = "/nfs/data/projects/{0}/analysis/assembly/genes_and_pathways/uniprotkb_HPRT_AZA6MP_related.faa.{0}.70cov_30id.blastp.out".format(prj_name)
    # fp = "/nfs/data/projects/{0}/analysis/assembly/genes_and_pathways/uniprotkb_IMPDH_AZA6MP_related.faa.{0}.70cov_30id.blastp.out".format(prj_name)
    # fp = "/nfs/data/projects/{0}/analysis/assembly/genes_and_pathways/Xdh_AZA6MP_related.faa.{0}.70cov_30id.blastp.out".format(prj_name)
    # fp = "/nfs/data/projects/{0}/analysis/assembly/genes_and_pathways/TPMT_AZA6MP_related.faa.{0}.70cov_30id.blastp.out".format(prj_name)

    # fp = "/nfs/data/projects/Healthy_500FG/analysis/assembly/genes_and_pathways/P_vulgatus_3_oxo_5_alpha_steroid_4_dehydrogenae.faa.Healthy_500FG.70cov_30id.blastp.out"
    # fp = "/nfs/data/projects/Healthy_500FG/analysis/assembly/genes_and_pathways/P_vulgatus_beta-glucuronidae.faa.Healthy_500FG.70cov_30id.blastp.out"
    # fp = "/nfs/data/projects/IBD_Elinav_2022/analysis/assembly/cdhit/cohort_prot90.faa.IBD_Elinav_2022.cov50_id80.tblastn.out"
    fp = "/nfs/data/projects/IBD_HMP2/analysis/assembly/genes_and_pathways/hep_hag_uniprot_May14.faa.IIBD_HMP2.80cov_30id_eval_0.01.tblastn.out"
    blast_res = parse_fmt6_blast_results(fp)
    # # merge all gc_ids from the blast_res
    # gc_id_lst = []
    # for q in blast_res:
    # 	gc_id_lst = gc_id_lst+blast_res[q]client
    # gc_id_lst = list(set(gc_id_lst))

    # print(len(gc_id_lst),"unique gc_id loaded")

    # # merge all gc_ids from the blastp_res
    # prot90_id_lst = []
    # for q in blast_res:
    # 	prot90_id_lst = prot90_id_lst+blast_res[q]
    # prot90_id_lst = list(set(prot90_id_lst))

    # print(len(prot90_id_lst),"unique prot90_id loaded")
    # # matched to gc_id
    # prot_id2gc_id_dic = _prot_id2gc_id(prot90_id_lst,prj_name,client)
    # # return {prot_id:[gc_id_lst]}
    # gc_id_lst = []
    # for prot_id in prot_id2gc_id_dic:
    # 	gc_id_lst = gc_id_lst+prot_id2gc_id_dic[prot_id]
    # gc_id_lst = list(set(gc_id_lst))
    # print(len(gc_id_lst),"unique gc_id loaded")

    gc_id_lst = []
    for q in blast_res:
        gc_id_lst = gc_id_lst + blast_res[q]
    print(len(gc_id_lst), "unique gc_id loaded")

    # test save to RPKM tsv file
    prj_name = "IBD_HMP2"
    sfp = fp + ".rpkm"
    save_gc_RPKM_tsv(prj_name, gc_id_lst, client, sfp)

    # test save to msp_annotation tsv file
    sfp = fp + ".msp_ann"
    save_gc_taxa_ann_tsv(prj_name, gc_id_lst, client, sfp)

    # test save to FG tsv file
    sfp = fp + ".FG_ann"
    save_gc_FG_tsv(prj_name, gc_id_lst, client, sfp)

    ######################################################################
    # test sequence extraction:
    ######################################################################
    # prj_name = "IBD_Elinav_2022"
    # gc_id_lst = ["Germany__6219__k141_9418::47::43375::43857::-"]
    # sfp = "/nfs/data/projects/IBD_Elinav_2022/analysis/assembly/genes_and_pathways/kp_GPX4__{0}.fna".format(prj_name)

    # prj_name = "IBD_HMP2"
    # gc_id_lst = ["cohort_merged__CSM79HIB_k105_28270::8::8300::8782::+"]
    # sfp = "/nfs/data/projects/IBD_Elinav_2022/analysis/assembly/genes_and_pathways/kp_GPX4__{0}.fna".format(prj_name)

    # prj_name = "IBD_HallAB_2017"
    # gc_id_lst = ["cohort_merged__SRR5650032__k99_53560::6::6133::6615::+"]
    # sfp = "/nfs/data/projects/IBD_Elinav_2022/analysis/assembly/genes_and_pathways/kp_GPX4__{0}.fna".format(prj_name)

    # prj_name = "IBD_HeQ_2017"
    # gc_id_lst = ["cohort_merged__ERR1620279__k99_543::8::8225::8707::+"]
    # sfp = "/nfs/data/projects/IBD_Elinav_2022/analysis/assembly/genes_and_pathways/kp_GPX4__{0}.fna".format(prj_name)

    # prj_name = "IBD_IjazUZ_2017"
    # gc_id_lst = ["cohort_merged__ERR1776102__k141_32268::8::8300::8782::+"]
    # sfp = "/nfs/data/projects/IBD_Elinav_2022/analysis/assembly/genes_and_pathways/kp_GPX4__{0}.fna".format(prj_name)

    # prj_name = "IBD_LewisJD_2015"
    # gc_id_lst = ["cohort_merged__SRR2145578__k99_20008::8::8177::8659::+","cohort_merged__SRR2145635__k99_1910::5::2811::3293::-"]
    # sfp = "/nfs/data/projects/IBD_Elinav_2022/analysis/assembly/genes_and_pathways/kp_GPX4__{0}.fna".format(prj_name)

    # prj_name = "CRC_YangY_2021"
    # gc_id_lst = ["cohort_merged__SRR16124185__k141_35005::7::7719::8201::+"]
    # sfp = "/nfs/data/projects/IBD_Elinav_2022/analysis/assembly/genes_and_pathways/kp_GPX4__{0}.fna".format(prj_name)

    # prj_name = "CRC_YachidaS_2019"
    # gc_id_lst = ["cohort_merged__DRR127628__k141_149657::3::1505::1987::-",
    # 			"cohort_merged__DRR127567__k141_44292::53::46904::47386::-",
    # 			"cohort_merged__DRR127557__k141_3223::4::2896::3378::-"]
    # sfp = "/nfs/data/projects/IBD_Elinav_2022/analysis/assembly/genes_and_pathways/kp_GPX4__{0}.fna".format(prj_name)
    # save_gc_seq_fna(prj_name,gc_id_lst,client,sfp)

    #


#
