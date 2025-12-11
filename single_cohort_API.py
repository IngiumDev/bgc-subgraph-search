import os, glob, sys
from collections import defaultdict

from pymongo import MongoClient
from pymongo.collection import Collection as pymongoCollection
from typing import List, Tuple, Dict, Optional
import warnings
import re

from utils import count_prots_per_pfam, matches_fg_query

MAX_BUFFER_SIZE = 5000  # split the input search list to avoid

### helper functions

"""
Search the values in search_lst against the search_key field in the target_collection
extract the values in the ret_key field
return {search_value:[target_val]}
"""


def search_protein_by_FG(
    queryFG: str,
    DB_name: str,
    client: MongoClient,
    col_name: str = "protein_pfam_domain",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):
    """

    :param queryFG:
        FG query string in the new regular-expression-like format, describing
        a Pfam domain architecture.
        Example: "*, (PF05658 | PF05662){1,20}, *, PF03895"

    :param DB_name:
        Name of the MongoDB database containing the Pfam and FG collections
        for the project to be searched.

    :param client:
        An active `MongoClient` instance connected to the MongoDB server.

    :param col_name:
        Name of the collection that stores per-domain Pfam annotations
        (`protein_pfam_domain`-style schema: one document per Pfam hit).

    :param max_buffer_size:
        Maximum number of entries to process in a single batch when
        splitting large queries or result sets (not strictly enforced in
        the first draft implementation, but reserved for future scaling).

    :return:
        A list of `protein_id` strings whose FG matches `queryFG`.
    """
    # Check if query is valid

    # extract mandatory domains
    # TODO: Replace with actual name
    mandatory_domains = get_mandatory_domains(queryFG)
    if len(mandatory_domains) == 0:
        raise ValueError("No mandatory domains found in the FG query.")
    # TODO: Tests

    # Query Database for mandatory domains

    prots_per_pfam = count_prots_per_pfam(mandatory_domains) # Finalize parameters

    # Get mandatory domain with least proteins

    min_pfam = min(prots_per_pfam, key=prots_per_pfam.get)



    # Query all the proteins that contain this PFAM domain
    collection = client[DB_name][col_name]

    # Extract all protein IDs that have this PFAM domain
    protein_ids = collection.distinct("protein_id", {"Pfam_id": min_pfam})


    matching_protein_ids = []
    compact_protein_collection = client[DB_name]["FG_interproscan_Pfam"]

    # TODO: Refactor to only parse and compile the query once
    for protein_id in protein_ids:
        # TODO: Optimize this by only querying once for all proteins
        # Query the compact protein collection to get the PFAM domain string for this protein_id, It's in the FG field and they are separeted by :::
        document = compact_protein_collection.find_one({
            "protein_id": protein_id,
            "FG": 1,
            "_id": 0
        }) # Only want the FG, what if there is none?
        if document is None:
            raise ValueError(f"No document found for protein_id: {protein_id}")

        protein_Pfam_domains = document["FG"]
        protein_Pfam_domains_string = protein_Pfam_domains.replace(":::", "")

        if matches_fg_query(queryFG, protein_Pfam_domains_string):
            matching_protein_ids.append(protein_id)


    return matching_protein_ids

        # print(protein_id)


    # Init matching result list i.e. matchingProteinIDs

    # Iterate through each protein ID and create PFAM string from this and then use regular regex to match

    # return a list of protein_id
    pass

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

if __name__ == "__main__":
    # testing
    db_host = str(sys.argv[1])
    db_port = int(sys.argv[2])
    user_name = str(sys.argv[3])
    user_pd = str(sys.argv[4])

    # search_protein_by_FG(TODO)

