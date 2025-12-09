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


## TODO:
def search_protein_by_FG(
    queryFG: str,
    DB_name: str,
    client: MongoClient,
    col_name: str = "protein_pfam_domain",
    max_buffer_size: int = MAX_BUFFER_SIZE,
):

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
