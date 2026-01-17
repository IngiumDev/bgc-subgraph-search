import os, glob, sys
from collections import defaultdict
import argparse

from pymongo import MongoClient
from pymongo.collection import Collection as pymongoCollection
from typing import List, Tuple, Dict, Optional
import warnings
import re
import logging

from utils import count_prots_per_pfam, matches_fg_query, get_mandatory_domains

MAX_BUFFER_SIZE = 5000  # split the input search list to avoid

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s: %(message)s"))
logger.addHandler(handler)

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
    # Log starting query and log the parameters for the method including queryFG, DB_name, col_name, max_buffer_size
    logger.info(
        f"Starting FG search with parameters: queryFG={queryFG}, DB_name={DB_name}, col_name={col_name}, max_buffer_size={max_buffer_size}")

    # TODO: Check if query is valid

    # extract mandatory domains
    logger.info("Extracting mandatory domains from FG query.")
    mandatory_domains = get_mandatory_domains(queryFG)
    if len(mandatory_domains) == 0:
        raise ValueError("No mandatory domains found in the FG query.")
    # TODO: Tests
    logger.info(f"Mandatory domains extracted: {mandatory_domains}")

    # Query Database for mandatory domains
    logger.info("Counting proteins per mandatory Pfam domain.")
    prots_per_pfam = count_prots_per_pfam(mandatory_domains, client)  # Finalize parameters
    logger.info(f"Proteins per Pfam domain: {prots_per_pfam}")

    # Get mandatory domain with least proteins

    min_pfam = min(prots_per_pfam, key=prots_per_pfam.get)
    logger.info(f"Pfam domain with least proteins: {min_pfam} ({prots_per_pfam[min_pfam]} proteins)")

    # Query all the proteins that contain this PFAM domain
    collection = client[DB_name][col_name]

    # Extract all protein IDs that have this PFAM domain
    logger.info(f"Querying proteins with Pfam domain: {min_pfam}")
    protein_ids = collection.distinct("protein_id", {"Pfam_id": min_pfam})
    logger.info(f"Number of proteins with Pfam domain {min_pfam}: {len(protein_ids)}")

    matching_protein_ids = []
    compact_protein_collection = client[DB_name]["FG_interproscan_Pfam"]

    FG_query = queryFG.replace(" *", ".*").replace(", ", "").replace(" ", "").replace(",", "")
    if FG_query.startswith("*"):
      FG_query = f".{FG_query}" # handle edge case
      
    logger.info(f"Compiled FG query regex: {FG_query}")
    pattern = re.compile(FG_query)
    logger.info("Searching for matching proteins based on FG query using batched $in queries.")

    # Process protein_ids in batches to reduce the number of queries.
    batch_size = 100
    total = len(protein_ids)
    for start in range(0, total, batch_size):
        end = min(start + batch_size, total)
        batch = protein_ids[start:end]
        logger.info(f"Processing proteins {start + 1}-{end}/{total} (batch size {len(batch)})")

        # Fetch all documents for this batch with a single $in query
        cursor = compact_protein_collection.find(
            {"protein_id": {"$in": batch}}, {"protein_id": 1, "FG": 1, "_id": 0}
        )

        # Map returned documents by protein_id for quick lookup
        found_docs = {}
        for doc in cursor:
            pid = doc.get("protein_id")
            if pid is None:
                raise ValueError("Document missing protein_id field during batch fetch")
            found_docs[pid] = doc

        # If any protein_id in the batch is missing from the returned docs, raise same error as before
        missing = [pid for pid in batch if pid not in found_docs]
        if missing:
            raise ValueError(f"No document found for protein_id(s): {', '.join(missing)}")

        # Iterate in original order to preserve ordering and logging
        for idx_in_batch, protein_id in enumerate(batch, start=start):
            document = found_docs[protein_id]
            protein_Pfam_domains = document.get("FG")
            protein_Pfam_domains_string = protein_Pfam_domains.replace(":::", "")

            logger.debug(f"Protein PFAM domain string: {protein_Pfam_domains_string}")
            if matches_fg_query(pattern, protein_Pfam_domains_string):
                logger.debug(f"Protein {protein_id} matches the FG query.")
                matching_protein_ids.append(protein_id)

    logger.info(f"Total matching proteins found: {len(matching_protein_ids)}")

    # NOTE: Here I added what still has to be done
    # right now we have query -> prot_ids (in theory we could have this for several queries) 
    # in the rest of get_FG2gc_id this is what happenes (we need to implement this):
    # FG2prot_id_dic = get_FG2prot_id(FG_lst,DB_name,client,col_name,max_buffer_size=MAX_BUFFER_SIZE) # old func which we replace (both return p_ids)
    #
    # ret = {}
    # for FG in FG2prot_id_dic:
    # 	cur_prot_id_lst = FG2prot_id_dic[FG] # all prots matching query
    #     # log: 
    #     # FG2prot_id_dic (first entry)
    #     # PF02606: ['Germany__4249__k119_38096::141::143966::145006::-', 'Germany__4249__k119_11459::3::3403::4437::+', ... ]
    # 	ret[FG] = []
    # 	cur_prot_id2gc_id_dic = _prot_id2gc_id(cur_prot_id_lst,DB_name,client,max_buffer_size=MAX_BUFFER_SIZE) # {prot_id:[gc_id_lst]}
    #     # log: 
    #     # cur_prot_id2gc_id_dic
    #     # Germany__4249__k119_38096::141::143966::145006::- ['Germany__4249__k119_38096::141::143966::145006::-', 'US__11562__k141_79879::2::843::1766::+', 'Israel__10088__k141_85188::6::4068::4940::+']
    # 	for cur_prot_id in cur_prot_id2gc_id_dic:
    # 		ret[FG] = ret[FG] + cur_prot_id2gc_id_dic[cur_prot_id]
    # 	
    # return ret

    return matching_protein_ids

    # print(protein_id)

    # Init matching result list i.e. matchingProteinIDs

    # Iterate through each protein ID and create PFAM string from this and then use regular regex to match

    # return a list of protein_id

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



def search_protein_by_FG_old(
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
    # Log starting query and log the parameters for the method including queryFG, DB_name, col_name, max_buffer_size
    logger.info(
        f"Starting FG search with parameters: queryFG={queryFG}, DB_name={DB_name}, col_name={col_name}, max_buffer_size={max_buffer_size}")

    # TODO: Check if query is valid

    # extract mandatory domains
    logger.info("Extracting mandatory domains from FG query.")
    mandatory_domains = get_mandatory_domains(queryFG)
    if len(mandatory_domains) == 0:
        raise ValueError("No mandatory domains found in the FG query.")
    # TODO: Tests
    logger.info(f"Mandatory domains extracted: {mandatory_domains}")

    # Query Database for mandatory domains
    logger.info("Counting proteins per mandatory Pfam domain.")
    prots_per_pfam = count_prots_per_pfam(mandatory_domains, client)  # Finalize parameters
    logger.info(f"Proteins per Pfam domain: {prots_per_pfam}")

    # Get mandatory domain with least proteins

    min_pfam = min(prots_per_pfam, key=prots_per_pfam.get)
    logger.info(f"Pfam domain with least proteins: {min_pfam} ({prots_per_pfam[min_pfam]} proteins)")

    # Query all the proteins that contain this PFAM domain
    collection = client[DB_name][col_name]

    # Extract all protein IDs that have this PFAM domain
    logger.info(f"Querying proteins with Pfam domain: {min_pfam}")
    protein_ids = collection.distinct("protein_id", {"Pfam_id": min_pfam})
    logger.info(f"Number of proteins with Pfam domain {min_pfam}: {len(protein_ids)}")

    matching_protein_ids = []
    compact_protein_collection = client[DB_name]["FG_interproscan_Pfam"]

    pattern = compile_fg_query_pattern(queryFG)
    
    logger.info("Searching for matching proteins based on FG query.")
    for idx, protein_id in enumerate(protein_ids):
        # TODO: Optimize this by only querying once for all proteins
        # Query the compact protein collection to get the PFAM domain string for this protein_id, It's in the FG field and they are separeted by :::
        if idx % 100 == 0:
            logger.info(f"Processing protein_id {idx + 1}/{len(protein_ids)}: {protein_id}")
        document = compact_protein_collection.find_one({
            "protein_id": protein_id},
            {"FG": 1,
             "_id": 0
             })  # Only want the FG, what if there is none?
        if document is None:
            raise ValueError(f"No document found for protein_id: {protein_id}")

        protein_Pfam_domains = document["FG"]

        protein_Pfam_domains_string = protein_Pfam_domains.replace(":::", "")

        logger.debug(f"Protein PFAM domain string: {protein_Pfam_domains_string}")
        if matches_fg_query(pattern, protein_Pfam_domains_string):
            logger.debug(f"Protein {protein_id} matches the FG query.")
            matching_protein_ids.append(protein_id)

    logger.info(f"Total matching proteins found: {len(matching_protein_ids)}")

    # NOTE: Here I added what still has to be done
    # right now we have query -> prot_ids (in theory we could have this for several queries) 
    # in the rest of get_FG2gc_id this is what happenes (we need to implement this):
    # FG2prot_id_dic = get_FG2prot_id(FG_lst,DB_name,client,col_name,max_buffer_size=MAX_BUFFER_SIZE) # old func which we replace (both return p_ids)
    #
    # ret = {}
    # for FG in FG2prot_id_dic:
    # 	cur_prot_id_lst = FG2prot_id_dic[FG] # all prots matching query
    #     # log: 
    #     # FG2prot_id_dic (first entry)
    #     # PF02606: ['Germany__4249__k119_38096::141::143966::145006::-', 'Germany__4249__k119_11459::3::3403::4437::+', ... ]
    # 	ret[FG] = []
    # 	cur_prot_id2gc_id_dic = _prot_id2gc_id(cur_prot_id_lst,DB_name,client,max_buffer_size=MAX_BUFFER_SIZE) # {prot_id:[gc_id_lst]}
    #     # log: 
    #     # cur_prot_id2gc_id_dic
    #     # Germany__4249__k119_38096::141::143966::145006::- ['Germany__4249__k119_38096::141::143966::145006::-', 'US__11562__k141_79879::2::843::1766::+', 'Israel__10088__k141_85188::6::4068::4940::+']
    # 	for cur_prot_id in cur_prot_id2gc_id_dic:
    # 		ret[FG] = ret[FG] + cur_prot_id2gc_id_dic[cur_prot_id]
    # 	
    # return ret

    return matching_protein_ids
#

def compile_fg_query_pattern(queryFG: str) -> re.Pattern[str]:
	"""Convert an FG query string into a compiled regex pattern."""
	cleaned_query = (queryFG or "").strip()
	if not cleaned_query:
		raise ValueError("FG query cannot be empty.")

	fg_regex = (
		cleaned_query.replace(" *", ".*")
		.replace(", ", "")
		.replace(" ", "")
		.replace(",", "")
	)

	if fg_regex.startswith("*"):
		fg_regex = f".{fg_regex}"

	try:
		return re.compile(fg_regex)
	except re.error as exc:
		raise ValueError(
			f"Failed to compile FG regex from query '{queryFG}': {exc}"
		) from exc


if __name__ == "__main__":
    # testing
    db_host = os.environ["MONGO_URI"]
    db_port = int(os.environ["MONGO_PORT"])
    user_pd = os.environ["MONGO_PASSWORD"]
    user_name = os.environ["MONGO_USER"]
    client = MongoClient(host=db_host, port=db_port, username=user_name, password=user_pd, authSource='admin')
    client.list_database_names()  # check available databases
    print(client[os.environ["MONGO_DB"]].list_collection_names())
    # add arg parser
    # if the --use-old flag is provided, use the old function
    # add aregment for queryFG use the example one as default
    

    parser = argparse.ArgumentParser(description="Search proteins by FG query.")
    parser.add_argument("--use-old", action="store_true", help="Use the old function")
    parser.add_argument("--queryFG", type=str, default="*, (PF05658 | PF05662){1,20}, *, PF03895", help="FG query string")
    args = parser.parse_args()
    
    if args.use_old:
        res = search_protein_by_FG_old(
            queryFG=args.queryFG,
            DB_name=os.environ["MONGO_DB"],
            client=client,
        )
    else:
        res = search_protein_by_FG(
            queryFG=args.queryFG,
            DB_name=os.environ["MONGO_DB"],
            client=client,
        )
    print(res)

