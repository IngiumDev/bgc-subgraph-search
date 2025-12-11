from pymongo import MongoClient
import re
from collections import defaultdict

import numpy as np
import os
import sys


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
        host=str(db_host),
        port=int(db_port),
        username=str(user_name),
        password=str(user_pd),
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

        for doc in results:
            pfam_dict[p_id].add(doc["protein_id"])

    res = {k: len(v) for k, v in pfam_dict.items()}  # convert to pfam: count
    print(res)
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
    query_parts = [part.strip() for part in FG_query.split(',')]
    mandatory_domains = []

    pfam_pattern = re.compile(r'PF\d{4,6}') 

    for part in query_parts:
        # a) Exclude arbitrary wildcards.
        if part in ['.', '.*', '*']:
            continue
        
        # b) Exclude parts containing the 'OR' operator '|'.
        if '|' in part:
            continue
            
        # c) Check if the part contains a specific Pfam ID.
        match = pfam_pattern.search(part)
        if match:
            pfam_id = match.group(0)
            
            is_mandatory = True
            
            # Extract the part following the Pfam ID (the quantifier).
            quantifier_part = part[match.end():]
            
            if quantifier_part:
                if '?' in quantifier_part or '*' in quantifier_part:
                    is_mandatory = False
                    
                # Check for explicit quantifiers like {n,m}
                quantifier_match = re.search(r'\{(\d+),?\d*\}', quantifier_part)
                if quantifier_match:
                    # Check the minimum count 'n'. If n=0, it's not mandatory.
                    min_count = int(quantifier_match.group(1))
                    if min_count == 0:
                        is_mandatory = False 

            if is_mandatory:
                mandatory_domains.append(pfam_id)
        
    unique_mandatory_domains = list(set(mandatory_domains))

    return unique_mandatory_domains


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
    print("\n--- Running b: get_mandatory_domains ---")
    mandatory_list = get_mandatory_domains(FG_query)
    print("Mandatory Domains:", mandatory_list)
    # c)
    # FG_candidate = (
    #     "PF03895PF03895PF03895PF05658PF05662PF05662PF03895PF03895PF03895PF0023"
    # )
    # s = matches_fg_query(FG_query, FG_candidate)
    # print(s)
