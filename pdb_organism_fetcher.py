import requests
import time
import json
import argparse

# Configuration
INPUT_FILE = "extracted_data.txt"
RCSB_ENTRY_SUMMARY_API_URL_TEMPLATE = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_ENTITY_API_URL_TEMPLATE = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"

REQUEST_TIMEOUT = 15
API_CALL_DELAY = 0.25 
PRINT_FULL_API_RESPONSE = False # Controlled by command line arg

def _parse_organisms_from_entity_data_obj(entity_data_obj):
    """Helper to parse organism names from a JSON object representing a single entity's data."""
    found_organisms_set = set()
    entity_had_source_info_block = False

    source_organisms_list = entity_data_obj.get('rcsb_entity_source_organism')
    if source_organisms_list:
        entity_had_source_info_block = True
        for org_entry in source_organisms_list:
            if org_entry:
                ncbi_name = org_entry.get('ncbi_scientific_name')
                if ncbi_name:
                    found_organisms_set.add(ncbi_name)
                elif not ncbi_name: 
                    sci_name_fallback = org_entry.get('scientific_name')
                    if sci_name_fallback:
                        found_organisms_set.add(sci_name_fallback)
    
    entity_src_gen_list = entity_data_obj.get('entity_src_gen')
    if entity_src_gen_list:
        entity_had_source_info_block = True
        for src_gen_entry in entity_src_gen_list:
            if src_gen_entry:
                name_options = [
                    src_gen_entry.get('pdbx_gene_src_scientific_name'),
                    src_gen_entry.get('gene_src_scientific_name'), 
                    src_gen_entry.get('gene_src_common_name')
                ]
                for name_opt in name_options:
                    if name_opt:
                        found_organisms_set.add(name_opt)
                        break 
            
    if not entity_had_source_info_block and 'struct_ref' in entity_data_obj:
        struct_ref_list = entity_data_obj.get('struct_ref')
        if struct_ref_list:
            for ref in struct_ref_list:
                if ref.get('db_name') == 'UNP' and ref.get('pdbx_organism_scientific'):
                    found_organisms_set.add(ref.get('pdbx_organism_scientific'))
    return found_organisms_set

def get_organisms_for_pdb_entry(pdb_id_upper):
    """
    Fetches all unique organism/source names for a given PDB entry.
    Uses Strategy 1 (entry summary) and falls back to Strategy 2 (iterating all entities) if needed.
    Returns a string of comma-separated organism names or an error/status message.
    """
    entry_summary_url = RCSB_ENTRY_SUMMARY_API_URL_TEMPLATE.format(pdb_id=pdb_id_upper)
    all_entry_organisms_s1 = set()
    strategy1_response_data = None
    strategy1_failed_specifically = False
    strategy1_http_error_status = None
    s1_failure_reason = "S1: Unknown reason"


    try: # Strategy 1
        response_s1 = requests.get(entry_summary_url, timeout=REQUEST_TIMEOUT)
        response_s1.raise_for_status()
        strategy1_response_data = response_s1.json()

        if PRINT_FULL_API_RESPONSE:
            print(f"\n--- Full API Response (Strategy 1: Entry Summary for {pdb_id_upper}) ---")
            print(json.dumps(strategy1_response_data, indent=2))

        polymer_entities_s1 = strategy1_response_data.get('polymer_entities')
        if polymer_entities_s1:
            for entity in polymer_entities_s1:
                orgs_from_entity = _parse_organisms_from_entity_data_obj(entity)
                all_entry_organisms_s1.update(orgs_from_entity)
            if all_entry_organisms_s1: return ", ".join(sorted(list(all_entry_organisms_s1)))
            s1_failure_reason = "S1: 'polymer_entities' present but no orgs found"
            strategy1_failed_specifically = True
        else: 
            top_level_orgs_s1 = _parse_organisms_from_entity_data_obj(strategy1_response_data)
            if top_level_orgs_s1: return ", ".join(sorted(list(top_level_orgs_s1)))
            s1_failure_reason = "S1: No 'polymer_entities' key and no top-level orgs"
            strategy1_failed_specifically = True
    except requests.exceptions.HTTPError as http_err:
        strategy1_http_error_status = http_err.response.status_code if http_err.response else 'N/A'
        s1_failure_reason = f"S1: HTTP error {strategy1_http_error_status}"
        if strategy1_http_error_status == 404: return f"PDB ID {pdb_id_upper} not found (404)"
    except Exception as e:
        s1_failure_reason = f"S1: Exception ({str(e)[:30]})"
        if PRINT_FULL_API_RESPONSE: print(f"DEBUG: S1 Exception: {e}")

    if PRINT_FULL_API_RESPONSE or args.test_entry : # Print S1 failure only in debug/test modes
        print(f"  Strategy 1 for {pdb_id_upper} inconclusive ('{s1_failure_reason}'). Trying Strategy 2 (fallback)...")

    # --- Strategy 2: Iterate all polymer entities for the PDB ID ---
    all_entry_organisms_s2 = set()
    polymer_entity_ids_list = []
    
    if strategy1_response_data:
        container_ids = strategy1_response_data.get('rcsb_entry_container_identifiers')
        if container_ids: polymer_entity_ids_list = container_ids.get('polymer_entity_ids', [])
    
    if not polymer_entity_ids_list and strategy1_http_error_status != 404:
        if PRINT_FULL_API_RESPONSE: print(f"  S2: Re-fetching {entry_summary_url} for entity ID list...")
        try:
            response_s2_entry = requests.get(entry_summary_url, timeout=REQUEST_TIMEOUT)
            response_s2_entry.raise_for_status()
            entry_data_for_s2 = response_s2_entry.json()
            if PRINT_FULL_API_RESPONSE:
                print(f"\n--- Full API Response (S2 re-fetch Entry Summary for {pdb_id_upper}) ---")
                print(json.dumps(entry_data_for_s2, indent=2))
            container_ids = entry_data_for_s2.get('rcsb_entry_container_identifiers')
            if container_ids: polymer_entity_ids_list = container_ids.get('polymer_entity_ids', [])
        except Exception as e:
            return f"{s1_failure_reason}, S2: Failed to re-fetch entry for entity list ({str(e)[:50]})"

    if not polymer_entity_ids_list:
        return f"{s1_failure_reason}, S2: No polymer_entity_ids found for PDB {pdb_id_upper}"

    if PRINT_FULL_API_RESPONSE: 
        print(f"  S2: Found {len(polymer_entity_ids_list)} polymer entity IDs: {polymer_entity_ids_list}. Querying each...")
    
    for entity_id in polymer_entity_ids_list:
        entity_api_url = RCSB_ENTITY_API_URL_TEMPLATE.format(pdb_id=pdb_id_upper, entity_id=entity_id)
        try:
            response_entity_s2 = requests.get(entity_api_url, timeout=REQUEST_TIMEOUT)
            response_entity_s2.raise_for_status()
            entity_data_s2 = response_entity_s2.json()
            if PRINT_FULL_API_RESPONSE:
                print(f"\n--- Full API Response (S2: Entity {entity_id} for {pdb_id_upper}) ---")
                print(json.dumps(entity_data_s2, indent=2))
            
            orgs_from_this_entity = _parse_organisms_from_entity_data_obj(entity_data_s2)
            all_entry_organisms_s2.update(orgs_from_this_entity)
            time.sleep(API_CALL_DELAY)
        except Exception as e_s2:
            if PRINT_FULL_API_RESPONSE:
                status_s2 = e_s2.response.status_code if hasattr(e_s2, 'response') and e_s2.response else 'N/A'
                print(f"    S2: Error for entity {entity_id} (status {status_s2}): {str(e_s2)[:50]}. Skipping.")

    if all_entry_organisms_s2:
        return ", ".join(sorted(list(all_entry_organisms_s2)))
    else:
        return f"{s1_failure_reason}, S2: No organisms found after checking all entities."


def main():
    parser = argparse.ArgumentParser(description="Fetch organism data for PDB entries.")
    parser.add_argument('--test_entry', type=str, default=None,
                        help="Process only this specific entry (e.g., 1A9N_B or 1A9N) from the input file.")
    parser.add_argument('--print_api_response', action='store_true', help="Print full API JSON and extra debug info.")
    global args # Make args globally accessible for PRINT_FULL_API_RESPONSE checks in helper
    args = parser.parse_args()


    global PRINT_FULL_API_RESPONSE
    PRINT_FULL_API_RESPONSE = args.print_api_response
    if PRINT_FULL_API_RESPONSE: print("INFO: Verbose logging and full API responses will be printed.\n")
    if args.test_entry: print(f"INFO: --test_entry '{args.test_entry}' active.\n")

    results_log = []
    
    target_pdb_id_from_arg = None
    if args.test_entry:
        try:
            target_pdb_id_from_arg = args.test_entry.split('_', 1)[0].upper() if '_' in args.test_entry else args.test_entry.upper()
        except: #pylint: disable=bare-except
            print(f"Warning: Could not parse PDB ID from --test_entry '{args.test_entry}'.")
            target_pdb_id_from_arg = args.test_entry.upper() # Try using it as is

    if args.test_entry:
        if not target_pdb_id_from_arg:
             print(f"ERROR: Invalid --test_entry format or empty PDB ID: '{args.test_entry}'.")
             return
        
        input_line_for_logging = args.test_entry # Use the argument itself for logging
        
        organism_str = get_organisms_for_pdb_entry(target_pdb_id_from_arg)
        is_homo_sapiens = "homo sapiens" in organism_str.lower() if isinstance(organism_str, str) else False
        
        marker = " <-- HOMO SAPIENS" if is_homo_sapiens else ""
        print(f"{input_line_for_logging}: {organism_str}{marker}")
        # No extensive summary needed for single test_entry run
        return # Exit after processing the single test entry

    # Normal mode (process all entries from file)
    processed_pdb_ids_organisms = {} 
    unique_pdb_ids_processed_count = 0
    try:
        with open(INPUT_FILE, 'r') as f:
            for line_content in f:
                line_content = line_content.strip()
                if not line_content.startswith(">>"): continue

                pdb_chain_original_from_file = line_content[2:].strip()
                pdb_id_val = "ParseError"
                
                try:
                    pdb_id_val = pdb_chain_original_from_file.split('_', 1)[0] if '_' in pdb_chain_original_from_file else pdb_chain_original_from_file
                except ValueError: # Should not happen with split if '_' is checked
                    pass # pdb_id_val remains ParseError
                
                if pdb_id_val == "ParseError":
                    print(f"{pdb_chain_original_from_file}: Parse Error (PDB_Chain format)")
                    results_log.append({"input": pdb_chain_original_from_file, "org": "Parse Error", "hs": False})
                    continue
                
                pdb_id_upper = pdb_id_val.upper()
                organism_str = "Not Fetched Yet"

                if pdb_id_upper in processed_pdb_ids_organisms:
                    organism_str = processed_pdb_ids_organisms[pdb_id_upper]
                    if PRINT_FULL_API_RESPONSE: # Only print this if debugging
                         print(f"  PDB ID {pdb_id_upper} (from {pdb_chain_original_from_file}) already processed. Using cached: {organism_str}")
                else:
                    if PRINT_FULL_API_RESPONSE:
                         print(f"Fetching organism(s) for PDB Entry: {pdb_id_upper} (from input: {pdb_chain_original_from_file})...")
                    organism_str = get_organisms_for_pdb_entry(pdb_id_upper)
                    processed_pdb_ids_organisms[pdb_id_upper] = organism_str
                    unique_pdb_ids_processed_count +=1
                
                is_homo_sapiens = "homo sapiens" in organism_str.lower() if isinstance(organism_str, str) else False
                results_log.append({"input": pdb_chain_original_from_file, "org": organism_str, "hs": is_homo_sapiens})
                
                marker = " <-- HOMO SAPIENS" if is_homo_sapiens else ""
                print(f"{pdb_chain_original_from_file}: {organism_str}{marker}")

    except FileNotFoundError: print(f"ERROR: File '{INPUT_FILE}' not found."); return
    except Exception as e: print(f"ERROR (batch mode): {e}"); import traceback; traceback.print_exc(); return

    # --- Summary for normal mode ---
    if not args.test_entry:
        print("\n--- Summary ---")
        homo_sapiens_pdb_ids = set()
        if not results_log:
            print("No PDB entries were processed from the file.")
        else:
            # Create a summary based on unique PDB IDs from the cache
            for pdb_key_sorted in sorted(processed_pdb_ids_organisms.keys()):
                org_data = processed_pdb_ids_organisms[pdb_key_sorted]
                is_hs = "homo sapiens" in org_data.lower() if isinstance(org_data, str) else False
                # This part of summary can be removed if too verbose, as results are printed per line
                # marker = " <-- HOMO SAPIENS" if is_hs else ""
                # print(f"PDB Entry {pdb_key_sorted}: {org_data}{marker}") 
                if is_hs:
                    homo_sapiens_pdb_ids.add(pdb_key_sorted)
        
        print(f"\nTotal unique PDB IDs processed from file: {unique_pdb_ids_processed_count}")
        print(f"Unique PDB IDs where Homo sapiens was found: {len(homo_sapiens_pdb_ids)}")

if __name__ == "__main__":
    main()