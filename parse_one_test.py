import requests
import json
import sys # Keep sys for exit codes indicating success/failure

# --- Configuration ---
PDB_ID = "1FJE"
CHAIN_ID = "B" # Author Asymmetric ID

# --- Functions ---

def get_api_data(api_url, query_type="Generic"):
    """Fetches data from a specified RCSB API endpoint."""
    try:
        response = requests.get(api_url, timeout=30)
        response.raise_for_status()
        data = response.json()
        return data
    except requests.exceptions.HTTPError as http_err:
        # Optionally print a user-facing error, but for now just return None
        # print(f"Error fetching {query_type} data: {http_err}", file=sys.stderr)
        pass
    except requests.exceptions.RequestException as req_err:
        # print(f"Network error fetching {query_type} data: {req_err}", file=sys.stderr)
        pass
    except requests.exceptions.JSONDecodeError:
        # print(f"Error decoding JSON for {query_type} data.", file=sys.stderr)
        pass
    except Exception as e:
        # print(f"An unexpected error occurred fetching {query_type} data: {e}", file=sys.stderr)
        pass
    return None


def get_entity_id_and_auth_start(instance_data):
    """Extracts the entity_id and the starting author sequence number from instance data."""
    if not instance_data:
        return None, None

    entity_id = None
    auth_start_num = None

    try:
        container_ids = instance_data.get('rcsb_polymer_entity_instance_container_identifiers', {})
        entity_id = container_ids.get('entity_id')
        mapping = container_ids.get('auth_to_entity_poly_seq_mapping', [])

        if not entity_id:
            return None, None

        if not mapping:
            # Return entity_id but None for auth_start if mapping is missing
            return entity_id, None
        else:
            try:
                auth_start_str = mapping[0]
                auth_start_num = int(auth_start_str)
            except (IndexError, ValueError, TypeError):
                 # Return entity_id but None for auth_start if conversion fails
                return entity_id, None

    except Exception:
        # Catch any other unexpected errors during parsing
        return None, None

    return entity_id, auth_start_num


def extract_uniprot_id_and_uniprot_start(entity_data):
    """Extracts UniProt ID and the corresponding UniProt start number (for entity residue 1) from entity data."""
    if not entity_data:
        return None, None

    uniprot_id = None
    uniprot_start_num = None

    # 1. Extract UniProt ID
    try:
        identifiers = entity_data.get('rcsb_polymer_entity_container_identifiers', {}).get('reference_sequence_identifiers', [])
        if identifiers:
            for identifier in identifiers:
                if identifier.get('database_name') == 'UniProt':
                    uniprot_id = identifier.get('database_accession')
                    if uniprot_id:
                        break
    except Exception:
         # Ignore errors here, try getting alignment anyway
         pass

    # 2. Extract UniProt Start Number from Alignment (for entity_beg_seq_id == 1)
    try:
        entity_align_list = entity_data.get('rcsb_polymer_entity_align', [])
        if entity_align_list:
            found_start_num = False
            for alignment in entity_align_list:
                if alignment.get('reference_database_name') == 'UniProt':
                    regions = alignment.get('aligned_regions', [])
                    if regions:
                        for region in regions:
                            if region.get('entity_beg_seq_id') == 1:
                                ref_start_str = region.get('ref_beg_seq_id')
                                if ref_start_str is not None:
                                    try:
                                        uniprot_start_num = int(ref_start_str)
                                        found_start_num = True
                                        break
                                    except (ValueError, TypeError):
                                        # If conversion fails, keep uniprot_start_num as None
                                        uniprot_start_num = None
                                        found_start_num = False # Mark as not found
                                        break # Stop trying regions if conversion failed
                        if found_start_num:
                            break # Stop trying alignments if start found
            # If loop finishes without finding, uniprot_start_num remains None

    except Exception:
        # Ignore errors during alignment parsing, uniprot_start_num will remain None
        pass

    return uniprot_id, uniprot_start_num


# --- Main Execution ---
if __name__ == "__main__":
    final_uniprot_id = None
    final_shift = None
    entity_id = None
    auth_start = None
    uniprot_start = None

    # 1. Get instance data
    instance_api_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{PDB_ID}/{CHAIN_ID}"
    instance_data = get_api_data(instance_api_url, query_type="Instance")

    if instance_data:
        # 2. Extract entity_id and author start number from instance data
        entity_id, auth_start = get_entity_id_and_auth_start(instance_data)

    # 3. If entity_id found, get entity data
    if entity_id:
        entity_api_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{PDB_ID}/{entity_id}"
        entity_data = get_api_data(entity_api_url, query_type="Entity")

        if entity_data:
            # 4. Extract UniProt ID and UniProt start number from entity data
            final_uniprot_id, uniprot_start = extract_uniprot_id_and_uniprot_start(entity_data)

    # 5. Calculate Shift if both start numbers were found
    if auth_start is not None and uniprot_start is not None:
        try:
            final_shift = int(uniprot_start) - int(auth_start)
        except (ValueError, TypeError):
            final_shift = None # Mark shift as undetermined if calculation fails
    else:
        final_shift = None

    # --- Final Output ---
    print("\n" + "="*30)
    print(f"Results for PDB {PDB_ID}, Chain {CHAIN_ID} (using API):")
    print("-"*30)
    if final_uniprot_id:
        print(f"Associated UniProt ID: {final_uniprot_id}")
    else:
        print("Associated UniProt ID: Not Found via API")

    if final_shift is not None:
        print(f"Numbering Shift (UniProt - Auth): {final_shift}")
    else:
        print("Numbering Shift: Could Not Determine via API")
    print("="*30)

    # Exit with error code if data wasn't found
    if final_uniprot_id is None or final_shift is None:
        sys.exit(1)
    else:
        sys.exit(0)