import requests
import re
from collections import defaultdict
import sys
import json
import os

def parse_journal_file(journal_file_path):
    """
    Parse the journal file to create a mapping between short UniProt IDs and full UniProt IDs.
    
    Args:
        journal_file_path: Path to the journal.pcbi.1010859.s010.txt file
        
    Returns:
        A dictionary mapping PDB_chain ("1G2E_A") to full UniProt ID ("P26378_RRM2")
    """
    pdb_chain_to_full_uniprot = {}
    
    try:
        with open(journal_file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                # Parse the line (format: P26378_RRM2_1G2E_A_127-195_139-207_123-222_B)
                parts = line.split('_')
                if len(parts) < 4:
                    continue
                
                # Extract the full UniProt ID ("P26378_RRM2")
                full_uniprot_id = f"{parts[0]}_{parts[1]}"
                
                # Extract the PDB ID and chain
                pdb_id = parts[2]
                chain_id = parts[3]
                
                # Create a key for the PDB chain ("1G2E_A")
                pdb_chain_key = f"{pdb_id}_{chain_id}"
                
                # Store the mapping
                pdb_chain_to_full_uniprot[pdb_chain_key] = full_uniprot_id
        
        return pdb_chain_to_full_uniprot
    
    except Exception as e:
        print(f"Error parsing journal file: {e}")
        return {}

def get_uniprot_id(pdb_id, chain_id):
    """Get UniProt ID for a given PDB ID and chain using multiple methods."""
    # Method 1: RCSB GraphQL API
    url = "https://data.rcsb.org/graphql"
    query = """
    {
      entry(entry_id: "%s") {
        polymer_entities {
          entity_poly {
            pdbx_strand_id
            rcsb_entity_polymer_type
          }
          rcsb_polymer_entity_container_identifiers {
            auth_asym_ids
            uniprot_ids
          }
        }
      }
    }
    """ % pdb_id
    
    try:
        response = requests.post(url, json={'query': query})
        if response.status_code == 200:
            data = response.json()
            if "data" in data and "entry" in data["data"] and data["data"]["entry"] is not None:
                for entity in data["data"]["entry"]["polymer_entities"]:
                    identifiers = entity.get("rcsb_polymer_entity_container_identifiers", {})
                    asym_ids = identifiers.get("auth_asym_ids", [])
                    
                    if chain_id in asym_ids and "uniprot_ids" in identifiers and identifiers["uniprot_ids"]:
                        return identifiers["uniprot_ids"][0]
    except Exception as e:
        print(f"Error with GraphQL query for {pdb_id}: {e}")
    
    # Method 2: PDBe SIFTS API
    try:
        url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            if pdb_id in data:
                mappings = data[pdb_id].get("UniProt", {})
                for uniprot_id, mapping_data in mappings.items():
                    for segment in mapping_data["mappings"]:
                        if segment["chain_id"] == chain_id:
                            return uniprot_id
    except Exception as e:
        print(f"Error with PDBe SIFTS API for {pdb_id}: {e}")
    
    # Method 3: Legacy XML API
    try:
        url = f"https://www.rcsb.org/pdb/rest/describeMol?structureId={pdb_id}"
        response = requests.get(url)
        if response.status_code == 200:
            try:
                import xml.etree.ElementTree as ET
                root = ET.fromstring(response.content)
                for polymer in root.findall(".//polymer"):
                    chains = polymer.get("chains", "")
                    if chain_id in chains.split(","):
                        for db_ref in polymer.findall(".//dbRef"):
                            if db_ref.get("dbSource") == "UniProt":
                                return db_ref.get("dbAccession", "Unknown")
            except Exception as e:
                print(f"Error parsing XML for {pdb_id}: {e}")
    except Exception as e:
        print(f"Error with legacy XML API for {pdb_id}: {e}")
    
    return "Unknown"

def get_uniprot_sequence(pdb_id, chain_id):
    """
    Get the UniProt sequence for a specific chain from RCSB PDB.
    This improved function specifically targets the sequence shown highlighted in blue on the website.
    """
    # First try to use the RCSB API to get entity information and UniProt mapping
    try:
        # Get entity ID for the chain
        url = f"https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{pdb_id}/{chain_id}"
        print(f"Fetching entity data from {url}")
        
        response = requests.get(url)
        if response.status_code != 200:
            print(f"Warning: Could not fetch entity data. Status code: {response.status_code}")
            entity_id = None
        else:
            data = response.json()
            entity_id = data.get('rcsb_polymer_entity_instance', {}).get('entity_id')
            print(f"Found entity_id: {entity_id}")
        
        if entity_id:
            # Get sequence from entity
            entity_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
            entity_response = requests.get(entity_url)
            
            if entity_response.status_code == 200:
                entity_data = entity_response.json()
                uniprot_ids = entity_data.get('rcsb_polymer_entity', {}).get('rcsb_uniprot_accession', [])
                
                if uniprot_ids:
                    uniprot_id = uniprot_ids[0]
                    print(f"Found UniProt ID: {uniprot_id}")
                    
                    # Get the aligned sequence (this is specifically what's shown in the web UI)
                    align_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity_align/{pdb_id}/{entity_id}"
                    align_response = requests.get(align_url)
                    
                    if align_response.status_code == 200:
                        align_data = align_response.json()
                        alignments = align_data.get('polymer_entity_align', [])
                        
                        for alignment in alignments:
                            ref_db_name = alignment.get('ref_db_name')
                            ref_db_id = alignment.get('ref_db_id')
                            
                            if ref_db_name == 'UniProt' and ref_db_id == uniprot_id:
                                # This is the UniProt alignment we want
                                aligned_regions = alignment.get('aligned_regions', [])
                                
                                if aligned_regions:
                                    # There might be multiple regions - get all of them
                                    segments = []
                                    
                                    for region in aligned_regions:
                                        ref_seq = region.get('ref_sequence')
                                        if ref_seq:
                                            segments.append(ref_seq)
                                    
                                    if segments:
                                        return ''.join(segments)
    except Exception as e:
        print(f"Error using API approach: {str(e)}")

    # Try to get the FASTA file which has UniProt info
    print("API approach failed, trying to scrape the FASTA file...")
    try:
        fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/display"
        fasta_response = requests.get(fasta_url)
        
        if fasta_response.status_code == 200:
            fasta_text = fasta_response.text
            
            # Extract UniProt ID from FASTA
            uniprot_id = None
            for line in fasta_text.split('\n'):
                if line.startswith('>') and f"|Chain {chain_id}|" in line:
                    # This is the header for our chain
                    print(f"Found FASTA header: {line}")
                    
                    # Check for UniProt mention in comment
                    match = re.search(r'UniProt: ([A-Z0-9]+)', line)
                    if match:
                        uniprot_id = match.group(1)
                        print(f"Found UniProt ID from FASTA: {uniprot_id}")
            
            if uniprot_id is None:
                # Try to get UniProt ID using other methods
                uniprot_id = get_uniprot_id(pdb_id, chain_id)
                
            if uniprot_id and uniprot_id != "Unknown":
                # Use SIFTS mapping to get the sequence
                sifts_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
                sifts_response = requests.get(sifts_url)
                
                if sifts_response.status_code == 200:
                    sifts_data = sifts_response.json()
                    pdb_entry = sifts_data.get(pdb_id.lower(), {})
                    
                    if pdb_entry and 'UniProt' in pdb_entry and uniprot_id in pdb_entry['UniProt']:
                        mappings = pdb_entry['UniProt'][uniprot_id]['mappings']
                        
                        for mapping in mappings:
                            if mapping['chain_id'] == chain_id:
                                unp_start = mapping.get('unp_start')
                                unp_end = mapping.get('unp_end')
                                
                                if unp_start is not None and unp_end is not None:
                                    # Get the UniProt sequence
                                    uniprot_url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
                                    uniprot_response = requests.get(uniprot_url)
                                    
                                    if uniprot_response.status_code == 200:
                                        fasta_lines = uniprot_response.text.strip().split('\n')
                                        # Skip the header line
                                        uniprot_sequence = ''.join(fasta_lines[1:])
                                        
                                        # Extract the mapped region
                                        return uniprot_sequence[unp_start-1:unp_end]
    except Exception as e:
        print(f"Error with FASTA approach: {str(e)}")
    
    # Fallback to old method as last resort
    return get_domain_sequence(pdb_id, chain_id)

def get_domain_sequence(pdb_id, chain_id):
    """Get the domain sequence (RRM) as shown in blue on the website."""
    try:
        url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/polymer_coverage/{pdb_id}"
        response = requests.get(url)
        
        if response.status_code == 200:
            data = response.json()
            if pdb_id in data:
                for entity in data[pdb_id]:
                    for chain in entity.get("chains", []):
                        if chain["chain_id"] == chain_id:
                            observed = chain.get("observed", [])
                            if observed:
                                seq_url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_id}"
                                seq_response = requests.get(seq_url)
                                
                                if seq_response.status_code == 200:
                                    seq_data = seq_response.json()
                                    
                                    for mol in seq_data.get(pdb_id, []):
                                        for chain_data in mol.get("chains", []):
                                            if chain_data["chain_id"] == chain_id:
                                                seq = mol.get("sequence", "")
                                                
                                                domain_seq = ""
                                                for obs_range in observed:
                                                    start = obs_range["start"]["residue_number"] - 1  # Convert to 0-indexed
                                                    end = obs_range["end"]["residue_number"]
                                                    domain_seq += seq[start:end]
                                                
                                                return domain_seq
    except Exception as e:
        print(f"Error getting domain sequence from PDBe: {e}")
    
    try:
        url = f"https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{pdb_id}/{chain_id}"
        response = requests.get(url)
        
        if response.status_code == 200:
            data = response.json()
            if "rcsb_polymer_instance_feature" in data:
                for feature in data["rcsb_polymer_instance_feature"]:
                    if feature.get("type") == "RRM_DOM":
                        start = feature.get("feature_positions", [{}])[0].get("beg_seq_id", 0)
                        end = feature.get("feature_positions", [{}])[0].get("end_seq_id", 0)
                        
                        if start > 0 and end > 0:
                            entity_id = data.get("rcsb_polymer_entity_instance_container_identifiers", {}).get("entity_id")
                            
                            if entity_id:
                                entity_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
                                entity_response = requests.get(entity_url)
                                
                                if entity_response.status_code == 200:
                                    entity_data = entity_response.json()
                                    seq = entity_data.get("entity_poly", {}).get("pdbx_seq_one_letter_code", "")
                                    
                                    if seq:
                                        return seq[start-1:end]
    except Exception as e:
        print(f"Error getting domain sequence from RCSB: {e}")
    
    try:
        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        pdb_response = requests.get(pdb_url)
        
        if pdb_response.status_code == 200:
            pdb_text = pdb_response.text
            
            # Get residue ranges from the PDB file
            residues = {}
            for line in pdb_text.splitlines():
                if line.startswith("ATOM"):
                    curr_chain = line[21].strip()
                    if curr_chain != chain_id:
                        continue
                    
                    try:
                        res_num = int(line[22:26].strip())
                        res_name = line[17:20].strip()
                        
                        # Map three-letter code to one-letter code
                        three_to_one = {
                            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
                            "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
                            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
                            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
                        }
                        
                        one_letter = three_to_one.get(res_name, "X")
                        residues[res_num] = one_letter
                    except:
                        continue
            
            if residues:
                # Create sequence from residues
                min_res = min(residues.keys())
                max_res = max(residues.keys())
                sequence = ''.join([residues.get(i, "") for i in range(min_res, max_res + 1)])
                return sequence
    except Exception as e:
        print(f"Error extracting sequence from PDB file: {e}")
    
    # Try to scrape the actual sequence from RCSB website
    try:
        url = f"https://www.rcsb.org/structure/{pdb_id}"
        response = requests.get(url)
        
        if response.status_code == 200:
            html = response.text
            
            sequence_match = re.search(r'sequenceInfo\s*=\s*(\[.*?\]);', html, re.DOTALL)
            if sequence_match:
                sequence_data = sequence_match.group(1)
                try:
                    import json
                    sequences = json.loads(sequence_data)
                    for seq in sequences:
                        if seq.get("chainId") == chain_id:
                            return seq.get("sequence", "")
                except:
                    pass
    except Exception as e:
        print(f"Error scraping sequence from RCSB website: {e}")
    
    return ""

def extract_residue_number_range(pdb_id, chain_id):
    """Extract the residue number range directly from the structure in PDB file."""
    try:
        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(pdb_url)
        
        if response.status_code != 200:
            return None, None
        
        pdb_text = response.text
        
        # Find all residue numbers in the specified chain
        residue_numbers = []
        for line in pdb_text.splitlines():
            if line.startswith("ATOM"):
                current_chain = line[21].strip()
                if current_chain != chain_id:
                    continue
                
                try:
                    res_num = int(line[22:26].strip())
                    residue_numbers.append(res_num)
                except ValueError:
                    continue
        
        if residue_numbers:
            return min(residue_numbers), max(residue_numbers)
        
        return None, None
    except Exception as e:
        print(f"Error extracting residue range: {e}")
        return None, None

def get_api_data(api_url, query_type="Generic"):
    """Fetches data from a specified RCSB API endpoint."""
    try:
        response = requests.get(api_url, timeout=30)
        response.raise_for_status()
        data = response.json()
        return data
    except requests.exceptions.HTTPError as http_err:
        pass
    except requests.exceptions.RequestException as req_err:
        pass
    except requests.exceptions.JSONDecodeError:
        pass
    except Exception as e:
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
            return entity_id, None
        else:
            try:
                auth_start_str = mapping[0]
                auth_start_num = int(auth_start_str)
            except (IndexError, ValueError, TypeError):
                return entity_id, None

    except Exception:
        return None, None

    return entity_id, auth_start_num


def extract_uniprot_id_and_uniprot_start(entity_data):
    """Extracts UniProt ID and the corresponding UniProt start number (for entity residue 1) from entity data."""
    if not entity_data:
        return None, None

    uniprot_id = None
    uniprot_start_num = None

    try:
        identifiers = entity_data.get('rcsb_polymer_entity_container_identifiers', {}).get('reference_sequence_identifiers', [])
        if identifiers:
            for identifier in identifiers:
                if identifier.get('database_name') == 'UniProt':
                    uniprot_id = identifier.get('database_accession')
                    if uniprot_id:
                        break
    except Exception:
        pass

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
                                        uniprot_start_num = None
                                        found_start_num = False
                                        break
                        if found_start_num:
                            break

    except Exception:
        pass

    return uniprot_id, uniprot_start_num


def calculate_uniprot_shift(pdb_id, chain_id):
    """Calculate the numbering shift between UniProt and Author residue numbering."""
    final_shift = None
    
    # 1. Get instance data
    instance_api_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{pdb_id}/{chain_id}"
    instance_data = get_api_data(instance_api_url, query_type="Instance")

    if instance_data:
        # 2. Extract entity_id and author start number from instance data
        entity_id, auth_start = get_entity_id_and_auth_start(instance_data)

        # 3. If entity_id found, get entity data
        if entity_id:
            entity_api_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
            entity_data = get_api_data(entity_api_url, query_type="Entity")

            if entity_data:
                # 4. Extract UniProt ID and UniProt start number from entity data
                _, uniprot_start = extract_uniprot_id_and_uniprot_start(entity_data)

                # 5. Calculate Shift if both start numbers were found
                if auth_start is not None and uniprot_start is not None:
                    try:
                        final_shift = int(uniprot_start) - int(auth_start)
                    except (ValueError, TypeError):
                        final_shift = None

    return final_shift

def process_data_file(input_file, protein_cds_file, output_file="extracted_data.txt", journal_file="journal.pcbi.1010859.s010.txt"):
    """Process the input data file and protein/RNA sequences to extract and write results.
    
    Args:
        input_file: File with PDB chains and connection info
        protein_cds_file: FASTA file with protein and RNA sequences
        output_file: Output file to write results
        journal_file: Path to the journal file with full UniProt IDs
    """
    connections = defaultdict(list)
    
    protein_cache = {}
    
    print(f"Parsing journal file {journal_file} for full UniProt IDs...")
    pdb_chain_to_full_uniprot = parse_journal_file(journal_file)
    print(f"Found {len(pdb_chain_to_full_uniprot)} PDB chain to full UniProt ID mappings")
    
    short_to_full_uniprot = {}
    for pdb_chain, full_uniprot in pdb_chain_to_full_uniprot.items():
        short_uniprot = full_uniprot.split('_')[0]
        short_to_full_uniprot[short_uniprot] = full_uniprot
    
    print(f"Created mapping for {len(short_to_full_uniprot)} UniProt IDs")
    
    # Read input file with PDB chains
    try:
        with open(input_file, 'r') as f:
            input_lines = [line.strip() for line in f if line.strip()]
        
        print(f"Processing {len(input_lines)} lines from {input_file}...")
        
        valid_lines = [line for line in input_lines if ',' in line]
        if not valid_lines:
            print(f"Error: No valid input lines found in {input_file}. Each line should be in format: PDB_CHAIN,POSITION_AA_...")
            return
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        return
    except Exception as e:
        print(f"Error reading input file: {e}")
        return
    
    # Read protein_cds.fasta file to map proteins to RNA sequences
    protein_to_rna = {}
    try:
        with open(protein_cds_file, 'r') as f:
            fasta_content = f.read()
        
        print(f"Reading protein and RNA sequences from {protein_cds_file}...")
        
        fasta_entries = fasta_content.strip().split('>')
        
        for entry in fasta_entries:
            if not entry:
                continue
            
            entry_lines = entry.strip().split('\n')
            header = entry_lines[0].strip()
            sequence = ''.join(entry_lines[1:]).strip()
            
            if not header or not sequence:
                continue
            
            if header.startswith('NM_') or header.startswith('XM_'):
                continue
                
            protein_id = header
            protein_to_rna[protein_id] = {"protein": sequence, "rna": None}
        
        protein_rna_pairs = []
        
        for entry in fasta_entries:
            if not entry:
                continue
            
            entry_lines = entry.strip().split('\n')
            header = entry_lines[0].strip()
            sequence = ''.join(entry_lines[1:]).strip()
            
            if not header or not sequence:
                continue
            
            if header.startswith('NM_') or header.startswith('XM_'):
                protein_rna_pairs.append({"rna_id": header, "rna_seq": sequence})
        
        if len(protein_rna_pairs) == len(protein_to_rna):
            for i, (protein_id, data) in enumerate(protein_to_rna.items()):
                if i < len(protein_rna_pairs):
                    data["rna"] = protein_rna_pairs[i]["rna_seq"]
        else:
            for rna_data in protein_rna_pairs:
                for protein_id, data in protein_to_rna.items():
                    if data["rna"] is None:
                        data["rna"] = rna_data["rna_seq"]
                        break
        
        print(f"Found {len(protein_to_rna)} protein sequences in FASTA file")
    except FileNotFoundError:
        print(f"Error: Protein CDS file '{protein_cds_file}' not found.")
        return
    except Exception as e:
        print(f"Error reading protein CDS file: {e}")
        return
    
    unique_pdb_chains = set()
    for line in input_lines:
        parts = line.split(',')
        if len(parts) >= 1:
            unique_pdb_chains.add(parts[0])

    total_chains = len(unique_pdb_chains)
    processed_chains = 0
    
    for line in input_lines:
        parts = line.split(',')
        if len(parts) != 2:
            continue
        
        pdb_chain = parts[0]
        connection_part = parts[1]
        
        pdb_match = re.match(r'(.+)_(.)', pdb_chain)
        if not pdb_match:
            print(f"Invalid PDB ID format: {pdb_chain}")
            continue
        
        pdb_id = pdb_match.group(1)
        chain_id = pdb_match.group(2)
        
        connection_match = re.match(r'(\d+)_([A-Z])_', connection_part)
        if not connection_match:
            print(f"Invalid connection format: {connection_part}")
            continue
        
        position = int(connection_match.group(1))
        expected_aa = connection_match.group(2)
        
        connections[pdb_chain].append(position)
        
        # Get protein data if not already cached
        if pdb_chain not in protein_cache:
            processed_chains += 1
            progress_percent = (processed_chains / total_chains) * 100
            print(f"\rProgress: {progress_percent:.1f}% - Processing {pdb_chain} ({processed_chains}/{total_chains})", end="", flush=True)
            
            uniprot_id = get_uniprot_id(pdb_id, chain_id)
            
            # Use the new improved method to get the UniProt sequence
            sequence = get_uniprot_sequence(pdb_id, chain_id)
            
            # find the corresponding nucleotide sequence
            nucleotide_seq = ""
            if uniprot_id != "Unknown" and sequence and len(sequence) > 0:
                if uniprot_id in protein_to_rna:
                    full_protein = protein_to_rna[uniprot_id]["protein"]
                    full_rna = protein_to_rna[uniprot_id]["rna"]
                    
                    if full_protein and full_rna:
                        subseq_pos = full_protein.find(sequence)
                        
                        if subseq_pos != -1:
                            rna_start = subseq_pos * 3
                            rna_end = (subseq_pos + len(sequence)) * 3
                            
                            if rna_start < len(full_rna) and rna_end <= len(full_rna):
                                nucleotide_seq = full_rna[rna_start:rna_end]
                        else:
                            for window_size in range(min(20, len(sequence)), 5, -1):
                                for i in range(len(sequence) - window_size + 1):
                                    subseq = sequence[i:i+window_size]
                                    subseq_pos = full_protein.find(subseq)
                                    
                                    if subseq_pos != -1:
                                        offset = i
                                        full_start = subseq_pos - offset
                                        
                                        if full_start >= 0:
                                            rna_start = full_start * 3
                                            rna_end = (full_start + len(sequence)) * 3
                                            
                                            if rna_start < len(full_rna) and rna_end <= len(full_rna):
                                                nucleotide_seq = full_rna[rna_start:rna_end]
                                                break
                                
                                if nucleotide_seq:
                                    break
            
            protein_cache[pdb_chain] = (sequence, uniprot_id, nucleotide_seq)
    
    print()
    
    # Write results to output file
    try:
        with open(output_file, 'w') as f:
            for pdb_chain in sorted(connections.keys()):
                sequence, uniprot_id, nucleotide_seq = protein_cache.get(pdb_chain, ("", "Unknown", ""))
                
                pdb_match = re.match(r'(.+)_(.)', pdb_chain)
                if pdb_match:
                    pdb_id = pdb_match.group(1)
                    chain_id = pdb_match.group(2)
                    pdb_chain_key = f"{pdb_id}_{chain_id}"
                    
                    # Calculate the UniProt shift
                    uniprot_shift = calculate_uniprot_shift(pdb_id, chain_id)
                    shift_info = f"Numbering Shift (UniProt - Auth): {uniprot_shift}" if uniprot_shift is not None else "Numbering Shift: Could Not Determine"
                    
                    if pdb_chain_key in pdb_chain_to_full_uniprot:
                        full_uniprot_id = pdb_chain_to_full_uniprot[pdb_chain_key]
                        uniprot_id = full_uniprot_id
                    elif uniprot_id in short_to_full_uniprot:
                        full_uniprot_id = short_to_full_uniprot[uniprot_id]
                        uniprot_id = full_uniprot_id
                
                unique_connections = sorted(set(connections[pdb_chain]))
                
                f.write(f">>{pdb_chain}\n")
                f.write(f"UniProt ID: {uniprot_id}\n")
                f.write(f"Protein Sequence: {sequence}\n")
                f.write(f"Nucleotide Sequence: {nucleotide_seq}\n")
                f.write(f"{shift_info}\n")
                f.write(f"List of connections: {', '.join(map(str, unique_connections))}\n\n")
        
        print(f"Data processing complete. Results written to {output_file}")
    except Exception as e:
        print(f"Error writing to output file: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python extract_data.py journal.pcbi.1010859.s011.txt proten_cds.fasta [output_file.txt] [journal_file.txt]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    protein_cds_file = sys.argv[2]
    
    # Optionally allow specifying the output file and journal file
    output_file = "extracted_data.txt"
    journal_file = "journal.pcbi.1010859.s010.txt"
    
    if len(sys.argv) > 3:
        output_file = sys.argv[3]
    
    if len(sys.argv) > 4:
        journal_file = sys.argv[4]
    
    process_data_file(input_file, protein_cds_file, output_file, journal_file)