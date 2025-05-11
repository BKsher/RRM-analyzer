import requests
from bs4 import BeautifulSoup
import argparse
import re
import json

def get_uniprot_sequence(pdb_id, chain):
    """
    Extract the UniProt sequence for a specific chain from RCSB PDB.
    This function tries several approaches to get the sequence.
    """
    
    # First try to use the RCSB API to get entity information and UniProt mapping
    try:
        # Get entity ID for the chain
        url = f"https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{pdb_id}/{chain}"
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

    # If the API approach fails, try screen scraping
    print("API approach failed, trying to scrape the FASTA file...")
    
    try:
        # Try to get the FASTA file which has UniProt info
        fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/display"
        fasta_response = requests.get(fasta_url)
        
        if fasta_response.status_code == 200:
            fasta_text = fasta_response.text
            
            # Extract UniProt ID from FASTA
            uniprot_id = None
            for line in fasta_text.split('\n'):
                if line.startswith('>') and f"|Chain {chain}|" in line:
                    # This is the header for our chain
                    print(f"Found FASTA header: {line}")
                    
                    # Check for UniProt mention in comment
                    match = re.search(r'UniProt: ([A-Z0-9]+)', line)
                    if match:
                        uniprot_id = match.group(1)
                        print(f"Found UniProt ID from FASTA: {uniprot_id}")
                    else:
                        # Try to get UniProt ID from chain info
                        entity_url = f"https://data.rcsb.org/graphql?query={{polymer_entity_instances(entry_id:%22{pdb_id}%22){{rcsb_id,rcsb_polymer_entity_instance_container_identifiers{{auth_asym_id}},polymer_entity{{rcsb_polymer_entity_container_identifiers{{uniprot_ids}}}}}}}}"
                        entity_response = requests.get(entity_url)
                        
                        if entity_response.status_code == 200:
                            graphql_data = entity_response.json()
                            instances = graphql_data.get('data', {}).get('polymer_entity_instances', [])
                            
                            for instance in instances:
                                asym_id = instance.get('rcsb_polymer_entity_instance_container_identifiers', {}).get('auth_asym_id')
                                
                                if asym_id == chain:
                                    uniprot_ids = instance.get('polymer_entity', {}).get('rcsb_polymer_entity_container_identifiers', {}).get('uniprot_ids', [])
                                    
                                    if uniprot_ids:
                                        uniprot_id = uniprot_ids[0]
                                        print(f"Found UniProt ID from GraphQL: {uniprot_id}")
                                        break
            
            if uniprot_id:
                # Get UniProt sequence
                uniprot_url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
                uniprot_response = requests.get(uniprot_url)
                
                if uniprot_response.status_code == 200:
                    fasta_lines = uniprot_response.text.strip().split('\n')
                    # Skip the header line
                    uniprot_sequence = ''.join(fasta_lines[1:])
                    
                    # Now we need to align this with the PDB sequence to get the exact displayed region
                    # This is complex, but we can use a simpler approach for now
                    print(f"Got full UniProt sequence with length {len(uniprot_sequence)}")
                    
                    # Try to get PDB to UniProt mapping
                    mapping_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
                    mapping_response = requests.get(mapping_url)
                    
                    if mapping_response.status_code == 200:
                        mapping_data = mapping_response.json()
                        pdb_mappings = mapping_data.get(pdb_id.lower(), {}).get('UniProt', {}).get(uniprot_id, {}).get('mappings', [])
                        
                        for mapping in pdb_mappings:
                            chain_id = mapping.get('chain_id')
                            
                            if chain_id == chain:
                                uniprot_start = mapping.get('unp_start')
                                uniprot_end = mapping.get('unp_end')
                                
                                if uniprot_start is not None and uniprot_end is not None:
                                    print(f"Found mapping from PDB to UniProt: {uniprot_start}-{uniprot_end}")
                                    # Convert to 0-based indexing
                                    uniprot_start = int(uniprot_start) - 1
                                    uniprot_end = int(uniprot_end)
                                    
                                    return uniprot_sequence[uniprot_start:uniprot_end]
                    
                    # If we couldn't get the exact mapping, return the full sequence
                    return uniprot_sequence
    except Exception as e:
        print(f"Error with FASTA approach: {str(e)}")
    
    # As a last resort, try to use the SIFTS mapping
    print("Trying SIFTS mapping...")
    try:
        sifts_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
        sifts_response = requests.get(sifts_url)
        
        if sifts_response.status_code == 200:
            sifts_data = sifts_response.json()
            pdb_entry = sifts_data.get(pdb_id.lower(), {})
            
            # Look through all UniProt mappings
            for uniprot_id, uniprot_info in pdb_entry.get('UniProt', {}).items():
                mappings = uniprot_info.get('mappings', [])
                
                for mapping in mappings:
                    chain_id = mapping.get('chain_id')
                    
                    if chain_id == chain:
                        print(f"Found mapping for chain {chain} to UniProt {uniprot_id}")
                        
                        # Get the UniProt sequence
                        uniprot_url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
                        uniprot_response = requests.get(uniprot_url)
                        
                        if uniprot_response.status_code == 200:
                            fasta_lines = uniprot_response.text.strip().split('\n')
                            # Skip the header line
                            uniprot_sequence = ''.join(fasta_lines[1:])
                            
                            # Get the mapped region
                            uniprot_start = mapping.get('unp_start')
                            uniprot_end = mapping.get('unp_end')
                            
                            if uniprot_start is not None and uniprot_end is not None:
                                print(f"Mapped region: {uniprot_start}-{uniprot_end}")
                                # Convert to 0-based indexing
                                uniprot_start = int(uniprot_start) - 1
                                uniprot_end = int(uniprot_end)
                                
                                return uniprot_sequence[uniprot_start:uniprot_end]
                            
                            # If no specific region is mapped, return the full sequence
                            return uniprot_sequence
    except Exception as e:
        print(f"Error with SIFTS approach: {str(e)}")
    
    print("Could not extract UniProt sequence")
    return None

def main():
    parser = argparse.ArgumentParser(description="Extract UniProt sequence for a specific chain")
    parser.add_argument("pdb_id", help="PDB ID (e.g., 5M8I)")
    parser.add_argument("chain", help="Chain identifier (e.g., A)")
    args = parser.parse_args()
    
    sequence = get_uniprot_sequence(args.pdb_id, args.chain)
    
    if sequence:
        print("\nExtracted UniProt sequence:")
        print(sequence)
    else:
        print("\nFailed to extract UniProt sequence")

if __name__ == "__main__":
    main()