#!/usr/bin/env python3

import random
import numpy as np
import os
import re
import sys
import argparse
import matplotlib
matplotlib.use('Agg')  # Use Agg backend for non-interactive plotting
import matplotlib.pyplot as plt
from pathlib import Path
import traceback

# --- IMPORTS AND CONSTANTS FOR ORGANISM FETCHER ---
import requests
import time
import json

# Configuration for organism_fetcher
RCSB_ENTRY_SUMMARY_API_URL_TEMPLATE = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_ENTITY_API_URL_TEMPLATE = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
REQUEST_TIMEOUT = 15
API_CALL_DELAY = 0.25 # Respect API call delays
# ----------------------------------------------------

# --- RRM Analyzer Configuration ---
TESTING = False  # SET THIS TO True FOR A QUICK TEST ON THE FIRST RRM ENTRY
WINDOW_SIZE = 5
DEFAULT_OUTPUT_DIR = "random_distribution_plots"
WINDOWS_SUBDIR = "windows"


# --- ORGANISM FETCHER FUNCTIONS ---
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
    Returns a string of comma-separated organism names or an error/status message.
    """
    entry_summary_url = RCSB_ENTRY_SUMMARY_API_URL_TEMPLATE.format(pdb_id=pdb_id_upper)
    all_entry_organisms_s1 = set()
    strategy1_response_data = None
    strategy1_http_error_status = None
    s1_failure_reason = "S1: Unknown reason"

    try: # Strategy 1
        response_s1 = requests.get(entry_summary_url, timeout=REQUEST_TIMEOUT)
        time.sleep(API_CALL_DELAY)
        response_s1.raise_for_status()
        strategy1_response_data = response_s1.json()

        polymer_entities_s1 = strategy1_response_data.get('polymer_entities')
        if polymer_entities_s1:
            for entity in polymer_entities_s1:
                orgs_from_entity = _parse_organisms_from_entity_data_obj(entity)
                all_entry_organisms_s1.update(orgs_from_entity)
            if all_entry_organisms_s1: return ", ".join(sorted(list(all_entry_organisms_s1)))
            s1_failure_reason = "S1: 'polymer_entities' present but no orgs found"
        else:
            top_level_orgs_s1 = _parse_organisms_from_entity_data_obj(strategy1_response_data)
            if top_level_orgs_s1: return ", ".join(sorted(list(top_level_orgs_s1)))
            s1_failure_reason = "S1: No 'polymer_entities' key and no top-level orgs"
    except requests.exceptions.HTTPError as http_err:
        strategy1_http_error_status = http_err.response.status_code if http_err.response else 'N/A'
        s1_failure_reason = f"S1: HTTP error {strategy1_http_error_status}"
        if strategy1_http_error_status == 404: return f"PDB ID {pdb_id_upper} not found (404)"
    except Exception as e:
        s1_failure_reason = f"S1: Exception ({str(e)[:30]})"

    all_entry_organisms_s2 = set()
    polymer_entity_ids_list = []

    if strategy1_response_data:
        container_ids = strategy1_response_data.get('rcsb_entry_container_identifiers')
        if container_ids: polymer_entity_ids_list = container_ids.get('polymer_entity_ids', [])

    if not polymer_entity_ids_list and strategy1_http_error_status != 404:
        try:
            response_s2_entry = requests.get(entry_summary_url, timeout=REQUEST_TIMEOUT)
            time.sleep(API_CALL_DELAY)
            response_s2_entry.raise_for_status()
            entry_data_for_s2 = response_s2_entry.json()
            container_ids = entry_data_for_s2.get('rcsb_entry_container_identifiers')
            if container_ids: polymer_entity_ids_list = container_ids.get('polymer_entity_ids', [])
        except Exception as e:
            return f"{s1_failure_reason}, S2: Failed to re-fetch entry for entity list ({str(e)[:50]})"

    if not polymer_entity_ids_list:
        return f"{s1_failure_reason}, S2: No polymer_entity_ids found for PDB {pdb_id_upper}"

    for entity_id in polymer_entity_ids_list:
        entity_api_url = RCSB_ENTITY_API_URL_TEMPLATE.format(pdb_id=pdb_id_upper, entity_id=entity_id)
        try:
            response_entity_s2 = requests.get(entity_api_url, timeout=REQUEST_TIMEOUT)
            time.sleep(API_CALL_DELAY)
            response_entity_s2.raise_for_status()
            entity_data_s2 = response_entity_s2.json()
            orgs_from_this_entity = _parse_organisms_from_entity_data_obj(entity_data_s2)
            all_entry_organisms_s2.update(orgs_from_this_entity)
        except Exception:
            pass # Silently skip entities that error out in S2

    if all_entry_organisms_s2:
        return ", ".join(sorted(list(all_entry_organisms_s2)))
    else:
        return f"{s1_failure_reason}, S2: No organisms found after checking all entities."
# -------------------------------------------------------------------


# --- RRMDataset Class ---
class RRMDataset:
    """Class for handling RRM dataset files."""
    def __init__(self, file_path):
        self.file_path = file_path
        self.entries = []
        self.parse_file()

    def parse_file(self):
        try:
            with open(self.file_path, 'r') as f:
                content = f.read()
            raw_entries = re.split(r'>>([^\n]+)\n', content)
            if not raw_entries[0].strip():
                raw_entries = raw_entries[1:]

            for i in range(0, len(raw_entries), 2):
                if i + 1 < len(raw_entries):
                    header_line = raw_entries[i].strip()
                    data_block = raw_entries[i+1]
                    
                    pdb_match_header = re.match(r'([0-9A-Za-z]{4})_([A-Za-z0-9]+)', header_line)
                    if not pdb_match_header:
                        print(f"Warning: Could not parse PDB ID and Chain from header: {header_line}")
                        continue
                    
                    pdb_id = pdb_match_header.group(1)
                    chain_id = pdb_match_header.group(2)
                    score_file_key = f"{pdb_id}_{chain_id}"

                    try:
                        entry = self._extract_entry_data(header_line, data_block, pdb_id, chain_id, score_file_key)
                        if entry:
                            self.entries.append(entry)
                    except Exception as e:
                        print(f"Error parsing entry with header: {header_line}\nDetails: {str(e)}\n{traceback.format_exc()}")
            print(f"Successfully parsed {len(self.entries)} entries from dataset.")
        except Exception as e:
            print(f"Error reading or parsing dataset file: {str(e)}\n{traceback.format_exc()}")

    def _extract_entry_data(self, header_line, data, pdb_id_from_header, chain_id_from_header, score_file_key):
        uniprot_match = re.search(r'UniProt ID: ([^\n]+)', data)
        uniprot_id_full = uniprot_match.group(1).strip() if uniprot_match else "UnknownUniprot"
        
        rrm_num_match_uniprot = re.search(r'_RRM(\d+)', uniprot_id_full)
        rrm_num = rrm_num_match_uniprot.group(1) if rrm_num_match_uniprot else None

        if rrm_num:
            rrm_identifier = f"{uniprot_id_full.split('_RRM')[0]}_RRM{rrm_num}"
        else:
            rrm_identifier = uniprot_id_full

        nuc_match = re.search(r'Nucleotide Sequence: ([^\n]+)', data)
        nuc_seq = nuc_match.group(1).strip() if nuc_match else ""
        rna_seq = self._dna_to_rna(nuc_seq)

        prot_match = re.search(r'Protein Sequence: ([^\n]+)', data)
        prot_seq = prot_match.group(1).strip() if prot_match else ""

        shift = 0
        shift_match = re.search(r'Numbering Shift.*?:([^\n]+)', data)
        if shift_match:
            shift_text = shift_match.group(1).strip()
            num_match = re.search(r'(-?\d+)', shift_text)
            if num_match: shift = int(num_match.group(1))

        connections = []
        connections_match = re.search(r'List of connections: ([^\n]+)', data)
        if connections_match:
            connections_str = connections_match.group(1).strip()
            try:
                connections = [int(x.strip()) for x in connections_str.split(',') if x.strip()]
            except ValueError:
                try:
                    connections = [int(x.strip()) for x in connections_str.split() if x.strip()]
                except ValueError:
                    print(f"Warning: Could not parse connections for {score_file_key}: '{connections_str}'")
        
        # These connections are assumed to be protein-level indices.
        # Shift is applied to these protein-level indices.
        shifted_connections = sorted([c + shift for c in connections])

        return {
            'raw_header': header_line,
            'pdb_id': pdb_id_from_header, 
            'chain_id': chain_id_from_header,
            'score_file_key': score_file_key,
            'uniprot_id_full': uniprot_id_full,
            'rrm_identifier': rrm_identifier,
            'rna_sequence': rna_seq,
            'protein_sequence': prot_seq,
            'numbering_shift': shift,
            'connections': shifted_connections # These are protein-level connection indices
        }

    def _dna_to_rna(self, sequence):
        if not sequence: return ""
        rna = sequence.upper().replace('T', 'U')
        return re.sub(r'[^ACGU]', '', rna)

# --- ScoreLoader Class ---
class ScoreLoader:
    def __init__(self, windows_dir_base):
        self.windows_dir_base = Path(windows_dir_base)
        if not self.windows_dir_base.is_dir():
            script_dir = Path(__file__).parent
            self.windows_dir_base = script_dir / windows_dir_base
            if not self.windows_dir_base.is_dir():
                 raise FileNotFoundError(f"Windows scores directory '{windows_dir_base}' not found.")
        print(f"ScoreLoader: Using windows scores directory: {self.windows_dir_base.resolve()}")
        self.score_cache = {}
        self.nucleotide_map = {'A': 0, 'C': 1, 'G': 2, 'U': 3}

    def _window_to_index(self, window_sequence):
        if len(window_sequence) != WINDOW_SIZE:
            raise ValueError(f"Window '{window_sequence}' must be length {WINDOW_SIZE}.")
        index = 0
        for i, nucleotide in enumerate(window_sequence):
            if nucleotide not in self.nucleotide_map:
                raise ValueError(f"Invalid nucleotide '{nucleotide}' in window '{window_sequence}'.")
            index += self.nucleotide_map[nucleotide] * (4 ** (WINDOW_SIZE - 1 - i))
        return index

    def _load_score_file(self, score_file_key):
        if score_file_key in self.score_cache:
            return self.score_cache[score_file_key]
        score_file_path = self.windows_dir_base / f"{score_file_key}.txt"
        if not score_file_path.exists(): return None
        try:
            with open(score_file_path, 'r') as f:
                scores = [float(line.strip()) for line in f if line.strip()]
            if len(scores) != 4**WINDOW_SIZE:
                print(f"Warning: Score file {score_file_path} has {len(scores)} lines, expected {4**WINDOW_SIZE}.")
                return None
            self.score_cache[score_file_key] = scores
            return scores
        except Exception as e:
            print(f"Error loading score file {score_file_path}: {e}")
            return None

    def get_score(self, score_file_key, window_sequence):
        all_scores = self._load_score_file(score_file_key)
        if all_scores is None: return 0.0
        try:
            index = self._window_to_index(window_sequence.upper())
            return all_scores[index] if 0 <= index < len(all_scores) else 0.0
        except ValueError: return 0.0 
        except Exception: return 0.0


# --- ConnectionAnalyzer Class ---
class ConnectionAnalyzer:
    def __init__(self, score_loader_instance, all_dataset_entries):
        self.score_loader = score_loader_instance
        self.all_entries = all_dataset_entries

    def identify_connection_blocks(self, connections):
        if not connections: return []
        connections = sorted(list(set(connections)))
        if not connections: return []
        blocks, current_block = [], [connections[0]]
        for i in range(1, len(connections)):
            if connections[i] - connections[i-1] == 1:
                current_block.append(connections[i])
            else:
                blocks.append(current_block)
                current_block = [connections[i]]
        blocks.append(current_block)
        return blocks

    def get_rna_positions_for_connection(self, protein_connection_idx, rna_seq_length):
        """Maps a protein connection index to a triplet of RNA indices."""
        rna_start = protein_connection_idx * 3
        rna_end = protein_connection_idx * 3 + 2
        return (rna_start, rna_end) if rna_end < rna_seq_length else None

    def get_rna_windows_for_block(self, block_of_protein_connections, rna_seq_length, window_size=WINDOW_SIZE):
        """
        For a block of protein_connection_indices, determine RNA windows.
        Each protein connection maps to an RNA triplet (standard 1:3 mapping).
        Returns a sorted list of unique RNA window start indices.
        """
        involved_rna_indices = set() # Set of individual RNA nucleotide indices
        for protein_conn_idx in block_of_protein_connections:
            pos_triplet = self.get_rna_positions_for_connection(protein_conn_idx, rna_seq_length)
            if pos_triplet: 
                involved_rna_indices.update(range(pos_triplet[0], pos_triplet[1] + 1))
        
        if not involved_rna_indices: return []
        
        # Find all window start positions whose windows overlap with involved_rna_indices
        overlapping_window_starts = set()
        min_rna_idx_affected = min(involved_rna_indices)
        max_rna_idx_affected = max(involved_rna_indices)

        # A window starting at 's' covers [s, s + window_size - 1]
        # Smallest 's' such that s + window_size - 1 >= min_rna_idx_affected
        #   => s >= min_rna_idx_affected - window_size + 1
        # Largest 's' such that s <= max_rna_idx_affected
        start_range_min = max(0, min_rna_idx_affected - window_size + 1)
        start_range_max = min(rna_seq_length - window_size, max_rna_idx_affected) # Max start for a full window

        for s in range(start_range_min, start_range_max + 1):
            # Window is [s, s + window_size -1]
            current_window_indices_set = set(range(s, s + window_size))
            if not current_window_indices_set.isdisjoint(involved_rna_indices):
                overlapping_window_starts.add(s)
        return sorted(list(overlapping_window_starts))

    def calculate_connection_score(self, score_file_key_for_scores, rna_sequence_to_scan,
                                   protein_connections_list, # These are always protein-level indices
                                   window_size=WINDOW_SIZE,
                                   use_direct_rna_mapping_for_original=False, # If True, protein_idx -> rna_idx (1:1)
                                   max_windows_to_score=None): # Max unique windows to sum scores for
        if not rna_sequence_to_scan or len(rna_sequence_to_scan) < window_size:
            return 0.0, [], []
        
        # protein_connections_list are protein residue indices.
        # blocks are contiguous blocks of these protein residue indices.
        blocks = self.identify_connection_blocks(protein_connections_list)
        if not blocks:
            return 0.0, [], []

        # block_sizes are lengths of protein residue blocks.
        block_sizes = [len(b) for b in blocks]
        
        potential_unique_window_starts = set()

        for protein_block in blocks: # protein_block is a list of protein residue indices
            rna_window_starts_for_this_block = []
            if use_direct_rna_mapping_for_original:
                # Map protein_idx directly to rna_idx (1:1)
                # This path is for original score calculation as per request.
                involved_direct_rna_indices = set()
                for protein_idx in protein_block:
                    rna_idx = protein_idx # Direct mapping
                    if 0 <= rna_idx < len(rna_sequence_to_scan):
                        involved_direct_rna_indices.add(rna_idx)
                
                if involved_direct_rna_indices:
                    min_rna_idx_affected = min(involved_direct_rna_indices)
                    max_rna_idx_affected = max(involved_direct_rna_indices)
                    
                    start_range_min = max(0, min_rna_idx_affected - window_size + 1)
                    start_range_max = min(len(rna_sequence_to_scan) - window_size, max_rna_idx_affected)

                    current_block_overlapping_starts = set()
                    for s in range(start_range_min, start_range_max + 1):
                        current_window_indices_set = set(range(s, s + window_size))
                        if not current_window_indices_set.isdisjoint(involved_direct_rna_indices):
                            current_block_overlapping_starts.add(s)
                    rna_window_starts_for_this_block = sorted(list(current_block_overlapping_starts))
            else:
                # Standard: protein_idx maps to RNA triplet (1:3). Used for random scores.
                rna_window_starts_for_this_block = self.get_rna_windows_for_block(
                    protein_block, len(rna_sequence_to_scan), window_size
                )
            
            for rna_start_idx in rna_window_starts_for_this_block:
                potential_unique_window_starts.add(rna_start_idx)

        # Now score these unique windows, up to the limit
        total_score = 0.0
        all_scored_rna_window_starts = set() # Tracks unique window starts actually scored
        
        # If max_windows_to_score is 0, score nothing.
        if max_windows_to_score == 0:
             return 0.0, block_sizes, []

        sorted_potential_starts = sorted(list(potential_unique_window_starts))

        for rna_start_idx in sorted_potential_starts:
            # Check if we have reached the limit on the number of unique windows to score
            if max_windows_to_score is not None and len(all_scored_rna_window_starts) >= max_windows_to_score:
                break 

            window_seq = rna_sequence_to_scan[rna_start_idx : rna_start_idx + window_size]
            # Ensure window is valid (should be, if generation logic is correct)
            if len(window_seq) == window_size:
                total_score += self.score_loader.get_score(score_file_key_for_scores, window_seq)
                all_scored_rna_window_starts.add(rna_start_idx)
        
        return total_score, block_sizes, sorted(list(all_scored_rna_window_starts))


    def generate_random_connections(self, block_sizes, rna_seq_length, min_gap_between_rrm_blocks=1):
        # This function generates protein-level connection indices.
        # block_sizes are lengths of protein residue blocks.
        num_rrm_protein_slots = rna_seq_length // 3
        if num_rrm_protein_slots <= 0 or not block_sizes: return []
        
        k = len(block_sizes) # Number of blocks
        total_rrm_residues_in_blocks = sum(block_sizes)
        
        # Minimum slots needed for all blocks and minimum gaps between them
        min_required_slots = total_rrm_residues_in_blocks + max(0, k - 1) * min_gap_between_rrm_blocks
        if min_required_slots > num_rrm_protein_slots: return [] # Not enough space

        # Remaining slots to distribute as additional gaps
        extra_slots = num_rrm_protein_slots - min_required_slots
        num_gap_bins = k + 1 # k blocks means k+1 possible gap locations (before first, between, after last)
        
        additional_gaps = [0] * num_gap_bins
        if extra_slots > 0:
            # Distribute extra_slots into num_gap_bins using "stars and bars" sampling
            # Generate num_gap_bins-1 random dividers in a range of extra_slots + num_gap_bins -1 positions
            dividers = sorted(random.sample(range(extra_slots + num_gap_bins - 1), num_gap_bins - 1))
            last_div = -1
            for i, d in enumerate(dividers):
                additional_gaps[i] = d - last_div - 1
                last_div = d
            additional_gaps[num_gap_bins - 1] = (extra_slots + num_gap_bins - 1) - last_div - 1
        
        # Combine minimum gaps with additional gaps
        final_gaps = [0] * num_gap_bins
        final_gaps[0] = additional_gaps[0] # Gap before the first block
        for i in range(k -1): # Gaps between k blocks (k-1 such gaps)
             final_gaps[i+1] = additional_gaps[i+1] + min_gap_between_rrm_blocks
        if k > 0 : # Space after last block (if there are blocks)
            final_gaps[k] = additional_gaps[k] # The k-th gap (index k) is after the k-th block

        random_connections_protein_indices = []
        current_protein_pos = final_gaps[0] # Start after the first gap
        
        shuffled_block_sizes = list(block_sizes) # Operate on a copy
        random.shuffle(shuffled_block_sizes)

        for i, block_len in enumerate(shuffled_block_sizes):
            random_connections_protein_indices.extend(range(current_protein_pos, current_protein_pos + block_len))
            current_protein_pos += block_len # Move past the current block
            if i < k - 1: # If not the last block, add the next inter-block gap
                current_protein_pos += final_gaps[i+1]
        
        return sorted(list(set(random_connections_protein_indices)))

    def calculate_random_distribution(self, current_rrm_entry, n_iterations=1000):
        # Original score calculation:
        # protein_connections are from the entry.
        # use_direct_rna_mapping_for_original=True means protein_idx -> rna_idx (1:1 mapping).
        original_score, original_protein_block_sizes, original_scored_rna_window_starts = \
            self.calculate_connection_score(
                current_rrm_entry['score_file_key'],
                current_rrm_entry['rna_sequence'], 
                current_rrm_entry['connections'], # These are protein-level connection indices
                WINDOW_SIZE,
                use_direct_rna_mapping_for_original=True # Request 1: Use 1:1 mapping for original
            )
        
        # Number of unique windows scored for the original entry. This will be the limit for random.
        num_windows_in_original = len(original_scored_rna_window_starts)

        if not original_protein_block_sizes: # No blocks means no connections to randomize
            return original_score, [], [], 0 # No random scores if no original blocks

        random_scores, successful_iterations = [], 0
        for _ in range(n_iterations):
            target_entry = random.choice(self.all_entries)
            if not target_entry['rna_sequence'] or len(target_entry['rna_sequence']) < WINDOW_SIZE:
                continue
            
            # Generate random protein-level connections using original (protein) block sizes
            rand_protein_conns = self.generate_random_connections(
                original_protein_block_sizes, # These are protein block sizes
                len(target_entry['rna_sequence'])
            )
            if not rand_protein_conns:
                continue

            # Random score calculation:
            # rand_protein_conns are protein-level.
            # use_direct_rna_mapping_for_original=False means standard protein_idx -> RNA triplet mapping.
            # max_windows_to_score limits the number of scored windows.
            score_iter, _, _ = self.calculate_connection_score(
                target_entry['score_file_key'],
                target_entry['rna_sequence'],
                rand_protein_conns, # These are protein-level connection indices
                WINDOW_SIZE,
                use_direct_rna_mapping_for_original=False, # Standard 1:3 mapping for random
                max_windows_to_score=num_windows_in_original # Request 2: Limit scored windows
            )
            random_scores.append(score_iter)
            successful_iterations += 1
            
        return original_score, random_scores, original_protein_block_sizes, successful_iterations

# --- ConnectionDistributionPlot Class ---
class ConnectionDistributionPlot:
    def __init__(self, output_dir_path):
        self.output_dir = Path(output_dir_path)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def plot_distribution(self, original_score, random_scores, rrm_entry, 
                          total_attempted_iterations, successful_iterations_count, 
                          block_sizes_info="", organism_str="N/A"):
        fig, ax = plt.subplots(figsize=(12, 8))
        entry_label = f"{rrm_entry['score_file_key']} (UniProt: {rrm_entry['uniprot_id_full']})"
        org_display = organism_str[:70] + "..." if len(organism_str) > 70 else organism_str
        
        title_l1 = f'Connection Score Distribution for {entry_label}'
        title_l2 = f'Organism(s): {org_display}'
        title_l3 = f'({successful_iterations_count}/{total_attempted_iterations} Successful Random Shuffles vs. Original)'
        full_title = f"{title_l1}\n{title_l2}\n{title_l3}"
        if block_sizes_info and block_sizes_info != "N/A":
             full_title += f"\nOriginal Protein Block Sizes: {block_sizes_info}" # Clarified block sizes
        ax.set_title(full_title, fontsize=10)

        if not random_scores:
            ax.text(0.5, 0.5, "No successful random shuffles or no random scores generated.", ha='center', va='center')
            if original_score is not None:
                 ax.axvline(x=original_score, color='red', ls='--', lw=2, label=f'Original: {original_score:.3f}')
                 ax.legend()
        else:
            ax.hist(random_scores, bins=30, alpha=0.75, color='skyblue', ec='black', label=f'Random Scores ({successful_iterations_count})')
            ax.axvline(x=original_score, color='red', ls='--', lw=2, label=f'Original: {original_score:.3f}')
            mean_rand, std_rand = np.mean(random_scores), np.std(random_scores)
            ax.axvline(x=mean_rand, color='darkblue', ls=':', lw=1.5, label=f'Random Mean: {mean_rand:.3f}')
            
            count_ext = sum(1 for x in random_scores if x >= original_score)
            p_val = (count_ext + 1) / (len(random_scores) + 1) # (N_extreme + 1) / (N_random + 1)
            z_sc = (original_score - mean_rand) / std_rand if std_rand > 1e-9 else 0 # Avoid division by zero
            stats_text = (f'Original: {original_score:.3f}\nRandom Mean: {mean_rand:.3f}\n'
                          f'Random StdDev: {std_rand:.3f}\nZ-score: {z_sc:.2f}\nP-value: {p_val:.4f}')
            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=9, va='top', 
                    bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))
            ax.legend()

        ax.set_xlabel('Sum of Connection Window Scores'); ax.set_ylabel('Frequency')
        ax.grid(True, ls='--', alpha=0.7)
        plt.tight_layout(rect=[0, 0, 1, 0.93]) # Adjust layout to prevent title overlap
        plot_filename = self.output_dir / f"{rrm_entry['score_file_key']}_distribution_org.png"
        try: plt.savefig(plot_filename, dpi=150)
        except Exception as e: print(f"Error saving plot {plot_filename}: {e}")
        plt.close(fig)

    def plot_p_value_summary(self, results_list, plot_title_suffix="", filename_suffix=""):
        if not results_list:
            print(f"No results for P-value summary plot {plot_title_suffix}.")
            return
        # Filter out results where p_value might be None (e.g., if random_scores was empty)
        valid_plot_results = [res for res in results_list if res.get('p_value') is not None]
        if not valid_plot_results:
            print(f"No valid P-values for summary plot {plot_title_suffix}.")
            return

        sorted_results = sorted(valid_plot_results, key=lambda x: x['p_value'])
        rrm_labels = [f"{res.get('uniprot_id_full', 'N/A')} ({res.get('id', 'N/A')})" for res in sorted_results]
        p_values = [res['p_value'] for res in sorted_results]

        fig, ax = plt.subplots(figsize=(max(12, len(rrm_labels) * 0.4), 8))
        bar_colors = ['red' if p < 0.05 else 'blue' for p in p_values]
        bars = ax.bar(rrm_labels, p_values, color=bar_colors, ec='black')
        ax.axhline(y=0.05, color='green', ls='--', label='P=0.05 Threshold')

        if len(rrm_labels) <= 30: # Add text only for a reasonable number of bars
            for i, bar in enumerate(bars):
                yval = bar.get_height()
                txt_col = 'white' if bar_colors[i] == 'red' and yval > 0.02 else 'black'
                plt.text(bar.get_x() + bar.get_width()/2.0, yval + 0.002, f'{yval:.3f}', 
                         ha='center', va='bottom', fontsize=7, rotation=90, color=txt_col)

        ax.set_xlabel('RRM Entry: UniProt ID (PDBID_Chain)'); ax.set_ylabel('P-value (Original >= Random)')
        ax.set_title(f'P-value Summary {plot_title_suffix} ({len(rrm_labels)} RRMs, Sorted)')
        plt.xticks(rotation=90, ha='right', fontsize=8 if len(rrm_labels) > 20 else 9)
        # Adjust legend position based on number of labels to prevent overlap with x-axis labels
        legend_y_offset = -0.35 if len(rrm_labels) > 25 else (-0.25 if len(rrm_labels) > 15 else -0.15)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, legend_y_offset), 
                  fancybox=True, shadow=True, ncol=1)
        ax.grid(True, axis='y', ls='--', alpha=0.6)
        upper_y_limit = max(0.06, (max(p_values) if p_values else 0.06) * 1.15) # Ensure 0.05 line is visible
        plt.ylim(0, upper_y_limit)
        plt.tight_layout(rect=[0, 0.05, 1, 0.97]) # Adjust bottom margin for rotated labels
        summary_plot_filename = self.output_dir / f"p_value_summary_colored{filename_suffix}.png"
        try: plt.savefig(summary_plot_filename, dpi=150)
        except Exception as e: print(f"Error saving P-value summary: {e}")
        plt.close(fig)

# --- Main script logic ---
def main():
    parser = argparse.ArgumentParser(description='RRM Connection Score Analyzer.')
    parser.add_argument('dataset_path', type=str, help='Path to the extracted_data.txt file.')
    parser.add_argument('--output_dir', type=str, default=DEFAULT_OUTPUT_DIR,
                        help=f'Directory to save output plots (default: {DEFAULT_OUTPUT_DIR}).')
    parser.add_argument('--windows_scores_dir', type=str, default=WINDOWS_SUBDIR,
                        help=f'Directory containing RRM score files (default: {WINDOWS_SUBDIR}).')
    parser.add_argument('--iterations', type=int, default=5000, # Default changed to 5000 from 1000
                        help='Number of random iterations for distribution (default: 5000).')
    args = parser.parse_args()

    print("--- RRM Connection Score Analyzer (Organism Fetching Enabled) ---")
    if TESTING: print("!!! TESTING MODE ENABLED: Processing only the first RRM entry. !!!")

    rrm_data = RRMDataset(args.dataset_path)
    if not rrm_data.entries: print("No entries in dataset. Exiting."); return
    
    try: score_loader = ScoreLoader(args.windows_scores_dir)
    except FileNotFoundError as e: print(f"Error: {e}"); return
    
    conn_analyzer = ConnectionAnalyzer(score_loader, rrm_data.entries)
    plotter = ConnectionDistributionPlot(args.output_dir)
    
    organism_cache = {}
    all_results_for_summary = []
    entries_to_process = [rrm_data.entries[0]] if TESTING and rrm_data.entries else rrm_data.entries

    for i, rrm_entry in enumerate(entries_to_process):
        print(f"\nProcessing RRM ({i+1}/{len(entries_to_process)}): {rrm_entry['score_file_key']}")
        
        pdb_id_for_org = rrm_entry['pdb_id'].upper()
        organism_str = organism_cache.get(pdb_id_for_org)
        if organism_str is None:
            print(f"  Fetching organism for PDB ID: {pdb_id_for_org}...")
            organism_str = get_organisms_for_pdb_entry(pdb_id_for_org)
            organism_cache[pdb_id_for_org] = organism_str
            print(f"  Organism (fetched): {organism_str}")
        else:
            print(f"  Organism (cached): {organism_str}")
        is_homo_sapiens = "homo sapiens" in organism_str.lower() if isinstance(organism_str, str) else False

        current_entry_result = {
            'id': rrm_entry['score_file_key'],
            'uniprot_id_full': rrm_entry['uniprot_id_full'],
            'organism': organism_str,
            'is_homo_sapiens': is_homo_sapiens,
            'original_score': None, 'p_value': None, 'z_score': None,
            'valid_for_summary': False # Default to False
        }

        if not rrm_entry['connections'] or not rrm_entry['rna_sequence'] or len(rrm_entry['rna_sequence']) < WINDOW_SIZE:
            reason = "No connections" if not rrm_entry['connections'] else ("RNA too short/missing" if not rrm_entry['rna_sequence'] or len(rrm_entry['rna_sequence']) < WINDOW_SIZE else "Unknown reason")
            print(f"  Skipping {rrm_entry['score_file_key']}: {reason}.")
            all_results_for_summary.append(current_entry_result) # Append with default/error values
            continue
        
        try:
            original_score, random_scores, protein_block_sizes, successful_shuffles = \
                conn_analyzer.calculate_random_distribution(rrm_entry, args.iterations)
            
            current_entry_result['original_score'] = original_score
            print(f"  Original Score: {original_score:.4f}")

            if random_scores and successful_shuffles > 0 : # Ensure there are scores and iterations
                mean_rand, std_rand = np.mean(random_scores), np.std(random_scores)
                print(f"  {len(random_scores)} random scores from {successful_shuffles}/{args.iterations} shuffles. Mean: {mean_rand:.4f}, StdDev: {std_rand:.4f}")
                
                count_ext = sum(1 for x in random_scores if x >= original_score)
                p_val = (count_ext + 1) / (len(random_scores) + 1) # Proper p-value for permutation test
                current_entry_result['p_value'] = p_val
                
                if std_rand > 1e-9: # Avoid division by zero or tiny std dev
                    current_entry_result['z_score'] = (original_score - mean_rand) / std_rand
                else:
                    current_entry_result['z_score'] = 0.0 if abs(original_score - mean_rand) < 1e-9 else float('inf') * np.sign(original_score - mean_rand)

                # Validate for summary based on percentage difference - this seems application specific, keeping it
                perc_diff = abs(original_score - mean_rand) / abs(mean_rand) * 100 if abs(mean_rand) > 1e-9 else (0 if abs(original_score) < 1e-9 else float('inf'))
                print(f"  P-value: {p_val:.4f}, Z-score: {current_entry_result['z_score']:.2f}")
                print(f"  Percentage difference from random mean: {perc_diff:.2f}%")
                
                # Condition for validity in summary plot - you might adjust this
                if perc_diff <= 50.0 or successful_shuffles < 100 : # Example: if few shuffles, maybe still include
                    current_entry_result['valid_for_summary'] = True
                else:
                    print(f"  Skipping from P-value summary: Difference > 50% and sufficient shuffles.")

            elif args.iterations > 0:
                print(f"  No successful random shuffles or no random scores generated out of {args.iterations} attempts.")
                current_entry_result['p_value'] = None # Cannot calculate p-value
                current_entry_result['z_score'] = None
            
            # Always append, even if some calculations failed; 'valid_for_summary' controls actual plotting inclusion
            all_results_for_summary.append(current_entry_result)
            
            plotter.plot_distribution(original_score, random_scores, rrm_entry, args.iterations, 
                                      successful_shuffles, str(protein_block_sizes or "N/A"), organism_str)
            print(f"  Saved distribution plot for {rrm_entry['score_file_key']}")

        except Exception as e:
            print(f"  ERROR processing {rrm_entry['score_file_key']}: {e}\n{traceback.format_exc()}")
            # Append with default/error values if an exception occurs mid-processing for this entry
            current_entry_result['error_message'] = str(e) # Store error for debugging if needed
            all_results_for_summary.append(current_entry_result)


    if not (TESTING and len(entries_to_process) == 1):
        valid_results_general = [r for r in all_results_for_summary if r['valid_for_summary'] and r.get('p_value') is not None]
        if valid_results_general:
            print("\nGenerating general P-value summary plot...")
            plotter.plot_p_value_summary(valid_results_general, plot_title_suffix="(All Valid)", filename_suffix="_all_valid")
        else: print("No RRMs met criteria for the general P-value summary.")

        valid_results_hs = [r for r in valid_results_general if r['is_homo_sapiens']] # Already filtered by valid_for_summary
        if valid_results_hs:
            print("\nGenerating Homo sapiens specific P-value summary plot...")
            plotter.plot_p_value_summary(valid_results_hs, plot_title_suffix="(Homo sapiens)", filename_suffix="_homo_sapiens")
        elif valid_results_general : # Only print this if there were general results but no HS ones
            print("No Homo sapiens RRMs met criteria for P-value summary.")
            
    print("\n--- Analysis Complete ---")

if __name__ == "__main__":
    main()