import re
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use Agg backend for non-interactive plotting
import matplotlib.pyplot as plt
from pathlib import Path
import random
import argparse
import traceback

# Global configuration - can be overridden by args
WINDOW_SIZE = 5
NUCLEOTIDE_MAP = {'A': 0, 'C': 1, 'G': 2, 'U': 3}

class RRMDataEntry:
    """Holds data for a single RRM entry and its analysis results."""
    def __init__(self, raw_header, pdb_id, chain_id, uniprot_id_full, rrm_identifier,
                 rna_sequence, protein_sequence, numbering_shift, connections,
                 score_file_key):
        self.raw_header = raw_header
        self.pdb_id = pdb_id
        self.chain_id = chain_id
        self.uniprot_id_full = uniprot_id_full # This is the one to use, e.g., "P0A7Z4_RRM1" or "P0A7Z4"
        self.rrm_identifier = rrm_identifier  # More specific identifier like P0A7Z4_RRM1 or PDBID_CHAIN_RRMx
        self.rna_sequence = rna_sequence
        self.protein_sequence = protein_sequence
        self.numbering_shift = numbering_shift # Shift for connection indices
        self.connections = connections  # RRM-level indices, already shifted
        self.score_file_key = score_file_key  # e.g., 1A9N_B, used for score file lookup

        # Analysis results to be filled later
        self.original_score = None
        self.num_windows_in_original_score = 0  # This is 'N'
        self.random_score_distribution = []
        self.p_value = None
        self.z_score = None

class DatasetParser:
    """Parses the RRM dataset file."""
    def __init__(self, file_path):
        self.file_path = file_path
        self.entries = []

    def _dna_to_rna(self, sequence):
        if not sequence: return ""
        rna = sequence.upper().replace('T', 'U')
        return re.sub(r'[^ACGU]', '', rna) # Keep only valid RNA characters

    def parse_file(self):
        print(f"Parsing dataset file: {self.file_path}")
        try:
            with open(self.file_path, 'r') as f:
                content = f.read()
            
            raw_entries_parts = re.split(r'(>>[^\n]+\n)', content)
            
            if not raw_entries_parts[0].strip(): 
                raw_entries_parts = raw_entries_parts[1:]

            parsed_count = 0
            for i in range(0, len(raw_entries_parts), 2): 
                if i + 1 < len(raw_entries_parts):
                    header_line_full = raw_entries_parts[i].strip() 
                    data_block = raw_entries_parts[i+1]
                    
                    header_match = re.match(r'>>([0-9A-Za-z]{4})_([A-Za-z0-9]+)', header_line_full)
                    if not header_match:
                        print(f"Warning: Could not parse PDB ID and Chain from header line: '{header_line_full}'. Skipping.")
                        continue
                    
                    pdb_id = header_match.group(1)
                    chain_id = header_match.group(2)
                    score_file_key = f"{pdb_id}_{chain_id}"

                    uniprot_match = re.search(r'UniProt ID: ([^\n]+)', data_block)
                    uniprot_id_full = uniprot_match.group(1).strip() if uniprot_match else "UnknownUniprot"
                    
                    rrm_num_match_uniprot = re.search(r'_RRM(\d+)', uniprot_id_full)
                    rrm_num_from_uniprot = rrm_num_match_uniprot.group(1) if rrm_num_match_uniprot else None

                    # Construct rrm_identifier (used internally if needed, but uniprot_id_full is primary for labels)
                    if rrm_num_from_uniprot:
                        base_uniprot = uniprot_id_full.split('_RRM')[0]
                        rrm_identifier_derived = f"{base_uniprot}_RRM{rrm_num_from_uniprot}"
                    else:
                        rrm_identifier_match_header = re.search(r'([0-9A-Za-z]{4}_[A-Za-z0-9]+(_RRM\d+)?)', header_line_full)
                        if rrm_identifier_match_header and '_RRM' in rrm_identifier_match_header.group(1):
                             rrm_identifier_derived = rrm_identifier_match_header.group(1)
                        else: 
                            rrm_identifier_derived = uniprot_id_full if uniprot_id_full != "UnknownUniprot" else score_file_key

                    nuc_match = re.search(r'Nucleotide Sequence: ([^\n]+)', data_block)
                    nuc_seq = nuc_match.group(1).strip() if nuc_match else ""
                    rna_seq = self._dna_to_rna(nuc_seq)

                    prot_match = re.search(r'Protein Sequence: ([^\n]+)', data_block)
                    prot_seq = prot_match.group(1).strip() if prot_match else ""

                    shift = 0
                    shift_match = re.search(r'Numbering Shift.*?:([^\n]+)', data_block)
                    if shift_match:
                        shift_text = shift_match.group(1).strip()
                        num_match_shift = re.search(r'(-?\d+)', shift_text)
                        if num_match_shift: shift = int(num_match_shift.group(1))

                    connections_list = []
                    connections_match = re.search(r'List of connections: ([^\n]+)', data_block)
                    if connections_match:
                        connections_str = connections_match.group(1).strip()
                        try:
                            connections_list = [int(x.strip()) for x in connections_str.split(',') if x.strip()]
                        except ValueError:
                            try:
                                connections_list = [int(x.strip()) for x in connections_str.split() if x.strip()]
                            except ValueError:
                                print(f"Warning: Could not parse connections for {score_file_key}: '{connections_str}'")
                    
                    shifted_connections = sorted([c + shift for c in connections_list])

                    entry = RRMDataEntry(
                        raw_header=header_line_full,
                        pdb_id=pdb_id,
                        chain_id=chain_id,
                        uniprot_id_full=uniprot_id_full,
                        rrm_identifier=rrm_identifier_derived,
                        rna_sequence=rna_seq,
                        protein_sequence=prot_seq,
                        numbering_shift=shift,
                        connections=shifted_connections,
                        score_file_key=score_file_key
                    )
                    self.entries.append(entry)
                    parsed_count += 1
            print(f"Successfully parsed {parsed_count} entries from dataset.")
        except Exception as e:
            print(f"Error reading or parsing dataset file: {str(e)}\n{traceback.format_exc()}")
        return self.entries

class ScoreCalculator:
    """Calculates original and random scores based on pre-computed window scores."""
    def __init__(self, windows_dir_base):
        self.windows_dir_base = Path(windows_dir_base)
        if not self.windows_dir_base.is_dir():
            script_dir = Path(__file__).resolve().parent
            self.windows_dir_base = script_dir / windows_dir_base
            if not self.windows_dir_base.is_dir():
                 raise FileNotFoundError(f"Windows scores directory '{windows_dir_base}' not found at {self.windows_dir_base.resolve()} or directly.")
        print(f"ScoreCalculator: Using windows scores directory: {self.windows_dir_base.resolve()}")
        self.score_cache = {} 

    def _window_to_index(self, window_sequence):
        global WINDOW_SIZE, NUCLEOTIDE_MAP 
        if len(window_sequence) != WINDOW_SIZE:
            return -1  
        index = 0
        for i, nucleotide in enumerate(window_sequence):
            if nucleotide not in NUCLEOTIDE_MAP:
                return -1  
            index += NUCLEOTIDE_MAP[nucleotide] * (4 ** (WINDOW_SIZE - 1 - i))
        return index

    def _load_score_file(self, score_file_key):
        global WINDOW_SIZE
        if score_file_key in self.score_cache:
            return self.score_cache[score_file_key]
        
        score_file_path = self.windows_dir_base / f"{score_file_key}.txt"
        if not score_file_path.exists():
            print(f"Warning: Score file not found: {score_file_path}")
            self.score_cache[score_file_key] = None 
            return None
        
        try:
            with open(score_file_path, 'r') as f:
                scores = [float(line.strip()) for line in f if line.strip()]
            
            expected_lines = 4**WINDOW_SIZE
            if len(scores) != expected_lines:
                print(f"Warning: Score file {score_file_path} has {len(scores)} lines, expected {expected_lines}. Scores will not be used for {score_file_key}.")
                self.score_cache[score_file_key] = None
                return None
            
            self.score_cache[score_file_key] = scores
            return scores
        except Exception as e:
            print(f"Error loading score file {score_file_path}: {e}")
            self.score_cache[score_file_key] = None
            return None

    def get_score_for_window(self, score_file_key, window_sequence):
        all_scores = self._load_score_file(score_file_key)
        if all_scores is None:
            return 0.0 

        index = self._window_to_index(window_sequence.upper())
        if index == -1: 
            return 0.0
        
        if 0 <= index < len(all_scores):
            return all_scores[index]
        else: 
            print(f"Warning: Index {index} out of bounds for score file {score_file_key} (len {len(all_scores)}).")
            return 0.0

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

    def get_rna_positions_for_rrm_connection(self, rrm_protein_connection_idx, rna_seq_length):
        rna_start = rrm_protein_connection_idx * 3
        rna_end = rna_start + 2 
        if rna_end < rna_seq_length:
            return (rna_start, rna_end)
        return None

    def get_rna_window_starts_for_block(self, block_of_rrm_connections, rna_seq_length):
        global WINDOW_SIZE
        involved_rna_indices = set()
        for rrm_conn_idx in block_of_rrm_connections:
            rna_triplet_pos = self.get_rna_positions_for_rrm_connection(rrm_conn_idx, rna_seq_length)
            if rna_triplet_pos:
                involved_rna_indices.update(range(rna_triplet_pos[0], rna_triplet_pos[1] + 1))
        
        if not involved_rna_indices: return set()

        overlapping_window_starts = set()
        for s in range(rna_seq_length - WINDOW_SIZE + 1): 
            current_window_indices = set(range(s, s + WINDOW_SIZE))
            if not current_window_indices.isdisjoint(involved_rna_indices):
                overlapping_window_starts.add(s)
        return overlapping_window_starts

    def calculate_original_score(self, entry: RRMDataEntry):
        global WINDOW_SIZE
        if not entry.rna_sequence or len(entry.rna_sequence) < WINDOW_SIZE:
            print(f"RNA sequence for {entry.score_file_key} is too short or missing.")
            entry.original_score = 0.0
            entry.num_windows_in_original_score = 0
            return
        if not entry.connections:
            print(f"No connections for {entry.score_file_key}.")
            entry.original_score = 0.0
            entry.num_windows_in_original_score = 0
            return

        connection_blocks = self.identify_connection_blocks(entry.connections)
        if not connection_blocks:
            entry.original_score = 0.0
            entry.num_windows_in_original_score = 0
            return

        all_contributing_rna_window_starts = set()
        for block in connection_blocks:
            window_starts_for_this_block = self.get_rna_window_starts_for_block(block, len(entry.rna_sequence))
            all_contributing_rna_window_starts.update(window_starts_for_this_block)

        original_score_sum = 0.0
        if not all_contributing_rna_window_starts:
            entry.original_score = 0.0
            entry.num_windows_in_original_score = 0
            return
            
        for rna_start_idx in sorted(list(all_contributing_rna_window_starts)):
            window_seq = entry.rna_sequence[rna_start_idx : rna_start_idx + WINDOW_SIZE]
            if len(window_seq) == WINDOW_SIZE: 
                score = self.get_score_for_window(entry.score_file_key, window_seq)
                original_score_sum += score
        
        entry.original_score = original_score_sum
        entry.num_windows_in_original_score = len(all_contributing_rna_window_starts)

    def generate_random_score_distribution(self, score_file_key, N, num_iterations):
        if N == 0: return [] 

        all_scores_for_pdb = self._load_score_file(score_file_key)
        if all_scores_for_pdb is None:
            print(f"Cannot generate random distribution for {score_file_key}: score file not loaded/valid.")
            return []
        
        if N > len(all_scores_for_pdb):
             print(f"Warning: N ({N}) for {score_file_key} is greater than available scores ({len(all_scores_for_pdb)}). Clamping N.")
             N = len(all_scores_for_pdb)
             if N == 0 : return []

        random_sums = []
        for _ in range(num_iterations):
            sampled_scores = random.sample(all_scores_for_pdb, N)
            current_sum = sum(sampled_scores)
            random_sums.append(current_sum)
        return random_sums

class Plotter:
    """Handles creation and saving of plots."""
    def __init__(self, output_dir_base):
        self.output_dir = Path(output_dir_base)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Plotter: Saving plots to {self.output_dir.resolve()}")

    def plot_individual_distribution(self, entry: RRMDataEntry, num_iterations: int):
        if entry.original_score is None or not entry.random_score_distribution:
            print(f"Skipping plot for {entry.score_file_key}: missing original score or random distribution.")
            return

        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Use "uniprot_id_full (pdb_id_chain_id)" for display and filename
        display_label = f"{entry.uniprot_id_full} ({entry.pdb_id}_{entry.chain_id})"
        
        title = (f"Score Distribution for {display_label}\n"
                 f"Original Score (N={entry.num_windows_in_original_score} windows) vs. {num_iterations} Random Samples of N windows")
        ax.set_title(title, fontsize=12)

        ax.hist(entry.random_score_distribution, bins=50, alpha=0.75, label='Random Score Sums', color='skyblue', edgecolor='black')
        ax.axvline(entry.original_score, color='red', linestyle='dashed', linewidth=2, label=f'Original Score: {entry.original_score:.3f}')

        mean_random = np.mean(entry.random_score_distribution)
        std_random = np.std(entry.random_score_distribution)
        ax.axvline(mean_random, color='green', linestyle='dotted', linewidth=2, label=f'Random Mean: {mean_random:.3f}')

        count_ge_original = sum(1 for rs in entry.random_score_distribution if rs >= entry.original_score)
        entry.p_value = (count_ge_original + 1) / (len(entry.random_score_distribution) + 1) 

        if std_random > 1e-9: 
            entry.z_score = (entry.original_score - mean_random) / std_random
        else:
            entry.z_score = 0 if abs(entry.original_score - mean_random) < 1e-9 else float('inf') * np.sign(entry.original_score - mean_random)
            
        stats_text = (f"P-value (orig >= rand): {entry.p_value:.4f}\n"
                      f"Z-score: {entry.z_score:.2f}\n"
                      f"Random Mean: {mean_random:.3f}\n"
                      f"Random StdDev: {std_random:.3f}")
        
        x_text_pos = 0.05 if entry.original_score > mean_random else 0.95 
        ha_text = 'left' if entry.original_score > mean_random else 'right'

        ax.text(x_text_pos, 0.95, stats_text, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', horizontalalignment=ha_text,
                bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))

        ax.set_xlabel(f'Sum of Scores for N={entry.num_windows_in_original_score} Windows')
        ax.set_ylabel('Frequency')
        ax.legend()
        ax.grid(True, linestyle='--', alpha=0.6)
        
        plt.tight_layout()
        # Sanitize display_label for filename
        safe_filename_base = re.sub(r'[^\w\-_.\(\) ]', '_', display_label) # Allow spaces, parens
        safe_filename_base = re.sub(r'\s+', '_', safe_filename_base) # Replace spaces with underscore
        plot_filename = self.output_dir / f"{safe_filename_base}_distribution.png"
        plt.savefig(plot_filename)
        plt.close(fig)

    def plot_p_value_summary(self, entries_list, filename="p_value_summary.png"):
        valid_entries = [e for e in entries_list if e.p_value is not None]
        if not valid_entries:
            print("No valid entries with P-values to plot for summary.")
            return

        sorted_entries = sorted(valid_entries, key=lambda e: e.p_value)
        
        # Use "uniprot_id_full (pdb_id_chain_id)" for labels
        labels = [f"{e.uniprot_id_full} ({e.pdb_id}_{e.chain_id})" for e in sorted_entries]
        p_values = [e.p_value for e in sorted_entries]

        fig, ax = plt.subplots(figsize=(max(12, len(labels) * 0.35), 9)) 
        
        bar_colors = ['red' if p < 0.05 else 'blue' for p in p_values]
        bars = ax.bar(labels, p_values, color=bar_colors, edgecolor='black')
        ax.axhline(y=0.05, color='green', linestyle='--', label='P=0.05 Threshold')

        ax.set_xlabel('UniProt ID (PDB ID_Chain)', fontsize=10) # Updated x-axis label description
        ax.set_ylabel('P-value (Original Score >= Random Score)', fontsize=10)
        ax.set_title(f'P-value Summary ({len(labels)} RRMs, Sorted by P-value)', fontsize=12)
        
        plt.xticks(rotation=90, ha='right', fontsize=max(5, 9 - len(labels)//15)) 
        ax.legend()
        ax.grid(True, axis='y', linestyle='--', alpha=0.6)
        
        upper_y_limit = max(0.07, min(1.0, np.max(p_values) * 1.2 if p_values else 0.07))
        ax.set_ylim(0, upper_y_limit)

        plt.tight_layout(pad=0.5) 
        summary_plot_filename = self.output_dir / filename
        plt.savefig(summary_plot_filename)
        plt.close(fig)
        print(f"Saved P-value summary plot: {summary_plot_filename}")

def main():
    parser = argparse.ArgumentParser(description="RRM Connection Score Analyzer (1024 random window method)")
    parser.add_argument('dataset_path', type=str, help='Path to the extracted_data.txt file.')
    parser.add_argument('--windows_dir', type=str, default='windows',
                        help='Directory containing PDBID_CHAIN.txt score files (default: windows).')
    parser.add_argument('--output_dir', type=str, default='1024_random_distribution',
                        help='Directory to save output plots (default: 1024_random_distribution).')
    parser.add_argument('--iterations', type=int, default=5000,
                        help='Number of random iterations for distribution (default: 5000).')
    parser.add_argument('--window_size_global', type=int, default=5, choices=[3,4,5],
                        help='Global RNA window size (default: 5). Must match score files structure.')

    args = parser.parse_args()

    global WINDOW_SIZE 
    WINDOW_SIZE = args.window_size_global
    print(f"--- RRM Connection Score Analyzer ---")
    print(f"Using global RNA WINDOW_SIZE: {WINDOW_SIZE}")
    print(f"Dataset: {args.dataset_path}")
    print(f"Windows scores directory: {args.windows_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Iterations for random sampling: {args.iterations}")

    dataset_parser = DatasetParser(args.dataset_path)
    rrm_entries = dataset_parser.parse_file()

    if not rrm_entries:
        print("No RRM entries parsed from the dataset. Exiting.")
        return

    try:
        score_calculator = ScoreCalculator(args.windows_dir)
    except FileNotFoundError as e:
        print(f"Critical Error: {e}")
        print("Please ensure the windows scores directory exists and is correctly specified.")
        return

    plotter = Plotter(args.output_dir)
    processed_entries_for_summary = []

    for i, entry in enumerate(rrm_entries):
        # Construct the primary display/file label here to ensure consistency
        primary_label = f"{entry.uniprot_id_full} ({entry.pdb_id}_{entry.chain_id})"
        print(f"\nProcessing entry {i+1}/{len(rrm_entries)}: {primary_label}")


        score_calculator.calculate_original_score(entry)
        
        if entry.num_windows_in_original_score == 0:
            print(f"  Original Score: {entry.original_score:.4f} (N=0 unique RNA windows). Skipping random distribution.")
            entry.p_value = None 
            processed_entries_for_summary.append(entry)
            continue
        
        print(f"  Original Score: {entry.original_score:.4f} (from N={entry.num_windows_in_original_score} unique RNA windows)")

        entry.random_score_distribution = score_calculator.generate_random_score_distribution(
            entry.score_file_key, entry.num_windows_in_original_score, args.iterations
        )

        if not entry.random_score_distribution:
            print(f"  Failed to generate random score distribution for {entry.score_file_key}.")
            entry.p_value = None
            processed_entries_for_summary.append(entry)
            continue
            
        plotter.plot_individual_distribution(entry, args.iterations)
        print(f"  Individual plot saved. P-value: {entry.p_value:.4f}, Z-score: {entry.z_score:.2f}")
        processed_entries_for_summary.append(entry)

    print("\n--- All entries processed ---")
    if processed_entries_for_summary:
        plotter.plot_p_value_summary(processed_entries_for_summary)
    else:
        print("No entries were processed to generate a P-value summary.")

    print(f"\nAnalysis complete. Outputs are in '{Path(args.output_dir).resolve()}'")

if __name__ == "__main__":
    main()