import random
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import traceback

# --- ConnectionAnalyzer Class (Algorithm Update Here) ---
class ConnectionAnalyzer:
    """Class for analyzing RRM-RNA connections and calculating their scores."""

    def __init__(self, scorer, dataset):
        self.scorer = scorer
        self.dataset = dataset
        self.cache = {}

    def identify_connection_blocks(self, connections):
        if not connections:
            return []

        connections = sorted(list(set(connections))) # Ensure unique and sorted
        if not connections:
            return []

        blocks = []
        current_block = [connections[0]]

        for i in range(1, len(connections)):
            if connections[i] - connections[i-1] == 1 or (0.3 <= connections[i] - connections[i-1] <= 0.4):
                current_block.append(connections[i])
            else:
                blocks.append(current_block)
                current_block = [connections[i]]

        blocks.append(current_block) # Add the last block

        return blocks

    def get_rna_positions_for_connection(self, connection, rna_seq_length):
        rna_start = round(connection * 3)
        rna_end = rna_start + 2

        if rna_end >= rna_seq_length:
            # print(f"DEBUG: RNA positions for connection {connection} out of bounds ({rna_seq_length}).")
            return None

        return (rna_start, rna_end)

    def get_rna_windows_for_block(self, block, rna_seq_length, window_size=5):
        windows = [] # Stores (start_pos, end_pos) of RNA windows
        all_rna_indices_in_block = [] # Collect all individual RNA nts covered by connections in this block
        for conn_rrm_idx in block:
            # Convert RRM connection index to RNA triplet start/end
            rna_triplet_pos = self.get_rna_positions_for_connection(conn_rrm_idx, rna_seq_length)
            if rna_triplet_pos:
                all_rna_indices_in_block.extend(range(rna_triplet_pos[0], rna_triplet_pos[1] + 1))
        
        if not all_rna_indices_in_block:
            return []

        all_rna_indices_in_block = sorted(list(set(all_rna_indices_in_block))) # Unique, sorted RNA indices

        min_rna_idx_in_block = all_rna_indices_in_block[0]
        max_rna_idx_in_block = all_rna_indices_in_block[-1]

        
        possible_starts = range(max(0, min_rna_idx_in_block - window_size + 1), 
                                min(rna_seq_length - window_size + 1, max_rna_idx_in_block + 1))

        generated_windows_starts = set()
        for s_start in possible_starts:
            s_end = s_start + window_size - 1
            # Check for overlap
            window_indices_set = set(range(s_start, s_end + 1))
            if not window_indices_set.isdisjoint(all_rna_indices_in_block):
                if s_start not in generated_windows_starts:
                    windows.append((s_start, s_end)) # Storing start and end of the window
                    generated_windows_starts.add(s_start)
        
        return sorted(list(set(windows)))


    def cache_all_window_scores(self, rrm_id, rna_sequence, window_size=5):
        # More robust key, but careful with very long sequences if memory is an issue
        cache_key = f"{rrm_id}_{rna_sequence[:20]}_{rna_sequence[-20:]}_{len(rna_sequence)}_{window_size}" 

        if cache_key in self.cache:
            return self.cache[cache_key]

        result = self.scorer.run(rrm_id, rna_sequence, window_size)

        if not result['success'] or not result.get('adjusted_scores'): # Check for adjusted_scores specifically
            # print(f"Warning: Failed to get scores for {rrm_id} on sequence. RRMScorer error: {result.get('error')}")
            windows, scores, positions = self.scorer.generate_synthetic_scores(rna_sequence, window_size)
            if not scores:
                # print(f"ERROR: Could not get or generate scores for {rrm_id}. Returning empty dict.")
                self.cache[cache_key] = {} # Cache empty result to avoid re-computation
                return {} 
            # print(f"Warning: Using synthetic scores for {rrm_id}")
            min_score = min(scores) if scores else 0
            adjusted_scores = [score - min_score for score in scores]
            # Use positions from synthetic scores
        else:
            adjusted_scores = result['adjusted_scores']
            positions = result['positions']
            # windows = result['windows'] # Not directly used for the position_to_score map

        position_to_score = {}
        for pos, adj_score in zip(positions, adjusted_scores):
            if pos not in position_to_score:
                 position_to_score[pos] = adj_score

        self.cache[cache_key] = position_to_score
        return position_to_score

    def find_subsequence_start(self, full_sequence, subsequence):
        if isinstance(full_sequence, list): # Should not happen if RNA strings are used
            full_sequence = ''.join(str(item) for item in full_sequence)
        if isinstance(subsequence, list):
            subsequence = ''.join(str(item) for item in subsequence)

        pos = full_sequence.find(subsequence)
        if pos != -1:
            return pos

        for length in [50, 40, 30, 20, 15, 10]:
             if len(subsequence) > length and len(full_sequence) > length : # Ensure full_sequence is also long enough
                  prefix = subsequence[:length]
                  pos_prefix = full_sequence.find(prefix)
                  if pos_prefix != -1:
                      # print(f"DEBUG: Found subsequence using prefix match (len={length}) at {pos_prefix}")
                      return pos_prefix 

                  suffix = subsequence[-length:]
                  pos_suffix = full_sequence.rfind(suffix) # Use rfind for suffix to get last occurrence
                  if pos_suffix != -1:
                      start_pos_from_suffix = pos_suffix - (len(subsequence) - length)
                      # print(f"DEBUG: Found subsequence using suffix match (len={length}) at {pos_suffix}, calculated start {start_pos_from_suffix}")
                      if start_pos_from_suffix >= 0:
                           return start_pos_from_suffix
        # print(f"DEBUG: Subsequence not found in full sequence.")
        return -1


    def calculate_connection_score(self, rrm_id, rna_sequence, connections, window_size=5):
        """
        Calculates the sum of adjusted scores for RNA windows corresponding to RRM connections.
        'connections' are RRM indices.
        """
        position_to_score = self.cache_all_window_scores(rrm_id, rna_sequence, window_size)

        if not position_to_score:
            # print(f"Warning: No scores available for {rrm_id} on sequence '{rna_sequence[:20]}...'. Returning score 0.")
            return 0.0, [], []

        rrm_blocks = self.identify_connection_blocks(connections) # These are blocks of RRM indices
        if not rrm_blocks:
             return 0.0, [], []

        total_block_set_score = 0.0
        all_scored_rna_window_starts_for_blocks = set()
        
        # This will store (start_rna_pos, end_rna_pos) for each window identified for any block
        all_rna_windows_coords_for_all_blocks = [] 

        for rrm_block in rrm_blocks: # rrm_block is a list of RRM indices
            # Get all unique RNA scoring windows (start_pos, end_pos) that overlap this RRM block
            rna_windows_for_this_block = self.get_rna_windows_for_block(rrm_block, len(rna_sequence), window_size)
            all_rna_windows_coords_for_all_blocks.extend(rna_windows_for_this_block)
            
            current_block_score = 0.0
            for rna_win_start, rna_win_end in rna_windows_for_this_block:
                # Score is associated with the start position of the RNA window
                if rna_win_start not in all_scored_rna_window_starts_for_blocks: # Ensure each window scored once across all blocks
                    score_val = position_to_score.get(rna_win_start)
                    if score_val is not None:
                        current_block_score += score_val
                        all_scored_rna_window_starts_for_blocks.add(rna_win_start)
                    # else:
                        # print(f"DEBUG: No score found for RNA window starting at {rna_win_start} for RRM block {rrm_block}")
            total_block_set_score += current_block_score
        
        block_sizes = [len(b) for b in rrm_blocks] # Sizes of RRM blocks

        # Return unique sorted list of (start, end) RNA window coordinates involved
        unique_rna_windows_coords = sorted(list(set(all_rna_windows_coords_for_all_blocks)))

        return total_block_set_score, block_sizes, unique_rna_windows_coords

    def generate_random_connections(self, block_sizes, rna_seq_length, min_gap_rrm=1):
        """
        Generate random non-overlapping RRM connection positions based on given block_sizes.
        Returns a sorted list of RRM indices.
        min_gap_rrm is the minimum number of RRM units between blocks.
        """
        n_rrm_slots_available = rna_seq_length // 3
        if n_rrm_slots_available <= 0:
            # print("Warning: RNA sequence too short for any RRM positions in generate_random_connections.")
            return []

        if not block_sizes:
            return []

        k_num_blocks = len(block_sizes)
        
        # Content length: sum of sizes of all blocks
        total_rrm_units_in_blocks = sum(block_sizes)
        
        min_rrm_slots_for_gaps = max(0, k_num_blocks - 1) * min_gap_rrm
        min_total_rrm_slots_needed = total_rrm_units_in_blocks + min_rrm_slots_for_gaps

        if min_total_rrm_slots_needed > n_rrm_slots_available:
            # print(f"Warning: Not enough RRM slots ({n_rrm_slots_available}) to place blocks (total len {total_rrm_units_in_blocks}) with gaps (total {min_rrm_slots_for_gaps}). Required: {min_total_rrm_slots_needed}. Trying with 0 gap.")
            min_gap_rrm = 0 # Try with no gap
            min_rrm_slots_for_gaps = 0
            min_total_rrm_slots_needed = total_rrm_units_in_blocks
            if min_total_rrm_slots_needed > n_rrm_slots_available:
                # print("Error: Cannot even place blocks without gaps.")
                return []

        n_distributable_extra_rrm_slots = n_rrm_slots_available - min_total_rrm_slots_needed
        
        if k_num_blocks == 0: # No blocks, no connections
            return []
        
        bar_positions = sorted(random.sample(range(n_distributable_extra_rrm_slots + k_num_blocks), k_num_blocks))
        
        # Calculate sizes of each of the k_num_blocks+1 gap segments
        gap_segment_sizes = []
        last_bar_pos = -1
        for bar_pos in bar_positions:
            gap_segment_sizes.append(bar_pos - last_bar_pos - 1)
            last_bar_pos = bar_pos
        gap_segment_sizes.append((n_distributable_extra_rrm_slots + k_num_blocks -1) - last_bar_pos)

        # Construct the final list of RRM connections
        random_connections_list = []
        current_rrm_position = gap_segment_sizes[0] # Start after the first gap segment
        
        # Shuffle block_sizes so their placement order is random
        shuffled_block_sizes = list(block_sizes)
        random.shuffle(shuffled_block_sizes)

        for i in range(k_num_blocks):
            block_len = shuffled_block_sizes[i]
            # Add RRM indices for the current block
            for rrm_idx_in_block in range(block_len):
                if current_rrm_position + rrm_idx_in_block < n_rrm_slots_available:
                     random_connections_list.append(current_rrm_position + rrm_idx_in_block)
            
            current_rrm_position += block_len # Move past the current block
            if i < k_num_blocks - 1: # If not the last block, add mandatory gap and then the next gap segment
                current_rrm_position += min_gap_rrm + gap_segment_sizes[i+1]
        
        return sorted(random_connections_list)


    def calculate_random_distribution(self, rrm_id, rna_sequence_context, original_connections_rrm, n_iterations=1000, window_size=5):
        """
        Calculate the original connection score and a distribution of random scores.
        original_connections_rrm: RRM indices relative to the start of rna_sequence_context.
        rna_sequence_context: The RNA sequence string in which connections are defined and random placements occur.
        """
        # print(f"\nDEBUG ConnectionAnalyzer: Calculating distribution for RRM ID: {rrm_id}")
        # print(f"DEBUG ConnectionAnalyzer: Context RNA Length: {len(rna_sequence_context)}, Window Size: {window_size}")
        # print(f"DEBUG ConnectionAnalyzer: Input Original RRM Connections: {original_connections_rrm}")

        orig_score, orig_block_rrm_sizes, block_rna_window_coords = self.calculate_connection_score(
            rrm_id, rna_sequence_context, original_connections_rrm, window_size
        )
        # print(f"DEBUG ConnectionAnalyzer: Original RRM block sizes: {orig_block_rrm_sizes}")
        # print(f"DEBUG ConnectionAnalyzer: Original connection score: {orig_score:.4f}")
        # print(f"DEBUG ConnectionAnalyzer: RNA windows for original blocks: {block_rna_window_coords}")


        if not orig_block_rrm_sizes: # No blocks means nothing to randomize based on
             # print("Warning: No RRM connection blocks identified from original connections. Cannot generate random distribution.")
             return orig_score, [], block_rna_window_coords


        # --- Calculate Random Distribution ---
        random_scores = []
        
        # Pre-cache scores for the rna_sequence_context
        _ = self.cache_all_window_scores(rrm_id, rna_sequence_context, window_size)

        # print(f"DEBUG ConnectionAnalyzer: Generating {n_iterations} random placements within {len(rna_sequence_context)} nt RNA ({len(rna_sequence_context) // 3} RRM slots)...")
        for i in range(n_iterations):
            # Generate random RRM connection positions based on original RRM block sizes
            # within the provided rna_sequence_context length
            random_rrm_connections = self.generate_random_connections(orig_block_rrm_sizes, len(rna_sequence_context))

            if not random_rrm_connections: # Should ideally not happen if orig_block_rrm_sizes is valid
                random_scores.append(0.0) 
                continue

            # Calculate score for this random placement using the same rna_sequence_context
            rand_score, _, _ = self.calculate_connection_score(
                rrm_id, rna_sequence_context, random_rrm_connections, window_size
            )
            random_scores.append(rand_score)


        return orig_score, random_scores, block_rna_window_coords


# --- ConnectionDistributionPlot Class ---
class ConnectionDistributionPlot:
    """Class for creating distribution plots of connection scores."""

    def __init__(self, figure, ax, canvas): # canvas can be None for batch plotting
        self.figure = figure
        self.ax = ax
        self.canvas = canvas

    def plot_distribution(self, original_score, random_scores, entry_info, iterations=1000):
        """
        Create a distribution plot of random scores with the original score highlighted.
        Calculates and displays a consistent one-sided (right-tailed) p-value.
        entry_info can be an entry dict or a string RRM ID.
        """
        if self.ax is None or self.figure is None : # If used in batch and fig/ax not pre-created
            self.figure, self.ax = plt.subplots(figsize=(10,7)) # Default size

        self.figure.clear() # Clear if reusing figure object
        self.ax = self.figure.add_subplot(111)


        entry_label = "Unknown Entry"
        if isinstance(entry_info, dict): # If full entry object
            entry_label = f'{entry_info.get("rrm_id", "N/A")} (PDB: {entry_info.get("pdb_id", "N/A")}_{entry_info.get("chain_id", "N/A")})'
        elif isinstance(entry_info, str): # If just RRM ID string
            entry_label = entry_info
        
        title = f'Connection Score Distribution for {entry_label}\n' \
                f'({iterations} Random Shuffles vs. Original)'
        self.ax.set_title(title, fontsize=12)


        if not random_scores or len(random_scores) == 0:
            self.ax.text(0.5, 0.5, "No random scores generated to plot.", ha='center', va='center', fontsize=12)
            # print(f"Warning: No random scores provided for {entry_label} in plot_distribution. Plotting original score only.")
            if original_score is not None:
                 self.ax.axvline(x=original_score, color='red', linestyle='--', linewidth=2,
                             label=f'Original: {original_score:.2f}')
                 self.ax.legend()
            self.ax.set_xlabel('Sum of Scores for Connection Block Windows')
            self.ax.set_ylabel('Frequency')
            if self.canvas: 
                 self.canvas.draw()
            return 

        n_total_random = len(random_scores)
        mean_random = np.mean(random_scores)
        std_random = np.std(random_scores)
        
        # Percentile: % of random scores <= original_score
        percentile = sum(1 for x in random_scores if x <= original_score) / n_total_random * 100 if n_total_random > 0 else 0

        z_score = (original_score - mean_random) / std_random if std_random > 0 else 0
        
        # P-value (one-sided, right-tailed): (count(random >= original) + 1) / (N_random + 1)
        # This assesses if the original score is significantly HIGHER than random.
        count_equal_or_greater = sum(1 for x in random_scores if x >= original_score)
        p_value_right_tailed = (count_equal_or_greater + 1) / (n_total_random + 1) if n_total_random > 0 else 1.0


        self.ax.hist(random_scores, bins=min(30, max(10, n_total_random // 20 if n_total_random > 20 else 10)), # Dynamic bins
                     alpha=0.75, color='skyblue', edgecolor='black', linewidth=0.5, density=False) # density=False for frequency

        self.ax.axvline(x=original_score, color='red', linestyle='--', linewidth=2,
                        label=f'Original: {original_score:.2f}')

        self.ax.axvline(x=mean_random, color='darkblue', linestyle='-', linewidth=1.5, alpha=0.8,
                        label=f'Random Mean: {mean_random:.2f} (SD: {std_random:.2f})')

        self.ax.set_xlabel('Sum of Scores for Connection Block Windows')
        self.ax.set_ylabel('Frequency')
        
        legend_text = (
            f'Original: {original_score:.2f}\n'
            f'Random Mean: {mean_random:.2f} (SD: {std_random:.2f})\n'
            f'Z-score: {z_score:.2f}\n'
            f'P-value (orig >= random): {p_value_right_tailed:.4f}\n' 
            f'Percentile (orig >= %random): {percentile:.1f}%'
        )
        # Adding legend outside, then annotation box inside for clarity
        handles, labels = self.ax.get_legend_handles_labels()
        # self.ax.legend(handles, labels, loc='upper right', bbox_to_anchor=(1, 0.9)) # Simple legend
        
        hist_max_y = self.ax.get_ylim()[1]
        
        # Decide horizontal alignment based on original score position relative to mean
        x_annot_pos = 0.03
        ha_annot = 'left'
        if original_score > mean_random and (original_score - mean_random) > 0.5 * std_random : # If original is significantly to the right
             pass # Default left is fine
        elif mean_random > original_score and (mean_random - original_score) > 0.5 * std_random: # If original is significantly to the left
             x_annot_pos = 0.97
             ha_annot = 'right'


        self.ax.annotate(legend_text,
                         xy=(x_annot_pos, 0.95), xycoords='axes fraction',
                         horizontalalignment=ha_annot,
                         verticalalignment='top',
                         fontsize=8, 
                         bbox=dict(boxstyle="round,pad=0.4", fc="ivory", alpha=0.7, ec="gray"))


        self.ax.grid(axis='y', linestyle='--', alpha=0.5) 

        self.figure.tight_layout(pad=1.5) 
        if self.canvas: 
            self.canvas.draw()

# --- Helper functions for connection analysis that were in main or global scope ---
def random_normal_bounded(mean=1.0, std_dev=0.05, min_val=0.9, max_val=1.1):
    # import numpy as np # Already imported
    while True:
        value = np.random.normal(mean, std_dev)
        if min_val <= value <= max_val:
            return value

def calculate_and_display_pvalue(original_score, random_scores): # Renamed for clarity if used outside plotting
    """
    Calculates p-value. Original script had a slightly different logic for display.
    This version will return right-tailed and left-tailed for flexibility.
    """
    if not random_scores: return None, None
    
    n = len(random_scores)
    p_greater_equal = (sum(1 for x in random_scores if x >= original_score) + 1) / (n + 1)

    p_less_equal = (sum(1 for x in random_scores if x <= original_score) + 1) / (n + 1)
    
    return p_greater_equal, p_less_equal


# --- Batch processing for connection analysis ---
def batch_connection_analysis_processing(dataset, scorer, conn_analyzer_instance, iterations, window_size, output_dir, generate_plots=False):
    """
    Performs connection analysis for entries in batch mode.
    Saves results to CSV and optionally generates individual distribution plots.
    """
    entries_with_connections = [entry for entry in dataset.entries if entry.get('connections')]
    if not entries_with_connections:
        # print("No entries with connection data found for batch connection analysis.")
        return []

    # print(f"\nRunning Batch Connection Analysis for {len(entries_with_connections)} entries...")
    connection_results_summary = []

    for i, entry in enumerate(entries_with_connections):
        rrm_id = entry['rrm_id']
        display_label = f"{entry['rrm_id']} (PDB: {entry['pdb_id']}_{entry['chain_id']})"

        context_rna = entry['rna_sequence'].upper() # Default to dataset RNA
        connections_for_analysis = entry['connections'] # Assumed relative to context_rna

        base_uniprot_id = entry['uniprot_id'].split('_')[0]
        protein_cds_dict = getattr(scorer, 'protein_cds_dict', None)

        # This logic mirrors what might be in GUI's _run_connection_analysis for context determination
        if protein_cds_dict and base_uniprot_id in protein_cds_dict:
            full_cds_rna_raw = protein_cds_dict[base_uniprot_id]
            if full_cds_rna_raw:
                full_cds_rna = dataset._dna_to_rna(full_cds_rna_raw).upper()
                dataset_rna_original = entry['rna_sequence'].upper()
                
                start_pos_rna_in_full = conn_analyzer_instance.find_subsequence_start(full_cds_rna, dataset_rna_original)

                if start_pos_rna_in_full != -1:
                    context_rna = full_cds_rna
                    rrm_shift = start_pos_rna_in_full // 3
                    connections_for_analysis = [c + rrm_shift for c in entry['connections']]


        try:
            # Pass the determined context_rna and connections_for_analysis
            orig_score, random_scores, _ = conn_analyzer_instance.calculate_random_distribution(
                rrm_id, context_rna, connections_for_analysis, iterations, window_size
            )

            current_entry_results = {
                'rrm_id': rrm_id,
                'display_label': display_label,
                'original_score': orig_score, # Using the directly calculated original score
                'random_mean': None, 'random_std': None, 'z_score': None,
                'p_value_greater_equal': None, 'p_value_less_equal': None, 'percentile': None
            }

            if random_scores:
                n_total = len(random_scores)
                mean_random = np.mean(random_scores)
                std_random = np.std(random_scores)
                percentile = sum(1 for x in random_scores if x <= orig_score) / n_total * 100 if n_total > 0 else 0.0
                z_score = (orig_score - mean_random) / std_random if std_random > 0 else 0.0
                
                p_val_ge, p_val_le = calculate_and_display_pvalue(orig_score, random_scores)

                current_entry_results.update({
                    'random_mean': mean_random, 'random_std': std_random, 'z_score': z_score,
                    'p_value_greater_equal': p_val_ge, 'p_value_less_equal': p_val_le, 'percentile': percentile
                })
                
                # print(f"  Original Score: {orig_score:.4f}") # Using orig_score
                # print(f"  Random Mean: {mean_random:.4f}, SD: {std_random:.4f}")
                # print(f"  Z-score: {z_score:.2f}, P-val (orig>=rnd): {p_val_ge:.4f}, P-val (orig<=rnd): {p_val_le:.4f}, Percentile: {percentile:.1f}%")

                if generate_plots:
                    plot_output_path = os.path.join(output_dir, f"{entry['rrm_id'].replace('/', '_')}_connection_distribution.png")
                    # Create a new figure and axes for each plot to avoid issues
                    fig_dist, ax_dist = plt.subplots(figsize=(10, 6))
                    dist_plotter = ConnectionDistributionPlot(fig_dist, ax_dist, None) # No canvas for batch
                    dist_plotter.plot_distribution(orig_score, random_scores, entry, iterations) # Pass full entry
                    try:
                        fig_dist.savefig(plot_output_path, dpi=300, bbox_inches='tight')
                    except Exception as e_save:
                        pass
                    plt.close(fig_dist) # Close the figure
            
            connection_results_summary.append(current_entry_results)

        except Exception as e_proc:
            # print(f"  ERROR processing {rrm_id} for connection analysis: {str(e_proc)}")
            # print(traceback.format_exc())
            connection_results_summary.append({
                'rrm_id': rrm_id, 'display_label': display_label, 'original_score': None, 
                'error': str(e_proc)
            })
    
    return connection_results_summary


def create_additional_summary_plots(connection_results_summary, output_dir):
    """
    Generates summary plots from batch connection analysis results.
    connection_results_summary is a list of dicts from batch_connection_analysis_processing.
    """
    filtered_results = [r for r in connection_results_summary 
                        if r.get('original_score') is not None and 
                           r.get('random_mean') is not None and # Ensure random stats were computed
                           r.get('p_value_greater_equal') is not None and 
                           r.get('z_score') is not None]


    if not filtered_results:
        # print("No valid entries found for additional connection summary plots after filtering.")
        return

    # Sort by original_score for the first plot
    sorted_by_orig_score = sorted(filtered_results, key=lambda x: x['original_score'])
    rrm_ids_orig_sort = [r['display_label'] for r in sorted_by_orig_score] # Use display_label
    original_scores = [r['original_score'] for r in sorted_by_orig_score]
    random_means = [r['random_mean'] for r in sorted_by_orig_score]

    # Using p_value_greater_equal as the primary p-value for plotting significance
    p_values_ge = [r['p_value_greater_equal'] for r in sorted_by_orig_score] 
    z_scores = [r['z_score'] for r in sorted_by_orig_score]


    # 1. Original vs Mean Shuffle Values Plot
    if len(rrm_ids_orig_sort) > 0:
        fig1, ax1 = plt.subplots(figsize=(max(14, len(rrm_ids_orig_sort)*0.5), 8)) # Dynamic width
        x_indices = np.arange(len(rrm_ids_orig_sort))
        width = 0.35
        ax1.bar(x_indices - width/2, original_scores, width, label='Original Scores', color='coral', edgecolor='black')
        ax1.bar(x_indices + width/2, random_means, width, label='Mean Shuffle Scores', color='skyblue', edgecolor='black')
        ax1.set_xlabel('RRM Entry')
        ax1.set_ylabel('Connection Score Sum')
        ax1.set_title('Original Scores vs Mean Shuffle Scores', fontsize=14)
        ax1.set_xticks(x_indices)
        ax1.set_xticklabels(rrm_ids_orig_sort, rotation=90, ha="right", fontsize=max(6, 10 - len(rrm_ids_orig_sort)//10))
        ax1.legend()
        ax1.grid(axis='y', linestyle='--', alpha=0.7)
        fig1.tight_layout(pad=1.5)
        output_path1 = os.path.join(output_dir, "conn_original_vs_mean_scores.png")
        try:
            fig1.savefig(output_path1, dpi=300, bbox_inches='tight')
        except Exception as e: # print(f"Error saving plot {output_path1}: {e}")
            pass
        plt.close(fig1)
        # print(f"Saved original vs mean connection scores plot to {output_path1}")

    # 2. P-values Plot (p_value_greater_equal)
    # Sort by p-value for this plot
    p_val_sorted_results = sorted(filtered_results, key=lambda x: x['p_value_greater_equal'])
    p_rrm_ids = [r['display_label'] for r in p_val_sorted_results]
    p_values_plot = [r['p_value_greater_equal'] for r in p_val_sorted_results]

    if len(p_rrm_ids) > 0:
        fig2, ax2 = plt.subplots(figsize=(max(14, len(p_rrm_ids)*0.5), 8))
        bars_p = ax2.bar(range(len(p_rrm_ids)), p_values_plot, color='lightgreen', edgecolor='black')
        for i, bar in enumerate(bars_p):
            if p_values_plot[i] <= 0.05: # Highlight significant p-values
                bar.set_color('salmon')
        ax2.set_xlabel('RRM Entry')
        ax2.set_ylabel('P-value (Original >= Random)')
        ax2.set_title('Connection Score P-values (orig >= random)', fontsize=14)
        ax2.set_xticks(range(len(p_rrm_ids)))
        ax2.set_xticklabels(p_rrm_ids, rotation=90, ha="right", fontsize=max(6, 10 - len(p_rrm_ids)//10))
        ax2.axhline(y=0.05, color='blue', linestyle='--', alpha=0.7, label='Significance Threshold (p=0.05)')
        ax2.legend()
        ax2.grid(axis='y', linestyle='--', alpha=0.7)
        ax2.set_ylim(0, max(0.1, min(1.0, np.percentile(p_values_plot, 99) * 1.2 if p_values_plot else 0.1))) # Dynamic y-axis for p-values
        fig2.tight_layout(pad=1.5)
        output_path2 = os.path.join(output_dir, "conn_p_values_summary.png")
        try:
            fig2.savefig(output_path2, dpi=300, bbox_inches='tight')
        except Exception as e: # print(f"Error saving plot {output_path2}: {e}")
            pass
        plt.close(fig2)
        # print(f"Saved connection p-values plot to {output_path2}")

    # 3. Z-scores Plot
    # Sort by Z-score for this plot (e.g., highest Z-score first)
    z_score_sorted_results = sorted(filtered_results, key=lambda x: x['z_score'], reverse=True)
    z_rrm_ids = [r['display_label'] for r in z_score_sorted_results]
    z_scores_plot = [r['z_score'] for r in z_score_sorted_results]
    
    if len(z_rrm_ids) > 0:
        mean_z_val = np.mean(z_scores_plot) if z_scores_plot else 0
        fig3, ax3 = plt.subplots(figsize=(max(14, len(z_rrm_ids)*0.5), 8))
        bars_z = ax3.bar(range(len(z_rrm_ids)), z_scores_plot, color='lightblue', edgecolor='black')
        for i, bar in enumerate(bars_z): # Highlight significant Z-scores
            if abs(z_scores_plot[i]) >= 1.96: # Corresponds to p=0.05 two-tailed
                bar.set_color('salmon')
        ax3.set_xlabel('RRM Entry')
        ax3.set_ylabel('Z-score')
        ax3.set_title('Connection Score Z-scores', fontsize=14)
        ax3.set_xticks(range(len(z_rrm_ids)))
        ax3.set_xticklabels(z_rrm_ids, rotation=90, ha="right", fontsize=max(6, 10 - len(z_rrm_ids)//10))
        ax3.axhline(y=1.96, color='blue', linestyle='--', alpha=0.7, label='Significance Threshold (|Z| >= 1.96)')
        ax3.axhline(y=-1.96, color='blue', linestyle='--', alpha=0.7)
        ax3.axhline(y=mean_z_val, color='green', linestyle='-', alpha=0.7, label=f'Mean Z-score ({mean_z_val:.2f})')
        ax3.legend()
        ax3.grid(axis='y', linestyle='--', alpha=0.7)
        fig3.tight_layout(pad=1.5)
        output_path3 = os.path.join(output_dir, "conn_z_scores_summary.png")
        try:
            fig3.savefig(output_path3, dpi=300, bbox_inches='tight')
        except Exception as e: # print(f"Error saving plot {output_path3}: {e}")
            pass
        plt.close(fig3)
        # print(f"Saved connection Z-scores plot to {output_path3}")