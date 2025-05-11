import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import traceback

def calculate_moving_average(positions, scores, window_size=63):

    if len(positions) < window_size:
        mean_value = np.mean(scores)
        return [mean_value] * len(positions)

    position_to_score = {pos: score for pos, score in zip(positions, scores)}

    moving_avgs = []
    half_window = window_size // 2

    for i, pos in enumerate(positions):
        window_positions = list(range(pos - half_window, pos + half_window + 1))
        window_scores = [position_to_score.get(p, None) for p in window_positions]
        valid_scores = [s for s in window_scores if s is not None]

        if valid_scores:
            moving_avgs.append(np.mean(valid_scores))
        else:
            # Fallback if no scores in window, use the original score at that position
            # This might need adjustment based on desired behavior if all surrounding scores are None
            if i < len(scores):
                 moving_avgs.append(scores[i])
            else: # Should not happen if positions and scores are aligned
                 moving_avgs.append(0)


    return moving_avgs


def prepare_binding_plot_data(entry, scorer, dataset_obj, window_size, protein_cds_dict_ref):
    """
    Runs RRMScorer for both dataset RNA and full protein_cds RNA (if available).
    Returns data necessary for plotting, including scores, positions, and mean scores.
    This function is designed to be called by the GUI or batch processes.
    """
    dataset_rna = entry['rna_sequence'].upper()
    base_uniprot_id = entry['uniprot_id'].split('_')[0]

    # Data for the main plot (either full CDS or dataset if full CDS fails/missing)
    main_plot_windows = []
    main_plot_scores = []
    main_plot_adjusted_scores = []
    main_plot_positions = []
    main_plot_is_synthetic = False

    # Mean scores
    dataset_mean_score = None
    full_cds_exclusive_mean_score = None # This is the "full_mean_score" from the original logic for the full CDS
    
    # Boundary information
    start_boundary = None
    end_boundary = None
    boundary_found = False

    # Attempt to process full protein CDS RNA first
    if protein_cds_dict_ref and base_uniprot_id in protein_cds_dict_ref:
        protein_cds_rna_raw = protein_cds_dict_ref[base_uniprot_id]
        protein_cds_rna = dataset_obj._dna_to_rna(protein_cds_rna_raw).upper()

        protein_cds_result = scorer.run(entry['rrm_id'], protein_cds_rna, window_size)

        if protein_cds_result['success']:
            main_plot_windows = protein_cds_result['windows']
            main_plot_scores = protein_cds_result['scores']
            main_plot_adjusted_scores = protein_cds_result['adjusted_scores']
            main_plot_positions = protein_cds_result['positions']
            main_plot_is_synthetic = False

            if main_plot_adjusted_scores:
                full_cds_exclusive_mean_score = np.mean(main_plot_adjusted_scores)

            # Calculate dataset_mean_score based on scoring the dataset_rna separately
            dataset_result_for_mean = scorer.run(entry['rrm_id'], dataset_rna, window_size)
            if dataset_result_for_mean['success'] and dataset_result_for_mean['adjusted_scores']:
                dataset_mean_score = np.mean(dataset_result_for_mean['adjusted_scores'])
            
            # Determine boundaries
            dataset_in_full_rna_start = protein_cds_rna.find(dataset_rna)
            if dataset_in_full_rna_start != -1:
                start_boundary = dataset_in_full_rna_start
                end_boundary = dataset_in_full_rna_start + len(dataset_rna) - window_size # end of the last window fully within dataset
                boundary_found = True

            plot_data = {
                'entry': entry,
                'windows': main_plot_windows,
                'scores': main_plot_scores,
                'adjusted_scores': main_plot_adjusted_scores,
                'positions': main_plot_positions,
                'is_synthetic': main_plot_is_synthetic,
                'dataset_mean_score': dataset_mean_score,
                'full_cds_exclusive_mean_score': full_cds_exclusive_mean_score,
                'start_boundary': start_boundary,
                'end_boundary': end_boundary,
                'boundary_found': boundary_found,
                'error_message': None,
                'using_protein_cds': True
            }
            # print(f"DEBUG binding_analyzer: Successfully processed protein CDS for {entry['rrm_id']}")
            return plot_data
        else:
            # print(f"DEBUG binding_analyzer: RRMScorer failed on protein_cds_rna for {entry['rrm_id']}. Falling back.")
            pass # Fall through to dataset RNA processing

    # Fallback or if protein_cds_rna was not available/processed
    # print(f"DEBUG binding_analyzer: Processing dataset RNA only for {entry['rrm_id']}")
    dataset_result = scorer.run(entry['rrm_id'], dataset_rna, window_size)
    error_message_for_dataset_fallback = None

    if dataset_result['success']:
        main_plot_windows = dataset_result['windows']
        main_plot_scores = dataset_result['scores']
        main_plot_adjusted_scores = dataset_result['adjusted_scores']
        main_plot_positions = dataset_result['positions']
        main_plot_is_synthetic = False
        if main_plot_adjusted_scores:
            dataset_mean_score = np.mean(main_plot_adjusted_scores)
    else:
        # print(f"DEBUG binding_analyzer: RRMScorer failed on dataset_rna for {entry['rrm_id']}. Generating synthetic.")
        main_plot_windows, main_plot_scores, main_plot_positions = scorer.generate_synthetic_scores(dataset_rna, window_size)
        min_score_syn = min(main_plot_scores) if main_plot_scores else 0
        main_plot_adjusted_scores = [score - min_score_syn for score in main_plot_scores]
        main_plot_is_synthetic = True
        dataset_mean_score = None # Or np.mean if you want mean of synthetic
        error_message_for_dataset_fallback = f"RRMScorer failed: {dataset_result.get('error', 'Unknown error')}. Using synthetic data for dataset."

    # For dataset-only plot, boundaries are effectively start/end of the dataset sequence plot itself
    if main_plot_positions:
        start_boundary = min(main_plot_positions) if main_plot_positions else 0
        # end_boundary is the start of the last window for the dataset_rna
        end_boundary = max(main_plot_positions) if main_plot_positions else len(dataset_rna) - window_size
        boundary_found = True # Boundaries are implicitly the dataset itself
    
    plot_data = {
        'entry': entry,
        'windows': main_plot_windows,
        'scores': main_plot_scores,
        'adjusted_scores': main_plot_adjusted_scores,
        'positions': main_plot_positions,
        'is_synthetic': main_plot_is_synthetic,
        'dataset_mean_score': dataset_mean_score,
        'full_cds_exclusive_mean_score': None, # Not applicable here
        'start_boundary': start_boundary,
        'end_boundary': end_boundary,
        'boundary_found': boundary_found, # True if dataset positions exist
        'error_message': error_message_for_dataset_fallback,
        'using_protein_cds': False
    }
    # print(f"DEBUG binding_analyzer: Processed dataset RNA for {entry['rrm_id']}")
    return plot_data


def _calculate_exclusive_mean(dataset_rna, protein_cds_rna, protein_cds_scores_adj, protein_cds_positions, window_size):
    """
    Helper to calculate mean score of protein_cds_rna excluding the dataset_rna portion.
    This was a method in RRMVisualizerGUI, now a helper here.
    Assumes protein_cds_scores_adj and protein_cds_positions are from the full CDS run.
    """
    try:
        dataset_in_full_start = protein_cds_rna.find(dataset_rna)

        if dataset_in_full_start == -1:
            # print(f"DEBUG binding_analyzer._calculate_exclusive_mean: Dataset RNA not found in protein CDS RNA.")
            return None # Dataset not found within the protein_cds_rna

        dataset_window_start_in_cds = dataset_in_full_start
        dataset_window_end_in_cds_exclusive = dataset_in_full_start + len(dataset_rna) - window_size
        
        exclusive_scores = []
        for i, pos in enumerate(protein_cds_positions):
            if not (dataset_window_start_in_cds <= pos <= dataset_window_end_in_cds_exclusive):
                exclusive_scores.append(protein_cds_scores_adj[i])
            # else:
                # print(f"DEBUG: Pos {pos} is within dataset range [{dataset_window_start_in_cds}, {dataset_window_end_in_cds_exclusive}], excluding from exclusive mean.")


        if exclusive_scores:
            return np.mean(exclusive_scores)
        else:
            return None
    except Exception as e:
        return None


# --- Batch processing for binding scores ---
def batch_process_binding_scores_comparison(dataset, scorer, window_size=5, output_dir="output"):
    """
    Process all entries with full RNA comparison in batch mode.
    Calculate and plot the difference between dataset mean score and full RNA mean score.
    """
    results = []
    differences = []


    os.makedirs(output_dir, exist_ok=True)
    
    if not hasattr(scorer, 'protein_cds_dict') or not scorer.protein_cds_dict:
        pass


    for i, entry in enumerate(dataset.entries):
        rrm_id = entry['rrm_id']
        display_label = f"{rrm_id} (PDB: {entry['pdb_id']}_{entry['chain_id']})"

        # print(f"Processing {i+1}/{len(dataset.entries)}: {display_label}")

        dataset_rna = entry['rna_sequence'].upper()
        dataset_result = scorer.run(rrm_id, dataset_rna, window_size)

        if not dataset_result['success'] or not dataset_result['adjusted_scores']:
            # print(f"  Failed to score dataset RNA or no adjusted scores: {dataset_result.get('error', 'Unknown error')}")
            continue
        
        dataset_mean_score = np.mean(dataset_result['adjusted_scores'])

        base_uniprot_id = entry['uniprot_id'].split('_')[0]
        
        # Check for full RNA specific to this entry
        protein_cds_rna_raw = None
        if hasattr(scorer, 'protein_cds_dict') and scorer.protein_cds_dict and base_uniprot_id in scorer.protein_cds_dict:
             protein_cds_rna_raw = scorer.protein_cds_dict[base_uniprot_id]

        if not protein_cds_rna_raw:
            # print(f"  No full RNA found for {base_uniprot_id} for entry {rrm_id}, skipping full comparison for this entry.")
            results.append({
                'rrm_id': rrm_id,
                'dataset_mean': dataset_mean_score,
                'full_rna_mean': None, # Mark as not available
                'difference': None      # Mark as not available
            })
            continue

        protein_cds_rna = dataset._dna_to_rna(protein_cds_rna_raw).upper()
        protein_cds_result = scorer.run(rrm_id, protein_cds_rna, window_size)

        if not protein_cds_result['success'] or not protein_cds_result['adjusted_scores']:
            # print(f"  Failed to score full RNA or no adjusted scores for {rrm_id}: {protein_cds_result.get('error', 'Unknown error')}")
            results.append({
                'rrm_id': rrm_id,
                'dataset_mean': dataset_mean_score,
                'full_rna_mean': None, # Mark as failed
                'difference': None
            })
            continue

        full_rna_mean_score = np.mean(protein_cds_result['adjusted_scores'])
        difference = dataset_mean_score - full_rna_mean_score

        # print(f"  Dataset Mean: {dataset_mean_score:.2f}")
        # print(f"  Full RNA Mean: {full_rna_mean_score:.2f}")
        # print(f"  Difference: {difference:.2f}")

        current_result = {
            'rrm_id': rrm_id,
            'dataset_mean': dataset_mean_score,
            'full_rna_mean': full_rna_mean_score,
            'difference': difference
        }
        results.append(current_result)
        if difference is not None: # Only add valid differences
            differences.append(difference)

        # Plotting for individual comparison
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.bar(['Dataset RNA', 'Full RNA'], [dataset_mean_score, full_rna_mean_score],
               color=['coral', 'skyblue'], edgecolor='black')
        ax.set_ylabel('Mean Binding Score (Adjusted)')
        ax.set_title(f'Mean Binding Score Comparison for {display_label}') # Use display_label
        ax.text(0, dataset_mean_score + abs(dataset_mean_score*0.05) + 0.01, f'{dataset_mean_score:.2f}', ha='center', va='bottom')
        ax.text(1, full_rna_mean_score + abs(full_rna_mean_score*0.05) + 0.01, f'{full_rna_mean_score:.2f}', ha='center', va='bottom')
        
        max_val_for_text = max(dataset_mean_score, full_rna_mean_score)
        ax.text(0.5, max_val_for_text + abs(max_val_for_text*0.1) + 0.02, # Adjusted y for difference text
                f'Difference: {difference:.2f}', ha='center', va='bottom',
                bbox=dict(facecolor='yellow', alpha=0.3))

        fig.tight_layout()
        output_path = os.path.join(output_dir, f"{rrm_id}_binding_score_comparison.png")
        try:
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
        except Exception as e:
            # print(f"Error saving plot {output_path}: {e}")
            pass
        plt.close(fig)
        # print(f"  Saved comparison plot to {output_path}")

    if results:
        import csv
        output_csv = os.path.join(output_dir, "binding_score_comparison_results.csv")
        fieldnames = ['rrm_id', 'dataset_mean', 'full_rna_mean', 'difference']
        
        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        # print(f"\nSaved summary results to {output_csv}")

        if differences: # Only plot if there are valid differences
            fig_dist, ax_dist = plt.subplots(figsize=(10, 6))
            ax_dist.hist(differences, bins=15, color='skyblue', edgecolor='black', alpha=0.7)
            ax_dist.axvline(x=0, color='red', linestyle='--', linewidth=2, label='No Difference')
            mean_difference = np.mean(differences)
            ax_dist.axvline(x=mean_difference, color='green', linestyle='-', linewidth=2,
                            label=f'Mean Difference: {mean_difference:.2f}')
            ax_dist.set_xlabel('Difference (Dataset Mean - Full RNA Mean)')
            ax_dist.set_ylabel('Frequency')
            ax_dist.set_title('Distribution of Binding Score Differences')
            ax_dist.legend()
            fig_dist.tight_layout()
            dist_output_path = os.path.join(output_dir, "binding_score_difference_distribution.png")
            try:
                fig_dist.savefig(dist_output_path, dpi=300, bbox_inches='tight')
            except Exception as e:
                # print(f"Error saving difference distribution plot: {e}")
                pass
            plt.close(fig_dist)
            # print(f"Saved difference distribution plot to {dist_output_path}")
        # else:
            # print("No valid differences to plot for the distribution.")


    return results

# General batch process
def batch_process_score_generation(dataset, scorer, window_size=5):
    """Process all entries in batch mode to generate scores."""
    results = []

    for i, entry in enumerate(dataset.entries):
        display_label = f"{entry['rrm_id']} (PDB: {entry['pdb_id']}_{entry['chain_id']})"

        result = scorer.run(entry['rrm_id'], entry['rna_sequence'], window_size)

        if result['success']:
            # print(f"  Success! Found {len(result['windows'])} windows with scores")
            results.append({
                'rrm_id': entry['rrm_id'],
                'pdb_id': entry['pdb_id'],
                'chain_id': entry['chain_id'],
                'display_label': display_label,
                'windows': result['windows'],
                'scores': result['scores'],
                'adjusted_scores': result['adjusted_scores'],
                'positions': result['positions'],
                'connections': entry.get('connections', []),
                'success': True
            })
        else:
            # print(f"  Failed: {result.get('error', 'Unknown error')}")
            results.append({
                'rrm_id': entry['rrm_id'],
                'pdb_id': entry['pdb_id'],
                'chain_id': entry['chain_id'],
                'display_label': display_label,
                'error': result.get('error', 'Unknown error'),
                'success': False
            })
    return results