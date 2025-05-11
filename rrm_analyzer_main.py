import sys
import os
import argparse
import traceback
import tkinter as tk
from core_logic import RRMDataset, RRMScorer
import binding_analyzer
import connection_analyzer
import gui


def main():
    parser = argparse.ArgumentParser(description='RRM Binding Score Analyzer')
    parser.add_argument('dataset_path', help='Path to the dataset file')
    parser.add_argument('output_dir', help='Directory for output files')
    parser.add_argument('--window-size', type=int, default=5, choices=[3, 5], help='Window size (3 or 5)')
    parser.add_argument('--gui', action='store_true', help='Run in interactive GUI mode (default if no batch op specified)')
    
    # Batch operation flags
    parser.add_argument('--batch-score-generation', action='store_true',
                        help='Batch: Generate and save RRMScorer scores for all dataset entries.')
    parser.add_argument('--batch-binding-comparison', action='store_true',
                        help='Batch: Generate comparison plots of dataset vs full RNA binding scores.')
    parser.add_argument('--batch-connection-analysis', action='store_true',
                        help='Batch: Run connection block analysis, save results and summary plots.')
    parser.add_argument('--batch-connection-dist-plots', action='store_true',
                        help='Batch (requires --batch-connection-analysis): Also generate individual distribution plots for each RRM during connection analysis.')

    # Parameters for batch connection analysis
    parser.add_argument('--iterations', type=int, default=1000,
                      help='Number of random iterations for connection analysis (default: 1000)')
    # parser.add_argument('--highlight-connections', action='store_true', # This is a GUI feature primarily
    #                   help='Highlight connection positions in plots (batch mode only where applicable)')


    args = parser.parse_args()

    print("RRM Binding Score Analyzer")
    print(f"Dataset: {args.dataset_path}")
    print(f"Output directory: {args.output_dir}")
    print(f"Window size: {args.window_size}")

    if not os.path.exists(args.dataset_path):
        print(f"Error: Dataset file not found: {args.dataset_path}")
        return 1

    # Check for rrm_rna_wrapper.py early, as RRMScorer depends on it.
    wrapper_path = "rrm_rna_wrapper.py" # RRMScorer looks for this itself
    if not os.path.exists(wrapper_path):
        print(f"Error: RRMScorer wrapper script not found: {wrapper_path}")
        print(f"Make sure '{wrapper_path}' exists in the current directory.")
        return 1

    os.makedirs(args.output_dir, exist_ok=True)

    try:
        dataset = RRMDataset(args.dataset_path)
        if not dataset.entries:
            print("Error: No valid entries found in dataset")
            return 1

        scorer = RRMScorer(args.output_dir) # Initializes RRMScorer
        if not scorer.available:
            print("Error: RRMScorer is not available or not working. Please check its setup.")

        scorer.entries = dataset.entries # For _generate_complete_windows in RRMScorer

        # Determine if any batch operation is requested
        is_batch_mode = args.batch_score_generation or \
                        args.batch_binding_comparison or \
                        args.batch_connection_analysis

        needs_protein_cds = args.batch_binding_comparison or \
                              args.batch_connection_analysis or \
                              args.gui

        if needs_protein_cds:
            if os.path.exists("proten_cds.fasta"):
                 scorer.protein_cds_dict = scorer.load_protein_cds_fasta("proten_cds.fasta")
                 if not scorer.protein_cds_dict:
                      print("Warning: protein_cds.fasta found but could not be loaded or is empty. Full sequence analysis might be affected.")
            else:
                 print("Warning: proten_cds.fasta not found. Full sequence analysis will not be available.")


        if is_batch_mode:
            print("Running in Batch Mode...")
            if args.batch_score_generation:
                print("\n--- Starting Batch Score Generation ---")
                binding_analyzer.batch_process_score_generation(dataset, scorer, args.window_size)
                print("--- Finished Batch Score Generation ---")

            if args.batch_binding_comparison:
                if not hasattr(scorer, 'protein_cds_dict') or not scorer.protein_cds_dict:
                    print("Error: protein_cds.fasta is required for --batch-binding-comparison but was not loaded successfully.")
                else:
                    print("\n--- Starting Batch Binding Score Comparison ---")
                    binding_analyzer.batch_process_binding_scores_comparison(
                        dataset, scorer, args.window_size, args.output_dir
                    )
                    print("--- Finished Batch Binding Score Comparison ---")

            if args.batch_connection_analysis:
                print("\n--- Starting Batch Connection Analysis ---")
                # ConnectionAnalyzer needs scorer and dataset
                conn_analyzer_instance = connection_analyzer.ConnectionAnalyzer(scorer, dataset)
                
                connection_summary = connection_analyzer.batch_connection_analysis_processing(
                    dataset, scorer, conn_analyzer_instance,
                    args.iterations, args.window_size, args.output_dir,
                    generate_plots=args.batch_connection_dist_plots # Pass flag for individual plots
                )
                if connection_summary: # If any results were processed
                    connection_analyzer.create_additional_summary_plots(connection_summary, args.output_dir)
                print("--- Finished Batch Connection Analysis ---")
            
            print("\nBatch processing complete.")

        elif args.gui: # Explicit GUI mode
            print("Launching Interactive GUI Mode...")
            root = tk.Tk()
            # Pass the imported modules themselves to the GUI
            app = gui.RRMVisualizerGUI(root, dataset, scorer, 
                                       binding_analyzer, 
                                       connection_analyzer)
            root.mainloop()
        else:
            print("No specific batch operation selected. Defaulting to Interactive GUI Mode...")
            root = tk.Tk()
            app = gui.RRMVisualizerGUI(root, dataset, scorer, 
                                       binding_analyzer, 
                                       connection_analyzer)
            root.mainloop()

    except Exception as e:
        print(f"An unexpected error occurred in main: {str(e)}")
        print(traceback.format_exc())
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())