import os
import numpy as np
import sys 
import traceback
import threading
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.lines import Line2D 
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

class RRMVisualizerGUI:
    """GUI for RRM score visualization."""

    def __init__(self, root, dataset, scorer, binding_analyzer_module, connection_analyzer_module):
        self.root = root
        self.dataset = dataset 
        self.scorer = scorer   
        self.binding_analyzer = binding_analyzer_module
        self.connection_analyzer_module = connection_analyzer_module 

        self.cache = {} 
        self.conn_analyzer_instance = self.connection_analyzer_module.ConnectionAnalyzer(self.scorer, self.dataset)

        root.title("RRM Binding Score Analyzer")
        root.geometry("1200x800")

        try:
            root.state('zoomed')
        except tk.TclError: 
            try:
                root.attributes('-zoomed', True)
            except tk.TclError:
                print("Could not zoom window automatically.")

        self.main_frame = ttk.Frame(root, padding=10)
        self.main_frame.pack(fill=tk.BOTH, expand=True)

        self.main_frame.columnconfigure(0, weight=1)
        self.main_frame.rowconfigure(0, weight=0) 
        self.main_frame.rowconfigure(1, weight=3) 
        self.main_frame.rowconfigure(2, weight=1) 

        self._create_control_panel()
        self._create_plot_area()
        self._create_info_panel()

        if self.dataset.entries:
            self._update_info_display() 

    def _create_control_panel(self):
        self.control_frame = ttk.LabelFrame(self.main_frame, text="Controls", padding=10)
        self.control_frame.grid(row=0, column=0, sticky="ew", pady=5)
        
        col = 0
        ttk.Label(self.control_frame, text="Select RRM:").grid(row=0, column=col, sticky=tk.W, pady=5, padx=5)
        col+=1
        self.combined_ids_map = {} 
        self.display_ids = []

        for entry in self.dataset.entries:
            combined_id_val = entry['combined_id']
            rrm_id_val = entry['rrm_id']
            display_text = f"{rrm_id_val} ({entry['pdb_id']}_{entry['chain_id']})"
            self.display_ids.append(display_text)
            self.combined_ids_map[display_text] = (combined_id_val, rrm_id_val)

        self.display_ids.sort() 

        self.rrm_var = tk.StringVar()
        self.rrm_dropdown = ttk.Combobox(self.control_frame, textvariable=self.rrm_var, values=self.display_ids, width=35, state="readonly")
        self.rrm_dropdown.grid(row=0, column=col, sticky=tk.W, pady=5, padx=5)
        if self.display_ids:
            self.rrm_dropdown.current(0)
        self.rrm_dropdown.bind("<<ComboboxSelected>>", self._update_info_display)
        col+=1

        ttk.Label(self.control_frame, text="Window Size:").grid(row=0, column=col, sticky=tk.W, pady=5, padx=(10, 5))
        col+=1
        self.window_size_var = tk.IntVar(value=5)
        window_size_frame = ttk.Frame(self.control_frame)
        window_size_frame.grid(row=0, column=col, sticky=tk.W, pady=5)
        ttk.Radiobutton(window_size_frame, text="3", variable=self.window_size_var, value=3, command=self.clear_cache_on_param_change).pack(side=tk.LEFT)
        ttk.Radiobutton(window_size_frame, text="5", variable=self.window_size_var, value=5, command=self.clear_cache_on_param_change).pack(side=tk.LEFT)
        col+=1
        
        self.highlight_connections_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            self.control_frame, text="Highlight Connections", variable=self.highlight_connections_var
        ).grid(row=0, column=col, sticky=tk.W, pady=5, padx=10)
        col+=1

        ttk.Label(self.control_frame, text="Plot Type:").grid(row=0, column=col, sticky=tk.W, pady=5, padx=(10,5))
        col+=1
        self.plot_type_var = tk.StringVar(value="binding_score")
        plot_type_frame = ttk.Frame(self.control_frame)
        plot_type_frame.grid(row=0, column=col, sticky=tk.W, pady=5)
        ttk.Radiobutton(plot_type_frame, text="Binding Score",
                    variable=self.plot_type_var, value="binding_score").pack(side=tk.LEFT)
        ttk.Radiobutton(plot_type_frame, text="Connection Dist.",
                    variable=self.plot_type_var, value="connection_dist").pack(side=tk.LEFT)
        col+=1

        ttk.Label(self.control_frame, text="Iterations:").grid(row=0, column=col, sticky=tk.W, pady=5, padx=(10,0))
        col+=1
        self.iterations_var = tk.StringVar(value="1000") 
        iterations_entry = ttk.Entry(self.control_frame, textvariable=self.iterations_var, width=7)
        iterations_entry.grid(row=0, column=col, sticky=tk.W, pady=5, padx=(0,5))
        col+=1

        self.plot_button = ttk.Button(self.control_frame, text="Plot", command=self.plot_button_dispatcher)
        self.plot_button.grid(row=0, column=col, sticky=tk.W, pady=5, padx=10)
        col+=1

        self.status_var = tk.StringVar(value="Ready")
        self.status_label = ttk.Label(self.control_frame, textvariable=self.status_var, font=("Arial", 10, "italic"))
        self.status_label.grid(row=0, column=col, sticky="ew", pady=5, padx=10)
        self.control_frame.columnconfigure(col, weight=1) 

    def plot_button_dispatcher(self):
        selected_display_text = self.rrm_var.get()
        if not selected_display_text:
            messagebox.showerror("Selection Error", "Please select an RRM.")
            return

        self.clear_cache_for_current_rrm() 

        plot_type = self.plot_type_var.get()
        if plot_type == "binding_score":
            self.run_binding_score_analysis_and_plot()
        elif plot_type == "connection_dist":
            self.run_connection_distribution_analysis_and_plot()
        else:
            messagebox.showerror("Error", "Unknown plot type selected.")
            self.status_var.set("Ready")

    def _create_plot_area(self):
        self.plot_frame = ttk.LabelFrame(self.main_frame, text="Plot Area", padding=10) 
        self.plot_frame.grid(row=1, column=0, sticky="nsew", pady=5)

        self.figure = plt.Figure(figsize=(12, 6), dpi=100) # Slightly reduced height for better layout
        self.ax = self.figure.add_subplot(111)
        self.ax.text(0.5, 0.5, "Select an RRM and click 'Plot' to visualize",
                     ha='center', va='center', fontsize=12)
        self.ax.set_xticks([])
        self.ax.set_yticks([])

        self.canvas = FigureCanvasTkAgg(self.figure, self.plot_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=True)
        
        plot_notes_frame = ttk.Frame(self.plot_frame)
        plot_notes_frame.pack(fill=tk.X, pady=(5,0))

        self.plot_note_label = ttk.Label(plot_notes_frame,
                                         text="X-axis shows RNA sequence position indices (start of window).",
                                         font=("Arial", 8, "italic"))
        self.plot_note_label.pack(side=tk.LEFT, padx=5)

        self.connection_legend_label = ttk.Label(plot_notes_frame,
                                                 text="", 
                                                 font=("Arial", 8, "italic"),
                                                 foreground="red")
        self.connection_legend_label.pack(side=tk.RIGHT, padx=5)
        self.connection_legend_label.pack_forget() 

    def _create_info_panel(self):
        self.info_frame = ttk.LabelFrame(self.main_frame, text="Sequence & Analysis Information", padding=10)
        self.info_frame.grid(row=2, column=0, sticky="nsew", pady=5)

        self.info_frame.columnconfigure(0, weight=1)
        self.info_frame.rowconfigure(0, weight=1) 

        text_frame = ttk.Frame(self.info_frame)
        text_frame.grid(row=0, column=0, sticky="nsew")
        text_frame.columnconfigure(0, weight=1)
        text_frame.rowconfigure(0, weight=1)

        self.info_text = tk.Text(text_frame, wrap=tk.WORD, height=10, relief=tk.SUNKEN, borderwidth=1)
        v_scrollbar = ttk.Scrollbar(text_frame, orient="vertical", command=self.info_text.yview)
        self.info_text.configure(yscrollcommand=v_scrollbar.set)
        h_scrollbar = ttk.Scrollbar(text_frame, orient="horizontal", command=self.info_text.xview)
        self.info_text.configure(xscrollcommand=h_scrollbar.set)

        self.info_text.grid(row=0, column=0, sticky="nsew")
        v_scrollbar.grid(row=0, column=1, sticky="ns")
        h_scrollbar.grid(row=1, column=0, sticky="ew") 

    def _update_info_display(self, event=None):
        selected_display = self.rrm_var.get()
        if not selected_display or selected_display not in self.combined_ids_map:
            self.info_text.config(state=tk.NORMAL)
            self.info_text.delete(1.0, tk.END)
            self.info_text.insert(tk.END, "Please select an RRM.")
            self.info_text.config(state=tk.DISABLED)
            return

        combined_id_val, _ = self.combined_ids_map[selected_display]
        entry = next((e for e in self.dataset.entries if e['combined_id'] == combined_id_val), None)
        if not entry:
            self.info_text.config(state=tk.NORMAL)
            self.info_text.delete(1.0, tk.END)
            self.info_text.insert(tk.END, "Selected RRM entry not found in dataset.")
            self.info_text.config(state=tk.DISABLED)
            return

        self.info_text.config(state=tk.NORMAL)
        self.info_text.delete(1.0, tk.END)

        self.info_text.insert(tk.END, f"RRM ID: {entry['rrm_id']}\n")
        self.info_text.insert(tk.END, f"UniProt ID: {entry['uniprot_id']}\n")
        self.info_text.insert(tk.END, f"PDB ID: {entry['pdb_id']} (Chain {entry['chain_id']})\n")

        rna_seq = entry['rna_sequence']
        formatted_seq = "Dataset RNA Sequence:\n" # Clarify it's the dataset RNA
        chunk_size = 70 
        for i in range(0, len(rna_seq), chunk_size):
            formatted_seq += rna_seq[i:i+chunk_size] + "\n"
        self.info_text.insert(tk.END, formatted_seq)

        seq_len = len(rna_seq)
        window_size = self.window_size_var.get()
        num_windows = max(0, seq_len - window_size + 1)

        self.info_text.insert(tk.END, f"\nDataset RNA Length: {seq_len} nucleotides\n")
        self.info_text.insert(tk.END, f"Number of {window_size}-nt Windows in Dataset RNA: {num_windows}\n")

        if 'connections' in entry and entry['connections']:
            conn_str = ', '.join(map(str, entry['connections']))
            if len(conn_str) > 100: conn_str = conn_str[:100] + "..." 
            self.info_text.insert(tk.END, f"\nConnection RRM Indices (UniProt num.): {conn_str}\n")
        
        self.info_text.config(state=tk.DISABLED)

    def _display_top_scores_in_info(self, windows, scores, adjusted_scores, positions, mean_score=None, context_label=""):
        if not windows or not scores:
            return
        
        self.info_text.config(state=tk.NORMAL)
        scores_for_sorting = adjusted_scores if adjusted_scores else scores
        window_data = {}
        for i, window_seq in enumerate(windows):
            current_score = scores_for_sorting[i]
            current_pos = positions[i] if positions and i < len(positions) else -1
            
            if window_seq not in window_data or current_score > window_data[window_seq]['score']:
                 window_data[window_seq] = {'score': current_score, 'positions': [current_pos], 'original_score': scores[i]}
            elif window_seq in window_data and current_score == window_data[window_seq]['score']:
                 if current_pos not in window_data[window_seq]['positions']:
                      window_data[window_seq]['positions'].append(current_pos)

        sorted_windows = sorted(window_data.items(), key=lambda x: x[1]['score'], reverse=True)
        top_5 = sorted_windows[:5]

        self.info_text.insert(tk.END, f"\n\nTop 5 Scoring Windows ({context_label}):\n")
        for i, (window, data) in enumerate(top_5, 1):
            pos_str = "N/A"
            if data['positions'] and all(p >= 0 for p in data['positions']):
                pos_str = f"pos {', '.join(str(p) for p in sorted(data['positions']))}"
            
            score_display = f"{data['original_score']:.2f}"
            if data['original_score'] != data['score']: 
                score_display += f" (Adj: {data['score']:.2f})"
            self.info_text.insert(tk.END, f"{i}. {window}: {score_display} - {pos_str}\n")

        if mean_score is not None:
            self.info_text.insert(tk.END, f"\nMean Adjusted Score ({context_label}): {mean_score:.2f}\n")
        self.info_text.config(state=tk.DISABLED)

    def run_binding_score_analysis_and_plot(self):
        selected_display = self.rrm_var.get()
        combined_id_val, rrm_id_val = self.combined_ids_map[selected_display]
        entry = next((e for e in self.dataset.entries if e['combined_id'] == combined_id_val), None)
        
        if not entry:
            messagebox.showerror("Error", "Selected RRM entry not found.")
            self.status_var.set("Ready")
            return

        window_size = self.window_size_var.get()
        self.status_var.set(f"Processing binding scores for {rrm_id_val}...")
        self.root.update_idletasks()

        cache_key = f"{combined_id_val}_binding_{window_size}"
        if cache_key in self.cache:
            plot_data = self.cache[cache_key]
            self._display_binding_score_plot(plot_data) # Directly use cached data
            self._populate_binding_info_panel(plot_data)
            self.status_var.set("Ready")
            return

        threading.Thread(target=self._thread_run_binding_analysis, args=(entry, window_size, cache_key)).start()

    def _thread_run_binding_analysis(self, entry, window_size, cache_key):
        try:
            plot_data = self.binding_analyzer.prepare_binding_plot_data(
                entry, self.scorer, self.dataset, window_size, 
                getattr(self.scorer, 'protein_cds_dict', None) 
            )
            self.cache[cache_key] = plot_data
            self.root.after(0, self._finalize_binding_plot, plot_data) 
        except Exception as e:
            error_details = traceback.format_exc()
            self.root.after(0, messagebox.showerror, "Binding Analysis Error", f"An error occurred: {e}\n{error_details}")
            self.root.after(0, self.status_var.set, "Error")

    def _finalize_binding_plot(self, plot_data):
        if plot_data.get('error_message') and not plot_data.get('is_synthetic'): # Only show warning if not already synthetic
            messagebox.showwarning("Binding Score Info", plot_data['error_message'])
        
        self._display_binding_score_plot(plot_data)
        self._populate_binding_info_panel(plot_data)
        self.status_var.set("Ready")

    def _display_binding_score_plot(self, plot_data_dict): # Renamed to avoid conflict
        """Displays the RRM binding score plot using data from prepare_binding_plot_data."""
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        
        entry = plot_data_dict['entry']
        raw_windows = plot_data_dict['windows'] # Raw list of window sequences
        raw_scores = plot_data_dict['scores'] # Original scores from RRMScorer
        raw_adjusted_scores = plot_data_dict['adjusted_scores'] # Min-subtracted scores
        raw_positions = plot_data_dict['positions'] # Start positions of windows

        is_synthetic = plot_data_dict['is_synthetic']
        dataset_mean_score = plot_data_dict['dataset_mean_score']
        full_cds_mean_score = plot_data_dict['full_cds_exclusive_mean_score'] 
        
        # Boundary information from plot_data_dict
        # These are RNA nucleotide positions in the context of the plotted sequence
        dataset_start_rna_in_context = plot_data_dict['start_boundary'] 
        # 'end_boundary' from prepare_binding_plot_data was start of last window.
        # We need actual end for line/span.
        dataset_end_rna_in_context = None
        if dataset_start_rna_in_context is not None and plot_data_dict['using_protein_cds']:
            dataset_rna_len = len(entry['rna_sequence'])
            dataset_end_rna_in_context = dataset_start_rna_in_context + dataset_rna_len -1 # Inclusive end
        elif not plot_data_dict['using_protein_cds'] and raw_positions: # Dataset only context
            dataset_start_rna_in_context = min(raw_positions) if raw_positions else 0
            dataset_end_rna_in_context = max(raw_positions) + self.window_size_var.get() - 1 # End of last window plotted

        boundary_found = plot_data_dict['boundary_found']

        if not raw_windows or not raw_adjusted_scores:
            self.ax.text(0.5, 0.5, "No binding score data to display.", ha='center', va='center')
            self.canvas.draw()
            return

        # --- Data preparation for plotting (mimic original _display_plot) ---
        # We need a dense representation for x-axis if scores are sparse
        # but bars are plotted only where scores exist.
        
        # These are the positions where we actually have scores
        valid_plot_positions = []
        valid_plot_adjusted_scores = []

        if raw_positions and all(p >= 0 for p in raw_positions):
            # Filter out potential duplicate positions if scorer produced them, take first score
            # This step also ensures alignment between positions and scores for plotting
            temp_pos_score_map = {}
            for i, pos in enumerate(raw_positions):
                if pos not in temp_pos_score_map: # Keep first score for a position
                    temp_pos_score_map[pos] = raw_adjusted_scores[i]
            
            sorted_unique_positions = sorted(temp_pos_score_map.keys())
            valid_plot_positions = sorted_unique_positions
            valid_plot_adjusted_scores = [temp_pos_score_map[p] for p in sorted_unique_positions]
        else: # Fallback if positions are invalid (e.g. -1)
            valid_plot_positions = list(range(len(raw_windows)))
            valid_plot_adjusted_scores = raw_adjusted_scores
        
        if not valid_plot_positions: # Still no data to plot
             self.ax.text(0.5, 0.5, "No valid positions to plot.", ha='center', va='center')
             self.canvas.draw()
             return

        min_render_pos = min(valid_plot_positions)
        max_render_pos = max(valid_plot_positions)
        
        # Colormap setup
        min_adj_score_val = min(valid_plot_adjusted_scores) if valid_plot_adjusted_scores else 0
        max_adj_score_val = max(valid_plot_adjusted_scores) if valid_plot_adjusted_scores else 1
        score_range_for_norm = max_adj_score_val - min_adj_score_val
        
        norm = plt.Normalize(0, 1) # Normalize based on 0-1 for consistency after min-subtraction and range scaling
        if score_range_for_norm > 0:
            normalized_for_color = [(s - min_adj_score_val) / score_range_for_norm for s in valid_plot_adjusted_scores]
        else: # All scores are same
            normalized_for_color = [0.5] * len(valid_plot_adjusted_scores)

        cmap = plt.cm.get_cmap('YlGnBu') 
        bar_colors = [cmap(nc) for nc in normalized_for_color]

        bars = self.ax.bar(valid_plot_positions, valid_plot_adjusted_scores, color=bar_colors, width=1.0)

        # Moving average - calculated on the 'valid_plot_positions' and their scores
        if len(valid_plot_positions) > 3:
             ma_scores = self.binding_analyzer.calculate_moving_average(valid_plot_positions, valid_plot_adjusted_scores, window_size=63)
             smooth_line, = self.ax.plot(valid_plot_positions, ma_scores, color='red', linewidth=2.0, alpha=0.9, zorder=10) # Ensure MA on top

        # Boundary lines and shading (Corrected logic)
        boundary_line1_obj, boundary_line2_obj = None, None
        if boundary_found and dataset_start_rna_in_context is not None and dataset_end_rna_in_context is not None:
            boundary_line1_obj = self.ax.axvline(x=dataset_start_rna_in_context, color='blue', linestyle='--', linewidth=1.5, alpha=0.7)
            # For line, use inclusive end; for span, use exclusive end if needed by axvspan
            boundary_line2_obj = self.ax.axvline(x=dataset_end_rna_in_context, color='blue', linestyle='--', linewidth=1.5, alpha=0.7)
            # Span from start to length of dataset_rna
            if plot_data_dict['using_protein_cds']:
                self.ax.axvspan(dataset_start_rna_in_context, dataset_start_rna_in_context + len(entry['rna_sequence']), 
                                alpha=0.10, color='blue', zorder=0) # Span the full dataset length
        elif not plot_data_dict['using_protein_cds']: # Fallback boundaries for dataset-only view
            fb_start = min_render_pos
            fb_end = max_render_pos + self.window_size_var.get() -1 # End of last window
            boundary_line1_obj = self.ax.axvline(x=fb_start, color='dimgray', linestyle=':', linewidth=1.2)
            boundary_line2_obj = self.ax.axvline(x=fb_end, color='dimgray', linestyle=':', linewidth=1.2)
        
        # Highlight connections
        connections_rrm_indices = entry.get('connections', []) # These are UniProt RRM residue indices
        any_connections_highlighted = False
        if connections_rrm_indices and self.highlight_connections_var.get():
            window_len = self.window_size_var.get()
            for rrm_idx_uniprot in connections_rrm_indices:
                # Map RRM UniProt index to RNA nucleotide triplet in the *current plot context*
                # If plot is full CDS, rrm_idx_uniprot is used directly for mapping.
                # If plot is dataset_rna, and if connections are relative to PDB structure that aligns
                # to dataset_rna, then this mapping needs care. BUT entry['connections'] are already UniProt shifted.
                # So, these RRM indices are absolute w.r.t. the start of the protein's CDS.
                
                conn_rna_start_abs = rrm_idx_uniprot * 3
                conn_rna_end_abs = conn_rna_start_abs + 2

                # Now check overlap with plotted RNA windows (valid_plot_positions)
                for i, win_plot_start_pos in enumerate(valid_plot_positions):
                    win_plot_end_pos = win_plot_start_pos + window_len - 1
                    if max(win_plot_start_pos, conn_rna_start_abs) <= min(win_plot_end_pos, conn_rna_end_abs):
                        if i < len(bars):
                            bars[i].set_color('red') # Set entire bar to red
                            bars[i].set_zorder(5) # Ensure highlighted bars are prominent
                            any_connections_highlighted = True
        
        title_str = f'RRM Binding Scores for {entry["rrm_id"]} ({entry["pdb_id"]}_{entry["chain_id"]})'
        if is_synthetic: title_str += ' (SYNTHETIC DATA)'
        context_title = "Dataset RNA context"
        if plot_data_dict['using_protein_cds']: context_title = "Full CDS RNA context"
        title_str += f' ({context_title})'
            
        self.ax.set_xlabel('Position in RNA Sequence (Window Start Index)')
        self.ax.set_ylabel('RRM Binding Score (Adjusted)')
        self.ax.set_title(title_str, fontsize=11)

        # Legend
        legend_elements = []
        if 'smooth_line' in locals() and smooth_line: 
            legend_elements.append(Line2D([0], [0], color='red', lw=2, label='63-point Moving Average'))
        if dataset_mean_score is not None:
            legend_elements.append(Line2D([0], [0], color='green', linestyle='--', label=f'Dataset Mean: {dataset_mean_score:.2f}'))
        if full_cds_mean_score is not None:
             legend_elements.append(Line2D([0], [0], color='purple', linestyle='-.', label=f'Full CDS (excl. dataset): {full_cds_mean_score:.2f}'))
        
        boundary_label = "Dataset Boundaries"
        if not plot_data_dict['using_protein_cds']: boundary_label = "Plotted Range Boundaries"
        elif not boundary_found : boundary_label = "Dataset Boundaries (Estimated)"
            
        if boundary_line1_obj : # Add boundary to legend
            legend_elements.append(Line2D([0], [0], color='blue' if plot_data_dict['using_protein_cds'] and boundary_found else 'dimgray', 
                                          linestyle='--' if plot_data_dict['using_protein_cds'] and boundary_found else ':', label=boundary_label))

        if any_connections_highlighted:
            self.connection_legend_label.config(text="Red bars: Connection Position Overlap")
            self.connection_legend_label.pack(side=tk.RIGHT, padx=5)
            legend_elements.append(Line2D([0], [0], color='red', marker='s', linestyle='None', markersize=7, label='Connection Position'))
        else:
            self.connection_legend_label.pack_forget() 

        if legend_elements:
            self.ax.legend(handles=legend_elements, loc='best', fontsize=9)

        # X-axis ticks (restored original logic)
        if len(valid_plot_positions) > 0:
            # min_render_pos and max_render_pos already calculated
            plot_actual_range = max_render_pos - min_render_pos
            
            interval = 10 # Default
            if plot_actual_range > 1000: interval = 100
            elif plot_actual_range > 500: interval = 50
            elif plot_actual_range > 200: interval = 20
            elif plot_actual_range > 100: interval = 10
            elif plot_actual_range > 20: interval = 5
            elif plot_actual_range > 10: interval = 2
            else: interval = 1
            
            start_tick = (min_render_pos // interval) * interval
            # Ensure ticks cover the full range of plotted data
            tick_positions = list(range(start_tick, max_render_pos + interval, interval)) # Go up to or past max_render_pos
            # Filter ticks to be within a reasonable visual range if start_tick is far off
            tick_positions = [t for t in tick_positions if t >= min_render_pos - interval and t <= max_render_pos + interval]


            self.ax.set_xticks(tick_positions)
            self.ax.set_xticklabels([str(p) for p in tick_positions], rotation=45, ha='right', fontsize=8)
            self.ax.set_xlim(min_render_pos - 0.01*plot_actual_range, max_render_pos + 0.01*plot_actual_range) # Add small padding

            minor_interval = interval // 5 if interval >= 5 else 1
            if minor_interval > 0 and plot_actual_range > 10 : # Add minor ticks for larger ranges
                 minor_ticks = list(range(min_render_pos, max_render_pos + 1, minor_interval))
                 self.ax.set_xticks(minor_ticks, minor=True)

        self.ax.grid(True, axis='y', linestyle=':', alpha=0.6)
        self.ax.grid(True, axis='x', linestyle=':', alpha=0.3, which='major') # Grid for major x-ticks
        
        # Colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm) # norm should be 0-1 for normalized_for_color
        sm.set_array([]) # Necessary for scalar mappable
        cbar = self.figure.colorbar(sm, ax=self.ax, label='Normalized Binding Strength (0-1)', aspect=30, pad=0.03)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label('Normalized Binding Strength (0-1)', size=9)


        self.figure.tight_layout(rect=[0, 0, 0.97, 1]) # Adjust rect to make space for colorbar
        self.canvas.draw()

    def _populate_binding_info_panel(self, plot_data_dict):
        self.info_text.config(state=tk.NORMAL)
        self._update_info_display() 
        self.info_text.config(state=tk.NORMAL) 

        entry = plot_data_dict['entry']
        self.info_text.insert(tk.END, "\n--- Binding Score Analysis ---\n")
        if plot_data_dict['is_synthetic']:
            self.info_text.insert(tk.END, "NOTE: Displaying SYNTHETIC binding scores.\n")
        if plot_data_dict.get('error_message') and not plot_data_dict.get('is_synthetic'):
             self.info_text.insert(tk.END, f"Warning: {plot_data_dict['error_message']}\n")

        context = "Full CDS" if plot_data_dict['using_protein_cds'] else "Dataset RNA"
        self.info_text.insert(tk.END, f"Analysis Context: {context}\n")

        mean_to_display = plot_data_dict['full_cds_exclusive_mean_score'] if plot_data_dict['using_protein_cds'] else plot_data_dict['dataset_mean_score']
        self._display_top_scores_in_info(
            plot_data_dict['windows'], 
            plot_data_dict['scores'], 
            plot_data_dict['adjusted_scores'], 
            plot_data_dict['positions'],
            mean_to_display, # Display mean relevant to current context
            context_label=context
        )

        if plot_data_dict['dataset_mean_score'] is not None:
            self.info_text.insert(tk.END, f"Dataset RNA Mean Adjusted Score: {plot_data_dict['dataset_mean_score']:.2f}\n")
        if plot_data_dict['full_cds_exclusive_mean_score'] is not None:
            self.info_text.insert(tk.END, f"Full CDS Mean Adjusted Score: {plot_data_dict['full_cds_exclusive_mean_score']:.2f}\n")
        
        if plot_data_dict['dataset_mean_score'] is not None and plot_data_dict['full_cds_exclusive_mean_score'] is not None:
            diff = plot_data_dict['dataset_mean_score'] - plot_data_dict['full_cds_exclusive_mean_score']
            self.info_text.insert(tk.END, f"Difference (Dataset Mean - Full CDS Mean): {diff:.2f}\n")
            if diff > 0.1 : self.info_text.insert(tk.END, "  Interpretation: Dataset region shows higher avg. affinity than rest of CDS.\n")
            elif diff < -0.1: self.info_text.insert(tk.END, "  Interpretation: Dataset region shows lower avg. affinity than rest of CDS.\n")
            else: self.info_text.insert(tk.END, "  Interpretation: Dataset region shows similar avg. affinity to rest of CDS.\n")

        self.info_text.config(state=tk.DISABLED)

    def run_connection_distribution_analysis_and_plot(self):
        selected_display = self.rrm_var.get()
        combined_id_val, rrm_id_val = self.combined_ids_map[selected_display]
        entry = next((e for e in self.dataset.entries if e['combined_id'] == combined_id_val), None)

        if not entry:
            messagebox.showerror("Error", "Selected RRM entry not found.")
            self.status_var.set("Ready")
            return

        if not entry.get('connections'):
            messagebox.showwarning("No Connections", f"RRM {rrm_id_val} has no connection data.")
            self.status_var.set("Ready")
            self.figure.clear() # Clear plot if no data
            self.ax = self.figure.add_subplot(111)
            self.ax.text(0.5,0.5, f"No connection data for {rrm_id_val}.", ha='center', va='center')
            self.canvas.draw()
            self._populate_connection_info_panel({'entry': entry, 'error_message': "No connection data."})
            return
            
        try:
            iterations = int(self.iterations_var.get())
            if not (10 <= iterations <= 20000): 
                raise ValueError("Iterations must be between 10 and 20000.")
        except ValueError as e_iter:
            messagebox.showerror("Invalid Input", f"Iterations error: {e_iter}")
            self.status_var.set("Ready")
            return

        window_size = self.window_size_var.get()
        self.status_var.set(f"Analyzing connections for {rrm_id_val}...")
        self.root.update_idletasks()

        cache_key = f"{combined_id_val}_connection_{window_size}_{iterations}"
        if cache_key in self.cache:
            conn_plot_data = self.cache[cache_key]
            self._finalize_connection_plot(conn_plot_data) # Use finalize to show plot and info
            return # Status will be set by finalize

        threading.Thread(target=self._thread_run_connection_analysis, args=(entry, window_size, iterations, cache_key)).start()

    def _thread_run_connection_analysis(self, entry, window_size, iterations, cache_key):
        try:
            rrm_id = entry['rrm_id']
            original_connections_rrm = entry['connections'] 
            
            context_rna = entry['rna_sequence'].upper() 
            connections_in_context = original_connections_rrm 
            using_full_cds_context = False
            dataset_rna_for_context_check = entry['rna_sequence'].upper()


            base_uniprot_id = entry['uniprot_id'].split('_')[0]
            protein_cds_dict = getattr(self.scorer, 'protein_cds_dict', None)

            if protein_cds_dict and base_uniprot_id in protein_cds_dict:
                full_cds_rna_raw = protein_cds_dict[base_uniprot_id]
                if full_cds_rna_raw:
                    full_cds_rna = self.dataset._dna_to_rna(full_cds_rna_raw).upper()
                    start_pos_rna_in_full = self.conn_analyzer_instance.find_subsequence_start(full_cds_rna, dataset_rna_for_context_check)

                    if start_pos_rna_in_full != -1:
                        context_rna = full_cds_rna
                        rrm_shift = start_pos_rna_in_full // 3
                        connections_in_context = [c + rrm_shift for c in original_connections_rrm]
                        using_full_cds_context = True
            
            orig_score, random_scores, block_rna_window_coords = \
                self.conn_analyzer_instance.calculate_random_distribution(
                    rrm_id, context_rna, connections_in_context, iterations, window_size
            )
            
            conn_plot_data = {
                'entry': entry, 'original_score': orig_score, 'random_scores': random_scores,
                'iterations': iterations, 'block_rna_window_coords': block_rna_window_coords, 
                'context_rna_length': len(context_rna),
                'connections_in_context': connections_in_context, 
                'using_full_cds_context': using_full_cds_context, 'error_message': None
            }
            self.cache[cache_key] = conn_plot_data
            self.root.after(0, self._finalize_connection_plot, conn_plot_data)

        except Exception as e:
            error_details = traceback.format_exc()
            conn_plot_data = {'entry': entry, 'error_message': f"{e}\n{error_details}"} 
            self.root.after(0, self._finalize_connection_plot, conn_plot_data)

    def _finalize_connection_plot(self, conn_plot_data):
        entry_for_info = conn_plot_data.get('entry')
        if conn_plot_data.get('error_message'):
            messagebox.showerror("Connection Analysis Error", 
                                 f"An error occurred: {conn_plot_data['error_message']}")
            self.status_var.set("Error")
            if entry_for_info : self._populate_connection_info_panel(conn_plot_data) 
            return

        if not conn_plot_data.get('random_scores'): # Check if random_scores is empty or None
            messagebox.showwarning("Connection Analysis Info", "No random scores generated (original score might be valid but distribution cannot be plotted).")
            self.figure.clear()
            self.ax = self.figure.add_subplot(111)
            self.ax.text(0.5,0.5, "No random scores for distribution plot.", ha='center', va='center')
            if 'original_score' in conn_plot_data: # If we have original score, at least show that info
                 self.ax.set_title(f"Original Connection Score: {conn_plot_data['original_score']:.2f} (No random distribution)")
            self.canvas.draw()
        else:
            self._display_connection_distribution_plot(conn_plot_data)
        
        if entry_for_info : self._populate_connection_info_panel(conn_plot_data)
        self.status_var.set("Ready")

    def _display_connection_distribution_plot(self, conn_plot_data):
        dist_plotter = self.connection_analyzer_module.ConnectionDistributionPlot(
            self.figure, self.ax, self.canvas
        )
        dist_plotter.plot_distribution(
            conn_plot_data['original_score'],
            conn_plot_data['random_scores'],
            conn_plot_data['entry'], 
            conn_plot_data['iterations']
        )
        self.connection_legend_label.pack_forget() 

    def _populate_connection_info_panel(self, conn_plot_data):
        self.info_text.config(state=tk.NORMAL)
        self._update_info_display() 
        self.info_text.config(state=tk.NORMAL) 

        entry = conn_plot_data['entry']
        self.info_text.insert(tk.END, "\n--- Connection Analysis ---\n")

        if conn_plot_data.get('error_message') and not conn_plot_data.get('original_score'): # Only if analysis truly failed
            self.info_text.insert(tk.END, f"Error during analysis: {conn_plot_data['error_message']}\n")
            self.info_text.config(state=tk.DISABLED)
            return
            
        context_label = "Full CDS" if conn_plot_data.get('using_full_cds_context') else "Dataset RNA"
        self.info_text.insert(tk.END, f"Analysis Context: {context_label} (length {conn_plot_data.get('context_rna_length', 'N/A')} nt)\n")
        
        conns_in_ctx_str = str(conn_plot_data.get('connections_in_context', 'N/A'))
        self.info_text.insert(tk.END, f"RRM Connections Used (in context): {conns_in_ctx_str[:100]}{'...' if len(conns_in_ctx_str)>100 else ''}\n")

        original_score = conn_plot_data.get('original_score')
        if original_score is not None:
            self.info_text.insert(tk.END, f"Original Connection Score Sum: {original_score:.4f}\n")
        else:
            self.info_text.insert(tk.END, "Original Connection Score Sum: Not calculated.\n")


        random_scores = conn_plot_data.get('random_scores', [])
        if random_scores and original_score is not None : # Ensure original_score is available for comparison
            mean_random = np.mean(random_scores)
            std_random = np.std(random_scores)
            n_total = len(random_scores)
            percentile = sum(1 for x in random_scores if x <= original_score) / n_total * 100 if n_total > 0 else 0.0
            z_score = (original_score - mean_random) / std_random if std_random > 0 else 0.0
            
            p_val_ge, p_val_le = self.connection_analyzer_module.calculate_and_display_pvalue(original_score, random_scores)

            self.info_text.insert(tk.END, f"Random Shuffles Mean: {mean_random:.4f}, SD: {std_random:.4f}\n")
            self.info_text.insert(tk.END, f"Z-score: {z_score:.2f}\n")
            if p_val_ge is not None: self.info_text.insert(tk.END, f"P-value (Original >= Random): {p_val_ge:.4f}\n")
            if p_val_le is not None: self.info_text.insert(tk.END, f"P-value (Original <= Random): {p_val_le:.4f}\n")
            self.info_text.insert(tk.END, f"Percentile (Original score is >= {percentile:.1f}% of random scores)\n")
        elif original_score is not None: # Original score calculated, but no random scores
            self.info_text.insert(tk.END, "Random score distribution not generated or N/A.\n")


        rrm_blocks = self.conn_analyzer_instance.identify_connection_blocks(conn_plot_data.get('connections_in_context', []))
        self.info_text.insert(tk.END, f"\nIdentified RRM Connection Blocks ({len(rrm_blocks)} based on connections in context):\n")
        for i, rrm_block in enumerate(rrm_blocks[:5]): # Show first 5 blocks
             self.info_text.insert(tk.END, f"  Block {i+1} (RRM indices): {str(rrm_block)[:80]}{'...' if len(str(rrm_block))>80 else ''}\n")
        if len(rrm_blocks) > 5: self.info_text.insert(tk.END, "  ... and more blocks.\n")


        block_rna_windows = conn_plot_data.get('block_rna_window_coords', []) 
        if block_rna_windows:
            self.info_text.insert(tk.END, f"\nRNA Windows Scored for Original Connections ({len(block_rna_windows)} unique windows):\n")
            for i, (start_rna, end_rna) in enumerate(block_rna_windows[:5]): 
                self.info_text.insert(tk.END, f"  RNA Window {i+1}: Start {start_rna}, End {end_rna}\n")
            if len(block_rna_windows) > 5:
                self.info_text.insert(tk.END, "  ... and more RNA windows.\n")

        self.info_text.config(state=tk.DISABLED)

    def clear_cache_on_param_change(self, event=None):
        self.cache.clear() 

    def clear_cache_for_current_rrm(self):
        selected_display = self.rrm_var.get()
        if not selected_display or selected_display not in self.combined_ids_map:
            return

        combined_id_val, _ = self.combined_ids_map[selected_display]
        window_size = self.window_size_var.get()
        plot_type = self.plot_type_var.get()
        cache_key_pattern = None

        if plot_type == "binding_score":
            cache_key_pattern = f"{combined_id_val}_binding_{window_size}"
        elif plot_type == "connection_dist":
            iterations = self.iterations_var.get() 
            cache_key_pattern = f"{combined_id_val}_connection_{window_size}_{iterations}"
        
        if cache_key_pattern and cache_key_pattern in self.cache:
            del self.cache[cache_key_pattern]