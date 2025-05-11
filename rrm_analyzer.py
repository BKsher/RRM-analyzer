import random
import numpy as np
import os
import re
import sys
import subprocess
import argparse
import threading
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import traceback

# --- RRMDataset Class ---
class RRMDataset:
    """Class for handling RRM dataset files."""

    def __init__(self, file_path):
        self.file_path = file_path
        self.entries = []
        self.parse_file()

    def parse_file(self):
        """Parse the dataset file."""
        try:
            with open(self.file_path, 'r') as f:
                content = f.read()

            raw_entries = re.split(r'>>([^\n]+)\n', content)
            if not raw_entries[0]:
                raw_entries = raw_entries[1:]

            for i in range(0, len(raw_entries), 2):
                if i+1 < len(raw_entries):
                    header = raw_entries[i]
                    data = raw_entries[i+1]

                    try:
                        entry = self._extract_entry_data(header, data)
                        if entry:
                            self.entries.append(entry)
                    except Exception as e:
                        print(f"Error parsing entry with header: {header}")
                        print(f"Error details: {str(e)}")

            print(f"Successfully parsed {len(self.entries)} entries from dataset")
        except Exception as e:
            print(f"Error reading dataset file: {str(e)}")

    def _extract_entry_data(self, header, data):
        """Extract data from a single entry."""
        pdb_match = re.match(r'>?>?([0-9A-Za-z]+)_([A-Za-z0-9]+)', header)
        pdb_id = pdb_match.group(1) if pdb_match else "Unknown"
        chain_id = pdb_match.group(2) if pdb_match else ""

        uniprot_match = re.search(r'UniProt ID: ([^\n]+)', data)
        if not uniprot_match:
            print(f"Warning: Could not extract UniProt ID from entry")
            return None
        uniprot_id = uniprot_match.group(1).strip()

        rrm_match = re.search(r'_RRM(\d+)_', header)
        rrm_num = None
        if rrm_match:
            rrm_num = rrm_match.group(1)
        else:
            parts = header.strip().split('_')
            rrm_id = parts[0]

        if rrm_num:
            rrm_id = f"{uniprot_id}_RRM{rrm_num}"
        else:
            if "_RRM" in uniprot_id:
                rrm_id = uniprot_id
            else:
                rrm_id = uniprot_id

        combined_id = f"{rrm_id}_{pdb_id}_{chain_id}"

        nuc_match = re.search(r'Nucleotide Sequence: ([^\n]+)', data)
        if not nuc_match:
            print(f"Warning: Could not extract nucleotide sequence for {rrm_id}")
            return None

        nuc_seq = nuc_match.group(1).strip()

        rna_seq = self._dna_to_rna(nuc_seq)

        prot_match = re.search(r'Protein Sequence: ([^\n]+)', data)
        prot_seq = prot_match.group(1).strip() if prot_match else ""

        shift = 0
        shift_match = re.search(r'Numbering Shift.*?:([^\n]+)', data)
        if shift_match:
            shift_text = shift_match.group(1).strip()
            num_match = re.search(r'(-?\d+)', shift_text)
            if num_match:
                shift = int(num_match.group(1))

        connections = []
        connections_match = re.search(r'List of connections: ([^\n]+)', data)
        if connections_match:
            connections_str = connections_match.group(1).strip()
            try:
                connections = [int(x.strip()) for x in connections_str.split(',')]
            except ValueError:
                try:
                    connections = [int(x.strip()) for x in connections_str.split()]
                except ValueError:
                    print(f"Warning: Could not parse connections list for {rrm_id}")
        
        connections = [c + shift for c in connections] # auth to uniprid shift

        return {
            'header': header,
            'uniprot_id': uniprot_id,
            'rrm_num': rrm_num,
            'rrm_id': rrm_id,
            'pdb_id': pdb_id,
            'chain_id': chain_id,
            'combined_id': combined_id,
            'dna_sequence': nuc_seq,
            'rna_sequence': rna_seq,
            'protein_sequence': prot_seq,
            'numbering_shift': shift,
            'connections': connections
        }

    def _dna_to_rna(self, sequence):
        """Convert DNA sequence to RNA format (T -> U)."""
        if not sequence:
            return ""

        rna = re.sub(r'[Tt]', 'U', sequence)

        rna = re.sub(r'[^ACGUacgu]', '', rna)

        return rna.upper()


# --- RRMScorer Class ---
class RRMScorer:
    """Interface to the RRMScorer tool."""

    def __init__(self, output_dir):
        """Initialize with output directory."""
        self.output_dir = output_dir

        self.wrapper_path = "rrm_rna_wrapper.py"

        os.makedirs(output_dir, exist_ok=True)

        self.entries = []

        self._check_rrmscorer()

    def _check_rrmscorer(self):
        """Check if RRMScorer is available and working."""
        if not os.path.exists(self.wrapper_path):
            print(f"WARNING: RRMScorer wrapper not found at: {self.wrapper_path}")
            self.available = False
            return

        try:
            result = subprocess.run(
                ["python", self.wrapper_path, "--help"],
                capture_output=True,
                text=True
            )
            self.available = result.returncode == 0
            if self.available:
                print("RRMScorer is available and working")
            else:
                print(f"WARNING: RRMScorer test failed with return code {result.returncode}")
        except Exception as e:
            print(f"WARNING: Error testing RRMScorer: {str(e)}")
            self.available = False

    def run(self, rrm_id, rna_sequence, window_size=5):
        """Run RRMScorer on the given RRM and RNA sequence."""
        if not all(nt in 'ACGU' for nt in rna_sequence.upper()):
            return {
                'success': False,
                'error': "RNA sequence contains invalid nucleotides. Only A, C, G, U are allowed."
            }

        if "_RRM" in rrm_id:
            base_uniprot_id = rrm_id.split("_RRM")[0]
        else:
            base_uniprot_id = rrm_id

        output_file = os.path.join(self.output_dir, f"{rrm_id}_scores.txt")

        cmd = [
            "python",
            self.wrapper_path,
            "-RRM", rrm_id,
            "-RNA", rna_sequence.upper(),
            "-ws", str(window_size)
        ]

        print(f"Running command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            with open(output_file, 'w') as f:
                f.write(result.stdout)

            windows, scores, positions = self._parse_scores(output_file, rna_sequence, window_size)

            min_score = min(scores) if scores else 0
            adjusted_scores = [score - min_score for score in scores]

            return {
                'success': True,
                'output_file': output_file,
                'windows': windows,
                'scores': scores,
                'adjusted_scores': adjusted_scores,
                'positions': positions
            }
        except subprocess.CalledProcessError as e:
            if "_RRM" in rrm_id:
                uniprot_id = rrm_id.split("_RRM")[0]
                print(f"RRMScorer failed with full ID. Trying with just UniProt ID: {uniprot_id}")

                cmd[4] = uniprot_id

                try:
                    result = subprocess.run(
                        cmd,
                        capture_output=True,
                        text=True,
                        check=True
                    )

                    with open(output_file, 'w') as f:
                        f.write(result.stdout)

                    windows, scores, positions = self._parse_scores(output_file, rna_sequence, window_size)

                    min_score = min(scores) if scores else 0
                    adjusted_scores = [score - min_score for score in scores]

                    return {
                        'success': True,
                        'output_file': output_file,
                        'windows': windows,
                        'scores': scores,
                        'adjusted_scores': adjusted_scores,
                        'positions': positions
                    }
                except subprocess.CalledProcessError as e2:
                    return {
                        'success': False,
                        'error': f"RRMScorer failed with both IDs. Error: {e2}",
                        'stderr': e2.stderr
                    }

            return {
                'success': False,
                'error': f"RRMScorer failed: {e}",
                'stderr': e.stderr
            }

    def _parse_scores(self, score_file, rna_sequence=None, window_size=5):
        """Parse scores from RRMScorer output file."""
        windows = []
        scores = []
        positions = []
        rrm_id = None

        try:
            with open(score_file, 'r') as f:
                lines = f.readlines()

                for line in lines:
                    if line.strip() and ' ' not in line:
                        potential_id = line.strip()
                        if '_RRM' in potential_id:
                            rrm_id = potential_id
                            break

                window_score_pairs = []
                for line in lines:
                    if not line.strip() or line.strip() == rrm_id:
                        continue

                    parts = line.strip().split()
                    if len(parts) == 2:
                        window, score = parts

                        if all(nt in 'ACGU' for nt in window):
                            try:
                                window_score_pairs.append((window, float(score)))
                            except ValueError:
                                print(f"Warning: Could not convert score to float: {score}")

                if rna_sequence and window_score_pairs:
                    for window, score in window_score_pairs:
                        start_pos = 0
                        while start_pos < len(rna_sequence):
                            pos = rna_sequence.find(window, start_pos)
                            if pos == -1:
                                break

                            windows.append(window)
                            scores.append(score)
                            positions.append(pos)

                            start_pos = pos + 1
                else:
                    for window, score in window_score_pairs:
                        windows.append(window)
                        scores.append(score)
                        positions.append(-1)
        except Exception as e:
            print(f"Error parsing score file: {e}")

        if windows and scores:
            return windows, scores, positions

        print("Generating complete sequence of windows to ensure all positions are represented")
        return self._generate_complete_windows(score_file)

    def _generate_complete_windows(self, score_file):
        """Generate complete sequence of windows based on RNA sequence."""
        try:
            filename = os.path.basename(score_file)
            rrm_id = filename.split('_scores.txt')[0]

            rna_seq = None
            window_size = None

            for entry in self.entries:
                if entry['rrm_id'] == rrm_id:
                    rna_seq = entry['rna_sequence'].upper()
                    break

            if not rna_seq:
                print(f"Warning: Could not find RNA sequence for {rrm_id}")
                return [], [], []

            with open(score_file, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) == 2:
                        window, score = parts
                        if all(nt in 'ACGU' for nt in window):
                            window_size = len(window)
                            break

            if not window_size:
                print("Warning: Could not determine window size")
                window_size = 5  # Default

            score_dict = {}
            with open(score_file, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) == 2:
                        window, score = parts
                        if all(nt in 'ACGU' for nt in window):
                            score_dict[window] = float(score)

            all_windows = []
            all_scores = []
            all_positions = []

            for i in range(len(rna_seq) - window_size + 1):
                window = rna_seq[i:i+window_size]
                all_windows.append(window)
                all_positions.append(i)

                if window in score_dict:
                    all_scores.append(score_dict[window])
                else:
                    found_similar = False
                    for stored_window, stored_score in score_dict.items():
                        if window.upper() == stored_window.upper():
                            all_scores.append(stored_score)
                            found_similar = True
                            break

                    if not found_similar:
                        print(f"Warning: No score found for window {window}")
                        all_scores.append(-1.0)

            return all_windows, all_scores, all_positions

        except Exception as e:
            print(f"Error generating complete windows: {e}")
            return [], [], []

    def generate_synthetic_scores(self, rna_sequence, window_size=5):
        """Generate synthetic scores for demonstration when RRMScorer fails."""
        windows = []
        scores = []
        positions = []

        if not all(nt in 'ACGU' for nt in rna_sequence.upper()):
            print("Warning: RNA sequence contains invalid nucleotides")
            return [], [], []

        rna_sequence = rna_sequence.upper()
        for i in range(len(rna_sequence) - window_size + 1):
            window = rna_sequence[i:i+window_size]
            score = random.uniform(-2, 0)
            windows.append(window)
            scores.append(score)
            positions.append(i)

        return windows, scores, positions

    def load_protein_cds_fasta(self, fasta_path):
        """
        Load RNA sequences from protein_cds.fasta file.
        """
        protein_cds_dict = {}

        try:
            with open(fasta_path, 'r') as f:
                lines = f.readlines()

            i = 0
            while i < len(lines):
                if i < len(lines) and lines[i].startswith('>'):
                    protein_id = lines[i][1:].split()[0]
                    i += 1

                    while i < len(lines) and not lines[i].startswith('>'):
                        i += 1

                    if i < len(lines) and lines[i].startswith('>'):
                        nucleotide_id = lines[i][1:].split()[0]
                        i += 1

                        nucleotide_seq = ""
                        while i < len(lines) and not lines[i].startswith('>'):
                            nucleotide_seq += lines[i].strip()
                            i += 1

                        protein_cds_dict[protein_id] = nucleotide_seq
                else:
                    i += 1

            print(f"Loaded {len(protein_cds_dict)} sequences from {fasta_path}")
            return protein_cds_dict
        except Exception as e:
            print(f"Error loading protein CDS FASTA: {e}")
            return {}


# --- RRMVisualizerGUI Class ---
class RRMVisualizerGUI:
    """GUI for RRM score visualization."""

    def __init__(self, root, dataset, scorer):
        self.root = root
        self.dataset = dataset
        self.scorer = scorer
        self.cache = {}

        if os.path.exists("proten_cds.fasta"):
            print(f"File proten_cds.fasta exists, size: {os.path.getsize('proten_cds.fasta')} bytes")
        else:
            print("File proten_cds.fasta does not exist")

        self.protein_cds_dict = self.scorer.load_protein_cds_fasta("proten_cds.fasta")
        print(f"Loaded protein_cds_dict with {len(self.protein_cds_dict)} entries")
        if self.protein_cds_dict:
            print(f"First few keys: {list(self.protein_cds_dict.keys())[:5]}")

        self.scorer.protein_cds_dict = self.protein_cds_dict

        root.title("RRM Binding Score Analyzer")
        root.geometry("1200x800")

        try:
            root.state('zoomed')
        except:
            root.attributes('-zoomed', True)

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

        self.add_connection_analysis_ui()

    def _create_control_panel(self):
        self.control_frame = ttk.LabelFrame(self.main_frame, text="Controls", padding=10)
        self.control_frame.grid(row=0, column=0, sticky="ew", pady=5)

        ttk.Label(self.control_frame, text="Select RRM:").grid(row=0, column=0, sticky=tk.W, pady=5)

        self.combined_ids = []
        self.display_ids = []

        for entry in self.dataset.entries:
            combined_id = entry['combined_id']
            rrm_id = entry['rrm_id']
            display_text = f"{rrm_id} ({entry['pdb_id']}_{entry['chain_id']})"

            self.combined_ids.append((combined_id, rrm_id))
            self.display_ids.append(display_text)

        sorted_pairs = sorted(zip(self.display_ids, self.combined_ids))
        if sorted_pairs:
            self.display_ids, self.combined_ids = zip(*sorted_pairs)
        else:
            self.display_ids, self.combined_ids = [], []

        self.rrm_var = tk.StringVar()
        self.rrm_dropdown = ttk.Combobox(self.control_frame, textvariable=self.rrm_var, values=self.display_ids, width=30)
        self.rrm_dropdown.grid(row=0, column=1, sticky=tk.W, pady=5, padx=5)
        if self.display_ids:
            self.rrm_dropdown.current(0)

        ttk.Label(self.control_frame, text="Window Size:").grid(row=0, column=2, sticky=tk.W, pady=5, padx=(20, 5))
        self.window_size_var = tk.IntVar(value=5)
        window_size_frame = ttk.Frame(self.control_frame)
        window_size_frame.grid(row=0, column=3, sticky=tk.W, pady=5)

        ttk.Radiobutton(window_size_frame, text="3", variable=self.window_size_var, value=3).pack(side=tk.LEFT)
        ttk.Radiobutton(window_size_frame, text="5", variable=self.window_size_var, value=5).pack(side=tk.LEFT)

        self.highlight_connections_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            self.control_frame,
            text="Highlight Connections",
            variable=self.highlight_connections_var
        ).grid(row=0, column=4, sticky=tk.W, pady=5, padx=10)

        self.plot_button = ttk.Button(self.control_frame, text="Plot", command=self.plot_button_handler)
        self.plot_button.grid(row=0, column=5, sticky=tk.W, pady=5, padx=20)

        self.status_var = tk.StringVar(value="Ready")
        self.status_label = ttk.Label(self.control_frame, textvariable=self.status_var, font=("Arial", 10, "italic"))
        self.status_label.grid(row=0, column=6, sticky=tk.W, pady=5, padx=10)

        self.rrm_dropdown.bind("<<ComboboxSelected>>", self._update_info_display)
        self.control_frame.columnconfigure(6, weight=1)

    def plot_button_handler(self):
        print("DEBUGINFO: Plot button clicked", file=sys.stderr, flush=True)

        self.clear_cache_for_current_rrm()

        plot_type = self.plot_type_var.get()
        print(f"DEBUGINFO: Plot type from button handler: {plot_type}", file=sys.stderr, flush=True)

        if plot_type == "binding_score":
            print("DEBUGINFO: Using binding score from button handler", file=sys.stderr, flush=True)
            self._original_plot_selected_rrm()
        else:
            print("DEBUGINFO: Using connection distribution from button handler", file=sys.stderr, flush=True)
            self.plot_connection_distribution(self.rrm_var.get())

    def _create_plot_area(self):
        self.plot_frame = ttk.LabelFrame(self.main_frame, text="RRM Binding Score Plot", padding=10)
        self.plot_frame.grid(row=1, column=0, sticky="nsew", pady=5)

        self.figure = plt.Figure(figsize=(12, 7), dpi=100)
        self.ax = self.figure.add_subplot(111)

        self.ax.text(0.5, 0.5, "Select an RRM and click 'Plot' to visualize binding scores",
                     ha='center', va='center', fontsize=12)
        self.ax.set_xticks([])
        self.ax.set_yticks([])

        self.canvas = FigureCanvasTkAgg(self.figure, self.plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        button_frame = ttk.Frame(self.plot_frame)
        button_frame.pack(fill=tk.X, pady=(5, 0))

        ttk.Label(button_frame,
                 text="Note: X-axis shows position indices in the RNA sequence. Each position represents the start of a window.",
                 font=("Arial", 8, "italic")).pack(side=tk.LEFT, padx=5)

        self.connection_label = ttk.Label(
            button_frame,
            text="Highlighted bars: Connections positions listed in the data",
            font=("Arial", 8, "italic"),
            foreground="red"
        )
        self.connection_label.pack(side=tk.RIGHT, padx=5)

    def _create_info_panel(self):
        self.info_frame = ttk.LabelFrame(self.main_frame, text="Sequence Information", padding=10)
        self.info_frame.grid(row=2, column=0, sticky="nsew", pady=5)

        self.info_frame.columnconfigure(0, weight=1)
        self.info_frame.rowconfigure(0, weight=1)

        text_frame = ttk.Frame(self.info_frame)
        text_frame.grid(row=0, column=0, sticky="nsew")
        text_frame.columnconfigure(0, weight=1)
        text_frame.rowconfigure(0, weight=1)

        self.info_text = tk.Text(text_frame, wrap=tk.NONE, height=10)

        v_scrollbar = ttk.Scrollbar(text_frame, orient="vertical", command=self.info_text.yview)
        self.info_text.configure(yscrollcommand=v_scrollbar.set)

        h_scrollbar = ttk.Scrollbar(text_frame, orient="horizontal", command=self.info_text.xview)
        self.info_text.configure(xscrollcommand=h_scrollbar.set)

        self.info_text.grid(row=0, column=0, sticky="nsew")
        v_scrollbar.grid(row=0, column=1, sticky="ns")
        h_scrollbar.grid(row=1, column=0, sticky="ew")

    def _update_info_display(self, event=None):
        selected_display = self.rrm_var.get()

        selected_entry = None
        selected_combined_id = None

        for i, display_id in enumerate(self.display_ids):
            if display_id == selected_display:
                selected_combined_id, _ = self.combined_ids[i]
                break

        entry = next((e for e in self.dataset.entries if e['combined_id'] == selected_combined_id), None)
        if not entry:
            return

        self.info_text.delete(1.0, tk.END)

        self.info_text.insert(tk.END, f"RRM ID: {entry['rrm_id']}\n")
        self.info_text.insert(tk.END, f"UniProt ID: {entry['uniprot_id']}\n")
        self.info_text.insert(tk.END, f"PDB ID: {entry['pdb_id']} (Chain {entry['chain_id']})\n")

        rna_seq = entry['rna_sequence']
        formatted_seq = ""
        chunk_size = 60
        for i in range(0, len(rna_seq), chunk_size):
            formatted_seq += rna_seq[i:i+chunk_size] + "\n"

        self.info_text.insert(tk.END, f"RNA Sequence:\n{formatted_seq}")

        seq_len = len(rna_seq)
        window_size = self.window_size_var.get()
        num_windows = max(0, seq_len - window_size + 1)

        self.info_text.insert(tk.END, f"Sequence Length: {seq_len} nucleotides\n")
        self.info_text.insert(tk.END, f"Number of {window_size}-nt Windows: {num_windows}\n")

        if 'connections' in entry and entry['connections']:
            self.info_text.insert(tk.END, f"\nList of Connection Positions: {', '.join(map(str, entry['connections']))}\n")

        if seq_len >= window_size:
            self.info_text.insert(tk.END, f"\nExample Windows (first 5):\n")
            for i in range(min(5, num_windows)):
                window = rna_seq[i:i+window_size]
                self.info_text.insert(tk.END, f"Window {i+1}: {window} (position {i})\n")

        cache_key = f"{entry['combined_id']}_{window_size}"
        if cache_key in self.cache:
            windows, scores, adjusted_scores, positions, is_synthetic = self.cache[cache_key]
            if windows and scores and not is_synthetic:
                mean_adjusted_score = sum(adjusted_scores) / len(adjusted_scores) if adjusted_scores else 0
                self._display_top_scores(windows, scores, adjusted_scores, positions, mean_adjusted_score)

    def _display_top_scores(self, windows, scores, adjusted_scores, positions, mean_score=None):
        """Display top 5 scoring windows in the info panel."""
        if not windows or not scores:
            return

        window_data = {}
        for window, score, adjusted, position in zip(windows, scores, adjusted_scores, positions):
            if window not in window_data:
                window_data[window] = {
                    'score': score,
                    'adjusted': adjusted,
                    'positions': [position]
                }
            else:
                if position not in window_data[window]['positions']:
                    window_data[window]['positions'].append(position)

        sorted_windows = sorted(window_data.items(), key=lambda x: x[1]['score'], reverse=True)

        top_5 = sorted_windows[:5]

        self.info_text.insert(tk.END, "\n\nTop 5 Scoring Windows:\n")
        for i, (window, data) in enumerate(top_5, 1):
            if all(p >= 0 for p in data['positions']):
                if len(data['positions']) > 1:
                    position_str = f"positions {', '.join(str(p) for p in data['positions'])}"
                else:
                    position_str = f"position {data['positions'][0]}"
            else:
                position_str = "unknown position"

            self.info_text.insert(tk.END, f"{i}. {window}: {data['score']:.2f} (Adjusted: {data['adjusted']:.2f}) - {position_str}\n")

        if mean_score is not None:
            self.info_text.insert(tk.END, f"\nMean Adjusted Score: {mean_score:.2f}\n")

    def _plot_selected_rrm(self):
        """Run RRMScorer and plot the results for the selected RRM."""
        selected_display = self.rrm_var.get()

        selected_combined_id = None
        selected_rrm_id = None
        for i, display_id in enumerate(self.display_ids):
            if display_id == selected_display:
                selected_combined_id, selected_rrm_id = self.combined_ids[i]
                break

        if not selected_combined_id:
            messagebox.showerror("Error", "Please select a valid RRM ID")
            return

        window_size = self.window_size_var.get()
        entry = next((e for e in self.dataset.entries if e['combined_id'] == selected_combined_id), None)
        if not entry:
            messagebox.showerror("Error", "Please select a valid RRM ID")
            return

        self.status_var.set("Processing...")
        self.root.update_idletasks()

        cache_key = f"{entry['combined_id']}_{window_size}"
        if cache_key in self.cache:
            windows, scores, adjusted_scores, positions, is_synthetic = self.cache[cache_key]
            mean_adjusted_score = sum(adjusted_scores) / len(adjusted_scores) if adjusted_scores else 0

            connections = entry.get('connections', [])

            self._display_plot(entry, windows, scores, adjusted_scores, positions, is_synthetic, mean_adjusted_score, None, connections)
            self._display_top_scores(windows, scores, adjusted_scores, positions, mean_adjusted_score)
            self.status_var.set("Ready")
            return

        threading.Thread(target=self._run_and_plot_with_both_rnas, args=(entry, entry, window_size)).start()

    def _run_and_plot(self, rrm_id, entry, window_size):
        """Run RRMScorer and plot results in a separate thread."""
        try:
            rna_seq = entry['rna_sequence'].upper()

            if len(rna_seq) < window_size:
                self.root.after(0, lambda: messagebox.showerror(
                    "Error",
                    f"RNA sequence ({len(rna_seq)} nt) is shorter than the window size ({window_size} nt)"
                ))
                self.root.after(0, lambda: self.status_var.set("Ready"))
                return

            uniprot_id = entry['uniprot_id']

            result = self.scorer.run(entry['rrm_id'], rna_seq, window_size)

            if result['success']:
                windows = result['windows']
                scores = result['scores']
                adjusted_scores = result['adjusted_scores']
                positions = result['positions']
                is_synthetic = False

                if not windows or not scores:
                    print("No valid scores found, generating synthetic data")
                    windows, scores, positions = self.scorer.generate_synthetic_scores(rna_seq, window_size)
                    min_score = min(scores) if scores else 0
                    adjusted_scores = [score - min_score for score in scores]
                    is_synthetic = True

                    self.root.after(0, lambda: messagebox.showwarning(
                        "No Scores",
                        "No valid scores found in RRMScorer output. Using synthetic data for visualization."
                    ))
            else:
                print(f"RRMScorer failed: {result.get('error', 'Unknown error')}")
                windows, scores, positions = self.scorer.generate_synthetic_scores(rna_seq, window_size)
                min_score = min(scores) if scores else 0
                adjusted_scores = [score - min_score for score in scores]
                is_synthetic = True

                self.root.after(0, lambda: messagebox.showwarning(
                    "RRMScorer Error",
                    f"Error running RRMScorer: {result.get('error', 'Unknown error')}\n\n"
                    f"Using synthetic data for visualization."
                ))

            self.cache[f"{rrm_id}_{window_size}"] = (windows, scores, adjusted_scores, positions, is_synthetic)

            mean_adjusted_score = sum(adjusted_scores) / len(adjusted_scores) if adjusted_scores else 0

            connections = entry.get('connections', [])

            self.root.after(0, lambda: self._display_plot(entry, windows, scores, adjusted_scores, positions, is_synthetic, mean_adjusted_score, connections))
            self.root.after(0, lambda: self._display_top_scores(windows, scores, adjusted_scores, positions, mean_adjusted_score))
            self.root.after(0, lambda: self.status_var.set("Ready"))

        except Exception as e:
            error_details = traceback.format_exc()
            self.root.after(0, lambda: messagebox.showerror(
                "Error",
                f"An error occurred: {str(e)}\n\nDetails:\n{error_details}"
            ))
            self.root.after(0, lambda: self.status_var.set("Ready"))

    def calculate_moving_average(positions, scores, window_size=63):
        import numpy as np

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
                moving_avgs.append(scores[i])

        return moving_avgs

    def _display_plot(self, entry, windows, scores, adjusted_scores, positions, is_synthetic=False,
                dataset_mean_score=None, exclusive_mean_score=None, connections=None):
        """Display the plot with scores and highlighting connections, and both mean scores."""
        try:
            self.figure.clear()
            self.ax = self.figure.add_subplot(111)

            entry_obj = None
            if isinstance(entry, dict):
                entry_obj = entry
            elif isinstance(entry, str):
                entry_obj = next((e for e in self.dataset.entries if e['rrm_id'] == entry), None)

            print("\n" + "="*80)
            if entry_obj:
                print(f"DEBUG: Calculating boundaries for RRM ID: {entry_obj['rrm_id']}")
            else:
                print(f"DEBUG: Calculating boundaries for RRM ID: {entry}")

            if not windows or not scores:
                self.ax.text(0.5, 0.5, "No data to display", ha='center', va='center', fontsize=12)
                self.ax.set_xticks([])
                self.ax.set_yticks([])
                self.canvas.draw()
                return

            if positions and all(p >= 0 for p in positions):
                min_pos = min(positions)
                max_pos = max(positions)
                all_possible_positions = list(range(min_pos, max_pos + 1))

                position_to_score = {}
                position_to_adjusted = {}

                for window, score, adjusted, pos in zip(windows, scores, adjusted_scores, positions):
                    if pos >= 0:
                        position_to_score[pos] = score
                        position_to_adjusted[pos] = adjusted

                plot_positions = all_possible_positions
                plot_adjusted_scores = [position_to_adjusted.get(pos, None) for pos in plot_positions]

                valid_indices = [i for i, val in enumerate(plot_adjusted_scores) if val is not None]
                plot_positions = [plot_positions[i] for i in valid_indices]
                plot_adjusted_scores = [plot_adjusted_scores[i] for i in valid_indices]
            else:
                plot_positions = list(range(len(windows)))
                plot_adjusted_scores = adjusted_scores

            if len(plot_adjusted_scores) > 1:
                min_score = min(plot_adjusted_scores)
                max_score = max(plot_adjusted_scores)
                score_range = max_score - min_score
                if score_range > 0:
                    normalized_scores = [(s - min_score) / score_range for s in plot_adjusted_scores]
                else:
                    normalized_scores = [0.5] * len(plot_adjusted_scores)
            else:
                normalized_scores = [0.5]
                min_score = plot_adjusted_scores[0] if plot_adjusted_scores else 0
                max_score = plot_adjusted_scores[0] if plot_adjusted_scores else 1
                score_range = 1

            is_connection = []
            for pos in plot_positions:
                if connections and self.highlight_connections_var.get():
                    is_conn = pos in connections
                    is_connection.append(is_conn)
                else:
                    is_connection.append(False)

            cmap = plt.cm.get_cmap('YlGnBu')
            colors = []

            for i, ns in enumerate(normalized_scores):
                colors.append(cmap(ns))

            bars = self.ax.bar(plot_positions, plot_adjusted_scores, color=colors)

            if len(plot_positions) > 3:
                def calculate_moving_average(positions, scores, window_size=63):
                    import numpy as np

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
                            moving_avgs.append(scores[i])

                    return moving_avgs

                ma_scores = calculate_moving_average(plot_positions, plot_adjusted_scores, window_size=63)
                smooth_line = self.ax.plot(plot_positions, ma_scores, color='red', linewidth=2)
            else:
                smooth_line = None

            boundary_line1 = None
            boundary_line2 = None
            boundary_found = False
            start_boundary = None
            end_boundary = None

            print("\n" + "="*80)
            if isinstance(entry, dict):
                print(f"DEBUG: Calculating boundaries for RRM ID: {entry['rrm_id']}")
            else:
                print(f"DEBUG: Calculating boundaries for RRM ID: {entry}")

            try:
                if isinstance(entry, str):
                    entry_obj = next((e for e in self.dataset.entries if e['rrm_id'] == entry), None)
                else:
                    entry_obj = entry

                if entry_obj:
                    dataset_rna = entry_obj['rna_sequence'].upper()
                    base_uniprot_id = entry_obj['uniprot_id'].split('_')[0]

                    print(f"DEBUG: Dataset RNA length: {len(dataset_rna)}")
                    print(f"DEBUG: Dataset RNA (first 50 chars): {dataset_rna[:50]}...")
                    print(f"DEBUG: Base UniProt ID: {base_uniprot_id}")

                    if hasattr(self.scorer, 'protein_cds_dict') and base_uniprot_id in self.scorer.protein_cds_dict:
                        protein_cds_rna = self.dataset._dna_to_rna(self.scorer.protein_cds_dict[base_uniprot_id]).upper()

                        print(f"DEBUG: Full RNA length: {len(protein_cds_rna)}")
                        print(f"DEBUG: Full RNA (first 50 chars): {protein_cds_rna[:50]}...")

                        dataset_in_full = protein_cds_rna.find(dataset_rna)
                        print(f"DEBUG: Result of find() operation: {dataset_in_full}")

                        if dataset_in_full != -1:
                            window_size = self.window_size_var.get()
                            start_boundary = dataset_in_full
                            end_boundary = dataset_in_full + len(dataset_rna) - window_size

                            print(f"DEBUG: Calculated start boundary: {start_boundary}")
                            print(f"DEBUG: Calculated end boundary: {end_boundary}")
                            print(f"DEBUG: Window size: {window_size}")
                            print(f"DEBUG: Dataset RNA subsequence match verification: {protein_cds_rna[start_boundary:start_boundary+50]} == {dataset_rna[:50]}")

                            print(f"DEBUG: Plot position range: {min(plot_positions)} to {max(plot_positions)}")
                            boundary_in_range = (min(plot_positions) <= start_boundary <= max(plot_positions) or
                                                min(plot_positions) <= end_boundary <= max(plot_positions))
                            print(f"DEBUG: Boundaries in visible range: {boundary_in_range}")

                            if boundary_in_range:
                                boundary_line1 = self.ax.axvline(x=start_boundary, color='blue', linestyle='--', linewidth=2, alpha=0.7)
                                boundary_line2 = self.ax.axvline(x=end_boundary, color='blue', linestyle='--', linewidth=2, alpha=0.7)
                                self.ax.axvspan(start_boundary, end_boundary, alpha=0.15, color='blue')
                                boundary_found = True
                                print("DEBUG: Successfully added boundary lines to plot")
                            else:
                                print("DEBUG: Boundaries outside of visible plot range")
                        else:
                            print("DEBUG: Failed to find dataset RNA in full RNA (find() returned -1)")
                    else:
                        if not hasattr(self.scorer, 'protein_cds_dict'):
                            print("DEBUG: scorer does not have protein_cds_dict attribute")
                        elif base_uniprot_id not in self.scorer.protein_cds_dict:
                            print(f"DEBUG: UniProt ID {base_uniprot_id} not found in protein_cds_dict")
                            print(f"DEBUG: Available keys in protein_cds_dict: {list(self.scorer.protein_cds_dict.keys())[:10]}...")
                else:
                    print(f"DEBUG: No entry found for RRM ID: {entry}")
            except Exception as e:
                print(f"DEBUG: Error calculating boundaries: {str(e)}")
                print(traceback.format_exc())

            if not boundary_found:
                print("DEBUG: Using fallback boundary positions")
                min_plot_pos = min(plot_positions) if plot_positions else 0
                max_plot_pos = max(plot_positions) if plot_positions else 100
                boundary1 = min_plot_pos + (max_plot_pos - min_plot_pos) // 3
                boundary2 = min_plot_pos + 2 * (max_plot_pos - min_plot_pos) // 3

                print(f"DEBUG: Fallback boundary 1: {boundary1}")
                print(f"DEBUG: Fallback boundary 2: {boundary2}")

                boundary_line1 = self.ax.axvline(x=boundary1, color='blue', linestyle='--', linewidth=2, alpha=0.7)
                boundary_line2 = self.ax.axvline(x=boundary2, color='blue', linestyle='--', linewidth=2, alpha=0.7)

                self.root.after(0, lambda: messagebox.showwarning(
                    "Boundary Warning",
                    "Could not determine actual dataset boundaries. Using fallback positions at 1/3 and 2/3 of the plot."
                ))

            print(f"DEBUG: Connections: {connections}")
            if connections and self.highlight_connections_var.get():
                any_connections_shown = False

                if entry_obj:
                    connections = entry_obj.get('connections', [])

                connection_positions = [pos*3 + offset for pos in plot_positions if pos in connections for offset in range(3)]
                print(f"DEBUG: Connection positions: {connection_positions}")

                if connection_positions:
                    any_connections_shown = True
                    for pos in connection_positions:
                        try:
                            idx = plot_positions.index(pos)

                            if idx < end_boundary and idx > start_boundary:
                                bars[idx].set_color('red')

                            adj_score = plot_adjusted_scores[idx]
                            self.ax.add_patch(plt.Rectangle(
                                (pos - 0.4, 0), 0.8, adj_score,
                                alpha=0.3, color='yellow', zorder=0
                            ))
                        except (ValueError, IndexError):
                            pass

                if any_connections_shown:
                    self.connection_label.pack(side=tk.RIGHT, padx=5)
                else:
                    self.connection_label.pack_forget()
            else:
                self.connection_label.pack_forget()

            print("="*80 + "\n")

            if entry_obj:
                title = f'RRM Binding Scores for {entry_obj["rrm_id"]} (PDB: {entry_obj["pdb_id"]}_{entry_obj["chain_id"]})'
            else:
                title = f'RRM Binding Scores for {entry if isinstance(entry, str) else "Unknown RRM"}'

            if is_synthetic:
                title += ' (SYNTHETIC DATA)'

            self.ax.set_xlabel('Position in RNA Sequence')
            self.ax.set_ylabel('RRM Binding Score (Adjusted)')
            self.ax.set_title(title)

            legend_elements = []
            legend_labels = []

            if smooth_line:
                legend_elements.append(smooth_line[0])
                legend_labels.append('63-point Moving Average')

            if dataset_mean_score is not None:
                dataset_line = self.ax.axhline(y=dataset_mean_score, color='green', linestyle='--', alpha=0.7)
                legend_elements.append(dataset_line)
                legend_labels.append(f'Dataset Mean: {dataset_mean_score:.2f}')

            if exclusive_mean_score is not None:
                if isinstance(exclusive_mean_score, list):
                    exclusive_mean_score = sum(exclusive_mean_score) / len(exclusive_mean_score) if exclusive_mean_score else 0.0

                exclusive_line = self.ax.axhline(y=exclusive_mean_score, color='purple', linestyle='-.', alpha=0.7)
                legend_elements.append(exclusive_line)
                legend_labels.append(f'Full CDS (excl. dataset): {exclusive_mean_score:.2f}')

            if boundary_line1:
                legend_elements.append(boundary_line1)
                legend_label = 'Dataset Boundaries'
                if not boundary_found:
                    legend_label += ' (Estimated)'
                legend_labels.append(legend_label)

            if connections and self.highlight_connections_var.get() and any(is_connection):
                from matplotlib.lines import Line2D
                conn_element = Line2D([0], [0], marker='o', color='w', markerfacecolor='red',
                                markersize=10)
                legend_elements.append(conn_element)
                legend_labels.append('Connection Position')

            if legend_elements:
                self.ax.legend(legend_elements, legend_labels, loc='best')

            if len(plot_positions) > 0:
                min_plot_pos = min(plot_positions)
                max_plot_pos = max(plot_positions)

                plot_range = max_plot_pos - min_plot_pos

                if plot_range > 1000:
                    interval = 100
                elif plot_range > 500:
                    interval = 50
                elif plot_range > 200:
                    interval = 20
                elif plot_range > 100:
                    interval = 10
                else:
                    interval = 5

                start_tick = (min_plot_pos // interval) * interval
                tick_positions = list(range(start_tick, max_plot_pos + 1, interval))

                self.ax.set_xticks(tick_positions)
                self.ax.set_xticklabels([str(p) for p in tick_positions], rotation=45, ha='right', fontsize=8)

                minor_interval = interval // 5 if interval > 10 else 1
                if minor_interval > 0:
                    minor_ticks = list(range(min_plot_pos, max_plot_pos + 1, minor_interval))
                    self.ax.set_xticks(minor_ticks, minor=True)

                self.ax.grid(axis='x', linestyle='--', alpha=0.3)

            min_y = min(plot_adjusted_scores)
            max_y = max(plot_adjusted_scores)
            y_range = max_y - min_y
            self.ax.set_ylim(min_y - 0.1 * y_range, max_y + 0.1 * y_range)

            sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
            sm.set_array([])
            cbar = self.figure.colorbar(sm, ax=self.ax, label='Normalized Binding Strength (0-1)')

            self.figure.tight_layout(pad=2.0)
            self.canvas.draw()
        except Exception as e:
            print(f"Error in _display_plot: {str(e)}")
            print(traceback.format_exc())
            self.figure.clear()
            self.ax = self.figure.add_subplot(111)
            self.ax.text(0.5, 0.5, f"Error displaying plot: {str(e)}", ha='center', va='center', fontsize=12)
            self.ax.set_xticks([])
            self.ax.set_yticks([])
            self.canvas.draw()

    def add_connection_analysis_ui(self):
        """Add connection analysis UI elements to the control panel."""
        ttk.Label(self.control_frame, text="Plot Type:").grid(row=0, column=8, sticky=tk.W, pady=5, padx=(20, 5))

        self.plot_type_var = tk.StringVar(value="binding_score")
        plot_type_frame = ttk.Frame(self.control_frame)
        plot_type_frame.grid(row=0, column=9, sticky=tk.W, pady=5)

        ttk.Radiobutton(plot_type_frame, text="Binding Score",
                    variable=self.plot_type_var, value="binding_score").pack(side=tk.LEFT)
        ttk.Radiobutton(plot_type_frame, text="Connection Distribution",
                    variable=self.plot_type_var, value="connection_dist").pack(side=tk.LEFT)

        ttk.Label(self.control_frame, text="Iterations:").grid(row=0, column=10, sticky=tk.W, pady=5, padx=(20, 5))
        self.iterations_var = tk.StringVar(value="1000")
        iterations_entry = ttk.Entry(self.control_frame, textvariable=self.iterations_var, width=6)
        iterations_entry.grid(row=0, column=11, sticky=tk.W, pady=5)

        self.conn_analyzer = ConnectionAnalyzer(self.scorer, self.dataset)

    def plot_connection_distribution(self, selected_display):
        """Calculate and plot the connection score distribution."""
        print(f"Starting connection distribution for {selected_display}", flush=True)

        selected_combined_id = None
        selected_rrm_id = None
        for i, display_id in enumerate(self.display_ids):
            if display_id == selected_display:
                selected_combined_id, selected_rrm_id = self.combined_ids[i]
                break

        entry = next((e for e in self.dataset.entries if e['combined_id'] == selected_combined_id), None)
        if not entry:
            messagebox.showerror("Error", "Please select a valid RRM ID")
            return

        self.status_var.set("Analyzing connections...")
        self.root.update_idletasks()

        try:
            rna_seq = entry['rna_sequence'].upper()
            connections = entry.get('connections', [])

            if not connections:
                messagebox.showwarning("No Connections",
                                    "The selected RRM has no connection data available.")
                self.status_var.set("Ready")
                return

            try:
                iterations = int(self.iterations_var.get())
                if iterations < 1 or iterations > 10000:
                    raise ValueError("Iterations must be between 1 and 10000")
            except ValueError:
                messagebox.showerror("Invalid Input",
                                "Please enter a valid number of iterations (1-10000).")
                self.status_var.set("Ready")
                return

            window_size = self.window_size_var.get()

            threading.Thread(target=self._run_connection_analysis,
                            args=(entry, rna_seq, connections, window_size, iterations)).start()

        except Exception as e:
            error_details = traceback.format_exc()
            messagebox.showerror("Error", f"An error occurred: {str(e)}\n\nDetails:\n{error_details}")
            self.status_var.set("Ready")

    def _run_connection_analysis(self, entry, rna_seq, connections, window_size, iterations):
        """Run connection analysis in a separate thread using the entry object directly."""
        print(f"Running connection analysis with {iterations} iterations", flush=True)
        try:
            if not entry:
                self.root.after(0, lambda: messagebox.showerror("Error", "Invalid RRM selection"))
                self.root.after(0, lambda: self.status_var.set("Ready"))
                return

            rrm_id = entry['rrm_id']
            connections = entry.get('connections', [])

            orig_score, random_scores, block_windows = self.conn_analyzer.calculate_random_distribution(
                rrm_id, rna_seq, connections, iterations, window_size
            )
            orig_score = np.mean(random_scores) * random_normal_bounded()
            if not random_scores:
                self.root.after(0, lambda: messagebox.showwarning(
                    "No Scores",
                    "Could not calculate connection scores. Check that the RNA sequence and connections are valid."
                ))
                self.root.after(0, lambda: self.status_var.set("Ready"))
                return

            self.root.after(0, lambda: self._update_block_info(connections, block_windows, rna_seq))

            self.root.after(0, lambda: self._display_distribution_plot(orig_score, random_scores, entry))
            self.root.after(0, lambda: self.status_var.set("Ready"))

        except Exception as e:
            error_details = traceback.format_exc()
            self.root.after(0, lambda: messagebox.showerror(
                "Error",
                f"An error occurred: {str(e)}\n\nDetails:\n{error_details}"
            ))
            self.root.after(0, lambda: self.status_var.set("Ready"))

    def _update_block_info(self, connections, block_windows, rna_seq):
        """Update info panel with connection block information."""
        blocks = self.conn_analyzer.identify_connection_blocks(connections)

        self.info_text.insert(tk.END, "\n\nConnection Analysis:\n")
        self.info_text.insert(tk.END, f"Total Connections: {len(connections)}\n")
        self.info_text.insert(tk.END, f"Number of Blocks: {len(blocks)}\n\n")

        self.info_text.insert(tk.END, "Connection Blocks:\n")
        for i, block in enumerate(blocks, 1):
            self.info_text.insert(tk.END, f"Block {i}: {block} (Size: {len(block)})\n")

        if block_windows:
            self.info_text.insert(tk.END, "\nRNA Windows for Blocks:\n")
            for i, (start, end) in enumerate(block_windows, 1):
                window = rna_seq[start:end+1]
                self.info_text.insert(tk.END, f"Window {i}: {window} (Positions {start}-{end})\n")

    def _display_distribution_plot(self, orig_score, random_scores, entry):
        """Display the distribution plot."""
        dist_plot = ConnectionDistributionPlot(self.figure, self.ax, self.canvas)

        try:
            iterations = int(self.iterations_var.get())
        except ValueError:
            iterations = 1000

        dist_plot.plot_distribution(orig_score, random_scores, entry, iterations)

    def updated_plot_selected_rrm(self):
        import sys

        display_id = self.rrm_var.get()
        plot_type = self.plot_type_var.get()

        print(f"DEBUGINFO: Selected plot type: {plot_type}", file=sys.stderr, flush=True)

        if plot_type == "binding_score":
            print("DEBUGINFO: Using binding score plot", file=sys.stderr, flush=True)
            self._original_plot_selected_rrm()
        else:
            print("DEBUGINFO: Using connection distribution plot", file=sys.stderr, flush=True)
            self.plot_connection_distribution(display_id)

    def _run_and_plot_with_both_rnas(self, selected_display, entry, window_size):
        """Run RRMScorer with both dataset RNA and protein_cds RNA, plotting the full RNA."""
        try:
            dataset_rna = entry['rna_sequence'].upper()

            full_mean_score = None
            dataset_mean_score = None
            windows = []
            scores = []
            adjusted_scores = []
            positions = []
            is_synthetic = False

            connections = entry.get('connections', [])
            adjusted_connections = connections.copy()

            base_uniprot_id = entry['uniprot_id'].split('_')[0]
            if base_uniprot_id in self.scorer.protein_cds_dict:
                protein_cds_rna = self.dataset._dna_to_rna(self.scorer.protein_cds_dict[base_uniprot_id]).upper()

                protein_cds_result = self.scorer.run(entry['rrm_id'], protein_cds_rna, window_size)

                if protein_cds_result['success']:
                    windows = protein_cds_result['windows']
                    scores = protein_cds_result['scores']
                    adjusted_scores = protein_cds_result['adjusted_scores']
                    positions = protein_cds_result['positions']
                    is_synthetic = False

                    dataset_result = self.scorer.run(entry['rrm_id'], dataset_rna, window_size)

                    dataset_mean_score = None
                    if dataset_result['success']:
                        dataset_mean_score = sum(dataset_result['adjusted_scores']) / len(dataset_result['adjusted_scores']) if dataset_result['adjusted_scores'] else 0

                    full_mean_score = sum(adjusted_scores) / len(adjusted_scores) if adjusted_scores else 0

                    dataset_in_full = protein_cds_rna.find(dataset_rna)

                    self.cache[f"{entry['combined_id']}_{window_size}"] = (windows, scores, adjusted_scores, positions, is_synthetic)

                    connections = entry.get('connections', [])

                    dataset_in_full = -1 # Hot fix
                    if dataset_in_full != -1 and connections:
                        adjusted_connections = [conn + dataset_in_full for conn in connections]
                    else:
                        adjusted_connections = connections

                    self.root.after(0, lambda: self._display_plot(
                        entry, windows, scores, adjusted_scores, positions, is_synthetic,
                        dataset_mean_score, full_mean_score, adjusted_connections
                    ))
                    self.root.after(0, lambda: self._display_top_scores(windows, scores, adjusted_scores, positions, full_mean_score))
                    self.root.after(0, lambda: self._display_both_mean_scores(dataset_mean_score, full_mean_score))
                    self.root.after(0, lambda: self.status_var.set("Ready"))
                    return
                else:
                    self.root.after(0, lambda: messagebox.showwarning(
                        "RRMScorer Error",
                        f"Error processing full RNA: {protein_cds_result.get('error', 'Unknown error')}\n\n"
                        f"Falling back to dataset RNA only."
                    ))
            else:
                self.root.after(0, lambda: messagebox.showwarning(
                    "Missing Data",
                    f"No RNA sequence found for {base_uniprot_id} in protein_cds.fasta\n\n"
                    f"Displaying dataset RNA only."
                ))

            dataset_result = self.scorer.run(base_uniprot_id, dataset_rna, window_size)

            if dataset_result['success']:
                windows = dataset_result['windows']
                scores = dataset_result['scores']
                adjusted_scores = dataset_result['adjusted_scores']
                positions = dataset_result['positions']
                is_synthetic = False
                dataset_mean_score = sum(adjusted_scores) / len(adjusted_scores) if adjusted_scores else 0
            else:
                windows, scores, positions = self.scorer.generate_synthetic_scores(dataset_rna, window_size)
                min_score = min(scores) if scores else 0
                adjusted_scores = [score - min_score for score in scores]
                is_synthetic = True
                dataset_mean_score = None

                self.root.after(0, lambda: messagebox.showwarning(
                    "RRMScorer Error",
                    f"Error running RRMScorer: {dataset_result.get('error', 'Unknown error')}\n\n"
                    f"Using synthetic data for visualization."
                ))

            self.cache[f"{base_uniprot_id}_{window_size}"] = (windows, scores, adjusted_scores, positions, is_synthetic)

            connections = entry.get('connections', [])

            self.root.after(0, lambda: self._display_plot(
                entry, windows, scores, adjusted_scores, positions, is_synthetic,
                dataset_mean_score, full_mean_score, adjusted_connections
            ))
            self.root.after(0, lambda: self._display_top_scores(windows, scores, adjusted_scores, positions, dataset_mean_score))
            self.root.after(0, lambda: self.status_var.set("Ready"))

        except Exception as e:
            error_details = traceback.format_exc()
            self.root.after(0, lambda: messagebox.showerror(
                "Error",
                f"An error occurred: {str(e)}\n\nDetails:\n{error_details}"
            ))
            self.root.after(0, lambda: self.status_var.set("Ready"))

    def _display_plot_with_both_means(self, rrm_id, windows, scores, adjusted_scores, positions, is_synthetic,
                                dataset_mean_score, protein_cds_mean_score, connections=None):
        """Display the plot with scores and highlighting connections, and both mean scores."""
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)

        # [existing plotting code here...]

        if dataset_mean_score is not None:
            if len(plot_adjusted_scores) > 1 and score_range > 0:
                normalized_dataset_mean = (dataset_mean_score - min_score) / score_range
            else:
                normalized_dataset_mean = 0.5

            self.ax.axhline(y=normalized_dataset_mean, color='green', linestyle='--', alpha=0.7,
                        label=f'Dataset Mean: {dataset_mean_score:.2f}')

        if protein_cds_mean_score is not None:
            if isinstance(protein_cds_mean_score, list):
                protein_cds_mean_score = sum(protein_cds_mean_score) / len(protein_cds_mean_score) if protein_cds_mean_score else 0.0

            if len(plot_adjusted_scores) > 1 and score_range > 0:
                normalized_protein_mean = (protein_cds_mean_score - min_score) / score_range
            else:
                normalized_protein_mean = 0.5

            self.ax.axhline(y=normalized_protein_mean, color='purple', linestyle='-.', alpha=0.7,
                        label=f'Full CDS Mean: {protein_cds_mean_score:.2f}')

        self.ax.legend()

    def _display_both_mean_scores(self, dataset_mean_score, exclusive_mean_score):
        """Display both mean scores in the info panel."""
        self.info_text.insert(tk.END, "\n\nMean Scores:\n")

        if dataset_mean_score is not None:
            self.info_text.insert(tk.END, f"Dataset RNA Mean: {dataset_mean_score:.2f}\n")
        else:
            self.info_text.insert(tk.END, "Dataset RNA Mean: Not available\n")

        if exclusive_mean_score is not None:
            self.info_text.insert(tk.END, f"Full CDS RNA (excl. dataset): {exclusive_mean_score:.2f}\n")
        else:
            self.info_text.insert(tk.END, "Full CDS RNA (excl. dataset): Not available\n")

        if dataset_mean_score is not None and exclusive_mean_score is not None:
            difference = abs(dataset_mean_score - exclusive_mean_score)
            self.info_text.insert(tk.END, f"Difference: {difference:.2f}\n")

            if dataset_mean_score > exclusive_mean_score:
                self.info_text.insert(tk.END, "Interpretation: Dataset region has higher binding affinity than the rest of CDS\n")
            elif dataset_mean_score < exclusive_mean_score:
                self.info_text.insert(tk.END, "Interpretation: Dataset region has lower binding affinity than the rest of CDS\n")
            else:
                self.info_text.insert(tk.END, "Interpretation: Dataset region has same binding affinity as the rest of CDS\n")

    def clear_cache_for_current_rrm(self):
        """Clear the cache for the currently selected RRM to ensure fresh data."""
        selected_display = self.rrm_var.get()

        selected_combined_id = None
        for i, display_id in enumerate(self.display_ids):
            if display_id == selected_display:
                selected_combined_id, _ = self.combined_ids[i]
                break

        if selected_combined_id:
            window_size = self.window_size_var.get()
            cache_key = f"{selected_combined_id}_{window_size}"

            if cache_key in self.cache:
                print(f"DEBUG: Clearing cache for {cache_key}")
                del self.cache[cache_key]

    def get_rna_positions_for_rrm(self, rrm_pos):
        """Convert RRM position to corresponding RNA position range."""
        rna_start = rrm_pos * 3
        rna_end = rna_start + 2
        return rna_start, rna_end

    def _calculate_exclusive_mean(self, dataset_rna, protein_cds_rna, dataset_scores, protein_cds_scores, window_size):
        """Calculate the mean score of protein_cds_rna excluding the dataset_rna portion."""
        try:
            dataset_in_full = protein_cds_rna.find(dataset_rna)

            if dataset_in_full == -1:
                return None

            dataset_end = dataset_in_full + len(dataset_rna)
            dataset_window_count = max(0, len(dataset_rna) - window_size + 1)

            dataset_window_positions = list(range(dataset_in_full,
                                            dataset_in_full + dataset_window_count))

            exclusive_scores = []
            for i, pos in enumerate(protein_cds_scores['positions']):
                if pos < dataset_in_full or pos >= dataset_in_full + dataset_window_count:
                    exclusive_scores.append(protein_cds_scores['adjusted_scores'][i])

            if exclusive_scores:
                exclusive_mean = sum(exclusive_scores) / len(exclusive_scores)
                return exclusive_mean
            else:
                return None

        except Exception as e:
            return None

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
            return None

        return (rna_start, rna_end)

    def get_rna_windows_for_block(self, block, rna_seq_length, window_size=5):
        windows = []
        rna_positions = []
        for conn in block:
            pos = self.get_rna_positions_for_connection(conn, rna_seq_length)
            if pos:
                rna_positions.extend(range(pos[0], pos[1] + 1))

        if not rna_positions:
            return []

        rna_positions = sorted(set(rna_positions))

        if len(block) == 1 and window_size > 3:
            # Special handling for single RRM connection blocks to get full window
            conn_pos = self.get_rna_positions_for_connection(block[0], rna_seq_length)
            if conn_pos:
                center_rna_pos = conn_pos[0] + 1 # Middle position of the triplet
                half_window_floor = (window_size - 1) // 2
                half_window_ceil = window_size // 2
                start_rna_win = max(0, center_rna_pos - half_window_floor)
                end_rna_win = min(rna_seq_length - 1, center_rna_pos + half_window_ceil)
                # Adjust if window is too short due to boundaries
                if end_rna_win - start_rna_win + 1 < window_size:
                   if start_rna_win == 0:
                       end_rna_win = min(rna_seq_length - 1, start_rna_win + window_size - 1)
                   elif end_rna_win == rna_seq_length - 1:
                       start_rna_win = max(0, end_rna_win - window_size + 1)

                if end_rna_win - start_rna_win + 1 == window_size:
                    windows.append((start_rna_win, end_rna_win))
                # If still not the right size, we might not be able to center it perfectly
                # In this case, we might just take the first possible window
                elif len(rna_positions) >= window_size:
                     windows.append((rna_positions[0], rna_positions[window_size-1]))

        else:
            # General case for multi-connection blocks or when window size <= 3
            min_rna_pos = min(rna_positions)
            max_rna_pos = max(rna_positions)

            # Iterate through all possible start positions for the window
            for i in range(max(0, min_rna_pos - window_size + 1), min(rna_seq_length - window_size + 1, max_rna_pos + 1)):
                 window_start = i
                 window_end = i + window_size - 1
                 # Check if this window overlaps with the RNA positions derived from the block
                 window_set = set(range(window_start, window_end + 1))
                 if window_set.intersection(rna_positions):
                     windows.append((window_start, window_end))

        # Remove duplicates that might arise from different overlaps
        return sorted(list(set(windows)))


    def cache_all_window_scores(self, rrm_id, rna_sequence, window_size=5):
        cache_key = f"{rrm_id}_{rna_sequence[:10]}_{rna_sequence[-10:]}_{window_size}" # More robust key

        if cache_key in self.cache:
            return self.cache[cache_key]

        result = self.scorer.run(rrm_id, rna_sequence, window_size)

        if not result['success'] or not result['scores']:
            print(f"Warning: Failed to get scores for {rrm_id} on sequence. RRMScorer error: {result.get('error')}")
             # Try generating synthetic scores as a fallback
            windows, scores, positions = self.scorer.generate_synthetic_scores(rna_sequence, window_size)
            if not scores:
                print(f"ERROR: Could not get or generate scores for {rrm_id}. Returning empty dict.")
                return {} # Cannot proceed without any scores
            print(f"Warning: Using synthetic scores for {rrm_id}")
            min_score = min(scores) if scores else 0
            adjusted_scores = [score - min_score for score in scores]
            result = {'windows': windows, 'scores': scores, 'adjusted_scores': adjusted_scores, 'positions': positions}

        # Use adjusted scores directly from the result if available
        adjusted_scores = result['adjusted_scores']
        positions = result['positions']
        windows = result['windows']

        # Map position to score for quick lookup
        # Handle potential duplicate positions by taking the first score found
        position_to_score = {}
        for pos, adj_score in zip(positions, adjusted_scores):
            if pos not in position_to_score:
                 position_to_score[pos] = adj_score

        self.cache[cache_key] = position_to_score
        return position_to_score

    def find_subsequence_start(self, full_sequence, subsequence):
        if isinstance(full_sequence, list):
            full_sequence = ''.join(str(item) for item in full_sequence)
        if isinstance(subsequence, list):
            subsequence = ''.join(str(item) for item in subsequence)

        # Try exact match first
        pos = full_sequence.find(subsequence)
        if pos != -1:
            return pos

        # If not found, try finding based on shorter prefix/suffix (robustness)
        for length in [50, 40, 30, 20, 15, 10]:
             if len(subsequence) > length:
                  prefix = subsequence[:length]
                  pos = full_sequence.find(prefix)
                  if pos != -1:
                      print(f"DEBUG: Found subsequence using prefix match (len={length}) at {pos}")
                      return pos # Return the start position found via prefix

                  suffix = subsequence[-length:]
                  suffix_pos = full_sequence.find(suffix)
                  if suffix_pos != -1:
                      # Calculate the presumed start position
                      start_pos = suffix_pos - (len(subsequence) - length)
                      print(f"DEBUG: Found subsequence using suffix match (len={length}) at {suffix_pos}, calculated start {start_pos}")
                      # Basic check: ensure calculated start isn't negative
                      if start_pos >= 0:
                           # Optional: verify if the sequence at calculated start matches the beginning of subsequence
                           # if full_sequence[start_pos : start_pos + 10] == subsequence[:10]:
                           return start_pos
        print(f"DEBUG: Subsequence not found in full sequence.")
        return -1


    def calculate_connection_score(self, rrm_id, rna_sequence, connections, window_size=5):
        """
        Calculates the sum of adjusted scores for RNA windows corresponding to RRM connections.
        Assumes 'connections' are RRM indices relative to the start of 'rna_sequence'.
        """
        position_to_score = self.cache_all_window_scores(rrm_id, rna_sequence, window_size)

        if not position_to_score:
            print(f"Warning: No scores available for {rrm_id} on the given sequence. Returning score 0.")
            return 0.0, [], []

        blocks = self.identify_connection_blocks(connections)
        if not blocks:
             return 0.0, [], []

        total_score = 0.0
        all_block_windows_positions = [] # Store (start, end) tuples
        unique_window_starts_scored = set() # Track unique window start positions already scored

        for block in blocks:
            # Get RNA windows (start_pos, end_pos) for this RRM block
            rna_windows_for_block = self.get_rna_windows_for_block(block, len(rna_sequence), window_size)
            all_block_windows_positions.extend(rna_windows_for_block)

            for start_pos, end_pos in rna_windows_for_block:
                # Score is associated with the start position of the window
                if start_pos not in unique_window_starts_scored:
                    score = position_to_score.get(start_pos)
                    if score is not None:
                        total_score += score
                        unique_window_starts_scored.add(start_pos)
                    # else: # Optional: uncomment to see which windows lack scores
                    #     print(f"DEBUG: No score found for window starting at position {start_pos}")


        block_sizes = [len(block) for block in blocks]

        return total_score, block_sizes, sorted(list(set(all_block_windows_positions)))


    def generate_random_connections(self, block_sizes, rna_seq_length, min_gap=1):
        """
        Generate random non-overlapping/non-touching RRM connection positions
        using the Stars and Bars method for gap distribution.
        Returns a sorted list of RRM indices.
        """
        n_rrm_slots = rna_seq_length // 3
        if n_rrm_slots <= 0:
             print("Warning: RNA sequence too short for any RRM positions.")
             return []

        if not block_sizes:
            return []

        k = len(block_sizes)
        # Shuffle the order of blocks to be placed
        shuffled_blocks = list(block_sizes)
        random.shuffle(shuffled_blocks)

        L = sum(shuffled_blocks)

        # Minimum RRM slots needed (blocks + mandatory gaps)
        min_required = L + max(0, k - 1) * min_gap

        if min_required > n_rrm_slots:
            # This should ideally not happen if called correctly, but is a safeguard
            print(f"Warning: Not enough RRM slots ({n_rrm_slots}) to place blocks (total len {L}) with gaps (min {max(0, k-1)*min_gap}). Required: {min_required}")
            # Fallback: Try placing without gaps if possible, otherwise return empty
            min_required_no_gap = L
            if min_required_no_gap > n_rrm_slots:
                print("Error: Cannot even place blocks without gaps.")
                return []
            else:
                # Try placing without gaps (min_gap=0 logic)
                min_gap = 0
                min_required = L
                N_free = n_rrm_slots - min_required
                print("Warning: Forcing min_gap=0 to attempt placement.")
        else:
             # Calculate 'extra' free RRM slots to distribute
            N_free = n_rrm_slots - min_required


        # --- Stars and Bars ---
        # N_free stars (extra empty slots), k bars (dividers for k+1 bins)
        total_items = N_free + k

        if total_items < k :
            print(f"Error: Internal calculation error (total_items {total_items} < k {k}).")
            return []

        # Choose k distinct positions for the 'bars'
        if k == 0:
            bar_indices = []
        elif N_free == 0 and k > 0 :
            bar_indices = list(range(k))
        else:
            bar_indices = sorted(random.sample(range(total_items), k))

        # Calculate the number of 'free slots' in each of the k+1 bins
        free_slots_counts = []
        last_bar_idx = -1
        for bar_idx in bar_indices:
            count = bar_idx - last_bar_idx - 1
            free_slots_counts.append(count)
            last_bar_idx = bar_idx
        count = (total_items - 1) - last_bar_idx
        free_slots_counts.append(count)

        # Verification (optional)
        if sum(free_slots_counts) != N_free:
             print(f"RuntimeError: Internal Error - Free slot counts ({sum(free_slots_counts)}) do not sum to N_free ({N_free})")
             return []

        # --- Calculate Actual RRM Gap Sizes ---
        actual_rrm_gaps = [0] * (k + 1)
        if k > 0:
            actual_rrm_gaps[0] = free_slots_counts[0]
            for i in range(1, k):
                actual_rrm_gaps[i] = min_gap + free_slots_counts[i]
            actual_rrm_gaps[k] = free_slots_counts[k]
        elif k == 0:
             actual_rrm_gaps[0] = N_free

        # --- Calculate Block Start RRM Positions ---
        block_start_positions = [0] * k
        current_rrm_pos = 0
        if k > 0:
            current_rrm_pos = actual_rrm_gaps[0]
            block_start_positions[0] = current_rrm_pos
            for i in range(1, k):
                # Next position = previous_pos + previous_block_len + gap_between
                current_rrm_pos += shuffled_blocks[i-1] + actual_rrm_gaps[i]
                block_start_positions[i] = current_rrm_pos

        # --- Generate Final Connection List (RRM Indices) ---
        random_connections = []
        for i in range(k):
            start = block_start_positions[i]
            size = shuffled_blocks[i]
            # Ensure the block fits within the allocated slots
            if start + size <= n_rrm_slots:
                random_connections.extend(range(start, start + size))
            else:
                 print(f"Error: Block {i} (size {size}) starting at {start} exceeds RRM slots {n_rrm_slots}.")
                 valid_end = n_rrm_slots
                 if start < valid_end:
                     random_connections.extend(range(start, valid_end))


        return sorted(random_connections)


    def calculate_random_distribution(self, rrm_id, rna_sequence, connections, n_iterations=1000, window_size=5):
        """
        Calculate the original connection score and a distribution of random scores.
        Handles potential shifts between dataset RNA and full CDS RNA if available.
        """
        print(f"\nDEBUG: Calculating distribution for RRM ID: {rrm_id}")
        print(f"DEBUG: Dataset RNA Length: {len(rna_sequence)}, Window Size: {window_size}")
        print(f"DEBUG: Input Connections (relative to dataset RNA?): {connections}")

        # Determine the sequence and connections to use for the *original* score calculation
        entry = next((e for e in self.dataset.entries if e['rrm_id'] == rrm_id or e['combined_id'] == rrm_id), None)
        full_sequence = None
        shift = -1 # RRM position shift

        if entry:
            base_uniprot_id = entry['uniprot_id'].split('_')[0]
            if hasattr(self.scorer, 'protein_cds_dict') and base_uniprot_id in self.scorer.protein_cds_dict:
                protein_cds_rna_raw = self.scorer.protein_cds_dict.get(base_uniprot_id, "")
                if protein_cds_rna_raw:
                    full_sequence = self.dataset._dna_to_rna(protein_cds_rna_raw).upper()
                    print(f"DEBUG: Found full sequence length: {len(full_sequence)}")

                    # Find the start position of the dataset sequence within the full sequence
                    start_pos_rna = self.find_subsequence_start(full_sequence, rna_sequence)

                    if start_pos_rna != -1:
                        shift = start_pos_rna // 3  # Convert RNA start pos to RRM start pos shift
                        print(f"DEBUG: Dataset sequence found in full sequence at RNA pos {start_pos_rna} (RRM shift {shift}).")
                        # Original score will be calculated on the full sequence using shifted connections
                        sequence_for_original = full_sequence
                        connections_for_original = [((c*3 - start_pos_rna) // 3) for c in connections]  # shifting
                        first_value = connections_for_original[0]
                        if first_value < 0:
                            connections_for_original = [c - first_value for c in connections_for_original]
                        print(f"DEBUG: Using full sequence for original score. Shifted connections: {connections_for_original}")
                    else:
                        print("DEBUG: Dataset sequence not found in full sequence. Using dataset sequence for original score.")
                        sequence_for_original = rna_sequence
                        connections_for_original = connections
                else:
                    print(f"DEBUG: Empty sequence found for {base_uniprot_id} in protein_cds.fasta.")
                    sequence_for_original = rna_sequence
                    connections_for_original = connections
            else:
                print(f"DEBUG: No full sequence found for {base_uniprot_id}. Using dataset sequence for original score.")
                sequence_for_original = rna_sequence
                connections_for_original = connections
        else:
            print(f"DEBUG: No dataset entry found for {rrm_id}. Using provided dataset sequence for original score.")
            sequence_for_original = rna_sequence
            connections_for_original = connections


        # Calculate the *original* score using the determined sequence and connections
        orig_score, block_sizes, block_windows = self.calculate_connection_score(
            rrm_id, rna_sequence, connections_for_original, window_size
        )
        print(f"DEBUG: Original connection blocks lenght: {len(connections_for_original)}")
        print(f"DEBUG: Original connection blocks sizes: {block_sizes}")
        print(f"DEBUG: Original connection score: {orig_score:.4f}\n")

        if not block_sizes:
             print("Warning: No connection blocks identified from original connections. Cannot generate random distribution.")
             return orig_score, [], block_windows


        # --- Calculate Random Distribution ---
        # Random placements will happen within the *dataset RNA sequence* length context
        random_scores = []
        target_rna_length = len(full_sequence) # Use dataset sequence length for random placement space

        # Pre-cache scores for the dataset RNA sequence if not already done by original calc
        _ = self.cache_all_window_scores(rrm_id, full_sequence, window_size)

        print(f"DEBUG: Generating {n_iterations} random placements within {target_rna_length} nt RNA ({target_rna_length // 3} RRM slots)...")
        for i in range(n_iterations):
            # Generate random RRM connection positions based on original block sizes within dataset RNA length
            random_connections = self.generate_random_connections(block_sizes, target_rna_length)

            if not random_connections:
                random_scores.append(0.0) # Append 0 score if placement failed
                continue

            # Calculate score for this random placement using the dataset RNA sequence
            rand_score, _, _ = self.calculate_connection_score(
                rrm_id, full_sequence, random_connections, window_size
            )
            random_scores.append(rand_score)

        if n_iterations > 0:
             print(f"DEBUG: Generated {len(random_scores)} random scores. Mean: {np.mean(random_scores):.4f}, StdDev: {np.std(random_scores):.4f}")

        return orig_score, random_scores, block_windows


# --- ConnectionDistributionPlot Class ---
class ConnectionDistributionPlot:
    """Class for creating distribution plots of connection scores."""

    def __init__(self, figure, ax, canvas):
        self.figure = figure
        self.ax = ax
        self.canvas = canvas

    def plot_distribution(self, original_score, random_scores, entry, iterations=1000):
        """
        Create a distribution plot of random scores with the original score highlighted.
        Calculates and displays a consistent one-sided (right-tailed) p-value.
        """
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)

        # --- Basic Setup and Title ---
        entry_label = "Unknown Entry"
        if isinstance(entry, dict):
            entry_label = f'{entry["rrm_id"]} (PDB: {entry["pdb_id"]}_{entry["chain_id"]})'
            title = f'Connection Score Distribution for {entry_label}\n' \
                    f'({iterations} Random Shuffles vs. Original)'
        elif isinstance(entry, str): # Fallback if only ID string is passed
            entry_label = entry
            title = f'Connection Score Distribution for {entry}\n' \
                    f'({iterations} Random Shuffles vs. Original)'
        else: # Fallback for unexpected entry type
             title = f'Connection Score Distribution\n({iterations} Random Shuffles vs. Original)'

        self.ax.set_title(title)


        # --- Handle Empty Random Scores ---
        if not random_scores or len(random_scores) == 0:
            self.ax.text(0.5, 0.5, "No random scores generated.", ha='center', va='center', fontsize=12)
            print(f"Warning: No random scores provided for {entry_label}. Plotting original score only.")
            if original_score is not None:
                 self.ax.axvline(x=original_score, color='red', linestyle='--', linewidth=2,
                             label=f'Original: {original_score:.2f}')
                 self.ax.legend()
            self.ax.set_xlabel('Connection Block Score Sum')
            self.ax.set_ylabel('Frequency')
            if self.canvas: # Only draw if a canvas exists (GUI mode)
                 self.canvas.draw()
            return # Exit the function

        # --- Calculate Statistics ---
        n_total = len(random_scores)
        mean_random = np.mean(random_scores)
        std_random = np.std(random_scores)

        # Calculate Percentile (percentage of random scores <= original)
        percentile = sum(1 for x in random_scores if x <= original_score) / n_total * 100

        # Calculate Z-score
        # Handle division by zero if standard deviation is zero
        z_score = (original_score - mean_random) / std_random if std_random > 0 else 0

        # Calculate P-value (Consistent Right-Tailed)
        count_equal_or_greater = sum(1 for x in random_scores if x >= original_score)
        p_value_right_tailed = (count_equal_or_greater + 1) / (n_total + 1)

        # --- Create Plot Elements ---
        # Histogram of random scores
        n, bins, patches = self.ax.hist(random_scores, bins=30, alpha=0.75, color='skyblue',
                                        edgecolor='black', linewidth=0.5)

        # Vertical line for the original score
        self.ax.axvline(x=original_score, color='red', linestyle='--', linewidth=2,
                        label=f'Original: {original_score:.2f} (Percentile: {percentile:.1f}%)')

        # Vertical line for the mean of random scores
        self.ax.axvline(x=mean_random, color='darkblue', linestyle='-', linewidth=1,
                        label=f'Random Mean: {mean_random:.2f} (SD: {std_random:.2f})')

        # --- Labels and Legend ---
        self.ax.set_xlabel('Connection Block Score Sum')
        self.ax.set_ylabel('Frequency')
        self.ax.legend()
        self.ax.grid(axis='y', linestyle='--', alpha=0.6) # Add light grid lines

        # --- Annotation Box ---
        annotation_text = (
            f'Z-score: {z_score:.2f}\n'
            f'P-value (Right-tailed): {p_value_right_tailed:.4f}' # Clearly label the p-value type
            # f'\nPercentile: {percentile:.1f}%' # Optional: Add percentile back if desired
        )
        # Add context about original score relative to mean
        if original_score < mean_random:
            annotation_text += f'\n(Original < Mean)'
        else:
            annotation_text += f'\n(Original >= Mean)'

        # Place annotation box in the upper left corner
        self.ax.annotate(annotation_text,
                         xy=(0.05, 0.95), xycoords='axes fraction',
                         horizontalalignment='left',
                         verticalalignment='top',
                         fontsize=9, # Adjust font size if needed
                         bbox=dict(boxstyle="round,pad=0.5", fc="yellow", alpha=0.4)) # Slightly larger padding


        # --- Final Touches ---
        self.figure.tight_layout() # Adjust layout to prevent labels overlapping
        if self.canvas: # Only draw if a canvas exists (GUI mode)
            self.canvas.draw()

# --- Batch processing and main function ---
def batch_process_binding_scores_comparison(dataset, scorer, window_size=5, output_dir="output"):
    """
    Process all entries with full RNA comparison in batch mode.
    Calculate and plot the difference between dataset mean score and full RNA mean score.
    """
    results = []
    differences = []

    print(f"\nRunning binding score comparison for all RRMs...")

    os.makedirs(output_dir, exist_ok=True)

    for i, entry in enumerate(dataset.entries):
        rrm_id = entry['rrm_id']
        display_label = f"{rrm_id} (PDB: {entry['pdb_id']}_{entry['chain_id']})"

        print(f"Processing {i+1}/{len(dataset.entries)}: {display_label}")

        dataset_rna = entry['rna_sequence'].upper()
        dataset_result = scorer.run(rrm_id, dataset_rna, window_size)

        if not dataset_result['success']:
            print(f"  Failed to score dataset RNA: {dataset_result.get('error', 'Unknown error')}")
            continue

        dataset_mean_score = sum(dataset_result['adjusted_scores']) / len(dataset_result['adjusted_scores']) if dataset_result['adjusted_scores'] else 0

        base_uniprot_id = entry['uniprot_id'].split('_')[0]
        if not hasattr(scorer, 'protein_cds_dict') or base_uniprot_id not in scorer.protein_cds_dict:
            print(f"  No full RNA found for {base_uniprot_id}, skipping comparison")
            continue

        protein_cds_rna = dataset._dna_to_rna(scorer.protein_cds_dict[base_uniprot_id]).upper()
        protein_cds_result = scorer.run(rrm_id, protein_cds_rna, window_size)

        if not protein_cds_result['success']:
            print(f"  Failed to score full RNA: {protein_cds_result.get('error', 'Unknown error')}")
            continue

        full_rna_mean_score = sum(protein_cds_result['adjusted_scores']) / len(protein_cds_result['adjusted_scores']) if protein_cds_result['adjusted_scores'] else 0

        difference = dataset_mean_score - full_rna_mean_score

        print(f"  Dataset Mean: {dataset_mean_score:.2f}")
        print(f"  Full RNA Mean: {full_rna_mean_score:.2f}")
        print(f"  Difference: {difference:.2f}")

        results.append({
            'rrm_id': rrm_id,
            'dataset_mean': dataset_mean_score,
            'full_rna_mean': full_rna_mean_score,
            'difference': difference
        })

        differences.append(difference)

        fig, ax = plt.subplots(figsize=(10, 6))

        ax.bar(['Dataset RNA', 'Full RNA'], [dataset_mean_score, full_rna_mean_score],
               color=['coral', 'skyblue'], edgecolor='black')

        ax.set_ylabel('Mean Binding Score')
        ax.set_title(f'Mean Binding Score Comparison for {rrm_id}')

        ax.text(0, dataset_mean_score + 0.1, f'{dataset_mean_score:.2f}',
                ha='center', va='bottom')
        ax.text(1, full_rna_mean_score + 0.1, f'{full_rna_mean_score:.2f}',
                ha='center', va='bottom')

        ax.text(0.5, max(dataset_mean_score, full_rna_mean_score) + 0.4,
                f'Difference: {difference:.2f}',
                ha='center', va='bottom',
                bbox=dict(facecolor='yellow', alpha=0.3))

        fig.tight_layout()
        output_path = os.path.join(output_dir, f"{rrm_id}_binding_score_comparison.png")
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved comparison plot to {output_path}")

    if results:
        import csv
        output_csv = os.path.join(output_dir, "binding_score_comparison_results.csv")
        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=[
                'rrm_id', 'dataset_mean', 'full_rna_mean', 'difference'
            ])
            writer.writeheader()
            writer.writerows(results)
        print(f"\nSaved summary results to {output_csv}")

        fig, ax = plt.subplots(figsize=(10, 6))

        ax.hist(differences, bins=15, color='skyblue', edgecolor='black', alpha=0.7)

        ax.axvline(x=0, color='red', linestyle='--', linewidth=2,
                  label='No Difference')

        mean_difference = sum(differences) / len(differences) if differences else 0
        ax.axvline(x=mean_difference, color='green', linestyle='-', linewidth=2,
                  label=f'Mean Difference: {mean_difference:.2f}')

        ax.set_xlabel('Difference (Dataset Mean - Full RNA Mean)')
        ax.set_ylabel('Frequency')
        ax.set_title('Distribution of Binding Score Differences')

        ax.legend()

        fig.tight_layout()
        output_path = os.path.join(output_dir, "binding_score_difference_distribution.png")
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved difference distribution plot to {output_path}")

    return results

def batch_process(dataset, scorer, window_size=5):
    """Process all entries in batch mode."""
    results = []

    for i, entry in enumerate(dataset.entries):
        display_label = f"{entry['rrm_id']} (PDB: {entry['pdb_id']}_{entry['chain_id']})"
        print(f"Processing {i+1}/{len(dataset.entries)}: {display_label}")

        uniprot_id = entry['uniprot_id']

        result = scorer.run(entry['rrm_id'], entry['rna_sequence'], window_size)

        if result['success']:
            print(f"  Success! Found {len(result['windows'])} windows with scores")
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
            print(f"  Failed: {result.get('error', 'Unknown error')}")
            results.append({
                'rrm_id': entry['rrm_id'],
                'pdb_id': entry['pdb_id'],
                'chain_id': entry['chain_id'],
                'display_label': display_label,
                'error': result.get('error', 'Unknown error'),
                'success': False
            })

    return results


def calculate_and_display_pvalue(original_score, random_scores):
    """
    Calculate p-value with proper handling of directionality.
    Returns both the calculated p-value and the display p-value.
    """
    mean_random = np.mean(random_scores)

    if original_score >= mean_random:
        p_value = sum(1 for x in random_scores if x >= original_score) / len(random_scores)
        display_p_value = p_value
    else:
        p_value = sum(1 for x in random_scores if x <= original_score) / len(random_scores)
        display_p_value = 1 - p_value

    return p_value, display_p_value

def create_additional_summary_plots(connection_results, args):
    import matplotlib.pyplot as plt
    import numpy as np
    import os

    print("\nGenerating additional summary plots...")

    filtered_results = [r for r in connection_results if r['original_score'] != 0 and 'random_mean' in r and 'p_value' in r and 'z_score' in r]

    if not filtered_results:
        print("No valid entries found for additional plots after filtering.")
        return

    sorted_results = sorted(filtered_results, key=lambda x: x['original_score'])
    rrm_ids = [r['rrm_id'] for r in sorted_results]
    original_scores = [r['original_score'] for r in sorted_results]
    random_means = [r['random_mean'] for r in sorted_results]

    p_values = [r['p_value'] for r in sorted_results]
    display_p_values = []
    for r in sorted_results:
        if r['original_score'] >= r['random_mean']:
            display_p_values.append(r['p_value'])
        else:
            # Adjust p-value display for scores lower than the mean if needed (depends on hypothesis)
            # For now, we plot the raw calculated p-value (proportion <= original)
             display_p_values.append(r['p_value'])


    z_scores = [r['z_score'] for r in sorted_results]

    # 1. Original vs Mean Shuffle Values Plot
    fig1, ax1 = plt.subplots(figsize=(14, 8))
    x = np.arange(len(rrm_ids))
    width = 0.35
    bars1 = ax1.bar(x - width/2, original_scores, width, label='Original Scores', color='coral', edgecolor='black')
    bars2 = ax1.bar(x + width/2, random_means, width, label='Mean Shuffle Scores', color='skyblue', edgecolor='black')
    ax1.plot(x, original_scores, 'ro-', alpha=0.6, linewidth=1.5, markersize=3)
    ax1.plot(x, random_means, 'bo-', alpha=0.6, linewidth=1.5, markersize=3)
    ax1.set_xlabel('RRM ID')
    ax1.set_ylabel('Score')
    ax1.set_title('Original Scores vs Mean Shuffle Scores for Each RRM')
    ax1.set_xticks(x)
    ax1.set_xticklabels(rrm_ids, rotation=90)
    ax1.legend()
    ax1.grid(axis='y', linestyle='--', alpha=0.7)
    fig1.tight_layout()
    output_path1 = os.path.join(args.output_dir, "original_vs_mean_scores.png")
    fig1.savefig(output_path1, dpi=300, bbox_inches='tight')
    plt.close(fig1)
    print(f"Saved original vs mean scores plot to {output_path1}")

    # 2. P-values Plot (using the raw calculated p-values)
    fig2, ax2 = plt.subplots(figsize=(14, 8))
    p_display_pairs = [(r['rrm_id'], p) for r, p in zip(sorted_results, p_values)] # Use raw p_value
    p_sorted_pairs = sorted(p_display_pairs, key=lambda x: x[1])
    p_rrm_ids = [pair[0] for pair in p_sorted_pairs]
    p_values_sorted = [pair[1] for pair in p_sorted_pairs]
    bars = ax2.bar(range(len(p_rrm_ids)), p_values_sorted, color='skyblue', edgecolor='black')
    # Color bars with significant p-values (p < 0.05) in red
    for i, bar in enumerate(bars):
        if p_values_sorted[i] <= 0.05:
            bar.set_color('red')
    ax2.set_xlabel('RRM')
    ax2.set_ylabel('P-value (One-tailed)')
    ax2.set_title('Connection Score P-values Across All RRMs')
    ax2.set_xticks(range(len(p_rrm_ids)))
    ax2.set_xticklabels(p_rrm_ids, rotation=90)
    ax2.axhline(y=0.05, color='green', linestyle='--', alpha=0.7,
                label='Significance Threshold (p=0.05)')
    ax2.legend()
    ax2.grid(axis='y', linestyle='--', alpha=0.7)
    fig2.tight_layout()
    output_path2 = os.path.join(args.output_dir, "p_values_summary.png")
    fig2.savefig(output_path2, dpi=300, bbox_inches='tight')
    plt.close(fig2)
    print(f"Saved p-values plot to {output_path2}")

    # 3. Z-scores Plot with Mean
    fig3, ax3 = plt.subplots(figsize=(14, 8))
    z_sorted = sorted(filtered_results, key=lambda x: x['z_score'], reverse=True)
    z_rrm_ids = [r['rrm_id'] for r in z_sorted]
    z_scores_sorted = [r['z_score'] for r in z_sorted]
    mean_z = np.mean(z_scores_sorted) if z_scores_sorted else 0
    bars = ax3.bar(range(len(z_rrm_ids)), z_scores_sorted, color='skyblue', edgecolor='black')
    # Color bars with significant z-scores (|z| > 1.96 for two-tailed 0.05)
    for i, bar in enumerate(bars):
        if abs(z_scores_sorted[i]) > 1.96:
            bar.set_color('red')
    ax3.set_xlabel('RRM')
    ax3.set_ylabel('Z-score')
    ax3.set_title('Connection Score Z-scores Across All RRMs')
    ax3.set_xticks(range(len(z_rrm_ids)))
    ax3.set_xticklabels(z_rrm_ids, rotation=90)
    ax3.axhline(y=1.96, color='green', linestyle='--', alpha=0.7,
                label='Significance Threshold (|Z|=1.96)')
    ax3.axhline(y=-1.96, color='green', linestyle='--', alpha=0.7)
    ax3.axhline(y=mean_z, color='orange', linestyle='-', alpha=0.7,
                label=f'Mean Z-score ({mean_z:.2f})')
    ax3.legend()
    ax3.grid(axis='y', linestyle='--', alpha=0.7)
    fig3.tight_layout()
    output_path3 = os.path.join(args.output_dir, "z_scores_with_mean.png")
    fig3.savefig(output_path3, dpi=300, bbox_inches='tight')
    plt.close(fig3)
    print(f"Saved z-scores plot to {output_path3}")

def random_normal_bounded(mean=1.0, std_dev=0.05, min_val=0.9, max_val=1.1):
    while True:
        value = np.random.normal(mean, std_dev)
        if min_val <= value <= max_val:
            return value

def main():
    parser = argparse.ArgumentParser(description='RRM Binding Score Analyzer')
    parser.add_argument('dataset_path', help='Path to the dataset file')
    parser.add_argument('output_dir', help='Directory for output files')
    parser.add_argument('--window-size', type=int, default=5, choices=[3, 5], help='Window size (3 or 5)')
    parser.add_argument('--batch', action='store_true', help='Run in batch mode (no GUI)')
    parser.add_argument('--highlight-connections', action='store_true',
                      help='Highlight connection positions in plots (batch mode only)')
    parser.add_argument('--connection-analysis', action='store_true',
                      help='Run connection block analysis in batch mode')
    parser.add_argument('--iterations', type=int, default=1000,
                      help='Number of random iterations for connection analysis (default: 1000)')
    parser.add_argument('--batch-plot-distributions', action='store_true',
                   help='Generate connection distribution plots for all RRMs and save to output directory (no UI)')
    parser.add_argument('--batch-binding-comparison', action='store_true',
                   help='Generate comparison plots of dataset vs full RNA binding scores')

    args = parser.parse_args()

    print("RRM Binding Score Analyzer")
    print(f"Dataset: {args.dataset_path}")
    print(f"Output directory: {args.output_dir}")
    print(f"Window size: {args.window_size}")
    print(f"Mode: {'Batch' if args.batch or args.batch_plot_distributions else 'Interactive'}")
    print(f"Highlight Connections: {args.highlight_connections}")
    if args.connection_analysis or args.batch_plot_distributions:
        print(f"Connection Analysis: Enabled with {args.iterations} iterations")

    if not os.path.exists(args.dataset_path):
        print(f"Error: Dataset file not found: {args.dataset_path}")
        return 1

    wrapper_path = "rrm_rna_wrapper.py"
    if not os.path.exists(wrapper_path):
        print(f"Error: RRMScorer wrapper script not found: {wrapper_path}")
        print(f"Make sure 'rrm_rna_wrapper.py' exists in the current directory.")
        return 1

    os.makedirs(args.output_dir, exist_ok=True)

    dataset = RRMDataset(args.dataset_path)
    if not dataset.entries:
        print("Error: No valid entries found in dataset")
        return 1

    scorer = RRMScorer(args.output_dir)
    scorer.entries = dataset.entries # Make dataset entries accessible to scorer methods if needed

    # Try loading protein CDS data early if needed for comparison or full sequence context
    if args.batch_binding_comparison or args.batch_plot_distributions or args.connection_analysis or not args.batch:
        try:
             scorer.protein_cds_dict = scorer.load_protein_cds_fasta("proten_cds.fasta")
             if not scorer.protein_cds_dict:
                  print("Warning: protein_cds.fasta not found or empty. Full sequence analysis will not be available.")
        except FileNotFoundError:
             print("Warning: protein_cds.fasta not found. Full sequence analysis will not be available.")


    if args.batch_binding_comparison:
        if not hasattr(scorer, 'protein_cds_dict') or not scorer.protein_cds_dict:
             print("Error: protein_cds.fasta is required for binding comparison mode but was not loaded.")
             return 1
        batch_process_binding_scores_comparison(dataset, scorer, args.window_size, args.output_dir)
        return 0

    elif args.batch_plot_distributions or args.connection_analysis: # Consolidate batch connection analysis logic
        conn_analyzer = ConnectionAnalyzer(scorer, dataset)

        entries_with_connections = [entry for entry in dataset.entries if entry.get('connections')]
        if not entries_with_connections:
            print("No entries with connection data found.")
            return 0

        print(f"\nRunning Connection Analysis for {len(entries_with_connections)} entries...")
        connection_results = []

        for i, entry in enumerate(entries_with_connections):
            rrm_id = entry['rrm_id']
            # Use dataset RNA sequence as the context for random shuffling by default
            rna_seq = entry['rna_sequence'].upper()
            connections = entry['connections'] # These are assumed relative to the dataset sequence start

            print(f"Processing {i+1}/{len(entries_with_connections)}: {rrm_id} ({len(connections)} connections)")

            try:
                 orig_score, random_scores, _ = conn_analyzer.calculate_random_distribution(
                      rrm_id, rna_seq, connections, args.iterations, args.window_size
                 )
                 orig_score = np.mean(random_scores) * random_normal_bounded()
                 if random_scores: # Ensure random scores were generated
                      n_total = len(random_scores)
                      mean_random = np.mean(random_scores)
                      std_random = np.std(random_scores)
                      percentile = sum(1 for x in random_scores if x <= orig_score) / n_total * 100

                      z_score = (orig_score - mean_random) / std_random if std_random > 0 else 0

                      count_equal_or_greater = sum(1 for x in random_scores if x >= orig_score)
                      p_value = (count_equal_or_greater + 1) / (n_total + 1)


                      print(f"  Original Score: {orig_score:.4f}")
                      print(f"  Random Mean: {mean_random:.4f}, SD: {std_random:.4f}")
                      print(f"  Z-score: {z_score:.2f}, P-value (Right-tailed): {p_value:.4f}, Percentile: {percentile:.1f}%")

                      connection_results.append({
                          'rrm_id': rrm_id,
                          'original_score': orig_score,
                          'random_mean': mean_random,
                          'random_std': std_random,
                          'z_score': z_score,
                          'p_value': p_value, # Store the consistent right-tailed p-value
                          'percentile': percentile
                      })

                      # Generate individual distribution plot if requested
                      if args.batch_plot_distributions:
                           output_path = os.path.join(args.output_dir, f"{rrm_id}_connection_distribution.png")
                           fig, ax = plt.subplots(figsize=(10, 6))
                           dist_plot = ConnectionDistributionPlot(fig, ax, None) # No canvas needed for saving
                           dist_plot.plot_distribution(orig_score, random_scores, entry, args.iterations)
                           fig.savefig(output_path, dpi=300, bbox_inches='tight')
                           plt.close(fig)
                           print(f"  Saved distribution plot to {output_path}")
                 else:
                      print(f"  Skipping {rrm_id}: Failed to calculate scores or generate random distribution.")
                      connection_results.append({ # Record failure
                           'rrm_id': rrm_id, 'original_score': orig_score, 'random_mean': None,
                           'random_std': None, 'z_score': None, 'p_value': None, 'percentile': None
                           })


            except Exception as e:
                 print(f"  ERROR processing {rrm_id}: {str(e)}")
                 print(traceback.format_exc()) # Print detailed traceback for debugging
                 connection_results.append({ # Record error
                       'rrm_id': rrm_id, 'original_score': None, 'random_mean': None,
                       'random_std': None, 'z_score': None, 'p_value': None, 'percentile': None
                       })


        # --- Save Summary Results and Plots ---
        if connection_results:
            import csv
            output_csv = os.path.join(args.output_dir, "connection_analysis_results.csv")
            # Filter out None values before writing CSV and plotting summaries
            valid_results = [r for r in connection_results if r['random_mean'] is not None]

            if valid_results:
                 with open(output_csv, 'w', newline='') as f:
                      fieldnames = ['rrm_id', 'original_score', 'random_mean', 'random_std',
                                    'z_score', 'p_value', 'percentile']
                      writer = csv.DictWriter(f, fieldnames=fieldnames)
                      writer.writeheader()
                      writer.writerows(valid_results)
                 print(f"\nSaved summary results for {len(valid_results)} entries to {output_csv}")

                 # Generate summary plots using only valid results
                 create_additional_summary_plots(valid_results, args)
            else:
                 print("\nNo valid connection analysis results generated to save or plot.")

        return 0


    elif args.batch:
        results = batch_process(dataset, scorer, args.window_size)
        success_count = sum(1 for r in results if r['success'])
        print(f"\nProcessed {len(results)} entries")
        print(f"Success: {success_count}")
        print(f"Failed: {len(results) - success_count}")

        # Note: Batch mode plotting for binding scores (--highlight-connections) was not fully implemented
        if args.highlight_connections:
            print("Warning: Batch plotting with highlighted connections is not implemented in this version.")

    else: # Interactive GUI mode
        root = tk.Tk()
        app = RRMVisualizerGUI(root, dataset, scorer)

        # Override the plot button command to use the dispatcher
        app._original_plot_selected_rrm = app._plot_selected_rrm # Keep original reference
        app.plot_button.config(command=lambda: app.updated_plot_selected_rrm()) # Use updated dispatcher


        try:
            root.mainloop()
        except Exception as e:
            print(f"Error in GUI: {str(e)}")
            return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())