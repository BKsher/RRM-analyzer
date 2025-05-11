import random
import os
import re
import subprocess

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