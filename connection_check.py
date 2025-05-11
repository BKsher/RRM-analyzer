import re
import sys

warnings_list = []

def log_warning(message):
    """Prints a warning to stderr and adds it to the global list."""
    print(f"Warning: {message}", file=sys.stderr)
    warnings_list.append(f"Warning: {message}")

def log_error(message):
    """Prints an error to stderr, adds it to the global list, and potentially raises exception."""
    print(f"Error: {message}", file=sys.stderr)
    warnings_list.append(f"Error: {message}")

def parse_journal_file(filename="journal.pcbi.1010859.s016.fasta"):
    """
    Parses the journal file to extract unique RRM IDs and their associated RNA/DNA subsequences.
    Converts 'U' to 'T', removes '-', and stores the cleaned subsequence.
    """
    unique_rrms = {}
    print(f"Parsing journal file: {filename}...")
    line_num = 0
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()

        if len(lines) % 2 != 0:
             log_warning(f"File '{filename}' has an odd number of lines ({len(lines)}). Assuming FASTA-like pairs, last line might be ignored if it's a header.")

        i = 0
        while i < len(lines) -1 : 
            line_num = i + 1
            header_line = lines[i].strip()
            sequence_line = lines[i+1].strip()

            header = None
            sequence = None
            if header_line.startswith('>'):
                 header = header_line
                 sequence = sequence_line
                 i += 2
            elif sequence_line.startswith('>'):
                 log_warning(f"Possible swapped header/sequence at lines {line_num}-{line_num+1} in '{filename}'. Attempting to correct.")
                 header = sequence_line
                 sequence = header_line
                 i += 2
            elif i == 0 and lines[1].strip().startswith('>'):
                 log_warning(f"First line of '{filename}' doesn't look like a header, but second does. Skipping first line.")
                 i += 1
                 continue
            else:
                 log_warning(f"Cannot reliably determine FASTA pair structure at lines {line_num}-{line_num+1} in '{filename}'. Skipping entry.")
                 i += 2
                 continue

            header_content = header[1:] if header.startswith('>') else header
            match = re.match(r'([A-Z0-9]+_RRM\d+)', header_content)
            if not match:
                 log_warning(f"Could not extract RRM ID (format Pxxxx_RRMx) from header: '{header}' (Line {line_num}). Skipping entry.")
                 continue

            rrm_id = match.group(1)

            if rrm_id not in unique_rrms:
                 cleaned_sequence = sequence.upper().replace('U', 'T').replace('-', '')

                 if not cleaned_sequence:
                      log_warning(f"Sequence for {rrm_id} (Line {line_num+1}) is empty after cleaning dashes/Us. Skipping.")
                      continue
                 if not all(c in 'ATGC' for c in cleaned_sequence):
                      log_warning(f"Sequence for {rrm_id} (Line {line_num+1}) contains invalid characters after cleaning ('{sequence}' -> '{cleaned_sequence}'). Skipping.")
                      continue

                 unique_rrms[rrm_id] = cleaned_sequence

    except FileNotFoundError:
        log_error(f"File not found - {filename}")
        return None
    except Exception as e:
        log_error(f"An unexpected error occurred while parsing '{filename}' around line {line_num}: {e}")
        return None

    print(f"Found {len(unique_rrms)} unique RRM IDs with valid subsequences in {filename}.")
    return unique_rrms


def parse_extracted_data(filename="extracted_data.txt"):
    """
    Parses the extracted_data.txt file to map RRM IDs to their
    full nucleotide sequence and list of connection indices.
    """
    rrm_details = {}
    current_rrm_id = None
    current_pdb_id = None
    buffer = {}
    line_num = 0
    print(f"Parsing extracted data file: {filename}...")
    try:
        with open(filename, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()

                is_pdb_line = re.match(r'^[A-Z0-9]+_[A-Z0-9]$', line)
                end_of_block = not line or is_pdb_line

                if end_of_block and buffer:
                    rrm_id_to_store = buffer.get('rrm_id')
                    pdb_id_context = buffer.get('pdb_id', 'Unknown PDB')

                    if not rrm_id_to_store:
                         log_warning(f"Block ending near line {line_num} in '{filename}' (PDB: {pdb_id_context}) is missing 'UniProt ID:'.")
                    elif buffer.get('nucleotide_seq') is None:
                         log_warning(f"Block for RRM ID '{rrm_id_to_store}' (PDB: {pdb_id_context}) in '{filename}' is missing 'Nucleotide Sequence:'.")
                    elif buffer.get('connections') is None:
                        log_warning(f"Block for RRM ID '{rrm_id_to_store}' (PDB: {pdb_id_context}) in '{filename}' is missing 'List of connections:'.")
                    else:
                         if rrm_id_to_store not in rrm_details:
                              rrm_details[rrm_id_to_store] = {
                                   'nucleotide_seq': buffer['nucleotide_seq'],
                                   'connections': buffer['connections']
                              }

                    buffer = {}

                if is_pdb_line:
                     buffer['pdb_id'] = line

                if line.startswith("UniProt ID:"):
                    buffer['rrm_id'] = line.split(":", 1)[1].strip()
                elif line.startswith("Nucleotide Sequence:"):
                    buffer['nucleotide_seq'] = line.split(":", 1)[1].strip().upper()
                elif line.startswith("List of connections:"):
                    connection_str = line.split(":", 1)[1].strip()
                    connections_list = []
                    if connection_str:
                         try:
                             raw_parts = connection_str.split(',')
                             connections_list = []
                             for part in raw_parts:
                                 part = part.strip()
                                 if part.isdigit():
                                     connections_list.append(int(part))
                                 elif part:
                                     log_warning(f"Non-numeric value '{part}' found in connections list for RRM '{buffer.get('rrm_id', 'Unknown')}' (PDB: {buffer.get('pdb_id', 'Unknown')}) near line {line_num}. It will be ignored.")
                             buffer['connections'] = connections_list
                         except Exception as e:
                             log_warning(f"Could not parse connections string '{connection_str}' for RRM '{buffer.get('rrm_id', 'Unknown')}' (PDB: {buffer.get('pdb_id', 'Unknown')}) near line {line_num}: {e}. Assigning empty list.")
                             buffer['connections'] = []
                    else:
                         buffer['connections'] = []


            if buffer:
                 rrm_id_to_store = buffer.get('rrm_id')
                 pdb_id_context = buffer.get('pdb_id', 'Unknown PDB')
                 if not rrm_id_to_store:
                     log_warning(f"Last block in '{filename}' (PDB: {pdb_id_context}) is missing 'UniProt ID:'.")
                 elif buffer.get('nucleotide_seq') is None:
                     log_warning(f"Last block for RRM ID '{rrm_id_to_store}' (PDB: {pdb_id_context}) in '{filename}' is missing 'Nucleotide Sequence:'.")
                 elif buffer.get('connections') is None:
                      log_warning(f"Last block for RRM ID '{rrm_id_to_store}' (PDB: {pdb_id_context}) in '{filename}' is missing 'List of connections:'.")
                 else:
                     if rrm_id_to_store not in rrm_details:
                          rrm_details[rrm_id_to_store] = {
                               'nucleotide_seq': buffer['nucleotide_seq'],
                               'connections': buffer['connections']
                          }

    except FileNotFoundError:
        log_error(f"File not found - {filename}")
        return None
    except Exception as e:
        log_error(f"An unexpected error occurred while parsing '{filename}' around line {line_num}: {e}")
        return None

    print(f"Processed data for {len(rrm_details)} unique RRM IDs from {filename}.")
    return rrm_details


def check_interactions(unique_rrms, rrm_details):
    """
    Checks interactions based on cleaned subsequences and categorizes results.
    Categories: 'Interaction', 'Found_NoInteraction', 'NotFound', 'Data Missing'.
    Returns a list of tuples: (rrm_id, category, interacting_connection_or_None).
    """
    results = []
    if not unique_rrms:
        log_error("No unique RRM data provided from journal file. Cannot perform interaction check.")
        return results

    print("Checking interactions...")
    total_rrms = len(unique_rrms)
    processed_count = 0

    for rrm_id, dna_subsequence in unique_rrms.items():
        processed_count += 1
        category = "NotFound"
        interacting_connection = None
        subsequence_found_at_all = False
        interaction_found = False
        status_msg = "Processing"

        if rrm_id not in rrm_details:
            category = "Data Missing"
            status_msg = "Data Missing"
        else:
            details = rrm_details[rrm_id]
            nucleotide_seq = details.get('nucleotide_seq')
            protein_connections = details.get('connections')

            if not nucleotide_seq:
                 log_warning(f"Nucleotide sequence for {rrm_id} is missing or empty in extracted data. Cannot check interaction.")
                 category = "Data Missing"
                 status_msg = "Seq Missing"
            elif not protein_connections:
                 status_msg = "No Connections"
            else:
                 status_msg = "Checking..."
                 connection_nucleotide_indices = set()
                 for prot_idx in protein_connections:
                     if prot_idx <= 0: continue
                     start_nucleotide_idx = 3 * (prot_idx - 1)
                     if start_nucleotide_idx < len(nucleotide_seq):
                         connection_nucleotide_indices.add(start_nucleotide_idx)
                     if start_nucleotide_idx + 1 < len(nucleotide_seq):
                         connection_nucleotide_indices.add(start_nucleotide_idx + 1)
                     if start_nucleotide_idx + 2 < len(nucleotide_seq):
                         connection_nucleotide_indices.add(start_nucleotide_idx + 2)

            if nucleotide_seq:
                start_index = 0
                while True:
                    found_pos = nucleotide_seq.find(dna_subsequence, start_index)
                    if found_pos == -1:
                        break

                    subsequence_found_at_all = True

                    if connection_nucleotide_indices:
                        subsequence_indices = set(range(found_pos, found_pos + len(dna_subsequence)))
                        overlap = subsequence_indices.intersection(connection_nucleotide_indices)

                        if overlap:
                            interaction_found = True
                            status_msg = "Interaction Found"
                            first_overlap_nucleotide_idx = min(overlap)
                            interacting_connection = (first_overlap_nucleotide_idx // 3) + 1
                            break

                    start_index = found_pos + 1

            if category != "Data Missing":
                 if interaction_found:
                      category = "Interaction"
                 elif subsequence_found_at_all:
                      category = "Found_NoInteraction"
                      status_msg = "Found (No Interaction)"
                 else:
                      category = "NotFound"
                      status_msg = "Not Found"

        results.append((rrm_id, category, interacting_connection))

        progress = processed_count / total_rrms
        percent = int(progress * 100)
        bar_len = 20
        bar = '#' * int(progress * bar_len) + '-' * (bar_len - int(progress * bar_len))
        progress_line = f"\rProgress: [{bar}] {processed_count}/{total_rrms} ({percent}%) - {rrm_id} ({status_msg})"
        sys.stdout.write(progress_line.ljust(80))
        sys.stdout.flush()

    print()
    print("Interaction check complete.")
    return results

def write_results(results, filename="results.txt"):
    """Writes the interaction results into three columns and appends warnings."""
    print(f"Writing results to {filename}...")

    col1_not_found = []
    col2_interaction = []
    col3_found_no_interaction = []
    missing_data_rrms = []

    for rrm_id, category, connection in results:
        if category == "NotFound":
            col1_not_found.append(rrm_id)
        elif category == "Interaction":
            col2_interaction.append(f"{rrm_id} (Conn: {connection})")
        elif category == "Found_NoInteraction":
            col3_found_no_interaction.append(rrm_id)
        elif category == "Data Missing":
            missing_data_rrms.append(rrm_id)
            if not any(rrm_id in w for w in warnings_list):
                 warnings_list.append(f"Warning: Data for RRM ID '{rrm_id}' was missing in extracted_data.txt during interaction check.")


    max_rows = max(len(col1_not_found), len(col2_interaction), len(col3_found_no_interaction))

    try:
        with open(filename, 'w') as f:
            f.write("Subsequence Not Found\tInteraction Found\tFound (No Interaction)\n")
            f.write("---\t---\t---\n")

            for i in range(max_rows):
                c1 = col1_not_found[i] if i < len(col1_not_found) else ""
                c2 = col2_interaction[i] if i < len(col2_interaction) else ""
                c3 = col3_found_no_interaction[i] if i < len(col3_found_no_interaction) else ""
                f.write(f"{c1}\t{c2}\t{c3}\n")

            if warnings_list:
                f.write("\n\n---\n")
                f.write("Warnings and Errors Log:\n")
                f.write("---\n")
                for warning in warnings_list:
                    f.write(f"{warning}\n")
            else:
                 f.write("\n\n---\nNo warnings or errors were logged during execution.\n---\n")


        print("Results successfully written.")
    except Exception as e:
        log_error(f"Failed to write results to {filename}: {e}")
        print(f"Error writing results to {filename}. Check stderr and warnings log.", file=sys.stderr)

if __name__ == "__main__":
    print("--- RRM Interaction Check Script ---")

    unique_rrms_data = parse_journal_file()

    rrm_connection_details = parse_extracted_data()

    if unique_rrms_data:
        interaction_results = check_interactions(unique_rrms_data, rrm_connection_details)

        if interaction_results is not None:
            write_results(interaction_results)
        else:
            log_error("Interaction check failed to produce results.")
            print("\nInteraction check failed. Results file will not be written.", file=sys.stderr)
            if warnings_list:
                 try:
                     with open("results.txt.warnings_only", 'w') as f:
                           f.write("---\nWarnings and Errors Log (Interaction Check Failed):\n---\n")
                           for warning in warnings_list:
                                f.write(f"{warning}\n")
                     print("Warnings log saved to results.txt.warnings_only", file=sys.stderr)
                 except Exception as e:
                      print(f"Could not even write warnings log: {e}", file=sys.stderr)


    elif not unique_rrms_data and warnings_list:
         print("\nJournal file parsing failed or yielded no valid data. Writing warnings log.")
         try:
            with open("results.txt", 'w') as f:
                f.write("---\nJournal File Parsing Failed - Warnings and Errors Log:\n---\n")
                for warning in warnings_list:
                    f.write(f"{warning}\n")
            print("Warnings log written to results.txt")
         except Exception as e:
             print(f"Could not write warnings log to results.txt: {e}", file=sys.stderr)
    else:
         print("\nCould not proceed: No valid data obtained from journal file and no specific errors logged.")