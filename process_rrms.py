import os
import subprocess
import itertools
import traceback # For detailed error messages

# --- Configuration ---
EXTRACTED_DATA_FILE = "extracted_data.txt"
RRM_RNA_WRAPPER_SCRIPT = "rrm_rna_wrapper.py"
OUTPUT_FOLDER = "windows"
RNA_BASES = ['A', 'C', 'G', 'U']
BLOCK_SIZE = 5
NUM_UNIQUE_BLOCKS = 4**BLOCK_SIZE  # 1024

# --- Testing Flag ---
TESTING_MODE = False # <<< CHANGE THIS TO False FOR FULL RUN

# --- Step 1: Parse extracted_data.txt ---
def parse_data_file(filepath):
    """
    Parses the extracted_data.txt file to get PDB ID with chain and UniProt ID.
    """
    rrm_info_list = []
    try:
        with open(filepath, 'r') as f:
            current_pdb_chain = None
            current_uniprot_id = None
            for line_idx, line_content in enumerate(f):
                line = line_content.strip()
                if line.startswith(">>"):
                    if current_pdb_chain and current_uniprot_id:
                        rrm_info_list.append({
                            "pdb_chain": current_pdb_chain,
                            "uniprot_id": current_uniprot_id
                        })
                    elif current_pdb_chain and not current_uniprot_id:
                        print(f"Warning in {filepath}: PDB entry '{current_pdb_chain}' (around line {line_idx+1}) seems to be missing a UniProt ID before a new '>>' entry.")

                    current_pdb_chain = line[2:].strip()
                    current_uniprot_id = None
                elif line.startswith("UniProt ID:"):
                    current_uniprot_id = line.split(":", 1)[1].strip()
            
            if current_pdb_chain and current_uniprot_id:
                rrm_info_list.append({
                    "pdb_chain": current_pdb_chain,
                    "uniprot_id": current_uniprot_id
                })
            elif current_pdb_chain and not current_uniprot_id:
                 print(f"Warning in {filepath}: Last PDB entry '{current_pdb_chain}' is missing a UniProt ID.")
    except FileNotFoundError:
        print(f"Error: File '{filepath}' not found.")
        return None
    return rrm_info_list

# --- Step 2: Generate the long RNA sequence ---
def generate_rna_sequence_and_blocks():
    """
    Generates 1024 unique 5-mer blocks and concatenates them.
    Returns the full RNA string and the list of individual blocks.
    """
    if 4**BLOCK_SIZE != NUM_UNIQUE_BLOCKS:
        print(f"Warning: Mismatch in expected unique blocks. 4^{BLOCK_SIZE} = {4**BLOCK_SIZE}, NUM_UNIQUE_BLOCKS = {NUM_UNIQUE_BLOCKS}")

    unique_blocks_list = []
    for p in itertools.product(RNA_BASES, repeat=BLOCK_SIZE):
        unique_blocks_list.append("".join(p))
        if len(unique_blocks_list) >= NUM_UNIQUE_BLOCKS:
            break
    
    if len(unique_blocks_list) < NUM_UNIQUE_BLOCKS:
        print(f"Error: Could only generate {len(unique_blocks_list)} unique blocks, expected {NUM_UNIQUE_BLOCKS}.")
        return None, None

    long_rna_sequence = "".join(unique_blocks_list)
    
    print(f"Generated {len(unique_blocks_list)} unique RNA blocks of length {BLOCK_SIZE}.")
    print(f"Total RNA sequence length: {len(long_rna_sequence)} ({NUM_UNIQUE_BLOCKS} * {BLOCK_SIZE}).")
    return long_rna_sequence, unique_blocks_list

# --- Main execution ---
def main():
    print("Starting RRM processing script...")
    if TESTING_MODE:
        print("!!! TESTING MODE ENABLED: Only the first RRM entry will be processed. !!!")
        print("!!! To process all entries, set TESTING_MODE = False in the script. !!!")

    rrm_data = parse_data_file(EXTRACTED_DATA_FILE)
    if not rrm_data:
        print("Exiting due to parsing error or no data from extracted_data.txt.")
        return

    print(f"Found {len(rrm_data)} RRM entries to process from '{EXTRACTED_DATA_FILE}'.")
    display_limit = 5
    for i, entry in enumerate(rrm_data):
        if i < display_limit:
            print(f"  {i+1}. PDB: {entry['pdb_chain']}, UniProt: {entry['uniprot_id']}")
        elif i == display_limit:
            print(f"  ... and {len(rrm_data) - display_limit} more.")
            break
    
    long_rna, rna_blocks_ordered = generate_rna_sequence_and_blocks()
    if not long_rna:
        print("Exiting due to RNA generation error.")
        return

    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)
        print(f"Created output directory: '{OUTPUT_FOLDER}'")

    processed_count = 0
    for rrm_idx, rrm_entry in enumerate(rrm_data):
        if TESTING_MODE and rrm_idx > 0:
            print("\nTesting mode: Reached limit. Stopping further RRM processing.")
            break
        
        pdb_chain_id = rrm_entry['pdb_chain']
        uniprot_id = rrm_entry['uniprot_id']
        
        print(f"\nProcessing RRM {rrm_idx + 1}/{len(rrm_data)}: {pdb_chain_id} (UniProt: {uniprot_id})...")

        command = [
            "python", RRM_RNA_WRAPPER_SCRIPT,
            "-RRM", uniprot_id,
            "-RNA", long_rna,
            "-ws", str(BLOCK_SIZE)
        ]
        
        process = None
        try:
            print(f"  Running command: {' '.join(command)}")
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1, universal_newlines=True)
            
            # It's generally better to use communicate for Popen
            stdout, stderr = process.communicate(timeout=600) # 10 minutes timeout

            if process.returncode != 0:
                print(f"  Error: '{RRM_RNA_WRAPPER_SCRIPT}' for {pdb_chain_id} (UniProt: {uniprot_id}) exited with code {process.returncode}.")
                if stderr: print(f"  Stderr from wrapper:\n{stderr.strip()}")
                if stdout: print(f"  Stdout from wrapper (if any on error):\n{stdout.strip()}")
                continue

            # --- MODIFIED SCORE PARSING ---
            raw_output_lines = stdout.strip().split('\n')
            all_scores = []
            lines_not_parsed_as_score = 0
            
            for line_num, s_line in enumerate(raw_output_lines):
                s_val_stripped = s_line.strip()
                if not s_val_stripped: # Skip empty lines
                    continue
                
                parts = s_val_stripped.split() # Split by whitespace
                if not parts: # Should not happen if s_val_stripped is not empty
                    lines_not_parsed_as_score += 1
                    if lines_not_parsed_as_score <= 5:
                        print(f"  Debug: Encountered an empty line after stripping (line #{line_num+1}): '{s_line}'. Skipping.")
                    continue

                # Assume score is the last part
                score_candidate = parts[-1]
                try:
                    score_float = float(score_candidate)
                    all_scores.append(score_float)
                except ValueError:
                    lines_not_parsed_as_score += 1
                    # Print only a few initial warnings to avoid flooding the console
                    if lines_not_parsed_as_score <= 5:
                         print(f"  Warning: Could not convert last part ('{score_candidate}') of stdout line #{line_num+1} to float. Full line: '{s_val_stripped}'. Skipping this line.")
            
            if lines_not_parsed_as_score > 0:
                print(f"  Info: Total lines from wrapper output not parsed as scores for {pdb_chain_id}: {lines_not_parsed_as_score}")
            # --- END OF MODIFIED SCORE PARSING ---

            expected_total_scores = len(long_rna) - BLOCK_SIZE + 1
            if len(all_scores) != expected_total_scores:
                print(f"  Warning: Expected {expected_total_scores} numerical scores from wrapper for {pdb_chain_id}, but extracted {len(all_scores)} scores.")
                print(f"    This could be due to the wrapper not producing one score per sliding window, or parsing issues.")

            selected_scores = []
            if all_scores: # Only proceed if we have some scores
                # Select every BLOCK_SIZE-th score (0th, 5th, 10th, etc.)
                for i in range(0, len(all_scores), BLOCK_SIZE):
                    selected_scores.append(all_scores[i])
                    if len(selected_scores) >= NUM_UNIQUE_BLOCKS:
                        selected_scores = selected_scores[:NUM_UNIQUE_BLOCKS] # Ensure we don't exceed
                        break
            
            if len(selected_scores) != NUM_UNIQUE_BLOCKS:
                print(f"  Warning: Expected {NUM_UNIQUE_BLOCKS} scores for our unique blocks for {pdb_chain_id}, but selected {len(selected_scores)}.")
                print(f"    Total scores extracted: {len(all_scores)}. Total lines not parsed: {lines_not_parsed_as_score}.")
                if not selected_scores:
                    print(f"  No valid scores selected for {pdb_chain_id}. Skipping file save.")
                    continue

            output_filename = os.path.join(OUTPUT_FOLDER, f"{pdb_chain_id}.txt")
            with open(output_filename, 'w') as outfile:
                for score in selected_scores:
                    outfile.write(f"{score}\n")
            
            print(f"  Successfully processed and saved {len(selected_scores)} scores to '{output_filename}'.")
            processed_count += 1

        except subprocess.TimeoutExpired:
            print(f"  Error: Command for {pdb_chain_id} (UniProt: {uniprot_id}) timed out after 10 minutes.")
            if process:
                process.kill()
                # process.communicate() # Ensure pipes are cleared after kill
        except FileNotFoundError:
            print(f"  Critical Error: Script '{RRM_RNA_WRAPPER_SCRIPT}' not found. Make sure it's in the correct path and Python can execute it.")
            print("  Aborting further processing.")
            break
        except Exception as e:
            print(f"  An unexpected error occurred while processing {pdb_chain_id} (UniProt: {uniprot_id}): {e}")
            print("  Traceback:")
            traceback.print_exc()

    print(f"\nScript finished. Processed {processed_count} RRM entries.")
    if TESTING_MODE and processed_count == 1 and len(rrm_data) > 1 :
         print(f"Note: Processing was limited to the first RRM entry due to TESTING_MODE. {len(rrm_data) - 1} entries were skipped.")
    elif TESTING_MODE and processed_count == 0 and len(rrm_data) > 0:
         print(f"Note: TESTING_MODE was on, but the first RRM entry was not successfully processed.")
    elif not TESTING_MODE:
        print(f"All {len(rrm_data)} RRM entries were scheduled for processing.")


if __name__ == "__main__":
    main()