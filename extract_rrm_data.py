import os
import re

def parse_rrm_id_file(file_path):
    """Parse a file with RRM IDs and their ranges (P26378_RRM2_1G2E_A_127-195_139-207_123-222_B)"""
    rrm_info = {}
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('_')
            if len(parts) < 6:
                print(f"Warning: Skipping malformed line: {line}")
                continue
            
            # Extract protein ID and RRM number
            prot_id = f"{parts[0]}_{parts[1]}"
            
            # Find ranges (they contain '-')
            ranges = [part for part in parts if '-' in part]
            if len(ranges) < 2:
                print(f"Warning: Could not find enough ranges in line: {line}")
                continue
            
            # Extract second range
            second_range = ranges[1]
            try:
                start, end = map(int, second_range.split('-'))
                rrm_info[prot_id] = {
                    'full_id': line,
                    'range': (start, end)
                }
            except ValueError:
                print(f"Warning: Invalid range format in line: {line}")
                continue
    
    return rrm_info

def parse_fasta_content(content):
    """Parse FASTA content with alternating protein/nucleotide entries"""
    entries = {}
    current_id = None
    current_seq = []
    
    for line in content.split('\n'):
        line = line.strip()
        if not line:
            continue
        
        if line.startswith('>'):
            # Save previous sequence if it exists
            if current_id:
                entries[current_id] = ''.join(current_seq)
            
            # Start new sequence
            current_id = line[1:].split()[0]  # Remove '>' and take first word
            current_seq = []
        else:
            current_seq.append(line)
    
    # Save the last sequence
    if current_id:
        entries[current_id] = ''.join(current_seq)
    
    # Separate protein and nucleotide sequences
    protein_entries = {}
    nucleotide_entries = {}
    
    for seq_id, sequence in entries.items():
        if seq_id.startswith('NM_') or seq_id.startswith('ENST'):
            nucleotide_entries[seq_id] = sequence
        else:
            protein_entries[seq_id] = sequence
    
    # Map protein IDs to nucleotide IDs (assuming they alternate in the file)
    protein_to_nucleotide = {}
    protein_ids = list(protein_entries.keys())
    nucleotide_ids = list(nucleotide_entries.keys())
    
    for i in range(min(len(protein_ids), len(nucleotide_ids))):
        protein_to_nucleotide[protein_ids[i]] = nucleotide_ids[i]
    
    return protein_entries, nucleotide_entries, protein_to_nucleotide

def extract_sequences(rrm_info, protein_entries, nucleotide_entries, protein_to_nucleotide):
    """Extract protein and nucleotide subsequences based on RRM info"""
    results = []
    
    for rrm_id, info in rrm_info.items():
        base_id = rrm_id.split('_')[0]  # e.g., P26378
        
        # Find protein sequence
        prot_seq = protein_entries.get(base_id)
        if not prot_seq:
            print(f"Warning: No protein sequence found for {base_id}")
            continue
        
        # Extract protein subsequence
        start, end = info['range']
        if 1 <= start <= len(prot_seq) and 1 <= end <= len(prot_seq):
            prot_subseq = prot_seq[start-1:end]
        else:
            print(f"Warning: Range {start}-{end} is out of bounds for protein {base_id} (length: {len(prot_seq)})")
            continue
        
        # Find nucleotide sequence
        nuc_id = protein_to_nucleotide.get(base_id)
        nuc_seq = nucleotide_entries.get(nuc_id) if nuc_id else None
        
        if nuc_seq:
            # Calculate nucleotide range based on protein range
            nuc_start = start * 3 + 1
            nuc_end = end * 3 + 3
            
            if 1 <= nuc_start <= len(nuc_seq) and 1 <= nuc_end <= len(nuc_seq):
                nuc_subseq = nuc_seq[nuc_start-1:nuc_end]
                nuc_range = f"{nuc_start}-{nuc_end}"
            else:
                print(f"Warning: Nucleotide range {nuc_start}-{nuc_end} is out of bounds for {nuc_id} (length: {len(nuc_seq)})")
                nuc_subseq = "Not available"
                nuc_range = "Not available"
        else:
            nuc_id = "Not available"
            nuc_subseq = "Not available"
            nuc_range = "Not available"
        
        # Create formatted entry
        entry = {
            'full_id': info['full_id'],
            'uniprot_id': base_id,
            'ccds_id': nuc_id if nuc_id and nuc_id.startswith('CCDS') else None,
            'protein_range': f"{start}-{end}",
            'nucleotide_range': nuc_range,
            'protein_sequence': prot_subseq,
            'nucleotide_sequence': nuc_subseq
        }
        
        results.append(entry)
    
    return results

def write_output(results, output_file):
    """Write results to output file in the specified format"""
    with open(output_file, 'w') as out:
        for entry in results:
            out.write(f">>{entry['full_id']}\n")
            out.write(f"UniProt ID: {entry['uniprot_id']}\n")
            
            if entry['ccds_id']:
                out.write(f"CCDS ID: {entry['ccds_id']}\n")
            
            out.write(f"Protein Range: {entry['protein_range']}\n")
            out.write(f"Nucleotide Range: {entry['nucleotide_range']}\n")
            out.write(f"Protein Sequence: {entry['protein_sequence']}\n")
            out.write(f"Nucleotide Sequence: {entry['nucleotide_sequence']}\n")
            out.write("\n")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Extract RRM sequences from FASTA files')
    parser.add_argument('--rrm-info', required=True, help='File with RRM IDs and ranges')
    parser.add_argument('--fasta', required=True, help='FASTA file with protein and nucleotide sequences')
    parser.add_argument('--output', default='extracted_sequences.txt', help='Output file')
    
    args = parser.parse_args()
    
    # Parse RRM info
    rrm_info = parse_rrm_id_file(args.rrm_info)
    print(f"Parsed {len(rrm_info)} RRM entries from {args.rrm_info}")
    
    # Read FASTA file
    with open(args.fasta, 'r') as f:
        fasta_content = f.read()
    
    # Parse FASTA content
    protein_entries, nucleotide_entries, protein_to_nucleotide = parse_fasta_content(fasta_content)
    print(f"Found {len(protein_entries)} protein sequences and {len(nucleotide_entries)} nucleotide sequences")
    
    # Extract sequences
    results = extract_sequences(rrm_info, protein_entries, nucleotide_entries, protein_to_nucleotide)
    print(f"Extracted {len(results)} sequences")
    
    # Write output
    write_output(results, args.output)
    print(f"Results written to {args.output}")

if __name__ == "__main__":
    main()
