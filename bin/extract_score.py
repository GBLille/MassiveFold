#!/usr/bin/env python3

import os
import pickle
import csv
import argparse

def extract_scores(pkl_file):
    """
    Extract scores from a .pickle file.

    Parameters:
        pkl_file (str): Path to the .pickle file.

    Returns:
        dict: Extracted scores.
    """
    with open(pkl_file, 'rb') as pkl:
        data = pickle.load(pkl)
    keys = list(data.keys())
    
    if 'iptm' in keys: 
        return {'ranking_ptm': data['ptm'], 'ranking_iptm': data['iptm'], 'ranking_debug': data['ranking_confidence']}
    else:
        return {'ranking_ptm': data['ptm'], 'ranking_debug': data['ranking_confidence']}

def extract_and_write_to_csv( id_batch, sequence_name, batch_start, batch_end, batch_model, dir):
    """
    Extract scores from pickle files in a directory and write them to a CSV with additional metadata.

    Parameters:
        directory (str): Path to the directory containing .pickle files.
        output_csv (str): Path to the output CSV file.
        id_batch (str): The batch ID.
        sequence_name (str): The sequence name.
        batch_start (int): The starting index of the batch.
        batch_end (int): The ending index of the batch.
        batch_model (str): The model version.
    """
    # Initialize the CSV file and write the headers
    with open(id_batch+"output.csv", mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['File', 'ID_Batch', 'Sequence_Name', 'Batch_Start', 'Batch_End', 'Batch_Model', 'Ranking_PTM', 'Ranking_IPTM', 'Ranking_Debug'])
        
        # Iterate over .pickle files and extract scores
        for filename in os.listdir(dir):
            if filename.endswith(".pickle"):
                file_path = os.path.join(dir, filename)
                scores = extract_scores(file_path)
                
                # Extract the float values
                ranking_ptm = float(scores['ranking_ptm'])
                ranking_iptm = float(scores['ranking_iptm']) if 'ranking_iptm' in scores else None
                ranking_debug = float(scores['ranking_debug'])
                
                # Write the row to CSV
                writer.writerow([file_path, id_batch, sequence_name, batch_start, batch_end, batch_model, ranking_ptm, ranking_iptm, ranking_debug])
                print(f"Processed {file_path}: PTM={ranking_ptm}, IPTM={ranking_iptm}, Debug={ranking_debug}")


def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Extract scores from .pickle files and save them to a CSV.')
    
    # Add batch and sequence parameters
    parser.add_argument('id_batch', type=str, help='Batch ID (e.g., "batch_001")')
    parser.add_argument('sequence_name', type=str, help='Sequence name (e.g., "sequence_xyz")')
    parser.add_argument('batch_start', type=int, help='Start index of the batch')
    parser.add_argument('batch_end', type=int, help='End index of the batch')
    parser.add_argument('batch_model', type=str, help='Model version (e.g., "alphafold2_multimer_v1_model_5")')
    parser.add_argument('dir', type=str, help='Directory (if not ".")')

    # Parse the arguments
    args = parser.parse_args()

    # Call the function to extract scores and write them to the CSV
    extract_and_write_to_csv( args.id_batch, args.sequence_name, args.batch_start, args.batch_end, args.batch_model, args.dir)

if __name__ == "__main__":
    main()
