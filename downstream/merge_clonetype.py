import os
import csv

input_folders = [
    './P03/P03-T-I',
    './P04/P04-T-I',
    './P05/P05-T-I',
    './P06/P06-T-I',
]
output_file = './T/fil_clone_celltype-I.tsv'

target_files = []
for folder in input_folders:
    # Construct the complete file path (case-insensitive)
    target_path = os.path.join(folder, "fil_clone_celltype.tsv")
    
    # Verify if the file exists
    if os.path.isfile(target_path):
        target_files.append(target_path)
    else:
        print(f"Warning: Path {target_path} does not exist, skipping")

# Merge files
with open(output_file, 'w', newline='', encoding='utf-8') as outfile:
    writer = csv.writer(outfile)
    header_written = False
    
    for file_path in target_files:
        try:
            with open(file_path, 'r', encoding='utf-8') as infile:
                reader = csv.reader(infile)
                header = next(reader, None)
                
                # Write header (only once)
                if not header_written and header:
                    writer.writerow(header)
                    header_written = True
                
                # Write data rows
                for row in reader:
                    writer.writerow(row)
            print(f"Merged file: {file_path}")
        except Exception as e:
            print(f"Error occurred while processing file {file_path}: {str(e)}")

print(f"\nMerge completed! Successfully merged {len(target_files)} files to {output_file}")
