def merge_tcr_timepoints(patient_id, input_dir, output_file):
    """Merge TCR clonotype abundance data across multiple timepoints for a single patient."""
    
    timepoints = ['I', 'II', 'III', 'IV']
    all_data = []
    
    for tp in timepoints:
        file_path = os.path.join(input_dir, patient_id, f"{patient_id}-T-{tp}", "tcr_clonotype_abundance.csv")
        
        if not os.path.exists(file_path):
            print(f"Warning: Timepoint {tp} file does not exist, skipping: {file_path}")
            continue
            
        try:
            df = pd.read_csv(file_path)
            df['Timepoint'] = tp
            
            required_cols = ['Clonotype', 'Count', 'Abundance (%)']
            if not all(col in df.columns for col in required_cols):
                missing = [col for col in required_cols if col not in df.columns]
                print(f"File {file_path} missing required columns: {missing}, skipping")
                continue
                
            all_data.append(df)
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue
    
    if not all_data:
        raise ValueError("No valid data available for merging")
    
    merged = pd.concat(all_data, ignore_index=True)
    
    # Create pivot tables for counts and abundances
    pivot_count = merged.pivot(
        index='Clonotype',
        columns='Timepoint',
        values='Count'
    ).add_prefix('Count_')
    
    pivot_abundance = merged.pivot(
        index='Clonotype',
        columns='Timepoint',
        values='Abundance (%)'
    ).add_prefix('Abundance_')
    
    # Combine the pivot tables
    result = pd.concat([pivot_count, pivot_abundance], axis=1)
    
    # Reset index to make Clonotype a column
    result = result.reset_index()
    
    # Order columns by timepoint
    timepoint_order = ['I', 'II', 'III', 'IV']
    ordered_columns = ['Clonotype']
    
    for tp in timepoint_order:
        ordered_columns.extend([
            f'Count_{tp}',
            f'Abundance_{tp}'
        ])
    
    # Reorder columns (some timepoints may be missing)
    result = result[[col for col in ordered_columns if col in result.columns]]

    # Save the merged result
    result.to_csv(output_file, index=False)
    print(f"Merging completed, results saved to: {output_file}")
    print(f"Total number of clonotypes: {len(result)}")
    
    return result
