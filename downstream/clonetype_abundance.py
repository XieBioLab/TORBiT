import pandas as pd

def calculate_tcr_clonotype_abundance(input_file, output_file):
    try:
        data = pd.read_csv(input_file, sep='\t')
    except Exception as e:
        print(f"Unable to read input file: {e}")
        return
    
    required_columns = ['V', 'J', 'CDR3aa']  
    for col in required_columns:
        if col not in data.columns:
            print(f"Missing required column: {col}")
            return

    data['Clonotype'] = (
        data['V'] + "_" +
        data['J'] + "_" +
        data['CDR3aa']
    )
    
    print("Clonotype data sample:")
    print(data[['Clonotype']].head())
    
    clonotype_counts = data['Clonotype'].value_counts()
    if clonotype_counts.empty:
        print("No clonotype counts were collected")
        return
    
    total_clonotypes = clonotype_counts.sum()
    
    clonotype_abundance = clonotype_counts.reset_index()
    clonotype_abundance.columns = ['Clonotype', 'Count']
    clonotype_abundance['Abundance (%)'] = (clonotype_abundance['Count'] / total_clonotypes) * 100
    
    clonotype_abundance.to_csv(output_file, index=False)
    print(f"Clonotype abundance has been saved to file: {output_file}")

