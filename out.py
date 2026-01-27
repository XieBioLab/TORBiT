from collections import defaultdict
import csv
import argparse
import re

def process_file(input_file, output_file, error_file):
    output_records = []
    error_records = []
    id_counter = defaultdict(int)
    seen_ids = set()

    with open(input_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            current_record = {}
            
            original_id = row.get('sequence_id', '')
            base_id = original_id.split('_')[0]
            
            if base_id in seen_ids:
                id_counter[base_id] += 1
                unique_id = f"{base_id}_{id_counter[base_id]}"
            else:
                seen_ids.add(base_id)
                unique_id = base_id
            
            current_record['cell_id'] = unique_id
            
            cdr3_nt = row.get('junction', '')
            cdr3_aa = row.get('junction_aa', '')
            if not cdr3_nt or not cdr3_aa:
                error_records.append(create_error_record(
                    current_record, '', '', '', "Missing CDR3 sequence"
                ))
                continue
            
            raw_v_gene = row.get('v_call', '')
            chain_type = detect_chain_type(raw_v_gene)
            
            raw_j_gene = row.get('j_call', '')
            raw_d_gene = row.get('d_call', '') if chain_type in ['VB','VD'] else ''
            
            cleaned_v = clean_v_gene(raw_v_gene, chain_type, raw_j_gene)
            
            validation_result, error_msg = validate_chain_strict(cleaned_v, raw_j_gene, chain_type)
            if not validation_result:
                error_records.append(create_error_record(
                    current_record, cleaned_v, raw_j_gene, chain_type, error_msg
                ))
                continue

            if not is_valid_allele(cleaned_v) or not is_valid_allele(raw_j_gene):
                error_records.append(create_error_record(
                    current_record, cleaned_v, raw_j_gene, chain_type, "Invalid allele format"
                ))
                continue

            current_record.update({
                'V': cleaned_v,
                'D': raw_d_gene,
                'J': raw_j_gene,
                'Chain_type': chain_type,
                'Productive': 'Yes' if row.get('productive', '').strip().upper() == 'T' else 'No',
                'CDR3_nt': cdr3_nt,
                'CDR3_aa': cdr3_aa,
                'locus': row.get('locus', '')
            })
            
            output_records.append(current_record)

    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, 
            fieldnames=['cell_id','locus','V','D','J','Chain_type','Productive','CDR3_nt','CDR3_aa'],
            delimiter='\t')
        writer.writeheader()
        writer.writerows(output_records)

    with open(error_file, 'w', newline='') as errfile:
        writer = csv.DictWriter(errfile, 
            fieldnames=['cell_id','V','J','Chain_type','error_message'],
            delimiter='\t')
        writer.writeheader()
        writer.writerows(error_records)

def clean_v_gene(raw_v: str, chain_type: str, j_gene: str) -> str:

    if '/' not in raw_v:
        return raw_v
    
    j_prefix = j_gene.split('J')[0] if 'J' in j_gene else ''
    if j_prefix:
        parts = raw_v.split('/')
        for part in parts:
            if part.startswith(j_prefix):
                return combine_allele(part, raw_v)
    
    expected_prefix = {
        'VA': 'TRAV',
        'VB': 'TRBV',
        'VG': 'TRGV',
        'VD': 'TRDV'
    }.get(chain_type, '')
    
    if expected_prefix:
        parts = raw_v.split('/')
        for part in parts:
            if part.startswith(expected_prefix):
                return combine_allele(part, raw_v)
    
    main_part = raw_v.split('/')[0]
    return combine_allele(main_part, raw_v)

def combine_allele(main_part: str, original: str) -> str:
    allele_match = re.search(r'(\*[\d\-]+)$', original)
    if allele_match:
        return main_part.split('*')[0] + allele_match.group(1)
    return main_part

def is_valid_allele(gene: str) -> bool:
    if not gene:
        return True
    if re.search(r'[^A-Za-z0-9*\-]', gene):
        return False
    if '*' in gene:
        parts = gene.split('*')
        if len(parts) != 2:
            return False
        if not parts[1] or not parts[1][0].isdigit():
            return False
    return True

def detect_chain_type(v_gene: str) -> str:
    if v_gene.startswith('TRBV'): return 'VB'
    if v_gene.startswith('TRAV'): return 'VA'
    if v_gene.startswith('TRGV'): return 'VG'
    if v_gene.startswith('TRDV'): return 'VD'
    return 'Unknown'

def validate_chain_strict(v, j, chain_type) -> (bool, str):
    if not v or not j:
        return False, "Empty gene"
    
    expected = {
        'VB': ('TRBV', 'TRBJ'),
        'VA': ('TRAV', 'TRAJ'),
        'VG': ('TRGV', 'TRGJ'),
        'VD': ('TRDV', 'TRDJ')
    }.get(chain_type, ('',''))
    
    if not expected[0] or not expected[1]:
        return False, f"Unknown chain type: {chain_type}"
    
    v_valid = v.startswith(expected[0])
    j_valid = j.startswith(expected[1])
    
    if not v_valid or not j_valid:
        error_msg = []
        if not v_valid:
            error_msg.append(f"V gene {v} doesn't match chain type {chain_type}")
        if not j_valid:
            error_msg.append(f"J gene {j} doesn't match chain type {chain_type}")
        return False, "; ".join(error_msg)
    
    return True, ""

def create_error_record(record, v, j, chain_type, msg=''):
    return {
        'cell_id': record.get('cell_id', 'N/A'),
        'V': v,
        'J': j,
        'Chain_type': chain_type,
        'error_message': msg
    }

