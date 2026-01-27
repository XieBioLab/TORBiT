import gzip
import os
from collections import defaultdict

def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def parse_fq(fq1_file, barcode_len, umi_len=None, fq2_file=None):
    if fq2_file:
        f = {}  
        r = {}  
        fq1_dict = {}  

        with open_file(fq1_file) as fq1:
            for line_num, line in enumerate(fq1):
                line = line.strip()
                if line_num % 4 == 0:  # ID行
                    seq_id = line  
                elif line_num % 4 == 1: 
                    seq = line
                elif line_num % 4 == 2: 
                    continue
                elif line_num % 4 == 3: 
                    qual = line

                    barcode = seq[:barcode_len]
                    if umi_len is not None:
                        umi = seq[barcode_len:barcode_len + umi_len]
                        if ':0:' in seq_id:  
                            parts = seq_id.rsplit(':0:', 1)
                            new_id = f"{parts[0]} 1:{barcode}:{umi}"
                        elif ':' in seq_id:
                            parts = seq_id.rsplit(':', 1)
                            new_id = f"{parts[0]} 1:{barcode}:{umi}"
                        else:
                            new_id = f"{seq_id} 1:{barcode}:{umi}"
                    else:
                        umi = None
                        if ':0:' in seq_id:
                            parts = seq_id.rsplit(':0:', 1)
                            new_id = f"{parts[0]} 1:{barcode}"
                        elif ':' in seq_id:
                            parts = seq_id.rsplit(':', 1)
                            new_id = f"{parts[0]} 1:{barcode}"
                        else:
                            new_id = f"{seq_id} 1:{barcode}"

                    fq1_dict[seq_id] = (barcode, umi)

                    if new_id not in f:
                        f[new_id] = []
                    remaining_seq = seq[barcode_len + (umi_len if umi_len else 0):]
                    remaining_qual = qual[barcode_len + (umi_len if umi_len else 0):]
                    f[new_id].append((remaining_seq, remaining_qual))

        with open_file(fq2_file) as fq2:
            for line_num, line in enumerate(fq2):
                line = line.strip()
                if line_num % 4 == 0:  
                    seq_id = line  
                elif line_num % 4 == 1:  
                    seq = line
                elif line_num % 4 == 2:  
                    continue
                elif line_num % 4 == 3:  
                    qual = line

                    if seq_id in fq1_dict:
                        barcode, umi = fq1_dict[seq_id]
                        if umi is not None:
                            if ':0:' in seq_id:
                                parts = seq_id.rsplit(':0:', 1)
                                new_id = f"{parts[0]} 2:{barcode}:{umi}"
                            elif ':' in seq_id:
                                parts = seq_id.rsplit(':', 1)
                                new_id = f"{parts[0]} 2:{barcode}:{umi}"
                            else:
                                new_id = f"{seq_id} 2:{barcode}:{umi}"
                        else:
                            if ':0:' in seq_id:
                                parts = seq_id.rsplit(':0:', 1)
                                new_id = f"{parts[0]} 2:{barcode}"
                            elif ':' in seq_id:
                                parts = seq_id.rsplit(':', 1)
                                new_id = f"{parts[0]} 2:{barcode}"
                            else:
                                new_id = f"{seq_id} 2:{barcode}"
                    else:
                        print(f"Warning: {seq_id} not found in fq1, skipping.")
                        continue

                    if new_id not in r:
                        r[new_id] = []
                    r[new_id].append((seq, qual))

        return f, r
    else:
        f = {}
        with open_file(fq1_file) as fq1:
            for line_num, line in enumerate(fq1):
                line = line.strip()
                if line_num % 4 == 0:
                    seq_id = line  
                elif line_num % 4 == 1:
                    seq = line
                elif line_num % 4 == 2:
                    continue
                elif line_num % 4 == 3:
                    qual = line

                    barcode = seq[:barcode_len]
                    if umi_len is not None:
                        umi = seq[barcode_len:barcode_len + umi_len]
                        if ':0:' in seq_id:
                            parts = seq_id.rsplit(':0:', 1)
                            new_id = f"{parts[0]} 1:{barcode}:{umi}"
                        elif ':' in seq_id:
                            parts = seq_id.rsplit(':', 1)
                            new_id = f"{parts[0]} 1:{barcode}:{umi}"
                        else:
                            new_id = f"{seq_id} 1:{barcode}:{umi}"
                    else:
                        if ':0:' in seq_id:
                            parts = seq_id.rsplit(':0:', 1)
                            new_id = f"{parts[0]} 1:{barcode}"
                        elif ':' in seq_id:
                            parts = seq_id.rsplit(':', 1)
                            new_id = f"{parts[0]} 1:{barcode}"
                        else:
                            new_id = f"{seq_id} 1:{barcode}"

                    if new_id not in f:
                        f[new_id] = []
                    remaining_seq = seq[barcode_len + (umi_len if umi_len else 0):]
                    remaining_qual = qual[barcode_len + (umi_len if umi_len else 0):]
                    f[new_id].append((remaining_seq, remaining_qual))
        return f

def cluster_by_barcode_umi(f, r=None, umi_len=None):
    clustered = defaultdict(list)
    if r:
        for new_id in f:
            if umi_len is not None:
                parts = new_id.split(':')
                barcode = parts[-2]
                umi = parts[-1]
                barcode_umi = f"{barcode}_{umi}"
            else:
                barcode_umi = new_id.split(':')[-1]
            clustered[barcode_umi].append((new_id, 'f', f[new_id]))

        for new_id in r:
            if umi_len is not None:
                parts = new_id.split(':')
                barcode = parts[-2]
                umi = parts[-1]
                barcode_umi = f"{barcode}_{umi}"
            else:
                barcode_umi = new_id.split(':')[-1]
            clustered[barcode_umi].append((new_id, 'r', r[new_id]))
        return clustered
    else:
        for new_id in f:
            if umi_len is not None:
                parts = new_id.split(':')
                barcode = parts[-2]
                umi = parts[-1]
                barcode_umi = f"{barcode}_{umi}"
            else:
                barcode_umi = new_id.split(':')[-1]
            clustered[barcode_umi].append((new_id, 'f', f[new_id]))
        return clustered

def write_clustered_fq(clustered, output_fq1, output_fq2=None):
    if output_fq2:
        os.makedirs(os.path.dirname(output_fq1), exist_ok=True)
        os.makedirs(os.path.dirname(output_fq2), exist_ok=True)

        with open(output_fq1, 'w') as fq1_out, open(output_fq2, 'w') as fq2_out:
            for barcode_umi, sequences in clustered.items():
                for seq_info in sequences:
                    new_id, direction, seq_qual_list = seq_info
                    for seq, qual in seq_qual_list:
                        if direction == 'f':
                            fq1_out.write(f"{new_id}\n{seq}\n+\n{qual}\n")
                        else:
                            fq2_out.write(f"{new_id}\n{seq}\n+\n{qual}\n")
    else:
        os.makedirs(os.path.dirname(output_fq1), exist_ok=True)
        with open(output_fq1, 'w') as fq1_out:
            for barcode_umi, sequences in clustered.items():
                for seq_info in sequences:
                    new_id, direction, seq_qual_list = seq_info
                    for seq, qual in seq_qual_list:
                        fq1_out.write(f"{new_id}\n{seq}\n+\n{qual}\n")

