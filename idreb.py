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
        f = {}  # 正向序列字典
        r = {}  # 反向序列字典
        fq1_dict = {}  # 用于存储fq1中的原始ID和其对应的barcode, umi

        # 解析fq1文件 (正向序列)
        with open_file(fq1_file) as fq1:
            for line_num, line in enumerate(fq1):
                line = line.strip()
                if line_num % 4 == 0:  # ID行
                    seq_id = line  # 保留完整的原始ID
                elif line_num % 4 == 1:  # 序列行
                    seq = line
                elif line_num % 4 == 2:  # '+'行
                    continue
                elif line_num % 4 == 3:  # 质量分数行
                    qual = line

                    # 提取barcode和UMI
                    barcode = seq[:barcode_len]
                    if umi_len is not None:
                        umi = seq[barcode_len:barcode_len + umi_len]
                        # 新的ID格式: 保留原始ID结构，添加链信息(1:)和barcode:umi
                        if ':0:' in seq_id:  # 匹配类似 2:N:0: 的结构
                            parts = seq_id.rsplit(':0:', 1)
                            new_id = f"{parts[0]} 1:{barcode}:{umi}"
                        elif ':' in seq_id:
                            parts = seq_id.rsplit(':', 1)
                            new_id = f"{parts[0]} 1:{barcode}:{umi}"
                        else:
                            new_id = f"{seq_id} 1:{barcode}:{umi}"
                    else:
                        umi = None
                        # 新的ID格式: 保留原始ID结构，添加链信息(1:)和barcode
                        if ':0:' in seq_id:
                            parts = seq_id.rsplit(':0:', 1)
                            new_id = f"{parts[0]} 1:{barcode}"
                        elif ':' in seq_id:
                            parts = seq_id.rsplit(':', 1)
                            new_id = f"{parts[0]} 1:{barcode}"
                        else:
                            new_id = f"{seq_id} 1:{barcode}"

                    # 存储到fq1_dict，供fq2使用
                    fq1_dict[seq_id] = (barcode, umi)

                    # 存储到正向序列字典
                    if new_id not in f:
                        f[new_id] = []
                    # 移除barcode和umi后的序列和quality
                    remaining_seq = seq[barcode_len + (umi_len if umi_len else 0):]
                    remaining_qual = qual[barcode_len + (umi_len if umi_len else 0):]
                    f[new_id].append((remaining_seq, remaining_qual))

        # 解析fq2文件 (反向序列)
        with open_file(fq2_file) as fq2:
            for line_num, line in enumerate(fq2):
                line = line.strip()
                if line_num % 4 == 0:  # ID行
                    seq_id = line  # 保留完整的原始ID
                elif line_num % 4 == 1:  # 序列行
                    seq = line
                elif line_num % 4 == 2:  # '+'行
                    continue
                elif line_num % 4 == 3:  # 质量分数行
                    qual = line

                    # 从fq1_dict获取barcode和UMI
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

                    # 存储到反向序列字典
                    if new_id not in r:
                        r[new_id] = []
                    # 反向序列不需要移除barcode和umi
                    r[new_id].append((seq, qual))

        return f, r
    else:
        f = {}
        with open_file(fq1_file) as fq1:
            for line_num, line in enumerate(fq1):
                line = line.strip()
                if line_num % 4 == 0:
                    seq_id = line  # 保留完整的原始ID
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
                    # 移除barcode和umi后的序列和quality
                    remaining_seq = seq[barcode_len + (umi_len if umi_len else 0):]
                    remaining_qual = qual[barcode_len + (umi_len if umi_len else 0):]
                    f[new_id].append((remaining_seq, remaining_qual))
        return f

def cluster_by_barcode_umi(f, r=None, umi_len=None):
    clustered = defaultdict(list)
    if r:
        for new_id in f:
            # 提取barcode和UMI部分
            if umi_len is not None:
                # 从类似 @...:barcode:umi 中提取 barcode 和 umi
                parts = new_id.split(':')
                barcode = parts[-2]
                umi = parts[-1]
                barcode_umi = f"{barcode}_{umi}"
            else:
                # 从类似 @...:barcode 中提取 barcode
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

# def main():
#     input_dir = "/data/zhanqh/soloseq_test"
#     fq1_file = os.path.join(input_dir, "soloseq_nn_f1.fastq")
#     fq2_file = os.path.join(input_dir, "soloseq_nn_r2.fastq")

#     output_fq1 = os.path.join(input_dir, "zclustered1.fastq")
#     output_fq2 = os.path.join(input_dir, "zclustered2.fastq")

#     barcode_len = 16
#     #umi_len = 12  # 或设置为实际的UMI长度，例如 12

#     print("Parsing FASTQ files...")
#     if fq2_file:
#         f, r = parse_fq(fq1_file, barcode_len, umi_len, fq2_file)
#     else:
#         f = parse_fq(fq1_file, barcode_len, umi_len)

#     print("Clustering by barcode and UMI...")
#     clustered = cluster_by_barcode_umi(f, r, umi_len)

#     print("Writing clustered FASTQ files...")
#     write_clustered_fq(clustered, output_fq1, output_fq2)

#     print("Processing completed. Output files saved to:")
#     print(f"- Forward reads: {output_fq1}")
#     if fq2_file:
#         print(f"- Reverse reads: {output_fq2}")

# if __name__ == "__main__":
#     main()