import os
import csv

# 定义四个输入文件夹路径和输出文件路径
input_folders = [
    '/data/zhanqh/new_sampleTCR/P03/P03-T-I',
    '/data/zhanqh/new_sampleTCR/P04/P04-T-I',
    '/data/zhanqh/new_sampleTCR/P05/P05-T-I',
    '/data/zhanqh/new_sampleTCR/P06/P06-T-I',
    #'/data/zhanqh/data_all/GSA/HRA005546/HRR13731/HRR1373107/tribit',
    #'/data/zhanqh/data_all/GSA/HRA005546/HRR13731/HRR1373108/tribit'
]
output_file = '/data/zhanqh/new_sampleTCR/MSI/T/fil_clone_celltype-I.tsv'

# 获取所有指定文件名路径
target_files = []
for folder in input_folders:
    # 构建目标文件完整路径（不区分大小写）
    target_path = os.path.join(folder, "fil_clone_celltype.tsv")
    
    # 验证文件是否存在
    if os.path.isfile(target_path):
        target_files.append(target_path)
    else:
        print(f"警告: 路径 {target_path} 不存在，已跳过")

# 合并文件
with open(output_file, 'w', newline='', encoding='utf-8') as outfile:
    writer = csv.writer(outfile)
    header_written = False
    
    for file_path in target_files:
        try:
            with open(file_path, 'r', encoding='utf-8') as infile:
                reader = csv.reader(infile)
                header = next(reader, None)
                
                # 写入表头（仅首次）
                if not header_written and header:
                    writer.writerow(header)
                    header_written = True
                
                # 写入数据行
                for row in reader:
                    writer.writerow(row)
            print(f"已合并文件: {file_path}")
        except Exception as e:
            print(f"处理文件 {file_path} 时发生错误: {str(e)}")

print(f"\n合并完成！共合并 {len(target_files)} 个文件到 {output_file}")