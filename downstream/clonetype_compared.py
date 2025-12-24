import csv
from collections import Counter
import os
from openpyxl import load_workbook, Workbook
from collections import defaultdict
def analyze_tcr_usage_to_excel(csv_file, output_excel):
    # 初始化计数器
    alpha_v_gene_count = Counter()
    alpha_j_gene_count = Counter()
    beta_v_gene_count = Counter()
    beta_j_gene_count = Counter()
    
    # 读取 CSV 文件
    with open(csv_file, 'r', encoding='utf-8-sig') as file:  # 使用utf-8-sig处理可能的BOM头
        reader = csv.DictReader(file, delimiter='\t')  # 假设是制表符分隔
        
        # 打印列名用于调试
        print("检测到的列名:", reader.fieldnames)
        
        for row in reader:
            try:
                chain_type = row["Chain_type"].strip()  # 添加strip()去除可能的空格
                v_gene = row["V"].strip()
                j_gene = row["J"].strip()
                
                # 根据链类型统计 V 和 J 基因使用次数
                if chain_type == "VA":  # α链
                    for v in v_gene.split(","):
                        alpha_v_gene_count[v.strip()] += 1
                    for j in j_gene.split(","):
                        alpha_j_gene_count[j.strip()] += 1
                elif chain_type == "VB":  # β链
                    for v in v_gene.split(","):
                        beta_v_gene_count[v.strip()] += 1
                    for j in j_gene.split(","):
                        beta_j_gene_count[j.strip()] += 1
            except KeyError as e:
                print(f"行 {reader.line_num} 缺少关键列: {e}")
                continue
    
    # 创建 Excel 文件
    wb = Workbook()
    
    # 添加工作表
    sheets_data = [
        ("αV", alpha_v_gene_count),
        ("αJ", alpha_j_gene_count),
        ("βV", beta_v_gene_count),
        ("βJ", beta_j_gene_count)
    ]
    
    for sheet_name, counter in sheets_data:
        sheet = wb.create_sheet(sheet_name)
        sheet.append([f"{sheet_name[0]} Gene", "Count"])
        for gene, count in counter.most_common():  # 按计数排序
            sheet.append([gene, count])
    
    # 删除默认工作表
    if "Sheet" in wb.sheetnames:
        wb.remove(wb["Sheet"])
    
    # 保存 Excel 文件
    wb.save(output_excel)
    print(f"统计结果已保存到 {output_excel}")





def merge_excel_files(main_dir, output_file):
    # 定义子目录列表
    sub_dirs = ["P01-T-I", "P01-T-II", "P01-T-III", "P01-T-IV"]
    # 定义要处理的sheet名称
    sheet_names = ["αV", "αJ", "βV", "βJ"]
    
    # 创建新的工作簿
    merged_wb = Workbook()
    # 删除默认创建的sheet
    while merged_wb.sheetnames:
        merged_wb.remove(merged_wb[merged_wb.sheetnames[0]])
    
    # 为每个sheet创建数据结构
    sheet_data = {name: defaultdict(dict) for name in sheet_names}
    sample_names = [d.split("-")[-1] for d in sub_dirs]  # 获取样本标识如"I", "II"等
    
    # 首先收集所有sheet的数据
    for sheet_name in sheet_names:
        for sample_idx, sub_dir in enumerate(sub_dirs):
            file_path = os.path.join(main_dir, sub_dir, "out.xlsx")
            
            if not os.path.exists(file_path):
                print(f"警告: 文件 {file_path} 不存在，跳过")
                continue
            
            try:
                wb = load_workbook(file_path)
                if sheet_name not in wb.sheetnames:
                    print(f"警告: {file_path} 中不存在sheet {sheet_name}，跳过")
                    continue
                
                ws = wb[sheet_name]
                # 假设第一列是基因名，第二列是计数
                for row in ws.iter_rows(min_row=2, values_only=True):  # 跳过标题行
                    if len(row) >= 2 and row[0] is not None:
                        gene_name = str(row[0]).strip()
                        count = row[1] if row[1] is not None else 0
                        sheet_data[sheet_name][gene_name][sample_idx] = count
            except Exception as e:
                print(f"处理 {file_path} 的 {sheet_name} 时出错: {e}")
                continue
    
    # 为每个sheet创建合并后的数据
    for sheet_name in sheet_names:
        merged_ws = merged_wb.create_sheet(sheet_name)
        
        # 写入表头
        header = ["Gene"]
        for name in sample_names:
            header.append(name)
            header.append("")  # 空白间隔列
        merged_ws.append(header[:-1])  # 去掉最后一个多余的空白列
        
        # 写入数据（按基因名排序）
        for gene in sorted(sheet_data[sheet_name].keys()):
            row_data = [gene]
            for sample_idx in range(len(sub_dirs)):
                if sample_idx in sheet_data[sheet_name][gene]:
                    row_data.append(sheet_data[sheet_name][gene][sample_idx])
                else:
                    row_data.append(0)  # 该样本中没有此基因则填0
                row_data.append("")  # 空白间隔列
            merged_ws.append(row_data[:-1])
    
    # 保存合并后的文件
    merged_wb.save(output_file)
    print(f"所有文件已按sheet合并到 {output_file}，相同基因已对齐")
