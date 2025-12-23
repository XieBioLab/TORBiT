import csv
from collections import Counter
from openpyxl import Workbook

def analyze_tcr_usage_to_excel(csv_file, output_excel):
    # 初始化计数器
    alpha_v_gene_count = Counter()
    alpha_j_gene_count = Counter()
    beta_v_gene_count = Counter()
    beta_j_gene_count = Counter()
    
    # 读取 CSV 文件
    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            chain_type = row["Chain_type"]
            v_gene = row["V"]
            j_gene = row["J"]
            
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
    
    # 创建 Excel 文件
    wb = Workbook()
    
    # 添加 αV 表
    alpha_v_sheet = wb.create_sheet("αV")
    alpha_v_sheet.append(["V Gene", "Count"])
    for v_gene, count in alpha_v_gene_count.items():
        alpha_v_sheet.append([v_gene, count])
    
    # 添加 αJ 表
    alpha_j_sheet = wb.create_sheet("αJ")
    alpha_j_sheet.append(["J Gene", "Count"])
    for j_gene, count in alpha_j_gene_count.items():
        alpha_j_sheet.append([j_gene, count])
    
    # 添加 βV 表
    beta_v_sheet = wb.create_sheet("βV")
    beta_v_sheet.append(["V Gene", "Count"])
    for v_gene, count in beta_v_gene_count.items():
        beta_v_sheet.append([v_gene, count])
    
    # 添加 βJ 表
    beta_j_sheet = wb.create_sheet("βJ")
    beta_j_sheet.append(["J Gene", "Count"])
    for j_gene, count in beta_j_gene_count.items():
        beta_j_sheet.append([j_gene, count])
    
    # 删除默认工作表
    if "Sheet" in wb.sheetnames:
        wb.remove(wb["Sheet"])
    
    # 保存 Excel 文件
    wb.save(output_excel)
    print(f"统计结果已保存到 {output_excel}")

if __name__ == "__main__":
    csv_file = "/data/zhanqh/sampleTCR/P01/P01-T-I/output.csv"  # 替换为你的 CSV 文件路径
    output_excel = "/data/zhanqh/sampleTCR/P01/P01-T-I/tcr_usage_analysis.xlsx"  # 输出 Excel 文件路径
    analyze_tcr_usage_to_excel(csv_file, output_excel)
