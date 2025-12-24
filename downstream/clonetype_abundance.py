import pandas as pd

def calculate_tcr_clonotype_abundance(input_file, output_file):
    # 尝试读取数据
    try:
        data = pd.read_csv(input_file, sep='\t')
    except Exception as e:
        print(f"无法读取输入文件: {e}")
        return
    
    # 检查数据列是否齐全
    required_columns = ['V', 'J', 'CDR3aa']  # 删除了 'Productive' 列
    for col in required_columns:
        if col not in data.columns:
            print(f"缺少必要的列: {col}")
            return

    # 定义克隆型：组合 (V, J, CDR3_aa)
    data['Clonotype'] = (
        data['V'] + "_" +
        data['J'] + "_" +
        data['CDR3aa']
    )
    
    # 打印克隆型定义后的数据（调试信息）
    print("克隆型数据示例:")
    print(data[['Clonotype']].head())
    
    # 统计克隆型丰度
    clonotype_counts = data['Clonotype'].value_counts()
    if clonotype_counts.empty:
        print("没有统计到任何克隆型数据。")
        return
    
    total_clonotypes = clonotype_counts.sum()
    
    # 转换为 DataFrame 并计算丰度百分比
    clonotype_abundance = clonotype_counts.reset_index()
    clonotype_abundance.columns = ['Clonotype', 'Count']
    clonotype_abundance['Abundance (%)'] = (clonotype_abundance['Count'] / total_clonotypes) * 100
    
    # 保存结果到文件
    clonotype_abundance.to_csv(output_file, index=False)
    print(f"克隆型丰度已保存到文件: {output_file}")

