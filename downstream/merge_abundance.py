import pandas as pd
import os

def merge_tcr_timepoints(patient_id, input_dir, output_file):
    """
    合并同一患者不同时间点的TCR克隆型数据，每个克隆型一行，各时间点的Count和Abundance作为列
    
    参数:
        patient_id: 患者ID (如 'P01')
        input_dir: 基础目录路径 (如 '/data/zhanqh/sampleTCR')
        output_file: 输出文件路径
    """
    # 定义四个时间点
    timepoints = ['I', 'II', 'III', 'IV']
    
    # 存储所有数据框
    all_data = []
    
    for tp in timepoints:
        file_path = os.path.join(input_dir, patient_id, f"{patient_id}-T-{tp}", "tcr_clonotype_abundance.csv")
        
        if not os.path.exists(file_path):
            print(f"警告: 时间点 {tp} 文件不存在，跳过: {file_path}")
            continue
            
        try:
            # 读取数据（严格匹配你的列名）
            df = pd.read_csv(file_path)
            
            # 添加时间点标识
            df['Timepoint'] = tp
            
            # 验证必要列是否存在
            required_cols = ['Clonotype', 'Count', 'Abundance (%)']
            if not all(col in df.columns for col in required_cols):
                missing = [col for col in required_cols if col not in df.columns]
                print(f"文件 {file_path} 缺少必要列: {missing}，跳过")
                continue
                
            all_data.append(df)
            
        except Exception as e:
            print(f"处理 {file_path} 时出错: {e}")
            continue
    
    if not all_data:
        raise ValueError("没有有效数据可供合并")
    
    # 合并所有数据
    merged = pd.concat(all_data, ignore_index=True)
    
    # 透视数据（确保使用你的准确列名）
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
    
    # 合并结果
    result = pd.concat([pivot_count, pivot_abundance], axis=1)
    
    # 重置索引并重排列
    result = result.reset_index()
    
    # 按时间点顺序排列列
    timepoint_order = ['I', 'II', 'III', 'IV']
    ordered_columns = ['Clonotype']
    
    for tp in timepoint_order:
        ordered_columns.extend([
            f'Count_{tp}',
            f'Abundance_{tp}'
        ])
    
    # 只保留实际存在的列
    result = result[[col for col in ordered_columns if col in result.columns]]
    
    # 保存结果
    result.to_csv(output_file, index=False)
    print(f"合并完成，结果保存到: {output_file}")
    print(f"总克隆型数量: {len(result)}")
    
    return result

