import os
import logging
import shutil
import subprocess
from collections import defaultdict
from multiprocessing import Lock
import psutil
import multiprocessing
import csv
# 初始化日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# 线程安全锁
lock = Lock()

# def get_barcode_umi(r1_files):
#     """带详细调试日志的解析函数"""
#     if isinstance(r1_files, str):
#         r1_files = [r1_files]

#     barcode_umi_pairs = []
#     is_bulk = False
    
#     for r1_file in r1_files:
#         with open(r1_file, 'r') as r1:
#             for i, line in enumerate(r1):
#                 if i % 4 != 0:  # 只处理header行
#                     continue
                    
#                 header = line.strip()
#                 # 调试日志
#                 if i < 12:  # 打印前3条记录的header
#                     logging.debug(f"Header样本 {i//4}: {header[:100]}...")
                
#                 if " " not in header:
#                     is_bulk = True
#                     logging.warning(f"检测到bulk格式header: {header[:100]}...")
#                     return [], True
                
#                 try:
#                     parts = header.split(" ")
#                     barcode_umi = parts[1].split(":")
#                     barcode = barcode_umi[1] if barcode_umi else None
#                     umi = barcode_umi[2] if len(barcode_umi) > 1 else None
#                     if barcode:
#                         barcode_umi_pairs.append((barcode, umi))
#                 except Exception as e:
#                     logging.warning(f"解析失败: {header[:100]}... 错误: {str(e)}")
    
#     logging.info(f"成功解析 {len(barcode_umi_pairs)} 条barcode记录")
#     return barcode_umi_pairs, False
from collections import defaultdict
import logging

def get_barcode_umi(r1_files):
    """带详细调试日志的解析函数"""
    if isinstance(r1_files, str):
        r1_files = [r1_files]

    barcode_umi_counts = defaultdict(int)  # 使用defaultdict统计每个组合的出现次数
    is_bulk = False
    
    for r1_file in r1_files:
        with open(r1_file, 'r') as r1:
            for i, line in enumerate(r1):
                if i % 4 != 0:  # 只处理header行
                    continue
                    
                header = line.strip()
                # 调试日志
                if i < 12:  # 打印前3条记录的header
                    logging.debug(f"Header样本 {i//4}: {header[:100]}...")
                
                if " " not in header:
                    is_bulk = True
                    logging.warning(f"检测到bulk格式header: {header[:100]}...")
                    return list(barcode_umi_counts.keys()), True  # 返回当前已解析的唯一组合和True
                
                try:
                    parts = header.split(" ")
                    barcode_umi = parts[1].split(":")
                    
                    # 判断是否有umi
                    if len(barcode_umi) == 2:  # 只有一个:，只有barcode
                        barcode = barcode_umi[1]
                        umi = None
                    elif len(barcode_umi) == 3:  # 有两个:，第一个后面的为barcode，第二个后面的为umi
                        barcode = barcode_umi[1]
                        umi = barcode_umi[2]
                    else:
                        raise ValueError("不支持的header格式")
                    
                    if barcode:
                        barcode_umi = (barcode, umi) if umi is not None else (barcode, None)
                        barcode_umi_counts[barcode_umi] += 1  # 统计组合的出现次数
                except Exception as e:
                    logging.warning(f"解析失败: {header[:100]}... 错误: {str(e)}")
    
    logging.info(f"成功解析 {len(barcode_umi_counts)} 条唯一的barcode记录")
    return list(barcode_umi_counts.keys()), False
def monitor_memory():
    mem = psutil.virtual_memory()
    if mem.percent > 90:
        logging.warning(f"High memory usage detected: {mem.percent}%")

def run_trinity(subset_r1, subset_r2=None, output_dir=None, i=None):
    """改进的Trinity运行函数"""
    os.makedirs(output_dir, exist_ok=True)
    base_cmd = [
        'Trinity',
        '--seqType', 'fq',
        '--max_memory', '100G',
        '--CPU', '6',
        '--output', output_dir,
        '--no_version_check'
    ]

    if subset_r2:
        base_cmd += ['--left', subset_r1, '--right', subset_r2]
    else:
        base_cmd += ['--single', subset_r1]
    
    result = subprocess.run(base_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    return result.returncode == 0

def run_catt(f1, f2=None, output_path=None, species="hs", i=None):
    """运行 CATT (TCR 组装工具)"""
    # 确保输出目录存在
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
    catt_path = "/data/zhanqh/catt/CATT/catt"
    out_file = os.path.join(output_path, 'clones.csv')
    
    # 构建基本命令
    base_cmd = [
        'python', catt_path,
        '--chain',"TRA"
        '--f1', f1,
        '--species', species,
        '-o', out_file
    ]
    
    if f2:
        base_cmd += ['--f2', f2]
    
    logging.info(f"Executing CATT command: {' '.join(base_cmd)}")  # 记录完整命令
    
    try:
        result = subprocess.run(
            base_cmd, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            text=True,
            check=True  # 这将自动在非零返回码时引发异常
        )
        logging.info(f"CATT stdout: {result.stdout[:200]}...")  # 截断输出
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"CATT failed with return code {e.returncode}")
        logging.error(f"Error output: {e.stderr}")
        return False
    except Exception as e:
        logging.error(f"Unexpected error running CATT: {str(e)}")
        return False

def sc_assembly(subset, r1_file, output_dir, output_fa, r2_file=None):
    for barcode, umi in subset:
        # 生成唯一标识符
        output_id = f"{barcode}_{umi}" if umi is not None else barcode
        
        # 创建临时工作目录
        temp_work_dir = os.path.join(output_dir, output_id)
        os.makedirs(temp_work_dir, exist_ok=True)
        
        # 创建 Trinity 专用输出目录
        trinity_output_dir = os.path.join(temp_work_dir, "trinity")
        os.makedirs(trinity_output_dir, exist_ok=True)
        
        # 临时文件路径
        temp_r1 = os.path.join(temp_work_dir, 'temp_r1.fq')
        if r2_file:
            temp_r2 = os.path.join(temp_work_dir, 'temp_r2.fq')
        
        sequence_count = 0
        try:
            # 按照是否有 r2_file 选择不同的读取逻辑
            if r2_file is not None:
                with open(temp_r1, 'w') as out_r1, open(temp_r2, 'w') as out_r2:
                    with open(r1_file, 'r') as r1, open(r2_file, 'r') as r2:
                        while True:
                            # 读取双端 FASTQ 的 4 行记录
                            r1_header = r1.readline().strip()
                            r2_header = r2.readline().strip()
                            if not r1_header:
                                break
                            r1_sequence = r1.readline().strip()
                            r2_sequence = r2.readline().strip()
                            r1_plus = r1.readline().strip()
                            r2_plus = r2.readline().strip()
                            r1_quality = r1.readline().strip()
                            r2_quality = r2.readline().strip()

                            # 新的barcode和umi提取逻辑
                            key_pos = r1_header.find(" 1:")
                            if key_pos == -1:
                                continue  # 跳过不符合格式的记录
                            
                            barcode_umi_part = r1_header[key_pos + 3:]  # 跳过" 1:"
                            parts = barcode_umi_part.split(':')
                            
                            r1_barcode = None
                            r1_umi = None
                            
                            if len(parts) >= 1:
                                r1_barcode = parts[0]
                                if len(parts) >= 2:
                                    r1_umi = parts[1].split()[0]

                            # 检查是否匹配当前barcode和umi（如果有）
                            if umi is not None:
                                match_condition = (r1_barcode == barcode and r1_umi == umi)
                            else:
                                match_condition = (r1_barcode == barcode)

                            if match_condition:
                                out_r1.write(f"{r1_header}\n{r1_sequence}\n{r1_plus}\n{r1_quality}\n")
                                out_r2.write(f"{r2_header}\n{r2_sequence}\n{r2_plus}\n{r2_quality}\n")
                                sequence_count += 1
            else:
                with open(temp_r1, 'w') as out_r1:
                    with open(r1_file, 'r') as r1:
                        while True:
                            r1_header = r1.readline().strip()
                            if not r1_header:
                                break
                            r1_sequence = r1.readline().strip()
                            r1_plus = r1.readline().strip()
                            r1_quality = r1.readline().strip()

                            # 新的barcode和umi提取逻辑
                            key_pos = r1_header.find(" 1:")
                            if key_pos == -1:
                                continue  # 跳过不符合格式的记录
                            
                            barcode_umi_part = r1_header[key_pos + 3:]  # 跳过" 1:"
                            parts = barcode_umi_part.split(':')
                            
                            r1_barcode = None
                            r1_umi = None
                            
                            if len(parts) >= 1:
                                r1_barcode = parts[0]
                                if len(parts) >= 2:
                                    r1_umi = parts[1].split()[0]

                            # 检查是否匹配当前barcode和umi（如果有）
                            if umi is not None:
                                match_condition = (r1_barcode == barcode and r1_umi == umi)
                            else:
                                match_condition = (r1_barcode == barcode)

                            if match_condition:
                                out_r1.write(f"{r1_header}\n{r1_sequence}\n{r1_plus}\n{r1_quality}\n")
                                sequence_count += 1
        except Exception as e:
            logging.error(f"Error processing FASTQ files for {output_id}: {e}")
            continue

        # 当仅有一条序列时，直接写入最终文件
        if sequence_count == 1:
            try:
                with lock:
                    with open(output_fa, 'a') as fa_out, open(temp_r1, 'r') as temp_r1:
                        for line in temp_r1:
                            if line.startswith('@'):
                                fa_out.write(f">{output_id}\n")
                            elif not line.startswith(('+', '@', '!')):
                                fa_out.write(line)
                                break
            except Exception as e:
                logging.error(f"Error writing single sequence to output file {output_fa}: {e}")
            #尝试删除临时目录
            try:
                shutil.rmtree(temp_work_dir)
            except Exception as e:
                logging.error(f"Error cleaning up temporary directory for {output_id}: {e}")
            continue

        # 运行Trinity
        success = False
        if r2_file:
            success = run_trinity(temp_r1, temp_r2, trinity_output_dir)
        else:
            success = run_trinity(temp_r1, None, trinity_output_dir)
        
        if success:
            # 改进的Trinity输出文件搜索机制
            trinity_output = None
            # 首先检查Trinity默认输出位置
            paths = os.path.join(trinity_output_dir, 'Trinity.fasta'),  # 默认输出文件名

            # 尝试所有可能的路径
            for path in paths:
                if os.path.exists(path):
                    trinity_output = path
                    break
            
            if trinity_output:
                try:
                    with lock:
                        with open(trinity_output, 'r') as f_in, open(output_fa, 'a') as f_out:
                            for line in f_in:
                                if line.startswith('>'):
                                    f_out.write(f">{output_id}\n")
                                else:
                                    f_out.write(line)
                except Exception as e:
                    logging.error(f"Error writing output for {output_id}: {e}")
                
        # 尝试删除临时目录
        try:
            shutil.rmtree(temp_work_dir)
        except Exception as e:
            logging.error(f"Error cleaning up temporary directory for {output_id}: {e}")

def rename_contig_ids(output_fa):
    """后处理函数，重新命名contigs ID为output_id_id_数字格式"""
    temp_file = output_fa + '.tmp'
    barcode_counter = {}
    
    try:
        with open(output_fa, 'r') as f_in, open(temp_file, 'w') as f_out:
            for line in f_in:
                if line.startswith('>'):
                    # 提取原始output_id
                    original_id = line[1:].strip()
                    # 如果是单端数据，可能包含换行符
                    original_id = original_id.split('\n')[0]
                    
                    # 判断是否有umi部分（检查下划线数量）
                    parts = original_id.split('_')
                    if len(parts) >= 2:
                        # 有umi的情况：格式为barcode_umi
                        barcode_part = parts[0]
                        remaining_part = '_'.join(parts[1:])  # 保留umi部分
                    else:
                        # 没有umi的情况：整个就是barcode
                        barcode_part = original_id
                        remaining_part = ""
                    
                    # 更新计数器
                    if barcode_part in barcode_counter:
                        barcode_counter[barcode_part] += 1
                    else:
                        barcode_counter[barcode_part] = 1
                    
                    # 构建新ID
                    if remaining_part:
                        new_id = f"{barcode_part}_{barcode_counter[barcode_part]}"
                    else:
                        new_id = f"{barcode_part}_{barcode_counter[barcode_part]}"
                    
                    f_out.write(f">{new_id}\n")
                else:
                    f_out.write(line)
        
        # 替换原文件
        shutil.move(temp_file, output_fa)
    except Exception as e:
        logging.error(f"Error renaming contig IDs: {e}")
        if os.path.exists(temp_file):
            os.remove(temp_file)


def convert_catt_to_fasta(input_csv, output_fa):
    """将CATT的CSV输出转换为FASTA格式"""
    try:
        with open(input_csv, 'r') as f_in, open(output_fa, 'w') as f_out:
            reader = csv.DictReader(f_in)
            for i, row in enumerate(reader):
                # 生成序列ID
                seq_id = f"contig_{i+1}"
                # 获取CDR3序列
                cdr3_seq = row.get('cdr3_aa', '')
                # 获取核酸序列
                nt_seq = row.get('cdr3_nt', '')
                
                # 优先写入核酸序列
                if nt_seq:
                    f_out.write(f">{seq_id}\n{nt_seq}\n")
                elif cdr3_seq:
                    f_out.write(f">{seq_id}_translated\n{cdr3_seq}\n")
    except Exception as e:
        logging.error(f"Error converting CATT output: {e}")
        raise

def parallel_batch(r1_file, threads, output_dir, r2_file=None, output_fa=None):
    """简化版本，完全依赖自动检测"""
    # 初始化输出
    os.makedirs(output_dir, exist_ok=True)
    output_fa = output_fa or os.path.join(output_dir, 'assembled_contigs.fa')
    with open(output_fa, 'w'): pass

    # 自动检测数据类型
    barcode_umi_list, is_bulk = get_barcode_umi(r1_file)
    mode = 'paired-end' if r2_file else 'single-end'
    logging.info(f"自动检测结果: {'bulk' if is_bulk else 'single-cell'} {mode}")

    if is_bulk:
        # Bulk处理流程
        catt_output = os.path.join(output_dir, 'catt_output')
        os.makedirs(catt_output, exist_ok=True)
        
        run_catt(
            f1=r1_file,
            f2=r2_file,
            output_path=catt_output,
            species="hs"

        )
    else:
        # 单细胞处理流程
        if not barcode_umi_list:
            raise ValueError("自动检测到单细胞数据但未解析出barcode")
        else:
            total_barcodes = len(barcode_umi_list)
            batches = [barcode_umi_list[i:i + threads] for i in range(0, total_barcodes, threads)]

            # 定义 mode 变量
            mode = 'paired-end' if r2_file is not None else 'single-end'
            for batch_index, batch in enumerate(batches):
                logging.info(f"Processing {mode} batch {batch_index + 1}/{len(batches)}")
                monitor_memory()
                subsets = [batch[i::threads] for i in range(threads)]
                with multiprocessing.Pool(processes=threads) as pool:
                    if r2_file is not None:
                        pool.starmap(
                            sc_assembly,
                            [(subset, r1_file, output_dir, output_fa, r2_file) for subset in subsets]
                        )
                    else:
                        pool.starmap(
                            sc_assembly,
                            [(subset, r1_file, output_dir, output_fa) for subset in subsets]
                        )

    rename_contig_ids(output_fa)
def main():
    r1_file = '/data/zhanqh/data_all/GSA/HRR1373633/trinity/bclustered1.fastq'
    r2_file = '/data/zhanqh/data_all/GSA/HRR1373633/trinity/bclustered2.fastq'
    output_dir = '/data/zhanqh/data_all/GSA/HRR1373633/trinity/TRinity'
    threads = 8
    parallel_batch(r1_file, threads, output_dir, r2_file, None)

if __name__ == "__main__":
    main()