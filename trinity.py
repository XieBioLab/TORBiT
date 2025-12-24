import os
import logging
import shutil
import subprocess
from collections import defaultdict
from multiprocessing import Lock
import psutil
import multiprocessing

# 初始化日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# 线程安全锁
lock = Lock()

def get_barcode_umi(r1_files):
    if isinstance(r1_files, str):  # 如果是单个文件路径
        r1_files = [r1_files]

    barcode_umi_counts = defaultdict(int)
    for r1_file in r1_files:
        with open(r1_file, 'r') as r1:
            while True:
                r1_header = r1.readline().strip()
                if not r1_header:
                    break
                
                barcode = None
                umi = None
                
                # 新的识别逻辑：定位到" 1:"的位置
                key_pos = r1_header.find(" ")
                if key_pos == -1:
                    logging.warning(f"Header does not contain ' ' pattern: {r1_header}")
                    # 跳过后续三行
                    r1.readline()  # 序列行
                    r1.readline()  # '+' 行
                    r1.readline()  # 质量值行
                    continue
                
                # 提取" 1:"后面的部分
                barcode_umi_part = r1_header[key_pos + 3:]  # 跳过" 1:"
                
                # 分割冒号获取barcode和umi
                parts = barcode_umi_part.split(':')
                
                # 第一个冒号前的内容是barcode
                if len(parts) >= 1:
                    barcode = parts[0]
                    
                    # 如果有第二个冒号，则后面的内容是umi
                    if len(parts) >= 2:
                        umi = parts[1].split()[0]  # 去除可能存在的后续空格和内容
                
                if barcode is not None:
                    barcode_umi = (barcode, umi) if umi is not None else (barcode, None)
                    barcode_umi_counts[barcode_umi] += 1
                else:
                    logging.warning(f"Could not extract barcode from header: {r1_header}")
                
                # 跳过后续三行
                r1.readline()  # 序列行
                r1.readline()  # '+' 行
                r1.readline()  # 质量值行

    return list(barcode_umi_counts.keys())

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

def parallel_batch(r1_file, threads, output_dir, r2_file=None, has_barcode_umi=bool, output_fa=None):
    # 确保output_dir存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 设置输出文件路径（在output_dir目录下）
    if not output_fa:
        output_fa = os.path.join(output_dir, 'assembled_contigs.fa')
    
    # 初始化输出文件（清空或创建）
    with open(output_fa, 'w') as f:
        pass  # 只是创建/清空文件

    if has_barcode_umi:
        # 如果存在 barcode 和 umi，执行原逻辑
        barcode_umi_list = get_barcode_umi(r1_file)
        total_barcodes = len(barcode_umi_list)
        
        if total_barcodes == 0:
            logging.error("No barcode information found in the input file")
            return
            
        # 检查是否有umi信息
        has_umi = any(umi is not None for (barcode, umi) in barcode_umi_list)
        logging.info(f"Found {total_barcodes} barcode{' and umi' if has_umi else ''} combinations")
        
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

    else:
        # 如果不存在 barcode 和 umi，直接处理 bulk 数据
        mode = 'paired-end' if r2_file is not None else 'single-end'
        logging.info(f"Processing {mode} bulk data using {threads} threads")

        # 创建输出目录
        bulk_output_dir = os.path.join(output_dir, 'trinity_output')
        os.makedirs(bulk_output_dir, exist_ok=True)

        # 使用多进程处理 bulk 数据
        try:
            logging.info(f"Starting Trinity assembly for bulk data")
            if r2_file:
                logging.info(f"Running Trinity in paired-end mode")
                success = run_trinity(r1_file, r2_file, bulk_output_dir)
            else:
                logging.info(f"Running Trinity in single-end mode")
                success = run_trinity(r1_file, None, bulk_output_dir)
            
            if not success:
                logging.error("Trinity failed to run for bulk data")
                return
            
            # 检查Trinity输出文件
            trinity_output = os.path.join(bulk_output_dir, 'Trinity.fasta')
            if not os.path.exists(trinity_output):
                # 尝试查找其他可能的输出位置
                found = False
                for root, dirs, files in os.walk(bulk_output_dir):
                    if 'Trinity.fasta' in files:
                        trinity_output = os.path.join(root, 'Trinity.fasta')
                        found = True
                        break
                
                if not found:
                    logging.error("Could not find Trinity output file")
                    return

            # 直接将Trinity输出复制到最终文件
            shutil.copy(trinity_output, output_fa)
            logging.info(f"Bulk assembly completed successfully")

        except Exception as e:
            logging.error(f"Error processing bulk data: {e}")
            return

    # 在所有处理完成后，重命名contig IDs
    rename_contig_ids(output_fa)
    logging.info(f"Assembly completed for {mode} data.")

