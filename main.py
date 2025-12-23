import os
import pandas as pd
import timestamp
import argparse
import tempfile
import Alignment
import filtering
import idreb
import trinity
import annotation
import out
import time  # 导入time模块

tprint = timestamp.print_timestamp

def main():
    # Record program start time
    start_time = time.time()

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Process FASTQ files, perform BWA alignment, and separate forward/reverse strands",
        formatter_class=argparse.RawTextHelpFormatter,  # Allow line breaks and formatting
        epilog="""
Usage Examples:
    Single-end sequencing: python script.py -rf ref.fa -r1 read1.fq -o output_dir -i index_dir -t 8
    Paired-end sequencing: python script.py -rf ref.fa -r1 read1.fq -r2 read2.fq -pe -o output_dir -i index_dir -t 8
    Single-cell mode: python script.py -rf ref.fa -r1 read1.fq -sc -bc "1:16" -umi 17:26 -o output_dir -i index_dir -t 8

Detailed Description:
    This program processes FASTQ files, performs sequence alignment, and supports single-cell data analysis.
    Supports both single-end and paired-end sequencing data, with barcode and UMI extraction capabilities.
        """
    )
    
    # Group parameters by functionality (optional, for better clarity)
    parser.add_argument_group('Required Parameters')
    parser.add_argument('-rf', required=True, 
                       help="Path to reference genome file")
    parser.add_argument('-r1', required=True,
                       help="Path to first read file (read1.fq)")
    parser.add_argument('-o', type=str, required=True,
                       help="Output directory for SAM files")
    parser.add_argument('-i', type=str, required=True,
                       help="Directory to store Bowtie2 index files")
    parser.add_argument('-t', type=int, required=True,
                       help="Number of threads to use")
    
    parser.add_argument_group('Paired-end Sequencing Parameters')
    parser.add_argument('-pe', action='store_true',
                       help='Data is paired-end sequencing')
    parser.add_argument('-r2',
                       help="Path to second read file (read2.fq), required for paired-end sequencing")
    
    parser.add_argument_group('Single-cell Analysis Parameters')
    parser.add_argument('-sc', action='store_true',
                       help="Enable single-cell mode (requires barcode and UMI parameters)")
    parser.add_argument('-bc',
                       help="Barcode start and end positions, format 'start:end' (required for single-cell mode)")
    parser.add_argument('-umi', type=int, default=0,
                       help="UMI start and end positions, format like 17:26 (default: %(default)s)")
    
    parser.add_argument_group('Other Parameters')
    parser.add_argument('-l', default='alignment.log',
                       help="Path to log file (default: %(default)s)")
    parser.add_argument("-IMGT", type=str, default="/data/zhanqh/tribit/ref/human_IMGT_T.fa",
                       help="Path to IMGT reference genome (default: %(default)s)")
    
    # Add custom detailed help parameter
    parser.add_argument('--help-full', action='store_true',
                       help='Display detailed help information and usage examples')
    
    args = parser.parse_args()
    
    # Handle custom help
    if args.help_full:
        print("""
==============================================
FASTQ Processing and Alignment Tool - Detailed Help
==============================================

Functional Description:
    This tool is specifically designed for processing high-throughput sequencing data, providing a complete analysis pipeline:
    1. FASTQ file quality control and preprocessing
    2. Sequence alignment to reference genome (BWA/Bowtie2)
    3. Single-cell data analysis support (barcode/UMI extraction)
    4. Forward/reverse strand separation

Parameter Detailed Explanation:
    [Required Parameters]
    -rf <file>   : Path to reference genome FASTA file
    -r1 <file>   : Path to Read1 FASTQ file
    -o  <dir>    : Output directory
    -i  <dir>    : Bowtie2 index directory
    -t  <int>    : Number of parallel threads (recommended: 4-16)
    
    [Paired-end Sequencing]
    -pe          : Specify as paired-end sequencing data
    -r2 <file>   : Path to Read2 FASTQ file (use with -pe)
    
    [Single-cell Mode]
    -sc          : Enable single-cell analysis mode
    -bc <pos>    : Barcode position, format "start:end", e.g., "1:16"
    -umi <pos>   : UMI position, format "start:end", e.g., "17:26"
    
    [Advanced Options]
    -l  <file>   : Path to log file
    -IMGT <file> : Path to IMGT reference file

Typical Workflow:
    1. Prepare reference genome and index
    2. Run alignment: python script.py -rf ref.fa -r1 reads.fq -o output -i index -t 8
    3. Check log file: tail -f alignment.log
    4. Result files are located in output directory

Important Notes:
    - Ensure sufficient disk space for intermediate files
    - Paired-end sequencing requires both -pe and -r2 parameters
    - Single-cell mode requires correct barcode and UMI position specification
    - Recommended to use -t parameter to fully utilize multi-core CPU

Version Information: v1.0.0
Author: Your Name
Contact: email@example.com
        """)
        sys.exit(0)
    # 新增检查：如果指定-pe，则必须提供-r2参数
    if args.pe and not args.r2:
        parser.error("When --pe is specified, -r2 is required!")

    fq1 = args.r1
    fq2 = args.r2
    rf = args.rf
    index_dir = args.i
    outdir = args.o
    threads = args.t
    is_paired = args.pe
    has_barcode_umi = args.sc
    imgt = args.IMGT
    if args.bc:
        barcode_len = int(args.bc)
    else:
        barcode_len = None  # 或默认值

    if args.umi:
        umi_len = int(args.umi)
    else:
        umi_len = None  # 或默认值

    # 生成或使用已有索引
    success, index_prefix = Alignment.GBI(rf, index_dir)
    if success:
                # 创建新的子文件夹（例如：'alignment_results'）
        #subdir = os.path.join(outdir, "my_out")
        os.makedirs(outdir, exist_ok=True)  # 如果文件夹不存在则创建

        # 更新路径以保存 SAM 文件
        sam_output_path = os.path.join(outdir, "output.sam")

        # 在新的子文件夹中保存文件
        tprint(f"Saving alignment result to {sam_output_path}")

        # 调用run_bwa函数，传递具体的SAM文件路径
        sam_output_path = Alignment.run_bwa(index_prefix, 
                                            fq1, 
                                            fq2, 
                                            output_sam=sam_output_path)

    if sam_output_path:
        tprint(f"Alignment was completed, result was saved to {sam_output_path}")
        # 检查输出路径下的以 .sam 为后缀的文件
        input_sam = [f for f in os.listdir(outdir) if f.endswith('.sam')]
        tprint("The sam_file was found, filtering and translation is beginning...")
        if input_sam:
            input_sam_path = os.path.join(outdir, input_sam[0])  # 取第一个 .sam 文件

            if args.pe and fq2:#双端
                if args.bc:
                    # 创建临时文件来接收输出
                    with tempfile.NamedTemporaryFile(delete=False, suffix='_1.fastq') as temp_fq1, \
                        tempfile.NamedTemporaryFile(delete=False, suffix='_2.fastq') as temp_fq2:

                        # 调用过滤函数，输出到临时文件
                        filtering.filter_bwa_sam(input_sam_path, temp_fq1.name, temp_fq2.name)

                        # 读取临时文件内容
                        with open(temp_fq1.name, 'r') as f1, open(temp_fq2.name, 'r') as f2:
                            tfq1 = f1.readlines()  # 使用 readlines() 按行读取
                            tfq2 = f2.readlines()
                        tprint("Filtering was complete! ID rebuilding and clustering...")
                        if args.umi:
                            # 重建序列id
                            f, r = idreb.parse_fq(temp_fq1.name, barcode_len, umi_len, temp_fq2.name)

                            # 对具有相同形态的序列进行聚类
                            clustered = idreb.cluster_by_barcode_umi(f, r, umi_len)
                            tprint("Clustering done! Fastqs are written...")
                        else:
                            f, r = idreb.parse_fq(temp_fq1.name, barcode_len, umi_len, temp_fq2.name)

                            # 对具有相同形态的序列进行聚类
                            clustered = idreb.cluster_by_barcode_umi(f, r, umi_len)
                            tprint("Clustering done! Fastqs are written...")
                        # 写出新的双端文件                
                        output_fq1_path = os.path.join(outdir, "clustered1.fastq")
                        output_fq2_path = os.path.join(outdir, "clustered2.fastq")
                        idreb.write_clustered_fq(clustered, output_fq1_path, output_fq2_path)
                        tprint("Fastq files were created successfully! Assembly will be initiated...")

                        # 调用组装模块
                        if os.path.exists(output_fq1_path) and os.path.exists(output_fq2_path):
                            trinity.parallel_batch(output_fq1_path, threads, outdir, output_fq2_path, has_barcode_umi= has_barcode_umi)
                            # 删除SAM文件
                            os.remove(input_sam_path)
                            tprint(f"Removed SAM file: {input_sam_path}")
                        else:
                            tprint("Error: clustered1.fastq or clustered2.fastq not found in output directory.")
                else:
                    # 双端 bulk RNA-seq 数据，跳过聚类步骤
                    # 定义输出文件路径
                    output_fq1_path = os.path.join(outdir, "clustered1.fastq")
                    output_fq2_path = os.path.join(outdir, "clustered2.fastq")
                    
                    # 调用过滤函数，直接写入输出文件
                    filtering.filter_bwa_sam(input_sam_path, output_fq1_path, output_fq2_path)
                    
                    # 调用组装模块
                    if os.path.exists(output_fq1_path) and os.path.exists(output_fq2_path):
                        tprint("Paired-end bulk RNA-seq data detected, skipping clustering and directly calling assembly...")
                        trinity.parallel_batch(output_fq1_path, threads, outdir, output_fq2_path, has_barcode_umi= has_barcode_umi)
                        # 删除SAM文件
                        os.remove(input_sam_path)
                        tprint(f"Removed SAM file: {input_sam_path}")
                    else:
                        tprint("Error: clustered1.fastq or clustered2.fastq not found in output directory.")

            else:#单端
                if args.bc: 
                    if args.umi:#单细胞
                        with tempfile.NamedTemporaryFile(delete=False, suffix='_s.fastq') as temp_fq1:
                                                    # 调用过滤函数，输出到临时文件
                            filtering.filter_bwa_sam(input_sam_path, temp_fq1.name)
                            # 对具有相同形态的序列进行聚类
                                        # 解析 FASTQ 文件并获取 f
                            f = idreb.parse_fq(temp_fq1.name, barcode_len, umi_len)
                            clustered = idreb.cluster_by_barcode_umi(f)
                            tprint("Clustering done! Fastqs are written...")
                            output_fq1_path = os.path.join(outdir, "clustered1.fastq")
                            idreb.write_clustered_fq(clustered, output_fq1_path)
                            tprint("Fastq files were created successfully! Assembly will be initiated...")

                            # 调用组装模块
                            if os.path.exists(output_fq1_path):
                                tprint(f"Using clustered fastq files: {output_fq1_path}")
                                trinity.parallel_batch(output_fq1_path, threads, outdir, r2_file=None, has_barcode_umi=has_barcode_umi)
                                # 删除SAM文件
                                os.remove(input_sam_path)
                                tprint(f"Removed SAM file: {input_sam_path}")
                    else:
                        with tempfile.NamedTemporaryFile(delete=False, suffix='_s.fastq') as temp_fq1:
                                                    # 调用过滤函数，输出到临时文件
                            filtering.filter_bwa_sam(input_sam_path, temp_fq1.name)
                            # 对具有相同形态的序列进行聚类
                                        # 解析 FASTQ 文件并获取 f
                            f = idreb.parse_fq(temp_fq1.name, barcode_len, umi_len)
                            clustered = idreb.cluster_by_barcode_umi(f)
                            tprint("Clustering done! Fastqs are written...")
                            output_fq1_path = os.path.join(outdir, "clustered1.fastq")
                            idreb.write_clustered_fq(clustered, output_fq1_path)
                            tprint("Fastq files were created successfully! Assembly will be initiated...")

                            # 调用组装模块
                            if os.path.exists(output_fq1_path):
                                tprint(f"Using clustered fastq files: {output_fq1_path}")
                                trinity.parallel_batch(output_fq1_path, threads, outdir, r2_file=None, has_barcode_umi=has_barcode_umi)
                                # 删除SAM文件
                                os.remove(input_sam_path)
                                tprint(f"Removed SAM file: {input_sam_path}")
                else:
                    # 单端 bulk RNA-seq 数据，跳过聚类步骤
                    # 定义输出文件路径
                    output_fq1_path = os.path.join(outdir, "clustered1.fastq")
                    
                    # 调用过滤函数，直接写入输出文件
                    filtering.filter_bwa_sam(input_sam_path, output_fq1_path)
                    
                    # 调用组装模块
                    if os.path.exists(output_fq1_path):
                        # 单端bulkRNAseq数据直接调用组装
                        tprint("Single-end data detected, skipping clustering and directly calling assembly...")
                        trinity.parallel_batch(output_fq1_path, threads, outdir, r2_file=None, has_barcode_umi= False)
                        # 删除SAM文件
                        os.remove(input_sam_path)
                        tprint(f"Removed SAM file: {input_sam_path}")
            # 检查输出路径下是否存在 .fa 文件
            fa_file = [i for i in os.listdir(outdir) if i.endswith('.fa')]
            if fa_file:
                tprint(f"Found {len(fa_file)} .fa.")
                fa_file_path = os.path.join(outdir, fa_file[0]) 
                annotation.annotation(
                        imgt,
                        fa_file_path,
                        outdir,
                        threads
                    )
            else: 
                tprint(f"fasta_file is not found ,igblast error")
        else:
            tprint("No SAM files found for filtering.")
    else:
        tprint("Alignment failed.")

    # 记录程序结束时间并计算总耗时
    end_time = time.time()
    total_time = end_time - start_time
    tprint(f"Total execution time: {total_time:.2f} seconds.")

if __name__ == "__main__":
    main()
