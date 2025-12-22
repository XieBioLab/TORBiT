# 简介

TORBiT（TCR Reconstruction-Based Immunophenotyping Tool）是一个专为从单细胞转录组测序数据中重建 T 细胞受体序列并分析 TCR 库动态的计算流程。该流程整合了多种生物信息学工具与方法，实现了基于scRNA-seq 原始测序数据，在单细胞水平上完整 TCR α/β/γ/δ 链序列的重建。

# 克隆仓库

git clone https://github.com/yourusername/fastq-processor.git
cd TORBiT

# 安装依赖

pip install -r requirements.txt(未编写)

# 安装生物信息学工具依赖

mamba install -c bioconda bwa samtools Trinity TRUST4 # 如果使用conda

#### **运行示例**

# 单端测序

python main.py -rf ./ref/hg_tcr.fa -r1 ./test_data/soloseq_f1.fastq -o ./results -i ./ref/ -t 8

# 双端测序

python main.py -rf  ./ref/hg_tcr.fa -r1 ./test_data/soloseq_f1.fastq -r2 ./test_data/soloseq_r2.fastq -pe -o ./results -i ./ref/ -t 8

# 单细胞模式(单端)

python main.py -rf ./ref/hg_tcr.fa -r1 ./test_data/soloseq_f1.fastq -sc -bc 16 -o ./results -i ./ref/ -t 8

### 7. **输出文件说明**

results/
├── aligned.sam            # 比对结果SAM文件
├── forward_strand.sam     # 正链比对结果
├── reverse_strand.sam     # 反链比对结果
├── barcode_umi.csv        # 单细胞模式：barcode和UMI信息
├── alignment.log          # 运行日志
└── stats_report.txt       # 统计报告
