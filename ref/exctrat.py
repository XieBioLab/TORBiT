from Bio import SeqIO
import re

def extract_tcr_genes(input_fa, output_fa, keywords):
    pattern = re.compile("|".join(keywords), re.IGNORECASE)
    with open(output_fa, "w") as out_file:
        for record in SeqIO.parse(input_fa, "fasta"):
            if pattern.search(record.description):
                # 关键修改：添加 wrap=None 保持单行
                SeqIO.write(record, out_file, "fasta-2line")  # 或使用以下自定义写入
                # 或者手动写入（完全控制格式）：
                # out_file.write(f">{record.description}\n{str(record.seq)}\n")
    print(f"筛选完成！结果已保存到 {output_fa}")



# 用法示例
if __name__ == "__main__":
    input_fa = "/data/zhanqh/py_pre-test/ref_all/hg38_tcr.fa"          # 输入文件路径
    output_fa = "/data/zhanqh/py_pre-test/ref_all/hg38n_tcr.fa"    # 输出文件路径
    keywords = [
        "TRAV", "TRBV", "TRDV", "TRGV",  # TCR V 基因
        "TRAJ", "TRBJ", "TRDJ", "TRGJ",  # TCR J 基因
        "TRAC", "TRBC", "TRDC", "TRGC",  # TCR C 基因
    ]
extract_tcr_genes(input, output_fa, keywords)