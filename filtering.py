import os
import subprocess

def filter_bwa_sam(input_sam, output_fq1, output_fq2=None):
    if not os.path.exists(input_sam):
        print(f"输入文件 {input_sam} 不存在，无法继续处理。")
        return None

    try:
        # 获取 samtools 的完整路径
        samtools_path = "./samtools-1.3/samtools"  # 修改为你的 samtools 路径

        if output_fq2:  # 双端模式
            # 添加 -h 参数保留头部
            command = f"{samtools_path} view -h -F 12 {input_sam} | {samtools_path} fastq - -1 {output_fq1} -2 {output_fq2}"
        else:  # 单端模式
            # 单端模式下，只生成一个 FASTQ 文件
            command = f"{samtools_path} view -h -F 4 {input_sam} | {samtools_path} fastq - > {output_fq1}"

        # 执行命令
        result = subprocess.run(command, 
                                shell=True, 
                                check=True, 
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE, 
                                text=True)

        print("command output:", result.stdout)
        print("error information:", result.stderr)

        if output_fq2:  # 双端模式
            return output_fq1, output_fq2
        else:  # 单端模式
            return output_fq1

    except subprocess.CalledProcessError as e:
        print(f"Command error: {e.stderr}")
        return None, None if output_fq2 else None
    except Exception as e:
        print(f"An error occurred during the execution process: {e}")
        return None, None if output_fq2 else None
    

# def main():
#     filter_bwa_sam(input_sam="/data/zhanqh/soloseq_test/trinity/sample_Aligned.out.sam", 
#                    output_fq1="/data/zhanqh/soloseq_test/trinity/sample1.fastq",
#                    output_fq2="/data/zhanqh/soloseq_test/trinity/sample2.fastq")
# if __name__ == "__main__":
#     main()
