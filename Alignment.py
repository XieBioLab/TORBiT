import subprocess
import os

def GBI(ref_file, index_dir):
    # 确保索引目录存在
    if not os.path.exists(index_dir):
        os.makedirs(index_dir)

    index_prefix = os.path.join(index_dir, os.path.basename(ref_file))

    # 检查索引文件是否已存在
    required_files = [".bwt", ".pac", ".ann", ".amb", ".sa"]
    if all(os.path.exists(index_prefix + ext) for ext in required_files):
        print("Index already, alignment begin。")
        return True, index_prefix  # 返回两个值

    # 生成索引
    try:
        subprocess.run(["bwa", "index", ref_file], check=True)
        print("Index generation successfully。")
        return True, index_prefix  # 返回两个值
    except subprocess.CalledProcessError as e:
        print(f"Index generation failed: {e}")
        return False, None  # 返回两个值


def run_bwa(index_prefix, fq1_file, fq2_file=None, output_sam="output.sam", threads=8, filter_unmapped=True):
    try:
        command = [
            "bwa", "mem",
            "-t", str(threads),
            "-M",
            index_prefix,
            fq1_file
        ]
        if fq2_file:
            command.append(fq2_file)

        # 根据 filter_unmapped 参数决定是否过滤未比对读段
        if filter_unmapped:
            # 使用管道将 BWA 输出传递给 samtools 过滤
            bwa_process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL
            )
            
            samtools_command = [
                "samtools", "view",
                "-F", "4",  # 过滤未比对读段
                "-h",       # 保留头信息
                "-o", output_sam,
                "-"         # 从标准输入读取
            ]
            
            subprocess.run(
                samtools_command,
                stdin=bwa_process.stdout,
                check=True
            )
        else:
            # 原始逻辑：直接输出所有读段
            with open(output_sam, "w") as sam_file:
                subprocess.run(
                    command,
                    stdout=sam_file,
                    stderr=subprocess.DEVNULL,
                    check=True
                )
        return output_sam
    except subprocess.CalledProcessError as e:
        print(f"Alignment failed: {e}")
        return None
    
