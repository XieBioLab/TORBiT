import os
import subprocess

def filter_bwa_sam(input_sam, output_fq1, output_fq2=None):
    if not os.path.exists(input_sam):
        print(f"The input file {input_sam} not exist, processing cannot continue.")
        return None

    try:
        samtools_path = "./samtools-1.3/samtools" 

        if output_fq2:  
            command = f"{samtools_path} view -h -F 12 {input_sam} | {samtools_path} fastq - -1 {output_fq1} -2 {output_fq2}"
        else:  
            command = f"{samtools_path} view -h -F 4 {input_sam} | {samtools_path} fastq - > {output_fq1}"

        result = subprocess.run(command, 
                                shell=True, 
                                check=True, 
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE, 
                                text=True)

        print("command output:", result.stdout)
        print("error information:", result.stderr)

        if output_fq2:  
            return output_fq1, output_fq2
        else:  
            return output_fq1

    except subprocess.CalledProcessError as e:
        print(f"Command error: {e.stderr}")
        return None, None if output_fq2 else None
    except Exception as e:
        print(f"An error occurred during the execution process: {e}")
        return None, None if output_fq2 else None
    


