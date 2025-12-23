import subprocess
import os
from pathlib import Path

def annotation(ref, input, output, threads = int):
    dir = "/data/zhanqh/miniforge3/envs/torbit/bin/"
    # Prepare path
    annot = os.path.join(output, "annot.txt")
    cdr3 = os.path.join(output, "_cdr3.out")
    abundance = os.path.join(output, "clone.csv")
    clonetype = os.path.join(dir, "trust-simplerep.pl")

    try:
        print("Running annotator...")
        # etter way to handle output redirection
        with open(annot, 'w') as outfile:
            subprocess.run([
                'annotator',
                '-f', ref,
                '-a', input,
                '--fasta',
                '-t', str(threads),
                '--needReverseComplement',
                '--noImpute',
                '--outputCDR3File',
                '-o',output
            ], check=True, stdout=outfile, stderr=subprocess.PIPE)
        
        print(f"Annotation completed. Results saved to {annot}")
        
    except subprocess.CalledProcessError as e:
        print(f"Error running annotator: {e}")
        print(f"Stderr: {e.stderr.decode('utf-8')}")
        return
    
    # Step 2: Check if CDR3 output was created
    if not os.path.exists(cdr3):
        print(f"Error: Expected CDR3 output file {cdr3} not found")
        return
    
    # Step 3: Run the Perl script to generate final report
    try:
        print("Generating final report with trust-simplerep.pl...")
        perl_cmd = [
            'perl',
            clonetype,
            cdr3,
            '>', abundance  # Again, need alternative for subprocess
        ]
        
        # Better way to handle output redirection
        with open(abundance, 'w') as outfile:
            subprocess.run([
                'perl',
                clonetype,
                cdr3
            ], check=True, stdout=outfile, stderr=subprocess.PIPE)
        
        print(f"Report generation completed. Final report saved to {abundance}")
        
    except subprocess.CalledProcessError as e:
        print(f"Error running trust-simplerep.pl: {e}")
        print(f"Stderr: {e.stderr.decode('utf-8')}")
    except Exception as e:
        print(f"Unexpected error during report generation: {e}")


