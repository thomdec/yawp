import allel
import numpy as np
import pandas as pd
from warnings import filterwarnings
import sys
import argparse
import os

import func

## TO DO
#   - Add region parameters! + only read in the samples supplied in the population file! 
#   - Make more verbose: ADD TEXT INFO output AND TIME
#   - Add a multiprocessing option (across chromosomes? or accross windows?)
#   - TRY not to use the allel dependency (for reading in VCF --> the bottleneck!)

def main():

    filterwarnings('ignore')
    sys.tracebacklimit = 0

    yawp_image = \
    " __ __   ____  __    __  ____  \n" \
    "|  |  | /    ||  |__|  ||    \ \n" \
    "|  |  ||  o  ||  |  |  ||  o  )\n" \
    "|  ~  ||     ||  |  |  ||   _/ \n" \
    "|___, ||  _  ||  `  '  ||  |   \n" \
    "|     ||  |  | \      / |  |   \n" \
    "|____/ |__|__|  \_/\_/  |__|   \n" \
                                
    help = "\nYAWP: computes common popgen stats"

    parser = argparse.ArgumentParser(description=yawp_image + help, prog='tool', formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog,max_help_position=40))

    parser._action_groups.pop()
    required = parser.add_argument_group('Required parameters')
    optional = parser.add_argument_group('Optional parameters')

    #required
    required.add_argument('-v', '--vcf', type=str, help='VCF file', required=True)
    required.add_argument('-b', '--multicallable_bed', type=str, help='Multicallable bed array (can be gzip compressed)', required=True)
    required.add_argument('-p', '--populations', type=str, help='Populations file.\nCsv file with id and population in each row: [id,pop]', required=True)
    required.add_argument('-o', '--output', type=str, help='Output name. (.tsv)', required=True)

    # optional
    #parser.add_argument('-s', '--statistics', type=str, help='pi, dxy, tajD', required=True)
    optional.add_argument('-w', '--window_size', type=int, help='Window size. If not supplied, will return global estimates', required=False)
    optional.add_argument('-m', '--mask', type=str, help='Mask (e.g 4D sites). Must be in the BED format', required=False) #should I rename mask to filter?
    optional.add_argument('-c', '--chunk_size', type=int, help='Chunk size [Default = 1_000_000]', default=1_000_000, required=False) # choose default wisely
    optional.add_argument('-k', '--keep_tmp', help='Keep temporary VCF and multicallable.bed created when using a mask', action='store_true')

    args = parser.parse_args()

    # check parameters
    print("[YAWP] Checking command line arguments")

    if os.path.exists(args.vcf) is not True:
        raise Exception(f"[YAWP] ERROR: The specified VCF {args.vcf} does not exist") 

    if not os.path.exists(args.vcf + ".csi") | os.path.exists(args.vcf + ".tbi"): 
        raise Exception('[YAWP] ERROR: The vcf is not indexed. The vcf can be indexed with "bcftools index" or "tabix -p vcf"') 

    # + check if bgzip'd + check if tabix indexed! (Should do that in prep as well! ) CHECKKKK !!! 
    if os.path.exists(args.multicallable_bed) is not True:
        raise Exception(f"[YAWP] ERROR: The specified multicallable bed file {args.multicallable_bed} does not exist")
    
    if os.access(os.path.dirname("./" + args.output), os.W_OK) is not True:
        raise Exception(f"[YAWP] ERROR: Cannot write output to {args.output}")

    if os.path.exists(args.populations) is not True:
        raise Exception(f"[YAWP] ERROR: The specified populations file {args.populations} does not exist")
    
    pop_df = pd.read_csv(args.populations, names=["id", "pop"])
    vcf_header = allel.read_vcf_headers(args.vcf)

    if not np.all(np.isin(pop_df.id, vcf_header.samples)):
        raise Exception(f"[YAWP] ERROR: The following sample in the population file could not be found in the VCF:\n{', '.join(pop_df.id[~np.isin(pop_df.id, vcf_header.samples)])}")

    pops = {}
    unique_pops = np.unique(pop_df["pop"]).astype(str)
    for pop in unique_pops:
        ids = pop_df[pop_df['pop'] == pop]['id']
        pops[pop] = np.where(np.isin(vcf_header.samples, ids))[0]
        
    if args.mask is not None:

        if args.window_size is not None:
            raise Exception("[YAWP] ERROR: Masks are not supported with windows at the moment! Come back soon")

        print("[YAWP]: Masking the VCF and the multicallable bed file")
        if os.path.exists(args.mask) is not True:
            raise Exception(f"[YAWP] ERROR: The specified mask {args.mask} does not exist") 

        vcf, multicallable_bed = func.mask(args.vcf, args.multicallable_bed, args.mask)

        if os.stat(multicallable_bed).st_size == 28:
            raise Exception("[YAWP] ERROR: The filtered multicallable bed and vcf could not be created. Check the mask bed file")

    else: 

        vcf, multicallable_bed = args.vcf, args.multicallable_bed

    if args.window_size is not None:
        
        if not os.path.exists(args.multicallable_bed + ".tbi"): # or + "tabix index! "
            raise Exception('[YAWP] ERROR: The multicallable bed is not indexed. It can be indexed with or "tabix -p bed"') 
         
        print("[YAWP] YAWPING")
        func.compute_windows(args.vcf, args.multicallable_bed, pops, args.window_size, args.chunk_size, args.output)
        print("[YAWP] Successfully ran YAWP in windows")

    else:

        print("[YAWP] Global estimation started")
        func.compute_global(vcf, multicallable_bed, pops, args.output)
        print("[YAWP] Successfully ran YAWP global")

    if args.mask is not None: 
        
        if args.keep_tmp:
            
            print(f"[YAWP] Keeping tmp files: {multicallable_bed} & {vcf}")
        
        else:

            func.run_command(f"rm {multicallable_bed} && rm {multicallable_bed}.tbi && rm {vcf} && rm {vcf}.csi")

if __name__ == "__main__":
            
   main()