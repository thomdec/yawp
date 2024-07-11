import subprocess
import sys
from glob import glob
from uuid import uuid4
import argparse
import os

def main():
    sys.tracebacklimit = 0

    help = "\nCreate a multicallable array from the directory with sample callable files and a VCF with the FILTER \nfield populated with PASS and fails or directly from the gimbleprep temporary directory)"

    parser = argparse.ArgumentParser(description=help, prog='tool', formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog,max_help_position=40))

    parser.add_argument('-g', '--gimbleprep_dir', type=str, help='Gimbleprep directory', required=False)
    parser.add_argument('-v', '--vcf', type=str, help='VCF with "PASS" and fails in the FILTER field', required=False)
    parser.add_argument('-c', '--callabe_dir', type=str, help='Directory with the sample callable bed files', required=False)
    parser.add_argument('-o', '--output', type=str, help='Output name of the multicallable array (.bed.gz)', required=True)

    args = parser.parse_args()

    if args.gimbleprep_dir is not None:

        if os.path.exists(args.gimbleprep_dir) is not True:
            raise Exception(f"[ForgeMulticallable] ERROR: The specified directory with callable sites {args.gimbleprep_dir} does not exist")
        
        if os.access(os.path.dirname("./" + args.output), os.W_OK) is not True:
            raise Exception(f"[ForgeMulticallable] ERROR: Cannot write output to {args.output}")

        vcf = (args.gimbleprep_dir + "/" if not args.gimbleprep_dir.endswith("/") else args.gimbleprep_dir) + "vcf.filtered.vcf.gz"

        if os.path.exists(vcf) is not True:
            raise Exception(f"[ForgeMulticallable] ERROR: The vcf.filtered.vcf.gz could not be found in the tmp_gimble directory")
        
        tmp_multiinter_bed = f"tmp_multiinter_bed_{uuid4().hex}.bed.gz"
        tmp_fail_bed = f"tmp_fail_bed_{uuid4().hex}.bed.gz"

        _stdout = subprocess.check_output("bcftools query -l " + vcf, shell=True, text=True)
        sample_ids = [sample_id for sample_id in _stdout.split("\n") if sample_id]
        callable_beds = list(glob((args.gimbleprep_dir + "/" if not args.gimbleprep_dir.endswith("/") else args.gimbleprep_dir) + "*.callable.bed" ))
        ordered_beds = [[x for x in callable_beds if sample_ids[i] in x][0] for i in range(len(sample_ids))]

        # creates the multi callable bed file
        run_command(f"bedtools multiinter -i {' '.join(ordered_beds)} -names {' '.join(sample_ids)} | cut -f 1-3,6- | bgzip > {tmp_multiinter_bed}")

        # creates a bed file with failed sites in the VCF
        run_command(f"""bcftools view {args.vcf} -i "%FILTER!='PASS'" | bcftools query -f '%CHROM\t%POS0\t%END\t%FILTER\n' | bgzip > {tmp_fail_bed}""")

        # remove VCF failed sites from the multi callable array
        run_command(f"""bedtools subtract -a {tmp_multiinter_bed} -b {tmp_fail_bed} | bgzip > {args.output}""")

        # index the multicallable array
        run_command(f"""tabix -p bed {args.output}""")

        run_command(f"rm {tmp_fail_bed} && rm {tmp_multiinter_bed}")

        print("Successfully created the multicallable bed file")

    else:

        if os.path.exists(args.vcf) is not True:
            raise Exception(f"[ForgeMulticallable] ERROR: The specified VCF {args.vcf} does not exist")
        
        if os.path.exists(args.callabe_dir) is not True:
            raise Exception(f"[ForgeMulticallable] ERROR: The specified directory with callable sites {args.callabe_dir} does not exist")
        
        if os.access(os.path.dirname("./" + args.output), os.W_OK) is not True:
            raise Exception(f"[ForgeMulticallable] ERROR: Cannot write output to {args.output}")

        tmp_multiinter_bed = f"tmp_multiinter_bed_{uuid4().hex}.bed.gz"
        tmp_fail_bed = f"tmp_fail_bed_{uuid4().hex}.bed.gz"

        _stdout = subprocess.check_output("bcftools query -l " + args.vcf, shell=True, text=True)
        sample_ids = [sample_id for sample_id in _stdout.split("\n") if sample_id]
        callable_beds = list(glob((args.callabe_dir + "/" if not args.callabe_dir.endswith("/") else args.callabe_dir) + "*.callable.bed" ))
        ordered_beds = [[x for x in callable_beds if sample_ids[i] in x][0] for i in range(len(sample_ids))]

        # creates the multi callable bed file
        run_command(f"bedtools multiinter -i {' '.join(ordered_beds)} -names {' '.join(sample_ids)} | cut -f 1-3,6- | bgzip > {tmp_multiinter_bed}")

        # creates a bed file with failed sites in the VCF
        run_command(f"""bcftools view {args.vcf} -i "%FILTER!='PASS'" | bcftools query -f '%CHROM\t%POS0\t%END\t%FILTER\n' | bgzip > {tmp_fail_bed}""")

        # remove VCF failed sites from the multi callable array
        run_command(f"""bedtools subtract -a {tmp_multiinter_bed} -b {tmp_fail_bed} | bgzip > {args.output}""")

        # index the multicallable array
        run_command(f"""tabix -p bed {args.output}""")

        run_command(f"rm {tmp_fail_bed} && rm {tmp_multiinter_bed}")

        print("Successfully created the multicallable bed file")

# run command in func.py
def run_command(args):
    try:
        call = subprocess.run(args, text=True, capture_output=True, shell=True)
        call.check_returncode()
    except subprocess.CalledProcessError as cpe:
        if call.stdout:
            sys.exit('[ForgeMulticallable] The command following command failed with error code %r:\n[X] => %s\n[X] (STDOUT): %r\n[X] (STDERR): %r' % (cpe.returncode, cpe.cmd, cpe.stdout, cpe.stderr))
        sys.exit('[ForgeMulticallable] The command following command failed with error code %r:\n[X] => %s\n[X] (STDERR): %r' % (cpe.returncode, cpe.cmd, cpe.stderr.rstrip("\n")))


if __name__ == "__main__":
            
   main()