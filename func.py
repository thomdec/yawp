import allel
import numpy as np
import pandas as pd
import subprocess
from scipy.special import comb
from itertools import combinations, product, groupby
from uuid import uuid4
from io import StringIO
import sys

# To do:
#   - Rolling windows and SNP count windows
#   - Add a site missingness and genotype missingness measure! (To facilitate post-yawp filtering)

def diff_within(g, pop):
    """Returns per site allelic differences within a population"""
    ac_pop = g[:, pop].count_alleles()
    an_pop = np.sum(ac_pop, axis=1)
    n_pairs = an_pop * (an_pop - 1) / 2
    n_same = np.sum(ac_pop * (ac_pop - 1) / 2, axis=1)
    n_diff = np.sum(n_pairs - n_same)
    return(n_diff)

def diff_between(g, pop1, pop2):
    """Returns per site allelic differences between two populations"""
    ac_pop1 = g[:, pop1].count_alleles(max_allele=3)
    ac_pop2 = g[:, pop2].count_alleles(max_allele=3)
    an_pop1 = np.sum(ac_pop1, axis=1)
    an_pop2 = np.sum(ac_pop2, axis=1)
    n_pairs = an_pop1 * an_pop2 
    n_same = np.sum(ac_pop1 * ac_pop2, axis=1)
    n_diff = np.sum(n_pairs - n_same)
    return(n_diff)

def comp_within(callable_array, start_end, pop):
    """Returns the total number of nucleotide comparison within a population"""
    sum_cov_pop = np.sum(callable_array[:, pop], axis=1)
    n_comp = np.sum(comb(2 * sum_cov_pop, 2) * (start_end[:, 1] - start_end[:, 0]))
    return(n_comp)

def comp_between(callable_array, start_end, pop1, pop2):
    """Returns the total number of nucleotide comparison between two population"""
    sum_cov_pop1 = np.sum(callable_array[:, pop1], axis=1)
    sum_cov_pop2 = np.sum(callable_array[:, pop2], axis=1)
    n_comp = np.sum(4 * sum_cov_pop1 * sum_cov_pop2 * (start_end[:, 1] - start_end[:, 0]))
    return(n_comp)

def calc_pi(g, callable_array, start_end, pop):
    """Returns nucleotide diversity within a population"""

    n_diff = diff_within(g, pop)
    n_comp = comp_within(callable_array, start_end, pop)
    pi = n_diff / n_comp
    return(pi, int(n_diff), int(n_comp))

def calc_dxy(g, callable_array, start_end, pop1, pop2):
    """"Returns nucleotide divergence between two populations"""

    n_diff = diff_between(g, pop1, pop2)
    n_comp = comp_between(callable_array, start_end, pop1, pop2)
    dxy = n_diff / n_comp
    return(dxy, int(n_diff), int(n_comp))

def mask(vcf, multicallabe_bed, mask):
    """Mask a vcf and and bed file using bcftools -filter and bedtools intersect. Returns masked files paths"""
    
    masked_vcf = f"tmp.masked.{uuid4().hex}.vcf.gz"
    masked_multicallable = f"tmp.masked.multicallable.{uuid4().hex}.bed.gz"

    run_command(f"""bcftools filter -R {mask} {vcf} -Oz -o {masked_vcf}""")
    run_command(f"""bcftools index {masked_vcf}""")

    run_command(f"""bedtools intersect -a {multicallabe_bed} -b {mask} | bgzip > {masked_multicallable}""")
    run_command(f"""tabix -p bed {masked_multicallable}""")

    return masked_vcf, masked_multicallable

def create_chunked_windows(window_size, chunk_size, chr_len):
    """Create windows in chunks"""

    window_starts = np.array([*range(0, int(chr_len), window_size)])
    window_ends = np.array([*range(0 + window_size, int(chr_len) + window_size, window_size)])
    window_ends[-1] = int(chr_len)
    windows = np.stack((window_starts, window_ends)).T

    chunk_idx = np.ceil(window_ends / chunk_size).astype(int)
    chunk_starts = np.array([window_starts[chunk_idx==i].min() for i, _ in groupby(chunk_idx)])
    chunk_ends = np.array([window_ends[chunk_idx==i].max() for i, _ in groupby(chunk_idx)])
    chunks = np.stack((chunk_starts, chunk_ends)).T

    window_idx = [np.where(chunk_idx==i)[0] for i, _ in groupby(chunk_idx)]

    return windows, chunks, window_idx

def parse_bed_chunk(bed, region):
    """Parse a specific region of a bed file, using tabix"""

    cmd = f"tabix {bed} {region} | cut -f 2-"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    output = StringIO(result.stdout)
    array = np.loadtxt(output, dtype=int)
    return array

def intersect(bed_array, start, end):
    """Intersect bed array with a defined window"""

    intersect_mask = (bed_array[:, 0] < end) & (bed_array[:, 1] > start)
    intersect_array = bed_array[intersect_mask]
    if not intersect_array[0,0] == start:
        intersect_array[0,0] = start
    if not intersect_array[-1, 1] == end:
        intersect_array[-1, 1] = end

    return(intersect_array)

def compute_global(vcf, multicallabe_bed, pops, output):
    """Computes global estimates of pi and dxy. Writes to output files"""

    pi_file = open(f'{output}.pi.tsv', 'w')
    dxy_file = open(f'{output}.dxy.tsv', 'w')

    pi_file.write('pop' + '\t' +  'pi' + '\t' +  'n_diff' + '\t' +  'n_comp' + '\n')
    dxy_file.write('pop1' + '\t' +  'pop2' + '\t' +  'dxy' + '\t' + 'n_diff' + '\t' +  'n_comp' + '\n')

    callset = allel.read_vcf(vcf, fields=['calldata/GT', 'variants/POS', 'samples'])
    g = allel.GenotypeArray(allel.GenotypeDaskArray(callset['calldata/GT']))
    samples = callset['samples']

    multicallable_array = np.loadtxt(multicallabe_bed, usecols=(list(range(1, len(samples) + 3))), dtype=int)

    for pop in pops:
        pi, n_diff, n_comp = calc_pi(g, multicallable_array[:, 2:], multicallable_array[:, :2], pops[pop])
        pi_file.write(pop + '\t' +  str(pi) + '\t' +  str(n_diff) + '\t' +  str(n_comp) + '\n')

    for pop1, pop2 in combinations(pops, 2):
        dxy, n_diff, n_comp  = calc_dxy(g, multicallable_array[:, 2:], multicallable_array[:, :2], pops[pop1], pops[pop2])
        dxy_file.write(pop1 + '\t' +  pop2 + '\t' +  str(dxy) + '\t' +  str(n_diff) + '\t' +  str(n_comp) + '\n')

    pi_file.close()
    dxy_file.close()


## check if the chunk is empty!!! check if the window is empty --> n_comp = 0
def compute_windows(vcf, multicallable_bed, pops, window_size, chunk_size, output):
    """Computes windowed estimates of pi and dxy. Writes to output files"""

    pi_file = open(f'{output}.pi.tsv', 'w')
    dxy_file = open(f'{output}.dxy.tsv', 'w')

    pi_file.write('chromosome' + '\t' + 'start' + '\t' +  'end' + '\t' +  'pop' + '\t' +  'pi' + '\t' +  'n_diff' + '\t' +  'n_comp' + '\n')
    dxy_file.write('chromosome' + '\t' + 'start' + '\t' +  'end' + '\t' +  'pop1' + '\t' +  'pop2' + '\t' +  'dxy' + '\t' + 'n_diff' + '\t' +  'n_comp' + '\n')

    data = subprocess.check_output("bcftools index -s " + vcf, shell=True, text=True)
    chromosomes = pd.DataFrame([x.split('\t') for x in data[:-1].split('\n')], columns=["id", "chr_size", "vcf_line"])

    for iteration, chr in chromosomes.iterrows():

        windows, chunks, window_idx = create_chunked_windows(window_size, chunk_size, chr.chr_size)

        for i, chunk in enumerate(chunks):

            region_string = f"{chr.id}:{chunk[0]}-{chunk[1]}"

            # load genotype chunk
            callset = allel.read_vcf(vcf, region=region_string, fields=['calldata/GT', 'variants/POS', 'samples'])
            pos = callset['variants/POS']
            g = allel.GenotypeArray(allel.GenotypeDaskArray(callset['calldata/GT']))

            # load multicallable chunk
            multicallable_array = parse_bed_chunk(multicallable_bed, region_string)

            for window in windows[window_idx[i]]:
                mask_window = (pos >= window[0]) & (pos <= window[1])
                win_multicallable_array = intersect(multicallable_array, window[0], window[1])

                for pop in pops:
                    win_pi, n_diff, n_comp = calc_pi(g[mask_window], win_multicallable_array[:, 2:], win_multicallable_array[:, :2], pops[pop])
                    pi_file.write(str(chr.id) + '\t' + str(window[0]) + '\t' +  str(window[1]) + '\t' +  pop + '\t' +  str(win_pi) + '\t' +  str(n_diff) + '\t' +  str(n_comp) + '\n')

                for pop1, pop2 in combinations(pops, 2):
                    win_dxy, n_diff, n_comp  = calc_dxy(g[mask_window], win_multicallable_array[:, 2:], win_multicallable_array[:, :2], pops[pop1], pops[pop2])
                    dxy_file.write(str(chr.id) + '\t' + str(window[0]) + '\t' +  str(window[1]) + '\t' +  pop1 + '\t' +  pop2 + '\t' +  str(win_dxy) + '\t' +  str(n_diff) + '\t' +  str(n_comp) + '\n')

    pi_file.close()
    dxy_file.close()

def run_command(args):
    """Run a command as a subprocess"""

    try:
        call = subprocess.run(args, text=True, capture_output=True, shell=True)
        call.check_returncode()
    except subprocess.CalledProcessError as cpe:
        if call.stdout:
            sys.exit('[YAWP] The command following command failed with error code %r:\n[X] => %s\n[X] (STDOUT): %r\n[X] (STDERR): %r' % (cpe.returncode, cpe.cmd, cpe.stdout, cpe.stderr))
        sys.exit('[YAWP] The command following command failed with error code %r:\n[X] => %s\n[X] (STDERR): %r' % (cpe.returncode, cpe.cmd, cpe.stderr.rstrip("\n")))

