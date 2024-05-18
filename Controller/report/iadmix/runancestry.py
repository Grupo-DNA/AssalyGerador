#! /usr/bin/env python2.7
# author Vikas Bansal vibansal@ucsd.edu

# TODO
# 1. print # of sites and average read-depth per site in calculateGLL code
# 3. option to flip alleles for - strand SNPs
# 4. option to delete intermediate ancestry files created
# 5. option to match SNPs using rsid or chr:position
# 6. option to read VCF file with one or more individuals
# 7. option to read a file with list of BAM files for sequencing experiment
# 8. option to specify --delta threshold for parsimonious ancestry estimation
# 9. along with admixture estimates: also print delta for each value...
# 10. read plink .bed files binary format for genotypes
# 11. option to do ancestry for only specific samples in ped file

import glob
import math
import os
import string
import sys
import time
from math import log
from optparse import OptionParser
from os import PathLike
from pathlib import Path
from subprocess import DEVNULL, call
from typing import Union

POOLSIZE = 2
INCLUDE_STRAND_AMBIG = 0
CORES = 1
WINDOW = 0
PARSIMONY = 1
HWE_CHECK = 0
OUTPUT_GLL_ONLY = 0
ADD_chr_names = False
SORT_SNPS = 1  # sort the 'forGLL' file before running calculateGLL code, added 07/07/2017

##################################################################################################################################

DIRNAME = Path(__file__).resolve().parent


def get_file_encoding(filepath: str):
    """Chooses a possible encoding for .csv files
    Only between UTF-8 and Latin-1 (ISO 8859-1)

    Args:
        filepath (str): .csv filepath.

    Returns:
        str: Encoding found.
    """
    for encoding in ('utf-8', 'latin-1'):
        try:
            with open(filepath, "r", newline="", encoding=encoding) as handle:
                handle.read()
            return encoding
        except UnicodeDecodeError:
            pass
    raise ValueError(f"Failed to identify encoding for file '{filepath}'")


def output_GLL_inputfile(afmatrixfile, outfile):
    File1 = open(outfile, 'w')
    snp_list = []
    with open(afmatrixfile) as File:
        for line in File:
            if line[0] != '#':
                snp = line.strip().split()
                if not ADD_chr_names:
                    print('SNP', snp[0].strip('chr'), snp[1], snp[3],
                          snp[4], snp[3] + '/' + snp[4], 100, file=File1)
                elif 'chr' in snp[0]:
                    print('SNP', snp[0], snp[1], snp[3], snp[4],
                          snp[3] + '/' + snp[4], 100, file=File1)

    File1.close()
    if SORT_SNPS == 1:
        call('sort -k 2,2 -k 3,3g ' + outfile +
             ' > ' + outfile + '.sorted', shell=True)
        call('mv ' + outfile + '.sorted ' + outfile, shell=True)


def make_ancestry_inputfile_simple(afmatrixfile, GLLfile, outfile):
    File1 = open(outfile, 'w')
    GLLs = {}

    encod = get_file_encoding(GLLfile)
    with open(GLLfile, encoding=encod) as File:
        for line in File:
            snp = line.strip().split()
            GLLs[(snp[POOLSIZE+1].strip('chr'), int(snp[POOLSIZE+2]))
                 ] = ' '.join(snp[0:POOLSIZE+1])

    encod = get_file_encoding(afmatrixfile)
    with open(afmatrixfile, encoding=encod) as File:
        for line in File:
            snp = line.strip().split()
            if line[0] == '#':
                print('#GLL', end=' ', file=File1)
                for i in range(5, len(snp)):
                    print(snp[i], end=' ', file=File1)
                print('\n', end=' ', file=File1)
            else:
                try:
                    ll = GLLs[(snp[0].strip('chr'), int(snp[1]))]
                    print(ll, end=' ', file=File1)
                    for i in range(5, len(snp)):
                        print(snp[i], end=' ', file=File1)
                    print('\n', end=' ', file=File1)
                except KeyError:
                    pass
    File1.close()

# function that processes bam file for ancestry calculation


def ancestry_pipeline_bam(matrixfile, bamfile, outfile):

    print("making input file for getting genotype likelihoods", file=sys.stderr)
    output_GLL_inputfile(matrixfile, outfile+'.forGLL')

    print("calculating genotype likelihoods for bam file",
          bamfile, "poolsize is: ", POOLSIZE, file=sys.stderr)
    if OUTPUT_GLL_ONLY == 1:
        call([f"\"{DIRNAME!s}/calculateGLL\" --allsites 1 -p {POOLSIZE!r} --bam " +
             bamfile + " --variants " + outfile + ".forGLL" + " --mmq 30 > " + outfile + ".GLL"], shell=True)
    else:
        call([f"\"{DIRNAME!s}/calculateGLL\" -p {POOLSIZE!r} --bam " + bamfile +
             " --variants " + outfile + ".forGLL" + " --mmq 30 > " + outfile + ".GLL"], shell=True)

        print("making input file for ancestry calculations", file=sys.stderr)
        make_ancestry_inputfile_simple(
            matrixfile, outfile + ".GLL", outfile + ".ancestry.input")
        # python runancestry.py ancestry-data/hapmap3.allchroms.shared.matrix.hg19 ancestry-data/DM00231.bwamem.MD.GLL > ancestry-data/DM00231.bwamem.MD.input
        print("ancestry admixture calculations for", bamfile,
              "poolsize is ", POOLSIZE, file=sys.stderr)
        call([f"time \"{DIRNAME!s}/ANCESTRY\" -p {POOLSIZE!r} " +
             "--pr 1 -i " + outfile + ".ancestry.input > " + outfile + '.ancestry.out'], shell=True)
        # call(["rm -f " + bamfile + ".GLL"],shell=True);
        # call(["rm -f " + bamfile + ".forGLL"],shell=True);
        # call(["rm -f " + bamfile + ".ancestry.input"],shell=True);


##################################################################################################################################

# add support for plink ped files for multiple samples or processing a file with list of bamfile names
# function that reads in PLINK ped and map files and determines ancestry for each individual in ped file
def make_ancestry_inputfile_plink(pedfile, mapfile, afmatrixfile, outfile, modulo_CORE):

    RC = {}
    RC['A'] = 'T'
    RC['T'] = 'A'
    RC['C'] = 'G'
    RC['G'] = 'C'
    gll_low = -4
    gll_high = 0
    gll_equal = -0.477121

    # read the allele frequency matrix file into a hashtable once for all samples
    print("reading allele frequency file", afmatrixfile, file=sys.stderr)

    lines = 0
    AFtable = {}
    headerline = []
    with open(afmatrixfile) as File:
        for line in File:
            snp = line.strip().split()
            if line[0] == '#' or lines == 0:
                headerline = snp
            elif INCLUDE_STRAND_AMBIG == 1:
                AFtable[snp[2]] = snp
            elif snp[3] == 'A' and snp[4] == 'T':
                continue
            elif snp[3] == 'T' and snp[4] == 'A':
                continue
            elif snp[3] == 'C' and snp[4] == 'G':
                continue
            elif snp[3] == 'G' and snp[4] == 'C':
                continue
            else:
                AFtable[snp[2]] = snp
            lines += 1

    # read the map file once so that we have SNP rsid for each SNP
    print("reading mapfile", mapfile, file=sys.stderr)
    try:
        File = open(mapfile, 'r')
    except IOError:
        print("map file", mapfile, "not found", file=sys.stderr)
        return -1
    SNPlist = []
    for line in File:
        snp = line.strip().split()
        if snp[0][0] == '#':
            continue  # BUG if first line has '# chrom', make sure the lines match up to ped file... 10/04/2014
        SNPlist.append([snp[1], snp[0], snp[2]])
    File.close()
    #######

    # read the ped file line by line and calculate ancestry for each individual
    File = open(pedfile, 'r')
    samples = 0
    for line in File:
        if CORES > 1 and samples % CORES != modulo_CORE:
            samples += 1
            continue
        GLLs = {}
        genotypes = line.strip().split()
        sampleid = genotypes[1]
        outfilename = outfile+'.'+sampleid + '.input'
        File1 = open(outfilename, 'w')

        # Family ID      Individual ID   FatherID MotherID Sex  phenotype #SIM_1 SIM_1 0 0 1 -9 C C C C A

        print('#GLL', end=' ', file=File1)
        for i in range(5, len(headerline)):
            print(headerline[i], end=' ', file=File1)
        print('\n', end=' ', file=File1)

        overlapping_markers = 0
        missing_markers = 0
        for i in range(6, len(genotypes), 2):
            s = i-6
            s //= 2
            try:
                snp = AFtable[SNPlist[s][0]]
                flag = 0
                if genotypes[i] == snp[3] and genotypes[i+1] == snp[3]:
                    print(gll_high, gll_low, gll_low, end=' ', file=File1)
                elif genotypes[i] == snp[4] and genotypes[i+1] == snp[4]:
                    print(gll_low, gll_low, gll_high, end=' ', file=File1)
                elif genotypes[i] == snp[3] and genotypes[i+1] == snp[4]:
                    print(gll_low, gll_high, gll_low, end=' ', file=File1)
                elif genotypes[i+1] == snp[3] and genotypes[i] == snp[4]:
                    print(gll_low, gll_high, gll_low, end=' ', file=File1)
                elif genotypes[i+1] == RC[snp[3]] and genotypes[i] == RC[snp[4]]:
                    print(gll_low, gll_high, gll_low, end=' ', file=File1)
                elif genotypes[i] == RC[snp[3]] and genotypes[i+1] == RC[snp[4]]:
                    print(gll_low, gll_high, gll_low, end=' ', file=File1)
                elif genotypes[i] == RC[snp[4]] and genotypes[i+1] == RC[snp[4]]:
                    print(gll_low, gll_low, gll_high, end=' ', file=File1)
                elif genotypes[i] == RC[snp[3]] and genotypes[i+1] == RC[snp[3]]:
                    print(gll_high, gll_low, gll_low, end=' ', file=File1)
                else:
                    flag = 1
                    missing_markers += 1
                    # print i,'missing',snp,SNPlist[s],genotypes[i],genotypes[i+1]
                if flag == 0:
                    for j in range(5, len(snp)):
                        print(snp[j], end=' ', file=File1)
                    print('\n', end=' ', file=File1)
                    overlapping_markers += 1
            except KeyError:
                # print >>sys.stderr, 'variant not found',SNPlist[s][0];
                pass

        File1.close()
        print("\n\nancestry admixture calculations for individual:", samples+1,
              sampleid, 'using', overlapping_markers, 'markers', file=sys.stderr)
        # print >>sys.stderr, "missing",missing_markers,len(genotypes)/2;

        if WINDOW == 0:
            call([f"\"{DIRNAME!s}/ANCESTRY\" -p 2 --pr {PARSIMONY!r} --HWE " +
                 repr(HWE_CHECK) + " -i " + outfilename + " > " + outfilename + ".ancestry"], shell=True)
        else:
            print("calling BFGS method in windows", file=sys.stderr)
            call([f"time \"{DIRNAME!s}/ANCESTRY.1\" -p 2 --pr {PARSIMONY!r} " +
                  "-i " + outfilename + " > " + outfilename + ".ancestry"], shell=True)
        # call(["rm -f " + outfilename],shell=True);
        # call(["rm -f " + outfilename +  " " + outfilename + ".ancestry"],shell=True);
        samples += 1
        # sys.exit();
        """
        """
    File.close()

####################################################################################################################################################


# genotype file is rsid AG pairs
def make_ancestry_inputfile_rsid(afmatrixfile, GLLfile, outfile):

    RC = {}
    RC['A'] = 'T'
    RC['T'] = 'A'
    RC['C'] = 'G'
    RC['G'] = 'C'
    # gll_low = -3.602; gll_high = -0.000021715; gll_equal = -0.477121255; ## assuming error rate = 1/20000 or 50 errors per million genotypes
    gll_low = -4
    gll_high = -0.000086868
    gll_equal = -0.477121
    # gll_low = -9; gll_high = 0.00000; gll_equal = -0.477121;
    column = 1

    # read genotype file and store genotypes in hashtable
    File1 = open(outfile, 'w')
    GLLs = {}
    encod = get_file_encoding(GLLfile)
    with open(GLLfile, encoding=encod) as File:
        for line in File:
            snp = line.rstrip().split()
            if snp[column] == 'NN' or snp[column] == 'N/N' or snp[column] == 'NoCall' or snp[column] == '--':
                continue
            if snp[column][1] == '/':
                genotype = snp[column][0] + snp[column][2]
            else:
                genotype = snp[column]
            GLLs[snp[0]] = genotype

    lines = 0
    valid_snps = 0

    encod = get_file_encoding(afmatrixfile)
    with open(afmatrixfile) as File:
        for line in File:
            snp = line.rstrip().split()
            try:
                snp[3], snp[4]
            except IndexError:
                print(snp)
                raise
            if line[0] == '#' or lines == 0:
                print('#GLL', end=' ', file=File1)
                for i in range(5, len(snp)):
                    print(snp[i], end=' ', file=File1)
                print('\n', end=' ', file=File1)
            # ignore A/T and C/G snps that are strand ambiguous
            elif INCLUDE_STRAND_AMBIG == 0 and snp[3] == 'A' and snp[4] == 'T':
                continue
            elif INCLUDE_STRAND_AMBIG == 0 and snp[3] == 'T' and snp[4] == 'A':
                continue
            elif INCLUDE_STRAND_AMBIG == 0 and snp[3] == 'C' and snp[4] == 'G':
                continue
            elif INCLUDE_STRAND_AMBIG == 0 and snp[3] == 'G' and snp[4] == 'C':
                continue
            # the two alleles are [ACTG], ignore indel D/I alleles
            elif snp[3] in RC and snp[4] in RC:
                try:
                    GENOTYPE = GLLs[snp[2]]
                    flag = 0
                    if GENOTYPE[0] == snp[3] and GENOTYPE[1] == snp[3]:
                        print(gll_high, gll_low, gll_low, end=' ', file=File1)
                    elif GENOTYPE[0] == snp[4] and GENOTYPE[1] == snp[4]:
                        print(gll_low, gll_low, gll_high, end=' ', file=File1)
                    elif GENOTYPE[0] == snp[3] and GENOTYPE[1] == snp[4]:
                        print(gll_low, gll_high, gll_low, end=' ', file=File1)
                    elif GENOTYPE[1] == snp[3] and GENOTYPE[0] == snp[4]:
                        print(gll_low, gll_high, gll_low, end=' ', file=File1)
                    elif GENOTYPE[1] == RC[snp[3]] and GENOTYPE[0] == RC[snp[4]]:
                        print(gll_low, gll_high, gll_low, end=' ', file=File1)
                    elif GENOTYPE[0] == RC[snp[3]] and GENOTYPE[1] == RC[snp[4]]:
                        print(gll_low, gll_high, gll_low, end=' ', file=File1)
                    elif GENOTYPE[0] == RC[snp[4]] and GENOTYPE[1] == RC[snp[4]]:
                        print(gll_low, gll_low, gll_high, end=' ', file=File1)
                    elif GENOTYPE[0] == RC[snp[3]] and GENOTYPE[1] == RC[snp[3]]:
                        print(gll_high, gll_low, gll_low, end=' ', file=File1)
                    else:
                        # print >>File1, gll_equal,gll_equal,gll_equal,
                        # print >>sys.stderr, GENOTYPE,snp;
                        flag = 1
                    if flag == 0:
                        for i in range(5, len(snp)):
                            print(snp[i], end=' ', file=File1)
                        print('\n', end=' ', file=File1)
                        valid_snps += 1
                        # print >>sys.stderr, snp[2];
                except KeyError:
                    pass
            lines += 1
    # print("ancestry genotype file has", valid_snps, "markers", file=sys.stderr)

    File1.close()


def ancestry_pipeline_genotypes(matrixfile: Union[PathLike[str], str], genotype_file, outfile):
    # print("making input file for ancestry calculations", file=sys.stderr)
    make_ancestry_inputfile_rsid(
        matrixfile, genotype_file, genotype_file + ".ancestry.input")
    # print("ancestry admixture calculations", genotype_file, file=sys.stderr)
    if WINDOW == 0:
        out = call([f"time '{DIRNAME!s}/ANCESTRY' --LRT 7.68 -p 2 --pr 1 -i " +
                    f"'{genotype_file}.ancestry.input' >  '{outfile}'"], shell=True, stdout=DEVNULL, stderr=DEVNULL)
    else:
        print("calling BFGS method in windows", file=sys.stderr)
        out = call([f"time \"{DIRNAME!s}/ANCESTRY.1\" -p 2 --pr 1 -i " +
                    f"\"{genotype_file}.ancestry.input\" > \"{outfile}\""])
    # call(["time " + executable_dir_path + "/ANCESTRY -p 2 --pr 1 -i " + genotype_file + ".ancestry.input > " + outfile],shell=True);
    # call(["time ./ANCESTRY --ea 1 -p 2 --pr 1 -i " + bamfile + ".ancestry.input > " + outfile],shell=True);
    # call(["rm -f " + bamfile + ".ancestry.input"],shell=True);


###########################################################################################################################


parser = OptionParser()
# parser.add_option("--type",dest="type",type="string",help="type of input file: bam/geno/genotype/VCF/ped",default="");
parser.add_option("-b", "--bam", dest="BAMfile", type="string",
                  help="input bam file name", default="")
parser.add_option("--geno", dest="genofile", type="string",
                  help="input file name for simple genotype file: rsid AG on each line", default="")
parser.add_option("--vcf", "--VCF", dest="vcffile", type="string",
                  help="input file name for VCF file", default="")
parser.add_option("--plink", dest="pedfile", type="string",
                  help="input file name for plink file in ped/map format (provide prefix of ped file name without the .ped)", default="")
parser.add_option("-f", "--freq", dest="AFfile", type="string",
                  help="allele frequency file", default="")
parser.add_option("-o", "--out", dest="outfile", type="string",
                  help="prefix of output file names", default="iADMIX.out")
parser.add_option("-p", "--poolsize", dest="POOLSIZE", type="int",
                  help="pool size for non-diploid samples, default = 2", default=2)
parser.add_option("--pr", "--parsimony", dest="pr", type="int",
                  help="parsimonious ancestry estimation 0/1, default = 1", default=1)
parser.add_option("--hwe", "--HWE", dest="HWE", type="int",
                  help="HWE deviation estimation 0/1, default = 0", default=0)
parser.add_option("--strand", dest="include_strand_ambiguous", type="int",
                  help="<0/1> include SNPs that are strand ambiguous (A/T, C/G): default = 0 (not included)", default=0)
# parser.add_option("--path", dest="path", type="string",
#                   help="path to directory with executables ANCESTRY and calculateGLL", default=".")
parser.add_option("-c", "--cores", dest="CORES", type="int",
                  help="number of cores for processing in parallel", default=1)
parser.add_option("-m", "--modcore", dest="MOD_CORE",
                  type="int", help="0/1/2/3..../CORES-1", default=0)
parser.add_option("-w", "--windows", dest="WINDOW",
                  type="int", help="0/1", default=0)
# only output genotypes likelihoods, no ancestry calculation...
parser.add_option("-g", "--gll", dest="GLL_ONLY",
                  type="int", help="0/1", default=0)
parser.add_option("--addchr", dest="ADD_chr_names",
                  action="store_true", default=False)
(options, args) = parser.parse_args()
POOLSIZE = options.POOLSIZE
CORES = options.CORES
HWE_CHECK = options.HWE
OUTPUT_GLL_ONLY = options.GLL_ONLY
WINDOW = options.WINDOW
ADD_chr_names = options.ADD_chr_names
INCLUDE_STRAND_AMBIG = options.include_strand_ambiguous
PARSIMONY = options.pr

# file = open(sys.argv[1],'rb'); while True: byte = file.read(1); print ord(byte); sys.exit();

if __name__ == "__main__":
    if (options.BAMfile == "" and options.genofile == "" and options.pedfile == "") or options.AFfile == "":

        print("\n####program iAdmix for estimation of ancestry using population allele frequencies###\n", file=sys.stderr)
        print("python runancestry.py --bam sample.bam --geno sample.genotypes --freq allelefrequency_file --out outputfile --poolsize 2", file=sys.stderr)
        print("\nExample: python runancestry.py --freq hapmap3.allchroms.shared.matrix --geno HGDP01254.genotypes --out HGDP12054.ancestry", file=sys.stderr)
        print("Example: python runancestry.py --freq hapmap3.allchroms.shared.matrix --bam HGDP01254.sorted.bam --out HGDP12054.ancestry", file=sys.stderr)
        print("Example: python runancestry.py --freq hapmap3.allchroms.shared.matrix --plink HGDP01254.genotypes (ped/map) --out HGDP12054.ancestry\n\n", file=sys.stderr)
        parser.print_help()
        sys.exit()
    elif options.pedfile != "":
        make_ancestry_inputfile_plink(options.pedfile+'.ped', options.pedfile +
                                      '.map', options.AFfile, options.outfile, options.MOD_CORE)

    elif options.BAMfile != "":
        if not (DIRNAME / "ANCESTRY").is_file():
            print("ANCESTRY executable does not exist\n", file=sys.stderr)
            sys.exit()
        if not (DIRNAME / "calculateGLL").is_file():
            print("calculateGLL executable does not exist\n", file=sys.stderr)
            sys.exit()
        ancestry_pipeline_bam(options.AFfile, options.BAMfile,
                              options.outfile)
    elif options.genofile != "":
        if not (DIRNAME / "ANCESTRY").is_file():
            print("ANCESTRY executable does not exist\n", file=sys.stderr)
            sys.exit()
        ancestry_pipeline_genotypes(
            options.AFfile, options.genofile, options.outfile)


########################################################################################################################################
