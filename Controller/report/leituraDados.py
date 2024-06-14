import os
import re
from tqdm import tqdm
from Controller.controller.logger import get_logger
DNA_REGEX = re.compile(r"[ACTGD-]{2}")

def get_file_encoding(filepath):
    """Detects the encoding of the given file."""
    for encoding in ('utf-8', 'latin-1'):
        try:
            with open(filepath, "r", encoding=encoding) as file:
                file.read()
            return encoding
        except UnicodeDecodeError:
            continue
    raise Exception("Unable to determine file encoding.")

def get_file_delimiter(filepath, encoding=None):
    """Detects the delimiter used in the given file."""
    if encoding is None:
        encoding = get_file_encoding(filepath)
    with open(filepath, "r", encoding=encoding) as file:
        first_line = file.readline().strip()
        if "\t" in first_line:
            return "\t"
        elif "," in first_line:
            return ","
        elif ";" in first_line:
            return ";"
    raise Exception("Unable to determine file delimiter.")

def get_callback_fromline(line, delimiter):
    """Returns a callback function to extract genotype data based on the apparent format."""
    def callback1(value, delim):
        rsid, _, (a1, a2) = value.rstrip("\n").split(delim)
        if rsid == ".":
            return None, None, None
        return rsid, a1, a2

    def callback2(value, delim):
        rsid, chrom, _, alleles = value.rstrip("\n").split(delim)
        a1 = alleles[0] if len(alleles) > 0 else None
        a2 = alleles[1] if len(alleles) == 2 else None
        if rsid == "." or len(alleles) > 2:
            return None, None, None
        if a2 is None and chrom in ["MT", "X", "Y"]:
            a2 = "-"
        return rsid, a1, a2

    def callback3(value, delim):
        rsid, _, a1, a2, _ = value.rstrip("\n").split(delim)
        if rsid == ".":
            return None, None, None
        return rsid, a1, a2

    def callback4(value, delim):
        rsid, _, _, a1, a2, _ = value.rstrip("\n").split(delim)
        if rsid == ".":
            return None, None, None
        return rsid, a1, a2

    splitted = line.rstrip("\n").split(delimiter)
    if len(splitted) == 3 and len(splitted[-1]) == 2:
        return callback1
    elif len(splitted) == 4 and len(splitted[-1]) == 2:
        return callback2
    elif len(splitted) == 5 and len(splitted[2]) == 1 and len(splitted[3]) == 1:
        return callback3
    elif len(splitted) == 6 and len(splitted[3]) == 1 and len(splitted[4]) == 1:
        return callback4
    else:
        raise Exception("Unexpected data format")

def read_SNPs(ID):
    print("Lendo dados brutos...\n")
    snp_path = os.path.join("../Controller", "DataFiles", "Files", "SNPs.txt")
    snp = []

    # Read SNPs file
    with open(snp_path, "r") as file:
        snp = [line.strip() for line in file.readlines()]

    # Identify the correct raw data file
    raw_data_path = os.path.join("../Controller", "DataFiles", "Brutos", f"{ID}.txt")
    if not os.path.exists(raw_data_path):
        raw_data_path = os.path.join("../Controller", "DataFiles", "Brutos", f"{ID}_23andMe.txt")
    if not os.path.exists(raw_data_path):
        print(f"Dado Bruto {ID} não encontrado.")
        return -1

    # Get file encoding and delimiter
    encoding = get_file_encoding(raw_data_path)
    delimiter = get_file_delimiter(raw_data_path, encoding)

    # Initialize the progress bar
    total_lines = sum(1 for _ in open(raw_data_path, encoding=encoding))
    progress_bar = tqdm(total=total_lines, desc="Processing raw data")

    alelos_dict = {}

    with open(raw_data_path, encoding=encoding) as file:
        # Skip header
        for line in file:
            if re.match(r'rs\d+', line):
                callback = get_callback_fromline(line, delimiter)
                break
            progress_bar.update(1)

        # Read and process the file line by line
        for line in file:
            if not line.strip() or line.startswith("#") or not re.match(r'rs\d+', line):
                continue
            rsid, a1, a2 = callback(line, delimiter)
            if a1 not in ["A", "T", "C", "G", "D", "I"] or a2 not in ["A", "T", "C", "G", "D", "I"]:
                a1 = "-"
                a2 = "-"
            alelos_dict[rsid] = f"{a1}{a2}"
            progress_bar.update(1)

    progress_bar.close()

    # Update SNP list with data from raw file
    for i in range(len(snp)):
        sp = snp[i].split("\t")
        rsid = sp[0]
        if rsid in alelos_dict:
            snp[i] = f"{rsid}\t{alelos_dict[rsid]}\t{sp[2]}\t{sp[3]}"
        else:
            snp[i] = f"{rsid}\t--\t{sp[2]}\t{sp[3]}"

    print("Dados brutos lidos\n")
    return snp
