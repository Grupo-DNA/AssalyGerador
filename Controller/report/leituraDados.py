import os
import re
from tqdm import tqdm
from Controller.controller.auxiliar_functions import get_logger

DNA_REGEX = re.compile(r"[ACTGD-]{2}")

def get_callback_fromline(line: str, delimiter: str):
    """Retorna uma função de extração dos dados genotípicos,
    baseada na formatação aparente dos mesmos."""
    logger = get_logger()

    def callback1(value: str, delim: str):
        # rsid, chrom, pos, genotype
        parts = value.rstrip("\n").split(delim)
        if len(parts) < 4:
            return None, None, None
        rsid, chrom, pos, alleles = parts
        if rsid == "." or len(alleles) != 2:
            return None, None, None
        a1 = alleles[0]
        a2 = alleles[1]
        return rsid, a1, a2

    splitted = line.rstrip("\n").split(delimiter)

    # Assume the format is: rsid chrom pos genotype
    if len(splitted) == 4:
        return callback1
    else:
        logger.error('Erro tipo 3 (dados brutos/normscore)\n\tDado bruto em formato inadequado')
        raise Exception()

def get_file_encoding(filepath: str):
    """Chooses a possible encoding for .csv files
    Only between UTF-8 and Latin-1 (ISO 8859-1)
    Args:
        filepath (str): .csv filepath.
    Returns:
        str: Encoding found.
    """
    logger = get_logger()
    for encoding in ('utf-8', 'latin-1'):
        try:
            with open(filepath, "r", newline="", encoding=encoding) as handle:
                handle.read()
            return encoding
        except UnicodeDecodeError:
            pass
    logger.error('Erro tipo 3 (dados brutos/normscore)\n\tDado bruto em formato inadequado (incapaz de determinar o encoding)')
    raise Exception()

def get_file_delimiter(filepath: str, encoding: str = None):
    """Guess .csv file delimiter.
    Args:
        filepath (PathLike): .csv filepath.
        encoding (str, Optional): .csv file encoding. Defaults to None.
    Returns:
        str: Delimiter chosen between "\\t", "," and ";".
    """
    if encoding is None:
        encoding = get_file_encoding(filepath)
    delims = ["\t", ",", ";"]
    with open(filepath, "r", newline="", encoding=encoding) as handle:
        sniffer = handle.readline().rstrip("\n")
        by_order = sorted(
            delims, key=sniffer.count,
            reverse=True)
        return by_order[0]

def read_SNPs(ID):
    print("Lendo dados brutos...\n")
    
    # Endereço do arquivo SNPs.txt
    snps_file_path = os.path.join("../Controller", "DataFiles", "Files", "SNPs.txt")
    print(snps_file_path)
    
    # Lendo SNPs
    snp_list = []
    with open(snps_file_path, "r") as file:
        for line in file:
            line = line.strip()
            snp_list.append(line)
    
    # Determinando o caminho do arquivo bruto
    raw_data_path = os.path.join("../Controller", "DataFiles", "Brutos", f"{ID}.txt")
    if not os.path.exists(raw_data_path):
        raw_data_path = os.path.join("Brutos", f"{ID}_23andMe.txt")
    if not os.path.exists(raw_data_path):
        print(f"Dado Bruto {ID} não encontrado.")
        return -1
    
    # Lendo o arquivo bruto
    encoding = get_file_encoding(raw_data_path)
    delimiter = get_file_delimiter(raw_data_path, encoding=encoding)
    with open(raw_data_path, "r", encoding=encoding) as file:
        raw_data_lines = file.readlines()
    
    # Determinando o callback apropriado
    callback = get_callback_fromline(raw_data_lines[0], delimiter)
    
    # Criando um dicionário de genótipos a partir do arquivo bruto
    genotype_dict = {}
    for line in raw_data_lines:
        if not line.strip() or line.startswith("#") or not re.match(r'rs\d+', line):
            continue
        rsid, a1, a2 = callback(line, delimiter)
        if rsid:
            genotype_dict[rsid] = f"{a1}{a2}"
    
    # Atualizando a lista de SNPs com genótipos lidos do arquivo bruto
    for i in range(len(snp_list)):
        snp_fields = snp_list[i].split("\t")
        rsid = snp_fields[0]
        if rsid in genotype_dict:
            genotype = genotype_dict[rsid]
            snp_list[i] = f"{rsid}\t{genotype}\t{snp_fields[2]}\t{snp_fields[3]}"
        else:
            snp_list[i] = f"{rsid}\t--\t{snp_fields[2]}\t{snp_fields[3]}"
    
    print("Dados brutos lidos\n")
    return snp_list
