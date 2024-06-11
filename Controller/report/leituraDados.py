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
	endereco = os.path.join("../Controller", "DataFiles", "Files", "SNPs.txt")
	snp = []
	with open(endereco, "r") as file:
		reading = file.readlines()
		for line in reading:
			line = line.replace("\n", "")
			snp.append(line)
	endereco = os.path.join("../Controller", "DataFiles", "Brutos", f"{ID}.txt")
	if os.path.exists(endereco) == False:
		endereco = os.path.join("../Controller", "DataFiles", "Brutos", f"{ID}_23andMe.txt")
	if os.path.exists(endereco) == False:
		print(f"Dado Bruto {ID} não encontrado.")
		return -1
	with open(endereco, "r") as file:
		reading = file.readlines()
		for i in range(len(snp)):
			check = 0
			sp = snp[i].split("\t")
			for line in reading:
				line = line.replace("\n", "")
				bruto = line.split("\t")
				if sp[0] == bruto[0]:
					snp[i] = sp[0]+"\t"+bruto[3]+"\t"+sp[2]+"\t"+sp[3]
					check = 1
					break
			if check == 0:
				snp[i] = sp[0]+"\t"+"--"+"\t"+sp[2]+"\t"+sp[3]
	print("Dados brutos lidos\n")
	return snp