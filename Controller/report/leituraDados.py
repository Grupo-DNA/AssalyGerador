import os
import re
from tqdm import tqdm
from google.oauth2 import service_account
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload

# Arquivo SNPs.txt vai ser montado dentro do programa com os dados de cada paciente
def read_SNPs(ID):
	print("Lendo dados brutos...\n")
	endereco = os.path.join("../Controller", "DataFiles", "Files", "SNPs.txt")
	print(endereco)
	snp = []
	with open(endereco, "r") as file:
		reading = file.readlines()
		for line in reading:
			line = line.replace("\n", "")
			snp.append(line)
	endereco = os.path.join("../Controller", "DataFiles", "Brutos", f"{ID}.txt")
	if os.path.exists(endereco) == False:
		endereco = os.path.join("Brutos", f"{ID}_23andMe.txt")
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
