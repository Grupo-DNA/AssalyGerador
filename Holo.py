from laudos import *
from pypdf import PdfWriter

with open("Gerar.txt", "r") as handle:
	reading = handle.readlines()
	for line in reading:
		if line[0] == "#" or line == "\n":
			continue
		line = line.replace("\n", "")
		ID = line
		outpdf = PdfWriter()
		snp = read_SNPs(ID)
		holobionte(outpdf, snp)
		#rotas_resultados(outpdf)
		laudo_salvar(outpdf, f"Test{ID}")