from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from reportlab.platypus import Paragraph
from io import BytesIO
from Controller.controller.config import CONFIG
from pypdf import PageObject, PdfReader
from datetime import date
import os
import sys
import math
from ancestralidade import AncestralidadeMaker
from base import Paciente
from reportlab.lib.colors import HexColor
from Controller.controller.auxiliar_functions import *

def style(canv, name):
	style = CONFIG.styles[name]
	canv.setFont(style.fontName, style.fontSize)
	canv.setFillColor(style.textColor)

def dicio_descri():
	endereco = os.path.join("Files", "Descri.txt")
	dicio = {}
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			sp = line.replace("\n", "").split("\t")
			dicio[sp[0]] = sp[1]
	return dicio

def switch_highest(high):
	if high == "europeus":
		return "européia"
	if high == "africanos":
		return "africana"
	if high == "norte-americanos":
		return "norte-americana"
	if high == "sul-americanos":
		return "sul-americana"
	if high == "asiáticos":
		return "asiática"
	if high == "oceânicos":
		return "oceânica"

# roundown arredonda para casas de 10
def roundown(number):
	if math.floor(number/10)*10 > 91:
		return 99
	else:
		return math.floor(number/10)*10

def switch_column(col):
	if col == 0:
		return 36
	elif col == 1:
		return 222
	elif col == 2:
		return 408
	else:
		return 0

def draw_box_rota(col,y,w,h,h2,name,canv,snp):
	if "2" in name or "1" in name:
		name = name[0:len(name)-1]
	x = switch_column(col)-9
	canv.setStrokeColorRGB(0.56,0.75,0.82)
	canv.rect(x,y,w,h, stroke=1, fill=0)
	if "baixo risco" in name.lower():
		num = 1
	elif "médio risco" in name.lower():
		num = 2
	elif "alto risco" in name.lower():
		num = 3
	elif "sem efeito" in name.lower():
		num = 4
	else:
		num = find_color(name,snp)
	#Green
	if num == 1:
		rect_style = CONFIG.styles["ancestralidade.rect-green"]
	#Yellow
	elif num == 2:
		rect_style = CONFIG.styles["ancestralidade.rect-yellow"]
	#Red
	elif num == 3:
		rect_style = CONFIG.styles["ancestralidade.rect-red"]
	#Gray
	elif num == 4:
		rect_style = CONFIG.styles["ancestralidade.rect-grey"]
	canv.setFillColor(rect_style.textColor)
	canv.rect(x,y+h-h2,w,h2, stroke=0, fill=1)
	canv.setFillColorRGB(255,255,255)
	canv.setFont(rect_style.fontName, rect_style.fontSize)
	if len(name) > 25:
		canv.setFont(rect_style.fontName, 8)
		name = name.split(" ", 1)
		canv.drawCentredString((w/2)+x, y+h-h2+15, name[0].upper())
		canv.drawCentredString((w/2)+x, y+h-h2+5, name[1].upper())
	else:
		canv.drawCentredString((w/2)+x, y+h-h2+8, name.upper())

def count_genes_size(name):
	print("name misterioso", name)
	ini = 0
	if "1" in name:
		name = name[0:len(name)-1]
	elif "2" in name:
		name = name[0:len(name)-1]
		ini = 22
	count = 0
	check = 0
	endereco = os.path.join("Files", "Lista.txt")
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			sep = line.replace("\n", "").split("\t")
			if check == 1:
				if sep[0] == "":
					break
				elif sep[0][0] == 'r' and sep[0][1] == 's':
					count += 1
			elif sep[0] == name:
				check = 1
	if count > 22 and ini == 0:
		return (22*30)
	return (count*30)-(ini*30)

def draw_rs_rotas(outpdf, snp, canv, name, column, posy):
	posy += count_genes_size(name)-28
	count = 0
	if "1" in name:
		name = name[0:len(name)-1]
	elif "2" in name:
		name = name[0:len(name)-1]
		count = 22
	posx = switch_column(column)
	alto = []
	medio = []
	baixo = []
	semef = []
	check = 0
	endereco = os.path.join("Files", "Lista.txt")
	end2 = os.path.join("Files", "Efeitos.txt")
	conjunto_genes_unicos = set()
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			sep = line.replace("\n", "").split("\t")
			if check > 0:
				if sep[0] == "":
					break
				elif sep[0][0] == 'r' and sep[0][1] == 's':
					for i in range(len(snp)):
						seu_gene = snp[i].split("\t")
						if sep[0] == seu_gene[0] and seu_gene[0] not in conjunto_genes_unicos:
							conjunto_genes_unicos.add(seu_gene[0])
							if seu_gene[1] == "--":
								semef.append(f"{seu_gene[0]}\t{seu_gene[1]}\t{seu_gene[3]}")
							else:
								with open(end2, "r") as f:
									reading = f.readlines()
									for l in reading:
										sep2 = l.replace("\n", "").split("\t")
										if sep2[0] == name and sep2[2] == sep[0] and sep2[3] == seu_gene[1]:
											if sep2[4] == "0":
												baixo.append(f"{seu_gene[0]}\t{seu_gene[1]}\t{seu_gene[3]}")
											elif sep2[4] == "1":
												medio.append(f"{seu_gene[0]}\t{seu_gene[1]}\t{seu_gene[3]}")
											elif sep2[4] == "2":
												alto.append(f"{seu_gene[0]}\t{seu_gene[1]}\t{seu_gene[3]}")
											break
			elif sep[0] == name:
				check = 1
	style = CONFIG.styles["ancestralidade.gene"]
	gene_style = CONFIG.styles["ancestralidade.text-red"]
	counter = 0
	for i in range(len(alto)):
		if count <= 0 and counter < 22:
			sp = alto[i].split("\t")
			canv.setFont(style.fontName, style.fontSize)
			canv.setFillColor(style.textColor)
			canv.drawString(posx, posy, sp[0])
			canv.setFont(gene_style.fontName, gene_style.fontSize)
			canv.setFillColor(gene_style.textColor)
			canv.drawString(posx, posy+10, sp[2])
			canv.drawString(posx+130, posy+10, sp[1])
			posy -= 30
			counter += 1
		count -= 1
	gene_style = CONFIG.styles["ancestralidade.text-yellow"]
	for i in range(len(medio)):
		if count <= 0 and counter < 22:
			sp = medio[i].split("\t")
			canv.setFont(style.fontName, style.fontSize)
			canv.setFillColor(style.textColor)
			canv.drawString(posx, posy, sp[0])
			canv.setFont(gene_style.fontName, gene_style.fontSize)
			canv.setFillColor(gene_style.textColor)
			canv.drawString(posx, posy+10, sp[2])
			canv.drawString(posx+130, posy+10, sp[1])
			posy -= 30
			counter += 1
		count -= 1
	gene_style = CONFIG.styles["ancestralidade.text-green"]
	for i in range(len(baixo)):
		if count <= 0 and counter < 22:
			sp = baixo[i].split("\t")
			canv.setFont(style.fontName, style.fontSize)
			canv.setFillColor(style.textColor)
			canv.drawString(posx, posy, sp[0])
			canv.setFont(gene_style.fontName, gene_style.fontSize)
			canv.setFillColor(gene_style.textColor)
			canv.drawString(posx, posy+10, sp[2])
			canv.drawString(posx+130, posy+10, sp[1])
			posy -= 30
			counter += 1
		count -= 1
	gene_style = CONFIG.styles["ancestralidade.text-grey"]
	for i in range(len(semef)):
		if count <= 0 and counter < 22:
			sp = semef[i].split("\t")
			canv.setFont(style.fontName, style.fontSize)
			canv.setFillColor(style.textColor)
			canv.drawString(posx, posy, sp[0])
			canv.setFont(gene_style.fontName, gene_style.fontSize)
			canv.setFillColor(gene_style.textColor)
			canv.drawString(posx, posy+10, sp[2])
			canv.drawString(posx+130, posy+10, sp[1])
			posy -= 30
			counter += 1
		count -= 1

# calculo PRS rotas
def find_color(name, snp):
    lista_file_path = os.path.join("Files", "Lista.txt")
    efeitos_file_path = os.path.join("Files", "Efeitos.txt")
    
    # Mapeando o nome "Glicação" para "Insulina e Glicose", se aplicável
    if name == "Glicação":
        name = "Insulina e Glicose"
    
    found_name = False
    calc = 0
    maxval = 0
    
    # Abrindo o arquivo Lista.txt e iterando sobre suas linhas
    with open(lista_file_path, "r") as lista_file:
        for line in lista_file:
            fields = line.strip().split("\t")
            # Verificando se encontramos o nome desejado
            if not found_name:
                if fields[0] == name:
                    found_name = True
                continue
            if fields[0] == "":
                break  # Interrompe o loop se encontrar uma linha vazia
            if fields[0][0] == 'r' and fields[0][1] == 's':
                for snp_entry in snp:
                    snp_fields = snp_entry.split("\t")
                    if fields[0] == snp_fields[0] and snp_fields[1] != "--":
                        with open(efeitos_file_path, "r") as efeitos_file:
                            for efeito_line in efeitos_file:
                                efeito_fields = efeito_line.strip().split("\t")
                                if efeito_fields[0] == name and efeito_fields[2] == fields[0] and efeito_fields[3] == snp_fields[1]:
                                    calc += int(efeito_fields[4])
                                    maxval += 2
    if maxval == 0:
        return 4
    color = calc / maxval
    if color < 0.3:
        return 1  # Verde
    elif color < 0.7:
        return 2  # Amarelo
    else:
        return 3  # Vermelho
	
def read_SNPs(service, ID):
    print("Lendo dados brutos...\n")
    print("Procurando arquivo no drive...")

    file_id = get_file_id(service, ID)
    if not file_id:
        print(f"Dado Bruto {ID} não encontrado no Drive.")
        return -1

    filename = f"bruto_{ID}.txt"
    download_file(service, file_id, filename)

    snp_file_path = os.path.join("Files", "SNPs.txt")
    snps = []
    with open(snp_file_path, "r") as file:
        snps = [line.strip() for line in file.readlines()]

    if not os.path.exists(filename):
        print(f"Dado Bruto {ID} não encontrado.")
        return -1

    delimiter = get_file_delimiter(filename)
    total_lines = sum(1 for _ in open(filename, encoding="utf-8"))
    progress_bar = tqdm(total=total_lines)
    snp_dict = {line.split("\t")[0]: line.split("\t")[1:] for line in snps}

    with open(filename, encoding="utf-8") as file:
        for line in file:
            if re.match(r'rs\d+', line):
                break
            progress_bar.update(1)

        callback = get_callback_fromline(line, delimiter)
        for line in file:
            if not line.strip() or line.startswith("#") or not re.match(r'rs\d+', line):
                continue
            rsid, a1, a2 = callback(line, delimiter)
            if rsid in snp_dict:
                snp_dict[rsid] = [a1, a2] + snp_dict[rsid][1:]
            progress_bar.update(1)
    
    progress_bar.close()
    os.remove(filename)

    snps = [f"{rsid}\t{data[0]}\t{data[1]}\t{data[2]}" if data[0] != "--" else f"{rsid}\t--\t{data[1]}\t{data[2]}" for rsid, data in snp_dict.items()]

    print("Dados brutos lidos\n")
    return snps

def make_rect_color(name, canv, posx, posy, snp, width):
	num=find_color(name,snp)
	#Green
	if num == 1:
		rect_style = CONFIG.styles["ancestralidade.rect-green"]
		canv.setFillColor(rect_style.textColor)
		if not width > 0:
			canv.rect(posx,posy-7,283,23, stroke=0, fill=1)
		else:
			canv.rect(posx,posy-7,width,23, stroke=0, fill=1)
	#Yellow
	if num == 2:
		rect_style = CONFIG.styles["ancestralidade.rect-yellow"]
		canv.setFillColor(rect_style.textColor)
		if not width > 0:
			canv.rect(posx,posy-7,414,23, stroke=0, fill=1)
		else:
			canv.rect(posx,posy-7,width,23, stroke=0, fill=1)
	#Red
	if num == 3:
		rect_style = CONFIG.styles["ancestralidade.rect-red"]
		canv.setFillColor(rect_style.textColor)
		if not width > 0:
			canv.rect(posx,posy-7,544,23, stroke=0, fill=1)
		else:
			canv.rect(posx,posy-7,width,23, stroke=0, fill=1)
	#Gray
	if num == 4:
		rect_style = CONFIG.styles["ancestralidade.rect-grey"]
		canv.setFillColor(rect_style.textColor)
		if not width > 0:
			canv.rect(posx,posy-7,283,23, stroke=0, fill=1)
		else:
			canv.rect(posx,posy-7,width,23, stroke=0, fill=1)
	style = CONFIG.styles["ancestralidade.visao-geral"]
	canv.setFont(style.fontName, style.fontSize)

def laudo_capa(outpdf, name):
	endereco = os.path.join(CONFIG.template, "capa.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	style = CONFIG.styles["ancestralidade.capa-nome"]
	c.setFillColor(style.textColor)
	c.setFont(style.fontName, style.fontSize)
	if len(name) > 35:
		c.setFont(style.fontName, style.fontSize-3)
		c.drawString(48, 128, name)
	else:
		c.drawString(48, 125, name)
	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)
	endereco = os.path.join(CONFIG.template, "sumario.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
	outpdf.add_page(pagina)
	endereco = os.path.join(CONFIG.template, "nossa_filosofia.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
	outpdf.add_page(pagina)
	print("Capa gerada\n")

def ancestralidade(outpdf, ID, name):
	endereco = os.path.join(CONFIG.template, "ancestralidade.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)

	style = CONFIG.styles["ancestralidade.ances-name"]
	c.setFont(style.fontName, style.fontSize)
	c.setFillColor(style.textColor)

	txt = name.split(" ")[0].capitalize() + ", de acordo com seu mapeamento, você possui"
	c.drawString(100, 710, txt)

	endereco = os.path.join("Brutos", f"{ID}.txt")
	if os.path.exists(endereco) == False:
		endereco = os.path.join("Brutos", f"{ID}_23andMe.txt")
	ances = AncestralidadeMaker(Paciente(endereco))

	# .highest devolve a população e o valor percentual da maior população do paciente
	high, high_pct = ances.highest()

	# roundown arredonda para casas de 10
	high_pct = roundown(high_pct)

	# ajeita o nome da ancestralidade que será mostrado no pdf
	high = switch_highest(high)

	par = Paragraph(f">{high_pct}%", style=CONFIG.styles["ancestralidade.ances-number"])
	par.wrap(300, 100)
	par.drawOn(c, 80, 670)

	style = CONFIG.styles["ancestralidade.ances-high-top"]
	c.setFont(style.fontName, style.fontSize)
	c.setFillColor(style.textColor)
	c.drawString(270, 650, "de ancestralidade")

	style = CONFIG.styles["ancestralidade.ances-high-bot"]
	c.setFont(style.fontName, style.fontSize)
	c.setFillColor(style.textColor)
	c.drawString(270, 620, high)


	# cria a lista de valores de populações e porcentagens do paciente
	values = []
	for val in ances.ancestries:
		txt = str(val.population).split(".")[1]
		numb = str(val.value).split(".")
		pct = numb[0] + "." + numb[1][0]
		values.append([txt, float(pct)])

	ancestry = []
	# Soma populações de mesmo ID em um mesmo item
	for item in values:
		check = 0
		for val in ancestry:
			if ances.display_name(item[0]) in ances.display_name(val[0]):
				val[1] += item[1]
				check = 1
		if check == 0:
			ancestry.append([item[0], item[1]])

	ancestry = sorted(ancestry, key=lambda x: x[1], reverse = True)

	posy = 517
	for item in ancestry:
		style = CONFIG.styles["ancestralidade.ances-text"]
		c.setFont(style.fontName, style.fontSize)
		c.drawString(57, posy, str(ances.display_name(item[0])))
		style = CONFIG.styles["ancestralidade.ances-pct"]
		c.setFont(style.fontName, style.fontSize)
		if item[1] > 99:
			c.drawString(425, posy, "99%")
		else:
			c.drawString(425, posy, f"{item[1]}%")
		posy-= 16
	c.drawString(425, posy, "<1.00%")
	style = CONFIG.styles["ancestralidade.ances-text"]
	c.setFont(style.fontName, style.fontSize)
	c.drawString(57, posy, "Outros")

	pins = []
	# checa se as 2 maiores populações são de continentes diferentes
	for item in ancestry:
		if ances.pin_coordinates(ances.display_name(item[0])) not in pins:
			pins.append(ances.pin_coordinates(ances.display_name(item[0])))
		if len(pins) > 1:
			break

	endereco = os.path.join(CONFIG.template, "Circle Big.png")
	c.drawImage(endereco, pins[0][0], pins[0][1], width=39, height=39, mask='auto')
	if len(pins) > 1:
		endereco = os.path.join(CONFIG.template, "Circle Small.png")
		c.drawImage(endereco, pins[1][0], pins[1][1], width=25, height=25, mask='auto')

	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)
	print("Ancestralidade gerada\n")

def make_petal_color(outpdf, canv, name, posx, posy, pos, x, y, snp, big):
	num = find_color(name, snp)
	print('valor do num:',num)
	endereco=None
	if big == 0:
		big = ""
	else:
		big = "Big"
	#print("color", num)
	if num == 1:
		#print(f"Pet{big}G{pos}.png")
		endereco = os.path.join(os.path.join("Files", "Holobionte"), f"Pet{big}G{pos}.png")
	elif num == 2:
		#print(f"Pet{big}Y{pos}.png")
		endereco = os.path.join(os.path.join("Files", "Holobionte"), f"Pet{big}Y{pos}.png")
	elif num == 3:
		#print(f"Pet{big}R{pos}.png")
		endereco = os.path.join(os.path.join("Files", "Holobionte"), f"Pet{big}R{pos}.png")
	else:
		#print("\n\nERRO Pétalas\n\n")
		print('valor inesperados')
	canv.drawImage(endereco, posx, posy, width=x, height=y, mask='auto')

def find_impactful(outpdf, canv, name, snp, x, y):
	endereco = os.path.join("Files", "Lista.txt")
	end2 = os.path.join("Files", "Efeitos.txt")
	check = 0
	alto = []
	if name == "Glicação":
		name = "Insulina e Glicose"
	if name == "Intestinal":
		name = "Saúde Intestinal"
	with open(endereco, "r") as file:
		reading = file.readlines()
		for line in reading:
			sep = line.replace("\n", "").split("\t")
			if check == 1:
				if sep[0] == "":
					break
				elif sep[0][0] == 'r' and sep[0][1] == 's':
					for i in range(len(snp)):
						sp = snp[i].split("\t")
						with open(end2, "r") as f:
							r = f.readlines()
							for l in r:
								sep2 = l.replace("\n", "").split("\t")
								if sep2[0] == name and sep2[2] == sep[0] and sep2[3] == sp[1]:
									if sep[0] == sp[0] and sep2[4] == "2" and sp[3] not in alto:
										if ", " in sp[3]:
											genes = sp[3].split(", ")
											for gene in genes:
												alto.append(gene)
										else:
											alto.append(sp[3])
										break
			elif sep[0] == name:
				check = 1

	style = CONFIG.styles["ancestralidade.gene"]
	canv.setFont(style.fontName, style.fontSize)
	canv.setFillColor(style.textColor)
	col = 0
	line = 1
	tp = ""
	if len(alto) == 0:
		canv.drawString(x, y-10*line, "Sem genes de")
		canv.drawString(x, y-20*line, "alto risco")
	while 0 < len(alto):
		if line >= 9:
			print("Texto do Holobionte saindo do espaço definido. O programa será terminado.")
			exit(-5)
		if col == 0:
			tp = alto.pop(0)
		elif len(tp)+2+len(alto[0]) > 18:
			pass
		else:
			tp = f"{tp}, {alto.pop(0)}"
		col += 1
		if col == 3:
			canv.drawString(x, y-10*line, tp)
			col = 0
			line += 1
	if col != 0:
		canv.drawString(x, y-10*line, tp)
	style = CONFIG.styles["ancestralidade.holo"]
	canv.setFont(style.fontName, style.fontSize)
	canv.setFillColor(style.textColor)

def find_impactful_visao_geral(outpdf, canv, name, snp, x, y):
	endereco = os.path.join("Files", "Lista.txt")
	end2 = os.path.join("Files", "Efeitos.txt")
	check = 0
	alto = []
	with open(endereco, "r") as file:
		reading = file.readlines()
		for line in reading:
			sep = line.replace("\n", "").split("\t")
			if check == 1:
				if sep[0] == "":
					break
				elif sep[0][0] == 'r' and sep[0][1] == 's':
					for i in range(len(snp)):
						sp = snp[i].split("\t")
						with open(end2, "r") as f:
							r = f.readlines()
							for l in r:
								sep2 = l.replace("\n", "").split("\t")
								if sep2[0] == name and sep2[2] == sep[0] and sep2[3] == sp[1]:
									if sep[0] == sp[0] and sep2[4] == "2" and sp[3] not in alto:
										if ", " in sp[3]:
											genes = sp[3].split(", ")
											for gene in genes:
												alto.append(gene)
										else:
											alto.append(sp[3])
										break
			elif sep[0] == name:
				check = 1

	style = CONFIG.styles["ancestralidade.visao-geral-risco"]
	canv.setFont(style.fontName, style.fontSize)
	canv.setFillColor(style.textColor)
	#canv.setFillColorRGB(0,0,0)
	tp = ""
	x += 7
	y += 3
	if len(alto) == 0:
		canv.drawString(x, y, "Sem genes de alto risco")
	while 0 < len(alto):
			if tp != "":
				tp = f"{tp}, {alto.pop(0)}"
			else:
				tp = alto.pop(0)
	canv.drawString(x, y, tp)
	style = CONFIG.styles["ancestralidade.holo"]
	canv.setFont(style.fontName, style.fontSize)
	canv.setFillColor(style.textColor)

def find_impactful_blue(outpdf, canv, name, snp, x, y):
	endereco = os.path.join("Files", "Lista.txt")
	end2 = os.path.join("Files", "Efeitos.txt")
	if name == "Movimento":
		names = ["Danos e Lesões", "Desempenho em Esportes de Endurance", "Adaptabilidade esportiva", "Desempenho em Força e Potência", "Elevação da Frequência Cardíaca", "Colágeno e Articulações"]
	elif name == "Neurológico":
		names = ["Transtornos de Humor", "Enxaqueca", "Memória e Atenção", "Ciclo Circadiano", "Doenças Neurodegenerativas", "Canabinóides", "Desempenho Motor", "Comportamentos Adictos", "Emoções e Comportamento"]
	elif name == "Cardiovascular":
		names = ["Saúde Cardiovascular", "Dislipidemias", "Hematologia", "Pressão Arterial Sistêmica"]
	alto = []
	for item in names:
		print(item)
		check = 0
		with open(endereco, "r") as file:
			reading = file.readlines()
			for line in reading:
				sep = line.replace("\n", "").split("\t")
				if check == 1:
					if sep[0] == "":
						break
					elif sep[0][0] == 'r' and sep[0][1] == 's':
						print(sep[0])
						for i in range(len(snp)):
							sp = snp[i].split("\t")
							with open(end2, "r") as f:
								r = f.readlines()
								for l in r:
									sep2 = l.replace("\n", "").split("\t")
									if sep2[0] == item and sep2[2] == sep[0] and sep2[3] == sp[1]:
										if sep[0] == sp[0] and sep2[4] == "2" and sp[3] not in alto:
											if ", " in sp[3]:
												genes = sp[3].split(", ")
												for gene in genes:
													alto.append(gene)
											else:
												alto.append(sp[3])
				elif sep[0] == item:
					check = 1
	style = CONFIG.styles["ancestralidade.gene"]
	canv.setFont(style.fontName, style.fontSize)
	canv.setFillColor(style.textColor)
	col = 0
	line = 1
	tp = ""
	if len(alto) == 0:
		canv.drawString(x, y-10*line, "Sem genes de")
		canv.drawString(x, y-20*line, "alto risco")
	while 0 < len(alto):
		if line >= 9:
			print("Texto do Holobionte saindo do espaço definido. O programa será terminado.")
			exit(-5)
		if col == 0:
			tp = alto.pop(0)
		elif len(tp)+2+len(alto[0]) > 18:
			pass
		else:
			tp = f"{tp}, {alto.pop(0)}"
		col += 1
		if col == 3:
			canv.drawString(x, y-10*line, tp)
			col = 0
			line += 1
	if col != 0:
		canv.drawString(x, y-10*line, tp)
	style = CONFIG.styles["ancestralidade.holo"]
	canv.setFont(style.fontName, style.fontSize)
	canv.setFillColor(style.textColor)

def holobionte(outpdf, snp):
	endereco = os.path.join(CONFIG.template, "holobionte.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)

	x = 1373 * 0.4
	y = 1197 * 0.4
	posx = (595-x)/2
	posy = (872-y)/2
	names = ["Inflamação", "Glicação", "Biogênese Mitocondrial", "Detoxificação", "Imunidade", "Metilação"]
	names_outer = ["Alzheimer", "Parkinson", "Esclerose", "Depressão", "Esquizofrenia", "Ansiedade"]
	count = 0
	
	style = CONFIG.styles["ancestralidade.holo-pet-big"]
	c.setFont(style.fontName, style.fontSize)
	c.setFillColor(style.textColor)

	# for i in range(1,5):
	# 	endereco = os.path.join(os.path.join("Files", "Holobionte"), f"PetBig{i}.png")
	# 	c.drawImage(endereco, posx, posy, width=x, height=y, mask='auto')

	for i in range(6):
		endereco = os.path.join(os.path.join("Files", "Holobionte"), f"Ball{i}.png")
		c.drawImage(endereco, posx, posy, width=x, height=y, mask='auto')
		make_petal_color(outpdf, c, names[i], posx, posy, i+1, x, y, snp, 0)
		make_petal_color(outpdf, c, names_outer[i], posx, posy, i+1, x, y, snp, 1)

	endereco = os.path.join(os.path.join("Files", "Holobionte"), "PetalSmall.png")
	c.drawImage(endereco, posx, posy, width=x, height=y, mask='auto')

	# c.drawString(460, 443, "Intestinal")
	# find_impactful(outpdf, c, "Intestinal", snp, 460, 443)
	
	# c.setFont(style.fontName, style.fontSize)
	# c.setFillColor(style.textColor)
	# c.drawString(355, 610, "Cardiovascular")
	# find_impactful_blue(outpdf, c, "Cardiovascular", snp, 355, 610)
	
	# c.setFont(style.fontName, style.fontSize)
	# c.setFillColor(style.textColor)
	# c.drawString(160, 610, "Neurológico")
	# find_impactful_blue(outpdf, c, "Neurológico", snp, 160, 610)
	
	# c.setFont(style.fontName, style.fontSize)
	# c.setFillColor(style.textColor)
	# c.drawString(70, 443, "Movimento")
	# find_impactful_blue(outpdf, c, "Movimento", snp, 70, 443)

	style = CONFIG.styles["ancestralidade.holo"]
	c.setFont(style.fontName, style.fontSize)
	c.setFillColor(style.textColor)

	# Escrevendo nas pétalas internas e externas
	# Código comentado escreve os genes de risco de cada pétala externa, bastando cadastrar cada rota corretamente nos arquivos
	count = 0
	c.drawString(272, 570, names[count])
	find_impactful(outpdf, c, names[count], snp, 272, 570)
	c.drawString(65, 458, names_outer[count])
	find_impactful(outpdf, c, names_outer[count], snp, 65, 458)
	count += 1
	c.drawString(365, 513, names[count])
	find_impactful(outpdf, c, names[count], snp, 365, 513)
	c.drawString(155, 620, names_outer[count])
	find_impactful(outpdf, c, names_outer[count], snp, 155, 620)
	count += 1
	space = names[count].split(" ", 1)
	c.drawString(365, 388, space[0])
	c.drawString(365, 378, space[1])
	find_impactful(outpdf, c, names[count], snp, 365, 378)
	c.drawString(365, 620, names_outer[count])
	find_impactful(outpdf, c, names_outer[count], snp, 365, 620)
	count += 1
	c.drawString(260, 315, names[count])
	find_impactful(outpdf, c, names[count], snp, 260, 315)
	c.drawString(465, 458, names_outer[count])
	find_impactful(outpdf, c, names_outer[count], snp, 465, 458)
	count += 1
	c.drawString(170, 388, names[count])
	find_impactful(outpdf, c, names[count], snp, 170, 388)
	c.drawString(365, 305, names_outer[count])
	find_impactful(outpdf, c, names_outer[count], snp, 365, 305)
	count += 1
	c.drawString(165, 513, names[count])
	find_impactful(outpdf, c, names[count], snp, 165, 513)
	c.drawString(155, 305, names_outer[count])
	find_impactful(outpdf, c, names_outer[count], snp, 155, 305)
	endereco = os.path.join(os.path.join("Files", "Holobionte"), "homem.png")
	c.drawImage(endereco, ((595-x*0.9)/2), ((842-y*0.9)/2), width=x*0.9, height=y*0.9, mask='auto')

	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)
	print("Holobionte gerado\n")

def rotas_nutrientes(outpdf,snp):
	endereco = os.path.join(CONFIG.template, "rotas-nutrientes.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	endereco = os.path.join("Files", "Rotas.txt")
	found = 0
	size = 750
	col = 0
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			item = line.replace("\n", "")
			if found == 1 and item == "":
				break
			if found == 1:
				h = count_genes_size(item)+28
				# se daria overflow na coluna, gera nova coluna
				if size-h-30 <= 0:
					size = 750
					col += 1
				# contabiliza a legenda para a última coluna
				if col == 2 and size-h-222 <= 0:
					size = 750
					col += 1
				# se daria overflow na pagina, gera nova pagina
				if col > 2:
					col = 0
					c.save()
					inpdf = PdfReader(packet)
					for pageidx in range(len(template.pages)):
						pagina: PageObject = template.pages[pageidx]
						pagina.merge_page(inpdf.pages[pageidx])
					outpdf.add_page(pagina)
					endereco = os.path.join(CONFIG.template, "rotas-nutrientes.pdf")
					template = PdfReader(open(endereco, "rb"), strict=False)
					packet = BytesIO()
					c = canvas.Canvas(packet, pagesize=A4)
				size -= h + 10
				draw_rs_rotas(outpdf, snp, c, item, col, size+5)
				draw_box_rota(col,size,170,h,24,item,c,snp)
			if item == "Necessidade de Nutrientes":
				found = 1
	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)

def rotas_dietas(outpdf,snp):
	endereco = os.path.join(CONFIG.template, "rotas-dietas.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	endereco = os.path.join("Files", "Rotas.txt")
	found = 0
	size = 750
	col = 0
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			item = line.replace("\n", "")
			if found == 1 and item == "":
				break
			if found == 1:
				h = count_genes_size(item)+28
				if size-h-30 <= 0:
					size = 750
					col += 1
				# contabiliza a legenda para a última coluna
				if col == 2 and size-h-222 <= 0:
					size = 750
					col += 1
				if col > 2:
					col = 0
					c.save()
					inpdf = PdfReader(packet)
					for pageidx in range(len(template.pages)):
						pagina: PageObject = template.pages[pageidx]
						pagina.merge_page(inpdf.pages[pageidx])
					outpdf.add_page(pagina)
					endereco = os.path.join(CONFIG.template, "rotas-dietas.pdf")
					template = PdfReader(open(endereco, "rb"), strict=False)
					packet = BytesIO()
					c = canvas.Canvas(packet, pagesize=A4)
				size -= h + 10
				draw_rs_rotas(outpdf, snp, c, item, col, size+5)
				draw_box_rota(col,size,170,h,24,item,c,snp)
			if item == "Resposta à Dietas":
				found = 1
	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)

def rotas_sistemico(outpdf,snp,sex):
	endereco = os.path.join(CONFIG.template, "rotas-sistemico.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	endereco = os.path.join("Files", "Rotas.txt")
	found = 0
	size = 750
	col = 0
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			item = line.replace("\n", "")
			if found == 1 and item == "":
				break
			if found == 1:
				if item == 'Saúde Mamária' or item == 'Endometriose':
					if 'm' in sex or 'M' in sex:
						continue
				h = count_genes_size(item)+28
				if size-h-30 <= 0:
					size = 750
					col += 1
				# contabiliza a legenda para a última coluna
				if col == 2 and size-h-222 <= 0:
					size = 750
					col += 1
				if col > 2:
					col = 0
					c.save()
					inpdf = PdfReader(packet)
					for pageidx in range(len(template.pages)):
						pagina: PageObject = template.pages[pageidx]
						pagina.merge_page(inpdf.pages[pageidx])
					outpdf.add_page(pagina)
					endereco = os.path.join(CONFIG.template, "rotas-sistemico.pdf")
					template = PdfReader(open(endereco, "rb"), strict=False)
					packet = BytesIO()
					c = canvas.Canvas(packet, pagesize=A4)
				size -= h + 10
				draw_rs_rotas(outpdf, snp, c, item, col, size+5)
				draw_box_rota(col,size,170,h,24,item,c,snp)
			if item == "Sistêmico":
				found = 1
	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)

def rotas_atividades(outpdf,snp):
	endereco = os.path.join(CONFIG.template, "rotas-atividades.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	endereco = os.path.join("Files", "Rotas.txt")
	found = 0
	size = 750
	col = 0
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			item = line.replace("\n", "")
			if found == 1 and item == "":
				break
			if found == 1:
				h = count_genes_size(item)+28
				if size-h-30 <= 0:
					size = 750
					col += 1
				# contabiliza a legenda para a última coluna
				if col == 2 and size-h-222 <= 0:
					size = 750
					col += 1
				if col > 2:
					col = 0
					c.save()
					inpdf = PdfReader(packet)
					for pageidx in range(len(template.pages)):
						pagina: PageObject = template.pages[pageidx]
						pagina.merge_page(inpdf.pages[pageidx])
					outpdf.add_page(pagina)
					endereco = os.path.join(CONFIG.template, "rotas-atividades.pdf")
					template = PdfReader(open(endereco, "rb"), strict=False)
					packet = BytesIO()
					c = canvas.Canvas(packet, pagesize=A4)
				size -= h + 10
				draw_rs_rotas(outpdf, snp, c, item, col, size+5)
				draw_box_rota(col,size,170,h,24,item,c,snp)
			if item == "Atividades Físicas":
				found = 1
	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)

def rotas_cardio(outpdf,snp):
	endereco = os.path.join(CONFIG.template, "rotas-cardio.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	endereco = os.path.join("Files", "Rotas.txt")
	found = 0
	size = 750
	col = 0
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			item = line.replace("\n", "")
			if found == 1 and item == "":
				break
			if found == 1:
				h = count_genes_size(item)+28
				if size-h-30 <= 0:
					size = 750
					col += 1
				# contabiliza a legenda para a última coluna
				if col == 2 and size-h-222 <= 0:
					size = 750
					col += 1
				if col > 2:
					col = 0
					c.save()
					inpdf = PdfReader(packet)
					for pageidx in range(len(template.pages)):
						pagina: PageObject = template.pages[pageidx]
						pagina.merge_page(inpdf.pages[pageidx])
					outpdf.add_page(pagina)
					endereco = os.path.join(CONFIG.template, "rotas-cardio.pdf")
					template = PdfReader(open(endereco, "rb"), strict=False)
					packet = BytesIO()
					c = canvas.Canvas(packet, pagesize=A4)
				size -= h + 10
				draw_rs_rotas(outpdf, snp, c, item, col, size+5)
				draw_box_rota(col,size,170,h,24,item,c,snp)
			if item == "Saúde Cardiovascular":
				found = 1
	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)

def rotas_saudemental(outpdf,snp):
	endereco = os.path.join(CONFIG.template, "rotas-saudemental.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	endereco = os.path.join("Files", "Rotas.txt")
	found = 0
	size = 750
	col = 0
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			item = line.replace("\n", "")
			if found == 1 and item == "":
				break
			if found == 1:
				h = count_genes_size(item)+28
				if size-h-30 <= 0:
					size = 750
					col += 1
				# contabiliza a legenda para a última coluna
				if col == 2 and size-h-222 <= 0:
					size = 750
					col += 1
				if col > 2:
					col = 0
					c.save()
					inpdf = PdfReader(packet)
					for pageidx in range(len(template.pages)):
						pagina: PageObject = template.pages[pageidx]
						pagina.merge_page(inpdf.pages[pageidx])
					outpdf.add_page(pagina)
					endereco = os.path.join(CONFIG.template, "rotas-saudemental.pdf")
					template = PdfReader(open(endereco, "rb"), strict=False)
					packet = BytesIO()
					c = canvas.Canvas(packet, pagesize=A4)
				size -= h + 10
				draw_rs_rotas(outpdf, snp, c, item, col, size+5)
				draw_box_rota(col,size,170,h,24,item,c,snp)
			if item == "Saúde Mental":
				found = 1
	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)

def rotas_energia(outpdf,snp):
	endereco = os.path.join(CONFIG.template, "rotas-energia.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	endereco = os.path.join("Files", "Rotas.txt")
	found = 0
	size = 750
	col = 0
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			item = line.replace("\n", "")
			if found == 1 and item == "":
				break
			if found == 1:
				h = count_genes_size(item)+28
				if size-h-30 <= 0:
					size = 750
					col += 1
				# contabiliza a legenda para a última coluna
				if col == 2 and size-h-222 <= 0:
					size = 750
					col += 1
				if col > 2:
					col = 0
					c.save()
					inpdf = PdfReader(packet)
					for pageidx in range(len(template.pages)):
						pagina: PageObject = template.pages[pageidx]
						pagina.merge_page(inpdf.pages[pageidx])
					outpdf.add_page(pagina)
					endereco = os.path.join(CONFIG.template, "rotas-energia.pdf")
					template = PdfReader(open(endereco, "rb"), strict=False)
					packet = BytesIO()
					c = canvas.Canvas(packet, pagesize=A4)
				size -= h + 10
				draw_rs_rotas(outpdf, snp, c, item, col, size+5)
				draw_box_rota(col,size,170,h,24,item,c,snp)
			if item == "Energia e Metabolismo":
				found = 1
	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)

def rotas_sinalizacao(outpdf,snp):
	endereco = os.path.join(CONFIG.template, "rotas-sinalizacao.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	endereco = os.path.join("Files", "Rotas.txt")
	found = 0
	size = 750
	col = 0
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			item = line.replace("\n", "")
			if found == 1 and item == "":
				break
			if found == 1:
				h = count_genes_size(item)+28
				if size-h-30 <= 0:
					size = 750
					col += 1
				# contabiliza a legenda para a última coluna
				if col == 2 and size-h-222 <= 0:
					size = 750
					col += 1
				if col > 2:
					col = 0
					c.save()
					inpdf = PdfReader(packet)
					for pageidx in range(len(template.pages)):
						pagina: PageObject = template.pages[pageidx]
						pagina.merge_page(inpdf.pages[pageidx])
					outpdf.add_page(pagina)
					endereco = os.path.join(CONFIG.template, "rotas-sinalizacao.pdf")
					template = PdfReader(open(endereco, "rb"), strict=False)
					packet = BytesIO()
					c = canvas.Canvas(packet, pagesize=A4)
				size -= h + 10
				draw_rs_rotas(outpdf, snp, c, item, col, size+5)
				draw_box_rota(col,size,170,h,24,item,c,snp)
			if item == "Sinalização Celular":
				found = 1
	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)

def rotas_resultados(outpdf):
	endereco = os.path.join(CONFIG.template, "capa-rotas.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
	outpdf.add_page(pagina)

def rotas_descri(outpdf):
	endereco = os.path.join(CONFIG.template, "capa-rotas2.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
	outpdf.add_page(pagina)

def draw_blue_box(x,y,w,h,canv):
	canv.setStrokeColorRGB(0.56,0.75,0.82)
	canv.setLineWidth(0.5)
	canv.roundRect(x,y,w,h, 4, stroke=1, fill=0)

def draw_risco_visao(canv):
	x = 83
	y = 727
	w = 132
	h = 18
	rect_style = CONFIG.styles["ancestralidade.rect-green"]
	canv.setFillColor(rect_style.textColor)
	canv.rect(x,y,w,h, stroke=0, fill=1)
	style(canv, "ancestralidade.ances-text")
	canv.drawCentredString(x+w/2, y+5, "Baixo Risco")
	x+=15+w
	rect_style = CONFIG.styles["ancestralidade.rect-yellow"]
	canv.setFillColor(rect_style.textColor)
	canv.rect(x,y,w,h, stroke=0, fill=1)
	style(canv, "ancestralidade.ances-text")
	canv.drawCentredString(x+w/2, y+5, "Risco Moderado")
	x+=15+w
	rect_style = CONFIG.styles["ancestralidade.rect-red"]
	canv.setFillColor(rect_style.textColor)
	canv.rect(x,y,w,h, stroke=0, fill=1)
	style(canv, "ancestralidade.ances-text")
	canv.drawCentredString(x+w/2, y+5, "Alto Risco")

def visao_geral_test(outpdf, snp, sex):
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	draw_risco_visao(c)
	posx = 25
	h = 650
	endereco = os.path.join("Files", "Rotas.txt")
	title = ""
	names = []
	
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			item = line.replace("\n", "")
			if item == "":
				if h < 100:
					print("hzinho", h)
					c.save()
					inpdf = PdfReader(packet)
					endereco = os.path.join(CONFIG.template, "rotas.pdf")
					template = PdfReader(open(endereco, "rb"), strict=False)
					for pageidx in range(len(template.pages)):
						pagina: PageObject = template.pages[pageidx]
						pagina.merge_page(inpdf.pages[pageidx])
					outpdf.add_page(pagina)
					packet = BytesIO()
					c = canvas.Canvas(packet, pagesize=A4)
					draw_risco_visao(c)
					h = 650
				par = Paragraph(title, style=CONFIG.styles["ancestralidade.visao-geral-title"])
				par.wrap(564, 20)
				par.drawOn(c, posx, h+28)

				counter = 0

				if 'm' in sex or 'M' in sex:
					if "Saúde Mamária" in names:
						names.remove("Saúde Mamária")
					if "Endometriose" in names:
						names.remove("Endometriose")
				while len(names) > 0:
					if h < 50:
						size = 30+counter*48
						draw_blue_box(posx-6,h+18,564,size,c)
						c.save()
						inpdf = PdfReader(packet)
						endereco = os.path.join(CONFIG.template, "rotas.pdf")
						template = PdfReader(open(endereco, "rb"), strict=False)
						for pageidx in range(len(template.pages)):
							pagina: PageObject = template.pages[pageidx]
							pagina.merge_page(inpdf.pages[pageidx])
						outpdf.add_page(pagina)
						packet = BytesIO()
						c = canvas.Canvas(packet, pagesize=A4)
						draw_risco_visao(c)
						h = 650
						par = Paragraph(title, style=CONFIG.styles["ancestralidade.visao-geral-title"])
						par.wrap(564, 20)
						par.drawOn(c, posx, h+28)
						counter = 0
					rota = names.pop(0)
					make_rect_color(rota,c,posx,h,snp,0)
					c.setFillColorRGB(255,255,255)
					c.drawString(32, h, rota)
					find_impactful_visao_geral(outpdf, c, rota, snp, posx, h-24)
					h -= 48
					counter += 1

				#size = 30+counter*48-check_title*27
				size = 30+counter*48
				draw_blue_box(posx-6,h+18,564,size,c)
				h -= 45
				print(title)
				title = ""
				names = []
			elif title == "":
				title = item
			else:
				names.append(item)

	c.save()
	inpdf = PdfReader(packet)
	endereco = os.path.join(CONFIG.template, "rotas.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)

def visao_geral(outpdf, snp, sex):
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	draw_risco_visao(c)
	posx = 25
	h = 650
	endereco = os.path.join("Files", "Rotas.txt")
	title = ""
	counter = 0

	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			item = line.replace("\n", "")
			if title == "":
				if counter != 0:
					draw_blue_box(posx-6,h+11,564,36+counter*48,c)
					h -= 45
					if h < 100:
						c.save()
						try:
							# Abra o arquivo PDF usando PdfReader
							inpdf = PdfReader(packet)
							
							# Obtenha o valor de startxref
							startxref = inpdf.startxref
							print(f"startxref: {startxref}")
						except Exception as e:
							print(f"Erro ao processar o PDF: {e}")

						#inpdf = PdfReader(packet)
						endereco = os.path.join(CONFIG.template, "rotas.pdf")
						template = PdfReader(open(endereco, "rb"), strict=False)
						for pageidx in range(len(template.pages)):
							pagina: PageObject = template.pages[pageidx]
							pagina.merge_page(inpdf.pages[pageidx])
						outpdf.add_page(pagina)
						c = canvas.Canvas(packet, pagesize=A4)
						draw_risco_visao(c)
						posx = 25
						h = 650
				counter = 0
				title = item
				par = Paragraph(title, style=CONFIG.styles["ancestralidade.visao-geral-title"])
				par.wrap(564, 20)
				par.drawOn(c, posx, h+28)
			elif item == "":
				title = ""
			else:
				name = item
				if name == 'Saúde Mamária' or name == 'Endometriose':
					if 'm' in sex or 'M' in sex:
						continue
				make_rect_color(name,c,posx,h,snp,0)
				c.setFillColorRGB(255,255,255)
				c.drawString(32, h, name)
				make_rect_color(name,c,posx,h-24,snp,0)
				c.setFillColorRGB(255,255,255)
				find_impactful_visao_geral(outpdf, c, name, snp, posx, h-24)
				h -= 48
				counter += 1
		draw_blue_box(posx-6,h+11,564,36+counter*48,c)
		h -= 45
	c.save()
	inpdf = PdfReader(packet)
	endereco = os.path.join(CONFIG.template, "rotas.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)

def sum_genes(outpdf):
	endereco = os.path.join(CONFIG.template, "capa-sumario.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
	outpdf.add_page(pagina)

def contatos(outpdf):
	endereco = os.path.join(CONFIG.template, "verso.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
	outpdf.add_page(pagina)

def laudo_salvar(outpdf, name):
	today = date.today().strftime("%m-%d")
	endereco = os.path.join("Laudos Gerados", str(today))
	if not os.path.isdir(endereco):
		os.makedirs(endereco)
	endereco = os.path.join(endereco, f"Assaly {name}.pdf")
	with open(endereco, "wb") as file:
		outpdf.write(file)

def make_descri_box(name, canv, snp, posy, height):
	canv.setStrokeColorRGB(0.56,0.75,0.82)
	height += 47
	canv.rect(27,posy-height+6,537,height, stroke=1, fill=0)
	make_rect_color(name,canv,26.5,posy,snp,537.5)
	canv.setFillColorRGB(255,255,255)
	canv.drawString(41, posy, name)

def write_descri(name, canv, posy, dicio):

	par = Paragraph(dicio[name], style=CONFIG.styles["ancestralidade.text-regular"])
	w, h = par.wrap(513, 1)
	lines = h/9
	lines -= 2
	par.drawOn(canv, 39, posy - 40 - lines*9)
	return h

def make_descri_generic(names,c,snp,dicio):
	posy = 701
	posy_step = 90

	for i in range(len(names)):
		h = write_descri(names[i], c, posy, dicio)
		make_descri_box(names[i], c, snp, posy, h)
		posy -= posy_step+h-15

def descri_sinalizacao(outpdf,snp,dicio):
	endereco = os.path.join(CONFIG.template, "descricao-sinalizacao.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	found = 0
	names = []
	endereco = os.path.join("Files", "Rotas.txt")
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			item = line.replace("\n", "")
			if found == 1 and item == "":
				break
			if found == 1:
				names.append(item)
			if item == "Sinalização Celular":
				found = 1
	if len(names) > 7:
		make_descri_generic(names[:7], c, snp, dicio)
		c.save()
		inpdf = PdfReader(packet)
		for pageidx in range(len(template.pages)):
			pagina: PageObject = template.pages[pageidx]
			pagina.merge_page(inpdf.pages[pageidx])
		outpdf.add_page(pagina)
		endereco = os.path.join(CONFIG.template, "descricao-sinalizacao.pdf")
		template = PdfReader(open(endereco, "rb"), strict=False)
		packet = BytesIO()
		c = canvas.Canvas(packet, pagesize=A4)
		make_descri_generic(names[7:], c, snp, dicio)
	else:
		make_descri_generic(names, c, snp, dicio)
	
	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)
	print("Descrição Sinalização gerada\n")

def descri_sistemico(outpdf,snp,dicio,sex):
	endereco = os.path.join(CONFIG.template, "descricao-sistemico.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	found = 0
	names = []
	endereco = os.path.join("Files", "Rotas.txt")
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			item = line.replace("\n", "")
			if found == 1 and item == "":
				break
			if found == 1:
				if item == 'Saúde Mamária' or item == 'Endometriose':
					if 'm' in sex or 'M' in sex:
						continue
				names.append(item)
			if item == "Sistêmico":
				found = 1
	if len(names) > 7:
		make_descri_generic(names[:7], c, snp, dicio)
		c.save()
		inpdf = PdfReader(packet)
		for pageidx in range(len(template.pages)):
			pagina: PageObject = template.pages[pageidx]
			pagina.merge_page(inpdf.pages[pageidx])
		outpdf.add_page(pagina)
		endereco = os.path.join(CONFIG.template, "descricao-sistemico.pdf")
		template = PdfReader(open(endereco, "rb"), strict=False)
		packet = BytesIO()
		c = canvas.Canvas(packet, pagesize=A4)
		make_descri_generic(names[7:], c, snp, dicio)
	else:
		make_descri_generic(names, c, snp, dicio)
	
	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)
	print("Descrição Sistêmico gerada\n")

# def descri_sistemico2(outpdf,snp,sex,dicio):
# 	endereco = os.path.join(CONFIG.template, "descricao-sistemico.pdf")
# 	packet = BytesIO()
# 	c = canvas.Canvas(packet, pagesize=A4)
# 	names = ["Longevidade", "Capacidade Auditiva"]
# 	if 'f' in sex or 'F' in sex:
# 		names.append("Endometriose")
# 		names.append("Saúde Mamária")
# 	template = PdfReader(open(endereco, "rb"), strict=False)
	
# 	make_descri_generic(names, c, snp, dicio)

# 	c.save()
# 	inpdf = PdfReader(packet)
# 	for pageidx in range(len(template.pages)):
# 		pagina: PageObject = template.pages[pageidx]
# 		pagina.merge_page(inpdf.pages[pageidx])
# 	outpdf.add_page(pagina)
# 	print("Descrição Sistêmico 2 gerada\n")

def descri_saudemental(outpdf,snp,dicio):
	found = 0
	names = []
	endereco = os.path.join("Files", "Rotas.txt")
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			item = line.replace("\n", "")
			if found == 1 and item == "":
				break
			if found == 1:
				names.append(item)
			if item == "Saúde Mental":
				found = 1
	if len(names) > 7:		
		tam = len(names)
		low = 0
		high = 0

		while (tam) > 0:
			endereco = os.path.join(CONFIG.template, "descricao-saudemental.pdf")
			template = PdfReader(open(endereco, "rb"), strict=False)
			packet = BytesIO()
			c = canvas.Canvas(packet, pagesize=A4)
			min_heigh = min(tam, 7)
			if(min_heigh == 7):
				high += 7
			else : 
				high += min_heigh

			make_descri_generic(names[low:high], c, snp, dicio)
			tam -= min_heigh
			low = high
			c.save()
			inpdf = PdfReader(packet)
			for pageidx in range(len(template.pages)):
				pagina: PageObject = template.pages[pageidx]
				pagina.merge_page(inpdf.pages[pageidx])
			outpdf.add_page(pagina)


		#make_descri_generic(names[7:], c, snp, dicio)
		
		#inpdf = PdfReader(packet)
		#for pageidx in range(len(template.pages)):
			#pagina: PageObject = template.pages[pageidx]
			#pagina.merge_page(inpdf.pages[pageidx])
		#outpdf.add_page(pagina)
		#endereco = os.path.join(CONFIG.template, "descricao-saudemental.pdf")
		#template = PdfReader(open(endereco, "rb"), strict=False)
		#packet = BytesIO()
		#c = canvas.Canvas(packet, pagesize=A4)
		#outpdf.add_page(pagina)
		#make_descri_generic(names[14:], c, snp, dicio)

	else:
		make_descri_generic(names, c, snp, dicio)
	
	#c.save()
	#inpdf = PdfReader(packet)
	#for pageidx in range(len(template.pages)):
	#	pagina: PageObject = template.pages[pageidx]
	#	pagina.merge_page(inpdf.pages[pageidx])
	#outpdf.add_page(pagina)
	print("Descrição Saúde Mental gerada\n")

# def descri_saudemental2(outpdf,snp,dicio):
# 	endereco = os.path.join(CONFIG.template, "descricao-saudemental.pdf")
# 	template = PdfReader(open(endereco, "rb"), strict=False)
# 	packet = BytesIO()
# 	c = canvas.Canvas(packet, pagesize=A4)
# 	names = ["Desempenho Motor", "Resposta ao Tratamento para Depressão", "Canabinóides"]
	
# 	make_descri_generic(names, c, snp, dicio)

# 	c.save()
# 	inpdf = PdfReader(packet)
# 	for pageidx in range(len(template.pages)):
# 		pagina: PageObject = template.pages[pageidx]
# 		pagina.merge_page(inpdf.pages[pageidx])
# 	outpdf.add_page(pagina)
# 	print("Descrição Saúde Mental 2 gerada\n")

def descri_cardio(outpdf,snp,dicio):
	endereco = os.path.join(CONFIG.template, "descricao-cardio.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	found = 0
	names = []
	endereco = os.path.join("Files", "Rotas.txt")
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			item = line.replace("\n", "")
			if found == 1 and item == "":
				break
			if found == 1:
				names.append(item)
			if item == "Saúde Cardiovascular":
				found = 1
	if len(names) > 7:
		make_descri_generic(names[:7], c, snp, dicio)
		c.save()
		inpdf = PdfReader(packet)
		for pageidx in range(len(template.pages)):
			pagina: PageObject = template.pages[pageidx]
			pagina.merge_page(inpdf.pages[pageidx])
		outpdf.add_page(pagina)
		endereco = os.path.join(CONFIG.template, "descricao-cardio.pdf")
		template = PdfReader(open(endereco, "rb"), strict=False)
		packet = BytesIO()
		c = canvas.Canvas(packet, pagesize=A4)
		make_descri_generic(names[7:], c, snp, dicio)
	else:
		make_descri_generic(names, c, snp, dicio)
	
	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)
	print("Descrição Cardiovascular gerada\n")

def descri_energia(outpdf,snp,dicio):
	endereco = os.path.join(CONFIG.template, "descricao-energia.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	found = 0
	names = []
	endereco = os.path.join("Files", "Rotas.txt")
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			item = line.replace("\n", "")
			if found == 1 and item == "":
				break
			if found == 1:
				names.append(item)
			if item == "Energia e Metabolismo":
				found = 1
	if len(names) > 7:
		make_descri_generic(names[:7], c, snp, dicio)
		c.save()
		inpdf = PdfReader(packet)
		for pageidx in range(len(template.pages)):
			pagina: PageObject = template.pages[pageidx]
			pagina.merge_page(inpdf.pages[pageidx])
		outpdf.add_page(pagina)
		endereco = os.path.join(CONFIG.template, "descricao-energia.pdf")
		template = PdfReader(open(endereco, "rb"), strict=False)
		packet = BytesIO()
		c = canvas.Canvas(packet, pagesize=A4)
		make_descri_generic(names[7:], c, snp, dicio)
	else:
		make_descri_generic(names, c, snp, dicio)
	
	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)
	print("Descrição Energia gerada\n")

def descri_atividades(outpdf,snp,dicio):
	endereco = os.path.join(CONFIG.template, "descricao-atividades.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	found = 0
	names = []
	endereco = os.path.join("Files", "Rotas.txt")
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			item = line.replace("\n", "")
			if found == 1 and item == "":
				break
			if found == 1:
				names.append(item)
			if item == "Atividades Físicas":
				found = 1
	if len(names) > 7:
		make_descri_generic(names[:7], c, snp, dicio)
		c.save()
		inpdf = PdfReader(packet)
		for pageidx in range(len(template.pages)):
			pagina: PageObject = template.pages[pageidx]
			pagina.merge_page(inpdf.pages[pageidx])
		outpdf.add_page(pagina)
		endereco = os.path.join(CONFIG.template, "descricao-atividades.pdf")
		template = PdfReader(open(endereco, "rb"), strict=False)
		packet = BytesIO()
		c = canvas.Canvas(packet, pagesize=A4)
		make_descri_generic(names[7:], c, snp, dicio)
	else:
		make_descri_generic(names, c, snp, dicio)
	
	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)
	print("Descrição Atividades gerada\n")

def descri_nutrientes(outpdf,snp,dicio):
	endereco = os.path.join(CONFIG.template, "descricao-nutrientes.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	found = 0
	names = []
	endereco = os.path.join("Files", "Rotas.txt")
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			item = line.replace("\n", "")
			if found == 1 and item == "":
				break
			if found == 1:
				names.append(item)
			if item == "Necessidade de Nutrientes":
				found = 1
	if len(names) > 7:
		make_descri_generic(names[:7], c, snp, dicio)
		c.save()
		inpdf = PdfReader(packet)
		for pageidx in range(len(template.pages)):
			pagina: PageObject = template.pages[pageidx]
			pagina.merge_page(inpdf.pages[pageidx])
		outpdf.add_page(pagina)
		endereco = os.path.join(CONFIG.template, "descricao-nutrientes.pdf")
		template = PdfReader(open(endereco, "rb"), strict=False)
		packet = BytesIO()
		c = canvas.Canvas(packet, pagesize=A4)
		make_descri_generic(names[7:], c, snp, dicio)
	else:
		make_descri_generic(names, c, snp, dicio)
	
	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)
	print("Descrição Nutrientes gerada\n")

def descri_dietas(outpdf,snp,dicio):
	endereco = os.path.join(CONFIG.template, "descricao-dietas.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	found = 0
	names = []
	endereco = os.path.join("Files", "Rotas.txt")
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			item = line.replace("\n", "")
			if found == 1 and item == "":
				break
			if found == 1:
				names.append(item)
			if item == "Resposta à Dietas":
				found = 1
	if len(names) > 7:
		make_descri_generic(names[:7], c, snp, dicio)
		c.save()
		inpdf = PdfReader(packet)
		for pageidx in range(len(template.pages)):
			pagina: PageObject = template.pages[pageidx]
			pagina.merge_page(inpdf.pages[pageidx])
		outpdf.add_page(pagina)
		endereco = os.path.join(CONFIG.template, "descricao-dietas.pdf")
		template = PdfReader(open(endereco, "rb"), strict=False)
		packet = BytesIO()
		c = canvas.Canvas(packet, pagesize=A4)
		make_descri_generic(names[7:], c, snp, dicio)
	else:
		make_descri_generic(names, c, snp, dicio)
	
	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)
	print("Descrição Dietas gerada\n")

def gene_efeitos(outpdf,snp):
	end2 = os.path.join("Files", "Efeitos.txt")
	baixo = []
	medio = []
	alto = []
	semef = []
	print(snp)
	print(len(snp))
	print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
	for i in range(len(snp)):
		val = int(0)
		sp = snp[i].split("\t")
		if sp[1] == "--":
			semef.append(f"{sp[0]}\t{sp[1]}")
		else:
			with open(end2, "r") as f:
				reading = f.readlines()
				for l in reading:
					sep2 = l.replace("\n", "").split("\t")
					print(sep2)
					if sep2[2] == sp[0] and sep2[3] == sp[1]:
						if int(sep2[4]) > int(val):
							val = int(sep2[4])
			if str(val) == "0":
				baixo.append(f"{sp[0]}\t{sp[1]}\t{sp[3]}")
			elif str(val) == "1":
				medio.append(f"{sp[0]}\t{sp[1]}\t{sp[3]}")
			elif str(val) == "2":
				alto.append(f"{sp[0]}\t{sp[1]}\t{sp[3]}")

	endereco = os.path.join(CONFIG.template, "genes-por-efeito.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	column = 0
	style = CONFIG.styles["ancestralidade.gene"]
	
	gene_style = CONFIG.styles["ancestralidade.text-red"]
	rep = 0
	while rep < len(alto):
		if column == 3:
			column = 0
			c.save()
			inpdf = PdfReader(packet)
			for pageidx in range(len(template.pages)):
				pagina: PageObject = template.pages[pageidx]
				pagina.merge_page(inpdf.pages[pageidx])
			outpdf.add_page(pagina)
			endereco = os.path.join(CONFIG.template, "genes-por-efeito.pdf")
			template = PdfReader(open(endereco, "rb"), strict=False)
			packet = BytesIO()
			c = canvas.Canvas(packet, pagesize=A4)
		posx = switch_column(column)
		posy = 692
		count = 0
		while count < 22 and (count+rep) < len(alto):
			sp = alto[count+rep].split("\t")
			c.setFont(style.fontName, style.fontSize)
			c.setFillColor(style.textColor)
			c.drawString(posx, posy, sp[0])
			for i in range(len(snp)):
				seu_gene = snp[i].split("\t")
				if sp[0] == seu_gene[0]:
					c.setFont(gene_style.fontName, gene_style.fontSize)
					c.setFillColor(gene_style.textColor)
					c.drawString(posx, posy+10, seu_gene[3])
					c.drawString(posx+130, posy+10, seu_gene[1])
			count += 1
			posy -= 30
		size = count*30+28
		draw_box_rota(column,740-size,170,size,24,"ALTO RISCO",c,snp)
		column += 1
		rep += count

	gene_style = CONFIG.styles["ancestralidade.text-yellow"]
	rep = 0
	while rep < len(medio):
		if column == 3:
			column = 0
			c.save()
			inpdf = PdfReader(packet)
			for pageidx in range(len(template.pages)):
				pagina: PageObject = template.pages[pageidx]
				pagina.merge_page(inpdf.pages[pageidx])
			outpdf.add_page(pagina)
			endereco = os.path.join(CONFIG.template, "genes-por-efeito.pdf")
			template = PdfReader(open(endereco, "rb"), strict=False)
			packet = BytesIO()
			c = canvas.Canvas(packet, pagesize=A4)
		posx = switch_column(column)
		posy = 692
		count = 0
		while count < 22 and (count+rep) < len(medio):
			sp = medio[count+rep].split("\t")
			c.setFont(style.fontName, style.fontSize)
			c.setFillColor(style.textColor)
			c.drawString(posx, posy, sp[0])
			for i in range(len(snp)):
				seu_gene = snp[i].split("\t")
				if sp[0] == seu_gene[0]:
					c.setFont(gene_style.fontName, gene_style.fontSize)
					c.setFillColor(gene_style.textColor)
					c.drawString(posx, posy+10, seu_gene[3])
					c.drawString(posx+130, posy+10, seu_gene[1])
			count += 1
			posy -= 30
		size = count*30+28
		draw_box_rota(column,740-size,170,size,24,"MÉDIO RISCO",c,snp)
		column += 1
		rep += count

	gene_style = CONFIG.styles["ancestralidade.text-green"]
	rep = 0
	while rep < len(baixo):
		if column == 3:
			column = 0
			c.save()
			inpdf = PdfReader(packet)
			for pageidx in range(len(template.pages)):
				pagina: PageObject = template.pages[pageidx]
				pagina.merge_page(inpdf.pages[pageidx])
			outpdf.add_page(pagina)
			endereco = os.path.join(CONFIG.template, "genes-por-efeito.pdf")
			template = PdfReader(open(endereco, "rb"), strict=False)
			packet = BytesIO()
			c = canvas.Canvas(packet, pagesize=A4)
		posx = switch_column(column)
		posy = 692
		count = 0
		while count < 22 and (count+rep) < len(baixo):
			sp = baixo[count+rep].split("\t")
			c.setFont(style.fontName, style.fontSize)
			c.setFillColor(style.textColor)
			c.drawString(posx, posy, sp[0])
			for i in range(len(snp)):
				seu_gene = snp[i].split("\t")
				if sp[0] == seu_gene[0]:
					c.setFont(gene_style.fontName, gene_style.fontSize)
					c.setFillColor(gene_style.textColor)
					c.drawString(posx, posy+10, seu_gene[3])
					c.drawString(posx+130, posy+10, seu_gene[1])
			count += 1
			posy -= 30
		size = count*30+28
		draw_box_rota(column,740-size,170,size,24,"BAIXO RISCO",c,snp)
		column += 1
		rep += count

	gene_style = CONFIG.styles["ancestralidade.text-grey"]
	rep = 0
	while rep < len(semef):
		if column == 3:
			column = 0
			c.save()
			inpdf = PdfReader(packet)
			for pageidx in range(len(template.pages)):
				pagina: PageObject = template.pages[pageidx]
				pagina.merge_page(inpdf.pages[pageidx])
			outpdf.add_page(pagina)
			endereco = os.path.join(CONFIG.template, "genes-por-efeito.pdf")
			template = PdfReader(open(endereco, "rb"), strict=False)
			packet = BytesIO()
			c = canvas.Canvas(packet, pagesize=A4)
		posx = switch_column(column)
		posy = 692
		count = 0
		while count < 22 and (count+rep) < len(semef):
			sp = semef[count+rep].split("\t")
			c.setFont(style.fontName, style.fontSize)
			c.setFillColor(style.textColor)
			c.drawString(posx, posy, sp[0])
			for i in range(len(snp)):
				seu_gene = snp[i].split("\t")
				if sp[0] == seu_gene[0]:
					c.setFont(gene_style.fontName, gene_style.fontSize)
					c.setFillColor(gene_style.textColor)
					c.drawString(posx, posy+10, seu_gene[3])
					c.drawString(posx+130, posy+10, seu_gene[1])
			count += 1
			posy -= 30
		size = count*30+28
		draw_box_rota(column,740-size,170,size,24,"SEM EFEITO",c,snp)
		column += 1
		rep += count

	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)
	print("Sumário de Genes gerado\n")