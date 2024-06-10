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
from Controller.report.ancestralidade import AncestralidadeMaker
from Controller.report.base import Paciente
from reportlab.lib.colors import HexColor

def style(canv, name):
	style = CONFIG.styles[name]
	canv.setFont(style.fontName, style.fontSize)
	canv.setFillColor(style.textColor)

def dicio_descri():
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Descri.txt")
	dicio = {}
	with open(endereco, "r") as file:
		r = file.readlines()
		for line in r:
			sp = line.replace("\n", "").split("\t")
			if len(sp) == 2:
				dicio[sp[0]] = sp[1]
			else:
				print(f"Line format error: {line}")
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
		print('box rota, com snp - ',snp)
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
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Lista.txt")
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
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Lista.txt")
	end2 = os.path.join("../Controller", "DataFiles", "Files", "Efeitos.txt")
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

def read_file(file_path):
	with open(file_path, "r") as file:
		return file.readlines()
	
def find_color(trait_name, snp_data):
	print(snp_data)
	sleep
	
	
	print('Trait being processed:', trait_name)

	lista_file_path = os.path.join("../Controller", "DataFiles", "Files", "Lista.txt")
	effects_file_path = os.path.join('../Controller', "DataFiles", "Files", "Efeitos.txt")
	routes_file_path = os.path.join('../Controller', "DataFiles", "Files", "Rotas.txt")

	routes_dict = load_routes(routes_file_path)

	if trait_name == "Glicação":
		trait_name = "Insulina e Glicose"

	related_traits = get_route_characteristics(trait_name, routes_dict)
	if len(related_traits) < 2:
		related_traits = [trait_name]
	
	print('Related traits for', trait_name, ':', related_traits)
	
	total_effect = 0
	max_possible_effect = 0
	
	try:
		lista_lines = read_file(lista_file_path)
		effects_lines = read_file(effects_file_path)

		effects_dict = {}
		for line in effects_lines:
			line = line.strip()
			if not line:
				continue
			fields = line.split("\t")
			if len(fields) >= 5:
				trait = fields[0]
				snp_id = fields[2]
				genotype = fields[3]
				effect = int(fields[4])
				if trait not in effects_dict:
					effects_dict[trait] = {}
				if snp_id not in effects_dict[trait]:
					effects_dict[trait][snp_id] = {}
				effects_dict[trait][snp_id][genotype] = effect

				#ID=23177081028823
	except Exception as e:
		print('Error processing effects dictionary:', e)

	for trait in related_traits:
		print('Processing trait:', trait)
		found_trait = False

		with open(lista_file_path, "r") as lista_file:
			for line in lista_file:
				fields = line.strip().split("\t")

				if not found_trait:
					if fields[0] == trait:
						found_trait = True
						print(f'Trait {trait} found in Lista.txt')
					continue

				if fields[0] == "":
					break

				if fields[0].startswith('rs'):
					snp_id = fields[0]
					for snp_entry in snp_data:
						snp_fields = snp_entry.split("\t")
						if snp_id == snp_fields[0] and snp_fields[1] != "--":
							genotype = snp_fields[1]
							if (trait in effects_dict and snp_id in effects_dict[trait]
									and genotype in effects_dict[trait][snp_id]):
								effect = effects_dict[trait][snp_id][genotype]
								total_effect += effect
								max_possible_effect += 2
								print(f'Match found: Trait: {trait}, SNP: {snp_id}, Genotype: {genotype}, Effect: {effect}, Total Effect: {total_effect}, Max Possible Effect: {max_possible_effect}')
							else:
								print(f'SNP {snp_id} with genotype {genotype} not found in effects_dict for trait {trait}')
						else:
							print(f'SNP {snp_id} not matched with SNP entry or genotype is "--"')
				else:
					print(f'Line does not start with rs, ignoring: {fields[0]}')

	if max_possible_effect == 0:
		print("Max possible effect is 0, returning default color code 4")
		return 4  # Default color code if no effects found
	
	effect_ratio = total_effect / max_possible_effect
	print(f'Effect ratio calculated: {effect_ratio}')
	if effect_ratio < 0.3:
		print("Returning color code 1 (Green)")
		return 1  # Green
	elif effect_ratio < 0.7:
		print("Returning color code 2 (Yellow)")
		return 2  # Yellow
	else:
		print("Returning color code 3 (Red)")
		return 3  # Red
	
def read_file(file_path):

	with open(file_path, "r") as file:
		return file.readlines()

def load_routes(routes_file_path):
	"""
	Carrega as rotas e suas características a partir de um arquivo e as organiza em um dicionário.

	:param routes_file_path: Caminho para o arquivo que contém as rotas e suas características.
	:return: Dicionário com rotas como chaves e listas de características como valores.
	"""
	routes_dict = {}
	try:
		with open(routes_file_path, 'r') as file:
			lines = file.readlines()
		
		current_route = None

		for line in lines:
			line = line.strip()
			if not line:  # Ignorar linhas vazias
				current_route = None
				continue
			
			if current_route is None:
				# Detectar nova rota
				current_route = line
				routes_dict[current_route] = []
			else:
				# Adicionar característica à rota atual
				routes_dict[current_route].append(line)

	except FileNotFoundError:
		print(f"File not found: {routes_file_path}")
	except Exception as e:
		print(f"An error occurred: {e}")
	
	return routes_dict

def get_route_characteristics(name, routes_dict):
	if name == "":
		pass
	"""
	Obtém as características associadas a uma rota ou retorna o nome da característica se for uma característica.

	:param name: Nome da rota ou característica.
	:param routes_dict: Dicionário com rotas como chaves e listas de características como valores.
	:return: Lista de características associadas à rota ou a própria característica.
	"""
	if name in routes_dict:
		# Se o nome estiver nas rotas, retornar todas as características associadas a essa rota
		return routes_dict[name]
	else:
		# Caso contrário, retornar o nome como uma característica única
		return [name]

def make_rect_color(name, canv, posx, posy, snp, width):
	print('make rect color')
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

	endereco = os.path.join("../Controller", "DataFiles", "Brutos", f"{ID}.txt")
	if os.path.exists(endereco) == False:
		endereco = os.path.join("../Controller", "DataFiles", "Brutos", f"{ID}_23andMe.txt")
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
	if name == "":
		pass
	print('NOME QUE ESTA NA PETAL COLOR', name)
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
		endereco = os.path.join(os.path.join("../Controller", "DataFiles","Files", "Holobionte"), f"Pet{big}G{pos}.png")
	elif num == 2:
		print(f"Pet{big}Y{pos}.png")
		endereco = os.path.join(os.path.join("../Controller", "DataFiles", "Files", "Holobionte"), f"Pet{big}Y{pos}.png")
	elif num == 3:
		#print(f"Pet{big}R{pos}.png")
		endereco = os.path.join(os.path.join("../Controller", "DataFiles", "Files", "Holobionte"), f"Pet{big}R{pos}.png")
	else:
		#print("\n\nERRO Pétalas\n\n")
		print('valor inesperados')
	if name == " ":
		endereco = os.path.join(os.path.join("../Controller", "DataFiles", "Files", "Holobionte"), f"PetG{pos}.png")

	canv.drawImage(endereco, posx, posy, width=x, height=y, mask='auto')


def find_impactful(outpdf, canv, name, snp, x, y):
	if name == " ":
		pass
	print(f"Chamando find_impactful para {name} nas coordenadas ({x}, {y})")
	routes_file_path = os.path.join("../Controller", "DataFiles", "Files", "Rotas.txt")
	routes_dict = load_routes(routes_file_path)
	characteristics = get_route_characteristics(name, routes_dict)

	# Limitar a três características
	shown_characteristics = characteristics[:3]
	remaining_characteristics_count = max(0, len(characteristics) - 3)

	# Ajustar o tamanho da fonte diretamente
	canv.setFont("Helvetica", 7)  # Fonte Helvetica com tamanho 8
	canv.setFillColorRGB(0, 0, 0)  # Cor preta

	# Ajustar a posição inicial para um pouco abaixo do nome da rota
	y -= 12
	x+= 8
	print(f"Posição inicial ajustada para características: ({x}, {y})")

	if not shown_characteristics:
		print("Nenhuma característica encontrada")
		canv.drawString(x, y, "Sem características")
	else:
		for characteristic in shown_characteristics:
			print(f"Escrevendo característica '{characteristic}' nas coordenadas ({x}, {y})")
			canv.drawString(x, y, characteristic)
			y -= 12  # Ajustar a posição y para a próxima linha
	
	canv.setFont("Helvetica", 6)  # Fonte Helvetica com tamanho 8
	canv.setFillColorRGB(0, 0, 0)  # Cor preta

	if remaining_characteristics_count > 0:
		y -= 1
		x+= 10
		remaining_text = f"+ {remaining_characteristics_count}"
		print(f"Escrevendo informação adicional '{remaining_text}' nas coordenadas ({x}, {y})")
		canv.drawString(x, y, remaining_text)

	canv.setFont("Helvetica-Bold", 10)
	canv.setFillColorRGB(0.2, 0.2, 0.2)  # Cor cinza escuro para mais contraste


def find_impactful_visao_geral(outpdf, canv, name, snp, x, y):
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Lista.txt")
	end2 = os.path.join("../Controller", "DataFiles", "Files", "Efeitos.txt")
	check = 0
	alto = []
	with open(endereco, "r") as file:
		reading = file.readlines()
		for line in reading:
			sep = line.replace("\n", "").split("\t")
			if check == 1:
				if sep[0] == "":
					break
				elif sep[0].startswith('rs'):
					for snp_entry in snp:
						sp = snp_entry.split("\t")
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
	# canv.setFillColorRGB(0,0,0)

	x += 10
	y += 5
	if len(alto) == 0:
		canv.drawString(x, y, "Sem genes de alto risco")
	else:
		tp = ", ".join(alto)
		canv.drawString(x, y, tp)

	style = CONFIG.styles["ancestralidade.holo"]
	canv.setFont(style.fontName, style.fontSize)
	canv.setFillColor(style.textColor)

def find_impactful_blue(outpdf, canv, name, snp, x, y):
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Lista.txt")
	end2 = os.path.join("../Controller", "DataFiles", "Files", "Efeitos.txt")
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
	
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Template", "holobionte.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)

	c.setFont("Helvetica-Bold", 10)
	c.setFillColorRGB(0.2, 0.2, 0.2)  # Cor cinza escuro para mais contraste

	x = 1373 * 0.4
	y = 1197 * 0.4
	posx = (595 - x) / 2
	posy = (872 - y) / 2
	names = [" ", "Saúde Mental", "Doenças Crônicas", " ", "Envelhecimento", "Neuro"]
	names_outer = ["Atividades Físicas", "Saúde Cardiovascular", "Nutrição", 'Sinalização Celular', "Metabolismo", "Sistêmico"]

	num_names = len(names)
	num_names_outer = len(names_outer)

	for i in range(num_names):
		ball_path = os.path.join("../Controller", "DataFiles", "Files", "Holobionte", f"Ball{i}.png")
		print(f"Desenhando imagem: {ball_path}")
		c.drawImage(ball_path, posx, posy, width=x, height=y, mask='auto')
		make_petal_color(outpdf, c, names[i], posx, posy, i + 1, x, y, snp, 0)

		if i < num_names_outer:
			make_petal_color(outpdf, c, names_outer[i], posx, posy, i + 1, x, y, snp, 1)

	petal_small_path = os.path.join("../Controller", "DataFiles", "Files", "Holobionte", "PetalSmall.png")
	c.drawImage(petal_small_path, posx, posy, width=x, height=y, mask='auto')

	# Write names on petals with the original positions
	names_positions = [
		(278, 580),  # COR AZUL
		(375, 522),  # Saúde Mental
		(366, 409),  # Doença Crônica
		(295, 324),  # COR AZUL
		(170, 404),  # Envelhecimento
		(177, 520)  # Neuro
	]

	names_outer_positions = [
		(53, 458),   # Fisica
		(149, 620),  # Cardio
		(344, 625),  # Nutri
		(453, 458),  # Celular
		(359, 305),  # Metabolismo
		(144, 310)   # Sistemico
	]

	for count in range(num_names):
		posx_name, posy_name = names_positions[count]
		print(f"Escrevendo nome '{names[count]}' nas coordenadas ({posx_name}, {posy_name})")
		posy_name = posy_name - 14
		posx_name = posx_name - 14

		c.drawString(posx_name, posy_name, names[count])
		find_impactful(outpdf, c, names[count], snp, posx_name, posy_name)  # Ajustar posição y

		if count < num_names_outer:
			posx_outer, posy_outer = names_outer_positions[count]
			print(f"Escrevendo nome '{names_outer[count]}' nas coordenadas ({posx_outer}, {posy_outer})")
			c.drawString(posx_outer, posy_outer, names_outer[count])
			find_impactful(outpdf, c, names_outer[count], snp, posx_outer, posy_outer)

	homem_image_path = os.path.join("../Controller", "DataFiles", "Files", "Holobionte", "homem.png")
	c.drawImage(homem_image_path, (595 - x * 0.9) / 2, (842 - y * 0.9) / 2, width=x * 0.9, height=y * 0.9, mask='auto')

	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)

	print("Holobionte gerado\n")

def rotas_nutrientes(outpdf,snp):
	endereco = os.path.join(CONFIG.template, "rotas-nutrientes.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Rotas.txt")
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
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Rotas.txt")
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
	endereco = os.path.join("../Controller","DataFiles", "Files", "Rotas.txt")
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
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Rotas.txt")
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
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Rotas.txt")
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
	endereco = os.path.join("../Controller","DataFiles", "Files", "Rotas.txt")
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
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Rotas.txt")
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
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Rotas.txt")
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
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Rotas.txt")
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
					print('rota',rota)
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
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Rotas.txt")
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
				#make_rect_color(name,c,posx,h,snp,0)
				c.setFillColorRGB(255,255,255)
				c.drawString(32, h, name)
				#make_rect_color(name,c,posx,h-24,snp,0)
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
	#make_rect_color(name,canv,26.5,posy,snp,537.5)
	canv.setFillColorRGB(255,255,255)
	canv.drawString(41, posy, name)

def write_descri(name, canv, posy, dicio):
	print(dicio[name])
	
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
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Rotas.txt")
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
	endereco = os.path.join("../Controller","DataFiles", "Files", "Rotas.txt")
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
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Rotas.txt")
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
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Rotas.txt")
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
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Rotas.txt")
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
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Rotas.txt")
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
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Rotas.txt")
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
	endereco = os.path.join("../Controller", "DataFiles", "Files", "Rotas.txt")
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

def gene_efeitos(outpdf, snp):
	def read_effects_file(filepath):
		with open(filepath, "r") as file:
			return [line.strip().split("\t") for line in file.readlines()]

	def classify_snps(snp_list, effects):
		low_risk = []
		medium_risk = []
		high_risk = []
		no_effect = []

		for snp in snp_list:
			sp = snp.split("\t")
			if len(sp) < 2 or sp[1] == "--":
				no_effect.append(snp)
				continue

			max_val = 0
			for effect in effects:
				if len(effect) < 5:
					continue
				if effect[2] == sp[0] and effect[3] == sp[1]:
					max_val = max(max_val, int(effect[4]))

			if max_val == 0:
				low_risk.append(snp)
			elif max_val == 1:
				medium_risk.append(snp)
			elif max_val == 2:
				high_risk.append(snp)

		return low_risk, medium_risk, high_risk, no_effect

	def draw_genes(c, genes, title, style, gene_style, column):
		posx = switch_column(column)
		posy = 692
		count = 0
		for gene_info in genes:
			if count >= 22:
				break
			sp = gene_info.split("\t")
			if len(sp) < 3:
				continue
			c.setFont(style.fontName, style.fontSize)
			c.setFillColor(style.textColor)
			c.drawString(posx, posy, sp[0])
			c.setFont(gene_style.fontName, gene_style.fontSize)
			c.setFillColor(gene_style.textColor)
			c.drawString(posx, posy + 10, sp[3])
			c.drawString(posx + 130, posy + 10, sp[1])
			count += 1
			posy -= 30

		size = count * 30 + 28
		draw_box_rota(column, 740 - size, 170, size, 24, title, c, snp)
		return count

	end2 = os.path.join("../Controller", "DataFiles", "Files", "Efeitos.txt")
	effects = read_effects_file(end2)
	
	low_risk, medium_risk, high_risk, no_effect = classify_snps(snp, effects)

	endereco = os.path.join(CONFIG.template, "genes-por-efeito.pdf")
	template = PdfReader(open(endereco, "rb"), strict=False)
	packet = BytesIO()
	c = canvas.Canvas(packet, pagesize=A4)
	column = 0
	style = CONFIG.styles["ancestralidade.gene"]

	for risk_group, title, gene_style in [
		(high_risk, "ALTO RISCO", CONFIG.styles["ancestralidade.text-red"]),
		(medium_risk, "MÉDIO RISCO", CONFIG.styles["ancestralidade.text-yellow"]),
		(low_risk, "BAIXO RISCO", CONFIG.styles["ancestralidade.text-green"]),
		(no_effect, "SEM EFEITO", CONFIG.styles["ancestralidade.text-grey"])
	]:
		rep = 0
		while rep < len(risk_group):
			if column == 3:
				column = 0
				c.save()
				inpdf = PdfReader(packet)
				for pageidx in range(len(template.pages)):
					pagina: PageObject = template.pages[pageidx]
					pagina.merge_page(inpdf.pages[pageidx])
				outpdf.add_page(pagina)
				packet = BytesIO()
				c = canvas.Canvas(packet, pagesize=A4)

			posx = switch_column(column)
			posy = 692
			count = draw_genes(c, risk_group[rep:rep + 22], title, style, gene_style, column)
			column += 1
			rep += count

	c.save()
	inpdf = PdfReader(packet)
	for pageidx in range(len(template.pages)):
		pagina: PageObject = template.pages[pageidx]
		pagina.merge_page(inpdf.pages[pageidx])
	outpdf.add_page(pagina)
	print("Sumário de Genes gerado\n")
