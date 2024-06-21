import os
import sys
sys.path.append('../')
from pypdf import PdfWriter
from Controller.report.laudos import *
from google.auth.transport.requests import Request
from google.oauth2.credentials import Credentials
from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError
from Controller.report.leituraDados import read_SNPs
from Controller.controller.auxiliar_functions import *

# If modifying these scopes, delete the file token.json.
SCOPES = ['https://www.googleapis.com/auth/spreadsheets.readonly']

# The ID and range of a sample spreadsheet.
SAMPLE_SPREADSHEET_ID = '1bj0Hn9Oo0wl2TkwDFYGBEzm-RK1WDcOEo9d4uOTnfiU'
#SAMPLE_RANGE_NAME = 'Usuários próprios BD Club!B3:E'
SAMPLE_RANGE_NAME = 'Usuários BD Club 2.0!B3:E'

# Autorização Planilha
creds = None
if os.path.exists('token.json'):
	creds = Credentials.from_authorized_user_file('token.json', SCOPES)
if not creds or not creds.valid:
	if creds and creds.expired and creds.refresh_token:
		creds.refresh(Request())
	else:
		flow = InstalledAppFlow.from_client_secrets_file(
			'credentials.json', SCOPES)
		creds = flow.run_local_server(port=0)
	with open('token.json', 'w') as token:
		token.write(creds.to_json())

try:
	service = build('sheets', 'v4', credentials=creds)
	sheet = service.spreadsheets()
	result = sheet.values().get(spreadsheetId=SAMPLE_SPREADSHEET_ID, range=SAMPLE_RANGE_NAME).execute()
	values = result.get('values', [])

except HttpError as err:
	print(err)

with open("gerar.txt", "r") as handle:
	reading = handle.readlines()
	for line in reading:
		# Pula linhas que comecem com # ou estejam em branco
		if line[0] == "#" or line == "\n":
			continue
		line = line.replace("\n", "")
		for row in values:
			if row[0] in line:
				name = row[1].strip().title()
				sex = row[3]
				ID = row[0]
				print(f"Gerando laudo {ID} - {name}\n")
				
				snp = read_SNPs(ID)
				
				if snp == -1:
					break
				outpdf = PdfWriter()
				dicio = dicio_descri()
				laudo_capa(outpdf, name)
				ancestralidade(outpdf, ID, name)
				holobionte(outpdf, snp)
				print('holobionte gerada')
				rotas_resultados(outpdf)
				visao_geral_test(outpdf, snp, sex)
				rotas_nutrientes(outpdf, snp)
				rotas_sistemico(outpdf, snp, sex)
				rotas_atividades(outpdf, snp)
				rotas_cardio(outpdf, snp)
				rotas_saudemental(outpdf, snp)
				rotas_energia(outpdf, snp)
				rotas_sinalizacao(outpdf, snp)
				rotas_descri(outpdf)
				print('Rotas Geradas')
				descri_sistemico(outpdf, snp, dicio, sex)
				descri_saudemental(outpdf, snp, dicio)
				descri_cardio(outpdf, snp, dicio)
				descri_Osteoarticular(outpdf,snp,dicio)
				descri_energia(outpdf, snp, dicio)
				descri_atividades(outpdf, snp, dicio)
				descri_nutrientes(outpdf, snp, dicio)
				
				print('Descricao gerada')
				sum_genes(outpdf)
				print('A chamar funcao de genes')
				gene_efeitos(outpdf, snp)
				contatos(outpdf)
				laudo_salvar(outpdf, name)
				print(f"Laudo Assaly {ID} - {name} gerado com sucesso\n\n")
				break