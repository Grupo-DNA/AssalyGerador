from googleapiclient.discovery import build
from google.auth.transport.requests import Request
from google.oauth2.credentials import Credentials
from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.http import MediaIoBaseDownload, MediaFileUpload
from googleapiclient.errors import HttpError
import io, os, pandas as pd, random, math
from PIL import Image
import logging
from datetime import datetime
from time import time, sleep
from Controller.controller.logger import logger_str
import json
import ast
import csv

# If modifying these scopes, delete the file token.json.
SCOPES = ['https://www.googleapis.com/auth/spreadsheets.readonly','https://www.googleapis.com/auth/drive']
DEFAULT_LOGGER_NAME = "principal"
DB_LINK = '1bj0Hn9Oo0wl2TkwDFYGBEzm-RK1WDcOEo9d4uOTnfiU' # ID da planilha de "relação de IDs" (parte do link)
USERS_SHEET = 'Usuários BD Club 2.0' # Id da planilha da Club_relação-de-IDs-integração
NORMSCORE_LINK = '108KQpjPqOkWeOwH8z_2HST8EeXWWKJe-FFo2IpktFSI' # Id da planilha de Normscore Integrado
CUSTOMIZATION_SHEET = "Personalização"

def login():
    """Login for usage of the Drive v3 API.
    """
    """
    login faz autenticação de segurança para acessar arquivos pelo google api. É preciso um arquivo com as credenciais 
        ou um arquivo já com o token.
    A partir das credenciais, a própria função gera um arquivo de token para ser usado e acessado futuramente.

    :param: None
    :return: creds - credenciais de login para o google API
    """
    try:

        creds = None
        if os.path.exists('../Controller/AUTH/token.json'):
            
                creds = Credentials.from_authorized_user_file('../Controller/AUTH/token.json', SCOPES)
        # If there are no (valid) credentials available, let the user log in.
        if not creds or not creds.valid:
            if creds and creds.expired and creds.refresh_token:
                creds.refresh(Request())
            else:
                flow = InstalledAppFlow.from_client_secrets_file(
                    '../Controller/AUTH/credentials.json', SCOPES)
                creds = flow.run_local_server(port=0)

            # Save the credentials for the next run
            with open('../Controller/AUTH/token.json', 'w') as token:
                token.write(creds.to_json())
        return creds
    
    
    except:
        raise Exception('Erro TIPO 0 (infraestrutura)\n\tFalha no Login com Google Drive.')

def download_image(service, file_id, file_nickname = "", cache = True, identifier_to_filename = False):
    if identifier_to_filename:
        image_path = f'../Controller/DataFiles/imagens_drive/{file_nickname}_{file_id}.png'
    else:
        image_path = f'../Controller/DataFiles/imagens_drive/{file_id}.png'
    if not os.path.isfile(image_path) or not cache:
        file_download(service, file_id, image_path)
    return image_path

def get_file_id(drive_service, id_sample):
    # pasta no drive de dados brutos
    file = (
        drive_service.files()
        .list(
            q=f"mimeType='text/plain' and name contains '{id_sample}' and trashed = false",
            corpora = "drive",
            driveId = "0ACNfyK39P99DUk9PVA",
            includeItemsFromAllDrives=True, 
            supportsAllDrives=True,
            fields="files(id, name)"
        )
        .execute()
    )
    results = file.get("files", [])
    
    if not results:# folder do not exist yet
        return None
    else:
        return results[0].get("id")
    
def file_download(service, file_id, file_name):
    """
    file_download does the download of a file based in its id

    :param service: drive api client
    :param file_id: id of the file for download
    :param file_name: name for save the file in local machine
    :return: None
    """
    request = service.files().get_media(fileId=file_id)
    fh = io.BytesIO()
        
    # Initialise a downloader object to download the file
    downloader = MediaIoBaseDownload(fh, request)
    done = False
    while not done:
        status, done = downloader.next_chunk()
        #print(F'Download {int(status.progress() * 100)}.')
        
    # Write the received data to the file
    with open(os.path.join(file_name), 'wb') as f:
        f.write(fh.getvalue())
        f.close()
    return