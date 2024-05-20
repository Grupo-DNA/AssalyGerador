import os
import re
from tqdm import tqdm
from google.oauth2 import service_account
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload
from Controller.controller.auxiliar_functions import get_logger

DNA_REGEX = re.compile(r"[ACTGD-]{2}")

def get_credentials():
    SCOPES = ['https://www.googleapis.com/auth/drive.readonly']
    SERVICE_ACCOUNT_FILE = 'path/to/service_account.json'
    
    credentials = service_account.Credentials.from_service_account_file(
        SERVICE_ACCOUNT_FILE, scopes=SCOPES)
    return credentials

def get_drive_service():
    credentials = get_credentials()
    return build('drive', 'v3', credentials=credentials)

def get_callback_fromline(line: str, delimiter: str):
    logger = get_logger()
    def callback1(value: str, delim: str):
        rsid, _, (a1, a2) = value.rstrip("\n").split(delim)
        if rsid == ".":
            return None, None, None
        return rsid, a1, a2
    def callback2(value: str, delim: str):
        rsid, chrom, _, alelos = value.rstrip("\n").split(delim)
        a1 = None
        a2 = None
        if len(alelos) > 0:
            a1 = alelos[0]
        if len(alelos) == 2:
            a2 = alelos[1]
        if rsid == "." or len(alelos) > 2:
            return None, None, None
        if a2 == None and chrom in ["MT", "X", "Y"]:
            a2 = "-"
        return rsid, a1, a2
    def callback3(value: str, delim: str):
        rsid, _, a1, a2, _ = value.rstrip("\n").split(delim)
        if rsid == ".":
            return None, None, None
        return rsid, a1, a2
    def callback4(value: str, delim: str):
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
        logger.error('Erro tipo 3 (dados brutos/normscore)\n\tDado bruto em formato inadequado')
        raise Exception()

def get_file_encoding(filepath: str):
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
    if encoding is None:
        encoding = get_file_encoding(filepath)
    delims = ["\t", ",", ";"]
    with open(filepath, "r", newline="", encoding=encoding) as handle:
        sniffer = handle.readline().rstrip("\n")
        by_order = sorted(
            delims, key=sniffer.count,
            reverse=True)
        return by_order[0]

def download_file(service, file_id, filename):
    request = service.files().get_media(fileId=file_id)
    fh = open(filename, "wb")
    downloader = MediaIoBaseDownload(fh, request)
    done = False
    while done is False:
        status, done = downloader.next_chunk()
        if status:
            print("Download %d%%." % int(status.progress() * 100))
    fh.close()