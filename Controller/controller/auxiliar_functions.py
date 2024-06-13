import os
import re
from tqdm import tqdm
from google.oauth2 import service_account
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload
from Controller.controller.logger import *
from Controller.controller.logger import logger_str
from datetime import datetime
DNA_REGEX = re.compile(r"[ACTGD-]{2}")


def config_logger(name = DEFAULT_LOGGER_NAME):
    
    logger = logging.getLogger(name)
    global DEFAULT_LOGGER_NAME 
    DEFAULT_LOGGER_NAME = name

    log_formatter = logging.Formatter('%(asctime)s %(levelname)s %(funcName)s(linha: %(lineno)d de %(module)s.py) %(message)s', 
                                  datefmt='%d/%m/%Y %H:%M:%S')
    now = datetime.now().strftime('%d_%m_%Y-%H_%M_%S')
    logFile = f"../logs/log-{now}.txt"

    #Setup File handler
    file_handler = logging.FileHandler(logFile)
    file_handler.setFormatter(log_formatter)
    file_handler.setLevel(logging.INFO)

    #Setup Stream Handler (i.e. console)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(log_formatter)
    stream_handler.setLevel(logging.INFO)

    #Get our logger
    logger.setLevel(logging.INFO)

    #Add both Handlers
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    #Write some Data
    logger.info(f"Início do programa: {now}")

def get_logger():
    return logging.getLogger(DEFAULT_LOGGER_NAME)

def get_data(service, ID, df_norm):
        """Retorna os dados do arquivo de genótipos.
        Args:
            file: Caminho do arquivo de entrada.
        Returns:
            dict[str, SNP]: SNPs obtidos e seus genótipos na forma rsid: SNP.
        """
        logger = get_logger()
        # pasta no drive de dados brutos
        file_id = get_file_id(service, ID) #TODO: log caso file_id = None (n ache codigo bruto)
        if file_id == None:
            logger.error("Erro tipo 1 (usuário/cliente)\n\tFalta de dados brutos")
            raise Exception()
        filename = f"bruto_{ID}.txt"
        file_download(service, file_id, filename)
        path_input_ancestralidade = f'../Controller/report/iadmix/input/input_{ID}.txt'# caminho e nome do arquivo de input para a ancestralidade
        alelos_dict = {}  # Dicionário para mapear rsids para alelos
        with open(filename, encoding="utf-8") as handle:# abre o arquivo de brutos para leitura
            
            delimiter = get_file_delimiter(filename, encoding="utf-8")
            total_lines = sum(1 for _ in open(f"./{filename}", encoding="utf-8"))
            progress_bar = tqdm(total=total_lines)

            with open(path_input_ancestralidade, 'w') as file:# abre o arquivo de input da ancestralidade para escrever dados
                for line in handle:
                    if re.match(r'rs\d+', line):
                        # Agora que encontramos uma linha que começa com "rsXXXX", podemos parar de pular o cabeçalho
                        break
                    # Atualizar a barra de progresso
                    progress_bar.update(1)

                callback = get_callback_fromline(line, delimiter)
                for line in handle:
                    if not line.strip() or line.startswith("#") or not(re.match(r'rs\d+', line)):
                        continue
                    rsids, a1, a2 = callback(line, delimiter)
                    #TODO:
                    if(a1 not in ["A","T","C","G","D","I"] or a2 not in ["A","T","C","G","D","I"]):
                        a1 = "-"
                        a2 = "-"
                    if "," in rsids:
                        rsids = rsids.split(",")
                        for rsid in rsids:
                            rsid = rsid.strip()
                            if rsid and DNA_REGEX.fullmatch(f"{a1}{a2}"):
                                alelos_dict[rsid] = f"{a1}{a2}"
                                if (a1 != "-" or a2 != "-"):
                                    file.write(f"{rsid}\t{a1}{a2}\n")
                    elif rsids and DNA_REGEX.fullmatch(f"{a1}{a2}"):
                        alelos_dict[rsids] = f"{a1}{a2}"
                        if (a1 != "-" or a2 != "-"):
                            file.write(f"{rsids}\t{a1}{a2}\n")
                    # Atualizar a barra de progresso
                    progress_bar.update(1)
        
        # Fechar a barra de progresso quando terminar
        
        progress_bar.close()
        os.remove(filename) #remove arquivo temporario
        # Atualize a coluna "alelos" do DataFrame com os valores do dicionário
        df_norm['alelos'] = df_norm['rsid'].map(alelos_dict)
        df_norm['alelos'].fillna('--', inplace=True)
        string_erro = "Erro tipo 3 (dados brutos/normscore)\n\t" + "Alelos de efeito e risco inconsistentes entre dados brutos e normscore!"
        flag_erro = False
        for _, row in df_norm.iterrows():
                interest_pair = row["effect-allele"]
                read_pair = row["alelos"]
                if read_pair != "--" and (read_pair[0] not in interest_pair or read_pair[1] not in interest_pair):
                    flag_erro = True
                    string_erro += f"\n\tSNP lido {row['rsid']} não coincide com o par cadastrado: Leu {read_pair[0]} e " \
                          + f"{read_pair[1]}, esperava {interest_pair[0]} ou {interest_pair[2]}"
        if(flag_erro == True):
            logger.error(string_erro)
            
        # Substitua os valores vazios na coluna 'alelos' por '--'
        return df_norm

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
        print('Erro tipo 3 (dados brutos/normscore)\n\tDado bruto em formato inadequado')
        raise Exception()

def extract_specific_fields(transformed_data):
    extracted_data = []
    for item in transformed_data:
        # Acessamos os campos específicos do dicionário
        extracted_data.append({
            'Trilha': item['Trilha'],
            'Snp': item['Snps'],
            'Gene': item['Gene']
        })
    return extracted_data

def transform_data(sheet_data):
    headers = sheet_data[0]
    transformed_data = []
    
    # Em seguida, iteramos sobre cada linha de dados (ignorando os cabeçalhos)
    for row in sheet_data[1:]:
        # Criamos um dicionário para cada linha, associando os cabeçalhos aos valores correspondentes
        data_dict = {headers[i]: row[i] for i in range(len(headers))}
        transformed_data.append(data_dict)
    
    return transformed_data


def get_file_encoding(filepath: str):
    
    for encoding in ('utf-8', 'latin-1'):
        try:
            with open(filepath, "r", newline="", encoding=encoding) as handle:
                handle.read()
            return encoding
        except UnicodeDecodeError:
            pass
    print('Erro tipo 3 (dados brutos/normscore)\n\tDado bruto em formato inadequado (incapaz de determinar o encoding)')
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