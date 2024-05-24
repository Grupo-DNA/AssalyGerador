import logging

# Este arquivo contem as funcoes relacionadas ao log de erros

DEFAULT_LOGGER_NAME = "principal"

class logger_str:
    '''
        Contem todas as mensagens que podem aparecer no log de erros
    '''
    def __init__(self):
        pass

    def start_generation(self,ID:str,date_time:str) -> str:
        return ("===========================================================</p>" + 
                f"<ul><b>Geração do laudo {ID} começou em {date_time}</b></ul>")
    
    def drive_error(self,msg_err:str) -> str:
        return ("<font color='red'>Erro tipo 0 (infraestrutura):</font>Ao tentar conectar-se ao drive:" 
                + msg_err + "</p>")
    
    def fallen_internet(self) -> str:
        return ("<font color='red'>Erro tipo 0 (infraestrutura):</font> Falta de conexão à Internet! Tentando novamente...</p>")
    
    def abort(self,ID:str) -> str:
        return (f"<u><font color='darkred'>Abortando a geração do laudo {ID}!</u></p>")
    
    def small_prs(self, ID:str) -> str:
        return (f"<font color='yellow'>laudo {ID} com baixa compatibilidade.</font></p>")
    
    def time_spent(self,ID:str,timef:float,time_start:float) -> str:
        return (f"<u>Tempo gasto gerando o laudo para {ID}: {(timef - time_start)//60} minutos e {((timef - time_start)%60):.2f} segundos</u></p>")
    
    def statistics(self,laudos_gerados:list[str], n_laudos_gerados:int, n_total:int, 
                   laudos_baixa_compatibilidade:list[str], n_baixa_compatibilidade:int,
                   laudos_falhos:list[str],percem_duplo_traco:float,num_snip_repetido:int,
                   snips_repetidos_genotipos_diferentes:int) -> str:
        return (f"===========================================================</p>"+
        f"<ul><font color='darkblue'>Número de gerações de laudo bem-sucedidas: {len(laudos_gerados)} ({100*n_laudos_gerados/n_total}%)</font></ul>" +
        f"<ul><font color='darkblue'><u>IDs dos laudos bem-sucedidos: {laudos_gerados}</u></font></ul>" +
        f"<ul><font color='darkblue'>Número de laudos gerados mas com baixa compatibilidade: {len(laudos_baixa_compatibilidade)} ({100*n_baixa_compatibilidade/n_total}%)</font></ul>" +
        f"<ul><font color='darkblue'><u>IDs dos laudos gerados com baixa compatibilidade: {laudos_baixa_compatibilidade}</u></font></ul>" + 
        f"<ul><font color='darkblue'>Número de gerações de laudo que falharam: {len(laudos_falhos)} ({100*len(laudos_falhos)/(len(laudos_gerados)+len(laudos_falhos))}%)</font></ul>" +
        f"<ul><font color='darkblue'><u>IDs dos laudos falhos: {laudos_falhos})</u></font></ul>"+
        f"<ul><font color='darkblue'><u>Porcentagem de alelos '--': {100*(percem_duplo_traco):.2f}%)</u></font></ul>"+
        f"<ul><font color='darkblue'><u>Número de SNIPS repetidos: {num_snip_repetido})</u></font></ul>"+
        f"<ul><font color='darkblue'><u>Número de SNIPS repetidos que apresentam genótipos diferentes: {snips_repetidos_genotipos_diferentes})</u></font></ul>")
    
    def bad_data_formatter(self) -> str:
        return ("<font color = 'red'>Erro tipo 3 (dados brutos/normscore):</font> Dado bruto em formato inadequado.</p>")
    
    def missing_data(self) -> str:
        return ("<font color = 'red'>Erro tipo 1 (usuário/cliente):</font> Falta de dados brutos</p>")
    
    def wrong_sheet_link(self) -> str:
        return ("<font color='red'>Erro tipo 0 (infraestrutura):</font> Link para planilha inválido! Abortando geração!</p>")
    
    def sheet_access_failed(self, msg_err:str) -> str:
        return ("<font color='red'>Erro tipo 0 (infraestrutura):</font> Ao tentar acessar planilha:" + msg_err + "</p>")
    
    def permission_denied(self) -> str:
        return ("<font color='red'>Erro TIPO 0 (infraestrutura):</font> Falha no Login com Google Drive.</p>")
    
    def missing_header(self, header:str) -> str:
        return (f"<font color='red'>Erro TIPO 0 (infraestrutura):</font> Planilha das trilhas do laudo com problema de formatação (Cabeçalho não possui coluna '{header}').</p>")
    
    def access_failed(self, IDlaudo:str, msg_err:str) -> str:
        return (f"<font color='red'>Erro tipo 0 (infraestrutura):</font> Ao tentar acessar planilha {IDlaudo}: " + msg_err + "</p>")
    
    def missing_carac(self,missing_carac_list:list) -> str:
        return (f"font color='red'>Erro tipo 2 (personalização):</font> Falta de informações de texto das características: {missing_carac_list}.</p>")
    
    def normscoreID_not_found(self,header:str) -> str:
        return (f"<font color='red'>Erro TIPO 3 (Normscore): {header} não encontrado na planilha.</p>")
    
    def normscoresheet_not_found(self, msg_err:str) -> str:
        return ("<font color='red'>Erro tipo 0 (infraestrutura):</font> Ao tentar acessar planilha Normscore Global: " + msg_err + "</p>")

    def sheet_not_found(self, USERS_SHEET:str) -> str:
        return (f"<font color='red'>Erro TIPO 0 (infraestrutura):</font> Aba de {USERS_SHEET} não encontrada no BD. Verifique o nome.</p>")

    def ID_not_found(self) -> str:
        return ("<font color='red'>Erro TIPO 1 (cliente):</font> ID Club não encontrado na planilha.</p>")
    
    def duplicatedID(self) -> str:
        return ("<font color='red'>Erro TIPO 1 (cliente):</font> ID Club duplicado na planilha.</p>")
    
    def carac_not_generated(self, num_decartadas:int,total_caract:int) -> str:
        return (f"Entre todas as características previstas para a esse laudo {(num_decartadas*100)/total_caract:.2f}% não conseguiram ser geradas")

    def empty_file(self) -> str:
        return ('Erro TIPO 0: arquivo está vazio.')
    
    def laudo_not_found(self, IDlaudo, ID) -> str:
        return (f"Formato de laudo {IDlaudo} não encontrado na planilha para o {ID}.</p>")

def config_logger(now:str, name = DEFAULT_LOGGER_NAME) -> None:
    '''
        Inicializa o log e define seu formato e arquivo onde sera escrito.

        Args:
        - now: Tempo do início da execução do gerador, com data e horário.
        - name: Nome dado ao logger.
    '''    

    logger = logging.getLogger(name)
    global DEFAULT_LOGGER_NAME 
    DEFAULT_LOGGER_NAME = name

    log_formatter = logging.Formatter('<p>%(asctime)s %(levelname)s %(funcName)s(linha: %(lineno)d de %(module)s.py) %(message)s', 
                                  datefmt='%d/%m/%Y %H:%M:%S')
    logFile = f"./logs/log-{now}.txt"

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
    logger.info(f"Início do programa: {now}</p>")

# Para criar loggers dentro de outras funcoes 
def get_logger():
    '''
        Recebe um Logger já configurado. 
        
        Usada para não criar várias configurações loggers durante a execução
    '''
    return logging.getLogger(DEFAULT_LOGGER_NAME)