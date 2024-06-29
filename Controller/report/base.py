# -*- coding: utf-8 -*-
# packages/base.py

"""Módulo que contém os modelos de dados."""

import logging
import re
from datetime import date
from decimal import Decimal
from enum import Enum, Flag, auto
from functools import cached_property, singledispatchmethod
from math import prod, sqrt
from os import PathLike, remove
from pathlib import Path
from statistics import NormalDist
from typing import Iterable, Optional, TypedDict, Union

import pandas as pd
from reportlab.lib.colors import Color, HexColor

if __name__ in ("__main__", "base"):
    from Controller.controller.config import CONFIG
    from Controller.controller.exceptions import NoSNPFoundError, RestrictionError
    from Controller.controller.helpers import (csv_reader, eval_rule, get_file_delimiter,get_temp_filename, truncate)
    from iadmix.runancestry import ancestry_pipeline_genotypes
else:
    from Controller.controller.config import CONFIG
    from Controller.controller.exceptions import NoSNPFoundError, RestrictionError
    from Controller.controller.helpers import (csv_reader, eval_rule, get_file_delimiter,get_temp_filename, truncate)
    from .iadmix.runancestry import ancestry_pipeline_genotypes


logger = logging.getLogger("base")

DNA_REGEX = re.compile(r"[ACTGD]{2}")
cdf = NormalDist(0.0, 1.0).cdf
ZScoreSnpData = TypedDict('ZScoreSnpData', {
    "rsid": str,
    "effect_allele": str,
    "major_allele": str,
    "MAF": float,
    "effect_size": float
})


class Sexo(Flag):
    """Sexo do paciente."""
    M = auto()
    F = auto()
    INDEFINIDO = F | M
    AMBOS = INDEFINIDO
    MASCULINO = M
    FEMININO = F

    def __repr__(self):
        return f"{self.__class__.__name__}.{self.name.upper()}"  # type: ignore

    @classmethod
    def _missing_(cls, value: str):
        for item in cls:
            if str(item.name).lower() == value.lower():
                return item
        if value == "?":
            return cls.INDEFINIDO

    def format(self, txt: str):
        """Return an adapted string based on gender specific internal flags.

        ## Supported flags:

            `[oa]`\n
            `[se M]...[se F]...[senão]...` (order independent switches)
        """
        oa = {
            Sexo.MASCULINO: "o",
            Sexo.FEMININO: "a",
            Sexo.INDEFINIDO: "o(a)"
        }
        if "[oa]" in txt.lower():
            txt = txt.replace("[oa]", oa[self]).replace(
                "[OA]", oa[self].upper())

        anyflag = r"\[se(?: ?n[ãa]o| F| M)\]"
        flag = {
            Sexo.M: r"\[se M\]([^[]+)",
            Sexo.F: r"\[se F\]([^[]+)",
            Sexo.INDEFINIDO: r"\[se ?n[ãa]o\]([^[]+)",
        }[self]
        if match := re.search(anyflag, txt, re.IGNORECASE):
            if match := re.search(flag, txt, re.IGNORECASE):
                txt = match.group(1)
            else:
                raise ValueError(f"Tag {flag!r} não encontrada em {txt!r}")

        return txt


class SuperPopulation(Enum):
    """Código de super-população para cálculo de scores poligênicos.

    Atributos:
    ---------
        zscore_filename: str
            nome do arquivo contendo os dados para calcular z-scores.
        normscore_filename: str
            nome do arquivo contendo os dados para calcular score normalmente.
    """
    EUR = auto()
    EAS = auto()
    AMR = auto()
    AFR = auto()
    SAS = auto()
    GLOBAL = auto()

    @property
    def zscore_filename(self):
        """Nome do arquivo com os SNPs de Z-Score."""
        return f"SNPs_zscore_{self.name.upper()}.tsv"

    @property
    def normscore_filename(self):
        """Nome do arquivo com os SNPs de score oligogênico."""
        return f"SNPs_normscore_{self.name.upper()}.tsv"

    def __repr__(self):
        return f"{self.__class__.__name__}.{self.name.upper()}"


class Continente(Enum):
    """Continentes"""
    AMERICA_DO_SUL = auto()
    AMERICA_DO_NORTE = auto()
    EUROPA = auto()
    AFRICA = auto()
    ASIA = auto()
    OCEANIA = auto()


class Population(Enum):
    """Códigos de populações para Ancestralidade.

    Atributos:
    ----------
        display_name : str
            Nome que aparece sobre a barra do gráfico de barras.
        bar_color : str
            Cor da barra do gráfico de barras.
        continente : Continente
            Continente usado na escolha da cor do pin e do texto do resumo.
        superpop_code : SuperPopulation
            Código da superpopulação, usado para selecionar os pesos dos SNPs.
        pin_coordinates : tuple[float, float]
            Coordenadas do pin no PDF.
        pin_file : Path
            Arquivo com a imagem do pin.

    Métodos:
    ----------
        from_display_name(name: str) -> Self@Population
            Recupera a população com base em seu nome de display (evita duplicatas)
    """
    TSI = auto()
    CEU = auto()
    MKK = auto()
    YRI = auto()
    JPT = auto()
    CHB = auto()
    LWK = auto()
    CHD = auto()
    DEFAULT = auto()
    OUTROS = DEFAULT

    @classmethod
    def _missing_(cls, value=None):
        if value is None:
            return Population.DEFAULT
        elif isinstance(value, str):
            val = value.upper()
            for member in cls:
                if member.name == val:
                    return member

        raise KeyError(value)

    @classmethod
    def from_display_name(cls, name: str):
        """Recupera a população com base em seu nome de display (evita duplicatas)

        Args:
            name (str): Texto de Population.display_name

        Raises:
            ValueError: Se o texto não corresponder a nenhuma das populações.
        """
        for pop in cls:
            if pop.display_name == name:
                return pop

        raise ValueError(f"Population with display_name {name!r} not found")

    @property
    def bar_color(self):
        """Cor da barra do gráfico de barras."""
        DEFAULT_BAR_COLOR = "#F3B2A4"

        EUROPA_BAR_COLOR = "#92296D"
        AMERICA_SUL_BAR_COLOR = "#7D8CBD"
        AMERICA_NORTE_BAR_COLOR = "#65B1CF"
        AFRICA_BAR_COLOR = "#BD8DBD"
        ASIA_BAR_COLOR = "#EB5C61"
        bar_colors = {
            Population.TSI: EUROPA_BAR_COLOR,
            Population.CEU: EUROPA_BAR_COLOR,
            Population.MKK: AFRICA_BAR_COLOR,
            Population.YRI: AFRICA_BAR_COLOR,
            Population.JPT: ASIA_BAR_COLOR,
            Population.CHB: ASIA_BAR_COLOR,
            Population.LWK: AFRICA_BAR_COLOR,
            Population.CHD: ASIA_BAR_COLOR,
            Population.OUTROS: "#797474"
        }

        return bar_colors.get(self, DEFAULT_BAR_COLOR)

    @property
    def continente(self):
        """Continente usado na escolha da cor do pin e do texto do resumo."""
        continentes = {
            Continente.AMERICA_DO_SUL: (),
            Continente.AMERICA_DO_NORTE: (),
            Continente.EUROPA: (Population.TSI, Population.CEU),
            Continente.AFRICA: (Population.MKK, Population.YRI, Population.LWK),
            Continente.ASIA: (Population.CHB, Population.JPT, Population.CHD),
            Continente.OCEANIA: ()
        }
        for continente, pops in continentes.items():
            if self in pops:
                return continente

        raise KeyError(f"Continente não encontrado para {self}")

    @property
    def superpop_code(self):
        """Código da superpopulação, usado para selecionar os pesos dos SNPs."""
        codes = {
            SuperPopulation.AMR: (),
            SuperPopulation.EUR: (Population.TSI, Population.CEU),
            SuperPopulation.AFR: (Population.MKK, Population.YRI, Population.LWK),
            SuperPopulation.EAS: (Population.CHB, Population.JPT, Population.CHD),
            SuperPopulation.SAS: ()
        }
        for continente, pops in codes.items():
            if self in pops:
                return continente

        raise KeyError(f"Código de super população não encontrado para {self}")

    @property
    def pin_coordinates(self):
        """Coordenadas do pin no PDF em mm a partir de topleft."""
        pin_coordinates_dict = {
            Population.TSI: (111.0, 106.0),
            Population.MKK: (105.0, 125.0),
            Population.CEU: (102.5, 99.5),
            Population.CHB: (160.0, 103.0),
            Population.JPT: (180.0, 103.0),
            Population.YRI: (118.0, 130.0),
            Population.LWK: (118.0, 135.0),
            Population.CHD: (165.0, 110.0)
        }

        return pin_coordinates_dict[self]

    @property
    def pin_file(self):
        """Arquivo com a imagem do pin."""
        DEFAULT_PIN_FILE = "pin-demais.png"
        pin_files = {
            Continente.AMERICA_DO_SUL: "pin-america-sul.png",
            Continente.AMERICA_DO_NORTE: "pin-america-norte.png",
            Continente.EUROPA: "pin-europa.png",
            Continente.AFRICA: "pin-africa.png",
            Continente.ASIA: "pin-japao-oceania.png",
            Continente.OCEANIA: "pin-japao-oceania.png"
        }

        pin_filename = pin_files.get(self.continente, DEFAULT_PIN_FILE)
        return CONFIG.pins / pin_filename

    @property
    def display_name(self):
        """Nome da população para o gráfico de barras."""
        display_names = {
            Population.TSI: "Itália, Europa (TSI)",
            Population.CEU: "Europa do Norte e Ocidental (CEU)",
            Population.MKK: "África",
            Population.CHB: "China (CHB)",
            Population.JPT: "Japão (JPT)",
            Population.YRI: "África",
            Population.LWK: "África",
            Population.CHD: "Ásia (CHD)",
            Population.DEFAULT: "Outros"
        }
        return display_names[self]


class SNP:
    """Classe representando um SNP de uma pessoa.

    Atributos
    ---------
        rsid : str
            Código rsid do SNP
        genotipo : str
            Genótipo daquele SNP na forma dos dois alelos em ordem alfabética.
        gene: str | None
            Nome do gene que corresponde ao SNP.

    Métodos
    -------
        get_infile_data(filepath: PathLike[str] | str) -> dict[str, SNP]
            Método estático público. Lê um arquivo de genótipos.
        get_callback_from_line(line: str, delimiter: str) -> Callable[[str, str], tuple[str, str, str]]
            Método estático interno. Identifica a função de extração de dados do arquivo.
    """

    def __init__(self, rsid: str, a1: str, a2: str, gene: str = None):
        self.rsid = rsid
        self.genotipo = "".join(sorted(a1.upper()+a2.upper()))
        self.gene = gene

    def __hash__(self):
        return hash(self.rsid)

    def __eq__(self, other):
        if isinstance(other, str):
            return self.rsid == other
        elif isinstance(other, type(self)):
            return self.rsid == other.rsid

        raise NotImplementedError(
            f"__eq__ not implemented between SNP and {type(other).__name__}")

    @staticmethod
    def get_infile_data2(filepath: Union[PathLike[str], str]):
        """Retorna os dados do arquivo de genótipos.

        Args:
            filepath (Union[PathLike[str], str]): Caminho do arquivo de entrada.

        Returns:
            dict[str, SNP]: SNPs obtidos e seus genótipos na forma rsid: SNP.
        """
        file = Path(filepath)

        snps: dict[str, "SNP"] = {}
        with open(file, encoding="utf-8") as handle:
            header = handle.readline()
            delimiter = "\t" if "\t" in header else ";"

            for line in handle:
                splitted = line.strip("\n").split(delimiter)
                if len(splitted) == 3 and len(splitted[-1]) == 2:
                    rsid, _, (a1, a2) = splitted
                elif len(splitted) == 4 and len(splitted[-1]) == 2:
                    rsid, _, _, (a1, a2) = splitted
                elif len(splitted) == 5 and len(splitted[2]) == 1 and len(splitted[3]) == 1:
                    rsid, _, a1, a2, _ = splitted
                elif len(splitted) == 6 and len(splitted[3]) == 1 and len(splitted[4]) == 1:
                    rsid, _, _, a1, a2, _ = splitted
                else:
                    raise ValueError(
                        f"Falha na leitura da seguinte linha do arquivo {filepath}:", line)

                if rsid and DNA_REGEX.fullmatch(f"{a1}{a2}"):
                    snps[rsid] = SNP(rsid, a1, a2)

        return snps

    def get_callback_fromline(line: str, delimiter: str):
        """Retorna uma função de extração dos dados genotípicos,
        baseada na formatação aparente dos mesmos."""
        
        #TODO: Consertar todos os callbacks para não quebrar caso encontrem cromossomos X, Y ou MT 
        #(apenas um valor). Além disso, verificar se é mulher e recebeu apenas um valor em X
        #Checar na função que chama se é valor defeituoso e tratar; Retornar cromossomo sempre?
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
            if a2 == None and chrom in ["MT","X","Y"]:
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

    @classmethod
    def get_infile_data(cls, filepath: Union[PathLike[str], str]):
        """Retorna os dados do arquivo de genótipos.

        Args:
            filepath (Union[PathLike[str], str]): Caminho do arquivo de entrada.

        Returns:
            dict[str, SNP]: SNPs obtidos e seus genótipos na forma rsid: SNP.
        """
        file = Path(filepath)

        snps: dict[str, SNP] = {}
        with open(file, encoding="utf-8") as handle:
            header = handle.readline()
            delimiter = get_file_delimiter(file, encoding="utf-8")

            line = handle.readline()
            callback = cls.get_callback_fromline(line, delimiter)

            rsid, a1, a2 = callback(line, delimiter)
            if rsid and DNA_REGEX.fullmatch(f"{a1}{a2}"):
                snps[rsid] = SNP(rsid, a1, a2)  # type: ignore

            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                rsid, a1, a2 = callback(line, delimiter)

                if rsid and DNA_REGEX.fullmatch(f"{a1}{a2}"):
                    snps[rsid] = SNP(rsid, a1, a2)  # type: ignore

        return snps

    def __repr__(self):
        return f"SNP({self.rsid}, {self.genotipo})"


class Paciente:
    """
    Classe representando um paciente.

    Atributos
    ----------
        userid : str
            ID único do paciente, usado para identificar duplicatas.
        filepath : PathLike[str]
            Caminho do arquivo contendo os genótipos do paciente.
        nome : str
            Nome completo do paciente. (case sensitive)
        sexo : Sexo
            Sexo do paciente.
        super_population: SuperPopulation
            Código da superpopulação do paciente.
        data_nascimento: date
            Data de nascimento do paciente.
        documento: str
            Documento de identificação do paciente.
        names: list[str]
            Lista com o nome e sobrenomes do paciente.
        iniciais: str
            Iniciais do paciente. (apenas primeiro e último nomes)
        snps: dict[str, SNP]
            SNPs do paciente contendo rsid e genótipo.
        ancestralidades: list[Ancestralidade]
            Lista de ancestralidades do paciente. (cached)
        idade: int
            Idade do paciente em anos. (readonly)
    """

    def __init__(self, filepath: PathLike[str], sexo: Sexo = Sexo.INDEFINIDO, super_population: SuperPopulation = SuperPopulation.GLOBAL, data_nascimento: Optional[date] = None, documento: Optional[str] = None):
        """Construtor de Paciente.

        Args:
        -----
            userid (str): ID do paciente.
            filepath (PathLike[str]): Caminho completo do arquivo com os genótipos do paciente.
            nome (str): Nome completo do paciente (case sensitive).
            sexo (Sexo, optional): Sexo do paciente. Defaults to Sexo.INDEFINIDO.
            super_population (SuperPopulation, optional): Código da super população do paciente. Defaults to SuperPopulation.GLOBAL.
            data_nascimento (date, optional): Data de nascimento do paciente. Defaults to None.
            documento (str, optional): Documento (RG, CPF) do paciente. Defaults to None.

        Raises:
        -------
            ValueError: Se o paciente não possuir sobrenome (incapaz de definir iniciais).
            ValueError: Se o arquivo de entrada não for fornecido ou não existir.
        """
        self.filepath = Path(filepath)
        self.sexo = sexo
        self.super_population = super_population
        self.data_nascimento = data_nascimento
        self.documento = documento

        if not filepath:
            raise ValueError("Arquivo de entrada não fornecido")
        if not self.filepath.is_file():
            raise FileNotFoundError(
                f"Arquivo de entrada {str(self.filepath)!r} não encontrado")

        self.snps = SNP.get_infile_data(self.filepath)

    @cached_property
    def ancestralidades(self):
        """Retorna a lista com as populações ancestrais estimadas do paciente.
        Realiza o cálculo apenas uma vez, podendo ser sobrescrita.

        Returns:
            list[Ancestralidade]: Lista com ancestralidades
                contendo o código populacional e a porcentagem.
        """
        return Ancestralidade.calc_ancestralidades(self.snps, Path(self.filepath).stem)

    @property
    def idade(self):
        """Idade em anos do paciente."""
        if self.data_nascimento is None:
            raise ValueError(
                "Impossível calcular a idade sem uma data de nascimento.")
        return date.today().year - self.data_nascimento.year

    def format(self, text: str):
        txt = self.sexo.format(text).strip()
        if "[" in txt or "]" in txt:
            logger.error(
                f"Falha ao interpretar alguma bandeira do texto {txt!r}")
        return txt


class Ancestralidade:
    """Classe representando a porcentagem de ancestralidade de um paciente.

    Atributos
    ---------
        population : Population
            População de origem. (TSI, YRI, ...)
        value : Decimal
            Valor da porcentagem de ancestralidade entre 0 e 100

    Métodos
    --------
        get_highest_continente(ancestries: list[Ancestralidade]) -> tuple[Continente, Decimal]
            Retorna o continente com maior porcentagem relativa de ancestralidades e seu valor.
        calc_ancestralidades(snps: dict[str, SNP], filename: str) -> list[Ancestralidade]
            Obtém as porcentagens relativas das ancestralidades.
    """

    def __init__(self, population: Union[Population, None] = Population.DEFAULT, value: Decimal = None):
        if value is None:
            raise ValueError(
                "Algum valor deve ser fornecido para a ancestralidade")

        self.population = Population.DEFAULT if population is None else population
        self.value = value

    def __repr__(self):
        return f"{self.__class__.__name__}({self.population}, {self.value})"

    @staticmethod
    def get_highest_continente(ancestries: list["Ancestralidade"]):
        """Recupera o continente com maior porcentagem relativa de ancestralidades e seu valor."""
        if not ancestries:
            raise ValueError(
                "Pelo menos uma ancestralidade deve ser fornecida")

        continentes = {cont: Decimal(0) for cont in Continente}

        for ancestry in ancestries:
            continentes[ancestry.population.continente] += ancestry.value

        return sorted(continentes.items(), reverse=True, key=lambda pair: pair[1])[0]

    @staticmethod
    def calc_ancestralidades(snps: dict[str, SNP], filename: str):
        """Usa o IADMIX para calcular as porcentagens relativas de ancestralidade."""
        if not snps:
            raise ValueError(
                "Não há SNPs suficientes para calcular a ancestralidade")


        tmp_infile = get_temp_filename(
            root=CONFIG.iadmix_dir/"input", prefix=f"tmp-{filename}")
        tmp_outfile = f"{tmp_infile}.out"

        tmp_infile_path = f"{CONFIG.iadmix_dir}/input/{tmp_infile}"
        tmp_outfile_path = f"{CONFIG.iadmix_dir}/output/{tmp_outfile}"

        e = None
        try:
            with open(tmp_infile_path, "w", encoding="utf-8") as handle:
                for snp in snps.values():
                    handle.write(f"{snp.rsid}\t{snp.genotipo}\n")

            # run iadmix here
            ancestry_pipeline_genotypes(
                CONFIG.iadmix_freq_file, tmp_infile_path, tmp_outfile_path)

            with open(tmp_outfile_path, encoding="utf-8") as handle:
                lines = handle.readlines()

                if "FINAL_NZ_PROPS" not in lines[-1]:
                    raise ValueError("Falha ao ler o output do iadmix")

                *ancests_, _ = lines[-1].split()
                ancests: list[tuple[str, str]] = [
                    tuple(anc.split(":")) for anc in ancests_]
                ancestries = [
                    Ancestralidade(Population[pop], 100*Decimal(value))
                    for pop, value in ancests]
        finally:
            try:
                remove(tmp_infile_path)
                remove(tmp_outfile_path)
            except FileNotFoundError:
                pass

        if (soma := sum(anc.value for anc in ancestries)) > Decimal("100.01"):
            raise ValueError(
                f"As ancestralidades deram soma superior a 100.01% ({soma})")

        return ancestries


class CaracteristicaBase:
    """Classe reunindo as funções básicas de uma característica

    ...

    Atributos
    ---------
        id_ : str
            ID único da característica, usado para identificar duplicatas.
        paciente : Paciente
            Paciente que fornece os dados para a característica.
        caderno: str
            ID do caderno da característica.
            (sempre termina em 'dna', mas não impede a busca dos arquivos)
        is_z_score: bool
            Se a característica deve ou não ter seu score absoluto calculado pelo Z-Score.
            (impede o cálculo da freqência populacional do genótipo)
        total_snp_count : int
            Quantidade de SNPs total da característica.
        abs_score : float
            Valor do score absoluto entre 0 e 100.
        snps: dict[str, SNP]
            SNPs usados para o cálculo do score absoluto.
        resumo_var: str
            Resumo que aparece sob o item "Seus Resultados"
        resumo_var_card: str
            Resumo que aparece sob o item "Seus Resultados" nos cards de Resumo
        display_abs_score: float
            Valor exibido na barra do score absoluto. (de 0 a 100)
        qual_score: str
            Texto qualitativo exibido dentro da barra do score absoluto.
        abs_color : Color | None
            Se presente, força a cor da barra do score absoluto.
        pop_geno_user: Decimal | None
            Frequência populacional do genótipo do paciente.
            (apenas se não for Z-Score)
        pop_freqs: list[tuple[str, float]] | None
            Frequências populacionais dos genótipos.
            (mostra os 3 genótipos se houver apenas um SNP,
            ou apenas duas opções se mais de um,
            apenas se não for Z-Score)
        index_title: str
            Título no índice da característica.
            (permite o agrupamento no índice)
        confiabilidade: int
            Nível de confiabilidade dos dados da característica (4 ou 5)
        resumo_const: str
            Resumo que aparece sob o item "Seus Genótipos"
        icone : PathLike[str]
            Caminho do arquivo com o ícone da caraterística
        referencias : list[str]
            Referências bibliográficas da característica.

    Métodos
    ---------
        get_dados_var()
            Calcula, se necessário, os dados variáveis de settings.VARDATA_FILE
        get_resumo_var()
            Retorna um resumo que depende do score absoluto.
        get_abs_score(snps: dict[str, SNP], pop: SuperPopulation) -> tuple[list[str], float]
            Retorna os SNPs usados no cálculo e o valor do score absoluto em percentagem.
        get_abs_z_score(...)
            Calcula o score absoluto a partir do Z-Score.
        get_z_score(...)
            Calcula o Z-Score.
        get_abs_norm_score(...)
            Calcula o score absoluto a partir dos dados curados.
        get_restricoes(...)
            Obtém as restrições da característica.
            (por enquanto apenas Sexo do paciente)
        restrict(paciente: Paciente) -> bool
            Define se a característica deve ser impedida conforme as restrições.
        get_pop_freqs(...) -> tuple[None, None] | tuple[Decimal, list[tuple[str, float]]]
            Obtém os dados de frequência populacional dos genótipos.
            (usa apenas os dados curados por enquanto)
        get_dados_const() -> tuple[str, int, str]
            Obtém os dados constantes de settings.CARACTERISTICAS_FILE
    """

    def __init__(self, id_: str, paciente: Paciente, caderno: str, is_z_score: bool):
        self.id_ = id_.lower()
        self.caderno = Caderno.make_valid_id(caderno)
        self.paciente = paciente
        self.is_zscore = is_z_score

        self.restricoes: dict[str, Sexo] = self.get_restricoes(
            self.caderno, id_)

        if self.caderno not in CONFIG.cadernos_validos:
            raise ValueError(f"Caderno {self.caderno!r} inválido")

        if self.restrict(paciente):
            raise RestrictionError(
                f"Restrição de {self.id_!r} impediu o uso dessa característica")

        self.total_snp_count, snp_list, self.abs_score = self.get_abs_score(
            paciente.snps, paciente.super_population)
        self.snp_count = len(snp_list)

        self.snps = {rsid: paciente.snps[rsid] for rsid in snp_list}

        self.resumo_var, self.resumo_var_card = self.get_resumo_var()
        self.display_abs_score, self.qual_score, self.abs_color = self.get_dados_var()
        self.index_title, self.confiabilidade, self.resumo_const, self.referencias = self.get_dados_const()
        self.icone = self.get_icone_file()

        self.pop_geno_user, self.pop_freqs = self.get_pop_freqs(
            self.snps, paciente.super_population)

        total_s = "s" if self.total_snp_count != 1 else ""
        s = "s" if self.snp_count != 1 else ""
        msg = f"{self.id_}: {self.snp_count} de {self.total_snp_count} SNP{total_s} encontrado{s}"

        if self.snp_count != self.total_snp_count:
            logger.info(msg)
        else:
            logger.debug(msg)

    def __hash__(self):
        return hash(self.id_)

    def get_abs_score(self, snps: dict[str, SNP], pop: SuperPopulation):
        """Retorna o score absoluto (entre 0 e 100)"""
        total_count, seen, score = self.get_abs_z_score(
            snps, pop) if self.is_zscore else self.get_abs_norm_score(snps, pop)

        return total_count, seen, score

    def get_abs_z_score(self, snps: dict[str, SNP], superpop: SuperPopulation):
        """Calcula o score absoluto a partir do z-score."""
        snp_count, rsids, z_score = self.get_z_score(snps, superpop)

        abs_score = 100*cdf(z_score)
        return snp_count, rsids, abs_score

    def get_z_score(self, snps: dict[str, SNP], superpop: SuperPopulation):
        """Calcula o Z-Score da característica."""
        # O arquivo já deve ser válido
        pop_file = CONFIG.scores / superpop.zscore_filename

        if not pop_file.is_file():
            raise ValueError(
                f"Arquivo {superpop.zscore_filename} não encontrado")

        df = pd.read_table(pop_file, index_col=0,
                           decimal=".", encoding="utf-8")

        df: pd.DataFrame = df.loc[self.id_]  # type: ignore
        filtered_df = df[
            ~df.select_dtypes("object")
            .fillna('')
            .astype(str)
            .add('|')
            .sum(axis=1)
            .str.contains('?', regex=False)
        ]  # remove os SNPs zoados

        rsids: list[str] = []
        missing_effect_info_snps = []
        missing_genotype_snps = []
        pop_sd_snps = []
        scores = []
        snp_count = len(filtered_df)

        for _, snp in filtered_df.iterrows():  # type: ignore
            snp: ZScoreSnpData
            rsid = snp["rsid"]
            eff_allele, _, neff_allele = snp["effect_allele"].upper()
            _, _, min_allele = snp["major_allele"].upper()
            min_allele_freq = snp["MAF"]
            effect_size = snp["effect_size"]

            if rsid not in snps:
                missing_genotype_snps.append(rsid)
                continue
            elif not set(snps[rsid].genotipo).issubset({eff_allele, neff_allele}):
                missing_effect_info_snps.append(rsid)
                continue

            rsids.append(rsid)
            eff_allele_freq = min_allele_freq if eff_allele == min_allele else 1-min_allele_freq
            effect_allele_count = snps[rsid].genotipo.count(eff_allele)

            frac_0 = (1-eff_allele_freq)**2
            frac_2 = eff_allele_freq**2
            frac_1 = 1-frac_0-frac_2

            mean = effect_size*(frac_1 + 2*frac_2)
            sigma = frac_0*mean**2 + frac_1 * \
                (effect_size-mean)**2 + frac_2*(2*effect_size-mean)**2
            pop_sd_snps.append(sigma)

            score = effect_size*effect_allele_count - 2*eff_allele_freq*effect_size
            scores.append(score)

        pop_sd = sqrt(sum(pop_sd_snps))

        if missing_effect_info_snps:
            logger.warning(
                f"SNPs no arquivo de entrada que não tem equivalência de alelos no arquivo de config ({len(missing_effect_info_snps)})")
            # logger.warn(*missing_effect_info_snps, sep=", ")
        if missing_genotype_snps:
            logger.warning(
                f"SNPs ausentes no arquivo de entrada ({len(missing_genotype_snps)})")
            # logger.warn(*missing_genotype_snps, sep=', ')

        if not scores:
            raise ValueError(
                f"Não há SNPs suficientes para calcular o Z-score de {self.id_!r}")

        z_score = round(sum(scores)/pop_sd, 4)

        return snp_count, rsids, z_score

    def restrict(self, paciente: Paciente) -> bool:
        """Se a característica não pode ser calculada conforme o paciente."""
        should_stop = False

        rsexo = self.restricoes.get("sexo")
        if rsexo is not None and rsexo not in paciente.sexo:
            should_stop = True

        return should_stop

    @staticmethod
    def get_restricoes(caderno: str, caracteristica: str):
        """Retorna as restrições da característica."""
        out_rests = {}
        for line in CONFIG.cadernos_data:
            cad_id, char_id, restrict = line["caderno"].lower(
            ), line["caract"].lower(), line["restrict"].lower()

            if cad_id == caderno and char_id == caracteristica and restrict:
                in_rests = restrict.split(";")
                for rest in in_rests:
                    key, val = rest.split("=")
                    if key == "sexo":
                        out_rests["sexo"] = Sexo[val.upper()]
        return out_rests

    def get_dados_var(self):
        """Retorna o valor quantitativo e qualitativo do score absoluto.
        Esses valores não dependem do valor individual de cada SNP nem é ciente dos mesmos."""
        for line in CONFIG.vardata_data:
            char_id, regra = line["id da característica"], line["regra"]
            quant, qual, color = line["valor da porcentagem"], line["valor qualitativo"], line["cor da barra"]

            if char_id.lower() != self.id_:
                continue

            if color:
                try:
                    color_ = HexColor(color)
                except ValueError as e:
                    logger.error(
                        f"Falha na leitura da cor da barra da característica {self.id_}: {e.args[0]}")
                    color_ = None
            else:
                color_ = None

            if eval_rule(regra, self.abs_score):
                return float(eval_rule(quant, self.abs_score)), qual, color_

        # default values
        # logger.debug(f"Defaulting values")
        qual = CONFIG.abs_chart_display_text(self.abs_score)
        display = CONFIG.abs_chart_display_value(self.abs_score)
        return display, qual, None

    def get_resumo_var(self):
        """Retorna resumos que dependem do score absoluto.
        Esse resumo não depende do valor individual de cada SNP nem é ciente dos mesmos."""
        for line in CONFIG.resumos_data:
            char_id, regra, resumo = line["id da característica"], line["regra"], line["resumo"]
            resumo_card = line['resumo card'].strip() or resumo
            if char_id.lower() != self.id_:
                continue

            if eval_rule(regra, self.abs_score):
                return self.paciente.format(resumo), self.paciente.format(resumo_card)

        raise ValueError(
            f"Resumo da característica {self.id_!r} não definido quando score vale {self.abs_score:.2f}")

    def get_pop_freqs(self, snps: dict[str, SNP], superpop: SuperPopulation):
        """Indica a qual grupo o paciente pertence, e separa em porcentagem os grupos genotípicos.
        Assume que as freqs. genotípicas são independentes."""
        if self.is_zscore:
            return None, None

        pop_file = CONFIG.scores / superpop.normscore_filename

        freqs_geno: dict[str, list[tuple[str, float]]] = {}

        for line in csv_reader(pop_file, lower_keys=True):
            caract_id, rsid = line['caract'], line['rsid']
            eff, _, non_eff = line['effect_allele'].upper()
            f1, f2, f3 = map(Decimal, line['frequencies'].split("/"))

            if caract_id.lower() != self.id_ or rsid not in snps:
                continue

            EE = eff*2
            NE = "".join(sorted(eff+non_eff))
            NN = non_eff*2

            freqs_geno[rsid] = [
                (EE, float(f1)), (NE, float(f2)), (NN, float(f3))]

        if not freqs_geno:
            raise ValueError(
                f"Não há SNPs suficientes para calcular as frequências genotípicas de {self.id_!r}")

        if len(freqs_geno) == 1:
            rsid, data = freqs_geno.popitem()
            value = [val for geno, val in data if geno ==
                     snps[rsid].genotipo][0]
            graph_data = [(f"{round(val)}%<br/>{geno}", val)
                          for geno, val in data]
            return Decimal(value), graph_data

        else:
            # logger.debug(freqs_geno)
            vals: list[float] = []
            for rsid, data in freqs_geno.items():
                val = next(float(val) for geno, val in data if geno ==
                           snps[rsid].genotipo) / 100
                vals.append(val)

            value = 100*prod(vals)
            graph_data: list[tuple[str, float]] = [
                (f"{round(value)}%<br/>Seu genótipo", value), ("", 100-value)]
            return Decimal(round(value)), graph_data

    def get_dados_const(self):
        """Retorna os dados que não dependem do score da característica:
            - título no índice;
            - confiabilidade;
            - resumo constante;
            - referencia.
        """
        for line in CONFIG.caracteristicas_data:
            char_id, index_title = line['característica'], line['título no índice']
            resumo = line['resumo constante']
            referencias = [ref.strip()
                           for ref in line["referências"].split("|") if ref.strip()]

            try:
                confiabilidade = int(line['confiabilidade'])
            except ValueError as e:
                logger.critical(f"Offending line: {line!r}")
                raise ValueError(
                    "Falha ao obter o valor da confiabilidade") from e

            if char_id == self.id_:
                return index_title, confiabilidade, resumo, referencias

        raise ValueError(
            f"Dados constantes da característica {self.id_} não encontrados")

    def get_icone_file(self):
        """Retorna o caminho do arquivo com o ícone da característica."""
        file = CONFIG.icones / f"{self.id_}.png"
        default = CONFIG.icones / "missing.png"

        if not file.is_file():
            logger.warning(f"Arquivo com o ícone de {self.id_} não encontrado")
            return default

        return file


class Caracteristica(CaracteristicaBase):
    """
    Classe representando uma característica qualquer
    (diabetes, Intolerância ao glutén, resistência atlética...).

    ...

    Atributos (não herdados)
    ---------
        template: str
            Caminho do arquivo com o modelo da característica.

    Métodos (não herdados)
    --------
        get_template_file()
            Retorna o caminho do arquivo com o template da característica.
        get_efeito_file(...) -> PathLike[str] | None
            Obtém o caminho dos desenhos do efeito do SNP (as 3 bolinhas)
    """

    def __init__(self, id_: str, paciente: Paciente, caderno: str, is_z_score: bool):
        super().__init__(id_, paciente, caderno, is_z_score)

        # logger.debug(repr(self))

        self.template = self.get_template_file()

    def __repr__(self):
        return f"Caracteristica({self.id_}, {self.paciente!r}, {self.is_zscore})"

    def __str__(self):
        return f"Característica {self.id_!r} com {self.snp_count}/{self.total_snp_count} SNPs"

    def get_template_file(self):
        """Retorna o caminho do arquivo com o template da característica."""
        temp1 = CONFIG.templates / self.caderno.upper() / f"{self.id_}.pdf"
        temp2 = CONFIG.templates / self.caderno.upper().removesuffix("DNA") / \
            f"{self.id_}.pdf"

        if temp1.is_file():
            return temp1
        elif temp2.is_file():
            return temp2

        display_name = Caderno.get_display_name(self.caderno)
        raise ValueError(
            f"Template de {self.id_!r} não encontrado no caderno {display_name!r}")

    def get_efeito_file(self, snp: SNP, color: Color):
        """Retorna o caminho do arquivo com a imagem do efeito de um dado SNP"""
        if snp.genotipo == "--" or self.is_zscore:
            return

        pop_file = CONFIG.scores / self.paciente.super_population.normscore_filename

        for line in csv_reader(pop_file, lower_keys=True):
            caract_id, rsid = line['caract'], line['rsid']
            eff_allele, begin = line['effect_allele'], line['scale_begin']
            if caract_id.lower() == self.id_ and rsid == snp.rsid:
                break
        else:
            raise ValueError(f"Falha ao obter o efeito de {snp}")
        effect = snp.genotipo.count(eff_allele[0].upper())
        level = truncate(effect + int(begin), 0, 3)

        DEFAULT_FILE = CONFIG.efeitos / \
            f"{self.caderno.upper()}-{level}.png"

        # hexcode = color.hexval().upper()
        # filepath = MEDIA_EFEITOS_DIR / f"{hexcode}-{level}.png"
        # if filepath.is_file():
        #     return filepath

        # logger.warn(f"Falha ao obter o arquivo de efeito, usando o padrão do caderno {self.caderno}")
        # logger.info(f"{hexcode=}, {snp=}, {level=}\n")
        return DEFAULT_FILE


class FunCaracteristica(CaracteristicaBase):
    """
    Classe representando uma característica divertida qualquer
    (consumo de refrigerante, ...).
    """

    def __init__(self, id_: str, paciente: Paciente, caderno: str = "fundna", is_z_score: bool = False):
        super().__init__(id_, paciente, "fundna", False)

    def __repr__(self):
        return f"Caracteristica({self.id_}, {self.paciente!r}, {self.is_zscore})"


class Caderno:
    """Classe representando um caderno qualquer
    (nutri, beauty, mint...).

    ...

    Atributos
    ---------
        caracteristicas : list[Caracteristica]
            Características que compoem o caderno, já na ordem.
        id_ : str
            ID do caderno.
        caderno : str
            Nome do caderno.
        paciente : Paciente
            Paciente a quem o caderno se refere.
        template_index : Path
            Caminho do arquivo de template do índice do caderno.
        extra_pages : tuple[list[Path], ...]
            Caminhos dos arquivos (uma página cada idealmente):
            antes do índice, entre o índice e as características, e após as características.

    Métodos
    --------
        from_id(caderno: str, paciente: Paciente, verbose=False) -> tuple[Caderno, list[str]]
            Factory method. Cria um caderno com todas as características e lista as características faltantes.
        append(char: Caracteristica) -> None
            Método público. Adiciona uma característica ao final.
        insert(idx: int, char: Caracteristica) -> None
            Método público. Adiciona uma característica na posição idx.
        extend(chars: Iterable[Caracteristica]) -> None
            Método público. Adiciona as características ao final.
        remove(char: Caracteristica | str) -> None
            Método público. Remove a característica.
        clear() -> None
            Método público. Remove todas as características.
        make_valid_id(caderno: str) -> str
            Método estático público. Converte em um ID de caderno válido.
        get_display_name(caderno: str) -> str
            Método estático público. Retorna o nome do caderno.
    """

    def __init__(self, caracteristicas: list[Union[Caracteristica, FunCaracteristica]], caderno: str, paciente: Paciente):
        self.caracteristicas = caracteristicas
        self.id_ = self.make_valid_id(caderno)
        self.caderno = self.get_display_name(self.id_)
        self.paciente = paciente

        self_dir = list(CONFIG.templates.glob(
            f"{self.id_.upper().removesuffix('DNA')}*"))[0]
        templates = list(self_dir.glob("Índice.pdf"))
        if len(templates) == 1:
            self.template_index = templates[0]
        elif len(templates) > 1:
            raise ValueError(
                f"Há mais de um template do índice disponível para o caderno {self.caderno}: {templates}")
        else:
            raise ValueError(
                f"Template do índice não encontrado para o caderno {self.caderno}")

        pages_before = sorted(self_dir.glob(
            "Pré-Índice/*"))
        pages_after = sorted(self_dir.glob(
            "Pós-Índice/*"))
        pages_after_after = sorted(self_dir.glob(
            "Pós-Características/*"))
        self.extra_pages = (pages_before, pages_after, pages_after_after)

        if not self.caracteristicas:
            return

        if any(char.caderno != self.id_ for char in self.caracteristicas):
            raise ValueError(
                f"Alguma característica possui caderno diferente do esperado ({self.caderno!r})")
        if any(char.paciente != self.paciente for char in self.caracteristicas):
            raise ValueError(
                f"Alguma característica possui paciente diferente do esperado ({self.paciente})")

    @classmethod
    def from_id(cls, caderno: str, paciente: Paciente):
        """Cria um caderno com o máximo possível de características na ordem do arquivo settings.CADERNOS_FILE

        Args:
            caderno (str): ID do caderno a ser criado.
            paciente (Paciente): Paciente cujos dados permitem a criação do caderno.

        Raises:
            ValueError: Se o paciente não tiver SNPs
                (não é possível em teoria pois indica um erro anterior)

        Returns:
            tupl[Caderno, list[str]]: Caderno gerado e uma lista com os ids das características faltantes.
        """
        if not paciente.snps:
            raise ValueError("NENHUM SNP ENCONTRADO")

        caderno = cls.make_valid_id(caderno)

        char_type = FunCaracteristica if cls.make_valid_id(
            "fun") == caderno else Caracteristica

        chars: list[Union[Caracteristica, FunCaracteristica]] = []
        true_missing_chars: list[str] = []
        restricted_missing_chars: list[str] = []
        for line in CONFIG.cadernos_data:
            cad_id, char_id, is_z_score = line["caderno"].lower(
            ), line["caract"].lower(), line["is_zscore"]
            if cad_id == caderno:
                try:
                    char = char_type(
                        char_id, paciente, caderno, is_z_score=bool(int(is_z_score)))
                except NoSNPFoundError as e:
                    true_missing_chars.append(char_id)
                    logger.warning(f"{char_id}: {e}")
                except RestrictionError as e:
                    restricted_missing_chars.append(char_id)
                    logger.warning(f"{char_id}: {e}")
                else:
                    chars.append(char)

        return cls(chars, caderno, paciente), true_missing_chars, restricted_missing_chars

    def insert(self, idx: int, char: Caracteristica) -> None:
        """Adiciona a característica no índice da lista de características (0-based)."""
        return self.caracteristicas.insert(idx, char)

    def append(self, char: Caracteristica) -> None:
        """Adiciona a característica ao final da lista de características."""
        return self.caracteristicas.append(char)

    def extend(self, chars: Iterable[Caracteristica]) -> None:
        """Adiciona características ao final da lista de características."""
        return self.caracteristicas.extend(chars)

    @singledispatchmethod
    def remove(self, char):
        """Remove a característica"""
        raise NotImplementedError("Cannot remove unknown char")

    @remove.register
    def _(self, char: str):
        for caract in self.caracteristicas:
            if caract.id_ == char:
                self.caracteristicas.remove(caract)
                return

        raise ValueError(f"Failed to remove {char}")

    @remove.register
    def _(self, char: Caracteristica) -> None:
        """Remove a característica da lista de características."""
        return self.caracteristicas.remove(char)

    def clear(self) -> None:
        """Esvazia a lista de características."""
        return self.caracteristicas.clear()

    def __len__(self):
        return len(self.caracteristicas)

    def __bool__(self):
        return bool(self.caracteristicas)

    def __repr__(self):
        s = 's' if len(self.caracteristicas) != 1 else ''
        return f"Caderno {self.caderno} com {len(self.caracteristicas)} característica{s}"

    @staticmethod
    def make_valid_id(caderno: str):
        """Tenta converter o valor em um ID válido."""
        caderno = caderno.lower().removesuffix("dna") + "dna"
        if caderno not in CONFIG.cadernos_validos:
            raise ValueError(f"ID de caderno {caderno!r} inválido")

        return caderno

    @staticmethod
    def get_display_name(caderno_id: str):
        """Recupera o nome do caderno a partir do ID."""
        for line in CONFIG.cadernos_data:
            if line['caderno'].lower() == caderno_id:
                return line['caderno']

        raise ValueError(f"Caderno {caderno_id!r} não encontrado")