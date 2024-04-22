# -*- coding: utf-8 -*-
# packages/helpers.py

"""Módulo que contém funções auxiliares sem relação com os dados"""

import sys
import argparse
import csv
import logging
import operator
import os
import re
from os import PathLike
from pathlib import Path
from typing import Any, Callable, Type, TypeVar, Union

from reportlab.lib.utils import ImageReader
from reportlab.platypus.flowables import Image
from typing_extensions import TypeGuard

logger = logging.getLogger("helpers")

# This corrects the type signature for pylance
_TFunc = TypeVar("_TFunc", bound=Callable[..., Any])


def get_file_encoding(filepath: Union[PathLike[str], str]):
    """Chooses a possible encoding for .csv files
    Only between UTF-8 and Latin-1 (ISO 8859-1)

    Args:
        filepath (str): .csv filepath.

    Returns:
        str: Encoding found.
    """
    for encoding in ('utf-8', 'latin-1'):
        try:
            with open(filepath, "r", newline="", encoding=encoding) as handle:
                handle.read()
            return encoding
        except UnicodeDecodeError:
            pass
    raise ValueError(f"Failed to identify encoding for file '{filepath}'")


def get_file_delimiter(filepath: Union[PathLike[str], str], encoding: str = None):
    """Guess .csv file delimiter.

    Args:
        filepath (PathLike): .csv filepath.
        encoding (str, Optional): .csv file encoding. Defaults to None.

    Returns:
        str: Delimiter chosen between "\\t", "," and ";".
    """
    if encoding is None:
        encoding = get_file_encoding(filepath)

    delims = ["\t", ",", ";"]
    with open(filepath, "r", newline="", encoding=encoding) as handle:
        sniffer = handle.readline().rstrip("\n")
        by_order = sorted(
            delims, key=sniffer.count,
            reverse=True)

        return by_order[0]


def csv_reader(filepath: Path, remove_trailing_whitespace=True, lower_keys=False, lower_values=False, guess_type=False):
    """Gerador DictReader para arquivos.

    Args:
        filepath (str): Filepath.
        remove_trailing_whitespace (bool, optional): Se o arquivo deve ser normalizado.
            Defaults to True.
        lower_keys (bool, optional): Se os cabeçalhos devem ser normalizados.
            Defaults to False.
        lower_values (bool, optional): Se os valores devem ser normalizados.
            Defaults to False.

    Yields:
        dict: Linhas do arquivo.
    """
    encod = get_file_encoding(filepath)
    delim = get_file_delimiter(filepath, encod)
    with open(filepath, "r", newline="", encoding=encod) as handle:
        reader = csv.DictReader(handle, delimiter=delim)
        for line in reader:
            if not line:
                continue
            elif not any([v and str(v).strip() for v in line.values()]):
                continue

            line = {
                str(key).strip(): "" if value is None else str(value)
                for key, value in line.items()
                if key and str(key).strip()
            }

            if remove_trailing_whitespace:
                line = {
                    key: value.strip()
                    for key, value in line.items()
                }

            if lower_keys:
                line = {
                    key.lower(): value
                    for key, value in line.items()
                }

            if lower_values:
                line = {
                    key: str(value).lower()
                    for key, value in line.items()
                }

            yield line


def truncate(n: int, min_: int, max_: int):
    """Returns a value between min_ and max_ (inclusive)"""
    return max(min(n, max_), min_)


def get_temp_filename(root: Path = None, prefix: str = "tmp"):
    """Create a filename that doesn't exist."""
    if root is None:
        root = Path(__file__).resolve().parent

    root = Path(root)

    if not (root / f"{prefix}.txt").is_file():
        return f"{prefix}.txt"

    i = 0
    while (tmpfile := root / f"{prefix}-{i}.txt").is_file():
        i += 1

    return tmpfile.name


def get_image(path: Union[PathLike[str], str], width: float = None, height: float = None, **kwargs):
    """Retorna uma imagem corretamente escalonada a partir do arquivo path.
    Se apenas um dentre width e height for fornecido,
    a imagem é escalonada visando manter sua proporção original.

    Args:
        path (PathLike[str]): Caminho do arquivo da imagem (não lembro quais formatos, provv jpg, png).
        width (float, optional): Largura da imagem em unidades de PDF. Defaults to None.
        height (float, optional): Altura da imagem em unidades de PDF. Defaults to None.

    Raises:
        FileNotFoundError: Se a imagem não existe.

    Returns:
        Image: Imagem gerada a partir do arquivo.
    """
    if width is not None and height is None:
        img = ImageReader(str(path))
        iw, ih = img.getSize()
        aspect = ih/float(iw)
        return Image(path, width=width, height=(width*aspect), **kwargs)
    elif width is None and height is not None:
        img = ImageReader(str(path))
        iw, ih = img.getSize()
        aspect = iw/float(ih)
        return Image(path, width=(height*aspect), height=height, **kwargs)
    else:
        return Image(path, width=width, height=height, **kwargs)


_T = TypeVar("_T")


def is_list_of(val: list[Any], type_: Type[_T]) -> TypeGuard[list[_T]]:
    """Confirma se val é uma lista cujos items são todos do tipo type_
    """
    return all(isinstance(element, type_) for element in val)


ops: dict[str, Callable[[Any, Any], bool]] = {
    "<": operator.lt,
    ">": operator.gt,
    ">=": operator.ge,
    "<=": operator.le,
    "==": operator.eq,
    "!=": operator.ne,
    "=": operator.eq,
    "<>": operator.ne,
    "+": operator.add,
    "-": operator.sub,
    "*": operator.mul,
    "/": operator.truediv,
}

comps: dict[str, Callable[[Any, Any], bool]] = {
    "AND": operator.and_,
    "E": operator.and_,
    "OR": operator.or_,
    "OU": operator.or_
}


def eval_rule(rule: str, value: float):
    """Calcula a regra conforme pedido.

    Args:
        rule (str): Regra que descreve uma sequência de operações matemáticas.
        value (float): Valor alimentado para a variável MANTER dentro da regra.

    Raises:
        ValueError: Se a função falhar em resolver a regra para apenas um token.

    Returns:
        str | float | bool: Valor final após o cálculo da regra.
    """
    rule = rule.upper().replace(" ", "")

    if any(rule.startswith(op) for op in ops):
        expression = f"{value}{rule}".replace("MANTER", str(value))
    else:
        expression = rule.replace("MANTER", str(value))
    tokens: list[Union[str, float, bool]] = [
        tk for tk in re.split(r"(?: |([-<>=+*/]+|AND|E|OR|OU))", expression)
        if tk is not None]
    for idx, tk in enumerate(tokens):
        try:
            flt = float(tk)
        except ValueError:
            flt = tk

        tokens[idx] = flt

    while any(tk in ops for tk in tokens):
        for idx, tk in enumerate(tokens):
            if tk in ops:
                try:
                    booler = ops[str(tk)](tokens[idx-1], tokens[idx+1])
                except TypeError:
                    logger.critical(f"Failed to parse {rule=}, {tokens=}")
                    raise
                tokens[idx-1:idx+2] = [booler]
                break
    while any(tk in comps for tk in tokens):
        for idx, tk in enumerate(tokens):
            if tk in comps:
                booler = comps[str(tk)](tokens[idx-1], tokens[idx+1])
                tokens[idx-1:idx+2] = [booler]
                break

    if len(tokens) == 1:
        return tokens[0]

    raise ValueError(f"Something went wrong while evaluating rule: {tokens}")


def extant_file(x: str):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.isfile(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError(f"{x} does not exist")
    return Path(x)


def set_default_subparser(self: argparse.ArgumentParser, subparser_name: str, recv_args: list[str] = None):
    """
    see http://stackoverflow.com/questions/5176691/argparse-how-to-specify-a-default-subcommand
    """
    # pylint: disable=protected-access

    subparser_found = False
    if not sys.argv[1:]:
        sys.argv.append("--help")
        return

    for arg in sys.argv[1:]:
        if arg in ['-h', '--help', '--version', '-v']:
            # global help or version if no subparser
            break
    else:
        if self._subparsers is None:
            return
        for x in self._subparsers._actions:
            if not isinstance(x, argparse._SubParsersAction):
                continue
            for sp_name in x._name_parser_map.keys():
                if sp_name in sys.argv[1:]:
                    subparser_found = True
        if not subparser_found:
            # insert default in first position, this implies no
            # global options without a sub_parsers specified
            if recv_args is None:
                sys.argv.insert(1, subparser_name)
            else:
                recv_args.insert(0, subparser_name)


def filter_both(lst: list[_T], pred: Callable[[_T], bool] = None) -> tuple[list[_T], list[_T]]:
    if pred is None:
        pred = bool
    
    truy, falsy = [], []
    for val in lst:
        if pred(val):
            truy.append(val)
        else:
            falsy.append(val)
    return truy, falsy


def partition(lst: list[_T], n: int) -> list[list[_T]]:
    """Separates a list equally in n sublists, with bias toward first item."""
    pcholder = [None]*len(lst)
    division = len(lst)/n
    scaffold = [pcholder[round(division*i):round(division*(i+1))]
                for i in range(n)]
    scaffold.sort(key=len, reverse=True)
    idx = 0
    out: list[list] = []
    for col in scaffold:
        out.append([])
        for _ in col:
            out[-1].append(lst[idx])
            idx += 1
    return out
