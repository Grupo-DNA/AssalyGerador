import logging
import logging.config
import platform
import pprint
import re
import os
import pandas as pd
from functools import cached_property
from pathlib import Path
from typing import Any, TypedDict, Union, cast
from typing_extensions import TypeGuard

import yaml
from reportlab.lib.colors import Color, HexColor, toColorOrNone
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT, TA_RIGHT
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.pdfbase import pdfmetrics, ttfonts

if __name__ == "packages.config":
    from Controller.controller.exceptions import ConfigError
else:
    from Controller.controller.exceptions import ConfigError


_user_configpath = Path(__file__).resolve().parent.parent / "config.yaml"
_global_configpath = Path(__file__).resolve().parent / "DataFiles" /"Files"/ "config.yaml"
_default_configpath = Path(__file__).resolve().parent.parent / "DataFiles" / "Files"/ "Defaults" / "config.yaml"
_DIRNAME = Path(__file__).resolve().parent.parent


logger = logging.getLogger("config")

'''
def constructor_path(loader, node):
    """
    paths:
        root: &ROOT !path [path, to, root]
        patha: !path [*ROOT, a]
        pathb: !path [*ROOT, b]
    {paths:
      {root: /full/path/to/root,
       patha: /full/path/to/root/a,
       pathb: /full/path/to/root/b}
    }
    """
    seq = loader.construct_sequence(node)
    out_seq: list[str] = []
    for i in seq:  # this might deal with yaml shenanigans that I've no seen yet
        if isinstance(i, list):
            out_seq.extend(map(str, i))
        else:
            out_seq.append(str(i))

    out = Path(*out_seq)
    if out.is_absolute() and out.exists():
        return out
    nout = _DIRNAME / out
    if nout.is_absolute() and nout.exists():
        return nout
    else:
        raise FileNotFoundError(out)
        '''


def constructor_color(loader, node):
    txt = str(loader.construct_scalar(node))
    try:
        return HexColor(f"#{txt}")
    except ValueError:
        return toColorOrNone(txt)


def constructor_flatten(loader, node):
    seq = loader.construct_sequence(node)
    nseq = []
    for sub in seq:
        if isinstance(sub, list):
            nseq.extend(sub)
        else:
            nseq.append(sub)
    return nseq


def constructor_paragraph_style(loader, node):
    values: dict[str, Any] = loader.construct_mapping(node)
    name = values.pop("name")
    return ParagraphStyle(name, **values)


def constructor_text_alignment(loader, node) -> int:
    value: str = loader.construct_scalar(node)
    return {
        "CENTER": TA_CENTER,
        "RIGHT": TA_RIGHT,
        "LEFT": TA_LEFT,
        "JUSTIFY": TA_JUSTIFY
    }[value.upper()]

def constructor_path(loader, node):
    """
    paths:
        root: &ROOT !path [path, to, root]
        patha: !path [*ROOT, a]
        pathb: !path [*ROOT, b]
    {paths:
      {root: /full/path/to/root,
       patha: /full/path/to/root/a,
       pathb: /full/path/to/root/b}
    }
    """
    
    seq = loader.construct_sequence(node)
    out_seq: list[str] = []
    for i in seq:  # this might deal with yaml shenanigans that I've not seen yet
        if isinstance(i, list):
            out_seq.extend(map(str, i))
        else:
            out_seq.append(str(i))
    out = Path(*out_seq)
    print('CAMINHO:',out)
    if out.is_absolute() and out.exists():
        return out
    nout = _DIRNAME / out
    print(_DIRNAME)
    if nout.is_absolute() and nout.exists():
        return nout
    else:
        raise FileNotFoundError(out)

color_pattern = re.compile(r"^#?[a-fA-F0-9]{6}$")
yaml.SafeLoader.add_implicit_resolver("!color", color_pattern, None)
yaml.SafeLoader.add_constructor("!path", constructor_path)
yaml.SafeLoader.add_constructor("!color", constructor_color)
yaml.SafeLoader.add_constructor("!text_style", constructor_paragraph_style)
yaml.SafeLoader.add_constructor("!alignment", constructor_text_alignment)
yaml.SafeLoader.add_constructor("!aln", constructor_text_alignment)
yaml.SafeLoader.add_constructor("!flatten", constructor_flatten)


def get_config_data():
    with open(_default_configpath) as handle:
        def_data = handle.read()

    if _user_configpath.is_file():
        print("user config file found")
        with open(_user_configpath) as handle:
            usr_data = handle.read()
    elif _global_configpath.is_file():
        print("global config file found")
        with open(_global_configpath) as handle:
            usr_data = handle.read()
    else:
        usr_data = ""

    # this is so references in usr files
    # can find their definitions in the default one
    data = def_data + "\n" + usr_data
    config_data = yaml.safe_load(def_data)

    return config_data


_config_data = get_config_data()


def get_config(*args: str, test_mode: bool = False) -> Any:
    largs = list(map(str.lower, args))
    if test_mode:
        largs = ["test_mode", *largs]

    value = _config_data
    for arg in largs:
        for key in value.keys():
            if arg == key.lower():
                value = value[key]
                break
        else:
            break
    else:
        return value

    configpath = "/".join(map(str, largs))
    raise ConfigError(f"Config {configpath!r} not found!")


class _Config:
    def __init__(self):
        """This method is defined to remind you that this is not a static class"""
        self._test_mode = False
        logging.config.dictConfig(self.logging)
        self.register_fonts()
        self._styles = self.register_styles()

    def register_fonts(self):
        # Os arquivos .ttf podem estar em dois lugares:
        # .local/share/fonts
        # Files\Constants\fonts
        try:
            pdfmetrics.registerFont(ttfonts.TTFont("Montserrat", "Montserrat-Black.ttf"))
        except:
            endereco = os.path.join("../Controller", "DataFiles", "Files", "Constants", "Fonts", "Montserrat-Black.ttf")
            print('Endereco buscado:',endereco)
            pdfmetrics.registerFont(ttfonts.TTFont("Montserrat", endereco))
        
        try:
            pdfmetrics.registerFont(ttfonts.TTFont("Montserrat-Bold", "Montserrat-Bold.ttf"))
        except:
            endereco = os.path.join("../Controller", "DataFiles", "Files", "Constants", "Fonts", "Montserrat-Black.ttf")
            pdfmetrics.registerFont(ttfonts.TTFont("Montserrat-Bold", endereco))

        try:
            pdfmetrics.registerFont(ttfonts.TTFont("Montserrat-SemiBold", "Montserrat-SemiBold.ttf"))
        except:
            endereco = os.path.join("../Controller", "DataFiles", "Files", "Constants", "Fonts", "Montserrat-Black.ttf")
            pdfmetrics.registerFont(ttfonts.TTFont("Montserrat-SemiBold", endereco))
        
        try:
            pdfmetrics.registerFont(ttfonts.TTFont("Montserrat-Extra-Bold", "Montserrat-ExtraBold.ttf"))
        except:
            endereco = os.path.join("../Controller", "DataFiles", "Files", "Constants", "Fonts", "Montserrat-Black.ttf")
            pdfmetrics.registerFont(ttfonts.TTFont("Montserrat Extra Bold", endereco))
        
        try:
            pdfmetrics.registerFont(ttfonts.TTFont("Montserrat-Light", "Montserrat-Light.ttf"))
        except:
            endereco = os.path.join("../Controller", "DataFiles", "Files", "Constants", "Fonts", "Montserrat-Light.ttf")
            pdfmetrics.registerFont(ttfonts.TTFont("Montserrat Light", endereco))
        
        try:
            pdfmetrics.registerFont(ttfonts.TTFont("Montserrat-Regular", "Montserrat-Regular.ttf"))
        except:
            endereco = os.path.join("../Controller", "DataFiles", "Files", "Constants", "Fonts", "Montserrat-Black.ttf")
            pdfmetrics.registerFont(ttfonts.TTFont("Montserrat-Regular", endereco))
        
        try:
            pdfmetrics.registerFont(ttfonts.TTFont("Montserrat-Medium", "Montserrat-Medium.ttf"))
        except:
            endereco = os.path.join("../Controller", "DataFiles", "Files", "Constants", "Fonts", "Montserrat-Black.ttf")
            pdfmetrics.registerFont(ttfonts.TTFont("Montserrat-Medium", endereco))

    def register_styles(self):
        def _is_font_dict(data) -> TypeGuard[dict]:
            return isinstance(data, dict) and "fontName" in data
        
        def _register_styles(data: Union[dict, list], root: str = ""):
            # logger.debug(f"Received {root=!r} and {data=!r}")
            if _is_font_dict(data):
                style = ParagraphStyle(root.lower(), **data)
                styles.add(style)
            elif isinstance(data, dict):
                for key, value in data.items():
                    name = f"{root}.{key}" if root else key
                    _register_styles(value, name.lower())
            elif isinstance(data, list):
                for idx, value in enumerate(data):
                    name = f"{root}{idx}" if root else str(idx)
                    _register_styles(value, name.lower())
            elif isinstance(data, ParagraphStyle):
                styles.add(data)
            elif data is None:
                return
            else:
                raise TypeError(
                    f"Argument data of type {type(data).__name__} is not a valid paragraph style config: {data}")

        styles = getSampleStyleSheet()
        _register_styles(get_config(
            "styles", "ancestralidade"), "ancestralidade")

        return styles

    def get_text_mode_enabled(self):
        return self._test_mode

    def set_test_mode_enabled(self, value: bool):
        self._test_mode = value

    @property
    def styles(self):
        return self._styles

    @property
    def parceria(self):
        return "padrão"

    @cached_property
    def templates(self) -> Path:
        return get_config("paths", "templates")

    @cached_property
    def template(self) -> Path:
        return get_config("paths", "template")

    @cached_property
    def scores(self) -> Path:
        # TODO: refatorar para pegar pelas populações, provv algo do tipo:
        #    scores(self, z-score, superpop)
        return get_config("paths", "scores")

    @cached_property
    def media(self) -> Path:
        return get_config("paths", "media", "media")

    @cached_property
    def efeitos(self) -> Path:
        return get_config("paths", "media", "efeitos")

    @cached_property
    def confiabilidade(self) -> Path:
        return get_config("paths", "media", "confiabilidade")

    @cached_property
    def icones(self) -> Path:
        return get_config("paths", "media", "icones")

    @cached_property
    def pins(self) -> Path:
        return get_config("paths", "media", "pins")

    @cached_property
    def iadmix_dir(self) -> Path:
        return get_config("paths", "ancestralidade", "iadmix")

    @cached_property
    def iadmix_freq_file(self) -> Path:
        return get_config("paths", "ancestralidade", "iadmix freq file")

    @property
    def logging(self) -> dict:
        return get_config("logging")

    @cached_property
    def tqdm_loggers(self):
        return [logging.getLogger(name) for name in get_config("tqdm")]

    def page_color(self, caderno: str) -> Color:
        return get_config("styles", "cadernos", caderno, "page color")

    @property
    def geno_freq_chart_startangle(self) -> float:
        return get_config("styles", "charts", "pizza frequência populacional", "startAngle")

    @property
    def abs_chart_innerRadiusFraction(self) -> float:
        return get_config("styles", "charts", "gráfico do risco absoluto", "innerRadiusFraction")

    @property
    def has_resumos(self) -> bool:
        return get_config("laudos", self.parceria, "cadernos", "possui resumo?")

    @property
    def has_ancestralidade(self) -> bool:
        # TODO: implement switch based CLI
        return True


CONFIG = _Config()
get_config()