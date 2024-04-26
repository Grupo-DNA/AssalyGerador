# -*- coding: utf-8 -*-
# packages/ancestralidade.py

"""Script que permite a criação de uma página PDF contendo dados da ancestralidade."""

if __name__ == "packages.ancestralidade":
    from .base import Ancestralidade, Continente, Paciente, Population
else:
    from base import Ancestralidade, Continente, Paciente, Population

class AncestralidadeMaker:
    """
    Classe reunindo funções de criação de uma característica do laudo.

    """

    def __init__(self, paciente: Paciente):
        self.paciente = paciente
        self.ancestries = paciente.ancestralidades
    
    def pin_coordinates(self, pin):
        pin_coordinates_dict = {
            "Itália, Europa (TSI)": (288, 280),
            "África": (302, 248),
            "Europa do Norte e Ocidental (CEU)": (273, 319),
            "China (CHB)": (405, 267),
            "Japão (JPT)": (438, 280),
            "Ásia (CHD)": (378, 304)
        }
        return pin_coordinates_dict[pin]

    def display_name(self, place):
        """Nome da população para o gráfico de barras."""
        display_names = {
            "TSI": "Itália, Europa (TSI)",
            "MKK": "África",
            "CEU": "Europa do Norte e Ocidental (CEU)",
            "CHB": "China (CHB)",
            "JPT": "Japão (JPT)",
            "YRI": "África",
            "LWK": "África",
            "CHD": "Ásia (CHD)"
        }
        return display_names[place]

    def highest(self):
        cont_name, cont_pct = Ancestralidade.get_highest_continente(
            self.ancestries)
        display_conts = {
            Continente.EUROPA: "europeus",
            Continente.AFRICA: "africanos",
            Continente.AMERICA_DO_NORTE: "norte-americanos",
            Continente.AMERICA_DO_SUL: "sul-americanos",
            Continente.ASIA: "asiáticos",
            Continente.OCEANIA: "oceânicos"
        }
        return display_conts[cont_name], cont_pct