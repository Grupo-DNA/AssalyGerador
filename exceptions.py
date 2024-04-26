# -*- coding: utf-8 -*-
# packages/exceptions.py

"""Módulo que contém os erros para uso interno do programa"""


class Error(Exception):
    """Base class for this program exceptions."""


class WrongNumberOfArguments(Error, TypeError):
    """Trying to pass wrong number of arguments to a function."""


class WrongTypeOfArguments(Error, TypeError):
    """Trying to pass arguments of wrong type to a function."""


class LoginError(Error, ValueError):
    """Trying to login with invalid credentials."""


class RestrictionError(Error, ValueError):
    """Trying to create a Caracteristica that's unavailable for a Paciente."""


class NoSNPFoundError(Error, ValueError):
    """Trying to create a Caracteristica from a genofile that doesn't have any SNP for it."""


class ConfigError(Error):
    """Trying to get a mal-configured value."""

class InvalidBookError(Exception):
    def __init__(self, message):
        super().__init__(message=None)

class InvalidIdError(Exception):
    def __init__(self, message=None):
        super().__init__(message)
