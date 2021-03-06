#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  errors.py
#  
#  Copyright 2017 notna <notna@apparat.org>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

class ChemError(Exception):pass

class IncompleteFormulaError(Exception):pass
class InvalidFormulaError(Exception):pass
class InvalidMultiplier(Exception):pass
class MultipleMoleculesError(Exception):pass
class BaseAtomOutOfRangeError(Exception):pass
class FormulaTooLargeError(Exception):pass

class NotAnAtomError(Exception):pass

class BindingError(Exception):pass
class AlreadyBoundError(BindingError):pass
class NotEnoughBindingsError(BindingError):pass
class NotBoundError(BindingError):pass

class InvalidPrefixError(Exception):pass
class InvalidSuffixError(Exception):pass

class CyclicMoleculeError(Exception):pass
class InvalidGroupError(Exception):pass

class UnsupportedFeatureError(NotImplementedError):pass
class UnsupportedElementError(NotImplementedError):pass
class UnsupportedGroupError(NotImplementedError):pass
class UnsupportedFormulaTypeError(NotImplementedError):pass
class UnsupportedBindingError(NotImplementedError):pass

class InternalError(Exception):pass

class SMILESError(Exception):pass
class SMILESSyntaxError(SMILESError):pass
class UnsupportedSMILESFeatureError(NotImplementedError):pass
