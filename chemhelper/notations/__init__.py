#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  __init__.py
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

class BaseNotation(object):
    def checkValid(self):
        raise NotImplementedError("%s cannot be checked for validity"%self.__class__.__name__)
    def countAtoms(self):
        raise NotImplementedError("Atoms in %s cannot be counted"%self.__class__.__name__)
    
    def asStructuralFormula(self):
        raise NotImplementedError("%s cannot be converted to a Structural Formula"%self.__class__.__name__)
    
    def asMolecularFormula(self):
        raise NotImplementedError("%s cannot be converted to a Molecular Formula"%self.__class__.__name__)
    
    def asCondensedFormula(self):
        raise NotImplementedError("%s cannot be converted to IUPAC Nomenclature"%self.__class__.__name__)
    
    def asIUPACName(self):
        raise NotImplementedError("%s cannot be converted to IUPAC Nomenclature"%self.__class__.__name__)
    
    def saveAsSMILES(self,fname):
        data = self.dumpAsSMILES()
        with open(fname,"w") as f:
            f.write(data)
        return data
    
    def saveAsInChI(self,fname):
        data = self.dumpAsInChI()
        with open(fname,"w") as f:
            f.write(data)
        return data
    
    def dumpAsSMILES(self):
        raise NotImplementedError("%s cannot be exported as a SMILES File"%self.__class__.__name__)
    
    def dumpAsInChI(self):
        raise NotImplementedError("%s cannot be exported as a InChI File"%self.__class__.__name__)

from . import structural
#from . import molecular
from . import iupac
from . import condensed
