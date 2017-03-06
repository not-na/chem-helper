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
    
    # Conversion Methods
    def asStructuralFormula(self):
        raise NotImplementedError("%s cannot be converted to a Structural Formula"%self.__class__.__name__)
    
    def asMolecularFormula(self):
        raise NotImplementedError("%s cannot be converted to a Molecular Formula"%self.__class__.__name__)
    
    def asCondensedFormula(self):
        raise NotImplementedError("%s cannot be converted to IUPAC Nomenclature"%self.__class__.__name__)
    
    def asIUPACName(self):
        raise NotImplementedError("%s cannot be converted to IUPAC Nomenclature"%self.__class__.__name__)
    
    # Save to File Methods
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
    
    # Save to String Methods
    def dumpAsSMILES(self):
        raise NotImplementedError("%s cannot be exported as a SMILES File"%self.__class__.__name__)
    
    def dumpAsInChI(self):
        raise NotImplementedError("%s cannot be exported as a InChI File"%self.__class__.__name__)
    
    # Load from File Methods
    @classmethod
    def loadFromSMILES(cls,fname):
        with open(fname,"r") as f:
            data = f.read()
        return cls.loadsFromSMILES(data)
    
    @classmethod
    def loadFromInChI(cls,fname):
        with open(fname,"r") as f:
            data = f.read()
        return cls.loadsFromInChI(data)
    
    # Load from String Methods
    @classmethod
    def loadsFromSMILES(cls,data):
        raise NotImplementedError("%s cannot be loaded from SMILES"%cls.__name__)
    
    @classmethod
    def loadsFromInChI(cls,data):
        raise NotImplementedError("%s cannot be loaded from InChI"%cls.__name__)
    
    # General Dump Methods
    def dump(self,fname,mime=None):
        if mime is None:
            mime = self.mimeTypeFromFilename(fname)
        
        if mime=="chemical/x-daylight-smiles":
            return self.saveAsSMILES(fname)
        elif mime=="chemical/x-inchi":
            return self.saveAsInChI(fname)
        else:
            raise ValueError("Invalid Mimetype '%s'"%mime)
    def dumps(self,mime):
        if mime=="chemical/x-daylight-smiles":
            return self.dumpAsSMILES()
        elif mime=="chemical/x-inchi":
            return self.dumpAsInChI()
        else:
            raise ValueError("Invalid Mimetype '%s'"%mime)
    
    # General Load Methods
    @classmethod
    def load(cls,fname,mime=None):
        if mime is None:
            mime = cls.mimeTypeFromFilename(fname)
        
        if mime=="chemical/x-daylight-smiles":
            return cls.loadFromSMILES(fname)
        elif mime=="chemical/x-inchi":
            return cls.loadFromInChI(fname)
        else:
            raise ValueError("Invalid Mimetype '%'"%mime)
    @classmethod
    def loads(cls,data,mime):
        if mime=="chemical/x-daylight-smiles":
            return cls.loadsFromSMILES(data)
        elif mime=="chemical/x-inchi":
            return cls.loadsFromInChI(data)
        else:
            raise ValueError("Invalid Mimetype '%s'"%mime)
    
    @staticmethod
    def mimeTypeFromFilename(fname):
        if fname is None:
            raise ValueError("Filename may not be None")
        elif fname.endswith("smi") or fname.endswith("smiles"):
            # SMILES
            mimetype = "chemical/x-daylight-smiles"
        elif fname.endswith("inchi"):
            # InChI
            mimetype = "chemical/x-inchi"
        else:
            # Defaults to InChI
            mimetype = "chemical/x-inchi"
    

from . import structural
#from . import molecular
from . import iupac
from . import condensed
