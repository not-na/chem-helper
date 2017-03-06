#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  condensed.py
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

import re

# Matches formulas of the form CH3-(CH2)n-CH3 where the dashes, parentheses and n are optional
RE_SIMPLE_FORMULA = re.compile(r"""
                                   CH3                                          # First CH3 group marks the start
                                   -?                                           # Optional dash
                                   (\()?CH2(?(1)\)|)                            # Matches middle part, with parentheses optional
                                   ([1-9])?(?(2)[0-9]*|)                        # Optional multiplier
                                   -?                                           # Optional dash
                                   CH3                                          # Last CH3 group marks the end
                                   \Z                                           # Make sure it ends there
                                   """,re.VERBOSE)
# Concatenated: CH3-?(\()+CH2(?(1)\)|)([1-9])+(?(2)[0-9]*|)-?CH3\Z

from . import BaseNotation
from . import structural, iupac
from .. import errors

class CondensedMolecularNotation(BaseNotation):
    # Formula of the form CH3(CH2)3CH3 for Pentane
    def __init__(self,formula):
        self.formula = formula
    
    # Conversion Methods
    def asCondensedFormula(self):
        return self
    
    def asStructuralFormula(self):
        struct = structural.StructuralNotation()
        
        c_amount = 0
        if RE_SIMPLE_FORMULA.match(self.formula):
            # Simple mode, e.g. straight chain alkane
            
            # Removes extra CH3 groups and dashes from start and end
            f = self.formula.rstrip("3").rstrip("CH-").lstrip("CH3").lstrip("-")
            # Formula should now be of form (?CH2)?[0-9]*
            
            # Extracts number from the end of the molecule
            n_s = ""
            print("Looping %s"%f)
            while f[-1] in "0123456789":
                if len(f[:-1])<3:
                    break # prevents the 2 from CH2 from being parsed
                n_s = f[-1]+n_s
                f = f[:-1]
            
            print("Parsing %s"%n_s)
            # Parse number
            if n_s == "":
                # For Propane, e.g. CH2-CH3-CH2
                n = 3
            else:
                n = int(n_s)+2 # +2 due to stripped CH3 ends
            
            c_amount+=n
            
            # Add base carbons to struct
            carbons = []
            for i in range(n):
                carbons.append(struct.addCarbon(name="C%s"%(i+1)))
            
            # Bind them together
            for c in carbons:
                if carbons.index(c)==0:
                    continue
                c.bindToAtom(carbons[carbons.index(c)-1])
        elif self.formula.startswith("CH3") and self.formula.endswith("CH3"):
            # Most other formula types
            raise errors.UnsupportedFormulaTypeError("Cannot convert between non-simple condensed formulas and structural formulas")
        elif self.formula.endswith("OH"):
            raise errors.UnsupportedFormulaTypeError("Alcohols represented as condensed formulas cannot be converted yet")
        else:
            raise errors.UnsupportedFormulaTypeError("The given formula cannot be converted yet")
        
        r = struct.countAtoms()
        if c_amount!=r.get("C",0):
            raise errors.InternalError("Mismatch between expected and gotten carbon atoms: exp=%s versus tot=%s"%(c_amount,r["C"]))
        
        struct.fillWithHydrogen()
        return struct
    
    def asIUPACName(self):
        return self.asStructuralFormula().asIUPACName()
    
    # Save to String Methods
    def dumpAsSMILES(self):
        return self.asStructuralFormula().dumpAsSMILES()
    
    def dumpAsInChI(self):
        return self.asStructuralFormula().dumpAsInChI()
    
    # Load from String Methods
    @classmethod
    def loadsFromSMILES(cls,data):
        return structural.StructuralNotation.loadsFromSMILES(data).asCondensedFormula()
    
    @classmethod
    def loadsFromInChI(cls,data):
        return structural.StructuralNotation.loadsFromInChI(data).asCondensedFormula()
    
    # Magic Methods
    def __repr__(self):
        return "<CondensedMolecularNotation(formula='%s')>"%self.formula
    def __str__(self):
        return self.formula
    
    def __eq__(self,other):
        return self.__class__==other.__class__ and \
               self.formula==other.formula
