#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  elements.py
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

class Atom(object):
    atomtype = "Atom"
    symbol = "-"
    max_bindings = 0
    isotope = None # Only specify if needed
    def __init__(self,structure,pos=None,name=""):
        self.structure = structure
        
        self.pos = pos
        
        self.name = name
        
        self.bindings = {}
        
        self.num_bindings = 0
    
    def bindToAtom(self,other,bindings=1):
        if self is other:
            # Prevents bugs if atom is bound to itself
            raise errors.BindingError("Cannot bind an atom to itself")
        elif other in self.bindings:
            # Already bound, prevents weird bindings
            # TODO: merge bindings instead
            raise errors.AlreadyBoundError("Atoms are already bound to each other")
        elif self.num_bindings+bindings>self.max_bindings:
            # Not enough bindings are available to bind to this (self) atom
            raise errors.NotEnoughBindingsError("Not enough bindings available to bind from this atom")
        elif other.num_bindings+bindings>other.max_bindings:
            # Not enough bindings are available to bind to this (other) atom
            raise errors.NotEnoughBindingsError("Not enough bindings available to bind to this atom")
        
        self.bindings[other]=bindings
        other.bindings[self]=bindings
        self.num_bindings+=bindings
        other.num_bindings+=bindings
    
    def fillWithHydrogen(self):
        # TODO: Add some smart positioning for hydrogen "childs"
        while self.num_bindings<self.max_bindings:
            h = Hydrogen(self.structure)
            self.structure.addAtom(h)
            self.bindToAtom(h)
    
    def __repr__(self):
        if self.name != "":
            return "<Atom(symbol='%s',bindings=%s,name='%s')>"%(self.symbol,self.num_bindings,self.name)
        else:
            return "<Atom(symbol='%s',bindings=%s) at %s>"%(self.symbol,self.num_bindings,hex(id(self)))
    
    def __lt__(self,other):
        if isinstance(other,Atom):
            return self.name<other.name # Allows for sorting
        raise TypeError("Cannot compare %s to %s"%(self.__class__.__name__,other.__class__.__name__))
    
    def __gt__(self,other):
        if isinstance(other,Atom):
            return self.name>other.name # Allows for sorting
        raise TypeError("Cannot compare %s to %s"%(self.__class__.__name__,other.__class__.__name__))

class Carbon(Atom):
    atomtype = "Carbon"
    symbol = "C"
    max_bindings = 4

class Hydrogen(Atom):
    atomtype = "Hydrogen"
    symbol = "H"
    max_bindings = 1

class Oxygen(Atom):
    atomtype = "Oxygen"
    symbol = "O"
    max_bindings = 2
    # TODO: implement special render with "shields" for oxygen only

class Nitrogen(Atom):
    atomtype = "Nitrogen"
    symbol = "N"
    max_bindings = 3

class Sulfur(Atom):
    atomtype = "Sulfur"
    symbol = "S"
    max_bindings = 2

class Phosphorus(Atom):
    atomtype = "Phosporus"
    symbol = "P"
    max_bindings = 3

class Fluorine(Atom):
    atomtype = "Fluorine"
    symbol = "F"
    max_bindings = 1

class Chlorine(Atom):
    atomtype = "Chlorine"
    symbol = "Cl"
    max_bindings = 1

class Bromine(Atom):
    atomtype = "Bromine"
    symbol = "Br"
    max_bindings = 1

class Iodine(Atom):
    atomtype = "Iodine"
    symbol = "I"
    max_bindings = 1

class Boron(Atom):
    atomtype = "Boron"
    symbol = "B"
    max_bindings = 3

ELEMENTS = {
    "C":Carbon,
    "H":Hydrogen,
    "O":Oxygen,
    "N":Nitrogen,
    "S":Sulfur,
    "P":Phosphorus,
    "F":Fluorine,
    "Cl":Chlorine,
    "Br":Bromine,
    "I":Iodine,
    "B":Boron,
    }

# List of all IUPAC-Accepted Elements as of the 7th of March 2017 (07.03.2017)
# All elements are in order of their Number, e.g. from top-left to bottom-right

## Changelog (dates in DD.MM.YYYY format)
# 07.03.2017 by notna: first type-down of all elements, source ptable.com and handout from school
## End Changelog, for further information view the git history on https://github.com/not-na/chem-helper
ALL_ELEMENTS = [
    # First Period
    "H","He",
    # Second Period
    "Li","Be","B","C","N","O","F","Ne",
    # Third Period
    "Na","Mg","Al","Si","P","S","Cl","Ar",
    # Fourth Period
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
    # Fifth Period
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
    # Sixth Period (including Lanthanoids)
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
    # Seventh Period (including Actinoids)
    "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
    ]
