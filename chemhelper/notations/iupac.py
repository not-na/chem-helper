#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  iupac.py
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

import bidict

from . import BaseNotation
from .. import errors

SPECIAL_ALKANE_PREFIXES = bidict.bidict({
    1   :"meth",
    2   :"eth",
    3   :"prop",
    4   :"but",
    })

SPECIAL_PREFIXES = bidict.bidict({
    # First few are special
    1   :"mono",
    2   :"di",
    # Nonane, undecane and tridecane are also special
    11  :"undeca",
    13  :"trideca",
    #20  :"eicos", # uncomment to enable alternative spelling
    # Icosane is recommended by the IUPAC, it is thus used here
    })

BASE_PREFIXES_1 = bidict.bidict({
    0:"", # causes the algorithm to ignore the digit
    1:"hen",
    2:"do",
    3:"tri",
    4:"tetra",
    5:"penta",
    6:"hexa",
    7:"hepta",
    8:"octa",
    9:"nona"
    })

BASE_PREFIXES_10 = {
    0:"",
    1:"deca",
    2:"icosa",
    3:"triaconta",
    4:"tetraconta",
    5:"pentaconta",
    6:"hexaconta",
    7:"heptaconta",
    8:"octaconta",
    9:"nonaconta",
    }

BASE_PREFIXES_100 = {
    0:"",
    1:"hecta",
    2:"dicta",
    3:"tricta",
    4:"tetracta",
    5:"pentacta",
    6:"hexacta",
    7:"heptacta",
    8:"octacta",
    9:"nonacta",
    }

BASE_PREFIXES_1000 = {
    0:"",# probably unneccessary
    1:"kilia",
    2:"dilia",
    3:"trilia",
    4:"tetralia",
    5:"pentalia",
    6:"hexalia",
    7:"heptalia",
    8:"octalia",
    9:"nonalia",
    }

class IUPACNotation(BaseNotation):
    def __init__(self,name=""):
        self.name = name
    
    def checkValid(self):
        pass
    def countAtoms(self):
        pass
    
    def asIUPACName(self):
        return self
    
    def __repr__(self):
        return "<IUPACNotation(name='%s')>"%self.name
    def __str__(self):
        return self.name

def getNumericalPrefix(n):
    # Returns the numerical prefix associated with the given number
    # Note that this function was designed for use with alkane, this may influence some prefixes, e.g. 1=meth etc.
    
    if n in SPECIAL_PREFIXES:
        return SPECIAL_PREFIXES[n]
    elif n<10:
        # Single digit
        prefix = BASE_PREFIXES_1[n]#.rstrip("a")
        return prefix
    elif n<10000:
        out = ""
        
        # last digit
        # multiplier 1, ones
        digit1 = int(str(n)[-1])
        pre1 = BASE_PREFIXES_1[digit1]
        out+=pre1
        
        # 2nd digit from the right
        # multiplier 10, tens
        digit2 = int(str(n)[-2])
        pre2 = BASE_PREFIXES_10[digit2]
        if digit2==2 and len(out)>0 and out[-1] in ["a","e","i","o","y"]: # ends with a vowel
            pre2 = pre2.lstrip("i")
        out+=pre2
        
        # 3rd digit from the right
        # multiplier 100, hundreds
        if n>=100:
            digit3 = int(str(n)[-3])
            pre3 = BASE_PREFIXES_100[digit3]
            out+=pre3
        
        # 4th digit from the right
        # multiplier 1000, thousands
        if n>=1000:
            digit4 = int(str(n)[-4])
            pre4 = BASE_PREFIXES_1000[digit4]
            out+=pre4
        
        return out#.rstrip("a")
    else:
        raise NotImplementedError("Numerical Prefixes for numbers higher than 9999 are not implemented")

def parseNumericalPrefix(s_in):
    s_in = s_in.lower()
    if s_in in SPECIAL_PREFIXES.inv:
        return SPECIAL_PREFIXES.inv[s_in]
    
    if s_in=="":
        raise errors.InvalidPrefixError("Prefix cannot be an empty string")
    
    n = 0
    
    # Thousands
    if s_in.endswith("lia"): # suffix for thousands
        for val,prefix in BASE_PREFIXES_1000.items():
            if s_in.endswith(prefix) and prefix!="":
                n+=val*1000
                s_in = s_in[:-len(prefix)] # cuts off the thousands
    
    # Hundreds
    if s_in.endswith("cta"):
        for val,prefix in BASE_PREFIXES_100.items():
            if s_in.endswith(prefix) and prefix!="":
                n+=val*100
                s_in = s_in[:-len(prefix)] # cuts off the hundreds
    
    # Tens
    if s_in.endswith("aconta") or s_in.endswith("deca") or s_in.endswith("icosa"):
        for val,prefix in BASE_PREFIXES_10.items():
            if s_in.endswith(prefix) and prefix!="":
                n+=val*10
                s_in = s_in[:-len(prefix)] # cuts off the tens
    if s_in.endswith("cosa"):
        n+=20
        s_in = s_in[:-4]
    if s_in.endswith("tr"): # fix for tricosa
        s_in+="i"
    
    # Ones
    for val,prefix in BASE_PREFIXES_1.items():
        if s_in.endswith(prefix) and prefix!="":
            n+=val
            s_in = s_in[:-len(prefix)] # cuts off the ones
    
    if s_in!="":
        raise errors.InvalidPrefixError("Prefix could not be fully parsed, %s remained"%s_in)
    
    return n

def getAlkanePrefix(n):
    if n in SPECIAL_ALKANE_PREFIXES:
        return SPECIAL_ALKANE_PREFIXES[n]
    return getNumericalPrefix(n).rstrip("a")

def parseAlkanePrefix(s_in):
    s_in = s_in.lower()
    if s_in in SPECIAL_ALKANE_PREFIXES.inv:
        return SPECIAL_ALKANE_PREFIXES.inv[s_in]
    return parseNumericalPrefix(s_in)

def getAlkylMultiplierPrefix(n):
    if n==1:
        return ""
    return getNumericalPrefix(n)

def parseAlkylMultiplierPrefix(s_in):
    return parseNumericalPrefix(s_in)
