#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  test_structural2iupac.py
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

import pytest

import chemhelper

from conftest import basic_alkane

ELEMENT_STR = "{element}[{count}]"

test_cases_i2s_branched_alkane = [
    # Some basic tests
    ["4-Ethylheptane",              7,  9,  20],
    ["4-Ethyloctane",               8, 10,  22],
    # Would use 2,3-dimethyloctane, but the current algo sometimes chooses C 2-methyl #1 instead of C1
    # Technically both are valid, but this behavior is too complicated to simply test, needs non-parametrized tests
    ["3,4-Dimethyloctane",          8, 10,  22],
    
    # Various other tests
    # Checks sort order of both positions and alkyl groups
    ["4,5-Diethyl-3,4-dimethyloctane", 8, 14, 30],
    ["3,4,5,6,7,8,9-Heptamethylundecane", 11, 18, 38],
    ["3,3,4,4,5,5-Hexamethylheptane", 7, 13, 28],
    
    # Tests for ambiguous replacements, see http://www.acdlabs.com/iupac/nomenclature/93/r93_55.htm#r_0_1_4_2
    # Checks if trisdecyl is parsed correctly
    ["12,13,14-Trisdecyltriacontane", 30, 60, 122],
    # Checks if pentakisdecyl is parsed correctly
    ["12,13,14,15,16-Pentakisdecylhectane", 100, 150, 302],
    
    ]

test_cases_autocorrect = [
    # Tests autocorrection
    
    # Basic case correction
    ["ethane","Ethane",                         "C[2]H[6]"],
    ["METHANOL","Methanol",                     "CH[4]O"],
    ["3,4-dImEThylheXANoL","3,4-Dimethylhexanol","C[8]H[18]O"],
    
    # Autofill of Group positions
    ["Aminomethane","1-Aminomethane",           "CH[5]N"],
    ["Tetraaminomethane","1,1,1,1-Tetraaminomethane","CH[8]N[4]"],
    ["Pentaiodoethane","1,1,1,2,2-Pentaiodoethane","C[2]HI[5]"],
    ["Tetrafluorobutane","1,2,3,4-Tetrafluorobutane","C[4]H[6]F[4]"],
    
    # Autofill of mixed Group positions
    ["2-Amino-fluorobutane","2-Amino-1-fluorobutane","C[4]H[10]FN"],
    
    # Ordering of Groups
    ["6-Butyl-20-aminohectane","20-Amino-6-butylhectane","C[104]H[211]N"],
    
    # Correction of main chain with implausible side chains
    ["1-Methylmethane","Ethane",                  "C[2]H[6]"],
    ["1,1,1,1-Tetramethylmethane","2,2-Dimethylpropane","C[5]H[12]"],
    ["2-Ethyloctane","3-Methylnonane",          "C[10]H[22]"],
    
    # Simplification of unneccessary locants
    ["Ethan-1-ol","Ethanol",                    "C[2]H[6]O"],
    ]

def test_s2i_branched_alkane_1():
    # Note that the backbone is not the straight line
    #
    # C-C-C-C-C-C
    #       |
    #       C
    #       |
    #     C-C
    # Numbering:
    # 1-2-3-4-5-6
    #       |
    #       7-8-9
    # Right chain is 1-2-3-4-7-8-9
    
    struct = chemhelper.notations.structural.StructuralNotation()
    
    carbons = []
    for i in range(1,10):
        carbons.append(struct.addCarbon(name="C%s"%i))
    
    c1,c2,c3,c4,c5,c6,c7,c8,c9 = carbons
    
    c1.bindToAtom(c2)
    c2.bindToAtom(c3)
    c3.bindToAtom(c4)
    c4.bindToAtom(c5)
    c4.bindToAtom(c7)
    c5.bindToAtom(c6)
    c7.bindToAtom(c8)
    c8.bindToAtom(c9)
    
    struct.fillWithHydrogen()
    
    assert struct.checkValid()==[]
    
    assert struct.countAtoms()=={"C":9,"H":20}
    
    iupacname = struct.asIUPACName()
    
    assert iupacname.name == "4-Ethylheptane"

def test_s2i_branched_alkane_2():
    # Other example with reverse chain
    #
    # C-C-C-C-C-C-C
    #     |
    #     C-C-C-C
    #
    # This example demonstrates that simply recursively searching from one end will not work
    # Numbering:
    # 1-2-3-4-5-6-7
    #     |
    #     8-9-10
    # Right chain is 7-6-5-4-3-8-9-10
    
    struct = chemhelper.notations.structural.StructuralNotation()
    
    carbons = []
    for i in range(1,11):
        carbons.append(struct.addCarbon(name="C%s"%i))
    
    c1,c2,c3,c4,c5,c6,c7,c8,c9,c10 = carbons
    
    c1.bindToAtom(c2)
    c2.bindToAtom(c3)
    c3.bindToAtom(c4)
    c4.bindToAtom(c5)
    c5.bindToAtom(c6)
    c6.bindToAtom(c7)
    c3.bindToAtom(c8)
    c8.bindToAtom(c9)
    c9.bindToAtom(c10)
    
    struct.fillWithHydrogen()
    
    assert struct.checkValid()==[]
    
    assert struct.countAtoms()=={"C":10,"H":22}
    
    iupacname = struct.asIUPACName()
    
    assert iupacname.name == "4-Ethyloctane"

@pytest.mark.parametrize(("name", "chainlen","exp_c","exp_h"), test_cases_i2s_branched_alkane)
def test_i2s_branched_alkane(name,chainlen,exp_c,exp_h):
    iupacname = chemhelper.notations.iupac.IUPACNotation(name)
    struct = iupacname.asStructuralFormula()
    
    assert struct.countAtoms()=={"C":exp_c,"H":exp_h}
    
    backbone = [c.name for c in struct.getCarbonBackbone()]
    
    backbone_expected = ["C%s"%(i+1) for i in range(chainlen)]
    
    print(backbone)
    print(backbone_expected)
    
    assert len(backbone)==len(backbone_expected)
    assert backbone == backbone_expected or \
           backbone == list(reversed(backbone_expected))
    
    assert struct.asIUPACName()==iupacname

@pytest.mark.parametrize(("orig","exp","sum_formula"), test_cases_autocorrect)
def test_autocorrect(orig,exp,sum_formula):
    orig_i = chemhelper.notations.iupac.IUPACNotation(orig)
    exp_i = chemhelper.notations.iupac.IUPACNotation(exp)
    
    orig_s = orig_i.asStructuralFormula()
    assert orig_s.asIUPACName()==exp_i
    
    orig_sum = orig_s.getSumFormula(ELEMENT_STR)
    assert orig_sum == sum_formula
    
    exp_s = exp_i.asStructuralFormula()
    assert exp_s.asIUPACName()==exp_i
    
    exp_sum = exp_s.getSumFormula(ELEMENT_STR)
    assert exp_sum == sum_formula
