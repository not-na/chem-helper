#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  test_iupac2smiles.py
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

# Some constants for better readability
NONYL = "(CCCCCCCCC)"
DECYL = "(CCCCCCCCCC)"

test_cases_i2s_branched_alkane = [
    # Some Basic Tests
    ["3,4-dimethyloctane","CCC(C)C(C)CCCC"], # Tests proper flipping
    ["4,4-dimethyloctane","CCCC(C)(C)CCCC"], # Tests proper flipping and two chains on the same base
    ["4-ethylheptane","CCCC(CC)CCC"], # Tests Ethyl Groups
    ["3-ethyl-3-methylhexane","CCC(C)(CC)CCC"], # Tests order of sub-chains
    # NOTE: 3-ethyl-3-methyl-4,4-dipropylnonane was used originally, but it was ambigous and could be converted to multiple SMILES
    ["4-ethyl-4-methyl-6,6-dipropyldodecane","CCCC(C)(CC)CC(CCC)(CCC)CCCCCC"], # More sub-chains
    
    # Large Molecule Tests
    ["10,11-dinonyltricosane","CCCCCCCCCC"+NONYL+"C"+NONYL+"CCCCCCCCCCCC"], # Tests nonyl, tricosane (23) and long SMILES
    ["12,12,13,13,14,14-hexakisdecyltriacontane",("C"*12)+(DECYL*2)+"C"+(DECYL*2)+"C"+(DECYL*2)+("C"*16)], # Tests hexakis- prefix and three double-sidechains next to eachother
    ["Hectane","C"*100], # Tests very large molecules
    ["Dictane","C"*200], # More large molecules
    
    # Basic Alkanol Tests
    ["Ethanol","C(O)C"],
    ["Propanol","C(O)CC"],
    ["Hexanol","C(O)CCCCC"],
    
    # Alkanol Tests with specified position
    ["Propan-2-ol","CC(O)C"],
    
    # Currently only works in one direction
    #["Hexane-1,2-diol","C(O)C(O)CCCC"]
    ]

test_cases_s2i_special = [
    # Tests that use special features only supported by loading, not saving
    ["[C][C][C][C][C]([C])[C]([C])[C][C]","3,4-dimethyloctane"], # Atoms in Brackets TODO: dont add hydrogen if atom in brackets
    ["[CH3][CH2][CH2][CH2][CH1]([CH3])[CH1]([CH3])[CH2][CH3]","3,4-dimethyloctane"], # Atoms in brackets with explicit Hydrogen
    ["CC[CH](C)(C)CC","3,3-dimethylpentane"], # Mixed brackets/no brackets and single hydrogen without number
    ["CC(CC)C(CCCC)C","3,4-dimethyloctane"], # Start from branch
    ["CCCCC(C) C(C) CC\n","3,4-dimethyloctane"], # Space/Newline ignored
    ]

@pytest.mark.parametrize(("name","smiles"),test_cases_i2s_branched_alkane)
def test_i2s_branched_alkane(name,smiles):
    iupacname = chemhelper.notations.iupac.IUPACNotation(name)
    struct = iupacname.asStructuralFormula()
    
    print(smiles,struct.dumpAsSMILES())
    print(iupacname.loadsFromSMILES(smiles))
    
    # Check conversion from struct to smiles
    assert struct.dumpAsSMILES()==smiles
    # Check conversion from smiles to name
    assert iupacname.loadsFromSMILES(smiles) == iupacname

@pytest.mark.parametrize(("smiles","name"),test_cases_s2i_special)
def test_s2i_special(smiles,name):
    assert name==chemhelper.notations.iupac.IUPACNotation.loadsFromSMILES(smiles).name

def test_s2i_doubleletterelements():
    # Tests double-letter elements
    
    # Probably not a valid structure, but good for a testcase
    struct = chemhelper.notations.structural.StructuralNotation.loadsFromSMILES("CC(Br)C(Cl)C")
    assert struct.countAtoms()=={"C":4,"H":8,"Br":1,"Cl":1}

def main(args):
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
