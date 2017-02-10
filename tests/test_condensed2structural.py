#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  test_condensed2structural.py
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

test_cases_c2s_basic_alkane = [
    # Tests for simple mode
    ["CH3(CH2)2CH3",            4,      4,      10],
    ["CH3-(CH2)2-CH3",          4,      4,      10],
    ["CH3CH2CH3",               3,      3,       8],
    ["CH3-CH2-CH3",             3,      3,       8],
    ["CH3(CH2)123CH3",        125,    125,     252],
    ["CH3-(CH2)123-CH3",      125,    125,     252],
    ]

@pytest.mark.parametrize(("formula", "chainlen","exp_c","exp_h"), test_cases_c2s_basic_alkane)
def test_c2s_basic_alkane(formula,chainlen,exp_c,exp_h):
    formula = chemhelper.notations.condensed.CondensedMolecularNotation(formula)
    struct = formula.asStructuralFormula()
    
    assert struct.countAtoms()=={"C":exp_c,"H":exp_h}
    
    backbone = [c.name for c in struct.getCarbonBackbone()]
    
    backbone_expected = ["C%s"%(i+1) for i in range(chainlen)]
    
    print(backbone)
    print(backbone_expected)
    
    assert len(backbone)==len(backbone_expected)
    assert backbone == backbone_expected or \
           backbone == list(reversed(backbone_expected))


def main(args):
    cond = chemhelper.notations.condensed.CondensedMolecularNotation("CH3(CH2)8CH3")
    print(cond)
    print(cond.asStructuralFormula())
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
