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

def test_branched_alkane_1():
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
    
    assert iupacname.name == "4-ethylheptane"

def test_branched_alkane_2():
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
    
    assert iupacname.name == "4-ethyloctane"


def main(args):
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
    # Should be 4-ethyloctane
    
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
    
    print("IUPAC Name:")
    print(struct.asIUPACName().name)
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
