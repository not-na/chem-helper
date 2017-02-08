#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  test_structural_backbone.py
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

test_cases_simple_alkane = [
    # n,length of chain
    (1,1),
    (2,2),
    (3,3),
    (4,4),
    (5,5),
    (6,6),
    (10,10),
    (20,20),
    (30,30),
    ]

test_cases_long_alkane = [
    # n,length of chain
    (100,100),
    (500,500),
    (1000,1000),
    ]

@pytest.mark.parametrize(("n", "expected"), test_cases_simple_alkane)
def test_simple_alkane_chainlen(n,expected):
    assert expected==len(basic_alkane(n).getCarbonBackbone())


@pytest.mark.parametrize(("n", "expected"), test_cases_long_alkane)
def test_long_alkane_chainlen(n,expected):
    assert expected==len(basic_alkane(n).getCarbonBackbone())


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
    
    backbone = struct.getCarbonBackbone()
    
    assert backbone == [c1,c2,c3,c4,c7,c8,c9] or \
           backbone == [c9,c8,c7,c4,c3,c2,c1]

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
    
    backbone = struct.getCarbonBackbone()
    
    assert backbone == [c7,c6,c5,c4,c3,c8,c9,c10] or \
           backbone == [c10,c9,c8,c3,c4,c5,c6,c7]

def main(args):
    struct = chemhelper.notations.structural.StructuralNotation()
    
    c1 = struct.addCarbon(name="C1")
    c2 = struct.addCarbon(name="C2")
    c3 = struct.addCarbon(name="C3")
    c4 = struct.addCarbon(name="C4")
    
    c1.bindToAtom(c2)
    c2.bindToAtom(c3)
    c3.bindToAtom(c4)
    
    struct.fillWithHydrogen()
    
    print("Invalid atoms:")
    print(struct.checkValid())
    
    print("Atom count:")
    print(struct.countAtoms())
    
    print("Carbon backbone:")
    r = struct.getCarbonBackbone()
    print([i.name for i in r])
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
