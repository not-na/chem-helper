#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  structural_basic.py
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

# Test cases

@pytest.mark.parametrize(("n", "expected_c", "expected_h"), [
    (   1,   1,   4),
    (   2,   2,   6),
    (   3,   3,   8),
    (   4,   4,  10),
    (  10,  10,  22),
    ( 100, 100, 202),
    (1000,1000,2002),
])
def test_count(n,expected_c,expected_h):
    atoms = basic_alkane(n).countAtoms()
    
    assert atoms["C"]==expected_c
    assert atoms["H"]==expected_h
    assert sorted(atoms.keys())==["C","H"] # sort needed to prevent random order of keys

def test_empty_count():
    struct = chemhelper.notations.structural.StructuralNotation()
    
    assert struct.checkValid()==[]
    
    struct.fillWithHydrogen()
    
    assert struct.checkValid()==[]
    
    assert struct.countAtoms()=={}

# Interactive mode

def main(args):
    struct = chemhelper.notations.structural.StructuralNotation()
    
    c1 = struct.addCarbon()
    c2 = struct.addCarbon()
    c3 = struct.addCarbon()
    
    c1.bindToAtom(c2)
    c2.bindToAtom(c3)
    
    struct.fillWithHydrogen()
    print(struct.checkValid())
    
    print(struct.countAtoms())
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
