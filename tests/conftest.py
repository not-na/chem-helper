#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  conftest.py
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

# Fixtures

def basic_alkane(n,do_hydrogen=True):
    struct = chemhelper.notations.structural.StructuralNotation()
    
    carbons = []
    for i in range(1,n+1):
        carbons.append(struct.addCarbon(name="C%s"%i))
    
    for c in carbons:
        if carbons.index(c)==0:
            continue
        c.bindToAtom(carbons[carbons.index(c)-1])
    
    if do_hydrogen:
        struct.fillWithHydrogen()
    
    assert struct.checkValid()==[]
    
    return struct
