#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  iupac_numerical.py
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

import traceback

import pytest

import chemhelper

SUFFIX = ""

test_cases_general = [
    [1,"mono"],
    [2,"di"],
    [3,"tri"],
    [4,"tetra"],
    [5,"penta"],
    [6,"hexa"],
    [7,"hepta"],
    [8,"octa"],
    [9,"nona"],
    [10,"deca"],
    [11,"undeca"],
    [12,"dodeca"],
    [13,"trideca"],
    [14,"tetradeca"],
    [15,"pentadeca"],
    [16,"hexadeca"],
    [17,"heptadeca"],
    [18,"octadeca"],
    [19,"nonadeca"],
    [20,"icosa"],
    [21,"henicosa"],
    [22,"docosa"],
    [23,"tricosa"],
    [24,"tetracosa"],
    [25,"pentacosa"],
    [26,"hexacosa"],
    [27,"heptacosa"],
    [28,"octacosa"],
    [29,"nonacosa"],
    [30,"triaconta"],
    [31,"hentriaconta"],
    [32,"dotriaconta"],
    [33,"tritriaconta"],
    [40,"tetraconta"],
    [50,"pentaconta"],
    [60,"hexaconta"],
    [70,"heptaconta"],
    [80,"octaconta"],
    [90,"nonaconta"],
    [100,"hecta"],
    [132,"dotriacontahecta"],
    # from http://www.acdlabs.com/iupac/nomenclature/93/r93_317.htm
    # note that the -kis has been removed
    [231,"hentriacontadicta"],
    # from http://www.chem.qmul.ac.uk/iupac/misc/numb.html#1
    [468,"octahexacontatetracta"]
    ]

test_cases_alkane = [
    [1,"meth"],
    [2,"eth"],
    [3,"prop"],
    [4,"but"],
    ]

# Simple tests
@pytest.mark.parametrize(("n", "expected"), test_cases_alkane)
def test_basic_alkane_n2name(n,expected):
    assert expected==chemhelper.notations.iupac.getAlkanePrefix(n)

@pytest.mark.parametrize(("n", "expected"), test_cases_general)
def test_basic_general_n2name(n,expected):
    assert expected==chemhelper.notations.iupac.getNumericalPrefix(n)


@pytest.mark.parametrize(("n", "name"), test_cases_alkane)
def test_reverse_alkane_name2n(n,name):
    assert n==chemhelper.notations.iupac.parseAlkanePrefix(name)

@pytest.mark.parametrize(("n", "name"), test_cases_general)
def test_reverse_general_name2n(n,name):
    assert n==chemhelper.notations.iupac.parseNumericalPrefix(name)

# Tests with string prepended

@pytest.mark.parametrize(("n", "name"), test_cases_alkane)
def test_reverse_alkane_name2n_prefixed(n,name):
    assert n==chemhelper.notations.iupac.parseAlkanePrefix("abcdef"+name,True)[0]

@pytest.mark.parametrize(("n", "name"), test_cases_general)
def test_reverse_general_name2n_prefixed(n,name):
    assert n==chemhelper.notations.iupac.parseNumericalPrefix("abcdef"+name,True)[0]

# TODO: clean this up and make it an external tool
def main(args):
    global SUFFIX
    
    if len(args)==2:
        SUFFIX = args[1]
        print("Found commandline argument for suffix")
    if SUFFIX!="":
        print("Automatically adding suffix '-%s'"%SUFFIX)
    
    # Test cases
    print("Running automated test cases...")
    print("General")
    for n,expected in test_cases_general:
        result = chemhelper.notations.iupac.getNumericalPrefix(n)
        if result!=expected:
            print("Unexpected result:")
        print("n=%s, expected %s, got %s"%(n,expected,result))
    print("Alkane special")
    for n,expected in test_cases_alkane:
        result = chemhelper.notations.iupac.getAlkanePrefix(n)
        if result!=expected:
            print("Unexpected result:")
        print("n=%s, expected %s%s, got %s%s"%(n,expected,SUFFIX,result,SUFFIX))
    
    # Interactive mode
    
    print("Enter number to convert, press enter without number to quit")
    
    while True:
        i = input("n:")
        try:
            if i=="":
                break
            elif len(i)>4:
                print("Only numbers up to 9999 are supported")
                continue
            elif int(i)<=0:
                print("Only numbers higher than zero are allowed")
                continue
            
            print("n=%s, got %s"%(int(i),chemhelper.notations.iupac.getAlkanePrefix(int(i)).title()+SUFFIX))
        except Exception:
            print("Got an error while trying to convert number:")
            traceback.print_exc()
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
