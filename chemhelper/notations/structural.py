#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  structural.py
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

import sys
import time

from collections import defaultdict
from functools import reduce

from . import BaseNotation
from . import iupac
from .. import errors
from .. import elements
from ..elements import Atom, Carbon, Hydrogen, Oxygen, Nitrogen, Sulfur, Phosphorus, Fluorine, Chlorine, Bromine, Iodine, Boron

iupac.structural = sys.modules["chemhelper.notations.structural"] # to avoid circular dependency

class StructuralNotation(BaseNotation):
    def __init__(self):
        self.atoms = set()
    
    # Structure Modification Methods
    def addAtom(self,atom):
        self.atoms.add(atom)
    
    def addCarbon(self,pos=None,name=""):
        a = Carbon(self,pos,name)
        self.atoms.add(a)
        return a
    def addHydrogen(self,pos=None,name=""):
        a = Hydrogen(self,pos,name)
        self.atoms.add(a)
        return a
    def addOxygen(self,pos=None,name=""):
        a = Oxygen(self,pos,name)
        self.atoms.add(a)
        return a
    def addNitrogen(self,pos=None,name=""):
        a = Nitrogen(self,pos,name)
        self.atoms.add(a)
        return a
    def addSulfur(self,pos=None,name=""):
        a = Sulfur(self,pos,name)
        self.atoms.add(a)
        return a
    def addPhosphorus(self,pos=None,name=""):
        a = Phosphorus(self,pos,name)
        self.atoms.add(a)
        return a
    def addFluorine(self,pos=None,name=""):
        a = Fluorine(self,pos,name)
        self.atoms.add(a)
        return a
    def addChlorine(self,pos=None,name=""):
        a = Chlorine(self,pos,name)
        self.atoms.add(a)
        return a
    def addBromine(self,pos=None,name=""):
        a = Bromine(self,pos,name)
        self.atoms.add(a)
        return a
    def addIodine(self,pos=None,name=""):
        a = Iodine(self,pos,name)
        self.atoms.add(a)
        return a
    def addBoron(self,pos=None,name=""):
        a = Boron(self,pos,name)
        self.atoms.add(a)
        return a
    
    def countAtoms(self):
        count = {}
        for atom in self.atoms:
            if atom.symbol not in count:
                count[atom.symbol]=0
            count[atom.symbol]+=1
        return count
    
    def getSumFormula(self,element_str="{element}<sub>{count}</sub>",atomOrder=["C","H","O","F","Cl","Br","I"]):
        sum_formula = ""
        count = self.countAtoms()
        elements = set()
        
        # Generate Sum Formula
        for element in atomOrder:
            if count.get(element,0)==0:
                continue # No atoms of this type here, ignore it
            elif count[element]==1:
                sum_formula+=element
                elements.add(element)
            else:
                sum_formula+=element_str.format(
                    element=element,
                    count=count[element],
                )
                elements.add(element)
        
        # Check that all atoms have been accounted for
        for element,n in count.items():
            if n==0:
                continue
            elif element not in elements:
                raise errors.UnsupportedElementError("Element %s is currently not supported by sum formulas"%element)
        
        return sum_formula
    
    def fillWithHydrogen(self):
        for atom in set(self.atoms):
            atom.fillWithHydrogen()
    def checkValid(self):
        out = []
        for atom in self.atoms:
            totbinds = sum(atom.bindings.values())
            if totbinds>atom.max_bindings:
                # TODO: add support for ions
                out.append((atom,">"))
            elif totbinds<atom.max_bindings:
                out.append((atom,"<"))
            else:
                continue
        return out
    
    # Conversion Methods
    def asStructuralFormula(self):
        return self
    
    # Old version of algorithm
    # Most sub-routines and sub-algorithms have been ported over to the more flexible newer algorithm
    # There is no real reason to use this old algorithm
    """
    def asIUPACName_OLD(self):
        if self.checkValid()!=[]:
            raise errors.IncompleteFormulaError("At least %s atoms are invalid, cannot convert if not valid"%len(self.checkValid()))
        
        # Check for methane, special
        if self.countAtoms()=={"C":1,"H":4}:
            return iupac.IUPACNotation("Methane")
        elif self.countAtoms().get("C",0)==0:
            # Prevents bugs further down
            return iupac.IUPACNotation("")
        
        # Check that all atoms are connected to eachother, to prevent bugs with multiple molecules in one formula
        self.checkConnected(True)
        
        # find longest carbon chain
        backbone = self.getCarbonBackbone()
        if len(backbone)>9999:
            raise errors.FormulaTooLargeError("Backbone is %s atoms long, only up to 9999 supported"%len(backbone))
        
        # List of (position,type,extradata)
        groups = []
        
        # Parse branches
        n = 0
        for c in backbone:
            # Go through each atom of the backbone and count the number
            n+=1
            for neighbour in c.bindings:
                if neighbour in backbone:
                    # Neighbour is part of the backbone
                    continue
                elif neighbour.symbol=="H":
                    # Hydrogen is (currently) ignored, as it is not relevant
                    continue
                elif neighbour.symbol=="C":
                    # Found a branch
                    # Follow it and measure its length
                    grouptype,extradata = self.analyzeBranch(backbone,c,neighbour)
                    groups.append([n,grouptype,extradata])
                elif neighbour.symbol=="O":
                    # Found an oxygen side-branch
                    # Check if it is a Hydroxy Group by checking the binding
                    if c.bindings[neighbour] == 1:
                        # Hydroxy Group
                        h = 0
                        for n2 in neighbour.bindings:
                            if n2.symbol=="H":
                                h+=1
                        if h==1:
                            grouptype = "hydroxyl"
                            extradata = {"c":c,"n":n}
                            groups.append([n,grouptype,extradata])
                        else:
                            # Not bound to a hydrogen on the other end, not supported
                            raise errors.UnsupportedGroupError("Non-Hydroxy Oxygen based Groups are currently not supported")
                    elif c.bindings[neighbour] == 2:
                        # Keto Group
                        raise errors.UnsupportedGroupError("Keto Groups are currently not supported")
                    else:
                        # Not chemically possible
                        raise errors.InvalidFormulaError("Triple bindings are not possible for oxygen atoms")
                elif neighbour.symbol=="F":
                    # Found a Fluoro group
                    # No further checking required, since Fluorine only takes one binding
                    grouptype = "fluoro"
                    extradata = {"c":c,"n":n}
                    groups.append([n,grouptype,extradata])
                else:
                    raise errors.UnsupportedElementError("Element '%s' (%s) is not currently supported"%(neighbour.symbol,neighbour.atomtype))
        max_n = len(backbone)+1 # needed for an off-by-one bug
        
        unflip_sum = sum([group[0] for group in groups])
        flip_sum = sum([max_n-group[0] for group in groups])
        
        if unflip_sum<flip_sum:
            pass # Do nothing, unflipped yields the lowest numbers
        elif flip_sum<unflip_sum:
            n_groups = []
            for n,grouptype,extradata in groups:
                n_groups.append([max_n-n,grouptype,extradata])
            groups = n_groups
        
        #print("Flip: %s flip: %s unflip: %s"%(flip_sum<unflip_sum,flip_sum,unflip_sum))
        
        # Dict of n:list of groups
        alkyl_groups = {}
        
        hydroxyl_groups = {}
        
        # Group together the alkyl groups
        for n,grouptype,extradata in groups:
            if grouptype=="alkyl":
                if extradata["n"] not in alkyl_groups:
                    alkyl_groups[extradata["n"]]=[]
                alkyl_groups[extradata["n"]].append([n,grouptype,extradata])
            elif grouptype=="hydroxyl":
                if extradata["n"] not in hydroxyl_groups:
                    hydroxyl_groups[extradata["n"]]=[]
                else:
                    # there was already one group at this carbon
                    # erlenmeyer rule prevents this
                    raise errors.InvalidFormulaError("Cannot have more than one hydroxy group per carbon")
                hydroxyl_groups[extradata["n"]].append([n,grouptype,extradata])
            
            else:
                # Should only happen if functional group detection has been added for a group that is not yet supported here
                raise errors.UnsupportedGroupError("Groups of type '%s' are not yet supported"%grouptype)
        
        # Dict of prefix:data
        alkyl_prefixname = {}
        
        # Create the alkyl prefix strings
        for n,groups in alkyl_groups.items():
            d = {}
            d["n"]=n
            d["groups"]=groups
            d["positions"]=[]
            d["prefix"]=iupac.getAlkanePrefix(n)+"yl"
            for g_n,g_grouptype,g_extradata in groups:
                d["positions"].append(g_n)
            alkyl_prefixname[d["prefix"]]=d
        
        prefixnames_sorted = sorted(alkyl_prefixname.keys())
        
        prefixes = []
        for prefixname in prefixnames_sorted:
            d = alkyl_prefixname[prefixname]
            name_out = iupac.getAlkylMultiplierPrefix(len(d["groups"]),d["n"])+prefixname
            num_prefix = ",".join([str(i) for i in sorted(d["positions"])])
            name_out = num_prefix+"-"+name_out
            prefixes.append(name_out)
        alkyl_prefix = "-".join(prefixes)
        
        # Hydroxyl Groups/Alkanol
        # Already grouped together in hydroxyl_groups
        suffixes = []
        count = 0
        for n,groups in hydroxyl_groups.items():
            for g_n,g_gt,g_ed in groups:
                suffixes.append(n)
                count+=1
        suffixes = sorted(suffixes)
        if count==0:
            # No hydroxyl groups, just use -ane as the suffix
            hydroxy_suffix="ane"
        elif count==1:
            # Just one group, use -ol but no multiplier
            if suffixes[0]==1 or suffixes[0]==len(backbone):
                # at the start of the molecule, no need to specify
                # TODO: verify this
                hydroxy_suffix="anol"
            else:
                hydroxy_suffix="an-%s-ol"%suffixes[0]
        else:
            # Possibly many different suffixes
            # TODO: implement this
            raise errors.UnsupportedFeatureError("Multiple hydroxy groups cannot yet be named")
        
        # Create the base name
        base_name = iupac.getAlkanePrefix(len(backbone))+hydroxy_suffix
        
        # Combine it
        out = alkyl_prefix+base_name
        
        if not out[0].isdigit():
            # Only capitalizes first character, otherwise -ol suffixes would also be capitalized
            out=out[0].upper()+out[1:]
        
        out = iupac.IUPACNotation(out)
        return out
    """
    
    def asIUPACName(self,advanced=False):
        # Parsing is done in multiple stages
        # 1. Validate the molecule and check for edge cases
        # 2. Find carbon backbone
        # 3. Parse all branches into functional groups
        # 4. Determine if the ordering must be flipped
        # 5. Group all groups
        # 6. Seperate Prefix and Suffix Groups
        # 7. Create the prefix names for all groups of groups
        # 8. Create the suffix names for all groups of groups
        # 9. Merge prefix, base name, and suffix
        
        # Rules that are being followed:
        # 1.1: Unbranched Alkane base names
        # 2.2: Numbering of multiple side-chains, TODO: may sometimes be buggy
        # 2.3: Ordering of multiple side-chains of different nature
        # 2.5a: Multiplying prefixes in front of identical groups
        # 102.1: Halogen Derivate naming using prefixes
        # 201.1: Alcohols using -ol suffix
        
        # Rules that are partially being followed:
        # 1.2: Alkyl groups naming, only supported as functional group and not standalone
        # 2.1: Basic branches, but iso- and neo- prefixes are not supported
        
        # Rules that may be implemented in the future:
        # 2.25: Numbering of branched Alkyl groups
        # 2.4: Ordering of side chains in equivalent positions
        # 2.6: Main chain selection in case of multiple same-length candidates
        # 3.1: Double bond in main chain suffix -ene
        # 3.2: Triple bond in main chain suffix -yne
        # 3.3: Both double and triple bond in main chain
        # 3.4: Main chain selection based on max amount of double and triple bonds
        # 3.5: Endings of radicals based alkenes/alkynes
        # 3.6: Main chain selection in radicals
        # 103.1: Halogens using Radicofunctional names like "Methyl chloride"
        # 105.1: Naming of Halogens replacing all Hydrogen using per- prefix
        # 108.1: Halogen Derivate trivial names
        # 108.2: Halogen Derivate inorganic nomenclature
        # 201.2: Alcohols using hydroxy- prefix
        # 201.3: Alcohols using Radicofunctional names like "Methyl alcohol"
        # 201.4: Alcohol trivial names
        
        data = {}
        self.s2i_stage1(data)
        self.s2i_stage2(data)
        self.s2i_stage3(data)
        self.s2i_stage4(data)
        self.s2i_stage5(data)
        self.s2i_stage6(data)
        self.s2i_stage7(data)
        self.s2i_stage8(data)
        self.s2i_stage9(data)
        
        if not advanced:
            return data["out"]
        else:
            return data["out"],data
    
    def s2i_stage1(self,data):
        # Stage 1
        # 1. Validate the molecule and check for edge cases
        if self.checkValid()!=[]:
            raise errors.IncompleteFormulaError("At least %s atoms are invalid, cannot convert if not valid"%len(self.checkValid()))
        
        # Check for methane, special
        if self.countAtoms()=={"C":1,"H":4}:
            return iupac.IUPACNotation("Methane")
        elif self.countAtoms().get("C",0)==0:
            # Prevents bugs further down
            return iupac.IUPACNotation("")
        
        # Check that all atoms are connected to eachother, to prevent bugs with multiple molecules in one formula
        self.checkConnected(True)
    def s2i_stage2(self,data):
        # Stage 2
        # 2. Find carbon backbone
        
        # find longest carbon chain
        data["backbone"] = self.getCarbonBackbone()
        if len(data["backbone"])>9999:
            raise errors.FormulaTooLargeError("Backbone is %s atoms long, only up to 9999 supported"%len(data["backbone"]))
    def s2i_stage3(self,data):
        # Stage 3
        # 3. Parse all branches into functional groups
        
        # List of [base_n,type,data]
        groups = []
        
        # Parses branches
        n = 0
        for c in data["backbone"]:
            # Go through each atom of the backbone and count the number
            n+=1
            for neighbour in c.bindings:
                if neighbour in data["backbone"]:
                    # Neighbour is part of the backbone
                    continue
                elif neighbour.symbol=="H":
                    # Hydrogen is (currently) ignored, as it is not relevant
                    continue
                elif neighbour.symbol=="C":
                    # Found a carbon side-branch
                    # Follow it and measure its length
                    if c.bindings[neighbour] == 1:
                        # Alkyl Group
                        grouptype,extradata = self.analyzeBranch(data["backbone"],c,neighbour)
                        groups.append([n,grouptype,extradata])
                    else:
                        raise errors.UnsupportedGroupError("Only single bonds are supported between carbon atoms")
                elif neighbour.symbol=="O":
                    # Found an oxygen side-branch
                    # Check if it is a Hydroxy Group by checking the binding
                    if c.bindings[neighbour] == 1:
                        # Hydroxy Group
                        h = 0
                        for n2 in neighbour.bindings:
                            if n2.symbol=="H":
                                h+=1
                        if h==1:
                            grouptype = "hydroxyl"
                            extradata = {"c":c,"n":n}
                            groups.append([n,grouptype,extradata])
                        else:
                            # Not bound to a hydrogen on the other end, not supported
                            raise errors.UnsupportedGroupError("Non-Hydroxy Oxygen based Groups are currently not supported")
                    elif c.bindings[neighbour] == 2:
                        # Keto Group
                        raise errors.UnsupportedGroupError("Keto Groups are currently not supported")
                    else:
                        # Not chemically possible
                        raise errors.InvalidFormulaError("Triple bindings are not possible for oxygen atoms")
                elif neighbour.symbol=="F":
                    # Found a Fluoro group
                    # No further checking required, since Fluorine only takes one binding
                    grouptype = "fluoro"
                    extradata = {"c":c,"n":n}
                    groups.append([n,grouptype,extradata])
                elif neighbour.symbol=="Cl":
                    # Found a Chloro group
                    # No further checking required, since Chlorine only takes one binding
                    grouptype = "chloro"
                    extradata = {"c":c,"n":n}
                    groups.append([n,grouptype,extradata])
                elif neighbour.symbol=="Br":
                    # Found a Bromo group
                    # No further checking required, since Bromine only takes one binding
                    grouptype = "bromo"
                    extradata = {"c":c,"n":n}
                    groups.append([n,grouptype,extradata])
                elif neighbour.symbol=="I":
                    # Found a Iodo group
                    # No further checking required, since Iodine only takes one binding
                    grouptype = "iodo"
                    extradata = {"c":c,"n":n}
                    groups.append([n,grouptype,extradata])
                else:
                    # May happen if an unsupported element is loaded via a SMILES File
                    raise errors.UnsupportedElementError("Element '%s' (%s) is not currently supported"%(neighbour.symbol,neighbour.atomtype))
        
        data["f_groups"] = groups
    def s2i_stage4(self,data):
        # Stage 4
        # 4. Determine if the ordering must be flipped
        
        max_n = len(data["backbone"])+1 # needed for an off-by-one bug
        
        unflip_sum = sum([group[0] for group in data["f_groups"]])
        flip_sum = sum([max_n-group[0] for group in data["f_groups"]])
        
        if unflip_sum<flip_sum:
            pass # Do nothing, unflipped yields the lowest numbers
        elif flip_sum<unflip_sum:
            n_groups = []
            for n,grouptype,extradata in data["f_groups"]:
                n_groups.append([max_n-n,grouptype,extradata])
            data["f_groups"] = n_groups
    def s2i_stage5(self,data):
        # Stage 5
        # 5. Group all groups
        
        # Dict of name:list of groups
        # If name is an integer, the group is an alkyl group of the given length
        g_groups = {}
        
        for n,grouptype,edata in data["f_groups"]:
            if grouptype=="alkyl":
                # Alkyl groups are treated specially, since there are different types
                if edata["n"] not in g_groups:
                    g_groups[edata["n"]]=[]
                g_groups[edata["n"]].append([n,grouptype,edata])
            else:
                # All other groups get grouped by their grouptype
                if grouptype not in g_groups:
                    g_groups[grouptype]=[]
                g_groups[grouptype].append([n,grouptype,edata])
        
        data["g_groups"] = g_groups
    def s2i_stage6(self,data):
        # Stage 6
        # 6. Seperate Prefix and Suffix Groups
        
        data["s_groups"] = {}
        data["p_groups"] = {}
        
        for name,groups in data["g_groups"].items():
            if name=="hydroxyl":
                # Currently always as a suffix
                data["s_groups"][name]=groups
            elif name=="fluoro":
                # Currently always as a prefix
                data["p_groups"][name]=groups
            elif name=="chloro":
                # Currently always as a prefix
                data["p_groups"][name]=groups
            elif name=="bromo":
                # Currently always as a prefix
                data["p_groups"][name]=groups
            elif name=="iodo":
                # Currently always as a prefix
                data["p_groups"][name]=groups
            elif isinstance(name,int): # Alkyls use numbers to indicate the length
                # Always as a prefix
                data["p_groups"][name]=groups
            else:
                # Should only happen if a functional group has been added in Stage 3 but not yet here
                raise errors.UnsupportedGroupError("Groups of type '%s' are not yet supported"%name)
    def s2i_stage7(self,data):
        # Stage 7
        # 7. Create the prefix names for all groups of groups
        
        # list of [prefix,sortkey]
        prefixes = []
        
        # Generate the prefixes and appropriate sortkeys
        for name,groups in data["p_groups"].items():
            if isinstance(name,int): # Alkyls use numbers to indicate the length
                base = iupac.getAlkylMultiplierPrefix(len(groups),name)+iupac.getAlkanePrefix(name)+"yl"
                # Sorts the sub-groups by the position of the base carbon and then joins those positions together using commata
                num = ",".join([str(i[0]) for i in sorted(groups,key=lambda g:g[0])])
                prefix = num+"-"+base
                prefixes.append([prefix,base])
            elif name in ["hydroxyl","fluoro","chloro","bromo","iodo"]:
                # Currently, all hydroxyl groups are added as suffixes
                # This causes this to not actually be used for hydroxyl groups
                if name=="hydroxyl":
                    name=name[:-1] # strips the l
                base = iupac.getAlkylMultiplierPrefix(len(groups))+name
                # Sorts the sub-groups by the position of the base carbon and then joins those positions together using commata
                num = ",".join([str(i[0]) for i in sorted(groups,key=lambda g:g[0])])
                prefix = num+"-"+base
                prefixes.append([prefix,base])
            else:
                raise errors.UnsupportedGroupError("Groups of type %s are not supported in prefixes"%name)
        
        # Sort the prefixes by the sortkey
        # Basically equivalent to alphabetical sorting and ignoring the numbers and dashes in front of the prefixes
        prefixes = [p[0] for p in sorted(prefixes,key=lambda p:p[1])]
        
        data["prefixes"]=prefixes
        data["prefix_name"]="-".join(prefixes)
    def s2i_stage8(self,data):
        # Stage 8
        # 8. Create the suffix names for all groups of groups
        
        # TODO: make this more flexible
        # Currently only works with hydroxy groups as alkanols
        
        # Extract a list of base positions
        suffixes = [g[0] for g in data["s_groups"].get("hydroxyl",[])]
        
        for name,groups in data["s_groups"].items():
            if name!="hydroxyl":
                raise errors.UnsupportedGroupError("Groups of type %s are not supported in suffixes"%name)
        
        suffixes = sorted(suffixes)
        if len(suffixes)==0:
            # No hydroxyl groups, just use -ane as the suffix
            suffix="ane"
        elif len(suffixes)==1:
            # Just one group, use -ol but no multiplier
            if suffixes[0]==1 or suffixes[0]==len(data["backbone"]):
                # at the start or end of the molecule, no need to specify
                # TODO: verify this
                suffix="anol"
            else:
                suffix="an-%s-ol"%suffixes[0]
        else:
            # Possibly many different suffixes
            # TODO: implement this
            #raise errors.UnsupportedFeatureError("Multiple hydroxy groups cannot yet be named")
            base = iupac.getAlkylMultiplierPrefix(len(suffixes))
            num = ",".join([str(i) for i in suffixes])
            suffix = "an-"+num+"-"+base+"ol"
        
        data["suffixes"]=suffixes
        data["suffix_name"]=suffix
    def s2i_stage9(self,data):
        # Stage 9
        # 9. Merge prefix, base name, and suffix
        
        data["base_name"]=iupac.getAlkanePrefix(len(data["backbone"]))
        
        data["merged_name"] = data["prefix_name"]+data["base_name"]+data["suffix_name"]
        
        #if not data["merged_name"][0].isdigit():
        #    # Only capitalizes first character, otherwise -ol suffixes would also be capitalized
        #    data["merged_name"]=data["merged_name"][0].upper()+data["merged_name"][1:]
        # Capitalize first alphabetical character
        for char in data["merged_name"]:
            if char in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ":
                i = list(data["merged_name"]).index(char) # This works because we know it is the first character of its kind
                data["merged_name"] = data["merged_name"][:i]+char.upper()+data["merged_name"][i+1:]
                break
        
        data["out"] = iupac.IUPACNotation(data["merged_name"])
    
    def checkConnected(self,raise_error=False):
        if len(self.atoms)==0:
            return True # Prevents crashes when trying to find a starting node
        visited = set()
        
        # Simply walks the graph of all reachable atoms
        # Stores all reachable atoms
        stack = set([list(self.atoms)[0]])
        while len(stack)>0:
            atom = stack.pop()
            visited.add(atom)
            
            for neighbour in atom.bindings:
                if neighbour not in visited:
                    stack.add(neighbour)
        
        # The number of reachable atoms should equal the number of all atoms
        # If not, there are some unreachable atoms
        if len(visited)!=len(self.atoms):
            if raise_error:
                raise errors.MultipleMoleculesError("Multiple molecules in one formula detected")
            else:
                return False
        else:
            return True
    
    # Save to String Methods
    def dumpAsSMILES(self):
        if self.checkValid()!=[]:
            raise errors.IncompleteFormulaError("At least %s atoms are invalid, cannot convert if not valid"%len(self.checkValid()))
        
        # Check for methane, special
        if self.countAtoms()=={"C":1,"H":4}:
            return "C"
        
        # find longest carbon chain
        backbone = self.getCarbonBackbone()
        
        # Dict of position: [(position,type,extradata),*]
        groups = {}
        
        # Parse branches
        n = 0
        for c in backbone:
            # Go through each atom of the backbone and count the number
            n+=1
            for neighbour in c.bindings:
                if neighbour in backbone:
                    # Neighbour is part of the backbone
                    continue
                elif neighbour.symbol=="H":
                    # Hydrogen is (currently) ignored, as it is not relevant
                    continue
                elif neighbour.symbol=="C":
                    # Found a branch
                    # Follow it and measure its length
                    grouptype,extradata = self.analyzeBranch(backbone,c,neighbour)
                    group = [n,grouptype,extradata]
                    if n not in groups:
                        groups[n]=[]
                    groups[n].append(group)
                elif neighbour.symbol=="O":
                    # Only hydroxy groups are currently supported
                    h = 0
                    for n2 in neighbour.bindings:
                        if n2.symbol=="H":
                            h+=1
                    if h==1:
                        # Is a hydroxy group
                        group = [n,"hydroxy",None]
                        if n not in groups:
                            groups[n]=[]
                        groups[n].append(group)
                    else:
                        raise errors.UnsupportedGroupError("Ketogroups are currently not supported")
                elif neighbour.symbol in ["F","Cl","Br","I"]:
                    # Halogen derivative
                    group = [n,"halogen",neighbour.symbol]
                    if n not in groups:
                        groups[n]=[]
                    groups[n].append(group)
                else:
                    raise errors.UnsupportedElementError("Element '%s' is not currently supported"%neighbour.symbol)
        
        max_n = len(backbone)+1 # needed for an off-by-one bug
        
        g_flip = {}
        for i in groups.values():
            for n,grouptype,extradata in i:
                n = max_n-n
                if n not in g_flip:
                    g_flip[n] = []
                g_flip[n].append([n,grouptype,extradata])
        
        unflip_sum = sum([sum([j[0] for j in i]) for i in groups.values()])
        flip_sum = sum([sum([j[0] for j in i]) for i in g_flip.values()])
        
        if unflip_sum<flip_sum:
            pass # do nothing, unflipped yields lowest locants
        elif unflip_sum==flip_sum:
            # special case
            # TODO: implement properly using rule 2.4
            # currently, just don't flip it
            pass
        else:
            # flip, for lower locants
            groups = g_flip
        
        # Calculates sum of flipping the molecule numbering versus not flipping it
        #unflip_sum = sum([sum([g[0] for g in glist]) for glist in groups.values()])
        #flip_sum = sum([sum([max_n-g[0] for g in glist]) for glist in groups.values()])
        
        # Uncomment if problems with flipping code arise
        #print("Flip Data:")
        #for kv in groups.items():
        #    print("%s: %s"%kv)
        #print("flip: %s unflip: %s"%(flip_sum,unflip_sum))
        
        #if unflip_sum<flip_sum:
        #    pass # Do nothing, unflipped yields the lowest numbers
        #elif flip_sum<unflip_sum:
        #    n_groups = {}
        #    for i in groups.values():
        #        for n,grouptype,extradata in i:
        #            new_n = max_n-n
        #            if new_n not in n_groups:
        #                n_groups[new_n] = []
        #            n_groups[new_n].append([new_n,grouptype,extradata])
        #    groups = n_groups
        
        #print("Final Groups:")
        #for n,g in groups.items():
        #    print("n=%s:"%n)
        #    for group in g:
        #        print("\tGroup with len=%s:"%group[2]["n"])
        #        print("\t\t%s"%group)
        #        print("\t\tAtoms:")
        #        for atom in group[2]["atoms"]:
        #            print("\t\t\tAtom %s connected to %s"%(atom.name,[a.name for a in atom.bindings if atom.symbol=="C"]))
        #        print("\t\tEnd Group")
        #print("End Groups")
        
        # Compile output
        out = ""
        n = 0
        for c in backbone:
            # Add the base carbon
            n+=1
            out += "C"
            
            if n in groups:
                # If there are side chains
                
                #print("Side chains at %s:"%n)
                
                groups_sorted = sorted(groups[n],key=(lambda group: (group[2]["n"] if group[1]=="alkyl" else 0)))
                
                #for g in groups_sorted:
                #    print("\t%s"%g)
                #print("End side chains %s"%n)
                
                # Add all groups
                for _,gtype,gdata in groups_sorted:
                    if gtype == "alkyl":
                        # Add parentheses containing the group
                        out += "("+("C"*gdata["n"])+")"
                    elif gtype == "hydroxy":
                        # Add parentheses containing the hydroxy group
                        out += "(O)"
                    elif gtype == "halogen":
                        # Add parentheses containing the halogen group
                        out += "(%s)"%gdata
                    else:
                        raise errors.InvalidGroupError("Unknown group type '%s'"%gtype)
        
        return out
    
    # Load from String Methods
    @classmethod
    def loadsFromSMILES(cls,data):
        out = cls()
        
        d_list = list(data)
        c_index = 0
        binding_type = 1
        prev = None
        stack = []
        while d_list!=[]:
            # Parses the token
            # If it is a parentheses, the branch will be started/finished
            # If not, the atom will be parsed and created
            if d_list[0] in "\n ":
                # Ignores spaces and newlines
                d_list.pop(0)
                c_index+=1
                continue
            elif d_list[0] in ">":
                raise errors.UnsupportedSMILESFeatureError("SMILES Reactions are not supported")
            elif d_list[0] in "@123456789%/\\.":
                # Features in order:
                # Chirality 1x ,Cyclic Structures 10x,Directional Bonds 2x,Disconnected Structures 1x
                raise errors.UnsupportedSMILESFeatureError("SMILES Feature at %s is not yet supported"%c_index)
            elif d_list[0] in "-=#:":
                # Check if we are at the end of a group, no bindings may be specified there
                if len(d_list)<=1:
                    raise errors.SMILESSyntaxError("Dangling binding type at the end of input")
                elif d_list[1]==")":
                    raise errors.SMILESSyntaxError("Tried to specify binding type at the end of a group")
                elif binding_type!=1:
                    # Does not catch multiple single bond specifications
                    raise errors.SMILESSyntaxError("Multiple binding specification at char %s"%(c_index+1))
                if d_list[0]=="-":
                    binding_type=1
                elif d_list[0]=="=":
                    binding_type=2
                elif d_list[0]=="#":
                    binding_type=3
                elif d_list[0]==":":
                    raise errors.UnsupportedSMILESFeatureError("Aromatic bindings are not yet supported")
                d_list.pop(0)
                c_index+=1
            elif d_list[0]=="(":
                # Start of branch
                if c_index==0:
                    raise errors.SMILESSyntaxError("Cannot start branch at the beginning of the molecule")
                
                # Pushes the current prev as a backup to the stack
                stack.append(prev)
                
                # Pop the parenthesis and skip to the next cycle
                d_list.pop(0)
                c_index+=1
                continue
            elif d_list[0]==")":
                # End of branch
                if stack==[]:
                    raise errors.SMILESSyntaxError("Found closing parenthesis outside of a branch")
                
                # Retrieves the prev of the parent chain from the stack
                prev = stack.pop()
                
                # Pop the parenthesis and skip to the next cycle
                d_list.pop(0)
                c_index+=1
                continue
            elif d_list[0]=="[":
                # Bracketized Atom
                
                c_start = c_index
                
                # Collects the entire atom definition within brackets
                # Note that this might fail if there is a bracket mismatch
                a = ""
                d_list.pop(0) # Pops off the starting bracket
                while True:
                    if len(d_list)==0:
                        raise errors.SMILESSyntaxError("Missing closing bracket for bracket starting at char %s"%c_start)
                    elif d_list[0]=="]":
                        d_list.pop(0)
                        c_index+=1
                        break
                    elif d_list[0]=="[":
                        raise errors.SMILESSyntaxError("Doubly-opened bracket at %s for start bracket %s"%(c_index+1,c_start))
                    
                    c = d_list.pop(0)
                    c_index+=1
                    a+=c
                
                # Gets the element
                element = None
                for e in elements.ALL_ELEMENTS:
                    if a.startswith(e):
                        element = e
                if element is None:
                    raise errors.SMILESSyntaxError("Unknown element with Symbol '%s' at char %s"%(a[:2],c_start+2))
                
                # Creates the atom
                if element in elements.ELEMENTS.keys():
                    # Known element
                    atom = elements.ELEMENTS[element](out,name="%s from char %s"%(element,c_start))
                else:
                    raise errors.UnsupportedElementError("Unsupported element %s")
                    # TODO: add support for arbitrary elements
                out.addAtom(atom)
                
                # Adds specified hydrogen
                if a.startswith("H"):
                    # Assumes that at maximum 9 hydrogen will be added
                    if len(a)<2:
                        amount = 1 # If there is no amount specified, the default is one
                    elif a[1] not in "0123456789":
                        raise errors.SMILESSyntaxError("H for attached Hydrogen found, but following character was not a number")
                    else:
                        amount = int(a[1])
                    for i in range(amount):
                        h = Hydrogen(out,name="H #%s of char %s"%(i+1,c_start))
                        out.addAtom(a)
                        atom.bindToAtom(h)
                
                # TODO: add support for isotopes and charge
            elif len(d_list)>=2 and d_list[0]+d_list[1] in ["Cl","Br"]:
                # Double-letter elements
                element = d_list.pop(0)+d_list.pop(0)
                c_index+=2
                
                atom = elements.ELEMENTS[element](out,name="%s from char %s"%(element,c_index))
                out.addAtom(atom)
            elif d_list[0] in "BCNOPSFI":
                # Single-letter elements
                element = d_list.pop(0)
                c_index+=1
                
                atom = elements.ELEMENTS[element](out,name="%s from char %s"%(element,c_index))
                out.addAtom(atom)
            elif d_list[0]=="]":
                # Extraneous Closing Bracket
                raise errors.SMILESSyntaxError("Extraneous closing bracket at char %s"%(c_index+1))
            else:
                raise errors.SMILESSyntaxError("Invalid token at char %s"%(c_index+1))
            
            # Connect to the previous atom
            if prev is not None:
                if binding_type!=1:
                    raise errors.UnsupportedSMILESFeatureError("Multiple bindings are currently not supported")
                prev.bindToAtom(atom,binding_type)
                binding_type=1
            prev = atom
        if stack != []:
            raise errors.SMILESSyntaxError("%s parentheses have not been closed at the end, starting with %s"%(len(stack),stack[-1]))
        
        out.fillWithHydrogen()
        return out
    
    # Helper Methods and backbone extraction Algorithm for asIUPACName and dumpAsSMILES
    def analyzeBranch(self,backbone,c,start):
        d = {}
        d["atoms"]=[]
        
        grouptype=""
        if start.symbol=="H":
            raise errors.InvalidGroupError("Invalid group starting with Hydrogen")
        elif start.symbol=="C":
            grouptype="alkyl" # E.g. methyl, ethyl etc.
        else:
            raise errors.UnsupportedElementError("Element '%s' is not currently supported as a branch-starting element"%start.symbol)
        
        if grouptype=="alkyl":
            n = 0 # number of carbon atoms found in this branch
            
            # List of atoms
            stack = [(start)]
            
            visited = set()
            
            # Parse the group
            while len(stack)>0:
                atom = stack.pop()
                visited.add(atom)
                if atom.symbol=="H":
                    continue # Ignore hydrogen
                elif atom.symbol=="C":
                    n+=1
                    d["atoms"].append(atom)
                    double_branch = False
                    lc = None
                    for neighbour in atom.bindings:
                        # TODO: detect if there is a branch on the branch
                        if neighbour == c:
                            continue # Stops us from accidentally going back to the backbone
                        elif neighbour in visited:
                            continue # Stops infinite loops
                        elif neighbour in backbone:
                            # Happens if the group leads back to the backbone
                            raise errors.CyclicMoleculeError("Molecule is not acyclic")
                        elif neighbour.symbol=="H":
                            continue
                        elif neighbour.symbol=="C":
                            if double_branch:
                                # Triggers if more than one valid carbon to go to is detected on a single side-chain carbon
                                print("ERROR: Prev:")
                                print(lc)
                                print("Cur:")
                                print(neighbour)
                                raise errors.UnsupportedFeatureError("Double branch detected, not supported")
                            double_branch = True
                            lc = neighbour
                            stack.append((neighbour))
                        else:
                            # Prevents atoms of other elements being silently replaced with hydrogen
                            raise errors.UnsupportedFeatureError("Element '%s' is not currently supported in alkyl groups"%neighbour.symbol)
                else:
                    raise errors.UnsupportedElementError("Element '%s' is not currently supported in alkyl groups"%atom.symbol)
            
            # Interpret results
            if len(d["atoms"])!=n:
                raise ValueError("Inconsistent results for alkyl group (len(atoms)!=n)")
            d["n"]=n
            # name of the group will be determined after the grouping of the groups
        else:
            raise errors.InvalidGroupError("Unknown group type '%s'"%grouptype)
        
        return grouptype,d
    
    ## Begin Carbon-Backbone extraction Algorithm
    
    def getCarbonBackbone(self):
        # Example
        # Note that the backbone is not the straight line
        #
        # C-C-C-C-C-C
        #       |
        #       C
        #       |
        #     C-C
        # 
        # Other example with reverse chain
        #
        # C-C-C-C-C-C-C
        #     |
        #     C-C-C-C
        #
        # This example demonstrates that simply recursively searching from one end will not work
        
        # Check for methane, special
        if self.countAtoms()["C"]==1:
            for atom in self.atoms:
                if atom.symbol=="C":
                    return [atom]
        elif self.countAtoms()["C"]==0:
            raise errors.UnsupportedFeatureError("Molecules without a carbon backbone are currently not supported")
        
        # List of all dead-end carbon atoms
        ends = []
        for atom in self.atoms:
            #print("ATOM %s"%atom)
            # for every carbon atom
            if atom.symbol=="C":
                #print("Carbon")
                f = 0
                # go through all bound atoms
                for other in atom.bindings.keys():
                    #print("\tNeighbour %s"%atom)
                    if other.symbol=="C":
                        #print("\tNeighbouring Carbon")
                        f +=1
                if f<=1:
                    #print("END %s"%atom)
                    # if it is connected to one or fewer carbons, it is a dead end
                    ends.append(atom)
        
        # Check that there are at least two endpoints
        # Note that this would not catch most cyclo-molecules
        if len(ends)<2:
            raise errors.CyclicMoleculeError("Molecule is not acyclic")
        
        return self._longestChain(ends,)
    
    def _longestChain(self,ends):
        max_c = self.countAtoms().get("C",0)
        
        # List of (startpoint,dag)
        dags = []
        
        # Convert the Undirected Acyclic Graph representation of the molecule into a Directed Acyclic Graph
        # Multiple graphs are output, one for each starting point
        for end in ends:
            dag = {}
            stack = [(None,end)]
            visited = set()
            while len(stack)>0:
                parent,node = stack.pop()
                if node in visited:
                    raise errors.CyclicMoleculeError("Molecule is not acyclic")
                visited.add(node)
                
                dag[node]=[]
                if parent is not None:
                    dag[parent].append(node)
                
                for neighbour in node.bindings:
                    if node.bindings[neighbour]!=1:
                        raise errors.UnsupportedBindingError("%s-Binds are currently not supported"%node.bindings[neighbour])
                    elif neighbour.symbol!="C":
                        continue
                    elif neighbour in visited:
                        continue
                    stack.append((node,neighbour))
            dags.append([end,dag])
        
        longest_chain = []
        
        # Go through each endpoint and find the longest path from it
        for end,dag in dags:
            l_chain = self._processDAG(ends,end,dag)
            if len(l_chain)>len(longest_chain):
                longest_chain = l_chain
            if len(longest_chain)>=max_c:
                # Stop if the found chain is longer than the amount of carbon atoms in the molecule, cannot get longer
                # Also serves as an additional safety mechanism in case of cyclo molecules
                break # cannot get any longer
        
        return longest_chain
    
    def _processDAG(self,ends,end,dag):
        # Convert the DAG consisting of a dict with node:list of adjacent to dict with node:set of adjacent
        dag_s = {}
        for key,value in dag.items():
            dag_s[key]=set(value)
        topodag = self.toposort2(dag_s)
        
        # Needed for longest_path() as it required a pred dict
        # Normally provided by networkx, but not used here so created manually
        # TODO: create this graph only once, not for every DAG
        pred = {}
        for node in dag:
            pred[node]={}
            for neighbour in node.bindings:
                if neighbour.symbol!="C":
                    continue
                if topodag.index(neighbour)<topodag.index(node):
                    pred[node][neighbour]=None
        
        # Actually compute the longest path
        out = self.longest_path(dag,topodag,pred)
        return out
    
    # Longest path algorithm based on http://stackoverflow.com/a/17997977/3490549
    def longest_path(self,dag,topodag,pred):
        dist = {} # stores [node, distance] pair
        for node in topodag:
            # pairs of dist,node for all incoming edges
            pairs = [(dist[v][0]+1,v) for v in pred[node]] 
            if pairs:
                dist[node] = max(pairs)
            else:
                dist[node] = (0, node)
        node,(length,_)  = max(dist.items(), key=lambda x:x[1])
        path = []
        while length > 0:
            path.append(node)
            length,node = dist[node]
        return list(reversed(path))
    
    # toposort2 based on http://rosettacode.org/wiki/Topological_sort#Python
    def toposort2(self,data):
        out = []
        
        for k, v in data.items():
            v.discard(k) # Ignore self dependencies
        extra_items_in_deps = reduce(set.union, data.values()) - set(data.keys())
        data.update({item:set() for item in extra_items_in_deps})
        while True:
            ordered = set(item for item,dep in data.items() if not dep)
            if not ordered:
                break
            out.extend(sorted(ordered))
            data = {item: (dep - ordered) for item,dep in data.items()
                    if item not in ordered}
        
        if data:
            raise errors.CyclicMoleculeError("Molecule is not acyclic")
        return out
    
    ## End Carbon-Backbone extraction algorithm
