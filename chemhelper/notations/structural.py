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
from ..elements import Atom, Carbon, Hydrogen, Oxygen

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
    
    def countAtoms(self):
        count = {}
        for atom in self.atoms:
            if atom.symbol not in count:
                count[atom.symbol]=0
            count[atom.symbol]+=1
        return count
    
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
    
    def asIUPACName(self):
        if self.checkValid()!=[]:
            raise errors.IncompleteFormulaError("At least %s atoms are invalid, cannot convert if not valid"%len(self.checkValid()))
        
        # Check for methane, special
        if self.countAtoms()=={"C":1,"H":4}:
            return iupac.IUPACNotation("Methane")
        
        # find longest carbon chain
        backbone = self.getCarbonBackbone()
        
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
                else:
                    raise errors.UnsupportedElementError("Element '%s' is not currently supported"%neighbour.symbol)
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
        
        # Group together the alkyl groups
        for n,grouptype,extradata in groups:
            if grouptype=="alkyl":
                if extradata["n"] not in alkyl_groups:
                    alkyl_groups[extradata["n"]]=[]
                alkyl_groups[extradata["n"]].append([n,grouptype,extradata])
            else:
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
        
        # Create the base name
        base_name = iupac.getAlkanePrefix(len(backbone))+"ane"
        
        # Combine it
        out = alkyl_prefix+base_name
        
        if not out[0].isdigit():
            out = out.title()
        
        out = iupac.IUPACNotation(out)
        return out
    
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
                else:
                    raise errors.UnsupportedElementError("Element '%s' is not currently supported"%neighbour.symbol)
        
        # Compile output
        out = ""
        n = 0
        for c in backbone:
            # Add the base carbon
            n+=1
            out += "C"
            
            # Add all groups
            for _,gtype,gdata in groups.get(n,[]):
                if gtype == "alkyl":
                    # Add parentheses containing the group
                    out += "("+("C"*gdata["n"])+")"
                else:
                    raise errors.InvalidGroupError("Unknown group type '%s'"%gtype)
        
        return out
    
    # Load from String Methods
    # TODO
    
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
                            stack.append((neighbour))
                        # else is missing, other elements handled when they are parsed in the queue
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
        if self.countAtoms()=={"C":1,"H":4}:
            for atom in self.atoms:
                if atom.symbol=="C":
                    return [atom]
        
        # List of all dead-end carbon atoms
        ends = []
        for atom in self.atoms:
            if atom.symbol=="C":
                f = 0
                for other in atom.bindings.keys():
                    if other.symbol=="C":
                        f +=1
                if f<=1:
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
