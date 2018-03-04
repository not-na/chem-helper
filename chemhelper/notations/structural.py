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
                else:
                    raise errors.UnsupportedElementError("Element '%s' is not currently supported"%neighbour.symbol)
        
        max_n = len(backbone)+1 # needed for an off-by-one bug
        
        # Calculates sum of flipping the molecule numbering versus not flipping it
        unflip_sum = sum([sum([g[0] for g in glist]) for glist in groups.values()])
        flip_sum = sum([sum([max_n-g[0] for g in glist]) for glist in groups.values()])
        
        # Uncomment if problems with flipping code arise
        #print("Flip Data:")
        #for kv in groups.items():
        #    print("%s: %s"%kv)
        #print("flip: %s unflip: %s"%(flip_sum,unflip_sum))
        
        if unflip_sum<flip_sum:
            pass # Do nothing, unflipped yields the lowest numbers
        elif flip_sum<unflip_sum:
            n_groups = {}
            for i in groups.values():
                for n,grouptype,extradata in i:
                    new_n = max_n-n
                    if new_n not in n_groups:
                        n_groups[new_n] = []
                    n_groups[new_n].append([new_n,grouptype,extradata])
            groups = n_groups
        
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
                                raise errors.UnsupportedFeatureError("Double branch detected, not supported")
                            double_branch = True
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
