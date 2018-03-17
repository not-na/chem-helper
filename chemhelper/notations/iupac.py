#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  iupac.py
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

import bidict

from . import BaseNotation
from .. import errors

SPECIAL_ALKANE_PREFIXES = bidict.bidict({
    1   :"meth",
    2   :"eth",
    3   :"prop",
    4   :"but",
    })

SPECIAL_PREFIXES = bidict.bidict({
    # First few are special
    1   :"mono",
    2   :"di",
    # Nonane, undecane and tridecane are also special
    11  :"undeca",
    13  :"trideca",
    #20  :"eicos", # uncomment to enable alternative spelling
    # Icosane is recommended by the IUPAC, it is thus used here
    })

BASE_PREFIXES_1 = bidict.bidict({
    0:"", # causes the algorithm to ignore the digit
    1:"hen",
    2:"do",
    3:"tri",
    4:"tetra",
    5:"penta",
    6:"hexa",
    7:"hepta",
    8:"octa",
    9:"nona"
    })

BASE_PREFIXES_10 = {
    0:"",
    1:"deca",
    2:"icosa",
    3:"triaconta",
    4:"tetraconta",
    5:"pentaconta",
    6:"hexaconta",
    7:"heptaconta",
    8:"octaconta",
    9:"nonaconta",
    }

BASE_PREFIXES_100 = {
    0:"",
    1:"hecta",
    2:"dicta",
    3:"tricta",
    4:"tetracta",
    5:"pentacta",
    6:"hexacta",
    7:"heptacta",
    8:"octacta",
    9:"nonacta",
    }

BASE_PREFIXES_1000 = {
    0:"",# probably unneccessary
    1:"kilia",
    2:"dilia",
    3:"trilia",
    4:"tetralia",
    5:"pentalia",
    6:"hexalia",
    7:"heptalia",
    8:"octalia",
    9:"nonalia",
    }

class IUPACNotation(BaseNotation):
    def __init__(self,name=""):
        self.name = name
        
        self.fgroups = {}
        
        self.init_fgroups()
    def init_fgroups(self):
        self.fgroups["alkyl"]=self.fg_alkyl
        self.fgroups["hydroxyl"]=self.fg_hydroxyl
        self.fgroups["amino"]=self.fg_amino
        self.fgroups["fluoro"]=self.fg_fluoro
        self.fgroups["chloro"]=self.fg_chloro
        self.fgroups["bromo"]=self.fg_bromo
        self.fgroups["iodo"]=self.fg_iodo
    
    # Conversion Methods
    def asIUPACName(self):
        return self
    
    # Old version of algorithm
    # Most sub-routines and sub-algorithms have been ported over to the more flexible newer algorithm
    # There is no real reason to use this old algorithm
    """
    def asStructuralFormula_OLD(self):
        struct = structural.StructuralNotation()
        
        mtype = ""
        if self.name.endswith("ane"):
            mtype = "alkane"
        elif self.name.endswith("ene"):
            raise errors.UnsupportedFormulaTypeError("Alkenes are not yet supported")
        elif self.name.endswith("yne"):
            raise errors.UnsupportedFormulaTypeError("Alkynes are not yet supported")
        elif self.name.endswith("ol"):
            #raise errors.UnsupportedFormulaTypeError("Alcohols/Alkanols are not yet supported")
            mtype = "alkanol"
        elif self.name=="":
            # Prevents bugs further down
            return struct
        else:
            raise errors.UnsupportedFormulaTypeError("Unsupported formula type")
        
        if mtype=="alkane":
            # Pre-process the formula and split into prefix and base alkane
            formula = self.name[:-3] # Removes the -ane from the end
            n,prefix = parseAlkanePrefix(formula,True)
            
            c_amount = n
            
            # Add base carbons to struct
            carbons = []
            for i in range(n):
                carbons.append(struct.addCarbon(name="C%s"%(i+1)))
            
            # Bind them together
            for c in carbons:
                if carbons.index(c)==0:
                    continue
                c.bindToAtom(carbons[carbons.index(c)-1])
            
            # To prevent unneccessary processing if it is a simple alkane
            if prefix!="":
                s_prefix = iter(prefix.split("-"))
                
                # Go through every group of numerical prefix-name
                for _i in s_prefix:
                    num_prefix,prefixname = _i,next(s_prefix)
                    
                    prefixname = prefixname[:-2] # removes the -yl from the end
                    alkyl_n,multprefix = parseAlkanePrefix(prefixname,True)
                    multiplier = parseAlkylMultiplierPrefix(multprefix)
                    pname = getAlkanePrefix(alkyl_n)+"yl"
                    
                    positions = [int(i) for i in num_prefix.split(",")]
                    
                    if multiplier!=len(positions):
                        raise errors.InvalidMultiplier("Multiplier %s announces %s positions, but %s found in alkyl group %s"%(multprefix,multiplier,len(positions),num_prefix+"-"+prefixname+"yl"))
                    
                    # Add each alkyl group to struct
                    for position in positions:
                        if position>len(carbons):
                            raise errors.BaseAtomOutOfRangeError("Tried to place sidechain at carbon #%s, but there are only %s"%(position,len(carbons)))
                        parent = carbons[position-1]
                        # Add every atom of this alkyl group at the position
                        for i2 in range(alkyl_n):
                            name = "C %s-%s #%s"%(position,pname,i2+1)
                            c = struct.addCarbon(name=name)
                            c.bindToAtom(parent)
                            parent = c
                            c_amount+=1
            
            # TODO: add sanity check by checking amount of carbon
            #if c_amount!=...
        elif mtype=="alkanol":
            # Pre-process the formula and split into prefix and base alkane
            formula = self.name[:-2] # Removes the -ol from the end
            
            # Possible formats:
            # <prefix>anol
            # <prefix>an-<n>-ol
            # <prefix>ane<npre>ol
            # <prefix>ane<npre>-<n>,<n>-ol
            # Note that this stores the information for later processing
            if formula.endswith("an"):
                # Unspecified Hydroxy-Group position
                alc_n=1
                alc_pos=None
                formula=formula[:-2]
            elif formula.endswith("-"):
                # Specified Position, but may be 1 or more groups
                pass
            else:
                # Multiple groups, but unspecified position
                alc_n,formula=parseNumericalPrefix(formula,True)
                # the "ane" suffix remains, remove it
                formula = formula[:-3]
                alc_pos=None
            
            
            n,prefix = parseAlkanePrefix(formula,True)
            
            c_amount = n
            
            # Add base carbons to struct
            carbons = []
            for i in range(n):
                carbons.append(struct.addCarbon(name="C%s"%(i+1)))
            
            # Bind them together
            for c in carbons:
                if carbons.index(c)==0:
                    continue
                c.bindToAtom(carbons[carbons.index(c)-1])
            
            # To prevent unneccessary processing if it is a simple alkane
            if prefix!="":
                s_prefix = iter(prefix.split("-"))
                
                # Go through every group of numerical prefix-name
                for _i in s_prefix:
                    try:
                        num_prefix,prefixname = _i,next(s_prefix)
                    except StopIteration:
                        break # break the loop
                    
                    prefixname = prefixname[:-2] # removes the -yl from the end
                    alkyl_n,multprefix = parseAlkanePrefix(prefixname,True)
                    multiplier = parseAlkylMultiplierPrefix(multprefix)
                    pname = getAlkanePrefix(alkyl_n)+"yl"
                    
                    positions = [int(i) for i in num_prefix.split(",")]
                    
                    if multiplier!=len(positions):
                        raise errors.InvalidMultiplier("Multiplier %s announces %s positions, but %s found in alkyl group %s"%(multprefix,multiplier,len(positions),num_prefix+"-"+prefixname+"yl"))
                    
                    # Add each alkyl group to struct
                    for position in positions:
                        if position>len(carbons):
                            raise errors.BaseAtomOutOfRangeError("Tried to place sidechain at carbon #%s, but there are only %s"%(position,len(carbons)))
                        parent = carbons[position-1]
                        # Add every atom of this alkyl group at the position
                        for i2 in range(alkyl_n):
                            name = "C %s-%s #%s"%(position,pname,i2+1)
                            c = struct.addCarbon(name=name)
                            c.bindToAtom(parent)
                            parent = c
                            c_amount+=1
            
            # Now, add the hydroxy groups
            # alc_pos is a list of numbers/positions, 1 based
            if alc_pos is None:
                if alc_n==1:
                    # at one end of the molecule
                    alc_pos = [1]
                elif alc_n==2:
                    # at opposing ends of the molecule
                    alc_pos = [1,n]
                else:
                    # alternating closing in from both ends
                    alc_pos = []
                    for i in range(alc_n):
                        pass
                        # TODO: generate alc_pos
            for p in alc_pos:
                name = "O %s"%p
                o = struct.addOxygen(name=name)
                carbons[p-1].bindToAtom(o)
            
            # TODO: add sanity check by checking amount of carbon
            #if c_amount!=...
        else:
            raise errors.UnsupportedFormulaTypeError("Unsupported formula type '%s'"%mtype)
        
        # Fill everything up with hydrogen
        struct.fillWithHydrogen()
        return struct
    """
    def asStructuralFormula(self):
        # Parsing is done in multiple stages
        # 1. Split the main name in prefixes, main chain and suffix
        # 2. Parse Prefixes and Suffixes into functional groups
        # 3. Create the functional groups, but leave them disconnected from the main chain
        # 4. Create the main chain
        # 5. Connect the functional groups to the main chain
        # 6. Fill the molecule with hydrogen
        
        # Rules that are being followed:
        # 1.1: Unbranched Alkane base names
        # 2.2: Numbering of multiple side-chains, TODO: may sometimes be buggy
        # 2.3: Ordering of multiple side-chains of different nature
        # 2.5a: Multiplying prefixes in front of identical groups
        # 102.1: Halogen Derivate naming using prefixes
        # 201.1: Alcohols using -ol suffix
        # 811.3: Amines using amino- if not the principal group
        
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
        # 811.1: Amines as part of rings
        # 811.2: Amines as generic names using nitrogen
        # 811.4: Amine Radicals and trivial names
        # 812.1: Monoamines using -amine and trivial names
        
        data = {}
        self.i2s_stage1(data)
        self.i2s_stage2(data)
        self.i2s_stage3(data)
        self.i2s_stage4(data)
        self.i2s_stage5(data)
        self.i2s_stage6(data)
        
        return data["struct"]
    
    def i2s_stage1(self,data):
        # Stage 1
        # Split the main name in prefixes, main chain and suffix
        # Name format:
        # 2,2-dimethyl-3,3,4,5-tetraethyl-4,5-dihydroxyheptan-6-ol
        # Output data:
        # prefixes: list of 2-tuples: [[(2,2),"dimethyl"],[(3,3,4,5),"tetraethyl"],[(4,5),"dihydroxy"]]
        # main_chain_length: int of main chain length: 7
        # suffixes: list of strings: [[(6,),"ol"]]
        
        # First, split by "an" from the right
        # Will potentially leave an extra e at the start
        if "an" not in self.name:
            # Check to make sure we are actually dealing with an alkane
            raise errors.UnsupportedFormulaTypeError("Cannot convert non alkane-based names")
        pre_main,suffix = self.name.rsplit("an",maxsplit=1)
        # Remove the extra e and dash at the start of the suffix, if any
        #print(suffix)
        suffix = suffix.lstrip("e")
        suffix = suffix.lstrip("-")
        print("Suffix: %s"%suffix)
        
        # Parse main chain length out of the remaining string
        # Note that parseAlkanePrefix will automatically lower-case all outputs
        data["main_chain_length"],prefix = parseAlkanePrefix(pre_main,True)
        
        # Now, we parse the prefix
        # Create the output prefix list
        data["prefixes"] = []
        
        print("Prefix: %s"%prefix)
        s_prefix = iter(prefix.split("-"))
        for _i in s_prefix:
            try:
                num_prefix,prefixname = _i,next(s_prefix)
            except StopIteration:
                # May occur if there is an odd number of dashes
                #raise errors.InvalidPrefixError("Invalid prefix with an odd number of dashes")
                break
            positions = [int(i) for i in num_prefix.split(",")]
            data["prefixes"].append([positions,prefixname])
        
        # Lastly, the suffix is parsed
        # The suffix can be similar to the prefix, but positions are not required
        
        # Create the output suffix list
        data["suffixes"] = []
        
        if suffix=="":
            return # skip
        s_suffix = iter(suffix.split("-"))
        for _i in s_suffix:
            if _i.strip("1234567890,")=="":
                # Suffix with explicit positions
                try:
                    num_suffix,suffixname = _i,next(s_suffix)
                except StopIteration:
                    # May occur if there is an odd number of dashes
                    #raise errors.InvalidSuffixError("Invalid suffix with an odd number of dashes")
                    break
                
                # TODO: support suffixes without positions
                print("Suffix: %s|%s"%(num_suffix,suffixname))
                print()
                print(num_suffix.split(","))
                positions = [int(i) for i in num_suffix.split(",")]
                data["suffixes"].append([positions,suffixname])
            else:
                # Suffix without explicit positions, only multipliers
                positions = None
                suffixname = _i
                data["suffixes"].append([positions,suffixname])
    def i2s_stage2(self,data):
        # Stage 2
        # Parse Prefixes and Suffixes into functional groups
        
        # Initialize required data structures
        data["fgroups"] = []
        # Each fg is a dict:
        # type: typename: "alkyl"
        # base: int of base carbon: 2
        # Rest of dict is dependent on type
        
        # Convert the prefix to functional groups
        for positions,prefixname in data["prefixes"]:
            if prefixname.endswith("hydroxy"):
                raise errors.UnsupportedFeatureError("Hydroxyl groups are currently not supported in prefixes")
            elif prefixname.endswith("formyl"):
                raise errors.UnsupportedFeatureError("Aldehydes are currently not supported")
            elif prefixname.endswith("oxo"):
                raise errors.UnsupportedFeatureError("Ketones are currently not supported")
            elif prefixname.endswith("carboxy"):
                raise errors.UnsupportedFeatureError("Carboxylic Acids are currently not supported")
            elif prefixname.endswith("amino"):
                # Amino group
                # Based on rule C-811.3
                multprefix = prefixname[:-5] # removes the amino from the end
                multiplier = parseAlkylMultiplierPrefix(multprefix)
                
                # Checks if positions were given, if not then generate them
                if positions is None:
                    if multiplier>data["main_chain_length"]:
                        raise errors.InvalidMultiplier("Multiplier %s of amino Group announces %s positions, but only %s available"%(
                                multprefix,multiplier,data["main_chain_length"],
                                ))
                    raise errors.UnsupportedFeatureError("Auto-placement of Amino-Groups is currently not supported")
                
                # Validates the multiplier
                if multiplier!=len(positions):
                    raise errors.InvalidMultiplier("Multiplier %s announces %s positions, but %s found in amino group %s"%(
                            multprefix,multiplier,len(positions),",".join([str(i) for i in positions])+"-"+prefixname,
                            ))
                
                # Go through every functional group specified by this suffix
                for position in positions:
                    # Store the result
                    fg = {
                        "type":"amino",
                        "base":position,
                    }
                    data["fgroups"].append(fg)
            elif prefixname.endswith("fluoro"):
                # Fluorine Group
                multprefix = prefixname[:-6] # removes the fluoro from the end
                multiplier = parseAlkylMultiplierPrefix(multprefix)
                
                if positions is None:
                    if multiplier>data["main_chain_length"]:
                        raise errors.InvalidMultiplier("Multiplier %s of fluoro Group announces %s positions, but only %s available"%(
                                multprefix,multiplier,data["main_chain_length"]
                                ))
                    positions = []
                    # positions were not specified, need to be generated
                    # Alternating from both ends and at opposite ends
                    # Note that this may not always be chemically correct
                    # TODO: make the placement smarter
                    # When the first/last atom are already specified "away" and are full, this algorithm may cause a NotEnoughBindingsError
                    for i in range(multiplier):
                        if i%2==0: # left
                            positions.append(int(i/2)+1)
                        else: # right
                            positions.append(data["main_chain_length"]-(int(i/2)+1))
                
                if multiplier!=len(positions):
                    raise errors.InvalidMultiplier("Multiplier %s announces %s positions, but %s found in fluoro group %s"%(
                            multprefix,multiplier,len(positions),",".join([str(i) for i in positions])+"-"+prefixname
                            ))
                
                # Go through every functional group specified by this suffix
                for position in positions:
                    # Store the result
                    fg = {
                        "type":"fluoro",
                        "base":position,
                        }
                    data["fgroups"].append(fg)
            elif prefixname.endswith("chloro"):
                # Chlorine Group
                multprefix = prefixname[:-6] # removes the chloro from the end
                multiplier = parseAlkylMultiplierPrefix(multprefix)
                
                if positions is None:
                    if multiplier>data["main_chain_length"]:
                        raise errors.InvalidMultiplier("Multiplier %s of chloro Group announces %s positions, but only %s available"%(
                                multprefix,multiplier,data["main_chain_length"]
                                ))
                    positions = []
                    # positions were not specified, need to be generated
                    # Alternating from both ends and at opposite ends
                    # Note that this may not always be chemically correct
                    for i in range(multiplier):
                        if i%2==0: # left
                            positions.append(int(i/2)+1)
                        else: # right
                            positions.append(data["main_chain_length"]-(int(i/2)+1))
                
                if multiplier!=len(positions):
                    raise errors.InvalidMultiplier("Multiplier %s announces %s positions, but %s found in chloro group %s"%(
                            multprefix,multiplier,len(positions),",".join([str(i) for i in positions])+"-"+prefixname
                            ))
                
                # Go through every functional group specified by this suffix
                for position in positions:
                    # Store the result
                    fg = {
                        "type":"chloro",
                        "base":position,
                        }
                    data["fgroups"].append(fg)
            elif prefixname.endswith("bromo"):
                # Bromine Group
                multprefix = prefixname[:-5] # removes the bromo from the end
                multiplier = parseAlkylMultiplierPrefix(multprefix)
                
                if positions is None:
                    if multiplier>data["main_chain_length"]:
                        raise errors.InvalidMultiplier("Multiplier %s of bromo Group announces %s positions, but only %s available"%(
                                multprefix,multiplier,data["main_chain_length"]
                                ))
                    positions = []
                    # positions were not specified, need to be generated
                    # Alternating from both ends and at opposite ends
                    # Note that this may not always be chemically correct
                    for i in range(multiplier):
                        if i%2==0: # left
                            positions.append(int(i/2)+1)
                        else: # right
                            positions.append(data["main_chain_length"]-(int(i/2)+1))
                
                if multiplier!=len(positions):
                    raise errors.InvalidMultiplier("Multiplier %s announces %s positions, but %s found in bromo group %s"%(
                            multprefix,multiplier,len(positions),",".join([str(i) for i in positions])+"-"+prefixname
                            ))
                
                # Go through every functional group specified by this suffix
                for position in positions:
                    # Store the result
                    fg = {
                        "type":"bromo",
                        "base":position,
                        }
                    data["fgroups"].append(fg)
            elif prefixname.endswith("iodo"):
                # Iodine Group
                multprefix = prefixname[:-4] # removes the iodo from the end
                multiplier = parseAlkylMultiplierPrefix(multprefix)
                
                if positions is None:
                    if multiplier>data["main_chain_length"]:
                        raise errors.InvalidMultiplier("Multiplier %s of iodo Group announces %s positions, but only %s available"%(
                                multprefix,multiplier,data["main_chain_length"]
                                ))
                    positions = []
                    # positions were not specified, need to be generated
                    # Alternating from both ends and at opposite ends
                    # Note that this may not always be chemically correct
                    for i in range(multiplier):
                        if i%2==0: # left
                            positions.append(int(i/2)+1)
                        else: # right
                            positions.append(data["main_chain_length"]-(int(i/2)+1))
                
                if multiplier!=len(positions):
                    raise errors.InvalidMultiplier("Multiplier %s announces %s positions, but %s found in iodo group %s"%(
                            multprefix,multiplier,len(positions),",".join([str(i) for i in positions])+"-"+prefixname
                            ))
                
                # Go through every functional group specified by this suffix
                for position in positions:
                    # Store the result
                    fg = {
                        "type":"iodo",
                        "base":position,
                        }
                    data["fgroups"].append(fg)
            else:
                # Just assume it is an alkyl
                try:
                    prefixname = prefixname[:-2] # removes the -yl from the end
                    alkyl_n,multprefix = parseAlkanePrefix(prefixname,True)
                    multiplier = parseAlkylMultiplierPrefix(multprefix)
                    pname = getAlkanePrefix(alkyl_n)+"yl"
                    
                    if multiplier!=len(positions):
                        raise errors.InvalidMultiplier("Multiplier %s announces %s positions, but %s found in alkyl group %s"%(
                                multprefix,multiplier,len(positions),",".join(positions)+"-"+prefixname+"yl"
                                ))
                    
                    # Go through every functional group specified by this prefix
                    for position in positions:
                        # Store the result
                        if position>data["main_chain_length"]:
                            raise errors.BaseAtomOutOfRangeError("Tried to add %syl to Carbon #%s, but main chain is only %s Carbons long"%(
                                    prefixname,position,data["main_chain_length"]
                                    ))
                        fg = {
                            "type":"alkyl",
                            "base":position,
                            "alkyl_length":alkyl_n,
                            "alkyl_name":pname,
                            }
                        data["fgroups"].append(fg)
                except errors.InvalidMultiplier:
                    raise
                except Exception:
                    raise errors.InvalidPrefixError("Invalid prefix %s"%prefixname)
        
        # Convert the suffix to functional groups
        for positions,suffixname in data["suffixes"]:
            if suffixname.endswith("ol"):
                # TODO: support a larger variety of syntaxes
                multprefix = suffixname[:-2] # removes the -ol from the end
                multiplier = parseAlkylMultiplierPrefix(multprefix)
                
                if positions is None:
                    if multiplier>data["main_chain_length"]:
                        raise errors.InvalidMultiplier("Multiplier %s announces %s positions, but only %s available"%(
                                multprefix,multiplier,data["main_chain_length"]
                                ))
                    positions = []
                    # positions were not specified, need to be generated
                    # Alternating from both ends and at opposite ends
                    for i in range(multiplier):
                        if i%2==0: # left
                            positions.append(int(i/2)+1)
                        else: # right
                            positions.append(data["main_chain_length"]-(int(i/2)+1))
                
                if multiplier!=len(positions):
                    raise errors.InvalidMultiplier("Multiplier %s announces %s positions, but %s found in alkanol group %s"%(
                            multprefix,multiplier,len(positions),",".join([str(i) for i in positions])+"-"+suffixname
                            ))
                
                # Go through every functional group specified by this suffix
                for position in positions:
                    # Store the result
                    fg = {
                        "type":"hydroxyl",
                        "base":position,
                        }
                    data["fgroups"].append(fg)
            elif suffixname.endswith("al"):
                raise errors.UnsupportedFeatureError("Aldehydes are currently not supported")
            elif suffixname.endswith("one"):
                raise errors.UnsupportedFeatureError("Ketones are currently not supported")
            else:
                if suffixname=="":
                    raise errors.InvalidFormulaError("Expected suffix after dash, but end-of-formula found")
                raise errors.UnsupportedFeatureError("Suffix type '%s' is currently not supported"%suffixname)
    def i2s_stage3(self,data):
        # Stage 3
        # Create the functional groups, but leave them disconnected from the main chain
        
        #  First stage that needs the structural formula object
        data["struct"] = structural.StructuralNotation()
        
        # Iterate through all functional groups and call the appropriate constructor
        # The constructors do not return anything, but they will modify the fg dict and add the bondinfo key
        for fg in data["fgroups"]:
            if fg["type"] in self.fgroups:
                f = self.fgroups[fg["type"]]
                # TODO: add exception catching to functional group creation
                f(fg,data)
            else:
                raise errors.UnsupportedGroupError("Functional groups of type '%s' cannot be converted to structures yet"%fg["type"])
    def i2s_stage4(self,data):
        # Stage 4
        # Create the main chain
        
        # Add carbons to struct
        carbons = []
        for i in range(data["main_chain_length"]):
            carbons.append(data["struct"].addCarbon(name="C%s"%(i+1)))
        
        # Connect them together, from left to right
        for c in carbons:
            if carbons.index(c)==0:
                continue
            c.bindToAtom(carbons[carbons.index(c)-1])
        
        data["carbons"]=carbons
    def i2s_stage5(self,data):
        # Stage 5
        # Connect the functional groups to the main chain
        for fg in data["fgroups"]:
            # Base Carbon Number, Connecting atom, Bond count, Bond data base, Bond data connecting
            base,conn,n,bdata,cdata = fg["bondinfo"]
            data["carbons"][base-1].bindToAtom(conn,n,bdata,cdata)
    def i2s_stage6(self,data):
        # Stage 6
        # Fill the molecule with hydrogen
        data["struct"].fillWithHydrogen()
    
    def fg_alkyl(self,fg,data):
        # Creates an alkyl group of the given length, but does not yet connect it with the main chain
        prev = None
        first = None
        for i in range(fg["alkyl_length"]):
            name = "C %s-%s #%s"%(fg["base"],fg["alkyl_name"],i+1)
            c = data["struct"].addCarbon(name=name)
            if prev is not None:
                c.bindToAtom(prev,sdata={"reason":"%s Group generated from IUPAC Name"%fg["alkyl_name"]})
            else:
                first = c
            prev = c
        
        # Create bondinfo
        base = fg["base"]
        conn = first
        n = 1
        bdata = {"reason":"%s Group generated from IUPAC Name "%fg["alkyl_name"]}
        cdata = bdata
        fg["bondinfo"] = base,conn,n,bdata,cdata
    def fg_hydroxyl(self,fg,data):
        # Creates an hydroxyl group, but does not yet connect it with the main chain
        # Create oxygen
        o = data["struct"].addOxygen(name="O@%s"%fg["base"])
        # Create Hydrogen
        h = data["struct"].addHydrogen(name="H@%s"%fg["base"])
        o.bindToAtom(h,sdata={"reason":"Hydroxyl Group generated from IUPAC Name"})
        
        # Create bondinfo
        base = fg["base"]
        conn = o
        n = 1
        bdata = {"reason":"Hydroxyl Group generated from IUPAC Name"}
        cdata = bdata
        fg["bondinfo"] = base,conn,n,bdata,cdata
    def fg_amino(self,fg,data):
        # Creates an amino group, but does not yet connect it with the main chain
        
        # Create Nitrogen
        n = data["struct"].addNitrogen(name="N@%s"%fg["base"])
        
        # Create Hydrogen 1
        h1 = data["struct"].addHydrogen(name="H1@%s"%fg["base"])
        n.bindToAtom(h1,sdata={"reason":"Amino Group generated from IUPAC Name"})
        
        # Create Hydrogen 2
        h2 = data["struct"].addHydrogen(name="H2@%s"%fg["base"])
        n.bindToAtom(h2,sdata={"reason":"Amino Group generated from IUPAC Name"})
        
        # Create bondinfo
        base = fg["base"]
        conn = n
        n = 1
        bdata = {"reason":"Amino Group generated from IUPAC Name"}
        cdata = bdata
        fg["bondinfo"] = base,conn,n,bdata,cdata
    def fg_fluoro(self,fg,data):
        # Creates a fluoro group, but does not yet connect it with the main chain
        # Create Fluorine
        f = data["struct"].addFluorine(name="F@%s"%fg["base"])
        
        # Create bondinfo
        base = fg["base"]
        conn = f
        n = 1
        bdata = {"reason":"Fluoro Group generated from IUPAC Name"}
        cdata = bdata
        fg["bondinfo"] = base,conn,n,bdata,cdata
    def fg_chloro(self,fg,data):
        # Creates a chloro group, but does not yet connect it with the main chain
        # Create Chlorine
        c = data["struct"].addChlorine(name="Cl@%s"%fg["base"])
        
        # Create bondinfo
        base = fg["base"]
        conn = c
        n = 1
        bdata = {"reason":"Chloro Group generated from IUPAC Name"}
        cdata = bdata
        fg["bondinfo"] = base,conn,n,bdata,cdata
    def fg_bromo(self,fg,data):
        # Creates a bromo group, but does not yet connect it with the main chain
        # Create Bromine
        b = data["struct"].addBromine(name="Br@%s"%fg["base"])
        
        # Create bondinfo
        base = fg["base"]
        conn = b
        n = 1
        bdata = {"reason":"Bromo Group generated from IUPAC Name"}
        cdata = bdata
        fg["bondinfo"] = base,conn,n,bdata,cdata
    def fg_iodo(self,fg,data):
        # Creates a iodo group, but does not yet connect it with the main chain
        # Create Iodine
        i = data["struct"].addIodine(name="I@%s"%fg["base"])
        
        # Create bondinfo
        base = fg["base"]
        conn = i
        n = 1
        bdata = {"reason":"Iodo Group generated from IUPAC Name"}
        cdata = bdata
        fg["bondinfo"] = base,conn,n,bdata,cdata
    
    # Save to String Methods
    def dumpAsSMILES(self):
        return self.asStructuralFormula().dumpAsSMILES()
    
    def dumpAsInChI(self):
        return self.asStructuralFormula().dumpAsInChI()
    
    # Load from String Methods
    @classmethod
    def loadsFromSMILES(cls,data):
        return structural.StructuralNotation.loadsFromSMILES(data).asIUPACName()
    
    @classmethod
    def loadsFromInChI(cls,data):
        return structural.StructuralNotation.loadsFromInChI(data).asIUPACName()
    
    # Magic Methods
    def __repr__(self):
        return "<IUPACNotation(name='%s')>"%self.name
    def __str__(self):
        return self.name
    
    def __eq__(self,other):
        return self.__class__==other.__class__ and \
               self.name==other.name

def getNumericalPrefix(n):
    # Returns the numerical prefix associated with the given number
    # Note that this function was designed for use with alkane, this may influence some prefixes, e.g. 1=meth etc.
    
    if n in SPECIAL_PREFIXES:
        return SPECIAL_PREFIXES[n]
    elif n<10:
        # Single digit
        prefix = BASE_PREFIXES_1[n]#.rstrip("a")
        return prefix
    elif n<10000:
        out = ""
        
        # last digit
        # multiplier 1, ones
        digit1 = int(str(n)[-1])
        pre1 = BASE_PREFIXES_1[digit1]
        out+=pre1
        
        # 2nd digit from the right
        # multiplier 10, tens
        digit2 = int(str(n)[-2])
        pre2 = BASE_PREFIXES_10[digit2]
        if digit2==2 and len(out)>0 and out[-1] in ["a","e","i","o","y"]: # ends with a vowel
            pre2 = pre2.lstrip("i")
        out+=pre2
        
        # 3rd digit from the right
        # multiplier 100, hundreds
        if n>=100:
            digit3 = int(str(n)[-3])
            pre3 = BASE_PREFIXES_100[digit3]
            out+=pre3
        
        # 4th digit from the right
        # multiplier 1000, thousands
        if n>=1000:
            digit4 = int(str(n)[-4])
            pre4 = BASE_PREFIXES_1000[digit4]
            out+=pre4
        
        return out#.rstrip("a")
    else:
        raise NotImplementedError("Numerical Prefixes for numbers higher than 9999 are not implemented")

def parseNumericalPrefix(s_in,return_leftover=False):
    s_in = s_in.lower()
    if not return_leftover:
        if s_in in SPECIAL_PREFIXES.inv:
            return SPECIAL_PREFIXES.inv[s_in]
    elif return_leftover:
        for n,prefix in SPECIAL_PREFIXES.items():
            if s_in.endswith(prefix) and prefix!="":
                s_in = s_in[:-len(prefix)]
                return n,s_in
    
    if s_in=="":
        raise errors.InvalidPrefixError("Prefix cannot be an empty string")
    
    n = 0
    
    # Thousands
    if s_in.endswith("lia"): # suffix for thousands
        for val,prefix in BASE_PREFIXES_1000.items():
            if s_in.endswith(prefix) and prefix!="":
                n+=val*1000
                s_in = s_in[:-len(prefix)] # cuts off the thousands
    
    # Hundreds
    if s_in.endswith("cta"):
        for val,prefix in BASE_PREFIXES_100.items():
            if s_in.endswith(prefix) and prefix!="":
                n+=val*100
                s_in = s_in[:-len(prefix)] # cuts off the hundreds
    
    # Tens
    if s_in.endswith("aconta") or s_in.endswith("deca") or s_in.endswith("icosa"):
        for val,prefix in BASE_PREFIXES_10.items():
            if s_in.endswith(prefix) and prefix!="":
                n+=val*10
                s_in = s_in[:-len(prefix)] # cuts off the tens
    if s_in.endswith("cosa"):
        n+=20
        s_in = s_in[:-4]
    if s_in.endswith("tr"): # fix for tricosa
        s_in+="i"
    
    # Ones
    for val,prefix in BASE_PREFIXES_1.items():
        if s_in.endswith(prefix) and prefix!="":
            n+=val
            s_in = s_in[:-len(prefix)] # cuts off the ones
    
    if not return_leftover and s_in!="":
        raise errors.InvalidPrefixError("Prefix could not be fully parsed, %s remained"%s_in)
    
    if not return_leftover:
        return n
    elif return_leftover:
        return n,s_in

def getAlkanePrefix(n):
    if n in SPECIAL_ALKANE_PREFIXES:
        return SPECIAL_ALKANE_PREFIXES[n]
    return getNumericalPrefix(n).rstrip("a")

def parseAlkanePrefix(s_in,return_leftover=False):
    s_in = s_in.lower()
    if not return_leftover:
        if s_in in SPECIAL_ALKANE_PREFIXES.inv:
            return SPECIAL_ALKANE_PREFIXES.inv[s_in]
    elif return_leftover:
        for n,prefix in SPECIAL_ALKANE_PREFIXES.items():
            if s_in.endswith(prefix) and prefix!="":
                s_in = s_in[:-len(prefix)]
                return n,s_in
    return parseNumericalPrefix(s_in+"a",return_leftover)

def getAlkylMultiplierPrefix(n,n2=0):
    if n==1:
        return ""
    out = getNumericalPrefix(n)
    if n2 != 0:
        if parseAlkanePrefix(out+getAlkanePrefix(n2),True)[0]!=n2:
            # Happens e.g. with n=3 and n2=10
            # tri- and dec- may not be able to be distinguished, so -kis is added to the end
            # Exceptions apply for bis/tris
            # See http://www.acdlabs.com/iupac/nomenclature/93/r93_55.htm#r_0_1_4_2
            if out=="bi":
                out="bis"
            elif out=="tri":
                out="tris"
            else:
                out+="kis"
    return out

def parseAlkylMultiplierPrefix(s_in,return_leftover=False):
    if s_in=="":
        if return_leftover:
            return 1,""
        else:
            return 1
    elif s_in.endswith("bis") or s_in.endswith("tris"):
        s_in = s_in[:-1] # remove extra -s
    elif s_in.endswith("kis"):
        s_in = s_in[:-3] # remove extra -kis
    return parseNumericalPrefix(s_in,return_leftover)
