# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 16:54:46 2023

@author: linus
"""

import time
import libsbml
import os
"""Reaction class represents instances of reaction with all their used informations
each reaction is initialized with the following parameters: name, reactants, products and their affiliated
stoichiometric parameters"""
"""Reaction class represents instances of reaction with all their used informations
each reaction is initialized with the following parameters: name, reactants, products and their affiliated
stoichiometric parameters"""
class Reaction:
    """Creates an instance of a Reaction and appends this instance to the listOfReactions"""
    def __init__(self, name, reactants, products, extract_stoich_rea=[], extract_stoich_prod=[]):
        
        self.closed=True                    # properties of closure
        
        self.defined_name=name
        self.reac_stoich=extract_stoich_rea
        self.prod_stoich=extract_stoich_prod
        self.listOfReactants=reactants
        self.listOfProducts=products
        self.always=(len(self.listOfReactants)==0) 
        self.reversible=False
        #listOfReactions.append(self)
    def __repr__(self):
        #output=self.defined_name + ":"+ self.listOfReactants+" "+self.reac_stoich+"->"+self.listOfProducts+":"+self.prod_stoich
        output=self.defined_name
        return(output)

#Reaktion-Klasse bildet Instanzen von minimalen Reaktionsräumen einer Reaktion und trägt den Namen dieser.
#spezielles set von reaktionen, der dem Abschluss einer bestimmten Reaktion entspricht
#single-reaction-closure
#elementarer reactionsabschluss
#elementary reaction closure represents instance of the minimal space of a reaction in regard to triggering other reactions.

def generate_closure_for_species(ele,LOR):
    x=ERC(LOR,solospecies=ele)
    return(x)

def generate_closure_for_reaction(ele,LOR):
    x=ERC(LOR,reaction=ele)
    return(x)


    
def create_closures(LOR, species=None, Timer_in_sec=False):
    time_start=False
    dict_output={}
    if type(Timer_in_sec)==int:
        time_start=time.process_time()
    
    if species==None:
        for ele in LOR:
            dict_output[ele.defined_name]=generate_closure_for_reaction(ele,LOR)
            if time_start:
                if time.process_time() > Timer_in_sec + time_start:
                    return("error")
       # if type(Timer_in_sec)==int:    
        #    dict_output[LOR[0].defined_name].timer=time.process_time()-time_start 
    else:
        for specie in species:
            dict_output[specie]=generate_closure_for_species(specie,LOR)
            if time_start:
                if time.process_time() > Timer_in_sec + time_start:
                    return("error")    
    return(dict_output)
    
class ERC:
    # initialization includes creating the list of species and the list of reactions of the ERC and adding the appropriate
    # items of the starting reaction / solospecies
    def __init__(self, lOR_solve, reaction=[], solospecies=[]):
        if reaction !=[]:
            self.defined_name=reaction.defined_name
#             ERC_dict[self.defined_name]=self
            self.reactions=[reaction]
            self.species=set()
            #self.timer=None
            for element in reaction.listOfReactants:
                self.species.add(element)
            for element in reaction.listOfProducts:
                self.species.add(element)
            self.eRC_aufstellung(lOR_solve)
            
        else:
            if solospecies==[]:
                self.species=set()
                self.defined_name="empty"
            else:
                self.defined_name=solospecies
                self.species={solospecies}
            self.reactions=[]
#             ERC_dict[self.defined_name]=self
            self.eRC_aufstellung(lOR_solve)
            
#function iterates over all reactions and checks if a reaction, which is not yet part of the ERC, is supported
#by the species set
#if this is the case, the reaction and its products are added to the ERC, and the check of the reactions continues
#if the species set is checked for every reaction without adding species, the function terminates
  
    def eRC_aufstellung(self, lOR_solve):
        itera= [(i, False) for i in range(len(lOR_solve))]
        checker={key: value for (key, value) in (itera)}

        last_match=(len(lOR_solve)-1)


 #       for checkreaction_index in range(len(lOR_solve)):
  #          filtered_ERC_reactions=lOR_solve.copy()
   #         if lOR_solve[checkreaction_index].always==True:
    #            del filtered_ERC_reactions[checkreaction_index]
        
        #infinite loop until broken by return command
        while True:
            #iterate over all reactions
            for checkreaction_index in range(len(lOR_solve)):
                #only check reactions, which are not a part already
                if not lOR_solve[checkreaction_index] in self.reactions:
                    #reaction supported
                    
                    if set(lOR_solve[checkreaction_index].listOfReactants).issubset(self.species):
                        #add reaction and species to ERC
                        
                        self.reactions.append(lOR_solve[checkreaction_index])
                        self.species.update(lOR_solve[checkreaction_index].listOfProducts)
                        
                        last_match=checkreaction_index
                        #last match is safed and loop is continued, so the loop only return after a whole loop over all reactions 
                        continue
                    
                    #also try for reversible reactions with swapped list of reactions and products
                    elif lOR_solve[checkreaction_index].reversible==True:
                        if set(lOR_solve[checkreaction_index].listOfProducts).issubset(self.species):
                            self.reactions.append(lOR_solve[checkreaction_index])
                            self.species.update(lOR_solve[checkreaction_index].listOfReactants)
                            last_match=checkreaction_index
                            continue
                    
                if checkreaction_index==last_match:
                    return
        
#Es gibt 2 Möglichkeiten, die Reaktionen in das Programm einzuspeisen. Wird kein Dateipfad mit einer SBML-Datei gegeben, 
#greift es die im Programmtext manuell gegebenen Reaktionen ab.
#uses the 2 options to get the reaction network. If no path to an SBML-file is given, it takes the  manually alterable
#reaction network down below
def getReaction(path="example", consider_reverse=True, consider_constant=False, consider_init_ammount=False,alternative_reverse=False):
    reverse_thing=[]
    listOfReactions_gR = []
    if path=="example":
        #listOfReactions_gR.append(Reaction("r1",["b1"],["b1" , "a1"],[1],[1,1]))
        #listOfReactions_gR.append(Reaction("r2",["b2"],[ "b2" , "a2"],[1],[1,1]))
        #listOfReactions_gR.append(Reaction("r3",["b3"],[ "b3" , "a3"],[1],[1,1]))
        #listOfReactions_gR.append(Reaction("r4",["b2","b3"],["b3"],[1,1],[2]))
        #listOfReactions_gR.append(Reaction("r5",["b1","a2"],[],[1,1],[]))
        #listOfReactions_gR.append(Reaction("r6",["b1","b3"],["b2"],[1,1],[2]))
        #listOfReactions_gR.append(Reaction("r7",["a1","a2"],[],[1,1],[]))
        #listOfReactions_gR.append(Reaction("r8",["a1","a2","a3"],["d1"],[1,1,1],[1]))
        #listOfReactions_gR.append(Reaction("r9",["d1","b1"],["b1","a2"],[1,1],[2,1]))
        
        
        #listOfReactions_gR.append(Reaction("r10",["b3","a1","a3"],["b3"],[1,1,1],[2]))
        listOfReactions_gR.append(Reaction("r1",[],["h"],[],[1]))
        listOfReactions_gR.append(Reaction("r2",["h","v"],["in","v"],[1,1],[1,1]))
        listOfReactions_gR.append(Reaction("r3",["h","v"],["h","vb"],[1,1],[1,1]))
        listOfReactions_gR.append(Reaction("r4",["vb"],["m"],[1],[1]))
        listOfReactions_gR.append(Reaction("r5",["p"],["v"],[1],[1]))
        listOfReactions_gR.append(Reaction("r6",["m"],["m" , "s"],[1],[1,1]))
        listOfReactions_gR.append(Reaction("r7",["s"],[],[1],[]))
        listOfReactions_gR.append(Reaction("r8",["in","s"],["s"],[1,1],[1]))
        listOfReactions_gR.append(Reaction("r9",["v","s"],["s"],[1,1],[1]))
        
    else:
#SBML-Reader from LibSBML-Package, is initialized
        reader = libsbml.SBMLReader()
        
#various errors are captured
        if not os.path.isfile(path):
            print("no file in path found")
            return()
        if reader == None:
            print("no object created")
        
        doc_extract = reader.readSBMLFromFile(path)
        if doc_extract.getNumErrors() > 0:
            if doc_extract.getError(0).getErrorId() == libsbml.XMLFileUnreadable:
                print("XMLFileUnreadable")
            elif doc_extract.getError(0).getErrorId() == libsbml.XMLFileOperationError:
                print("XMLFileOperationError")
                
                
#in the following section, functions and classes of the libsbml-package are used:
#model is generated:

        model_extr=doc_extract.getModel()
        listOfSpecies_gR=[]
        #model_name= getSBMLmodName(model_extr)
        if model_extr is None:
            return(None, None)
        model_name="no name"
        model_name= model_extr.getName()
#reactionlist is extracted from model:
        listrea_extract=model_extr.getListOfReactions()
#species list is extracted from model:
        reverse=[]
        #for schleife für jede einzelne Reaktion
        for i in range (len(listrea_extract)):
            name=listrea_extract.get(i).getId()
            
#sbml.ListOfReactions.get(x) extracts SBML.Reactions with the attributes ListOfReactants and ListOfProducts            
            extract_rea_list=listrea_extract.get(i).getListOfReactants()
            extract_pro_list=listrea_extract.get(i).getListOfProducts()
            extract_rea_list_species=[]
            extract_pro_list_species=[]
            extract_rea_list_species_stoich=[]
            extract_pro_list_species_stoich=[]
#from these lists we get the species with their respective stoichiometric factor            
            for j in range (len(extract_rea_list)):
                extract_rea_list_species.append(extract_rea_list[j].getSpecies())
                extract_rea_list_species_stoich.append(extract_rea_list[j].getStoichiometry())
                
            for j in range (len(extract_pro_list)):
                extract_pro_list_species.append(extract_pro_list[j].getSpecies())
                extract_pro_list_species_stoich.append(extract_pro_list[j].getStoichiometry())
#reactions are then initialized with the self implemented reaction class       
            if listrea_extract.get(i).getReversible():
                reverse.append(name)
                
            x= Reaction(name,extract_rea_list_species,extract_pro_list_species,extract_rea_list_species_stoich,extract_pro_list_species_stoich)
            listOfReactions_gR.append(x)
        if consider_reverse:
            for i in listOfReactions_gR:
    #             print(i.defined)
                if i.defined_name in reverse:
                    has_rea=0
                    for j in listOfReactions_gR:
                        if i==j:
                            continue
                        
                        if set(j.listOfReactants)==set(i.listOfProducts):
         
                            if set(j.listOfProducts)==set(i.listOfReactants):
                                first_works=1
                                for k in range(len(i.listOfReactants)):
                                    first_works=0
                                    ind=j.listOfProducts.index(i.listOfReactants[k])
                                    if i.reac_stoich[k]==j.prod_stoich[ind]:
                                        first_works=1
                                    else: break
                                if first_works==1:
                                    if len(i.listOfProducts)==0:
                                        has_rea=1
                                    for k in range(len(i.listOfProducts)):
                                        ind=j.listOfReactants.index(i.listOfProducts[k])
                                        if i.prod_stoich[k]==j.reac_stoich[ind]:
                                            has_rea=1
                                            print("skip_reverse")
                    if not has_rea:
                        if alternative_reverse==True:
                            i.reversible=True
                        else:
                            listOfReactions_gR.append(Reaction(i.defined_name+str("_reverse"),i.listOfProducts,i.listOfReactants,i.prod_stoich,i.reac_stoich))
                        #reverse_thing.append(i.defined_name)
                        #i.reversible=True
#        if show_reactions:
 #           for reaction123 in listOfReactions_gR:
  #              print(reaction123.defined_name, end=":")
   #             print(reaction123.listOfReactants, end=" ")
    #            print(reaction123.reac_stoich, end="->")
     #           print(reaction123.listOfProducts, end=":")
      #          print(reaction123.prod_stoich)
       # for reaction in listOfReactions_gR:
        #        if reaction.always==True:
         #           for ele in reaction.listOfProducts:
          #              add_rea=True
           #             for reaction2 in listOfReactions_gR:
            #                if reaction2.listOfProducts==[]:
             #                   if ele in reaction2.listOfReactants:
              #                      add_rea=False
               #         if add_rea:
                #            print("home")
                 #           listOfReactions_gR.append(Reaction(str(reaction.defined_name+str("outfl")+str(ele)),[str(ele)],[],[1],[0]))
        listspec_extract=model_extr.getListOfSpecies()
        inflows=[]
        for reaction in listOfReactions_gR:
            if reaction.always==True:
                inflows.extend(reaction.listOfProducts)
        dicto={}    
        for i in range(len(listspec_extract)):
            name_species=listspec_extract.get(i).getId()
            kaleta_name_species=listspec_extract.get(i).getName()
            dicto[name_species]=kaleta_name_species
            if consider_init_ammount and listspec_extract.get(i).getInitialAmount()>0:
                print("init_ammount_inflow")
                name_reaction=str(name_species)+str("init_ammount_inflow")
                listOfReactions_gR.append(Reaction(name_reaction,[str(name_species)],[str(name_species)],[1],[2]))
                listOfReactions_gR.append(Reaction(str(name_reaction+str(2)),[str(name_species)],[],[1],[0]))
       #     if consider_constant and listspec_extract.get(i).getConstant()==True:
        #        print("_constant_inflow")
         #       print(name_species)
          #      name_reaction=str(name_species)+str("_constant_inflow")
           #     for reaction in listOfReactions_gR:
            #        if reaction.always==True:
             #           inflows.extend(reaction.listOfProducts)
              #  if name_species not in inflows:
               #     listOfReactions_gR.append(Reaction(name_reaction,[],[str(name_species)],[0],[1]))
                 #   if reaction.listOfProducts==[] and reaction:
                #            listOfReactions_gR.append(Reaction(str(name_reaction+str(2)),[str(name_species)],[],[1],[0]))
        #print(dicto)
    #return(listOfReactions_gR,dicto)
    return(listOfReactions_gR)
 #   if give_reverse:
  #      return(listOfReactions_gR,reverse_thing)
   # else:
    #    return(listOfReactions_gR)    

#species can be given manually as a subset of the species of the reaction network, the default option is the set of species
#it then reduces the reaction network to the only usable ones and checks for closedness
def generate_Network_to_be_analyzed(listOfReactions,add_reaction=None, species=None, change_inflow=False,remove_species_from_reactions=False,remove_reaction=False):
    """functions makes changes in the reaction network. Either by reducing the list to a subset of reactions,
    given a subset of species or by changing reactions, which happen in every compartment, to be distributable."""
    listOfReactions_sS=listOfReactions[:]
    setOfSpecies_sS=set()
    if isinstance(add_reaction, Reaction):
        listOfReactions_sS.append(add_reaction)
#if there is no species given, its gained by iterating over all reactions and updating the set
    if species==None:
        for reaction in listOfReactions_sS:
            setOfSpecies_sS.update(reaction.listOfReactants)  
            setOfSpecies_sS.update(reaction.listOfProducts)

    else:
#checks for correct data type of the species set
        if not type(species) is set:
            print("Speciesset muss Set sein")
            return()
        setOfSpecies_sS = species 
        listOfReactions_iterator= listOfReactions_sS.copy()
#since we iterate over the reactions, which we might remove while the iteration, we create a copy of the list of reactions
#to ensure a correct iteration
        for reaction in listOfReactions_iterator:
            if set(reaction.listOfReactants).issubset(setOfSpecies_sS):
                if set(reaction.listOfProducts).issubset(setOfSpecies_sS):
                    reaction.closed=True
                    
                else:
                    reaction.closed=False
                    #closed attribute for reactions, which generate products, which are not in the species set
                    
            else:
                listOfReactions_sS.remove(reaction)
#reactions, which are impossible to be supoorted are removed
    if change_inflow:
#         enz_number=0
        for reaction in listOfReactions:
            if reaction.always:
#                 reaction.listOfReactants.append("made_up_enzyme_"+str(enz_number))
#                 reaction.reac_stoich.append(1)
#                 reaction.listOfProducts.append("made_up_enzyme_"+str(enz_number))
#                 reaction.prod_stoich.append(1)
#                 enz_number+=1
                reaction.listOfReactants.extend(reaction.listOfProducts)
                reaction.reac_stoich.extend([1]* len(reaction.listOfReactants))
#                 for i in range(len(reaction.prod_stoich)):
#                     reaction.prod_stoich[i]=reaction.prod_stoich[i]
                reaction.prod_stoich=[reaction.prod_stoich[i]+1 for i in range(len(reaction.prod_stoich))]
                                             
                                                                
                        
    if len(listOfReactions_sS)==0:
        raise ValueError('no reactions left in reaction_network')
        
    if bool(remove_reaction)==True:    
        for i in range(len(listOfReactions)):
            if listOfReactions[i].defined_name==remove_reaction:
                del listOfReactions[i]
                break
#     if change_null==True:
    if bool(remove_species_from_reactions)==True:
        for reaction in listOfReactions:
            try:
                
                i= reaction.listOfReactants.index(remove_species_from_reactions)
                if len(reaction.listOfReactants)==1:
                    reaction.always=True
                del reaction.listOfReactants[i]
                del reaction.reac_stoich[i]
            except ValueError:
                pass
            try:
                i= reaction.listOfProducts.index(remove_species_from_reactions)
                
                del reaction.listOfProducts[i]
                del reaction.prod_stoich[i]
            except ValueError:
                pass
    return(listOfReactions_sS,setOfSpecies_sS)

def erc_to_matrix(erc_dict):
    i=len(erc_dict)
    bool_matrix= [[False] * i for j in range(i)]
    j=0
    index_to_reac={}
    reac_to_index={}
    for key in erc_dict:
        index_to_reac[j]=key
        reac_to_index[key]=j
        j+=1
    for key in erc_dict:
        for i in erc_dict[key].reactions:
            bool_matrix[reac_to_index[key]][reac_to_index[i.defined_name]]=True
    return(bool_matrix, index_to_reac)

def sort_second(elem):
    return elem[1]

def ERC_meets_transitivity(ERC_dict):
    bool_ma,index_to_reac= erc_to_matrix(ERC_dict)
    remover=[]
    for i in range(len(bool_ma)):
        for j in range(len(bool_ma)):
            if bool_ma[i][j]and i!=j:
                for k in range(len(bool_ma)):
                    if (bool_ma[j][k]and j!=k):
                        bool_ma[i][k] = False
                        if i!=k:
                            remover.append([i,k])
    remover.sort(key=sort_second)
    remover.reverse()
    for ele in remover: 
        ERC_dict[index_to_reac[ele[0]]].reactions

        oko=ERC_dict[index_to_reac[ele[0]]].reactions
        for i in oko:
            if i.defined_name==index_to_reac[ele[1]]:
                oko2=ERC_dict[index_to_reac[ele[0]]].reactions.index(i)
                del ERC_dict[index_to_reac[ele[0]]].reactions[oko2]
    
    return(ERC_dict)