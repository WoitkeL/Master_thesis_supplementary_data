# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:08:37 2023

@author: linus
"""
import pulp
import regex as re
def get_compartments(lOR,true_reactions, extraspecies=""):
    
    #function to delete subsets
    def delete_smaller_sets(candidate_set):
        candidate_set_remover=set()
        for element in candidate_set:
            for element2 in candidate_set:
                if element!=element2:
                    if element.issubset(element2):
                        candidate_set_remover.add(element)
                        break
        candidate_set.difference_update(candidate_set_remover)
        return(candidate_set)
    
    #function to split speciesset for a not proper reaction
    def splitt_speciesset(element, reaction):
        extender=set()
        for reactant in reaction.listOfReactants:
            #for each reaction we remove each reactant once and add the
            #remaining set to the extender, which is then returned
            extender_element=set(element.copy())
            extender_element.remove(reactant)
            extender.add(frozenset(extender_element))
        return(extender)
    
    
    lOR_gc=lOR[:]
    LOR_dictionary={}
    setOfSpecies_gc=set()
    candidate_set=set()
    inactive_reactions=set()
    true_reactions_reac=set()

    #create LOR-dictionary to link name to reaction-class-objekt
    for reaction in lOR_gc:
        LOR_dictionary[reaction.defined_name]=reaction
    
    #set species
    for reaction in true_reactions:
        setOfSpecies_gc.update(LOR_dictionary[reaction].listOfReactants)  
        setOfSpecies_gc.update(LOR_dictionary[reaction].listOfProducts)
    setOfSpecies_gc.update(extraspecies)
    
    #remove uninvolved reactions
    for reaction in lOR_gc:
        if not set(reaction.listOfReactants).issubset(setOfSpecies_gc):
            lOR_gc.remove(reaction)
            continue

    #create list of active and inactive reactions
    for reaction in lOR_gc:
        if reaction.defined_name in true_reactions:
            true_reactions_reac.add(reaction)
        else:
            inactive_reactions.add(reaction)
    
    '''dittrich anderer name MRC: heaviest compartments, all species so you cant add any whitout acitvating inactive reac'''
    #start MRCs with set of all species as candidate set
    candidate_set.add(frozenset(setOfSpecies_gc.copy()))
    counter=0
    #iterate over every inactive reaction to separate the reactants
    for inactive_reaction in inactive_reactions:
        counter+=1      
        candidate_set_copy=candidate_set.copy()
        #iterate over all candidate sets to seperate reactants of inactive reactions
        for element_set in candidate_set_copy:
            if set(inactive_reaction.listOfReactants).issubset(element_set):

                #add all smaller candidates, which obey the absence of the support of the current inactive reaction   
                candidate_set.update(splitt_speciesset(element_set,inactive_reaction))
                
                #delete all old candidates, which contain the support of the current inactive reaction
                candidate_set.remove(element_set)
 

        if counter%5==0:
            candidate_set= delete_smaller_sets(candidate_set)

    
    #delete subsets in set
    candidate_set= delete_smaller_sets(candidate_set)
    #check for closed
    counter=0

    
    #iterate over every active reaction to check for closedness of every compartment
    nochange=False
    while nochange==False:
        nochange_current=True
        for reaction in true_reactions_reac:
            
            counter+=1
            candidate_set_copy=candidate_set.copy()
            for element_set in candidate_set_copy:    
                if set(reaction.listOfReactants).issubset(element_set):
    
                    if not set(reaction.listOfProducts).issubset(element_set):
                        #add all smaller candidates, which obey the absence of the support of the current non-closed reaction 
                        candidate_set.update(splitt_speciesset(element_set, reaction))
                        nochange_current=False
                        #delete all old candidates, which contain the support of the current non-closed reaction 
                        candidate_set.remove(element_set)
        if nochange_current==True:
            nochange=True
      #  if counter%50==0:
       #     candidate_set= delete_smaller_sets(candidate_set)

            
    candidate_set= delete_smaller_sets(candidate_set)

    #setup for LP_solver
    candidate_reactions_dictionary={}
    reaction_compartment_dictionary={}
    species_compartment_dictionary={}
    for reaction in true_reactions:
        reaction_compartment_dictionary[reaction]=[]
    for species in setOfSpecies_gc:
        species_compartment_dictionary[species]=[]
    for element in candidate_set:
        reactionlist=[]
        for reaction in true_reactions_reac:
            if set(reaction.listOfReactants).issubset(element):
                reactionlist.extend(reaction.defined_name)
                reaction_compartment_dictionary[reaction.defined_name].append(element)
        for species in element:
            species_compartment_dictionary[species].append(element)

    return((setOfSpecies_gc, candidate_set, reaction_compartment_dictionary, species_compartment_dictionary, true_reactions))


def get_min_compartments(setOfSpecies,candidate_set, reac_dic, species_dic,true_reactions):

    iterator=set()
    compartment_shortform={}
    compartment_shortform_return={}
    i=0
    #setup of dictionaries used constraint construction
    for ding in candidate_set:
        #maps an index to each of the candidates in both directions
        compartment_shortform[i]= ding
        compartment_shortform_return[str(ding)]= i
        i+=1
        iterator.add(str(ding))
        
    # setup variables including compartments as instance of existence of a speciesset
    #also setup of reaction and species variable, which have to be present in at least one active compartment
    Reac=pulp.LpVariable.dicts("reaction_value",true_reactions,lowBound=0, cat="Integer")
    species=pulp.LpVariable.dicts("species_value", frozenset(setOfSpecies),lowBound=0,cat="Integer")
    compartments=pulp.LpVariable.dicts("compartment", compartment_shortform.keys(),lowBound=0, upBound=1,cat="Integer" )
    
    #define lp objective function for minimization of the number of compartments active
    solve_min_comp=pulp.LpProblem("Solver Min Compartment",pulp.LpMinimize)
    objective_function=pulp.lpSum(compartments.values())  
    solve_min_comp+=objective_function  
    
    #setup of reaction and species in relation to their compartments, in which the occur
    for reaction in true_reactions:
        solve_min_comp+=pulp.lpSum(compartments[compartment_shortform_return[str(b)]] for b in reac_dic[reaction])>=1
    for single_species in setOfSpecies:
        solve_min_comp+=pulp.lpSum(compartments[compartment_shortform_return[str(b)]] for b in species_dic[single_species])>=1
#     print(solve_min_comp)
    solve_min_comp.solve()
#     print(pulp.LpStatus[solve_min_comp.status])

    
            
            
    #variable extraction and output
    solutionlist=[]
    species_over_list=[]
    pattern1 = re.compile("compartment")
    pattern2 = re.compile("species_value")
#     print(true_reactions)
    for variable in solve_min_comp.variables():
#        print("{} = {}".format(variable.name,variable.varValue))
        if pattern1.match(variable.name) and variable.varValue==1:
            solutionlist.append(compartment_shortform[int(variable.name[12:])])
        if variable.varValue==None:
            print(variable.name)
            continue
        if variable.varValue>=2 and pattern2.match(variable.name):
            species_over_list.append(variable.name[14:])
#     print(species_over_list)
    return(solutionlist)

"######################################################################################################################"
