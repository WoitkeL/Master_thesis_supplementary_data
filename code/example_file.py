# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 16:40:35 2023

@author: linus
"""

from analyse_class import Analyse
from Hasse import analyze_DO,check_ORG_LOR,check_ORG_LOR_expanded
from reaction_ERC import Reaction
import time 
from LP_compartments import get_compartments, get_min_compartments
from iterate_over_DB_multi import get_species_closure

example1=Analyse()
# or initializing with path to sbml file
#path=.../...xml
#example1=Analyse(path)

#you can change the handling of reverse reactions in the following ways:
#example1=Analyse(path,consider_reverse=False)

#or make it one directional: 
#but note that you need to use naive algorithm to draw the hasse diagram
#example1=Analyse(path,alternative_reverse=True)

#you can also create inflow reactions for constant species or speices with a concentration higher than 0:
#example1=Analyse(path,consider_constant=True,consider_init_ammount=True)


#options to change reaction network
#example1.change_network(remove_species_from_reactions="EmptySet")
#example1.change_network(change_inflow=True)

#print reaction network
example1.print_RN()



#get list of ERCs
example1.print_ERC()

#get list of species
example1.get_species()
print(example1.species)



#calculating DOS and SARs respectively

example1.SARs()

#solve for SARs with overproduction
#example1.SARs(second_optimize=True,use_gurobi=False)

#example1.DOs()

#print(example1.SAR_solution)
#print(example1.DO_solution)


# calculating all minimal compartments
#mrc numbers can give number of MRCs before using get_min_compartments()
# mrc_numbers=[]

# for i in example1.SAR_solution:
#     bol,boos=check_ORG_LOR_expanded(i[0],example1.reaction_network)
#     if not bol:
#         solution_gc=get_compartments(example1.reaction_network,i[0])
#         mrc_numbers.append(len(solution_gc[1]))
#         compartments=get_min_compartments(*solution_gc)
#         print(compartments)
        
#draw functions here you can see all optional parameters:
#example1.draw_hasse(solution="SARs",use_naive=False,second_value=False,show_species=False,show_compartments=False,show_new=False,shortform=False)
#this command line is the same as
example1.draw_hasse()

#use this to draw DOs:
#example1.draw_hasse(solution="DOs",use_naive=True,show_compartments=False, second_value=False)

#used to get partial data of exel file
#analyze_dict= analyze_DO(example1.SAR_solution, example1.reaction_network)
#print(analyze_dict)


