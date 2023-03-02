# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 16:41:31 2023

@author: linus
"""
import os
from libsbml import SBMLReader
from reaction_ERC import getReaction, generate_Network_to_be_analyzed, create_closures
from setup_LP_DO_final import setup_LP_Var
from Hasse import create_Hasse3, analyze_DO, create_Hasse
from LP_compartments import get_min_compartments
class Analyse:
    """class to handle all properties of a reaction network analysis. 
    It saves any information in class attributes and pipes them into their advancing/ processing functions. 
    """ 
    def __init__(self, source="example",consider_reverse=True, consider_constant=False, consider_init_ammount=False,alternative_reverse=False):
        if source=="example":
            self.source="from reaction network"
        if type(source)==str and source!="example" :
            self.source=source
            if not os.path.isfile(source):
                raise ValueError('no file in path found')
        options={consider_reverse:consider_reverse, consider_constant:consider_constant, consider_init_ammount:consider_init_ammount}
        #self.reaction_network=getReaction(source,options)
        #print("checko")
        #print(alternative_reverse)
        self.reaction_network=getReaction(source,alternative_reverse=alternative_reverse,consider_reverse=consider_reverse)
#    def get_species(self):
 #       setOfSpecies_sS=set()
  #      for reaction in self.reaction_network:
   #         setOfSpecies_sS.update(reaction.listOfReactants)  
    #        setOfSpecies_sS.update(reaction.listOfProducts)
     #   self.species=setOfSpecies_sS
        
    def get_species(self):
        if self.source!="from reaction network":
            setOfSpecies_sS=[]
            reader = SBMLReader()
            doc_extract = reader.readSBMLFromFile(self.source)
            model_extr=doc_extract.getModel()
            listspec_extract=model_extr.getListOfSpecies()
            for i in range(len(listspec_extract)):
                name_species=listspec_extract.get(i).getId()
                setOfSpecies_sS.append(name_species)

        else:
            setOfSpecies_sS=set()
            for reaction in self.reaction_network:
                setOfSpecies_sS.update(reaction.listOfReactants)  
                setOfSpecies_sS.update(reaction.listOfProducts)
        self.species=setOfSpecies_sS
        
    def get_name(self):
        if self.source!="from reaction network":
            reader = SBMLReader()
            doc_extract = reader.readSBMLFromFile(self.source)
            #model_name="error name"
            model_extr=doc_extract.getModel()
            model_name=model_extr.getName()
            self.name=model_name
        else:
            self.name="no name"
        
    def change_network(self, **options):
        """functions makes changes in the reaction network. uses following arguments:
        species = (set_x) -> reduces the list to a subset of reactions,to exclude not supportable reactions
        change_inflow = True: -> changes reactions, which happen in every compartment, to be distributable.
        add_reaction = xxxx -> adds a reaction """
        if "add_reactions" in options:
            pass
        self.reaction_network, self.species= generate_Network_to_be_analyzed(self.reaction_network, **options)
#     def reduce_network(self, species=None ,change_inflow=False):
#         self.reaction_network, self.species= generate_Network_to_be_analyzed(self.reaction_network, setspec=species, change_inflow=change_inflow)

            
    def SARs(self, **options):
        solutiondict={}
        solutiondict =setup_LP_Var(self.reaction_network,"SAR", **options)
        try:
            if solutiondict["error"]=="timeout_ERC":
                self.abort_run=True
                return
        except:
            pass
        self.SAR_solution=solutiondict["output_list"]
        self.ERC_dict=solutiondict["ERC_dict"]
        if "error" in solutiondict:
            print("error:"+str(solutiondict["error"]))
            self.SAR_solution_wrong=str(solutiondict["error"])
            #solutiondict['termination']
            if str(solutiondict["error"])=="Solution_timeout":
                self.info_so_far=[int(solutiondict["number_SAR_O"]),int(solutiondict["number_SAR_DO_only"])]
        if "constraints" in solutiondict:
            self.number_constr=solutiondict["constraints"]
        if "time_ERC" in solutiondict:
            self.time_ERC= solutiondict["time_ERC"]
            self.time_solve= solutiondict["time_solve"]

    def DOs(self, see_constraints=False):
        solutiondict =setup_LP_Var(self.reaction_network,"DO",see_constraints=see_constraints)
        self.DO_solution=solutiondict["output_list"]
        self.ERC_dict=solutiondict["ERC_dict"]
        
    def largest_DO(self, see_constraints=False):
        solutiondict =setup_LP_Var(self.reaction_network,"DO",see_constraints=see_constraints)
        self.largest_DO=solutiondict["output_list"]
        self.ERC_dict=solutiondict["ERC_dict"]
    
#     def draw_hasse(self, solution_set, output_style="pdf"):
#  DOs=False, shortform=False, show_species=False, 
#  show_new=False,show_compartments="length", second_value=False, analyze_only=False
    def draw_hasse(self, solution="SARs",use_naive=False, **options):
        
        draw_argument=""
        #try:
        if solution=="SARs":
            #draw_argument=getattr(self, "SAR_solution")
            draw_argument=self.SAR_solution
        elif solution=="DOs":
            options["DOs"]=True
            draw_argument=self.DO_solution
            #draw_argument=getattr(self, "DO_solution")
        else:
            print("incorrect job input")
        #create_Hasse3(draw_argument, self.reaction_network, **options)
        if use_naive:
            create_Hasse(draw_argument, self.reaction_network,**options)
        else:
            create_Hasse3(draw_argument, self.reaction_network, **options)
        #except AttributeError:
        if False:  
            print("Attribute not found")
            print(solution)
            print(".SARs exists: "+ str(hasattr(self, "SAR_solution")))
            print(".DOs exists: "+ str(hasattr(self,"DO_solution")))
#         create_Hasse3(list_of_solutions,LOR, DOs=False, shortform=False, show_species=False,
#         show_new=False,show_compartments="length", second_value=False, analyze_only=False):
   
#     def help_hasse(self):
#         print("use __class__.draw_hasse() to draw")
    def evaluate_solution(self, solution_set="default", output_style="exel"):
        """asd"""
        x=0
        try:
            bla=getattr(self,solution_set)
            x=analyze_DO(bla, self.reaction_network)
        except KeyError:
            print()
        return(x)
            
    def print_RN(self):
        line_List=[["reaction","reactants","products"]]
#         line_List=[]
        padding1=len("reaction")
        padding2=len("reactants")
        padding3=len("products")
        for reaction in self.reaction_network:
            left_alignment=""
            right_alignment=""
            for ele in range(len(reaction.listOfReactants)):
                left_alignment+=str(reaction.reac_stoich[ele]) + " "+ reaction.listOfReactants[ele] + "  "
                if ele<len(reaction.listOfReactants)-1:
                    left_alignment+="+ "
            for ele in range(len(reaction.listOfProducts)):
                right_alignment+=str(reaction.prod_stoich[ele]) + " "+ reaction.listOfProducts[ele] + "  "
                if ele<len(reaction.listOfProducts)-1:
                    right_alignment+="+ "
                    
            padding1 = max(padding1, len(reaction.defined_name))
            padding2 = max(padding2, len(left_alignment))
            padding3 = max(padding3, len(right_alignment))
            line_List.append([reaction.defined_name, left_alignment, right_alignment])
        for a,b,c in line_List:
            print(f'{a:<{padding1}}:  {b:<{padding2}}->  {c:<{padding3}}')
            
    def print_ERC(self):
        try:
            self.ERCs= self.ERCs
        except AttributeError:
            self.ERCs= create_closures(self.reaction_network, Timer_in_sec=120)
        if type(self.ERCs)==str:
                print("timeout ERC")
                return
        print("ERCs:", end=" ")
        [print(key,':',[ele.defined_name for ele in self.ERCs[key].reactions ]) for key in self.ERCs.keys()]
        
        
    def print_ERC_len(self):
        try:
            self.ERCs= self.ERCs
        except AttributeError:
            self.ERCs= create_closures(self.reaction_network, Timer_in_sec=120)
            if type(self.ERCs)==str:
                print("timeout ERC")
                return
#         print("ERCs:")
#         print([len(self.ERCs[key].reactions) for key in self.ERCs.keys()])
        output_array=[len(self.ERCs[key].reactions)-1 for key in self.ERCs.keys()]
#         ERC_counter=0
#         for key in self.ERCs.keys():
#             ERC_counter+=len(self.ERCs[key].reactions) 
#         avg_ERC=ERC_counter/len(self.ERCs)
#         return(avg_ERC)
        return(output_array)
    def get_comp_of_solution_x(self, index, solution="SARs", extra_species=""):
        get_min_compartments(self.reaction_network, getattr(self, solution), **extra_species)
#     def __repr__(self):

#         print(self.__dict__.keys())
#         return()
    #use generate_Network_to_be_analyzed with species subset or just normally to get species set
