'''
@authors: 
Janina Estañol Ruscalleda
Melitza Roxanne Hesseling
Francisco Rafael Castilla Patino
'''
from Sequence import Sequence
from Evolution import Evolution
from ToolsToWorkWithSequences import ToolsToWorkWithSequences
from DistanceMatrix import DistanceMatrix
from DistanceBasedAlgorithms import DistanceBasedAlgorithms


def main():

    # Initialise ancestral sequence
    sequence_ancestral = Sequence("SpeciesA", "ACTGACTGACTGACTGACTGACTGACTGACTGACTG")
    # Create transition matrix
    # transition_probability = {"A":{"G":0.02,"C":0.02,"T":0.02,"A":0.939,"_":0.0,"I":0.001}, 
    #                           "C":{"G":0.02,"C":0.939,"T":0.02,"A":0.02,"_":0.0,"I":0.001}, 
    #                           "G":{"G":0.939,"C":0.02,"T":0.02,"A":0.02,"_":0.0,"I":0.001},
    #                           "T":{"G":0.02,"C":0.02,"T":0.939,"A":0.02,"_":0.0,"I":0.001},
    #                           "_":{"G":0.002,"C":0.002,"T":0.002,"A":0.002,"_":0.992},
    #                           "I":{"G":0.25,"C":0.25,"T":0.25,"A":0.25}}
    # transition_probability = {"A":{"G":0.02,"C":0.02,"T":0.02,"A":0.939,"_":0.001,"I":0.0}, 
    #                           "C":{"G":0.02,"C":0.939,"T":0.02,"A":0.02,"_":0.001,"I":0.0}, 
    #                           "G":{"G":0.939,"C":0.02,"T":0.02,"A":0.02,"_":0.001,"I":0.0},
    #                           "T":{"G":0.02,"C":0.02,"T":0.939,"A":0.02,"_":0.001,"I":0.0},
    #                           "_":{"G":0.002,"C":0.002,"T":0.002,"A":0.002,"_":0.992},
    #                           "I":{"G":0.25,"C":0.25,"T":0.25,"A":0.25}}
    ins_rate = 0.0005  #beta
    del_rate = 0.0005  #alpha
    transition_probability = {"A":{"G":0.0,"C":0.0,"T":0.0,"A":(1-ins_rate-del_rate),"_":del_rate,"I":ins_rate}, 
                              "C":{"G":0.0,"C":(1-ins_rate-del_rate),"T":0.0,"A":0.0,"_":del_rate,"I":ins_rate}, 
                              "G":{"G":(1-ins_rate-del_rate),"C":0.0,"T":0.0,"A":0.0,"_":del_rate,"I":ins_rate},
                              "T":{"G":0.0,"C":0.0,"T":(1-ins_rate-del_rate),"A":0.0,"_":del_rate,"I":ins_rate},
                              "_":{"G":0.0,"C":0.0,"T":0.0,"A":0.0,"_":1,"I":0.0},
                              "I":{"G":0.25,"C":0.25,"T":0.25,"A":0.25}}
    # Create evolution class
    evolution = Evolution(sequence_ancestral, transition_probability)
    
    # Split ancestral species into SpeciesA and SpeciesC
    evolution.split_species_in_two("SpeciesA", "SpeciesC")
    # # Evolve all species for 100 generations
    evolution.evolve(100)
    # # Split SpeciesA into SpeciesA and SpeciesB
    evolution.split_species_in_two("SpeciesA", "SpeciesB")
    # # Evolve all species for 50 generations
    evolution.evolve(100)
    # # Split SpeciesC into SpeciesC and SpeciesD
    evolution.split_species_in_two("SpeciesC", "SpeciesD")
    # # Evolve all species for 150 generations
    evolution.evolve(100)

    # Get the list of all species and sort it
    specieslist = evolution.get_list_of_species_name()
    specieslist.sort()

    # Initialising ToolsToWorkWithSequences class
    stats = ToolsToWorkWithSequences()

    # Print the sequence and nucleotide frequency for each species in the list
    for species in specieslist:
        print(evolution.get_sequence_species(species))


    # Creating distance matrix from four sequences
    my_distance_matrix = DistanceMatrix(specieslist)
    # Iterating over all combinations of species to compute distances
    for i in range(len(specieslist)-1):
        for j in range(i+1, len(specieslist)):
            # computing distance between the two sequences
            dist = stats.observed_pairwise_nucleotide_distance(evolution.get_sequence_species(specieslist[i]), evolution.get_sequence_species(specieslist[j]))
            # adding distance to the matrix
            my_distance_matrix.add_value_i_j(i, j, dist)


    dba = DistanceBasedAlgorithms(my_distance_matrix)    
    print("UPGMA:", dba.UPGMA())
    print("NJ:", dba.NJ())


if __name__ == "__main__":
    main ()