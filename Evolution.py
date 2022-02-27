'''

@author: olao
'''
from Sequence import Sequence
from Alias_Vose import RandomMultinomial

class Evolution(object):
    '''
    classdocs
    '''
    
    def __init__(self, ancestral_sequence, transition_matrix):
        '''
        Constructor
        '''
        if not isinstance(ancestral_sequence, Sequence):
            raise TypeError("ancestral_sequence must be a Sequence object!")
        self.insertion={}
        self.insertion[ancestral_sequence.get_name()]=[0]*ancestral_sequence.sequence_length()
        self.evolving_sequences = {}
        self.evolving_sequences[ancestral_sequence.get_name()] = ancestral_sequence
        self.transition_matrix = transition_matrix
        self.random_transition = {}
        for key in self.transition_matrix:
            self.random_transition[key] = RandomMultinomial(list(transition_matrix[key].values()))
        
    
    def get_sequence_species(self,name):
        return self.evolving_sequences[name]
    
    
    def get_list_of_species_name(self):
        return list(self.evolving_sequences.keys())
    
    
    def split_species_in_two(self,name_of_species_that_splits, new_name_of_species):
        new_sequence = self.evolving_sequences[name_of_species_that_splits].copy(new_name_of_species)
        self.evolving_sequences[new_sequence.get_name()] = new_sequence
    
    
    '''
     you give this function an evolution time and it returns a mutated sequence
    
    '''
    def evolve(self, generations):
        for gen in range(generations):                                          #the number of generations
            for species in self.evolving_sequences:                             #the number of sequences
                sequence_species = self.evolving_sequences[species]             #name the species
                self.insertion[species]=[0]*sequence_species.sequence_length()  #inicializate the diccionary for key specie and list of of insertion
                for n in range(sequence_species.sequence_length()):             ##look at each nucleotide in the sequence
                    nucleotide_at_position_n = sequence_species.nucleotide_at_position(n)                   #look the nucleotid in position n
                    propose_change = self.random_transition[nucleotide_at_position_n].sample()              
                    nucleotide_propose_change = list(self.transition_matrix[nucleotide_at_position_n].keys())[propose_change]   #the nucleotid tha can change
                    print(nucleotide_at_position_n,nucleotide_propose_change)                   #print the chanche
                    if(nucleotide_propose_change=="I"):                         #fi have a insertion
                        self.insertion[species][n]=1                            #add 1 in a indice
                    else:                                                           
                         sequence_species.mutate_nucleotide_at_position(n,nucleotide_propose_change) #make the mutatio

            for species in self.evolving_sequences:  
                sequence_species = self.evolving_sequences[species] 
                leng=len(self.insertion[species])# the leen of specie
                #print(leng)
               #va mirado posicion por posicion
                for n in range(leng):
                    if (self.insertion[species][n]==1): #iff have a insertion
                        propose_change = self.random_transition["I"].sample()   
                        nucleotide_propose_change = list(self.transition_matrix["I"].keys())[propose_change] #the new nucleotid
                        sequence_species.mutate_nucleatide_insertion(n,nucleotide_propose_change)       #make the insertion
                        
                        for species_2 in self.evolving_sequences:                       #make the deletion of others species
                            if(species!=species_2):
                                sequence_species = self.evolving_sequences[species_2]       
                                sequence_species.add_deletion(n)
                                self.insertion[species_2].insert(n,0)

                    


                        

                        

                        

            
def main():
    sequ="ACTGACTGACTGACTGACTGACTGACTGACTGACTG"
    sequence_ancestral = Sequence("Ancestral", sequ)
    transition_probability = {"A":{"G":0.04,"C":0.04,"T":0.04,"A":0.81,"_":0.02,"I":0.01}, "C":{"G":0.04,"C":0.81,"T":0.04,"A":0.04,"_":0.02,"I":0.01}, "G":{"G":0.81,"C":0.04,"T":0.04,"A":0.04,"_":0.02,"I":0.01},"T":{"G":0.04,"C":0.04,"T":0.81,"A":0.04,"_":0.02,"I":0.01},"_":{"G":0.04,"C":0.04,"T":0.04,"A":0.04,"_":0.81,"I":0.01},"I":{"G":0.25,"C":0.25,"T":0.25,"A":0.25}}
    evolution = Evolution(sequence_ancestral, transition_probability)
    
    evolution.split_species_in_two("Ancestral", "Species2")
    
    evolution.evolve(10)
    
    print(evolution.get_sequence_species("Ancestral"))
    print(evolution.get_sequence_species("Species2"))  
    print(len(sequ))
    print(evolution.get_sequence_species("Ancestral").sequence_length())  

    
if __name__ == "__main__":
    main ()         
        
        