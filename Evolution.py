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
    
    
    def evolve(self, generations):
        for gen in range(generations):
            for species in self.evolving_sequences:
                sequence_species = self.evolving_sequences[species]
                for n in range(sequence_species.sequence_length()):
                    nucleotide_at_position_n = sequence_species.nucleotide_at_position(n)
                    propose_change = self.random_transition[nucleotide_at_position_n].sample()
                    nucleotide_propose_change = list(self.transition_matrix[nucleotide_at_position_n].keys())[propose_change]
                    print(nucleotide_at_position_n,nucleotide_propose_change)
                    sequence_species.mutate_nucleotide_at_position(n,nucleotide_propose_change)
                    if(nucleotide_propose_change=="I"):
                        if(n!=sequence_species.sequence_length()-1):
                            n+=1               

            
def main():
    sequence_ancestral = Sequence("Ancestral", "ACTGACTGACTGACTGACTGACTGACTGACTGACTG")
    transition_probability = {"A":{"G":0.029,"C":0.029,"T":0.029,"A":0.88,"D":0.002,"I":0.002}, "C":{"G":0.029,"C":0.88,"T":0.029,"A":0.029,"D":0.002,"I":0.002}, "G":{"G":0.88,"C":0.029,"T":0.029,"A":0.029,"D":0.002,"I":0.002},"T":{"G":0.029,"C":0.029,"T":0.88,"A":0.029,"D":0.002,"I":0.002},"D":{"G":0.029,"C":0.029,"T":0.029,"A":0.029,"D":0.88,"I":0.002},"I":{"G":0.029,"C":0.029,"T":0.029,"A":0.029,"D":0.002,"I":0.88}}
    evolution = Evolution(sequence_ancestral, transition_probability)
    
    evolution.split_species_in_two("Ancestral", "Species2")
    
    evolution.evolve(1000)
    
    print(evolution.get_sequence_species("Ancestral"))
    print(evolution.get_sequence_species("Species2"))    

    
if __name__ == "__main__":
    main ()         
        
        