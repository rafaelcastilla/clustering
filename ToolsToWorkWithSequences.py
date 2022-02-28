'''
@author: Melitza Hesseling
'''
from Sequence import Sequence


class ToolsToWorkWithSequences(object):
    '''
    Some methods to analyse sequences
    '''

    def __init__(self):
        '''
        Constructor
        '''
    
    def nucleotide_statistics(self, sequence):
        '''
        Method to calculate the percentage of the A, C, T, G found in the sequence
        '''
        # if the input for sequence is not of the class 'Sequence', raise an error
        if not isinstance(sequence, Sequence):
            raise TypeError("Chosen sequence must be a Sequence object!")
        
        # Initializing dictionary for all nucleotieds
        dict_stats = {}
        # Iterate over all positions in the sequence
        for i in range(sequence.sequence_length()):
            # Nucleotide at position i
            nucl = sequence.nucleotide_at_position(i)
            # If the nucleotide is already present in the dictionary, increase its value by 1/sequencelength
            if nucl in dict_stats:
                dict_stats[nucl] += 1.0/sequence.sequence_length()
            # Otherwise add the nucleotide to the dictionary and set its value to 1/sequencelength
            else:
                dict_stats[nucl] = 1.0/sequence.sequence_length()
        # Return a dictionary with the frequencies of the nucleotides
        return dict_stats

    
    def observed_pairwise_nucleotide_distance(self, S1, S2):
        '''
        Takes two sequences and returns the number of nucleotides that for the same 
        position are different in the two sequences
        '''
        # if the input for sequence is not of the class 'Sequence', raise an error
        if not isinstance(S1, Sequence) or not isinstance(S2, Sequence):
            raise TypeError("Chosen sequences must be a Sequence objects!")
        
        # initialising counter
        counter = 0
        # iterating over all the nucleotide positions in the sequences
        for i in range(S1.sequence_length()):
            # if the nucleotide at that position is not the same between the sequences, add one to the counter
            if S1.nucleotide_at_position(i) != S2.nucleotide_at_position(i):
                counter +=1
        # return the number of differences between the sequences
        return counter