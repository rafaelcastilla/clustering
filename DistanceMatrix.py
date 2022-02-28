'''
@author: olao, edited by Melitza Hesseling
'''
import copy


class DistanceMatrix(object):
    '''
    Class to store a symmetric distance matrix
    Requires a list with the name of the species
    '''

    def __init__(self, name_of_species):
        '''
        Constructor
        '''
        self.name_of_species = name_of_species
        # create the list of lists
        self.m = []
        # create for each row the column
        for r in range(len(name_of_species)):
            # create a column of lenght = # of species
            c = [0]*len(name_of_species)
            # append the colum to the list of lists
            self.m.append(c)
    


    def add_value_i_j(self, i, j, d):
        '''
        Add a value in position i,j (same for j,i)
        '''
        # raise expeption if i or j is out of range
        if (i > len(self.name_of_species)) or (j > len(self.name_of_species)):
            raise ValueError("Out of range")

        # setting the distance for the elements i,j and j,i
        self.m[i][j] = d
        self.m[j][i] = d
        
    

    def get_value_i_j(self,i,j):
        '''
        Get the value at position i,j
        '''
        return self.m[i][j]
    
    

    def change_name_species_i(self,i,new_name):
        '''
        Set a new name at species i
        '''
        self.name_of_species[i] = new_name

    

    def n_rows(self):
        '''
        Get the number of rows
        '''
        return len(self.m)
    
    

    def remove_species_i(self,i):
        '''
        Remove the row and column at position i
        '''
        # take the i'th species out of the specieslist
        self.name_of_species.pop(i)
        # take the species row out of the matrix
        self.m.pop(i)
        # iterate over all colums
        for j in range(len(self.m)):
            # take the i'th element out
            self.m[j].pop(i)
                
    

    def get_name_of_species(self):
        '''
         the name of the species
        '''
        return self.name_of_species



    def copy(self):
        '''
        Generate a copy of this distance matrix
        '''
        # creating a matrix with the amount of elements we need
        d_cop = DistanceMatrix(copy.deepcopy(self.name_of_species))
        # adding the values from the original matrix to our copy
        for i in range(len(self.name_of_species)-1):
            for j in range(i+1, len(self.name_of_species)):
                d_cop.add_value_i_j(i, j, self.get_value_i_j(i,j))
        # returning the copy of our original matrix
        return d_cop