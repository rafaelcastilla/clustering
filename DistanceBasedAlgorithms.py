'''
@author: olao, edited by Melitza Hesseling
'''
from DistanceMatrix import DistanceMatrix
from Sequence import Sequence
from Evolution import Evolution
from ToolsToWorkWithSequences import ToolsToWorkWithSequences


class DistanceBasedAlgorithms(object):
    '''
    A class used to construct a rooted evolution tree from a distance matrix
    '''

    def __init__(self, d_matrix):
        '''
        Constructor
        Takes a DistanceMatrix object as argument
        '''
        self.d_matrix = d_matrix
        
    
    def UPGMA(self):
        '''
        Uses UGMA to compute the evolution tree
        '''
        # creating a copy of the matrix, so that we can still work with the original matrix later on
        new_matrix = self.d_matrix.copy()
        # iterating until the matrix has schrunk into a 1x1 matrix
        while new_matrix.n_rows() > 1:
            # initializing the species indices and the minimal distance
            sp1 = 0
            sp2 = 1
            D_min = new_matrix.get_value_i_j(sp1, sp2)
            # looping over all elements in upper triangle (because of symmetry) of the matrix to find the smallest element
            for i in range(new_matrix.n_rows()-1):
                for j in range(i+1, new_matrix.n_rows()):
                    # if the distance found in the element is smaller than the current minimal distance
                    if D_min > new_matrix.get_value_i_j(i, j):
                        # set new minimal distance and indices of the species
                        D_min = new_matrix.get_value_i_j(i, j)
                        sp1 = i
                        sp2 = j

            # creating the new name for the node: (Sp(i):Dmin/2, Sp(j):Dmin/2) (Newick notation)
            new_name = "("+str(new_matrix.get_name_of_species()[sp1])+":"+str(D_min/2)+", "+str(new_matrix.get_name_of_species()[sp2])+":"+str(D_min/2)+")"

            # iterating over column/row of of matrix to change distances for sp1 to new calculated value
            for i in range(new_matrix.n_rows()):
                if i != sp1:
                    new_value = (new_matrix.get_value_i_j(i, sp1)+new_matrix.get_value_i_j(i ,sp2))/2
                    new_matrix.add_value_i_j(i, sp1, new_value)

            # change name of sp1 one to new name:
            new_matrix.change_name_species_i(sp1, new_name)

            # remove the column and row of sp2
            new_matrix.remove_species_i(sp2)

        # return the name (Newick notation) of the species of the last remaining element in the matrix
        return new_matrix.get_name_of_species()[0]



    '''
    Implement the NJ approach
    '''
    def NJ(self):
        '''
        Uses Neighbour-Joining to compute the evolution tree
        '''
        # creating a copy of the matrix, so that we can still work with the original matrix later on
        new_matrix = self.d_matrix.copy()
        # iterate until we have a 2x2 matrix
        # This is because the NJ algorithm does not return a rooted tree, that we do need for the Newick notation
        # So after this loop we'll connect the last two elements with picking one element as a root
        # Also, the NJ algorithm needs at least 3 elements to perform the limb lenght
        while new_matrix.n_rows() != 2:
            # creating an empty list with the total distance
            sum_by_row = [0]*new_matrix.n_rows()
            # looping over all elements in upper triangle (because of symmetry)
            for i in range(new_matrix.n_rows()-1):
                for j in range(i+1,new_matrix.n_rows()):
                    # getting and adding the elements value to the total distance
                    d = new_matrix.get_value_i_j(i,j)
                    sum_by_row[i] += d
                    sum_by_row[j] += d
            
            # creating the NJ matrix
            Q = DistanceMatrix(new_matrix.get_name_of_species())
            # looping over all elements in upper triangle (because of symmetry)
            for i in range(Q.n_rows()-1):
                for j in range(i+1,Q.n_rows()):
                    # compute NJ element value
                    q_i_j = (Q.n_rows()-2)*new_matrix.get_value_i_j(i,j) - sum_by_row[i] - sum_by_row[j]
                    # add this value to the matrix
                    Q.add_value_i_j(i, j, q_i_j)
                    
            #check lowest value of Q
            # initialize the element and minimal distance
            lowest_i = 0
            lowest_j = 1
            min_dist = Q.get_value_i_j(lowest_i, lowest_j)
            # looping over all elements in upper triangle (because of symmetry)
            for i in range(Q.n_rows()-1):
                for j in range(i+1,Q.n_rows()):
                    # if the distance found in the element is smaller than the current minimal distance
                    if(min_dist > Q.get_value_i_j(i, j)):
                        # set new minimal distance and indices of the species
                        min_dist =  Q.get_value_i_j(i, j)
                        lowest_i = i
                        lowest_j = j
            
            # get value of lowest element in Q from original matrix
            d_lowest_i_lowest_j = new_matrix.get_value_i_j(lowest_i,lowest_j)

            # calculate limblength(i) and limblenght(j)
            d_lowest_i_u = 0.5*d_lowest_i_lowest_j + 0.5*(sum_by_row[lowest_i]-sum_by_row[lowest_j])/(Q.n_rows()-2)
            d_lowest_j_u = d_lowest_i_lowest_j - d_lowest_i_u

            # calculate average distance from this done to each other node and change distances for lowest_i to new calculated value
            for k in range(new_matrix.n_rows()):
                avg = (new_matrix.get_value_i_j(lowest_i,k) + new_matrix.get_value_i_j(lowest_j,k) - d_lowest_i_lowest_j)/2.0
                new_matrix.add_value_i_j(lowest_i, k, avg)
            
            # make sure that the distince between node and itself is 0
            new_matrix.add_value_i_j(lowest_i,lowest_i,0)            
            # change name of lowest_i one to new name            
            new_matrix.change_name_species_i(lowest_i, ("(" + new_matrix.get_name_of_species()[lowest_i] + ":" + str(d_lowest_i_u) + ", " + new_matrix.get_name_of_species()[lowest_j] + ":" + str(d_lowest_j_u) + ")" ))
            # remove the column and row of lowest_j
            new_matrix.remove_species_i(lowest_j)          

        # Last iteration. Newick needs a rooted tree, even if nj is not rooted!
        # create a root between the last two elements. Set the length between the last element and the root as the final distance between the last element and all the others.
        # for completeness, set the distance of the node comprissing all the other elements to the final node to 0
        # getting the distance between last two nodes
        d_lowest_i_u = new_matrix.get_value_i_j(0,1)
        # change the name of the first species to final Newic notation 
        new_matrix.change_name_species_i(0, ("(" + new_matrix.get_name_of_species()[0] + ":" + str(0.0) + ", " + new_matrix.get_name_of_species()[1] + ":" + str(d_lowest_i_u) + ")" ))
        # remove the second species
        new_matrix.remove_species_i(1)  
        # return the name (Newick notation) of the species of the last remaining element in the matrix  
        return new_matrix.get_name_of_species()[0]




def main():
    
    # Example from midterm
    species = ["A","B","(CD)","E","F"]
    my_distance_matrix = DistanceMatrix(species)
    dist = [[0,6, 29, 24, 30],[6, 0,  31, 26, 28],[29, 31, 0,  32, 15],[24, 26, 32, 0,  30],[30, 28, 15, 30,  0]]
    
    for i in range(len(dist)-1):
        for j in range(i+1,len(dist)):
            d = dist[i][j]
            my_distance_matrix.add_value_i_j(i, j, d)
    
    dba = DistanceBasedAlgorithms(my_distance_matrix)    
    print(dba.UPGMA())


    # Example from wikipedia
    species = ["A","B","C","D","E"]
    my_distance_matrix = DistanceMatrix(species)
    dist = [[0,5,9,9,8],[5,0,10,10,9],[9,10,0,8,7],[9,10,8,0,3],[8,9,7,3,0]]
    
    for i in range(len(dist)-1):
        for j in range(i+1,len(dist)):
            d = dist[i][j]
            my_distance_matrix.add_value_i_j(i, j, d)
    
    dba = DistanceBasedAlgorithms(my_distance_matrix)    
    print(dba.UPGMA())
    print(dba.NJ())


if __name__ == "__main__":
    main ()