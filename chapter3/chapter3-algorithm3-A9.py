"""
--------------
Attachment A.9
--------------

Algorithm 3 - Search for linear-equivalence classes

REFERENCES: Chapter 3 of the thesis
            [2] - Christof Beirle and Gregor Leander - New instances of quadratic APN functions
            [3] - Christof Beirle, Marcus Brinkmann and Gregor Leander - Linearly self-equivalent APN permutations in small dimention

REQUIREMENTS:  - Python (version 3.12.2 was used), Python Software Foundation, https://www.python.org
               - Sage, the Sage MAthematical Software System (version 9.0 was used), Sage Developers, https://www.sagemath.org
			   - Cython Module for Python, Cython Developers, https://cython.readthedocs.io/en/latest/index.html
"""

# package for math over finite field
import galois
# package with procedure "product" which can generate all possible combinations of given elements
import itertools
# package sage.all is imported later in the algorithm

# defining finite field of 2 elements
GF = galois.GF(2**1)


# function for generating all possible companion matrices for given size
# INPUT:  dimension (size) of matrices that we want to generate
# OUTPUT: list of matricies in which each element is companion matrix (of some unknown order)
def Generate_Companion_Matrices(size):    
    # we save all possible combination of 0 and 1 (elements from F_2) for coefficients b_0,...,b_size
    all_coefficients = list(itertools.product([0, 1], repeat=size))
    # array for saving generated companion matrices
    companion_matrices = []
    # for every combination of zeros and ones in the array "coefficients" we will create "matrix" which can be a companion matrix
    for coeff_tuple in all_coefficients:
        # making polynomial from the coefficients in "coeff_tuple", we take the coefficients
        coef = list(coeff_tuple)
        # reverse the list because the methode "galois.Poly" makes the polynomial in a way that the element on 0-th position is the coefficient for x^len(coef)
        coef.reverse()
        # the element "coef[-1]" is coeeficient b_0 in the polynomial for which we are generating the companion matrix
        # but since we are interested only in matrices from GL, we need the b_0 to be equal to 1
        if (coef[-1] == 1):
            # at first we create zero matrix
            matrix = [[0]*size for _ in range(size)]
            # we set the last column of the matrix according to the "coeff_tuple"
            for i in range(size):
                matrix[i][-1] = coeff_tuple[i]
            # we add ones on the subdiagonal of the matrix
            for i in range(1, size):
                matrix[i][i-1] = 1
            companion_matrices.append(matrix)
    return companion_matrices


# function which verifies that the given sequence is not decreasing
# INPUT:  list of numbers
# OUTPUT: 0/1 
def Is_The_Sequence_Not_Decreasing(sequence):
    if (len(sequence)==1):
        return(1)
    else:
        for i in range(len(sequence)-1):
            if(sequence[i+1]<sequence[i]):
                return(0)
        return(1)


# function which for given array of polynomials "polynomials" verify if the condition of division of polynomials from the definition
# of RFC matrix is satisfied
# INPUT:  list of numbers
# OUTPUT: 0/1
def Division_Condition_Is_Satisfied(polynomials):
    zero_polynomial = galois.Poly([0], field=GF)
    length = len(polynomials)
    if length==1:
        return(1)
    else:
        for i in range(length-1):
            pol1 = polynomials[i]
            pol2 = polynomials[i+1]
            if (pol2 % pol1)!=zero_polynomial:
                return(0)
        return(1)


# we want to generate RCF matrices of prime order, therefore we need that each companion matrix is of the prime order
# thus if we search for prime order that exists in for every dimension of companion matrices in the RCF
# the dimensions are defined in the "comb", therefore we get dimensions from the "comb" and we find orders that appears in each "orders[comb[dimension]]"
# problematic is the identity matrix of dimension one which is only possible companion matrix (invertible) of dimension one, this matrix is of the order 1
# its the only companion matrix of order one, therefore in any "orders[comb[i]]" for i>1 the order "1" wont appear
# but if we multiply this matrix by itself "k" times, then we still get companion matrix, therefore we can think of this matrix as of any order
# therefore, if "comb" includes one, we put order one in the "output_orders" even though that the order "1" does not appear in any "orders[comb[i]]" for i>1
# INPUT:  - list "comb" which includes dimensions for the blocks in the RCF (blocks are companion matrices)
#         - list "orders" which includes all orders of companion matrices, on the i-th position are orders of the comapnion matrices of dimension "i"
# OUTPUT: - list of orders which apers in every 
def Select_Possible_Orders(orders,comb):
    output_orders = []
    # we split the funciton in two cases, first one is that the dimension of companion matrix is equal to the dimension of the RCF matrix, thus we are interested
    # only in the orders for the given dimension
    if len(comb) == 1:
        selected_orders = orders[comb[0]]        
        for order in selected_orders:
            is_element_in_output_orders = 0
            for element in output_orders:
                if element == order:
                    is_element_in_output_orders = 1
            if is_element_in_output_orders ==0:
                output_orders.append(order)
        return(output_orders)
    else:
        # first we deal with the ones in the "comb"
        position_in_comb = 0
        if comb[0] == 1:
            output_orders.append(1)
            condition = 1
            while (condition == 1):
                condition = 0
                position_in_comb = position_in_comb + 1
                if (position_in_comb < len(comb)):
                    if (comb[position_in_comb]==1):
                        condition = 1
        # if in the comb is any other number then one
        if (position_in_comb < len(comb)):
            # we select the first dimension in "comb" which is not one
            # we take all orders for this dimension and look up to the rest of dimension in the "comb" (if there are such other dimensions)
            selected_orders = orders[comb[position_in_comb]]
            for position_in_orders in range(len(selected_orders)):
                    selected_prime = selected_orders[position_in_orders]
                    prime_was_selected_before = 0 
                    for element in output_orders:
                        if element == selected_prime:
                            prime_was_selected_before = 1
                    if (prime_was_selected_before == 0):
                        the_prime_is_in_orders = [0]
                        for i in range(position_in_comb+1,len(comb)):
                            prime_is_there = 0
                            for element in orders[comb[i]]:
                                if element == selected_prime:
                                    prime_is_there = 1
                            the_prime_is_in_orders.append(prime_is_there)
                        if (sum(the_prime_is_in_orders) == (len(comb)-(position_in_comb+1))):
                            output_orders.append(selected_prime)
    return(output_orders)


# function which selects only those matrices that are of the given order            
# INPUT:  - number "order" which is the desired order of selected matrices
#         - list "matrices" which includes all matrices from which we need to select
#         - list "list_of_orders" which includes orders of matrices: order of matrix matrices[i] is equal to list_of_orders[i]
# OUTPUT: - list of matrices of given order 
def Select_Matrices_Based_On_Order(order,matrices,list_of_orders):
    output_matrices = []
    # the function can get as input order != 1 and in the list "matrices" only identity matrix of dimension 1. This matrix can be considered of any order for the algorithm
    # therefore we need to put it as an output
    if ((len(matrices[0]) == 1) & (order !=1)):
        output_matrices.append(matrices[0])
        return(output_matrices)
    for position in range(len(matrices)):
        if order==list_of_orders[position]:
            output_matrices.append(matrices[position])
    return(output_matrices)


# function which generates all RCF matrices of prime order of given dimension           
# INPUT:  - number "n" which is dimension of which we want the RCF matrices
# OUTPUT: - list of matrices of given order 
def Block_Diagonal_Matrices(n):
    output_matrices = []
    output_polynomials = []
    output_orders = []
    companion_matrices_of_given_size = []
    companion_matrices = []
    companion_matrices_orders = []
    # we append empty list because we want to have matrices of dimension one on the position [1]    
    companion_matrices.append([])
    companion_matrices_orders.append([])    
    #generate all possible combinations of dimensions of the companion matrices in the RFC
    combinations = Generate_Combinations(n)
    #for each size from 1 to n we generate all companion matrices of prime order
    for size in range(1,n+1):
        # for given dimension "size" we generate all companion matrices
        companion_matrices_of_given_size = Generate_Companion_Matrices(size)
        # from these matrices we select only those which have prime order
        (companion_matrices_of_prime_order, orders_of_matrices_for_given_size) = Select_Matrices_With_Prime_Order(size,companion_matrices_of_given_size)
        companion_matrices.append(companion_matrices_of_prime_order)
        companion_matrices_orders.append(orders_of_matrices_for_given_size)
    # for each combination "comb" we will try to make RCF matrix
    for comb in combinations:        
        # we select orders of matrices of dimension in "comb "
        distinct_orders = Select_Possible_Orders(companion_matrices_orders,comb)
        # for each order in the distinct order we will try to make RCF matrix
        for order in distinct_orders:
            companion_combinations=[]        
            # for each dimension we select matrices of order "order"
            for dimension in comb:
                matrices_for_dimension = []
                selected_matrices = Select_Matrices_Based_On_Order(order,companion_matrices[dimension],companion_matrices_orders[dimension])
                if selected_matrices != []:
                    for element in selected_matrices:
                        matrices_for_dimension.append(element)
                    companion_combinations.append(matrices_for_dimension)
            # if we found for each dimension some matrix
            if (len(companion_combinations)==len(comb)):
                # here we try all combination of companion matrices to make RCF matrix
                for companions in itertools.product(*companion_combinations):
                    # array for creating RCF matrix
                    RCF_matrix_in_block_form = []
                    # variable for indexing where to put each block (companion matrix) in the RCF matrix
                    place = 0
                    # array of polynomials which are taken from the last column of each companion matrix, it is needed for verifying the condition of divisibility of the polynomials from the definition of RCF
                    polynomials = []
                    # for each companion matrix
                    for companion in companions:
                        coef = []
                        size = len(companion)
                        # create rows in the RCF matrix with the block of companion matrix 
                        block = [[0]*place + row +[0]*(n-place-size) for row in companion]
                        # taking coefficients from the last column of the companion matrix
                        for i in range(size):
                            coef.insert(0,companion[i][-1])
                        # creating polynomial from the coefficients and append it into the array of all polynomials
                        coef.insert(0,1)
                        polynomials.append(galois.Poly(coef, field=GF))
                        # add the block into the output RCF matrix
                        RCF_matrix_in_block_form.extend(block)
                        # move the starting position for the next block
                        place = place + len(companion)
                    # verify that the matrix satisfy the condition of divisibility of polynomials from the definition of RFC, if so, append this matrix into the output array                    
                    if(Division_Condition_Is_Satisfied(polynomials)==1):
                        output_polynomials.append(polynomials)
                        output_orders.append(order)
                        output_matrices.append(RCF_matrix_in_block_form)              
    return(output_matrices,output_polynomials,output_orders)


# function for converting list of matrices (which are in the form of list) to the list of matrices as elements from GF
# INPUT:  - list "matrices" of matrices
# OUTPUT: - list of matrices as GF elements
def Convert_To_Finite_Field(matrices):
    finite_field_matrix = []
    for element in matrices:
        finite_field_matrix.append(GF(element))
    return finite_field_matrix


# function determing if given number is prime (or one) or not
# INPUT:  - number
# OUTPUT: - 0/1
def Is_Number_Prime_Or_One(number):
    if (number == 1):
        return(1)
    if number > 1:
        for i in range(2, int(number/2)+1):
            if (number % i) == 0:
                return(0)
                break
        else:
            return(1)
    else:
        return(0)


# function for creating identity matrix of given dimension as an element from GF
# INPUT:  - dimension
# OUTPUT: - matrix as an element from GF
def Create_Identity_Matrix(dimension):
    output_matrix = []
    for i in range(dimension):
        output_matrix.append([0]*dimension)
    for i in range(dimension):
        output_matrix[i][i]=1
    output_matrix=galois.GF(2**1)(output_matrix)
    return(output_matrix)


# function which generates number of invertibel matrices in given dimesion
# such number is order of group GL(dimension)(F_2)
# INPUT:  - dimension
# OUTPUT: - order of group GL
def Number_Of_Invertible_Matricecs_For_Given_Dimension(dimension):
    number = 1
    for i in range(dimension):
        number = number*(2**dimension-2**i)
    return(number)


# function which selects from the given list of matrices those, which are of prime order
# INPUT:  - dimension
# OUTPUT: - order of group GL
def Select_Matrices_With_Prime_Order(dimension,matrices):
    dimension = len(matrices[0])
    order_of_group_GL = Number_Of_Invertible_Matricecs_For_Given_Dimension(dimension)
    identity_matrix = Create_Identity_Matrix(dimension)
    output_matrices = []
    output_orders = []
    counter = 0
    for matrix in matrices:
        matrix_in_list = []
        matrix_in_list.append(matrix)
        matrix_in_GF = GF(matrix)
        multiplication_matrix = matrix_in_GF
        order = 1
        limit = order_of_group_GL/2+1
        while ((multiplication_matrix != identity_matrix).any()) & (order < limit):
            multiplication_matrix = multiplication_matrix@matrix_in_GF
            order = order +1
            if order == limit:
                print("ERROR")
        if (Is_Number_Prime_Or_One(order)==1):
            counter = counter +1
            output_matrices.append(matrix)
            output_orders.append(order)
    return(output_matrices, output_orders)


# function generates all possible combinations of integers such that:
# - the sum of integers is equal to "n"
# - the numbers are not in descending order
# INPUT:  - desired sum of elements in combination "n"
#         - number "max_number"  
# OUTPUT: - order of group GL
def Generate_Combinations(n, max_number=None):
    if max_number is None:
        max_number = n
    if n == 0:
        return [[]]
    if n < 0 or max_number == 0:
        return []
    combination = []
    for i in range(1, min(max_number, n) + 1):
        remaining_combinations = Generate_Combinations(n - i, i)
        for comb in remaining_combinations:
            result = [i]+comb
            result.reverse()
            result.sort()
            combination.append(result)
    return combination


# function which removes identity matrix (and coresponding order and polynomials) from the given lists 
# INPUT:  - list "matrices" which contains RCF matrices in the form of GF
#         - list "orders" which contains orders of the matrices in the list "matrices"
#         - list "polynomials" which contains polynomials which corespond to the companion matrices which forms the given matrices
# OUTPUT: - lists of matrices, orders and polynomials
def Remove_Identity_Matrix(matrices, orders, polynomials):
    output_matrices = []
    output_orders = []
    output_polynomials = []  
    for i in range(len(matrices)):
        if (orders[i]!=1):
            output_matrices.append(matrices[i])
            output_orders.append(orders[i])                                   
            output_polynomials.append(polynomials[i])
    return(output_matrices,output_orders,output_polynomials)


# function which prints found RCF matrices with their orders and coresponding polynomials to the file "chapter3-output.txt" 
# INPUT:  - list "matrices" which contains RCF matrices in the form of GF
#         - list "orders" which contains orders of the matrices in the list "matrices"
#         - list "polynomials" which contains polynomials which corespond to the companion matrices which forms the given matrices
def Print_Founded_RCF_Matrices_Of_Prime_Order(matrices, orders, polynomials,dimension):
    with open('chapter3-output.txt','w') as output_file:
        output_file.write('For dimension ')
        output_file.write(str(dimension))
        output_file.write(' we found ')
        output_file.write(str(len(matrices)))
        output_file.write(' RCF matrices. They are: \n')
    for i in range(len(matrices)):
        with open('chapter3-output.txt','a') as output_file:
            output_file.write('Order of the RCF matrix: ')
            output_file.write(str( orders[i]))
            output_file.write('\n')
            output_file.write('Polynomials of comapnions matrices: ')
            output_file.write(str( polynomials[i]))
            output_file.write('\n')
            for row in matrices[i]:            
                output_file.write(str(row))
                output_file.write('\n')
            output_file.write('\n')


# function which generates all possible tuples of RCF matrices such that the matrices in tuple have to be of the same order
# function also generates tuples of polynomials which corespond to the tuples of matrices
# INPUT:  - list "matrices" which contains RCF matrices in the form of GF
#         - list "orders" which contains orders of the matrices in the list "matrices"
#         - list "polynomials" which contains polynomials which corespond to the companion matrices which forms the given matrice
# OUTPUT: - lists of tuples of matrices, orders and list of tuples of polynomials
def Generate_All_Tuples_For_Each_Possible_Order (matrices, orders, polynomials):
    all_prime_orders_tuples=[]    
    exclusive_prime_orders = []
    polynomials_for_tuples = []
    help_list_for_polynomials = []
    for j in range(len(polynomials)):
        help_list_for_polynomials.append(j)
    # we generate tuples seperatly for each "order" in "orders"
    for order in orders:
        if order>1:
            if ((order in exclusive_prime_orders)==0):
                matrices_of_order = []
                help_polynomials_for_order = []
                polynomials_for_order = []
                help_polynomials_for_tuples = []
                exclusive_prime_orders.append(order)
                for position in range( len(orders)):
                    if order == orders[position]:
                        matrices_of_order.append(matrices[position])
                        help_polynomials_for_order.append(help_list_for_polynomials[position])
                all_prime_orders_tuples.append(Tuples(matrices_of_order,2))
                help_polynomials_for_tuples.append(Tuples(help_polynomials_for_order,2).list())
                for element in help_polynomials_for_tuples[0]:
                    first_polynomial = polynomials[element[0]]
                    second_polynomial = polynomials[element[1]]
                    polynomials_for_order.append([first_polynomial,second_polynomial])
                polynomials_for_tuples.append(polynomials_for_order)
    return(all_prime_orders_tuples, exclusive_prime_orders,polynomials_for_tuples)


# function which selects tuples of RCF matrices up to power-similarity
# INPUT:  - list "matrices" which contains RCF matrices in the form of GF
#         - list "orders" which contains orders of the matrices in the list "matrices"
#         - list "polynomials" which contains polynomials which corespond to the companion matrices which forms the given matrice
# OUTPUT: - lists of tuples of matrices, orders and list of tuples of polynomials up to power-similarity
def Select_Tuples_Of_RCF_Matrices_Up_To_Power_Similarity(all_prime_orders_tuples, exclusive_prime_orders, polynomials_for_tuples):
    not_power_similar_RCF = []
    coresponding_polynomials = []
    coresponding_orders = []
    for position in range(len(exclusive_prime_orders)):
        order = exclusive_prime_orders[position]
        all_tuples_of_given_order = all_prime_orders_tuples[position]
        polynomial_tuples_of_given_order = polynomials_for_tuples[position]
        # list for note if a tuple was power-equivalent to some previously examined tuples
        is_power_equivalent = [0]*len(all_tuples_of_given_order)
        for first_position in range(len(all_tuples_of_given_order)):
            if (is_power_equivalent[first_position]==0):
                for second_position in range(first_position+1,len(all_tuples_of_given_order)):
                    if (is_power_equivalent[second_position]==0):
                        i = 1
                        found = 0
                        while ((i<order) & (found == 0)):
                            if((all_tuples_of_given_order[first_position][0]**i).is_similar(all_tuples_of_given_order[second_position][0])&
                            (all_tuples_of_given_order[first_position][1]**i).is_similar(all_tuples_of_given_order[second_position][1])):
                                is_power_equivalent[second_position] = 1
                                found = 1
                            i = i+1
                not_power_similar_RCF.append([all_tuples_of_given_order[first_position][0],all_tuples_of_given_order[first_position][1]])
                coresponding_orders.append(order)
                coresponding_polynomials.append([polynomial_tuples_of_given_order[first_position][0], polynomial_tuples_of_given_order[first_position][1]])
    return(not_power_similar_RCF, coresponding_orders, coresponding_polynomials)


# function which prints found RCF matrices with their orders and coresponding polynomials to the file "chapter3-output.txt" 
# INPUT:  - list "matrices" which contains RCF matrices in the form of GF
#         - list "orders" which contains orders of the matrices in the list "matrices"
#         - list "polynomials" which contains polynomials which corespond to the companion matrices which forms the given matrices
def Print_Tuples(matrices, orders, polynomials, RCF_matrices, orders_of_RCF_matrices, polynomials_for_RCF_matrices, dimension):
    identity_matrix = Create_Identity_Matrix(dimension)
    with open('chapter3-output.txt','a') as output_file:
        output_file.write('-------------------------------- \n')
        output_file.write('We have found ')
        output_file.write(str(len(matrices)+len(RCF_matrices)+len(RCF_matrices)))
        output_file.write(' equivalence classes. \n')
        output_file.write('We have found ')
        output_file.write(str(len(matrices)))
        output_file.write(' tuples such that ord(A)=ord(B)=p for some prime p. \n')
        output_file.write('We have found ')
        output_file.write(str(len(RCF_matrices)))
        output_file.write(' tuples such that ord(A)=p for some prime p and B=I_')
        output_file.write(str(dimension))
        output_file.write('. \n')
        output_file.write('We have found ')
        output_file.write(str(len(RCF_matrices)))
        output_file.write(' tuples such that A=I_')
        output_file.write(str(dimension))
        output_file.write(' and ord(B)=p for some prime p. \n')
        output_file.write('They are: \n')
        output_file.write('\n')
    with open('chapter3-output.txt','a') as output_file:
        output_file.write('-------------------------------------------------------------------------------- \n')
        output_file.write('Equivalence classes where ord(A)=ord(B) is prime. \n')
        output_file.write('\n')
    for i in range(len(matrices)):
        with open('chapter3-output.txt','a') as output_file:
            output_file.write('Order of the tuple: ')
            output_file.write(str( orders[i]))
            output_file.write('\n')
            output_file.write('Polynomials of companion matrices: \n')
            output_file.write(str( polynomials[i][0]))
            output_file.write('\n')
            output_file.write(str( polynomials[i][1]))
            output_file.write('\n')
            for row in matrices[i][0]:            
                output_file.write(str(row))
                output_file.write('\n')
            output_file.write('\n')
            for row in matrices[i][1]:            
                output_file.write(str(row))
                output_file.write('\n')
            output_file.write('\n')
    with open('chapter3-output.txt','a') as output_file:
        output_file.write('-------------------------------------------------------------------------------- \n')
        output_file.write('Equivalence classes where ord(A) is prime and B is identity matrix. \n')
        output_file.write('\n')
    for i in range(len(RCF_matrices)):
        with open('chapter3-output.txt','a') as output_file:
            output_file.write('Order of the matrix A: ')
            output_file.write(str( orders_of_RCF_matrices[i]))
            output_file.write('\n')
            output_file.write('Polynomials of companion matrices in matrix A: \n')
            output_file.write(str( polynomials_for_RCF_matrices[i]))
            output_file.write('\n')
            for row in RCF_matrices[i]:            
                output_file.write(str(row))
                output_file.write('\n')
            output_file.write('\n')
            for row in identity_matrix:            
                output_file.write(str(row))
                output_file.write('\n')
            output_file.write('\n')
    with open('chapter3-output.txt','a') as output_file:
        output_file.write('-------------------------------------------------------------------------------- \n')
        output_file.write('Equivalence classes where A is identity matrix and ord(B) is prime. \n')
        output_file.write('\n')
    for i in range(len(RCF_matrices)):
        with open('chapter3-output.txt','a') as output_file:
            output_file.write('Order of the matrix B: ')
            output_file.write(str( orders_of_RCF_matrices[i]))
            output_file.write('\n')
            output_file.write('Polynomials of companion matrices in matrix B: \n')
            output_file.write(str( polynomials_for_RCF_matrices[i]))
            output_file.write('\n')
            for row in identity_matrix:            
                output_file.write(str(row))
                output_file.write('\n')
            output_file.write('\n')
            for row in RCF_matrices[i]:            
                output_file.write(str(row))
                output_file.write('\n')
            output_file.write('\n')



# --------------------------------------------------------------------------------------------------------
# -------------------- ALGORITHM 3 - algorithm for finding linear-equivalence classes --------------------
# --------------------------------------------------------------------------------------------------------

# THIS VARIABLE CAN BE CHANGE TO EXAMINE MATRICES OF DIMENSION "n"
n=7




# create all RCF matrices of prime order
RCF_matrices, polynomials_for_RCF, orders_of_RCF = Block_Diagonal_Matrices(n)
# convert them into GF elements
finite_field_matrices = Convert_To_Finite_Field(RCF_matrices)
# remove identity matrix
RCF_matrices, orders_of_RCF, polynomials_for_RCF = Remove_Identity_Matrix(RCF_matrices, orders_of_RCF, polynomials_for_RCF)
#print found matrices
Print_Founded_RCF_Matrices_Of_Prime_Order(RCF_matrices,orders_of_RCF,polynomials_for_RCF,n)

# import sage.all (we cannot import it earlier because during the programming of this algorithm the sage.all wasn't used in the begging, thus if we import it in the beganing,
# the algorithm won't work because there is would be a collison)
from sage.all import *

# transform matrices into the GF elements
RCF_matrices_GF = []
for matice in RCF_matrices:
    RCF_matrices_GF.append(matrix(GF(2),matice))


# we generate all possible tuples of matrices (RCF(A), RCF(B)) such that ord(RCF(A))=ord(RCF(B)) = prime 
all_prime_orders_tuples, exclusive_prime_orders, polynomials_for_tuples = Generate_All_Tuples_For_Each_Possible_Order(RCF_matrices_GF, orders_of_RCF, polynomials_for_RCF)
# we select tuples up to power-similarity
RCF_matrices_up_to_power_similarity, orders_of_RCF_matrices_up_to_power_similarity, corresponding_polynomials = Select_Tuples_Of_RCF_Matrices_Up_To_Power_Similarity(all_prime_orders_tuples, exclusive_prime_orders, polynomials_for_tuples)             
# we print tuples of RCF matrices up to power-similarity and similarity      
Print_Tuples(RCF_matrices_up_to_power_similarity, orders_of_RCF_matrices_up_to_power_similarity, corresponding_polynomials, RCF_matrices_GF, orders_of_RCF, polynomials_for_RCF, len(RCF_matrices[0]))

