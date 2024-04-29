"""
Attachment A.11

Algotithm 4 - algorithm for finding a APN function G:F_2^(n-1) -> F_2^(n-1) which represents the APN function from the trim spectrum of a give function F (if such APN trim exists) 

REFERENCES: Chapter 4 of the thesis
            [1] - Christof Beirle, Gregor Leander and Leo Perrin - Trims and extensions of quadratic APN functions
            [8] - LPP-crypto/sboxu: Tools for studying S-boxes (https://github.com/lpp-crypto/sboxU)
"""

from sboxU import *
import numpy as np

# value which cannot be obtained in some following functions
init_value = -1


# function which computes the Inner_product of two numbers
# INPUT:  - integer "a"
#         - integer "b"
# OUTPUT: - integer "result"
def Inner_Product(a,b):
    # l&x is "l" AND "x", which sets bit to 1 if and only if both bits are 1
    # .bits() convert the integer into bits
    # .count(1) return number of "1"
    number_of_bits = Integer(a & b).bits().count(1)
    # result is from \F_2
    result = number_of_bits%2
    return result


# function which computes the trim spectrum for a given function F
# the function follows the theory in the thesis which is based on theory provided in [1]
# INPUT:  - list "F" which is a look-up table of a function F
#         - integer "n", which is dimension of the function F
# OUTPUT: - list "trims" which contains all functions from the trim spectrum of F
#         - list "parameters" which contains parameters (alpha, beta, gamma, epsilon) of a trim  
#              - trim of a function on the "i"-th position (function trims[i]) corresponds to the parameters on the "i"-th position (parameters[i])
def Finding_Trim_Spectrum(F,n):
    epsilon=0
    trim_spectrum=[]
    perameters = []
    for H in range(2):
        for alpha in range(1,2**n):
            if (H == 1):
                epsilon = 1
                while (Inner_Product(epsilon,alpha)==0):
                    epsilon = epsilon+1
            for beta in range(1,2**n):
                gamma = 1
                while (Inner_Product(beta,gamma)==0):
                    gamma=gamma+1
                trim = [0]*(2**n)
                for x in range(2**n):
                    if(Inner_Product(alpha,x)==0):
                        if(Inner_Product(gamma,F[x^epsilon])==0):
                            trim[x]=int(F[x^epsilon])
                        else:
                            trim[x]=int(F[x^epsilon]^beta)
                    if(Inner_Product(alpha,x)==1):
                        trim[x]=init_value
                trim_spectrum.append(trim)
                perameters.append([alpha, beta, gamma, epsilon])
    return trim_spectrum, perameters


# function which computes DDT for a trim of a function
# this function is specific for trim, because the DDT needs to be evaluated with the respect to the fact that trim is a function from "alpha^\perp" into "gamma^\perp",
# which means that the trim is not defined on all elements from F_2^n
# values for which the trim is not defined are marked with "init_value" in the look-up table of the trim
# INPUT: - list "lut" which is a look-up table for the trim of a function 
#        - integer "n", which is dimension of the orignal function F
# OUTPUT: two dimensional list "DDT"
def Compute_DDT_For_Trim(lut,n):
    ddt = np.zeros((2**n, 2**n), dtype=int)
    for x in range(2**n):
        for a in range(1,2**n):
            if((lut[x]!=init_value) & (lut[x^a]!=init_value)):
                y = lut[x] ^ lut[x^a]
                ddt[a][y] += 1
    return ddt


# function dicides if the input trim of a function is APN or not
# this function is specific for trim, because the DDT needs to be evaluated with the respect to the fact that trim is a function from "alpha^\perp" into "gamma^\perp",
# which means that the trim is not defined on all elements from F_2^n
# INPUT: - list "function" which is a look-up table for the trim of a function 
#        - integer "n", which is dimension of the orignal function F
# OUTPUT: True/False
def Is_Function_APN_For_Trim(function,n):
    ddt = Compute_DDT_For_Trim(function,n)
    for i in range(1, 2**n):
        for j in range(2**n):
            if ddt[i][j] > 2:
                return False
    return True


# function which creates list of binary notation of number 1,2,...,2**n-1
# INPUT:  integer n
# OUTPUT: list of binary notation of integers 1,2,...,2*n-1, all binary notations have the same length equals to "n"
def Create_Arrays_of_Binary_Notation_For_Inputs(n):
    binary_array = []
    for i in range (1,2**n):
        binary_string = bin(i)[2:]
        bit_array = [int(bit) for bit in binary_string]
        # set the length of all binary notations to the same length equal to "n"
        for _ in range(n-len(bit_array)):
            bit_array.insert(0,0)
        binary_array.append(bit_array)
    return(binary_array)


# function which creates list of all elements from the set "alpha^\perp"
# INPUT:  - integer alpha
#         - integer n
# OUTPUT: list of all elements from the set "alpha^\perp"
def Find_All_Elements_of_Hyperplane(alpha,dimension):
    output_list = []
    # from the definition of "alpha^\perp", elements in this set are those for which <alpha,element> is equal to zero
    for i in range(2**dimension):
        if(Inner_Product(alpha,i)==0):
            output_list.append(i)
    return(output_list)


# function dicides if the given sequence is linear independent (the sequence contains "dimension" number of elements)
# INPUT:  - list "sequence" of elements we want to test for linear independency
#         - integer "dimension"
# OUTPUT: True/False
def Is_Sequence_Linear_Independent(sequence, dimension):
    combinations_of_coefficients = Create_Arrays_of_Binary_Notation_For_Inputs(dimension)
    for combination in combinations_of_coefficients:
        s = 0
        for i in range(len(combination)):
            s = s^(combination[i]*sequence[i])
        if (s == 0):
            return(False)
    return(True)


# recursive function for finding bases of hyperplane given by "alpha^\perp", where "alpha^\perp" is set of dimension "dimension"
# INPUT:  - list "sequence" of elements we want to test for linear independency
#         - integer n
# OUTPUT: True/False
def Find_Basis_of_Hyperplane(elements, dimension):
    if(dimension == 1):
        return([elements[1]])
    else:
        basis_of_lower_dimension = Find_Basis_of_Hyperplane(elements, dimension-1)
        for element in elements:
            if (not(element in basis_of_lower_dimension)):
                candidate_for_basis = basis_of_lower_dimension + [element]
                if (Is_Sequence_Linear_Independent(candidate_for_basis,dimension) == True):
                    return(candidate_for_basis)
                

# function which finds the isomorphic mapping between F_2^dimension and hyperplane of dimension "dimension"
# the mapping is represented as a two dimensional list, where element on position "mapping[0][i]"" (element from F_2^dimension) maps onto element "mapping[1][i]" (element from hyperplane)
# INPUT:  - list "basis" of elements which are the basis of the hyperplane of dimension "dimension"
#         - integer "dimension"
# OUTPUT: list "mapping"
def Find_Mapping(basis,dimension):
    # list combinations_of_coefficients consists of all possible combinations of zeros and ones in dimension "dimension" except for combination [0,0, ... , 0]
    # using this combinations (after adding combination [0,0, ..., 0]) and elements from basis, we can generate all elements from the hyperplane
    combinations_of_coefficients = Create_Arrays_of_Binary_Notation_For_Inputs(dimension)
    # adding combination [0,0, ..., 0]
    combinations_of_coefficients.insert(0,[0]*dimension)
    mapping = [[0]*(2**dimension),[0]*(2**dimension)]
    for i in range(2**dimension):
        mapping[0][i]=i
        for j in range(dimension):
            mapping[1][i]=(mapping[1][i]) ^ (combinations_of_coefficients[i][j]*basis[j])
    return(mapping)


# function which finds look-up table for a function G:\F_2^(n-1) -> \F_2^(n-1)  which is based on the trim of a function F: \F_2^n -> \F_2^n
# INPUT:  - integer "alpha" which represents the hyperplane "alpha^\perp"
#         - integer "gamma" which represents the hyperplane "gamma^\perp"
#         - list "trim" which represents the look-up table for the trim of a function
#             - note that the list is of the length 2**(n+1) and contains "init_value" on position where the trim is not defined
#         - integer "m" which represents the dimension of the hyperplanes "alpha^\perp" and "gamma^\perp"
# OUTPUT: list "lut" which is the look-up of the desired function
def Find_Function_G(alpha,gamma,trim,m):
    n = m+1
    basis_of_alpha_hyperplane = Find_Basis_of_Hyperplane(Find_All_Elements_of_Hyperplane(alpha,n),m)
    basis_of_gamma_hyperplane = Find_Basis_of_Hyperplane(Find_All_Elements_of_Hyperplane(gamma,n),m)
    mapping_for_alpha = Find_Mapping(basis_of_alpha_hyperplane,m)
    mapping_for_gamma = Find_Mapping(basis_of_gamma_hyperplane,m)
    lut = [0]*(2**m)
    for i in range(2**m):
        # "i" (element from \F_2^(n-1)) -> (via mapping \varphi) -> "mapping_for_alpha[1][i]" (element from "alpha^\perp")
        # "mapping_for_alpha[1][i]" (element from "alpha^\perp") -> (via trim of a function) -> "trim[mapping_for_alpha[1][i]]" which is element in "gamma^\perp"
        output_of_trim=trim[mapping_for_alpha[1][i]]
        index = 0
        #while ((mapping_for_gamma[1][index]!=output_of_trim) & (index<len(mapping_for_gamma[1])-1)):
        while ((mapping_for_gamma[1][index]!=output_of_trim) & (index<2**m-1)):
            index = index + 1
        # "trim[mapping_for_alpha[1][i]]" which is element in "gamma^\perp" -> (via mapping \pi) -> "mapping_for_gamma[0][index]"" (which is equal to integer "index"),
        # where "index" is integer such that "mapping_for_gamma[1][index]" == "trim[mapping_for_alpha[1][i]]"
        lut[i]=index
    return(lut)



# function which saves all founded values from the algorithm into the .txt file
# INPUT:  - list "trims" which contains all APN trims in a form of look-up table which contains value "init_value" on 2**(n-1) positions
#         - list "parameters" that contains "alpha", "beta", "gamma", "epsilon" which are parameters for which we can get APN trim from the input function
#         - list "functions_with_APN_trims" that contains all input functions, which have APN function in their trim spectrum
#         - list "functions_without_APN_trims" that contains all input functions, which do not have any APN function in their trim spectrum
#         - list "functions_G" that contains functions G:F_2^(n-1) -> F_2^(n-1) which are obtained from the APN trims
#         - integer "dimension" which is a dimension of input functions to the algorithm
#         - string "output_name" which is a name of the output file where this function will save the results
# OUTPUT: file "output_name"
def Print_Founded_Functions_Into_File(trims,parameters,functions_with_APN_trims,functions_without_APN_trims,functions_G,dimension,output_name):
    with open(output_name,'w') as output_file:
        output_file.write('The algorithm got as an input ')
        output_file.write(str(len(functions_without_APN_trims)+len(functions_with_APN_trims)))
        output_file.write(' APN functions. These functions were from dimension ')
        output_file.write(str(dimension))
        output_file.write('.\n')
        output_file.write('The number of functions with APN functions in their trim spectrum is: ')
        output_file.write(str(len(functions_with_APN_trims)))
        output_file.write('\n')
        output_file.write('The number of functions without APN functions in their trim spectrum is: ')
        output_file.write(str(len(functions_without_APN_trims)))
        output_file.write('\n')
        output_file.write('The list of APN functions from trim spectrum:\n')
    for i in range(len(trims)):
        with open(output_name,'a') as output_file:
            output_file.write(str(trims[i]))
            output_file.write(',')
            output_file.write('\n')

    with open(output_name,'a') as output_file:
        output_file.write('\n')
        output_file.write('The input function, APN trim and parameters: ')
        output_file.write('\n')
    for i in range(len(trims)):
            with open(output_name,'a') as output_file:
                output_file.write(str(i+1))
                output_file.write('\n')
                output_file.write(str(functions_with_APN_trims[i]))
                output_file.write('\n')
                output_file.write(str(trims[i]))
                output_file.write('\n')
                output_file.write(str(parameters[i]))
                output_file.write('\n')

    with open(output_name,'a') as output_file:
        output_file.write('\n')
        output_file.write('The input functions without APN function in their trim spectrum: ')
        output_file.write('\n')
        output_file.write('The number of such functions: ')
        output_file.write(str(len(functions_without_APN_trims)))
        output_file.write('\n')
        output_file.write('The list of such functions: ')
        output_file.write('\n')
    for i in range(len(functions_without_APN_trims)):          
        with open(output_name,'a') as output_file:
            output_file.write(str(functions_without_APN_trims[i]))
            output_file.write(',')
            output_file.write('\n')

    with open(output_name,'a') as output_file:
        output_file.write('\n')
        output_file.write('The input functions with APN function in their trim spectrum: ')
        output_file.write('\n')
        output_file.write('The number of such functions: ')
        output_file.write(str(len(functions_with_APN_trims)))
        output_file.write('\n')
        output_file.write('The list of such functions: ')
        output_file.write('\n')
    for i in range(len(functions_with_APN_trims)):          
        with open(output_name,'a') as output_file:
            output_file.write(str(functions_with_APN_trims[i]))
            output_file.write(',')
            output_file.write('\n')

    with open(output_name,'a') as output_file:
        output_file.write('\n')
        output_file.write('The functions "G" from \F_2^(n-1) to \F_2^(n-1) obtained from the APN trims.')
        output_file.write('\n')
        output_file.write('The number of such functions: ')
        output_file.write(str(len(functions_G)))
        output_file.write('\n')
        output_file.write('The list of such functions: ')
        output_file.write('\n')
    for i in range(len(functions_G)):          
        with open(output_name,'a') as output_file:
            output_file.write(str(functions_G[i]))
            output_file.write(',')
            output_file.write('\n')

    with open(output_name,'a') as output_file:
        output_file.write('\n')
        output_file.write('The parameters for the APN trims of a function.')
        output_file.write('\n')
        output_file.write('The list of such parameters: ')
        output_file.write('\n')
    for i in range(len(parameters)):          
        with open(output_name,'a') as output_file:
            output_file.write(str(parameters[i]))
            output_file.write(',')
            output_file.write('\n')



# ------------------------------------------------------------------------------------------------------------------
# ----------------------------- ALGORITHM 4 - for finding all APN trims of a function ------------------------------
# ------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
# 1st PART - choosing input file with APN functions

"""
n = 7
load("sboxU/known_functions/sevenBitAPN.py")
functions = all_quadratics()





n = 6
load("FINAL_PROGRAMS/chapter4-results/chapter4-algorithm4-A11-output-n7.py")
functions = classified_functions_G()


n = 6
# representatives of classes 10 and 13 which does not appear in trim spectra of any 7-bit quadratic APN function
functions = [
[0, 17, 62, 28, 35, 19, 38, 37, 14, 19, 34, 12, 58, 6, 45, 34, 41, 10, 55, 39, 15, 13, 42, 27, 55, 24, 59, 39, 6, 8, 49, 12, 53, 2, 61, 57, 8, 30, 59, 30, 28, 39, 6, 14, 54, 44, 23, 62, 56, 61, 16, 38, 0, 36, 19, 4, 1, 8, 59, 1, 46, 6, 47, 52],
[9, 46, 15, 38, 19, 8, 43, 62, 26, 56, 58, 22, 6, 24, 24, 8, 3, 8, 3, 6, 56, 15, 6, 63, 63, 49, 25, 25, 2, 48, 26, 38, 20, 63, 9, 44, 37, 50, 6, 31, 18, 60, 41, 9, 37, 55, 32, 60, 17, 22, 10, 3, 1, 58, 36, 17, 56, 58, 5, 9, 46, 16, 45, 29]
]





n = 5
load("FINAL_PROGRAMS/chapter4-results/chapter4-algorithm4-A11-output-n6.py")
functions = classified_functions_G()


"""



output_name = 'chapter4-algorithm4-A11-output.txt'

#------------------------------------------------------------------------
# 2nd PART - setting up variables which are used throughout the algorithm

ct = 0
list_of_functions_that_have_APN_in_trims_spectrum = []
list_of_functions_that_do_not_have_APN_in_trims_spectrum = []
list_of_APN_trims = []
list_of_parameters = []
list_of_functions_G = []


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3rd PART - the algorithm itself; 
#          - in "i"-th run the algorithm analise "i"-th input function
#                  - it computes "trim_spectrum" for "i"-th function and analise if any APN function is in the "trim_spectrum"
#                  - if so, it saves the input function, trim, parameters and compute the function G:\F_2^(n-1) -> \F_2^(n-1)
#                  - if none of the functions in the trim spectrum is APN, it saves the input function into other list

for i in range(len(functions)):
    k=0
    trim_is_apn = 0
    trim_spectrum,parameters = Finding_Trim_Spectrum(functions[i],n)
    # while loop which analise all functions from the "trim_spectrum" until it finds APN function (if such function exists)
    while((k<len(trim_spectrum)) & (trim_is_apn == 0)):
        trim_is_apn = Is_Function_APN_For_Trim(trim_spectrum[k],n)
        if trim_is_apn==True:
            ct = ct+1
            list_of_functions_that_have_APN_in_trims_spectrum.append(functions[i])
            list_of_APN_trims.append(trim_spectrum[k])
            list_of_parameters.append(parameters[k])
            list_of_functions_G.append(Find_Function_G(parameters[k][0],parameters[k][2],trim_spectrum[k], n-1))
        k=k+1
    print(i+1, "vysledek, je tam APN funkce? ", trim_is_apn, "pocet APN funkci: ",ct)
    if(trim_is_apn==0):
        list_of_functions_that_do_not_have_APN_in_trims_spectrum.append(functions[i])

Print_Founded_Functions_Into_File(list_of_APN_trims,list_of_parameters,list_of_functions_that_have_APN_in_trims_spectrum,list_of_functions_that_do_not_have_APN_in_trims_spectrum,list_of_functions_G,n,output_name)














"""

# function which was used during programming for verifying that the basis can generate all elements from the given set "elements" which contains
# all elements from the hyperplane
# INPUT:  - list "basis" of elements which are the basis of the hyperplane of dimension "dimension"
#         - list "elements" which contains all elements from the hyperplane
#         - integer "dimension"
# OUTPUT: list "mapping"
def Can_Basis_Really_Generate_All_Elements(basis,elements,dimension):
    combinations_of_coefficients = Create_Arrays_of_Binary_Notation_For_Inputs(dimension)
    was_element_founded = [0]*len(elements)
    index = 0
    # since the list "combinations_of_coefficients" does not contain combination "[0, 0, ..., 0]", we need solve this case separately
    while ((elements[index]!=0) & (index<len(elements)-1)):
        index = index + 1
    if (elements[index] == 0):
        was_element_founded[index]=1
    for combination in combinations_of_coefficients:
        s = 0
        for i in range(len(combination)):
            s = s^(combination[i]*basis[i])
        index = 0
        while ((elements[index]!=s) & (index<len(elements)-1)):
            index = index + 1
        if (elements[index] == s):
            was_element_founded[index]=1
    if (sum(was_element_founded) == len(elements)):
        return (1)
    else:
        return (0)


"""