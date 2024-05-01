"""
--------------
Attachment A.7
--------------

Algotithm 2 - modified - algorithm for sorting finded quadratic APN functions up to EA-equivalence using EA-invariants and function from [12]

REFERENCES: section 2.2 of the thesis
            subsection 2.3.3 of the thesis
            [14] - LPP-crypto/sboxu: Tools for studying S-boxes (https://github.com/lpp-crypto/sboxU)
            [2] - Christof Beirle and Gregor Leander - New instances of quadratic APN functions
            [3] - Christof Beirle, Marcus Brinkmann and Gregor Leander - Linearly self-equivalent APN permutations in small dimention
            [12] - Anne Canteaut, Alain Couvreur and L´eo Perrin - Recovering or Testing Extended-Affine Equivalence
                   - result of the paper is algorithm for finding matrices for EA-equivalence, which can be used more efficiently on the ortho-derivatives
                   - the source codes can be find: https://github.com/alaincouvreur/EA_equivalence_for_quadratic_functions
                   - IMPORTANT NOTE: for correct running it is neccessary to correct few mistakes in the source code of "Equivalence.py":
                        - line 335: "def full_tuple(F,G,A,B,m,n):"
                        - line 410: "return (None, None, None, None, None, None,0)"
                        - line 506: "abC = full_tuple(F, G, A_cand, B_cand,m,n)"

REQUIREMENTS:  - Python (version 3.12.2 was used), Python Software Foundation, https://www.python.org
               - Sage, the Sage MAthematical Software System (version 9.0 was used), Sage Developers, https://www.sagemath.org
			   - Cython Module for Python, Cython Developers, https://cython.readthedocs.io/en/latest/index.html
"""


# package "sboxU" contains usefull method for computing ortho-derivative of a vectorial Boolean function
from sboxU import *
# package sagemath
from sage.all import *
# value which cannot be obtained in some following functions
init_value = -1
# python code from [12] which contains main algorithm "get_equivalence" which can 
load("Equivalence.py")



"""
Based on which foud group of quadratic APN functions from the subsection 2.3.2 (attachment A.5) we want to examine we can ucomment desired following lines.



n = 2
load("chapter2/results-A5/chapter2-algorithm1-A3-output-n2.py")
functions = all_found_functions()


n = 3
load("chapter2/results-A5/chapter2-algorithm1-A3-output-n3.py")
functions = all_found_functions()



n = 4
load("chapter2/results-A5/chapter2-algorithm1-A4-output-n4.py")
functions = all_found_functions()



n = 5
load("chapter2/results-A5/chapter2-algorithm1-A4-output-n7.py")
functions = all_found_functions()



n = 6
load("chapter2/results-A5/chapter2-algorithm1-A4-output-n6.py")
functions = all_found_functions()



n = 7
load("chapter2/results-A5/chapter2-algorithm1-A4-output-n7.py")
functions = all_found_functions()


"""






# function for converting vectorial boolean function from the form of the look-up table in the polynomials
# INPUT:  look-up table of the boolean function
# OUTPUT: list of polynomials representing coordinate functions
def From_LUT_To_Polynomial(lut,n):
    # define the polynomial ring
    R = PolynomialRing(GF(2), 'x', n)
    x = R.gens()
    # initialize the output list 
    polynomials = [0 for _ in range(n)]
    # convert the boolean function to polynomial
    for i, output in enumerate(lut):
        # Convert i and output to binary
        inputs = bin(i)[2:].zfill(n)
        outputs = bin(output)[2:].zfill(n)
        for j, output in enumerate(outputs):
            if output == '1':
                term = 1
                for k, input in enumerate(inputs):
                    if input == '0':
                        term *= (1 - x[k])
                    else:
                        term *= x[k]
                polynomials[j] += term
    return polynomials


# function which checks if the given list "array" contains "init_value" at any position
# this function is mentioned in [1] but not in pseudo-code
# INPUT:  list which could contain "init_value" at any position
# OUTPUT: 0 if the input list contains "init_value" at any position
def Is_Complete(array):
    global init_value
    for i in range(0,len(array)):
        if(array[i] == init_value):
            return(0)
    return(1)


# function that returns the lowest index "i" where the given list contains the value "init_value"
# this function is mentioned in [1] but not in pseudo-code
# INPUT: list which could contain "init_value" at any position
# OUTPUT: index of the position where the value "init_value" is in the input list
def Next_Free_Position(array):
    global init_value
    for i in range(0,len(array)):
        if(array[i] == init_value):
            return(i)


# function prints all the found functions (classified by EA-equivalence classes) to the file "chapter2-algorithm2-output.txt"
# INPUT:  functions in the form of a look-up classified by EA-equivalence classes
# OUTPUT: None
def Print_Founded_Functions_Into_File(source,dimension,output_name):
    with open(output_name,'w') as output_file:
        output_file.write('The algorithm run with parametr for dimension equal to: ')
        output_file.write(str(dimension))
        output_file.write('\n')
        output_file.write('The number of found EA-equivalent classes: ')
        output_file.write(str(len(source)))
        output_file.write('\n')
        output_file.write('The found EA-equivalent classes are:\n')
    for i in range(len(source)):
        with open(output_name,'a') as output_file:
            output_file.write('The class number: ')
            output_file.write(str(i+1))
            output_file.write('\n')
        for j in range(len(source[i])):          
            with open(output_name,'a') as output_file:
                output_file.write(str(source[i][j]))
                output_file.write(',')
                output_file.write('\n')


# function that creates a list of binary notation of numbers 0,1,...,2**n-1
# INPUT:  integer "n"
# OUTPUT: list of binary notations of integers 0,...,2*n-1
def Create_Arrays_of_Binary_Notation_For_Inputs(n):
    binary_array = []
    for i in range (2**n):
        binary_string = bin(i)[2:]
        # Use list comprehension to create the bit array
        bit_array = [int(bit) for bit in binary_string]
        for _ in range(n-len(bit_array)):
            bit_array.insert(0,0)
        binary_array.append(bit_array)
    return(binary_array)


# function which convert bit notation in list to the integer
# INPUT: list with binary notation
# OUTPUT: integer
def Convert_Array_Of_Bits_To_Integer(list_for_convert):
    result = 0
    for ele in list_for_convert:
        result = (result << 1) | ele
    return(result)


# function which creates look-up table of the given function
# INPUT: - function in the form of a list consisting of coordinate functions, where each such coordinate function is in the form of a polynomial
#        - dimension n
# OUTPUT: look-up table
def Create_LUT_Table(G,n):
    binary_array = Create_Arrays_of_Binary_Notation_For_Inputs(n)
    G = Matrix(n,1,G)
    lut = []
    R = PolynomialRing(GF(2), 'x', n)
    zero_polynomial = R(0)
    zero = Matrix(1,1,zero_polynomial)(binary_array[0])
    one_polynomial = R(1)
    one = Matrix(1,1,one_polynomial)(binary_array[0])
    for i in range(2**n):
        array_of_bits = []
        for position in range(n):
            if (zero[0] == G(binary_array[i])[position]):
                array_of_bits.append(0)
            if (one[0] == G(binary_array[i])[position]):
                array_of_bits.append(1)
        lut.append(Convert_Array_Of_Bits_To_Integer(array_of_bits))
    return(lut)


# function to compute the Hamming weight of an integer which represents a vector from F_2^n
# INPUT:  integer "number"
# OUTPUT: number of ones in binary notation of "number"
def Hamming_Weight(number):
    one_count = 0
    number_in_bits = str(bin(number))
    for i in (number_in_bits):
        if (i == "1"):
            one_count+=1
    return(one_count)


# function that calculates the inner product of two numbers
# INPUT:  - integer "a"
#         - integer "b"
# OUTPUT: - integer "result"
def Inner_Product(a,b):
    # l&x is "l" AND "x", which sets bit to 1 if and only if both bits are 1
    # .bits() convert the integer into bits
    # .count(1) return number of "1"
    number_of_bits = Integer(a & b).bits().count(1)
    # result is from F_2
    result = number_of_bits%2
    return result


import numpy as np       
# function that calculates the extended Walsh spectrum for the given function "F" 
# INPUT:  - list "F", which is the look-up table of function F
#         - integer "n", which is the dimension of the input space of the function "F"
# OUTPUT: - list "extended_w_sp"
def Extended_Walsh_Spectrum(F,n):
    extended_w_sp = [0]*(2**n+1)
    for alpha in range(0,2**n):
        for beta in range(1,2**n):
            walsh_transform = 0
            for x in range(0,2**n):
                walsh_transform = walsh_transform + (-1)**((Inner_Product(beta,F[x])+Inner_Product(alpha,x))%2)
            extended_w_sp[abs(walsh_transform)] +=1
    return(extended_w_sp)


# function which, for given functions "F" and "G", compares their extended Walsh spectra
# INPUT:  - list "F", which is the look-up table of function F
#         - list "G", which is the look-up table of function G
#         - integer "n", which is the dimension of the input space of the functions "F" and "G"
# OUTPUT: - 0/1
def Have_Functions_Same_Ext_W_Spectrum(F,G,n):
    ext_w_sp_F = Extended_Walsh_Spectrum(F,n)
    ext_w_sp_G = Extended_Walsh_Spectrum(G,n)
    for i in range(2**n+1):
        if (ext_w_sp_F[i]!=ext_w_sp_G[i]):
            return(0)
    return(1)    
    

# function that calculates DDT for the given function "F"
# INPUT:  - list "F", which is the look-up table of function F
#         - integer "n", which is the dimension of the input space of the function "F"
# OUTPUT: - list "ddt"
def Compute_DDT(F,n):
    ddt = np.zeros((2**n, 2**n), dtype=int)
    for x in range(2**n):
        for alpha in range(1,2**n):
            beta = F[x] ^ F[x^alpha]
            ddt[alpha][beta] += 1
    return ddt


# function that calculates the differential spectrum for the given function "F"
# INPUT:  - list "F", which is the look-up table of function F
#         - integer "n", which is the dimension of the input space of the function "F"
# OUTPUT: - list "diff_spectrum"
def Compute_Diff_Spectrum(F,n):
    diff_spectrum = [0]*(2**n)
    ddt = Compute_DDT(F,n)
    for alpha in range(2**n):
        for beta in range(2**n):
            diff_spectrum[ddt[alpha][beta]]+=1
    return(diff_spectrum)


# function which, for given functions "F" and "G", compares their differential spectra
# INPUT:  - list "F", which is the look-up table of function F
#         - list "G", which is the look-up table of function G
#         - integer "n", which is the dimension of the input space of the functions "F" and "G"
# OUTPUT: - 0/1
def Have_Functions_Same_Diff_Spectrum(F,G,n):
    diff_sp_F = Compute_Diff_Spectrum(F,n)
    diff_sp_G = Compute_Diff_Spectrum(G,n)
    for i in range(2**n):
        if (diff_sp_F[i]!=diff_sp_G[i]):
            return(0)
    return(1)    


# function that compares extended Walsh spectra and differential spectra of given functions "F" and "G"
# if extended Walsh spectra and differential spectra are the same for both of the functions "F" and "G", then these functions can be EA-equivalent
# if extended Walsh spectra or differential spectra differs, then functions "F" and "G" cannot be EA-equivalent
# function using EA-invariants from section 2.2 of the thesis
# INPUT:  - list "F", which is the look-up table of function F
#         - list "G", which is the look-up table of function G
#         - integer "n", which is the dimension of input space of the functions "F" and "G"
# OUTPUT: - 0/1
def Can_Functions_Be_EA_equivalent(F,G,n):
    if (Have_Functions_Same_Diff_Spectrum(F,G,n))==0:
        return(0)
    if(Have_Functions_Same_Ext_W_Spectrum(F,G,n))==0:
        return(0)
    return(1)


# function that puts given functions "functions" into EA-equivalence classes using EA-invariants from section 2.2 of the thesis
# the function is described in more detail in the thesis
# INPUT:  - integer "n", which is the dimension of the input functions
#         - list "functions" containing the input functions in the form of a look-up table
# OUTPUT: list of input functions sorted into EA-classes
def Put_Functions_Into_EA_Classes(n,functions):
    EA_classes = []
    is_in_some_class = [init_value]*(len(functions))
    while(Is_Complete(is_in_some_class)==0):
        x = Next_Free_Position(is_in_some_class)
        EA_class = []
        EA_class.append(functions[x])
        is_in_some_class[x] = 1
        F = functions[x]
        # function ortho_derivative() from [12] and included in [14] 
        F_ortho = ortho_derivative(F)
        for i in range(x+1,len(functions)):
            if(is_in_some_class[i]==init_value):
                G = functions[i]
                G_ortho = ortho_derivative(G)
                if (Can_Functions_Be_EA_equivalent(F_ortho,G_ortho,n) & (Can_Functions_Be_EA_equivalent(F,G,n))):
                    A,B,C,D,E,H,I = get_equivalence(Matrix(n,1,From_LUT_To_Polynomial(F,n)), Matrix(n,1,From_LUT_To_Polynomial(G,n)),10,None,False)
                    if(A!=None):
                        EA_class.append(functions[i])
                        is_in_some_class[i] = 1
        EA_classes.append(EA_class.copy())
    return(EA_classes)








# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------- ALGORITHM 2 - modified - for classifying functions up to EA-equivalence using EA-invariants and function from [12] ------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

output_name = 'chapter2-algorithm2-A7-output-n'+str(n)+'.txt'
classes = Put_Functions_Into_EA_Classes(n,functions)
print('Number of EA-equivalence classes: ',len(classes))
Print_Founded_Functions_Into_File(classes,n,output_name)













              




"""
Some other functions that can help with handling of functions in the form of coordinate functions that are in the form of polynomials.
These functions have benn used when programming the functions in this Python file.


# function to evaluate if the condition "a\preceq b" is true for the vectors "a" and "b"
# INPUT:  vectors "a" and "b" in the form of integers 
# OUTPUT: 1 if the condition "a\preceq b" is true, otherwise 0
def Vector_Ordering(a,b):
    binary_string = bin(a)[2:]
    a_in_bits = [int(bit) for bit in binary_string]
    binary_string = bin(b)[2:]
    b_in_bits = [int(bit) for bit in binary_string]
    # length of both bit lists have to be same, thus the short one need to be extended with zeros
    while(len(a_in_bits) > len(b_in_bits)):
        b_in_bits.insert(0,0)
    while(len(a_in_bits) < len(b_in_bits)):
        a_in_bits.insert(0,0)
    for i in range(0,len(a_in_bits)):
        if(a_in_bits[i] != 0):
            if(b_in_bits[i] == 0):
                return(0)
    return(1)


# function decides if the given function "F" is APN or not
# INPUT:  - list "F", which is look-up table of function F
#         - integer "n", which is the dimension of input space of the functions "F"
# OUTPUT: - 0/1
def Is_Function_APN(F,n):
    ddt = Compute_DDT(F,n)
    for i in range(1, 2**n):
        for j in range(0, 2**n):
            if ddt[i][j] > 2:
                return False
    return True


# function decides if the given function "F" is quadratic
# INPUT:  - list "F", which is look-up table of function F
#         - integer "n", which is the dimension of input space of the functions "F"
# OUTPUT: - 0/1
def Is_Function_Quadratic(F,n):
    coef = [0]*(2**n)
    for u in range(2**n):
        if (Hamming_Weight(u)>2):
            au = 0
            for x in range(2**n):
                if (Vector_Ordering(x,u)==1):
                    au = au ^ F[x]
            coef[u]=au
    for i in range(2**n):
        if (coef[i]!=0):
            return(0)
    return (1)


# function which creates list of binary notation of number 0,1,...,2**n-1
# INPUT:  integer n
# OUTPUT: list of binary notation of integers 0,...,2*n-1, all binary notations have the same length equals to "n"
def Create_Arrays_of_Binary_Notation_For_Inputs(n):
    binary_array = []
    for i in range (2**n):
        binary_string = bin(i)[2:]
        bit_array = [int(bit) for bit in binary_string]
        for _ in range(n-len(bit_array)):
            bit_array.insert(0,0)
        binary_array.append(bit_array)
    return(binary_array)


# function which convert bit notation in list to the integer
# INPUT: list with binary notation
# OUTPUT: integer
def Convert_Array_Of_Bits_To_Integer(list_for_convert):
    result = 0
    for ele in list_for_convert:
        result = (result << 1) | ele
    return(result)


# function which creates look-up table of given function
# INPUT: - function in the form of list which consists of coordinate functions where each such coordinate function is of the form of polynomial
#        - dimension n
# OUTPUT: look-up table
def Create_LUT_Table(G,n):
    binary_array = Create_Arrays_of_Binary_Notation_For_Inputs(n)
    G = Matrix(n,1,G)
    lut = []
    R = PolynomialRing(GF(2), 'x', n)
    zero_polynomial = R(0)
    zero = Matrix(1,1,zero_polynomial)(binary_array[0])
    one_polynomial = R(1)
    one = Matrix(1,1,one_polynomial)(binary_array[0])
    for i in range(2**n):
        array_of_bits = []
        for position in range(n):
            if (zero[0] == G(binary_array[i])[position]):
                array_of_bits.append(0)
            if (one[0] == G(binary_array[i])[position]):
                array_of_bits.append(1)
        lut.append(Convert_Array_Of_Bits_To_Integer(array_of_bits))
    return(lut)


"""

