"""
---------------
Attachment A.13
---------------

Algorithm 5 - algorithm for classification of n+1-bit quadratic APN functions with maximum liunearity up to EA-equivalence 

REFERENCES: [5] - Christof Beirle, Gregor Leander and Leo Perrin - Trims and extensions of quadratic APN functions
			[14] - LPP-crypto/sboxu: Tools for studying S-boxes (https://github.com/lpp-crypto/sboxU)
			[A.5] - Attachment A.5 to the thesis

REQUIREMENTS:  - Python (version 3.12.2 was used), Python Software Foundation, https://www.python.org
               - Sage, the Sage MAthematical Software System (version 9.0 was used), Sage Developers, https://www.sagemath.org
			   - Cython Module for Python, Cython Developers, https://cython.readthedocs.io/en/latest/index.html
"""

# using tools provided by [14]
from sboxU import *


# this function is taken from [Attachment B, 5]
# function returns bit on "i"-th position in "a"
# INPUT:  - integer a
#         - integer i
# OUTPUT: integer 0 or 1
def bit(a,i):
	return((a>>i)%2)


# this function is taken from [Attachment B, 5]
# function which calculates inner product for two binary vectors represented as integers
# INPUT:  - integer l
#         - integer x
# OUTPUT: integer 0 or 1
def innerProduct(l,x):
	return Integer(l&x).bits().count(1)%2


# this function is taken from [Attachment B, 5]
# function which converts the matrix A to 2-dimensional array according to the theory from Subsection 5.4.1 in the thesis
# INPUT:  - list "A"
# OUTPUT: matrix T
def matrixToSbox(A):
	T = []
	V = VectorSpace(GF(2),A.nrows())
	for v in V:
		T. append(ZZ(list(A*v),base=2))
	return T


# this function is taken from [Attachment B, 5] and described more preciously in the thesis
# function for computing M (for computing L) with respect to the input (quadratic APN vectorial Boolean function G, ortho-derivative of G and function \ell)
# INPUT:  - integer l
#         - integer x
# OUTPUT: integer 0 or 1
def genSystem(G,ortho,l):
	M = []
	for alpha in range(2**n)[1::]:
		if(innerProduct(l,alpha)==0):
			M_alpha = []
			for bit_pi in range (n):
				for bit_alpha in range (n):
					M_alpha.append(bit(alpha,bit_alpha)*bit(ortho[alpha],bit_pi))
			M.append(M_alpha)
	return(matrix(GF(2),M))


# function for computing the function T in dimension n+1 with rspect to the input (quadratic APN vectorial Boolean function G, ortho-derivative of G and function \ell)
# this function is taken from [Attachment B, 5] and described more preciously in the thesis
# INPUT:  - list "G" which is a look-up table for function G
#		  - list "ortho" which is a look-up table for function ortho-derivative of a function G
#         - integer "l" which represents binary notation of a vector, which represents linear mapping \ell 
# OUTPUT: list "T"
def isExtendable(G,ortho,l):
	M = genSystem(G,ortho,l)
	v1 = vector(GF(2),[1]*(M.nrows()))
	try :
		L_tilde = M.solve_right(v1)
	except (ValueError) :
		return []
	if (M.right_kernel().dimension()>(2*n)):
		print ( "There might be more solutions! \n")
	#otherwise there is only on L up to Gamma−equivalence
	L = matrixToSbox (matrix(GF( 2 ) , n , n , list (L_tilde ) ) )
	T = [ ]
	for i in range (2**n) :
		T.append (int (G[i] << 1))
	for i in range (2**n) :
		T.append(int(((G[i]^L[i])<<1)^innerProduct(l,i)))
	return T


# this function is taken from [Attachment B, 5] and described more preciously in the thesis
# function returns a list of APN extensions of G
# INPUT:  - list "G" which is a look-up table for function G
# OUTPUT: list "sol"
def extensions (G) :
	sol = [ ]
	pi = ortho_derivative(G)
	for l in range (2**n) :
		T = isExtendable (G, pi , l )
		if (not (T == [])):
			sol.append(T)
	return sol


# function prints all the functions that are found in the file "chapter5-algorithm-output.txt.txt"
# INPUT: - list "source" containing the functions in the form of a look-up table
#        - integer "dimension"
# OUTPUT: None
def Print_Founded_Functions_Into_File(source,dimension):
	with open('chapter5-algorithm-output.txt','w') as output_file:
		output_file.write('The algorithm run with parametr for dimension equal to: ')
		output_file.write(str(dimension))
		output_file.write('\n')
		output_file.write('The number of found functions: ')
		output_file.write(str(len(source)))
		output_file.write('\n')
		output_file.write('The found functions are:\n')
	for i in range(len(source)):
		with open('chapter5-algorithm-output.txt','a') as output_file:
			output_file.write(str(source[i][0]))
			output_file.write(',')
			output_file.write('\n')




# ------------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------- Algorithm 5 - Classification of n+1-bit quadratic APN functions with maximum liunearity up to EA-equivalence --------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------
# the algorithm is taken from [Attachment B, 5]

"""
Based on which foud group of quadratic APN functions from the subsection 2.3.2 (attachment A.5) or from [14] we want to examine we can ucomment desired following lines.



n = 2
# the look-up table in "functions" are associated with the representatives of EA-equivalence class found in [A.5], those representatives are first functions of the EA-equivalence class 
functions = [[0, 1, 0, 2]]



n = 3
functions = [[0, 4, 3, 1, 6, 6, 0, 6]]



n = 4
functions = [[0, 0, 11, 0, 3, 0, 12, 4, 0, 1, 2, 8, 5, 7, 3, 10]]



n = 5
functions = [
[0, 12, 3, 18, 8, 24, 18, 31, 3, 12, 22, 4, 0, 19, 12, 2, 27, 30, 13, 21, 0, 25, 15, 11, 8, 14, 8, 19, 24, 2, 1, 6],
[0, 12, 3, 18, 8, 24, 18, 31, 3, 12, 22, 4, 0, 19, 12, 2, 27, 30, 13, 21, 0, 25, 15, 11, 22, 16, 22, 13, 6, 28, 31, 24]
]



n = 6
functions = [
[0, 0, 3, 25, 23, 48, 47, 18, 13, 3, 7, 19, 22, 63, 39, 20, 50, 61, 14, 27, 19, 59, 20, 38, 13, 12, 56, 35, 32, 6, 46, 18, 62, 53, 59, 42, 10, 38, 52, 2, 12, 9, 0, 31, 52, 22, 3, 59, 48, 52, 10, 20, 50, 17, 51, 10, 48, 58, 3, 19, 62, 19, 54, 1],
[0, 0, 3, 25, 23, 48, 47, 18, 13, 3, 7, 19, 22, 63, 39, 20, 50, 61, 14, 27, 19, 59, 20, 38, 13, 12, 56, 35, 32, 6, 46, 18, 62, 53, 59, 42, 10, 38, 52, 2, 21, 16, 25, 6, 45, 15, 26, 34, 6, 2, 60, 34, 4, 39, 5, 60, 31, 21, 44, 60, 17, 60, 25, 46],
[0, 0, 3, 25, 23, 48, 47, 18, 13, 3, 7, 19, 22, 63, 39, 20, 50, 61, 14, 27, 19, 59, 20, 38, 13, 12, 56, 35, 32, 6, 46, 18, 62, 53, 59, 42, 10, 38, 52, 2, 35, 38, 47, 48, 27, 57, 44, 20, 48, 52, 10, 20, 50, 17, 51, 10, 31, 21, 44, 60, 17, 60, 25, 46],
[0, 0, 3, 25, 23, 48, 47, 18, 13, 3, 7, 19, 22, 63, 39, 20, 50, 61, 14, 27, 19, 59, 20, 38, 13, 12, 56, 35, 32, 6, 46, 18, 62, 53, 59, 42, 48, 28, 14, 56, 17, 20, 29, 2, 19, 49, 36, 28, 2, 6, 56, 38, 58, 25, 59, 2, 31, 21, 44, 60, 43, 6, 35, 20],
[0, 0, 3, 25, 23, 48, 47, 18, 13, 3, 7, 19, 22, 63, 39, 20, 50, 61, 14, 27, 19, 59, 20, 38, 13, 12, 56, 35, 32, 6, 46, 18, 62, 53, 59, 42, 49, 29, 15, 57, 7, 2, 11, 20, 4, 38, 51, 11, 21, 17, 47, 49, 44, 15, 45, 20, 30, 20, 45, 61, 43, 6, 35, 20],
[0, 0, 3, 25, 23, 48, 47, 18, 13, 3, 7, 19, 22, 63, 39, 20, 50, 61, 14, 27, 19, 59, 20, 38, 13, 12, 56, 35, 32, 6, 46, 18, 62, 53, 59, 42, 55, 27, 9, 63, 30, 27, 18, 13, 27, 57, 44, 20, 13, 9, 55, 41, 50, 17, 51, 10, 31, 21, 44, 60, 44, 1, 36, 19],
[0, 0, 26, 30, 53, 22, 53, 18, 20, 34, 63, 13, 58, 47, 11, 26, 6, 8, 12, 6, 33, 12, 49, 24, 5, 61, 62, 2, 57, 34, 24, 7, 38, 31, 18, 47, 29, 7, 51, 45, 49, 62, 52, 63, 17, 61, 14, 38, 47, 24, 11, 56, 6, 18, 56, 40, 47, 46, 58, 63, 29, 63, 18, 52],
[0, 0, 49, 8, 28, 20, 54, 7, 31, 28, 11, 49, 36, 47, 43, 25, 18, 16, 4, 63, 52, 62, 57, 10, 36, 37, 23, 47, 37, 44, 13, 61, 19, 58, 29, 13, 34, 3, 55, 47, 63, 21, 20, 7, 41, 11, 25, 2, 60, 23, 21, 7, 55, 20, 5, 31, 57, 17, 53, 36, 21, 53, 2, 27],
[0, 8, 25, 24, 0, 15, 19, 21, 28, 8, 54, 43, 9, 26, 41, 51, 1, 1, 44, 37, 22, 17, 49, 63, 53, 41, 43, 62, 55, 44, 35, 49, 36, 30, 42, 25, 45, 16, 41, 29, 7, 33, 58, 21, 27, 58, 44, 4, 20, 38, 46, 21, 10, 63, 58, 6, 31, 49, 22, 49, 20, 61, 23, 55],
[0, 17, 62, 28, 35, 19, 38, 37, 14, 19, 34, 12, 58, 6, 45, 34, 41, 10, 55, 39, 15, 13, 42, 27, 55, 24, 59, 39, 6, 8, 49, 12, 53, 2, 61, 57, 8, 30, 59, 30, 28, 39, 6, 14, 54, 44, 23, 62, 56, 61, 16, 38, 0, 36, 19, 4, 1, 8, 59, 1, 46, 6, 47, 52],
[0, 26, 9, 21, 38, 13, 2, 47, 16, 24, 32, 46, 19, 42, 14, 49, 15, 61, 10, 62, 0, 3, 40, 45, 32, 0, 28, 58, 10, 27, 27, 12, 13, 36, 32, 15, 12, 20, 12, 18, 32, 27, 52, 9, 4, 14, 61, 49, 13, 12, 44, 43, 37, 21, 41, 31, 31, 12, 7, 18, 18, 48, 39, 3],
[0, 36, 12, 17, 58, 5, 28, 26, 6, 40, 24, 15, 14, 59, 58, 54, 46, 5, 30, 12, 51, 3, 41, 32, 10, 43, 40, 48, 37, 31, 45, 46, 12, 6, 49, 2, 44, 61, 59, 19, 20, 20, 59, 2, 6, 29, 3, 33, 26, 31, 27, 39, 29, 3, 54, 17, 32, 47, 51, 5, 21, 1, 44, 1],
[9, 46, 15, 38, 19, 8, 43, 62, 26, 56, 58, 22, 6, 24, 24, 8, 3, 8, 3, 6, 56, 15, 6, 63, 63, 49, 25, 25, 2, 48, 26, 38, 20, 63, 9, 44, 37, 50, 6, 31, 18, 60, 41, 9, 37, 55, 32, 60, 17, 22, 10, 3, 1, 58, 36, 17, 56, 58, 5, 9, 46, 16, 45, 29]
]



n = 7
# using tools provided by [14]
load("sboxU/known_functions/sevenBitAPN.py")
functions = all_quadratics()
"""


solutions = [ ]
for i in range(len(functions)):
	print(i)
	T = extensions (functions[i])
	if (not(T==[ ])) :
		solutions.append(T)
print (solutions)
print(len(solutions))
Print_Founded_Functions_Into_File(solutions,n)
