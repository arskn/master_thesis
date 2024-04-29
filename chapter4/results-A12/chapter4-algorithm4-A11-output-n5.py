# Attachment A.12 for n=5 to 4

# All APN functions found using the algorithm A.11 for 2 EA-inequivalent quadratic APN functions G listed in "chapter4-algorithm4-A11-output-n6.py".


"""
REFERENCES: Chapter 4 of the thesis
            Chapter 2 of the thesis
            [1] - Christof Beirle, Gregor Leander and Leo Perrin - Trims and extensions of quadratic APN functions
            [8] - LPP-crypto/sboxu: Tools for studying S-boxes (https://github.com/lpp-crypto/sboxU)
"""

"""
The number of input functions that have an APN trim in their trim spectrum: 2
The number of input functions which have no APN trim in their trim spectrum: 0
The dimension of the functions G: 4
The number of quadratic functions among the founded APN trims: 2 (using the observation from Chapter 2 of the thesis)
The number of EA-equivalence classes of functions G: 1 (using the "Algorithm 2 - modified")
"""



"""
The file contains the following functions that return a list of functions:
- classified_functions_G()
        - function returns representatives of EA-equivalence classes of functions G

        
Remaining results:
The file contains the following functions that return a list of functions:
- functions_G()
        - function returns list of all functions G:F_2^4 -> F_2^4

- parameters_of_APN_trims()
        - parameters alpha, beta, gamma, epsilon for which we could find APN trims

- functions_with_APN_trim()
        - list of input functions (from the file "chapter4-algorithm4-output-n7_to_n5-reduced.py") with APN functions in their trim spectrum

- functions_without_APN_trim()
        - list of input functions (from the file "chapter4-algorithm4-output-n7_to_n5-reduced.py") without APN functions in their trim spectrum
        
- APN_trims()
        - list of trims of a functions as they were obtained during execution of the algorithm (list containing 2^5 values, of which 2^4 are equal to init_value=-1)
"""


# representatives of EA-equivalence classes
def classified_functions_G():
    return[
[0, 4, 0, 13, 12, 10, 9, 6, 1, 9, 15, 14, 5, 15, 14, 13]
]


# functions G which are obtained from the APN trims
def functions_G():
    return[
[0, 4, 0, 13, 12, 10, 9, 6, 1, 9, 15, 14, 5, 15, 14, 13],
[0, 8, 5, 0, 4, 7, 5, 11, 4, 0, 11, 2, 13, 2, 6, 4]
]

# parameters \alpha, \beta, \gamma, \epsilon for which we could find APN trim
def parameters_of_APN_trims():
    return[
[1, 1, 1, 0],
[1, 1, 1, 0]
    ]

# input functions with APN functions in their trim spectrum
def functions_with_APN_trim():
    return[
[0, 4, 16, 12, 6, 19, 20, 25, 0, 18, 31, 21, 18, 17, 15, 20, 5, 24, 19, 22, 10, 6, 30, 10, 22, 29, 15, 28, 13, 23, 22, 20],
[0, 4, 4, 24, 5, 4, 22, 15, 20, 18, 27, 5, 21, 22, 13, 22, 2, 20, 0, 14, 29, 14, 8, 3, 0, 20, 9, 5, 27, 10, 5, 12]
]

# input functions without APN functions in their trim spectrum
def functions_without_APN_trim():
    return []

# APN trims
def APN_trims():
    return[
[0, -1, 16, -1, 6, -1, 20, -1, 0, -1, 30, -1, 18, -1, 14, -1, 4, -1, 18, -1, 10, -1, 30, -1, 22, -1, 14, -1, 12, -1, 22, -1],
[0, -1, 4, -1, 4, -1, 22, -1, 20, -1, 26, -1, 20, -1, 12, -1, 2, -1, 0, -1, 28, -1, 8, -1, 0, -1, 8, -1, 26, -1, 4, -1]
]
