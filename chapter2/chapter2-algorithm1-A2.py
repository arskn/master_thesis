"""
Attachment A2

Algorithm 1 - algorithm for finding any quadratic APN vectorial Boolean function in a given dimension n, with a time limit and a number of iterations

REFERENCES: subsection 2.3.1
			[8] - LPP-crypto/sboxu: Tools for studying S-boxes (https://github.com/lpp-crypto/sboxU)
			[9] - Christof Beirle and Gregor Leander - New instances of quadratic APN functions
	    	[10] - Christof Beirle, Marcus Brinkmann and Gregor Leander - Linearly self-equivalent APN permutations in small dimention
"""


# "random" is used to generate random values in the functions
import random
# package "sboxU" contains useful method Integer [8]
from sboxU import Integer
# package to access the system time
import time


# function to create a list "sbox"
# INPUT:  dimension n
# OUTPUT: list of "init_value" repeated 2**n times
def Generate_Sbox(n) :
	global init_value
	sbox = [init_value]*(2**n)
	return(sbox)


# function to generate a 2-dimensional list "P", containing 2**n permutations
# INPUT:  dimension n
# OUTPUT: list of 2**n permutations on 2**n elements
def Generate_P(n):
	# random.seed()
	P = []
	element_of_P = []	
	for j in range(0,2**n):
		element_of_P.append(j)
	for i in range(0,2**n):
		random.shuffle(element_of_P)
		P.append(element_of_P.copy())
	return(P)


# function to generate a 2-dimensional list "DDT" to check that the APN condition is still satisfied during the execution of the algorithm
# INPUT:  dimension n
# OUTPUT: list of 2**n lists of 2**n elements which are initialised at 0
def Generate_DDT(n):
	DDT = []
	element_of_DDT=[0]*(2**n)
	for i in range(0,2**n):
		DDT.append(element_of_DDT.copy())
	return(DDT)


# function for generating the list "CTR", which is used to check if all the vectors "x\preceq u" have been included in the xor for the coefficient a_u of ANF
# INPUT:  dimension n
# OUTPUT: two lists, each initialised at 0 repeated 2**n times
def Generate_SUM_and_CTR(n):
	SUM = [0]*(2**n)
	CTR = [0]*(2**n)
	return(SUM, CTR)


# function to compute Hamming weight of an integer which representing a vector from F_2^n
# INPUT:  integer "number"
# OUTPUT: number of ones in binary notation of "number"
def Hamming_Weight(number):
	one_count = 0
	number_in_bits = str(bin(number))
	for i in (number_in_bits):
		if (i == "1"):
			one_count+=1
	return(one_count)


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


# function to add information to the DDT list while the algorithm is running
# this function is taken from [10] and is described in more detail in the thesis
# INPUT: integer "x" representing the input of the function
# OUTPUT: 1 if the function is still APN, otherwise 0
def Add_DDT_Information(x):
	global sbox
	global DDT
	global init_value
	for alpha in range(1,2**n):
		if (sbox[x^alpha] != init_value):
			DDT[alpha][sbox[x]^sbox[x^alpha]] = DDT[alpha][sbox[x]^sbox[x^alpha]] + 2
			if (DDT[alpha][sbox[x]^sbox[x^alpha]] > 2):
				return(0)
	return(1)


# function to remove information from the DDT list while the algorithm is running
# this function is taken from [10]
# INPUT: integer "x" representing the input of the function
# OUTPUT: 0
def Remove_DDT_Information(x):
	global sbox
	global DDT
	global init_value
	for alpha in range(1,2**n):
		if(sbox[x^alpha] != init_value):
			DDT[alpha][sbox[x]^sbox[x^alpha]] = DDT[alpha][sbox[x]^sbox[x^alpha]] - 2
			if (DDT[alpha][sbox[x]^sbox[x^alpha]] == 2):
				return(0)
	return(1)


# function that checks if the condition for coefficients in ANF is not violated (Observation 14 from the thesis)
# this function is taken from [9] and is described in more detail in the thesis
# INPUT: integer "x" representing the input of the function
# OUTPUT: 1 if the function is still quadratic, otherwise 0
def Add_Degree_Information(x):
	global CTR
	global SUM
	global sbox
	for u in range(1,2**n):
		if(Hamming_Weight(u) > 2):
			if(Vector_Ordering(x,u) == 1):
				CTR[u]+=1
				SUM[u] = SUM[u]^sbox[x]
				if(CTR[u] == 2**Hamming_Weight(u)):
					if(SUM[u] != 0):
						return(0)
	return(1)


# function that removes sbox[u] from the xor for the coefficients of the ANF during the step from "depth" to "depth-1" in the "nextVal" function
# this function is taken from [9]
# INPUT: integer "x" representing the input of the function
# OUTPUT: 0/1
def Remove_Degree_Information(x):
	global CTR
	global SUM
	global sbox
	for u in range(1,2**n):
		if (Hamming_Weight(u) > 2):
			if(Vector_Ordering(x,u) == 1):
				CTR[u]-=1
				SUM[u] = SUM[u]^sbox[x]
				if(CTR[u] == (2**Hamming_Weight(u)-1)):
					if (SUM[u] != sbox[x]):
						return(0)
	return (1)


# function tries, if the point "x" with value "sbox[x]" does not violate the APN condition
# this function is taken from [9] and is described in more detail in the thesis
# INPUT: integer "x" representing the input of the function
# OUTPUT: 1 if the "x" was successfully added, otherwise 0
def Add_Point(x):
	if (Add_DDT_Information(x)):
		return (Add_Degree_Information(x))
	return(0)


# function that removes the points "x" with value "sbox[x]" 
# this function is taken from [9]
# INPUT: integer "x" representing the input to the function
# OUTPUT: None
def Remove_Point(x):
	if (Remove_DDT_Information(x)):
		Remove_Degree_Information(x)


# function that checks if the given list "array" contains "init_value" at any position
# this function is mentioned in [9]
# INPUT: list "array"
# OUTPUT: 0 if the list contains the value "init_value", otherwise 1
def Is_Complete(array):
	global init_value
	for i in range(0,len(array)):
		if(array[i] == init_value):
			return(0)
	return(1)

# function to check that the algorithm run time does not exceed the time limit set by user (the user input time limit for each run of the algorithm)
# INPUT: time which we want to check
# OUTPUT: 0 if the time limit has been exceeded, otherwise 1
def Is_Time_Limit_Over(time):
	global time_limit
	global time_start
	if(float(time-time_start)>time_limit):
		return(1)
	return(0)


# function to find a quadratic APN function by setting a new value for "x" to "sbox[x]" and checking, if the function can still be quadratic and APN
# function includes the restriction provided by the "Is_Time_Limit_Over" function
# this function is taken from [9] and described in more detial in the thesis; however, it has been modified, because in [9] there were mystakes that made the function unexecutable
# INPUT: integer "x"
# OUTPUT: None
def Next_Val(x):
	global output_sbox
	global sbox
	global P
	global init_value
	global time_limit
	if (Is_Complete(sbox)):
		output_sbox.append(sbox.copy())
		return
	for z in range(0,2**n): 
		y = P[x][z]
		sbox[x] = y
		b = Add_Point(x)
		if(b == 1):
			Next_Val(x+1)
		if (output_sbox != []):
			return	
		Remove_Point(x)
		sbox[x] = init_value
		if(Is_Time_Limit_Over(time.time())):
			return	
	return




# ----------------------------------------------------------------------------------
# ----------------------------- ALGORITHM 1 - modified -----------------------------
# ----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# 1st PART - dialog with the user to obtain the parameters for the algorithm

# dialog with the user to obtain the following information: dimension n, for how long (in seconds) will the tree search run, how many runs will be made
print('This algorithm is used to find a quadratic APN function. First select the following parameter:')

print('1) Enter the the dimension in which you want to search for the quadratic APN vectorial boolean function. (integer)')
n = int(input())

print('2) Enter the maximum number of seconds a run will take (positive real number).')
time_limit = float(input())

print('3) Enter the number of runs to be made (each run will take ',time_limit,' seconds).')
number_of_runs = int(input())

#------------------------------------------------------------------------
# 2nd PART - setting up variables which are used throughout the algorithm

output_sbox = []
init_value = -1

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3rd PART - the algorithm itself; first the algorithm sets up the variables to be used in the i-th run, and then the recursive tree search is executed by nextVal(0)

run = 0
time_start_of_whole = time.time()
# the algorithm checks if a function has been found and if the number of iterations has not been exceeded
while ((output_sbox == []) & (run<number_of_runs)):
	run = run + 1
	time_start = time.time()
	sbox = Generate_Sbox(n)
	rnd = random.random()
	random.seed(rnd)
	P = Generate_P(n)
	DDT = Generate_DDT(n)
	SUM, CTR = Generate_SUM_and_CTR(n)
	Next_Val(0)
print(rnd)

if (output_sbox != []):
	print('The quadratic APN vectorial Boolean function has been found. The lookup table is:')
	print(output_sbox)
else:
	print('No quadratic APN vectorial Boolean function has been found.')
