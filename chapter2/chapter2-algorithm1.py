"""
REFERENCES: [9] - Christof Beirle and Gregor Leander - New instances of quadratic APN functions
	    	[10] - Christof Beirle, Marcus Brinkmann and Gregor Leander - Linearly self-equivalent APN permutations in small dimention
"""


# "random" serves for generating random values in the functions
import random
# package "sboxU" contains usefull method for computing ortho-derivative of a vectorial Boolean function
from sboxU import Integer, ortho_derivative



# function for generating list "sbox"
def Generate_Sbox(n) :
	global init_value
	sbox = [init_value]*(2**n)
	return(sbox)


# function for generating 2-dimensional list "P", which includes 2**n permutations
def Generate_P(n):
	P = []
	random.seed(1) #TODO REMOVE
	element_of_P = []	
	for j in range(0,2**n):
		element_of_P.append(j)
	for i in range(0,2**n):
		random.shuffle(element_of_P)
		P.append(element_of_P.copy())
	return(P)


# function for generating 2-dimensional list "DDT", for checking if the condition of APN is still satisfied during the run of the algorithm 
def Generate_DDT(n):
	DDT = []
	element_of_DDT=[0]*(2**n)
	for i in range(0,2**n):
		DDT.append(element_of_DDT.copy())
	return(DDT)


# function for generating list "CTR", which servers for checking if all vectors "x\preceq u" were included in the xor for the coefficient a_u of ANF
def Generate_SUM_and_CTR(n):
	SUM = [0]*(2**n)
	CTR = [0]*(2**n)
	return(SUM, CTR)


# function for computing Hamming weight of integer which represents a vector from F_2^n
def Hamming_Weight(number):
	one_count = 0
	number_in_bits = str(bin(number))
	for i in (number_in_bits):
		if (i == "1"):
			one_count+=1
	return(one_count)


# function for evaluating if the condion "a\preceq b"" is satisfied for the vectors "a" and "b"
def Vector_Ordering(a,b):
	a_in_bits = str(bin(a))
	b_in_bits = str(bin(b))
	# length of both bit lists have to be same, thus the short one need to be extended with zeros
	while(len(a_in_bits) > len(b_in_bits)):
		b_in_bits = b_in_bits[:2]+ '0'+ b_in_bits[2:]
	while(len(a_in_bits) < len(b_in_bits)):
		a_in_bits = a_in_bits[:2]+ '0'+ a_in_bits[2:]
	for i in range(0,len(a_in_bits)):
		if(a_in_bits[i] != 0):
			if(b_in_bits[i] == 0):
				return(0)
	return(1)


# function for add information into DDT list during the run of the algorithm
# this function is taken from [10] and described more preciously in the thesis
def Add_DDT_Information(x):
	global sbox
	global DDT
	global init_value
	for alpha in range(1,2**n):
		if(Hamming_Weight(alpha)%2 == 0):
			if (sbox[x^alpha] != init_value):
				DDT[alpha][sbox[x]^sbox[x^alpha]] = DDT[alpha][sbox[x]^sbox[x^alpha]] + 2
				if (DDT[alpha][sbox[x]^sbox[x^alpha]] > 2):
					return(0)
	return(1)


# function for remove information from DDT list during the run of the algorithm
# this function is taken from [10]
def Remove_DDT_Information(x):
	global sbox
	global DDT
	global init_value
	for alpha in range(1,2**n):
		if(Hamming_Weight(alpha)%2 == 0):
			if(sbox[x^alpha] != init_value):
				DDT[alpha][sbox[x]^sbox[x^alpha]] = DDT[alpha][sbox[x]^sbox[x^alpha]] - 2
				if (DDT[alpha][sbox[x]^sbox[x^alpha]] == 2):
					return(0)


# function which is checking, if the condition for coefficients in ANF is not violated
# this function is taken from [9] and described more preciously in the thesis
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


# function which removes sbox[u] from the xor for coefficients of the ANF during the step from "depth" into "depth-1" in the function "nextVal" 
# this function is taken from [9]
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
# this function is taken from [9] and described more preciously in the thesis
def Add_Point(x):
	if (Add_DDT_Information(x)):
		return (Add_Degree_Information(x))
	return(0)


# function which removes the points "x" with value "sbox[x]" 
# this function is taken from [9]
def Remove_Point(x):
	if (Remove_DDT_Information(x)):
		Remove_Degree_Information(x)


# function which checks if the given list "array" contains "init_value" at any position
# this function is mentioned in [9]
def Is_Complete(array):
	global init_value
	for i in range(0,len(array)):
		if(array[i] == init_value):
			return(0)
	return(1)


# function for finding APN function by setting new value to "sbox[x]" for "x" and checking, if the function can still be APN
# this function is taken from [9] and described more preciously in the thesis; altough it has benn modified, because in [9] were mystakes, which made the function unrunnable 
def Next_Val(x):
	global output_sboxes
	global sbox
	global P
	global init_value
	global time_limit
	if (Is_Complete(sbox)):
		output_sboxes.append(sbox.copy())
		print("is complete")
		return
	for z in range(0,2**n): 
		y = P[x][z]
		sbox[x] = y
		b = Add_Point(x)
		if(b == 1):
			Next_Val(x+1)
		if (output_sboxes != []):
			return	
		Remove_Point(x)
		sbox[x] = init_value	
	return


# function for computing the inner product of two integeres (which represent the vector from F_2^n)
def Inner_Product(l,x):
	return Integer(l&x).bits().count(1)%2







# ------------------------------------------------------------------------
# ----------------------------- ALGORITHM 1 ------------------------------
# ------------------------------------------------------------------------	
#-------------------------------------------------------------------------
# 1st PART - dialog with user to get parameters for the algorithm

# dialog with user to get the following information: dimension n, for how long (in seconds) will the tree search run, how many runs will be made
print('This algorithm serves for finding a quadratic APN function. At first, choose the following parameter:')

print('Enter the the dimension in which will be the tree search made. (integer)')
n = int(input())

#------------------------------------------------------------------------
# 2nd PART - setting up variables which are used for the whole algorithm

output_sboxes = []
init_value = -1

#-------------------------------------------------------------------------------------------------------------------------------------------------------
# 3rd PART - algorithm itself; at first the algorithm set up variables used in the i-th run and then the recursive tree search is executed by nextVal(0)
sbox = Generate_Sbox(n)
P = Generate_P(n)
DDT = Generate_DDT(n)
SUM, CTR = Generate_SUM_and_CTR(n)
Next_Val(0)

print('The quadratic APN vectorial Boolean function was found. The lookup table is:')
print(output_sboxes)



















