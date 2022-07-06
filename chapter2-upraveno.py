"""
REFERENCES: [1] - Christof Beirle and Gregor Leander - New instances of quadratic APN functions
			[2] - Christof Beirle, Marcus Brinkmann and Gregor Leander - Linearly self-equivalent APN permutations in small dimention
"""


# "random" serves for generating random values in the functions
import random
# package "sboxU" contains usefull method for computing ortho-derivative of a vectorial Boolean function
from sboxU import *
# package for access to the system time
import time
# TODO - remove the following package
# from sage.all import GF, PolynomialRing


# function for generating list "sbox"
def generateSbox(n) :
	global init_value
	sbox = [init_value]*(2**n)
	return(sbox)


#function for generating 2-dimensional list "P", which includes 2**n permutations
def generateP(n):
	P = []
	#random.seed(1) #TODO REMOVE
	element_of_P = []	
	for j in range(0,2**n):
		element_of_P.append(j)
	for i in range(0,2**n):
		random.shuffle(element_of_P)
		P.append(element_of_P.copy())
	return(P)


# function for generating 2-dimensional list "DDT", for checking if the condition of APN is still satisfied during the run of the algorithm 
def generateDDT(n):
	DDT = []
	element_of_DDT=[0]*(2**n)
	for i in range(0,2**n):
		DDT.append(element_of_DDT.copy())
	return(DDT)


# function for generating list "CTR", which servers for checking if all vectors "x\preceq u" were included in the xor for the coefficient a_u of ANF
def generateSUMandCTR(n):
	SUM = [0]*(2**n)
	CTR = [0]*(2**n)
	return(SUM, CTR)


# function for computing Hamming weight of integer which represents a vector from F_2^n
def hammingWeight(number):
	one_count = 0
	number_in_bits = str(bin(number))
	for i in (number_in_bits):
		if (i == "1"):
			one_count+=1
	return(one_count)


# function for evaluating if the condion "a\preceq b"" is satisfied for the vectors "a" and "b"
def vectorOrdering(a,b):
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
# this function is taken from [2] and described more preciously in the thesis
def addDDTInformation(x):
	global sbox
	global DDT
	global init_value
	for alpha in range(1,2**n):
		if(hammingWeight(alpha)%2 == 0):
			if (sbox[x^alpha] != init_value):
				DDT[alpha][sbox[x]^sbox[x^alpha]] = DDT[alpha][sbox[x]^sbox[x^alpha]] + 2
				if (DDT[alpha][sbox[x]^sbox[x^alpha]] > 2):
					return(0)
	return(1)


# function for remove information from DDT list during the run of the algorithm
# this function is taken from [2] and described more preciously in the thesis
def removeDDTInformation(x):
	global sbox
	global DDT
	global init_value
	for alpha in range(1,2**n):
		if(hammingWeight(alpha)%2 == 0):
			if(sbox[x^alpha] != init_value):
				DDT[alpha][sbox[x]^sbox[x^alpha]] = DDT[alpha][sbox[x]^sbox[x^alpha]] - 2
				if (DDT[alpha][sbox[x]^sbox[x^alpha]] == 2):
					return(0)


# function which is checking, if the condition for coefficients in ANF is not violated
# this function is taken from [1] and described more preciously in the thesis
def addDegreeInformation(x):
	global CTR
	global SUM
	global sbox
	for u in range(1,2**n):
		if(hammingWeight(u) > 2):
			if(vectorOrdering(x,u) == 1):
				CTR[u]+=1
				SUM[u] = SUM[u]^sbox[x]
				if(CTR[u] == 2**hammingWeight(u)):
					if(SUM[u] != 0):
						return(0)
	return(1)


# function which removes sbox[u] from the xor for coefficients of the ANF during the step from "depth" into "depth-1" in the function "nextVal" 
# this function is taken from [1] and described more preciously in the thesis
def removeDegreeInformation(x):
	global CTR
	global SUM
	global sbox
	for u in range(1,2**n):
		if (hammingWeight(u) > 2):
			if(vectorOrdering(x,u) == 1):
				CTR[u]-=1
				SUM[u] = SUM[u]^sbox[x]
				if(CTR[u] == (2**hammingWeight(u)-1)):
					if (SUM[u] != sbox[x]):
						return(0)
	return (1)


# function tries, if the point "x" with value "sbox[x]" does not violate the APN condition
# this function is taken from [1] and described more preciously in the thesis
def addPoint(x):
	if (addDDTInformation(x)):
		return (addDegreeInformation(x))
	return(0)


# function which removes the points "x" with value "sbox[x]" 
# this function is taken from [1] and described more preciously in the thesis
def removePoint(x):
	if (removeDDTInformation(x)):
		removeDegreeInformation(x)


# function which checks if the given list "array" contains "init_value" at any position
# this function is mentioned in [1] but not in pseudo-code
def isComplete(array):
	global init_value
	for i in range(0,len(array)):
		if(array[i] == init_value):
			return(0)
	return(1)


# function which returns the lowest index "i" on  which the given list contains the value "init_value"
# this function is mentioned in [1] but not in pseudo-code
def nextFreePosition(array):
	global init_value
	for i in range(0,len(array)):
		if(array[i] == init_value):
			return(i)


# function which checks, if the run time of the algorithm does not exceed the time limit provided by user (the user input time limit for each run of the algorithm)
def isTimeLimitOver(time):
	global time_limit
	global time_start
	if(float(time-time_start)>time_limit):
		return(1)
	return(0)



# function which removes all duplicates from the input list "source"; the reduced list is the "output"
def removeDuplicates(source):
	output = []
	last_tested = []
	source.sort()
	i = 0
	while (i<len(source)):
		if (last_tested!=source[i]):
			last_tested = source[i]
			output.append(source[i])
			i = i+1
		else:
			i=i+1
	return(output)


# function which randomly selects a function from each EA-equivalence classes and print them into the file "chapter2-output.txt"
def printRandomElementFromEachEAClassIntoFile(source):
	with open('chapter2-output.txt','w') as output_file:
		output_file.write('Representative of each EA-equivalence class:\n')
	for i in range(0,len(source)):
		index = random.randint(0,len(source[i])-1)
		with open('chapter2-output.txt','a') as output_file:
			output_file.write(str(source[i][index]))
			output_file.write('\n')


# function for finding APN function by setting new value to "sbox[x]" for "x" and checking, if the function can still be APN
# this function is taken from [1] and described more preciously in the thesis; altough it has benn modified, because in [1] were mystakes, which made the function  
def nextVal(x):
	global output_sboxes
	global sbox
	global P
	global init_value
	global time_limit
	if (x==2**n):
		output_sboxes.append(sbox.copy())
		return
	for z in range(0,2**n): 
		y = P[x][z]
		sbox[x] = y
		b = addPoint(x)
		if(b == 1):
			nextVal(x+1)		
		removePoint(x)
		sbox[x] = init_value
		if(isTimeLimitOver(time.time())):
			return


# function for computing the inner product of two integeres (which represent the vector from F_2^n)
def inner_product(l,x):
	return Integer(l&x).bits().count(1)%2


# function which checks if two given functions are EA-equivalent or not; the method is base on the fact that Walsh spectrum is an invariant for EA-equivalence
# the function is described more preciously in the thesis
def areFunctionsEAEquivalent(F,G):
	pi_F = ortho_derivative(F)
	pi_G = ortho_derivative(G)
	#TODO - remove the if
	if(len(pi_F)!=len(pi_G)):
		print('rozdilna_delka')
	for u in range (0,2**(n)):	
		for v in range(0,2**(n)):	
			WS_F = 0
			WS_G = 0		
			for x in range(0,len(pi_F)):
				exponent_F = inner_product(u,pi_F[x])^inner_product(v,x)				
				WS_F = WS_F + (-1)**exponent_F				
			for x in range(0,len(pi_G)):
				exponent_G = inner_product(u,pi_G[x])^inner_product(v,x)
				WS_G = WS_G + (-1)**exponent_G
			if(WS_F != WS_G):
				return(0)
	return(1)


# function which sort given functions "functions" into EA-equivalence classes
# the function is described more preciously in the thesis
def putFunctionsIntoEAClasses(functions):
	EA_classes = []
	is_in_some_class = [init_value]*(len(functions))
	while(isComplete(is_in_some_class)==0):
		x = nextFreePosition(is_in_some_class)
		EA_class = []
		EA_class.append(functions[x])
		is_in_some_class[x] = 1
		for i in range(x+1,len(functions)):
			if(is_in_some_class[i]==init_value):
				if(areFunctionsEAEquivalent(functions[x],functions[i])):
					EA_class.append(functions[i])
					is_in_some_class[i] = 1
		EA_classes.append(EA_class.copy())
	return(EA_classes)	






# -------------------------------------------------------------	
# -------------------- ALGORITHM APNSearch --------------------
# -------------------------------------------------------------	


#----------------------------------------------------------------
# 1st PART - dialog with user to get parameters for the algorithm

# dialog with user to get the following information: dimension n, for how long (in seconds) will the tree search run, how many runs will be made
print('This algorithm serves for finding APN functions. At first, choose the following parameters:')

print('1) Enter the the dimension in which will be the tree search made. (integer between 1 and 8).')
n = int(input())

print('2) Enter for how many seconds will one run lasts at most (integer).')
time_limit = int(input())

print('3) Enter how many runs will be made (each will take ',time_limit,' seconds).')
number_of_runs = int(input())


#------------------------------------------------------------------------
# 2nd PART - setting up variables which are used for the whole algorithm

output_sboxes = []
init_value = -1


#-------------------------------------------------------------------------------------------------------------------------------------------------------
# 3rd PART - algorithm itself; at first the algorithm set up variables used in the i-th run and then the recursive tree search is executed by nextVal(0)

for i in range(0,number_of_runs):
	time_start = time.time()
	sbox = generateSbox(n)
	P = generateP(n)
	DDT = generateDDT(n)
	SUM, CTR = generateSUMandCTR(n)
	nextVal(0)


#---------------------------------------------------------------------------------------------------------------------------------------------
# 4th PART - removing duplicates from the array output_sboxes; this needed to be done only in case when the number of runs is greater than 1;
# that is because in each run, the search is defined by P, which contains only permutations

time_APN_found = time.time()
print("Algorithm found ", len(output_sboxes) ," APN functions for n=",n,'.\n')
if (number_of_runs == 1):
	print('There was only one run, therefore there cannot be any duplicates. The run took:', time_APN_found - time_start ,'seconds.')
	output_sboxes_reduced = output_sboxes
else:
	print('Since the array P is random generated, it is possible, that some APN functions appeared multiple times. Therefore we remove duplicates.\n')
	output_sboxes_reduced = removeDuplicates(output_sboxes)
	time_sort_done = time.time()	
	print('After removing duplicates, we have ', len(output_sboxes_reduced), 'unique APN fucntions. It took:', time_sort_done-time_APN_found ,'seconds. \n')


#-------------------------------------------------------------------------------------------------------------------------
# 5th PART - sort founded APN functions into EA-equivalence classes and print one random element from each class into .txt file 
time_ea_sort_starts = time.time()
print("Now the algorithm starts to sort these functions into EA-equivalence classes.\n")
EA_classes = putFunctionsIntoEAClasses(output_sboxes_reduced)
printRandomElementFromEachEAClassIntoFile(EA_classes)
time_ea_sort_done = time.time()
print("All APN functions founded in the run of algorithm were sort into ", len(EA_classes) ," EA-equivalence classes. It took:", time_ea_sort_done-time_ea_sort_starts, "seconds.")



















