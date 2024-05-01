# Attachment A5

# All founded functions using algorithm "chapter2-algorithm1-A3.py" for n=2 from the subsection 2.3.2 of the thesis.

"""
The algorithm run with parametr for dimension equal to: 2
The number of founded functions: 192
The found functions are:
"""

def all_found_functions():
    return [[0, 1, 0, 2],
[0, 1, 0, 0],
[0, 1, 0, 3],
[0, 1, 1, 2],
[0, 1, 1, 1],
[0, 1, 1, 3],
[0, 1, 2, 2],
[0, 1, 2, 1],
[0, 1, 2, 0],
[0, 1, 3, 1],
[0, 1, 3, 0],
[0, 1, 3, 3],
[0, 3, 0, 2],
[0, 3, 0, 1],
[0, 3, 0, 0],
[0, 3, 1, 1],
[0, 3, 1, 0],
[0, 3, 1, 3],
[0, 3, 2, 2],
[0, 3, 2, 0],
[0, 3, 2, 3],
[0, 3, 3, 2],
[0, 3, 3, 1],
[0, 3, 3, 3],
[0, 2, 0, 1],
[0, 2, 0, 0],
[0, 2, 0, 3],
[0, 2, 1, 2],
[0, 2, 1, 1],
[0, 2, 1, 0],
[0, 2, 2, 2],
[0, 2, 2, 1],
[0, 2, 2, 3],
[0, 2, 3, 2],
[0, 2, 3, 0],
[0, 2, 3, 3],
[0, 0, 0, 2],
[0, 0, 0, 1],
[0, 0, 0, 3],
[0, 0, 1, 2],
[0, 0, 1, 0],
[0, 0, 1, 3],
[0, 0, 2, 1],
[0, 0, 2, 0],
[0, 0, 2, 3],
[0, 0, 3, 2],
[0, 0, 3, 1],
[0, 0, 3, 0],
[3, 1, 0, 1],
[3, 1, 0, 0],
[3, 1, 0, 3],
[3, 1, 1, 2],
[3, 1, 1, 1],
[3, 1, 1, 0],
[3, 1, 2, 2],
[3, 1, 2, 1],
[3, 1, 2, 3],
[3, 1, 3, 2],
[3, 1, 3, 0],
[3, 1, 3, 3],
[3, 3, 0, 2],
[3, 3, 0, 1],
[3, 3, 0, 3],
[3, 3, 1, 2],
[3, 3, 1, 0],
[3, 3, 1, 3],
[3, 3, 2, 1],
[3, 3, 2, 0],
[3, 3, 2, 3],
[3, 3, 3, 2],
[3, 3, 3, 1],
[3, 3, 3, 0],
[3, 2, 0, 2],
[3, 2, 0, 0],
[3, 2, 0, 3],
[3, 2, 1, 2],
[3, 2, 1, 1],
[3, 2, 1, 3],
[3, 2, 2, 2],
[3, 2, 2, 1],
[3, 2, 2, 0],
[3, 2, 3, 1],
[3, 2, 3, 0],
[3, 2, 3, 3],
[3, 0, 0, 2],
[3, 0, 0, 1],
[3, 0, 0, 0],
[3, 0, 1, 1],
[3, 0, 1, 0],
[3, 0, 1, 3],
[3, 0, 2, 2],
[3, 0, 2, 0],
[3, 0, 2, 3],
[3, 0, 3, 2],
[3, 0, 3, 1],
[3, 0, 3, 3],
[2, 1, 0, 2],
[2, 1, 0, 1],
[2, 1, 0, 0],
[2, 1, 1, 1],
[2, 1, 1, 0],
[2, 1, 1, 3],
[2, 1, 2, 2],
[2, 1, 2, 0],
[2, 1, 2, 3],
[2, 1, 3, 2],
[2, 1, 3, 1],
[2, 1, 3, 3],
[2, 3, 0, 2],
[2, 3, 0, 0],
[2, 3, 0, 3],
[2, 3, 1, 2],
[2, 3, 1, 1],
[2, 3, 1, 3],
[2, 3, 2, 2],
[2, 3, 2, 1],
[2, 3, 2, 0],
[2, 3, 3, 1],
[2, 3, 3, 0],
[2, 3, 3, 3],
[2, 2, 0, 2],
[2, 2, 0, 1],
[2, 2, 0, 3],
[2, 2, 1, 2],
[2, 2, 1, 0],
[2, 2, 1, 3],
[2, 2, 2, 1],
[2, 2, 2, 0],
[2, 2, 2, 3],
[2, 2, 3, 2],
[2, 2, 3, 1],
[2, 2, 3, 0],
[2, 0, 0, 1],
[2, 0, 0, 0],
[2, 0, 0, 3],
[2, 0, 1, 2],
[2, 0, 1, 1],
[2, 0, 1, 0],
[2, 0, 2, 2],
[2, 0, 2, 1],
[2, 0, 2, 3],
[2, 0, 3, 2],
[2, 0, 3, 0],
[2, 0, 3, 3],
[1, 1, 0, 2],
[1, 1, 0, 1],
[1, 1, 0, 3],
[1, 1, 1, 2],
[1, 1, 1, 0],
[1, 1, 1, 3],
[1, 1, 2, 1],
[1, 1, 2, 0],
[1, 1, 2, 3],
[1, 1, 3, 2],
[1, 1, 3, 1],
[1, 1, 3, 0],
[1, 3, 0, 1],
[1, 3, 0, 0],
[1, 3, 0, 3],
[1, 3, 1, 2],
[1, 3, 1, 1],
[1, 3, 1, 0],
[1, 3, 2, 2],
[1, 3, 2, 1],
[1, 3, 2, 3],
[1, 3, 3, 2],
[1, 3, 3, 0],
[1, 3, 3, 3],
[1, 2, 0, 2],
[1, 2, 0, 1],
[1, 2, 0, 0],
[1, 2, 1, 1],
[1, 2, 1, 0],
[1, 2, 1, 3],
[1, 2, 2, 2],
[1, 2, 2, 0],
[1, 2, 2, 3],
[1, 2, 3, 2],
[1, 2, 3, 1],
[1, 2, 3, 3],
[1, 0, 0, 2],
[1, 0, 0, 0],
[1, 0, 0, 3],
[1, 0, 1, 2],
[1, 0, 1, 1],
[1, 0, 1, 3],
[1, 0, 2, 2],
[1, 0, 2, 1],
[1, 0, 2, 0],
[1, 0, 3, 1],
[1, 0, 3, 0],
[1, 0, 3, 3]]