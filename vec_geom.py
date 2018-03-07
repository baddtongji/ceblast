import math
from numpy import inner, array, sqrt, arccos

def dotproduct(v1, v2):
    """
    dot product between two vector
    
    Paramter:
    v1: array of number
    v2: array of number

    Return: 
    number
    """
    return inner(array(v1), array(v2))

def length(v):
    """
    the norm of vector v, aka, |v|
    """
    return sqrt(dotproduct(v, v))

def angle(v1, v2):
    """
    the angle between two vector.

    Paramter:
    v1: array of number
    v2: array of number

    Return: 
    float
    """
    return arccos(dotproduct(v1, v2) / (length(v1) * length(v2)))
  
