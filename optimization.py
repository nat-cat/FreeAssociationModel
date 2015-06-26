#Herein lies a list of optimization functions of various types

from scipy.optimize import fmin as simplex #uses the Nelder-Mead simplex method to minimize function inputs
from scipy.optimize import fmin_tnc as tnc #uses the truncated Newton method and allows us to bound the answer
import numpy
import time
import random
from math import *


function1 = lambda x: (4*x**2+15.6*x+8.2) #first function to minimize: 1D quadratic
guess1 = 0.5 #just some guess as to the answer


function2 = lambda x: -(x[0]-x[0]**2+x[1]-x[1]**2) #second function to minimize: 2D quadratic
guess2 = [0.2,0.7] #just some guess as to the answer


na = 10.0 #number of monomers in polymer A
nb = 100.0 #number of monomers in polymer B

xT = (sqrt(na*nb)-nb)/(na-nb) #critical phi value

chiCrit = ((sqrt(na)+sqrt(nb))**2)/(2*na*nb) #critical chi value


def f(x,chi):  #free energy density at a given phi value
	if x<1e-4: #if phi is basically 0, remove x*log(x) because it's 0
		return (1-x)/nb*log(1-x)+chi*x*(1-x)
	elif x>.9999: #if phi is basically 1, remove (1-x)*log(1-x) because it's 0
		return x/na*log(x)+chi*x*(1-x)
	else:
		return x/na*log(x)+(1-x)/nb*log(1-x)+chi*x*(1-x)

def fP(x,chi):  #derivative of free energy density at a given phi value
	return (1/na)*(log(x)+1)-(1/nb)*(log(1-x)-1)+chi*(1-2*x)

def fT(x,chi): #takes a list of 2 phi values, combined free energy equation with phi_1 and phi_2
	return (xT-x[1])/(x[0]-x[1])*f(x[0],chi)+(xT-x[0])/(x[1]-x[0])*f(x[1],chi)

def fTprime(x,chi): #returns gradient of the combined free energy equation wrt phi_1 and phi_2 in list format
	'''[derivative of f(phi1),derivative of f(phi2)]'''
	return [(xT-x[1])*((x[1]-x[0])*fP(x[0],chi)+f(x[0],chi)-f(x[1],chi))/(x[1]-x[0])**2,
			(xT-x[0])*((x[0]-x[1])*fP(x[1],chi)-f(x[0],chi)+f(x[1],chi))/(x[1]-x[0])**2]

def s(x): #spinodal curve at each value of phi
	return (1/(na*x)+1/(nb*x))/2

def binodalPoint(chiValue): #returns list of phi_1 and phi_2 to get that chi value
	chi = chiValue
	result = tnc(lambda x: fT(x,chi),phiGuess,fprime = lambda x: fTprime(x,chi),bounds=([0,xT],[xT,1]),approx_grad=0)
	return result[0] #[phi_1,phi_2]

phiGuess = [xT-.01,xT+.01] #just some guess as to the values of phi_1 and phi_2 at a given chi
multiplier = 10.0**(1/50.0)
# chiList = numpy.arange(chiCrit,10.0*chiCrit,chiCrit*multiplier)
# chiList = list(chiList)


def binodalCurve(): #creates a list, coords, with 50 sublists consisting of [phi_1,phi_2,chi] for 50 chi values
	coords = []
	chiList = []
	phiList = []
	chi = chiCrit
	dataPoints = 50 #how many datapoints you want
	multiplier = 10.0**(1.0/dataPoints)
	while chi < 10.0*chiCrit:
		result = binodalPoint(chi)
		phi1 = result[0]
		phi2 = result[1]
		coords.append([phi1,chi])
		coords.append([phi2,chi])
		chi *= multiplier
	#sort with low list and high list, flip low list at the end
	coords = sorted(coords,key=lambda l:l[0], reverse=False) #sorts the coordinates by phi in order from smallest phi to largest for graphing
	return coords #coords is a list of coordinates for binodal



print "binodal: ",binodalCurve()
# # print 'for Quadratic 1D:'
# # print simplex(function1,guess1)
# # print '----------------------------------------------------'
# # print 'for Quadratic 2D:'
# # print simplex(function2,guess2)
# # print '----------------------------------------------------'
# print 'for Binodal Point at Chi = 0.1:'
# print binodalPoint(.1)
# print "chiCrit = ",chiCrit
# print "phiCrit = ",xT
# print '----------------------------------------------------'
# print 'for Binodal Curve:'
# # coords,phiList = binodalCurve()
# for coord in coords:
# 	print coord
# print "chiCrit = ",chiCrit
# print "phiCrit = ",xT
# print '----------------------------------------------------'