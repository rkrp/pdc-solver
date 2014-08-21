#!/usr/bin/env python

from pprint import pprint
from matplotlib import pyplot as plt
from scipy.signal import argrelextrema
from scipy import interpolate
from math import factorial
import numpy as np
import csv

#Magic nums for field numbers
_PV_ = 2
_SP_ = 3
_PB_ = 5
_IT_ = 6
_DT_ = 7

class Dataset:
	def __init__(self):
		self.ss = 0
		self.ss_time = 0
		self.pv_ = list()
		self.sp_ = list()
		self.x_ = list()
		self.pb = 0
		self.it = 0
		self.dt = 0
		
def csv_import(datasets):
	data = list()

	for dsnum in range(1, datasets + 1):
		objdata = Dataset()
		pv = list()
		sp = list()
		
		with open('Input/data' + str(dsnum)  + '.csv') as fp:
			reader = csv.reader(fp)
			rowskip = True
			initvals = True
			for rec in reader:
				if(rowskip):
					rowskip = False
					continue
				else:
					pv.append(float(rec[_PV_]))
					sp.append(float(rec[_SP_]))
					if(initvals):
						objdata.pb = float(rec[_PB_])
						objdata.it = float(rec[_IT_])
						objdata.dt = float(rec[_DT_])
					
		#Zero adjustment with SP
		x = zero_adjustment(sp)

		objdata.x_ = x
		objdata.sp_ = sp 
		objdata.pv_ = pv
		data.append(objdata)
	return data
	
def zero_adjustment(sp):
	length = len(sp)
	breakpoint = 0
	for i in range(length):
		if(sp[i] != sp[i + 1]):
			breakpoint = i
			break
	return range(-breakpoint, length - breakpoint)

def cmp_y(a, b, tol=1e-1):
	return (abs(a - b)) <= tol

def steady_state(sp):
	#To be worked on later
	return max(sp)

def get_all_divs(x_, pv_, ss):
	divs = list()
	sameband = False
	for x, pv in zip(x_, pv_):
		if(cmp_y(pv, ss)):
			if not(sameband):
				divs.append(x)
				sameband = True
		else:
			sameband = False
	return divs

def steady_state_time(divs):
	return divs[-1]

def fluctuations(pv, ss, divs):
	flucs = list()
	op = 'max'
	i = 0
	divvalue = float(len(pv)) / float(divs[-1])
	
	while i < len(divs) - 1:
		try:
			start = divs[i] * divvalue
			end = divs[i+1] * divvalue
			
			if op == 'max':
				flucs.append(max(pv[start:end]) - ss)
				op = 'min'
			else:
				flucs.append(min(pv[start:end]) - ss)
				op = 'max'
			i += 1
		except IndexError:
			return flucs
	return flucs
	
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
	# Implementation of Savitzky-Golay Algorithm for smoothing the curve
	# without which determination of maxima and minima was a nightmare
	# Credits to scipy
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window+1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
    
def main():
	data = csv_import(5)
	splot = 320
	i = 0
	
	for dt in data:
		print 'Parsing Dataset', i + 1
		
		#~ f=interpolate.UnivariateSpline(data[i]['X'], data[i]['PV'])
		#~ g=interpolate.interp1d(data[i]['X'], data[i]['PV'])
		
		#~ pvsmooth = savitzky_golay(np.array(dt.pv_), 31, 3)
		pvsmooth = dt.pv_
		print len(pvsmooth)
		pvfunc = interpolate.interp1d(dt.x_, pvsmooth)

		xnew = np.arange(dt.x_[0], dt.x_[-1], 1e-1)
		pvinterp = pvfunc(xnew)
		
		ss = steady_state(dt.sp_)
		#~ divs = get_all_divs(data[i]['X'], pvsmooth, ss)
		#~ flucs = fluctuations(pvsmooth, ss, divs)
		divs = get_all_divs(xnew, pvinterp, ss)
		flucs = fluctuations(pvinterp, ss, divs)
		print 'Divs:', divs
		print 'Flucs:', flucs
		print ''
		
		splot += 1
		
		plt.subplot(splot)
		#~ plt.plot(data[i]['X'], pvsmooth, 'g-')
		plt.xlabel('Time')
		plt.ylabel('SP & PV')

		textdata = '$PB$ = ' + str(dt.pb) + '\n$IT$ = ' + str(dt.it) + '\n$DT$ = ' + str(dt.dt)
		plt.text(len(dt.x_) * 0.7, ss * 0.9, textdata)
		plt.plot(xnew, pvinterp, 'g-')
		plt.plot(dt.x_, dt.sp_, 'b-')
		plt.grid()
		

		#~ plt.plot(data[i]['X'], data[i]['PV'], data[i]['X'], data[i]['SP'])
	plt.show()
	

if __name__ == '__main__':
	main()
