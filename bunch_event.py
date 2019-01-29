from __future__ import division
import os
import numpy as np
import numpy.ma as ma
from astropy.io import fits,ascii
import math
import sys
from matplotlib import pyplot as plt
import scipy.stats.distributions as ps
import collections


f, axarr = plt.subplots(2,1,sharex=True) #this '2' is the number of vertical panels
f.subplots_adjust(hspace=0)

#~ f, ax = plt.subplots(4,1,sharex=True) #this '2' is the number of vertical panels
#~ f.subplots_adjust(hspace=0)


for Q_n in range(1,5):
	filename1 = "/home/user/Ajay/tifr/bunchfile/FRB_search/FRB160608/AS1G05_248T01_9000000486_03761cztM0_level2_bunch.fits"
	evt_file1 = fits.open(filename1)
	
	bunch_time   	= evt_file1[Q_n].data['TIME']
	bunch_numevt 	= evt_file1[Q_n].data['NumEvent']
	bunch_det1 		= evt_file1[Q_n].data['DetId1']
	bunch_det2 		= evt_file1[Q_n].data['DetId2']
	bunch_det3 		= evt_file1[Q_n].data['DetId3']
	bunch_det4 		= evt_file1[Q_n].data['DetId4']
	buntime_dfs 	= evt_file1[Q_n].data['Time_dfs']
	buntime_dsl 	= evt_file1[Q_n].data['Time_dsl']
	
	index_bunch_heavy = [bunch_numevt==63][0]
	
	bunch_time   	= bunch_time[index_bunch_heavy]
	bunch_numevt 	= bunch_numevt[index_bunch_heavy]
	bunch_det1 		= bunch_det1[index_bunch_heavy]
	bunch_det2 		= bunch_det2[index_bunch_heavy]
	bunch_det3 		= bunch_det3[index_bunch_heavy]
	bunch_det4 		= bunch_det4[index_bunch_heavy]
	buntime_dfs 	= buntime_dfs[index_bunch_heavy]
	buntime_dsl 	= buntime_dsl[index_bunch_heavy]
	

	#~ event_stamp_time = np.array([])
	#~ event_stamp_det1 = np.array([])
	#~ event_stamp_det2 = np.array([])
	#~ event_stamp_det3 = np.array([])
	#~ event_stamp_det4 = np.array([])
	
	#~ bunch_time_mod   = bunch_time
	limit = len(bunch_time)
	bunch_time_mod   = np.array([0.0]*len(bunch_time))
	bunch_numevt_mod = np.array([0]*len(bunch_time))
	
	bunch_det_mod = np.array([[0 for x in range(16)] for y in range(limit)] )
	

	
	i=0
	bunch_counter=0
	
	
	while i < limit-1:
		
		bunch_time_mod[bunch_counter]	 				= bunch_time[i]
	
		bunch_numevt_mod[bunch_counter] 				= bunch_numevt[i]
		bunch_det_mod[bunch_counter][bunch_det1[i]] 		= bunch_det_mod[bunch_counter][bunch_det1[i]]+1
		bunch_det_mod[bunch_counter][bunch_det2[i]] 		= bunch_det_mod[bunch_counter][bunch_det2[i]]+1
		bunch_det_mod[bunch_counter][bunch_det3[i]] 		= bunch_det_mod[bunch_counter][bunch_det3[i]]+1
		if bunch_numevt[i]>3:
			bunch_det_mod[bunch_counter][bunch_det4[i]] 		= bunch_det_mod[bunch_counter][bunch_det4[i]]+1
		
		for j in range(i+1,limit): 
			if ( ((bunch_time[j]-20.0*float(buntime_dfs[j])/1000000.0)	-	(bunch_time[j-1] + 20.0*float(buntime_dsl[j-1])/1000000.0) )        < 30.0/1000000.0                ):
				
				bunch_numevt_mod[bunch_counter] 				= bunch_numevt_mod[bunch_counter] + bunch_numevt[j]
				bunch_det_mod[bunch_counter][bunch_det1[j]] 	= bunch_det_mod[bunch_counter][bunch_det1[j]]+1
				bunch_det_mod[bunch_counter][bunch_det2[j]] 	= bunch_det_mod[bunch_counter][bunch_det2[j]]+1
				bunch_det_mod[bunch_counter][bunch_det3[j]] 	= bunch_det_mod[bunch_counter][bunch_det3[j]]+1
				if bunch_numevt[j]>3:
					bunch_det_mod[bunch_counter][bunch_det4[j]] 		= bunch_det_mod[bunch_counter][bunch_det4[j]]+1
		
			else: 
				break
		#~ print (str(i)+"\t"+str(limit)+"\n")		
		#~ print (bunch_det_mod[bunch_counter][12])
		i=j
		bunch_counter=bunch_counter+1
		
	index_heavy = [bunch_numevt_mod>100][0]	
	heavy_bunch_time_mod = bunch_time_mod[index_heavy]
	heavy_bunch_numevt_mod = bunch_numevt_mod[index_heavy]
	
	heavy_bunch_time_mod = np.around(heavy_bunch_time_mod, decimals=5, out=None)
	
	#~ bunch_test_time = bunch_time[index_heavy]
	#~ print (len(heavy_bunch_time_mod))
	heavy_bunch_det_mod = bunch_det_mod[index_heavy]
	#~ print( heavy_bunch_det_mod)
	#~ print (heavy_bunch_time_mod[10])
	data = [heavy_bunch_time_mod, heavy_bunch_numevt_mod]
	#~ print (data)
	#~ heavy_bunch_det_mod = np.transpose(heavy_bunch_det_mod)
	data = np.transpose(data)
	#~ print (np.shape(data)                       )
	#~ print (np.shape(heavy_bunch_det_mod)        )
	data = np.concatenate((data, heavy_bunch_det_mod), 1)
	#~ print (data)
	#~ ascii.write([heavy_bunch_time_mod, heavy_bunch_numevt_mod], 'test_'+str(Q_n)+'.txt', names=['TIME', 'NUMBER_OF_EVENTS'])
	ascii.write(data, 'test_'+str(Q_n)+'.txt', overwrite=True)
	
	#~ outfile = open('test_'+str(Q_n)+'.txt', "w")
	
	#~ for k in range(0,len(heavy_bunch_time_mod)):
		#~ outfile.write(str(heavy_bunch_time_mod[k])+"\t"+str(heavy_bunch_numevt_mod[k])+"\n")

#~ infile1  = ascii.read('test_1.txt',header_start = None, delimiter = '\s', data_start = 1)
infile1  = np.loadtxt('test_1.txt', delimiter = ' ', skiprows =1)
time_1 	 = np.array(infile1[:,[0]])
det_Q1	 = np.array(infile1[:,2:18])
infile2  = np.loadtxt('test_2.txt', delimiter = ' ', skiprows =1)
time_2 	 = np.array(infile2[:,[0]])
det_Q2	 = np.array(infile2[:,2:18])
infile3  = np.loadtxt('test_3.txt', delimiter = ' ', skiprows =1)
time_3 	 = np.array(infile3[:,[0]])
det_Q3	 = np.array(infile3[:,2:18])
infile4  = np.loadtxt('test_4.txt', delimiter = ' ', skiprows =1)
time_4 	 = np.array(infile4[:,[0]])
det_Q4	 = np.array(infile4[:,2:18])

#~ print (time_1)

#~ print (np.shape(time_1))
#~ print (np.shape(det_Q1))
#~ print (det_Q1)
#~ exit()

#~ common_index12 = np.isin(time_1, time_2, assume_unique=False, invert=False)
common_index12 = np.in1d(time_1, time_2, assume_unique=False, invert=False)
common_index13 = np.in1d(time_1, time_3, assume_unique=False, invert=False)
common_index14 = np.in1d(time_1, time_4, assume_unique=False, invert=False)
common_index1234	=	 np.multiply(common_index12,np.multiply(common_index13,common_index14))
common_index21 = np.in1d(time_2, time_1, assume_unique=False, invert=False)
common_index23 = np.in1d(time_2, time_3, assume_unique=False, invert=False)
common_index24 = np.in1d(time_2, time_4, assume_unique=False, invert=False)
common_index2134	=	 np.multiply(common_index21,np.multiply(common_index23,common_index24))
common_index31 = np.in1d(time_3, time_1, assume_unique=False, invert=False)
common_index32 = np.in1d(time_3, time_2, assume_unique=False, invert=False)
common_index34 = np.in1d(time_3, time_4, assume_unique=False, invert=False)
common_index3124	=	 np.multiply(common_index31,np.multiply(common_index32,common_index34))
common_index41 = np.in1d(time_4, time_1, assume_unique=False, invert=False)
common_index42 = np.in1d(time_4, time_2, assume_unique=False, invert=False)
common_index43 = np.in1d(time_4, time_3, assume_unique=False, invert=False)
common_index4123	=	 np.multiply(common_index41,np.multiply(common_index42,common_index43))

#~ print (common_index12)
print(np.sum(common_index12))
print(np.sum(common_index13))
print(np.sum(common_index14))
print(np.sum(common_index1234))
#~ print(np.sum(common_index2134))
#~ print(np.sum(common_index3124))
#~ print(np.sum(common_index4123))

time_1_common = time_1[common_index1234]
time_2_common = time_2[common_index2134]
time_3_common = time_3[common_index3124]
time_4_common = time_4[common_index4123]
#~ ascii.write(np.transpose(time_1_common), 'Time1_common.txt', overwrite=True)
det_Q1_common = det_Q1[common_index1234]
det_Q2_common = det_Q2[common_index2134]
det_Q3_common = det_Q3[common_index3124]
det_Q4_common = det_Q4[common_index4123]
det_Q1_non_common = det_Q1[~common_index1234]
det_Q2_non_common = det_Q2[~common_index2134]
det_Q3_non_common = det_Q3[~common_index3124]
det_Q4_non_common = det_Q4[~common_index4123]


counter_det_Q1_common = np.count_nonzero(np.array(det_Q1_common),axis=1)
counter_det_Q2_common = np.count_nonzero(np.array(det_Q2_common),axis=1)
counter_det_Q3_common = np.count_nonzero(np.array(det_Q3_common),axis=1)
counter_det_Q4_common = np.count_nonzero(np.array(det_Q4_common),axis=1)
counter_det_Q1_non_common = np.count_nonzero(np.array(det_Q1_non_common),axis=1)
counter_det_Q2_non_common = np.count_nonzero(np.array(det_Q2_non_common),axis=1)
counter_det_Q3_non_common = np.count_nonzero(np.array(det_Q3_non_common),axis=1)
counter_det_Q4_non_common = np.count_nonzero(np.array(det_Q4_non_common),axis=1)

counter_det_total_common=total = np.array([0.0]*len(counter_det_Q1_common))
for i in range(0,len(time_1_common)):
	if time_1_common[i] == time_2_common[i] == time_3_common[i] == time_4_common[i]:
		counter_det_total_common[i] = counter_det_Q1_common[i] + counter_det_Q2_common[i] + counter_det_Q3_common[i] + counter_det_Q4_common[i]

counter_det_total_non_common = np.append(counter_det_Q1_non_common, np.append(counter_det_Q2_non_common, np.append(counter_det_Q3_non_common,counter_det_Q4_non_common)))
	
print (counter_det_total_common)
print (counter_det_total_non_common)

bindata = np.arange(0,65,1)
non_com	,	bin_e	=	np.histogram(counter_det_total_non_common,bins=bindata)
com 	,	bin_e	=	np.histogram(counter_det_total_common,bins=bindata)
print (non_com)
bin_e = np.delete(bin_e, len(bin_e)-1)

plt.rc('axes', linewidth=2)
plt.rc('font', family='serif', serif='cm10')
#~ plt.rc('text', usetex=True)
#~ plt.rcParams['text.latex.preamble'] = [r'\boldmath']
P = np.pi # Dear old pi!
padding     = 3 # The padding of the axes labels.
size_font = 9 # The fontsize in the images.

#~ plt.plot(bin_e, non_com)
#~ plt.plot(bin_e, com)
axarr[0].plot(bin_e, non_com,  label = "Non Simultaneous ")
axarr[1].plot(bin_e, com,  label = "Simultaneous ")
#~ axarr[i-1].axvline(trig_time+1.36,color='red')
#~ axarr[i-1].axvspan(trig_time-19.3,trig_time-19.0 ,color='orange')
#~ plt.axvspan(trig_time-1.0, trig_time+1.0, alpha=0.5, color='green')
axarr[0].set_ylabel('No of bunches (bunch size >100) ', fontsize = size_font )
axarr[1].set_xlabel('No of detectors triggered per bunch  (bunch_size > 100)', fontsize = size_font )
axarr[0].legend()
axarr[1].legend()
plt.tight_layout()
#~ plt.show()	
plt.savefig("bunch_event_statistics.pdf")



trig_time = 203053983.088 - 2.0
time2frb = np.abs(time_1_common - trig_time)

outfile = open("test.txt","w")
for j in range(0,len(time2frb)):
	outfile.write(str(time2frb[j])+"\t"+ str(det_Q1_common[j])+"\t"+ str(det_Q2_common[j])+"\n")

print (min(time2frb))

	
		
'''
UT_0 	= heavy_bunch_time_mod[0]
UT_last = heavy_bunch_time_mod[len(heavy_bunch_time_mod)-1]
bindata = np.arange(UT_0,UT_last,0.001)
lc_main , bin_e1 = np.histogram(time_temp, bins=bindata)	
bin_e1 = np.delete(bin_e1 , len(bin_e1)-1)

plt.rc('axes', linewidth=2)
plt.rc('font', family='serif', serif='cm10')
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\boldmath']

P = np.pi # Dear old pi!
padding     = 3 # The padding of the axes labels.
size_font = 9 # The fontsize in the images.

axarr[i-1].plot(bin_e1, lc_main,  label = "Quadrant "+str(i-1))
#~ axarr[i-1].axvline(trig_time+1.36,color='red')
#~ axarr[i-1].axvspan(trig_time-19.3,trig_time-19.0 ,color='orange')
#~ plt.axvspan(trig_time-1.0, trig_time+1.0, alpha=0.5, color='green')
axarr[i-1].legend()
	'''
	
'''
axarr[1].set_ylabel('Counts', fontsize = size_font - 2 )
axarr[3].set_xlabel('TIME (s)', fontsize = size_font - 2 )
#~ plt.tight_layout()
plt.show()	
'''	
