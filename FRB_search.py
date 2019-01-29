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

f, axarr = plt.subplots(4,1,sharex=True) #this '2' is the number of vertical panels
f.subplots_adjust(hspace=0)

f, ax = plt.subplots(4,1,sharex=True) #this '2' is the number of vertical panels
f.subplots_adjust(hspace=0)

for i in range(1,5):
	filename1 = "/home/user/Ajay/tifr/bunchfile/AS1A02_005T01_9000000948cztM0_level2_bunch.fits"
	evt_file1 = fits.open(filename1)
	time_evt1 = evt_file1[i].data['TIME']
	bunch_numevt = evt_file1[i].data['NumEvent']
	
	det1 = evt_file1[i].data['DetId1']
	det2 = evt_file1[i].data['DetId2']
	det3 = evt_file1[i].data['DetId3']
	det4 = evt_file1[i].data['DetId4']

	
	
	
	
	#~ filename2 = "/home/cztworks/ajay/CZTI_DATA/FRB_search/FRB180311/AS1G08_077T01_9000001970_13248cztM0_level2_quad_clean.evt"
	#~ evt_file2 = fits.open(filename2)
	#~ time_evt2 = evt_file2[i].data['TIME']
	#~ veto_time  = evt_file2[5].data['Time']
	#~ veto_counts = evt_file2[5].data['VetoSpec']
	#~ veto_quadid = evt_file2[5].data['QuadID']
	#~ veto_counts	= veto_counts[ :, 24:107]
	#~ veto_counts = veto_counts.transpose()
	#~ index_quad = np.where(veto_quadid==i-1)[0]
	#~ veto_counts = veto_counts[:, index_quad]
	#~ veto_time = veto_time[index_quad]
	
	
	time_del = 100.0

	
	#~ FRB160608
	#~ trig_time = 203053983.088 - 2.0

	
	
	#~ FRB170107
	trig_time = 221515511.000-2.0
	
	#~ FRB180725
	#~ trig_time = 270237586.115-2.0
	#~ UT_0 = 270237586.115-2.0 - time_del
	#~ UT_last = 270237586.115-2.0+ time_del

	UT_0      = trig_time - time_del
	UT_last   = trig_time + time_del
	#comment the following if veto not required
	#~ ind_time		=	np.where( (  UT_0-200.0< veto_time) & ( UT_last+200.0> veto_time ) )[0]
	#~ veto_time		=	veto_time[ind_time]
	#~ veto_counts		=	veto_counts[:, ind_time]
	#~ ax[i-1].plot( veto_time, np.sum( veto_counts, axis = 0 ) )
	#~ ax[i-1].axvline(trig_time,color='red')
	#~ plt.show()
	
	
	#~ print veto_time 
	#~ print veto_counts
	#~ UT_0 = int(time_evt1[0])
	#~ UT_last = int(time_evt1[len(time_evt1)-1])
	


	bindata = np.arange(UT_0,UT_last,0.01)
	
	index_time = [ (time_evt1 > UT_0) & (time_evt1 < UT_last) ][0]
	bunch_numevt_temp = bunch_numevt[index_time]
	time_evt1_temp = time_evt1[index_time]
	det1=det1[index_time]
	det2=det2[index_time]
	det3=det3[index_time]
	det4=det4[index_time]
	#~ print time_evt1_temp
	
	time_temp = np.array([])
	for j in range(0,len(time_evt1_temp)):
		#~ if trig_time-1.55 <time_evt1_temp[j]< trig_time-1.40:
		if trig_time-3.0 <time_evt1_temp[j]< trig_time+3.0:
			print (str(time_evt1_temp[j])+"\t"+str(bunch_numevt_temp[j])+"\t"+str(det1[j])+"\t"+str(det2[j])+"\t"+str(det3[j])+"\t"+str(det4[j])+"\n")
		for k in range(0,bunch_numevt_temp[j]):
			
			time_temp = np.append(time_temp, time_evt1_temp[j])
	

	#~ print time_temp 
	
	lc_main , bin_e1 = np.histogram(time_temp, bins=bindata)
	
	#~ UT_0 = int(time_evt2[0])
	#~ UT_last = int(time_evt1[len(time_evt2)-1])
	#~ bindata = np.arange(UT_0,UT_last,0.001)
	
	
	#~ index_time = [ (time_evt2 > trig_time-5.0) & (time_evt2 < trig_time+5.0) ][0]
	#~ time_evt2_temp = time_evt2[index_time]
	
	bindata = np.arange(UT_0,UT_last,0.001)
	#~ lc_clean, bin_e2 = np.histogram(time_evt2_temp, bins=bindata)
	bin_e1 = np.delete(bin_e1 , len(bin_e1)-1)
	#~ bin_e2 = np.delete(bin_e2 , len(bin_e2)-1)
	
	
	
	
	
	
	
	plt.rc('axes', linewidth=2)
	plt.rc('font', family='serif', serif='cm10')
	#~ plt.rc('text', usetex=True)
	#~ plt.rcParams['text.latex.preamble'] = [r'\boldmath']
	
	P = np.pi # Dear old pi!
	padding     = 3 # The padding of the axes labels.
	size_font = 9 # The fontsize in the images.
	

	
	
	
	
	
	axarr[i-1].plot(bin_e1, lc_main,  label = "Quadrant "+str(i-1))
	
	#~ plt.plot(bin_e2, lc_clean)
	#~ axarr[i-1].axvline(trig_time+1.36,color='red')
	
	axarr[i-1].axvline(trig_time,color='red')
	#~ axarr[i-1].axvline(trig_time-4.0,color='red')
	
	#~ axarr[i-1].axvline(trig_time-1.55,color='green')
	
	axarr[i-1].axvspan(trig_time-3.0,trig_time+3.0 ,color='orange')
	#~ plt.axvspan(trig_time-1.0, trig_time+1.0, alpha=0.5, color='green')
	
	axarr[i-1].legend()
	
	#~ axarr[i-1].locator_params(axis='y', nbins=3)

axarr[1].set_ylabel('Counts', fontsize = size_font - 2 )
axarr[3].set_xlabel('TIME (s)', fontsize = size_font - 2 )
#~ plt.tight_layout()
plt.show()


	
