import numpy as np
import matplotlib.pyplot as plt
import time
import timeit

nmin,nmax,dn=4,100,4
numpts=np.arange(nmin,nmax+dn,dn,dtype=np.int_) 
t_dc=np.zeros(int(nmax/nmin)) # array storing time taken using direct computation
t_fft=np.zeros(int(nmax/nmin)) # array stroing time taken using Numpy library
for n in numpts:
	arr=np.arange(n,dtype=np.complex_)
	farr=np.zeros(n,dtype=np.complex_) # array for DFT using direct computation
	
	# computing time taken using direct cmputation
	t0=time.time()
	
	for i in range(0,n):
		for k in range(0,n):
			farr[i]+=arr[k]*np.exp(-1j*2*np.pi*i*k/n)
		farr[i]*=1/np.sqrt(n)
	t1=time.time()
	t_dc[int(n/nmin-1)]=t1-t0
	
	# computing time taken using Numpy library
	t2=time.time()
	fftarr=np.fft.fft(arr,norm='ortho') # array for DFT using Numpy library
	t3=time.time()
	t_fft[int(n/nmin-1)]=t3-t2

 # making plot for time vs numpts
plt.plot(numpts,t_dc,'k',label='direc computation')
plt.plot(numpts,t_fft,'r',label='Numpy library')
plt.yscale('log')
plt.xlabel('$n$')
plt.ylabel('$t$')
plt.legend(shadow=True,fontsize=8)
plt.title('Comparing time taken for DFT using direct computation and Numpy libraries',fontsize=8)
plt.savefig('q#5.png',dpi=500)
plt.show()
