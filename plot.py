import numpy as np
import matplotlib.pyplot as plt

# sinc function
def fx(x):
	return np.sin(x)/x
	
# Fourier transform of sinc function
def fk(k):
	if(k>=-1 and 1>=k):
		a=1
	else:
		a= 0
	return a*np.sqrt(np.pi/2)
fft=np.loadtxt('q#1.dat')
fftw=np.loadtxt('q#2.dat')
freq=fft[:,0]
yfft=fft[:,1]
yfftw=fftw[:,1]

numpts=256 # no. of sampling points 
xmin,xmax=np.array([-50,50]) # interval in x-space
dx=(xmax-xmin)/(numpts-1) # sampling rate

fkk=np.zeros(numpts)
for i in range(0,numpts):
	fkk[i]=fk(freq[i])
xarr=np.zeros(numpts)
yarr=np.zeros(numpts)

for i in range (0,numpts):
	xarr[i]=xmin+i*dx
	if(xarr[i]==0):
		yarr[i]=1
	yarr[i]=fx(xarr[i])

idx=np.argsort(freq)

fig=plt.figure()
ax1=fig.add_subplot(121)
ax2=fig.add_subplot(122)
ax1.plot(xarr,yarr,'y')
ax1.set_xlabel('$x$')
ax1.set_ylabel('$f(x)$')
ax1.set_title('sinc function',fontsize=8)
ax2.plot(freq[idx],yfft[idx],'c',label='numeric-python')
ax2.plot(freq[idx],yfftw[idx],'k',label='numeric-fftw')
ax2.plot(freq[idx],fkk[idx],'y',label='analytic')
ax2.set_xlabel('$k$')
ax2.set_ylabel('$f(k)$')
ax2.legend(fontsize=5,shadow=True,loc='best')
ax2.yaxis.tick_right()
ax2.set_title('Fourier transform using Numpy & FFTW library',fontsize=8)
plt.savefig('q#2.png',dpi=500)
plt.show()
