import numpy as np
import matplotlib.pyplot as plt

def fx(x):
	return np.sin(x)/x
def fk(k):
	if(k>=-1 and 1>=k):
		a=1
	else:
		a= 0
	return a*np.sqrt(np.pi/2)
	
numpts=256 # no. of sampling points 
xmin,xmax=np.array([-50,50]) # interval in x-space
dx=(xmax-xmin)/(numpts-1) # sampling rate

xarr=np.zeros(numpts)
yarr=np.zeros(numpts)

for i in range (0,numpts):
	xarr[i]=xmin+i*dx
	if(xarr[i]==0):
		yarr[i]=1
	yarr[i]=fx(xarr[i])
yfft=np.fft.fft(yarr,norm='ortho')
freq=np.fft.fftfreq(numpts,dx)
freq=2*np.pi*freq
factor=np.exp(-1j*freq*xmin)
yfft*=dx*np.sqrt(numpts/(2*np.pi))*factor

fkk=np.zeros(numpts)
for i in range(0,numpts):
	fkk[i]=fk(freq[i])

idx=np.argsort(freq)

fig=plt.figure()
ax1=fig.add_subplot(121)
ax2=fig.add_subplot(122)
ax1.plot(xarr,yarr,'r')
ax1.set_xlabel('$x$')
ax1.set_ylabel('$f(x)$')
ax2.plot(freq[idx],yfft[idx].real,'k',label='numeric')
ax2.plot(freq[idx],fkk[idx],'y',label='analytic')
ax2.set_xlabel('$k$')
ax2.set_ylabel('$f(k)$')
ax2.legend(fontsize=8,shadow=True)
ax2.yaxis.tick_right()
plt.suptitle('Fourier transform with 256 sampling points b/w [-50,50]',fontsize=8)
plt.savefig('q#1.png',dpi=500)
plt.show()
