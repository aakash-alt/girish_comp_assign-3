# constant function, f(x)=1


import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.ticker as ticker

def fx(x):
	return 1
	
numpts=1024 # no. of sampling points 
xmin,xmax=np.array([-1000,1000]) # interval in x-space
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

idx=np.argsort(freq)

fig=plt.figure()
ax1=fig.add_subplot(121)
ax2=fig.add_subplot(122)
ax1.plot(xarr,yarr,'y')
ax1.set_xlabel('$x$')
ax1.set_ylabel('$f(x)$')
ax1.set_title('constant function, $f(x)$=1',fontsize=8)
ax2.plot(freq[idx],yfft[idx].real,'k')
ax2.set_xlabel('$k$')
ax2.set_ylabel('$f(k)$')
ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.3f'))
ax2.yaxis.tick_right()
ax2.set_title('FT with 1024 sampling points b/w [-1000,1000]',fontsize=8)
plt.savefig('q#6.png',dpi=500)
plt.show()
