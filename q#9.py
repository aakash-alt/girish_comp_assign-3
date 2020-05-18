import numpy as np
import matplotlib.pyplot as plt

def gx(x):
	x=np.array(x)
	if (x>-1 and x<1):
		return 1
	else:
		return 0

def hx(x):
	x=np.array(x)
	if (x>-1 and x<1):
		return 1
	else:
		return 0
def fx(x):
	x=np.array(x)
	if(x>-2 and x<0):
		return x+2
	if(x>0 and x<2):
		return 2-x
	else:
		return 0

xmin,xmax=-5,5
numpts=256
dx=(xmax-xmin)/(numpts-1)

freq=np.fft.fftfreq(numpts,dx)
freq*=2*np.pi

xarr=np.zeros(numpts)
garr=np.zeros(numpts)
harr=np.zeros(numpts)
farr=np.zeros(numpts)

for i in range(0,numpts):
	xarr[i]=xmin+i*dx
	garr[i]=gx(xarr[i])
	harr[i]=hx(xarr[i])

garr=np.fft.fft(garr,norm='ortho')
harr=np.fft.fft(harr,norm='ortho')
harr=harr*garr
harr=np.fft.ifft(harr,norm='ortho')
harr*=dx*np.sqrt(numpts)
Xarr=np.zeros(numpts)
n1=int(numpts/2)
for i in range(0,n1):
	Xarr[i]=xmin+(n1-i-1)*dx
	farr[i]=fx(Xarr[i])
for i in range(n1,numpts):
	Xarr[i]=xmax+(n1-i)*dx
	farr[i]=fx(Xarr[i])

idx=np.argsort(Xarr)

plt.plot(Xarr[idx],harr[idx].real,'k',label='numeric-$f(x)$')
plt.plot(Xarr[idx],farr[idx],'y',label='analytic-$f(x)$')
plt.legend(shadow=True,fontsize=8)
plt.xlabel('$x$')
plt.ylabel('$f(x)$')
plt.title('Convolution of $g(x)$ & $h(x)$ using 256 sample points in range [-5,5]',fontsize=8)
plt.savefig('q#9.png',dpi=500)
plt.show()

