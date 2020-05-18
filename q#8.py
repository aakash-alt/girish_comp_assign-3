import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# defining f(x,y)
def fxy(x,y):
	return np.exp(-(x**2+y**2))

# defining FT of f(x,y)
def FTf(kx,ky):
	return (1/2)*np.exp(-(kx**2+ky**2)/4)
numpts=256
xmin,xmax=np.array([-50,50])
ymin,ymax=np.array([-50,50])
dx=(xmax-xmin)/(numpts-1)
dy=(ymax-ymin)/(numpts-1)

xarr=np.zeros(numpts,dtype=np.complex_)
yarr=np.zeros(numpts,dtype=np.complex_)
for i in range(0,numpts):
	xarr[i]=xmin+i*dx
	yarr[i]=ymin+i*dy

# sampling data
sampledata=np.ones((numpts,numpts),dtype=np.complex_)
for i in range(0,numpts):
	for j in range(0,numpts):
		sampledata[i,j]=fxy(xarr[i],yarr[j])
# DFT of sampled data (sampledata), and storing in itself
sampledata=np.fft.fft2(sampledata,norm='ortho')
kx=np.fft.fftfreq(numpts,dx)
ky=np.fft.fftfreq(numpts,dy)
kx*=2*np.pi
ky*=2*np.pi
for i in range(0,numpts):
	for j in range(0,numpts):
		sampledata[i][j]=dx*dy*(numpts/(2*np.pi))*(np.exp(-1j*kx[i]*xmin+ -1j*ky[j]*ymin))*sampledata[i][j]
		
# creating arrays for plotting
x=np.zeros(numpts*numpts)
y=np.zeros(numpts*numpts)
fkN=np.zeros(numpts*numpts)
fkA=np.zeros(numpts*numpts)

for i in range(0,numpts):
	for j in range(0,numpts):
		x[j+i*numpts]=kx[i]
		y[j+i*numpts]=ky[j]
		fkN[i+j*numpts]=sampledata[i][j].real
		fkA[i+j*numpts]=FTf(kx[i],ky[j])
		

fig=plt.figure()
ax1=fig.add_subplot(121,projection='3d')
ax2=fig.add_subplot(122,projection='3d')
ax1.plot3D(x,y,fkA,'.y',label='analytic')
ax1.plot3D(x,y,fkN,'.c',label='numeric')
ax2.plot3D(x,y,abs(fkN-fkA),'.c',label='numeric')
ax1.legend(fontsize=5,loc='center left')
ax1.set_xlabel('$kx$')
ax1.set_ylabel('$ky$')
ax2.set_xlabel('$kx$')
ax2.set_ylabel('$ky$')
ax1.set_title('FT of 2D-Gassian with 512 sample points in range [-50,50]',fontsize=6)
ax2.set_title('error plot',fontsize=6)
plt.savefig('q#8.png',dpi=500)
plt.show()

