import numpy as np
import matplotlib.pyplot as plt

data=np.loadtxt('noise.txt')


# plotting the measurement
plt.plot(data,'k')
plt.title('noise',fontsize=8)
plt.savefig('noise.png',dpi=500)
plt.show()


# plotting the DFT
numpts=len(data)
freq=np.fft.fftfreq(numpts)
dft=np.fft.fft(data,norm='ortho')
dft=abs(dft )
idx=np.argsort(freq)
plt.plot(freq[idx],dft[idx],'k')
plt.xlabel('$k$')
plt.ylabel('DFT')
plt.title('DFT of noise measurement',fontsize=8)
plt.savefig('noise-dft.png',dpi=500)
plt.show()


# plotting the power spectrum (using periodogram)
ps=dft*np.conj(dft)
ps=ps.real/numpts
plt.plot(freq[idx],ps[idx],'k')
plt.xlabel('$k$')
plt.ylabel('power-spectrum')
plt.title('power spectrum using peridogrma',fontsize=8)
plt.savefig('power-spectrum-periodogram.png',dpi=500)
plt.show()


# plotting the binned power spectrum (Bartlett's method)
nbin=10
avpspec=np.zeros(51,dtype=np.complex_) # size of each bin
data=np.array_split(data[0:510],nbin)
for i in range(nbin):
	freq=2*np.pi*np.fft.fftfreq(len(data[i]),d=1) # '1' is the step size
	dft=np.fft.fft(data[i],norm='ortho')
	pspec=dft*np.conj(dft)/len(data[i])
	avpspec=avpspec+pspec
avpspec*=1/nbin
idx2=np.argsort(freq)
plt.plot(freq[idx2],avpspec[idx2].real,'k')
plt.title('power psectrum using Bartlett method',fontsize=8)
plt.xlabel('$k$')
plt.ylabel('power-spectrum')
plt.savefig('power-spectrum-Bartlett.png',dpi=500)
plt.show()
