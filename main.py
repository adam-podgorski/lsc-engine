from datetime import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

date = dt.now()
P = 102.2e3
longitude = 150 
latitude = -34
coordinates = np.array([longitude,latitude])
zenith = 4 # Arbitrary for now.

# From the Simple Bird Spectral Model
class Spectrum():

	def __init__(self,date,coordinates,zenith):
		# Time
		self.hour = date.hour
		self.day = date.timetuple().tm_yday
		# Position
		self.longitude = coordinates[0]
		self.latitude = coordinates[1]
		# Sun position
		self.zenith = np.deg2rad(zenith) # Radian output
		# Extraterrestrial Solar Irradiance
		self.atmospheric_data = self.readTextFiles()
		self.wavelength = self.atmospheric_data['lambda'].values #microns
		# Sec Rayleigh Scattering	
		# Relative air mass
		self.P0 = 101.325e3
		self.M = 1/(np.cos(np.deg2rad(zenith))+0.15*np.power(93.885-self.zenith,-1.253))
		# Pressure corrected air mass
		self.M_prime = self.M*P/self.P0
		# Sec Aerosol Scattering
		self.beta_aerosol = 0.27
		# Sec H2O Scattering
		self.w_vert = 1.42 # cm
		# Sec O3 Scattering
		self.h_o = 22e3
		self.o3 = self.o3()

	def beamIrrad(self):
		h0 = self.atmospheric_data['etsi']
		dES = self.earthSunDistance()
		tr = self.rayleighT(self.M_prime)
		ta = self.aerosolT(self.M)
		tw = self.h2oT(self.M)
		to = self.o3T(self.m_o())
		tu = self.unifT(M_prime)
		return h0*dES*tr*ta*tw*to*tu*np.cos(self.zenith)

	def readTextFiles(self):
		df = pd.read_csv('./spectra/etsi.txt', sep="\t", header=None)
		df.columns = ['lambda','etsi']
		df['a_w'] = pd.read_csv('./spectra/a_w.txt',sep="\t",header=None).iloc[:,1]
		df['a_o'] = pd.read_csv('./spectra/a_o.txt',sep="\t",header=None).iloc[:,1]
		df['a_u'] = pd.read_csv('./spectra/a_u.txt',sep="\t",header=None).iloc[:,1]
		return df

	def earthSunDistance(self):
		psi = 2*np.pi*(self.day-1)/365
		d1 = 1.00011+0.034221*np.cos(psi)+0.00128*np.sin(psi)
		d2 = 0.000719*np.cos(2*psi)+0.000077*np.sin(2*psi)
		return d1+d2

	def rayleighT(self,M):
		lambda_squared = self.wavelength**2
		numer = -1*M
		denom = lambda_squared**2*(115.6406-1.3366/lambda_squared)
		return np.exp(numer/denom)

	def aerosolT(self,M):
		turbidity = self.turbidity()
		return np.exp(-1*turbidity*M)

	def turbidity(self):
		x = self.wavelength
		alpha1 = 1.0274
		alpha2 = 1.2060
		alpha = np.piecewise(x,[x<0.5,x>=0.5],[alpha1,alpha2])
		turbidity = self.beta_aerosol*np.power(x,-1*alpha)
		return turbidity

	def h2oT(self,M):
		a_w = self.atmospheric_data['a_w'].values
		exp1 = -0.2385*a_w*self.w_vert*self.M
		exp2 = np.power(1+20.07*a_w*self.w_vert*self.M,0.45)
		return np.exp(exp1/exp2)

	def m_o(self):
		h_o = self.h_o
		numerator = 1+h_o/6370
		denominator = np.sqrt(np.cos(self.zenith)**2+2*h_o/6370)
		return numerator/denominator

	def o3(self):
		varphi = self.latitude
		lambda_L = self.longitude
		beta = 1.28
		a = 150
		c = 40
		d = 0.9856
		e = self.day
		f = -30
		g = 20
		h = 3
		i = 20
		j = 235
		sin1 = np.sin(np.deg2rad(d*(e+f)))
		sin2 = np.sin(np.deg2rad(h*(lambda_L+i)))
		sin3 = np.sin(np.deg2rad(beta*varphi))
		return sin1
		#return j+(a+c*sin1+g*sin2)*sin3**2

	def o3T(self,M):
		a_o = self.atmospheric_data['a_o'].values
		return np.exp(-1*a_o*self.o3*M)

	def unifT(self,M):
		a_u = self.atmospheric_data['a_u'].values
		exp1 = -1.41*a_u*M
		exp2 = pow(1+118.93*a_u*M,0.45)
		return np.exp(exp1/exp2)

	def scatteredIrrad(self):
		i_r = self.rayleighIrrad()
		i_a = self.aerosolIrrad()
		i_g = self.groundIrrad()
		return i_r+i_a+i_g

	def rayleighIrrad(self):
		h0 = self.atmospheric_data['etsi']
		dES = self.earthSunDistance()
		tr = self.rayleighT(self.M_prime)
		taa = self.aerosolTa(self.M)
		tw = self.h2oT(self.M)
		to = self.o3T(self.m_o())
		tu = self.unifT(M_prime)
		return h0*dES*cos(self.zenith)*to*tu*tw*taa*(1-np.power(tr,0.95))*0.5

	def aerosolIrrad(self):
		h0 = self.atmospheric_data['etsi']
		dES = self.earthSunDistance()
		fs = self.f_s()
		tr = self.rayleighT(self.M_prime)
		taa = self.aerosolTa(self.M)
		tas = self.aerosolTs(self.M)
		tw = self.h2oT(self.M)
		to = self.o3T(self.m_o())
		tu = self.unifT(M_prime)
		part1 = h0*dES*cos(self.zenith)*to*tu*tw*taa
		part2 = np.power(tr,1.5)*(1-tas)(fs)
		return part1*part2

	def groundIrrad(self,rg=0.2):
		cosZ = np.cos(self.zenith)
		i_dir = self.beamIrrad()
		i_rayl = self.rayleighIrrad()
		i_aero = self.aerosolIrrad()
		rs = self.rs()
		numer = rs*rg*(i_dir*cosZ+i_rayl+i_aero)
		denom = 1 - rs*rg
		return numer/denom

	def rs(self):
		M=1.8
		fs = self.f_prime_s()
		tr = self.rayleighT(M)
		taa = self.aerosolTa(M)
		tas = self.aerosolTs(M)
		tw = self.h2oT(M)
		to = self.o3T(M)
		part1 = to*tw*taa
		part2 = 0.5*(1-tr)+(1-f_s)*tr*(1-tas)
		return part1*part2


	def aerosolTs(self,M):
		omega = self.omega()
		tau_a = self.turbidity()
		return np.exp(-1*omega*tau_a*M)

	def aerosolTa(self,M):
		omega = self.omega()
		tau_a = self.turbidity()
		return np.exp(-1*(1-omega)*tau_a*M)

	def f_s(self):
		afs = self.afs()
		bfs = self.bfs()
		cosZ = np.cos(self.zenith)
		return 1-0.5*np.exp((afs+bfs*cosZ)*cosZ)

	def f_prime_s(self):
		afs = self.afs()
		bfs = self.bfs()
		return 1-0.5*np.exp((afs+bfs/1.8)/1.8)

	def afs(self):
		alg = self.alg()
		return alg*(1.459+alg*(0.1595+alg*0.4129))

	def bfs(self):
		alg = self.alg()
		return alg*(0.0783+alg*(-0.3824-alg*0.5874))

	def alg(self):
		cos_theta
		return np.log(1-cos_theta)

	def omega(self):
		x = self.wavelength
		omega_04 = 0.945
		omega_prime = 0.095
		return omega_04*np.exp(-omega_prime*np.log(x/0.4)**2)

	def scatteredIrradCorr(self):
		x = self.wavelength
		C = np.piecewise(x,[x<=0.45,x>0.45],[lambda x: np.power(x+0.55,1.8),1])
		return C*self.scatteredIrrad()

def plot_spectra(date,coordinates,zenith):
	spectrum = Spectrum(date,coordinates,zenith)
	plt.plot(spectrum.wavelength,spectrum.atmospheric_data['etsi'])
	plt.plot(spectrum.wavelength,spectrum.beamIrrad())
	plt.show()

def plot_absorption_coefficients(date,coordinates,zenith):
	spectrum = Spectrum(date,coordinates,zenith)
	a_w = spectrum.atmospheric_data['a_w']
	a_o = spectrum.atmospheric_data['a_o']
	a_u = spectrum.atmospheric_data['a_u']
	#plt.plot(spectrum.wavelength,a_w)
	plt.plot(spectrum.wavelength,a_o)
	#plt.plot(spectrum.wavelength,a_u)
	plt.show()

def plot_transmissions(date,coordinates,zenith):
	spectrum = Spectrum(date,coordinates,zenith)
	#plt.plot(spectrum.wavelength,spectrum.rayleighT())
	#plt.plot(spectrum.wavelength,spectrum.aerosolT())
	#plt.plot(spectrum.wavelength,spectrum.h2oT())
	plt.plot(spectrum.wavelength,spectrum.o3T())
	#plt.plot(spectrum.wavelength,spectrum.unifT())
	plt.show()

#plot_absorption_coefficients()
#plot_transmissions()
#plot_spectra()

