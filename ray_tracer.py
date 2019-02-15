import lsc_geometry as lsc
import launch_photon as launch
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import numpy as np
import spectrum as spec

def lscParameters():
	l_x = 6.0
	l_y = 2.0
	l_z = 0.4
	lsc_dimensions = np.array([l_x,l_y,l_z])
	lsc_definition = [
	    (0,0,0), (0,l_y,0), (l_x,0,0), (0,0,l_z)
	]
	return lsc_dimensions,lsc_definition

def plotSpectrumCDF(fig):
	excel_file = ".\\spectra\\21_06_2019.xlsx"
	df = spec.loadSpectralData(excel_file)
	ax1 = fig.add_subplot(211)
	cdf = spec.plotCDF(df,ax1)
	plt.title('CDF')
	return df

def runLSCSim(fig,lsc_dimensions,lsc_definition,cdf):
	ax2 = fig.add_subplot(212, projection='3d')
	lsc.plotLSC(lsc_definition,fig,ax2)
	launch_parameters, boundary_pt = launch.launchPhoton(lsc_dimensions,cdf)
	launch_pt = launch_parameters[0:3]
	launch.plotIncidentRay(lsc_dimensions,launch_pt,boundary_pt,ax2)
	launch.printLaunchParameters(launch_parameters)

fig = plt.figure()
lsc_parameters = lscParameters()
spectrum = plotSpectrumCDF(fig)
runLSCSim(fig,lsc_parameters[0],lsc_parameters[1],spec.spectralCDF(spectrum))
plt.show()