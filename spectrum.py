# Load Spectrum Data
import numpy as np
import pandas as pd 
import scipy.constants as const

def getPhotonFlux(row,cols):
	wavelength = row[cols[0]]
	irradiance = row[cols[1]]
	photon_energy = const.Planck*const.c/wavelength
	photon_flux = irradiance/photon_energy
	return photon_flux

def addPhotonFluxCol(df):
	photon_flux = 'Photon Flux (photons/m2/nm)'
	cols = list(df)
	df[photon_flux] = df.apply(lambda row: getPhotonFlux(row,list(df)), axis=1)

def photonFluxCumhist(df):
	cumsum_photon_flux = 'Cumulative Photon Flux (%)'
	photon_flux = list(df)[2]
	total_photon_flux = df[photon_flux].sum()
	df[cumsum_photon_flux] = df[photon_flux].cumsum(axis=0)/total_photon_flux

def photonFluxPDF(df):
	pdf_photon_flux = 'Photon Flux PDF (%)'
	photon_flux = list(df)[2]
	total_photon_flux = df[photon_flux].sum()
	df[pdf_photon_flux] = df[photon_flux]/total_photon_flux

def addColumns(df):
	addPhotonFluxCol(df)
	photonFluxPDF(df)
	photonFluxCumhist(df)

def loadSpectralData(excel_file):
	sheet = "Spectral irradiance"
	col_nums = [0,7]
	df = pd.read_excel(excel_file,sheetname=sheet,usecols=col_nums)
	addColumns(df)
	return df

def spectralCDF(df):
	return np.array([df[list(df)[0]],df[list(df)[3]]])

def plotPDF(df,ax):

	ax.plot(df[list(df)[0]],df[list(df)[3]])

def plotCDF(df,ax):

	ax.plot(df[list(df)[0]],df[list(df)[4]])
