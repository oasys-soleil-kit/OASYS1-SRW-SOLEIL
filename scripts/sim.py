#!/env/python
# -*- coding: utf-8 -*-
#############################################################################
# Edge Radiation from Dipole
# Raul Freitas & Sergio Lordano
# Based on: SRWLIB Example#1: 
# Modified by Ferenc Borondics
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
import sys
#sys.path.append('/Users/raul.freitas/SRW_Dev/env/work/srw_python/')

from srwlib import *
from uti_plot import *
import os
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, c, e
#import my_package as mypkg

# Convert wavelength [µm] to energy [eV]
def wl2eV(wl):
    return 1E6 * h * c / e / wl

#********************** Auxiliary function to read tabulated 3D Magnetic Field data from ASCII file:
def AuxReadInMagFld3D(filePath, sCom):
	f = open(filePath, 'r')
	f.readline() #1st line: just pass
	global xStart,xStep,xNp,yStart,yStep,yNp,zStart,zStep,zNp
	xStart = float(f.readline().split(sCom, 2)[1]) #2nd line: initial X position [m]; it will not actually be used
	xStep = float(f.readline().split(sCom, 2)[1]) #3rd line: step vs X [m]
	xNp = int(f.readline().split(sCom, 2)[1]) #4th line: number of points vs X
	yStart = float(f.readline().split(sCom, 2)[1]) #5th line: initial Y position [m]; it will not actually be used
	yStep = float(f.readline().split(sCom, 2)[1]) #6th line: step vs Y [m]
	yNp = int(f.readline().split(sCom, 2)[1]) #7th line: number of points vs Y
	zStart = float(f.readline().split(sCom, 2)[1]) #8th line: initial Z position [m]; it will not actually be used
	zStep = float(f.readline().split(sCom, 2)[1]) #9th line: step vs Z [m]
	zNp = int(f.readline().split(sCom, 2)[1]) #10th line: number of points vs Z
	totNp = xNp*yNp*zNp
	locArBx = array('d', [0]*totNp)
	locArBy = array('d', [0]*totNp)
	locArBz = array('d', [0]*totNp)
	for i in range(totNp):
		curLineParts = f.readline().split('\t')
		locArBx[i] = float(curLineParts[0])
		locArBy[i] = float(curLineParts[1])
		locArBz[i] = float(curLineParts[2])
	f.close()
	xRange = xStep
	if xNp > 1: xRange = (xNp - 1)*xStep
	yRange = yStep
	if yNp > 1: yRange = (yNp - 1)*yStep
	zRange = zStep
	if zNp > 1: zRange = (zNp - 1)*zStep

	print('')
	print('locArBx', np.sum(locArBx))
	print('locArBy', np.sum(locArBy))
	print('locArBz', np.sum(locArBz))
	print('xnp', xNp) # int
	print(('ynp', yNp)) # int
	print(('znp', zNp)) # int
	print(('xrange', xRange)) # float
	print(('yrange', yRange)) # float
	print(('zrange', zRange)) # float
	print('......')

	return SRWLMagFld3D(locArBx, locArBy, locArBz, xNp, yNp, zNp, xRange, yRange, zRange, 1)

# borondics ********************** Auxiliary function to generate tabulated 3D Magnetic Field data
# can be used for looping wavelengths, optimization, etc.
def myBend(H, Xc, s, n, xstart, xend, npoints):
	'''
	H - magnetic field [T]
	Xc - center of straight section [mm]
	s - length of bend
	n - transition steepness
	xstart, xend, npoints 
	'''
	global xStart,xStep,xNp,yStart,yStep,yNp,zStart,zStep,zNp
	dummy = 0.0
	xStart = dummy/1000 # convert to [m]
	xStep = dummy
	xNp = 1
	yStart = dummy/1000 # convert to [m]
	yStep = dummy 
	yNp = 1
	zStart = xstart/1000 # convert to [m]
	zStep = (xend - xstart) / (npoints - 1)
	zNp = npoints
	
	x = np.linspace(xstart, xend, npoints)
	Bx = np.zeros(npoints)
	By = H * (((((x+Xc)**n)/(s**n + (x+Xc)**n))-1) + ((((x-Xc)**n)/(s**n + (x-Xc)**n))-1))
	Bz = np.zeros(npoints)

	return SRWLMagFld3D(array('d', Bx), array('d', By), array('d', Bz), 1, 1, npoints, 0.0, 0.0, 16.0, 1)
	

startTime = time.time()

# ======== Folders and File names ========================================================================================= #
#**********************Input Parameters:
strExDataFolderName = 'Input_field' #example data sub-folder name
arFldInFileNames = ['B1-B1_mag_fld_original.dat'] #3D Magnetic Field data file names
strTrajOutFileName = 'ex01_res_traj_off.dat' #file name for output trajectory data
prefix = 'B1-B1'

# ========================================================================================================================= #
# TRAJECTORY CALCULATION
# ========================================================================================================================= #
# ********************** Defines particle object
part = SRWLParticle()

# Initial Transverse Coordinates (initial Longitudinal Coordinate will be defined later on) [m]
part.x = -119.72390066034988e-3 #-119.34276573016813e-3			  
part.y = 0.0
# Initial Transverse Velocities
part.xp = 48.089011080417869e-3					  
part.yp = 0
part.gamma = 3/0.51099890221e-03	#Relative Energy
part.relE0 = 1					  #Electron Rest Mass
part.nq = -1						#Electron Charge

# Number of Points for Trajectory calculation
npTraj = 16001 

# Magnetic Field Interpolation Method, to be entered into 3D field structures below 
#	(to be used e.g. for trajectory calculation):
# 1- bi-linear (3D), 2- bi-quadratic (3D), 3- bi-cubic (3D), 4- 1D cubic spline (longitudinal) + 2D bi-cubic
fieldInterpMeth = 4 #2 

#General Precision parameters for Trajectory calculation:
#[0]: integration method No:
	#1- fourth-order Runge-Kutta (precision is driven by number of points)
	#2- fifth-order Runge-Kutta
#[1],[2],[3],[4],[5]: absolute precision values for X[m],X'[rad],Y[m],Y'[rad],Z[m] (yet to be tested!!) - 
# 	to be taken into account only for R-K fifth order or higher
#[6]: tolerance (default = 1) for R-K fifth order or higher
#[7]: max. number of auto-steps for R-K fifth order or higher (default = 5000)
arPrecPar = [1] 

#********************** Defining Magnetic Field:
xcID = 0 #Transverse Coordinates of ID Center [m]
ycID = 0
zcID = 0 #Longitudinal Coordinate of ID Center [m]

magFldCnt = SRWLMagFldC() #Container
magFldCnt.allocate(1) #Magnetic Field consists of 1 part

#MagFld =SRWLMagFldM(0.6, 1, 'n', 2.0)

print('   Reading magnetic field data from files ... ', end='')
for i in range(1): # handle this nicer for the myBend setup. We need a dictionary in the beginning
	filePath = os.path.join(os.getcwd(), strExDataFolderName, arFldInFileNames[i])
# 	# original magFldCnt.arMagFld[i] = AuxReadInMagFld3D(filePath, '#')
	magFldCnt.arMagFld[i] = myBend(1.72, 12405/2, 1050/2, 60, -8000, 8000, npTraj)
	magFldCnt.arMagFld[i].interp = fieldInterpMeth
	magFldCnt.arXc[i] = xcID
	magFldCnt.arYc[i] = ycID
# 	plt.plot(magFldCnt.arMagFld[i].arBy)
# 	plt.show()

print('done')

magFldCnt.arMagFld[0].nRep = 1 #Central part of ID

magFldCnt.arZc[0] = zcID
part.z = zcID -0.5*magFldCnt.arMagFld[0].rz

#********************** Trajectory structure, where the results will be stored
partTraj = SRWLPrtTrj()
partTraj.partInitCond = part
#partTraj.allocate(npTraj)
partTraj.allocate(npTraj, True)
partTraj.ctStart = 0 #Start Time for the calculation
#partTraj.ctEnd = (numPer + 2)*per + magFldCnt.arMagFld[0].rz + magFldCnt.arMagFld[2].rz #End Time
partTraj.ctEnd = magFldCnt.arMagFld[0].rz


if(1): # (1) calculates and plots trajectory, (0) do not plot trajectory
	print('   Performing calculation ... ', end='') 
	partTraj = srwl.CalcPartTraj(partTraj, magFldCnt, arPrecPar) # ACTUAL FUNCTION THAT PERFORMS TRAJECTORY CALCULATION
	print('done')
	#**********************Saving results to a file
	print('   Saving trajectory data to a file ... ', end='')
	#partTraj.save_ascii(os.path.join(os.getcwd(), 'SRW_traj_from_FAC_field_half.txt'))
	print('done')
	#**********************Plotting results
	print('   Plotting the results (blocks script execution; close any graph windows to proceed) ... ', end='')
	ctMesh = [partTraj.ctStart, partTraj.ctEnd, partTraj.np]
	for i in range(partTraj.np):
		partTraj.arX[i] *= 1000
		partTraj.arXp[i] *= 1000
	   
# 	uti_plot1d(partTraj.arX, ctMesh, ['ct [m]', 'Horizontal Position [mm]'])	
# 	uti_plot1d(partTraj.arXp, ctMesh, ['ct [m]', 'Horizontal Angles [mrad]'])
# 	uti_plot_show() #show all graphs (and block execution)
	print('done')

# 	# borondics - plotting directly (could be used for savefig or save data and then plot)
# 	f = plt.figure(figsize=(12,7))
# 	plt.plot(np.linspace(ctMesh[0], ctMesh[1], ctMesh[2]), partTraj.arX, 'r')
# 	plt.tight_layout()
# 	plt.show()
	
# ========================================================================================================================= #
# WAVEFRONT SAMPLING
# ========================================================================================================================= #

if(1): # (1) Calculate Wavefront, (0) do not run wavefront 
	
	#********************* Radiation Sampling
	wl = 10						 # Wavelength [µm]
	center_energy = wl2eV(wl)  # Energy [eV]
	print('   Wavelength = {0:.2f} um'.format(wl))
	OBS_h = 50/1000				   # Observation window width in (m)
	OBS_v = 50/1000				   # Observation window heigth in (m)
	OBS_z = 6.5					# Longitudinal observation position (m)
	OBS_hCEN = 0.000					# Observation window center position (m)

	
	#********************* Electron Beam
	beam = SRWLPartBeam()
	beam.Iavg = 500/1000		#Average current [A]
	print('   Current = {0:.2f} mA'.format(beam.Iavg*1000))
	beam.partStatMom1.x = 0#-119.72390066034988e-3				  #initial transverse positions [m]
	beam.partStatMom1.y = 0.0
	beam.partStatMom1.z = 0#-8.0				   #initial longitudinal positions (set in the middle of undulator)
	beam.partStatMom1.xp = 0#48.089011080417869e-3				  #initial relative transverse velocities
	beam.partStatMom1.yp = 0.0
	beam.partStatMom1.gamma = 3.0/0.5109989e-03 #3./0.51099890221e-03 #relative energy
	sigEperE = 0.00085 #relative RMS energy spread
	sigX = 64.90e-06 #horizontal RMS size of e-beam [m]
	sigXp = 3.80e-06 #horizontal RMS angular divergence [rad]
	sigY = 3.40e-06 #vertical RMS size of e-beam [m]
	sigYp = 0.70e-06 #vertical RMS angular divergence [rad]
	#2nd order stat. moments:
	beam.arStatMom2[0] = sigX*sigX #<(x-<x>)^2> 
	beam.arStatMom2[1] = 0 #<(x-<x>)(x'-<x'>)>
	beam.arStatMom2[2] = sigXp*sigXp #<(x'-<x'>)^2> 
	beam.arStatMom2[3] = sigY*sigY #<(y-<y>)^2>
	beam.arStatMom2[4] = 0 #<(y-<y>)(y'-<y'>)>
	beam.arStatMom2[5] = sigYp*sigYp #<(y'-<y'>)^2>
	beam.arStatMom2[10] = sigEperE*sigEperE #<(E-<E>)^2>/<E
	
	
	#********************* Wavefront initial structure
	wfr = SRWLWfr()						 #For intensity distribution at fixed photon energy
	wfr.allocate(1, 400, 400)			   #Numbers of points vs Photon Energy, Horizontal and Vertical number of pixels
	wfr.mesh.zStart = OBS_z				 #Longitudinal Position [m] at which SR has to be calculated
	wfr.mesh.eStart = center_energy		 #Initial Photon Energy [eV]
	wfr.mesh.eFin = center_energy		   #Final Photon Energy [eV]
	wfr.mesh.xStart = OBS_hCEN-OBS_h/2	  #Initial Horizontal Position [m]
	wfr.mesh.xFin = OBS_hCEN+OBS_h/2		#Final Horizontal Position [m]
	wfr.mesh.yStart = -OBS_v/2			  #Initial Vertical Position [m]
	wfr.mesh.yFin = OBS_v/2				 #Final Vertical Position [m]
	wfr.partBeam = beam
	#meshInitPartCoh = deepcopy(wfr.mesh)
	
	
	#***********Precision for Trajectory
	meth = 2					#SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
	relPrec = 0.1			  #relative precision
	zStartInteg = zStart		#-9.0 #longitudinal position to start integration (effective if < zEndInteg)
	zEndInteg = OBS_z-0.1	 #9.0 #longitudinal position to finish integration (effective if > zStartInteg)
	npTraj = zNp			  #Number of points for trajectory calculation 
	useTermin = 1			   #Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)
	sampFactNxNyForProp = 0	 #sampling factor for adjusting nx, ny (effective if > 0)
	arPrecParR = [meth, relPrec, zStartInteg, zEndInteg, npTraj, useTermin, sampFactNxNyForProp]
	
	print('   Performing Electric Field (wavefront at fixed photon energy) calculation ... ', end='')
	srwl.CalcElecFieldSR(wfr, 0, magFldCnt, arPrecParR) # Calculate by magnetic field 
#	srwl.CalcElecFieldSR(wfr, partTraj, 0, arPrecParR) # Calculate by electron trajectory 
	print('done')



	# ===================================================================================================================== #
	# WAVEFRONT PROPAGATION
	# ===================================================================================================================== #
	
	if(0): # (1) Propagate through thin lenses, (2) do not propagate wavefront
		
#		Drift_10m_Lens = SRWLOptD(5.0)
		
#		Ap = SRWLOptA('r', 'a', 40.0e-3, 40.0e-3)		
		
		d = 12
		Lens = SRWLOptL(6, 6)	
		Drift_Lens_obs = SRWLOptD(d)		
#		ppDrift_10m_Lens =	[ 0,  0, 1.5,  1,  0, 2.0, 2.0, 2.0, 2.0,  0,  0,   0]
#		ppAp =				[ 0,  0, 1.5,  0,  0, 1.0, 4.0, 1.0, 4.0,  0,  0,   0]
#		OptElement=		  [[1],[2],[3], [4],[5],[6], [7], [8], [9], [10],[11],[12]]   Example	 
		ppLens =			  [ 0,  0, 1.5,  0,  0, 1.0, 8.0, 1.0, 8.0,  0,  0,   0]
		ppDrift_Lens_obs =	[ 0,  0, 1.5,  0,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]
		ppFin  =			  [ 0,  0, 1.5,  0,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]
		
		#Wavefront Propagation Parameters:
		#[1]: Auto-Resize (1) or not (0) Before propagation
		#[2]: Auto-Resize (1) or not (0) After propagation
		#[3]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
		#[4]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
		#[5]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
		#[6]: Horizontal Range modification factor at Resizing (1. means no modification)
		#[7]: Horizontal Resolution modification factor at Resizing
		#[8]: Vertical Range modification factor at Resizing
		#[9]: Vertical Resolution modification factor at Resizing
		#[10]: Type of wavefront Shift before Resizing (not yet implemented)
		#[11]: New Horizontal wavefront Center position after Shift (not yet implemented)
		#[12]: New Vertical wavefront Center position after Shift (not yet implemented)



		#******************* Observation plane		
		
		# At Lens
#		optBL = SRWLOptC([Drift_10m_Lens],[ppDrift_10m_Lens])		
		# After Lens
#		optBL = SRWLOptC([Lens],[ppLens]) 
		# At OBS
#		optBL = SRWLOptC([Lens, Drift_Lens_obs],[ppLens, ppDrift_Lens_obs]) 
		
		# After Aperture		
#		optBL = SRWLOptC([Ap],[ppAp])		
		# After Aperture and Lens		
#		optBL = SRWLOptC([Ap, Lens, Drift_Lens_obs],[ppAp, ppLens, ppDrift_Lens_obs]) # Container of optical elements and propagation parameters
		optBL = SRWLOptC([Lens, Drift_Lens_obs],[ppLens, ppDrift_Lens_obs]) # Container of optical elements and propagation parameters
#		optBL = SRWLOptC([Drift_Lens_obs],[ppDrift_Lens_obs]) # Container of optical elements and propagation parameters
#		optBL = SRWLOptC([Lens, Drift_Lens_obs],[ppLens, ppDrift_Lens_obs, ppFin]) # Container of optical elements and propagation parameters
		
		
		print('   Simulating Electric Field Wavefront Propagation ... ', end='')
		srwl.PropagElecField(wfr, optBL) # ACTUAL FUNCTION THAT PROPAGATES WAVEFRONT
		print('done')
		# ===================================================================================================================== #


		 
	print('   Extracting Intensity from calculated Electric Field ... ', end='')
	mesh1 = deepcopy(wfr.mesh)
	arI1 = array('f', [0]*mesh1.nx*mesh1.ny) #"flat" 2D array to take intensity data
	srwl.CalcIntFromElecField(arI1, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0) # Intensity from electric field (Wfr_name, Disp_Name, Polariz., FuncOf, Extr, ObsEi , ConstHorPos, ConstVerPos, NewDisp)
#EXP int CALL srwlCalcIntFromElecField(char* pInt, SRWLWfr* pWfr, char pol, char intType, char depType, double e, double x, double y);
# * @param [out] pInt pointer to resulting Intensity (array)
# * @param [in] pWfr pointer to pre-calculated Wavefront structure
# * @param [in] pol polarization component to extract: 
# *			 0- Linear Horizontal; 
# *			 1- Linear Vertical; 
# *			 2- Linear 45 degrees; 
# *			 3- Linear 135 degrees;
# *			 4- Circular Right; 
# *			 5- Circular Left; 
# *			 6- Total
# * @param [in] intType "type" of a characteristic to be extracted: 
# *			 0- "Single-Electron" Intensity; 
# *			 1- "Multi-Electron" Intensity; 
# *			 2- "Single-Electron" Flux; 
# *			 3- "Multi-Electron" Flux; 
# *			 4- "Single-Electron" Radiation Phase; 
# *			 5- Re(E): Real part of Single-Electron Electric Field;
# *			 6- Im(E): Imaginary part of Single-Electron Electric Field;
# *			 7- "Single-Electron" Intensity, integrated over Time or Photon Energy (i.e. Fluence);
# * @param [in] depType type of dependence to extract: 
# *			 0- vs e (photon energy or time);
# *			 1- vs x (horizontal position or angle);
# *			 2- vs y (vertical position or angle);
# *			 3- vs x&y (horizontal and vertical positions or angles);
# *			 4- vs e&x (photon energy or time and horizontal position or angle);
# *			 5- vs e&y (photon energy or time and vertical position or angle);
# *			 6- vs e&x&y (photon energy or time, horizontal and vertical positions or angles);
# * @param [in] e photon energy (to keep fixed)
# * @param [in] x horizontal position (to keep fixed)
# * @param [in] y vertical position (to keep fixed)
	print('done')

	
	print("nx vs ny = {0} x {1}".format(wfr.mesh.nx, wfr.mesh.ny)) 
	print('Total time = {0:.1f} minutes'.format((time.time() - startTime)/60.0))

	# ===================================================================================================================== #
	# PLOTS
	# ===================================================================================================================== #

	if(1):
		uti_plot2d(arI1, [1000*wfr.mesh.xStart, 1000*wfr.mesh.xFin, wfr.mesh.nx], [1000*wfr.mesh.yStart, 1000*wfr.mesh.yFin, wfr.mesh.ny], ['Horizontal Position [mm]', 'Vertical Position [mm]', 'wfr, arI1 - Intensity at ' + str(wfr.mesh.eStart) + ' eV'])

	#	 print(len(arI1))
	#  # borondics - plotting directly
	# 	import colorcet as cc
	# 	a = np.array(arI1).reshape([400, 400])
	# 	f = plt.figure(figsize=(10,10))
	# 	plt.imshow(a, cmap=cc.cm['linear_kry_5_95_c72'])
	# 	plt.tight_layout()
	# 	plt.savefig('a.png', dpi=300)
	#	 plt.show()
	
		arI2x = array('f', [0]*wfr.mesh.nx) #array to take 1D intensity data (vs X)
		srwl.CalcIntFromElecField(arI2x, wfr, 6, 0, 1, wfr.mesh.eStart, 0, 0)
	#	I2x_fwhm = 1000*mypkg.get_fwhm(np.linspace(1000*wfr.mesh.xStart, 1000*wfr.mesh.xFin, wfr.mesh.nx), arI2x)
		uti_plot1d(arI2x, [1000*wfr.mesh.xStart, 1000*wfr.mesh.xFin, wfr.mesh.nx], ['Horizontal Position [mm]', 'Intensity [ph/s/.1%bw/mm^2]', 'Intensity at ' + str(wfr.mesh.eStart) + ' eV'])
	
		arI2y = array('f', [0]*wfr.mesh.ny) #array to take 1D intensity data (vs Y)
		srwl.CalcIntFromElecField(arI2y, wfr, 6, 0, 2, wfr.mesh.eStart, 0, 0)
	#	I2y_fwhm = 1000*mypkg.get_fwhm(np.linspace(1000*wfr.mesh.yStart, 1000*wfr.mesh.yFin, wfr.mesh.ny), arI2y)
		uti_plot1d(arI2y, [1000*wfr.mesh.yStart, 1000*wfr.mesh.yFin, wfr.mesh.ny], ['Vertical Position [mm]', 'Intensity [ph/s/.1%bw/mm^2]', 'Intensity at ' + str(wfr.mesh.eStart) + ' eV'])

		if(0):
			int_x = np.array(arI2x)
			int_y = np.array(arI2y)
			x_coord = np.linspace(wfr.mesh.xStart, wfr.mesh.xFin, wfr.mesh.nx)
			y_coord = np.linspace(wfr.mesh.yStart, wfr.mesh.yFin, wfr.mesh.ny)
			np.savetxt(prefix+'_IX.txt', np.array([x_coord, int_x]).transpose(), fmt='%.6e', delimiter='\t')
			np.savetxt(prefix+'_IY.txt', np.array([y_coord, int_y]).transpose(), fmt='%.6e', delimiter='\t')

		if(0):
			int_2D = np.zeros((len(y_coord)+1, len(x_coord)+1))
			int_2D[0, 1:] = x_coord
			int_2D[1:, 0] = y_coord
			int_2D[1:,1:] = np.array(arI1).reshape((len(y_coord), len(x_coord)))
			np.savetxt(prefix+'_2D.txt', int_2D, fmt='%.6e', delimiter='\t')
		
		if(0):
			ElecFldX = np.array(wfr.arEx)
			ElecFldY = np.array(wfr.arEy)
			np.savetxt(prefix+'_EX.txt', ElecFldX, fmt='%.6e', delimiter='\t')
			np.savetxt(prefix+'_EY.txt', ElecFldY, fmt='%.6e', delimiter='\t')
		
		uti_plot_show() #show all graphs (and block execution) 
	
		#**********************Plotting trajectory








