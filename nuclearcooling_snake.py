import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Differences from coaxial steam generator:
# 	The geometry of this reactor are tubes going from z=0 to z=h with
# a slope of dxdz at all times, where x is any direction perpendicular
# to z (the direction of x can change). This models, for example, tubes
# spiraling upwards with constant slope or approximates tubes zig-zagging
# upwards.
# 	Total length of one tube is no longer equal to h, variable l denotes
# the total length of one tube. l is calculated from h and dxdz, where
# dxdz is a parameter that needs to be provided. It is assumed that the
# slope of the tubes is approximately constant and that dxdz >> 1.
# 	Axh has been changed, its now calculated from l and not h. The 
# function getNusseltPb has been changed to the formula for cross flow.
# Lead flow cross sectional area corrected, sloped tubes take up more 
# of the cross sectional area.

## PHYSICAL PARAMETERS ##

# Thermal conductivity of tube material [W/(m*K)]
kss = 20
# Desired thermal power [W]
Q = 14e6
# Height of reactor [m]
h = 0.64
# Diameter of reactor [m]
D = 0.4
# Number of tubes
n = 18
# Inner diameter of tubes [m]
di = 11*0.001
# Thickness of tubes [m]
dd = 1*0.001
# "Slope" of the tubes, the x direction is any direction perpendicular to z.
# For example, dxdz=0 gives coaxial flow (will be inaccurate because the nusselt
# number assumes perpendicular flow), dxdz=1 gives 45 degree tubes, and so on.
# The value of dxdz should be >> 1.
dxdz = 22
# Lead inflow temperature [C]
T0Pb = 550
# Water inflow temperature (steam if above 342.11C, otherwise liquid) [C]
T0W = 100
# Desired temperature of water at outflow [C]
T1W = 500
# Desired temperature of lead at outflow [C]
T1Pb = 500
# Lowest alowed lead temperature
TminPb = 350


## SIMULATION PARAMETERS ##

# Number of length elements
N = 200
# Amount times to compute temperatures, for each iteration the calculation becomes more accurate
cycles = 20
# If True, print useful data about the solution after the simulation has completed
printData = True
# If True, will plot every intermediate state, otherwise it will just plot the final, most accurate, state
plotIntermediates = False
# If True will print a progress bar as it's simulating
progressBar = True
# If True will plot the "convegence error", a low error means the solution has converged
# (not necessarily to the right solution). A systematic error means cycles should be increased,
# a noisy error means that the amount of cycles is fine.
plotError = True


## PHYSICAL CONSTANTS ##

# pi
pi = np.pi
# Density of lead at inflow [kg/m^3]
rhoPb = 10.678*1000
# Density of water at inflow [kg/m^3]
rhoW = 965.19736410518
# Enthalpy of water at boiling point [J/kg]
HWb0 = 1609.829098436e3
# Enthalpy of steam at boiling point [J/kg]
HWb1 = 2611.4108223058e3


## CALCULATED CONSTANTS ##

# Step size [m]
dz = h/(N - 1)
# Radius of reactor [m]
R = D/2
# Inner radius of tubes [m]
ri = di/2
# Outer diameter of tubes [m]
do = di+2*dd
# Outer radius of tubes [m]
ro = do/2
# Change of l with respect to z
dldz = (1 + dxdz**2)**0.5
# Total length of one tube
l = dldz*h
# Total heat exchange area [m^2]
Ahx = do*pi*l*n
# Water flow total cross sectional area [m^2]
Aw = ri*ri*pi*n
# Lead flow cross sectional area [m^2]
Apb = R*R*pi - ro*ro*pi*n*dldz
# Difference in water enthalpy between outflow and inflow
deltaHW = 0
# Difference in lead enthalpy between outflow and inflow
deltaHPb = 0
# Water mass flow rate [kg/s]
#mDotW = uW*Aw*rhoW
mDotW = 0
# Lead mass flow rate [kg/s]
#mDotPb = uPb*Apb*rhoPb
mDotPb = 0
# Velocity of water at inflow [m/s]
uW = 0
# Velocity of lead at inflow [m/s]
uPb = 0

def updateConstants():
	global dz, R, ri, do, ro, dldz, l, Ahx, Aw, Apb, deltaHPb, deltaHW, mDotPb, mDotW, uPb, uW
	# Step size [m]
	dz = h/(N - 1)
	# Radius of reactor [m]
	R = D/2
	# Inner radius of tubes [m]
	ri = di/2
	# Outer diameter of tubes [m]
	do = di+2*dd
	# Outer radius of tubes [m]
	ro = do/2
	# Change of l with respect to z
	dldz = (1 + dxdz**2)**0.5
	# Total length of one tube
	l = dldz*h
	# Total heat exchange area [m^2]
	Ahx = do*pi*l*n
	# Water flow total cross sectional area [m^2]
	Aw = ri*ri*pi*n
	# Lead flow cross sectional area [m^2]
	Apb = R*R*pi - ro*ro*pi*n*dldz
	# Difference in water enthalpy between outflow and inflow
	deltaHW = getHW(T1W)-getHW(T0W)
	# Difference in lead enthalpy between outflow and inflow
	deltaHPb = getHPb(T1Pb)-getHPb(T0Pb)
	# Water mass flow rate [kg/s]
	#mDotW = uW*Aw*rhoW
	mDotW = Q/deltaHW
	# Lead mass flow rate [kg/s]
	#mDotPb = uPb*Apb*rhoPb
	mDotPb = Q/deltaHPb
	# Velocity of water at inflow [m/s]
	uW = mDotW/Aw/rhoW
	# Velocity of lead at inflow [m/s]
	uPb = mDotPb/Apb/rhoPb

def lerp(x, xList, yList):
	"""Returns an approximation of y(x), xList and yList should contain sampled values from y(x).
	xList has to contain x-values in increasing order."""
	if len(xList) != len(yList) or len(xList) == 0:
		print("Error: dimensions of lists dont match or equals zero.")
		return None
	#Binary search, find the neighbouring x-values
	i0 = 0
	i1 = len(xList) - 1
	x0 = xList[i0]
	x1 = xList[i1]
	if x > x1 or x < x0:
		print("Error: given value is out of bounds: " + str(x))
		return None
	while(i1 - i0 != 1):
		iMid = int((i0+i1)/2)
		xMid = xList[iMid]
		if xMid > x:
			i1 = iMid
		else:
			i0 = iMid
	xminus = xList[i0]
	xplus = xList[i1]
	p = (x-xminus)/(xplus-xminus)
	return yList[i0]*(1-p) + yList[i1]*p
	
#Sampled water temperatures [C]
TWsamp = [
	0, 					50, 				100, 				150, 
	200, 				250, 				280, 				300,
	320, 				330, 				337, 				342.11,
	342.12, 			342.2, 				345,				350, 
	360,				370, 				380, 				400, 
	440, 				480, 				520, 				580]
#Sampled specific enthalpy values [J/kg]
HWsamp = [
	15.069416858762e3, 	222.22555353748e3, 	430.32079609592e3, 	641.34032434486e3, 
	858.1170925424e3, 	1086.0355962798e3, 	1232.7888617023e3, 	1338.0632609535e3,
	1453.846049747e3, 	1518.639370056e3, 	1568.8433826561e3, 	1609.7439498365e3,
	1609.829098436e3, 	2611.4108223058e3, 	2644.4722229326e3, 	2692.9998154054e3, 
	2769.5620385581e3,	2831.4049104836e3, 	2884.6070536591e3, 	2975.5476758035e3,
	3124.5831873807e3, 	3251.7604839088e3, 	3367.7855983754e3, 	3530.7530092137e3]
#Sampled isobar specific heat capacity [kJ/kgK]
CWsamp = [
	4.1501246032742, 	4.1465669400481, 	4.1837678599079, 	4.265896612902, 
	4.4219269256657, 	4.7324546776592, 	5.0833000211663, 	5.4760165314084,
	6.183875433228, 	6.8289169597204, 	7.5793510257999, 	8.5136552680353,
	8.5160651264762, 	12.941502699273,	10.850428622717, 	8.7885126098705, 
	6.7740330631997,	5.6867522089636, 	5.0002808341009, 	4.1777776193013,
	3.386350240844, 	3.0126365014173, 	2.807448480426, 	2.644747678729]
#Sampled density at 15 MPA [kg/m^3]
rhoWsamp = [
	1007.2953533373, 	994.42906842479, 	965.19736410518, 	925.02855025451, 
	874.51144816007, 	811.02408541955, 	763.5659317279, 	725.55327522717, 
	678.75763644934	, 	649.61853654555, 	625.1842416931, 	603.73619073597, 
	603.68980109742, 	96.641787510605, 	92.58910388291, 	87.102683889095, 
	79.476556981261, 	74.112339548535, 	69.982275774721, 	63.811881187966, 
	55.664233554794, 	50.190066624199, 	46.090496592759, 	41.418643769644]
#Sampled (dynamic) viscosity at 15 MPA [Pa s]
muWsamp = [
	0.0017569222255541,	0.00054957548270974,0.00028572774719823, 0.00018609826343681, 
	0.00013761834458872,0.00010910372865988,9.6236453459579E-5, 8.8333362476973E-5,
	8.0272953163543E-5, 7.5857009417524E-5, 7.2406568846584E-5, 6.9527604710848E-5, 
	6.9521509303783E-5, 2.2793338904066E-5, 2.2822905858684E-5, 2.2935282504602E-5, 
	2.3263415562689E-5, 2.3651679469595E-5, 2.4066651122984E-5, 2.4930042731232E-5, 
	2.6687629398083E-5, 2.8429697881328E-5, 3.0139941628777E-5, 3.2640892162702E-5]
#Sampled thermal conductivity at 15 MPA [W/mK]
kWsamp = [
	0.56930910448506, 	0.65051792015857, 	0.68721569070846, 	0.69177924539095, 
	0.67486517874436, 	0.63476052513029, 	0.59568357949575, 	0.56143592682432, 
	0.51992444478752, 	0.49636008973223, 	0.47845802606602, 	0.46406300236298, 
	0.46403325075079, 	0.11521314787269, 	0.10844451493895, 	0.1007924810207, 
	0.092185222216949, 	0.0872445138624, 	0.084117453069393, 	0.080681703274648, 
	0.079006161973628, 	0.080461037513768, 	0.083484892747535, 	0.089621105028219]

## WATER PROPERTIES ##

def getHW(T):
	"""Specific enthalpy of water [J/kg] at given temperature [C]"""
	return lerp(T, TWsamp, HWsamp)

def getTW(H):
	"""Temperature of water [C] at given specific enthalpy [J/kg]"""
	return lerp(H, HWsamp, TWsamp)

def getCW(H):
	"""Isobar specific heat capacity [J/kgK] at given specific enthalpy [J/kg]"""
	return lerp(H, HWsamp, CWsamp)*1000

def getRhoW(H):
	"""Density of water [kg/m^3] at given specific enthalpy [J/kg]"""
	return lerp(H, HWsamp, rhoWsamp)

def getkW(H):
	"""Thermal conductivity of water [W/mK] at given specific enthalpy [J/kg]"""
	return lerp(H, HWsamp, kWsamp)

def getMuW(H):
	"""Dynamic viscosity of water [Pa s] at given specific enthalpy [J/kg]"""
	return lerp(H, HWsamp, muWsamp)

def getWaterVel(H):
	"""Velocity of water inside tubes [m/s] at given specific enthalpy [J/kg]"""
	return mDotW/(getRhoW(H)*Aw)

def getPrandtlW(H):
	"""Returns the prandtl number for water at given specific enthalpy [J/kg]"""
	return getMuW(H)*getCW(H)/getkW(H)

def getReynoldW(H):
	"""Returns Reynolds number for water at given specific enthalpy [J/kg]"""
	return di*abs(getWaterVel(H))*getRhoW(H)/getMuW(H)

def getNusseltW(H):
	"""Returns the nusselt number for water at given specific enthalpy [J/kg]"""
	if H <= HWb0:
		return 0.023*getReynoldW(H)**0.8*getPrandtlW(H)**0.4
	elif H < HWb1:
		return getNusseltW(HWb0)*5
	else:
		return 0.0157*getReynoldW(H)**0.84*getPrandtlW(H)**0.33

def getConvHtW(H):
	"""Returns the convective heat transfer coefficient, h, [W/m^2K] for water at given specific enthalpy [J/kg]"""
	return getNusseltW(H)*getkW(H)/di

## LEAD PROPERTIES ##

def getHPb(T):
	"""Specific enthalpy of lead [J/kg] at given temperature [C]"""
	return 145*T #from dmitry

def getTPb(H):
	"""Temperature of lead [C] at given specific enthalpy [J/kg]"""
	return H/145 #from dmitry

def getLambdaPb(T):
	"""Thermal conductivity of molten lead [W/mK] at given temperature [C]"""
	return 9.2 + 0.011*(T + 273.15) #from dmitry

def getAlphaPb(T):
	"""Thermal diffusivity of molten lead [m/s^2] at given temperature [C]"""
	return (3.408 + 0.0112*(T + 273.15))*1e-6 #from dmitry

def getPecletPb(T):
	"""Returns the peclet number for lead at given temperature [C]"""
	return do*abs(uPb)/getAlphaPb(T)

def getNusseltPb(T):
	"""Returns the nusselt number for lead at given temperature [C]"""
	#return 6 + 0.006*getPecletPb(T) # Lead flows parallel to tubes
	return 4.003 + 0.228*getPecletPb(T)**0.67 # Cross flow, tubes perpendicular to lead flow

def getConvHtPb(T):
	"""Returns the convective heat transfer coefficient, h, [W/m^2K] for lead at given temperature [C]"""
	return getNusseltPb(T)*getLambdaPb(T)/do

# Constant factors used when calculating Reff
cf1 = 1/(2*pi*ri)
cf2 = np.log(ro/ri)/(2*pi*kss)
cf3 = 1/(2*pi*ro)
def getReff(HW, TPb):
	"""Returns the overall thermal resistance coefficient, Reff, [mK/W] for given water and lead temp [C]"""
	return cf1/getConvHtW(HW) + cf2 + cf3/getConvHtPb(TPb)

def getQW(Hw):
	"""Returns the energy flow of water for given enthalpy [MW]"""
	return mDotW*Hw*1e-6
	
def getQPb(Hpb):
	"""Returns the energy flow of lead for given enthalpy [MW]"""
	return mDotPb*Hpb*1e-6

def plotEnthalpyData():
	_, (ax1, ax2) = plt.subplots(2, 1, True)
	ax1.plot(TWsamp, HWsamp)
	ax1.set_title("Water")
	plt.xlabel("Temperature [Celcius]")
	ax1.set_ylabel("Specific enthalpy [J/kg]")
	ax2.plot(TWsamp, CWsamp)
	ax2.set_ylabel("Specific heat capacity [kJ/kgK]")
	plt.show()

def printSolutionData(Hw, Hpb):
	print("Data for simulation with N = {0}, and {1} cycles".format(N, cycles))
	print("Reactor dimensions: height = {0} m, diameter = {1} m".format(h,D))
	print("Fraction of volume consisting of water and tube material: {0}".format(ro*ro/(R*R)*dldz*n))
	print("TUBE DATA:")
	print("Dimensions (1 tube):  length: {0:.1f}  m,  outer diameter:    {1:.1f} mm,      thickness: {2:.1f} mm".format(l, do*1000, dd*1000))
	print("Geometry:             Number of tubes:     {0}, slope of tubes dx/dz: {1:.1f}".format(n,dxdz))
	print("Material properties:  Conductivity: {0} W/mK".format(kss))
	print("")
	Qw = getQW(Hw[N-1]) - getQW(Hw[0])
	Qpb = getQPb(Hpb[N-1]) - getQPb(Hpb[0])
	print("")
	print("SOLUTION DETAILS                WATER        LEAD         PHYSICAL DIMENSION")
	print("Flow areas                    {0:.3f}       {1:.3f}               m^2       ".format(Aw, Apb))
	print("Inflow Temp                   {0:.3f}       {1:.3f}              deg. C     ".format(getTW(Hw[0]), getTPb(Hpb[-1])))
	print("Outflow Temp                  {0:.3f}       {1:.3f}              deg. C     ".format(getTW(Hw[-1]), getTPb(Hpb[0])))
	print("Mass flow (pos z)             {0:.3f}       {1:.3f}               kg/s      ".format(mDotW, mDotPb))
	print("Inflow vel. (pos z)           {0:.3f}       {1:.3f}               m/s       ".format(uW, uPb))
	print("Thermal pwr gain              {0:.3f}       {1:.3f}               MW        ".format(Qw, Qpb))
	print("Solution discrepancy (lower is better): 10^{0:.0f} W".format(np.log10(checkSolution(Hw, Hpb))))
##### SIMULATION #####

# Derivative of water enthalpy with respect to z
def dHWdz(TPb, TW, HW):
	return (TPb-TW)*dldz/(mDotW*getReff(HW, TPb))*n
# Derivative of lead enthalpy with respect to z
def dHPbdz(TPb, TW, HW):
	return (TW-TPb)*dldz/(mDotPb*getReff(HW, TPb))*n

# Main function
def simulate(printSol = printData, prgBar = progressBar):
	if prgBar:
		print('[', end='')
	# Water specific enthalpy at inflow (z=0)
	H0W = getHW(T0W)
	# Lead specific enthalpy at inflow (z=h)
	H0Pb = getHPb(T0Pb)
	#Initial guesses
	Hw = [H0W]*N
	Hpb = [H0Pb]*N
	Tw = [T0W]*N
	Tpb = [T0Pb]*N
	
	HWmin = 16000
	HWmax = 3530000
	
	tempsW = [Tw.copy()]
	tempsPb = [Tpb.copy()]
	for j in range(cycles):
		lastDeltaTw = 0
		lastDeltaTpb = 0
		for i in range(1, N):
			# Provide an estimation for the temperature (will not
			# change the solution, just allow it to converge faster)
			T = Tw[i] + lastDeltaTw
			Hw[i] = Hw[i-1] + dHWdz(Tpb[i], T, Hw[i])*dz
			# Cannot allow enthalpy outside of sampled values
			Hw[i] = np.max((HWmin, np.min((HWmax, Hw[i]))))
			# Record how much T changed, the assumtion is that the
			# next T value will change by roughly the same amount
			lastTw = Tw[i]
			Tw[i] = getTW(Hw[i])
			lastDeltaTw = Tw[i] - lastTw
		for i in range(N-2, -1, -1):
			T = Tpb[i+1] + lastDeltaTpb
			Hpb[i] = Hpb[i+1] - dHPbdz(T, Tw[i+1], Hw[i+1])*dz
			lastTpb = Tpb[i]
			Tpb[i] = getTPb(Hpb[i])
			lastDeltaTpb = Tpb[i] - lastTpb
		if plotIntermediates:
			tempsW += [Tw.copy()]
			tempsPb += [Tpb.copy()]
		if prgBar:
			if j%(cycles/19) <= 1:
				print('=', end='', flush=True)
	if prgBar:
		print(']')

	if printSol:
		printSolutionData(Hw, Hpb)
	
	z = np.linspace(0, h, N)
	if plotIntermediates:
		for j in range(len(tempsW)-1):
			plt.plot(z, tempsPb[j])
			plt.plot(z, tempsW[j])
	if printSol:
		plt.plot(z, Tpb, label="Final lead temp [C]")
		plt.plot(z, Tw, label="Final water temp [C]")
		# plt.plot(z, Hpb, label="Final lead enthalpy [J/kg]")
		# plt.plot(z, Hw, label="Final water enthalpy [J/kg]")
		# Qw = [mDotW*(H - Hw[0])*1e-6 for H in Hw]
		# Qpb = [mDotPb*(H - Hpb[0])*1e-6 for H in Hpb]
		# plt.plot(z, Qpb, label="Final lead energy flow [MW]")
		# plt.plot(z, Qw, label="Final water energy flow [MW]")
		plt.legend()
		plt.xlabel("z")
		plt.ylabel("Temperature [C]")
		plt.title("Water inflow at z=0, lead inflow at z=h=" + str(h))
		plt.show()
	return Tw[N-1], Hw, Hpb, Tw, Tpb

def checkSolution(Hw, Hpb):
	"""Checks how well the given solution satisfied the differential equaitons. Returns a 
	number [W] that is a measure of well the functions satisfy the differential 
	equations, lower is better."""
	Tw = [getTW(H) for H in Hw]
	Tpb = [getTPb(H) for H in Hpb]
	discrepancyQW = []
	discrepancyQPb = []
	for i in range(1, N):
		deltaHW = Hw[i] - Hw[i-1]
		discrepancy = deltaHW - dHWdz(Tpb[i], Tw[i], Hw[i])*dz
		discrepancyQW.append(getQW(discrepancy)*1e6)
	for i in range(1, N):
		deltaHPb = Hpb[i] - Hpb[i-1]
		discrepancy = deltaHPb - dHPbdz(Tpb[i], Tw[i], Hw[i])*dz
		discrepancyQPb.append(getQPb(discrepancy)*1e6)
	if plotError:
		z = np.linspace(0, h, N-1)
		plt.plot(z, discrepancyQW, label="water")
		plt.plot(z, discrepancyQPb, label="lead")
		plt.legend()
		plt.show()
	return (np.sum([abs(d) for d in discrepancyQW]) + np.sum([abs(d) for d in discrepancyQPb]))/N

def solutionIsValid(Hw, Hpb, Tw, Tpb):
	tol = 0.001
	if Aw <= 0 or Apb <= 0 or ro*ro/(R*R)*dldz*n > 0.5:
		return False, 'A'
	if abs(Tw[-1]/T1W - 1) > tol:
		return False, 'B'
	if abs(Tpb[0]/T1Pb - 1) > tol:
		return False, 'C'
	if uW < 1 or uW > 3 or uPb < -3 or uPb > -1:
		return False, 'D'
	return True, None

def searchFordxdz(prgBar = progressBar, **kwargs):
	global dxdz
	iterations = -1
	tol = 1e-10
	argIter = kwargs.get('iterations')
	argTol = kwargs.get('tol')
	if argIter != None:
		iterations = argIter
	if argTol != None:
		tol = argTol
	i = 1
	while True:
		updateConstants()
		TW, Hw, Hpb, Tw, Tpb = simulate(False, prgBar)
		diff = TW - T1W
		dxdz -= diff/15
		print(diff)
		if iterations == -1:
			if abs(diff) < tol or i > 30:
				break
		else:
			if i >= iterations:
				break
		i += 1
	return solutionIsValid(Hw, Hpb, Tw, Tpb)

def searchFordxdzSmart(prgBar = progressBar, **kwargs):
	global dxdz
	iterations = -1
	tol = 1e-10
	argIter = kwargs.get('iterations')
	argTol = kwargs.get('tol')
	if argIter != None:
		iterations = argIter
	if argTol != None:
		tol = argTol
	i = 1
	# Desired specific enthalpy of water at outflow
	H1W = getHW(T1W)
	while True:
		updateConstants()
		_, Hw, Hpb, Tw, Tpb = simulate(False, prgBar)
		deltaZ = (Hw[-1] - H1W)/dHWdz(Tpb[-1], Tw[-1], Hw[-1])
		diff = deltaZ*dldz/h
		if (diff > dxdz/2):
			dxdz = dxdz*0.9
		else:
			dxdz -= diff
		print(dHWdz(Tpb[-1], Tw[-1], Hw[-1]))
		print(deltaZ)
		print(diff)
		if iterations == -1:
			if abs(diff) < tol or i > 30:
				break
		else:
			if i >= iterations:
				break
		i += 1
	return solutionIsValid(Hw, Hpb, Tw, Tpb)

def searchFordxdzNewton(prgBar = progressBar, **kwargs):
	global dxdz
	iterations = -1
	tol = 1e-10
	argIter = kwargs.get('iterations')
	argTol = kwargs.get('tol')
	if argIter != None:
		iterations = argIter
	if argTol != None:
		tol = argTol
	i = 1
	# Desired specific enthalpy of water at outflow
	H1W = getHW(T1W)
	while True:
		updateConstants()
		_, Hw, Hpb, Tw, Tpb = simulate(False, prgBar)
		dxdz += 0.0001
		updateConstants()
		_, dHw, dHpb, dTw, dTpb = simulate(False, prgBar)
		dxdz -= 0.0001
		derivativeAprx = (dHw[-1]-Hw[-1])/0.0001
		diff = (Hw[-1] - H1W)/derivativeAprx
		if (diff > dxdz/2):
			dxdz = dxdz*0.9
		else:
			dxdz -= diff
		# dxdz -= diff
		print(diff)
		if iterations == -1:
			if abs(diff) < tol or i > 30:
				break
		else:
			if i >= iterations:
				break
		i += 1
	return solutionIsValid(Hw, Hpb, Tw, Tpb)

def parameterAnalysis():
	global D, di, n
	dim = 5
	lst = []
	vn = []
	vdi = []
	vD = []
	vval = []
	verr = []
	for i in range(dim):
		n = 14 + i*4
		for j in range(dim):
			di = 8*0.001 + j*0.002
			for k in range(dim):
				D = 0.35 + k*0.2/(dim - 1)
				if i == 0 and j == 0 and k == 0:
					iterations = -1#10
				else:
					iterations = -1#5
				valid, errType = searchFordxdzNewton(True, iterations=iterations, tol=1e-5)
				progress = round((i*dim**2 + j*dim + k)*100/(dim**3), 2)
				print(progress, "%")
				lst.append((n, di, D, valid, errType))
				vn.append(n)
				vdi.append(di)
				vD.append(D)
				vval.append(valid)
				verr.append(errType)
	print("-- RAW DATA OUTPUT --")	
	print(lst)
	fig = plt.figure(figsize=(6, 6))
	ax3 = fig.add_subplot(111, projection='3d')
	plt.tight_layout()
	col = ['g' if b else 'r' for b in vval]
	size = [500 if b else 50 for b in vval]
	shape = ['.' if e==None else '$'+e+'$' for e in verr]
	for n_, di_, D_, col_, size_, shape_ in zip(vn, vdi, vD, col, size, shape):
		ax3.scatter([n_], [di_], [D_], c=col_, s=size_, alpha=0.5, marker=shape_)
	plt.xlabel('n')
	plt.ylabel('di')
	ax3.set_zlabel('D')
	plt.show()




# valid, type = searchFordxdzNewton(tolerance=1e-10)
# print(dxdz)
# updateConstants()
# simulate(True)
# print(valid)

parameterAnalysis()

# Hinterval = np.linspace(HWsamp[0], HWsamp[len(HWsamp)-1], 1000)
# PW = [getReff(H, 540) for H in Hinterval]
# plt.plot(Hinterval, PW)
# plt.show()

# plotEnthalpyData()
