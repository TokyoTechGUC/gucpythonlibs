#!/gs/bs/tga-guc-lab/dependencies/dependencies_intel/conda/envs/guconda/bin/python
import mpi4py.MPI
import numpy as np
import netCDF4 as nc
from wrf import getvar
from glob import glob
from wrf import ALL_TIMES
import os
import wrf
from scipy import optimize
import os.path
from os import path

def heat_index_wrf(ifil):
	outfile = ifil.replace('.nc','')+'_heatparams.nc'
	if path.exists(outfile):
		print(f'warning: {outfile} already exists. Skipping')
		return 0
	print(f"Calculating climate stresses for {ifil}")
	df = nc.Dataset(ifil,'r')
	test = df.variables['T2'][:]
	UTCI_approx,MRT = UTCI(df)
	WBGT_approx,Tg,Tnwb,Tpwb = WBGT(df)
	UTCI_approx = np.asarray(UTCI_approx)
	MRT = np.asarray(MRT)
	WBGT_approx = np.asarray(WBGT_approx)
	Tg = np.asarray(Tg)
	Tnwb = np.asarray(Tnwb)
	Tpwb = np.asarray(Tpwb)
	HI = np.asarray(HeatIndex(df))
	appTemp = np.asarray(apparentTemp(df))
	effTemp = np.asarray(effectiveTemp(df))
	humidex = np.asarray(humidex_calc(df))
	discomInd = np.asarray(discomInd_calc(df))
	print(f"Calculation for {ifil} finished. Now writing to {outfile}.")
	with nc.Dataset(outfile,'w',format='NETCDF4') as dst:
		dst.setncatts(df.__dict__)
		for name,dimension in df.dimensions.items():
			if name in ['Time','DateStrLen','west_east','south_north','west_east_stag','south_north_stag']:
				dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
		for name,variable in df.variables.items():
			if name in ['XTIME']:
				x = dst.createVariable(name, variable.datatype, variable.dimensions,zlib=True,complevel=5,shuffle=True)
				dst[name].setncatts(df[name].__dict__)
				dst[name][:] = df[name][:]
			if name in ['XLAT']:
				x = dst.createVariable(name, variable.datatype, variable.dimensions,zlib=True,complevel=5,shuffle=True)
				dst[name].setncatts(df[name].__dict__)
				dst[name][:,:,:] = df[name][:,:,:]
			if name in ['XLONG']:
				x = dst.createVariable(name, variable.datatype, variable.dimensions,zlib=True,complevel=5,shuffle=True)
				dst[name].setncatts(df[name].__dict__)
				dst[name][:,:,:] = df[name][:,:,:]
			if name in ['T2']:
				x = dst.createVariable(name, variable.datatype, variable.dimensions,zlib=True,complevel=5,shuffle=True)
				dst[name].setncatts(df[name].__dict__)
				dst[name][:,:,:] = df[name][:,:,:]
				x = dst.createVariable('Tpwb',variable.datatype,variable.dimensions,zlib=True,complevel=5,shuffle=True)
				dst['Tpwb'].setncatts(df[name].__dict__)
				dst['Tpwb'].setncattr('description','Psychrometric Wet-Bulb Temperature')
				dst['Tpwb'].setncattr('units','K')
				dst['Tpwb'][:,:,:] = Tpwb[:,:,:]
				del(Tpwb)
				x = dst.createVariable('Tnwb',variable.datatype,variable.dimensions,zlib=True,complevel=5,shuffle=True)
				dst['Tnwb'].setncatts(df[name].__dict__)
				dst['Tnwb'].setncattr('description','Natural Wet-Bulb Temperature')
				dst['Tnwb'].setncattr('units','K')
				dst['Tnwb'][:,:,:] = Tnwb[:,:,:]
				del(Tnwb)
				x = dst.createVariable('Tg',variable.datatype,variable.dimensions,zlib=True,complevel=5,shuffle=True)
				dst['Tg'].setncatts(df[name].__dict__)
				dst['Tg'].setncattr('description','Globe Temperature')
				dst['Tg'].setncattr('units','K')
				dst['Tg'][:,:,:] = Tg[:,:,:]
				del(Tg)
				x = dst.createVariable('WBGT',variable.datatype,variable.dimensions,zlib=True,complevel=5,shuffle=True)
				dst['WBGT'].setncatts(df[name].__dict__)
				dst['WBGT'].setncattr('description','Wet-Bulb Globe Temperature')
				dst['WBGT'].setncattr('units','K')
				dst['WBGT'][:,:,:] = WBGT_approx[:,:,:]
				del(WBGT_approx)
				x = dst.createVariable('Tmrt',variable.datatype,variable.dimensions,zlib=True,complevel=5,shuffle=True)
				dst['Tmrt'].setncatts(df[name].__dict__)
				dst['Tmrt'].setncattr('description','Mean radiant temperature 2m')
				dst['Tmrt'].setncattr('units','K')
				dst['Tmrt'][:,:,:] = MRT[:,:,:]
				del(MRT)
				x = dst.createVariable('UTCI',variable.datatype,variable.dimensions,zlib=True,complevel=5,shuffle=True)
				dst['UTCI'].setncatts(df[name].__dict__)
				dst['UTCI'].setncattr('description','Universal Thermal Comfort Index 2m')
				dst['UTCI'].setncattr('units','degC')
				dst['UTCI'][:,:,:] = UTCI_approx[:,:,:]
				del(UTCI_approx)
				x = dst.createVariable('RH2',variable.datatype,variable.dimensions,zlib=True,complevel=5,shuffle=True)
				dst['RH2'].setncatts(df[name].__dict__)
				dst['RH2'].setncattr('description','Relative Humidity 2m')
				dst['RH2'].setncattr('units','%')
				RH = np.asarray(getvar(df,'rh2',timeidx=ALL_TIMES)[:])
				if len(RH.shape) < 3:
					RH = RH.reshape(1,RH.shape[0],RH.shape[1])
				dst['RH2'][:,:,:] = RH
				x = dst.createVariable('WNDSPD',variable.datatype,variable.dimensions,zlib=True,complevel=5,shuffle=True)
				dst['WNDSPD'].setncatts(df[name].__dict__)
				dst['WNDSPD'].setncattr('description','Wind speed 10m')
				dst['WNDSPD'].setncattr('units','m/s')
				dst['WNDSPD'][:,:,:] = np.sqrt(df.variables['U10'][:]**2.0+df.variables['V10'][:]**2.0)[:,:,:]
				x = dst.createVariable('HI',variable.datatype,variable.dimensions,zlib=True,complevel=5,shuffle=True)
				dst['HI'].setncatts(df[name].__dict__)
				dst['HI'].setncattr('description','Heat Index (T2,RH2)')
				dst['HI'].setncattr('units','degF')
				dst['HI'][:,:,:] = HI[:,:,:]
				x = dst.createVariable('appTemp',variable.datatype,variable.dimensions,zlib=True,complevel=5,shuffle=True)
				dst['appTemp'].setncatts(df[name].__dict__)
				dst['appTemp'].setncattr('description','Apparent Temperature (T2,RH2,Wind10)')
				dst['appTemp'].setncattr('units','degC')
				dst['appTemp'][:,:,:] = appTemp[:,:,:]
				x = dst.createVariable('effTemp',variable.datatype,variable.dimensions,zlib=True,complevel=5,shuffle=True)
				dst['effTemp'].setncatts(df[name].__dict__)
				dst['effTemp'].setncattr('description','Effective Temperature (T2,RH2,Wind10)')
				dst['effTemp'].setncattr('units','degC')
				dst['effTemp'][:,:,:] = effTemp[:,:,:]
				x = dst.createVariable('humidex',variable.datatype,variable.dimensions,zlib=True,complevel=5,shuffle=True)
				dst['humidex'].setncatts(df[name].__dict__)
				dst['humidex'].setncattr('description','Humidex (T2,RH2)')
				dst['humidex'].setncattr('units','degC')
				dst['humidex'][:,:,:] = humidex[:,:,:]
				x = dst.createVariable('discomInd',variable.datatype,variable.dimensions,zlib=True,complevel=5,shuffle=True)
				dst['discomInd'].setncatts(df[name].__dict__)
				dst['discomInd'].setncattr('description','Discomfort Index (T2,RH2)')
				dst['discomInd'].setncattr('units','degC')
				dst['discomInd'][:,:,:] = discomInd[:,:,:]
	return print(f"{outfile} created.")

def discomInd_calc(df):
	tas = df.variables['T2'][:]-273.15	
	hurs = np.asarray(getvar(df,'rh2',timeidx=ALL_TIMES)[:])
	return tas - 0.55*(1.0-0.01*hurs)*(tas-14.5)

def humidex_calc(df):
	tas = df.variables['T2'][:]-273.15	
	hurs = np.asarray(getvar(df,'rh2',timeidx=ALL_TIMES)[:])
	if len(hurs.shape) < 3:
		hurs = hurs.reshape(1,hurs.shape[0],hurs.shape[1])
	vp = tashurs2vap(tas,hurs)
	return tas+5.0/9.0*(vp-10.0)
	
	
def effectiveTemp(df):
	tas = df.variables['T2'][:]-273.15	
	wind = np.sqrt(df.variables['U10'][:]**2.0+df.variables['V10'][:]**2.0)
	hurs = np.asarray(getvar(df,'rh2',timeidx=ALL_TIMES)[:])
	if len(hurs.shape) < 3:
		hurs = hurs.reshape(1,hurs.shape[0],hurs.shape[1])
	return  37 - (37-tas)/(0.68 - 0.0014*hurs + 1/(1.76+ 1.4*(wind**0.75))) - 0.29 * tas * (1 - 0.01*hurs)
	
def apparentTemp(df):
	#Note that this returns degC
	tas = df.variables['T2'][:]-273.15	
	wind = np.sqrt(df.variables['U10'][:]**2.0+df.variables['V10'][:]**2.0)
	hurs = np.asarray(getvar(df,'rh2',timeidx=ALL_TIMES)[:])
	if len(hurs.shape) < 3:
		hurs = hurs.reshape(1,hurs.shape[0],hurs.shape[1])
	vp = tashurs2vap(tas,hurs)
	return tas + 0.33*vp - 0.7*wind - 4.0

def tashurs2vap(tas,hurs):
	c1 = 0.06107
	a1 = 17.368
	b1 = 2388.3
	c2 = 0.06108
	a2 = 17.856
	b2 = 2455.2
	T0 = 0.0
	tas = tas*10.0
	iceMask = tas<T0
	waterMask = tas>=T0
	vp = np.empty_like(tas)
	vp[waterMask] = hurs[waterMask]*c1*np.exp((a1*tas[waterMask])/(b1+tas[waterMask]))
	vp[iceMask] = hurs[iceMask] * c2 *np.exp((a2*tas[iceMask])/(b2+tas[iceMask]))
	return vp
	
def UTCI(df):
	# ehPa <- Water Vapour Pressure in hPa; RH <- Relative humidity in percent
	'''
	Excerpt from the original script (ensure units are consistent):
	      print '(a,f10.2)','  The air temperature in degC is          : ',Ta
	      print '(a,f10.2)','  The mean radiant temperature in degC is : ',tmrt
	      print '(a,f10.2)','  The wind velocity (10 m) in m/s is   : ',va10m
	      print '(a,f10.2)','  The water vapour pressure in hPa is   : ',ehPa
	      print '(a,f10.2)','  The relative humidity in percent is   : ',RH
	'''
	Ta = df.variables['T2'][:]
	va = np.sqrt(df.variables['U10'][:]**2.0+df.variables['V10'][:]**2.0)
	# Calculate MRT
	Tmrt = MRT(df)

	# Calculate UCTI
	D_Tmrt = Tmrt-Ta #This is in degC
	RH = np.asarray(getvar(df,'rh2',timeidx=ALL_TIMES)[:])
	if len(RH.shape) < 3:
		RH = RH.reshape(1,RH.shape[0],RH.shape[1])
	ehPa = es_calc(Ta)*RH/100.0
	Ta = Ta - 273.15
	Pa = ehPa/10.0 #This is vapour pressure in kPa
	UTCI_approx=Ta+(6.07562052E-01)+(-2.27712343E-02)*Ta+(8.06470249E-04)*Ta*Ta+(-1.54271372E-04)*Ta*Ta*Ta+(-3.24651735E-06)*Ta*Ta*Ta*Ta+(7.32602852E-08)*Ta*Ta*Ta*Ta*Ta+(1.35959073E-09)*Ta*Ta*Ta*Ta*Ta*Ta+(-2.25836520E+00)*va+(8.80326035E-02)*Ta*va+(2.16844454E-03)*Ta*Ta*va+(-1.53347087E-05)*Ta*Ta*Ta*va+(-5.72983704E-07)*Ta*Ta*Ta*Ta*va+(-2.55090145E-09)*Ta*Ta*Ta*Ta*Ta*va+(-7.51269505E-01)*va*va+(-4.08350271E-03)*Ta*va*va+(-5.21670675E-05)*Ta*Ta*va*va+(1.94544667E-06)*Ta*Ta*Ta*va*va+(1.14099531E-08)*Ta*Ta*Ta*Ta*va*va+(1.58137256E-01)*va*va*va+(-6.57263143E-05)*Ta*va*va*va+(2.22697524E-07)*Ta*Ta*va*va*va+(-4.16117031E-08)*Ta*Ta*Ta*va*va*va+(-1.27762753E-02)*va*va*va*va+(9.66891875E-06)*Ta*va*va*va*va+(2.52785852E-09)*Ta*Ta*va*va*va*va+(4.56306672E-04)*va*va*va*va*va+(-1.74202546E-07)*Ta*va*va*va*va*va+(-5.91491269E-06)*va*va*va*va*va*va+(3.98374029E-01)*D_Tmrt+(1.83945314E-04)*Ta*D_Tmrt+(-1.73754510E-04)*Ta*Ta*D_Tmrt+(-7.60781159E-07)*Ta*Ta*Ta*D_Tmrt+(3.77830287E-08)*Ta*Ta*Ta*Ta*D_Tmrt+(5.43079673E-10)*Ta*Ta*Ta*Ta*Ta*D_Tmrt+(-2.00518269E-02)*va*D_Tmrt+(8.92859837E-04)*Ta*va*D_Tmrt+(3.45433048E-06)*Ta*Ta*va*D_Tmrt+(-3.77925774E-07)*Ta*Ta*Ta*va*D_Tmrt+(-1.69699377E-09)*Ta*Ta*Ta*Ta*va*D_Tmrt+(1.69992415E-04)*va*va*D_Tmrt+(-4.99204314E-05)*Ta*va*va*D_Tmrt+(2.47417178E-07)*Ta*Ta*va*va*D_Tmrt+(1.07596466E-08)*Ta*Ta*Ta*va*va*D_Tmrt+(8.49242932E-05)*va*va*va*D_Tmrt+(1.35191328E-06)*Ta*va*va*va*D_Tmrt+(-6.21531254E-09)*Ta*Ta*va*va*va*D_Tmrt+(-4.99410301E-06)*va*va*va*va*D_Tmrt+(-1.89489258E-08)*Ta*va*va*va*va*D_Tmrt+(8.15300114E-08)*va*va*va*va*va*D_Tmrt+(7.55043090E-04)*D_Tmrt*D_Tmrt+(-5.65095215E-05)*Ta*D_Tmrt*D_Tmrt+(-4.52166564E-07)*Ta*Ta*D_Tmrt*D_Tmrt+(2.46688878E-08)*Ta*Ta*Ta*D_Tmrt*D_Tmrt+(2.42674348E-10)*Ta*Ta*Ta*Ta*D_Tmrt*D_Tmrt+(1.54547250E-04)*va*D_Tmrt*D_Tmrt+(5.24110970E-06)*Ta*va*D_Tmrt*D_Tmrt+(-8.75874982E-08)*Ta*Ta*va*D_Tmrt*D_Tmrt+(-1.50743064E-09)*Ta*Ta*Ta*va*D_Tmrt*D_Tmrt+(-1.56236307E-05)*va*va*D_Tmrt*D_Tmrt+(-1.33895614E-07)*Ta*va*va*D_Tmrt*D_Tmrt+(2.49709824E-09)*Ta*Ta*va*va*D_Tmrt*D_Tmrt+(6.51711721E-07)*va*va*va*D_Tmrt*D_Tmrt+(1.94960053E-09)*Ta*va*va*va*D_Tmrt*D_Tmrt+(-1.00361113E-08)*va*va*va*va*D_Tmrt*D_Tmrt+(-1.21206673E-05)*D_Tmrt*D_Tmrt*D_Tmrt+(-2.18203660E-07)*Ta*D_Tmrt*D_Tmrt*D_Tmrt+(7.51269482E-09)*Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt+(9.79063848E-11)*Ta*Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt+(1.25006734E-06)*va*D_Tmrt*D_Tmrt*D_Tmrt+(-1.81584736E-09)*Ta*va*D_Tmrt*D_Tmrt*D_Tmrt+(-3.52197671E-10)*Ta*Ta*va*D_Tmrt*D_Tmrt*D_Tmrt+(-3.36514630E-08)*va*va*D_Tmrt*D_Tmrt*D_Tmrt+(1.35908359E-10)*Ta*va*va*D_Tmrt*D_Tmrt*D_Tmrt+(4.17032620E-10)*va*va*va*D_Tmrt*D_Tmrt*D_Tmrt+(-1.30369025E-09)*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt+(4.13908461E-10)*Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt+(9.22652254E-12)*Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt+(-5.08220384E-09)*va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt+(-2.24730961E-11)*Ta*va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt+(1.17139133E-10)*va*va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt+(6.62154879E-10)*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt+(4.03863260E-13)*Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt+(1.95087203E-12)*va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt+(-4.73602469E-12)*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt+(5.12733497E+00)*Pa+(-3.12788561E-01)*Ta*Pa+(-1.96701861E-02)*Ta*Ta*Pa+(9.99690870E-04)*Ta*Ta*Ta*Pa+(9.51738512E-06)*Ta*Ta*Ta*Ta*Pa+(-4.66426341E-07)*Ta*Ta*Ta*Ta*Ta*Pa+(5.48050612E-01)*va*Pa+(-3.30552823E-03)*Ta*va*Pa+(-1.64119440E-03)*Ta*Ta*va*Pa+(-5.16670694E-06)*Ta*Ta*Ta*va*Pa+(9.52692432E-07)*Ta*Ta*Ta*Ta*va*Pa+(-4.29223622E-02)*va*va*Pa+(5.00845667E-03)*Ta*va*va*Pa+(1.00601257E-06)*Ta*Ta*va*va*Pa+(-1.81748644E-06)*Ta*Ta*Ta*va*va*Pa+(-1.25813502E-03)*va*va*va*Pa+(-1.79330391E-04)*Ta*va*va*va*Pa+(2.34994441E-06)*Ta*Ta*va*va*va*Pa+(1.29735808E-04)*va*va*va*va*Pa+(1.29064870E-06)*Ta*va*va*va*va*Pa+(-2.28558686E-06)*va*va*va*va*va*Pa+(-3.69476348E-02)*D_Tmrt*Pa+(1.62325322E-03)*Ta*D_Tmrt*Pa+(-3.14279680E-05)*Ta*Ta*D_Tmrt*Pa+(2.59835559E-06)*Ta*Ta*Ta*D_Tmrt*Pa+(-4.77136523E-08)*Ta*Ta*Ta*Ta*D_Tmrt*Pa+(8.64203390E-03)*va*D_Tmrt*Pa+(-6.87405181E-04)*Ta*va*D_Tmrt*Pa+(-9.13863872E-06)*Ta*Ta*va*D_Tmrt*Pa+(5.15916806E-07)*Ta*Ta*Ta*va*D_Tmrt*Pa+(-3.59217476E-05)*va*va*D_Tmrt*Pa+(3.28696511E-05)*Ta*va*va*D_Tmrt*Pa+(-7.10542454E-07)*Ta*Ta*va*va*D_Tmrt*Pa+(-1.24382300E-05)*va*va*va*D_Tmrt*Pa+(-7.38584400E-09)*Ta*va*va*va*D_Tmrt*Pa+(2.20609296E-07)*va*va*va*va*D_Tmrt*Pa+(-7.32469180E-04)*D_Tmrt*D_Tmrt*Pa+(-1.87381964E-05)*Ta*D_Tmrt*D_Tmrt*Pa+(4.80925239E-06)*Ta*Ta*D_Tmrt*D_Tmrt*Pa+(-8.75492040E-08)*Ta*Ta*Ta*D_Tmrt*D_Tmrt*Pa+(2.77862930E-05)*va*D_Tmrt*D_Tmrt*Pa+(-5.06004592E-06)*Ta*va*D_Tmrt*D_Tmrt*Pa+(1.14325367E-07)*Ta*Ta*va*D_Tmrt*D_Tmrt*Pa+(2.53016723E-06)*va*va*D_Tmrt*D_Tmrt*Pa+(-1.72857035E-08)*Ta*va*va*D_Tmrt*D_Tmrt*Pa+(-3.95079398E-08)*va*va*va*D_Tmrt*D_Tmrt*Pa+(-3.59413173E-07)*D_Tmrt*D_Tmrt*D_Tmrt*Pa+(7.04388046E-07)*Ta*D_Tmrt*D_Tmrt*D_Tmrt*Pa+(-1.89309167E-08)*Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt*Pa+(-4.79768731E-07)*va*D_Tmrt*D_Tmrt*D_Tmrt*Pa+(7.96079978E-09)*Ta*va*D_Tmrt*D_Tmrt*D_Tmrt*Pa+(1.62897058E-09)*va*va*D_Tmrt*D_Tmrt*D_Tmrt*Pa+(3.94367674E-08)*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa+(-1.18566247E-09)*Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa+(3.34678041E-10)*va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa+(-1.15606447E-10)*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa+(-2.80626406E+00)*Pa*Pa+(5.48712484E-01)*Ta*Pa*Pa+(-3.99428410E-03)*Ta*Ta*Pa*Pa+(-9.54009191E-04)*Ta*Ta*Ta*Pa*Pa+(1.93090978E-05)*Ta*Ta*Ta*Ta*Pa*Pa+(-3.08806365E-01)*va*Pa*Pa+(1.16952364E-02)*Ta*va*Pa*Pa+(4.95271903E-04)*Ta*Ta*va*Pa*Pa+(-1.90710882E-05)*Ta*Ta*Ta*va*Pa*Pa+(2.10787756E-03)*va*va*Pa*Pa+(-6.98445738E-04)*Ta*va*va*Pa*Pa+(2.30109073E-05)*Ta*Ta*va*va*Pa*Pa+(4.17856590E-04)*va*va*va*Pa*Pa+(-1.27043871E-05)*Ta*va*va*va*Pa*Pa+(-3.04620472E-06)*va*va*va*va*Pa*Pa+(5.14507424E-02)*D_Tmrt*Pa*Pa+(-4.32510997E-03)*Ta*D_Tmrt*Pa*Pa+(8.99281156E-05)*Ta*Ta*D_Tmrt*Pa*Pa+(-7.14663943E-07)*Ta*Ta*Ta*D_Tmrt*Pa*Pa+(-2.66016305E-04)*va*D_Tmrt*Pa*Pa+(2.63789586E-04)*Ta*va*D_Tmrt*Pa*Pa+(-7.01199003E-06)*Ta*Ta*va*D_Tmrt*Pa*Pa+(-1.06823306E-04)*va*va*D_Tmrt*Pa*Pa+(3.61341136E-06)*Ta*va*va*D_Tmrt*Pa*Pa+(2.29748967E-07)*va*va*va*D_Tmrt*Pa*Pa+(3.04788893E-04)*D_Tmrt*D_Tmrt*Pa*Pa+(-6.42070836E-05)*Ta*D_Tmrt*D_Tmrt*Pa*Pa+(1.16257971E-06)*Ta*Ta*D_Tmrt*D_Tmrt*Pa*Pa+(7.68023384E-06)*va*D_Tmrt*D_Tmrt*Pa*Pa+(-5.47446896E-07)*Ta*va*D_Tmrt*D_Tmrt*Pa*Pa+(-3.59937910E-08)*va*va*D_Tmrt*D_Tmrt*Pa*Pa+(-4.36497725E-06)*D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa+(1.68737969E-07)*Ta*D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa+(2.67489271E-08)*va*D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa+(3.23926897E-09)*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa+(-3.53874123E-02)*Pa*Pa*Pa+(-2.21201190E-01)*Ta*Pa*Pa*Pa+(1.55126038E-02)*Ta*Ta*Pa*Pa*Pa+(-2.63917279E-04)*Ta*Ta*Ta*Pa*Pa*Pa+(4.53433455E-02)*va*Pa*Pa*Pa+(-4.32943862E-03)*Ta*va*Pa*Pa*Pa+(1.45389826E-04)*Ta*Ta*va*Pa*Pa*Pa+(2.17508610E-04)*va*va*Pa*Pa*Pa+(-6.66724702E-05)*Ta*va*va*Pa*Pa*Pa+(3.33217140E-05)*va*va*va*Pa*Pa*Pa+(-2.26921615E-03)*D_Tmrt*Pa*Pa*Pa+(3.80261982E-04)*Ta*D_Tmrt*Pa*Pa*Pa+(-5.45314314E-09)*Ta*Ta*D_Tmrt*Pa*Pa*Pa+(-7.96355448E-04)*va*D_Tmrt*Pa*Pa*Pa+(2.53458034E-05)*Ta*va*D_Tmrt*Pa*Pa*Pa+(-6.31223658E-06)*va*va*D_Tmrt*Pa*Pa*Pa+(3.02122035E-04)*D_Tmrt*D_Tmrt*Pa*Pa*Pa+(-4.77403547E-06)*Ta*D_Tmrt*D_Tmrt*Pa*Pa*Pa+(1.73825715E-06)*va*D_Tmrt*D_Tmrt*Pa*Pa*Pa+(-4.09087898E-07)*D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa*Pa+(6.14155345E-01)*Pa*Pa*Pa*Pa+(-6.16755931E-02)*Ta*Pa*Pa*Pa*Pa+(1.33374846E-03)*Ta*Ta*Pa*Pa*Pa*Pa+(3.55375387E-03)*va*Pa*Pa*Pa*Pa+(-5.13027851E-04)*Ta*va*Pa*Pa*Pa*Pa+(1.02449757E-04)*va*va*Pa*Pa*Pa*Pa+(-1.48526421E-03)*D_Tmrt*Pa*Pa*Pa*Pa+(-4.11469183E-05)*Ta*D_Tmrt*Pa*Pa*Pa*Pa+(-6.80434415E-06)*va*D_Tmrt*Pa*Pa*Pa*Pa+(-9.77675906E-06)*D_Tmrt*D_Tmrt*Pa*Pa*Pa*Pa+(8.82773108E-02)*Pa*Pa*Pa*Pa*Pa+(-3.01859306E-03)*Ta*Pa*Pa*Pa*Pa*Pa+(1.04452989E-03)*va*Pa*Pa*Pa*Pa*Pa+(2.47090539E-04)*D_Tmrt*Pa*Pa*Pa*Pa*Pa+(1.48348065E-03)*Pa*Pa*Pa*Pa*Pa*Pa
	'''
	Thresholds:
		Ta .ge. -50 .and. TA .le. 50
		Tmrt-Ta .ge. -30 .and. Tmrt-Ta .le. 70
		va10m .ge. 0.5 .and. va10m .le. 30
		ehPa .le. 50 .and. RH .le. 100 .and. ehPa .ge. 0 .and. RH .ge. 0
		Note: ehPa <- Water Vapour Pressure in hPa; RH <- Relative humidity in percent
	'''
	UTCI_approx[(Ta<=-50.)&(Ta>=50.0)] = np.nan
	UTCI_approx[(D_Tmrt<=-30.0)&(D_Tmrt>=70.0)] = np.nan
	UTCI_approx[(va<=0.5)&(va>=30.0)] = np.nan
	UTCI_approx[(ehPa>=50.0)&(RH>=100.0)&(ehPa<=0.0)&(RH<=0.0)] = np.nan
	land = df.variables['LANDMASK'][:]
	UTCI_approx[land==0]=np.nan
	return UTCI_approx,Tmrt

def es_calc(tk):
	g = np.asarray([-2.8365744E3,-6.028076559E3,1.954263612E1,-2.737830188E-2,1.6261698E-5,7.0229056E-10,-1.8680009E-13,2.715030],dtype=np.float64)
	es = 0.01*np.exp(g[7]*np.log(tk) + g[0]*tk**(-2.0)
		+g[1]*tk**(-1.0) + g[2]*tk**(0.0)
		+g[3]*tk**(1.0) + g[4]*tk**(2.0)
		+g[5]*tk**(3.0) + g[6]*tk**(4.0))
	return es

def MRT(df):
	'''
	Reference: https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/gdj3.102
	arccos(COSZEN) -- Zenith Angle
	SE = Solar elevation angle = 90.0 - arccos(COSZEN)
	fp = surface projection factor
	'''
	SE = 90.0 - np.degrees(np.arccos(df.variables['COSZEN'][:]))
	fp = 0.308*np.cos(np.radians(SE*(0.998-SE*SE/50000.0)))
	fa = 0.5
	emissp = 0.97
	sigma = 5.67E-8
	A_ir = 0.7
	SRAD = df.variables['SWDOWN'][:]*(df.variables['ALBEDO'][:]+(fp/fa))
	MRT = ((fa/sigma)*(df.variables['GLW'][:]+sigma*df.variables['EMISS'][:]*df.variables['TSK'][:]**4.0
		+(A_ir/emissp)*SRAD))**0.25
	return MRT
	
	### WBGT Calculation ####
	'''
	Depending on the time of the day determined by solar radiation,
	the algorithm for WBGT will change.
	'''

def HeatIndex(df):
	#Calculate the heat index from temperature and relative humidity.
	tas = df.variables['T2'][:] - 273.15
	tasf = tas*1.8+32.0
	hurs = np.asarray(getvar(df,'rh2',timeidx=ALL_TIMES)[:])
	if len(hurs.shape) < 3:
		hurs = hurs.reshape(1,hurs.shape[0],hurs.shape[1])
	a = -42.379
	b = 2.04901523
	c = 10.14333127
	d = -0.22475541
	e = -6.83783E-3 
	f = -5.481717E-2
	g = 1.22874E-3
	h = 8.5282E-4
	i = -1.99e-6
	result = a + b * tasf + c * hurs + d * tasf * hurs + e * tasf**2 + f * hurs**2 + g * tasf**2 * hurs + h * tasf * hurs**2 + i * tasf**2 * hurs**2
	i = (tasf<80.)
	result[i] = 0.0
	i = (hurs>85.0)&(tasf<=87.0)
	result[i] = result[i] + ((hurs[i]-85)/10) * ((87-tasf[i])/5)
	i = (hurs<13.0)&(tasf<=112.0)
	result[i] = result[i] - ((13-hurs[i])/4)*np.sqrt((17-np.abs(tasf[i]-95))/17)
	return result
	
def WBGT(df,MinWindSpeed=0.1,propDirect=0.8,irad=1.0):
	tas = df.variables['T2'][:]	
	wind = np.sqrt(df.variables['U10'][:]**2.0+df.variables['V10'][:]**2.0)
	radiation = df.variables['SWDOWN'][:]
	Pair = df.variables['PSFC'][:]*0.01
	RH = np.asarray(getvar(df,'rh2',timeidx=ALL_TIMES)[:])
	if len(RH.shape) < 3:
		RH = RH.reshape(1,RH.shape[0],RH.shape[1])
	dewp = np.asarray(getvar(df,'td2',timeidx=ALL_TIMES,units='K')[:])
	if len(dewp.shape) < 3:
		dewp = dewp.reshape(1,dewp.shape[0],dewp.shape[1])
	cza  = df.variables['COSZEN'][:]
	SurfAlbedo = df.variables['ALBEDO'][:]
	land = df.variables['LANDMASK'][:]
	return np.vectorize(wbgt_indiv)(land,tas,wind,radiation,Pair,RH,dewp,cza,SurfAlbedo,MinWindSpeed,propDirect,irad)

def wbgt_indiv(land,tas,wind,radiation,Pair,RH,dewp,cza,SurfAlbedo,MinWindSpeed,propDirect,irad):
	if land == 0:
		return np.nan,np.nan,np.nan,np.nan
	if radiation > 0.:
		Tg = globe_temp(tas,RH,Pair,wind,MinWindSpeed,radiation,propDirect,cza,SurfAlbedo)
		Tnwb = natglob_temp(tas,dewp,RH,Pair,wind,MinWindSpeed,radiation,propDirect,cza,irad,SurfAlbedo)
		Tpwb = np.empty_like(tas)
		WBGT = 0.7*Tnwb + 0.2*Tg + 0.1*tas
	else:
		Tg = np.empty_like(tas)
		Tnwb = np.empty_like(tas)
		Tpwb = Bernard(tas-273.15,dewp-273.15)
		WBGT = 0.67*Tpwb+0.33*tas
	return WBGT,Tg,Tnwb,Tpwb

def Bernard(tas,dewp):
	c1 = 6.106
	c2 = 17.27
	c3 = 237.3
	c4 = 1556.0
	c5 = 1.484
	c6 = 1010.0
	ed = c1*np.exp((c2*dewp)/(c3+dewp))
	def func_to_opt_Bernard(TT):
		return np.abs(c4*ed-c5*ed*TT-c4*c1*np.exp((c2*TT)/(c3+TT))+c5*c1*np.exp((c2*TT)/(c3+TT))*TT+c6*(tas-TT))
	lowbound = np.min(dewp)-2.5
	upbound = np.max(tas)+2.5
	if(lowbound>upbound):lowbound = upbound - 5.0
	Tpwb = optimize.minimize_scalar(func_to_opt_Bernard,method='bounded',bounds=(lowbound,upbound)).x
	return Tpwb+273.15

def globe_temp(Tair,relh,Pair,wind,MinWindSpeed,radiation,propDirect,cza,SurfAlbedo=0.4):
	stefanb = 0.000000056696
	cp = 1003.5
	mair = 28.97
	rgas = 8314.34
	rair = rgas/mair
	Pr = cp/(cp + (1.25*rair))
	emisglob = 0.95 #Emissivity
	albglob = 0.05 #Albedo
	diamglob = 0.05 #50mm Diameter globe
	emissfc = 0.999
	albsfc = SurfAlbedo
	Tsfc = np.copy(Tair)
	RH = relh*0.01
	def func_to_opt_tglob(Tglobe_prev):
		Tref = 0.5*(Tglobe_prev + Tair)
		h = h_sphere_in_air(Tref, Pair, wind, MinWindSpeed, diamglob)
		Tglobe = (0.5*(emis_atm(Tair,RH)*Tair**4.0+emissfc*Tsfc**4.0)-h/(emisglob*stefanb)*(Tglobe_prev-Tair)+radiation/(2.0*emisglob*stefanb)*(1.0-albglob)*(propDirect*(1.0/(2.0*cza)-1.0)+1.0+albsfc))**0.25
		return np.abs(Tglobe-Tglobe_prev)
	#Tg = optimize.brent(func_to_opt_tglob,brack=(np.min(Tair-2.0),np.max(Tair+35.0)))
	Tg = optimize.minimize_scalar(func_to_opt_tglob,method='bounded',bounds=(np.min(Tair-2.0),np.max(Tair+10.0))).x
	return Tg
	#Note tas must be in Kelvin, relh must be in decimals. cza is conside of zenith angle.

def natglob_temp(Tair,dewp,RH,Pair,wind,MinWindSpeed,radiation,propDirect,cza,irad=1.0,SurfAlbedo=0.4):
	stefanb = 0.000000056696
	cp = 1003.5
	mair = 28.97
	mh2o = 18.015
	rgas = 8314.34
	rair = rgas/mair
	ratio = cp*mair/mh2o
	Pr = cp/(cp + (1.25*rair))
	#Wick constants
	emiswick = 0.95
	albwick = 0.4
	diamwick = 0.007
	lenwick = 0.0254
	#Globe consants
	emisglobe = 0.95
	albglobe = 0.05
	diamglobe = 0.0508
	#Surface constants
	emissfc = 0.999
	albsfc = SurfAlbedo
	
	#Note dewp,Tair must be Kelvin. RH must be in frac.
	RH = RH * 0.01
	eair = RH*esat(Tair)
	emisatm = emis_atm(Tair,RH)
	Tsfc = np.copy(Tair)
	density = Pair*100.0/(rair*Tair)

	def func_to_opt_natglob(Twb_prev):
		Tref = 0.5*(Twb_prev+Tair)
		Fatm = stefanb*emiswick*(0.5*(emisatm*Tair**4.0+emissfc*Tsfc**4.0) - Twb_prev**4.0) + (1.0-albwick)*radiation*((1.0-propDirect)*(1+0.25*diamwick/lenwick)+((np.tan(np.arccos(cza))/3.1416)+0.25*diamwick/lenwick)*propDirect + albsfc)
		Sc = viscosity(Tair)/(density*diffusivity(Tref,Pair))
		h = h_cylinder_in_air(Twb_prev,Pair,wind,MinWindSpeed,diamwick)
		ewick = esat(Twb_prev)
		evap = h_evap(Twb_prev)
		Twb = Tair - evap/ratio*(ewick-eair)/(Pair-ewick)*(Pr/Sc)**0.56+Fatm/h*irad
		return np.abs(Twb-Twb_prev)
	#Tnwb = optimize.brent(func_to_opt_natglob,brack=(np.min(dewp-2.0),np.max(Tair+1.0)))
	lowbound = np.min(dewp)-2.5
	upbound = np.max(Tair)+2.5
	if(lowbound>upbound):lowbound = upbound - 5.0
	Tnwb = optimize.minimize_scalar(func_to_opt_natglob,method='bounded',bounds=(lowbound,upbound)).x
	return Tnwb

def diffusivity(Tk,Pair):
	pcrit13 = (36.4*218.0)**(1.0/3.0)
	tcrit512 = (132.0*647.3)**(5.0/12.0)
	Tcrit12 = (132.0*647.3)**0.5
	Mmix = (1.0/28.97+1.0/18.015)**0.5
	return 0.000364*(Tk/Tcrit12)**2.334*pcrit13*tcrit512*Mmix/(Pair/1013.25)*0.0001

def therm_cond(Tk):
	mair = 28.97
	rgas = 8314.34
	rair = rgas/mair
	cp = 1003.5
	return (cp+1.25*rair)*viscosity(Tk)

def h_cylinder_in_air(Tk,Pair,speed,minspeed,diamwick):
	mair = 28.97
	rgas = 8314.34
	rair = rgas/mair
	cp = 1003.5
	Pr = cp/(cp + (1.25*rair))
	thermcon = therm_cond(Tk)

	density = Pair*100.0/(rair*Tk)
	if speed<minspeed:
		speed=minspeed
	
	Re = speed*density*diamwick/viscosity(Tk)
	Nu = 0.281*Re**0.6*Pr**0.44
	return Nu*thermcon/diamwick

def h_evap(Tk):
	return (313.15-Tk)/30.0*(-71100.0) + 2407300.0

def viscosity(Tk):
	omega = (Tk/97.0-2.9)/0.4*(-0.034)+1.048
	return 0.0000026693*(28.97*Tk)**0.5/(3.617**2.0*omega)


def thermal_cond(Tk):
	mair = 28.97
	rgas = 8314.34
	rair = rgas/mair
	cp = 1003.5
	return (cp+1.25*rair)*viscosity(Tk)

def h_sphere_in_air(Tk, Pair, wind, minspeed, diamglobe):
	mair = 28.97
	rgas = 8314.34
	rair = rgas/mair
	cp = 1003.5
	Pr = cp/(cp+1.25*rair)
	density = Pair*100.0/(rair*Tk)
	if wind<minspeed:
		wind = minspeed
	Re = wind*density*diamglobe/viscosity(Tk)
	Nu = 2.0+0.6*Re**0.5*Pr**0.3333
  	# Convective heat tranfer coefficient for flow around a sphere, W/(m2 K)
	return Nu*thermal_cond(Tk)/diamglobe
	
def esat(Tk):
	#Tk is in Kelvin.
	esat = 6.1121*np.exp(17.502*(Tk-273.15)/(Tk-32.18))
	return 1.004*esat # Contains correction for moist air, if pressure is not available; or pressure > 800 mb

def emis_atm(Tk,RH):
	#Note RH must be in fraction
	es = RH*esat(Tk)
	return 0.575*es**0.143
	
def find_closest(A,target):
	#A must be sorted
	idx = A.searchsorted(target)
	idx = np.clip(idx,1,len(A)-1)
	left = A[idx-1]
	right = A[idx]
	idx -= target -left < right - target
	return idx

if __name__=="__main__":
	filelist = glob('./**/wrfout*00',recursive=True)
	#highdom=max([int(ifil.split('_')[-3].replace('d','')) for ifil in filelist])
	#filelist = glob(f'./**/wrfout*d{highdom:02d}*',recursive=True)
	rank = mpi4py.MPI.COMM_WORLD.Get_rank()
	size = mpi4py.MPI.COMM_WORLD.Get_size()
	if rank==0:
		print(f"Found {len(filelist)} wrfout files.")
	for i,ifil in enumerate(filelist):
		if i%size!=rank: continue
		heat_index_wrf(ifil)
