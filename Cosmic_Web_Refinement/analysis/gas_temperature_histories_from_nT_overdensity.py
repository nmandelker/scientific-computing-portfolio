from gadget import *
from gadget_subfind import *
import const as c
import matplotlib as mpl
import os, sys

######################################################
### Input parameters ###

exclude = 5	# 0: Include all cells;	1: Exclude allcels in FoF groups; 2: Exclude all cells in Subfind subhalos; 3: Exclude all cells within Rvir of central halos if FoF groups; 4: Exclude halos in projection
min_mass  = 1.0e7	# Minimum halo mas for exclusion rules
min_Npart = 32		# Minimum number of DM particles in FoF groups that are excluded if exclude==5
xmax = 100.0								# proper kpc --> Integration depth is +/- xmax from midplane of sheet
renorm = 0	# 1 to normalize each histogram to its fraction of the total mass, 0 to normalize each histogram to unity
Z_from_ZF1 = 1 # Divide volume of each simulation based on projected metallicity distribution in ZF1 to make sure to probe same volume, and be consistent with covering fractions
Dweight = 0	   # 1 for density-weighted average metallicity, 0 for mass-weighted average metallicity
Pthresh = 0		# Only include cells with pressure > Pthresh, to roughly eliminate pre-shock material outside of the sheet. Set to 0 for no threshold

min_mass /= 1.0e10	# halo masses are given in units of 10^{10} solar masses

ZFs   = np.array([1, 0, 2, 4])
NZF = 4
#ZFs   = np.array([1, 0, 2])
#NZF = 3

SnapNumSelect = 91	# Snapshot where cells are selected for history tracking, based on T and n thresholds below
Tmin_box = 4.6
Tmax_box = 5.1
nmin_box = -3.8
nmax_box = -2.6

#SnapNums  = np.arange(78,110)
SnapNums  = np.arange(82,110)
Nsnap = np.size(SnapNums)
redshift_list = np.zeros([Nsnap])

zL   = 1800	# comoving h^{-1} kpc, for sheet boundary (+/-zL)
yL   = 1000	# comoving h^{-1} kpc, for sheet boundary (+/-yL)
imgcent = [0.,0.,0.]
pixres = 2.0			# comiving h^-1 kpc

Zmaxs = np.array([-3.0])							# Bins of metalicity with these boundaries
NZ_bins = np.size(Zmaxs)

######################################################

UnitLength_in_cm         = 3.085678e21	# kpc/h
UnitMass_in_g            = 1.989e43	# 1.0e10 Msun/h
UnitVelocity_in_cm_per_s = 1.0e5	# 1 km/sec
UnitTime_in_s       = UnitLength_in_cm / UnitVelocity_in_cm_per_s 
UnitDensity_in_cgs  = UnitMass_in_g    / UnitLength_in_cm**3
UnitPressure_in_cgs = UnitMass_in_g    / UnitLength_in_cm    / UnitTime_in_s**2
UnitEnergy_in_cgs   = UnitMass_in_g    * UnitLength_in_cm**2 / UnitTime_in_s**2

gravity = c.G                        #cm^3/(g s^2)
msolar = c.msol                      #gram
mproton = c.amu                      #gram
kB = c.KB                            #erg/K
kpcincm = c.parsec*1e3
Xh = 0.76
gamma_gas = 5.0/3.0
Omegab = 0.0486
HUBBLE = 3.2407789e-18 # in h/sec
sec_in_Gyr = 60.0*60.0*24.0*365.25*1.0e9

UnitNumberDensity_in_cgs = UnitDensity_in_cgs/mproton	# ntot in cm^{-3}

######################################################

#xsize   = 5.3
#ysize   = 5.3
#fsize   = 17
#fsize2  = 19
#fweight = 150
#ax1_pos  = [0.16, 0.12, 0.80, 0.80]

xsize   = 19.1
ysize   = 5.3
fsize   = 18
fsize2  = 22
fweight = 150
ax1_pos  = [0.08, 0.12, 0.23, 0.80]
ax2_pos  = [0.37, 0.12, 0.23, 0.80]
ax3_pos  = [0.66, 0.12, 0.23, 0.80]

xstr = r"$z$"
ystr_mass = r"${\rm log}(M_{\rm gas}\:{\rm [M_{\odot}]})$"
ystr_temp = r"${\rm log}(T\:{\rm [K]})$"
ystr_dens = r"${\rm log}(n_{\rm H}\:{\rm [cm^{-3}]})$"

xlaby = -0.08
ylabx = -0.12
tit_pos = [0.5, 1.005]

xtick         = 1.0 + np.arange(0,15)*0.5
xminorLocator = MultipleLocator(0.1)
if NZF>3:
        xlim = [2.5, 6.3]
else:
        xlim = [2.0,6.3]

ytick_mass =  6.0 + np.arange(0,13)*0.5
yminorLocator_mass = MultipleLocator(0.1)
ylim_mass = [7.0,  11.0]

ytick_temp =  3.0 + np.arange(0,10)*0.5
yminorLocator_temp = MultipleLocator(0.1)
ylim_temp = [ 3.5,  7.0]

ytick_dens = -6.0 + np.arange(0,9)*1.0
yminorLocator_dens = MultipleLocator(0.2)
ylim_dens = [-5.5,  0.5]

######################################################

FigDir   = "../plots/FilB-Zoom4"
if not os.path.exists(FigDir):
        os.mkdir(FigDir)
FigDir   = "%s/for_letter" % (FigDir)
if not os.path.exists(FigDir):
        os.mkdir(FigDir)
FigDir   = "%s/Temperature_histories_from_nT_overdensity" % (FigDir)
if not os.path.exists(FigDir):
        os.mkdir(FigDir)
if Z_from_ZF1==1:
	FigDir   = "%s/Z_from_ZF1" % (FigDir)
	if not os.path.exists(FigDir):
		os.mkdir(FigDir)
	if Dweight==1:
		FigDir   = "%s/Dweight" % (FigDir)
		if not os.path.exists(FigDir):
			os.mkdir(FigDir)
FigDir   = FigDir+"/Tmin_"+'{0:.1f}'.format(Tmin_box)+"_Tmax_"+'{0:.1f}'.format(Tmax_box)+"_nmin_"+'{0:.1f}'.format(nmin_box)+"_nmax_"+'{0:.1f}'.format(nmax_box)
if not os.path.exists(FigDir):
	os.mkdir(FigDir)

######################################################

if renorm==0:
	if exclude<5:
		fig_tit_low_Z = "%s/selected_at_snap_%03d_xmax_%03d_exclude_%01d_Mmin_%01d_Pmin_%02d_low_Z.png" % (FigDir,SnapNumSelect,xmax,exclude,np.log10(min_mass)+10,Pthresh)
		fig_tit_high_Z = "%s/selected_at_snap_%03d_xmax_%03d_exclude_%01d_Mmin_%01d_Pmin_%02d_high_Z.png" % (FigDir,SnapNumSelect,xmax,exclude,np.log10(min_mass)+10,Pthresh)
	elif exclude==5:
		fig_tit_low_Z = "%s/selected_at_snap_%03d_xmax_%03d_exclude_%01d_Nmin_%02d_Pmin_%02d_low_Z.png" % (FigDir,SnapNumSelect,xmax,exclude,min_Npart,Pthresh)
		fig_tit_high_Z = "%s/selected_at_snap_%03d_xmax_%03d_exclude_%01d_Nmin_%02d_Pmin_%02d_high_Z.png" % (FigDir,SnapNumSelect,xmax,exclude,min_Npart,Pthresh)
elif renorm==1:
	if exclude<5:
		fig_tit_low_Z = "%s/selected_at_snap_%03d_xmax_%03d_exclude_%01d_Mmin_%01d_Pmin_%02d_low_Z_renormalized.png" % (FigDir,SnapNumSelect,xmax,exclude,min_Npart,Pthresh)
		fig_tit_high_Z = "%s/selected_at_snap_%03d_xmax_%03d_exclude_%01d_Mmin_%01d_Pmin_%02d_high_Z_renormalized.png" % (FigDir,SnapNumSelect,xmax,exclude,min_Npart,Pthresh)
	elif exclude==5:
		fig_tit_low_Z = "%s/selected_at_snap_%03d_xmax_%03d_exclude_%01d_Nmin_%02d_Pmin_%02d_low_Z_renormalized.png" % (FigDir,SnapNumSelect,xmax,exclude,min_Npart,Pthresh)
		fig_tit_high_Z = "%s/selected_at_snap_%03d_xmax_%03d_exclude_%01d_Nmin_%02d_Pmin_%02d_high_Z_renormalized.png" % (FigDir,SnapNumSelect,xmax,exclude,min_Npart,Pthresh)

ttl_high_Z = r"$"+'{0:.1f}'.format(Zmaxs[0])+r"{\rm <[Z]}$"
ttl_low_Z  = r"${\rm [Z]<}"+'{0:.1f}'.format(Zmaxs[0])+"$"

######################################################

print("Selecting cells from snapshot ",SnapNumSelect)
print("with Log(T) between ",Tmin_box," and",Tmax_box)
print("and Log(n_H) between ",nmin_box," and",nmax_box)

SnapNum = SnapNumSelect
Ncell_metal_poor = np.zeros(NZF)
Ncell_metal_rich = np.zeros(NZF)
for ZF_ind in np.arange(0,NZF):
	ZF = ZFs[ZF_ind]
	if ZF_ind==0:	# I rearanged the ZFs array to have ZF1 first, in order to use its metallicity distribution for all ZFs. This helps store the output in the correct order
		ZF_ind2 = 1
	elif ZF_ind==1:
		ZF_ind2 = 0
	else:
		ZF_ind2 = ZF_ind 
	print(" ")
	print('ZF, ZF_ind, ZF_ind2= ',ZF, ZF_ind, ZF_ind2)

	SimLabel = 'FilB-Zoom%01d' % (ZF)
	in_filename = "halo-positions-%s.dat" % (SimLabel)
	SimDir   = "/home/mandelnr/basement/%s/output" % (SimLabel)
	SnapBase = 'snap_'

	posdata  = np.loadtxt(in_filename)
	RevSnap2 = posdata[0,:]
	aexp_vec2 = posdata[1,:]
	x1_vec    = posdata[2,:]
	y1_vec    = posdata[3,:]
	z1_vec    = posdata[4,:]
	x2_vec    = posdata[5,:]
	y2_vec    = posdata[6,:]
	z2_vec    = posdata[7,:]
	ind1_vec  = posdata[8,:]
	ind2_vec  = posdata[9,:]

	if ZF==4:
		SnapNum2 = SnapNum - 78
	else:
		SnapNum2 = SnapNum

	s  = gadget_readsnap( SnapNum2, SimDir, SnapBase, loadonlytype=[0], loadonly=['pos', 'rho', 'mass', 'gmet', 'u', 'ne', 'gz', 'id'], lazy_load=True )
	sf = load_subfind( SnapNum2, dir=SimDir+'/' )
	redshift = s.redshift
	BoxHalf = 0.5 * s.boxsize
	aexp = s.time
	h0 = s.hubbleparam
	BoxHalf *= aexp/h0
	BoxSize = 2.0*BoxHalf
	boxx = 2.0*zL			# comiving h^-1 kpc
	boxy = 2.0*yL			# comiving h^-1 kpc
	boxz = 2.0*xmax*(h0/aexp)	# comiving h^-1 kpc

	### Gas property arrays ###
	gas_id = s.data['id']
	Metal  = s.data['gz'] / 0.0127	# Solar units
	Mass   = 1.0*s.mass
	ugas   = s.u * UnitEnergy_in_cgs / UnitMass_in_g			# internal energy per unit mass in erg/gr
	MeanWeight = 4.0 / ( 3.0*Xh + 1.0 + 4.0*Xh*s.ne ) * mproton
	Temp   = ugas * MeanWeight/kB * (gamma_gas-1.0)				# temperature in Kelvin
	Pres   = Temp * s.rho * UnitDensity_in_cgs / MeanWeight		# P/kB [K cm^{-3}]
	Dens   = s.rho * UnitDensity_in_cgs
	nHgas  = Dens * Xh / mproton
#	nHgas  = s.rho * s.gmet[:,0] * UnitNumberDensity_in_cgs	# nH  cm^{-3}
#	s.data['rho']   = s.rho * s.gmet[:,0] * UnitNumberDensity_in_cgs	# nH  cm^{-3}

	if Z_from_ZF1==1:
		s.data['gz']    /= 0.0127	# Solar units
		if Dweight==1:
			s.data['rho']   *= s.gmet[:,0] * UnitNumberDensity_in_cgs	# nH  cm^{-3}
			s.data['DZ']     = s.rho * s.gz
		else:
			s.data['MZ']     = s.mass * s.gz

	b2, = np.where(RevSnap2-SnapNum==0)
	if ZF<4:
		if b2.size>0:
			if(np.abs(aexp_vec2[b2]-aexp)>0.0001):
				print(SnapNum)
				print(aexp)
				print(aexp_vec2[b2])
				sys.exit("Error: Snapshot aexp Mismatch")
		else:
			print(SnapNum)
			sys.exit("Error: Current snapshot is not in file")

	print(SnapNum,RevSnap2[b2])
	x_first = float(x1_vec[b2])
	y_first = float(y1_vec[b2])
	z_first = float(z1_vec[b2])
	x_second = float(x2_vec[b2])
	y_second = float(y2_vec[b2])
	z_second = float(z2_vec[b2])
	ind_first  = int(ind1_vec[b2])
	ind_second = int(ind2_vec[b2])
	print(ind_first, ind_second)
	print(x_first, y_first, z_first)
	print(x_second, y_second, z_second)
	print(sf.data['fpos'][ind_first,0], sf.data['fpos'][ind_first,1], sf.data['fpos'][ind_first,2])
	print(sf.data['fpos'][ind_second,0], sf.data['fpos'][ind_second,1], sf.data['fpos'][ind_second,2])
	print(sf.data['fpos'][0,0], sf.data['fpos'][0,1], sf.data['fpos'][0,2])
	print(sf.data['fpos'][1,0], sf.data['fpos'][1,1], sf.data['fpos'][1,2])
	print(sf.data['fpos'][2,0], sf.data['fpos'][2,1], sf.data['fpos'][2,2])

	if ZF==4:
		if SnapNum==99:
			ind_first = 1

			x_first = sf.data['fpos'][ind_first,0]
			y_first = sf.data['fpos'][ind_first,1]
			z_first = sf.data['fpos'][ind_first,2]

			x_second = sf.data['fpos'][ind_second,0]
			y_second = sf.data['fpos'][ind_second,1]
			z_second = sf.data['fpos'][ind_second,2]

	# determine mid-point between two largest halos (Halos with most particles)
	xc = x_first - x_second
	yc = y_first - y_second
	zc = z_first - z_second
	
	# Accounts for periodic box
	if xc<-BoxHalf:
		xc += BoxSize
	elif xc>BoxHalf:
		xc -= BoxSize
	if yc<-BoxHalf:
		yc += BoxSize
	elif yc>BoxHalf:
		yc -= BoxSize
	if zc<-BoxHalf:
		zc += BoxSize
	elif zc>BoxHalf:
		zc -= BoxSize

	xc *= 0.5
	yc *= 0.5
	zc *= 0.5

	# Position of midpoint between halos
	xc += x_second
	yc += y_second
	zc += z_second

	# compute coordinates relative to the center of the filament   (coordinates of gas cells)
	xp = s.data['pos'][:,0] - xc
	yp = s.data['pos'][:,1] - yc
	zp = s.data['pos'][:,2] - zc

	# also for groups
	gxp = sf.data['fpos'][:,0] - xc
	gyp = sf.data['fpos'][:,1] - yc
	gzp = sf.data['fpos'][:,2] - zc
	gR = sf.data['frt2']
	#gM = sf.data['fmas']
	gM = sf.data['fmt2']
	Ng = sf.totngroups
	
	# in case we have a periodic wrap, let's take care of it (should actually not be necessary here)
	b1, = np.where(xp<-BoxHalf)
	if b1.size>0:
		xp[b1] += BoxSize
	b1, = np.where(xp>BoxHalf)
	if b1.size>0:
		xp[b1] -= BoxSize
	b1, = np.where(yp<-BoxHalf)
	if b1.size>0:
		yp[b1] += BoxSize
	b1, = np.where(yp>BoxHalf)
	if b1.size>0:
		yp[b1] -= BoxSize
	b1, = np.where(zp<-BoxHalf)
	if b1.size>0:
		zp[b1] += BoxSize
	b1, = np.where(zp>BoxHalf)
	if b1.size>0:
		zp[b1] -= BoxSize

	b1, = np.where(gxp<-BoxHalf)
	if b1.size>0:
		gxp[b1] += BoxSize
	b1, = np.where(gxp>BoxHalf)
	if b1.size>0:
		gxp[b1] -= BoxSize
	b1, = np.where(gyp<-BoxHalf)
	if b1.size>0:
		gyp[b1] += BoxSize
	b1, = np.where(gyp>BoxHalf)
	if b1.size>0:
		gyp[b1] -= BoxSize
	b1, = np.where(gzp<-BoxHalf)
	if b1.size>0:
		gzp[b1] += BoxSize
	b1, = np.where(gzp>BoxHalf)
	if b1.size>0:
		gzp[b1] -= BoxSize

	# now turn the coordinates into a frame were the filament is aligned with the z-axis
	# first need to find a suitable coordinate transformation
	Jx = x_first - x_second
	Jy = y_first - y_second
	Jz = z_first - z_second
	if Jx<-BoxHalf:
		Jx += BoxSize
	elif Jx>BoxHalf:
		Jx -= BoxSize
	if Jy<-BoxHalf:
		Jy += BoxSize
	elif Jy>BoxHalf:
		Jy -= BoxSize
	if Jz<-BoxHalf:
		Jz += BoxSize
	elif Jz>BoxHalf:
		Jz -= BoxSize

	JJ = np.sqrt(Jx**2 + Jy**2 + Jz**2)
	rr = JJ
	Jx /= JJ
	Jy /= JJ
	Jz /= JJ
	if Jx<0:
		Jx *= -1
		Jy *= -1
		Jz *= -1
	print('Jx,y,z',Jx,Jy,Jz)

	# determine x-direction (nx, ny, nz)
	nx = -Jy
	ny = Jx
	nz = 0
	if nx<0:
		nx*= -1
		ny*= -1
	n = np.sqrt(nx**2 + ny**2 + nz**2)
	nx /= n
	ny /= n
	nz /= n
	print('nx,y,z',nx,ny,nz)

	# determine y-direction (mx, my, mz)
	mx = Jy * nz - Jz * ny
	my = Jz * nx - Jx * nz
	mz = Jx * ny - Jy * nx
	print('mx,y,z',mx,my,mz)

	# Rotate the XY plane by 45 degrees to try and capture empirical pancake
	nxx = nx*np.cos(np.pi/4.) + mx*np.sin(np.pi/4.)
	nyy = ny*np.cos(np.pi/4.) + my*np.sin(np.pi/4.)
	nzz = nz*np.cos(np.pi/4.) + mz*np.sin(np.pi/4.)

	mxx = mx*np.cos(np.pi/4.) - nx*np.sin(np.pi/4.)
	myy = my*np.cos(np.pi/4.) - ny*np.sin(np.pi/4.)
	mzz = mz*np.cos(np.pi/4.) - nz*np.sin(np.pi/4.)

	nx = nxx
	ny = nyy
	nz = nzz
	mx = mxx
	my = myy
	mz = mzz

	# now turn the coordinates into a frame were the filament is aligned with the z-axis
	xx = xp*nx + yp*ny + zp*nz	# x coordinate of point in disk frame
	yy = xp*mx + yp*my + zp*mz	# y coordinate of point in disk frame
	zz = xp*Jx + yp*Jy + zp*Jz	# z coordinate of point in disk frame

	gxx = gxp*nx + gyp*ny + gzp*nz	# x coordinate of point in disk frame
	gyy = gxp*mx + gyp*my + gzp*mz	# y coordinate of point in disk frame
	gzz = gxp*Jx + gyp*Jy + gzp*Jz	# z coordinate of point in disk frame

	xx *= h0/aexp			 # comiving h^-1 kpc
	yy *= h0/aexp			 # comiving h^-1 kpc
	zz *= h0/aexp			 # comiving h^-1 kpc

	gxx *= h0/aexp			 # comiving h^-1 kpc
	gyy *= h0/aexp			 # comiving h^-1 kpc
	gzz *= h0/aexp			 # comiving h^-1 kpc
	gR  *= h0/aexp			 # comiving h^-1 kpc

	s.data['pos'][:,0] = xx			# comiving h^-1 kpc
	s.data['pos'][:,1] = yy			# comiving h^-1 kpc
	s.data['pos'][:,2] = zz			# comiving h^-1 kpc

	################################
	### Select sheet region ###

	b1, = np.where( (np.abs(zz)<=0.5*boxx) & (np.abs(yy)<=0.5*boxy) & (np.abs(xx)<=0.5*boxz) )
	print( 'Ncell with halos=',np.size(b1) )

	if exclude==1:		### Remove gas cells associated with FoF groups ###
		print( 'exclude 1' )
		Ng = sf.totngroups
		first_cell = 0
		Ntrack_g = (Ng - np.mod(Ng,10)) / 10
		for g_ind in np.arange(0,Ng):
			if( (np.mod(g_ind,Ntrack_g)==0) | (g_ind==0) ):
				print("g_ind, Ng= ", g_ind, Ng)

			Ncell = sf.data['flty'][g_ind,0]
			last_cell = first_cell+Ncell
			if sf.data['fmas'][g_ind]>=min_mass:
				Mass[first_cell:last_cell] = 0.0
			first_cell += Ncell
	elif exclude==5:		### Remove gas cells associated with FoF groups above minimal number of DM particles ###
		print( 'exclude 5' )
		Ng = sf.totngroups
		first_cell = 0
		Ntrack_g = (Ng - np.mod(Ng,10)) / 10
		for g_ind in np.arange(0,Ng):
			if( (np.mod(g_ind,Ntrack_g)==0) | (g_ind==0) ):
				print("g_ind, Ng= ", g_ind, Ng)

			Ncell = sf.data['flty'][g_ind,0]
			last_cell = first_cell+Ncell
			if sf.data['flty'][g_ind,1]>=min_Npart:
				Mass[first_cell:last_cell] = 0.0
			first_cell += Ncell
	elif exclude==2:	### Remove gas cells associated with subfind subhalos with M>min_mass ###
		print( 'exclude 2' )
		first_cell = 0
		for g_ind in np.arange(0,Ng):
			first_sub = sf.data['ffsh'][g_ind]
			if first_sub>-1:
				Nsub = sf.data['fnsh'][g_ind]
				for s_ind in np.arange(0,Nsub):
					Ncell = sf.data['slty'][first_sub+s_ind,0]
					last_cell = first_cell+Ncell
					if sf.data['smas'][first_sub+s_ind]>min_mass:
						Mass[first_cell:last_cell] = 0.0
					first_cell += Ncell
			first_cell = np.int(np.sum(sf.data['flty'][0:g_ind,0]))

	Temp  = np.log10(Temp)
	nHgas = np.log10(nHgas)
	b1, = np.where( (np.abs(zz)<=0.5*boxx) & (np.abs(yy)<=0.5*boxy) & (np.abs(xx)<=0.5*boxz) & (Mass>0.0) & (Temp>Tmin_box) & (Temp<Tmax_box) & (nHgas>nmin_box) & (nHgas<nmax_box) )
	ID1    = gas_id[b1]			# solar units
	mass1  = Mass[b1]			# 10^{10} solar mass
	metal1 = Metal[b1]			# solar units
	xx1    = xx[b1]				# comiving h^-1 kpc
	yy1    = yy[b1]				# comiving h^-1 kpc
	zz1    = zz[b1]				# comiving h^-1 kpc

	Ncell = np.size(mass1)
	print( 'Intermediate Ncell=', np.int(np.size(mass1)) )
	if exclude==3:
		print( 'exclude 3' )
		bg2, = np.where( ((gzz+gR)>=-0.5*boxx) & ((gzz-gR)<=0.5*boxx) & ((gyy+gR)>=-0.5*boxy) & ((gyy-gR)<=0.5*boxy) & ((gxx+gR)>=-0.5*boxz) & ((gxx-gR)<=0.5*boxz) & (gM>=min_mass) & (gR>0.0) )
		Ng2 = np.size(bg2)
		gx2 = gxx[bg2]
		gy2 = gyy[bg2]
		gz2 = gzz[bg2]
		gR2 = gR[bg2]
		print("Ng2=",Ng2)
		if Ng2>0:
			Ntrack_g = (Ng2 - np.mod(Ng2,10)) / 10
			for g_ind in np.arange(0,Ng2):
				if( (np.mod(g_ind,Ntrack_g)==0) | (g_ind==0) ):
					print("g_ind, Ng2= ", g_ind, Ng2)
				b3, = np.where( np.sqrt( (xx1-gx2[g_ind])**2 + (yy1-gy2[g_ind])**2 + (zz1-gz2[g_ind])**2 ) <= gR2[g_ind] )
				mass1[b3]  = 0.0
	elif exclude==4:
		print( 'exclude 4' )
		bg2, = np.where( ((gzz+gR)>=-0.5*boxx) & ((gzz-gR)<=0.5*boxx) & ((gyy+gR)>=-0.5*boxy) & ((gyy-gR)<=0.5*boxy) & ((gxx+gR)>=-0.5*boxz) & ((gxx-gR)<=0.5*boxz) & (gM>=min_mass) & (gR>0.0) )
		Ng2 = np.size(bg2)
		gx2 = gxx[bg2]
		gy2 = gyy[bg2]
		gz2 = gzz[bg2]
		gR2 = gR[bg2]
		print("Ng2=",Ng2)
		if Ng2>0:
			Ntrack_g = (Ng2 - np.mod(Ng2,10)) / 10
			for g_ind in np.arange(0,Ng2):
				if( (np.mod(g_ind,Ntrack_g)==0) | (g_ind==0) ):
					print("g_ind, Ng2= ", g_ind, Ng2)
				b3, = np.where( np.sqrt( (yy1-gy2[g_ind])**2 + (zz1-gz2[g_ind])**2 ) <= gR2[g_ind] )
				mass1[b3]  = 0.0

	b1, = np.where( mass1>0.0 )
	ID2    = ID1[b1]			# 10^{10} solar mass
	mass2  = mass1[b1]			# 10^{10} solar mass
	metal2 = metal1[b1]			# solar units
	xx2    = xx1[b1]			# comiving h^-1 kpc
	yy2    = yy1[b1]			# comiving h^-1 kpc
	zz2    = zz1[b1]			# comiving h^-1 kpc

	print( 'Final Ncell=', np.int(np.size(mass2)) ) 
	print( 'Log10 total mass=',  np.log10(np.sum(mass2) * 1.0e10) )
		
	if Z_from_ZF1==1:
		print( 'Z_from_ZF1' )
		Ncell = np.size(mass2)
		Ntrack = (Ncell - np.mod(Ncell,10)) / 10
		if ZF==1:
			boxx = 2.0*zL			# comiving h^-1 kpc
			boxy = 2.0*yL			# comiving h^-1 kpc
			boxz = 2.0*xmax*(h0/aexp)	# comiving h^-1 kpc
			Ncell_x = int(boxx/pixres)
			Ncell_y = int(boxy/pixres)
			if Dweight==0:
				Z_proj_Q = s.get_Aslice("MZ", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz,  axes = [2,1], proj = True, numthreads=32)
				Z_proj_W = s.get_Aslice("mass", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz,  axes = [2,1], proj = True, numthreads=32)
			elif Dweight==1:
				Z_proj_Q = s.get_Aslice("DZ", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz,  axes = [2,1], proj = True, numthreads=32)
				Z_proj_W = s.get_Aslice("rho", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz,  axes = [2,1], proj = True, numthreads=32)

			map_Z = np.transpose( Z_proj_Q['grid'] / Z_proj_W['grid'] )
			xpl   = Z_proj_Q['x']
			ypl   = Z_proj_Q['y']
		print("Ncell_x, Ncell_y, shape(map_Z) ", Ncell_x, Ncell_y, np.shape(map_Z))

		### Version 2 ###
		Ntrack = (Ncell_x - np.mod(Ncell_x,10)) / 10
		for x_ind in np.arange(0,Ncell_x):
			if( (np.mod(x_ind,Ntrack)==0) | (x_ind==0) ):
				print("x_ind of Ncell_x ",x_ind,Ncell_x)
			b3, = np.where( (zz2>=xpl[x_ind]) & (zz2<xpl[x_ind+1]) )
			yy3 = yy2[b3]
			for y_ind in np.arange(0,Ncell_y):
				b4, = np.where( (yy3>=ypl[y_ind]) & (yy3<ypl[y_ind+1]) )
				metal2[b3[b4]] = map_Z[y_ind,x_ind]

		print("min, max metal2 = ", np.min(metal2), np.max(metal2))

	if renorm==0:
		if exclude<5:
			cell_list_outname_low_Z = "%s/cell_IDs_at_snap_%03d_ZF%01d_xmax_%03d_exclude_%01d_Mmin_%01d_Pmin_%02d_low_Z.txt" % (FigDir,SnapNumSelect,ZF,xmax,exclude,np.log10(min_mass)+10,Pthresh)
			cell_list_outname_high_Z = "%s/cell_IDs_at_snap_%03d_ZF%01d_xmax_%03d_exclude_%01d_Mmin_%01d_Pmin_%02d_high_Z.txt" % (FigDir,SnapNumSelect,ZF,xmax,exclude,np.log10(min_mass)+10,Pthresh)
		elif exclude==5:
			cell_list_outname_low_Z = "%s/cell_IDs_at_snap_%03d_ZF%01d_xmax_%03d_exclude_%01d_Nmin_%02d_Pmin_%02d_low_Z.txt" % (FigDir,SnapNumSelect,ZF,xmax,exclude,min_Npart,Pthresh)
			cell_list_outname_high_Z = "%s/cell_IDs_at_snap_%03d_ZF%01d_xmax_%03d_exclude_%01d_Nmin_%02d_Pmin_%02d_high_Z.txt" % (FigDir,SnapNumSelect,ZF,xmax,exclude,min_Npart,Pthresh)
	elif renorm==1:
		if exclude<5:
			cell_list_outname_low_Z = "%s/cell_IDs_at_snap_%03d_ZF%01d_xmax_%03d_exclude_%01d_Mmin_%01d_Pmin_%02d_low_Z_renormalized.txt" % (FigDir,SnapNumSelect,ZF,xmax,exclude,np.log10(min_mass)+10,Pthresh)
			cell_list_outname_high_Z = "%s/cell_IDs_at_snap_%03d_ZF%01d_xmax_%03d_exclude_%01d_Mmin_%01d_Pmin_%02d_high_Z_renormalized.txt" % (FigDir,SnapNumSelect,ZF,xmax,exclude,np.log10(min_mass)+10,Pthresh)
		elif exclude==5:
			cell_list_outname_low_Z = "%s/cell_IDs_at_snap_%03d_ZF%01d_xmax_%03d_exclude_%01d_Nmin_%02d_Pmin_%02d_low_Z_renormalized.txt" % (FigDir,SnapNumSelect,ZF,xmax,exclude,min_Npart,Pthresh)
			cell_list_outname_high_Z = "%s/cell_IDs_at_snap_%03d_ZF%01d_xmax_%03d_exclude_%01d_Nmin_%02d_Pmin_%02d_high_Z_renormalized.txt" % (FigDir,SnapNumSelect,ZF,xmax,exclude,min_Npart,Pthresh)

	metal2 = np.log10(metal2)
	b1, = np.where( metal2 < Zmaxs[0] )
	Ncell_metal_poor[ZF_ind2] = np.size( b1 )
	if ZF==0:
		id_list_0 = ID2[b1]
		np.savetxt( cell_list_outname_low_Z, np.sort(id_list_0) )
	elif ZF==1:
		id_list_1 = ID2[b1]
		np.savetxt( cell_list_outname_low_Z, np.sort(id_list_1) )
	elif ZF==2:
		id_list_2 = ID2[b1]
		np.savetxt( cell_list_outname_low_Z, np.sort(id_list_2) )
	elif ZF==4:
		id_list_4 = ID2[b1]
		np.savetxt( cell_list_outname_low_Z, np.sort(id_list_4) )
	print("Wrote selected cell IDs to ", cell_list_outname_low_Z)

######################################################

print(" ")
print("Done selecting cells")
print("Writing list of cell IDs")


print("Tracking cell histories")
Temp_history = np.zeros([Nsnap, NZF, 2])
dens_history = np.zeros([Nsnap, NZF, 2])
mass_history = np.zeros([Nsnap, NZF, 4])


for ZF_ind in np.arange(0,NZF):
	ZF = ZFs[ZF_ind]
	if ZF_ind==0:	# I rearanged the ZFs array to have ZF1 first, in order to use its metallicity distribution for all ZFs. This helps store the output in the correct order
		ZF_ind2 = 1
	elif ZF_ind==1:
		ZF_ind2 = 0
	else:
		ZF_ind2 = ZF_ind 
	print('ZF, ZF_ind, ZF_ind2= ',ZF, ZF_ind, ZF_ind2)

	if ZF==0:
		id_list = np.sort(id_list_0)
	elif ZF==1:
		id_list = np.sort(id_list_1)
	elif ZF==2:
		id_list = np.sort(id_list_2)
	elif ZF==4:
		id_list = np.sort(id_list_4)

	SimLabel = 'FilB-Zoom%01d' % (ZF)
	SimDir   = "/home/mandelnr/basement/%s/output" % (SimLabel)
	SnapBase = 'snap_'
	
	OutDir   = "%s/gas_IDs_ZF%01d" % (FigDir,ZF)
	if not os.path.exists(OutDir):
		os.mkdir(OutDir)
	print("cell IDs per snapshot written in ", OutDir)

	for snap_ind in np.arange(0,Nsnap):
		SnapNum = SnapNums[snap_ind]
		print('ZF, SnapNum ', ZF, SnapNum)

		if ZF==4:
			SnapNum2 = SnapNum - 78
		else:
			SnapNum2 = SnapNum

		s  = gadget_readsnap( SnapNum2, SimDir, SnapBase, loadonlytype=[0], loadonly=['rho', 'mass', 'u', 'ne', 'id'], lazy_load=True )
		if ZF_ind==0:
			redshift_list[snap_ind] = s.redshift

		### Gas property arrays ###
		gas_id = s.data['id']
		Mass   = 1.0*s.mass
		ugas   = s.u * UnitEnergy_in_cgs / UnitMass_in_g			# internal energy per unit mass in erg/gr
		MeanWeight = 4.0 / ( 3.0*Xh + 1.0 + 4.0*Xh*s.ne ) * mproton
		Temp   = ugas * MeanWeight/kB * (gamma_gas-1.0)				# temperature in Kelvin
		nHgas  = s.rho * UnitDensity_in_cgs * Xh / mproton

		Temp  = np.log10(Temp)
		nHgas = np.log10(nHgas)

		Nlist = np.int(Ncell_metal_poor[ZF_ind2])
		cell_list_local = -1 * np.ones([Nlist], dtype=int)
		print("Nlist= ", Nlist)
		print(id_list)
		id_max = id_list[Nlist-1]	# This is a sorted list
		id_min = id_list[0]			# This is a sorted list
		print('max, min id_list= ',id_max, id_min)
		print('max, min gas_id= ',np.max(gas_id), np.min(gas_id))

		gas_id_sorted = np.array(np.sort(gas_id))
		gas_id_sorted_index_list = np.argsort(gas_id)
		Ngas_id = np.size(gas_id)
		b1, = np.where(gas_id_sorted>id_max)
		if( np.size(b1)>0 ):
			end = np.int(b1[0])
		else:
			end = np.int(Ngas_id)
		print("Ngas_id, end ", Ngas_id, end)

		#cell_list_outname = "%s/snap_%03d.txt" % (OutDir,SnapNum)
		#output_array = np.column_stack( [gas_id_sorted, gas_id_sorted_index_list, Mass[gas_id_sorted_index_list], Temp[gas_id_sorted_index_list], nHgas[gas_id_sorted_index_list]] )
		#np.savetxt(cell_list_outname, output_array)
		#print("wrote cell IDs to ", cell_list_outname)

		Ntrack = (Nlist - np.mod(Nlist,10)) / 10
		counter = np.int(0)
		last_match = np.int(0)
		local_end = np.int(end)
		#cell_ind = np.int(0)
		#while( (cell_ind<Nlist) & (counter<end) ):
		for cell_ind in np.arange(0,Nlist):
			if( (np.mod(cell_ind,Ntrack)==0) ):
				print("cell_ind of Nlist ", cell_ind, Nlist)
				print("counter of end ", counter, end)

			if(counter < end):
				local_end = np.int( np.min( [end, counter + id_list[cell_ind] - last_match] ) )
				b2, = np.where( gas_id_sorted[counter:local_end] == id_list[cell_ind] )
				if np.size(b2)==1:
					last_match = id_list[cell_ind]
					cell_list_local[cell_ind] = gas_id_sorted_index_list[counter+b2]
					counter += np.int(b2)+1
				elif np.size(b2)>1:
					print("b2= ", b2, " cell_ind= ", cell_ind, " id_list[cell_ind]= ", id_list[cell_ind], "gas_id[b2]= ", gas_id_sorted[b2])
					sys.exit("Error: more than 1 cell in snapshot matches a given listed cell ID")

		b3, = np.where(cell_list_local >= 0)
		print("size(b3)= ",np.size(b3))
		#print("b3")
		#print(b3)
		#print("cell_list_local[b3]")
		#print(cell_list_local[b3])
		if np.size(b3)>0:
			mass_history[snap_ind, ZF_ind2, 0] = np.sum( Mass[cell_list_local[b3]] )

			b4, = np.where( Temp[cell_list_local[b3]] < Tmin_box )
			mass_history[snap_ind, ZF_ind2, 1] = np.sum( Mass[cell_list_local[b3[b4]]] )

			b4, = np.where( (Temp[cell_list_local[b3]] >= Tmin_box) & (Temp[cell_list_local[b3]] < Tmax_box) )
			mass_history[snap_ind, ZF_ind2, 2] = np.sum( Mass[cell_list_local[b3[b4]]] )

			b4, = np.where( Temp[cell_list_local[b3]] >= Tmax_box )
			mass_history[snap_ind, ZF_ind2, 3] = np.sum( Mass[cell_list_local[b3[b4]]] )

			Temp_history[snap_ind, ZF_ind2, 0] = np.sum( Mass[cell_list_local[b3]] * Temp[cell_list_local[b3]] )     / mass_history[snap_ind, ZF_ind2, 0]
			dens_history[snap_ind, ZF_ind2, 0] = np.sum( Mass[cell_list_local[b3]] * nHgas[cell_list_local[b3]] )    / mass_history[snap_ind, ZF_ind2, 0]
			if np.size(b3)>1:
				Temp_history[snap_ind, ZF_ind2, 1] = np.sum( Mass[cell_list_local[b3]] * Temp[cell_list_local[b3]]**2 )  / mass_history[snap_ind, ZF_ind2, 0]
				dens_history[snap_ind, ZF_ind2, 1] = np.sum( Mass[cell_list_local[b3]] * nHgas[cell_list_local[b3]]**2 ) / mass_history[snap_ind, ZF_ind2, 0]

				Temp_history[snap_ind, ZF_ind2, 1] = np.sqrt( Temp_history[snap_ind, ZF_ind2, 1] - Temp_history[snap_ind, ZF_ind2, 0]**2 )
				dens_history[snap_ind, ZF_ind2, 1] = np.sqrt( dens_history[snap_ind, ZF_ind2, 1] - dens_history[snap_ind, ZF_ind2, 0]**2 )
	
	average_outname = "%s/selected_at_snap_%03d_average_evolution.txt" % (OutDir,SnapNumSelect)
	print("writing redshift evolution data to ", average_outname)
	output_array = np.column_stack( [redshift_list, mass_history[:, ZF_ind2, 0], mass_history[:, ZF_ind2, 1], mass_history[:, ZF_ind2, 2], mass_history[:, ZF_ind2, 3], Temp_history[:, ZF_ind2, 0], Temp_history[:, ZF_ind2, 1], dens_history[:, ZF_ind2, 0], dens_history[:, ZF_ind2, 1]] )
	np.savetxt(average_outname, output_array)

########################################################################################################################################################################################################################
########################################################################################################################################################################################################################

### Open figure object ###
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig = plt.figure(figsize = (xsize,ysize))

######################################################

ax1 = plt.axes(ax1_pos)

line1, = plt.plot( redshift_list, np.log10(mass_history[:,0,0])+10.0, linestyle='-', color='g', linewidth=2 )
line2, = plt.plot( redshift_list, np.log10(mass_history[:,1,0])+10.0, linestyle='-', color='b', linewidth=2 )
line3, = plt.plot( redshift_list, np.log10(mass_history[:,2,0])+10.0, linestyle='-', color='r', linewidth=2 )
if NZF>3:
        line4, = plt.plot( redshift_list, np.log10(mass_history[:,3,0])+10.0, linestyle='-', color='k', linewidth=2 )

plt.plot( redshift_list, np.log10(mass_history[:,0,1])+10.0, linestyle='--', color='g', linewidth=2 )
plt.plot( redshift_list, np.log10(mass_history[:,1,1])+10.0, linestyle='--', color='b', linewidth=2 )
plt.plot( redshift_list, np.log10(mass_history[:,2,1])+10.0, linestyle='--', color='r', linewidth=2 )
if NZF>3:
        plt.plot( redshift_list, np.log10(mass_history[:,3,1])+10.0, linestyle='--', color='k', linewidth=2 )

xlabel(xstr, fontname='times', fontsize = fsize2, fontweight = fweight)
ylabel(ystr_mass, fontname='times', fontsize = fsize2, fontweight = fweight)
ax1.xaxis.set_label_coords(0.5, xlaby)
ax1.yaxis.set_label_coords(ylabx, 0.5)

plt.xticks(xtick, fontname='arial', fontsize = fsize, fontweight = fweight)
ax1.xaxis.set_tick_params(which = 'major', direction = 'in', width = 1, length = 6, pad = 5, bottom = True, top = True)
ax1.xaxis.set_minor_locator(xminorLocator)
ax1.xaxis.set_tick_params(which = 'minor', direction = 'in', width = 1, length = 3, pad = 5, bottom = True, top = True)

plt.yticks(ytick_mass, fontname='arial', fontsize = fsize, fontweight = fweight)
ax1.yaxis.set_tick_params(which = 'major', direction = 'in', width = 1, length = 6, pad = 5, left = True, right = True)
ax1.yaxis.set_minor_locator(yminorLocator_mass)
ax1.yaxis.set_tick_params(which = 'minor', direction = 'in', width = 1, length = 3, pad = 5, left = True, right = True)

disp1 = r"${\rm ZF0.5}$"
disp2 = r"${\rm ZF1.0}$"
disp3 = r"${\rm ZF2.0}$"
disp4 = r"${\rm ZF4.0}$"
if NZF>3:
        plt.legend([line1, line2, line3, line4],[disp1, disp2, disp3, disp4], fontsize = 14, frameon=False, loc='upper left')
else:
        plt.legend([line1, line2, line3],[disp1, disp2, disp3], fontsize = 14, frameon=False, loc='upper left')

ax1.set_xlim(xlim)
ax1.set_ylim(ylim_mass)

######################################################

ax2 = plt.axes(ax2_pos)

plt.plot( redshift_list, Temp_history[:,0,0], linestyle='-', color='g', linewidth=2 )
plt.plot( redshift_list, Temp_history[:,1,0], linestyle='-', color='b', linewidth=2 )
plt.plot( redshift_list, Temp_history[:,2,0], linestyle='-', color='r', linewidth=2 )
if NZF>3:
        plt.plot( redshift_list, Temp_history[:,3,0], linestyle='-', color='k', linewidth=2 )

plt.plot( redshift_list, Temp_history[:,0,0]+Temp_history[:,0,1], linestyle='--', color='g', linewidth=1 )
plt.plot( redshift_list, Temp_history[:,1,0]+Temp_history[:,1,1], linestyle='--', color='b', linewidth=1 )
plt.plot( redshift_list, Temp_history[:,2,0]+Temp_history[:,2,1], linestyle='--', color='r', linewidth=1 )
if NZF>3:
        plt.plot( redshift_list, Temp_history[:,3,0]+Temp_history[:,3,1], linestyle='--', color='k', linewidth=1 )

plt.plot( redshift_list, Temp_history[:,0,0]-Temp_history[:,0,1], linestyle='--', color='g', linewidth=1 )
plt.plot( redshift_list, Temp_history[:,1,0]-Temp_history[:,1,1], linestyle='--', color='b', linewidth=1 )
plt.plot( redshift_list, Temp_history[:,2,0]-Temp_history[:,2,1], linestyle='--', color='r', linewidth=1 )
if NZF>3:
        plt.plot( redshift_list, Temp_history[:,3,0]-Temp_history[:,3,1], linestyle='--', color='k', linewidth=1 )

xlabel(xstr, fontname='times', fontsize = fsize2, fontweight = fweight)
ylabel(ystr_temp, fontname='times', fontsize = fsize2, fontweight = fweight)
ax2.xaxis.set_label_coords(0.5, xlaby)
ax2.yaxis.set_label_coords(ylabx, 0.5)

plt.xticks(xtick, fontname='arial', fontsize = fsize, fontweight = fweight)
ax2.xaxis.set_tick_params(which = 'major', direction = 'in', width = 1, length = 6, pad = 5, bottom = True, top = True)
ax2.xaxis.set_minor_locator(xminorLocator)
ax2.xaxis.set_tick_params(which = 'minor', direction = 'in', width = 1, length = 3, pad = 5, bottom = True, top = True)

plt.yticks(ytick_temp, fontname='arial', fontsize = fsize, fontweight = fweight)
ax2.yaxis.set_tick_params(which = 'major', direction = 'in', width = 1, length = 6, pad = 5, left = True, right = True)
ax2.yaxis.set_minor_locator(yminorLocator_temp)
ax2.yaxis.set_tick_params(which = 'minor', direction = 'in', width = 1, length = 3, pad = 5, left = True, right = True)

ax2.set_xlim(xlim)
ax2.set_ylim(ylim_temp)

######################################################

ax3 = plt.axes(ax3_pos)

plt.plot( redshift_list, dens_history[:,0,0], linestyle='-', color='g', linewidth=2 )
plt.plot( redshift_list, dens_history[:,1,0], linestyle='-', color='b', linewidth=2 )
plt.plot( redshift_list, dens_history[:,2,0], linestyle='-', color='r', linewidth=2 )
if NZF>3:
        plt.plot( redshift_list, dens_history[:,3,0], linestyle='-', color='k', linewidth=2 )

plt.plot( redshift_list, dens_history[:,0,0]+dens_history[:,0,1], linestyle='--', color='g', linewidth=1 )
plt.plot( redshift_list, dens_history[:,1,0]+dens_history[:,1,1], linestyle='--', color='b', linewidth=1 )
plt.plot( redshift_list, dens_history[:,2,0]+dens_history[:,2,1], linestyle='--', color='r', linewidth=1 )
if NZF>3:
        plt.plot( redshift_list, dens_history[:,3,0]+dens_history[:,3,1], linestyle='--', color='k', linewidth=1 )

plt.plot( redshift_list, dens_history[:,0,0]-dens_history[:,0,1], linestyle='--', color='g', linewidth=1 )
plt.plot( redshift_list, dens_history[:,1,0]-dens_history[:,1,1], linestyle='--', color='b', linewidth=1 )
plt.plot( redshift_list, dens_history[:,2,0]-dens_history[:,2,1], linestyle='--', color='r', linewidth=1 )
if NZF>3:
        plt.plot( redshift_list, dens_history[:,3,0]-dens_history[:,3,1], linestyle='--', color='k', linewidth=1 )

xlabel(xstr, fontname='times', fontsize = fsize2, fontweight = fweight)
ylabel(ystr_dens, fontname='times', fontsize = fsize2, fontweight = fweight)
ax3.xaxis.set_label_coords(0.5, xlaby)
ax3.yaxis.set_label_coords(ylabx, 0.5)

plt.xticks(xtick, fontname='arial', fontsize = fsize, fontweight = fweight)
ax3.xaxis.set_tick_params(which = 'major', direction = 'in', width = 1, length = 6, pad = 5, bottom = True, top = True)
ax3.xaxis.set_minor_locator(xminorLocator)
ax3.xaxis.set_tick_params(which = 'minor', direction = 'in', width = 1, length = 3, pad = 5, bottom = True, top = True)

plt.yticks(ytick_dens, fontname='arial', fontsize = fsize, fontweight = fweight)
ax3.yaxis.set_tick_params(which = 'major', direction = 'in', width = 1, length = 6, pad = 5, left = True, right = True)
ax3.yaxis.set_minor_locator(yminorLocator_dens)
ax3.yaxis.set_tick_params(which = 'minor', direction = 'in', width = 1, length = 3, pad = 5, left = True, right = True)

ax3.set_xlim(xlim)
ax3.set_ylim(ylim_dens)

######################################################

savefig(fig_tit_low_Z)
print(fig_tit_low_Z)
plt.close()

########################################################################################################################################################################################################################
########################################################################################################################################################################################################################
