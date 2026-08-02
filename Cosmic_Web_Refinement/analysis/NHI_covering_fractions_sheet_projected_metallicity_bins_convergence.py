from gadget import *
from gadget_subfind import *
import const as c
import matplotlib as mpl
import os, sys

######################################################
### Input parameters ###

exclude = 5	# 0: Include all cells;	1: Exclude allcels in FoF groups; 2: Exclude all cells in Subfind subhalos; 3: Exclude all cells within Rvir of central halos if FoF groups; 4: Exclude halos in projection
min_mass  = 1.0e9	# Minimum halo mas for exclusion rules
min_mass /= 1.0e10	# halo masses are given in units of 10^{10} solar masses
min_Npart = 32		# Minimum number of DM particles in FoF groups that are excluded if exclude==5
Dweight = 0	# 1 for density-weighted average metallicity, 0 for mass-weighted average metallicity
Z_from_ZF1 = 1 # Divide surface of each simulation based on metallicity distribution in ZF1 to make sure to probe same volume.

ZFs   = np.array([1, 0, 2, 4])
NZF = 4

SnapNums  = np.array([91, 99, 108])
int_reds  = np.array([5, 4, 3])
Nsnap = 3
redshift_list = np.zeros([Nsnap])

#xmaxs = np.array([30.0, 30.0, 60.0])		# proper kpc --> Integration depth is +/- xmax from midplane of sheet
xmax = 100.0								# proper kpc --> Integration depth is +/- xmax from midplane of sheet
zL   = 1800	# comoving h^{-1} kpc, for sheet boundary (+/-zL)
yL   = 1000	# comoving h^{-1} kpc, for sheet boundary (+/-yL)

#Zmaxs = np.array([-4.0, -3.0, -2.0, -1.0, 0.0])	# Bins of metalicity with these boundaries
#Zmaxs = np.array([-3.5, -2.5, -1.5, -0.5, 0.5])	# Bins of metalicity with these boundaries
#Zmaxs = np.array([-3.0, -1.5])						# Bins of metalicity with these boundaries
Zmaxs = np.array([-3.0])							# Bins of metalicity with these boundaries
NZ_bins = np.size(Zmaxs)

######################################################
### Constants and units ###

gravity = c.G                        #cm^3/(g s^2)
msolar = c.msol                      #gram
mproton = c.amu                      #gram
kB = c.KB                            #erg/K
kpcincm = c.parsec*1e3
Xh = 0.76
gamma_gas = 5.0/3.0
Omegab = 0.0486
Omega0 = 0.3089
HUBBLE = 3.2407789e-18 			# 100 km/s/Mpc in 1/sec, i.e. H0 in h/sec
cm_in_pc     = 3.085678e18
gr_in_Ms     = 1.989e33

UnitLength_in_cm         = 3.085678e21	# kpc/h
UnitMass_in_g            = 1.989e43	# 1.0e10 Msun/h
UnitVelocity_in_cm_per_s = 1.0e5	# 1 km/sec

UnitTime_in_s       = UnitLength_in_cm / UnitVelocity_in_cm_per_s 
UnitDensity_in_cgs  = UnitMass_in_g    / UnitLength_in_cm**3
UnitNumberDensity_in_cgs = UnitDensity_in_cgs / mproton			# ntot in cm^{-3}
UnitPressure_in_cgs = UnitMass_in_g    / UnitLength_in_cm    / UnitTime_in_s**2
UnitEnergy_in_cgs   = UnitMass_in_g    * UnitLength_in_cm**2 / UnitTime_in_s**2

######################################################
### NHI array and plot parameters ###

NHI_max  = 22.0
NHI_min  = 9.0
dNHI     = 0.04
Nbin_NHI = np.int( (NHI_max - NHI_min) / dNHI ) + 1
NHI_bins = NHI_min + np.arange(0,Nbin_NHI)*dNHI
NHI_CDF = np.zeros([Nbin_NHI, NZ_bins+2, NZF, Nsnap])
XHI     = np.zeros([NZ_bins+2, NZF, Nsnap, 2])

xtick_NHI         = np.arange(NHI_min,NHI_max+1.0,1.0)
xminorLocator_NHI = MultipleLocator(0.2)
xstr_NHI = r"${\rm Log}(N_{\rm HI}\:{\rm [cm^{-2}]})$"
xmin_NHI = 12.2
xmax_NHI = 20.0

ymin  = -3.0
ymax  = 0.0
dy    = 1.0
ytick = 10.0**(np.arange(ymin,ymax+dy,dy))
ystr = r"$f_{\rm A}(>N_{\rm HI})$"

xsize   = 15.9
ysize   = 5.3
fsize   = 18
fsize2  = 22
fweight = 150
ax1_pos  = [0.100, 0.12, 0.277, 0.80]
ax2_pos  = [0.377, 0.12, 0.277, 0.80]
ax3_pos  = [0.654, 0.12, 0.277, 0.80]

xlaby = -0.07
ylabx = -0.12
tit_pos = [0.5, 1.005]

######################################################
### Read in snapshot and rotate frame ###

FigDir   = "../plots/FilB-Zoom4"
if not os.path.exists(FigDir):
	os.mkdir(FigDir)
FigDir   = "%s/for_letter" % (FigDir)
if not os.path.exists(FigDir):
	os.mkdir(FigDir)
FigDir   = "%s/HI_covering_fraction_Z_bins" % (FigDir)
if not os.path.exists(FigDir):
	os.mkdir(FigDir)
if Z_from_ZF1==1:
	FigDir   = "%s/Z_from_ZF1" % (FigDir)
	if not os.path.exists(FigDir):
		os.mkdir(FigDir)

if Dweight==0:
	if exclude<5:
		fname = FigDir+"/projected_convergence_xmax_%03d_exclude_%01d_Mmin_%01d.png" % (np.int(xmax),exclude,np.log10(min_mass)+10)
	elif exclude==5:
		fname = FigDir+"/projected_convergence_xmax_%03d_exclude_%01d_Nmin_%02d.png" % (np.int(xmax),exclude,min_Npart)
elif Dweight==1:
	if exclude<5:
		fname = FigDir+"/projected_convergence_xmax_%03d_exclude_%01d_Mmin_%01d_Dweight.png" % (np.int(xmax),exclude,np.log10(min_mass)+10)
	elif exclude==5:
		fname = FigDir+"/projected_convergence_xmax_%03d_exclude_%01d_Nmin_%02d_Dweight.png" % (np.int(xmax),exclude,min_Npart)

print(fname)

for snap_ind in np.arange(0,Nsnap):

	SnapNum = SnapNums[snap_ind]
	int_red = int_reds[snap_ind]
	#boxz = xmaxs[snap_ind]
	print('SnapNum ',SnapNum)

	if Dweight==0:
		if exclude<5:
			out_txt_name = FigDir+"/projected_output_snap_%03d_xmax_%03d_exclude_%01d_Mmin_%01d.txt" % (SnapNum,np.int(xmax),exclude,np.log10(min_mass)+10)
			out_xhi_name = FigDir+"/XHI_snap_%03d_xmax_%03d_exclude_%01d_Mmin_%01d.txt" % (SnapNum,np.int(xmax),exclude,np.log10(min_mass)+10)
		elif exclude==5:
			out_txt_name = FigDir+"/projected_output_snap_%03d_xmax_%03d_exclude_%01d_Nmin_%02d.txt" % (SnapNum,np.int(xmax),exclude,min_Npart)
			out_xhi_name = FigDir+"/XHI_snap_%03d_xmax_%03d_exclude_%01d_Nmin_%02d.txt" % (SnapNum,np.int(xmax),exclude,min_Npart)
	elif Dweight==1:
		if exclude<5:
			out_txt_name = FigDir+"/projected_output_snap_%03d_xmax_%03d_exclude_%01d_Mmin_%01d_Dweight.txt" % (SnapNum,np.int(xmax),exclude,np.log10(min_mass)+10)
			out_xhi_name = FigDir+"/XHI_snap_%03d_xmax_%03d_exclude_%01d_Mmin_%01d.txt" % (SnapNum,np.int(xmax),exclude,np.log10(min_mass)+10)
		elif exclude==5:
			out_txt_name = FigDir+"/projected_output_snap_%03d_xmax_%03d_exclude_%01d_Nmin_%02d_Dweight.txt" % (SnapNum,np.int(xmax),exclude,min_Npart)
			out_xhi_name = FigDir+"/XHI_snap_%03d_xmax_%03d_exclude_%01d_Nmin_%02d.txt" % (SnapNum,np.int(xmax),exclude,min_Npart)

	print(out_txt_name)

	for ZF_ind in np.arange(0,NZF):
		ZF = ZFs[ZF_ind]
		if ZF_ind==0:	# I rearanged the ZFs array to have ZF1 first, in order to use its metallicity distribution for all ZFs. This helps store the output in the correct order
			ZF_ind2 = 1
		elif ZF_ind==1:
			ZF_ind2 = 0
		else:
			ZF_ind2 = ZF_ind 

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

		s  = gadget_readsnap( SnapNum2, SimDir, SnapBase, loadonlytype=[0], loadonly=['pos', 'rho', 'mass', 'gz', 'nh', 'gmet'], lazy_load=True )
		sf = load_subfind( SnapNum2, dir=SimDir+'/' )
		redshift = s.redshift
		if ZF==2:
			redshift_list[snap_ind] = redshift
		BoxHalf = 0.5 * s.boxsize
		aexp = s.time
		h0 = s.hubbleparam
		BoxHalf *= aexp/h0
		BoxSize = 2.0*BoxHalf
		RhoMeanBaryons = Omegab * 3.0 * (h0*HUBBLE)**2 / (8.0 * pi * gravity * (aexp**3))	# mean baryon density (physical gr/cm^3)

		Dens = s.rho * UnitDensity_in_cgs					# density in gr cm^{-3}
		s.data['rho'] = s.rho * s.gmet[:,0] * UnitNumberDensity_in_cgs	# nH  cm^{-3}
		s.data['rho_HI'] = s.data['rho'] * s.nh					# nHI cm^{-3}
		s.data['gz']    /= 0.0127	# Solar units
		s.data['MZ']     = s.mass * s.gz
		s.data['DZ']     = s.rho * s.gz

		if exclude==1:		### Remove gas cells associated with FoF groups ###
			print("exclude 1")
			Ng = sf.totngroups
			first_cell = 0
			Ntrack_g = (Ng - np.mod(Ng,10)) / 10
			for g_ind in np.arange(0,Ng):
				if( (np.mod(g_ind,Ntrack_g)==0) | (g_ind==0) ):
					print("g_ind, Ng= ", g_ind, Ng)

				Ncell = sf.data['flty'][g_ind,0]
				last_cell = first_cell+Ncell
				if sf.data['fmas'][g_ind]>=min_mass:
					s.data['rho_HI'][first_cell:last_cell] = 0.0
					s.data['rho'][first_cell:last_cell] = 0.0
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
					s.data['rho_HI'][first_cell:last_cell] = 0.0
					s.data['rho'][first_cell:last_cell] = 0.0
				first_cell += Ncell
		elif exclude==2:	### Remove gas cells associated with subfind subhalos with M>min_mass ###
			print("exclude 2")
			Ng = sf.totngroups
			first_cell = 0
			Ntrack_g = (Ng - np.mod(Ng,10)) / 10
			for g_ind in np.arange(0,Ng):
				first_sub = sf.data['ffsh'][g_ind]
				Nsub = sf.data['fnsh'][g_ind]
				if( (np.mod(g_ind,Ntrack_g)==0) | (g_ind==0) ):
					print("g_ind, Ng, first_sub, Nsub= ", g_ind, Ng, first_sub, Nsub)
				if first_sub>-1:
					for s_ind in np.arange(0,Nsub):
						Ncell = sf.data['slty'][first_sub+s_ind,0]
						last_cell = first_cell+Ncell
						if sf.data['smas'][first_sub+s_ind]>min_mass:
							s.data['rho_HI'][first_cell:last_cell] = 0.0
							s.data['rho'][first_cell:last_cell] = 0.0
						first_cell += Ncell
				first_cell = np.sum(sf.data['flty'][0:g_ind,0])

		######################################################

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

		xx *= h0/aexp # comiving h^-1 kpc
		yy *= h0/aexp # comiving h^-1 kpc
		zz *= h0/aexp # comiving h^-1 kpc

		################################

		if exclude==3:	### Remove FoF groups with M_{fof}>10^7 ###
			b1 = np.where( (np.abs(zz)<=zL) & (np.abs(yy)<=yL) & (np.abs(xx)<=xmax*h0/aexp) )	### Select sheet region ###
			Ncell1 = np.size(b1)
			print('Ncell1=',Ncell1)
			xx1   = xx[b1]*(aexp/h0)	# proper kpc
			yy1   = yy[b1]*(aexp/h0)	# proper kpc
			zz1   = zz[b1]*(aexp/h0)	# proper kpc
			gx1 = gxx
			gy1 = gyy
			gz1 = gzz

			bg2 = np.where( ((gz1+gR)>=-zL*(aexp/h0)) & ((gz1-gR)<=zL*(aexp/h0)) & ((gy1+gR)>=-yL*(aexp/h0)) & ((gy1-gR)<=yL*(aexp/h0)) & ((gx1+gR)>=-xmax) & ((gx1-gR)<=xmax) & (gM>=min_mass) & (gR>0.0) )
			Ng2 = np.size(bg2)
			if Ng2>0:
				gx2 = gx1[bg2]
				gy2 = gy1[bg2]
				gz2 = gz1[bg2]
				gR2 = gR[bg2]
				gM2 = gM[bg2]
				print("Ng2=",Ng2)
				Ntrack_g2 = (Ng2 - np.mod(Ng2,10)) / 10
				for g_ind in np.arange(0,Ng2):
					if( (np.mod(g_ind,Ntrack_g2)==0) | (g_ind==0) ):
						print("g_ind, Ng2= ",g_ind,Ng2)
					b3, = np.where( np.sqrt( (xx1-gx2[g_ind])**2 + (yy1-gy2[g_ind])**2 + (zz1-gz2[g_ind])**2 ) <= gR2[g_ind] )
					s.data['rho_HI'][b1[b3]] = 0.0

	#			min_dist = np.zeros(Ncell1)
	#			for cell_ind in np.arange(0,Ncell1):
	#				min_dist[cell_ind] = np.min( np.sqrt( (xx1[cell_ind]-gx2)**2 + (yy1[cell_ind]-gy2)**2 + (zz1[cell_ind]-gz2)**2 ) / gR2 )
	#			b3 = np.where( min_dist>1.0 )
	#			Ncell3 = np.size(b3)
	#			b4 = np.where( min_dist<1.0 )
	#			s.data['rho_HI'][b4] = 0.0
				b3 = np.where( (np.abs(zz)<=zL) & (np.abs(yy)<=yL) & (np.abs(xx)<=xmax*h0/aexp) & (s.data['rho_HI']>0.0) )
				Ncell3 = np.size(b3)
				print( 'Final Ncell no halos=',Ncell3 )
			else:
				print( 'Final Ncell no halos=',Ncell1 )

		######################

		s.data['pos'][:,0] = xx			# comiving h^-1 kpc
		s.data['pos'][:,1] = yy			# comiving h^-1 kpc
		s.data['pos'][:,2] = zz			# comiving h^-1 kpc

		### For computing column density map
		imgcent = [0.,0.,0.]
		pixres = 2.0			# comiving h^-1 kpc
		boxx = 2.0*zL			# comiving h^-1 kpc
		boxy = 2.0*yL			# comiving h^-1 kpc
		boxz = 2.0*xmax*(h0/aexp)	# comiving h^-1 kpc
		VT = boxx*boxy
		nx = int(boxx/pixres)
		ny = int(boxy/pixres)

		print( boxx, boxy, boxz, pixres )
		print( int(boxx/pixres), int(boxy/pixres) )

		################################
		### Create surface maps ###

		print( np.max(s.data['rho_HI']), np.min(s.data['rho_HI']) )
		NHI_proj = s.get_Aslice("rho_HI", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz,  axes = [2,1], proj = True, numthreads=32)
		map_NHI  = np.transpose( NHI_proj['grid'] ) * (pixres*aexp/h0) * UnitLength_in_cm	# NHI cm^{-2}
		xpl   = NHI_proj['x']		# comoving kpc/h
		ypl   = NHI_proj['y']		# comoving kpc/h
		print(np.max(xpl), np.min(xpl))
		print(np.max(ypl), np.min(ypl))
		b3 = np.where( map_NHI>0.0 )
		print( np.max(map_NHI), np.min(map_NHI), np.min(map_NHI[b3]) )

		NH_proj = s.get_Aslice("rho", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz,  axes = [2,1], proj = True, numthreads=32)
		map_NH  = np.transpose( NH_proj['grid'] ) * (pixres*aexp/h0) * UnitLength_in_cm	# NH cm^{-2}

		if Dweight==0:
			Z_proj_Q = s.get_Aslice("MZ", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz,  axes = [2,1], proj = True, numthreads=32)
			Z_proj_W = s.get_Aslice("mass", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz,  axes = [2,1], proj = True, numthreads=32)
		elif Dweight==1:
			Z_proj_Q = s.get_Aslice("DZ", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz,  axes = [2,1], proj = True, numthreads=32)
			Z_proj_W = s.get_Aslice("rho", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz,  axes = [2,1], proj = True, numthreads=32)

		map_Z    = np.transpose( Z_proj_Q['grid'] / Z_proj_W['grid'] )
		if Z_from_ZF1==1:
			if ZF==1:
				map_Z2 = np.zeros(np.shape(map_Z))
				map_Z2[:] = map_Z[:]
		else:
			map_Z2 = np.zeros(np.shape(map_Z))
			map_Z2[:] = map_Z[:]

		################################
		### Make histograms ###

		nx = int(boxx/pixres)
		ny = int(boxy/pixres)

		gxx *= (h0/aexp)		# comoving kpc/h
		gyy *= (h0/aexp)		# comoving kpc/h
		gzz *= (h0/aexp)		# comoving kpc/h
		gR  *= (h0/aexp)		# comoving kpc/h
		if exclude==4:	### Exclude halos in projection
			b2 = np.where( ((gzz+gR)>=-0.5*boxx) & ((gzz-gR)<=0.5*boxx) & ((gyy+gR)>=-0.5*boxy) & ((gyy-gR)<=0.5*boxy) & ((gxx+gR)>=-0.5*boxz) & ((gxx-gR)<=0.5*boxz) & (gM>=min_mass) & (gR>=0.0) )
			Ng2 = np.size(b2)
			print("Ng2=",Ng2)
			gz2 = gzz[b2]
			gy2 = gyy[b2]
			gR2 = gR[b2]
			gM2 = gM[b2]
			if Ng2>0:
				print("gR2=",np.min(gR2),np.max(gR2))
				print("gM2=",np.min(gM2),np.max(gM2))
		else:
			Ng2 = 0		# To avoid the loop condition

		Ntrack_x = (nx - np.mod(nx,10)) / 10
		Ntrack_y = (ny - np.mod(ny,10)) / 10
		for x_ind in np.arange(0,nx):
			for y_ind in np.arange(0,ny):
				if( ((np.mod(x_ind,Ntrack_x)==0) & (np.mod(y_ind,Ntrack_y)==0)) | ((x_ind==0) & (y_ind==0)) ):
					print("(x_ind,y_ind) of (nx,ny) ",x_ind,y_ind,nx,ny)
				g_ind = 0
				while g_ind<Ng2:	# Check if cell is within projected Rvir of any halo with Mvir>min_mass
					if( np.sqrt( (xpl[x_ind]-gz2[g_ind])**2 + (ypl[y_ind]-gy2[g_ind])**2 ) <= gR2[g_ind] ):
						g_ind = 2*Ng2+10
					else:
						g_ind += 1
				if g_ind==Ng2:		# Cell is not within projected Rvir of any halo with Mvir>min_mass
					if map_NHI[y_ind,x_ind]>0.0:
						NHI_1 = np.log10(map_NHI[y_ind,x_ind])
						iNHI = np.int(np.floor( (NHI_1 - NHI_min) / dNHI ))
						iNHI = np.max([iNHI,0])
						iNHI = np.min([iNHI,Nbin_NHI-1])
						NHI_CDF[0:iNHI, 0, ZF_ind2, snap_ind] += pixres**2
						XHI[0, ZF_ind2, snap_ind, 0] += map_NHI[y_ind,x_ind]
						XHI[0, ZF_ind2, snap_ind, 1] += map_NH[y_ind,x_ind]

						Z_ind = 0
						while Z_ind<NZ_bins:
							if( np.log10(map_Z2[y_ind,x_ind]) < Zmaxs[Z_ind] ):
								NHI_CDF[0:iNHI, Z_ind+1, ZF_ind2, snap_ind] += pixres**2
								XHI[Z_ind+1, ZF_ind2, snap_ind, 0] += map_NHI[y_ind,x_ind]
								XHI[Z_ind+1, ZF_ind2, snap_ind, 1] += map_NH[y_ind,x_ind]
								Z_ind = 2*NZ_bins + 10
							else:
								Z_ind += 1
						if Z_ind==NZ_bins:
							NHI_CDF[0:iNHI, Z_ind+1, ZF_ind2, snap_ind] += pixres**2
							XHI[Z_ind+1, ZF_ind2, snap_ind, 0] += map_NHI[y_ind,x_ind]
							XHI[Z_ind+1, ZF_ind2, snap_ind, 1] += map_NH[y_ind,x_ind]

		for Z_ind in np.arange(0,NZ_bins+2):
			NHI_CDF[:, Z_ind, ZF_ind2, snap_ind] /= VT
			XHI[Z_ind, ZF_ind2, snap_ind, 0] /= XHI[Z_ind, ZF_ind2, snap_ind, 1]

		gxx /= (h0/aexp)		# comoving kpc/h
		gyy /= (h0/aexp)		# comoving kpc/h
		gzz /= (h0/aexp)		# comoving kpc/h
		gR  /= (h0/aexp)		# comoving kpc/h

			########################
	output_array = np.column_stack([NHI_bins, NHI_CDF[:, 0, 0, snap_ind], NHI_CDF[:, 1, 0, snap_ind], NHI_CDF[:, 2, 0, snap_ind], 
								              NHI_CDF[:, 0, 1, snap_ind], NHI_CDF[:, 1, 1, snap_ind], NHI_CDF[:, 2, 1, snap_ind], 
											  NHI_CDF[:, 0, 2, snap_ind], NHI_CDF[:, 1, 2, snap_ind], NHI_CDF[:, 2, 2, snap_ind], 
											  NHI_CDF[:, 0, 3, snap_ind], NHI_CDF[:, 1, 3, snap_ind], NHI_CDF[:, 2, 3, snap_ind]] )
	np.savetxt(out_txt_name, output_array)

	xhi_array = np.zeros([4,3])
	xhi_array[0,:] = XHI[:, 0, snap_ind, 0]
	xhi_array[1,:] = XHI[:, 1, snap_ind, 0]
	xhi_array[2,:] = XHI[:, 2, snap_ind, 0]
	xhi_array[3,:] = XHI[:, 3, snap_ind, 0]
	np.savetxt(out_xhi_name, xhi_array)

######################################################
### Open figure object ###

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig = plt.figure(figsize = (xsize,ysize))

disp01 = r"${\rm ZF0.5}$"
disp02 = r"${\rm ZF1.0}$"
disp03 = r"${\rm ZF2.0}$"
disp04 = r"${\rm ZF4.0}$"

disp11 = r"${\rm total}$"
disp12 = r"${\rm [Z]}<"+'{0:.1f}'.format(Zmaxs[0])+"$"
disp13 = r"${\rm "+'{0:.1f}'.format(Zmaxs[0])+"<[Z]}$"

tit0 = r"$z="+'{0:.2f}'.format(redshift_list[0])+"$"
tit1 = r"$z="+'{0:.2f}'.format(redshift_list[1])+"$"
tit2 = r"$z="+'{0:.2f}'.format(redshift_list[2])+"$"

disp21 = r"${\rm log}(x_{\rm HI})="+'{0:.2f}'.format(np.log10(XHI[0,0,0,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[1,0,0,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[2,0,0,0]))+"$"
disp22 = r"${\rm log}(x_{\rm HI})="+'{0:.2f}'.format(np.log10(XHI[0,1,0,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[1,1,0,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[2,1,0,0]))+"$"
disp23 = r"${\rm log}(x_{\rm HI})="+'{0:.2f}'.format(np.log10(XHI[0,2,0,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[1,2,0,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[2,2,0,0]))+"$"
disp24 = r"${\rm log}(x_{\rm HI})="+'{0:.2f}'.format(np.log10(XHI[0,3,0,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[1,3,0,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[2,3,0,0]))+"$"

disp31 = r"${\rm log}(x_{\rm HI})="+'{0:.2f}'.format(np.log10(XHI[0,0,1,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[1,0,1,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[2,0,1,0]))+"$"
disp32 = r"${\rm log}(x_{\rm HI})="+'{0:.2f}'.format(np.log10(XHI[0,1,1,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[1,1,1,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[2,1,1,0]))+"$"
disp33 = r"${\rm log}(x_{\rm HI})="+'{0:.2f}'.format(np.log10(XHI[0,2,1,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[1,2,1,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[2,2,1,0]))+"$"
disp34 = r"${\rm log}(x_{\rm HI})="+'{0:.2f}'.format(np.log10(XHI[0,3,1,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[1,3,1,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[2,3,1,0]))+"$"

disp41 = r"${\rm log}(x_{\rm HI})="+'{0:.2f}'.format(np.log10(XHI[0,0,2,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[1,0,2,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[2,0,2,0]))+"$"
disp42 = r"${\rm log}(x_{\rm HI})="+'{0:.2f}'.format(np.log10(XHI[0,1,2,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[1,1,2,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[2,1,2,0]))+"$"
disp43 = r"${\rm log}(x_{\rm HI})="+'{0:.2f}'.format(np.log10(XHI[0,2,2,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[1,2,2,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[2,2,2,0]))+"$"
disp44 = r"${\rm log}(x_{\rm HI})="+'{0:.2f}'.format(np.log10(XHI[0,3,2,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[1,3,2,0]))+",\,"+'{0:.2f}'.format(np.log10(XHI[2,3,2,0]))+"$"

################################################################
### Make plots ###

ax1 = plt.axes(ax1_pos)

line01, = plt.plot( NHI_bins, NHI_CDF[:, 0, 0, 0], marker='None', linestyle='-', color='g', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 1, 0, 0], marker='None', linestyle='--', color='g', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 2, 0, 0], marker='None', linestyle=':', color='g', linewidth=2 )

line02, = plt.plot( NHI_bins, NHI_CDF[:, 0, 1, 0], marker='None', linestyle='-', color='b', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 1, 1, 0], marker='None', linestyle='--', color='b', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 2, 1, 0], marker='None', linestyle=':', color='b', linewidth=2 )

line03, = plt.plot( NHI_bins, NHI_CDF[:, 0, 2, 0], marker='None', linestyle='-', color='r', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 1, 2, 0], marker='None', linestyle='--', color='r', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 2, 2, 0], marker='None', linestyle=':', color='r', linewidth=2 )

line04, = plt.plot( NHI_bins, NHI_CDF[:, 0, 3, 0], marker='None', linestyle='-', color='k', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 1, 3, 0], marker='None', linestyle='--', color='k', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 2, 3, 0], marker='None', linestyle=':', color='k', linewidth=2 )

yvec = np.logspace(-5.0, 1.0, 30)
xvec = 17.2*np.ones(30)
plt.plot( xvec, yvec, marker='None', linestyle=':', color='k', linewidth=1 )
ax1.set_yscale('log')

xlabel(xstr_NHI, fontname='times', fontsize = fsize2, fontweight = fweight)
ylabel(ystr,     fontname='times', fontsize = fsize2, fontweight = fweight)
ax1.xaxis.set_label_coords(0.5, xlaby)
ax1.yaxis.set_label_coords(ylabx, 0.5)

plt.xticks(xtick_NHI,  fontname='arial', fontsize = fsize, fontweight = fweight)
ax1.xaxis.set_tick_params(which = 'major', direction = 'in', width = 1, length = 6, pad = 5, bottom = True, top = True)
ax1.xaxis.set_minor_locator(xminorLocator_NHI)
ax1.xaxis.set_tick_params(which = 'minor', direction = 'in', width = 1, length = 3, pad = 5, bottom = True, top = True)

plt.yticks(ytick, fontname='arial', fontsize = fsize, fontweight = fweight)
ax1.yaxis.set_tick_params(which = 'major', direction = 'in', width = 1, length = 6, pad = 5, left = True, right = True)
yminorLocator = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=10)
ax1.yaxis.set_minor_locator(yminorLocator)
ax1.yaxis.set_tick_params(which = 'minor', direction = 'in', width = 1, length = 3, pad = 5, left = True, right = True)

ttl = ax1.set_title(tit0, fontname='times', fontsize = fsize2, fontweight = fweight, loc='center')
ttl.set_position(tit_pos)

plt.legend([line01, line02, line03, line04],[disp01,disp02,disp03,disp04], fontsize = 14, frameon=False, loc='lower left')

plt.text(12.5, 0.0600, disp21, horizontalalignment='left', verticalalignment='top', fontname='times', fontsize = 14, transform = ax1.transData, color = 'g')
plt.text(12.5, 0.0363, disp22, horizontalalignment='left', verticalalignment='top', fontname='times', fontsize = 14, transform = ax1.transData, color = 'b')
plt.text(12.5, 0.0220, disp23, horizontalalignment='left', verticalalignment='top', fontname='times', fontsize = 14, transform = ax1.transData, color = 'r')
plt.text(12.5, 0.0133, disp24, horizontalalignment='left', verticalalignment='top', fontname='times', fontsize = 14, transform = ax1.transData, color = 'k')

ax1.set_xlim([xmin_NHI, xmax_NHI])
ax1.set_ylim([10.0**ymin, 10.0**ymax])

#######################

ax2 = plt.axes(ax2_pos)

plt.plot( NHI_bins, NHI_CDF[:, 0, 0, 1], marker='None', linestyle='-', color='g', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 1, 0, 1], marker='None', linestyle='--', color='g', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 2, 0, 1], marker='None', linestyle=':', color='g', linewidth=2 )

plt.plot( NHI_bins, NHI_CDF[:, 0, 1, 1], marker='None', linestyle='-', color='b', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 1, 1, 1], marker='None', linestyle='--', color='b', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 2, 1, 1], marker='None', linestyle=':', color='b', linewidth=2 )

plt.plot( NHI_bins, NHI_CDF[:, 0, 2, 1], marker='None', linestyle='-', color='r', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 1, 2, 1], marker='None', linestyle='--', color='r', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 2, 2, 1], marker='None', linestyle=':', color='r', linewidth=2 )

line11, = plt.plot( NHI_bins, NHI_CDF[:, 0, 3, 1], marker='None', linestyle='-', color='k', linewidth=2 )
line12, = plt.plot( NHI_bins, NHI_CDF[:, 1, 3, 1], marker='None', linestyle='--', color='k', linewidth=2 )
line13, = plt.plot( NHI_bins, NHI_CDF[:, 2, 3, 1], marker='None', linestyle=':', color='k', linewidth=2 )

yvec = np.logspace(-5.0, 1.0, 30)
xvec = 17.2*np.ones(30)
plt.plot( xvec, yvec, marker='None', linestyle=':', color='k', linewidth=1 )
ax2.set_yscale('log')

xlabel(xstr_NHI, fontname='times', fontsize = fsize2, fontweight = fweight)
ax2.xaxis.set_label_coords(0.5, xlaby)

plt.xticks(xtick_NHI,  fontname='arial', fontsize = fsize, fontweight = fweight)
ax2.xaxis.set_tick_params(which = 'major', direction = 'in', width = 1, length = 6, pad = 5, bottom = True, top = True)
ax2.xaxis.set_minor_locator(xminorLocator_NHI)
ax2.xaxis.set_tick_params(which = 'minor', direction = 'in', width = 1, length = 3, pad = 5, bottom = True, top = True)

plt.yticks(ytick, fontname='arial', fontsize = fsize, fontweight = fweight)
ax2.yaxis.set_tick_params(which = 'major', direction = 'in', width = 1, length = 6, pad = 5, left = True, right = True)
yminorLocator = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=10)
ax2.yaxis.set_minor_locator(yminorLocator)
ax2.yaxis.set_tick_params(which = 'minor', direction = 'in', width = 1, length = 3, pad = 5, left = True, right = True)

labels = [item.get_text() for item in ax2.get_yticklabels()]
empty_string_labels = ['']*len(labels)
ax2.set_yticklabels(empty_string_labels)

ttl = ax2.set_title(tit1, fontname='times', fontsize = fsize2, fontweight = fweight, loc='center')
ttl.set_position(tit_pos)

plt.legend([line11, line12, line13],[disp11,disp12,disp13], fontsize = 14, frameon=False, loc='lower left')

plt.text(12.5, 0.0600, disp31, horizontalalignment='left', verticalalignment='top', fontname='times', fontsize = 14, transform = ax2.transData, color = 'g')
plt.text(12.5, 0.0363, disp32, horizontalalignment='left', verticalalignment='top', fontname='times', fontsize = 14, transform = ax2.transData, color = 'b')
plt.text(12.5, 0.0220, disp33, horizontalalignment='left', verticalalignment='top', fontname='times', fontsize = 14, transform = ax2.transData, color = 'r')
plt.text(12.5, 0.0133, disp34, horizontalalignment='left', verticalalignment='top', fontname='times', fontsize = 14, transform = ax2.transData, color = 'k')

ax2.set_xlim([xmin_NHI, xmax_NHI])
ax2.set_ylim([10.0**ymin, 10.0**ymax])

#######################

ax3 = plt.axes(ax3_pos)

plt.plot( NHI_bins, NHI_CDF[:, 0, 0, 2], marker='None', linestyle='-', color='g', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 1, 0, 2], marker='None', linestyle='--', color='g', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 2, 0, 2], marker='None', linestyle=':', color='g', linewidth=2 )

plt.plot( NHI_bins, NHI_CDF[:, 0, 1, 2], marker='None', linestyle='-', color='b', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 1, 1, 2], marker='None', linestyle='--', color='b', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 2, 1, 2], marker='None', linestyle=':', color='b', linewidth=2 )

plt.plot( NHI_bins, NHI_CDF[:, 0, 2, 2], marker='None', linestyle='-', color='r', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 1, 2, 2], marker='None', linestyle='--', color='r', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 2, 2, 2], marker='None', linestyle=':', color='r', linewidth=2 )

plt.plot( NHI_bins, NHI_CDF[:, 0, 3, 2], marker='None', linestyle='-', color='k', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 1, 3, 2], marker='None', linestyle='--', color='k', linewidth=2 )
plt.plot( NHI_bins, NHI_CDF[:, 2, 3, 2], marker='None', linestyle=':', color='k', linewidth=2 )

yvec = np.logspace(-5.0, 1.0, 30)
xvec = 17.2*np.ones(30)
plt.plot( xvec, yvec, marker='None', linestyle=':', color='k', linewidth=1 )
ax3.set_yscale('log')

xlabel(xstr_NHI, fontname='times', fontsize = fsize2, fontweight = fweight)
ax3.xaxis.set_label_coords(0.5, xlaby)

plt.xticks(xtick_NHI,  fontname='arial', fontsize = fsize, fontweight = fweight)
ax3.xaxis.set_tick_params(which = 'major', direction = 'in', width = 1, length = 6, pad = 5, bottom = True, top = True)
ax3.xaxis.set_minor_locator(xminorLocator_NHI)
ax3.xaxis.set_tick_params(which = 'minor', direction = 'in', width = 1, length = 3, pad = 5, bottom = True, top = True)

plt.yticks(ytick, fontname='arial', fontsize = fsize, fontweight = fweight)
ax3.yaxis.set_tick_params(which = 'major', direction = 'in', width = 1, length = 6, pad = 5, left = True, right = True)
yminorLocator = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=10)
ax3.yaxis.set_minor_locator(yminorLocator)
ax3.yaxis.set_tick_params(which = 'minor', direction = 'in', width = 1, length = 3, pad = 5, left = True, right = True)

labels = [item.get_text() for item in ax3.get_yticklabels()]
empty_string_labels = ['']*len(labels)
ax3.set_yticklabels(empty_string_labels)

ttl = ax3.set_title(tit2, fontname='times', fontsize = fsize2, fontweight = fweight, loc='center')
ttl.set_position(tit_pos)

plt.text(12.5, 0.00600, disp41, horizontalalignment='left', verticalalignment='top', fontname='times', fontsize = 14, transform = ax3.transData, color = 'g')
plt.text(12.5, 0.00363, disp42, horizontalalignment='left', verticalalignment='top', fontname='times', fontsize = 14, transform = ax3.transData, color = 'b')
plt.text(12.5, 0.00220, disp43, horizontalalignment='left', verticalalignment='top', fontname='times', fontsize = 14, transform = ax3.transData, color = 'r')
plt.text(12.5, 0.00133, disp44, horizontalalignment='left', verticalalignment='top', fontname='times', fontsize = 14, transform = ax3.transData, color = 'k')

ax3.set_xlim([xmin_NHI, xmax_NHI])
ax3.set_ylim([10.0**ymin, 10.0**ymax])

#######################

savefig(fname)
print(fname)
plt.close()

########################
