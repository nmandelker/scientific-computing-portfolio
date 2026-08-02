from gadget import *
from gadget_subfind import *
import const as c
import matplotlib as mpl
import os, sys

SimLabel = 'FilB-Zoom2'
SnapNums = np.arange(78,126)

in_filename = "halo-positions-%s.dat" % (SimLabel)

SimDir   = "/home/mandelnr/basement/%s/output" % (SimLabel)
SnapBase = 'snap_'
FigDir   = "../plots/%s" % (SimLabel)
if not os.path.exists(FigDir):
        os.mkdir(FigDir)
FigDir   = "%s/Temp_gas/" % (FigDir)
if not os.path.exists(FigDir):
        os.mkdir(FigDir)

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

gravity = c.G                        #cm^3/(g s^2)
msolar = c.msol                      #gram
mproton = c.amu                      #gram
kB = c.KB                            #erg/K
kpcincm = c.parsec*1e3
Xh = 0.76
gamma_gas = 5.0/3.0

boxx = 3600.
boxy = 2000.
boxz = 400. 
pixres = 2. 
imgcent = [0,0,0]
xsize = 18.0
ysize = 5.2857
fsize = 18
fsize2 = 22
fweight = 150

for SnapNum in SnapNums:
	print(SnapNum)
	s  = gadget_readsnap( SnapNum, SimDir, SnapBase, loadonlytype=[0], loadonly=['pos', 'vol', 'rho', 'mass', 'ne', 'u'], lazy_load=True )
	sf = load_subfind( SnapNum, dir=SimDir+'/' )
	redshift = s.redshift
	BoxHalf = 0.5 * s.boxsize
	aexp = s.time
	h0 = s.hubbleparam
	BoxHalf *= aexp/h0
	BoxSize = 2.0*BoxHalf
	ugas = s.u
	negas = s.ne
	MeanWeight= 4.0/(3*Xh+1+4*Xh*negas) * mproton
	Temp = ugas * MeanWeight/kB * (gamma_gas-1) * 1.989e53 / (1e10 * msolar) # * UnitEnergy_in_cgs/ UnitMass_in_g;		Now we have temperature in Kelvin

	b2, = np.where(RevSnap2-SnapNum==0)
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

	# determine mid-point between two largest halos (Halos with most particles)
	xc = x_first - x_second
	yc = y_first - y_second
	zc = z_first - z_second
	print('xc,yc,zc=',xc, yc, zc)
	
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
	gM = sf.data['fmas']
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
	xx = xp*nx + yp*ny + zp*nz # x coordinate of point in disk frame
	yy = xp*mx + yp*my + zp*mz # y coordinate of point in disk frame
	zz = xp*Jx + yp*Jy + zp*Jz # z coordinate of point in disk frame

	gxx = gxp*nx + gyp*ny + gzp*nz # x coordinate of point in disk frame
	gyy = gxp*mx + gyp*my + gzp*mz # y coordinate of point in disk frame
	gzz = gxp*Jx + gyp*Jy + gzp*Jz # z coordinate of point in disk frame

	s.data['pos'][:,0] = xx
	s.data['pos'][:,1] = yy
	s.data['pos'][:,2] = zz
	s.data['MT']   = s.rho * Temp
	s.data['pos'] *= h0/aexp # comiving h^-1 kpc
	gxx *= h0/aexp # comiving h^-1 kpc
	gyy *= h0/aexp # comiving h^-1 kpc
	gzz *= h0/aexp # comiving h^-1 kpc
	gR  *= h0/aexp # comiving h^-1 kpc

	proj_xz_Q = s.get_Aslice("MT", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz, axes = [2,0], proj = True, numthreads=32)
	proj_xz_W = s.get_Aslice("rho", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz, axes = [2,0], proj = True, numthreads=32)
	mapp_xz = proj_xz_Q['grid'] / proj_xz_W['grid']

	proj_yz_Q = s.get_Aslice("MT", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz, axes = [2,1], proj = True, numthreads=32)
	proj_yz_W = s.get_Aslice("rho", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz, axes = [2,1], proj = True, numthreads=32)
	mapp_yz = proj_yz_Q['grid'] / proj_yz_W['grid']

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	fig = plt.figure(figsize = (xsize,ysize))
	ax1 = plt.axes([0.07, 0.18, 0.37, 0.70])
	plt.pcolormesh( proj_yz_Q['x'], proj_yz_Q['y'], np.transpose(mapp_yz), norm = matplotlib.colors.LogNorm(vmin=1.0e4, vmax=1.0e7), cmap = 'RdBu_r', rasterized = True )
	xlabel(r"$z\:[h^{-1}{\rm kpc\: comoving}]$", fontname='times', fontsize = fsize2, fontweight = fweight)
	ylabel(r"$y'\:[h^{-1}{\rm kpc\: comoving}]$", fontname='times', fontsize = fsize2, fontweight = fweight)
	plt.xticks([-1500, -1000, -500, 0, 500, 1000, 1500], fontname='arial', fontsize = fsize, fontweight = fweight) #[' ', u'\u2212''1000', ' ', 0, ' ', '1000', ' ']
	plt.yticks([-1000, -500, 0, 500, 1000], fontname='arial', fontsize = fsize, fontweight = fweight)
	ax1.xaxis.set_tick_params(which = 'major', direction = 'out', width = 1, length = 6, pad = 5)
	ax1.yaxis.set_tick_params(which = 'major', direction = 'out', width = 1, length = 6, pad = 5)
	ax1.xaxis.set_label_coords(0.5, -0.12)
	ax1.yaxis.set_label_coords(-0.12, 0.5)

	ttl = ax1.set_title(r"${\rm " + SimLabel + r",\:\:\: projected\:density-weighted\:T_{\rm gas},\:\:\:z=" + '{0:.2f}'.format(redshift) + "}$", fontname='times', fontsize = fsize2, fontweight = fweight, loc='center')
	ttl.set_position([0.5, 1.02])

	ax2 = plt.axes([0.52, 0.18, 0.37, 0.70])
	plt.pcolormesh( proj_xz_Q['x'], proj_xz_Q['y'], np.transpose(mapp_xz), norm = matplotlib.colors.LogNorm(vmin=1.0e4, vmax=1.0e7), cmap = 'RdBu_r', rasterized = True )
	xlabel(r"$z\:[h^{-1}{\rm kpc\: comoving}]$", fontname='times', fontsize = fsize2, fontweight = fweight)
	ylabel(r"$x'\:[h^{-1}{\rm kpc\: comoving}]$", fontname='times', fontsize = fsize2, fontweight = fweight)
	plt.xticks([-1500, -1000, -500, 0, 500, 1000, 1500], fontname='arial', fontsize = fsize, fontweight = fweight) #[' ', u'\u2212''1000', ' ', 0, ' ', '1000', ' ']
	plt.yticks([-1000, -500, 0, 500, 1000], fontname='arial', fontsize = fsize, fontweight = fweight)
	ax2.xaxis.set_tick_params(which = 'major', direction = 'out', width = 1, length = 6, pad = 5)
	ax2.yaxis.set_tick_params(which = 'major', direction = 'out', width = 1, length = 6, pad = 5)
	ax2.xaxis.set_label_coords(0.5, -0.12)
	ax2.yaxis.set_label_coords(-0.12, 0.5)

	axcbar = fig.add_axes([0.90, 0.19, 0.025, 0.6956])
	cbar = matplotlib.colorbar.ColorbarBase(ax = axcbar, norm = matplotlib.colors.LogNorm(vmin=1.0e4, vmax=1.0e7), cmap = 'RdBu_r', orientation='vertical')
	cbar.set_label(label=r"$T{\rm [K]}$", fontname='times', fontsize = fsize2, fontweight = fweight, rotation=270, labelpad=30)
	cbar.ax.tick_params(direction = 'out', labelsize = fsize, width = 1, length = 6, pad = 5)

        ax1.set_xlim([-boxx/2.0,boxx/2.0])
        ax1.set_ylim([-boxy/2.0,boxy/2.0])
        ax2.set_xlim([-boxx/2.0,boxx/2.0])
        ax2.set_ylim([-boxy/2.0,boxy/2.0])
	fname = FigDir+"snap_%03d.png" % (SnapNum)
	savefig(fname)
	print(fname)

	t=np.linspace(0.0,2.0*pi,100)
	for i in range(Ng):
		if gR[i]>0. and gM[i]>0.1 and gzz[i]+gR[i]>=-boxx/2.0 and gzz[i]-gR[i]<=boxx/2.0:
			if gyy[i]+gR[i]>=-boxy/2.0 and gyy[i]-gR[i]<=boxy/2.0 and gxx[i]>=-boxz/2.0 and gxx[i]<=boxz/2.0:
				ax1.plot(gzz[i]+gR[i]*np.cos(t), gyy[i]+gR[i]*np.sin(t), color='k', linestyle='-', linewidth=1.0)
			if gxx[i]+gR[i]>=-boxy/2.0 and gxx[i]-gR[i]<=boxy/2.0 and gyy[i]>=-boxz/2.0 and gyy[i]<=boxz/2.0:
				ax2.plot(gzz[i]+gR[i]*np.cos(t), gxx[i]+gR[i]*np.sin(t), color='k', linestyle='-', linewidth=1.0)

	ax1.set_xlim([-boxx/2.0,boxx/2.0])
	ax1.set_ylim([-boxy/2.0,boxy/2.0])
	ax2.set_xlim([-boxx/2.0,boxx/2.0])
	ax2.set_ylim([-boxy/2.0,boxy/2.0])
	fname = FigDir+"with_halos_snap_%03d.png" % (SnapNum)
	savefig(fname)
	print(fname)

	plt.close()


