from gadget import *
from gadget_subfind import *
import const as c
import matplotlib.pyplot as plt
import matplotlib as mpl
import os, sys
#from numba import jit
#os.environ['KMP_DUPLICATE_LIB_OK']='True'

ZF = 4
SnapNum = 99
if ZF==4:
    SnapNum2 = SnapNum - 78
else:
    SnapNum2 = SnapNum
print('ZF, SnapNum, SnapNum2 = ', ZF, SnapNum, SnapNum2)

exclude = 0         # 0: Include all cells;	1: Exclude all cells in FoF groups with M>min_mass; 2: Exclude all cells in Subfind subhalos with M>min_mass; 5: Exclude all cells in FoF groups with N_{dm,part}>min_Npart
min_mass  = 1.0e7	# Minimum halo mas for exclusion rules
min_Npart = 32		# Minimum number of DM particles in FoF groups that are excluded if exclude==5

h_over_Delta = 0.0 # 1.0  # If this is >0, this is SPH smoothing kernel (radius) in units of the effective cell size (diameter). IF THIS IS 0, USE NEAREST NEIGHBOUR GAS CELL INSTEAD OF SPH SMOOTHING

L_slice = 15.0 # kpc - half length of the filament slices

# Bins

#Delta_z      = 0.30              # kpc
#Delta_r_norm = 1.0 / 180.0       # in units of r_shock - this is between ~(0.17-0.45)kpc for the different slices.
#Delta_p      = 2.0*np.pi/180.0   # 2.0 degrees ~ 0.035 radians

#Delta_z      = 0.40              # kpc
#Delta_r_norm = 1.0 / 120.0       # in units of r_shock - this is between ~(0.25-0.75)kpc for the different slices.
#Delta_p      = 2.0*np.pi/120.0   # 3.0 degrees ~ 0.05 radians

#Delta_z      = 0.50              # kpc
#Delta_r_norm = 1.0 / 100.0       # in units of r_shock - this is between ~(0.3-0.8)kpc for the different slices.
#Delta_p      = 2.0*np.pi / 90.0  # 4.0 degrees ~ 0.07 radians

Delta_z      = 0.70              # kpc
Delta_r_norm = 1.0 / 80.0        # in units of r_shock - this is between ~(0.35-1.0)kpc for the different slices.
Delta_p      = 2.0*np.pi/70.0    # 5.1 degrees ~ 0.09 radians

Zmin      = -L_slice - 4.0*Delta_z    # kpc - added a buffer so we can do high order derivatives at the z boundaries
Zmax      =  L_slice + 4.0*Delta_z    # kpc - added a buffer so we can do high order derivatives at the z boundaries
Rmin_norm =  0.02     # r_shock
Rmax_norm =  1.1      # r_shock
Pmin      =  0.0      # radians
Pmax      =  2.0*np.pi - Delta_p

zbin_num   = int(np.floor((Zmax-Zmin)/Delta_z) + 1)
rbin_num   = int(np.floor((Rmax_norm-Rmin_norm)/Delta_r_norm) + 1)
pbin_num   = int(np.floor((Pmax-Pmin)/Delta_p) + 1)
Nbin       = rbin_num * pbin_num * zbin_num

z_bin      = Zmin      + np.arange(0,zbin_num)*Delta_z
r_bin_norm = Rmin_norm + np.arange(0,rbin_num)*Delta_r_norm
p_bin      = Pmin      + np.arange(0,pbin_num)*Delta_p

print('rbin_num, pbin_num, zbin_num = ', rbin_num, pbin_num, zbin_num)
print('Nbin = ', rbin_num * pbin_num * zbin_num, Nbin)

# The centres we've been using; slices_1: 1 - 6; slices_2: 7 - 10
## In these scripts I do slices 1, 2, 3, 4, and 5, 6 (vertical). Since they require extra rotations.
'''
# ZF3
slices_1 = [[-0.2593815263129492, -156.77116040128288, -278.5769794500865], 
[19.29960307468462, -180.02177342681858, 18.783482611124054], 
[11.942188924548926, -172.22860052014127, 55], [18.958866455674514, -159.23324755380037, 135], 
[36.1618831230089, -16.19509203477297, 330.38220545849816], [56.19692148069989, -64.94119207620751, 357.6376986767729]]

slices_2 = [[29.97807990484243, -29.468658288254147, 53.61269984845511], [31, -29, 77], [31, -29, 113], 
[29.40874213348795, -29.541646560122633, 204.96891537756946]]

r_shock_list = [array([76.74496644]), array([64.7852349]), array([44.18791946]), array([34.88590604]), 
array([40.86577181]), array([49.5033557]), array([35.55033557]), array([30.23489933]), array([31.56375839]), array([30.89932886])]
'''
# ZF4
slices_1 = np.array([np.array([  -0.42191654, -156.95545196, -278.57697945]),
       np.array([  18.77496127, -179.99887174,   18.78348261]),
       np.array([  11.69567357, -172.66565623,   55.        ]),
       np.array([  18.8371828 , -159.61706258,  135.        ]),
       np.array([  36.16379242,  -16.19509203,  330.4726014 ]),
       np.array([  55.09858281,  -64.94119208,  356.43792306])])
slices_2 = np.array([np.array([ 30.06170486, -29.4384368 ,  53.61269985]),
       np.array([ 30.49680243, -29.12472163,  77.        ]),
       np.array([ 30.99852922, -28.69353305, 113.        ]),
       np.array([ 30.02304343, -28.92417905, 204.96891538])])

r_shock_list = np.array([78.73825503355705, 61.463087248322154, 50.167785234899334, 41.530201342281885, 38.872483221476514, 
48.83892617449665, 40.2013422818792, 30.23489932885906, 29.57046979865772, 32.89261744966443])

SimLabel = 'FilB-Zoom%01d' % (ZF)
SimDir   = "/sci/labs/nmandelker/nmandelker/IPM_sims/%s/output" % (SimLabel)
SnapBase = 'snap_'
FigDir   = "../plots/FilB-Zoom4"
if not os.path.exists(FigDir):
    os.mkdir(FigDir)
FigDir   = "%s/Sam_SPH_method" % (FigDir)
if not os.path.exists(FigDir):
    os.mkdir(FigDir)
FigDir   = FigDir+"/Nr_"+'{:.0f}'.format(rbin_num)+"_Np_"+'{:.0f}'.format(pbin_num)+"_Nz_"+'{:.0f}'.format(zbin_num)+"_h_"+'{0:.1f}'.format(h_over_Delta)+"_exclude_"+'{0:.1f}'.format(exclude)
if not os.path.exists(FigDir):
    os.mkdir(FigDir)

in_filename = "halo-positions-%s.dat" % (SimLabel)

posdata   = np.loadtxt(in_filename)
RevSnap2  = posdata[0,:]
aexp_vec2 = posdata[1,:]
x1_vec    = posdata[2,:]
y1_vec    = posdata[3,:]
z1_vec    = posdata[4,:]
x2_vec    = posdata[5,:]
y2_vec    = posdata[6,:]
z2_vec    = posdata[7,:]
ind1_vec  = posdata[8,:]
ind2_vec  = posdata[9,:]

gravity    = c.G                        #cm^3/(g s^2)
msolar     = c.msol                     #gram
mproton    = c.amu                      #gram
kB         = c.KB                       #erg/K
kpcincm    = c.parsec*1e3
Xh         = 0.76
gamma_gas  = 5.0/3.0
Omegab     = 0.0486
HUBBLE     = 3.2407789e-18              # in h/sec
sec_in_Gyr = 60.0*60.0*24.0*365.25*1.0e9

UnitLength_in_cm         = 3.085678e21  # kpc
UnitMass_in_g            = 1.989e43     # 1.0e10 Msun
UnitVelocity_in_cm_per_s = 1.0e5        # 1 km/sec
UnitTime_in_s            = UnitLength_in_cm / UnitVelocity_in_cm_per_s               
UnitDensity_in_cgs       = UnitMass_in_g    / UnitLength_in_cm**3                    # in g/cm^3
UnitPressure_in_cgs      = UnitMass_in_g    / UnitLength_in_cm    / UnitTime_in_s**2 # in dyne/cm^2
UnitEnergy_in_cgs        = UnitMass_in_g    * UnitLength_in_cm**2 / UnitTime_in_s**2 # in ergs
UnitNumberDensity_in_cgs = UnitDensity_in_cgs / mproton			                     # in cm^{-3}

s  = gadget_readsnap( SnapNum2, SimDir, SnapBase, loadonlytype=[0], loadonly=['pos', 'vel', 'vol', 'rho', 'mass', 'u', 'ne', 'gmet', 'pot'], lazy_load=True )
sd = gadget_readsnap( SnapNum2, SimDir, SnapBase, loadonlytype=[1], loadonly=['pos', 'pot'], lazy_load=True ) # Dark Matter
ss = gadget_readsnap( SnapNum2, SimDir, SnapBase, loadonlytype=[4], loadonly=['pos', 'mass', 'pot'], lazy_load=True ) # Stars
sf = load_subfind( SnapNum2, dir=SimDir+'/' )

redshift = s.redshift
BoxHalf = 0.5 * s.boxsize   #ckpc/h
aexp = s.time
h0 = s.hubbleparam
BoxHalf *= aexp/h0          #pkpc
BoxSize = 2.0*BoxHalf       #pkpc

### Gas property arrays ###
ugas       = s.u * UnitEnergy_in_cgs / UnitMass_in_g			# internal energy per unit mass in erg/gr
MeanWeight = 4.0 / ( 3.0*Xh + 1.0 + 4.0*Xh*s.ne ) * mproton     # gram
Temp       = ugas * MeanWeight/kB * (gamma_gas-1.0)				# temperature in Kelvin
mass       = s.mass * UnitMass_in_g                             # gas cell mass in grams
rho        = s.rho * UnitDensity_in_cgs                         # gas cell total density in gr/cm^3
pot        = s.pot * UnitEnergy_in_cgs / UnitMass_in_g          # gas cell gravitational potential in (cm/sec)^2
pressure   = kB * Temp * rho / MeanWeight                       # gas thermal pressure in erg/cm^3
volume     = 1.0 * s.vol[:] # * UnitLength_in_cm**3             # gas cell volume in kpc^3

stmass  = ss.mass * UnitMass_in_g                               # stellar particle mass in grams
ZF_real = np.maximum(1.0*ZF, 0.5)
dmass   = 8.2e4 * 1.989e33 * (4.0/ZF_real)**3                   # dark matter particle mass in grams

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

# compute coordinates relative to the center of the filament (coordinates of gas cells)
xp = s.data['pos'][:,0] - xc
yp = s.data['pos'][:,1] - yc
zp = s.data['pos'][:,2] - zc

#Same thing, for the dark matter particles
dxp = sd.data['pos'][:,0] - xc
dyp = sd.data['pos'][:,1] - yc
dzp = sd.data['pos'][:,2] - zc

#Same thing, for stars; st stands for stars
stxp = ss.data['pos'][:,0] - xc
styp = ss.data['pos'][:,1] - yc
stzp = ss.data['pos'][:,2] - zc

vxp = 1.0 * s.data['vel'][:,0]
vyp = 1.0 * s.data['vel'][:,1]
vzp = 1.0 * s.data['vel'][:,2]

# Friends of Friends groups
gxp = sf.data['fpos'][:,0] - xc
gyp = sf.data['fpos'][:,1] - yc
gzp = sf.data['fpos'][:,2] - zc
gR  = sf.data['frt2']
gM  = sf.data['fmas']
Ng  = sf.totngroups

# subhalo groups
sxp = sf.data['spos'][:,0] - xc
syp = sf.data['spos'][:,1] - yc
szp = sf.data['spos'][:,2] - zc

vx_halos = sf.data['svel'][:,0]
vy_halos = sf.data['svel'][:,1]
vz_halos = sf.data['svel'][:,2]

# in case we have a periodic wrap, let's take care of it (should actually not be necessary here)
# Gas
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

# DM
b1, = np.where(dxp<-BoxHalf)
if b1.size>0:
    dxp[b1] += BoxSize
b1, = np.where(dxp>BoxHalf)
if b1.size>0:
    dxp[b1] -= BoxSize
b1, = np.where(dyp<-BoxHalf)
if b1.size>0:
    dyp[b1] += BoxSize
b1, = np.where(dyp>BoxHalf)
if b1.size>0:
    dyp[b1] -= BoxSize
b1, = np.where(dzp<-BoxHalf)
if b1.size>0:
    dzp[b1] += BoxSize
b1, = np.where(dzp>BoxHalf)
if b1.size>0:
    dzp[b1] -= BoxSize

# Stars
b1, = np.where(stxp<-BoxHalf)
if b1.size>0:
    stxp[b1] += BoxSize
b1, = np.where(stxp>BoxHalf)
if b1.size>0:
    stxp[b1] -= BoxSize
b1, = np.where(styp<-BoxHalf)
if b1.size>0:
    styp[b1] += BoxSize
b1, = np.where(styp>BoxHalf)
if b1.size>0:
    styp[b1] -= BoxSize
b1, = np.where(stzp<-BoxHalf)
if b1.size>0:
    stzp[b1] += BoxSize
b1, = np.where(stzp>BoxHalf)
if b1.size>0:
    stzp[b1] -= BoxSize

# Groups
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

# Subhalos
b1, = np.where(sxp<-BoxHalf)
if b1.size>0:
    sxp[b1] += BoxSize
b1, = np.where(sxp>BoxHalf)
if b1.size>0:
    sxp[b1] -= BoxSize
b1, = np.where(syp<-BoxHalf)
if b1.size>0:
    syp[b1] += BoxSize
b1, = np.where(syp>BoxHalf)
if b1.size>0:
    syp[b1] -= BoxSize
b1, = np.where(szp<-BoxHalf)
if b1.size>0:
    szp[b1] += BoxSize
b1, = np.where(szp>BoxHalf)
if b1.size>0:
    szp[b1] -= BoxSize

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
# No special meaning for 45 degrees, It's just a coincidence
nxx = nx*np.cos(np.pi/4.) + mx*np.sin(np.pi/4.)
nyy = ny*np.cos(np.pi/4.) + my*np.sin(np.pi/4.)
nzz = nz*np.cos(np.pi/4.) + mz*np.sin(np.pi/4.)

mxx = mx*np.cos(np.pi/4.) - nx*np.sin(np.pi/4.)
myy = my*np.cos(np.pi/4.) - ny*np.sin(np.pi/4.)
mzz = mz*np.cos(np.pi/4.) - nz*np.sin(np.pi/4.)

nx = 1.0 * nxx
ny = 1.0 * nyy
nz = 1.0 * nzz
mx = 1.0 * mxx
my = 1.0 * myy
mz = 1.0 * mzz

# Sanity checks
print('coordinate dot products, Slices 7-10', Jx * mx + Jy * my + Jz * mz,    Jx * nx + Jy * ny + Jz * nz,    mx * nx + my * ny + mz * nz)
print('coordinate amplitudes, Slices 7-10',   np.sqrt(Jx**2 + Jy**2 + Jz**2), np.sqrt(nx**2 + ny**2 + nz**2), np.sqrt(mx**2 + my**2 + mz**2))

# Now rotate all the vector quantities
xx7   = xp*nx   + yp*ny   + zp*nz	# x coordinate of point in disk frame, gas
yy7   = xp*mx   + yp*my   + zp*mz	# y coordinate of point in disk frame
zz7   = xp*Jx   + yp*Jy   + zp*Jz	# z coordinate of point in disk frame

dxx7  = dxp*nx  + dyp*ny  + dzp*nz	# x coordinate of point in disk frame, dm
dyy7  = dxp*mx  + dyp*my  + dzp*mz	# y coordinate of point in disk frame
dzz7  = dxp*Jx  + dyp*Jy  + dzp*Jz	# z coordinate of point in disk frame

stxx7 = stxp*nx + styp*ny + stzp*nz	# x coordinate of point in disk frame, stars
styy7 = stxp*mx + styp*my + stzp*mz	# y coordinate of point in disk frame
stzz7 = stxp*Jx + styp*Jy + stzp*Jz	# z coordinate of point in disk frame

gxx7  = gxp*nx  + gyp*ny  + gzp*nz	# x coordinate of point in disk frame, groups
gyy7  = gxp*mx  + gyp*my  + gzp*mz	# y coordinate of point in disk frame
gzz7  = gxp*Jx  + gyp*Jy  + gzp*Jz	# z coordinate of point in disk frame

sxx7  = sxp*nx  + syp*ny  + szp*nz	# x coordinate of point in disk frame, subhalos
syy7  = sxp*mx  + syp*my  + szp*mz	# y coordinate of point in disk frame
szz7  = sxp*Jx  + syp*Jy  + szp*Jz	# z coordinate of point in disk frame

vxx7  = vxp*nx  + vyp*ny  + vzp*nz	# x velocity coordinate of point in rotated frame
vyy7  = vxp*mx  + vyp*my  + vzp*mz	# y velocity coordinate of point in rotated frame
vzz7  = vxp*Jx  + vyp*Jy  + vzp*Jz	# z velocity coordinate of point in rotated frame

vxx7_halos = vx_halos*nx + vy_halos*ny + vz_halos*nz
vyy7_halos = vx_halos*mx + vy_halos*my + vz_halos*mz
vzz7_halos = vx_halos*Jx + vy_halos*Jy + vz_halos*Jz

# This is another rotation so that the longest filament is aligned with the z axis
# No need for slices 7, 8, 9, 10
mxx = mx*np.cos(np.pi/9.) + Jx*np.sin(np.pi/9.) # pi/9 = 20/180*pi = 20 degrees
myy = my*np.cos(np.pi/9.) + Jy*np.sin(np.pi/9.)
mzz = mz*np.cos(np.pi/9.) + Jz*np.sin(np.pi/9.)

Jxx = Jx*np.cos(np.pi/9.) - mx*np.sin(np.pi/9.)
Jyy = Jy*np.cos(np.pi/9.) - my*np.sin(np.pi/9.)
Jzz = Jz*np.cos(np.pi/9.) - mz*np.sin(np.pi/9.)

mx = 1.0 * mxx
my = 1.0 * myy
mz = 1.0 * mzz
Jx = 1.0 * Jxx
Jy = 1.0 * Jyy
Jz = 1.0 * Jzz

# Sanity checks
print('coordinate dot products, Slices 1-6', Jx * mx + Jy * my + Jz * mz,    Jx * nx + Jy * ny + Jz * nz,    mx * nx + my * ny + mz * nz)
print('coordinate amplitudes, Slices 1-6',   np.sqrt(Jx**2 + Jy**2 + Jz**2), np.sqrt(nx**2 + ny**2 + nz**2), np.sqrt(mx**2 + my**2 + mz**2))

# Now rotate all the vector quantities
xx1   = xp*nx   + yp*ny   + zp*nz	# x coordinate of point in disk frame, gas
yy1   = xp*mx   + yp*my   + zp*mz	# y coordinate of point in disk frame
zz1   = xp*Jx   + yp*Jy   + zp*Jz	# z coordinate of point in disk frame

dxx1  = dxp*nx  + dyp*ny  + dzp*nz	# x coordinate of point in disk frame, dm
dyy1  = dxp*mx  + dyp*my  + dzp*mz	# y coordinate of point in disk frame
dzz1  = dxp*Jx  + dyp*Jy  + dzp*Jz	# z coordinate of point in disk frame

stxx1 = stxp*nx + styp*ny + stzp*nz	# x coordinate of point in disk frame, stars
styy1 = stxp*mx + styp*my + stzp*mz	# y coordinate of point in disk frame
stzz1 = stxp*Jx + styp*Jy + stzp*Jz	# z coordinate of point in disk frame

gxx1  = gxp*nx  + gyp*ny  + gzp*nz	# x coordinate of point in disk frame, groups
gyy1  = gxp*mx  + gyp*my  + gzp*mz	# y coordinate of point in disk frame
gzz1  = gxp*Jx  + gyp*Jy  + gzp*Jz	# z coordinate of point in disk frame

sxx1  = sxp*nx  + syp*ny  + szp*nz	# x coordinate of point in disk frame, subhalos
syy1  = sxp*mx  + syp*my  + szp*mz	# y coordinate of point in disk frame
szz1  = sxp*Jx  + syp*Jy  + szp*Jz	# z coordinate of point in disk frame

vxx1  = vxp*nx  + vyp*ny  + vzp*nz	# x velocity coordinate of point in rotated frame
vyy1  = vxp*mx  + vyp*my  + vzp*mz	# y velocity coordinate of point in rotated frame
vzz1  = vxp*Jx  + vyp*Jy  + vzp*Jz	# z velocity coordinate of point in rotated frame

vxx1_halos = vx_halos*nx + vy_halos*ny + vz_halos*nz
vyy1_halos = vx_halos*mx + vy_halos*my + vz_halos*mz
vzz1_halos = vx_halos*Jx + vy_halos*Jy + vz_halos*Jz


################################
### Remove halos ###

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
			mass[first_cell:last_cell] = 0.0
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
			mass[first_cell:last_cell] = 0.0
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
					mass[first_cell:last_cell] = 0.0
				first_cell += Ncell
		first_cell = np.int(np.sum(sf.data['flty'][0:g_ind,0]))

'''
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
# Test the centres and shock radii by making a large-scale map #

s.data['pos'][:,0] = xx
s.data['pos'][:,1] = yy
s.data['pos'][:,2] = zz
s.data['pos'] *= h0/aexp # comiving h^-1 kpc
gxx *= h0/aexp
gyy *= h0/aexp
gzz *= h0/aexp
gR  *= h0/aexp
slices_1 *= h0/aexp
slices_2 *= h0/aexp
r_shock_list *= h0/aexp
slice_length = 30.0 * h0/aexp

s.data['vxp'] = s.rho * vxx
s.data['vyp'] = s.rho * vyy
s.data['vzp'] = s.rho * vzz
s.data['vxp2'] = s.rho * vxx**2
s.data['vyp2'] = s.rho * vyy**2
s.data['vzp2'] = s.rho * vzz**2

boxx = 3600.
boxy = 2000.
boxz = 800.
pixres = 2. 
imgcent = [0,0,0]
xsize = 9.2
ysize = 5.3
fsize = 18
fsize2 = 22
fweight = 150

ax1pos = [0.127, 0.18, 0.725, 0.70]
caxpos = [0.87, 0.18, 0.025, 0.6956]

xtick    = [-1500, -1000, -500, 0, 500, 1000, 1500]
ytick    = [-1000, -500, 0, 500, 1000]
dx_minor = 100
dy_minor = 100
xlaby    = -0.12
ylabx    = -0.12
tit_pos  = [0.2, 1.02]

xstr   = r"$z\:[h^{-1}{\rm kpc\: comoving}]$"
ystr   = r"$y'\:[h^{-1}{\rm kpc\: comoving}]$"
titstr = r"${\rm " + SimLabel + r",\:\:\: z=" + '{0:.2f}'.format(redshift) + "}$"

s.data['rho'] = s.rho * s.gmet[:,0] * UnitNumberDensity_in_cgs	# nH  cm^{-3}
proj_yz_Q = s.get_Aslice("rho", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz, axes = [2,1], proj = True, numthreads=32)
mapp_yz   = np.transpose( proj_yz_Q['grid'] ) * (pixres*aexp/h0)*UnitLength_in_cm	# NH cm^{-2}

xpl = proj_yz_Q['x']
ypl = proj_yz_Q['y']

ma = 1.0e22
mi = 3.0e18

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig = plt.figure(figsize = (xsize,ysize))
ax1 = plt.axes(ax1pos)
plt.pcolormesh( xpl, ypl, mapp_yz, norm = matplotlib.colors.LogNorm(vmin=mi, vmax=ma), cmap = 'nipy_spectral', rasterized = True )
xlabel(xstr, fontname='times', fontsize = fsize2, fontweight = fweight)
ylabel(ystr, fontname='times', fontsize = fsize2, fontweight = fweight)

plt.xticks(xtick, fontname='arial', fontsize = fsize, fontweight = fweight)
plt.yticks(ytick, fontname='arial', fontsize = fsize, fontweight = fweight)

ax1.xaxis.set_tick_params(which = 'major', direction = 'out', width = 1, length = 6, pad = 5, bottom = True, top = True)
ax1.yaxis.set_tick_params(which = 'major', direction = 'out', width = 1, length = 6, pad = 5, bottom = True, top = True)

xminorLocator = MultipleLocator(dx_minor)
ax1.xaxis.set_minor_locator(xminorLocator)
ax1.xaxis.set_tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5, bottom = True, top = True)

yminorLocator = MultipleLocator(dy_minor)
ax1.yaxis.set_minor_locator(yminorLocator)
ax1.yaxis.set_tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5, left = True, right = True)

ax1.xaxis.set_label_coords(0.5, xlaby)
ax1.yaxis.set_label_coords(ylabx, 0.5)

ttl = ax1.set_title(titstr, fontname='times', fontsize = fsize2, fontweight = fweight, loc='left')
ttl.set_position(tit_pos)

axcbar = fig.add_axes(caxpos)
cbar = matplotlib.colorbar.ColorbarBase(ax = axcbar, norm = matplotlib.colors.LogNorm(vmin=mi, vmax=ma), cmap = 'nipy_spectral', orientation='vertical', ticks = LogLocator(subs=range(10)))
cbar.set_label(label=r"$N_{\rm H}\:[{\rm cm^{-2}}]$", fontname='times', fontsize = fsize2, fontweight = fweight, rotation=270, labelpad=19)
cbar.ax.tick_params(which = 'major', direction = 'out', labelsize = fsize, width = 1, length = 6, pad = 5)
cbar.ax.tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5)

t=np.linspace(0.0,2.0*pi,100)
gr_cnt  = -1
for i in range(Ng):
	if gR[i]>0. and gM[i]>1.0 and gzz[i]+gR[i]>=-boxx/2.0 and gzz[i]-gR[i]<=boxx/2.0:
		if gyy[i]+gR[i]>=-boxy/2.0 and gyy[i]-gR[i]<=boxy/2.0 and gxx[i]>=-boxz/2.0 and gxx[i]<=boxz/2.0:
			ax1.plot(gzz[i]+gR[i]*np.cos(t), gyy[i]+gR[i]*np.sin(t), color='k', linestyle='-', linewidth=1.0)

for reg_ind in np.arange(0,4):
    xc1 = slices_1[reg_ind][0]
    yc1 = slices_1[reg_ind][1]
    zc1 = slices_1[reg_ind][2]
    rs1 = r_shock_list[reg_ind]

    tx = np.linspace(zc1-0.5*slice_length,zc1+0.5*slice_length,100)
    ty = (yc1-rs1) * np.ones(100)
    ax1.plot(tx, ty, color='k', linestyle='-', linewidth=1.0)
    ty = (yc1+rs1) * np.ones(100)
    ax1.plot(tx, ty, color='k', linestyle='-', linewidth=1.0)
    ty = np.linspace(yc1-rs1,yc1+rs1,100)
    tx = (zc1-0.5*slice_length) * np.ones(100)
    ax1.plot(tx, ty, color='k', linestyle='-', linewidth=1.0)
    tx = (zc1+0.5*slice_length) * np.ones(100)
    ax1.plot(tx, ty, color='k', linestyle='-', linewidth=1.0)
    plt.text(zc1, yc1-0.4*rs1, r"${\rm S%01d}$" % (reg_ind+1), horizontalalignment='center', fontname='times', fontsize = fsize, fontweight = fweight, transform = ax1.transData)

for reg_ind in np.arange(4,6):
    xc1 = slices_1[reg_ind][0]
    yc1 = slices_1[reg_ind][1]
    zc1 = slices_1[reg_ind][2]
    rs1 = r_shock_list[reg_ind]

    tx = np.linspace(yc1-0.5*slice_length,yc1+0.5*slice_length,100)
    ty = (zc1-rs1) * np.ones(100)
    ax1.plot(ty, tx, color='k', linestyle='-', linewidth=1.0)
    ty = (zc1+rs1) * np.ones(100)
    ax1.plot(ty, tx, color='k', linestyle='-', linewidth=1.0)
    ty = np.linspace(zc1-rs1,zc1+rs1,100)
    tx = (yc1-0.5*slice_length) * np.ones(100)
    ax1.plot(ty, tx, color='k', linestyle='-', linewidth=1.0)
    tx = (yc1+0.5*slice_length) * np.ones(100)
    ax1.plot(ty, tx, color='k', linestyle='-', linewidth=1.0)
    plt.text(zc1-0.4*rs1, yc1, r"${\rm S%01d}$" % (reg_ind+1), horizontalalignment='center', fontname='times', fontsize = fsize, fontweight = fweight, transform = ax1.transData)

for reg_ind in np.arange(0,4):
    xc1 = slices_2[reg_ind][0]
    yc1 = slices_2[reg_ind][1]
    zc1 = slices_2[reg_ind][2]
    rs1 = r_shock_list[reg_ind+6]

    tx = np.linspace(zc1-0.5*slice_length,zc1+0.5*slice_length,100)
    ty = (yc1-rs1) * np.ones(100)
    ax1.plot(tx, ty, color='k', linestyle='-', linewidth=1.0)
    ty = (yc1+rs1) * np.ones(100)
    ax1.plot(tx, ty, color='k', linestyle='-', linewidth=1.0)
    ty = np.linspace(yc1-rs1,yc1+rs1,100)
    tx = (zc1-0.5*slice_length) * np.ones(100)
    ax1.plot(tx, ty, color='k', linestyle='-', linewidth=1.0)
    tx = (zc1+0.5*slice_length) * np.ones(100)
    ax1.plot(tx, ty, color='k', linestyle='-', linewidth=1.0)
    plt.text(zc1, yc1-0.4*rs1, r"${\rm S%01d}$" % (reg_ind+1), horizontalalignment='center', fontname='times', fontsize = fsize, fontweight = fweight, transform = ax1.transData)

ax1.set_xlim([-boxx/2.0,boxx/2.0])
ax1.set_ylim([-boxy/2.0,boxy/2.0])
fname = FigDir+"/snap_%03d_sheet_regions_NH.png" % (SnapNum)
print(fname)
savefig(fname)
plt.close()
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
'''

'''
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
# Test the centres and shock radii by making face-on temperature and density maps of the slices #
# Also try recomputing the centers using shrinking cylinders, no FoF groups, and only cold and dense gas (T<T_thresh, nH>D_thresh) #

boxx = 200.
boxy = 200.
boxz = 30.
pixres = 0.5

imgcent = [0,0,0]
xsize = 6.3
ysize = 6.3
fsize = 18
fsize2 = 22
fweight = 150

ax1pos = [0.14, 0.14, 0.67, 0.67]
caxpos = [0.82, 0.14, 0.025, 0.67]

xtick    = [-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100]
ytick    = [-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100]
dx_minor = 5
dy_minor = 5
xlaby    = -0.12
ylabx    = -0.12
tit_pos  = [0.2, 1.02]

xstr   = r"$x\:[{\rm kpc}]$"
ystr   = r"$y\:[{\rm kpc}]$"
titstr = r"${\rm " + SimLabel + r",\:\:\: z=" + '{0:.2f}'.format(redshift) + "}$"

s.data['DWght'] = rho
s.data['DTemp'] = rho * Temp
Hyd_dens        = s.rho * s.gmet[:,0] * UnitNumberDensity_in_cgs	# nH  cm^{-3}
s.data['DDens'] = rho * Hyd_dens

# Temp_thresh = 5.0e4 # Use gas with T<Temp_thresh to compute center
Temp_thresh = 10.0**(4.5) # Use gas with T<Temp_thresh to compute center
Dens_thresh = 10.0**(-3.0) # Use gas with T<Temp_thresh to compute center
rad_max = 30.0 # Maximal (initial) radius for shrinking cylinders
rad_min = 3.0  # Minimal (final) radius for shrinking cylinders
min_Npart = 32 # Exclude gas in halos with >min_Npart DM particles when computing center

### Slices 1-4 ###
for slice_ind in np.arange(0,4):
    xc1 = slices_1[slice_ind][0]
    yc1 = slices_1[slice_ind][1]
    zc1 = slices_1[slice_ind][2]
    rs1 = r_shock_list[slice_ind]

    s.data['pos'][:,0] = xx1 - xc1
    s.data['pos'][:,1] = yy1 - yc1
    s.data['pos'][:,2] = zz1 - zc1

    xc2 = 0.0
    yc2 = 0.0
    zc2 = 0.0
    zz2 = np.abs( s.data['pos'][:,2] )
    rad_check = rad_max
    mass_no_halos = 1.0 * mass[:]
    first_cell = 0
    for g_ind in np.arange(0,Ng):
        Ncell = sf.data['flty'][g_ind,0]
        last_cell = first_cell+Ncell
        if sf.data['flty'][g_ind,1]>=min_Npart:
            mass_no_halos[first_cell:last_cell] = 0.0
        first_cell += Ncell
    while rad_check>=rad_min:
        xx2 = s.data['pos'][:,0] - xc2
        yy2 = s.data['pos'][:,1] - yc2
        rr2 = np.sqrt( xx2**2 + yy2**2 )
        #b, = np.where( (rr2<=rad_check) & (zz2 <= L_slice) & (Temp<=Temp_thresh) )
        b, = np.where( (rr2<=rad_check) & (zz2 <= L_slice) & (Temp<=Temp_thresh) & (Hyd_dens>=Dens_thresh) )
        xc2 = np.sum( mass_no_halos[b]*s.data['pos'][b,0] ) / np.sum( mass_no_halos[b] )
        yc2 = np.sum( mass_no_halos[b]*s.data['pos'][b,1] ) / np.sum( mass_no_halos[b] )
        rad_check /= 1.1
    
    rad_check *= 1.1
    xx2 = s.data['pos'][:,0] - xc2
    yy2 = s.data['pos'][:,1] - yc2
    rr2 = np.sqrt( xx2**2 + yy2**2 )
    # b, = np.where( (rr2<=rad_check) & (zz2 <= L_slice) & (Temp<=Temp_thresh) )
    b, = np.where( (rr2<=rad_check) & (zz2 <= L_slice) & (Temp<=Temp_thresh) & (Hyd_dens>=Dens_thresh) )
    vx_cen_2 = np.sum( mass_no_halos[b]*vxx1[b] ) / np.sum( mass_no_halos[b] )
    vy_cen_2 = np.sum( mass_no_halos[b]*vyy1[b] ) / np.sum( mass_no_halos[b] )
    vz_cen_2 = np.sum( mass_no_halos[b]*vzz1[b] ) / np.sum( mass_no_halos[b] )

    xx2 = s.data['pos'][:,0]
    yy2 = s.data['pos'][:,1]
    rr2 = np.sqrt( xx2**2 + yy2**2 )
    b, = np.where( (rr2<=1.0) & (zz2 <= L_slice) )
    vx_cen_1 = np.sum( mass[b]*vxx1[b] ) / np.sum( mass[b] )
    vy_cen_1 = np.sum( mass[b]*vyy1[b] ) / np.sum( mass[b] )
    vz_cen_1 = np.sum( mass[b]*vzz1[b] ) / np.sum( mass[b] )

    print( ' ' )
    print( 'SLICE NUMBER', slice_ind+1 )
    print( 'New Center (x, y) = ', xc2, yc2 )
    print( 'New Center Distance = ', np.sqrt(xc2**2 + yc2**2) )
    print( 'Old CoM Velocity (vx, vy, vz) = ', vx_cen_1, vy_cen_1, vz_cen_1 )
    print( 'New CoM Velocity (vx, vy, vz) = ', vx_cen_2, vy_cen_2, vz_cen_2 )
    print( 'CoM Velocity Difference = ', np.sqrt( (vx_cen_1-vx_cen_2)**2 + (vy_cen_1-vy_cen_2)**2 + (vz_cen_1-vz_cen_2)**2 ) )
    print( ' ' )

    ##############
    proj_xy_Q = s.get_Aslice("DTemp", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz, axes = [0,1], proj = True, numthreads=32)
    proj_xy_W = s.get_Aslice("DWght", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz, axes = [0,1], proj = True, numthreads=32)
    mapp_xy   = np.transpose(proj_xy_Q['grid'] / proj_xy_W['grid'])
    ma = 1.0e6
    mi = 1.0e4

    xpl = proj_xy_Q['x']
    ypl = proj_xy_Q['y']

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    fig = plt.figure(figsize = (xsize,ysize))
    ax1 = plt.axes(ax1pos)
    plt.pcolormesh( xpl, ypl, mapp_xy, norm = matplotlib.colors.LogNorm(vmin=mi, vmax=ma), cmap = 'RdBu_r', rasterized = True )
    xlabel(xstr, fontname='times', fontsize = fsize2, fontweight = fweight)
    ylabel(ystr, fontname='times', fontsize = fsize2, fontweight = fweight)

    plt.xticks(xtick, fontname='arial', fontsize = fsize, fontweight = fweight)
    plt.yticks(ytick, fontname='arial', fontsize = fsize, fontweight = fweight)

    ax1.xaxis.set_tick_params(which = 'major', direction = 'out', width = 1, length = 6, pad = 5, bottom = True, top = True)
    ax1.yaxis.set_tick_params(which = 'major', direction = 'out', width = 1, length = 6, pad = 5, bottom = True, top = True)

    xminorLocator = MultipleLocator(dx_minor)
    ax1.xaxis.set_minor_locator(xminorLocator)
    ax1.xaxis.set_tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5, bottom = True, top = True)

    yminorLocator = MultipleLocator(dy_minor)
    ax1.yaxis.set_minor_locator(yminorLocator)
    ax1.yaxis.set_tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5, left = True, right = True)

    ax1.xaxis.set_label_coords(0.5, xlaby)
    ax1.yaxis.set_label_coords(ylabx, 0.5)

    ttl = ax1.set_title(titstr, fontname='times', fontsize = fsize2, fontweight = fweight, loc='left')
    ttl.set_position(tit_pos)

    axcbar = fig.add_axes(caxpos)
    cbar = matplotlib.colorbar.ColorbarBase(ax = axcbar, norm = matplotlib.colors.LogNorm(vmin=mi, vmax=ma), cmap = 'RdBu_r', orientation='vertical', ticks = LogLocator(subs=range(10)))
    cbar.set_label(label=r"$T\:[{\rm K}]$", fontname='times', fontsize = fsize2, fontweight = fweight, rotation=270, labelpad=20)
    cbar.ax.tick_params(which = 'major', direction = 'out', labelsize = fsize, width = 1, length = 6, pad = 5)
    cbar.ax.tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5)

    t=np.linspace(0.0,2.0*pi,100)
    ax1.plot(rs1*np.cos(t), rs1*np.sin(t), color='w', linestyle='--', linewidth=1.0)
    ax1.plot(xc2, yc2, color='w', marker='o', markersize=10.0, markerfacecolor='w', markeredgecolor='w')
    for i in range(Ng):
        if gR[i]>0.:
            if (gxx1[i]-xc1)+gR[i]>=-boxx/2.0 and (gxx1[i]-xc1)-gR[i]<=boxx/2.0 and (gyy1[i]-yc1)+gR[i]>=-boxy/2.0 and (gyy1[i]-yc1)-gR[i]<=boxy/2.0 and (gzz1[i]-zc1)>=-boxz/2.0 and (gzz1[i]-zc1)<=boxz/2.0:
    	        if gM[i]>1.0:
                    ax1.plot(gxx1[i]-xc1+gR[i]*np.cos(t), gyy1[i]-yc1+gR[i]*np.sin(t), color='g', linestyle='-', linewidth=1.0)
    	        elif gM[i]>0.1:
                    ax1.plot(gxx1[i]-xc1+gR[i]*np.cos(t), gyy1[i]-yc1+gR[i]*np.sin(t), color='k', linestyle='-', linewidth=1.0)
    	        elif gM[i]>0.01:
                    ax1.plot(gxx1[i]-xc1+gR[i]*np.cos(t), gyy1[i]-yc1+gR[i]*np.sin(t), color='w', linestyle='-', linewidth=1.0)

    ax1.set_xlim([-boxx/2.0,boxx/2.0])
    ax1.set_ylim([-boxy/2.0,boxy/2.0])
    fname = FigDir+"/slice_%03d_temperature.png" % (slice_ind+1)
    print(fname)
    savefig(fname)
    plt.close()

    ##############
    proj_xy_Q = s.get_Aslice("DDens", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz, axes = [0,1], proj = True, numthreads=32)
    proj_xy_W = s.get_Aslice("DWght", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz, axes = [0,1], proj = True, numthreads=32)
    mapp_xy   = np.transpose(proj_xy_Q['grid'] / proj_xy_W['grid'])
    ma = 3.0e-2
    mi = 1.0e-5

    xpl = proj_xy_Q['x']
    ypl = proj_xy_Q['y']

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    fig = plt.figure(figsize = (xsize,ysize))
    ax1 = plt.axes(ax1pos)
    plt.pcolormesh( xpl, ypl, mapp_xy, norm = matplotlib.colors.LogNorm(vmin=mi, vmax=ma), cmap = 'RdBu_r', rasterized = True )
    xlabel(xstr, fontname='times', fontsize = fsize2, fontweight = fweight)
    ylabel(ystr, fontname='times', fontsize = fsize2, fontweight = fweight)

    plt.xticks(xtick, fontname='arial', fontsize = fsize, fontweight = fweight)
    plt.yticks(ytick, fontname='arial', fontsize = fsize, fontweight = fweight)

    ax1.xaxis.set_tick_params(which = 'major', direction = 'out', width = 1, length = 6, pad = 5, bottom = True, top = True)
    ax1.yaxis.set_tick_params(which = 'major', direction = 'out', width = 1, length = 6, pad = 5, bottom = True, top = True)

    xminorLocator = MultipleLocator(dx_minor)
    ax1.xaxis.set_minor_locator(xminorLocator)
    ax1.xaxis.set_tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5, bottom = True, top = True)

    yminorLocator = MultipleLocator(dy_minor)
    ax1.yaxis.set_minor_locator(yminorLocator)
    ax1.yaxis.set_tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5, left = True, right = True)

    ax1.xaxis.set_label_coords(0.5, xlaby)
    ax1.yaxis.set_label_coords(ylabx, 0.5)

    ttl = ax1.set_title(titstr, fontname='times', fontsize = fsize2, fontweight = fweight, loc='left')
    ttl.set_position(tit_pos)

    axcbar = fig.add_axes(caxpos)
    cbar = matplotlib.colorbar.ColorbarBase(ax = axcbar, norm = matplotlib.colors.LogNorm(vmin=mi, vmax=ma), cmap = 'RdBu_r', orientation='vertical', ticks = LogLocator(subs=range(10)))
    cbar.set_label(label=r"$n_{\rm H}\:[{\rm cm^{-3}}]$", fontname='times', fontsize = fsize2, fontweight = fweight, rotation=270, labelpad=20)
    cbar.ax.tick_params(which = 'major', direction = 'out', labelsize = fsize, width = 1, length = 6, pad = 5)
    cbar.ax.tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5)

    t=np.linspace(0.0,2.0*pi,100)
    ax1.plot(rs1*np.cos(t), rs1*np.sin(t), color='w', linestyle='--', linewidth=1.0)
    ax1.plot(xc2, yc2, color='w', marker='o', markersize=10.0, markerfacecolor='w', markeredgecolor='w')
    for i in range(Ng):
        if gR[i]>0.:
            if (gxx1[i]-xc1)+gR[i]>=-boxx/2.0 and (gxx1[i]-xc1)-gR[i]<=boxx/2.0 and (gyy1[i]-yc1)+gR[i]>=-boxy/2.0 and (gyy1[i]-yc1)-gR[i]<=boxy/2.0 and (gzz1[i]-zc1)>=-boxz/2.0 and (gzz1[i]-zc1)<=boxz/2.0:
    	        if gM[i]>1.0:
                    ax1.plot(gxx1[i]-xc1+gR[i]*np.cos(t), gyy1[i]-yc1+gR[i]*np.sin(t), color='g', linestyle='-', linewidth=1.0)
    	        elif gM[i]>0.1:
                    ax1.plot(gxx1[i]-xc1+gR[i]*np.cos(t), gyy1[i]-yc1+gR[i]*np.sin(t), color='k', linestyle='-', linewidth=1.0)
    	        elif gM[i]>0.01:
                    ax1.plot(gxx1[i]-xc1+gR[i]*np.cos(t), gyy1[i]-yc1+gR[i]*np.sin(t), color='w', linestyle='-', linewidth=1.0)

    ax1.set_xlim([-boxx/2.0,boxx/2.0])
    ax1.set_ylim([-boxy/2.0,boxy/2.0])
    fname = FigDir+"/slice_%03d_density.png" % (slice_ind+1)
    print(fname)
    savefig(fname)
    plt.close()

### Slices 5-6 ###
for slice_ind in np.arange(4,6):
    xc1 = slices_1[slice_ind][0]
    yc1 = slices_1[slice_ind][1]
    zc1 = slices_1[slice_ind][2]
    rs1 = r_shock_list[slice_ind]

    s.data['pos'][:,0] = xx1 - xc1
    s.data['pos'][:,1] = yy1 - yc1
    s.data['pos'][:,2] = zz1 - zc1

    xc2 = 0.0
    yc2 = 0.0
    zc2 = 0.0
    zz2 = np.abs( s.data['pos'][:,1] )
    rad_check = rad_max
    mass_no_halos = 1.0 * mass[:]
    first_cell = 0
    for g_ind in np.arange(0,Ng):
        Ncell = sf.data['flty'][g_ind,0]
        last_cell = first_cell+Ncell
        if sf.data['flty'][g_ind,1]>=min_Npart:
            mass_no_halos[first_cell:last_cell] = 0.0
        first_cell += Ncell
    while rad_check>=rad_min:
        xx2 = s.data['pos'][:,0] - xc2
        yy2 = s.data['pos'][:,2] - yc2
        rr2 = np.sqrt( xx2**2 + yy2**2 )
        # b, = np.where( (rr2<=rad_check) & (zz2 <= L_slice) & (Temp<=Temp_thresh) )
        b, = np.where( (rr2<=rad_check) & (zz2 <= L_slice) & (Temp<=Temp_thresh) & (Hyd_dens>=Dens_thresh) )
        xc2 = np.sum( mass_no_halos[b]*s.data['pos'][b,0] ) / np.sum( mass_no_halos[b] )
        yc2 = np.sum( mass_no_halos[b]*s.data['pos'][b,2] ) / np.sum( mass_no_halos[b] )
        rad_check /= 1.1
    
    rad_check *= 1.1
    xx2 = s.data['pos'][:,0] - xc2
    yy2 = s.data['pos'][:,2] - yc2
    rr2 = np.sqrt( xx2**2 + yy2**2 )
    # b, = np.where( (rr2<=rad_check) & (zz2 <= L_slice) & (Temp<=Temp_thresh) )
    b, = np.where( (rr2<=rad_check) & (zz2 <= L_slice) & (Temp<=Temp_thresh) & (Hyd_dens>=Dens_thresh) )
    vx_cen_2 = np.sum( mass_no_halos[b]*vxx1[b] ) / np.sum( mass_no_halos[b] )
    vy_cen_2 = np.sum( mass_no_halos[b]*vzz1[b] ) / np.sum( mass_no_halos[b] )
    vz_cen_2 = np.sum( mass_no_halos[b]*vyy1[b] ) / np.sum( mass_no_halos[b] )

    xx2 = s.data['pos'][:,0]
    yy2 = s.data['pos'][:,2]
    rr2 = np.sqrt( xx2**2 + yy2**2 )
    b, = np.where( (rr2<=1.0) & (zz2 <= L_slice) )
    vx_cen_1 = np.sum( mass[b]*vxx1[b] ) / np.sum( mass[b] )
    vy_cen_1 = np.sum( mass[b]*vzz1[b] ) / np.sum( mass[b] )
    vz_cen_1 = np.sum( mass[b]*vyy1[b] ) / np.sum( mass[b] )

    print( ' ' )
    print( 'SLICE NUMBER', slice_ind+1 )
    print( 'New Center (x, y) = ', xc2, yc2 )
    print( 'New Center Distance = ', np.sqrt(xc2**2 + yc2**2) )
    print( 'Old CoM Velocity (vx, vy, vz) = ', vx_cen_1, vy_cen_1, vz_cen_1 )
    print( 'New CoM Velocity (vx, vy, vz) = ', vx_cen_2, vy_cen_2, vz_cen_2 )
    print( 'CoM Velocity Difference = ', np.sqrt( (vx_cen_1-vx_cen_2)**2 + (vy_cen_1-vy_cen_2)**2 + (vz_cen_1-vz_cen_2)**2 ) )
    print( ' ' )

    ##############
    proj_xy_Q = s.get_Aslice("DTemp", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz, axes = [0,2], proj = True, numthreads=32)
    proj_xy_W = s.get_Aslice("DWght", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz, axes = [0,2], proj = True, numthreads=32)
    mapp_xy   = np.transpose(proj_xy_Q['grid'] / proj_xy_W['grid'])
    ma = 1.0e6
    mi = 1.0e4

    xpl = proj_xy_Q['x']
    ypl = proj_xy_Q['y']

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    fig = plt.figure(figsize = (xsize,ysize))
    ax1 = plt.axes(ax1pos)
    plt.pcolormesh( xpl, ypl, mapp_xy, norm = matplotlib.colors.LogNorm(vmin=mi, vmax=ma), cmap = 'RdBu_r', rasterized = True )
    xlabel(xstr, fontname='times', fontsize = fsize2, fontweight = fweight)
    ylabel(ystr, fontname='times', fontsize = fsize2, fontweight = fweight)

    plt.xticks(xtick, fontname='arial', fontsize = fsize, fontweight = fweight)
    plt.yticks(ytick, fontname='arial', fontsize = fsize, fontweight = fweight)

    ax1.xaxis.set_tick_params(which = 'major', direction = 'out', width = 1, length = 6, pad = 5, bottom = True, top = True)
    ax1.yaxis.set_tick_params(which = 'major', direction = 'out', width = 1, length = 6, pad = 5, bottom = True, top = True)

    xminorLocator = MultipleLocator(dx_minor)
    ax1.xaxis.set_minor_locator(xminorLocator)
    ax1.xaxis.set_tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5, bottom = True, top = True)

    yminorLocator = MultipleLocator(dy_minor)
    ax1.yaxis.set_minor_locator(yminorLocator)
    ax1.yaxis.set_tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5, left = True, right = True)

    ax1.xaxis.set_label_coords(0.5, xlaby)
    ax1.yaxis.set_label_coords(ylabx, 0.5)

    ttl = ax1.set_title(titstr, fontname='times', fontsize = fsize2, fontweight = fweight, loc='left')
    ttl.set_position(tit_pos)

    axcbar = fig.add_axes(caxpos)
    cbar = matplotlib.colorbar.ColorbarBase(ax = axcbar, norm = matplotlib.colors.LogNorm(vmin=mi, vmax=ma), cmap = 'RdBu_r', orientation='vertical', ticks = LogLocator(subs=range(10)))
    cbar.set_label(label=r"$T\:[{\rm K}]$", fontname='times', fontsize = fsize2, fontweight = fweight, rotation=270, labelpad=20)
    cbar.ax.tick_params(which = 'major', direction = 'out', labelsize = fsize, width = 1, length = 6, pad = 5)
    cbar.ax.tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5)

    t=np.linspace(0.0,2.0*pi,100)
    ax1.plot(rs1*np.cos(t), rs1*np.sin(t), color='w', linestyle='--', linewidth=1.0)
    ax1.plot(xc2, yc2, color='w', marker='o', markersize=10.0, markerfacecolor='w', markeredgecolor='w')
    for i in range(Ng):
        if gR[i]>0.:
            if (gxx1[i]-xc1)+gR[i]>=-boxx/2.0 and (gxx1[i]-xc1)-gR[i]<=boxx/2.0 and (gzz1[i]-zc1)+gR[i]>=-boxy/2.0 and (gzz1[i]-zc1)-gR[i]<=boxy/2.0 and (gyy1[i]-yc1)>=-boxz/2.0 and (gyy1[i]-yc1)<=boxz/2.0:
    	        if gM[i]>1.0:
                    ax1.plot(gxx1[i]-xc1+gR[i]*np.cos(t), gzz1[i]-zc1+gR[i]*np.sin(t), color='g', linestyle='-', linewidth=1.0)
    	        elif gM[i]>0.1:
                    ax1.plot(gxx1[i]-xc1+gR[i]*np.cos(t), gzz1[i]-zc1+gR[i]*np.sin(t), color='k', linestyle='-', linewidth=1.0)
    	        elif gM[i]>0.01:
                    ax1.plot(gxx1[i]-xc1+gR[i]*np.cos(t), gzz1[i]-zc1+gR[i]*np.sin(t), color='w', linestyle='-', linewidth=1.0)

    ax1.set_xlim([-boxx/2.0,boxx/2.0])
    ax1.set_ylim([-boxy/2.0,boxy/2.0])
    fname = FigDir+"/slice_%03d_temperature.png" % (slice_ind+1)
    print(fname)
    savefig(fname)
    plt.close()

    ##############
    proj_xy_Q = s.get_Aslice("DDens", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz, axes = [0,2], proj = True, numthreads=32)
    proj_xy_W = s.get_Aslice("DWght", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz, axes = [0,2], proj = True, numthreads=32)
    mapp_xy   = np.transpose(proj_xy_Q['grid'] / proj_xy_W['grid'])
    ma = 3.0e-2
    mi = 1.0e-5

    xpl = proj_xy_Q['x']
    ypl = proj_xy_Q['y']

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    fig = plt.figure(figsize = (xsize,ysize))
    ax1 = plt.axes(ax1pos)
    plt.pcolormesh( xpl, ypl, mapp_xy, norm = matplotlib.colors.LogNorm(vmin=mi, vmax=ma), cmap = 'RdBu_r', rasterized = True )
    xlabel(xstr, fontname='times', fontsize = fsize2, fontweight = fweight)
    ylabel(ystr, fontname='times', fontsize = fsize2, fontweight = fweight)

    plt.xticks(xtick, fontname='arial', fontsize = fsize, fontweight = fweight)
    plt.yticks(ytick, fontname='arial', fontsize = fsize, fontweight = fweight)

    ax1.xaxis.set_tick_params(which = 'major', direction = 'out', width = 1, length = 6, pad = 5, bottom = True, top = True)
    ax1.yaxis.set_tick_params(which = 'major', direction = 'out', width = 1, length = 6, pad = 5, bottom = True, top = True)

    xminorLocator = MultipleLocator(dx_minor)
    ax1.xaxis.set_minor_locator(xminorLocator)
    ax1.xaxis.set_tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5, bottom = True, top = True)

    yminorLocator = MultipleLocator(dy_minor)
    ax1.yaxis.set_minor_locator(yminorLocator)
    ax1.yaxis.set_tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5, left = True, right = True)

    ax1.xaxis.set_label_coords(0.5, xlaby)
    ax1.yaxis.set_label_coords(ylabx, 0.5)

    ttl = ax1.set_title(titstr, fontname='times', fontsize = fsize2, fontweight = fweight, loc='left')
    ttl.set_position(tit_pos)

    axcbar = fig.add_axes(caxpos)
    cbar = matplotlib.colorbar.ColorbarBase(ax = axcbar, norm = matplotlib.colors.LogNorm(vmin=mi, vmax=ma), cmap = 'RdBu_r', orientation='vertical', ticks = LogLocator(subs=range(10)))
    cbar.set_label(label=r"$n_{\rm H}\:[{\rm cm^{-3}}]$", fontname='times', fontsize = fsize2, fontweight = fweight, rotation=270, labelpad=20)
    cbar.ax.tick_params(which = 'major', direction = 'out', labelsize = fsize, width = 1, length = 6, pad = 5)
    cbar.ax.tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5)

    t=np.linspace(0.0,2.0*pi,100)
    ax1.plot(rs1*np.cos(t), rs1*np.sin(t), color='w', linestyle='--', linewidth=1.0)
    ax1.plot(xc2, yc2, color='w', marker='o', markersize=10.0, markerfacecolor='w', markeredgecolor='w')
    for i in range(Ng):
        if gR[i]>0.:
            if (gxx1[i]-xc1)+gR[i]>=-boxx/2.0 and (gxx1[i]-xc1)-gR[i]<=boxx/2.0 and (gzz1[i]-zc1)+gR[i]>=-boxy/2.0 and (gzz1[i]-zc1)-gR[i]<=boxy/2.0 and (gyy1[i]-yc1)>=-boxz/2.0 and (gyy1[i]-yc1)<=boxz/2.0:
    	        if gM[i]>1.0:
                    ax1.plot(gxx1[i]-xc1+gR[i]*np.cos(t), gzz1[i]-zc1+gR[i]*np.sin(t), color='g', linestyle='-', linewidth=1.0)
    	        elif gM[i]>0.1:
                    ax1.plot(gxx1[i]-xc1+gR[i]*np.cos(t), gzz1[i]-zc1+gR[i]*np.sin(t), color='k', linestyle='-', linewidth=1.0)
    	        elif gM[i]>0.01:
                    ax1.plot(gxx1[i]-xc1+gR[i]*np.cos(t), gzz1[i]-zc1+gR[i]*np.sin(t), color='w', linestyle='-', linewidth=1.0)

    ax1.set_xlim([-boxx/2.0,boxx/2.0])
    ax1.set_ylim([-boxy/2.0,boxy/2.0])
    fname = FigDir+"/slice_%03d_density.png" % (slice_ind+1)
    print(fname)
    savefig(fname)
    plt.close()

### Slices 7-10 ###
for slice_ind in np.arange(0,4):
    xc1 = slices_2[slice_ind][0]
    yc1 = slices_2[slice_ind][1]
    zc1 = slices_2[slice_ind][2]
    rs1 = r_shock_list[slice_ind+6]

    s.data['pos'][:,0] = xx7 - xc1
    s.data['pos'][:,1] = yy7 - yc1
    s.data['pos'][:,2] = zz7 - zc1

    xc2 = 0.0
    yc2 = 0.0
    zc2 = 0.0
    zz2 = np.abs( s.data['pos'][:,2] )
    rad_check = rad_max
    mass_no_halos = 1.0 * mass[:]
    first_cell = 0
    for g_ind in np.arange(0,Ng):
        Ncell = sf.data['flty'][g_ind,0]
        last_cell = first_cell+Ncell
        if sf.data['flty'][g_ind,1]>=min_Npart:
            mass_no_halos[first_cell:last_cell] = 0.0
        first_cell += Ncell
    while rad_check>=rad_min:
        xx2 = s.data['pos'][:,0] - xc2
        yy2 = s.data['pos'][:,1] - yc2
        rr2 = np.sqrt( xx2**2 + yy2**2 )
        # b, = np.where( (rr2<=rad_check) & (zz2 <= L_slice) & (Temp<=Temp_thresh) )
        b, = np.where( (rr2<=rad_check) & (zz2 <= L_slice) & (Temp<=Temp_thresh) & (Hyd_dens>=Dens_thresh) )
        xc2 = np.sum( mass_no_halos[b]*s.data['pos'][b,0] ) / np.sum( mass_no_halos[b] )
        yc2 = np.sum( mass_no_halos[b]*s.data['pos'][b,1] ) / np.sum( mass_no_halos[b] )
        rad_check /= 1.1
    
    rad_check *= 1.1
    xx2 = s.data['pos'][:,0] - xc2
    yy2 = s.data['pos'][:,1] - yc2
    rr2 = np.sqrt( xx2**2 + yy2**2 )
    # b, = np.where( (rr2<=rad_check) & (zz2 <= L_slice) & (Temp<=Temp_thresh) )
    b, = np.where( (rr2<=rad_check) & (zz2 <= L_slice) & (Temp<=Temp_thresh) & (Hyd_dens>=Dens_thresh) )
    vx_cen_2 = np.sum( mass_no_halos[b]*vxx7[b] ) / np.sum( mass_no_halos[b] )
    vy_cen_2 = np.sum( mass_no_halos[b]*vyy7[b] ) / np.sum( mass_no_halos[b] )
    vz_cen_2 = np.sum( mass_no_halos[b]*vzz7[b] ) / np.sum( mass_no_halos[b] )

    xx2 = s.data['pos'][:,0]
    yy2 = s.data['pos'][:,1]
    rr2 = np.sqrt( xx2**2 + yy2**2 )
    b, = np.where( (rr2<=1.0) & (zz2 <= L_slice) )
    vx_cen_1 = np.sum( mass[b]*vxx7[b] ) / np.sum( mass[b] )
    vy_cen_1 = np.sum( mass[b]*vyy7[b] ) / np.sum( mass[b] )
    vz_cen_1 = np.sum( mass[b]*vzz7[b] ) / np.sum( mass[b] )

    print( ' ' )
    print( 'SLICE NUMBER', slice_ind+7 )
    print( 'New Center (x, y) = ', xc2, yc2 )
    print( 'New Center Distance = ', np.sqrt(xc2**2 + yc2**2) )
    print( 'Old CoM Velocity (vx, vy, vz) = ', vx_cen_1, vy_cen_1, vz_cen_1 )
    print( 'New CoM Velocity (vx, vy, vz) = ', vx_cen_2, vy_cen_2, vz_cen_2 )
    print( 'CoM Velocity Difference = ', np.sqrt( (vx_cen_1-vx_cen_2)**2 + (vy_cen_1-vy_cen_2)**2 + (vz_cen_1-vz_cen_2)**2 ) )
    print( ' ' )

    ##############
    proj_xy_Q = s.get_Aslice("DTemp", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz, axes = [0,1], proj = True, numthreads=32)
    proj_xy_W = s.get_Aslice("DWght", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz, axes = [0,1], proj = True, numthreads=32)
    mapp_xy   = np.transpose(proj_xy_Q['grid'] / proj_xy_W['grid'])
    ma = 1.0e6
    mi = 1.0e4

    xpl = proj_xy_Q['x']
    ypl = proj_xy_Q['y']

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    fig = plt.figure(figsize = (xsize,ysize))
    ax1 = plt.axes(ax1pos)
    plt.pcolormesh( xpl, ypl, mapp_xy, norm = matplotlib.colors.LogNorm(vmin=mi, vmax=ma), cmap = 'RdBu_r', rasterized = True )
    xlabel(xstr, fontname='times', fontsize = fsize2, fontweight = fweight)
    ylabel(ystr, fontname='times', fontsize = fsize2, fontweight = fweight)

    plt.xticks(xtick, fontname='arial', fontsize = fsize, fontweight = fweight)
    plt.yticks(ytick, fontname='arial', fontsize = fsize, fontweight = fweight)

    ax1.xaxis.set_tick_params(which = 'major', direction = 'out', width = 1, length = 6, pad = 5, bottom = True, top = True)
    ax1.yaxis.set_tick_params(which = 'major', direction = 'out', width = 1, length = 6, pad = 5, bottom = True, top = True)

    xminorLocator = MultipleLocator(dx_minor)
    ax1.xaxis.set_minor_locator(xminorLocator)
    ax1.xaxis.set_tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5, bottom = True, top = True)

    yminorLocator = MultipleLocator(dy_minor)
    ax1.yaxis.set_minor_locator(yminorLocator)
    ax1.yaxis.set_tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5, left = True, right = True)

    ax1.xaxis.set_label_coords(0.5, xlaby)
    ax1.yaxis.set_label_coords(ylabx, 0.5)

    ttl = ax1.set_title(titstr, fontname='times', fontsize = fsize2, fontweight = fweight, loc='left')
    ttl.set_position(tit_pos)

    axcbar = fig.add_axes(caxpos)
    cbar = matplotlib.colorbar.ColorbarBase(ax = axcbar, norm = matplotlib.colors.LogNorm(vmin=mi, vmax=ma), cmap = 'RdBu_r', orientation='vertical', ticks = LogLocator(subs=range(10)))
    cbar.set_label(label=r"$T\:[{\rm K}]$", fontname='times', fontsize = fsize2, fontweight = fweight, rotation=270, labelpad=20)
    cbar.ax.tick_params(which = 'major', direction = 'out', labelsize = fsize, width = 1, length = 6, pad = 5)
    cbar.ax.tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5)

    t=np.linspace(0.0,2.0*pi,100)
    ax1.plot(rs1*np.cos(t), rs1*np.sin(t), color='w', linestyle='--', linewidth=1.0)
    ax1.plot(xc2, yc2, color='w', marker='o', markersize=10.0, markerfacecolor='w', markeredgecolor='w')
    for i in range(Ng):
        if gR[i]>0.:
            if (gxx7[i]-xc1)+gR[i]>=-boxx/2.0 and (gxx7[i]-xc1)-gR[i]<=boxx/2.0 and (gyy7[i]-yc1)+gR[i]>=-boxy/2.0 and (gyy7[i]-yc1)-gR[i]<=boxy/2.0 and (gzz7[i]-zc1)>=-boxz/2.0 and (gzz7[i]-zc1)<=boxz/2.0:
    	        if gM[i]>1.0:
                    ax1.plot(gxx7[i]-xc1+gR[i]*np.cos(t), gyy7[i]-yc1+gR[i]*np.sin(t), color='g', linestyle='-', linewidth=1.0)
    	        elif gM[i]>0.1:
                    ax1.plot(gxx7[i]-xc1+gR[i]*np.cos(t), gyy7[i]-yc1+gR[i]*np.sin(t), color='k', linestyle='-', linewidth=1.0)
    	        elif gM[i]>0.01:
                    ax1.plot(gxx7[i]-xc1+gR[i]*np.cos(t), gyy7[i]-yc1+gR[i]*np.sin(t), color='w', linestyle='-', linewidth=1.0)

    ax1.set_xlim([-boxx/2.0,boxx/2.0])
    ax1.set_ylim([-boxy/2.0,boxy/2.0])
    fname = FigDir+"/slice_%03d_temperature.png" % (slice_ind+7)
    print(fname)
    savefig(fname)
    plt.close()

    ##############
    proj_xy_Q = s.get_Aslice("DDens", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz, axes = [0,1], proj = True, numthreads=32)
    proj_xy_W = s.get_Aslice("DWght", box = [boxx,boxy], center = imgcent, nx = int(boxx/pixres), ny = int(boxy/pixres), nz = int(boxz/pixres), boxz = boxz, axes = [0,1], proj = True, numthreads=32)
    mapp_xy   = np.transpose(proj_xy_Q['grid'] / proj_xy_W['grid'])
    ma = 3.0e-2
    mi = 1.0e-5

    xpl = proj_xy_Q['x']
    ypl = proj_xy_Q['y']

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    fig = plt.figure(figsize = (xsize,ysize))
    ax1 = plt.axes(ax1pos)
    plt.pcolormesh( xpl, ypl, mapp_xy, norm = matplotlib.colors.LogNorm(vmin=mi, vmax=ma), cmap = 'RdBu_r', rasterized = True )
    xlabel(xstr, fontname='times', fontsize = fsize2, fontweight = fweight)
    ylabel(ystr, fontname='times', fontsize = fsize2, fontweight = fweight)

    plt.xticks(xtick, fontname='arial', fontsize = fsize, fontweight = fweight)
    plt.yticks(ytick, fontname='arial', fontsize = fsize, fontweight = fweight)

    ax1.xaxis.set_tick_params(which = 'major', direction = 'out', width = 1, length = 6, pad = 5, bottom = True, top = True)
    ax1.yaxis.set_tick_params(which = 'major', direction = 'out', width = 1, length = 6, pad = 5, bottom = True, top = True)

    xminorLocator = MultipleLocator(dx_minor)
    ax1.xaxis.set_minor_locator(xminorLocator)
    ax1.xaxis.set_tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5, bottom = True, top = True)

    yminorLocator = MultipleLocator(dy_minor)
    ax1.yaxis.set_minor_locator(yminorLocator)
    ax1.yaxis.set_tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5, left = True, right = True)

    ax1.xaxis.set_label_coords(0.5, xlaby)
    ax1.yaxis.set_label_coords(ylabx, 0.5)

    ttl = ax1.set_title(titstr, fontname='times', fontsize = fsize2, fontweight = fweight, loc='left')
    ttl.set_position(tit_pos)

    axcbar = fig.add_axes(caxpos)
    cbar = matplotlib.colorbar.ColorbarBase(ax = axcbar, norm = matplotlib.colors.LogNorm(vmin=mi, vmax=ma), cmap = 'RdBu_r', orientation='vertical', ticks = LogLocator(subs=range(10)))
    cbar.set_label(label=r"$n_{\rm H}\:[{\rm cm^{-3}}]$", fontname='times', fontsize = fsize2, fontweight = fweight, rotation=270, labelpad=20)
    cbar.ax.tick_params(which = 'major', direction = 'out', labelsize = fsize, width = 1, length = 6, pad = 5)
    cbar.ax.tick_params(which = 'minor', direction = 'out', width = 1, length = 3, pad = 5)

    t=np.linspace(0.0,2.0*pi,100)
    ax1.plot(rs1*np.cos(t), rs1*np.sin(t), color='w', linestyle='--', linewidth=1.0)
    ax1.plot(xc2, yc2, color='w', marker='o', markersize=10.0, markerfacecolor='w', markeredgecolor='w')
    for i in range(Ng):
        if gR[i]>0.:
            if (gxx7[i]-xc1)+gR[i]>=-boxx/2.0 and (gxx7[i]-xc1)-gR[i]<=boxx/2.0 and (gyy7[i]-yc1)+gR[i]>=-boxy/2.0 and (gyy7[i]-yc1)-gR[i]<=boxy/2.0 and (gzz7[i]-zc1)>=-boxz/2.0 and (gzz7[i]-zc1)<=boxz/2.0:
    	        if gM[i]>1.0:
                    ax1.plot(gxx7[i]-xc1+gR[i]*np.cos(t), gyy7[i]-yc1+gR[i]*np.sin(t), color='g', linestyle='-', linewidth=1.0)
    	        elif gM[i]>0.1:
                    ax1.plot(gxx7[i]-xc1+gR[i]*np.cos(t), gyy7[i]-yc1+gR[i]*np.sin(t), color='k', linestyle='-', linewidth=1.0)
    	        elif gM[i]>0.01:
                    ax1.plot(gxx7[i]-xc1+gR[i]*np.cos(t), gyy7[i]-yc1+gR[i]*np.sin(t), color='w', linestyle='-', linewidth=1.0)

    ax1.set_xlim([-boxx/2.0,boxx/2.0])
    ax1.set_ylim([-boxy/2.0,boxy/2.0])
    fname = FigDir+"/slice_%03d_density.png" % (slice_ind+7)
    print(fname)
    savefig(fname)
    plt.close()
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
'''

###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################

### Slices 1-10 ###
for slice_ind in np.arange(0,10):
    print(' ')
    print('SLICE NUMBER ', slice_ind+1)

#################################################################################################
### Initialize shock radius, center, and rest frame
### Consider recalibrating the centres using the centre of mass of cold gas, T<~5x10^4 K, in shrinking cylinders
### This will also require recalibrating the profiles and shock radii, but likely will not change these by much. However, it can influence the dynamics and pressure gradients
### Should also compute rest frame self-consistently when finding the center itself, in the final cylinder of the shrinking cylinders method, as opposed to what is done here

    r_shock = r_shock_list[slice_ind]
    Rmin    = Rmin_norm * r_shock
    Rmax    = Rmax_norm * r_shock
    Delta_r = Delta_r_norm * r_shock
    r_bin   = r_bin_norm * r_shock
    print('r_shock = ', r_shock, 'kpc')
    print('Delta_r, Delta_p, Delta_z = ', Delta_r, Delta_p, Delta_z)
    print('r_bin = ', r_bin)
    print('p_bin = ', p_bin)
    print('z_bin = ', z_bin)


    if slice_ind<6:
        xcen = slices_1[slice_ind][0]
        ycen = slices_1[slice_ind][1]
        zcen = slices_1[slice_ind][2]
    elif slice_ind<10:
        xcen = slices_2[slice_ind-6][0]
        ycen = slices_2[slice_ind-6][1]
        zcen = slices_2[slice_ind-6][2]
    print('xcen, ycen, zcen = ', xcen, ycen, zcen, 'kpc')

    if slice_ind<4:
        cyl_rad   = np.sqrt( (xx1 - xcen)**2   + (yy1 - ycen)**2   ) # Cylindrical radius for gas
        dcyl_rad  = np.sqrt( (dxx1 - xcen)**2  + (dyy1 - ycen)**2  ) # Cylindrical radius for dark matter
        stcyl_rad = np.sqrt( (stxx1 - xcen)**2 + (styy1 - ycen)**2 ) # Cylindrical radius for stars
        # This computes the center of mass velocity of the gas within a 1kpc cylinder around the selected center
        CoM_r, = np.where( (cyl_rad <= 1.0) & (np.abs(zz1 - zcen) <= L_slice) )
        CoM_vx = sum( vxx1[CoM_r].astype('float64')*mass[CoM_r].astype('float64') ) / sum( mass[CoM_r].astype('float64') )
        CoM_vy = sum( vyy1[CoM_r].astype('float64')*mass[CoM_r].astype('float64') ) / sum( mass[CoM_r].astype('float64') )
        CoM_vz = sum( vzz1[CoM_r].astype('float64')*mass[CoM_r].astype('float64') ) / sum( mass[CoM_r].astype('float64') )

    elif slice_ind<6:   # In slices 5 and 6, the filament axis is the y axis, so the filament radius is in the XZ plane
        cyl_rad   = np.sqrt( (xx1 - xcen)**2   + (zz1 - zcen)**2   ) # Cylindrical radius for gas
        dcyl_rad  = np.sqrt( (dxx1 - xcen)**2  + (dzz1 - zcen)**2  ) # Cylindrical radius for dark matter
        stcyl_rad = np.sqrt( (stxx1 - xcen)**2 + (stzz1 - zcen)**2 ) # Cylindrical radius for stars
        # This computes the center of mass velocity of the gas within a 1kpc cylinder around the selected center
        CoM_r, = np.where( (cyl_rad <= 1.0) & (np.abs(yy1 - ycen) <= L_slice) ) # In slices 5 and 6, the filament axis is the y axis
        CoM_vx = sum( vxx1[CoM_r].astype('float64')*mass[CoM_r].astype('float64') ) / sum( mass[CoM_r].astype('float64') )
        CoM_vy = sum( vzz1[CoM_r].astype('float64')*mass[CoM_r].astype('float64') ) / sum( mass[CoM_r].astype('float64') )
        CoM_vz = sum( vyy1[CoM_r].astype('float64')*mass[CoM_r].astype('float64') ) / sum( mass[CoM_r].astype('float64') )

    elif slice_ind<10:
        cyl_rad   = np.sqrt( (xx7 - xcen)**2   + (yy7 - ycen)**2   ) # Cylindrical radius for gas
        dcyl_rad  = np.sqrt( (dxx7 - xcen)**2  + (dyy7 - ycen)**2  ) # Cylindrical radius for dark matter
        stcyl_rad = np.sqrt( (stxx7 - xcen)**2 + (styy7 - ycen)**2 ) # Cylindrical radius for stars
        # This computes the center of mass velocity of the gas within a 1kpc cylinder around the selected center
        CoM_r, = np.where( (cyl_rad <= 1.0) & (np.abs(zz7 - zcen) <= L_slice) )
        CoM_vx = sum( vxx7[CoM_r].astype('float64')*mass[CoM_r].astype('float64') ) / sum( mass[CoM_r].astype('float64') )
        CoM_vy = sum( vyy7[CoM_r].astype('float64')*mass[CoM_r].astype('float64') ) / sum( mass[CoM_r].astype('float64') )
        CoM_vz = sum( vzz7[CoM_r].astype('float64')*mass[CoM_r].astype('float64') ) / sum( mass[CoM_r].astype('float64') )

    print('CoM_vx, CoM_vy, CoM_vz = ', CoM_vx, CoM_vy, CoM_vz, 'km/s')

    #################################################################################################
    # Only focus on the interesting region to save time; avoid looping over all particles

    if slice_ind<4:
        slice_region_1,   = np.where( (cyl_rad <= 1.1*Rmax)   & (np.abs(zz1 - zcen) <= 1.1*Zmax) )   # Take a buffer around the boundaries
        dslice_region_1,  = np.where( (dcyl_rad <= 1.1*Rmax)  & (np.abs(dzz1 - zcen) <= 1.1*Zmax) )  # Take a buffer around the boundaries
        stslice_region_1, = np.where( (stcyl_rad <= 1.1*Rmax) & (np.abs(stzz1 - zcen) <= 1.1*Zmax) ) # Take a buffer around the boundaries
        print( 'Ncell with halos=',np.size(slice_region_1) )
        slice_region_1,   = np.where( (cyl_rad <= 1.1*Rmax)   & (np.abs(zz1 - zcen) <= 1.1*Zmax) & (mass > 0.0) ) # No buffer at all (except the centre region, which is needed for the CoM velocities)
        print( 'Ncell no halos=',  np.size(slice_region_1) )
        # Select gas cells in desired region
        xx_1   = xx1[slice_region_1].astype('float64') - xcen   # kpc
        yy_1   = yy1[slice_region_1].astype('float64') - ycen   # kpc
        zz_1   = zz1[slice_region_1].astype('float64') - zcen   # kpc
        vx_1   = ( vxx1[slice_region_1].astype('float64') - CoM_vx ) * UnitVelocity_in_cm_per_s  # cm/s
        vy_1   = ( vyy1[slice_region_1].astype('float64') - CoM_vy ) * UnitVelocity_in_cm_per_s  # cm/s
        vz_1   = ( vzz1[slice_region_1].astype('float64') - CoM_vz ) * UnitVelocity_in_cm_per_s  # cm/s
        # dark matter
        dxx_1  = dxx1[dslice_region_1] - xcen # kpc
        dyy_1  = dyy1[dslice_region_1] - ycen # kpc
        dzz_1  = dzz1[dslice_region_1] - zcen # kpc
        # stars
        stxx_1 = stxx1[stslice_region_1] - xcen # kpc
        styy_1 = styy1[stslice_region_1] - ycen # kpc
        stzz_1 = stzz1[stslice_region_1] - zcen # kpc

    elif slice_ind<6: # In slices 5 and 6, the filament axis is the y axis
        slice_region_1,   = np.where( (cyl_rad <= 1.1*Rmax)   & (np.abs(yy1 - ycen) <= 1.1*Zmax) )   # Take a buffer around the boundaries
        dslice_region_1,  = np.where( (dcyl_rad <= 1.1*Rmax)  & (np.abs(dyy1 - ycen) <= 1.1*Zmax) )  # Take a buffer around the boundaries
        stslice_region_1, = np.where( (stcyl_rad <= 1.1*Rmax) & (np.abs(styy1 - ycen) <= 1.1*Zmax) ) # Take a buffer around the boundaries
        print( 'Ncell with halos=',np.size(slice_region_1) )
        slice_region_1,   = np.where( (cyl_rad <= 1.1*Rmax)   & (np.abs(yy1 - ycen) <= 1.1*Zmax) & (mass > 0.0) ) # No buffer at all (except the centre region, which is needed for the CoM velocities)
        print( 'Ncell no halos=',  np.size(slice_region_1) )
        # Select gas cells in desired region - In slices 5 and 6, the filament axis is the y axis, so we swap the Y and Z axes when setting the frame
        xx_1   = xx1[slice_region_1].astype('float64') - xcen   # kpc
        yy_1   = zz1[slice_region_1].astype('float64') - zcen   # kpc
        zz_1   = yy1[slice_region_1].astype('float64') - ycen   # kpc
        vx_1   = ( vxx1[slice_region_1].astype('float64') - CoM_vx ) * UnitVelocity_in_cm_per_s  # cm/s
        vy_1   = ( vzz1[slice_region_1].astype('float64') - CoM_vy ) * UnitVelocity_in_cm_per_s  # cm/s
        vz_1   = ( vyy1[slice_region_1].astype('float64') - CoM_vz ) * UnitVelocity_in_cm_per_s  # cm/s
        # dark matter - In slices 5 and 6, the filament axis is the y axis, so we swap the Y and Z axes when setting the frame
        dxx_1  = dxx1[dslice_region_1] - xcen # kpc
        dyy_1  = dzz1[dslice_region_1] - zcen # kpc
        dzz_1  = dyy1[dslice_region_1] - ycen # kpc
        # stars - In slices 5 and 6, the filament axis is the y axis, so we swap the Y and Z axes when setting the frame
        stxx_1 = stxx1[stslice_region_1] - xcen # kpc
        styy_1 = stzz1[stslice_region_1] - zcen
        stzz_1 = styy1[stslice_region_1] - ycen

    elif slice_ind<10:
        slice_region_1,   = np.where( (cyl_rad <= 1.1*Rmax)   & (np.abs(zz7 - zcen) <= 1.1*Zmax) ) # No buffer at all (except the centre region, which is needed for the CoM velocities)
        dslice_region_1,  = np.where( (dcyl_rad <= 1.1*Rmax)  & (np.abs(dzz7 - zcen) <= 1.1*Zmax) ) # No buffer at all (except the centre region, which is needed for the CoM velocities)
        stslice_region_1, = np.where( (stcyl_rad <= 1.1*Rmax) & (np.abs(stzz7 - zcen) <= 1.1*Zmax) ) # Make sure the boudary's vacinities are included included
        print( 'Ncell with halos=',np.size(slice_region_1) )
        slice_region_1,   = np.where( (cyl_rad <= 1.1*Rmax)   & (np.abs(zz7 - zcen) <= 1.1*Zmax) & (mass > 0.0) ) # No buffer at all (except the centre region, which is needed for the CoM velocities)
        print( 'Ncell no halos=',  np.size(slice_region_1) )
        # Select gas cells in desired region
        xx_1   = xx7[slice_region_1].astype('float64') - xcen   # kpc
        yy_1   = yy7[slice_region_1].astype('float64') - ycen   # kpc
        zz_1   = zz7[slice_region_1].astype('float64') - zcen   # kpc
        vx_1   = ( vxx7[slice_region_1].astype('float64') - CoM_vx ) * UnitVelocity_in_cm_per_s  # cm/s
        vy_1   = ( vyy7[slice_region_1].astype('float64') - CoM_vy ) * UnitVelocity_in_cm_per_s  # cm/s
        vz_1   = ( vzz7[slice_region_1].astype('float64') - CoM_vz ) * UnitVelocity_in_cm_per_s  # cm/s
        # dark matter
        dxx_1  = dxx7[dslice_region_1] - xcen    # kpc
        dyy_1  = dyy7[dslice_region_1] - ycen    # kpc
        dzz_1  = dzz7[dslice_region_1] - zcen    # kpc
        # stars
        stxx_1 = stxx7[stslice_region_1] - xcen # stars
        styy_1 = styy7[stslice_region_1] - ycen
        stzz_1 = stzz7[stslice_region_1] - zcen

    ncell_1 = np.size(slice_region_1)
    ntrack = (ncell_1 - np.mod(ncell_1,20)) / 20
    print('ncell_1 =', ncell_1)
    print('ntrack =', ntrack)

    dncell_1 = np.size(dslice_region_1)
    dntrack = (dncell_1 - np.mod(dncell_1,20)) / 20
    print('dncell_1 =', dncell_1)
    print('dntrack =', dntrack)

    stncell_1 = np.size(stslice_region_1)
    stntrack = (stncell_1 - np.mod(stncell_1,20)) / 20
    print('stncell_1 =', stncell_1)
    print('stntrack =', stntrack)
    
    rr_1       = cyl_rad[slice_region_1].astype('float64')    # kpc
    mass_1     = mass[slice_region_1].astype('float64')       # gr
    rho_1      = rho[slice_region_1].astype('float64')        # gr/cm^3
    pressure_1 = pressure[slice_region_1].astype('float64')   # erg/cm^3
    volume_1   = volume[slice_region_1].astype('float64')     # kpc^3
    pot_1      = pot[slice_region_1].astype('float64')        # (cm/sec)^2
    Delta_cell = volume_1**(1.0/3.0)                          # kpc - effective cell size (diameter)
    #utherm_1  = ugas[slice_region_1] * mass[slice_region_1]    # erg - Total thermal energy in the cell = Pressure*Volume/(Gamma_gas - 1) = 1.5*P*V
    utherm_1   = pressure_1 * volume_1 * UnitLength_in_cm**3
    mpot_1     = pot_1 * mass_1
    vpot_1     = pot_1 * volume_1 * UnitLength_in_cm**3
    px_1 = mass_1 * vx_1    # x momentum
    py_1 = mass_1 * vy_1    # y momentum
    pz_1 = mass_1 * vz_1    # z momentum
    drr_1    = dcyl_rad[dslice_region_1].astype('float64')    # dark matter
    strr_1   = stcyl_rad[stslice_region_1].astype('float64')  # stars
    stmass_1 = stmass[stslice_region_1].astype('float64')
    print('min,max mass =', np.min(mass_1), np.max(mass_1))

    hlist_1    = h_over_Delta * Delta_cell # In kpc
    print(' ')
    print('hmin, hmax, hmean, hmean_vol', np.min(hlist_1), np.max(hlist_1), np.mean(hlist_1), np.sum(hlist_1 * volume_1)/np.sum(volume_1))
    print('min(h/Delta), max(h/Delta), mean(h/Delta)', np.min(hlist_1/Delta_cell), np.max(hlist_1/Delta_cell), np.mean(hlist_1/Delta_cell))

    pp_1 = np.zeros(ncell_1)
    for i in range(ncell_1):
        if yy_1[i] >= 0.0 and xx_1[i] >= 0.0:
            pp_1[i] =             np.arcsin( np.abs(yy_1[i]) / rr_1[i] )
        elif yy_1[i] >= 0.0 and xx_1[i] < 0.0:
            pp_1[i] = 0.5*np.pi + np.arcsin( np.abs(xx_1[i]) / rr_1[i] )
        elif yy_1[i] < 0.0 and xx_1[i] < 0.0:
            pp_1[i] =     np.pi + np.arcsin( np.abs(yy_1[i]) / rr_1[i] )
        elif yy_1[i] < 0.0 and xx_1[i] >= 0.0:
            pp_1[i] = 1.5*np.pi + np.arcsin( np.abs(xx_1[i]) / rr_1[i] )

    r_grid = np.zeros([rbin_num, pbin_num, zbin_num]) # r value for each bin
    p_grid = np.zeros([rbin_num, pbin_num, zbin_num]) # phi value for each bin
    z_grid = np.zeros([rbin_num, pbin_num, zbin_num]) # Probably not useful, but I'm doing it
    for i in range(rbin_num):
        r_grid[i, :, :] = r_bin[i]
    for j in range(pbin_num):
        p_grid[:, j, :] = p_bin[j]
    for k in range(zbin_num):
        z_grid[:, :, k] = z_bin[k]

    x_grid   = r_grid*np.cos(p_grid)
    y_grid   = r_grid*np.sin(p_grid)
    ### I'm going to try without this extra array, to save memory. Really we just need r_grid... ###
    # vol_grid = (r_grid - 0.5*Delta_r) * Delta_r * Delta_p * Delta_z # kpc^3, implies bins from [r-dr, r]
    # vol_grid = r_grid * Delta_r * Delta_p * Delta_z # kpc^3, implies bins from [r-dr/2, r+dr/2]

    #################################################################################################
    # Now we can play

    h1 = 8.0 / np.pi
    h2 = 6.0 * h1
    h3 = 2.0 * h1
    print('h1, h2, h3 = ', h1, h2, h3)

    mass_dumped_2   = np.zeros((rbin_num, pbin_num, zbin_num))
    px_dumped_2     = np.zeros((rbin_num, pbin_num, zbin_num))
    py_dumped_2     = np.zeros((rbin_num, pbin_num, zbin_num))
    pz_dumped_2     = np.zeros((rbin_num, pbin_num, zbin_num))
    utherm_dumped_2 = np.zeros((rbin_num, pbin_num, zbin_num))
    mpot_dumped_2   = np.zeros((rbin_num, pbin_num, zbin_num))
    vpot_dumped_2   = np.zeros((rbin_num, pbin_num, zbin_num))

    counter         = np.zeros(ncell_1)
    test_ms         = np.zeros(ncell_1)
    test_px         = np.zeros(ncell_1)
    test_py         = np.zeros(ncell_1)
    test_pz         = np.zeros(ncell_1)
    test_ut         = np.zeros(ncell_1)
    test_mp         = np.zeros(ncell_1)
    test_vp         = np.zeros(ncell_1)
    w_test          = np.zeros(ncell_1)
    n_cell_not_deposited = 0

    if(h_over_Delta > 0.0):
        for m in range(ncell_1): # loop over all the cells
            if( (np.mod(m,ntrack)==0) | (m==0) ):
                print(' ')
                print("m, ncell_1= ",      m, ncell_1 )

            w_tot = 0.0
            u = np.sqrt( (xx_1[m] - x_grid)**2 + (yy_1[m] - y_grid)**2 + (zz_1[m] - z_grid)**2 ) / hlist_1[m]
            # w1 = ( h3 * (1.0 - u)**3 ) * vol_grid
            # w2 = ( h1 + h2 * (u - 1.0) * u**2 ) * vol_grid
            umin = np.min(u)
            if umin>=1:  # This is very strange - a gas cell did not get deposited into any bins!
                n_cell_not_deposited += 1
                if (rr_1[m] <= Rmax) & (rr_1[m] >= Rmin) & (np.abs(zz_1[m]) <= Zmax):   # We took a buffer when selecting the cells, so now need to verify that the cell is actually in the volume
                    d_ind = np.where(u == u.min())
                    mass_dumped_2[d_ind[0][0]][d_ind[1][0]][d_ind[2][0]] += mass_1[m]
                    px_dumped_2[d_ind[0][0]][d_ind[1][0]][d_ind[2][0]] += px_1[m]
                    py_dumped_2[d_ind[0][0]][d_ind[1][0]][d_ind[2][0]] += py_1[m]
                    pz_dumped_2[d_ind[0][0]][d_ind[1][0]][d_ind[2][0]] += pz_1[m]
                    utherm_dumped_2[d_ind[0][0]][d_ind[1][0]][d_ind[2][0]] += utherm_1[m]
                    mpot_dumped_2[d_ind[0][0]][d_ind[1][0]][d_ind[2][0]]   += mpot_1[m]
                    vpot_dumped_2[d_ind[0][0]][d_ind[1][0]][d_ind[2][0]]   += vpot_1[m]

            else:
                b1 = np.where( (u >= 0.5) & (u < 1) )
                # sum1 = np.sum( w1[b1] )
                # sum1 = np.sum( ( h3 * (1.0 - u[b1])**3 ) * vol_grid[b1] )
                sum1 = np.sum( ( h3 * (1.0 - u[b1])**3 ) * r_grid[b1] ) # vol_grid = r_grid*Delta_r*Delta_p*Delta_z, but all Delta's are constant and are normalized out anyway
                w_tot += sum1
                if umin<=0.5:
                    b2 = np.where( u < 0.5 )
                    # sum2 = np.sum( w2[b2] )
                    # sum2 = np.sum( ( h1 + h2 * (u[b2] - 1.0) * u[b2]**2 ) * vol_grid[b2] )
                    sum2 = np.sum( ( h1 + h2 * (u[b2] - 1.0) * u[b2]**2 ) * r_grid[b2] )  # vol_grid = r_grid*Delta_r*Delta_p*Delta_z, but all Delta's are constant and are normalized out anyway
                    w_tot += sum2
                else:
                    sum2 = 0.0
        
                # mass_dumped_2[b1]   += mass_1[m]   * w1[b1] / w_tot
                # px_dumped_2[b1]     += px_1[m]     * w1[b1] / w_tot
                # py_dumped_2[b1]     += py_1[m]     * w1[b1] / w_tot
                # pz_dumped_2[b1]     += pz_1[m]     * w1[b1] / w_tot
                # utherm_dumped_2[b1] += utherm_1[m] * w1[b1] / w_tot
                # mpot_dumped_2[b1]   += mpot_1[m]   * w1[b1] / w_tot
                # vpot_dumped_2[b1]   += vpot_1[m]   * w1[b1] / w_tot
                #mass_dumped_2[b1]   += mass_1[m]   * ( h3 * (1.0 - u[b1])**3 ) * vol_grid[b1] / w_tot
                #px_dumped_2[b1]     += px_1[m]     * ( h3 * (1.0 - u[b1])**3 ) * vol_grid[b1] / w_tot
                #py_dumped_2[b1]     += py_1[m]     * ( h3 * (1.0 - u[b1])**3 ) * vol_grid[b1] / w_tot
                #pz_dumped_2[b1]     += pz_1[m]     * ( h3 * (1.0 - u[b1])**3 ) * vol_grid[b1] / w_tot
                #utherm_dumped_2[b1] += utherm_1[m] * ( h3 * (1.0 - u[b1])**3 ) * vol_grid[b1] / w_tot
                #mpot_dumped_2[b1]   += mpot_1[m]   * ( h3 * (1.0 - u[b1])**3 ) * vol_grid[b1] / w_tot
                #vpot_dumped_2[b1]   += vpot_1[m]   * ( h3 * (1.0 - u[b1])**3 ) * vol_grid[b1] / w_tot
                mass_dumped_2[b1]   += mass_1[m]   * ( h3 * (1.0 - u[b1])**3 ) * r_grid[b1] / w_tot
                px_dumped_2[b1]     += px_1[m]     * ( h3 * (1.0 - u[b1])**3 ) * r_grid[b1] / w_tot
                py_dumped_2[b1]     += py_1[m]     * ( h3 * (1.0 - u[b1])**3 ) * r_grid[b1] / w_tot
                pz_dumped_2[b1]     += pz_1[m]     * ( h3 * (1.0 - u[b1])**3 ) * r_grid[b1] / w_tot
                utherm_dumped_2[b1] += utherm_1[m] * ( h3 * (1.0 - u[b1])**3 ) * r_grid[b1] / w_tot
                mpot_dumped_2[b1]   += mpot_1[m]   * ( h3 * (1.0 - u[b1])**3 ) * r_grid[b1] / w_tot
                vpot_dumped_2[b1]   += vpot_1[m]   * ( h3 * (1.0 - u[b1])**3 ) * r_grid[b1] / w_tot
                counter[m]    += np.size(b1)
                w_test[m]     += sum1 / w_tot
                test_ms[m]    += mass_1[m]   * sum1 / w_tot
                test_px[m]    += px_1[m]     * sum1 / w_tot
                test_py[m]    += py_1[m]     * sum1 / w_tot
                test_pz[m]    += pz_1[m]     * sum1 / w_tot
                test_ut[m]    += utherm_1[m] * sum1 / w_tot
                test_mp[m]    += mpot_1[m]   * sum1 / w_tot
                test_vp[m]    += vpot_1[m]   * sum1 / w_tot
                if umin<=0.5:
                    # mass_dumped_2[b2]   += mass_1[m]   * w2[b2] / w_tot
                    # px_dumped_2[b2]     += px_1[m]     * w2[b2] / w_tot
                    # py_dumped_2[b2]     += py_1[m]     * w2[b2] / w_tot
                    # pz_dumped_2[b2]     += pz_1[m]     * w2[b2] / w_tot
                    # utherm_dumped_2[b2] += utherm_1[m] * w2[b2] / w_tot
                    # mpot_dumped_2[b2]   += mpot_1[m]   * w2[b2] / w_tot
                    # vpot_dumped_2[b2]   += vpot_1[m]   * w2[b2] / w_tot
                    #mass_dumped_2[b2]   += mass_1[m]   * ( h1 + h2 * (u[b2] - 1.0) * u[b2]**2 ) * vol_grid[b2] / w_tot
                    #px_dumped_2[b2]     += px_1[m]     * ( h1 + h2 * (u[b2] - 1.0) * u[b2]**2 ) * vol_grid[b2] / w_tot
                    #py_dumped_2[b2]     += py_1[m]     * ( h1 + h2 * (u[b2] - 1.0) * u[b2]**2 ) * vol_grid[b2] / w_tot
                    #pz_dumped_2[b2]     += pz_1[m]     * ( h1 + h2 * (u[b2] - 1.0) * u[b2]**2 ) * vol_grid[b2] / w_tot
                    #utherm_dumped_2[b2] += utherm_1[m] * ( h1 + h2 * (u[b2] - 1.0) * u[b2]**2 ) * vol_grid[b2] / w_tot
                    #mpot_dumped_2[b2]   += mpot_1[m]   * ( h1 + h2 * (u[b2] - 1.0) * u[b2]**2 ) * vol_grid[b2] / w_tot
                    #vpot_dumped_2[b2]   += vpot_1[m]   * ( h1 + h2 * (u[b2] - 1.0) * u[b2]**2 ) * vol_grid[b2] / w_tot
                    mass_dumped_2[b2]   += mass_1[m]   * ( h1 + h2 * (u[b2] - 1.0) * u[b2]**2 ) * r_grid[b2] / w_tot
                    px_dumped_2[b2]     += px_1[m]     * ( h1 + h2 * (u[b2] - 1.0) * u[b2]**2 ) * r_grid[b2] / w_tot
                    py_dumped_2[b2]     += py_1[m]     * ( h1 + h2 * (u[b2] - 1.0) * u[b2]**2 ) * r_grid[b2] / w_tot
                    pz_dumped_2[b2]     += pz_1[m]     * ( h1 + h2 * (u[b2] - 1.0) * u[b2]**2 ) * r_grid[b2] / w_tot
                    utherm_dumped_2[b2] += utherm_1[m] * ( h1 + h2 * (u[b2] - 1.0) * u[b2]**2 ) * r_grid[b2] / w_tot
                    mpot_dumped_2[b2]   += mpot_1[m]   * ( h1 + h2 * (u[b2] - 1.0) * u[b2]**2 ) * r_grid[b2] / w_tot
                    vpot_dumped_2[b2]   += vpot_1[m]   * ( h1 + h2 * (u[b2] - 1.0) * u[b2]**2 ) * r_grid[b2] / w_tot
                    counter[m]    += np.size(b2)
                    w_test[m]     += sum2 / w_tot
                    test_ms[m]    += mass_1[m]   * sum2 / w_tot
                    test_px[m]    += px_1[m]     * sum2 / w_tot
                    test_py[m]    += py_1[m]     * sum2 / w_tot
                    test_pz[m]    += pz_1[m]     * sum2 / w_tot
                    test_ut[m]    += utherm_1[m] * sum2 / w_tot
                    test_mp[m]    += mpot_1[m]   * sum2 / w_tot
                    test_vp[m]    += vpot_1[m]   * sum2 / w_tot
        
        b1, = np.where(counter>0)
        print('w_test=',(np.abs(w_test[b1]-1.0)).max(), (np.abs(w_test[b1]-1.0)).min()) # Just a sanity check to see if it's very small
        print('test_ms=',(np.abs(test_ms[b1]  - mass_1[b1])   / mass_1[b1]).max(),   (np.abs(test_ms[b1]  - mass_1[b1])   / mass_1[b1]).min())   # Just a sanity check to see if it's very small
        print('test_px=',(np.abs((test_px[b1] - px_1[b1])     / px_1[b1])).max(),    (np.abs((test_px[b1] - px_1[b1])     / px_1[b1])).min())   # Just a sanity check to see if it's very small
        print('test_py=',(np.abs((test_py[b1] - py_1[b1])     / py_1[b1])).max(),    (np.abs((test_py[b1] - py_1[b1])     / py_1[b1])).min())   # Just a sanity check to see if it's very small
        print('test_pz=',(np.abs((test_pz[b1] - pz_1[b1])     / pz_1[b1])).max(),    (np.abs((test_pz[b1] - pz_1[b1])     / pz_1[b1])).min())   # Just a sanity check to see if it's very small
        print('test_ut=',(np.abs(test_ut[b1]  - utherm_1[b1]) / utherm_1[b1]).max(), (np.abs(test_ut[b1]  - utherm_1[b1]) / utherm_1[b1]).min())   # Just a sanity check to see if it's very small
        print('test_mp=',(np.abs(test_mp[b1]  - mpot_1[b1])   / mpot_1[b1]).max(),   (np.abs(test_mp[b1]  - mpot_1[b1])   / mpot_1[b1]).min())   # Just a sanity check to see if it's very small
        print('test_vp=',(np.abs(test_vp[b1]  - vpot_1[b1])   / vpot_1[b1]).max(),   (np.abs(test_vp[b1]  - vpot_1[b1])   / vpot_1[b1]).min())   # Just a sanity check to see if it's very small
        print('bins per cell (>0)=',counter[b1].max(), counter[b1].min(), counter[b1].mean()) # How many bins is each cell being deposited into?
        print('The number of cells that are not deposited into any bin is', n_cell_not_deposited, np.size(np.where(counter==0)))
        print('The fraction of cells that are not depostied into any bin is', 1.0*n_cell_not_deposited / (1.0*ncell_1))

        #bin_volumes_2 = vol_grid * UnitLength_in_cm**3
        #rho_dumped_2  = mass_dumped_2   / bin_volumes_2
        #pres_dumped_2 = utherm_dumped_2 / bin_volumes_2
        # rho_dumped_2    = np.zeros((rbin_num, pbin_num, zbin_num))
        # pres_dumped_2   = np.zeros((rbin_num, pbin_num, zbin_num))
        # vx_dumped_2     = np.zeros((rbin_num, pbin_num, zbin_num))
        # vy_dumped_2     = np.zeros((rbin_num, pbin_num, zbin_num))
        # vz_dumped_2     = np.zeros((rbin_num, pbin_num, zbin_num))
        # b1 = np.where(mass_dumped_2>0.0)
        # rho_dumped_2  = mass_dumped_2   / ( r_grid * Delta_r * Delta_p * Delta_z * UnitLength_in_cm**3 )
        # pres_dumped_2 = utherm_dumped_2 / ( r_grid * Delta_r * Delta_p * Delta_z * UnitLength_in_cm**3 )
        # vx_dumped_2[b1]   = px_dumped_2[b1]     / mass_dumped_2[b1]
        # vy_dumped_2[b1]   = py_dumped_2[b1]     / mass_dumped_2[b1]
        # vz_dumped_2[b1]   = pz_dumped_2[b1]     / mass_dumped_2[b1]

        b1 = np.where(mass_dumped_2>0.0)
        px_dumped_2[b1]     = px_dumped_2[b1]     / mass_dumped_2[b1]   # velocity, cm/s
        py_dumped_2[b1]     = py_dumped_2[b1]     / mass_dumped_2[b1]   # velocity, cm/s
        pz_dumped_2[b1]     = pz_dumped_2[b1]     / mass_dumped_2[b1]   # velocity, cm/s
        mpot_dumped_2[b1]   = mpot_dumped_2[b1]   / mass_dumped_2[b1]   # mass-weighted potential, (cm/s)^2
        mass_dumped_2[b1]   = mass_dumped_2[b1]   / ( r_grid[b1] * Delta_r * Delta_p * Delta_z * UnitLength_in_cm**3 )  # density, gr/cm^3
        utherm_dumped_2[b1] = utherm_dumped_2[b1] / ( r_grid[b1] * Delta_r * Delta_p * Delta_z * UnitLength_in_cm**3 )  # pressure, erg/cm^3
        vpot_dumped_2[b1]   = vpot_dumped_2[b1]   / ( r_grid[b1] * Delta_r * Delta_p * Delta_z * UnitLength_in_cm**3 )  # volume-weighted potential, (cm/s)^2
    else:
        print('h_over_Delta=',h_over_Delta)
        print('This means we are going with nearest neighbours!')

    # Now let's take care of the empty bins
    n_empty_bin = 0
    for i in range(rbin_num): # loop over all bins
        print("(i, rbin_num)=", i, rbin_num)
        for j in range(pbin_num):
            for k in range(zbin_num):
                if( mass_dumped_2[i, j, k] == 0.0 ):
                    n_empty_bin += 1
                    distance_cells_to_bin = np.sqrt((xx_1 - x_grid[i,j,k])**2 + (yy_1 - y_grid[i,j,k])**2 + (zz_1 - z_grid[i,j,k])**2)
                    closest_cell = np.where(distance_cells_to_bin == distance_cells_to_bin.min())
                    # rho_dumped_2[i, j, k]  = rho_1[closest_cell]
                    # vx_dumped_2[i, j, k]   = vx_1[closest_cell]
                    # vy_dumped_2[i, j, k]   = vy_1[closest_cell]
                    # vz_dumped_2[i, j, k]   = vz_1[closest_cell]
                    # pres_dumped_2[i, j, k] = pressure_1[closest_cell]
                    mass_dumped_2[i, j, k]   = rho_1[closest_cell]
                    px_dumped_2[i, j, k]     = vx_1[closest_cell]
                    py_dumped_2[i, j, k]     = vy_1[closest_cell]
                    pz_dumped_2[i, j, k]     = vz_1[closest_cell]
                    utherm_dumped_2[i, j, k] = pressure_1[closest_cell]
                    mpot_dumped_2[i, j, k]   = pot_1[closest_cell]
                    vpot_dumped_2[i, j, k]   = pot_1[closest_cell]

    #vr_dumped_2   = vx_dumped_2 * np.cos(p_grid) + vy_dumped_2 * np.sin(p_grid)
    #vphi_dumped_2 = vy_dumped_2 * np.cos(p_grid) - vx_dumped_2 * np.sin(p_grid) # Should be already in cm/s
    vr_dumped_2   = px_dumped_2 * np.cos(p_grid) + py_dumped_2 * np.sin(p_grid)
    vphi_dumped_2 = py_dumped_2 * np.cos(p_grid) - px_dumped_2 * np.sin(p_grid) # Should be already in cm/s
    print('The number and fraction of empty bins is ', n_empty_bin, 1.0*n_empty_bin / (1.0*Nbin))

    ### Let's take the total mass profile directly from the data without the SPH ###
    ### This is because we can't define self-consistent smoothing lengths for the dark matter and stellar particles like for the gas cells ###
    Mtot_bin_2 = np.zeros(rbin_num)
    Mgas_bin_2 = np.zeros(rbin_num)
    for i in range(rbin_num): # loop over all bins
        b,  = np.where( (rr_1 <= r_bin[i]) & (np.abs(zz_1) <= L_slice) )
        Mtot_bin_2[i] += np.sum(mass_1[b])
        Mgas_bin_2[i] += np.sum(mass_1[b])

        db, = np.where( (drr_1 <= r_bin[i]) & (np.abs(dzz_1) <= L_slice) )
        Mtot_bin_2[i] += dmass * np.size(db)

        stb, = np.where( (strr_1 <= r_bin[i]) & (np.abs(stzz_1) <= L_slice) )
        Mtot_bin_2[i] += np.sum(stmass_1[stb])

    #################################################################################################
    Dir = FigDir+"/slice_%03d/" % (slice_ind+1)
    if not os.path.exists(Dir):
        os.mkdir(Dir)
    #np.savetxt(Dir + 'rho_%03d.txt'  % (slice_ind+1), np.reshape(rho_dumped_2, Nbin))
    #np.savetxt(Dir + 'pres_%03d.txt' % (slice_ind+1), np.reshape(pres_dumped_2, Nbin))
    #np.savetxt(Dir + 'vr_%03d.txt'   % (slice_ind+1), np.reshape(vr_dumped_2, Nbin))
    #np.savetxt(Dir + 'vphi_%03d.txt' % (slice_ind+1), np.reshape(vphi_dumped_2, Nbin))
    #np.savetxt(Dir + 'vz_%03d.txt'   % (slice_ind+1), np.reshape(vz_dumped_2, Nbin))
    np.savetxt(Dir + 'rho_%03d.txt'  % (slice_ind+1), np.reshape(mass_dumped_2, Nbin))
    np.savetxt(Dir + 'pres_%03d.txt' % (slice_ind+1), np.reshape(utherm_dumped_2, Nbin))
    np.savetxt(Dir + 'mpot_%03d.txt' % (slice_ind+1), np.reshape(mpot_dumped_2, Nbin))
    np.savetxt(Dir + 'vpot_%03d.txt' % (slice_ind+1), np.reshape(vpot_dumped_2, Nbin))
    np.savetxt(Dir + 'vr_%03d.txt'   % (slice_ind+1), np.reshape(vr_dumped_2, Nbin))
    np.savetxt(Dir + 'vphi_%03d.txt' % (slice_ind+1), np.reshape(vphi_dumped_2, Nbin))
    np.savetxt(Dir + 'vz_%03d.txt'   % (slice_ind+1), np.reshape(pz_dumped_2, Nbin))
    np.savetxt(Dir + 'Mtot_%03d.txt' % (slice_ind+1), Mtot_bin_2)
    np.savetxt(Dir + 'Mgas_%03d.txt' % (slice_ind+1), Mgas_bin_2)