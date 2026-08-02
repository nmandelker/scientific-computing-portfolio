from gadget import *
from gadget_subfind import *
import const as c
import matplotlib as mpl
import os, sys

######################################################

ZF = 4

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

UnitNumberDensity_in_cgs = UnitDensity_in_cgs/mproton	# ntot in cm^{-3}

######################################################

xsize   = 5.3
ysize   = 5.3
fsize   = 17
fsize2  = 19
fweight = 150
ax1_pos  = [0.16, 0.12, 0.80, 0.80]

xstr = r"$a$"
ystr = r"${\rm log}(N_{\rm cell})$"

xlaby = -0.08
ylabx = -0.12
tit_pos = [0.5, 1.005]

xtick = 0.1 + np.arange(0,11)*0.02
ytick = 7.5 + np.arange(0,7)*0.5

xminorLocator = MultipleLocator(0.01)
yminorLocator = MultipleLocator(0.1)

xlim = [0.1, 0.28]
ylim = [7.5, 10.5]

######################################################

FigDir   = "../plots/FilB-Zoom4"
if not os.path.exists(FigDir):
        os.mkdir(FigDir)
FigDir   = "%s/for_letter" % (FigDir)
if not os.path.exists(FigDir):
        os.mkdir(FigDir)

SimLabel = 'FilB-Zoom%01d' % (ZF)
in_filename = "halo-positions-%s.dat" % (SimLabel)
SimDir   = "/home/mandelnr/basement/%s/output" % (SimLabel)
SnapBase = 'snap_'

Ncell_list = np.zeros([32,3])
aexp_list  = np.zeros([32])
for snap_ind in np.arange(0,32):
	print('snap_ind ',snap_ind)

	s  = gadget_readsnap( snap_ind, SimDir, SnapBase, loadonlytype=[0], loadonly=['pos', 'mass', 'vol'], lazy_load=True )
	sf = load_subfind( snap_ind, dir=SimDir+'/' )
	redshift = s.redshift
	BoxHalf = 0.5 * s.boxsize
	aexp = s.time
	h0 = s.hubbleparam
	BoxHalf *= aexp/h0
	BoxSize = 2.0*BoxHalf
	Ng = sf.totngroups
	Ncell_halos = np.zeros([Ng])

	aexp_list[snap_ind] = aexp
	Ncell_list[snap_ind,0] = np.size(s.data['vol'])

	print('Ncell ',Ncell_list[snap_ind,0])
	print('Ngroup ',Ng)

	for i in range(Ng):
		if sf.data['fmt2'][i]>0.01:
			Ncell_halos[i] = sf.data['flen'][i]
			Rv = sf.data['frt2'][i]
			Ncell_list[snap_ind,1] += (4.0*np.pi/3.0)*(Rv/0.2)**3
			Ncell_list[snap_ind,2] += (4.0*np.pi/3.0)*(Rv/0.1)**3

	Ncell_IPM = Ncell_list[snap_ind,0] - np.sum(Ncell_halos)
	Ncell_list[snap_ind,1] += Ncell_IPM
	Ncell_list[snap_ind,2] += Ncell_IPM

######################################################

### Open figure object ###
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig = plt.figure(figsize = (xsize,ysize))
fig_tit = "%s/Ncell_estimate.png" % (FigDir)

######################################################

ax1 = plt.axes(ax1_pos)
line1, = plt.plot( aexp_list, np.log10(Ncell_list[:,0]), marker='o', markerfacecolor='k', markeredgecolor='k', markersize=4, linestyle='-', color='k', linewidth=2 )
line2, = plt.plot( aexp_list, np.log10(Ncell_list[:,1]), marker='o', markerfacecolor='r', markeredgecolor='r', markersize=4, linestyle='-', color='r', linewidth=2 )
line3, = plt.plot( aexp_list, np.log10(Ncell_list[:,2]), marker='o', markerfacecolor='g', markeredgecolor='g', markersize=4, linestyle='-', color='g', linewidth=2 )

xlabel(xstr, fontname='times', fontsize = fsize2, fontweight = fweight)
ylabel(ystr, fontname='times', fontsize = fsize2, fontweight = fweight)
ax1.xaxis.set_label_coords(0.5, xlaby)
ax1.yaxis.set_label_coords(ylabx, 0.5)

ttl = ax1.set_title(r"${\rm ZF}"+format(ZF)+"$", fontname='times', fontsize = fsize2, fontweight = fweight, loc='center')
ttl.set_position(tit_pos)

plt.xticks(xtick, fontname='arial', fontsize = fsize, fontweight = fweight)
plt.yticks(ytick, fontname='arial', fontsize = fsize, fontweight = fweight)
ax1.xaxis.set_tick_params(which = 'major', direction = 'in', width = 1, length = 6, pad = 5, bottom = True, top = True)
ax1.yaxis.set_tick_params(which = 'major', direction = 'in', width = 1, length = 6, pad = 5, left = True, right = True)

ax1.xaxis.set_minor_locator(xminorLocator)
ax1.xaxis.set_tick_params(which = 'minor', direction = 'in', width = 1, length = 3, pad = 5, bottom = True, top = True)

ax1.yaxis.set_minor_locator(yminorLocator)
ax1.yaxis.set_tick_params(which = 'minor', direction = 'in', width = 1, length = 3, pad = 5, left = True, right = True)

disp1 = r"${\rm mass\:\:refinement}$"
disp2 = r"${\rm 200pc}$"
disp3 = r"${\rm 100pc}$"
plt.legend([line1, line2, line3],[disp1, disp2, disp3], fontsize = 14, frameon=False, loc='upper left')

ax1.set_xlim(xlim)
ax1.set_ylim(ylim)

savefig(fig_tit)
print(fig_tit)
plt.close()


