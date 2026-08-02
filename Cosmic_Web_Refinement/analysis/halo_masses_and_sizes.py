from gadget import *
from gadget_subfind import *
import const as c
import matplotlib as mpl
import os, sys

SimLabel = 'FilB-Zoom2'
SnapNums = np.arange(78,126)

in_filename = "halo-positions-%s.dat" % (SimLabel)

RevSnap  = SnapNums[::-1]
SimDir   = "/home/mandelnr/basement/%s/output" % (SimLabel)
SnapBase = 'snap_'
FigDir   = "../plots/%s" % (SimLabel)
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

aexp_vec = 0.0*RevSnap
M1_vec   = 0.0*RevSnap
M2_vec   = 0.0*RevSnap
R1_vec   = 0.0*RevSnap
R2_vec   = 0.0*RevSnap
D_vec    = 0.0*RevSnap

for idx, SnapNum in enumerate(RevSnap):
	print(idx,SnapNum)
	sf = load_subfind( SnapNum, dir=SimDir+'/' )
	aexp = sf.time

	b1, = np.where(RevSnap2-SnapNum==0)
	if b1.size>0:
		if(np.abs(aexp_vec2[b1]-aexp)>0.0001):
			print(SnapNum)
			print(aexp)
			print(aexp_vec2[b1])
			sys.exit("Error: Snapshot aexp Mismatch")
	else:
		print(SnapNum)
		sys.exit("Error: Current snapshot is not in file")

	print(SnapNum,RevSnap2[b1])
	x_first = x1_vec[b1]
	y_first = y1_vec[b1]
	z_first = z1_vec[b1]
	x_second = x2_vec[b1]
	y_second = y2_vec[b1]
	z_second = z2_vec[b1]
	ind_first  = int(ind1_vec[b1])
	ind_second = int(ind2_vec[b1])
	print(ind_first, ind_second)
	print(x_first, y_first, z_first)
	print(x_second, y_second, z_second)

	aexp_vec[idx] = aexp
	M1_vec[idx] = sf.data['fmt2'][ind_first]
	M2_vec[idx] = sf.data['fmt2'][ind_second]
	R1_vec[idx] = sf.data['frt2'][ind_first]
	R2_vec[idx] = sf.data['frt2'][ind_second]
	D_vec[idx]  = np.sqrt( (x_first - x_second)**2  + (y_first - y_second)**2  + (z_first - z_second)**2 ) 


print(aexp_vec)
print(M1_vec)
print(M2_vec)

#xsize = 9.0
#ysize = 5.2857
xsize = 5.7
ysize = 5.3
fsize = 18
fsize2 = 22
fweight = 150

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig = plt.figure(figsize = (xsize,ysize))
ax1 = plt.axes([0.18, 0.1, 0.77175, 0.83])
line1, = plt.plot( aexp_vec, np.log10(M1_vec)+10, marker='None', linestyle='-', color='k', linewidth=2 )#, label=r"${\rm M_1}" )
line2, = plt.plot( aexp_vec, np.log10(M2_vec)+10, marker='None', linestyle='-', color='r', linewidth=2 )#, label=r"${\rm M_2}" )

xlabel(r"$a$", fontname='times', fontsize = fsize2, fontweight = fweight)
ylabel(r"${\rm log}(M_{\rm v}\:{\rm [M_{\odot}]})$", fontname='times', fontsize = fsize2, fontweight = fweight)
plt.xticks([0.10, 0.15, 0.20, 0.25, 0.30, 0.35], fontname='arial', fontsize = fsize, fontweight = fweight) #[' ', u'\u2212''1000', ' ', 0, ' ', '1000', ' ']
plt.yticks([10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0], fontname='arial', fontsize = fsize, fontweight = fweight)

ax1.xaxis.set_tick_params(which = 'major', direction = 'in', width = 1, length = 6, pad = 5, bottom = True, top = True)
ax1.yaxis.set_tick_params(which = 'major', direction = 'in', width = 1, length = 6, pad = 5, left = True, right = True)
minorLocator = MultipleLocator(0.01)
ax1.xaxis.set_minor_locator(minorLocator)
ax1.xaxis.set_tick_params(which = 'minor', direction = 'in', width = 1, length = 3, pad = 5, bottom = True, top = True)
minorLocator = MultipleLocator(0.1)
ax1.yaxis.set_minor_locator(minorLocator)
ax1.yaxis.set_tick_params(which = 'minor', direction = 'in', width = 1, length = 3, pad = 5, left = True, right = True)

ax1.xaxis.set_label_coords(0.5, -0.06)
ax1.yaxis.set_label_coords(-0.12, 0.5)

ax1.set_title(r"${\rm " + SimLabel + "}$", fontname='times', fontsize = fsize2, fontweight = fweight, loc='center')

#plt.legend(handles=[line1, line2])
plt.legend([line1, line2],[r"${\rm Halo\:1}$",r"${\rm Halo\:2}$"], fontsize = fsize, frameon=False, loc='best')

ax1.set_xlim([0.11,0.35])
ax1.set_ylim([10.0,13.0])
fname = "%s/halo_masses.png" % (FigDir)
savefig(fname)
print(fname)
plt.close(fig)

#######
fig = plt.figure(figsize = (xsize,ysize))
ax1 = plt.axes([0.18, 0.1, 0.77175, 0.83])
line1, = plt.plot( aexp_vec, np.log10(R1_vec), marker='None', linestyle='-', color='k', linewidth=2 )
line2, = plt.plot( aexp_vec, np.log10(R2_vec), marker='None', linestyle='-', color='r', linewidth=2 )
line3, = plt.plot( aexp_vec, np.log10(D_vec),  marker='None', linestyle='-', color='b', linewidth=2 )

xlabel(r"$a$", fontname='times', fontsize = fsize2, fontweight = fweight)
ylabel(r"${\rm log}(R_{\rm v},~D\:{\rm [kpc]})$", fontname='times', fontsize = fsize2, fontweight = fweight)
plt.xticks([0.10, 0.15, 0.20, 0.25, 0.30, 0.35], fontname='arial', fontsize = fsize, fontweight = fweight) #[' ', u'\u2212''1000', ' ', 0, ' ', '1000', ' ']
plt.yticks([1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2], fontname='arial', fontsize = fsize, fontweight = fweight)

ax1.xaxis.set_tick_params(which = 'major', direction = 'in', width = 1, length = 6, pad = 5, bottom = True, top = True)
ax1.yaxis.set_tick_params(which = 'major', direction = 'in', width = 1, length = 6, pad = 5, left = True, right = True)
minorLocator = MultipleLocator(0.01)
ax1.xaxis.set_minor_locator(minorLocator)
ax1.xaxis.set_tick_params(which = 'minor', direction = 'in', width = 1, length = 3, pad = 5, bottom = True, top = True)
minorLocator = MultipleLocator(0.05)
ax1.yaxis.set_minor_locator(minorLocator)
ax1.yaxis.set_tick_params(which = 'minor', direction = 'in', width = 1, length = 3, pad = 5, left = True, right = True)

ax1.xaxis.set_label_coords(0.5, -0.06)
ax1.yaxis.set_label_coords(-0.12, 0.5)

ax1.set_title(r"${\rm " + SimLabel + "}$", fontname='times', fontsize = fsize2, fontweight = fweight, loc='center')

#plt.legend(handles=[line1, line2])
plt.legend([line1, line2, line3],[r"${\rm Halo\:1}$",r"${\rm Halo\:2}$",r"${\rm Distance}$"], fontsize = fsize, frameon=False, loc='best')

ax1.set_xlim([0.11,0.35])
ax1.set_ylim([1.0,3.2])
fname = "%s/halo_sizes.png" % (FigDir)
savefig(fname)
print(fname)
plt.close(fig)

out_filename = "halo_masses_and_sizes_%s.dat" % (SimLabel)
np.savetxt(out_filename, (aexp_vec, np.log10(M1_vec)+10, np.log10(M2_vec)+10, np.log10(R1_vec), np.log10(R2_vec), np.log10(D_vec)))
