t = 0:(1/5):3;
gr_ptv = []; gr_p = []; gr_vz = [];
levelmax = [12,13,14];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir = 'schemes and convergence/kh_l01d1m05_llf1_l12';
nSnapshotMax=16; N=4096; zoomx=0.25; zoomy=0.25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%theoretical growth rate, KHI time
[omega, kzplus, kzminus, rhominus, rhoplus, V] = calc_dispersion_ramses(0.1, 1, 0.5, 5/3);
theoretical_growth_rate = imag(omega);
t_kh = 1/theoretical_growth_rate;
t = t*t_kh;

%get numerical results for color, pressure, and velocity
isoval=0.5;
ptv = color_ptv_directory(dir,nSnapshotMax,N,zoomx,zoomy,isoval);
[p_span, p_abs_mean] = variable_directory(dir,nSnapshotMax,N,zoomx,zoomy,'p');
[vz_span, vz_abs_mean] = variable_directory(dir,nSnapshotMax,N,zoomx,zoomy,'vz');
[t_ptv,growth_rate_ptv] = momemntary_growth_rate(t,ptv);
[t_p,growth_rate_p] = momemntary_growth_rate(t,p_abs_mean);
[t_vz,growth_rate_vz] = momemntary_growth_rate(t,vz_abs_mean);

%get the t->0 growth rates
tmp = polyfit(t_ptv(1:5),growth_rate_ptv(1:5),1);
gr_ptv = [gr_ptv,tmp(2)];
tmp = polyfit(t_p(1:5),growth_rate_p(1:5),1);
gr_p = [gr_p,tmp(2)];
tmp = polyfit(t_vz(1:5),growth_rate_vz(1:5),1);
gr_vz = [gr_vz,tmp(2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir = 'schemes and convergence/kh_l01d1m05_llf1';
nSnapshotMax=16; N=8192; zoomx=0.25; zoomy=0.25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get numerical results for color, pressure, and velocity
isoval=0.5;
ptv = color_ptv_directory(dir,nSnapshotMax,N,zoomx,zoomy,isoval);
[p_span, p_abs_mean] = variable_directory(dir,nSnapshotMax,N,zoomx,zoomy,'p');
[vz_span, vz_abs_mean] = variable_directory(dir,nSnapshotMax,N,zoomx,zoomy,'vz');
[t_ptv,growth_rate_ptv] = momemntary_growth_rate(t,ptv);
[t_p,growth_rate_p] = momemntary_growth_rate(t,p_abs_mean);
[t_vz,growth_rate_vz] = momemntary_growth_rate(t,vz_abs_mean);

%get the t->0 growth rates
tmp = polyfit(t_ptv(1:5),growth_rate_ptv(1:5),1);
gr_ptv = [gr_ptv,tmp(2)];
tmp = polyfit(t_p(1:5),growth_rate_p(1:5),1);
gr_p = [gr_p,tmp(2)];
tmp = polyfit(t_vz(1:5),growth_rate_vz(1:5),1);
gr_vz = [gr_vz,tmp(2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir = 'schemes and convergence/kh_l01d1m05_llf1_l14';
nSnapshotMax=16; N=8192*2; zoomx=0.25; zoomy=0.25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get numerical results for color, pressure, and velocity
isoval=0.5;
ptv = color_ptv_directory(dir,nSnapshotMax,N,zoomx,zoomy,isoval);
[p_span, p_abs_mean] = variable_directory(dir,nSnapshotMax,N,zoomx,zoomy,'p');
[vz_span, vz_abs_mean] = variable_directory(dir,nSnapshotMax,N,zoomx,zoomy,'vz');
[t_ptv,growth_rate_ptv] = momemntary_growth_rate(t,ptv);
[t_p,growth_rate_p] = momemntary_growth_rate(t,p_abs_mean);
[t_vz,growth_rate_vz] = momemntary_growth_rate(t,vz_abs_mean);

%get the t->0 growth rates
tmp = polyfit(t_ptv(1:5),growth_rate_ptv(1:5),1);
gr_ptv = [gr_ptv,tmp(2)];
tmp = polyfit(t_p(1:5),growth_rate_p(1:5),1);
gr_p = [gr_p,tmp(2)];
tmp = polyfit(t_vz(1:5),growth_rate_vz(1:5),1);
gr_vz = [gr_vz,tmp(2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%limits and ticks
xmin = [12 12]; xmax = [14 14];
ymin = [0.8 0.8]; ymax = [1 1];
xtick = 1; ytick = 0.02;

%labels and fonts
xlbl = '$l_{max}$'; ylbl = 'Instantaneous Growth Rate at $t \to 0$ [$1/t_kh$]'; ttl = ''; %ttl = 'Nondimensional Growth Rate'; 
axfnt = 14; lblfnt = 14; ttlfnt = 14;
fnt = 'Times New Roman';

%make figure
myfig(xmin,xmax,ymin,ymax,xtick,ytick,xlbl,ylbl,ttl,axfnt,lblfnt,ttlfnt,1,0);

%plot
plot(levelmax,gr_ptv*t_kh,'-or','LineWidth',2); hold on;
plot(levelmax,gr_vz*t_kh,'-ob','LineWidth',2); hold on;
plot(levelmax,gr_p*t_kh,'-og','LineWidth',2);

set(gcf,'renderer','zbuffer');
print(gcf, '-depsc', 'convergence_resolution');