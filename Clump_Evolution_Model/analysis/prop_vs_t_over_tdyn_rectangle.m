function prop_vs_t_over_tdyn_rectangle(is, norm_is, nis, ind, gal, prop, tmax_window, Mmax_window, SFR_window, ...
    zform_window, log_smoothing, tform, tend, versus_dist, ...
    alpha, eta1, eps, mu, etas, tmig_td, fg_i)

% gal is the simulation index: V**
% gal==0 means stack all simulations together

%(prop = 1: ID)
%prop = 2: redshift
%prop = 3: Rc
%prop = 4: M_{gas}
%prop = 5: M_{star}
%prop = 6: M_{bar}
%prop = 7: f_{gas}
%prop = 8: f_{dm}
%prop = 9: {\Sigma}_{gas}
%prop = 10: {\Sigma}_{*}
%prop = 11: {\Sigma}_{bar}
%prop = 12: age
%prop = 13: z_{gas}
%prop = 14: z_{stars}
%prop = 15: SFR
%prop = 16: {\Sigma}_{SFR}
%prop = 17: sSFR
%prop = 18: t_{dep}
%prop = 19: d/Rd
%prop = 20: z/Hd
%prop = 21: d
%prop = 22: z
%prop = 23: mean_residual
%prop = 24: shape parameter
%prop = 25: dark matter contrast
%(prop = 26: Ex situ)
%prop = 28: tff
%prop = 29: tdyn_local
%prop = 30: tdyn_global
%prop = 48: time

dirname = './clump_evolution/';
mkdir(dirname);

Tdown = tmax_window(1);
Tup   = tmax_window(2);

Mdown = Mmax_window(1);
Mup   = Mmax_window(2);

SFRdown = SFR_window(1);
SFRup   = SFR_window(2);

zdown = zform_window(1);
zup   = zform_window(2);

if(gal~=0)
    dirname = strcat(dirname,'V',num2str(gal,'%02i'),'/');
    tit = strcat('$V',num2str(gal,'%02i'),'$');
    mkdir(dirname);
else
    dirname = strcat(dirname,'all/');
    tit = 'all';
    mkdir(dirname);
end
dirname = strcat(dirname,num2str(Mdown,'%3.1f'),'_M50_',num2str(Mup,'%3.1f'),'/');
dirname = strcat(dirname,num2str(zdown,'%3.1f'),'_zform_',num2str(zup,'%3.1f'),'/');
mkdir(dirname);

if(prop==3)
    is(:,3) = norm_is(:,3);
    log_prop = 0;
    norm_stack = 0;
    yl = 0.01;
    yu = 1;
    dy = 0.1;
    ytick = [0.01 0.1 1];
    ystr = {'0.01','0.1','1'};
    normyl = 0.1;
    normyu = 10;
    normytick = [0.1 1 10];
    normystr = {num2str(0.1,'%3.1f'),num2str(1,'%1.0f'),num2str(10,'%2.0f')};
    dirname = strcat(dirname,'Rc_over_Rd/');
    fname = 'norm_Rc';
    ylab = '$R_{\rm c}\:/\:R_{\rm d}$';
    ylab2 = '$R_{\rm c}\:\:/\:\:R_{\rm d}$';
    norm_ylab = ylab;
    norm_ylab2 = ylab2;
elseif(prop==4)
    log_prop = 0;
    norm_stack = 1;
    yl = 1e6;
    yu = 1e9;
    ytick = [1e6 1e7 1e8 1e9];
    ystr = {'1e6','1e7','1e8','1e9'};
    normyl = 0.02;
    normyu = 2;
    normytick = [0.02 0.1 1 2];
    normystr = {num2str(0.02,'%4.2f'),num2str(0.1,'%3.1f'),num2str(1,'%1.0f'),num2str(2,'%1.0f')};
    dirname = strcat(dirname,'gas_mass/');
    fname = 'Mgas';
    ylab = '$M_{\rm g}\:{\rm [M_{\odot}]}$';
    ylab2 = '$M_{\rm g}\:\:{\rm [M_{\odot}]}$';
    norm_ylab = '$M_{\rm g}\:/\:M(t_{\rm n})$';
    norm_ylab2 = '$M_{\rm g}\:\:/\:\:M(t_{\rm n})$';
elseif(prop==5)
    log_prop = 0;
    norm_stack = 1;
    yl = 1e6;
    yu = 1e9;
    ytick = [1e6 1e7 1e8 1e9];
    ystr = {'1e6','1e7','1e8','1e9'};
    normyl = 0.1;
    normyu = 10;
    normytick = [0.1 1 10];
    normystr = {num2str(0.1,'%3.1f'),num2str(1,'%1.0f'),num2str(10,'%2.0f')};
    dirname = strcat(dirname,'stellar_mass/');
    fname = 'Mstar';
    ylab = '$M_{\rm s}\:{\rm [M_{\odot}]}$';
    ylab2 = '$M_{\rm s}\:\:{\rm [M_{\odot}]}$';
    norm_ylab = '$M_{\rm s}\:/\:M(t_{\rm n})$';
    norm_ylab2 = '$M_{\rm s}\:\:/\:\:M(t_{\rm n})$';
elseif(prop==6)
    log_prop = 0;
    norm_stack = 1;
    yl = 1e7;
    yu = 1e9;
    ytick = [1e7 1e8 1e9];
    ystr = {'1e7','1e8','1e9'};
    normyl = 0.1;
    normyu = 10;
    normytick = [0.1 1 10];
    normystr = {num2str(0.1,'%3.1f'),num2str(1,'%1.0f'),num2str(10,'%2.0f')};
    dirname = strcat(dirname,'mass/');
    fname = 'mass';
    ylab = '$M\:{\rm [M_{\odot}]}$';
    ylab2 = '$M\:\:{\rm [M_{\odot}]}$';
    norm_ylab = '$M\:/\:M(t_{\rm n})$';
    norm_ylab2 = '$M\:\:/\:\:M(t_{\rm n})$';
elseif(prop==7)
    log_prop = 0;
    norm_stack = 0;
    yl = 0.01;
    yu = 1;
    dy = 0.1;
    ytick = [0.01 0.1 1];
    ystr = {num2str(0.01,'%4.2f'),num2str(0.1,'%3.1f'),num2str(1,'%1.0f')};
    dirname = strcat(dirname,'fgas/');
    fname = 'fgas';
    ylab = '$f_{\rm g}$';
    ylab2 = '$f_{\rm g}$';
    norm_ylab = '$f_{\rm g}\:/\:f(g,\,t_{\rm n})$';
    norm_ylab2 = '$f_{\rm g}\:\:/\:\:f(g,\,t_{\rm n})$';
elseif(prop==12)
    log_prop = 0;
    norm_stack = 1;
    yl = 10;
    yu = 2000;
    dy = 200;
    ytick = [10 100 1000 2000];
    ystr = {num2str(10,'%2.0f'),num2str(100,'%3.0f'),num2str(1000,'%4.0f'),num2str(2000,'%4.0f')};
    normyl = 0.3;
    normyu = 30;
    normytick = [0.3 1 10 30];
    normystr = {num2str(0.3,'%3.1f'),num2str(1,'%1.0f'),num2str(10,'%2.0f'),num2str(30,'%2.0f')};
    dirname = strcat(dirname,'age/');
    fname = 'age';
    ylab = '${\rm age_*}\:{\rm [Myr]}$';
    ylab2 = '${\rm age_*}\:\:{\rm [Myr\:]}$';
    norm_ylab = '${\rm age_*}\:/\:t_{\rm d}$';
    norm_ylab2 = '${\rm age_*}\:\:/\:\:t_{\rm d}$';
elseif(prop==15)
    log_prop = 0;
    norm_stack = 1;
    yl = 0.1;
    yu = 10;
    ytick = [0.1 1 10];
    ystr = {'0.1','1','10'};
    if(versus_dist==0)
        normyl = 0.001;
        normyu = 1;
        normytick = [0.001 0.01 0.1 1];
        normystr = {num2str(0.001,'%5.3f'),num2str(0.01,'%4.2f'),num2str(0.1,'%3.1f'),num2str(1,'%1.0f')};
    elseif(versus_dist==1)
        normyl = 0.002;
        normyu = 0.2;
        normytick = [0.002 0.01 0.1 0.2];
        normystr = {num2str(0.002,'%5.3f'),num2str(0.01,'%4.2f'),num2str(0.1,'%3.1f'),num2str(0.2,'%3.1f')};
    end
    dirname = strcat(dirname,'sfr/');
    fname = 'sfr';
    ylab = '$SFR\:{\rm [M_{\odot}\:yr^{-1}]}$';
    ylab2 = '$SFR\:\:{\rm [M_{\odot}\:yr^{-1}\:]}$';
    norm_ylab = '$SFR\:/\:[M(t_{\rm n})\:/\:t_{\rm d}]$';
    norm_ylab2 = '$SFR\:\:/\:\:[M(t_{\rm n})\:\:/\:\:t_{\rm d}]$';
elseif(prop==17)
    log_prop = 1;
    norm_stack = 1;
    yl = -1;
    yu = 1.4;
    dy = 0.2;
    normyl = -1;
    normyu = 1;
    normdy = 0.2;
    dirname = strcat(dirname,'sSFR/');
    fname = 'ssfr';
    ylab = '${\rm log(}sSFR\:{\rm [Gyr^{-1}])}$';
    ylab2 = '${\rm log(\:\:}sSFR\:\:{\rm [Gyr^{-1}\:])}$';
    norm_ylab = '${\rm log(}sSFR\:{\rm [normalized])}$';
    norm_ylab2 = '${\rm log(\:\:}sSFR\:\:{\rm [normalized])}$';
elseif(prop==177)
    is(:,177) = 1000./is(:,17);
    log_prop = 0;
    norm_stack = 1;
    yl = 10;
    yu = 2000;
    ytick = [10 100 1000 2000];
    ystr = {num2str(10,'%2.0f'),num2str(100,'%3.0f'),num2str(1000,'%4.0f'),num2str(2000,'%4.0f')};
    normyl = 0.1;
    normyu = 100;
    normytick = [0.1 1 10 100];
    normystr = {num2str(0.1,'%3.1f'),num2str(1,'%1.0f'),num2str(10,'%2.0f'),num2str(100,'%2.0f')};
    dirname = strcat(dirname,'tsf/');
    fname = 'tsf';
    ylab = '$t_{\rm sf}\:{\rm [Myr]}$';
    ylab2 = '$t_{\rm sf}\:\:{\rm [Myr\:]}$';
    norm_ylab = '$t_{\rm sf}\:/\:t_{\rm d}$';
    norm_ylab2 = '$t_{\rm sf}\:\:/\:\:t_{\rm d}$';
elseif(prop==18)
    is(:,18) = 1000.*is(:,18);
    log_prop = 0;
    norm_stack = 1;
    yl = 10;
    yu = 2000;
    ytick = [10 100 1000 2000];
    ystr = {num2str(10,'%2.0f'),num2str(100,'%3.0f'),num2str(1000,'%4.0f'),num2str(2000,'%4.0f')};
    normyl = 0.3;
    normyu = 30;
    normytick = [0.3 1 10 30];
    normystr = {num2str(0.3,'%3.1f'),num2str(1,'%1.0f'),num2str(10,'%2.0f'),num2str(30,'%2.0f')};
    dirname = strcat(dirname,'tdep/');
    fname = 'tdep';
    ylab = '$t_{\rm dep}\:{\rm [Myr]}$';
    ylab2 = '$t_{\rm dep}\:\:{\rm [Myr\:]}$';
    norm_ylab = '$t_{\rm dep}\:/\:t_{\rm d}$';
    norm_ylab2 = '$t_{\rm dep}\:\:/\:\:t_{\rm d}$';
elseif(prop==19)
    log_prop = 0;
    norm_stack = 0;
    yl = 0.1;
    yu = 2;
    dy = 0.2;
    ytick = [0.1 1 2];
    ystr = {num2str(0.1,'%3.1f'),num2str(1,'%1.0f'),num2str(2,'%1.0f')};
    dirname = strcat(dirname,'dist/');
    fname = 'dist';
    ylab = '$d\:/\:R_{\rm d}$';
    ylab2 = '$d\:\:/\:\:R_{\rm d}$';
    norm_ylab = ylab;
    norm_ylab2 = ylab2;
elseif(prop==21)
    log_prop = 1;
    norm_stack = 0;
    yl = -1;
    yu = 3;
%     yl = -0.4;
%     yu = 1.6;
    dy = 0.5;
    normyl = -2;
    normyu = 2;
%     normyl = -1;
%     normyu = 0.4;
    normdy = 0.5;
    dirname = strcat(dirname,'dist/');
    fname = 'dist';
    ylab = '${\rm log(}d\:{\rm [kpc])}$';
    ylab2 = '${\rm log(\:\:}d\:\:{\rm [kpc\:])}$';
    norm_ylab = '${\rm log(}d\:{\rm [normalized])}$';
    norm_ylab2 = '${\rm log(\:\:}d\:\:{\rm [normalized])}$';
elseif(prop==22)
    log_prop = 1;
    norm_stack = 0;
    yl = -1;
    yu = 3;
    dy = 0.5;
    normyl = -2;
    normyu = 2;
    normdy = 0.5;
    dirname = strcat(dirname,'height/');
    fname = 'height';
    ylab = '${\rm log(}h\:{\rm [kpc])}$';
    ylab2 = '${\rm log(\:\:}h\:\:{\rm [kpc\:])}$';
    norm_ylab = '${\rm log(}h\:{\rm [normalized])}$';
    norm_ylab2 = '${\rm log(\:\:}h\:\:{\rm [normalized])}$';
elseif(prop==44)
    is(:,44) = 2e6.*is(:,44).*is(:,30)./is(:,6);
    log_prop = 0;
    norm_stack = 0;
    yl = 0.02;
    yu = 2;
    dy = 0.1;
    ytick = [0.02 0.1 1 2];
    ystr = {num2str(0.02,'%4.2f'),num2str(0.1,'%3.1f'),num2str(1,'%1.0f'),num2str(2,'%1.0f')};
    dirname = strcat(dirname,'alpha_stars/');
    fname = 'alphas';
    ylab = '$\alpha_{\rm s}$';
    ylab2 = '$\alpha_{\rm s}$';
elseif(prop==49)
    log_prop = 1;
    norm_stack = 0;
    yl = 0;
    yu = 2;
    dy = 0.2;
    normyl = -1;
    normyu = 1;
    normdy = 0.2;
    dirname = strcat(dirname,'td_local_over_tff/');
    fname = 'td_over_local_tff';
    ylab = '${\rm log(}t_{\rm d,\,l}/t_{\rm ff})$';
    ylab2 = '${\rm log(\:\:}t_{\rm d,\,l}\:/t_{\rm ff}\:)$';
elseif(prop==50)
    log_prop = 1;
    norm_stack = 0;
    yl = 0;
    yu = 2;
    dy = 0.2;
    normyl = -1;
    normyu = 1;
    normdy = 0.2;
    dirname = strcat(dirname,'td_global_over_tff/');
    fname = 'td_global_over_tff';
    ylab = '${\rm log(}t_{\rm d,\,g}/t_{\rm ff})$';
    ylab2 = '${\rm log(\:\:}t_{\rm d,\,g}\:/t_{\rm ff}\:)$';
elseif(prop==65)
    is(:,65) = is(:,65).*is(:,50).*(is(:,4)+is(:,15).*30.*1e6)./is(:,4);
    log_prop = 0;
    norm_stack = 0;
    yl = 0.03;
    yu = 3;
    dy = 0.1;
    ytick = [0.03 0.1 1 3];
    ystr = {num2str(0.03,'%4.2f'),num2str(0.1,'%3.1f'),num2str(1,'%1.0f'),num2str(3,'%1.0f')};
    dirname = strcat(dirname,'epsd/');
    fname = 'epsd';
    ylab = '$\epsilon_{\rm d}$';
    ylab2 = '$\epsilon_{\rm d}$';
elseif(prop==66)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 0.6;
    dy = 0.2;
    normyl = -2;
    normyu = 0.6;
    normdy = 0.2;
    dirname = strcat(dirname,'mass_loading_Rc_Vr_0/');
    fname = 'etag';
    ylab = '${\rm log(}\eta)$';
    ylab2 = '${\rm log(\:\:}\eta\:)$';    
elseif(prop==67)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 0.6;
    dy = 0.2;
    normyl = -2;
    normyu = 0.6;
    normdy = 0.2;
    dirname = strcat(dirname,'mass_loading_Rc_Vr_Vesc/');
    fname = 'etag';
    ylab = '${\rm log(}\eta)$';
    ylab2 = '${\rm log(\:\:}\eta\:)$';
elseif(prop==68)
    log_prop = 0;
    norm_stack = 0;
    yl = 0.1;
    yu = 10;
    dy = 0.1;
    ytick = [0.1 1 10];
    ystr = {num2str(0.1,'%3.1f'),num2str(1,'%1.0f'),num2str(10,'%2.0f')};
    dirname = strcat(dirname,'mass_loading_Rc_Vr_0_V_Vesc/');
    fname = 'etag';
    ylab  = '$\eta$';
    ylab2 = '$\eta$';
elseif(prop==69)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 0.6;
    dy = 0.2;
    normyl = -2;
    normyu = 0.6;
    normdy = 0.2;
    dirname = strcat(dirname,'mass_loading_avg_Vr_0/');
    fname = 'etag';
    ylab = '${\rm log(}\eta)$';
    ylab2 = '${\rm log(\:\:}\eta\:)$';
elseif(prop==70)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 0.6;
    dy = 0.2;
    normyl = -2;
    normyu = 0.6;
    normdy = 0.2;
    dirname = strcat(dirname,'mass_loading_avg_Vr_Vesc/');
    fname = 'etag';
    ylab = '${\rm log(}\eta)$';
    ylab2 = '${\rm log(\:\:}\eta\:)$';
elseif(prop==71)
    is(:,71) = min(is(:,71),100);
    log_prop = 0;
    norm_stack = 0;
    yl = 0.1;
    yu = 10;
    dy = 0.1;
    ytick = [0.1 1 10];
    ystr = {num2str(0.1,'%3.1f'),num2str(1,'%1.0f'),num2str(10,'%2.0f')};
    dirname = strcat(dirname,'mass_loading_avg_Vr_0_V_Vesc/');
    fname = 'etag';
    ylab = '$\eta$';
    ylab2 = '$\eta$';
elseif(prop==72)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 0.6;
    dy = 0.2;
    normyl = -2;
    normyu = 0.6;
    normdy = 0.2;
    dirname = strcat(dirname,'mass_loading_Frederic/');
    fname = 'etag';
    ylab = '${\rm log(}\eta)$';
    ylab2 = '${\rm log(\:\:}\eta\:)$';
elseif(prop==73)
    log_prop = 0;
    norm_stack = 0;
    yl = 0.01;
    yu = 100;
    dy = 0.1;
    ytick = [0.01 0.1 1 10 100];
    ystr = {num2str(0.01,'%4.2f'),num2str(0.1,'%3.1f'),num2str(1,'%1.0f'),num2str(10,'%2.0f'),num2str(100,'%3.0f')};
    dirname = strcat(dirname,'eta_stars/');
    fname = 'etas';
    ylab = '${\rm log(}\eta_{\rm s})$';
    ylab2 = '${\rm log(\:\:}\eta_{\rm s}\:\:)$';
elseif(prop==74)
    is2 = -is(:,74);
    is(:,74) = max(is(:,74),-100);
    is(:,74) = min(is(:,74),100);
    is2 = max(is2,-100);
    is2 = min(is2,100);
    log_prop = 0;
    norm_stack = 0;
    yl = 0.01;
    yu = 10;
    dy = 0.1;
    ytick = [0.01 0.1 1 10];
    ystr = {num2str(0.01,'%4.2f'),num2str(0.1,'%3.1f'),num2str(1,'%1.0f'),num2str(10,'%2.0f')};
    dirname = strcat(dirname,'eta_stars_net/');
    fname = 'etas';
    ylab = '$|\eta_{\rm s}|$';
    ylab2 = '$|\eta_{\rm s}|$';
elseif(prop==80)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_avg_Vr_0/');
    fname = 'alpha';
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==81)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_Rc_Vr_Vesc_gas_cons/');
    fname = 'alpha';
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==82)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_Rc_Vr_0_V_Vesc_gas_cons/');
    fname = 'alpha';
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==83)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_avg_Vr_Vesc_gas_cons/');
    fname = 'alpha';
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==84)
    log_prop = 0;
    norm_stack = 0;
    yl = 0.02;
    yu = 2;
    dy = 0.1;
    ytick = [0.02 0.1 1 2];
    ystr = {num2str(0.02,'%4.2f'),num2str(0.1,'%3.1f'),num2str(1,'%1.0f'),num2str(2,'%1.0f')};
    dirname = strcat(dirname,'alpha_acc_avg_Vr_0_V_Vesc_gas_cons/');
    fname = 'alpha';
    ylab = '$\alpha$';
    ylab2 = '$\alpha$';
elseif(prop==85)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_Frederic_gas_cons/');
    fname = 'alpha';
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==92)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_Rc_Vr_Vesc_bar_cons/');
    fname = 'alpha';
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==93)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_Rc_Vr_0_V_Vesc_bar_cons/');
    fname = 'alpha';
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==94)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_avg_Vr_Vesc_bar_cons/');
    fname = 'alpha';
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==95)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_avg_Vr_0_V_Vesc_bar_cons/');
    fname = 'alpha';
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
elseif(prop==96)
    log_prop = 1;
    norm_stack = 0;
    yl = -2;
    yu = 2;
    dy = 0.4;
    normyl = -2;
    normyu = 2;
    normdy = 0.4;
    dirname = strcat(dirname,'alpha_acc_Frederic_bar_cons/');
    fname = 'alpha';
    ylab = '${\rm log(}\alpha)$';
    ylab2 = '${\rm log(\:\:}\alpha\:)$';
end
if(versus_dist==1)
    fname = strcat(fname,'_vs_dist');
end
mkdir(dirname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

col(:,1)  = [0 0 0];
col(:,2)  = [1 0 0];
col(:,3)  = [0 1 0];
col(:,4)  = [0 0 1];
col(:,5)  = [0 1 1];
col(:,6)  = [1 0 1];
col(:,7)  = [1 1 0];
col(:,8)  = [0.7 0.7 0.7];
col(:,9)  = [0.2 0.7 0.4];
col(:,10) = [0.7 0.4 0.2];
col(:,11)  = [0.1 0.2 0.8];
col(:,12)  = [0.8 0.1 0.5];
col(:,13) = [0.3 0.8 0.3];
Ncol = 13;

if(versus_dist==0)
    xlab = '$t\:/\:t_{\rm d}$';
    xlab2 = '$t\:\:/\:\:t_{\rm d}$';
    xl = 0.3;
    xu = 30;
    xtick = [0.03 0.3 1 10 30];
    xstr = {num2str(0.03,'%4.2f'),num2str(0.3,'%3.1f'),num2str(1,'%1.0f'),num2str(10,'%2.0f'),num2str(30,'%2.0f')};
    tmin = 0;
    tmax = 40;
    dt = 0.25;
    tbin = (tmin+dt/2):dt:(tmax-dt/2);
%    dt = 0.1;
%    tbin = 10.^[-1:0.1:1.8];
    Nt = length(tbin);
elseif(versus_dist==1)
    xlab = '$d\:/\:R_{\rm d}$';
    xlab2 = '$d\:\:/\:\:R_{\rm d}$';
    xl = 0.1;
    xu = 2;
    xtick = [0.1 1 2];
    xstr = {num2str(0.1,'%3.1f'),num2str(1,'%1.0f'),num2str(2,'%1.0f')};
    dt = 0.1;
    tbin = 10.^[-1:0.1:0.6];
    Nt = length(tbin);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(gal~=0)
    m = find(ind==gal);
    vec = (nis(m)+1):(nis(m+1)-1);
    vec2 = vec+1;
    
    b = find(is(vec,62)>=Tdown & is(vec,62)<Tup & is(vec,1)~=is(vec2,1));
    end_index = vec(b);
    if( is(nis(m+1),62)>=Tdown & is(nis(m+1),62)<Tup )
        end_index = [end_index, nis(m+1)];
    end
    
    b = find(is(vec2,62)>=Tdown & is(vec2,62)<Tup & is(vec2,1)~=is(vec,1));
    start_index = vec2(b);
    if( is(nis(m)+1,62)>=Tdown & is(nis(m)+1,62)<Tup )
        start_index = [nis(m)+1, start_index];
    end
else    
    end_index = [];
    start_index = [];
    for m=1:length(ind)
        vec = (nis(m)+1):(nis(m+1)-1);
        vec2 = vec+1;
        
        b = find(is(vec,62)>=Tdown & is(vec,62)<Tup & is(vec,1)~=is(vec2,1));
        end_index = [end_index,vec(b)];
        if( is(nis(m+1),62)>=Tdown & is(nis(m+1),62)<Tup )
            end_index = [end_index, nis(m+1)];
        end        
        
        if( is(nis(m)+1,62)>=Tdown & is(nis(m)+1,62)<Tup )
            start_index = [nis(m)+1, start_index];
        end
        b = find(is(vec2,62)>=Tdown & is(vec2,62)<Tup & is(vec2,1)~=is(vec,1));
        start_index = [start_index,vec2(b)];
    end
end
b = find(is(start_index,2)>=zdown & is(start_index,2)<zup);
start_index = start_index(b);
end_index = end_index(b);
Mavg = 0.*end_index;
SFRavg = 0.*end_index;
for i=1:length(start_index)
    tvec = start_index(i):end_index(i);
    b = find(abs(is(tvec,48)-50)<=10);
    if(isempty(b))
        b = find(abs(is(tvec,48)-50)==min(abs(is(tvec,48)-50)));
    end
    Mavg(i) = mean(log10(is(tvec(b),6)));
    b = find(abs(is(tvec,48)-200)<=20);
    if(isempty(b))
        b = find(abs(is(tvec,48)-200)==min(abs(is(tvec,48)-200)));
    end
    SFRavg(i) = mean(log10(is(tvec(b),15)));
end
b = find(Mavg>=Mdown & Mavg<Mup & SFRavg>=SFRdown & SFRavg<SFRup);
[gal, length(find(Mavg>=Mdown & Mavg<Mup)), length(b), median(Mavg), median(SFRavg)]
start_index = start_index(b);
end_index = end_index(b);
Mavg   = Mavg(b);
SFRavg = SFRavg(b);
[gal, median(Mavg), median(SFRavg)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( length(start_index) ~= length(end_index) )
    [0,0,length(start_index), length(end_index)]
    return
elseif( min(end_index-start_index)<=0 )
    [0,0,min(end_index-start_index)]
    return
else
    
    Nclump = length(start_index);
    Nindex = max(end_index-start_index) + 1;
    ydat  = 1e50.*ones(Nclump,Nindex);
    ydat2 = ydat;
    ydat3 = ydat;
    zdat  = ydat;
    norm_ydat = ydat;
    norm_zdat = ydat;
    clump_time = ydat;
    stack  = 1e50.*ones(3,Nt);
    stack2 = 1e50.*ones(3,Nt);
    min(is(end_index(1:Nclump),2))
    [Nclump, Nindex]
    if(versus_dist==0)
        if(log_smoothing==0)
            tH = 1/2;
        else
            tH = 0.1;
        end
    elseif(versus_dist==1)
        if(log_smoothing==0)
            tH = 0.3;
        else
            tH = 0.10;
        end
    end
    
    med_mass = zeros(1,Nclump);
    for i=1:Nclump
        % Define vectors
        j = start_index(i);
        k = end_index(i);
        tsfc = is(j,prop);
        Msi  = is(j,5);
        for n=j:k
            ydat(i,n-j+1) = is(n,prop);
            if(prop==74)
                zdat(i,n-j+1) = is2(n);
            elseif(prop==177)
                %zdat(i,n-j+1) = max(0.1*is(j,30), is(n,prop)-tsfc);
                zdat(i,n-j+1) = max(0.1*is(j,30), 1e-6.*(is(n,5)-Msi)/is(n,15));
            end
            ydat2(i,n-j+1) = is(n,6);   % Mc
            ydat3(i,n-j+1) = is(n,30);  % td_global
            if(versus_dist==0)
                clump_time(i,n-j+1) = max(0.1,is(n,48)./is(n,30));
            elseif(versus_dist==1)
                clump_time(i,n-j+1) = max(0.05,is(n,19));
            end
        end
        xdat = clump_time(i,1:(k-j+1));
        [min(ydat(i,1:(k-j+1))), max(ydat(i,1:(k-j+1)))]
        % Normalize
        if(norm_stack == 1 & versus_dist==0)
            if(prop==177 || prop==18 || prop==12)
                norm_ydat(i,1:(k-j+1)) = ydat(i,1:(k-j+1))./ydat3(i,1:(k-j+1));
                if(prop==177)
                    norm_zdat(i,1:(k-j+1)) = zdat(i,1:(k-j+1))./ydat3(i,1:(k-j+1));
                end
            else
                b = find(abs(xdat-tform)<=0.25); % Window of 0.25*tdyn
                if(isempty(b))
                    b = find(abs(xdat-tform)==min(abs(xdat-tform)));
                end
                if(prop==4 || prop==5)
                    avg_prop = mean(ydat2(i,b));
                    med_mass(i) = avg_prop;
                elseif(prop==15)
                    avg_prop = mean(ydat2(i,b)./(1e6.*ydat3(i,b)));
                else
                    avg_prop = mean(ydat(i,b));
                    if(prop==74)
                        avg_prop2 = mean(zdat(i,b));
                    end
                end
                norm_ydat(i,1:(k-j+1)) = ydat(i,1:(k-j+1))./avg_prop;
                if(prop==74)
                    norm_zdat(i,1:(k-j+1)) = zdat(i,1:(k-j+1))./avg_prop2;
                end
            end
        else
            norm_ydat(i,1:(k-j+1)) = ydat(i,1:(k-j+1));
            if(prop==74 || prop==177)
                norm_zdat(i,1:(k-j+1)) = zdat(i,1:(k-j+1));
            end
        end
        
        % Smooth
        if(prop<=20 || prop==177)
            norm_ydat(i,1:(k-j+1)) = log10(norm_ydat(i,1:(k-j+1)));
            if(prop==177)
                norm_zdat(i,1:(k-j+1)) = log10(norm_zdat(i,1:(k-j+1)));
            end
        end
        if(log_smoothing==1)
            xdat = log10(xdat);
        end
        temp = norm_ydat(i,1:(k-j+1));
        for n=1:(k-j+1)
            % Gaussian with HWHM=tH
            temp(n) = sum( norm_ydat(i,1:(k-j+1)).*0.5.^(((xdat-xdat(n))./tH).^2) ) ./ sum( 0.5.^(((xdat-xdat(n))./tH).^2) );
        end
        norm_ydat(i,1:(k-j+1)) = temp;
        if(prop==74 || prop==177)
            temp2 = norm_zdat(i,1:(k-j+1));
            for n=1:(k-j+1)
                % Gaussian with HWHM=tH
                temp2(n) = sum( norm_zdat(i,1:(k-j+1)).*0.5.^(((xdat-xdat(n))./tH).^2) ) ./ sum( 0.5.^(((xdat-xdat(n))./tH).^2) );
            end
            norm_zdat(i,1:(k-j+1)) = temp2;
        end
        if(prop<=20 || prop==177)
            norm_ydat(i,1:(k-j+1)) = 10.^norm_ydat(i,1:(k-j+1));
            if(prop==177)
                norm_zdat(i,1:(k-j+1)) = 10.^(norm_zdat(i,1:(k-j+1)));
            end
        end
        if(versus_dist==1)
            clump_time(i,1:(k-j+1)) = log10(clump_time(i,1:(k-j+1)));
        end
    end
    if(prop==4)
        [gal, log10(median(med_mass))]
    end
    clear temp xdat ydat2 ydat3
    
    if(versus_dist==1)
        tbin = log10(tbin);
    end
    for i=1:Nt
        b = find(clump_time>(tbin(i)-dt/2) & clump_time<=(tbin(i)+dt/2) & norm_ydat~=1e50);
        temp_vec = sort(norm_ydat(b));
        Ntemp_vec = length(b);
        if(Ntemp_vec>=3)
            Nlow  = floor(Ntemp_vec/6) + 1;
            Nhigh = floor(5*Ntemp_vec/6);
            
            stack(1,i) = median(temp_vec);
            stack(2,i) = temp_vec(Nlow);
            stack(3,i) = temp_vec(Nhigh);
        end
        if(prop==74 || prop==177)
            b = find(clump_time>(tbin(i)-dt/2) & clump_time<=(tbin(i)+dt/2) & norm_zdat~=1e50);
            temp_vec = sort(norm_zdat(b));
            Ntemp_vec = length(b);
            if(Ntemp_vec>=3)
                Nlow  = floor(Ntemp_vec/6) + 1;
                Nhigh = floor(5*Ntemp_vec/6);
                
                stack2(1,i) = median(temp_vec);
                stack2(2,i) = temp_vec(Nlow);
                stack2(3,i) = temp_vec(Nhigh);
            end
        end
    end
    if(versus_dist==1)
        tbin = 10.^tbin;
        for i=1:Nclump
            j = start_index(i);
            k = end_index(i);
            clump_time(i,1:(k-j+1)) = 10.^clump_time(i,1:(k-j+1));
        end
    end
    if(prop<=20 || prop==177)
        stack(:,1:Nt) = log10(stack(:,1:Nt));
        b = find(stack==50);
        stack(b) = 1e50;
        if(prop==177)
            stack2(:,1:Nt) = log10(stack2(:,1:Nt));
            b = find(stack2==50);
            stack2(b) = 1e50;
        end
    end
    clear temp_vec Ntemp_vec Nlow Nhigh
    b = find( stack(1,1:Nt)~=1e50 );
    temp = stack(1,b);
    if(log_smoothing==0)
        xdat = tbin(b);
    else
        xdat = log10(tbin(b));
    end
    for n=1:length(b)
        % Gaussian with HWHM=tH
        temp(n) = sum( stack(1,b).*0.5.^(((xdat-xdat(n))./tH).^2) ) ./ sum( 0.5.^(((xdat-xdat(n))./tH).^2) );
    end
    stack(1,b) = temp;
    if(prop<=20 || prop==177)
        stack(1,b) = 10.^stack(1,b);
    end
    if(prop==74 || prop==177)
        b = find( stack2(1,1:Nt)~=1e50 );
        temp = stack2(1,b);
        if(log_smoothing==0)
            xdat = tbin(b);
        else
            xdat = log10(tbin(b));
        end
        for n=1:length(b)
            % Gaussian with HWHM=tH
            temp(n) = sum( stack2(1,b).*0.5.^(((xdat-xdat(n))./tH).^2) ) ./ sum( 0.5.^(((xdat-xdat(n))./tH).^2) );
        end
        stack2(1,b) = temp;
        if(prop==177)
            stack2(1,b) = 10.^stack2(1,b);
        end
    end
    
    b = find( stack(2,1:Nt)~=1e50 );
    temp = stack(2,b);
    if(log_smoothing==0)
        xdat = tbin(b);
    else
        xdat = log10(tbin(b));
    end
    for n=1:length(b)
        % Gaussian with HWHM=tH
        temp(n) = sum( stack(2,b).*0.5.^(((xdat-xdat(n))./tH).^2) ) ./ sum( 0.5.^(((xdat-xdat(n))./tH).^2) );
    end
    stack(2,b) = temp;
    if(prop<=20 || prop==177)
        stack(2,b) = 10.^stack(2,b);
    end
    if(prop==74 || prop==177)
        b = find( stack2(2,1:Nt)~=1e50 );
        temp = stack2(2,b);
        if(log_smoothing==0)
            xdat = tbin(b);
        else
            xdat = log10(tbin(b));
        end
        for n=1:length(b)
            % Gaussian with HWHM=tH
            temp(n) = sum( stack2(2,b).*0.5.^(((xdat-xdat(n))./tH).^2) ) ./ sum( 0.5.^(((xdat-xdat(n))./tH).^2) );
        end
        stack2(2,b) = temp;
        if(prop==177)
            stack2(2,b) = 10.^stack2(2,b);
        end
    end
    
    b = find( stack(3,1:Nt)~=1e50 );
    temp = stack(3,b);
    if(log_smoothing==0)
        xdat = tbin(b);
    else
        xdat = log10(tbin(b));
    end
    for n=1:length(b)
        % Gaussian with HWHM=tH
        temp(n) = sum( stack(3,b).*0.5.^(((xdat-xdat(n))./tH).^2) ) ./ sum( 0.5.^(((xdat-xdat(n))./tH).^2) );
    end
    stack(3,b) = temp;
    if(prop<=20 || prop==177)
        stack(3,b) = 10.^stack(3,b);
    end
    if(prop==74 || prop==177)
        b = find( stack2(3,1:Nt)~=1e50 );
        temp = stack2(3,b);
        if(log_smoothing==0)
            xdat = tbin(b);
        else
            xdat = log10(tbin(b));
        end
        for n=1:length(b)
            % Gaussian with HWHM=tH
            temp(n) = sum( stack2(3,b).*0.5.^(((xdat-xdat(n))./tH).^2) ) ./ sum( 0.5.^(((xdat-xdat(n))./tH).^2) );
        end
        stack2(3,b) = temp;
        if(prop==177)
            stack2(3,b) = 10.^stack2(3,b);
        end
    end
    clear temp xdat
    if(prop==74)
        stack = max(1e-6, stack);
        stack2 = max(1e-6, stack2);
    end
    
    if(log_prop==1)
        b = find(stack ~= 1e50);
        b1 = find(stack(b)<=0);
        if(~isempty(b1))
            stack(b(b1)) = -10;
        end
        b1 = find(stack(b)>0);
        if(~isempty(b1))
            stack(b(b1)) = log10(stack(b(b1)));
        end
        
        b = find(ydat ~= 1e50);
        b1 = find(ydat(b)<=0);
        if(~isempty(b1))
            ydat(b(b1)) = -10;
        end
        b1 = find(ydat(b)>0);
        if(~isempty(b1))
            ydat(b(b1)) = log10(ydat(b(b1)));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filename1 = strcat(dirname,'INDIVIDUAL__',num2str(Tdown,'%04i'),'_Tmax_',num2str(Tup,'%04i'));
%     filename2 = strcat(dirname,'STACKED__',num2str(Tdown,'%04i'),'_Tmax_',num2str(Tup,'%04i'));
    filename2 = strcat(dirname,'V',num2str(gal,'%02i'),'_',fname);
    filename1e = strcat(filename1,'.eps');
    filename1j = strcat(filename1,'.jpg');
    filename2e = strcat(filename2,'.eps');
    filename2j = strcat(filename2,'.jpg');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
    set(figure1,'WindowStyle','docked');
    if(norm_stack == 0)
        normyl=yl;
        normdy=dy;
        normyu=yu;
        norm_ylab = ylab;
        norm_ylab2 = ylab2;
        if(prop==177 || prop==84 || prop==74 || prop==73 || prop==71 || prop==65 || prop==44 || prop==18 || prop==15 || prop==12 || prop==7 || prop==6 || prop==5 || prop==4 || prop==3)
            normytick = ytick;
            normystr = ystr;
        end
    end
    if(prop==177 || prop==84 || prop==74 || prop==73 || prop==71 || prop==65 || prop==44 || prop==18 || prop==15 || prop==12 ||prop==7 || prop==6 || prop==5 || prop==4 || prop==3)
        axes1 = axes('Parent',figure1,'YTick',normytick,...
            'YScale','log',...
            'YtickLabel',normystr,...
            'XTick',xtick,...
            'XScale','log',...
            'XtickLabel',xstr,...%,xl:dx:xu,...'XMinorTick','on',...
            'TickLength',[0.02 0.04],...
            'PlotBoxAspectRatio',[1.618 1 1],...
            'LineWidth',1.5,...
            'FontSize',16,...
            'FontName','Arial',...
            'Position',[0.13 0.14 0.775 0.815]);
    else
        axes1 = axes('Parent',figure1,'YTick',normyl:normdy:normyu,...
            'YMinorTick','on',...
            'XTick',xtick,...
            'XScale','log',...
            'XtickLabel',xstr,...%,xl:dx:xu,...'XMinorTick','on',...
            'TickLength',[0.02 0.04],...
            'PlotBoxAspectRatio',[1.618 1 1],...
            'LineWidth',1.5,...
            'FontSize',16,...
            'FontName','Arial',...
            'Position',[0.13 0.14 0.775 0.815]);
    end
    xlim(axes1,[xl xu]);
    if(prop==177 || prop==84 || prop==74 || prop==73 || prop==71 || prop==65 || prop==44 || prop==18 || prop==15 || prop==7 || prop==6 || prop==5 || prop==4 || prop==3)
        ylim(axes1,[normyl normyu]);
    else
        ylim(axes1,[normyl normyu]);
    end
    box(axes1,'on');
    % grid on
    hold(axes1,'all');    
    xlabel(xlab,'Interpreter','latex','FontSize',22,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.08 0]);    
    ylabel(norm_ylab,'Interpreter','latex','FontSize',22,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.1 0.5 0]);    
    title(tit,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 0.90 0]);    
    set(gcf,'renderer','painters')
    
    for i=1:Nclump
        j = start_index(i);
        k = end_index(i);
        
        m = mod(i,10);
        m = m+1;
        if(mod(floor((i-1)/10),2)==0)
            lin = '-';
        else
            lin = '--';
        end
        if(prop==74)
            plot(clump_time(i,1:(k-j+1)),norm_ydat(i,1:(k-j+1)),...
                'marker','none',...
                'linestyle',lin,'linewidth',3,'color',col(:,m));
            plot(clump_time(i,1:(k-j+1)),norm_zdat(i,1:(k-j+1)),...
                'marker','none',...
                'linestyle',lin,'linewidth',1,'color',col(:,m));
        else
            plot(clump_time(i,1:(k-j+1)),norm_ydat(i,1:(k-j+1)),...
                'marker','none',...
                'linestyle',lin,'linewidth',2,'color',col(:,m));
        end
%         plot(log10(is(end_index(i),45)),last_snap(i),...
%             'linestyle','none',...
%             'marker','square','MarkerSize',10,'MarkerFaceColor',col(:,m),'MarkerEdgeColor',col(:,m));
    end
    print(gcf,'-depsc',filename1e);
    
    set(gcf,'visible','off');
    xlabel(xlab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.08 0]);    
    ylabel(ylab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.1 0.5 0]);    
    title(tit,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 0.90 0]);
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1)+0.1 ti(2)+0.04 0.84-ti(3)-ti(1) 1-ti(4)-ti(2)]);
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3)+0.0 pos(4)+ti(2)+ti(4)+5]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3)+2.6 pos(4)+ti(2)+ti(4)+1]);
    print(gcf,'-djpeg',filename1j);
    close all
    fclose all;

    figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
    set(figure1,'WindowStyle','docked');
    if(norm_stack == 0)
        normyl=yl;
        normdy=dy;
        normyu=yu;
        norm_ylab = ylab;
        norm_ylab2 = ylab2;
        if(prop==177 || prop==84 || prop==74 || prop==73 || prop==71 || prop==65 || prop==44 || prop==19 || prop==18 || prop==15 || prop==12 || prop==7 || prop==6 || prop==5 || prop==4 || prop==3)
            normytick = ytick;
            normystr = ystr;
        end
    end
    if(prop==177 || prop==84 || prop==74 || prop==73 || prop==71 || prop==65 || prop==44 || prop==19 || prop==18 || prop==15 || prop==12 || prop==7 || prop==6 || prop==5 || prop==4 || prop==3)
        axes1 = axes('Parent',figure1,'YTick',normytick,...
            'YScale','log',...
            'YtickLabel',normystr,...
            'XTick',xtick,...
            'XScale','log',...
            'XtickLabel',xstr,...%,xl:dx:xu,...'XMinorTick','on',...
            'TickLength',[0.02 0.04],...
            'PlotBoxAspectRatio',[1.618 1 1],...
            'LineWidth',1.5,...
            'FontSize',16,...
            'FontName','Arial',...
            'Position',[0.1500 0.0500 0.8100 0.6500]);%[0.13 0.14 0.775 0.815]);
    else
        axes1 = axes('Parent',figure1,'YTick',normyl:normdy:normyu,...
            'YMinorTick','on',...
            'XTick',xtick,...
            'XScale','log',...
            'XtickLabel',xstr,...%,xl:dx:xu,...'XMinorTick','on',...
            'TickLength',[0.02 0.04],...
            'PlotBoxAspectRatio',[1.618 1 1],...
            'LineWidth',1.5,...
            'FontSize',16,...
            'FontName','Arial',...
            'Position',[0.13 0.14 0.775 0.815]);
    end
    xlim(axes1,[xl xu]);
    if(prop==177 || prop==84 || prop==74 || prop==73 || prop==71 || prop==65 || prop==44 || prop==19 || prop==18 || prop==15 || prop==12 || prop==7 || prop==6 || prop==5 || prop==4 || prop==3)
        ylim(axes1,[normyl normyu]);
    else
        ylim(axes1,[normyl normyu]);
    end
    box(axes1,'on');
    % grid on
    hold(axes1,'all');    
    xlabel(xlab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.05 0]);    
    ylabel(norm_ylab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.11 0.5 0]);
    title(tit,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.90 0.85 0]);
    set(gcf,'renderer','painters')
    
    if(prop==74)
        b = find(stack(1,:) ~= 1e50 & tbin<tend);
        l(1) = plot(tbin(b),stack(1,b),'marker','none',...
            'linestyle','-','linewidth',3,'color','k','DisplayName','${\rm Sims,\:\eta_{\rm s}>0}$');
        b = find(stack(2,:) ~= 1e50 & tbin<=30);
        plot(tbin(b),stack(2,b),'marker','none',...
            'linestyle','--','linewidth',2,'color','k');
        b = find(stack(3,:) ~= 1e50 & tbin<=30);
        plot(tbin(b),stack(3,b),'marker','none',...
            'linestyle','--','linewidth',2,'color','k');    
        b = find(stack2(1,:) ~= 1e50 & tbin<tend);
        
        l(2) = plot(tbin(b),stack2(1,b),'marker','none',...
            'linestyle','-','linewidth',3,'color','g','DisplayName','${\rm Sims,\:\eta_{\rm s}<0}$');
        b = find(stack2(2,:) ~= 1e50 & tbin<=30);
        plot(tbin(b),stack2(2,b),'marker','none',...
            'linestyle','--','linewidth',2,'color','g');
        b = find(stack2(3,:) ~= 1e50 & tbin<=30);
        plot(tbin(b),stack2(3,b),'marker','none',...
            'linestyle','--','linewidth',2,'color','g');
    elseif(prop==177)
        b = find(stack(1,:) ~= 1e50 & tbin<tend);
        l(1) = plot(tbin(b),stack(1,b),'marker','none',...
            'linestyle','-','linewidth',3,'color','k','DisplayName','${\rm Sims,\:t_{\rm sf}}$');
        b = find(stack(2,:) ~= 1e50 & tbin<=30);
        plot(tbin(b),stack(2,b),'marker','none',...
            'linestyle','--','linewidth',2,'color','k');
        b = find(stack(3,:) ~= 1e50 & tbin<=30);
        plot(tbin(b),stack(3,b),'marker','none',...
            'linestyle','--','linewidth',2,'color','k');
        
        b = find(stack2(1,:) ~= 1e50 & tbin<tend);        
        l(2) = plot(tbin(b),stack2(1,b),'marker','none',...
            'linestyle','-.','linewidth',2,'color','m','DisplayName','${\rm Sims,\:t_{\rm sfc}}$');
    else
        b = find(stack(1,:) ~= 1e50 & tbin<tend);
        l(1) = plot(tbin(b),stack(1,b),'marker','none',...
            'linestyle','-','linewidth',3,'color','k','DisplayName','${\rm Sims}$');
        b = find(stack(2,:) ~= 1e50 & tbin<=30);
        plot(tbin(b),stack(2,b),'marker','none',...
            'linestyle','--','linewidth',2,'color','k');
        b = find(stack(3,:) ~= 1e50 & tbin<=30);
        plot(tbin(b),stack(3,b),'marker','none',...
            'linestyle','--','linewidth',2,'color','k');
    end
    
    if( (prop==4 || prop==5 || prop==6 || prop==7 || prop==12 || prop==15 || prop==17 || prop==18 || prop==177 ...
            || prop==65 || prop==71 || prop==74 || prop==84) & versus_dist==0 )
        [T,Y]=ode45(@(t,y)clump_evolution_exact(t, y, alpha, eta1, eps, mu, etas, 1, tmig_td),[0 300],[fg_i 1-fg_i 1]);
        Mg = Y(:,1);
        Ms = Y(:,2);
        R = Y(:,3);
        time = T(:,1);
        Mb = Mg + Ms;
        fg = (Mg./Mb);
        fs = 1-fg;
        SFR = eps.*Mg;  %SFR / (Mci/tdyn) = SFR*tdyn/Mci
        sSFR = eps.*(Mg./Ms);  %SFR / (Mci/tdyn) / (Ms/Mci) = SFR*tdyn/Ms
        tdep = ones(size(Ms))./eps;
        tsf = tdep.*(Ms./Mg);
        b = find(abs(time-tform) == min(abs(time-tform)));
        b1 = find(time<tend);
        if(prop==4)
            l(2) = plot(time(b1),Mg(b1)./Mb(b),'marker','none','linestyle','-','color','r','linewidth',2,'DisplayName','${\rm Model}$');
        elseif(prop==5)
            l(2) = plot(time(b1),Ms(b1)./Mb(b),'marker','none','linestyle','-','color','r','linewidth',2,'DisplayName','${\rm Model}$');
        elseif(prop==6)
            l(2) = plot(time(b1),Mb(b1)./Mb(b),'marker','none','linestyle','-','color','r','linewidth',2,'DisplayName','${\rm Model}$');
        elseif(prop==7)
            l(2) = plot(time(b1),fg(b1),'marker','none','linestyle','-','color','r','linewidth',2,'DisplayName','${\rm Model}$');
        elseif(prop==15)
            l(2) = plot(time(b1),SFR(b1)./Mb(b),'marker','none','linestyle','-','color','r','linewidth',2,'DisplayName','${\rm Model}$');
        elseif(prop==17)
            l(2) = plot(time(b1),log10(sSFR(b1)./sSFR(b)),'marker','none','linestyle','-','color','r','linewidth',2,'DisplayName','${\rm Model}$');
        elseif(prop==18)
            l(2) = plot(time(b1),tdep(b1),'marker','none','linestyle','-','color','r','linewidth',2,'DisplayName','${\rm Model}$');
        elseif(prop==177)
            l(3) = plot(time(b1),tsf(b1),'marker','none','linestyle','-','color','r','linewidth',2,'DisplayName','${\rm Model}$');
        elseif(prop==65)
            l(2) = plot(time(b1),eps.*ones(size(b1)),'marker','none','linestyle','-','color','r','linewidth',2,'DisplayName','${\rm Model}$');
        elseif(prop==71)
            l(2) = plot(time(b1),eta1.*ones(size(b1)),'marker','none','linestyle','-','color','r','linewidth',2,'DisplayName','${\rm Model}$');
        elseif(prop==74)
            if(etas>=0)
                dname = '${\rm Model,\:\eta_{\rm s}>0}$';
                l(3) = plot(time(b1),abs(etas).*ones(size(b1)),'marker','none',...
                    'linestyle','-','color','r','linewidth',2,'DisplayName',dname);
            else
                dname = '${\rm Model,\:\eta_{\rm s}<0}$';
                l(3) = plot(time(b1),abs(etas).*ones(size(b1)),'marker','none',...
                    'linestyle',':','color','r','linewidth',2,'DisplayName',dname);
            end
        elseif(prop==84)
            l(2) = plot(time(b1),alpha.*ones(size(b1)),'marker','none','linestyle','-','color','r','linewidth',2,'DisplayName','${\rm Model}$');
        end
        if(prop~=177 & prop~=3 & prop~=5 & prop~=12 & prop~=15 & prop~=19 & prop~=65 & prop~=44 & prop~=71 & prop~=74 & prop~=84)
            legend1 = legend(gca,l(1:2));
            set(legend1,...
                'Position',[0.225, 0.2, 0.16, 0.10],'FontSize',16,...
                'Interpreter','latex');
        elseif(prop==5||prop==15||prop==65||prop==71||prop==84)
            legend1 = legend(gca,l(1:2));
            set(legend1,...
                'Position',[0.78, 0.2, 0.16, 0.10],'FontSize',16,...
                'Interpreter','latex');
        elseif(prop==74 & gal==19)
            legend1 = legend(gca,l(1:3));
            set(legend1,...
                'Position',[0.74, 0.754, 0.16, 0.15],'FontSize',16,...
                'Interpreter','latex');
        elseif(prop==74 & gal==7)
            legend1 = legend(gca,l(1:3));
            set(legend1,...
                'Position',[0.25, 0.74, 0.16, 0.15],'FontSize',16,...
                'Interpreter','latex');
        elseif(prop==177)
            x = 0.01:1:1000;
            y = 0.01:1:1000;
            l(4) = plot(x,y,'marker','none','linestyle',':','linewidth',2,'color','k','DisplayName','$t_{\rm sf}=t$');
            legend1 = legend(gca,l(1:4));
            set(legend1,...
                'Position',[0.7, 0.22, 0.16, 0.2],'FontSize',16,...
                'Interpreter','latex');
        elseif(prop==12 & norm_stack==1)
            x = 0.01:1:1000;
            y = 0.01:1:1000;
            l(2) = plot(x,y,'marker','none','linestyle',':','linewidth',2,'color','b','DisplayName','${\rm age_*}=t$');
            legend1 = legend(gca,l(1:2));
            set(legend1,...
                'Position',[0.72, 0.2, 0.16, 0.10],'FontSize',16,...
                'Interpreter','latex');
        end
        legend boxoff
        
        if(prop==5)
            annotation(figure1,'textbox',[0.23 0.53 0.10 0.05],...
                'String',strcat('${\rm \mu}\:\:=',num2str(mu,'%4.2f'),'$'),...
                'FontSize',16,...
                'FontName','Times',...
                'Interpreter','Latex',...
                'FitBoxToText','off',...
                'LineStyle','none',...
                'BackgroundColor',[1 1 1]);
            annotation(figure1,'textbox',[0.23 0.59 0.10 0.05],...
                'String',strcat('${\rm f_{g0}}=',num2str(fg_i,'%4.2f'),'$'),...
                'FontSize',16,...
                'FontName','Times',...
                'Interpreter','Latex',...
                'FitBoxToText','off',...
                'LineStyle','none',...
                'BackgroundColor',[1 1 1]);
            annotation(figure1,'textbox',[0.23 0.70 0.10 0.05],...
                'String',strcat('${\rm \epsilon_{d}}=',num2str(eps,'%4.2f'),'$'),...
                'FontSize',18,...
                'FontName','Times',...
                'Interpreter','Latex',...
                'FitBoxToText','off',...
                'LineStyle','none',...
                'BackgroundColor',[1 1 1]);
            annotation(figure1,'textbox',[0.23 0.76 0.10 0.05],...
                'String',strcat('${\rm \alpha}\:=',num2str(alpha,'%4.2f'),'$'),...
                'FontSize',18,...
                'FontName','Times',...
                'Interpreter','Latex',...
                'FitBoxToText','off',...
                'LineStyle','none',...
                'BackgroundColor',[1 1 1]);
            annotation(figure1,'textbox',[0.23 0.82 0.10 0.05],...
                'String',strcat('${\rm \eta_{\rm s}}=',num2str(etas,'%4.2f'),'$'),...
                'FontSize',18,...
                'FontName','Times',...
                'Interpreter','Latex',...
                'FitBoxToText','off',...
                'LineStyle','none',...
                'BackgroundColor',[1 1 1]);
            annotation(figure1,'textbox',[0.23 0.88 0.10 0.05],...
                'String',strcat('${\rm \eta}\:=',num2str(eta1,'%4.2f'),'$'),...
                'FontSize',18,...
                'FontName','Times',...
                'Interpreter','Latex',...
                'FitBoxToText','off',...
                'LineStyle','none',...
                'BackgroundColor',[1 1 1]);
        end
    end
    print(gcf,'-depsc',filename2e);
    
    set(gcf,'visible','off');
    xlabel(xlab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.05 0]);    
    ylabel(norm_ylab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.11 0.5 0]);
    if(prop==74 & gal==19)
        title(tit,'Interpreter','latex','FontSize',20,...
            'FontName','Times New Roman','units','normalized',...
            'position',[0.4 0.85 0]);
    elseif(prop==74 & gal==7)
        title(tit,'Interpreter','latex','FontSize',20,...
            'FontName','Times New Roman','units','normalized',...
            'position',[0.6 0.85 0]);
    elseif(prop==7)
        title(tit,'Interpreter','latex','FontSize',20,...
            'FontName','Times New Roman','units','normalized',...
            'position',[0.9 0.85 0]);
    else
        title(tit,'Interpreter','latex','FontSize',20,...
            'FontName','Times New Roman','units','normalized',...
            'position',[0.5 0.85 0]);
    end
    
%     ax = gca;

%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset;
%     left = outerpos(1) + ti(1) + 0.02;
%     bottom = outerpos(2) + ti(2) + 0.02;
%     ax_width = outerpos(3) - ti(1) - ti(3) - 0.04;
%     ax_height = outerpos(4) - ti(2) - ti(4) - 0.06;
%     ax.Position = [left bottom ax_width ax_height]

%     ax.Position = [0.1500 0.0500 0.8100 0.9400];
    
%     ti = get(gca,'TightInset');
%     set(gca,'Position',[ti(1)+0.1 ti(2)+0.04 0.84-ti(3)-ti(1) 1-ti(4)-ti(2)]);
%     set(gca,'units','centimeters')
%     pos = get(gca,'Position');
%     set(gcf, 'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3)+0.0 pos(4)+ti(2)+ti(4)+0]);
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3)+2.6 pos(4)+ti(2)+ti(4)+2.6]);
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 16.18 10]);
    set(gca, 'Position',[0.1500 0.0800 0.8100 0.9400]);

    print(gcf,'-djpeg',filename2j);
    print(gcf,'-depsc',filename2e);
    close all
    fclose all;
end