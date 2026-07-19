function [b2,b3,e2,e3,migrate2,migrate3] = clump_evolution(is2, is3, mass_thresh, shape_thresh, ...
    res_thresh, fgas_thresh, dist_thresh, height_thresh, time_thresh, min_snap, ...
    xdat, ydat, stack, mean_type, xnorm, ynorm, legend_location, plot_type, xaxis, yaxis)
% stack == 'no stack' -> don't stack (show individual clumps)
% stack == 'linear stack' -> stack clump data
% stack == 'log stack' -> stack clump data in log space
% *norm == '* absolute' -> absolute values for * (x or y)
% *norm == '* normalized' -> values for * (x or y)
% plot_type == 'no error' -> no error bars on plot
% plot_type == 'error' -> error bars on plot
% mean_type == 'mean' -> use mean and standard deviation when stacking
% mean_type == 'median' -> use median and 67% scatter when stacking
% *axis == '* linear scale' -> linear values with linear scale on * (x or y) axis
% *axis == '* log scale' -> linear values with log scale on * (x or y) axis
% *axis == '* log values' ->  log values on * (x or y) axis

% dist_thresh is the minimum distance a clump is found from the center
% compared to the maximum distance. We may only want migrating clumps.
%prop = 2: redshift / aexp
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
%prop = 18: {\tau}^{-1}
%prop = 19: r/Rd
%prop = 20: z/Hd
%prop = 21: r
%prop = 22: z
%prop = 23: mean_residual
%prop = 24: shape parameter
%prop = 25: dark matter contrast
%prop = 26: classification
%prop = 27: ex situ
%prop = 28: t_ff
%prop = 29: t_d local
%prop = 30: M_dot_gas in from Rc
%prop = 31: M_dot_gas in from 1.5Rc
%prop = 32: M_dot_gas in from 2Rc
%prop = 33: M_dot_gas out from Rc
%prop = 34: M_dot_gas out from 1.5Rc
%prop = 35: M_dot_gas out from 2Rc
%prop = 36: M_dot_stars in
%prop = 37: M_dot_stars out
%prop = 38: M_dot_stars formed
%prop = 39: time since formation
%prop = 40: t_d global
%prop = 41: mean M_dot_gas in
%prop = 42: mean M_dot_gas out
%prop = 43: eta
%prop = 44: eps_ff
%prop = 45: td_local / t_ff
%prop = 46: alpha
%prop = 47: eta_stars
%prop = 48: Mdot_gas out / M_clump
%prop = 49: Mdot_gas in / M_clump
%prop = 50: net Mdot_stars out / M_clump
%prop = 51: Mdot_stars in / M_clump
%prop = 52: Mdot_stars out / M_clump


n2 = length(is2(:,1));
n3 = length(is3(:,1));
begin2 = [];
begin3 = [];
end2 = [];
end3 = [];
m = min_snap - 1;

for i=1:n2-m
    if( is2(i,39)==0 & is2(i+1:i+m,39)>0 & ...
    max(log10(is2(i,6)),log10(is2(i+1,6)))>=mass_thresh(1) & ...
    max(log10(is2(i,6)),log10(is2(i+1,6)))<mass_thresh(2) & ...
    max(log10(is2(i,23)),log10(is2(i+1,23)))>=res_thresh(1) & ...
    max(log10(is2(i,23)),log10(is2(i+1,23)))<res_thresh(2) & ...
    max(is2(i,24),is2(i+1,24))>=shape_thresh(1) & ...
    max(is2(i,24),is2(i+1,24))<shape_thresh(2) & ...
    max(log10(is2(i,7)),log10(is2(i+1,7)))>=fgas_thresh(1) & ... 
    max(log10(is2(i,7)),log10(is2(i+1,7)))<fgas_thresh(2) )
        begin2 = [begin2, i];
        j=i+1;
        while(is2(j,39)>0)
            j = j+1;
        end
        end2 = [end2, j-1];
    end
end

for i=1:n3-m
    if( is3(i,39)==0 & is3(i+1:i+m,39)>0 & ...
    max(log10(is3(i,6)),log10(is3(i+1,6)))>=mass_thresh(1) & ...
    max(log10(is3(i,6)),log10(is3(i+1,6)))<mass_thresh(2) & ...
    max(log10(is3(i,23)),log10(is3(i+1,23)))>=res_thresh(1) & ...
    max(log10(is3(i,23)),log10(is3(i+1,23)))<res_thresh(2) & ...
    max(is3(i,24),is3(i+1,24))>=shape_thresh(1) & ...
    max(is3(i,24),is3(i+1,24))<shape_thresh(2) & ...
    max(log10(is3(i,7)),log10(is3(i+1,7)))>=fgas_thresh(1) & ... 
    max(log10(is3(i,7)),log10(is3(i+1,7)))<fgas_thresh(2) )
        begin3 = [begin3, i];
        j=i+1;
        while(is3(j,39)>0)
            j = j+1;
        end
        end3 = [end3, j-1];
    end
end

if(xdat==2)
    titx = 'a';
    xl = [-1 -1];
    xu = [0 0];
    dx = 0.1;
elseif(xdat==3)
    titx = 'R_{\rm c}\:{\rm [kpc]}';
    xl = [-1 -1];
    xu = [0 0];
    dx = 0.1;
elseif(xdat==4)
    titx = 'M_{\rm c,\:gas}\:{\rm [M_{\odot}]}';
    xl = [2 2];
    xu = [9 9];
    dx = 1;
elseif(xdat==5)
    titx = 'M_{\rm c,\:star}\:{\rm [M_{\odot}]}';
    xl = [2 2];
    xu = [9 9];
    dx = 1;
elseif(xdat==6)
    titx = 'M_{\rm c,\:bar}\:{\rm [M_{\odot}]}';
    xl = [6 6];
    xu = [10 10];
    dx = 1;
elseif(xdat==7)
    titx = 'f_{\rm gas}';
    xl = [-6 -6];
    xu = [0 0];
    dx = 1;
elseif(xdat==8)
    titx = 'f_{\rm dm}';
    xl = [-6 -6];
    xu = [0 0];
    dx = 1;
elseif(xdat==9)
    titx = '\Sigma_{\rm c,\:gas}\:{\rm [M_{\odot}\:pc^{-2}]}';
    xl = [0 0];
    xu = [9 9];
    dx = 1;
elseif(xdat==10)
    titx = '\Sigma_{\rm c,\:star}\:{\rm [M_{\odot}\:pc^{-2}]}';
    xl = [2 2];
    xu = [9 9];
    dx = 1;
elseif(xdat==11)
    titx = '\Sigma_{\rm c,\:bar}\:{\rm [M_{\odot}\:pc^{-2}]}';
    xl = [6 6];
    xu = [10 10];
    dx = 1;
elseif(xdat==12)
    titx = 'age\:{\rm [Myr]}';
    xl = [0.5 0.5];
    xu = [4 4];
    dx = 0.5;
elseif(xdat==13)
    titx = '${\rm log(O\:/\:H)_{gas}}\:+\:12';
    xl = [7.5 7.5];
    xu = [9.5 9.5];
    dx = 0.2;
elseif(xdat==14)
    titx = '${\rm log(O\:/\:H)_{star}}\:+\:12';
    xl = [7.5 7.5];
    xu = [9.5 9.5];
    dx = 0.2;
elseif(xdat==15)
    titx = '$SFR_{\rm c}\:{\rm [M_{\odot}\:yr^{-1}]}';
    xl = [-6 -6];
    xu = [2 2];
    dx = 1;
elseif(xdat==16)
    titx = '{\Sigma}_{\rm SFR}\:{\rm [M_{\odot}\:yr^{-1}\:kpc^{-2}]}';
    xl = [-4 -4];
    xu = [3 3];
    dx = 1;
elseif(xdat==17)
    titx = '$sSFR_{\rm c}\:{\rm [Gyr^{-1}]}';
    xl = [-3 -3];
    xu = [2 2];
    dx = 0.5;
elseif(xdat==18)
    titx = '{\tau}^{-1}\:{\rm [Gyr^{-1}]}';
    xl = [-3 -3];
    xu = [2 2];
    dx = 1;
elseif(xdat==19)
    titx = 'd\:/\:R_{\rm d}';
    xl = [-2 -2];
    xu = [0.5 0.5];
    dx = 0.5;
elseif(xdat==20)
    titx = 'h\:/\:H_{\rm d}';
    xl = [-3 -3];
    xu = [1 1];
    dx = 0.5;
elseif(xdat==21)
    titx = 'd\:{\rm [kpc]}';
    xl = [-0.5 -0.5];
    xu = [2 2];
    dx = 0.5;
elseif(xdat==22)
    titx = 'h\:{\rm [kpc]}';
    xl = [-2 -2];
    xu = [2 2];
    dx = 0.5;
elseif(xdat==23)
    titx ='\delta_{\rm c}';
    xl = [1 1];
    xu = [3 3];
    dx = 0.2;
elseif(xdat==24)
    titx = 'S_{\rm c}';
    xl = [-2 -2];
    xu = [0 0];
    dx = 0.2;
elseif(xdat==25)
    titx ='\delta_{\rm dm, \:c}';
    xl = [1 1];
    xu = [3 3];
    dx = 0.2;
elseif(xdat==28)
    titx ='t_{\rm ff}{\rm [Myr]}';
    xl = [0 0];
    xu = [2 2];
    dx = 0.2;
elseif(xdat==29)
    titx ='t_{\rm d, \: local}{\rm [Myr]}';
    xl = [0 0];
    xu = [2 2];
    dx = 0.2;
elseif(xdat==30)
    titx = '{\dot {M}}_{\rm in, \:Rc}\:{\rm [M_{\odot}\:yr^{-1}]}';
    xl = [-6 -6];
    xu = [2 2];
    dx = 1;
elseif(xdat==31)
    titx = '{\dot {M}}_{\rm in, \:1.5Rc}\:{\rm [M_{\odot}\:yr^{-1}]}';
    xl = [-6 -6];
    xu = [2 2];
    dx = 1;
elseif(xdat==32)
    titx = '{\dot {M}}_{\rm in, \:2Rc}\:{\rm [M_{\odot}\:yr^{-1}]}';
    xl = [-6 -6];
    xu = [2 2];
    dx = 1;
elseif(xdat==33)
    titx = '{\dot {M}}_{\rm out, \:Rc}\:{\rm [M_{\odot}\:yr^{-1}]}';
    xl = [-6 -6];
    xu = [2 2];
    dx = 1;
elseif(xdat==34)
    titx = '{\dot {M}}_{\rm out, \:1.5Rc}\:{\rm [M_{\odot}\:yr^{-1}]}';
    xl = [-6 -6];
    xu = [2 2];
    dx = 1;
elseif(xdat==35)
    titx = '{\dot {M}}_{\rm out, \:2Rc}\:{\rm [M_{\odot}\:yr^{-1}]}';
    xl = [-6 -6];
    xu = [2 2];
    dx = 1;
elseif(xdat==36)
    titx = '{\dot {M}}_{\rm in, stars}\:{\rm [M_{\odot}\:yr^{-1}]}';
    xl = [-6 -6];
    xu = [2 2];
    dx = 1;
elseif(xdat==37)
    titx = '{\dot {M}}_{\rm out, stars}\:{\rm [M_{\odot}\:yr^{-1}]}';
    xl = [-6 -6];
    xu = [2 2];
    dx = 1;
elseif(xdat==38)
    titx = '{\dot {M}}_{\rm formed, stars}\:{\rm [M_{\odot}\:yr^{-1}]}';
    xl = [-6 -6];
    xu = [2 2];
    dx = 1;
elseif(xdat==39)
    titx = 't\:{\rm [Myr]}';
    xl=[1 1];
    xu=[3.5 3.5];
    dx=0.5;
elseif(xdat==40)
    titx ='t_{\rm d, \: global}{\rm [Myr]}$';
    xl = [0 0];
    xu = [2 2];
    dx = 0.2;
elseif(xdat==41)
    titx = '{\dot {M}}_{\rm in, avg}\:{\rm [M_{\odot}\:yr^{-1}]}';
    xl = [-6 -6];
    xu = [2 2];
    dx = 1;
elseif(xdat==42)
    titx = '{\dot {M}}_{\rm out, avg}\:{\rm [M_{\odot}\:yr^{-1}]}';
    xl = [-6 -6];
    xu = [2 2];
    dx = 1;
elseif(xdat==43)
    titx = '\eta';
    xl = [-2 -2];
    xu = [2 2];
    dx = 0.5;
elseif(xdat==44)
    titx = '\epsilon_{\rm ff}';
    xl = [-3 -3];
    xu = [0 0];
    dx = 0.3;
elseif(xdat==45)
    titx = 't_{\rm d}\:/\:t_{\rm ff}';
    xl = [0 0];
    xu = [2 2];
    dx = 0.2;
elseif(xdat==46)
    titx = '\alpha{\rm )';
    xl = [-3 -3];
    xu = [0 0];
    dx = 0.3;
elseif(xdat==47)
    titx = '\eta_{\rm s}';
    xl = [-3 -3];
    xu = [4 4];
    dx = 1;
elseif(xdat==48)
    titx = '{\dot {M}}_{\rm out,\:gas} / M_{\rm c}\:{\rm [Gyr^{-1}]}';
    xl = [-3 -3];
    xu = [2 2];
    dx = 0.5;
elseif(xdat==49)
    titx = '{\dot {M}}_{\rm in,\:gas} / M_{\rm c}\:{\rm [Gyr^{-1}]}';
    xl = [-3 -3];
    xu = [2 2];
    dx = 0.5;
elseif(xdat==50)
    titx = '{\dot {M}}_{\rm out, \:stars, \:net} / M_{\rm c}\:{\rm [Gyr^{-1}]}';
    xl = [-3 -3];
    xu = [2 2];
    dx = 1;
elseif(xdat==51)
    titx = '{\dot {M}}_{\rm in, \:stars} / M_{\rm c}\:{\rm [Gyr^{-1}]}';
    xl = [-3 -3];
    xu = [2 2];
    dx = 1;
elseif(xdat==52)
    titx = '{\dot {M}}_{\rm out, \:stars} / M_{\rm c}\:{\rm [Gyr^{-1}]}';
    xl = [-3 -3];
    xu = [2 2];
    dx = 1;
end
if(strcmp(xaxis,'x log values'))
    if(xdat==13 | xdat==14)
        titx = strcat('$',titx);
    else
        titx = strcat('${\rm log(}',titx,'{\rm )}');
    end
elseif(strcmp(xaxis,'x linear scale'))
    titx = strcat('$',titx);
    if(xdat == 2)
        xl = [0.2 0.2];
        xu = [0.5 0.5];
        dx = 0.05;
    elseif(xdat==39)
        xl = [0 0];
        xu = [1200 1200];
        dx = 200;
    elseif(xdat==47 | xdat==50)
        xl = [-5 -5];
        xu = [5 5];
        dx = 1;
    else
        xl = [0 0];
        xu = 10.^xu;
        dx = xu(2)./10;
    end
elseif(strcmp(xaxis,'x log scale'))
    titx = strcat('$',titx);
    if(xdat == 2)
        xtick = [0.2:0.05:0.5];
        xl = [0.2 0.2];
        xu = [0.5 0.5];
    elseif(xdat==39)
        xtick = [10 100 1000];
        xl = [10 10];
        xu = [3000 3000];
    elseif(xdat==47 | xdat==50)
        xtick = [0.001 0.01 0.1 1 10 100];
        xl = [0.001 0.001];
        xu = [100 100];
    else
        xtick = 10.^[xl:xu];
        xl = 10.^xl;
        xu = 10.^xu;
    end
    dx = (xu(2)-xl(2))/10;
end
if(strcmp(xnorm,'x normalized'))
    if(xdat==39)
        titx = '$t\:/\:t_{\rm d}$';
        if(strcmp(xaxis,'x log values'))
            xl = [-1 -1];
            xu = [2.2 2.2];
            dx = 0.4;
        elseif(strcmp(xaxis,'x linear scale'))
            xl = [0 0];
            xu = [50 50];
            dx = 5;
        elseif(strcmp(xaxis,'x log scale'))
            xl = [0.1 0.1];
            xu = [20 20];
            xtick = [0.1 1 10 100];
        end
    else
        titx = strcat(titx,'\:{\rm (normalized)}$');
        if(strcmp(xaxis,'x log values'))
            xl = [-1 -1];
            xu = [1 1];
            dx = 0.2;
        elseif(strcmp(xaxis,'x linear scale'))
            xl = [0 0];
            xu = [10 10];
            dx = 1;
        elseif(strcmp(xaxis,'x log scale'))
            xl = [0.1 0.1];
            xu = [10 10];
            xtick = [0.1 1 10];
        end
    end
elseif(strcmp(xnorm,'x absolute'))
    titx = strcat(titx,'$');
end

if(ydat==2)
    tity = 'a';
    yl = [-1 -1];
    yu = [0 0];
    dy = 0.1;
elseif(ydat==3)
    tity = 'R_{\rm c}\:{\rm [kpc]}';
    yl = [-1 -1];
    yu = [0 0];
    dy = 0.1;
elseif(ydat==4)
    tity = 'M_{\rm c,\:gas}\:{\rm [M_{\odot}]}';
    yl = [2 2];
    yu = [9 9];
    dy = 1;
elseif(ydat==5)
    tity = 'M_{\rm c,\:star}\:{\rm [M_{\odot}]}';
    yl = [2 2];
    yu = [9 9];
    dy = 1;
elseif(ydat==6)
    %tity = 'M_{\rm c,\:bar}\:{\rm [M_{\odot}]}';
    tity = 'M_{\rm c,\:bar}';
    yl = [6 6];
    yu = [10 10];
    dy = 1;
elseif(ydat==7)
    tity = 'f_{\rm gas}';
    yl = [-6 -6];
    yu = [0 0];
    dy = 1;
elseif(ydat==8)
    tity = 'f_{\rm dm}';
    yl = [-6 -6];
    yu = [0 0];
    dy = 1;
elseif(ydat==9)
    tity = '\Sigma_{\rm c,\:gas}\:{\rm [M_{\odot}\:pc^{-2}]}';
    yl = [0 0];
    yu = [9 9];
    dy = 1;
elseif(ydat==10)
    tity = '\Sigma_{\rm c,\:star}\:{\rm [M_{\odot}\:pc^{-2}]}';
    yl = [2 2];
    yu = [9 9];
    dy = 1;
elseif(ydat==11)
    tity = '\Sigma_{\rm c,\:bar}\:{\rm [M_{\odot}\:pc^{-2}]}';
    yl = [6 6];
    yu = [10 10];
    dy = 1;
elseif(ydat==12)
    tity = 'age\:{\rm [Myr]}';
    yl = [0.5 0.5];
    yu = [4 4];
    dy = 0.5;
elseif(ydat==13)
    tity = '${\rm log(O\:/\:H)_{gas}}\:+\:12';
    yl = [7.5 7.5];
    yu = [9.5 9.5];
    dy = 0.2;
elseif(ydat==14)
    tity = '${\rm log(O\:/\:H)_{star}}\:+\:12';
    yl = [7.5 7.5];
    yu = [9.5 9.5];
    dy = 0.2;
elseif(ydat==15)
    tity = 'SFR_{\rm c}\:{\rm [M_{\odot}\:yr^{-1}]}';
    yl = [-6 -6];
    yu = [2 2];
    dy = 1;
elseif(ydat==16)
    tity = '{\Sigma}_{\rm SFR}\:{\rm [M_{\odot}\:yr^{-1}\:kpc^{-2}]}';
    yl = [-4 -4];
    yu = [3 3];
    dy = 1;
elseif(ydat==17)
    tity = '$sSFR_{\rm c}\:{\rm [Gyr^{-1}]}';
    yl = [-3 -3];
    yu = [2 2];
    dy = 0.5;
elseif(ydat==18)
    tity = '{\tau}^{-1}\:{\rm [Gyr^{-1}]}';
    yl = [-3 -3];
    yu = [2 2];
    dy = 1;
elseif(ydat==19)
    tity = 'd\:/\:R_{\rm d}';
    yl = [-2 -2];
    yu = [0.5 0.5];
    dy = 0.5;
elseif(ydat==20)
    tity = 'h\:/\:H_{\rm d}';
    yl = [-3 -3];
    yu = [1 1];
    dy = 0.5;
elseif(ydat==21)
    tity = 'd\:{\rm [kpc]}';
    yl = [-0.5 -0.5];
    yu = [2 2];
    dy = 0.5;
elseif(ydat==22)
    tity = 'h\:{\rm [kpc]}';
    yl = [-2 -2];
    yu = [2 2];
    dy = 0.5;
elseif(ydat==23)
    tity ='\delta_{\rm c}';
    yl = [1 1];
    yu = [3 3];
    dy = 0.2;
elseif(ydat==24)
    tity = 'S_{\rm c}';
    yl = [-2 -2];
    yu = [0 0];
    dy = 0.2;
elseif(ydat==25)
    tity ='\delta_{\rm dm, \:c}';
    yl = [1 1];
    yu = [3 3];
    dy = 0.2;
elseif(ydat==28)
    tity ='t_{\rm ff}{\rm [Myr]}';
    yl = [0 0];
    yu = [2 2];
    dy = 0.2;
elseif(ydat==29)
    tity ='t_{\rm d, \: local}{\rm [Myr]}';
    yl = [0 0];
    yu = [2 2];
    dy = 0.2;
elseif(ydat==30)
    tity = '{\dot {M}}_{\rm in, \:Rc}\:{\rm [M_{\odot}\:yr^{-1}]}';
    yl = [-6 -6];
    yu = [2 2];
    dy = 1;
elseif(ydat==31)
    tity = '{\dot {M}}_{\rm in, \:1.5Rc}\:{\rm [M_{\odot}\:yr^{-1}]}';
    yl = [-6 -6];
    yu = [2 2];
    dy = 1;
elseif(ydat==32)
    tity = '{\dot {M}}_{\rm in, \:2Rc}\:{\rm [M_{\odot}\:yr^{-1}]}';
    yl = [-6 -6];
    yu = [2 2];
    dy = 1;
elseif(ydat==33)
    tity = '{\dot {M}}_{\rm out, \:Rc}\:{\rm [M_{\odot}\:yr^{-1}]}';
    yl = [-6 -6];
    yu = [2 2];
    dy = 1;
elseif(ydat==34)
    tity = '{\dot {M}}_{\rm out, \:1.5Rc}\:{\rm [M_{\odot}\:yr^{-1}]}';
    yl = [-6 -6];
    yu = [2 2];
    dy = 1;
elseif(ydat==35)
    tity = '{\dot {M}}_{\rm out, \:2Rc}\:{\rm [M_{\odot}\:yr^{-1}]}';
    yl = [-6 -6];
    yu = [2 2];
    dy = 1;
elseif(ydat==36)
    tity = '{\dot {M}}_{\rm in, stars}\:{\rm [M_{\odot}\:yr^{-1}]}';
    yl = [-6 -6];
    yu = [2 2];
    dy = 1;
elseif(ydat==37)
    tity = '{\dot {M}}_{\rm out, stars}\:{\rm [M_{\odot}\:yr^{-1}]}';
    yl = [-6 -6];
    yu = [2 2];
    dy = 1;
elseif(ydat==38)
    tity = '{\dot {M}}_{\rm formed, stars}\:{\rm [M_{\odot}\:yr^{-1}]}';
    yl = [-6 -6];
    yu = [2 2];
    dy = 1;
elseif(ydat==39)
    tity = 't\:{\rm [Myr]}';
    yl=[1 1];
    yu=[3.5 3.5];
    dy=0.5;
elseif(ydat==40)
    tity ='t_{\rm d, \: global}{\rm [Myr]}$';
    yl = [0 0];
    yu = [2 2];
    dy = 0.2;
elseif(ydat==41)
    tity = '{\dot {M}}_{\rm in, avg}\:{\rm [M_{\odot}\:yr^{-1}]}';
    yl = [-6 -6];
    yu = [2 2];
    dy = 1;
elseif(ydat==42)
    tity = '{\dot {M}}_{\rm out, avg}\:{\rm [M_{\odot}\:yr^{-1}]}';
    yl = [-6 -6];
    yu = [2 2];
    dy = 1;
elseif(ydat==43)
    tity = '\eta';
    yl = [-2 -2];
    yu = [2 2];
    dy = 0.5;
elseif(ydat==44)
    tity = '\epsilon_{\rm ff}';
    yl = [-3 -3];
    yu = [0 0];
    dy = 0.3;
elseif(ydat==45)
    tity = 't_{\rm d}\:/\:t_{\rm ff}';
    yl = [0 0];
    yu = [2 2];
    dy = 0.2;
elseif(ydat==46)
    tity = '\alpha';
    yl = [-2 -2];
    yu = [1 1];
    dy = 0.3;
elseif(ydat==47)
    tity = '\eta_{\rm s}';
    yl = [-3 -3];
    yu = [4 4];
    dy = 1;
elseif(ydat==48)
    tity = '{\dot {M}}_{\rm out,\:gas} / M_{\rm c}\:{\rm [Gyr^{-1}]}';
    yl = [-3 -3];
    yu = [2 2];
    dy = 0.5;
elseif(ydat==49)
    tity = '{\dot {M}}_{\rm in,\:gas} / M_{\rm c}\:{\rm [Gyr^{-1}]}';
    yl = [-3 -3];
    yu = [2 2];
    dy = 0.5;
elseif(ydat==50)
    tity = '{\dot {M}}_{\rm out, \:stars, \:net} / M_{\rm c}\:{\rm [Gyr^{-1}]}';
    yl = [-3 -3];
    yu = [2 2];
    dy = 1;
elseif(ydat==51)
    tity = '{\dot {M}}_{\rm in, \:stars} / M_{\rm c}\:{\rm [Gyr^{-1}]}';
    yl = [-3 -3];
    yu = [2 2];
    dy = 1;
elseif(ydat==52)
    tity = '{\dot {M}}_{\rm out, \:stars} / M_{\rm c}\:{\rm [Gyr^{-1}]}';
    yl = [-3 -3];
    yu = [2 2];
    dy = 1;
end
if(strcmp(yaxis,'y log values'))
    if(ydat==13 | ydat==14)
        tity = strcat('$',tity);
    else
        tity = strcat('${\rm log(}',tity,'{\rm )}');
    end
elseif(strcmp(yaxis,'y linear scale'))
    tity = strcat('$',tity);
    if(ydat == 2)
        yl = [0.2 0.2];
        yu = [0.5 0.5];
        dy = 0.05;
    elseif(ydat==39)
        yl = [0 0];
        yu = [2400 2400];
        dy = 400;
    elseif(ydat==47 | ydat==50)
        yl = [-5 -5];
        yu = [5 5];
        dy = 1;
    else
        yl = [0 0];
        yu = 10.^yu;
        dy = yu(2)./10;
    end
elseif(strcmp(yaxis,'y log scale'))
    tity = strcat('$',tity);
    if(ydat == 2)
        ytick = [0.2:0.05:0.5];
        yl = [0.2 0.2];
        yu = [0.5 0.5];
    elseif(ydat==39)
        ytick = [10 100 1000];
        yl = [10 10];
        yu = [3000 3000];
    elseif(ydat==47 | ydat==50)
        ytick = [0.001 0.01 0.1 1 10 100];
        yl = [0.001 0.001];
        yu = [100 100];
    else
        ytick = 10.^[yl:yu];
        yl = 10.^yl;
        yu = 10.^yu;
    end
    dy = (yu(2)-yl(2))/10;
end
if(strcmp(ynorm,'y normalized'))
    tity = strcat(tity,'\:{\rm (normalized)}$');
    if(strcmp(yaxis,'y log values'))
        if(ydat==7)
            yl = [-3 -3];
            dy = 0.5;
        else
            yl = [-1 -1];
            dy = 0.2;
        end
        yu = [1 1];
    elseif(strcmp(yaxis,'y linear scale'))
        yl = [0 0];
        yu = [10 10];
        dy = 1;
    elseif(strcmp(yaxis,'y log scale'))
        if(ydat==7)
            yl = [0.001 0.001];
            ytick = [0.001 0.01 0.1 1 10];
        else
            yl = [0.1 0.1];
            ytick = [0.1 1 10];
        end
        yu = [10 10];
    end
elseif(strcmp(ynorm,'y absolute'))
    tity = strcat(tity,'$');
end

figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
set(figure1,'WindowStyle','docked');

% Create axes
axes1 = axes('Parent',figure1,'YTick',[yl(1):dy:yu(1)],...
    'YMinorTick','on',...
    'XTick',[xl(1):dx:xu(1)],...
    'XMinorTick','on',...
    'TickLength',[0.025 0.05],...
    'Position',[0.12 0.12 0.86 0.82],...
    'PlotBoxAspectRatio',[1 1 1],...
    'LineWidth',1.5,...
    'FontSize',16,...
    'FontName','Arial');

if(strcmp(yaxis,'y log scale'))
    set(gca,'YScale','log','YTick',ytick);
end
if(strcmp(xaxis,'x log scale'))
    set(gca,'XScale','log','XTick',xtick);
end

xlim(axes1,[xl(2) xu(2)]);
ylim(axes1,[yl(2) yu(2)]);
box(axes1,'on');
hold(axes1,'all');
%grid on

% Create xlabel
xlabel(titx,'Interpreter','latex','FontSize',24,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 -0.06 0]);

% Create ylabel
ylabel(tity,'Interpreter','latex','FontSize',24,...
    'FontName','Times New Roman','units','normalized',...
    'position',[-0.06 0.515 0]);

b2 = begin2;
e2 = end2;
dat2 = is2;
dat2(:,39) = 1000.*(dat2(:,39));
aexp2 = 1./(1+dat2(b2,2));
t2 = 950./( ((1./aexp2)./7).^1.5 ) - 950./( ((1./(aexp2-0.01))./7).^1.5 );
dat2(:,2) = 1./(1+dat2(:,2));

b3 = begin3;
e3 = end3;
dat3 = is3;
dat3(:,39) = 1000.*(dat3(:,39));
aexp3 = 1./(1+dat3(b3,2));
t3 = 950./( ((1./aexp3)./7).^1.5 ) - 950./( ((1./(aexp3-0.01))./7).^1.5 );
dat3(:,2) = 1./(1+dat3(:,2));

%dat2(b2,39) = 0.5.*t2;
dat2(b2,39) = min(t2,dat2(b2,12));
%dat3(b3,39) = 0.5.*t3;
dat3(b3,39) = min(t3,dat3(b3,12));

for i=1:length(b2)
    dat2(b2(i)+1:e2(i),39) = dat2(b2(i)+1:e2(i),39) + dat2(b2(i),39);    
end
dat2(:,41) = (dat2(:,30)+dat2(:,31)+dat2(:,32))./3;
dat2(:,42) = (dat2(:,33)+dat2(:,34)+dat2(:,35))./3;
for i=1:length(dat2(:,17))
    if(dat2(i,17)>0.01)
        dat2(i,43) = dat2(i,42)./dat2(i,15);
        dat2(i,44) = abs(1e6.*dat2(i,28).*dat2(i,15)./dat2(i,4));
    else
        dat2(i,43) = 0;
        dat2(i,44) = 0;
    end
end
%dat2(:,45) = abs(dat2(:,29)./dat2(:,28));       %% local t_d
dat2(:,45) = abs(dat2(:,40)./dat2(:,28));       %% global t_d
dat2(:,46:52) = 0;
del2 = dat2(:,41);
for i=1:length(b2)
    dat2(b2(i),46) = ( 2 .* 1e6 .* abs(dat2(b2(i),40)) ./ dat2(b2(i),6) ) .* dat2(b2(i),41);
    for j = (b2(i)+1):e2(i)
        del = (dat2(j,6) - dat2(j-1,6)) ./ ( (dat2(j,39) - dat2(j-1,39)) .* 1e9 ) - ...
            ( dat2(j,36) - dat2(j,37) - 0.5.*(dat2(j,42) + dat2(j-1,42)) );
        if(del>0)
            del2(j) = del;
        end
        dat2(j,46) = ( 2 .* 1e6 .* ( 0.5.*(abs(dat2(j,40)) + abs(dat2(j-1,40))) ) ./ ...
            ( 0.5.*(dat2(j,6) + dat2(j-1,6)) ) ) .* del2(j);
        if( dat2(j,17) > 0.01 )
            dat2(j,47) = ( dat2(j,37) - dat2(j,36) ) ./ dat2(j,38);
        end
    end 
end
dat2(:,48) = 1e9 .* dat2(:,42)./dat2(:,6);
%dat2(:,49) = 1e9 .* dat2(:,41)./dat2(:,6);
dat2(:,49) = 1e9 .* del2(:)./dat2(:,6);
for i=1:length(b2)
    dat2(b2(i)+1:e2(i),50) = 1e9 .* (dat2(b2(i)+1:e2(i),37) - dat2(b2(i)+1:e2(i),36)) ./ ...
        (dat2(b2(i)+1:e2(i),6));
    dat2(b2(i)+1:e2(i),51) = 1e9 .* dat2(b2(i)+1:e2(i),36) ./ dat2(b2(i)+1:e2(i),6);
    dat2(b2(i)+1:e2(i),52) = 1e9 .* dat2(b2(i)+1:e2(i),37) ./ dat2(b2(i)+1:e2(i),6);
end

for i=1:length(b3)
    dat3(b3(i)+1:e3(i),39) = dat3(b3(i)+1:e3(i),39) + dat3(b3(i),39);
end
dat3(:,41) = (dat3(:,30)+dat3(:,31)+dat3(:,32))./3;
dat3(:,42) = (dat3(:,33)+dat3(:,34)+dat3(:,35))./3;
for i=1:length(dat3(:,17))
    if(dat3(i,17)>0.01)
        dat3(i,43) = dat3(i,42)./dat3(i,15);
        dat3(i,44) = abs(1e6.*dat3(i,28).*dat3(i,15)./dat3(i,4));
    else
        dat3(i,43) = 0;
        dat3(i,44) = 0;
    end
end
%dat3(:,45) = abs(dat3(:,29)./dat3(:,28));      %% local t_d
dat3(:,45) = abs(dat3(:,40)./dat3(:,28));       %% global t_d
dat3(:,46:52) = 0;
del3 = dat3(:,41);
for i=1:length(b3)
    dat3(b3(i),46) = ( 2 .* 1e6 .* abs(dat3(b3(i),40)) ./ dat3(b3(i),6) ) .* dat3(b3(i),41);
    for j = (b3(i)+1):e3(i)
        del = (dat3(j,6) - dat3(j-1,6)) ./ ( (dat3(j,39) - dat3(j-1,39)) .* 1e9 ) - ...
            ( dat3(j,36) - dat3(j,37) - 0.5.*(dat3(j,42) + dat3(j-1,42)) );
        if(del>0)
            del3(j) = del;
        end
        dat3(j,46) = ( 2 .* 1e6 .* ( 0.5.*(abs(dat3(j,40)) + abs(dat3(j-1,40))) ) ./ ...
            ( 0.5.*(dat3(j,6) + dat3(j-1,6)) ) ) .* del3(j);
        if( dat3(j,17) > 0.01 )
            dat3(j,47) = ( dat3(j,37) - dat3(j,36) ) ./ dat3(j,38);
        end
    end 
end
dat3(:,48) = 1e9 .* dat3(:,42)./dat3(:,6);
%dat3(:,49) = 1e9 .* dat3(:,41)./dat3(:,6);
dat3(:,49) = 1e9 .* del3(:)./dat3(:,6);
for i=1:length(b3)
    dat3(b3(i)+1:e3(i),50) = 1e9 .* (dat3(b3(i)+1:e3(i),37) - dat3(b3(i)+1:e3(i),36)) ./ ...
        (dat3(b3(i)+1:e3(i),6));
    dat3(b3(i)+1:e3(i),51) = 1e9 .* dat3(b3(i)+1:e3(i),36) ./ dat3(b3(i)+1:e3(i),6);
    dat3(b3(i)+1:e3(i),52) = 1e9 .* dat3(b3(i)+1:e3(i),37) ./ dat3(b3(i)+1:e3(i),6);
end

nmigrate2 = 0;
migrate2 = [];
for i=1:length(b2)
    if(e2(i)-b2(i)==1)
        cond = dat2(e2(i),21)./dat2(b2(i),21);
    elseif(e2(i)-b2(i)>1)
        cond = mean(dat2(e2(i)-1:e2(i),21))./mean(dat2(b2(i):b2(i)+1,21));
    end
    if(cond > dist_thresh(1) & cond <= dist_thresh(2) & ...
       max(abs(dat2(b2(i):e2(i),20))) > height_thresh(1) & ...
       max(abs(dat2(b2(i):e2(i),20))) <= height_thresh(2) & ...
       max(dat2(b2(i):e2(i),39)./dat2(b2(i):e2(i),40)) > time_thresh(1) & ...
       max(dat2(b2(i):e2(i),39)./dat2(b2(i):e2(i),40)) <= time_thresh(2))
        nmigrate2 = nmigrate2 + 1;
        migrate2(nmigrate2) = i;
    end
end
nmigrate3 = 0;
migrate3 = [];
for i=1:length(b3)
    if(e3(i)-b3(i)==1)
        cond = dat3(e3(i),21)./dat3(b3(i),21);
    elseif(e3(i)-b3(i)>1)
        cond = mean(dat3(e3(i)-1:e3(i),21))./mean(dat3(b3(i):b3(i)+1,21));
    end
    if(cond > dist_thresh(1) & cond <= dist_thresh(2) & ...
       max(abs(dat3(b3(i):e3(i),20))) > height_thresh(1) & ...
       max(abs(dat3(b3(i):e3(i),20))) <= height_thresh(2) & ...
       max(dat3(b3(i):e3(i),39)./dat3(b3(i):e3(i),40)) > time_thresh(1) & ...
       max(dat3(b3(i):e3(i),39)./dat3(b3(i):e3(i),40)) <= time_thresh(2))
        nmigrate3 = nmigrate3 + 1;
        migrate3(nmigrate3) = i;
    end
end            
            
if(strcmp(stack,'no stack'))
    if(res_thresh(1)==0 & shape_thresh(1)==0)
        if(mass_thresh(2)>=9.5)
            tit = strcat('$No \: \: RP \: \: \: ','{\rm log(}M_{\rm c}{\rm )}>',...
                num2str(mass_thresh(1),'%3.1f'),'$');
        else
            tit = strcat('$No \: \: RP \: \: \: ',num2str(mass_thresh(2),'%3.1f'),...
                '>{\rm log(}M_{\rm c}{\rm )}>',num2str(mass_thresh(1),'%3.1f'),'$');
        end
    else
        if(mass_thresh(2)>=9.5)
            tit = strcat('$No \: RP \: \: \: ',...
                '{\rm log(}M_{\rm c}{\rm )}>',num2str(mass_thresh(1),'%3.1f'),'\: \: \:',...
                'S_{\rm c}>',num2str(shape_thresh(1),'%3.1f'),'\: \: \:',...
                '{\rm log(}\delta_{\rm c}{\rm )}>',num2str(res_thresh(1),'%3.1f'),'$');
        else
            tit = strcat('$No \: RP \: \: \: ',num2str(mass_thresh(2),'%3.1f'),...
                '>{\rm log(}M_{\rm c}{\rm )}>',num2str(mass_thresh(1),'%3.1f'),'\: \: \:',...
                'S_{\rm c}>',num2str(shape_thresh(1),'%3.1f'),'\: \: \:',...
                '{\rm log(}\delta_{\rm c}{\rm )}>',num2str(res_thresh(1),'%3.1f'),'$');
        end
    end
    title(tit,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 1.005 0]);%,'BackgroundColor',[.7 .9 .7]);
    colOrd = get(gca, 'ColorOrder');
    for i=1:nmigrate2
        n = migrate2(i);
        j = mod(i,length(colOrd))+1;
        normx = mean(dat2(b2(n):b2(n)+1,xdat));
        normy = mean(dat2(b2(n):b2(n)+1,ydat));
%        normx = dat2(b2(n),xdat);
%        normy = dat2(b2(n),ydat);
        if(strcmp(xnorm,'x normalized'))
            if(strcmp(xaxis,'x linear scale') | strcmp(xaxis,'x log scale') | xdat==13 | xdat==14)
                xp = dat2(b2(n):e2(n),xdat) ./ normx;
                if(xdat == 39)
                    xp = dat2(b2(n):e2(n),39) ./ dat2(b2(n):e2(n),40);
                end
            else
                xp = log10( dat2(b2(n):e2(n),xdat) ./ normx );
                if(xdat == 39)
                    xp = log10( dat2(b2(n):e2(n),39) ./ dat2(b2(n):e2(n),40) );
                end
            end
        elseif(strcmp(xnorm,'x absolute'))
            if(strcmp(xaxis,'x linear scale') | strcmp(xaxis,'x log scale') | xdat==13 | xdat==14)
                xp = dat2(b2(n):e2(n),xdat);
            else
                xp = log10( dat2(b2(n):e2(n),xdat) );
            end
        end
        if(strcmp(ynorm,'y normalized'))
            if(strcmp(yaxis,'y linear scale') | strcmp(yaxis,'y log scale') | ydat==13 | ydat==14)
                yp = dat2(b2(n):e2(n),ydat) ./normy;
            else
                yp = log10( dat2(b2(n):e2(n),ydat) ./ normy );
            end
        elseif(strcmp(ynorm,'y absolute'))
            if(strcmp(yaxis,'y linear scale') | strcmp(yaxis,'y log scale') | ydat==13 | ydat==14)
                yp = dat2(b2(n):e2(n),ydat);
            else
                yp = log10( dat2(b2(n):e2(n),ydat) );
            end
        end
        plot(xp, yp, ...
            'linestyle','-','linewidth',1,'color',colOrd(j,:),...
            'marker','o','markersize',3,...
            'markerfacecolor',colOrd(j,:),'markeredgecolor',colOrd(j,:));
    end
    if(max(xdat,ydat)==39 & min(xdat,ydat)==12)
        if(strcmp(xaxis,'x log values'))
            xlim([1 4]);
            limx = [1 4];
            xtick = [1:0.5:4];
        elseif(strcmp(xaxis,'x linear scale'))
            xlim([0 3e3]);
            limx = [0 3e3];
            xtick = [0:100:1000];            
        elseif(strcmp(xaxis,'x log scale'))
            xlim([10 5e3]);
            limx = [10 5e3];
            xtick = [10 100 1000 5000];
        end
        if(strcmp(yaxis,'y log values'))
            ylim([1 4]);
            limy = [1 4];
            ytick = [1:0.5:4];
        elseif(strcmp(yaxis,'y linear scale'))
            ylim([0 3e3]);
            limy = [0 3e3];
            ytick = [0:100:1000];            
        elseif(strcmp(yaxis,'y log scale'))
            ylim([10 5e3]);
            limy = [10 5e3];
            ytick = [10 100 1000 5000];
        end
        set(gca,'XTick',xtick,'YTick',ytick);
        p = linspace(limx(1),limx(2),10);
        plot(p,p,'linewidth',3,'linestyle','--','color','k','marker','none')
    end

    figure2 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
    set(figure2,'WindowStyle','docked');
    axes2 = axes('Parent',figure2,'YTick',[yl(1):dy:yu(1)],...
        'YMinorTick','on',...
        'XTick',[xl(1):dx:xu(1)],...
        'XMinorTick','on',...
        'TickLength',[0.025 0.05],...
        'Position',[0.12 0.12 0.86 0.82],...
        'PlotBoxAspectRatio',[1 1 1],...
        'LineWidth',1.5,...
        'FontSize',16,...
        'FontName','Arial');
    if(strcmp(yaxis,'y log scale'))
        set(gca,'YScale','log','YTick',ytick);
    end
    if(strcmp(xaxis,'x log scale'))
        set(gca,'XScale','log','XTick',xtick);
    end
    
    xlim(axes2,[xl(2) xu(2)]);
    ylim(axes2,[yl(2) yu(2)]);
    box(axes2,'on');
    hold(axes2,'all');
    %grid on
    xlabel(titx,'Interpreter','latex','FontSize',24,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.06 0]);
    ylabel(tity,'Interpreter','latex','FontSize',24,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.06 0.515 0]);
    
    if(res_thresh(1)==0 & shape_thresh(1)==0)
        if(mass_thresh(2)>=9.5)
            tit = strcat('$RP \: \: \: ',...
                '{\rm log(}M_{\rm c}{\rm )}>',num2str(mass_thresh(1),'%3.1f'),'$');
        else
            tit = strcat('$RP \: \: \: ',num2str(mass_thresh(2),'%3.1f'),...
                '>{\rm log(}M_{\rm c}{\rm )}>',num2str(mass_thresh(1),'%3.1f'),'$');
        end
    else
        if(mass_thresh(2)>=9.5)
            tit = strcat('$RP \: \: \: ',...
                '{\rm log(}M_{\rm c}{\rm )}>',num2str(mass_thresh(1),'%3.1f'),'\: \: \:',...
                'S_{\rm c}>',num2str(shape_thresh(1),'%3.1f'),'\: \: \:',...
                '{\rm log(}\delta_{\rm c}{\rm )}>',num2str(res_thresh(1),'%3.1f'),'$');
        else
            tit = strcat('$RP \: \: \: ',num2str(mass_thresh(2),'%3.1f'),...
                '>{\rm log(}M_{\rm c}{\rm )}>',num2str(mass_thresh(1),'%3.1f'),'\: \: \:',...
                'S_{\rm c}>',num2str(shape_thresh(1),'%3.1f'),'\: \: \:',...
                '{\rm log(}\delta_{\rm c}{\rm )}>',num2str(res_thresh(1),'%3.1f'),'$');
        end
    end
    title(tit,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 1.005 0]);%,'BackgroundColor',[.7 .9 .7]);
    colOrd = get(gca, 'ColorOrder');
    for i=1:nmigrate3
        n = migrate3(i);
        j = mod(i,length(colOrd))+1;
        normx = mean(dat3(b3(n):b3(n)+1,xdat));
        normy = mean(dat3(b3(n):b3(n)+1,ydat));
%         normx = dat3(b3(n),xdat);
%         normy = dat3(b3(n),ydat);
        if(strcmp(xnorm,'x normalized'))
            if(strcmp(xaxis,'x linear scale') | strcmp(xaxis,'x log scale') | xdat==13 | xdat==14)
                xp = dat3(b3(n):e3(n),xdat) ./ normx;
                if(xdat == 39)
                    xp = dat3(b3(n):e3(n),39) ./ dat3(b3(n):e3(n),40);
                end
            else
                xp = log10( dat3(b3(n):e3(n),xdat) ./ normx );
                if(xdat == 39)
                    xp = log10( dat3(b3(n):e3(n),39) ./ dat3(b3(n):e3(n),40) );
                end
            end
        elseif(strcmp(xnorm,'x absolute'))
            if(strcmp(xaxis,'x linear scale') | strcmp(xaxis,'x log scale') | xdat==13 | xdat==14)
                xp = dat3(b3(n):e3(n),xdat);
            else
                xp = log10( dat3(b3(n):e3(n),xdat) );
            end
        end
        if(strcmp(ynorm,'y normalized'))
            if(strcmp(yaxis,'y linear scale') | strcmp(yaxis,'y log scale') | ydat==13 | ydat==14)
                yp = dat3(b3(n):e3(n),ydat) ./ normy;
            else
                yp = log10( dat3(b3(n):e3(n),ydat) ./ normy );
            end
        elseif(strcmp(ynorm,'y absolute'))
            if(strcmp(yaxis,'y linear scale') | strcmp(yaxis,'y log scale') | ydat==13 | ydat==14)
                yp = dat3(b3(n):e3(n),ydat);
            else
                yp = log10( dat3(b3(n):e3(n),ydat) );
            end
        end
        plot(xp, yp, ...
            'linestyle','-','linewidth',1,'color',colOrd(j,:),...
            'marker','o','markersize',3,...
            'markerfacecolor',colOrd(j,:),'markeredgecolor',colOrd(j,:));
    end
    if(max(xdat,ydat)==39 & min(xdat,ydat)==12)
        if(strcmp(xaxis,'x log values'))
            xlim([1 4]);
            xtick = [1:0.5:4];
        elseif(strcmp(xaxis,'x linear scale'))
            xlim([0 3e3]);
            xtick = [0:100:1000];            
        elseif(strcmp(xaxis,'x log scale'))
            xlim([10 5e3]);
            xtick = [10 100 1000 5000];
        end
        if(strcmp(yaxis,'y log values'))
            ylim([1 4]);
            ytick = [1:0.5:4];
        elseif(strcmp(yaxis,'y linear scale'))
            ylim([0 3e3]);
            ytick = [0:100:1000];            
        elseif(strcmp(yaxis,'y log scale'))
            ylim([10 5e3]);
            ytick = [10 100 1000 5000];
        end
        set(gca,'XTick',xtick,'YTick',ytick);
        p = linspace(limx(1),limx(2),10);
        plot(p,p,'linewidth',3,'linestyle','--','color','k','marker','none')
    end
    [length(b2),nmigrate2,length(b3),nmigrate3]
    
else
    if(res_thresh(1)==0 & shape_thresh(1)==0)
        if(mass_thresh(2)>=9.5)
            tit = strcat('$',...
                '{\rm log(}M_{\rm c}{\rm )}>',num2str(mass_thresh(1),'%3.1f'),'$');
        else
            tit = strcat('$',num2str(mass_thresh(2),'%3.1f'),...
                '>{\rm log(}M_{\rm c}{\rm )}>',num2str(mass_thresh(1),'%3.1f'),'$');
        end
    else
        if(mass_thresh(2)>=9.5)
            tit = strcat('$',...
                '{\rm log(}M_{\rm c}{\rm )}>',num2str(mass_thresh(1),'%3.1f'),'\: \: \:',...
                'S_{\rm c}>',num2str(shape_thresh(1),'%3.1f'),'\: \: \:',...
                '{\rm log(}\delta_{\rm c}{\rm )}>',num2str(res_thresh(1),'%3.1f'),'$');
        else
            tit = strcat('$',num2str(mass_thresh(2),'%3.1f'),...
                '>{\rm log(}M_{\rm c}{\rm )}>',num2str(mass_thresh(1),'%3.1f'),'\: \: \:',...
                'S_{\rm c}>',num2str(shape_thresh(1),'%3.1f'),'\: \: \:',...
                '{\rm log(}\delta_{\rm c}{\rm )}>',num2str(res_thresh(1),'%3.1f'),'$');
        end
    end
    title(tit,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 1.005 0]);%,'BackgroundColor',[.7 .9 .7]);
    if(xdat==39 & strcmp(xnorm,'x absolute'))
        %bin = log10([10,60,150:150:3000]);
        %bin = log10([10,60,150:150:900,1200:300:3000]);
        bin = log10([0,120,250,400,550,700,850,1000]);%,10.^[3.2:0.2:3.6]]);
        %bin = log10([0:200:1200]);
    elseif(xdat==39 & strcmp(xnorm,'x normalized'))
        bin = log10([0:2.5:10,15:5:30,40:10:100]);
    elseif(xdat==19)
        bin=[-2,-1:0.2:0.4];
    else
        bin = [xl(2):0.5.*dx:xu(2)];
    end
    
    x2 = zeros(length(bin)-1,1);
    x3 = zeros(length(bin)-1,1);
    y2 = zeros(length(bin)-1,1);
    y3 = zeros(length(bin)-1,1);
    s2 = zeros(length(bin)-1,1);
    s3 = zeros(length(bin)-1,1);
    n2 = zeros(length(bin)-1,1);
    n3 = zeros(length(bin)-1,1);
    xmed2 = zeros(length(bin)-1,1);
    xmed3 = zeros(length(bin)-1,1);
    ymed2 = zeros(length(bin)-1,1);
    ymed3 = zeros(length(bin)-1,1);
    xl2 = zeros(length(bin)-1,1);
    xl3 = zeros(length(bin)-1,1);
    xu2 = zeros(length(bin)-1,1);
    xu3 = zeros(length(bin)-1,1);    
    yl2 = zeros(length(bin)-1,1);
    yl3 = zeros(length(bin)-1,1);
    yu2 = zeros(length(bin)-1,1);
    yu3 = zeros(length(bin)-1,1);
    sortedx2 = zeros(length(bin)-1, sum(e2(1:nmigrate2)-b2(1:nmigrate2)+nmigrate2));
    sortedy2 = zeros(length(bin)-1, sum(e2(1:nmigrate2)-b2(1:nmigrate2)+nmigrate2));
    sortedx3 = zeros(length(bin)-1, sum(e3(1:nmigrate3)-b3(1:nmigrate3)+nmigrate3));
    sortedy3 = zeros(length(bin)-1, sum(e3(1:nmigrate3)-b3(1:nmigrate3)+nmigrate3));
    
    for j=1:nmigrate2
        n = migrate2(j);
        normx = mean(dat2(b2(n):b2(n)+1,xdat));
        normy = mean(dat2(b2(n):b2(n)+1,ydat));
%         normx = dat2(b2(n),xdat);
%         normy = dat2(b2(n),ydat);
        if(strcmp(xnorm,'x normalized'))
            if(xdat == 39)
                xp = dat2(b2(n):e2(n),39) ./ dat2(b2(n):e2(n),40);
            else
                xp = dat2(b2(n):e2(n),xdat) ./ normx;
            end
        elseif(strcmp(xnorm,'x absolute'))
            xp = dat2(b2(n):e2(n),xdat);
        end
        if(strcmp(ynorm,'y normalized'))
            yp = dat2(b2(n):e2(n),ydat) ./ normy;
        elseif(strcmp(ynorm,'y absolute'))
            yp = dat2(b2(n):e2(n),ydat);
        end
        for i=1:length(bin)-1
            if(strcmp(stack,'linear stack'))
                b = find(xp>0 & log10(xp)>=bin(i) & log10(xp)<bin(i+1));
                sortedx2(i, n2(i)+1:n2(i)+length(b)) = xp(b);
                sortedy2(i, n2(i)+1:n2(i)+length(b)) = yp(b);
                n2(i) = n2(i) + length(b);
            elseif(strcmp(stack,'log stack'))
                b = find(xp>1e-9 & yp>1e-9 & log10(xp)>=bin(i) & log10(xp)<bin(i+1));
                if(xdat==13 | xdat==14)
                    sortedx2(i, n2(i)+1:n2(i)+length(b)) = xp(b);
                    sortedy2(i, n2(i)+1:n2(i)+length(b)) = yp(b);
                else
                    sortedx2(i, n2(i)+1:n2(i)+length(b)) = log10(xp(b));
                    sortedy2(i, n2(i)+1:n2(i)+length(b)) = log10(yp(b));
                end
                n2(i) = n2(i) + length(b);
            end
            clear b
        end
    end
    for i=1:length(bin)-1
        x2(i) = mean(sortedx2(i,1:n2(i)));
        y2(i) = mean(sortedy2(i,1:n2(i)));
        s2(i) = std(sortedy2(i,1:n2(i)));
        xmed2(i) = median(sortedx2(i,1:n2(i)));
        ymed2(i) = median(sortedy2(i,1:n2(i)));
%        [x2(i), xmed2(i), y2(i), ymed2(i)]
%        sortedy2(i,1:n2(i))
        
        sortedx2(i,1:n2(i)) = sort(sortedx2(i,1:n2(i)));
        sortedy2(i,1:n2(i)) = sort(sortedy2(i,1:n2(i)));
%        sortedy2(i,1:n2(i))
        if(n2(i)<=1)
            xl2(i) = xmed2(i);
            yl2(i) = ymed2(i);
            xu2(i) = xmed2(i);
            yu2(i) = ymed2(i);
        elseif(n2(i)==2)
            xl2(i) = (5*sortedx2(i,1) + sortedx2(i,2))/6;
            yl2(i) = (5*sortedy2(i,1) + sortedy2(i,2))/6;
            xu2(i) = (sortedx2(i,1) + 5*sortedx2(i,2))/6;
            yu2(i) = (sortedy2(i,1) + 5*sortedy2(i,2))/6;
        elseif(n2(i)==3)
            xl2(i) = (2*sortedx2(i,1) + sortedx2(i,2))/3;
            yl2(i) = (2*sortedy2(i,1) + sortedy2(i,2))/3;
            xu2(i) = (sortedx2(i,2) + 2*sortedx2(i,3))/3;
            yu2(i) = (sortedy2(i,2) + 2*sortedy2(i,3))/3;
        elseif(n2(i)==4)
            xl2(i) = (sortedx2(i,1) + sortedx2(i,2))/2;
            yl2(i) = (sortedy2(i,1) + sortedy2(i,2))/2;
            xu2(i) = (sortedx2(i,3) + sortedx2(i,4))/2;
            yu2(i) = (sortedy2(i,3) + sortedy2(i,4))/2;
        elseif(n2(i)==5)
            xl2(i) = sortedx2(i,2);
            yl2(i) = sortedy2(i,2);
            xu2(i) = sortedx2(i,4);
            yu2(i) = sortedy2(i,4);
        elseif(n2(i)==6)
            xl2(i) = sortedx2(i,2);
            yl2(i) = sortedy2(i,2);
            xu2(i) = sortedx2(i,5);
            yu2(i) = sortedy2(i,5);
        elseif(n2(i)==7)
            xl2(i) = sortedx2(i,2);
            yl2(i) = sortedy2(i,2);
            xu2(i) = sortedx2(i,6);
            yu2(i) = sortedy2(i,6);
        elseif(n2(i)==8)
            xl2(i) = (sortedx2(i,2) + sortedx2(i,3))/2;
            yl2(i) = (sortedy2(i,2) + sortedy2(i,3))/2;
            xu2(i) = (sortedx2(i,6) + sortedx2(i,7))/2;
            yu2(i) = (sortedy2(i,6) + sortedy2(i,7))/2;
        elseif(n2(i)==9)
            xl2(i) = (sortedx2(i,2) + sortedx2(i,3))/2;
            yl2(i) = (sortedy2(i,2) + sortedy2(i,3))/2;
            xu2(i) = (sortedx2(i,7) + sortedx2(i,8))/2;
            yu2(i) = (sortedy2(i,7) + sortedy2(i,8))/2;
        elseif(n2(i)==10)
            xl2(i) = sortedx2(i,3);
            yl2(i) = sortedy2(i,3);
            xu2(i) = sortedx2(i,8);
            yu2(i) = sortedy2(i,8);
        elseif(n2(i)==11)
            xl2(i) = sortedx2(i,3);
            yl2(i) = sortedy2(i,3);
            xu2(i) = sortedx2(i,9);
            yu2(i) = sortedy2(i,9);
        elseif(n2(i)==12)
            xl2(i) = sortedx2(i,3);
            yl2(i) = sortedy2(i,3);
            xu2(i) = sortedx2(i,10);
            yu2(i) = sortedy2(i,10);
        elseif(n2(i)==13)
            xl2(i) = sortedx2(i,3);
            yl2(i) = sortedy2(i,3);
            xu2(i) = sortedx2(i,11);
            yu2(i) = sortedy2(i,11);
        elseif(n2(i)==14)
            xl2(i) = (sortedx2(i,3)+sortedx2(i,4))/2;
            yl2(i) = (sortedy2(i,3)+sortedy2(i,4))/2;
            xu2(i) = (sortedx2(i,11)+sortedx2(i,12))/2;
            yu2(i) = (sortedy2(i,11)+sortedy2(i,12))/2;
        elseif(n2(i)==15)
            xl2(i) = (sortedx2(i,3)+sortedx2(i,4))/2;
            yl2(i) = (sortedy2(i,3)+sortedy2(i,4))/2;
            xu2(i) = (sortedx2(i,12)+sortedx2(i,13))/2;
            yu2(i) = (sortedy2(i,12)+sortedy2(i,13))/2;
        elseif(n2(i)==16)
            xl2(i) = (sortedx2(i,3)+sortedx2(i,4))/2;
            yl2(i) = (sortedy2(i,3)+sortedy2(i,4))/2;
            xu2(i) = (sortedx2(i,13)+sortedx2(i,14))/2;
            yu2(i) = (sortedy2(i,13)+sortedy2(i,14))/2;
        elseif(n2(i)==17)
            xl2(i) = sortedx2(i,4);
            yl2(i) = sortedy2(i,4);
            xu2(i) = sortedx2(i,14);
            yu2(i) = sortedy2(i,14);
        elseif(n2(i)==18)
            xl2(i) = sortedx2(i,4);
            yl2(i) = sortedy2(i,4);
            xu2(i) = sortedx2(i,15);
            yu2(i) = sortedy2(i,15);
        elseif(n2(i)==19)
            xl2(i) = sortedx2(i,4);
            yl2(i) = sortedy2(i,4);
            xu2(i) = sortedx2(i,16);
            yu2(i) = sortedy2(i,16);
        elseif(n2(i)==20)
            xl2(i) = (sortedx2(i,4)+sortedx2(i,5))/2;
            yl2(i) = (sortedy2(i,4)+sortedy2(i,5))/2;
            xu2(i) = (sortedx2(i,16)+sortedx2(i,17))/2;
            yu2(i) = (sortedy2(i,16)+sortedy2(i,17))/2;
        elseif(n2(i)==21)
            xl2(i) = (sortedx2(i,4)+sortedx2(i,5))/2;
            yl2(i) = (sortedy2(i,4)+sortedy2(i,5))/2;
            xu2(i) = (sortedx2(i,17)+sortedx2(i,18))/2;
            yu2(i) = (sortedy2(i,17)+sortedy2(i,18))/2;
        elseif(n2(i)==22)
            xl2(i) = (sortedx2(i,4)+sortedx2(i,5))/2;
            yl2(i) = (sortedy2(i,4)+sortedy2(i,5))/2;
            xu2(i) = (sortedx2(i,18)+sortedx2(i,19))/2;
            yu2(i) = (sortedy2(i,18)+sortedy2(i,19))/2;
        elseif(n2(i)==23)
            xl2(i) = sortedx2(i,5);
            yl2(i) = sortedy2(i,5);
            xu2(i) = sortedx2(i,19);
            yu2(i) = sortedy2(i,19);
        else
            nl = ceil(n2(i)/6);
            nu = floor(5*n2(i)/6);
            xl2(i) = sortedx2(i,nl);
            yl2(i) = sortedy2(i,nl);
            xu2(i) = sortedx2(i,nu);
            yu2(i) = sortedy2(i,nu);
        end
        %[n2(i) nl nu yl2(i) yu2(i)]
    end        
    
    e3(1:nmigrate3)-b3(1:nmigrate3)+1
    e2(1:nmigrate2)-b2(1:nmigrate2)+1
    mean([e3(1:nmigrate3)-b3(1:nmigrate3)+1])
    mean([e2(1:nmigrate2)-b2(1:nmigrate2)+1])
    
    for j=1:nmigrate3
        n = migrate3(j);
        normx = mean(dat3(b3(n):b3(n)+1,xdat));
        normy = mean(dat3(b3(n):b3(n)+1,ydat));
%         normx = dat3(b3(n),xdat);
%         normy = dat3(b3(n),ydat);
        if(strcmp(xnorm,'x normalized'))
            if(xdat == 39)
                xp = dat3(b3(n):e3(n),39) ./ dat3(b3(n):e3(n),40);
            else
                xp = dat3(b3(n):e3(n),xdat) ./ normx;
            end
        elseif(strcmp(xnorm,'x absolute'))
            xp = dat3(b3(n):e3(n),xdat);
        end
        if(strcmp(ynorm,'y normalized'))
            yp = dat3(b3(n):e3(n),ydat) ./ normy;
        elseif(strcmp(ynorm,'y absolute'))
            yp = dat3(b3(n):e3(n),ydat);
        end
        for i=1:length(bin)-1
            if(strcmp(stack,'linear stack'))
                b = find(xp>0 & log10(xp)>=bin(i) & log10(xp)<bin(i+1));
                sortedx3(i, n3(i)+1:n3(i)+length(b)) = xp(b);
                sortedy3(i, n3(i)+1:n3(i)+length(b)) = yp(b);
                n3(i) = n3(i) + length(b);
            elseif(strcmp(stack,'log stack'))
                b = find(xp>1e-9 & yp>1e-9 & log10(xp)>=bin(i) & log10(xp)<bin(i+1));
                if(xdat==13 | xdat==14)
                    sortedx3(i, n3(i)+1:n3(i)+length(b)) = xp(b);
                    sortedy3(i, n3(i)+1:n3(i)+length(b)) = yp(b);
                else
                    sortedx3(i, n3(i)+1:n3(i)+length(b)) = log10(xp(b));
                    sortedy3(i, n3(i)+1:n3(i)+length(b)) = log10(yp(b));
                end
                n3(i) = n3(i) + length(b);
            end
            clear b
        end
    end
    for i=1:length(bin)-1
        x3(i) = mean(sortedx3(i,1:n3(i)));
        y3(i) = mean(sortedy3(i,1:n3(i)));
        s3(i) = std(sortedy3(i,1:n3(i)));
        xmed3(i) = median(sortedx3(i,1:n3(i)));
        ymed3(i) = median(sortedy3(i,1:n3(i)));
%        [x3(i), xmed3(i), y3(i), ymed3(i)]
%        sortedy3(i,1:n3(i))        

        sortedx3(i,1:n3(i)) = sort(sortedx3(i,1:n3(i)));
        sortedy3(i,1:n3(i)) = sort(sortedy3(i,1:n3(i)));
        if(n3(i)<=1)
            xl3(i) = xmed3(i);
            yl3(i) = ymed3(i);
            xu3(i) = xmed3(i);
            yu3(i) = ymed3(i);
        elseif(n3(i)==2)
            xl3(i) = (5*sortedx3(i,1) + sortedx3(i,2))/6;
            yl3(i) = (5*sortedy3(i,1) + sortedy3(i,2))/6;
            xu3(i) = (sortedx3(i,1) + 5*sortedx3(i,2))/6;
            yu3(i) = (sortedy3(i,1) + 5*sortedy3(i,2))/6;
        elseif(n3(i)==3)
            xl3(i) = (2*sortedx3(i,1) + sortedx3(i,2))/3;
            yl3(i) = (2*sortedy3(i,1) + sortedy3(i,2))/3;
            xu3(i) = (sortedx3(i,2) + 2*sortedx3(i,3))/3;
            yu3(i) = (sortedy3(i,2) + 2*sortedy3(i,3))/3;
        elseif(n3(i)==4)
            xl3(i) = (sortedx3(i,1) + sortedx3(i,2))/2;
            yl3(i) = (sortedy3(i,1) + sortedy3(i,2))/2;
            xu3(i) = (sortedx3(i,3) + sortedx3(i,4))/2;
            yu3(i) = (sortedy3(i,3) + sortedy3(i,4))/2;
        elseif(n3(i)==5)
            xl3(i) = sortedx3(i,2);
            yl3(i) = sortedy3(i,2);
            xu3(i) = sortedx3(i,4);
            yu3(i) = sortedy3(i,4);
        elseif(n3(i)==6)
            xl3(i) = sortedx3(i,2);
            yl3(i) = sortedy3(i,2);
            xu3(i) = sortedx3(i,5);
            yu3(i) = sortedy3(i,5);
        elseif(n3(i)==7)
            xl3(i) = sortedx3(i,2);
            yl3(i) = sortedy3(i,2);
            xu3(i) = sortedx3(i,6);
            yu3(i) = sortedy3(i,6);
        elseif(n3(i)==8)
            xl3(i) = (sortedx3(i,2) + sortedx3(i,3))/2;
            yl3(i) = (sortedy3(i,2) + sortedy3(i,3))/2;
            xu3(i) = (sortedx3(i,6) + sortedx3(i,7))/2;
            yu3(i) = (sortedy3(i,6) + sortedy3(i,7))/2;
        elseif(n3(i)==9)
            xl3(i) = (sortedx3(i,2) + sortedx3(i,3))/2;
            yl3(i) = (sortedy3(i,2) + sortedy3(i,3))/2;
            xu3(i) = (sortedx3(i,7) + sortedx3(i,8))/2;
            yu3(i) = (sortedy3(i,7) + sortedy3(i,8))/2;
        elseif(n3(i)==10)
            xl3(i) = sortedx3(i,3);
            yl3(i) = sortedy3(i,3);
            xu3(i) = sortedx3(i,8);
            yu3(i) = sortedy3(i,8);
        elseif(n3(i)==11)
            xl3(i) = sortedx3(i,3);
            yl3(i) = sortedy3(i,3);
            xu3(i) = sortedx3(i,9);
            yu3(i) = sortedy3(i,9);
        elseif(n3(i)==12)
            xl3(i) = sortedx3(i,3);
            yl3(i) = sortedy3(i,3);
            xu3(i) = sortedx3(i,10);
            yu3(i) = sortedy3(i,10);
        elseif(n3(i)==13)
            xl3(i) = sortedx3(i,3);
            yl3(i) = sortedy3(i,3);
            xu3(i) = sortedx3(i,11);
            yu3(i) = sortedy3(i,11);
        elseif(n3(i)==14)
            xl3(i) = (sortedx3(i,3)+sortedx3(i,4))/2;
            yl3(i) = (sortedy3(i,3)+sortedy3(i,4))/2;
            xu3(i) = (sortedx3(i,11)+sortedx3(i,12))/2;
            yu3(i) = (sortedy3(i,11)+sortedy3(i,12))/2;
        elseif(n3(i)==15)
            xl3(i) = (sortedx3(i,3)+sortedx3(i,4))/2;
            yl3(i) = (sortedy3(i,3)+sortedy3(i,4))/2;
            xu3(i) = (sortedx3(i,12)+sortedx3(i,13))/2;
            yu3(i) = (sortedy3(i,12)+sortedy3(i,13))/2;
        elseif(n3(i)==16)
            xl3(i) = (sortedx3(i,3)+sortedx3(i,4))/2;
            yl3(i) = (sortedy3(i,3)+sortedy3(i,4))/2;
            xu3(i) = (sortedx3(i,13)+sortedx3(i,14))/2;
            yu3(i) = (sortedy3(i,13)+sortedy3(i,14))/2;
        elseif(n3(i)==17)
            xl3(i) = sortedx3(i,4);
            yl3(i) = sortedy3(i,4);
            xu3(i) = sortedx3(i,14);
            yu3(i) = sortedy3(i,14);
        elseif(n3(i)==18)
            xl3(i) = sortedx3(i,4);
            yl3(i) = sortedy3(i,4);
            xu3(i) = sortedx3(i,15);
            yu3(i) = sortedy3(i,15);
        elseif(n3(i)==19)
            xl3(i) = sortedx3(i,4);
            yl3(i) = sortedy3(i,4);
            xu3(i) = sortedx3(i,16);
            yu3(i) = sortedy3(i,16);
        elseif(n3(i)==20)
            xl3(i) = (sortedx3(i,4)+sortedx3(i,5))/2;
            yl3(i) = (sortedy3(i,4)+sortedy3(i,5))/2;
            xu3(i) = (sortedx3(i,16)+sortedx3(i,17))/2;
            yu3(i) = (sortedy3(i,16)+sortedy3(i,17))/2;
        elseif(n3(i)==21)
            xl3(i) = (sortedx3(i,4)+sortedx3(i,5))/2;
            yl3(i) = (sortedy3(i,4)+sortedy3(i,5))/2;
            xu3(i) = (sortedx3(i,17)+sortedx3(i,18))/2;
            yu3(i) = (sortedy3(i,17)+sortedy3(i,18))/2;
        elseif(n3(i)==22)
            xl3(i) = (sortedx3(i,4)+sortedx3(i,5))/2;
            yl3(i) = (sortedy3(i,4)+sortedy3(i,5))/2;
            xu3(i) = (sortedx3(i,18)+sortedx3(i,19))/2;
            yu3(i) = (sortedy3(i,18)+sortedy3(i,19))/2;
        elseif(n3(i)==23)
            xl3(i) = sortedx3(i,5);
            yl3(i) = sortedy3(i,5);
            xu3(i) = sortedx3(i,19);
            yu3(i) = sortedy3(i,19);
        else
            nl = ceil(n3(i)/6);
            nu = floor(5*n3(i)/6);
            xl3(i) = sortedx3(i,nl);
            yl3(i) = sortedy3(i,nl);
            xu3(i) = sortedx3(i,nu);
            yu3(i) = sortedy3(i,nu);
        end
    end
    
    for i=2:length(bin)-2
        if(n2(i) == 0)
            x2(i) = 0.5.*( x2(i+1) + x2(i-1) );
            y2(i) = 0.5.*( y2(i+1) + y2(i-1) );
            xmed2(i) = 0.5.*( xmed2(i+1) + xmed2(i-1) );
            ymed2(i) = 0.5.*( ymed2(i+1) + ymed2(i-1) );
            xu2(i) = xmed2(i);
            yu2(i) = ymed2(i);
            xl2(i) = xmed2(i);
            yl2(i) = ymed2(i);
        end
        if(n3(i) == 0)
            x3(i) = 0.5.*( x3(i+1) + x3(i-1) );
            y3(i) = 0.5.*( y3(i+1) + y3(i-1) );
            xmed3(i) = 0.5.*( xmed3(i+1) + xmed3(i-1) );
            ymed3(i) = 0.5.*( ymed3(i+1) + ymed3(i-1) );
            xu3(i) = xmed3(i);
            yu3(i) = ymed3(i);
            xl3(i) = xmed3(i);
            yl3(i) = ymed3(i);
        end
    end
    
    if(strcmp(stack,'linear stack'))
        if(strcmp(xaxis,'x linear scale') | strcmp(xaxis,'x log scale') | xdat==13 | xdat==14)
            if(strcmp(mean_type,'mean'))
                x = x2;
                xx = x3;
            elseif(strcmp(mean_type,'median'))
                x = xmed2;
                xx = xmed3;
            end
        else
            if(strcmp(mean_type,'mean'))
                x = log10(x2);
                xx = log10(x3);
            elseif(strcmp(mean_type,'median'))
                x = log10(xmed2);
                xx = log10(xmed3);
            end
        end
        if(strcmp(yaxis,'y linear scale') | strcmp(yaxis,'y log scale') | ydat==13 | ydat==14)
            if(strcmp(mean_type,'mean'))
                y = y2;
                l = s2;
                u = s2;
                yy = y3;
                ll = s3;
                uu = s3;
            elseif(strcmp(mean_type,'median'))
                y = ymed2;
                l = ymed2 - yl2;
                u = yu2 - ymed2;
                yy = ymed3;
                ll = ymed3 - yl3;
                uu = yu3 - ymed3;
            end
        else
            if(strcmp(mean_type,'mean'))
                y = log10(y2);
                l = log10(y2./(y2 - s2));
                u = log10((y2 + s2)./y2);
                yy = log10(y3);
                ll = log10(y3./(y3 - s3));
                uu = log10((y3 + s3)./y3);
            elseif(strcmp(mean_type,'median'))
                y = log10(ymed2);
                l = log10(ymed2./(yl2));
                u = log10((yu2)./ymed2);
                yy = log10(ymed3);
                ll = log10(ymed3./(yl3));
                uu = log10((yu3)./ymed3);
            end
        end
        
        if(strcmp(plot_type,'no error'))
            l(2) = plot(x, y,'linewidth',3,'color','r','marker','square','markersize',8,...
                'MarkerFaceColor','r','MarkerEdgeColor','r','DisplayName','No RP');
            l(1) = plot(xx, yy,'linewidth',3,'color','b','marker','square','markersize',8,...
                'MarkerFaceColor','b','MarkerEdgeColor','b','DisplayName','RP');
        elseif(strcmp(plot_type,'error'))
            l(2) = errorbar(x, y, l, u, 'linewidth',3,'color','r','marker','square','markersize',8,...
                'MarkerFaceColor','r','MarkerEdgeColor','r','DisplayName','No RP');
            l(1) = errorbar(xx, yy, ll, uu, 'linewidth',3,'color','b','marker','square','markersize',8,...
                'MarkerFaceColor','b','MarkerEdgeColor','b','DisplayName','RP');
        end
        
    elseif(strcmp(stack,'log stack'))
        if(strcmp(xaxis,'x log values') | xdat==13 | xdat==14)
            if(strcmp(mean_type,'mean'))
                x = x2;
                xx = x3;
            elseif(strcmp(mean_type,'median'))
                x = xmed2;
                xx = xmed3;
            end
        else
            if(strcmp(mean_type,'mean'))
                x = 10.^x2;
                xx = 10.^x3;
            elseif(strcmp(mean_type,'median'))
                x = 10.^xmed2;
                xx = 10.^xmed3;
            end
        end
        if(strcmp(yaxis,'y log values') | ydat==13 | ydat==14)
            if(strcmp(mean_type,'mean'))
                y = y2;
                l = s2;
                u = s2;
                yy = y3;
                ll = s3;
                uu = s3;
            elseif(strcmp(mean_type,'median'))
                y = ymed2;
                l = ymed2 - yl2;
                u = yu2 - ymed2;
                yy = ymed3;
                ll = ymed3 - yl3;
                uu = yu3 - ymed3;
            end
        else
            if(strcmp(mean_type,'mean'))
                y = 10.^y2;
                l = 10.^y2 - 10.^(y2 - s2);
                u = 10.^(y2 + s2) - 10.^y2;
                yy = 10.^y3;
                ll = 10.^y3 - 10.^(y3 - s3);
                uu = 10.^(y3 + s3) - 10.^y3;
            elseif(strcmp(mean_type,'median'))
                y = 10.^ymed2;
                l = 10.^ymed2 - 10.^(yl2);
                u = 10.^(yu2) - 10.^ymed2;
                yy = 10.^ymed3;
                ll = 10.^ymed3 - 10.^(yl3);
                uu = 10.^(yu3) - 10.^ymed3;
            end
        end
        
        if(strcmp(plot_type,'no error'))
            l(2) = plot(x,y,'linewidth',3,'color','r','marker','square','markersize',8,...
                'MarkerFaceColor','r','MarkerEdgeColor','r','DisplayName','No RP');
            l(1) = plot(xx,yy,'linewidth',3,'color','b','marker','square','markersize',8,...
                'MarkerFaceColor','b','MarkerEdgeColor','b','DisplayName','RP');
        elseif(strcmp(plot_type,'error'))
            l(2) = errorbar(x, y, l, u, 'linewidth',3,'color','r','marker','square','markersize',8,...
                'MarkerFaceColor','r','MarkerEdgeColor','r','DisplayName','No RP');
            l(1) = errorbar(xx, yy, ll, uu, 'linewidth',3,'color','b','marker','square','markersize',8,...
                'MarkerFaceColor','b','MarkerEdgeColor','b','DisplayName','RP');
        end
    end
    
    legend1=legend(gca,l(1:2));
    if(legend_location == 'se')
        set(legend1,'edgecolor','w',...
            'Position',[0.7006 0.142 0.1773 0.1056],...
            'fontname','Times New Roman','fontsize',16);
    elseif(legend_location == 'sw')
        set(legend1,'edgecolor','w',...
            'Position',[0.224 0.142 0.1654 0.1056],...
            'fontname','Times New Roman','fontsize',16);
    elseif(legend_location == 'nw')
    set(legend1,'edgecolor','w',...
        'Position',[0.2204 0.8096 0.1773 0.1056],...
        'fontname','Times New Roman','fontsize',16);
    elseif(legend_location == 'ne')
    set(legend1,'edgecolor','w',...
        'Position',[0.7095 0.8122 0.1654 0.1056],...
        'fontname','Times New Roman','fontsize',16);
    end
    [length(b2),nmigrate2(1),length(b3),nmigrate3(1)]
    b2(migrate2);
    e2(migrate2);
    b3(migrate3);
    e3(migrate3);
    
    if(max(xdat,ydat)==39 & min(xdat,ydat)==12)
        plot([-50:50:3000],[-50:50:3000],'linewidth',3,'linestyle','--','color','k','marker','none')
    end
end