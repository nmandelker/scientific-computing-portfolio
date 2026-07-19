function clump_property_histograms(is3, norm_is3, nis3, ... 
    norm, ind3, clump_dens_thresh, norm_lifetime_thresh, mass_thresh, max_mass_thresh, ...
    type)
% type = 0/1/2/3 -> stellar mass / baryonic mass / SFR / maximal baryonic mass function (on X axis).
% norm = 0 / 1 -> X or X/Xd on X axis.
% Compares Gen2 and Gen3 mass functions with cuts in clump density and
% lifetime
% ONLY in situ clumps are used
% If common==1 then ONLY common snapshots to Gen2 and Gen3 are used

filename = '';

if(type==0)
    m = 5;
elseif(type==1)
    m = 6;
elseif(type==2)
    m = 15;
elseif(type==3)
    m = 49;
elseif(type==5)
    m = 48;
elseif(type==6)
    m = 11;
elseif(type==7)
    %m = 52;
    m = 48;
elseif(type==8)
    m = 41;
elseif(type==9)
    m = 61;
elseif(type==10)
    m = 24;
elseif(type==11)
    m = 23;
elseif(type==12)
    m = 10;
elseif(type==13)
    m = 46;
elseif(type==14)
    m = 62;
elseif(type==15)
    m = 63;
elseif(type==16)
    m = 64;
elseif(type==17)
    m = 78;
else
    disp('invalid type')
    type
    return    
end

if(norm == 1)
    x3 = norm_is3;
    if(type==0)
        titx1 = '${\rm log(}\:\:M_{\rm c,*}\:/\:M_{\rm d,*}\:{\rm )}$';
        titx2 = '${\rm log(}M_{\rm c,*}\:/\:M_{\rm d,*}{\rm )}$';
        bin = [-5:0.1:1];
        xl = -4;
        xu = 0.5;
        xd = 0.5;
    elseif(type==1)
        titx1 = '${\rm log(}\:\:M_{\rm c}\:/\:M_{\rm d}\:{\rm )}$';
        titx2 = '${\rm log(}M_{\rm c}\:/\:M_{\rm d}{\rm )}$';
        bin = [-5:0.1:1];
        xl = -4;
        xu = 0.5;
        xd = 0.5;
    elseif(type==2)
        titx1 = '${\rm log(}\:\:SFR_{\rm c}\:/\:SFR_{\rm d}\:{\rm )}$';
        titx2 = '${\rm log(}SFR_{\rm c}\:/\:SFR_{\rm d}{\rm )}$';
        bin = [-5:0.1:1];
        xl = -4;
        xu = 0.5;
        xd = 0.5;
    elseif(type==3)
        titx1 = '${\rm log(}\:\:M_{\rm c,max}\:/\:M_{\rm d}\:{\rm )}$';
        titx2 = '${\rm log(}M_{\rm c,max}\:/\:M_{\rm d}{\rm )}$';
        bin = [-5:0.1:1];
        xl = -4;
        xu = 0.5;
        xd = 0.5;
    elseif(type==5)
        titx1 = '${\rm log(}\:\:n_{\rm c}\:/\:n_{\rm d}\:{\rm )}$';
        titx2 = '${\rm log(}n_{\rm c}\:/\:n_{\rm d}{\rm )}$';
        bin = [-2:0.1:5];
        xl = -1;
        xu = 4.5;
        xd = 0.5;
    elseif(type==6)
        titx1 = '${\rm log(}\:\:\Sigma_{\rm c}\:/\:\Sigma_{\rm d}\:{\rm )}$';
        titx2 = '${\rm log(}\Sigma_{\rm c}\:/\:\Sigma_{\rm d}{\rm )}$';
        bin = [-4:0.1:3];
        xl = -3;
        xu = 2.5;
        xd = 0.5;
    elseif(type==7)
        titx1 = '${\rm log(}\:\:t_{\rm c,\,max}\:/\:t_{\rm ff}\:{\rm )}$';
        titx2 = '${\rm log(}t_{\rm c,\,max}\:/\:t_{\rm ff}{\rm )}$';
        bin = [-3:0.12:4];
        xl = -2;
        xu = 3;
        xd = 0.5;
    elseif(type==8)
        titx1 = '${\rm log(}\:\:t\:/\:t_{\rm ff}\:{\rm )}$';
        titx2 = '${\rm log(}t\:/\:t_{\rm ff}{\rm )}$';
        bin = [-3:0.15:4];
        xl = -2;
        xu = 3;
        xd = 0.5;
    elseif(type==9)
        titx1 = '${\rm log(}\:\:t_{\rm max}\:/\:t_{\rm d}\:{\rm )}$';
        titx2 = '${\rm log(}t_{\rm max}\:/\:t_{\rm d}{\rm )}$';
        bin = [-3:0.15:4];
        xl = -2;
        xu = 3;
        xd = 0.5;
    elseif(type==10)
        titx1 = '$S_{\rm c}$';
        titx2 = '$S_{\rm c}$';
        bin = [-0.1:0.025:1.1];
        xl = 0;
        xu = 1;
        xd = 0.1;
        x3 = 10.^x3;
    elseif(type==11)
        titx1 = '${\rm log(}\:\:{\bar {\delta}}_{\rm c}\:{\rm )}$';
        titx2 = '${\rm log(}{\bar {\delta}}_{\rm c}{\rm )}$';
        bin = [0.9:0.075:4];
        xl = 1;
        xu = 3;
        xd = 0.2;
    elseif(type==12)
        titx1 = '${\rm log(}\:\:\Sigma_{\rm *,c}\:/\:\Sigma_{\rm *,d}\:{\rm )}$';
        titx2 = '${\rm log(}\Sigma_{\rm *,c}\:/\:\Sigma_{\rm *,d}{\rm )}$';
        bin = [-4:0.1:3];
        xl = -3;
        xu = 2.5;
        xd = 0.5;
    elseif(type==13)
        titx1 = '${\rm log(}\:\:n_{\rm gas,\:c}\:/\:n_{\rm gas,\:d}\:{\rm )}$';
        titx2 = '${\rm log(}n_{\rm gas,\:c}\:/\:n_{\rm gas,\:d}{\rm )}$';
        bin = [-2:0.1:5];
        xl = -1;
        xu = 4.5;
        xd = 0.5;
    elseif(type==14)
        titx1 = '${\rm log(}\:\:V_{\rm circ,\:c}\:\:\:{\rm [km\:s^{-1}\:] )}$';
        titx2 = '${\rm log(}V_{\rm circ,\:c}\:{\rm [km\:s^{-1}] )}$';
        bin = [-2:0.1:6];
        xl = 0;
        xu = 3;
        xd = 0.2;
    elseif(type==15)
        titx1 = '${\rm log(}\epsilon_{\rm ff,\:c})$';
        titx2 = '${\rm log(}\epsilon_{\rm ff,\:c})$';
        bin = [-5:0.1:1];
        xl = -4;
        xu = 0.4;
        xd = 0.4;
    elseif(type==16)
        titx1 = '$\eta_{\rm c}$';
        titx2 = '$\eta_{\rm c}$';
        bin = [-1:0.2:12];
        xl = 0;
        xu = 10;
        xd = 1;
        x3 = 10.^x3;
    elseif(type==17)
        titx1 = '$\alpha_{\rm v}$';
        titx2 = '$\alpha_{\rm v}$';
        bin = [-10:0.2:10];
        xl = -0.5;
        xu = 3;
        xd = 0.5;
    end
elseif(norm == 0)
    x3 = is3;
    if(type==0)
        titx1 = '${\rm log(}\:\:M_{\rm c,*}\:{\rm [M_{\odot}\:] )}$';
        titx2 = '${\rm log(}M_{\rm c,*}\:{\rm [M_{\odot}] )}$';
        bin = [3:0.1:11];
        xl = 4.5;
        xu = 9.5;
        xd = 0.5;
    elseif(type==1)
        titx1 = '${\rm log(}\:\:M_{\rm c}\:{\rm [M_{\odot}\:] )}$';
        titx2 = '${\rm log(}M_{\rm c}\:{\rm [M_{\odot}] )}$';
        bin = [3:0.1:11];
        xl = 6;
        xu = 9.5;
        xd = 0.5;
    elseif(type==2)
        titx1 = '${\rm log(}\:\:SFR_{\rm c}\:{\rm [M_{\odot}\:yr^{-1}\:] )}$';
        titx2 = '${\rm log(}SFR_{\rm c}\:{\rm [M_{\odot}\:yr^{-1}] )}$';
        bin = [-5:0.1:2];
        xl = -4;
        xu = 1;
        xd = 0.5;
    elseif(type==3)
        titx1 = '${\rm log(}\:\:M_{\rm c,max}\:\:{\rm [M_{\odot}\:] )}$';
        titx2 = '${\rm log(}M_{\rm c,max}\:\:{\rm [M_{\odot}] )}$';
        bin = [3:0.1:11];
        xl = 6;
        xu = 9.5;
        xd = 0.5;
    elseif(type==5)
        titx1 = '${\rm log(}\:\:n_{\rm c}\:{\rm [cm^{-3}] )}$';
        titx2 = '${\rm log(}n_{\rm c}\:{\rm [cm^{-3}] )}$';
        bin = [-3:0.1:4];
        xl = -2;
        xu = 3.5;
        xd = 0.5;
    elseif(type==6)
        titx1 = '${\rm log(}\:\:\Sigma_{\rm c}\:{\rm [M_{\odot}\:pc^{-2}\:] )}$';
        titx2 = '${\rm log(}\Sigma_{\rm c}\:{\rm [M_{\odot}\:pc^{-2}] )}$';
        bin = [-2:0.1:5];
        xl = -1;
        xu = 4;
        xd = 0.5;
    elseif(type==7)
        titx1 = '${\rm log(}\:\:t_{\rm c,\,max}\:{\rm [Myr\:] )}$';
        titx2 = '${\rm log(}t_{\rm c,\,max}\:{\rm [Myr] )}$';
        bin = [-1.5:0.05:5];
        xl = -1;
        xu = 4;
        xd = 0.5;
    elseif(type==8)
        titx1 = '${\rm log(}\:\:t\:{\rm [Myr\:] )}$';
        titx2 = '${\rm log(}t\:{\rm [Myr] )}$';
        bin = [-0.5:0.1:5];
        xl = 0;
        xu = 4;
        xd = 0.5;
    elseif(type==9)
        titx1 = '${\rm log(}\:\:t\:/\:t_{\rm d}\:{\rm )}$';
        titx2 = '${\rm log(}t\:/\:t_{\rm d}{\rm )}$';
        bin = [-3:0.15:4];
        xl = -2;
        xu = 3;
        xd = 0.5;
    elseif(type==10)
        titx1 = '$S_{\rm c}$';
        titx2 = '$S_{\rm c}$';
        bin = [-0.1:0.025:1.1];
        xl = 0;
        xu = 1;
        xd = 0.1;
        x3 = 10.^x3;
    elseif(type==11)
        titx1 = '${\rm log(}\:\:{\bar {\delta}}_{\rm c}\:{\rm )}$';
        titx2 = '${\rm log(}{\bar {\delta}}_{\rm c}{\rm )}$';
        bin = [0.9:0.075:4];
        xl = 1;
        xu = 3;
        xd = 0.2;
    elseif(type==12)
        titx1 = '${\rm log(}\:\:\Sigma_{\rm *,c}\:{\rm [M_{\odot}\:pc^{-2}\:] )}$';
        titx2 = '${\rm log(}\Sigma_{\rm *,c}\:{\rm [M_{\odot}\:pc^{-2}] )}$';
        bin = [-2:0.1:5];
        xl = -1;
        xu = 4;
        xd = 0.5;
    elseif(type==13)
        titx1 = '${\rm log(}\:\:n_{\rm gas,\:c}\:{\rm [cm^{-3}] )}$';
        titx2 = '${\rm log(}n_{\rm gas,\:c}\:{\rm [cm^{-3}] )}$';
        bin = [-3:0.1:4];
        xl = -2;
        xu = 3.5;
        xd = 0.5;
    elseif(type==14)
        titx1 = '${\rm log(}\:\:V_{\rm circ,\:c}\:\:\:{\rm [km\:s^{-1}\:] )}$';
        titx2 = '${\rm log(}V_{\rm circ,\:c}\:{\rm [km\:s^{-1}] )}$';
        bin = [-2:0.1:6];
        xl = 0.2;
        xu = 2.6;
        xd = 0.2;
    elseif(type==15)
        titx1 = '${\rm log(}\epsilon_{\rm ff,\:c})$';
        titx2 = '${\rm log(}\epsilon_{\rm ff,\:c})$';
        bin = [-5:0.1:1];
        xl = -4;
        xu = 0.4;
        xd = 0.4;
    elseif(type==16)
        titx1 = '$\eta_{\rm c}$';
        titx2 = '$\eta_{\rm c}$';
        bin = [-1:0.2:12];
        xl = 0;
        xu = 10;
        xd = 1;
        x3 = 10.^x3;
    elseif(type==17)
        titx1 = '${\rm log(}\alpha_{\rm v})$';
        titx2 = '${\rm log(}\alpha_{\rm v})$';
        bin = [-10:0.2:10];
        xl = -1;
        xu = 3;
        xd = 0.5;
    end
end
dbin = (bin(2)-bin(1))/2;
multi = 0;
% dens thresh
if(clump_dens_thresh(2) <= 100)
    multi = 1;
    max_dens = clump_dens_thresh(2);
    tit = strcat('$\:',num2str(clump_dens_thresh(2),'%3i'),'cm^{-3}>n_{\rm c}');
else
    max_dens = 1e5;
    if(clump_dens_thresh(1) >= 1)
        tit = '$\:n_{\rm c}';
    else
        tit = '$\:';
    end
end
if(clump_dens_thresh(1) >= 1)
    multi = 1;
    min_dens = clump_dens_thresh(1);
    tit = strcat(tit,'>',num2str(clump_dens_thresh(1),'%3i'),'cm^{-3}');
else
    min_dens = 0;
end
if(min_dens>0)
    if(~strcmp(filename,''))
        filename = strcat(filename,'_');
    end
    filename = strcat(filename,'nmin_',num2str(min_dens,'%4.1f'));
end
if(max_dens<1e5)
    if(~strcmp(filename,''))
        filename = strcat(filename,'_');
    end
    filename = strcat(filename,'nmax_',num2str(max_dens,'%4.1f'));
end

% mass thresh
if(max_mass_thresh(2)>10 & max_mass_thresh(1)<6)
    if(mass_thresh(2) <= 10)
        multi = 1;
        max_mass = mass_thresh(2);
        tit = strcat(tit,'\:',num2str(mass_thresh(2),'%3.1f'),'>{\rm log}(M_{\rm c})');
    else
        max_mass = 1e5;
        if(mass_thresh(1) >= 6)
            tit = strcat(tit,'\:{\rm log}(M_{\rm c})');
        else
            tit = strcat(tit,'\:');
        end
    end
    if(mass_thresh(1) >= 6)
        multi = 1;
        min_mass = mass_thresh(1);
        tit = strcat(tit,'>',num2str(mass_thresh(1),'%3.1f'));
    else
        min_mass = 0;
    end
    if(min_mass>0)
        if(~strcmp(filename,''))
            filename = strcat(filename,'_');
        end
        filename = strcat(filename,'Mmin_',num2str(min_mass,'%3.1f'));
    end
    if(max_mass<1e5)
        if(~strcmp(filename,''))
            filename = strcat(filename,'_');
        end
        filename = strcat(filename,'Mmax_',num2str(max_mass,'%3.1f'));
    end
else
    min_mass = -1e9;
    max_mass = 1e9;
end
tit
% max mass thresh
if(max_mass_thresh(2) <= 10)
    multi = 1;
    max_max_mass = max_mass_thresh(2);
    tit = strcat(tit,'\:',num2str(max_mass_thresh(2),'%3.1f'),'>{\rm log}(M_{\rm c,\,max}\:)');
else
    max_max_mass = 1e5;
    if(max_mass_thresh(1) >= 6)
        tit = strcat(tit,'\:{\rm log}(M_{\rm c,\,max}\:)');
%     else
%         tit = '\:';
    end
end
if(max_mass_thresh(1) >= 6)
    multi = 1;
    min_max_mass = max_mass_thresh(1);
    tit = strcat(tit,'>',num2str(max_mass_thresh(1),'%3.1f'));
else
    min_max_mass = 0;
end
if(min_max_mass>0)
    if(~strcmp(filename,''))
        filename = strcat(filename,'_');
    end
    filename = strcat(filename,'MaxMmin_',num2str(min_max_mass,'%3.1f'));
end
if(max_max_mass<1e5)
    if(~strcmp(filename,''))
        filename = strcat(filename,'_');
    end
    filename = strcat(filename,'MaxMmax_',num2str(max_max_mass,'%3.1f'));
end
tit
% norm lifetime thresh (clump free fall time)
if(norm_lifetime_thresh(2) < 10000)
    multi = 1;
    if(~strcmp(tit,'\:'))
        tit = strcat(tit,',\:');
    end
    max_norm_lifetime = norm_lifetime_thresh(2);
    if(norm_lifetime_thresh(2) < 10)
        tit = strcat(tit,num2str(norm_lifetime_thresh(2),'%3.1f'),'>t_{\rm max}/t_{\rm ff}');
    elseif(norm_lifetime_thresh(2) < 100)
        tit = strcat(tit,num2str(norm_lifetime_thresh(2),'%4.1f'),'>t_{\rm max}/t_{\rm ff}');
    elseif(norm_lifetime_thresh(2) < 1000)
        tit = strcat(tit,num2str(norm_lifetime_thresh(2),'%5.1f'),'>t_{\rm max}/t_{\rm ff}');
    else
        tit = strcat(tit,num2str(norm_lifetime_thresh(2),'%6.1f'),'>t_{\rm max}/t_{\rm ff}');
    end
else
    max_norm_lifetime = 10000;
    if(norm_lifetime_thresh(1) > 0 & norm_lifetime_thresh(1)~=20)
        if(~strcmp(tit,'\:'))
            tit = strcat(tit,',\:');
        end
        tit = strcat(tit,'t_{\rm max}/t_{\rm ff}');
    end
end
if(norm_lifetime_thresh(1) > 0)
    multi = 1;
    min_norm_lifetime = norm_lifetime_thresh(1);
    if(norm_lifetime_thresh(1) == 20)
        tit = strcat(tit,'LLCs');
    elseif(norm_lifetime_thresh(1) >= 1000)
        tit = strcat(tit,'>',num2str(norm_lifetime_thresh(1),'%6.1f'));
    elseif(norm_lifetime_thresh(1) >= 100)
        tit = strcat(tit,'>',num2str(norm_lifetime_thresh(1),'%5.1f'));
    elseif(norm_lifetime_thresh(1) >= 10)
        tit = strcat(tit,'>',num2str(norm_lifetime_thresh(1),'%4.1f'));
    else
        tit = strcat(tit,'>',num2str(norm_lifetime_thresh(1),'%3.1f'));
    end
else
    min_norm_lifetime = -10;
end
if(min_norm_lifetime>0)
    if(~strcmp(filename,''))
        filename = strcat(filename,'_');
    end
    if(min_norm_lifetime>=1000)
        filename = strcat(filename,'NormTmin_',num2str(ceil(min_norm_lifetime),'%6.1f'));
    elseif(min_norm_lifetime>=100)
        filename = strcat(filename,'NormTmin_',num2str(ceil(min_norm_lifetime),'%5.1f'));
    elseif(min_norm_lifetime>=10)
        filename = strcat(filename,'NormTmin_',num2str(ceil(min_norm_lifetime),'%4.1f'));
    else
        filename = strcat(filename,'NormTmin_',num2str(ceil(min_norm_lifetime),'%3.1f'));
    end
end
if(max_norm_lifetime<10000)
    if(~strcmp(filename,''))
        filename = strcat(filename,'_');
    end
    if(max_norm_lifetime<10)
        filename = strcat(filename,'NormTmax_',num2str(ceil(max_norm_lifetime),'%3.1f'));
    elseif(max_norm_lifetime<100)
        filename = strcat(filename,'NormTmax_',num2str(ceil(max_norm_lifetime),'%4.1f'));
    elseif(max_norm_lifetime<1000)
        filename = strcat(filename,'NormTmax_',num2str(ceil(max_norm_lifetime),'%5.1f'));
    else
        filename = strcat(filename,'NormTmax_',num2str(ceil(max_norm_lifetime),'%6.1f'));
    end
end
tit
if(multi==0)
    filename = 'all';
end
tit = strcat(tit,'$')
tity = '${\bar {N}}_{\rm c}$';

yu = 10;
yl = 0.001;
ytick = [0.001 0.01 0.1 1 10];
ystr = {num2str(0.001,'%5.3f'),num2str(0.01,'%4.2f'),num2str(0.1,'%3.1f'),num2str(1,'%1.0f'),num2str(10,'%2.0f')};
yl = 0;
yu = 10;
ytick = yl:1:yu;

estr = {'','','','','',''};
if(type==7 | type==10 | type==11)
    xtick = xl:xd:xu;
else
    xtick = (xl+xd):xd:xu;
end

figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
set(figure1,'WindowStyle','docked');

axes1 = axes('Parent',figure1,'YTick',ytick,... 'YScale','log',... 'YtickLabel',ystr,...
    'YMinorTick','on',...
    'XTick',xtick,...
    'XMinorTick','on',...'
    'TickLength',[0.02 0.04],...
    'PlotBoxAspectRatio',[1 1 1],...
    'LineWidth',2,...
    'FontSize',16,...
    'FontName','Arial');

xlim(axes1,[xl xu]);
ylim(axes1,[yl yu]);
box(axes1,'on');
hold(axes1,'all');

% Create xlabel
xlab = xlabel(titx1,...
    'Interpreter','latex',...
    'FontSize',18,...
    'FontName','Times New Roman','units','normalized','position',[0.5 -0.06 0]);

% Create ylabel
ylab = ylabel(tity,...
    'Interpreter','latex','FontSize',18,...
    'FontName','Times New Roman',...
    'units','normalized','position',[-0.08 0.5 0]);

% Create title
title(tit,...
    'Interpreter','latex','FontSize',18,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 1.01 0]);%,'BackgroundColor',[.7 .9 .7]);

set(figure1,'renderer','painters')
% ti = get(gca,'TightInset');
% set(gca,'Position',[ti(1)+0.05 ti(2)+0.08 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
% set(gca,'units','centimeters')
% pos = get(gca,'Position');
% set(gcf, 'PaperUnits','centimeters');
% set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3)+1.6 pos(4)+ti(2)+ti(4)+3.2]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3)+5.4 pos(4)+ti(2)+ti(4)+3.2]);
    
ngal_3 = sum(length(unique(is3(nis3(1:end-1)+1:nis3(2:end),2))));
clump3 = 0.001.*ones(3,length(bin));
for i=1:length(bin)
%     if(type==5 | type==6 | type==17)
        clump3(1,i) = length(find(log10(x3(:,m))>=bin(i)-dbin & log10(x3(:,m))<bin(i)+dbin & ...
            log10(is3(:,6))>=min_mass & log10(is3(:,6))<max_mass & ...
            log10(is3(:,49))>=min_max_mass & log10(is3(:,49))<max_max_mass));
%     else
%         clump3(1,i) = length(find(log10(x3(:,m))>=bin(i)-dbin & log10(x3(:,m))<bin(i)+dbin));
%     end
    clump3(2,i) = length(find(log10(x3(:,m))>=bin(i)-dbin & log10(x3(:,m))<bin(i)+dbin & ...
        is3(:,48)>=min_dens & is3(:,48)<max_dens & ...
        norm_is3(:,52)>1e-6 & ...
        log10(is3(:,6))>=min_mass & log10(is3(:,6))<max_mass & ...
        log10(is3(:,49))>=min_max_mass & log10(is3(:,49))<max_max_mass));
    clump3(3,i) = length(find(log10(x3(:,m))>=bin(i)-dbin & log10(x3(:,m))<bin(i)+dbin & ...
        is3(:,48)>=min_dens & is3(:,48)<max_dens & ...
        norm_is3(:,52)>=min_norm_lifetime & norm_is3(:,52)<max_norm_lifetime & ...
        log10(is3(:,6))>=min_mass & log10(is3(:,6))<max_mass & ...
        log10(is3(:,49))>=min_max_mass & log10(is3(:,49))<max_max_mass));
end
[sum(clump3(1,:)), sum(clump3(2,:))]
clump3(1,:) = clump3(1,:) ./ ngal_3;
clump3(2,:) = max(0.001,clump3(2,:) ./ ngal_3);
clump3(3,:) = max(0.001,clump3(3,:) ./ ngal_3);

if(multi == 1)
    if(type~=7 & type~=8 & type ~=9 & min_norm_lifetime==20)
        l(1) = stairs(bin,clump3(1,:),'marker','none','linewidth',4,'color','r','linestyle',':','DisplayName','$RP,\:ZL+SL+LL$');%[.565,.22,.042]
        l(2) = stairs(bin,clump3(2,:),'marker','none','linewidth',3,'color','r','linestyle','--','DisplayName','$RP,\:SL+LL$');
        l(3) = stairs(bin,clump3(3,:),'marker','none','linewidth',3,'color','r','linestyle','-','DisplayName','$RP,\:LL$');
        leg_pos1 = [14.0 15.0 0.20 0.20];
        leg_pos2 = [6.5 15.0 0.20 0.20];
    else
%         l(1) = stairs(bin,clump3(1,:),'marker','none','linewidth',3,'color','r','linestyle','--','DisplayName','$RP,\:Non\:ZL$');
        l(1) = stairs(bin,clump3(3,:),'marker','none','linewidth',3,'color','r','linestyle','-','DisplayName',strcat('$RP,',tit));%[0.8706 0.4902 0]
    end
    if(type==7)
        %leg_pos = [12.4 14.2 0.20 0.20]; % good for max mass cut, right hand side
        leg_pos = [7 15.5 0.20 0.10]; % good for max mass cut, left hand side
    else
        leg_pos = [14.2 14.2 0.20 0.20];
    end
else
    l(1) = stairs(bin,clump3(1,:),'marker','none','linewidth',3,'color','r','linestyle','-','DisplayName','$RP$');
    leg_pos = [14.7 15.2 0.20 0.20];
end

if(type==10 & min_norm_lifetime==20)
    
    xlabel('Parent',axes1,titx2,'Interpreter','latex','FontSize',24,...
        'FontName','Times New Roman','units','normalized','position',[0.5 -0.06 0]);

    legend1=legend(gca,l(1:3));
    legend boxoff
    set(legend1,'edgecolor','w','Position',leg_pos1,...
        'fontname','Times New Roman','fontsize',14,...
        'Interpreter','Latex');
     ah2=axes('position',get(gca,'position'), 'visible','off');
     ti = get(gca,'TightInset');
    
    set(gca,'Position',[ti(1)+0.05 ti(2)+0.08 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3)+1.6 pos(4)+ti(2)+ti(4)+3.2]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3)+1.6 pos(4)+ti(2)+ti(4)+3.2]);
    
    legend2 = legend(ah2,l(4:6));
    legend boxoff
    set(legend2,'edgecolor','w','Position',leg_pos2,...
        'fontname','Times New Roman','fontsize',14,'interpreter','latex');
end

if(type==7)
    Nl = length(l);
    legend1=legend(gca,l(1:Nl));
    legend boxoff
    set(legend1,'edgecolor','w','Position',leg_pos,...
        'fontname','Times New Roman','fontsize',18,...
        'Interpreter','Latex');
    y=logspace(-6,4,11);
    x=log10(20).*y./y;
    plot(x,y,'linestyle','--','linewidth',2,'color','k','marker','none');
    % Create textbox
    annotation(gcf,'textbox',[0.40 0.72 0.09 0.0535],...
        'String','$SLCs$',...
        'Interpreter','Latex',...
        'FontSize',24,...
        'FontName','Times',...
        'FitBoxToText','off',...
        'LineStyle','none',...
        'BackgroundColor',[1 1 1]);
    annotation(gcf,'textbox',[0.67 0.72 0.09 0.0535],...
        'String','$LLCs$',...
        'Interpreter','Latex',...
        'FontSize',24,...
        'FontName','Times',...
        'FitBoxToText','off',...
        'LineStyle','none',...
        'BackgroundColor',[1 1 1]);
end

mkdir('./clump_property_histograms/');
mkdir(strcat('./clump_property_histograms/all_snapshots/'));
mkdir(strcat('./clump_property_histograms/all_snapshots/multiplot/'));
dirname = strcat('./clump_property_histograms/all_snapshots/multiplot/');
if(type==0)
    dirname = strcat(dirname,'Mstar/');
elseif(type==1)
    dirname = strcat(dirname,'Mbar/');
elseif(type==2)
    dirname = strcat(dirname,'SFR/');
elseif(type==3)
    dirname = strcat(dirname,'max_Mbar/');
elseif(type==5)
    dirname = strcat(dirname,'nbar/');
elseif(type==6)
    dirname = strcat(dirname,'Sigma/');
elseif(type==7)
    dirname = strcat(dirname,'lifetime/');
elseif(type==8)
    dirname = strcat(dirname,'time/');
elseif(type==9)
    dirname = strcat(dirname,'t_over_td/');
elseif(type==10)
    dirname = strcat(dirname,'shape/');
elseif(type==11)
    dirname = strcat(dirname,'residual/');
elseif(type==12)
    dirname = strcat(dirname,'Sigma_star/');
elseif(type==13)
    dirname = strcat(dirname,'ngas/');
elseif(type==14)
    dirname = strcat(dirname,'Vcirc/');
elseif(type==15)
    dirname = strcat(dirname,'eps_ff/');
elseif(type==16)
    dirname = strcat(dirname,'mass_loading/');
elseif(type==17)
    dirname = strcat(dirname,'alpha_vir/');
end
mkdir(dirname);
filename = strcat(dirname,filename);
if(norm==1)
    filename = strcat(filename,'_normalized');
end
filename1 = strcat(filename,'.jpg');
saveas(gcf,filename1);

set(ylab,'fontsize',24)
if(type==10 &  min_norm_lifetime==20)
    set(legend2,'fontsize',14)
    set(legend1,'fontsize',14)
end
filename1 = strcat(filename,'.eps');

xlabel('Parent',axes1,titx2,'Interpreter','latex','FontSize',24,...
    'FontName','Times New Roman','units','normalized','position',[0.5 -0.06 0]);
% if(type~=3 & type~=7)
%     legend hide
% end
set(gcf,'renderer','painters')
print(gcf, '-depsc', filename1);
close all
fclose all