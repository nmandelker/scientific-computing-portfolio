clear

KB = 1.38e-16;
mproton = 1.67e-24;
Myr = 60*60*24*365.25*1e6;
kpc = 3.086e21;
Zsolar = 0.02;

red = 2;    % for UVB
Rs = 3;     % kpc
Mb = 1;

ns    = 10.^([-4:0.01:0]);  % cm^{-3}
delta = 10.^([1:0.01:3]);
Ts    = 0.*ns;
mus   = 0.*ns;
Cs    = 0.*ns;
tsc   = 0.*ns;
Mtot  = 0.*delta;
alpha = 0.*delta;

Nd = length(delta);
Nn = length(ns);
Rs_crit_map    = zeros(Nn,Nd,2);

for i=1:Nd
    Mtot(i)  = ( sqrt(delta(i)) / (1+sqrt(delta(i))) ) * Mb;
    alpha(i) = 0.21 * ( 0.8*exp(-3*Mtot(i)^2) + 0.2 );
end

Zs = 0.03*Zsolar;  % absolute units
Zb = 0.10*Zsolar;  % absolute units
Zmix = sqrt(Zs*Zb);
tit1 = strcat('$M_{\rm b}=',num2str(Mb,'%3.1f'),',\,z=',num2str(red,'%1i'),...
    ',\,[Z_{\rm s}]=',num2str(log10(Zs/Zsolar),'%4.1f'),...
    ',\,[Z_{\rm b}]=',num2str(log10(Zb/Zsolar),'%4.1f'),'$');
for i=1:Nn
    [Ts(i),mus(i)] = find_cooling_equilibrium2(ns(i), red, Zs/Zsolar, 1);
    Cs(i) = sqrt( (5/3) * KB * Ts(i) / (mus(i) * mproton) ); % cm/s
    Cs(i) = Cs(i) * Myr / kpc; % kpc/Myr
    tsc(i)= 2 * Rs / Cs(i);    % Myr
end
for j=1:Nn
    for i=1:Nd
        nb = ns(j) / delta(i);
        Tb = Ts(j) * delta(i);
        Tmix = sqrt(Ts(j)*Tb);
        nmix = sqrt(ns(j)*nb);
        [nspec, mu, cool, heat, cool_tot, heat_tot] = cooling_equilibrium(Tmix, nmix, red, Zmix/Zsolar, 1, 1);
        tcool = ( 1.5*KB*Tmix / (nmix*(cool_tot-heat_tot)) ) / Myr;    % Myr        
        tdis = Rs / ( alpha(i) * Mb*Cs(j)*sqrt(delta(i)) );
        Rs_crit_map(j,i,1) = Rs .* tcool / tdis;
    end
end

Zs = 0.001*Zsolar;  % absolute units
Zb = 0.010*Zsolar;  % absolute units
Zmix = sqrt(Zs*Zb);
tit2 = strcat('$M_{\rm b}=',num2str(Mb,'%3.1f'),',\,z=',num2str(red,'%1i'),...
    ',\,[Z_{\rm s}]=',num2str(log10(Zs/Zsolar),'%4.1f'),...
    ',\,[Z_{\rm b}]=',num2str(log10(Zb/Zsolar),'%4.1f'),'$');
% tit2 = strcat('$[Z_{\rm s}]=',num2str(log10(Zs/Zsolar),'%4.1f'),...
%     ',\,[Z_{\rm b}]=',num2str(log10(Zb/Zsolar),'%4.1f'),'$');
for i=1:Nn
    [Ts(i),mus(i)] = find_cooling_equilibrium2(ns(i), red, Zs/Zsolar, 1);
    Cs(i) = sqrt( (5/3) * KB * Ts(i) / (mus(i) * mproton) ); % cm/s
    Cs(i) = Cs(i) * Myr / kpc; % kpc/Myr
    tsc(i)= 2 * Rs / Cs(i);    % Myr
end
for j=1:Nn
    for i=1:Nd
        nb = ns(j) / delta(i);
        Tb = Ts(j) * delta(i);
        Tmix = sqrt(Ts(j)*Tb);
        nmix = sqrt(ns(j)*nb);
        [nspec, mu, cool, heat, cool_tot, heat_tot] = cooling_equilibrium(Tmix, nmix, red, Zmix/Zsolar, 1, 1);
        tcool = ( 1.5*KB*Tmix / (nmix*(cool_tot-heat_tot)) ) / Myr;    % Myr        
        tdis = Rs / ( alpha(i) * Mb*Cs(j)*sqrt(delta(i)) );
        Rs_crit_map(j,i,2) = Rs .* tcool / tdis;
    end
end

cmap = [0 0 0;0.125 0.125 1;0.1071 0.1071 1;0.08929 0.08929 1;...
    0.07143 0.07143 1;0.05357 0.05357 1;0.03571 0.03571 1;...
    0.01786 0.01786 1;0 0 1;0 0.1111 1;0 0.2222 1;0 0.3333 1;0 0.4444 1;...
    0 0.5556 1;0 0.6667 1;0 0.7333 1;0 0.8 1;0 0.8667 1;0 0.9333 1;...
    0 1 1;0 1 0.975;0 1 0.95;0 1 0.925;0 1 0.9;0 1 0.7714;0 1 0.6429;...
    0 1 0.5143;0 1 0.3857;0 1 0.2571;0 1 0.1286;0 1 0;0.1111 1 0;...
    0.2222 1 0;0.3333 1 0;0.4444 1 0;0.5556 1 0;0.6667 1 0;0.7778 1 0;...
    0.8889 1 0;1 1 0;1 0.9821 0;1 0.9643 0;1 0.9464 0;1 0.9286 0;...
    1 0.898 0;1 0.8673 0;1 0.8367 0;1 0.8061 0;1 0.7755 0;1 0.7449 0;...
    1 0.7143 0;1 0.5952 0;1 0.4762 0;1 0.3571 0;1 0.2381 0;1 0.119 0;...
    1 0 0;0.9286 0 0;0.8571 0 0;0.7857 0 0;0.7143 0 0;0.6429 0 0;...
    0.5714 0 0;0.5 0 0];

xl = -4;
xu = 0;
dx = 0.5;
yl = 1;
yu = 3;
dy = 0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68],'Colormap',cmap);
set(figure1,'WindowStyle','docked');
axes1 = axes('Parent',figure1,'YTick',[yl:dy:yu],...
    'XTick',[xl:dx:xu-dx],...
    'PlotBoxAspectRatio',[1 1 1],...
    'FontSize',12,...
    'FontName','Times',...
    'Position',[0.12 0.13 0.79 0.80],...
    'CLim',[-3 2]);
box(axes1,'on');
axis square
xlim([xl xu]);
ylim([yl yu]);
xlabel('${\rm log}\:\:(n_{\rm H,s}\:{\rm [cm^{-3}]})$','Interpreter','latex','FontSize',16,'FontName','Times','units','normalized',...
    'position',[0.5,-0.08,0.0]);
ylabel('${\rm log}\:\:(\delta)$','Interpreter','latex','FontSize',16,'FontName','Times','units','normalized',...
    'position',[-0.1,0.5,0.0]);
title(tit1,'Interpreter','latex','FontSize',11,'FontName','Times New Roman','units','normalized',...    
    'position',[0.5 1.01 0]);
% bar=colorbar('peer',axes1,'FontSize',12,'FontName','Times',...
%     'YTick',[-3:0.5:2],'units','normalized','Position',[0.80,0.15,0.05,0.78]);
% set(get(bar,'title'),'String','${\rm log}(R_{\rm s,\,crit}\:{\rm [kpc]})$',...
%     'Interpreter','latex','Fontsize',16,...
%     'Rotation',270,'units','normalized','position',[2.7,0.5,9.16],'FontName','Times');
hold('all');
surf(log10(ns),log10(delta),log10(Rs_crit_map(:,:,1)'),'Parent',axes1,'LineStyle','None');
set(gcf,'renderer','zbuffer')
plot3(-3,2,6,'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',4)
plot3(-2,2,6,'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',4)
plot3(-1,2,6,'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',4)
plot3(-2,log10(30),6,'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',4)
M = contour(log10(ns),log10(delta),log10(Rs_crit_map(:,:,1)'),[log10(3) log10(3)],'linestyle','-','color','w','linewidth',2);
plot3(M(1,2:end),M(2,2:end),10.*ones(size(M(1,2:end))),'marker','none','linestyle','-','color','w','linewidth',2)

set(gcf,'visible','off');
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 20 9.5]);
set(gca, 'Position',[0.0420 0.1500 0.4740 0.7800]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axes2 = axes('Parent',figure1,'YTick',[ ],...
    'XTick',[xl:dx:xu],...
    'PlotBoxAspectRatio',[1 1 1],...
    'FontSize',12,...
    'FontName','Times',...
    'Position',[0.12 0.13 0.79 0.80],...
    'CLim',[-3 2]);
box(axes2,'on');
axis square
xlim([xl xu]);
ylim([yl yu]);
xlabel('${\rm log}\:\:(n_{\rm H,s}\:{\rm [cm^{-3}]})$','Interpreter','latex','FontSize',16,'FontName','Times','units','normalized',...
    'position',[0.5,-0.08,0.0]);
title(tit2,'Interpreter','latex','FontSize',11,'FontName','Times New Roman','units','normalized',...    
    'position',[0.5 1.01 0]);
bar=colorbar('peer',axes2,'FontSize',12,'FontName','Times',...
    'YTick',[-3:0.5:2],'units','normalized','Position',[0.85,0.15,0.03,0.78]);
set(get(bar,'title'),'String','${\rm log}(R_{\rm s,\,crit}\:{\rm [kpc]})$',...
    'Interpreter','latex','Fontsize',16,...
    'Rotation',270,'units','normalized','position',[2.7,0.5,9.16],'FontName','Times');
hold('all');
surf(log10(ns),log10(delta),log10(Rs_crit_map(:,:,2)'),'Parent',axes2,'LineStyle','None');
set(gcf,'renderer','zbuffer')
plot3(-3,2,6,'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',4)
plot3(-2,2,6,'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',4)
plot3(-1,2,6,'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',4)
plot3(-2,log10(30),6,'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',4)
M = contour(log10(ns),log10(delta),log10(Rs_crit_map(:,:,2)'),[log10(3) log10(3)],'linestyle','-','color','w','linewidth',2);
plot3(M(1,2:end),M(2,2:end),10.*ones(size(M(1,2:end))),'marker','none','linestyle','-','color','w','linewidth',2)

set(gcf,'visible','off');
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 20 9.5]);
set(gca, 'Position',[0.4200 0.1500 0.4740 0.7800]);

filename = strcat('./mixing_region_cooling/Rs_crit_panels');
print(gcf,'-djpeg',strcat(filename,'.jpg'));
print(gcf,'-depsc',strcat(filename,'.eps'));

close all
fclose all;