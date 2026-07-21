function solve_halo_decel_cooling_panels_NFW(c, beta, eta, Theta_h, M12, red, Zs0, Zb0)

z3  = (1+red)/3;

delta_vec = [30, 30, 30, 100, 100, 100, 300, 300, 300];% .* M12^(2/3) .* z3;
Rs_Rv_vec = [0.09, 0.29, 0.50, 0.05, 0.16, 0.28, 0.03, 0.09, 0.16];% .* M12^(-1/3) .* z3^(-1);
ns_vec    = [0.45, 0.45, 0.15, 1.50, 1.50, 0.50, 4.50, 4.50, 1.50];% .* M12^(2/3) .* z3^4;

Lnorm_vec = zeros(size(delta_vec));
B  = zeros(size(delta_vec));
n0 = zeros(size(delta_vec));
Ts = zeros(size(delta_vec));
Rs = zeros(size(delta_vec));
% Mb_vec = zeros(size(delta_vec));
% Theta_h_vec = zeros(size(delta_vec));

KB = 1.38e-16;
mproton = 1.67e-24;
Myr = 60*60*24*365.25*1e6;
kpc = 3.0856e21;
Zsolar = 0.02;

Mb = ( 200 / 185 ) * ( eta / sqrt(Theta_h) );
Zs = Zs0*Zsolar;     % absolute units
Zb = Zb0*Zsolar;   % absolute units
Zmix = sqrt(Zs*Zb); % absolute units

for i=1:length(delta_vec)    
    delta = delta_vec(i);
    n0(i) = ns_vec(i) * 0.01;         % Hydrogen number density of stream at Rv in cm^{-3}
    rhos = n0(i) .* 1.67e-24 ./ 0.76; % Stream density at Rv in gr/cm^3
    
    Rv = 100 .* M12.^(1/3) .* z3.^(-1.0);     % Virial radius in kpc, Dekel et al 2013
    Rs(i) = Rs_Rv_vec(i) .* Rv;               % Stream radius in kpc
    
    [Ts(i),mus] = find_cooling_equilibrium2(n0(i), red, Zs/Zsolar, 1);
    Tmin = 1.5*Ts(i);
    nmin = n0(i)/1.5;
    %%% Comment out to use prescribed Mb and/or theta_h values
%     Theta_s = Ts(i) / 1.5e4;
%     Theta_h_vec(i) = Theta_s / ( M12^(2/3) * z3 * (100/delta) );
%     Mb_vec(i) = ( 200 / 185 ) * ( eta / sqrt(Theta_h_vec(i)) );
%     Theta_h = Theta_h_vec(i);
%     Mb = Mb_vec(i);
    %%% Comment out to use UV background for Ts
    Theta_s = Theta_h * M12^(2/3) * z3 * (100/delta);
    Ts(i) = 1.5e4 * Theta_s;
    
    Cs_phys  = sqrt( (5/3) * KB * Ts(i) / (mus * mproton) ); % cm/s
    Cs_phys  = Cs_phys * Myr / kpc; % kpc/Myr
    tsc_phys = 2 * Rs(i) / Cs_phys; % Myr
    
    [nspec, mu, cool, heat, cool_tot, heat_tot] = cooling_equilibrium(Tmin, nmin, red, Zs/Zsolar, 1, 1);
    tcool = ( 1.5*KB*Tmin / (nmin*(cool_tot-heat_tot)) ) / Myr;    % Myr
    tcool_over_tsc = tcool / tsc_phys;
    tau_cool_norm = tcool_over_tsc/0.002;
    
    B(i) = 5 .* Rs_Rv_vec(i)^(-1) .* delta^(-1.5) .* Mb^(-1) .* tau_cool_norm^(-0.25);
    
    Lnorm_vec(i) = 2e40 .* M12^(5/3) .* z3^(-0.5) .* (Rs_Rv_vec(i)/0.16) .* ns_vec(i) .* ...
        (delta/100)^(-1.5) .* eta .* Mb^(-1) .* tau_cool_norm^(-0.25);    % Luminosity normalization in erg/sec
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xu1 = 1;
xl1 = 0.1;
dx1 = 0.1;

% yu1 = 3;
yu1 = 2.6;
yl1 = 0.8;
dy1 = 0.2;

% yu2 = 2;
% yl2 = 0.9;
% dy2 = 0.1;
yu2 = 1.4;
yl2 = 0.95;
dy2 = 0.05;

ytick3 = [1e39 1e40 1e41 1e42];
yu3 = 4*ytick3(end);
% yu3 = 3*ytick3(end);
yl3 = ytick3(1);

xtit = '$r/R_{\rm v}$';
ytit1 = '$|V|/V_{\rm v}$';
ytit2 = '$m/m_{\rm 0}$';
% ytit3 = '$\mathcal{L}_{\rm diss}(>r)\:{\rm [erg~s^{-1}]}$';
ytit3 = '$L_{\rm diss}(>r)\:{\rm [erg~s^{-1}]}$';

Mvexp = 12 + str2num(num2str(log10(M12),'%3.1f'));
tit = strcat('$M_{\rm v}=10^{',num2str(Mvexp,'%3.1f'),'}M_{\odot},\:z=',num2str(red,'%3.1f'),'$');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
set(figure1,'WindowStyle','docked');
set(gcf,'renderer','painters')

set(gcf,'visible','off');
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 30 9.5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axes1 = axes('Parent',figure1,...
    'YTick',yl1:dy1:yu1,...
    'YMinorTick','on',...
    'XTick',xl1:dx1:xu1,...
    'XMinorTick','on',...
    'TickLength',[0.02 0.04],...
    'PlotBoxAspectRatio',[1 1 1],...
    'LineWidth',1.5,...
    'FontSize',12,...
    'FontName','Arial',...
    'Position',[0.13 0.14 0.775 0.815]);
xlim(axes1,[xl1 xu1]);
ylim(axes1,[yl1 yu1]);
box(axes1,'on');
% grid on
hold(axes1,'all');
xlabel(xtit,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 -0.08 0]);
ylabel(ytit1,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[-0.11 0.5 0]);
title(tit,'Interpreter','latex','FontSize',13,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 0.87 0]);
annotation(figure1,'textbox',[0.12 0.76 0.10 0.05],...
    'String',strcat('$\eta=',num2str(eta,'%3.1f'),',\:\Theta_{\rm h}=',num2str(Theta_h,'%3.1f'),'$'),...
    'FontSize',14,...
    'FontName','Times',...
    'Interpreter','Latex',...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor',[1 1 1]);
annotation(figure1,'textbox',[0.135 0.68 0.10 0.05],...
    'String',strcat('$c=',num2str(c,'%2i'),',\:\beta=',num2str(beta,'%1i'),'$'),...
    'FontSize',14,...
    'FontName','Times',...
    'Interpreter','Latex',...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor',[1 1 1]);

[X,Y]=ode45(@(x,y)halo_decel_cooling_NFW(x, y, eta, 0, c, beta),[1 0.1],[eta^2, 1]);
l(1) = plot(X,sqrt(Y(:,1)),'linestyle','-','linewidth',2,'marker','none','color','k',...
    'DisplayName','${\rm free\:fall}$');
for i=1:length(delta_vec)
    if(i<=3)
%         dvec = strcat('\delta=',num2str(delta_vec(i),'%3i'),',\:\:\:');
        dvec = strcat('\delta=',num2str(delta_vec(i),'%3i'),',\,');
        col = [0, 1, 1];
    elseif(i<=6)
        dvec = strcat('\delta=',num2str(delta_vec(i),'%3i'),',\,');
        col = [1, 0, 0];
    elseif(i<=9)
        dvec = strcat('\delta=',num2str(delta_vec(i),'%3i'),',\,');
        col = [0, 0, 1];
    end
    if(mod(i,3)==1)
        ls = ':';
    elseif(mod(i,3)==2)
        ls = '-';
    elseif(mod(i,3)==0)
        ls = '--';
    end
%     rsvec = strcat('R_{\rm s}/R_{\rm v}=',num2str(Rs_Rv_vec(i),'%4.2f'),',\:');
    rsvec = strcat('R_{\rm sv}=',num2str(Rs_Rv_vec(i),'%4.2f'),',\,');
    nsvec = strcat('n_{\rm s,0.01}=',num2str(ns_vec(i),'%4.2f'));
    dname = strcat('$',dvec,rsvec,nsvec,'$');
    
    [X,Y]=ode45(@(x,y)halo_decel_cooling_NFW(x, y, eta, B(i), c, beta),[1 0.1],[eta^2, 1]);
    l(i+1) = plot(X,sqrt(Y(:,1)),'linestyle',ls,'linewidth',2,'marker','none','color',col,'DisplayName',dname);
end

legend1 = legend(axes1,l(1:4));
set(legend1,...
    'Position',[0.026, 0.17, 0.3, 0.15],'FontSize',10,'Interpreter','latex');
%     'Position',[0.024, 0.165, 0.3, 0.15],'FontSize',10,'Interpreter','latex');
legend boxoff

set(axes1, 'Position',[0.065 0.14 0.2467 0.78]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axes2 = axes('Parent',figure1,...
    'YTick',yl2:dy2:yu2,...
    'YMinorTick','on',...
    'XTick',xl1:dx1:xu1,...
    'XMinorTick','on',...
    'TickLength',[0.02 0.04],...
    'PlotBoxAspectRatio',[1 1 1],...
    'LineWidth',1.5,...
    'FontSize',12,...
    'FontName','Arial',...
    'Position',[0.13 0.14 0.775 0.815]);
xlim(axes2,[xl1 xu1]);
ylim(axes2,[yl2 yu2]);
box(axes2,'on');
% grid on
hold(axes2,'all');
xlabel(xtit,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 -0.08 0]);
ylabel(ytit2,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[-0.11 0.5 0]);

for i=1:length(delta_vec)
    if(i<=3)
%         dvec = strcat('\delta=',num2str(delta_vec(i),'%3i'),',\:\:\:');
        dvec = strcat('\delta=',num2str(delta_vec(i),'%3i'),',\,');
        col = [0, 1, 1];
    elseif(i<=6)
        dvec = strcat('\delta=',num2str(delta_vec(i),'%3i'),',\,');
        col = [1, 0, 0];
    elseif(i<=9)
        dvec = strcat('\delta=',num2str(delta_vec(i),'%3i'),',\,');
        col = [0, 0, 1];
    end
    if(mod(i,3)==1)
        ls = ':';
    elseif(mod(i,3)==2)
        ls = '-';
    elseif(mod(i,3)==0)
        ls = '--';
    end
%     rsvec = strcat('R_{\rm s}/R_{\rm v}=',num2str(Rs_Rv_vec(i),'%4.2f'),',\:');
    rsvec = strcat('R_{\rm sv}=',num2str(Rs_Rv_vec(i),'%4.2f'),',\,');
    nsvec = strcat('n_{\rm s,0.01}=',num2str(ns_vec(i),'%4.2f'));
    dname = strcat('$',dvec,rsvec,nsvec,'$');
    
    [X,Y]=ode45(@(x,y)halo_decel_cooling_NFW(x, y, eta, B(i), c, beta),[1 0.1],[eta^2,1]);
    l(i) = plot(X,Y(:,2),'linestyle',ls,'linewidth',2,'marker','none','color',col,'DisplayName',dname);
end

legend2 = legend(axes2,l(4:6));
set(legend2,...
    'Position',[0.34, 0.74, 0.3, 0.15],'FontSize',10,'Interpreter','latex');
legend boxoff

set(gca, 'Position',[0.365 0.14 0.2467 0.78]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axes3 = axes('Parent',figure1,...
    'YTick',ytick3,...
    'YScale','log',...
    'YMinorTick','on',...
    'XTick',xl1:dx1:xu1,...
    'XMinorTick','on',...
    'TickLength',[0.02 0.04],...
    'PlotBoxAspectRatio',[1 1 1],...
    'LineWidth',1.5,...
    'FontSize',12,...
    'FontName','Arial',...
    'Position',[0.13 0.14 0.775 0.815]);
xlim(axes3,[xl1 xu1]);
ylim(axes3,[yl3 yu3]);
box(axes3,'on');
% grid on
hold(axes3,'all');
xlabel(xtit,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 -0.08 0]);
ylabel(ytit3,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[-0.13 0.5 0]);

for i=1:length(delta_vec)
    if(i<=3)
%         dvec = strcat('\delta=',num2str(delta_vec(i),'%3i'),',\:\:\:');
        dvec = strcat('\delta=',num2str(delta_vec(i),'%3i'),',\,');
        col = [0, 1, 1];
    elseif(i<=6)
        dvec = strcat('\delta=',num2str(delta_vec(i),'%3i'),',\,');
        col = [1, 0, 0];
    elseif(i<=9)
        dvec = strcat('\delta=',num2str(delta_vec(i),'%3i'),',\,');
        col = [0, 0, 1];
    end
    if(mod(i,3)==1)
        ls = ':';
    elseif(mod(i,3)==2)
        ls = '-';
    elseif(mod(i,3)==0)
        ls = '--';
    end
%     rsvec = strcat('R_{\rm s}/R_{\rm v}=',num2str(Rs_Rv_vec(i),'%4.2f'),',\:');
    rsvec = strcat('R_{\rm sv}=',num2str(Rs_Rv_vec(i),'%4.2f'),',\,');
    nsvec = strcat('n_{\rm s,0.01}=',num2str(ns_vec(i),'%4.2f'));
    dname = strcat('$',dvec,rsvec,nsvec,'$');
    
    [X,Y]=ode45(@(x,y)halo_decel_cooling_NFW(x, y, eta, B(i), c, beta),[1 0.1],[eta^2,1]);
    EKdiss = abs( Y(:,1) );
    ETdiss = 1.8*eta^2/Mb^2 .* ones(size(EKdiss));
    
    % rs(r) derived from constant line mass
%     yint = ( ((5/3).*ETdiss + EKdiss) ./ (X.^(5*beta/8)) );    %yint = ( ( 1.8*eta^2/Mb_vec(i)^2 + EKdiss ) ./ (X.^1.25) );
    % rs(r) derived from self-consistent line mass
    yint = ( ((5/3).*ETdiss + EKdiss) ./ ( X.^(5*beta/8) .* Y(:,2).^(3/8) ) );    %yint = ( ( 1.8*eta^2/Mb_vec(i)^2 + EKdiss ) ./ (X.^1.25) );
    
    xq = linspace(X(1),X(end),1000);
    yq = interp1(X,log10(yint),xq);
    yq = 10.^yq;
    yplot = zeros(size(yq));
    for j=2:length(yq)
        yplot(j) = trapz(xq(j:-1:1),yq(j:-1:1));
    end
    yplot = Lnorm_vec(i) .* yplot;
    l(i) = plot(xq,yplot,'linestyle',ls,'linewidth',2,'marker','none','color',col,'DisplayName',dname);
end

legend3 = legend(axes3,l(7:9));
set(legend3,...
    'Position',[0.65, 0.74, 0.3, 0.15],'FontSize',10,'Interpreter','latex');
legend boxoff

set(gca, 'Position',[0.675 0.14 0.2467 0.78]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_direc = './decel_panels/';
mkdir(fig_direc);
fig_tit = strcat('toy_model_cooling_panels_c_',num2str(c,'%02i'),'_beta_',num2str(beta,'%1i'));
print(gcf,'-djpeg',strcat(fig_direc,fig_tit,'.jpg'));
print(gcf,'-depsc',strcat(fig_direc,fig_tit,'.eps'));
close all
fclose all;
