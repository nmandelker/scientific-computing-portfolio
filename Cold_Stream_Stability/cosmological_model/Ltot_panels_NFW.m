function Ltot_panels_NFW(c, beta, eta, Theta_h, Zs0, Zb0)

Mv = 10.^[11:0.1:14];
M12 = Mv./1e12;
red  = [1, 2, 4];
z3 = (1+red)./3;
Lmin = zeros(length(M12),length(red));
Lmax = zeros(length(M12),length(red));
Lfid = zeros(length(M12),length(red));

Vmin = zeros(length(M12),length(red));
Vmax = zeros(length(M12),length(red));
Vfid = zeros(length(M12),length(red));

Mumin = zeros(length(M12),length(red));
Mumax = zeros(length(M12),length(red));
Mufid = zeros(length(M12),length(red));

KB = 1.38e-16;
mproton = 1.67e-24;
Myr = 60*60*24*365.25*1e6;
kpc = 3.0856e21;
Zsolar = 0.02;

Mb = ( 200 / 185 ) * ( eta / sqrt(Theta_h) );
Zs = Zs0*Zsolar;     % absolute units
Zb = Zb0*Zsolar;   % absolute units
Zmix = sqrt(Zs*Zb); % absolute units

for nM=1:length(Mv)
    for nZ=1:length(red)
        
        Rv = 100 .* M12(nM).^(1/3) .* z3(nZ).^(-1.0);     % Virial radius in kpc, Dekel et al 2013
        
        delta_vec = [30, 30, 30, 100, 100, 100, 300, 300, 300] * M12(nM)^(2/3) * z3(nZ);
        Rs_Rv_vec = [0.09, 0.29, 0.50, 0.05, 0.16, 0.28, 0.03, 0.09, 0.16] * M12(nM)^(-1/3);
        ns_vec    = [0.45, 0.45, 0.15, 1.50, 1.50, 0.50, 4.50, 4.50, 1.50] * M12(nM)^(2/3) * z3(nZ)^4;

        Lnorm_vec = zeros(size(delta_vec));
        Vnorm_vec = zeros(size(delta_vec));
        Munorm_vec = zeros(size(delta_vec));
        
        for i=1:length(delta_vec)
            delta = delta_vec(i);
            n0 = ns_vec(i) * 0.01;         % Hydrogen number density of stream at Rv in cm^{-3}
            rhos = n0 .* 1.67e-24 ./ 0.76; % Stream density at Rv in gr/cm^3            
            Rs = Rs_Rv_vec(i) .* Rv;               % Stream radius in kpc
            
            [Ts,mus] = find_cooling_equilibrium2(n0, red(nZ), Zs/Zsolar, 1);
            T_mult = 1.5:0.1:5;
            tcool_vec = 0.*T_mult;
            for kT = 1:length(T_mult)
                Tmin = T_mult(kT) * Ts;
                nmin = n0 / T_mult(kT);
%                 nmin = n0;
                [nspec, mu, cool, heat, cool_tot, heat_tot] = ...
                    cooling_equilibrium(Tmin, nmin, red(nZ), Zs/Zsolar, 1, 1);
                tcool_vec(kT) = ( 1.5*KB*Tmin / (nmin*(cool_tot-heat_tot)) ) / Myr;    % Myr
            end
            b = find(tcool_vec>0);
            tcool = min(tcool_vec(b));
                
%             Tmin = 3*Ts;
% %             nmin = n0/1.5;
% %             Tmin = 3e4;
%             nmin = n0;
%             [nspec, mu, cool, heat, cool_tot, heat_tot] = ...
%                 cooling_equilibrium(Tmin, nmin, red(nZ), Zs/Zsolar, 1, 1);
%             tcool = ( 1.5*KB*Tmin / (nmin*(cool_tot-heat_tot)) ) / Myr;    % Myr
%             while(tcool<0)
%                 Tmin = 1.5*Tmin;
%                 nmin = nmin/1.5;
%                 [nspec, mu, cool, heat, cool_tot, heat_tot] = ...
%                     cooling_equilibrium(Tmin, nmin, red(nZ), Zs/Zsolar, 1, 1);
%                 tcool = ( 1.5*KB*Tmin / (nmin*(cool_tot-heat_tot)) ) / Myr;    % Myr
%             end
            
            %%% Comment out to use UV background for Ts
            Theta_s = Theta_h * M12(nM)^(2/3) * z3(nZ) * (100/delta);
            Ts = 1.5e4 * Theta_s;
            
            Cs_phys  = sqrt( (5/3) * KB * Ts / (mus * mproton) ); % cm/s
            Cs_phys  = Cs_phys * Myr / kpc; % kpc/Myr
            tsc_phys = 2 * Rs / Cs_phys; % Myr
            
            tcool_over_tsc = tcool / tsc_phys;
            tau_cool_norm = tcool_over_tsc/0.002;
            
            B = 5 .* Rs_Rv_vec(i)^(-1) .* delta^(-1.5) .* Mb^(-1) .* tau_cool_norm^(-0.25);
            % Luminosity normalization in erg/sec
            Lnorm = 2e40 .* M12(nM)^(5/3) .* z3(nZ)^(-0.5) .* (Rs_Rv_vec(i)/0.16) .* ...
                ns_vec(i) .* (delta/100)^(-1.5) .* eta .* Mb^(-1) .* tau_cool_norm^(-0.25);
            
            
            [X,Y]=ode45(@(x,y)halo_decel_cooling_NFW(x, y, eta, B, c, beta),[1 0.1],[eta^2,1]);
            EKdiss = abs( Y(:,1) );
            ETdiss = 1.8*eta^2/Mb^2 .* ones(size(EKdiss));
            
            % rs(r) derived from constant line mass
%             yint = ( ( (5/3)*ETdiss + EKdiss ) ./ (X.^(5*beta/8)) );
            % rs(r) derived from self-consistent line mass
            yint = ( ( (5/3)*ETdiss + EKdiss ) ./ ( X.^(5*beta/8) .* Y(:,2).^(3/8) ) );
            
            xq = linspace(X(1),X(end),1000);
            yq = interp1(X,log10(yint),xq);
            yq = 10.^yq;
            yplot = zeros(size(yq));
            for j=2:length(yq)
                yplot(j) = trapz(xq(j:-1:1),yq(j:-1:1));
            end
            yplot = Lnorm .* yplot;
            Lnorm_vec(i) = yplot(end);            
            Vnorm_vec(i) = sqrt(abs(Y(end,1)));
            Munorm_vec(i) = Y(end,2);
        end
        
        Lmin(nM,nZ) = min(Lnorm_vec);
        Lmax(nM,nZ) = max(Lnorm_vec);
        Lfid(nM,nZ) = Lnorm_vec(5);
        
        Vmin(nM,nZ) = min(Vnorm_vec);
        Vmax(nM,nZ) = max(Vnorm_vec);
        Vfid(nM,nZ) = Vnorm_vec(5);
        
        Mumin(nM,nZ) = min(Munorm_vec);
        Mumax(nM,nZ) = max(Munorm_vec);
        Mufid(nM,nZ) = Munorm_vec(5);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xtit = '$M_{\rm v}\:{\rm [M_{\odot}] }$';
xtick = [1e11, 1e12, 1e13, 1e14];
xstr = {'10^{11}', '10^{12}', '10^{13}', '10^{14}'};
xu = xtick(end);
xl = xtick(1);

ytit1 = '$|V|(0.1R_{\rm v})/V_{\rm v}$';
% yu1 = 3;
yu1 = 2.7;
yl1 = 2;
dy1 = 0.1;
ytick1 = yl1:dy1:yu1;

ytit2 = '$m(0.1R_{\rm v})/m_{\rm 0}$';
% yu2 = 4;
yu2 = 2;
yl2 = 0.9;
% dy2 = 0.5;
dy2 = 0.2;
ytick2 = yu2:-dy2:yl2;
ytick2 = ytick2(end:-1:1);

ytit3 = '$L_{\rm diss}(>0.1R_{\rm v})\:{\rm [erg~s^{-1}]}$';
ytick3 = [1e39 1e40 1e41 1e42 1e43 1e44 1e45];
yu3 = ytick3(end);
yl3 = ytick3(1);

tit = strcat('$\Theta_{\rm h}=',num2str(Theta_h,'%3.1f'),',\:\eta=',num2str(eta,'%3.1f'),'$');

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
    'YTick',ytick1,...
    'YMinorTick','on',...
    'XTick',xtick,...
    'XtickLabel',xstr,...
    'XMinorTick','on',...
    'XScale','log',...
    'TickLength',[0.02 0.04],...
    'PlotBoxAspectRatio',[1 1 1],...
    'LineWidth',1.5,...
    'FontSize',12,...
    'FontName','Arial',...
    'Position',[0.13 0.14 0.775 0.815]);
xlim(axes1,[xl xu]);
ylim(axes1,[yl1 yu1]);
box(axes1,'on');
% grid on
hold(axes1,'all');
xlabel(xtit,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 -0.09 0]);
ylabel(ytit1,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[-0.12 0.5 0]);

n = floor(length(Mv)/2);

l(2) = plot(Mv,Vfid(:,1),'marker','none','linestyle','-','linewidth',2,'color','g',...
    'DisplayName',strcat('$z=',num2str(red(1),'%4.2f'),'$'));
plot(Mv,Vmin(:,1),'marker','none','linestyle',':','linewidth',1,'color','g');
plot(Mv,Vmax(:,1),'marker','none','linestyle',':','linewidth',1,'color','g');
y = linspace(Vmin(n-2,1),Vmax(n-2,1),10);
x = Mv(n-2).*ones(size(y));
plot(x,y,'marker','none','linestyle','-','linewidth',1.5,'color','g');
plot(Mv(n-3:n-1),Vmin(n-3:n-1,1),'marker','none','linestyle','-','linewidth',1.5,'color','g');
plot(Mv(n-3:n-1),Vmax(n-3:n-1,1),'marker','none','linestyle','-','linewidth',1.5,'color','g');

l(3) = plot(Mv,Vfid(:,2),'marker','none','linestyle','-','linewidth',2,'color','b',...
    'DisplayName',strcat('$z=',num2str(red(2),'%4.2f'),'$'));
plot(Mv,Vmin(:,2),'marker','none','linestyle',':','linewidth',1,'color','b');
plot(Mv,Vmax(:,2),'marker','none','linestyle',':','linewidth',1,'color','b');
y = linspace(Vmin(n,2),Vmax(n,2),10);
x = Mv(n).*ones(size(y));
plot(x,y,'marker','none','linestyle','-','linewidth',1.5,'color','b');
plot(Mv(n-1:n+1),Vmin(n-1:n+1,2),'marker','none','linestyle','-','linewidth',1.5,'color','b');
plot(Mv(n-1:n+1),Vmax(n-1:n+1,2),'marker','none','linestyle','-','linewidth',1.5,'color','b');

l(4) = plot(Mv,Vfid(:,3),'marker','none','linestyle','-','linewidth',2,'color','r',...
    'DisplayName',strcat('$z=',num2str(red(3),'%4.2f'),'$'));
plot(Mv,Vmin(:,3),'marker','none','linestyle',':','linewidth',1,'color','r');
plot(Mv,Vmax(:,3),'marker','none','linestyle',':','linewidth',1,'color','r');
y = linspace(Vmin(n+2,3),Vmax(n+2,3),10);
x = Mv(n+2).*ones(size(y));
plot(x,y,'marker','none','linestyle','-','linewidth',1.5,'color','r');
plot(Mv(n+1:n+3),Vmin(n+1:n+3,3),'marker','none','linestyle','-','linewidth',1.5,'color','r');
plot(Mv(n+1:n+3),Vmax(n+1:n+3,3),'marker','none','linestyle','-','linewidth',1.5,'color','r');

% for i=1:(length(Mv)-2)
%     y = linspace(Vmin(i+1,1),Vmax(i+1,1),10);
%     x = Mv(i+1).*ones(size(y));
%     plot(x,y,'marker','none','linestyle','-','linewidth',0.5,'color','g');
%     
%     x = log10(linspace(Mv(i),Mv(i+2),10));
%     y = log10(Vmin(i,2)) + (log10(Vmax(i+2,2))-log10(Vmin(i,2)))./...
%         (log10(Mv(i+2))-log10(Mv(i))).*(x-log10(Mv(i)));
%     plot(10.^x,10.^y,'marker','none','linestyle','-','linewidth',0.5,'color','b');
%     
%     x = log10(linspace(Mv(i),Mv(i+2),10));
%     y = log10(Vmax(i,3)) + (log10(Vmin(i+2,3))-log10(Vmax(i,3)))./...
%         (log10(Mv(i+2))-log10(Mv(i))).*(x-log10(Mv(i)));
%     plot(10.^x,10.^y,'marker','none','linestyle','-','linewidth',0.5,'color','r');
% end

[X,Y]=ode45(@(x,y)halo_decel_cooling_NFW(x, y, eta, 0, c, beta),[1 0.1],[eta^2,1]);
l(1) = plot(Mv,sqrt(Y(end,1)).*ones(size(Mv)),'linestyle','-','linewidth',2,'marker','none','color','k',...
    'DisplayName','${\rm free\:fall}$');

legend1 = legend(axes1,l(1:4));
set(legend1,...
    'Location','SouthEast','FontSize',12,'Interpreter','latex');
legend boxoff

set(gca, 'Position',[0.065 0.14 0.2467 0.78]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axes2 = axes('Parent',figure1,...
    'YTick',ytick2,...    'YtickLabel',ystr,...
    'YMinorTick','on',...
    'XTick',xtick,...
    'XtickLabel',xstr,...
    'XMinorTick','on',...
    'XScale','log',...
    'TickLength',[0.02 0.04],...
    'PlotBoxAspectRatio',[1 1 1],...
    'LineWidth',1.5,...
    'FontSize',12,...
    'FontName','Arial',...
    'Position',[0.13 0.14 0.775 0.815]);
xlim(axes2,[xl xu]);
ylim(axes2,[yl2 yu2]);
box(axes2,'on');
% grid on
hold(axes2,'all');
xlabel(xtit,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 -0.09 0]);
ylabel(ytit2,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[-0.11 0.5 0]);

annotation(figure1,'textbox',[0.42 0.84 0.10 0.05],...
    'String',strcat('$\eta=',num2str(eta,'%3.1f'),',\:\Theta_{\rm h}=',num2str(Theta_h,'%3.1f'),'$'),...
    'FontSize',14,...
    'FontName','Times',...
    'Interpreter','Latex',...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor',[1 1 1]);
annotation(figure1,'textbox',[0.435 0.76 0.10 0.05],...
    'String',strcat('$c=',num2str(c,'%2i'),',\:\beta=',num2str(beta,'%1i'),'$'),...
    'FontSize',14,...
    'FontName','Times',...
    'Interpreter','Latex',...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor',[1 1 1]);

l(1) = plot(Mv,Mufid(:,1),'marker','none','linestyle','-','linewidth',2,'color','g',...
    'DisplayName',strcat('$z=',num2str(red(1),'%4.2f'),'$'));
plot(Mv,Mumin(:,1),'marker','none','linestyle',':','linewidth',1,'color','g');
plot(Mv,Mumax(:,1),'marker','none','linestyle',':','linewidth',1,'color','g');
y = linspace(Mumin(n-2,1),Mumax(n-2,1),10);
x = Mv(n-2).*ones(size(y));
plot(x,y,'marker','none','linestyle','-','linewidth',1.5,'color','g');
plot(Mv(n-3:n-1),Mumin(n-3:n-1,1),'marker','none','linestyle','-','linewidth',1.5,'color','g');
plot(Mv(n-3:n-1),Mumax(n-3:n-1,1),'marker','none','linestyle','-','linewidth',1.5,'color','g');

l(2) = plot(Mv,Mufid(:,2),'marker','none','linestyle','-','linewidth',2,'color','b',...
    'DisplayName',strcat('$z=',num2str(red(2),'%4.2f'),'$'));
plot(Mv,Mumin(:,2),'marker','none','linestyle',':','linewidth',1,'color','b');
plot(Mv,Mumax(:,2),'marker','none','linestyle',':','linewidth',1,'color','b');
y = linspace(Mumin(n,2),Mumax(n,2),10);
x = Mv(n).*ones(size(y));
plot(x,y,'marker','none','linestyle','-','linewidth',1.5,'color','b');
plot(Mv(n-1:n+1),Mumin(n-1:n+1,2),'marker','none','linestyle','-','linewidth',1.5,'color','b');
plot(Mv(n-1:n+1),Mumax(n-1:n+1,2),'marker','none','linestyle','-','linewidth',1.5,'color','b');

l(3) = plot(Mv,Mufid(:,3),'marker','none','linestyle','-','linewidth',2,'color','r',...
    'DisplayName',strcat('$z=',num2str(red(3),'%4.2f'),'$'));
plot(Mv,Mumin(:,3),'marker','none','linestyle',':','linewidth',1,'color','r');
plot(Mv,Mumax(:,3),'marker','none','linestyle',':','linewidth',1,'color','r');
y = linspace(Mumin(n+2,3),Mumax(n+2,3),10);
x = Mv(n+2).*ones(size(y));
plot(x,y,'marker','none','linestyle','-','linewidth',1.5,'color','r');
plot(Mv(n+1:n+3),Mumin(n+1:n+3,3),'marker','none','linestyle','-','linewidth',1.5,'color','r');
plot(Mv(n+1:n+3),Mumax(n+1:n+3,3),'marker','none','linestyle','-','linewidth',1.5,'color','r');

% for i=1:(length(Mv)-2)
%     y = linspace(Mumin(i+1,1),Mumax(i+1,1),10);
%     x = Mv(i+1).*ones(size(y));
%     plot(x,y,'marker','none','linestyle','-','linewidth',0.5,'color','g');
%     
%     x = log10(linspace(Mv(i),Mv(i+2),10));
%     y = log10(Mumin(i,2)) + (log10(Mumax(i+2,2))-log10(Mumin(i,2)))./...
%         (log10(Mv(i+2))-log10(Mv(i))).*(x-log10(Mv(i)));
%     plot(10.^x,10.^y,'marker','none','linestyle','-','linewidth',0.5,'color','b');
%     
%     x = log10(linspace(Mv(i),Mv(i+2),10));
%     y = log10(Mumax(i,3)) + (log10(Mumin(i+2,3))-log10(Mumax(i,3)))./...
%         (log10(Mv(i+2))-log10(Mv(i))).*(x-log10(Mv(i)));
%     plot(10.^x,10.^y,'marker','none','linestyle','-','linewidth',0.5,'color','r');
% end

% legend2 = legend(axes2,l(1:3));
% set(legend2,...
%     'Location','NorthEast','FontSize',12,'Interpreter','latex');
% legend boxoff

set(gca, 'Position',[0.365 0.14 0.2467 0.78]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axes3 = axes('Parent',figure1,...
    'YTick',ytick3,...    'YtickLabel',ystr,...
    'YMinorTick','on',...
    'YScale','log',...
    'XTick',xtick,...
    'XtickLabel',xstr,...
    'XMinorTick','on',...
    'XScale','log',...
    'TickLength',[0.02 0.04],...
    'PlotBoxAspectRatio',[1 1 1],...
    'LineWidth',1.5,...
    'FontSize',12,...
    'FontName','Arial',...
    'Position',[0.13 0.14 0.775 0.815]);
xlim(axes3,[xl xu]);
ylim(axes3,[yl3 yu3]);
box(axes3,'on');
% grid on
hold(axes3,'all');
xlabel(xtit,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 -0.09 0]);
ylabel(ytit3,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[-0.13 0.5 0]);

l(1) = plot(Mv,Lfid(:,1),'marker','none','linestyle','-','linewidth',2,'color','g',...
    'DisplayName',strcat('$z=',num2str(red(1),'%4.2f'),'$'));
plot(Mv,Lmin(:,1),'marker','none','linestyle',':','linewidth',1,'color','g');
plot(Mv,Lmax(:,1),'marker','none','linestyle',':','linewidth',1,'color','g');
y = linspace(Lmin(n-2,1),Lmax(n-2,1),10);
x = Mv(n-2).*ones(size(y));
plot(x,y,'marker','none','linestyle','-','linewidth',1.5,'color','g');
plot(Mv(n-3:n-1),Lmin(n-3:n-1,1),'marker','none','linestyle','-','linewidth',1.5,'color','g');
plot(Mv(n-3:n-1),Lmax(n-3:n-1,1),'marker','none','linestyle','-','linewidth',1.5,'color','g');

l(2) = plot(Mv,Lfid(:,2),'marker','none','linestyle','-','linewidth',2,'color','b',...
    'DisplayName',strcat('$z=',num2str(red(2),'%4.2f'),'$'));
plot(Mv,Lmin(:,2),'marker','none','linestyle',':','linewidth',1,'color','b');
plot(Mv,Lmax(:,2),'marker','none','linestyle',':','linewidth',1,'color','b');
y = linspace(Lmin(n,2),Lmax(n,2),10);
x = Mv(n).*ones(size(y));
plot(x,y,'marker','none','linestyle','-','linewidth',1.5,'color','b');
plot(Mv(n-1:n+1),Lmin(n-1:n+1,2),'marker','none','linestyle','-','linewidth',1.5,'color','b');
plot(Mv(n-1:n+1),Lmax(n-1:n+1,2),'marker','none','linestyle','-','linewidth',1.5,'color','b');

l(3) = plot(Mv,Lfid(:,3),'marker','none','linestyle','-','linewidth',2,'color','r',...
    'DisplayName',strcat('$z=',num2str(red(3),'%4.2f'),'$'));
plot(Mv,Lmin(:,3),'marker','none','linestyle',':','linewidth',1,'color','r');
plot(Mv,Lmax(:,3),'marker','none','linestyle',':','linewidth',1,'color','r');
y = linspace(Lmin(n+2,3),Lmax(n+2,3),10);
x = Mv(n+2).*ones(size(y));
plot(x,y,'marker','none','linestyle','-','linewidth',1.5,'color','r');
plot(Mv(n+1:n+3),Lmin(n+1:n+3,3),'marker','none','linestyle','-','linewidth',1.5,'color','r');
plot(Mv(n+1:n+3),Lmax(n+1:n+3,3),'marker','none','linestyle','-','linewidth',1.5,'color','r');

% for i=1:(length(Mv)-2)
%     y = linspace(Lmin(i+1,1),Lmax(i+1,1),10);
%     x = Mv(i+1).*ones(size(y));
%     plot(x,y,'marker','none','linestyle','-','linewidth',0.5,'color','g');
%     
%     x = log10(linspace(Mv(i),Mv(i+2),10));
%     y = log10(Lmin(i,2)) + (log10(Lmax(i+2,2))-log10(Lmin(i,2)))./...
%         (log10(Mv(i+2))-log10(Mv(i))).*(x-log10(Mv(i)));
%     plot(10.^x,10.^y,'marker','none','linestyle','-','linewidth',0.5,'color','b');
%     
%     x = log10(linspace(Mv(i),Mv(i+2),10));
%     y = log10(Lmax(i,3)) + (log10(Lmin(i+2,3))-log10(Lmax(i,3)))./...
%         (log10(Mv(i+2))-log10(Mv(i))).*(x-log10(Mv(i)));
%     plot(10.^x,10.^y,'marker','none','linestyle','-','linewidth',0.5,'color','r');
% end

% legend3 = legend(axes3,l(1:3));
% set(legend3,...
%     'Location','NorthWest','FontSize',12,'Interpreter','latex');
% legend boxoff

set(gca, 'Position',[0.675 0.14 0.2467 0.78]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_direc = './Ltot_panels/';
mkdir(fig_direc);
fig_tit = strcat('Ltot_panels_c_',num2str(c,'%02i'),'_beta_',num2str(beta,'%1i'));
print(gcf,'-djpeg',strcat(fig_direc,fig_tit,'.jpg'));
print(gcf,'-depsc',strcat(fig_direc,fig_tit,'.eps'));
close all
fclose all;