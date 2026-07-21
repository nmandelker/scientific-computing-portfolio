clear

Theta_h  = 1;
eta = 1;
Mb = ( 200 / 185 ) * ( eta / sqrt(Theta_h) );

Zs0 = 0.03;
Zb0 = 0.1;

Mv = 10.^[11:0.1:14];
M12 = Mv./1e12;
red  = [1, 2, 4];
z3 = (1+red)./3;
Lmin = zeros(length(M12),length(red));
Lmax = zeros(length(M12),length(red));
Lfid = zeros(length(M12),length(red));

KB = 1.38e-16;
mproton = 1.67e-24;
Myr = 60*60*24*365.25*1e6;
kpc = 3.0856e21;
Zsolar = 0.02;

for nM=1:length(Mv)
    for nZ=1:length(red)
        
        Rv = 100 .* M12(nM).^(1/3) .* z3(nZ).^(-1.0);     % Virial radius in kpc, Dekel et al 2013
        Vv = 200 .* 1e5 .* M12(nM).^(1/3) .* z3(nZ).^0.5; % Virial velocity in cm/sec, Dekel et al 2013
        
        delta_vec = [30, 30, 30, 100, 100, 100, 300, 300, 300] * M12(nM)^(2/3) * z3(nZ);
        Rs_Rv_vec = [0.09, 0.29, 0.50, 0.05, 0.16, 0.28, 0.03, 0.09, 0.16] * M12(nM)^(-1/3);
        ns_vec    = [0.45, 0.45, 0.15, 1.50, 1.50, 0.50, 4.50, 4.50, 1.50] * M12(nM)^(2/3) * z3(nZ)^4;

        Lnorm_vec = zeros(size(delta_vec));
        
        for i=1:length(delta_vec)
            delta = delta_vec(i);
            n0 = ns_vec(i) * 0.01;         % Hydrogen number density of stream at Rv in cm^{-3}
            rhos = n0 .* 1.67e-24 ./ 0.76; % Stream density at Rv in gr/cm^3            
            Rs = Rs_Rv_vec(i) .* Rv;               % Stream radius in kpc            
            m0 = pi .* (Rs*kpc)^2 .* (Rv*kpc) .* rhos; % initial stream mass in gr
            
            Zs = Zs0*Zsolar;     % absolute units
            Zb = Zb0*Zsolar;   % absolute units
            Zmix = sqrt(Zs*Zb); % absolute units
            
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
            
            
            [X,Y]=ode45(@(x,y)halo_decel_cooling(x, y, B),[1 0.1],[1,1]);
            EKdiss = abs( Y(:,1) + 2.*( 1 - log(X(:)) ) );
            yint = ( ( 1.8*eta^2/Mb^2 + EKdiss ) ./ (X.^1.25) );
            xq = linspace(X(1),X(end),1000);
            yq = interp1(X,log10(yint),xq);
            yq = 10.^yq;
            yplot = zeros(size(yq));
            for j=2:length(yq)
                yplot(j) = trapz(xq(j:-1:1),yq(j:-1:1));
            end
            yplot = Lnorm .* yplot;
            Lnorm_vec(i) = yplot(end);
        end
        
        Lmin(nM,nZ) = min(Lnorm_vec);
        Lmax(nM,nZ) = max(Lnorm_vec);
        Lfid(nM,nZ) = Lnorm_vec(5);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xtit = '$M_{\rm v}\:{\rm [M_{\odot}] }$';
xtick = [1e11, 1e12, 1e13, 1e14];
xstr = {'10^{11}', '10^{12}', '10^{13}', '10^{14}'};
xu = xtick(end);
xl = xtick(1);

ytit = '$L_{\rm diss}(>0.1R_{\rm v})\:{\rm [erg~s^{-1}]}$';
ytick = [1e40 1e41 1e42 1e43 1e44];
yu = 2*ytick(end);
yl = ytick(1);

tit = strcat('$\Theta_{\rm h}=',num2str(Theta_h,'%3.1f'),',\:\eta=',num2str(eta,'%3.1f'),'$');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
set(figure1,'WindowStyle','docked');
set(gcf,'renderer','painters')

set(gcf,'visible','off');
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 10 9.5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axes1 = axes('Parent',figure1,...
    'YTick',ytick,...    'YtickLabel',ystr,...
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
xlim(axes1,[xl xu]);
ylim(axes1,[yl yu]);
box(axes1,'on');
% grid on
hold(axes1,'all');
xlabel(xtit,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 -0.08 0]);
ylabel(ytit,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[-0.13 0.5 0]);

l(1) = plot(Mv,Lfid(:,1),'marker','none','linestyle','-','linewidth',2,'color','g',...
    'DisplayName',strcat('$z=',num2str(red(1),'%4.2f'),'$'));
plot(Mv,Lmin(:,1),'marker','none','linestyle','-','linewidth',1,'color','g');
plot(Mv,Lmax(:,1),'marker','none','linestyle','-','linewidth',1,'color','g');

l(2) = plot(Mv,Lfid(:,2),'marker','none','linestyle','-','linewidth',2,'color','b',...
    'DisplayName',strcat('$z=',num2str(red(2),'%4.2f'),'$'));
plot(Mv,Lmin(:,2),'marker','none','linestyle','-','linewidth',1,'color','b');
plot(Mv,Lmax(:,2),'marker','none','linestyle','-','linewidth',1,'color','b');

l(3) = plot(Mv,Lfid(:,3),'marker','none','linestyle','-','linewidth',2,'color','r',...
    'DisplayName',strcat('$z=',num2str(red(3),'%4.2f'),'$'));
plot(Mv,Lmin(:,3),'marker','none','linestyle','-','linewidth',1,'color','r');
plot(Mv,Lmax(:,3),'marker','none','linestyle','-','linewidth',1,'color','r');
for i=1:(length(Mv)-2)
    y = linspace(Lmin(i+1,1),Lmax(i+1,1),10);
    x = Mv(i+1).*ones(size(y));
    plot(x,y,'marker','none','linestyle','-','linewidth',0.5,'color','g');
    
    x = log10(linspace(Mv(i),Mv(i+2),10));
    y = log10(Lmin(i,2)) + (log10(Lmax(i+2,2))-log10(Lmin(i,2)))./...
        (log10(Mv(i+2))-log10(Mv(i))).*(x-log10(Mv(i)));
    plot(10.^x,10.^y,'marker','none','linestyle','-','linewidth',0.5,'color','b');
    
    x = log10(linspace(Mv(i),Mv(i+2),10));
    y = log10(Lmax(i,3)) + (log10(Lmin(i+2,3))-log10(Lmax(i,3)))./...
        (log10(Mv(i+2))-log10(Mv(i))).*(x-log10(Mv(i)));
    plot(10.^x,10.^y,'marker','none','linestyle','-','linewidth',0.5,'color','r');
end

legend1 = legend(axes1,l(1:3));
set(legend1,...
    'Location','NorthWest','FontSize',10,'Interpreter','latex');
legend boxoff

set(gca, 'Position',[0.20 0.14 0.74 0.78]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print(gcf,'-djpeg',strcat('./Ltot.jpg'));
print(gcf,'-depsc',strcat('./Ltot.eps'));
close all
fclose all;