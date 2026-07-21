
KB = 1.38e-16;
mproton = 1.67e-24;
Myr = 60*60*24*365.25*1e6;
kpc = 3.086e21;
Zsolar = 0.02;

% Zs = 0.001*Zsolar;  % absolute units
% Zb = 0.010*Zsolar;  % absolute units
Zs = 0.03*Zsolar;  % absolute units
Zb = 0.10*Zsolar;  % absolute units
Zmix = sqrt(Zs*Zb);

red = 2;    % for UVB
Rs0 = 3;     % kpc
nb  = 1e-4;  % cm^{-3}

delta_i = [10,  50,  50,  50];
delta_f = [100, 100, 500, 1000];
Nparam = length(delta_i);
Nd = 1000;

ns_i = 0.*delta_i;
ns_f = 0.*delta_i;
Ts_f = 0.*delta_i;
mus_f = 0.*delta_i;
T2s_f = 0.*delta_i;
Tb = 0.*delta_i;
mub = 0.*delta_i;
T2b = 0.*delta_i;

delta = zeros(Nd,Nparam);
ns    = zeros(Nd,Nparam);
T2s   = zeros(Nd,Nparam);
Ts    = zeros(Nd,Nparam);
mus   = zeros(Nd,Nparam);
Cs    = zeros(Nd,Nparam);
tsc   = zeros(Nd,Nparam);
tcool = zeros(Nd,Nparam);
Rs    = zeros(Nd,Nparam);
for j=1:Nparam
    ns_i(j)  = delta_i(j) * nb;
    ns_f(j)  = delta_f(j) * nb;
    [Ts_f(j),mus_f(j)] = find_cooling_equilibrium2(ns_f(j), red, Zs/Zsolar, 1);
    T2s_f(j) = Ts_f(j)/mus_f(j);
    T2b(j) = T2s_f(j) * delta_f(j);
    [nspec, mub(j), cool, heat, cool_tot, heat_tot] = cooling_equilibrium_T2_input(T2b(j), nb, red, Zb/Zsolar, 1, 1);
    Tb(j) = T2b(j) * mub(j);
    
    % units_density = 2.184210526e-28;
    % units_length  = 2.9568e23;
    % if( delta_f==100 )
    %     units_time = 2.0616635e16;
    % elseif( delta_f==500 )
    %     units_time = 9.8711656645e15;
    % end
    
    delta(:,j) = 10.^( linspace(log10(delta_i(j)),log10(delta_f(j)),Nd) );
    ns(:,j)    = nb .* delta(:,j);
    T2s(:,j)   = T2s_f(j) .* (delta_f(j) ./ delta(:,j));
    
    for i=1:Nd
        [nspec, mus(i,j), cool, heat, cool_tot, heat_tot] = cooling_equilibrium_T2_input(T2s(i,j), ns(i,j), red, Zs/Zsolar, 1, 1);
        Ts(i,j)    = T2s(i,j) * mus(i,j);
        Rs(i,j)    = Rs0 .* ( (ns(i,j)./ns_i(j)).^(-0.5) );
        Cs(i,j)    = sqrt( (5/3) * KB * T2s(i,j) / mproton ); % cm/s
        Cs(i,j)    = Cs(i,j) * Myr / kpc; % kpc/Myr
        tsc(i,j)   = 2 * Rs(i,j) / Cs(i,j);    % Myr
        tcool(i,j) = ( 1.5*KB*Ts(i,j) / (ns(i,j)*(cool_tot-heat_tot)) ) / Myr;    % Myr
    end
end

xtick = [1, 10, 100, 1000];
xstr  = {'1','10','100','1000'};
xl = xtick(1);
xu = xtick(end);

ytick = [0.01, 0.1, 1, 10, 100, 1000];
ystr  = {'0.01','0.1','1','10','100', '1000'};
yl = ytick(1);
yu = ytick(end);


% tit = strcat('$z=',num2str(red,'%1i'),',\,R_{\rm s,i}=',num2str(Rs0,'%3.1f'),'{\rm kpc},\,Z_{\rm s}=',num2str(Zs/Zsolar,'%4.2f'),...
%     'Z_{\odot},\,n_{\rm H,b}=',num2str(nb,'%3.1e'),'{\rm cm^{-3}}$');
tit = strcat('$R_{\rm s,i}=',num2str(Rs0,'%3.1f'),'{\rm kpc},\,Z_{\rm s}=',num2str(Zs/Zsolar,'%4.2f'),...
    'Z_{\odot},\,n_{\rm H,b}=',num2str(nb,'%3.1e'),'{\rm cm^{-3}}$');
xtit = '$\delta$';
ytit = '$R_{\rm s}\,/\,c_{\rm s}t_{\rm cool}$';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
set(figure1,'WindowStyle','docked');
set(gcf,'renderer','painters')

set(gcf,'visible','off');
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 10 9.5]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes1 = axes('Parent',figure1,'YTick',ytick,...
    'YMinorTick','on',...
    'YScale','log',...
    'YtickLabel',ystr,...
    'XTick',xtick,...
    'XMinorTick','on',...
    'XScale','log',...
    'XtickLabel',xstr,...
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
    'position',[0.5 -0.07 0]);
ylabel(ytit,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[-0.09 0.5 0]);
title(tit,'Interpreter','latex','FontSize',12,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 1.02 0]);

col(:,1) = [0,0,0];
col(:,2) = [1,0,0];
col(:,3) = [0,0,1];
col(:,4) = [0,1,0];
for j=1:Nparam
    DisplayName = strcat('$\delta_{\rm i}=',num2str(delta_i(j),'%4i'),...
        ',\,\delta_{\rm f}=',num2str(delta_f(j),'%4i'),',\,T_{\rm b}=',num2str(Tb(j),'%3.1e'),'{\rm K}$');
    l(j) = plot(delta(:,j),Rs(:,j)./(Cs(:,j).*tcool(:,j)),'linestyle','-','linewidth',2,'color',col(:,j),'DisplayName',DisplayName);
end
legend1 = legend(axes1,l(1:Nparam));
set(legend1,...
    'Location','SouthWest','FontSize',10,'Interpreter','latex');
legend boxoff

set(axes1, 'Position',[0.12 0.12 0.8 0.8]);

filename = strcat('./Rs_over_Cs_tcool/Rs0_',num2str(Rs0,'%3.1f'),'_Zs_',num2str(Zs/Zsolar,'%4.2f'),'_nb_',num2str(log10(nb),'%4.1f'));
print(gcf,'-djpeg',strcat(filename,'.jpg'));
% print(gcf,'-depsc',strcat(filename,'.eps'));