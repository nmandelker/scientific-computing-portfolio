clear

Thmin = 3/8;
Thmax = 1;
Thfid = 1;

Tsmin = 0.5;
Tsmax = sqrt(2.0);
Tsfid = 1;

Tmin = Thmin/Tsmax;
Tmax = Thmax/Tsmin;
% Tmin = 0.3;
% Tmax = 3.0;
Tfid = 1.0;

eta_min = 0.5;
eta_max = sqrt(2);
eta_fid = 1.0;
% eta_min = 1.0;
% eta_max = 1.0;

fmin = 1.0;
fmax = 3.0;
ffid = 1.0;

smin = 0.3;
smax = 3.0;
sfid = 1.0;

alpha_min = 0.5;
alpha_max = 1.0;
alpha_fid = 1.0;

Mb_min = 1.08 * eta_min / sqrt(Thmax);
Mb_max = 1.08 * eta_max / sqrt(Thmin);
Mb_fid = 1.08 * eta_fid / sqrt(Thfid);

Lmin = 0.5;
Lmax = 2.0;
Lfid = 1.0;

Zmin = 0.0;
Zmax = 0.1;
Zfid = 0.03;

Amin = sqrt(smin * fmin) * Lmin / ( sqrt(eta_max) * alpha_max * Mb_max * Tsmax );
Amax = sqrt(smax * fmax) * Lmax / ( sqrt(eta_min) * alpha_min * Mb_min * Tsmin );
Afid = sqrt(sfid * ffid) * Lfid / ( sqrt(eta_fid) * alpha_fid * Mb_fid * Tsfid );

Mv = 10.^[11:0.01:15];
M12 = Mv./1e12;
z = 1:0.01:6;
z3 = (1+z)./3;

delta_min = zeros(length(Mv),length(z));
delta_max = zeros(length(Mv),length(z));
delta_fid = zeros(length(Mv),length(z));
R_min     = zeros(length(Mv),length(z));
R_max     = zeros(length(Mv),length(z));
R_fid     = zeros(length(Mv),length(z));

M_stream  = zeros(3,length(z));

size(M_stream)
size(z)

for i=1:length(z)
    delta_min(:,i) = 100    .* z3(i)       .* (M12(:).^(2/3))                  .* Tmin;
    delta_max(:,i) = 100    .* z3(i)       .* (M12(:).^(2/3))                  .* Tmax;
    delta_fid(:,i) = 100    .* z3(i)       .* (M12(:).^(2/3))                  .* Tfid;
    
    R_min(:,i)     = 18.0   .* (z3(i).^2.5) .* (M12(:).^(1/3)) .* (delta_max(:,i)./100).^(-1.0) .* Amin;
    R_max(:,i)     = 18.0   .* (z3(i).^2.5) .* (M12(:).^(1/3)) .* (delta_min(:,i)./100).^(-1.0) .* Amax;
    R_fid(:,i)     = 18.0   .* (z3(i).^2.5) .* (M12(:).^(1/3)) .* (delta_fid(:,i)./100).^(-1.0) .* Afid;
    for j=1:length(Mv)
        Mtot_min = Mb_min .* sqrt(delta_min(j,i)) ./ ( 1 + sqrt(delta_min(j,i)) );
        Mtot_max = Mb_max .* sqrt(delta_max(j,i)) ./ ( 1 + sqrt(delta_max(j,i)) );
        Mtot_fid = Mb_fid .* sqrt(delta_fid(j,i)) ./ ( 1 + sqrt(delta_fid(j,i)) );
        
        alpha_min2 = 0.21 * ( 0.8*exp(-3*(Mtot_max^2)) + 0.2 ) / 0.1;
        alpha_max2 = 0.21 * ( 0.8*exp(-3*(Mtot_min^2)) + 0.2 ) / 0.1;
        alpha_fid2 = 0.21 * ( 0.8*exp(-3*(Mtot_fid^2)) + 0.2 ) / 0.1;
        
        Amin2 = Amin * alpha_max / alpha_max2;
        Amax2 = Amax * alpha_min / alpha_min2;
        Afid2 = Afid * alpha_fid / alpha_fid2;
        
        R_min(j,i) = (R_min(j,i) / Amin) .* Amin2;
        R_max(j,i) = (R_max(j,i) / Amax)  .* Amax2;
        R_fid(j,i) = (R_fid(j,i) / Afid)  .* Afid2;
    end
    
%     for j=1:length(Mv)
%         Tmix_min = Tsmin .* 1.5e4 .* sqrt( delta_min(j,i) );
%         Tmix_max = Tsmax .* 1.5e4 .* sqrt( delta_max(j,i) );
%         Tmix_fid = Tsfid .* 1.5e4 .* sqrt( delta_fid(j,i) );
%         Tvec = linspace(Tmix_min, Tmix_max, 20);
%         
%         nmix_min = 5.1e-5 .* (z3(i)^3) .* fmin .* sqrt( delta_min(j,i) );
%         nmix_max = 5.1e-5 .* (z3(i)^3) .* fmax .* sqrt( delta_max(j,i) );
%         nmix_fid = 5.1e-5 .* (z3(i)^3) .* ffid .* sqrt( delta_fid(j,i) );
%         nvec = linspace(nmix_min, nmix_max, 20);
        
%         Lmin2 = zeros(length(Tvec), length(nvec));
%         Lmax2 = zeros(length(Tvec), length(nvec));
%         for k=1:length(Tvec)
%             for l=1:length(nvec)
%                 [nspec, mu, cool, heat, cool_tot, heat_tot] = cooling_equilibrium(Tvec(k), nvec(l), ...
%                     z(i), Zfid, 1, 1);
%                 Lmin2(k,l) = abs((cool_tot-heat_tot) .* (10^(22.5)));
%                 
%                 [nspec, mu, cool, heat, cool_tot, heat_tot] = cooling_equilibrium(Tvec(k), nvec(l), ...
%                     z(i), Zfid, 1, 1);
%                 Lmax2(k,l) = abs((cool_tot-heat_tot) .* (10^(22.5)));
%             end
%         end
%         Amin2 = (Amin ./ Lmin) .* min(min(Lmin2));
%         Amax2 = (Amax ./ Lmax) .* max(max(Lmax2));
        
%         [nspec, mu, cool, heat, cool_tot, heat_tot] = cooling_equilibrium(Tmix_fid, nmix_fid, ...
%             z(i), Zmax, 1, 1);
%         Afid2 = (Afid ./ Lfid) .* (cool_tot-heat_tot)./(10.^(-22.5));
    
%         R_min(j,i)     = 18.0   .* (z3(i).^2.5) .* (M12(j).^(1/3)) .* (delta_max(j,i)./100).^(-1.0) .* Amin2;
%         R_max(j,i)     = 18.0   .* (z3(i).^2.5) .* (M12(j).^(1/3)) .* (delta_min(j,i)./100).^(-1.0) .* Amax2;
%         R_fid(j,i)     = 18.0   .* (z3(i).^2.5) .* (M12(j).^(1/3)) .* (delta_fid(j,i)./100).^(-1.0) .* Afid2;
%     end
    
    b = find( abs(R_fid(:,i)-15) == min( abs(R_fid(:,i)-15) ) );
    M_stream(1,i) = Mv(b(1));
    b = find( abs(R_fid(:,i)-20) == min( abs(R_fid(:,i)-20) ) );
    M_stream(2,i) = Mv(b(1));
    b = find( abs(R_fid(:,i)-25) == min( abs(R_fid(:,i)-25) ) );
    M_stream(3,i) = Mv(b(1));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ytit = '$M_{\rm stream}\:{\rm [M_{\odot}] }$';
ytick = [1e11, 1e12, 1e13, 1e14, 1e15];
ystr = {'10^{11}', '10^{12}', '10^{13}', '10^{14}', '10^{15}'};
yu = ytick(end);
yl = ytick(1);

xtit = '$z$';
xtick = 0:1:6;
xu = xtick(end);
xl = xtick(1);

figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
set(figure1,'WindowStyle','docked');
% set(gcf,'renderer','painters')
set(gcf,'renderer','zbuffer')

set(gcf,'visible','off');
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 10 9.5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes1 = axes('Parent',figure1,...
    'YTick',ytick,...
    'YtickLabel',ystr,...
    'YMinorTick','on',...
    'YScale','log',...
    'XTick',xtick,...
    'XMinorTick','on',...
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
    'position',[0.5 -0.09 0]);
ylabel(ytit,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[-0.13 0.5 0]);
% title(tit,'Interpreter','latex','FontSize',16,...
%     'FontName','Times New Roman','units','normalized',...
%     'position',[0.25 0.87 0]);

l(1) = plot(z,M_stream(1,:),'marker','none','linestyle','-','linewidth',2,'color','k',...
    'DisplayName',strcat('$R_{\rm stream}=15R_{\rm crit}$'));
l(2) = plot(z,M_stream(2,:),'marker','none','linestyle','-','linewidth',2,'color','b',...
    'DisplayName',strcat('$R_{\rm stream}=20R_{\rm crit}$'));
l(3) = plot(z,M_stream(3,:),'marker','none','linestyle','-','linewidth',2,'color','r',...
    'DisplayName',strcat('$R_{\rm stream}=25R_{\rm crit}$'));

legend1 = legend(axes1,l(1:3));
set(legend1,...
    'Location','SouthEast','FontSize',10,'Interpreter','latex');
legend boxoff

set(axes1, 'Position',[0.20 0.14 0.74 0.78]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print(gcf,'-djpeg',strcat('./Mstream_20Rcrit.jpg'));
print(gcf,'-depsc',strcat('./Mstream_20Rcrit.eps'));
close all
fclose all;