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

Mv = 10.^[11:0.05:14];
M12 = Mv./1e12;
z  = [0, 2, 4];
z3 = (1+z)./3;

delta_min = zeros(length(Mv),length(z));
delta_max = zeros(length(Mv),length(z));
delta_fid = zeros(length(Mv),length(z));
R_min     = zeros(length(Mv),length(z));
R_max     = zeros(length(Mv),length(z));
R_fid     = zeros(length(Mv),length(z));

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
%     
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
%         
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
%         
%         [nspec, mu, cool, heat, cool_tot, heat_tot] = cooling_equilibrium(Tmix_fid, nmix_fid, ...
%             z(i), Zfid, 1, 1);
%         Afid2 = (Afid ./ Lfid) .* (cool_tot-heat_tot)./(10.^(-22.5));
%     
%         R_min(j,i)     = 18.0   .* (z3(i).^2.5) .* (M12(j).^(1/3)) .* (delta_max(j,i)./100).^(-1.0) .* Amin2;
%         R_max(j,i)     = 18.0   .* (z3(i).^2.5) .* (M12(j).^(1/3)) .* (delta_min(j,i)./100).^(-1.0) .* Amax2;
%         R_fid(j,i)     = 18.0   .* (z3(i).^2.5) .* (M12(j).^(1/3)) .* (delta_fid(j,i)./100).^(-1.0) .* Afid2;
%     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xtit = '$M_{\rm v}\:{\rm [M_{\odot}] }$';
xtick = [1e11, 1e12, 1e13, 1e14];
xstr = {'10^{11}', '10^{12}', '10^{13}', '10^{14}'};
xu = xtick(end);
xl = xtick(1);

ytit = '$R_{\rm s}/R_{\rm s,crit}$';
ytick = [0.01, 0.1, 1, 10, 100, 1000, 1e4];
ystr = {'0.01', '0.1', '1', '10', '100', '1000', '10000'};
yu = ytick(end);
yl = ytick(1);

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
    'position',[0.5 -0.09 0]);
ylabel(ytit,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[-0.13 0.5 0]);
% title(tit,'Interpreter','latex','FontSize',16,...
%     'FontName','Times New Roman','units','normalized',...
%     'position',[0.25 0.87 0]);

% xshade1 = linspace(1.05.*Mv(1),0.95.*Mv(end),100);
xshade1 = Mv;
xshade = [xshade1, fliplr(xshade1)];
yshade1 = 1    .* ones(size(xshade1));
yshade2 = 1e-3 .* ones(size(xshade1));
yshade = [yshade1, fliplr(yshade2)];
fill(xshade,yshade,[0.7, 0.7, 0.7],'facealpha',0.7);

set(axes1,'YTick',ytick,'YtickLabel',ystr,'YMinorTick','on','YScale','log','TickLength',[0.02 0.04])
n = floor(length(Mv)/2);

l(1) = plot(Mv,R_fid(:,1),'marker','none','linestyle','-','linewidth',2,'color','g',...
    'DisplayName',strcat('$z=',num2str(z(1),'%4.2f'),'$'));
plot(Mv,R_min(:,1),'marker','none','linestyle',':','linewidth',1,'color','g');
plot(Mv,R_max(:,1),'marker','none','linestyle',':','linewidth',1,'color','g');
y = linspace(R_min(n-2,1),R_max(n-2,1),10);
x = Mv(n-2).*ones(size(y));
plot(x,y,'marker','none','linestyle','-','linewidth',1.5,'color','g');
plot(Mv(n-3:n-1),R_min(n-3:n-1,1),'marker','none','linestyle','-','linewidth',1.5,'color','g');
plot(Mv(n-3:n-1),R_max(n-3:n-1,1),'marker','none','linestyle','-','linewidth',1.5,'color','g');

l(2) = plot(Mv,R_fid(:,2),'marker','none','linestyle','-','linewidth',2,'color','b',...
    'DisplayName',strcat('$z=',num2str(z(2),'%4.2f'),'$'));
plot(Mv,R_min(:,2),'marker','none','linestyle',':','linewidth',1,'color','b');
plot(Mv,R_max(:,2),'marker','none','linestyle',':','linewidth',1,'color','b');
y = linspace(R_min(n,2),R_max(n,2),10);
x = Mv(n).*ones(size(y));
plot(x,y,'marker','none','linestyle','-','linewidth',1.5,'color','b');
plot(Mv(n-1:n+1),R_min(n-1:n+1,2),'marker','none','linestyle','-','linewidth',1.5,'color','b');
plot(Mv(n-1:n+1),R_max(n-1:n+1,2),'marker','none','linestyle','-','linewidth',1.5,'color','b');

l(3) = plot(Mv,R_fid(:,3),'marker','none','linestyle','-','linewidth',2,'color','r',...
    'DisplayName',strcat('$z=',num2str(z(3),'%4.2f'),'$'));
plot(Mv,R_min(:,3),'marker','none','linestyle',':','linewidth',1,'color','r');
plot(Mv,R_max(:,3),'marker','none','linestyle',':','linewidth',1,'color','r');
y = linspace(R_min(n+2,3),R_max(n+2,3),10);
x = Mv(n+2).*ones(size(y));
plot(x,y,'marker','none','linestyle','-','linewidth',1.5,'color','r');
plot(Mv(n+1:n+3),R_min(n+1:n+3,3),'marker','none','linestyle','-','linewidth',1.5,'color','r');
plot(Mv(n+1:n+3),R_max(n+1:n+3,3),'marker','none','linestyle','-','linewidth',1.5,'color','r');

% for i=1:(length(Mv)-2)
%     y = linspace(R_min(i+1,1),R_max(i+1,1),10);
%     x = Mv(i+1).*ones(size(y));
%     plot(x,y,'marker','none','linestyle','-','linewidth',0.5,'color','g');
%     
%     x = log10(linspace(Mv(i),Mv(i+2),10));
%     y = log10(R_min(i,2)) + (log10(R_max(i+2,2))-log10(R_min(i,2)))./...
%         (log10(Mv(i+2))-log10(Mv(i))).*(x-log10(Mv(i)));
%     plot(10.^x,10.^y,'marker','none','linestyle','-','linewidth',0.5,'color','b');
%     
%     x = log10(linspace(Mv(i),Mv(i+2),10));
%     y = log10(R_max(i,3)) + (log10(R_min(i+2,3))-log10(R_max(i,3)))./...
%         (log10(Mv(i+2))-log10(Mv(i))).*(x-log10(Mv(i)));
%     plot(10.^x,10.^y,'marker','none','linestyle','-','linewidth',0.5,'color','r');
% end

annotation(figure1,'textbox',[0.5 0.181 0.10 0.05],...
    'String',strcat('${\rm Stream\:Disruption}$'),...
    'FontSize',10,...
    'FontName','Times',...
    'Interpreter','Latex',...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor','none');

legend3 = legend(axes1,l(1:3));
set(legend3,...
    'Location','SouthWest','FontSize',10,'Interpreter','latex');
legend boxoff

set(axes1, 'Position',[0.20 0.14 0.74 0.78]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print(gcf,'-djpeg',strcat('./Rscrit.jpg'));
print(gcf,'-depsc',strcat('./Rscrit.eps'));
% print(gcf,'-djpeg',strcat(fig_direc,'net_lumis_',num2str(snap_num,'%03i'),'.jpg'));
%print(gcf,'-depsc',strcat(fig_direc,'net_lumis_',num2str(snap_num,'%03i'),'.eps'));
close all
fclose all;