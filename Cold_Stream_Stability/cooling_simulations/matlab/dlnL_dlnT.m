clear

T=10.^([3:0.01:9]);
n = 10^(-3.3);
red = 2.99;
Z = 1e-3;

NT = length(T);
cool_net = zeros(3,NT);
dlnLdlnT = zeros(3,NT);

for i=1:NT
    [nspec, mu, cool, heat, cool_tot, heat_tot] = cooling_equilibrium(T(i), n, red, Z, 0, 0);
    cool_net(1,i) = cool_tot;
    [nspec, mu, cool, heat, cool_tot, heat_tot] = cooling_equilibrium(T(i), n, red, Z, 1, 1);
    cool_net(2,i) = cool_tot;
    cool_net(3,i) = (cool_tot - heat_tot);
end
for i=1:NT-1
    for j=1:3
        if(cool_net(j,i+1)>0 & cool_net(j,i)>0)
            dlnLdlnT(j,i+1) = log10(cool_net(j,i+1)/cool_net(j,i)) / log10(T(i+1)/T(i));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tit = strcat('${\rm [Z]=-3.0},\:{\rm log}(n/{\rm cm^{-3}})=-3.3,\:z=',num2str(red,'%4.2f'),'$');

xtit  = '${\rm log}(T~{\rm [K]})$';
xtick = 4:0.2:6;
xstr  = {'4.0', '4.2', '4.4', '4.6', '4.8', '5.0', '5.2', '5.4', '5.6', '5.8', '6.0'};
xl    = 4.0;
xu    = 6.0;

ytit  = '$d{\rm log}~\Lambda/d{\rm log~T}$';
ytick = -4:1:4;
ystr  = {'-4.0', '-3.0', '-2.0', '-1.0', '0.0', '1.0', '2.0', '3.0', '4.0'};
yl    = -4.0;
yu    = 4.0;
% ytick = -3:0.5:3;
% ystr  = {'-3.0', '-2.5', '-2.0', '-1.5', '-1.0', '-0.5', '0.0', '0.5', '1.0', '1.5', '2.0', '2.5', '3.0'};
% yl    = -3.0;
% yu    = 3.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
set(figure1,'WindowStyle','docked');
set(gcf,'renderer','painters')

set(gcf,'visible','off');
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 10 9.5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axes1 = axes('Parent',figure1,...
    'YTick',ytick,...
    'YMinorTick','on',...
    'YtickLabel',ystr,...
    'XTick',xtick,...
    'XMinorTick','on',...
    'XtickLabel',xstr,...
    'TickLength',[0.02 0.04],...
    'PlotBoxAspectRatio',[1 1 1],...
    'LineWidth',1.0,...
    'FontSize',12,...
    'FontName','Times New Roman',...
    'Position',[0.13 0.14 0.775 0.815]);
xlim(axes1,[xl xu]);
ylim(axes1,[yl yu]);
axes1.TickLabelInterpreter = 'latex';
box(axes1,'on');
% grid on
hold(axes1,'all');
xlabel(xtit,'Interpreter','latex','FontSize',14,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 -0.08 0]);
ylabel(ytit,'Interpreter','latex','FontSize',14,...
    'FontName','Times New Roman','units','normalized',...
    'position',[-0.12 0.5 0]);
title(tit,'Interpreter','latex','FontSize',12,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 1.02 0]);

% l(1) = plot(log10(T),dlnLdlnT(1,:),'linestyle','-','linewidth',1.5,'marker','none','color','k','DisplayName','${\rm w/o\:PI,\:w/o\:PH}$');
l(1) = plot(log10(T),dlnLdlnT(2,:),'linestyle','-','linewidth',1.5,'marker','none','color','r','DisplayName','${\rm w/\:\:\:PI,\:w/o\:PH}$');
b = find(cool_net(3,:)>0);
l(2) = plot(log10(T(b)),dlnLdlnT(3,b),'linestyle','-','linewidth',1.5,'marker','none','color','b','DisplayName','${\rm w/\:\:\:PI,\:w/\:\:\:PH}$');

x = log10(T);
y = zeros(size(x));
plot(x,y,'linestyle','--','linewidth',1.0,'marker','none','color','k');
y = 2.0.*ones(size(x));
plot(x,y,'linestyle','--','linewidth',1.0,'marker','none','color','k');

legend1 = legend(axes1,l(1:2));
set(legend1,...
    'Location','SouthEast','FontSize',10,'Interpreter','latex');
%       'Position',[0.055, 0.73, 0.133, 0.149],'FontSize',11,'Interpreter','latex');
legend boxoff

set(axes1, 'Position',[0.18 0.13 0.74 0.78]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print(gcf,'-dpng',strcat('dlnL_dlnT_z',num2str(round(red),'%01i'),'.png'));
close all
fclose all;
