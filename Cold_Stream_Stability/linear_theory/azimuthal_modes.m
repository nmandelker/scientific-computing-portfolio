m = 5;
t = 0:0.0001:2*pi;
R0 = 0.5;
R = R0 + 0.1.*R0.*cos(m.*t);
x = R.*cos(t);
y = R.*sin(t);
x0 = R0.*cos(t);
y0 = R0.*sin(t);

% Create figure
figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
set(figure1,'WindowStyle','docked');

xt = (-1.2*R0):0.1:(1.2*R0);
yt = (-1.2*R0):0.1:(1.2*R0);
% Create axes
axes1 = axes('Parent',figure1,...
    'YTick',yt,...
    'XTick',xt,...
    'Position',[0.12 0.15 0.8 0.8],...
    'PlotBoxAspectRatio',[1 1 1],...
    'FontSize',16,...
    'FontName','Arial');

xlim([xt(1) xt(length(xt))]);
ylim([yt(1) yt(length(yt))]);
grid('off');
hold('all');

titx = '$x$';
tity = '$y$';
tit = strcat('$m = ',num2str(m,'%1i'),'$');

xlabel(titx,'Interpreter','latex','FontSize',24,...
    'units','normalized','position',[0.5 -0.07 0],...
    'FontName','Times New Roman');

ylabel(tity,'Interpreter','latex','FontSize',24,...
    'units','normalized','position',[-0.07 0.5 0],...
    'FontName','Times New Roman');

title(tit,'Interpreter','latex','FontSize',30,...
    'units','normalized','position',[0.5 0.5 0],...
    'FontName','Times New Roman');%,'BackGroundColor','w');

plot(x0,y0,'linestyle','--','color','k','marker','none','linewidth',3)
plot(x,y,'linestyle','-','color','k','marker','none','linewidth',3)

set(gcf,'renderer','painters')