clear 
LogDel = 0:0.005:4;
Mb = [0:0.01:5];
Del = 10.^LogDel;
nD = length(Del);
nM = length(Mb);
tkh_over_tshear_t = zeros(nD,nM);
tkh_over_tshear_s = zeros(nD,nM);

for j=1:nM
    for i=1:nD
        M = Mb(j);
        D = Del(i);
        Mtot = ( sqrt(D) / (1+sqrt(D)) ) * M;
        alpha = 0.21 * ( 0.8*exp(-3*Mtot^2) + 0.2 );
        tkh = (1+D) ./ sqrt(D);
        tshear_t = 1 ./ alpha;
        tshear_s = (1 + sqrt(D))/alpha;
        tkh_over_tshear_t(i,j) = tkh ./ tshear_t;
        tkh_over_tshear_s(i,j) = tkh ./ tshear_s;
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

xmin = 0;
xmax = 5;
dx = 0.5;
xtick = xmin:dx:xmax;
ymin = 0;
ymax = 4;
dy = 0.5;
ytick = ymin:dy:ymax;
cmin1 = -2.5;
cmax1 = 2.5;
dc1 = 0.5;
cmin2 = -2.5;
cmax2 = 2.5;
dc2 = 0.5;
titx = '$M_{\rm b}$';
tity = '${\rm log}\:\:(\delta)$';
titc_1 = '${\rm log}\:\:(t_{\rm KH}/t_{\rm shear,tot})$';
titc_2 = '${\rm log}\:\:(t_{\rm KH}/t_{\rm shear,stream})$';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68],'Colormap',cmap);
set(figure1,'WindowStyle','docked');

axes1 = axes('Parent',figure1,...
    'YTick',ytick,...
    'YMinorTick','on',...
    'XTick',xtick,...
    'XMinorTick','on',...
    'TickLength',[0.02 0.04],...
    'PlotBoxAspectRatio',[1 1 1],...
    'FontSize',18,...
    'FontName','Times',...
    'CLim',[cmin1 cmax1],...
    'Position',[0.13 0.14 0.775 0.815]);

cv = [cmin1:dc1:cmax1];
bar = colorbar('peer',axes1,'FontSize',18,'FontName','Times','YTick',cv);
set(get(bar,'Ylabel'),'String',titc_1,...
    'Interpreter','latex','Fontsize',22,...
    'Rotation',270,'units','normalized','position',[4.5,0.5,9.16],'FontName','Times');
        
xlim(axes1,[xmin xmax]);
ylim(axes1,[ymin ymax]);
ylabel(tity,'Interpreter','latex','FontSize',22,...
    'FontName','Times New Roman','units','normalized','position',[-0.12 0.5 0]);
xlabel(titx,'Interpreter','latex','FontSize',22,...
    'FontName','Times New Roman','units','normalized','position',[0.5 -0.10 0]);

box(axes1,'on');
view(2);
hold(axes1,'all');
set(gcf,'renderer','Zbuffer')
surf(Mb,LogDel,log10(tkh_over_tshear_t),'linestyle','none');

set(gcf,'visible','off');
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 18 15]);
set(gca, 'Position',[0.10 0.15 0.7500 0.786]);

figdir = './timescale_comparison/';
mkdir(figdir);
figname1 = strcat(figdir,'tkh_over_tshear_total.jpg');

print(gcf,'-djpeg',figname1);
close all
fclose all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68],'Colormap',cmap);
set(figure1,'WindowStyle','docked');

axes1 = axes('Parent',figure1,...
    'YTick',ytick,...
    'YMinorTick','on',...
    'XTick',xtick,...
    'XMinorTick','on',...
    'TickLength',[0.02 0.04],...
    'PlotBoxAspectRatio',[1 1 1],...
    'FontSize',18,...
    'FontName','Times',...
    'CLim',[cmin2 cmax2],...
    'Position',[0.13 0.14 0.775 0.815]);

cv = [cmin2:dc2:cmax2];
bar = colorbar('peer',axes1,'FontSize',18,'FontName','Times','YTick',cv);
set(get(bar,'Ylabel'),'String',titc_2,...
    'Interpreter','latex','Fontsize',22,...
    'Rotation',270,'units','normalized','position',[4.5,0.5,9.16],'FontName','Times');
        
xlim(axes1,[xmin xmax]);
ylim(axes1,[ymin ymax]);
ylabel(tity,'Interpreter','latex','FontSize',22,...
    'FontName','Times New Roman','units','normalized','position',[-0.12 0.5 0]);
xlabel(titx,'Interpreter','latex','FontSize',22,...
    'FontName','Times New Roman','units','normalized','position',[0.5 -0.10 0]);

box(axes1,'on');
view(2);
hold(axes1,'all');
set(gcf,'renderer','Zbuffer')
surf(Mb,LogDel,log10(tkh_over_tshear_s),'linestyle','none');

set(gcf,'visible','off');
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 18 15]);
set(gca, 'Position',[0.10 0.15 0.7500 0.786]);

figdir = './timescale_comparison/';
mkdir(figdir);
figname1 = strcat(figdir,'tkh_over_tshear_stream.jpg');

print(gcf,'-djpeg',figname1);
close all
fclose all;