clear

Thmin = 3/8;
Thmax = 1;

Tsmin = 0.5;
Tsmax = 2.0;

% Tmin = Thmin/Tsmax;
% Tmax = Thmax/Tsmin;
Tmin = 0.3;
Tmax = 3.0;
Tfid = 1.0;

% eta_min = 0.5;
% eta_max = sqrt(2);
eta_fid = 1.0;
eta_min = 1.0;
eta_max = 1.0;

fmin = 1.0;
fmax = 3.0;
ffid = 1.0;

smin = 0.3;
smax = 3.0;
sfid = 1.0;

Amin = sqrt(smin / (eta_max*fmax));
Amax = sqrt(smax / (eta_min*fmin));
Afid = sqrt(sfid / (eta_fid*ffid));

Mv = 10.^[11:0.05:14];
M12 = Mv./1e12;
z  = [1, 2, 4];
z3 = (1+z)./3;

delta_min = zeros(length(Mv),length(z));
delta_max = zeros(length(Mv),length(z));
delta_fid = zeros(length(Mv),length(z));
n_min     = zeros(length(Mv),length(z));
n_max     = zeros(length(Mv),length(z));
n_fid     = zeros(length(Mv),length(z));
R_min     = zeros(length(Mv),length(z));
R_max     = zeros(length(Mv),length(z));
R_fid     = zeros(length(Mv),length(z));

for i=1:length(z)
    delta_min(:,i) = 100    .* z3(i)       .* (M12.^(2/3))                  .* Tmin;
    delta_max(:,i) = 100    .* z3(i)       .* (M12.^(2/3))                  .* Tmax;
    delta_fid(:,i) = 100    .* z3(i)       .* (M12.^(2/3))                  .* Tfid;
    
    n_min(:,i)     = 5.1e-3 .* (z3(i)^3)   .* (delta_min(:,i)./100)         .* fmin;
    n_max(:,i)     = 5.1e-3 .* (z3(i)^3)   .* (delta_max(:,i)./100)         .* fmax;
    n_fid(:,i)     = 5.1e-3 .* (z3(i)^3)   .* (delta_fid(:,i)./100)         .* ffid;
    
    R_min(:,i)     = 0.16   .* (z3(i)^0.5) .* (delta_max(:,i)./100).^(-0.5) .* Amin;
    R_max(:,i)     = 0.16   .* (z3(i)^0.5) .* (delta_min(:,i)./100).^(-0.5) .* Amax;
    R_fid(:,i)     = 0.16   .* (z3(i)^0.5) .* (delta_fid(:,i)./100).^(-0.5) .* Afid;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xtit = '$M_{\rm v}\:{\rm [M_{\odot}] }$';
xtick = [1e11, 1e12, 1e13, 1e14];
xstr = {'10^{11}', '10^{12}', '10^{13}', '10^{14}'};
xu = xtick(end);
xl = xtick(1);

ytit1 = '$\delta$';
ytick1 = [1, 10, 100, 1000, 10000];
ystr1 = {'1', '10', '100', '1000', '10000'};
yu1 = ytick1(end);
yl1 = ytick1(1);

ytit2 = '$n_{\rm H,0}\:{\rm [cm^{-3}] }$';
ytick2 = [1e-4, 1e-3, 1e-2, 0.1, 1];
ystr2 = {'10^{-4}', '10^{-3}', '0.01', '0.1', '1'};
yu2 = ytick2(end);
yl2 = ytick2(1);

ytit3 = '$R_{\rm s}/R_{\rm v}$';
ytick3 = [0.01, 0.1, 1];
ystr3 = {'0.01', '0.1', '1'};
yu3 = ytick3(end);
yl3 = ytick3(1);

figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
set(figure1,'WindowStyle','docked');
set(gcf,'renderer','painters')

set(gcf,'visible','off');
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 30 9.5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes1 = axes('Parent',figure1,...
    'YTick',ytick1,...
    'YtickLabel',ystr1,...
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
ylim(axes1,[yl1 yu1]);
box(axes1,'on');
% grid on
hold(axes1,'all');
xlabel(xtit,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 -0.09 0]);
ylabel(ytit1,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[-0.13 0.5 0]);
% title(tit,'Interpreter','latex','FontSize',16,...
%     'FontName','Times New Roman','units','normalized',...
%     'position',[0.25 0.87 0]);

n = floor(length(Mv)/2);

l(1) = plot(Mv,delta_fid(:,1),'marker','none','linestyle','-','linewidth',2,'color','g',...
    'DisplayName',strcat('$z=',num2str(z(1),'%4.2f'),'$'));
plot(Mv,delta_min(:,1),'marker','none','linestyle',':','linewidth',1,'color','g');
plot(Mv,delta_max(:,1),'marker','none','linestyle',':','linewidth',1,'color','g');
y = linspace(delta_min(n-2,1),delta_max(n-2,1),10);
x = Mv(n-2).*ones(size(y));
plot(x,y,'marker','none','linestyle','-','linewidth',1.5,'color','g');
plot(Mv(n-3:n-1),delta_min(n-3:n-1,1),'marker','none','linestyle','-','linewidth',1.5,'color','g');
plot(Mv(n-3:n-1),delta_max(n-3:n-1,1),'marker','none','linestyle','-','linewidth',1.5,'color','g');

l(2) = plot(Mv,delta_fid(:,2),'marker','none','linestyle','-','linewidth',2,'color','b',...
    'DisplayName',strcat('$z=',num2str(z(2),'%4.2f'),'$'));
plot(Mv,delta_min(:,2),'marker','none','linestyle',':','linewidth',1,'color','b');
plot(Mv,delta_max(:,2),'marker','none','linestyle',':','linewidth',1,'color','b');
y = linspace(delta_min(n,2),delta_max(n,2),10);
x = Mv(n).*ones(size(y));
plot(x,y,'marker','none','linestyle','-','linewidth',1.5,'color','b');
plot(Mv(n-1:n+1),delta_min(n-1:n+1,2),'marker','none','linestyle','-','linewidth',1.5,'color','b');
plot(Mv(n-1:n+1),delta_max(n-1:n+1,2),'marker','none','linestyle','-','linewidth',1.5,'color','b');

l(3) = plot(Mv,delta_fid(:,3),'marker','none','linestyle','-','linewidth',2,'color','r',...
    'DisplayName',strcat('$z=',num2str(z(3),'%4.2f'),'$'));
plot(Mv,delta_min(:,3),'marker','none','linestyle',':','linewidth',1,'color','r');
plot(Mv,delta_max(:,3),'marker','none','linestyle',':','linewidth',1,'color','r');
y = linspace(delta_min(n+2,3),delta_max(n+2,3),10);
x = Mv(n+2).*ones(size(y));
plot(x,y,'marker','none','linestyle','-','linewidth',1.5,'color','r');
plot(Mv(n+1:n+3),delta_min(n+1:n+3,3),'marker','none','linestyle','-','linewidth',1.5,'color','r');
plot(Mv(n+1:n+3),delta_max(n+1:n+3,3),'marker','none','linestyle','-','linewidth',1.5,'color','r');

% for i=1:(length(Mv)-2)
%     y = linspace(delta_min(i+1,1),delta_max(i+1,1),10);
%     x = Mv(i+1).*ones(size(y));
%     plot(x,y,'marker','none','linestyle','-','linewidth',0.5,'color','g');
%     
%     x = log10(linspace(Mv(i),Mv(i+2),10));
%     y = log10(delta_min(i,2)) + (log10(delta_max(i+2,2))-log10(delta_min(i,2)))./...
%         (log10(Mv(i+2))-log10(Mv(i))).*(x-log10(Mv(i)));
%     plot(10.^x,10.^y,'marker','none','linestyle','-','linewidth',0.5,'color','b');
%     
%     x = log10(linspace(Mv(i),Mv(i+2),10));
%     y = log10(delta_max(i,3)) + (log10(delta_min(i+2,3))-log10(delta_max(i,3)))./...
%         (log10(Mv(i+2))-log10(Mv(i))).*(x-log10(Mv(i)));
%     plot(10.^x,10.^y,'marker','none','linestyle','-','linewidth',0.5,'color','r');
% end

legend1 = legend(axes1,l(1:3));
set(legend1,...
    'Location','SouthEast','FontSize',10,'Interpreter','latex');
legend boxoff

set(axes1, 'Position',[0.065 0.14 0.2467 0.78]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes2 = axes('Parent',figure1,...
    'YTick',ytick2,...
    'YtickLabel',ystr2,...
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
    'position',[-0.13 0.5 0]);
% title(tit,'Interpreter','latex','FontSize',16,...
%     'FontName','Times New Roman','units','normalized',...
%     'position',[0.25 0.87 0]);

l(1) = plot(Mv,n_fid(:,1),'marker','none','linestyle','-','linewidth',2,'color','g',...
    'DisplayName',strcat('$z=',num2str(z(1),'%4.2f'),'$'));
plot(Mv,n_min(:,1),'marker','none','linestyle',':','linewidth',1,'color','g');
plot(Mv,n_max(:,1),'marker','none','linestyle',':','linewidth',1,'color','g');
y = linspace(n_min(n-2,1),n_max(n-2,1),10);
x = Mv(n-2).*ones(size(y));
plot(x,y,'marker','none','linestyle','-','linewidth',1.5,'color','g');
plot(Mv(n-3:n-1),n_min(n-3:n-1,1),'marker','none','linestyle','-','linewidth',1.5,'color','g');
plot(Mv(n-3:n-1),n_max(n-3:n-1,1),'marker','none','linestyle','-','linewidth',1.5,'color','g');

l(2) = plot(Mv,n_fid(:,2),'marker','none','linestyle','-','linewidth',2,'color','b',...
    'DisplayName',strcat('$z=',num2str(z(2),'%4.2f'),'$'));
plot(Mv,n_min(:,2),'marker','none','linestyle',':','linewidth',1,'color','b');
plot(Mv,n_max(:,2),'marker','none','linestyle',':','linewidth',1,'color','b');
y = linspace(n_min(n,2),n_max(n,2),10);
x = Mv(n).*ones(size(y));
plot(x,y,'marker','none','linestyle','-','linewidth',1.5,'color','b');
plot(Mv(n-1:n+1),n_min(n-1:n+1,2),'marker','none','linestyle','-','linewidth',1.5,'color','b');
plot(Mv(n-1:n+1),n_max(n-1:n+1,2),'marker','none','linestyle','-','linewidth',1.5,'color','b');

l(3) = plot(Mv,n_fid(:,3),'marker','none','linestyle','-','linewidth',2,'color','r',...
    'DisplayName',strcat('$z=',num2str(z(3),'%4.2f'),'$'));
plot(Mv,n_min(:,3),'marker','none','linestyle','-','linewidth',1,'color','r');
plot(Mv,n_max(:,3),'marker','none','linestyle','-','linewidth',1,'color','r');
y = linspace(n_min(n+2,3),n_max(n+2,3),10);
x = Mv(n+2).*ones(size(y));
plot(x,y,'marker','none','linestyle','-','linewidth',1.5,'color','r');
plot(Mv(n+1:n+3),n_min(n+1:n+3,3),'marker','none','linestyle','-','linewidth',1.5,'color','r');
plot(Mv(n+1:n+3),n_max(n+1:n+3,3),'marker','none','linestyle','-','linewidth',1.5,'color','r');

% for i=1:(length(Mv)-2)
%     y = linspace(n_min(i+1,1),n_max(i+1,1),10);
%     x = Mv(i+1).*ones(size(y));
%     plot(x,y,'marker','none','linestyle','-','linewidth',0.5,'color','g');
%     
%     x = log10(linspace(Mv(i),Mv(i+2),10));
%     y = log10(n_min(i,2)) + (log10(n_max(i+2,2))-log10(n_min(i,2)))./...
%         (log10(Mv(i+2))-log10(Mv(i))).*(x-log10(Mv(i)));
%     plot(10.^x,10.^y,'marker','none','linestyle','-','linewidth',0.5,'color','b');
%     
%     x = log10(linspace(Mv(i),Mv(i+2),10));
%     y = log10(n_max(i,3)) + (log10(n_min(i+2,3))-log10(n_max(i,3)))./...
%         (log10(Mv(i+2))-log10(Mv(i))).*(x-log10(Mv(i)));
%     plot(10.^x,10.^y,'marker','none','linestyle','-','linewidth',0.5,'color','r');
% end

legend2 = legend(axes2,l(1:3));
set(legend2,...
    'Location','SouthEast','FontSize',10,'Interpreter','latex');
legend boxoff

set(axes2, 'Position',[0.375 0.14 0.2467 0.78]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes3 = axes('Parent',figure1,...
    'YTick',ytick3,...
    'YtickLabel',ystr3,...
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
% title(tit,'Interpreter','latex','FontSize',16,...
%     'FontName','Times New Roman','units','normalized',...
%     'position',[0.25 0.87 0]);

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
plot(Mv,R_min(:,3),'marker','none','linestyle','-','linewidth',1,'color','r');
plot(Mv,R_max(:,3),'marker','none','linestyle','-','linewidth',1,'color','r');
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

legend3 = legend(axes3,l(1:3));
set(legend3,...
    'Location','SouthWest','FontSize',10,'Interpreter','latex');
legend boxoff

set(axes3, 'Position',[0.685 0.14 0.2467 0.78]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print(gcf,'-djpeg',strcat('./params.jpg'));
print(gcf,'-depsc',strcat('./params.eps'));
% print(gcf,'-djpeg',strcat(fig_direc,'net_lumis_',num2str(snap_num,'%03i'),'.jpg'));
%print(gcf,'-depsc',strcat(fig_direc,'net_lumis_',num2str(snap_num,'%03i'),'.eps'));
close all
fclose all;