clear
zsolar = 0.03;
red = 2;
nH = 0.01;

fid=fopen('./cooling_simulation_output/cooling_00001.out');
ntemp=fread(fid,1,'int')
NRho=fread(fid,1,'int');
NT=fread(fid,1,'int');
ntemp=fread(fid,1,'int')
ntemp=fread(fid,1,'int')
dens=fread(fid,NRho,'float64');
ntemp=fread(fid,1,'int')
ntemp=fread(fid,1,'int')
T2=fread(fid,NT,'float64');
ntemp=fread(fid,1,'int')
ntemp=fread(fid,1,'int')
cool=fread(fid,NRho*NT,'float64');
ntemp=fread(fid,1,'int')
ntemp=fread(fid,1,'int')
heat=fread(fid,NRho*NT,'float64');
ntemp=fread(fid,1,'int')
ntemp=fread(fid,1,'int')
cool_com=fread(fid,NRho*NT,'float64');
ntemp=fread(fid,1,'int')
ntemp=fread(fid,1,'int')
heat_com=fread(fid,NRho*NT,'float64');
ntemp=fread(fid,1,'int')
ntemp=fread(fid,1,'int')
metal=fread(fid,NRho*NT,'float64');
ntemp=fread(fid,1,'int')
ntemp=fread(fid,1,'int')
cool_prime=fread(fid,NRho*NT,'float64');
ntemp=fread(fid,1,'int')
ntemp=fread(fid,1,'int')
heat_prime=fread(fid,NRho*NT,'float64');
ntemp=fread(fid,1,'int')
ntemp=fread(fid,1,'int')
cool_com_prime=fread(fid,NRho*NT,'float64');
ntemp=fread(fid,1,'int')
ntemp=fread(fid,1,'int')
heat_com_prime=fread(fid,NRho*NT,'float64');
ntemp=fread(fid,1,'int')
ntemp=fread(fid,1,'int')
metal_prime=fread(fid,NRho*NT,'float64');
ntemp=fread(fid,1,'int')
ntemp=fread(fid,1,'int')
mu=fread(fid,NRho*NT,'float64');
ntemp=fread(fid,1,'int')
ntemp=fread(fid,1,'int')
nspec=fread(fid,NRho*NT*6,'float64');

cool=reshape(cool,[NRho,NT]);
heat=reshape(heat,[NRho,NT]);
cool_com=reshape(cool_com,[NRho,NT]);
heat_com=reshape(heat_com,[NRho,NT]);
metal=reshape(metal,[NRho,NT]);
cool_prime=reshape(cool_prime,[NRho,NT]);
heat_prime=reshape(heat_prime,[NRho,NT]);
cool_com_prime=reshape(cool_com_prime,[NRho,NT]);
heat_com_prime=reshape(heat_com_prime,[NRho,NT]);
metal_prime=reshape(metal_prime,[NRho,NT]);
mu=reshape(mu,[NRho,NT]);
nspec=reshape(nspec,[NRho,NT,6]);
dens = 10.^dens';
T2 = 10.^T2';

b = find(abs(dens-nH)==min(abs(dens-nH)));
%b=82;
mu = mu(b,:);
T = mu.*T2;

Lcool = 10.^(cool(b,:));
Lheat = 10.^(heat(b,:));
Lmetal = zsolar.*(10.^(metal(b,:)));
Lcool_com = 10.^(cool_com(b,:))./dens(b);
Lheat_com = 10.^(heat_com(b,:))./dens(b);
nspec = 10.^(nspec(b,:,:))./dens(b);

Lcool1 = zeros(size(T));
Lheat1 = zeros(size(T));
Lcool_com1 = zeros(size(T));
Lheat_com1 = zeros(size(T));
Lmetal1 = zeros(size(T));
mu1 = zeros(size(T));
nspec1 = zeros(length(T),6);

for i=1:length(T)
%     [nspec1(i,:), mu1(i), cool1, heat1, cool_tot1, heat_tot1] = cooling_equilibrium(T(i), nH, red, zsolar, 1, 1);
    [nspec1(i,:), mu1(i), cool1, heat1, cool_tot1, heat_tot1] = cooling_equilibrium_T2_input(T2(i), nH, red, zsolar, 1, 1);
    Lcool1(i) = sum(cool1(1:13));
    Lheat1(i) = sum(heat1(1:3));
    Lcool_com1(i) = cool1(14);
    Lheat_com1(i) = heat1(4);
    Lmetal1(i) = sum(cool1(15:16));
end
nspec1 = nspec1./nH;

xtick = [1e3 1e4 1e5 1e6 1e7];
ytick = [1e-28 1e-27 1e-26 1e-25 1e-24 1e-23 1e-22 1e-21];
xstr = {'10^3','10^4','10^5','10^6','10^7'};
ystr = {'10^{-28}','10^{-27}','10^{-26}','10^{-25}','10^{-24}','10^{-23}','10^{-22}','10^{-21}'};
figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
set(figure1,'WindowStyle','docked');
axes1 = axes('Parent',figure1,'YTick',ytick,...
    'YScale','log',...
    'YtickLabel',ystr,...
    'XTick',xtick,...
    'XScale','log',...
    'XtickLabel',xstr,...%,xl:dx:xu,...'XMinorTick','on',...
    'TickLength',[0.02 0.04],...
    'PlotBoxAspectRatio',[1 1 1],...
    'LineWidth',1.5,...
    'FontSize',16,...
    'FontName','Arial',...
    'Position',[0.1076 0.1200 0.8724 0.8524]);%[0.13 0.14 0.775 0.815]);
xlim(axes1,[xtick(1) xtick(length(xtick))]);
ylim(axes1,[ytick(1) ytick(length(ytick))]);
box(axes1,'on');
hold(axes1,'all');
xlabel('$T\:{\rm [K]}$','Interpreter','latex','FontSize',22,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 -0.07 0]);
ylabel('$|\Lambda|\:{\rm [erg\,s^{-1}\,cm^3]}$','Interpreter','latex','FontSize',22,...
    'FontName','Times New Roman','units','normalized',...
    'position',[-0.11 0.5 0]);

plot(T2,Lcool,'linewidth',2,'color','k','linestyle','-');
plot(T2,Lheat,'linewidth',2,'color','r','linestyle','-');
plot(T2,Lmetal,'linewidth',2,'color','b','linestyle','-');
plot(T2,Lcool_com,'linewidth',2,'color','c','linestyle','-');
plot(T2,Lheat_com,'linewidth',2,'color','m','linestyle','-');
% plot(T,abs(Lcool+Lmetal+Lcool_com-Lheat-Lheat_com),'linewidth',2,'color','g','linestyle','-');

plot(T2,Lcool1,'linewidth',2,'color','k','linestyle','--');
plot(T2,Lheat1,'linewidth',2,'color','r','linestyle','--');
plot(T2,Lmetal1,'linewidth',2,'color','b','linestyle','--');
plot(T2,Lcool_com1,'linewidth',2,'color','c','linestyle','--');
plot(T2,Lheat_com1,'linewidth',2,'color','m','linestyle','--');
% plot(T,abs(Lcool1+Lmetal1+Lcool_com1-Lheat1-Lheat_com1),'linewidth',2,'color','g','linestyle','--');


xtick = [1e3 1e4 1e5 1e6 1e7];
ytick = 0.55:0.01:0.65;
xstr = {'10^3','10^4','10^5','10^6','10^7'};
figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
set(figure1,'WindowStyle','docked');
axes1 = axes('Parent',figure1,'YTick',ytick,...
    'YMinorTick','on',...
    'XTick',xtick,...
    'XScale','log',...
    'XtickLabel',xstr,...%,xl:dx:xu,...'XMinorTick','on',...
    'TickLength',[0.02 0.04],...
    'PlotBoxAspectRatio',[1 1 1],...
    'LineWidth',1.5,...
    'FontSize',16,...
    'FontName','Arial',...
    'Position',[0.1076 0.1200 0.8724 0.8524]);%[0.13 0.14 0.775 0.815]);
xlim(axes1,[xtick(1) xtick(length(xtick))]);
ylim(axes1,[ytick(1) ytick(length(ytick))]);
box(axes1,'on');
hold(axes1,'all');
xlabel('$T\:{\rm [K]}$','Interpreter','latex','FontSize',22,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 -0.07 0]);
ylabel('$\mu$','Interpreter','latex','FontSize',22,...
    'FontName','Times New Roman','units','normalized',...
    'position',[-0.11 0.5 0]);

plot(T2,mu,'linewidth',2,'color','k','linestyle','-');
plot(T2,mu1,'linewidth',2,'color','r','linestyle','--');


xtick = [1e3 1e4 1e5 1e6 1e7];
ytick = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1];
xstr = {'10^3','10^4','10^5','10^6','10^7'};
ystr = {'10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}','10^{1}'};
figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
set(figure1,'WindowStyle','docked');
axes1 = axes('Parent',figure1,'YTick',ytick,...
    'YScale','log',...
    'YtickLabel',ystr,...
    'XTick',xtick,...
    'XScale','log',...
    'XtickLabel',xstr,...%,xl:dx:xu,...'XMinorTick','on',...
    'TickLength',[0.02 0.04],...
    'PlotBoxAspectRatio',[1 1 1],...
    'LineWidth',1.5,...
    'FontSize',16,...
    'FontName','Arial',...
    'Position',[0.1076 0.1200 0.8724 0.8524]);%[0.13 0.14 0.775 0.815]);
xlim(axes1,[xtick(1) xtick(length(xtick))]);
ylim(axes1,[ytick(1) ytick(length(ytick))]);
box(axes1,'on');
hold(axes1,'all');
xlabel('$T\:{\rm [K]}$','Interpreter','latex','FontSize',22,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 -0.07 0]);
ylabel('$n_{\rm spec}\:/\:n_{\rm H}$','Interpreter','latex','FontSize',22,...
    'FontName','Times New Roman','units','normalized',...
    'position',[-0.11 0.5 0]);

plot(T2,nspec(1,:,1),'linewidth',2,'color','k','linestyle','-');
plot(T2,nspec(1,:,2),'linewidth',2,'color','r','linestyle','-');
plot(T2,nspec(1,:,3),'linewidth',2,'color','b','linestyle','-');
plot(T2,nspec(1,:,4),'linewidth',2,'color','c','linestyle','-');
plot(T2,nspec(1,:,5),'linewidth',2,'color','m','linestyle','-');
plot(T2,nspec(1,:,6),'linewidth',2,'color','g','linestyle','-');

plot(T2,nspec1(:,1),'linewidth',2,'color','k','linestyle','--');
plot(T2,nspec1(:,2),'linewidth',2,'color','r','linestyle','--');
plot(T2,nspec1(:,3),'linewidth',2,'color','b','linestyle','--');
plot(T2,nspec1(:,4),'linewidth',2,'color','c','linestyle','--');
plot(T2,nspec1(:,5),'linewidth',2,'color','m','linestyle','--');
plot(T2,nspec1(:,6),'linewidth',2,'color','g','linestyle','--');