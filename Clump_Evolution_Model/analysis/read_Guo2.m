function [Nclumps, Ngals] = read_Guo2(zbin, Mbin, subtraction, ind, is, norm_is, nis)
% zbin = [1.0 3.0];
% Mbin = [9.8 11.4];
% subtraction = 4;

mass_dist = false;
sfr_dist = false;
tsf_dist = false;
age_dist = false;
mass_age = false;
sfr_age = false;
tsf_age = false;
%%%%%%%%%%%%%%%%%%%%%%%%%
Tdown7 = 50;
Tup7 = 3000;
Tdown19 = 50;
Tup19 = 3000;
Mdown7 = 8.0;
Mup7 = 9.5;
Mdown19 = 8.0;
Mup19 = 9.5;
zdown7 = 1;
zup7 = 2.5;
zdown19 = 3;
zup19 = 5;

SFRdown = -1;
%%%%%%%%%%%%%%%%%%%%%%%%%

% Collect V07 clumps using same cuts as in clump evolution plots
m = find(ind==7);
vec = (nis(m)+1):(nis(m+1)-1);
vec2 = vec+1;
b = find(is(vec,62)>=Tdown7 & is(vec,62)<Tup7 & is(vec,1)~=is(vec2,1));
end_index7 = vec(b);
if( is(nis(m+1),62)>=Tdown7 & is(nis(m+1),62)<Tup7 )
    end_index7 = [end_index7, nis(m+1)];
end
b = find(is(vec2,62)>=Tdown7 & is(vec2,62)<Tup7 & is(vec2,1)~=is(vec,1));
start_index7 = vec2(b);
if( is(nis(m)+1,62)>=Tdown7 & is(nis(m)+1,62)<Tup7 )
    start_index7 = [nis(m)+1, start_index7];
end
b = find(is(start_index7,2)>=zdown7 & is(start_index7,2)<zup7);
start_index7 = start_index7(b);
end_index7 = end_index7(b);
Mavg = 0.*end_index7;
for i=1:length(start_index7)
    tvec = start_index7(i):end_index7(i);
    b = find(abs(is(tvec,48)-50)<=10);
    if(isempty(b))
        b = find(abs(is(tvec,48)-50)==min(abs(is(tvec,48)-50)));
    end
    Mavg(i) = mean(log10(is(tvec(b),6)));
end
b = find(Mavg>=Mdown7 & Mavg<Mup7);
start_index7 = start_index7(b);
end_index7 = end_index7(b);

if( length(start_index7) ~= length(end_index7) )
    [0,0,length(start_index7), length(end_index7)]
    return
elseif( min(end_index7-start_index7)<=0 )
    [0,0,min(end_index7-start_index7)]
    return
else
    Nclump = length(start_index7);
    clump_length = end_index7(1:Nclump) - start_index7(1:Nclump) + 1;
    Nmax = sum(clump_length);    
    clump_length = [0, clump_length];
    sim_age7   = zeros(1,Nmax);
    sim_mstar7 = -50.*ones(1,Nmax);
    sim_sfr7   = -50.*ones(1,Nmax);
    sim_sfr72  = -50.*ones(1,Nmax);
    sim_tsf7   = zeros(1,Nmax);
    sim_dist7  = zeros(1,Nmax);
    for i=1:Nclump
        j = start_index7(i);
        k = end_index7(i);
        begin = sum(clump_length(1:i));
        sim_age7(begin+1:begin+clump_length(i+1))   = is(j:k,12);
        sim_mstar7(begin+1:begin+clump_length(i+1)) = log10(norm_is(j:k,5));
        sim_sfr7(begin+1:begin+clump_length(i+1))   = log10(norm_is(j:k,15));
        sim_sfr72(begin+1:begin+clump_length(i+1))   = log10(is(j:k,15));
        sim_tsf7(begin+1:begin+clump_length(i+1))   = 1000./is(j:k,17);
        sim_dist7(begin+1:begin+clump_length(i+1))  = is(j:k,19);
    end
end

% Collect V19 clumps using same cuts as in clump evolution plots
m = find(ind==19);
vec = (nis(m)+1):(nis(m+1)-1);
vec2 = vec+1;
b = find(is(vec,62)>=Tdown19 & is(vec,62)<Tup19 & is(vec,1)~=is(vec2,1));
end_index19 = vec(b);
if( is(nis(m+1),62)>=Tdown19 & is(nis(m+1),62)<Tup19 )
    end_index19 = [end_index19, nis(m+1)];
end
b = find(is(vec2,62)>=Tdown19 & is(vec2,62)<Tup19 & is(vec2,1)~=is(vec,1));
start_index19 = vec2(b);
if( is(nis(m)+1,62)>=Tdown19 & is(nis(m)+1,62)<Tup19 )
    start_index19 = [nis(m)+1, start_index19];
end
b = find(is(start_index19,2)>=zdown19 & is(start_index19,2)<zup19);
start_index19 = start_index19(b);
end_index19 = end_index19(b);
Mavg = 0.*end_index19;
for i=1:length(start_index19)
    tvec = start_index19(i):end_index19(i);
    b = find(abs(is(tvec,48)-50)<=10);
    if(isempty(b))
        b = find(abs(is(tvec,48)-50)==min(abs(is(tvec,48)-50)));
    end
    Mavg(i) = mean(log10(is(tvec(b),6)));
end
b = find(Mavg>=Mdown19 & Mavg<Mup19);
start_index19 = start_index19(b);
end_index19 = end_index19(b);

if( length(start_index19) ~= length(end_index19) )
    [0,0,length(start_index19), length(end_index19)]
    return
elseif( min(end_index19-start_index19)<=0 )
    [0,0,min(end_index19-start_index19)]
    return
else
    Nclump = length(start_index19);
    clump_length = end_index19(1:Nclump) - start_index19(1:Nclump) + 1;
    Nmax = sum(clump_length);    
    clump_length = [0, clump_length];
    sim_age19   = zeros(1,Nmax);
    sim_mstar19 = -50.*ones(1,Nmax);
    sim_sfr19   = -50.*ones(1,Nmax);
    sim_sfr192  = -50.*ones(1,Nmax);
    sim_tsf19   = zeros(1,Nmax);
    sim_dist19  = zeros(1,Nmax);
    for i=1:Nclump
        j = start_index19(i);
        k = end_index19(i);
        begin = sum(clump_length(1:i));
        sim_age19(begin+1:begin+clump_length(i+1))   = is(j:k,12);
        sim_mstar19(begin+1:begin+clump_length(i+1)) = log10(norm_is(j:k,5));
        sim_sfr19(begin+1:begin+clump_length(i+1))   = log10(norm_is(j:k,15));
        sim_sfr192(begin+1:begin+clump_length(i+1))  = log10(is(j:k,15));
        sim_tsf19(begin+1:begin+clump_length(i+1))   = 1000./is(j:k,17);
        sim_dist19(begin+1:begin+clump_length(i+1))  = is(j:k,19);
    end
end

% b = find(ind==7);
% sim7 = is(nis(b)+1:nis(b+1),:);
% norm_sim7 = norm_is(nis(b)+1:nis(b+1),:);
% Mdisc7 = sim7(:,6)./norm_sim7(:,6);
% Msdisc7 = sim7(:,5)./norm_sim7(:,5);
% %b = find(sim7(:,2)>=zbin(1) & sim7(:,2)<=zbin(2) & log10(Mdisc7) >= Mbin(1) & log10(Mdisc7) <= Mbin(2));
% b = find(sim7(:,2)>=1 & sim7(:,2)<=3 & sim7(:,62)>10 & log10(sim7(:,6))>=8 & log10(sim7(:,15))>=-1);
% sim7 = sim7(b,:);
% norm_sim7  = norm_sim7(b,:);
% sim_age7   = sim7(:,12);
% %sim_mstar7 = log10(norm_sim7(:,5) .* Msdisc7(b) ./ Mdisc7(b));
% sim_mstar7 = log10(norm_sim7(:,5));
% sim_sfr7   = log10(norm_sim7(:,15));
% sim_tsf7   = 1000./sim7(:,17);
% sim_dist7  = sim7(:,19);
% 
% b = find(ind==19);
% sim19 = is(nis(b)+1:nis(b+1),:);
% norm_sim19 = norm_is(nis(b)+1:nis(b+1),:);
% Mdisc19 = sim19(:,6)./norm_sim19(:,6);
% Msdisc19 = sim19(:,5)./norm_sim19(:,5);
% %b = find(sim19(:,2)>=zbin(1) & sim19(:,2)<=zbin(2) & log10(Mdisc19) >= Mbin(1) & log10(Mdisc19) <= Mbin(2));
% b = find(sim19(:,2)>=3 & sim19(:,2)<=5 & sim19(:,62)>10 & log10(sim19(:,6))>=8 & log10(sim19(:,15))>=-1);
% sim19 = sim19(b,:);
% norm_sim19  = norm_sim19(b,:);
% sim_age19   = sim19(:,12);
% %sim_mstar19 = log10(norm_sim19(:,5) .* Msdisc19(b) ./ Mdisc19(b));
% sim_mstar19 = log10(norm_sim19(:,5));
% sim_sfr19   = log10(norm_sim19(:,15));
% sim_tsf19   = 1000./sim19(:,17);
% sim_dist19  = sim19(:,19);

a=load('Guo_clump_cat_edited.txt');

b=find(a(:,19)==0);
a0 = a(b,:);
b=find(a(:,19)==1);
a1 = a(b,:);
b=find(a(:,19)==2);
a2 = a(b,:);
b=find(a(:,19)==3);
a3 = a(b,:);
b=find(a(:,19)==4);
a4 = a(b,:);
b=find(a(:,19)==5);
a5 = a(b,:);
b=find(a(:,19)==6);
a6 = a(b,:);
[length(a0), length(a1), length(a2), length(a3), length(a4), length(a5), length(a6)]
clear a b

com0=[];
com1=[];
com2=[];
com3=[];
com4=[];
com5=[];
com6=[];
for i=1:length(a1)
    com5 = [com5, find(a5(:,1)==a1(i,1) & a5(:,14)==a1(i,14))];
end
for i=1:length(com5)
    com0 = [com0, find(a0(:,1)==a5(com5(i),1) & a0(:,14)==a5(com5(i),14))];
    com1 = [com1, find(a1(:,1)==a5(com5(i),1) & a1(:,14)==a5(com5(i),14))];
    com2 = [com2, find(a2(:,1)==a5(com5(i),1) & a2(:,14)==a5(com5(i),14))];
    com3 = [com3, find(a3(:,1)==a5(com5(i),1) & a3(:,14)==a5(com5(i),14))];
    com4 = [com4, find(a4(:,1)==a5(com5(i),1) & a4(:,14)==a5(com5(i),14))];
    com6 = [com6, find(a6(:,1)==a5(com5(i),1) & a6(:,14)==a5(com5(i),14))];
end
a0 = a0(com0,:);
a1 = a1(com1,:);
a2 = a2(com2,:);
a3 = a3(com3,:);
a4 = a4(com4,:);
a5 = a5(com5,:);
a6 = a6(com6,:);
clear com0 com1 com2 com3 com4 com5 com6 i
% [max(abs(a0(:,1)-a1(:,1))), max(abs(a0(:,1)-a2(:,1))), ...
%     max(abs(a0(:,1)-a3(:,1))), max(abs(a0(:,1)-a4(:,1))), ...
%     max(abs(a0(:,1)-a5(:,1))), max(abs(a0(:,1)-a6(:,1)))]
% 
% [max(abs(a0(:,14)-a1(:,14))), max(abs(a0(:,14)-a2(:,14))), ...
%     max(abs(a0(:,14)-a3(:,14))), max(abs(a0(:,14)-a4(:,14))), ...
%     max(abs(a0(:,14)-a5(:,14))), max(abs(a0(:,14)-a6(:,14)))]
N = length(a0);
a = [a0
    a1
    a2
    a3
    a4
    a5
    a6];
clear a1 a2 a3 a4 a5 a6

ind = (subtraction*N+1):((subtraction+1)*N);
a0 = a(ind,:);
b = find(a0(:,4)>=zbin(1) & a0(:,4)<=zbin(2) & a0(:,5)>=Mbin(1) & a0(:,5)<=Mbin(2));
a0 = a0(b,:);
Nclumps = length(b);
Ngals = length(unique(a0(:,1)));

Mclump = a0(:,41);                      %Log(normalized to galaxy
Mclump_norm = Mclump-a0(:,5);           %Log(normalized to galaxy
%SFR = a0(:,44);                        %SED SFR Log(normalized to galaxy
SFR = a0(:,62) - log10(1.8);            %Log UV SFR converted to Salpeter
b = find(SFR>-10);
[max(SFR(b)), min(SFR(b))]
SFR_norm = SFR-a0(:,7);                 %Log(normalized to galaxy)
tsf = 1e-6.*(10.^(Mclump-SFR));  %Myr
tsf_norm = Mclump_norm-SFR_norm;        %Log(normalized to galaxy
age = 1e3.*(10.^(a0(:,50)));            %Myr
tau = 1e3.*(10.^(a0(:,53)));            %Myr
dist = a0(:,38);                        %Normalized by galactic SMA
[length(find(Mclump<=-20)),length(find(SFR<=-20))]
dirname = './Guo_data/';
mkdir(dirname);
dirname = strcat(dirname,'redshift_',num2str(zbin(1),'%3.1f'),'_',num2str(zbin(2),'%3.1f'),...
    '_Mstar_',num2str(Mbin(1),'%3.1f'),'_',num2str(Mbin(2),'%3.1f'),'/');
mkdir(dirname);
dirname = strcat(dirname,'subtraction_',num2str(subtraction,'%1i'),'/');
mkdir(dirname);

median(age)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fgi_model = 0.97;
alpha_model = 0.3;
eta_model = 1.0;
eps_model = 0.1;
mu_model = 0.8;
etas_model = 0.3;
tmig_td_model = 5;
tdyn_model = 30;

[T1,Y1]=ode45(@(t1,y1)clump_evolution_exact(t1, y1, alpha_model, eta_model, ...
    eps_model, mu_model, etas_model, 1, tmig_td_model),[0 300],[fgi_model 1-fgi_model 1]);
time_model = tdyn_model.*T1(:,1);
Mg_model = Y1(:,1);
Ms_model = Y1(:,2);
%R_model = Y1(:,3);
%R_model = 1 - T1(:,1)./tmig_td_model;
R_model = 2.5.*exp(-T1(:,1)./tmig_td_model);
Mb_model = Mg_model + Ms_model;
fg_model = (Mg_model./Mb_model);
SFR_norm_model = eps_model.*Mg_model.*(10^8)./(1e6.*tdyn_model);
sSFR_norm_model = 5e8.*SFR_norm_model./((10^8).*Ms_model); % = 1000.*(Mg./Ms).*eps.*td_tff./tdyn;
tdep_model = (Mg_model.*(10^8))./(1e6.*SFR_norm_model); % = tdyn./(eps.*td_tff)
tsf_model = 1000./sSFR_norm_model; % = ( tdyn./(eps.*td_tff) ) .* ( Ms./Mg )

% R_model = 1.70.*R_model;
%Ms_model = log10(0.05.*Ms_model);              %tdyn=30Myr
%SFR_norm_model = log10(0.10.*SFR_norm_model);  %tdyn=30Myr
Ms_model = log10(0.075.*Ms_model);              %tdyn=50Myr
SFR_norm_model = log10(0.2.*SFR_norm_model);   %tdyn=50Myr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(tsf_dist) 
    xl = 0.15;
    xu = 5;
    xtick = [0.2 1 5];
    xstr = {num2str(0.2,'%3.1f'),num2str(1,'%01i'),num2str(5,'%01i')};
    yl = 1;
    yu = 2000;
    ytick = [1 10 100 1000 2000];
    ystr = {num2str(1,'%01i'),num2str(10,'%02i'),num2str(100,'%03i'),num2str(1000,'%04i'),num2str(2000,'%04i')};
    %ystr = {'10^0','10^1','10^2','10^3','10^4'};
    
    ylab  = '$t_{\rm SF}\:{\rm [Myr]}$';
    ylab2 = '$t_{\rm SF}\:\:{\rm [Myr\:]}$';
    xlab  = '$d\:/\:{\rm SMA}$';
    xlab2 = '$d\:\:/\:\:{\rm SMA}$';
    
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
        'Position',[0.13 0.14 0.775 0.815]);
    xlim(axes1,[xl xu]);
    ylim(axes1,[yl yu]);
    box(axes1,'on');
    % grid on
    hold(axes1,'all');
    xlabel(xlab,'Interpreter','latex','FontSize',22,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);
    ylabel(ylab,'Interpreter','latex','FontSize',22,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.11 0.5 0]);
    set(gcf,'renderer','painters')
    
    plot(dist, tsf, 'linestyle','none','marker','o','markerfacecolor','k','markeredgecolor','k','markersize',3)
    
    xbin = -0.5:0.10:0.4;
    dbin = (xbin(2)-xbin(1))/2;
    Nbin = length(xbin);
    tsfbin = zeros(Nbin,3);
    sim_tsfbin7  = zeros(Nbin,3);
    sim_tsfbin19 = zeros(Nbin,3);
    for i=1:Nbin
        b = find(SFR_norm>=-9 & log10(dist)>=xbin(i)-dbin & log10(dist)<xbin(i)+dbin);
        temp = sort(tsf(b));
        Ntemp = length(temp);
        if(Ntemp>10)
            N1 = ceil(0.16*Ntemp);
            N2 = ceil(0.84*Ntemp);
            tsfbin(i,1) = median(tsf(b));
            tsfbin(i,2) = median(temp(N1));
            tsfbin(i,3) = median(temp(N2));
        end
        b = find(sim_sfr72>=SFRdown & log10(sim_dist7)>=xbin(i)-dbin & log10(sim_dist7)<xbin(i)+dbin);
        temp = sort(sim_tsf7(b));
        Ntemp = length(temp);
        if(Ntemp>10)
            N1 = ceil(0.16*Ntemp);
            N2 = ceil(0.84*Ntemp);
            sim_tsfbin7(i,1) = median(sim_tsf7(b));
            sim_tsfbin7(i,2) = median(temp(N1));
            sim_tsfbin7(i,3) = median(temp(N2));
        end
        b = find(sim_sfr192>=SFRdown & log10(sim_dist19)>=xbin(i)-dbin & log10(sim_dist19)<xbin(i)+dbin);
        temp = sort(sim_tsf19(b));
        Ntemp = length(temp);
        if(Ntemp>10)
            N1 = ceil(0.16*Ntemp);
            N2 = ceil(0.84*Ntemp);
            sim_tsfbin19(i,1) = median(sim_tsf19(b));
            sim_tsfbin19(i,2) = median(temp(N1));
            sim_tsfbin19(i,3) = median(temp(N2));
        end
    end
    xbin = 10.^xbin;
    l(1) = plot(xbin, tsfbin(:,1), 'linestyle','-','marker','o',...
        'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8,...
        'color','k','linewidth',4,'DisplayName','median obs');
    l(4) = plot(R_model, tsf_model, 'linestyle','-','marker','none',...
        'color','c','linewidth',4,'DisplayName','model');
    l(2) = plot(xbin, sim_tsfbin7(:,1), 'linestyle','-','marker','o',...
        'MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',6,...
        'color','g','linewidth',3,'DisplayName','median V07');
    l(3) = plot(xbin, sim_tsfbin19(:,1), 'linestyle','-','marker','o',...
        'MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',6,...
        'color','m','linewidth',3,'DisplayName','median V19');
    plot(xbin, tsfbin(:,2), 'linestyle','--','marker','none','color','k','linewidth',4);
    plot(xbin, tsfbin(:,3), 'linestyle','--','marker','none','color','k','linewidth',4);
%     plot(xbin, sim_tsfbin7(:,2), 'linestyle','--','marker','none','color','c','linewidth',3);
%     plot(xbin, sim_tsfbin7(:,3), 'linestyle','--','marker','none','color','c','linewidth',3);
%     plot(xbin, sim_tsfbin19(:,2), 'linestyle','--','marker','none','color','m','linewidth',3);
%     plot(xbin, sim_tsfbin19(:,3), 'linestyle','--','marker','none','color','m','linewidth',3);
    legend1 = legend(gca,l(1:3));
    set(legend1,...
        'Position',[0.24, 0.135, 0.10, 0.15],'FontSize',14,...
        'Interpreter','tex');
    legend boxoff
    
    print(gcf,'-depsc',strcat(dirname,'tsf_dist.eps'));
    set(gcf,'visible','off');
    xlabel(xlab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);
    ylabel(ylab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.11 0.5 0]);
    
    ax = gca;
    %     outerpos = ax.OuterPosition;
    %     ti = ax.TightInset;
    %     left = outerpos(1) + ti(1) + 0.02;
    %     bottom = outerpos(2) + ti(2) + 0.02;
    %     ax_width = outerpos(3) - ti(1) - ti(3) - 0.04;
    %     ax_height = outerpos(4) - ti(2) - ti(4) - 0.06;
    %     ax.Position = [left bottom ax_width ax_height];
    ax.Position = [0.1076 0.1200 0.8724 0.8524];
    
    print(gcf,'-djpeg',strcat(dirname,'tsf_dist.jpg'));
    print(gcf,'-depsc',strcat(dirname,'tsf_dist.eps'));
    close all
    fclose all;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(age_dist)
    xl = 0.15;
    xu = 5;
    xtick = [0.2 1 5];
    xstr = {num2str(0.2,'%3.1f'),num2str(1,'%01i'),num2str(5,'%01i')};
    yl = 10;
    yu = 1000;
    ytick = [10 100 1000 10000];
    ystr = {num2str(10,'%02i'),num2str(100,'%03i'),num2str(1000,'%04i'),num2str(10000,'%05i')};
    %ystr = {'10^1','10^2','10^3','10^4'};
    
    ylab  = '${\rm age}\:{\rm [Myr]}$';
    ylab2 = '${\rm age}\:\:{\rm [Myr]}$';
    xlab  = '$d\:/\:{\rm SMA}$';
    xlab2 = '$d\:\:/\:\:{\rm SMA}$';
    
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
        'Position',[0.13 0.14 0.775 0.815]);
    xlim(axes1,[xl xu]);
    ylim(axes1,[yl yu]);
    box(axes1,'on');
    % grid on
    hold(axes1,'all');
    xlabel(xlab,'Interpreter','latex','FontSize',22,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);
    ylabel(ylab,'Interpreter','latex','FontSize',22,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.11 0.5 0]);
    set(gcf,'renderer','painters')
    
    plot(dist, age, 'linestyle','none','marker','o','markerfacecolor','k','markeredgecolor','k','markersize',3)
    
    xbin = -0.5:0.10:0.4;
    dbin = (xbin(2)-xbin(1))/2;
    Nbin = length(xbin);
    agebin = zeros(Nbin,3);
    sim_agebin7 = zeros(Nbin,3);
    sim_agebin19 = zeros(Nbin,3);
    for i=1:Nbin
        b = find(SFR_norm>=-9 & log10(dist)>=xbin(i)-dbin & log10(dist)<xbin(i)+dbin);
        temp = sort(age(b));
        Ntemp = length(temp);
        if(Ntemp>10)
            N1 = ceil(0.16*Ntemp);
            N2 = ceil(0.84*Ntemp);
            agebin(i,1) = median(age(b));
            agebin(i,2) = median(temp(N1));
            agebin(i,3) = median(temp(N2));
        end
        b = find(sim_sfr72>=SFRdown & log10(sim_dist7)>=xbin(i)-dbin & log10(sim_dist7)<xbin(i)+dbin);
        temp = sort(sim_age7(b));
        Ntemp = length(temp);
        if(Ntemp>10)
            N1 = ceil(0.16*Ntemp);
            N2 = ceil(0.84*Ntemp);
            sim_agebin7(i,1) = median(sim_age7(b));
            sim_agebin7(i,2) = median(temp(N1));
            sim_agebin7(i,3) = median(temp(N2));
        end
        b = find(sim_sfr192>=SFRdown & log10(sim_dist19)>=xbin(i)-dbin & log10(sim_dist19)<xbin(i)+dbin);
        temp = sort(sim_age19(b));
        Ntemp = length(temp);
        if(Ntemp>10)
            N1 = ceil(0.16*Ntemp);
            N2 = ceil(0.84*Ntemp);
            sim_agebin19(i,1) = median(sim_age19(b));
            sim_agebin19(i,2) = median(temp(N1));
            sim_agebin19(i,3) = median(temp(N2));
        end
    end
    xbin = 10.^xbin;
    l(1) = plot(xbin, agebin(:,1), 'linestyle','-','marker','o',...
        'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8,...
        'color','k','linewidth',4,'DisplayName','median obs');
    l(2) = plot(xbin, sim_agebin7(:,1), 'linestyle','-','marker','o',...
        'MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',6,...
        'color','g','linewidth',3,'DisplayName','median V07');
    l(3) = plot(xbin, sim_agebin19(:,1), 'linestyle','-','marker','o',...
        'MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',6,...
        'color','m','linewidth',3,'DisplayName','median V19');
    plot(xbin, agebin(:,2), 'linestyle','--','marker','none','color','k','linewidth',4);
    plot(xbin, agebin(:,3), 'linestyle','--','marker','none','color','k','linewidth',4);
%     plot(xbin, sim_agebin7(:,2), 'linestyle','--','marker','none','color','c','linewidth',3);
%     plot(xbin, sim_agebin7(:,3), 'linestyle','--','marker','none','color','c','linewidth',3);
%     plot(xbin, sim_agebin19(:,2), 'linestyle','--','marker','none','color','m','linewidth',3);
%     plot(xbin, sim_agebin19(:,3), 'linestyle','--','marker','none','color','m','linewidth',3);
    legend1 = legend(gca,l(1:3));
    set(legend1,...
        'Position',[0.24, 0.135, 0.10, 0.15],'FontSize',14,...
        'Interpreter','tex');
    legend boxoff
    
    print(gcf,'-depsc',strcat(dirname,'age_dist.eps'));
    set(gcf,'visible','off');
    xlabel(xlab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);
    ylabel(ylab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.11 0.5 0]);
    
    ax = gca;
    %     outerpos = ax.OuterPosition;
    %     ti = ax.TightInset;
    %     left = outerpos(1) + ti(1) + 0.02;
    %     bottom = outerpos(2) + ti(2) + 0.02;
    %     ax_width = outerpos(3) - ti(1) - ti(3) - 0.04;
    %     ax_height = outerpos(4) - ti(2) - ti(4) - 0.06;
    %     ax.Position = [left bottom ax_width ax_height];
    ax.Position = [0.1076 0.1200 0.8724 0.8524];
    
    print(gcf,'-djpeg',strcat(dirname,'age_dist.jpg'));
    print(gcf,'-depsc',strcat(dirname,'age_dist.eps'));
    close all
    fclose all;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(mass_dist)
    xl = 0.15;
    xu = 5;
    xtick = [0.2 1 5];
    xstr = {num2str(0.2,'%3.1f'),num2str(1,'%01i'),num2str(5,'%01i')};
    yl = 1e-4;
    yu = 1;
    ytick = [1e-4 1e-3 0.01 0.1 1];
    % ystr = {num2str(0.0001,'%6.4f'),num2str(0.001,'%5.3f'),num2str(0.01,'%4.2f'),num2str(0.1,'%3.1'),num2str(1,'%1.0f')};
    ystr = {'10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^0'};
    
    ylab  = '$M_*\:/\:M^{\rm gal}_*$';
    ylab2 = '$M_*\:\:/\:\:M^{\rm gal}_*$';
    xlab  = '$d\:/\:{\rm SMA}$';
    xlab2 = '$d\:\:/\:\:{\rm SMA}$';
    
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
        'Position',[0.13 0.14 0.775 0.815]);
    xlim(axes1,[xl xu]);
    ylim(axes1,[yl yu]);
    box(axes1,'on');
    % grid on
    hold(axes1,'all');
    xlabel(xlab,'Interpreter','latex','FontSize',22,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);
    ylabel(ylab,'Interpreter','latex','FontSize',22,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.11 0.5 0]);
    set(gcf,'renderer','painters')
    
    b=find(Mclump_norm>=-3.6);
    plot(dist(b), 10.^Mclump_norm(b), 'linestyle','none','marker','o','markerfacecolor','k','markeredgecolor','k','markersize',3)
    
    xbin = -0.5:0.10:0.4;
    dbin = (xbin(2)-xbin(1))/2;
    Nbin = length(xbin);
    Mclump_norm_bin = -50.*ones(Nbin,3);
    sim_Mclump_norm_bin7 = -50.*ones(Nbin,3);
    sim_Mclump_norm_bin19 = -50.*ones(Nbin,3);
    for i=1:Nbin
        b = find(Mclump_norm>=-9 & log10(dist)>=xbin(i)-dbin & log10(dist)<xbin(i)+dbin);
        temp = sort(Mclump_norm(b));
        Ntemp = length(temp);
        if(Ntemp>10)
            N1 = ceil(0.16*Ntemp);
            N2 = ceil(0.84*Ntemp);
            Mclump_norm_bin(i,1) = median(Mclump_norm(b));
            Mclump_norm_bin(i,2) = median(temp(N1));
            Mclump_norm_bin(i,3) = median(temp(N2));
        end
        b = find(sim_sfr72>=SFRdown & log10(sim_dist7)>=xbin(i)-dbin & log10(sim_dist7)<xbin(i)+dbin);
        temp = sort(sim_mstar7(b));
        Ntemp = length(temp);
        if(Ntemp>10)
            N1 = ceil(0.16*Ntemp);
            N2 = ceil(0.84*Ntemp);
            sim_Mclump_norm_bin7(i,1) = median(sim_mstar7(b));
            sim_Mclump_norm_bin7(i,2) = median(temp(N1));
            sim_Mclump_norm_bin7(i,3) = median(temp(N2));
        end
        b = find(sim_sfr192>=SFRdown & log10(sim_dist19)>=xbin(i)-dbin & log10(sim_dist19)<xbin(i)+dbin);
        temp = sort(sim_mstar19(b));
        Ntemp = length(temp);
        if(Ntemp>10)
            N1 = ceil(0.16*Ntemp);
            N2 = ceil(0.84*Ntemp);
            sim_Mclump_norm_bin19(i,1) = median(sim_mstar19(b));
            sim_Mclump_norm_bin19(i,2) = median(temp(N1));
            sim_Mclump_norm_bin19(i,3) = median(temp(N2));
        end
    end
    xbin = 10.^xbin;
    b1 = find(Mclump_norm_bin(:,1)>-50);
    b2 = find(sim_Mclump_norm_bin7(:,1)>-50);
    b3 = find(sim_Mclump_norm_bin19(:,1)>-50);
    l(1) = plot(xbin(b1), 10.^Mclump_norm_bin(b1,1), 'linestyle','-','marker','o',...
        'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8,...
        'color','k','linewidth',4,'DisplayName','median obs');
    l(4) = plot(R_model, 10.^Ms_model, 'linestyle','-','marker','none',...
        'color','c','linewidth',4,'DisplayName','model');
    l(2) = plot(xbin(b2), 10.*10.^sim_Mclump_norm_bin7(b2,1), 'linestyle','-','marker','o',...
        'MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',6,...
        'color','g','linewidth',3,'DisplayName','median V07 x10');
    l(3) = plot(xbin(b3), 10.*10.^sim_Mclump_norm_bin19(b3,1), 'linestyle','-','marker','o',...
        'MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',6,...
        'color','m','linewidth',3,'DisplayName','median V19 x10');
    plot(xbin(b1), 10.^Mclump_norm_bin(b1,2), 'linestyle','--','marker','none','color','k','linewidth',4);
    plot(xbin(b1), 10.^Mclump_norm_bin(b1,3), 'linestyle','--','marker','none','color','k','linewidth',4);
%     plot(xbin(b2), 10.^sim_Mclump_norm_bin7(b2,2), 'linestyle','--','marker','none','color','c','linewidth',3);
%     plot(xbin(b2), 10.^sim_Mclump_norm_bin7(b2,3), 'linestyle','--','marker','none','color','c','linewidth',3);
%     plot(xbin(b3), 10.^sim_Mclump_norm_bin19(b3,2), 'linestyle','--','marker','none','color','m','linewidth',3);
%     plot(xbin(b3), 10.^sim_Mclump_norm_bin19(b3,3), 'linestyle','--','marker','none','color','m','linewidth',3);
    legend1 = legend(gca,l(1:3));
    set(legend1,...
        'Position',[0.30, 0.135, 0.10, 0.15],'FontSize',14,...
        'Interpreter','tex');
    legend boxoff
    
    print(gcf,'-depsc',strcat(dirname,'Mstar_norm_dist.eps'));
    set(gcf,'visible','off');
    xlabel(xlab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);
    ylabel(ylab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.11 0.5 0]);
    
    ax = gca;
    %     outerpos = ax.OuterPosition;
    %     ti = ax.TightInset;
    %     left = outerpos(1) + ti(1) + 0.02;
    %     bottom = outerpos(2) + ti(2) + 0.02;
    %     ax_width = outerpos(3) - ti(1) - ti(3) - 0.04;
    %     ax_height = outerpos(4) - ti(2) - ti(4) - 0.06;
    %     ax.Position = [left bottom ax_width ax_height];
    ax.Position = [0.1076 0.1200 0.8724 0.8524];
    
    print(gcf,'-djpeg',strcat(dirname,'Mstar_norm_dist.jpg'));
    print(gcf,'-depsc',strcat(dirname,'Mstar_norm_dist.eps'));
    close all
    fclose all;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(sfr_dist)
    xl = 0.15;
    xu = 5;
    xtick = [0.2 1 5];
    xstr = {num2str(0.2,'%3.1f'),num2str(1,'%01i'),num2str(5,'%01i')};
    yl = 1e-3;
    yu = 10;
    ytick = [1e-3 0.01 0.1 1 10];
    % ystr = {num2str(0.0001,'%6.4f'),num2str(0.001,'%5.3f'),num2str(0.01,'%4.2f'),num2str(0.1,'%3.1'),num2str(1,'%1.0f')};
    ystr = {'10^{-3}','10^{-2}','10^{-1}','10^0','10^1'};
    
    ylab  = '${\rm SFR}\:/\:{\rm SFR^{gal}}$';
    ylab2 = '${\rm SFR}\:\:/\:\:{\rm SFR^{gal}}$';
    xlab  = '$d\:/\:{\rm SMA}$';
    xlab2 = '$d\:\:/\:\:{\rm SMA}$';
    
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
        'Position',[0.13 0.14 0.775 0.815]);
    xlim(axes1,[xl xu]);
    ylim(axes1,[yl yu]);
    box(axes1,'on');
    % grid on
    hold(axes1,'all');
    xlabel(xlab,'Interpreter','latex','FontSize',22,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);
    ylabel(ylab,'Interpreter','latex','FontSize',22,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.11 0.5 0]);
    set(gcf,'renderer','painters')
    
    b=find(SFR_norm>=-3.6);
    plot(dist(b), 10.^SFR_norm(b), 'linestyle','none','marker','o','markerfacecolor','k','markeredgecolor','k','markersize',3)
    
    xbin = -0.5:0.10:0.4;
    dbin = (xbin(2)-xbin(1))/2;
    Nbin = length(xbin);
    SFR_norm_bin = -50.*ones(Nbin,3);
    sim_SFR_norm_bin7 = -50.*ones(Nbin,3);
    sim_SFR_norm_bin19 = -50.*ones(Nbin,3);
    for i=1:Nbin
        b = find(SFR_norm>=-9 & log10(dist)>=xbin(i)-dbin & log10(dist)<xbin(i)+dbin);
        temp = sort(SFR_norm(b));
        Ntemp = length(temp);
        if(Ntemp>10)
            N1 = ceil(0.16*Ntemp);
            N2 = ceil(0.84*Ntemp);
            SFR_norm_bin(i,1) = median(SFR_norm(b));
            SFR_norm_bin(i,2) = median(temp(N1));
            SFR_norm_bin(i,3) = median(temp(N2));
        end
        b = find(sim_sfr72>=SFRdown & log10(sim_dist7)>=xbin(i)-dbin & log10(sim_dist7)<xbin(i)+dbin);
        temp = sort(sim_sfr7(b));
        Ntemp = length(temp);
        if(Ntemp>10)
            N1 = ceil(0.16*Ntemp);
            N2 = ceil(0.84*Ntemp);
            sim_SFR_norm_bin7(i,1) = median(sim_sfr7(b));
            sim_SFR_norm_bin7(i,2) = median(temp(N1));
            sim_SFR_norm_bin7(i,3) = median(temp(N2));
        end
        b = find(sim_sfr192>=SFRdown & log10(sim_dist19)>=xbin(i)-dbin & log10(sim_dist19)<xbin(i)+dbin);
        temp = sort(sim_sfr19(b));
        Ntemp = length(temp);
        if(Ntemp>10)
            N1 = ceil(0.16*Ntemp);
            N2 = ceil(0.84*Ntemp);
            sim_SFR_norm_bin19(i,1) = median(sim_sfr19(b));
            sim_SFR_norm_bin19(i,2) = median(temp(N1));
            sim_SFR_norm_bin19(i,3) = median(temp(N2));
        end
    end
    xbin = 10.^xbin;
    b1 = find(SFR_norm_bin(:,1)>-50);
    b2 = find(sim_SFR_norm_bin7(:,1)>-50);
    b3 = find(sim_SFR_norm_bin19(:,1)>-50);
    l(1) = plot(xbin(b1), 10.^SFR_norm_bin(b1,1), 'linestyle','-','marker','o',...
        'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8,...
        'color','k','linewidth',4,'DisplayName','median obs');
    l(2) = plot(xbin(b2), 3*10.^sim_SFR_norm_bin7(b2,1), 'linestyle','-','marker','o',...
        'MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',6,...
        'color','g','linewidth',3,'DisplayName','median V07 x4');
    l(3) = plot(xbin(b3), 3*10.^sim_SFR_norm_bin19(b3,1), 'linestyle','-','marker','o',...
        'MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',6,...
        'color','m','linewidth',3,'DisplayName','median V19 x4');
    plot(xbin(b1), 10.^SFR_norm_bin(b1,2), 'linestyle','--','marker','none','color','k','linewidth',3);
    plot(xbin(b1), 10.^SFR_norm_bin(b1,3), 'linestyle','--','marker','none','color','k','linewidth',3);
%     plot(xbin(b2), 10.^sim_SFR_norm_bin7(b2,2), 'linestyle','--','marker','none','color','c','linewidth',3);
%     plot(xbin(b2), 10.^sim_SFR_norm_bin7(b2,3), 'linestyle','--','marker','none','color','c','linewidth',3);
%     plot(xbin(b3), 10.^sim_SFR_norm_bin19(b3,2), 'linestyle','--','marker','none','color','m','linewidth',3);
%     plot(xbin(b3), 10.^sim_SFR_norm_bin19(b3,3), 'linestyle','--','marker','none','color','m','linewidth',3);
    legend1 = legend(gca,l(1:3));
    set(legend1,...
        'Position',[0.30, 0.135, 0.10, 0.15],'FontSize',14,...
        'Interpreter','tex');
    legend boxoff
    
    print(gcf,'-depsc',strcat(dirname,'sfr_norm_dist.eps'));
    set(gcf,'visible','off');
    xlabel(xlab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);
    ylabel(ylab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.11 0.5 0]);
    
    ax = gca;
    %     outerpos = ax.OuterPosition;
    %     ti = ax.TightInset;
    %     left = outerpos(1) + ti(1) + 0.02;
    %     bottom = outerpos(2) + ti(2) + 0.02;
    %     ax_width = outerpos(3) - ti(1) - ti(3) - 0.04;
    %     ax_height = outerpos(4) - ti(2) - ti(4) - 0.06;
    %     ax.Position = [left bottom ax_width ax_height];
    ax.Position = [0.1076 0.1200 0.8724 0.8524];
    
    print(gcf,'-djpeg',strcat(dirname,'sfr_norm_dist.jpg'));
    print(gcf,'-depsc',strcat(dirname,'sfr_norm_dist.eps'));
    close all
    fclose all;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(tsf_age)
    xl = 10;
    xu = 1000;
    xtick = [10 100 1000];
    xstr = {num2str(10,'%02i'),num2str(100,'%03i'),num2str(1000,'%04i')};
    yl = 1;
    yu = 1000;
    ytick = [1 10 100 1000 10000];
    ystr = {num2str(1,'%01i'),num2str(10,'%02i'),num2str(100,'%03i'),num2str(1000,'%04i'),num2str(10000,'%05i')};
    %ystr = {'10^0','10^1','10^2','10^3','10^4'};
    
    ylab  = '$t_{\rm SF}\:{\rm [Myr]}$';
    ylab2 = '$t_{\rm SF}\:\:{\rm [Myr\:]}$';
    xlab  = '${\rm age}\:{\rm [Myr]}$';
    xlab2 = '${\rm age}\:\:{\rm [Myr]}$';
    
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
        'Position',[0.13 0.14 0.775 0.815]);
    xlim(axes1,[xl xu]);
    ylim(axes1,[yl yu]);
    box(axes1,'on');
    % grid on
    hold(axes1,'all');
    xlabel(xlab,'Interpreter','latex','FontSize',22,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);
    ylabel(ylab,'Interpreter','latex','FontSize',22,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.11 0.5 0]);
    set(gcf,'renderer','painters')
    
    plot(age, tsf, 'linestyle','none','marker','o','markerfacecolor','k','markeredgecolor','k','markersize',3)
    
    xbin = 0.6:0.2:3.4;
    dbin = (xbin(2)-xbin(1))/2;
    Nbin = length(xbin);
    tsfbin = zeros(Nbin,3);
    for i=1:Nbin
        b = find(SFR_norm>=-9 & log10(age)>=xbin(i)-dbin & log10(age)<xbin(i)+dbin);
        temp = sort(tsf(b));
        Ntemp = length(temp);
        if(Ntemp>10)
            N1 = ceil(0.16*Ntemp);
            N2 = ceil(0.84*Ntemp);
            tsfbin(i,1) = median(tsf(b));
            tsfbin(i,2) = median(temp(N1));
            tsfbin(i,3) = median(temp(N2));
        end
    end
    l(1) = plot(10.^xbin, tsfbin(:,1), 'linestyle','-','marker','o',...
        'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8,...
        'color','k','linewidth',4,'DisplayName','median');
    plot(10.^xbin, tsfbin(:,2), 'linestyle','--','marker','none','color','k','linewidth',4);
    plot(10.^xbin, tsfbin(:,3), 'linestyle','--','marker','none','color','k','linewidth',4);
    
    l(2) = plot(time_model, tsf_model, 'linestyle','-','marker','none','color','c','linewidth',4,'DisplayName','model');
    legend1 = legend(gca,l(1:2));
    set(legend1,...
        'Location','NorthWest','FontSize',16,...
        'Interpreter','tex');
    legend boxoff
    
    
    annotation(figure1,'textbox',[0.75 0.15 0.10 0.05],...
        'String',strcat('${\rm \mu}\:\:=',num2str(mu_model,'%4.2f'),'$'),...
        'FontSize',14,...
        'FontName','Times',...
        'Interpreter','Latex',...
        'FitBoxToText','off',...
        'LineStyle','none',...
        'BackgroundColor',[1 1 1]);
    annotation(figure1,'textbox',[0.75 0.21 0.10 0.05],...
        'String',strcat('${\rm f_{g0}}=',num2str(fgi_model,'%4.2f'),'$'),...
        'FontSize',14,...
        'FontName','Times',...
        'Interpreter','Latex',...
        'FitBoxToText','off',...
        'LineStyle','none',...
        'BackgroundColor',[1 1 1]);
    annotation(figure1,'textbox',[0.45 0.15 0.10 0.05],...
        'String',strcat('${\rm \epsilon_{d}}=',num2str(eps_model,'%4.2f'),'$'),...
        'FontSize',16,...
        'FontName','Times',...
        'Interpreter','Latex',...
        'FitBoxToText','off',...
        'LineStyle','none',...
        'BackgroundColor',[1 1 1]);
    annotation(figure1,'textbox',[0.45 0.21 0.10 0.05],...
        'String',strcat('${\rm \alpha}\:=',num2str(alpha_model,'%4.2f'),'$'),...
        'FontSize',16,...
        'FontName','Times',...
        'Interpreter','Latex',...
        'FitBoxToText','off',...
        'LineStyle','none',...
        'BackgroundColor',[1 1 1]);
    annotation(figure1,'textbox',[0.26 0.15 0.10 0.05],...
        'String',strcat('${\rm \eta_{\rm s}}=',num2str(etas_model,'%4.2f'),'$'),...
        'FontSize',16,...
        'FontName','Times',...
        'Interpreter','Latex',...
        'FitBoxToText','off',...
        'LineStyle','none',...
        'BackgroundColor',[1 1 1]);
    annotation(figure1,'textbox',[0.26 0.21 0.10 0.05],...
        'String',strcat('${\rm \eta}\:=',num2str(eta_model,'%4.2f'),'$'),...
        'FontSize',16,...
        'FontName','Times',...
        'Interpreter','Latex',...
        'FitBoxToText','off',...
        'LineStyle','none',...
        'BackgroundColor',[1 1 1]);
    
    print(gcf,'-depsc',strcat(dirname,'tsf_age.eps'));
    set(gcf,'visible','off');
    xlabel(xlab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);
    ylabel(ylab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.11 0.5 0]);
    
    ax = gca;
    %     outerpos = ax.OuterPosition;
    %     ti = ax.TightInset;
    %     left = outerpos(1) + ti(1) + 0.02;
    %     bottom = outerpos(2) + ti(2) + 0.02;
    %     ax_width = outerpos(3) - ti(1) - ti(3) - 0.04;
    %     ax_height = outerpos(4) - ti(2) - ti(4) - 0.06;
    %     ax.Position = [left bottom ax_width ax_height];
    ax.Position = [0.1076 0.1200 0.8724 0.8524];
    
    print(gcf,'-djpeg',strcat(dirname,'tsf_age.jpg'));
    print(gcf,'-depsc',strcat(dirname,'tsf_age.eps'));
    close all
    fclose all;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(mass_age)
    xl = 10;
    xu = 1000;
    xtick = [10 100 1000];
    xstr = {num2str(10,'%02i'),num2str(100,'%03i'),num2str(1000,'%04i')};
    yl = 1e-4;
    yu = 1;
    ytick = [1e-4 1e-3 0.01 0.1 1];
    % ystr = {num2str(0.0001,'%6.4f'),num2str(0.001,'%5.3f'),num2str(0.01,'%4.2f'),num2str(0.1,'%3.1'),num2str(1,'%1.0f')};
    ystr = {'10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^0'};
    
    ylab  = '$M_*\:/\:M^{\rm gal}_*$';
    ylab2 = '$M_*\:\:/\:\:M^{\rm gal}_*$';
    xlab  = '${\rm age}\:{\rm [Myr]}$';
    xlab2 = '${\rm age}\:\:{\rm [Myr]}$';
    
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
        'Position',[0.13 0.14 0.775 0.815]);
    xlim(axes1,[xl xu]);
    ylim(axes1,[yl yu]);
    box(axes1,'on');
    % grid on
    hold(axes1,'all');
    xlabel(xlab,'Interpreter','latex','FontSize',22,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);
    ylabel(ylab,'Interpreter','latex','FontSize',22,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.11 0.5 0]);
    set(gcf,'renderer','painters')
    
    b=find(Mclump_norm>=-3.6);
    plot(age(b), 10.^Mclump_norm(b), 'linestyle','none','marker','o','markerfacecolor','k','markeredgecolor','k','markersize',3)
    
    xbin = 0.6:0.2:3.4;
    dbin = (xbin(2)-xbin(1))/2;
    Nbin = length(xbin);
    Mclump_norm_bin = zeros(Nbin,3);
    for i=1:Nbin
        b = find(Mclump_norm>=-9 & log10(age)>=xbin(i)-dbin & log10(age)<xbin(i)+dbin);
        temp = sort(Mclump_norm(b));
        Ntemp = length(temp);
        if(Ntemp>10)
            N1 = ceil(0.16*Ntemp);
            N2 = ceil(0.84*Ntemp);
            Mclump_norm_bin(i,1) = median(Mclump_norm(b));
            Mclump_norm_bin(i,2) = median(temp(N1));
            Mclump_norm_bin(i,3) = median(temp(N2));
        end
    end
    l(1) = plot(10.^xbin, 10.^Mclump_norm_bin(:,1), 'linestyle','-','marker','o',...
        'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8,...
        'color','k','linewidth',4,'DisplayName','median');
    plot(10.^xbin, 10.^Mclump_norm_bin(:,2), 'linestyle','--','marker','none','color','k','linewidth',4);
    plot(10.^xbin, 10.^Mclump_norm_bin(:,3), 'linestyle','--','marker','none','color','k','linewidth',4);
    
    l(2) = plot(time_model, 10.^Ms_model, 'linestyle','-','marker','none','color','c','linewidth',4,'DisplayName','model');
    legend1 = legend(gca,l(1:2));
    set(legend1,...
        'Location','NorthWest','FontSize',16,...
        'Interpreter','tex');
    legend boxoff
    
    print(gcf,'-depsc',strcat(dirname,'Mstar_norm_age.eps'));
    set(gcf,'visible','off');
    xlabel(xlab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);
    ylabel(ylab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.11 0.5 0]);
    
    ax = gca;
    %     outerpos = ax.OuterPosition;
    %     ti = ax.TightInset;
    %     left = outerpos(1) + ti(1) + 0.02;
    %     bottom = outerpos(2) + ti(2) + 0.02;
    %     ax_width = outerpos(3) - ti(1) - ti(3) - 0.04;
    %     ax_height = outerpos(4) - ti(2) - ti(4) - 0.06;
    %     ax.Position = [left bottom ax_width ax_height];
    ax.Position = [0.1076 0.1200 0.8724 0.8524];
    
    print(gcf,'-djpeg',strcat(dirname,'Mstar_norm_age.jpg'));
    print(gcf,'-depsc',strcat(dirname,'Mstar_norm_age.eps'));
    close all
    fclose all;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(sfr_age)
    xl = 10;
    xu = 1000;
    xtick = [10 100 1000];
    xstr = {num2str(10,'%02i'),num2str(100,'%03i'),num2str(1000,'%04i')};
    yl = 1e-3;
    yu = 10;
    ytick = [1e-3 0.01 0.1 1 10];
    % ystr = {num2str(0.0001,'%6.4f'),num2str(0.001,'%5.3f'),num2str(0.01,'%4.2f'),num2str(0.1,'%3.1'),num2str(1,'%1.0f')};
    ystr = {'10^{-3}','10^{-2}','10^{-1}','10^0','10^1'};
    
    ylab  = '${\rm SFR}\:/\:{\rm SFR^{gal}}$';
    ylab2 = '${\rm SFR}\:\:/\:\:{\rm SFR^{gal}}$';
    xlab  = '${\rm age}\:{\rm [Myr]}$';
    xlab2 = '${\rm age}\:\:{\rm [Myr]}$';
    
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
        'Position',[0.13 0.14 0.775 0.815]);
    xlim(axes1,[xl xu]);
    ylim(axes1,[yl yu]);
    box(axes1,'on');
    % grid on
    hold(axes1,'all');
    xlabel(xlab,'Interpreter','latex','FontSize',22,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);
    ylabel(ylab,'Interpreter','latex','FontSize',22,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.11 0.5 0]);
    set(gcf,'renderer','painters')
    
    b=find(SFR_norm>=-3.6);
    plot(age(b), 10.^SFR_norm(b), 'linestyle','none','marker','o','markerfacecolor','k','markeredgecolor','k','markersize',3)
    
    xbin = 0.6:0.2:3.4;
    dbin = (xbin(2)-xbin(1))/2;
    Nbin = length(xbin);
    SFR_norm_bin = zeros(Nbin,3);
    for i=1:Nbin
        b = find(SFR_norm>=-9 & log10(age)>=xbin(i)-dbin & log10(age)<xbin(i)+dbin);
        temp = sort(SFR_norm(b));
        Ntemp = length(temp);
        if(Ntemp>10)
            N1 = ceil(0.16*Ntemp);
            N2 = ceil(0.84*Ntemp);
            SFR_norm_bin(i,1) = median(SFR_norm(b));
            SFR_norm_bin(i,2) = median(temp(N1));
            SFR_norm_bin(i,3) = median(temp(N2));
        end
    end
    l(1) = plot(10.^xbin, 10.^SFR_norm_bin(:,1), 'linestyle','-','marker','o',...
        'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8,...
        'color','k','linewidth',4,'DisplayName','median');
    plot(10.^xbin, 10.^SFR_norm_bin(:,2), 'linestyle','--','marker','none','color','k','linewidth',4);
    plot(10.^xbin, 10.^SFR_norm_bin(:,3), 'linestyle','--','marker','none','color','k','linewidth',4);
    l(2) = plot(time_model, 10.^SFR_norm_model, 'linestyle','-','marker','none','color','c','linewidth',4,'DisplayName','model');
    legend1 = legend(gca,l(1:2));
    set(legend1,...
        'Position',[0.20 0.12 0.15 0.15],'FontSize',16,...
        'Interpreter','tex');
    legend boxoff
    
    print(gcf,'-depsc',strcat(dirname,'sfr_norm_age.eps'));
    set(gcf,'visible','off');
    xlabel(xlab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 -0.07 0]);
    ylabel(ylab2,'Interpreter','latex','FontSize',20,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.11 0.5 0]);
    
    ax = gca;
    %     outerpos = ax.OuterPosition;
    %     ti = ax.TightInset;
    %     left = outerpos(1) + ti(1) + 0.02;
    %     bottom = outerpos(2) + ti(2) + 0.02;
    %     ax_width = outerpos(3) - ti(1) - ti(3) - 0.04;
    %     ax_height = outerpos(4) - ti(2) - ti(4) - 0.06;
    %     ax.Position = [left bottom ax_width ax_height];
    ax.Position = [0.1076 0.1200 0.8724 0.8524];
    
    print(gcf,'-djpeg',strcat(dirname,'sfr_norm_age.jpg'));
    print(gcf,'-depsc',strcat(dirname,'sfr_norm_age.eps'));
    close all
    fclose all;
end