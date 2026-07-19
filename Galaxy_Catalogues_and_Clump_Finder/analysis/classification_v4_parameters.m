function classification_v4_parameters(bulge, insitu, exsitu, norm_is, norm_es, ...
    norm_bulge, norm, is_only, tit, type)
if(is_only == 0)
    redshift = [bulge(:,2)',exsitu(:,2)',insitu(:,2)']';
    rad = log10(1000.*[bulge(:,3)',exsitu(:,3)',insitu(:,3)']');
    mgas = log10([bulge(:,4)',exsitu(:,4)',insitu(:,4)']');
    norm_mgas = log10([norm_bulge(:,4)',norm_es(:,4)',norm_is(:,4)']');
    mstar = log10([bulge(:,5)',exsitu(:,5)',insitu(:,5)']');
    norm_mstar = log10([norm_bulge(:,5)',norm_es(:,5)',norm_is(:,5)']');
    mass = log10([bulge(:,6)',exsitu(:,6)',insitu(:,6)']');
    norm_mass = log10([norm_bulge(:,6)',norm_es(:,6)',norm_is(:,6)']');
    fgas = [bulge(:,7)',exsitu(:,7)',insitu(:,7)']';
    norm_fgas = [norm_bulge(:,7)',norm_es(:,7)',norm_is(:,7)']';
    Sig_g = log10([bulge(:,9)',exsitu(:,9)',insitu(:,9)']');
    Sig_s = log10([bulge(:,10)',exsitu(:,10)',insitu(:,10)']');
    norm_Sig_s = log10([norm_bulge(:,10)',norm_es(:,10)',norm_is(:,10)']');
    age = log10([bulge(:,12)',exsitu(:,12)',insitu(:,12)']');
    norm_age = log10([norm_bulge(:,12)',norm_es(:,12)',norm_is(:,12)']');
    zgas = [bulge(:,13)',exsitu(:,13)',insitu(:,13)']';
    norm_zgas = [norm_bulge(:,13)',norm_es(:,13)',norm_is(:,13)']';
    zstars = [bulge(:,14)',exsitu(:,14)',insitu(:,14)']';
    SFR = [bulge(:,15)',exsitu(:,15)',insitu(:,15)']';
    norm_SFR = [norm_bulge(:,15)',norm_es(:,15)',norm_is(:,15)']';
    Sig_SFR = [bulge(:,16)',exsitu(:,16)',insitu(:,16)']';
    norm_Sig_SFR = [norm_bulge(:,16)',norm_es(:,16)',norm_is(:,16)']';
    sSFR = [bulge(:,17)',exsitu(:,17)',insitu(:,17)']';
    norm_sSFR = [norm_bulge(:,17)',norm_es(:,17)',norm_is(:,17)']';
    tdep = [bulge(:,18)',exsitu(:,18)',insitu(:,18)']';
    norm_dist = [bulge(:,19)',exsitu(:,19)',insitu(:,19)']';
    norm_height = abs([bulge(:,20)',exsitu(:,20)',insitu(:,20)'])';
    dist = [bulge(:,21)',exsitu(:,21)',insitu(:,21)']';
    height = abs([bulge(:,22)',exsitu(:,22)',insitu(:,22)'])';
    del = log10([bulge(:,23)',exsitu(:,23)',insitu(:,23)']');    
    eta = [bulge(:,24)',exsitu(:,24)',insitu(:,24)']';
    del_DM = log10([bulge(:,27)',exsitu(:,27)',insitu(:,27)']');
    tff = [bulge(:,29)',exsitu(:,29)',insitu(:,29)']';
    tdyn = [bulge(:,30)',exsitu(:,30)',insitu(:,30)']';
    tdyn_g = [bulge(:,31)',exsitu(:,31)',insitu(:,31)']';
elseif(is_only == 1)
    redshift = insitu(:,2);
    rad = log10(1000.*insitu(:,3));
    mgas = log10(insitu(:,4));
    norm_mgas = log10(norm_is(:,4));
    mstar = log10(insitu(:,5));
    norm_mstar = log10(norm_is(:,5));
    mass = log10(insitu(:,6));
    norm_mass = log10(norm_is(:,6));
    fgas = insitu(:,7);
    norm_fgas = norm_is(:,7);
    Sig_g = log10(insitu(:,9));
    Sig_s = log10(insitu(:,10));    
    norm_Sig_s = log10(norm_is(:,10));
    age = log10(insitu(:,12));
    norm_age = log10(norm_is(:,12));
    zgas = insitu(:,13);
    norm_zgas = norm_is(:,13);
    zstars = insitu(:,14);
    SFR = insitu(:,15);
    norm_SFR = norm_is(:,15);
    Sig_SFR = insitu(:,16);
    norm_Sig_SFR = norm_is(:,16);
    sSFR = insitu(:,17);
    norm_sSFR = norm_is(:,17);
    tdep = insitu(:,18);
    norm_dist = insitu(:,19);
    norm_height = abs(insitu(:,20));
    dist = insitu(:,21);
    height = abs(insitu(:,22));
    del = log10(insitu(:,23));
    eta = insitu(:,24);
    del_DM = log10(insitu(:,27));
    tff = insitu(:,29);
    tdyn = insitu(:,30);
    tdyn_g = insitu(:,31);
end
Nc = length(eta);

b=find(fgas<=10^(-1.3));
fgas(b)=10^(-1.3);
b=find(norm_fgas<=10^(-1.9));
norm_fgas(b)=10^(-1.9);

b=find(SFR<=1e-4);
SFR(b)=10^(-3.8);
norm_SFR(b)=10^(-3.8);

b=find(Sig_SFR<=10^(-3.8));
Sig_SFR(b)=10^(-3.8);
norm_Sig_SFR(b)=10^(-3.8);

b=find(sSFR<=1e-2);
sSFR(b)=10^(-1.8);
norm_sSFR(b)=10^(-1.8);

b=find(tdep<=1e-2);
tdep(b)=(1e-2);
tdep = log10(1./tdep);

b=find(isinf(del_DM));
del_DM(b)=-0.9;
b=find(isnan(del_DM));
del_DM(b)=-0.9;

b=find(isinf(mstar) | mstar<=3.1);
mstar(b)=3.1;
sSFR(b)=10^(-1.8);
norm_sSFR(b)=10^(-1.8);
age(b) = 1.5;
norm_age(b)=-0.9;
b=find(isinf(Sig_s) | Sig_s<=0.1);
Sig_s(b) = 0.1;
norm_Sig_s(b) = 10^(-3.8);

b=find( isinf(mgas) | mgas<=3.1 );
mgas(b)=3.1;
b=find(isinf(Sig_g) | Sig_g<=0.1);
Sig_g(b) = 0.1;
b=find(isinf(zgas) | zgas<=7.85);
zgas(b) = 7.85;
b=find(isinf(norm_zgas) | norm_zgas<=10^(-0.45));
norm_zgas(b) = 10^(-0.45);

td_ff = log10(abs(tdyn)./tff);
tdg_ff = log10(abs(tdyn_g)./tff);

fgas=log10(fgas);
norm_fgas=log10(norm_fgas);

SFR=log10(SFR);
norm_SFR=log10(norm_SFR);

Sig_SFR=log10(Sig_SFR);
norm_Sig_SFR=log10(norm_Sig_SFR);

sSFR=log10(sSFR);
norm_sSFR=log10(norm_sSFR);

norm_zgas=log10(norm_zgas);

Rd = 1000.*dist./norm_dist;
norm_rad = rad - log10(Rd);

if(norm == 0)
    x = mass;
    binx = 5.8:0.03:11.2;
    xt = 6:0.5:11;
    y = rad;
    biny = 1:0.03:3.6;
    yt = 1.6:0.2:3.4;
else
    binx = -4.8:0.03:1.2;
    x = norm_mass;
    xt = -4:0.5:0.5;
    y = norm_rad;
    biny = -4.8:0.03:1.2;
    yt = -2.5:0.5:-0.5;
end
nx = length(binx);
ny = length(biny);
n = -50.*ones(ny,nx);
num = zeros(ny,nx);
for i=1:nx-1
    for j=1:ny-1
        b = find(x>=binx(i) & x<binx(i+1) & y>=biny(j) & y<biny(j+1));
        if(~isempty(b))
            num(j,i) = length(b)/Nc;
            if(type==1)
                n(j,i) = log10(length(b)/Nc);
            elseif(type==2)
                n(j,i) = median(redshift(b));
            elseif(type==3)
                n(j,i) = median(del_DM(b));
            elseif(type==4)
                n(j,i) = median(del(b));
            elseif(type==5)
                n(j,i) = median(eta(b));            
            elseif(type==6)
                n(j,i) = median(fgas(b));
            elseif(type==7)
                n(j,i) = median(dist(b));
            elseif(type==8)
                n(j,i) = median(height(b));
            elseif(type==9)
                n(j,i) = median(norm_dist(b));
            elseif(type==10)
                n(j,i) = median(norm_height(b));
            elseif(type==11)
                n(j,i) = median(SFR(b));
            elseif(type==12)
                n(j,i) = median(mstar(b));
            elseif(type==14)
                n(j,i) = median(mgas(b));
            elseif(type==15)
                n(j,i) = median(age(b));
            elseif(type==16)
                n(j,i) = median(sSFR(b));
            elseif(type==17)
                n(j,i) = median(zgas(b));
            elseif(type==18)
                n(j,i) = median(zstars(b));
            elseif(type==19)
                n(j,i) = median(Sig_g(b));
            elseif(type==20)
                n(j,i) = median(Sig_s(b));
            elseif(type==21)
                n(j,i) = median(tdep(b));
            elseif(type==22)
                n(j,i) = median(td_ff(b));
            elseif(type==23)
                n(j,i) = median(tdg_ff(b));
            elseif(type==25)
                n(j,i) = median(norm_SFR(b));
            elseif(type==26)
                n(j,i) = median(max_norm_dist(b));
            elseif(type==27)
                n(j,i) = median(max_norm_height(b));
            elseif(type==28)
                n(j,i) = median(Sig_SFR(b));
            elseif(type==29)
                n(j,i) = median(norm_Sig_SFR(b));
            elseif(type==30)
                n(j,i) = median(norm_Sig_s(b));
            elseif(type==32)
                n(j,i) = median(norm_fgas(b));
            elseif(type==33)
                n(j,i) = median(norm_sSFR(b));
            elseif(type==34)
                n(j,i) = median(norm_zgas(b));
            elseif(type==35)
                n(j,i) = median(norm_age(b));
            end
        end
    end
end
if(type==1)
    titc = '${\rm log(N\:/\:N_{tot})}$';
    cl = [-6 -2 0.4];
elseif(type==2)
    titc = '$z$';
    cl = [0 10 1];
elseif(type==3)
    titc = '${\rm log(1 + \delta_{DM})}$';
    cl = [-1 1 0.2];
elseif(type==4)
    titc = '${\rm log({\bar{\delta}})}$';
    cl = [0.9 2.3 0.2];
elseif(type==5)
    titc = '${\rm S_{\rm c}}$';
    cl = [0 1 0.1];
elseif(type==6)
    titc = '${\rm log(}\:\:f_{\rm gas}\:{\rm )}$';
    cl = [-1.4 0 0.1];
elseif(type==7)
	titc = '$dist\:{\rm [kpc]}$';
    cl = [0 10 1];
elseif(type==8)
	titc = '$h\:{\rm [kpc]}$';
    cl = [0 7 1];
elseif(type==9)
    titc = '$d\:/\:R_{\rm d}$';
    cl = [0 2 0.2];
elseif(type==10)
    titc = '$h\:/\:H_{\rm d}$';
    cl = [0 3.5 0.5];
elseif(type==11)
    titc = '${\rm log(}\:\:SFR\:{\rm [\:\:M_{\odot}\:\:yr^{-1}])}$';
    cl = [-4 1 0.5];
elseif(type==12)
    titc = '${\rm log(}\:\:M_{\rm *}\:{\rm [M_{\odot}])}$';
    cl = [3 9 0.5];
elseif(type==14)
    titc = '${\rm log(}\:\:M_{\rm gas}\:{\rm [M_{\odot}])}$';
    cl = [3 9 0.5];
elseif(type==15)
    titc = '${\rm log(}\:\:age_{\rm *}\:{\rm [Myr])}$';
    cl = [1.4 3.4 0.2];
elseif(type==16)
    titc = '${\rm log(}\:\:sSFR\:{\rm [Gyr^{-1}\:])}$';
    cl = [-1.6 1.6 0.2];
elseif(type==17)
    titc = '$Z_{\rm gas}\:\:{\rm [log(O/H)+12]}$';
    cl = [7.8 9 0.1];
elseif(type==18)
    titc = '$Z_{\rm stars}\:\:{\rm [log(O/H)+12]}$';
    cl = [7.8 9 0.1];
elseif(type==19)
    titc = '${\rm log(}\:\:\Sigma_{\rm gas}\:\:{\rm [M_{\odot}\:pc^{-2}\:])}$';
    cl = [0 3 0.5];
elseif(type==20)
    titc = '${\rm log(}\:\:\Sigma_{\rm stars}\:\:{\rm [M_{\odot}\:pc^{-2}\:])}$';
    cl = [0 3 0.5];
elseif(type==21)
    titc = '${\rm log(}\:\:t_{\rm dep}\:{\rm [Gyr\:])}$';
    cl = [-2 2 0.5];
elseif(type==22)
    titc = '${\rm log(}\:\:t_{\rm d}\:/\:t_{\rm ff}{\rm )}$';
    cl = [-1 1 0.2];
elseif(type==23)
    titc = '${\rm log(}\:\:t_{\rm d,\:global}\:\:\:\:/\:t_{\rm ff}{\rm )}$';
    cl = [-1 1 0.2];
elseif(type==25)
    titc = '${\rm log(}\:\:SFR_{\rm c}\:/\:SFR_{\rm d}{\rm )}$';
    cl = [-4 0 0.5];
elseif(type==26)
    titc = '${\rm max(}\:\:\:d\:/\:R_{\rm d}\:\:{\rm )}$';
    cl = [0 2 0.2];
elseif(type==27)
    titc = '${\rm max(}\:\:\:h\:/\:H_{\rm d}\:\:{\rm )}$';
    cl = [0 3.5 0.5];
elseif(type==28)
    titc = '${\rm log(}\:\:\Sigma_{\rm SFR}\:\:{\rm [\:M_{\odot}\:\:yr^{-1}\:\:kpc^{-2}\:\:])}$';
    cl = [-4 1 0.5];
elseif(type==29)
    titc = '${\rm log(}\:\:\Sigma_{\rm SFR,\:c}\:\:\:/\:\Sigma_{\rm SFR,\:d}\:\:{\rm )}$';
    cl = [-4 2 0.5];
elseif(type==30)
    titc = '${\rm log(}\:\:\Sigma_{\rm *,\:c}\:\:\:/\:\Sigma_{\rm *,\:d}\:\:{\rm )}$';
    cl = [-4 3 1];
elseif(type==32)
    titc = '${\rm log(}\:\:f_{\rm gas,c}\:/\:f_{\rm gas,d}\:\:{\rm )}$';
    cl = [-2 0.8 0.4];
elseif(type==33)
    titc = '${\rm log(}\:\:sSFR_{\rm c}\:/\:sSFR_{\rm d}\:\:{\rm )}$';
    cl = [-1.6 1.6 0.4];
elseif(type==34)
    titc = '${\rm log(}\:\:Z_{\rm gas,c}\:/\:Z_{\rm gas,d}\:\:{\rm )}$';
    cl = [-0.5 0.5 0.1];
elseif(type==35)
    titc = '${\rm log(}\:\:age_{\rm *,c}\:/\:age_{\rm *,d}\:\:{\rm )}$';
    cl = [-1 0.6 0.2];
end
            
% Create figure
figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68],...
    'Colormap',[0 0 0;0 0 0.25;0 0 0.5;0 0 0.75;0 0 1;0 0.07143 0.9929;0 0.1429 0.9857;0 0.2143 0.9786;0 0.2857 0.9714;0 0.3571 0.9643;0 0.4286 0.9571;0 0.5 0.95;0 0.5714 0.9429;0 0.6429 0.9357;0 0.7143 0.9286;0 0.7857 0.9214;0 0.8571 0.9143;0 0.9286 0.9071;0 1 0.9;0 1 0.825;0 1 0.75;0 1 0.675;0 1 0.6;0 1 0.525;0 1 0.45;0 1 0.375;0 1 0.3;0 1 0.225;0 1 0.15;0 1 0.075;0 1 0;0.15 1 0;0.3 1 0;0.45 1 0;0.6 1 0;0.75 1 0;0.9 1 0;0.92 1 0;0.94 1 0;0.96 1 0;0.98 1 0;1 1 0;1 0.9352 0;1 0.8704 0;1 0.8056 0;1 0.7407 0;1 0.6759 0;1 0.6111 0;1 0.584 0;1 0.5568 0;1 0.5296 0;1 0.5025 0;1 0.4753 0;1 0.3565 0;1 0.2377 0;1 0.1188 0;1 0 0;1 0.2 0.2;1 0.4 0.4;1 0.6 0.6;1 0.8 0.8;1 0.8667 0.8667;1 0.9333 0.9333;1 1 1]);
set(figure1,'WindowStyle','docked','Visible','off');

% Create axes
axes1 = axes('Parent',figure1,...
    'YTick',yt,...
    'XTick',xt,...
    'Position',[0.12 0.15 0.8 0.8],...
    'PlotBoxAspectRatio',[1 1 1],...
    'FontSize',16,...
    'FontName','Arial',...
    'CLim',[cl(1) cl(2)]);
% Uncomment the following line to preserve the X-limits of the axes
 xlim([xt(1) xt(length(xt))]);
% Uncomment the following line to preserve the Y-limits of the axes
 ylim([yt(1) yt(length(yt))]);
grid('off');
hold('all');

% Create surf
surf(binx,biny,n,'Parent',axes1,'linestyle','none');

% Create xlabel
if(norm==0)
    titx = '${\rm Log(}\:\:\:M_{\rm c}\:{\rm [M_{\odot}])}$';
else
    titx = '${\rm Log(}\:\:\:M_{\rm c}\:/\:M_{\rm d}{\rm )}$';
end
xlabel(titx,'Interpreter','latex',...
    'units','normalized','position',[0.5 -0.08 0],...
    'FontSize',18,'FontName','Times New Roman');

% Create ylabel
if(norm==0)
    tity = '${\rm Log(}\:\:\:R_{\rm c}\:{\rm [pc])}$';
else
    tity = '${\rm Log(}\:\:\:R_{\rm c}\:/\:R_{\rm d}{\rm )}$';
end
ylabel(tity,'Interpreter','latex','FontSize',18,...
    'units','normalized','position',[-0.09 0.5 0],...
    'FontName','Times New Roman');

% Create title
if(~strcmp(tit,''))
%     title(tit,...
%         'Interpreter','latex','FontSize',12,'units','normalized','position',[0.5 0.915 0],...
%         'FontName','Times New Roman','BackGroundColor','w');
    annotation(figure1,'textbox',[0.27 0.896 0.40 0.05],...
        'String',tit,...
        'Interpreter','latex',...
        'FontSize',12,...
        'FontName','Times',...
        'FitBoxToText','off',...
        'LineStyle','none',...
        'BackgroundColor',[1 1 1]);
end

% Create colorbar
bar=colorbar('peer',axes1,'FontSize',16,'FontName','Arial',...
    'CLim',[1 64],'YTick',[cl(1):cl(3):cl(2)]);
set(get(bar,'Ylabel'),'String',titc,...
    'Interpreter','Latex','Fontsize',18,'Rotation',270,'position',[10.5,(cl(1)+cl(2))/2,9.16],'FontName','Times');

if(norm==0)
    % Create plot
    x1 = (xt(1)-1):0.5:(xt(length(xt))+1);
    y1 = (1/3).*x1 - (1/3).*log10(0.0003363 * 4*pi/3);
    y2 = (1/3).*x1 - (1/3).*log10(0.003363 * 4*pi/3);
    y3 = (1/3).*x1 - (1/3).*log10(0.033630 * 4*pi/3);
    y4 = (1/3).*x1 - (1/3).*log10(0.336300 * 4*pi/3);
    y5 = (1/3).*x1 - (1/3).*log10(3.363000 * 4*pi/3);
    y6 = (1/3).*x1 - (1/3).*log10(33.63000 * 4*pi/3);

    plot3(x1,y1,50.*x1./x1,'Parent',axes1,'LineWidth',3,'Color',[1 0 0]);
    plot3(x1,y2,50.*x1./x1,'Parent',axes1,'LineWidth',3,'Color',[1 0 0]);
    plot3(x1,y3,50.*x1./x1,'Parent',axes1,'LineWidth',3,'Color',[1 0 0]);
    plot3(x1,y4,50.*x1./x1,'Parent',axes1,'LineWidth',3,'Color',[1 0 0]);
    plot3(x1,y5,50.*x1./x1,'Parent',axes1,'LineWidth',3,'Color',[1 0 0]);
    plot3(x1,y6,50.*x1./x1,'Parent',axes1,'LineWidth',3,'Color',[1 0 0]);

%     % Create plot
%     y1 = (1/2).*x1 - (1/2).*log10(1.00 * pi);
%     y2 = (1/2).*x1 - (1/2).*log10(10.0 * pi);
%     y3 = (1/2).*x1 - (1/2).*log10(100. * pi);
%     y4 = (1/2).*x1 - (1/2).*log10(1000 * pi);
% 
%     plot3(x1,y1,50.*x1./x1,'Parent',axes1,'LineWidth',3,'Color','c');
%     plot3(x1,y2,50.*x1./x1,'Parent',axes1,'LineWidth',3,'Color','c');
%     plot3(x1,y3,50.*x1./x1,'Parent',axes1,'LineWidth',3,'Color','c');
%     plot3(x1,y4,50.*x1./x1,'Parent',axes1,'LineWidth',3,'Color','c');
elseif(norm==1)
    % Create plot
    x1 = (xt(1)-1):0.49:(xt(length(xt))+10);
    y1 = (1/3).*x1 - (1/3).*log10(1);
    y2 = (1/3).*x1 - (1/3).*log10(10);
    y3 = (1/3).*x1 - (1/3).*log10(100);
    y4 = (1/3).*x1 - (1/3).*log10(1000);

    plot3(x1,y1,50.*x1./x1,'Parent',axes1,'LineWidth',3,'Color','r');
    plot3(x1,y2,50.*x1./x1,'Parent',axes1,'LineWidth',3,'Color','r');
    plot3(x1,y3,50.*x1./x1,'Parent',axes1,'LineWidth',3,'Color','r');
    plot3(x1,y4,50.*x1./x1,'Parent',axes1,'LineWidth',3,'Color','r');

%     y1 = (1/2).*x1 - (1/2).*log10(0.1);
%     y2 = (1/2).*x1 - (1/2).*log10(1);
%     y3 = (1/2).*x1 - (1/2).*log10(10);
%     y4 = (1/2).*x1 - (1/2).*log10(100);
%     
%     plot3(x1,y1,50.*x1./x1,'Parent',axes1,'LineWidth',3,'Color','c');
%     plot3(x1,y2,50.*x1./x1,'Parent',axes1,'LineWidth',3,'Color','c');
%     plot3(x1,y3,50.*x1./x1,'Parent',axes1,'LineWidth',3,'Color','c');
%     plot3(x1,y4,50.*x1./x1,'Parent',axes1,'LineWidth',3,'Color','c');
end

for i=2:(ny-1)
    bx = find(num(i,2:nx-1)>0);
    if(~isempty(bx))
        num(i,bx) = ( num(i+1,bx+1) + num(i+1,bx) + num(i+1,bx-1) + ...
            num(i,bx+1) + num(i,bx) + num(i,bx-1) + ...
            num(i-1,bx+1) + num(i-1,bx) + num(i-1,bx-1) ) / 9;
    end
    clear bx
end
bx = find(num(1,2:nx-1)>0);
num(1,bx) = ( num(2,bx+1) + num(2,bx) + num(2,bx-1) + num(1,bx+1) + num(1,bx) + num(1,bx-1) ) / 6;
bx = find(num(ny,2:nx-1)>0);
num(ny,bx) = ( num(ny,bx+1) + num(ny,bx) + num(ny,bx-1) + num(ny-1,bx+1) + num(ny-1,bx) + num(ny-1,bx-1) ) / 6;
by = find(num(2:ny-1,1)>0);
num(by,1) = ( num(by+1,2) + num(by,2) + num(by-1,2) + num(by+1,1) + num(by,1) + num(by-1,1) ) / 6;
by = find(num(2:ny-1,nx)>0);
num(by,nx) = ( num(by+1,nx) + num(by,nx) + num(by-1,nx) + num(by+1,nx-1) + num(by,nx-1) + num(by-1,nx-1) ) / 6;

test = (max(max(num))+0.000001):-0.000001:(min(min(num))-0.000001);
test2 = 0.*test;
for i=1:length(test)
    b=find(num>=test(i));
    test2(i) = sum(num(b));
end
test2 = test2-0.5;
b = find(abs(test2)==min(abs(test2)));
v = [100+test(b(1)) 100+test(b(1))];
[c3,h] = contour3(binx,biny,100+num,v);
set(h,'linewidth',3,'edgecolor','w');
set(gcf,'renderer','painters')
