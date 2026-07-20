function [Lambda, Re_Omega, Im_Omega, Qb, Qs, Pmax] = P1_cyl(M, Del, m_mode, n_mode, Lam)

Mc = sqrt(Del)*M;
Rho_out = 1;
Rho_in = Del*Rho_out;
P_out = 1;
P_in = 1;
G_out = 5/3;
G_in = 5/3;
C_out = sqrt(G_out*P_out/Rho_out);
C_in = sqrt(G_in*P_in/Rho_in);
V_out = 0;
V_in = M * C_out;

Rs = 0.00625;

Length = 250001;
dirname = strcat('./M_',num2str(M,'%3.1f'),'_del_',num2str(Del,'%3i'),'_modes/data/');
filelist1 = dir(dirname);
filelist = filelist1(3:length(filelist1));
n = length(filelist);
X = zeros(1,Length);
ReP = zeros(1,Length);
ImP = zeros(1,Length);
for j=1:n
    n1=length(filelist(j).name);
    m2 = str2num(filelist(j).name(n1-4));
    if(filelist(j).name(n1-8)=='n')
        n2 = str2num(filelist(j).name(n1-7));
    else
        n2 = str2num(filelist(j).name(n1-8:n1-7));
    end
    if(n2==n_mode & m2==m_mode)
        if(filelist(j).name(1)=='x')
            a1 = load(strcat(dirname,filelist(j).name));
            Na1 = min(Length,length(a1));
            X(1:Na1) = abs(a1(1:Na1));
        elseif(filelist(j).name(1)=='p' & filelist(j).name(5)=='r')
            a1 = load(strcat(dirname,filelist(j).name));
            Na1 = min(Length,length(a1));
            ReP(1:Na1) = abs(a1(1:Na1));
        elseif(filelist(j).name(1)=='p' & filelist(j).name(5)=='i')
            a1 = load(strcat(dirname,filelist(j).name));
            Na1 = min(Length,length(a1));
            ImP(1:Na1) = abs(a1(1:Na1));
        end
    end
end
clear a1 Na1

if(Lam==0)  % input 0 wavelength means look for maximum growth rate of mode
    Om = X.*ImP;
    max_mode = max(Om);
    b = find(Om==max_mode);
    if(length(b)>1)
        b = b(1);
    end
    K = X(b)/Rs;
    Re_Omega = ReP(b)*K*V_in;
    Im_Omega = ImP(b)*K*V_in;
else
    wavelength = 2*pi./(X);  % in units of a
    min_err = min(abs(wavelength-Lam));
    b = find(abs(wavelength-Lam)==min_err);
    if(length(b)>1)
        b = b(1);
    end
    K = X(b)/Rs;
    Re_Omega = ReP(b)*K*V_in;
    Im_Omega = ImP(b)*K*V_in;
end
Omega = Re_Omega + 1i*Im_Omega;

Qb = sqrt(K^2 - (Omega/C_out - K*V_out/C_out)^2);
Qs = sqrt(K^2 - (Omega/C_in  - K*V_in/C_in)^2);

Lambda = 2*pi/K;
A = 0.05;
if(m_mode==0)
    H_up = -1.0 .* A .* Qb .* ( besseli(0,Qs.*Rs)./besselk(0,Qb.*Rs) ) .* ...
        besselk(1,Qb.*Rs) ./ ( Omega^2 * Rho_out );
else
    H_up = -1.0 .* A .* Qb .* ( besseli(m_mode,Qs.*Rs)./besselk(m_mode,Qb.*Rs) ) .* ...
        0.5.*(besselk(m_mode+1,Qb.*Rs)+besselk(m_mode-1,Qb.*Rs)) ./ ( Omega^2 * Rho_out );
end
abs(H_up) / Lambda
D = exp(1i.*m_mode.*pi).*A;

Ngrid = 401;
cmap = [0 0 0;0.05602 0.05672 0.1007;0.112 0.1134 0.2014;0.1681 0.1702 0.302;...
    0.2241 0.2269 0.4027;0.1821 0.1844 0.5147;0.1401 0.1418 0.6267;...
    0.09804 0.09927 0.7387;0.05602 0.05672 0.8507;0.04202 0.04254 0.888;...
    0.02801 0.02836 0.9254;0.014 0.01418 0.9627;0 0 1;0 0.1143 1;0 0.2286 1;...
    0 0.3429 1;0 0.4571 1;0 0.5714 1;0 0.6857 1;0 0.8 1;0 0.84 1;0 0.88 1;...
    0 0.92 1;0 0.96 1;0 1 1;0 1 0.8333;0 1 0.6667;0 1 0.5;0 1 0.3333;0 1 0.1667;...
    0 1 0;0.125 1 0;0.25 1 0;0.375 1 0;0.5 1 0;0.625 1 0;0.75 1 0;0.875 1 0;...
    1 1 0;1 0.9688 0;1 0.9375 0;1 0.9063 0;1 0.875 0;1 0.7969 0;1 0.7188 0;...
    1 0.6406 0;1 0.5625 0;1 0.5313 0;1 0.5 0;1 0.4688 0;1 0.4375 0;1 0.3646 0;...
    1 0.2917 0;1 0.2188 0;1 0.1458 0;1 0.09722 0;1 0.04861 0;1 0 0;0.9167 0 0;...
    0.8333 0 0;0.75 0 0;0.6667 0 0;0.5833 0 0;0.5 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xz plane
[x,y] = meshgrid(linspace(0,10,Ngrid),linspace(-5,5,Ngrid));
P1 = zeros(Ngrid,Ngrid);
r = abs(y);

b = find(y>=1);
P1(b) = A .* ( besselk(m_mode,Qb.*r(b).*Rs) ./ besselk(m_mode,Qb.*Rs) );
b = find(y<1 & y>=0);
P1(b) = A .* ( besseli(m_mode,Qs.*r(b).*Rs) ./ besseli(m_mode,Qs.*Rs) );
b = find(y<0 & y>-1);
P1(b) = D .* ( besseli(m_mode,Qs.*r(b).*Rs) ./ besseli(m_mode,Qs.*Rs) );
b = find(y<=-1);
P1(b) = D .* ( besselk(m_mode,Qb.*r(b).*Rs) ./ besselk(m_mode,Qb.*Rs) );

P1 = P1.*exp(1i.*K.*x.*Rs);
Pmax = max(max(real(P1)))
P1 = real(P1./A);

nice_fig([0 0],[10 10],[-5 -5],[5 5],1,1,...
    '$z\:/\:R_{\rm s}$','$x\:/\:R_{\rm s}$','',...
    16,20,18,1,0,0);
set(gcf,'ColorMap',cmap);
surf(x,y,P1,'linestyle','none');
view(2);
caxis([-1 1]);
bar=colorbar('peer',gca,'FontSize',16,'FontName','Arial',...
	'CLim',[1 64],'YTick',[-1:0.2:1]);
set(get(bar,'Ylabel'),'String','P_1 / A',...
	'Interpreter','tex','Fontsize',24,'Rotation',270,'position',...
[10.5,0,9.16],'FontName','Times');

x = -1:0.0001:11;
y = 1 + real(0.1.*exp(1i*K.*x.*Rs));
z = 30 * ones(length(y),1);
plot3(x,y,z,'marker','none','linestyle','-','linewidth',2,'color','w');
y = -1 - real(D./A.*0.1.*exp(1i*K.*x.*Rs));
plot3(x,y,z,'marker','none','linestyle','-','linewidth',2,'color','w');

% Create title
tit = strcat('$M_{\rm b}=',num2str(M,'%3.1f'),'\:\: \delta=',num2str(Del,'%3i'),...
    '\:\: \lambda=',num2str(Lambda/Rs,'%3.1f'),'R_{\rm s}\:\: m=',num2str(m_mode,'%1i'),...
    '\:\: n=',num2str(n_mode,'%1i'),'\:\: y=0$');
title(tit,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 1.01 0.1]);%,'BackgroundColor',[1 1 1]);
set(gcf,'renderer','Zbuffer')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xy plane
[x,y] = meshgrid(linspace(-5,5,Ngrid),linspace(-5,5,Ngrid));
r = sqrt(x.^2+y.^2);
phi = r;
b = find(x==0 & y==0);
phi(b) = 0;
b = find(x==0 & y>0);
phi(b) = pi/2;
b = find(x==0 & y<0);
phi(b) = 3*pi/2;
b = find(x>0 & y==0);
phi(b) = 0;
b = find(x<0 & y==0);
phi(b) = pi;
b = find(x>0 & y>0);
phi(b) = atan(abs(y(b)./x(b)));
b = find(x<0 & y>0);
phi(b) = pi - atan(abs(y(b)./x(b)));
b = find(x<0 & y<0);
phi(b) = pi + atan(abs(y(b)./x(b)));
b = find(x>0 & y<0);
phi(b) = 2*pi - atan(abs(y(b)./x(b)));

P1 = zeros(Ngrid,Ngrid);

b = find(r>=1);
P1(b) = A.*(besseli(m_mode,Qs.*Rs)./besselk(m_mode,Qb.*Rs)).*besselk(m_mode,Qb.*r(b).*Rs);
b = find(r<1);
P1(b) = A.*besseli(m_mode,Qs.*r(b).*Rs);

P1 = P1.*exp(1i.*m_mode.*phi);
Pmax = max(max(real(P1)));
P1 = real(P1./Pmax);

nice_fig([-5 -5],[5 5],[-5 -5],[5 5],1,1,...
    '$x\:/\:R_{\rm s}$','$y\:/\:R_{\rm s}$','',...
    16,20,18,1,0,0);
set(gcf,'ColorMap',cmap);
surf(x,y,P1,'linestyle','none');
view(2);
caxis([-1 1]);
bar=colorbar('peer',gca,'FontSize',16,'FontName','Arial',...
	'CLim',[1 64],'YTick',[-1:0.2:1]);
set(get(bar,'Ylabel'),'String','P_1 / A',...
	'Interpreter','tex','Fontsize',24,'Rotation',270,'position',...
[10.5,0,9.16],'FontName','Times');

t = 0:0.0001:(2*pi);
R = 1 + 0.1.*cos(m_mode.*t);
x = R.*cos(t);
y = R.*sin(t);
z = 30 * ones(length(y),1);
plot3(x,y,z,'marker','none','linestyle','-','linewidth',2,'color','w');

% Create title
tit = strcat('$M_{\rm b}=',num2str(M,'%3.1f'),'\:\: \delta=',num2str(Del,'%3i'),...
    '\:\: \lambda=',num2str(Lambda/Rs,'%3.1f'),'R_{\rm s}\:\: m=',num2str(m_mode,'%1i'),...
    '\:\: n=',num2str(n_mode,'%1i'),'\:\: z=0$');
title(tit,'Interpreter','latex','FontSize',16,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 1.01 0.1]);%,'BackgroundColor',[1 1 1]);
set(gcf,'renderer','Zbuffer')