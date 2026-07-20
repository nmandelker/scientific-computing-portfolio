function cyl_mode_structure(Del, Mh, m_mode, write_tit)

% In directory
Num_m_modes = 6;
Length = 250001;

Mc = sqrt(Del)*Mh;
Mtot = 1/(1/Mh+1/Mc)
Nm = length(m_mode);

if(Del>=1)
    Delform = '%3i';
else
    Delform = '%3.1f';
end    
dirname1 = strcat('./M_',num2str(Mh,'%3.1f'),'_del_',num2str(Del,Delform),'_modes/');
dirname = strcat(dirname1,'data/');
tit = strcat('$\delta=',num2str(Del,Delform),...
    '\:\:\: M_{\rm b}=',num2str(Mh,'%3.1f'),...
    '\:\:\: M_{\rm s}=',num2str(Mc,'%3.1f'),...
    '\:\:\: M_{\rm tot}\simeq',num2str(Mtot,'%4.2f'),'$');

filelist1 = dir(dirname);
filelist = filelist1(3:length(filelist1));
n = length(filelist)
nmode = n/(Num_m_modes*3);
X = zeros(Nm,nmode,Length);
ReP = zeros(Nm,nmode,Length);
ImP = zeros(Nm,nmode,Length);
Na = zeros(Nm,nmode,1);
for i=1:nmode
    k = i-1;
    for ii=1:Nm
        kk = m_mode(ii);
        for j=1:n
            n1=length(filelist(j).name);
            m2 = str2num(filelist(j).name(n1-4));
            if(filelist(j).name(n1-8)=='n')
                n2 = str2num(filelist(j).name(n1-7));
            else
                n2 = str2num(filelist(j).name(n1-8:n1-7));
            end
            if(n2==k & m2==kk)
                if(filelist(j).name(1)=='x')
                    a = load(strcat(dirname,filelist(j).name));
                    Na(ii,i) = min(Length,length(a));
                    X(ii,i,1:Na(ii,i)) = abs(a(1:Na(ii,i)));
                elseif(filelist(j).name(1)=='p' & filelist(j).name(5)=='i')
                    a = load(strcat(dirname,filelist(j).name));
                    Na(ii,i) = min(Length,length(a));
                    ImP(ii,i,1:Na(ii,i)) = abs(a(1:Na(ii,i)));
                elseif(filelist(j).name(1)=='p' & filelist(j).name(5)=='r')
                    a = load(strcat(dirname,filelist(j).name));
                    Na(ii,i) = min(Length,length(a));
                    ReP(ii,i,1:Na(ii,i)) = abs(a(1:Na(ii,i)));
                end
            end
        end
    end
end

col1(:,1) = [0, 0, 0];
col1(:,2) = [0.4, 0.4, 0.4];
col1(:,3) = [0.7, 0.7, 0.7];
col1(:,4) = [0.9, 0.9, 0.9];
col1(:,5) = [0.8 0.5 0.4];

col(:,1) = [0.2 0.2 0.8];
col(:,2) = [1 0 0];
col(:,3) = [0 1 0];
col(:,4) = [0 0 1];
col(:,5) = [0 0 0];
col(:,6) = [1 0 1];
col(:,7) = [0 1 1];
col(:,8) = [0.8 0.8 0.2];
col(:,9) = [0.4 0.2 0.8];
col(:,10) = [0.8 0.4 0.2];
col(:,11) = [0.2 0.8 0.5];
ncol = 11;

xlow = 0.01;
xup = 100;

x_incomp = linspace(xlow,xup,100000);
incomp_growth = zeros(Nm,100000);
incomp_real = zeros(Nm,100000);
for i=1:Nm
    k = m_mode(i);
    Im = besseli(k,x_incomp);
    Km = besselk(k,x_incomp);
    if(k==0)
        ImPrime =  besseli(1,x_incomp);
        KmPrime = -besselk(1,x_incomp);
    else
        ImPrime =  0.5 .* ( besseli(k+1,x_incomp) + besseli(k-1,x_incomp) );
        KmPrime = -0.5 .* ( besselk(k+1,x_incomp) + besselk(k-1,x_incomp) );
    end
    bes_root = 1 ./ sqrt( Del .* Im .* KmPrime ./ (ImPrime .* Km) );
    Phi_incomp_1 = 1 ./ ( 1 - bes_root ); 
    Phi_incomp_2 = 1 ./ ( 1 + bes_root );
    if( imag(Phi_incomp_1) > imag(Phi_incomp_2) )
        Phi_incomp = Phi_incomp_1;
    else
        Phi_incomp = Phi_incomp_2;
    end
    incomp_growth(i,:) = 2 .* x_incomp .* Mc .* imag(Phi_incomp);
    incomp_real(i,:)   = 2 .* x_incomp .* Mc .* real(Phi_incomp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ylow = 0.01;
yup = 100;
tick = [0.01 0.1 1 10 100];
str = {num2str(0.01,'%4.2f'),num2str(0.1,'%3.1f'),...
    num2str(1,'%1.0f'),num2str(10,'%2.0f'),num2str(100,'%3.0f')};

figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
set(figure1,'WindowStyle','docked');
axes1 = axes('Parent',figure1,'YTick',tick,...
    'YScale','log',...
    'YtickLabel',str,...
    'YMinorTick','on',...'
    'XTick',tick,...
    'XScale','log',...
    'XtickLabel',str,...
    'XMinorTick','on',...'
    'TickLength',[0.02 0.04],...
    'PlotBoxAspectRatio',[1 1 1],...
    'LineWidth',1.5,...
    'FontSize',16,...
    'FontName','Arial');

xlim(axes1,[xlow xup]);
ylim(axes1,[ylow yup]);
box(axes1,'on');
hold(axes1,'all');

% Create xlabel
xlabel('$K\:=\:kR_{\rm s}\:=\:2\pi R_{\rm s} \:/\:\lambda$',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FontName','Times New Roman','units','normalized','position',[0.5 -0.06 0]);

% Create ylabel
ylabel('$\omega_{_{\rm I}}\cdot 2R_{\rm s}/c_{\rm s}\:=\:t_{\rm sc}\:/\:t_{\rm KH}$','Interpreter','latex','FontSize',24,...
    'FontName','Times New Roman',...
    'units','normalized','position',[-0.085 0.5 0]);

% Create title
if(write_tit)
    title(tit,...
        'Interpreter','latex','FontSize',18,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 1.01 0]);
end

set(figure1,'renderer','painters')

for k=2:nmode
    xplot=reshape(X(1,k,1:Na(i,k))  ,[Na(1,k),1]);
    yplot=reshape(ImP(1,k,1:Na(i,k)),[Na(1,k),1]);
    plot(xplot,2.*Mc.*xplot.*yplot,'linestyle','-','linewidth',2,...
        'color','r','marker','none');
end
l(1) = plot(x_incomp,incomp_growth(1,:),...
    'linestyle','--','linewidth',2,...
    'color','k','marker','none','DisplayName','$incomp$');
xplot=reshape(X(1,1,1:Na(1,1))  ,[Na(1,1),1]);
yplot=reshape(ImP(1,1,1:Na(1,1)),[Na(1,1),1]);
l(2) = plot(xplot,2.*Mc.*xplot.*yplot,'linestyle','-','linewidth',2,...
    'color','k','marker','none','DisplayName','$n=0$');
l(3) = plot(1e6.*xplot,1e6.*2.*Mc.*xplot.*yplot,'linestyle','-','linewidth',2,...
    'color','r','marker','none','DisplayName','$n=1-20$');
if(Nm>1)
    for i=2:Nm
%        k = m_mode(i);
%        j = mod(i,ncol)+1;
%        c = col(:,j);
        for k=2:nmode
            xplot=reshape(X(i,k,1:Na(i,k))  ,[Na(i,k),1]);
            yplot=reshape(ImP(i,k,1:Na(i,k)),[Na(i,k),1]);
            plot(xplot,2.*Mc.*xplot.*yplot,'linestyle','-','linewidth',2,...
                'color','r','marker','none');
        end
        plot(x_incomp,incomp_growth(1,:),...
            'linestyle','--','linewidth',2,...
            'color','k','marker','none');
        xplot=reshape(X(i,1,1:Na(i,1))  ,[Na(i,1),1]);
        yplot=reshape(ImP(i,1,1:Na(i,1)),[Na(i,1),1]);
        plot(xplot,2.*Mc.*xplot.*yplot,'linestyle','-','linewidth',2,...
            'color','k','marker','none');
    end
end

y1 = logspace(-4,3,8);
x1 = 2*pi.*y1./y1;
x2 = pi.*y1./y1;
x3 = 4*pi.*y1./y1;
plot(x1,y1,'marker','none','linestyle',':','linewidth',2,'color','k');
plot(x2,y1,'marker','none','linestyle',':','linewidth',2,'color','k');
plot(x3,y1,'marker','none','linestyle',':','linewidth',2,'color','k');

x4 = 1:1:1000;
y4 = log(4.*Mh.*((sqrt(Del)/(1+sqrt(Del))).^2).*x4);
y5 = y4 - log(y4);
y6 = y4 - log(y5);
l(4) = plot(x4,y6,'marker','none','linestyle','-.','linewidth',4,'color','k',...
    'DisplayName','$\omega_{\rm I,\:3}$');

legend1 = legend(axes1,l(1:4));
set(legend1,'edgecolor','w','Units','Normalized',...
    'fontname','Times New Roman','fontsize',14,'Interpreter','Latex',...
    'Position',[0.18, 0.65+0.05, 0.30, 0.2]);
legend boxoff

m_string = strcat('m',{' '},'=',{' '},num2str(m_mode(1),'%1i'));
if(Nm>1)
    for i=2:Nm
        m_string = strcat(m_string,',',{' '},num2str(m_mode(i),'%1i'));
    end
end

annotation(gcf,'textbox',[0.47 0.84 0.11 0.06],...
    'String',m_string,...
    'FontSize',18,...
    'FontName','Times',...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ylow = 0.1;
yup = 1000;
xlow = 0.01;
xup = 100;
xtick = [0.01 0.1 1 10 100];
ytick = [0.1 1 10 100 1000];
xstr = {num2str(0.01,'%4.2f'),num2str(0.1,'%3.1f'),...
    num2str(1,'%1.0f'),num2str(10,'%2.0f'),num2str(100,'%3.0f')};
ystr = {num2str(0.1,'%3.1f'),num2str(1,'%1.0f'),...
    num2str(10,'%2.0f'),num2str(100,'%3.0f'),num2str(1000,'%4.0f')};

figure2 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
set(figure2,'WindowStyle','docked');
axes2 = axes('Parent',figure2,'YTick',ytick,...
    'YScale','log',...
    'YtickLabel',ystr,...
    'YMinorTick','on',...'
    'XTick',xtick,...
    'XScale','log',...
    'XtickLabel',xstr,...
    'XMinorTick','on',...'
    'TickLength',[0.02 0.04],...
    'PlotBoxAspectRatio',[1 1 1],...
    'LineWidth',1.5,...
    'FontSize',16,...
    'FontName','Arial');

xlim(axes2,[xlow xup]);
ylim(axes2,[ylow yup]);
box(axes2,'on');
hold(axes2,'all');

% Create xlabel
xlabel('$K\:=\:kR_{\rm s}\:=\:2\pi R_{\rm s} \:/\:\lambda$',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FontName','Times New Roman','units','normalized','position',[0.5 -0.06 0]);

% Create ylabel
ylabel('$\omega_{_{\rm R}}\cdot 2R_{\rm s}/c_{\rm s}\:=\:2\pi\,t_{\rm sc}\:/\:t_{\rm period}$','Interpreter','latex','FontSize',24,...
    'FontName','Times New Roman',...
    'units','normalized','position',[-0.085 0.5 0]);

% Create title
if(write_tit)
    title(tit,...
        'Interpreter','latex','FontSize',18,...
        'FontName','Times New Roman','units','normalized',...
        'position',[0.5 1.01 0]);
end

set(figure2,'renderer','painters')

for k=2:nmode
    xplot=reshape(X(1,k,1:Na(i,k))  ,[Na(1,k),1]);
    yplot=reshape(ReP(1,k,1:Na(i,k)),[Na(1,k),1]);
    plot(xplot,2.*Mc.*xplot.*yplot,'linestyle','-','linewidth',2,...
        'color','r','marker','none');
end
p(1) = plot(x_incomp,incomp_real(1,:),...
    'linestyle','--','linewidth',2,...
    'color','k','marker','none','DisplayName','$incomp$');
xplot=reshape(X(1,1,1:Na(1,1))  ,[Na(1,1),1]);
yplot=reshape(ReP(1,1,1:Na(1,1)),[Na(1,1),1]);
p(2) = plot(xplot,2.*Mc.*xplot.*yplot,'linestyle','-','linewidth',2,...
    'color','k','marker','none','DisplayName','$n=0$');
p(3) = plot(1e6.*xplot,1e6.*2.*Mc.*xplot.*yplot,'linestyle','-','linewidth',2,...
    'color','r','marker','none','DisplayName','$n=1-20$');
if(Nm>1)
    for i=2:Nm
%        k = m_mode(i);
%        j = mod(i,ncol)+1;
%        c = col(:,j);
        for k=2:nmode
            xplot=reshape(X(i,k,1:Na(i,k))  ,[Na(i,k),1]);
            yplot=reshape(ReP(i,k,1:Na(i,k)),[Na(i,k),1]);
            plot(xplot,2.*Mc.*xplot.*yplot,'linestyle','-','linewidth',2,...
                'color','r','marker','none');
        end
        plot(x_incomp,incomp_real(1,:),...
            'linestyle','--','linewidth',2,...
            'color','k','marker','none');
        xplot=reshape(X(i,1,1:Na(i,1))  ,[Na(i,1),1]);
        yplot=reshape(ReP(i,1,1:Na(i,1)),[Na(i,1),1]);
        plot(xplot,2.*Mc.*xplot.*yplot,'linestyle','-','linewidth',2,...
            'color','k','marker','none');
    end
end

y1 = logspace(-4,3,8);
x1 = 2*pi.*y1./y1;
x2 = pi.*y1./y1;
x3 = 4*pi.*y1./y1;
plot(x1,y1,'marker','none','linestyle',':','linewidth',2,'color','k');
plot(x2,y1,'marker','none','linestyle',':','linewidth',2,'color','k');
plot(x3,y1,'marker','none','linestyle',':','linewidth',2,'color','k');

legend2 = legend(axes2,p(1:3));
set(legend2,'edgecolor','w','Units','Normalized',...
    'fontname','Times New Roman','fontsize',14,'Interpreter','Latex',...
    'Position',[0.18, 0.65+0.1, 0.30, 0.15]);
legend boxoff

annotation(gcf,'textbox',[0.47 0.84 0.11 0.06],...
    'String',m_string,...
    'FontSize',18,...
    'FontName','Times',...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor','w');