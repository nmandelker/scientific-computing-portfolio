function KHI_phase_diagrams_analytic(Diam_over_Rvir, lambda_over_Diam)

Mh = [0:0.01:3];
LogDel = 0:0.005:2;
Del = 10.^LogDel;

Tvir_over_Tkh = zeros(length(Del),length(Mh));
for i=1:length(Mh)
    for j=1:length(Del)
        Mc = sqrt(Del(j)).*Mh(i);
        Mtot = ( sqrt(Del(j))./(1+sqrt(Del(j))) ).*Mh(i);
        Tsc_over_Tvir = Diam_over_Rvir .* Mc;
        if( Mtot < 1 )
            p = [Mh(i)^2, -2*Mh(i)^2, Mh(i)^2-1-1/Del(j), 2, -1]; % Quartic polynomial for compressible sheet
            Im_Phi = max(imag(roots(p)));
            Tvir_over_Tkh(j,i) = (2*pi*Im_Phi) / (lambda_over_Diam) / (Diam_over_Rvir);
        else
            p = [Mh(i)^2, -2*Mh(i)^2, Mh(i)^2-1-1/Del(j), 2, -1]; % Quartic polynomial for compressible sheet
            Im_Phi = max(imag(roots(p)));
            Tvir_over_Tkh1 = (2*pi*Im_Phi) / (lambda_over_Diam) / (Diam_over_Rvir);
            
            y = 4.*Mtot .* ( sqrt(Del(j))/(1+sqrt(Del(j))) ) .* pi/lambda_over_Diam;
            Tvir_over_Tkh2 = (log(y) - log(log(y)-log(log(y)))) ./ Tsc_over_Tvir;
            Tvir_over_Tkh2 = Tvir_over_Tkh2 .* (1 - Tsc_over_Tvir);
            
            Tvir_over_Tkh(j,i) = max([1e-20, Tvir_over_Tkh1, Tvir_over_Tkh2]);
        end
    end
end

% Mh2 = 0.75:0.01:2.25;
% Del2 = 10:0.1:110;
% Ts_Over_Tk = zeros(length(Del2),length(Mh2));
% for i=1:length(Mh2)
%     for j=1:length(Del2)
%         Mc = sqrt(Del2(j)).*Mh2(i);
%         Ts_Over_Tv = Diam_over_Rvir .* Mc;
%         
%         p = [Mh2(i)^2, -2*Mh2(i)^2, Mh2(i)^2-1-1/Del2(j), 2, -1]; % Quartic polynomial for compressible sheet
%         Im_Phi = max(imag(roots(p)));
%         Tv_Over_Tk = (2*pi*Im_Phi) / (lambda_over_Diam) / (Diam_over_Rvir);
%         Ts_Over_Tk(j,i) = Ts_Over_Tv * Tv_Over_Tk;
%     end
% end
% coherence = 110.*Mh2./Mh2;
% for i=1:length(Mh2)
%     if(min(Ts_Over_Tk(:,i))<=1 & max(Ts_Over_Tk(:,i))>=1)
%         b = find(abs(Ts_Over_Tk(:,i)-1) == min(abs(Ts_Over_Tk(:,i)-1)), 1, 'first');
%         coherence(i) = Del2(b);
%     elseif(min(Ts_Over_Tk(:,i))<=1 & max(Ts_Over_Tk(:,i))<=1)
%         coherence(i) = 10;
%     elseif(min(Ts_Over_Tk(:,i))>=1 & max(Ts_Over_Tk(:,i))>=1)
%         coherence(i) = 110;
%     end
% end

cmap = [0 0 0;0 0 0.25;0 0 0.5;0 0 0.75;0 0 1;0 0.07143 0.9929;0 0.1429 0.9857;
    0 0.2143 0.9786;0 0.2857 0.9714;0 0.3571 0.9643;0 0.4286 0.9571;0 0.5 0.95;
    0 0.5714 0.9429;0 0.6429 0.9357;0 0.7143 0.9286;0 0.7857 0.9214;0 0.8571 0.9143;
    0 0.9286 0.9071;0 1 0.9;0 1 0.825;0 1 0.75;0 1 0.675;0 1 0.6;0 1 0.525;0 1 0.45;
    0 1 0.375;0 1 0.3;0 1 0.225;0 1 0.15;0 1 0.075;0 1 0;0.15 1 0;0.3 1 0;0.45 1 0;
    0.6 1 0;0.75 1 0;0.9 1 0;0.92 1 0;0.94 1 0;0.96 1 0;0.98 1 0;1 1 0;1 0.9352 0;
    1 0.8704 0;1 0.8056 0;1 0.7407 0;1 0.6759 0;1 0.6111 0;1 0.584 0;1 0.5568 0;
    1 0.5296 0;1 0.5025 0;1 0.4753 0;1 0.3565 0;1 0.2377 0;1 0.1188 0;1 0 0;1 0.2 0.2;
    1 0.4 0.4;1 0.6 0.6;1 0.8 0.8;1 0.8667 0.8667;1 0.9333 0.9333;1 1 1];

% Create figure
figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68],'Colormap',cmap);
set(figure1,'WindowStyle','docked');

xt = 0:0.5:3;
yt = 0:0.2:2;
cl = [0 2.5 0.5];
% Create axes
axes1 = axes('Parent',figure1,...
    'YTick',yt,...
    'XTick',xt,...
    'Position',[0.12 0.15 0.8 0.8],...
    'PlotBoxAspectRatio',[1 1 1],...
    'FontSize',16,...
    'FontName','Arial',...
    'CLim',[cl(1) cl(2)]);

xlim([xt(1) xt(length(xt))]);
ylim([yt(1) yt(length(yt))]);
grid('off');
hold('all');

%surf(Mh,Del,Tkh_over_Tvir,'Parent',axes1,'linestyle','none');
surf(Mh,LogDel,log10(Tvir_over_Tkh),'Parent',axes1,'linestyle','none');
titx = '$M_{\rm b}$';
tity = '${\rm log}\:\:(\delta)$';
titc = '${\rm log}\:\:(N_{\rm e\:\:foldings}\:\:\:\:)$';
tit = strcat('$2R_{\rm s}/R_{\rm vir} = ',num2str(Diam_over_Rvir,'%4.2f'),'\:\:\:\: \lambda/2R_{\rm s} = ',num2str(lambda_over_Diam,'%4.2f'),'$');

xlabel(titx,'Interpreter','latex','FontSize',22,...
    'units','normalized','position',[0.5 -0.08 0],...
    'FontName','Times New Roman');

ylabel(tity,'Interpreter','latex','FontSize',22,...
    'units','normalized','position',[-0.08 0.5 0],...
    'FontName','Times New Roman');

title(tit,'Interpreter','latex','FontSize',18,...
    'units','normalized','position',[0.41 1.01 0],...
    'FontName','Times New Roman');%,'BackGroundColor','w');

% Create colorbar
bar=colorbar('peer',axes1,'FontSize',16,'FontName','Arial',...
    'CLim',[1 64],'YTick',[cl(1):cl(3):cl(2)]);
set(get(bar,'Ylabel'),'String',titc,...
    'Interpreter','Latex','Fontsize',18,'Rotation',270,'position',[10.0,(cl(1)+cl(2))/2,9.16],'FontName','Times');

set(gcf,'renderer','zbuffer')