function Tsc_over_Tvir = compressible_sheet_growth_times(a_over_Rvir, lambda_over_a, Mh, Del);

mach = 'tot';

Del = [0.01:0.1:100];
Mb = [0:0.03:3];
Ms = sqrt(max(Del)).*Mb;
Mtot = Mb;

nD = length(Del);
nM = length(Mb);
R = zeros(nD,nM);
if(strcmp(mach,'b'))
    titx = '$M_{\rm b}$';
    for j=1:nM
        for i=1:nD
            p = [Mb(j)^2, -2*Mb(j)^2, Mb(j)^2-1-1/Del(i), 2, -1];
            R(i,j) = max(imag(roots(p)));
        end
    end
elseif(strcmp(mach,'s'))
    titx = '$M_{\rm s}$';
    for j=1:nM
        for i=1:nD
            p = [(Ms(j)/sqrt(Del(i)))^2, -2*(Ms(j)/sqrt(Del(i)))^2, (Ms(j)/sqrt(Del(i)))^2-1-1/Del(i), 2, -1];
            R(i,j) = max(imag(roots(p)));
        end
    end
elseif(strcmp(mach,'tot'))
    titx = '$M_{\rm tot}$';
    for j=1:nM
        for i=1:nD
            Mt = (1+sqrt(Del(i)))/(sqrt(Del(i)));
            p = [(Mt*Mtot(j))^2, -2*(Mt*Mtot(j))^2, (Mt*Mtot(j))^2-1-1/Del(i), 2, -1];
            R(i,j) = max(imag(roots(p)));
        end
    end
end

cmap = [1 1 1;0.125 0.125 1;0.1071 0.1071 1;0.08929 0.08929 1;...
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

xl = [0 0];
if(strcmp(mach,'s'))
    xu = [30 30];
    dx = 3;
else
    xu = [3 3];
    dx = 0.5;
end
yl = [1 0];
yu = [100 100];
dy = 10;
nice_fig(xl,xu,yl,yu,dx,dy,...
    titx,'$\delta$','$Compressible\:\:Vortex\:\:Sheet$',...
    16,20,18,1,1,0);
set(gca,'YTick',[1 10 20 30 40 50 60 70 80 90 100]);
set(gcf,'ColorMap',cmap);
if(strcmp(mach,'b'))
    surf(Mb,Del,R,'linestyle','none');
elseif(strcmp(mach,'s'))
    surf(Ms,Del,R,'linestyle','none');
elseif(strcmp(mach,'tot'))
    surf(Mtot,Del,R,'linestyle','none');
end
view(2);
bar=colorbar('peer',gca,'FontSize',16,'FontName','Arial',...
	'CLim',[1 64],'YTick',[0:0.1:0.5]);
caxis([0 0.5]);
set(get(bar,'Ylabel'),'String','${\rm Im}\left(\Phi\right)$',...
	'Interpreter','Latex','Fontsize',18,'Rotation',270,'position',...
[10,0.25,9.16],'FontName','Times');
