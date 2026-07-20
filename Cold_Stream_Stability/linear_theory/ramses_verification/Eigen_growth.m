function Eigen_growth(Nrad)

delt = [0.000157, 0.000334, 0.00149, 0.00204, 0.000236];
Omega_I = [318.5297967047999, 149.6194816042138, 33.595959445774078, 24.420383427003003, 211.7125140159174];
leg_tit1 = '$M_{\rm b}=1.5,\:\:\delta=1$';
leg_tit2 = '$M_{\rm b}=1.5,\:\:\delta=10\:\:{\rm S}$';
leg_tit3 = '$M_{\rm b}=1.5,\:\:\delta=10\:\:{\rm B}$';
leg_tit4 = '$M_{\rm b}=1.5,\:\:\delta=100$';
leg_tit5 = '$M_{\rm b}=5.0,\:\:\delta=1$';
leg_tit6 = '${\rm exp}(t/t_{\rm KH})$';
np = zeros(1,length(delt));
nvz = zeros(1,length(delt));
grid_direc = 'D:/shared_backup/Kelvin_Helmholtz/ramses/eigenmode_runs/grids/';

time = zeros(100,length(delt));
P1 = zeros(100,length(delt));
P1_2 = zeros(100,length(delt));
VZ = zeros(100,length(delt));
VZ_2 = zeros(100,length(delt));

time2 = linspace(0,5,300);
KHI = exp(time2);

for i=1:length(delt)
    if(i==1)
        direc = strcat(grid_direc,'M_1.5_Del_1/L15_sig4/');
    elseif(i==2)
        direc = strcat(grid_direc,'M_1.5_Del_10_surface/L15_sig4/');
    elseif(i==3)
        direc = strcat(grid_direc,'M_1.5_Del_10_body/L15_sig16/');
    elseif(i==4)
        direc = strcat(grid_direc,'M_1.5_Del_100/L15_sig4/');
    elseif(i==5)
        direc = strcat(grid_direc,'M_5.0_Del_1/L15_sig4/');
    end

    % Read file names and identify data
    filelist1 = dir(direc);
    filelist = filelist1(3:length(filelist1));
    n = length(filelist);
    
    for j=1:n
        if(filelist(j).name(5)=='p')
            np(i) = np(i)+1;
            pres_file(np(i)) = filelist(j);
        end
        if(filelist(j).name(5:6)=='vy')
            nvz(i) = nvz(i)+1;
            vz_file(nvz(i)) = filelist(j);
        end
    end
    np(i) = min([np(i),100]);
    nvz(i) = min([nvz(i),100]);
    
    if(np(i)>1)
        time(1:np(i),i) = Omega_I(i).*[0:delt(i):(np(i)-1)*delt(i)];
    elseif(nvz(i)>1)
        time(1:nvz(i),i) = Omega_I(i).*[0:delt(i):(nvz(i)-1)*delt(i)];
    end
    
    % Read in pressure
    for j=1:np(i)
        [i,j,1]
        test = pres_file(j).name;
        test2 = strcat(direc,test);
        fid = fopen(test2, 'rb');
        ntemp = fread(fid,1, 'int');
        nx = fread(fid,1, 'int');
        ny = fread(fid,1, 'int');
        ntemp = fread(fid,1, 'int');
        rtemp = fread(fid,1, 'float32');
        pres = fread(fid, nx*ny, 'float32');
        rtemp = fread(fid,1, 'float64');
        xmin = fread(fid, 1, 'float64');
        xmax = fread(fid, 1, 'float64');
        rtemp = fread(fid,1, 'float64');
        ymin = fread(fid, 1, 'float64');
        ymax = fread(fid, 1, 'float64');
        fclose(fid);
        
        % Check smaller box %
        pres = reshape(pres,[ny,nx]);
        pres = permute(pres,[2,1]);
        if(j==1)
            [x,y] = meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));
            rad = 1/160;
            b = find(abs(y-0.5)<=Nrad*rad);
        end
        %P1(j,i) = mean(abs(pres(b)-1));
        %P1_2(j,i) = sqrt(mean((pres(b)-1).^2));
        P1(j,i) = sum(abs(pres(b)-1));
        P1_2(j,i) = sqrt(sum((pres(b)-1).^2));
    end
    P1(:,i) = P1(:,i)./P1(1,i);
    P1_2(:,i) = P1_2(:,i)./P1_2(1,i);
    
    % Read in velocity
    for j=1:nvz(i)
        [i,j,2]
        test = vz_file(j).name;
        test2 = strcat(direc,test);
        fid = fopen(test2, 'rb');
        ntemp = fread(fid,1, 'int');
        nx = fread(fid,1, 'int');
        ny = fread(fid,1, 'int');
        ntemp = fread(fid,1, 'int');
        rtemp = fread(fid,1, 'float32');
        vel = fread(fid, nx*ny, 'float32');
        rtemp = fread(fid,1, 'float64');
        xmin = fread(fid, 1, 'float64');
        xmax = fread(fid, 1, 'float64');
        rtemp = fread(fid,1, 'float64');
        ymin = fread(fid, 1, 'float64');
        ymax = fread(fid, 1, 'float64');
        fclose(fid);
        
        % Check smaller box %
        vel = reshape(vel,[ny,nx]);
        vel = permute(vel,[2,1]);
        if(j==1)
            [x,y] = meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));
            rad = 1/160;
            b = find(abs(y-0.5)<=Nrad*rad);
        end
        %VZ(j,i) = mean(abs(vel(b)));
        %VZ_2(j,i) = sqrt(mean((vel(b)).^2));
        VZ(j,i) = sum(abs(vel(b)));
        VZ_2(j,i) = sqrt(sum((vel(b)).^2));
    end
    VZ(:,i) = VZ(:,i)./VZ(1,i);
    VZ_2(:,i) = VZ_2(:,i)./VZ_2(1,i);
end

xl = [0 0];
xu = [5 5];
xd = 0.5;
yl = [0.5 0.5];
yu = [100 100];
ytick = [0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100];
ystr = {num2str(0.5,'%3.1f'),' ',' ',' ',' ',num2str(1,'%1.0f'),' ',' ',' ',' ',' ',' ',' ',' ',num2str(10,'%2.0f'),' ',' ',' ',' ',' ',' ',' ',' ',num2str(100,'%3.0f')};
xlab = '$t/t_{\rm KH}$';
ylab = '$\left <P_1\right>/\left<P_{\rm 1,\,t=0}\right>$';
%tit = strcat('$Eigenmodes$');
tit = '';

nice_fig2(xl,xu,yl,yu,xd,ytick,ystr,xlab,ylab,tit,16,24,18,1,0,0)
l(6) = plot(time2, KHI, 'linestyle', '--', 'linewidth', 2, 'color', 'k', 'marker', 'none', 'DisplayName', leg_tit6);
l(1) = plot(time(:,1), P1(:,1), 'linestyle', '-', 'linewidth', 2, 'color', 'g', 'marker', 'none', 'DisplayName', leg_tit1);
l(2) = plot(time(:,2), P1(:,2), 'linestyle', '-', 'linewidth', 2, 'color', 'c', 'marker', 'none', 'DisplayName', leg_tit2);
l(3) = plot(time(:,3), P1(:,3), 'linestyle', '-', 'linewidth', 2, 'color', 'b', 'marker', 'none', 'DisplayName', leg_tit3);
l(4) = plot(time(:,4), P1(:,4), 'linestyle', '-', 'linewidth', 2, 'color', 'm', 'marker', 'none', 'DisplayName', leg_tit4);
l(5) = plot(time(:,5), P1(:,5), 'linestyle', '-', 'linewidth', 2, 'color', 'r', 'marker', 'none', 'DisplayName', leg_tit5);
legend1 = legend(gca,l(1:6));
set(legend1,'edgecolor','w','Units','Normalized',...
    'fontname','Times New Roman','fontsize',14,'Interpreter','Latex',...
    'Position',[0.525, 0.127, 0.342, 0.317]);
set(gcf,'renderer','painters')

nice_fig2(xl,xu,yl,yu,xd,ytick,ystr,xlab,ylab,tit,16,24,18,1,0,0)
l(6) = plot(time2, KHI, 'linestyle', '--', 'linewidth', 2, 'color', 'k', 'marker', 'none', 'DisplayName', leg_tit6);
l(1) = plot(time(:,1), P1_2(:,1), 'linestyle', '-', 'linewidth', 2, 'color', 'g', 'marker', 'none', 'DisplayName', leg_tit1);
l(2) = plot(time(:,2), P1_2(:,2), 'linestyle', '-', 'linewidth', 2, 'color', 'c', 'marker', 'none', 'DisplayName', leg_tit2);
l(3) = plot(time(:,3), P1_2(:,3), 'linestyle', '-', 'linewidth', 2, 'color', 'b', 'marker', 'none', 'DisplayName', leg_tit3);
l(4) = plot(time(:,4), P1_2(:,4), 'linestyle', '-', 'linewidth', 2, 'color', 'm', 'marker', 'none', 'DisplayName', leg_tit4);
l(5) = plot(time(:,5), P1_2(:,5), 'linestyle', '-', 'linewidth', 2, 'color', 'r', 'marker', 'none', 'DisplayName', leg_tit5);
legend1 = legend(gca,l(1:6));
set(legend1,'edgecolor','w','Units','Normalized',...
    'fontname','Times New Roman','fontsize',14,'Interpreter','Latex',...
    'Position',[0.525, 0.127, 0.342, 0.317]);
set(gcf,'renderer','painters')

ylab = '$\left <u_{\rm x}\right>/\left<u_{\rm x,\,t=0}\right>$';
nice_fig2(xl,xu,yl,yu,xd,ytick,ystr,xlab,ylab,tit,16,24,18,1,0,0)
l(6) = plot(time2, KHI, 'linestyle', '--', 'linewidth', 2, 'color', 'k', 'marker', 'none', 'DisplayName', leg_tit6);
l(1) = plot(time(:,1), VZ(:,1), 'linestyle', '-', 'linewidth', 2, 'color', 'g', 'marker', 'none', 'DisplayName', leg_tit1);
l(2) = plot(time(:,2), VZ(:,2), 'linestyle', '-', 'linewidth', 2, 'color', 'c', 'marker', 'none', 'DisplayName', leg_tit2);
l(3) = plot(time(:,3), VZ(:,3), 'linestyle', '-', 'linewidth', 2, 'color', 'b', 'marker', 'none', 'DisplayName', leg_tit3);
l(4) = plot(time(:,4), VZ(:,4), 'linestyle', '-', 'linewidth', 2, 'color', 'm', 'marker', 'none', 'DisplayName', leg_tit4);
l(5) = plot(time(:,5), VZ(:,5), 'linestyle', '-', 'linewidth', 2, 'color', 'r', 'marker', 'none', 'DisplayName', leg_tit5);
legend1 = legend(gca,l(1:6));
set(legend1,'edgecolor','w','Units','Normalized',...
    'fontname','Times New Roman','fontsize',14,'Interpreter','Latex',...
    'Position',[0.525, 0.127, 0.342, 0.317]);
set(gcf,'renderer','painters')

nice_fig2(xl,xu,yl,yu,xd,ytick,ystr,xlab,ylab,tit,16,24,18,1,0,0)
l(6) = plot(time2, KHI, 'linestyle', '--', 'linewidth', 2, 'color', 'k', 'marker', 'none', 'DisplayName', leg_tit6);
l(1) = plot(time(:,1), VZ_2(:,1), 'linestyle', '-', 'linewidth', 2, 'color', 'g', 'marker', 'none', 'DisplayName', leg_tit1);
l(2) = plot(time(:,2), VZ_2(:,2), 'linestyle', '-', 'linewidth', 2, 'color', 'c', 'marker', 'none', 'DisplayName', leg_tit2);
l(3) = plot(time(:,3), VZ_2(:,3), 'linestyle', '-', 'linewidth', 2, 'color', 'b', 'marker', 'none', 'DisplayName', leg_tit3);
l(4) = plot(time(:,4), VZ_2(:,4), 'linestyle', '-', 'linewidth', 2, 'color', 'm', 'marker', 'none', 'DisplayName', leg_tit4);
l(5) = plot(time(:,5), VZ_2(:,5), 'linestyle', '-', 'linewidth', 2, 'color', 'r', 'marker', 'none', 'DisplayName', leg_tit5);
legend1 = legend(gca,l(1:6));
set(legend1,'edgecolor','w','Units','Normalized',...
    'fontname','Times New Roman','fontsize',14,'Interpreter','Latex',...
    'Position',[0.525, 0.127, 0.342, 0.317]);
set(gcf,'renderer','painters')

% filename = strcat(fig_direc,num2str(delt*(i-1),'%6.4f'),'.jpg');
% saveas(gca,filename);
