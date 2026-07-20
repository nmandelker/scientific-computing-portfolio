function convergence(Del, type, Nrad)

col(1) = 'g';
col(2) = 'c';
col(3) = 'b';
col(4) = 'r';
grid_direc = 'D:/shared_backup/Kelvin_Helmholtz/ramses/eigenmode_runs/grids/';

if(Del==1 | Del==101)
    L = [13, 14, 15, 16];
    if(Del==1)
        grid_direc = strcat(grid_direc,'M_1.5_Del_1/');
        tit = '$M_{\rm b}=1.5\:\:\:\delta=1';
        delt = 0.000157;
        Omega_I = 318.5297967047999;
    elseif(Del==101)
        grid_direc = strcat(grid_direc,'M_1.5_Del_10_surface/');
        tit = '$M_{\rm b}=1.5\:\:\:\delta=10\:\:\:surface';
        delt = 0.000334;
        Omega_I = 149.6194816042138;
    end
    if(type==1)
        sig = [4, 4, 4, 4];
        tit = strcat(tit,'\:\:\:\sigma/\Delta=4$');
    elseif(type==2)
        sig = [4, 8, 16, 32];
        tit = strcat(tit,'\:\:\:\lambda/\sigma=25$');
    elseif(type==3)
        L = [13, 14, 15];
        sig = [2, 4, 8];
        tit = strcat(tit,'\:\:\:\lambda/\sigma=50$');
    elseif(type==4)
        L = [13, 14, 15];
        sig = [1, 2, 4];
        tit = strcat(tit,'\:\:\:\lambda/\sigma=100$');
    elseif(type==5)
        L = [14, 15];
        sig = [1, 2];
        tit = strcat(tit,'\:\:\:\lambda/\sigma=200$');
    end
elseif(Del==102 | Del==100)
    L = [13, 14, 15];
    if(Del==102)
        grid_direc = strcat(grid_direc,'M_1.5_Del_10_body/');
        tit = '$M_{\rm b}=1.5\:\:\:\delta=10\:\:\:body';
        delt = 0.00149;
        Omega_I = 33.595959445774078;
        Omega_I2 = 149.6194816042138;
    elseif(Del==100)
        grid_direc = strcat(grid_direc,'M_1.5_Del_100/');
        tit = '$M_{\rm b}=1.5\:\:\:\delta=100';
        delt = 0.00204;
        Omega_I = 24.420383427003003;
    end
    if(type==1)
        sig = [4, 4, 4];
        tit = strcat(tit,'\:\:\:\sigma/\Delta=4$');
    elseif(type==2)
        sig = [4, 8, 16];
        tit = strcat(tit,'\:\:\:\lambda/\sigma=25$');
    elseif(type==3)
        L = [13, 14, 15];
        sig = [2, 4, 8];
        tit = strcat(tit,'\:\:\:\lambda/\sigma=50$');
    elseif(type==4)
        L = [13, 14, 15];
        sig = [1, 2, 4];
        tit = strcat(tit,'\:\:\:\lambda/\sigma=100$');
    elseif(type==5)
        L = [14, 15];
        sig = [1, 2];
        tit = strcat(tit,'\:\:\:\lambda/\sigma=200$');
    end
end
Ndirec = length(L);

np = zeros(1,Ndirec);
nvz = zeros(1,Ndirec);
nh = zeros(1,Ndirec);

time = zeros(100,Ndirec);
P1 = zeros(100,Ndirec);
P1_2 = zeros(100,Ndirec);
VZ = zeros(100,Ndirec);
VZ_2 = zeros(100,Ndirec);
ptv = zeros(100,Ndirec);

time2 = linspace(0,5,300);
KHI = exp(time2);
if(Del==102)
    KHI2 = exp((Omega_I2/Omega_I).*time2);
end

for i=1:Ndirec
    grid_direc1 = strcat(grid_direc,'L',num2str(L(i),'%2i'),'_sig',num2str(sig(i),'%2i'),'/')

    % Read file names and identify data
    filelist1 = dir(grid_direc1);
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
        if(filelist(j).name(5:6)=='co')
            nh(i) = nh(i)+1;
            color_file(nh(i)) = filelist(j);
        end
    end
    np(i) = min([np(i),100]);
    nvz(i) = min([nvz(i),100]);
    nh(i) = min([nh(i),100]);
    
    if(np(i)>1)
        time(1:np(i),i)  = Omega_I.*[0:delt:(np(i)-1)*delt];
    elseif(nvz(i)>1)
        time(1:nvz(i),i) = Omega_I.*[0:delt:(nvz(i)-1)*delt];
    elseif(nh(i)>1)
        time(1:nh(i),i)  = Omega_I.*[0:delt:(nh(i)-1)*delt];
    end
    
    % Read in pressure
    for j=1:np(i)
        [i,j,1]
        test = pres_file(j).name;
        test2 = strcat(grid_direc1,test);
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
        P1(j,i) = mean(abs(pres(b)-1));
        P1_2(j,i) = sqrt(mean((pres(b)-1).^2));
    end
    P1(:,i) = P1(:,i)./P1(1,i);
    P1_2(:,i) = P1_2(:,i)./P1_2(1,i);
    clear pres x y
    
    % Read in velocity
    for j=1:nvz(i)
        [i,j,2]
        test = vz_file(j).name;
        test2 = strcat(grid_direc1,test);
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
        VZ(j,i) = mean(abs(vel(b)));
        VZ_2(j,i) = sqrt(mean((vel(b)).^2));
    end
    VZ(:,i) = VZ(:,i)./VZ(1,i);
    VZ_2(:,i) = VZ_2(:,i)./VZ_2(1,i);
    clear vel x y
    
    % Read in color
    for j=1:nh(i)
        [i,j,3]
        test = color_file(j).name;
        test2 = strcat(grid_direc1,test);
        fid = fopen(test2, 'rb');
        ntemp = fread(fid,1, 'int');
        nx = fread(fid,1, 'int');
        ny = fread(fid,1, 'int');
        ntemp = fread(fid,1, 'int');
        rtemp = fread(fid,1, 'float32');
        color = fread(fid, nx*ny, 'float32');
        rtemp = fread(fid,1, 'float64');
        xmin = fread(fid, 1, 'float64');
        xmax = fread(fid, 1, 'float64');
        rtemp = fread(fid,1, 'float64');
        ymin = fread(fid, 1, 'float64');
        ymax = fread(fid, 1, 'float64');
        fclose(fid);
        
        color = reshape(color,[ny,nx]);
        color = permute(color,[2 1]);
        x = linspace(xmin,xmax,nx);
        y = linspace(ymin,ymax,ny);        
        c = contourc(x, y, color, [0.5 0.5]);        
        %if the contour isn't continuous, extract only (x,y) data of all segments
        cc = [];
        fin = c(2,1)+1; %number of points in the next segment+1
        while fin<size(c,2)
            cc = [cc,c(1:2,2:fin)];
            c = c(:,fin+1:end);
            fin = c(2,1)+1;
        end
        cc = [cc,c(1:2,2:fin)];        
        b = find( cc(2,:) >= ((ymin+ymax)/2) );
        yup = max(cc(2,b));
        ydown = min(cc(2,b));
        ptv(j,i) = 0.5*(yup-ydown);
    end        
    clear color x y
    b = find(ptv(:,i)>0);
    ptv(:,i) = ptv(:,i)./ptv(b(1),i);
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

nice_fig2(xl,xu,yl,yu,xd,ytick,ystr,xlab,ylab,tit,16,24,18,1,0,0)
if(Del==102)
    l(Ndirec+2) = plot(time2+0.5, KHI2, 'linestyle', ':', 'linewidth', 2, 'color', 'k', ...
        'marker', 'none', 'DisplayName', '${\rm exp}(t/t_{\rm KH,\,surface})$');
end
l(Ndirec+1) = plot(time2, KHI, 'linestyle', '--', 'linewidth', 2, 'color', 'k', ...
    'marker', 'none', 'DisplayName', '${\rm exp}(t/t_{\rm KH})$');
for i=1:Ndirec
    l(i) = plot(time(:,i), P1(:,i), 'linestyle', '-', 'linewidth', 2, ...
        'color', col(i), 'marker', 'none', ...
        'DisplayName', strcat('$\lambda/\Delta=',num2str( floor(102.5*2^(L(i)-13))),'$' ));
end
if(Del==102)
    Ndirec=Ndirec+1;
end
legend1 = legend(gca,l(1:Ndirec+1));
set(legend1,'edgecolor','w','Units','Normalized',...
    'fontname','Times New Roman','fontsize',14,'Interpreter','Latex',...
    'Position',[0.60, 0.15, 0.24, 0.12+0.06*(Ndirec-1)]);
set(gcf,'renderer','painters')
if(Del==102)
    Ndirec=Ndirec-1;
end

nice_fig2(xl,xu,yl,yu,xd,ytick,ystr,xlab,ylab,tit,16,24,18,1,0,0)
if(Del==102)
    l(Ndirec+2) = plot(time2+0.5, KHI2, 'linestyle', ':', 'linewidth', 2, 'color', 'k', ...
        'marker', 'none', 'DisplayName', '${\rm exp}(t/t_{\rm KH,\,surface})$');
end
l(Ndirec+1) = plot(time2, KHI, 'linestyle', '--', 'linewidth', 2, 'color', 'k', ...
    'marker', 'none', 'DisplayName', '${\rm exp}(t/t_{\rm KH})$');
for i=1:Ndirec
    l(i) = plot(time(:,i), P1_2(:,i), 'linestyle', '-', 'linewidth', 2, ...
        'color', col(i), 'marker', 'none', ...
        'DisplayName', strcat('$\lambda/\Delta=',num2str( floor(102.5*2^(L(i)-13))),'$' ));
end
if(Del==102)
    Ndirec=Ndirec+1;
end
legend1 = legend(gca,l(1:Ndirec+1));
set(legend1,'edgecolor','w','Units','Normalized',...
    'fontname','Times New Roman','fontsize',14,'Interpreter','Latex',...
    'Position',[0.60, 0.15, 0.24, 0.12+0.06*(Ndirec-1)]);
set(gcf,'renderer','painters')
if(Del==102)
    Ndirec=Ndirec-1;
end

ylab = '$\left <u_{\rm x}\right>/\left<u_{\rm x,\,t=0}\right>$';
nice_fig2(xl,xu,yl,yu,xd,ytick,ystr,xlab,ylab,tit,16,24,18,1,0,0)
if(Del==102)
    l(Ndirec+2) = plot(time2+0.5, KHI2, 'linestyle', ':', 'linewidth', 2, 'color', 'k', ...
        'marker', 'none', 'DisplayName', '${\rm exp}(t/t_{\rm KH,\,surface})$');
end
l(Ndirec+1) = plot(time2, KHI, 'linestyle', '--', 'linewidth', 2, 'color', 'k', ...
    'marker', 'none', 'DisplayName', '${\rm exp}(t/t_{\rm KH})$');
for i=1:Ndirec
    l(i) = plot(time(:,i), VZ(:,i), 'linestyle', '-', 'linewidth', 2, ...
        'color', col(i), 'marker', 'none', ...
        'DisplayName', strcat('$\lambda/\Delta=',num2str( floor(102.5*2^(L(i)-13))),'$' ));
end
if(Del==102)
    Ndirec=Ndirec+1;
end
legend1 = legend(gca,l(1:Ndirec+1));
set(legend1,'edgecolor','w','Units','Normalized',...
    'fontname','Times New Roman','fontsize',14,'Interpreter','Latex',...
    'Position',[0.60, 0.15, 0.24, 0.12+0.06*(Ndirec-1)]);
set(gcf,'renderer','painters')
if(Del==102)
    Ndirec=Ndirec-1;
end

nice_fig2(xl,xu,yl,yu,xd,ytick,ystr,xlab,ylab,tit,16,24,18,1,0,0)
if(Del==102)
    l(Ndirec+2) = plot(time2+0.5, KHI2, 'linestyle', ':', 'linewidth', 2, 'color', 'k', ...
        'marker', 'none', 'DisplayName', '${\rm exp}(t/t_{\rm KH,\,surface})$');
end
l(Ndirec+1) = plot(time2, KHI, 'linestyle', '--', 'linewidth', 2, 'color', 'k', ...
    'marker', 'none', 'DisplayName', '${\rm exp}(t/t_{\rm KH})$');
for i=1:Ndirec
    l(i) = plot(time(:,i), VZ_2(:,i), 'linestyle', '-', 'linewidth', 2, ...
        'color', col(i), 'marker', 'none', ...
        'DisplayName', strcat('$\lambda/\Delta=',num2str( floor(102.5*2^(L(i)-13))),'$' ));
end
if(Del==102)
    Ndirec=Ndirec+1;
end
legend1 = legend(gca,l(1:Ndirec+1));
set(legend1,'edgecolor','w','Units','Normalized',...
    'fontname','Times New Roman','fontsize',14,'Interpreter','Latex',...
    'Position',[0.60, 0.15, 0.24, 0.12+0.06*(Ndirec-1)]);
set(gcf,'renderer','painters')
if(Del==102)
    Ndirec=Ndirec-1;
end

ylab = '$h/h_{t=0}$';
nice_fig2(xl,xu,yl,yu,xd,ytick,ystr,xlab,ylab,tit,16,24,18,1,0,0)
if(Del==102)
    l(Ndirec+2) = plot(time2+0.5, KHI2, 'linestyle', ':', 'linewidth', 2, 'color', 'k', ...
        'marker', 'none', 'DisplayName', '${\rm exp}(t/t_{\rm KH,\,surface})$');
end
l(Ndirec+1) = plot(time2, KHI, 'linestyle', '--', 'linewidth', 2, 'color', 'k', ...
    'marker', 'none', 'DisplayName', '${\rm exp}(t/t_{\rm KH})$');
for i=1:Ndirec
    l(i) = plot(time(:,i), ptv(:,i), 'linestyle', '-', 'linewidth', 2, ...
        'color', col(i), 'marker', 'none', ...
        'DisplayName', strcat('$\lambda/\Delta=',num2str( floor(102.5*2^(L(i)-13))),'$' ));
end
if(Del==102)
    Ndirec=Ndirec+1;
end
legend1 = legend(gca,l(1:Ndirec+1));
set(legend1,'edgecolor','w','Units','Normalized',...
    'fontname','Times New Roman','fontsize',14,'Interpreter','Latex',...
    'Position',[0.60, 0.15, 0.24, 0.12+0.06*(Ndirec-1)]);
set(gcf,'renderer','painters')
if(Del==102)
    Ndirec=Ndirec-1;
end

% filename = strcat(fig_direc,num2str(delt*(i-1),'%6.4f'),'.jpg');
% saveas(gca,filename);
