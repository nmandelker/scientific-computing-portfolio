function [dens,xvec,yvec,vx,vy,p]=read_ramses_2d(direc,delt,n0,pl,rho1,rho2,lambda,urel,pres,diam,gama)
% Take every n0-ith snapshot
% If pl=0, then don't plot
% If pl=1, then plot snapshot sequence, density only
% If pl=2, then plot snapshot sequence, with density and velocity fields

cs = sqrt(gama * pres * rho1);
tcross = diam / cs;
tkh = (rho1 + rho2) * lambda / sqrt(rho1 * rho2 * urel.^2);
% Read file names and identify data
filelist1 = dir(direc);
filelist = filelist1(3:length(filelist1));
n = length(filelist)
nd = 0;
np = 0;
nvx = 0;
nvy = 0;
for i=1:n
    if(filelist(i).name(6)=='e')
        nd = nd+1;
        dens_file(nd) = filelist(i);
    elseif(filelist(i).name(6)=='r')
        np = np+1;
        pres_file(np) = filelist(i);
    elseif(filelist(i).name(6)=='x')
        nvx = nvx+1;
        vx_file(nvx) = filelist(i);
    elseif(filelist(i).name(6)=='y')
        nvy = nvy+1;
        vy_file(nvy) = filelist(i);
    end
end
[nd nvx nvy np]
n1 = floor(nd/n0)

% Read in density
test = dens_file(1).name;
test2 = strcat(direc,test);
fid = fopen(test2, 'rb');
ntemp = fread(fid,1, 'int');
nx = fread(fid,1, 'int');
ny = fread(fid,1, 'int');
dens(:,:)=zeros(nx*ny,n1);
ntemp = fread(fid,1, 'int');
dens(:,1) = fread(fid, nx*ny, 'float32');
rtemp = fread(fid,1, 'float32');
rtemp = fread(fid,1, 'float64');
xmin = fread(fid, 1, 'float64');
xmax = fread(fid, 1, 'float64');
rtemp = fread(fid,1, 'float64');
ymin = fread(fid, 1, 'float64');
ymax = fread(fid, 1, 'float64');
fclose(fid);
    
for i=2:n1
    i
    %l1 = length(test);
    %time(i,1:5) = test(l1-8:l1-4);
    test = dens_file(1+n0*(i-1)).name;
    test2 = strcat(direc,test);
    fid = fopen(test2, 'rb');
    ntemp = fread(fid,1, 'int');
    nx = fread(fid,1, 'int');
    ny = fread(fid,1, 'int');
    ntemp = fread(fid,1, 'int');
    dens(:,i) = fread(fid, nx*ny, 'float32');
    fclose(fid);
end
dens = reshape(dens,[ny,nx,n1]);
dens = permute(dens,[2 1 3]);
xvec = linspace(xmin,xmax,nx);
yvec = linspace(ymin,ymax,ny);

% Read in velocities
if(nargout>3)
    vx(:,:)=zeros(nx*ny,n1);
    vy(:,:)=zeros(nx*ny,n1);
    
    if(nvx>0)
        for i=1:n1
            i
            test = vx_file(1+n0*(i-1)).name;
            test2 = strcat(direc,test);
            fid = fopen(test2, 'rb');
            ntemp = fread(fid,1, 'int');
            nx = fread(fid,1, 'int');
            ny = fread(fid,1, 'int');
            ntemp = fread(fid,1, 'int');
            vx(:,i) = fread(fid, nx*ny, 'float32');
            fclose(fid);
        end
    end

    if(nvy>0)
        for i=1:n1
            i
            test = vy_file(1+n0*(i-1)).name;    
            test2 = strcat(direc,test);
            fid = fopen(test2, 'rb');
            ntemp = fread(fid,1, 'int');
            nx = fread(fid,1, 'int');
            ny = fread(fid,1, 'int');
            ntemp = fread(fid,1, 'int');
            vy(:,i) = fread(fid, nx*ny, 'float32');
            fclose(fid);
        end
    end
    vx = reshape(vx,[ny,nx,n1]);
    vy = reshape(vy,[ny,nx,n1]);
    vx = permute(vx,[2 1 3]);
    vy = permute(vy,[2 1 3]);
end

% Read in pressure
if(nargout>4)
    p(:,:)=zeros(nx*ny,n1);
    if(np>0)
        for i=1:n1
            i
            test = pres_file(1+n0*(i-1)).name;
            test2 = strcat(direc,test);
            fid = fopen(test2, 'rb');
            ntemp = fread(fid,1, 'int');
            nx = fread(fid,1, 'int');
            ny = fread(fid,1, 'int');
            ntemp = fread(fid,1, 'int');
            p(:,i) = fread(fid, nx*ny, 'float32');
            fclose(fid);
        end
    end
    p = reshape(p,[ny,nx,n1]);
    p = permute(p,[2 1 3]);
end

if(pl==1)
    plot_ramses_2d(direc,delt,1,dens,xvec,yvec,tkh,tcross);
elseif(pl==2)
    plot_ramses_2d(direc,delt,1,dens,xvec,yvec,tkh,tcross,vx,vy);
end

