function turbulence_density_maps(gal, version, Sg, Vr, Vp, Vz, sr, sp, sz, scm, scv, ssm, ssv, st, socm, socv,...
    cSg, cVr, cVp, cVz, csr, csp, csz, csc, css, cst, csoc, clumps)
% gal should be galaxy number '**' with no ./ or /
%
% Sg, Vr, Vp, Vz, sr, sp, sz, sc, ss, st and soc are 0 or 1, if you want 
% to plot gas surface density, radial velocity, azimuthal velocity, 
% vertical velocity, radial dispersion, azimuthal dispersion, vetical 
% dispersion, compressive turbulence amplitude, solenoidal dispersion 
% amplitude, total dispersion amplitude, solenoidal to compressive ratio
%
% All tracers are plotted face. 
%
% cSg, cVr, cVp, cVz, csr, csp, csz, csc, css, cst and csoc are vectors 
% with 3 entries: Lower and upper limits for color bar and tick-step.

ft = 'float64';
if(version==0)
    dirname0 = './outputs/3D_top_hat_FW_0900pc/NOT_subtract_rot/vol_weight_local/NOT_smoothed/';
elseif(version==1)
    dirname0 = './outputs/3D_top_hat_FW_0900pc/NOT_subtract_rot/mass_weight_local/NOT_smoothed/';
elseif(version==2)
    dirname0 = './outputs/3D_top_hat_FW_0900pc/subtract_rot/mass_weight_local/NOT_smoothed/';
elseif(version==3)
    dirname0 = './outputs/3D_top_hat_FW_0900pc/subtract_rot/mass_weight_local/smoothed_BEFORE/';
elseif(version==4)
    dirname0 = './outputs/3D_Gaussian_FW_0900pc/subtract_rot/mass_weight_local/smoothed_BEFORE/';
elseif(version==5)
    dirname0 = './outputs/2D_Gaussian_FW_0900pc/subtract_rot/mass_weight_local/smoothed_BEFORE/';
elseif(version==6)
    dirname0 = './outputs/2D_Gaussian_FW_0600pc/subtract_rot/mass_weight_local/smoothed_BEFORE/';
elseif(version==7)
    dirname0 = './outputs/2D_Gaussian_FW_0300pc/subtract_rot/mass_weight_local/smoothed_BEFORE/';
elseif(version==8)
    dirname0 = './outputs/grid_res_200pc/2D_Gaussian_FW_0600pc/subtract_rot/mass_weight_local/smoothed_BEFORE/';
elseif(version==15)
    dirname0 = './outputs/thin_slices/';
elseif(version==18)
	dirname0 = './outputs/grid_res_070pc/3D_top_hat_FW_0500pc/NOT_subtract_rot/NOT_subtract_loc/mass_weight_local/decomp_vel/NOT_smoothed/';
elseif(version==19)
    dirname0 = './outputs/grid_res_070pc/3D_top_hat_FW_0500pc/NOT_subtract_rot/subtract_loc/mass_weight_local/decomp_vel/NOT_smoothed/';
elseif(version==20)
    dirname0 = './outputs/grid_res_070pc/3D_top_hat_FW_0500pc/subtract_rot/NOT_subtract_loc/mass_weight_local/decomp_vel/NOT_smoothed/';
elseif(version==21)
    dirname0 = './outputs/grid_res_070pc/3D_top_hat_FW_0500pc/subtract_rot/subtract_loc/mass_weight_local/decomp_vel/NOT_smoothed/';
elseif(version==22)
    dirname0 = './outputs/grid_res_250pc/3D_top_hat_FW_0500pc/NOT_subtract_rot/NOT_subtract_loc/mass_weight_local/decomp_vel/NOT_smoothed/';
elseif(version==23)
    dirname0 = './outputs/grid_res_250pc/3D_top_hat_FW_0500pc/NOT_subtract_rot/subtract_loc/mass_weight_local/decomp_vel/NOT_smoothed/';
elseif(version==24)
    dirname0 = './outputs/grid_res_250pc/3D_top_hat_FW_0500pc/subtract_rot/NOT_subtract_loc/mass_weight_local/decomp_vel/NOT_smoothed/';
elseif(version==25)
    dirname0 = './outputs/grid_res_250pc/3D_top_hat_FW_0500pc/subtract_rot/subtract_loc/mass_weight_local/decomp_vel/NOT_smoothed/';
end
dirname0 = strcat(dirname0,'/VELA_v2_',num2str(gal,'%02i'),'/binary_grid_outputs/');

dirname2 = strcat('./outputs/figs/version',num2str(version,'%02i'),'/surface_plots/');
mkdir(dirname2);

dirname2 = strcat(dirname2,'VELA_v2_',num2str(gal,'%02i'),'/');
mkdir(dirname2);

if(clumps == 0)
    if(Sg == 1)
        fSg = strcat(dirname2,'Sigma_gas/');
        mkdir(fSg);
    end
    if(Vr == 1)
        fVr = strcat(dirname2,'Vr/');
        mkdir(fVr);
    end
    if(Vp == 1)
        fVp = strcat(dirname2,'Vphi/');
        mkdir(fVp);
    end
    if(Vz == 1)
        fVz = strcat(dirname2,'Vz/');
        mkdir(fVz);
    end
    if(sr == 1)
        fsr = strcat(dirname2,'sigma_r/');
        mkdir(fsr);
    end
    if(sp == 1)
        fsp = strcat(dirname2,'sigma_phi/');
        mkdir(fsp);
    end
    if(sz == 1)
        fsz = strcat(dirname2,'sigma_z/');
        mkdir(fsz);
    end
    if(scm == 1)
        fscm = strcat(dirname2,'sigma_c_energy/');
        mkdir(fscm);
    end
    if(scv == 1)
        fscv = strcat(dirname2,'sigma_c_momentum/');
        mkdir(fscv);
    end
    if(ssm == 1)
        fssm = strcat(dirname2,'sigma_s_energy/');
        mkdir(fssm);  
    end      
    if(ssv == 1)
        fssv = strcat(dirname2,'sigma_s_momentum/');
        mkdir(fssv);
    end
    if(st == 1)
        fst = strcat(dirname2,'sigma_tot/');
        mkdir(fst);
    end
    if(socm == 1)
        fsocm = strcat(dirname2,'C_over_S_energy/');
        mkdir(fsocm);
    end
    if(socv == 1)
        fsocv = strcat(dirname2,'C_over_S_momentum/');
        mkdir(fsocv);
    end
elseif(clumps == 1)
    if(Sg == 1)
        fSg = strcat(dirname2,'Sigma_gas_with_clumps/');
        mkdir(fSg);
    end
    if(Vr == 1)
        fVr = strcat(dirname2,'Vr_with_clumps/');
        mkdir(fVr);
    end
    if(Vp == 1)
        fVp = strcat(dirname2,'Vphi_with_clumps/');
        mkdir(fVp);
    end
    if(Vz == 1)
        fVz = strcat(dirname2,'Vz_with_clumps/');
        mkdir(fVz);
    end
    if(sr == 1)
        fsr = strcat(dirname2,'sigma_r_with_clumps/');
        mkdir(fsr);
    end
    if(sp == 1)
        fsp = strcat(dirname2,'sigma_phi_with_clumps/');
        mkdir(fsp);
    end
    if(sz == 1)
        fsz = strcat(dirname2,'sigma_z_with_clumps/');
        mkdir(fsz);
    end
    if(scm == 1)
        fscm = strcat(dirname2,'sigma_c_mass_with_clumps/');
        mkdir(fscm);
    end
    if(scv == 1)
        fscv = strcat(dirname2,'sigma_c_volume_with_clumps/');
        mkdir(fscv);
    end
    if(ssm == 1)
        fssm = strcat(dirname2,'sigma_s_mass_with_clumps/');
        mkdir(fssm);
    end
    if(ssv == 1)
        fssv = strcat(dirname2,'sigma_s_volume_with_clumps/');
        mkdir(fssv);
    end
    if(st == 1)
        fst = strcat(dirname2,'sigma_tot_with_clumps/');
        mkdir(fst);
    end
    if(socm == 1)
        fsocm = strcat(dirname2,'compressive_to_solenoidal_ratio_mass_with_clumps/');
        mkdir(fsocm);
    end
    if(socv == 1)
        fsocv = strcat(dirname2,'compressive_to_solenoidal_ratio_volume_with_clumps/');
        mkdir(fsocv);
    end
end

filelist1 = dir(dirname0);
ndir = length(filelist1)-2;
filelist = filelist1(3:ndir+2);
cmap = [0 0 0;0 0 0.25;0 0 0.5;0 0 0.75;0 0 1;0 0.06478 0.9935;...
    0 0.1296 0.987;0 0.1943 0.9806;0 0.2591 0.9741;0 0.3239 0.9676;...
    0 0.3887 0.9611;0 0.4534 0.9547;0 0.5182 0.9482;0 0.583 0.9417;...
    0 0.6478 0.9352;0 0.7126 0.9287;0 0.7773 0.9223;0 0.8421 0.9158;...
    0 0.8816 0.9118;0 0.9211 0.9079;0 0.9605 0.9039;0 1 0.9;0 1 0.795;...
    0 1 0.69;0 1 0.585;0 1 0.48;0 1 0.375;0 1 0.27;0 1 0.2025;0 1 0.135;...
    0 1 0.0675;0 1 0;0.1 1 0;0.2 1 0;0.3 1 0;0.4 1 0;0.5 1 0;0.6 1 0;...
    0.7 1 0;0.8 1 0;0.9 1 0;0.9333 1 0;0.9667 1 0;1 1 0;1 0.9222 0;...
    1 0.8444 0;1 0.7667 0;1 0.6889 0;1 0.6111 0;1 0.5771 0;1 0.5432 0;...
    1 0.5093 0;1 0.4753 0;1 0.3565 0;1 0.2377 0;1 0.1188 0;1 0 0;...
    1 0.2 0.2;1 0.4 0.4;1 0.6 0.6;1 0.8 0.8;1 0.8667 0.8667;...
    1 0.9333 0.9333;1 1 1];
for i=1:ndir
	filename = strcat(dirname0,filelist(i).name);
	proceed = false;
% 	filelist(i).name
% 	filelist(i).name(1)
% 	filelist(i).name(7)
	if(Sg==1 & filelist(i).name(1)=='S' & filelist(i).name(7)=='g')
        proceed = true;
        colg = 'w';
        cl = [cSg(1) cSg(2)];
        cvec = [cl(1):cSg(3):cl(2)];
        clab = '${\rm log(}\:\:{\Sigma}_{\rm gas}{\rm \:)[\:M_{\odot}\:pc^{-2}\:]}$';
        filename2 = strcat(fSg,'a',filename(length(filename)-8:length(filename)-4),'.jpg');
        fid = fopen(filename, 'rb','ieee-be');
        ntemp = fread(fid,1, 'int')
        n = fread(fid, 1, 'int')
        n = fread(fid, 1, 'int')
        s1 = fread(fid, 1, 'float32')
        s2 = fread(fid, 1, 'float32')
        face = fread(fid, n * n, ft);
        fclose(fid);
        n1 = sqrt(length(face));
        [n, n1]
        n = n1;
        face = reshape(face, [n, n]);
        face = permute(face, [2, 1]);
        b = find(face == 0);
        face(b) = 1e-9;
        face = log10(face);
        [min(min(face)), max(max(face))]
    elseif(Vr==1 & filelist(i).name(1)=='V' & filelist(i).name(2)=='r')
        proceed = true;
        colg = 'k';
        cl = [cVr(1) cVr(2)];
        cvec = [cl(1):cVr(3):cl(2)];
        clab = '${\rm V_r\:[\:km\:sec^{-1}\:]}$';
        filename2 = strcat(fVr,'a',filename(length(filename)-8:length(filename)-4),'.jpg');
        fid = fopen(filename, 'rb','ieee-be');
        ntemp = fread(fid,1, 'int')
        n = fread(fid, 1, 'int')
        n = fread(fid, 1, 'int')
        s1 = fread(fid, 1, 'float32')
        s2 = fread(fid, 1, 'float32')
        face = fread(fid, n * n, ft);
        fclose(fid);
        n1 = sqrt(length(face));
        [n, n1]
        n = n1;
        face = reshape(face, [n, n]);
        face = permute(face, [2, 1]);
    elseif(Vp==1 & filelist(i).name(1)=='V' & filelist(i).name(2)=='p')
        proceed = true;
        colg = 'k';
        cl = [cVp(1) cVp(2)];
        cvec = [cl(1):cVp(3):cl(2)];
        clab = '${\rm V_{\phi}\:[\:km\:sec^{-1}\:]}$';
        filename2 = strcat(fVp,'a',filename(length(filename)-8:length(filename)-4),'.jpg');
        fid = fopen(filename, 'rb','ieee-be');
        ntemp = fread(fid,1, 'int')
        n = fread(fid, 1, 'int')
        n = fread(fid, 1, 'int')
        s1 = fread(fid, 1, 'float32')
        s2 = fread(fid, 1, 'float32')
        face = fread(fid, n * n, ft);
        fclose(fid);
        n1 = sqrt(length(face));
        [n, n1]
        n = n1;
        face = reshape(face, [n, n]);
        face = permute(face, [2, 1]);
    elseif(Vz==1 & filelist(i).name(1)=='V' & filelist(i).name(2)=='z')
        proceed = true;
        colg = 'k';
        cl = [cVz(1) cVz(2)];
        cvec = [cl(1):cVz(3):cl(2)];
        clab = '${\rm V_z\:[\:km\:sec^{-1}\:]}$';
        filename2 = strcat(fVz,'a',filename(length(filename)-8:length(filename)-4),'.jpg');
        fid = fopen(filename, 'rb','ieee-be');
        ntemp = fread(fid,1, 'int')
        n = fread(fid, 1, 'int')
        n = fread(fid, 1, 'int')
        s1 = fread(fid, 1, 'float32')
        s2 = fread(fid, 1, 'float32')
        face = fread(fid, n * n, ft);
        fclose(fid);
        n1 = sqrt(length(face));
        [n, n1]
        n = n1;
        face = reshape(face, [n, n]);
        face = permute(face, [2, 1]);
    elseif(sr==1 & filelist(i).name(1)=='s' & filelist(i).name(7)=='r')
        proceed = true;
        colg = 'k';
        cl = [csr(1) csr(2)];
        cvec = [cl(1):csr(3):cl(2)];
        clab = '${\rm \sigma_r\:[\:km\:sec^{-1}\:]}$';
        filename2 = strcat(fsr,'a',filename(length(filename)-8:length(filename)-4),'.jpg');
        fid = fopen(filename, 'rb','ieee-be');
        ntemp = fread(fid,1, 'int')
        n = fread(fid, 1, 'int')
        n = fread(fid, 1, 'int')
        s1 = fread(fid, 1, 'float32')
        s2 = fread(fid, 1, 'float32')
        face = fread(fid, n * n, ft);
        fclose(fid);
        n1 = sqrt(length(face));
        [n, n1]
        n = n1;
        face = reshape(face, [n, n]);
        face = permute(face, [2, 1]);
    elseif(sp==1 & filelist(i).name(1)=='s' & filelist(i).name(7)=='p')
        proceed = true;
        colg = 'k';
        cl = [csp(1) csp(2)];
        cvec = [cl(1):csp(3):cl(2)];
        clab = '${\rm \sigma_{\phi}\:[\:km\:sec^{-1}\:]}$';
        filename2 = strcat(fsp,'a',filename(length(filename)-8:length(filename)-4),'.jpg');
        fid = fopen(filename, 'rb','ieee-be');
        ntemp = fread(fid,1, 'int')
        n = fread(fid, 1, 'int')
        n = fread(fid, 1, 'int')
        s1 = fread(fid, 1, 'float32')
        s2 = fread(fid, 1, 'float32')
        face = fread(fid, n * n, ft);
        fclose(fid);
        n1 = sqrt(length(face));
        [n, n1]
        n = n1;
        face = reshape(face, [n, n]);
        face = permute(face, [2, 1]);
    elseif(sz==1 & filelist(i).name(1)=='s' & filelist(i).name(7)=='z')
        proceed = true;
        colg = 'k';
        cl = [csz(1) csz(2)];
        cvec = [cl(1):csz(3):cl(2)];
        clab = '${\rm \sigma_z\:[\:km\:sec^{-1}\:]}$';
        filename2 = strcat(fsz,'a',filename(length(filename)-8:length(filename)-4),'.jpg');
        fid = fopen(filename, 'rb','ieee-be');
        ntemp = fread(fid,1, 'int')
        n = fread(fid, 1, 'int')
        n = fread(fid, 1, 'int')
        s1 = fread(fid, 1, 'float32')
        s2 = fread(fid, 1, 'float32')
        face = fread(fid, n * n, ft);
        fclose(fid);
        n1 = sqrt(length(face));
        [n, n1]
        n = n1;
        face = reshape(face, [n, n]);
        face = permute(face, [2, 1]);
    elseif(scm==1 & filelist(i).name(1)=='s' & filelist(i).name(9)=='m' & filelist(i).name(7)=='c')
        proceed = true;
        colg = 'w';
        cl = [csc(1) csc(2)];
        cvec = [cl(1):csc(3):cl(2)];
        clab = '${\rm \sigma_c\:[\:km\:sec^{-1}\:]}$';
        filename2 = strcat(fscm,'a',filename(length(filename)-8:length(filename)-4),'.jpg');
        fid = fopen(filename, 'rb','ieee-be');
        ntemp = fread(fid,1, 'int')
        n = fread(fid, 1, 'int')
        n = fread(fid, 1, 'int')
        s1 = fread(fid, 1, 'float32')
        s2 = fread(fid, 1, 'float32')
        face = fread(fid, n * n, ft);
        fclose(fid);
        n1 = sqrt(length(face));
        [n, n1]
        n = n1;
        face = reshape(face, [n, n]);
        face = permute(face, [2, 1]);
        face = sqrt(face);
    elseif(scv==1 & filelist(i).name(1)=='s' & filelist(i).name(9)=='v' & filelist(i).name(7)=='c')
        proceed = true;
        colg = 'w';
        cl = [csc(1) csc(2)];
        cvec = [cl(1):csc(3):cl(2)];
        clab = '${\rm \sigma_c\:[\:km\:sec^{-1}\:]}$';
        filename2 = strcat(fscv,'a',filename(length(filename)-8:length(filename)-4),'.jpg');
        fid = fopen(filename, 'rb','ieee-be');
        ntemp = fread(fid,1, 'int')
        n = fread(fid, 1, 'int')
        n = fread(fid, 1, 'int')
        s1 = fread(fid, 1, 'float32')
        s2 = fread(fid, 1, 'float32')
        face = fread(fid, n * n, ft);
        fclose(fid);
        n1 = sqrt(length(face));
        [n, n1]
        n = n1;
        face = reshape(face, [n, n]);
        face = permute(face, [2, 1]);
        face = sqrt(face);
    elseif(ssm==1 & filelist(i).name(1)=='s' & filelist(i).name(9)=='m'  & filelist(i).name(7)=='s')
        proceed = true;
        colg = 'w';
        cl = [css(1) css(2)];
        cvec = [cl(1):css(3):cl(2)];
        clab = '${\rm \sigma_s\:[\:km\:sec^{-1}\:]}$';
        filename2 = strcat(fssm,'a',filename(length(filename)-8:length(filename)-4),'.jpg');
        fid = fopen(filename, 'rb','ieee-be');
        ntemp = fread(fid,1, 'int')
        n = fread(fid, 1, 'int')
        n = fread(fid, 1, 'int')
        s1 = fread(fid, 1, 'float32')
        s2 = fread(fid, 1, 'float32')
        face = fread(fid, n * n, ft);
        fclose(fid);
        n1 = sqrt(length(face));
        [n, n1]
        n = n1;
        face = reshape(face, [n, n]);
        face = permute(face, [2, 1]);
        face = sqrt(face);
    elseif(ssv==1 & filelist(i).name(1)=='s' & filelist(i).name(9)=='v'  & filelist(i).name(7)=='s')
        proceed = true;
        colg = 'w';
        cl = [css(1) css(2)];
        cvec = [cl(1):css(3):cl(2)];
        clab = '${\rm \sigma_s\:[\:km\:sec^{-1}\:]}$';
        filename2 = strcat(fssv,'a',filename(length(filename)-8:length(filename)-4),'.jpg');
        fid = fopen(filename, 'rb','ieee-be');
        ntemp = fread(fid,1, 'int')
        n = fread(fid, 1, 'int')
        n = fread(fid, 1, 'int')
        s1 = fread(fid, 1, 'float32')
        s2 = fread(fid, 1, 'float32')
        face = fread(fid, n * n, ft);
        fclose(fid);
        n1 = sqrt(length(face));
        [n, n1]
        n = n1;
        face = reshape(face, [n, n]);
        face = permute(face, [2, 1]);
        face = sqrt(face);
    elseif(st==1 & filelist(i).name(1)=='s' & filelist(i).name(7)=='t')
        proceed = true;
        colg = 'w';
        cl = [cst(1) cst(2)];
        cvec = [cl(1):cst(3):cl(2)];
        clab = '${\rm \sigma_{tot}\:[\:km\:sec^{-1}\:]}$';
        filename2 = strcat(fst,'a',filename(length(filename)-8:length(filename)-4),'.jpg');
        fid = fopen(filename, 'rb','ieee-be');
        ntemp = fread(fid,1, 'int')
        n = fread(fid, 1, 'int')
        n = fread(fid, 1, 'int')
        s1 = fread(fid, 1, 'float32')
        s2 = fread(fid, 1, 'float32')
        face = fread(fid, n * n, ft);
        fclose(fid);
        n1 = sqrt(length(face));
        [n, n1]
        n = n1;
        face = reshape(face, [n, n]);
        face = permute(face, [2, 1]);
    elseif((socm==1 | socv==1) & filelist(i).name(1)=='C')
        if(socm==1 & filelist(i).name(10)=='m')
            proceed = true;
            colg = 'k';
            cl = [csoc(1) csoc(2)];
            cvec = [cl(1):csoc(3):cl(2)];
            clab = '${\rm log(\:\:2\:\times\:\left(\sigma_c\:/\:\sigma_s \right)\:\:^2\:\:)}$';
            filename2 = strcat(fsocm,'a',filename(length(filename)-8:length(filename)-4),'.jpg');
            fid = fopen(filename, 'rb','ieee-be');
            ntemp = fread(fid,1, 'int')
            n = fread(fid, 1, 'int')
            n = fread(fid, 1, 'int')
            s1 = fread(fid, 1, 'float32')
            s2 = fread(fid, 1, 'float32')
            face = fread(fid, n * n, ft);
            fclose(fid);
            n1 = sqrt(length(face));
            [n, n1]
            n = n1;
            face = reshape(face, [n, n]);
            face = permute(face, [2, 1]);
            b = find(face == 0);
            face(b) = 1e-9;
            face = log10(2.*(face));
        elseif(socv==1 & filelist(i).name(10)=='v')
            proceed = true;
            colg = 'k';
            cl = [csoc(1) csoc(2)];
            cvec = [cl(1):csoc(3):cl(2)];
            clab = '${\rm log(\:\:2\:\times\:\left(\sigma_c\:/\:\sigma_s \right)\:\:^2\:\:)}$';
            filename2 = strcat(fsocv,'a',filename(length(filename)-8:length(filename)-4),'.jpg');
            fid = fopen(filename, 'rb','ieee-be');
            ntemp = fread(fid,1, 'int')
            n = fread(fid, 1, 'int')
            n = fread(fid, 1, 'int')
            s1 = fread(fid, 1, 'float32')
            s2 = fread(fid, 1, 'float32')
            face = fread(fid, n * n, ft);
            fclose(fid);
            n1 = sqrt(length(face));
            [n, n1]
            n = n1;
            face = reshape(face, [n, n]);
            face = permute(face, [2, 1]);
            b = find(face == 0);
            face(b) = 1e-9;
            face = log10(2.*(face));
        end
    end
	if(proceed == true)
        s1 = floor(s1);        
        x = linspace(-s1,s1,n);
        y = x;
        if(s1<8)
            d=1;
        elseif(s1==15 | s1==18 | s1==21)
            d=3;
        elseif(s1<16)
            d=2;
        elseif(s1>=16)
            d=4;
        end

        xlab = '$x\:[kpc]$';
        ylab = '$y\:[kpc]$';
        figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68],...
            'Colormap',cmap,'Visible','Off');
        set(figure1,'WindowStyle','docked');

	% Create axes
        axes1 = axes('Parent',figure1,'YTick',[-s1:d:s1],...
        'XTick',[-s1:d:s1],...
        'PlotBoxAspectRatio',[1 1 1],...
        'FontSize',16,...
        'FontName','Times',...
        'CLim',[cl(1) cl(2)]);
        box(axes1,'on');
        axis square
         
	% Create textbox
    annotation(figure1,'textbox',[0.1675 0.8605 0.12 0.064],...
        'String',strcat('a=',filename(length(filename)-8:length(filename)-4)),...
        'FontSize',16,...
        'FontName','Times',...
        'FitBoxToText','off',...
        'LineStyle','none',...
        'BackgroundColor',[1 1 1]);

        ylabel(ylab,'Interpreter','latex','FontSize',16,'FontName','Times','units','normalized',...
            'position',[-0.05,0.5,0.0]);
        xlabel(xlab,'Interpreter','latex','FontSize',16,'FontName','Times','units','normalized',...
            'position',[0.5,-0.065,0.0]);
        bar=colorbar('peer',axes1,'FontSize',16,'FontName','Times',...
            'CLim',[1 64],'YTick',cvec);
        set(get(bar,'Ylabel'),'String',clab,'Interpreter','latex','Fontsize',16,...
            'Rotation',270,'position',[10,(cl(1)+cl(2))/2,9.16],'FontName','Times');
        xlim([-s1 s1]);
        ylim([-s1 s1]);
        hold('all');
        [max(max(imag(face))),max(max(real(face))),max(max(imag(face)./real(face)))]
        surf(x,y,face,'Parent',axes1,'LineStyle','none');
                
        x1 = [-s1-1.5:1:s1+1.5];
        y1 = (s1/4).*x1./x1;
        z1 = 500.*x1./x1;
        plot3(x1,y1,z1,'linewidth',1,'color',colg);
        plot3(x1,2.*y1,z1,'linewidth',1,'color',colg);
        plot3(x1,3.*y1,z1,'linewidth',1,'color',colg);
        plot3(x1,-y1,z1,'linewidth',1,'color',colg);
        plot3(x1,-2.*y1,z1,'linewidth',1,'color',colg);
        plot3(x1,-3.*y1,z1,'linewidth',1,'color',colg);
        plot3(x1,0.*y1,z1,'linewidth',1,'color',colg);
        
        y1 = [-s1-1.5:1:s1+1.5];
        x1 = (s1/4).*y1./y1;
        z1 = 500.*y1./y1;
        plot3(x1,y1,z1,'linewidth',1,'color',colg);
        plot3(2.*x1,y1,z1,'linewidth',1,'color',colg);
        plot3(3.*x1,y1,z1,'linewidth',1,'color',colg);
        plot3(-x1,y1,z1,'linewidth',1,'color',colg);
        plot3(-2.*x1,y1,z1,'linewidth',1,'color',colg);
        plot3(-3.*x1,y1,z1,'linewidth',1,'color',colg);
        plot3(0.*x1,y1,z1,'linewidth',1,'color',colg);
        
%         if(clumps == 1)
%             filenamec = strcat('./',dirname1, ...
%                 '_non_spherical/binary_grid_outputs/face_on_surface_densities_Halpha_t',...
%                 filename(length(filename)-marker+1:length(filename)-4),'.bin')
%             fidc = fopen(filenamec, 'rb','ieee-be');
%             ntemp = fread(fidc,1, 'int')
%             nc = fread(fidc,2, 'int')
%             sc = fread(fidc,2, 'float64')
%             Rdc = fread(fidc,1, 'float32')
%             rawc = fread(fidc, nc(1) * nc(2), 'float32');
%             rtemp = fread(fidc,2, 'float32')
%             smooth1c = fread(fidc, nc(1) * nc(2), 'float32');
%             rtemp = fread(fidc,2, 'float32')
%             smooth2c = fread(fidc, nc(1) * nc(2), 'float32');
%             rtemp = fread(fidc,2, 'float32')
%             diffc = fread(fidc, nc(1) * nc(2), 'float32');
%             rtemp = fread(fidc,2, 'float32')
%             nclump = fread(fidc, 1, 'int');
%             idclump = fread(fidc, nclump, 'int');
%             spec = fread(fidc, nclump, 'int');
%             es = fread(fidc, nclump, 'int');
%             xclump = fread(fidc, nclump, 'float32');
%             yclump = fread(fidc, nclump, 'float32');
%             rclump = fread(fidc, nclump, 'float64');
%             mclump = fread(fidc, nclump, 'float64');
%             fclose(fidc);
%             clear rawc smooth1c smooth2c diffc
%             theta = linspace(0,2*pi,100);
%             for j=1:nclump
%                 if( log10(mclump(j)) >= 7 )
%                     if(log10(mclump(j)) >= 9)
%                         col = 'k';
%                     elseif(log10(mclump(j)) >= 8)
%                         col = 'r';
%                     else
%                         col = 'w';
%                     end
%                     if(es(j)==1)
%                         x1 = xclump(j) + rclump(j).*cos(theta);
%                         y1 = yclump(j) + rclump(j).*sin(theta);
%                         z1 = 50.*x1./x1;
%                         plot3(x1,y1,z1,'linewidth',2,'color',col);
%                     elseif(es(j)==0)
%                         x1 = linspace(xclump(j)-rclump(j),xclump(j)+rclump(j),10);
%                         y1 = (yclump(j)+rclump(j)).*x1./x1;
%                         z1=50.*x1./x1;
%                         x2 = linspace(xclump(j)-rclump(j),xclump(j)+rclump(j),10);
%                         y2 = (yclump(j)-rclump(j)).*x2./x2;
%                         z2=50.*x2./x2;
%                         y3 = linspace(yclump(j)-rclump(j),yclump(j)+rclump(j),10);
%                         x3 = (xclump(j)-rclump(j)).*y3./y3;
%                         z3=50.*x3./x3;
%                         y4 = linspace(yclump(j)-rclump(j),yclump(j)+rclump(j),10);
%                         x4 = (xclump(j)+rclump(j)).*y4./y4;
%                         z4=50.*x4./x4;
%                         plot3(x1,y1,z1,'linewidth',2,'color',col);
%                         plot3(x2,y2,z2,'linewidth',2,'color',col);
%                         plot3(x3,y3,z3,'linewidth',2,'color',col);
%                         plot3(x4,y4,z4,'linewidth',2,'color',col);
%                     end
%                 end
%             end
%             clear xclump yclump rclump spec es idclump nclump
%         end
%         
         set(gcf,'renderer','zbuffer')
         saveas(gcf,filename2);
         close all
         fclose all
        
        
%         Rad_d = 5;
%         x = linspace(-s1,s1,n);
%         y = x;
%         if(Rad_d<8)
%             d=1;
%         elseif(Rad_d==15 | Rad_d==18 | Rad_d==21)
%             d=3;
%         elseif(Rad_d<16)
%             d=2;
%         elseif(Rad_d>=16)
%             d=4;
%         end
% 
%         xlab = '$x\:[kpc]$';
%         ylab = '$y\:[kpc]$';
%         figure2 = figure('PaperType','a4letter','PaperSize',[20.98 29.68],...
%             'Colormap',cmap,'Visible','Off');
%         set(figure2,'WindowStyle','docked');
% 
% 	% Create axes
%         axes2 = axes('Parent',figure2,'YTick',[-Rad_d:d:Rad_d],...
%         'XTick',[-Rad_d:d:Rad_d],...
%         'PlotBoxAspectRatio',[1 1 1],...
%         'FontSize',16,...
%         'FontName','Times',...
%         'CLim',[cl(1) cl(2)]);
%         box(axes2,'on');
%         axis square
%          
% 	% Create textbox
%     annotation(figure1,'textbox',[0.1675 0.8605 0.12 0.064],...
%         'String',strcat('a=',filename(length(filename)-8:length(filename)-4)),...
%         'FontSize',16,...
%         'FontName','Times',...
%         'FitBoxToText','off',...
%         'LineStyle','none',...
%         'BackgroundColor',[1 1 1]);
%         filename2 = strcat(filename2(1:length(filename2)-9),'zoom_a',filename(length(filename)-8:length(filename)-4),'.jpg');
%     
%         ylabel(ylab,'Interpreter','latex','FontSize',16,'FontName','Times','units','normalized',...
%             'position',[-0.05,0.5,0.0]);
%         xlabel(xlab,'Interpreter','latex','FontSize',16,'FontName','Times','units','normalized',...
%             'position',[0.5,-0.065,0.0]);
%         bar=colorbar('peer',axes2,'FontSize',16,'FontName','Times',...
%             'CLim',[1 64],'YTick',cvec);
%         set(get(bar,'Ylabel'),'String',clab,'Interpreter','latex','Fontsize',16,...
%             'Rotation',270,'position',[10,(cl(1)+cl(2))/2,9.16],'FontName','Times');
%         xlim([-Rad_d Rad_d]);
%         ylim([-Rad_d Rad_d]);
%         hold('all');
%         [max(max(imag(face))),max(max(real(face))),max(max(imag(face)./real(face)))]
%         surf(x,y,face,'Parent',axes2,'LineStyle','none');
%         
% %         theta = linspace(0,2*pi,100);
% %         x1 = Rd.*cos(theta)./1000;
% %         y1 = Rd.*sin(theta)./1000;
% %         z1 = 500.*x1./x1;
% %         plot3(x1,y1,z1,'Parent',axes2,'linewidth',2,'color',colg,'linestyle','--');
% 
%         x1 = [-Rad_d-1.5:1:Rad_d+1.5];
%         y1 = (Rad_d/4).*x1./x1;
%         z1 = 500.*x1./x1;
%         plot3(x1,y1,z1,'Parent',axes2,'linewidth',1,'color',colg);
%         plot3(x1,2.*y1,z1,'Parent',axes2,'linewidth',1,'color',colg);
%         plot3(x1,3.*y1,z1,'Parent',axes2,'linewidth',1,'color',colg);
%         plot3(x1,-y1,z1,'Parent',axes2,'linewidth',1,'color',colg);
%         plot3(x1,-2.*y1,z1,'Parent',axes2,'linewidth',1,'color',colg);
%         plot3(x1,-3.*y1,z1,'Parent',axes2,'linewidth',1,'color',colg);
%         plot3(x1,0.*y1,z1,'Parent',axes2,'linewidth',1,'color',colg);
%         
%         y1 = [-Rad_d-1.5:1:Rad_d+1.5];
%         x1 = (Rad_d/4).*y1./y1;
%         z1 = 500.*y1./y1;
%         plot3(x1,y1,z1,'Parent',axes2,'linewidth',1,'color',colg);
%         plot3(2.*x1,y1,z1,'Parent',axes2,'linewidth',1,'color',colg);
%         plot3(3.*x1,y1,z1,'Parent',axes2,'linewidth',1,'color',colg);
%         plot3(-x1,y1,z1,'Parent',axes2,'linewidth',1,'color',colg);
%         plot3(-2.*x1,y1,z1,'Parent',axes2,'linewidth',1,'color',colg);
%         plot3(-3.*x1,y1,z1,'Parent',axes2,'linewidth',1,'color',colg);
%         plot3(0.*x1,y1,z1,'Parent',axes2,'linewidth',1,'color',colg);
% %         
% %         if(clumps == 1)
% %             filenamec = strcat('./',dirname1, ...
% %                 '_non_spherical/binary_grid_outputs/face_on_surface_densities_Halpha_t',...
% %                 filename(length(filename)-marker+1:length(filename)-4),'.bin')
% %             fidc = fopen(filenamec, 'rb','ieee-be');
% %             ntemp = fread(fidc,1, 'int')
% %             nc = fread(fidc,2, 'int')
% %             sc = fread(fidc,2, 'float64')
% %             Rdc = fread(fidc,1, 'float32')
% %             rawc = fread(fidc, nc(1) * nc(2), 'float32');
% %             rtemp = fread(fidc,2, 'float32')
% %             smooth1c = fread(fidc, nc(1) * nc(2), 'float32');
% %             rtemp = fread(fidc,2, 'float32')
% %             smooth2c = fread(fidc, nc(1) * nc(2), 'float32');
% %             rtemp = fread(fidc,2, 'float32')
% %             diffc = fread(fidc, nc(1) * nc(2), 'float32');
% %             rtemp = fread(fidc,2, 'float32')
% %             nclump = fread(fidc, 1, 'int');
% %             idclump = fread(fidc, nclump, 'int');
% %             spec = fread(fidc, nclump, 'int');
% %             es = fread(fidc, nclump, 'int');
% %             xclump = fread(fidc, nclump, 'float32');
% %             yclump = fread(fidc, nclump, 'float32');
% %             rclump = fread(fidc, nclump, 'float64');
% %             mclump = fread(fidc, nclump, 'float64');
% %             fclose(fidc);
% %             clear rawc smooth1c smooth2c diffc
% %             theta = linspace(0,2*pi,100);
% %             for j=1:nclump
% %                 if( log10(mclump(j)) >= 7 )
% %                     if(log10(mclump(j)) >= 9)
% %                         col = 'k';
% %                     elseif(log10(mclump(j)) >= 8)
% %                         col = 'r';
% %                     else
% %                         col = 'w';
% %                     end
% %                     if(es(j)==1)
% %                         x1 = xclump(j) + rclump(j).*cos(theta);
% %                         y1 = yclump(j) + rclump(j).*sin(theta);
% %                         z1 = 50.*x1./x1;
% %                         plot3(x1,y1,z1,'linewidth',2,'color',col);
% %                     elseif(es(j)==0)
% %                         x1 = linspace(xclump(j)-rclump(j),xclump(j)+rclump(j),10);
% %                         y1 = (yclump(j)+rclump(j)).*x1./x1;
% %                         z1=50.*x1./x1;
% %                         x2 = linspace(xclump(j)-rclump(j),xclump(j)+rclump(j),10);
% %                         y2 = (yclump(j)-rclump(j)).*x2./x2;
% %                         z2=50.*x2./x2;
% %                         y3 = linspace(yclump(j)-rclump(j),yclump(j)+rclump(j),10);
% %                         x3 = (xclump(j)-rclump(j)).*y3./y3;
% %                         z3=50.*x3./x3;
% %                         y4 = linspace(yclump(j)-rclump(j),yclump(j)+rclump(j),10);
% %                         x4 = (xclump(j)+rclump(j)).*y4./y4;
% %                         z4=50.*x4./x4;
% %                         plot3(x1,y1,z1,'linewidth',2,'color',col);
% %                         plot3(x2,y2,z2,'linewidth',2,'color',col);
% %                         plot3(x3,y3,z3,'linewidth',2,'color',col);
% %                         plot3(x4,y4,z4,'linewidth',2,'color',col);
% %                     end
% %                 end
% %             end
% %             clear xclump yclump rclump spec es idclump nclump
% %         end
% %         
%          set(gcf,'renderer','zbuffer')
%          saveas(gcf,filename2);
%          
%          clear face
%          close all;
%          fclose all;
     end
 end
