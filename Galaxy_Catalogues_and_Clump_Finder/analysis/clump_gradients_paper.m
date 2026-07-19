function clump_gradients_paper(in_situ, norm_in_situ, ex_situ, norm_ex_situ, nis, nes, ...
    prop, norm, gen, gal, combined, mass_min, norm_mass_min, norm_sfr_min, ...
    h_max, norm_lifetime_thresh, Md_thresh, z_thresh, include_zero_lifetime)

%(prop = 1: ID)
%prop = 2: redshift
%prop = 3: Rc
%prop = 4: M_{gas}
%prop = 5: M_{star}
%prop = 6: M_{bar}
%prop = 7: f_{gas}
%prop = 8: f_{dm}
%prop = 9: {\Sigma}_{gas}
%prop = 10: {\Sigma}_{*}
%prop = 11: {\Sigma}_{bar}
%prop = 12: age
%prop = 13: z_{gas}
%prop = 14: z_{stars}
%prop = 15: SFR
%prop = 16: {\Sigma}_{SFR}
%prop = 17: sSFR
%prop = 18: t_{dep}
%prop = 19: d/Rd
%prop = 20: z/Hd
%prop = 21: d
%prop = 22: z
%prop = 23: mean_residual
%prop = 24(,25,26): shape parameter
%prop = 27: dark matter contrast
%(prop = 28: Ex situ)
%prop = 29: tff
%prop = 30: tdyn_local
%prop = 31: tdyn_global
%prop = 41: life time
%prop = 63: eps_ff
%prop = 64: mass loading
%prop = 65: net in over out
%prop = 67: in over out
%prop = 100: t/td (global)

% Choose gal
if(gen==2)
    ind = [1:3, 5:16, 19, 21, 23, 25:35];
    gal_pref = 'VELA_';
%     tit = '${\rm No\:RP}\:\:';
elseif(gen==3)
    ind = [1:17, 19:35];
    gal_pref = 'VELA_v2_';
%     tit = '${\rm RP}\:\:';
end
tit = '$';
% if(Md_thresh(1)>7.5 | Md_thresh(2)<11)
%     tit = strcat(tit,num2str(Md_thresh(1),'%3.1f'),'<{\rm log}\:\:(\:M_{\rm d}\:\:)<',num2str(Md_thresh(2),'%3.1f'),'\:\:\:');
% end
% if(z_thresh(1)>1 | z_thresh(2)<4)
%     tit = strcat(tit,num2str(z_thresh(1),'%3.1f'),'<z<',num2str(z_thresh(2),'%3.1f'),'\:\:\:');
% end
if(gal>0)
    gal_pref = strcat(gal_pref,num2str(gal,'%02i'));
    bgal = find(ind==gal);
    if(nis(bgal+1)>nis(bgal))
        bis1 = nis(bgal)+1:nis(bgal+1);
        bes1 = nes(bgal)+1:nes(bgal+1);
    else
        bis1 = [];
        bes1 = [];
    end
else
    gal_pref = strcat(gal_pref,'stacked');
    bis1 = 1:length(in_situ);
    bes1 = 1:length(ex_situ);
end

if( ~isempty(bis1) | ~isempty(bes1) )
    if(prop==3)
        prop_list = 'Rc';
    elseif(prop==4)
        prop_list = 'Mgas';
    elseif(prop==5)
        prop_list = 'Mstar';
    elseif(prop==6)
        prop_list = 'Mbar';
    elseif(prop==7)
        prop_list = 'fgas';
    elseif(prop==8)
        prop_list = 'fdm';
    elseif(prop==9)
        prop_list = 'SigGas';
    elseif(prop==10)
        prop_list = 'SigStar';
    elseif(prop==11)
        prop_list = 'SigBar';
    elseif(prop==12)
        prop_list = 'AgeStar';
    elseif(prop==13)
        prop_list = 'ZGas';
    elseif(prop==14)
        prop_list = 'ZStar';
    elseif(prop==15)
        prop_list = 'SFR';
    elseif(prop==16)
        prop_list = 'SigSFR';
    elseif(prop==17)
        prop_list = 'sSFR';
    elseif(prop==18)
        prop_list = 'tdep';
    elseif(prop==19)
        prop_list = 'norm_dist';
    elseif(prop==20)
        prop_list = 'norm_height';
    elseif(prop==21)
        prop_list = 'dist';
    elseif(prop==22)
        prop_list = 'height';
    elseif(prop==23)
        prop_list = 'residual';
    elseif(prop==24)
        prop_list = 'shape1';
    elseif(prop==27)
        prop_list = 'DM_residual';
    elseif(prop==41)
        prop_list = 'life_time';
    elseif(prop==42)
        prop_list = 'td_tff_local';
    elseif(prop==43)
        prop_list = 'td_tff_global';
    elseif(prop==44)
        prop_list = 'Mgas_in';
    elseif(prop==45)
        prop_list = 'Mgas_out';
    elseif(prop==46)
        prop_list = 'n_gas';
    elseif(prop==47)
        prop_list = 'n_star';
    elseif(prop==48)
        prop_list = 'n_bar';
    elseif(prop==49)
        prop_list = 'max_mass';
    elseif(prop==63)
        prop_list = 'eps_ff';
    elseif(prop==64)
        prop_list = 'mass_loading';
    elseif(prop==65)
        prop_list = 'net_in_over_out';
    elseif(prop==67)
        prop_list = 'out_over_in';
    elseif(prop==100)
        prop_list = 't_over_td';
    end
    
    mkdir('./gradients_paper/');
    if(combined == 1)
        sample = 'common_sample';
    elseif(combined == 0)
        sample = 'total_sample';
    end
    % Normalized or absolute
    if(norm==1)
        is_dat = norm_in_situ(bis1,:);
        es_dat = norm_ex_situ(bes1,:);
    elseif(norm==0)
        is_dat = in_situ(bis1,:);
        es_dat = ex_situ(bes1,:);
    end
    norm_is_dat = norm_in_situ(bis1,:);
    norm_es_dat = norm_ex_situ(bes1,:);
    
    Mdstr = strcat(num2str(10*Md_thresh(1),'%3i'),'_Md_',num2str(10*Md_thresh(2),'%3i'));
    zstr = strcat(num2str(10*z_thresh(1),'%2i'),'_z_',num2str(10*z_thresh(2),'%2i'));
    Tstr = strcat('age_thresh_',num2str(norm_lifetime_thresh(1),'%1i'),'Tff');
    mkdir(strcat('./gradients_paper/',sample,'/'));
    mkdir(strcat('./gradients_paper/',sample,'/',gal_pref,'/'));
    mkdir(strcat('./gradients_paper/',sample,'/',gal_pref,'/',Tstr,'__',Mdstr,'__',zstr,'/'));
    mkdir(strcat('./gradients_paper/',sample,'/',gal_pref,'/',Tstr,'__',Mdstr,'__',zstr,'/',prop_list,'/'));
    filename = strcat('./gradients_paper/',sample,'/',gal_pref,'/',Tstr,'__',Mdstr,'__',zstr,'/',prop_list,'/');
    if(norm==1)
        filename = strcat('./gradients_paper/',sample,'/',gal_pref,'/',Tstr,'__',Mdstr,'__',zstr,'/',prop_list,'/Normalized__');
    end
    str = '';
    
    % normalized thresholds take priority
    if(norm_mass_min>-1e9)
        mass_min = -1e9;
        min_m = num2str(norm_mass_min,'%4.2f');
        str = strcat(str,'NM_a_',min_m);
        tit = strcat(tit,'log(M_{\rm c}/M_{\rm d})\:\:>',num2str(norm_mass_min,'%4.2f'),'\:\:\:');
    else
        if(mass_min>-1e9)
            min_m = num2str(mass_min,'%3.1f');
            tit = strcat(tit,'log(M_{\rm c})\:>',num2str(mass_min,'%3.1f'),'\:\:\:');
            str = strcat(str,'M_a_',min_m);
        end
    end
    
    if(norm_sfr_min>-1e9)
        norm_SFR_min = num2str(norm_sfr_min,'%4.2f');
        tit = strcat(tit,'SFR_{\rm c}\:/\:SFR_{\rm d}\:\:\:>',num2str(norm_sfr_min,'%4.2f'),'\:\:\:');
        if(~strcmp(str,''))
            str = strcat(str,'_');
        end
        str = strcat(str,'NS_a_',norm_SFR_min);
    end
    
    if(h_max<10)
        max_h = num2str(h_max,'%3.1f');
        tit = strcat(tit,num2str(h_max,'%3.1f'),'\:\:>h\:/\:H_{\rm d}');
        if(~strcmp(str,''))
            str = strcat(str,'_');
        end
        str = strcat(str,'H_b_',max_h);
    end
    if(strcmp(str,''))
        filename = strcat(filename,'all');
    else
        filename = strcat(filename,str);
    end
    tit = strcat(tit,'$');
   
    % Fix zero SFR
    b0 = find(is_dat(:,15)>1e-6);
    b  = find(is_dat(:,15)<=1e-6);
    if(~isempty(b))
        if(isempty(b0))
            b0 = b;
        end
        is_dat(b,15)      = min(is_dat(b0,15));
        is_dat(b,16)      = min(is_dat(b0,16));
        is_dat(b,17)      = min(is_dat(b0,17));
        is_dat(b,18)      = min(is_dat(b0,18));
        norm_is_dat(b,15) = min(norm_is_dat(b0,15));
        norm_is_dat(b,16) = min(norm_is_dat(b0,16));
        norm_is_dat(b,17) = min(norm_is_dat(b0,17));
        norm_is_dat(b,18) = min(norm_is_dat(b0,18));        
    end
    
    b0 = find(es_dat(:,15)>1e-6);
    b  = find(es_dat(:,15)<=1e-6);
    if(~isempty(b))
        if(isempty(b0))
            b0 = b;
        end
        es_dat(b,15)      = min(es_dat(b0,15));
        es_dat(b,16)      = min(es_dat(b0,16));
        es_dat(b,17)      = min(es_dat(b0,17));
        es_dat(b,18)      = min(es_dat(b0,18));
        norm_es_dat(b,15) = min(norm_es_dat(b0,15));
        norm_es_dat(b,16) = min(norm_es_dat(b0,16));
        norm_es_dat(b,17) = min(norm_es_dat(b0,17));
        norm_es_dat(b,18) = min(norm_es_dat(b0,18));
    end
    
    % Fix zero gas
    b0 = find(is_dat(:,4)>1e-6);
    b = find(is_dat(:,4)<=1e-6);
    if(~isempty(b))
        if(isempty(b0))
            b0 = b;
        end
        is_dat(b,4) = min(is_dat(b0,4));
        is_dat(b,7) = min(is_dat(b0,7));
        is_dat(b,9) = min(is_dat(b0,9));
        norm_is_dat(b,4) = min(norm_is_dat(b0,4));
        norm_is_dat(b,7) = min(norm_is_dat(b0,7));
        norm_is_dat(b,9) = min(norm_is_dat(b0,9));
    end
    b0 = find(es_dat(:,4)>1e-6);
    b = find(es_dat(:,4)<=1e-6);
    if(~isempty(b))
        if(isempty(b0))
            b0 = b;
        end
        es_dat(b,4) = min(es_dat(b0,4));
        es_dat(b,7) = min(es_dat(b0,7));
        es_dat(b,9) = min(es_dat(b0,9));
        norm_es_dat(b,4) = min(norm_es_dat(b0,4));
        norm_es_dat(b,7) = min(norm_es_dat(b0,7));
        norm_es_dat(b,9) = min(norm_es_dat(b0,9));
    end
    
    %%%% Thresholds
    b = find( log10(in_situ(bis1,6)./norm_in_situ(bis1,6))<=Md_thresh(2) & log10(in_situ(bis1,6)./norm_in_situ(bis1,6))>Md_thresh(1) & ...
        in_situ(bis1,2)<=z_thresh(2) & in_situ(bis1,2)>=z_thresh(1) & ...
        log10(in_situ(bis1,49))>=mass_min & ...
        log10(norm_in_situ(bis1,49))>=norm_mass_min & ...
        norm_in_situ(bis1,15)>=norm_sfr_min & ...
        in_situ(bis1,51)<=h_max);% &...
        %in_situ(bis1,23)>=20 );
        %abs(in_situ(bis1,42)./sqrt(3/(4*pi)))>1 );
    is_dat = is_dat(b,:);
    norm_is_dat = norm_is_dat(b,:);
    
    b = find( log10(ex_situ(bes1,6)./norm_ex_situ(bes1,6))<=Md_thresh(2) & log10(ex_situ(bes1,6)./norm_ex_situ(bes1,6))>Md_thresh(1) & ...
        ex_situ(bes1,2)<=z_thresh(2) & ex_situ(bes1,2)>=z_thresh(1) & ...
        log10(ex_situ(bes1,6))>=mass_min & ...
        log10(norm_ex_situ(bes1,6))>=norm_mass_min & ...
        norm_ex_situ(bes1,15)>=norm_sfr_min & ...
        ex_situ(bes1,51)<=h_max);% &...
        %ex_situ(bes1,23)>=20);
        %abs(ex_situ(bes1,42)./sqrt(3/(4*pi)))>1 );
    es_dat = es_dat(b,:);
    norm_es_dat = norm_es_dat(b,:);

    
    % Set things up
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    b=find(norm_is_dat(:,52)>=norm_lifetime_thresh);
    b1=find(norm_is_dat(:,52)<norm_lifetime_thresh & norm_is_dat(:,52)>1e-6);
    b2=find(norm_is_dat(:,52)<1e-6);
    
    if(prop==3)
        if(norm==0)
            Y2=log10(1000.*is_dat(b1,prop));    %SLCs
            Y1=log10(1000.*es_dat(:,prop));
            Y3=log10(1000.*is_dat(b,prop));     %LLCs
            Y4=log10(1000.*is_dat(b2,prop));     %ZLCs
            ylab1='${\rm log(}\:\:R_{\rm c}\:{\rm [pc])}$';
            ylab2='${\rm log(}R_{\rm c}\:{\rm [pc])}$';
            yu0=[3.2,3.2];
            yl0=[1.6,1.6];
            dy = 0.2;
        else
            Y2=log10(is_dat(b1,prop));
            Y1=log10(es_dat(:,prop));
            Y3=log10(is_dat(b,prop));
            Y4=log10(is_dat(b2,prop));
            ylab1='${\rm log(}\:\:R_{\rm c}\:/\:R_{\rm d}{\rm )}$';
            ylab2='${\rm log(}R_{\rm c}\:/\:R_{\rm d}{\rm )}$';
            yu0=[-0.6,-0.6];
            yl0=[-2.0,-2.0];
            dy = 0.2;
        end
    end
    
    if(prop==4)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        if(norm==0)
            ylab1='${\rm log(}\:\:M_{\rm gas,\:c}\:{\rm [M_{\odot}])}$';
            ylab2='${\rm log(}M_{\rm gas,\:c}\:{\rm [M_{\odot}])}$';
            yu0=[10,10];
            yl0=[6,6];
            dy = 0.5;
        else
            ylab1='${\rm log(}\:\:M_{\rm gas,\:c}\:/\:M_{\rm gas,\:d}{\rm )}$';
            ylab2='${\rm log(}M_{\rm gas,\:c}\:/\:M_{\rm gas,\:d}{\rm )}$';
            yu0=[0,0];
            yl0=[-4,-4];
            dy = 0.5;
        end
    end
    
    if(prop==5)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        if(norm==0)
            ylab1='${\rm log(}\:\:M_{\rm *,\:c}\:{\rm [M_{\odot}])}$';
            ylab2='${\rm log(}M_{\rm *,\:c}\:{\rm [M_{\odot}])}$';
            yu0=[10 10];
            yl0=[6 6];
            dy = 0.5;
        else
            ylab1='${\rm log(}\:\:M_{\rm *,\:c}\:/\:M_{\rm *,\:d}{\rm )}$';
            ylab2='${\rm log(}M_{\rm *,\:c}\:/\:M_{\rm *,\:d}{\rm )}$';
            yu0=[0,0];
            yl0=[-4,-4];
            dy = 0.5;
        end
    end
    
    if(prop==6)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        if(norm==0)
            ylab1='${\rm log(}\:\:M_{\rm c}\:{\rm [M_{\odot}])}$';
            ylab2='${\rm log(}M_{\rm c}\:{\rm [M_{\odot}])}$';
            yu0=[10,10];
            yl0=[6,6];
            dy = 0.5;
        else
            ylab1='${\rm log(}\:\:M_{\rm c}\:/\:M_{\rm d}{\rm )}$';
            ylab2='${\rm log(}M_{\rm c}\:/\:M_{\rm d}{\rm )}$';
            yu0=[-1, -1];
            yl0=[-4,-4];
            dy = 0.5;
        end
    end
    
    if(prop==7)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        if(norm==0)
            ylab1='${\rm log(}\:\:f_{\rm gas,\:c}\:\:{\rm )}$';
            ylab2='${\rm log(}f_{\rm gas,\:c}\:{\rm )}$';
            yu0=[0,0];
            yl0=[-2,-2];
            dy = 0.2;
        else
            ylab1='${\rm log(}\:\:f_{\rm gas,\:c}\:/\:f_{\rm gas,\:d}\:\:{\rm )}$';
            ylab2='${\rm log(}f_{\rm gas,\:c}\:/\:f_{\rm gas,\:d}\:{\rm )}$';
            yu0=[0.8,0.8];
            yl0=[-1.2,-1.2];
            dy = 0.2;
        end
    end
    
    if(prop==9)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        if(norm==0)
            ylab1='${\rm log(}\:\:{\Sigma}_{\rm gas,\:c}\:{\rm [M_{\odot}\:pc^{-2}])}$';
            ylab2='${\rm log(}{\Sigma}_{\rm gas,\:c}\:{\rm [M_{\odot}\:pc^{-2}])}$';
            yu0=[4,4];
            yl0=[-0.5,-0.5];
            dy = 0.5;
        else
            ylab1='${\rm log(}\:\:{\Sigma}_{\rm gas,\:c}\:/\:{\Sigma}_{\rm gas,\:d}{\rm )}$';
            ylab2='${\rm log(}{\Sigma}_{\rm gas,\:c}\:/\:{\Sigma}_{\rm gas,\:d}{\rm )}$';
            yu0=[3,3];
            yl0=[-1.5,-1.5];
            dy = 0.5;
        end
    end
    
    if(prop==10)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        if(norm==0)
            ylab1='${\rm log(}\:\:{\Sigma}_{\rm *,\:c}\:{\rm [M_{\odot}\:pc^{-2}])}$';
            ylab2='${\rm log(}{\Sigma}_{\rm *,\:c}\:{\rm [M_{\odot}\:pc^{-2}])}$';
            yu0=[4,4];
            yl0=[-0.5,-0.5];
            dy = 0.5;
        else
            ylab1='${\rm log(}\:\:{\Sigma}_{\rm *,\:c}\:/\:{\Sigma}_{\rm *,\:d}{\rm )}$';
            ylab2='${\rm log(}{\Sigma}_{\rm *,\:c}\:/\:{\Sigma}_{\rm *,\:d}{\rm )}$';
            yu0=[3,3];
            yl0=[-1.5,-1.5];
            dy = 0.5;
        end
    end
    
    if(prop==11)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        if(norm==0)
            ylab1='${\rm log(}\:\:{\Sigma}_{\rm c}\:{\rm [M_{\odot}\:pc^{-2}])}$';
            ylab2='${\rm log(}{\Sigma}_{\rm c}\:{\rm [M_{\odot}\:pc^{-2}])}$';
            yu0=[4,4];
            yl0=[-0.5,-0.5];
            dy = 0.5;
        else
            ylab1='${\rm log(}\:\:{\Sigma}_{\rm c}\:/\:{\Sigma}_{\rm d}{\rm )}$';
            ylab2='${\rm log(}{\Sigma}_{\rm c}\:/\:{\Sigma}_{\rm d}{\rm )}$';
            yu0=[3,3];
            yl0=[-1.5,-1.5];
            dy = 0.5;
        end
    end
    
    if(prop==12)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        if(norm==0)
            ylab1='${\rm log(}\:\:age_{\rm c}\:{\rm [Myr])}$';
            ylab2='${\rm log(}age_{\rm c}\:{\rm [Myr])}$';
            yu0=[4,4];
            yl0=[1,1];
            dy = 0.5;
        else
            ylab1='${\rm log(}\:\:age_{\rm c}\:/\:age_{\rm d}{\rm )}$';
            ylab2='${\rm log(}age_{\rm c}\:/\:age_{\rm d}{\rm )}$';
            yu0=[1,1];
            yl0=[-2,-2];
            dy = 0.5;
        end
    end
    
    if(prop==15)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        if(norm==0)
            ylab1='${\rm log(}\:\:SFR_{\rm c}\:{\rm [M_{\odot}\:yr^{-1}])}$';
            ylab2='${\rm log(}SFR_{\rm c}\:{\rm [M_{\odot}\:yr^{-1}])}$';
            yu0=[0,0];
            yl0=[-4,-4];
            dy = 0.5;
        else
            ylab1='${\rm log(}\:\:SFR_{\rm c}\:/\:SFR_{\rm d}{\rm )}$';
            ylab2='${\rm log(}SFR_{\rm c}\:/\:SFR_{\rm d}{\rm )}$';
            yu0=[0,0];
            yl0=[-4,-4];
            dy = 0.5;
        end
    end
    
    if(prop==16)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        if(norm==0)
            ylab1='${\rm log(}\:\:{\Sigma}_{\rm SFR,\:c}\:\:{\rm [M_{\odot}\:\:yr^{-1}\:\:kpc^{-2}]\:)}$';
            ylab2='${\rm log(}{\Sigma}_{\rm SFR,\:c}\:\:{\rm [M_{\odot}\:yr^{-1}\:kpc^{-2}])}$';
            yu0=[2.4,2.5];
            yl0=[-4,-4];
            dy = 0.5;
        else
            ylab1='${\rm log(}\:\:{\Sigma}_{\rm SFR,\:c}\:\:/\:{\Sigma}_{\rm SFR,\:d}\:\:{\rm )}$';
            ylab2='${\rm log(}{\Sigma}_{\rm SFR,\:c}\:/\:{\Sigma}_{\rm SFR,\:d}{\rm )}$';
            yu0=[2.5,2.5];
            yl0=[-1.5,-1.5];
            dy = 0.5;
        end
    end
    
    if(prop==17)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        if(norm==0)
            ylab1='${\rm log(}\:\:sSFR_{\rm c}\:{\rm [Gyr^{-1}])}$';
            ylab2='${\rm log(}sSFR_{\rm c}\:{\rm [Gyr^{-1}])}$';
            yu0=[2,2];
            yl0=[-2,-2];
            dy = 0.5;
        else
            ylab1='${\rm log(}\:\:sSFR_{\rm c}\:/\:sSFR_{\rm d}{\rm )}$';
            ylab2='${\rm log(}sSFR_{\rm c}\:/\:sSFR_{\rm d}{\rm )}$';
            yu0=[2,2];
            yl0=[-2,-2];
            dy = 0.5;
        end
    end
    
    if(prop==18)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        if(norm==0)
            ylab1='${\rm log(}\:\:t_{\rm dep,\:c}\:{\rm [Gyr])}$';
            ylab2='${\rm log(}t_{\rm dep,\:c}\:{\rm [Gyr])}$';
            yu0=[1.6,2];
            yl0=[-3.3,-3.5];
            dy = 0.5;
        else
            ylab1='${\rm log(}\:\:t_{\rm dep,\:c}\:/\:t_{\rm dep,\:d}\:\:{\rm )}$';
            ylab2='${\rm log(}t_{\rm dep,\:c}\:/\:t_{\rm dep,\:d}{\rm )}$';
            yu0=[2.5,2.5];
            yl0=[-3,-3];
            dy = 0.5;
        end
    end
    
    if(prop==13)
        if(norm==0)
            Y2=is_dat(b1,prop);
            Y1=es_dat(:,prop);
            Y3=is_dat(b,prop);
            Y4=is_dat(b2,prop);
            ylab1='${\rm log(O\:/\:H)_{gas}}\:\:\:+\:\:\:12$';
            ylab2='${\rm log(O\:/\:H)_{gas}}\:+\:12$';
            yu0=[9.4,9.4];
            yl0=[7.4,7.4];
            dy = 0.2;
        else
            Y2=log10(is_dat(b1,prop));
            Y1=log10(es_dat(:,prop));
            Y3=log10(is_dat(b,prop));
            Y4=log10(is_dat(b2,prop));
            ylab1='${\rm log}\:(\:Z_{\rm gas,\:c}\:\:/\:\:Z_{\rm gas,\:d}\:)$';
            ylab2='${\rm log}(Z_{\rm gas,\:c}\:/\:Z_{\rm gas,\:d})$';
            yu0=[1,1];
            yl0=[-1,-1];
            dy = 0.2;
        end
    end
    
    if(prop==14)
        Y2=is_dat(b1,prop);
        Y1=es_dat(:,prop);
        Y3=is_dat(b,prop);
        Y4=is_dat(b2,prop);
        if(norm==0)
            ylab1='${\rm log(O\:/\:H)_{stars}}\:\:\:+\:\:\:12$';
            ylab2='${\rm log(O\:/\:H)_{stars}}\:+\:12$';
            yu0=[9.4,9.4];
            yl0=[7.4,7.4];
            dy = 0.2;
        else
            Y2=log10(is_dat(b1,prop));
            Y1=log10(es_dat(:,prop));
            Y3=log10(is_dat(b,prop));
            Y4=log10(is_dat(b2,prop));
            ylab1='${\rm log}\:(\:Z_{\rm stars,\:c}\:\:/\:\:Z_{\rm stars,\:d}\:)$';
            ylab2='${\rm log}(Z_{\rm stars,\:c}\:/\:Z_{\rm stars,\:d})$';
            yu0=[2,2];
            yl0=[-2,-2];
            dy = 0.5;
        end
    end
    
    if(prop==20)
        if(norm==0)
            Y2=abs(is_dat(b1,prop));
            Y1=abs(es_dat(:,prop));
            Y3=abs(is_dat(b,prop));
            Y4=abs(is_dat(b2,prop));
            ylab1='$h\:/\:H_{\rm d}$';
            ylab2='$h\:/\:H_{\rm d}$';
            yu0=[10,10];
            yl0=[0,0];
            dy = 1;
        else
            Y2=log10(abs(is_dat(b1,prop)));
            Y1=log10(abs(es_dat(:,prop)));
            Y3=log10(abs(is_dat(b,prop)));
            Y4=log10(abs(is_dat(b2,prop)));
            ylab1='${\rm log(}\:\:h\:/\:H_{\rm d}{\rm )}$';
            ylab2='${\rm log(}h\:/\:H_{\rm d}{\rm )}$';
            yu0=[2,2];
            yl0=[-2,-2];
            dy = 0.5;
        end
    end
    
    if(prop==23)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        ylab1='${\rm log}(\:\:{\bar {\delta}}_{\rm c})$';
        ylab2='${\rm log}{\bar {\delta}}_{\rm c})$';
        yu0=[4 4];
        yl0=[0.5 0.5];
        dy = 0.5;
    end
    
    if(prop==24)
        Y2=is_dat(b1,prop);
        Y1=es_dat(:,prop);
        Y3=is_dat(b,prop);
        Y4=is_dat(b2,prop);
        ylab1='$S_{\rm c}$';
        ylab2='$S_{\rm c}$';
        yu0=[1 1];
        yl0=[0 0];
        dy = 0.1;
    end
    
    if(prop==41)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        if(norm==0)
            ylab1='${\rm log(}\:\:t_{\rm c}\:{\rm [Myr])}$';
            ylab2='${\rm log(}t_{\rm c}\:{\rm [Myr])}$';
            yu0=[3.7,4];
            yl0=[1,1];
            dy = 0.5;
        else
            ylab1='${\rm log(}\:\:t_{\rm c}\:/\:t_{\rm ff,\:c}{\rm )}$';
            ylab2='${\rm log(}t_{\rm c}\:/\:t_{\rm ff,\:c}{\rm )}$';
            yu0=[3,3];
            yl0=[-3,-3];
            dy = 0.5;
        end
    end    
    
    if(prop==42)
        Y2=log10( is_dat(b1,prop)./sqrt(3/(4*pi)) );
        Y1=log10( es_dat(:,prop)./sqrt(3/(4*pi)) );
        Y3=log10( is_dat(b,prop)./sqrt(3/(4*pi)) );
        Y4=log10( is_dat(b2,prop)./sqrt(3/(4*pi)) );
        ylab1='${\rm log(}\:\:t_{\rm d,\,l}\:/\:t_{\rm ff,\:c}\:)$';
        ylab2='${\rm log(}t_{\rm d,\,l}\:/\:t_{\rm ff,\:c})$';
        yu0=[2,2];
        yl0=[-2,-2];
        dy = 0.5;
    end
    
    if(prop==43)
        Y2=log10(abs( is_dat(b1,prop)./sqrt(3/(4*pi)) ));
        Y1=log10(abs( es_dat(:,prop)./sqrt(3/(4*pi)) ));
        Y3=log10(abs( is_dat(b,prop)./sqrt(3/(4*pi)) ));
        Y4=log10(abs( is_dat(b2,prop)./sqrt(3/(4*pi)) ));
        ylab1='${\rm log(}\:\:t_{\rm d,\,g}\:/\:t_{\rm ff,\:c}\:)$';
        ylab2='${\rm log(}t_{\rm d,\,g}\:/\:t_{\rm ff,\:c})$';
        yu0=[2,2];
        yl0=[-2,-2];
        dy = 0.5;
    end
    
    if(prop==46)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        if(norm==0)
            ylab1='${\rm log(}\:\:n_{\rm gas,\:c}\:\:{\rm [cm^{-3}]\:)}$';
            ylab2='${\rm log(}n_{\rm gas,\:c}\:{\rm [cm^{-3}])}$';
            yu0=[3,3];
            yl0=[-1,-1];
            dy = 0.5;
        else
            ylab1='${\rm log(}\:\:n_{\rm gas,\:c}\:\:/\:n_{\rm gas,d}\:{\rm )}$';
            ylab2='${\rm log(}n_{\rm gas,\:c}\:/\:n_{\rm gas,d}{\rm )}$';
            yu0=[4,4];
            yl0=[0,0];
            dy = 0.5;
        end
    end
    
    if(prop==47)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        if(norm==0)
            ylab1='${\rm log(}\:\:n_{\rm *,\:c}\:\:{\rm [cm^{-3}]\:)}$';
            ylab2='${\rm log(}n_{\rm *,\:c}\:{\rm [cm^{-3}]\:)}$';
            yu0=[3,3];
            yl0=[-1,-1];
            dy = 0.5;
        else
            ylab1='${\rm log(}\:\:n_{\rm *,\:c}\:\:{\rm [cm^{-3}]\:)}$';
            ylab2='${\rm log(}n_{\rm *,\:c}\:{\rm [cm^{-3}])}$';
            yu0=[4,4];
            yl0=[0,0];
            dy = 0.5;
        end
    end
    
    if(prop==48)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        if(norm==0)
            ylab1='${\rm log(}\:\:n_{\rm bar,\:c}\:\:{\rm [cm^{-3}]\:)}$';
            ylab2='${\rm log(}\:n_{\rm bar,\:c}\:{\rm [cm^{-3}])}$';
            yu0=[3,3];
            yl0=[-1,-1];
            dy = 0.5;
        else
            ylab1='${\rm log(}\:\:n_{\rm bar,\:c}\:\:/\:n_{\rm bar,d}\:{\rm )}$';
            ylab2='${\rm log(}n_{\rm bar,\:c}\:/\:n_{\rm bar,d}{\rm )}$';
            yu0=[4,4];
            yl0=[0,0];
            dy = 0.5;
        end
    end
    
    if(prop==49)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        if(norm==0)
            ylab1='${\rm log(}\:\:{\rm max}\:\:(M_{\rm c})\:{\rm [M_{\odot}])}$';
            ylab2='${\rm log(}{\rm max}\:(M_{\rm c})\:{\rm [M_{\odot}])}$';
            yu0=[10,10];
            yl0=[5.5,5.5];
            dy = 0.5;
        else
            ylab1='${\rm log(}\:\:{\rm max}\:\:(M_{\rm c}\:/\:M_{\rm d}\:){\rm )}$';
            ylab2='${\rm log(}{\rm max}\:(M_{\rm c}\:/\:M_{\rm d}){\rm )}$';
            yu0=[0.5,0.5];
            yl0=[-4.5,-4.5];
            dy = 0.5;
        end
    end    
    
    if(prop==63)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        ylab1='${\rm log}\:(\:\epsilon_{\rm ff,\,c}\:)$';
        ylab2='${\rm log}(\epsilon_{\rm ff,\,c})$';
        yu0=[0.5 0.5];
        yl0=[-4 -4];
        dy = 0.5;
    end
    
    if(prop==64)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        ylab1='${\rm log}\:(\:\eta_{\rm c}\:)$';
        ylab2='${\rm log}(\eta_{\rm c})$';
        yu0=[1.5 1.5];
        yl0=[-1.5 -1.5];
        dy = 0.3;
    end
    
    if(prop==65)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        ylab1='${\rm log}\:[\:{\dot {M}}_{\rm in}\:/\:({\dot {M}}_{\rm out}+SFR)]$';
        ylab2='${\rm log}[{\dot {M}}_{\rm in}\:/\:({\dot {M}}_{\rm out}+SFR)]$';
        yu0=[1.5 1.5];
        yl0=[-1.5 -1.5];
        dy = 0.3;
    end    
    
    if(prop==67)
        Y2=log10(is_dat(b1,prop));
        Y1=log10(es_dat(:,prop));
        Y3=log10(is_dat(b,prop));
        Y4=log10(is_dat(b2,prop));
        ylab1='${\rm log}\:(\:{\dot {M}}_{\rm out}\:/\:{\dot {M}}_{\rm in})$';
        ylab2='${\rm log}({\dot {M}}_{\rm out}\:/\:{\dot {M}}_{\rm in})$';
        yu0=[1.5 1.5];
        yl0=[-1.5 -1.5];
        dy = 0.3;
    end
    
    if(prop==100)
        Y2=log10( is_dat(b1,41)./abs(is_dat(b1,31)) );
        Y1=log10( es_dat(:,41)./abs(es_dat(:,31)) );
        Y3=log10( is_dat(b,41)./abs(is_dat(b,31)) );
        Y4=log10( is_dat(b2,41)./abs(is_dat(b2,31)) );
        ylab1='${\rm log(}\:\:t_{\rm c}\:/\:t_{\rm dyn}\:)$';
        ylab2='${\rm log(}t_{\rm c}\:/\:t_{\rm dyn})$';
        yu0=[2.5,2.5];
        yl0=[-0.5,-0.5];
        dy = 0.5;
    end
    
    if(prop==6 | prop==7 | prop==12)
        yl0(2) = yl0(1)+dy;
    end
    
    if(isempty(Y1))
        Y1 = 1e-6;
    end
    if(isempty(Y2))
        Y2 = 1e-6;
    end
    if(isempty(Y3))
        Y3 = 1e-6;
    end
    
    X1=log10(es_dat(:,19));    % ex situ clumps
    X2=log10(is_dat(b1,19));   % SLC
    X3=log10(is_dat(b,19));    % LLC
    X4=log10(is_dat(b2,19));   % ZLC
    xlab1='${\rm log(}\:\:\:d\:/\:R_{\rm d}{\rm )}$';
    xlab2='${\rm log(}d\:/\:R_{\rm d}{\rm )}$';
    xl0 = [-1.4 -1.4];
    xu0 = [0.6 0.6];
    dx = 0.2;
    if(isempty(X1))
        X1 = 1e-6;
    end
    if(isempty(X2))
        X2 = 1e-6;
    end
    if(isempty(X3))
        X3 = 1e-6;
    end
    if(isempty(X4))
        X4 = 1e-6;
    end
    dbin = 0.15;
    leg_pos = [0.07 0.17 0.2 0.12];
    leg_pos2 = [0.15 0.155 0.2 0.07];
    leg_pos3 = [7.00 3.30 0.2 0.03];
    
        %%% Only All Medians and 1-sigma scatter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%     nice_fig(xl0,xu0,yl0,yu0,dx,dy,xlab1,ylab1,tit,16,20,14,1,0,1)    
%     
%     xleft  = -2:dbin:1;
%     xright = xleft+dbin;
% %     N_per_bin = 30;
% %     xsort = sort(X3);
% %     Nsort = ceil(length(xsort)/N_per_bin);
% %     xleft = zeros(1,Nsort);
% %     xright = zeros(1,Nsort);
% %     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
% %     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
% 
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X1<=xright(i) & X1>xleft(i));   % ESCs
%         if(prop==20)
%             is4=sort(abs(Y1(b4)));
%         else
%             is4=sort(Y1(b4));
%         end
%         ycen(i) = median(is4);
%         if(length(is4)>=7)
%             nl(i) = ceil(length(is4)/6);
%             nu(i) = ceil(5*length(is4)/6);            
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(4) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','b',...
%         'marker','^','markerfacecolor','b','markeredgecolor','b','markersize',6,'DisplayName','$ex\:situ$');
%     plot(xcen(b4),yu(b4),'linestyle','--','linewidth',2,'color','b','marker','^',...
%         'markerfacecolor','b','markeredgecolor','b','markersize',4);
%     plot(xcen(b4),yl(b4),'linestyle','--','linewidth',2,'color','b','marker','^',...
%         'markerfacecolor','b','markeredgecolor','b','markersize',4);
%     
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X4<=xright(i) & X4>xleft(i));   % ZLCs
%         if(prop==20)
%             is4=sort(abs(Y4(b4)));
%         else
%             is4=sort(Y4(b4));
%         end
%         ycen(i) = median(is4);
%         if(length(is4)>=7)
%             nl(i) = ceil(length(is4)/6);
%             nu(i) = ceil(5*length(is4)/6);            
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(3) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','g',...
%         'marker','d','markerfacecolor','g','markeredgecolor','g','markersize',6,'DisplayName','$ZLCs$');
%     plot(xcen(b4),yu(b4),'linestyle','--','linewidth',2,'color','g','marker','d',...
%         'markerfacecolor','g','markeredgecolor','g','markersize',4);
%     plot(xcen(b4),yl(b4),'linestyle','--','linewidth',2,'color','g','marker','d',...
%         'markerfacecolor','g','markeredgecolor','g','markersize',4);
%     
%     xleft  = -2:dbin:1;
%     xright = xleft+dbin;
% %     N_per_bin = 30;
% %     xsort = sort(X2);
% %     Nsort = ceil(length(xsort)/N_per_bin);
% %     xleft = zeros(1,Nsort);
% %     xright = zeros(1,Nsort);
% %     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
% %     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
%     
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X2<=xright(i) & X2>xleft(i));   % SLCs
%         if(prop==20)
%             is4=sort(abs(Y2(b4)));
%         else
%             is4=sort(Y2(b4));
%         end
%         ycen(i) = median(is4);
%         if(length(is4)>=7)
%             nl(i) = ceil(length(is4)/6);
%             nu(i) = ceil(5*length(is4)/6);            
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(2) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','r',...
%         'marker','s','markerfacecolor','r','markeredgecolor','r','markersize',6,'DisplayName','$SLCs$');
%     plot(xcen(b4),yu(b4),'linestyle','--','linewidth',2,'color','r','marker','s',...
%         'markerfacecolor','r','markeredgecolor','r','markersize',4);
%     plot(xcen(b4),yl(b4),'linestyle','--','linewidth',2,'color','r','marker','s',...
%         'markerfacecolor','r','markeredgecolor','r','markersize',4);
%         
%     xleft  = -2:dbin:1;
%     xright = xleft+dbin;
% %     N_per_bin = 30;
% %     xsort = sort(X3);
% %     Nsort = ceil(length(xsort)/N_per_bin);
% %     xleft = zeros(1,Nsort);
% %     xright = zeros(1,Nsort);
% %     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
% %     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
%     
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X3<=xright(i) & X3>xleft(i));   % LLCs
%         if(prop==20)
%             is4=sort(abs(Y3(b4)));
%         else
%             is4=sort(Y3(b4));
%         end
%         ycen(i) = median(is4);
%         if(length(is4)>=7)
%             nl(i) = ceil(length(is4)/6);
%             nu(i) = ceil(5*length(is4)/6);            
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(1) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','k',...
%         'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',6,'DisplayName','$LLCs$');
%     plot(xcen(b4),yu(b4),'linestyle','--','linewidth',2,'color','k','marker','o',...
%         'markerfacecolor','k','markeredgecolor','k','markersize',4);
%     plot(xcen(b4),yl(b4),'linestyle','--','linewidth',2,'color','k','marker','o',...
%         'markerfacecolor','k','markeredgecolor','k','markersize',4);
%     
%     if(prop==13 | prop==14)
%         % Solar metallicity
%         x=[-2:0.3:4];
%         y=8.69.*x./x;
%         plot(x,y,'Parent',gca,'marker','none','color','k','linewidth',2,'linestyle','-.')
%     end
%     
%     % Create legend    
%     legend1 = legend(gca,l(1:4));
%     set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
%         'Location','SouthWest','FontSize',14,...
%         'Interpreter','latex');
%     legend boxoff
%     
%     filename1 = strcat(filename,'_all_medians.jpg');
%     print(gcf,'-djpeg',filename1);
%     
%     xlabel(xlab2,'Interpreter','latex','FontSize',24,...
%     'FontName','Times New Roman','units','normalized',...
%     'position',[0.5 -0.06 0]);
%     ylabel(ylab2,'Interpreter','latex','FontSize',24,...
%         'FontName','Times New Roman','units','normalized',...
%         'position',[-0.08 0.5 0]);
%     title('');
%     set(legend1,'FontSize',18,'units','normalized','position',leg_pos);
%     filename1 = strcat(filename,'_all_medians.eps');
%     print(gcf,'-depsc',filename1);
%     close all
%     fclose all
    
            %%% LLC, SLC, ESC points + all medians
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    nice_fig(xl0,xu0,yl0,yu0,dx,dy,xlab1,ylab1,tit,16,20,14,1,0,1)  
    % Create plot
    ll(3) = plot(X1,Y1,'Parent',gca,'MarkerFaceColor',[0.2 0.8 0.8],...
        'MarkerEdgeColor',[0.2 0.8 0.8],...
        'MarkerSize',5,...
        'Marker','^',...
        'LineStyle','none', 'DisplayName',' '); 
    ll(2) = plot(X2,Y2,'Parent',gca,'MarkerFaceColor',[0.8 0.8 0.2],...
        'MarkerEdgeColor',[0.8 0.8 0.2],...
        'MarkerSize',5,...
        'Marker','s',...
        'LineStyle','none', 'DisplayName',' ');  
    ll(1) = plot(X3,Y3,'Parent',gca,'MarkerFaceColor',[0.5 0.5 0.5],...
        'MarkerEdgeColor',[0.5 0.5 0.5],...
        'MarkerSize',5,...
        'Marker','o',...
        'LineStyle','none', 'DisplayName',' ');
    
    xleft  = -2:dbin:1;
    xright = xleft+dbin;
%     N_per_bin = 30;
%     xsort = sort(X3);
%     Nsort = ceil(length(xsort)/N_per_bin);
%     xleft = zeros(1,Nsort);
%     xright = zeros(1,Nsort);
%     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
%     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];

    xcen   = (xleft+xright)/2;
    ycen   = -1e9.*ones(size(xcen));
    nl   = ycen;
    nu   = ycen;
    yl   = ycen;
    yu   = ycen;
    for i=1:length(xcen)
        b4=find(X1<=xright(i) & X1>xleft(i));   % ESCs
        if(prop==20)
            es4=sort(abs(Y1(b4)));
        else
            es4=sort(Y1(b4));
        end
        ycen(i) = median(es4);
        if(length(es4)>=5)
            nl(i) = ceil(length(es4)/6);
            nu(i) = ceil(5*length(es4)/6);            
            yu(i)=es4(nu(i));
            yl(i)=es4(nl(i));
        end
    end
    b4 = find(yu ~= -1e9);
    l(3) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','b',...
        'marker','^','markerfacecolor','b','markeredgecolor','b','markersize',8,'DisplayName','$ex\:situ$');
%     plot(xcen(b4),yu(b4),'linestyle','--','linewidth',2,'color','b','marker','^',...
%         'markerfacecolor','b','markeredgecolor','b','markersize',4);
%     plot(xcen(b4),yl(b4),'linestyle','--','linewidth',2,'color','b','marker','^',...
%         'markerfacecolor','b','markeredgecolor','b','markersize',4);
    
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X4<=xright(i) & X4>xleft(i));   % ZLCs
%         if(prop==20)
%             is4=sort(abs(Y4(b4)));
%         else
%             is4=sort(Y4(b4));
%         end
%         ycen(i) = median(is4);
%         if(length(is4)>=7)
%             nl(i) = ceil(length(is4)/6);
%             nu(i) = ceil(5*length(is4)/6);            
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(3) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','g',...
%         'marker','d','markerfacecolor','g','markeredgecolor','g','markersize',6,'DisplayName','$ZLCs$');
%     plot(xcen(b4),yu(b4),'linestyle','--','linewidth',2,'color','g','marker','d',...
%         'markerfacecolor','g','markeredgecolor','g','markersize',4);
%     plot(xcen(b4),yl(b4),'linestyle','--','linewidth',2,'color','g','marker','d',...
%         'markerfacecolor','g','markeredgecolor','g','markersize',4);
%     
    xleft  = -2:dbin:1;
    xright = xleft+dbin;
%     N_per_bin = 30;
%     xsort = sort(X2);
%     Nsort = ceil(length(xsort)/N_per_bin);
%     xleft = zeros(1,Nsort);
%     xright = zeros(1,Nsort);
%     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
%     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
    
    xcen   = (xleft+xright)/2;
    ycen   = -1e9.*ones(size(xcen));
    nl   = ycen;
    nu   = ycen;
    yl   = ycen;
    yu   = ycen;
    for i=1:length(xcen)
        b4=find(X2<=xright(i) & X2>xleft(i));   % SLCs
        if(prop==20)
            is4=sort(abs(Y2(b4)));
        else
            is4=sort(Y2(b4));
        end
        ycen(i) = median(is4);
        if(length(is4)>=7)
            nl(i) = ceil(length(is4)/6);
            nu(i) = ceil(5*length(is4)/6);            
            yu(i)=is4(nu(i));
            yl(i)=is4(nl(i));
        end
    end
    b4 = find(yu ~= -1e9);
    l(2) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','r',...
        'marker','s','markerfacecolor','r','markeredgecolor','r','markersize',8,'DisplayName','$SLCs$');
%     plot(xcen(b4),yu(b4),'linestyle','--','linewidth',2,'color','r','marker','s',...
%         'markerfacecolor','r','markeredgecolor','r','markersize',4);
%     plot(xcen(b4),yl(b4),'linestyle','--','linewidth',2,'color','r','marker','s',...
%         'markerfacecolor','r','markeredgecolor','r','markersize',4);
        
    xleft  = -2:dbin:1;
    xright = xleft+dbin;
%     N_per_bin = 30;
%     xsort = sort(X3);
%     Nsort = ceil(length(xsort)/N_per_bin);
%     xleft = zeros(1,Nsort);
%     xright = zeros(1,Nsort);
%     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
%     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
    
    xcen   = (xleft+xright)/2;
    ycen   = -1e9.*ones(size(xcen));
    nl   = ycen;
    nu   = ycen;
    yl   = ycen;
    yu   = ycen;
    for i=1:length(xcen)
        b4=find(X3<=xright(i) & X3>xleft(i));   % LLCs
        if(prop==20)
            is4=sort(abs(Y3(b4)));
        else
            is4=sort(Y3(b4));
        end
        ycen(i) = median(is4);
        if(length(is4)>=7)
            nl(i) = ceil(length(is4)/6);
            nu(i) = ceil(5*length(is4)/6);            
            yu(i)=is4(nu(i));
            yl(i)=is4(nl(i));
        end
    end
    b4 = find(yu ~= -1e9);
    l(1) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','k',...
        'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',8,'DisplayName','$LLCs$');
%     plot(xcen(b4),yu(b4),'linestyle','--','linewidth',2,'color','k','marker','o',...
%         'markerfacecolor','k','markeredgecolor','k','markersize',4);
%     plot(xcen(b4),yl(b4),'linestyle','--','linewidth',2,'color','k','marker','o',...
%         'markerfacecolor','k','markeredgecolor','k','markersize',4);
    
    if(prop==13 | prop==14)
        % Solar metallicity
        x=[-2:0.3:4];
        y=8.69.*x./x;
        plot(x,y,'Parent',gca,'marker','none','color','k','linewidth',2,'linestyle','-.')
    end
    
    % Create legend    
    legend1 = legend(gca,l(1:3));
    set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
        'Location','SouthWest','FontSize',14,...
        'Interpreter','latex');
    legend boxoff
    
    filename1 = strcat(filename,'_data_points.jpg');
    print(gcf,'-djpeg',filename1);
    
    xlabel(xlab2,'Interpreter','latex','FontSize',24,...
    'FontName','Times New Roman','units','normalized',...
    'position',[0.5 -0.06 0]);
    ylabel(ylab2,'Interpreter','latex','FontSize',24,...
        'FontName','Times New Roman','units','normalized',...
        'position',[-0.08 0.5 0]);
    title('');
    set(legend1,'FontSize',18,'units','normalized','position',leg_pos2);
    
    ah2=axes('position',get(gca,'position'), 'visible','off');
    ti = get(gca,'TightInset');    
    set(gca,'Position',[ti(1)+0.05 ti(2)+0.08 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3)+1.6 pos(4)+ti(2)+ti(4)+3.2]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3)+1.6 pos(4)+ti(2)+ti(4)+2.0]);    
    legend2 = legend(ah2,ll(1:3));
    legend boxoff
    set(legend2,'edgecolor','w','Position',leg_pos3,...
        'fontname','Times New Roman','fontsize',18,'interpreter','latex');
    
    filename1 = strcat(filename,'_data_points.eps');
    print(gcf,'-depsc',filename1);
    close all
    fclose all
    
    %%% Long Lived
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%     nice_fig(xl0,xu0,yl0,yu0,dx,dy,xlab1,ylab1,tit,16,20,14,1,0,1)    
%     % Create plot
%     l(1)=plot(X3,Y3,'Parent',gca,'MarkerFaceColor',[0.5 0.5 0.5],...
%         'MarkerEdgeColor',[0.5 0.5 0.5],...
%         'MarkerSize',3,...
%         'Marker','o',...
%         'LineStyle','none');
%     xleft  = -2:dbin:1;
%     xright = xleft+dbin;
% %     N_per_bin = 30;
% %     xsort = sort(X3);
% %     Nsort = ceil(length(xsort)/N_per_bin);
% %     xleft = zeros(1,Nsort);
% %     xright = zeros(1,Nsort);
% %     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
% %     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
%     
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X4<=xright(i) & X4>xleft(i));   % ZLCs
%         if(prop==20)
%             is4=sort(abs(Y4(b4)));
%         else
%             is4=sort(Y4(b4));
%         end
%         ycen(i) = median(is4);
%         nl(i) = ceil(length(is4)/6);
%         nu(i) = ceil(5*length(is4)/6);
%         if(nl(i)>0 & nu(i)>nl(i))
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(4) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','g',...
%         'marker','none','DisplayName','$ZLCs$');
%     
%     xleft  = -2:dbin:1;
%     xright = xleft+dbin;
% %     N_per_bin = 30;
% %     xsort = sort(X2);
% %     Nsort = ceil(length(xsort)/N_per_bin);
% %     xleft = zeros(1,Nsort);
% %     xright = zeros(1,Nsort);
% %     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
% %     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
%     
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X2<=xright(i) & X2>xleft(i));   % SLCs
%         if(prop==20)
%             is4=sort(abs(Y2(b4)));
%         else
%             is4=sort(Y2(b4));
%         end
%         ycen(i) = median(is4);
%         nl(i) = ceil(length(is4)/6);
%         nu(i) = ceil(5*length(is4)/6);
%         if(nl(i)>0 & nu(i)>nl(i))
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(3) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','r',...
%         'marker','none','DisplayName','$SLCs$');
%         
%     xleft  = -2:dbin:1;
%     xright = xleft+dbin;
% %     N_per_bin = 30;
% %     xsort = sort(X3);
% %     Nsort = ceil(length(xsort)/N_per_bin);
% %     xleft = zeros(1,Nsort);
% %     xright = zeros(1,Nsort);
% %     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
% %     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
%     
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X3<=xright(i) & X3>xleft(i));   % LLCs
%         if(prop==20)
%             is4=sort(abs(Y3(b4)));
%         else
%             is4=sort(Y3(b4));
%         end
%         ycen(i) = median(is4);
%         nl(i) = ceil(length(is4)/6);
%         nu(i) = ceil(5*length(is4)/6);
%         if(nl(i)>0 & nu(i)>nl(i))
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(2) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','k',...
%         'marker','none','DisplayName','$LLCs$');
% 
%     if(prop==13 | prop==14)
%         % Solar metallicity
%         x=[-2:0.3:4];
%         y=8.69.*x./x;
%         plot(x,y,'Parent',gca,'marker','none','color','k','linewidth',2,'linestyle','-.')
%     end
%     
%     % Create legend    
%     legend1 = legend(gca,l(2:4));
%     set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
%         'Location','SouthWest','FontSize',14,...
%         'Interpreter','latex');
%     legend boxoff
%     
%     filename1 = strcat(filename,'_LLC.jpg');
%     print(gcf,'-djpeg',filename1);
%     
%     xlabel(xlab2,'Interpreter','latex','FontSize',24,...
%         'FontName','Times New Roman','units','normalized',...
%         'position',[0.5 -0.07 0]);
%     ylabel(ylab2,'Interpreter','latex','FontSize',24,...
%         'FontName','Times New Roman','units','normalized',...
%         'position',[-0.09 0.5 0]);
%     title('');
%     set(legend1,'FontSize',18,'Location','SouthWest');
%     filename1 = strcat(filename,'_LLC.eps');
%     print(gcf,'-depsc',filename1);
%     close all
%     fclose all
%     
%     %%% Short Lived
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%     nice_fig(xl0,xu0,yl0,yu0,dx,dy,xlab1,ylab1,tit,16,20,14,1,0,1)    
%     % Create plot
%     l(1)=plot(X2,Y2,'Parent',gca,'MarkerFaceColor',[0.8706 0.4902 0],...
%         'MarkerEdgeColor',[0.8706 0.4902 0],...
%         'MarkerSize',3,...
%         'Marker','o',...
%         'LineStyle','none');
%     
%     xleft  = -2:dbin:1;
%     xright = xleft+dbin;
% %     N_per_bin = 30;
% %     xsort = sort(X3);
% %     Nsort = ceil(length(xsort)/N_per_bin);
% %     xleft = zeros(1,Nsort);
% %     xright = zeros(1,Nsort);
% %     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
% %     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
%     
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X4<=xright(i) & X4>xleft(i));   % ZLCs
%         if(prop==20)
%             is4=sort(abs(Y4(b4)));
%         else
%             is4=sort(Y4(b4));
%         end
%         ycen(i) = median(is4);
%         if(length(is4)>=7)
%             nl(i) = ceil(length(is4)/6);
%             nu(i) = ceil(5*length(is4)/6);            
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(4) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','g',...
%         'marker','none','DisplayName','$ZLCs$');
%     
%     xleft  = -2:dbin:1;
%     xright = xleft+dbin;
% %     N_per_bin = 30;
% %     xsort = sort(X2);
% %     Nsort = ceil(length(xsort)/N_per_bin);
% %     xleft = zeros(1,Nsort);
% %     xright = zeros(1,Nsort);
% %     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
% %     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
%     
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X2<=xright(i) & X2>xleft(i));   % SLCs
%         if(prop==20)
%             is4=sort(abs(Y2(b4)));
%         else
%             is4=sort(Y2(b4));
%         end
%         ycen(i) = median(is4);
%         if(length(is4)>=7)
%             nl(i) = ceil(length(is4)/6);
%             nu(i) = ceil(5*length(is4)/6);            
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(3) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','r',...
%         'marker','none','DisplayName','$SLCs$');
%         
%     xleft  = -2:dbin:1;
%     xright = xleft+dbin;
% %     N_per_bin = 30;
% %     xsort = sort(X3);
% %     Nsort = ceil(length(xsort)/N_per_bin);
% %     xleft = zeros(1,Nsort);
% %     xright = zeros(1,Nsort);
% %     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
% %     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
%     
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X3<=xright(i) & X3>xleft(i));   % LLCs
%         if(prop==20)
%             is4=sort(abs(Y3(b4)));
%         else
%             is4=sort(Y3(b4));
%         end
%         ycen(i) = median(is4);
%         if(length(is4)>=7)
%             nl(i) = ceil(length(is4)/6);
%             nu(i) = ceil(5*length(is4)/6);            
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(2) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','k',...
%         'marker','none','DisplayName','$LLCs$');
% 
%     if(prop==13 | prop==14)
%         % Solar metallicity
%         x=[-2:0.3:4];
%         y=8.69.*x./x;
%         plot(x,y,'Parent',gca,'marker','none','color','k','linewidth',2,'linestyle','-.')
%     end
%     
%     % Create legend    
%     legend1 = legend(gca,l(2:4));
%     set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
%         'Location','SouthWest','FontSize',14,...
%         'Interpreter','latex');
%     legend boxoff
%     
%     filename1 = strcat(filename,'_SLC.jpg');
%     print(gcf,'-djpeg',filename1);
%     
%     xlabel(xlab2,'Interpreter','latex','FontSize',24,...
%         'FontName','Times New Roman','units','normalized',...
%         'position',[0.5 -0.07 0]);
%     ylabel(ylab2,'Interpreter','latex','FontSize',24,...
%         'FontName','Times New Roman','units','normalized',...
%         'position',[-0.09 0.5 0]);
%     title('');
%     set(legend1,'FontSize',18,'Location','SouthWest');
%     filename1 = strcat(filename,'_SLC.eps');
%     print(gcf,'-depsc',filename1);
%     close all
%     fclose all    
%     
%     %%% Zero Lived
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%     nice_fig(xl0,xu0,yl0,yu0,dx,dy,xlab1,ylab1,tit,16,20,14,1,0,1)    
%     % Create plot
%     l(1)=plot(X4,Y4,'Parent',gca,'MarkerFaceColor','b',...
%         'MarkerEdgeColor','b',...
%         'MarkerSize',3,...
%         'Marker','o',...
%         'LineStyle','none');
%     
%     xleft  = -2:dbin:1;
%     xright = xleft+dbin;
% %     N_per_bin = 30;
% %     xsort = sort(X3);
% %     Nsort = ceil(length(xsort)/N_per_bin);
% %     xleft = zeros(1,Nsort);
% %     xright = zeros(1,Nsort);
% %     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
% %     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
%     
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X4<=xright(i) & X4>xleft(i));   % ZLCs
%         if(prop==20)
%             is4=sort(abs(Y4(b4)));
%         else
%             is4=sort(Y4(b4));
%         end
%         ycen(i) = median(is4);
%         nl(i) = ceil(length(is4)/6);
%         nu(i) = ceil(5*length(is4)/6);
%         if(nl(i)>0 & nu(i)>nl(i))
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(4) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','g',...
%         'marker','none','DisplayName','$ZLCs$');
%     
%     xleft  = -2:dbin:1;
%     xright = xleft+dbin;
% %     N_per_bin = 30;
% %     xsort = sort(X2);
% %     Nsort = ceil(length(xsort)/N_per_bin);
% %     xleft = zeros(1,Nsort);
% %     xright = zeros(1,Nsort);
% %     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
% %     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
%     
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X2<=xright(i) & X2>xleft(i));   % SLCs
%         if(prop==20)
%             is4=sort(abs(Y2(b4)));
%         else
%             is4=sort(Y2(b4));
%         end
%         ycen(i) = median(is4);
%         nl(i) = ceil(length(is4)/6);
%         nu(i) = ceil(5*length(is4)/6);
%         if(nl(i)>0 & nu(i)>nl(i))
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(3) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','r',...
%         'marker','none','DisplayName','$SLCs$');
%         
%     xleft  = -2:dbin:1;
%     xright = xleft+dbin;
% %     N_per_bin = 30;
% %     xsort = sort(X3);
% %     Nsort = ceil(length(xsort)/N_per_bin);
% %     xleft = zeros(1,Nsort);
% %     xright = zeros(1,Nsort);
% %     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
% %     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
%     
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X3<=xright(i) & X3>xleft(i));   % LLCs
%         if(prop==20)
%             is4=sort(abs(Y3(b4)));
%         else
%             is4=sort(Y3(b4));
%         end
%         ycen(i) = median(is4);
%         nl(i) = ceil(length(is4)/6);
%         nu(i) = ceil(5*length(is4)/6);
%         if(nl(i)>0 & nu(i)>nl(i))
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(2) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','k',...
%         'marker','none','DisplayName','$LLCs$');
% 
%     if(prop==13 | prop==14)
%         % Solar metallicity
%         x=[-2:0.3:4];
%         y=8.69.*x./x;
%         plot(x,y,'Parent',gca,'marker','none','color','k','linewidth',2,'linestyle','-.')
%     end
%     
%     % Create legend    
%     legend1 = legend(gca,l(2:4));
%     set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
%         'Location','SouthWest','FontSize',14,...
%         'Interpreter','latex');
%     legend boxoff
%     
%     filename1 = strcat(filename,'_ZLC.jpg');
%     print(gcf,'-djpeg',filename1);
%     
%     xlabel(xlab2,'Interpreter','latex','FontSize',24,...
%         'FontName','Times New Roman','units','normalized',...
%         'position',[0.5 -0.07 0]);
%     ylabel(ylab2,'Interpreter','latex','FontSize',24,...
%         'FontName','Times New Roman','units','normalized',...
%         'position',[-0.09 0.5 0]);
%     title('');
%     set(legend1,'FontSize',18,'Location','SouthWest');
%     filename1 = strcat(filename,'_ZLC.eps');
%     print(gcf,'-depsc',filename1);
%     close all
%     fclose all
%     
%         %%% all
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%     nice_fig(xl0,xu0,yl0,yu0,dx,dy,xlab1,ylab1,tit,16,20,14,1,0,1)    
%     % Create plot
%     plot(X4,Y4,'Parent',gca,'MarkerFaceColor','b',...
%         'MarkerEdgeColor','b',...
%         'MarkerSize',3,...
%         'Marker','o',...
%         'LineStyle','none'); 
%     plot(X2,Y2,'Parent',gca,'MarkerFaceColor',[0.8706 0.4902 0],...
%         'MarkerEdgeColor',[0.8706 0.4902 0],...
%         'MarkerSize',3,...
%         'Marker','o',...
%         'LineStyle','none');  
%     plot(X3,Y3,'Parent',gca,'MarkerFaceColor',[0.5 0.5 0.5],...
%         'MarkerEdgeColor',[0.5 0.5 0.5],...
%         'MarkerSize',3,...
%         'Marker','o',...
%         'LineStyle','none');
%     
%     xleft  = -2:dbin:1;
%     xright = xleft+dbin;
% %     N_per_bin = 30;
% %     xsort = sort(X3);
% %     Nsort = ceil(length(xsort)/N_per_bin);
% %     xleft = zeros(1,Nsort);
% %     xright = zeros(1,Nsort);
% %     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
% %     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
%     
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X4<=xright(i) & X4>xleft(i));   % ZLCs
%         if(prop==20)
%             is4=sort(abs(Y4(b4)));
%         else
%             is4=sort(Y4(b4));
%         end
%         ycen(i) = median(is4);
%         nl(i) = ceil(length(is4)/6);
%         nu(i) = ceil(5*length(is4)/6);
%         if(nl(i)>0 & nu(i)>nl(i))
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(4) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','g',...
%         'marker','none','DisplayName','$ZLCs$');
%     
%     xleft  = -2:dbin:1;
%     xright = xleft+dbin;
% %     N_per_bin = 30;
% %     xsort = sort(X2);
% %     Nsort = ceil(length(xsort)/N_per_bin);
% %     xleft = zeros(1,Nsort);
% %     xright = zeros(1,Nsort);
% %     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
% %     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
%     
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X2<=xright(i) & X2>xleft(i));   % SLCs
%         if(prop==20)
%             is4=sort(abs(Y2(b4)));
%         else
%             is4=sort(Y2(b4));
%         end
%         ycen(i) = median(is4);
%         nl(i) = ceil(length(is4)/6);
%         nu(i) = ceil(5*length(is4)/6);
%         if(nl(i)>0 & nu(i)>nl(i))
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(2) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','r',...
%         'marker','none','DisplayName','$SLCs$');
%         
%     xleft  = -2:dbin:1;
%     xright = xleft+dbin;
% %     N_per_bin = 30;
% %     xsort = sort(X3);
% %     Nsort = ceil(length(xsort)/N_per_bin);
% %     xleft = zeros(1,Nsort);
% %     xright = zeros(1,Nsort);
% %     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
% %     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
%     
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X3<=xright(i) & X3>xleft(i));   % LLCs
%         if(prop==20)
%             is4=sort(abs(Y3(b4)));
%         else
%             is4=sort(Y3(b4));
%         end
%         ycen(i) = median(is4);
%         nl(i) = ceil(length(is4)/6);
%         nu(i) = ceil(5*length(is4)/6);
%         if(nl(i)>0 & nu(i)>nl(i))
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(3) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','k',...
%         'marker','none','DisplayName','$LLCs$');
%     
%     if(prop==13 | prop==14)
%         % Solar metallicity
%         x=[-2:0.3:4];
%         y=8.69.*x./x;
%         plot(x,y,'Parent',gca,'marker','none','color','k','linewidth',2,'linestyle','-.')
%     end
%     
%     % Create legend    
%     legend1 = legend(gca,l(2:4));
%     set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
%         'Location','SouthWest','FontSize',14,...
%         'Interpreter','latex');
%     legend boxoff
%     
%     filename1 = strcat(filename,'_all.jpg');
%     print(gcf,'-djpeg',filename1);
%     
%     xlabel(xlab2,'Interpreter','latex','FontSize',24,...
%         'FontName','Times New Roman','units','normalized',...
%         'position',[0.5 -0.06 0]);
%     ylabel(ylab2,'Interpreter','latex','FontSize',24,...
%         'FontName','Times New Roman','units','normalized',...
%         'position',[-0.09 0.5 0]);
%     title('');
%     set(legend1,'FontSize',18,'Location','SouthWest');
%     filename1 = strcat(filename,'_all.eps');
%     print(gcf,'-depsc',filename1);
%     close all
%     fclose all
    
    %%% Two Medians
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%     nice_fig(xl0,xu0,yl0,yu0,dx,dy,xlab1,ylab1,tit,16,20,14,1,0,1)    
%     
%     xleft  = -2:dbin:1;
%     xright = xleft+dbin;
% %     N_per_bin = 30;
% %     xsort = sort(X2);
% %     Nsort = ceil(length(xsort)/N_per_bin);
% %     xleft = zeros(1,Nsort);
% %     xright = zeros(1,Nsort);
% %     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
% %     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
%     
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X2<=xright(i) & X2>xleft(i));   % SLCs
%         if(prop==20)
%             is4=sort(abs(Y2(b4)));
%         else
%             is4=sort(Y2(b4));
%         end
%         ycen(i) = median(is4);
%         if(length(is4)>=7)
%             nl(i) = ceil(length(is4)/6);
%             nu(i) = ceil(5*length(is4)/6);            
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(1) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','r',...
%         'marker','o','markerfacecolor','r','markeredgecolor','r','markersize',6,'DisplayName','$SLCs$');
%     plot(xcen(b4),yu(b4),'linestyle','--','linewidth',2,'color','r','marker','o',...
%         'markerfacecolor','r','markeredgecolor','r','markersize',4);
%     plot(xcen(b4),yl(b4),'linestyle','--','linewidth',2,'color','r','marker','o',...
%         'markerfacecolor','r','markeredgecolor','r','markersize',4);
%         
%     xleft  = -2:dbin:1;
%     xright = xleft+dbin;
% %     N_per_bin = 30;
% %     xsort = sort(X3);
% %     Nsort = ceil(length(xsort)/N_per_bin);
% %     xleft = zeros(1,Nsort);
% %     xright = zeros(1,Nsort);
% %     xleft(1:Nsort) = xsort(1+N_per_bin.*(0:Nsort-1));
% %     xright(1:Nsort) = [xsort(N_per_bin.*(1:Nsort-1))',xsort(length(xsort))'];
%     
%     xcen   = (xleft+xright)/2;
%     ycen   = -1e9.*ones(size(xcen));
%     nl   = ycen;
%     nu   = ycen;
%     yl   = ycen;
%     yu   = ycen;
%     for i=1:length(xcen)
%         b4=find(X3<=xright(i) & X3>xleft(i));   % LLCs
%         if(prop==20)
%             is4=sort(abs(Y3(b4)));
%         else
%             is4=sort(Y3(b4));
%         end
%         ycen(i) = median(is4);
%         if(length(is4)>=7)
%             nl(i) = ceil(length(is4)/6);
%             nu(i) = ceil(5*length(is4)/6);            
%             yu(i)=is4(nu(i));
%             yl(i)=is4(nl(i));
%         end
%     end
%     b4 = find(yu ~= -1e9);
%     l(2) = plot(xcen(b4),ycen(b4),'linestyle','-','linewidth',3,'color','k',...
%         'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',6,'DisplayName','$LLCs$');
%     plot(xcen(b4),yu(b4),'linestyle','--','linewidth',2,'color','k','marker','o',...
%         'markerfacecolor','k','markeredgecolor','k','markersize',4);
%     plot(xcen(b4),yl(b4),'linestyle','--','linewidth',2,'color','k','marker','o',...
%         'markerfacecolor','k','markeredgecolor','k','markersize',4);
%     
%     if(prop==13 | prop==14)
%         % Solar metallicity
%         x=[-2:0.3:4];
%         y=8.69.*x./x;
%         plot(x,y,'Parent',gca,'marker','none','color','k','linewidth',2,'linestyle','-.')
%     end
%     
%     % Create legend    
%     legend1 = legend(gca,l(1:2));
%     set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
%         'Location','SouthWest','FontSize',14,...
%         'Interpreter','latex');
%     legend boxoff
%     
%     filename1 = strcat(filename,'_two_medians.jpg');
%     print(gcf,'-djpeg',filename1);
%     
%     xlabel(xlab2,'Interpreter','latex','FontSize',24,...
%     'FontName','Times New Roman','units','normalized',...
%     'position',[0.5 -0.06 0]);
%     ylabel(ylab2,'Interpreter','latex','FontSize',24,...
%         'FontName','Times New Roman','units','normalized',...
%         'position',[-0.08 0.5 0]);
%     title('');
%     set(legend1,'FontSize',18,'units','normalized','position',leg_pos2);
%     filename1 = strcat(filename,'_two_medians.eps');
%     print(gcf,'-depsc',filename1);
%     close all
%     fclose all
end