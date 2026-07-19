function launch_classification_v4_parameters(insitu, exsitu, bulge_clump, norm_insitu, ...
    norm_exsitu, norm_bulge_clump, nis, nes, nbulge, gal, Delt, Fw, type)

% xtype = 1 --> Mc, Mc/Md        (baryons, baryons disc)
% xtype = 2 --> Ms, Ms/Msd       (stars,   stars disc)
% xtype = 3 --> Mg, Mg/Mgd       (gas,     gas disc)
% xtype = 4 --> Mc, Mc/Mcold,d   (baryons, cold disc)
% xtype = 5 --> Ms, Ms/Ms,sphere (stars,   stars sphere)

if(type==1)
    type_tit = 'number';
elseif(type==2)
    type_tit = 'redshift';
elseif(type==3)
    type_tit = 'DM_residual';
elseif(type==4)
    type_tit = 'residual';
elseif(type==5)
    type_tit = 'shape';
elseif(type==6)
    type_tit = 'gas_fraction';
elseif(type==7)
    type_tit = 'distance';
elseif(type==8)
    type_tit = 'height';
elseif(type==9)
    type_tit = 'normalized_distance';
elseif(type==10)
    type_tit = 'normalized_height';
elseif(type==11)
    type_tit = 'SFR';
elseif(type==12)
    type_tit = 'Mstar';
elseif(type==14)
    type_tit = 'Mgas';
elseif(type==15)
    type_tit = 'stellar_age';
elseif(type==16)
    type_tit = 'sSFR';
elseif(type==17)
    type_tit = 'gas_metallicity';
elseif(type==18)
    type_tit = 'stellar_metallicity';
elseif(type==19)
    type_tit = 'gas_surface_density';
elseif(type==20)
    type_tit = 'stellar_surface_density';
elseif(type==21)
    type_tit = 'depletion_time';
elseif(type==22)
    type_tit = 'local_dynamical_over_free_fall';
elseif(type==23)
    type_tit = 'global_dynamical_over_free_fall';
elseif(type==24)
    type_tit = 'time_since_formation';
elseif(type==25)
    type_tit = 'normalized_SFR';
elseif(type==26)
    type_tit = 'max_normalized_dist';
elseif(type==27)
    type_tit = 'max_normalized_height';
elseif(type==28)
    type_tit = 'Sigma_SFR';
elseif(type==29)
    type_tit = 'normalized_Sigma_SFR';
elseif(type==30)
    type_tit = 'normalized_stellar_surface_density';
elseif(type==31)
    type_tit = 't_over_tff';
elseif(type==32)
    type_tit = 'normalized_gas_fraction';
elseif(type==33)
    type_tit = 'normalized_ssfr';
elseif(type==34)
    type_tit = 'normalized_gas_metalicity';
elseif(type==35)
    type_tit = 'normalized_stellar_age';
elseif(type==36)
    type_tit = 'normalized_SFR_in_sphere';
elseif(type==37)
    type_tit = 'normalized_gas_fraction_sphere';
elseif(type==38)
    type_tit = 'normalized_ssfr_sphere';
elseif(type==39)
    type_tit = 'normalized_gas_metalicity_sphere';
elseif(type==40)
    type_tit = 'normalized_stellar_age_sphere';
end

mkdir('./mass_size_plane/');
gal_pref = strcat('VELA_v2_',num2str(gal,'%02i'));
is = [];
es = [];
bulge = [];
norm_is = [];
norm_es = [];
norm_bulge = [];
Nis = 0;
Nes = 0;
Nbulge = 0;

for i=1:length(Delt)
    if(nis(i+1)>nis(i))
        is = insitu(nis(i)+1:nis(i+1),:);
        norm_is = norm_insitu(nis(i)+1:nis(i+1),:);
        Nis = length(is(:,2));
    end
    if(nes(i+1)>nes(i))
        es = exsitu(nes(i)+1:nes(i+1),:);
        norm_es = norm_exsitu(nes(i)+1:nes(i+1),:);
        Nes = length(es(:,2));
    end
    if(nbulge(i+1)>nbulge(i))
        bulge = bulge_clump(nbulge(i)+1:nbulge(i+1),:);
        norm_bulge = norm_bulge_clump(nbulge(i)+1:nbulge(i+1),:);
        Nbulge = length(bulge(:,2));
    end
    if(Nis>0)
        mkdir(strcat('./mass_size_plane/param_check/',gal_pref));
        mkdir(strcat('./mass_size_plane/param_check/',gal_pref,'/',type_tit));
        mkdir(strcat('./mass_size_plane/param_check/',gal_pref,'/',type_tit,'/IS_only'));
        tit1 = strcat('$F_{\rm W}=',num2str(Fw(i),'%3.1f'),'kpc\:\:\delta_{\rm min}=',num2str(Delt(i),'%02i'),'\:\:N_{\rm c}=',num2str(Nis,'%5i'),'$');
        
        classification_v4_parameters(bulge, is, es, norm_is, norm_es, norm_bulge, 0, 1, tit1, type);
        filename = strcat('./mass_size_plane/param_check/',gal_pref,'/',type_tit,'/IS_only/Fw_',num2str(10*Fw(i),'%02i'),'_Dt_',num2str(Delt(i),'%02i'),'.jpg');
        print(gcf,'-djpeg',filename);
        close all
        fclose all
        classification_v4_parameters(bulge, is, es, norm_is, norm_es, norm_bulge, 1, 1, tit1, type);
        filename = strcat('./mass_size_plane/param_check/',gal_pref,'/',type_tit,'/IS_only/normalized_Fw_',num2str(10*Fw(i),'%02i'),'_Dt_',num2str(Delt(i),'%02i'),'.jpg');
        print(gcf,'-djpeg',filename);
        close all
        fclose all
     end
    if(Nis+Nes+Nbulge>0)
        mkdir(strcat('./mass_size_plane/param_check/',gal_pref));
        mkdir(strcat('./mass_size_plane/param_check/',gal_pref,'/',type_tit));
        mkdir(strcat('./mass_size_plane/param_check/',gal_pref,'/',type_tit,'/All'));
        tit1 = strcat('$F_{\rm W}=',num2str(Fw(i),'%3.1f'),'kpc\:\:\delta_{\rm min}=',num2str(Delt(i),'%02i'),'\:\:N_{\rm c}=',num2str(Nis,'%5i'));
        
        classification_v4_parameters(bulge, is, es, norm_is, norm_es, norm_bulge, 0, 0, tit1, type);
        filename = strcat('./mass_size_plane/param_check/',gal_pref,'/',type_tit,'/All/Fw_',num2str(10*Fw(i),'%02i'),'_Dt_',num2str(Delt(i),'%02i'),'.jpg');
        print(gcf,'-djpeg',filename);
        close all
        fclose all
        classification_v4_parameters(bulge, is, es, norm_is, norm_es, norm_bulge, 1, 0, tit1, type);
        filename = strcat('./mass_size_plane/param_check/',gal_pref,'/',type_tit,'/All/normalized_Fw_',num2str(10*Fw(i),'%02i'),'_Dt_',num2str(Delt(i),'%02i'),'.jpg');
        print(gcf,'-djpeg',filename);
        close all
        fclose all
    end
end