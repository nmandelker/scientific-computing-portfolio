ind2 = [1:3, 5:16, 19, 21, 23, 25:35];
ind3 = [1:17, 19:35];
ind_common = ind2;
%ind_common = [1:3, 5:16, 19, 21, 23, 25:34];
N_common = length(ind_common);
Amin_common = zeros(N_common,1);
Amax_common = zeros(N_common,1);

disc_cat2     = [];
Mstar_cat2    = [];
sphere_cat2   = [];
Rvir_015_cat2 = [];
N_cat2        = [0];
disc_cat3     = [];
Mstar_cat3    = [];
sphere_cat3   = [];
Rvir_015_cat3 = [];
N_cat3        = [0];

cat2_common   = [];
N_cat2_common = zeros(1,length(ind2));
cat3_common   = [];
N_cat3_common = zeros(1,length(ind3));

for i=1:length(ind2)
    A1 = importdata(strcat('./galaxy_catalogues/','VELA',num2str(ind2(i),'%02i'),'/galaxy_catalogue/Nir_disc_cat.txt'));
    N_cat2 = [N_cat2, N_cat2(i) + A1(1)];
    disc_cat2 = [disc_cat2, reshape(A1(2:length(A1)),[(length(A1)-1)/A1(1), A1(1)])];
    
    B1 = importdata(strcat('./galaxy_catalogues/','VELA',num2str(ind2(i),'%02i'),'/galaxy_catalogue/Mstar.txt'));
    Mstar_cat2 = [Mstar_cat2, reshape(B1(2:length(B1)),[(length(B1)-1)/B1(1), B1(1)])];
    
    D1 = importdata(strcat('./galaxy_catalogues/','VELA',num2str(ind2(i),'%02i'),'/galaxy_catalogue/Nir_015_Rvir_cat.txt'));
    Rvir_015_cat2 = [Rvir_015_cat2, reshape(D1(2:length(D1)),[(length(D1)-1)/D1(1), D1(1)])];
    
    C1 = importdata(strcat('./galaxy_catalogues/','VELA',num2str(ind2(i),'%02i'),'/galaxy_catalogue/Nir_spherical_galaxy_cat.txt'));
    sphere_cat2 = [sphere_cat2, reshape(C1(2:length(C1)),[(length(C1)-1)/C1(1), C1(1)])];
    clear A1 B1 C1 D1
end
for i=1:length(ind3)    
    A1 = importdata(strcat('./galaxy_catalogues/','VELA_v2_',num2str(ind3(i),'%02i'),'/galaxy_catalogue/Nir_disc_cat.txt'));
    N_cat3 = [N_cat3, N_cat3(i) + A1(1)];
    disc_cat3 = [disc_cat3, reshape(A1(2:length(A1)),[(length(A1)-1)/A1(1), A1(1)])];
    
    B1 = importdata(strcat('./galaxy_catalogues/','VELA_v2_',num2str(ind3(i),'%02i'),'/galaxy_catalogue/Mstar.txt'));
    Mstar_cat3 = [Mstar_cat3, reshape(B1(2:length(B1)),[(length(B1)-1)/B1(1), B1(1)])];
    
    D1 = importdata(strcat('./galaxy_catalogues/','VELA_v2_',num2str(ind3(i),'%02i'),'/galaxy_catalogue/Nir_015_Rvir_cat.txt'));
    Rvir_015_cat3 = [Rvir_015_cat3, reshape(D1(2:length(D1)),[(length(D1)-1)/D1(1), D1(1)])];
    
    C1 = importdata(strcat('./galaxy_catalogues/','VELA_v2_',num2str(ind3(i),'%02i'),'/galaxy_catalogue/Nir_spherical_galaxy_cat.txt'));
    sphere_cat3 = [sphere_cat3, reshape(C1(2:length(C1)),[(length(C1)-1)/C1(1), C1(1)])];
    clear A1 B1 C1 D1
end
disc_cat2     = disc_cat2';
Mstar_cat2    = Mstar_cat2';
sphere_cat2   = sphere_cat2';
Rvir_015_cat2 = Rvir_015_cat2';
disc_cat3     = disc_cat3';
Mstar_cat3    = Mstar_cat3';
sphere_cat3   = sphere_cat3';
Rvir_015_cat3 = Rvir_015_cat3';

for i=1:length(ind_common)
    j = ind_common(i);    
    b2 = find(ind2==j);
    Amin2 = disc_cat2(N_cat2(b2)+1,1);
    Amax2 = disc_cat2(N_cat2(b2+1),1);
    b3 = find(ind3==j);
    Amin3 = disc_cat3(N_cat3(b3)+1,1);
    Amax3 = disc_cat3(N_cat3(b3+1),1);
    Amin_common(i) = max(Amin2,Amin3);
    Amax_common(i) = min(Amax2,Amax3);
    
    % Now fix Amin/Amax for bugs
    if(j==1)    % VELA 01
        Amin_common(i) = max(0.42,Amin_common(i));   %Before this, gen2 and gen3 trace different progenitors of a merger
    end
    if(j==3)    % VELA 03
        Amin_common(i) = max(0.18,Amin_common(i));   %I don't understand why, but my clump catalogues are missin a=0.17 for gen3
    end
    if(j==9)    % VELA 09
        Amax_common(i) = min(0.40,Amax_common(i));   %Bug after this. Rvir decreases and gas disc dissapears
    end
    if(j==11)    % VELA 11
        Amax_common(i) = min(0.46,Amax_common(i));   %Bug after this
    end
    if(j==12)    % VELA 12
        Amax_common(i) = min(0.44,Amax_common(i));   %Bug after this
    end
    if(j==14)    % VELA 14
        Amax_common(i) = min(0.41,Amax_common(i));   %Bug after this
    end
    if(j==34)    % VELA 34
        Amin_common(i) = max(0.20,Amin_common(i));   %Before this, gen2 and gen3 trace different progenitors of a merger
    end
    if(j==35)    % VELA 34
        Amin_common(i) = max(0.12,Amin_common(i));   %I don't understand why, but my clump catalogues start from a=0.12
    end
    
    a3  = disc_cat3(N_cat3(b3)+1:N_cat3(b3+1),1);
    a2  = disc_cat2(N_cat2(b2)+1:N_cat2(b2+1),1);
    for k=(N_cat2(b2)+1):N_cat2(b2+1)
        adisc2 = disc_cat2(k,1);
        if(adisc2<(Amax_common(i)+0.005) & adisc2>(Amin_common(i)-0.005))
            if(min(abs(a3(:)-adisc2))<0.005)
                N_cat2_common(b2) = N_cat2_common(b2) + 1;
                cat2_common = [cat2_common', k]';
                N_cat3_common(b3) =  N_cat3_common(b3) + 1;
                bcommon = find(abs(a3(:)-adisc2)<0.005);
                cat3_common = [cat3_common', N_cat3(b3)+bcommon']';                
            end
        end
    end
end

clear b2 b3 i j k adisc2 a2 a3 Amax2 Amax3 Amin2 Amin3 bcommon N_common Amax_common Amin_common