prop_list(1) = {'id'};
prop_list(2) = {'z'};
prop_list(3) = {'Rc'};
prop_list(4) = {'Mgas'};
prop_list(5) = {'Mstar'};
prop_list(6) = {'Mbar'};
prop_list(7) = {'fgas'};
prop_list(8) = {'fdm'};
prop_list(9) = {'SigGas'};
prop_list(10) = {'SigStar'};
prop_list(11) = {'SigBar'};
prop_list(12) = {'AgeStar'};
prop_list(13) = {'ZGas'};
prop_list(14) = {'ZStar'};
prop_list(15) = {'SFR'};
prop_list(16) = {'SigSFR'};
prop_list(17) = {'sSFR'};
prop_list(18) = {'tdep'};
prop_list(19) = {'norm_dist'};
prop_list(20) = {'norm_height'};
prop_list(21) = {'dist'};
prop_list(22) = {'height'};
prop_list(23) = {'residual'};
prop_list(24) = {'shape1'};
prop_list(25) = {'shape2'};
prop_list(26) = {'shape3'};
prop_list(27) = {'DM_residual'};
prop_list(28) = {'Es'};
prop_list(29) = {'tff'};
prop_list(30) = {'td_local'};
prop_list(31) = {'td_global'};
prop_list(32) = {'Mgas_in_1_0_Rc'};
prop_list(33) = {'Mgas_in_1_5_Rc'};
prop_list(34) = {'Mgas_in_2_0_Rc'};
prop_list(35) = {'Mgas_out_1_0_Rc'};
prop_list(36) = {'Mgas_out_1_5_Rc'};
prop_list(37) = {'Mgas_out_2_0_Rc'};
prop_list(38) = {'Mstar_in'};
prop_list(39) = {'Mstar_out'};
prop_list(40) = {'Mstar_formed'};
prop_list(41) = {'time'};

unit_list(1) = {'1'};
unit_list(2) = {'1'};
unit_list(3) = {'kpc'};
unit_list(4) = {'Msun'};
unit_list(5) = {'Msun'};
unit_list(6) = {'Msun'};
unit_list(7) = {'1'};
unit_list(8) = {'1'};
unit_list(9) = {'Msun/pc/pc'};
unit_list(10) = {'Msun/pc/pc'};
unit_list(11) = {'Msun/pc/pc'};
unit_list(12) = {'Myr'};
unit_list(13) = {'Log O/H + 12'};
unit_list(14) = {'Log O/H + 12'};
unit_list(15) = {'Msun/yr'};
unit_list(16) = {'Msun/yr/kpc/kpc'};
unit_list(17) = {'1/Gyr'};
unit_list(18) = {'Gyr'};
unit_list(19) = {'1'};
unit_list(20) = {'1'};
unit_list(21) = {'kpc'};
unit_list(22) = {'kpc'};
unit_list(23) = {'1'};
unit_list(24) = {'1'};
unit_list(25) = {'1'};
unit_list(26) = {'1'};
unit_list(27) = {'1'};
unit_list(28) = {'1'};
unit_list(29) = {'Myr'};
unit_list(30) = {'Myr'};
unit_list(31) = {'Myr'};
unit_list(32) = {'Msun/yr'};
unit_list(33) = {'Msun/yr'};
unit_list(34) = {'Msun/yr'};
unit_list(35) = {'Msun/yr'};
unit_list(36) = {'Msun/yr'};
unit_list(37) = {'Msun/yr'};
unit_list(38) = {'Msun/yr'};
unit_list(39) = {'Msun/yr'};
unit_list(40) = {'Msun/yr'};
unit_list(41) = {'Myr'};

% local td over tff
is2(:,42)         = is2(:,30)./is2(:,29);
es2(:,42)         = es2(:,30)./es2(:,29);
bulge2(:,42)      = bulge2(:,30)./bulge2(:,29);
norm_is2(:,42)    = norm_is2(:,30)./norm_is2(:,29);
norm_es2(:,42)    = norm_es2(:,30)./norm_es2(:,29);
norm_bulge2(:,42) = norm_bulge2(:,30)./norm_bulge2(:,29);
is3(:,42)         = is3(:,30)./is3(:,29);
es3(:,42)         = es3(:,30)./es3(:,29);
bulge3(:,42)      = bulge3(:,30)./bulge3(:,29);
norm_is3(:,42)    = norm_is3(:,30)./norm_is3(:,29);
norm_es3(:,42)    = norm_es3(:,30)./norm_es3(:,29);
norm_bulge3(:,42) = norm_bulge3(:,30)./norm_bulge3(:,29);
prop_list(42)     = {'td local over tff'};
unit_list(42)     = {'1'};

% global td over tff
is2(:,43)         = is2(:,31)./is2(:,29);
es2(:,43)         = es2(:,31)./es2(:,29);
bulge2(:,43)      = bulge2(:,31)./bulge2(:,29);
norm_is2(:,43)    = norm_is2(:,31)./norm_is2(:,29);
norm_es2(:,43)    = norm_es2(:,31)./norm_es2(:,29);
norm_bulge2(:,43) = norm_bulge2(:,31)./norm_bulge2(:,29);
is3(:,43)         = is3(:,31)./is3(:,29);
es3(:,43)         = es3(:,31)./es3(:,29);
bulge3(:,43)      = bulge3(:,31)./bulge3(:,29);
norm_is3(:,43)    = norm_is3(:,31)./norm_is3(:,29);
norm_es3(:,43)    = norm_es3(:,31)./norm_es3(:,29);
norm_bulge3(:,43) = norm_bulge3(:,31)./norm_bulge3(:,29);
prop_list(43)     = {'td global over tff'};
unit_list(43)     = {'1'};

% average Mgas in
is2(:,44)         = (is2(:,32)        +is2(:,33)        +is2(:,34))         ./ 3;
es2(:,44)         = (es2(:,32)        +es2(:,33)        +es2(:,34))         ./ 3;
bulge2(:,44)      = (bulge2(:,32)     +bulge2(:,33)     +bulge2(:,34))      ./ 3;
norm_is2(:,44)    = (norm_is2(:,32)   +norm_is2(:,33)   +norm_is2(:,34))    ./ 3;
norm_es2(:,44)    = (norm_es2(:,32)   +norm_es2(:,33)   +norm_es2(:,34))    ./ 3;
norm_bulge2(:,44) = (norm_bulge2(:,32)+norm_bulge2(:,33)+norm_bulge2(:,34)) ./ 3;
is3(:,44)         = (is3(:,32)        +is3(:,33)        +is3(:,34))         ./ 3;
es3(:,44)         = (es3(:,32)        +es3(:,33)        +es3(:,34))         ./ 3;
bulge3(:,44)      = (bulge3(:,32)     +bulge3(:,33)     +bulge3(:,34))      ./ 3;
norm_is3(:,44)    = (norm_is3(:,32)   +norm_is3(:,33)   +norm_is3(:,34))    ./ 3;
norm_es3(:,44)    = (norm_es3(:,32)   +norm_es3(:,33)   +norm_es3(:,34))    ./ 3;
norm_bulge3(:,44) = (norm_bulge3(:,32)+norm_bulge3(:,33)+norm_bulge3(:,34)) ./ 3;
prop_list(44)     = {'average Mgas in'};
unit_list(44)     = {'Msun/yr'};

% average Mgas out
is2(:,45)         = (is2(:,35)        +is2(:,36)        +is2(:,37))         ./ 3;
es2(:,45)         = (es2(:,35)        +es2(:,36)        +es2(:,37))         ./ 3;
bulge2(:,45)      = (bulge2(:,35)     +bulge2(:,36)     +bulge2(:,37))      ./ 3;
norm_is2(:,45)    = (norm_is2(:,35)   +norm_is2(:,36)   +norm_is2(:,37))    ./ 3;
norm_es2(:,45)    = (norm_es2(:,35)   +norm_es2(:,36)   +norm_es2(:,37))    ./ 3;
norm_bulge2(:,45) = (norm_bulge2(:,35)+norm_bulge2(:,36)+norm_bulge2(:,37)) ./ 3;
is3(:,45)         = (is3(:,35)        +is3(:,36)        +is3(:,37))         ./ 3;
es3(:,45)         = (es3(:,35)        +es3(:,36)        +es3(:,37))         ./ 3;
bulge3(:,45)      = (bulge3(:,35)     +bulge3(:,36)     +bulge3(:,37))      ./ 3;
norm_is3(:,45)    = (norm_is3(:,35)   +norm_is3(:,36)   +norm_is3(:,37))    ./ 3;
norm_es3(:,45)    = (norm_es3(:,35)   +norm_es3(:,36)   +norm_es3(:,37))    ./ 3;
norm_bulge3(:,45) = (norm_bulge3(:,35)+norm_bulge3(:,36)+norm_bulge3(:,37)) ./ 3;
prop_list(45)     = {'average Mgas out'};
unit_list(45)     = {'Msun/yr'};

% spherical gas density
cm3 = 1./0.03363;
is2(:,46)         = cm3.*is2(:,4)    ./ ((4*pi/3).*((1000.*is2(:,3)).^3));
es2(:,46)         = cm3.*es2(:,4)    ./ ((4*pi/3).*((1000.*es2(:,3)).^3));
bulge2(:,46)      = cm3.*bulge2(:,4) ./ ((4*pi/3).*((1000.*bulge2(:,3)).^3));
norm_is2(:,46)    = norm_is2(:,4)    ./ (norm_is2(:,3).^3);
norm_es2(:,46)    = norm_es2(:,4)    ./ (norm_es2(:,3).^3);
norm_bulge2(:,46) = norm_bulge2(:,4) ./ (norm_bulge2(:,3).^3);
is3(:,46)         = cm3.*is3(:,4)    ./ ((4*pi/3).*((1000.*is3(:,3)).^3));
es3(:,46)         = cm3.*es3(:,4)    ./ ((4*pi/3).*((1000.*es3(:,3)).^3));
bulge3(:,46)      = cm3.*bulge3(:,4) ./ ((4*pi/3).*((1000.*bulge3(:,3)).^3));
norm_is3(:,46)    = norm_is3(:,4)    ./ (norm_is3(:,3).^3);
norm_es3(:,46)    = norm_es3(:,4)    ./ (norm_es3(:,3).^3);
norm_bulge3(:,46) = norm_bulge3(:,4) ./ (norm_bulge3(:,3).^3);
prop_list(46)     = {'Gas density'};
unit_list(46)     = {'1/cm3'};

% spherical star density
cm3 = 1./0.03363;
is2(:,47)         = cm3.*is2(:,5)    ./ ((4*pi/3).*((1000.*is2(:,3)).^3));
es2(:,47)         = cm3.*es2(:,5)    ./ ((4*pi/3).*((1000.*es2(:,3)).^3));
bulge2(:,47)      = cm3.*bulge2(:,5) ./ ((4*pi/3).*((1000.*bulge2(:,3)).^3));
norm_is2(:,47)    = norm_is2(:,5)    ./ (norm_is2(:,3).^3);
norm_es2(:,47)    = norm_es2(:,5)    ./ (norm_es2(:,3).^3);
norm_bulge2(:,47) = norm_bulge2(:,5) ./ (norm_bulge2(:,3).^3);
is3(:,47)         = cm3.*is3(:,5)    ./ ((4*pi/3).*((1000.*is3(:,3)).^3));
es3(:,47)         = cm3.*es3(:,5)    ./ ((4*pi/3).*((1000.*es3(:,3)).^3));
bulge3(:,47)      = cm3.*bulge3(:,5) ./ ((4*pi/3).*((1000.*bulge3(:,3)).^3));
norm_is3(:,47)    = norm_is3(:,5)    ./ (norm_is3(:,3).^3);
norm_es3(:,47)    = norm_es3(:,5)    ./ (norm_es3(:,3).^3);
norm_bulge3(:,47) = norm_bulge3(:,5) ./ (norm_bulge3(:,3).^3);
prop_list(47)     = {'Star density'};
unit_list(47)     = {'1/cm3'};

% spherical baryonic density
cm3 = 1./0.03363;
is2(:,48)         = cm3.*is2(:,6)    ./ ((4*pi/3).*((1000.*is2(:,3)).^3));
es2(:,48)         = cm3.*es2(:,6)    ./ ((4*pi/3).*((1000.*es2(:,3)).^3));
bulge2(:,48)      = cm3.*bulge2(:,6) ./ ((4*pi/3).*((1000.*bulge2(:,3)).^3));
norm_is2(:,48)    = norm_is2(:,6)    ./ (norm_is2(:,3).^3);
norm_es2(:,48)    = norm_es2(:,6)    ./ (norm_es2(:,3).^3);
norm_bulge2(:,48) = norm_bulge2(:,6) ./ (norm_bulge2(:,3).^3);
is3(:,48)         = cm3.*is3(:,6)    ./ ((4*pi/3).*((1000.*is3(:,3)).^3));
es3(:,48)         = cm3.*es3(:,6)    ./ ((4*pi/3).*((1000.*es3(:,3)).^3));
bulge3(:,48)      = cm3.*bulge3(:,6) ./ ((4*pi/3).*((1000.*bulge3(:,3)).^3));
norm_is3(:,48)    = norm_is3(:,6)    ./ (norm_is3(:,3).^3);
norm_es3(:,48)    = norm_es3(:,6)    ./ (norm_es3(:,3).^3);
norm_bulge3(:,48) = norm_bulge3(:,6) ./ (norm_bulge3(:,3).^3);
prop_list(48)     = {'Bar density'};
unit_list(48)     = {'1/cm3'};

% maximal mass, normalized distance, normalized height and lifetime
is2(:,49:78)=0;
es2(:,49:78)=0;
bulge2(:,49:78)=0;
is3(:,49:78)=0;
es3(:,49:78)=0;
bulge3(:,49:78)=0;
norm_is2(:,49:78)=0;
norm_es2(:,49:78)=0;
norm_bulge2(:,49:78)=0;
norm_is3(:,49:78)=0;
norm_es3(:,49:78)=0;
norm_bulge3(:,49:78)=0;
avg_tff = zeros(1,2);

i=1;
while(i<length(is2))
    if(is2(i,1) == is2(i+1,1) & is2(i,2) > is2(i+1,2))   % multiple snapshot clump
        j = i + 1;
        avg_tff(1) = is2(i,6).*is2(i,29);
        avg_tff(2) = is2(i,6);
        while(is2(i,1)==is2(j,1) & is2(i,2) > is2(j,2) & j<length(is2))
            j = j + 1;
            avg_tff(1) = avg_tff(1) + is2(j,6).*is2(j,29);
            avg_tff(2) = avg_tff(2) + is2(j,6);
        end
        if(j==length(is2) & is2(i,1)==is2(j,1) & is2(i,2) > is2(j,2))
            j = j + 1;
            avg_tff(1) = avg_tff(1) + is2(j,6).*is2(j,29);
            avg_tff(2) = avg_tff(2) + is2(j,6);
        end
        avg_tff(1) = avg_tff(1) ./ avg_tff(2);
        is2(i:(j-1),49) = max(is2(i:(j-1),6));
        is2(i:(j-1),66) = max(is2(i:(j-1),15));
        is2(i:(j-1),50) = max(is2(i:(j-1),19));
        is2(i:(j-1),51) = max(is2(i:(j-1),20));
        is2(i:(j-1),52) = is2(j-1,41);
        norm_is2(i:(j-1),49) = max(norm_is2(i:(j-1),6));
        norm_is2(i:(j-1),66) = max(norm_is2(i:(j-1),15));
        norm_is2(i:(j-1),50) = max(norm_is2(i:(j-1),19));
        norm_is2(i:(j-1),51) = max(norm_is2(i:(j-1),20));
        norm_is2(i:(j-1),52) = is2(j-1,41)./avg_tff(1);
        i = j;
    else
        is2(i,49) = is2(i,6);
        is2(i,66) = is2(i,15);
        is2(i,50) = is2(i,19);
        is2(i,51) = is2(i,20);
        is2(i,52) = is2(i,41);
        norm_is2(i,49) = norm_is2(i,6);
        norm_is2(i,66) = norm_is2(i,15);
        norm_is2(i,50) = norm_is2(i,19);
        norm_is2(i,51) = norm_is2(i,20);
        norm_is2(i,52) = norm_is2(i,41);
        i = i + 1;
    end
end
if(i==length(is2))
    is2(length(is2),49) = is2(length(is2),6);
    is2(length(is2),66) = is2(length(is2),15);
    is2(length(is2),50) = is2(length(is2),19);
    is2(length(is2),51) = is2(length(is2),20);
    is2(length(is2),52) = is2(length(is2),41);
    norm_is2(length(is2),49) = norm_is2(length(is2),6);
    norm_is2(length(is2),66) = norm_is2(length(is2),15);
    norm_is2(length(is2),50) = norm_is2(length(is2),19);
    norm_is2(length(is2),51) = norm_is2(length(is2),20);
    norm_is2(length(is2),52) = norm_is2(length(is2),41);
end
[max(norm_is2(:,41)),min(norm_is2(:,41)),length(find(isnan(norm_is2(:,41))))]
[max(norm_is2(:,52)),min(norm_is2(:,52)),length(find(isnan(norm_is2(:,52))))]
[max(is2(:,41)),min(is2(:,41)),length(find(isnan(is2(:,41))))]
[max(is2(:,52)),min(is2(:,52)),length(find(isnan(is2(:,52))))]

i=1;
while(i<length(es2))
    if(es2(i,1) == es2(i+1,1) & es2(i,2) > es2(i+1,2))   % multiple snapshot clump
        j = i + 1;
        avg_tff(1) = es2(i,6).*es2(i,29);
        avg_tff(2) = es2(i,6);
        while(es2(i,1)==es2(j,1) & es2(i,2) > es2(j,2) & j<length(es2))
            j = j + 1;
            avg_tff(1) = avg_tff(1) + es2(j,6).*es2(j,29);
            avg_tff(2) = avg_tff(2) + es2(j,6);
        end
        if(j==length(es2) & es2(i,1)==es2(j,1) & es2(i,2) > es2(j,2))
            j = j + 1;
            avg_tff(1) = avg_tff(1) + es2(j,6).*es2(j,29);
            avg_tff(2) = avg_tff(2) + es2(j,6);
        end
        avg_tff(1) = avg_tff(1) ./ avg_tff(2);
        es2(i:(j-1),49) = max(es2(i:(j-1),6));
        es2(i:(j-1),66) = max(es2(i:(j-1),15));
        es2(i:(j-1),50) = max(es2(i:(j-1),19));
        es2(i:(j-1),51) = max(es2(i:(j-1),20));
        es2(i:(j-1),52) = es2(j-1,41);
        norm_es2(i:(j-1),49) = max(norm_es2(i:(j-1),6));
        norm_es2(i:(j-1),66) = max(norm_es2(i:(j-1),15));
        norm_es2(i:(j-1),50) = max(norm_es2(i:(j-1),19));
        norm_es2(i:(j-1),51) = max(norm_es2(i:(j-1),20));
        norm_es2(i:(j-1),52) = es2(j-1,41)./avg_tff(1);
        i = j;
    else
        es2(i,49) = es2(i,6);
        es2(i,66) = es2(i,15);
        es2(i,50) = es2(i,19);
        es2(i,51) = es2(i,20);
        es2(i,52) = es2(i,41);
        norm_es2(i,49) = norm_es2(i,6);
        norm_es2(i,66) = norm_es2(i,15);
        norm_es2(i,50) = norm_es2(i,19);
        norm_es2(i,51) = norm_es2(i,20);
        norm_es2(i,52) = norm_es2(i,41);
        i = i + 1;
    end
end
if(i==length(es2))
    es2(length(es2),49) = es2(length(es2),6);
    es2(length(es2),66) = es2(length(es2),15);
    es2(length(es2),50) = es2(length(es2),19);
    es2(length(es2),51) = es2(length(es2),20);
    es2(length(es2),52) = es2(length(es2),41);
    norm_es2(length(es2),49) = norm_es2(length(es2),6);
    norm_es2(length(es2),66) = norm_es2(length(es2),15);
    norm_es2(length(es2),50) = norm_es2(length(es2),19);
    norm_es2(length(es2),51) = norm_es2(length(es2),20);
    norm_es2(length(es2),52) = norm_es2(length(es2),41);
end

i=1;
while(i<length(bulge2))
    if(bulge2(i,1) == bulge2(i+1,1) & bulge2(i,2) > bulge2(i+1,2))   % multiple snapshot clump
        j = i + 1;
        avg_tff(1) = bulge2(i,6).*bulge2(i,29);
        avg_tff(2) = bulge2(i,6);
        while(bulge2(i,1)==bulge2(j,1) & bulge2(i,2) > bulge2(j,2) & j<length(bulge2))
            j = j + 1;
            avg_tff(1) = avg_tff(1) + bulge2(j,6).*bulge2(j,29);
            avg_tff(2) = avg_tff(2) + bulge2(j,6);
        end
        if(j==length(bulge2) & bulge2(i,1)==bulge2(j,1) & bulge2(i,2) > bulge2(j,2))
            j = j + 1;
            avg_tff(1) = avg_tff(1) + bulge2(j,6).*bulge2(j,29);
            avg_tff(2) = avg_tff(2) + bulge2(j,6);
        end
        avg_tff(1) = avg_tff(1) ./ avg_tff(2);
        bulge2(i:(j-1),49) = max(bulge2(i:(j-1),6));
        bulge2(i:(j-1),66) = max(bulge2(i:(j-1),15));
        bulge2(i:(j-1),50) = max(bulge2(i:(j-1),19));
        bulge2(i:(j-1),51) = max(bulge2(i:(j-1),20));
        bulge2(i:(j-1),52) = bulge2(j-1,41);
        norm_bulge2(i:(j-1),49) = max(norm_bulge2(i:(j-1),6));
        norm_bulge2(i:(j-1),66) = max(norm_bulge2(i:(j-1),15));
        norm_bulge2(i:(j-1),50) = max(norm_bulge2(i:(j-1),19));
        norm_bulge2(i:(j-1),51) = max(norm_bulge2(i:(j-1),20));
        norm_bulge2(i:(j-1),52) = bulge2(j-1,41)./avg_tff(1);
        i = j;
    else
        bulge2(i,49) = bulge2(i,6);
        bulge2(i,66) = bulge2(i,15);
        bulge2(i,50) = bulge2(i,19);
        bulge2(i,51) = bulge2(i,20);
        bulge2(i,52) = bulge2(i,41);
        norm_bulge2(i,49) = norm_bulge2(i,6);
        norm_bulge2(i,66) = norm_bulge2(i,15);
        norm_bulge2(i,50) = norm_bulge2(i,19);
        norm_bulge2(i,51) = norm_bulge2(i,20);
        norm_bulge2(i,52) = norm_bulge2(i,41);
        i = i + 1;
    end
end
if(i==length(bulge2))
    bulge2(length(bulge2),49) = bulge2(length(bulge2),6);
    bulge2(length(bulge2),66) = bulge2(length(bulge2),15);
    bulge2(length(bulge2),50) = bulge2(length(bulge2),19);
    bulge2(length(bulge2),51) = bulge2(length(bulge2),20);
    bulge2(length(bulge2),52) = bulge2(length(bulge2),41);
    norm_bulge2(length(bulge2),49) = norm_bulge2(length(bulge2),6);
    norm_bulge2(length(bulge2),66) = norm_bulge2(length(bulge2),15);
    norm_bulge2(length(bulge2),50) = norm_bulge2(length(bulge2),19);
    norm_bulge2(length(bulge2),51) = norm_bulge2(length(bulge2),20);
    norm_bulge2(length(bulge2),52) = norm_bulge2(length(bulge2),41);
end

i=1;
while(i<length(is3))
    if(is3(i,1) == is3(i+1,1) & is3(i,2) > is3(i+1,2))   % multiple snapshot clump
        j = i + 1;
        avg_tff(1) = is3(i,6).*is3(i,29);
        avg_tff(2) = is3(i,6);
        while(is3(i,1)==is3(j,1) & is3(i,2) > is3(j,2) & j<length(is3))
            j = j + 1;
            avg_tff(1) = avg_tff(1) + is3(j,6).*is3(j,29);
            avg_tff(2) = avg_tff(2) + is3(j,6);
        end
        if(j==length(is3) & is3(i,1)==is3(j,1) & is3(i,2) > is3(j,2))
            j = j + 1;
            avg_tff(1) = avg_tff(1) + is3(j,6).*is3(j,29);
            avg_tff(2) = avg_tff(2) + is3(j,6);
        end
        avg_tff(1) = avg_tff(1) ./ avg_tff(2);
        is3(i:(j-1),49) = max(is3(i:(j-1),6));
        is3(i:(j-1),66) = max(is3(i:(j-1),15));
        is3(i:(j-1),50) = max(is3(i:(j-1),19));
        is3(i:(j-1),51) = max(is3(i:(j-1),20));
        is3(i:(j-1),52) = is3(j-1,41);
        norm_is3(i:(j-1),49) = max(norm_is3(i:(j-1),6));
        norm_is3(i:(j-1),66) = max(norm_is3(i:(j-1),15));
        norm_is3(i:(j-1),50) = max(norm_is3(i:(j-1),19));
        norm_is3(i:(j-1),51) = max(norm_is3(i:(j-1),20));
        norm_is3(i:(j-1),52) = is3(j-1,41)./avg_tff(1);
        i = j;
    else
        is3(i,49) = is3(i,6);
        is3(i,66) = is3(i,15);
        is3(i,50) = is3(i,19);
        is3(i,51) = is3(i,20);
        is3(i,52) = is3(i,41);
        norm_is3(i,49) = norm_is3(i,6);
        norm_is3(i,66) = norm_is3(i,15);
        norm_is3(i,50) = norm_is3(i,19);
        norm_is3(i,51) = norm_is3(i,20);
        norm_is3(i,52) = norm_is3(i,41);
        i = i + 1;
    end
end
if(i==length(is3))
    is3(length(is3),49) = is3(length(is3),6);
    is3(length(is3),66) = is3(length(is3),15);
    is3(length(is3),50) = is3(length(is3),19);
    is3(length(is3),51) = is3(length(is3),20);
    is3(length(is3),52) = is3(length(is3),41);
    norm_is3(length(is3),49) = norm_is3(length(is3),6);
    norm_is3(length(is3),66) = norm_is3(length(is3),15);
    norm_is3(length(is3),50) = norm_is3(length(is3),19);
    norm_is3(length(is3),51) = norm_is3(length(is3),20);
    norm_is3(length(is3),52) = norm_is3(length(is3),41);
end
[max(norm_is3(:,41)),min(norm_is3(:,41)),length(find(isnan(norm_is3(:,41))))]
[max(norm_is3(:,52)),min(norm_is3(:,52)),length(find(isnan(norm_is3(:,52))))]
[max(is3(:,41)),min(is3(:,41)),length(find(isnan(is3(:,41))))]
[max(is3(:,52)),min(is3(:,52)),length(find(isnan(is3(:,52))))]

i=1;
while(i<length(es3))
    if(es3(i,1) == es3(i+1,1) & es3(i,2) > es3(i+1,2))   % multiple snapshot clump
        j = i + 1;
        avg_tff(1) = es3(i,6).*es3(i,29);
        avg_tff(2) = es3(i,6);
        while(es3(i,1)==es3(j,1) & es3(i,2) > es3(j,2) & j<length(es3))
            j = j + 1;
            avg_tff(1) = avg_tff(1) + es3(j,6).*es3(j,29);
            avg_tff(2) = avg_tff(2) + es3(j,6);
        end
        if(j==length(es3) & es3(i,1)==es3(j,1) & es3(i,2) > es3(j,2))
            j = j + 1;
            avg_tff(1) = avg_tff(1) + es3(j,6).*es3(j,29);
            avg_tff(2) = avg_tff(2) + es3(j,6);
        end
        avg_tff(1) = avg_tff(1) ./ avg_tff(2);
        es3(i:(j-1),49) = max(es3(i:(j-1),6));
        es3(i:(j-1),66) = max(es3(i:(j-1),15));
        es3(i:(j-1),50) = max(es3(i:(j-1),19));
        es3(i:(j-1),51) = max(es3(i:(j-1),20));
        es3(i:(j-1),52) = es3(j-1,41);
        norm_es3(i:(j-1),49) = max(norm_es3(i:(j-1),6));
        norm_es3(i:(j-1),66) = max(norm_es3(i:(j-1),15));
        norm_es3(i:(j-1),50) = max(norm_es3(i:(j-1),19));
        norm_es3(i:(j-1),51) = max(norm_es3(i:(j-1),20));
        norm_es3(i:(j-1),52) = es3(j-1,41)./avg_tff(1);
        i = j;
    else
        es3(i,49) = es3(i,6);
        es3(i,66) = es3(i,15);
        es3(i,50) = es3(i,19);
        es3(i,51) = es3(i,20);
        es3(i,52) = es3(i,41);
        norm_es3(i,49) = norm_es3(i,6);
        norm_es3(i,66) = norm_es3(i,15);
        norm_es3(i,50) = norm_es3(i,19);
        norm_es3(i,51) = norm_es3(i,20);
        norm_es3(i,52) = norm_es3(i,41);
        i = i + 1;
    end
end
if(i==length(es3))
    es3(length(es3),49) = es3(length(es3),6);
    es3(length(es3),66) = es3(length(es3),15);
    es3(length(es3),50) = es3(length(es3),19);
    es3(length(es3),51) = es3(length(es3),20);
    es3(length(es3),52) = es3(length(es3),41);
    norm_es3(length(es3),49) = norm_es3(length(es3),6);
    norm_es3(length(es3),66) = norm_es3(length(es3),15);
    norm_es3(length(es3),50) = norm_es3(length(es3),19);
    norm_es3(length(es3),51) = norm_es3(length(es3),20);
    norm_es3(length(es3),52) = norm_es3(length(es3),41);
end

i=1;
while(i<length(bulge3))
    if(bulge3(i,1) == bulge3(i+1,1) & bulge3(i,2) > bulge3(i+1,2))   % multiple snapshot clump
        j = i + 1;
        avg_tff(1) = bulge3(i,6).*bulge3(i,29);
        avg_tff(2) = bulge3(i,6);
        while(bulge3(i,1)==bulge3(j,1) & bulge3(i,2) > bulge3(j,2) & j<length(bulge3))
            j = j + 1;
            avg_tff(1) = avg_tff(1) + bulge3(j,6).*bulge3(j,29);
            avg_tff(2) = avg_tff(2) + bulge3(j,6);
        end
        if(j==length(bulge3) & bulge3(i,1)==bulge3(j,1) & bulge3(i,2) > bulge3(j,2))
            j = j + 1;
            avg_tff(1) = avg_tff(1) + bulge3(j,6).*bulge3(j,29);
            avg_tff(2) = avg_tff(2) + bulge3(j,6);
        end
        avg_tff(1) = avg_tff(1) ./ avg_tff(2);
        bulge3(i:(j-1),49) = max(bulge3(i:(j-1),6));
        bulge3(i:(j-1),66) = max(bulge3(i:(j-1),15));
        bulge3(i:(j-1),50) = max(bulge3(i:(j-1),19));
        bulge3(i:(j-1),51) = max(bulge3(i:(j-1),20));
        bulge3(i:(j-1),52) = bulge3(j-1,41);
        norm_bulge3(i:(j-1),49) = max(norm_bulge3(i:(j-1),6));
        norm_bulge3(i:(j-1),66) = max(norm_bulge3(i:(j-1),15));
        norm_bulge3(i:(j-1),50) = max(norm_bulge3(i:(j-1),19));
        norm_bulge3(i:(j-1),51) = max(norm_bulge3(i:(j-1),20));
        norm_bulge3(i:(j-1),52) = bulge3(j-1,41)./avg_tff(1);
        i = j;
    else
        bulge3(i,49) = bulge3(i,6);
        bulge3(i,66) = bulge3(i,15);
        bulge3(i,50) = bulge3(i,19);
        bulge3(i,51) = bulge3(i,20);
        bulge3(i,52) = bulge3(i,41);
        norm_bulge3(i,49) = norm_bulge3(i,6);
        norm_bulge3(i,66) = norm_bulge3(i,15);
        norm_bulge3(i,50) = norm_bulge3(i,19);
        norm_bulge3(i,51) = norm_bulge3(i,20);
        norm_bulge3(i,52) = norm_bulge3(i,41);
        i = i + 1;
    end
end
if(i==length(bulge3))
    bulge3(length(bulge3),49) = bulge3(length(bulge3),6);
    bulge3(length(bulge3),66) = bulge3(length(bulge3),15);
    bulge3(length(bulge3),50) = bulge3(length(bulge3),19);
    bulge3(length(bulge3),51) = bulge3(length(bulge3),20);
    bulge3(length(bulge3),52) = bulge3(length(bulge3),41);
    norm_bulge3(length(bulge3),49) = norm_bulge3(length(bulge3),6);
    norm_bulge3(length(bulge3),66) = norm_bulge3(length(bulge3),15);
    norm_bulge3(length(bulge3),50) = norm_bulge3(length(bulge3),19);
    norm_bulge3(length(bulge3),51) = norm_bulge3(length(bulge3),20);
    norm_bulge3(length(bulge3),52) = norm_bulge3(length(bulge3),41);
end

prop_list(49) = {'Max Mass'};
prop_list(50) = {'Max normalized distance'};
prop_list(51) = {'Max normalized height'};
prop_list(52) = {'Max age'};
unit_list(49) = {'Msun'};
unit_list(50) = {'1'};
unit_list(51) = {'1'};
unit_list(52) = {'Myr'};

clear cm3
clear i j

% Add age of stars in first snapshot to lifetime
for i=1:(length(nis3)-1)
    i
    j=1;
    while(j<=(nis3(i+1)-nis3(i)))
        if(is3(nis3(i)+j,41)==0)
            z = disc3(ndisc3(i)+1:ndisc3(i+1),1);
            [z1, I] = min( abs( z - is3(nis3(i)+j,2) ) );
            z1 = z(I);
            if(I>1)
                z2 = z(I-1);
            else
                z2 = z(I+1);
            end
            t1 = 0.95./( ((1+z1)/7)^1.5 );
            t2 = 0.95./( ((1+z2)/7)^1.5 );
            delt = 1000.*abs(t1-t2); %Myr
            if(is3(nis3(i)+j,12)>1e-6 & is3(nis3(i)+j,12)<delt)
                k = j+1;
                k1 = k;
                while(k<=(nis3(i+1)-nis3(i)))
                    if(is3(nis3(i)+k,1)==is3(nis3(i)+j,1) & is3(nis3(i)+k,2)<is3(nis3(i)+j,2))
                        k = k+1;
                        k1 = k;
                    else
                        k = 2*(nis3(i+1)-nis3(i)) + 10;
                    end
                end
                k = k1-1;
                if(k==j)
                    is3(nis3(i)+j,41) = is3(nis3(i)+j,41) + is3(nis3(i)+j,12);
                    is3(nis3(i)+j,52) = is3(nis3(i)+j,52) + is3(nis3(i)+j,12);
                    norm_is3(nis3(i)+j,41) = is3(nis3(i)+j,41)./is3(nis3(i)+j,29);
                    norm_is3(nis3(i)+j,52) = is3(nis3(i)+j,52)./is3(nis3(i)+j,29);
                else
                    temp = ( is3(nis3(i)+k,52)./norm_is3(nis3(i)+k,52) );
                    is3(nis3(i)+j:nis3(i)+k,41) = is3(nis3(i)+j:nis3(i)+k,41) + is3(nis3(i)+j,12);
                    is3(nis3(i)+j:nis3(i)+k,52) = is3(nis3(i)+j:nis3(i)+k,52) + is3(nis3(i)+j,12);
                    norm_is3(nis3(i)+j:nis3(i)+k,41) = is3(nis3(i)+j:nis3(i)+k,41)./is3(nis3(i)+j:nis3(i)+k,29);
                    norm_is3(nis3(i)+j:nis3(i)+k,52) = is3(nis3(i)+j:nis3(i)+k,52)./temp;
                end
                j = k+1;
            else
                j = j+1;
            end
        else
            j = j+1;
        end
    end
end
clear i j k k1 I z z1 z2 t1 t2 delt

% Add age of stars in first snapshot to lifetime
for i=1:(length(nis2)-1)
    i
    j=1;
    while(j<=(nis2(i+1)-nis2(i)))
        if(is2(nis2(i)+j,41)==0)
            z = disc2(ndisc2(i)+1:ndisc2(i+1),1);
            [z1, I] = min( abs( z - is2(nis2(i)+j,2) ) );
            z1 = z(I);
            if(I>1)
                z2 = z(I-1);
            else
                z2 = z(I+1);
            end
            t1 = 0.95./( ((1+z1)/7)^1.5 );
            t2 = 0.95./( ((1+z2)/7)^1.5 );
            delt = 1000.*abs(t1-t2); %Myr
            if(is2(nis2(i)+j,12)>1e-6 & is2(nis2(i)+j,12)<delt)
                k = j+1;
                k1 = k;
                while(k<=(nis2(i+1)-nis2(i)))
                    if(is2(nis2(i)+k,41)>0)
                        k = k+1;
                        k1 = k;
                    else
                        k = 2*(nis2(i+1)-nis2(i)) + 10;
                    end
                end
                k = k1-1;
                if(k==j)
                    is2(nis2(i)+j,41) = is2(nis2(i)+j,41) + is2(nis2(i)+j,12);
                    is2(nis2(i)+j,52) = is2(nis2(i)+j,52) + is2(nis2(i)+j,12);
                    norm_is2(nis2(i)+j,41) = is2(nis2(i)+j,41)./is2(nis2(i)+j,29);
                    norm_is2(nis2(i)+j,52) = is2(nis2(i)+j,52)./is2(nis2(i)+j,29);
                else
                    temp = ( is2(nis2(i)+j+1,52)./norm_is2(nis2(i)+j+1,52) );
                    is2(nis2(i)+j:nis2(i)+k,41) = is2(nis2(i)+j:nis2(i)+k,41) + is2(nis2(i)+j,12);
                    is2(nis2(i)+j:nis2(i)+k,52) = is2(nis2(i)+j:nis2(i)+k,52) + is2(nis2(i)+j,12);
                    norm_is2(nis2(i)+j:nis2(i)+k,41) = is2(nis2(i)+j:nis2(i)+k,41)./is2(nis2(i)+j:nis2(i)+k,29);
                    norm_is2(nis2(i)+j:nis2(i)+k,52) = is2(nis2(i)+j:nis2(i)+k,52)./temp;
                end
                j = k+1;
            else
                j = j+1;
            end
        else
            j = j+1;
        end
    end
end
clear i j k k1 I z z1 z2 t1 t2 delt temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prop_list(53) = {'Mc / Mcold disc'};
prop_list(54) = {'Mstar / Mstar sphere'};
prop_list(55) = {'SFR / SFR sphere'};
prop_list(56) = {'Zgas / Zgas sphere'};
prop_list(57) = {'Zstar / Zstar sphere'};
prop_list(58) = {'age / age sphere'};
prop_list(59) = {'sSFR / sSFR sphere'};
prop_list(60) = {'fgas / fgas sphere'};

unit_list(53) = {'Msun'};
unit_list(54) = {'Msun'};
unit_list(55) = {'Msun/yr'};
unit_list(56) = {'Log O/H + 12'};
unit_list(57) = {'Log O/H + 12'};
unit_list(58) = {'Myr'};
unit_list(59) = {'1/Gyr'};
unit_list(60) = {'1'};

is2(:,53)    = is2(:,6);
es2(:,53)    = es2(:,6);
bulge2(:,53) = bulge2(:,6);
is2(:,54)    = is2(:,5);
es2(:,54)    = es2(:,5);
bulge2(:,54) = bulge2(:,5);
is2(:,55)    = is2(:,15);
es2(:,55)    = es2(:,15);
bulge2(:,55) = bulge2(:,15);
is2(:,56)    = is2(:,13);
es2(:,56)    = es2(:,13);
bulge2(:,56) = bulge2(:,13);
is2(:,57)    = is2(:,14);
es2(:,57)    = es2(:,14);
bulge2(:,57) = bulge2(:,14);
is2(:,58)    = is2(:,12);
es2(:,58)    = es2(:,12);
bulge2(:,58) = bulge2(:,12);
is2(:,59)    = is2(:,17);
es2(:,59)    = es2(:,17);
bulge2(:,59) = bulge2(:,17);
is2(:,60)    = is2(:,7);
es2(:,60)    = es2(:,7);
bulge2(:,60) = bulge2(:,7);

is3(:,53)    = is3(:,6);
es3(:,53)    = es3(:,6);
bulge3(:,53) = bulge3(:,6);
is3(:,54)    = is3(:,5);
es3(:,54)    = es3(:,5);
bulge3(:,54) = bulge3(:,5);
is3(:,55)    = is3(:,15);
es3(:,55)    = es3(:,15);
bulge3(:,55) = bulge3(:,15);
is3(:,56)    = is3(:,13);
es3(:,56)    = es3(:,13);
bulge3(:,56) = bulge3(:,13);
is3(:,57)    = is3(:,14);
es3(:,57)    = es3(:,14);
bulge3(:,57) = bulge3(:,14);
is3(:,58)    = is3(:,12);
es3(:,58)    = es3(:,12);
bulge3(:,58) = bulge3(:,12);
is3(:,59)    = is3(:,17);
es3(:,59)    = es3(:,17);
bulge3(:,59) = bulge3(:,17);
is3(:,60)    = is3(:,7);
es3(:,60)    = es3(:,7);
bulge3(:,60) = bulge3(:,7);

norm_is2(:,53)    = norm_is2(:,6);
norm_es2(:,53)    = norm_es2(:,6);
norm_bulge2(:,53) = norm_bulge2(:,6);
norm_is2(:,54)    = norm_is2(:,5);
norm_es2(:,54)    = norm_es2(:,5);
norm_bulge2(:,54) = norm_bulge2(:,5);
norm_is2(:,55)    = norm_is2(:,15);
norm_es2(:,55)    = norm_es2(:,15);
norm_bulge2(:,55) = norm_bulge2(:,15);
norm_is2(:,56)    = norm_is2(:,13);
norm_es2(:,56)    = norm_es2(:,13);
norm_bulge2(:,56) = norm_bulge2(:,13);
norm_is2(:,57)    = norm_is2(:,14);
norm_es2(:,57)    = norm_es2(:,14);
norm_bulge2(:,57) = norm_bulge2(:,14);
norm_is2(:,58)    = norm_is2(:,12);
norm_es2(:,58)    = norm_es2(:,12);
norm_bulge2(:,58) = norm_bulge2(:,12);
norm_is2(:,59)    = norm_is2(:,17);
norm_es2(:,59)    = norm_es2(:,17);
norm_bulge2(:,59) = norm_bulge2(:,17);
norm_is2(:,60)    = norm_is2(:,7);
norm_es2(:,60)    = norm_es2(:,7);
norm_bulge2(:,60) = norm_bulge2(:,7);

norm_is3(:,53)    = norm_is3(:,6);
norm_es3(:,53)    = norm_es3(:,6);
norm_bulge3(:,53) = norm_bulge3(:,6);
norm_is3(:,54)    = norm_is3(:,5);
norm_es3(:,54)    = norm_es3(:,5);
norm_bulge3(:,54) = norm_bulge3(:,5);
norm_is3(:,55)    = norm_is3(:,15);
norm_es3(:,55)    = norm_es3(:,15);
norm_bulge3(:,55) = norm_bulge3(:,15);
norm_is3(:,56)    = norm_is3(:,13);
norm_es3(:,56)    = norm_es3(:,13);
norm_bulge3(:,56) = norm_bulge3(:,13);
norm_is3(:,57)    = norm_is3(:,14);
norm_es3(:,57)    = norm_es3(:,14);
norm_bulge3(:,57) = norm_bulge3(:,14);
norm_is3(:,58)    = norm_is3(:,12);
norm_es3(:,58)    = norm_es3(:,12);
norm_bulge3(:,58) = norm_bulge3(:,12);
norm_is3(:,59)    = norm_is3(:,17);
norm_es3(:,59)    = norm_es3(:,17);
norm_bulge3(:,59) = norm_bulge3(:,17);
norm_is3(:,60)    = norm_is3(:,7);
norm_es3(:,60)    = norm_es3(:,7);
norm_bulge3(:,60) = norm_bulge3(:,7);

for i=1:length(ind2)
    b = (N_cat2(i)+1):N_cat2(i+1);
    a = disc_cat2(b,1);
    for j=(nis2(i)+1):nis2(i+1)
        z2 = is2(j,2);
        a2 = 1./(1+z2);
        k = find( abs(a-a2)<=0.003 );
        k = k(1);
        norm_is2(j,53) = norm_is2(j,53) .*( disc_cat2(b(k),14)  + disc_cat2(b(k),16) )./ disc_cat2(b(k),15);
        norm_is2(j,54) = norm_is2(j,54) .*  disc_cat2(b(k),16) ./ sphere_cat2(b(k),11);
        norm_is2(j,55) = norm_is2(j,55) .*  disc_cat2(b(k),19) ./ sphere_cat2(b(k),14);
        norm_is2(j,56) = norm_is2(j,56) .*  disc_cat2(b(k),21) ./ sphere_cat2(b(k),16);
        norm_is2(j,57) = norm_is2(j,57) .*  disc_cat2(b(k),22) ./ sphere_cat2(b(k),17);
        norm_is2(j,58) = norm_is2(j,58) .*  disc_cat2(b(k),20) ./ sphere_cat2(b(k),15);
        norm_is2(j,60) = norm_is2(j,60) .*( disc_cat2(b(k),14) ./(disc_cat2(b(k),14)+disc_cat2(b(k),16)) ) ./ ...
            ( sphere_cat2(b(k),9)./(sphere_cat2(b(k),9)+sphere_cat2(b(k),11)) );
    end
    for j=(nes2(i)+1):nes2(i+1)
        z2 = es2(j,2);
        a2 = 1./(1+z2);
        k = find( abs(a-a2)<=0.003 );
        k = k(1);
        norm_es2(j,53) = norm_es2(j,53) .*( disc_cat2(b(k),14)  + disc_cat2(b(k),16) )./ disc_cat2(b(k),15);
        norm_es2(j,54) = norm_es2(j,54) .*  disc_cat2(b(k),16) ./ sphere_cat2(b(k),11);
        norm_es2(j,55) = norm_es2(j,55) .*  disc_cat2(b(k),19) ./ sphere_cat2(b(k),14);
        norm_es2(j,56) = norm_es2(j,56) .*  disc_cat2(b(k),21) ./ sphere_cat2(b(k),16);
        norm_es2(j,57) = norm_es2(j,57) .*  disc_cat2(b(k),22) ./ sphere_cat2(b(k),17);
        norm_es2(j,58) = norm_es2(j,58) .*  disc_cat2(b(k),20) ./ sphere_cat2(b(k),15);
        norm_es2(j,60) = norm_es2(j,60) .*( disc_cat2(b(k),14) ./(disc_cat2(b(k),14)+disc_cat2(b(k),16)) ) ./ ...
            ( sphere_cat2(b(k),9)./(sphere_cat2(b(k),9)+sphere_cat2(b(k),11)) );
    end
    for j=(nbulge2(i)+1):nbulge2(i+1)
        z2 = bulge2(j,2);
        a2 = 1./(1+z2);
        k = find( abs(a-a2)<=0.003 );
        k = k(1);
        norm_bulge2(j,53) = norm_bulge2(j,53) .*( disc_cat2(b(k),14)  + disc_cat2(b(k),16) )./ disc_cat2(b(k),15);
        norm_bulge2(j,54) = norm_bulge2(j,54) .*  disc_cat2(b(k),16) ./ sphere_cat2(b(k),11);
        norm_bulge2(j,55) = norm_bulge2(j,55) .*  disc_cat2(b(k),19) ./ sphere_cat2(b(k),14);
        norm_bulge2(j,56) = norm_bulge2(j,56) .*  disc_cat2(b(k),21) ./ sphere_cat2(b(k),16);
        norm_bulge2(j,57) = norm_bulge2(j,57) .*  disc_cat2(b(k),22) ./ sphere_cat2(b(k),17);
        norm_bulge2(j,58) = norm_bulge2(j,58) .*  disc_cat2(b(k),20) ./ sphere_cat2(b(k),15);
        norm_bulge2(j,60) = norm_bulge2(j,60) .*( disc_cat2(b(k),14) ./(disc_cat2(b(k),14)+disc_cat2(b(k),16)) ) ./ ...
            ( sphere_cat2(b(k),9)./(sphere_cat2(b(k),9)+sphere_cat2(b(k),11)) );
    end
end
norm_is2(:,59)    = norm_is2(:,54)   ./norm_is2(:,55);
norm_es2(:,59)    = norm_es2(:,54)   ./norm_es2(:,55);
norm_bulge2(:,59) = norm_bulge2(:,54)./norm_bulge2(:,55);
for i=1:length(ind3)
    b = (N_cat3(i)+1):N_cat3(i+1);
    a = disc_cat3(b,1);
    for j=(nis3(i)+1):nis3(i+1)
        z3 = is3(j,2);
        a3 = 1./(1+z3);
        k = find( abs(a-a3)<=0.003 );
        k = k(1);
        norm_is3(j,53) = norm_is3(j,53) .*( disc_cat3(b(k),14)  + disc_cat3(b(k),16) )./ disc_cat3(b(k),15);
        norm_is3(j,54) = norm_is3(j,54) .*  disc_cat3(b(k),16) ./ sphere_cat3(b(k),11);
        norm_is3(j,55) = norm_is3(j,55) .*  disc_cat3(b(k),19) ./ sphere_cat3(b(k),14);
        norm_is3(j,56) = norm_is3(j,56) .*  disc_cat3(b(k),21) ./ sphere_cat3(b(k),16);
        norm_is3(j,57) = norm_is3(j,57) .*  disc_cat3(b(k),22) ./ sphere_cat3(b(k),17);
        norm_is3(j,58) = norm_is3(j,58) .*  disc_cat3(b(k),20) ./ sphere_cat3(b(k),15);
        norm_is3(j,60) = norm_is3(j,60) .*( disc_cat3(b(k),14) ./(disc_cat3(b(k),14)+disc_cat3(b(k),16)) ) ./ ...
            ( sphere_cat3(b(k),9)./(sphere_cat3(b(k),9)+sphere_cat3(b(k),11)) );
    end
    for j=(nes3(i)+1):nes3(i+1)
        z3 = es3(j,2);
        a3 = 1./(1+z3);
        k = find( abs(a-a3)<=0.003 );
        k = k(1);
        norm_es3(j,53) = norm_es3(j,53) .*( disc_cat3(b(k),14)  + disc_cat3(b(k),16) )./ disc_cat3(b(k),15);
        norm_es3(j,54) = norm_es3(j,54) .*  disc_cat3(b(k),16) ./ sphere_cat3(b(k),11);
        norm_es3(j,55) = norm_es3(j,55) .*  disc_cat3(b(k),19) ./ sphere_cat3(b(k),14);
        norm_es3(j,56) = norm_es3(j,56) .*  disc_cat3(b(k),21) ./ sphere_cat3(b(k),16);
        norm_es3(j,57) = norm_es3(j,57) .*  disc_cat3(b(k),22) ./ sphere_cat3(b(k),17);
        norm_es3(j,58) = norm_es3(j,58) .*  disc_cat3(b(k),20) ./ sphere_cat3(b(k),15);
        norm_es3(j,60) = norm_es3(j,60) .*( disc_cat3(b(k),14) ./(disc_cat3(b(k),14)+disc_cat3(b(k),16)) ) ./ ...
            ( sphere_cat3(b(k),9)./(sphere_cat3(b(k),9)+sphere_cat3(b(k),11)) );
    end
    for j=(nbulge3(i)+1):nbulge3(i+1)
        z3 = bulge3(j,2);
        a3 = 1./(1+z3);
        k = find( abs(a-a3)<=0.003 );
        k = k(1);
        norm_bulge3(j,53) = norm_bulge3(j,53) .*( disc_cat3(b(k),14)  + disc_cat3(b(k),16) )./ disc_cat3(b(k),15);
        norm_bulge3(j,54) = norm_bulge3(j,54) .*  disc_cat3(b(k),16) ./ sphere_cat3(b(k),11);
        norm_bulge3(j,55) = norm_bulge3(j,55) .*  disc_cat3(b(k),19) ./ sphere_cat3(b(k),14);
        norm_bulge3(j,56) = norm_bulge3(j,56) .*  disc_cat3(b(k),21) ./ sphere_cat3(b(k),16);
        norm_bulge3(j,57) = norm_bulge3(j,57) .*  disc_cat3(b(k),22) ./ sphere_cat3(b(k),17);
        norm_bulge3(j,58) = norm_bulge3(j,58) .*  disc_cat3(b(k),20) ./ sphere_cat3(b(k),15);
        norm_bulge3(j,60) = norm_bulge3(j,60) .*( disc_cat3(b(k),14) ./(disc_cat3(b(k),14)+disc_cat3(b(k),16)) ) ./ ...
            ( sphere_cat3(b(k),9)./(sphere_cat3(b(k),9)+sphere_cat3(b(k),11)) );
    end
end
norm_is3(:,59)    = norm_is3(:,54)   ./norm_is3(:,55);
norm_es3(:,59)    = norm_es3(:,54)   ./norm_es3(:,55);
norm_bulge3(:,59) = norm_bulge3(:,54)./norm_bulge3(:,55);
clear a a2 a3 b z2 z3 i j k

% t over global td
is2(:,61)         = is2(:,41)./is2(:,31);
es2(:,61)         = es2(:,41)./es2(:,31);
bulge2(:,61)      = bulge2(:,41)./bulge2(:,31);
norm_is2(:,61)    = is2(:,52)./is2(:,31);
norm_es2(:,61)    = es2(:,52)./es2(:,31);
norm_bulge2(:,61) = bulge2(:,52)./bulge2(:,31);
is3(:,61)         = is3(:,41)./is3(:,31);
es3(:,61)         = es3(:,41)./es3(:,31);
bulge3(:,61)      = bulge3(:,41)./bulge3(:,31);
norm_is3(:,61)    = is3(:,52)./is3(:,31);
norm_es3(:,61)    = es3(:,52)./es3(:,31);
norm_bulge3(:,61) = bulge3(:,52)./bulge3(:,31);
prop_list(61)     = {'td local over tff'};
unit_list(61)     = {'1'};

% Vcirc
is2(:,62)         = sqrt(4.3e-3.*is2(:,6)./(1000.*is2(:,3)));
es2(:,62)         = sqrt(4.3e-3.*es2(:,6)./(1000.*es2(:,3)));
bulge2(:,62)      = sqrt(4.3e-3.*bulge2(:,6)./(1000.*bulge2(:,3)));
norm_is2(:,62)    = sqrt(4.3e-3.*is2(:,6)./(1000.*is2(:,3)));
norm_es2(:,62)    = sqrt(4.3e-3.*es2(:,6)./(1000.*es2(:,3)));
norm_bulge2(:,62) = sqrt(4.3e-3.*bulge2(:,6)./(1000.*bulge2(:,3)));
is3(:,62)         = sqrt(4.3e-3.*is3(:,6)./(1000.*is3(:,3)));
es3(:,62)         = sqrt(4.3e-3.*es3(:,6)./(1000.*es3(:,3)));
bulge3(:,62)      = sqrt(4.3e-3.*bulge3(:,6)./(1000.*bulge3(:,3)));
norm_is3(:,62)    = sqrt(4.3e-3.*is3(:,6)./(1000.*is3(:,3)));
norm_es3(:,62)    = sqrt(4.3e-3.*es3(:,6)./(1000.*es3(:,3)));
norm_bulge3(:,62) = sqrt(4.3e-3.*bulge3(:,6)./(1000.*bulge3(:,3)));
prop_list(62)     = {'Vcirc'};
unit_list(62)     = {'km/s'};

% eps_ff
is2(:,63)         = is2(:,15)./( (is2(:,4)+is2(:,15).*60.*1e6)./(1e6.*is2(:,29)) );
es2(:,63)         = es2(:,15)./( (es2(:,4)+es2(:,15).*60.*1e6)./(1e6.*es2(:,29)) );
bulge2(:,63)      = bulge2(:,15)./( (bulge2(:,4)+bulge2(:,15).*60.*1e6)./(1e6.*bulge2(:,29)) );
norm_is2(:,63)    = norm_is2(:,15)./( norm_is2(:,4)./(1./is2(:,43)) );
norm_es2(:,63)    = norm_es2(:,15)./( norm_es2(:,4)./(1./es2(:,43)) );
norm_bulge2(:,63) = norm_bulge2(:,15)./( norm_bulge2(:,4)./(1./bulge2(:,43)) );
is3(:,63)         = is3(:,15)./( (is3(:,4)+is3(:,15).*60.*1e6)./(1e6.*is3(:,29)) );
es3(:,63)         = es3(:,15)./( (es3(:,4)+es3(:,15).*60.*1e6)./(1e6.*es3(:,29)) );
bulge3(:,63)      = bulge3(:,15)./( (bulge3(:,4)+bulge3(:,15).*60.*1e6)./(1e6.*bulge3(:,29)) );
norm_is3(:,63)    = norm_is3(:,15)./( norm_is3(:,4)./(1./is3(:,43)) );
norm_es3(:,63)    = norm_es3(:,15)./( norm_es3(:,4)./(1./es3(:,43)) );
norm_bulge3(:,63) = norm_bulge3(:,15)./( norm_bulge3(:,4)./(1./bulge3(:,43)) );
prop_list(63)     = {'EpsFF'};
unit_list(63)     = {'1'};

% mass_loading
is2(:,64)         = is2(:,45)./is2(:,15);
es2(:,64)         = es2(:,45)./es2(:,15);
bulge2(:,64)      = bulge2(:,45)./bulge2(:,15);
norm_is2(:,64)    = is2(:,45)./is2(:,15);
norm_es2(:,64)    = es2(:,45)./es2(:,15);
norm_bulge2(:,64) = bulge2(:,45)./bulge2(:,15);
is3(:,64)         = is3(:,45)./is3(:,15);
es3(:,64)         = es3(:,45)./es3(:,15);
bulge3(:,64)      = bulge3(:,45)./bulge3(:,15);
norm_is3(:,64)    = is3(:,45)./is3(:,15);
norm_es3(:,64)    = es3(:,45)./es3(:,15);
norm_bulge3(:,64) = bulge3(:,45)./bulge3(:,15);
prop_list(64)     = {'Mass Loading'};
unit_list(64)     = {'1'};

% net in_over_out
is2(:,65)         = is2(:,44)./(is2(:,45)+is2(:,15));
es2(:,65)         = es2(:,44)./(es2(:,45)+es2(:,15));
bulge2(:,65)      = bulge2(:,44)./(bulge2(:,45)+bulge2(:,15));
norm_is2(:,65)    = is2(:,44)./(is2(:,45)+is2(:,15));
norm_es2(:,65)    = es2(:,44)./(es2(:,45)+es2(:,15));
norm_bulge2(:,65) = bulge2(:,44)./(bulge2(:,45)+bulge2(:,15));
is3(:,65)         = is3(:,44)./(is3(:,45)+is3(:,15));
es3(:,65)         = es3(:,44)./(es3(:,45)+es3(:,15));
bulge3(:,65)      = bulge3(:,44)./(bulge3(:,45)+bulge3(:,15));
norm_is3(:,65)    = is3(:,44)./(is3(:,45)+is3(:,15));
norm_es3(:,65)    = es3(:,44)./(es3(:,45)+es3(:,15));
norm_bulge3(:,65) = bulge3(:,44)./(bulge3(:,45)+bulge3(:,15));
prop_list(65)     = {'net in over out'};
unit_list(65)     = {'1'};

prop_list(66) = {'Max SFR'};
unit_list(66) = {'Msun/yr'};

% out_over_in
is2(:,67)         = is2(:,45)./(is2(:,44));
es2(:,67)         = es2(:,45)./(es2(:,44));
bulge2(:,67)      = bulge2(:,45)./(bulge2(:,44));
norm_is2(:,67)    = is2(:,45)./(is2(:,44));
norm_es2(:,67)    = es2(:,45)./(es2(:,44));
norm_bulge2(:,67) = bulge2(:,45)./(bulge2(:,44));
is3(:,67)         = is3(:,45)./(is3(:,44));
es3(:,67)         = es3(:,45)./(es3(:,44));
bulge3(:,67)      = bulge3(:,45)./(bulge3(:,44));
norm_is3(:,67)    = is3(:,45)./(is3(:,44));
norm_es3(:,67)    = es3(:,45)./(es3(:,44));
norm_bulge3(:,67) = bulge3(:,45)./(bulge3(:,44));
prop_list(67)     = {'out over in'};
unit_list(67)     = {'1'};


% mass_loading
is2(:,68)         = is2(:,35)./is2(:,15);
es2(:,68)         = es2(:,35)./es2(:,15);
bulge2(:,68)      = bulge2(:,35)./bulge2(:,15);
norm_is2(:,68)    = is2(:,35)./is2(:,15);
norm_es2(:,68)    = es2(:,35)./es2(:,15);
norm_bulge2(:,68) = bulge2(:,35)./bulge2(:,15);
is3(:,68)         = is3(:,35)./is3(:,15);
es3(:,68)         = es3(:,35)./es3(:,15);
bulge3(:,68)      = bulge3(:,35)./bulge3(:,15);
norm_is3(:,68)    = is3(:,35)./is3(:,15);
norm_es3(:,68)    = es3(:,35)./es3(:,15);
norm_bulge3(:,68) = bulge3(:,35)./bulge3(:,15);
prop_list(68)     = {'Mass Loading Rc'};
unit_list(68)     = {'1'};

% mass_loading
is2(:,69)         = is2(:,36)./is2(:,15);
es2(:,69)         = es2(:,36)./es2(:,15);
bulge2(:,69)      = bulge2(:,36)./bulge2(:,15);
norm_is2(:,69)    = is2(:,36)./is2(:,15);
norm_es2(:,69)    = es2(:,36)./es2(:,15);
norm_bulge2(:,69) = bulge2(:,36)./bulge2(:,15);
is3(:,69)         = is3(:,36)./is3(:,15);
es3(:,69)         = es3(:,36)./es3(:,15);
bulge3(:,69)      = bulge3(:,36)./bulge3(:,15);
norm_is3(:,69)    = is3(:,36)./is3(:,15);
norm_es3(:,69)    = es3(:,36)./es3(:,15);
norm_bulge3(:,69) = bulge3(:,36)./bulge3(:,15);
prop_list(69)     = {'Mass Loading 1.5Rc'};
unit_list(69)     = {'1'};

% mass_loading
is2(:,70)         = is2(:,37)./is2(:,15);
es2(:,70)         = es2(:,37)./es2(:,15);
bulge2(:,70)      = bulge2(:,37)./bulge2(:,15);
norm_is2(:,70)    = is2(:,37)./is2(:,15);
norm_es2(:,70)    = es2(:,37)./es2(:,15);
norm_bulge2(:,70) = bulge2(:,37)./bulge2(:,15);
is3(:,70)         = is3(:,37)./is3(:,15);
es3(:,70)         = es3(:,37)./es3(:,15);
bulge3(:,70)      = bulge3(:,37)./bulge3(:,15);
norm_is3(:,70)    = is3(:,37)./is3(:,15);
norm_es3(:,70)    = es3(:,37)./es3(:,15);
norm_bulge3(:,70) = bulge3(:,37)./bulge3(:,15);
prop_list(70)     = {'Mass Loading 2Rc'};
unit_list(70)     = {'1'};

% out_over_in
is2(:,71)         = is2(:,35)./(is2(:,32));
es2(:,71)         = es2(:,35)./(es2(:,32));
bulge2(:,71)      = bulge2(:,35)./(bulge2(:,32));
norm_is2(:,71)    = is2(:,35)./(is2(:,32));
norm_es2(:,71)    = es2(:,35)./(es2(:,32));
norm_bulge2(:,71) = bulge2(:,35)./(bulge2(:,32));
is3(:,71)         = is3(:,35)./(is3(:,32));
es3(:,71)         = es3(:,35)./(es3(:,32));
bulge3(:,71)      = bulge3(:,35)./(bulge3(:,32));
norm_is3(:,71)    = is3(:,35)./(is3(:,32));
norm_es3(:,71)    = es3(:,35)./(es3(:,32));
norm_bulge3(:,71) = bulge3(:,35)./(bulge3(:,32));
prop_list(71)     = {'out over in Rc'};
unit_list(71)     = {'1'};

% out_over_in
is2(:,72)         = is2(:,36)./(is2(:,33));
es2(:,72)         = es2(:,36)./(es2(:,33));
bulge2(:,72)      = bulge2(:,36)./(bulge2(:,33));
norm_is2(:,72)    = is2(:,36)./(is2(:,33));
norm_es2(:,72)    = es2(:,36)./(es2(:,33));
norm_bulge2(:,72) = bulge2(:,36)./(bulge2(:,33));
is3(:,72)         = is3(:,36)./(is3(:,33));
es3(:,72)         = es3(:,36)./(es3(:,33));
bulge3(:,72)      = bulge3(:,36)./(bulge3(:,33));
norm_is3(:,72)    = is3(:,36)./(is3(:,33));
norm_es3(:,72)    = es3(:,36)./(es3(:,33));
norm_bulge3(:,72) = bulge3(:,36)./(bulge3(:,33));
prop_list(72)     = {'out over in 1.5Rc'};
unit_list(72)     = {'1'};

% out_over_in
is2(:,73)         = is2(:,37)./(is2(:,34));
es2(:,73)         = es2(:,37)./(es2(:,34));
bulge2(:,73)      = bulge2(:,37)./(bulge2(:,34));
norm_is2(:,73)    = is2(:,37)./(is2(:,34));
norm_es2(:,73)    = es2(:,37)./(es2(:,34));
norm_bulge2(:,73) = bulge2(:,37)./(bulge2(:,34));
is3(:,73)         = is3(:,37)./(is3(:,34));
es3(:,73)         = es3(:,37)./(es3(:,34));
bulge3(:,73)      = bulge3(:,37)./(bulge3(:,34));
norm_is3(:,73)    = is3(:,37)./(is3(:,34));
norm_es3(:,73)    = es3(:,37)./(es3(:,34));
norm_bulge3(:,73) = bulge3(:,37)./(bulge3(:,34));
prop_list(73)     = {'out over in 2Rc'};
unit_list(73)     = {'1'};

% mass_loading_Vesc
is2(:,74)         = is2_Mgas_out_Vesc(:,5)./is2(:,15);
es2(:,74)         = es2_Mgas_out_Vesc(:,5)./es2(:,15);
bulge2(:,74)      = bulge2_Mgas_out_Vesc(:,5)./bulge2(:,15);
norm_is2(:,74)    = is2_Mgas_out_Vesc(:,5)./is2(:,15);
norm_es2(:,74)    = es2_Mgas_out_Vesc(:,5)./es2(:,15);
norm_bulge2(:,74) = bulge2_Mgas_out_Vesc(:,5)./bulge2(:,15);
is3(:,74)         = is3_Mgas_out_Vesc(:,5)./is3(:,15);
es3(:,74)         = es3_Mgas_out_Vesc(:,5)./es3(:,15);
bulge3(:,74)      = bulge3_Mgas_out_Vesc(:,5)./bulge3(:,15);
norm_is3(:,74)    = is3_Mgas_out_Vesc(:,5)./is3(:,15);
norm_es3(:,74)    = es3_Mgas_out_Vesc(:,5)./es3(:,15);
norm_bulge3(:,74) = bulge3_Mgas_out_Vesc(:,5)./bulge3(:,15);
prop_list(74)     = {'Mass Loading Vesc Rc'};
unit_list(74)     = {'1'};

% mass_loading_Vesc
is2(:,75)         = is2_Mgas_out_Vesc(:,6)./is2(:,15);
es2(:,75)         = es2_Mgas_out_Vesc(:,6)./es2(:,15);
bulge2(:,75)      = bulge2_Mgas_out_Vesc(:,6)./bulge2(:,15);
norm_is2(:,75)    = is2_Mgas_out_Vesc(:,6)./is2(:,15);
norm_es2(:,75)    = es2_Mgas_out_Vesc(:,6)./es2(:,15);
norm_bulge2(:,75) = bulge2_Mgas_out_Vesc(:,6)./bulge2(:,15);
is3(:,75)         = is3_Mgas_out_Vesc(:,6)./is3(:,15);
es3(:,75)         = es3_Mgas_out_Vesc(:,6)./es3(:,15);
bulge3(:,75)      = bulge3_Mgas_out_Vesc(:,6)./bulge3(:,15);
norm_is3(:,75)    = is3_Mgas_out_Vesc(:,6)./is3(:,15);
norm_es3(:,75)    = es3_Mgas_out_Vesc(:,6)./es3(:,15);
norm_bulge3(:,75) = bulge3_Mgas_out_Vesc(:,6)./bulge3(:,15);
prop_list(75)     = {'Mass Loading Vesc 1.5Rc'};
unit_list(75)     = {'1'};

% mass_loading_Vesc
is2(:,76)         = is2_Mgas_out_Vesc(:,7)./is2(:,15);
es2(:,76)         = es2_Mgas_out_Vesc(:,7)./es2(:,15);
bulge2(:,76)      = bulge2_Mgas_out_Vesc(:,7)./bulge2(:,15);
norm_is2(:,76)    = is2_Mgas_out_Vesc(:,7)./is2(:,15);
norm_es2(:,76)    = es2_Mgas_out_Vesc(:,7)./es2(:,15);
norm_bulge2(:,76) = bulge2_Mgas_out_Vesc(:,7)./bulge2(:,15);
is3(:,76)         = is3_Mgas_out_Vesc(:,7)./is3(:,15);
es3(:,76)         = es3_Mgas_out_Vesc(:,7)./es3(:,15);
bulge3(:,76)      = bulge3_Mgas_out_Vesc(:,7)./bulge3(:,15);
norm_is3(:,76)    = is3_Mgas_out_Vesc(:,7)./is3(:,15);
norm_es3(:,76)    = es3_Mgas_out_Vesc(:,7)./es3(:,15);
norm_bulge3(:,76) = bulge3_Mgas_out_Vesc(:,7)./bulge3(:,15);
prop_list(76)     = {'Mass Loading Vesc 2Rc'};
unit_list(76)     = {'1'};

clear is2_Mgas_out_Vesc es2_Mgas_out_Vesc bulge2_Mgas_out_Vesc
clear is3_Mgas_out_Vesc es3_Mgas_out_Vesc bulge3_Mgas_out_Vesc

% mass_loading_Vesc
is2(:,77)         = ( is2(:,74) + is2(:,75) + is2(:,76) ) ./ 3;
es2(:,77)         = ( es2(:,74) + es2(:,75) + es2(:,76) ) ./ 3;
bulge2(:,77)      = ( bulge2(:,74) + bulge2(:,75) + bulge2(:,76) ) ./ 3;
norm_is2(:,77)    = ( is2(:,74) + is2(:,75) + is2(:,76) ) ./ 3;
norm_es2(:,77)    = ( es2(:,74) + es2(:,75) + es2(:,76) ) ./ 3;
norm_bulge2(:,77) = ( bulge2(:,74) + bulge2(:,75) + bulge2(:,76) ) ./ 3;
is3(:,77)         = ( is3(:,74) + is3(:,75) + is3(:,76) ) ./ 3;
es3(:,77)         = ( es3(:,74) + es3(:,75) + es3(:,76) ) ./ 3;
bulge3(:,77)      = ( bulge3(:,74) + bulge3(:,75) + bulge3(:,76) ) ./ 3;
norm_is3(:,77)    = ( is3(:,74) + is3(:,75) + is3(:,76) ) ./ 3;
norm_es3(:,77)    = ( es3(:,74) + es3(:,75) + es3(:,76) ) ./ 3;
norm_bulge3(:,77) = ( bulge3(:,74) + bulge3(:,75) + bulge3(:,76) ) ./ 3;
prop_list(77)     = {'Mass Loading Vesc average'};
unit_list(77)     = {'1'};

% alpha_vir
is2(:,78)         = abs(is2_alpha_vir);
es2(:,78)         = abs(es2_alpha_vir);
bulge2(:,78)      = abs(bulge2_alpha_vir);
norm_is2(:,78)    = abs(is2_alpha_vir);
norm_es2(:,78)    = abs(es2_alpha_vir);
norm_bulge2(:,78) = abs(bulge2_alpha_vir);
is3(:,78)         = abs(is3_alpha_vir);
es3(:,78)         = abs(es3_alpha_vir);
bulge3(:,78)      = abs(bulge3_alpha_vir);
norm_is3(:,78)    = abs(is3_alpha_vir);
norm_es3(:,78)    = abs(es3_alpha_vir);
norm_bulge3(:,78) = abs(bulge3_alpha_vir);
prop_list(78)     = {'alpha vir'};
unit_list(78)     = {'1'};

clear is2_alpha_vir es2_alpha_vir bulge2_alpha_vir
clear is3_alpha_vir es3_alpha_vir bulge3_alpha_vir
