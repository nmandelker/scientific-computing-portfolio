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
prop_list(24) = {'shape'};
prop_list(25) = {'DM_residual'};
prop_list(26) = {'Es'};
prop_list(27) = {'merger'};
prop_list(28) = {'tff'};
prop_list(29) = {'td_local'};
prop_list(30) = {'td_global'};
prop_list(31) = {'Mgas_in_1_0_Rc'};
prop_list(32) = {'Mgas_in_1_5_Rc'};
prop_list(33) = {'Mgas_in_2_0_Rc'};
prop_list(34) = {'Mgas_out_1_0_Rc'};
prop_list(35) = {'Mgas_out_1_5_Rc'};
prop_list(36) = {'Mgas_out_2_0_Rc'};
prop_list(37) = {'Mgas_out_Vesc_1_0_Rc'};
prop_list(38) = {'Mgas_out_Vesc_1_5_Rc'};
prop_list(39) = {'Mgas_out_Vesc_2_0_Rc'};
prop_list(40) = {'Mgas_out_Frederic'};
prop_list(41) = {'Mgas_out_Vesc_v2_1_0_Rc'};
prop_list(42) = {'Mgas_out_Vesc_v2_1_5_Rc'};
prop_list(43) = {'Mgas_out_Vesc_v2_2_0_Rc'};
prop_list(44) = {'Mstar_in'};
prop_list(45) = {'Mstar_out'};
prop_list(46) = {'Mstar_formed'};
prop_list(47) = {'alpha_vir'};
prop_list(48) = {'time'};

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
unit_list(28) = {'Myr'};
unit_list(29) = {'Myr'};
unit_list(30) = {'Myr'};
unit_list(31) = {'Msun/yr'};
unit_list(32) = {'Msun/yr'};
unit_list(33) = {'Msun/yr'};
unit_list(34) = {'Msun/yr'};
unit_list(35) = {'Msun/yr'};
unit_list(36) = {'Msun/yr'};
unit_list(37) = {'Msun/yr'};
unit_list(38) = {'Msun/yr'};
unit_list(39) = {'Msun/yr'};
unit_list(40) = {'Msun/yr'};
unit_list(41) = {'Msun/yr'};
unit_list(42) = {'Msun/yr'};
unit_list(43) = {'Msun/yr'};
unit_list(44) = {'Msun/yr'};
unit_list(45) = {'Msun/yr'};
unit_list(46) = {'Msun/yr'};
unit_list(47) = {'1'};
unit_list(48) = {'Myr'};

% local td over tff
is(:,49)         = is(:,29)./is(:,28);
es(:,49)         = es(:,29)./es(:,28);
bulge(:,49)      = bulge(:,29)./bulge(:,28);
norm_is(:,49)    = norm_is(:,29)./norm_is(:,28);
norm_es(:,49)    = norm_es(:,29)./norm_es(:,28);
norm_bulge(:,49) = norm_bulge(:,29)./norm_bulge(:,28);
prop_list(49)     = {'td local over tff'};
unit_list(49)     = {'1'};

% global td over tff
is(:,50)         = is(:,30)./is(:,28);
es(:,50)         = es(:,30)./es(:,28);
bulge(:,50)      = bulge(:,30)./bulge(:,28);
norm_is(:,50)    = norm_is(:,30)./norm_is(:,28);
norm_es(:,50)    = norm_es(:,30)./norm_es(:,28);
norm_bulge(:,50) = norm_bulge(:,30)./norm_bulge(:,28);
prop_list(50)     = {'td global over tff'};
unit_list(50)     = {'1'};

% average Mgas in Vr<0
is(:,51)         = (is(:,31)        +is(:,32)        +is(:,33))         ./ 3;
es(:,51)         = (es(:,31)        +es(:,32)        +es(:,33))         ./ 3;
bulge(:,51)      = (bulge(:,31)     +bulge(:,32)     +bulge(:,33))      ./ 3;
norm_is(:,51)    = (norm_is(:,31)   +norm_is(:,32)   +norm_is(:,33))    ./ 3;
norm_es(:,51)    = (norm_es(:,31)   +norm_es(:,32)   +norm_es(:,33))    ./ 3;
norm_bulge(:,51) = (norm_bulge(:,31)+norm_bulge(:,32)+norm_bulge(:,33)) ./ 3;
prop_list(51)     = {'average Mgas in Vr<0'};
unit_list(51)     = {'Msun/yr'};

% average Mgas out Vr>0
is(:,52)         = (is(:,34)        +is(:,35)        +is(:,36))         ./ 3;
es(:,52)         = (es(:,34)        +es(:,35)        +es(:,36))         ./ 3;
bulge(:,52)      = (bulge(:,34)     +bulge(:,35)     +bulge(:,36))      ./ 3;
norm_is(:,52)    = (norm_is(:,34)   +norm_is(:,35)   +norm_is(:,36))    ./ 3;
norm_es(:,52)    = (norm_es(:,34)   +norm_es(:,35)   +norm_es(:,36))    ./ 3;
norm_bulge(:,52) = (norm_bulge(:,34)+norm_bulge(:,35)+norm_bulge(:,36)) ./ 3;
prop_list(52)     = {'average Mgas out Vr>0'};
unit_list(52)     = {'Msun/yr'};

% average Mgas out Vr>Vesc
is(:,53)         = (is(:,37)        +is(:,38)        +is(:,39))         ./ 3;
es(:,53)         = (es(:,37)        +es(:,38)        +es(:,39))         ./ 3;
bulge(:,53)      = (bulge(:,37)     +bulge(:,38)     +bulge(:,39))      ./ 3;
norm_is(:,53)    = (norm_is(:,37)   +norm_is(:,38)   +norm_is(:,39))    ./ 3;
norm_es(:,53)    = (norm_es(:,37)   +norm_es(:,38)   +norm_es(:,39))    ./ 3;
norm_bulge(:,53) = (norm_bulge(:,37)+norm_bulge(:,38)+norm_bulge(:,39)) ./ 3;
prop_list(53)     = {'average Mgas out Vr>Vesc'};
unit_list(53)     = {'Msun/yr'};

% average Mgas out Vr>0, V>Vesc
is(:,54)         = (is(:,41)        +is(:,42)        +is(:,43))         ./ 3;
es(:,54)         = (es(:,41)        +es(:,42)        +es(:,43))         ./ 3;
bulge(:,54)      = (bulge(:,41)     +bulge(:,42)     +bulge(:,43))      ./ 3;
norm_is(:,54)    = (norm_is(:,41)   +norm_is(:,42)   +norm_is(:,43))    ./ 3;
norm_es(:,54)    = (norm_es(:,41)   +norm_es(:,42)   +norm_es(:,43))    ./ 3;
norm_bulge(:,54) = (norm_bulge(:,41)+norm_bulge(:,42)+norm_bulge(:,43)) ./ 3;
prop_list(54)     = {'average Mgas out Vr>0 V>Vesc'};
unit_list(54)     = {'Msun/yr'};

% spherical gas density
cm3 = 1./0.03363;
is(:,55)         = cm3.*is(:,4)    ./ ((4*pi/3).*((1000.*is(:,3)).^3));
es(:,55)         = cm3.*es(:,4)    ./ ((4*pi/3).*((1000.*es(:,3)).^3));
bulge(:,55)      = cm3.*bulge(:,4) ./ ((4*pi/3).*((1000.*bulge(:,3)).^3));
norm_is(:,55)    = norm_is(:,4)    ./ (norm_is(:,3).^3);
norm_es(:,55)    = norm_es(:,4)    ./ (norm_es(:,3).^3);
norm_bulge(:,55) = norm_bulge(:,4) ./ (norm_bulge(:,3).^3);
prop_list(55)     = {'Gas density'};
unit_list(55)     = {'1/cm3'};

% spherical star density
is(:,56)         = cm3.*is(:,5)    ./ ((4*pi/3).*((1000.*is(:,3)).^3));
es(:,56)         = cm3.*es(:,5)    ./ ((4*pi/3).*((1000.*es(:,3)).^3));
bulge(:,56)      = cm3.*bulge(:,5) ./ ((4*pi/3).*((1000.*bulge(:,3)).^3));
norm_is(:,56)    = norm_is(:,5)    ./ (norm_is(:,3).^3);
norm_es(:,56)    = norm_es(:,5)    ./ (norm_es(:,3).^3);
norm_bulge(:,56) = norm_bulge(:,5) ./ (norm_bulge(:,3).^3);
prop_list(56)     = {'Star density'};
unit_list(56)     = {'1/cm3'};

% spherical baryonic density
is(:,57)         = cm3.*is(:,6)    ./ ((4*pi/3).*((1000.*is(:,3)).^3));
es(:,57)         = cm3.*es(:,6)    ./ ((4*pi/3).*((1000.*es(:,3)).^3));
bulge(:,57)      = cm3.*bulge(:,6) ./ ((4*pi/3).*((1000.*bulge(:,3)).^3));
norm_is(:,57)    = norm_is(:,6)    ./ (norm_is(:,3).^3);
norm_es(:,57)    = norm_es(:,6)    ./ (norm_es(:,3).^3);
norm_bulge(:,57) = norm_bulge(:,6) ./ (norm_bulge(:,3).^3);
prop_list(57)     = {'Bar density'};
unit_list(57)     = {'1/cm3'};
clear cm3

% maximal mass, normalized distance, normalized height and lifetime
is(:,58:96)=0;
es(:,58:96)=0;
bulge(:,58:96)=0;
norm_is(:,58:96)=0;
norm_es(:,58:96)=0;
norm_bulge(:,58:96)=0;
avg_tff = zeros(1,2);

i=1;
while(i<length(is))
    if(is(i,1) == is(i+1,1) & is(i,2) > is(i+1,2))   % multiple snapshot clump
        j = i + 1;
        avg_tff(1) = is(i,6).*is(i,28);
        avg_tff(2) = is(i,6);
        while(is(i,1)==is(j,1) & is(i,2) > is(j,2) & j<length(is))
            j = j + 1;
            avg_tff(1) = avg_tff(1) + is(j,6).*is(j,28);
            avg_tff(2) = avg_tff(2) + is(j,6);
        end
        if(j==length(is) & is(i,1)==is(j,1) & is(i,2) > is(j,2))
            j = j + 1;
            avg_tff(1) = avg_tff(1) + is(j,6).*is(j,28);
            avg_tff(2) = avg_tff(2) + is(j,6);
        end
        avg_tff(1) = avg_tff(1) ./ avg_tff(2);
        is(i:(j-1),58) = max(is(i:(j-1),6));
        is(i:(j-1),59) = max(is(i:(j-1),15));
        is(i:(j-1),60) = max(is(i:(j-1),19));
        is(i:(j-1),61) = max(is(i:(j-1),20));
        is(i:(j-1),62) = is(j-1,48);
        norm_is(i:(j-1),58) = max(norm_is(i:(j-1),6));
        norm_is(i:(j-1),59) = max(norm_is(i:(j-1),15));
        norm_is(i:(j-1),60) = max(norm_is(i:(j-1),19));
        norm_is(i:(j-1),61) = max(norm_is(i:(j-1),20));
        norm_is(i:(j-1),62) = is(j-1,48)./avg_tff(1);
        i = j;
    else
        is(i,58) = is(i,6);
        is(i,59) = is(i,15);
        is(i,60) = is(i,19);
        is(i,61) = is(i,20);
        is(i,62) = is(i,48);
        norm_is(i,58) = norm_is(i,6);
        norm_is(i,59) = norm_is(i,15);
        norm_is(i,60) = norm_is(i,19);
        norm_is(i,61) = norm_is(i,20);
        norm_is(i,62) = norm_is(i,48);
        i = i + 1;
    end
end
if(i==length(is))
    is(length(is),58) = is(length(is),6);
    is(length(is),59) = is(length(is),15);
    is(length(is),60) = is(length(is),19);
    is(length(is),61) = is(length(is),20);
    is(length(is),62) = is(length(is),48);
    norm_is(length(is),58) = norm_is(length(is),6);
    norm_is(length(is),59) = norm_is(length(is),15);
    norm_is(length(is),60) = norm_is(length(is),19);
    norm_is(length(is),61) = norm_is(length(is),20);
    norm_is(length(is),62) = norm_is(length(is),48);
end
[max(norm_is(:,48)),min(norm_is(:,48)),length(find(isnan(norm_is(:,48))))]
[max(norm_is(:,62)),min(norm_is(:,62)),length(find(isnan(norm_is(:,62))))]
[max(is(:,48)),min(is(:,48)),length(find(isnan(is(:,48))))]
[max(is(:,62)),min(is(:,62)),length(find(isnan(is(:,62))))]

i=1;
while(i<length(es))
    if(es(i,1) == es(i+1,1) & es(i,2) > es(i+1,2))   % multiple snapshot clump
        j = i + 1;
        avg_tff(1) = es(i,6).*es(i,28);
        avg_tff(2) = es(i,6);
        while(es(i,1)==es(j,1) & es(i,2) > es(j,2) & j<length(es))
            j = j + 1;
            avg_tff(1) = avg_tff(1) + es(j,6).*es(j,28);
            avg_tff(2) = avg_tff(2) + es(j,6);
        end
        if(j==length(es) & es(i,1)==es(j,1) & es(i,2) > es(j,2))
            j = j + 1;
            avg_tff(1) = avg_tff(1) + es(j,6).*es(j,28);
            avg_tff(2) = avg_tff(2) + es(j,6);
        end
        avg_tff(1) = avg_tff(1) ./ avg_tff(2);
        es(i:(j-1),58) = max(es(i:(j-1),6));
        es(i:(j-1),59) = max(es(i:(j-1),15));
        es(i:(j-1),60) = max(es(i:(j-1),19));
        es(i:(j-1),61) = max(es(i:(j-1),20));
        es(i:(j-1),62) = es(j-1,48);
        norm_es(i:(j-1),58) = max(norm_es(i:(j-1),6));
        norm_es(i:(j-1),59) = max(norm_es(i:(j-1),15));
        norm_es(i:(j-1),60) = max(norm_es(i:(j-1),19));
        norm_es(i:(j-1),61) = max(norm_es(i:(j-1),20));
        norm_es(i:(j-1),62) = es(j-1,48)./avg_tff(1);
        i = j;
    else
        es(i,58) = es(i,6);
        es(i,59) = es(i,15);
        es(i,60) = es(i,19);
        es(i,61) = es(i,20);
        es(i,62) = es(i,48);
        norm_es(i,58) = norm_es(i,6);
        norm_es(i,59) = norm_es(i,15);
        norm_es(i,60) = norm_es(i,19);
        norm_es(i,61) = norm_es(i,20);
        norm_es(i,62) = norm_es(i,48);
        i = i + 1;
    end
end
if(i==length(es))
    es(length(es),58) = es(length(es),6);
    es(length(es),59) = es(length(es),15);
    es(length(es),60) = es(length(es),19);
    es(length(es),61) = es(length(es),20);
    es(length(es),62) = es(length(es),48);
    norm_es(length(es),58) = norm_es(length(es),6);
    norm_es(length(es),59) = norm_es(length(es),15);
    norm_es(length(es),60) = norm_es(length(es),19);
    norm_es(length(es),61) = norm_es(length(es),20);
    norm_es(length(es),62) = norm_es(length(es),48);
end

i=1;
while(i<length(bulge))
    if(bulge(i,1) == bulge(i+1,1) & bulge(i,2) > bulge(i+1,2))   % multiple snapshot clump
        j = i + 1;
        avg_tff(1) = bulge(i,6).*bulge(i,28);
        avg_tff(2) = bulge(i,6);
        while(bulge(i,1)==bulge(j,1) & bulge(i,2) > bulge(j,2) & j<length(bulge))
            j = j + 1;
            avg_tff(1) = avg_tff(1) + bulge(j,6).*bulge(j,28);
            avg_tff(2) = avg_tff(2) + bulge(j,6);
        end
        if(j==length(bulge) & bulge(i,1)==bulge(j,1) & bulge(i,2) > bulge(j,2))
            j = j + 1;
            avg_tff(1) = avg_tff(1) + bulge(j,6).*bulge(j,28);
            avg_tff(2) = avg_tff(2) + bulge(j,6);
        end
        avg_tff(1) = avg_tff(1) ./ avg_tff(2);
        bulge(i:(j-1),58) = max(bulge(i:(j-1),6));
        bulge(i:(j-1),59) = max(bulge(i:(j-1),15));
        bulge(i:(j-1),60) = max(bulge(i:(j-1),19));
        bulge(i:(j-1),61) = max(bulge(i:(j-1),20));
        bulge(i:(j-1),62) = bulge(j-1,48);
        norm_bulge(i:(j-1),58) = max(norm_bulge(i:(j-1),6));
        norm_bulge(i:(j-1),59) = max(norm_bulge(i:(j-1),15));
        norm_bulge(i:(j-1),60) = max(norm_bulge(i:(j-1),19));
        norm_bulge(i:(j-1),61) = max(norm_bulge(i:(j-1),20));
        norm_bulge(i:(j-1),62) = bulge(j-1,48)./avg_tff(1);
        i = j;
    else
        bulge(i,58) = bulge(i,6);
        bulge(i,59) = bulge(i,15);
        bulge(i,60) = bulge(i,19);
        bulge(i,61) = bulge(i,20);
        bulge(i,62) = bulge(i,48);
        norm_bulge(i,58) = norm_bulge(i,6);
        norm_bulge(i,59) = norm_bulge(i,15);
        norm_bulge(i,60) = norm_bulge(i,19);
        norm_bulge(i,61) = norm_bulge(i,20);
        norm_bulge(i,62) = norm_bulge(i,48);
        i = i + 1;
    end
end
if(i==length(bulge))
    bulge(length(bulge),58) = bulge(length(bulge),6);
    bulge(length(bulge),59) = bulge(length(bulge),15);
    bulge(length(bulge),60) = bulge(length(bulge),19);
    bulge(length(bulge),61) = bulge(length(bulge),20);
    bulge(length(bulge),62) = bulge(length(bulge),48);
    norm_bulge(length(bulge),58) = norm_bulge(length(bulge),6);
    norm_bulge(length(bulge),59) = norm_bulge(length(bulge),15);
    norm_bulge(length(bulge),60) = norm_bulge(length(bulge),19);
    norm_bulge(length(bulge),61) = norm_bulge(length(bulge),20);
    norm_bulge(length(bulge),62) = norm_bulge(length(bulge),48);
end

prop_list(58) = {'Max Mass'};
prop_list(59) = {'Max SFR'};
prop_list(60) = {'Max normalized distance'};
prop_list(61) = {'Max normalized height'};
prop_list(62) = {'Max time'};
unit_list(58) = {'Msun'};
unit_list(59) = {'Msun/yr'};
unit_list(60) = {'1'};
unit_list(61) = {'1'};
unit_list(62) = {'Myr'};

clear avg_tff
clear i j

% Add age of stars in first snapshot to lifetime
for i=1:(length(nis)-1)
    i
    z = is(nis(i)+1:nis(i+1),2);
    j=1;
    while(j<=(nis(i+1)-nis(i)))
        if(is(nis(i)+j,48)==0)
            z1 = is(nis(i)+j,2);
            if(z1==max(z))
                b = find(z<z1);
                z2 = max(z(b));
            else
                b = find(z>z1);
                z2 = min(z(b));
            end
            t1 = 0.95./( ((1+z1)/7)^1.5 );
            t2 = 0.95./( ((1+z2)/7)^1.5 );
            delt = 1000.*abs(t1-t2); %Myr
            if(is(nis(i)+j,12)>1e-6 & is(nis(i)+j,12)<delt)
                k = j+1;
                k1 = k;
                while(k<=(nis(i+1)-nis(i)))
                    if(is(nis(i)+k,1)==is(nis(i)+j,1) & is(nis(i)+k,2)<is(nis(i)+j,2))
                        k = k+1;
                        k1 = k;
                    else
                        k = 2*(nis(i+1)-nis(i)) + 10;
                    end
                end
                k = k1-1;
                if(k==j)
                    is(nis(i)+j,48) = is(nis(i)+j,48) + is(nis(i)+j,12);
                    is(nis(i)+j,62) = is(nis(i)+j,62) + is(nis(i)+j,12);
                    norm_is(nis(i)+j,48) = is(nis(i)+j,48)./is(nis(i)+j,28);
                    norm_is(nis(i)+j,62) = is(nis(i)+j,62)./is(nis(i)+j,28);
                else
                    temp = ( is(nis(i)+k,62)./norm_is(nis(i)+k,62) );
                    is(nis(i)+j:nis(i)+k,48) = is(nis(i)+j:nis(i)+k,48) + is(nis(i)+j,12);
                    is(nis(i)+j:nis(i)+k,62) = is(nis(i)+j:nis(i)+k,62) + is(nis(i)+j,12);
                    norm_is(nis(i)+j:nis(i)+k,48) = is(nis(i)+j:nis(i)+k,48)./is(nis(i)+j:nis(i)+k,28);
                    norm_is(nis(i)+j:nis(i)+k,62) = is(nis(i)+j:nis(i)+k,62)./temp;
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
clear i j k k1 z z1 z2 t1 t2 delt b temp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% t over global td
is(:,63)         = is(:,48)./is(:,30);
es(:,63)         = es(:,48)./es(:,30);
bulge(:,63)      = bulge(:,48)./bulge(:,30);
norm_is(:,63)    = is(:,62)./is(:,30);
norm_es(:,63)    = es(:,62)./es(:,30);
norm_bulge(:,63) = bulge(:,62)./bulge(:,30);
prop_list(63)     = {'t over global td'};
unit_list(63)     = {'1'};

% Vcirc
is(:,64)         = sqrt(4.3e-3.*is(:,6)./(1000.*is(:,3)));
es(:,64)         = sqrt(4.3e-3.*es(:,6)./(1000.*es(:,3)));
bulge(:,64)      = sqrt(4.3e-3.*bulge(:,6)./(1000.*bulge(:,3)));
norm_is(:,64)    = sqrt(4.3e-3.*is(:,6)./(1000.*is(:,3)));
norm_es(:,64)    = sqrt(4.3e-3.*es(:,6)./(1000.*es(:,3)));
norm_bulge(:,64) = sqrt(4.3e-3.*bulge(:,6)./(1000.*bulge(:,3)));
prop_list(64)     = {'Vcirc'};
unit_list(64)     = {'km/s'};

% eps_ff
is(:,65)         = is(:,15)./( (is(:,4)+is(:,15).*30.*1e6)./(1e6.*is(:,28)) );
es(:,65)         = es(:,15)./( (es(:,4)+es(:,15).*30.*1e6)./(1e6.*es(:,28)) );
bulge(:,65)      = bulge(:,15)./( (bulge(:,4)+bulge(:,15).*30.*1e6)./(1e6.*bulge(:,28)) );
norm_is(:,65)    = norm_is(:,15)./( norm_is(:,4)./(1./is(:,47)) );
norm_es(:,65)    = norm_es(:,15)./( norm_es(:,4)./(1./es(:,47)) );
norm_bulge(:,65) = norm_bulge(:,15)./( norm_bulge(:,4)./(1./bulge(:,47)) );
prop_list(65)     = {'EpsFF'};
unit_list(65)     = {'1'};

% mass_loading Rc Vr>0
is(:,66)         = is(:,34)./is(:,15);
es(:,66)         = es(:,34)./es(:,15);
bulge(:,66)      = bulge(:,34)./bulge(:,15);
norm_is(:,66)    = is(:,34)./is(:,15);
norm_es(:,66)    = es(:,34)./es(:,15);
norm_bulge(:,66) = bulge(:,34)./bulge(:,15);
prop_list(66)     = {'Mass Loading Rc Vr>0'};
unit_list(66)     = {'1'};

% mass_loading Rc Vr>Vesc
is(:,67)         = is(:,37)./is(:,15);
es(:,67)         = es(:,37)./es(:,15);
bulge(:,67)      = bulge(:,37)./bulge(:,15);
norm_is(:,67)    = is(:,37)./is(:,15);
norm_es(:,67)    = es(:,37)./es(:,15);
norm_bulge(:,67) = bulge(:,37)./bulge(:,15);
prop_list(67)     = {'Mass Loading Rc Vr>Vesc'};
unit_list(67)     = {'1'};

% mass_loading Rc Vr>0 V>Vesc
is(:,68)         = is(:,41)./is(:,15);
es(:,68)         = es(:,41)./es(:,15);
bulge(:,68)      = bulge(:,41)./bulge(:,15);
norm_is(:,68)    = is(:,41)./is(:,15);
norm_es(:,68)    = es(:,41)./es(:,15);
norm_bulge(:,68) = bulge(:,41)./bulge(:,15);
prop_list(68)     = {'Mass Loading Rc Vr>0 V>Vesc'};
unit_list(68)     = {'1'};

% mass_loading avg Vr>0
is(:,69)         = is(:,52)./is(:,15);
es(:,69)         = es(:,52)./es(:,15);
bulge(:,69)      = bulge(:,52)./bulge(:,15);
norm_is(:,69)    = is(:,52)./is(:,15);
norm_es(:,69)    = es(:,52)./es(:,15);
norm_bulge(:,69) = bulge(:,52)./bulge(:,15);
prop_list(69)     = {'Mass Loading avg Vr>0'};
unit_list(69)     = {'1'};

% mass_loading avg Vr>Vesc
is(:,70)         = is(:,53)./is(:,15);
es(:,70)         = es(:,53)./es(:,15);
bulge(:,70)      = bulge(:,53)./bulge(:,15);
norm_is(:,70)    = is(:,53)./is(:,15);
norm_es(:,70)    = es(:,53)./es(:,15);
norm_bulge(:,70) = bulge(:,53)./bulge(:,15);
prop_list(70)     = {'Mass Loading avg Vr>Vesc'};
unit_list(70)     = {'1'};

% mass_loading avg Vr>0 V>Vesc
is(:,71)         = is(:,54)./is(:,15);
es(:,71)         = es(:,54)./es(:,15);
bulge(:,71)      = bulge(:,54)./bulge(:,15);
norm_is(:,71)    = is(:,54)./is(:,15);
norm_es(:,71)    = es(:,54)./es(:,15);
norm_bulge(:,71) = bulge(:,54)./bulge(:,15);
prop_list(71)     = {'Mass Loading avg Vr>0 V>Vesc'};
unit_list(71)     = {'1'};

% mass_loading Frederic
is(:,72)         = is(:,40)./is(:,15);
es(:,72)         = es(:,40)./es(:,15);
bulge(:,72)      = bulge(:,40)./bulge(:,15);
norm_is(:,72)    = is(:,40)./is(:,15);
norm_es(:,72)    = es(:,40)./es(:,15);
norm_bulge(:,72) = bulge(:,40)./bulge(:,15);
prop_list(72)     = {'Mass Loading Frederic'};
unit_list(72)     = {'1'};

% eta_stars
is(:,73)         = is(:,45)./(is(:,15));
es(:,73)         = es(:,45)./(es(:,15));
bulge(:,73)      = bulge(:,45)./(bulge(:,15));
norm_is(:,73)    = is(:,45)./(is(:,15));
norm_es(:,73)    = es(:,45)./(es(:,15));
norm_bulge(:,73) = bulge(:,45)./(bulge(:,15));
prop_list(73)     = {'eta stars'};
unit_list(73)     = {'1'};

% eta_stars_net
b = find( is(:,15) > 1e-6 );
is(b,74)         = (is(b,45)-is(b,44))./(is(b,15));
b = find( es(:,15) > 1e-6 );
es(b,74)         = (es(b,45)-es(b,44))./(es(b,15));
b = find( bulge(:,15) > 1e-6 );
bulge(b,74)      = (bulge(b,45)-bulge(b,44))./(bulge(b,15));
norm_is(:,74)    = is(:,74);
norm_es(:,74)    = es(:,74);
norm_bulge(:,74) = bulge(:,74);
prop_list(74)     = {'eta stars net'};
unit_list(74)     = {'1'};
clear b

% Inflow gas mass conservation: Rc Vr>Vesc, Rc Vr>0 V>Vesc, avg Vr>Vesc, avg Vr>0 V>Vesc, Frederic
vec = 1:(length(is)-1);
b = find(is(vec,1)==is(vec+1,1) & is(vec,2)>is(vec+1,2) & is(vec,48)<is(vec+1,48));
delt = (is(vec(b+1),48) - is(vec(b),48)) .* 1e6;
delm = is(vec(b+1),4) - is(vec(b),4);
mu = 0.8;
is(vec(b),75) = max( 1e-6, delm./delt + mu.*is(vec(b),15) + is(vec(b),37) );
is(vec(b),76) = max( 1e-6, delm./delt + mu.*is(vec(b),15) + is(vec(b),41) );
is(vec(b),77) = max( 1e-6, delm./delt + mu.*is(vec(b),15) + is(vec(b),53) );
is(vec(b),78) = max( 1e-6, delm./delt + mu.*is(vec(b),15) + is(vec(b),54) );
is(vec(b),79) = max( 1e-6, delm./delt + mu.*is(vec(b),15) + is(vec(b),40) );
clear vec b delt delm mu

vec = 1:(length(es)-1);
b = find(es(vec,1)==es(vec+1,1) & es(vec,2)>es(vec+1,2) & es(vec,48)<es(vec+1,48));
delt = (es(vec(b+1),48) - es(vec(b),48)) .* 1e6;
delm = es(vec(b+1),4) - es(vec(b),4);
mu = 0.8;
es(vec(b),75) = max( 1e-6, delm./delt + mu.*es(vec(b),15) + es(vec(b),37) );
es(vec(b),76) = max( 1e-6, delm./delt + mu.*es(vec(b),15) + es(vec(b),41) );
es(vec(b),77) = max( 1e-6, delm./delt + mu.*es(vec(b),15) + es(vec(b),53) );
es(vec(b),78) = max( 1e-6, delm./delt + mu.*es(vec(b),15) + es(vec(b),54) );
es(vec(b),79) = max( 1e-6, delm./delt + mu.*es(vec(b),15) + es(vec(b),40) );
clear vec b delt delm mu

vec = 1:(length(bulge)-1);
b = find(bulge(vec,1)==bulge(vec+1,1) & bulge(vec,2)>bulge(vec+1,2) & bulge(vec,48)<bulge(vec+1,48));
delt = (bulge(vec(b+1),48) - bulge(vec(b),48)) .* 1e6;
delm = bulge(vec(b+1),4) - bulge(vec(b),4);
mu = 0.8;
bulge(vec(b),75) = max( 1e-6, delm./delt + mu.*bulge(vec(b),15) + bulge(vec(b),37) );
bulge(vec(b),76) = max( 1e-6, delm./delt + mu.*bulge(vec(b),15) + bulge(vec(b),41) );
bulge(vec(b),77) = max( 1e-6, delm./delt + mu.*bulge(vec(b),15) + bulge(vec(b),53) );
bulge(vec(b),78) = max( 1e-6, delm./delt + mu.*bulge(vec(b),15) + bulge(vec(b),54) );
bulge(vec(b),79) = max( 1e-6, delm./delt + mu.*bulge(vec(b),15) + bulge(vec(b),40) );
clear vec b delt delm mu

norm_is(:,75:79)    = is(:,75:79);
norm_es(:,75:79)    = es(:,75:79);
norm_bulge(:,75:79) = bulge(:,75:79);

prop_list(75)     = {'Gas inflow gas mass conservation Rc Vr>Vesc'};
prop_list(76)     = {'Gas inflow gas mass conservation Rc Vr>0 V>Vesc'};
prop_list(77)     = {'Gas inflow gas mass conservation avg Vr>Vesc'};
prop_list(78)     = {'Gas inflow gas mass conservation avg Vr>0 V>Vesc'};
prop_list(79)     = {'Gas inflow gas mass conservation Frederic'};
unit_list(75)     = {'Msun/yr'};
unit_list(76)     = {'Msun/yr'};
unit_list(77)     = {'Msun/yr'};
unit_list(78)     = {'Msun/yr'};
unit_list(79)     = {'Msun/yr'};

% alpha_accretion
is(:,80)         = max( 1e-2, 2e6.*is(:,51).*is(:,30)./(is(:,6)));
es(:,80)         = max( 1e-2, 2e6.*es(:,51).*es(:,30)./(es(:,6)));
bulge(:,80)      = max( 1e-2, 2e6.*bulge(:,51).*bulge(:,30)./(bulge(:,6)));
prop_list(80)     = {'alpha accretion avg Min Vr<0'};
unit_list(80)     = {'1'};

is(:,81)         = max( 1e-2, 2e6.*is(:,75).*is(:,30)./(is(:,6)));
es(:,81)         = max( 1e-2, 2e6.*es(:,75).*es(:,30)./(es(:,6)));
bulge(:,81)      = max( 1e-2, 2e6.*bulge(:,75).*bulge(:,30)./(bulge(:,6)));
prop_list(81)     = {'alpha accretion gas mass conserv Rc Vr>Vesc'};
unit_list(81)     = {'1'};

is(:,82)         = max( 1e-2, 2e6.*is(:,76).*is(:,30)./(is(:,6)));
es(:,82)         = max( 1e-2, 2e6.*es(:,76).*es(:,30)./(es(:,6)));
bulge(:,82)      = max( 1e-2, 2e6.*bulge(:,76).*bulge(:,30)./(bulge(:,6)));
prop_list(82)     = {'alpha accretion gas mass conserv Rc Vr>0 V>Vesc'};
unit_list(82)     = {'1'};

is(:,83)         = max( 1e-2, 2e6.*is(:,77).*is(:,30)./(is(:,6)));
es(:,83)         = max( 1e-2, 2e6.*es(:,77).*es(:,30)./(es(:,6)));
bulge(:,83)      = max( 1e-2, 2e6.*bulge(:,77).*bulge(:,30)./(bulge(:,6)));
prop_list(83)     = {'alpha accretion gas mass conserv avg Vr>Vesc'};
unit_list(83)     = {'1'};

is(:,84)         = max( 1e-2, 2e6.*is(:,78).*is(:,30)./(is(:,6)));
es(:,84)         = max( 1e-2, 2e6.*es(:,78).*es(:,30)./(es(:,6)));
bulge(:,84)      = max( 1e-2, 2e6.*bulge(:,78).*bulge(:,30)./(bulge(:,6)));
prop_list(84)     = {'alpha accretion gas mass conserv avg Vr>0 V>Vesc'};
unit_list(84)     = {'1'};

is(:,85)         = max( 1e-2, 2e6.*is(:,79).*is(:,30)./(is(:,6)));
es(:,85)         = max( 1e-2, 2e6.*es(:,79).*es(:,30)./(es(:,6)));
bulge(:,85)      = max( 1e-2, 2e6.*bulge(:,79).*bulge(:,30)./(bulge(:,6)));
prop_list(85)     = {'alpha accretion gas mass conserv Frederic'};
unit_list(85)     = {'1'};

norm_is(:,80:85)    = is(:,80:85);
norm_es(:,80:85)    = es(:,80:85);
norm_bulge(:,80:85) = bulge(:,80:85);

% Inflow mass conservation: Rc Vr>Vesc, Rc Vr>0 V>Vesc, avg Vr>Vesc, avg Vr>0 V>Vesc, Frederic
vec = 1:(length(is)-1);
b = find(is(vec,1)==is(vec+1,1) & is(vec,2)>is(vec+1,2) & is(vec,48)<is(vec+1,48));
delt = (is(vec(b+1),48) - is(vec(b),48)) .* 1e6;
delm = is(vec(b+1),6) - is(vec(b),6);
is(vec(b),86) = max(1e-6, delm./delt + is(vec(b),37) + is(vec(b+1),45) - is(vec(b+1),44));
is(vec(b),87) = max(1e-6, delm./delt + is(vec(b),41) + is(vec(b+1),45) - is(vec(b+1),44));
is(vec(b),88) = max(1e-6, delm./delt + is(vec(b),53) + is(vec(b+1),45) - is(vec(b+1),44));
is(vec(b),89) = max(1e-6, delm./delt + is(vec(b),54) + is(vec(b+1),45) - is(vec(b+1),44));
is(vec(b),90) = max(1e-6, delm./delt + is(vec(b),40) + is(vec(b+1),45) - is(vec(b+1),44));
clear vec b delt delm

vec = 1:(length(es)-1);
b = find(es(vec,1)==es(vec+1,1) & es(vec,2)>es(vec+1,2) & es(vec,48)<es(vec+1,48));
delt = (es(vec(b+1),48) - es(vec(b),48)) .* 1e6;
delm = es(vec(b+1),6) - es(vec(b),6);
es(vec(b),86) = max(1e-6, delm./delt + es(vec(b),37) + es(vec(b+1),45) - es(vec(b+1),44));
es(vec(b),87) = max(1e-6, delm./delt + es(vec(b),41) + es(vec(b+1),45) - es(vec(b+1),44));
es(vec(b),88) = max(1e-6, delm./delt + es(vec(b),53) + es(vec(b+1),45) - es(vec(b+1),44));
es(vec(b),89) = max(1e-6, delm./delt + es(vec(b),54) + es(vec(b+1),45) - es(vec(b+1),44));
es(vec(b),90) = max(1e-6, delm./delt + es(vec(b),40) + es(vec(b+1),45) - es(vec(b+1),44));
clear vec b delt delm

vec = 1:(length(bulge)-1);
b = find(bulge(vec,1)==bulge(vec+1,1) & bulge(vec,2)>bulge(vec+1,2) & bulge(vec,48)<bulge(vec+1,48));
delt = (bulge(vec(b+1),48) - bulge(vec(b),48)) .* 1e6;
delm = bulge(vec(b+1),6) - bulge(vec(b),6);
bulge(vec(b),86) = max(1e-6, delm./delt + bulge(vec(b),37) + bulge(vec(b+1),45) - bulge(vec(b+1),44));
bulge(vec(b),87) = max(1e-6, delm./delt + bulge(vec(b),41) + bulge(vec(b+1),45) - bulge(vec(b+1),44));
bulge(vec(b),88) = max(1e-6, delm./delt + bulge(vec(b),53) + bulge(vec(b+1),45) - bulge(vec(b+1),44));
bulge(vec(b),89) = max(1e-6, delm./delt + bulge(vec(b),54) + bulge(vec(b+1),45) - bulge(vec(b+1),44));
bulge(vec(b),90) = max(1e-6, delm./delt + bulge(vec(b),40) + bulge(vec(b+1),45) - bulge(vec(b+1),44));
clear vec b delt delm

norm_is(:,86:90)    = is(:,86:90);
norm_es(:,86:90)    = es(:,86:90);
norm_bulge(:,86:90) = bulge(:,86:90);

prop_list(86)     = {'Gas inflow baryonic mass conservation Rc Vr>Vesc'};
prop_list(87)     = {'Gas inflow baryonic mass conservation Rc Vr>0 V>Vesc'};
prop_list(88)     = {'Gas inflow baryonic mass conservation avg Vr>Vesc'};
prop_list(89)     = {'Gas inflow baryonic mass conservation avg Vr>0 V>Vesc'};
prop_list(90)     = {'Gas inflow baryonic mass conservation Frederic'};
unit_list(86)     = {'Msun/yr'};
unit_list(87)     = {'Msun/yr'};
unit_list(88)     = {'Msun/yr'};
unit_list(89)     = {'Msun/yr'};
unit_list(90)     = {'Msun/yr'};

% alpha_accretion
is(:,91)         = max( 1e-2, 2e6.*is(:,51).*is(:,30)./(is(:,6)));
es(:,91)         = max( 1e-2, 2e6.*es(:,51).*es(:,30)./(es(:,6)));
bulge(:,91)      = max( 1e-2, 2e6.*bulge(:,51).*bulge(:,30)./(bulge(:,6)));
prop_list(91)     = {'alpha accretion avg Min Vr<0'};
unit_list(91)     = {'1'};

is(:,92)         = max( 1e-2, 2e6.*is(:,86).*is(:,30)./(is(:,6)));
es(:,92)         = max( 1e-2, 2e6.*es(:,86).*es(:,30)./(es(:,6)));
bulge(:,92)      = max( 1e-2, 2e6.*bulge(:,86).*bulge(:,30)./(bulge(:,6)));
prop_list(92)     = {'alpha accretion baryonic mass conserv Rc Vr>Vesc'};
unit_list(92)     = {'1'};

is(:,93)         = max( 1e-2, 2e6.*is(:,87).*is(:,30)./(is(:,6)));
es(:,93)         = max( 1e-2, 2e6.*es(:,87).*es(:,30)./(es(:,6)));
bulge(:,93)      = max( 1e-2, 2e6.*bulge(:,87).*bulge(:,30)./(bulge(:,6)));
prop_list(93)     = {'alpha accretion baryonic mass conserv Rc Vr>0 V>Vesc'};
unit_list(93)     = {'1'};

is(:,94)         = max( 1e-2, 2e6.*is(:,88).*is(:,30)./(is(:,6)));
es(:,94)         = max( 1e-2, 2e6.*es(:,88).*es(:,30)./(es(:,6)));
bulge(:,94)      = max( 1e-2, 2e6.*bulge(:,88).*bulge(:,30)./(bulge(:,6)));
prop_list(94)     = {'alpha accretion baryonic mass conserv avg Vr>Vesc'};
unit_list(94)     = {'1'};

is(:,95)         = max( 1e-2, 2e6.*is(:,89).*is(:,30)./(is(:,6)));
es(:,95)         = max( 1e-2, 2e6.*es(:,89).*es(:,30)./(es(:,6)));
bulge(:,95)      = max( 1e-2, 2e6.*bulge(:,89).*bulge(:,30)./(bulge(:,6)));
prop_list(95)     = {'alpha accretion baryonic mass conserv avg Vr>Vesc'};
unit_list(95)     = {'1'};

is(:,96)         = max( 1e-2, 2e6.*is(:,90).*is(:,30)./(is(:,6)));
es(:,96)         = max( 1e-2, 2e6.*es(:,90).*es(:,30)./(es(:,6)));
bulge(:,96)      = max( 1e-2, 2e6.*bulge(:,90).*bulge(:,30)./(bulge(:,6)));
prop_list(96)     = {'alpha accretion baryonic mass conserv Frederic'};
unit_list(96)     = {'1'};

norm_is(:,91:96)    = is(:,91:96);
norm_es(:,91:96)    = es(:,91:96);
norm_bulge(:,91:96) = bulge(:,91:96);


