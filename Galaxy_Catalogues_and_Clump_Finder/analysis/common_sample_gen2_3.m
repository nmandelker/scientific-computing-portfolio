function [is2_common, es2_common, bulge2_common, disc2_common, nis2_common, nes2_common, ...
    nbulge2_common, ndisc2_common, is3_common, es3_common, bulge3_common, ...
    disc3_common, nis3_common, nes3_common, nbulge3_common, ndisc3_common, ...
    ind_common, amax_common, amin_common] = ...
common_sample_gen2_3(is2, es2, bulge2, disc2, nis2, nes2, nbulge2, ndisc2, ...
    is3, es3, bulge3, disc3, nis3, nes3, nbulge3, ndisc3)

ind2 = [1:3, 5:16, 19, 21, 23, 25:35];
ind3 = [1:17, 19:35];
ind_common = ind2;
%ind_common = [1:3, 5:16, 19, 21, 23, 25:34];
N_common = length(ind_common);
amin_common = zeros(N_common,1);
amax_common = zeros(N_common,1);

is2_common = [];
es2_common = [];
bulge2_common = [];
disc2_common = [];
is3_common = [];
es3_common = [];
bulge3_common = [];
disc3_common = [];
nis2_common = zeros(1,length(ind2));
nes2_common = zeros(1,length(ind2));
nbulge2_common = zeros(1,length(ind2));
ndisc2_common = zeros(1,length(ind2));
nis3_common = zeros(1,length(ind3));
nes3_common = zeros(1,length(ind3));
nbulge3_common = zeros(1,length(ind3));
ndisc3_common = zeros(1,length(ind3));

for i=1:N_common
    j = ind_common(i);
    b2 = find(ind2==j);
    zmax2 = disc2(ndisc2(b2)+1,1);
    zmin2 = disc2(ndisc2(b2+1),1);
    amin2 = 1./(1+zmax2);
    amax2 = 1./(1+zmin2);
    b3 = find(ind3==j);
    zmax3 = disc3(ndisc3(b3)+1,1);
    zmin3 = disc3(ndisc3(b3+1),1);
    amin3 = 1./(1+zmax3);
    amax3 = 1./(1+zmin3);
    amin_common(i) = max(amin2,amin3);
    amax_common(i) = min(amax2,amax3);
    
    % Now fix amin/amax for bugs
    if(j==1)    % VELA 01
        amin_common(i) = max(0.42,amin_common(i));   %Before this, gen2 and gen3 trace different progenitors of a merger
    end
    if(j==3)    % VELA 03
        amin_common(i) = max(0.18,amin_common(i));   %I don't understand why, but my clump catalogues are missin a=0.17 for gen3
    end
    if(j==9)    % VELA 09
        amax_common(i) = min(0.40,amax_common(i));   %Bug after this. Rvir decreases and gas disc dissapears
    end
    if(j==11)    % VELA 11
        amax_common(i) = min(0.46,amax_common(i));   %Bug after this
    end
    if(j==12)    % VELA 12
        amax_common(i) = min(0.44,amax_common(i));   %Bug after this
    end
    if(j==14)    % VELA 14
        amax_common(i) = min(0.41,amax_common(i));   %Bug after this
    end
    if(j==34)    % VELA 34
        amin_common(i) = max(0.20,amin_common(i));   %Before this, gen2 and gen3 trace different progenitors of a merger
    end
    if(j==35)    % VELA 35
        amin_common(i) = max(0.12,amin_common(i));   %I don't understand why, but my clump catalogues start from a=0.12
    end
    
    adisc3  = 1./(1+disc3(ndisc3(b3)+1:ndisc3(b3+1),1));
    ais3    = 1./(1+is3(nis3(b3)+1:nis3(b3+1),2));
    aes3    = 1./(1+es3(nes3(b3)+1:nes3(b3+1),2));
    abulge3 = 1./(1+bulge3(nbulge3(b3)+1:nbulge3(b3+1),2));
    ais2    = 1./(1+is2(nis2(b2)+1:nis2(b2+1),2));
    aes2    = 1./(1+es2(nes2(b2)+1:nes2(b2+1),2));
    abulge2 = 1./(1+bulge2(nbulge2(b2)+1:nbulge2(b2+1),2));
    for k=(ndisc2(b2)+1):ndisc2(b2+1)
        adisc2 = 1/(1+disc2(k,1));
        if(adisc2<(amax_common(i)+0.005) & adisc2>(amin_common(i)-0.005))
            if(min(abs(adisc3(:)-adisc2))<0.005)
                ndisc2_common(b2) = ndisc2_common(b2) + 1;
                disc2_common = [disc2_common', k]';
                ndisc3_common(b3) = ndisc3_common(b3) + 1;
                bcommon = find(abs(adisc3(:)-adisc2)<0.005);
                disc3_common = [disc3_common', ndisc3(b3)+bcommon']';
                                
                bcommon = find(abs(ais2(:)-adisc2)<0.005);
                nis2_common(b2) = nis2_common(b2) + length(bcommon);
                is2_common = sort([is2_common', nis2(b2)+sort(bcommon)'])';
                bcommon = find(abs(aes2(:)-adisc2)<0.005);
                nes2_common(b2) = nes2_common(b2) + length(bcommon);
                es2_common = sort([es2_common', nes2(b2)+sort(bcommon)'])';
                bcommon = find(abs(abulge2(:)-adisc2)<0.005);
                nbulge2_common(b2) = nbulge2_common(b2) + length(bcommon);
                bulge2_common = sort([bulge2_common', nbulge2(b2)+sort(bcommon)'])';
                
                bcommon = find(abs(ais3(:)-adisc2)<0.005);
                nis3_common(b3) = nis3_common(b3) + length(bcommon);
                is3_common = sort([is3_common', nis3(b3)+sort(bcommon)'])';
                bcommon = find(abs(aes3(:)-adisc2)<0.005);
                nes3_common(b3) = nes3_common(b3) + length(bcommon);
                es3_common = sort([es3_common', nes3(b3)+sort(bcommon)'])';
                bcommon = find(abs(abulge3(:)-adisc2)<0.005);
                nbulge3_common(b3) = nbulge3_common(b3) + length(bcommon);
                bulge3_common = sort([bulge3_common', nbulge3(b3)+sort(bcommon)'])';
            end
        end
    end
end
amax_common = amax_common';
amin_common = amin_common';