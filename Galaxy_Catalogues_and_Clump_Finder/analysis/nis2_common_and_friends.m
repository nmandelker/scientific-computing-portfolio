nis2_common_net = nis2_common;
nes2_common_net = nes2_common;
nbulge2_common_net = nbulge2_common;
ndisc2_common_net = ndisc2_common;
N_cat2_common_net = N_cat2_common;

nis3_common_net = nis3_common;
nes3_common_net = nes3_common;
nbulge3_common_net = nbulge3_common;
ndisc3_common_net = ndisc3_common;
N_cat3_common_net = N_cat3_common;

nis2_common = zeros(1,30);
nes2_common = zeros(1,30);
nbulge2_common = zeros(1,30);
ndisc2_common = zeros(1,30);
N_cat2_common = zeros(1,30);

nis3_common = zeros(1,35);
nes3_common = zeros(1,35);
nbulge3_common = zeros(1,35);
ndisc3_common = zeros(1,35);
N_cat3_common = zeros(1,35);

for i=2:30
    nis2_common(i) = sum(nis2_common_net(1:i-1));
    nes2_common(i) = sum(nes2_common_net(1:i-1));
    nbulge2_common(i) = sum(nbulge2_common_net(1:i-1));
    ndisc2_common(i) = sum(ndisc2_common_net(1:i-1));
    N_cat2_common(i) = sum(N_cat2_common_net(1:i-1));
end
for i=2:35
    nis3_common(i) = sum(nis3_common_net(1:i-1));
    nes3_common(i) = sum(nes3_common_net(1:i-1));
    nbulge3_common(i) = sum(nbulge3_common_net(1:i-1));
    ndisc3_common(i) = sum(ndisc3_common_net(1:i-1));
    N_cat3_common(i) = sum(N_cat3_common_net(1:i-1));
end
clear i