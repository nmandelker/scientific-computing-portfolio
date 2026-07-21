function [nspec, mu, cool, heat, cool_tot, heat_tot] = cooling_equilibrium(T, nH, z, zmet, Madau, J0_multiply)

aexp = 1/(1+z);
if(aexp < 1/9)
    J0simple = 0;
elseif(aexp < 1/4)
    J0simple = 4.*aexp;
elseif(aexp < 1/3)
    J0simple = 1;
else
    J0simple = 1/((3*aexp)^3);
end
%Parameters
aspec = 0;
J0_Theuns = 1e-21   .* J0_multiply .* J0simple;
normfacJ0 = 0.74627 .* J0_multiply;

X = 0.76;
Y = 1-X;
nHe = nH.*Y./(4.*X);

dumfac_ion = 2;
dumfac_rec = 0.75;
T3 = T./1e3;
T5 = T./1e5;
T6 = T./1e6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RATES
%%%%%%%

if(Madau == 1)
    J0_min = 1e-25;
    %Photo-Ionization rates (Madau)
    taux_rad(1) = max( normfacJ0.*exp( -31.04 + 2.795.*z - 0.5589.*(z.^2) ), ...
        1.26e10.*J0_min./(3.0+aspec) );
    taux_rad(2) = max( normfacJ0.*exp( -31.08 + 2.822.*z - 0.5664.*(z.^2) ), ...
        1.48e10.*J0_min.*(0.553.^aspec).*( 1.66./(2.05+aspec) - 0.66./(3.05+aspec) ) );
    taux_rad(3) = max( normfacJ0.*exp( -34.30 + 1.826.*z - 0.3899.*(z.^2) ), ...
        3.34e9.* J0_min.*(0.249.^aspec)./(3.0+aspec) );
    
    %Radiative heating rates (Madau)
    h_rad(1) = max( normfacJ0.*exp( -56.62 + 2.788.*z - 0.5594.*(z.^2) ), ...
        ( 2.91e-1.*J0_min./(2.0+aspec) ) ./ ( 3.0+aspec ) );
    h_rad(2) = max( normfacJ0.*exp( -56.06 + 2.800.*z - 0.5532.*(z.^2) ), ...
        ( 5.84e-1.*J0_min.*(0.553.^aspec) ) .* ( 1.66./(1.05+aspec) - 2.32./(2.05+aspec) + 0.66./(3.05+aspec) ) );
    h_rad(3) = max( normfacJ0.*exp( -58.67 + 1.888.*z - 0.3947.*(z.^2) ), ...
        ( 2.92e-1.*J0_min.*(0.249.^aspec)./(2.0+aspec) ) ./ (3.0+aspec) );
else
    %Photo-Ionization rates (Theuns)  
    taux_rad(1) = 1.26e10.*J0_Theuns./(3.0+aspec);
    taux_rad(2) = 1.48e10.*J0_Theuns.*(0.553.^aspec).*( 1.66./(2.05+aspec) - 0.66./(3.05+aspec) );
    taux_rad(3) = 3.34e9.* J0_Theuns.*(0.249.^aspec)./(3.0+aspec);
    
    %Radiative heating rates (Theuns)
    h_rad(1) = ( 2.91e-1.*J0_Theuns./(2.0+aspec) ) ./ ( 3.0+aspec );
    h_rad(2) = ( 5.84e-1.*J0_Theuns.*(0.553.^aspec) ) .* ...
        ( 1.66./(1.05+aspec) - 2.32./(2.05+aspec) + 0.66./(3.05+aspec) );
    h_rad(3) = ( 2.92e-1.*J0_Theuns.*(0.249.^aspec)./(2.0+aspec) ) ./ (3.0+aspec);
end

%Collisional Ionization rates
taux_ion(1) = dumfac_ion.*(5.85e-11.*sqrt(T)./(1.+sqrt(T5))).*exp(-157809.1./T);
taux_ion(2) = dumfac_ion.*(2.38e-11.*sqrt(T)./(1.+sqrt(T5))).*exp(-285335.4./T);
taux_ion(3) = dumfac_ion.*(5.68e-12.*sqrt(T)./(1.+sqrt(T5))).*exp(-631515.0./T);

%Dielectric recombination
taux_die = 1.9e-3.*(T.^(-1.5)).*exp(-470000./T).*(1+0.3.*exp(-94000./T));

%Recombination rates
taux_rec(1) = dumfac_rec.*8.40e-11./(sqrt(T).*(T3^0.2).*(1.+T6.^0.7));
taux_rec(2) = 1.50e-10./(T.^0.6353) + taux_die;
taux_rec(3) = 3.36e-10./(sqrt(T).*(T3.^0.2).*(1.+T6.^0.7));

%Bremsstrahlung cooling rate
cool_bre(1) = 1.42e-27 .* sqrt(T) .* (1.1 + 0.34.*exp( -((5.5-log10(T)).^2)./3 ));
cool_bre(2) = 1.42e-27 .* sqrt(T) .* (1.1 + 0.34.*exp( -((5.5-log10(T)).^2)./3 ));
cool_bre(3) = 5.68e-27 .* sqrt(T) .* (1.1 + 0.34.*exp( -((5.5-log10(T)).^2)./3 ));

%Line cooling rate
cool_exc(1) = ( 7.50e-19 ./  (1.+sqrt(T5))               ) .* exp(-118348./T);
cool_exc(2) = ( 9.10e-27 ./ ((1.+sqrt(T5)).*(T.^0.1687)) ) .* exp(-13179./T);
cool_exc(3) = ( 5.54e-17 ./ ((1.+sqrt(T5)).*(T.^0.3970)) ) .* exp(-473638./T);

%Recombination cooling rate
cool_rec(1) = 8.70e-27 .* sqrt(T) ./ ( (T3.^0.2) .* (1.+T6.^0.7) );
cool_rec(2) = 1.55e-26 .* T.^0.3647;
cool_rec (3)= 3.48e-26 .* sqrt(T) ./ ( (T3.^0.2) .* (1.+T6.^0.7) );

%Dielectric cooling
cool_die = 1.24e-13 .* (T.^(-1.5)) .* exp(-470000./T) .* (1 + 0.3.*exp(-94000./T));

%Ion cooling
cool_ion(1) = dumfac_ion .* 1.27e-21 .* sqrt(T) .* exp(-157809.1 ./ T) ./ (1 + sqrt(T5));
cool_ion(2) = dumfac_ion .* 9.38e-22 .* sqrt(T) .* exp(-285335.4 ./ T) ./ (1 + sqrt(T5));
cool_ion(3) = dumfac_ion .* 4.95e-22 .* sqrt(T) .* exp(-631515.0 ./ T) ./ (1 + sqrt(T5));

%Compton cooling
cool_com = 5.406e-36 .* T .* ((1+z).^4);

%Compton heating
heat_com = 5.406e-36 .* 2.726 .* ((1+z).^5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Chemical equilibrium
ne = nH;
err_ne = 1;
while(err_ne > 1e-8)
    t_ion2(1) = taux_ion(1) + taux_rad(1)./max(ne, 1e-15.*nH);
    t_ion2(2) = taux_ion(2) + taux_rad(2)./max(ne, 1e-15.*nH);
    t_ion2(3) = taux_ion(3) + taux_rad(3)./max(ne, 1e-15.*nH);
    
    nHI  = taux_rec(1)./(t_ion2(1) + taux_rec(1)).*nH;
    nHII = t_ion2(1)  ./(t_ion2(1) + taux_rec(1)).*nH;
    
    x1 = (taux_rec(3).*taux_rec(2) + t_ion2(2).*taux_rec(3) + t_ion2(3).*t_ion2(2));
    nHeIII = (t_ion2(3)  .*t_ion2(2)  ./x1).*nHe;
    nHeII  = (t_ion2(2)  .*taux_rec(3)./x1).*nHe;
    nHeI   = (taux_rec(3).*taux_rec(2)./x1).*nHe;
    
    err_ne = abs( (ne - (nHII + nHeII + 2.*nHeIII)) ./ nH );
    ne = 0.5.*ne + 0.5.*(nHII + nHeII + 2.*nHeIII);
end
ntot     = ne + nHI + nHII + nHeI + nHeII + nHeIII;
mu       = nH./(X.*ntot);
nspec(1)=ne;
nspec(2)=nHI;
nspec(3)=nHII;
nspec(4)=nHeI;
nspec(5)=nHeII;
nspec(6)=nHeIII;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Metal cooling
temperature_cc07 = [3.9684, 4.0187, 4.0690, 4.1194, 4.1697, 4.2200, 4.2703, ...
                    4.3206, 4.3709, 4.4212, 4.4716, 4.5219, 4.5722, 4.6225, ...
                    4.6728, 4.7231, 4.7734, 4.8238, 4.8741, 4.9244, 4.9747, ...
                    5.0250, 5.0753, 5.1256, 5.1760, 5.2263, 5.2766, 5.3269, ...
                    5.3772, 5.4275, 5.4778, 5.5282, 5.5785, 5.6288, 5.6791, ...
                    5.7294, 5.7797, 5.8300, 5.8804, 5.9307, 5.9810, 6.0313, ...
                    6.0816, 6.1319, 6.1822, 6.2326, 6.2829, 6.3332, 6.3835, ...
                    6.4338, 6.4841, 6.5345, 6.5848, 6.6351, 6.6854, 6.7357, ...
                    6.7860, 6.8363, 6.8867, 6.9370, 6.9873, 7.0376, 7.0879, ...
                    7.1382, 7.1885, 7.2388, 7.2892, 7.3395, 7.3898, 7.4401, ...
                    7.4904, 7.5407, 7.5911, 7.6414, 7.6917, 7.7420, 7.7923, ...
                    7.8426, 7.8929, 7.9433, 7.9936, 8.0439, 8.0942, 8.1445, ...
                    8.1948, 8.2451, 8.2955, 8.3458, 8.3961, 8.4464, 8.4967];

excess_cooling_cc07 = [-24.9949, -24.7270, -24.0473, -23.0713, -22.2907, -21.8917, -21.8058, ...
                       -21.8501, -21.9142, -21.9553, -21.9644, -21.9491, -21.9134, -21.8559, ...
                       -21.7797, -21.6863, -21.5791, -21.4648, -21.3640, -21.2995, -21.2691, ...
                       -21.2658, -21.2838, -21.2985, -21.2941, -21.2845, -21.2809, -21.2748, ...
                       -21.2727, -21.3198, -21.4505, -21.5921, -21.6724, -21.6963, -21.6925, ...
                       -21.6892, -21.7142, -21.7595, -21.7779, -21.7674, -21.7541, -21.7532, ...
                       -21.7679, -21.7866, -21.8052, -21.8291, -21.8716, -21.9316, -22.0055, ...
                       -22.0800, -22.1600, -22.2375, -22.3126, -22.3701, -22.4125, -22.4353, ...
                       -22.4462, -22.4450, -22.4406, -22.4337, -22.4310, -22.4300, -22.4356, ...
                       -22.4455, -22.4631, -22.4856, -22.5147, -22.5444, -22.5718, -22.5904, ...
                       -22.6004, -22.5979, -22.5885, -22.5728, -22.5554, -22.5350, -22.5159, ...
                       -22.4955, -22.4781, -22.4600, -22.4452, -22.4262, -22.4089, -22.3900, ...
                       -22.3722, -22.3529, -22.3339, -22.3137, -22.2936, -22.2729, -22.2521 ];

excess_prime_cc07 = [ 2.0037,  4.7267, 12.2283, 13.5820,  9.8755,  4.8379,  1.8046, ...
                      1.4574,  1.8086,  2.0685,  2.2012,  2.2250,  2.2060,  2.1605, ...
                      2.1121,  2.0335,  1.9254,  1.7861,  1.5357,  1.1784,  0.7628, ...
                      0.1500, -0.1401,  0.1272,  0.3884,  0.2761,  0.1707,  0.2279, ...
                     -0.2417, -1.7802, -3.0381, -2.3511, -0.9864, -0.0989,  0.1854, ...
                     -0.1282, -0.8028, -0.7363, -0.0093,  0.3132,  0.1894, -0.1526, ...
                     -0.3663, -0.3873, -0.3993, -0.6790, -1.0615, -1.4633, -1.5687, ...
                     -1.7183, -1.7313, -1.8324, -1.5909, -1.3199, -0.8634, -0.5542, ...
                     -0.1961, -0.0552,  0.0646, -0.0109, -0.0662, -0.2539, -0.3869, ...
                     -0.6379, -0.8404, -1.1662, -1.3930, -1.6136, -1.5706, -1.4266, ...
                     -1.0460, -0.7244, -0.3006, -0.1300,  0.1491,  0.0972,  0.2463, ...
                      0.0252,  0.1079, -0.1893, -0.1033, -0.3547, -0.2393, -0.4280, ...
                     -0.2735, -0.3670, -0.2033, -0.2261, -0.0821, -0.0754,  0.0634];

z_courty = [0.00000, 0.04912, 0.10060, 0.15470, 0.21140, 0.27090, 0.33330, 0.39880, ...
            0.46750, 0.53960, 0.61520, 0.69450, 0.77780, 0.86510, 0.95670, 1.05300, ...
            1.15400, 1.25900, 1.37000, 1.48700, 1.60900, 1.73700, 1.87100, 2.01300, ...
            2.16000, 2.31600, 2.47900, 2.64900, 2.82900, 3.01700, 3.21400, 3.42100, ...
            3.63800, 3.86600, 4.10500, 4.35600, 4.61900, 4.89500, 5.18400, 5.48800, ...
            5.80700, 6.14100, 6.49200, 6.85900, 7.24600, 7.65000, 8.07500, 8.52100, ...
            8.98900, 9.50000];
% z_courty = linspace(0,9.5,50);

phi_courty = [0.0499886, 0.0582622, 0.0678333, 0.0788739, 0.0915889, 0.1061913, 0.1229119, ...
              0.1419961, 0.1637082, 0.1883230, 0.2161014, 0.2473183, 0.2822266, 0.3210551, ...
              0.3639784, 0.4111301, 0.4623273, 0.5172858, 0.5752659, 0.6351540, 0.6950232, ...
              0.7529284, 0.8063160, 0.8520859, 0.8920522, 0.9305764, 0.9682031, 1.0058810, ...
              1.0444020, 1.0848160, 1.1282190, 1.1745120, 1.2226670, 1.2723200, 1.3231350, ...
              1.3743020, 1.4247480, 1.4730590, 1.5174060, 1.5552610, 1.5833640, 1.5976390, ...
              1.5925270, 1.5613110, 1.4949610, 1.3813710, 1.2041510, 0.9403100, 0.5555344, ...
              0.0000000];
          
Ncc07 = length(temperature_cc07);
Ncourty = length(z_courty);
c1 = 0.4;
c2 = 10.0;
TT0 = 1e5;
TTC = 1d6;
alpha1 = 0.15;
f_courty = 1;
TT = T;
lTT = log10(TT);

% This is a simple model to take into account the ionization background 
% on metal cooling (calibrated using CLOUDY).
if(Madau == 1)
    if (z <= 0.0 | z >= z_courty(50))
        ux=0.0;
    else
        iZ = 1 + round(z./z_courty(50)*49);
        iZ = min(iZ, 49);
        iZ = max(iZ, 1);
        delta_z = z_courty(iZ+1) - z_courty(iZ);
        ux = 1e-4 .* ( phi_courty(iZ+1).*(z-z_courty(iZ))./delta_z + ...
            phi_courty(iZ).*(z_courty(iZ+1)-z)./delta_z ) ./ nH;
    end
else
    ux = 1e-4 .* J0simple ./ nH;
end
g_courty       =   c1.*(TT/TT0).^alpha1         + c2.*exp(-TTC/TT);
g_courty_prime = ( c1.*alpha1.*(TT/TT0).^alpha1 + c2.*exp(-TTC/TT).*TTC/TT ) ./ TT;
f_courty       = 1  ./ (1 + ux./g_courty);
f_courty_prime = ( ux ./ (g_courty.*((1+ux./g_courty).^2)) ) .* ( g_courty_prime ./ g_courty );

if(lTT >= temperature_cc07(91))
    metal_tot1   = 1e-100;
    metal_tot2   = 1e-100;
    metal_prime1 = 0;
    metal_prime2 = 0;
elseif(lTT >= 1.0)
    lcool1       = -100;
    lcool1_prime = 0;
    if(lTT >= temperature_cc07(1))
        iT = 1 + round( 90.0 .* (lTT-temperature_cc07(1)) ./ (temperature_cc07(91)-temperature_cc07(1)) );
        iT = min(iT, 90);
        iT = max(iT, 1);
        deltaT = temperature_cc07(iT+1) - temperature_cc07(iT);
        lcool1 = excess_cooling_cc07(iT+1) .* (lTT-temperature_cc07(iT))   ./ deltaT + ...
                 excess_cooling_cc07(iT)   .* (temperature_cc07(iT+1)-lTT) ./ deltaT;
        lcool1_prime = excess_prime_cc07(iT+1) .* (lTT-temperature_cc07(iT))   ./ deltaT + ...
                       excess_prime_cc07(iT)   .* (temperature_cc07(iT+1)-lTT) ./ deltaT;
    end
% Fine structure cooling from infrared lines
    lcool2       = -31.522879 + 2.0*lTT - 20.0./TT - TT.*4.342944e-5;
    lcool2_prime = 2 + (20./TT - TT.*4.342944e-5).*log(10);
% Total metal cooling and temperature derivative
    metal_tot1   =    10.^lcool1;
    metal_tot2   =    10.^lcool2;
    metal_prime1 = ( (10.^lcool1).*lcool1_prime ) ./ (metal_tot1+metal_tot2);
    metal_prime2 = ( (10.^lcool2).*lcool2_prime ) ./ (metal_tot1+metal_tot2);
    metal_prime1 = metal_prime1.*f_courty + (metal_tot1+metal_tot2).*f_courty_prime;
    metal_prime2 = metal_prime2.*f_courty + (metal_tot1+metal_tot2).*f_courty_prime;
    metal_tot1   = metal_tot1.*f_courty;
    metal_tot2   = metal_tot2.*f_courty;
else
    metal_tot1 = 1e-100;
    metal_tot2 = 1e-100;
    metal_prime1 = 0;
    metal_prime2 = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%net heating and cooling
%%%%%%%%%%%%%%%%%%%%%%%%
%Bremstrahlung
cb1 = cool_bre(1) * ne * nHII   / nH^2;
cb2 = cool_bre(2) * ne * nHeII  / nH^2;
cb3 = cool_bre(3) * ne * nHeIII / nH^2;

%Ionization cooling
ci1 = cool_ion(1) * ne *nHI   / nH^2;
ci2 = cool_ion(2) * ne *nHeI  / nH^2;
ci3 = cool_ion(3) * ne *nHeII / nH^2;

%Recombination cooling
cr1 = cool_rec(1) * ne * nHII   / nH^2;
cr2 = cool_rec(2) * ne * nHeII  / nH^2;
cr3 = cool_rec(3) * ne * nHeIII / nH^2;

%Dielectric recombination cooling
cd  = cool_die * ne * nHeII / nH^2;

%Line cooling
ce1 = cool_exc(1) *ne * nHI   / nH^2;
ce2 = cool_exc(2) *ne * nHeI  / nH^2;
ce3 = cool_exc(3) *ne * nHeII / nH^2;

%Radiative heating
ch1 = h_rad(1) * nHI   / nH^2;
ch2 = h_rad(2) * nHeI  / nH^2;
ch3 = h_rad(3) * nHeII / nH^2;

%Compton cooling
coc = cool_com .* ne ./ nH^2;

%Compton heating
coh = heat_com .* ne ./ nH^2;

%Metal cooling
cM1 = metal_tot1 .* zmet;
cM2 = metal_tot2 .* zmet;

%Total cooling and heating rates
heat_tot = ch1+ch2+ch3 + coh;
cool_tot = cb1+cb2+cb3 + ci1+ci2+ci3 + cr1+cr2+cr3 + cd + ce1+ce2+ce3 + coc + cM1 + cM2;

heat = [ch1 ch2 ch3 coh];
cool = [cb1 cb2 cb3 ci1 ci2 ci3 cr1 cr2 cr3 cd ce1 ce2 ce3 coc cM1 cM2];
