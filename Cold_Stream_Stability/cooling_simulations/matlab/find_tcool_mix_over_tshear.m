function [tcool_over_tsc, tcool_over_tdis, tcool_over_tdis2, tcool_over_tkh, Rs_crit_map] = find_tcool_mix_over_tshear(Mb, delta, ns, Rs)

KB = 1.38e-16;
mproton = 1.67e-24;
Myr = 60*60*24*365.25*1e6;
kpc = 3.086e21;
Zsolar = 0.02;

red = 2;    % for UVB
Zs = 0.03*Zsolar;  % absolute units
Zb = 0.10*Zsolar;  % absolute units
Zmix = sqrt(Zs*Zb);

[Ts,mus] = find_cooling_equilibrium2(ns, red, Zs/Zsolar, 1);
Cs = sqrt( (5/3) * KB * Ts / (mus * mproton) ); % cm/s
Cs = Cs * Myr / kpc; % kpc/Myr
tsc= 2 * Rs / Cs;    % Myr

Mtot  = ( sqrt(delta) / (1+sqrt(delta)) ) * Mb;
alpha = 0.21 * ( 0.8*exp(-3*Mtot^2) + 0.2 );

nb = ns / delta;
Tb = Ts * delta;
Tmix = sqrt(Ts*Tb);
nmix = sqrt(ns*nb);
[nspec, mu, cool, heat, cool_tot, heat_tot] = cooling_equilibrium(Tmix, nmix, red, Zmix/Zsolar, 1, 1);
tcool = ( 1.5*KB*Tmix / (nmix*(cool_tot-heat_tot)) ) / Myr;    % Myr
tcool_over_tsc = tcool / tsc;
        
tdis = Rs / ( alpha * Mb*Cs*sqrt(delta) );
tcool_over_tdis = tcool / tdis;
Rs_crit_map = Rs .* tcool / tdis;

tdis2 = (1+sqrt(delta)) * Rs / ( alpha * Mb*Cs*sqrt(delta) );
tcool_over_tdis2 = tcool / tdis2;
Rs_crit_map = Rs .* tcool / tdis2;

tkh = sqrt(delta) * Rs / ( Mb*Cs*sqrt(delta) );
tcool_over_tkh = tcool / tkh;
Rs_crit_map = Rs .* tcool / tkh;