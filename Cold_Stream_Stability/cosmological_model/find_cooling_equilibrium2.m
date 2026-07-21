function [T,mu] = find_cooling_equilibrium2(nH, z, zmet, Madau)

T2_left  = 1e2;
T2_right = 1e9;
err_T2 = 1;

boost = 1;
% boost = exp(-nH/0.01);

while(err_T2 > 1e-10)
    T2 = 0.5 .* (T2_left + T2_right);
%     [nspec, mu, cool, heat, cool_tot, heat_tot] = cooling_equilibrium(T2, nH, z, zmet, Madau, 1);
    [nspec, mu, cool, heat, cool_tot, heat_tot] = cooling_equilibrium_T2_input(T2, nH, z, zmet, Madau, boost);
    diff = (cool_tot-heat_tot);
    if(diff<0)
        T2_left  = 0.5 .* (T2_left + T2_right);
    else
        T2_right = 0.5 .* (T2_left + T2_right);
    end
    err_T2 = abs(T2_right - T2_left) ./ T2_left;
end

% T = 0.5 .* (T2_left + T2_right);
T = 0.5 .* (T2_left + T2_right) * mu;
mu = mu;
