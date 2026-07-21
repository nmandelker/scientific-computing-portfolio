function [J21] = find_cooling_equilibrium(T, nH, z, zmet, Madau)

J0_left  = 1e-20;
J0_right = 1e20;
err_J0 = 1;
while(err_J0 > 1e-8)
    J0_multiply = 0.5 .* (J0_left + J0_right);
    [nspec, mu, cool, heat, cool_tot, heat_tot] = cooling_equilibrium(T, nH, z, zmet, Madau, J0_multiply);
    diff = (cool_tot-heat_tot);
    if(diff>0)
        J0_left  = 0.5 .* (J0_left + J0_right);
    else
        J0_right = 0.5 .* (J0_left + J0_right);
    end
    err_J0 = abs(J0_right - J0_left) ./ J0_left;
end

J21 = 0.5 .* (J0_left + J0_right);
