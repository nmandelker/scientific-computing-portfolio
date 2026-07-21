function dy = halo_decel_cooling_NFW(x, y, eta, tau, c, beta)

fc = log(1+c) - c/(1+c);
alpha = 5*beta/8;
cx = c*x;
fcx = log(1+cx) - cx/(1+cx);

dy=zeros(2,1);
% derive rs(r) from constant line mass
% dy(1) = 2 * eta * tau * y(1)^(0.5) / ( y(2) * (x^alpha) ) - 2 * (fcx/fc) / x^2;
% dy(2) = -eta * tau / ( y(1)^(0.5) * (x^alpha) );
% derive rs(r) from self-consistent line mass
dy(1) = 2 * eta * tau * y(1)^(0.5) / ( y(2)^(11/8) * (x^alpha) ) - 2 * (fcx/fc) / x^2;
dy(2) = -eta * tau / ( y(1)^(0.5) * y(2)^(3/8) * (x^alpha) );