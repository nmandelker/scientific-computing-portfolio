function dy = clump_evolution_exact(t, y, alpha, eta, eps, mu, etas, td_tff, tmig_td)
% First index is \phi_g and second is \phi_s
dy=zeros(3,1);
dy(1)= 0.5 .* alpha .* ( y(1)+y(2) ) - td_tff * (mu + eta) .* eps .* y(1) ; 
dy(2)= td_tff * (mu - etas) .* eps .* y(1);
dy(3)= -( y(1) + y(2) ) ./ (tmig_td);

%[T Y]=ode45(@(t,y)bathtub_exact(t,y,1,0,10^5,0.02,5,0.54),[16^(-1.5)*17.5
%17.5],[10^5 0])