%prop = [4 5 6 7 15 21 22 49 50 65 66 67 68 69 70 71 72 73 74 80 81 82 83 84 85 92 93 94 95 96];
% alpha, eta1, eps, mu, etas, tmig_td, fg_i
% 4 5 6 7 15 177 18 65 71 74 84
prop = [21 22];
for i=1:length(prop)
     prop_vs_t_over_tdyn_rectangle(is, norm_is, nis, ind, 7, prop(i), [500 3000], [5 9.5], [-9,-3], [1 2.5], 1, 2, 30, 0, ...
         0.25, 1.0, 0.5, 0.8, 0.15, 20, 0.9)
%      prop_vs_t_over_tdyn_rectangle(is, norm_is, nis, ind, 7, prop(i), [50 3000], [8.0 9.5], [-1 9], [1 2.5], 1, 2, 30, 0, ...
%          0.25, 1.0, 0.5, 0.8, 0.15, 20, 0.9)
   % alpha~0.1-0.4(0.2) at t~100-600Myr and ~0.07-0.2(0.13) at t>600Myr
   % eta~0.55-2.5(1.1)  at t~100-400Myr and ~0.35-1(0.6)    at t>400Myr
   % eps~0.03-0.04(0.035)
   % td/tff~3.5-8.9(5.5)
    
    %prop_vs_t_over_tdyn2(is, nis, ind, 8, prop(i), [100 3000], [8.5 9.5], [1 2], 70, 80)
    
%    prop_vs_t_over_tdyn_rectangle(is, norm_is, nis, ind, 19, prop(i), [500 3000], [5.0 9.5], [-9 -3], [0.5 8.5], 1, 2, 17, 0, ...
%       0.42, 1.7, 0.15, 0.8, -0.6, 20, 0.8)
%    prop_vs_t_over_tdyn_rectangle(is, norm_is, nis, ind, 19, prop(i), [50 3000], [8.0 9.5], [-1, 9], [3 5], 1, 2, 17, 0, ...
%       0.42, 1.7, 0.15, 0.8, -0.6, 20, 0.8)
   % alpha~0.15-0.7(0.4)
   % eta~0.8-2.2(1.6)
   % eps~0.025-0.04(0.03-0.035)
   % td/tff~1.5-2.5(2)
end