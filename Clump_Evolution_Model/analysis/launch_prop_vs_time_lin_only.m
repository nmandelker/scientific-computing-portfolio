prop = [4 5 6 7 15 21 22 49 50 65 66 67 68 69 70 71 72 73 74 80 81 82 83 84 85 92 93 94 95 96];
for i=1:length(prop)
    prop_vs_time_lin_only(is, nis, ind, 7, prop(i), [100 3000], [8.5 9.5], [1 2.5], 50, 80)
    
    prop_vs_time_lin_only(is, nis, ind, 8, prop(i), [100 3000], [8.5 9.5], [1 2], 70, 80)
    
    prop_vs_time_lin_only(is, nis, ind, 19, prop(i), [50 3000], [8.0 9.5], [3 5], 20, 30)
end