%type = 1:12;
type = [6 11];
for i=1:length(type)
    specific_rates(is, nis, ind, 7, type(i), [100 3000], [8.5 9.5], [1 2.5], 50, 1)
    
    %specific_rates(is, nis, ind, 8, type(i), [50 3000], [8 8.5], [1 2], 50, 1)
    
    %specific_rates(is, nis, ind, 19, type(i), [50 3000], [8.0 9.5], [3 5], 20, 0)
end