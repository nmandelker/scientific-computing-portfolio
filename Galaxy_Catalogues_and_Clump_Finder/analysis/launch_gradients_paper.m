close all
fclose all
clear m Md_thresh z_thresh T_thresh prop
m = 1e9;
%T_thresh = [1 2 3];

%prop = [3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20, 23, 24, 41, 42, 43, 48, 100];
prop = [13];
T_thresh = 20;
z_thresh = [0.7 2.5];
Md_thresh = [10 11];
% z_thresh = [0.7, 2
%     0.7 2.5];
% Md_thresh = [9 10
%     10 11
%     9 11];

for i=1:length(prop)
    for k=1:length(T_thresh)
        for l=1:length(Md_thresh(:,1))
            for ll=1:length(z_thresh(:,1))
%                 if(prop(i)==100)
%                     norm = 0;
%                 else
%                     norm = 1;
%                 end

%                 clump_gradients_paper(is3, norm_is3, es3, norm_es3, nis3, nes3, ...
%                     prop(i), 1, 3, 0, 0, ...
%                     7.5, -m, 0.001, m, T_thresh(k), Md_thresh(l,:), z_thresh(ll,:),1);
%                 
%                 clump_gradients_paper(is3, norm_is3, es3, norm_es3, nis3, nes3, ...
%                     prop(i), 0, 3, 0, 0, ...
%                     7.5, -m, 0.001, m, T_thresh(k), Md_thresh(l,:), z_thresh(ll,:),1);
%                 
                clump_gradients_paper(is3, norm_is3, es3, norm_es3, nis3, nes3, ...
                    prop(i), 1, 3, 0, 0, ...
                    7.0, -m, 0.001, m, T_thresh(k), Md_thresh(l,:), z_thresh(ll,:),1);
                
                clump_gradients_paper(is3, norm_is3, es3, norm_es3, nis3, nes3, ...
                    prop(i), 0, 3, 0, 0, ...
                    7.0, -m, 0.001, m, T_thresh(k), Md_thresh(l,:), z_thresh(ll,:),1);

                clump_gradients_paper(is3, norm_is3, es3, norm_es3, nis3, nes3, ...
                    prop(i), 1, 3, 0, 0, ...
                    7.0, -m, 10^(-1.5), m, T_thresh(k), Md_thresh(l,:), z_thresh(ll,:),1);
                
                clump_gradients_paper(is3, norm_is3, es3, norm_es3, nis3, nes3, ...
                    prop(i), 0, 3, 0, 0, ...
                    7.0, -m, 10^(-1.5), m, T_thresh(k), Md_thresh(l,:), z_thresh(ll,:),1);
%                 
%                 clump_gradients_paper(is3, norm_is3, es3, norm_es3, nis3, nes3, ...
%                     prop(i), 1, 3, 0, 0, ...
%                     -m, -m, 10^(-1.5), m, T_thresh(k), Md_thresh(l,:), z_thresh(ll,:),1);
%                 
%                 clump_gradients_paper(is3, norm_is3, es3, norm_es3, nis3, nes3, ...
%                     prop(i), 0, 3, 0, 0, ...
%                     -m, -m, 10^(-1.5), m, T_thresh(k), Md_thresh(l,:), z_thresh(ll,:),1);
            end
        end
    end
end

clear prop norm T_thresh zthresh m