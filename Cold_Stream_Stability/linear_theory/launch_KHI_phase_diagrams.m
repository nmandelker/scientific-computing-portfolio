function [a_over_Rvir, lambda_over_a, Tvir_over_Tkh, coherence] = launch_KHI_phase_diagrams()

Mh = 0.750:0.100:2.250;
Del = 20:2:100;
Mh2 = 0.75:0.01:2.25;

a_over_Rvir = [0.01 0.1];
lambda_over_a = [0.1 1 2];

Tvir_over_Tkh = zeros(length(Del), length(Mh), length(a_over_Rvir), length(lambda_over_a));
coherence = zeros(length(Mh2), length(a_over_Rvir), length(lambda_over_a));
for i=1:length(a_over_Rvir)
    for j=1:length(lambda_over_a)
        [a_over_Rvir(i), lambda_over_a(j)]
        
        [Tvir_over_Tkh(:,:,i,j), coherence(:,i,j)] = KHI_phase_diagrams(a_over_Rvir(i), lambda_over_a(j));
%         filename = strcat('./figures/a_over_Rvir_',num2str(a_over_Rvir(i),'%4.2f'),'__lambda_over_a_',num2str(lambda_over_a(j),'%4.2f'),'.jpg');
%         saveas(gcf,filename);
%         close all
%         fclose all
    end
end