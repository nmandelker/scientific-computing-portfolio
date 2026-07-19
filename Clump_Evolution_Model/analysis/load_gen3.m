function [is3, es3, bulge3, nis3, nes3, nbulge3, norm_is3, norm_es3, norm_bulge3] = load_gen3()

is_1=load('./V07/in_situ.out');
is_2=load('./V08/in_situ.out');
is_3=load('./V11/in_situ.out');
is_4=load('./V12/in_situ.out');
is_5=load('./V14/in_situ.out');
is_6=load('./V19/in_situ.out');
is_7=load('./V25/in_situ.out');
is_8=load('./V26/in_situ.out');
is_9=load('./V27/in_situ.out');

nis_1 = size(is_1,1);
nis_2 = size(is_2,1);
nis_3 = size(is_3,1);
nis_4 = size(is_4,1);
nis_5 = size(is_5,1);
nis_6 = size(is_6,1);
nis_7 = size(is_7,1);
nis_8 = size(is_8,1);
nis_9 = size(is_9,1);

is3 = [is_1
    is_2
    is_3
    is_4
    is_5
    is_6
    is_7
    is_8
    is_9];
nt=[nis_1 nis_2 nis_3 nis_4 nis_5 nis_6 nis_7 nis_8 nis_9];
nis3=[0 nis_1];
for i=2:length(nt)
    nis3(i+1) = nis3(i) + nt(i);
end
clear is_1 is_2 is_3 is_4 is_5 is_6 is_7 is_8 is_9
clear nis_1 nis_2 nis_3 nis_4 nis_5 nis_6 nis_7 nis_8 nis_9

es_1=load('./V07/ex_situ.out');
es_2=load('./V08/ex_situ.out');
es_3=load('./V11/ex_situ.out');
es_4=load('./V12/ex_situ.out');
es_5=load('./V14/ex_situ.out');
es_6=load('./V19/ex_situ.out');
es_7=load('./V25/ex_situ.out');
es_8=load('./V26/ex_situ.out');
es_9=load('./V27/ex_situ.out');

nes_1 = size(es_1,1);
nes_2 = size(es_2,1);
nes_3 = size(es_3,1);
nes_4 = size(es_4,1);
nes_5 = size(es_5,1);
nes_6 = size(es_6,1);
nes_7 = size(es_7,1);
nes_8 = size(es_8,1);
nes_9 = size(es_9,1);

es3 = [es_1
    es_2
    es_3
    es_4
    es_5
    es_6
    es_7
    es_8
    es_9];

nt=[nes_1 nes_2 nes_3 nes_4 nes_5 nes_6 nes_7 nes_8 nes_9];
nes3=[0 nes_1];
for i=2:length(nt)
    nes3(i+1) = nes3(i) + nt(i);
end
clear es_1 es_2 es_3 es_4 es_5 es_6 es_7 es_8 es_9
clear nes_1 nes_2 nes_3 nes_4 nes_5 nes_6 nes_7 nes_8 nes_9

bulge_1=load('./V07/bulge.out');
bulge_2=load('./V08/bulge.out');
bulge_3=load('./V11/bulge.out');
bulge_4=load('./V12/bulge.out');
bulge_5=load('./V14/bulge.out');
bulge_6=load('./V19/bulge.out');
bulge_7=load('./V25/bulge.out');
bulge_8=load('./V26/bulge.out');
bulge_9=load('./V27/bulge.out');

nbulge_1 = size(bulge_1,1);
nbulge_2 = size(bulge_2,1);
nbulge_3 = size(bulge_3,1);
nbulge_4 = size(bulge_4,1);
nbulge_5 = size(bulge_5,1);
nbulge_6 = size(bulge_6,1);
nbulge_7 = size(bulge_7,1);
nbulge_8 = size(bulge_8,1);
nbulge_9 = size(bulge_9,1);

bulge3 = [bulge_1
    bulge_2
    bulge_3
    bulge_4
    bulge_5
    bulge_6
    bulge_7
    bulge_8
    bulge_9];

nt=[nbulge_1 nbulge_2 nbulge_3 nbulge_4 nbulge_5 nbulge_6 nbulge_7 nbulge_8 nbulge_9];
nbulge3=[0 nbulge_1];
for i=2:length(nt)
    nbulge3(i+1) = nbulge3(i) + nt(i);
end
clear bulge_1 bulge_2 bulge_3 bulge_4 bulge_5 bulge_6 bulge_7 bulge_8 bulge_9 
clear nbulge_1 nbulge_2 nbulge_3 nbulge_4 nbulge_5 nbulge_6 nbulge_7 nbulge_8 nbulge_9

norm_is_1=load('./V07/normalized_in_situ.out');
norm_is_2=load('./V08/normalized_in_situ.out');
norm_is_3=load('./V11/normalized_in_situ.out');
norm_is_4=load('./V12/normalized_in_situ.out');
norm_is_5=load('./V14/normalized_in_situ.out');
norm_is_6=load('./V19/normalized_in_situ.out');
norm_is_7=load('./V25/normalized_in_situ.out');
norm_is_8=load('./V26/normalized_in_situ.out');
norm_is_9=load('./V27/normalized_in_situ.out');

norm_is3 = [norm_is_1
    norm_is_2
    norm_is_3
    norm_is_4
    norm_is_5
    norm_is_6
    norm_is_7
    norm_is_8
    norm_is_9];

clear norm_is_1 norm_is_2 norm_is_3 norm_is_4 norm_is_5 norm_is_6 norm_is_7 norm_is_8 norm_is_9

norm_es_1=load('./V07/normalized_ex_situ.out');
norm_es_2=load('./V08/normalized_ex_situ.out');
norm_es_3=load('./V11/normalized_ex_situ.out');
norm_es_4=load('./V12/normalized_ex_situ.out');
norm_es_5=load('./V14/normalized_ex_situ.out');
norm_es_6=load('./V19/normalized_ex_situ.out');
norm_es_7=load('./V25/normalized_ex_situ.out');
norm_es_8=load('./V26/normalized_ex_situ.out');
norm_es_9=load('./V27/normalized_ex_situ.out');

norm_es3 = [norm_es_1
    norm_es_2
    norm_es_3
    norm_es_4
    norm_es_5
    norm_es_6
    norm_es_7
    norm_es_8
    norm_es_9];

clear norm_es_1 norm_es_2 norm_es_3 norm_es_4 norm_es_5 norm_es_6 norm_es_7 norm_es_8 norm_es_9

norm_bulge_1=load('./V07/normalized_bulge.out');
norm_bulge_2=load('./V08/normalized_bulge.out');
norm_bulge_3=load('./V11/normalized_bulge.out');
norm_bulge_4=load('./V12/normalized_bulge.out');
norm_bulge_5=load('./V14/normalized_bulge.out');
norm_bulge_6=load('./V19/normalized_bulge.out');
norm_bulge_7=load('./V25/normalized_bulge.out');
norm_bulge_8=load('./V26/normalized_bulge.out');
norm_bulge_9=load('./V27/normalized_bulge.out');

norm_bulge3 = [norm_bulge_1
    norm_bulge_2
    norm_bulge_3
    norm_bulge_4
    norm_bulge_5
    norm_bulge_6
    norm_bulge_7
    norm_bulge_8
    norm_bulge_9];

clear norm_bulge_1 norm_bulge_2 norm_bulge_3 norm_bulge_4 norm_bulge_5 norm_bulge_6 norm_bulge_7 norm_bulge_8 norm_bulge_9
clear nt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = find(is3(:,5)<1e-6); % No stars
is3(b,5)  = 1e-6; % correct stellar mass
is3(b,10) = 1e-6; % correct stellar surface density
is3(b,12) = 1e-6; % correct stellar age
is3(b,14) = 1e-6; % correct stellar metalicity
is3(b,15) = 1e-6; % correct SFR
is3(b,16) = 1e-6; % correct Sigma SFR
is3(b,17) = 1e-6; % correct sSFR
is3(b,18) = 1e-6; % correct tdep
norm_is3(b,5)  = 1e-6; % correct stellar mass
norm_is3(b,10) = 1e-6; % correct stellar surface density
norm_is3(b,12) = 1e-6; % correct stellar age
norm_is3(b,14) = 1e-6; % correct stellar metalicity
norm_is3(b,15) = 1e-6; % correct SFR
norm_is3(b,16) = 1e-6; % correct Sigma SFR
norm_is3(b,17) = 1e-6; % correct sSFR
norm_is3(b,18) = 1e-6; % correct tdep

b = find(is3(:,15)<1e-6); % No SFR
is3(b,15) = 1e-6; % correct SFR
is3(b,16) = 1e-6; % correct Sigma SFR
is3(b,17) = 1e-6; % correct sSFR
is3(b,18) = 1e-6; % correct tdep
norm_is3(b,15) = 1e-6; % correct SFR
norm_is3(b,16) = 1e-6; % correct Sigma SFR
norm_is3(b,17) = 1e-6; % correct sSFR
norm_is3(b,18) = 1e-6; % correct tdep

b = find(is3(:,4)<1e-6); % No gas
is3(b,4)  = 1e-6; % correct gas mass
is3(b,7)  = 1e-6; % correct gas fraction
is3(b,13) = 1e-6; % correct gas metalicity
is3(b,18) = 1e-6; % correct tdep
norm_is3(b,4)  = 1e-6; % correct gas mass
norm_is3(b,7)  = 1e-6; % correct gas fraction
norm_is3(b,13) = 1e-6; % correct gas metalicity
norm_is3(b,18) = 1e-6; % correct tdep

b = find(isnan(is3(:,29))); % Bad local t_dyn
is3(b,29) = 1e10;
norm_is3(b,29) = 1e10;
b = find(isnan(is3(:,30))); % Bad global t_dyn
is3(b,30) = 1e10;
norm_is3(b,30) = 1e10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = find(es3(:,5)<1e-6); % No stars
es3(b,5)  = 1e-6; % correct stellar mass
es3(b,10) = 1e-6; % correct stellar surface density
es3(b,12) = 1e-6; % correct stellar age
es3(b,14) = 1e-6; % correct stellar metalicity
es3(b,15) = 1e-6; % correct SFR
es3(b,16) = 1e-6; % correct Sigma SFR
es3(b,17) = 1e-6; % correct sSFR
es3(b,18) = 1e-6; % correct tdep
norm_es3(b,5)  = 1e-6; % correct stellar mass
norm_es3(b,10) = 1e-6; % correct stellar surface density
norm_es3(b,12) = 1e-6; % correct stellar age
norm_es3(b,14) = 1e-6; % correct stellar metalicity
norm_es3(b,15) = 1e-6; % correct SFR
norm_es3(b,16) = 1e-6; % correct Sigma SFR
norm_es3(b,17) = 1e-6; % correct sSFR
norm_es3(b,18) = 1e-6; % correct tdep

b = find(es3(:,15)<1e-6); % No SFR
es3(b,15) = 1e-6; % correct SFR
es3(b,16) = 1e-6; % correct Sigma SFR
es3(b,17) = 1e-6; % correct sSFR
es3(b,18) = 1e-6; % correct tdep
norm_es3(b,15) = 1e-6; % correct SFR
norm_es3(b,16) = 1e-6; % correct Sigma SFR
norm_es3(b,17) = 1e-6; % correct sSFR
norm_es3(b,18) = 1e-6; % correct tdep

b = find(es3(:,4)<1e-6); % No gas
es3(b,4)  = 1e-6; % correct gas mass
es3(b,7)  = 1e-6; % correct gas fraction
es3(b,13) = 1e-6; % correct gas metalicity
es3(b,18) = 1e-6; % correct tdep
norm_es3(b,4)  = 1e-6; % correct gas mass
norm_es3(b,7)  = 1e-6; % correct gas fraction
norm_es3(b,13) = 1e-6; % correct gas metalicity
norm_es3(b,18) = 1e-6; % correct tdep

b = find(isnan(es3(:,29))); % Bad local t_dyn
es3(b,29) = 1e10;
norm_es3(b,29) = 1e10;
b = find(isnan(es3(:,30))); % Bad global t_dyn
es3(b,30) = 1e10;
norm_es3(b,30) = 1e10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = find(bulge3(:,5)<1e-6); % No stars
bulge3(b,5)  = 1e-6; % correct stellar mass
bulge3(b,10) = 1e-6; % correct stellar surface density
bulge3(b,12) = 1e-6; % correct stellar age
bulge3(b,14) = 1e-6; % correct stellar metalicity
bulge3(b,15) = 1e-6; % correct SFR
bulge3(b,16) = 1e-6; % correct Sigma SFR
bulge3(b,17) = 1e-6; % correct sSFR
bulge3(b,18) = 1e-6; % correct tdep
norm_bulge3(b,5)  = 1e-6; % correct stellar mass
norm_bulge3(b,10) = 1e-6; % correct stellar surface density
norm_bulge3(b,12) = 1e-6; % correct stellar age
norm_bulge3(b,14) = 1e-6; % correct stellar metalicity
norm_bulge3(b,15) = 1e-6; % correct SFR
norm_bulge3(b,16) = 1e-6; % correct Sigma SFR
norm_bulge3(b,17) = 1e-6; % correct sSFR
norm_bulge3(b,18) = 1e-6; % correct tdep

b = find(bulge3(:,15)<1e-6); % No SFR
bulge3(b,15) = 1e-6; % correct SFR
bulge3(b,16) = 1e-6; % correct Sigma SFR
bulge3(b,17) = 1e-6; % correct sSFR
bulge3(b,18) = 1e-6; % correct tdep
norm_bulge3(b,15) = 1e-6; % correct SFR
norm_bulge3(b,16) = 1e-6; % correct Sigma SFR
norm_bulge3(b,17) = 1e-6; % correct sSFR
norm_bulge3(b,18) = 1e-6; % correct tdep

b = find(bulge3(:,4)<1e-6); % No gas
bulge3(b,4)  = 1e-6; % correct gas mass
bulge3(b,7)  = 1e-6; % correct gas fraction
bulge3(b,13) = 1e-6; % correct gas metalicity
bulge3(b,18) = 1e-6; % correct tdep
norm_bulge3(b,4)  = 1e-6; % correct gas mass
norm_bulge3(b,7)  = 1e-6; % correct gas fraction
norm_bulge3(b,13) = 1e-6; % correct gas metalicity
norm_bulge3(b,18) = 1e-6; % correct tdep

b = find(isnan(bulge3(:,29))); % Bad local t_dyn
bulge3(b,29) = 1e10;
norm_bulge3(b,29) = 1e10;
b = find(isnan(bulge3(:,30))); % Bad global t_dyn
bulge3(b,30) = 1e10;
norm_bulge3(b,30) = 1e10;

clear b i j

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalized radius
norm_rad_is3 = is3(:,21)./is3(:,19); %kpc
norm_rad_es3 = es3(:,21)./es3(:,19); %kpc
norm_rad_bulge3 = bulge3(:,21)./bulge3(:,19); %kpc
norm_is3(:,3) = norm_is3(:,3)./norm_rad_is3;
norm_es3(:,3) = norm_es3(:,3)./norm_rad_es3;
norm_bulge3(:,3) = norm_bulge3(:,3)./norm_rad_bulge3;

clear norm_rad_is3 norm_rad_es3 norm_rad_bulge3

% depletion time
is3(:,18) = 1./is3(:,18);
es3(:,18) = 1./es3(:,18);
bulge3(:,18) = 1./bulge3(:,18);
norm_is3(:,18) = 1./norm_is3(:,18);
norm_es3(:,18) = 1./norm_es3(:,18);
norm_bulge3(:,18) = 1./norm_bulge3(:,18);

% absolute value height
is3(:,20) = abs(is3(:,20));
es3(:,20) = abs(es3(:,20));
bulge3(:,20) = abs(bulge3(:,20));
norm_is3(:,20) = abs(norm_is3(:,20));
norm_es3(:,20) = abs(norm_es3(:,20));
norm_bulge3(:,20) = abs(norm_bulge3(:,20));
is3(:,22) = abs(is3(:,22));
es3(:,22) = abs(es3(:,22));
bulge3(:,22) = abs(bulge3(:,22));
norm_is3(:,22) = abs(norm_is3(:,22));
norm_es3(:,22) = abs(norm_es3(:,22));
norm_bulge3(:,22) = abs(norm_bulge3(:,22));

% Correct free fall time
is3(:,28) = is3(:,28).*sqrt(3*pi/32);
es3(:,28) = es3(:,28).*sqrt(3*pi/32);
bulge3(:,28) = bulge3(:,28).*sqrt(3*pi/32);
norm_is3(:,28) = norm_is3(:,28).*sqrt(3*pi/32);
norm_es3(:,28) = norm_es3(:,28).*sqrt(3*pi/32);
norm_bulge3(:,28) = norm_bulge3(:,28).*sqrt(3*pi/32);

% absolute value dynamical times
is3(:,29) = abs(is3(:,29));
es3(:,29) = abs(es3(:,29));
bulge3(:,29) = abs(bulge3(:,29));
norm_is3(:,29) = abs(norm_is3(:,29));
norm_es3(:,29) = abs(norm_es3(:,29));
norm_bulge3(:,29) = abs(norm_bulge3(:,29));

is3(:,30) = abs(is3(:,30));
es3(:,30) = abs(es3(:,30));
bulge3(:,30) = abs(bulge3(:,30));
norm_is3(:,30) = abs(norm_is3(:,30));
norm_es3(:,30) = abs(norm_es3(:,30));
norm_bulge3(:,30) = abs(norm_bulge3(:,30));

% Fix Mstars_in, Mstars_out and Mstars_formed
b=find( isinf(is3(:,44)) | isnan(is3(:,44)) );
is3(b,44) = 0;
norm_is3(b,44)=0;
b=find( isinf(es3(:,44)) | isnan(es3(:,44)));
es3(b,44) = 0;
norm_es3(b,44)=0;
b=find(isinf(bulge3(:,44)) | isnan(bulge3(:,44)));
bulge3(b,44) = 0;
norm_bulge3(b,44)=0;

b=find( isinf(is3(:,45)) | isnan(is3(:,45)) );
is3(b,45) = 0;
norm_is3(b,45)=0;
b=find( isinf(es3(:,45)) | isnan(es3(:,45)));
es3(b,45) = 0;
norm_es3(b,45)=0;
b=find(isinf(bulge3(:,45)) | isnan(bulge3(:,45)));
bulge3(b,45) = 0;
norm_bulge3(b,45)=0;

b=find( isinf(is3(:,46)) | isnan(is3(:,46)) );
is3(b,46) = 0;
norm_is3(b,46)=0;
b=find( isinf(es3(:,46)) | isnan(es3(:,46)));
es3(b,46) = 0;
norm_es3(b,46)=0;
b=find(isinf(bulge3(:,46)) | isnan(bulge3(:,46)));
bulge3(b,46) = 0;
norm_bulge3(b,46)=0;

% lifetime in Myr, normalized by free fall time
is3(:,48) = 1000.*is3(:,48);
es3(:,48) = 1000.*es3(:,48);
bulge3(:,48) = 1000.*bulge3(:,48);
norm_is3(:,48) = is3(:,48)./is3(:,28);
norm_es3(:,48) = es3(:,48)./es3(:,28);
norm_bulge3(:,48) = bulge3(:,48)./bulge3(:,28);
[max(norm_is3(:,48)),min(norm_is3(:,48)),length(find(isnan(norm_is3(:,48))))]