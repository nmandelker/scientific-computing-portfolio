A1 = importdata(strcat('./Behroozi/smmr/c_smmr_z1.00_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat'),' ',1);
Behroozi_1 = A1.data(:,:);
clear A1

A1 = importdata(strcat('./Behroozi/smmr/c_smmr_z2.00_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat'),' ',1);
Behroozi_2 = A1.data(:,:);
clear A1

A1 = importdata(strcat('./Behroozi/smmr/c_smmr_z3.00_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat'),' ',1);
Behroozi_3 = A1.data(:,:);
clear A1

A1 = importdata(strcat('./Behroozi/smmr/c_smmr_z4.00_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat'),' ',1);
Behroozi_4 = A1.data(:,:);
clear A1