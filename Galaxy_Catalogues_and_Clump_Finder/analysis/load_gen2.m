function [is2, es2, bulge2, nis2, nes2, nbulge2, norm_is2, norm_es2, norm_bulge2, disc2, ndisc2] = ...
    load_gen2()

is_1=load('./gen2/same_clumps/VELA01_in_situ.out');
is_2=load('./gen2/same_clumps/VELA02_in_situ.out');
is_3=load('./gen2/same_clumps/VELA03_in_situ.out');
is_4=load('./gen2/same_clumps/VELA05_in_situ.out');
is_5=load('./gen2/same_clumps/VELA06_in_situ.out');
is_6=load('./gen2/same_clumps/VELA07_in_situ.out');
is_7=load('./gen2/same_clumps/VELA08_in_situ.out');
is_8=load('./gen2/same_clumps/VELA09_in_situ.out');
is_9=load('./gen2/same_clumps/VELA10_in_situ.out');
is_10=load('./gen2/same_clumps/VELA11_in_situ.out');
is_11=load('./gen2/same_clumps/VELA12_in_situ.out');
is_12=load('./gen2/same_clumps/VELA13_in_situ.out');
is_13=load('./gen2/same_clumps/VELA14_in_situ.out');
is_14=load('./gen2/same_clumps/VELA15_in_situ.out');
is_15=load('./gen2/same_clumps/VELA16_in_situ.out');
is_16=load('./gen2/same_clumps/VELA19_in_situ.out');
is_17=load('./gen2/same_clumps/VELA21_in_situ.out');
is_18=load('./gen2/same_clumps/VELA23_in_situ.out');
is_19=load('./gen2/same_clumps/VELA25_in_situ.out');
is_20=load('./gen2/same_clumps/VELA26_in_situ.out');
is_21=load('./gen2/same_clumps/VELA27_in_situ.out');
is_22=load('./gen2/same_clumps/VELA28_in_situ.out');
is_23=load('./gen2/same_clumps/VELA29_in_situ.out');
is_24=load('./gen2/same_clumps/VELA30_in_situ.out');
is_25=load('./gen2/same_clumps/VELA31_in_situ.out');
is_26=load('./gen2/same_clumps/VELA32_in_situ.out');
is_27=load('./gen2/same_clumps/VELA33_in_situ.out');
is_28=load('./gen2/same_clumps/VELA34_in_situ.out');
is_29=load('./gen2/same_clumps/VELA35_in_situ.out');

nis_1 = size(is_1,1);
nis_2 = size(is_2,1);
nis_3 = size(is_3,1);
nis_4 = size(is_4,1);
nis_5 = size(is_5,1);
nis_6 = size(is_6,1);
nis_7 = size(is_7,1);
nis_8 = size(is_8,1);
nis_9 = size(is_9,1);
nis_10 = size(is_10,1);
nis_11 = size(is_11,1);
nis_12 = size(is_12,1);
nis_13 = size(is_13,1);
nis_14 = size(is_14,1);
nis_15 = size(is_15,1);
nis_16 = size(is_16,1);
nis_17 = size(is_17,1);
nis_18 = size(is_18,1);
nis_19 = size(is_19,1);
nis_20 = size(is_20,1);
nis_21 = size(is_21,1);
nis_22 = size(is_22,1);
nis_23 = size(is_23,1);
nis_24 = size(is_24,1);
nis_25 = size(is_25,1);
nis_26 = size(is_26,1);
nis_27 = size(is_27,1);
nis_28 = size(is_28,1);
nis_29 = size(is_29,1);

is2 = [is_1
    is_2
    is_3
    is_4
    is_5
    is_6
    is_7
    is_8
    is_9
    is_10
    is_11
    is_12
    is_13
    is_14
    is_15
    is_16
    is_17
    is_18
    is_19
    is_20
    is_21
    is_22
    is_23
    is_24
    is_25
    is_26
    is_27
    is_28
    is_29];
nt=[nis_1 nis_2 nis_3 nis_4 nis_5 nis_6 nis_7 nis_8 nis_9 nis_10 nis_11 ...
    nis_12 nis_13 nis_14 nis_15 nis_16 nis_17 nis_18 nis_19 nis_20 ...
    nis_21 nis_22 nis_23 nis_24 nis_25 nis_26 nis_27 nis_28 nis_29];
nis2=[0 nis_1];
for i=2:length(nt)
    nis2(i+1) = nis2(i) + nt(i);
end
clear is_1 is_2 is_3 is_4 is_5 is_6 is_7 is_8 is_9 is_10 is_11 is_12 is_13 is_14 is_15
clear is_16 is_17 is_18 is_19 is_20 is_21 is_22 is_23 is_24 is_25 is_26 is_27 is_28 is_29
clear nis_1 nis_2 nis_3 nis_4 nis_5 nis_6 nis_7 nis_8 nis_9 nis_10 nis_11 nis_12 nis_13
clear nis_14 nis_15 nis_16 nis_17 nis_18 nis_19 nis_20 nis_21 nis_22 nis_23 nis_24 nis_25 
clear nis_26 nis_27 nis_28 nis_29

es_1=load('./gen2/same_clumps/VELA01_ex_situ.out');
es_2=load('./gen2/same_clumps/VELA02_ex_situ.out');
es_3=load('./gen2/same_clumps/VELA03_ex_situ.out');
es_4=load('./gen2/same_clumps/VELA05_ex_situ.out');
es_5=load('./gen2/same_clumps/VELA06_ex_situ.out');
es_6=load('./gen2/same_clumps/VELA07_ex_situ.out');
es_7=load('./gen2/same_clumps/VELA08_ex_situ.out');
es_8=load('./gen2/same_clumps/VELA09_ex_situ.out');
es_9=load('./gen2/same_clumps/VELA10_ex_situ.out');
es_10=load('./gen2/same_clumps/VELA11_ex_situ.out');
es_11=load('./gen2/same_clumps/VELA12_ex_situ.out');
es_12=load('./gen2/same_clumps/VELA13_ex_situ.out');
es_13=load('./gen2/same_clumps/VELA14_ex_situ.out');
es_14=load('./gen2/same_clumps/VELA15_ex_situ.out');
es_15=load('./gen2/same_clumps/VELA16_ex_situ.out');
es_16=load('./gen2/same_clumps/VELA19_ex_situ.out');
es_17=load('./gen2/same_clumps/VELA21_ex_situ.out');
es_18=load('./gen2/same_clumps/VELA23_ex_situ.out');
es_19=load('./gen2/same_clumps/VELA25_ex_situ.out');
es_20=load('./gen2/same_clumps/VELA26_ex_situ.out');
es_21=load('./gen2/same_clumps/VELA27_ex_situ.out');
es_22=load('./gen2/same_clumps/VELA28_ex_situ.out');
es_23=load('./gen2/same_clumps/VELA29_ex_situ.out');
es_24=load('./gen2/same_clumps/VELA30_ex_situ.out');
es_25=load('./gen2/same_clumps/VELA31_ex_situ.out');
es_26=load('./gen2/same_clumps/VELA32_ex_situ.out');
es_27=load('./gen2/same_clumps/VELA33_ex_situ.out');
es_28=load('./gen2/same_clumps/VELA34_ex_situ.out');
es_29=load('./gen2/same_clumps/VELA35_ex_situ.out');

nes_1 = size(es_1,1);
nes_2 = size(es_2,1);
nes_3 = size(es_3,1);
nes_4 = size(es_4,1);
nes_5 = size(es_5,1);
nes_6 = size(es_6,1);
nes_7 = size(es_7,1);
nes_8 = size(es_8,1);
nes_9 = size(es_9,1);
nes_10 = size(es_10,1);
nes_11 = size(es_11,1);
nes_12 = size(es_12,1);
nes_13 = size(es_13,1);
nes_14 = size(es_14,1);
nes_15 = size(es_15,1);
nes_16 = size(es_16,1);
nes_17 = size(es_17,1);
nes_18 = size(es_18,1);
nes_19 = size(es_19,1);
nes_20 = size(es_20,1);
nes_21 = size(es_21,1);
nes_22 = size(es_22,1);
nes_23 = size(es_23,1);
nes_24 = size(es_24,1);
nes_25 = size(es_25,1);
nes_26 = size(es_26,1);
nes_27 = size(es_27,1);
nes_28 = size(es_28,1);
nes_29 = size(es_29,1);

es2 = [es_1
    es_2
    es_3
    es_4
    es_5
    es_6
    es_7
    es_8
    es_9
    es_10
    es_11
    es_12
    es_13
    es_14
    es_15
    es_16
    es_17
    es_18
    es_19
    es_20
    es_21
    es_22
    es_23
    es_24
    es_25
    es_26
    es_27
    es_28
    es_29];
nt=[nes_1 nes_2 nes_3 nes_4 nes_5 nes_6 nes_7 nes_8 nes_9 nes_10 nes_11 ...
    nes_12 nes_13 nes_14 nes_15 nes_16 nes_17 nes_18 nes_19 nes_20 nes_21 ...
    nes_22 nes_23 nes_24 nes_25 nes_26 nes_27 nes_28 nes_29];
nes2=[0 nes_1];
for i=2:length(nt)
    nes2(i+1) = nes2(i) + nt(i);
end
clear es_1 es_2 es_3 es_4 es_5 es_6 es_7 es_8 es_9 es_10 es_11 es_12 es_13 es_14 es_15
clear es_16 es_17 es_18 es_19 es_20 es_21 es_22 es_23 es_24 es_25 es_26 es_27 es_28 es_29
clear nes_1 nes_2 nes_3 nes_4 nes_5 nes_6 nes_7 nes_8 nes_9 nes_10 nes_11 nes_12 nes_13
clear nes_14 nes_15 nes_16 nes_17 nes_18 nes_19 nes_20 nes_21 nes_22 nes_23 nes_24 nes_25 
clear nes_26 nes_27 nes_28 nes_29

bulge_1=load('./gen2/same_clumps/VELA01_bulge.out');
bulge_2=load('./gen2/same_clumps/VELA02_bulge.out');
bulge_3=load('./gen2/same_clumps/VELA03_bulge.out');
bulge_4=load('./gen2/same_clumps/VELA05_bulge.out');
bulge_5=load('./gen2/same_clumps/VELA06_bulge.out');
bulge_6=load('./gen2/same_clumps/VELA07_bulge.out');
bulge_7=load('./gen2/same_clumps/VELA08_bulge.out');
bulge_8=load('./gen2/same_clumps/VELA09_bulge.out');
bulge_9=load('./gen2/same_clumps/VELA10_bulge.out');
bulge_10=load('./gen2/same_clumps/VELA11_bulge.out');
bulge_11=load('./gen2/same_clumps/VELA12_bulge.out');
bulge_12=load('./gen2/same_clumps/VELA13_bulge.out');
bulge_13=load('./gen2/same_clumps/VELA14_bulge.out');
bulge_14=load('./gen2/same_clumps/VELA15_bulge.out');
bulge_15=load('./gen2/same_clumps/VELA16_bulge.out');
bulge_16=load('./gen2/same_clumps/VELA19_bulge.out');
bulge_17=load('./gen2/same_clumps/VELA21_bulge.out');
bulge_18=load('./gen2/same_clumps/VELA23_bulge.out');
bulge_19=load('./gen2/same_clumps/VELA25_bulge.out');
bulge_20=load('./gen2/same_clumps/VELA26_bulge.out');
bulge_21=load('./gen2/same_clumps/VELA27_bulge.out');
bulge_22=load('./gen2/same_clumps/VELA28_bulge.out');
bulge_23=load('./gen2/same_clumps/VELA29_bulge.out');
bulge_24=load('./gen2/same_clumps/VELA30_bulge.out');
bulge_25=load('./gen2/same_clumps/VELA31_bulge.out');
bulge_26=load('./gen2/same_clumps/VELA32_bulge.out');
bulge_27=load('./gen2/same_clumps/VELA33_bulge.out');
bulge_28=load('./gen2/same_clumps/VELA34_bulge.out');
bulge_29=load('./gen2/same_clumps/VELA35_bulge.out');

nbulge_1 = size(bulge_1,1);
nbulge_2 = size(bulge_2,1);
nbulge_3 = size(bulge_3,1);
nbulge_4 = size(bulge_4,1);
nbulge_5 = size(bulge_5,1);
nbulge_6 = size(bulge_6,1);
nbulge_7 = size(bulge_7,1);
nbulge_8 = size(bulge_8,1);
nbulge_9 = size(bulge_9,1);
nbulge_10 = size(bulge_10,1);
nbulge_11 = size(bulge_11,1);
nbulge_12 = size(bulge_12,1);
nbulge_13 = size(bulge_13,1);
nbulge_14 = size(bulge_14,1);
nbulge_15 = size(bulge_15,1);
nbulge_16 = size(bulge_16,1);
nbulge_17 = size(bulge_17,1);
nbulge_18 = size(bulge_18,1);
nbulge_19 = size(bulge_19,1);
nbulge_20 = size(bulge_20,1);
nbulge_21 = size(bulge_21,1);
nbulge_22 = size(bulge_22,1);
nbulge_23 = size(bulge_23,1);
nbulge_24 = size(bulge_24,1);
nbulge_25 = size(bulge_25,1);
nbulge_26 = size(bulge_26,1);
nbulge_27 = size(bulge_27,1);
nbulge_28 = size(bulge_28,1);
nbulge_29 = size(bulge_29,1);

bulge2 = [bulge_1
    bulge_2
    bulge_3
    bulge_4
    bulge_5
    bulge_6
    bulge_7
    bulge_8
    bulge_9
    bulge_10
    bulge_11
    bulge_12
    bulge_13
    bulge_14
    bulge_15
    bulge_16
    bulge_17
    bulge_18
    bulge_19
    bulge_20
    bulge_21
    bulge_22
    bulge_23
    bulge_24
    bulge_25
    bulge_26
    bulge_27
    bulge_28
    bulge_29];
nt=[nbulge_1 nbulge_2 nbulge_3 nbulge_4 nbulge_5 nbulge_6 nbulge_7 nbulge_8 ...
    nbulge_9 nbulge_10 nbulge_11 nbulge_12 nbulge_13 nbulge_14 nbulge_15 ...
    nbulge_16 nbulge_17 nbulge_18 nbulge_19 nbulge_20 nbulge_21 nbulge_22 ...
    nbulge_23 nbulge_24 nbulge_25 nbulge_26 nbulge_27 nbulge_28 nbulge_29];
nbulge2=[0 nbulge_1];
for i=2:length(nt)
    nbulge2(i+1) = nbulge2(i) + nt(i);
end
clear bulge_1 bulge_2 bulge_3 bulge_4 bulge_5 bulge_6 bulge_7 bulge_8 bulge_9 
clear bulge_10 bulge_11 bulge_12 bulge_13 bulge_14 bulge_15 bulge_16 bulge_17
clear bulge_18 bulge_19 bulge_20 bulge_21 bulge_22 bulge_23 bulge_24 bulge_25 
clear bulge_26 bulge_27 bulge_28 bulge_29
clear nbulge_1 nbulge_2 nbulge_3 nbulge_4 nbulge_5 nbulge_6 nbulge_7 nbulge_8 
clear nbulge_9 nbulge_10 nbulge_11 nbulge_12 nbulge_13 nbulge_14 nbulge_15
clear nbulge_16 nbulge_17 nbulge_18 nbulge_19 nbulge_20 nbulge_21 nbulge_22 
clear nbulge_23 nbulge_24 nbulge_25 nbulge_26 nbulge_27 nbulge_28 nbulge_29

norm_is_1=load('./gen2/same_clumps/VELA01_normalized_in_situ.out');
norm_is_2=load('./gen2/same_clumps/VELA02_normalized_in_situ.out');
norm_is_3=load('./gen2/same_clumps/VELA03_normalized_in_situ.out');
norm_is_4=load('./gen2/same_clumps/VELA05_normalized_in_situ.out');
norm_is_5=load('./gen2/same_clumps/VELA06_normalized_in_situ.out');
norm_is_6=load('./gen2/same_clumps/VELA07_normalized_in_situ.out');
norm_is_7=load('./gen2/same_clumps/VELA08_normalized_in_situ.out');
norm_is_8=load('./gen2/same_clumps/VELA09_normalized_in_situ.out');
norm_is_9=load('./gen2/same_clumps/VELA10_normalized_in_situ.out');
norm_is_10=load('./gen2/same_clumps/VELA11_normalized_in_situ.out');
norm_is_11=load('./gen2/same_clumps/VELA12_normalized_in_situ.out');
norm_is_12=load('./gen2/same_clumps/VELA13_normalized_in_situ.out');
norm_is_13=load('./gen2/same_clumps/VELA14_normalized_in_situ.out');
norm_is_14=load('./gen2/same_clumps/VELA15_normalized_in_situ.out');
norm_is_15=load('./gen2/same_clumps/VELA16_normalized_in_situ.out');
norm_is_16=load('./gen2/same_clumps/VELA19_normalized_in_situ.out');
norm_is_17=load('./gen2/same_clumps/VELA21_normalized_in_situ.out');
norm_is_18=load('./gen2/same_clumps/VELA23_normalized_in_situ.out');
norm_is_19=load('./gen2/same_clumps/VELA25_normalized_in_situ.out');
norm_is_20=load('./gen2/same_clumps/VELA26_normalized_in_situ.out');
norm_is_21=load('./gen2/same_clumps/VELA27_normalized_in_situ.out');
norm_is_22=load('./gen2/same_clumps/VELA28_normalized_in_situ.out');
norm_is_23=load('./gen2/same_clumps/VELA29_normalized_in_situ.out');
norm_is_24=load('./gen2/same_clumps/VELA30_normalized_in_situ.out');
norm_is_25=load('./gen2/same_clumps/VELA31_normalized_in_situ.out');
norm_is_26=load('./gen2/same_clumps/VELA32_normalized_in_situ.out');
norm_is_27=load('./gen2/same_clumps/VELA33_normalized_in_situ.out');
norm_is_28=load('./gen2/same_clumps/VELA34_normalized_in_situ.out');
norm_is_29=load('./gen2/same_clumps/VELA35_normalized_in_situ.out');

norm_is2 = [norm_is_1
    norm_is_2
    norm_is_3
    norm_is_4
    norm_is_5
    norm_is_6
    norm_is_7
    norm_is_8
    norm_is_9
    norm_is_10
    norm_is_11
    norm_is_12
    norm_is_13
    norm_is_14
    norm_is_15
    norm_is_16
    norm_is_17
    norm_is_18
    norm_is_19
    norm_is_20
    norm_is_21
    norm_is_22
    norm_is_23
    norm_is_24
    norm_is_25
    norm_is_26
    norm_is_27
    norm_is_28
    norm_is_29];
clear norm_is_1 norm_is_2 norm_is_3 norm_is_4 norm_is_5 norm_is_6 norm_is_7
clear norm_is_8 norm_is_9 norm_is_10 norm_is_11 norm_is_12 norm_is_13
clear norm_is_14 norm_is_15 norm_is_16 norm_is_17 norm_is_18 norm_is_19
clear norm_is_20 norm_is_21 norm_is_22 norm_is_23 norm_is_24 norm_is_25 
clear norm_is_26 norm_is_27 norm_is_28 norm_is_29

norm_es_1=load('./gen2/same_clumps/VELA01_normalized_ex_situ.out');
norm_es_2=load('./gen2/same_clumps/VELA02_normalized_ex_situ.out');
norm_es_3=load('./gen2/same_clumps/VELA03_normalized_ex_situ.out');
norm_es_4=load('./gen2/same_clumps/VELA05_normalized_ex_situ.out');
norm_es_5=load('./gen2/same_clumps/VELA06_normalized_ex_situ.out');
norm_es_6=load('./gen2/same_clumps/VELA07_normalized_ex_situ.out');
norm_es_7=load('./gen2/same_clumps/VELA08_normalized_ex_situ.out');
norm_es_8=load('./gen2/same_clumps/VELA09_normalized_ex_situ.out');
norm_es_9=load('./gen2/same_clumps/VELA10_normalized_ex_situ.out');
norm_es_10=load('./gen2/same_clumps/VELA11_normalized_ex_situ.out');
norm_es_11=load('./gen2/same_clumps/VELA12_normalized_ex_situ.out');
norm_es_12=load('./gen2/same_clumps/VELA13_normalized_ex_situ.out');
norm_es_13=load('./gen2/same_clumps/VELA14_normalized_ex_situ.out');
norm_es_14=load('./gen2/same_clumps/VELA15_normalized_ex_situ.out');
norm_es_15=load('./gen2/same_clumps/VELA16_normalized_ex_situ.out');
norm_es_16=load('./gen2/same_clumps/VELA19_normalized_ex_situ.out');
norm_es_17=load('./gen2/same_clumps/VELA21_normalized_ex_situ.out');
norm_es_18=load('./gen2/same_clumps/VELA23_normalized_ex_situ.out');
norm_es_19=load('./gen2/same_clumps/VELA25_normalized_ex_situ.out');
norm_es_20=load('./gen2/same_clumps/VELA26_normalized_ex_situ.out');
norm_es_21=load('./gen2/same_clumps/VELA27_normalized_ex_situ.out');
norm_es_22=load('./gen2/same_clumps/VELA28_normalized_ex_situ.out');
norm_es_23=load('./gen2/same_clumps/VELA29_normalized_ex_situ.out');
norm_es_24=load('./gen2/same_clumps/VELA30_normalized_ex_situ.out');
norm_es_25=load('./gen2/same_clumps/VELA31_normalized_ex_situ.out');
norm_es_26=load('./gen2/same_clumps/VELA32_normalized_ex_situ.out');
norm_es_27=load('./gen2/same_clumps/VELA33_normalized_ex_situ.out');
norm_es_28=load('./gen2/same_clumps/VELA34_normalized_ex_situ.out');
norm_es_29=load('./gen2/same_clumps/VELA35_normalized_ex_situ.out');

norm_es2 = [norm_es_1
    norm_es_2
    norm_es_3
    norm_es_4
    norm_es_5
    norm_es_6
    norm_es_7
    norm_es_8
    norm_es_9
    norm_es_10
    norm_es_11
    norm_es_12
    norm_es_13
    norm_es_14
    norm_es_15
    norm_es_16
    norm_es_17
    norm_es_18
    norm_es_19
    norm_es_20
    norm_es_21
    norm_es_22
    norm_es_23
    norm_es_24
    norm_es_25
    norm_es_26
    norm_es_27
    norm_es_28
    norm_es_29];
clear norm_es_1 norm_es_2 norm_es_3 norm_es_4 norm_es_5 norm_es_6 norm_es_7
clear norm_es_8 norm_es_9 norm_es_10 norm_es_11 norm_es_12 norm_es_13
clear norm_es_14 norm_es_15 norm_es_16 norm_es_17 norm_es_18 norm_es_19
clear norm_es_20 norm_es_21 norm_es_22 norm_es_23 norm_es_24 norm_es_25 
clear norm_es_26 norm_es_27 norm_es_28 norm_es_29

norm_bulge_1=load('./gen2/same_clumps/VELA01_normalized_bulge.out');
norm_bulge_2=load('./gen2/same_clumps/VELA02_normalized_bulge.out');
norm_bulge_3=load('./gen2/same_clumps/VELA03_normalized_bulge.out');
norm_bulge_4=load('./gen2/same_clumps/VELA05_normalized_bulge.out');
norm_bulge_5=load('./gen2/same_clumps/VELA06_normalized_bulge.out');
norm_bulge_6=load('./gen2/same_clumps/VELA07_normalized_bulge.out');
norm_bulge_7=load('./gen2/same_clumps/VELA08_normalized_bulge.out');
norm_bulge_8=load('./gen2/same_clumps/VELA09_normalized_bulge.out');
norm_bulge_9=load('./gen2/same_clumps/VELA10_normalized_bulge.out');
norm_bulge_10=load('./gen2/same_clumps/VELA11_normalized_bulge.out');
norm_bulge_11=load('./gen2/same_clumps/VELA12_normalized_bulge.out');
norm_bulge_12=load('./gen2/same_clumps/VELA13_normalized_bulge.out');
norm_bulge_13=load('./gen2/same_clumps/VELA14_normalized_bulge.out');
norm_bulge_14=load('./gen2/same_clumps/VELA15_normalized_bulge.out');
norm_bulge_15=load('./gen2/same_clumps/VELA16_normalized_bulge.out');
norm_bulge_16=load('./gen2/same_clumps/VELA19_normalized_bulge.out');
norm_bulge_17=load('./gen2/same_clumps/VELA21_normalized_bulge.out');
norm_bulge_18=load('./gen2/same_clumps/VELA23_normalized_bulge.out');
norm_bulge_19=load('./gen2/same_clumps/VELA25_normalized_bulge.out');
norm_bulge_20=load('./gen2/same_clumps/VELA26_normalized_bulge.out');
norm_bulge_21=load('./gen2/same_clumps/VELA27_normalized_bulge.out');
norm_bulge_22=load('./gen2/same_clumps/VELA28_normalized_bulge.out');
norm_bulge_23=load('./gen2/same_clumps/VELA29_normalized_bulge.out');
norm_bulge_24=load('./gen2/same_clumps/VELA30_normalized_bulge.out');
norm_bulge_25=load('./gen2/same_clumps/VELA31_normalized_bulge.out');
norm_bulge_26=load('./gen2/same_clumps/VELA32_normalized_bulge.out');
norm_bulge_27=load('./gen2/same_clumps/VELA33_normalized_bulge.out');
norm_bulge_28=load('./gen2/same_clumps/VELA34_normalized_bulge.out');
norm_bulge_29=load('./gen2/same_clumps/VELA35_normalized_bulge.out');

norm_bulge2 = [norm_bulge_1
    norm_bulge_2
    norm_bulge_3
    norm_bulge_4
    norm_bulge_5
    norm_bulge_6
    norm_bulge_7
    norm_bulge_8
    norm_bulge_9
    norm_bulge_10
    norm_bulge_11
    norm_bulge_12
    norm_bulge_13
    norm_bulge_14
    norm_bulge_15
    norm_bulge_16
    norm_bulge_17
    norm_bulge_18
    norm_bulge_19
    norm_bulge_20
    norm_bulge_21
    norm_bulge_22
    norm_bulge_23
    norm_bulge_24
    norm_bulge_25
    norm_bulge_26
    norm_bulge_27
    norm_bulge_28
    norm_bulge_29];
clear norm_bulge_1 norm_bulge_2 norm_bulge_3 norm_bulge_4 norm_bulge_5 norm_bulge_6 norm_bulge_7
clear norm_bulge_8 norm_bulge_9 norm_bulge_10 norm_bulge_11 norm_bulge_12 norm_bulge_13
clear norm_bulge_14 norm_bulge_15 norm_bulge_16 norm_bulge_17 norm_bulge_18 norm_bulge_19
clear norm_bulge_20 norm_bulge_21 norm_bulge_22 norm_bulge_23 norm_bulge_24 norm_bulge_25 
clear norm_bulge_26 norm_bulge_27 norm_bulge_28 norm_bulge_29
clear nt

disc2 = load('./gen2/disc_sizes_and_clumps_HcSp.out');
ndisc2 = zeros(1,length(nis2));
j = 2;
for i=2:max(size(disc2))
    if(disc2(i,1)>disc2(i-1,1))
        ndisc2(j) = i-1;
        j = j+1;
    end
end
ndisc2(length(ndisc2)) = max(size(disc2));

b = find(isnan(disc2)); % Will happen if disc has SFR=0
disc2(b) = 1e-6;
b = find(isinf(disc2)); % Will happen if disc has SFR=0
disc2(b) = 1e-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = find(is2(:,5)<1e-6); % No stars
is2(b,5)  = 1e-6; % correct stellar mass
is2(b,10) = 1e-6; % correct stellar surface density
is2(b,12) = 1e-6; % correct stellar age
is2(b,14) = 1e-6; % correct stellar metalicity
is2(b,15) = 1e-6; % correct SFR
is2(b,16) = 1e-6; % correct Sigma SFR
is2(b,17) = 1e-6; % correct sSFR
is2(b,18) = 1e-6; % correct tdep
norm_is2(b,5)  = 1e-6; % correct stellar mass
norm_is2(b,10) = 1e-6; % correct stellar surface density
norm_is2(b,12) = 1e-6; % correct stellar age
norm_is2(b,14) = 1e-6; % correct stellar metalicity
norm_is2(b,15) = 1e-6; % correct SFR
norm_is2(b,16) = 1e-6; % correct Sigma SFR
norm_is2(b,17) = 1e-6; % correct sSFR
norm_is2(b,18) = 1e-6; % correct tdep

b = find(is2(:,15)<1e-6); % No SFR
is2(b,15) = 1e-6; % correct SFR
is2(b,16) = 1e-6; % correct Sigma SFR
is2(b,17) = 1e-6; % correct sSFR
is2(b,18) = 1e-6; % correct tdep
norm_is2(b,15) = 1e-6; % correct SFR
norm_is2(b,16) = 1e-6; % correct Sigma SFR
norm_is2(b,17) = 1e-6; % correct sSFR
norm_is2(b,18) = 1e-6; % correct tdep

b = find(is2(:,4)<1e-6); % No gas
is2(b,4)  = 1e-6; % correct gas mass
is2(b,7)  = 1e-6; % correct gas fraction
is2(b,13) = 1e-6; % correct gas metalicity
is2(b,18) = 1e-6; % correct tdep
norm_is2(b,4)  = 1e-6; % correct gas mass
norm_is2(b,7)  = 1e-6; % correct gas fraction
norm_is2(b,13) = 1e-6; % correct gas metalicity
norm_is2(b,18) = 1e-6; % correct tdep

b = find(isnan(is2(:,30))); % Bad local t_dyn
is2(b,30) = 1e10;
norm_is2(b,30) = 1e10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = find(es2(:,5)<1e-6); % No stars
es2(b,5)  = 1e-6; % correct stellar mass
es2(b,10) = 1e-6; % correct stellar surface density
es2(b,12) = 1e-6; % correct stellar age
es2(b,14) = 1e-6; % correct stellar metalicity
es2(b,15) = 1e-6; % correct SFR
es2(b,16) = 1e-6; % correct Sigma SFR
es2(b,17) = 1e-6; % correct sSFR
es2(b,18) = 1e-6; % correct tdep
norm_es2(b,5)  = 1e-6; % correct stellar mass
norm_es2(b,10) = 1e-6; % correct stellar surface density
norm_es2(b,12) = 1e-6; % correct stellar age
norm_es2(b,14) = 1e-6; % correct stellar metalicity
norm_es2(b,15) = 1e-6; % correct SFR
norm_es2(b,16) = 1e-6; % correct Sigma SFR
norm_es2(b,17) = 1e-6; % correct sSFR
norm_es2(b,18) = 1e-6; % correct tdep

b = find(es2(:,15)<1e-6); % No SFR
es2(b,15) = 1e-6; % correct SFR
es2(b,16) = 1e-6; % correct Sigma SFR
es2(b,17) = 1e-6; % correct sSFR
es2(b,18) = 1e-6; % correct tdep
norm_es2(b,15) = 1e-6; % correct SFR
norm_es2(b,16) = 1e-6; % correct Sigma SFR
norm_es2(b,17) = 1e-6; % correct sSFR
norm_es2(b,18) = 1e-6; % correct tdep

b = find(es2(:,4)<1e-6); % No gas
es2(b,4)  = 1e-6; % correct gas mass
es2(b,7)  = 1e-6; % correct gas fraction
es2(b,13) = 1e-6; % correct gas metalicity
es2(b,18) = 1e-6; % correct tdep
norm_es2(b,4)  = 1e-6; % correct gas mass
norm_es2(b,7)  = 1e-6; % correct gas fraction
norm_es2(b,13) = 1e-6; % correct gas metalicity
norm_es2(b,18) = 1e-6; % correct tdep

b = find(isnan(es2(:,30))); % Bad local t_dyn
es2(b,30) = 1e10;
norm_es2(b,30) = 1e10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = find(bulge2(:,5)<1e-6); % No stars
bulge2(b,5)  = 1e-6; % correct stellar mass
bulge2(b,10) = 1e-6; % correct stellar surface density
bulge2(b,12) = 1e-6; % correct stellar age
bulge2(b,14) = 1e-6; % correct stellar metalicity
bulge2(b,15) = 1e-6; % correct SFR
bulge2(b,16) = 1e-6; % correct Sigma SFR
bulge2(b,17) = 1e-6; % correct sSFR
bulge2(b,18) = 1e-6; % correct tdep
norm_bulge2(b,5)  = 1e-6; % correct stellar mass
norm_bulge2(b,10) = 1e-6; % correct stellar surface density
norm_bulge2(b,12) = 1e-6; % correct stellar age
norm_bulge2(b,14) = 1e-6; % correct stellar metalicity
norm_bulge2(b,15) = 1e-6; % correct SFR
norm_bulge2(b,16) = 1e-6; % correct Sigma SFR
norm_bulge2(b,17) = 1e-6; % correct sSFR
norm_bulge2(b,18) = 1e-6; % correct tdep

b = find(bulge2(:,15)<1e-6); % No SFR
bulge2(b,15) = 1e-6; % correct SFR
bulge2(b,16) = 1e-6; % correct Sigma SFR
bulge2(b,17) = 1e-6; % correct sSFR
bulge2(b,18) = 1e-6; % correct tdep
norm_bulge2(b,15) = 1e-6; % correct SFR
norm_bulge2(b,16) = 1e-6; % correct Sigma SFR
norm_bulge2(b,17) = 1e-6; % correct sSFR
norm_bulge2(b,18) = 1e-6; % correct tdep

b = find(bulge2(:,4)<1e-6); % No gas
bulge2(b,4)  = 1e-6; % correct gas mass
bulge2(b,7)  = 1e-6; % correct gas fraction
bulge2(b,13) = 1e-6; % correct gas metalicity
bulge2(b,18) = 1e-6; % correct tdep
norm_bulge2(b,4)  = 1e-6; % correct gas mass
norm_bulge2(b,7)  = 1e-6; % correct gas fraction
norm_bulge2(b,13) = 1e-6; % correct gas metalicity
norm_bulge2(b,18) = 1e-6; % correct tdep

b = find(isnan(bulge2(:,30))); % Bad local t_dyn
bulge2(b,30) = 1e10;
norm_bulge2(b,30) = 1e10;

clear b i j

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalized radius
norm_rad_is2 = is2(:,21)./is2(:,19); %kpc
norm_rad_es2 = es2(:,21)./es2(:,19); %kpc
norm_rad_bulge2 = bulge2(:,21)./bulge2(:,19); %kpc
norm_is2(:,3) = norm_is2(:,3)./norm_rad_is2;
norm_es2(:,3) = norm_es2(:,3)./norm_rad_es2;
norm_bulge2(:,3) = norm_bulge2(:,3)./norm_rad_bulge2;

clear norm_rad_is2 norm_rad_es2 norm_rad_bulge2

% depletion time
is2(:,18) = 1./is2(:,18);
es2(:,18) = 1./es2(:,18);
bulge2(:,18) = 1./bulge2(:,18);
norm_is2(:,18) = 1./norm_is2(:,18);
norm_es2(:,18) = 1./norm_es2(:,18);
norm_bulge2(:,18) = 1./norm_bulge2(:,18);

% absolute value height
is2(:,20) = abs(is2(:,20));
es2(:,20) = abs(es2(:,20));
bulge2(:,20) = abs(bulge2(:,20));
norm_is2(:,20) = abs(norm_is2(:,20));
norm_es2(:,20) = abs(norm_es2(:,20));
norm_bulge2(:,20) = abs(norm_bulge2(:,20));

% Correct free fall time
is2(:,29) = is2(:,29).*sqrt(3*pi/32);
es2(:,29) = es2(:,29).*sqrt(3*pi/32);
bulge2(:,29) = bulge2(:,29).*sqrt(3*pi/32);
norm_is2(:,29) = norm_is2(:,29).*sqrt(3*pi/32);
norm_es2(:,29) = norm_es2(:,29).*sqrt(3*pi/32);
norm_bulge2(:,29) = norm_bulge2(:,29).*sqrt(3*pi/32);

% absolute value dynamical times
is2(:,30) = abs(is2(:,30));
es2(:,30) = abs(es2(:,30));
bulge2(:,30) = abs(bulge2(:,30));
norm_is2(:,30) = abs(norm_is2(:,30));
norm_es2(:,30) = abs(norm_es2(:,30));
norm_bulge2(:,30) = abs(norm_bulge2(:,30));

is2(:,31) = abs(is2(:,31));
es2(:,31) = abs(es2(:,31));
bulge2(:,31) = abs(bulge2(:,31));
norm_is2(:,31) = abs(norm_is2(:,31));
norm_es2(:,31) = abs(norm_es2(:,31));
norm_bulge2(:,31) = abs(norm_bulge2(:,31));

% lifetime in Myr, normalized by free fall time
is2(:,41) = 1000.*is2(:,41);
es2(:,41) = 1000.*es2(:,41);
bulge2(:,41) = 1000.*bulge2(:,41);
norm_is2(:,41) = is2(:,41)./is2(:,29);
norm_es2(:,41) = es2(:,41)./es2(:,29);
norm_bulge2(:,41) = bulge2(:,41)./bulge2(:,29);
[max(norm_is2(:,41)),min(norm_is2(:,41)),length(find(isnan(norm_is2(:,41))))]
