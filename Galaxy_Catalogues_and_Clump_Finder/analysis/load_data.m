clear

[is2, es2, bulge2, nis2, nes2, nbulge2, norm_is2, norm_es2, norm_bulge2, disc2, ndisc2] = load_gen2();
[is3, es3, bulge3, nis3, nes3, nbulge3, norm_is3, norm_es3, norm_bulge3, disc3, ndisc3] = load_gen3();

gen2_vs_gen3_galaxy_discs_spheres;

post_referee_properties;

add_properties;

[is2_common, es2_common, bulge2_common, disc2_common, nis2_common, nes2_common, ...
    nbulge2_common, ndisc2_common, is3_common, es3_common, bulge3_common, ...
    disc3_common, nis3_common, nes3_common, nbulge3_common, ndisc3_common, ...
    ind_common, amax_common, amin_common] = ...
common_sample_gen2_3(is2, es2, bulge2, disc2, nis2, nes2, nbulge2, ndisc2, ...
    is3, es3, bulge3, disc3, nis3, nes3, nbulge3, ndisc3);

nis2_common_and_friends;

load_Behroozi;