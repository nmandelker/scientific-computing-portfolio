program main
	implicit none
	character(len=20),allocatable :: gal_name(:)
	character(len=20) :: gen_string
	character(len=256) :: input_arg
	character(len=10000) :: filename1, filename2, filename3, filename4, filename5, filename6, filename7, filename8, filename9, filename10
	integer :: i, j, k, Nsnapshot, Nsimulation, nbinsitu, nbexsitu, nmed
	real(4) :: aexp, rcom(3), vcom(3), Ldisc(3), Rdisc, Hdisc, Mgas_disc, Mstar_disc, M_Es_star_disc, SFR_disc, age_disc, metgas_disc, metstars_disc, gen

	if(iargc().ge.1) then
		call getarg(1,input_arg)
		read(input_arg(1:3),*) gen
	else
		gen = 2.0
	endif
	print *, gen

	if(gen .eq. 1.0) then
		write(gen_string,'(a)') 'gen1'
		Nsimulation = 30
		allocate( gal_name(Nsimulation) )
		gal_name(1) = 'MW1'
		gal_name(2) = 'MW2'
		gal_name(3) = 'MW3'
		gal_name(4) = 'MW4'
		gal_name(5) = 'MW5'
		gal_name(6) = 'MW6'
		gal_name(7) = 'MW7'
		gal_name(8) = 'MW8'
		gal_name(9) = 'MW9'
		gal_name(10) = 'MW10'
		gal_name(11) = 'MW11'
		gal_name(12) = 'VL01'
		gal_name(13) = 'VL02'
		gal_name(14) = 'VL03'
		gal_name(15) = 'VL04'
		gal_name(16) = 'VL05'
		gal_name(17) = 'VL06'
		gal_name(18) = 'VL07'
		gal_name(19) = 'VL08'
		gal_name(20) = 'VL09'
		gal_name(21) = 'VL10'
		gal_name(22) = 'VL11'
		gal_name(23) = 'VL12'
		gal_name(24) = 'SFG1'
		gal_name(25) = 'SFG2'
		gal_name(26) = 'SFG4'
		gal_name(27) = 'SFG5'
		gal_name(28) = 'SFG7'
		gal_name(29) = 'SFG8'
		gal_name(30) = 'SFG9'

	elseif(gen .eq. 1.5) then
		write(gen_string,'(a)') 'gen15'
		Nsimulation = 3
		allocate( gal_name(Nsimulation) )
		gal_name(1) = 'SFG1_SF_NEW'
		gal_name(2) = 'SFG4_SF_NEW'
		gal_name(3) = 'SFG5_SF_NEW'

	elseif(gen .eq. 2.0) then
		write(gen_string,'(a)') 'gen2'
		Nsimulation = 29
		allocate( gal_name(Nsimulation) )
		gal_name(1) = 'VELA01'
		gal_name(2) = 'VELA02'
		gal_name(3) = 'VELA03'
		gal_name(4) = 'VELA05'
		gal_name(5) = 'VELA06'
		gal_name(6) = 'VELA07'
		gal_name(7) = 'VELA08'
		gal_name(8) = 'VELA09'
		gal_name(9) = 'VELA10'
		gal_name(10) = 'VELA11'
		gal_name(11) = 'VELA12'
		gal_name(12) = 'VELA13'
		gal_name(13) = 'VELA14'
		gal_name(14) = 'VELA15'
		gal_name(15) = 'VELA16'
		gal_name(16) = 'VELA19'
		gal_name(17) = 'VELA21'
		gal_name(18) = 'VELA23'
		gal_name(19) = 'VELA25'
		gal_name(20) = 'VELA26'
		gal_name(21) = 'VELA27'
		gal_name(22) = 'VELA28'
		gal_name(23) = 'VELA29'
		gal_name(24) = 'VELA30'
		gal_name(25) = 'VELA31'
		gal_name(26) = 'VELA32'
		gal_name(27) = 'VELA33'
		gal_name(28) = 'VELA34'
		gal_name(29) = 'VELA35'

	elseif(gen .eq. 3.0) then
		write(gen_string,'(a)') 'gen3'
		Nsimulation = 34
		allocate( gal_name(Nsimulation) )
		gal_name(1) = 'VELA_v2_01'
		gal_name(2) = 'VELA_v2_02'
		gal_name(3) = 'VELA_v2_03'
		gal_name(4) = 'VELA_v2_04'
		gal_name(5) = 'VELA_v2_05'
		gal_name(6) = 'VELA_v2_06'
		gal_name(7) = 'VELA_v2_07'
		gal_name(8) = 'VELA_v2_08'
		gal_name(9) = 'VELA_v2_09'
		gal_name(10) = 'VELA_v2_10'
		gal_name(11) = 'VELA_v2_11'
		gal_name(12) = 'VELA_v2_12'
		gal_name(13) = 'VELA_v2_13'
		gal_name(14) = 'VELA_v2_14'
		gal_name(15) = 'VELA_v2_15'
		gal_name(16) = 'VELA_v2_16'
		gal_name(17) = 'VELA_v2_17'
		gal_name(18) = 'VELA_v2_19'
		gal_name(19) = 'VELA_v2_20'
		gal_name(20) = 'VELA_v2_21'
		gal_name(21) = 'VELA_v2_22'
		gal_name(22) = 'VELA_v2_23'
		gal_name(23) = 'VELA_v2_24'
		gal_name(24) = 'VELA_v2_25'
		gal_name(25) = 'VELA_v2_26'
		gal_name(26) = 'VELA_v2_27'
		gal_name(27) = 'VELA_v2_28'
		gal_name(28) = 'VELA_v2_29'
		gal_name(29) = 'VELA_v2_30'
		gal_name(30) = 'VELA_v2_31'
		gal_name(31) = 'VELA_v2_32'
		gal_name(32) = 'VELA_v2_33'
		gal_name(33) = 'VELA_v2_34'
		gal_name(34) = 'VELA_v2_35'

	elseif(gen .ne. 1.0 .and. gen .ne. 1.5 .and. gen .ne. 2.0 .and. gen .ne. 3.0 ) then
		print *, 'Error: What generation simulation do you want?'
		print *, 'gen=',gen
		stop
	end if

	write(filename1,'(a,a)') 'mkdir -p ./',trim(gen_string)
	call system(filename1)

	write(filename1,'(a)') 'cat'
	write(filename2,'(a)') 'cat'
	write(filename3,'(a)') 'cat'
	write(filename4,'(a)') 'cat'
	write(filename5,'(a)') 'cat'
	write(filename6,'(a)') 'cat'
	write(filename7,'(a)') 'cat'
	do i=1,Nsimulation
		write(filename1,'(a,a,a,a)') trim(filename1),' ./',trim(gal_name(i)),'_thresh_10_wide_FWHM_2.50/Matlab_friendly_data/bulge.out'
		write(filename2,'(a,a,a,a)') trim(filename2),' ./',trim(gal_name(i)),'_thresh_10_wide_FWHM_2.50/Matlab_friendly_data/in_situ.out'
		write(filename3,'(a,a,a,a)') trim(filename3),' ./',trim(gal_name(i)),'_thresh_10_wide_FWHM_2.50/Matlab_friendly_data/ex_situ.out'
		write(filename4,'(a,a,a,a)') trim(filename4),' ./',trim(gal_name(i)),'_thresh_10_wide_FWHM_2.50/Matlab_friendly_data/normalized_bulge.out'
		write(filename5,'(a,a,a,a)') trim(filename5),' ./',trim(gal_name(i)),'_thresh_10_wide_FWHM_2.50/Matlab_friendly_data/normalized_in_situ.out'
		write(filename6,'(a,a,a,a)') trim(filename6),' ./',trim(gal_name(i)),'_thresh_10_wide_FWHM_2.50/Matlab_friendly_data/normalized_ex_situ.out'
		write(filename7,'(a,a,a,a)') trim(filename7),' ./',trim(gal_name(i)),'_thresh_10_wide_FWHM_2.50/galaxy_catalogue/disc_sizes_and_clumps_HcSp.out'
	end do
	write(filename1,'(a,a,a,a)') trim(filename1),' > ./',trim(gen_string),'/bulge.out'
	write(filename2,'(a,a,a,a)') trim(filename2),' > ./',trim(gen_string),'/in_situ.out'
	write(filename3,'(a,a,a,a)') trim(filename3),' > ./',trim(gen_string),'/ex_situ.out'
	write(filename4,'(a,a,a,a)') trim(filename4),' > ./',trim(gen_string),'/normalized_bulge.out'
	write(filename5,'(a,a,a,a)') trim(filename5),' > ./',trim(gen_string),'/normalized_in_situ.out'
	write(filename6,'(a,a,a,a)') trim(filename6),' > ./',trim(gen_string),'/normalized_ex_situ.out'
	write(filename7,'(a,a,a,a)') trim(filename7),' > ./',trim(gen_string),'/disc_sizes_and_clumps_HcSp.out'
	print*, trim(filename1)
	print*, ''
	print*, trim(filename2)
	print*, ''
	print*, trim(filename3)
	print*, ''
	print*, trim(filename4)
	print*, ''
	print*, trim(filename5)
	print*, ''
	print*, trim(filename6)
	print*, ''
	print*, trim(filename7)
	print*, ''
	call system(filename1)
	call system(filename2)
	call system(filename3)
	call system(filename4)
	call system(filename5)
	call system(filename6)
	call system(filename7)
end program main

