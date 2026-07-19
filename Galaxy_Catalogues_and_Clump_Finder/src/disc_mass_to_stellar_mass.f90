!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program main
	implicit none
	character(len=20),allocatable :: gal_name(:)
	character(len=20) :: gen_string
	character(len=10000) :: filename, input_arg
	integer :: i, j, k, Nsimulation, Nsnapshot1, Nsnapshot2, Nsnapshot_tot
	real(4) :: aexp, rcom(3), vcom(3), Rdisc, spher_mgas, spher_mcold, spher_mstars, spher_exsitu_stellar_mass, spher_mdm, spher_SFR, spher_age, spher_metgas, spher_metstars
	real(4) :: aexp2, rcom2(3), vcom2(3), ang_mom(3), Lmag, Rdisc2, Hdisc, mean_mgas, mean_mcold, mean_mstars, exsitu_stellar_mass, mean_mdm
	real(4) :: mean_SFR, mean_age, mean_metgas, mean_metstars
	real(4) :: gen, redshift


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

	write(filename,'(a,a,a)') './',trim(gen_string),'/disc_mass_to_stellar_mass.txt'
	open(unit=12,file=filename,form='formatted')

	Nsnapshot_tot = 0
	do i=1,Nsimulation
		write(filename,'(a,a,a)') '/BIGDATA/nirm/galaxy_inputs/',trim(gal_name(i)),'/galaxy_catalogue/Nir_spherical_galaxy_cat.txt'
		open(unit=13,file=filename,form='formatted')
		read(13,'(1x,i3)') Nsnapshot1
		write(filename,'(a,a,a)') '/BIGDATA/nirm/galaxy_inputs/',trim(gal_name(i)),'/galaxy_catalogue/Nir_disc_cat.txt'
		open(unit=14,file=filename,form='formatted')
		read(14,'(1x,i3)') Nsnapshot2

		if(Nsnapshot1 .ne. Nsnapshot2) then
			print *, 'Error: Problem with number of snapshots'
			print *, 'gal name=',gal_name(i)
			print *, 'Nsnapshot1=',Nsnapshot1,'Nsnapshot2=',Nsnapshot2
			stop
		end if

		do j=1,Nsnapshot1
			read(13,'(1x,f5.3,6(1x,es12.5),1x,f7.3,7(1x,es12.5),2(1x,f7.3))') aexp, rcom(:), vcom(:), Rdisc, spher_mgas, spher_mcold, spher_mstars, spher_exsitu_stellar_mass, spher_mdm, spher_SFR, spher_age, spher_metgas, spher_metstars
			read(14,'(1x,f5.3,10(1x,es12.5),2(1x,f7.3),7(1x,es12.5),2(1x,f7.3))') aexp2, rcom2(:), vcom2(:), ang_mom(:), Lmag, Rdisc2, Hdisc, mean_mgas, mean_mcold, mean_mstars, exsitu_stellar_mass, mean_mdm, mean_SFR, mean_age, mean_metgas, mean_metstars

			if(aexp .ne. aexp2) then
				print *, 'Error: Problem with expansion factor'
				print *, 'gal name=',gal_name(i)
				print *, 'aexp=',aexp,'aexp2=',aexp2
				stop
			end if

			if(rcom(1) .ne. rcom2(1) .or. rcom(2) .ne. rcom2(2) .or. rcom(3) .ne. rcom2(3) ) then
				print *, 'Error: Problem with COM'
				print *, 'gal name=',gal_name(i)
				print *, 'rcom=',rcom(:),'rcom2=',rcom2(:)
				stop
			end if

			if(vcom(1) .ne. vcom2(1) .or. vcom(2) .ne. vcom2(2) .or. vcom(3) .ne. vcom2(3) ) then
				print *, 'Error: Problem with COM velocity'
				print *, 'gal name=',gal_name(i)
				print *, 'vcom=',vcom(:),'vcom2=',vcom2(:)
				stop
			end if

			if(Rdisc .ne. Rdisc2) then
				print *, 'Error: Problem with radius'
				print *, 'gal name=',gal_name(i)
				print *, 'Rd=',Rdisc,'Rd2=',Rdisc2
				stop
			end if

			redshift = 1.0_4 / (aexp) - 1.0_4
			write(12,'(1x,a,1x,f5.3,2(1x,f7.3),4(1x,es12.5))') trim(gal_name(i)),redshift, Rdisc, Hdisc, mean_mgas, mean_mstars, spher_mgas, spher_mstars
			Nsnapshot_tot = Nsnapshot_tot + 1
		end do
		close(unit=13)
		close(unit=14)
	end do
	print *, Nsnapshot_tot
	close(unit=12)
end program main



