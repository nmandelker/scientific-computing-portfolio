program main
	implicit none
	integer,parameter :: Neps_energy = 2
	integer :: i, Nsnapshot
	real(4),allocatable :: tsnap(:)
	real(4) :: tsnap2
	real(8) :: stream_momentum(24), Hb(2), Hs(2), stream_radius(2), clumping(Neps_energy+1), pure_frac(Neps_energy), Vc(3)
	real(8) :: sigma_gauss, stream_masses(4), stream_volumes(4), turb_vel(6), turb_smooth_volume(20), turb_smooth_mass(20)
	real(8) :: clumping_dense_cold(6), energies(9), energy_fluxes(6), luminosity_vec(2)
	character(len=256) :: filename, format_string, input_arg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(iargc().ge.1) then
		call getarg(1,input_arg)
		read(input_arg,*) sigma_gauss
	else
		sigma_gauss = 0.25
	end if
	print *, 'sigma_gauss'
	print *, sigma_gauss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	open(unit=16,file='./output/time.txt')
	read(16,*) Nsnapshot
	print *, Nsnapshot
	allocate( tsnap(Nsnapshot) )
	do i=1,Nsnapshot
		read(16,'(F6.4)') tsnap(i)
	end do
	close(unit=16)

	stream_momentum(:)     = 0.0_8
	Hb(:)                  = 0.0_8
	Hs(:)                  = 0.0_8
	stream_radius(:)       = 0.0_8
	clumping(:)            = 0.0_8
	pure_frac(:)           = 0.0_8
	Vc(:)                  = 0.0_8
	stream_masses(:)       = 0.0_8
	turb_vel(:)            = 0.0_8
	clumping_dense_cold(:) = 0.0_8
	energies(:)            = 0.0_8
	energy_fluxes(:)       = 0.0_8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	write(filename,'(a)') './stream_analysis/momentum.txt'
	open(unit=20,file=filename,form='formatted')
	filename = ''

	write(filename,'(a)') './stream_analysis/thickness.txt'
	open(unit=21,file=filename,form='formatted')
	filename = ''

	write(filename,'(a)') './stream_analysis/clumping.txt'
	open(unit=22,file=filename,form='formatted')
	filename = ''

	write(filename,'(a)') './stream_analysis/pure_frac.txt'
	open(unit=23,file=filename,form='formatted')
	filename = ''

	write(filename,'(a)') './stream_analysis/Vc.txt'
	open(unit=24,file=filename,form='formatted')
	filename = ''

	write(filename,'(a)') './stream_analysis/stream_masses.txt'
	open(unit=25,file=filename,form='formatted')
	filename = ''

	write(filename,'(a)') './stream_analysis/turbulence_from_profile.txt'
	open(unit=26,file=filename,form='formatted')
	filename = ''

	write(filename,'(a)') './stream_analysis/clumping_factor_dense_cold.txt'
	open(unit=27,file=filename,form='formatted')
	filename = ''

	write(filename,'(a)') './stream_analysis/energies_from_profile.txt'
	open(unit=28,file=filename,form='formatted')
	filename = ''

	write(filename,'(a)') './stream_analysis/energy_fluxes.txt'
	open(unit=29,file=filename,form='formatted')
	filename = ''

	write(filename,'(a)') './stream_analysis/stream_volumes.txt'
	open(unit=43,file=filename,form='formatted')
	filename = ''

	write(filename,'(a)') './stream_analysis/turbulence_from_smoothing_volume_weighted.txt'
	open(unit=44,file=filename,form='formatted')
	filename = ''

	write(filename,'(a)') './stream_analysis/turbulence_from_smoothing_mass_weighted.txt'
	open(unit=45,file=filename,form='formatted')
	filename = ''

	write(filename,'(a)') './stream_analysis/Luminosity.txt'
	open(unit=46,file=filename,form='formatted')
	filename = ''
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	do i=1,Nsnapshot
		print *, ''
		print *, 'i, Nsnapshot=',i,Nsnapshot
		print *, 't=',tsnap(i)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		write(filename,'(a,I5.5,a)') './stream_analysis/momentum/momentum_',i,'.txt'
		open(unit=30,file=filename,form='formatted')
		read(30,'(F6.4,24(1x,ES12.5))') tsnap2, stream_momentum(1:22)
		if(tsnap2 .ne. tsnap(i)) then
			print *, 'ERROR with tsnap in momentum:', i, tsnap(i), tsnap2
		end if
		close(unit=30)
		filename = ''
		write(20,'(F6.4,24(1x,ES12.5))') tsnap(i), stream_momentum(1:24)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		write(filename,'(a,I5.5,a)') './stream_analysis/thickness/thickness_',i,'.txt'
		open(unit=31,file=filename,form='formatted')
		read(31,'(F6.4,6(1x,ES12.5))') tsnap2, Hs(1), Hb(1), stream_radius(1), Hs(2), Hb(2), stream_radius(2)
		if(tsnap2 .ne. tsnap(i)) then
			print *, 'ERROR with tsnap in thickness:', i, tsnap(i), tsnap2
		end if
		close(unit=31)
		filename = ''
		write(21,'(F6.4,6(1x,ES12.5))') tsnap(i), Hs(1), Hb(1), stream_radius(1), Hs(2), Hb(2), stream_radius(2)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		write(filename,'(a,I5.5,a)') './stream_analysis/clumping/clumping_',i,'.txt'
		open(unit=32,file=filename,form='formatted')
		write(format_string,'(a,i1,a)') '(F6.4,',Neps_energy+1,'(1x,ES12.5))'
		read(32,trim(format_string)) tsnap2, clumping(1:(Neps_energy+1))
		if(tsnap2 .ne. tsnap(i)) then
			print *, 'ERROR with tsnap in clumping:', i, tsnap(i), tsnap2
		end if
		close(unit=32)
		filename = ''
		write(22,trim(format_string)) tsnap(i), clumping(1:(Neps_energy+1))
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		write(filename,'(a,I5.5,a)') './stream_analysis/pure_frac/pure_frac_',i,'.txt'
		open(unit=33,file=filename,form='formatted')
		write(format_string,'(a,i1,a)') '(F6.4,',Neps_energy,'(1x,ES12.5))'
		read(33,trim(format_string)) tsnap2, pure_frac(1:Neps_energy)
		if(tsnap2 .ne. tsnap(i)) then
			print *, 'ERROR with tsnap in pure_frac:', i, tsnap(i), tsnap2
		end if
		close(unit=33)
		filename = ''
		write(23,trim(format_string)) tsnap(i), pure_frac(1:Neps_energy)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		write(filename,'(a,I5.5,a)') './stream_analysis/Vc/Vc_',i,'.txt'
		open(unit=34,file=filename,form='formatted')
		read(34,'(F6.4,3(1x,ES12.5))') tsnap2, Vc(1:3)
		if(tsnap2 .ne. tsnap(i)) then
			print *, 'ERROR with tsnap in Vc:', i, tsnap(i), tsnap2
		end if
		close(unit=34)
		filename = ''
		write(24,'(F6.4,3(1x,ES12.5))') tsnap(i), Vc(1:3)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		write(filename,'(a,I5.5,a)') './stream_analysis/stream_masses/stream_masses_',i,'.txt'
		open(unit=35,file=filename,form='formatted')
		read(35,'(F6.4,4(1x,ES12.5))') tsnap2, stream_masses(1:4)
		if(tsnap2 .ne. tsnap(i)) then
			print *, 'ERROR with tsnap in stream_masses:', i, tsnap(i), tsnap2
		end if
		close(unit=35)
		filename = ''
		write(25,'(F6.4,4(1x,ES12.5))') tsnap(i), stream_masses(1:4)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		write(filename,'(a,I5.5,a)') './stream_analysis/turbulence_from_profile/turbulence_from_profile_',i,'.txt'
		open(unit=36,file=filename,form='formatted')
		read(36,'(F6.4,6(1x,ES12.5))') tsnap2, turb_vel(1:6)
		if(tsnap2 .ne. tsnap(i)) then
			print *, 'ERROR with tsnap in turbulence_from_profile:', i, tsnap(i), tsnap2
		end if
		close(unit=36)
		filename = ''
		write(26,'(F6.4,6(1x,ES12.5))') tsnap(i), turb_vel(1:6)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		write(filename,'(a,I5.5,a)') './stream_analysis/clumping_factor_dense_cold/clumping_factor_dense_cold_',i,'.txt'
		open(unit=37,file=filename,form='formatted')
		read(37,'(F6.4,6(1x,ES12.5))') tsnap2, clumping_dense_cold(1:6)
		if(tsnap2 .ne. tsnap(i)) then
			print *, 'ERROR with tsnap in clumping_factor_dense_cold:', i, tsnap(i), tsnap2
		end if
		close(unit=37)
		filename = ''
		write(27,'(F6.4,6(1x,ES12.5))') tsnap(i), clumping_dense_cold(1:6)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		write(filename,'(a,I5.5,a)') './stream_analysis/energies_from_profile/energies_from_profile_',i,'.txt'
		open(unit=38,file=filename,form='formatted')
		read(38,'(F6.4,9(1x,ES12.5))') tsnap2, energies(1:9)
		if(tsnap2 .ne. tsnap(i)) then
			print *, 'ERROR with tsnap in energies_from_profile_:', i, tsnap(i), tsnap2
		end if
		close(unit=38)
		filename = ''
		write(28,'(F6.4,9(1x,ES12.5))') tsnap(i), energies(1:9)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		write(filename,'(a,I5.5,a)') './stream_analysis/energy_fluxes/energy_fluxes_',i,'.txt'
		open(unit=39,file=filename,form='formatted')
		read(39,'(F6.4,6(1x,ES12.5))') tsnap2, energy_fluxes(1:6)
		if(tsnap2 .ne. tsnap(i)) then
			print *, 'ERROR with tsnap in energy_fluxes:', i, tsnap(i), tsnap2
		end if
		close(unit=39)
		filename = ''
		write(29,'(F6.4,6(1x,ES12.5))') tsnap(i), energy_fluxes(1:6)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		write(filename,'(a,I5.5,a)') './stream_analysis/stream_volumes/stream_volumes_',i,'.txt'
		open(unit=53,file=filename,form='formatted')
		read(53,'(F6.4,4(1x,ES12.5))') tsnap2, stream_volumes(1:4)
		if(tsnap2 .ne. tsnap(i)) then
			print *, 'ERROR with tsnap in stream_volmes:', i, tsnap(i), tsnap2
		end if
		close(unit=53)
		filename = ''
		write(43,'(F6.4,4(1x,ES12.5))') tsnap(i), stream_volumes(1:4)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(i.le.300) then
			write(filename,'(a,f4.2,a,I5.5,a)') './stream_analysis/turbulence_from_smoothing/sigma_',sigma_gauss,'/volume_weighted_',i,'.txt'
			open(unit=54,file=filename,form='formatted')
			read(54,'(F6.4,20(1x,ES12.5))') tsnap2, turb_smooth_volume(1:20)
			if(tsnap2 .ne. tsnap(i)) then
				print *, 'ERROR with tsnap in volume weighted smoothed turbulence:', i, tsnap(i), tsnap2
			end if
			close(unit=54)
			filename = ''
			write(44,'(F6.4,20(1x,ES12.5))') tsnap(i), turb_smooth_volume(1:20)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			write(filename,'(a,f4.2,a,I5.5,a)') './stream_analysis/turbulence_from_smoothing/sigma_',sigma_gauss,'/mass_weighted_',i,'.txt'
			open(unit=55,file=filename,form='formatted')
			read(55,'(F6.4,20(1x,ES12.5))') tsnap2, turb_smooth_mass(1:20)
			if(tsnap2 .ne. tsnap(i)) then
				print *, 'ERROR with tsnap in mass weighted smoothed turbulence:', i, tsnap(i), tsnap2
			end if
			close(unit=55)
			filename = ''
			write(45,'(F6.4,20(1x,ES12.5))') tsnap(i), turb_smooth_mass(1:20)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			write(filename,'(a,I5.5,a)') './stream_analysis/Luminosity/Luminosity',i,'.txt'
			open(unit=56,file=filename,form='formatted')
			read(56,'(F6.4,2(1x,ES12.5))') tsnap2, luminosity_vec(1:2)
			if(tsnap2 .ne. tsnap(i)) then
				print *, 'ERROR with tsnap in Luminosity:', i, tsnap(i), tsnap2
			end if
			close(unit=56)
			filename = ''
			write(46,'(F6.4,2(1x,ES12.5))') tsnap(i), luminosity_vec(1:2)
		end if

	end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	deallocate( tsnap )
	close(unit=20)
	close(unit=21)
	close(unit=22)
	close(unit=23)
	close(unit=24)
	close(unit=25)
	close(unit=26)
	close(unit=27)
	close(unit=28)
	close(unit=29)
	close(unit=43)
	close(unit=44)
	close(unit=45)
	close(unit=46)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program main

