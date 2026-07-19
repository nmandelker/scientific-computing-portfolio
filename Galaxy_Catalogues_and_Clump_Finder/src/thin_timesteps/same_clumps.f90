module globvar
	implicit none
	real(4),allocatable :: is_redshift(:,:), is_rad(:,:), is_gas_mass(:,:), is_star_mass(:,:), is_mass(:,:)
	real(4),allocatable :: is_gas_frac(:,:), is_dm_frac(:,:), is_gas_sig(:,:), is_star_sig(:,:), is_sig(:,:)
	real(4),allocatable :: is_age(:,:), is_gas_met(:,:), is_star_met(:,:), is_SFR(:,:), is_sig_SFR(:,:), is_sSFR(:,:), is_tau(:,:)
	real(4),allocatable :: is_d_Rd(:,:), is_z_Hd(:,:), is_d(:,:), is_z(:,:), is_residual(:,:), is_shape(:,:), is_dm(:,:)
	real(4),allocatable :: is_tff(:,:), is_td(:,:), is_td_global(:,:), is_mgas_in(:,:,:), is_mgas_out(:,:,:) 
	real(4),allocatable :: is_mstars_in(:,:), is_mstars_out(:,:), is_mstars_formed(:,:), is_alpha(:,:)
	integer,allocatable :: is_es(:,:), is_merger(:,:), is_id(:,:), is_list(:)

	real(4),allocatable :: es_redshift(:,:), es_rad(:,:), es_gas_mass(:,:), es_star_mass(:,:), es_mass(:,:)
	real(4),allocatable :: es_gas_frac(:,:), es_dm_frac(:,:), es_gas_sig(:,:), es_star_sig(:,:), es_sig(:,:)
	real(4),allocatable :: es_age(:,:), es_gas_met(:,:), es_star_met(:,:), es_SFR(:,:), es_sig_SFR(:,:), es_sSFR(:,:), es_tau(:,:)
	real(4),allocatable :: es_d_Rd(:,:), es_z_Hd(:,:), es_d(:,:), es_z(:,:), es_residual(:,:), es_shape(:,:), es_dm(:,:)
	real(4),allocatable :: es_tff(:,:), es_td(:,:), es_td_global(:,:), es_mgas_in(:,:,:), es_mgas_out(:,:,:) 
	real(4),allocatable :: es_mstars_in(:,:), es_mstars_out(:,:), es_mstars_formed(:,:), es_alpha(:,:)
	integer,allocatable :: es_es(:,:), es_merger(:,:), es_id(:,:), es_list(:)

	real(4),allocatable :: bulge_redshift(:,:), bulge_rad(:,:), bulge_gas_mass(:,:), bulge_star_mass(:,:), bulge_mass(:,:)
	real(4),allocatable :: bulge_gas_frac(:,:), bulge_dm_frac(:,:), bulge_gas_sig(:,:), bulge_star_sig(:,:), bulge_sig(:,:)
	real(4),allocatable :: bulge_age(:,:), bulge_gas_met(:,:), bulge_star_met(:,:), bulge_SFR(:,:), bulge_sig_SFR(:,:), bulge_sSFR(:,:), bulge_tau(:,:)
	real(4),allocatable :: bulge_d_Rd(:,:), bulge_z_Hd(:,:), bulge_d(:,:), bulge_z(:,:), bulge_residual(:,:), bulge_shape(:,:), bulge_dm(:,:)
	real(4),allocatable :: bulge_tff(:,:), bulge_td(:,:), bulge_td_global(:,:), bulge_mgas_in(:,:,:), bulge_mgas_out(:,:,:) 
	real(4),allocatable :: bulge_mstars_in(:,:), bulge_mstars_out(:,:), bulge_mstars_formed(:,:), bulge_alpha(:,:)
	integer,allocatable :: bulge_es(:,:), bulge_merger(:,:), bulge_id(:,:), bulge_list(:)
end module globvar

program main
use globvar
	implicit none
	integer :: i, j, k, l, m, n, ind, nis, nes, nbulge, Nsimulation, nex_in, nbulge_in, nbulge_ex
	integer,allocatable :: ex_in(:), bulge_in(:), bulge_ex(:)
	real(4) :: time, time2, dt, dens_thresh, wide
	character(len=256) :: filename, dirname, input_arg
	character(len=20),allocatable :: gal_name(:)

!!!!!!!!!! Get smoothing scale and residual threshold !!!!!!!!!!
	if(iargc().ge.1) then
		call getarg(1,input_arg)
		read(input_arg,*) dens_thresh
		if(iargc().ge.2) then
			call getarg(2,input_arg)
			read(input_arg,*) wide
		else
			wide = 2500.0_4
		end if
	else
		dens_thresh = 10.0_4
		wide = 2500.0_4
	end if
	print *, 'dens_thresh, wide_FWHM [pc]'
	print *, dens_thresh, wide

!!!!!!!!!! What simulations will we be looking at today? !!!!!!!!!!
	write(filename,'(a)') './same_clumps_input.dat'
	open(unit=15,file=filename,form='formatted')
	read(15,*) Nsimulation
	print *, Nsimulation
	allocate( gal_name(Nsimulation), stat=i ) 	!!! Deallocated at end of program
	if(i.ne.0) then
		print *, 'error in allocation gal_name. stat= ', i
		stop
	end if
	do i=1,Nsimulation
		read(15,*) gal_name(i)
	end do
	close(unit=15)

	!!!!!!!!!! Loop over all simulations !!!!!!!!!!
	do k=1,Nsimulation
		print *, 'Simulation:', k, trim(gal_name(k))
		write(dirname,'(a,a,a,i2.2,a,f4.2)') './',trim(gal_name(k)),'_thresh_',floor(dens_thresh),'_wide_FWHM_',wide/1000.0_4

		write(filename,'(a,a,a)') 'mkdir -p ./',trim(dirname),'/Matlab_friendly_data/same_clumps'
		call system(trim(filename))

		write(filename,'(a,a)') trim(dirname),'/Matlab_friendly_data/in_situ.out'
		open(unit=12,file=filename,status="old",action="read",form='formatted')
		nis = 0
		do
			read(12,*,end=6)
			nis = nis + 1
		end do
 6		rewind(unit=12)
		if (nis .gt. 0) then
			allocate( is_redshift(nis,2), is_rad(nis,2), is_gas_mass(nis,2), is_star_mass(nis,2), is_mass(nis,2) )
			allocate( is_gas_frac(nis,2), is_dm_frac(nis,2), is_gas_sig(nis,2), is_star_sig(nis,2), is_sig(nis,2) )
			allocate( is_age(nis,2), is_gas_met(nis,2), is_star_met(nis,2), is_SFR(nis,2), is_sig_SFR(nis,2), is_sSFR(nis,2), is_tau(nis,2) )
			allocate( is_d_Rd(nis,2), is_z_Hd(nis,2), is_d(nis,2), is_z(nis,2), is_residual(nis,2), is_shape(nis,2), is_dm(nis,2) )
			allocate( is_tff(nis,2), is_td(nis,2), is_td_global(nis,2), is_mgas_in(nis,2,3), is_mgas_out(nis,2,10) )
			allocate( is_mstars_in(nis,2), is_mstars_out(nis,2), is_mstars_formed(nis,2), is_alpha(nis,2) )
			allocate( is_es(nis,2), is_merger(nis,2), is_id(nis,2), is_list(nis) )
			is_list(:) = 1
			do i = 1, nis
				read(12,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,20(1x,es10.2))') is_id(i,1), is_redshift(i,1), is_rad(i,1), is_gas_mass(i,1), &
				& is_star_mass(i,1), is_mass(i,1), is_gas_frac(i,1), is_dm_frac(i,1), is_gas_sig(i,1), is_star_sig(i,1), is_sig(i,1), is_age(i,1), &
				& is_gas_met(i,1), is_star_met(i,1), is_SFR(i,1), is_sig_SFR(i,1), is_sSFR(i,1), is_tau(i,1), is_d_Rd(i,1), is_z_Hd(i,1), is_d(i,1), is_z(i,1), &
				& is_residual(i,1), is_shape(i,1), is_dm(i,1), is_es(i,1), is_merger(i,1), is_tff(i,1), is_td(i,1), is_td_global(i,1), &
				& is_mgas_in(i,1,1), is_mgas_in(i,1,2), is_mgas_in(i,1,3), &
				& is_mgas_out(i,1,1), is_mgas_out(i,1,2), is_mgas_out(i,1,3), is_mgas_out(i,1,4), is_mgas_out(i,1,5), is_mgas_out(i,1,6), is_mgas_out(i,1,7), & 
				& is_mgas_out(i,1,8), is_mgas_out(i,1,9), is_mgas_out(i,1,10), &
				& is_mstars_in(i,1), is_mstars_out(i,1), is_mstars_formed(i,1), is_alpha(i,1)
			end do
		else
			allocate( is_id(1,2), is_redshift(1,2), is_list(1) )
			is_id(:,:) = 0
			is_redshift(:,:) = 0.0_4
			is_list(:) = 0
		end if
		close(unit=12)
		print *, 'nis=', nis

		if( nis .gt. 0) then
			write(filename,'(a,a)') trim(dirname),'/Matlab_friendly_data/normalized_in_situ.out'
			open(unit=13,file=filename,status="old",action="read",form='formatted')
			do i = 1, nis
				read(13,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,20(1x,es10.2))') is_id(i,2), is_redshift(i,2), is_rad(i,2), is_gas_mass(i,2), &
				& is_star_mass(i,2), is_mass(i,2), is_gas_frac(i,2), is_dm_frac(i,2), is_gas_sig(i,2), is_star_sig(i,2), is_sig(i,2), is_age(i,2), &
				& is_gas_met(i,2), is_star_met(i,2), is_SFR(i,2), is_sig_SFR(i,2), is_sSFR(i,2), is_tau(i,2), is_d_Rd(i,2), is_z_Hd(i,2), is_d(i,2), is_z(i,2), &
				& is_residual(i,2), is_shape(i,2), is_dm(i,2), is_es(i,2), is_merger(i,2), is_tff(i,2), is_td(i,2), is_td_global(i,2), &
				& is_mgas_in(i,2,1), is_mgas_in(i,2,2), is_mgas_in(i,2,3), &
				& is_mgas_out(i,2,1), is_mgas_out(i,2,2), is_mgas_out(i,2,3), is_mgas_out(i,2,4), is_mgas_out(i,2,5), is_mgas_out(i,2,6), is_mgas_out(i,2,7), &
				& is_mgas_out(i,2,8), is_mgas_out(i,2,9), is_mgas_out(i,2,10), &  
				& is_mstars_in(i,2), is_mstars_out(i,2), is_mstars_formed(i,2), is_alpha(i,2)
			end do
			close(unit=13)
		end if

		write(filename,'(a,a)') trim(dirname),'/Matlab_friendly_data/ex_situ.out'
		open(unit=14,file=filename,status="old",action="read",form='formatted')
		nes = 0
		do
			read(14,*,end=7)
			nes = nes + 1
		end do
 7		rewind(unit=14)
		if( nes .gt. 0 ) then
			allocate( es_redshift(nes,2), es_rad(nes,2), es_gas_mass(nes,2), es_star_mass(nes,2), es_mass(nes,2) )
			allocate( es_gas_frac(nes,2), es_dm_frac(nes,2), es_gas_sig(nes,2), es_star_sig(nes,2), es_sig(nes,2) )
			allocate( es_age(nes,2), es_gas_met(nes,2), es_star_met(nes,2), es_SFR(nes,2), es_sig_SFR(nes,2), es_sSFR(nes,2), es_tau(nes,2) )
			allocate( es_d_Rd(nes,2), es_z_Hd(nes,2), es_d(nes,2), es_z(nes,2), es_residual(nes,2), es_shape(nes,2), es_dm(nes,2) )
			allocate( es_tff(nes,2), es_td(nes,2), es_td_global(nes,2), es_mgas_in(nes,2,3), es_mgas_out(nes,2,10) )
			allocate( es_mstars_in(nes,2), es_mstars_out(nes,2), es_mstars_formed(nes,2), es_alpha(nes,2) )
			allocate( es_es(nes,2), es_merger(nes,2), es_id(nes,2), es_list(nes) )
			es_list(:) = 1
			do i = 1, nes
				read(14,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,20(1x,es10.2))') es_id(i,1), es_redshift(i,1), es_rad(i,1), es_gas_mass(i,1), &
				& es_star_mass(i,1), es_mass(i,1), es_gas_frac(i,1), es_dm_frac(i,1), es_gas_sig(i,1), es_star_sig(i,1), es_sig(i,1), es_age(i,1), &
				& es_gas_met(i,1), es_star_met(i,1), es_SFR(i,1), es_sig_SFR(i,1), es_sSFR(i,1), es_tau(i,1), es_d_Rd(i,1), es_z_Hd(i,1), es_d(i,1), es_z(i,1), &
				& es_residual(i,1), es_shape(i,1), es_dm(i,1), es_es(i,1), es_merger(i,1), es_tff(i,1), es_td(i,1), es_td_global(i,1), &
				& es_mgas_in(i,1,1), es_mgas_in(i,1,2), es_mgas_in(i,1,3), &
				& es_mgas_out(i,1,1), es_mgas_out(i,1,2), es_mgas_out(i,1,3), es_mgas_out(i,1,4), es_mgas_out(i,1,5), es_mgas_out(i,1,6), es_mgas_out(i,1,7), &
				& es_mgas_out(i,1,8), es_mgas_out(i,1,9), es_mgas_out(i,1,10), &  
				& es_mstars_in(i,1), es_mstars_out(i,1), es_mstars_formed(i,1), es_alpha(i,1)
			end do
		else
			allocate( es_id(1,2), es_redshift(1,2), es_list(1) )
			es_id(:,:) = 0
			es_redshift(:,:) = 0.0_4
			es_list(:) = 0
		end if
		close(unit=14)
		print *, 'nes=', nes

		if( nes .gt. 0 ) then
			write(filename,'(a,a)') trim(dirname),'/Matlab_friendly_data/normalized_ex_situ.out'
			open(unit=15,file=filename,status="old",action="read",form='formatted')
			do i = 1, nes
				read(15,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,20(1x,es10.2))') es_id(i,2), es_redshift(i,2), es_rad(i,2), es_gas_mass(i,2), &
				& es_star_mass(i,2), es_mass(i,2), es_gas_frac(i,2), es_dm_frac(i,2), es_gas_sig(i,2), es_star_sig(i,2), es_sig(i,2), es_age(i,2), &
				& es_gas_met(i,2), es_star_met(i,2), es_SFR(i,2), es_sig_SFR(i,2), es_sSFR(i,2), es_tau(i,2), es_d_Rd(i,2), es_z_Hd(i,2), es_d(i,2), es_z(i,2), &
				& es_residual(i,2), es_shape(i,2), es_dm(i,2), es_es(i,2), es_merger(i,2), es_tff(i,2), es_td(i,2), es_td_global(i,2), &
				& es_mgas_in(i,2,1), es_mgas_in(i,2,2), es_mgas_in(i,2,3), &
				& es_mgas_out(i,2,1), es_mgas_out(i,2,2), es_mgas_out(i,2,3), es_mgas_out(i,2,4), es_mgas_out(i,2,5), es_mgas_out(i,2,6), es_mgas_out(i,2,7), &
				& es_mgas_out(i,2,8), es_mgas_out(i,2,9), es_mgas_out(i,2,10), &  
				& es_mstars_in(i,2), es_mstars_out(i,2), es_mstars_formed(i,2), es_alpha(i,2)
			end do
			close(unit=15)
		end if

		write(filename,'(a,a)') trim(dirname),'/Matlab_friendly_data/bulge.out'
		open(unit=16,file=filename,status="old",action="read",form='formatted')
		nbulge = 0
		do
			read(16,*,end=8)
			nbulge = nbulge + 1
		end do
 8		rewind(unit=16)
		if( nbulge .gt. 0 ) then
			allocate( bulge_redshift(nbulge,2), bulge_rad(nbulge,2), bulge_gas_mass(nbulge,2), bulge_star_mass(nbulge,2), bulge_mass(nbulge,2) )
			allocate( bulge_gas_frac(nbulge,2), bulge_dm_frac(nbulge,2), bulge_gas_sig(nbulge,2), bulge_star_sig(nbulge,2), bulge_sig(nbulge,2) )
			allocate( bulge_age(nbulge,2), bulge_gas_met(nbulge,2), bulge_star_met(nbulge,2), bulge_SFR(nbulge,2), bulge_sig_SFR(nbulge,2) )
			allocate( bulge_sSFR(nbulge,2), bulge_tau(nbulge,2) )
			allocate( bulge_d_Rd(nbulge,2), bulge_z_Hd(nbulge,2), bulge_d(nbulge,2), bulge_z(nbulge,2), bulge_residual(nbulge,2), bulge_shape(nbulge,2), bulge_dm(nbulge,2) )
			allocate( bulge_tff(nbulge,2), bulge_td(nbulge,2), bulge_td_global(nbulge,2), bulge_mgas_in(nbulge,2,3), bulge_mgas_out(nbulge,2,10) )
			allocate( bulge_mstars_in(nbulge,2), bulge_mstars_out(nbulge,2), bulge_mstars_formed(nbulge,2), bulge_alpha(nbulge,2) )
			allocate( bulge_es(nbulge,2), bulge_merger(nbulge,2), bulge_id(nbulge,2), bulge_list(nbulge) )
			bulge_list(:) = 1
			do i = 1, nbulge
				read(16,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,20(1x,es10.2))') bulge_id(i,1), bulge_redshift(i,1), bulge_rad(i,1), bulge_gas_mass(i,1), &
				& bulge_star_mass(i,1), bulge_mass(i,1), bulge_gas_frac(i,1), bulge_dm_frac(i,1), bulge_gas_sig(i,1), bulge_star_sig(i,1), bulge_sig(i,1), bulge_age(i,1), &
				& bulge_gas_met(i,1), bulge_star_met(i,1), bulge_SFR(i,1), bulge_sig_SFR(i,1), bulge_sSFR(i,1), bulge_tau(i,1), bulge_d_Rd(i,1), bulge_z_Hd(i,1), &
				& bulge_d(i,1), bulge_z(i,1), bulge_residual(i,1), bulge_shape(i,1), bulge_dm(i,1), bulge_es(i,1), bulge_merger(i,1), bulge_tff(i,1), &
				& bulge_td(i,1), bulge_td_global(i,1), bulge_mgas_in(i,1,1), bulge_mgas_in(i,1,2), bulge_mgas_in(i,1,3), &
				& bulge_mgas_out(i,1,1), bulge_mgas_out(i,1,2), bulge_mgas_out(i,1,3), bulge_mgas_out(i,1,4), bulge_mgas_out(i,1,5), bulge_mgas_out(i,1,6), &
				& bulge_mgas_out(i,1,7), bulge_mgas_out(i,1,8), bulge_mgas_out(i,1,9), bulge_mgas_out(i,1,10), &
				& bulge_mstars_in(i,1), bulge_mstars_out(i,1), bulge_mstars_formed(i,1), bulge_alpha(i,1)
			end do
		else
			allocate( bulge_id(1,2), bulge_redshift(1,2), bulge_list(1) )
			bulge_id(:,:) = 0
			bulge_redshift(:,:) = 0.0_4
			bulge_list(:) = 0
		end if
		close(unit=16)
		print *, 'nbulge=', nbulge

		if( nbulge .gt. 0 ) then
			write(filename,'(a,a)') trim(dirname),'/Matlab_friendly_data/normalized_bulge.out'
			open(unit=17,file=filename,status="old",action="read",form='formatted')
			do i = 1, nbulge
				read(17,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,20(1x,es10.2))') bulge_id(i,2), bulge_redshift(i,2), bulge_rad(i,2), bulge_gas_mass(i,2), &
				& bulge_star_mass(i,2), bulge_mass(i,2), bulge_gas_frac(i,2), bulge_dm_frac(i,2), bulge_gas_sig(i,2), bulge_star_sig(i,2), bulge_sig(i,2), bulge_age(i,2), &
				& bulge_gas_met(i,2), bulge_star_met(i,2), bulge_SFR(i,2), bulge_sig_SFR(i,2), bulge_sSFR(i,2), bulge_tau(i,2), bulge_d_Rd(i,2), bulge_z_Hd(i,2), &
				& bulge_d(i,2), bulge_z(i,2), bulge_residual(i,2), bulge_shape(i,2), bulge_dm(i,2), bulge_es(i,2), bulge_merger(i,2), bulge_tff(i,2), &
				& bulge_td(i,2), bulge_td_global(i,2), bulge_mgas_in(i,2,1), bulge_mgas_in(i,2,2), bulge_mgas_in(i,2,3), &
				& bulge_mgas_out(i,2,1), bulge_mgas_out(i,2,2), bulge_mgas_out(i,2,3), bulge_mgas_out(i,2,4), bulge_mgas_out(i,2,5), bulge_mgas_out(i,2,6), & 
				& bulge_mgas_out(i,2,7), bulge_mgas_out(i,2,8), bulge_mgas_out(i,2,9), bulge_mgas_out(i,2,10), &
				& bulge_mstars_in(i,2), bulge_mstars_out(i,2), bulge_mstars_formed(i,2), bulge_alpha(i,2)
			end do
			close(unit=17)
		end if

		print *, 'verifying that all IDs in each snapshot are different'
		i = 1
		do while (i .lt. nis)
			do j=i+1,nis
				if( is_id(j,1) .eq. is_id(i,1) .and. is_redshift(j,1) .eq. is_redshift(i,1) ) then
					print *, 'there are 2 in situ clumps with ID', is_id(i,1)
					print *, 'at redshift', is_redshift(i,1)
					stop
				end if
			end do
			if( nes .gt. 0) then
				do j=1,nes
					if( es_id(j,1) .eq. is_id(i,1) .and. es_redshift(j,1) .eq. is_redshift(i,1) ) then
						print *, 'there are an in situ clump and an ex situ clump with ID', is_id(i,1)
						print *, 'at redshift', is_redshift(i,1)
						stop
					end if
				end do
			end if
			if( nbulge .gt. 0) then
				do j=1,nbulge
					if( bulge_id(j,1) .eq. is_id(i,1) .and. bulge_redshift(j,1) .eq. is_redshift(i,1) ) then
						print *, 'there are an in situ clump and a bulge clump with ID', is_id(i,1)
						print *, 'at redshift', is_redshift(i,1)
						stop
					end if
				end do
			end if
			i = i + 1
		end do
		i = 1
		do while (i .lt. nes)
			do j=i+1,nes
				if( es_id(j,1) .eq. es_id(i,1) .and. es_redshift(j,1) .eq. es_redshift(i,1) ) then
					print *, 'there are 2 ex situ clumps with ID', es_id(i,1)
					print *, 'at redshift', es_redshift(i,1)
					stop
				end if
			end do
			if( nbulge .gt. 0) then
				do j=1,nbulge
					if( bulge_id(j,1) .eq. es_id(i,1) .and. bulge_redshift(j,1) .eq. es_redshift(i,1) ) then
						print *, 'there are an ex situ clump and a bulge clump with ID', es_id(i,1)
						print *, 'at redshift', es_redshift(i,1)
						stop
					end if
				end do
			end if
			i = i + 1
		end do
		i = 1
		do while (i .lt. nbulge)
			do j=i+1,nbulge
				if( bulge_redshift(j,1) .eq. bulge_redshift(i,1) ) then
					print *, 'there are 2 bulge clumps at redshift', bulge_redshift(i,1)
					stop
				end if
			end do
			i = i + 1
		end do
		print *, 'looks like youre good to go!'

		print *, 'making files'
		write(filename,'(a,a,a,a,a)') './',trim(dirname),'/Matlab_friendly_data/same_clumps/in_situ.out'
		open(unit=18,file=filename,form='formatted')
		write(filename,'(a,a,a,a,a)') './',trim(dirname),'/Matlab_friendly_data/same_clumps/ex_situ.out'
		open(unit=19,file=filename,form='formatted')
		write(filename,'(a,a,a,a,a)') './',trim(dirname),'/Matlab_friendly_data/same_clumps/bulge.out'
		open(unit=20,file=filename,form='formatted')

		write(filename,'(a,a,a,a,a)') './',trim(dirname),'/Matlab_friendly_data/same_clumps/normalized_in_situ.out'
		open(unit=21,file=filename,form='formatted')
		write(filename,'(a,a,a,a,a)') './',trim(dirname),'/Matlab_friendly_data/same_clumps/normalized_ex_situ.out'
		open(unit=22,file=filename,form='formatted')
		write(filename,'(a,a,a,a,a)') './',trim(dirname),'/Matlab_friendly_data/same_clumps/normalized_bulge.out'
		open(unit=23,file=filename,form='formatted')

		allocate( ex_in(nis+nes+nbulge), bulge_in(nis+nes+nbulge), bulge_ex(nis+nes+nbulge) )
		ex_in(:) = 0
		bulge_in(:) = 0
		bulge_ex(:) = 0
		nex_in = 0
		nbulge_in = 0
		nbulge_ex = 0

		print *, 'in situ clumps'
		i = 1
		do while (i .le. nis)
			if( is_list(i) .ne. 0 ) then
				if( any( bulge_id(:,1).eq.is_id(i,1) ) ) then
					if( .not. any( bulge_redshift(:,1).eq.is_redshift(i,1) ) ) then
						nbulge_in = nbulge_in + 1
						bulge_in(nbulge_in) = i
					elseif( any(bulge_redshift(:,1).eq.is_redshift(i,1) ) ) then
						j = 1
						do while(j .lt. nbulge)
							if( bulge_redshift(j,1) .gt. is_redshift(i,1) ) then
								j = j + 1
							elseif( bulge_redshift(j,1) .eq. is_redshift(i,1) ) then
								n = bulge_id(j,1)
								bulge_id(j,:) = is_id(i,:)
								is_id(i,:) = n
								j = 2*nbulge + 10
								i = i - 1
							elseif( bulge_redshift(j,1) .lt. is_redshift(i,1) ) then
								print *, 'Something got fucked up when trying to correct bulge ID'
								print *, j,bulge_redshift(j,1),is_redshift(i,1)
								stop
							end if
						end do
					end if
				else if( any( es_id(:,1).eq.is_id(i,1) .and. es_redshift(:,1).gt.is_redshift(i,1) .and. es_list(:).eq.1 ) ) then
					nex_in = nex_in + 1
					ex_in(nex_in) = i
				else
					time = 0.95_4/(((1+is_redshift(i,1))/7.0_4)**1.5)  !!! Gyr
					dt = 0.0_4
					is_list(i) = 0
					write(18,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))') is_id(i,1), is_redshift(i,1), is_rad(i,1), is_gas_mass(i,1), &
					& is_star_mass(i,1), is_mass(i,1), is_gas_frac(i,1), is_dm_frac(i,1), is_gas_sig(i,1), is_star_sig(i,1), is_sig(i,1), is_age(i,1), &
					& is_gas_met(i,1), is_star_met(i,1), is_SFR(i,1), is_sig_SFR(i,1), is_sSFR(i,1), is_tau(i,1), is_d_Rd(i,1), is_z_Hd(i,1), is_d(i,1), is_z(i,1), &
					& is_residual(i,1), is_shape(i,1), is_dm(i,1), is_es(i,1), is_merger(i,1), is_tff(i,1), is_td(i,1), is_td_global(i,1), &
					& is_mgas_in(i,1,1), is_mgas_in(i,1,2), is_mgas_in(i,1,3), &
					& is_mgas_out(i,1,1), is_mgas_out(i,1,2), is_mgas_out(i,1,3), is_mgas_out(i,1,4), is_mgas_out(i,1,5), is_mgas_out(i,1,6), is_mgas_out(i,1,7), &
					& is_mgas_out(i,1,8), is_mgas_out(i,1,9), is_mgas_out(i,1,10), &
					& is_mstars_in(i,1), is_mstars_out(i,1), is_mstars_formed(i,1), is_alpha(i,1), dt
					write(21,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')is_id(i,2), is_redshift(i,2), is_rad(i,2), is_gas_mass(i,2), &
					& is_star_mass(i,2), is_mass(i,2), is_gas_frac(i,2), is_dm_frac(i,2), is_gas_sig(i,2), is_star_sig(i,2), is_sig(i,2), is_age(i,2), &
					& is_gas_met(i,2), is_star_met(i,2), is_SFR(i,2), is_sig_SFR(i,2), is_sSFR(i,2), is_tau(i,2), is_d_Rd(i,2), is_z_Hd(i,2), is_d(i,2), is_z(i,2), &
					& is_residual(i,2), is_shape(i,2), is_dm(i,2), is_es(i,2), is_merger(i,2), is_tff(i,2), is_td(i,2), is_td_global(i,2), &
					& is_mgas_in(i,2,1), is_mgas_in(i,2,2), is_mgas_in(i,2,3), &
					& is_mgas_out(i,2,1), is_mgas_out(i,2,2), is_mgas_out(i,2,3), is_mgas_out(i,2,4), is_mgas_out(i,2,5), is_mgas_out(i,2,6), is_mgas_out(i,2,7), &
					& is_mgas_out(i,2,8), is_mgas_out(i,2,9), is_mgas_out(i,2,10), &  
					& is_mstars_in(i,2), is_mstars_out(i,2), is_mstars_formed(i,2), is_alpha(i,2), dt

					ind = i
					if(i .lt. nis) then
						do j=i+1,nis
							if( is_id(j,1) .eq. is_id(i,1) ) then
								ind = j
								if( nes .gt. 0 ) then
									do l=1,nes
										if( es_id(l,1).eq.is_id(i,1) .and. es_redshift(l,1).lt.is_redshift(i,1) .and. es_redshift(l,1).gt.is_redshift(j,1) .and. es_list(l).eq.1 ) then
											dt = 0.95_4/(((1+es_redshift(l,1))/7.0_4)**1.5) - time      !!! Gyr
											es_list(l) = 0

											write(18,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(l,1), es_redshift(l,1), &
											& es_rad(l,1), es_gas_mass(l,1), es_star_mass(l,1), es_mass(l,1), es_gas_frac(l,1), & 
											& es_dm_frac(l,1), es_gas_sig(l,1), es_star_sig(l,1), es_sig(l,1), es_age(l,1), & 
											& es_gas_met(l,1), es_star_met(l,1), es_SFR(l,1), es_sig_SFR(l,1), es_sSFR(l,1), es_tau(l,1), &
											& es_d_Rd(l,1), es_z_Hd(l,1), es_d(l,1), es_z(l,1), &
											& es_residual(l,1), es_shape(l,1), es_dm(l,1), es_es(l,1), es_merger(l,1), es_tff(l,1), es_td(l,1), & 
											& es_td_global(l,1), es_mgas_in(l,1,1), es_mgas_in(l,1,2), es_mgas_in(l,1,3), & 
											& es_mgas_out(l,1,1), es_mgas_out(l,1,2), es_mgas_out(l,1,3), & 
											& es_mgas_out(l,1,4), es_mgas_out(l,1,5), es_mgas_out(l,1,6), es_mgas_out(l,1,7), &
											& es_mgas_out(l,1,8), es_mgas_out(l,1,9), es_mgas_out(l,1,10), &  
											& es_mstars_in(l,1), es_mstars_out(l,1), es_mstars_formed(l,1), es_alpha(l,1), dt
											write(21,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(l,2), es_redshift(l,2), & 
											& es_rad(l,2), es_gas_mass(l,2), es_star_mass(l,2), es_mass(l,2), es_gas_frac(l,2), & 
											& es_dm_frac(l,2), es_gas_sig(l,2), es_star_sig(l,2), es_sig(l,2), es_age(l,2), &
											& es_gas_met(l,2), es_star_met(l,2), es_SFR(l,2), es_sig_SFR(l,2), es_sSFR(l,2), es_tau(l,2), & 
											& es_d_Rd(l,2), es_z_Hd(l,2), es_d(l,2), es_z(l,2), & 
											& es_residual(l,2), es_shape(l,2), es_dm(l,2), es_es(l,2), es_merger(l,2), es_tff(l,2), es_td(l,2), & 
											& es_td_global(l,2), es_mgas_in(l,2,1), es_mgas_in(l,2,2), es_mgas_in(l,2,3), &
											& es_mgas_out(l,2,1), es_mgas_out(l,2,2), es_mgas_out(l,2,3), &
											& es_mgas_out(l,2,4), es_mgas_out(l,2,5), es_mgas_out(l,2,6), es_mgas_out(l,2,7), &
											& es_mgas_out(l,2,8), es_mgas_out(l,2,9), es_mgas_out(l,2,10), &
											& es_mstars_in(l,2), es_mstars_out(l,2), es_mstars_formed(l,2), es_alpha(l,2), dt
										end if
									end do
								end if

								dt = 0.95_4/(((1+is_redshift(j,1))/7.0_4)**1.5) - time      !!! Gyr
								is_list(j) = 0

								write(18,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')is_id(j,1), is_redshift(j,1), is_rad(j,1), & 
								& is_gas_mass(j,1), is_star_mass(j,1), is_mass(j,1), is_gas_frac(j,1), is_dm_frac(j,1), & 
								& is_gas_sig(j,1), is_star_sig(j,1), is_sig(j,1), is_age(j,1), is_gas_met(j,1), is_star_met(j,1), &
								& is_SFR(j,1), is_sig_SFR(j,1), is_sSFR(j,1), is_tau(j,1), is_d_Rd(j,1), is_z_Hd(j,1), is_d(j,1), is_z(j,1), &
								& is_residual(j,1), is_shape(j,1), is_dm(j,1), is_es(j,1), is_merger(j,1), is_tff(j,1), is_td(j,1), is_td_global(j,1), &
								& is_mgas_in(j,1,1), is_mgas_in(j,1,2), is_mgas_in(j,1,3), &
								& is_mgas_out(j,1,1), is_mgas_out(j,1,2), is_mgas_out(j,1,3), is_mgas_out(j,1,4), is_mgas_out(j,1,5), is_mgas_out(j,1,6), & 
								& is_mgas_out(j,1,7), is_mgas_out(j,1,8), is_mgas_out(j,1,9), is_mgas_out(j,1,10), &
								& is_mstars_in(j,1), is_mstars_out(j,1), is_mstars_formed(j,1), is_alpha(j,1), dt
								write(21,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')is_id(j,2), is_redshift(j,2), is_rad(j,2), & 
								& is_gas_mass(j,2), is_star_mass(j,2), is_mass(j,2), is_gas_frac(j,2), is_dm_frac(j,2), & 
								& is_gas_sig(j,2), is_star_sig(j,2), is_sig(j,2), is_age(j,2), is_gas_met(j,2), is_star_met(j,2), & 
								& is_SFR(j,2), is_sig_SFR(j,2), is_sSFR(j,2), is_tau(j,2), is_d_Rd(j,2), is_z_Hd(j,2), is_d(j,2), is_z(j,2), &
								& is_residual(j,2), is_shape(j,2), is_dm(j,2), is_es(j,2), is_merger(j,2), is_tff(j,2), is_td(j,2), is_td_global(j,2), &
								& is_mgas_in(j,2,1), is_mgas_in(j,2,2), is_mgas_in(j,2,3), &
								& is_mgas_out(j,2,1), is_mgas_out(j,2,2), is_mgas_out(j,2,3), is_mgas_out(j,2,4), is_mgas_out(j,2,5), is_mgas_out(j,2,6), & 
								& is_mgas_out(j,2,7), is_mgas_out(j,2,8), is_mgas_out(j,2,9), is_mgas_out(j,2,10), &
								& is_mstars_in(j,2), is_mstars_out(j,2), is_mstars_formed(j,2), is_alpha(j,2),dt
							end if
						end do
					end if
					if( nes .gt. 0 ) then
						do l=1,nes
							if( es_id(l,1).eq.is_id(ind,1) .and. es_redshift(l,1).lt.is_redshift(ind,1) .and. es_list(l).eq.1 ) then
								dt = 0.95_4/(((1+es_redshift(l,1))/7.0_4)**1.5) - time      !!! Gyr
								es_list(l) = 0

								write(18,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(l,1), es_redshift(l,1), es_rad(l,1), & 
								& es_gas_mass(l,1), es_star_mass(l,1), es_mass(l,1), es_gas_frac(l,1), es_dm_frac(l,1), & 
								& es_gas_sig(l,1), es_star_sig(l,1), es_sig(l,1), es_age(l,1), es_gas_met(l,1), es_star_met(l,1), & 
								& es_SFR(l,1), es_sig_SFR(l,1), es_sSFR(l,1), es_tau(l,1), es_d_Rd(l,1), es_z_Hd(l,1), es_d(l,1), es_z(l,1), &
								& es_residual(l,1), es_shape(l,1), es_dm(l,1), es_es(l,1), es_merger(l,1), es_tff(l,1), es_td(l,1), es_td_global(l,1), &
								& es_mgas_in(l,1,1), es_mgas_in(l,1,2), es_mgas_in(l,1,3), &
								& es_mgas_out(l,1,1), es_mgas_out(l,1,2), es_mgas_out(l,1,3), es_mgas_out(l,1,4), es_mgas_out(l,1,5), es_mgas_out(l,1,6), &
								& es_mgas_out(l,1,7), es_mgas_out(l,1,8), es_mgas_out(l,1,9), es_mgas_out(l,1,10), &
								& es_mstars_in(l,1), es_mstars_out(l,1), es_mstars_formed(l,1), es_alpha(l,1), dt
								write(21,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(l,2), es_redshift(l,2), es_rad(l,2), & 
								& es_gas_mass(l,2), es_star_mass(l,2), es_mass(l,2), es_gas_frac(l,2), es_dm_frac(l,2), & 
								& es_gas_sig(l,2), es_star_sig(l,2), es_sig(l,2), es_age(l,2), es_gas_met(l,2), es_star_met(l,2), & 
								& es_SFR(l,2), es_sig_SFR(l,2), es_sSFR(l,2), es_tau(l,2), es_d_Rd(l,2), es_z_Hd(l,2), es_d(l,2), es_z(l,2), &
								& es_residual(l,2), es_shape(l,2), es_dm(l,2), es_es(l,2), es_merger(l,2), es_tff(l,2), es_td(l,2), es_td_global(l,2), &
								& es_mgas_in(l,2,1), es_mgas_in(l,2,2), es_mgas_in(l,2,3), &
								& es_mgas_out(l,2,1), es_mgas_out(l,2,2), es_mgas_out(l,2,3), es_mgas_out(l,2,4), es_mgas_out(l,2,5), es_mgas_out(l,2,6), & 
								& es_mgas_out(l,2,7), es_mgas_out(l,2,8), es_mgas_out(l,2,9), es_mgas_out(l,2,10), &
								& es_mstars_in(l,2), es_mstars_out(l,2), es_mstars_formed(l,2), es_alpha(l,2), dt
							end if
						end do
					end if
				end if
			end if
			i = i + 1
		end do

		!!! Check for errors !!!
		do i=1,nis
			if( any(bulge_id(:,1).eq.is_id(i,1)) .and. .not. any(bulge_redshift(:,1).eq.is_redshift(i,1)) ) then
				if( is_list(i) .ne. 1) then
					print *, 'error in is_list 1', i
					stop
				end if
			elseif ( any( es_id(:,1).eq.is_id(i,1) .and. es_redshift(:,1).gt.is_redshift(i,1) .and. es_list(:).eq.1 ) ) then
				if( is_list(i) .ne. 1) then
					print *, 'error in is_list 2', i
					stop
				end if
			else
				if( is_list(i) .ne. 0) then
					print *, 'error in is_list 3', i
					stop
				end if
			end if
		end do

		do i=1,nes
			if ( any( is_id(:,1).eq.es_id(i,1) .and. is_redshift(:,1).gt.es_redshift(i,1) .and. is_list(:).eq.0 ) ) then
				if( es_list(i) .ne. 0) then
					print *, 'error in es_list 1', i
					stop
				end if
			elseif( es_list(i) .ne. 1) then
				print *, 'error in es_list 2', i
				stop
			end if
		end do

		do i=1,nbulge
			if( bulge_list(i) .ne. 1) then
				print *, 'error in bulge_list 1', i
				stop
			end if
		end do
		!!!!!!!!!!!!!!!!!!!!!!!!!!

		print *, 'ex situ clumps'
		i = 1
		do while (i .le. nes)
			if( es_list(i) .ne. 0 ) then
				if( any(bulge_id(:,1).eq.es_id(i,1)) .and. .not. any(bulge_redshift(:,1).eq.es_redshift(i,1)) ) then
					nbulge_ex = nbulge_ex + 1
					bulge_ex(nbulge_ex) = i
				else if( any(is_id(:,1).eq.es_id(i,1) .and. is_list(:).eq.1 ) ) then
					time = 0.95_4/(((1+es_redshift(i,1))/7.0_4)**1.5)  !!! Gyr
					dt = 0.0_4
					es_list(i) = 0

					write(19,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(i,1), es_redshift(i,1), es_rad(i,1), es_gas_mass(i,1), &
					& es_star_mass(i,1), es_mass(i,1), es_gas_frac(i,1), es_dm_frac(i,1), es_gas_sig(i,1), es_star_sig(i,1), es_sig(i,1), es_age(i,1), &
					& es_gas_met(i,1), es_star_met(i,1), es_SFR(i,1), es_sig_SFR(i,1), es_sSFR(i,1), es_tau(i,1), es_d_Rd(i,1), es_z_Hd(i,1), es_d(i,1), es_z(i,1), &
					& es_residual(i,1), es_shape(i,1), es_dm(i,1), es_es(i,1), es_merger(i,1), es_tff(i,1), es_td(i,1), es_td_global(i,1), &
					& es_mgas_in(i,1,1), es_mgas_in(i,1,2), es_mgas_in(i,1,3), &
					& es_mgas_out(i,1,1), es_mgas_out(i,1,2), es_mgas_out(i,1,3), es_mgas_out(i,1,4), es_mgas_out(i,1,5), es_mgas_out(i,1,6), es_mgas_out(i,1,7), &
					& es_mgas_out(i,1,8), es_mgas_out(i,1,9), es_mgas_out(i,1,10), &  
					& es_mstars_in(i,1), es_mstars_out(i,1), es_mstars_formed(i,1), es_alpha(i,1), dt
					write(22,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(i,2), es_redshift(i,2), es_rad(i,2), es_gas_mass(i,2), &
					& es_star_mass(i,2), es_mass(i,2), es_gas_frac(i,2), es_dm_frac(i,2), es_gas_sig(i,2), es_star_sig(i,2), es_sig(i,2), es_age(i,2), &
					& es_gas_met(i,2), es_star_met(i,2), es_SFR(i,2), es_sig_SFR(i,2), es_sSFR(i,2), es_tau(i,2), es_d_Rd(i,2), es_z_Hd(i,2), es_d(i,2), es_z(i,2), &
					& es_residual(i,2), es_shape(i,2), es_dm(i,2), es_es(i,2), es_merger(i,2), es_tff(i,2), es_td(i,2), es_td_global(i,2), &
					& es_mgas_in(i,2,1), es_mgas_in(i,2,2), es_mgas_in(i,2,3), &
					& es_mgas_out(i,2,1), es_mgas_out(i,2,2), es_mgas_out(i,2,3), es_mgas_out(i,2,4), es_mgas_out(i,2,5), es_mgas_out(i,2,6), es_mgas_out(i,2,7), &
					& es_mgas_out(i,2,8), es_mgas_out(i,2,9), es_mgas_out(i,2,10), &  
					& es_mstars_in(i,2), es_mstars_out(i,2), es_mstars_formed(i,2), es_alpha(i,2), dt

					ind = i
					if(i .lt. nes) then
						do l=i+1,nes
							if( es_id(l,1) .eq. es_id(i,1) ) then
								ind = l
								do j=1,nex_in
									m = ex_in(j)
									if( is_id(m,1) .eq. es_id(l,1) .and. is_list(m) .eq. 1 .and. is_redshift(m,1) .lt. es_redshift(i,1) .and. is_redshift(m,1) .gt. es_redshift(l,1) ) then
										dt = 0.95_4/(((1+is_redshift(m,1))/7.0_4)**1.5) - time      !!! Gyr
										is_list(m) = 0

										write(19,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')is_id(m,1), is_redshift(m,1), & 
										& is_rad(m,1), is_gas_mass(m,1), is_star_mass(m,1), is_mass(m,1), is_gas_frac(m,1), is_dm_frac(m,1), & 
										& is_gas_sig(m,1), is_star_sig(m,1), is_sig(m,1), is_age(m,1), is_gas_met(m,1), is_star_met(m,1), & 
										& is_SFR(m,1), is_sig_SFR(m,1), is_sSFR(m,1), is_tau(m,1), is_d_Rd(m,1), is_z_Hd(m,1), is_d(m,1), is_z(m,1), &
										& is_residual(m,1), is_shape(m,1), is_dm(m,1), is_es(m,1), is_merger(m,1), is_tff(m,1), is_td(m,1), & 
										& is_td_global(m,1), is_mgas_in(m,1,1), is_mgas_in(m,1,2), is_mgas_in(m,1,3), &
										& is_mgas_out(m,1,1), is_mgas_out(m,1,2), is_mgas_out(m,1,3), &
										& is_mgas_out(m,1,4), is_mgas_out(m,1,5), is_mgas_out(m,1,6), is_mgas_out(m,1,7), &
										& is_mgas_out(m,1,8), is_mgas_out(m,1,9), is_mgas_out(m,1,10), &
										& is_mstars_in(m,1), is_mstars_out(m,1), is_mstars_formed(m,1), is_alpha(m,1), dt
										write(22,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')is_id(m,2), is_redshift(m,2), & 
										& is_rad(m,2), is_gas_mass(m,2), is_star_mass(m,2), is_mass(m,2), is_gas_frac(m,2), is_dm_frac(m,2), &
										& is_gas_sig(m,2), is_star_sig(m,2), is_sig(m,2), is_age(m,2), is_gas_met(m,2), is_star_met(m,2), & 
										& is_SFR(m,2), is_sig_SFR(m,2), is_sSFR(m,2), is_tau(m,2), is_d_Rd(m,2), is_z_Hd(m,2), is_d(m,2), is_z(m,2), &
										& is_residual(m,2), is_shape(m,2), is_dm(m,2), is_es(m,2), is_merger(m,2), is_tff(m,2), is_td(m,2), & 
										& is_td_global(m,2), is_mgas_in(m,2,1), is_mgas_in(m,2,2), is_mgas_in(m,2,3), &
										& is_mgas_out(m,2,1), is_mgas_out(m,2,2), is_mgas_out(m,2,3), & 
										& is_mgas_out(m,2,4), is_mgas_out(m,2,5), is_mgas_out(m,2,6), is_mgas_out(m,2,7), &
										& is_mgas_out(m,2,8), is_mgas_out(m,2,9), is_mgas_out(m,2,10), & 
										& is_mstars_in(m,2), is_mstars_out(m,2), is_mstars_formed(m,2), is_alpha(m,2), dt
									end if
								end do

								dt = 0.95_4/(((1+es_redshift(l,1))/7.0_4)**1.5) - time      !!! Gyr
								es_list(l) = 0
								write(19,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(l,1), es_redshift(l,1), es_rad(l,1), & 
								& es_gas_mass(l,1), es_star_mass(l,1), es_mass(l,1), es_gas_frac(l,1), es_dm_frac(l,1), & 
								& es_gas_sig(l,1), es_star_sig(l,1), es_sig(l,1), es_age(l,1), es_gas_met(l,1), es_star_met(l,1), & 
								& es_SFR(l,1), es_sig_SFR(l,1), es_sSFR(l,1), es_tau(l,1), es_d_Rd(l,1), es_z_Hd(l,1), es_d(l,1), es_z(l,1), &
								& es_residual(l,1), es_shape(l,1), es_dm(l,1), es_es(l,1), es_merger(l,1), es_tff(l,1), es_td(l,1), es_td_global(l,1), &
								& es_mgas_in(l,1,1), es_mgas_in(l,1,2), es_mgas_in(l,1,3), &
								& es_mgas_out(l,1,1), es_mgas_out(l,1,2), es_mgas_out(l,1,3), es_mgas_out(l,1,4), es_mgas_out(l,1,5), es_mgas_out(l,1,6), & 
								& es_mgas_out(l,1,7), es_mgas_out(l,1,8), es_mgas_out(l,1,9), es_mgas_out(l,1,10), &
								& es_mstars_in(l,1), es_mstars_out(l,1), es_mstars_formed(l,1), es_alpha(l,1), dt
								write(22,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(l,2), es_redshift(l,2), es_rad(l,2), & 
								& es_gas_mass(l,2), es_star_mass(l,2), es_mass(l,2), es_gas_frac(l,2), es_dm_frac(l,2), & 
								& es_gas_sig(l,2), es_star_sig(l,2), es_sig(l,2), es_age(l,2), es_gas_met(l,2), es_star_met(l,2), & 
								& es_SFR(l,2), es_sig_SFR(l,2), es_sSFR(l,2), es_tau(l,2), es_d_Rd(l,2), es_z_Hd(l,2), es_d(l,2), es_z(l,2), &
								& es_residual(l,2), es_shape(l,2), es_dm(l,2), es_es(l,2), es_merger(l,2), es_tff(l,2), es_td(l,2), es_td_global(l,2), &
								& es_mgas_in(l,2,1), es_mgas_in(l,2,2), es_mgas_in(l,2,3), &
								& es_mgas_out(l,2,1), es_mgas_out(l,2,2), es_mgas_out(l,2,3), es_mgas_out(l,2,4), es_mgas_out(l,2,5), es_mgas_out(l,2,6), &
								& es_mgas_out(l,2,7), es_mgas_out(l,2,8), es_mgas_out(l,2,9), es_mgas_out(l,2,10), &
								& es_mstars_in(l,2), es_mstars_out(l,2), es_mstars_formed(l,2), es_alpha(l,2), dt
							end if
						end do
					end if
					do j=1,nex_in
						m = ex_in(j)
						if( is_id(m,1) .eq. es_id(ind,1) .and. is_list(m) .eq. 1 .and. is_redshift(m,1) .lt. es_redshift(ind,1) ) then
							dt = 0.95_4/(((1+is_redshift(m,1))/7.0_4)**1.5) - time      !!! Gyr
							is_list(m) = 0

							write(19,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')is_id(m,1), is_redshift(m,1), is_rad(m,1), is_gas_mass(m,1), &
							& is_star_mass(m,1), is_mass(m,1), is_gas_frac(m,1), is_dm_frac(m,1), is_gas_sig(m,1), is_star_sig(m,1), is_sig(m,1), & 
							& is_age(m,1), is_gas_met(m,1), is_star_met(m,1), is_SFR(m,1), is_sig_SFR(m,1), is_sSFR(m,1), is_tau(m,1), & 
							& is_d_Rd(m,1), is_z_Hd(m,1), is_d(m,1), is_z(m,1), is_residual(m,1), is_shape(m,1), is_dm(m,1), is_es(m,1), is_merger(m,1), & 
							& is_tff(m,1), is_td(m,1), is_td_global(m,1), is_mgas_in(m,1,1), is_mgas_in(m,1,2), is_mgas_in(m,1,3), is_mgas_out(m,1,1), &
							& is_mgas_out(m,1,2), is_mgas_out(m,1,3), is_mgas_out(m,1,4), is_mgas_out(m,1,5), is_mgas_out(m,1,6), is_mgas_out(m,1,7), &
							& is_mgas_out(m,1,8), is_mgas_out(m,1,9), is_mgas_out(m,1,10), &
							& is_mstars_in(m,1), is_mstars_out(m,1), is_mstars_formed(m,1), is_alpha(m,1), dt
							write(22,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')is_id(m,2), is_redshift(m,2), is_rad(m,2), is_gas_mass(m,2), &
							& is_star_mass(m,2), is_mass(m,2), is_gas_frac(m,2), is_dm_frac(m,2), is_gas_sig(m,2), is_star_sig(m,2), is_sig(m,2), & 
							& is_age(m,2), is_gas_met(m,2), is_star_met(m,2), is_SFR(m,2), is_sig_SFR(m,2), is_sSFR(m,2), is_tau(m,2), & 
							& is_d_Rd(m,2), is_z_Hd(m,2), is_d(m,2), is_z(m,2), is_residual(m,2), is_shape(m,2), is_dm(m,2), is_es(m,2), is_merger(m,2), & 
							& is_tff(m,2), is_td(m,2), is_td_global(m,2), is_mgas_in(m,2,1), is_mgas_in(m,2,2), is_mgas_in(m,2,3), is_mgas_out(m,2,1), &
							& is_mgas_out(m,2,2), is_mgas_out(m,2,3), is_mgas_out(m,2,4), is_mgas_out(m,2,5), is_mgas_out(m,2,6), is_mgas_out(m,2,7), &
							& is_mgas_out(m,2,8), is_mgas_out(m,2,9), is_mgas_out(m,2,10), &  
							& is_mstars_in(m,2), is_mstars_out(m,2), is_mstars_formed(m,2), is_alpha(m,2), dt
						end if
					end do
				else
					time = 0.95_4/(((1+es_redshift(i,1))/7.0_4)**1.5)  !!! Gyr
					dt = 0
					es_list(i) = 0

					write(19,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(i,1), es_redshift(i,1), es_rad(i,1), es_gas_mass(i,1), &
					& es_star_mass(i,1), es_mass(i,1), es_gas_frac(i,1), es_dm_frac(i,1), es_gas_sig(i,1), es_star_sig(i,1), es_sig(i,1), es_age(i,1), &
					& es_gas_met(i,1), es_star_met(i,1), es_SFR(i,1), es_sig_SFR(i,1), es_sSFR(i,1), es_tau(i,1), es_d_Rd(i,1), es_z_Hd(i,1), es_d(i,1), es_z(i,1), &
					& es_residual(i,1), es_shape(i,1), es_dm(i,1), es_es(i,1), es_merger(i,1), es_tff(i,1), es_td(i,1), es_td_global(i,1), &
					& es_mgas_in(i,1,1), es_mgas_in(i,1,2), es_mgas_in(i,1,3), &
					& es_mgas_out(i,1,1), es_mgas_out(i,1,2), es_mgas_out(i,1,3), es_mgas_out(i,1,4), es_mgas_out(i,1,5), es_mgas_out(i,1,6), es_mgas_out(i,1,7), &
					& es_mgas_out(i,1,8), es_mgas_out(i,1,9), es_mgas_out(i,1,10), &
					& es_mstars_in(i,1), es_mstars_out(i,1), es_mstars_formed(i,1), es_alpha(i,1), dt
					write(22,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(i,2), es_redshift(i,2), es_rad(i,2), es_gas_mass(i,2), &
					& es_star_mass(i,2), es_mass(i,2), es_gas_frac(i,2), es_dm_frac(i,2), es_gas_sig(i,2), es_star_sig(i,2), es_sig(i,2), es_age(i,2), &
					& es_gas_met(i,2), es_star_met(i,2), es_SFR(i,2), es_sig_SFR(i,2), es_sSFR(i,2), es_tau(i,2), es_d_Rd(i,2), es_z_Hd(i,2), es_d(i,2), es_z(i,2), &
					& es_residual(i,2), es_shape(i,2), es_dm(i,2), es_es(i,2), es_merger(i,2), es_tff(i,2), es_td(i,2), es_td_global(i,2), &
					& es_mgas_in(i,2,1), es_mgas_in(i,2,2), es_mgas_in(i,2,3), &
					& es_mgas_out(i,2,1), es_mgas_out(i,2,2), es_mgas_out(i,2,3), es_mgas_out(i,2,4), es_mgas_out(i,2,5), es_mgas_out(i,2,6), es_mgas_out(i,2,7), &
					& es_mgas_out(i,2,8), es_mgas_out(i,2,9), es_mgas_out(i,2,10), &  
					& es_mstars_in(i,2), es_mstars_out(i,2), es_mstars_formed(i,2), es_alpha(i,2), dt

					if(i .lt. nes) then
						do j=i+1,nes
							if( es_id(j,1) .eq. es_id(i,1) ) then
								dt = 0.95_4/(((1+es_redshift(j,1))/7.0_4)**1.5) - time      !!! Gyr
								es_list(j) = 0

								write(19,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(j,1), es_redshift(j,1), es_rad(j,1), & 
								& es_gas_mass(j,1), es_star_mass(j,1), es_mass(j,1), es_gas_frac(j,1), es_dm_frac(j,1), es_gas_sig(j,1), & 
								& es_star_sig(j,1), es_sig(j,1), es_age(j,1), es_gas_met(j,1), es_star_met(j,1), es_SFR(j,1), es_sig_SFR(j,1), & 
								& es_sSFR(j,1), es_tau(j,1), es_d_Rd(j,1), es_z_Hd(j,1), es_d(j,1), es_z(j,1), es_residual(j,1), es_shape(j,1), & 
								& es_dm(j,1), es_es(j,1), es_merger(j,1), es_tff(j,1), es_td(j,1), es_td_global(j,1), &
								& es_mgas_in(j,1,1), es_mgas_in(j,1,2), es_mgas_in(j,1,3), es_mgas_out(j,1,1), es_mgas_out(j,1,2), &
								& es_mgas_out(j,1,3), es_mgas_out(j,1,4), es_mgas_out(j,1,5), es_mgas_out(j,1,6), es_mgas_out(j,1,7), &
								& es_mgas_out(j,1,8), es_mgas_out(j,1,9), es_mgas_out(j,1,10), &  
								& es_mstars_in(j,1), es_mstars_out(j,1), es_mstars_formed(j,1), es_alpha(j,1), dt
								write(22,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(j,2), es_redshift(j,2), es_rad(j,2), & 
								& es_gas_mass(j,2), es_star_mass(j,2), es_mass(j,2), es_gas_frac(j,2), es_dm_frac(j,2), es_gas_sig(j,2), & 
								& es_star_sig(j,2), es_sig(j,2), es_age(j,2), es_gas_met(j,2), es_star_met(j,2), es_SFR(j,2), es_sig_SFR(j,2), & 
								& es_sSFR(j,2), es_tau(j,2), es_d_Rd(j,2), es_z_Hd(j,2), es_d(j,2), es_z(j,2), es_residual(j,2), es_shape(j,2), & 
								& es_dm(j,2), es_es(j,2), es_merger(j,2), es_tff(j,2), es_td(j,2), es_td_global(j,2), &
								& es_mgas_in(j,2,1), es_mgas_in(j,2,2), es_mgas_in(j,2,3), es_mgas_out(j,2,1), es_mgas_out(j,2,2), &
								& es_mgas_out(j,2,3), es_mgas_out(j,2,4), es_mgas_out(j,2,5), es_mgas_out(j,2,6), es_mgas_out(j,2,7), &
								& es_mgas_out(j,2,8), es_mgas_out(j,2,9), es_mgas_out(j,2,10), &
								& es_mstars_in(j,2), es_mstars_out(j,2), es_mstars_formed(j,2), es_alpha(j,2), dt
							end if
						end do
					end if
				end if
			end if
			i = i + 1
		end do

		!!! Check for errors !!!
		do i=1,nes
			if( any(bulge_id(:,1).eq.es_id(i,1)) .and. .not. any(bulge_redshift(:,1).eq.es_redshift(i,1)) ) then
				if( es_list(i) .ne. 1) then
					print *, 'error in es_list 3', i
					stop
				end if
			else
				if( es_list(i) .ne. 0) then
					print *, 'error in es_list 4', i
					stop
				end if
			end if
		end do          

		do i=1,nis
			if( any(bulge_id(:,1).eq.is_id(i,1)) .and. .not. any(bulge_redshift(:,1).eq.is_redshift(i,1)) ) then
				if( is_list(i) .ne. 1) then
					print *, 'error in is_list 4', i
					stop
				end if
			else
				if( is_list(i) .ne. 0) then
					print *, 'error in is_list 5', i
					stop
				end if
			end if
		end do

		do i=1,nbulge
			if( bulge_list(i) .ne. 1) then
				print *, 'error in bulge_list 2', i
				stop
			end if
		end do
		!!!!!!!!!!!!!!!!!!!!!!!!!!

		print *, 'bulge clumps'
		time = maxval(bulge_redshift(1:nbulge,1))
		if( nbulge_in .gt. 0 ) then
			time2 = maxval(is_redshift(bulge_in(1:nbulge_in),1), mask = is_list(bulge_in(1:nbulge_in)) .eq. 1)
			time = max(time, time2)
		end if
		if( nbulge_ex .gt. 0) then
			time2 = maxval(es_redshift(bulge_ex(1:nbulge_ex),1), mask = es_list(bulge_ex(1:nbulge_ex)) .eq. 1)
			time = max(time, time2)
		end if
		time = 0.95_4/(((1+time)/7.0_4)**1.5)  !!! Gyr

		do i=1,nbulge
			if( nbulge_in .gt. 0 ) then
				do j=1,nbulge_in
					m = bulge_in(j)
					if( is_list(m) .eq. 1 .and. is_redshift(m,1) .gt. bulge_redshift(i,1) ) then
						if( nbulge_ex .gt. 0 ) then
							do l=1,nbulge_ex
								n = bulge_ex(l)
								if( es_list(n) .eq. 1 .and. es_redshift(n,1) .gt. is_redshift(m,1) ) then
									dt = 0.95_4/(((1+es_redshift(n,1))/7.0_4)**1.5) - time      !!! Gyr
									es_list(n) = 0

									write(20,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(n,1), es_redshift(n,1), & 
									& es_rad(n,1), es_gas_mass(n,1), es_star_mass(n,1), es_mass(n,1), es_gas_frac(n,1), es_dm_frac(n,1), & 
									& es_gas_sig(n,1), es_star_sig(n,1), es_sig(n,1), es_age(n,1), es_gas_met(n,1), es_star_met(n,1), & 
									& es_SFR(n,1), es_sig_SFR(n,1), es_sSFR(n,1), es_tau(n,1), es_d_Rd(n,1), es_z_Hd(n,1), es_d(n,1), es_z(n,1), & 
									& es_residual(n,1), es_shape(n,1), es_dm(n,1), es_es(n,1), es_merger(n,1), es_tff(n,1), es_td(n,1), & 
									& es_td_global(n,1), es_mgas_in(n,1,1), es_mgas_in(n,1,2), es_mgas_in(n,1,3), &
									& es_mgas_out(n,1,1), es_mgas_out(n,1,2), es_mgas_out(n,1,3), & 
									& es_mgas_out(n,1,4), es_mgas_out(n,1,5), es_mgas_out(n,1,6), es_mgas_out(n,1,7), & 
									& es_mgas_out(n,1,8), es_mgas_out(n,1,9), es_mgas_out(n,1,10), & 
									& es_mstars_in(n,1), es_mstars_out(n,1), es_mstars_formed(n,1), es_alpha(n,1), dt
									write(23,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(n,2), es_redshift(n,2), & 
									& es_rad(n,2), es_gas_mass(n,2), es_star_mass(n,2), es_mass(n,2), es_gas_frac(n,2), es_dm_frac(n,2), & 
									& es_gas_sig(n,2), es_star_sig(n,2), es_sig(n,2), es_age(n,2), es_gas_met(n,2), es_star_met(n,2), & 
									& es_SFR(n,2), es_sig_SFR(n,2), es_sSFR(n,2), es_tau(n,2), es_d_Rd(n,2), es_z_Hd(n,2), es_d(n,2), es_z(n,2), & 
									& es_residual(n,2), es_shape(n,2), es_dm(n,2), es_es(n,2), es_merger(n,2), es_merger(n,2), es_tff(n,2), es_td(n,2), & 
									& es_td_global(n,2), es_mgas_in(n,2,1), es_mgas_in(n,2,2), es_mgas_in(n,2,3), & 
									& es_mgas_out(n,2,1), es_mgas_out(n,2,2), es_mgas_out(n,2,3), & 
									& es_mgas_out(n,2,4), es_mgas_out(n,2,5), es_mgas_out(n,2,6), es_mgas_out(n,2,7), & 
									& es_mgas_out(n,2,8), es_mgas_out(n,2,9), es_mgas_out(n,2,10), & 
									& es_mstars_in(n,2), es_mstars_out(n,2), es_mstars_formed(n,2), es_alpha(n,2), dt
								end if
							end do
						end if
						dt = 0.95_4/(((1+is_redshift(m,1))/7.0_4)**1.5) - time      !!! Gyr
						is_list(m) = 0

						write(20,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')is_id(m,1), is_redshift(m,1), is_rad(m,1), is_gas_mass(m,1), &
						& is_star_mass(m,1), is_mass(m,1), is_gas_frac(m,1), is_dm_frac(m,1), is_gas_sig(m,1), is_star_sig(m,1), is_sig(m,1), is_age(m,1), &
						& is_gas_met(m,1), is_star_met(m,1), is_SFR(m,1), is_sig_SFR(m,1), is_sSFR(m,1), is_tau(m,1), is_d_Rd(m,1), is_z_Hd(m,1), is_d(m,1), & 
						& is_z(m,1), is_residual(m,1), is_shape(m,1), is_dm(m,1), is_es(m,1), is_merger(m,1), is_tff(m,1), is_td(m,1), is_td_global(m,1), &
						& is_mgas_in(m,1,1), is_mgas_in(m,1,2), is_mgas_in(m,1,3), is_mgas_out(m,1,1), is_mgas_out(m,1,2), &
						& is_mgas_out(m,1,3), is_mgas_out(m,1,4), is_mgas_out(m,1,5), is_mgas_out(m,1,6), is_mgas_out(m,1,7), &
						& is_mgas_out(m,1,8), is_mgas_out(m,1,9), is_mgas_out(m,1,10), &  
						& is_mstars_in(m,1), is_mstars_out(m,1), is_mstars_formed(m,1), is_alpha(m,1), dt
						write(23,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')is_id(m,2), is_redshift(m,2), is_rad(m,2), is_gas_mass(m,2), &
						& is_star_mass(m,2), is_mass(m,2), is_gas_frac(m,2), is_dm_frac(m,2), is_gas_sig(m,2), is_star_sig(m,2), is_sig(m,2), is_age(m,2), &
						& is_gas_met(m,2), is_star_met(m,2), is_SFR(m,2), is_sig_SFR(m,2), is_sSFR(m,2), is_tau(m,2), is_d_Rd(m,2), is_z_Hd(m,2), is_d(m,2), & 
						& is_z(m,2), is_residual(m,2), is_shape(m,2), is_dm(m,2), is_es(m,2), is_merger(m,2), is_tff(m,2), is_td(m,2), is_td_global(m,2), &
						& is_mgas_in(m,2,1), is_mgas_in(m,2,2), is_mgas_in(m,2,3), is_mgas_out(m,2,1), is_mgas_out(m,2,2), &
						& is_mgas_out(m,2,3), is_mgas_out(m,2,4), is_mgas_out(m,2,5), is_mgas_out(m,2,6), is_mgas_out(m,2,7), &
						& is_mgas_out(m,2,8), is_mgas_out(m,2,9), is_mgas_out(m,2,10), &  
						& is_mstars_in(m,2), is_mstars_out(m,2), is_mstars_formed(m,2), is_alpha(m,2), dt
					end if
				end do
			else
				if( nbulge_ex .gt. 0 ) then
					do l=1,nbulge_ex
						n = bulge_ex(l)
						if( es_list(n) .eq. 1 .and. es_redshift(n,1) .gt. bulge_redshift(i,1) ) then
							dt = 0.95_4/(((1+es_redshift(n,1))/7.0_4)**1.5) - time      !!! Gyr
							es_list(n) = 0

							write(20,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(n,1), es_redshift(n,1), es_rad(n,1), es_gas_mass(n,1), &
							& es_star_mass(n,1), es_mass(n,1), es_gas_frac(n,1), es_dm_frac(n,1), es_gas_sig(n,1), es_star_sig(n,1), es_sig(n,1), es_age(n,1), &
							& es_gas_met(n,1), es_star_met(n,1), es_SFR(n,1), es_sig_SFR(n,1), es_sSFR(n,1), es_tau(n,1), es_d_Rd(n,1), es_z_Hd(n,1), es_d(n,1), &
							& es_z(n,1), es_residual(n,1), es_shape(n,1), es_dm(n,1), es_es(n,1), es_merger(n,1), es_tff(n,1), es_td(n,1), es_td_global(n,1), &
							& es_mgas_in(n,1,1), es_mgas_in(n,1,2), es_mgas_in(n,1,3), es_mgas_out(n,1,1), es_mgas_out(n,1,2), &
							& es_mgas_out(n,1,3), es_mgas_out(n,1,4), es_mgas_out(n,1,5), es_mgas_out(n,1,6), es_mgas_out(n,1,7), &
							& es_mgas_out(n,1,8), es_mgas_out(n,1,9), es_mgas_out(n,1,10), &  
							& es_mstars_in(n,1), es_mstars_out(n,1), es_mstars_formed(n,1), es_alpha(n,1), dt
							write(23,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(n,2), es_redshift(n,2), es_rad(n,2), es_gas_mass(n,2), &
							& es_star_mass(n,2), es_mass(n,2), es_gas_frac(n,2), es_dm_frac(n,2), es_gas_sig(n,2), es_star_sig(n,2), es_sig(n,2), es_age(n,2), &
							& es_gas_met(n,2), es_star_met(n,2), es_SFR(n,2), es_sig_SFR(n,2), es_sSFR(n,2), es_tau(n,2), es_d_Rd(n,2), es_z_Hd(n,2), es_d(n,2), &
							& es_z(n,2), es_residual(n,2), es_shape(n,2), es_dm(n,2), es_es(n,2), es_merger(n,2), es_tff(n,2), es_td(n,2), es_td_global(n,2), & 
							& es_mgas_in(n,2,1), es_mgas_in(n,2,2), es_mgas_in(n,2,3), es_mgas_out(n,2,1), es_mgas_out(n,2,2), &
							& es_mgas_out(n,2,3), es_mgas_out(n,2,4), es_mgas_out(n,2,5), es_mgas_out(n,2,6), es_mgas_out(n,2,7), &
							& es_mgas_out(n,2,8), es_mgas_out(n,2,9), es_mgas_out(n,2,10), &  
							& es_mstars_in(n,2), es_mstars_out(n,2), es_mstars_formed(n,2), es_alpha(n,2), dt
						end if
					end do
				end if
			end if
			dt = 0.95_4/(((1+bulge_redshift(i,1))/7.0_4)**1.5) - time      !!! Gyr
			bulge_list(i) = 0

			write(20,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')bulge_id(i,1), bulge_redshift(i,1), bulge_rad(i,1), bulge_gas_mass(i,1), &
			& bulge_star_mass(i,1), bulge_mass(i,1), bulge_gas_frac(i,1), bulge_dm_frac(i,1), bulge_gas_sig(i,1), bulge_star_sig(i,1), bulge_sig(i,1), & 
			& bulge_age(i,1), bulge_gas_met(i,1), bulge_star_met(i,1), bulge_SFR(i,1), bulge_sig_SFR(i,1), bulge_sSFR(i,1), bulge_tau(i,1), bulge_d_Rd(i,1), & 
			& bulge_z_Hd(i,1), bulge_d(i,1), bulge_z(i,1), bulge_residual(i,1), bulge_shape(i,1), bulge_dm(i,1), bulge_es(i,1), bulge_merger(i,1), bulge_tff(i,1), & 
			& bulge_td(i,1), bulge_td_global(i,1), bulge_mgas_in(i,1,1), bulge_mgas_in(i,1,2), bulge_mgas_in(i,1,3), & 
			& bulge_mgas_out(i,1,1), bulge_mgas_out(i,1,2), bulge_mgas_out(i,1,3), bulge_mgas_out(i,1,4), bulge_mgas_out(i,1,5), bulge_mgas_out(i,1,6), & 
			& bulge_mgas_out(i,1,7), bulge_mgas_out(i,1,8), bulge_mgas_out(i,1,9), bulge_mgas_out(i,1,10), &
			& bulge_mstars_in(i,1), bulge_mstars_out(i,1), bulge_mstars_formed(i,1), bulge_alpha(i,1), dt
			write(23,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')bulge_id(i,2), bulge_redshift(i,2), bulge_rad(i,2), bulge_gas_mass(i,2), &
			& bulge_star_mass(i,2), bulge_mass(i,2), bulge_gas_frac(i,2), bulge_dm_frac(i,2), bulge_gas_sig(i,2), bulge_star_sig(i,2), bulge_sig(i,2), & 
			& bulge_age(i,2), bulge_gas_met(i,2), bulge_star_met(i,2), bulge_SFR(i,2), bulge_sig_SFR(i,2), bulge_sSFR(i,2), bulge_tau(i,2), bulge_d_Rd(i,2), & 
			& bulge_z_Hd(i,2), bulge_d(i,2), bulge_z(i,2), bulge_residual(i,2), bulge_shape(i,2), bulge_dm(i,2), bulge_es(i,2), bulge_merger(i,2), bulge_tff(i,2), & 
			& bulge_td(i,2), bulge_td_global(i,2), bulge_mgas_in(i,2,1), bulge_mgas_in(i,2,2), bulge_mgas_in(i,2,3), & 
			& bulge_mgas_out(i,2,1), bulge_mgas_out(i,2,2), bulge_mgas_out(i,2,3), bulge_mgas_out(i,2,4), bulge_mgas_out(i,2,5), bulge_mgas_out(i,2,6), & 
			& bulge_mgas_out(i,2,7), bulge_mgas_out(i,2,8), bulge_mgas_out(i,2,9), bulge_mgas_out(i,2,10), &
			& bulge_mstars_in(i,2), bulge_mstars_out(i,2), bulge_mstars_formed(i,2), bulge_alpha(i,2), dt
		end do
		if( nbulge_in .gt. 0 ) then
			do j=1,nbulge_in
				m = bulge_in(j)
				if( is_list(m) .eq. 1 .and. is_redshift(m,1) .lt. bulge_redshift(nbulge,1) ) then
					if( nbulge_ex .gt. 0 ) then
						do l=1,nbulge_ex
							n = bulge_ex(l)
							if( es_list(n) .eq. 1 .and. es_redshift(n,1) .lt. bulge_redshift(nbulge,1) .and. es_redshift(n,1) .gt. is_redshift(m,1) ) then
								dt = 0.95_4/(((1+es_redshift(n,1))/7.0_4)**1.5) - time      !!! Gyr
								es_list(n) = 0

								write(20,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(n,1), es_redshift(n,1), es_rad(n,1), & 
								& es_gas_mass(n,1), es_star_mass(n,1), es_mass(n,1), es_gas_frac(n,1), es_dm_frac(n,1), es_gas_sig(n,1), & 
								& es_star_sig(n,1), es_sig(n,1), es_age(n,1), es_gas_met(n,1), es_star_met(n,1), es_SFR(n,1), es_sig_SFR(n,1), & 
								& es_sSFR(n,1), es_tau(n,1), es_d_Rd(n,1), es_z_Hd(n,1), es_d(n,1), es_z(n,1), es_residual(n,1), es_shape(n,1), & 
								& es_dm(n,1), es_es(n,1), es_merger(n,1), es_tff(n,1), es_td(n,1), es_td_global(n,1), es_mgas_in(n,1,1), es_mgas_in(n,1,2), & 
								& es_mgas_in(n,1,3), es_mgas_out(n,1,1), es_mgas_out(n,1,2), es_mgas_out(n,1,3), es_mgas_out(n,1,4), es_mgas_out(n,1,5), & 
								& es_mgas_out(n,1,6), es_mgas_out(n,1,7), es_mgas_out(n,1,8), es_mgas_out(n,1,9), es_mgas_out(n,1,10), &
								& es_mstars_in(n,1), es_mstars_out(n,1), es_mstars_formed(n,1), es_alpha(n,1), dt
								write(23,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(n,2), es_redshift(n,2), es_rad(n,2), & 
								& es_gas_mass(n,2), es_star_mass(n,2), es_mass(n,2), es_gas_frac(n,2), es_dm_frac(n,2), es_gas_sig(n,2), & 
								& es_star_sig(n,2), es_sig(n,2), es_age(n,2), es_gas_met(n,2), es_star_met(n,2), es_SFR(n,2), es_sig_SFR(n,2), & 
								& es_sSFR(n,2), es_tau(n,2), es_d_Rd(n,2), es_z_Hd(n,2), es_d(n,2), es_z(n,2), es_residual(n,2), es_shape(n,2), & 
								& es_dm(n,2), es_es(n,2), es_merger(n,2), es_tff(n,2), es_td(n,2), es_td_global(n,2), es_mgas_in(n,2,1), es_mgas_in(n,2,2), & 
								& es_mgas_in(n,2,3), es_mgas_out(n,2,1), es_mgas_out(n,2,2), es_mgas_out(n,2,3), es_mgas_out(n,2,4), es_mgas_out(n,2,5), & 
								& es_mgas_out(n,2,6), es_mgas_out(n,2,7), es_mgas_out(n,2,8), es_mgas_out(n,2,9), es_mgas_out(n,2,10), & 
								& es_mstars_in(n,2), es_mstars_out(n,2), es_mstars_formed(n,2), es_alpha(n,2), dt
							end if
						end do
					end if

					dt = 0.95_4/(((1+is_redshift(m,1))/7.0_4)**1.5) - time      !!! Gyr
					is_list(m) = 0

					write(20,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')is_id(m,1), is_redshift(m,1), is_rad(m,1), is_gas_mass(m,1), &
					& is_star_mass(m,1), is_mass(m,1), is_gas_frac(m,1), is_dm_frac(m,1), is_gas_sig(m,1), is_star_sig(m,1), is_sig(m,1), is_age(m,1), &
					& is_gas_met(m,1), is_star_met(m,1), is_SFR(m,1), is_sig_SFR(m,1), is_sSFR(m,1), is_tau(m,1), is_d_Rd(m,1), is_z_Hd(m,1), is_d(m,1), is_z(m,1), &
					& is_residual(m,1), is_shape(m,1), is_dm(m,1), is_es(m,1), is_merger(m,1), is_tff(m,1), is_td(m,1), is_td_global(m,1), &
					& is_mgas_in(m,1,1), is_mgas_in(m,1,2), is_mgas_in(m,1,3), is_mgas_out(m,1,1), is_mgas_out(m,1,2), &
					& is_mgas_out(m,1,3), is_mgas_out(m,1,4), is_mgas_out(m,1,5), is_mgas_out(m,1,6), is_mgas_out(m,1,7), &
					& is_mgas_out(m,1,8), is_mgas_out(m,1,9), is_mgas_out(m,1,10), &  
					& is_mstars_in(m,1), is_mstars_out(m,1), is_mstars_formed(m,1), is_alpha(m,1), dt
					write(23,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')is_id(m,2), is_redshift(m,2), is_rad(m,2), is_gas_mass(m,2), &
					& is_star_mass(m,2), is_mass(m,2), is_gas_frac(m,2), is_dm_frac(m,2), is_gas_sig(m,2), is_star_sig(m,2), is_sig(m,2), is_age(m,2), &
					& is_gas_met(m,2), is_star_met(m,2), is_SFR(m,2), is_sig_SFR(m,2), is_sSFR(m,2), is_tau(m,2), is_d_Rd(m,2), is_z_Hd(m,2), is_d(m,2), is_z(m,2), &
					& is_residual(m,2), is_shape(m,2), is_dm(m,2), is_es(m,2), is_merger(m,2), is_tff(m,2), is_td(m,2), is_td_global(m,2), &
					& is_mgas_in(m,2,1), is_mgas_in(m,2,2), is_mgas_in(m,2,3), is_mgas_out(m,2,1), is_mgas_out(m,2,2), &
					& is_mgas_out(m,2,3), is_mgas_out(m,2,4), is_mgas_out(m,2,5), is_mgas_out(m,2,6), is_mgas_out(m,2,7), &
					& is_mgas_out(m,2,8), is_mgas_out(m,2,9), is_mgas_out(m,2,10), &
					& is_mstars_in(m,2), is_mstars_out(m,2), is_mstars_formed(m,2), is_alpha(m,2), dt
				end if
			end do
		end if
		if( nbulge_ex .gt. 0 ) then
			do l=1,nbulge_ex
				n = bulge_ex(l)
				if( es_list(n) .eq. 1 .and. es_redshift(n,1) .lt. bulge_redshift(nbulge,1) ) then
					if( nbulge_in .gt. 0 ) then 
						if( es_redshift(n,1) .lt. minval(is_redshift(bulge_in(1:nbulge_in),1)) ) then
							dt = 0.95_4/(((1+es_redshift(n,1))/7.0_4)**1.5) - time      !!! Gyr
							es_list(n) = 0

							write(20,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(n,1), es_redshift(n,1), es_rad(n,1), es_gas_mass(n,1), &
							& es_star_mass(n,1), es_mass(n,1), es_gas_frac(n,1), es_dm_frac(n,1), es_gas_sig(n,1), es_star_sig(n,1), es_sig(n,1), es_age(n,1), &
							& es_gas_met(n,1), es_star_met(n,1), es_SFR(n,1), es_sig_SFR(n,1), es_sSFR(n,1), es_tau(n,1), es_d_Rd(n,1), es_z_Hd(n,1), es_d(n,1), &
							& es_z(n,1), es_residual(n,1), es_shape(n,1), es_dm(n,1), es_es(n,1), es_merger(n,1), es_tff(n,1), es_td(n,1), es_td_global(n,1), & 
							& es_mgas_in(n,1,1), es_mgas_in(n,1,2), es_mgas_in(n,1,3), es_mgas_out(n,1,1), es_mgas_out(n,1,2), & 
							& es_mgas_out(n,1,3), es_mgas_out(n,1,4), es_mgas_out(n,1,5), es_mgas_out(n,1,6), es_mgas_out(n,1,7), &
							& es_mgas_out(n,1,8), es_mgas_out(n,1,9), es_mgas_out(n,1,10), &
							& es_mstars_in(n,1), es_mstars_out(n,1), es_mstars_formed(n,1), es_alpha(n,1), dt
							write(23,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(n,2), es_redshift(n,2), es_rad(n,2), es_gas_mass(n,2), &
							& es_star_mass(n,2), es_mass(n,2), es_gas_frac(n,2), es_dm_frac(n,2), es_gas_sig(n,2), es_star_sig(n,2), es_sig(n,2), es_age(n,2), &
							& es_gas_met(n,2), es_star_met(n,2), es_SFR(n,2), es_sig_SFR(n,2), es_sSFR(n,2), es_tau(n,2), es_d_Rd(n,2), es_z_Hd(n,2), es_d(n,2), &
							& es_z(n,2), es_residual(n,2), es_shape(n,2), es_dm(n,2), es_es(n,2), es_merger(n,2), es_tff(n,2), es_td(n,2), es_td_global(n,2), & 
							& es_mgas_in(n,2,1), es_mgas_in(n,2,2), es_mgas_in(n,2,3), es_mgas_out(n,2,1), es_mgas_out(n,2,2), & 
							& es_mgas_out(n,2,3), es_mgas_out(n,2,4), es_mgas_out(n,2,5), es_mgas_out(n,2,6), es_mgas_out(n,2,7), &
							& es_mgas_out(n,2,8), es_mgas_out(n,2,9), es_mgas_out(n,2,10), &  
							& es_mstars_in(n,2), es_mstars_out(n,2), es_mstars_formed(n,2), es_alpha(n,2), dt
						end if
					else
						dt = 0.95_4/(((1+es_redshift(n,1))/7.0_4)**1.5) - time      !!! Gyr
						es_list(n) = 0

						write(20,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(n,1), es_redshift(n,1), es_rad(n,1), es_gas_mass(n,1), &
						& es_star_mass(n,1), es_mass(n,1), es_gas_frac(n,1), es_dm_frac(n,1), es_gas_sig(n,1), es_star_sig(n,1), es_sig(n,1), es_age(n,1), &
						& es_gas_met(n,1), es_star_met(n,1), es_SFR(n,1), es_sig_SFR(n,1), es_sSFR(n,1), es_tau(n,1), es_d_Rd(n,1), es_z_Hd(n,1), es_d(n,1), & 
						& es_z(n,1), es_residual(n,1), es_shape(n,1), es_dm(n,1), es_es(n,1), es_merger(n,1), es_tff(n,1), es_td(n,1), es_td_global(n,1), &
						& es_mgas_in(n,1,1), es_mgas_in(n,1,2), es_mgas_in(n,1,3), es_mgas_out(n,1,1), es_mgas_out(n,1,2), &
						& es_mgas_out(n,1,3), es_mgas_out(n,1,4), es_mgas_out(n,1,5), es_mgas_out(n,1,6), es_mgas_out(n,1,7), &
						& es_mgas_out(n,1,8), es_mgas_out(n,1,9), es_mgas_out(n,1,10), &
						& es_mstars_in(n,1), es_mstars_out(n,1), es_mstars_formed(n,1), es_alpha(n,1), dt
						write(23,'(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,21(1x,es10.2))')es_id(n,2), es_redshift(n,2), es_rad(n,2), es_gas_mass(n,2), &
						& es_star_mass(n,2), es_mass(n,2), es_gas_frac(n,2), es_dm_frac(n,2), es_gas_sig(n,2), es_star_sig(n,2), es_sig(n,2), es_age(n,2), &
						& es_gas_met(n,2), es_star_met(n,2), es_SFR(n,2), es_sig_SFR(n,2), es_sSFR(n,2), es_tau(n,2), es_d_Rd(n,2), es_z_Hd(n,2), es_d(n,2), & 
						& es_z(n,2), es_residual(n,2), es_shape(n,2), es_dm(n,2), es_es(n,2), es_merger(n,2), es_tff(n,2), es_td(n,2), es_td_global(n,2), & 
						& es_mgas_in(n,2,1), es_mgas_in(n,2,2), es_mgas_in(n,2,3), es_mgas_out(n,2,1), es_mgas_out(n,2,2), & 
						& es_mgas_out(n,2,3), es_mgas_out(n,2,4), es_mgas_out(n,2,5), es_mgas_out(n,2,6), es_mgas_out(n,2,7), &
						& es_mgas_out(n,2,8), es_mgas_out(n,2,9), es_mgas_out(n,2,10), &  
						& es_mstars_in(n,2), es_mstars_out(n,2), es_mstars_formed(n,2), es_alpha(n,2), dt
					end if
				end if
			end do
		end if

		!!! Check for errors !!!
		do i=1,nes
			if( es_list(i) .ne. 0) then
				print *, 'error in es_list 5', i
				stop
			end if
		end do          

		do i=1,nis
			if( is_list(i) .ne. 0) then
				print *, 'error in is_list 6', i
				stop
			end if
		end do          

		do i=1,nbulge
			if( bulge_list(i) .ne. 0) then
				print *, 'error in bulge_list 3', i
				stop
			end if
		end do
		!!!!!!!!!!!!!!!!!!!!!!!!!!

		close(unit=18)
		close(unit=19)
		close(unit=20)
		close(unit=21)
		close(unit=22)
		close(unit=23)
		if(nis .gt. 0) then
			deallocate( is_redshift, is_rad, is_gas_mass, is_star_mass, is_mass, is_gas_frac, is_dm_frac, is_gas_sig, is_star_sig, is_sig, is_age, is_gas_met, is_star_met, & 
			& is_SFR, is_sig_SFR, is_sSFR, is_tau, is_d_Rd, is_z_Hd, is_d, is_z, is_residual, is_shape, is_dm, is_tff, is_td, is_td_global, is_mgas_in, is_mgas_out, & 
			& is_mstars_in, is_mstars_out, is_mstars_formed, is_alpha, is_es, is_merger, is_id, is_list )
		else
			deallocate( is_id, is_redshift, is_list )
		end if

		if(nes .gt. 0) then
			deallocate( es_redshift, es_rad, es_gas_mass, es_star_mass, es_mass, es_gas_frac, es_dm_frac, es_gas_sig, es_star_sig, es_sig, es_age, es_gas_met, es_star_met, & 
			& es_SFR, es_sig_SFR, es_sSFR, es_tau, es_d_Rd, es_z_Hd, es_d, es_z, es_residual, es_shape, es_dm, es_tff, es_td, es_td_global, es_mgas_in, es_mgas_out, & 
			& es_mstars_in, es_mstars_out, es_mstars_formed, es_alpha, es_es, es_merger, es_id, es_list )
		else
			deallocate( es_id, es_redshift, es_list )
		end if

		if(nbulge .gt. 0) then
			deallocate( bulge_redshift, bulge_rad, bulge_gas_mass, bulge_star_mass, bulge_mass, bulge_gas_frac, bulge_dm_frac, bulge_gas_sig, bulge_star_sig, bulge_sig, &
			& bulge_age, bulge_gas_met, bulge_star_met, bulge_SFR, bulge_sig_SFR, bulge_sSFR, bulge_tau, bulge_d_Rd, bulge_z_Hd, bulge_d, bulge_z, bulge_residual, & 
			& bulge_shape, bulge_dm, bulge_tff, bulge_td, bulge_td_global, bulge_mgas_in, bulge_mgas_out, bulge_mstars_in, bulge_mstars_out, bulge_mstars_formed, &
			& bulge_alpha, bulge_es, bulge_merger, bulge_id, bulge_list )
		else
			deallocate( bulge_id, bulge_redshift, bulge_list )
		end if

            deallocate( ex_in, bulge_in, bulge_ex )
	end do
	deallocate( gal_name )
end program main




