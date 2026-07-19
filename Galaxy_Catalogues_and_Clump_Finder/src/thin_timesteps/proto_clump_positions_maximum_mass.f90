program main
	implicit none
	character(len=20),allocatable :: gal_name(:)                                !!! This is the name of the galaxy you want to search for clumps. Can have either 3 characters (MW3) or 4 characters (VL01, SFG1, MW10)
	character(len=256) :: filename, path_name, input_arg, dirname
	integer :: i, j, k, m, n, Nsnapshot, Nsimulation
	real(4) :: wide, dens_thresh
	real(4),allocatable :: Rvir(:), aexp(:), Mvir(:), Vvir(:), redshift(:)

	real(4),allocatable :: is_aexp(:), is_redshift(:), is_rad(:), is_gas_mass(:), is_star_mass(:), is_mass(:)
	real(4),allocatable :: is_gas_frac(:), is_dm_frac(:), is_gas_sig(:), is_star_sig(:), is_sig(:)
	real(4),allocatable :: is_age(:), is_gas_met(:), is_star_met(:), is_SFR(:), is_sig_SFR(:), is_sSFR(:), is_tau(:)
	real(4),allocatable :: is_d_Rd(:), is_z_Hd(:), is_d(:), is_z(:), is_residual(:), is_shape(:,:), is_dm(:)
	real(4),allocatable :: is_es(:), is_tff(:), is_td(:), is_td_global(:), is_mgas_in(:,:), is_mgas_out(:,:), is_mstars_in(:), is_mstars_out(:), is_mstars_formed(:), is_delt(:)
	integer,allocatable :: is_id(:), is_snap(:)

	real(4) :: Ld1(3), Ld2(3), Ld3(3), Rd, Hd, Rc, rcar_clump(3), vcar_clump(3), clump_cen(3), vcyl_clump(3), SFR, age, metg, mets, delt_clump_max
	real(8) :: Mgc, Msc, Mdmc, Mnc, Mclump_max
	integer :: nis, nes, es, ncell, nstar, ind, new, Ngroup, IDc, Merger, Nis_tot, aexp_loc(1), Nsnap_clump
	character(len=4) :: comm
	logical :: file_exist
	real(4) :: delt, xnew, ynew, znew, rnew, phinew, xold(2), yold(2), zold, rold, phiold, vx, vy, vz, vr, omega
	real(4),parameter :: pi = 3.141592654_4, pi2 = 2.0_4*pi, pi4_3 = (4.0_4 / 3.0_4)*pi

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
	write(filename,'(a)') './protoclump_input.dat'
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
	do j=1,Nsimulation
		print *, 'Simulation:', gal_name(j)
		write(dirname,'(a,a,i2.2,a,f4.2)') trim(gal_name(j)),'_thresh_',floor(dens_thresh),'_wide_FWHM_',wide/1000.0_4

		write(filename,'(a,a,a)') './',trim(dirname),'/Matlab_friendly_data/same_clumps/in_situ.out'
		open(unit=12,file=filename,status="old",action="read",form='formatted')
		Nis_tot = 0
		do
			read(12,*,end=6)
			Nis_tot = Nis_tot + 1
		end do
 6		rewind(unit=12)
		print *, 'Nclumps:',Nis_tot
		if( Nis_tot .gt. 0 ) then
			allocate( is_aexp(Nis_tot), is_redshift(Nis_tot), is_rad(Nis_tot), is_gas_mass(Nis_tot), is_star_mass(Nis_tot), is_mass(Nis_tot) )
			allocate( is_gas_frac(Nis_tot), is_dm_frac(Nis_tot), is_gas_sig(Nis_tot), is_star_sig(Nis_tot), is_sig(Nis_tot) )
			allocate( is_age(Nis_tot), is_gas_met(Nis_tot), is_star_met(Nis_tot), is_SFR(Nis_tot), is_sig_SFR(Nis_tot), is_sSFR(Nis_tot), is_tau(Nis_tot) )
			allocate( is_d_Rd(Nis_tot), is_z_Hd(Nis_tot), is_d(Nis_tot), is_z(Nis_tot), is_residual(Nis_tot), is_shape(Nis_tot,3), is_dm(Nis_tot) )
			allocate( is_es(Nis_tot), is_tff(Nis_tot), is_td(Nis_tot), is_td_global(Nis_tot), is_mgas_in(Nis_tot,3), is_mgas_out(Nis_tot,3), is_mstars_in(Nis_tot), is_mstars_out(Nis_tot), is_mstars_formed(Nis_tot), is_delt(nis_tot) )
			allocate( is_id(Nis_tot), is_snap(Nis_tot) )
			do i = 1, Nis_tot
				read(12,'(1x,i5.5,1x,f6.3,25(1x,es10.2),1x,i3,13(1x,es10.2))') is_id(i), is_redshift(i), is_rad(i), is_gas_mass(i), &
				& is_star_mass(i), is_mass(i), is_gas_frac(i), is_dm_frac(i), is_gas_sig(i), is_star_sig(i), is_sig(i), is_age(i), &
				& is_gas_met(i), is_star_met(i), is_SFR(i), is_sig_SFR(i), is_sSFR(i), is_tau(i), is_d_Rd(i), is_z_Hd(i), is_d(i), is_z(i), &
				& is_residual(i), is_shape(i,1), is_shape(i,2), is_shape(i,3), is_dm(i), is_es(i), is_tff(i), is_td(i), is_td_global(i), &
				& is_mgas_in(i,1), is_mgas_in(i,2), is_mgas_in(i,3), is_mgas_out(i,1), is_mgas_out(i,2), is_mgas_out(i,3), is_mstars_in(i), is_mstars_out(i), is_mstars_formed(i), is_delt(i)
			end do
			close(unit=12)

			!!!!!!!!!! How many snapshots in the simulation? !!!!!!!!!!
			Nsnapshot = 0
			write(filename,'(a,a,a)') '/BIGDATA/nirm/thin_time_slices/galaxy_inputs/',trim(gal_name(j)),'/galaxy_catalogue/Nir_halo_cat.txt'
			open(unit=16,file=filename,form='formatted')
			read(16,'(1x,i4)') Nsnapshot
			print *, 'Nsnapshot:',Nsnapshot
			call allocate_global(Nsnapshot)
			do i=1,Nsnapshot
				read(16,'(1x,f6.4,3(1x,es12.5))') aexp(i), Rvir(i), Mvir(i), Vvir(i)
			end do
			close(unit=16)
			redshift(1:Nsnapshot) = 1/aexp(1:Nsnapshot) - 1
			call open_files()

			do i=1,Nis_tot
				aexp_loc(:) = minloc( abs(is_redshift(i)-redshift(1:Nsnapshot)) )
				is_aexp(i) = aexp(aexp_loc(1))
				is_snap(i) = aexp_loc(1)
			end do

			n = 1
			do while( n .le. Nis_tot )
				print *, 'clump',n,'of',Nis_tot
				if( is_delt(n) .eq. 0.0_4 ) then		!!! Looks like a "new" clump
					Mclump_max = is_mass(n)
					Nsnap_clump = 1
					delt_clump_max = 0.0_4
					m = n+1
					do while( is_delt(m) .gt. 0.0_4 )
						Nsnap_clump = Nsnap_clump + 1
						delt_clump_max = is_delt(m)
						Mclump_max = max( Mclump_max, is_mass(m) )
						m = m+1
					end do					!!! At the end, m is the index of the next new clump, where m>n and delt(m)=0
					if( is_snap(n) .ge. 2 ) then		!!! look from second snapshot
						k = is_snap(n)
						delt = 0.95_4/(((1+redshift(k))/7.0_4)**1.5) - 0.95_4/(((1+redshift(k-1))/7.0_4)**1.5)	!!! Gyr
						write(filename,'(a,a,a,f6.4,a)') './',trim(dirname),'/clump_catalogue/Nir_clump_cat_a',aexp(k),'.txt'
						inquire(file=filename, exist=file_exist)
						if(file_exist) then
							open(unit=24,file=filename,form='formatted')
							read(24,'(3(1x,e12.5))') Ld1(:)
							read(24,'(3(1x,e12.5))') Ld2(:)
							read(24,'(3(1x,e12.5))') Ld3(:)
							read(24,'(2(1x,f7.3))') Rd, Hd
							read(24,'(3(1x,i3))') Ngroup, nis, nes
!							print ('(3(1x,e12.5))'), Ld1(:)
!							print ('(3(1x,e12.5))'), Ld2(:)
!							print ('(3(1x,e12.5))'), Ld3(:)
!							print ('(2(1x,f7.3))'), Rd, Hd
!							print ('(3(1x,i3))'), Ngroup, nis, nes
							if( Ngroup .gt. 0 ) then
								i = 1
								do while (i .le. Ngroup)
									read(24,'(1x,a4,2(1x,i5.5),1x,i5,1x,i1,3(1x,f7.3),3(1xes12.5),4(1xf7.3),7(1x,es12.5),1x,i5,1x,i10,2(1x,es12.5),2(1x,f7.3),1x,i3)') &
									& comm, IDc, Merger, ind, new, rcar_clump(:), vcar_clump(:), Rc, clump_cen(:), vcyl_clump(:), Mgc, Msc, Mdmc, Mnc, ncell, nstar, SFR, age, metg, mets, es
									if( IDc .eq. is_id(n) ) then
										if( new .eq. 1 .and. trim(comm) .eq. 'Is' ) then
											xnew = clump_cen(1)
											ynew = clump_cen(2)
											znew = clump_cen(3)
											rnew = sqrt( xnew**2 + ynew**2 )
											if(xnew.eq.0.0_4.and.ynew.eq.0.0_4) then
												phinew = 0.0_4
											elseif(xnew>=0.0_4.and.ynew>=0.0_4) then
												phinew = atan(ynew/xnew)
											elseif(xnew<0.0_4) then
												phinew = pi+atan(ynew/xnew)
											else
												phinew = pi2+atan(ynew/xnew)
											end if
											vr    = 1.023_4 * vcyl_clump(1)
											omega = 1.023_4 * vcyl_clump(2)/rnew
											vz    = 1.023_4 * vcyl_clump(3)
											vx    = 1.023_4 * ( vr*cos(phinew) - vcyl_clump(2)*sin(phinew) )
											vy    = 1.023_4 * ( vr*sin(phinew) + vcyl_clump(2)*cos(phinew) )

											xold(1) = xnew   - vx   *delt
											yold(1) = ynew   - vy   *delt
											zold    = znew   - vz   *delt
											rold    = rnew   - vr   *delt
											phiold  = phinew - omega*delt
											xold(2) = rold*cos(phiold)
											yold(2) = rold*sin(phiold)

											if( log10(Mclump_max) .le. 6.5 ) then
												write(18,'(1x,f6.4,1x,i5.5,3(1x,f7.3),4(1xes12.5),1x,f7.3,1x,i5,1x,i10,6(1x,es12.5),1x,i4,3x,f6.4,6(1x,f7.3))') aexp(k), IDc, xnew, ynew, znew, vr, omega*rnew, vz, &
												& delt, Rc, ncell, nstar, Mgc, Msc, SFR, age, Mclump_max, delt_clump_max, Nsnap_clump, aexp(k-1), xold(1), yold(1), zold, xold(2), yold(2), zold
											elseif( log10(Mclump_max) .le. 7.0 ) then
												write(19,'(1x,f6.4,1x,i5.5,3(1x,f7.3),4(1xes12.5),1x,f7.3,1x,i5,1x,i10,6(1x,es12.5),1x,i4,3x,f6.4,6(1x,f7.3))') aexp(k), IDc, xnew, ynew, znew, vr, omega*rnew, vz, &
												& delt, Rc, ncell, nstar, Mgc, Msc, SFR, age, Mclump_max, delt_clump_max, Nsnap_clump, aexp(k-1), xold(1), yold(1), zold, xold(2), yold(2), zold
											elseif( log10(Mclump_max) .le. 7.5 ) then
												write(20,'(1x,f6.4,1x,i5.5,3(1x,f7.3),4(1xes12.5),1x,f7.3,1x,i5,1x,i10,6(1x,es12.5),1x,i4,3x,f6.4,6(1x,f7.3))') aexp(k), IDc, xnew, ynew, znew, vr, omega*rnew, vz, &
												& delt, Rc, ncell, nstar, Mgc, Msc, SFR, age, Mclump_max, delt_clump_max, Nsnap_clump, aexp(k-1), xold(1), yold(1), zold, xold(2), yold(2), zold
											elseif( log10(Mclump_max) .le. 8.0 ) then
												write(21,'(1x,f6.4,1x,i5.5,3(1x,f7.3),4(1xes12.5),1x,f7.3,1x,i5,1x,i10,6(1x,es12.5),1x,i4,3x,f6.4,6(1x,f7.3))') aexp(k), IDc, xnew, ynew, znew, vr, omega*rnew, vz, &
												& delt, Rc, ncell, nstar, Mgc, Msc, SFR, age, Mclump_max, delt_clump_max, Nsnap_clump, aexp(k-1), xold(1), yold(1), zold, xold(2), yold(2), zold
											elseif( log10(Mclump_max) .le. 8.5 ) then
												write(22,'(1x,f6.4,1x,i5.5,3(1x,f7.3),4(1xes12.5),1x,f7.3,1x,i5,1x,i10,6(1x,es12.5),1x,i4,3x,f6.4,6(1x,f7.3))') aexp(k), IDc, xnew, ynew, znew, vr, omega*rnew, vz, &
												& delt, Rc, ncell, nstar, Mgc, Msc, SFR, age, Mclump_max, delt_clump_max, Nsnap_clump, aexp(k-1), xold(1), yold(1), zold, xold(2), yold(2), zold
											else
												write(23,'(1x,f6.4,1x,i5.5,3(1x,f7.3),4(1xes12.5),1x,f7.3,1x,i5,1x,i10,6(1x,es12.5),1x,i4,3x,f6.4,6(1x,f7.3))') aexp(k), IDc, xnew, ynew, znew, vr, omega*rnew, vz, &
												& delt, Rc, ncell, nstar, Mgc, Msc, SFR, age, Mclump_max, delt_clump_max, Nsnap_clump, aexp(k-1), xold(1), yold(1), zold, xold(2), yold(2), zold
											end if
										else
											print *, 'There is a problem'
											print *, 'The Matlab files say an in situ clump with delt=0'
											print *, 'but the catalogues have the clump as old or not in situ'
											print *, 'snapshot number, expansion factor, clump ID, delta t, new, comm'
											print *, is_snap(n), is_aexp(n), is_id(n), is_delt(n), new, trim(comm)
											stop
										end if
										i = 2*Ngroup + 5
									else
										i = i + 1
									end if
								end do
								if( i .eq. Ngroup + 1 ) then
									print *, 'There is a problem'
									print *, 'The clump does not seem to be in the catalogue matching aexpn taken from the Matlab file'
									print *, 'snapshot number, expansion factor, redshift, clump ID, delta t, Ngroup'
									print *, is_snap(n), is_aexp(n), is_redshift(n), is_id(n), is_delt(n), Ngroup
									stop
								end if
							else
								print *, 'There is a problem'
								print *, 'The clump from the Matlab file points to an empty catalogue'
								print *, 'snapshot number, expansion factor, redshift, clump ID, delta t, Ngroup'
								print *, is_snap(n), is_aexp(n), is_redshift(n), is_id(n), is_delt(n), Ngroup
								stop
							end if
							close(unit=24)
						else
							print *, 'There is a problem'
							print *, 'The clump from the Matlab file points to a catalogue that does not exist'
							print *, 'snapshot number, expansion factor, redshift, clump ID, delta t'
							print *, is_snap(n), is_aexp(n), is_redshift(n), is_id(n), is_delt(n)
							stop
						end if
					end if
					n = m					!!! At the end, m is the index of the next new clump, since delt(m)=0
				else
					print *, 'Something is wrong'
					print *, 'You are getting is_delt>0, when that shouldnt happen'
					print *, 'snapshot number, expansion factor, clump ID, delt'
					print *, is_snap(n), is_aexp(n), is_id(n), is_delt(n)
					stop
				end if
			end do
			deallocate( is_aexp, is_redshift, is_rad, is_gas_mass, is_star_mass, is_mass, is_gas_frac, is_dm_frac, is_gas_sig, is_star_sig, is_sig, is_age, is_gas_met, is_star_met, is_SFR, is_sig_SFR, is_sSFR, is_tau, &
			& is_d_Rd, is_z_Hd, is_d, is_z, is_residual, is_shape, is_dm, is_es, is_tff, is_td, is_td_global, is_mgas_in, is_mgas_out, is_mstars_in, is_mstars_out, is_mstars_formed, is_delt, is_id, is_snap )
			call deallocate_global()
			call close_files()
			print *, 'done with', trim(gal_name(j))
			print *, ''
		else
			print *, 'What do you want from me?'
			print *, 'There are no clumps in the Matlab file...'
			print *, ''
		end if
	end do
	deallocate( gal_name )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contains
	subroutine allocate_global(N)
!!!!!!!!!! Allocate global arrays per simulation !!!!!!!!!!
!!!!!!!!!! Deallocated in 'deallocate_global' !!!!!!!!!!
		implicit none
		integer,intent(in) :: N
		integer :: i

		allocate( aexp(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation aexp. stat= ', i
			stop
		end if
		allocate( Rvir(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Rvir. stat= ', i
			stop
		end if
		allocate( Mvir(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Mvir. stat= ', i
			stop
		end if
		allocate( Vvir(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Vvir. stat= ', i
			stop
		end if
		allocate( redshift(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation redshift. stat= ', i
			stop
		end if

	end subroutine allocate_global
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine deallocate_global()
!!! Deallocates all arrays from alloctae_global !!!
		implicit none

		if(allocated(aexp)) deallocate(aexp)
		if(allocated(Rvir)) deallocate(Rvir)
		if(allocated(Mvir)) deallocate(Mvir)
		if(allocated(Vvir)) deallocate(Vvir)
		if(allocated(redshift)) deallocate(redshift)

	end subroutine deallocate_global
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine open_files()
		implicit none

		write(filename,'(a,a,a)') 'mkdir -p ./',trim(dirname),'/proto_clumps'
		call system(trim(filename))

		write(filename,'(a,a,a)') 'mkdir -p ./',trim(dirname),'/proto_clumps/max_mass_ranking'
		call system(trim(filename))

		write(filename,'(a,a,a)') './',trim(dirname),'/proto_clumps/max_mass_ranking/log_Mc_60_65.txt'
		open(unit=18,file=filename,form='formatted')

		write(filename,'(a,a,a)') './',trim(dirname),'/proto_clumps/max_mass_ranking/log_Mc_65_70.txt'
		open(unit=19,file=filename,form='formatted')

		write(filename,'(a,a,a)') './',trim(dirname),'/proto_clumps/max_mass_ranking/log_Mc_70_75.txt'
		open(unit=20,file=filename,form='formatted')

		write(filename,'(a,a,a)') './',trim(dirname),'/proto_clumps/max_mass_ranking/log_Mc_75_80.txt'
		open(unit=21,file=filename,form='formatted')

		write(filename,'(a,a,a)') './',trim(dirname),'/proto_clumps/max_mass_ranking/log_Mc_80_85.txt'
		open(unit=22,file=filename,form='formatted')

		write(filename,'(a,a,a)') './',trim(dirname),'/proto_clumps/max_mass_ranking/log_Mc_85_00.txt'
		open(unit=23,file=filename,form='formatted')

		write(18, *) 'a1     ID      x       y       z      vr           vphi         vz           delta t       Rc     Ncells      N*   Mgc          Msc          SFRc         agec        max Mc       tmax        Nsnap  a2       x1     y1       z1      x2     y2       z2'
		write(19, *) 'a1     ID      x       y       z      vr           vphi         vz           delta t       Rc     Ncells      N*   Mgc          Msc          SFRc         agec        max Mc       tmax        Nsnap  a2       x1     y1       z1      x2     y2       z2'
		write(20, *) 'a1     ID      x       y       z      vr           vphi         vz           delta t       Rc     Ncells      N*   Mgc          Msc          SFRc         agec        max Mc       tmax        Nsnap  a2       x1     y1       z1      x2     y2       z2'
		write(21, *) 'a1     ID      x       y       z      vr           vphi         vz           delta t       Rc     Ncells      N*   Mgc          Msc          SFRc         agec        max Mc       tmax        Nsnap  a2       x1     y1       z1      x2     y2       z2'
		write(22, *) 'a1     ID      x       y       z      vr           vphi         vz           delta t       Rc     Ncells      N*   Mgc          Msc          SFRc         agec        max Mc       tmax        Nsnap  a2       x1     y1       z1      x2     y2       z2'
		write(23, *) 'a1     ID      x       y       z      vr           vphi         vz           delta t       Rc     Ncells      N*   Mgc          Msc          SFRc         agec        max Mc       tmax        Nsnap  a2       x1     y1       z1      x2     y2       z2'

	end subroutine open_files
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine close_files()
		implicit none

		close(unit=18)
		close(unit=19)
		close(unit=20)
		close(unit=21)
		close(unit=22)
		close(unit=23)

	end subroutine close_files

end program main

