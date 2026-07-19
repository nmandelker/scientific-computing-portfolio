program main
	implicit none
	character(len=20),allocatable :: gal_name(:)                                !!! This is the name of the galaxy you want to search for clumps. Can have either 3 characters (MW3) or 4 characters (VL01, SFG1, MW10)
	character(len=256) :: filename, path_name, input_arg, dirname
	integer :: i, j, k, Nsnapshot, Nsimulation
	real(4) :: wide, dens_thresh
	real(4),allocatable :: Rvir(:), aexp(:), Mvir(:), Vvir(:), redshift(:)

	real(4) :: Ld1(3), Ld2(3), Ld3(3), Rd, Hd, Rc, rcar_clump(3), vcar_clump(3), clump_cen(3), vcyl_clump(3), SFR, age, metg, mets
	real(8) :: Mgc, Msc, Mdmc, Mnc
	integer :: nis, nes, es, ncell, nstar, ind, new, Ngroup, Ntot(6), IDc, Merger
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
		Nsnapshot = 0

		!!!!!!!!!! How many snapshots in the simulation? !!!!!!!!!!
		write(filename,'(a,a,a)') '/BIGDATA/nirm/thin_time_slices/galaxy_inputs/',trim(gal_name(j)),'/galaxy_catalogue/Nir_halo_cat.txt'
		open(unit=16,file=filename,form='formatted')
		read(16,'(1x,i4)') Nsnapshot
		print *, Nsnapshot
		call allocate_global(Nsnapshot)
		do i=1,Nsnapshot
			read(16,'(1x,f6.4,3(1x,es12.5))') aexp(i), Rvir(i), Mvir(i), Vvir(i)
		end do
		close(unit=16)
		redshift(1:Nsnapshot) = 1/aexp(1:Nsnapshot) - 1
		call open_files()

		do k=2,Nsnapshot
			print *, gal_name(j), k, aexp(k)
			Ntot(:) = 0
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

!				print ('(3(1x,e12.5))'), Ld1(:)
!				print ('(3(1x,e12.5))'), Ld2(:)
!				print ('(3(1x,e12.5))'), Ld3(:)
!				print ('(2(1x,f7.3))'), Rd, Hd
!				print ('(3(1x,i3))'), Ngroup, nis, nes

				if( Ngroup .gt. 0 ) then
					do i=1,Ngroup
						read(24,'(1x,a4,2(1x,i5.5),1x,i5,1x,i1,3(1x,f7.3),3(1xes12.5),4(1xf7.3),7(1x,es12.5),1x,i5,1x,i10,2(1x,es12.5),2(1x,f7.3),1x,i3)') &
						& comm, IDc, Merger, ind, new, rcar_clump(:), vcar_clump(:), Rc, clump_cen(:), vcyl_clump(:), Mgc, Msc, Mdmc, Mnc, ncell, nstar, SFR, age, metg, mets, es

!						print ('(1x,a4,2(1x,i5.5),1x,i5,1x,i1,3(1x,f7.3),3(1xes12.5),4(1xf7.3),7(1x,es12.5),1x,i5,1x,i10,2(1x,es12.5),2(1x,f7.3),1x,i3)'), &
!						& comm, IDc, Merger, ind, new, rcar_clump(:), vcar_clump(:), Rc, clump_cen(:), vcyl_clump(:), Mgc, Msc, Mdmc, Mnc, ncell, nstar, SFR, age, metg, mets, es

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

							if( log10(Mgc + Msc) .le. 6.5 ) then
								write(18,'(1x,f6.4,1x,i5.5,3(1x,f7.3),4(1xes12.5),1x,f7.3,1x,i5,1x,i10,4(1x,es12.5),1x,f6.4,6(1x,f7.3))') aexp(k), IDc, xnew, ynew, znew, vr, omega*rnew, vz, delt, Rc, ncell, nstar, & 
								& Mgc, Msc, SFR, age, aexp(k-1), xold(1), yold(1), zold, xold(2), yold(2), zold
								Ntot(1) = Ntot(1) + 1
							elseif( log10(Mgc + Msc) .le. 7.0 ) then
								write(19,'(1x,f6.4,1x,i5.5,3(1x,f7.3),4(1xes12.5),1x,f7.3,1x,i5,1x,i10,4(1x,es12.5),1x,f6.4,6(1x,f7.3))') aexp(k), IDc, xnew, ynew, znew, vr, omega*rnew, vz, delt, Rc, ncell, nstar, & 
								& Mgc, Msc, SFR, age, aexp(k-1), xold(1), yold(1), zold, xold(2), yold(2), zold
								Ntot(2) = Ntot(2) + 1
							elseif( log10(Mgc + Msc) .le. 7.5 ) then
								write(20,'(1x,f6.4,1x,i5.5,3(1x,f7.3),4(1xes12.5),1x,f7.3,1x,i5,1x,i10,4(1x,es12.5),1x,f6.4,6(1x,f7.3))') aexp(k), IDc, xnew, ynew, znew, vr, omega*rnew, vz, delt, Rc, ncell, nstar, & 
								& Mgc, Msc, SFR, age, aexp(k-1), xold(1), yold(1), zold, xold(2), yold(2), zold
								Ntot(3) = Ntot(3) + 1
							elseif( log10(Mgc + Msc) .le. 8.0 ) then
								write(21,'(1x,f6.4,1x,i5.5,3(1x,f7.3),4(1xes12.5),1x,f7.3,1x,i5,1x,i10,4(1x,es12.5),1x,f6.4,6(1x,f7.3))') aexp(k), IDc, xnew, ynew, znew, vr, omega*rnew, vz, delt, Rc, ncell, nstar, & 
								& Mgc, Msc, SFR, age, aexp(k-1), xold(1), yold(1), zold, xold(2), yold(2), zold
								Ntot(4) = Ntot(4) + 1
							elseif( log10(Mgc + Msc) .le. 8.5 ) then
								write(22,'(1x,f6.4,1x,i5.5,3(1x,f7.3),4(1xes12.5),1x,f7.3,1x,i5,1x,i10,4(1x,es12.5),1x,f6.4,6(1x,f7.3))') aexp(k), IDc, xnew, ynew, znew, vr, omega*rnew, vz, delt, Rc, ncell, nstar, & 
								& Mgc, Msc, SFR, age, aexp(k-1), xold(1), yold(1), zold, xold(2), yold(2), zold
								Ntot(5) = Ntot(5) + 1
							else
								write(23,'(1x,f6.4,1x,i5.5,3(1x,f7.3),4(1xes12.5),1x,f7.3,1x,i5,1x,i10,4(1x,es12.5),1x,f6.4,6(1x,f7.3))') aexp(k), IDc, xnew, ynew, znew, vr, omega*rnew, vz, delt, Rc, ncell, nstar, & 
								& Mgc, Msc, SFR, age, aexp(k-1), xold(1), yold(1), zold, xold(2), yold(2), zold
								Ntot(6) = Ntot(6) + 1
							end if
						end if
					end do
					print *, 'nclump', Ngroup
					print *, 'Nnew', sum(Ntot(1:6))
					print *, 'by mass', Ntot(1:6)
				else
					print *, 'No clumps'
				end if
				close(unit=24)
			else
				print *, 'no file'
			end if
			print *, ''
		end do
		print *, 'done with', trim(gal_name(j))
		print *, ''
		call deallocate_global()
		call close_files()
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

		write(filename,'(a,a,a)') './',trim(dirname),'/proto_clumps/log_Mc_60_65.txt'
		open(unit=18,file=filename,form='formatted')

		write(filename,'(a,a,a)') './',trim(dirname),'/proto_clumps/log_Mc_65_70.txt'
		open(unit=19,file=filename,form='formatted')

		write(filename,'(a,a,a)') './',trim(dirname),'/proto_clumps/log_Mc_70_75.txt'
		open(unit=20,file=filename,form='formatted')

		write(filename,'(a,a,a)') './',trim(dirname),'/proto_clumps/log_Mc_75_80.txt'
		open(unit=21,file=filename,form='formatted')

		write(filename,'(a,a,a)') './',trim(dirname),'/proto_clumps/log_Mc_80_85.txt'
		open(unit=22,file=filename,form='formatted')

		write(filename,'(a,a,a)') './',trim(dirname),'/proto_clumps/log_Mc_85_00.txt'
		open(unit=23,file=filename,form='formatted')

		write(18, *) 'a1     ID      x       y       z      vr           vphi         vz           delta t       Rc     Ncells     N*   Mgc          Msc          SFRc          agec          a2       x1      y1      z1      x2      y2      z2'
		write(19, *) 'a1     ID      x       y       z      vr           vphi         vz           delta t       Rc     Ncells     N*   Mgc          Msc          SFRc          agec          a2       x1      y1      z1      x2      y2      z2'
		write(20, *) 'a1     ID      x       y       z      vr           vphi         vz           delta t       Rc     Ncells     N*   Mgc          Msc          SFRc          agec          a2       x1      y1      z1      x2      y2      z2'
		write(21, *) 'a1     ID      x       y       z      vr           vphi         vz           delta t       Rc     Ncells     N*   Mgc          Msc          SFRc          agec          a2       x1      y1      z1      x2      y2      z2'
		write(22, *) 'a1     ID      x       y       z      vr           vphi         vz           delta t       Rc     Ncells     N*   Mgc          Msc          SFRc          agec          a2       x1      y1      z1      x2      y2      z2'
		write(23, *) 'a1     ID      x       y       z      vr           vphi         vz           delta t       Rc     Ncells     N*   Mgc          Msc          SFRc          agec          a2       x1      y1      z1      x2      y2      z2'

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

