module parameters
!!! Parameters used in the code, which can be varied to change performance, memory usage, definitions, implimentations, etc
	implicit none

	integer,parameter :: Ngrid_max = 600							!!! Maximum number of cells per axes in the 2D surface maps
	real(4),parameter :: narrow(5) = (/ 0.067_4, 0.134_4, 0.268_4, 0.536_4, 1.072_4 /)	!!! Narrow gaussian FWHM in physical kpc (=2 cells if Ngrid=600)
	real(4),parameter :: cool_gas = 1.5e4							!!! Maximum temperature (in K) for "cool gas" map
	real(4),parameter :: hot_gas = 1.0e6							!!! Minimum temperature (in K) for "hot gas" map
	real(4),parameter :: young_stars = 0.1_4						!!! Maximum age (in Gyr) for "young stars" map
	real(4),parameter :: old_stars = 0.5_4							!!! Minimum age (in Gyr) for "old stars" map
	real(4),parameter :: grid_size(5) = (/ 10.0_4, 20.0_4, 40.0_4, 80.0_4, 160.0_4 /)	!!! (in kpc) 2D maps go from -grid_size to +grid_size
	logical,parameter :: include_gas = .true.
	logical,parameter :: include_stars = .true.
	logical,parameter :: include_dm = .true.
	logical,parameter :: move_center = .false.
	real(4),parameter :: pi = 3.141592654_4, pi2 = 2.0_4*pi, pi4_3 = (4.0_4 / 3.0_4)*pi
end module parameters
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module globvar
!!! Global variables and allocatable arrays used throughout the code
	implicit none

!!!	Arrays for gas, stars and DM data
!!!	---------------------------------
	real(4),allocatable :: xgas(:), ygas(:), zgas(:), vxgas(:), vygas(:), vzgas(:)
	real(4),allocatable :: density_gas(:), temperature_gas(:), cell_size_gas(:), SNI_gas(:), SNII_gas(:)
	integer,allocatable :: Ngas(:)

	real(8),allocatable :: xstars(:), ystars(:), zstars(:), vxstars(:), vystars(:), vzstars(:), mass_stars(:)
	real(4),allocatable :: age_stars(:)
	integer,allocatable :: idstars(:)
	integer,allocatable :: Nstars(:)

	real(8),allocatable :: xdm(:), ydm(:), zdm(:), vxdm(:), vydm(:), vzdm(:), mass_dm(:)
	integer,allocatable :: iddm(:)
	integer,allocatable :: Ndm(:)

	real(4) :: res								      !!! resolution of the uni-grid

!!!	Arrays for disc data to be read from input files
!!!	---------------------------------
	real(4),allocatable :: Rvir(:), aexp(:), redshift(:), rcom(:,:), vcom(:,:), Ldisc(:,:), Lmag(:), Rdisc(:), Hdisc(:)
	real(8),allocatable :: Mvir(:), Vvir(:), Mgas_disc(:), Mcold_disc(:), Mstar_disc(:), M_Es_star_disc(:), Mdm_disc(:)
	real(8),allocatable :: SFR_disc(:), age_disc(:), metgas_disc(:), metstars_disc(:)

end module globvar
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module read_binary
!!! Reads post-treated ART files of gas, stars and dark matter data
use parameters
use globvar
	implicit none

contains
	subroutine data_gas(filename, nstop)
		implicit none
		character (*),intent (in) :: filename
		integer,intent(in) :: nstop
		integer :: i

		open ( 12 , file = filename, form = 'unformatted' )
		cell_size_gas(:) = 1.d-10
		xgas(:) = 1.d-10
		ygas(:) = 1.d-10
		zgas(:) = 1.d-10
		vxgas(:) = 1.d-10
		vygas(:) = 1.d-10
		vzgas(:) = 1.d-10
		density_gas(:) = 1.d-10
		temperature_gas(:) = 1.d-10
		SNII_gas(:) = 1.d-10
		SNI_gas(:) = 1.d-10
		i = 1

		DO WHILE (i .le. nstop)
			read (12,end=6) cell_size_gas(i), xgas(i), ygas(i), zgas(i), Vxgas(i), Vygas(i), Vzgas(i), &
			& density_gas(i), temperature_gas(i), SNII_gas(i), SNI_gas(i)
			i = i + 1
		end do
 6		continue
		close (12)
		print *, i, nstop
	end subroutine data_gas
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine data_dm(filename, nstop)
		implicit none
		character (*),intent (in) :: filename
		integer,intent(in) :: nstop
		integer :: i

		open ( 13 , file = filename, form = 'unformatted' )
		iddm(:) = 0
		xdm(:) = 1.d-10
		ydm(:) = 1.d-10
		zdm(:) = 1.d-10
		vxdm(:) = 1.d-10
		vydm(:) = 1.d-10
		vzdm(:) = 1.d-10
		mass_dm(:) = 1.d-10
		i = 1
		DO WHILE (i .le. nstop)
			read (13,end=6) iddm(i), xdm(i), ydm(i), zdm(i), Vxdm(i), Vydm(i), Vzdm(i) ,mass_dm(i)
			i = i + 1
		end do
 6		continue
		close (13)
		print *, i, nstop
	end subroutine data_dm
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine data_stars(filename, nstop)
		implicit none
		character (*),intent (in) :: filename
		integer,intent(in) :: nstop
		integer :: i

		open ( 14 , file = filename, form = 'unformatted' )
		idstars(:) = 0
		xstars(:) = 1.d-10
		ystars(:) = 1.d-10
		zstars(:) = 1.d-10
		vxstars(:) = 1.d-10
		vystars(:) = 1.d-10
		vzstars(:) = 1.d-10
		mass_stars(:) = 1.d-10
		age_stars(:) = 1.d-10
		i = 1
		DO WHILE (i .le. nstop)
			read (14,end=6) idstars(i), xstars(i), ystars(i), zstars(i), Vxstars(i), &
			& Vystars(i), Vzstars(i), mass_stars(i), age_stars(i)
			i = i + 1
		end do
 6		continue
		close (14)
		print *, i, nstop
	end subroutine data_stars
end module read_binary
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module splitter
!!! Splits and refines cells for interpolation using oct-tree method.
!!! Also contains routines for rotating data to arbitrary frame
use parameters
use globvar
	implicit none

	real(4) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, vec(3), vec_vel(3)
	real(4),allocatable :: xprime(:), yprime(:), zprime(:), vxprime(:), vyprime(:), vzprime(:), rprime(:)
	real(4) :: saxis1(3), saxis2(3), saxis3(3), theta, phi
	real(4),parameter :: ocx(8) = (/ 1.0_4,  1.0_4, -1.0_4, -1.0_4,  1.0_4,  1.0_4, -1.0_4, -1.0_4 /)
	real(4),parameter :: ocy(8) = (/ 1.0_4, -1.0_4,  1.0_4, -1.0_4,  1.0_4, -1.0_4,  1.0_4, -1.0_4 /)
	real(4),parameter :: ocz(8) = (/ 1.0_4,  1.0_4,  1.0_4,  1.0_4, -1.0_4, -1.0_4, -1.0_4, -1.0_4 /)
contains
!!!	CREATE ORTHONORMAL CO-ORDINATE AXES
!!!	----------------------------------------
	subroutine axes(a1,a2,a3)
		implicit none
		real(4),intent(in) :: a3(3)
		real(4),intent(inout) :: a1(3),a2(3)

		theta = acos(a3(3))
		if(a3(1).eq.0.0_4.and.a3(2).eq.0.0_4) then
			phi = 0.0_4
		elseif(a3(1)>=0.0_4.and.a3(2)>=0.0_4) then
			phi = atan(a3(2)/a3(1))
		elseif(a3(1)<0.0_4) then
			phi = pi+atan(a3(2)/a3(1))
		else
			phi = pi2+atan(a3(2)/a3(1))
		end if
		a1(:) = (/ cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta) /)
		a2(:) = (/ -sin(phi),cos(phi),0.0_4 /)	
	end subroutine axes

!!!	DEALLOCATE PRIMED VALUES AND RE-ALLOCATE ACCORDINGLY
!!!	----------------------------------------
	subroutine allocation(nsplit)
		implicit none
		integer,intent(in) :: nsplit
		integer :: i

		call deallocate_primes()
		allocate( xprime(nsplit),stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of xprime. stat= ', i
			stop
		end if
		allocate( yprime(nsplit),stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of yprime. stat= ', i
			stop
		end if
		allocate( zprime(nsplit),stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of zprime. stat= ', i
			stop
		end if
		allocate( vxprime(nsplit),stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of vxprime. stat= ', i
			stop
		end if
		allocate( vyprime(nsplit),stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of vyprime. stat= ', i
			stop
		end if
		allocate( vzprime(nsplit),stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of vzprime. stat= ', i
			stop
		end if
		allocate( rprime(nsplit),stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of rprime. stat= ', i
			stop
		end if
	end subroutine allocation

!!!	DEALLOCATE PRIMED VALUES
!!!	----------------------------------------
	subroutine deallocate_primes()
		implicit none

		if(allocated(xprime))  deallocate(xprime)
		if(allocated(yprime))  deallocate(yprime)
		if(allocated(zprime))  deallocate(zprime)
		if(allocated(vxprime)) deallocate(vxprime)
		if(allocated(vyprime)) deallocate(vyprime)
		if(allocated(vzprime)) deallocate(vzprime)
		if(allocated(rprime))  deallocate(rprime)
	end subroutine deallocate_primes

!!!	ROTATE PARTICLE DATA - NO CELL SPLITTING
!!!	----------------------------------------
	subroutine split0(x, y, z, rcm, vx, vy, vz, vcm)
		implicit none
		real(4), intent(in) :: x, y, z, rcm(3)
		real(4), optional, intent(in) :: vx, vy, vz, vcm(3)

		call allocation(1)

		vec(:) = (/ x-rcm(1),y-rcm(2),z-rcm(3) /)				!!! kpc
		xprime(1) = dot_product(vec,saxis1)					!!! kpc
		yprime(1) = dot_product(vec,saxis2)					!!! kpc
		zprime(1) = dot_product(vec,saxis3)					!!! kpc
		rprime(1) = sqrt( xprime(1)**2 + yprime(1)**2 )				!!! kpc

		if(present(vx) .and. present(vy) .and. present(vz) .and. present(vcm) ) then
			vec_vel(:) = (/ vx-vcm(1),vy-vcm(2),vz-vcm(3) /)			!!! km/s
			vxprime(1) = dot_product(vec_vel,saxis1)				!!! km/s
			vyprime(1) = dot_product(vec_vel,saxis2)				!!! km/s
			vzprime(1) = dot_product(vec_vel,saxis3)				!!! km/s
		end if
	end subroutine split0

!!!	ROTATE HIGHEST DENSITY GAS DATA - SPLIT CELLS ONCE
!!!	--------------------------------------------------
	subroutine split1(x, y, z, rcm, cell_size, vx, vy, vz, vcm)

		implicit none
		real(4), intent(in) :: x, y, z, rcm(3), cell_size
		real(4), optional, intent(in) :: vx, vy, vz, vcm(3)
		integer :: k1

		call allocation(8)

		if(present(vx) .and. present(vy) .and. present(vz) .and. present(vcm) ) then
			do k1=1,8
				x1 = x + ocx(k1)*cell_size/4000.0_4				!!! kpc
				y1 = y + ocy(k1)*cell_size/4000.0_4				!!! kpc
				z1 = z + ocz(k1)*cell_size/4000.0_4				!!! kpc

				vec(:) = (/ x1-rcm(1),y1-rcm(2),z1-rcm(3) /)			!!! kpc
				xprime(k1) = dot_product(vec,saxis1)				!!! kpc
				yprime(k1) = dot_product(vec,saxis2)				!!! kpc
				zprime(k1) = dot_product(vec,saxis3)				!!! kpc
				rprime(k1) = sqrt( xprime(k1)**2 + yprime(k1)**2 )		!!! kpc

				vec_vel(:) = (/ vx-vcm(1),vy-vcm(2),vz-vcm(3) /)		!!! km/s
				vxprime(k1) = dot_product(vec_vel,saxis1)			!!! km/s
				vyprime(k1) = dot_product(vec_vel,saxis2)			!!! km/s
				vzprime(k1) = dot_product(vec_vel,saxis3)			!!! km/s
			end do
		else
			do k1=1,8
				x1 = x + ocx(k1)*cell_size/4000.0_4				!!! kpc
				y1 = y + ocy(k1)*cell_size/4000.0_4				!!! kpc
				z1 = z + ocz(k1)*cell_size/4000.0_4				!!! kpc

				vec(:) = (/ x1-rcm(1),y1-rcm(2),z1-rcm(3) /)			!!! kpc
				xprime(k1) = dot_product(vec,saxis1)				!!! kpc
				yprime(k1) = dot_product(vec,saxis2)				!!! kpc
				zprime(k1) = dot_product(vec,saxis3)				!!! kpc
				rprime(k1) = sqrt( xprime(k1)**2 + yprime(k1)**2 )		!!! kpc
			end do
		end if
	end subroutine split1

!!!	ROTATE LEVEL 2 GAS DATA - SPLIT CELLS TWICE
!!!	-------------------------------------------
	subroutine split2(x, y, z, rcm, cell_size, vx, vy, vz, vcm)
		implicit none
		real(4), intent(in) :: x, y, z, rcm(3), cell_size
		real(4), optional, intent(in) :: vx, vy, vz, vcm(3)
		integer :: k1,k2

		call allocation(64)

		if(present(vx) .and. present(vy) .and. present(vz) .and. present(vcm) ) then
			do k1=1,8
				x1 = x + ocx(k1)*cell_size/4000.0_4				!!! kpc
				y1 = y + ocy(k1)*cell_size/4000.0_4				!!! kpc
				z1 = z + ocz(k1)*cell_size/4000.0_4				!!! kpc
				do k2=1,8
					x2 = x1 + ocx(k2)*cell_size/8000.0_4			!!! kpc
					y2 = y1 + ocy(k2)*cell_size/8000.0_4			!!! kpc
					z2 = z1 + ocz(k2)*cell_size/8000.0_4			!!! kpc

					vec(:) = (/ x2-rcm(1),y2-rcm(2),z2-rcm(3) /)		!!! kpc
					xprime( k2+8*(k1-1) ) = dot_product(vec,saxis1)		!!! kpc
					yprime( k2+8*(k1-1) ) = dot_product(vec,saxis2)		!!! kpc
					zprime( k2+8*(k1-1) ) = dot_product(vec,saxis3)		!!! kpc
					rprime( k2+8*(k1-1) ) = sqrt( xprime( k2+8*(k1-1) )**2 + yprime( k2+8*(k1-1) )**2 )		!!! kpc

					vec_vel(:) = (/ vx-vcm(1),vy-vcm(2),vz-vcm(3) /)	!!! km/s
					vxprime( k2+8*(k1-1) ) = dot_product(vec_vel,saxis1)	!!! km/s
					vyprime( k2+8*(k1-1) ) = dot_product(vec_vel,saxis2)	!!! km/s
					vzprime( k2+8*(k1-1) ) = dot_product(vec_vel,saxis3)	!!! km/s
				end do
			end do
		else
			do k1=1,8
				x1 = x + ocx(k1)*cell_size/4000.0_4				!!! kpc
				y1 = y + ocy(k1)*cell_size/4000.0_4				!!! kpc
				z1 = z + ocz(k1)*cell_size/4000.0_4				!!! kpc
				do k2=1,8
					x2 = x1 + ocx(k2)*cell_size/8000.0_4			!!! kpc
					y2 = y1 + ocy(k2)*cell_size/8000.0_4			!!! kpc
					z2 = z1 + ocz(k2)*cell_size/8000.0_4			!!! kpc

					vec(:) = (/ x2-rcm(1),y2-rcm(2),z2-rcm(3) /)		!!! kpc
					xprime( k2+8*(k1-1) ) = dot_product(vec,saxis1)		!!! kpc
					yprime( k2+8*(k1-1) ) = dot_product(vec,saxis2)		!!! kpc
					zprime( k2+8*(k1-1) ) = dot_product(vec,saxis3)		!!! kpc
					rprime( k2+8*(k1-1) ) = sqrt( xprime( k2+8*(k1-1) )**2 + yprime( k2+8*(k1-1) )**2 )		!!! kpc
				end do
			end do
		end if
	end subroutine split2

!!!	ROTATE LEVEL 3 GAS DATA - SPLIT CELLS 3 TIMES
!!!	---------------------------------------------
	subroutine split3(x, y, z, rcm, cell_size, vx, vy, vz, vcm)
		implicit none
		real(4), intent(in) :: x, y, z, rcm(3), cell_size
		real(4), optional, intent(in) :: vx, vy, vz, vcm(3)
		integer :: k1,k2,k3

		call allocation(512)

		if(present(vx) .and. present(vy) .and. present(vz) .and. present(vcm) ) then
			do k1=1,8
				x1 = x + ocx(k1)*cell_size/4000.0_4									!!! kpc
				y1 = y + ocy(k1)*cell_size/4000.0_4									!!! kpc
				z1 = z + ocz(k1)*cell_size/4000.0_4									!!! kpc
				do k2=1,8
					x2 = x1 + ocx(k2)*cell_size/8000.0_4								!!! kpc
					y2 = y1 + ocy(k2)*cell_size/8000.0_4								!!! kpc
					z2 = z1 + ocz(k2)*cell_size/8000.0_4								!!! kpc
					do k3=1,8
						x3 = x2 + ocx(k3)*cell_size/16000.0_4							!!! kpc
						y3 = y2 + ocy(k3)*cell_size/16000.0_4							!!! kpc
						z3 = z2 + ocz(k3)*cell_size/16000.0_4							!!! kpc

						vec(:) = (/ x3-rcm(1),y3-rcm(2),z3-rcm(3) /)						!!! kpc
						xprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = dot_product(vec,saxis1)				!!! kpc
						yprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = dot_product(vec,saxis2)				!!! kpc
						zprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = dot_product(vec,saxis3)				!!! kpc
						rprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = sqrt( xprime(k3+8*((k2-1)+ 8*(k1-1)))**2 + &
						& yprime(k3+8*((k2-1)+ 8*(k1-1)))**2 )							!!! kpc

						vec_vel(:) = (/ vx-vcm(1),vy-vcm(2),vz-vcm(3) /)					!!! km/s
						vxprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = dot_product(vec_vel,saxis1)			!!! km/s
						vyprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = dot_product(vec_vel,saxis2)			!!! km/s
						vzprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = dot_product(vec_vel,saxis3)			!!! km/s
					end do
				end do
			end do
		else
			do k1=1,8
				x1 = x + ocx(k1)*cell_size/4000.0_4									!!! kpc
				y1 = y + ocy(k1)*cell_size/4000.0_4									!!! kpc
				z1 = z + ocz(k1)*cell_size/4000.0_4									!!! kpc
				do k2=1,8
					x2 = x1 + ocx(k2)*cell_size/8000.0_4								!!! kpc
					y2 = y1 + ocy(k2)*cell_size/8000.0_4								!!! kpc
					z2 = z1 + ocz(k2)*cell_size/8000.0_4								!!! kpc
					do k3=1,8
						x3 = x2 + ocx(k3)*cell_size/16000.0_4							!!! kpc
						y3 = y2 + ocy(k3)*cell_size/16000.0_4							!!! kpc
						z3 = z2 + ocz(k3)*cell_size/16000.0_4							!!! kpc

						vec(:) = (/ x3-rcm(1),y3-rcm(2),z3-rcm(3) /)						!!! kpc
						xprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = dot_product(vec,saxis1)				!!! kpc
						yprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = dot_product(vec,saxis2)				!!! kpc
						zprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = dot_product(vec,saxis3)				!!! kpc
						rprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = sqrt( xprime(k3+8*((k2-1)+ 8*(k1-1)))**2 + &
						& yprime(k3+8*((k2-1)+ 8*(k1-1)))**2 )							!!! kpc
					end do
				end do
			end do
		end if
	end subroutine split3

!!!	ROTATE LEVEL 4 GAS DATA - SPLIT CELLS 4 TIMES
!!!	---------------------------------------------
	subroutine split4(x, y, z, rcm, cell_size, vx, vy, vz, vcm)
		implicit none
		real(4), intent(in) :: x, y, z, rcm(3), cell_size
		real(4), optional, intent(in) :: vx, vy, vz, vcm(3)
		integer :: k1,k2,k3,k4

		call allocation(4096)

		if(present(vx) .and. present(vy) .and. present(vz) .and. present(vcm) ) then
			do k1=1,8
				x1 = x + ocx(k1)*cell_size/4000.0_4										!!! kpc
				y1 = y + ocy(k1)*cell_size/4000.0_4										!!! kpc
				z1 = z + ocz(k1)*cell_size/4000.0_4										!!! kpc
				do k2=1,8
					x2 = x1 + ocx(k2)*cell_size/8000.0_4									!!! kpc
					y2 = y1 + ocy(k2)*cell_size/8000.0_4									!!! kpc
					z2 = z1 + ocz(k2)*cell_size/8000.0_4									!!! kpc
					do k3=1,8
						x3 = x2 + ocx(k3)*cell_size/16000.0_4								!!! kpc
						y3 = y2 + ocy(k3)*cell_size/16000.0_4								!!! kpc
						z3 = z2 + ocz(k3)*cell_size/16000.0_4								!!! kpc
						do k4=1,8
							x4 = x3 + ocx(k4)*cell_size/32000.0_4							!!! kpc
							y4 = y3 + ocy(k4)*cell_size/32000.0_4							!!! kpc
							z4 = z3 + ocz(k4)*cell_size/32000.0_4							!!! kpc

							vec(:) = (/ x4-rcm(1),y4-rcm(2),z4-rcm(3) /)						!!! kpc
							xprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = dot_product(vec,saxis1)		!!! kpc
							yprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = dot_product(vec,saxis2)		!!! kpc
							zprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = dot_product(vec,saxis3)		!!! kpc
							rprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = &
							& sqrt( xprime(k4+8*((k3-1)+8*((k2-1)+8*(k1-1))))**2 + &
							& yprime(k4+8*((k3-1)+8*((k2-1)+8*(k1-1))))**2 )					!!! kpc

							vec_vel(:) = (/ vx-vcm(1),vy-vcm(2),vz-vcm(3) /)					!!! km/s
							vxprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = dot_product(vec_vel,saxis1)	!!! km/s
							vyprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = dot_product(vec_vel,saxis2)	!!! km/s
							vzprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = dot_product(vec_vel,saxis3)	!!! km/s
						end do
					end do
				end do
			end do
		else
			do k1=1,8
				x1 = x + ocx(k1)*cell_size/4000.0_4										!!! kpc
				y1 = y + ocy(k1)*cell_size/4000.0_4										!!! kpc
				z1 = z + ocz(k1)*cell_size/4000.0_4										!!! kpc
				do k2=1,8
					x2 = x1 + ocx(k2)*cell_size/8000.0_4									!!! kpc
					y2 = y1 + ocy(k2)*cell_size/8000.0_4									!!! kpc
					z2 = z1 + ocz(k2)*cell_size/8000.0_4									!!! kpc
					do k3=1,8
						x3 = x2 + ocx(k3)*cell_size/16000.0_4								!!! kpc
						y3 = y2 + ocy(k3)*cell_size/16000.0_4								!!! kpc
						z3 = z2 + ocz(k3)*cell_size/16000.0_4								!!! kpc
						do k4=1,8
							x4 = x3 + ocx(k4)*cell_size/32000.0_4							!!! kpc
							y4 = y3 + ocy(k4)*cell_size/32000.0_4							!!! kpc
							z4 = z3 + ocz(k4)*cell_size/32000.0_4							!!! kpc

							vec(:) = (/ x4-rcm(1),y4-rcm(2),z4-rcm(3) /)						!!! kpc
							xprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = dot_product(vec,saxis1)		!!! kpc
							yprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = dot_product(vec,saxis2)		!!! kpc
							zprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = dot_product(vec,saxis3)		!!! kpc
							rprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = &
							& sqrt( xprime(k4+8*((k3-1)+8*((k2-1)+8*(k1-1))))**2 + &
							& yprime(k4+8*((k3-1)+8*((k2-1)+8*(k1-1))))**2 )					!!! kpc
						end do
					end do
				end do
			end do
		end if
	end subroutine split4
!!!	ROTATE LEVEL 5 GAS DATA - SPLIT CELLS 5 TIMES
!!!	---------------------------------------------
	subroutine split5(x, y, z, rcm, cell_size, vx, vy, vz, vcm)
		implicit none
		real(4), intent(in) :: x, y, z, rcm(3), cell_size
		real(4), optional, intent(in) :: vx, vy, vz, vcm(3)
		integer :: k1,k2,k3,k4,k5

		call allocation(32768)

		if(present(vx) .and. present(vy) .and. present(vz) .and. present(vcm) ) then
			do k1=1,8
				x1 = x + ocx(k1)*cell_size/4000.0_4											!!! kpc
				y1 = y + ocy(k1)*cell_size/4000.0_4											!!! kpc
				z1 = z + ocz(k1)*cell_size/4000.0_4											!!! kpc
				do k2=1,8
					x2 = x1 + ocx(k2)*cell_size/8000.0_4										!!! kpc
					y2 = y1 + ocy(k2)*cell_size/8000.0_4										!!! kpc
					z2 = z1 + ocz(k2)*cell_size/8000.0_4										!!! kpc
					do k3=1,8
						x3 = x2 + ocx(k3)*cell_size/16000.0_4									!!! kpc
						y3 = y2 + ocy(k3)*cell_size/16000.0_4									!!! kpc
						z3 = z2 + ocz(k3)*cell_size/16000.0_4									!!! kpc
						do k4=1,8
							x4 = x3 + ocx(k4)*cell_size/32000.0_4								!!! kpc
							y4 = y3 + ocy(k4)*cell_size/32000.0_4								!!! kpc
							z4 = z3 + ocz(k4)*cell_size/32000.0_4								!!! kpc
							do k5=1,8
								x5 = x4 + ocx(k5)*cell_size/64000.0_4							!!! kpc
								y5 = y4 + ocy(k5)*cell_size/64000.0_4							!!! kpc
								z5 = z4 + ocz(k5)*cell_size/64000.0_4							!!! kpc

								vec(:) = (/ x5-rcm(1),y5-rcm(2),z5-rcm(3) /)						!!! kpc
								xprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec,saxis1)	!!! kpc
								yprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec,saxis2)	!!! kpc
								zprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec,saxis3)	!!! kpc
								rprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = &
								& sqrt( xprime( k5+ 8*((k4-1)+ 8*((k3-1)+8*((k2-1)+8*(k1-1)))) )**2 + &
								& yprime( k5+ 8*((k4-1)+ 8*((k3-1)+8*((k2-1)+8*(k1-1)))) )**2 )				!!! kpc

								vec_vel(:) = (/ vx-vcm(1),vy-vcm(2),vz-vcm(3) /)					    !!! km/s
								vxprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec_vel,saxis1)!!! km/s
								vyprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec_vel,saxis2)!!! km/s
								vzprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec_vel,saxis3)!!! km/s
							end do
						end do
					end do
				end do
			end do
		else
			do k1=1,8
				x1 = x + ocx(k1)*cell_size/4000.0_4											!!! kpc
				y1 = y + ocy(k1)*cell_size/4000.0_4											!!! kpc
				z1 = z + ocz(k1)*cell_size/4000.0_4											!!! kpc
				do k2=1,8
					x2 = x1 + ocx(k2)*cell_size/8000.0_4										!!! kpc
					y2 = y1 + ocy(k2)*cell_size/8000.0_4										!!! kpc
					z2 = z1 + ocz(k2)*cell_size/8000.0_4										!!! kpc
					do k3=1,8
						x3 = x2 + ocx(k3)*cell_size/16000.0_4									!!! kpc
						y3 = y2 + ocy(k3)*cell_size/16000.0_4									!!! kpc
						z3 = z2 + ocz(k3)*cell_size/16000.0_4									!!! kpc
						do k4=1,8
							x4 = x3 + ocx(k4)*cell_size/32000.0_4								!!! kpc
							y4 = y3 + ocy(k4)*cell_size/32000.0_4								!!! kpc
							z4 = z3 + ocz(k4)*cell_size/32000.0_4								!!! kpc
							do k5=1,8
								x5 = x4 + ocx(k5)*cell_size/64000.0_4							!!! kpc
								y5 = y4 + ocy(k5)*cell_size/64000.0_4							!!! kpc
								z5 = z4 + ocz(k5)*cell_size/64000.0_4							!!! kpc

								vec(:) = (/ x5-rcm(1),y5-rcm(2),z5-rcm(3) /)						!!! kpc
								xprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec,saxis1)	!!! kpc
								yprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec,saxis2)	!!! kpc
								zprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec,saxis3)	!!! kpc
								rprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = &
								& sqrt( xprime( k5+ 8*((k4-1)+ 8*((k3-1)+8*((k2-1)+8*(k1-1)))) )**2 + &
								& yprime( k5+ 8*((k4-1)+ 8*((k3-1)+8*((k2-1)+8*(k1-1)))) )**2 )				!!! kpc
							end do
						end do
					end do
				end do
			end do
		end if
	end subroutine split5
end module splitter

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module surface_maps
use parameters
use globvar
use splitter
	implicit none
!!! Grid variables
	integer :: ngrid, nR
	real(4) :: del
	real(4),allocatable :: gas_surface_density(:,:,:), cool_gas_surface_density(:,:,:), hot_gas_surface_density(:,:,:), stellar_surface_density(:,:,:)
	real(4),allocatable :: tot_met_mass(:,:,:)
	real(4),allocatable :: young_stellar_surface_density(:,:,:), old_stellar_surface_density(:,:,:), dm_surface_density(:,:,:)
	real(4),allocatable :: cold_surface_density(:,:,:), baryonic_surface_density(:,:,:), temp_grid(:,:)

contains
	subroutine calculate_grid_size(snapshot, grid_ind)
	! calculates sizes of uni-grids
	! There is a small grid of size grid_size^2 with number of cells nR^2
	! There is a larger grid which has a buffer in order to use fft and not worry about edge effects
	! All lengths are in kpc
	! 'snapshot' is integer index of snapshot from main loop
		implicit none
		integer,intent(in) :: snapshot, grid_ind
		real(4) :: xmax, buffer
		integer :: i

		i = grid_ind
		nR = min(ceiling( (2.0_4*grid_size(i))/res ), Ngrid_max)	!!! Grid goes from '-grid_size' to '+grid_size' and has 'nR' cells 
		if(MOD(nR,2).ne.0) then
			nR = nR + 1
		end if
		del = 2.0_4*grid_size(i) / sngl(nR)		!!! size of cell in grid [kpc]

		buffer = 10.0_4 * narrow(i)			!!! Must be kpc
		xmax = grid_size(i) + buffer			!!! Must be kpc

		ngrid = nR + ceiling(2.0_4*buffer/del)		!!! Grid goes from -xmax to +xmax and has 'ngrid' cells
		if(MOD(ngrid,2).ne.0) then
			ngrid = ngrid + 1			!!! Need even number of cells for fft
		end if
		print *, 'Rvir, Rdisc, Hdisc =', Rvir(snapshot), Rdisc(snapshot), Hdisc(snapshot)
		print *, 'grid_size, xmax =', grid_size(i), xmax
		print *, 'ngrid = ',ngrid, 'ngrid^2 = ',sngl(ngrid)**2
		print *, 'nR=', nR, 'nR^2 = ', sngl(nR)**2
		print *, 'res=', 1000.0_4*res, 'del=', 1000.0_4*del

	end subroutine calculate_grid_size
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine allocate_grids()
		implicit none
		integer :: i

		if( include_gas ) then
			print *, 'allocating gas arrays'
			allocate( gas_surface_density(ngrid,ngrid,3), stat=i )				!!! deallocated at end of subroutine !!!
			if(i.ne.0) then
				print *, 'error in allocation of gas_surface_density. stat= ', i
				stop
			end if
			gas_surface_density(:,:,:) = 1.d-10

			allocate( cool_gas_surface_density(ngrid,ngrid,3), stat=i )				!!! deallocated at end of subroutine !!!
			if(i.ne.0) then
				print *, 'error in allocation of cool_gas_surface_density. stat= ', i
				stop
			end if
			cool_gas_surface_density(:,:,:) = 1.d-10

			allocate( hot_gas_surface_density(ngrid,ngrid,3), stat=i )				!!! deallocated at end of subroutine !!!
			if(i.ne.0) then
				print *, 'error in allocation of hot_gas_surface_density. stat= ', i
				stop
			end if
			hot_gas_surface_density(:,:,:) = 1.d-10

			allocate( tot_met_mass(ngrid,ngrid,3), stat=i )				!!! deallocated at end of subroutine !!!
			if(i.ne.0) then
				print *, 'error in allocation of tot_met_mass. stat= ', i
				stop
			end if
			tot_met_mass(:,:,:) = 1.d-10

			if( include_stars ) then
				allocate( cold_surface_density(ngrid,ngrid,3), stat=i )				!!! deallocated at end of 'add_stars' subroutine !!!
				if(i.ne.0) then
					print *, 'error in allocation of cold_surface_density. stat= ', i
					stop
				end if
				cold_surface_density(:,:,:) = 1.d-10

				allocate( baryonic_surface_density(ngrid,ngrid,3), stat=i )			!!! deallocated at end of 'add_stars' subroutine !!!
				if(i.ne.0) then
					print *, 'error in allocation of baryonic_surface_density. stat= ', i
					stop
				end if
				baryonic_surface_density(:,:,:) = 1.d-10
			end if
		end if

		if( include_stars ) then
			print *, 'allocating star arrays'
			allocate( stellar_surface_density(ngrid, ngrid, 3), stat=i )		!!! deallocated at the end of this subroutine !!!
			if(i.ne.0) then
				print *, 'error in allocation stellar_surface_density. stat= ', i
				stop
			end if
			stellar_surface_density(:,:,:) = 1d-10

			allocate( young_stellar_surface_density(ngrid, ngrid, 3), stat=i )	!!! deallocated at the end of this subroutine !!!
			if(i.ne.0) then
				print *, 'error in allocation young_stellar_surface_density. stat= ', i
				stop
			end if
			young_stellar_surface_density(:,:,:) = 1d-10

			allocate( old_stellar_surface_density(ngrid, ngrid, 3), stat=i )	!!! deallocated at the end of this subroutine !!!
			if(i.ne.0) then
				print *, 'error in allocation young_stellar_surface_density. stat= ', i
				stop
			end if
			old_stellar_surface_density(:,:,:) = 1d-10
		end if

		if( include_dm ) then
			allocate( dm_surface_density(ngrid, ngrid, 3), stat=i )			!!! deallocated after collecting dm data !!!
			if(i.ne.0) then
				print *, 'error in allocation dm_surface_density. stat= ', i
				stop
			end if
			dm_surface_density(:,:,:) = 1d-10
		end if

		allocate( temp_grid(nR,nR), stat=i )						!!! deallocated at end of 'add_dm' subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation of temp_grid. stat= ', i
			stop
		end if
		temp_grid(:,:) = 1.d-10
	end subroutine allocate_grids
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine add_gas(snapshot, grid_ind)
	! Makes 2D grids for the gas
	! The z axis of the grid is defined by AM.
	! 'snapshot' is integer index of snapshot from main loop
		implicit none
		integer,intent(in) :: snapshot, grid_ind
		integer :: i, j, k, m, ip(8), jp(8), kp(8), Ngas10
		real(4) :: rcen(3), split, xgrid, ygrid, zgrid, rthresh1
		real(4),allocatable :: rgas(:)

		vxgas(:) = vxgas(:) - real(vcom(1,snapshot),8)
		vygas(:) = vygas(:) - real(vcom(2,snapshot),8)
		vzgas(:) = vzgas(:) - real(vcom(3,snapshot),8)
		if(move_center) then
			xgas(:) = xgas(:) - real(rcom(1,snapshot),8)
			ygas(:) = ygas(:) - real(rcom(2,snapshot),8)
			zgas(:) = zgas(:) - real(rcom(3,snapshot),8)
		end if
		rcen(:) = 0.0_4
		saxis3(:) = Ldisc(:,snapshot)
		call axes(saxis1,saxis2,saxis3)

		allocate( rgas(Ngas(snapshot)), stat=i )					!!! deallocated at the end of the subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation of rgas. stat= ', i
			stop
		end if
		rgas(:) = sqrt(xgas(:)**2 + ygas(:)**2 + zgas(:)**2) / del
		rthresh1 = sqrt(3.0_4)*sngl(ngrid/2)

		Ngas10 = (Ngas(snapshot) - mod(Ngas(snapshot),10)) / 10				!!! This just helps keep track of where I am in the loop
		do i=1,Ngas(snapshot)
			if ( mod(i,Ngas10) .eq. 0 ) then
				print*, 'i of Ngas',i,'of',Ngas(snapshot)
			end if
			if ( rgas(i) .le. rthresh1 ) then
				if(cell_size_gas(i)<1500.0_4*del) then
					call split0(xgas(i),ygas(i),zgas(i),rcen(:))
				else if(cell_size_gas(i)<3000.0_4*del) then
					call split1(xgas(i),ygas(i),zgas(i),rcen(:),cell_size_gas(i))
				else if(cell_size_gas(i)<5000.0_4*del) then
					call split2(xgas(i),ygas(i),zgas(i),rcen(:),cell_size_gas(i))
				else if(cell_size_gas(i)<9000.0_4*del) then
					call split3(xgas(i),ygas(i),zgas(i),rcen(:),cell_size_gas(i))
				else
					call split4(xgas(i),ygas(i),zgas(i),rcen(:),cell_size_gas(i))
				end if
				k = size(xprime(:))
				split = log(sngl(k)) / log(8.0_4)
				do j=1,k
					xgrid = xprime(j)/del + ngrid/2.0_4	!!! 'xprime'=0 --> 'xgrid'=ngrid/2, 'xprime'=-xmax --> xgrid=0, 'xprime'=xmax --> xgrid=ngrid
					ygrid = yprime(j)/del + ngrid/2.0_4 
					zgrid = zprime(j)/del + ngrid/2.0_4

					ip(1) = floor( xgrid )
					jp(1) = floor( ygrid )
					kp(1) = floor( zgrid )

					ip(2) = ip(1)
					jp(2) = jp(1)
					kp(2) = kp(1)+1

					ip(3) = ip(1)
					jp(3) = jp(1)+1
					kp(3) = kp(1)

					ip(4) = ip(1)+1
					jp(4) = jp(1)
					kp(4) = kp(1)

					ip(5) = ip(1)
					jp(5) = jp(1)+1
					kp(5) = kp(1)+1

					ip(6) = ip(1)+1
					jp(6) = jp(1)
					kp(6) = kp(1)+1

					ip(7) = ip(1)+1
					jp(7) = jp(1)+1
					kp(7) = kp(1)

					ip(8) = ip(1)+1
					jp(8) = jp(1)+1
					kp(8) = kp(1)+1

					do m=1,8
						if(ip(m) .ge. 1 .and. ip(m) .le. ngrid .and. jp(m) .ge. 1 .and. jp(m) .le. ngrid .and. &
						& kp(m) .ge. 1 .and. kp(m) .le. ngrid) then
							if(kp(m) .ge. ngrid/2-nR/2+1 .and. kp(m) .le. ngrid/2+nR/2) then
								if(m .eq. 1) then
									gas_surface_density(ip(m),jp(m),1) = gas_surface_density(ip(m),jp(m),1) + &
									& 0.03363_4 * density_gas(i) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(7)) * (ygrid-jp(7)) )

									tot_met_mass(ip(m),jp(m),1) = tot_met_mass(ip(m),jp(m),1) + &
									& 0.03363_4 * density_gas(i) * ( SNI_gas(i) + SNII_gas(i) ) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(7)) * (ygrid-jp(7)) )

									if( temperature_gas(i) .le. cool_gas ) then
										cool_gas_surface_density(ip(m),jp(m),1) = cool_gas_surface_density(ip(m),jp(m),1) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(7)) * (ygrid-jp(7)) )
									elseif( temperature_gas(i) .ge. hot_gas ) then
										hot_gas_surface_density(ip(m),jp(m),1) = hot_gas_surface_density(ip(m),jp(m),1) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(7)) * (ygrid-jp(7)) )
									end if
								elseif(m .eq. 3) then
						      			gas_surface_density(ip(m),jp(m),1) = gas_surface_density(ip(m),jp(m),1) + &
									& 0.03363_4 * density_gas(i) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(4)) * (ygrid-jp(4)) )

						      			tot_met_mass(ip(m),jp(m),1) = tot_met_mass(ip(m),jp(m),1) + &
									& 0.03363_4 * density_gas(i) * ( SNI_gas(i) + SNII_gas(i) ) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(4)) * (ygrid-jp(4)) )

									if( temperature_gas(i) .le. cool_gas ) then
							      			cool_gas_surface_density(ip(m),jp(m),1) = cool_gas_surface_density(ip(m),jp(m),1) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(4)) * (ygrid-jp(4)) )
									elseif( temperature_gas(i) .ge. hot_gas ) then
							      			hot_gas_surface_density(ip(m),jp(m),1) = hot_gas_surface_density(ip(m),jp(m),1) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(4)) * (ygrid-jp(4)) )
									end if
								elseif(m .eq. 4) then
						      			gas_surface_density(ip(m),jp(m),1) = gas_surface_density(ip(m),jp(m),1) + &
									& 0.03363_4 * density_gas(i) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(3)) * (ygrid-jp(3)) )

						      			tot_met_mass(ip(m),jp(m),1) = tot_met_mass(ip(m),jp(m),1) + &
									& 0.03363_4 * density_gas(i) * ( SNI_gas(i) + SNII_gas(i) ) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(3)) * (ygrid-jp(3)) )

									if( temperature_gas(i) .le. cool_gas ) then
							      			cool_gas_surface_density(ip(m),jp(m),1) = cool_gas_surface_density(ip(m),jp(m),1) + &
										& 0.03363_4*density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(3)) * (ygrid-jp(3)) )
									elseif( temperature_gas(i) .ge. hot_gas ) then
							      			hot_gas_surface_density(ip(m),jp(m),1) = hot_gas_surface_density(ip(m),jp(m),1) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(3)) * (ygrid-jp(3)) )
									end if
								elseif(m .eq. 7) then
									gas_surface_density(ip(m),jp(m),1) = gas_surface_density(ip(m),jp(m),1) + &
									& 0.03363_4 * density_gas(i) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(1))*(ygrid-jp(1)) )

						      			tot_met_mass(ip(m),jp(m),1) = tot_met_mass(ip(m),jp(m),1) + &
									& 0.03363_4 * density_gas(i) * ( SNI_gas(i) + SNII_gas(i) ) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(1)) * (ygrid-jp(1)) )

									if( temperature_gas(i) .le. cool_gas ) then
										cool_gas_surface_density(ip(m),jp(m),1) = cool_gas_surface_density(ip(m),jp(m),1) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(1)) * (ygrid-jp(1)) )
									elseif( temperature_gas(i) .ge. hot_gas ) then
										hot_gas_surface_density(ip(m),jp(m),1) = hot_gas_surface_density(ip(m),jp(m),1) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(1)) * (ygrid-jp(1)) )
									end if
								end if
							end if
							if(ip(m) .ge. ngrid/2-nR/2+1 .and. ip(m) .le. ngrid/2+nR/2) then
								if(m .eq. 1) then
									gas_surface_density(jp(m),kp(m),2) = gas_surface_density(jp(m),kp(m),2) + & 
									& 0.03363_4 * density_gas(i) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (ygrid-jp(5))*(zgrid-kp(5)) )

						      			tot_met_mass(jp(m),kp(m),2) = tot_met_mass(jp(m),kp(m),2) + &
									& 0.03363_4 * density_gas(i) * ( SNI_gas(i) + SNII_gas(i) ) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (ygrid-jp(5)) * (zgrid-kp(5)) )

									if( temperature_gas(i) .le. cool_gas ) then
										cool_gas_surface_density(jp(m),kp(m),2) = cool_gas_surface_density(jp(m),kp(m),2) + & 
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (ygrid-jp(5)) * (zgrid-kp(5)) )
									elseif( temperature_gas(i) .ge. hot_gas ) then
										hot_gas_surface_density(jp(m),kp(m),2) = hot_gas_surface_density(jp(m),kp(m),2) + & 
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (ygrid-jp(5)) * (zgrid-kp(5)) )
									end if
								elseif(m .eq. 2) then
									gas_surface_density(jp(m),kp(m),2) = gas_surface_density(jp(m),kp(m),2) + &
									& 0.03363_4 * density_gas(i) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (ygrid-jp(3)) * (zgrid-kp(3)) )

						      			tot_met_mass(jp(m),kp(m),2) = tot_met_mass(jp(m),kp(m),2) + &
									& 0.03363_4 * density_gas(i) * ( SNI_gas(i) + SNII_gas(i) ) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (ygrid-jp(3)) * (zgrid-kp(3)) )

									if( temperature_gas(i) .le. cool_gas ) then
										cool_gas_surface_density(jp(m),kp(m),2) = cool_gas_surface_density(jp(m),kp(m),2) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (ygrid-jp(3)) * (zgrid-kp(3)) )
									elseif( temperature_gas(i) .ge. hot_gas ) then
										hot_gas_surface_density(jp(m),kp(m),2) = hot_gas_surface_density(jp(m),kp(m),2) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (ygrid-jp(3)) * (zgrid-kp(3)) )
									end if
								elseif(m .eq. 3) then
									gas_surface_density(jp(m),kp(m),2) = gas_surface_density(jp(m),kp(m),2) + &
									& 0.03363_4 * density_gas(i) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (ygrid-jp(2)) * (zgrid-kp(2)) )

						      			tot_met_mass(jp(m),kp(m),2) = tot_met_mass(jp(m),kp(m),2) + &
									& 0.03363_4 * density_gas(i) * ( SNI_gas(i) + SNII_gas(i) ) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (ygrid-jp(2)) * (zgrid-kp(2)) )

									if( temperature_gas(i) .le. cool_gas ) then
										cool_gas_surface_density(jp(m),kp(m),2) = cool_gas_surface_density(jp(m),kp(m),2) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (ygrid-jp(2)) * (zgrid-kp(2)) )
									elseif( temperature_gas(i) .ge. hot_gas ) then
										hot_gas_surface_density(jp(m),kp(m),2) = hot_gas_surface_density(jp(m),kp(m),2) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (ygrid-jp(2)) * (zgrid-kp(2)) )
									end if
								elseif(m .eq. 5) then
									gas_surface_density(jp(m),kp(m),2) = gas_surface_density(jp(m),kp(m),2) + &
									& 0.03363_4 * density_gas(i) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (ygrid-jp(1)) * (zgrid-kp(1)) )

						      			tot_met_mass(jp(m),kp(m),2) = tot_met_mass(jp(m),kp(m),2) + &
									& 0.03363_4 * density_gas(i) * ( SNI_gas(i) + SNII_gas(i) ) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (ygrid-jp(1)) * (zgrid-kp(1)) )

									if( temperature_gas(i) .le. cool_gas ) then
										cool_gas_surface_density(jp(m),kp(m),2) = cool_gas_surface_density(jp(m),kp(m),2) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (ygrid-jp(1)) * (zgrid-kp(1)) )
									elseif( temperature_gas(i) .ge. hot_gas ) then
										hot_gas_surface_density(jp(m),kp(m),2) = hot_gas_surface_density(jp(m),kp(m),2) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (ygrid-jp(1)) * (zgrid-kp(1)) )
									end if
								end if
							end if
							if(jp(m) .ge. ngrid/2-nR/2+1 .and. jp(m) .le. ngrid/2+nR/2) then
								if(m .eq. 1) then
									gas_surface_density(ip(m),kp(m),3) = gas_surface_density(ip(m),kp(m),3) + &
									& 0.03363_4 * density_gas(i) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(6)) * (zgrid-kp(6)) )

						      			tot_met_mass(ip(m),kp(m),3) = tot_met_mass(ip(m),kp(m),3) + &
									& 0.03363_4 * density_gas(i) * ( SNI_gas(i) + SNII_gas(i) ) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(6)) * (zgrid-kp(6)) )

									if( temperature_gas(i) .le. cool_gas ) then
										cool_gas_surface_density(ip(m),kp(m),3) = cool_gas_surface_density(ip(m),kp(m),3) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(6)) * (zgrid-kp(6)) )
									elseif( temperature_gas(i) .ge. hot_gas ) then
										hot_gas_surface_density(ip(m),kp(m),3) = hot_gas_surface_density(ip(m),kp(m),3) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(6))*(zgrid-kp(6)) )
									end if
								elseif(m .eq. 2) then
									gas_surface_density(ip(m),kp(m),3) = gas_surface_density(ip(m),kp(m),3) + &
									& 0.03363_4 * density_gas(i) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(4)) * (zgrid-kp(4)) )

						      			tot_met_mass(ip(m),kp(m),3) = tot_met_mass(ip(m),kp(m),3) + &
									& 0.03363_4 * density_gas(i) * ( SNI_gas(i) + SNII_gas(i) ) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(4)) * (zgrid-kp(4)) )

									if( temperature_gas(i) .le. cool_gas ) then
										cool_gas_surface_density(ip(m),kp(m),3) = cool_gas_surface_density(ip(m),kp(m),3) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(4)) * (zgrid-kp(4)) )
									elseif( temperature_gas(i) .ge. hot_gas ) then
										hot_gas_surface_density(ip(m),kp(m),3) = hot_gas_surface_density(ip(m),kp(m),3) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(4)) * (zgrid-kp(4)) )
									end if
								elseif(m .eq. 4) then
									gas_surface_density(ip(m),kp(m),3) = gas_surface_density(ip(m),kp(m),3) + &
									& 0.03363_4 * density_gas(i) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(2)) * (zgrid-kp(2)) )

						      			tot_met_mass(ip(m),kp(m),3) = tot_met_mass(ip(m),kp(m),3) + &
									& 0.03363_4 * density_gas(i) * ( SNI_gas(i) + SNII_gas(i) ) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(2)) * (zgrid-kp(2)) )

									if( temperature_gas(i) .le. cool_gas ) then
										cool_gas_surface_density(ip(m),kp(m),3) = cool_gas_surface_density(ip(m),kp(m),3) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(2)) * (zgrid-kp(2)) )
									elseif( temperature_gas(i) .ge. hot_gas ) then
										hot_gas_surface_density(ip(m),kp(m),3) = hot_gas_surface_density(ip(m),kp(m),3) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(2)) * (zgrid-kp(2)) )
									end if
								elseif(m .eq. 6) then
									gas_surface_density(ip(m),kp(m),3) = gas_surface_density(ip(m),kp(m),3) + &
									& 0.03363_4 * density_gas(i) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(1)) * (zgrid-kp(1)) )

						      			tot_met_mass(ip(m),kp(m),3) = tot_met_mass(ip(m),kp(m),3) + &
									& 0.03363_4 * density_gas(i) * ( SNI_gas(i) + SNII_gas(i) ) * &
									& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(1)) * (zgrid-kp(1)) )

									if( temperature_gas(i) .le. cool_gas ) then
										cool_gas_surface_density(ip(m),kp(m),3) = cool_gas_surface_density(ip(m),kp(m),3) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(1)) * (zgrid-kp(1)) )
									elseif( temperature_gas(i) .ge. hot_gas ) then
										hot_gas_surface_density(ip(m),kp(m),3) = hot_gas_surface_density(ip(m),kp(m),3) + &
										& 0.03363_4 * density_gas(i) * &
										& ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (xgrid-ip(1))*(zgrid-kp(1)) )
									end if
								end if
							end if
						end if
					end do
				end do
				call deallocate_primes()
			end if
		end do
		print *, 'done with gas grid'
		deallocate( rgas )
		tot_met_mass(:,:,:) = tot_met_mass(:,:,:) / gas_surface_density					!!! unitless
		gas_surface_density(:,:,:) = gas_surface_density(:,:,:) / ((1000.0_4*del)**2)			!!! M_sun pc^{-2}
		cool_gas_surface_density(:,:,:) = cool_gas_surface_density(:,:,:) / ((1000.0_4*del)**2)		!!! M_sun pc^{-2}
		hot_gas_surface_density(:,:,:) = hot_gas_surface_density(:,:,:) / ((1000.0_4*del)**2)		!!! M_sun pc^{-2}
		if( include_stars ) then
			cold_surface_density(:,:,:) = cool_gas_surface_density(:,:,:)
			baryonic_surface_density(:,:,:) = gas_surface_density(:,:,:)
		end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		print *, 'XY gas surface density'
		call fft2d(narrow(grid_ind)/del, gas_surface_density(:,:,1))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  gas_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 1)
			end do
		end do
		write(20) nR, grid_size(grid_ind), Rvir(snapshot), Rdisc(snapshot), Hdisc(snapshot), temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'YZ gas surface density'
		call fft2d(narrow(grid_ind)/del, gas_surface_density(:,:,2))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  gas_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 2)
			end do
		end do
		write(20) 2, temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'XZ gas surface density'
		call fft2d(narrow(grid_ind)/del, gas_surface_density(:,:,3))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  gas_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 3)
			end do
		end do
		write(20) 3, temp_grid
		temp_grid(:,:) = 1.d-10
		deallocate( gas_surface_density )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		print *, 'XY cool gas surface density'
		call fft2d(narrow(grid_ind)/del, cool_gas_surface_density(:,:,1))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  cool_gas_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 1)
			end do
		end do
		write(27) nR, grid_size(grid_ind), Rvir(snapshot), Rdisc(snapshot), Hdisc(snapshot), temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'YZ cool gas surface density'
		call fft2d(narrow(grid_ind)/del, cool_gas_surface_density(:,:,2))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  cool_gas_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 2)
			end do
		end do
		write(27) 2, temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'XZ cool gas surface density'
		call fft2d(narrow(grid_ind)/del, cool_gas_surface_density(:,:,3))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  cool_gas_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 3)
			end do
		end do
		write(27) 3, temp_grid
		temp_grid(:,:) = 1.d-10
		deallocate( cool_gas_surface_density )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		print *, 'XY hot gas surface density'
		call fft2d(narrow(grid_ind)/del, hot_gas_surface_density(:,:,1))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  hot_gas_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 1)
			end do
		end do
		write(28) nR, grid_size(grid_ind), Rvir(snapshot), Rdisc(snapshot), Hdisc(snapshot), temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'YZ hot gas surface density'
		call fft2d(narrow(grid_ind)/del, hot_gas_surface_density(:,:,2))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  hot_gas_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 2)
			end do
		end do
		write(28) 2, temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'XZ hot gas surface density'
		call fft2d(narrow(grid_ind)/del, hot_gas_surface_density(:,:,3))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  hot_gas_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 3)
			end do
		end do
		write(28) 3, temp_grid
		temp_grid(:,:) = 1.d-10
		deallocate( hot_gas_surface_density )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		print *, 'XY total gas metallicity mass weighted'
		call fft2d(narrow(grid_ind)/del, tot_met_mass(:,:,1))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  tot_met_mass(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 1)
			end do
		end do
		write(29) nR, grid_size(grid_ind), Rvir(snapshot), Rdisc(snapshot), Hdisc(snapshot), temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'YZ total gas metallicity mass weighted'
		call fft2d(narrow(grid_ind)/del, tot_met_mass(:,:,2))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  tot_met_mass(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 2)
			end do
		end do
		write(29) 2, temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'XZ total gas metallicity mass weighted'
		call fft2d(narrow(grid_ind)/del, tot_met_mass(:,:,3))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  tot_met_mass(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 3)
			end do
		end do
		write(29) 3, temp_grid
		temp_grid(:,:) = 1.d-10
		deallocate( tot_met_mass )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end subroutine add_gas
!________________________________________________________________________________________________________

	subroutine add_stars(snapshot, grid_ind)
	! Makes 2D grids for stars, young stars, old stars, baryons and H_alpha
	! The z axis of the grid is defined by AM.
	! 'snapshot' is integer index of snapshot from main loop
		implicit none
		integer,intent(in) :: snapshot, grid_ind
		integer :: i, j, m, ip(8), jp(8), kp(8), ntrack
		real(4) :: rcen(3), xgrid, ygrid, zgrid, rthresh1
		real(4),allocatable :: rstars(:)

		print *, 'adding stars to grid'
		saxis3(:) = Ldisc(:,snapshot)
		call axes(saxis1,saxis2,saxis3)

		vxstars(:) = vxstars(:) - real(vcom(1,snapshot),8)
		vystars(:) = vystars(:) - real(vcom(2,snapshot),8)
		vzstars(:) = vzstars(:) - real(vcom(3,snapshot),8)
		if(move_center) then
			xstars(:) = xstars(:) - real(rcom(1,snapshot),8)
			ystars(:) = ystars(:) - real(rcom(2,snapshot),8)
			zstars(:) = zstars(:) - real(rcom(3,snapshot),8)
		end if
		rcen(:) = 0.0_4

		print *, 'going through all stars'
		allocate( rstars(Nstars(snapshot)), stat=i )				!!! deallocated at the end of the subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation rstars. stat= ', i
			stop
		end if
		rstars(:) = sqrt(sngl(xstars(:)**2 + ystars(:)**2 + zstars(:)**2)) / del
		rthresh1 = sqrt(3.0_4)*sngl(ngrid/2)

		ntrack = (Nstars(snapshot) - mod(Nstars(snapshot),10)) / 10		!!! This just helps keep track of where I am in the loop
		do i=1,Nstars(snapshot)
			if ( mod(i,ntrack) .eq. 0 ) then
				print*, 'i of Nstars',i,'of',Nstars(snapshot)
			end if

			if ( rstars(i) .le. rthresh1 ) then
				call split0( sngl(xstars(i)), sngl(ystars(i)), sngl(zstars(i)), rcen(:) )

				xgrid = xprime(1)/del + ngrid/2.0_4	!!! 'xprime'=0 --> 'xgrid'=ngrid/2, 'xprime'=-xmax --> xgrid=0, 'xprime'=xmax --> xgrid=ngrid
				ygrid = yprime(1)/del + ngrid/2.0_4 
				zgrid = zprime(1)/del + ngrid/2.0_4

				ip(1) = floor( xgrid )
				jp(1) = floor( ygrid )
				kp(1) = floor( zgrid )

				ip(2) = ip(1)
				jp(2) = jp(1)
				kp(2) = kp(1)+1

				ip(3) = ip(1)
				jp(3) = jp(1)+1
				kp(3) = kp(1)

				ip(4) = ip(1)+1
				jp(4) = jp(1)
				kp(4) = kp(1)

				ip(5) = ip(1)
				jp(5) = jp(1)+1
				kp(5) = kp(1)+1

				ip(6) = ip(1)+1
				jp(6) = jp(1)
				kp(6) = kp(1)+1

				ip(7) = ip(1)+1
				jp(7) = jp(1)+1
				kp(7) = kp(1)

				ip(8) = ip(1)+1
				jp(8) = jp(1)+1
				kp(8) = kp(1)+1

				do m=1,8
					if(ip(m) .ge. 1 .and. ip(m) .le. ngrid .and. jp(m) .ge. 1 .and. jp(m) .le. ngrid .and. &
					& kp(m) .ge. 1 .and. kp(m) .le. ngrid) then
						if(kp(m) .ge. ngrid/2-nR/2+1 .and. kp(m) .le. ngrid/2+nR/2) then
							if(m .eq. 1) then
								stellar_surface_density(ip(m),jp(m),1) = stellar_surface_density(ip(m),jp(m),1) + & 
								& sngl(mass_stars(i))*abs( (xgrid-ip(7))*(ygrid-jp(7)) )
								if(age_stars(i).le.young_stars) then
									young_stellar_surface_density(ip(m),jp(m),1) = &
									& young_stellar_surface_density(ip(m),jp(m),1) + &
									& sngl(mass_stars(i))*abs( (xgrid-ip(7))*(ygrid-jp(7)) )
								else if(age_stars(i).ge.old_stars) then
									old_stellar_surface_density(ip(m),jp(m),1) = &
									& old_stellar_surface_density(ip(m),jp(m),1) + &
									& sngl(mass_stars(i))*abs( (xgrid-ip(7))*(ygrid-jp(7)) )
								end if
							elseif(m .eq. 3) then
								stellar_surface_density(ip(m),jp(m),1) = stellar_surface_density(ip(m),jp(m),1) + &
								& sngl(mass_stars(i))*abs( (xgrid-ip(4))*(ygrid-jp(4)) )
								if(age_stars(i).le.young_stars) then
									young_stellar_surface_density(ip(m),jp(m),1) = &
									& young_stellar_surface_density(ip(m),jp(m),1) + & 
									& sngl(mass_stars(i))*abs( (xgrid-ip(4))*(ygrid-jp(4)) )
								else if(age_stars(i).ge.old_stars) then
									old_stellar_surface_density(ip(m),jp(m),1) = &
									& old_stellar_surface_density(ip(m),jp(m),1) + &
									& sngl(mass_stars(i))*abs( (xgrid-ip(4))*(ygrid-jp(4)) )
								end if
							elseif(m .eq. 4) then
								stellar_surface_density(ip(m),jp(m),1) = stellar_surface_density(ip(m),jp(m),1) + &
								& sngl(mass_stars(i))*abs( (xgrid-ip(3))*(ygrid-jp(3)) )
								if(age_stars(i).le.young_stars) then
									young_stellar_surface_density(ip(m),jp(m),1) = & 
									& young_stellar_surface_density(ip(m),jp(m),1) + &
									& sngl(mass_stars(i))*abs( (xgrid-ip(3))*(ygrid-jp(3)) )
								else if(age_stars(i).ge.old_stars) then
									old_stellar_surface_density(ip(m),jp(m),1) = & 
									& old_stellar_surface_density(ip(m),jp(m),1) + &
									& sngl(mass_stars(i))*abs( (xgrid-ip(3))*(ygrid-jp(3)) )
								end if
							elseif(m .eq. 7) then
								stellar_surface_density(ip(m),jp(m),1) = stellar_surface_density(ip(m),jp(m),1) + &
								& sngl(mass_stars(i))*abs( (xgrid-ip(1))*(ygrid-jp(1)) )
								if(age_stars(i).le.young_stars) then
									young_stellar_surface_density(ip(m),jp(m),1) = &
									& young_stellar_surface_density(ip(m),jp(m),1) + &
									& sngl(mass_stars(i))*abs( (xgrid-ip(1))*(ygrid-jp(1)) )
								else if(age_stars(i).ge.old_stars) then
									old_stellar_surface_density(ip(m),jp(m),1) = &
									& old_stellar_surface_density(ip(m),jp(m),1) + &
									& sngl(mass_stars(i))*abs( (xgrid-ip(1))*(ygrid-jp(1)) )
								end if
							end if
						end if
						if(ip(m) .ge. ngrid/2-nR/2+1 .and. ip(m) .le. ngrid/2+nR/2) then
							if(m .eq. 1) then
								stellar_surface_density(jp(m),kp(m),2) = stellar_surface_density(jp(m),kp(m),2) + &
								& sngl(mass_stars(i))*abs( (ygrid-jp(5))*(zgrid-kp(5)) )
								if(age_stars(i).le.young_stars) then
									young_stellar_surface_density(jp(m),kp(m),2) = &
									& young_stellar_surface_density(jp(m),kp(m),2) + &
									& sngl(mass_stars(i))*abs( (ygrid-jp(5))*(zgrid-kp(5)) )
								else if(age_stars(i).ge.old_stars) then
									old_stellar_surface_density(jp(m),kp(m),2) = &
									& old_stellar_surface_density(jp(m),kp(m),2) + & 
									& sngl(mass_stars(i))*abs( (ygrid-jp(5))*(zgrid-kp(5)) )
								end if
							elseif(m .eq. 2) then
								stellar_surface_density(jp(m),kp(m),2) = stellar_surface_density(jp(m),kp(m),2) + &
								& sngl(mass_stars(i))*abs( (ygrid-jp(3))*(zgrid-kp(3)) )
								if(age_stars(i).le.young_stars) then
									young_stellar_surface_density(jp(m),kp(m),2) = &
									& young_stellar_surface_density(jp(m),kp(m),2) + &
									& sngl(mass_stars(i))*abs( (ygrid-jp(3))*(zgrid-kp(3)) )
								else if(age_stars(i).ge.old_stars) then
									old_stellar_surface_density(jp(m),kp(m),2) = &
									& old_stellar_surface_density(jp(m),kp(m),2) + &
									& sngl(mass_stars(i))*abs( (ygrid-jp(3))*(zgrid-kp(3)) )
								end if
							elseif(m .eq. 3) then
								stellar_surface_density(jp(m),kp(m),2) = stellar_surface_density(jp(m),kp(m),2) + &
								& sngl(mass_stars(i))*abs( (ygrid-jp(2))*(zgrid-kp(2)) )
								if(age_stars(i).le.young_stars) then
									young_stellar_surface_density(jp(m),kp(m),2) = &
									& young_stellar_surface_density(jp(m),kp(m),2) + &
									& sngl(mass_stars(i))*abs( (ygrid-jp(2))*(zgrid-kp(2)) )
								else if(age_stars(i).ge.old_stars) then
									old_stellar_surface_density(jp(m),kp(m),2) = &
									& old_stellar_surface_density(jp(m),kp(m),2) + &
									& sngl(mass_stars(i))*abs( (ygrid-jp(2))*(zgrid-kp(2)) )
								end if
							elseif(m .eq. 5) then
								stellar_surface_density(jp(m),kp(m),2) = stellar_surface_density(jp(m),kp(m),2) + &
								& sngl(mass_stars(i))*abs( (ygrid-jp(1))*(zgrid-kp(1)) )
								if(age_stars(i).le.young_stars) then
									young_stellar_surface_density(jp(m),kp(m),2) = &
									& young_stellar_surface_density(jp(m),kp(m),2) + &
									& sngl(mass_stars(i))*abs( (ygrid-jp(1))*(zgrid-kp(1)) )
								else if(age_stars(i).ge.old_stars) then
									old_stellar_surface_density(jp(m),kp(m),2) = &
									& old_stellar_surface_density(jp(m),kp(m),2) + &
									& sngl(mass_stars(i))*abs( (ygrid-jp(1))*(zgrid-kp(1)) )
								end if
							end if
						end if
						if(jp(m) .ge. ngrid/2-nR/2+1 .and. jp(m) .le. ngrid/2+nR/2) then
							if(m .eq. 1) then
								stellar_surface_density(ip(m),kp(m),3) = stellar_surface_density(ip(m),kp(m),3) + &
								& sngl(mass_stars(i))*abs( (xgrid-ip(6))*(zgrid-kp(6)) )
								if(age_stars(i).le.young_stars) then
									young_stellar_surface_density(ip(m),kp(m),3) = &
									& young_stellar_surface_density(ip(m),kp(m),3) + &
									& sngl(mass_stars(i))*abs( (xgrid-ip(6))*(zgrid-kp(6)) )
								else if(age_stars(i).ge.old_stars) then
									old_stellar_surface_density(ip(m),kp(m),3) = & 
									& old_stellar_surface_density(ip(m),kp(m),3) + &
									& sngl(mass_stars(i))*abs( (xgrid-ip(6))*(zgrid-kp(6)) )
								end if
							elseif(m .eq. 2) then
								stellar_surface_density(ip(m),kp(m),3) = stellar_surface_density(ip(m),kp(m),3) + & 
								& sngl(mass_stars(i))*abs( (xgrid-ip(4))*(zgrid-kp(4)) )
								if(age_stars(i).le.young_stars) then
									young_stellar_surface_density(ip(m),kp(m),3) = &
									& young_stellar_surface_density(ip(m),kp(m),3) + &
									& sngl(mass_stars(i))*abs( (xgrid-ip(4))*(zgrid-kp(4)) )
								else if(age_stars(i).ge.old_stars) then
									old_stellar_surface_density(ip(m),kp(m),3) = &
									& old_stellar_surface_density(ip(m),kp(m),3) + &
									& sngl(mass_stars(i))*abs( (xgrid-ip(4))*(zgrid-kp(4)) )
								end if
							elseif(m .eq. 4) then
								stellar_surface_density(ip(m),kp(m),3) = stellar_surface_density(ip(m),kp(m),3) + &
								& sngl(mass_stars(i))*abs( (xgrid-ip(2))*(zgrid-kp(2)) )
								if(age_stars(i).le.young_stars) then
									young_stellar_surface_density(ip(m),kp(m),3) = &
									& young_stellar_surface_density(ip(m),kp(m),3) + &
									& sngl(mass_stars(i))*abs( (xgrid-ip(2))*(zgrid-kp(2)) )
								else if(age_stars(i).ge.old_stars) then
									old_stellar_surface_density(ip(m),kp(m),3) = &
									& old_stellar_surface_density(ip(m),kp(m),3) + &
									& sngl(mass_stars(i))*abs( (xgrid-ip(2))*(zgrid-kp(2)) )
								end if
							elseif(m .eq. 6) then
								stellar_surface_density(ip(m),kp(m),3) = stellar_surface_density(ip(m),kp(m),3) + &
								& sngl(mass_stars(i))*abs( (xgrid-ip(1))*(zgrid-kp(1)) )
								if(age_stars(i).le.young_stars) then
									young_stellar_surface_density(ip(m),kp(m),3) = &
									& young_stellar_surface_density(ip(m),kp(m),3) + &
									& sngl(mass_stars(i))*abs( (xgrid-ip(1))*(zgrid-kp(1)) )
								else if(age_stars(i).ge.old_stars) then
									old_stellar_surface_density(ip(m),kp(m),3) = &
									& old_stellar_surface_density(ip(m),kp(m),3) + &
									& sngl(mass_stars(i))*abs( (xgrid-ip(1))*(zgrid-kp(1)) )
								end if
							end if
						end if
					end if
				end do
				call deallocate_primes()
			end if
		end do
		print *, 'done with star grid'
		deallocate( rstars )

		stellar_surface_density(:,:,:) = stellar_surface_density(:,:,:) / ((1000.0_4*del)**2)	   	   	!!! M_sun pc^{-2}
		young_stellar_surface_density(:,:,:) = young_stellar_surface_density(:,:,:) / ((1000.0_4*del)**2)	!!! M_sun pc^{-2}
		old_stellar_surface_density(:,:,:) = old_stellar_surface_density(:,:,:) / ((1000.0_4*del)**2) 	 	!!! M_sun pc^{-2}
		cold_surface_density(:,:,:) = cold_surface_density(:,:,:) + young_stellar_surface_density(:,:,:)
		baryonic_surface_density(:,:,:) = baryonic_surface_density(:,:,:) + stellar_surface_density(:,:,:)

		print *, 'XY star surface density'
		call fft2d(narrow(grid_ind)/del, stellar_surface_density(:,:,1))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  stellar_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 1)
			end do
		end do
		write(21) nR, grid_size(grid_ind), Rvir(snapshot), Rdisc(snapshot), Hdisc(snapshot), temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'YZ star surface density'
		call fft2d(narrow(grid_ind)/del, stellar_surface_density(:,:,2))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  stellar_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 2)
			end do
		end do
		write(21) 2, temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'XZ star surface density'
		call fft2d(narrow(grid_ind)/del, stellar_surface_density(:,:,3))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  stellar_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 3)
			end do
		end do
		write(21) 3, temp_grid
		temp_grid(:,:) = 1.d-10
		deallocate(stellar_surface_density)

		print *, 'XY young star surface density'
		call fft2d(narrow(grid_ind)/del, young_stellar_surface_density(:,:,1))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  young_stellar_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 1)
			end do
		end do
		write(22) nR, grid_size(grid_ind), Rvir(snapshot), Rdisc(snapshot), Hdisc(snapshot), temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'YZ young star surface density'
		call fft2d(narrow(grid_ind)/del, young_stellar_surface_density(:,:,2))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  young_stellar_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 2)
			end do
		end do
		write(22) 2, temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'XZ young star surface density'
		call fft2d(narrow(grid_ind)/del, young_stellar_surface_density(:,:,3))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  young_stellar_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 3)
			end do
		end do
		write(22) 3, temp_grid
		temp_grid(:,:) = 1.d-10
		deallocate(young_stellar_surface_density)

		print *, 'XY old star surface density'
		call fft2d(narrow(grid_ind)/del, old_stellar_surface_density(:,:,1))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  old_stellar_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 1)
			end do
		end do
		write(23) nR, grid_size(grid_ind), Rvir(snapshot), Rdisc(snapshot), Hdisc(snapshot), temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'YZ old star surface density'
		call fft2d(narrow(grid_ind)/del, old_stellar_surface_density(:,:,2))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  old_stellar_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 2)
			end do
		end do
		write(23) 2, temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'XZ old star surface density'
		call fft2d(narrow(grid_ind)/del, old_stellar_surface_density(:,:,3))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  old_stellar_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 3)
			end do
		end do
		write(23) 3, temp_grid
		temp_grid(:,:) = 1.d-10
		deallocate(old_stellar_surface_density)

		print *, 'XY cold surface density'
		call fft2d(narrow(grid_ind)/del, cold_surface_density(:,:,1))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  cold_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 1)
			end do
		end do
		write(24) nR, grid_size(grid_ind), Rvir(snapshot), Rdisc(snapshot), Hdisc(snapshot), temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'YZ cold surface density'
		call fft2d(narrow(grid_ind)/del, cold_surface_density(:,:,2))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  cold_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 2)
			end do
		end do
		write(24) 2, temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'XZ cold surface density'
		call fft2d(narrow(grid_ind)/del, cold_surface_density(:,:,3))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  cold_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 3)
			end do
		end do
		write(24) 3, temp_grid
		temp_grid(:,:) = 1.d-10
		deallocate(cold_surface_density)

		print *, 'XY baryonic surface density'
		call fft2d(narrow(grid_ind)/del, baryonic_surface_density(:,:,1))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  baryonic_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 1)
			end do
		end do
		write(25) nR, grid_size(grid_ind), Rvir(snapshot), Rdisc(snapshot), Hdisc(snapshot), temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'YZ baryonic surface density'
		call fft2d(narrow(grid_ind)/del, baryonic_surface_density(:,:,2))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  baryonic_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 2)
			end do
		end do
		write(25) 2, temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'XZ baryonic surface density'
		call fft2d(narrow(grid_ind)/del, baryonic_surface_density(:,:,3))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  baryonic_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 3)
			end do
		end do
		write(25) 3, temp_grid
		temp_grid(:,:) = 1.d-10

		deallocate( baryonic_surface_density )
	end subroutine add_stars
!________________________________________________________________________________________________________

	subroutine add_dm(snapshot, grid_ind)
	! Makes 2D grids for stars, young stars, old stars, baryons and H_alpha
	! The z axis of the grid is defined by AM.
	! 'snapshot' is integer index of snapshot from main loop
		implicit none
		integer,intent(in) :: snapshot, grid_ind
		integer :: i, j, m, ip(8), jp(8), kp(8), ntrack
		real(4) :: rcen(3), xgrid, ygrid, zgrid, rthresh1
		real(4),allocatable :: rdm(:)

		print *, 'adding dm to grid'
		saxis3(:) = Ldisc(:,snapshot)
		call axes(saxis1,saxis2,saxis3)

		vxdm(:) = vxdm(:) - real(vcom(1,snapshot),8)
		vydm(:) = vydm(:) - real(vcom(2,snapshot),8)
		vzdm(:) = vzdm(:) - real(vcom(3,snapshot),8)
		if(move_center) then
			xdm(:) = xdm(:) - real(rcom(1,snapshot),8)
			ydm(:) = ydm(:) - real(rcom(2,snapshot),8)
			zdm(:) = zdm(:) - real(rcom(3,snapshot),8)
		end if
		rcen(:) = 0.0_4

		allocate( rdm(Ndm(snapshot)), stat=i )					!!! deallocated after collecting dm data !!!
		if(i.ne.0) then
			print *, 'error in allocation rdm. stat= ', i
			stop
		end if
		rdm(:) = sqrt(sngl(xdm(:)**2 + ydm(:)**2 + zdm(:)**2)) / del
		rthresh1 = sqrt(3.0_4)*sngl(ngrid/2)

		ntrack = (Ndm(snapshot) - mod(Ndm(snapshot),10)) / 10			!!! This just helps keep track of where I am in the loop
		do i=1,Ndm(snapshot)
			if ( mod(i,ntrack) .eq. 0 ) then
				print*, 'i of Ndm',i,'of',Ndm(snapshot)
			end if
			if ( rdm(i) .le. rthresh1 ) then
				call split0(sngl(xdm(i)), sngl(ydm(i)), sngl(zdm(i)), rcen)
				xgrid = xprime(1)/del + ngrid/2.0_4	!!! 'xprime'=0 --> 'xgrid'=ngrid/2, 'xprime'=-xmax --> xgrid=0, 'xprime'=xmax --> xgrid=ngrid
				ygrid = yprime(1)/del + ngrid/2.0_4 
				zgrid = zprime(1)/del + ngrid/2.0_4

				ip(1) = floor( xgrid )
				jp(1) = floor( ygrid )
				kp(1) = floor( zgrid )

				ip(2) = ip(1)
				jp(2) = jp(1)
				kp(2) = kp(1)+1

				ip(3) = ip(1)
				jp(3) = jp(1)+1
				kp(3) = kp(1)

				ip(4) = ip(1)+1
				jp(4) = jp(1)
				kp(4) = kp(1)

				ip(5) = ip(1)
				jp(5) = jp(1)+1
				kp(5) = kp(1)+1

				ip(6) = ip(1)+1
				jp(6) = jp(1)
				kp(6) = kp(1)+1

				ip(7) = ip(1)+1
				jp(7) = jp(1)+1
				kp(7) = kp(1)

				ip(8) = ip(1)+1
				jp(8) = jp(1)+1
				kp(8) = kp(1)+1

				do m=1,8
					if(ip(m) .ge. 1 .and. ip(m) .le. ngrid .and. jp(m) .ge. 1 .and. jp(m) .le. ngrid .and. &
					& kp(m) .ge. 1 .and. kp(m) .le. ngrid) then
						if(kp(m) .ge. ngrid/2-nR/2+1 .and. kp(m) .le. ngrid/2+nR/2) then
							if(m .eq. 1) then
								dm_surface_density(ip(m),jp(m),1) = dm_surface_density(ip(m),jp(m),1) + &
								& sngl(mass_dm(i))*abs( (xgrid-ip(7))*(ygrid-jp(7)) )
							elseif(m .eq. 3) then
								dm_surface_density(ip(m),jp(m),1) = dm_surface_density(ip(m),jp(m),1) + &
								& sngl(mass_dm(i))*abs( (xgrid-ip(4))*(ygrid-jp(4)) )
							elseif(m .eq. 4) then
								dm_surface_density(ip(m),jp(m),1) = dm_surface_density(ip(m),jp(m),1) + &
								& sngl(mass_dm(i))*abs( (xgrid-ip(3))*(ygrid-jp(3)) )
							elseif(m .eq. 7) then
								dm_surface_density(ip(m),jp(m),1) = dm_surface_density(ip(m),jp(m),1) + & 
								& sngl(mass_dm(i))*abs( (xgrid-ip(1))*(ygrid-jp(1)) )
							end if
						end if
						if(ip(m) .ge. ngrid/2-nR/2+1 .and. ip(m) .le. ngrid/2+nR/2) then
							if(m .eq. 1) then
								dm_surface_density(jp(m),kp(m),2) = dm_surface_density(jp(m),kp(m),2) + &
								& sngl(mass_dm(i))*abs( (ygrid-jp(5))*(zgrid-kp(5)) )
							elseif(m .eq. 2) then
								dm_surface_density(jp(m),kp(m),2) = dm_surface_density(jp(m),kp(m),2) + &
								& sngl(mass_dm(i))*abs( (ygrid-jp(3))*(zgrid-kp(3)) )
							elseif(m .eq. 3) then
								dm_surface_density(jp(m),kp(m),2) = dm_surface_density(jp(m),kp(m),2) + &
								& sngl(mass_dm(i))*abs( (ygrid-jp(2))*(zgrid-kp(2)) )
							elseif(m .eq. 5) then
								dm_surface_density(jp(m),kp(m),2) = dm_surface_density(jp(m),kp(m),2) + &
								& sngl(mass_dm(i))*abs( (ygrid-jp(1))*(zgrid-kp(1)) )
							end if
						end if
						if(jp(m) .ge. ngrid/2-nR/2+1 .and. jp(m) .le. ngrid/2+nR/2) then
							if(m .eq. 1) then
								dm_surface_density(ip(m),kp(m),3) = dm_surface_density(ip(m),kp(m),3) + &
								& sngl(mass_dm(i))*abs( (xgrid-ip(6))*(zgrid-kp(6)) )
							elseif(m .eq. 2) then
								dm_surface_density(ip(m),kp(m),3) = dm_surface_density(ip(m),kp(m),3) + &
								& sngl(mass_dm(i))*abs( (xgrid-ip(4))*(zgrid-kp(4)) )
							elseif(m .eq. 4) then
								dm_surface_density(ip(m),kp(m),3) = dm_surface_density(ip(m),kp(m),3) + &
								& sngl(mass_dm(i))*abs( (xgrid-ip(2))*(zgrid-kp(2)) )
							elseif(m .eq. 6) then
								dm_surface_density(ip(m),kp(m),3) = dm_surface_density(ip(m),kp(m),3) + &
								& sngl(mass_dm(i))*abs( (xgrid-ip(1))*(zgrid-kp(1)) )
							end if
						end if
					end if
				end do
			end if
		end do
		print *, 'done with dm grid'
		deallocate( rdm )

		dm_surface_density(:,:,:) = dm_surface_density(:,:,:) / ((1000.0_4*del)**2)
		print *, 'XY dm surface density'
		call fft2d(narrow(grid_ind)/del, dm_surface_density(:,:,1))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  dm_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 1)
			end do
		end do
		write(26) nR, grid_size(grid_ind), Rvir(snapshot), Rdisc(snapshot), Hdisc(snapshot), temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'YZ dm surface density'
		call fft2d(narrow(grid_ind)/del, dm_surface_density(:,:,2))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  dm_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 2)
			end do
		end do
		write(26) 2, temp_grid
		temp_grid(:,:) = 1.d-10
		print *, 'XZ dm surface density'
		call fft2d(narrow(grid_ind)/del, dm_surface_density(:,:,3))
		do j=1,nR
			do i=1,nR
				temp_grid(i,j) =  dm_surface_density(ngrid/2-nR/2+i, ngrid/2-nR/2+j, 3)
			end do
		end do
		write(26) 3, temp_grid
		deallocate( dm_surface_density )

      end subroutine add_dm
!________________________________________________________________________________________________________

	subroutine fft2d(Gwidth, grid)
	use MKL_DFTI
! Performs 2-d gaussian interpolation using fft
! Gwidth is FWHM given in number of cells.
! grid is input and output 2D field to be smoothed

		implicit none
		real(4),intent(in) :: Gwidth
		real(4),intent(inout) :: grid(ngrid, ngrid)
		complex(4),allocatable :: fft_in1d(:), gauss_1d(:)
		integer :: i, j, k, error, length(2)
		type(DFTI_DESCRIPTOR), POINTER :: fft_desc

		print *, 'entering fft'
		length = (/ ngrid,ngrid /)

		print *, 'defining fft_in1d matrix'
		allocate( fft_in1d(ngrid**2), stat=i )                      !!! deallocated at the end of this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation fft_in1d. stat= ', i
			stop
		end if
		print *, 'defining gaussian_1d matrix'
		allocate( gauss_1d(ngrid**2), stat=i )                      !!! deallocated at the end of this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation gauss_1d. stat= ', i
			stop
		end if

		do j=1,ngrid
			do i=1,ngrid
				fft_in1d( (j-1)*ngrid + i ) = grid(i,j)
				gauss_1d( (j-1)*ngrid + i ) = 0.5_4**(( ((1.0_4*(ngrid+1))/2.0_4 - i)**2 + &
								&	((1.0_4*(ngrid+1))/2.0_4 - j)**2 ) / ((Gwidth/2.0_4)**2) ) 
			end do
		end do
		gauss_1d(:) = gauss_1d(:) / sum(gauss_1d(:))
		print *, 'defined 1d vectors'

		error=DftiCreateDescriptor(fft_desc, DFTI_SINGLE, DFTI_COMPLEX, 2, length )
		print *, 'error11=',error
		error=DftiCommitDescriptor( fft_desc )
		print *, 'error12=',error
		error = DftiComputeForward(fft_desc, fft_in1d)
		print *, 'error13=',error
		error = DftiFreeDescriptor(fft_desc )
		print *, 'error14=',error

		error=DftiCreateDescriptor(fft_desc, DFTI_SINGLE, DFTI_COMPLEX, 2, length )
		print *, 'error21=',error
		error=DftiCommitDescriptor( fft_desc )
		print *, 'error22=',error
		error = DftiComputeForward(fft_desc, gauss_1d)
		print *, 'error23=',error
		error = DftiFreeDescriptor(fft_desc )
		print *, 'error24=',error

		print *, 'convolving'
		fft_in1d(:) = fft_in1d(:)*gauss_1d(:)
		print *, 'convolved'

		deallocate(gauss_1d)

		error=DftiCreateDescriptor(fft_desc, DFTI_SINGLE, DFTI_COMPLEX, 2, length )
		print *, 'error31=',error
		error=DftiCommitDescriptor( fft_desc )
		print *, 'error32=',error
		error = DftiComputeBackward(fft_desc, fft_in1d)
		print *, 'error33=',error
		error = DftiFreeDescriptor(fft_desc )
		print *, 'error34=',error

		print *, 'making output smoothed matrix'
		do j=1,nR/2
			do i=1,nR/2
				grid(ngrid/2-nR/2+i, ngrid/2-nR/2+j) = fft_in1d( (ngrid-nR/2+j-1)*ngrid + ngrid-nR/2+i )/(ngrid**2)
				grid(ngrid/2+i, ngrid/2+j) = fft_in1d( (j-1)*ngrid + i )/(ngrid**2)
				grid(ngrid/2-nR/2+i, ngrid/2+j) = fft_in1d( (j-1)*ngrid + ngrid-nR/2+i )/(ngrid**2)
				grid(ngrid/2+i, ngrid/2-nR/2+j) = fft_in1d( (ngrid-nR/2+j-1)*ngrid + i )/(ngrid**2)
			end do
		end do

		deallocate(fft_in1d)
		print *, 'done with smoothing'
	end subroutine fft2d
end module surface_maps

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program main
!#include "adef.h"
use parameters
use globvar
use read_binary
use splitter
use surface_maps
!use omp_lib

	implicit none
	character(len=256),allocatable :: gal_name(:)	!!! This is the name of the galaxy, e.g. MW3, VL01, SFG1, MW10, VELA07, VELA_v2_10, etc.
	character(len=4) :: Rvir_string
	character(len=5) :: grid_size_string
	character(len=7) :: comm
	character(len=20) :: snap_tag, snap_tag2
	character(len=256) :: filename, path_name, input_arg
	character(len=256),allocatable :: dm_file_name(:), gas_file_name(:), stars_file_name(:)
	integer :: i, j, k, l, Nsnapshot, Nsnapshot2, Nsimulation
	real(4) :: aexp2, Rv4

!!!!!!!!!! What simulations will we be looking at today? !!!!!!!!!!
	if(iargc().ge.1) then
		call getarg(1,input_arg)
		write(filename,'(a)') trim(input_arg)
	else
		write(filename,'(a)') './surface_densities_input.dat'
	end if
	print *, 'input file'
	print *, trim(filename)

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
		Nsnapshot = 0
		Nsnapshot2 = 0
		!!!!!!!!!! How many snapshots in the simulation? !!!!!!!!!!
		write(filename,'(a,a,a)') '/BIGDATA/nirm/galaxy_inputs/',trim(gal_name(j)),'/galaxy_catalogue/Nir_halo_cat.txt'
		open(unit=16,file=filename,form='formatted')
		read(16,'(1x,i3)') Nsnapshot
		print *, Nsnapshot
		call allocate_global(Nsnapshot)

		!!!!!!!!!! Get array sizes !!!!!!!!!!
		allocate( Ngas(Nsnapshot), Ndm(Nsnapshot), Nstars(Nsnapshot) )          !!! deallocated at the end of simulation loop !!!
		write(filename,'(a,a,a)') '/BIGDATA/nirm/galaxy_inputs/',trim(gal_name(j)),'/Nmax.txt'
		open(unit=17,file=filename,form='formatted')
		read(17,'(1x,i3)') Nsnapshot2
		if( Nsnapshot2 .ne. Nsnapshot ) then
			print *, 'DATA Nsnapshot - error with catalogue files. Inconsistent number of snapshots per galaxy'
			print *, 'Nsnapshot',Nsnapshot,'Nsnapshot2',Nsnapshot2
			stop
		end if
		do i=1,Nsnapshot
			read(17,'(1x,f5.3,3(1x,i))') aexp(i), Ngas(i), Nstars(i), Ndm(i)
		end do
		close(unit=17)
		print *, Ngas(Nsnapshot), Nstars(Nsnapshot), Ndm(Nsnapshot)

		write(filename,'(a,a)') 'mkdir -p ./',trim(gal_name(j))
		call system(trim(filename))

		!!!!!!!!!! Initialize the arrays per simulation !!!!!!!!!!
		do i=1,Nsnapshot
			read(16,'(1x,f5.3,3(1x,es12.5))') aexp2, Rvir(i), Mvir(i), Vvir(i)
			if( aexp2 .ne. aexp(i) ) then
				print *, 'DATA aexp - error with catalogue files. Inconsistent expansion factor'
				print *, 'snapshot #',i,'aexp',aexp(i),'aexp2',aexp2
				stop
			end if
		end do
		close(unit=16)
		redshift(:) = ( 1.0_4 / aexp(:) ) - 1.0_4

		write(filename,'(a,a,a)') '/BIGDATA/nirm/galaxy_inputs/',trim(gal_name(j)),'/galaxy_catalogue/Nir_disc_cat.txt'
		open(unit=18,file=filename,form='formatted')
		read(18,'(1x,i3)') Nsnapshot2
		if( Nsnapshot2 .ne. Nsnapshot ) then
			print *, 'DATA Nsnapshot - error with catalogue files. Inconsistent number of snapshots per galaxy'
			print *, 'Nsnapshot',Nsnapshot,'Nsnapshot2',Nsnapshot2
			stop
		end if
		do i=1,Nsnapshot
			read(18,'(1x,f5.3,10(1x,e12.5),2(1x,f7.3),7(1x,es12.5),2(1x,f7.3))') aexp2, rcom(1,i), rcom(2,i), rcom(3,i), &
			& vcom(1,i), vcom(2,i), vcom(3,i), Ldisc(1,i), Ldisc(2,i), Ldisc(3,i), Lmag(i), Rdisc(i), Hdisc(i), &
			& Mgas_disc(i), Mcold_disc(i), Mstar_disc(i), M_Es_star_disc(i), Mdm_disc(i), SFR_disc(i), age_disc(i), metgas_disc(i), metstars_disc(i)
			if( aexp2 .ne. aexp(i) ) then
				print *, 'DATA aexp - error with catalogue files. Inconsistent expansion factor'
				print *, 'snapshot #',i,'aexp',aexp(i),'aexp2',aexp2
				stop
			end if
		end do
		close(unit=18)

		!!!!!!!!!! Loop over all snapshots !!!!!!!!!!
		write(path_name,'(a,a,a,a)') '/BIGDATA/nirm/SIMULATIONS/',trim(gal_name(j)),'/',trim(gal_name(j))

		do k=1,Nsnapshot
			print *, gal_name(j), k, Nsnapshot
			Rv4 = 4.0_4 * Rvir(k)
			If (Rv4 .ge. 1000) then
				write(Rvir_string,'(i4.4)') int(Rv4)
			else
				write(Rvir_string,'(i3.3,a1)') int(Rv4),'.'
			end if
			write(snap_tag,'(a4,a1,f5.3,a4)') trim(Rvir_string),'a',aexp(k),'.dat'
			write(dm_file_name(k),'(a,a,a)') trim(path_name),'_D',trim(snap_tag)
			write(gas_file_name(k),'(a,a,a)') trim(path_name),'_GZ',trim(snap_tag)
			write(stars_file_name(k),'(a,a,a)') trim(path_name),'_S',trim(snap_tag)
			write(snap_tag2,'(a1,f5.3)') 'a',aexp(k)
			call allocate_input( k )

			print *, trim(dm_file_name(k))
			print *, trim(gas_file_name(k))
			print *, trim(stars_file_name(k))
			print *, Ngas(k), Nstars(k), Ndm(k)

			l = size(grid_size)
			if( l .ne. size(narrow) ) then
				print *, 'grid_size - narrow mis-match'
				stop
			end if
			print *, 'starting to make maps'
			do i=1,l
				if(grid_size(i).lt.10.0_4) then
					write(grid_size_string,'(f3.1)') grid_size(i)
				elseif(grid_size(i).lt.100.0_4) then
					write(grid_size_string,'(f4.1)') grid_size(i)
				else
					write(grid_size_string,'(f5.1)') grid_size(i)
				end if
				write(filename,'(a,a,a,a,a)') 'mkdir -p ./',trim(gal_name(j)),'/',trim(grid_size_string),'_kpc'
				call system(trim(filename))
				call open_files()
				print *, ''
				print *, 'Box size = +/-', grid_size(i),'kpc'
				call calculate_grid_size( k, i )
				call allocate_grids()
				if( include_gas ) then
					call add_gas( k, i )
				end if
				if( include_stars ) then
					call add_stars( k, i )
				end if
				if( include_dm ) then
					call add_dm( k, i )
				end if
				deallocate( temp_grid )
				call deallocate_grids()
				call close_files()
			end do
			call deallocate_data()
			print *, ''
			print *, ''
		end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		deallocate( Ngas, Nstars, Ndm )
		deallocate( dm_file_name, gas_file_name, stars_file_name )
		deallocate( rcom, vcom, Ldisc, Lmag, Rdisc, Hdisc, Mgas_disc, Mcold_disc, Mstar_disc, M_Es_star_disc, Mdm_disc, SFR_disc, age_disc, metgas_disc, metstars_disc )
		deallocate( aexp, Rvir, Mvir, Vvir, redshift )
		print *, ''
	end do
	deallocate( gal_name )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contains
	subroutine allocate_global(N)
!!!!!!!!!! Allocate global arrays per simulation !!!!!!!!!!
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
		allocate( rcom(3,N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation rcom. stat= ', i
			stop
		end if
		allocate( vcom(3,N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation vcom. stat= ', i
			stop
		end if
		allocate( Ldisc(3,N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Ldisc. stat= ', i
			stop
		end if
		allocate( Lmag(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Lmag. stat= ', i
			stop
		end if
		allocate( Rdisc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Rdisc. stat= ', i
			stop
		end if
		allocate( Hdisc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Hdisc. stat= ', i
			stop
		end if
		allocate( Mgas_disc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Mgas_disc. stat= ', i
			stop
		end if
		allocate( Mcold_disc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Mcold_disc. stat= ', i
			stop
		end if
		allocate( Mstar_disc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Mstar_disc. stat= ', i
			stop
		end if
		allocate( M_Es_star_disc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation M_Es_star_disc. stat= ', i
			stop
		end if
		allocate( Mdm_disc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Mdm_disc. stat= ', i
			stop
		end if
		allocate( SFR_disc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation SFR_disc. stat= ', i
			stop
		end if
		allocate( age_disc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation age_disc. stat= ', i
			stop
		end if
		allocate( metgas_disc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation metgas_disc. stat= ', i
			stop
		end if
		allocate( metstars_disc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation metstars_disc. stat= ', i
			stop
		end if
		allocate( dm_file_name(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation dm_file_name. stat= ', i
			stop
		end if
		allocate( gas_file_name(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation gas_file_name. stat= ', i
			stop
		end if
		allocate( stars_file_name(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation stars_file_name. stat= ', i
			stop
		end if

	end subroutine allocate_global
!________________________________________________________________________________________________________
	subroutine allocate_input(snap)
            implicit none
            integer,intent(in) :: snap
            integer :: i

		allocate( xgas(Ngas(snap)), stat=i )			!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation xgas. stat= ', i
			stop
		end if
		allocate( ygas(Ngas(snap)), stat=i )			!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation ygas. stat= ', i
			stop
		end if
		allocate( zgas(Ngas(snap)), stat=i )			!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation zgas. stat= ', i
			stop
		end if
		allocate( density_gas(Ngas(snap)), stat=i )		!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation density_gas. stat= ', i
			stop
		end if
		allocate( cell_size_gas(Ngas(snap)), stat=i )		!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation cell_size_gas. stat= ', i
			stop
		end if
		allocate( vxgas(Ngas(snap)), stat=i )			!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation vxgas. stat= ', i
			stop
		end if
		allocate( vygas(Ngas(snap)), stat=i )			!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation vygas. stat= ', i
			stop
		end if
		allocate( vzgas(Ngas(snap)), stat=i )			!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation vzgas. stat= ', i
			stop
		end if
		allocate( temperature_gas(Ngas(snap)), stat=i )		!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation temperature_gas. stat= ', i
			stop
		end if
		allocate( SNI_gas(Ngas(snap)), stat=i )		!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation SNI_gas. stat= ', i
			stop
		end if
		allocate( SNII_gas(Ngas(snap)), stat=i )		!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation SNII_gas. stat= ', i
			stop
		end if

		call data_gas( gas_file_name(snap), Ngas(snap) )
		res = minval( cell_size_gas(1:Ngas(snap)) ) / 1000.0_4	!!! in kpc

		allocate( xstars(Nstars(snap)), stat=i )		!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation xstars. stat= ', i
			stop
		end if
		allocate( ystars(Nstars(snap)), stat=i )		!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation ystars. stat= ', i
			stop
		end if
		allocate( zstars(Nstars(snap)), stat=i )		!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation zstars. stat= ', i
			stop
		end if
		allocate( vxstars(Nstars(snap)), stat=i )		!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation vxstars. stat= ', i
			stop
		end if
		allocate( vystars(Nstars(snap)), stat=i )		!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation vystars. stat= ', i
			stop
		end if
		allocate( vzstars(Nstars(snap)), stat=i )		!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation vzstars. stat= ', i
			stop
		end if
		allocate( mass_stars(Nstars(snap)), stat=i )		!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation mass_stars. stat= ', i
			stop
		end if
		allocate( idstars(Nstars(snap)), stat=i )		!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation idstars. stat= ', i
			stop
		end if
		allocate( age_stars(Nstars(snap)), stat=i )		!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation age_stars. stat= ', i
			stop
		end if

		call data_stars( stars_file_name(snap), Nstars(snap))
		deallocate( idstars )

		allocate( xdm(Ndm(snap)), stat=i )			!!! deallocated after collecting dm data !!!
		if(i.ne.0) then
			print *, 'error in allocation xdm. stat= ', i
			stop
		end if
		allocate( ydm(Ndm(snap)), stat=i )			!!! deallocated after collecting dm data !!!
		if(i.ne.0) then
			print *, 'error in allocation ydm. stat= ', i
			stop
		end if
		allocate( zdm(Ndm(snap)), stat=i )			!!! deallocated after collecting dm data !!!
		if(i.ne.0) then
			print *, 'error in allocation zdm. stat= ', i
			stop
		end if
		allocate( mass_dm(Ndm(snap)), stat=i )			!!! deallocated after collecting dm data !!!
		if(i.ne.0) then
			print *, 'error in allocation mass_dm. stat= ', i
			stop
		end if
		allocate( vxdm(Ndm(snap)), stat=i )			!!! deallocated immediately !!!
		if(i.ne.0) then
			print *, 'error in allocation vxdm. stat= ', i
			stop
		end if
		allocate( vydm(Ndm(snap)), stat=i )			!!! deallocated immediately !!!
		if(i.ne.0) then
			print *, 'error in allocation vydm. stat= ', i
			stop
		end if
		allocate( vzdm(Ndm(snap)), stat=i )			!!! deallocated immediately !!!
		if(i.ne.0) then
			print *, 'error in allocation vzdm. stat= ', i
			stop
		end if
		allocate( iddm(Ndm(snap)), stat=i )			!!! deallocated immediately !!!
		if(i.ne.0) then
			print *, 'error in allocation iddm. stat= ', i
			stop
		end if

		call data_dm(dm_file_name(snap), Ndm(snap))
		deallocate( iddm )

      end subroutine allocate_input
!________________________________________________________________________________________________________
	subroutine deallocate_data()
!!! This is a failsafe. Actually, everything should be deallocated by now. !!!
		implicit none

		if(allocated(xgas)) deallocate(xgas)
		if(allocated(ygas)) deallocate(ygas)
		if(allocated(zgas)) deallocate(zgas)
		if(allocated(vxgas)) deallocate(vxgas)
		if(allocated(vygas)) deallocate(vygas)
		if(allocated(vzgas)) deallocate(vzgas)
		if(allocated(density_gas)) deallocate(density_gas)
		if(allocated(cell_size_gas)) deallocate(cell_size_gas)
		if(allocated(temperature_gas)) deallocate(temperature_gas)
		if(allocated(SNII_gas)) deallocate(SNII_gas)
		if(allocated(SNI_gas)) deallocate(SNI_gas)

		if(allocated(xstars)) deallocate( xstars )
		if(allocated(ystars)) deallocate( ystars )
		if(allocated(zstars)) deallocate( zstars )
		if(allocated(vxstars)) deallocate( vxstars )
		if(allocated(vystars)) deallocate( vystars )
		if(allocated(vzstars)) deallocate( vzstars )
		if(allocated(mass_stars)) deallocate( mass_stars )
		if(allocated(age_stars)) deallocate( age_stars )
		if(allocated(idstars)) deallocate( idstars )

		if(allocated(xdm)) deallocate( xdm )
		if(allocated(ydm)) deallocate( ydm )
		if(allocated(zdm)) deallocate( zdm )
		if(allocated(vxdm)) deallocate( vxdm )
		if(allocated(vydm)) deallocate( vydm )
		if(allocated(vzdm)) deallocate( vzdm )
		if(allocated(mass_dm)) deallocate( mass_dm )
		if(allocated(iddm)) deallocate( iddm )

		call deallocate_primes()

	end subroutine deallocate_data

!________________________________________________________________________________________________________
	subroutine deallocate_grids()
!!! This is a failsafe. Actually, everything should be deallocated by now. !!!
		implicit none

		if(allocated(gas_surface_density)) deallocate(gas_surface_density)
		if(allocated(cool_gas_surface_density)) deallocate(cool_gas_surface_density)
		if(allocated(hot_gas_surface_density)) deallocate(hot_gas_surface_density)
		if(allocated(tot_met_mass)) deallocate(tot_met_mass)

		if(allocated(stellar_surface_density)) deallocate(stellar_surface_density)
		if(allocated(young_stellar_surface_density)) deallocate(young_stellar_surface_density)
		if(allocated(old_stellar_surface_density)) deallocate(old_stellar_surface_density)

		if(allocated(cold_surface_density)) deallocate(cold_surface_density)
		if(allocated(baryonic_surface_density)) deallocate(baryonic_surface_density)

		if(allocated(dm_surface_density)) deallocate(dm_surface_density)

		if(allocated(temp_grid)) deallocate(temp_grid)

		call deallocate_primes()

	end subroutine deallocate_grids
!________________________________________________________________________________________________________

	subroutine open_files()
		implicit none

		write(filename,'(a,a,a,a,a,a,a)') './',trim(gal_name(j)),'/',trim(grid_size_string),'_kpc/gas_surface_density_',trim(snap_tag2),'.bin'
		open(unit=20,file=filename,form='unformatted')

		write(filename,'(a,a,a,a,a,a,a)') './',trim(gal_name(j)),'/',trim(grid_size_string),'_kpc/stellar_surface_density_',trim(snap_tag2),'.bin'
		open(unit=21,file=filename,form='unformatted')

		write(filename,'(a,a,a,a,a,a,a)') './',trim(gal_name(j)),'/',trim(grid_size_string),'_kpc/young_stellar_surface_density_',trim(snap_tag2),'.bin'
		open(unit=22,file=filename,form='unformatted')

		write(filename,'(a,a,a,a,a,a,a)') './',trim(gal_name(j)),'/',trim(grid_size_string),'_kpc/old_stellar_surface_density_',trim(snap_tag2),'.bin'
		open(unit=23,file=filename,form='unformatted')

		write(filename,'(a,a,a,a,a,a,a)') './',trim(gal_name(j)),'/',trim(grid_size_string),'_kpc/cold_surface_density_',trim(snap_tag2),'.bin'
		open(unit=24,file=filename,form='unformatted')

		write(filename,'(a,a,a,a,a,a,a)') './',trim(gal_name(j)),'/',trim(grid_size_string),'_kpc/baryonic_surface_density_',trim(snap_tag2),'.bin'
		open(unit=25,file=filename,form='unformatted')

		write(filename,'(a,a,a,a,a,a,a)') './',trim(gal_name(j)),'/',trim(grid_size_string),'_kpc/dm_surface_density_',trim(snap_tag2),'.bin'
		open(unit=26,file=filename,form='unformatted')

		write(filename,'(a,a,a,a,a,a,a)') './',trim(gal_name(j)),'/',trim(grid_size_string),'_kpc/cool_gas_surface_density_',trim(snap_tag2),'.bin'
		open(unit=27,file=filename,form='unformatted')

		write(filename,'(a,a,a,a,a,a,a)') './',trim(gal_name(j)),'/',trim(grid_size_string),'_kpc/hot_gas_surface_density_',trim(snap_tag2),'.bin'
		open(unit=28,file=filename,form='unformatted')

		write(filename,'(a,a,a,a,a,a,a)') './',trim(gal_name(j)),'/',trim(grid_size_string),'_kpc/total_gas_metallicity_mass_weighted_',trim(snap_tag2),'.bin'
		open(unit=29,file=filename,form='unformatted')

	end subroutine open_files
!________________________________________________________________________________________________________

	subroutine close_files()
		implicit none

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

	end subroutine close_files
end program main

