module parameters
!!! Parameters used in the code, which can be varied to change performance, memory usage, definitions, implimentations, etc
	implicit none

	real(4),parameter :: cool_gas = 1.5e4							!!! Maximum temperature (in K) for "cool gas"
	real(4),parameter :: hot_gas = 1.0e6							!!! Minimum temperature (in K) for "hot gas"
	real(4),parameter :: pi = 3.141592654_4, pi2 = 2.0_4*pi, pi4_3 = (4.0_4 / 3.0_4)*pi
end module parameters
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module globvar
!!! Global variables and allocatable arrays used throughout the code
	implicit none

!!!	Arrays for gas, stars and DM data
!!!	---------------------------------
	real(4),allocatable :: xgas(:), ygas(:), zgas(:), vxgas(:), vygas(:), vzgas(:)
	real(4),allocatable :: density_gas(:), temperature_gas(:), cell_size_gas(:)
	integer,allocatable :: Ngas(:), Ndm(:), Nstars(:)
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
		i = 1

		DO WHILE (i .le. nstop)
			read (12,end=6) cell_size_gas(i), xgas(i), ygas(i), zgas(i), Vxgas(i), Vygas(i), Vzgas(i), &
			& density_gas(i), temperature_gas(i)
			i = i + 1
		end do
 6		continue
		close (12)
		print *, i, nstop
	end subroutine data_gas
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

program main
use parameters
use globvar
use read_binary
use splitter

	implicit none
	character(len=256),allocatable :: gal_name(:)	!!! This is the name of the galaxy, e.g. MW3, VL01, SFG1, MW10, VELA07, VELA_v2_10, etc.
	character(len=4) :: Rvir_string
	character(len=20) :: snap_tag, snap_tag2
	character(len=256) :: filename, path_name, input_arg
	character(len=256),allocatable :: gas_file_name(:)
	integer :: i, j, k, l, m, n, Nsnapshot, Nsnapshot2, Nsimulation, Nbin, Ngas10
	real(4) :: aexp2, Rv4, dbin, rcen(3), split, ld, rgas
	real(4),allocatable :: bin(:)
	real(8),allocatable :: dens_bin(:,:)

!!!!!!!!!! initialize binning !!!!!!!!!!
	Nbin = 201
	dbin = 0.05_4
	allocate( bin(Nbin), dens_bin(Nbin,8) )
	bin(1:Nbin) = -5.0_4 - dbin + (/ (i,i=1,Nbin) /)*dbin

!!!!!!!!!! What simulations will we be looking at today? !!!!!!!!!!
	if(iargc().ge.1) then
		call getarg(1,input_arg)
		write(filename,'(a)') trim(input_arg)
	else
		write(filename,'(a)') './density_pdf_input.dat'
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
			write(gas_file_name(k),'(a,a,a)') trim(path_name),'_G',trim(snap_tag)
			write(snap_tag2,'(a1,f5.3)') 'a',aexp(k)
			call allocate_input( k )
			print *, trim(gas_file_name(k))
			print *, Ngas(k)

			write(filename,'(a,a,a,a,a)') './',trim(gal_name(j)),'/density_pdf_',trim(snap_tag2),'.out'
			open(unit=20,file=filename)

			rcen(:) = 0.0_4
			saxis3(:) = Ldisc(:,k)
			call axes(saxis1,saxis2,saxis3)

			dens_bin(:,:) = 0.0_8
			Ngas10 = (Ngas(k) - mod(Ngas(k),10)) / 10				!!! This just helps keep track of where I am in the loop
			do i=1,Ngas(k)
				if ( mod(i,Ngas10) .eq. 0 ) then
					print*, 'i of Ngas',i,'of',Ngas(k)
				end if
				rgas = max( abs(xgas(i)), abs(ygas(i)), abs(zgas(i)) )
				if ( rgas .le. Rdisc(k) ) then
					if(cell_size_gas(i)<1.5_4*res) then
						call split0(xgas(i),ygas(i),zgas(i),rcen(:))
					else if(cell_size_gas(i)<3.0_4*res) then
						call split1(xgas(i),ygas(i),zgas(i),rcen(:),cell_size_gas(i))
					else if(cell_size_gas(i)<5.0_4*res) then
						call split2(xgas(i),ygas(i),zgas(i),rcen(:),cell_size_gas(i))
					else if(cell_size_gas(i)<9.0_4*res) then
						call split3(xgas(i),ygas(i),zgas(i),rcen(:),cell_size_gas(i))
					else
						call split4(xgas(i),ygas(i),zgas(i),rcen(:),cell_size_gas(i))
					end if
					l = size(xprime(:))
					split = log(sngl(l)) / log(8.0_4)
					do n=1,l
						if( abs(zprime(n)) .le. Hdisc(k) .and. rprime(n) .le. Rdisc(k) ) then
							do m=1,Nbin
								ld = log10(density_gas(i))
								if( ld .le. bin(m) + 0.5_4*dbin .and. ld .gt. bin(m) - 0.5_4*dbin ) then
									dens_bin(m,1) = dens_bin(m,1) + real( (cell_size_gas(i) / (2.0_4**split))**3, 8)**3					!!! volume weighting
									dens_bin(m,2) = dens_bin(m,2) + real( 0.03363_4 * density_gas(i) * ( cell_size_gas(i)/(2.0_4**split) )**3, 8)		!!! mass weighting
									if( temperature_gas(i) .le. cool_gas ) then
										dens_bin(m,3) = dens_bin(m,3) + real( (cell_size_gas(i) / (2.0_4**split))**3, 8)**3					!!! volume weighting
										dens_bin(m,4) = dens_bin(m,4) + real( 0.03363_4 * density_gas(i) * ( cell_size_gas(i)/(2.0_4**split) )**3, 8)		!!! mass weighting
									elseif( temperature_gas(i) .le. hot_gas ) then
										dens_bin(m,5) = dens_bin(m,5) + real( (cell_size_gas(i) / (2.0_4**split))**3, 8)**3					!!! volume weighting
										dens_bin(m,6) = dens_bin(m,6) + real( 0.03363_4 * density_gas(i) * ( cell_size_gas(i)/(2.0_4**split) )**3, 8)		!!! mass weighting
									else
										dens_bin(m,7) = dens_bin(m,7) + real( (cell_size_gas(i) / (2.0_4**split))**3, 8)**3					!!! volume weighting
										dens_bin(m,8) = dens_bin(m,8) + real( 0.03363_4 * density_gas(i) * ( cell_size_gas(i)/(2.0_4**split) )**3, 8)		!!! mass weighting
									end if
								end if
							end do
						end if
					end do
				end if
			end do
			do m=1,Nbin
				write(20,'(9(1x,es12.5))') bin(m), dens_bin(m,1), dens_bin(m,2), dens_bin(m,3), dens_bin(m,4), dens_bin(m,5), dens_bin(m,6), dens_bin(m,7), dens_bin(m,8)
			end do

			deallocate( xgas, ygas, zgas, cell_size_gas, density_gas, temperature_gas )
			call deallocate_data()
			close(unit=20)
			print *, ''
			print *, ''
		end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		deallocate( Ngas, Nstars, Ndm )
		deallocate( gas_file_name )
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
		allocate( gas_file_name(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation gas_file_name. stat= ', i
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

		call data_gas( gas_file_name(snap), Ngas(snap) )
		deallocate( vxgas, vygas, vzgas )
		res = minval( cell_size_gas(1:Ngas(snap)) )	!!! in pc

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
		call deallocate_primes()

	end subroutine deallocate_data
end program main

