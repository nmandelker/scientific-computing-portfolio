! Try inverse fftin the density and velocity as well to see if the vectors are better conserved

module parameters
!!! Parameters used in the code, which are read in from the file ART_decomposition_input.dat
	implicit none
	! Don't touch these unless you know what you are doing and have a good reason
	real(8),parameter :: G = 4.3E-6		!!! Units: Kpc * (km/sec)^2 * M_sun^(-1)
	real(8),parameter :: pi = 3.141592654

	real(8) :: ngrid_max			!!! Maximum number of possible cells in the uni-grid (CPU memory limited)
	logical :: check_CIC			!!! Verify mass conservation of CIC interpolation during grid creation
	logical :: move_center			!!! Correct galaxy center to center of mass of cold gas (calculated previously)
	logical :: move_velocity_center		!!! Correct galaxy center to center of mass velocity of cold gas (calculated previously)

	real(4) :: res				!!! [pc], resolution of base uni grid
	real(4) :: decomp_weight		!!! vector which is decomposed is rho^{g}*{\vec {v}}. If 0, pure velocity. If 1, momentum. If 0.5, Sqrt[Ek]*{\hat {v}}
	logical :: pre_smooth			!!! If true, smooth data with a Gaussian prior to decomposing
	integer :: pre_smooth_dim		!!! Dimension of Gaussian smoothing. 2D (like in Innoe et al 2015) or 3D
	real(4) :: pre_smooth_FWHM		!!! FWHM of Gaussian used for smoothing the fields, in pc ( sigma ~ 0.425*FWHM )
	integer :: decomp_type			!!! 1 for full field, 2 for rotation curve subtracted, 3 for local average subtracted
	real(4) :: type_2_FWHM			!!! FWHM of Gaussian used to smooth rotation curve before subtracting it
	integer :: Nvert			!!! Number of vertical bins to use when subtracting rotation curve (decomp_type=2)
	integer :: type_2_mass			!!! 1 (0) to subtract mass (volume) weighted rotation curve
	real(4) :: type_3_FWHM			!!! FWHM of Gaussian used to define "background" which is subtracted in decomp_type=3, given in pc

	! These set the scale of the region where everything is calculated, without the buffer
	logical :: use_fixed_R			!!! A fixed size for the disc radius
	real(4) :: Fixed_R			!!! Disc radius
	logical :: use_Rvir			!!! Scale disc radius with the virial radius
	real(4) :: R_over_Rvir			!!! Rd divided by Rvir
	logical :: use_old_R			!!! Use the values of Rd based on cold gas and young stars, from Nir_disc_cat.txt
	logical :: use_fixed_H			!!! A fixed size for the disc height
	real(4) :: Fixed_H			!!! Disc height
	logical :: use_R			!!! Scale disc height with disc radius
	real(4) :: H_over_R			!!! Hd divided by Rd
	logical :: use_old_H			!!! Use the values of Hd based on cold gas and young stars, from Nir_disc_cat.txt

	! This determines if a temperature threshold is used for selecting only "cold" gas
	logical :: use_temp_thresh		!!! limit study to gas below temperature threshold
	real(4) :: temp_thresh

	!!! for Rd = 20kpc, Hd = 4kpc and cells of 100 pc, there are 400 cells per diamter and 80 cells per thickness --> Should FFT a 1024*1024*256 ~ 2.7e8 grid
end module parameters
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module globvar
!!! Global variables and allocatable arrays used throughout the code
	implicit none

!!!	Arrays for gas, stars and DM data
!!!	---------------------------------
	real(4),allocatable :: xgas(:), ygas(:), zgas(:), vxgas(:), vygas(:), vzgas(:)
	real(4),allocatable :: density_gas(:), temperature_gas(:), cell_size_gas(:)
	integer,allocatable :: Ngas(:), Nstars(:), Ndm(:)

!!!	Arrays for disc data to be read from input files
!!!	---------------------------------
	real(4),allocatable :: Rvir(:), aexp(:), rcom(:,:), vcom(:,:), Ldisc(:,:), Lmag(:), Rdisc(:), Hdisc(:)
	real(8),allocatable :: Mvir(:), Vvir(:), Mgas_disc(:), Mcold_disc(:), Mstar_disc(:), M_Es_star_disc(:), Mdm_disc(:)
	real(8),allocatable :: SFR_disc(:), age_disc(:), metgas_disc(:), metstars_disc(:)

	real(4) :: grid_rad, grid_height	!!! physical size of grid
	real(4) :: FWHM_smooth			!!! [pc], FWHM of all smoothing windows - corrected for small discs
	real(4) :: eps4 = 1.e-9
	real(8) :: eps8 = 1.d-9

	character(len=2048) :: output_dirname
	character(len=20) ::snap_tag2
end module globvar
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module read_binary
!!! Reads post-treated ART files of gas, stars and dark matter data
use parameters
use globvar
	implicit none

contains
	subroutine data_gas(filename,nstop)
		implicit none
		character (*),intent (in) :: filename
		integer,intent(in) :: nstop
		integer :: i

		print *, 'reading gas data'
		print *, trim(filename), nstop
		open ( 12 , file = filename, form = 'unformatted' )
		cell_size_gas(:) = 1.e-6
		xgas(:) = 1.e-6
		ygas(:) = 1.e-6
		zgas(:) = 1.e-6
		vxgas(:) = 1.e-6
		vygas(:) = 1.e-6
		vzgas(:) = 1.e-6
		density_gas(:) = 1.e-6
		temperature_gas(:) = 1.e-6
		i = 1

		DO WHILE (i .le. nstop)
			read (12,end=6) cell_size_gas(i), xgas(i), ygas(i), zgas(i), vxgas(i), vygas(i), vzgas(i), &
			& density_gas(i), temperature_gas(i)
			i = i + 1
		end do
 6		continue
		close (12)
		print *, 'data_gas', i, nstop
	end subroutine data_gas
end module read_binary
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
			phi = 2.0_4*pi+atan(a3(2)/a3(1))
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
						rprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = sqrt( xprime(k3+8*((k2-1)+ 8*(k1-1)))**2 + yprime(k3+8*((k2-1)+ 8*(k1-1)))**2 )	!!! kpc

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
						rprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = sqrt( xprime(k3+8*((k2-1)+ 8*(k1-1)))**2 + yprime(k3+8*((k2-1)+ 8*(k1-1)))**2 )	!!! kpc
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
							rprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = sqrt( xprime(k4+8*((k3-1)+8*((k2-1)+8*(k1-1))))**2 + &
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
							rprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = sqrt( xprime(k4+8*((k3-1)+8*((k2-1)+8*(k1-1))))**2 + &
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

								vec(:) = (/ x5-rcm(1),y5-rcm(2),z5-rcm(3) /)							!!! kpc
								xprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec,saxis1)		!!! kpc
								yprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec,saxis2)		!!! kpc
								zprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec,saxis3)		!!! kpc
								rprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = sqrt( xprime( k5+ 8*((k4-1)+ 8*((k3-1)+8*((k2-1)+8*(k1-1)))) )**2 + &
									& yprime( k5+ 8*((k4-1)+ 8*((k3-1)+8*((k2-1)+8*(k1-1)))) )**2 )				!!! kpc

								vec_vel(:) = (/ vx-vcm(1),vy-vcm(2),vz-vcm(3) /)						!!! km/s
								vxprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec_vel,saxis1)	!!! km/s
								vyprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec_vel,saxis2)	!!! km/s
								vzprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec_vel,saxis3)	!!! km/s
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
								rprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = sqrt( xprime( k5+ 8*((k4-1)+ 8*((k3-1)+8*((k2-1)+8*(k1-1)))) )**2 + &
									& yprime( k5+ 8*((k4-1)+ 8*((k3-1)+8*((k2-1)+8*(k1-1)))) )**2 )			!!! kpc
							end do
						end do
					end do
				end do
			end do
		end if
	end subroutine split5
end module splitter
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module decomp
use parameters
use globvar
use read_binary
use splitter
	implicit none
!!! Grid variables
	integer :: ngrid, ngrid2, nR, nH, nbuff
	real(8),allocatable :: density_grid(:,:,:), vx_grid(:,:,:), vy_grid(:,:,:), vz_grid(:,:,:)
	real(8),allocatable :: sigcx(:,:,:), sigcy(:,:,:), sigsx(:,:,:)
	real(4) :: xmax, zmax, disc_rad, disc_height
	real(8) :: sigma_ratio(14)

contains
	subroutine grid_size(Rd, Hd, buffer, galex, snapshot)
	! calculates sizes of uni-grids and makes sure they will not be too large
	! There is a small grid of size 2Rd*2Rd*2Hd with number of cells nR*nR*nH
	! There is a larger grid which has a buffer in order to use fft and not worry about edge effects
	! Rd, Hd and buffer should be given in pc (res, which is a global parameter, is also in pc)
	! 'galex' = a0.xxx
	! 'snapshot' is integer index of snapshot from main loop
		implicit none
		real(4),intent(in) :: Rd, Hd, buffer
		real(8) :: grid_correc, rad_correc
		integer,intent(in) :: snapshot
		character(len=20),intent(in) :: galex

		xmax = Rd + buffer
		zmax = Hd + buffer
		grid_correc = 1.0_8							!!! makes grid smaller if it contains more than 1.5e9 cells
		rad_correc = 1.0_8							!!! makes grid smaller if it contains more than 1.5e9 cells

		ngrid = ceiling((2.0_4*xmax)/res)				!!! Grid goes from -xmax to +xmax and has cell size res
		if(MOD(ngrid,2).ne.0) then
			ngrid = ngrid - 1					!!! Need even number of cells for fft
		end if
		ngrid2 = ceiling((2.0_4*zmax)/res)
		if(MOD(ngrid2,2).ne.0) then
			ngrid2 = ngrid2 - 1
		end if
		print '(a23,1x,f4.1,2(1x,f5.2))', 'AMR res, Rdisc, Hdisc =',minval(cell_size_gas(1:Ngas(snapshot))), Rdisc(snapshot), Hdisc(snapshot)
		print '(a35,1x,f4.1,4(1x,f5.2))', 'grid res, Rgrid, xmax Hgrid, zmax =', res, Rd/1000.0_4, xmax/1000.0_4, Hd/1000.0_4, zmax/1000.0_4
		print '(a19,2(1x,i4),1x,es10.3)', 'ngrid, ngrid2, ntot',ngrid,ngrid2,real(ngrid,8) * real(ngrid,8) * real(ngrid2,8)

		if( real(ngrid,8) * real(ngrid,8) * real(ngrid2,8) > ngrid_max ) then		!!! For memory reasons, don't want more than ngrid_max cells in the grid
			grid_correc = ( real(ngrid,8) * real(ngrid,8) * real(ngrid2,8) / ngrid_max )**( 1.0_8 / 3.0_8 )
			rad_correc = grid_correc*( real(Rd,8) / ( real(Rd,8) + real(buffer,8)*(1.0_8 - grid_correc) ) )	!!! corrects smaller grid according to larger grid correction, making sure to keep buffer
			ngrid = ceiling( real(ngrid,8) / grid_correc )
			if( MOD(ngrid,2) .ne. 0 ) then
				ngrid = ngrid - 1
			end if
			ngrid2 = ceiling( real(ngrid2) / grid_correc )
			if( MOD(ngrid2,2) .ne. 0 ) then
				ngrid2 = ngrid2 - 1
			end if
			print *, 'grid dimensions were too large'
			print *, '(a19,2(1x,i4),1x,es10.3)','ngrid, ngrid2, ntot',ngrid,ngrid2,real(ngrid,8) * real(ngrid,8) * real(ngrid2,8)
			write(19,*) galex, grid_correc, rad_correc
			write(19,*) ''
		end if

		nR = ceiling( 2.0_8*(real(Rd / res,8) / rad_correc) )		!!! Grid goes from '-disc_dim*Rdisc' to '+disc_dim*Rdisc' and has 'nR' cells of size 'res'
		if(MOD(nR,2).ne.0) then
			nR = nR - 1
		end if
		nH = ceiling( 2.0_8*(real(Hd / res,8) / rad_correc) )
		if(MOD(nH,2).ne.0) then
			nH = nH - 1
		end if
		print '(a12,2(1x,i4),1x,es10.3)', 'nR, nH, ntot',nR,nH,real(nR,8)*real(nR,8)*real(nH,8)
	end subroutine grid_size
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine create_grid(snapshot)
	! Creates uniform grids for density and 3 velocity components
	! The z axis of the grid is defined by AM.
	! 'snapshot' is integer index of snapshot from main loop
		implicit none
		integer,intent(in) :: snapshot
		integer :: i, j, k, m, ip(8), jp(8), kp(8), ntrack
		integer :: ntot
		real(4) :: split, rthresh1, rcen(3), vcen(3), mass_temp, ng_half, ng2_half
		real(4),allocatable :: rgas(:)
		real(8) :: test1(4), test2(4)				!!! Used to test mass conservation while making the grid
		real(4) :: max_1 ,max_2, in_grid			!!! Used to test mass conservation while making the grid
		real(4) :: mind, maxvx, maxvy, maxvz

		ng_half  = sngl(ngrid+1)  / 2.0_4
		ng2_half = sngl(ngrid2+1) / 2.0_4
		ntot = ngrid*ngrid*ngrid2
		if(check_CIC) then						!!! Test mass conservation while making the grid
			test1 = 0.0_8
			test2 = 0.0_8
			max_1 = 0.0_8
			max_2 = 0.0_8
		end if
		saxis3(:) = Ldisc(:,snapshot)
		call axes(saxis1,saxis2,saxis3)
		rcen(:) = (/ 0.0_4, 0.0_4, 0.0_4 /)
		vcen(:) = (/ 0.0_4, 0.0_4, 0.0_4 /)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		print *, 'allocating initial grids'
		if( use_temp_thresh ) then
			print *, 'T_gas <',temp_thresh
		else
			print *, 'No temperature threshold'
		end if
		allocate( density_grid(ngrid,ngrid,ngrid2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of density_grid. stat= ', i
			stop
		end if
		allocate( vx_grid(ngrid,ngrid,ngrid2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of vx_grid. stat= ', i
			stop
		end if
		allocate( vy_grid(ngrid,ngrid,ngrid2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of vy_grid. stat= ', i
			stop
		end if
		allocate( vz_grid(ngrid,ngrid,ngrid2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of vz_grid. stat= ', i
			stop
		end if
		density_grid(:,:,:) = 0.0_8
		vx_grid(:,:,:) = 0.0_8
		vy_grid(:,:,:) = 0.0_8
		vz_grid(:,:,:) = 0.0_8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		print *, 'adding gas to grid'
		allocate( rgas(Ngas(snapshot)), stat=i )				!!! deallocated in this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation of rgas. stat= ', i
			stop
		end if
		rgas(:) = 1000.0_4 * sqrt(xgas(:)**2 + ygas(:)**2 + zgas(:)**2) / res	!!! grid units
		rthresh1 = sqrt(2.0_4*(sngl(ng_half)**2) + sngl(ng2_half)**2)
		mind  = minval(density_gas, mask = rgas .lt. rthresh1)
		maxvx = maxval(abs(vxgas), mask = rgas .lt. rthresh1)
		maxvy = maxval(abs(vygas), mask = rgas .lt. rthresh1)
		maxvz = maxval(abs(vzgas), mask = rgas .lt. rthresh1)
		print '(a16,4(1x,es12.5))', 'min dens, max v ', mind, maxvx, maxvy, maxvz
		mind = max( mind * (res**3), 1.e-6 )

		ntrack = (Ngas(snapshot) - mod(Ngas(snapshot),10)) / 10		!!! This just helps keep track of where I am in the loop
		do i=1,Ngas(snapshot)
			if ( i .eq. 1 .or. mod(i,ntrack) .eq. 0 ) then
				print*, 'i of Ngas',i,'of',Ngas(snapshot)
			end if
			if ( rgas(i) .le. rthresh1 ) then
				if( temperature_gas(i) .le. temp_thresh .or. .not. use_temp_thresh ) then
					if(cell_size_gas(i) .le. res) then
						call split0( xgas(i), ygas(i), zgas(i), rcen(:), vxgas(i), vygas(i), vzgas(i), vcen(:) )
					else if(cell_size_gas(i) .le. 2.0_4*res) then
						call split1( xgas(i), ygas(i), zgas(i), rcen(:), cell_size_gas(i), vxgas(i), vygas(i), vzgas(i), vcen(:) )
					else if(cell_size_gas(i) .le. 4.0_4*res) then
						call split2( xgas(i), ygas(i), zgas(i), rcen(:), cell_size_gas(i), vxgas(i), vygas(i), vzgas(i), vcen(:) )
					else if(cell_size_gas(i) .le. 8.0_4*res) then
						call split3( xgas(i), ygas(i), zgas(i), rcen(:), cell_size_gas(i), vxgas(i), vygas(i), vzgas(i), vcen(:) )
!					else if(cell_size_gas(i) .le. 16.0_4*res) then
!						call split4( xgas(i), ygas(i), zgas(i), rcen(:), cell_size_gas(i), vxgas(i), vygas(i), vzgas(i), vcen(:) )
					else
						call split4( xgas(i), ygas(i), zgas(i), rcen(:), cell_size_gas(i), vxgas(i), vygas(i), vzgas(i), vcen(:) )
					end if
					k = size(xprime(:))
					split = log(sngl(k)) / log(8.0_4)
					mass_temp = 0.03363_4*density_gas(i)*( ( cell_size_gas(i) / (2.0_4**split) )**3 )
					xprime = 1000.0_4*xprime/res + ng_half		!!! 'xprime' = {0, -xmax, xmax} --> {(ngrid+1)/2, 0.5, ngrid+0.5}
					yprime = 1000.0_4*yprime/res + ng_half
					zprime = 1000.0_4*zprime/res + ng2_half
					do j=1,k
						ip(1) = floor( xprime(j) )
						jp(1) = floor( yprime(j) )
						kp(1) = floor( zprime(j) )

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

						in_grid = 0
						do m=1,8
							if( ip(m) .ge. 1 .and. ip(m) .le. ngrid .and. jp(m) .ge. 1 .and. jp(m) .le. ngrid .and. kp(m) .ge. 1 .and. kp(m) .le. ngrid2) then
								density_grid(ip(m),jp(m),kp(m)) = density_grid(ip(m),jp(m),kp(m)) + mass_temp * &
									& abs( (xprime(j)-ip(9-m))*(yprime(j)-jp(9-m))*(zprime(j)-kp(9-m)) )	! mass in Solar masses

								vx_grid(ip(m),jp(m),kp(m)) = vx_grid(ip(m),jp(m),kp(m)) + mass_temp * vxprime(j) * &
									& abs( (xprime(j)-ip(9-m))*(yprime(j)-jp(9-m))*(zprime(j)-kp(9-m)) )	! x momentum in Solar masses times km sec^{-1}

								vy_grid(ip(m),jp(m),kp(m)) = vy_grid(ip(m),jp(m),kp(m)) + mass_temp * vyprime(j) * &
									& abs( (xprime(j)-ip(9-m))*(yprime(j)-jp(9-m))*(zprime(j)-kp(9-m)) )	! y momentum in Solar masses times km sec^{-1}

								vz_grid(ip(m),jp(m),kp(m)) = vz_grid(ip(m),jp(m),kp(m)) + mass_temp * vzprime(j) * &
									& abs( (xprime(j)-ip(9-m))*(yprime(j)-jp(9-m))*(zprime(j)-kp(9-m)) )	! z momentum in Solar masses times km sec^{-1}

								in_grid = in_grid + abs( (xprime(j)-ip(9-m))*(yprime(j)-jp(9-m))*(zprime(j)-kp(9-m)) )
							end if
						end do
						if(check_CIC) then
							if( in_grid .gt. 1.0001_4 ) then
								print *, 'SOMETHING IS WRONG WITH YOUR CIC NORMALIZATIO, AKA IN_GRID'
								print *, in_grid
								stop
							end if
							test1(1) = test1(1) + real(in_grid * mass_temp,8)
							test1(2) = test1(2) + real(in_grid * mass_temp * vxprime(j),8)
							test1(3) = test1(3) + real(in_grid * mass_temp * vyprime(j),8)
							test1(4) = test1(4) + real(in_grid * mass_temp * vzprime(j),8)
							if(density_gas(i) .ge. max_1) then
								max_1 = density_gas(i)
							end if
						end if
					end do
					call deallocate_primes()
				end if
			end if
		end do
		if(check_CIC) then
			test2(1) = sum(density_grid(:,:,:))
			test2(2) = sum(vx_grid(:,:,:))
			test2(3) = sum(vy_grid(:,:,:))
			test2(4) = sum(vz_grid(:,:,:))
			max_2 = maxval(density_grid(:,:,:)) / (0.03363_4*res**3)
			print '(a,2(1x,es12.5))', 'density maximum',max_1,max_2
			print '(a,3(1x,es12.5))', 'mass conservation',      test1(1),test2(1),abs( (test2(1)-test1(1))/test1(1) )
			print '(a,3(1x,es12.5))', 'x momentum conservation',test1(2),test2(2),abs( (test2(2)-test1(2))/test1(2) )
			print '(a,3(1x,es12.5))', 'y momentum conservation',test1(3),test2(3),abs( (test2(3)-test1(3))/test1(3) ) 
			print '(a,3(1x,es12.5))', 'z momentum conservation',test1(4),test2(4),abs( (test2(4)-test1(4))/test1(4) )
			if( abs((test2(1)-test1(1))/test1(1)) .gt. 0.01_8 ) then
				print *, 'SERIOUSLY?!? THATS A PRETTY BIG MASS CONSERVATION ERROR DUDE!'
				print *, abs((test2(1)-test1(1))/test1(1))
				stop
			end if
			if( abs((test2(2)-test1(2))/test1(2)) .gt. 0.01_8 ) then
				print *, 'SERIOUSLY?!? THATS A PRETTY BIG Px CONSERVATION ERROR DUDE!'
				print *, abs((test2(2)-test1(2))/test1(2))
				stop
			end if
			if( abs((test2(3)-test1(3))/test1(3)) .gt. 0.01_8 ) then
				print *, 'SERIOUSLY?!? THATS A PRETTY BIG Py CONSERVATION ERROR DUDE!'
				print *, abs((test2(3)-test1(3))/test1(3))
				stop
			end if
			if( abs((test2(4)-test1(4))/test1(4)) .gt. 0.01_8 ) then
				print *, 'SERIOUSLY?!? THATS A PRETTY BIG Pz CONSERVATION ERROR DUDE!'
				print *, abs((test2(4)-test1(4))/test1(4))
				stop
			end if
		end if
		print *, 'done making density and velocity grids'
		deallocate( rgas )

		print *, 'HOW MANY ZEROS 3D, after creating grids: mass, P{x,y,z}?'
		print '(2(1x,i),3(1x,es12.5))', count(mask=density_grid.lt.0.0_8), count(mask=density_grid.eq.0.0_8), minval(density_grid,mask=density_grid.ne.0.0_8), maxval(density_grid), sum(density_grid)/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(vx_grid).lt.0.0_8), count(mask=abs(vx_grid).eq.0.0_8), minval(abs(vx_grid),mask=abs(vx_grid).ne.0.0_8), maxval(abs(vx_grid)), sum(abs(vx_grid))/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(vy_grid).lt.0.0_4), count(mask=abs(vy_grid).eq.0.0_4), minval(abs(vy_grid),mask=abs(vy_grid).ne.0.0_4), maxval(abs(vy_grid)), sum(abs(vy_grid))/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(vz_grid).lt.0.0_4), count(mask=abs(vz_grid).eq.0.0_4), minval(abs(vz_grid),mask=abs(vz_grid).ne.0.0_4), maxval(abs(vz_grid)), sum(abs(vz_grid))/ntot
		print *, ''

		if( pre_smooth ) then
			print *, 'smoothing'
			print *, 'FWHM=', FWHM_smooth, 'cells=', FWHM_smooth*res, 'pc'
			if( pre_smooth_dim .eq. 3 ) then
				print *, 'You elected 3D smoothing'
				call fft_smooth_3D(FWHM_smooth, density_grid, ngrid, ngrid, ngrid2, 1, 2)
				call fft_smooth_3D(FWHM_smooth, vx_grid,      ngrid, ngrid, ngrid2, 1, 2)
				call fft_smooth_3D(FWHM_smooth, vy_grid,      ngrid, ngrid, ngrid2, 1, 2)
				call fft_smooth_3D(FWHM_smooth, vz_grid,      ngrid, ngrid, ngrid2, 1, 2)
			elseif( pre_smooth_dim .eq. 2 ) then
				print *, 'You elected 2D smoothing'
				print *, 'Like in Innoue et al 2015'
				do k=1,ngrid2
					call fft_smooth_2D(FWHM_smooth, density_grid(:,:,k), ngrid, ngrid, 1, 2)
					call fft_smooth_2D(FWHM_smooth, vx_grid(:,:,k),      ngrid, ngrid, 1, 2)
					call fft_smooth_2D(FWHM_smooth, vy_grid(:,:,k),      ngrid, ngrid, 1, 2)
					call fft_smooth_2D(FWHM_smooth, vz_grid(:,:,k),      ngrid, ngrid, 1, 2)
				end do
			end if
			print *, 'HOW MANY ZEROS 3D, after smoothing: mass, P{x,y,z}?'
			print '(2(1x,i),3(1x,es12.5))', count(mask=density_grid.lt.0.0_8), count(mask=density_grid.eq.0.0_8), minval(density_grid,mask=density_grid.ne.0.0_8), maxval(density_grid), sum(density_grid)/ntot
			print '(2(1x,i),3(1x,es12.5))', count(mask=abs(vx_grid).lt.0.0_8), count(mask=abs(vx_grid).eq.0.0_8), minval(abs(vx_grid),mask=abs(vx_grid).ne.0.0_8), maxval(abs(vx_grid)), sum(abs(vx_grid))/ntot
			print '(2(1x,i),3(1x,es12.5))', count(mask=abs(vy_grid).lt.0.0_8), count(mask=abs(vy_grid).eq.0.0_8), minval(abs(vy_grid),mask=abs(vy_grid).ne.0.0_8), maxval(abs(vy_grid)), sum(abs(vy_grid))/ntot
			print '(2(1x,i),3(1x,es12.5))', count(mask=abs(vz_grid).lt.0.0_8), count(mask=abs(vz_grid).eq.0.0_8), minval(abs(vz_grid),mask=abs(vz_grid).ne.0.0_8), maxval(abs(vz_grid)), sum(abs(vz_grid))/ntot
			print *, ''
		end if

		density_grid = max( density_grid, real(mind,8) )	! mass M_sun
		vx_grid = vx_grid / density_grid			! x velocity, km / sec
		vy_grid = vy_grid / density_grid			! y velocity, km / sec
		vz_grid = vz_grid / density_grid			! z velocity, km / sec
		density_grid = density_grid / (real(res,8)**3)		! mass density M_sun pc^{-3}

		print *, 'HOW MANY ZEROS 3D, after renormalizing grids: density, v{x,y,z}?'
		print '(2(1x,i),3(1x,es12.5))', count(mask=density_grid.lt.0.0_8), count(mask=density_grid.eq.0.0_8), minval(density_grid,mask=density_grid.ne.0.0_8), maxval(density_grid), sum(density_grid)/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(vx_grid).lt.0.0_8), count(mask=abs(vx_grid).eq.0.0_8), minval(abs(vx_grid),mask=abs(vx_grid).ne.0.0_8), maxval(abs(vx_grid)), sum(abs(vx_grid))/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(vy_grid).lt.0.0_8), count(mask=abs(vy_grid).eq.0.0_8), minval(abs(vy_grid),mask=abs(vy_grid).ne.0.0_8), maxval(abs(vy_grid)), sum(abs(vy_grid))/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(vz_grid).lt.0.0_8), count(mask=abs(vz_grid).eq.0.0_8), minval(abs(vz_grid),mask=abs(vz_grid).ne.0.0_8), maxval(abs(vz_grid)), sum(abs(vz_grid))/ntot
		print *, ''
		where( abs(vx_grid) .gt. maxvx ) vx_grid = dsign(real(maxvx,8), vx_grid)
		where( abs(vy_grid) .gt. maxvx ) vy_grid = dsign(real(maxvy,8), vy_grid)
		where( abs(vz_grid) .gt. maxvx ) vz_grid = dsign(real(maxvz,8), vz_grid)

	end subroutine create_grid
!_____________________________________________________________________________________________________________________________________________________________________
	subroutine rot_curve_surface_maps(snapshot)
	! Calculates smoothed 1D profiles of face on and edge on surface density and mass-weighted rotation curve
	! May subtract rotation curve from velocity grids
	! Creates projected face on maps of surface density and cylindrical components of velocity (mass weighted)
	! 'snapshot' is integer index of snapshot from main loop
		implicit none
		integer,intent(in) :: snapshot
		integer :: i, j, k, m, n, n1, n2, ng_half, ng2_half, nR_half, nH_half, Nrot, ntot
		integer :: lower_ind, upper_ind, lower_ind2, upper_ind2
		real(8) :: rad, rx, vrot_interp, vx, vy, vx2, vy2, dens, Hw, norm
		real(8),allocatable :: Sigma_gas(:,:), Vz(:,:,:), Vr(:,:,:), Vphi(:,:,:), vrot(:,:), rrot(:), test(:,:)
		integer,allocatable :: Nbin(:)
		character(len=256) :: filename

		print *, ''
		print *, 'making surface plots and getting rotation curve'

		ng_half  = ngrid/2
		ng2_half = ngrid2/2
		nR_half  = nR/2
		nH_half  = nH/2
		Nrot = ceiling( ngrid/sqrt(2.0_4) ) - 1 ! number of cells from center to corner in 2D

		lower_ind  = ng_half  - nR_half + 1
		lower_ind2 = ng2_half - nH_half + 1
		upper_ind  = ng_half  + nR_half
		upper_ind2 = ng2_half + nH_half

		disc_rad = 1000.0_4 * grid_rad
		disc_height = 1000.0_4 * grid_height

		print *, 'Edge on'
		allocate( Sigma_gas(nR,nR), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of Sigma_gas. stat= ', i
			stop
		end if
		allocate( vrot(Nrot,3+3*Nvert), stat=i )	! extra room for azimuthally averaged rotation curve in a square
		if(i.ne.0) then
			print *, 'error in allocation of vrot. stat= ', i
			stop
		end if
		vrot      = 0.0_8
		Sigma_gas(1:nR,1:nH) = sum( density_grid( lower_ind:upper_ind, lower_ind:upper_ind, lower_ind2:upper_ind2 ), dim=2 ) * real(res, 8)	!!! XZ plane M_sun pc^{-2}
		do m=1,nH_half
			vrot(m,1) = vrot(m,1) + sum( Sigma_gas(1:nR,nH_half+m) + Sigma_gas(1:nR,nH_half-m+1), dim=1 )
		end do
		vrot(1:nH_half,1) = vrot(1:nH_half,1) / (2.0_8*real(nR))	!!! Average surface density in the vertical bin, M_{sun} pc^{-2}
		write(filename,'(a,a,a,a)') trim(output_dirname),'/rotation_curves/',trim(snap_tag2),'_edge_on_profile.out'
		open(unit=21,file=filename)
		do i=1,nH_half
			write(21,'(2(1x,es12.5))') ( sngl(i) - 0.5_4 ) * res, vrot(i,1)
		end do
		write(21,'(2(1x,es12.5))') disc_height, 1000.0_4 * Hdisc(snapshot)
		close(unit=21)

		print *, 'Face on'
		allocate( Vz(nR,nR,2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of Vz. stat= ', i
			stop
		end if
		print *, 'Vz'
		allocate( Vr(nR,nR,2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of Vr. stat= ', i
			stop
		end if
		print *, 'Vr'
		allocate( Vphi(nR,nR,2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of Vphi. stat= ', i
			stop
		end if
		print *, 'Vphi'
		allocate( rrot(Nrot), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of rrot. stat= ', i
			stop
		end if
		print *, 'rrot'
		allocate( Nbin(Nrot), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of Nbin. stat= ', i
			stop
		end if
		print *, 'Nbin'
		vrot(:,:) = 0.0_8
		rrot(:) = 0.0_8
		Nbin = 0
		rx = real(ngrid+1)/2.0_8
		print *, 'getting curve'
		do j=1,ngrid
			do i=1,ngrid
				rad  = sqrt( (real(i)-rx)**2 + (real(j)-rx)**2 )
				dens = sum( density_grid(i,j,lower_ind2:upper_ind2) ) * real(res,8)
				vx   = sum( density_grid(i,j,lower_ind2:upper_ind2) * vx_grid(i,j,lower_ind2:upper_ind2) ) * real(res,8)
				vy   = sum( density_grid(i,j,lower_ind2:upper_ind2) * vy_grid(i,j,lower_ind2:upper_ind2) ) * real(res,8)
				vx2  = sum( vx_grid(i,j,lower_ind2:upper_ind2) ) / nH
				vy2  = sum( vy_grid(i,j,lower_ind2:upper_ind2) ) / nH
				if( i .ge. lower_ind .and. i .le. upper_ind .and. j .ge. lower_ind .and. j .le. upper_ind ) then
					k = i-lower_ind+1
					m = j-lower_ind+1
					Sigma_gas(k,m) = max( dens, 1.d-10 )
					Vz(k,m,1)   = sum( density_grid(i,j,lower_ind2:upper_ind2) * vz_grid(i,j,lower_ind2:upper_ind2) ) * real(res,8) / Sigma_gas(k,m)
					Vr(k,m,1)   = ( vx*(real(i)-rx)/rad + vy*(real(j)-rx)/rad ) / Sigma_gas(k,m)
					Vphi(k,m,1) = ( vx*(rx-real(j))/rad + vy*(real(i)-rx)/rad ) / Sigma_gas(k,m)
					Vz(k,m,2)   = sum( vz_grid(i,j,lower_ind2:upper_ind2) ) / nH
					Vr(k,m,2)   = vx2*(real(i)-rx)/rad + vy2*(real(j)-rx)/rad
					Vphi(k,m,2) = vx2*(rx-real(j))/rad + vy2*(real(i)-rx)/rad
				end if
				m = min( ceiling(rad), Nrot )
				vrot(m,1) = vrot(m,1) + dens
				vrot(m,2) = vrot(m,2) + vx*(rx-real(j))/rad + vy*(real(i)-rx)/rad
				vrot(m,3) = vrot(m,3) + vx2*(rx-real(j))/rad + vy2*(real(i)-rx)/rad
				rrot(m)   = rrot(m)   + rad * dens
				Nbin(m)   = Nbin(m)   + 1

				k = ceiling(sngl(nH)/sngl(Nvert))
				do n=1,Nvert
					n1 = min( lower_ind2+n*k-1, upper_ind2 )
					dens = sum( density_grid(i,j,(lower_ind2+(n-1)*k):n1) ) * real(res,8)
					vx   = sum( density_grid(i,j,(lower_ind2+(n-1)*k):n1) * vx_grid(i,j,(lower_ind2+(n-1)*k):n1) ) * real(res,8)
					vy   = sum( density_grid(i,j,(lower_ind2+(n-1)*k):n1) * vy_grid(i,j,(lower_ind2+(n-1)*k):n1) ) * real(res,8)
					vx2  = sum( vx_grid(i,j,(lower_ind2+(n-1)*k):n1) ) / (n1 - (lower_ind2+(n-1)*k) + 1)
					vy2  = sum( vy_grid(i,j,(lower_ind2+(n-1)*k):n1) ) / (n1 - (lower_ind2+(n-1)*k) + 1)
					vrot(m,3+3*n-2) = vrot(m,3+3*n-2) + dens
					vrot(m,3+3*n-1) = vrot(m,3+3*n-1) + vx*(rx-real(j))/rad + vy*(real(i)-rx)/rad
					vrot(m,3+3*n)   = vrot(m,3+3*n)   + vx2*(rx-real(j))/rad + vy2*(real(i)-rx)/rad
				end do
			end do
		end do
		print *, 'done'

		print *, 'HOW MANY ZEROS 2D: Sigma, V_{r,phi,z; mass, vol}?'
		ntot = nR*nR
		print '(2(1x,i),3(1x,es12.5))', count(mask=Sigma_gas.lt.0.0_8),count(mask=Sigma_gas.eq.0.0_8),minval(Sigma_gas,mask=Sigma_gas.ne.0.0_8),maxval(Sigma_gas),sum(Sigma_gas)/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(Vr(:,:,1)).lt.0.0_8), count(mask=abs(Vr(:,:,1)).eq.0.0_8), minval(abs(Vr(:,:,1)),mask=abs(Vr(:,:,1)).ne.0.0_8), maxval(abs(Vr(:,:,1))), sum(abs(Vr(:,:,1)))/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(Vphi(:,:,1)).lt.0.0_8), count(mask=abs(Vphi(:,:,1)).eq.0.0_8), minval(abs(Vphi(:,:,1)),mask=abs(Vphi(:,:,1)).ne.0.0_8),maxval(abs(Vphi(:,:,1))),sum(abs(Vphi(:,:,1)))/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(Vz(:,:,1)).lt.0.0_8), count(mask=abs(Vz(:,:,1)).eq.0.0_8), minval(abs(Vz(:,:,1)),mask=abs(Vz(:,:,1)).ne.0.0_8), maxval(abs(Vz(:,:,1))), sum(abs(Vz(:,:,1)))/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(Vr(:,:,2)).lt.0.0_8), count(mask=abs(Vr(:,:,2)).eq.0.0_8), minval(abs(Vr(:,:,2)),mask=abs(Vr(:,:,2)).ne.0.0_8), maxval(abs(Vr(:,:,2))), sum(abs(Vr(:,:,2)))/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(Vphi(:,:,2)).lt.0.0_8), count(mask=abs(Vphi(:,:,2)).eq.0.0_8), minval(abs(Vphi(:,:,2)),mask=abs(Vphi(:,:,2)).ne.0.0_8),maxval(abs(Vphi(:,:,2))),sum(abs(Vphi(:,:,2)))/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(Vz(:,:,2)).lt.0.0_8), count(mask=abs(Vz(:,:,2)).eq.0.0_8), minval(abs(Vz(:,:,2)),mask=abs(Vz(:,:,2)).ne.0.0_8), maxval(abs(Vz(:,:,2))), sum(abs(Vz(:,:,2)))/ntot
		print *,''
		print *, 'outputting 2d Sigma, Vphi, Vr and Vz'
		write(31) nR, nR, grid_rad, grid_rad, disc_rad, Sigma_gas, Vphi(:,:,1), Vr(:,:,1), Vz(:,:,1), Vphi(:,:,2), Vr(:,:,2), Vz(:,:,2)
		deallocate( Sigma_gas, Vz, Vr, Vphi )

		where( vrot(:,1) .gt. 0.0_8 ) vrot(:,2) = vrot(:,2)   / vrot(:,1)
		where( vrot(:,1) .gt. 0.0_8 ) rrot(:)   = rrot(:)*res / vrot(:,1)
		where( vrot(:,1) .gt. 0.0_8 ) vrot(:,1) = vrot(:,1)   / real(Nbin(:))
		where( vrot(:,1) .gt. 0.0_8 ) vrot(:,3) = vrot(:,3)   / real(Nbin(:))
		do n=1,Nvert
			where( vrot(:,3+3*n-2) .gt. 0.0_8 ) vrot(:,3+3*n-1) = vrot(:,3+3*n-1) / vrot(:,3+3*n-2)
			where( Nbin(:) .gt. 0 ) vrot(:,3+3*n) = vrot(:,3+3*n) / real(Nbin(:))
		end do
		if( type_2_FWHM .gt. 0.0_8 ) then
			Hw = 0.5_8 * type_2_FWHM
			allocate( test(Nrot,3+3*Nvert), stat=i )
			if(i.ne.0) then
				print *, 'error in allocation of Nbin. stat= ', i
				stop
			end if
			print *, 'smoothing rotation curve'
			test(:,:) = vrot(:,:)
			do j=1,Nrot
				norm = sum( 0.5_8**((abs(rrot(:)-rrot(j))/Hw)**2), mask = abs(rrot(:)-rrot(j)) .le. 5.0_8*Hw )
				do i=1,3
					vrot(j,i)  = sum( test(:,i)  * 0.5_8**((abs(rrot(:)-rrot(j))/Hw)**2), mask = abs(rrot(:)-rrot(j)) .le. 5.0_8*Hw ) / norm
				end do
				do i=1,Nvert
					vrot(j,3+3*i-1) = sum( test(:,3+3*i-1) * 0.5_8**((abs(rrot(:)-rrot(j))/Hw)**2), mask = abs(rrot(:)-rrot(j)) .le. 5.0_8*Hw ) / norm
					vrot(j,3+3*i)   = sum( test(:,3+3*i)  * 0.5_8**((abs(rrot(:)-rrot(j))/Hw)**2), mask = abs(rrot(:)-rrot(j)) .le. 5.0_8*Hw ) / norm
				end do
			end do
			deallocate( test )
		end if
		write(filename,'(a,a,a,a)') trim(output_dirname),'/rotation_curves/',trim(snap_tag2),'_face_on_profile.out'
		open(unit=21,file=filename)
		do i=1,Nrot
			write(21,'(4(1x,es12.5))') rrot(i), vrot(i,1), vrot(i,2), vrot(i,3)
		end do
		write(21,'(4(1x,es12.5))') disc_rad, 1000.0_4 * Rdisc(snapshot), 1300.0_4 * Rdisc(snapshot), 1000.0_4 * Hdisc(snapshot)
		close(unit=21)
		write(filename,'(a,a,a,a)') trim(output_dirname),'/rotation_curves/',trim(snap_tag2),'_rotation_curves.out'
		open(unit=22,file=filename)
		if( 2*Nvert+3 .lt. 9 ) then
			write(filename,'(a,i1,a)') '(',2*Nvert+3,'(1x,es12.5))'
		elseif( 2*Nvert+3 .lt. 99 ) then
			write(filename,'(a,i2,a)') '(',2*Nvert+3,'(1x,es12.5))'
		else
			write(filename,'(a,i3,a)') '(',2*Nvert+3,'(1x,es12.5))'
		end if
		do i=1,Nrot
			write(22,trim(filename)) res*(sngl(i)-0.5_4), vrot(i,2), vrot(i,3), vrot(i,3+3*(/ (j,j=1,Nvert) /)-1), vrot(i,3+3*(/ (j,j=1,Nvert) /))
		end do
		close(unit=22)
		rrot(:) = rrot(:) / res

		print *, '!!! disc rad, height [kpc] !!!'
		print *, disc_rad / 1000.0_4, disc_height / 1000.0_4

		if( decomp_type .eq. 2 ) then
			print *, 'Decomp type 2'
			print *, 'Subtracting the rotation curve from the velocity field'
			print *, 'max and min Vrot within and outside buffer'
			print '(4(1x,es12.5))', maxval(vrot(1:nR_half,2)),minval(vrot(1:nR_half,2)),maxval(vrot(nR_half+1:Nrot,2)),minval(vrot(nR_half+1:Nrot,2))
			print '(4(1x,es12.5))', maxval(vrot(1:nR_half,3)),minval(vrot(1:nR_half,3)),maxval(vrot(nR_half+1:Nrot,3)),minval(vrot(nR_half+1:Nrot,3))
			n1 = ceiling(sngl(nH)/sngl(Nvert))
			if( type_2_mass .eq. 1 ) then
				n2 = 1
			else
				n2 = 0
			end if
			do j=1,ngrid
				do i=1,ngrid
					rad = sqrt( (real(i)-rx)**2 + (real(j)-rx)**2 )
					m = min( ceiling(rad), Nrot )
					do k=1,ngrid2
						n = max( 1, ceiling(real(k - lower_ind2)/real(n1)) )
						n = min( n, Nvert )
						n = 3 + 3*n - n2
						if( m .eq. 1 .or. m .eq. Nrot ) then
							vrot_interp = vrot(m,n)
						else
							if( rad .lt. rrot(m) ) then
								vrot_interp = vrot(m-1,n) + ( rad - rrot(m-1) ) * ( vrot(m,n) - vrot(m-1,n) ) / ( rrot(m) - rrot(m-1) )
							else
								vrot_interp = vrot(m,n)   + ( rad - rrot(m) )   * ( vrot(m+1,n) - vrot(m,n) ) / ( rrot(m+1) - rrot(m) )
							end if
						end if
						vx_grid(i,j,k) = vx_grid(i,j,k) - vrot_interp * real(rx-j)/rad
						vy_grid(i,j,k) = vy_grid(i,j,k) - vrot_interp * real(i-rx)/rad
					end do
				end do
			end do
			print *, 'HOW MANY ZEROS 3D, after removing rotation: density, v{x,y,z}?'
			ntot = ngrid*ngrid*ngrid2
			print '(2(1x,i),3(1x,es12.5))', count(mask=density_grid.lt.0.0_8),count(mask=density_grid.eq.0.0_8),minval(density_grid,mask=density_grid.ne.0.0_8),maxval(density_grid),sum(density_grid)/ntot
			print '(2(1x,i),3(1x,es12.5))', count(mask=abs(vx_grid).lt.0.0_8),count(mask=abs(vx_grid).eq.0.0_8),minval(abs(vx_grid),mask=abs(vx_grid).ne.0.0_8),maxval(abs(vx_grid)),sum(abs(vx_grid))/ntot
			print '(2(1x,i),3(1x,es12.5))', count(mask=abs(vy_grid).lt.0.0_4),count(mask=abs(vy_grid).eq.0.0_4),minval(abs(vy_grid),mask=abs(vy_grid).ne.0.0_4),maxval(abs(vy_grid)),sum(abs(vy_grid))/ntot
			print '(2(1x,i),3(1x,es12.5))', count(mask=abs(vz_grid).lt.0.0_4),count(mask=abs(vz_grid).eq.0.0_4),minval(abs(vz_grid),mask=abs(vz_grid).ne.0.0_4),maxval(abs(vz_grid)),sum(abs(vz_grid))/ntot
			print *,''
		end if
		deallocate( vrot, rrot, Nbin )
	end subroutine rot_curve_surface_maps
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine turb_def()
	! Subtractes from the momentum gids the average momentum smoothed over a local volume with either a top hat or Gaussian filter, in 2D or 3D
	! This defines the turbulence grids
	! These may be then be smoothed, if selected
	! Then creates 2D face on projections of cylindrical turbulence components and total turbulence magnitude
		implicit none
		integer :: i, j, k, m
		integer :: lower_ind, upper_ind, lower_ind2, upper_ind2, ntot
		real(8) :: mind, rad, rx, vx, vy, vx2, vy2, dens
		real(8),allocatable :: sigz(:,:,:), sigr(:,:,:), sigp(:,:,:), sigt(:,:,:)
		character(len=256) :: filename

		lower_ind  = ngrid/2  - nR/2 + 1
		lower_ind2 = ngrid2/2 - nH/2 + 1
		upper_ind  = ngrid/2  + nR/2
		upper_ind2 = ngrid2/2 + nH/2

		if( decomp_type .eq. 3 ) then
			print *, 'Decomp type 3'
			print *, 'Subtracting from each cell the velocity averaged over some scale, Gaussian weighted'
			print *, 'Equivalent to filtering out large scales'
			print '(a5,i4.4,a9)', 'FWHM=', type_3_FWHM, 'pc, in 3D'
			mind = minval(density_grid, mask = density_grid .gt. 0.0_8)
			call fft_smooth_3D(type_3_FWHM / res, density_grid, ngrid, ngrid, ngrid2, 0, 2)
			call fft_smooth_3D(type_3_FWHM / res, vx_grid, ngrid, ngrid, ngrid2, 0, 2)
			call fft_smooth_3D(type_3_FWHM / res, vy_grid, ngrid, ngrid, ngrid2, 0, 2)
			call fft_smooth_3D(type_3_FWHM / res, vz_grid, ngrid, ngrid, ngrid2, 0, 2)
			density_grid = max(density_grid, mind)
			print *, 'Done with local subtraction'

			print *, 'HOW MANY ZEROS 3D, after filtering average: density, v_{x,y,z}?'
			ntot = ngrid*ngrid*ngrid2
			print '(2(1x,i),3(1x,es12.5))', count(mask=density_grid.lt.0.0_8),count(mask=density_grid.eq.0.0_8),minval(density_grid,mask=density_grid.ne.0.0_8),maxval(density_grid),sum(density_grid)/ntot
			print '(2(1x,i),3(1x,es12.5))', count(mask=abs(vx_grid).lt.0.0_8),count(mask=abs(vx_grid).eq.0.0_8),minval(abs(vx_grid),mask=abs(vx_grid).ne.0.0_8),maxval(abs(vx_grid)),sum(abs(vx_grid))/ntot
			print '(2(1x,i),3(1x,es12.5))', count(mask=abs(vy_grid).lt.0.0_8),count(mask=abs(vy_grid).eq.0.0_8),minval(abs(vy_grid),mask=abs(vy_grid).ne.0.0_8),maxval(abs(vy_grid)),sum(abs(vy_grid))/ntot
			print '(2(1x,i),3(1x,es12.5))', count(mask=abs(vz_grid).lt.0.0_8),count(mask=abs(vz_grid).eq.0.0_8),minval(abs(vz_grid),mask=abs(vz_grid).ne.0.0_8),maxval(abs(vz_grid)),sum(abs(vz_grid))/ntot
			print *,''
		end if

		allocate( sigz(nR,nR,2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of sigz. stat= ', i
			stop
		end if
		allocate( sigr(nR,nR,2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of sigr. stat= ', i
			stop
		end if
		allocate( sigp(nR,nR,2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of sigphi. stat= ', i
			stop
		end if
		allocate( sigt(nR,nR,2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of sigtot. stat= ', i
			stop
		end if

		print *, 'calculating 2d projections'
		rx = real(nR+1)/2.0_8
		do j=1,nR
			do i=1,nR
				rad = sqrt( (real(i)-rx)**2 + (real(j)-rx)**2 )
				k = lower_ind+i-1
				m = lower_ind+j-1
				dens = sum( density_grid(k,m,lower_ind2:upper_ind2) )
				dens = max(dens, 1.d-10)
				vx   = sum( density_grid(k,m,lower_ind2:upper_ind2) * vx_grid(k,m,lower_ind2:upper_ind2) ) / dens
				vy   = sum( density_grid(k,m,lower_ind2:upper_ind2) * vy_grid(k,m,lower_ind2:upper_ind2) ) / dens
				vx2  = sum( vx_grid(k,m,lower_ind2:upper_ind2) ) / nH
				vy2  = sum( vy_grid(k,m,lower_ind2:upper_ind2) ) / nH
				sigz(i,j,1) = sum( density_grid(k,m,lower_ind2:upper_ind2) * vz_grid(k,m,lower_ind2:upper_ind2) ) / dens
				sigr(i,j,1) = vx*(real(i)-rx)/rad + vy*(real(j)-rx)/rad
				sigp(i,j,1) = vx*(rx-real(j))/rad + vy*(real(i)-rx)/rad
				sigt(i,j,1) = sum( density_grid(k,m,lower_ind2:upper_ind2) * ( & 
					& vx_grid(k,m,lower_ind2:upper_ind2)**2 + vy_grid(k,m,lower_ind2:upper_ind2)**2 + vz_grid(k,m,lower_ind2:upper_ind2)**2 ) ) / dens
				sigt(i,j,1) = sqrt( sigt(i,j,1) )
				sigz(i,j,2) = sum( vz_grid(k,m,lower_ind2:upper_ind2) ) / nH
				sigr(i,j,2) = vx2*(real(i)-rx)/rad + vy2*(real(j)-rx)/rad
				sigp(i,j,2) = vx2*(rx-real(j))/rad + vy2*(real(i)-rx)/rad
				sigt(i,j,2) = sqrt( sum( vx_grid(k,m,lower_ind2:upper_ind2)**2 + vy_grid(k,m,lower_ind2:upper_ind2)**2 + vz_grid(k,m,lower_ind2:upper_ind2)**2 ) / nH )
			end do
		end do
		print *, 'HOW MANY ZEROS 2D, within buffer: sigma_{r,phi,z,tot; mass,vol}'
		ntot = nR*nR
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(sigr(:,:,1)).lt.0.0_8),count(mask=abs(sigr(:,:,1)).eq.0.0_8),minval(abs(sigr(:,:,1)),mask=abs(sigr(:,:,1)).ne.0.0_8),maxval(abs(sigr(:,:,1))),sum(abs(sigr(:,:,1)))/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(sigp(:,:,1)).lt.0.0_8),count(mask=abs(sigp(:,:,1)).eq.0.0_8),minval(abs(sigp(:,:,1)),mask=abs(sigp(:,:,1)).ne.0.0_8),maxval(abs(sigp(:,:,1))),sum(abs(sigp(:,:,1)))/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(sigz(:,:,1)).lt.0.0_8),count(mask=abs(sigz(:,:,1)).eq.0.0_8),minval(abs(sigz(:,:,1)),mask=abs(sigz(:,:,1)).ne.0.0_8),maxval(abs(sigz(:,:,1))),sum(abs(sigz(:,:,1)))/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(sigt(:,:,1)).lt.0.0_8),count(mask=abs(sigt(:,:,1)).eq.0.0_8),minval(abs(sigt(:,:,1)),mask=abs(sigt(:,:,1)).ne.0.0_8),maxval(abs(sigt(:,:,1))),sum(abs(sigt(:,:,1)))/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(sigr(:,:,2)).lt.0.0_8),count(mask=abs(sigr(:,:,2)).eq.0.0_8),minval(abs(sigr(:,:,2)),mask=abs(sigr(:,:,2)).ne.0.0_8),maxval(abs(sigr(:,:,2))),sum(abs(sigr(:,:,2)))/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(sigp(:,:,2)).lt.0.0_8),count(mask=abs(sigp(:,:,2)).eq.0.0_8),minval(abs(sigp(:,:,2)),mask=abs(sigp(:,:,2)).ne.0.0_8),maxval(abs(sigp(:,:,2))),sum(abs(sigp(:,:,2)))/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(sigz(:,:,2)).lt.0.0_8),count(mask=abs(sigz(:,:,2)).eq.0.0_8),minval(abs(sigz(:,:,2)),mask=abs(sigz(:,:,2)).ne.0.0_8),maxval(abs(sigz(:,:,2))),sum(abs(sigz(:,:,2)))/ntot
		print '(2(1x,i),3(1x,es12.5))', count(mask=abs(sigt(:,:,2)).lt.0.0_8),count(mask=abs(sigt(:,:,2)).eq.0.0_8),minval(abs(sigt(:,:,2)),mask=abs(sigt(:,:,2)).ne.0.0_8),maxval(abs(sigt(:,:,2))),sum(abs(sigt(:,:,2)))/ntot

		print *, 'outputting 2d turbulence maps'
		write(32) nR, nR, grid_rad, grid_rad, disc_rad, sigp(:,:,1), sigr(:,:,1), sigz(:,:,1), sigt(:,:,1), sigp(:,:,2), sigr(:,:,2), sigz(:,:,2), sigt(:,:,2)
		deallocate( sigz, sigr, sigp, sigt )

	end subroutine turb_def
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine turb_decomp(snapshot)
	! Final part
	! decomposes turbulence into compreesive and solenoidal components, calculates the global ratios, and creates projected 2D plots of them
	! snapshot is snapshot number from main loop, needed for io purposes
		implicit none
		integer,intent(in) :: snapshot
		integer :: i, j, k, m, m1, m2, m3
		integer :: lower_ind, upper_ind, lower_ind2, upper_ind2, nghalf, ng2half, ntot
		real(8) :: mnorm, nnorm, rad, vx, vy, vx2, vy2, dens, rx, ry, rz
		real(8),allocatable :: sigz(:,:,:), sigr(:,:,:), sigp(:,:,:), sigt(:,:,:), Sigc_2d(:,:,:), Sigs_2d(:,:,:), dot_2d(:,:,:), test(:,:)
		character(len=256) :: filename

		nghalf  = ngrid/2
		ng2half = ngrid2/2
		lower_ind  = nghalf  - nR/2 + 1
		lower_ind2 = ng2half - nH/2 + 1
		upper_ind  = nghalf  + nR/2
		upper_ind2 = ng2half + nH/2

		print *, ''
		print *, 'Starting decomposition'
		if( abs(decomp_weight) .gt. 0.1_8 ) then
			vx_grid = (density_grid**decomp_weight) * vx_grid
			vy_grid = (density_grid**decomp_weight) * vy_grid
			vz_grid = (density_grid**decomp_weight) * vz_grid
		end if

		call fft_decompose
		call local_decompose
		deallocate( vx_grid, vy_grid, vz_grid )

		print *, 'making surface plots of turbulence modes'
		allocate( Sigc_2d(nR,nR,4), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of Sigc_2d. stat= ', i
			stop
		end if
		allocate( Sigs_2d(nR,nR,2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of Sigs_2d. stat= ', i
			stop
		end if
		do j=1,nR
			do i=1,nR
				k = lower_ind+i-1
				m = lower_ind+j-1
				dens = sum( density_grid(k,m,lower_ind2:upper_ind2) )
				dens = max(dens, 1.d-10)
				Sigc_2d(i,j,1) = sum( density_grid(k,m,lower_ind2:upper_ind2) * sigcx(k,m,lower_ind2:upper_ind2) ) / dens		!!! ( km sec^{-1} )^2
				Sigc_2d(i,j,2) = sum( density_grid(k,m,lower_ind2:upper_ind2) * sigcy(k,m,lower_ind2:upper_ind2) ) / dens		!!! ( km sec^{-1} )^2
				Sigs_2d(i,j,1) = sum( density_grid(k,m,lower_ind2:upper_ind2) * sigsx(k,m,lower_ind2:upper_ind2) ) / dens		!!! ( km sec^{-1} )^2
				Sigc_2d(i,j,3) = sum( sigcx(k,m,lower_ind2:upper_ind2) ) / real(nH)							!!! ( km sec^{-1} )^2
				Sigc_2d(i,j,4) = sum( sigcy(k,m,lower_ind2:upper_ind2) ) / real(nH)							!!! ( km sec^{-1} )^2
				Sigs_2d(i,j,2) = sum( sigsx(k,m,lower_ind2:upper_ind2) ) / real(nH)							!!! ( km sec^{-1} )^2
			end do
		end do

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		print *, 'file 1'
		Sigs_2d(:,:,1) = sqrt(Sigs_2d(:,:,1))
		Sigc_2d(:,:,1) = sqrt(Sigc_2d(:,:,1))
		Sigc_2d(:,:,2) = sqrt(Sigc_2d(:,:,2))
		write(34) nR, nR, grid_rad, grid_rad, disc_rad, Sigs_2d(:,:,1), Sigc_2d(:,:,1), Sigc_2d(:,:,2)
		Sigs_2d(:,:,1) = Sigs_2d(:,:,1)**2
		Sigc_2d(:,:,1) = Sigc_2d(:,:,1)**2
		Sigc_2d(:,:,2) = Sigc_2d(:,:,2)**2

		print *, 'file 2'
		Sigs_2d(:,:,2) = sqrt(Sigs_2d(:,:,2))
		Sigc_2d(:,:,3) = sqrt(Sigc_2d(:,:,3))
		Sigc_2d(:,:,4) = sqrt(Sigc_2d(:,:,4))
		write(35) nR, nR, grid_rad, grid_rad, disc_rad, Sigs_2d(:,:,2), Sigc_2d(:,:,3), Sigc_2d(:,:,4)
		Sigs_2d(:,:,2) = Sigs_2d(:,:,2)**2
		Sigc_2d(:,:,3) = Sigc_2d(:,:,3)**2
		Sigc_2d(:,:,4) = Sigc_2d(:,:,4)**2

		print *, 'file 3'
		allocate( test(nR,nR), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of test. stat= ', i
			stop
		end if
		test = Sigs_2d(:,:,1)+Sigc_2d(:,:,1)+Sigc_2d(:,:,2)
		Sigs_2d(:,:,1) = Sigs_2d(:,:,1) / test
		Sigc_2d(:,:,1) = Sigc_2d(:,:,1) / test
		Sigc_2d(:,:,2) = Sigc_2d(:,:,2) / test
		write(36) nR, nR, grid_rad, grid_rad, disc_rad, Sigs_2d(:,:,1), Sigc_2d(:,:,1), Sigc_2d(:,:,2)
		Sigs_2d(:,:,1) = Sigs_2d(:,:,1) * test
		Sigc_2d(:,:,1) = Sigc_2d(:,:,1) * test
		Sigc_2d(:,:,2) = Sigc_2d(:,:,2) * test

		print *, 'file 4'
		test = Sigs_2d(:,:,2)+Sigc_2d(:,:,3)+Sigc_2d(:,:,4)
		Sigs_2d(:,:,2) = Sigs_2d(:,:,2) / test
		Sigc_2d(:,:,3) = Sigc_2d(:,:,3) / test
		Sigc_2d(:,:,4) = Sigc_2d(:,:,4) / test
		write(37) nR, nR, grid_rad, grid_rad, disc_rad, Sigs_2d(:,:,2), Sigc_2d(:,:,3), Sigc_2d(:,:,4)
		print *, 'done with 2d maps'
		Sigs_2d(:,:,2) = Sigs_2d(:,:,2) * test
		Sigc_2d(:,:,3) = Sigc_2d(:,:,3) * test
		Sigc_2d(:,:,4) = Sigc_2d(:,:,4) * test

		nnorm = sum( sigcx + sigcy + sigsx )
		mnorm = sum( density_grid * (sigcx + sigcy + sigsx) )
		sigma_ratio(3) = sum( sigsx ) / nnorm
		sigma_ratio(4) = sum( sigcx ) / nnorm
		sigma_ratio(5) = sum( sigcy ) / nnorm
		sigma_ratio(6) = sum( density_grid * sigsx ) / mnorm
		sigma_ratio(7) = sum( density_grid * sigcx ) / mnorm
		sigma_ratio(8) = sum( density_grid * sigcy ) / mnorm

		nnorm = 0.0_8
		mnorm = 0.0_8
		sigma_ratio(9:14) = 0.0_8
		rx = real(ngrid+1)/2.0_8
		ry = real(ngrid+1)/2.0_8
		rz = real(ngrid2+1)/2.0_8
		do k=1,nH
			do j=1,nR
				do i=1,nR
					m1 = lower_ind  + i - 1
					m2 = lower_ind  + j - 1
					m3 = lower_ind2 + k - 1
					if( sqrt( (real(m1)-rx)**2 + (real(m2)-ry)**2 ) * res .le. disc_Rad .and. abs( real(m3)-ng2half ) * res .le. disc_height ) then
						nnorm = nnorm + sigcx(m1,m2,m3) + sigcy(m1,m2,m3) + sigsx(m1,m2,m3)
						mnorm = mnorm + density_grid(m1,m2,m3) * ( sigcx(m1,m2,m3) + sigcy(m1,m2,m3) + sigsx(m1,m2,m3) )
						sigma_ratio(9)  = sigma_ratio(9)  + sigsx(m1,m2,m3)
						sigma_ratio(10) = sigma_ratio(10) + sigcx(m1,m2,m3)
						sigma_ratio(11) = sigma_ratio(11) + sigcy(m1,m2,m3)
						sigma_ratio(12) = sigma_ratio(12) + density_grid(m1,m2,m3) * sigsx(m1,m2,m3)
						sigma_ratio(13) = sigma_ratio(13) + density_grid(m1,m2,m3) * sigcx(m1,m2,m3)
						sigma_ratio(14) = sigma_ratio(14) + density_grid(m1,m2,m3) * sigcy(m1,m2,m3)
					end if
				end do
			end do
		end do
		sigma_ratio(9:11)  = sigma_ratio(9:11)  / nnorm
		sigma_ratio(12:14) = sigma_ratio(12:14) / mnorm

		print *, 'volume ratios, whole box: s, e, c'
		print '(4(1x,es12.5))', sigma_ratio(3), sigma_ratio(4), sigma_ratio(5)
		print *, 'mass ratios, whole box: s, e, c'
		print '(4(1x,es12.5))', sigma_ratio(6), sigma_ratio(7), sigma_ratio(8)
		print *, 'volume ratios, disc: s, e, c'
		print '(4(1x,es12.5))', sigma_ratio(9), sigma_ratio(10), sigma_ratio(11)
		print *, 'mass ratios, disc: s, e, c'
		print '(4(1x,es12.5))', sigma_ratio(12), sigma_ratio(13), sigma_ratio(14)
		deallocate( density_grid, sigcx, sigcy, sigsx )

		write(20,'(1x,f5.3,14(1x,es12.5))') aexp(snapshot), sigma_ratio(1:2), sigma_ratio(3:5), sigma_ratio(9:11), sigma_ratio(6:8), sigma_ratio(12:14)

		deallocate( Sigc_2d, Sigs_2d, test )
	end subroutine turb_decomp
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine fft_decompose
	use MKL_DFTI
	! decomposes turbulence into solenoidal and compressive modes in Fourier space

		implicit none
		real(8),allocatable :: kx(:,:,:), ky(:,:,:), kz(:,:,:), kmag(:,:,:), Skmag(:,:,:)
		real(8),allocatable :: kvec(:), PS_out(:,:)
		complex(8),allocatable :: fft_inx(:), fft_iny(:), fft_inz(:)
		complex(8),allocatable :: Gcx_1d(:),  Gcy_1d(:), Gcz_1d(:), Gsx_1d(:),  Gsy_1d(:), Gsz_1d(:)
		complex(8) :: im, Gk
		real(8) :: L, L2, dk1, dk2, normc, norms
		integer :: i, j, m, n, n1, n2, n3, nvec, error, length(3) 
		integer :: Nxy, Nxyz, nghalf, ng2half, lower_ind, upper_ind, lower_ind2, upper_ind2
		character(len=256) :: filename
		type(DFTI_DESCRIPTOR), POINTER :: fft_desc

		Nxy = ngrid*ngrid
		Nxyz = Nxy*ngrid2

		nghalf     = ngrid/2
		ng2half    = ngrid2/2
		lower_ind  = nghalf  - nR/2 + 1
		lower_ind2 = ng2half - nH/2 + 1
		upper_ind  = nghalf  + nR/2
		upper_ind2 = ng2half + nH/2

		im = (0.0_8, 1.0_8)
		L  = 2.0_8 * xmax
		L2 = 2.0_8 * zmax
		dk1 = 2.0_8*pi/L
		dk2 = 2.0_8*pi/L2
		length = (/ ngrid, ngrid, ngrid2 /)

		print *, 'allocating K grids'
		allocate( kx(ngrid,ngrid,ngrid2), ky(ngrid,ngrid,ngrid2), kz(ngrid,ngrid,ngrid2), kmag(ngrid,ngrid,ngrid2) )
		kx = 0.0_8
		ky = 0.0_8
		kz = 0.0_8
		kmag = 0.0_8
		do m=1,ngrid2
			do j=1,ngrid
				kx(1:nghalf, j, m)       = real( (/ (i,i=1,nghalf) /) - 1 )
				kx((nghalf+1):ngrid, j, m) = real( (/ (i,i=1,nghalf) /) - 1 - nghalf )
			end do
		end do
		do m=1,ngrid2
			do i=1,ngrid
				ky(i, 1:nghalf, m)       = real( (/ (j,j=1,nghalf) /) - 1 )
				ky(i, (nghalf+1):ngrid, m) = real( (/ (j,j=1,nghalf) /) - 1 - nghalf )
			end do
		end do
		do j=1,ngrid
			do i=1,ngrid
				kz(i, j, 1:ng2half)       = real( (/ (m,m=1,ng2half) /) - 1 )
				kz(i, j, (ng2half+1):ngrid2) = real( (/ (m,m=1,ng2half) /) - 1 - ng2half )
			end do
		end do
		kx = kx * dk1
		ky = ky * dk1
		kz = kz * dk2
		kmag = sqrt( kx**2 + ky**2 + kz**2 )

		allocate( fft_inx(Nxyz), fft_iny(Nxyz), fft_inz(Nxyz) )
		do m=1,ngrid2
			do j=1,ngrid
				do i=1,ngrid
					fft_inx( (m-1)*Nxy + (j-1)*ngrid + i ) = vx_grid(i,j,m)
					fft_iny( (m-1)*Nxy + (j-1)*ngrid + i ) = vy_grid(i,j,m)
					fft_inz( (m-1)*Nxy + (j-1)*ngrid + i ) = vz_grid(i,j,m)
				end do
			end do
		end do
		norms = sum(abs(fft_inx)**2 + abs(fft_iny)**2 + abs(fft_inz)**2)

		print *, 'Forward FFT-ing'
		error = DftiCreateDescriptor(fft_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc)
		error = DftiComputeForward(fft_desc, fft_inx)
		error = DftiFreeDescriptor(fft_desc)

		error = DftiCreateDescriptor(fft_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc)
		error = DftiComputeForward(fft_desc, fft_iny)
		error = DftiFreeDescriptor(fft_desc)

		error = DftiCreateDescriptor(fft_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc)
		error = DftiComputeForward(fft_desc, fft_inz)
		error = DftiFreeDescriptor(fft_desc)

		fft_inx = fft_inx / sqrt(real(Nxyz))
		fft_iny = fft_iny / sqrt(real(Nxyz))
		fft_inz = fft_inz / sqrt(real(Nxyz))
		normc = sum(abs(fft_inx)**2 + abs(fft_iny)**2 + abs(fft_inz)**2)
		print '(a9,3(1x,es12.5))', 'PARCEVAL ', norms, normc, norms/normc

		allocate( Skmag(ngrid,ngrid,ngrid2) )
		Skmag = kmag
		kx = sin(kx*res) / res
		ky = sin(ky*res) / res
		kz = sin(kz*res) / res
		kmag = sqrt( kx**2 + ky**2 + kz**2 )
		allocate( Gcx_1d(Nxyz), Gcy_1d(Nxyz), Gcz_1d(Nxyz), Gsx_1d(Nxyz), Gsy_1d(Nxyz), Gsz_1d(Nxyz) )
		! compressive part
		do m=1,ngrid2
			do j=1,ngrid
				do i=1,ngrid
					if( kmag(i,j,m) .eq. 0.0_8 ) then
						Gcx_1d = 0.0_8
						Gcy_1d = 0.0_8
						Gcz_1d = 0.0_8
					else
						n = (m-1)*Nxy + (j-1)*ngrid + i
						Gk = ( kx(i,j,m)*fft_inx(n) + ky(i,j,m)*fft_iny(n) + kz(i,j,m)*fft_inz(n) ) / max(kmag(i,j,m)**2, 1.d-16)
						Gcx_1d(n) = kx(i,j,m) * Gk
						Gcy_1d(n) = ky(i,j,m) * Gk
						Gcz_1d(n) = kz(i,j,m) * Gk
					end if
				end do
			end do
		end do
		deallocate( kx, ky, kz )
		! solenoidal part
		Gsx_1d = fft_inx - Gcx_1d
		Gsy_1d = fft_iny - Gcy_1d
		Gsz_1d = fft_inz - Gcz_1d
		kmag = Skmag
		deallocate( Skmag )

		print *, 'Power spectrum'
	!!! LEAVE COMMENTED TO HAVE UNIFORM BINNING IN K SPACE, WITH CELL SIZE CALIBRATED ON DISC DIMAETER !!!
	!!! UNCOMMENT TO HAVE NON-UNIFORM BINNING IN K SPACE, WITH JUMP AT DISC THICKNESS !!!
		n1 = ceiling( dk2/dk1 )	! = ngrid/ngrid2 = L/L2 = xmax/zmax
		n2 = ceiling( sqrt(3.0_8)*(pi/res) / dk2 ) + 1
		nvec = n1 + n2 - 1
		print *, 'n1, n2, nvec ',n1, n2, nvec
!		nvec = ceiling( sqrt(3.0_8)*(pi/res) / dk1 ) + 1

		allocate( kvec(nvec), PS_out(nvec,13) )

		kvec(1:n1)  = (/ (i,i=1,n1) /)*dk1
		kvec(n1:(n1+n2-1)) = (/ (i,i=1,n2) /)*dk2
		kvec(1:n1) = kvec(1:n1) - 0.5_8*dk1
		kvec((n1+1):nvec) = kvec((n1+1):nvec) - 0.5_8*dk2
!		kvec(1:nvec)  = (/ (i,i=1,nvec) /)*dk1 - 0.5_8*dk1
		PS_out = 0.0_8
		do m=1,ngrid2
			do j=1,ngrid
				do i=1,ngrid
					if( kmag(i,j,m) .ge. dk2 ) then
						n3 = min( floor( kmag(i,j,m)/dk2 ) + n1, nvec )
					else
						n3 = min( floor( kmag(i,j,m)/dk1 ) + 1, nvec )
					end if
					n = (m-1)*Nxy + (j-1)*ngrid + i
					PS_out(n3,1)  = PS_out(n3,1)  + abs(fft_inx(n))**2 + abs(fft_iny(n))**2 + abs(fft_inz(n))**2
					PS_out(n3,2)  = PS_out(n3,2)  + abs(Gsx_1d(n))**2  + abs(Gsy_1d(n))**2  + abs(Gsz_1d(n))**2
					PS_out(n3,3)  = PS_out(n3,3)  + abs(Gcx_1d(n))**2  + abs(Gcy_1d(n))**2  + abs(Gcz_1d(n))**2
					PS_out(n3,4)  = PS_out(n3,4)  + fft_inx(n)
					PS_out(n3,5)  = PS_out(n3,5)  + fft_iny(n)
					PS_out(n3,6)  = PS_out(n3,6)  + fft_inz(n)
					PS_out(n3,7)  = PS_out(n3,7)  + Gsx_1d(n)
					PS_out(n3,8)  = PS_out(n3,8)  + Gsy_1d(n)
					PS_out(n3,9)  = PS_out(n3,9)  + Gsz_1d(n)
					PS_out(n3,10) = PS_out(n3,10) + Gcx_1d(n)
					PS_out(n3,11) = PS_out(n3,11) + Gcy_1d(n)
					PS_out(n3,12) = PS_out(n3,12) + Gcz_1d(n)
					PS_out(n3,13) = PS_out(n3,13) + 1.0_8
				end do
			end do
		end do
		deallocate( fft_inx, fft_iny, fft_inz )
		write(filename,'(a,a,a,a)') trim(output_dirname),'/power_spectra/',trim(snap_tag2),'.out'
		open(unit=21,file=filename)
		do i=1,nvec
			write(21,'(7(1x,es12.5))') kvec(i), real(PS_out(i,1)), real(PS_out(i,2)), real(PS_out(i,3)), & 
				& abs(PS_out(i,4)/PS_out(i,13))**2  + abs(PS_out(i,5)/PS_out(i,13))**2  + abs(PS_out(i,6)/PS_out(i,13))**2, &
				& abs(PS_out(i,7)/PS_out(i,13))**2  + abs(PS_out(i,8)/PS_out(i,13))**2  + abs(PS_out(i,9)/PS_out(i,13))**2, & 
				& abs(PS_out(i,10)/PS_out(i,13))**2 + abs(PS_out(i,11)/PS_out(i,13))**2 + abs(PS_out(i,12)/PS_out(i,13))**2
		end do
		close(unit=21)
		print *, 'done'
		deallocate( kvec, PS_out )

		print *, 'Filter out largest and smallest modes'
		do m=1,ngrid2
			do j=1,ngrid
				do i=1,ngrid
					if( kmag(i,j,m) .lt. dk2 .or. kmag(i,j,m) .ge. pi/res ) then
						n = (m-1)*Nxy + (j-1)*ngrid + i
						Gsx_1d(n) = 0.0_8
						Gsy_1d(n) = 0.0_8
						Gsz_1d(n) = 0.0_8
						Gcx_1d(n) = 0.0_8
						Gcy_1d(n) = 0.0_8
						Gcz_1d(n) = 0.0_8
					end if
				end do
			end do
		end do
		norms = sum(abs(Gsx_1d)**2 + abs(Gsy_1d)**2 + abs(Gsz_1d)**2)
		normc = sum(abs(Gcx_1d)**2 + abs(Gcy_1d)**2 + abs(Gcz_1d)**2)
		deallocate( kmag )
		deallocate( Gcx_1d, Gcy_1d, Gcz_1d, Gsx_1d, Gsy_1d, Gsz_1d )

		print *, 'total power fourier space: t, s/t c/t'
		print '(3(1x,es12.5))', norms+normc, norms/(norms+normc), normc/(norms+normc)
		sigma_ratio(1) = norms/(norms+normc)
		sigma_ratio(2) = normc/(norms+normc)

	end subroutine fft_decompose
!_____________________________________________________________________________________________________________________________________________________________________
	subroutine local_decompose
	! decomposes turbulence into solenoidal and compressive modes locally in real space (using strain rate tensor)

		implicit none
		real(8) :: divxx, divxy, divxz, divyx, divyy, divyz, divzx, divzy, divzz
		integer :: i, j, m

		print *, 'allocating grids for compressive and solenoidal turbulence'
		allocate( sigcx(ngrid,ngrid,ngrid2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of sigcx. stat= ', i
			stop
		end if
		allocate( sigcy(ngrid,ngrid,ngrid2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of sigcy. stat= ', i
			stop
		end if
		allocate( sigsx(ngrid,ngrid,ngrid2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation of sigsx. stat= ', i
			stop
		end if
		print *, 'OK'
		sigcx = 1.d-10
		sigcy = 1.d-10
		sigsx = 1.d-10

		do m=2,ngrid2-1
			do j=2,ngrid-1
				do i=2,ngrid-1
					divxx = (vx_grid(i+1,j,m) - vx_grid(i-1,j,m)) / 2.0_8
					divxy = (vx_grid(i,j+1,m) - vx_grid(i,j-1,m)) / 2.0_8
					divxz = (vx_grid(i,j,m+1) - vx_grid(i,j,m-1)) / 2.0_8
					divyx = (vy_grid(i+1,j,m) - vy_grid(i-1,j,m)) / 2.0_8
					divyy = (vy_grid(i,j+1,m) - vy_grid(i,j-1,m)) / 2.0_8
					divyz = (vy_grid(i,j,m+1) - vy_grid(i,j,m-1)) / 2.0_8
					divzx = (vz_grid(i+1,j,m) - vz_grid(i-1,j,m)) / 2.0_8
					divzy = (vz_grid(i,j+1,m) - vz_grid(i,j-1,m)) / 2.0_8
					divzz = (vz_grid(i,j,m+1) - vz_grid(i,j,m-1)) / 2.0_8
					if( divxx + divyy + divzz .ge. 0.0_8 ) then
						sigcx(i,j,m) = (divxx + divyy + divzz)**2
					else
						sigcy(i,j,m) = (divxx + divyy + divzz)**2
					end if
					sigsx(i,j,m) = (divzy - divyz)**2 + (divxz - divzx)**2 + (divyx - divxy)**2
				end do
			end do
		end do
		print *, 'done'

	end subroutine local_decompose
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine fft_smooth_1D(Gwidth, grid, nx, smooth, window)
	use MKL_DFTI
! Convolves input 'grid' with window function of FWHM 'Gwidth', given in number of cells
! 'grid' is assumed to be 1 dimensional, and the convolution is done in 1 dimension
! 'nx' is the number of cells in 1D in the grid and MUST be an even number
! If 'smooth'=1, then the input grid is replaced with the smoothed (convolved) grid
! If 'smooth'=0, then the smoothed (convolved) grid is subtracted from the input grid
! If 'window'=1 then a top hat window function is used
! If 'window'=2 then a Gaussian window function is used

		implicit none
		real(4),intent(in) :: Gwidth
		integer,intent(in) :: nx, smooth, window
		real(8),intent(inout) :: grid(nx)
		complex(8),allocatable :: fft_in1d(:), gauss_1d(:)
		real(8) :: Hwidth, Hwidth2, rx
		integer :: i, error, gw, fw, thx(2), nxhalf
		type(DFTI_DESCRIPTOR), POINTER :: fft_desc1

		if( smooth .ne. 0 .and. smooth .ne. 1 ) then
			print *, 'fft_smooth_1D, bad input for smooth parameter'
			print *, 'It should be 1 to replace original grid with smooth grid'
			print *, 'It should be 0 to subtract smoothed grid from original grid'
			print *, 'But you input',smooth
			stop
		end if
		if( window .ne. 1 .and. window .ne. 2 ) then
			print *, 'fft_smooth_1D, bad input for window parameter'
			print *, 'It should be 1 to use top hat window function'
			print *, 'It should be 2 to use Gaussian window function'
			print *, 'But you input',window
			stop
		end if

		allocate( fft_in1d(nx), stat=i )                      !!! deallocated at the end of this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation fft_in1d. stat= ', i
			stop
		end if
		allocate( gauss_1d(nx), stat=i )                      !!! deallocated at the end of this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation gauss_1d. stat= ', i
			stop
		end if

		rx = sngl(nx+1)/2.0_8
		fft_in1d(:) = grid(:)
		if( window .eq. 2 ) then
			Hwidth2 = (real(Gwidth,8)/2.0_8)**2
			do i=1,nx
				gauss_1d(i) = 0.5_8**( ((real(i)-rx)**2) / Hwidth2 ) 
			end do
		elseif( window .eq. 1 ) then
			Hwidth = real(Gwidth,8)/2.0_8
			gauss_1d(:) = 0.0_8
			thx(1) = ceiling( rx - Hwidth )
			thx(2) = ceiling( rx + Hwidth ) - 1
			do i=thx(1),thx(2)
				gauss_1d(i) = 1.0_8
			end do
		end if
		gauss_1d(:) = gauss_1d(:) / sum(gauss_1d(:))

		error = DftiCreateDescriptor(fft_desc1, DFTI_DOUBLE, DFTI_COMPLEX, 1, nx )
		error = DftiCommitDescriptor(fft_desc1)
		error = DftiComputeForward(fft_desc1, fft_in1d)
		error = DftiFreeDescriptor(fft_desc1)

		error = DftiCreateDescriptor(fft_desc1, DFTI_DOUBLE, DFTI_COMPLEX, 1, nx )
		error = DftiCommitDescriptor(fft_desc1)
		error = DftiComputeForward(fft_desc1, gauss_1d)
		error = DftiFreeDescriptor(fft_desc1)

		fft_in1d(:) = fft_in1d(:)*gauss_1d(:)
		print *, 'convolved'
		deallocate(gauss_1d)

		error = DftiCreateDescriptor(fft_desc1, DFTI_DOUBLE, DFTI_COMPLEX, 1, nx )
		error = DftiCommitDescriptor(fft_desc1)
		error = DftiComputeBackward(fft_desc1, fft_in1d)
		error = DftiFreeDescriptor(fft_desc1)
		print *, 'max and min, real and imag parts of inverse FFT'
		print '(4(1x,es12.5))', maxval(real(fft_in1d)), minval(real(fft_in1d)), maxval(aimag(fft_in1d)), minval(aimag(fft_in1d))

		print *, 'making output smoothed matrix'
		gw = 1-smooth
		fw = 1-2*smooth
		fft_in1d = fft_in1d / nx
		nxhalf = nx/2
		do i=1,nx/2
			grid(i)        = gw*grid(i)        - fw*real(fft_in1d(nxhalf+i))
			grid(nxhalf+i) = gw*grid(nxhalf+i) - fw*real(fft_in1d(i))
		end do
		deallocate(fft_in1d)
		print *, 'done with 1D smoothing'
	end subroutine fft_smooth_1D
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine fft_smooth_2D(Gwidth, grid, nx, ny, smooth, window)
	use MKL_DFTI
! Convolves input 'grid' with window function of FWHM 'Gwidth', given in number of cells
! 'grid' is assumed to be 2 dimensional, and the convolution is done in 2 dimensions
! 'nx' and 'ny' are the number of cells in 1D in the grid, and MUST both be even
! If 'smooth'=1, then the input grid is replaced with the smoothed (convolved) grid
! If 'smooth'=0, then the smoothed (convolved) grid is subtracted from the input grid
! If 'window'=1 then a top hat window function is used
! If 'window'=2 then a Gaussian window function is used

		implicit none
		real(4),intent(in) :: Gwidth
		integer,intent(in) :: nx, ny, smooth, window
		real(8),intent(inout) :: grid(nx, ny)
		complex(8),allocatable :: fft_in1d(:), gauss_1d(:)
		real(8) :: Hwidth, Hwidth2, rx, ry
		integer :: i, j, error, length(2), nxy, gw, fw, thx(2), thy(2), nxhalf, nyhalf
		type(DFTI_DESCRIPTOR), POINTER :: fft_desc2

		if( smooth .ne. 0 .and. smooth .ne. 1 ) then
			print *, 'fft_smooth_2D, bad input for smooth parameter'
			print *, 'It should be 1 to replace original grid with smooth grid'
			print *, 'It should be 0 to subtract smoothed grid from original grid'
			print *, 'But you input',smooth
			stop
		end if
		if( window .ne. 1 .and. window .ne. 2 ) then
			print *, 'fft_smooth_2D, bad input for window parameter'
			print *, 'It should be 1 to use top hat window function'
			print *, 'It should be 2 to use Gaussian window function'
			print *, 'But you input',window
			stop
		end if
		length = (/ nx, ny /)
		nxy = nx*ny
		rx = real(nx+1)/2.0_8
		ry = real(ny+1)/2.0_8

		allocate( fft_in1d(nxy), stat=i )                      !!! deallocated at the end of this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation fft_in1d. stat= ', i
			stop
		end if
		allocate( gauss_1d(nxy), stat=i )                      !!! deallocated at the end of this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation gauss_1d. stat= ', i
			stop
		end if

		do j=1,ny
			do i=1,nx
				fft_in1d( (j-1)*nx + i ) = grid(i,j)
			end do
		end do
		if( window .eq. 2 ) then
			Hwidth2 = (real(Gwidth,8)/2.0_8)**2
			do j=1,ny
				do i=1,nx
					gauss_1d( (j-1)*nx + i ) = 0.5_8**( ((real(i)-rx)**2 + (real(j)-ry)**2) / Hwidth2 ) 
				end do
			end do
		elseif( window .eq. 1 ) then
			Hwidth = real(Gwidth,8)/2.0_8
			gauss_1d(:) = 0.0_8
			thx(1) = ceiling( rx - Hwidth )
			thy(1) = ceiling( ry - Hwidth )
			thx(2) = ceiling( rx + Hwidth ) - 1
			thy(2) = ceiling( ry + Hwidth ) - 1
			do j=thy(1),thy(2)
				do i=thx(1),thx(2)
					gauss_1d( (j-1)*nx + i ) = 1.0_8
				end do
			end do
		end if
		gauss_1d(:) = gauss_1d(:) / sum(gauss_1d(:))

		error = DftiCreateDescriptor(fft_desc2, DFTI_DOUBLE, DFTI_COMPLEX, 2, length )
		error = DftiCommitDescriptor(fft_desc2)
		error = DftiComputeForward(fft_desc2, fft_in1d)
		error = DftiFreeDescriptor(fft_desc2)

		error = DftiCreateDescriptor(fft_desc2, DFTI_DOUBLE, DFTI_COMPLEX, 2, length )
		error = DftiCommitDescriptor(fft_desc2)
		error = DftiComputeForward(fft_desc2, gauss_1d)
		error = DftiFreeDescriptor(fft_desc2)

		fft_in1d(:) = fft_in1d(:)*gauss_1d(:)
		print *, 'convolved'
		deallocate(gauss_1d)

		error = DftiCreateDescriptor(fft_desc2, DFTI_DOUBLE, DFTI_COMPLEX, 2, length )
		error = DftiCommitDescriptor(fft_desc2)
		error = DftiComputeBackward(fft_desc2, fft_in1d)
		error = DftiFreeDescriptor(fft_desc2)
		print *, 'max and min, real and imag parts of inverse FFT'
		print '(4(1x,es12.5))', maxval(real(fft_in1d)), minval(real(fft_in1d)), maxval(aimag(fft_in1d)), minval(aimag(fft_in1d))

		print *, 'making output smoothed matrix'
		gw = 1-smooth
		fw = 1-2*smooth
		fft_in1d = fft_in1d / nxy
		nxhalf = nx/2
		nyhalf = ny/2
		do j=1,nyhalf
			do i=1,nxhalf
				grid(i, j)               = gw*grid(i, j)               - fw*real(fft_in1d( (nyhalf+j-1)*nx + nxhalf+i ))
				grid(nxhalf+i, nyhalf+j) = gw*grid(nxhalf+i, nyhalf+j) - fw*real(fft_in1d( (j-1)*nx + i ))
				grid(i, nyhalf+j)        = gw*grid(i, nyhalf+j)        - fw*real(fft_in1d( (j-1)*nx + nxhalf+i ))
				grid(nxhalf+i, j)        = gw*grid(nxhalf+i, j)        - fw*real(fft_in1d( (nyhalf+j-1)*nx + i ))
			end do
		end do
		deallocate(fft_in1d)
		print *, 'done with 2D smoothing'
	end subroutine fft_smooth_2D
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine fft_smooth_3D(Gwidth, grid, nx, ny, nz, smooth, window)
	use MKL_DFTI
! Convolves input 'grid' with Gaussian of FWHM 'Gwidth', given in number of cells
! 'grid' is assumed to be 3 dimensional, and the convolution is done in 3 dimensions
! 'nx' 'ny' and 'nz' are the number of cells in 1D in the grid and MUST all be even
! If 'smooth'=1, then the input grid is replaced with the smoothed (convolved) grid
! If 'smooth'=0, then the smoothed (convolved) grid is subtracted from the input grid
! If 'window'=1 then a top hat window function is used
! If 'window'=2 then a Gaussian window function is used

		implicit none
		real(4),intent(in) :: Gwidth
		integer,intent(in) :: nx, ny, nz, smooth, window
		real(8),intent(inout) :: grid(nx, ny, nz)
		complex(8),allocatable :: fft_in1d(:), gauss_1d(:)
		real(8) :: Hwidth, Hwidth2, rx, ry, rz
		integer :: i, j, k, error, length(3), nxy, nxyz, gw, fw, thx(2), thy(2), thz(2), nxhalf, nyhalf, nzhalf
		type(DFTI_DESCRIPTOR), POINTER :: fft_desc3

		print *, 'entering fft'
		if( smooth .ne. 0 .and. smooth .ne. 1 ) then
			print *, 'fft_smooth_3D, bad input for smooth parameter'
			print *, 'It should be 1 to replace original grid with smooth grid'
			print *, 'It should be 0 to subtract smoothed grid from original grid'
			print *, 'But you input',smooth
			stop
		end if
		if( window .ne. 1 .and. window .ne. 2 ) then
			print *, 'fft_smooth_3D, bad input for window parameter'
			print *, 'It should be 1 to use top hat window function'
			print *, 'It should be 2 to use Gaussian window function'
			print *, 'But you input',window
			stop
		end if
		length = (/ nx,ny,nz /)
		nxy = nx*ny
		nxyz = nxy*nz
		rx = real(nx+1)/2.0_8
		ry = real(ny+1)/2.0_8
		rz = real(nz+1)/2.0_8

		allocate( fft_in1d(nxyz), stat=i )				!!! deallocated in this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation fft_in1d. stat= ', i
			stop
		end if
		allocate( gauss_1d(nxyz), stat=i )				!!! deallocated in this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation gauss_1d. stat= ', i
			stop
		end if

		do k=1,nz
			do j=1,ny
				do i=1,nx
					fft_in1d( (k-1)*nxy + (j-1)*nx + i ) = grid(i,j,k)
				end do
			end do
		end do
		if( window .eq. 2 ) then
			Hwidth2  = (real(Gwidth,8)/2.0_8)**2
			do k=1,nz
				do j=1,ny
					do i=1,nx
						gauss_1d( (k-1)*nxy + (j-1)*nx + i ) = 0.5_4**( ((real(i)-rx)**2 + (real(j)-ry)**2 + (real(k)-rz)**2) / Hwidth2 )
					end do
				end do
			end do
		elseif( window .eq. 1 ) then
			Hwidth = real(Gwidth,8)/2.0_8
			gauss_1d(:) = 0.0_8
			thx(1) = ceiling( rx - Hwidth )
			thy(1) = ceiling( ry - Hwidth )
			thz(1) = ceiling( rz - Hwidth )
			thx(2) = ceiling( rx + Hwidth ) - 1
			thy(2) = ceiling( ry + Hwidth ) - 1
			thz(2) = ceiling( rz + Hwidth ) - 1
			do k=thz(1),thz(2)
				do j=thy(1),thy(2)
					do i=thx(1),thx(2)
						gauss_1d( (k-1)*nxy + (j-1)*nx + i ) = 1.0_8
					end do
				end do
			end do
		end if
		gauss_1d(:) = gauss_1d(:) / sum(gauss_1d(:))

		error = DftiCreateDescriptor(fft_desc3, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc3)
		error = DftiComputeForward(fft_desc3, fft_in1d)
		error = DftiFreeDescriptor(fft_desc3)

		error = DftiCreateDescriptor(fft_desc3, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc3)
		error = DftiComputeForward(fft_desc3, gauss_1d)
		error = DftiFreeDescriptor(fft_desc3)

		fft_in1d(:) = fft_in1d(:)*gauss_1d(:)
		print *, 'convolved'
		deallocate(gauss_1d)

		error = DftiCreateDescriptor(fft_desc3, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc3)
		error = DftiComputeBackward(fft_desc3, fft_in1d)
		error = DftiFreeDescriptor(fft_desc3)
		print *, 'max and min, real and imag parts of inverse FFT'
		print '(4(1x,es12.5))', maxval(real(fft_in1d)), minval(real(fft_in1d)), maxval(aimag(fft_in1d)), minval(aimag(fft_in1d))

		print *, 'making output smoothed matrix'
		gw = 1-smooth
		fw = 1-2*smooth
		fft_in1d = fft_in1d / nxyz
		print *, 'maximum and minimum real parts of inverse FFT'
		print '(2(1x,es12.5))', maxval(real(fft_in1d)), minval(real(fft_in1d))
		print *, 'maximum and minimum imaginary parts of inverse FFT'
		print '(2(1x,es12.5))', maxval(aimag(fft_in1d)), minval(aimag(fft_in1d))
		nxhalf = nx/2
		nyhalf = ny/2
		nzhalf = nz/2
		do k=1,nzhalf
			do j=1,nyhalf
				do i=1,nxhalf
					grid(i, j, k)                      = gw*grid(i, j, k)                      - fw*real(fft_in1d( (nzhalf+k-1)*nxy + (nyhalf+j-1)*nx + nxhalf+i ))
					grid(nxhalf+i, nyhalf+j, nzhalf+k) = gw*grid(nxhalf+i, nyhalf+j, nzhalf+k) - fw*real(fft_in1d( (k-1)*nxy + (j-1)*nx + i ))
					grid(i, j, nzhalf+k)               = gw*grid(i, j, nzhalf+k)               - fw*real(fft_in1d( (k-1)*nxy + (nyhalf+j-1)*nx + nxhalf+i ))
					grid(nxhalf+i, nyhalf+j, k)        = gw*grid(nxhalf+i, nyhalf+j, k)        - fw*real(fft_in1d( (nzhalf+k-1)*nxy + (j-1)*nx + i ))
					grid(i, nyhalf+j, k)               = gw*grid(i, nyhalf+j, k)               - fw*real(fft_in1d( (nzhalf+k-1)*nxy + (j-1)*nx + nxhalf+i ))
					grid(nxhalf+i, j, nzhalf+k)        = gw*grid(nxhalf+i, j, nzhalf+k)        - fw*real(fft_in1d( (k-1)*nxy + (nyhalf+j-1)*nx + i ))
					grid(nxhalf+i, j, k)               = gw*grid(nxhalf+i, j, k)               - fw*real(fft_in1d( (nzhalf+k-1)*nxy + (nyhalf+j-1)*nx + i ))
					grid(i, nyhalf+j, nzhalf+k)        = gw*grid(i, nyhalf+j, nzhalf+k)        - fw*real(fft_in1d( (k-1)*nxy + (j-1)*nx + nxhalf+i ))
				end do
			end do
		end do
		deallocate(fft_in1d)
		print *, 'done with 3D smoothing'
	end subroutine fft_smooth_3D
end module decomp
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program main
!#include "adef.h"
use parameters
use globvar
use read_binary
use splitter
use decomp
!use omp_lib

	implicit none
	character(len=20),allocatable :: gal_name(:)	!!! This is the name of the galaxy, e.g. MW3, VL01, SFG1, MW10, VELA07, VELA_v2_10, etc.
	character(len=4) :: Rvir_string
	character(len=20) :: snap_tag
	character(len=7) :: comm
	character(len=256) :: filename, path_name, input_arg
	character(len=256),allocatable :: gas_file_name(:)
	integer :: i, j, k, l, m, Nsnapshot, Nsnapshot2, Nsimulation
	real(4) :: aexp2, Rv4, buff_size

	call read_param()

!!!!!!!!!! What simulations will we be looking at today? !!!!!!!!!!
	if(iargc().ge.1) then
		call getarg(1,input_arg)
		write(filename,'(a)') trim(input_arg)
	else
		write(filename,'(a)') './input_output/turbulence_decomp_input.dat'
	end if
	print *, 'input filename'
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
		print *, ''
		print *, 'sim',j,'of',Nsimulation
		print *, gal_name(j)
		Nsnapshot = 0
		Nsnapshot2 = 0
		!!!!!!!!!! Get array sizes !!!!!!!!!!
		write(filename,'(a,a,a)') '/BIGDATA/nirm/galaxy_inputs/',trim(gal_name(j)),'/Nmax.txt'
		open(unit=16,file=filename,form='formatted')
		read(16,'(1x,i3)') Nsnapshot
		print *, 'has ',Nsnapshot,' snapshots'
		allocate( aexp(Nsnapshot), Ngas(Nsnapshot), Ndm(Nsnapshot), Nstars(Nsnapshot) )          !!! deallocated at the end of simulation loop !!!
		do i=1,Nsnapshot
			read(16,'(1x,f5.3,3(1x,i))') aexp(i), Ngas(i), Nstars(i), Ndm(i)
		end do
		close(unit=16)
		print *, 'first and last aexp, Ngas, Nstars, Ndm'
		print *, aexp(1), Ngas(1), Nstars(1), Ndm(1)
		print *, aexp(Nsnapshot), Ngas(Nsnapshot), Nstars(Nsnapshot), Ndm(Nsnapshot)

		call allocate_global(Nsnapshot)
		!!!!!!!!!! Initialize the arrays per simulation !!!!!!!!!!
		write(filename,'(a,a,a)') '/BIGDATA/nirm/galaxy_inputs/',trim(gal_name(j)),'/galaxy_catalogue/Nir_halo_cat.txt'
		open(unit=17,file=filename,form='formatted')
		read(17,'(1x,i3)') Nsnapshot2
		if( Nsnapshot2 .ne. Nsnapshot ) then
			print *, 'DATA Nsnapshot - error with halo cat or Nmax. Inconsistent number of snapshots per galaxy'
			print *, 'Nsnapshot',Nsnapshot,'Nsnapshot2',Nsnapshot2
			stop
		end if
		do i=1,Nsnapshot
			read(17,'(1x,f5.3,3(1x,es12.5))') aexp2, Rvir(i), Mvir(i), Vvir(i)
			if( aexp2 .ne. aexp(i) ) then
				print *, 'DATA aexp - error with halo cat or Nmax. Inconsistent expansion factor'
				print *, 'snapshot #',i,'aexp',aexp(i),'aexp2',aexp2
				stop
			end if
		end do
		close(unit=17)

		write(filename,'(a,a,a)') '/BIGDATA/nirm/galaxy_inputs/',trim(gal_name(j)),'/galaxy_catalogue/Nir_disc_cat.txt'
		open(unit=18,file=filename,form='formatted')
		read(18,'(1x,i3)') Nsnapshot2
		if( Nsnapshot2 .ne. Nsnapshot ) then
			print *, 'DATA Nsnapshot - error with disc cat. Inconsistent number of snapshots per galaxy'
			print *, 'Nsnapshot',Nsnapshot,'Nsnapshot2',Nsnapshot2
			stop
		end if
		do i=1,Nsnapshot
			read(18,'(1x,f5.3,10(1x,e12.5),2(1x,f7.3),7(1x,es12.5),2(1x,f7.3))') aexp2, rcom(1,i), rcom(2,i), rcom(3,i), &
			& vcom(1,i), vcom(2,i), vcom(3,i), Ldisc(1,i), Ldisc(2,i), Ldisc(3,i), Lmag(i), Rdisc(i), Hdisc(i), &
			& Mgas_disc(i), Mcold_disc(i), Mstar_disc(i), M_Es_star_disc(i), Mdm_disc(i), SFR_disc(i), age_disc(i), metgas_disc(i), metstars_disc(i)
			if( aexp2 .ne. aexp(i) ) then
				print *, 'DATA aexp - error with disc cat. Inconsistent expansion factor'
				print *, 'snapshot #',i,'aexp',aexp(i),'aexp2',aexp2
				stop
			end if
		end do
		close(unit=18)

		call open_files()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Loop over all snapshots !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Allocate gas data and compute grid size !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			call allocate_gas( k )
			if(use_fixed_R) then
				grid_rad = Fixed_R					!!! SHOULD BE IN KPC
			elseif(use_Rvir) then
				grid_rad = R_over_Rvir * Rvir(k)			!!! SHOULD BE IN KPC
			elseif(use_old_R) then
				grid_rad = Rdisc(k)					!!! SHOULD BE IN KPC
			end if
			if(use_fixed_H) then
				grid_height = Fixed_H					!!! SHOULD BE IN KPC
			elseif(use_R) then
				grid_height = H_over_R * grid_rad			!!! SHOULD BE IN KPC
			elseif(use_old_H) then
				grid_height = Hdisc(k)					!!! SHOULD BE IN KPC
			end if
			grid_rad = max(grid_rad, 2.0_4)					!!! SHOULD BE IN KPC
			grid_rad = min(grid_rad, Rvir(k))				!!! SHOULD BE IN KPC
			grid_height = max(grid_height, 1.0_4)				!!! SHOULD BE IN KPC
			grid_height = min(grid_height, grid_rad)			!!! SHOULD BE IN KPC
			FWHM_smooth = min( pre_smooth_FWHM, 500.0_4*grid_rad )		!!! 0.5 the disc radius, in pc
			print *, 'radius, height, FWHM'
			print *, grid_rad, grid_height, FWHM_smooth/1000.0_4

!			buff_size = 5.0_4*res
			buff_size = 0.0_4
			if( pre_smooth ) then
				buff_size = max(buff_size, 1.5_4 * FWHM_smooth)
			end if
			if( decomp_type .eq. 3 ) then
				buff_size = max(buff_size, 1.5_4 * type_3_FWHM)
			end if
			FWHM_smooth = FWHM_smooth / res

			write(filename,'(a,a,a,a)') trim(output_dirname),'/binary_grid_outputs/',trim(snap_tag2),'_raw_maps.bin'
			open(unit=31,file=filename,form='unformatted')
			write(filename,'(a,a,a,a)') trim(output_dirname),'/binary_grid_outputs/',trim(snap_tag2),'_turbulence.bin'
			open(unit=32,file=filename,form='unformatted')
			write(filename,'(a,a,a,a)') trim(output_dirname),'/binary_grid_outputs/',trim(snap_tag2),'_decomposed_turbulence_mass.bin'
			open(unit=34,file=filename,form='unformatted')
			write(filename,'(a,a,a,a)') trim(output_dirname),'/binary_grid_outputs/',trim(snap_tag2),'_decomposed_turbulence_volume.bin'
			open(unit=35,file=filename,form='unformatted')
			write(filename,'(a,a,a,a)') trim(output_dirname),'/binary_grid_outputs/',trim(snap_tag2),'_mode_ratios_mass.bin'
			open(unit=36,file=filename,form='unformatted')
			write(filename,'(a,a,a,a)') trim(output_dirname),'/binary_grid_outputs/',trim(snap_tag2),'_mode_ratios_volume.bin'
			open(unit=37,file=filename,form='unformatted')

			call grid_size( grid_rad * 1000.0_4, grid_height * 1000.0_4, buff_size, snap_tag2, k )	!!! Input dimensions in pc
			call create_grid( k )
			call rot_curve_surface_maps( k )
			call turb_def()
			call turb_decomp( k )
			call deallocate_all()

			close(unit=31)
			close(unit=32)
			close(unit=34)
			close(unit=35)
			close(unit=36)
			close(unit=37)
			print *, ''
		end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		print *, 'done with', gal_name(j)
		print *, Nsimulation-j,' simulations left'

		deallocate( aexp, Ngas, Nstars, Ndm )
		deallocate( gas_file_name )
		deallocate( rcom, vcom, Ldisc, Lmag, Rdisc, Hdisc )
		deallocate( Mgas_disc, Mcold_disc, Mstar_disc, M_Es_star_disc, Mdm_disc, SFR_disc, age_disc, metgas_disc, metstars_disc )
		deallocate( Rvir, Mvir, Vvir )
		call close_files()
		print *, ''
	end do
	deallocate( gal_name )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contains
	subroutine read_param()
!!!!!!!!!! Read in parameter file and check for inconsistencies !!!!!!!!!!
		implicit none
		character(len=24) :: s
		integer :: j

		write(filename,'(a)') './parameter_input.dat'
		open(unit=14,file=filename)
		read(14,'(a24,E9.2E2)') s,ngrid_max
		print '(a24,E9.2E2)', s,ngrid_max
		read(14,'(a24,I1)') s,j
		print '(a24,I1)', s,j
		if(j .eq. 1) then
			check_CIC = .true.
		elseif(j .eq. 0) then
			check_CIC = .false.
		else
			print *, 'Problem with check_CIC in parameter input file'
			print *, 'It should be equal to 0 or 1, but instead it is:'
			print *, j
			stop
		end if
		read(14,'(a24,I1)') s,j
		print '(a24,I1)', s,j
		if(j .eq. 1) then
			move_center = .true.
		elseif(j .eq. 0) then
			move_center = .false.
		else
			print *, 'Problem with move_center in parameter input file'
			print *, 'It should be equal to 0 or 1, but instead it is:'
			print *, j
			stop
		end if
		read(14,'(a24,I1)') s,j
		print '(a24,I1)', s,j
		if(j .eq. 1) then
			move_velocity_center = .true.
		elseif(j .eq. 0) then
			move_velocity_center = .false.
		else
			print *, 'Problem with move_velocity_center in parameter input file'
			print *, 'It should be equal to 0 or 1, but instead it is:'
			print *, j
			stop
		end if
		read(14,'(a24,F6.1)') s,res
		print '(a24,F6.1)', s,res
		read(14,'(a24,F4.1)') s,decomp_weight
		print '(a24,F4.1)', s,decomp_weight
		read(14,'(a24,I1)') s,j
		print '(a24,I1)', s,j
		if(j .eq. 1) then
			pre_smooth = .true.
		elseif(j .eq. 0) then
			pre_smooth = .false.
		else
			print *, 'Problem with pre_smooth in parameter input file'
			print *, 'It should be equal to 0 or 1, but instead it is:'
			print *, j
			stop
		end if
		read(14,'(a24,I1)') s,pre_smooth_dim
		print '(a24,I1)', s,pre_smooth_dim
		if(pre_smooth_dim .ne. 2 .and. pre_smooth_dim .ne. 3) then
			print *, 'Problem with pre_smooth_dim in parameter input file'
			print *, 'It should be equal to 2 or 3, but instead it is:'
			print *, pre_smooth_dim
			stop
		end if
		read(14,'(a24,F6.1)') s,pre_smooth_FWHM
		print '(a24,F6.1)', s,pre_smooth_FWHM
		read(14,'(a24,I1)') s,decomp_type
		print '(a24,I1)', s,decomp_type
		if( decomp_type .ne. 1 .and. decomp_type .ne. 2 .and. decomp_type .ne. 3 ) then
			print *, 'Problem with decomp_type in parameter input file'
			print *, 'It should be equal to 1 or 2 or 3, but instead it is:'
			print *, decomp_type
			stop
		end if
		read(14,'(a24,F6.1)') s,type_2_FWHM
		print '(a24,F6.1)', s,type_2_FWHM
		read(14,'(a24,I3)') s,Nvert
		print '(a24,I3)', s,Nvert
		read(14,'(a24,I1)') s,type_2_mass
		print '(a24,I1)', s,type_2_mass
		if( type_2_mass .ne. 0 .and. type_2_mass .ne. 1 ) then
			print *, 'Problem with type_2_mass in parameter input file'
			print *, 'It should be equal to 0 or 1, but instead it is:'
			print *, type_2_mass
			stop
		end if
		read(14,'(a24,F6.1)') s,type_3_FWHM
		print '(a24,F6.1)', s,type_3_FWHM
		read(14,'(a24,I1)') s,j
		print '(a24,I1)', s,j
		if(j .eq. 1) then
			use_fixed_R = .true.
		elseif(j .eq. 0) then
			use_fixed_R = .false.
		else
			print *, 'Problem with use_fixed_R in parameter input file'
			print *, 'It should be equal to 0 or 1, but instead it is:'
			print *, j
			stop
		end if
		read(14,'(a24,F4.1)') s,Fixed_R
		print '(a24,F4.1)', s,Fixed_R
		read(14,'(a24,I1)') s,j
		print '(a24,I1)', s,j
		if(j .eq. 1) then
			use_Rvir = .true.
		elseif(j .eq. 0) then
			use_Rvir = .false.
		else
			print *, 'Problem with use_Rvir in parameter input file'
			print *, 'It should be equal to 0 or 1, but instead it is:'
			print *, j
			stop
		end if
		read(14,'(a24,F4.2)') s,R_over_Rvir
		print '(a24,F4.2)', s,R_over_Rvir
		read(14,'(a24,I1)') s,j
		print '(a24,I1)', s,j
		if(j .eq. 1) then
			use_old_R = .true.
		elseif(j .eq. 0) then
			use_old_R = .false.
		else
			print *, 'Problem with use_old_R in parameter input file'
			print *, 'It should be equal to 0 or 1, but instead it is:'
			print *, j
			stop
		end if
		read(14,'(a24,I1)') s,j
		print '(a24,I1)', s,j
		if(j .eq. 1) then
			use_fixed_H = .true.
		elseif(j .eq. 0) then
			use_fixed_H = .false.
		else
			print *, 'Problem with use_fixed_H in parameter input file'
			print *, 'It should be equal to 0 or 1, but instead it is:'
			print *, j
			stop
		end if
		read(14,'(a24,F4.1)') s,Fixed_H
		print '(a24,F4.1)', s,Fixed_H
		read(14,'(a24,I1)') s,j
		print '(a24,I1)', s,j
		if(j .eq. 1) then
			use_R = .true.
		elseif(j .eq. 0) then
			use_R = .false.
		else
			print *, 'Problem with use_R in parameter input file'
			print *, 'It should be equal to 0 or 1, but instead it is:'
			print *, j
			stop
		end if
		read(14,'(a24,F4.2)') s,H_over_R
		print '(a24,F4.2)', s,H_over_R
		read(14,'(a24,I1)') s,j
		print '(a24,I1)', s,j
		if(j .eq. 1) then
			use_old_H = .true.
		elseif(j .eq. 0) then
			use_old_H = .false.
		else
			print *, 'Problem with use_old_H in parameter input file'
			print *, 'It should be equal to 0 or 1, but instead it is:'
			print *, j
			stop
		end if
		read(14,'(a24,I1)') s,j
		print '(a24,I1)', s,j
		if(j .eq. 1) then
			use_temp_thresh = .true.
		elseif(j .eq. 0) then
			use_temp_thresh = .false.
		else
			print *, 'Problem with use_temp_thresh in parameter input file'
			print *, 'It should be equal to 0 or 1, but instead it is:'
			print *, j
			stop
		end if
		read(14,'(a24,E9.2E2)') s,temp_thresh
		print '(a24,E9.2E2)', s,temp_thresh
		close(unit=14)
		if( use_fixed_R ) then
			if( use_old_R .or. use_Rvir ) then
				print *, 'Multiple different definitions for disc R'
				print *, 'MAKE UP YOUR MIND!!'
				print *, use_fixed_R, use_Rvir, use_old_R
				stop
			end if
		else
			if( use_Rvir ) then
				if( use_old_R ) then
					print *, 'Multiple different definitions for disc R'
					print *, 'MAKE UP YOUR MIND!!'
					print *, use_fixed_R, use_Rvir, use_old_R
					stop
				end if
			else
				if( .not. use_old_R ) then
					print *, 'No disc R definition was chosen'
					print *, 'MAKE UP YOUR MIND!!'
					print *, use_fixed_R, use_Rvir, use_old_R
					stop
				end if
			end if
		end if
		if( use_fixed_H ) then
			if( use_old_H .or. use_R ) then
				print *, 'Multiple different definitions for disc H'
				print *, 'MAKE UP YOUR MIND!!'
				print *, use_fixed_H, use_R, use_old_H
				stop
			end if
		else
			if( use_R ) then
				if( use_old_H ) then
					print *, 'Multiple different definitions for disc H'
					print *, 'MAKE UP YOUR MIND!!'
					print *, use_fixed_H, use_R, use_old_H
					stop
				end if
			else
				if( .not. use_old_H ) then
					print *, 'No disc H definition was chosen'
					print *, 'MAKE UP YOUR MIND!!'
					print *, use_fixed_H, use_R, use_old_H
					stop
				end if
			end if
		end if
	end subroutine read_param
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine allocate_global(N)
!!!!!!!!!! Allocate global arrays per simulation !!!!!!!!!!
		implicit none
		integer,intent(in) :: N
		integer :: i

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
!_____________________________________________________________________________________________________________________________________________________________________
	subroutine allocate_gas(snap)
!!!!!!!!!! Allocate global arrays per simulation !!!!!!!!!!
		implicit none
		integer,intent(in) :: snap
		integer :: i, nold

		allocate( xgas(Ngas(snap)), stat=i )      	!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation xgas. stat= ', i
			stop
		end if
		allocate( ygas(Ngas(snap)), stat=i )      	!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation ygas. stat= ', i
			stop
		end if
		allocate( zgas(Ngas(snap)), stat=i )      	!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation zgas. stat= ', i
			stop
		end if
		allocate( vxgas(Ngas(snap)), stat=i )     	!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation vxgas. stat= ', i
			stop
		end if
		allocate( vygas(Ngas(snap)), stat=i )     	!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation vygas. stat= ', i
			stop
		end if
		allocate( vzgas(Ngas(snap)), stat=i )     	!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation vzgas. stat= ', i
			stop
		end if
		allocate( density_gas(Ngas(snap)), stat=i )	!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation density_gas. stat= ', i
			stop
		end if
		allocate( cell_size_gas(Ngas(snap)), stat=i )	!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation cell_size_gas. stat= ', i
			stop
		end if
		allocate( temperature_gas(Ngas(snap)), stat=i )	!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation temperature_gas. stat= ', i
			stop
		end if

		call data_gas( gas_file_name(snap), Ngas(snap) )

		!!! Move gas to disc rest frame, then don't have to worry about it later !!!
		if(move_center) then
			xgas(:) = xgas(:) - rcom(1,snap)
			ygas(:) = ygas(:) - rcom(2,snap)
			zgas(:) = zgas(:) - rcom(3,snap)
		end if
		if(move_velocity_center) then
			vxgas(:) = vxgas(:) - vcom(1,snap)
			vygas(:) = vygas(:) - vcom(2,snap)
			vzgas(:) = vzgas(:) - vcom(3,snap)
		end if

	end subroutine allocate_gas
!_____________________________________________________________________________________________________________________________________________________________________
	subroutine deallocate_all()
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

	end subroutine deallocate_all
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine open_files()
		implicit none

		write(output_dirname,'(a)') './outputs'
		write(filename,'(a,a)') 'mkdir -p ',trim(output_dirname)
		call system(trim(filename))

		write(output_dirname,'(a,a,i1)') trim(output_dirname),'/V',decomp_type
		if( decomp_type .eq. 2 ) then
			write(output_dirname,'(a,a,i4.4,a,i3.3)') trim(output_dirname),'_',floor(type_2_FWHM),'pc_Nvert_',Nvert
			if( type_2_mass .eq. 1 ) then
				write(output_dirname,'(a,a)') trim(output_dirname),'_M'
			elseif( type_2_mass .eq. 0 ) then
				write(output_dirname,'(a,a)') trim(output_dirname),'_V'
			end if
		elseif( decomp_type .eq. 3 ) then
			write(output_dirname,'(a,a,i4.4,a)') trim(output_dirname),'_',floor(type_3_FWHM),'pc'
		end if
		write(filename,'(a,a)') 'mkdir -p ',trim(output_dirname)
		call system(trim(filename))

		write(output_dirname,'(a,a,f3.1)') trim(output_dirname),'/weight_',decomp_weight
		write(filename,'(a,a)') 'mkdir -p ',trim(output_dirname)
		call system(trim(filename))

		write(output_dirname,'(a,a,i3.3,a)') trim(output_dirname),'/res_',floor(res),'pc'
		write(filename,'(a,a)') 'mkdir -p ',trim(output_dirname)
		call system(trim(filename))

		if( pre_smooth ) then
			write(output_dirname,'(a,a,i1,a,i4.4,a)') trim(output_dirname),'/',pre_smooth_dim,'D_smooth_',floor(pre_smooth_FWHM),'pc'
		else
			write(output_dirname,'(a,a)') trim(output_dirname),'/no_pre_smooth'
		end if
		write(filename,'(a,a)') 'mkdir -p ',trim(output_dirname)
		call system(trim(filename))

		if(.not. use_temp_thresh) then
			temp_thresh = 1.E9
		end if

		if(use_fixed_R) then
			if( Fixed_R .ge. 10.0_4 ) then
				write(output_dirname,'(a,a,f4.1,a)') trim(output_dirname),'/',Fixed_R,'kpc_x_'
			else
				write(output_dirname,'(a,a,f3.1,a)') trim(output_dirname),'/',Fixed_R,'kpc_x_'
			end if
		elseif(use_old_R) then
			write(output_dirname,'(a,a)') trim(output_dirname),'/Rd_x_'
		elseif(use_Rvir) then
			write(output_dirname,'(a,a,f4.2,a)') trim(output_dirname),'/',R_over_Rvir,'Rv_x_'
		end if
		if(use_fixed_H) then
			write(output_dirname,'(a,f4.2,a,f3.1)') trim(output_dirname),Fixed_H,'kpc_Tmax_',log10(temp_thresh)
		elseif(use_old_H) then
			write(output_dirname,'(a,a,f3.1)') trim(output_dirname),'Hd_Tmax_',log10(temp_thresh)
		elseif(use_R) then
			write(output_dirname,'(a,f4.2,a,f3.1)') trim(output_dirname),H_over_R,'R_Tmax_',log10(temp_thresh)
		end if
		write(filename,'(a,a)') 'mkdir -p ',trim(output_dirname)
		call system(trim(filename))

		write(output_dirname,'(a,a,a)') trim(output_dirname),'/',trim(gal_name(j))
		write(filename,'(a,a)') 'mkdir -p ',trim(output_dirname)
		call system(trim(filename))

		write(filename,'(a,a,a)') 'mkdir -p ',trim(output_dirname),'/binary_grid_outputs'
		call system(trim(filename))

		write(filename,'(a,a,a)') 'mkdir -p ',trim(output_dirname),'/rotation_curves'
		call system(trim(filename))

		write(filename,'(a,a,a)') 'mkdir -p ',trim(output_dirname),'/power_spectra'
		call system(trim(filename))

		write(filename,'(a,a)') trim(output_dirname),'/too_large_grids.out'
		open(unit=19,file=filename)

		write(filename,'(a,a)') trim(output_dirname),'/compressive_over_solenoidal.out'
		open(unit=20,file=filename)

	end subroutine open_files
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine close_files()
		implicit none

		close(unit=19)
		close(unit=20)

	end subroutine close_files
end program main

