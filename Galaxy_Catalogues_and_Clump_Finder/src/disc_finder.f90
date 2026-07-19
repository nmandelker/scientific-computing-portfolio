module parameters
	implicit none
	integer,parameter :: Nmax = 55000000
	real(8),parameter :: fj = 0.7_8					!!! j_z fraction for disc stars

	real(8),parameter :: fM = 0.85_8				!!! Mass fraction for disc dimensions
	real(4),parameter :: fRvir = 0.15_4				!!! fraction of Rvir within which everything is initially calculated
	real(4),parameter :: initial_height = 1.0_4			!!! Initial height of disc for computing Rd, in [kpc]
	real(4),parameter :: fRd = 1.0_4				!!! fraction of Rd within which disc thickness is calculated

	real(4),parameter :: max_T = 1.5e4				!!! maximum temperature for "cold" gas in [K]
	real(4),parameter :: max_age = 0.1_4				!!! maximum age for "young" stars in [Gyr]
	logical,parameter :: use_stars = .true.				!!! disc defined with cold gas AND young stars
	logical,parameter :: move_center = .false.			!!! move to rest frame of tracer at each iteration

	real(4),parameter :: dtmax = 0.2_4				!!! For SFR calculation, Myr
	real(4),parameter :: tmax_i(3) = (/ 40.0_4, 30.0_4, 80.0_4 /)	!!! For SFR calculation, Myr
	real(4),parameter :: tmax_f(3) = (/ 80.0_4, 60.0_4, 120.0_4 /)	!!! For SFR calculation, Myr

	real(8),parameter :: A0 = 0.05_8				!!! For stellar mass loss
	real(8),parameter :: T0 = 0.005_8				!!! For stellar mass loss, in Gyr

	real(8),parameter :: G=4.3d-6					!!! Kpc*(km/sec)^2/M_sun
	real(4),parameter :: pi4_3=4.1887902047864_4, pi2=1.5_4*pi4_3, pi=0.5_4*pi2

!!!	Arrays for gas, stars and DM data
!!!	---------------------------------
	real(4),allocatable :: xgas(:), ygas(:), zgas(:), vxgas(:), vygas(:), vzgas(:)
	real(4),allocatable :: density_gas(:), temperature_gas(:), cell_size_gas(:), SNI_gas(:), SNII_gas(:)
	real(4) :: res
	integer :: Ngas

	real(8),allocatable :: xstars(:), ystars(:), zstars(:), vxstars(:), vystars(:), vzstars(:), mass_stars(:), initial_mass_stars(:)
	real(4),allocatable :: age_stars(:),SNI_stars(:),SNII_stars(:)
	integer,allocatable :: idstars(:), insitustars(:,:)
	integer :: Nstars

	real(8),allocatable :: xdm(:), ydm(:), zdm(:), vxdm(:), vydm(:), vzdm(:), mass_dm(:)
	integer,allocatable :: iddm(:)
	integer :: Ndm

!!!	Disc data
!!!	---------------------------------
	real(8) :: mean_mgas, mean_mcold, mean_mstars, mean_mdm, mean_mbar, mean_fgas, mean_sigma_gas, mean_sigma_stars, mean_sigma
	real(8) :: mean_SFR, mean_Sig_SFR, mean_SSFR, mean_tau, mean_age, mean_metgas, mean_metstars
	real(8) :: spher_mgas, spher_mcold, spher_mstars, spher_initial_mstars, spher_mdm, spher_mbar, spher_fgas, spher_sigma_gas, spher_sigma_stars,  spher_sigma
	real(8) :: spher_SFR, spher_Sig_SFR, spher_SSFR, spher_tau, spher_age, spher_metgas, spher_metstars
	real(8) :: exsitu_stellar_mass, spher_exsitu_stellar_mass
	real(8) :: Mgas_015Rv, Mstars_015Rv, initial_Mstars_015Rv, Mcold_gas_015Rv, Myoung_stars_015Rv, Mdm_015Rv, SFR_015Rv
	real(8) :: Mgas_1kpc, Mstars_1kpc, initial_Mstars_1kpc, Mcold_gas_1kpc, Myoung_stars_1kpc, Mdm_1kpc, SFR_1kpc
	real(8) :: initial_stellar_mass_disc, initial_stellar_mass_sphere, analytic_stellar_mass_disc, analytic_stellar_mass_sphere
	real(8) :: initial_stellar_mass_disc_es, initial_stellar_mass_sphere_es, analytic_stellar_mass_disc_es, analytic_stellar_mass_sphere_es
	real(8) :: SFR_test_disc(3), SFR_test_sphere(3)
	real(8) :: Mvir, Vvir
	real(4) :: Rvir
	real(8) :: Ms10, Ms0_1, MsV, IMs10, IMs0_1, IMsV, Mg10, Mg0_1, MgV, Md10, Md0_1, MdV

end module parameters
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module read_binary
use parameters
	implicit none
contains
	subroutine data_gas(filename)
		implicit none
		character (*),intent (in) :: filename
	
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
		Ngas=1

		DO WHILE (Ngas.lt.Nmax)
			read (12,end=6) cell_size_gas(Ngas),xgas(Ngas),ygas(Ngas),zgas(Ngas),Vxgas(Ngas),Vygas(Ngas),Vzgas(Ngas),&
			& density_gas(Ngas),temperature_gas(Ngas),SNII_gas(Ngas),SNI_gas(Ngas)

			Ngas=Ngas+1
		end do
 6   		continue
		close (12)
		if( Ngas.eq.Nmax ) then
			print *, 'DATA GAS - error reading data: too many cells. Change Nmax.'
			print *, 'Ngas=',Ngas
			stop
		end if
		Ngas=Ngas-1	
		print *,Ngas
	end subroutine data_gas
!________________________________________________________________________________________
	subroutine data_dm(filename)
		implicit none
		character (*),intent (in) :: filename
				
		open ( 13 , file = filename, form = 'unformatted' )
		iddm(:) = 0
		xdm(:) = 1.d-10
		ydm(:) = 1.d-10
		zdm(:) = 1.d-10
		vxdm(:) = 1.d-10
		vydm(:) = 1.d-10
		vzdm(:) = 1.d-10
		mass_dm(:) = 1.d-10
		Ndm=1
		DO WHILE (Ndm.lt.Nmax)
			read (13,end=6) iddm(Ndm),xdm(Ndm),ydm(Ndm),zdm(Ndm),Vxdm(Ndm),Vydm(Ndm),Vzdm(Ndm),mass_dm(Ndm)
			Ndm=Ndm+1
		end do
 6   		continue
		close (13)		
		if( Ndm.eq.Nmax ) then
			print *, 'DATA DM - error reading data: too many DM particles. Change Nmax.'
			print *, 'Ndm=',Ndm
			stop
		end if
		Ndm = Ndm-1
		print *,Ndm
	end subroutine data_dm
!________________________________________________________________________________________

	subroutine data_stars(filename, initial)
		implicit none
		character (*),intent (in) :: filename
		integer,intent(in) :: initial
				
		open ( 14 , file = filename, form = 'unformatted' )
		Nstars = 1
		if(initial .eq. 0) then
			idstars(:) = 0
			xstars(:) = 0.0_8
			ystars(:) = 0.0_8
			zstars(:) = 0.0_8
			vxstars(:) = 0.0_8
			vystars(:) = 0.0_8
			vzstars(:) = 0.0_8
			mass_stars(:) = 0.0_8
			age_stars(:) = 0.0_4
			DO WHILE (Nstars .le. Nmax)
				read (14,end=6) idstars(Nstars), xstars(Nstars), ystars(Nstars), zstars(Nstars), vxstars(Nstars), &
				& vystars(Nstars), vzstars(Nstars), mass_stars(Nstars), age_stars(Nstars)
				Nstars = Nstars + 1
			end do
		else
			idstars(:) = 0
			xstars(:) = 0.0_8
			ystars(:) = 0.0_8
			zstars(:) = 0.0_8
			vxstars(:) = 0.0_8
			vystars(:) = 0.0_8
			vzstars(:) = 0.0_8
			initial_mass_stars(:) = 0.0_8
			age_stars(:) = 0.0_4
			SNII_stars(:) = 0.0_4
			SNI_stars(:) = 0.0_4

			DO WHILE (Nstars .le. Nmax)
				read (14,end=6) idstars(Nstars), xstars(Nstars), ystars(Nstars), zstars(Nstars), vxstars(Nstars), &
				& vystars(Nstars), vzstars(Nstars), initial_mass_stars(Nstars), age_stars(Nstars), SNII_stars(Nstars), SNI_stars(Nstars)
				Nstars = Nstars + 1
			end do
		end if
 6   		continue
		close (14)
		if( Nstars.eq.Nmax ) then
			print *, 'DATA STARS - error reading data: too many star particles. Change Nmax.'
			print *, 'Nstars=',Nstars
			stop
		end if
		Nstars = Nstars-1
		print *,Nstars
	end subroutine data_stars
end module read_binary
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module splitter
use parameters
	implicit none
		real(4) :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,vec(3),vec_vel(3)
		real(4),allocatable :: xprime(:),yprime(:),zprime(:),vxprime(:),vyprime(:),vzprime(:),rprime(:)
		real(4) :: saxis1(3),saxis2(3),saxis3(3),theta,phi
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

		if(allocated(xprime)) deallocate(xprime)
		if(allocated(yprime)) deallocate(yprime)
		if(allocated(zprime)) deallocate(zprime)
		if(allocated(vxprime)) deallocate(vxprime)
		if(allocated(vyprime)) deallocate(vyprime)
		if(allocated(vzprime)) deallocate(vzprime)
		if(allocated(rprime)) deallocate(rprime)
		allocate( xprime(nsplit),yprime(nsplit),zprime(nsplit),vxprime(nsplit),vyprime(nsplit),vzprime(nsplit),rprime(nsplit) )
	end subroutine allocation

!!!	DEALLOCATE PRIMED VALUES
!!!	----------------------------------------
	subroutine deallocate_primes()
	implicit none

		if(allocated(xprime)) deallocate(xprime)
		if(allocated(yprime)) deallocate(yprime)
		if(allocated(zprime)) deallocate(zprime)
		if(allocated(vxprime)) deallocate(vxprime)
		if(allocated(vyprime)) deallocate(vyprime)
		if(allocated(vzprime)) deallocate(vzprime)
		if(allocated(rprime)) deallocate(rprime)
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module disc_properties
use parameters
use read_binary
use splitter
	implicit none
	real(4) :: Rhalf_spher, vcom(3), rcom(3), ang_mom(3), Rdisc, Hdisc
	real(4) :: ang_mom2(3), Rdisc2, Hdisc2
	real(8) :: L_mag
	integer :: closest_star(1)

contains
	subroutine half_mass_radius_spherical(Rvir)
!!! Calculates the spherical half mass radius for gas colder than max_T (and maybe also stars younger than max_age) within fRvir
		implicit none
		real(4),intent(in) :: Rvir
		real(4),allocatable :: rad(:)
		real(8),allocatable :: mass(:)
		real(8) :: mtot
		integer :: i,j,n

		print *, 'calculating half mass radius'
		n = ceiling(fRvir * Rvir / 0.05_4)	!!! 50 pc bins
		allocate( rad(n+1), mass(n+1) )
		mass(:) = 0.0_8
		rad(1:n) = (/ (i,i=1,n) /) * (fRvir * Rvir) / n
		rad(n+1) = 2.0_4*rad(n)

		do i=1,Ngas
			if( temperature_gas(i) .le. max_T ) then
				if( sqrt(xgas(i)**2 + ygas(i)**2 + zgas(i)**2) .le. rad(n) ) then
					j=1
					do while( sqrt(xgas(i)**2 + ygas(i)**2 + zgas(i)**2) .gt. rad(j) )
						j = j+1
					end do
					mass(j:n+1) = mass(j:n+1) + 0.03363_8 * real(density_gas(i) * (cell_size_gas(i)**3), 8)
				end if
			end if
		end do
		if( use_stars ) then
			do i=1,Nstars
				if( age_stars(i) .le. max_age ) then
					if( sngl(sqrt(xstars(i)**2 + ystars(i)**2 + zstars(i)**2)) .le. rad(n) ) then
						j=1
						do while( sngl(sqrt(xstars(i)**2 + ystars(i)**2 + zstars(i)**2)) .gt. rad(j) )
							j = j+1
						end do
						mass(j:n+1) = mass(j:n+1) + mass_stars(i)
					end if
				end if
			end do
		end if
		mtot = mass(n) / 2.0_8	!!! half the mass

		j=2
		do while( mass(j) .lt. mtot .and. j .lt. n )
			j = j+1
		end do
!		Rhalf_spher = rad(j-1) + ( (rad(j)-rad(j-1))/(mass(j)-mass(j-1)) )*(mtot-mass(j-1))
		Rhalf_spher = 10**( log10(rad(j-1)) + ( log10(rad(j)/rad(j-1)) / sngl(log10(mass(j)/mass(j-1))) ) * sngl(log10(mtot/mass(j-1))) )
		if( Rhalf_spher .lt. rad(j-1) .or. Rhalf_spher .gt. rad(j) ) then
			Rhalf_spher = 0.5_4*(rad(j-1)+rad(j))
		end if

		write(16,*) 'Rhalf_spher = ',Rhalf_spher
		write(16,*) ''

		deallocate( rad, mass )
	end subroutine half_mass_radius_spherical
!________________________________________________________________________________________

	subroutine center(radius, height, zed_prime)
! Calculate COM and COM velocity within specified volume.
! If 'height'=0, volume is sphere of 'radius'.
! If 'height'>0, volume is cylinder of radius 'radius' and half thickness equal to 'height', with axis 'zed_prime'.

		implicit none
		real(4),intent(in) :: radius, height,zed_prime(3)
		real(4) :: x_prime, y_prime, z_prime
		real(8) :: mass, rcm(3), vcm(3)
		integer :: i

		mass   = 0.0_8
		rcm(:) = 0.0_8
		vcm(:) = 0.0_8

		saxis3(:) = zed_prime(:)
		call axes(saxis1,saxis2,saxis3)

		do i=1,Ngas
			if( temperature_gas(i) .le. max_T ) then
				if ( height .eq. 0.0_4 ) then
					if( sqrt(xgas(i)**2 + ygas(i)**2 + zgas(i)**2) .le. radius ) then
						mass   = mass   + 0.03363_8 * real(density_gas(i) * (cell_size_gas(i)**3), 8)
						rcm(1) = rcm(1) + 0.03363_8 * real(density_gas(i) * (cell_size_gas(i)**3) * xgas(i), 8)
						rcm(2) = rcm(2) + 0.03363_8 * real(density_gas(i) * (cell_size_gas(i)**3) * ygas(i), 8)
						rcm(3) = rcm(3) + 0.03363_8 * real(density_gas(i) * (cell_size_gas(i)**3) * zgas(i), 8)
						vcm(1) = vcm(1) + 0.03363_8 * real(density_gas(i) * (cell_size_gas(i)**3) * vxgas(i), 8)
						vcm(2) = vcm(2) + 0.03363_8 * real(density_gas(i) * (cell_size_gas(i)**3) * vygas(i), 8)
						vcm(3) = vcm(3) + 0.03363_8 * real(density_gas(i) * (cell_size_gas(i)**3) * vzgas(i), 8)
					end if
				else
					x_prime = xgas(i)*saxis1(1) + ygas(i)*saxis1(2) + zgas(i)*saxis1(3)
					y_prime = xgas(i)*saxis2(1) + ygas(i)*saxis2(2) + zgas(i)*saxis2(3)
					z_prime = xgas(i)*saxis3(1) + ygas(i)*saxis3(2) + zgas(i)*saxis3(3)
					if( sqrt(x_prime**2 + y_prime**2) .le. radius .and. abs(z_prime) .le. height ) then
						mass   = mass   + 0.03363_8 * real(density_gas(i) * (cell_size_gas(i)**3), 8)
						rcm(1) = rcm(1) + 0.03363_8 * real(density_gas(i) * (cell_size_gas(i)**3) * xgas(i), 8)
						rcm(2) = rcm(2) + 0.03363_8 * real(density_gas(i) * (cell_size_gas(i)**3) * ygas(i), 8)
						rcm(3) = rcm(3) + 0.03363_8 * real(density_gas(i) * (cell_size_gas(i)**3) * zgas(i), 8)
						vcm(1) = vcm(1) + 0.03363_8 * real(density_gas(i) * (cell_size_gas(i)**3) * vxgas(i), 8)
						vcm(2) = vcm(2) + 0.03363_8 * real(density_gas(i) * (cell_size_gas(i)**3) * vygas(i), 8)
						vcm(3) = vcm(3) + 0.03363_8 * real(density_gas(i) * (cell_size_gas(i)**3) * vzgas(i), 8)
					end if
				end if
			end if
		end do

		if( use_stars ) then
			do i=1,Nstars
				if( age_stars(i) .le. max_age ) then
					if (height.eq.0) then
						if( sngl(sqrt(xstars(i)**2 + ystars(i)**2 + zstars(i)**2)) .le. radius ) then
							mass = mass + mass_stars(i)
							rcm(1) = rcm(1) + mass_stars(i)*xstars(i)
							rcm(2) = rcm(2) + mass_stars(i)*ystars(i)
							rcm(3) = rcm(3) + mass_stars(i)*zstars(i)
							vcm(1) = vcm(1) + mass_stars(i)*vxstars(i)
							vcm(2) = vcm(2) + mass_stars(i)*vystars(i)
							vcm(3) = vcm(3) + mass_stars(i)*vzstars(i)
						end if
					else
						x_prime = sngl(xstars(i))*saxis1(1) + sngl(ystars(i))*saxis1(2) + sngl(zstars(i))*saxis1(3)
						y_prime = sngl(xstars(i))*saxis2(1) + sngl(ystars(i))*saxis2(2) + sngl(zstars(i))*saxis2(3)
						z_prime = sngl(xstars(i))*saxis3(1) + sngl(ystars(i))*saxis3(2) + sngl(zstars(i))*saxis3(3)
						if( sqrt(x_prime**2 + y_prime**2) .le. radius .and. abs(z_prime) .le. height ) then
							mass   = mass   + mass_stars(i)
							rcm(1) = rcm(1) + mass_stars(i)*xstars(i)
							rcm(2) = rcm(2) + mass_stars(i)*ystars(i)
							rcm(3) = rcm(3) + mass_stars(i)*zstars(i)
							vcm(1) = vcm(1) + mass_stars(i)*vxstars(i)
							vcm(2) = vcm(2) + mass_stars(i)*vystars(i)
							vcm(3) = vcm(3) + mass_stars(i)*vzstars(i)
						end if
					end if
				end if
			end do
		end if

		rcm(:) = rcm(:) / mass
		vcm(:) = vcm(:) / mass
		rcom(:) = rcm(:)
		vcom(:) = vcm(:)
		write(16,*) 'rcm = ',rcm
		write(16,*) 'vcm = ',vcm
	end subroutine center
!________________________________________________________________________________________

	subroutine disc_spin(radius, height, zed_prime)
! Calculate angular momentum of cold gas (and young stars) within specified volume.
! If 'height'=0, volume is sphere of 'radius'.
! If 'height'>0, volume is cylinder of radius 'radius' and half thickness equal to 'height', with axis 'zed_prime'.

		implicit none
		real(4),intent(in) :: radius, height, zed_prime(3)
		real(4) :: rcen(3), vcen(3)
		real(8) :: ang_mom_temp(3), mass_temp
		integer :: i

		print *, 'calculating angular momentum'
		ang_mom_temp(:) = 0.0_8
		mass_temp = 0.0_8
		saxis3(:) = zed_prime
		call axes(saxis1,saxis2,saxis3)
		if( move_center ) then
			rcen(:) = rcom(:)
			vcen(:) = vcom(:)
		else
			rcen(:) = 0.0_4
			vcen(:) = 0.0_4
		end if

		do i=1,Ngas
			if( temperature_gas(i) .le. max_T ) then
		      		if( height .eq. 0.0_4 ) then
					if ( sqrt(xgas(i)**2 + ygas(i)**2 + zgas(i)**2) .le. radius ) then
						ang_mom_temp(1) = ang_mom_temp(1) + 0.03363_8 * real(density_gas(i) * ( (cell_size_gas(i))**3 ) * ( (ygas(i)-rcen(2))*(vzgas(i)-vcen(3)) - (zgas(i)-rcen(3))*(vygas(i)-vcen(2)) ), 8)
						ang_mom_temp(2) = ang_mom_temp(2) + 0.03363_8 * real(density_gas(i) * ( (cell_size_gas(i))**3 ) * ( (zgas(i)-rcen(3))*(vxgas(i)-vcen(1)) - (xgas(i)-rcen(1))*(vzgas(i)-vcen(3)) ), 8)
						ang_mom_temp(3) = ang_mom_temp(3) + 0.03363_8 * real(density_gas(i) * ( (cell_size_gas(i))**3 ) * ( (xgas(i)-rcen(1))*(vygas(i)-vcen(2)) - (ygas(i)-rcen(2))*(vxgas(i)-vcen(1)) ), 8)
						mass_temp       = mass_temp       + 0.03363_8 * real(density_gas(i) * ( (cell_size_gas(i))**3 ), 8)
					end if
      				else
					call split0( xgas(i), ygas(i), zgas(i), rcen(:) )
					if( rprime(1) .le. radius .and. abs(zprime(1)) .le. height ) then
						ang_mom_temp(1) = ang_mom_temp(1) + 0.03363_8 * real(density_gas(i) * ( (cell_size_gas(i))**3 ) * ( (ygas(i)-rcen(2))*(vzgas(i)-vcen(3)) - (zgas(i)-rcen(3))*(vygas(i)-vcen(2)) ), 8)
						ang_mom_temp(2) = ang_mom_temp(2) + 0.03363_8 * real(density_gas(i) * ( (cell_size_gas(i))**3 ) * ( (zgas(i)-rcen(3))*(vxgas(i)-vcen(1)) - (xgas(i)-rcen(1))*(vzgas(i)-vcen(3)) ), 8)
						ang_mom_temp(3) = ang_mom_temp(3) + 0.03363_8 * real(density_gas(i) * ( (cell_size_gas(i))**3 ) * ( (xgas(i)-rcen(1))*(vygas(i)-vcen(2)) - (ygas(i)-rcen(2))*(vxgas(i)-vcen(1)) ), 8)
						mass_temp       = mass_temp       + 0.03363_8 * real(density_gas(i) * ( (cell_size_gas(i))**3 ), 8)
					end if
				end if
			end if
		end do
		if( use_stars ) then
			do i=1,Nstars
				if( age_stars(i) .le. max_age ) then
			      		if( height .eq. 0.0_4 ) then
						if ( sngl(sqrt(xstars(i)**2 + ystars(i)**2 + zstars(i)**2)) .le. radius ) then
							ang_mom_temp(1) = ang_mom_temp(1) + mass_stars(i) * ( (ystars(i)-real(rcen(2),8))*(vzstars(i)-real(vcen(3),8)) - (zstars(i)-real(rcen(3),8))*(vystars(i)-real(vcen(2),8)) )
							ang_mom_temp(2) = ang_mom_temp(2) + mass_stars(i) * ( (zstars(i)-real(rcen(3),8))*(vxstars(i)-real(vcen(1),8)) - (xstars(i)-real(rcen(1),8))*(vzstars(i)-real(vcen(3),8)) )
							ang_mom_temp(3) = ang_mom_temp(3) + mass_stars(i) * ( (xstars(i)-real(rcen(1),8))*(vystars(i)-real(vcen(2),8)) - (ystars(i)-real(rcen(2),8))*(vxstars(i)-real(vcen(1),8)) )
							mass_temp       = mass_temp       + mass_stars(i)
						end if
      					else
						call split0( sngl(xstars(i)), sngl(ystars(i)), sngl(zstars(i)), rcen(:) )
						if( rprime(1) .le. radius .and. abs(zprime(1)) .le. height ) then
							ang_mom_temp(1) = ang_mom_temp(1) + mass_stars(i) * ( (ystars(i)-real(rcen(2),8))*(vzstars(i)-real(vcen(3),8)) - (zstars(i)-real(rcen(3),8))*(vystars(i)-real(vcen(2),8)) )
							ang_mom_temp(2) = ang_mom_temp(2) + mass_stars(i) * ( (zstars(i)-real(rcen(3),8))*(vxstars(i)-real(vcen(1),8)) - (xstars(i)-real(rcen(1),8))*(vzstars(i)-real(vcen(3),8)) )
							ang_mom_temp(3) = ang_mom_temp(3) + mass_stars(i) * ( (xstars(i)-real(rcen(1),8))*(vystars(i)-real(vcen(2),8)) - (ystars(i)-real(rcen(2),8))*(vxstars(i)-real(vcen(1),8)) )
							mass_temp       = mass_temp       + mass_stars(i)
						end if
					end if
				end if
			end do
		end if

		L_mag = sqrt(dot_product(ang_mom_temp,ang_mom_temp))
		ang_mom_temp(:) = ang_mom_temp(:) / L_mag
		ang_mom(:) = sngl( ang_mom_temp(:) )
		L_mag = L_mag / mass_temp
		write(16,*) 'Angular momentum = ',ang_mom(:)
	end subroutine disc_spin
!________________________________________________________________________________________

	recursive subroutine disc_radius(counter)
! Calculates disc dimensions in cylindrical coordinates where z-axis is disc AM axis.
! Radius is defined as radius containing 85% of the mass of cold gas (and young stars) within a cylinder with R = 0.15*Rvir and H = +/- 1kpc. (Rvir is a global variable)
! Thickness is defined as thickness containing 85% of the mass of cold gas (and young stars) within a cylinder with R = Rdisc and H = +/- Rdisc.
! Calculation is recursive, and updates rest frame velocity as well as angular momentum of cold gas in the new volume 
! and then re-computes dimensions in new frame.

		implicit none
		integer,intent(in) :: counter
		real(4),allocatable :: rad(:), zed(:)
		real(8),allocatable :: mass(:), mass_edge(:)
		real(4) :: height, condition, rcen(3), vcen(3), split
		real(8) :: Mdisc_edge, Mdisc_face
		integer :: i, j, k, m, n1, n2, n, Ngas_list, Nstar_list, ind
		integer,allocatable :: gas_list(:), star_list(:)

		print *, 'calculating disc radius'
		if(allocated(rad)) deallocate(rad)
		if(allocated(zed)) deallocate(zed)
		if(allocated(mass)) deallocate(mass)
		if(allocated(mass_edge)) deallocate(mass_edge)
		if(allocated(gas_list)) deallocate(gas_list)
		if(allocated(star_list)) deallocate(star_list)

		n1 = ceiling(fRvir * Rvir / 0.05_4)	!!! 50 pc bins
		n = n1+1
		print *, 'n1=',n1
		allocate ( rad(n), zed(n), mass(n), mass_edge(n) )

		rad(1:n1) = (/ (i,i=1,n1) /) * (fRvir * Rvir) / n1
		rad(n) = 2.0_4*rad(n1)
		mass(:) = 0.0_8
		mass_edge(:) = 0.0_8
		zed(:) =  0.0_4

		saxis3(:) = ang_mom(:)
		call axes(saxis1,saxis2,saxis3)
!		height = min(Hdisc2,1.0_4)
		height = initial_height
		if( move_center ) then
			rcen(:) = rcom(:)
			vcen(:) = vcom(:)
		else
			rcen(:) = 0.0_4
			vcen(:) = 0.0_4
		end if

		allocate( gas_list(Ngas) )
		gas_list(:) = 0
		Ngas_list = 0
		do i=1,Ngas
			if( temperature_gas(i) .le. max_T ) then
				if ( xgas(i)**2 + ygas(i)**2 + zgas(i)**2 .le. 3.0_4*(rad(n1)**2) ) then
					Ngas_list = Ngas_list+1
					gas_list(Ngas_list) = i
					if(cell_size_gas(i) .lt. 1.5_4*res) then
						call split0(xgas(i), ygas(i), zgas(i), rcen(:))
					else if(cell_size_gas(i) .lt. 3.0_4*res) then
						call split1(xgas(i), ygas(i), zgas(i), rcen(:), cell_size_gas(i))
					else
						call split2(xgas(i), ygas(i), zgas(i), rcen(:), cell_size_gas(i))
					end if
					k = size(xprime(:))
					split = log(sngl(k))/log(8.0_4)
					do j=1,k
						if( abs(zprime(j)) .le. height .and. rprime(j) .le. rad(n1) ) then
							m = 1
							do while( rprime(j) .gt. rad(m) )
								m = m + 1
							end do
							mass(m:n) = mass(m:n) + 0.03363_8 * real(density_gas(i)*( (cell_size_gas(i)/(2.0_4**split))**3 ), 8)
						end if
					end do
				end if
				call deallocate_primes()
			end if
		end do

		if( use_stars ) then
			allocate( star_list(Nstars) )
			star_list(:) = 0
			Nstar_list = 0
			do i=1,Nstars
				if( age_stars(i) .le. max_age ) then
					if ( sngl(xstars(i)**2 + ystars(i)**2 + zstars(i)**2) .le. 3.0_4*(rad(n1)**2) ) then
						Nstar_list = Nstar_list+1
						star_list(Nstar_list) = i
						call split0(sngl(xstars(i)), sngl(ystars(i)), sngl(zstars(i)), rcen(:))
						if( abs(zprime(1)) .le. height .and. rprime(1) .le. rad(n1) ) then
							m = 1
							do while( rprime(1) .gt. rad(m) )
								m = m + 1
							end do
							mass(m:n) = mass(m:n) + mass_stars(i)
						end if
					end if
					call deallocate_primes()
				end if
			end do
		end if

		Mdisc_face = fM * mass(n1)
		j=2
		do while( mass(j) .lt. Mdisc_face .and. j .lt. n1 )
			j = j+1
		end do

!		Rdisc = rad(j-1) + ( (rad(j)-rad(j-1))/(mass(j)-mass(j-1)) )*(Mdisc_face-mass(j-1))
		Rdisc = 10**( log10(rad(j-1)) + ( log10(rad(j)/rad(j-1)) / log10(sngl(mass(j)/mass(j-1))) ) * log10(sngl(Mdisc_face/mass(j-1))) )
		if( Rdisc .lt. rad(j-1) .or. Rdisc .gt. rad(j) ) then
			Rdisc = 0.5_4*(rad(j-1)+rad(j))
		end if
		Rdisc = max(Rdisc, 0.05_4)			!!! minimum allowed disc radius is 50pc
!		Rdisc = max(Rdisc, 1.0_4)			!!! minimum allowed disc radius is 1 kpc
		Rdisc = min(Rdisc, fRvir*Rvir)			!!! maximum allowed disc radius is fRvir*Rvir

		n2 = ceiling(fRd * Rdisc / 0.05_4)	!!! 50 pc bins
		print *, 'n2=',n2
		zed(1:n2) = (/ (i,i=1,n2) /) * (fRd * Rdisc) / n2
		zed(n2+1:n) = 2.0_4*zed(n2)

		do ind=1,Ngas_list
			i = gas_list(ind)
			if ( xgas(i)**2 + ygas(i)**2 + zgas(i)**2 .le. 3.0_4*(Rdisc**2) ) then
				if( cell_size_gas(i) .lt. 1.5_4*res ) then
					call split0(xgas(i), ygas(i), zgas(i), rcen(:))
				else if( cell_size_gas(i) .lt. 3.0_4*res) then
					call split1(xgas(i), ygas(i), zgas(i), rcen(:), cell_size_gas(i))
				else
					call split2(xgas(i), ygas(i), zgas(i), rcen(:), cell_size_gas(i))
				end if
				k = size(xprime(:))
				split = log(sngl(k))/log(8.0_4)
				do j=1,k
					if(rprime(j) .le. Rdisc .and. abs(zprime(j)) .le. zed(n2)) then
						m=1
						do while( abs(zprime(j)) .gt. zed(m))
							m = m+1
						end do
						mass_edge(m:n) = mass_edge(m:n) + 0.03363_8 * real(density_gas(i) * ((cell_size_gas(i)/(2.0_4**split))**3 ), 8)
					end if      
				end do
				call deallocate_primes()
			end if
		end do

		if( use_stars ) then
			do ind=1,Nstar_list
				i = star_list(ind)
				if ( sngl(xstars(i)**2 + ystars(i)**2 + zstars(i)**2) .le. 3.0_4*(Rdisc**2) ) then
					call split0(sngl(xstars(i)), sngl(ystars(i)), sngl(zstars(i)), rcen(:))
					if(rprime(1) .le. Rdisc .and. abs(zprime(1)) .le. zed(n2)) then
						m=1
						do while( abs(zprime(1)) .gt. zed(m))
							m = m+1
						end do
						mass_edge(m:n) = mass_edge(m:n) + mass_stars(i)
					end if      
					call deallocate_primes()
				end if
			end do
		end if

		Mdisc_edge = fM * mass_edge(n2)
		j = 2
		do while(mass_edge(j) .lt. Mdisc_edge .and. j .lt. n2)
			j = j + 1
		end do
!		Hdisc = zed(j-1) + ( (zed(j)-zed(j-1))/(mass_edge(j)-mass_edge(j-1)) )*(Mdisc_edge-mass_edge(j-1))
		Hdisc = 10**( log10(zed(j-1)) + ( log10(zed(j)/zed(j-1)) / log10(sngl(mass_edge(j)/mass_edge(j-1))) ) * log10(sngl(Mdisc_edge/mass_edge(j-1))) )
		if( Hdisc .lt. zed(j-1) .or. Hdisc .gt. zed(j) ) then
			Hdisc = 0.5_4*(zed(j-1)+zed(j))
		end if
		Hdisc = max(Hdisc, 0.05_4)			!!! minimum allowed disc height is 50pc
!		Hdisc = max(Hdisc, 0.5_4)			!!! minimum allowed disc height is 0.5 kpc
		Hdisc = min(Hdisc, fRd*Rdisc)			!!! maximum allowed disc height is fRd*Rdisc

		write(16,*) 'Disc Radius = ',Rdisc
		write(16,*) 'Disc Gas Thickness = ',Hdisc
		write(16,*) '--------------------------------------------'

		deallocate( gas_list, star_list )
		condition = max( abs(Rdisc-Rdisc2)/Rdisc, abs(Hdisc-Hdisc2)/Hdisc, abs((ang_mom(1)-ang_mom2(1))/ang_mom(1)), abs((ang_mom(2)-ang_mom2(2))/ang_mom(2)), abs((ang_mom(3)-ang_mom2(3))/ang_mom(3)) )
		if( condition .gt. 0.05 .and. counter .lt. 5) then
			Rdisc2 = Rdisc
			Hdisc2 = Hdisc
			ang_mom2(:) = ang_mom(:)
			deallocate( rad, zed, mass, mass_edge )
			call center(Rdisc2, Hdisc2, ang_mom2(:))
			call disc_spin(Rdisc2, Hdisc2, ang_mom2(:))
			call disc_radius(counter+1)
		else
			write(16,*) 'counter Rdisc = ',counter
			write(16,*) 'r,    m(r),    z,    m(z)'
			do j=1,max(n1,n2)
				write(16,*) rad(j), mass(j), zed(j), mass_edge(j)
			end do
			write(16,*) ''
			deallocate( rad, zed, mass, mass_edge )
		end if
	end subroutine disc_radius
!________________________________________________________________________________________

	subroutine mean_disc_values(radius, height, zed_prime, rad_spher, dm_file_name, stars_file_name, stars_z_file_name)
! Calculates mean baryonic properties for disc.
! Disc is cylinder with radius 'radius' and height 'height' with axis 'zed_prime'
! Spherical properties are then calculated in sphere of radius 'rad_spher'
	implicit none
	character(len=256) :: stars_file_name, stars_z_file_name, dm_file_name
	real(4),intent(in) :: radius, height, zed_prime(3), rad_spher
	integer :: i, j, k, ntmax(3)
	real(4) :: split, rcen(3), vcen(3), Ltmax(3), loop_rad
	real(8) :: added_mass, added_initial_mass
	real(8),allocatable :: SFR_tmax_disc(:,:), SFR_tmax_sphere(:,:), SFR_tmax_1kpc(:), SFR_tmax_015Rv(:)
	real(4),allocatable :: tmax(:,:)

     	mean_mgas = 0.0_8
     	mean_mcold = 0.0_8
     	mean_mstars = 0.0_8
     	mean_mdm = 0.0_8
     	mean_SFR = 0.0_8
     	mean_age = 0.0_8
     	mean_metgas = 0.0_8
     	mean_metstars = 0.0_8
	exsitu_stellar_mass = 0.0_8

     	spher_mgas = 0.0_8
     	spher_mcold = 0.0_8
     	spher_mstars = 0.0_8
     	spher_initial_mstars = 0.0_8
     	spher_mdm = 0.0_8
     	spher_SFR = 0.0_8
     	spher_age = 0.0_8
     	spher_metgas = 0.0_8
     	spher_metstars = 0.0_8
	spher_exsitu_stellar_mass = 0.0_8

	initial_stellar_mass_disc = 0.0_8
	initial_stellar_mass_sphere = 0.0_8
	analytic_stellar_mass_disc = 0.0_8
	analytic_stellar_mass_sphere = 0.0_8
	initial_stellar_mass_disc_es = 0.0_8
	initial_stellar_mass_sphere_es = 0.0_8
	analytic_stellar_mass_disc_es = 0.0_8
	analytic_stellar_mass_sphere_es = 0.0_8
	SFR_test_disc(:) = 0.0_8
	SFR_test_sphere(:) = 0.0_8

	Mvir = 0.0_8
	Vvir = 0.0_8
	MgV = 0.0_8
	MsV = 0.0_8
	IMsV = 0.0_8
	MdV = 0.0_8
	Mg0_1 = 0.0_8
	Ms0_1 = 0.0_8
	IMs0_1 = 0.0_8
	Md0_1 = 0.0_8
	Mg10 = 0.0_8
	Ms10 = 0.0_8
	IMs10 = 0.0_8
	Md10 = 0.0_8

	Mgas_015Rv = 0.0_8
	Mstars_015Rv = 0.0_8
	initial_Mstars_015Rv = 0.0_8
	Mcold_gas_015Rv = 0.0_8
	Myoung_stars_015Rv = 0.0_8
	Mdm_015Rv = 0.0_8
	SFR_015Rv = 0.0_8

	Mgas_1kpc = 0.0_8
	Mstars_1kpc = 0.0_8
	initial_Mstars_1kpc = 0.0_8
	Mcold_gas_1kpc = 0.0_8
	Myoung_stars_1kpc = 0.0_8
	Mdm_1kpc = 0.0_8
	SFR_1kpc = 0.0_8

	if( move_center ) then
		rcen(:) = rcom(:)
		vcen(:) = vcom(:)
	else
		rcen(:) = 0.0_4
		vcen(:) = 0.0_4
	end if

	Ltmax(:) = tmax_f(:) - tmax_i(:)
	ntmax(:) = floor(Ltmax(:) / dtmax) + 1
	allocate( SFR_tmax_disc(3,maxval(ntmax(1:3))), SFR_tmax_sphere(3,maxval(ntmax(1:3))), SFR_tmax_1kpc(ntmax(1)), SFR_tmax_015Rv(ntmax(1)), tmax(3,maxval(ntmax(1:3))) )
	do j=1,3
		tmax(j,1:ntmax(j)) = (tmax_i(j) - dtmax) + (/ (i,i=1,ntmax(j)) /)*dtmax	!!! Myr
	end do
	tmax(:,:) = tmax(:,:) / 1000.0_4	!!! Gyr
	SFR_tmax_disc(:,:) = 0.0_8
	SFR_tmax_sphere(:,:) = 0.0_8
	SFR_tmax_1kpc(:) = 0.0_8
	SFR_tmax_015Rv(:) = 0.0_8

	saxis3(:) = zed_prime(:)
	call axes(saxis1,saxis2,saxis3)
	do i=1,Ngas
		loop_rad = sqrt(xgas(i)**2 + ygas(i)**2 + zgas(i)**2)

		if ( loop_rad .le. sqrt(3.0_4)*radius ) then
			if(cell_size_gas(i) < 1.5_4*res) then
				call split0( xgas(i), ygas(i), zgas(i), rcen(:) )
			else if(cell_size_gas(i) < 3.0_4*res) then
				call split1( xgas(i), ygas(i), zgas(i), rcen(:), cell_size_gas(i) )
			else 
				call split2( xgas(i), ygas(i), zgas(i), rcen(:), cell_size_gas(i))
			end if
			k = size(xprime(:))
			split = log(sngl(k)) / log(8.0_4)
			added_mass = 0.03363_8 * real(density_gas(i) * ( cell_size_gas(i) / (2.0_4**split) )**3, 8)
			do j=1,k
				if( rprime(j) .le. radius .and. abs(zprime(j)) .le. height ) then
					mean_mgas   = mean_mgas   + added_mass
					mean_metgas = mean_metgas + added_mass * real( SNII_gas(i), 8 )
					if( temperature_gas(i) .lt. max_T ) then
						mean_mcold = mean_mcold + added_mass
					end if
				end if
			end do
		end if

		if( loop_rad .le. Rvir ) then
			added_mass = 0.03363_8 * real( density_gas(i) * (cell_size_gas(i)**3), 8 )
			MgV = MgV + added_mass
			if( loop_rad .le. fRvir*Rvir ) then
				Mgas_015Rv = Mgas_015Rv + added_mass
				if( temperature_gas(i) .le. max_T ) then
					Mcold_gas_015Rv = Mcold_gas_015Rv + added_mass
				end if
				if( loop_rad .le. 0.1_4*Rvir ) then
					Mg0_1 = Mg0_1 + added_mass
				end if
				if( loop_rad .le. rad_spher ) then
					spher_mgas   = spher_mgas   + added_mass
					spher_metgas = spher_metgas + added_mass * real( SNII_gas(i), 8 )
					if( temperature_gas(i) .le. max_T ) then
						spher_mcold = spher_mcold + added_mass
					end if
				end if
			end if
		end if

		if( loop_rad .le. 10.0_4 ) then
			added_mass = 0.03363_8 * real( density_gas(i) * (cell_size_gas(i)**3), 8 )
			Mg10 = Mg10 + added_mass
			if( loop_rad .le. 1.0_4 ) then
				Mgas_1kpc = Mgas_1kpc + added_mass
				if( temperature_gas(i) .le. max_T ) then
					Mcold_gas_1kpc = Mcold_gas_1kpc + added_mass
				end if
			end if
		end if
	end do
	Mvir = Mvir + MgV
	spher_metgas =  spher_metgas * (1.d12) / (0.755_8 * 16.0_8 * 2.0_8)	!!! [O/H] + 12
	mean_metgas =  mean_metgas * (1.d12) / (0.755_8 * 16.0_8 * 2.0_8)	!!! [O/H] + 12
	deallocate( xgas, ygas, zgas, vxgas, vygas, vzgas, density_gas, temperature_gas, cell_size_gas, SNII_gas )

	allocate( xdm(Nmax), ydm(Nmax), zdm(Nmax), mass_dm(Nmax), vxdm(Nmax), vydm(Nmax), vzdm(Nmax), iddm(Nmax) )
	call data_dm(dm_file_name)
	do i=1,Ndm
		loop_rad = sngl( sqrt(xdm(i)**2 + ydm(i)**2 + zdm(i)**2) )
		added_mass = mass_dm(i)
		if ( loop_rad .le. sqrt(3.0_4)*radius ) then
			call split0( sngl(xdm(i)), sngl(ydm(i)), sngl(zdm(i)), rcen(:), sngl(vxdm(i)), sngl(vydm(i)), sngl(vzdm(i)), vcen(:) )
			if( rprime(1) .le. radius .and. abs(zprime(1)) .le. height ) then
				if( xprime(1)*vyprime(1) - yprime(1)*vxprime(1) .gt. fj * sqrt( vxprime(1)**2 + vyprime(1)**2 + vzprime(1)**2 ) * sqrt( xprime(1)**2 + yprime(1)**2 + zprime(1)**2 ) ) then
					mean_mdm = mean_mdm + added_mass
				end if
			end if
		end if

		if( loop_rad .le. Rvir ) then
			MdV = MdV + added_mass
			if( loop_rad .le. fRvir*Rvir ) then
				Mdm_015Rv = Mdm_015Rv + added_mass
				if( loop_rad .le. 0.1_4*Rvir ) then
					Md0_1 = Md0_1 + added_mass
				end if
				if( loop_rad .le. rad_spher ) then
					spher_mdm = spher_mdm + added_mass
				end if
			end if
		end if

		if( loop_rad .le. 10.0_4 ) then
			Md10 = Md10 + added_mass
			if( loop_rad .le. 1.0_4 ) then
				Mdm_1kpc = Mdm_1kpc + added_mass
			end if
		end if			
	end do
	Mvir = Mvir + MdV
	deallocate( vxdm, vydm, vzdm, iddm, xdm, ydm, zdm, mass_dm )

	if( .not. use_stars ) then
		allocate( xstars(Nmax), ystars(Nmax), zstars(Nmax), vxstars(Nmax), vystars(Nmax), vzstars(Nmax), idstars(Nmax), &
			& mass_stars(Nmax), initial_mass_stars(Nmax), age_stars(Nmax), SNII_stars(Nmax) )
		allocate( SNI_stars(Nmax) )
		call data_stars(stars_file_name,0)
		call data_stars(stars_z_file_name,1)
		deallocate( SNI_stars )
	end if

	do i=1,Nstars
		loop_rad = sngl( sqrt(xstars(i)**2 + ystars(i)**2 + zstars(i)**2) )
		added_mass = mass_stars(i)
		added_initial_mass = initial_mass_stars(i)

		call split0( sngl(xstars(i)), sngl(ystars(i)), sngl(zstars(i)), rcen(:), sngl(vxstars(i)), sngl(vystars(i)), sngl(vzstars(i)), vcen(:) )
		if(insitustars(idstars(i),1).eq.0) then
			if( rprime(1) .le. radius .and. abs(zprime(1)) .le. height ) then
				insitustars(idstars(i),1) = 1
			else
				insitustars(idstars(i),1) = 2
			end if
		end if
		if(insitustars(idstars(i),2).eq.0) then
			if( loop_rad .le. rad_spher ) then
				insitustars(idstars(i),2) = 1
			else
				insitustars(idstars(i),2) = 2
			end if
		end if

		if( rprime(1) .le. radius .and. abs(zprime(1)) .le. height ) then
			if( xprime(1)*vyprime(1) - yprime(1)*vxprime(1) .gt. fj * sqrt( vxprime(1)**2 + vyprime(1)**2 + vzprime(1)**2 ) * sqrt( xprime(1)**2 + yprime(1)**2 + zprime(1)**2 ) ) then
				mean_mstars   = mean_mstars   + added_mass
				mean_age      = mean_age      + added_mass * real(age_stars(i),8)
				mean_metstars = mean_metstars + added_mass * real(SNII_stars(i),8)
				initial_stellar_mass_disc  = initial_stellar_mass_disc  + initial_mass_stars(i)
				analytic_stellar_mass_disc = analytic_stellar_mass_disc + initial_mass_stars(i) * ( 1.0_8 - A0*log(1.0_8 + real(age_stars(i),8)/T0) )
				if( age_stars(i) .le. max_age ) then
					mean_mcold = mean_mcold + added_mass
				end if
				do j=1,3
					if( age_stars(i) .le. tmax(j,ntmax(j)) ) then
						k = 1
						do while(k .le. ntmax(j))
							if( age_stars(i) .le. tmax(j,k) ) then
								SFR_tmax_disc(j,k:ntmax(j)) = SFR_tmax_disc(j,k:ntmax(j)) + added_initial_mass
								k = 2*ntmax(j) + 5
							else
								k = k+1
							end if
						end do
					end if
				end do
				if( insitustars(idstars(i),1) .eq. 2 ) then
					exsitu_stellar_mass           = exsitu_stellar_mass + added_mass
					initial_stellar_mass_disc_es  = initial_stellar_mass_disc_es  + added_initial_mass
					analytic_stellar_mass_disc_es = analytic_stellar_mass_disc_es + added_initial_mass * ( 1.0_8 - A0*log(1.0_8 + real(age_stars(i),8)/T0) )
				end if
			end if
		end if

		if( loop_rad .le. Rvir ) then
			MsV = MsV + added_mass
			IMsV = IMsV + added_initial_mass
			if( loop_rad .le. fRvir*Rvir ) then
				Mstars_015Rv = Mstars_015Rv + added_mass
				initial_Mstars_015Rv = initial_Mstars_015Rv + added_initial_mass
				if( age_stars(i) .le. max_age ) then
					Myoung_stars_015Rv = Myoung_stars_015Rv + added_mass
					if( age_stars(i) .le. tmax(1,ntmax(1)) ) then
						k = 1
						do while(k .le. ntmax(1))
							if( age_stars(i) .le. tmax(1,k) ) then
								SFR_tmax_015Rv(k:ntmax(1)) = SFR_tmax_015Rv(k:ntmax(1)) + added_initial_mass
								k = 2*ntmax(1) + 5
							else
								k = k+1
							end if
						end do
					end if
				end if
				if( loop_rad .le. 0.1_4*Rvir ) then
					Ms0_1 = Ms0_1 + added_mass
					IMs0_1 = IMs0_1 + added_initial_mass
				end if
				if( loop_rad .le. rad_spher ) then
					spher_mstars   = spher_mstars   + added_mass
					spher_initial_mstars   = spher_initial_mstars   + added_initial_mass
					spher_age      = spher_age      + added_mass * real(age_stars(i),8)
					spher_metstars = spher_metstars + added_mass * real(SNII_stars(i),8)
					initial_stellar_mass_sphere  = initial_stellar_mass_sphere  + added_initial_mass
					analytic_stellar_mass_sphere = analytic_stellar_mass_sphere + added_initial_mass * ( 1.0_8 - A0*log(1.0_8 + real(age_stars(i),8)/T0) )
					if( age_stars(i) .le. max_age ) then
						spher_mcold = spher_mcold + added_mass
					end if
					do j=1,3
						if( age_stars(i) .le. tmax(j,ntmax(j)) ) then
							k = 1
							do while(k .le. ntmax(j))
								if( age_stars(i) .le. tmax(j,k) ) then
									SFR_tmax_sphere(j,k:ntmax(j)) = SFR_tmax_sphere(j,k:ntmax(j)) + added_initial_mass
									k = 2*ntmax(j) + 5
								else
									k = k+1
								end if
							end do
						end if
					end do
					if( insitustars(idstars(i),2) .eq. 2 ) then
						spher_exsitu_stellar_mass       = spher_exsitu_stellar_mass + added_mass
						initial_stellar_mass_sphere_es  = initial_stellar_mass_sphere_es  + added_initial_mass
						analytic_stellar_mass_sphere_es = analytic_stellar_mass_sphere_es + added_initial_mass * ( 1.0_8 - A0*log(1.0_8 + real(age_stars(i),8)/T0) )
					end if
				end if
			end if
		end if

		if( loop_rad .le. 10.0_4 ) then
			Ms10 = Ms10 + added_mass
			IMs10 = IMs10 + added_initial_mass
			if( loop_rad .le. 1.0_4 ) then
				Mstars_1kpc = Mstars_1kpc + added_mass
				initial_Mstars_1kpc = initial_Mstars_1kpc + added_initial_mass
				if( age_stars(i) .le. max_age ) then
					Myoung_stars_1kpc = Myoung_stars_1kpc + added_mass
					if( age_stars(i) .le. tmax(1,ntmax(1)) ) then
						k = 1
						do while(k .le. ntmax(1))
							if( age_stars(i) .le. tmax(1,k) ) then
								SFR_tmax_1kpc(k:ntmax(1)) = SFR_tmax_1kpc(k:ntmax(1)) + added_initial_mass
								k = 2*ntmax(1) + 5
							else
								k = k+1
							end if
						end do
					end if
				end if

			end if
		end if
	end do
	Mvir = Mvir + MsV
	Vvir = sqrt( (G*Mvir) / real(Rvir,8) )
	mean_metstars  = mean_metstars  * (1.d12) / (0.755_8 * 16.0_8 * 2.0_8)
	spher_metstars = spher_metstars * (1.d12) / (0.755_8 * 16.0_8 * 2.0_8)
	mean_age  = mean_age  * 1000.0_8
	spher_age = spher_age * 1000.0_8
	do j=1,3
		do k=1,ntmax(j)
			SFR_tmax_disc(j,k) = SFR_tmax_disc(j,k) / ( 1.d9 * real(tmax(j,k),8) )
			SFR_tmax_sphere(j,k) = SFR_tmax_sphere(j,k) / ( 1.d9 * real(tmax(j,k),8) )
		end do
		SFR_test_disc(j) = sum(SFR_tmax_disc(j,1:ntmax(j))) / (1.0_8 * ntmax(j))
		SFR_test_sphere(j) = sum(SFR_tmax_sphere(j,1:ntmax(j))) / (1.0_8 * ntmax(j))
	end do
	mean_SFR = SFR_test_disc(1)
	spher_SFR = SFR_test_sphere(1)
	SFR_tmax_1kpc(1:ntmax(1))  = SFR_tmax_1kpc(1:ntmax(1))  / (1.d9 * real( tmax(1,1:ntmax(1)), 8))
	SFR_tmax_015RV(1:ntmax(1)) = SFR_tmax_015RV(1:ntmax(1)) / (1.d9 * real( tmax(1,1:ntmax(1)), 8))
	SFR_1kpc  = sum(SFR_tmax_1kpc(1:ntmax(1)))  / (1.0_8 * ntmax(1))
	SFR_015RV = sum(SFR_tmax_015RV(1:ntmax(1))) / (1.0_8 * ntmax(1))
	deallocate( SFR_tmax_disc, SFR_tmax_sphere, SFR_tmax_1kpc, SFR_tmax_015Rv, tmax )

	closest_star = minloc(sqrt( xstars(1:Nstars)**2 + ystars(1:Nstars)**2 + zstars(1:Nstars)**2 ))
	print *, 'closest star', closest_star(1), idstars(closest_star(1))
	print *, 'closest star position', xstars(closest_star(1)), ystars(closest_star(1)), zstars(closest_star(1))

	mean_mbar = mean_mstars + mean_mgas
	mean_age = mean_age / mean_mstars
	mean_fgas = mean_mgas / mean_mbar
	mean_sigma = mean_mbar / real(pi * ( (1000.0_4 * radius)**2 ), 8)			!!! M_sun pc^-2
	mean_sigma_stars = mean_mstars / real(pi * ( (1000.0_4 * radius)**2 ), 8)		!!! M_sun pc^-2
	mean_sigma_gas = mean_mgas / real(pi * ( (1000.0_4 * radius)**2 ), 8)			!!! M_sun pc^-2
	mean_Sig_SFR = mean_SFR / real(pi * ( radius**2 ), 8)					!!! M_sun yr^-1 kpc^-2
	mean_SSFR = 1.d9 * mean_SFR / mean_mstars						!!! Gyr^-1
	mean_tau  = 1.d9 * mean_SFR / mean_mgas							!!! Gyr^-1
 	mean_metgas = log10( mean_metgas / mean_mgas )						!!! Log[O/H] + 12
	mean_metstars = log10( mean_metstars / mean_mstars )					!!! Log[O/H] + 12

	spher_mbar = spher_mstars + spher_mgas
	spher_age = spher_age / spher_mstars
	spher_fgas = spher_mgas / spher_mbar
	spher_sigma = spher_mbar / real(pi * ( (1000.0_4 * rad_spher)**2 ), 8)			!!! M_sun pc^-2
	spher_sigma_stars = spher_mstars / real(pi * ( (1000.0_4 * rad_spher)**2 ), 8)		!!! M_sun pc^-2
	spher_sigma_gas = spher_mgas / real(pi * ( (1000.0_4 * rad_spher)**2 ), 8)		!!! M_sun pc^-2
	spher_Sig_SFR = spher_SFR / real(pi * ( rad_spher**2 ), 8)				!!! M_sun yr^-1 kpc^-2
	spher_SSFR = 1.d9 * spher_SFR / spher_mstars						!!! Gyr^-1
	spher_tau  = 1.d9 * spher_SFR / spher_mgas						!!! Gyr^-1
 	spher_metgas = log10( spher_metgas / spher_mgas )					!!! Log[O/H] + 12
	spher_metstars = log10( spher_metstars / spher_mstars )					!!! Log[O/H] + 12

      end subroutine mean_disc_values
end module disc_properties

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program main
!#include "adef.h"
	use parameters
	use read_binary
	use splitter
	use disc_properties
!	use omp_lib
	implicit none
	integer :: i, j, k, Nsnapshot
	character(len=256),allocatable :: gas_file_name(:), stars_file_name(:), stars_z_file_name(:), dm_file_name(:)
	character(len=32),allocatable :: sim(:)
	character(len=256) :: filename, file_temp, input_arg
	character(len=6) :: snap_tag
	character(len=5) :: snap_tag2

	if(iargc().ge.1) then
		call getarg(1,input_arg)
		write(filename,'(a)') trim(input_arg)
	else
		write(filename,'(a)') './disc_input.dat'
	end if
	print *, 'input file'
	print *, trim(filename)

	open(unit=16,file=trim(filename))
	read(16,*) Nsnapshot
	print *, Nsnapshot
	allocate( gas_file_name(Nsnapshot), stars_file_name(Nsnapshot), stars_z_file_name(Nsnapshot), dm_file_name(Nsnapshot), sim(Nsnapshot) )
	do k=1,Nsnapshot
		read(16,*) gas_file_name(k), sim(k)

		write(file_temp,'(a256)') gas_file_name(k)
		i = scan(file_temp,'Z')
		write(stars_file_name(k),*) trim(file_temp(1:i-2)), 'S', trim(file_temp(i+1:i+15))
		write(stars_z_file_name(k),*) trim(file_temp(1:i-2)), 'SZ', trim(file_temp(i+1:i+15))
		write(dm_file_name(k),*) trim(file_temp(1:i-2)), 'D', trim(file_temp(i+1:i+15))
	end do
	close(unit=16)
	print *, trim(gas_file_name(1))
	print *, trim(dm_file_name(1))
	print *, trim(stars_file_name(1))
	print *, trim(stars_z_file_name(1))
	print *, trim(gas_file_name(Nsnapshot))
	print *, trim(dm_file_name(Nsnapshot))
	print *, trim(stars_file_name(Nsnapshot))
	print *, trim(stars_z_file_name(Nsnapshot))

	do k=1,Nsnapshot
		Rvir = 0.0_4

		write(file_temp,'(a256)') gas_file_name(k)
		i = scan(file_temp,'Z')
		read(file_temp(i+1:i+4),*) Rvir
		Rvir = Rvir / 4.0_4

		i = scan(file_temp,'.')
		if(Rvir < 250.0_4) then
			write(snap_tag,'(a6)') trim(adjustl(file_temp(i+1:i+6)))		!!! a0.***
			write(snap_tag2,'(a5)') trim(adjustl(file_temp(i+2:i+6)))		!!! 0.***
		else
			write(snap_tag,'(a6)') trim(adjustl(file_temp(i-2:i+3)))		!!! a0.***
			write(snap_tag2,'(a5)') trim(adjustl(file_temp(i-1:i+3)))		!!! 0.***
		end if
		print *, 'snap_tag=',snap_tag

		j = k
		do while(j .ge. 1 .and. trim(sim(j)) .eq. trim(sim(k)))
			j = j - 1
			if(j .eq. 0) exit
		end do
		i = k
		do while(i .le. Nsnapshot .and. trim(sim(i)) .eq. trim(sim(k)))
			i = i + 1
			if(i .eq. Nsnapshot + 1) exit
		end do
		print *, 'k',k
		print *, 'j',j
		print *, 'i',i
		if(j.eq.k-1) then
			call open_files(i,j,k)

			if(allocated(insitustars)) deallocate(insitustars)
			allocate( insitustars(Nmax/2,2) )
			insitustars(:,:) = 0
		end if
		print *, 'alocating gas arrays'
		allocate( xgas(Nmax), ygas(Nmax), zgas(Nmax), density_gas(Nmax), cell_size_gas(Nmax), SNII_gas(Nmax) )
		allocate( vxgas(Nmax), vygas(Nmax), vzgas(Nmax), temperature_gas(Nmax) )
		allocate( SNI_gas(Nmax) )
		print *, 'reading gas data'
		call data_gas(gas_file_name(k))
		deallocate( SNI_gas )
		res = minval( cell_size_gas(1:Ngas) )

		if( use_stars ) then
			print *, 'alocating star arrays'
			allocate( xstars(Nmax), ystars(Nmax), zstars(Nmax), vxstars(Nmax), vystars(Nmax), vzstars(Nmax), idstars(Nmax), mass_stars(Nmax), initial_mass_stars(Nmax), age_stars(Nmax), SNII_stars(Nmax) )
			allocate( SNI_stars(Nmax) )
			print *, 'reading star data'
			call data_stars(stars_file_name(k),0)
			call data_stars(stars_z_file_name(k),1)
			deallocate( SNI_stars )
		end if

		write(16,*) trim(snap_tag)
		write(16,*) '------'
		call half_mass_radius_spherical( Rvir )
		call center( Rhalf_spher, 0.0_4, (/ 0.0_4, 0.0_4, 1.0_4 /) )
		call disc_spin( Rhalf_spher, 0.0_4, (/ 0.0_4, 0.0_4, 1.0_4 /) )
		Rdisc2 = Rhalf_spher
		Hdisc2 = Rhalf_spher
		ang_mom2(:) = 0.0_4
		call disc_radius( 1 )
		write(16,*) ''
		write(16,*) ''
		write(16,*) ''

		call mean_disc_values( Rdisc, Hdisc, ang_mom(:), Rdisc, dm_file_name(k), stars_file_name(k), stars_z_file_name(k) )

		write(15,'(1x,a5,3(1x,i))') trim(snap_tag2), Ngas, Nstars, Ndm

		write(17,'(18es12.2)') Rdisc, Hdisc, mean_mgas, mean_mcold, mean_mstars, exsitu_stellar_mass, &
		& mean_mbar, mean_fgas, mean_sigma_gas, mean_sigma_stars, mean_sigma, mean_age, mean_metgas, &
		& mean_metstars, mean_SFR, mean_Sig_SFR, mean_SSFR, mean_tau										!!! /mean_disc_values.out

		write(18,'(1x,a5,3(1x,es12.5))') trim(snap_tag2), Rvir, Mvir, Vvir									!!! /Nir_halo_cat.txt
		write(20,'(1x,a5,1x,i12,9(1x,e12.5),2(1x,f7.3))') trim(snap_tag2), idstars(closest_star(1)), xstars(closest_star(1)), &  		!!! '_simplified_disc_cat.txt'
		& ystars(closest_star(1)), zstars(closest_star(1)), vxstars(closest_star(1)), vystars(closest_star(1)), & 
		& vzstars(closest_star(1)), ang_mom(:), Rdisc, Hdisc
		write(22,'(1x,a5,10(1x,es12.5),2(1x,f7.3),7(1x,es12.5),2(1x,f7.3))') trim(snap_tag2), rcom(:), vcom(:), ang_mom(:), L_mag, &		!!! /Nir_disc_cat.txt'
		&  Rdisc, Hdisc, mean_mgas, mean_mcold, mean_mstars, exsitu_stellar_mass, mean_mdm, mean_SFR, mean_age, mean_metgas, mean_metstars
		write(23,'(1x,a5,6(1x,es12.5),1x,f7.3,7(1x,es12.5),2(1x,f7.3))') trim(snap_tag2), rcom(:), vcom(:), Rdisc, &
		& spher_mgas, spher_mcold, spher_mstars, spher_exsitu_stellar_mass, spher_mdm, spher_SFR, spher_age, spher_metgas, spher_metstars	!!!/Nir_spherical_galaxy_cat.txt'

		write(24,'(1x,a5,2(1x,f7.3),3(1x,es12.5))') trim(snap_tag2), Rdisc, Hdisc, SFR_test_disc(1:3)						!!! disc SFR convergence

		write(25,'(1x,a5,1x,f7.3,3(1x,es12.5))') trim(snap_tag2), Rdisc, SFR_test_sphere(1:3)							!!! sphere SFR convergence

		write(26,'(1x,a5,2(1x,f7.3),10(1x,es12.5))') trim(snap_tag2), Rdisc, Hdisc, mean_mstars, initial_stellar_mass_disc, analytic_stellar_mass_disc, mean_age, &
							& exsitu_stellar_mass, initial_stellar_mass_disc_es, analytic_stellar_mass_disc_es, &
							& mean_mstars - exsitu_stellar_mass, initial_stellar_mass_disc - initial_stellar_mass_disc_es, analytic_stellar_mass_disc - analytic_stellar_mass_disc_es			!!! stellar mass loss disc

		write(27,'(1x,a5,1x,f7.3,10(1x,es12.5))') trim(snap_tag2), Rdisc, spher_mstars, initial_stellar_mass_sphere, analytic_stellar_mass_sphere, spher_age, &
							& spher_exsitu_stellar_mass, initial_stellar_mass_sphere_es, analytic_stellar_mass_sphere_es, &
							& spher_mstars - spher_exsitu_stellar_mass, initial_stellar_mass_sphere - initial_stellar_mass_sphere_es, analytic_stellar_mass_sphere - analytic_stellar_mass_sphere_es	!!! stellar mass loss sphere

		write(28,'(1x,a5,2(1x,f7.3),17(1x,es12.5))') trim(snap_tag2), Rvir, Rdisc, Mvir, Mgv, Msv, IMsv, Mdv, Mg0_1, Ms0_1, IMs0_1, Md0_1, Mg10, Ms10, IMs10, Md10, spher_mgas, spher_mstars, spher_initial_mstars, spher_mdm
		write(29,'(1x,a5,8(1x,es12.5))') trim(snap_tag2), 0.15_4*Rvir, Mgas_015Rv, Mstars_015Rv, initial_Mstars_015Rv, Mcold_gas_015Rv, Myoung_stars_015Rv, Mdm_015Rv, SFR_015Rv
		write(30,'(1x,a5,7(1x,es12.5))') trim(snap_tag2), Mgas_1kpc, Mstars_1kpc, initial_Mstars_1kpc, Mcold_gas_1kpc, Myoung_stars_1kpc, Mdm_1kpc, SFR_1kpc

		deallocate( xstars, ystars, zstars, vxstars, vystars, vzstars, idstars, mass_stars, initial_mass_stars )
		call deallocate_all()
		if(i.eq.k+1) then
			if(allocated(insitustars)) deallocate(insitustars)
			call close_files()
		end if
		print *, ''
	end do

	deallocate( gas_file_name, stars_file_name, stars_z_file_name, dm_file_name, sim )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contains
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
		if(allocated(SNI_gas)) deallocate(SNI_gas)
		if(allocated(SNII_gas)) deallocate(SNII_gas)

		if(allocated(xstars)) deallocate( xstars )
		if(allocated(ystars)) deallocate( ystars )
		if(allocated(zstars)) deallocate( zstars )
		if(allocated(vxstars)) deallocate( vxstars )
		if(allocated(vystars)) deallocate( vystars )
		if(allocated(vzstars)) deallocate( vzstars )
		if(allocated(mass_stars)) deallocate( mass_stars )
		if(allocated(initial_mass_stars)) deallocate( initial_mass_stars )
		if(allocated(age_stars)) deallocate( age_stars )
		if(allocated(idstars)) deallocate( idstars )
		if(allocated(SNI_stars)) deallocate( SNI_stars )
		if(allocated(SNII_stars)) deallocate( SNII_stars )

		if(allocated(xdm)) deallocate( xdm )
		if(allocated(ydm)) deallocate( ydm )
		if(allocated(zdm)) deallocate( zdm )
		if(allocated(vxdm)) deallocate( vxdm )
		if(allocated(vydm)) deallocate( vydm )
		if(allocated(vzdm)) deallocate( vzdm )
		if(allocated(mass_dm)) deallocate( mass_dm )
		if(allocated(iddm)) deallocate( iddm )

		if(allocated(xprime)) deallocate(xprime)
		if(allocated(yprime)) deallocate(yprime)
		if(allocated(zprime)) deallocate(zprime)
		if(allocated(vxprime)) deallocate(vxprime)
		if(allocated(vyprime)) deallocate(vyprime)
		if(allocated(vzprime)) deallocate(vzprime)
		if(allocated(rprime)) deallocate(rprime)

	end subroutine deallocate_all
!________________________________________________________________________________________

	subroutine open_files(i1,j1,k1)
		implicit none
		integer,intent(in) :: i1,j1,k1

		print *, 'enter open files2'
		write(filename,'(a,a)') 'mkdir -p ./',trim(sim(k1))
		call system(filename)
		filename = ''

		write(filename,'(a,a,a)') 'mkdir -p ./',trim(sim(k1)),'/galaxy_catalogue'
		call system(filename)
		filename = ''

		write(filename,'(a,a,a)') './',trim(sim(k1)),'/Nmax.txt'
		open(unit=15,file=filename,form='formatted')
		filename = ''
		write(15,'(1x,i3)') i1-j1-1

		write(filename,'(a,a,a)') './',trim(sim(k1)),'/disc_recursion_check.out'
		open(unit=16,file=filename,form='formatted')
		filename = ''

		write(filename,'(a,a,a)') './',trim(sim(k1)),'/mean_disc_values.out'
		open(unit=17,file=filename,form='formatted')
		filename = ''
		write(17,*) 'UNITS:'
		write(17,*) '[mass] = M_sun'
		write(17,*) '[age] = Myr'
		write(17,*) '[Rdisc / Hdisc] = kpc'
		write(17,*) '[metalicity] = log([O/H]) + 12'
		write(17,*) '[Sigma] = M_sun / pc^2'
		write(17,*) '[SFR] = M_sun / yr'
		write(17,*) '[sig_SFR] = M_sun / yr / kpc'
		write(17,*) '[SSFR] = 1 / Gyr '
		write(17,*) '[tau] = 1 / Gyr 	tau = SFR/M_gas = gas consumption time'
		write(17,*) 'gal         ','Rdisc       ','Hdisc       ','Mg          ','Mc          ','Ms          ','Es_Ms       ','Mb          ','fg          ',&
		& 'Sigg        	','Sigs        ','Sigb        ','age         ','zg          ','zs          ','SFR          ','Sig_SFR     ',& 
		& 'SSFR        ','tau         '
		write(17,*) ''

		write(filename,'(a,a,a)') './',trim(sim(k1)),'/galaxy_catalogue/Nir_halo_cat.txt'
		open(unit=18,file=filename)
		filename = ''
		write(18,'(1x,i3)') i1-j1-1

		write(filename,'(a,a,a)') './',trim(sim(k1)),'/galaxy_catalogue/Nir_simplified_disc_cat.txt'
		open(unit=20,file=filename)
		filename = ''
		write(20,'(1x,i3)') i1-j1-1

		write(filename,'(a,a,a)') './',trim(sim(k1)),'/galaxy_catalogue/Nir_disc_cat.txt'
		open(unit=22,file=filename)
		filename = ''
		write(22,'(1x,i3)') i1-j1-1

		write(filename,'(a,a,a)') './',trim(sim(k1)),'/galaxy_catalogue/Nir_spherical_galaxy_cat.txt'
		open(unit=23,file=filename)
		filename = ''
		write(23,'(1x,i3)') i1-j1-1

		write(filename,'(a,a,a)') './',trim(sim(k1)),'/galaxy_catalogue/disc_SFR_convergence.txt'
		open(unit=24,file=filename)
		filename = ''
		write(24,'(1x,i3)') i1-j1-1

		write(filename,'(a,a,a)') './',trim(sim(k1)),'/galaxy_catalogue/sphere_SFR_convergence.txt'
		open(unit=25,file=filename)
		filename = ''
		write(25,'(1x,i3)') i1-j1-1

		write(filename,'(a,a,a)') './',trim(sim(k1)),'/galaxy_catalogue/disc_mass_loss.txt'
		open(unit=26,file=filename)
		filename = ''
		write(26,'(1x,i3)') i1-j1-1

		write(filename,'(a,a,a)') './',trim(sim(k1)),'/galaxy_catalogue/sphere_mass_loss.txt'
		open(unit=27,file=filename)
		filename = ''
		write(27,'(1x,i3)') i1-j1-1

		write(filename,'(a,a,a)') './',trim(sim(k1)),'/galaxy_catalogue/Mstar.txt'
		open(unit=28,file=filename)
		filename = ''
		write(28,'(1x,i3)') i1-j1-1

		write(filename,'(a,a,a)') './',trim(sim(k1)),'/galaxy_catalogue/Nir_015_Rvir_cat.txt'
		open(unit=29,file=filename,form='formatted')
		write(29,'(1x,i3)') i1-j1-1

		write(filename,'(a,a,a)') './',trim(sim(k1)),'/galaxy_catalogue/Nir_1kpc_cat.txt'
		open(unit=30,file=filename,form='formatted')
		write(30,'(1x,i3)') i1-j1-1

		print *, 'exit open files2'
	end subroutine open_files
!________________________________________________________________________________________

	subroutine close_files()
		implicit none

		close(unit=15)
		close(unit=16)
		close(unit=17)
		close(unit=18)
		close(unit=20)
		close(unit=22)
		close(unit=23)
		close(unit=24)
		close(unit=25)
		close(unit=26)
		close(unit=27)
		close(unit=28)
		close(unit=29)
		close(unit=30)
	end subroutine close_files
end program main

