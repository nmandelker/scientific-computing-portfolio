! 2 possible changes to make:
! ---------------------------
! (1) Linked List for stellar particles instead of cell by cell search. What about DM particles?
! (2) Only cut clumps by mass in main clump_prop routine using total baryonic mass, instead of using tracer mass when identifying clumps

module parameters
!!! Parameters used in the code, which can be varied to change performance, memory usage, definitions, implimentations, etc
	implicit none

	real(8),parameter :: ngrid_max = 1.6d9				!!! Maximum number of possible cells in the uni-grid (CPU memory limited)
	integer,parameter :: Nclumps_max = 500				!!! Maximum number of clumps per snapshot

	logical,parameter :: h_equals_r = .true.      			!!! Set Hd = Rd for clump finder
	real(4),parameter :: disc_dim = 2.0_4				!!! Region searched for clumps: r' < disc_dim * Rd, |z'| < disc_dim * Hd
	real(4),parameter :: min_box_size = 5.0_4			!!! Region searched for clumps is at least +/- 5kpc in each direction.

	logical,parameter :: check_CIC = .true.				!!! Verify mass conservation of CIC interpolation during grid creation
	logical,parameter :: output_3D_smooth = .false.			!!! Output (large) binary files containing 3D smoothed maps for Matlab
	logical,parameter :: output_residual = .false.			!!! Output (large) binary files containing 3D residual maps for Matlab
	logical,parameter :: move_center_velocity = .false.		!!! Correct velocities to center of mass velocity of cold gas within disc cylinder (calculated previously)
	logical,parameter :: move_center = .false.			!!! Correct positions to center of mass of cold gas within disc cylinder (calculated previously)

	logical,parameter :: fix_res = .true.				!!! Pre-fixed size for uni-grid. Otherwise, scale byAMR resolution
	real(4),parameter :: fixed_res = 70.0_4				!!! If fix_res, then this is the resolution of the uni-grid in pc
	real(4),parameter :: res_multiplier = 2.0_4			!!! Otherwise, this is the resolution of uni-grid in units of AMR resolution

	logical,parameter :: dens_cell_thresh = .false.			!!! Use minimum (narrow) density to count cell - avoids having too many spurious cells in low density regions
	real(4),parameter :: min_dens_cell = 0.03363_4 * 0.5_4		!!! Minimum allowable density of cell to count for clump finder. 0.5cm^{-3} in M_{sun} pc^{-3}

	logical,parameter :: use_narrow = .false.			!!! Smooth initial CiC grid by a narrow Gaussian. Otherwise, just compare the CiC grid to the wide Gaussian.
	real(4),parameter :: narrow = 70.0_4				!!! Narrow gaussian FWHM in physical pc (for 2xAMR cell)

	real(4) :: wide							!!! Wide gaussian FWHM in physical pc				 --> input when running
	real(4) :: dens_thresh						!!! Minimum delta(rho)/rho threshold for detecting clumps	 --> input when running

	integer,parameter :: Ntracer = 2				!!! Number of tracers to use for clump finding

	logical,parameter :: track_clumps = .true.			!!! Track clump stars to follow their histories
	logical,parameter :: by_mass = .false.				!!! Track clump based on mass of same stars vs stellar mass of ancestor
	logical,parameter :: by_num = .true.				!!! Track clump based on number of same stars vs number of stars in ancestor (if by accident both are true, this is default)
	real(8),parameter :: same_mass_thresh = 0.5_8			!!! Threshold in same stars mass over ancestor clump stellar mass if 'by_mass'=.true.
	real(4),parameter :: same_num_thresh = 0.25_4			!!! Threshold in number of same stars over number of stars in ancestor clump if 'by_num'=.true.
	integer,parameter :: min_star_track = 10			!!! minimum number of stars in a clump to reliably track it
	real(4),parameter :: Ndyn_track = 1.e6				!!! Number of dynamical times (td_global) to look back in time when tracking clumps. (For unlimitted time, make 1.e6 --> >10^{12} years)

	real(8),parameter :: fj = 0.7_8					!!! Disc stars have jz/jmax > fj
	real(8),parameter :: dm_thresh = (10.0_8)**(0.9_8)		!!! rho_c/rho_bg for dark matter for in-situ / ex-situ clump

	real(4),parameter :: tmax_i = 20.0_4				!!! For SFR calculation, in Myr
	real(4),parameter :: tmax_f = 40.0_4				!!! For SFR calculation, in Myr
	real(4),parameter :: dtmax = 0.2_4				!!! For SFR calculation, in Myr

	real(4),parameter :: max_T = 1.5e4				!!! Maximum temperature ( in K) for "cold" gas
	integer,parameter :: vol_thresh_abs = 8				!!! Minimum number of cells in clump
	real(8),parameter :: min_mass(3) = (/ 1.d-4, 1.d-3, 1.d-2 /)	!!! Different baryonic mass thresholds for clumps, in units of disc mass
	real(8),parameter :: min_mass_abs = 1.d6			!!! Baryonic mass threshold for clumps, in solar masses
	real(8),parameter :: min_dens_abs = 0.0_4			!!! Baryonic volume density threshold for clumps, in solar mass per cubic parsec (0.03363 ~ 1cm^{-3})

	real(8),parameter :: G = 4.3d-6					!!! Units: Kpc * (km/sec)^2 * M_sun^(-1)
	real(4),parameter :: pi = 3.141592654_4, pi2 = 2.0_4*pi, pi4_3 = (4.0_4 / 3.0_4)*pi

end module parameters
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module globvar
!!! Global variables and allocatable arrays used throughout the code
	implicit none

!!!	Arrays for gas, stars and DM data
!!!	---------------------------------
	character(len=256) :: dirname					!!! output directory name
	character(len=256),allocatable :: dm_file_name(:), gas_file_name(:), stars_z_file_name(:), stars_file_name(:)

	real(4),allocatable :: xgas(:), ygas(:), zgas(:), vxgas(:), vygas(:), vzgas(:)
	real(4),allocatable :: density_gas(:), temperature_gas(:), cell_size_gas(:), SNI_gas(:), SNII_gas(:)
	integer,allocatable :: Ngas(:)

	real(8),allocatable :: xstars(:), ystars(:), zstars(:), vxstars(:), vystars(:), vzstars(:), mass_stars(:), initial_mass_stars(:)
	real(4),allocatable :: age_stars(:), SNI_stars(:), SNII_stars(:)
	integer,allocatable :: idstars(:), insitustars(:), insitustars2(:)				!!! insitustars(i) = 1 / 2 if star(i) was formed inside / outside of the disc (initialized at 0)
	integer,allocatable :: Nstars(:)

	real(8),allocatable :: xdm(:), ydm(:), zdm(:), vxdm(:), vydm(:), vzdm(:), mass_dm(:)
	integer,allocatable :: iddm(:)
	integer,allocatable :: Ndm(:)

	real(4) :: resolution, res							!!! AMR resolution and uni-grid resolution
	character(len=512),allocatable :: clump_stars_filename(:)			!!! files containing lists of stellar particles in the clumps / background in each snapshot per simulation
	character(len=512),allocatable :: final_clump_stars_filename(:,:)		!!! files containing lists of stellar particles in the clumps / background in each snapshot per simulation


!!!	Arrays for disc data to be read from input files
!!!	---------------------------------
	real(4),allocatable :: Rvir(:), aexp(:), redshift(:), rcom(:,:), vcom(:,:), Ldisc(:,:), Lmag(:), Rdisc(:), Hdisc(:), Hdisc2(:)	!!! Hdisc = disc height,   Hdisc2 = Rdisc .OR. Hdisc, depending on "h_equals_r"
	real(8),allocatable :: Mvir(:), Vvir(:), Mgas_disc(:), Mcold_disc(:), Mstar_disc(:), M_Es_star_disc(:), Mdm_disc(:), Mbar_disc(:), fgas_disc(:)
	real(8),allocatable :: sigma_gas_disc(:), sigma_stars_disc(:), sigma_bar_disc(:)
	real(8),allocatable :: SFR_disc(:),  Sig_SFR_disc(:), SSFR_disc(:), tau_disc(:)
	real(8),allocatable :: age_disc(:), metgas_disc(:), metstars_disc(:)

!!!	Arrays for tracking clump histories
!!!	---------------------------------
	integer,allocatable :: ngroup_hist(:), id_hist(:,:), clump_merger_hist(:,:), nstar_hist(:,:)
	real(8),allocatable :: clump_mass_hist(:,:)
	real(8),allocatable :: spher_mass_test(:,:)

!!!	Arrays for grid tracers
!!!	---------------------------------
	real(4),allocatable :: temp_thresh(:), age_thresh(:)
	real(4) :: box_size_r, box_size_h
	character(len=20),allocatable :: tracer_name(:)
	integer :: max_ngroup

end module globvar
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module read_binary
!!! Reads post-treated ART files of gas, stars and dark matter data
use parameters
use globvar
	implicit none
	integer :: Ngas_list, Nstars_list
	integer,allocatable :: gas_list(:), star_list(:)

contains
	subroutine allocate_gas(snap)
!!!!!!!!!! Allocate gas arrays for snapshot !!!!!!!!!!
!!!!!!!!!! Deallocated in 'deallocate_gas_stars_dm' !!!!!!!!!!
		implicit none
		integer,intent(in) :: snap
		integer :: i

		print *, 'Ngas',Ngas(snap)
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
		allocate( SNII_gas(Ngas(snap)), stat=i )  	!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation SNII_gas. stat= ', i
			stop
		end if
		allocate( SNI_gas(Ngas(snap)), stat=i )		!!! deallocated immediately !!!
		if(i.ne.0) then
			print *, 'error in allocation SNI_gas. stat= ', i
			stop
		end if
		call data_gas( gas_file_name(snap), 1, Ngas(snap) )

		resolution = minval( cell_size_gas(1:Ngas(snap)) )
		if( fix_res ) then
			print *, 'YOU SELECTED A FIXED RESOLUTION FOR THE UNI-GRID'
			print *, 'GOOD FOR YOU!'
			res = fixed_res
		else
			print *, 'YOU SELECTED TO SCALE THE UNI-GRID WITH THE AMR GRID'
			print *, 'ALRIGHT THEN!'
			if(resolution .lt. 35.0) then			!!! Keep uni-grid resolution between 35 - 70 pc
				res = res_multiplier * resolution
			else
				res = resolution
			end if
		end if

		!!! Move gas to disc rest frame, then don't have to worry about it later !!!
		if(move_center_velocity) then
			vxgas(:) = vxgas(:) - vcom(1,snap)
			vygas(:) = vygas(:) - vcom(2,snap)
			vzgas(:) = vzgas(:) - vcom(3,snap)
		end if
		if(move_center) then
			xgas(:) = xgas(:) - rcom(1,snap)
			ygas(:) = ygas(:) - rcom(2,snap)
			zgas(:) = zgas(:) - rcom(3,snap)
		end if
	end subroutine allocate_gas
!___________________________________________________________________________________________________________________________________________________________________________________________________________________
	subroutine allocate_stars(snap, typ)
!!!!!!!!!! Allocate stellar arrays for snapshot !!!!!!!!!!
!!!!!!!!!! Deallocated in 'deallocate_gas_stars_dm' !!!!!!!!!!
		implicit none
		integer,intent(in) :: snap, typ
		integer :: i, m, n

		print *, 'Nstars',Nstars(snap)
		allocate( idstars(Nstars(snap)), stat=i )	!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation idstars. stat= ', i
			stop
		end if
		allocate( xstars(Nstars(snap)), stat=i )      	!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation xstars. stat= ', i
			stop
		end if
		allocate( ystars(Nstars(snap)), stat=i )      	!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation ystars. stat= ', i
			stop
		end if
		allocate( zstars(Nstars(snap)), stat=i )      	!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation zstars. stat= ', i
			stop
		end if
		allocate( vxstars(Nstars(snap)), stat=i )     	!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation vxstars. stat= ', i
			stop
		end if
		allocate( vystars(Nstars(snap)), stat=i )     	!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation vystars. stat= ', i
			stop
		end if
		allocate( vzstars(Nstars(snap)), stat=i )     	!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation vzstars. stat= ', i
			stop
		end if
		if( typ .eq. 1 .or. typ .eq. 3 ) then
			allocate( mass_stars(Nstars(snap)), stat=i )	!!! deallocated at the end !!!
			if(i.ne.0) then
				print *, 'error in allocation mass_stars. stat= ', i
				stop
			end if
		end if
		allocate( age_stars(Nstars(snap)), stat=i )	!!! deallocated at the end !!!
		if(i.ne.0) then
			print *, 'error in allocation age_stars. stat= ', i
			stop
		end if
		if( typ .eq. 2 .or. typ .eq. 3 ) then
			allocate( initial_mass_stars(Nstars(snap)), stat=i )	!!! deallocated at the end !!!
			if(i.ne.0) then
				print *, 'error in allocation initial_mass_stars. stat= ', i
				stop
			end if
			allocate( SNII_stars(Nstars(snap)), stat=i )  	!!! deallocated at the end !!!
			if(i.ne.0) then
				print *, 'error in allocation SNII_stars. stat= ', i
				stop
			end if
			allocate( SNI_stars(Nstars(snap)), stat=i )		!!! deallocated immediately !!!
			if(i.ne.0) then
				print *, 'error in allocation SNI_stars. stat= ', i
				stop
			end if
		end if
		if( typ .eq. 1 .or. typ .eq. 3 ) call data_stars( stars_file_name(snap), 0, Nstars(snap))
		if( typ .eq. 2 .or. typ .eq. 3 ) call data_stars( stars_z_file_name(snap), 1, Nstars(snap))

		!!! Move stars to disc rest frame, then don't have to worry about it later !!!
		if(move_center_velocity) then
			vxstars(:) = vxstars(:) - real(vcom(1,snap),8)
			vystars(:) = vystars(:) - real(vcom(2,snap),8)
			vzstars(:) = vzstars(:) - real(vcom(3,snap),8)
		end if
		if(move_center) then
			xstars(:) = xstars(:) - real(rcom(1,snap),8)
			ystars(:) = ystars(:) - real(rcom(2,snap),8)
			zstars(:) = zstars(:) - real(rcom(3,snap),8)
		end if

		n = size(insitustars)
		m = maxval(idstars(1:Nstars(snap)))
		if( m .gt. n ) then
			print *, 'fixing length of insitustars array'
			print *, n, m
			allocate( insitustars2(n), stat=i )
			if(i.ne.0) then
				print *, 'error in allocation insitustars2. stat= ', i
				stop
			end if
			insitustars2(1:n) =insitustars(1:n) 
			deallocate( insitustars )

			allocate( insitustars(m), stat=i )
			if(i.ne.0) then
				print *, 'error in allocation insitustars_2. stat= ', i
				stop
			end if
			insitustars(1:n) = insitustars2(1:n)
			insitustars(n+1:m) = 0
			deallocate( insitustars2 )
		end if

	end subroutine allocate_stars
!___________________________________________________________________________________________________________________________________________________________________________________________________________________
	subroutine allocate_dm(snap)
!!!!!!!!!! Allocate dm arrays for snapshot !!!!!!!!!!
!!!!!!!!!! Deallocated in 'deallocate_gas_stars_dm' !!!!!!!!!!
		implicit none
		integer,intent(in) :: snap
		integer :: i

		allocate( xdm(Ndm(snap)), stat=i )									!!! deallocated after collecting dm data !!!
		if(i.ne.0) then
			print *, 'error in allocation xdm. stat= ', i
			stop
		end if
		allocate( ydm(Ndm(snap)), stat=i )									!!! deallocated after collecting dm data !!!
		if(i.ne.0) then
			print *, 'error in allocation ydm. stat= ', i
			stop
		end if
		allocate( zdm(Ndm(snap)), stat=i )									!!! deallocated after collecting dm data !!!
		if(i.ne.0) then
			print *, 'error in allocation zdm. stat= ', i
			stop
		end if
		allocate( mass_dm(Ndm(snap)), stat=i )								!!! deallocated after collecting dm data !!!
		if(i.ne.0) then
			print *, 'error in allocation mass_dm. stat= ', i
			stop
		end if
		allocate( vxdm(Ndm(snap)), stat=i )									!!! deallocated immediately !!!
		if(i.ne.0) then
			print *, 'error in allocation vxdm. stat= ', i
			stop
		end if
		allocate( vydm(Ndm(snap)), stat=i )									!!! deallocated immediately !!!
		if(i.ne.0) then
			print *, 'error in allocation vydm. stat= ', i
			stop
		end if
		allocate( vzdm(Ndm(snap)), stat=i )									!!! deallocated immediately !!!
		if(i.ne.0) then
			print *, 'error in allocation vzdm. stat= ', i
			stop
		end if
		allocate( iddm(Ndm(snap)), stat=i )									!!! deallocated immediately !!!
		if(i.ne.0) then
			print *, 'error in allocation iddm. stat= ', i
			stop
		end if

		call data_dm(dm_file_name(snap), Ndm(snap))

		!!! Move stars to disc rest frame, then don't have to worry about it later !!!
		if(move_center_velocity) then
			vxdm(:) = vxdm(:) - real(vcom(1,snap),8)
			vydm(:) = vydm(:) - real(vcom(2,snap),8)
			vzdm(:) = vzdm(:) - real(vcom(3,snap),8)
		end if
		if(move_center) then
			xdm(:) = xdm(:) - real(rcom(1,snap),8)
			ydm(:) = ydm(:) - real(rcom(2,snap),8)
			zdm(:) = zdm(:) - real(rcom(3,snap),8)
		end if

	end subroutine allocate_dm
!___________________________________________________________________________________________________________________________________________________________________________________________________________________
	subroutine data_gas(filename,metals,nstop)
		implicit none
		character (*),intent (in) :: filename
		integer,intent(in) :: metals, nstop
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
		if(metals.eq.1) then
			SNII_gas(:) = 1.d-10
			SNI_gas(:) = 1.d-10
		end if
		i = 1

		DO WHILE (i .le. nstop)
			if (metals.eq.1) then
				read (12,end=6) cell_size_gas(i), xgas(i), ygas(i), zgas(i), vxgas(i), vygas(i), vzgas(i), &
				& density_gas(i), temperature_gas(i), SNII_gas(i), SNI_gas(i)
			else
				read (12,end=6) cell_size_gas(i), xgas(i), ygas(i), zgas(i), vxgas(i), vygas(i), vzgas(i), &
				& density_gas(i), temperature_gas(i)
			end if
			i = i + 1
		end do
 6		continue
		close (12)
		print *, i, nstop
	end subroutine data_gas
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine data_dm(filename,nstop)
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
			read (13,end=6) iddm(i), xdm(i), ydm(i), zdm(i), vxdm(i), vydm(i), vzdm(i), mass_dm(i)
			i = i + 1
		end do
 6		continue
		close (13)
		print *, i, nstop
	end subroutine data_dm
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine data_stars(filename, metals, nstop)
		implicit none
		character (*),intent (in) :: filename
		integer,intent(in) :: metals, nstop
		integer :: i
				
		open ( 14 , file = filename, form = 'unformatted' )
		if(metals .eq. 0) then
			idstars(:) = 0
			xstars(:) = 0.0_8
			ystars(:) = 0.0_8
			zstars(:) = 0.0_8
			vxstars(:) = 0.0_8
			vystars(:) = 0.0_8
			vzstars(:) = 0.0_8
			mass_stars(:) = 0.0_8
			age_stars(:) = 0.0_4
			i = 1
			DO WHILE (i .le. nstop)
				read (14,end=6) idstars(i), xstars(i), ystars(i), zstars(i), Vxstars(i), &
				& Vystars(i), Vzstars(i), mass_stars(i), age_stars(i)
				i = i + 1
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
			i = 1
			DO WHILE (i .le. nstop)
				read (14,end=6) idstars(i), xstars(i), ystars(i), zstars(i), Vxstars(i), &
				& Vystars(i), Vzstars(i), initial_mass_stars(i), age_stars(i), SNII_stars(i), SNI_stars(i)
				i = i + 1
			end do
		end if
 6   		continue
		close (14)
		print *, i, nstop
	end subroutine data_stars
!___________________________________________________________________________________________________________________________________________________________________________________________________________________
	subroutine deallocate_gas_stars_dm()
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

		if(allocated(gas_list)) deallocate(gas_list)

		if(allocated(xstars)) deallocate( xstars )
		if(allocated(ystars)) deallocate( ystars )
		if(allocated(zstars)) deallocate( zstars )
		if(allocated(vxstars)) deallocate( vxstars )
		if(allocated(vystars)) deallocate( vystars )
		if(allocated(vzstars)) deallocate( vzstars )
		if(allocated(mass_stars)) deallocate( mass_stars )
		if(allocated(initial_mass_stars)) deallocate( initial_mass_stars )
		if(allocated(idstars)) deallocate( idstars )
		if(allocated(age_stars)) deallocate( age_stars )
		if(allocated(SNII_stars)) deallocate( SNII_stars )
		if(allocated(SNI_stars)) deallocate( SNI_stars )

		if(allocated(star_list)) deallocate(star_list)

		if(allocated(xdm)) deallocate( xdm )
		if(allocated(ydm)) deallocate( ydm )
		if(allocated(zdm)) deallocate( zdm )
		if(allocated(vxdm)) deallocate( vxdm )
		if(allocated(vydm)) deallocate( vydm )
		if(allocated(vzdm)) deallocate( vzdm )
		if(allocated(mass_dm)) deallocate( mass_dm )
		if(allocated(iddm)) deallocate( iddm )
	end subroutine deallocate_gas_stars_dm
end module read_binary
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module tools
	implicit none

contains
	subroutine jacobi(a,n,np,d,v,nrot)
		implicit none	
		integer,intent(in) :: n, np
		integer,intent(inout) :: nrot
		real(8),intent(inout) :: a(np,np), d(np), v(np,np)
		integer :: i, ip, iq, j
		real(8) :: c, gJ, h, s, sm, t, tau, thetaJ, tresh, b(n), z(n)

		do ip=1,n
			do iq=1,n
				v(ip,iq) = 0.0_8
			end do
			v(ip,ip) = 1.0_8
		end do
		do ip=1,n
			b(ip) = a(ip,ip)
			d(ip) = b(ip)
			z(ip) = 0.0_8
		end do
		nrot=0
		do i=1,50
			sm = 0.0_8
			do ip=1,n-1
				do iq=ip+1,n
					sm = sm + abs(a(ip,iq))
				end do
			end do
			if(sm.eq.0.0_8) return
			if(i.lt.4) then
				tresh=1.0d-6*sm/n**2
			else
				tresh=0.0_8
			end if
			do ip=1,n-1
				do iq=ip+1,n
					gJ=100.0_8*abs(a(ip,iq))
					if( i .gt. 4 .and. gJ .le. 1.0d-8*abs(d(ip)) .and. gJ .le. 1.0d-8*abs(d(iq)) ) then
						a(ip,iq)=0.0_8
					else if(abs(a(ip,iq)).gt.tresh) then
						h=d(iq)-d(ip)
						if(gJ .le. 1.0d-8*abs(h)) then
							t=a(ip,iq)/h
						else
							thetaJ=0.5_8*h/a(ip,iq)
							t=1.0_8/(abs(thetaJ)+sqrt(1.0_8+thetaJ**2))
							if(thetaJ.lt.0.0_8) then
								t=-t
							end if
						end if
						c=1.0_8/sqrt(1+t**2)
						s=t*c
						tau=s/(1.0_8+c)
						h=t*a(ip,iq)
						z(ip)=z(ip)-h
						z(iq)=z(iq)+h
						d(ip)=d(ip)-h
						d(iq)=d(iq)+h
						a(ip,iq)=0.0_8
						do j=1,ip-1
							gJ=a(j,ip)
							h=a(j,iq)
							a(j,ip)=gJ-s*(h+gJ*tau)
							a(j,iq)=h+s*(gJ-h*tau)
						end do
						do j=ip+1,iq-1
							gJ=a(ip,j)
							h=a(j,iq)
							a(ip,j)=gJ-s*(h+gJ*tau)
							a(j,iq)=h+s*(gJ-h*tau)
						end do
						do j=iq+1,n
							gJ=a(ip,j)
							h=a(iq,j)
							a(ip,j)=gJ-s*(h+gJ*tau)
							a(iq,j)=h+s*(gJ-h*tau)
						end do
						do j=1,n
							gJ=v(j,ip)
							h=v(j,iq)
							v(j,ip)=gJ-s*(h+gJ*tau)
							v(j,iq)=h+s*(gJ-h*tau)
						end do
						nrot=nrot+1
					endif
				end do
			end do
			do ip=1,n
				b(ip)=b(ip)+z(ip)
				d(ip)=b(ip)
				z(ip)=0.0_8
			end do
		end do
		pause 'too many iterations in jacobi'
		return
	end subroutine jacobi
end module tools
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module clump_finder
use parameters
use globvar
use read_binary
use splitter
	implicit none
!!! Grid variables
	integer :: ngrid, ngrid2, nR, nH
	real(4),allocatable :: density_grid(:,:,:), smoothed_density1(:,:,:), smoothed_density_diff(:,:,:)
	real(8) :: grid_correc, rad_correc
	integer :: Nhigh_diff1, Nhigh_diff2
	integer,allocatable :: high_diff_pos1(:), high_diff_pos2(:), grid_2_cell(:,:,:), cell_2_clump(:)
	real(4),allocatable :: high_diff_res1(:), high_diff_res2(:)
!!! Clump / Group variables
	integer :: ngroup, ncell_tot, max_counter, global_counter
	integer,allocatable :: igroup(:), jgroup(:), kgroup(:), pos_group(:), ncell(:)
	real(4),allocatable :: res_group(:), max_res(:), mean_residual(:)
contains
	subroutine grid_size(Rd, Hd, buffer, galex, snapshot)
	! calculates sizes of uni-grids and makes sure they will not be too large
	! There is a small grid of size 2Rd*2Rd*2Hd with number of cells nR*nR*nH
	! There is a larger grid which has a buffer in order to use fft and not worry about edge effects
	! Rd, Hd and buffer should be given in pc (res is also in pc)
	! 'galex' = a0.xxx
	! 'snapshot' is integer index of snapshot from main loop
		implicit none
		real(4),intent(in) :: Rd, Hd, buffer
		integer,intent(in) :: snapshot
		character(len=20),intent(in) :: galex
		real(4) :: xmax, zmax

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
		print '(a34,2(1x,f4.1),4(1x,f5.2))', 'res, del, Rdisc, xmax, Hdisc, zmax',resolution, res, Rdisc(snapshot), xmax/1000.0_4, Hdisc(snapshot), zmax/1000.0_4
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
			write(23,*) galex, grid_correc, rad_correc
			write(23,*) ''
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
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine create_grid(snapshot, temp_thresh, age_thresh, trace_num)
	! Creates uniform grids for density, with cell size = 'res_multiplier' * AMR resolution at that snapshot.
	! Gas colder than temp_thresh (in degrees kelvin) and stars younger than age_thresh (in Gyr)
	! The z axis of the grid is defined by AM.
	! 'snapshot' is integer index of snapshot from main loop
	! 'trace_num' is the index of the tracer currently used. For the first tracer, i.e. the first grid made per snapshot, we make short lists for gas cells and stellar particles within the grid
		implicit none
		integer,intent(in) :: snapshot, trace_num
		real(4),intent(in) :: temp_thresh, age_thresh
		integer :: i, j, k, m, ip(8), jp(8), kp(8), ntrack, n0, n1, n2, n3, n4
		real(4) :: split, xgrid, ygrid, zgrid, rthresh1, rthresh2
		real(4),allocatable :: rgas(:), rstars(:)
		real(8) :: test1, test2						!!! Used to test mass conservation while making the grid
		real(4) :: max_1, max_2						!!! Used to test mass conservation while making the grid
		integer :: star_count						!!! Used to test mass conservation while making the grid

		if(check_CIC) then						!!! Test mass conservation while making the grid
			test1 = 0.0_8
			test2 = 0.0_8
			max_1 = 0.0_4
			max_2 = 0.0_4
		end if
		saxis3(:) = Ldisc(:,snapshot)
		call axes(saxis1,saxis2,saxis3)

		rthresh1 = sqrt(2.0_4*(sngl(ngrid/2)**2) + sngl(ngrid2/2)**2)
		rthresh2 = sqrt(2.0_4*(sngl(nR/2)**2) + sngl(nH/2)**2)

		print *, 'allocating initial grid'
		allocate( density_grid(ngrid,ngrid,ngrid2), stat=i )					!!! deallocated in subroutine 'gaussian_smoothing' after smooth maps are made and smaller un-smoothed grid is kept !!!
		if(i.ne.0) then
			print *, 'error in allocation of density_grid. stat= ', i
			stop
		end if
		density_grid(:,:,:) = 1.d-10

		print *, 'adding gas to grid'
		print *, 'Ngas=',Ngas(snapshot)
		allocate( rgas(Ngas(snapshot)), stat=i )						!!! deallocated in this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation of rgas. stat= ', i
			stop
		end if
		rgas(:) = 1000.0_4 * sqrt(xgas(:)**2 + ygas(:)**2 + zgas(:)**2) / res

		n0 = 0
		n1 = 0
		n2 = 0
		n3 = 0
		n4 = 0
		ntrack = (Ngas(snapshot) - mod(Ngas(snapshot),10)) / 10					!!! This just helps keep track of where I am in the loop
		do i=1,Ngas(snapshot)
			if ( mod(i,ntrack) .eq. 0 .or. i .eq. 1 ) then
				print*, 'i of Ngas',i,'of',Ngas(snapshot)
			end if
			if ( rgas(i) .le. rthresh1 ) then
				if( temperature_gas(i) .le. temp_thresh ) then
					if(cell_size_gas(i)<1.5_4*res) then
						call split0( xgas(i), ygas(i), zgas(i), (/ 0.0_4, 0.0_4, 0.0_4 /) )
						n0 = n0+1
					else if(cell_size_gas(i)<3.0_4*res) then
						call split1( xgas(i), ygas(i), zgas(i), (/ 0.0_4, 0.0_4, 0.0_4 /), cell_size_gas(i))
						n1 = n1+1
					else if(cell_size_gas(i)<5.0_4*res) then
						call split2( xgas(i), ygas(i), zgas(i), (/ 0.0_4, 0.0_4, 0.0_4 /), cell_size_gas(i))
						n2 = n2+1
					else if(cell_size_gas(i)<9.0_4*res) then
						call split3( xgas(i), ygas(i), zgas(i), (/ 0.0_4, 0.0_4, 0.0_4 /), cell_size_gas(i))
						n3 = n3+1
					else
						call split4( xgas(i), ygas(i), zgas(i), (/ 0.0_4, 0.0_4, 0.0_4 /), cell_size_gas(i))
						n4 = n4+1
					end if
					k = size(xprime(:))
					split = log(sngl(k)) / log(8.0_4)
					do j=1,k
						xgrid = 1000.0_4*xprime(j)/res + ngrid/2.0_4		!!! 'xprime'=0 --> 'xgrid'=ngrid/2, 'xprime'=-xmax --> xgrid=0, 'xprime'=xmax --> xgrid=ngrid
						ygrid = 1000.0_4*yprime(j)/res + ngrid/2.0_4 
						zgrid = 1000.0_4*zprime(j)/res + ngrid2/2.0_4

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

						star_count = 0
						do m=1,8
							if(ip(m) .ge. 1 .and. ip(m) .le. ngrid .and. jp(m) .ge. 1 .and. jp(m) .le. ngrid .and. kp(m) .ge. 1 .and. kp(m) .le. ngrid2) then
								density_grid(ip(m),jp(m),kp(m)) = density_grid(ip(m),jp(m),kp(m)) + 0.03363_4*density_gas(i)*&
								& ( ( cell_size_gas(i) / (2.0_4**split) )**3 )*abs( (xgrid-ip(9-m))*(ygrid-jp(9-m))*(zgrid-kp(9-m)) )
								if(check_CIC ) star_count = star_count + 1
							end if
						end do
						if(check_CIC) then
							test1 = test1 + real(0.03363_4*density_gas(i)*( ( cell_size_gas(i) / (2.0_4**split) )**3 ),8) * (1.0_8*star_count/8.0_8)
							if(density_gas(i) .ge. max_1) then
								max_1 = density_gas(i)
							end if
						end if
					end do
					call deallocate_primes()
				end if
				if( trace_num .eq. 1 ) then
					if ( rgas(i) .le. rthresh2 ) then
						Ngas_list = Ngas_list+1
						gas_list(Ngas_list) = i
					end if
				end if
			end if
		end do
		print *, 'done adding gas to grid'
		print *, 'n0, n1, n2, n3, n4'
		print *, n0, n1, n2, n3, n4
		deallocate( rgas )
		if(check_CIC .and. temp_thresh .gt. 100.0_4) then
			test2 = sum(real(density_grid(:,:,:),8))
			max_2 = maxval(density_grid(:,:,:)) / (0.03363_4*res**3)
			print '(a17,3(1x,es11.4))', 'mass conservation',test1,test2,(test2-test1)/test1
			print '(a15,2(1x,es11.4))', 'density maximum',max_1,max_2
			if( abs((test2-test1)/test1) .gt. 0.01_8 ) then
				print *, 'SERIOUSLY?!? THATS A PRETTY BIG MASS CONSERVATION ERROR DUDE!'
				print *, (test1-test2)/test1
				stop
			end if
		end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		print *, 'adding stars to grid'
		allocate( rstars(Nstars(snapshot)), stat=i )					!!! deallocated in this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation rstars. stat= ', i
			stop
		end if
		rstars(:) = 1000.0_4 * sqrt(sngl(xstars(:)**2 + ystars(:)**2 + zstars(:)**2)) / res

		ntrack = (Nstars(snapshot) - mod(Nstars(snapshot),10)) / 10			!!! This just helps keep track of where I am in the loop
		do i=1,Nstars(snapshot)
			if ( mod(i,ntrack) .eq. 0 ) then
				print*, 'i of Nstars',i,'of',Nstars(snapshot)
			end if

			if ( rstars(i) .gt. rthresh1 ) then
				if(insitustars(idstars(i)) .eq. 0) insitustars(idstars(i)) = 2
			else
				call split0( sngl(xstars(i)), sngl(ystars(i)), sngl(zstars(i)), (/ 0.0_4, 0.0_4, 0.0_4 /) )
				if(insitustars(idstars(i)) .eq. 0) then
					if( abs(zprime(1)) .le. disc_dim*Hdisc(snapshot) .and. rprime(1) .le. disc_dim*Rdisc(snapshot) ) then
						insitustars(idstars(i)) = 1
					else
						insitustars(idstars(i)) = 2
					end if
				end if
				if(age_stars(i).le.age_thresh) then
					xgrid = 1000.0_4*xprime(1)/res + ngrid/2.0_4	!!! 'xprime'=0 --> 'xgrid'=ngrid/2, 'xprime'=-xmax --> xgrid=0, 'xprime'=xmax --> xgrid=ngrid
					ygrid = 1000.0_4*yprime(1)/res + ngrid/2.0_4 
					zgrid = 1000.0_4*zprime(1)/res + ngrid2/2.0_4

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

					star_count = 0
					do m=1,8
						if(ip(m) .ge. 1 .and. ip(m) .le. ngrid .and. jp(m) .ge. 1 .and. jp(m) .le. ngrid .and. kp(m) .ge. 1 .and. kp(m) .le. ngrid2) then
							density_grid(ip(m),jp(m),kp(m)) = density_grid(ip(m),jp(m),kp(m)) + sngl(mass_stars(i))* &
							& abs( (xgrid-ip(9-m))*(ygrid-jp(9-m))*(zgrid-kp(9-m)) )
							if(check_CIC ) star_count = star_count + 1
						end if
					end do
					if(check_CIC ) then
						test1 = test1 + mass_stars(i) * (1.0_8*star_count/8.0_8)
					end if
				end if

				call deallocate_primes()
				if( trace_num .eq. 1 ) then
					if ( rstars(i) .le. rthresh2 ) then
						Nstars_list = Nstars_list+1
						star_list(Nstars_list) = i
					end if
				end if
			end if
		end do
		print *, 'done adding stars to grid'
		deallocate( rstars )
		if(check_CIC .and. age_thresh .gt. 0.01_4) then
			test2 = sum(real(density_grid(:,:,:),8))
			print '(a24,3(1x,es11.4))', 'mass conservation part 2',test1,test2,(test2-test1)/test1
			if( abs((test2-test1)/test1) .gt. 0.01_8 ) then
				print *, 'SERIOUSLY?!? THATS A PRETTY BIG MASS CONSERVATION ERROR DUDE!'
				print *, (test1-test2)/test1
				stop
			end if
		end if

		density_grid(:,:,:) = density_grid(:,:,:) / (res**3)				!!! M_sun pc^{-3}

	end subroutine create_grid
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine gaussian_smoothing(FWHM1, FWHM2, snapshot, trace_num, gal)
! Smooths gas density grid with spherical Gaussain filters with Full-Width-Half-Maxima of FWHM1 and FWHM2.
! FWHM are given in number of cells.
! 'snapshot' is the current snapshot, used to read the disc dimensions 'Rdisc' and 'Hdisc'.
! 'trace_num' is the index of the tracer currently used, in order to place the output (surface densities, etc) in the correct file

		implicit none
		real(4),intent(in)  :: FWHM1, FWHM2
		integer,intent(in)  :: snapshot, trace_num
		character(len=20)   :: gal	!!! snap_tag2 -->  a0.***
		real(4),allocatable :: surface_density_raw(:,:), surface_density1(:,:), surface_density2(:,:), surface_density_diff(:,:)
		character(len=256)  :: filename
		integer :: i, j, k, m

		print *, 'allocating smoothed grid'
		allocate( smoothed_density1(nR,nR,nH), stat=i)				!!! deallocated in subroutine 'group_finder' after clumps are found and clump data is stored !!!
		if(i.ne.0) then
			print *, 'error in allocation smoothed_density1. stat= ', i
			stop
		end if
		if( use_narrow ) then
			print *, 'YOU SELECTED NARROW SMOOTHING'
			print *, 'SO THATS WHAT WERE DOING NOW!'
			call fft(FWHM1,1)
		else
			print *, 'YOU SELECTED NO NARROW SMOOTHING'
			print *, 'SO WERE JUST GOING TO PUT THE RAW GRID INTO THE SMOOTHED GRID'
			smoothed_density1(:,:,:) = density_grid((ngrid/2)-(nR/2)+1:(ngrid/2)+(nR/2),(ngrid/2)-(nR/2)+1:(ngrid/2)+(nR/2),(ngrid2/2)-(nH/2)+1:(ngrid2/2)+(nH/2))
		end if
		print *, 'density maximum=', maxval(density_grid(:,:,:)), maxval(smoothed_density1(:,:,:))

		print *, 'calculating raw and narrow surface density'			!!! Creates 2D face on and edge-on maps for Matlab, within (+-2Rd)*(+-2Rd)*(+-2Hd)
		allocate( surface_density_raw(nR,nR), stat=i )				!!! deallocated in this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation surface_density_raw. stat= ', i
			stop
		end if
		allocate( surface_density1(nR,nR), stat=i )				!!! deallocated in this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation surface_density1. stat= ', i
			stop
		end if
		allocate( surface_density2(nR,nR), stat=i )				!!! deallocated in this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation surface_density2. stat= ', i
			stop
		end if
		allocate( surface_density_diff(nR,nR), stat=i )				!!! deallocated in this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation surface_density_diff. stat= ', i
			stop
		end if

		surface_density1(:,:) = sum(smoothed_density1,dim=3)*res
		surface_density_raw(:,:) = sum(density_grid((ngrid/2)-(nR/2)+1:(ngrid/2)+(nR/2),(ngrid/2)-(nR/2)+1:(ngrid/2)+(nR/2),(ngrid2/2)-(nH/2)+1:(ngrid2/2)+(nH/2)),dim=3)*res

		print *, 'allocating another grid'
		allocate( smoothed_density_diff(nR,nR,nH),stat=i )                      !!! deallocated in subroutine 'group_finder' after clumps are found and clump data is stored !!!
		if(i.ne.0) then
			print *, 'error in allocation smoothed_density_diff. stat= ', i
			stop
		end if

		print *, 'I DIDNT EVEN GIVE YOU A CHOICE THIS TIME'
		print *, 'WERE DOING THE WIDE SMOOTHING'
		call fft(FWHM2,2)
		deallocate( density_grid )

		print *, 'calculating wide surface density'
		surface_density2(:,:) = sum(smoothed_density_diff,dim=3)*res

		if(output_3D_smooth) then
			print *, 'writing 3d smoothed density'
			m = 40 + Ntracer + trace_num
			write(filename,'(a,a,a,a,a,a,a)') './',trim(dirname),'/binary_grid_outputs/smoothed_density_',trim(tracer_name(trace_num)),'_',trim(gal),'.bin'
			open(unit=m,file=filename,form='unformatted')
			write(m) nR, nR, nH, box_size_r, box_size_r, box_size_h, smoothed_density1
			write(m) 2, smoothed_density_diff
			close(unit=m)
			print *, 'Ok'
		end if

		do k=1,nH
			do j=1,nR
				do i=1,nR
					smoothed_density_diff(i,j,k) = ( min(smoothed_density1(i,j,k), 1e6) / max(smoothed_density_diff(i,j,k), 1e-6) ) - 1.0_4	!!! Now, smoothed_density_diff contains the residuals !!!
					if( dens_cell_thresh ) then
						if( smoothed_density1(i,j,k) .lt. min_dens_cell .and. smoothed_density_diff(i,j,k) .ge. dens_thresh ) then
							smoothed_density_diff(i,j,k) = -1.0_4 * smoothed_density_diff(i,j,k)	!!! low density cells have negative residuals so they can be recognized later
						end if
					end if
				end do
			end do
		end do

		deallocate( smoothed_density1 )

		if(output_residual) then
		!!! Output 3D residuals to Matlab !!!
			m = 40 + 2*Ntracer + trace_num
			write(filename,'(a,a,a,a,a,a,a)') './',trim(dirname),'/binary_grid_outputs/subtracted_density_',trim(tracer_name(trace_num)),'_',trim(gal),'.bin'
			open(unit=m,file=filename,form='unformatted')
			write(m) nR, nR, nH, box_size_r, box_size_r, box_size_h, abs(smoothed_density_diff)
			close(unit=m)
		end if

		surface_density_diff(:,:) =  ( surface_density1(:,:) / surface_density2(:,:) ) - 1.0_4 
		print *, 'writing face on surface densities'
		m = 40 + trace_num
		write(m) nR, nR, real(box_size_r,8), real(box_size_r,8), Rdisc(snapshot), surface_density_raw(:,:)
		write(m) surface_density1(:,:)
		write(m) surface_density2(:,:)
		write(m) surface_density_diff(:,:)
		print *, 'Ok'

		deallocate( surface_density_raw, surface_density1, surface_density2, surface_density_diff )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SURFACE DENSITIES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		print *, 'finished!'
	end subroutine gaussian_smoothing
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine fft(Gwidth1, mat)
	use MKL_DFTI
! Performs 3-d gaussian interpolation using fft
! Gwidth1 is FWHM given in number of AMR cells.
! If mat=1 then output result into smooth_density1
! If mat=2 then output result into smooth_density_diff

		implicit none
		real(4),intent(in) :: Gwidth1
		integer,intent(in) :: mat
		complex(4),allocatable :: fft_in1d(:), gauss_1d(:)
		integer :: i, j, k, error, length(3)
		type(DFTI_DESCRIPTOR), POINTER :: fft_desc

		print *, 'entering fft'
		length = (/ ngrid,ngrid,ngrid2 /)

		print *, 'defining fft_in1d matrix'
		allocate( fft_in1d(ngrid2*ngrid**2), stat=i )				!!! deallocated in this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation fft_in1d. stat= ', i
			stop
		end if
		print *, 'defining gaussian_1d matrix'
		allocate( gauss_1d(ngrid2*ngrid**2), stat=i )				!!! deallocated in this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation gauss_1d. stat= ', i
			stop
		end if

		do k=1,ngrid2
			do j=1,ngrid
				do i=1,ngrid
					fft_in1d( (k-1)*(ngrid**2) + (j-1)*ngrid + i ) = density_grid(i,j,k)
					gauss_1d( (k-1)*(ngrid**2) + (j-1)*ngrid + i ) = 0.5_4**(( ((1.0_4*(ngrid+1))/2.0_4 - i)**2 + &
												&  ((1.0_4*(ngrid+1))/2.0_4 - j)**2 + &
												&  ((1.0_4*(ngrid2+1))/2.0_4 - k)**2 ) / ((Gwidth1/2.0_4)**2) ) 
				end do
			end do
		end do
		gauss_1d(:) = gauss_1d(:) / sum(gauss_1d(:))
		print *, 'defined 1d vectors'

		error=DftiCreateDescriptor(fft_desc, DFTI_SINGLE, DFTI_COMPLEX, 3, length )
		print *, 'error11=',error
		error=DftiCommitDescriptor( fft_desc )
		print *, 'error12=',error
		error = DftiComputeForward(fft_desc, fft_in1d)
		print *, 'error13=',error
		error = DftiFreeDescriptor(fft_desc )
		print *, 'error14=',error

		error=DftiCreateDescriptor(fft_desc, DFTI_SINGLE, DFTI_COMPLEX, 3, length )
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

		error=DftiCreateDescriptor(fft_desc, DFTI_SINGLE, DFTI_COMPLEX, 3, length )
		print *, 'error31=',error
		error=DftiCommitDescriptor( fft_desc )
		print *, 'error32=',error
		error = DftiComputeBackward(fft_desc, fft_in1d)
		print *, 'error33=',error
		error = DftiFreeDescriptor(fft_desc )
		print *, 'error34=',error

		print *, 'making output smoothed matrix'
		if(mat .eq. 1) then
			do k=1,nH/2
				do j=1,nR/2
					do i=1,nR/2
						smoothed_density1(i,j,k) = fft_in1d( (ngrid2-nH/2+k-1)*(ngrid**2) + (ngrid-nR/2+j-1)*ngrid + ngrid-nR/2+i )/(ngrid2*ngrid**2)
						smoothed_density1(nR/2+i,nR/2+j,nH/2+k) = fft_in1d( (k-1)*(ngrid**2) + (j-1)*ngrid + i )/(ngrid2*ngrid**2)
						smoothed_density1(i,j,nH/2+k) = fft_in1d( (k-1)*(ngrid**2) + (ngrid-nR/2+j-1)*ngrid + ngrid-nR/2+i )/(ngrid2*ngrid**2)
						smoothed_density1(nR/2+i,nR/2+j,k) = fft_in1d( (ngrid2-nH/2+k-1)*(ngrid**2) + (j-1)*ngrid + i )/(ngrid2*ngrid**2)
						smoothed_density1(i,nR/2+j,k) = fft_in1d( (ngrid2-nH/2+k-1)*(ngrid**2) + (j-1)*ngrid + ngrid-nR/2+i )/(ngrid2*ngrid**2)
						smoothed_density1(nR/2+i,j,nH/2+k) = fft_in1d( (k-1)*(ngrid**2) + (ngrid-nR/2+j-1)*ngrid + i )/(ngrid2*ngrid**2)
						smoothed_density1(nR/2+i,j,k) = fft_in1d( (ngrid2-nH/2+k-1)*(ngrid**2) + (ngrid-nR/2+j-1)*ngrid + i )/(ngrid2*ngrid**2)
						smoothed_density1(i,nR/2+j,nH/2+k) = fft_in1d( (k-1)*(ngrid**2) + (j-1)*ngrid + ngrid-nR/2+i )/(ngrid2*ngrid**2)
					end do
				end do
			end do
		elseif(mat .eq. 2) then
			do k=1,nH/2
				do j=1,nR/2
					do i=1,nR/2
						smoothed_density_diff(i,j,k) = fft_in1d( (ngrid2-nH/2+k-1)*(ngrid**2) + (ngrid-nR/2+j-1)*ngrid + ngrid-nR/2+i )/(ngrid2*ngrid**2)
						smoothed_density_diff(nR/2+i,nR/2+j,nH/2+k) = fft_in1d( (k-1)*(ngrid**2) + (j-1)*ngrid + i )/(ngrid2*ngrid**2)
						smoothed_density_diff(i,j,nH/2+k) = fft_in1d( (k-1)*(ngrid**2) + (ngrid-nR/2+j-1)*ngrid + ngrid-nR/2+i )/(ngrid2*ngrid**2)
						smoothed_density_diff(nR/2+i,nR/2+j,k) = fft_in1d( (ngrid2-nH/2+k-1)*(ngrid**2) + (j-1)*ngrid + i )/(ngrid2*ngrid**2)
						smoothed_density_diff(i,nR/2+j,k) = fft_in1d( (ngrid2-nH/2+k-1)*(ngrid**2) + (j-1)*ngrid + ngrid-nR/2+i )/(ngrid2*ngrid**2)
						smoothed_density_diff(nR/2+i,j,nH/2+k) = fft_in1d( (k-1)*(ngrid**2) + (ngrid-nR/2+j-1)*ngrid + i )/(ngrid2*ngrid**2)
						smoothed_density_diff(nR/2+i,j,k) = fft_in1d( (ngrid2-nH/2+k-1)*(ngrid**2) + (ngrid-nR/2+j-1)*ngrid + i )/(ngrid2*ngrid**2)
						smoothed_density_diff(i,nR/2+j,nH/2+k) = fft_in1d( (k-1)*(ngrid**2) + (j-1)*ngrid + ngrid-nR/2+i )/(ngrid2*ngrid**2)
					end do
				end do
			end do
		else
			print *, 'fft - error in definig output grid:'
			print *, 'mat should be 1 or 2'
			print *, 'mat=',mat
			stop
		end if

		deallocate(fft_in1d)
		print *, 'done with smoothing'
	end subroutine fft
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine group_finder()
! Takes 3-d residuals above threshold and finds clumps
! Removes clumps that have too few cells or that are at volume boundary

		implicit none
		real(4),allocatable :: max_res2(:), mean_residual2(:)
		integer,allocatable :: igroup2(:), jgroup2(:), kgroup2(:), ncell2(:)
		integer :: i, j, k, l, n, ntot, ngroup_good

		print*, 'finding clumps'
		ntot = 0
		ntot = count(mask = abs(smoothed_density_diff(:,:,:)) .ge. dens_thresh)
		print *, 'max residual=', maxval(abs(smoothed_density_diff(:,:,:)))
		print *, 'number of cells above resudual thresh=', ntot
		print *, 'allocating clump finding arrays'

		allocate( pos_group(ntot),stat=i )			!!! deallocated in this subroutine, after finding initial clump candidates !!!
		if(i.ne.0) then
			print *, 'error in allocation pos_group. stat= ', i
			stop
		end if
		allocate( res_group(ntot),stat=i )			!!! deallocated in this subroutine, after finding initial clump candidates !!!
		if(i.ne.0) then
			print *, 'error in allocation res_group. stat= ', i
			stop
		end if
		allocate( ncell(0:Nclumps_max),stat=i )			!!! deallocated  in this subroutine, after cutting clumps at boundary !!!
		if(i.ne.0) then
			print *, 'error in allocation ncell. stat= ', i
			stop
		end if

		pos_group(:) = 0
		res_group(:) = 0.0_4
		ncell(:) = 0
		ngroup = 0
		ncell_tot = 0
		global_counter= 0
		max_counter = 0

		print*, 'Grouping clumps together'
		n = size(ncell) - 1				!!! Because it goes from 0 to ngroup_max
		print *, 'max number of allowed groups',n
		do k=1,nH
			do j=1,nR
				do i=1,nR
					if( smoothed_density_diff(i,j,k) .ge. dens_thresh ) then	!!! In the main run through, only look at high density cells --> will not collect clumps containing only low density cells, since all their residuals are negative
						ngroup = ngroup + 1
						call collect_neigh(i, j, k, 0)

						!!! remove clumps whose volume is less than the vol_thresh_abs parameter and whose mass is less than the mass threshold !!!
						if( ncell_tot - ncell(ngroup-1) .lt. vol_thresh_abs ) then
							pos_group( ncell(ngroup-1)+1:ncell_tot ) = 0
							res_group( ncell(ngroup-1)+1:ncell_tot ) = 0.0_4
							ncell_tot = ncell(ngroup-1)
							ngroup = ngroup - 1
						else
							ncell(ngroup) = ncell_tot		!!! Index of the last cell belonging to the clump
						end if
						if( ngroup .eq. n ) then
							print *, 'Too many groups. Allocating room for another 100'

							allocate( ncell2(0:n+100),stat=l )				!!! deallocated  in this subroutine, after cutting clumps at boundary !!!
							if(l.ne.0) then
								print *, 'error in allocation ncell2. stat= ', l
								stop
							end if

							print *, 'old',n,'new',size(ncell2)-1
							ncell2(0:n) = ncell(0:n)
							deallocate(ncell)

							allocate( ncell(0:n+100),stat=l )				!!! deallocated  in this subroutine, after cutting clumps at boundary !!!
							if(l.ne.0) then
								print *, 'error in allocation ncell_2. stat= ', l
								stop
							end if

							ncell(0:n) = ncell2(0:n)
							ncell(n+1:n+100) = 0
							deallocate(ncell2)
							n = size(ncell) - 1				!!! Because it goes from 0 to ngroup_max
							print *, 'max number of allowed groups',n
						end if
					end if
				end do
			end do
		print *, k
		end do
		print *, 'ngroup=', ngroup
		print *, 'ncell_tot=', ncell_tot
		print *, 'global_counter', global_counter
		print *, 'max_counter', max_counter
		print *, ''
		deallocate( smoothed_density_diff )
		max_ngroup = max( max_ngroup, ngroup )

		if(ngroup .eq. 0) then
			print *, 'no groups were found at all'
			print *, 'ngroup=', ngroup
			deallocate( pos_group, res_group, ncell )
		else
			allocate( ncell2(0:ngroup),stat=i )
			if(i.ne.0) then
				print *, 'error in allocation ncell2. stat= ', i
				stop
			end if
			ncell2(0:ngroup) = ncell(0:ngroup)
			deallocate( ncell )

			allocate( igroup2(ncell_tot),stat=i )			!!! deallocated in this subroutine, after cutting clumps at boundary !!!
			if(i.ne.0) then
				print *, 'error in allocation igroup2. stat= ', i
				stop
			end if
			allocate( jgroup2(ncell_tot),stat=i )			!!! deallocated in this subroutine, after cutting clumps at boundary !!!
			if(i.ne.0) then
				print *, 'error in allocation jgroup2. stat= ', i
				stop
			end if
			allocate( kgroup2(ncell_tot),stat=i )			!!! deallocated in this subroutine, after cutting clumps at boundary !!!
			if(i.ne.0) then
				print *, 'error in allocation kgroup2. stat= ', i
				stop
			end if
			allocate( mean_residual2(ngroup), stat=i )		!!! deallocated in this subroutine, after cutting clumps at boundary !!!
			if(i.ne.0) then
				print *, 'error in allocation mean_residual2. stat= ', i
				stop
			end if
			allocate( max_res2(ngroup), stat=i )			!!! deallocated in this subroutine, after cutting clumps at boundary !!!
			if(i.ne.0) then
				print *, 'error in allocation max_res2. stat= ', i
				stop
			end if

			do i=1,ncell_tot
				kgroup2(i) = ceiling( sngl(pos_group(i)) / nR**2 )
				jgroup2(i) = ceiling( sngl(pos_group(i) - (nR**2) * (kgroup2(i)-1)) / nR )
				igroup2(i) = pos_group(i) - (nR**2) * (kgroup2(i)-1) - nR * (jgroup2(i)-1)
				if( pos_group(i) .ne. (kgroup2(i)-1)*(nR**2) + (jgroup2(i)-1)*nR + igroup2(i) ) then
					print *, 'YOU FUCKED UP WITH POS_GROUP!'
					stop
				end if
			end do
			deallocate(pos_group)
			do i=1,ngroup
				mean_residual2(i) = sum(res_group(ncell2(i-1)+1:ncell2(i))) / (ncell2(i)-ncell2(i-1))
				max_res2(i) = maxval(res_group(ncell2(i-1)+1:ncell2(i)))
			end do
			deallocate(res_group)

			print *, 'max ncell',maxval( ncell2(1:ngroup)-ncell2(0:ngroup-1) )
			print *, 'min ncell',minval( ncell2(1:ngroup)-ncell2(0:ngroup-1) )
			print *, 'cutting clumps by location'
			print *, 'i, good, ncells, max res, mean res, minloc, maxloc, nR, nH'
			ngroup_good = 0
			do i=1,ngroup
				n = ncell2(i)-ncell2(i-1)
				if( maxval(igroup2(ncell2(i-1)+1:ncell2(i))) .lt. nR .and. maxval(jgroup2(ncell2(i-1)+1:ncell2(i))) .lt. nR .and. maxval(kgroup2(ncell2(i-1)+1:ncell2(i))) .lt. nH .and. &
				  & minval(igroup2(ncell2(i-1)+1:ncell2(i))) .gt. 1  .and. minval(jgroup2(ncell2(i-1)+1:ncell2(i))) .gt. 1  .and. minval(kgroup2(ncell2(i-1)+1:ncell2(i))) .gt. 1 ) then
					ngroup_good = ngroup_good + 1
					if(ngroup_good .ne. i) then
						igroup2( ncell2(ngroup_good-1)+1:ncell2(ngroup_good-1)+n ) = igroup2( ncell2(i-1)+1:ncell2(i) )
						jgroup2( ncell2(ngroup_good-1)+1:ncell2(ngroup_good-1)+n ) = jgroup2( ncell2(i-1)+1:ncell2(i) )
						kgroup2( ncell2(ngroup_good-1)+1:ncell2(ngroup_good-1)+n ) = kgroup2( ncell2(i-1)+1:ncell2(i) )
						mean_residual2(ngroup_good) = mean_residual2(i)
						max_res2(ngroup_good) = max_res2(i)
						ncell2(ngroup_good) = ncell2(ngroup_good-1)+n
					end if
				end if
				print '(2(1x,i5),1x,i7,2(1x,es10.3),6(1x,i3),2(1x,i4))', i, ngroup_good, n, max_res2(i), mean_residual2(i), & 
				  & min(minval(igroup2(ncell2(i-1)+1:ncell2(i))), minval(jgroup2(ncell2(i-1)+1:ncell2(i))), minval(kgroup2(ncell2(i-1)+1:ncell2(i)))), &
				  & max(maxval(igroup2(ncell2(i-1)+1:ncell2(i))), maxval(jgroup2(ncell2(i-1)+1:ncell2(i))), maxval(kgroup2(ncell2(i-1)+1:ncell2(i)))), nR, nH
			end do
			print *, 'ngroup, ngroup_good', ngroup, ngroup_good
			print *, 'ncell_tot, ncell_tot_good', ncell_tot, ncell2(ngroup_good)
			ngroup = ngroup_good
			ncell_tot = ncell2(ngroup_good)
			print *, 'max ncell',maxval( ncell2(1:ngroup)-ncell2(0:ngroup-1) )
			print *, 'min ncell',minval( ncell2(1:ngroup)-ncell2(0:ngroup-1) )
			print *, ''

			if(ngroup.gt.0) then
				allocate( ncell(0:ngroup),stat=i )			!!! deallocated in 'deallocate_clump_prop' !!!
				if(i.ne.0) then
					print *, 'error in allocation kgroup. stat= ', i
					stop
				end if
				allocate( igroup(ncell_tot),stat=i )			!!! deallocated in 'deallocate_clump_prop' !!!
				if(i.ne.0) then
				print *, 'error in allocation igroup. stat= ', i
				stop
				end if
				allocate( jgroup(ncell_tot),stat=i )			!!! deallocated in 'deallocate_clump_prop' !!!
				if(i.ne.0) then
					print *, 'error in allocation jgroup. stat= ', i
					stop
				end if
				allocate( kgroup(ncell_tot),stat=i )			!!! deallocated in 'deallocate_clump_prop' !!!
				if(i.ne.0) then
					print *, 'error in allocation kgroup. stat= ', i
					stop
				end if
				allocate( mean_residual(ngroup), stat=i )		!!! deallocated in 'deallocate_clump_prop' !!!
				if(i.ne.0) then
					print *, 'error in allocation mean_residual. stat= ', i
					stop
				end if
				allocate( max_res(ngroup), stat=i )			!!! deallocated in 'deallocate_clump_prop' !!!
				if(i.ne.0) then
					print *, 'error in allocation max_res. stat= ', i
					stop
				end if

				ncell(0:ngroup) = ncell2(0:ngroup)
				igroup(1:ncell_tot) = igroup2(1:ncell_tot)
				jgroup(1:ncell_tot) = jgroup2(1:ncell_tot)
				kgroup(1:ncell_tot) = kgroup2(1:ncell_tot)
				mean_residual(1:ngroup) = mean_residual2(1:ngroup)
				max_res(1:ngroup) = max_res2(1:ngroup)
				deallocate( ncell2, igroup2, jgroup2, kgroup2, mean_residual2, max_res2 )

				print *, 'preparing grid that points to clumps'
				allocate( grid_2_cell(nR,nR,nH), stat=i)		!!! deallocated in 'clump_prop' after placing gas, stars and dark matter in clump cells !!!
				if(i.ne.0) then
					print *, 'error in allocation grid_2_cell. stat= ', i
					stop
				end if
				allocate( cell_2_clump(ncell_tot), stat=i)			!!! deallocated in 'clump_prop' after placing gas, stars and dark matter in clump cells !!!
				if(i.ne.0) then
					print *, 'error in allocation grid_2_cell. stat= ', i
					stop
				end if
				grid_2_cell(:,:,:) = 0
				cell_2_clump(:) = 0
				do i=1,ncell_tot
					grid_2_cell(igroup(i),jgroup(i),kgroup(i)) = i
				end do
				do i=1,ngroup
					cell_2_clump(ncell(i-1)+1:ncell(i)) = i
				end do
			else
				print *, 'you killed all the clumps :-('
				deallocate( ncell2, igroup2, jgroup2, kgroup2, mean_residual2, max_res2 )
			end if
		end if
		print *, 'finished with group finder'
		print *, ''
	end subroutine group_finder
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	recursive subroutine collect_neigh(i, j, k, counter)
		implicit none
		integer,intent(in) :: i, j, k, counter
		integer :: counter1

		counter1 = counter + 1
		global_counter = global_counter + 1
		max_counter = max(counter1, max_counter)

		ncell_tot = ncell_tot + 1
		pos_group(ncell_tot) = (k-1)*(nR*nR) + (j-1)*nR + i
		res_group(ncell_tot) = abs(smoothed_density_diff(i,j,k))
		smoothed_density_diff(i,j,k) = 0.0_4

		!!! If you've gotten this far, then you must have at least 1 high density cell. Now attach all cells, even low density ones.
		if( i-1 .ge. 1 ) then
			if( abs(smoothed_density_diff(i-1,j,k)) .ge. dens_thresh ) then
				call collect_neigh( i-1, j, k, counter1 )
			end if
		end if

		if( i+1 .le. nR ) then
			if( abs(smoothed_density_diff(i+1,j,k)) .ge. dens_thresh ) then
				call collect_neigh( i+1, j, k, counter1 )
			end if
		end if

		if( j-1 .ge. 1 ) then
			if( abs(smoothed_density_diff(i,j-1,k)) .ge. dens_thresh ) then
				call collect_neigh( i, j-1, k, counter1 )
			end if
		end if

		if( j+1 .le. nR ) then
			if( abs(smoothed_density_diff(i,j+1,k)) .ge. dens_thresh ) then
				call collect_neigh( i, j+1, k, counter1 )
			end if
		end if

		if( k-1 .ge. 1 ) then
			if( abs(smoothed_density_diff(i,j,k-1)) .ge. dens_thresh ) then
				call collect_neigh( i, j, k-1, counter1 )
			end if
		end if

		if( k+1 .le. nH ) then
			if( abs(smoothed_density_diff(i,j,k+1)) .ge. dens_thresh ) then
				call collect_neigh( i, j, k+1, counter1 )
			end if
		end if
	end subroutine collect_neigh
end module clump_finder
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module clump_properties
use parameters
use globvar
use tools
use read_binary
use splitter
use clump_finder
	implicit none
!!!!!!!!!!!!!!!!!!!!!!!!! FOR CLUMPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real(8),allocatable :: dens_group(:), vcar_clump(:,:)
	real(8),allocatable :: clump_gas_mass(:), clump_star_mass(:), clump_dm_mass(:), mexsitu(:), clump_age(:), clump_gas_met(:), clump_star_met(:), clump_SFR(:)
	integer,allocatable :: nstar_group(:,:), star_group(:,:)
	integer :: nstar_total

	real(8),allocatable :: rcar_clump(:,:), rcar_clump_disc(:,:), rcar_clump_density_peak(:,:), rcar_clump_disc_density_peak(:,:), vcyl_clump(:,:), group_rad(:,:), max_dens(:), inert_eigen(:,:), eta_group(:,:)
	real(8),allocatable :: clump_mass(:), clump_dm_dens(:), clump_gas_frac(:), clump_dm_frac(:),  dm_back(:)
	real(8),allocatable :: clump_gas_sig(:), clump_star_sig(:), clump_sig(:), clump_Sig_SFR(:), clump_SSFR(:), clump_tau(:), clump_dist(:), clump_height(:)
	real(8),allocatable :: tff(:), td(:), Mgas_in(:,:), Mgas_out(:,:), Mstars_in(:), Mstars_out(:), Mstars_formed(:)
	real(8),allocatable :: v_clump_frame(:,:), v_sq_clump_frame(:,:), sigma_clump_frame(:,:), mass_clump_frame(:), alpha_vir(:)
	real(4),allocatable :: vgas_back(:,:), vgas_sq_back(:,:), sigma_gas_back(:,:), mgas_back(:)
	integer,allocatable :: clump_id(:), exsitu(:), new(:), ncell_2d(:), group_center_dens(:,:)
	character(len=15),allocatable :: clump_comment(:)

	integer :: ngroup_total
	integer :: nbinsitu(3), nbexsitu_rot(3), nbexsitu_nonrot(3), nbbulge(3)
	real(8) :: mdisc_insitu(3), mdisc_exsitu_rot(3), mdisc_exsitu_nonrot(3), mdisc_bulge(3), SFRdisc_insitu(3), SFRdisc_exsitu_rot(3), SFRdisc_exsitu_nonrot(3), SFRdisc_bulge(3)
contains
	subroutine allocate_clump_prop1()
		implicit none
		integer :: i

		allocate( dens_group(ncell_tot), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation dens_group. stat= ', i
			stop
		end if
		allocate( vcar_clump(ngroup,3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation vcar_clump. stat= ', i
			stop
		end if
		dens_group(:) = 0.0_8
		vcar_clump(:,:) = 0.0_8

		allocate( clump_gas_mass(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_gas_mass. stat= ', i
			stop
		end if
		allocate( clump_star_mass(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_star_mass. stat= ', i
			stop
		end if
		allocate( clump_dm_mass(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_dm_mass. stat= ', i
			stop
		end if
		allocate( mexsitu(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation mexsitu. stat= ', i
			stop
		end if
		allocate( clump_age(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_age. stat= ', i
			stop
		end if
		allocate( clump_gas_met(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_gas_met. stat= ', i
			stop
		end if
		allocate( clump_star_met(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_star_met. stat= ', i
			stop
		end if
		allocate( clump_SFR(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_SFR. stat= ', i
			stop
		end if
		clump_gas_mass(:) = 0.0_8
		clump_star_mass(:) = 0.0_8
		clump_dm_mass(:) = 0.0_8
		mexsitu(:) = 0.0_8
		clump_age(:) = 0.0_8
		clump_gas_met(:) = 0.0_8
		clump_star_met(:) = 0.0_8
		clump_SFR(:) = 0.0_8

		allocate( nstar_group(ngroup,2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation nstar_group. stat= ', i
			stop
		end if
		allocate( star_group(Nstars_list,2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation star_group. stat= ', i
			stop
		end if
		nstar_group(:,:) = 0
		star_group(:,:) = 0

	end subroutine allocate_clump_prop1
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine fix_allocation_clump_prop()
		implicit none
		integer :: i
		real(8),allocatable :: dens_group2(:), vcar_clump2(:,:)
		real(8),allocatable :: clump_gas_mass2(:), clump_star_mass2(:), clump_dm_mass2(:), mexsitu2(:), clump_age2(:), clump_gas_met2(:), clump_star_met2(:), clump_SFR2(:)
		integer,allocatable :: nstar_group2(:,:)
		integer,allocatable :: igroup2(:), jgroup2(:), kgroup2(:), ncell2(:)
		real(4),allocatable :: max_res2(:), mean_residual2(:)

		allocate( igroup2(ncell_tot), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation igroup2. stat= ', i
			stop
		end if
		allocate( jgroup2(ncell_tot), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation jgroup2. stat= ', i
			stop
		end if
		allocate( kgroup2(ncell_tot), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation kgroup2. stat= ', i
			stop
		end if
		allocate( dens_group2(ncell_tot), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation dens_group2. stat= ', i
			stop
		end if

		allocate( nstar_group2(ngroup,2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation nstar_group2. stat= ', i
			stop
		end if
		allocate( ncell2(0:ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation ncell2. stat= ', i
			stop
		end if


		allocate( vcar_clump2(ngroup,3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation vcar_clump2. stat= ', i
			stop
		end if
		allocate( clump_gas_mass2(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_gas_mass2. stat= ', i
			stop
		end if
		allocate( clump_star_mass2(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_star_mass2. stat= ', i
			stop
		end if
		allocate( clump_dm_mass2(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_dm_mass2. stat= ', i
			stop
		end if
		allocate( mexsitu2(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation mexsitu2. stat= ', i
			stop
		end if
		allocate( clump_age2(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_age2. stat= ', i
			stop
		end if
		allocate( clump_gas_met2(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_gas_met2. stat= ', i
			stop
		end if
		allocate( clump_star_met2(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_star_met2. stat= ', i
			stop
		end if
		allocate( clump_SFR2(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_SFR2. stat= ', i
			stop
		end if
		allocate( max_res2(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation max_res2. stat= ', i
			stop
		end if
		allocate( mean_residual2(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation mean_residual2. stat= ', i
			stop
		end if

		dens_group2(1:ncell_tot) = dens_group(1:ncell_tot)
		igroup2(1:ncell_tot) = igroup(1:ncell_tot)
		jgroup2(1:ncell_tot) = jgroup(1:ncell_tot)
		kgroup2(1:ncell_tot) = kgroup(1:ncell_tot)

		nstar_group2(1:ngroup,:) = nstar_group(1:ngroup,:)
		ncell2(0:ngroup) = ncell(0:ngroup)

		vcar_clump2(1:ngroup,:) = vcar_clump(1:ngroup,:)
		clump_gas_mass2(1:ngroup) = clump_gas_mass(1:ngroup)
		clump_star_mass2(1:ngroup) = clump_star_mass(1:ngroup)
		clump_dm_mass2(1:ngroup) = clump_dm_mass(1:ngroup)
		mexsitu2(1:ngroup) = mexsitu(1:ngroup)
		clump_age2(1:ngroup) = clump_age(1:ngroup)
		clump_gas_met2(1:ngroup) = clump_gas_met(1:ngroup)
		clump_star_met2(1:ngroup) = clump_star_met(1:ngroup)
		clump_SFR2(1:ngroup) = clump_SFR(1:ngroup)
		max_res2(1:ngroup) = max_res(1:ngroup)
		mean_residual2(1:ngroup) = mean_residual(1:ngroup)

		deallocate( dens_group, vcar_clump, clump_gas_mass, clump_star_mass, clump_dm_mass, mexsitu, clump_age, clump_gas_met, clump_star_met, clump_SFR, nstar_group )
		deallocate( igroup, jgroup, kgroup, ncell, max_res, mean_residual )

		allocate( dens_group(ncell_tot), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation dens_group2. stat= ', i
			stop
		end if
		allocate( igroup(ncell_tot), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation igroup2. stat= ', i
			stop
		end if
		allocate( jgroup(ncell_tot), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation jgroup2. stat= ', i
			stop
		end if
		allocate( kgroup(ncell_tot), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation kgroup2. stat= ', i
			stop
		end if

		allocate( nstar_group(ngroup,2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation nstar_group2. stat= ', i
			stop
		end if
		allocate( ncell(0:ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation ncell2. stat= ', i
			stop
		end if

		allocate( vcar_clump(ngroup,3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation vcar_clump2. stat= ', i
			stop
		end if
		allocate( clump_gas_mass(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_gas_mass2. stat= ', i
			stop
		end if
		allocate( clump_star_mass(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_star_mass2. stat= ', i
			stop
		end if
		allocate( clump_dm_mass(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_dm_mass2. stat= ', i
			stop
		end if
		allocate( mexsitu(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation mexsitu2. stat= ', i
			stop
		end if
		allocate( clump_age(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_age2. stat= ', i
			stop
		end if
		allocate( clump_gas_met(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_gas_met2. stat= ', i
			stop
		end if
		allocate( clump_star_met(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_star_met2. stat= ', i
			stop
		end if
		allocate( clump_SFR(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_SFR2. stat= ', i
			stop
		end if
		allocate( max_res(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation max_res2. stat= ', i
			stop
		end if
		allocate( mean_residual(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation mean_residual2. stat= ', i
			stop
		end if

		dens_group(1:ncell_tot) = dens_group2(1:ncell_tot)
		igroup(1:ncell_tot) = igroup2(1:ncell_tot)
		jgroup(1:ncell_tot) = jgroup2(1:ncell_tot)
		kgroup(1:ncell_tot) = kgroup2(1:ncell_tot)

		nstar_group(1:ngroup,:) = nstar_group2(1:ngroup,:)
		ncell(0:ngroup) = ncell2(0:ngroup)

		vcar_clump(1:ngroup,:) = vcar_clump2(1:ngroup,:)
		clump_gas_mass(1:ngroup) = clump_gas_mass2(1:ngroup)
		clump_star_mass(1:ngroup) = clump_star_mass2(1:ngroup)
		clump_dm_mass(1:ngroup) = clump_dm_mass2(1:ngroup)
		mexsitu(1:ngroup) = mexsitu2(1:ngroup)
		clump_age(1:ngroup) = clump_age2(1:ngroup)
		clump_gas_met(1:ngroup) = clump_gas_met2(1:ngroup)
		clump_star_met(1:ngroup) = clump_star_met2(1:ngroup)
		clump_SFR(1:ngroup) = clump_SFR2(1:ngroup)
		max_res(1:ngroup) = max_res2(1:ngroup)
		mean_residual(1:ngroup) = mean_residual2(1:ngroup)

		deallocate( dens_group2, vcar_clump2, clump_gas_mass2, clump_star_mass2, clump_dm_mass2, mexsitu2, clump_age2, clump_gas_met2, clump_star_met2, clump_SFR2, nstar_group2 )
		deallocate( igroup2, jgroup2, kgroup2, ncell2, max_res2, mean_residual2 )

	end subroutine fix_allocation_clump_prop
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine allocate_clump_prop2()
		implicit none
		integer :: i

		allocate( rcar_clump(ngroup,3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation rcar_clump. stat= ', i
			stop
		end if
		allocate( rcar_clump_disc(ngroup,3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation rcar_clump_disc. stat= ', i
			stop
		end if
		allocate( rcar_clump_density_peak(ngroup,3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation rcar_clump_density_peak. stat= ', i
			stop
		end if
		allocate( rcar_clump_disc_density_peak(ngroup,3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation rcar_clump_disc_density_peak. stat= ', i
			stop
		end if
		allocate( vcyl_clump(ngroup,3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation vcyl_clump. stat= ', i
			stop
		end if
		allocate( group_rad(ngroup,4), stat=i )			!!! deallocated in 'deallocate_clump_prop' !!!
		if(i.ne.0) then
			print *, 'error in allocation group_rad. stat= ', i
			stop
		end if
		allocate( max_dens(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation max_dens. stat= ', i
			stop
		end if
		allocate( group_center_dens(ngroup,3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation group_center_dens. stat= ', i
			stop
		end if
		allocate( inert_eigen(ngroup,3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation inert_eigen. stat= ', i
			stop
		end if
		allocate( eta_group(ngroup,3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation eta_group. stat= ', i
			stop
		end if
		rcar_clump(:,:) = 0.0_8
		rcar_clump_disc(:,:) = 0.0_8
		rcar_clump_density_peak(:,:) = 0.0_8
		rcar_clump_disc_density_peak(:,:) = 0.0_8
		vcyl_clump(:,:) = 0.0_8
		group_rad(:,:) = 0.0_8
		max_dens(:) = 0.0_8
		group_center_dens(:,:) = 0
		inert_eigen(:,:) = 0.0_8
		eta_group(:,:) = 0.0_8

		allocate( clump_mass(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_mass. stat= ', i
			stop
		end if
		allocate( clump_dm_dens(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_dm_dens. stat= ', i
			stop
		end if
		allocate( clump_gas_frac(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_gas_frac. stat= ', i
			stop
		end if
		allocate( clump_dm_frac(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_dm_frac. stat= ', i
			stop
		end if
		allocate( dm_back(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation dm_back. stat= ', i
			stop
		end if
		allocate( spher_mass_test(ngroup,2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation spher_mass_test. stat= ', i
			stop
		end if
		clump_mass(:) = 0.0_8
		clump_dm_dens(:) = 0.0_8
		clump_gas_frac(:) = 0.0_8
		clump_dm_frac(:) = 0.0_8
		dm_back(:) = 0.0_8
		spher_mass_test(:,:) = 0.0_8

		allocate( clump_gas_sig(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_gas_sig. stat= ', i
			stop
		end if
		allocate( clump_star_sig(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_star_sig. stat= ', i
			stop
		end if
		allocate( clump_sig(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_sig. stat= ', i
			stop
		end if
		allocate( clump_Sig_SFR(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_Sig_SFR. stat= ', i
			stop
		end if
		allocate( clump_SSFR(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_SSFR. stat= ', i
			stop
		end if
		allocate( clump_tau(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_tau. stat= ', i
			stop
		end if
		allocate( clump_dist(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_dist. stat= ', i
			stop
		end if
		allocate( clump_height(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_height. stat= ', i
			stop
		end if
		clump_gas_sig(:) = 0.0_8
		clump_star_sig(:) = 0.0_8
		clump_sig(:) = 0.0_8
		clump_Sig_SFR(:) = 0.0_8
		clump_SSFR(:) = 0.0_8
		clump_tau(:) = 0.0_8
		clump_dist(:) = 0.0_8
		clump_height(:) = 0.0_8

		allocate( tff(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation tff. stat= ', i
			stop
		end if
		allocate( td(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation td. stat= ', i
			stop
		end if
		allocate( Mgas_in(ngroup,3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Mgas_in. stat= ', i
			stop
		end if
		allocate( Mgas_out(ngroup,10), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Mgas_out. stat= ', i
			stop
		end if
		allocate( Mstars_in(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Mstars_in. stat= ', i
			stop
		end if
		allocate( Mstars_out(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Mstars_out. stat= ', i
			stop
		end if
		allocate( Mstars_formed(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Mstars_formed. stat= ', i
			stop
		end if
		tff(:) = 0.0_8
		td(:) = 0.0_8
		Mgas_in(:,:) = 0.0_8
		Mgas_out(:,:) = 0.0_8
		Mstars_in(:) = 0.0_8
		Mstars_out(:) = 0.0_8
		Mstars_formed(:) = 0.0_8

		allocate( v_clump_frame(ngroup,3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation v_clump_frame. stat= ', i
			stop
		end if
		allocate( v_sq_clump_frame(ngroup,3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation v_sq_clump_frame. stat= ', i
			stop
		end if
		allocate( sigma_clump_frame(ngroup,3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation sigma_clump_frame. stat= ', i
			stop
		end if
		allocate( mass_clump_frame(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation mass_clump_frame. stat= ', i
			stop
		end if
		allocate( alpha_vir(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation mass_clump_frame. stat= ', i
			stop
		end if
		v_clump_frame(:,:) = 0.0_8
		v_sq_clump_frame(:,:) = 0.0_8
		sigma_clump_frame(:,:) = 0.0_8
		mass_clump_frame(:) = 0.0_8
		alpha_vir(:) = 0.0_8

		allocate( vgas_back(ngroup,3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation vgas_back. stat= ', i
			stop
		end if
		allocate( vgas_sq_back(ngroup,3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation vgas_sq_back. stat= ', i
			stop
		end if
		allocate( sigma_gas_back(ngroup,3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation sigma_gas_back. stat= ', i
			stop
		end if
		allocate( mgas_back(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation mgas_back. stat= ', i
			stop
		end if
		vgas_back(:,:) = 0.0_4
		vgas_sq_back(:,:) = 0.0_4
		sigma_gas_back(:,:) = 0.0_4
		mgas_back(:) = 0.0_4

		allocate( clump_id(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_id. stat= ', i
			stop
		end if
		allocate( exsitu(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation exsitu. stat= ', i
			stop
		end if
		allocate( new(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation new. stat= ', i
			stop
		end if
		allocate( ncell_2d(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation ncell_2d. stat= ', i
			stop
		end if
		clump_id(:) = 0
		exsitu(:) = 0
		new(:) = 0
		ncell_2d(:) = 1

		allocate( clump_comment(ngroup), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_comment. stat= ', i
			stop
		end if
		clump_comment(:) = ''

		nbinsitu(:) = 0
		nbexsitu_rot(:) = 0
		nbexsitu_nonrot(:) = 0
		nbbulge(:) = 0
		mdisc_insitu(:) = 0.0_8
		mdisc_exsitu_rot(:) = 0.0_8
		mdisc_exsitu_nonrot(:) = 0.0_8
		mdisc_bulge(:) = 0.0_8
		SFRdisc_insitu(:) = 0.0_8
		SFRdisc_exsitu_rot(:) = 0.0_8
		SFRdisc_exsitu_nonrot(:) = 0.0_8
		SFRdisc_bulge(:) = 0.0_8

	end subroutine allocate_clump_prop2
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine deallocate_clump_prop()
		implicit none
		print *, 'entering deallocate clump prop'
		if(allocated(igroup))              deallocate(igroup)
		if(allocated(jgroup))              deallocate(jgroup)
		if(allocated(kgroup))              deallocate(kgroup)
		if(allocated(ncell))               deallocate(ncell)
		if(allocated(max_res))             deallocate(max_res)
		if(allocated(mean_residual))       deallocate(mean_residual)
		if(allocated(group_rad))           deallocate(group_rad)

		if(allocated(dens_group))          deallocate(dens_group)
		if(allocated(vcar_clump))          deallocate(vcar_clump)

		if(allocated(clump_gas_mass))      deallocate(clump_gas_mass)
		if(allocated(clump_star_mass))     deallocate(clump_star_mass)
		if(allocated(clump_dm_mass))       deallocate(clump_dm_mass)
		if(allocated(mexsitu))             deallocate(mexsitu)
		if(allocated(clump_age))           deallocate(clump_age)
		if(allocated(clump_gas_met))       deallocate(clump_gas_met)
		if(allocated(clump_star_met))      deallocate(clump_star_met)
		if(allocated(clump_SFR))           deallocate(clump_SFR)

		if(allocated(nstar_group))         deallocate(nstar_group)
		if(allocated(star_group))          deallocate(star_group)

		if(allocated(rcar_clump))          deallocate(rcar_clump)
		if(allocated(rcar_clump_disc))     deallocate(rcar_clump_disc)
		if(allocated(rcar_clump_density_peak))      deallocate(rcar_clump_density_peak)
		if(allocated(rcar_clump_disc_density_peak)) deallocate(rcar_clump_disc_density_peak)
		if(allocated(vcyl_clump))          deallocate(vcyl_clump)
		if(allocated(max_dens))            deallocate(max_dens)
		if(allocated(group_center_dens))   deallocate(group_center_dens)
		if(allocated(inert_eigen))         deallocate(inert_eigen)
		if(allocated(eta_group))           deallocate(eta_group)

		if(allocated(clump_mass))          deallocate(clump_mass)
		if(allocated(clump_dm_dens))       deallocate(clump_dm_dens)
		if(allocated(clump_gas_frac))      deallocate(clump_gas_frac)
		if(allocated(clump_dm_frac))       deallocate(clump_dm_frac)
		if(allocated(dm_back))             deallocate(dm_back)
		if(allocated(spher_mass_test))     deallocate(spher_mass_test)

		if(allocated(clump_gas_sig))       deallocate(clump_gas_sig)
		if(allocated(clump_star_sig))      deallocate(clump_star_sig)
		if(allocated(clump_sig))           deallocate(clump_sig)
		if(allocated(clump_Sig_SFR))       deallocate(clump_Sig_SFR)
		if(allocated(clump_ssfr))          deallocate(clump_ssfr)
		if(allocated(clump_tau))           deallocate(clump_tau)
		if(allocated(clump_dist))          deallocate(clump_dist)
		if(allocated(clump_height))        deallocate(clump_height)

		if(allocated(tff))                 deallocate(tff)
		if(allocated(td))                  deallocate(td)
		if(allocated(Mgas_in))             deallocate(Mgas_in)
		if(allocated(Mgas_out))            deallocate(Mgas_out)
		if(allocated(Mstars_in))           deallocate(Mstars_in)
		if(allocated(Mstars_out))          deallocate(Mstars_out)
		if(allocated(Mstars_formed))       deallocate(Mstars_formed)

		if(allocated(v_clump_frame))       deallocate(v_clump_frame)
		if(allocated(v_sq_clump_frame))    deallocate(v_sq_clump_frame)
		if(allocated(sigma_clump_frame))   deallocate(sigma_clump_frame)
		if(allocated(mass_clump_frame))    deallocate(mass_clump_frame)
		if(allocated(alpha_vir))           deallocate(alpha_vir)

		if(allocated(vgas_back))           deallocate(vgas_back)
		if(allocated(vgas_sq_back))        deallocate(vgas_sq_back)
		if(allocated(sigma_gas_back))      deallocate(sigma_gas_back)
		if(allocated(mgas_back))           deallocate(mgas_back)

		if(allocated(clump_id))            deallocate(clump_id)
		if(allocated(exsitu))              deallocate(exsitu)
		if(allocated(new))                 deallocate(new)
		if(allocated(ncell_2d))            deallocate(ncell_2d)

		if(allocated(clump_comment))       deallocate(clump_comment)
		print *, 'exiting deallocate clump prop'
	end subroutine deallocate_clump_prop
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine clump_prop(gal, snapshot)
		implicit none
		character(len=20),intent(in) :: gal	!!! snap_tag2 -->  a0.***
		integer,intent(in) :: snapshot		!!! k --> snapshot number
		character(len=4) :: comm
		character(len=15) :: temporary_comment
		character(len=512) :: temp_name, filename
		integer,allocatable :: vec_loc(:), closest_star1(:,:)
		integer,allocatable :: IND_star(:), ID_star(:), IS_star(:)
		integer,allocatable :: IND_star2(:), ID_star2(:), IS_star2(:)
		integer :: i, j, k, l, l1, m, n, n1, n2
		integer :: temp_loc(1), ntrack, ngroup_good, nrot, star_in_clump(4), ntmax
		integer :: nstar_clump, nstar_clump2, same_stars_num, for_merger(4), bad_clump1, bad_clump2, bad_clump_id, bulge_ind(2)
		integer :: ip(8), jp(8), kp(8)
		real(4) :: split, rthresh, delt, vgas_back_global, rgas_back_global, mgas_back_global, td_global, vtemp, Ltmax
		real(4) :: xgrid, ygrid, zgrid
		real(4),allocatable :: tmax(:)
		real(8),allocatable :: temp_vec(:), mass_bins(:), rad_bins(:), SFR_tmax(:,:)
		real(8),allocatable :: xstar(:),  ystar(:),  zstar(:),  vxstar(:),  vystar(:),  vzstar(:),  mass_star(:),  initial_mass_star(:),  age_star(:),  met_star(:)
		real(8),allocatable :: xstar2(:), ystar2(:), zstar2(:), vxstar2(:), vystar2(:), vzstar2(:), mass_star2(:), initial_mass_star2(:), age_star2(:), met_star2(:)  
		real(8) :: inert(3,3), inert2(3,3), v(3,3), temp_var, temp_var2,  aI, bI, cI, same_stars_mass, mstar_in_clump(4), Vesc1, Vesc2, Vesc3, mass_temporary, rthresh8
		logical :: go_on

		print *, 'allocating initial set of clump properties'
		call allocate_clump_prop1()
		saxis3(:) = Ldisc(:,snapshot)
		call axes(saxis1,saxis2,saxis3)
		write(22,'(3(1x,e12.5))') saxis1
		write(22,'(3(1x,e12.5))') saxis2
		write(22,'(3(1x,e12.5))') saxis3
		write(22,'(2(1x,f7.3))') Rdisc(snapshot),Hdisc(snapshot)

		Ltmax = tmax_f - tmax_i
		ntmax = floor(Ltmax / dtmax) + 1
		allocate( SFR_tmax(ngroup,ntmax), tmax(ntmax) )
		tmax(1:ntmax) = (tmax_i - dtmax) + (/ (i,i=1,ntmax) /)*dtmax	!!! Myr
		tmax(:) = tmax(:) / 1000.0_4	!!! Gyr
		SFR_tmax(:,:) = 0.0_8

		!!! clump gas properties !!!
		print *, 'clump gas properties'
		ntrack = (Ngas_list - mod(Ngas_list,20)) / 20						!!! This just helps keep track of where I am in the loop
		do i=1,Ngas_list
			if( i .eq. 1 ) print *, 'i of Ngas_list',i,'of',Ngas_list
			if( mod(i,ntrack) .eq. 0 ) then
				print*, 'i of Ngas_list',i,'of',Ngas_list
			end if
			m = gas_list(i)
			if(cell_size_gas(m) < 1.5_4*res) then
				call split0( xgas(m), ygas(m), zgas(m), (/ 0.0_4, 0.0_4, 0.0_4 /) )
			else if(cell_size_gas(m) < 3.0_4*res) then
				call split1( xgas(m), ygas(m), zgas(m), (/ 0.0_4, 0.0_4, 0.0_4 /), cell_size_gas(m) )
			else if(cell_size_gas(m) < 5.0_4*res) then
				call split2( xgas(m), ygas(m), zgas(m), (/ 0.0_4, 0.0_4, 0.0_4 /), cell_size_gas(m) )
			else if(cell_size_gas(m) < 9.0_4*res) then
				call split3( xgas(m), ygas(m), zgas(m), (/ 0.0_4, 0.0_4, 0.0_4 /), cell_size_gas(m) )
			else
				call split4( xgas(m), ygas(m), zgas(m), (/ 0.0_4, 0.0_4, 0.0_4 /), cell_size_gas(m) )
			end if
			l1 = size(xprime(:))
			split = log(1.0_4*l1)/log(8.0_4)
			do l=1,l1
				if(abs(xprime(l)) .le. box_size_r .and. abs(yprime(l)) .le. box_size_r .and. abs(zprime(l)) .le. box_size_h) then
					xgrid = 1000.0_4*xprime(l)/res + nR/2.0_4		!!! 'xprime'=0 --> 'xgrid'=ngrid/2, 'xprime'=-xmax --> xgrid=0, 'xprime'=xmax --> xgrid=ngrid
					ygrid = 1000.0_4*yprime(l)/res + nR/2.0_4 
					zgrid = 1000.0_4*zprime(l)/res + nH/2.0_4

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
					do j=1,8
						if(ip(j) .ge. 1 .and. ip(j) .le. nR .and. jp(j) .ge. 1 .and. jp(j) .le. nR .and. kp(j) .ge. 1 .and. kp(j) .le. nH) then
							n = grid_2_cell(ip(j), jp(j), kp(j))
							if(n .ne. 0) then
								n1 = cell_2_clump(n)
								dens_group(n) = dens_group(n) + 0.03363_8*real( density_gas(m) * ( (cell_size_gas(m) / (2.0_4**split))**3 ) * &
									& abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8)		!!! M_{sun}
								vcar_clump(n1,1)   = vcar_clump(n1,1)   + 0.03363_8*real( density_gas(m) * ( (cell_size_gas(m) / (2.0_4**split))**3 ) * &
									& abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ) * vxgas(m),8)	!!! M_{sun} * km/s
								vcar_clump(n1,2)   = vcar_clump(n1,2)   + 0.03363_8*real( density_gas(m) * ( (cell_size_gas(m) / (2.0_4**split))**3 ) * &
									& abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ) * vygas(m),8)	!!! M_{sun} * km/s
								vcar_clump(n1,3)   = vcar_clump(n1,3)   + 0.03363_8*real( density_gas(m) * ( (cell_size_gas(m) / (2.0_4**split))**3 ) * &
									& abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ) * vzgas(m),8)	!!! M_{sun} * km/s
								clump_gas_met(n1)  = clump_gas_met(n1)  + 0.03363_8*real( density_gas(m) * ( (cell_size_gas(m) / (2.0_4**split))**3 ) * &
									& abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8) * (1.d12)*real(SNII_gas(m),8) / (0.755_8*16.0_8*2.0_8)	!!! M_{sun} * 10^{12}*O/H
							end if
						end if
					end do
				end if
			end do
			call deallocate_primes()
		end do
		do i=1,ngroup
			clump_gas_mass(i) = sum( dens_group(ncell(i-1)+1:ncell(i)) )			!!! M_{sun}
			if( clump_gas_mass(i) .gt. 0.0_8 ) then
				clump_gas_met(i) = log10( clump_gas_met(i) / clump_gas_mass(i) )	!!! Log(O/H)+12
			else
				clump_gas_met(i) = 0.0_8						!!! Log(O/H)+12
			end if
		end do

		!!! clump star properties !!!
		print *, 'clump star properties'
		nstar_total = 0
		ntrack = (Nstars_list - mod(Nstars_list,20)) / 20					!!! This just helps keep track of where I am in the loop
		do i=1,Nstars_list
			star_in_clump(:) = 0
			mstar_in_clump(:) = 0.0_8
			if( i .eq. 1 ) print *, 'i of Nstars_list',i,'of',Nstars_list
			if( mod(i,ntrack) .eq. 0 ) then
				print*, 'i of Nstars_list',i,'of',Nstars_list
			end if
			m = star_list(i)
			call split0( sngl(xstars(m)), sngl(ystars(m)), sngl(zstars(m)), (/ 0.0_4, 0.0_4, 0.0_4 /) )
			if(abs(xprime(1)) .le. box_size_r .and. abs(yprime(1)) .le. box_size_r .and. abs(zprime(1)) .le. box_size_h) then
				!!! 'xprime'=0 --> 'xgrid' = nR/2, 'xprime' = - disc_dim * Rdisc --> xgrid = 0, 'xprime' = disc_dim * Rdisc --> xgrid = nR
				xgrid = 1000.0_4*xprime(1)/res + nR/2.0_4
				ygrid = 1000.0_4*yprime(1)/res + nR/2.0_4 
				zgrid = 1000.0_4*zprime(1)/res + nH/2.0_4

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
				do j=1,8
					if(ip(j) .ge. 1 .and. ip(j) .le. nR .and. jp(j) .ge. 1 .and. jp(j) .le. nR .and. kp(j) .ge. 1 .and. kp(j) .le. nH) then
						n = grid_2_cell(ip(j), jp(j), kp(j))
						if(n .ne. 0) then
							n1 = cell_2_clump(n)
							dens_group(n)       = dens_group(n)        + mass_stars(m) * real( abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8)
							clump_star_mass(n1) = clump_star_mass(n1)  + mass_stars(m) * real( abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8)
							vcar_clump(n1,1)    = vcar_clump(n1,1)     + mass_stars(m) * real( abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8) * vxstars(m)
							vcar_clump(n1,2)    = vcar_clump(n1,2)     + mass_stars(m) * real( abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8) * vystars(m)
							vcar_clump(n1,3)    = vcar_clump(n1,3)     + mass_stars(m) * real( abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8) * vzstars(m)
							clump_age(n1)       = clump_age(n1)        + mass_stars(m) * real( abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8) * 1000.0_8 * real(age_stars(m), 8)
							clump_star_met(n1)  = clump_star_met(n1)   + mass_stars(m) * real( abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8) * &
								& (1.d12)*real(SNII_stars(m), 8) / (0.755_8*16.0_8*2.0_8)
							if( insitustars(idstars(m)) .ne. 1 ) then
								mexsitu(n1) = mexsitu(n1)          + mass_stars(m) * real( abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8)
							end if
							if( age_stars(m) .le. tmax(ntmax) ) then
								n2 = 1
								do while(n2 .le. ntmax)
									if( age_stars(m) .le. tmax(n2) ) then
										SFR_tmax(n1,n2:ntmax) = SFR_tmax(n1,n2:ntmax) + initial_mass_stars(m) * real( abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8)
										n2 = ntmax + 10
									else
										n2 = n2 + 1
									end if
								end do
							end if
			!!! Place star in clump for list of stellar particles within the clump. Possible that 2 clumps get part of the star's mass due to CiC. Only list the star if the clump has more than half the particle mass
							if( star_in_clump(1) .eq. 0 ) then
								star_in_clump(1) = n1
								mstar_in_clump(1) = mstar_in_clump(1) + real( abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8)
							elseif( star_in_clump(1) .eq. n1 ) then
								mstar_in_clump(1) = mstar_in_clump(1) + real( abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8)
							else
								if( star_in_clump(2) .eq. 0 ) then
									star_in_clump(2) = n1
									mstar_in_clump(2) = mstar_in_clump(2) + real( abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8)
								elseif( star_in_clump(2) .eq. n1 ) then		!!! Possible to have 2 clumps getting 2 cells each
									mstar_in_clump(2) = mstar_in_clump(2) + real( abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8)
								else
									if( star_in_clump(3) .eq. 0 ) then
										star_in_clump(3) = n1
										mstar_in_clump(3) = mstar_in_clump(3) + real( abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8)
									elseif( star_in_clump(3) .eq. n1 ) then		!!! Possible to have 3 clumps getting 1 cell each
										mstar_in_clump(3) = mstar_in_clump(3) + real( abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8)
									else
										if( star_in_clump(4) .eq. 0 ) then
											star_in_clump(4) = n1
											mstar_in_clump(4) = mstar_in_clump(4) + real( abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8)
										elseif( star_in_clump(4) .eq. n1 ) then		!!! Possible to have 4 clumps getting 1 cell each
											mstar_in_clump(4) = mstar_in_clump(4) + real( abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8)
										else
											print *, 'something is wrong placing the stars in the clumps'
											print *, 'It looks like the same star connects to 5 DIFFERENT clumps'
											print *, star_in_clump(1),star_in_clump(2),star_in_clump(3),star_in_clump(4),n1
											stop
										end if
									end if
								end if
							end if
						end if
					end if
				end do
				temp_loc(:) = maxloc( mstar_in_clump(1:4) )
				if( mstar_in_clump(temp_loc(1)) .ge. 0.5_8 ) then
					n1 = star_in_clump(temp_loc(1))
					nstar_total = nstar_total + 1

					nstar_group(n1,1) = nstar_group(n1,1) + 1	!!! Counts the number of stars in the group
					star_group(nstar_total,1) = m			!!! Lists all stars that are in any group, by their index in the main stellar arrays

					star_group(nstar_total,2) = nstar_group(n1,2)	!!! Pointing to the previous star that was added to the group, or 0 if there were no stars added yet
					nstar_group(n1,2) = nstar_total			!!! Index pointing to the latest star to be added to the group
				end if
			!!! The idea is that at the end, nstar_group(n1,2) will point to the last star added to the group. Then the list will take you back to all previous stars in the group, until you reach 0
			end if
			call deallocate_primes()
		end do
		do i=1,ngroup
			if( clump_star_mass(i) .gt. 0.0_8 ) then
				clump_star_met(i) = log10( clump_star_met(i) / clump_star_mass(i) )		!!! Log(O/H)+12
				clump_age(i) = clump_age(i) / clump_star_mass(i)				!!! Myr
			else
				clump_star_met(i) = 0.0_8 							!!! Log(O/H)+12
				clump_age(i) = 0.0_8 								!!! Myr
			end if
			if( clump_gas_mass(i) + clump_star_mass(i) .gt. 0.0_8 ) then
				vcar_clump(i,:) = vcar_clump(i,:) / (clump_gas_mass(i) + clump_star_mass(i))	!!! center of mass velocity, box frame, km/s
			else
				vcar_clump(i,:) = 0.0_8
			end if
			if( abs(sum(dens_group(ncell(i-1)+1:ncell(i))) - (clump_gas_mass(i) + clump_star_mass(i))) / sum(dens_group(ncell(i-1)+1:ncell(i))) .gt. 0.001_8 ) then
				print *, 'something is wrong with the clump masses'
				print *, sum(dens_group(ncell(i-1)+1:ncell(i))), (clump_gas_mass(i) + clump_star_mass(i))
				stop
			end if
			SFR_tmax(i,1:ntmax) = SFR_tmax(i,1:ntmax) / ( 1.d9 * real(tmax(1:ntmax),8) )
			clump_SFR(i) = sum(SFR_tmax(i,1:ntmax)) / (1.0_8 * ntmax)
		end do
		deallocate( SFR_tmax, tmax )
		!!! clump dark matter properties !!!
		print *, 'clump dark matter properties'
		call allocate_dm(snapshot)
		deallocate( vxdm, vydm, vzdm, iddm )

		ntrack = (Ndm(snapshot) - mod(Ndm(snapshot),20)) / 20								!!! This just helps keep track of where I am in the loop
		do i=1,Ndm(snapshot)
			if( i .eq. 1 ) print*, 'i of Ndm',i,'of',Ndm(snapshot)
			if( mod(i,ntrack) .eq. 0 ) then
				print*, 'i of Ndm',i,'of',Ndm(snapshot)
			end if
			if( xdm(i)**2 + ydm(i)**2 + zdm(i)**2 .le. real(2.0_4 * (box_size_r)**2 + (box_size_h)**2, 8) ) then
				call split0( sngl(xdm(i)), sngl(ydm(i)), sngl(zdm(i)), (/ 0.0_4, 0.0_4, 0.0_4 /) )
				if(abs(xprime(1)) .le. box_size_r .and. abs(yprime(1)) .le. box_size_r .and. abs(zprime(1)) .le. box_size_h) then
					!!! 'xprime'=0 --> 'xgrid' = nR/2, 'xprime' = - disc_dim * Rdisc --> xgrid = 0, 'xprime' = disc_dim * Rdisc --> xgrid = nR
					xgrid = 1000.0_4*xprime(1)/res + nR/2.0_4
					ygrid = 1000.0_4*yprime(1)/res + nR/2.0_4 
					zgrid = 1000.0_4*zprime(1)/res + nH/2.0_4

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
					do j=1,8
						if(ip(j) .ge. 1 .and. ip(j) .le. nR .and. jp(j) .ge. 1 .and. jp(j) .le. nR .and. kp(j) .ge. 1 .and. kp(j) .le. nH) then
							n = grid_2_cell(ip(j), jp(j), kp(j))
							if(n .ne. 0) then
								n1 = cell_2_clump(n)
								clump_dm_mass(n1) = clump_dm_mass(n1) + mass_dm(i) * real( abs( (xgrid-ip(9-j))*(ygrid-jp(9-j))*(zgrid-kp(9-j)) ), 8)
							end if
						end if
					end do
				end if
			end if
		end do
		bulge_ind(:) = 0
		do k=0,1
			do j=0,1
				do i=0,1
					n = grid_2_cell(nR/2+i, nR/2+j, nH/2+k)
					if(n .ne. 0) then
						if(bulge_ind(1) .eq. 0) then
							bulge_ind(1:2) = cell_2_clump(n)
						elseif(bulge_ind(1) .ne. cell_2_clump(n)) then
							bulge_ind(2) = cell_2_clump(n)
						end if
					end if
				end do
			end do
		end do
		deallocate( grid_2_cell, cell_2_clump )
		if( bulge_ind(1) .ne. bulge_ind(2) ) then	!!! Either they're both 0, so they're equal; or one is 0 and one is a number, and the number must be (1); or they're both numbers
			if( clump_gas_mass(bulge_ind(2)) + clump_star_mass(bulge_ind(2)) .gt. clump_gas_mass(bulge_ind(1)) + clump_star_mass(bulge_ind(1)) ) then
				bulge_ind(1) = bulge_ind(2)
			end if
			n = ncell(bulge_ind(1)) - ncell(bulge_ind(1)-1)
			if( clump_gas_mass(bulge_ind(1)) + clump_star_mass(bulge_ind(1)) .lt. min_mass_abs .or. &
			& ( clump_gas_mass(bulge_ind(1)) + clump_star_mass(bulge_ind(1)) ) / real(n * (res**(3.0_4)), 8) .lt. min_dens_abs ) then
				bulge_ind(1) = 0
			end if
		end if

		!!! cut clumps by mass !!!
		print *, 'cutting clumps by baryonic mass'
		if( bulge_ind(1) .gt. 0 ) then
			print *, 'bulge_ind',bulge_ind(1), clump_gas_mass(bulge_ind(1)) + clump_star_mass(bulge_ind(1))
		else
			print *, 'bulge_ind',bulge_ind(1)
		end if
		print *, 'min, max Mbar', minval( clump_gas_mass(1:ngroup) + clump_star_mass(1:ngroup)), maxval(clump_gas_mass(1:ngroup) + clump_star_mass(1:ngroup))
		print *, 'min, max normalized Mbar', minval( (clump_gas_mass(1:ngroup) + clump_star_mass(1:ngroup)) / Mbar_disc(snapshot) ), maxval( (clump_gas_mass(1:ngroup) + clump_star_mass(1:ngroup)) / Mbar_disc(snapshot) )
		print *, 'min, max Ncell', minval( ncell(1:ngroup)-ncell(0:ngroup-1) ), maxval( ncell(1:ngroup)-ncell(0:ngroup-1) )
		print *, 'min, max Nstars', minval( nstar_group(1:ngroup,1) ), maxval(nstar_group(1:ngroup,1) )
		print *, 'min, max density', minval( (clump_gas_mass(1:ngroup) + clump_star_mass(1:ngroup)) / real((ncell(1:ngroup)-ncell(0:ngroup-1)) * (res**(3.0_4)), 8) / 0.03363_8 ), &
					&    maxval( (clump_gas_mass(1:ngroup) + clump_star_mass(1:ngroup)) / real((ncell(1:ngroup)-ncell(0:ngroup-1)) * (res**(3.0_4)), 8) / 0.03363_8 )
		print *, '-----------------------------------------------------------'
		print *, 'i, good, ncells, nstars, mclump, normalized mclump, dark matter mass, baryon density'
		ngroup_good = 0
		do i=1,ngroup
			n = ncell(i)-ncell(i-1)
			if( clump_gas_mass(i) + clump_star_mass(i) .ge. min_mass_abs .and. & 
			& ( clump_gas_mass(i) + clump_star_mass(i) ) / real(n * (res**(3.0_4)), 8) .ge. min_dens_abs ) then
				ngroup_good = ngroup_good + 1
				if(ngroup_good .ne. i) then
					igroup( ncell(ngroup_good-1)+1:ncell(ngroup_good-1)+n ) = igroup( ncell(i-1)+1:ncell(i) )
					jgroup( ncell(ngroup_good-1)+1:ncell(ngroup_good-1)+n ) = jgroup( ncell(i-1)+1:ncell(i) )
					kgroup( ncell(ngroup_good-1)+1:ncell(ngroup_good-1)+n ) = kgroup( ncell(i-1)+1:ncell(i) )
					dens_group( ncell(ngroup_good-1)+1:ncell(ngroup_good-1)+n ) = dens_group( ncell(i-1)+1:ncell(i) )
					ncell(ngroup_good) = ncell(ngroup_good-1)+n

					mean_residual(ngroup_good) = mean_residual(i)
					max_res(ngroup_good) = max_res(i)
					vcar_clump(ngroup_good,1:3) = vcar_clump(i,1:3)
					clump_gas_mass(ngroup_good) = clump_gas_mass(i)
					clump_star_mass(ngroup_good) = clump_star_mass(i)
					clump_dm_mass(ngroup_good) = clump_dm_mass(i)
					mexsitu(ngroup_good) = mexsitu(i)
					clump_age(ngroup_good) = clump_age(i)
					clump_gas_met(ngroup_good) = clump_gas_met(i)
					clump_star_met(ngroup_good) = clump_star_met(i)
					clump_SFR(ngroup_good) = clump_SFR(i)

					nstar_group(ngroup_good,1:2) = nstar_group(i,1:2)

					if( bulge_ind(1) .eq. i ) then
						bulge_ind(1) = ngroup_good
					end if
				end if
			end if
			print '(2(1x,i5),2(1x,i7),4(1x,es10.3))', i, ngroup_good, n, nstar_group(i,1), clump_gas_mass(i)+clump_star_mass(i), (clump_gas_mass(i)+clump_star_mass(i))/Mbar_disc(snapshot), clump_dm_mass(i), &
			& ( clump_gas_mass(i) + clump_star_mass(i) ) / real(n * (res**(3.0_4)), 8) / 0.03363_8
		end do
		print *, '-----------------------------------------------------------'
		print *, 'ngroup, ngroup_good', ngroup, ngroup_good
		print *, 'ncell_tot, ncell_tot_good', ncell_tot, ncell(ngroup_good)
		ngroup = ngroup_good
		ncell_tot = ncell(ngroup_good)
		if( bulge_ind(1) .gt. 0 ) then
			print *, 'bulge_ind',bulge_ind(1), clump_gas_mass(bulge_ind(1)) + clump_star_mass(bulge_ind(1))
		else
			print *, 'bulge_ind',bulge_ind(1)
		end if

		if(ngroup.gt.0) then
			print *, 'min, max Mbar', minval(clump_gas_mass(1:ngroup) + clump_star_mass(1:ngroup)), maxval(clump_gas_mass(1:ngroup) + clump_star_mass(1:ngroup))
			print *, 'min, max normalized Mbar', minval( (clump_gas_mass(1:ngroup) + clump_star_mass(1:ngroup)) / Mbar_disc(snapshot) ), maxval( (clump_gas_mass(1:ngroup) + clump_star_mass(1:ngroup)) / Mbar_disc(snapshot) )
			print *, 'min, max Ncell', minval( ncell(1:ngroup)-ncell(0:ngroup-1) ), maxval( ncell(1:ngroup)-ncell(0:ngroup-1) )
			print *, 'min, max Nstars', minval( nstar_group(1:ngroup,1) ), maxval( nstar_group(1:ngroup,1) )
			print *, 'min, max density', minval( (clump_gas_mass(1:ngroup) + clump_star_mass(1:ngroup)) / real((ncell(1:ngroup)-ncell(0:ngroup-1)) * (res**(3.0_4)), 8) / 0.03363_8 ), &
						&    maxval( (clump_gas_mass(1:ngroup) + clump_star_mass(1:ngroup)) / real((ncell(1:ngroup)-ncell(0:ngroup-1)) * (res**(3.0_4)), 8) / 0.03363_8 )
			print *, ''
			print *, 'writing initial stellar particle files'
			allocate( clump_stars_filename(ngroup), stat=i )         !!! deallocated at end of snapshot loop !!!
			if(i.ne.0) then
				print *, 'error in allocation clump_stars_filename. stat= ', i
				stop
			end if
			do i=1,ngroup
				write(clump_stars_filename(i),'(a,a,a,a,a,i5.5,a)') './',trim(dirname),'/clump_stars/',trim(gal),'_initial_clump',i,'.out'
				open(unit=81,file=clump_stars_filename(i))
				write(81,'(i10)') nstar_group(i,1)
				if(nstar_group(i,1).gt.0) then
					n = 0
					l = nstar_group(i,2)
					do while(l .gt. 0)
						m = star_group(l,1)
						write(81,'( 2(1x,i10),1x,i1,10(1x,es12.5) )') m, idstars(m), insitustars(idstars(m)), &
						& xstars(m), ystars(m), zstars(m), vxstars(m), vystars(m), vzstars(m), mass_stars(m), initial_mass_stars(m), age_stars(m), SNII_stars(m)
						l = star_group(l,2)
						n = n + 1
					end do
					if( n .ne. nstar_group(i,1) ) then
						print *, 'The stellar linked list is not behaving as its supposed to'
						print *, nstar_group(i,1), n
						stop
					end if
				end if
				close(unit=81)
			end do
			deallocate( star_group )
			print *, 'allocating rest of clump properties'
			call fix_allocation_clump_prop()
			call allocate_clump_prop2()
			print *, '2d face on projection'
			do i=1,ngroup
				do j=ncell(i-1)+2,ncell(i)
					if( minval( abs(igroup(ncell(i-1)+1:j-1)-igroup(j)) + abs(jgroup(ncell(i-1)+1:j-1)-jgroup(j)) ) .gt. 0 ) then
						ncell_2d(i) = ncell_2d(i) + 1
					end if
				end do
				print *, 'ncell_3d',ncell(i)-ncell(i-1),'ncell_2d',ncell_2d(i)
			end do
	
			print *, 'masses, fractions, densities, centers, Rmax'
			do i=1,ngroup
				n = ncell(i) - ncell(i-1)										!!! # of cells in the clump
				clump_mass(i)     = clump_gas_mass(i)   + clump_star_mass(i)						!!! M_{sun}
				clump_dm_dens(i)  = clump_dm_mass(i)    / ( n           * (real(res,8)/1000.0_8)**3 )			!!! M_{sun} kpc^{-3}
				clump_gas_sig(i)  = clump_gas_mass(i)   / ( ncell_2d(i) * real(res,8)**2 )				!!! M_{sun} pc^{-2}
				clump_star_sig(i) = clump_star_mass(i)  / ( ncell_2d(i) * real(res,8)**2 )				!!! M_{sun} pc^{-2}
				clump_sig(i)      = clump_mass(i)       / ( ncell_2d(i) * real(res,8)**2 )				!!! M_{sun} pc^{-2}
				clump_Sig_SFR(i)  = clump_SFR(i)        / ( ncell_2d(i) * (real(res,8)/1000.0_8)**2 )			!!! M_{sun} yr^{-1} kpc^{-2}

				if( clump_gas_mass(i) .gt. 0.0_8 ) then
					clump_tau(i) = 1.d9 * clump_SFR(i) / clump_gas_mass(i)						!!! Gyr^{-1}
				else
					clump_tau(i) = 0.0_8
				end if

				if( clump_star_mass(i) .gt. 0.0_8 ) then
					clump_SSFR(i) = 1.d9 * clump_SFR(i) / clump_star_mass(i)					!!! Gyr^{-1}
				else
					clump_SSFR(i) = 0.0_8
				end if

				if( clump_mass(i) .gt. 0.0_8 ) then
					clump_gas_frac(i) = clump_gas_mass(i)   / clump_mass(i)
					!!! center of mass, disc frame, kpc
					rcar_clump_disc(i,1) = sum( dens_group(ncell(i-1)+1:ncell(i)) * real(igroup(ncell(i-1)+1:ncell(i)) - nR/2, 8) ) * (real(res,8)/1000.0_8) / clump_mass(i)
					rcar_clump_disc(i,2) = sum( dens_group(ncell(i-1)+1:ncell(i)) * real(jgroup(ncell(i-1)+1:ncell(i)) - nR/2, 8) ) * (real(res,8)/1000.0_8) / clump_mass(i)
					rcar_clump_disc(i,3) = sum( dens_group(ncell(i-1)+1:ncell(i)) * real(kgroup(ncell(i-1)+1:ncell(i)) - nH/2, 8) ) * (real(res,8)/1000.0_8) / clump_mass(i)
				else
					clump_gas_frac(i) = 0.0_8 
					rcar_clump_disc(i,1:3) = 0.0_8 
				end if

				if( clump_mass(i) + clump_dm_mass(i) .gt. 0.0_8 ) then
					clump_dm_frac(i) = clump_dm_mass(i)    / ( clump_mass(i) + clump_dm_mass(i) )
				else
					clump_dm_frac(i) = 0.0_8
				end if

				max_dens(i) = maxval( dens_group(ncell(i-1)+1:ncell(i)) ) / (real(res,8)**3)				!!! M_{sun} pc^{-3}
				temp_loc(:) = maxloc( dens_group(ncell(i-1)+1:ncell(i)) )
				group_center_dens(i,1) = igroup(ncell(i-1)+temp_loc(1))							!!! density peak
				group_center_dens(i,2) = jgroup(ncell(i-1)+temp_loc(1))							!!! density peak
				group_center_dens(i,3) = kgroup(ncell(i-1)+temp_loc(1))							!!! density peak
				rcar_clump_disc_density_peak(i,1) = real(group_center_dens(i,1) - nR/2, 8) * (real(res,8)/1000.0_8)
				rcar_clump_disc_density_peak(i,2) = real(group_center_dens(i,2) - nR/2, 8) * (real(res,8)/1000.0_8)
				rcar_clump_disc_density_peak(i,3) = real(group_center_dens(i,3) - nH/2, 8) * (real(res,8)/1000.0_8)
				group_rad(i,1) = maxval(sqrt( real(igroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,1),8)**2 + real(jgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,2),8)**2 + &
					& real(kgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,3),8)**2 ))*real(res,8)/1000.0_8	!!! kpc
				group_rad(i,3) = real((res/1000.0_4)*( (n/pi4_3)**(1.0_4/3.0_4) ),8)					!!! kpc
				!!! There are 4 different definitions for radius. 1 is maximal distance from clump density peak. 3 is sphere with same volume of clump
				!!! 3 IS THE ONLY ONE THAT WILL COUNT !!!
			end do
			!!! center of mass, box frame, kpc
			rcar_clump(:,1)              =  cos(theta)*cos(phi)*rcar_clump_disc(:,1)              - sin(phi)*rcar_clump_disc(:,2)              + sin(theta)*cos(phi)*rcar_clump_disc(:,3)	
			rcar_clump(:,2)              =  cos(theta)*sin(phi)*rcar_clump_disc(:,1)              + cos(phi)*rcar_clump_disc(:,2)              + sin(theta)*sin(phi)*rcar_clump_disc(:,3)	
			rcar_clump(:,3)              = -sin(theta)*         rcar_clump_disc(:,1)                                                           + cos(theta)*         rcar_clump_disc(:,3)

			rcar_clump_density_peak(:,1) =  cos(theta)*cos(phi)*rcar_clump_disc_density_peak(:,1) - sin(phi)*rcar_clump_disc_density_peak(:,2) + sin(theta)*cos(phi)*rcar_clump_disc_density_peak(:,3)	
			rcar_clump_density_peak(:,2) =  cos(theta)*sin(phi)*rcar_clump_disc_density_peak(:,1) + cos(phi)*rcar_clump_disc_density_peak(:,2) + sin(theta)*sin(phi)*rcar_clump_disc_density_peak(:,3)	
			rcar_clump_density_peak(:,3) = -sin(theta)*         rcar_clump_disc_density_peak(:,1)                                              + cos(theta)*         rcar_clump_disc_density_peak(:,3)
	
			vcyl_clump(:,1) =  cos(theta)*cos(phi)*vcar_clump(:,1) + cos(theta)*sin(phi)*vcar_clump(:,2) - sin(theta)*vcar_clump(:,3)	! Vx',  km/s
			vcyl_clump(:,2) = -sin(phi)*vcar_clump(:,1) + cos(phi)*vcar_clump(:,2)								! Vy',  km/s
			vcyl_clump(:,3) = ( vcyl_clump(:,1)*rcar_clump_disc(:,1) + vcyl_clump(:,2)*rcar_clump_disc(:,2)) / sqrt(rcar_clump_disc(:,1)**2 + rcar_clump_disc(:,2)**2)																! Vr,   km/s
			vcyl_clump(:,2) = (-vcyl_clump(:,1)*rcar_clump_disc(:,2) + vcyl_clump(:,2)*rcar_clump_disc(:,1)) / sqrt(rcar_clump_disc(:,1)**2 + rcar_clump_disc(:,2)**2)																! Vphi, km/s
			vcyl_clump(:,1) = vcyl_clump(:,3)												! Vr,   km/s
			vcyl_clump(:,3) =  sin(theta)*cos(phi)*vcar_clump(:,1) + sin(theta)*sin(phi)*vcar_clump(:,2) + cos(theta)*vcar_clump(:,3)	! Vz',  km/s
			print *, 'identify bulge'
			if( bulge_ind(1) .gt. 0 ) then
				clump_comment(bulge_ind(1)) = ' bulge'
			end if
			print *, 'normalized distance and height'
			do i=1,ngroup
!				if( maxval(abs(rcar_clump(i,:))) .lt. 2.0_8 * real(res,8) / 1000.0_8 ) then
!					clump_comment(i) = ' bulge'
!				end if
				clump_dist(i) = sqrt(rcar_clump_disc(i,1)**2 + rcar_clump_disc(i,2)**2) / real(Rdisc(snapshot), 8)
				clump_height(i) = rcar_clump_disc(i,3) / real(Hdisc(snapshot), 8)
			end do
			print *, 'R90, inertia tensor and shape parameter'
			n = ceiling( 1000.0_8*maxval(group_rad(1:ngroup,1)) / res )	!!! maximal Rmax in units of grid resolution
			allocate( rad_bins(n), stat=i )					!!! deallocated in this subroutine !!!
			if(i.ne.0) then
				print *, 'error in allocation mass_bins. stat= ', i
				stop
			end if
			allocate( mass_bins(n), stat=i )				!!! deallocated in this subroutine !!!
			if(i.ne.0) then
				print *, 'error in allocation mass_bins. stat= ', i
				stop
			end if
			rad_bins(:) = (/ (i,i=1,n) /) * real(res,8)			!!! pc

			write(filename,'(a,a,a,a,a)') './',trim(dirname),'/extra_clump_parameters/',trim(gal),'.out'
			open(unit=40,file=filename)
			temp_var2 = real(res,8)/1000.0_8
			do i=1,ngroup
				mass_bins(:) = 0.0_8
				do j=ncell(i-1)+1,ncell(i)
					temp_var = sqrt( real(igroup(j)-group_center_dens(i,1), 8)**2 + real(jgroup(j)-group_center_dens(i,2), 8)**2 + real(kgroup(j)-group_center_dens(i,3), 8)**2 ) * real(res,8) !pc
					k=1
					do while(temp_var .ge. rad_bins(k) .and. k .lt. n)
						k = k+1
					end do
					mass_bins(k:n) = mass_bins(k:n) + dens_group(j)
				end do
				k=2
				do while(mass_bins(k) .le. 0.90_8*mass_bins(n) .and. k .lt. n)
					k = k+1
				end do
				group_rad(i,2) = ( rad_bins(k-1) + ( ( rad_bins(k)-rad_bins(k-1) )/( mass_bins(k)-mass_bins(k-1) ) ) * ( 0.90_8*mass_bins(n) - mass_bins(k-1) ) ) / 1000.0_8	!!! kpc
				if (group_rad(i,2) .gt. rad_bins(k)/1000.0_8) then
					group_rad(i,2) = (0.5_4*(rad_bins(k) + rad_bins(k-1))) / 1000.0_8
				end if	!!! There are 4 different definitions for radius. 2 is the radius of a sphere centered on density peak containing 90% of clump baryonic mass

				!!! Calculate Inertia Tensor !!!
				inert(:,:) = 0.0_8
				v(:,:) = 0.0_8
				inert_eigen(:,:) = 0.0_8
				inert(1,1) =  sum( dens_group(ncell(i-1)+1:ncell(i)) * (temp_var2**2)*real( (jgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,2))**2 + (kgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,3))**2, 8), &
					& mask = sqrt( real((igroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,1))**2, 8) + real((jgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,2))**2, 8) + &
					& real((kgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,3))**2, 8) ) * temp_var2 .le. group_rad(i,2) )

				inert(2,2) =  sum( dens_group(ncell(i-1)+1:ncell(i)) * (temp_var2**2)*real( (igroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,1))**2 + (kgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,3))**2, 8), &
					& mask = sqrt( real((igroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,1))**2, 8) + real((jgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,2))**2, 8) + &
					& real((kgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,3))**2, 8) ) * temp_var2 .le. group_rad(i,2) )

				inert(3,3) =  sum( dens_group(ncell(i-1)+1:ncell(i)) * (temp_var2**2)*real( (igroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,1))**2 + (jgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,2))**2, 8), &
					& mask = sqrt( real((igroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,1))**2, 8) + real((jgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,2))**2, 8) + &
					& real((kgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,3))**2, 8) ) * temp_var2 .le. group_rad(i,2) )

				inert(1,2) = -sum( dens_group(ncell(i-1)+1:ncell(i)) * (temp_var2**2)*real( (igroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,1))    * (jgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,2)), 8), &
					& mask = sqrt( real((igroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,1))**2, 8) + real((jgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,2))**2, 8) + &
					& real((kgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,3))**2, 8) ) * temp_var2 .le. group_rad(i,2) )

				inert(1,3) = -sum( dens_group(ncell(i-1)+1:ncell(i)) * (temp_var2**2)*real( (igroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,1))    * (kgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,3)), 8), &
					& mask = sqrt( real((igroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,1))**2, 8) + real((jgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,2))**2, 8) + &
					& real((kgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,3))**2, 8) ) * temp_var2 .le. group_rad(i,2) )

				inert(2,3) = -sum( dens_group(ncell(i-1)+1:ncell(i)) * (temp_var2**2)*real( (jgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,2))    * (kgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,3)), 8), &
					& mask = sqrt( real((igroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,1))**2, 8) + real((jgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,2))**2, 8) + &
					& real((kgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,3))**2, 8) ) * temp_var2 .le. group_rad(i,2) )
				inert(2,1) = inert(1,2)
				inert(3,1) = inert(1,3)
				inert(3,2) = inert(2,3)
				inert2(:,:) = inert(:,:)
				call jacobi(inert2, 3, 3, inert_eigen(i,:), v, nrot)
	
				!!! Order Inertia Tensor Eigenvalues !!!
				temp_var = maxval(inert_eigen(i,:))
				temp_loc = maxloc(inert_eigen(i,:))
				inert_eigen(i,temp_loc(1)) = inert_eigen(i,1)
				inert_eigen(i,1) = temp_var
				temp_var = maxval(inert_eigen(i,2:3))
				temp_loc = maxloc(inert_eigen(i,2:3))
				inert_eigen(i,temp_loc(1)+1) = inert_eigen(i,2)
				inert_eigen(i,2) = temp_var
	
				!!! Axes assuming homogeneous ellipsoid !!!
				cI = sqrt( (2.5_8/clump_mass(i)) * (inert_eigen(i,2) - inert_eigen(i,1) + inert_eigen(i,3)) )
				bI = sqrt( (2.5_8/clump_mass(i)) * (inert_eigen(i,3) - inert_eigen(i,2) + inert_eigen(i,1)) )
				aI = sqrt( (2.5_8/clump_mass(i)) * (inert_eigen(i,1) - inert_eigen(i,3) + inert_eigen(i,2)) )
		
				inert_eigen(i,:) = inert_eigen(i,:) / sqrt(dot_product(inert_eigen(i,:), inert_eigen(i,:)))
				eta_group(i,1) = inert_eigen(i,3)/inert_eigen(i,1)
				group_rad(i,4) = (aI * bI * cI)**(1.0_8/3.0_8)	!!! There are 4 different definitions for radius. 4 is radius of sphere with same volume as ellipsoid having these axes	
				write(40,'(1x,i3,2(1x,i5),15(1x,es12.4))') i, ncell(i)-ncell(i-1), ncell_2d(i), max_dens(i), max_res(i), mean_residual(i), inert_eigen(i,1:3), &
				& eta_group(i,1), aI, bI, cI, cI/aI, group_rad(i,1:4)
				do j=1,n
					write(40,*) rad_bins(j)/1000, mass_bins(j), mass_bins(j)/mass_bins(n)
				end do
				write(40,*) ''
			end do
			deallocate(rad_bins, mass_bins, max_dens, max_res)
			close(unit=40)

			print *, 'Alternate shape parameters'
			do i=1,ngroup
				!!! Calculate Inertia Tensor !!!
				inert(:,:) = 0.0_8
				v(:,:) = 0.0_8
				inert_eigen(:,:) = 0.0_8
				inert(1,1) =  sum( dens_group(ncell(i-1)+1:ncell(i)) * (temp_var2**2)*real( (jgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,2))**2 + (kgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,3))**2, 8) )
				inert(2,2) =  sum( dens_group(ncell(i-1)+1:ncell(i)) * (temp_var2**2)*real( (igroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,1))**2 + (kgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,3))**2, 8) )
				inert(3,3) =  sum( dens_group(ncell(i-1)+1:ncell(i)) * (temp_var2**2)*real( (igroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,1))**2 + (jgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,2))**2, 8) )
				inert(1,2) = -sum( dens_group(ncell(i-1)+1:ncell(i)) * (temp_var2**2)*real( (igroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,1))*(jgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,2)), 8) )
				inert(1,3) = -sum( dens_group(ncell(i-1)+1:ncell(i)) * (temp_var2**2)*real( (igroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,1))*(kgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,3)), 8) )
				inert(2,3) = -sum( dens_group(ncell(i-1)+1:ncell(i)) * (temp_var2**2)*real( (jgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,2))*(kgroup(ncell(i-1)+1:ncell(i))-group_center_dens(i,3)), 8) )
				inert(2,1) = inert(1,2)
				inert(3,1) = inert(1,3)
				inert(3,2) = inert(2,3)
				inert2(:,:) = inert(:,:)
				call jacobi(inert2, 3, 3, inert_eigen(i,:), v, nrot)
				eta_group(i,2) = minval(inert_eigen(i,:)) / maxval(inert_eigen(i,:))
			end do
			do i=1,ngroup
				!!! Calculate Inertia Tensor !!!
				inert(:,:) = 0.0_8
				v(:,:) = 0.0_8
				inert_eigen(:,:) = 0.0_8
				inert(1,1) =  sum( dens_group(ncell(i-1)+1:ncell(i)) * ( (temp_var2*(jgroup(ncell(i-1)+1:ncell(i))-nR/2)-rcar_clump_disc(i,2))**2 + (temp_var2*(kgroup(ncell(i-1)+1:ncell(i))-nH/2)-rcar_clump_disc(i,3))**2 ) )
				inert(2,2) =  sum( dens_group(ncell(i-1)+1:ncell(i)) * ( (temp_var2*(igroup(ncell(i-1)+1:ncell(i))-nR/2)-rcar_clump_disc(i,1))**2 + (temp_var2*(kgroup(ncell(i-1)+1:ncell(i))-nH/2)-rcar_clump_disc(i,3))**2 ) )
				inert(3,3) =  sum( dens_group(ncell(i-1)+1:ncell(i)) * ( (temp_var2*(igroup(ncell(i-1)+1:ncell(i))-nR/2)-rcar_clump_disc(i,1))**2 + (temp_var2*(jgroup(ncell(i-1)+1:ncell(i))-nR/2)-rcar_clump_disc(i,2))**2 ) )
				inert(1,2) = -sum( dens_group(ncell(i-1)+1:ncell(i)) * ( (temp_var2*(igroup(ncell(i-1)+1:ncell(i))-nR/2)-rcar_clump_disc(i,1)) * (temp_var2*(jgroup(ncell(i-1)+1:ncell(i))-nR/2)-rcar_clump_disc(i,2)) ) )
				inert(1,3) = -sum( dens_group(ncell(i-1)+1:ncell(i)) * ( (temp_var2*(igroup(ncell(i-1)+1:ncell(i))-nR/2)-rcar_clump_disc(i,1)) * (temp_var2*(kgroup(ncell(i-1)+1:ncell(i))-nH/2)-rcar_clump_disc(i,3)) ) )
				inert(2,3) = -sum( dens_group(ncell(i-1)+1:ncell(i)) * ( (temp_var2*(jgroup(ncell(i-1)+1:ncell(i))-nR/2)-rcar_clump_disc(i,2)) * (temp_var2*(kgroup(ncell(i-1)+1:ncell(i))-nH/2)-rcar_clump_disc(i,3)) ) )
				inert(2,1) = inert(1,2)
				inert(3,1) = inert(1,3)
				inert(3,2) = inert(2,3)
				inert2(:,:) = inert(:,:)
				call jacobi(inert2, 3, 3, inert_eigen(i,:), v, nrot)
				eta_group(i,3) = minval(inert_eigen(i,:)) / maxval(inert_eigen(i,:))
			end do
			deallocate( igroup, jgroup, kgroup, dens_group, group_center_dens, inert_eigen  )
	
			print *, 'stellar spherical mass test'
			ntrack = (Nstars_list - mod(Nstars_list,20)) / 20						!!! This just helps keep track of where I am in the loop
			do i=1,Nstars_list
				if(i .eq. 1) print *, 'i of Nstars_list',i,'of',Nstars_list
				if ( mod(i,ntrack) .eq. 0 ) then
					print*, 'i of Nstars_list',i,'of',Nstars_list
				end if
				m = star_list(i)
				do j=1,ngroup
					if( sqrt( (xstars(m)-rcar_clump(j,1))**2 + (ystars(m)-rcar_clump(j,2))**2 + (zstars(m)-rcar_clump(j,3))**2 ) .le. group_rad(j,3) ) then
						spher_mass_test(j,1)  = spher_mass_test(j,1)  + mass_stars(m)
					end if
					if( sqrt( (xstars(m)-rcar_clump(j,1))**2 + (ystars(m)-rcar_clump(j,2))**2 + (zstars(m)-rcar_clump(j,3))**2 ) .le. group_rad(j,4) ) then
						spher_mass_test(j,2) = spher_mass_test(j,2) + mass_stars(m)
					end if
					if( sqrt( (xstars(m)-rcar_clump_density_peak(j,1))**2 + (ystars(m)-rcar_clump_density_peak(j,2))**2 + (zstars(m)-rcar_clump_density_peak(j,3))**2 ) .le. group_rad(j,3) ) then
						v_clump_frame(j,1)    = v_clump_frame(j,1)    + mass_stars(m)*vxstars(m)
						v_clump_frame(j,2)    = v_clump_frame(j,2)    + mass_stars(m)*vystars(m)
						v_clump_frame(j,3)    = v_clump_frame(j,3)    + mass_stars(m)*vzstars(m)
						v_sq_clump_frame(j,1) = v_sq_clump_frame(j,1) + mass_stars(m)*vxstars(m)**2
						v_sq_clump_frame(j,2) = v_sq_clump_frame(j,2) + mass_stars(m)*vystars(m)**2
						v_sq_clump_frame(j,3) = v_sq_clump_frame(j,3) + mass_stars(m)*vzstars(m)**2
						mass_clump_frame(j)   = mass_clump_frame(j)   + mass_stars(m)
					end if
				end do
			end do
			deallocate( xstars, ystars, zstars, vxstars, vystars, vzstars, mass_stars, initial_mass_stars, idstars, age_stars, SNII_stars )
			deallocate( star_list )
		
			print *, 'gas background and spherical mass test'
			vgas_back_global = 0.0_4
			rgas_back_global = 0.0_4
			mgas_back_global = 0.0_4
			ntrack = (Ngas_list - mod(Ngas_list,20)) / 20						!!! This just helps keep track of where I am in the loop
			do i=1,Ngas_list
				if(i .eq. 1) print *, 'i of Ngas_list',i,'of',Ngas_list
				if ( mod(i,ntrack) .eq. 0 ) then
					print*, 'i of Ngas_list',i,'of',Ngas_list
				end if
				m = gas_list(i)
				mass_temporary = 0.03363_8*real(density_gas(m)*((cell_size_gas(m))**3),8)

				call split0( xgas(m), ygas(m), zgas(m), (/ 0.0_4, 0.0_4, 0.0_4 /), vxgas(m), vygas(m), vzgas(m), (/ 0.0_4, 0.0_4, 0.0_4 /) )
				do j=1,ngroup
					if( sqrt( (real(xgas(m),8)-rcar_clump(j,1))**2 + (real(ygas(m),8)-rcar_clump(j,2))**2 + (real(zgas(m),8)-rcar_clump(j,3))**2 ) .le. group_rad(j,3) ) then
						spher_mass_test(j,1) = spher_mass_test(j,1) + mass_temporary
					end if
					if( sqrt( (real(xgas(m),8)-rcar_clump(j,1))**2 + (real(ygas(m),8)-rcar_clump(j,2))**2 + (real(zgas(m),8)-rcar_clump(j,3))**2 ) .le. group_rad(j,4) ) then
						spher_mass_test(j,2) = spher_mass_test(j,2) + mass_temporary
					end if
					if( sqrt( (real(xgas(m),8)-rcar_clump_density_peak(j,1))**2 + (real(ygas(m),8)-rcar_clump_density_peak(j,2))**2 + (real(zgas(m),8)-rcar_clump_density_peak(j,3))**2 ) .le. group_rad(j,3) ) then
						v_clump_frame(j,1)    = v_clump_frame(j,1)    + mass_temporary*real(vxgas(m),8)
						v_clump_frame(j,2)    = v_clump_frame(j,2)    + mass_temporary*real(vygas(m),8)
						v_clump_frame(j,3)    = v_clump_frame(j,3)    + mass_temporary*real(vzgas(m),8)
						v_sq_clump_frame(j,1) = v_sq_clump_frame(j,1) + mass_temporary*real(vxgas(m)**2,8)
						v_sq_clump_frame(j,2) = v_sq_clump_frame(j,2) + mass_temporary*real(vygas(m)**2,8)
						v_sq_clump_frame(j,3) = v_sq_clump_frame(j,3) + mass_temporary*real(vzgas(m)**2,8)
						mass_clump_frame(j)   = mass_clump_frame(j)   + mass_temporary
					end if
					!!! Ring of total width 1kpc centered on the clump, and total thickness 1 kpc centered on the midplane
					if( rprime(1) .ge. sngl(sqrt( rcar_clump_disc(j,1)**2 + rcar_clump_disc(j,2)**2 )) - 0.5_4 .and. &
					&   rprime(1) .le. sngl(sqrt( rcar_clump_disc(j,1)**2 + rcar_clump_disc(j,2)**2 )) + 0.5_4 .and. &
					&   abs(zprime(1)) .le. 0.5_4 .and. temperature_gas(m) .lt. max_T ) then
						mgas_back(j)      = mgas_back(j)      + sngl(mass_temporary)
						vgas_back(j,1)    = vgas_back(j,1)    + sngl(mass_temporary) *    (vxprime(1)*xprime(1) + vyprime(1)*yprime(1))/rprime(1)
						vgas_back(j,2)    = vgas_back(j,2)    + sngl(mass_temporary) *    (vyprime(1)*xprime(1) - vxprime(1)*yprime(1))/rprime(1)
						vgas_back(j,3)    = vgas_back(j,3)    + sngl(mass_temporary) *     vzprime(1)
						vgas_sq_back(j,1) = vgas_sq_back(j,1) + sngl(mass_temporary) * ( ((vxprime(1)*xprime(1) + vyprime(1)*yprime(1))/rprime(1))**2 )
						vgas_sq_back(j,2) = vgas_sq_back(j,2) + sngl(mass_temporary) * ( ((vyprime(1)*xprime(1) - vxprime(1)*yprime(1))/rprime(1))**2 )
						vgas_sq_back(j,3) = vgas_sq_back(j,3) + sngl(mass_temporary) * (   vzprime(1)**2 )
					end if
				end do
				if( rprime(1) .ge. 0.5_4 * Rdisc(snapshot) .and. rprime(1) .le. Rdisc(snapshot) .and. abs(zprime(1)) .le. 0.5_4 .and. temperature_gas(m) .lt. max_T ) then
					mgas_back_global = mgas_back_global + sngl(mass_temporary)
					vgas_back_global = vgas_back_global + sngl(mass_temporary) * (vyprime(1)*xprime(1) - vxprime(1)*yprime(1))/rprime(1)
					rgas_back_global = rgas_back_global + sngl(mass_temporary) * rprime(1)
				end if
				call deallocate_primes()
			end do
			vgas_back_global = vgas_back_global / mgas_back_global
			rgas_back_global = rgas_back_global / mgas_back_global
			td_global = ( rgas_back_global / (1.023_4 * vgas_back_global) ) * 1000.0_8				!!! kpc / (kpc/Gyr) * 1000 --> Myr
			do i=1,ngroup
				vgas_back(i,:) = vgas_back(i,:) / mgas_back(i)
				vgas_sq_back(i,:) = vgas_sq_back(i,:) / mgas_back(i)
				sigma_gas_back(i,:) = sqrt( vgas_sq_back(i,:) - ((vgas_back(i,:))**2) )
	
				td(i)     = ( sqrt(rcar_clump_disc(i,1)**2 + rcar_clump_disc(i,2)**2) / real(1.023_4 * vgas_back(i,2),8) ) * 1000.0_8	!!! kpc / (kpc/Gyr) * 1000 --> Myr
				tff(i) = sqrt( real(pi4_3,8) * group_rad(i,3)**3 / ( (1.023_8**2)*G * clump_mass(i) ) ) * 1000.0_8	!!! sqrt( kpc^3 / ( kpc * ((kpc/Gyr)^2) * M_{sun}^{-1} * M_{sun} ) ) * 1000 --	> Myr
			end do
			do i=1,ngroup
				v_clump_frame(i,:)     = v_clump_frame(i,:)    / mass_clump_frame(i)
				v_sq_clump_frame(i,:)  = v_sq_clump_frame(i,:) / mass_clump_frame(i)
				sigma_clump_frame(i,1) = v_sq_clump_frame(i,1) - v_clump_frame(i,1)**2
				sigma_clump_frame(i,2) = v_sq_clump_frame(i,2) - v_clump_frame(i,2)**2
				sigma_clump_frame(i,3) = v_sq_clump_frame(i,3) - v_clump_frame(i,3)**2
			end do

			print *, 'dark matter background'
			ntrack = (Ndm(snapshot) - mod(Ndm(snapshot),20)) / 20								!!! This just helps keep track of where I am in the loop
			do i=1,Ndm(snapshot)
				if( i .eq. 1 ) print*, 'i of Ndm',i,'of',Ndm(snapshot)
				if( mod(i,ntrack) .eq. 0 ) then
					print*, 'i of Ndm',i,'of',Ndm(snapshot)
				end if
				do j=1,ngroup
					if( sqrt( (xdm(i)-rcar_clump_density_peak(j,1))**2 + (ydm(i)-rcar_clump_density_peak(j,2))**2 + (zdm(i)-rcar_clump_density_peak(j,3))**2 ) .le. group_rad(j,3) ) then
						mass_clump_frame(j) = mass_clump_frame(j) + mass_dm(i)
					end if
					if( sqrt( xdm(i)**2 + ydm(i)**2 + zdm(i)**2 ) .ge. sqrt( rcar_clump(j,1)**2 + rcar_clump(j,2)**2  + rcar_clump(j,3)**2 ) - 0.5_8 .and. &
					&   sqrt( xdm(i)**2 + ydm(i)**2 + zdm(i)**2 ) .le. sqrt( rcar_clump(j,1)**2 + rcar_clump(j,2)**2  + rcar_clump(j,3)**2 ) + 0.5_8) then
						dm_back(j) = dm_back(j) + mass_dm(i)
					end if
				end do
			end do
			deallocate( xdm, ydm, zdm, mass_dm )
			do i=1,ngroup
				dm_back(i) = dm_back(i) / ( real(pi4_3,8) * ( (sqrt(rcar_clump(i,1)**2 + rcar_clump(i,2)**2  + rcar_clump(i,3)**2) + 0.5_8)**3 - &
					& (sqrt(rcar_clump(i,1)**2 + rcar_clump(i,2)**2  + rcar_clump(i,3)**2) - 0.5_8)**3 ) )		!!! M_{sun} kpc^{-3}

				alpha_vir(i) = ( 5.0_8 / 3.0_8 ) * ( sigma_clump_frame(i,1) + sigma_clump_frame(i,2) + sigma_clump_frame(i,3) ) * group_rad(i,3) / ( G * mass_clump_frame(i) )
			end do

			print *, 'Gas inflow and outflow'
			ntrack = (Ngas_list - mod(Ngas_list,20)) / 20						!!! This just helps keep track of where I am in the loop
			do i=1,Ngas_list
				if(i .eq. 1) print *, 'i of Ngas_list',i,'of',Ngas_list
				if ( mod(i,ntrack) .eq. 0 ) then
					print*, 'i of Ngas_list',i,'of',Ngas_list
				end if
				m = gas_list(i)
				do j=1,ngroup
					Vesc1 = sqrt( 2.0_8 * G * mass_clump_frame(j) / group_rad(j,3) )
					Vesc2 = sqrt( 2.0_8 * G * mass_clump_frame(j) / (1.5_4 * group_rad(j,3)) )
					Vesc3 = sqrt( 2.0_8 * G * mass_clump_frame(j) / (2.0_4 * group_rad(j,3)) )
					rthresh = sqrt( (xgas(m)-sngl(rcar_clump_density_peak(j,1)))**2 + (ygas(m)-sngl(rcar_clump_density_peak(j,2)))**2 + (zgas(m)-sngl(rcar_clump_density_peak(j,3)))**2 )
					if ( rthresh .le. max( 3.0_4*sngl(group_rad(j,3)), 1.5_4 )) then
						if(cell_size_gas(m) < 0.75_4*res) then
							call split0( xgas(m), ygas(m), zgas(m), sngl(rcar_clump_density_peak(j,:)), vxgas(m), vygas(m), vzgas(m), sngl(v_clump_frame(j,:)) )
						else if(cell_size_gas(m) < 1.5_4*res) then
							call split1( xgas(m), ygas(m), zgas(m), sngl(rcar_clump_density_peak(j,:)), cell_size_gas(m), vxgas(m), vygas(m), vzgas(m), sngl(v_clump_frame(j,:)) )
						else if(cell_size_gas(m) < 3.0_4*res) then
							call split2( xgas(m), ygas(m), zgas(m), sngl(rcar_clump_density_peak(j,:)), cell_size_gas(m), vxgas(m), vygas(m), vzgas(m), sngl(v_clump_frame(j,:)) )
						else if(cell_size_gas(m) < 5.0_4*res) then
							call split3( xgas(m), ygas(m), zgas(m), sngl(rcar_clump_density_peak(j,:)), cell_size_gas(m), vxgas(m), vygas(m), vzgas(m), sngl(v_clump_frame(j,:)) )
						else
							call split4( xgas(m), ygas(m), zgas(m), sngl(rcar_clump_density_peak(j,:)), cell_size_gas(m), vxgas(m), vygas(m), vzgas(m), sngl(v_clump_frame(j,:)) )
						end if
						l1 = size(xprime(:))
						split = log(1.0_4*l1)/log(8.0_4)
						mass_temporary = 0.03363_8*real(density_gas(m)*((cell_size_gas(m)/(2.0_4**split))**3),8)
						do l=1,l1
							if( sqrt(rprime(l)**2 + zprime(l)**2) .ge. sngl(group_rad(j,3)) - (res/1000.0_4) .and. &
							&   sqrt(rprime(l)**2 + zprime(l)**2) .le. sngl(group_rad(j,3)) + (res/1000.0_4) ) then
								vtemp = (vxprime(l)*xprime(l) + vyprime(l)*yprime(l) + vzprime(l)*zprime(l)) / sqrt(rprime(l)**2+zprime(l)**2)
								if(vtemp .ge. Vesc1 ) then
									Mgas_out(j,8) = Mgas_out(j,8) + mass_temporary*1.023_8*real(abs(vtemp),8)	!!! M_sun * kpc/Gyr
									Mgas_out(j,4) = Mgas_out(j,4) + mass_temporary*1.023_8*real(abs(vtemp),8)	!!! M_sun * kpc/Gyr
									Mgas_out(j,1) = Mgas_out(j,1) + mass_temporary*1.023_8*real(abs(vtemp),8)	!!! M_sun * kpc/Gyr
								elseif(vtemp .ge. 0.0_4) then
									Mgas_out(j,1) = Mgas_out(j,1) + mass_temporary*1.023_8*real(abs(vtemp),8)	!!! M_sun * kpc/Gyr
									if( sqrt( vxprime(l)**2 + vyprime(l)**2 + vzprime(l)**2 ) .ge. Vesc1 ) then
										Mgas_out(j,8) = Mgas_out(j,8) + mass_temporary*1.023_8*real(abs(vtemp),8)	!!! M_sun * kpc/Gyr
									end if
								elseif(vtemp .lt. 0.0_4) then
									Mgas_in(j,1)  = Mgas_in(j,1)  + mass_temporary*1.023_8*real(abs(vtemp),8)	!!! M_sun * kpc/Gyr
								end if
							end if
							if( sqrt(rprime(l)**2 + zprime(l)**2) .ge. 1.5_4 * sngl(group_rad(j,3)) - (res/1000.0_4) .and. &
							&   sqrt(rprime(l)**2 + zprime(l)**2) .le. 1.5_4 * sngl(group_rad(j,3)) + (res/1000.0_4)  ) then
								vtemp = (vxprime(l)*xprime(l) + vyprime(l)*yprime(l) + vzprime(l)*zprime(l)) / sqrt(rprime(l)**2+zprime(l)**2)
								if(vtemp .ge. Vesc2 ) then
									Mgas_out(j,9) = Mgas_out(j,9) + mass_temporary*1.023_8*real(abs(vtemp),8)	!!! M_sun * kpc/Gyr
									Mgas_out(j,5) = Mgas_out(j,5) + mass_temporary*1.023_8*real(abs(vtemp),8)	!!! M_sun * kpc/Gyr
									Mgas_out(j,2) = Mgas_out(j,2) + mass_temporary*1.023_8*real(abs(vtemp),8)	!!! M_sun * kpc/Gyr
								elseif(vtemp .ge. 0.0_4) then
									Mgas_out(j,2) = Mgas_out(j,2) + mass_temporary*1.023_8*real(abs(vtemp),8)	!!! M_sun * kpc/Gyr
									if( sqrt( vxprime(l)**2 + vyprime(l)**2 + vzprime(l)**2 ) .ge. Vesc2 ) then
										Mgas_out(j,9) = Mgas_out(j,9) + mass_temporary*1.023_8*real(abs(vtemp),8)	!!! M_sun * kpc/Gyr
									end if
								elseif(vtemp .lt. 0.0_4) then
									Mgas_in(j,2)  = Mgas_in(j,2)  + mass_temporary*1.023_8*real(abs(vtemp),8)	!!! M_sun * kpc/Gyr
								end if
							end if
							if( sqrt(rprime(l)**2 + zprime(l)**2) .ge. 2.0_4 * sngl(group_rad(j,3)) - (res/1000.0_4) .and. &
							&   sqrt(rprime(l)**2 + zprime(l)**2) .le. 2.0_4 * sngl(group_rad(j,3)) + (res/1000.0_4)  ) then
								vtemp = (vxprime(l)*xprime(l) + vyprime(l)*yprime(l) + vzprime(l)*zprime(l)) / sqrt(rprime(l)**2+zprime(l)**2)
								if(vtemp .ge. Vesc3) then
									Mgas_out(j,10) = Mgas_out(j,10) + mass_temporary*1.023_8*real(abs(vtemp),8)	!!! M_sun * kpc/Gyr
									Mgas_out(j,6)  = Mgas_out(j,6)  + mass_temporary*1.023_8*real(abs(vtemp),8)	!!! M_sun * kpc/Gyr
									Mgas_out(j,3)  = Mgas_out(j,3)  + mass_temporary*1.023_8*real(abs(vtemp),8)	!!! M_sun * kpc/Gyr
								elseif(vtemp .ge. 0.0_4) then
									Mgas_out(j,3) = Mgas_out(j,3) + mass_temporary*1.023_8*real(abs(vtemp),8)	!!! M_sun * kpc/Gyr
									if( sqrt( vxprime(l)**2 + vyprime(l)**2 + vzprime(l)**2 ) .ge. Vesc3 ) then
										Mgas_out(j,10) = Mgas_out(j,10) + mass_temporary*1.023_8*real(abs(vtemp),8)	!!! M_sun * kpc/Gyr
									end if
								elseif(vtemp .lt. 0.0_4 ) then
									Mgas_in(j,3) = Mgas_in(j,3)   + mass_temporary*1.023_8*real(abs(vtemp),8)
								end if
							end if
							if( abs(xprime(l)) .le. 0.5_4 .and. abs(yprime(l)) .le. 0.5_4 ) then
								vtemp = vzprime(l)
								if( zprime(l) .ge. 1.0_4 - (res/1000.0_4) .and. zprime(l) .le. 1.0_4 + (res/1000.0_4)  ) then
									if( vtemp .gt. 0.0_4 ) then
										Mgas_out(j,7) = Mgas_out(j,7) + mass_temporary*1.023_8*real(abs(vtemp),8) ! M_sun * kpc/Gyr
									end if
								elseif( zprime(l) .le. -1.0_4 + (res/1000.0_4) .and. zprime(l) .ge. -1.0_4 - (res/1000.0_4)  ) then
									if( vtemp .lt. 0.0_4 ) then
										Mgas_out(j,7) = Mgas_out(j,7) + mass_temporary*1.023_8*real(abs(vtemp),8) ! M_sun * kpc/Gyr
									end if
								end if
							end if
						end do
					end if					
				end do
			end do
			do i=1,ngroup
				Mgas_in(i,:)  = Mgas_in(i,:)  / ( 1.d9 * 2.0_4 * (res/1000.0_4) )		!!! M_sun / yr
				Mgas_out(i,:) = Mgas_out(i,:) / ( 1.d9 * 2.0_4 * (res/1000.0_4) )		!!! M_sun / yr
			end do
			deallocate( xgas, ygas, zgas, vxgas, vygas, vzgas, density_gas, cell_size_gas, temperature_gas, SNII_gas )
			deallocate( gas_list )
			saxis3(:) = Ldisc(:,snapshot)
			call axes(saxis1,saxis2,saxis3)
	
			print *, 'merger trees and ex-situ clumps'
			allocate( closest_star1(ngroup,1), stat=i )		!!! deallocated at end of subroutine !!!
			if(i.ne.0) then
				print *, 'error in allocation closest_star. stat= ', i
				stop
			end if
			do i=1,ngroup
				nstar_clump = 0
				if( nstar_group(i,1) > 0 ) then
					open(unit=83,file=clump_stars_filename(i))			!!! stars in clump
					read(83,'(i10)') nstar_clump
					allocate( ID_star(nstar_clump), stat=j )			!!! deallocated after merger trees !!!
					if(j.ne.0) then
						print *, 'error in allocation ID_star. stat= ', j
						stop
					end if
					allocate( mass_star(nstar_clump), stat=j )			!!! deallocated after merger trees !!!
					if(j.ne.0) then
						print *, 'error in allocation mass_star. stat= ', j
						stop
					end if
					allocate( age_star(nstar_clump), stat=j )			!!! deallocated after merger trees !!!
					if(j.ne.0) then
						print *, 'error in allocation age_star. stat= ', j
						stop
					end if
					allocate( IND_star(nstar_clump), stat=j )			!!! deallocated after reading !!!
					if(j.ne.0) then
						print *, 'error in allocation IND_star. stat= ', j
						stop
					end if
					allocate( IS_star(nstar_clump), stat=j )			!!! deallocated after reading !!!
					if(j.ne.0) then
						print *, 'error in allocation IS_star. stat= ', j
						stop
					end if
					allocate( xstar(nstar_clump), stat=j )				!!! deallocated after reading !!!
					if(j.ne.0) then
						print *, 'error in allocation xstar. stat= ', j
						stop
					end if
					allocate( ystar(nstar_clump), stat=j )				!!! deallocated after reading !!!
					if(j.ne.0) then
						print *, 'error in allocation ystar. stat= ', j
						stop
					end if
					allocate( zstar(nstar_clump), stat=j )				!!! deallocated after reading !!!
					if(j.ne.0) then
						print *, 'error in allocation zstar. stat= ', j
						stop
					end if
					allocate( vxstar(nstar_clump), stat=j )				!!! deallocated after reading !!!
					if(j.ne.0) then
						print *, 'error in allocation vxstar. stat= ', j
						stop
					end if
					allocate( vystar(nstar_clump), stat=j )				!!! deallocated after reading !!!
					if(j.ne.0) then
						print *, 'error in allocation vystar. stat= ', j
						stop
					end if
					allocate( vzstar(nstar_clump), stat=j )				!!! deallocated after reading !!!
					if(j.ne.0) then
						print *, 'error in allocation vzstar. stat= ', j
						stop
					end if
					allocate( met_star(nstar_clump), stat=j )			!!! deallocated after reading !!!
					if(j.ne.0) then
						print *, 'error in allocation met_star. stat= ', j
						stop
					end if
					allocate( initial_mass_star(nstar_clump), stat=j )		!!! deallocated after reading !!!
					if(j.ne.0) then
						print *, 'error in allocation initial_mass_star. stat= ', j
						stop
					end if
	
					do j=1,nstar_clump
						read(83,'( 2(1x,i10),1x,i1,10(1x,es12.5) )') IND_star(j), ID_star(j), IS_star(j), xstar(j), ystar(j), zstar(j), vxstar(j), &
						& vystar(j), vzstar(j), mass_star(j), initial_mass_star(j), age_star(j), met_star(j)
					end do
					close(unit=83)
					closest_star1(i,:) = minloc( (xstar(:) - rcar_clump(i,1))**2 + (ystar(:) - rcar_clump(i,2))**2 + (zstar(:) - rcar_clump(i,3))**2 )
					closest_star1(i,1) = ID_star(closest_star1(i,1))
					deallocate( IND_star, IS_star, xstar, ystar, zstar, vxstar, vystar, vzstar, initial_mass_star, met_star )
	
					!!! look for merger trees !!!
					if( track_clumps ) then
						print *, 'YOU SELECTED TO TRACK CLUMP HISTORIES!'
						print *, 'GOOD LUCK!'
						if( snapshot .gt. 1 .and. clump_comment(i) .ne. ' bulge' .and. nstar_group(i,1) .ge. min_star_track ) then
							j = snapshot-1
							delt = 0.0_4
							do while(j .ge. 1 .and. 1000.0_4 * delt .lt. Ndyn_track * td_global)						!!! disc diameter crossing time for clump tracking
								delt = 0.95_4/(((1+redshift(snapshot))/7.0_4)**1.5) - 0.95_4/(((1+redshift(j))/7.0_4)**1.5)	!!! Gyr
								for_merger(:) = 0
								do k=1,ngroup_hist(j)
									if( nstar_hist(j,k) .ge. min_star_track ) then
										nstar_clump2 = 0
										same_stars_mass = 0.0_8
										same_stars_num = 0
										open(unit=84,file=final_clump_stars_filename(j,k))		!!! stars in clump
										read(84,'(i10)') nstar_clump2
										allocate( ID_star2(nstar_clump2), stat=l )			!!! deallocated after merger trees !!!
										if(l.ne.0) then
											print *, 'error in allocation ID_star2. stat= ', l
											stop
										end if
										allocate( mass_star2(nstar_clump2), stat=l )    		!!! deallocated after merger trees !!!
										if(l.ne.0) then
											print *, 'error in allocation mass_star2. stat= ', l
											stop
										end if
										allocate( IND_star2(nstar_clump2), stat=l )			!!! deallocated after reading !!!
										if(l.ne.0) then
											print *, 'error in allocation IND_star2. stat= ', l
											stop
										end if
										allocate( IS_star2(nstar_clump2), stat=l )			!!! deallocated after reading !!!
										if(l.ne.0) then
											print *, 'error in allocation IS_star2. stat= ', l
											stop
										end if
										allocate( xstar2(nstar_clump2), stat=l )		      	!!! deallocated after reading !!!
										if(l.ne.0) then
											print *, 'error in allocation xstar2. stat= ', l
											stop
										end if
										allocate( ystar2(nstar_clump2), stat=l )		      	!!! deallocated after reading !!!
										if(l.ne.0) then
											print *, 'error in allocation ystar2. stat= ', l
											stop
										end if
										allocate( zstar2(nstar_clump2), stat=l )		      	!!! deallocated after reading !!!
										if(l.ne.0) then
											print *, 'error in allocation zstar2. stat= ', l
											stop
										end if
										allocate( vxstar2(nstar_clump2), stat=l )		      	!!! deallocated after reading !!!
										if(l.ne.0) then
											print *, 'error in allocation vxstar2. stat= ', l
											stop
										end if
										allocate( vystar2(nstar_clump2), stat=l )		      	!!! deallocated after reading !!!
										if(l.ne.0) then
											print *, 'error in allocation vystar2. stat= ', l
											stop
										end if
										allocate( vzstar2(nstar_clump2), stat=l )		      	!!! deallocated after reading !!!
										if(l.ne.0) then
											print *, 'error in allocation vzstar2. stat= ', l
											stop
										end if
										allocate( age_star2(nstar_clump2), stat=l )			!!! deallocated after reading !!!
										if(l.ne.0) then
											print *, 'error in allocation age_star2. stat= ', l
											stop
										end if
										allocate( met_star2(nstar_clump2), stat=l )			!!! deallocated after reading !!!
										if(l.ne.0) then
											print *, 'error in allocation met_star2. stat= ', l
											stop
										end if
										allocate( initial_mass_star2(nstar_clump2), stat=l )		!!! deallocated after reading !!!
										if(l.ne.0) then
											print *, 'error in allocation initial_mass_star2. stat= ', l
											stop
										end if
	
										do l=1,nstar_clump2
											read(84,'( 2(1x,i10),1x,i1,10(1x,es12.5) )') IND_star2(l), ID_star2(l), IS_star2(l), xstar2(l), &
											& ystar2(l), zstar2(l), vxstar2(l), vystar2(l), vzstar2(l), mass_star2(l), initial_mass_star2(l), age_star2(l), met_star2(l)
										end do
										close(unit=84)
										deallocate( IND_star2, IS_star2, xstar2, ystar2, zstar2, vxstar2, vystar2, vzstar2, initial_mass_star2, age_star2, met_star2 )
							!!!!!!!! This is definately simpler, but I feel requires MANY more operations and must be slower, despite being vectorized !!!!!!!
!										do l=1,nstar_clump2
!											if( any( ID_star(:) .eq. ID_star2(l) ) ) then
!												same_stars_mass = same_stars_mass + mass_star2(l)
!												same_stars_num = same_stars_num + 1
!											end if
!										end do
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
										n1 = 1
										n2 = 1
										do while( n2 .le. nstar_clump2 )
											do while( ID_star(n1) .gt. ID_star2(n2) .and. n1 .lt. nstar_clump )	!!! stars are written from the linked list in order of decreasing ID
												n1 = n1+1
											end do
											if( n1 .eq. nstar_clump ) then
												if( ID_star(n1) .eq. ID_star2(n2) ) then
													same_stars_mass = same_stars_mass + mass_star2(n2)
													same_stars_num = same_stars_num + 1
												end if
												n2 = nstar_clump2 + 10
											else
												if( ID_star(n1) .eq. ID_star2(n2) ) then
													same_stars_mass = same_stars_mass + mass_star2(n2)
													same_stars_num = same_stars_num + 1
												end if
												n2 = n2+1
											end if
										end do
										go_on = .false.
										if( same_stars_num .ge. min_star_track ) then
											if( by_num ) then
												if( sngl(same_stars_num) .gt. same_num_thresh * sngl(nstar_clump2) ) go_on = .true.
											else if( by_mass ) then
												if( same_stars_mass .gt. same_mass_thresh * clump_mass_hist(j,k) ) go_on = .true.
											end if
										end if
										if( go_on ) then	!!! Old clump donated many stars to new clump --> Found possible parent !!!
											if( clump_id(i) .eq. 0 ) then			!!! Have not yet assigned any parents !!!
												clump_id(i) = id_hist(j,k)
												for_merger(1) = j
												for_merger(2) = k
												Mstars_in(i) = 0.0_8
												Mstars_out(i) = 0.0_8
												Mstars_formed(i) = 0.0_8
!												do l=1,nstar_clump2
!													if( .not. any(ID_star(:) .eq. ID_star2(l)) ) then
!														Mstars_out(i) = Mstars_out(i) + mass_star2(l)
!													end if
!												end do
												n1 = 1
												n2 = 1
												do while( n2 .le. nstar_clump2 )
													if( ID_star2(n2) .gt. ID_star(n1) ) then		!!! no match
														Mstars_out(i) = Mstars_out(i) + mass_star2(n2)
														n2 = n2+1
													elseif( ID_star2(n2) .eq. ID_star(n1) ) then		!!! match
														n2 = n2+1
													elseif( ID_star2(n2) .lt. ID_star(n1) ) then		!!! search for match
														do while( ID_star2(n2) .lt. ID_star(n1) .and. n1 .lt. nstar_clump )
															n1 = n1+1
														end do
														if( ID_star2(n2) .gt. ID_star(n1) ) then	!!! didn't find a match
															Mstars_out(i) = Mstars_out(i) + mass_star2(n2)
															n2 = n2+1
														elseif( ID_star2(n2) .eq. ID_star(n1) ) then	!!! found a match
															n2 = n2+1
														elseif( ID_star2(n2) .lt. ID_star(n1) ) then	!!! n1 must have reached nstar_clump -> no more matches
															Mstars_out(i) = Mstars_out(i) + sum( mass_star2(n2:nstar_clump2) )
															n2 = nstar_clump2 + 10
														end if
													end if
												end do
												Mstars_out(i) = Mstars_out(i) / ( 1.d9 * real(delt,8) )					!!! M_{sun} yr^{-1}
	
!												do l=1,nstar_clump
!													if( .not. any(ID_star2(:) .eq. ID_star(l)) )then
!														if( age_star(l) .le. real(delt,8) ) then
!															Mstars_formed(i) = Mstars_formed(i) + mass_star(l)
!														else
!															Mstars_in(i) = Mstars_in(i) + mass_star(l)
!														end if
!													end if
!												end do
												n1 = 1
												n2 = 1
												do while( n1 .le. nstar_clump )
													if( ID_star(n1) .gt. ID_star2(n2) ) then		!!! no match
														if( age_star(n1) .le. real(delt,8) ) then
															Mstars_formed(i) = Mstars_formed(i) + mass_star(n1)
														else
															Mstars_in(i) = Mstars_in(i) + mass_star(n1)
														end if
														n1 = n1+1
													elseif( ID_star(n1) .eq. ID_star2(n2) ) then		!!! match
														n1 = n1+1
													elseif( ID_star(n1) .lt. ID_star2(n2) ) then		!!! search for match
														do while( ID_star(n1) .lt. ID_star2(n2) .and. n2 .lt. nstar_clump2 )
															n2 = n2+1
														end do
														if( ID_star(n1) .gt. ID_star2(n2) ) then	!!! didn't find a match
															if( age_star(n1) .le. real(delt,8) ) then
																Mstars_formed(i) = Mstars_formed(i) + mass_star(n1)
															else
																Mstars_in(i) = Mstars_in(i) + mass_star(n1)
															end if
															n1 = n1+1
														elseif( ID_star(n1) .eq. ID_star2(n2) ) then	!!! found a match
															n1 = n1+1
														elseif( ID_star(n1) .lt. ID_star2(n2) ) then	!!! n2 must have reached nstar_clump2 -> no more matches
															Mstars_formed(i) = Mstars_formed(i) + sum( mass_star(n1:nstar_clump), &
																& mask = age_star(n1:nstar_clump) .le. real(delt,8) )
															Mstars_in(i)     = Mstars_in(i)     + sum( mass_star(n1:nstar_clump), &
																& mask = age_star(n1:nstar_clump) .gt. real(delt,8) )
															n1 = nstar_clump + 10
														end if
													end if
												end do
												Mstars_in(i)     = Mstars_in(i)     / ( 1.d9 * real(delt,8) )				!!! M_{sun} yr^{-1}
												Mstars_formed(i) = Mstars_formed(i) / ( 1.d9 * real(delt,8) )				!!! M_{sun} yr^{-1}
											else						!!! Have already found a parent --> We have a merger of more than 1 parent !!!
												if( clump_mass_hist(j,k) .gt. clump_mass_hist(for_merger(1),for_merger(2)) ) then
											!!! If the new parent candidate is more massive than the previous !!!
													clump_id(i) = id_hist(j,k)
													for_merger(3) = for_merger(1)
													for_merger(4) = for_merger(2)
													for_merger(1) = j
													for_merger(2) = k
													Mstars_in(i) = 0.0_8
													Mstars_out(i) = 0.0_8
													Mstars_formed(i) = 0.0_8
!													do l=1,nstar_clump2
!														if( .not. any(ID_star(:) .eq. ID_star2(l)) ) then
!															Mstars_out(i) = Mstars_out(i) + mass_star2(l)
!														end if
!													end do
													n1 = 1
													n2 = 1
													do while( n2 .le. nstar_clump2 )
														if( ID_star2(n2) .gt. ID_star(n1) ) then		!!! no match
															Mstars_out(i) = Mstars_out(i) + mass_star2(n2)
															n2 = n2+1
														elseif( ID_star2(n2) .eq. ID_star(n1) ) then		!!! match
															n2 = n2+1
														elseif( ID_star2(n2) .lt. ID_star(n1) ) then		!!! search for match
															do while( ID_star2(n2) .lt. ID_star(n1) .and. n1 .lt. nstar_clump )
																n1 = n1+1
															end do
															if( ID_star2(n2) .gt. ID_star(n1) ) then	!!! didn't find a match
																Mstars_out(i) = Mstars_out(i) + mass_star2(n2)
																n2 = n2+1
															elseif( ID_star2(n2) .eq. ID_star(n1) ) then	!!! found a match
																n2 = n2+1
															elseif( ID_star2(n2) .lt. ID_star(n1) ) then	!!! n1 must have reached nstar_clump -> no match
																Mstars_out(i) = Mstars_out(i) + sum( mass_star2(n2:nstar_clump2) )
																n2 = nstar_clump2 + 10
															end if
														end if
													end do
													Mstars_out(i) = Mstars_out(i) / ( 1.d9 * real(delt,8) )				!!! M_{sun} yr^{-1}
	
!													do l=1,nstar_clump
!														if( .not. any(ID_star2(:) .eq. ID_star(l)) )then
!															if( age_star(l) .le. real(delt,8) ) then
!																Mstars_formed(i) = Mstars_formed(i) + mass_star(l)
!															else
!																Mstars_in(i) = Mstars_in(i) + mass_star(l)
!															end if
!														end if
!													end do
													n1 = 1
													n2 = 1
													do while( n1 .le. nstar_clump )
														if( ID_star(n1) .gt. ID_star2(n2) ) then		!!! no match
															if( age_star(n1) .le. real(delt,8) ) then
																Mstars_formed(i) = Mstars_formed(i) + mass_star(n1)
															else
																Mstars_in(i) = Mstars_in(i) + mass_star(n1)
															end if
															n1 = n1+1
														elseif( ID_star(n1) .eq. ID_star2(n2) ) then		!!! match
															n1 = n1+1
														elseif( ID_star(n1) .lt. ID_star2(n2) ) then		!!! search for match
															do while( ID_star(n1) .lt. ID_star2(n2) .and. n2 .lt. nstar_clump2 )
																n2 = n2+1
															end do
															if( ID_star(n1) .gt. ID_star2(n2) ) then	!!! didn't find a match
																if( age_star(n1) .le. real(delt,8) ) then
																	Mstars_formed(i) = Mstars_formed(i) + mass_star(n1)
																else
																	Mstars_in(i) = Mstars_in(i) + mass_star(n1)
																end if
																n1 = n1+1
															elseif( ID_star(n1) .eq. ID_star2(n2) ) then	!!! found a match
																n1 = n1+1
															elseif( ID_star(n1) .lt. ID_star2(n2) ) then	!!! n1 must have reached nstar_clump -> no match
																Mstars_formed(i) = Mstars_formed(i) + sum( mass_star(n1:nstar_clump), &
																	& mask = age_star(n1:nstar_clump) .le. real(delt,8) )
																Mstars_in(i)     = Mstars_in(i)     + sum( mass_star(n1:nstar_clump), &
																	& mask = age_star(n1:nstar_clump) .gt. real(delt,8) )
																n1 = nstar_clump + 10
															end if
														end if
													end do
													Mstars_in(i)     = Mstars_in(i)     / ( 1.d9 * real(delt,8) )				!!! M_{sun} yr^{-1}
													Mstars_formed(i) = Mstars_formed(i) / ( 1.d9 * real(delt,8) )				!!! M_{sun} yr^{-1}
												else								!!! New candidate is less massive than previous one !!!
													if( for_merger(3) .eq. 0 .and. for_merger(4) .eq.0 ) then	!!! Is there already a merger candidate? !!!
														for_merger(3) = j
														for_merger(4) = k
													else								!!! There is another secondary candidate !!!
														if( clump_mass_hist(j,k) .gt. clump_mass_hist(for_merger(3),for_merger(4)) ) then
														!!! Which secondary candidate is more massive? !!!
															for_merger(3) = j
															for_merger(4) = k
														end if
													end if
												end if
											end if
										end if
										deallocate( ID_star2, mass_star2 )
									end if
								end do
								if(clump_id(i).eq.0) then
									j = j-1
								else					!!! the clump has been given a used ID
									if( i .gt. 1 ) then		!!! Make sure that the ID given to this clump is not occupied by another clump in this same snapshot
										if( any(clump_id(1:i-1) .eq. clump_id(i)) ) then	!!! Have any of the previous clumps in this snapshot been given the same ID?
											print *, '!!!!!!!!!!'
											print *, 'WARNING: 2 CLUMPS TAGGED WITH THE SAME ID'
											bad_clump_id = clump_id(i)
											l = 1
											do while (l .le. i-1)				!!! Which one?
												if(clump_id(l) .eq. bad_clump_id) then
													if( clump_mass(l) .ge. clump_mass(i) ) then	!!! Which is more massive?
														bad_clump1 = l
														bad_clump2 = i
													else
														bad_clump1 = i
														bad_clump2 = l
													end if
													l = i + 5
												else
													l = l + 1
												end if
											end do
											print *, 'initial clump',bad_clump1,'and initial clump',bad_clump2
											print *, 'both have final ID',bad_clump_id
											do l = 1,ngroup_hist(j)
												if(id_hist(j,l) .eq. bad_clump_id) then
													if( clump_merger_hist(j,l) .ne. 0 ) then
														print *, 'The 2 clumps were previously tagged as a merger'
														print *, 'The new IDs are',bad_clump_id,'and',clump_merger_hist(j,l)
														clump_id(bad_clump2) = clump_merger_hist(j,l)
													else
														print *, 'The 2 clumps were 1 clump that split'
														print *, 'Opening a new ID'
														print *, 'The new IDs are',bad_clump_id,'and',ngroup_total + 1
														ngroup_total = ngroup_total + 1
														clump_id(bad_clump2) = ngroup_total
														new(bad_clump2) = 1
														Mstars_out(bad_clump2) = 0.0_8
														Mstars_in(bad_clump2) = 0.0_8
														Mstars_formed(bad_clump2) = 0.0_8
														clump_merger_hist(snapshot,bad_clump2) = 0
														write(clump_comment(bad_clump2),'(a,i5.5)') ' S',clump_id(bad_clump1)
														if(bad_clump2 .gt. bad_clump1) then	!!! the smaller clump, which is now new, is the later one
															for_merger(3) = 0
															for_merger(4) = 0
														end if
													end if
												end if
											end do
											print *, '!!!!!!!!!!'
										else
											if( for_merger(3) .ne. 0 ) then
												write(clump_comment(i),'(a,i5.5)') ' M',id_hist(for_merger(3),for_merger(4))
												clump_merger_hist(snapshot,i) = id_hist(for_merger(3),for_merger(4))
											end if
										end if
									else
										if( for_merger(3) .ne. 0 ) then
											write(clump_comment(i),'(a,i5.5)') ' M',id_hist(for_merger(3),for_merger(4))
											clump_merger_hist(snapshot,i) = id_hist(for_merger(3),for_merger(4))
										end if
									end if
									j = 0
								end if
							end do
						end if
					else
						print *, 'YOU SELECTED TO SKIP TRACKING CLUMP HISTORIES!'
						print *, 'THAT WAS FAST!'
					end if
					deallocate( ID_star, mass_star, age_star )
				end if
				if( track_clumps ) then
					print *, 'found merger tree for',i,'of',ngroup
				end if
	
				if(clump_id(i).eq.0) then
					ngroup_total = ngroup_total + 1
					clump_id(i) = ngroup_total
					new(i) = 1
					Mstars_out(i) = 0.0_8
					Mstars_in(i) = 0.0_8
					Mstars_formed(i) = 0.0_8
				end if
	
				!!! Determine Ex Situ !!!
				if( clump_comment(i) .ne. ' bulge' ) then
					if( mexsitu(i) .gt. 0.5_8 * clump_mass(i) ) then
						if( (clump_dm_dens(i) / dm_back(i)) .gt. dm_thresh) then
							if( abs( sngl(vcyl_clump(i,1)) - vgas_back(i,1) ) .gt. 2.0_4 * sigma_gas_back(i,1) .or. &
							&   abs( sngl(vcyl_clump(i,2)) - vgas_back(i,2) ) .gt. 2.0_4 * sigma_gas_back(i,2) .or. &
							&   abs( sngl(vcyl_clump(i,3)) - vgas_back(i,3) ) .gt. 2.0_4 * sigma_gas_back(i,3) ) then
!							if( ( vcyl_clump(i,2) * sqrt(rcar_clump_disc(i,1)**2 + rcar_clump_disc(i,2)**2) ) / ( sqrt(vcyl_clump(i,1)**2 + vcyl_clump(i,2)**2 + vcyl_clump(i,3)**2) * &
!							& sqrt( rcar_clump(i,1)**2 + rcar_clump(i,2)**2  + rcar_clump(i,3)**2) ) .lt. fj ) then
								exsitu(i) = 123
								if(clump_comment(i).eq.'') then
									write(clump_comment(i),'(a)') ' Es123'
								else
									write(temporary_comment,'(a)') trim(clump_comment(i))
									write(clump_comment(i),'(a,a)') ' Es123', trim(temporary_comment)
								end if
							else
								exsitu(i) = 12
								if(clump_comment(i).eq.'') then
									write(clump_comment(i),'(a)') ' Es12'
								else
									write(temporary_comment,'(a)') trim(clump_comment(i))
									write(clump_comment(i),'(a,a)') ' Es12', trim(temporary_comment)
								end if
							end if
						else
							if( abs( sngl(vcyl_clump(i,1)) - vgas_back(i,1) ) .gt. 2.0_4 * sigma_gas_back(i,1) .or. &
							&   abs( sngl(vcyl_clump(i,2)) - vgas_back(i,2) ) .gt. 2.0_4 * sigma_gas_back(i,2) .or. &
							&   abs( sngl(vcyl_clump(i,3)) - vgas_back(i,3) ) .gt. 2.0_4 * sigma_gas_back(i,3) ) then
!							if( (vcyl_clump(i,2)*sqrt(rcar_clump_disc(i,1)**2 + rcar_clump_disc(i,2)**2)) / ( sqrt(vcyl_clump(i,1)**2 + vcyl_clump(i,2)**2 + vcyl_clump(i,3)**2) * &
!							& sqrt( rcar_clump(i,1)**2 + rcar_clump(i,2)**2  + rcar_clump(i,3)**2) ) .lt. fj ) then
								exsitu(i) = 13
								if(clump_comment(i).eq.'') then
									write(clump_comment(i),'(a)') ' Es13'
								else
									write(temporary_comment,'(a)') trim(clump_comment(i))
									write(clump_comment(i),'(a,a)') ' Es13', trim(temporary_comment)
								end if
							else
								exsitu(i) = 1
								if(clump_comment(i).eq.'') then
									write(clump_comment(i),'(a)') ' Es1'
								else
									write(temporary_comment,'(a)') trim(clump_comment(i))
									write(clump_comment(i),'(a,a)') ' Es1', trim(temporary_comment)
								end if
							end if
						end if
					else
						if( (clump_dm_dens(i) / dm_back(i)) .gt. dm_thresh) then
							if( abs( sngl(vcyl_clump(i,1)) - vgas_back(i,1) ) .gt. 2.0_4 * sigma_gas_back(i,1) .or. &
							&   abs( sngl(vcyl_clump(i,2)) - vgas_back(i,2) ) .gt. 2.0_4 * sigma_gas_back(i,2) .or. &
							&   abs( sngl(vcyl_clump(i,3)) - vgas_back(i,3) ) .gt. 2.0_4 * sigma_gas_back(i,3) ) then
!							if( (vcyl_clump(i,2)*sqrt(rcar_clump_disc(i,1)**2 + rcar_clump_disc(i,2)**2)) / ( sqrt(vcyl_clump(i,1)**2 + vcyl_clump(i,2)**2 + vcyl_clump(i,3)**2) * &
!							& sqrt( rcar_clump(i,1)**2 + rcar_clump(i,2)**2  + rcar_clump(i,3)**2) ) .lt. fj ) then
								exsitu(i) = 23
								if(clump_comment(i).eq.'') then
									write(clump_comment(i),'(a)') ' Es23'
								else
									write(temporary_comment,'(a)') trim(clump_comment(i))
									write(clump_comment(i),'(a,a)') ' Es23', trim(temporary_comment)
								end if
							else
								exsitu(i) = 2
								if(clump_comment(i).eq.'') then
									write(clump_comment(i),'(a)') ' Es2'
								else
									write(temporary_comment,'(a)') trim(clump_comment(i))
									write(clump_comment(i),'(a,a)') ' Es2', trim(temporary_comment)
								end if
							end if
						else
							if( abs( sngl(vcyl_clump(i,1)) - vgas_back(i,1) ) .gt. 2.0_4 * sigma_gas_back(i,1) .or. &
							&   abs( sngl(vcyl_clump(i,2)) - vgas_back(i,2) ) .gt. 2.0_4 * sigma_gas_back(i,2) .or. &
							&   abs( sngl(vcyl_clump(i,3)) - vgas_back(i,3) ) .gt. 2.0_4 * sigma_gas_back(i,3) ) then
!							if( (vcyl_clump(i,2)*sqrt(rcar_clump_disc(i,1)**2 + rcar_clump_disc(i,2)**2)) / ( sqrt(vcyl_clump(i,1)**2 + vcyl_clump(i,2)**2 + vcyl_clump(i,3)**2) * &
!							& sqrt( rcar_clump(i,1)**2 + rcar_clump(i,2)**2  + rcar_clump(i,3)**2) ) .lt. fj ) then
								exsitu(i) = 3
								if(clump_comment(i).eq.'') then
									write(clump_comment(i),'(a)') ' Es3'
								else
									write(temporary_comment,'(a)') trim(clump_comment(i))
									write(clump_comment(i),'(a,a)') ' Es3', trim(temporary_comment)
								end if
							end if
						end if
					end if
				end if
				print *, 'determined ex-situ for',i,'of',ngroup
			end do

			!!! Adjust size of clump history arrays if needed !!!
			n1 = size(id_hist, 1)
			n2 = size(id_hist, 2)
			if(ngroup .gt. n2) then
				print *, 'updating sizes of history arrays'
				print *, n1, n2, ngroup
				call update_history_arrays(n1, n2, ngroup)
			end if

			!!! Adjust clump_stars filenames to match clump IDs !!!
			print *, 'fixing file names to match IDs'
			do i=1,ngroup
				ngroup_hist(snapshot) = ngroup
				id_hist(snapshot,i) = clump_id(i)
				nstar_hist(snapshot,i) = nstar_group(i,1)
				clump_mass_hist(snapshot,i) = clump_mass(i)

				write(final_clump_stars_filename(snapshot,i),'(a,a,a,a,a,i5.5,a)') './',trim(dirname),'/clump_stars/',trim(gal),'_final_clump',clump_id(i),'.out'
				write(temp_name,'(a,a,a,a)') 'cp ',trim(clump_stars_filename(i)),' ',trim(final_clump_stars_filename(snapshot,i))
				print *, i
				print *, trim(temp_name)
				call system(trim(temp_name))
			end do
			deallocate( clump_stars_filename )
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
			print *, 'disc data'
			do i=1,ngroup
				do j=1,3
					if( clump_mass(i) / Mbar_disc(snapshot) .gt. min_mass(j) ) then
						if( exsitu(i) .eq. 2 .or. exsitu(i) .eq. 12 ) then
							nbexsitu_rot(j) = nbexsitu_rot(j) + 1
							mdisc_exsitu_rot(j) = mdisc_exsitu_rot(j) + clump_mass(i)
							SFRdisc_exsitu_rot(j) = SFRdisc_exsitu_rot(j) + clump_SFR(i)
						end if
						if( exsitu(i) .eq. 23 .or. exsitu(i) .eq. 123 ) then
							nbexsitu_nonrot(j) = nbexsitu_nonrot(j) + 1
							mdisc_exsitu_nonrot(j) = mdisc_exsitu_nonrot(j) + clump_mass(i)
							SFRdisc_exsitu_nonrot(j) = SFRdisc_exsitu_nonrot(j) + clump_SFR(i)
						end if
						if( ( exsitu(i) .eq. 0 .or. exsitu(i) .eq. 1 .or. exsitu(i) .eq. 3 .or. exsitu(i) .eq. 13 ) .and. clump_comment(i) .ne.' bulge' ) then
							nbinsitu(j) = nbinsitu(j) + 1
							mdisc_insitu(j) = mdisc_insitu(j) + clump_mass(i)
							SFRdisc_insitu(j) = SFRdisc_insitu(j) + clump_SFR(i)
						end if
						if( clump_comment(i) .eq. ' bulge' ) then
							nbbulge(j) = nbbulge(j) + 1
							mdisc_bulge(j) = mdisc_bulge(j) + clump_mass(i)
							SFRdisc_bulge(j) = SFRdisc_bulge(j) + clump_SFR(i)
						end if
					end if
				end do
			end do
			print *, 'done disc data'
			print *, 'mdisc ex situ rot',      mdisc_exsitu_rot(:)
			mdisc_exsitu_rot(:)      = mdisc_exsitu_rot(:)      / Mbar_disc(snapshot)
			print *, 'mdisc ex situ nonrot',   mdisc_exsitu_nonrot(:)
			mdisc_exsitu_nonrot(:)   = mdisc_exsitu_nonrot(:)   / Mbar_disc(snapshot)
			print *, 'mdisc in situ',          mdisc_insitu(:)
			mdisc_insitu(:)          = mdisc_insitu(:)          / Mbar_disc(snapshot)
			print *, 'mdisc bulge',            mdisc_bulge(:)
			mdisc_bulge(:)           = mdisc_bulge(:)           / Mbar_disc(snapshot)
			print *, 'sfrdisc ex situ rot',    SFRdisc_exsitu_rot(:)
			SFRdisc_exsitu_rot(:)    = SFRdisc_exsitu_rot(:)    / SFR_disc(snapshot)
			print *, 'sfrdisc ex situ nonrot', SFRdisc_exsitu_nonrot(:)
			SFRdisc_exsitu_nonrot(:) = SFRdisc_exsitu_nonrot(:) / SFR_disc(snapshot)
			print *, 'SFRdisc in situ',        SFRdisc_insitu(:)
			SFRdisc_insitu(:)        = SFRdisc_insitu(:)        / SFR_disc(snapshot)
			print *, 'SFRdisc bulge',          SFRdisc_bulge(:)
			SFRdisc_bulge(:)         = SFRdisc_bulge(:)         / SFR_disc(snapshot)
			print *, 'Mbar disc', Mbar_disc(snapshot)
			print *, 'SFR disc', SFR_disc(snapshot)
			print *, 'ngroup',ngroup
			print *, 'nbinsitu',nbinsitu(:)
			print *, 'nbexsitu rot',nbexsitu_rot(:)
			print *, 'nbexsitu nonrot',nbexsitu_nonrot(:)
			write(22,'(3(1x,i3))') ngroup, nbinsitu(1), nbexsitu_rot(1)+nbexsitu_nonrot(1)
	
			print *, 'output'
			!!! Output to clump files !!!
			do i=1,ngroup
				write(32,'(1x,f6.3,1x,i5.5,18es12.2,a11)') redshift(snapshot), clump_id(i), group_rad(i,3), clump_gas_mass(i), clump_star_mass(i), clump_mass(i), &
				& clump_gas_frac(i), clump_dm_frac(i), clump_gas_sig(i), clump_star_sig(i), clump_sig(i), clump_age(i), &
				& clump_gas_met(i), clump_star_met(i), clump_SFR(i), clump_Sig_SFR(i), clump_SSFR(i), clump_tau(i), &
				& clump_dist(i), clump_height(i), trim(clump_comment(i))

				write(33,'(1x,i5.5,1x,f6.3,8(1x,es10.2))') clump_id(i), redshift(snapshot), eta_group(i,1:3), clump_mass(i), group_rad(i,3), spher_mass_test(i,1), &
				& group_rad(i,4), spher_mass_test(i,2)

				if( clump_comment(i) .ne. ' bulge' .and. exsitu(i) .ne. 2 .and. exsitu(i) .ne. 12 .and. exsitu(i) .ne. 23 .and. exsitu(i) .ne. 123  ) then
					j = 26
					k = 29
				end if
				if( clump_comment(i) .ne. ' bulge' .and. exsitu(i) .ne. 0 .and. exsitu(i) .ne. 1 .and. exsitu(i) .ne. 3 .and. exsitu(i) .ne. 13 ) then
					j = 27
					k = 30
				end if
				if( clump_comment(i) .eq. ' bulge' ) then
					j = 28
					k = 31
				end if

				write(j,'(1x,i5.5,1x,f6.3,25(1x,es10.2),1x,i3,20(1x,es10.2))') clump_id(i), redshift(snapshot), group_rad(i,3), & 
				& clump_gas_mass(i), clump_star_mass(i), clump_mass(i), &
				& clump_gas_frac(i), clump_dm_frac(i), clump_gas_sig(i), clump_star_sig(i), clump_sig(i), clump_age(i), &
				& clump_gas_met(i), clump_star_met(i), clump_SFR(i), clump_Sig_SFR(i), clump_SSFR(i), clump_tau(i), &
				& clump_dist(i), clump_height(i), Rdisc(snapshot)*clump_dist(i), Hdisc(snapshot)*clump_height(i), mean_residual(i), eta_group(i,1:3), & 
				& clump_dm_dens(i) / dm_back(i), exsitu(i), tff(i), td(i), td_global, Mgas_in(i,:), Mgas_out(i,:), Mstars_in(i), Mstars_out(i), Mstars_formed(i), alpha_vir(i)

				write(k,'(1x,i5.5,1x,f6.3,25(1x,es10.2),1x,i3,20(1x,es10.2))') clump_id(i), redshift(snapshot), group_rad(i,3), clump_gas_mass(i) / Mgas_disc(snapshot), &
				& clump_star_mass(i) / Mstar_disc(snapshot), clump_mass(i) / Mbar_disc(snapshot), clump_gas_frac(i) / fgas_disc(snapshot), clump_dm_frac(i)/1.0_8, &
				& clump_gas_sig(i) / sigma_gas_disc(snapshot), clump_star_sig(i) / sigma_stars_disc(snapshot), clump_sig(i) / sigma_bar_disc(snapshot), &
				& clump_age(i) / age_disc(snapshot), 10.0_8**( clump_gas_met(i) - metgas_disc(snapshot) ), 10.0_8**( clump_star_met(i) - metstars_disc(snapshot) ), &
				& clump_SFR(i) / SFR_disc(snapshot), clump_Sig_SFR(i) / Sig_SFR_disc(snapshot), clump_SSFR(i) / SSFR_disc(snapshot), clump_tau(i) / tau_disc(snapshot), &
				& clump_dist(i), clump_height(i), Rdisc(snapshot) * clump_dist(i), Hdisc(snapshot) * clump_height(i), mean_residual(i), eta_group(i,1:3), & 
				& clump_dm_dens(i) / dm_back(i), exsitu(i), tff(i), td(i), td_global, Mgas_in(i,:), Mgas_out(i,:), Mstars_in(i), Mstars_out(i), Mstars_formed(i), alpha_vir(i)
			end do
			write(32,*) ''
			print *, 'wrote Matlab friendly files'

			!!! output clump centers to Matlab binary files !!!
			do i=1,Ntracer
				m = 40+i
				write(m) ngroup, clump_id(1:ngroup), new(1:ngroup), exsitu(1:ngroup), rcar_clump_disc(1:ngroup,1), rcar_clump_disc(1:ngroup,2), group_rad(1:ngroup,3), clump_mass(1:ngroup), eta_group(1:ngroup,1)
				close(unit=m)
			end do
			print *, 'wrote binary grid files'

			!!! output to clump catalogue !!!
			allocate( temp_vec(ngroup), stat=i )				!!! deallocated after output !!!
			if(i.ne.0) then
				print *, 'error in allocation temp_vec. stat= ', i
				stop
			end if
			allocate( vec_loc(ngroup), stat=i )				!!! deallocated after output !!!
			if(i.ne.0) then
				print *, 'error in allocation vec_loc. stat= ', i
				stop
			end if
	
			temp_vec(:) = clump_mass(1:ngroup)
			vec_loc(:) = (/ (i,i=1,ngroup) /)
			j=1
			do while(j.le.ngroup)
				temp_loc = maxloc(temp_vec)
				i = vec_loc(temp_loc(1))
				if( clump_comment(i) .ne. ' bulge' .and. exsitu(i) .ne. 2 .and. exsitu(i) .ne. 12 .and. exsitu(i) .ne. 23 .and. exsitu(i) .ne. 123 ) then
					comm = 'Is  '
				elseif( clump_comment(i) .ne. ' bulge' .and. exsitu(i) .ne. 0 .and. exsitu(i) .ne. 1 .and. exsitu(i) .ne. 3 .and. exsitu(i) .ne. 13 ) then
					comm = 'Es  '
				elseif( clump_comment(i) .eq. ' bulge' ) then
					comm = 'Bulg'
				end if
				write(22,'(1x,a4,2(1x,i5.5),1x,i5,1x,i1,3(1x,f7.3),3(1xes12.5),4(1xf7.3),7(1x,es12.5),1x,i5,1x,i10,2(1x,es12.5),2(1x,f7.3),1x,i3)') &
				& comm, clump_id(i), clump_merger_hist(snapshot,i), j, new(i), rcar_clump(i,:), vcar_clump(i,:), group_rad(i,3), rcar_clump_disc(i,:), vcyl_clump(i,:), &
				& clump_gas_mass(i), clump_star_mass(i), clump_dm_mass(i), clump_mass(i) / Mbar_disc(snapshot), ncell(i)-ncell(i-1), nstar_group(i,1), clump_SFR(i), clump_age(i), clump_gas_met(i), &
				& clump_star_met(i), exsitu(i)
				temp_vec(temp_loc(1)) = 0.0_8
	
				write(21,'(1x,i5.5,1x,i10,1x,e12.5)') clump_id(i), closest_star1(i,1), group_rad(i,3)
	
				j = j+1
			end do
			deallocate( temp_vec, vec_loc )
			print *, 'wrote clump catalogues'

			!!! output to disc tables !!!
			write(19,'(9es12.4,13(1x,i3),24(1x,es12.4))') redshift(snapshot), Rdisc(snapshot), Hdisc(snapshot), Rdisc(snapshot) / Rvir(snapshot), &
				& Mgas_disc(snapshot), Mstar_disc(snapshot), Mbar_disc(snapshot), SFR_disc(snapshot), SSFR_disc(snapshot), ngroup, & 
				& nbbulge(:), nbinsitu(:), nbexsitu_rot(:), nbexsitu_nonrot(:), mdisc_bulge(:), mdisc_insitu(:), mdisc_exsitu_rot(:), mdisc_exsitu_nonrot(:), & 
				& SFRdisc_bulge(:), SFRdisc_insitu(:), SFRdisc_exsitu_rot(:), SFRdisc_exsitu_nonrot(:)
			print *, 'wrote clump-disc data'
		else
			print *, 'Mass killed all the clumps :-('

			write(22,'(3(1x,i3))') 0, 0, 0

			write(21,'(1x,i5.5,1x,i10,1x,e12.5)') 0, 0, 0.0_8

			do i=1,Ntracer
				m = 40+i
				write(m) 0
				close(unit=m)
			end do

			write(19,'(9es12.4,13(1x,i3),24(1x,es12.4))') redshift(snapshot), Rdisc(snapshot), Hdisc(snapshot), Rdisc(snapshot) / Rvir(snapshot), &
				& Mgas_disc(snapshot), Mstar_disc(snapshot), Mbar_disc(snapshot), SFR_disc(snapshot), SSFR_disc(snapshot), 0, & 
				& 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, & 
				& 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4
		end if
	end subroutine clump_prop
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine update_history_arrays(m1, m2, m3)
		implicit none
		integer,intent(in) :: m1, m2, m3
		integer :: i, j, k
		integer,allocatable :: id_hist2(:,:), nstar_hist2(:,:), clump_merger_hist2(:,:)
		real(8),allocatable :: clump_mass_hist2(:,:)
		character(len=512),allocatable :: final_clump_stars_filename2(:,:)

		allocate( id_hist2(m1,m2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation id_hist2. stat= ', i
			stop
		end if
		allocate( nstar_hist2(m1,m2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation nstar_hist2. stat= ', i
			stop
		end if
		allocate( clump_merger_hist2(m1,m2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_merger_hist2. stat= ', i
			stop
		end if
		allocate( clump_mass_hist2(m1,m2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_mass_hist2. stat= ', i
			stop
		end if
		allocate( final_clump_stars_filename2(m1,m2), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation final_clump_stars_filename2. stat= ', i
			stop
		end if
		id_hist2(1:m1,1:m2) = id_hist(1:m1,1:m2)
		nstar_hist2(1:m1,1:m2) = nstar_hist(1:m1,1:m2)
		clump_merger_hist2(1:m1,1:m2) = clump_merger_hist(1:m1,1:m2)
		clump_mass_hist2(1:m1,1:m2) = clump_mass_hist(1:m1,1:m2)
		do j=1,m1
			do k=1,m2
				write(final_clump_stars_filename2(j,k),'(a)') trim(final_clump_stars_filename(j,k))
			end do
		end do
		deallocate( id_hist, nstar_hist, clump_merger_hist, clump_mass_hist, final_clump_stars_filename )

		allocate( id_hist(m1,m3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation id_hist_2. stat= ', i
			stop
		end if
		allocate( nstar_hist(m1,m3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation nstar_hist_2. stat= ', i
			stop
		end if
		allocate( clump_merger_hist(m1,m3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_merger_hist_2. stat= ', i
			stop
		end if
		allocate( clump_mass_hist(m1,m3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_mass_hist_2. stat= ', i
			stop
		end if
		allocate( final_clump_stars_filename(m1,m3), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation final_clump_stars_filename_2. stat= ', i
			stop
		end if
		id_hist(1:m1,1:m2) = id_hist2(1:m1,1:m2)
		nstar_hist(1:m1,1:m2) = nstar_hist2(1:m1,1:m2)
		clump_merger_hist(1:m1,1:m2) = clump_merger_hist2(1:m1,1:m2)
		clump_mass_hist(1:m1,1:m2) = clump_mass_hist2(1:m1,1:m2)
		do j=1,m1
			do k=1,m2
				write(final_clump_stars_filename(j,k),'(a)') trim(final_clump_stars_filename2(j,k))
			end do
		end do
		id_hist(1:m1, m2+1:m3) = 0
		nstar_hist(1:m1, m2+1:m3) = 0
		clump_merger_hist(1:m1, m2+1:m3) = 0
		clump_mass_hist(1:m1, m2+1:m3) = 0.0_8
		final_clump_stars_filename(1:m1, m2+1:m3) = ''

		deallocate( id_hist2, nstar_hist2, clump_merger_hist2, clump_mass_hist2, final_clump_stars_filename2 )
	end subroutine update_history_arrays
end module clump_properties

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program main
!#include "adef.h"
use parameters
use globvar
use read_binary
use splitter
use clump_finder
use clump_properties
!use omp_lib

	implicit none
	character(len=20),allocatable :: gal_name(:)                                !!! This is the name of the galaxy you want to search for clumps. Can have either 3 characters (MW3) or 4 characters (VL01, SFG1, MW10)
	character(len=4) :: Rvir_string
	character(len=20) :: snap_tag, snap_tag2
	character(len=256) :: filename, path_name, input_arg
	integer :: h, i, j, k, l, m, n, Nsnapshot, Nsnapshot2, Nsimulation, snap_i, snap_f
	real(4) :: aexp2, Rv4, wide1

	real(4),allocatable :: aexp_min(:), aexp_max(:)
	integer,allocatable :: appen(:)

	real(4) :: Ltemp1(3), Ltemp2(3), Ltemp3(3), Rtemp, Htemp, rctemp, xctemp, yctemp, zctemp, vcar_clump_temp(3), clump_cen_temp(3), vcyl_clump_temp(3), SFR_temp, age_temp, met_temp(2)
	real(8) :: clump_gas_mass_temp, clump_star_mass_temp, clump_dm_mass_temp, clump_normalized_mass_temp
	integer :: nistemp, nestemp, estemp, ncell_temp, nstar_temp
	character(len=8) :: comm_temp
	logical :: file_exist

	!!!!!!!!!! Define tracers !!!!!!!!!!
	allocate( temp_thresh(Ntracer), age_thresh(Ntracer), tracer_name(Ntracer) )
	temp_thresh(1) = 1.5e4	!!! Kelvin
	age_thresh(1) = 0.1	!!! Gyr
	write(tracer_name(1),'(a)') 'Halpha'

	temp_thresh(2) = 1.0	!!! Kelvin
	age_thresh(2) = 15.0	!!! Gyr
	write(tracer_name(2),'(a)') 'stars'

!!!!!!!!!! Get smoothing scale and residual threshold !!!!!!!!!!
	if(iargc().ge.1) then
		call getarg(1,input_arg)
		write(filename,'(a)') trim(input_arg)
		if(iargc().ge.2) then
			call getarg(2,input_arg)
			read(input_arg,*) dens_thresh
			if(iargc().ge.3) then
				call getarg(3,input_arg)
				read(input_arg,*) wide1
			else
				wide1 = 2500.0_4
			end if
		else
			dens_thresh = 10.0_4
			wide1 = 2500.0_4
		end if
	else
		write(filename,'(a)') './input_output/clump_finder_input.dat'
		dens_thresh = 10.0_4
		wide1 = 2500.0_4
	end if
	print *, 'input filename'
	print *, trim(filename)
	print *, 'dens_thresh, wide_FWHM [pc]'
	print *, dens_thresh, wide1

!!!!!!!!!! What simulations will we be looking at today? !!!!!!!!!!
	open(unit=15,file=filename,form='formatted')
	read(15,*) Nsimulation
	print *, Nsimulation
	allocate( gal_name(Nsimulation), stat=i ) 	!!! Deallocated at end of program
	if(i.ne.0) then
		print *, 'error in allocation gal_name. stat= ', i
		stop
	end if
	allocate( aexp_min(Nsimulation), stat=i ) 	!!! Deallocated at end of program
	if(i.ne.0) then
		print *, 'error in allocation aexp_min. stat= ', i
		stop
	end if
	allocate( aexp_max(Nsimulation), stat=i ) 	!!! Deallocated at end of program
	if(i.ne.0) then
		print *, 'error in allocation aexp_max. stat= ', i
		stop
	end if
	allocate( appen(Nsimulation), stat=i )		!!! Deallocated at end of program
	if(i.ne.0) then
		print *, 'error in allocation appen. stat= ', i
		stop
	end if

	do i=1,Nsimulation
		read(15,*) gal_name(i), aexp_min(i), aexp_max(i), appen(i)
	end do
	close(unit=15)

	!!!!!!!!!! Loop over all simulations !!!!!!!!!!
	do j=1,Nsimulation
		print *, 'Simulation:', gal_name(j)
		print *, 'aexp_min:', aexp_min(j)
		print *, 'aexp_max:', aexp_max(j)
		print *, 'apend?:', appen(j)

		write(dirname,'(a,a,i2.2,a,f4.2)') trim(gal_name(j)),'_thresh_',floor(dens_thresh),'_wide_FWHM_',wide1/1000.0_4
		write(filename,'(a,a)') 'mkdir -p ./',trim(dirname)
		call system(trim(filename))

		Nsnapshot = 0
		Nsnapshot2 = 0
		!!!!!!!!!! How many snapshots in the simulation? !!!!!!!!!!
		write(filename,'(a,a,a)') '/BIGDATA/nirm/galaxy_inputs/',trim(gal_name(j)),'/galaxy_catalogue/Nir_halo_cat.txt'
		open(unit=16,file=filename,form='formatted')
		read(16,'(1x,i3)') Nsnapshot
		print *, Nsnapshot

		!!!!!!!!!! Get array sizes !!!!!!!!!!
		allocate( aexp(Nsnapshot), Ngas(Nsnapshot), Ndm(Nsnapshot), Nstars(Nsnapshot) )          !!! deallocated at end of simulation loop !!!
		write(filename,'(a,a,a)') '/BIGDATA/nirm/galaxy_inputs/',trim(gal_name(j)),'/Nmax.txt'
		open(unit=20,file=filename,form='formatted')
		read(20,'(1x,i3)') Nsnapshot2
		if( Nsnapshot2 .ne. Nsnapshot ) then
			print *, 'DATA Nsnapshot - error with catalogue files. Inconsistent number of snapshots per galaxy'
			print *, 'Nsnapshot',Nsnapshot,'Nsnapshot2',Nsnapshot2
			stop
		end if
		do i=1,Nsnapshot
			read(20,'(1x,f5.3,3(1x,i))') aexp(i), Ngas(i), Nstars(i), Ndm(i)
		end do
		close(unit=20)
		print *, Ngas(Nsnapshot), Nstars(Nsnapshot), Ndm(Nsnapshot)
		call allocate_global(Nsnapshot)

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
		open(unit=17,file=filename,form='formatted')
		read(17,'(1x,i3)') Nsnapshot2
		if( Nsnapshot2 .ne. Nsnapshot ) then
			print *, 'DATA Nsnapshot - error with catalogue files. Inconsistent number of snapshots per galaxy'
			print *, 'Nsnapshot',Nsnapshot,'Nsnapshot2',Nsnapshot2
			stop
		end if
		do i=1,Nsnapshot
			read(17,'(1x,f5.3,10(1x,e12.5),2(1x,f7.3),7(1x,es12.5),2(1x,f7.3))') aexp2, rcom(1,i), rcom(2,i), rcom(3,i), vcom(1,i), vcom(2,i), vcom(3,i), &
			& Ldisc(1,i), Ldisc(2,i), Ldisc(3,i), Lmag(i), Rdisc(i), Hdisc(i), Mgas_disc(i), Mcold_disc(i), Mstar_disc(i), M_Es_star_disc(i), Mdm_disc(i), SFR_disc(i), age_disc(i), metgas_disc(i), metstars_disc(i)
			if( aexp2 .ne. aexp(i) ) then
				print *, 'DATA aexp - error with catalogue files. Inconsistent expansion factor'
				print *, 'snapshot #',i,'aexp',aexp(i),'aexp2',aexp2
				stop
			end if
		end do
		close(unit=17)
		if(h_equals_r) then
			Hdisc2(:) = Rdisc(:)
		else
			Hdisc2(:) = Hdisc(:)
		end if
		Mbar_disc(:) = Mstar_disc(:) + Mgas_disc(:)
		fgas_disc(:) = Mgas_disc(:) / Mbar_disc(:)
		sigma_bar_disc(:) = Mbar_disc(:) / real(pi * ( (1000.0_4 * Rdisc(:))**2 ), 8)			!!! M_sun pc^-2
		sigma_stars_disc(:) = Mstar_disc(:) / real(pi * ( (1000.0_4 * Rdisc(:))**2 ), 8)		!!! M_sun pc^-2
		sigma_gas_disc(:) = Mgas_disc(:) / real(pi * ( (1000.0_4 * Rdisc(:))**2 ), 8)			!!! M_sun pc^-2
		Sig_SFR_disc(:) = SFR_disc(:) / real(pi * ( Rdisc(:)**2 ), 8)					!!! M_sun yr^-1 kpc^-2
		SSFR_disc(:) = 1.d9 * SFR_disc(:) / Mstar_disc(:)						!!! Gyr^-1
		tau_disc(:)  = 1.d9 * SFR_disc(:) / Mgas_disc(:)						!!! Gyr^-1

		insitustars(:) = 0
		final_clump_stars_filename(:,:) = ''
		ngroup_hist(:) = 0
		id_hist(:,:) = 0
		clump_merger_hist(:,:) = 0
		nstar_hist(:,:) = 0
		clump_mass_hist(:,:) = 0.0_8
		ngroup_total = 0

		k = 1
		do while(k.le.Nsnapshot)
			if( aexp(k) .lt. aexp_min(j) ) then
				k = k + 1
			else
				snap_i = k
				k = 2*Nsnapshot + 5
			end if
		end do
		k = snap_i
		do while(k.le.Nsnapshot)
			if( aexp(k) .lt. aexp_max(j) .and. k .lt. Nsnapshot ) then
				k = k + 1
			else
				snap_f = k
				k = 2*Nsnapshot + 5
			end if
		end do
		print *, 'snap_i:', snap_i
		print *, 'snap_f:', snap_f
		write(path_name,'(a,a,a,a)') '/BIGDATA/nirm/SIMULATIONS/',trim(gal_name(j)),'/',trim(gal_name(j))

		if(appen(j) .eq. 1 .and. snap_i .gt.1 ) then
			print *, '---------------------------------'
			print *, 'appending: reading old clump data'
			print *, 'number of snapshots to read',snap_i-1
		!!! The important thing is to get the 6 history arrays, including final_clump_stars_filename, and the insitustars array !!!
			do k=1,snap_i-1
				print *, gal_name(j), k, snap_i-1, aexp(k)
				print *, 'Ngas, Nstars, Ndm'
				print *, Ngas(k), Nstars(k), Ndm(k)
				Rv4 = 4.0_4 * Rvir(k)
				If (Rv4 .ge. 1000) then
					write(Rvir_string,'(i4.4)') int(Rv4)
				else
					write(Rvir_string,'(i3.3,a1)') int(Rv4),'.'
				end if
				write(snap_tag,'(a4,a1,f5.3,a4)') trim(Rvir_string),'a',aexp(k),'.dat'
				write(dm_file_name(k),'(a,a,a)') trim(path_name),'_D',trim(snap_tag)
				write(gas_file_name(k),'(a,a,a)') trim(path_name),'_GZ',trim(snap_tag)
				write(stars_z_file_name(k),'(a,a,a)') trim(path_name),'_SZ',trim(snap_tag)
				write(stars_file_name(k),'(a,a,a)') trim(path_name),'_S',trim(snap_tag)
				write(snap_tag2,'(a1,f5.3)') 'a',aexp(k)
				print *, trim(gas_file_name(k))
				print *, trim(stars_z_file_name(k))
				print *, trim(stars_file_name(k))
				print *, trim(dm_file_name(k))
				print *, trim(snap_tag),'','', trim(snap_tag2)

				print *, 'collecting insitu stars data'
				call allocate_stars( k, 1 )
				deallocate( age_stars, mass_stars, vzstars, vystars, vxstars )

				n = (Nstars(k) - mod(Nstars(k),10)) / 10			!!! This just helps keep track of where I am in the loop
				do i=1,Nstars(k)
					if ( mod(i,n) .eq. 0 ) then
						print*, 'i of Nstars',i,'of',Nstars(k)
					end if
					call split0( sngl(xstars(i)), sngl(ystars(i)), sngl(zstars(i)), (/ 0.0_4, 0.0_4, 0.0_4 /) )
					if(insitustars(idstars(i)) .eq. 0) then
						if( abs(zprime(1)) .le. disc_dim*Hdisc(k) .and. rprime(1) .le. disc_dim*Rdisc(k) ) then
							insitustars(idstars(i)) = 1
						else
							insitustars(idstars(i)) = 2
						end if
					end if
				end do
				deallocate( zstars, ystars, xstars, idstars )

				print *, 'collecting clump history data'
				write(filename,'(a,a,a,f5.3,a)') './',trim(dirname),'/clump_catalogue/Nir_clump_cat_a',aexp(k),'.txt'
				inquire(file=filename, exist=file_exist)
				if(file_exist) then
					open(unit=24,file=filename,form='formatted')
					read(24,'(3(1x,e12.5))') Ltemp1(:)
					read(24,'(3(1x,e12.5))') Ltemp2(:)
					read(24,'(3(1x,e12.5))') Ltemp3(:)
					read(24,'(2(1x,f7.3))') Rtemp, Htemp
					read(24,'(3(1x,i3))') ngroup_hist(k), nistemp, nestemp

					print *, 'nclump(',k,')',ngroup_hist(k)

!		                        print ('(3(1x,e12.5))'), Ltemp1(:)
!        		                print ('(3(1x,e12.5))'), Ltemp2(:)
!        		                print ('(3(1x,e12.5))'), Ltemp3(:)
!        		                print ('(2(1x,f7.3))'), Rtemp, Htemp
!        		                print ('(3(1x,i3))'), ngroup_hist(k), nistemp, nestemp
					if( ngroup_hist(k) .gt. 0 ) then
						do i=1,ngroup_hist(k)
							read(24,'(1x,a4,2(1x,i5.5),1x,i5,1x,i1,3(1x,f7.3),3(1xes12.5),4(1xf7.3),7(1x,es12.5),1x,i5,1x,i10,2(1x,es12.5),2(1x,f7.3),1x,i3)') &
							& comm_temp, id_hist(k,i), clump_merger_hist(k,i), m, h, xctemp, yctemp, zctemp, vcar_clump_temp(:), rctemp, clump_cen_temp(:), vcyl_clump_temp(:), &
							& clump_gas_mass_temp, clump_star_mass_temp, clump_dm_mass_temp, clump_normalized_mass_temp, ncell_temp, nstar_hist(k,i), &
							& SFR_temp, age_temp, met_temp(:), estemp
!				                        print('(1x,a8,1x,i2,3(1x,f7.3),3(1xes12.5),4(1xf7.3),7(1x,es12.5),1x,i5,1x,i10,2(1x,es12.5),2(1x,f7.3),1x,i3)'), &
!							& comm, l, xclump(i), yclump(i), zclump(i), vcar_clump(:), rclump(i), rcar_clump_disc(:), vcyl_clump(:), &
!							& clump_gas_mass_temp, clump_star_mass_temp, clump_dm_mass_temp, clump_normalized_mass_temp, ncell_temp, nstar_temp, &
!							& SFR_temp, age_temp, met_temp(:), estemp

							clump_mass_hist(k,i) = clump_gas_mass_temp + clump_star_mass_temp
							write(final_clump_stars_filename(k,i),'(a,a,a,f5.3,a,i5.5,a)') './',trim(dirname),'/clump_stars/a',aexp(k),'_final_clump',id_hist(k,i),'.out'

							print *, 'clump id', id_hist(k,i)
							print *, 'clump merger', clump_merger_hist(k,i)
							print *, 'nstar',nstar_hist(k,i)
							print *, 'clump mass', clump_mass_hist(k,i)
							print *, trim(final_clump_stars_filename(k,i))
						end do
					end if
					close(unit=24)
				end if
				print *, ''
			end do
			ngroup_total = maxval(id_hist(1:snap_i-1,:))
			print *, 'nclump_total:',ngroup_total
			print *, '---------------------------------'
		end if
		call open_files( j )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Loop over all snapshots !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		max_ngroup = 0
		do k=snap_i,snap_f
			if( Rdisc(k) .lt. 2.0_4 * (wide1 / 1000.0_4) ) then
				wide = 0.5_4 * Rdisc(k) * 1000.0_4
			else
				wide = wide1
			end if

			box_size_r = max(disc_dim * Rdisc(k), min_box_size)
			box_size_h = max(disc_dim * Hdisc2(k), min_box_size)

			print *, gal_name(j), k, snap_f, aexp(k)
			print *, 'Rdisc', Rdisc(k), 'Hdisc', Hdisc(k),'kpc'
			print *, 'box size is +/-', box_size_r, box_size_h,'kpc'
			print *, 'wide FWHM is', wide/1000.0_4,'kpc'

			print *, 'Ngas, Nstars, Ndm'
			print *, Ngas(k), Nstars(k), Ndm(k)
			Rv4 = 4.0_4 * Rvir(k)
			If (Rv4 .ge. 1000) then
				write(Rvir_string,'(i4.4)') int(Rv4)
			else
				write(Rvir_string,'(i3.3,a1)') int(Rv4),'.'
			end if
			write(snap_tag,'(a4,a1,f5.3,a4)') trim(Rvir_string),'a',aexp(k),'.dat'
			write(dm_file_name(k),'(a,a,a)') trim(path_name),'_D',trim(snap_tag)
			write(gas_file_name(k),'(a,a,a)') trim(path_name),'_GZ',trim(snap_tag)
			write(stars_z_file_name(k),'(a,a,a)') trim(path_name),'_SZ',trim(snap_tag)
			write(stars_file_name(k),'(a,a,a)') trim(path_name),'_S',trim(snap_tag)
			write(snap_tag2,'(a1,f5.3)') 'a',aexp(k)
			print *, trim(gas_file_name(k))
			print *, trim(stars_z_file_name(k))
			print *, trim(stars_file_name(k))
			print *, trim(dm_file_name(k))
			print *, trim(snap_tag),'','', trim(snap_tag2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Allocate gas and star data and compute grid size !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			call allocate_gas( k )
			deallocate( SNI_gas )
			call grid_size( box_size_r * 1000.0_4, box_size_h * 1000.0_4, 1.5_4 * wide, snap_tag2, k )	!!! Input dimensions in pc
			box_size_r = box_size_r / sngl(rad_correc)
			box_size_h = box_size_h / sngl(rad_correc)

			call allocate_stars( k, 3 )
			deallocate( SNI_stars )

			!!! used for smaller lists, only of stars within the galaxy !!!
			allocate( gas_list(Ngas(k)), stat=i )							!!! deallocated in 'deallocate all' at the end of the snapshot loop !!!
			if(i.ne.0) then
				print *, 'error in allocation of gas_list. stat= ', i
				stop
			end if
			allocate( star_list(Nstars(k)), stat=i )						!!! deallocated in 'clump_stars' after stars are written to files !!!
			if(i.ne.0) then
				print *, 'error in allocation of star_list. stat= ', i
				stop
			end if
			gas_list(:) = 0
			Ngas_list = 0
			star_list(:) = 0
			Nstars_list = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Open files, make grids and find clumps !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			print *, ''
			do i=1,Ntracer
				m = 40+i
				write(filename,'(a,a,a,a,a,a,a)') './',trim(dirname),'/binary_grid_outputs/face_on_surface_densities_',trim(tracer_name(i)),'_',trim(snap_tag2),'.bin'
				open(unit=m,file=filename,form='unformatted')
			end do
			Nhigh_diff1 = 0
			Nhigh_diff2 = 0
			do i=1,Ntracer
				print *, 'tracer', i, trim(tracer_name(i))
				print *, 'going in:'
				print *, 'Nhigh_diff1',Nhigh_diff1,'Nhigh_diff2',Nhigh_diff2
				call create_grid( k, temp_thresh(i), age_thresh(i), i )
				call gaussian_smoothing( narrow/res, wide/res, k, i, snap_tag2 )

				if( i .gt. 1 ) then
					print *, 'number of high residual cells in tracer',i-1,'was',Nhigh_diff2
					if( Nhigh_diff2 .gt. 0 ) then
						do h=1,Nhigh_diff2
							n = ceiling( sngl(high_diff_pos2(h)) / nR**2 )
							m = ceiling( sngl(high_diff_pos2(h) - (nR**2) * (n-1)) / nR )
							l = high_diff_pos2(h) - (nR**2) * (n-1) - nR * (m-1)
							if( high_diff_pos2(h) .ne. (n-1)*(nR*nR) + (m-1)*nR + l ) then
								print *, 'YOU FUCKED UP WITH HIGH_DIFF_POS2!'
								stop
							end if
							if( smoothed_density_diff(l,m,n) .gt. 0.0_4 .or. high_diff_res2(h) .gt. 0.0_4 ) then		!!! One of the tracers, at least, had high density
								smoothed_density_diff(l,m,n) = max( abs(smoothed_density_diff(l,m,n)), abs(high_diff_res2(h)) )
							else												!!! They both had low density
								smoothed_density_diff(l,m,n) = min( smoothed_density_diff(l,m,n), high_diff_res2(h) )
							end if							
						end do
						deallocate( high_diff_pos2, high_diff_res2 )
						Nhigh_diff2 = 0
					end if
				end if
				if( i .lt. Ntracer ) then
					allocate(high_diff_pos1(nR*nR*nH), stat=m )
					if(m.ne.0) then
						print *, 'error in allocation high_diff_pos1. stat= ', m
						stop
					end if
					allocate(high_diff_res1(nR*nR*nH), stat=m )
					if(m.ne.0) then
						print *, 'error in allocation high_diff_pos1. stat= ', m
						stop
					end if
					high_diff_pos1(:) = 0
					Nhigh_diff1 = 0
					do n=1,nH
						do m=1,nR
							do l=1,nR
								if(abs(smoothed_density_diff(l,m,n)) .ge. dens_thresh) then				!!! To capture high residual cells with low density as well
									Nhigh_diff1 = Nhigh_diff1 + 1
									high_diff_pos1(Nhigh_diff1) = (n-1)*(nR*nR) + (m-1)*nR + l
									high_diff_res1(Nhigh_diff1) = smoothed_density_diff(l,m,n)			!!! Keep information about the density itself --> high / low
								end if
							end do
						end do
					end do
					deallocate( smoothed_density_diff )
					print *, 'number of high residual cells in tracer',i,'is',Nhigh_diff1
					if( Nhigh_diff1 .gt. 0 ) then
						Nhigh_diff2 = Nhigh_diff1
						allocate(high_diff_pos2(Nhigh_diff2), stat=m )
						if(m.ne.0) then
							print *, 'error in allocation high_diff_pos1. stat= ', m
							stop
						end if
						allocate(high_diff_res2(Nhigh_diff2), stat=m )
						if(m.ne.0) then
							print *, 'error in allocation high_diff_pos1. stat= ', m
							stop
						end if
						high_diff_pos2 = high_diff_pos1(1:Nhigh_diff2)
						high_diff_res2 = high_diff_res1(1:Nhigh_diff2)
						deallocate( high_diff_pos1, high_diff_res1 )
						Nhigh_diff1 = 0
					end if
				end if
			end do
			if( Ntracer .gt. 1 .and. output_residual) then
			!!! Output 3D residuals to Matlab !!!
				m = 40 + 3*Ntracer + 1
				write(filename,'(a,a,a,a,a)') './',trim(dirname),'/binary_grid_outputs/subtracted_density_combined_',trim(snap_tag2),'.bin'
				open(unit=m,file=filename,form='unformatted')
				write(m) nR, nR, nH, box_size_r, box_size_r, box_size_h, abs(smoothed_density_diff)
				close(unit=m)
			end if
			call group_finder()

			write(filename,'(a,a,a,a,a)') './',trim(dirname),'/clump_catalogue/Nir_simplified_clump_cat_',trim(snap_tag2),'.txt'
			open(unit=21,file=filename)
			write(filename,'(a,a,a,a,a)') './',trim(dirname),'/clump_catalogue/Nir_clump_cat_',trim(snap_tag2),'.txt'
			open(unit=22,file=filename)
			if( ngroup .eq. 0 ) then
				print *, 'no clumps at all'

				saxis3(:) = Ldisc(:,k)
				call axes(saxis1,saxis2,saxis3)
				write(22,'(3(1x,e12.5))') saxis1
				write(22,'(3(1x,e12.5))') saxis2
				write(22,'(3(1x,e12.5))') saxis3
				write(22,'(2(1x,f7.3))') Rdisc(k),Hdisc(k)
				write(22,'(3(1x,i3))') 0, 0, 0

				write(21,'(1x,i5.5,1x,i10,1x,e12.5)') 0, 0, 0.0_8

				write(18,'(1x,f5.3,10(1x,e12.5),2(1x,f7.3),2(1x,i3),6(1x,es12.5),2(1x,f7.3))') aexp(k), rcom(:,k), vcom(:,k), Ldisc(:,k), Lmag(k), &
				& Rdisc(k), Hdisc(k), 0, 0, Mgas_disc(k), Mcold_disc(k), Mstar_disc(k), M_Es_star_disc(k), Mdm_disc(k), SFR_disc(k), age_disc(k), metgas_disc(k), metstars_disc(k)

				do i=1,Ntracer
					m = 40+i
					write(m) 0
				end do

				write(19,'(9es12.4,13(1x,i3),24(1x,es12.4))') redshift(k), Rdisc(k), Hdisc(k), Rdisc(k) / Rvir(k), Mgas_disc(k), Mstar_disc(k), Mbar_disc(k), &		!!! disc_sizes_and_clumps
					& SFR_disc(k), SSFR_disc(k), 0, & 
					& 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, & 
					& 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4, 0.0_4
				call deallocate_clump_prop()
			else
				print *, 'calculating clump properties'
				call clump_prop( snap_tag2, k )
				print *, 'finshed clump prop'
				write(18,'(1x,f5.3,10(1x,e12.5),2(1x,f7.3),2(1x,i3),7(1x,es12.5),2(1x,f7.3))') aexp(k), rcom(:,k), vcom(:,k), Ldisc(:,k), Lmag(k), &	!!! disc_cat2
				& Rdisc(k), Hdisc(k), nbinsitu(1), nbexsitu_rot(1)+nbexsitu_nonrot(1), Mgas_disc(k), Mcold_disc(k), Mstar_disc(k), M_Es_star_disc(k), Mdm_disc(k), &
				& SFR_disc(k), age_disc(k), metgas_disc(k), metstars_disc(k)
				print *, 'wrote to disc cat2'
				call deallocate_clump_prop()
			end if
			close(unit=22)
			close(unit=21)

			call deallocate_gas_stars_dm()
			call deallocate_primes()
			print *, ''
		end do
		print *, 'done'

		call deallocate_global()
		call close_files()
		deallocate( aexp, Ngas, Ndm, Nstars )
		print *, ''

		print *, 'Deleting old clump-star files'
		write(filename,'(a,a,a)') 'rm ./',trim(dirname),'/clump_stars/*_initial_*'
		print *, trim(filename)
		call system(trim(filename))
		print *, 'max_ngroup', max_ngroup
	end do
	deallocate( gal_name, aexp_min, aexp_max, appen )
	deallocate( temp_thresh, age_thresh, tracer_name )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contains
	subroutine allocate_global(N)
!!!!!!!!!! Allocate global arrays per simulation !!!!!!!!!!
!!!!!!!!!! Deallocated in 'deallocate_global' !!!!!!!!!!
		implicit none
		integer,intent(in) :: N
		integer :: i

		print *, 'starting allocate global'
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
		allocate( Hdisc2(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Hdisc2. stat= ', i
			stop
		end if
		allocate( Mbar_disc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Mbar_disc. stat= ', i
			stop
		end if
		allocate( fgas_disc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation fgas_disc. stat= ', i
			stop
		end if
		allocate( sigma_bar_disc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation sigma_bar_disc. stat= ', i
			stop
		end if
		allocate( sigma_stars_disc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation sigma_stars_disc. stat= ', i
			stop
		end if
		allocate( sigma_gas_disc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation sigma_gas_disc. stat= ', i
			stop
		end if
		allocate( Sig_SFR_disc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation Sig_SFR_disc. stat= ', i
			stop
		end if
		allocate( SSFR_disc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation SSFR_disc. stat= ', i
			stop
		end if
		allocate( tau_disc(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation tau_disc. stat= ', i
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
		allocate( stars_z_file_name(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation stars_z_file_name. stat= ', i
			stop
		end if
		allocate( stars_file_name(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation stars_file_name. stat= ', i
			stop
		end if
		allocate( insitustars(Nstars(N)), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation insitustars. stat= ', i
			stop
		end if
		allocate( final_clump_stars_filename(N,Nclumps_max), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation final_clump_stars_filename. stat= ', i
			stop
		end if
		allocate( ngroup_hist(N), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation ngroup_hist. stat= ', i
			stop
		end if
		allocate( id_hist(N,Nclumps_max), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation id_hist. stat= ', i
			stop
		end if
		allocate( clump_merger_hist(N,Nclumps_max), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_merger_hist. stat= ', i
			stop
		end if
		allocate( nstar_hist(N,Nclumps_max), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation nstar_hist. stat= ', i
			stop
		end if
		allocate( clump_mass_hist(N,Nclumps_max), stat=i )
		if(i.ne.0) then
			print *, 'error in allocation clump_mass_hist. stat= ', i
			stop
		end if
		print *, 'finished allocate global'

	end subroutine allocate_global
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine deallocate_global()
!!! Deallocates all arrays from alloctae_global !!!
		implicit none

		if(allocated(ngroup_hist)) deallocate(ngroup_hist)
		if(allocated(id_hist)) deallocate(id_hist)
		if(allocated(clump_merger_hist)) deallocate(clump_merger_hist)
		if(allocated(nstar_hist)) deallocate(nstar_hist)
		if(allocated(clump_mass_hist)) deallocate(clump_mass_hist)
		if(allocated(final_clump_stars_filename)) deallocate(final_clump_stars_filename)
		if(allocated(insitustars)) deallocate(insitustars)

		if(allocated(dm_file_name)) deallocate(dm_file_name)
		if(allocated(gas_file_name)) deallocate(gas_file_name)
		if(allocated(stars_z_file_name)) deallocate(stars_z_file_name)
		if(allocated(stars_file_name)) deallocate(stars_file_name)

		if(allocated(Mbar_disc)) deallocate(Mbar_disc)
		if(allocated(fgas_disc)) deallocate(fgas_disc)
		if(allocated(sigma_bar_disc)) deallocate(sigma_bar_disc)
		if(allocated(sigma_stars_disc)) deallocate(sigma_stars_disc)
		if(allocated(sigma_gas_disc)) deallocate(sigma_gas_disc)
		if(allocated(Sig_SFR_disc)) deallocate(Sig_SFR_disc)
		if(allocated(SSFR_disc)) deallocate(SSFR_disc)
		if(allocated(tau_disc)) deallocate(tau_disc)
		if(allocated(rcom)) deallocate(rcom)
		if(allocated(vcom)) deallocate(vcom)
		if(allocated(Ldisc)) deallocate(Ldisc)
		if(allocated(Lmag)) deallocate(Lmag)
		if(allocated(Rdisc)) deallocate(Rdisc)
		if(allocated(Hdisc)) deallocate(Hdisc)
		if(allocated(Hdisc2)) deallocate(Hdisc2)
		if(allocated(Mgas_disc)) deallocate(Mgas_disc)
		if(allocated(Mcold_disc)) deallocate(Mcold_disc)
		if(allocated(Mstar_disc)) deallocate(Mstar_disc)
		if(allocated(M_Es_star_disc)) deallocate(M_Es_star_disc)
		if(allocated(Mdm_disc)) deallocate(Mdm_disc)
		if(allocated(SFR_disc)) deallocate(SFR_disc)
		if(allocated(age_disc)) deallocate(age_disc)
		if(allocated(metgas_disc)) deallocate(metgas_disc)
		if(allocated(metstars_disc)) deallocate(metstars_disc)

		if(allocated(Rvir)) deallocate(Rvir)
		if(allocated(Mvir)) deallocate(Mvir)
		if(allocated(Vvir)) deallocate(Vvir)
		if(allocated(redshift)) deallocate(redshift)

	end subroutine deallocate_global
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine open_files( sim_num )
		implicit none
		integer,intent(in) :: sim_num

		if( appen(sim_num) .eq. 0 ) then
			write(filename,'(a,a,a)') 'rm -r ./',trim(dirname),'/*'
			call system(trim(filename))

			write(filename,'(a,a,a)') 'mkdir -p ./',trim(dirname),'/galaxy_catalogue'
			call system(trim(filename))

			write(filename,'(a,a,a)') 'mkdir -p ./',trim(dirname),'/binary_grid_outputs'
			call system(trim(filename))

			write(filename,'(a,a,a)') 'mkdir -p ./',trim(dirname),'/clump_catalogue'
			call system(trim(filename))

			write(filename,'(a,a,a)') 'mkdir -p ./',trim(dirname),'/clump_stars'
			call system(trim(filename))

			write(filename,'(a,a,a)') 'mkdir -p ./',trim(dirname),'/extra_clump_parameters'
			call system(trim(filename))

			write(filename,'(a,a,a)') 'mkdir -p ./',trim(dirname),'/Matlab_friendly_data'
			call system(trim(filename))

			write(filename,'(a,a,a)') './',trim(dirname),'/galaxy_catalogue/Nir_disc_cat2.txt'
			open(unit=18,file=filename,form='formatted')
			write(18,'(1x,i3)') snap_f - snap_i + 1

			write(filename,'(a,a,a)') './',trim(dirname),'/galaxy_catalogue/disc_sizes_and_clumps_HcSp.out'
			open(unit=19,file=filename,form='formatted')

			write(filename,'(a,a,a)') './',trim(dirname),'/galaxy_catalogue/too_large_grids.out'
			open(unit=23,file=filename)

			write(filename,'(a,a,a)') './',trim(dirname),'/Matlab_friendly_data/in_situ.out'
			open(unit=26,file=filename,form='formatted')

			write(filename,'(a,a,a)') './',trim(dirname),'/Matlab_friendly_data/ex_situ.out'
			open(unit=27,file=filename,form='formatted')

			write(filename,'(a,a,a)') './',trim(dirname),'/Matlab_friendly_data/bulge.out'
			open(unit=28,file=filename,form='formatted')

			write(filename,'(a,a,a)') './',trim(dirname),'/Matlab_friendly_data/normalized_in_situ.out'
			open(unit=29,file=filename,form='formatted')

			write(filename,'(a,a,a)') './',trim(dirname),'/Matlab_friendly_data/normalized_ex_situ.out'
			open(unit=30,file=filename,form='formatted')

			write(filename,'(a,a,a)') './',trim(dirname),'/Matlab_friendly_data/normalized_bulge.out'
			open(unit=31,file=filename,form='formatted')

			write(filename,'(a,a,a)') './',trim(dirname),'/extra_clump_parameters/clump_table.out'
			open(unit=32,file=filename,form='formatted')
			write(32,*) 'UNITS:'
			write(32,*) '[mass] = M_sun'
			write(32,*) '[age] = Myr'
			write(32,*) '[dist] = kpc'
			write(32,*) '[metalicity] = log([O/H]) + 12'
			write(32,*) '[Sigma] = M_sun / pc^2'
			write(32,*) '[SFR] = M_sun / yr'
			write(32,*) '[sig_SFR] = M_sun / yr / kpc'
			write(32,*) '[SSFR] = 1 / Gyr '
			write(32,*) '[tau] = 1 / Gyr 	tau = SFR/M_gas = gas consumption time'
			write(32,*) 'gal   ','id    ','Rad         ','Mg          ','Ms          ','Mb          ','fg          ','fdm         ','Sigg        ','Sigs        ','Sigb        ','age         ',&
			& 'zg          ','zs          ','SFR          ','Sig_SFR     ','SSFR        ','tau         ','r/Rd        ','z/Hd     ','spec','com'
			write(32,*) ''

			write(filename,'(a,a,a)') './',trim(dirname),'/Matlab_friendly_data/mass_comparison.out'
			open(unit=33,file=filename,form='formatted')
		else
			write(filename,'(a,a,a)') './',trim(dirname),'/galaxy_catalogue/Nir_disc_cat2.txt'
			open(unit=18,file=filename,form='formatted',STATUS='OLD',POSITION='APPEND')

			write(filename,'(a,a,a)') './',trim(dirname),'/galaxy_catalogue/disc_sizes_and_clumps_HcSp.out'
			open(unit=19,file=filename,form='formatted',STATUS='OLD',POSITION='APPEND')

			write(filename,'(a,a,a)') './',trim(dirname),'/galaxy_catalogue/too_large_grids.out'
			open(unit=23,file=filename,STATUS='OLD',POSITION='APPEND')

			write(filename,'(a,a,a)') './',trim(dirname),'/Matlab_friendly_data/in_situ.out'
			open(unit=26,file=filename,form='formatted',STATUS='OLD',POSITION='APPEND')

			write(filename,'(a,a,a)') './',trim(dirname),'/Matlab_friendly_data/ex_situ.out'
			open(unit=27,file=filename,form='formatted',STATUS='OLD',POSITION='APPEND')

			write(filename,'(a,a,a)') './',trim(dirname),'/Matlab_friendly_data/bulge.out'
			open(unit=28,file=filename,form='formatted',STATUS='OLD',POSITION='APPEND')

			write(filename,'(a,a,a)') './',trim(dirname),'/Matlab_friendly_data/normalized_in_situ.out'
			open(unit=29,file=filename,form='formatted',STATUS='OLD',POSITION='APPEND')

			write(filename,'(a,a,a)') './',trim(dirname),'/Matlab_friendly_data/normalized_ex_situ.out'
			open(unit=30,file=filename,form='formatted',STATUS='OLD',POSITION='APPEND')

			write(filename,'(a,a,a)') './',trim(dirname),'/Matlab_friendly_data/normalized_bulge.out'
			open(unit=31,file=filename,form='formatted',STATUS='OLD',POSITION='APPEND')

			write(filename,'(a,a,a)') './',trim(dirname),'/extra_clump_parameters/clump_table.out'
			open(unit=32,file=filename,form='formatted',STATUS='OLD',POSITION='APPEND')

			write(filename,'(a,a,a)') './',trim(dirname),'/Matlab_friendly_data/mass_comparison.out'
			open(unit=33,file=filename,form='formatted',STATUS='OLD',POSITION='APPEND')
		end if

	end subroutine open_files
!___________________________________________________________________________________________________________________________________________________________________________________________________________________

	subroutine close_files()
		implicit none

		close(unit=18)
		close(unit=19)
		close(unit=23)
		close(unit=26)
		close(unit=27)
		close(unit=28)
		close(unit=29)
		close(unit=30)
		close(unit=31)
		close(unit=32)
		close(unit=33)

	end subroutine close_files

end program main


