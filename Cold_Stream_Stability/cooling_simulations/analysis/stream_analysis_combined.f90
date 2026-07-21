module parameters
	implicit none
	integer,parameter :: Nmax = 118000000
	integer,parameter :: Neps_energy = 2
	real(8),parameter :: eps_energy(2) = (/ 0.1_8, 0.01_8 /)	! Internal energy error for "pure" stream fluid; see PaperII
	real(8),parameter :: eps_visual = 0.04	!0.02				! Colour threshold for visual thickness; see PaperII
	real(8),parameter :: pi=3.141592654_8, pi2=2.0_8*pi
	real(8),parameter :: gamma_s = 5.0_8/3.0_8

!!!	Arrays for gas, stars and DM data
!!!	---------------------------------
	character(len=256),allocatable :: ART_file_name(:)
	real(4),allocatable :: cell_size_gas(:), ygas(:), zgas(:), rgas(:), vxgas(:), vrgas(:), vpgas(:), density_gas(:), pressure_gas(:), colour_gas(:)
	real(8) :: res
	integer :: Ngas

!!!	Stream data
!!!	---------------------------------
	real(4),allocatable :: tsnap(:)
	real(8),allocatable :: stream_momentum(:,:), Hb(:,:), Hs(:,:), stream_radius(:,:), clumping(:,:), pure_frac(:,:), Vc(:,:)
	real(8),allocatable :: stream_masses(:,:), stream_volumes(:,:), turb_vel(:,:), clumping_dense_cold(:,:), energies(:,:), energy_fluxes(:,:), R_fluxes(:)
	real(8) :: dR_flux
	integer :: Nflux
	real(8) :: Rs_init, Vs_init, Vb_init, Vc_init, Cs_init, Cb_init, Rhos_init, Rhob_init, Press_init, Presb_init, Temps_init, Tempb_init, delta_init

end module parameters
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module read_binary
use parameters
	implicit none
contains
	subroutine data_gas(filename)
		implicit none
		character (*),intent (in) :: filename
		real(4) :: xgas, vygas, vzgas
	
	    open ( 12 , file = filename, form = 'unformatted', convert = 'big_endian' )
		cell_size_gas(:) = 0.0_4
		ygas(:) = 0.0_4
		zgas(:) = 0.0_4
		rgas(:) = 0.0_4
		vxgas(:) = 0.0_4
		vrgas(:) = 0.0_4
		vpgas(:) = 0.0_4
		density_gas(:) = 0.0_4
		pressure_gas(:) = 0.0_4
		colour_gas(:) = 0.0_4
		Ngas=1

		DO WHILE (Ngas.lt.Nmax)
			read (12,end=6) cell_size_gas(Ngas), xgas, ygas(Ngas), zgas(Ngas), vxgas(Ngas), vygas, vzgas,&
			& density_gas(Ngas), pressure_gas(Ngas), colour_gas(Ngas)
			rgas(Ngas)  = sqrt( (ygas(Ngas)-0.5_4)**2 + (zgas(Ngas)-0.5_4)**2 )
			vrgas(Ngas) = ( vygas*(ygas(Ngas)-0.5_4)  + vzgas*(zgas(Ngas)-0.5_4) ) / max(rgas(Ngas), 1.e-6)
			vpgas(Ngas) = ( vzgas*(ygas(Ngas)-0.5_4)  - vygas*(zgas(Ngas)-0.5_4) ) / max(rgas(Ngas), 1.e-6)
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
		print *, 'Ngas', Ngas

	end subroutine data_gas
end module read_binary
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module stream_properties
use parameters
use read_binary
	implicit none
contains
	subroutine analysis(snap_num)

		implicit none
		integer,intent(in) :: snap_num
		character(len=256) :: filename, format_string
		integer :: i, j, k, kk, l, m, Nbin, Nbin_PDF, Nmask, Ntrack, j_in, j_out
		real(8) :: Rmax, vol_gas, mass_gas, mass_stream, mass_bg, temp_gas, PDF_min, PDF_max, PDF_res, dv, dvr
		real(8) :: mass_stream_tot, mass_bg_tot, mass_cold_tot, mass_dense_tot
		real(8) :: dens_avg_dense_cold(9), dens_squared_avg_dense_cold(6)
		real(8),allocatable :: Rbin(:), Vbin(:,:), Vrbin(:,:), Dbin(:,:), Tbin(:,:), Cbin(:,:), Pbin(:,:), volume_bin(:), sig2_profiles(:,:)
		real(8),allocatable :: PDF_bin(:), density_pdf(:,:), temperature_pdf(:,:)
		real(8),allocatable :: col_thresh_energy(:), pure_mass(:), pure_dens_avg(:), pure_dens_squared_avg(:)

		print *, 'MAIN ROUTINE'
		print *, 'SNAP NUM=', snap_num
		stream_momentum(:,snap_num)     = 0.0_8
		Hb(:,snap_num)                  = 0.0_8
		Hs(:,snap_num)                  = 0.0_8
		stream_radius(:,snap_num)       = 0.0_8
		clumping(:,snap_num)            = 0.0_8
		pure_frac(:,snap_num)           = 0.0_8
		Vc(:,snap_num)                  = 0.0_8
		stream_masses(:,snap_num)       = 0.0_8
		stream_volumes(:,snap_num)      = 0.0_8
		turb_vel(:,snap_num)            = 0.0_8
		clumping_dense_cold(:,snap_num) = 0.0_8
		energies(:,snap_num)            = 0.0_8
		energy_fluxes(:,snap_num)       = 0.0_8

		call data_gas(trim(ART_file_name(snap_num)))

		res = real(minval( cell_size_gas(1:Ngas) ), 8)
		Nbin = nint(0.5_8 / res)
		print *, 'res=',  res
		print *, 'Lmax=', log(1/res)/log(2.0_8)
		print *, 'Nbin=', Nbin

		print *, '!!! rgas', minval(rgas), maxval(rgas), '!!!'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		print *, 'Calculating Profiles'
		allocate( Rbin(Nbin), Vbin(2,Nbin), Vrbin(2,Nbin), Dbin(2,Nbin), Tbin(2,Nbin), Cbin(2,Nbin), Pbin(2,Nbin), volume_bin(Nbin) )
		Rbin = (/ (i,i=1,Nbin) /) * res - res/2.0_8
		Vbin = 0.0_8
		Vrbin = 0.0_8
		Dbin = 0.0_8
		Tbin = 0.0_8
		Cbin = 0.0_8
		Pbin = 0.0_8
		volume_bin = 0.0_8
		Rmax = sngl(Nbin) * res
		print *, 'Rmax=Nbin*res=', Rmax

		Ntrack = (Ngas - mod(Ngas,10)) / 10
		do i=1,Ngas
			if ( i .eq. 1 .or.  mod(i,Ntrack) .eq. 0 ) then
				print*, 'i of Ngas',i,'of',Ngas
			end if
			if( real(rgas(i),8) .le. Rmax ) then
				j = ceiling(real(rgas(i),8)/res)
				j = min(j, Nbin)
				j = max(j, 1)
				vol_gas     = real(cell_size_gas(i)**3, 8)
				mass_gas    = real(density_gas(i), 8) * vol_gas
				mass_stream = real(colour_gas(i), 8)  * mass_gas

				volume_bin(j) = volume_bin(j) + vol_gas

				Cbin(1,j)  = Cbin(1,j)  + real(colour_gas(i),8) * vol_gas
				Cbin(2,j)  = Cbin(2,j)  + real(colour_gas(i),8) * mass_gas
				Dbin(1,j)  = Dbin(1,j)  + mass_gas
				Dbin(2,j)  = Dbin(2,j)  + mass_stream

				Vbin(1,j)  = Vbin(1,j)  + real(vxgas(i),8) * mass_gas
				Vbin(2,j)  = Vbin(2,j)  + real(vxgas(i),8) * mass_stream
				Vrbin(1,j) = Vrbin(1,j) + real(vrgas(i),8) * mass_gas
				Vrbin(2,j) = Vrbin(2,j) + real(vrgas(i),8) * mass_stream

				Tbin(1,j)  = Tbin(1,j)  + real(pressure_gas(i),8) * vol_gas
				Tbin(2,j)  = Tbin(2,j)  + real(pressure_gas(i),8) * real(colour_gas(i),8) * vol_gas

			end if
		end do
		print *, 'Normalizing Profiles'
		write(filename,'(a,F6.4,a)') './stream_analysis/profiles/t',tsnap(snap_num),'.txt'
		open(unit=30,file=filename,form='formatted')
		filename = ''
		do i=1,Nbin
			if( volume_bin(i) .gt. 0.0_8 ) then
				Pbin(1,i) = Tbin(1,i) / volume_bin(i)

				Tbin(1,i) = Tbin(1,i) / Dbin(1,i)
				Tbin(2,i) = Tbin(2,i) / max( Dbin(2,i), 1.d-9 )

				Vbin(1,i) = Vbin(1,i) / Dbin(1,i)
				Vbin(2,i) = Vbin(2,i) / max( Dbin(2,i), 1.d-9 )
				
				Vrbin(1,i) = Vrbin(1,i) / Dbin(1,i)
				Vrbin(2,i) = Vrbin(2,i) / max( Dbin(2,i), 1.d-9 )

				Dbin(2,i) = Dbin(2,i) / max( Cbin(1,i), 1.d-9 )	! This is not the total stream mass divided by the total volume, but rather the weighted-average stream mass among cells that contain stream mass
				Cbin(2,i) = Cbin(2,i) / Dbin(1,i)

				Dbin(1,i) = Dbin(1,i) / volume_bin(i)
				Cbin(1,i) = Cbin(1,i) / volume_bin(i)
			end if
			write(30,'(13(1x,Es12.5))') Rbin(i), Cbin(1,i), Cbin(2,i), Dbin(1,i), Dbin(2,i), Tbin(1,i), Tbin(2,i), Vbin(1,i), Vbin(2,i), Vrbin(1,i), Vrbin(2,i), Pbin(1,i), Pbin(2,i)
		end do
		close(unit=30)
		Rbin = Rbin + res/2.0_8
		if( snap_num .eq. 1 ) then
			print *, 'Initial Conditions'
			j = 1
			do while( j .le. Nbin )
				if( Cbin(1,j) .gt. 0.5_8 ) then
					j = j + 1
				else
					Rs_init = Rbin(j) - res
					j = 2*Nbin + 5
				end if
			end do
			Vs_init = Vbin(1,1)
			Vb_init = Vbin(1,Nbin)
			Cs_init = sqrt( gamma_s * Tbin(1,1) )
			Cb_init = sqrt( gamma_s * Tbin(1,Nbin) )
			Rhos_init = Dbin(1,1)
			Rhob_init = Dbin(1,Nbin)
			Press_init = Pbin(1,1)
			Presb_init = Pbin(1,Nbin)
			Temps_init = Tbin(1,1)
			Tempb_init = Tbin(1,Nbin)
			delta_init = Rhos_init / Rhob_init
			Vc_init = Vs_init * sqrt(delta_init) / ( 1 + sqrt(delta_init) )
			print *, 'Rs', Rs_init, 1.0_8/Rs_init
			print *, 'Ds', Rhos_init
			print *, 'Db', Rhob_init
			print *, 'Vs', Vs_init
			print *, 'Vb', Vb_init
			print *, 'Cs', Cs_init
			print *, 'Cb', Cb_init
			print *, 'Ps', Press_init
			print *, 'Pb', Presb_init
			print *, 'Ms', Vs_init / Cs_init
			print *, 'Mb', Vs_init / Cb_init
			print *, 'delta', delta_init, (Cb_init / Cs_init)**2
			print *, 'Vc', Vc_init
			print *, 'Ts', Temps_init
			print *, 'Tb', Tempb_init
		end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		print *, 'Calculating Hs, Hb, Rs'
		print *, 'Rs_init=', Rs_init
		j = 1
		do while( j .le. Nbin )
			if( Cbin(1,j) .gt. 1.0_8 - eps_visual ) then
				j = j + 1
			else
				Hs(1,snap_num) = max( 0.0_8, Rs_init - Rbin(j) + res )
				j_in = j
				j = 2*Nbin + 5
			end if
		end do
		j = 1
		do while( j .le. Nbin )
			if( Cbin(2,j) .gt. 1.0_8 - eps_visual ) then
				j = j + 1
			else
				Hs(2,snap_num) = max( 0.0_8, Rs_init - Rbin(j) + res )
				j = 2*Nbin + 5
			end if
		end do
		j = Nbin-1
		do while( j .ge. 1 )
			if( Cbin(1,j) .lt. eps_visual ) then
				j = j - 1
			else
				Hb(1,snap_num) = max( 0.0_8, Rbin(j+1) - res - Rs_init )
				j_out = j
				j = 0
			end if
		end do
		j = Nbin-1
		do while( j .ge. 1 )
			if( Cbin(2,j) .lt. eps_visual ) then
				j = j - 1
			else
				Hb(2,snap_num) = max( 0.0_8, Rbin(j+1) - res - Rs_init )
				j = 0
			end if
		end do
		j = 1
		do while( j .le. Nbin )
			if( Cbin(1,j) .gt. 0.5_8 ) then
				j = j + 1
			else
				stream_radius(1,snap_num) = Rbin(j) - res
				j = 2*Nbin + 5
			end if
		end do
		j = 1
		do while( j .le. Nbin )
			if( Cbin(2,j) .gt. 0.5_8 ) then
				j = j + 1
			else
				stream_radius(2,snap_num) = Rbin(j) - res
				j = 2*Nbin + 5
			end if
		end do
		print *, 'Rs0=', Rs_init
		print *, 'Rs/Rs0=', stream_radius(1,snap_num) / Rs_init
		print *, '!!! rgas', minval(rgas), maxval(rgas), '!!!'
		print *, '!!! rgas/Rs0: ', minval(rgas/Rs_init), maxval(rgas/Rs_init), '!!!'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!! Pure fraction and clumping factor !!!
		print *, 'Initializing pure fraction and clumping factor'
		allocate( col_thresh_energy(Neps_energy+1), pure_mass(Neps_energy+1), pure_dens_avg(Neps_energy+1), pure_dens_squared_avg(Neps_energy+1) )
		col_thresh_energy(:) = 0.0_8
		pure_mass(:) = 0.0_8
		pure_dens_avg(:) = 0.0_8
		pure_dens_squared_avg(:) = 0.0_8
		do j=1,Neps_energy
			col_thresh_energy(j+1) = 1.0_8 - eps_energy(j) / delta_init
		end do
		print *, 'N epsilon energy=', Neps_energy
		print *, 'delta=', delta_init
		print *, 'epsilon energy=', eps_energy(:)
		print *, 'col thresh energy=', col_thresh_energy(:)

		!!! PDFs !!!
		print *, 'Initializing PDFs'
		Nbin_PDF = 300
		PDF_min = -3.0_8
		PDF_max = 3.0_8
		PDF_res = (PDF_max - PDF_min)/real(Nbin_PDF,8)
		allocate( PDF_bin(Nbin_PDF), density_pdf(Neps_energy+1,Nbin_PDF), temperature_pdf(Neps_energy+1,Nbin_PDF) )
		PDF_bin(:) = PDF_min + (/ (i,i=1,Nbin_PDF) /) * PDF_res - 0.5_8 * PDF_res
		density_pdf(:,:) = 0.0_8
		temperature_pdf(:,:) = 0.0_8
		Nmask = 0

		!!! Turbulence and dense/cold clumping factor !!!
		print *, 'Initializing turbulence and dense/cold clumping factor'
		allocate( sig2_profiles(8,Nbin) )
		sig2_profiles = 0.0_8
		dens_avg_dense_cold(:) = 0.0_8
		dens_squared_avg_dense_cold(:) = 0.0_8

		!!! Energies and Energy Fluxes !!!
		print *, 'Initializing eneries and energy fluxes'
		if( snap_num .eq. 1 ) then
			R_fluxes(:) = R_fluxes(:) * Rs_init
			dR_flux = dR_flux * Rs_init
			Nflux = size(R_fluxes)
		end if
		print*, 'Nflux=',    Nflux
		print*, 'R_fluxes=', R_fluxes(:)
		print*, 'dR_flux=',  dR_flux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		do i=1,Ngas
			if ( i .eq. 1 .or.  mod(i,Ntrack) .eq. 0 ) then
				print*, 'i of Ngas',i,'of',Ngas
			end if

			!!! Stream and background masses and volumes per cell !!!
			vol_gas     = real(cell_size_gas(i)**3, 8)
			mass_gas    = real(density_gas(i), 8) * vol_gas
			mass_stream = real(colour_gas(i), 8)  * mass_gas
			mass_bg     = mass_gas - mass_stream
			temp_gas    = real(pressure_gas(i), 8) / real(density_gas(i), 8)

			!!! Pure fraction in bins of energy error (see Padnos et al 2018) !!!
			do j=1,Neps_energy+1
				if( real(colour_gas(i), 8) .ge. col_thresh_energy(j) ) then
					pure_mass(j)             = pure_mass(j)             + mass_stream
					pure_dens_avg(j)         = pure_dens_avg(j)         + mass_stream * real(density_gas(i), 8)
					pure_dens_squared_avg(j) = pure_dens_squared_avg(j) + mass_stream * real(density_gas(i)**2, 8)
				end if
			end do

			!!! Masses !!!
			stream_masses(1,snap_num) = stream_masses(1,snap_num) + mass_gas
			stream_masses(2,snap_num) = stream_masses(2,snap_num) + mass_stream
			if( real(density_gas(i),8) .ge. Rhos_init/3.0_8 ) then
				stream_masses(3,snap_num) = stream_masses(3,snap_num) + mass_gas
			end if
			if( temp_gas .le. Temps_init*3.0_8 ) then
				stream_masses(4,snap_num) = stream_masses(4,snap_num) + mass_gas
			end if

			!!! Volumes !!!
			stream_volumes(1,snap_num) = stream_volumes(1,snap_num) + vol_gas
			stream_volumes(2,snap_num) = stream_volumes(2,snap_num) + vol_gas * real(colour_gas(i), 8)
			if( real(density_gas(i),8) .ge. Rhos_init/3.0_8 ) then
				stream_volumes(3,snap_num) = stream_volumes(3,snap_num) + vol_gas
			end if
			if( temp_gas .le. Temps_init*3.0_8 ) then
				stream_volumes(4,snap_num) = stream_volumes(4,snap_num) + vol_gas
			end if

			if( real(rgas(i),8) .le. Rmax ) then
				!!! Convective velocity in cylindrical bin and two orthogonal planar bins (see Mandelker et al 2019) !!!
				if( real(rgas(i),8) .le. Rs_init+Hb(1,snap_num) .and. real(rgas(i),8) .ge. Rs_init-Hs(1,snap_num) ) then
					Vc(1,snap_num) = Vc(1,snap_num) + mass_gas * real(vxgas(i), 8)
					Vc(2,snap_num) = Vc(2,snap_num) + mass_gas
					if( zgas(i) .le. 0.501_4 .and. zgas(i) .ge. 0.499_4 ) then
						Vc(3,snap_num) = Vc(3,snap_num) + mass_gas * real(vxgas(i), 8)
						Vc(4,snap_num) = Vc(4,snap_num) + mass_gas
					end if
					if( ygas(i) .le. 0.501_4 .and. ygas(i) .ge. 0.499_4 ) then
						Vc(5,snap_num) = Vc(5,snap_num) + mass_gas * real(vxgas(i), 8)
						Vc(6,snap_num) = Vc(6,snap_num) + mass_gas
					end if
				end if

				!!! (1) Total momentum profile, (2) center of mass stream velocity, (3) center of mass background velocity !!!
				j = ceiling(real(rgas(i),8)/Rs_init)
				j = min(j, 20)
				j = max(j, 1)
				do kk=j,20
					if(kk .eq. 20) then
						Nmask = Nmask + 1
					end if
					stream_momentum(kk,snap_num) = stream_momentum(kk,snap_num) + mass_gas * real(vxgas(i), 8)
				end do
				stream_momentum(21,snap_num) = stream_momentum(21,snap_num) + mass_stream * real(vxgas(i), 8)
				mass_stream_tot              = mass_stream_tot              + mass_stream
				stream_momentum(22,snap_num) = stream_momentum(22,snap_num) + mass_bg     * real(vxgas(i), 8)
				mass_bg_tot                  = mass_bg_tot                  + mass_bg
				if( temp_gas .le. Temps_init*3.0_8 ) then
					stream_momentum(23,snap_num) = stream_momentum(23,snap_num) + mass_gas * real(vxgas(i), 8)
					mass_cold_tot                = mass_cold_tot                + mass_gas
				end if
				if( real(density_gas(i),8) .ge. Rhos_init/3.0_8 ) then
					stream_momentum(24,snap_num) = stream_momentum(24,snap_num) + mass_gas * real(vxgas(i), 8)
					mass_dense_tot               = mass_dense_tot               + mass_gas
				end if

				!!! Density and temperature PDFs !!!
				do k=1,Neps_energy+1
					if( real(colour_gas(i), 8) .ge. col_thresh_energy(k) ) then
						l = ceiling( (log10(real(density_gas(i),8)/Rhos_init) - PDF_min) / PDF_res )
						l = min(l, Nbin_PDF)
						l = max(l, 1)
						density_pdf(k,l) = density_pdf(k,l) + mass_stream
						m = ceiling( (log10(temp_gas/Temps_init) - PDF_min) / PDF_res )
						m = min(m, Nbin_PDF)
						m = max(m, 1)
						temperature_pdf(k,m) = temperature_pdf(k,m) + mass_stream
					end if
				end do
			end if

			!!! Turbulent velocity from profile, and dense/cold clumping factor !!!
			if( real(rgas(i),8) .le. Rs_init+Hb(1,snap_num) ) then
				j = ceiling(real(rgas(i),8)/res)
				j = min(j, Nbin)
				j = max(j, 1)
				sig2_profiles(1,j) = sig2_profiles(1,j) + mass_gas    * ( (real(vxgas(i),8)-Vbin(1,j))**2 + (real(vrgas(i),8)-Vrbin(1,j))**2 + real(vpgas(i),8)**2 )
				sig2_profiles(2,j) = sig2_profiles(2,j) + mass_stream * ( (real(vxgas(i),8)-Vbin(2,j))**2 + (real(vrgas(i),8)-Vrbin(2,j))**2 + real(vpgas(i),8)**2 )
				sig2_profiles(3,j) = sig2_profiles(3,j) + mass_gas
				sig2_profiles(4,j) = sig2_profiles(4,j) + mass_stream

				if( j .le. j_out ) then
					dens_avg_dense_cold(1)         = dens_avg_dense_cold(1)         + mass_gas * real(density_gas(i), 8)
					dens_squared_avg_dense_cold(1) = dens_squared_avg_dense_cold(1) + mass_gas * real(density_gas(i), 8)**2

					dens_avg_dense_cold(4)         = dens_avg_dense_cold(4)         + vol_gas * real(density_gas(i), 8)
					dens_squared_avg_dense_cold(4) = dens_squared_avg_dense_cold(4) + vol_gas * real(density_gas(i), 8)**2

					dens_avg_dense_cold(7)         = dens_avg_dense_cold(7)         + vol_gas
				end if
				if( real(density_gas(i),8) .ge. Rhos_init/3.0_8 ) then
					sig2_profiles(5,j)  = sig2_profiles(5,j) + mass_gas * ( (real(vxgas(i),8)-Vbin(1,j))**2 + (real(vrgas(i),8)-Vrbin(1,j))**2 + real(vpgas(i),8)**2 )
					sig2_profiles(7,j)  = sig2_profiles(7,j) + mass_gas
					if( j .le. j_out ) then
						dens_avg_dense_cold(2)         = dens_avg_dense_cold(2)         + mass_gas * real(density_gas(i), 8)
						dens_squared_avg_dense_cold(2) = dens_squared_avg_dense_cold(2) + mass_gas * real(density_gas(i), 8)**2

						dens_avg_dense_cold(5)         = dens_avg_dense_cold(5)         + vol_gas * real(density_gas(i), 8)
						dens_squared_avg_dense_cold(5) = dens_squared_avg_dense_cold(5) + vol_gas * real(density_gas(i), 8)**2

						dens_avg_dense_cold(8)         = dens_avg_dense_cold(8)         + vol_gas
					end if
				end if
				if( temp_gas .le. Temps_init*3.0_8 ) then
					sig2_profiles(6,j) = sig2_profiles(6,j) + mass_gas * ( (real(vxgas(i),8)-Vbin(1,j))**2 + (real(vrgas(i),8)-Vrbin(1,j))**2 + real(vpgas(i),8)**2 )
					sig2_profiles(8,j) = sig2_profiles(8,j) + mass_gas
					if( j .le. j_out ) then
						dens_avg_dense_cold(3)         = dens_avg_dense_cold(3)         + mass_gas * real(density_gas(i), 8)
						dens_squared_avg_dense_cold(3) = dens_squared_avg_dense_cold(3) + mass_gas * real(density_gas(i), 8)**2
						
						dens_avg_dense_cold(6)         = dens_avg_dense_cold(6)         + vol_gas * real(density_gas(i), 8)
						dens_squared_avg_dense_cold(6) = dens_squared_avg_dense_cold(6) + vol_gas * real(density_gas(i), 8)**2

						dens_avg_dense_cold(9)         = dens_avg_dense_cold(9)         + vol_gas
					end if
				end if
			end if

			!!! Energies and Energy Fluxes!!!
			j = ceiling(real(rgas(i),8)/res)
			j = min(j, Nbin)
			j = max(j, 1)
			dv  = real(vxgas(i),8)-Vbin(1,j)
			dvr = real(vrgas(i),8)-Vrbin(1,j)
!			print*, 'Vbulk, dv, mass_gas=',Vbulk,dv,mass_gas
			energies(1,snap_num)  = energies(1,snap_num)  + 0.5_8 * mass_gas * Vbin(1,j)**2
			energies(2,snap_num)  = energies(2,snap_num)  + mass_gas * Vbin(1,j) * dv
			energies(3,snap_num)  = energies(3,snap_num)  + 0.5_8 * mass_gas * dv**2
			energies(4,snap_num)  = energies(4,snap_num)  + 0.5_8 * mass_gas * Vrbin(1,j)**2
			energies(5,snap_num)  = energies(5,snap_num)  + mass_gas * Vrbin(1,j) * dvr
			energies(6,snap_num)  = energies(6,snap_num)  + 0.5_8 * mass_gas * dvr**2
			energies(7,snap_num)  = energies(7,snap_num)  + 0.5_8 * mass_gas * real(vpgas(i),8)**2
			energies(8,snap_num)  = energies(8,snap_num)  + ( 1.0_8/(gamma_s-1.0_8) ) * real(colour_gas(i),8) * real(pressure_gas(i),8) * vol_gas
			energies(9,snap_num)  = energies(9,snap_num)  + ( 1.0_8/(gamma_s-1.0_8) ) * (1.0_8 - real(colour_gas(i),8)) * real(pressure_gas(i),8) * vol_gas
!			if( vrgas(i) .gt. 0.0_4 ) then
			j = 1
			do while ( j .le. Nflux )
				if( real(rgas(i),8) .ge. R_fluxes(j) - 0.5_8*dR_flux .and. real(rgas(i),8) .le. R_fluxes(j) + 0.5_8*dR_flux ) then
					energy_fluxes(j,snap_num) = energy_fluxes(j,snap_num) + ( 0.5_8 * mass_gas * ( real(vxgas(i),8)**2 + real(vrgas(i),8)**2 + real(vpgas(i),8)**2 ) + & 
						&  ( 1.0_8/(gamma_s-1.0_8) ) * real(pressure_gas(i),8) * vol_gas ) * real(vrgas(i),8) 
					j = 2*Nflux + 1
				else
					j = j + 1
				end if
			end do
!			end if

		end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!! Pure fraction in bins of energy error (see Padnos et al 2018) !!!
		do j=1,Neps_energy+1
			pure_dens_squared_avg(j) = pure_dens_squared_avg(j) / max(pure_mass(j), 1.d-10)
			pure_dens_avg(j)         = pure_dens_avg(j) / max(pure_mass(j), 1.d-10)
			clumping(j,snap_num)     = pure_dens_squared_avg(j) / max( pure_dens_avg(j)**2, 1.d-10 )
		end do
		do j=1,Neps_energy
			pure_frac(j,snap_num) = pure_mass(j+1) / max(pure_mass(1), 1.d-10)
		end do
		print *, 'pure frac=', pure_frac(:,snap_num)

		!!! Masses !!!
		print '(a, 4(1x,ES12.4))', 'M/M0 = ', stream_masses(1:4,snap_num)/stream_masses(1:4,1)

		!!! Convective velocity in cylindrical bin and two orthogonal planar bin (see Mandelker et al 2019) !!!
		Vc(1,snap_num) = Vc(1,snap_num) / max(Vc(2,snap_num), 1.d-10)
		Vc(3,snap_num) = Vc(3,snap_num) / max(Vc(4,snap_num), 1.d-10)
		Vc(5,snap_num) = Vc(5,snap_num) / max(Vc(6,snap_num), 1.d-10)

		!!! (1) Total momentum profile, (2) center of mass stream velocity, (3) center of mass background velocity !!!
		stream_momentum(22,snap_num) = ( stream_momentum(21,snap_num) + stream_momentum(22,snap_num) ) / ( mass_stream_tot + mass_bg_tot )
		stream_momentum(21,snap_num) = stream_momentum(21,snap_num) / mass_stream_tot
!		stream_momentum(22,snap_num) = stream_momentum(22,snap_num) / mass_bg_tot
		stream_momentum(23,snap_num) = stream_momentum(23,snap_num) / mass_cold_tot
		stream_momentum(24,snap_num) = stream_momentum(24,snap_num) / mass_dense_tot
		print *, 'Ncells in momentum bin 20:', Nmask
		print *, 'momentum conservation', stream_momentum(20,snap_num) / stream_momentum(20,1)

		!!! Density and temperature PDFs !!!
		write(filename,'(a,F6.4,a)') './stream_analysis/PDFs/t',tsnap(snap_num),'.txt'
		open(unit=31,file=filename,form='formatted')
		filename = ''
		do k=1,Neps_energy+1
			density_pdf(k,:)     = density_pdf(k,:)     / max( sum( density_pdf(k,:) * PDF_res ), 1.d-10 )
			temperature_pdf(k,:) = temperature_pdf(k,:) / max( sum( temperature_pdf(k,:) * PDF_res ), 1.d-10 )
		end do
		k = Neps_energy+1
		write(format_string,'(a,i1,a)') '(',1+2*k,'(1x,ES12.5))'
		do i=1,Nbin_PDF
			write(31,trim(format_string)) PDF_bin(i), density_pdf(1:k,i), temperature_pdf(1:k,i)
		end do
		print *, 'PDF sums:', sum(  density_pdf(1,:) * PDF_res ), sum(  temperature_pdf(1,:) * PDF_res )
		close(unit=31)

		!!! Turbulent velocity from profile, and dense/cold clumping factor !!!
		print '(a,1x,I4,1x,a,1x,I4)', 'j_in=', j_in, ', j_out', j_out
		turb_vel(1,snap_num) = sum( sig2_profiles(1,j_in:j_out) ) / max(sum( sig2_profiles(3,j_in:j_out) ), 1.d-10)	! total gas
		turb_vel(2,snap_num) = sum( sig2_profiles(2,j_in:j_out) ) / max(sum( sig2_profiles(4,j_in:j_out) ), 1.d-10)	! stream gas
		turb_vel(3,snap_num) = sum( sig2_profiles(1,1:j_in) )     / max(sum( sig2_profiles(3,1:j_in) ), 1.d-10)		! total gas
		turb_vel(4,snap_num) = sum( sig2_profiles(2,1:j_in) )     / max(sum( sig2_profiles(4,1:j_in) ), 1.d-10)		! stream gas
		turb_vel(5,snap_num) = sum( sig2_profiles(5,j_in:j_out) ) / max(sum( sig2_profiles(7,j_in:j_out) ), 1.d-10)	! dense gas only
		turb_vel(6,snap_num) = sum( sig2_profiles(6,j_in:j_out) ) / max(sum( sig2_profiles(8,j_in:j_out) ), 1.d-10)	! cold gas only

		turb_vel(1:6,snap_num) = sqrt( turb_vel(1:6,snap_num) )

		print '(a,1x,ES12.5,1x,a,1x,ES12.5)', 'hs/Rs=', Hs(1,snap_num)/Rs_init, ', hb/Rs=', Hb(1,snap_num)/Rs_init
		print '(a, 4(1x,ES12.4))', 'turb{1,2}/Vc, turb{3,4}/V = ', turb_vel(1:2,snap_num)/Vc_init, turb_vel(3:4,snap_num)/Vs_init

		dens_avg_dense_cold(1) = dens_avg_dense_cold(1) / max( sum( sig2_profiles(3,1:j_out) ), 1.d-10 )	! mass weighted, all
		dens_avg_dense_cold(2) = dens_avg_dense_cold(2) / max( sum( sig2_profiles(7,1:j_out) ), 1.d-10 )	! mass weighted, dense
		dens_avg_dense_cold(3) = dens_avg_dense_cold(3) / max( sum( sig2_profiles(8,1:j_out) ), 1.d-10 )	! mass weighted, cold
		dens_avg_dense_cold(4) = dens_avg_dense_cold(4) / max( dens_avg_dense_cold(7), 1.d-10 )				! volume weighted, all
		dens_avg_dense_cold(5) = dens_avg_dense_cold(5) / max( dens_avg_dense_cold(8), 1.d-10 )				! volume weighted, dense
		dens_avg_dense_cold(6) = dens_avg_dense_cold(6) / max( dens_avg_dense_cold(9), 1.d-10 )				! volume weighted, cold

		dens_squared_avg_dense_cold(1) = dens_squared_avg_dense_cold(1) / max(sum( sig2_profiles(3,1:j_out) ), 1.d-10)	! mass weighted all
		dens_squared_avg_dense_cold(2) = dens_squared_avg_dense_cold(2) / max(sum( sig2_profiles(7,1:j_out) ), 1.d-10)	! mass weighted dense
		dens_squared_avg_dense_cold(3) = dens_squared_avg_dense_cold(3) / max(sum( sig2_profiles(8,1:j_out) ), 1.d-10)	! mass weighted cold
		dens_squared_avg_dense_cold(4) = dens_squared_avg_dense_cold(4) / max( dens_avg_dense_cold(7), 1.d-10 )			! volume weighted, all
		dens_squared_avg_dense_cold(5) = dens_squared_avg_dense_cold(5) / max( dens_avg_dense_cold(8), 1.d-10 )			! volume weighted, dense
		dens_squared_avg_dense_cold(6) = dens_squared_avg_dense_cold(6) / max( dens_avg_dense_cold(9), 1.d-10 )			! volume weighted, cold

		clumping_dense_cold(1,snap_num) = dens_squared_avg_dense_cold(1) / max((dens_avg_dense_cold(1)**2), 1.d-10)
		clumping_dense_cold(2,snap_num) = dens_squared_avg_dense_cold(2) / max((dens_avg_dense_cold(2)**2), 1.d-10)
		clumping_dense_cold(3,snap_num) = dens_squared_avg_dense_cold(3) / max((dens_avg_dense_cold(3)**2), 1.d-10)
		clumping_dense_cold(4,snap_num) = dens_squared_avg_dense_cold(4) / max((dens_avg_dense_cold(4)**2), 1.d-10)
		clumping_dense_cold(5,snap_num) = dens_squared_avg_dense_cold(5) / max((dens_avg_dense_cold(5)**2), 1.d-10)
		clumping_dense_cold(6,snap_num) = dens_squared_avg_dense_cold(6) / max((dens_avg_dense_cold(6)**2), 1.d-10)

		!!! Energies and Energy Fluxes!!!
		energy_fluxes(:,snap_num) = energy_fluxes(:,snap_num) / dR_flux
		print *, 'Etot=',sum(energies(:,1)),sum(energies(:,snap_num))
		print *, 'Eflux (2,15)Rs=',energy_fluxes(1,snap_num),energy_fluxes(Nflux,snap_num)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		deallocate( Rbin, Vbin, Vrbin, Dbin, Tbin, Cbin, Pbin, volume_bin )
		deallocate( col_thresh_energy, pure_mass, pure_dens_avg, pure_dens_squared_avg )
		deallocate( PDF_bin, density_pdf, temperature_pdf )
		deallocate( sig2_profiles )

	end subroutine analysis
end module stream_properties

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program main
	use parameters
	use read_binary
	use stream_properties
	implicit none
	integer :: i, j, k, Nsnapshot, Nsnap
	character(len=256) :: filename, format_string, input_arg

	if(iargc().ge.1) then
		call getarg(1,input_arg)
		read(input_arg,*) Nsnap
	else
		Nsnap = 1
	end if
	print *, 'Nsnap'
	print *, Nsnap

	open(unit=16,file='./output/time.txt')
	read(16,*) Nsnapshot
	print *, Nsnapshot
	allocate( ART_file_name(Nsnapshot), tsnap(Nsnapshot) )
	do k=1,Nsnapshot
		read(16,'(F6.4)') tsnap(k)
		write(ART_file_name(k),'(a,F6.4,a)') './output/ART_format_t',tsnap(k),'.dat'
	end do
	close(unit=16)
	print *, trim(ART_file_name(1)), tsnap(1)
	print *, trim(ART_file_name(Nsnap)), tsnap(Nsnap)
	call open_files(Nsnap)

	print *, 'alocating gas arrays'
	allocate( cell_size_gas(Nmax), ygas(Nmax), zgas(Nmax), rgas(Nmax), vxgas(Nmax), vrgas(Nmax), vpgas(Nmax), density_gas(Nmax), pressure_gas(Nmax), colour_gas(Nmax) )

	allocate( stream_momentum(24,Nsnapshot), Hb(2,Nsnapshot), Hs(2,Nsnapshot), stream_radius(2,Nsnapshot) )
	stream_momentum(:,:)     = 0.0_8
	Hb(:,:)                  = 0.0_8
	Hs(:,:)                  = 0.0_8
	stream_radius(:,:)       = 0.0_8

	allocate( clumping(Neps_energy+1,Nsnapshot), pure_frac(Neps_energy,Nsnapshot), Vc(6,Nsnapshot) )
	clumping(:,:)            = 0.0_8
	pure_frac(:,:)           = 0.0_8
	Vc(:,:)                  = 0.0_8

	allocate( stream_masses(4,Nsnapshot), stream_volumes(4,Nsnapshot), turb_vel(6,Nsnapshot), clumping_dense_cold(6,Nsnapshot) )
	stream_masses(:,:)       = 0.0_8
	stream_volumes(:,:)      = 0.0_8
	turb_vel(:,:)            = 0.0_8
	clumping_dense_cold(:,:) = 0.0_8

	allocate( energies(9,Nsnapshot), energy_fluxes(6,Nsnapshot), R_fluxes(6) )
	energies(:,:) = 0.0_8
	energy_fluxes(:,:) = 0.0_8
	R_fluxes(1) = 2.0_8	! In units of Rs
	R_fluxes(2) = 3.0_8
	R_fluxes(3) = 4.0_8
	R_fluxes(4) = 5.0_8
	R_fluxes(5) = 10.0_8
	R_fluxes(6) = 15.0_8
	dR_flux     = 0.5_8

	Rs_init = 0.0_8
	Vs_init = 0.0_8
	Vb_init = 0.0_8 
	Cs_init = 0.0_8
	Cb_init = 0.0_8
	Rhos_init = 0.0_8
	Rhob_init = 0.0_8
	Press_init = 0.0_8
	Presb_init = 0.0_8

	print *, ''
	print *, 't=',tsnap(1)
	call analysis( 1 )
	print *, 't=',tsnap(Nsnap)
	call analysis( Nsnap )
	print *, 'Writing Output'

	write(20,'(F6.4,24(1x,ES12.5))') tsnap(Nsnap), stream_momentum(:,Nsnap)

	write(21,'(F6.4,6(1x,ES12.5))') tsnap(Nsnap), Hs(1,Nsnap), Hb(1,Nsnap), stream_radius(1,Nsnap), Hs(2,Nsnap), Hb(2,Nsnap), stream_radius(2,Nsnap)

	write(format_string,'(a,i1,a)') '(F6.4,',Neps_energy+1,'(1x,ES12.5))'
	write(22,trim(format_string)) tsnap(Nsnap), clumping(:,Nsnap)

	write(format_string,'(a,i1,a)') '(F6.4,',Neps_energy,'(1x,ES12.5))'
	write(23,trim(format_string)) tsnap(Nsnap), pure_frac(:,Nsnap)

	write(24,'(F6.4,3(1x,ES12.5))') tsnap(Nsnap), Vc(1,Nsnap), Vc(3,Nsnap), Vc(5,Nsnap)

	write(25,'(F6.4,4(1x,ES12.5))') tsnap(Nsnap), stream_masses(1:4,Nsnap)

	write(43,'(F6.4,4(1x,ES12.5))') tsnap(Nsnap), stream_volumes(1:4,Nsnap)

	write(26,'(F6.4,6(1x,ES12.5))') tsnap(Nsnap), turb_vel(1:6,Nsnap)

	write(27,'(F6.4,6(1x,ES12.5))') tsnap(Nsnap), clumping_dense_cold(1:6,Nsnap)	! mass weighted all, dense, cold; volume weighted all, dense, cold

	write(28,'(F6.4,9(1x,ES12.5))') tsnap(Nsnap), energies(1:9,Nsnap)

	write(29,'(F6.4,6(1x,ES12.5))') tsnap(Nsnap), energy_fluxes(1:6,Nsnap)
	call close_files()
	deallocate( stream_momentum, Hb, Hs, stream_radius, clumping, pure_frac, Vc )
	deallocate( stream_masses, stream_volumes, turb_vel, clumping_dense_cold )
	deallocate( energies, energy_fluxes, R_fluxes )

	deallocate( cell_size_gas ,ygas, zgas, rgas, vxgas, vrgas, vpgas, density_gas, pressure_gas, colour_gas )
	deallocate( ART_file_name, tsnap )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contains
	subroutine open_files(snap_num)
		implicit none
		integer,intent(in) :: snap_num

		print *, 'enter open files'
		write(filename,'(a)') 'mkdir -p ./stream_analysis'
		call system(filename)
		filename = ''

		write(filename,'(a)') 'mkdir -p ./stream_analysis/profiles'
		call system(filename)
		filename = ''

		write(filename,'(a)') 'mkdir -p ./stream_analysis/PDFs'
		call system(filename)
		filename = ''

		write(filename,'(a)') 'mkdir -p ./stream_analysis/momentum'
		call system(filename)
		filename = ''

		write(filename,'(a)') 'mkdir -p ./stream_analysis/thickness'
		call system(filename)
		filename = ''

		write(filename,'(a)') 'mkdir -p ./stream_analysis/clumping'
		call system(filename)
		filename = ''

		write(filename,'(a)') 'mkdir -p ./stream_analysis/pure_frac'
		call system(filename)
		filename = ''

		write(filename,'(a)') 'mkdir -p ./stream_analysis/Vc'
		call system(filename)
		filename = ''

		write(filename,'(a)') 'mkdir -p ./stream_analysis/stream_masses'
		call system(filename)
		filename = ''

		write(filename,'(a)') 'mkdir -p ./stream_analysis/stream_volumes'
		call system(filename)
		filename = ''

		write(filename,'(a)') 'mkdir -p ./stream_analysis/turbulence_from_profile'
		call system(filename)
		filename = ''

		write(filename,'(a)') 'mkdir -p ./stream_analysis/clumping_factor_dense_cold'
		call system(filename)
		filename = ''

		write(filename,'(a)') 'mkdir -p ./stream_analysis/energies_from_profile'
		call system(filename)
		filename = ''

		write(filename,'(a)') 'mkdir -p ./stream_analysis/energy_fluxes'
		call system(filename)
		filename = ''

		write(filename,'(a,I5.5,a)') './stream_analysis/momentum/momentum_',snap_num,'.txt'
		open(unit=20,file=filename,form='formatted')
		filename = ''

		write(filename,'(a,I5.5,a)') './stream_analysis/thickness/thickness_',snap_num,'.txt'
		open(unit=21,file=filename,form='formatted')
		filename = ''

		write(filename,'(a,I5.5,a)') './stream_analysis/clumping/clumping_',snap_num,'.txt'
		open(unit=22,file=filename,form='formatted')
		filename = ''

		write(filename,'(a,I5.5,a)') './stream_analysis/pure_frac/pure_frac_',snap_num,'.txt'
		open(unit=23,file=filename,form='formatted')
		filename = ''

		write(filename,'(a,I5.5,a)') './stream_analysis/Vc/Vc_',snap_num,'.txt'
		open(unit=24,file=filename,form='formatted')
		filename = ''

		write(filename,'(a,I5.5,a)') './stream_analysis/stream_masses/stream_masses_',snap_num,'.txt'
		open(unit=25,file=filename,form='formatted')
		filename = ''

		write(filename,'(a,I5.5,a)') './stream_analysis/turbulence_from_profile/turbulence_from_profile_',snap_num,'.txt'
		open(unit=26,file=filename,form='formatted')
		filename = ''

		write(filename,'(a,I5.5,a)') './stream_analysis/clumping_factor_dense_cold/clumping_factor_dense_cold_',snap_num,'.txt'
		open(unit=27,file=filename,form='formatted')
		filename = ''

		write(filename,'(a,I5.5,a)') './stream_analysis/energies_from_profile/energies_from_profile_',snap_num,'.txt'
		open(unit=28,file=filename,form='formatted')
		filename = ''

		write(filename,'(a,I5.5,a)') './stream_analysis/energy_fluxes/energy_fluxes_',snap_num,'.txt'
		open(unit=29,file=filename,form='formatted')
		filename = ''

		write(filename,'(a,I5.5,a)') './stream_analysis/stream_volumes/stream_volumes_',snap_num,'.txt'
		open(unit=43,file=filename,form='formatted')
		filename = ''

		print *, 'exit open files'
	end subroutine open_files
!________________________________________________________________________________________

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
		close(unit=43)
	end subroutine close_files

end program main

