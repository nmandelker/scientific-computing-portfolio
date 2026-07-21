module parameters
	implicit none
	integer,parameter :: Nmax = 118000000
	real(8),parameter :: eps_visual = 0.04	!0.02		! Colour threshold for visual thickness; see PaperII
	real(8),parameter :: pi=3.141592654_8, pi2=2.0_8*pi

!!!	Simulation parameters which are not output anywhere unfortunately...
!!!	---------------------------------
	real(8),parameter :: gamma_s  = 5.0_8/3.0_8			! Adiabatic index of the gas in the simulations
	real(8),parameter :: Zs_init  = 0.0006				! Initial metalicity in stream     in absolute units (Z_{solar}=0.02)
	real(8),parameter :: Zb_init  = 0.002				! Initial metalicity in background in absolute units (Z_{solar}=0.02)
	real(8),parameter :: redshift = 2.0				! Redshift of Haardt and Madau UVB used in the simulations

!!!	Physical constants
!!!	---------------------------------
	real(8),parameter :: Z_solar    = 0.02				! Solar metalicity in code units
	real(8),parameter :: XH         = 0.76				! Hydrogen mass fraction
	real(8),parameter :: KB         = 1.38062d-16
	real(8),parameter :: mproton    = 1.66d-24
	real(8),parameter :: Myr_in_sec = 3.1536d13
	real(8),parameter :: kpc_in_cm  = 3.086d21
	real(8),parameter :: Msun_in_gr = 1.98892d33

!!!	Arrays for gas data
!!!	---------------------------------
	character(len=256),allocatable :: ART_file_name(:)
	real(4),allocatable :: cell_size_gas(:), density_gas(:), pressure_gas(:), colour_gas(:), met_gas(:), rgas(:)
	real(8) :: res
	integer :: Ngas

!!!	Stream data
!!!	---------------------------------
	real(4),allocatable :: tsnap(:)
	real(8) :: Rs_init, Vs_init, Vb_init, Cs_init, Cb_init, Rhos_init, Rhob_init, Press_init, Presb_init, delta_init, Ts_init, Tb_init

!!!	Cooling table data
!!!	---------------------------------
	integer :: NDens, NTemp
	real(8) :: T2max, T1max
	real(8),allocatable :: luminosity_vec(:,:)
	real(8),allocatable :: dens_tab(:), T2_tab(:), mu_tab(:,:), nspec_tab(:,:,:)
	real(8),allocatable :: cool_tab(:,:), heat_tab(:,:), cool_com_tab(:,:), heat_com_tab(:,:), metal_tab(:,:)
	real(8),allocatable :: cool_prime_tab(:,:), heat_prime_tab(:,:), cool_com_prime_tab(:,:), heat_com_prime_tab(:,:), metal_prime_tab(:,:)

!!! 	Cooling interpolation parameters
!!!	---------------------------------
	integer :: i_nH, i_T2
	real(8) :: dlog_nH, facH, w1H, w2H
	real(8) :: dlog_T2, facT, w1T, w2T
	real(8) :: h, h2, h3, yy, yy2, yy3
	real(8) :: fa, fb, fprimea, fprimeb, alpha, beta, gama
	real(8) :: cool1, cool_prime1, heat1, heat_prime1, cool_com1, cool_com_prime1, heat_com1, heat_com_prime1, metal1, metal_prime1
	real(8) :: mu1_1, mu1_2, mu1, ne1, ne2, nHI1, nHI2
	real(8) :: Lambda_net, Lambda_tot, Lambda_LyA, Lambda_prime

!!!	Units
!!!	---------------------------------
	real(8) :: unit_length_cgs, unit_density_cgs, unit_time_cgs, unit_velocity_cgs, unit_mass_cgs, unit_energy_cgs


end module parameters
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module read_binary
use parameters
	implicit none
contains
	subroutine data_gas(filename)
		implicit none
		character (*),intent (in) :: filename
		real(4) :: xgas, ygas, zgas, vxgas, vygas, vzgas
	
	    open ( 12 , file = filename, form = 'unformatted', convert = 'big_endian' )
		cell_size_gas(:) = 0.0_4
		density_gas(:)   = 0.0_4
		pressure_gas(:)  = 0.0_4
		colour_gas(:)    = 0.0_4
		met_gas(:)       = 0.0_4
		rgas(:)       = 0.0_4
		Ngas=1

		DO WHILE (Ngas.lt.Nmax)
			read (12,end=6) cell_size_gas(Ngas), xgas, ygas, zgas, vxgas, vygas, vzgas,&
			& density_gas(Ngas), pressure_gas(Ngas), colour_gas(Ngas)
			rgas(Ngas) = sqrt( (ygas-0.5_4)**2 + (zgas-0.5_4)**2 )
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

		met_gas(1:Ngas) = colour_gas(1:Ngas) * sngl(Zs_init) + (1.0_4 - colour_gas(1:Ngas)) * sngl(Zb_init)

	end subroutine data_gas
end module read_binary
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module compute_cooling_rates
use parameters
use read_binary
	implicit none
	real(8),allocatable :: dens_vec(:), temp_vec(:), emis_vec(:,:), vol_vec(:), mass_grid(:,:), emis_grid(:,:,:)
	real(8) :: Dmin, Dmax, Tmin, Tmax, d_dens_vec, d_temp_vec
	integer :: Ndens_vec, Ntemp_vec
contains
	subroutine analysis(snap_num)
		implicit none
		integer,intent(in) :: snap_num
		character(len=256) :: filename
		integer :: i, j, k, Nbin, Ntrack, N_bad_T
		real(8) :: q_LyA, h_nu_LyA
		real(8) :: nH_gas, T2_gas, ne_gas, nHI_gas, met_solar_gas, T_gas, vol_gas, mass_gas, Rmax
		real(8),allocatable :: profiles(:,:), Rbin(:), Tbin(:), Ebin(:,:), Dbin(:), Mbin(:), Vbin(:), Cbin(:)
		logical :: file_exist

		print *, 'MAIN ROUTINE:'
		print *, 'SNAP NUM=', snap_num

		write(filename,'(a,I5.5,a)') './stream_analysis/phase_diagrams/mass_weighted_',snap_num,'.dat'
		open(unit=30,file=filename,form='unformatted')
		filename = ''
		write(filename,'(a,I5.5,a)') './stream_analysis/phase_diagrams/luminosity_weighted_',snap_num,'.dat'
		open(unit=31,file=filename,form='unformatted')
		filename = ''
		write(filename,'(a,I5.5,a)') './stream_analysis/phase_diagrams/emissivity_vs_temperature_',snap_num,'.txt'
		open(unit=32,file=filename)
		filename = ''
		write(filename,'(a,I5.5,a)') './stream_analysis/profiles/cooling_',snap_num,'.txt'
		open(unit=33,file=filename)
		filename = ''

		mass_grid(:,:)   = 0.0_8
		emis_grid(:,:,:) = 0.0_8
		emis_vec(:,:)    = 0.0_8
		vol_vec(:)       = 0.0_8

		write(filename,'(a,F6.4,a)') './stream_analysis/profiles/t',tsnap(1),'.txt'
		inquire(file=filename, exist=file_exist)
		if(file_exist) then
			print *, 'Profiles exist for t=0'
			print *, 'Reading ICs from there'
			print *, 'First read in gas data from current snapshot to get cell size and number of bins'

			print *, 't=',tsnap(snap_num)
			call data_gas(trim(ART_file_name(snap_num)))
			res = real(minval( cell_size_gas(1:Ngas) ), 8)
			Nbin = nint(0.5_8 / res)
			print *, 'res=',  res
			print *, 'Lmax=', log(1.0_8/res)/log(2.0_8)
			print *, 'Nbin=', Nbin

			allocate( profiles(13,Nbin) )
			open(unit=17,file=filename)
			do j=1,Nbin
				read(17,'(13(1x,Es12.5))') profiles(1:13,j)
			end do
			close(unit=17)
			
			j = 1
			do while( j .le. Nbin )
				if( profiles(2,j) .gt. 0.5_8 ) then
					j = j + 1
				else
					Rs_init = profiles(1,j) - res
					j = 2*Nbin + 5
				end if
			end do
			Vs_init    = profiles(8,1)
			Vb_init    = profiles(8,Nbin)
			Cs_init    = sqrt( gamma_s * profiles(6,1) )
			Cb_init    = sqrt( gamma_s * profiles(6,Nbin) )
			Rhos_init  = profiles(4,1)
			Rhob_init  = profiles(4,Nbin)
			Ts_init = profiles(6,1)
			Tb_init = profiles(6,Nbin)
			Press_init = profiles(6,1) * profiles(4,1)
			Presb_init = profiles(6,Nbin) * profiles(4,Nbin)
			delta_init = Rhos_init / Rhob_init

			deallocate( profiles )
		else
			print *, 'Profiles do not exist for t=0'
			print *, 'Must compute ICs from scratch'
			print *, 'Read in gas data from snap=1 to generate profiles and compute ICs'

			print *, 't=',tsnap(1)
			call data_gas(trim(ART_file_name(1)))
			res = real(minval( cell_size_gas(1:Ngas) ), 8)
			Nbin = nint(0.5_8 / res)
			print *, 'res=',  res
			print *, 'Lmax=', log(1.0_8/res)/log(2.0_8)
			print *, 'Nbin=', Nbin
			
			print *, 'Calculating Profiles'
			allocate( Rbin(Nbin), Vbin(Nbin), Dbin(Nbin), Tbin(Nbin), Cbin(Nbin) )
			Rbin = (/ (i,i=1,Nbin) /) * res
			Vbin = 0.0_8
			Dbin = 0.0_8
			Tbin = 0.0_8
			Cbin = 0.0_8
			Rmax = real(Nbin,8) * res
			print *, 'Rmax=Nbin*res=', Rmax

			Ntrack = (Ngas - mod(Ngas,10)) / 10
			do i=1,Ngas
				if ( i .eq. 1 .or.  mod(i,Ntrack) .eq. 0 ) then
					print*, 'i of Ngas',i,'of',Ngas
				end if
				if( rgas(i) .le. Rmax ) then
					j = ceiling(rgas(i) / res)
					j = min(j, Nbin)
					j = max(j, 1)
					vol_gas = real(cell_size_gas(i)**3, 8)
					mass_gas = real(density_gas(i), 8) * vol_gas

					Vbin(j) = Vbin(j) + vol_gas

					Cbin(j)  = Cbin(j)  + real(colour_gas(i),8) * vol_gas
					Dbin(j)  = Dbin(j)  + mass_gas
					Tbin(j)  = Tbin(j)  + real(pressure_gas(i),8) * vol_gas
				end if
			end do
			print *, 'Normalizing Profiles'
			do j=1,Nbin
				if( Vbin(j) .gt. 0.0_8 ) then
					Tbin(j) = Tbin(j) / Dbin(j)
					Dbin(j) = Dbin(j) / Vbin(j)
					Cbin(j) = Cbin(j) / Vbin(j)
				end if
			end do
			print *, 'Initial Conditions'
			j = 1
			do while( j .le. Nbin )
				if( Cbin(j) .gt. 0.5_8 ) then
					j = j + 1
				else
					Rs_init = Rbin(j) - res
					j = 2*Nbin + 5
				end if
			end do
			Cs_init = sqrt( gamma_s * Tbin(1) )
			Cb_init = sqrt( gamma_s * Tbin(Nbin) )
			Rhos_init = Dbin(1)
			Rhob_init = Dbin(Nbin)
			Ts_init = Tbin(1)
			Tb_init = Tbin(Nbin)
			Press_init = Ts_init * Rhos_init
			Presb_init = Tb_init * Rhob_init
			delta_init = Rhos_init / Rhob_init

			deallocate( Rbin, Vbin, Dbin, Tbin, Cbin )
			
			print *, 'Now read in gas data from current snapshot'
			print *, 't=',tsnap(snap_num)
			call data_gas(trim(ART_file_name(snap_num)))
			res = real(minval( cell_size_gas(1:Ngas) ), 8)
			Nbin = nint(0.5_8 / res)
			print *, 'res=',  res
			print *, 'Lmax=', log(1.0_8/res)/log(2.0_8)
			print *, 'Nbin=', Nbin
		end if

		allocate( Rbin(Nbin), Tbin(Nbin), Ebin(Nbin,3), Dbin(Nbin), Mbin(Nbin), Vbin(Nbin) )
		Rbin = (/ (i,i=1,Nbin) /) * res - res/2.0_8
		Tbin(:)   = 0.0_8
		Ebin(:,:) = 0.0_8
		Dbin(:)   = 0.0_8
		Mbin(:)   = 0.0_8
		Vbin(:)   = 0.0_8
		Rmax = real(Nbin) * real(res,8)

		Ntrack = (Ngas - mod(Ngas,10)) / 10
		N_bad_T = 0
		do i=1,Ngas
			if ( i .eq. 1 .or.  mod(i,Ntrack) .eq. 0 ) then
				print *, ''
				print*, 'i of Ngas',i,'of',Ngas
			end if

			nH_gas = real(density_gas(i), 8) * unit_density_cgs * XH / mproton 				! Hydrogen number density in cm^{-3}
			if( nH_gas .lt. 1.d-10 .or. nH_gas .ne. nH_gas .or. nH_gas-1. .eq. nH_gas ) then
				print '(a,3(1x,Es12.5))', 'nH_gas problem ', rgas(i)/Rs_init, nH_gas, density_gas(i)
				nH_gas = 1.d-10
			end if

			T2_gas = ( real(pressure_gas(i),8) * mproton / real(density_gas(i), 8) ) * ( unit_energy_cgs / unit_mass_cgs ) / KB		! T/mu in Kelvin
			if( T2_gas .lt. 1.0_8 .or. T2_gas .ne. T2_gas .or. T2_gas-1. .eq. T2_gas ) then
				print '(a,4(1x,Es12.5))', 'T2_gas problem ', rgas(i)/Rs_init, T2_gas, pressure_gas(i), density_gas(i)
				T2_gas = 1.0_8
				N_bad_T = N_bad_T + 1
			end if

			met_solar_gas = met_gas(i) / Z_solar										! Gas metalicity in solar units

			facH = max( log10(nH_gas), dens_tab(1) )
			facH = min( facH, dens_tab(NDens) )
			i_nH = max( nint( (facH-dens_tab(1)) * dlog_nH ) + 1, 1 )
			i_nH = min( i_nH, NDens-1 )
			w1H  = ( dens_tab(i_nH+1) - facH ) * dlog_nH
			w2H  = ( facH - dens_tab(i_nH) ) * dlog_nH

			facT = max( log10(T2_gas), T2_tab(1) )
			facT = min( facT, T2_tab(NTemp) )
			i_T2 = max( nint((facT-T2_tab(1))*dlog_T2) + 1, 1 )
			i_T2 = min( i_T2, NTemp-1 )
			yy   = facT - T2_tab(i_T2)
			yy2  = yy**2
			yy3  = yy**3

			!!! H and He cooling !!!
			fa = cool_tab(i_nH,i_T2  )*w1H + cool_tab(i_nH+1,i_T2  )*w2H
			fb = cool_tab(i_nH,i_T2+1)*w1H + cool_tab(i_nH+1,i_T2+1)*w2H
			fprimea = cool_prime_tab(i_nH,i_T2  )*w1H + cool_prime_tab(i_nH+1,i_T2  )*w2H
			fprimeb = cool_prime_tab(i_nH,i_T2+1)*w1H + cool_prime_tab(i_nH+1,i_T2+1)*w2H
			alpha = fprimea
			beta  = 3.0_8 * (fb-fa) / h2 - ( 2.0_8 * fprimea + fprimeb ) / h
			gama  = ( fprimea + fprimeb ) / h2 - 2.0_8 * ( fb - fa ) / h3
			cool1       = 10.0_8 ** ( fa + alpha*yy + beta*yy2 + gama*yy3 )
			cool_prime1 = ( cool1 / T2_gas ) * (alpha + 2.0_8*beta*yy + 3.0_8*gama*yy2)

			!!! radiative heating !!!
			fa = heat_tab(i_nH,i_T2  )*w1H + heat_tab(i_nH+1,i_T2  )*w2H
			fb = heat_tab(i_nH,i_T2+1)*w1H + heat_tab(i_nH+1,i_T2+1)*w2H
			fprimea = heat_prime_tab(i_nH,i_T2  )*w1H + heat_prime_tab(i_nH+1,i_T2  )*w2H
			fprimeb = heat_prime_tab(i_nH,i_T2+1)*w1H + heat_prime_tab(i_nH+1,i_T2+1)*w2H
			alpha = fprimea
			beta  = 3.0_8 * (fb-fa) / h2 - ( 2.0_8 * fprimea + fprimeb ) / h
			gama  = ( fprimea + fprimeb ) / h2 - 2.0_8 * ( fb - fa ) / h3
			heat1       = 10.0_8 ** ( fa + alpha*yy + beta*yy2 + gama*yy3 )
			heat_prime1 = ( heat1 / T2_gas ) * (alpha + 2.0_8*beta*yy + 3.0_8*gama*yy2)

			!!! Compton cooling !!!
			fa = cool_com_tab(i_nH,i_T2  )*w1H + cool_com_tab(i_nH+1,i_T2  )*w2H
			fb = cool_com_tab(i_nH,i_T2+1)*w1H + cool_com_tab(i_nH+1,i_T2+1)*w2H
			fprimea = cool_com_prime_tab(i_nH,i_T2  )*w1H + cool_com_prime_tab(i_nH+1,i_T2  )*w2H
			fprimeb = cool_com_prime_tab(i_nH,i_T2+1)*w1H + cool_com_prime_tab(i_nH+1,i_T2+1)*w2H
			alpha = fprimea
			beta  = 3.0_8 * (fb-fa) / h2 - ( 2.0_8 * fprimea + fprimeb ) / h
			gama  = ( fprimea + fprimeb ) / h2 - 2.0_8 * ( fb - fa ) / h3
			cool_com1       = 10.0_8 ** ( fa + alpha*yy + beta*yy2 + gama*yy3 )
			cool_com_prime1 = ( cool_com1 / T2_gas ) * (alpha + 2.0_8*beta*yy + 3.0_8*gama*yy2)

			!!! Compton heating !!!
			fa = heat_com_tab(i_nH,i_T2  )*w1H + heat_com_tab(i_nH+1,i_T2  )*w2H
			fb = heat_com_tab(i_nH,i_T2+1)*w1H + heat_com_tab(i_nH+1,i_T2+1)*w2H
			fprimea = heat_com_prime_tab(i_nH,i_T2  )*w1H + heat_com_prime_tab(i_nH+1,i_T2  )*w2H
			fprimeb = heat_com_prime_tab(i_nH,i_T2+1)*w1H + heat_com_prime_tab(i_nH+1,i_T2+1)*w2H
			alpha = fprimea
			beta  = 3.0_8 * (fb-fa) / h2 - ( 2.0_8 * fprimea + fprimeb ) / h
			gama  = ( fprimea + fprimeb ) / h2 - 2.0_8 * ( fb - fa ) / h3
			heat_com1       = 10.0_8 ** ( fa + alpha*yy + beta*yy2 + gama*yy3 )
			heat_com_prime1 = ( heat_com1 / T2_gas ) * (alpha + 2.0_8*beta*yy + 3.0_8*gama*yy2)

			!!! Metal line cooling !!!
			fa = metal_tab(i_nH,i_T2  )*w1H + metal_tab(i_nH+1,i_T2  )*w2H
			fb = metal_tab(i_nH,i_T2+1)*w1H + metal_tab(i_nH+1,i_T2+1)*w2H
			fprimea = metal_prime_tab(i_nH,i_T2  )*w1H + metal_prime_tab(i_nH+1,i_T2  )*w2H
			fprimeb = metal_prime_tab(i_nH,i_T2+1)*w1H + metal_prime_tab(i_nH+1,i_T2+1)*w2H
			alpha = fprimea
			beta  = 3.0_8 * (fb-fa) / h2 - ( 2.0_8 * fprimea + fprimeb ) / h
			gama  = ( fprimea + fprimeb ) / h2 - 2.0_8 * ( fb - fa ) / h3
			metal1       = 10.0_8 ** ( fa + alpha*yy + beta*yy2 + gama*yy3 )
			metal_prime1 = ( metal1 / T2_gas ) * (alpha + 2.0_8*beta*yy + 3.0_8*gama*yy2)

			!!! Total cooling !!!
			Lambda_tot   = cool1 + met_solar_gas * metal1 + cool_com1 / nH_gas
			Lambda_tot   = Lambda_tot * nH_gas**2		!!! net cooling rate (emissivity) per unit volume in erg/sec/cm^3
			if( Lambda_tot .ne. Lambda_tot .or. Lambda_tot-1. .eq. Lambda_tot ) then
				print '(a,2(1x,Es12.5))', 'Lambda_tot problem ', rgas(i)/Rs_init, Lambda_tot
				print '(a,5(1x,Es12.5))', 'components ', cool1, metal1, heat1, cool_com1, heat_com1
				print '(a,3(1x,Es12.5))', 'nH, T2, Zmet=', nH_gas, T2_gas, met_solar_gas
			end if

			!!! Net cooling !!!
			Lambda_net   = cool1 + met_solar_gas * metal1 - heat1 + ( cool_com1 - heat_com1 ) / nH_gas
			Lambda_net   = Lambda_net * nH_gas**2		!!! net cooling rate (emissivity) per unit volume in erg/sec/cm^3
			if( Lambda_net .ne. Lambda_net .or. Lambda_net-1. .eq. Lambda_net ) then
				print '(a,2(1x,Es12.5))', 'Lambda_net problem ', rgas(i)/Rs_init, Lambda_net
				print '(a,5(1x,Es12.5))', 'components ', cool1, metal1, heat1, cool_com1, heat_com1
				print '(a,3(1x,Es12.5))', 'nH, T2, Zmet=', nH_gas, T2_gas, met_solar_gas
			end if

			!!! Mean molecular weight !!!
			mu1_1 = mu_tab(i_nH  ,i_T2) + ( ( mu_tab(i_nH  ,i_T2+1) - mu_tab(i_nH  ,i_T2) ) / ( T2_tab(i_T2+1) - T2_tab(i_T2) ) ) * ( facT - T2_tab(i_T2) )
			mu1_2 = mu_tab(i_nH+1,i_T2) + ( ( mu_tab(i_nH+1,i_T2+1) - mu_tab(i_nH+1,i_T2) ) / ( T2_tab(i_T2+1) - T2_tab(i_T2) ) ) * ( facT - T2_tab(i_T2) )
			mu1   = mu1_1 + ( ( mu1_2 - mu1_1 ) / ( dens_tab(i_nH+1) - dens_tab(i_nH) ) ) * ( facH - dens_tab(i_nH) )
			if( mu1 .lt. 0.1_8 .or. mu1 .gt. 10.0_8 .or. mu1 .ne. mu1 .or. mu1-1 .eq. mu1 ) then
				print '(a,4(1x,Es12.5))', 'mu issue ', rgas(i)/Rs_init, mu1, mu1_1, mu1_2
				print '(a,3(1x,Es12.5))', 'nH, T2, Zmet=', nH_gas, T2_gas, met_solar_gas
				mu1 = mu_tab(i_nH, i_T2)
			end if

			!!! Gas temperature !!!
			T_gas = T2_gas*mu1				!!! gas temperature in Kelvin
			if( T_gas .ne. T_gas .or. T_gas-1. .eq. T_gas ) then
				print '(a,2(1x,Es12.5))', 'Tgas problem at ', rgas(i)/Rs_init, T_gas
				print '(a,3(1x,Es12.5))', 'nH, T2, Zmet=', nH_gas, T2_gas, met_solar_gas
			end if

			!!! Electron density !!!
			ne1 = nspec_tab(i_nH  ,i_T2, 1) + ( ( nspec_tab(i_nH  ,i_T2+1, 1) - nspec_tab(i_nH  ,i_T2, 1) ) / ( T2_tab(i_T2+1) - T2_tab(i_T2) ) ) * ( facT - T2_tab(i_T2) )
			ne2 = nspec_tab(i_nH+1,i_T2, 1) + ( ( nspec_tab(i_nH+1,i_T2+1, 1) - nspec_tab(i_nH+1,i_T2, 1) ) / ( T2_tab(i_T2+1) - T2_tab(i_T2) ) ) * ( facT - T2_tab(i_T2) )
			ne_gas   = ne1 + ( ( ne2 - ne1 ) / ( dens_tab(i_nH+1) - dens_tab(i_nH) ) ) * ( facH - dens_tab(i_nH) )
			if( ne_gas .ne. ne_gas .or. ne_gas-1 .eq. ne_gas ) then
				print '(a,4(1x,Es12.5))', 'ne issue ', rgas(i)/Rs_init, ne_gas, ne1, ne2
				print '(a,3(1x,Es12.5))', 'nH, T2, Zmet=', nH_gas, T2_gas, met_solar_gas
				ne_gas = nspec_tab(i_nH, i_T2, 1)
			end if
			ne_gas = 10.0**ne_gas

			!!! HI density !!!
			nHI1 = nspec_tab(i_nH  ,i_T2, 2) + ( ( nspec_tab(i_nH  ,i_T2+1, 2) - nspec_tab(i_nH  ,i_T2, 2) ) / ( T2_tab(i_T2+1) - T2_tab(i_T2) ) ) * ( facT - T2_tab(i_T2) )
			nHI2 = nspec_tab(i_nH+1,i_T2, 2) + ( ( nspec_tab(i_nH+1,i_T2+1, 2) - nspec_tab(i_nH+1,i_T2, 2) ) / ( T2_tab(i_T2+1) - T2_tab(i_T2) ) ) * ( facT - T2_tab(i_T2) )
			nHI_gas   = nHI1 + ( ( nHI2 - nHI1 ) / ( dens_tab(i_nH+1) - dens_tab(i_nH) ) ) * ( facH - dens_tab(i_nH) )
			if( nHI_gas .ne. nHI_gas .or. nHI_gas-1 .eq. nHI_gas ) then
				print '(a,4(1x,Es12.5))', 'nHI issue ', rgas(i)/Rs_init, nHI_gas, nHI1, nHI2
				print '(a,3(1x,Es12.5))', 'nH, T2, Zmet=', nH_gas, T2_gas, met_solar_gas
				nHI_gas = nspec_tab(i_nH, i_T2, 2)
			end if
			nHI_gas = 10.0**nHI_gas

			!!! Lyman Alpha cooling - This assumes only excitation cooling based on Goerdt et al 2010 !!!
			h_nu_LyA = 1.63d-11					!!! In [erg]. See Goerdt et al 2010, section 4
			q_LyA    = ( 2.41d-6/sqrt(T_gas) ) * ( (T_gas/1.d4)**0.22 ) * exp( -1.0 * h_nu_LyA / (KB*T_gas) )	!!! in [cm^3 s^-1]. See Goerdt et al 2010, section 4
			Lambda_LyA = ne_gas * nHI_gas * q_LyA * h_nu_LyA	!!! in [erg/s/cm^3]. See Goerdt et al 2010, section 4

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			vol_gas  = real(cell_size_gas(i)**3, 8) * unit_length_cgs**3		!!! Cell volume in cm^3
			mass_gas = real(density_gas(i), 8) * unit_density_cgs * vol_gas		!!! Cell mass in gr
			mass_gas = mass_gas / Msun_in_gr					!!! Cell mass in Msun

			j = ceiling( (log10(nH_gas) - Dmin) / d_dens_vec )
			j = max(j,1)
			j = min(j,Ndens_vec)

			k = ceiling( (log10(T_gas) - Tmin) / d_temp_vec )
			k = max(k, 1)
			k = min(k, Ntemp_vec)

			mass_grid(j,k)   = mass_grid(j,k)   + mass_gas
			emis_grid(j,k,1) = emis_grid(j,k,1) + Lambda_net * vol_gas			!!! Net (cooling-heating) emissivity from cell in erg/sec
			emis_grid(j,k,2) = emis_grid(j,k,2) + Lambda_tot * vol_gas			!!! Total cooling emissivity from cell in erg/sec
			emis_grid(j,k,3) = emis_grid(j,k,3) + Lambda_LyA * vol_gas			!!! Approximate Lyman alpha emissivity from cell in erg/sec

			if( T2_gas .lt. T2max ) then
				luminosity_vec(1,snap_num) = luminosity_vec(1,snap_num) + Lambda_net * vol_gas
				luminosity_vec(2,snap_num) = luminosity_vec(2,snap_num) + Lambda_tot * vol_gas
			end if

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!! Renormalize properties to have them of order unity for summation
			vol_gas  = real(cell_size_gas(i)**3, 8)		!!! Cell volume in code units
			mass_gas = real(density_gas(i), 8) * vol_gas	!!! Cell mass in code units
			Lambda_net = Lambda_net * 1.d27
			Lambda_tot = Lambda_tot * 1.d27
			Lambda_LyA = Lambda_LyA * 1.d27

			vol_vec(k)    = vol_vec(k)    + vol_gas
			emis_vec(k,1) = emis_vec(k,1) + Lambda_net * vol_gas
			emis_vec(k,2) = emis_vec(k,2) + Lambda_tot * vol_gas
			emis_vec(k,3) = emis_vec(k,3) + Lambda_LyA * vol_gas

			j = ceiling( real(rgas(i),8) / res )
			j = max(j,1)
			j = min(j,Nbin)
			Vbin(j) = Vbin(j) + vol_gas
			Mbin(j) = Mbin(j) + mass_gas
			Ebin(j,1) = Ebin(j,1) + Lambda_net * vol_gas
			Ebin(j,2) = Ebin(j,2) + Lambda_tot * vol_gas
			Ebin(j,3) = Ebin(j,3) + Lambda_LyA * vol_gas
			Dbin(j) = Dbin(j) + nH_gas     * vol_gas
			Tbin(j) = Tbin(j) + T_gas      * mass_gas
		end do
		print *, ''
		print *, 'Number of bad T cells ', N_bad_T

		mass_grid = mass_grid / (d_dens_vec * d_temp_vec)
		emis_grid = emis_grid / (d_dens_vec * d_temp_vec)
		write(30) Ndens_vec, Ntemp_vec, dens_vec, temp_vec, mass_grid
		write(31) Ndens_vec, Ntemp_vec, dens_vec, temp_vec, emis_grid(:,:,1), emis_grid(:,:,2), emis_grid(:,:,3)

		emis_vec = emis_vec / 1.d27		!!! Renormalize emissivity
		do k=1,Ntemp_vec
			if( vol_vec(k) .gt. 0.0_8 ) then
				emis_vec(k,:) = emis_vec(k,:) / vol_vec(k)
				write(32,'(4(1x,Es12.5))') temp_vec(k), emis_vec(k,1), emis_vec(k,2), emis_vec(k,3)
			end if
		end do

		Ebin = Ebin / 1.d27			!!! Renormalize emissivity
		do j=1,Nbin
!			print '(8(1x,Es12.5))', Rbin(j), Vbin(j), Mbin(j), Tbin(j), Dbin(j), Ebin(j,1), Ebin(j,2), Ebin(j,3)
			if( Vbin(j) .gt. 0.0_8 ) then
				Ebin(j,:) = Ebin(j,:) / Vbin(j)
				Dbin(j) = Dbin(j) / Vbin(j)
				Tbin(j) = Tbin(j) / Mbin(j)
				write(33,'(6(1x,Es12.5))') Rbin(j), Tbin(j), Dbin(j), Ebin(j,1), Ebin(j,2), Ebin(j,3)
			end if
		end do

		close(unit=30)
		close(unit=31)
		close(unit=32)
		close(unit=33)

		deallocate( Rbin, Tbin, Ebin, Dbin, Mbin, Vbin )
	end subroutine analysis
end module compute_cooling_rates
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program main
	use parameters
	use read_binary
	use compute_cooling_rates
	implicit none
	integer :: i, j, k, Nsnapshot, Nsnap
	character(len=256) :: filename, temp_label, input_arg

	if(iargc().ge.1) then
		call getarg(1,input_arg)
		read(input_arg,*) Nsnap
		if(iargc().ge.2) then
			call getarg(2,input_arg)
			read(input_arg,*) T1max		! sigma for smoothing in units of Rs
		else
			T1max = 6.d5
		end if
	else
		Nsnap = 1
		T1max = 8.d5
	end if
	T2max = T1max / 0.5885_8
	print *, 'Nsnap, T1max, T2max'
	print *, Nsnap, T1max, T2max

	!!! Read in cooling table !!!
	write(filename,'(a)') './output/cooling.out'
	open(unit=16,file=filename,form='unformatted')
	read(16) NDens, NTemp
	allocate( dens_tab(NDens), T2_tab(NTemp), mu_tab(NDens,NTemp), nspec_tab(NDens,NTemp,6) )
	allocate( cool_tab(NDens,NTemp), heat_tab(NDens,NTemp), cool_com_tab(NDens,NTemp), heat_com_tab(NDens,NTemp), metal_tab(NDens,NTemp) )
	allocate( cool_prime_tab(NDens,NTemp), heat_prime_tab(NDens,NTemp), cool_com_prime_tab(NDens,NTemp), heat_com_prime_tab(NDens,NTemp), metal_prime_tab(NDens,NTemp) )
	read(16) dens_tab
	read(16) T2_tab
	read(16) cool_tab
	read(16) heat_tab
	read(16) cool_com_tab
	read(16) heat_com_tab
	read(16) metal_tab
	read(16) cool_prime_tab
	read(16) heat_prime_tab
	read(16) cool_com_prime_tab
	read(16) heat_com_prime_tab
	read(16) metal_prime_tab
	read(16) mu_tab
	read(16) nspec_tab
	close(unit=16)
	print *, 'NDens=',NDens
	print *, 'NTemp=',NTemp
	print *, ' '
	print *, 'Dens_tab=',dens_tab
	print *, ' '
	print *, 'T2_tab=',T2_tab

	dlog_nH = 1.0_8*(NDens-1) / ( dens_tab(NDens) - dens_tab(1) )
	dlog_T2 = 1.0_8*(NTemp-1) / ( T2_tab(NTemp)   - T2_tab(1) )
	h       = 1.0/dlog_T2
	h2      = h**2
	h3      = h**3

	!!! Read in unit normalizations !!!
	write(filename,'(a,i5.5,a,i5.5,a)') './output/output_00001/info_00001.txt'
	open(unit=17,file=filename,form='formatted')
	read(17,*)
	read(17,*)
	read(17,*)
	read(17,*)
	read(17,*)
	read(17,*)
	read(17,*)
	read(17,*)
	read(17,*)
	read(17,*)
	read(17,*)
	read(17,*)
	read(17,*)
	read(17,*)
	read(17,*)
	read(17,'(A13,E23.15)')temp_label, unit_length_cgs
	read(17,'(A13,E23.15)')temp_label, unit_density_cgs
	read(17,'(A13,E23.15)')temp_label, unit_time_cgs
	close(unit=17)
	print *, 'Unit Length  in cgs:', unit_length_cgs
	print *, 'Unit Density in cgs:', unit_density_cgs
	print *, 'Unit Time    in cgs:', unit_time_cgs

	unit_velocity_cgs = unit_length_cgs / unit_time_cgs
	unit_mass_cgs     = unit_density_cgs * unit_length_cgs**3
	unit_energy_cgs   = unit_mass_cgs * unit_velocity_cgs**2

	print *, 'Unit Velocity  in cgs:', unit_velocity_cgs
	print *, 'Unit Mass      in cgs:', unit_mass_cgs
	print *, 'Unit Energy    in cgs:', unit_energy_cgs

	!!! Read in snapshot list !!!
	open(unit=18,file='./output/time.txt')
	read(18,*) Nsnapshot
	print *, Nsnapshot
	allocate( ART_file_name(Nsnapshot), tsnap(Nsnapshot) )
	do k=1,Nsnapshot
		read(18,'(F6.4)') tsnap(k)
		write(ART_file_name(k),'(a,F6.4,a)') './output/ART_format_t',tsnap(k),'.dat'
	end do
	close(unit=18)
	print *, trim(ART_file_name(1))
	print *, trim(ART_file_name(Nsnapshot))

	!!! Open output directory for cooling rates as a function of gas densty and temperature !!!
	write(filename,'(a)') 'mkdir -p ./stream_analysis/phase_diagrams'
	call system(filename)
	filename = ''

	!!! Initialize initial stream properties for sanity check !!!
	Rs_init = 0.0_8
	Vs_init = 0.0_8
	Vb_init = 0.0_8
	Cs_init = 0.0_8
	Cb_init = 0.0_8
	Rhos_init = 0.0_8
	Rhob_init = 0.0_8
	Press_init = 0.0_8
	Presb_init = 0.0_8
	Ts_init = 0.0_8
	Tb_init = 0.0_8

	!!! Set up grids for mass and emissivity weighted distributions !!!
	Dmin = -5.5
	Dmax = -0.5
	Ndens_vec = 500
	d_dens_vec = (Dmax-Dmin)/real(Ndens_vec)

	Tmin = 3.0
	Tmax = 7.0
	Ntemp_vec = 400
	d_temp_vec = (Tmax-Tmin)/real(Ntemp_vec)

	allocate( dens_vec(Ndens_vec), temp_vec(Ntemp_vec) )
	dens_vec = Dmin + (/ (i,i=1,Ndens_vec) /) * d_dens_vec - 0.5_8 * d_dens_vec
	temp_vec = Tmin + (/ (i,i=1,Ntemp_vec) /) * d_temp_vec - 0.5_8 * d_temp_vec
	allocate( mass_grid(Ndens_vec, Ntemp_vec), emis_grid(Ndens_vec, Ntemp_vec, 3) )
	allocate( vol_vec(Ntemp_vec),  emis_vec(Ntemp_vec,3) )

	print *, 'alocating gas arrays'
	allocate( cell_size_gas(Nmax), density_gas(Nmax), pressure_gas(Nmax), colour_gas(Nmax), met_gas(Nmax), rgas(Nmax) )

	allocate( luminosity_vec(2,Nsnapshot) )
	luminosity_vec(:,:) = 0.0_8

	print *, ''
	print *, 'Entering main routine'
	print *, 't=',tsnap(Nsnap)
	call analysis( Nsnap )

	print *, 'Writing Output'
	write(filename,'(a)') 'mkdir -p ./stream_analysis/Luminosity'
	call system(filename)
	filename = ''
	write(filename,'(a,I5.5,a)') './stream_analysis/Luminosity/Luminosity',Nsnap,'.txt'
	open(unit=40,file=filename,form='formatted')
	filename = ''
	write(40,'(F6.4,2(1x,ES12.5))') tsnap(Nsnap), luminosity_vec(1:2,Nsnap)
	close(unit=40)

	deallocate( luminosity_vec )
	deallocate( cell_size_gas, density_gas, pressure_gas, colour_gas, met_gas, rgas )
	deallocate( vol_vec, emis_vec )
	deallocate( mass_grid, emis_grid )
	deallocate( dens_vec, temp_vec )
	deallocate( dens_tab, T2_tab, mu_tab, nspec_tab )
	deallocate( cool_tab, heat_tab, cool_com_tab, heat_com_tab, metal_tab )
	deallocate( cool_prime_tab, heat_prime_tab, cool_com_prime_tab, heat_com_prime_tab, metal_prime_tab )

end program main
