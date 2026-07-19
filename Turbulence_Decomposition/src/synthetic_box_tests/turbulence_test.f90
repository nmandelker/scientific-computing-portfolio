module parameters
!!! Parameters used in the code, which are read in from the file parameter_input.dat
	implicit none
	real(8) :: pi

	!!! Variables for turbulence generation (from Krumholz code)
	integer :: Nseed	! Seed for random number generator
	integer :: kmin		! Minimum wavenumber, in units of 2*pi/Nx0 (smallest wavenumber increment)
	integer :: kmax		! Maximum wavenumber, in units of 2*pi/Nz0 (largest wavenumber increment)
	integer :: Nx0		! Number of grid cells in the x direction when initializing turbulence
	integer :: Ny0		! Number of grid cells in the y direction when initializing turbulence
	integer :: Nz0		! Number of grid cells in the z direction when initializing turbulence
	real(8) :: alpha	! Negative of power law slope. Power ~ k^-alpha. Supersonic turbulence is near alpha=2. Driving over a narrow band of two modes is often done with alpha=0
	real(8) :: f_solenoidal	! Determines the volume fraction of solenoidal power relative to the total (At low wave numbers, this is sensitve to the choice of radom seed). If <0 or >1, the initial Gaussian field is unchanged. 
	logical :: use_sine	! When defining compressive modes, use sin(2pik/N) rather than k (discreet formula)
	real(8) :: Vrot		! Normalization of uniform rotation velocity to be added to turbulence, in units if turbulent velocity amplitude. If Vrot<=0, then no rotation is added.
	real(8) :: beta	! Slope of rotation curve, if added

	!!! Variables for decomposition method
	logical :: pre_smooth			!!! smooth data with a Gaussian prior to decomposing
	integer :: pre_smooth_dim		!!! Dimension of Gaussian smoothing. 2D (like in Innoe et al 2015) or 3D
	integer :: pre_smooth_FWHM		!!! FWHM of Gaussian used for smoothing the fields, in number of cells ( sigma ~ 0.425*FWHM )
	integer :: decomp_type			!!! 1 for full field, 2 for rotation curve subtracted, 3 for subtracting smoother field
	integer :: type_3_FWHM			!!! FWHM of Gaussian used for defining field to be decomposed if decomp_type=3, in number of cells
	logical :: fourier_decomp		!!! perform Helmholtz decomposition in Fourier space. Otherwise use strain rate tensor in real space
	integer :: local_version		!!! 1 for diagonal of strain rate tensor, 2 for full divergence times {\vec {r}}

	!!! Variables for volume of region within which turbulence is defined and decomposed
	real(8) :: sphererad	! If >0, then perturbations are set to zero outside spherical region, and the perturbation field is shifted and renormalized to keep the center of mass velocity at zero and the variance at unity; the spherical region cut out is centered at the center of the perturbation cube, and has a radius given by the value of this parameter, with sphererad = 1 corresponding to the spherical region going all the way to the edge of the perturbation cube
	real(8) :: cylheigh	! If >0, then perturbations are set to zero outside cylindrical region, and the perturbation field is shifted and renormalized to keep the center of mass velocity at zero and the variance at unity; the cylindrical region cut out is centered at the center of the perturbation cube, and has a height given by the value of this parameter, with cylheigh = 1 corresponding to the height of the region going all the way to the edge of the perturbation cube. The radius of the cylinder is given by sphererad
	logical :: Ebox		! If extracting sphere or cylinder, rather than set cells outside geometric region within the large grid to 0, create smaller grid which won't have any 0s.

	integer :: Nxf		! Number of grid cells in the x direction to decompose, AFTER initializing the turbulence in a cubic grid
	integer :: Nyf		! Number of grid cells in the y direction to decompose, AFTER initializing the turbulence in a cubic grid
	integer :: Nzf		! Number of grid cells in the z direction to decompose, AFTER initializing the turbulence in a cubic grid

	character(len=2048) :: output_dirname
end module parameters

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module generate
use parameters
use mkl_vsl_type
use mkl_vsl
use MKL_DFTI
	implicit none
	complex(8),allocatable :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
	real(8),allocatable :: kx(:,:,:), ky(:,:,:), kz(:,:,:), kmag(:,:,:)
	real(8),allocatable :: pertx(:,:,:), perty(:,:,:), pertz(:,:,:)
	type(VSL_STREAM_STATE) :: stream_state_desc

contains
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine init_perturbations
		implicit none
		real(8),allocatable :: phasex(:,:,:), phasey(:,:,:), phasez(:,:,:)
		integer :: i, j, m, n, Nturb, error, Nxh, Nyh, Nzh
		complex(8) :: im
		real(8) :: dk, ran_temp(1)
		logical,allocatable :: kmask(:,:,:)

		print *, 'initiating perturbations'
		im = (0.0_8, 1.0_8)
		kx = 0.0_8
		ky = 0.0_8
		kz = 0.0_8
		kmag = 0.0_8
		! Perform fft k-ordering convention shifts
		! Nx0, Ny0 and Nz0 should all be even
		Nxh = Nx0/2
		Nyh = Ny0/2
		Nzh = Nz0/2
		do m=1,Nz0
			do j=1,Ny0
				kx(1:Nxh, j, m)       = real( (/ (i,i=1,Nxh) /) - 1 )
				kx((Nxh+1):Nx0, j, m) = real( (/ (i,i=1,Nxh) /) - 1 - Nxh )
			end do
		end do
		do m=1,Nz0
			do i=1,Nx0
				ky(i, 1:Nyh, m)       = real( (/ (j,j=1,Nyh) /) - 1 )
				ky(i, (Nyh+1):Ny0, m) = real( (/ (j,j=1,Nyh) /) - 1 - Nyh )
			end do
		end do
		do j=1,Ny0
			do i=1,Nx0
				kz(i, j, 1:Nzh)       = real( (/ (m,m=1,Nzh) /) - 1 )
				kz(i, j, (Nzh+1):Nz0) = real( (/ (m,m=1,Nzh) /) - 1 - Nzh )
			end do
		end do
		allocate( kmask(Nx0,Ny0,Nz0), phasex(Nx0,Ny0,Nz0), phasey(Nx0,Ny0,Nz0), phasez(Nx0,Ny0,Nz0) )
		kx = kx * 2.0_8*pi/real(Nx0)
		ky = ky * 2.0_8*pi/real(Ny0)
		kz = kz * 2.0_8*pi/real(Nz0)
		kmag = sqrt( kx**2 + ky**2 + kz**2 )

		dk = 2.0_8*pi/real(Nx0)		! Nx0 is necessarily the largest, so 2pi/Nx0 is the smallest dk interval available
		kmask = kmag .ge. real(kmin)*dk .and. kmag .lt. real(kmax+1)*dk
		Nturb = count( kmask )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		phasex = 0.0_8
		phasey = 0.0_8
		phasez = 0.0_8
		fx = 0.0_8
		fy = 0.0_8
		fz = 0.0_8
		print *, 'generating random numbers'
		do m=1,Nz0
			do j=1,Ny0
				do i=1,Nx0
					if( kmask(i,j,m) ) then
						error = vdrnguniform( VSL_RNG_METHOD_UNIFORM_STD, stream_state_desc, 1, phasex(i,j,m), 0.0_8, 1.0_8 )
						error = vdrnguniform( VSL_RNG_METHOD_UNIFORM_STD, stream_state_desc, 1, phasey(i,j,m), 0.0_8, 1.0_8 )
						error = vdrnguniform( VSL_RNG_METHOD_UNIFORM_STD, stream_state_desc, 1, phasez(i,j,m), 0.0_8, 1.0_8 )
						error = vdrnggaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream_state_desc, 1, ran_temp, 0.0_8, 1.0_8 )
						fx(i,j,m) = ran_temp(1)
						error = vdrnggaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream_state_desc, 1, ran_temp, 0.0_8, 1.0_8 )
						fy(i,j,m) = ran_temp(1)
						error = vdrnggaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream_state_desc, 1, ran_temp, 0.0_8, 1.0_8 )
						fz(i,j,m) = ran_temp(1)
					end if
				end do
			end do
		end do
		print *, 'done'
		phasex = 2.0_8 * pi * phasex
		phasey = 2.0_8 * pi * phasey
		phasez = 2.0_8 * pi * phasez
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Take away two powers to balance out the fact that the volume element dV \propto k^2 dk, meaning there are more gridpoints at higher k.
		where( kmask ) fx = fx * (kmag**(-0.5*(alpha+2)))
		where( kmask ) fy = fy * (kmag**(-0.5*(alpha+2)))
		where( kmask ) fz = fz * (kmag**(-0.5*(alpha+2)))

		! Add in phases
		where( kmask ) fx = cos(phasex)*fx + im*sin(phasex)*fx
		where( kmask ) fx = cos(phasey)*fy + im*sin(phasey)*fy
		where( kmask ) fx = cos(phasez)*fz + im*sin(phasez)*fz

		deallocate( kmask, phasex, phasey, phasez )
	end subroutine init_perturbations
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine make_perturbations
		implicit none
		complex(8),allocatable :: fxs(:,:,:), fys(:,:,:), fzs(:,:,:), fxc(:,:,:), fyc(:,:,:), fzc(:,:,:)
		complex(8),allocatable :: fx1d(:), fy1d(:), fz1d(:)
		complex(8),allocatable :: fxs1d(:), fys1d(:), fzs1d(:), fxc1d(:), fyc1d(:), fzc1d(:), PS_out(:,:)
		real(8) :: norms, normc, rescale
		real(8),allocatable :: pertsx(:,:,:), pertsy(:,:,:), pertsz(:,:,:), pertcx(:,:,:), pertcy(:,:,:), pertcz(:,:,:)
		real(8),allocatable :: test(:,:,:), Skmag(:,:,:), kvec(:)
		real(8) :: rad, rx, ry, dk1, dk2, dk3
		real(8),allocatable :: Vx(:,:), Vy(:,:), Vz(:,:), Vr(:,:), Vphi(:,:)
		logical,allocatable :: kmask(:,:,:)
		integer :: i, j, m, n1, n2, n3, nvec, n4, length(3), error
		character(len=256) :: filename
		type(DFTI_DESCRIPTOR), POINTER :: fft_desc

		allocate( kx(Nx0,Ny0,Nz0), ky(Nx0,Ny0,Nz0), kz(Nx0,Ny0,Nz0), kmag(Nx0,Ny0,Nz0) )
		allocate( fx(Nx0,Ny0,Nz0), fy(Nx0,Ny0,Nz0), fz(Nx0,Ny0,Nz0) )
		error = vslnewstream( stream_state_desc, VSL_BRNG_MT19937, Nseed )
		call init_perturbations
		if( use_sine) then
			allocate( Skmag(Nx0,Ny0,Nz0) )
			Skmag = kmag
			kx = sin(kx)
			ky = sin(ky)
			kz = sin(kz)
			kmag = sqrt( kx**2 + ky**2 + kz**2 )
		end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		print *, 'starting decomposition of initial field'
		allocate( fxs(Nx0,Ny0,Nz0), fys(Nx0,Ny0,Nz0), fzs(Nx0,Ny0,Nz0), fxc(Nx0,Ny0,Nz0), fyc(Nx0,Ny0,Nz0), fzc(Nx0,Ny0,Nz0) )
		fxs = 0.0_8
		fys = 0.0_8
		fzs = 0.0_8
		fxc = 0.0_8
		fyc = 0.0_8
		fzc = 0.0_8

		! compressive part
		fxc = kx*( kx*fx + ky*fy + kz*fz )/max(kmag**2, 1.d-16)
		fyc = ky*( kx*fx + ky*fy + kz*fz )/max(kmag**2, 1.d-16)
		fzc = kz*( kx*fx + ky*fy + kz*fz )/max(kmag**2, 1.d-16)
		where( kmag .eq. 0.0_8 ) fxc = 0.0_8
		where( kmag .eq. 0.0_8 ) fyc = 0.0_8
		where( kmag .eq. 0.0_8 ) fzc = 0.0_8
		normc = sqrt( sum(abs(fxc)**2 + abs(fyc)**2 + abs(fzc)**2) )
		deallocate( kx, ky, kz )

		! solenoidal part
		fxs = fx - fxc
		fys = fy - fyc
		fzs = fz - fzc
		where( kmag .eq. 0.0_8 ) fxs = 0.0_8
		where( kmag .eq. 0.0_8 ) fys = 0.0_8
		where( kmag .eq. 0.0_8 ) fzs = 0.0_8
		norms = sqrt( sum(abs(fxs)**2 + abs(fys)**2 + abs(fzs)**2) )

		! close random number stream and deallocate
		error = vsldeletestream( stream_state_desc )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		print *, '!!! Output Initial Power Spectra !!!'
		if( use_sine ) then 
			kmag = Skmag
			deallocate( Skmag )
		end if

		dk1 = 2.0_8*pi/Nx0
		dk2 = 2.0_8*pi/Ny0
		dk3 = 2.0_8*pi/Nz0
		n1 = ceiling( dk2/dk1 )
		n2 = ceiling( dk3/dk2 )
		n3 = ceiling( sqrt(3.0_8)*pi / dk3 ) + 1
		nvec = n1 + n2 + n3 - 2
		print *, 'n1, n2, n3, nvec ',n1, n2, n3, nvec
		allocate( kvec(nvec), PS_out(nvec,13) )
		kvec(1:n1)  = (/ (i,i=1,n1) /) * dk1
		kvec(n1:(n1+n2-1)) = (/ (i,i=1,n2) /) * dk2
		kvec((n1+n2-1):nvec) = (/ (i,i=1,n3) /) * dk3
		PS_out = 0.0_8
		do m=1,Nz0
			do j=1,Ny0
				do i=1,Nx0
					if( kmag(i,j,m) .ge. dk3 ) then
						n4 = min( floor( kmag(i,j,m)/dk3 ) + n1 + n2 - 1, nvec )
					elseif( kmag(i,j,m) .ge. dk2 ) then
						n4 = floor( kmag(i,j,m)/dk2 ) + n1
					else
						n4 = floor( kmag(i,j,m)/dk1 ) + 1
					end if
					PS_out(n4,1)  = PS_out(n4,1)  + abs(fx(i,j,m))**2  + abs(fy(i,j,m))**2  + abs(fz(i,j,m))**2
					PS_out(n4,2)  = PS_out(n4,2)  + abs(fxs(i,j,m))**2 + abs(fys(i,j,m))**2 + abs(fzs(i,j,m))**2
					PS_out(n4,3)  = PS_out(n4,3)  + abs(fxc(i,j,m))**2 + abs(fyc(i,j,m))**2 + abs(fzc(i,j,m))**2
					PS_out(n4,4)  = PS_out(n4,4)  + fx(i,j,m)
					PS_out(n4,5)  = PS_out(n4,5)  + fy(i,j,m)
					PS_out(n4,6)  = PS_out(n4,6)  + fz(i,j,m)
					PS_out(n4,7)  = PS_out(n4,7)  + fxs(i,j,m)
					PS_out(n4,8)  = PS_out(n4,8)  + fys(i,j,m)
					PS_out(n4,9)  = PS_out(n4,9)  + fzs(i,j,m)
					PS_out(n4,10) = PS_out(n4,10) + fxc(i,j,m)
					PS_out(n4,11) = PS_out(n4,11) + fyc(i,j,m)
					PS_out(n4,12) = PS_out(n4,12) + fzc(i,j,m)
					PS_out(n4,13) = PS_out(n4,13) + 1.0_8
				end do
			end do
		end do
		deallocate( kmag )
		write(filename,'(a,a)') trim(output_dirname),'/initial_power_spectra.out'
		open(unit=20,file=filename)
		do i=1,n1
			write(20,'(7(1x,es12.5))') kvec(i)-0.5_8*dk1, real(PS_out(i,1)), real(PS_out(i,2)), real(PS_out(i,3)), & 
				& abs(PS_out(i,4)/PS_out(i,13))**2  + abs(PS_out(i,5)/PS_out(i,13))**2  + abs(PS_out(i,6)/PS_out(i,13))**2, &
				& abs(PS_out(i,7)/PS_out(i,13))**2  + abs(PS_out(i,8)/PS_out(i,13))**2  + abs(PS_out(i,9)/PS_out(i,13))**2, & 
				& abs(PS_out(i,10)/PS_out(i,13))**2 + abs(PS_out(i,11)/PS_out(i,13))**2 + abs(PS_out(i,12)/PS_out(i,13))**2
		end do
		if( n1+n2-1 .gt. 1 ) then
			do i=n1+1,n1+n2-1
				write(20,'(7(1x,es12.5))') kvec(i)-0.5_8*dk2, real(PS_out(i,1)), real(PS_out(i,2)), real(PS_out(i,3)), & 
					& abs(PS_out(i,4)/PS_out(i,13))**2  + abs(PS_out(i,5)/PS_out(i,13))**2  + abs(PS_out(i,6)/PS_out(i,13))**2, &
					& abs(PS_out(i,7)/PS_out(i,13))**2  + abs(PS_out(i,8)/PS_out(i,13))**2  + abs(PS_out(i,9)/PS_out(i,13))**2, & 
					& abs(PS_out(i,10)/PS_out(i,13))**2 + abs(PS_out(i,11)/PS_out(i,13))**2 + abs(PS_out(i,12)/PS_out(i,13))**2
			end do
		end if
		do i=n1+n2,nvec
			write(20,'(7(1x,es12.5))') kvec(i)-0.5_8*dk3, real(PS_out(i,1)), real(PS_out(i,2)), real(PS_out(i,3)), & 
				& abs(PS_out(i,4)/PS_out(i,13))**2  + abs(PS_out(i,5)/PS_out(i,13))**2  + abs(PS_out(i,6)/PS_out(i,13))**2, &
				& abs(PS_out(i,7)/PS_out(i,13))**2  + abs(PS_out(i,8)/PS_out(i,13))**2  + abs(PS_out(i,9)/PS_out(i,13))**2, & 
				& abs(PS_out(i,10)/PS_out(i,13))**2 + abs(PS_out(i,11)/PS_out(i,13))**2 + abs(PS_out(i,12)/PS_out(i,13))**2
		end do
		close(unit=20)
		print *, 'done'
		deallocate( kvec, PS_out )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! total fourier vector
		print *, 'total fourier vector'
		if( f_solenoidal .ge. 0.0_8 .and. f_solenoidal .le. 1.0_8 ) then
			fx = f_solenoidal*fxs + (1.0_8-f_solenoidal)*fxc
			fy = f_solenoidal*fys + (1.0_8-f_solenoidal)*fyc
			fz = f_solenoidal*fzs + (1.0_8-f_solenoidal)*fzc
		end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! back to real space
		print *, 'back to real space'
		n1 = Nx0*Ny0*Nz0
		allocate( fx1d(n1), fy1d(n1), fz1d(n1), fxs1d(n1), fys1d(n1), fzs1d(n1), fxc1d(n1), fyc1d(n1), fzc1d(n1) )
		do m=1,Nz0
			do j=1,Ny0
				do i=1,Nx0
					n1 = (m-1)*Nx0*Ny0 + (j-1)*Nx0 + i
					fx1d(n1) = fx(i,j,m)
					fy1d(n1) = fy(i,j,m)
					fz1d(n1) = fz(i,j,m)
					fxs1d(n1) = fxs(i,j,m)
					fys1d(n1) = fys(i,j,m)
					fzs1d(n1) = fzs(i,j,m)
					fxc1d(n1) = fxc(i,j,m)
					fyc1d(n1) = fyc(i,j,m)
					fzc1d(n1) = fzc(i,j,m)
				end do
			end do
		end do
		deallocate( fx, fy, fz, fxs, fys, fzs, fxc, fyc, fzc )

		print *, 'inverse fft'
		length = (/ Nx0, Ny0, Nz0 /)

		print *, 'fx'
		error = DftiCreateDescriptor(fft_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc)
		error = DftiComputeBackward(fft_desc, fx1d)
		error = DftiFreeDescriptor(fft_desc)
		print '(4(1x,es12.5))', maxval(real(fx1d)), minval(real(fx1d)), maxval(aimag(fx1d)), minval(aimag(fx1d))

		print *, 'fy'
		error = DftiCreateDescriptor(fft_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc)
		error = DftiComputeBackward(fft_desc, fy1d)
		error = DftiFreeDescriptor(fft_desc)
		print '(4(1x,es12.5))', maxval(real(fy1d)), minval(real(fy1d)), maxval(aimag(fy1d)), minval(aimag(fy1d))

		print *, 'fz'
		error = DftiCreateDescriptor(fft_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc)
		error = DftiComputeBackward(fft_desc, fz1d)
		error = DftiFreeDescriptor(fft_desc)
		print '(4(1x,es12.5))', maxval(real(fz1d)), minval(real(fz1d)), maxval(aimag(fz1d)), minval(aimag(fz1d))

		print *, 'fxs'
		error = DftiCreateDescriptor(fft_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc)
		error = DftiComputeBackward(fft_desc, fxs1d)
		error = DftiFreeDescriptor(fft_desc)
		print '(4(1x,es12.5))', maxval(real(fxs1d)), minval(real(fxs1d)), maxval(aimag(fxs1d)), minval(aimag(fxs1d))

		print *, 'fys'
		error = DftiCreateDescriptor(fft_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc)
		error = DftiComputeBackward(fft_desc, fys1d)
		error = DftiFreeDescriptor(fft_desc)
		print '(4(1x,es12.5))', maxval(real(fys1d)), minval(real(fys1d)), maxval(aimag(fys1d)), minval(aimag(fys1d))

		print *, 'fzs'
		error = DftiCreateDescriptor(fft_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc)
		error = DftiComputeBackward(fft_desc, fzs1d)
		error = DftiFreeDescriptor(fft_desc)
		print '(4(1x,es12.5))', maxval(real(fzs1d)), minval(real(fzs1d)), maxval(aimag(fzs1d)), minval(aimag(fzs1d))

		print *, 'fxc'
		error = DftiCreateDescriptor(fft_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc)
		error = DftiComputeBackward(fft_desc, fxc1d)
		error = DftiFreeDescriptor(fft_desc)
		print '(4(1x,es12.5))', maxval(real(fxc1d)), minval(real(fxc1d)), maxval(aimag(fxc1d)), minval(aimag(fxc1d))

		print *, 'fyc'
		error = DftiCreateDescriptor(fft_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc)
		error = DftiComputeBackward(fft_desc, fyc1d)
		error = DftiFreeDescriptor(fft_desc)
		print '(4(1x,es12.5))', maxval(real(fyc1d)), minval(real(fyc1d)), maxval(aimag(fyc1d)), minval(aimag(fyc1d))

		print *, 'fzc'
		error = DftiCreateDescriptor(fft_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc)
		error = DftiComputeBackward(fft_desc, fzc1d)
		error = DftiFreeDescriptor(fft_desc)
		print '(4(1x,es12.5))', maxval(real(fzc1d)), minval(real(fzc1d)), maxval(aimag(fzc1d)), minval(aimag(fzc1d))

		allocate( pertx(Nx0,Ny0,Nz0), perty(Nx0,Ny0,Nz0), pertz(Nx0,Ny0,Nz0), pertsx(Nx0,Ny0,Nz0), pertsy(Nx0,Ny0,Nz0), pertsz(Nx0,Ny0,Nz0), & 
			& pertcx(Nx0,Ny0,Nz0), pertcy(Nx0,Ny0,Nz0), pertcz(Nx0,Ny0,Nz0) )
		do m=1,Nz0
			do j=1,Ny0
				do i=1,Nx0
					n1 = (m-1)*Nx0*Ny0 + (j-1)*Nx0 + i
					pertx(i,j,m)  = real(fx1d(n1))
					perty(i,j,m)  = real(fy1d(n1))
					pertz(i,j,m)  = real(fz1d(n1))
					pertsx(i,j,m) = real(fxs1d(n1))
					pertsy(i,j,m) = real(fys1d(n1))
					pertsz(i,j,m) = real(fzs1d(n1))
					pertcx(i,j,m) = real(fxc1d(n1))
					pertcy(i,j,m) = real(fyc1d(n1))
					pertcz(i,j,m) = real(fzc1d(n1))
				end do
			end do
		end do
		deallocate( fx1d, fy1d, fz1d, fxs1d, fys1d, fzs1d, fxc1d, fyc1d, fzc1d )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		print *, '!!!!!'
		print *, 'TESTS'
		allocate( test(Nx0-2, Ny0-2, Nz0-2) )
		m = (Nx0-2)*(Ny0-2)*(Nz0-2)
		print *, 'divergence of solonoidal mode'
		test = 0.5_8 * (  pertsx(3:Nx0,     2:(Ny0-1), 2:(Nz0-1)) - pertsx(1:(Nx0-2), 2:(Ny0-1), 2:(Nz0-1)) + &
				& pertsy(2:(Nx0-1), 3:Ny0,     2:(Nz0-1)) - pertsy(2:(Nx0-1), 1:(Ny0-2), 2:(Nz0-1)) + &
				& pertsz(2:(Nx0-1), 2:(Ny0-1), 3:Nz0)     - pertsz(2:(Nx0-1), 2:(Ny0-1), 1:(Nz0-2)) )
		print '(3(1x,es12.5))', maxval(abs(test)), minval(abs(test)), sum(abs(test))/m
		print *, 'divergence of compressive mode'
		test = 0.5_8 * (  pertcx(3:Nx0,     2:(Ny0-1), 2:(Nz0-1)) - pertcx(1:(Nx0-2), 2:(Ny0-1), 2:(Nz0-1)) + &
				& pertcy(2:(Nx0-1), 3:Ny0,     2:(Nz0-1)) - pertcy(2:(Nx0-1), 1:(Ny0-2), 2:(Nz0-1)) + &
				& pertcz(2:(Nx0-1), 2:(Ny0-1), 3:Nz0)     - pertcz(2:(Nx0-1), 2:(Ny0-1), 1:(Nz0-2)) )
		print '(3(1x,es12.5))', maxval(abs(test)), minval(abs(test)), sum(abs(test))/m
		print *, 'rotor of solonoidal mode'
		test = 0.5_8 * sqrt( &
				&  ( (pertsz(2:(Nx0-1), 3:Ny0, 2:(Nz0-1)) - pertsz(2:(Nx0-1), 1:(Ny0-2), 2:(Nz0-1))) - & 
				&    (pertsy(2:(Nx0-1), 2:(Ny0-1), 3:Nz0) - pertsy(2:(Nx0-1), 2:(Ny0-1), 1:(Nz0-2))) )**2 + &
				&  ( (pertsx(2:(Nx0-1), 2:(Ny0-1), 3:Nz0) - pertsx(2:(Nx0-1), 2:(Ny0-1), 1:(Nz0-2))) - & 
				&    (pertsz(3:Nx0, 2:(Ny0-1), 2:(Nz0-1)) - pertsz(1:(Nx0-2), 2:(Ny0-1), 2:(Nz0-1))) )**2 + &
				&  ( (pertsy(3:Nx0, 2:(Ny0-1), 2:(Nz0-1)) - pertsy(1:(Nx0-2), 2:(Ny0-1), 2:(Nz0-1))) - & 
				&    (pertsx(2:(Nx0-1), 3:Ny0, 2:(Nz0-1)) - pertsx(2:(Nx0-1), 1:(Ny0-2), 2:(Nz0-1))) )**2   &
				&  )
		print '(3(1x,es12.5))', maxval(abs(test)), minval(abs(test)), sum(abs(test))/m
		print *, 'rotor of compressive mode'
		test = 0.5_8 * sqrt( &
				&  ( (pertcz(2:(Nx0-1), 3:Ny0, 2:(Nz0-1)) - pertcz(2:(Nx0-1), 1:(Ny0-2), 2:(Nz0-1))) - & 
				&    (pertcy(2:(Nx0-1), 2:(Ny0-1), 3:Nz0) - pertcy(2:(Nx0-1), 2:(Ny0-1), 1:(Nz0-2))) )**2 + &
				&  ( (pertcx(2:(Nx0-1), 2:(Ny0-1), 3:Nz0) - pertcx(2:(Nx0-1), 2:(Ny0-1), 1:(Nz0-2))) - & 
				&    (pertcz(3:Nx0, 2:(Ny0-1), 2:(Nz0-1)) - pertcz(1:(Nx0-2), 2:(Ny0-1), 2:(Nz0-1))) )**2 + &
				&  ( (pertcy(3:Nx0, 2:(Ny0-1), 2:(Nz0-1)) - pertcy(1:(Nx0-2), 2:(Ny0-1), 2:(Nz0-1))) - & 
				&    (pertcx(2:(Nx0-1), 3:Ny0, 2:(Nz0-1)) - pertcx(2:(Nx0-1), 1:(Ny0-2), 2:(Nz0-1))) )**2   &
				&  )
		print '(3(1x,es12.5))', maxval(abs(test)), minval(abs(test)), sum(abs(test))/m
		deallocate( test )

		allocate( test(Nx0, Ny0, Nz0) )
		m = Nx0*Ny0*Nz0
		if( f_solenoidal .gt. 1.0_8 .or. f_solenoidal .lt. 0.0_8 ) then
			print *, 'difference magnitude from full field'
			test = sqrt( (pertx-pertsx-pertcx)**2 + (perty-pertsy-pertcy)**2 + (pertz-pertsz-pertcz)**2 )
			print '(3(1x,es12.5))', maxval(abs(test)), minval(abs(test)), sum(abs(test))/m
		end if

		print *, 'dot product'
		test = 2.0_8*( pertcx*pertsx + pertcy*pertsy + pertcz*pertsz )
		print '(4(1x,es12.5))', maxval(abs(test)), minval(abs(test)), sum(abs(test))/m, sum(test)/m

		print *, 'total power fourier space: t, s/t c/t'
		print '(3(1x,es12.5))', norms**2+normc**2, norms**2 / (norms**2+normc**2), normc**2 / (norms**2+normc**2)
		print *, 'total power real space: t, s/t c/t, dot/t'
		norms = sum(pertsx**2 + pertsy**2 + pertsz**2) / m
		normc = sum(pertcx**2 + pertcy**2 + pertcz**2) / m
		print '(4(1x,es12.5))', norms+normc, norms/(norms+normc), normc/(norms+normc), sum(test)/(norms+normc)
		print *, 'average power ratios real space: s/t c/t, dot/t'
		test = pertsx**2 + pertsy**2 + pertsz**2 + pertcx**2 + pertcy**2 + pertcz**2
		print '(3(1x,es12.5))',   sum( (pertsx**2 + pertsy**2 + pertsz**2)/test )/m, sum( (pertcx**2 + pertcy**2 + pertcz**2)/test )/m, & 
					& sum( 2.0_8*(pertcx*pertsx + pertcy*pertsy + pertcz*pertsz)/test )/m

		deallocate( test) 
		deallocate( pertsx, pertsy, pertsz, pertcx, pertcy, pertcz )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Subtract off COM (assuming uniform density)
		print '(a3,3(1x,es12.5))', 'COM', sum(pertx)/m, sum(perty)/m, sum(pertz)/m
		pertx = pertx - ( sum(pertx)/m )
		perty = perty - ( sum(perty)/m )
		pertz = pertz - ( sum(pertz)/m )

		! Scale RMS of perturbation cube to unity
		rescale = sqrt( sum(pertx**2 + perty**2 + pertz**2)/m )
		print '(a3,1x,es12.5)', 'RMS', rescale
		pertx  = pertx / rescale
		perty  = perty / rescale
		pertz  = pertz / rescale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		print *, 'Making surface maps of initial turbulence components in cylindrical coordinates'
		allocate( Vx(Nx0,Ny0), Vy(Nx0,Ny0), Vz(Nx0,Ny0), Vr(Nx0,Ny0), Vphi(Nx0,Ny0) )
		Vx = sum(pertx, dim=3) / real(Nz0)
		Vy = sum(perty, dim=3) / real(Nz0)
		Vz = sum(pertz, dim=3) / real(Nz0)
		rx = real(Nx0+1)/2.0_8
		ry = real(Ny0+1)/2.0_8
		do j=1,Ny0
			do i=1,Nx0
				rad = sqrt( real(i-rx)**2 + real(j-ry)**2 )
				Vr(i,j)   = Vx(i,j)*real(i-rx)/rad + Vy(i,j)*real(j-ry)/rad
				Vphi(i,j) = Vx(i,j)*real(ry-j)/rad + Vy(i,j)*real(i-rx)/rad
			end do
		end do
		Vx = sqrt( sum(pertx**2 + perty**2 + pertz**2, dim=3) / real(Nz0) )
		write(filename,'(a,a)') trim(output_dirname),'/initial_turbulence.bin'
		open(unit=20,file=filename,form='unformatted')
		write(20) Nx0, Ny0, Vr, Vphi, Vz, Vx
		close(unit=20)
		deallocate( Vx, Vy, Vz, Vr, Vphi )

		if( sphererad .gt. 0.0_8 ) then
			if( cylheigh .gt. 0.0_8 ) then
				call extract_cyl()
			else
				call extract_sphere()
			end if
		else
			Nxf = Nx0
			Nyf = Ny0
			Nzf = Nz0
		end if

		if( Vrot .gt. 0.0_8 ) then
			call add_rotation()
		end if
	end subroutine make_perturbations
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine extract_sphere
		implicit none
		integer :: i, j, m, Navgx, Navgy, Navgz
		real(8) :: rx, ry, rz, avgx, avgy, avgz
		real(8),allocatable :: rad(:,:,:), tempx(:,:,:), tempy(:,:,:), tempz(:,:,:)
		logical,allocatable :: rmask(:,:,:)

		if( Ebox ) then
			Nxf = floor( sphererad*real(Nz0) )
			if( MOD(Nxf,2) .ne. 0 ) then
				Nxf = Nxf - 1
			end if
			Nyf = Nxf
			Nzf = Nxf
			allocate( tempx(Nxf, Nyf, Nzf), tempy(Nxf, Nyf, Nzf), tempz(Nxf, Nyf, Nzf) )

			Navgx = Nx0/2
			Navgy = Ny0/2
			Navgz = Nz0/2
			i = Nxf/2
			j = Nyf/2
			m = Nzf/2
			tempx = pertx(Navgx-i+1:Navgx+i, Navgy-j+1:Navgy+j, Navgz-m+1:Navgz+m)
			tempy = perty(Navgx-i+1:Navgx+i, Navgy-j+1:Navgy+j, Navgz-m+1:Navgz+m)
			tempz = pertz(Navgx-i+1:Navgx+i, Navgy-j+1:Navgy+j, Navgz-m+1:Navgz+m)
			deallocate( pertx, perty, pertz )
			allocate( pertx(Nxf, Nyf, Nzf), perty(Nxf, Nyf, Nzf), pertz(Nxf, Nyf, Nzf) )
			pertx = tempx
			perty = tempy
			pertz = tempz
			deallocate( tempx, tempy, tempz )

			Navgx = Nxf*Nyf*Nzf
			avgx = sum(pertx) / Navgx
			avgy = sum(perty) / Navgx
			avgz = sum(pertz) / Navgx
			pertx = pertx - avgx
			perty = perty - avgy
			pertz = pertz - avgz

			avgx = sqrt( sum(pertx**2 + perty**2 + pertz**2) / Navgx )
			pertx = pertx / avgx
			perty = perty / avgx
			pertz = pertz / avgx
		else
			rx = real(Nx0+1)/2.0_8
			ry = real(Ny0+1)/2.0_8
			rz = real(Nz0+1)/2.0_8
			allocate( rad(Nx0,Ny0,Nz0), rmask(Nx0,Ny0,Nz0) )
			do m=1,Nz0
				do j=1,Ny0
					do i=1,Nx0
						rad(i,j,m) = sqrt( (real(i)-rx)**2 + (real(j)-ry)**2 + (real(m)-rz)**2 )
					end do
				end do
			end do
			rz = real(Nz0)/2.0_8
			rmask = rad .le. sphererad*rz
			deallocate( rad )

			where( .not. rmask ) pertx = 0.0_8
			where( .not. rmask ) perty = 0.0_8
			where( .not. rmask ) pertz = 0.0_8

			Navgx = count( rmask )
			avgx = sum(pertx, rmask) / Navgx
			avgy = sum(perty, rmask) / Navgx
			avgz = sum(pertz, rmask) / Navgx
			where( rmask ) pertx = pertx - avgx
			where( rmask ) perty = perty - avgy
			where( rmask ) pertz = pertz - avgz

			avgx = sqrt( sum(pertx**2 + perty**2 + pertz**2, rmask) / Navgx )
			where( rmask ) pertx = pertx / avgx
			where( rmask ) perty = perty / avgx
			where( rmask ) pertz = pertz / avgx

			deallocate( rmask )
			Nxf = Nx0
			Nyf = Ny0
			Nzf = Nz0
		end if
		print *, 'Nx,Ny,Nz',Nxf,Nyf,Nzf

	end subroutine extract_sphere
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine extract_cyl
		implicit none
		integer :: i, j, m, Navgx, Navgy, Navgz
		real(8) :: rx, ry, rz, avgx, avgy, avgz
		real(8),allocatable :: rad(:,:,:), height(:,:,:), tempx(:,:,:), tempy(:,:,:), tempz(:,:,:)
		logical,allocatable :: rmask(:,:,:)

		if( Ebox ) then
			Nxf = floor( sphererad*real(Nx0) )
			if( MOD(Nxf,2) .ne. 0 ) then
				Nxf = Nxf - 1
			end if
			Nyf = Nxf
			Nzf = floor( cylheigh*real(Nz0) )
			if( MOD(Nzf,2) .ne. 0 ) then
				Nzf = Nzf - 1
			end if
			allocate( tempx(Nxf, Nyf, Nzf), tempy(Nxf, Nyf, Nzf), tempz(Nxf, Nyf, Nzf) )

			Navgx = Nx0/2
			Navgy = Ny0/2
			Navgz = Nz0/2
			i = Nxf/2
			j = Nyf/2
			m = Nzf/2
			tempx = pertx(Navgx-i+1:Navgx+i, Navgy-j+1:Navgy+j, Navgz-m+1:Navgz+m)
			tempy = perty(Navgx-i+1:Navgx+i, Navgy-j+1:Navgy+j, Navgz-m+1:Navgz+m)
			tempz = pertz(Navgx-i+1:Navgx+i, Navgy-j+1:Navgy+j, Navgz-m+1:Navgz+m)
			deallocate( pertx, perty, pertz )
			allocate( pertx(Nxf, Nyf, Nzf), perty(Nxf, Nyf, Nzf), pertz(Nxf, Nyf, Nzf) )
			pertx = tempx
			perty = tempy
			pertz = tempz
			deallocate( tempx, tempy, tempz )

			Navgx = Nxf*Nyf*Nzf
			avgx = sum(pertx) / Navgx
			avgy = sum(perty) / Navgx
			avgz = sum(pertz) / Navgx
			pertx = pertx - avgx
			perty = perty - avgy
			pertz = pertz - avgz

			avgx = sqrt( sum(pertx**2 + perty**2 + pertz**2) / Navgx )
			pertx = pertx / avgx
			perty = perty / avgx
			pertz = pertz / avgx
		else
			rx = real(Nx0+1)/2.0_8
			ry = real(Ny0+1)/2.0_8
			rz = real(Nz0+1)/2.0_8
			allocate( rad(Nx0,Ny0,Nz0), height(Nx0,Ny0,Nz0), rmask(Nx0,Ny0,Nz0) )
			do m=1,Nz0
				do j=1,Ny0
					do i=1,Nx0
						rad(i,j,m) = sqrt( (real(i)-rx)**2 + (real(j)-ry)**2 )
						height(i,j,m) = abs( real(m)-rz )
					end do
				end do
			end do
			ry = real(Ny0)/2.0_8
			rz = real(Nz0)/2.0_8
			rmask = rad .le. sphererad*ry .and. height .le. cylheigh*rz
			deallocate( rad, height )

			where( .not. rmask ) pertx = 0.0_8
			where( .not. rmask ) perty = 0.0_8
			where( .not. rmask ) pertz = 0.0_8

			Navgx = count( rmask )
			avgx = sum(pertx, rmask) / Navgx
			avgy = sum(perty, rmask) / Navgx
			avgz = sum(pertz, rmask) / Navgx
			where( rmask ) pertx = pertx - avgx
			where( rmask ) perty = perty - avgy
			where( rmask ) pertz = pertz - avgz

			avgx = sqrt( sum(pertx**2 + perty**2 + pertz**2, rmask) / Navgx )
			where( rmask ) pertx = pertx / avgx
			where( rmask ) perty = perty / avgx
			where( rmask ) pertz = pertz / avgx

			deallocate( rmask )
			Nxf = Nx0
			Nyf = Ny0
			Nzf = Nz0
		end if
		print *, 'Nx,Ny,Nz',Nxf,Nyf,Nzf

	end subroutine extract_cyl
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine add_rotation
		implicit none
		integer :: i, j, m
		real(8) :: rad, rad2, rx, ry, rz

		rx = real(Nxf+1)/2.0_8
		ry = real(Nyf+1)/2.0_8
		rz = real(Nzf+1)/2.0_8
		do j=1,Nyf
			do i=1,Nxf
				rad = sqrt( (real(i)-rx)**2 + (real(j)-ry)**2 )
				do m=1,Nzf
					rad2 = abs(real(m)-rz) / (0.5_8 * real(Nzf))
					pertx(i,j,m) = pertx(i,j,m) + Vrot * (rad**(beta-1)) * real(ry-j) / ( cosh(rad2) )**2
					perty(i,j,m) = perty(i,j,m) + Vrot * (rad**(beta-1)) * real(i-rx) / ( cosh(rad2) )**2
				end do
			end do
		end do
	end subroutine add_rotation
end module generate

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module decomp
use parameters
use generate
	implicit none
!!! Grid variables
	real(8),allocatable :: sigcx(:,:,:), sigcy(:,:,:), sigcz(:,:,:)
	real(8),allocatable :: sigsx(:,:,:), sigsy(:,:,:), sigsz(:,:,:)
	real(8) :: sigma_ratio(8)

contains
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine turb_def()
		implicit none
		integer :: i, j, k, m, n, n1, Nrot, Nvert
		real(8) :: rad, rx, ry, vrot_interp, rm, vxt, vyt
		real(8),allocatable :: Vx(:,:), Vy(:,:), Vz(:,:), Vr(:,:), Vphi(:,:), vrot_curve(:,:), rrot(:)
		integer,allocatable :: Nrot_curve(:)
		character(len=256) :: filename

		print *, 'PREPARING GRIDS FOR DECOMPOSITION ACCORDING TO PARAMETER FILE'

		if( pre_smooth ) then
			print '(a13)', 'pre smoothing'
			print '(a5,i4.4,a11,i1,a1)', 'FWHM=', pre_smooth_FWHM, ' cells, in ', pre_smooth_dim, 'D'
			if( pre_smooth_dim .eq. 3 ) then
				call fft_smooth_3D(pre_smooth_FWHM, pertx, Nxf, Nyf, Nzf, 1, 2)
				call fft_smooth_3D(pre_smooth_FWHM, perty, Nxf, Nyf, Nzf, 1, 2)
				call fft_smooth_3D(pre_smooth_FWHM, pertz, Nxf, Nyf, Nzf, 1, 2)
			elseif( pre_smooth_dim .eq. 2 ) then
				do m=1,Nzf
					call fft_smooth_2D(pre_smooth_FWHM, pertx(:,:,m), Nxf, Nyf, 1, 2)
					call fft_smooth_2D(pre_smooth_FWHM, perty(:,:,m), Nxf, Nyf, 1, 2)
					call fft_smooth_2D(pre_smooth_FWHM, pertz(:,:,m), Nxf, Nyf, 1, 2)
				end do
			end if
		end if

		if( decomp_type .eq. 2 ) then
			print *, 'subtracting rotation curve'
			print *, 'rotation calculated from 2d map integrated along Z axis'
			Nrot = ceiling( sqrt(real(Nxf/2)**2 + real(Nyf/2)**2) ) - 1
			rx = real(Nxf+1)/2.0_8
			ry = real(Nyf+1)/2.0_8

			Nvert = 5
			allocate( Vx(Nxf,Nyf), Vy(Nxf,Nyf), vrot_curve(Nrot,Nvert), rrot(Nrot), Nrot_curve(Nrot) )
			if( sphererad .gt. 0.0_8 .and. Ebox .eq. 0 ) then
				if( cylheigh .gt. 0.0_8 ) then
					rad = real(Nzf) * cylheigh
				else
					rad = real(Nzf) * sphererad
				end if
			else
				rad = real(Nzf)
			end if
			Vx = sum(pertx, dim=3) / rad
			Vy = sum(perty, dim=3) / rad
			vrot_curve = 0.0_8
			rrot = 0.0_8
			Nrot_curve = 0
			k = ceiling(sngl(Nzf)/sngl(Nvert))
			do j=1,Nyf
				do i=1,Nxf
					rad = sqrt( real(i-rx)**2 + real(j-ry)**2 )
					m = min( ceiling(rad), Nrot )
					do n=1,Nvert
						n1 = min(Nzf, n*k)
						vxt = sum( pertx(i,j,(1+(n-1)*k):n1) ) / real(n1 - (n-1)*k)
						vyt = sum( perty(i,j,(1+(n-1)*k):n1) ) / real(n1 - (n-1)*k)
						vrot_curve(m,n) = vrot_curve(m,n) + vxt*real(ry-j)/rad + vyt*real(i-rx)/rad
					end do
!					rrot(m) = rrot(m) + Vx(i,j)*real(ry-j) + Vy(i,j)*real(i-rx)
					rrot(m) = rrot(m) + rad
					Nrot_curve(m) = Nrot_curve(m) + 1
				end do
			end do
!			where( Nrot_curve .gt. 0 ) rrot = rrot / vrot_curve
			where( Nrot_curve .gt. 0 ) rrot = rrot / Nrot_curve
			do i=1,Nvert
				where( Nrot_curve .gt. 0 ) vrot_curve(:,i) = vrot_curve(:,i) / Nrot_curve(:)
			end do
			write(filename,'(a,a)') trim(output_dirname),'/rotation_curve.out'
			open(unit=20,file=filename)
			if( Nvert .lt. 9 ) then
				write(filename,'(a,i1,a)') '(',Nvert+1,'(1x,es12.5))'
			elseif( Nvert .lt. 99 ) then
				write(filename,'(a,i2,a)') '(',Nvert+1,'(1x,es12.5))'
			else
				write(filename,'(a,i3,a)') '(',Nvert+1,'(1x,es12.5))'
			end if
			do i=1,Nrot
				write(20,trim(filename)) rrot(i), vrot_curve(1:Nvert,i)
			end do
			close(unit=20)
			n1 = ceiling(sngl(Nzf)/sngl(Nvert))
			do j=1,Nyf
				do i=1,Nxf
					rad = sqrt( real(i-rx)**2 + real(j-ry)**2 )
					m = min( ceiling(rad), Nrot )
					do k=1,Nzf
						n = max( 1, ceiling(real(k)/real(n1)) )
						n = min( n, Nvert )
						if( m .eq. 1 .or. m .eq. Nrot ) then
							vrot_interp = vrot_curve(m,n)
						else
							rm = real(m)
							if( rad .lt. rrot(m) ) then
								vrot_interp = log10( vrot_curve(m-1,n) ) + log10( rad/rrot(m-1) ) * log10( vrot_curve(m,n)/vrot_curve(m-1,n) ) / &
									&     log10( rrot(m)/rrot(m-1) )
							else
								vrot_interp = log10( vrot_curve(m,n) ) + log10( rad/rrot(m) ) * log10( vrot_curve(m+1,n)/vrot_curve(m,n) ) / &
									&     log10( rrot(m+1)/rrot(m) )
							end if
							vrot_interp = 10**(vrot_interp)
						end if
						pertx(i,j,k) = pertx(i,j,k) - vrot_interp * real(ry-j)/rad
						perty(i,j,k) = perty(i,j,k) - vrot_interp * real(i-rx)/rad
					end do
				end do
			end do
			deallocate( Vx, Vy, vrot_curve, rrot, Nrot_curve )
		elseif( decomp_type .eq. 3 ) then
			print *, 'Subtracting from each cell the velocity averaged over some scale, Gaussian weighted'
			print *, 'Equivalent to filtering out large scales'
			print '(a5,i4.4,a13)', 'FWHM=', type_3_FWHM, ' cells, in 3D'
			call fft_smooth_3D(type_3_FWHM, pertx, Nxf, Nyf, Nzf, 0, 2)
			call fft_smooth_3D(type_3_FWHM, perty, Nxf, Nyf, Nzf, 0, 2)
			call fft_smooth_3D(type_3_FWHM, pertz, Nxf, Nyf, Nzf, 0, 2)
		end if

		print *, 'Making surface maps of pre-decomp turbulence components in cylindrical coordinates'
		allocate( Vx(Nxf,Nyf), Vy(Nxf,Nyf), Vz(Nxf,Nyf), Vr(Nxf,Nyf), Vphi(Nxf,Nyf) )
		if( sphererad .gt. 0.0_8 .and. Ebox .eq. 0 ) then
			if( cylheigh .gt. 0.0_8 ) then
				rad = real(Nzf) * cylheigh
			else
				rad = real(Nzf) * sphererad
			end if
		else
			rad = real(Nzf)
		end if
		Vx = sum(pertx, dim=3) / rad
		Vy = sum(perty, dim=3) / rad
		Vz = sum(pertz, dim=3) / rad
		rx = real(Nxf+1)/2.0_8
		ry = real(Nyf+1)/2.0_8
		do j=1,Nyf
			do i=1,Nxf
				rad = sqrt( real(i-rx)**2 + real(j-ry)**2 )
				Vr(i,j)   = Vx(i,j)*real(i-rx)/rad + Vy(i,j)*real(j-ry)/rad
				Vphi(i,j) = Vx(i,j)*real(ry-j)/rad + Vy(i,j)*real(i-rx)/rad
			end do
		end do
		Vx = sqrt( sum(pertx**2 + perty**2 + pertz**2, dim=3) / rad )
		write(filename,'(a,a)') trim(output_dirname),'/pre_decomp_turbulence.bin'
		open(unit=20,file=filename,form='unformatted')
		write(20) Nxf, Nyf, Vr, Vphi, Vz, Vx
		close(unit=20)
		deallocate( Vx, Vy, Vz, Vr, Vphi )

	end subroutine turb_def
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine turb_decomp()
	! Final part
	! decomposes turbulence into compreesive and solenoidal components, calculates the global ratios, and creates projected 2D plots of them
	! snapshot is snapshot number from main loop, needed for io purposes
		implicit none
		integer :: i, j, k, m
		real(8) :: Ei, fsol, rad
		real(8),allocatable :: Sigc_2d(:,:,:), Sigs_2d(:,:), dot_2d(:,:), test(:,:,:)
		character(len=256) :: filename

		! Define the field that will be decomposed
		call turb_def()
		print *, 'Starting decomposition'
		Ei = sum(pertx**2 + perty**2 + pertz**2) / (Nxf*Nyf*Nzf)
		print '(a3,1x,es12.5)', 'RMS', sqrt(Ei)

		if( fourier_decomp ) then
			call fft_decompose

			allocate( test(Nxf, Nyf, Nzf) )
			test = 0.0_8
			pertx = 0.0_8
			perty = 0.0_8
			test(2:(Nxf-1), 2:(Nyf-1), 2:(Nzf-1)) = sigcx(3:Nxf, 2:(Nyf-1), 2:(Nzf-1)) - sigcx(1:(Nxf-2), 2:(Nyf-1), 2:(Nzf-1)) + &
							&       sigcy(2:(Nxf-1), 3:Nyf, 2:(Nzf-1)) - sigcy(2:(Nxf-1), 1:(Nyf-2), 2:(Nzf-1)) + &
							&       sigcz(2:(Nxf-1), 2:(Nyf-1), 3:Nzf) - sigcz(2:(Nxf-1), 2:(Nyf-1), 1:(Nzf-2))
			where( test .ge. 0.0_8 ) pertx = sigcx**2 + sigcy**2 + sigcz**2
			where( test .lt. 0.0_8 ) perty = sigcx**2 + sigcy**2 + sigcz**2
			sigcz = 2.0_8 * (sigcx*sigsx + sigcy*sigsy + sigcz*sigsz)
			sigcx = pertx
			sigcy = perty
			sigsx = sigsx**2 + sigsy**2 + sigsz**2
			deallocate( test )
		else
			call local_decompose
		end if
		deallocate( pertx, perty, pertz, sigsy, sigsz )
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		print *, 'making surface plots of turbulence modes'
		allocate( Sigc_2d(Nxf,Nyf,2), Sigs_2d(Nxf,Nyf), dot_2d(Nxf,Nyf) )
		if( sphererad .gt. 0.0_8 .and. Ebox .eq. 0 ) then
			if( cylheigh .gt. 0.0_8 ) then
				rad = real(Nzf) * cylheigh
			else
				rad = real(Nzf) * sphererad
			end if
		else
			rad = real(Nzf)
		end if
		Sigc_2d(:,:,1) = sum( sigcx, dim=3 ) / rad
		Sigc_2d(:,:,2) = sum( sigcy, dim=3 ) / rad
		Sigs_2d(:,:)   = sum( sigsx, dim=3 ) / rad
		dot_2d(:,:)    = sum( sigcz, dim=3 ) / rad

		write(filename,'(a,a)') trim(output_dirname),'/decomposed_turbulence.bin'
		open(unit=20,file=filename,form='unformatted')
		write(20) Nxf, Nyf, sqrt( Sigc_2d(:,:,1) ), sqrt( Sigc_2d(:,:,2) ), sqrt( sigs_2d )
		close(unit=20)
		write(filename,'(a,a)') trim(output_dirname),'/mode_ratios.bin'
		open(unit=20,file=filename,form='unformatted')
		write(20) Nxf, Nyf, Sigc_2d(:,:,1) / Sigs_2d, Sigc_2d(:,:,2) / Sigs_2d, dot_2d / ( Sigc_2d(:,:,1) + Sigc_2d(:,:,2) + Sigs_2d + dot_2d )
		close(unit=20)
		deallocate( Sigc_2d, Sigs_2d, dot_2d )

		print *, 'conservation: w/ and wo/ cross term'
		print '(2(1x,es12.5))', sum( sigcx + sigcy + sigsx + sigcz ) / (Ei*Nxf*Nyf*Nzf), sum( sigcx + sigcy + sigsx ) / (Ei*Nxf*Nyf*Nzf)

		print *, 'total power real space: s/t, c/t, e/t, dot/t'
		Ei = sum(sigcx + sigcy + sigsx + sigcz)
		if( f_solenoidal .ge. 0.0_8 .and. f_solenoidal .le. 1.0_8 ) then
			fsol = (2.0_8*f_solenoidal**2) / ( 2.0_8*f_solenoidal**2 + (1.0_8-f_solenoidal)**2 )
		else
			fsol = 2.0_8 / 3.0_8
		end if
		print '(5(1x,es12.5))', sum(sigsx)/Ei, sum(sigcy)/Ei, sum(sigcx)/Ei, sum(sigcz)/Ei, fsol
		allocate( test(Nxf, Nyf, Nzf) )
		test = sigcx + sigcy + sigsx + sigcz
		m  = count( mask = test .gt. 1.d-6 )
		print '(5(1x,es12.5))',   sum(sigsx/test, mask=test .gt. 1.d-6)/m, sum(sigcy/test, mask=test .gt. 1.d-6)/m, sum(sigcx/test, mask=test .gt. 1.d-6)/m, & 
					& sum(sigcz/test, mask=test .gt. 1.d-6)/m, fsol

		deallocate( sigcx, sigcy, sigcz, sigsx, test )

	end subroutine turb_decomp
!_____________________________________________________________________________________________________________________________________________________________________
	subroutine fft_decompose
	use MKL_DFTI
	! decomposes turbulence into solenoidal and compressive modes in Fourier space

		implicit none
		real(8),allocatable :: Skmag(:,:,:), test(:,:,:), test1d(:), kvec(:)
		complex(8),allocatable :: fft_inx(:), fft_iny(:), fft_inz(:), PS_out(:,:)
		complex(8),allocatable :: Gcx_1d(:),  Gcy_1d(:), Gcz_1d(:), Gsx_1d(:),  Gsy_1d(:), Gsz_1d(:)
		complex(8) :: im, Gk
		real(8) :: normc, norms, dk1, dk2, dk3
		integer :: i, j, m, n, n1, n2, n3, n4, nvec, error, length(3), Nxy, Nxyz, Nxh, Nyh, Nzh
		character(len=256) :: filename
		type(DFTI_DESCRIPTOR), POINTER :: fft_desc

		print *, 'FFT Decompose'
		Nxy = Nxf*Nyf
		Nxyz = Nxy*Nzf
		im = (0.0_8, 1.0_8)
		length = (/ Nxf, Nyf, Nzf /)

		allocate( kx(Nxf,Nyf,Nzf), ky(Nxf,Nyf,Nzf), kz(Nxf,Nyf,Nzf), kmag(Nxf,Nyf,Nzf) )
		kx = 0.0_8
		ky = 0.0_8
		kz = 0.0_8
		kmag = 0.0_8
		Nxh = Nxf/2
		Nyh = Nyf/2
		Nzh = Nzf/2
		do m=1,Nzf
			do j=1,Nyf
				kx(1:Nxh, j, m)       = real( (/ (i,i=1,Nxh) /) - 1 )
				kx((Nxh+1):Nxf, j, m) = real( (/ (i,i=1,Nxh) /) - 1 - Nxh )
			end do
		end do
		do m=1,Nzf
			do i=1,Nxf
				ky(i, 1:Nyh, m)       = real( (/ (j,j=1,Nyh) /) - 1 )
				ky(i, (Nyh+1):Nyf, m) = real( (/ (j,j=1,Nyh) /) - 1 - Nyh )
			end do
		end do
		do j=1,Nyf
			do i=1,Nxf
				kz(i, j, 1:Nzh)       = real( (/ (m,m=1,Nzh) /) - 1 )
				kz(i, j, (Nzh+1):Nzf) = real( (/ (m,m=1,Nzh) /) - 1 - Nzh )
			end do
		end do
		kx = kx * 2.0_8*pi/Nxf
		ky = ky * 2.0_8*pi/Nyf
		kz = kz * 2.0_8*pi/Nzf
		kmag = sqrt( kx**2 + ky**2 + kz**2 )

		allocate( fft_inx(Nxyz), fft_iny(Nxyz), fft_inz(Nxyz) )
		do m=1,Nzf
			do j=1,Nyf
				do i=1,Nxf
					fft_inx( (m-1)*Nxy + (j-1)*Nxf + i ) = pertx(i,j,m)
					fft_iny( (m-1)*Nxy + (j-1)*Nxf + i ) = perty(i,j,m)
					fft_inz( (m-1)*Nxy + (j-1)*Nxf + i ) = pertz(i,j,m)
				end do
			end do
		end do
		norms = sum(abs(fft_inx)**2 + abs(fft_iny)**2 + abs(fft_inz)**2)

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

		!!! Filter out modes that were not initiated, i.e. larger than the disc scale height and smaller than ~2 cells
		dk1 = 2.0_8*pi/real(Nx0)
		do m=1,Nzf
			do j=1,Nyf
				do i=1,Nxf
					if( kmag(i,j,m) .lt. real(kmin)*dk1 .or. kmag(i,j,m) .ge. real(kmax+1)*dk1 ) then
						n = (m-1)*Nxy + (j-1)*Nxf + i
						fft_inx(n) = 0.0_8
						fft_iny(n) = 0.0_8
						fft_inz(n) = 0.0_8
					end if
				end do
			end do
		end do

		if( use_sine ) then
			allocate( Skmag(Nxf,Nyf,Nzf) )
			Skmag = kmag
			kx = sin(kx)
			ky = sin(ky)
			kz = sin(kz)
			kmag = sqrt( kx**2 + ky**2 + kz**2 )
		end if
		allocate( Gcx_1d(Nxyz), Gcy_1d(Nxyz), Gcz_1d(Nxyz), Gsx_1d(Nxyz), Gsy_1d(Nxyz), Gsz_1d(Nxyz) )
		! compressive part
		do m=1,Nzf
			do j=1,Nyf
				do i=1,Nxf
					if( kmag(i,j,m) .eq. 0.0_8 ) then
						Gcx_1d = 0.0_8
						Gcy_1d = 0.0_8
						Gcz_1d = 0.0_8
					else
						n = (m-1)*Nxy + (j-1)*Nxf + i
						Gk = ( kx(i,j,m)*fft_inx(n) + ky(i,j,m)*fft_iny(n) + kz(i,j,m)*fft_inz(n) ) / max(kmag(i,j,m)**2, 1.d-16)
						Gcx_1d(n) = kx(i,j,m) * Gk
						Gcy_1d(n) = ky(i,j,m) * Gk
						Gcz_1d(n) = kz(i,j,m) * Gk
					end if
				end do
			end do
		end do
		normc = sqrt( sum(abs(Gcx_1d)**2 + abs(Gcy_1d)**2 + abs(Gcz_1d)**2) )
		deallocate( kx, ky, kz )

		! solenoidal part
		Gsx_1d = fft_inx - Gcx_1d
		Gsy_1d = fft_iny - Gcy_1d
		Gsz_1d = fft_inz - Gcz_1d
		norms = sqrt( sum(abs(Gsx_1d)**2 + abs(Gsy_1d)**2 + abs(Gsz_1d)**2) )

		if( use_sine ) then 
			kmag = Skmag
			deallocate( Skmag )
		end if

		print *, 'Final power spectrum'
		dk1 = 2.0_8*pi/Nxf
		dk2 = 2.0_8*pi/Nyf
		dk3 = 2.0_8*pi/Nzf
		n1 = ceiling( dk2/dk1 )
		n2 = ceiling( dk3/dk2 )
		n3 = ceiling( sqrt(3.0_8)*pi / dk3 ) + 1
		nvec = n1 + n2 + n3 - 2
		print *, 'n1, n2, n3, nvec ',n1, n2, n3, nvec
		allocate( kvec(nvec), PS_out(nvec,13) )
		kvec(1:n1)  = (/ (i,i=1,n1) /) * dk1
		kvec(n1:(n1+n2-1)) = (/ (i,i=1,n2) /) * dk2
		kvec((n1+n2-1):nvec) = (/ (i,i=1,n3) /) * dk3
		PS_out = 0.0_8
		do m=1,Nzf
			do j=1,Nyf
				do i=1,Nxf
					if( kmag(i,j,m) .ge. dk3 ) then
						n4 = min( floor( kmag(i,j,m)/dk3 ) + n1 + n2 - 1, nvec )
					elseif( kmag(i,j,m) .ge. dk2 ) then
						n4 = floor( kmag(i,j,m)/dk2 ) + n1
					else
						n4 = floor( kmag(i,j,m)/dk1 ) + 1
					end if
					n = (m-1)*Nxy + (j-1)*Nxf + i
					PS_out(n4,1)  = PS_out(n4,1)  + abs(fft_inx(n))**2 + abs(fft_iny(n))**2 + abs(fft_inz(n))**2
					PS_out(n4,2)  = PS_out(n4,2)  + abs(Gsx_1d(n))**2  + abs(Gsy_1d(n))**2  + abs(Gsz_1d(n))**2
					PS_out(n4,3)  = PS_out(n4,3)  + abs(Gcx_1d(n))**2  + abs(Gcy_1d(n))**2  + abs(Gcz_1d(n))**2
					PS_out(n4,4)  = PS_out(n4,4)  + fft_inx(n)
					PS_out(n4,5)  = PS_out(n4,5)  + fft_iny(n)
					PS_out(n4,6)  = PS_out(n4,6)  + fft_inz(n)
					PS_out(n4,7)  = PS_out(n4,7)  + Gsx_1d(n)
					PS_out(n4,8)  = PS_out(n4,8)  + Gsy_1d(n)
					PS_out(n4,9)  = PS_out(n4,9)  + Gsz_1d(n)
					PS_out(n4,10) = PS_out(n4,10) + Gcx_1d(n)
					PS_out(n4,11) = PS_out(n4,11) + Gcy_1d(n)
					PS_out(n4,12) = PS_out(n4,12) + Gcz_1d(n)
					PS_out(n4,13) = PS_out(n4,13) + 1.0_8
				end do
			end do
		end do
		deallocate( kmag )
		write(filename,'(a,a)') trim(output_dirname),'/final_power_spectra.out'
		open(unit=20,file=filename)
		do i=1,n1
			write(20,'(7(1x,es12.5))') kvec(i)-0.5_8*dk1, real(PS_out(i,1)), real(PS_out(i,2)), real(PS_out(i,3)), & 
				& abs(PS_out(i,4)/PS_out(i,13))**2  + abs(PS_out(i,5)/PS_out(i,13))**2  + abs(PS_out(i,6)/PS_out(i,13))**2, &
				& abs(PS_out(i,7)/PS_out(i,13))**2  + abs(PS_out(i,8)/PS_out(i,13))**2  + abs(PS_out(i,9)/PS_out(i,13))**2, & 
				& abs(PS_out(i,10)/PS_out(i,13))**2 + abs(PS_out(i,11)/PS_out(i,13))**2 + abs(PS_out(i,12)/PS_out(i,13))**2
		end do
		if( n1+n2-1 .gt. 1 ) then
			do i=n1+1,n1+n2-1
				write(20,'(7(1x,es12.5))') kvec(i)-0.5_8*dk2, real(PS_out(i,1)), real(PS_out(i,2)), real(PS_out(i,3)), & 
					& abs(PS_out(i,4)/PS_out(i,13))**2  + abs(PS_out(i,5)/PS_out(i,13))**2  + abs(PS_out(i,6)/PS_out(i,13))**2, &
					& abs(PS_out(i,7)/PS_out(i,13))**2  + abs(PS_out(i,8)/PS_out(i,13))**2  + abs(PS_out(i,9)/PS_out(i,13))**2, & 
					& abs(PS_out(i,10)/PS_out(i,13))**2 + abs(PS_out(i,11)/PS_out(i,13))**2 + abs(PS_out(i,12)/PS_out(i,13))**2
			end do
		end if
		do i=n1+n2,nvec
			write(20,'(7(1x,es12.5))') kvec(i)-0.5_8*dk3, real(PS_out(i,1)), real(PS_out(i,2)), real(PS_out(i,3)), & 
				& abs(PS_out(i,4)/PS_out(i,13))**2  + abs(PS_out(i,5)/PS_out(i,13))**2  + abs(PS_out(i,6)/PS_out(i,13))**2, &
				& abs(PS_out(i,7)/PS_out(i,13))**2  + abs(PS_out(i,8)/PS_out(i,13))**2  + abs(PS_out(i,9)/PS_out(i,13))**2, & 
				& abs(PS_out(i,10)/PS_out(i,13))**2 + abs(PS_out(i,11)/PS_out(i,13))**2 + abs(PS_out(i,12)/PS_out(i,13))**2
		end do
		close(unit=20)
		print *, 'done'
		deallocate( kvec, PS_out )

		print *, 'Gcx_1d'
		error = DftiCreateDescriptor(fft_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc)
		error = DftiComputeBackward(fft_desc, Gcx_1d)
		error = DftiFreeDescriptor(fft_desc)
		print '(4(1x,es12.5))', maxval(real(Gcx_1d)), minval(real(Gcx_1d)), maxval(aimag(Gcx_1d)), minval(aimag(Gcx_1d))

		print *, 'Gcy_1d'
		error = DftiCreateDescriptor(fft_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc)
		error = DftiComputeBackward(fft_desc, Gcy_1d)
		error = DftiFreeDescriptor(fft_desc)
		print '(4(1x,es12.5))', maxval(real(Gcy_1d)), minval(real(Gcy_1d)), maxval(aimag(Gcy_1d)), minval(aimag(Gcy_1d))

		print *, 'Gcz_1d'
		error = DftiCreateDescriptor(fft_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc)
		error = DftiComputeBackward(fft_desc, Gcz_1d)
		error = DftiFreeDescriptor(fft_desc)
		print '(4(1x,es12.5))', maxval(real(Gcz_1d)), minval(real(Gcz_1d)), maxval(aimag(Gcz_1d)), minval(aimag(Gcz_1d))

		print *, 'Gsx_1d'
		error = DftiCreateDescriptor(fft_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc)
		error = DftiComputeBackward(fft_desc, Gsx_1d)
		error = DftiFreeDescriptor(fft_desc)
		print '(4(1x,es12.5))', maxval(real(Gsx_1d)), minval(real(Gsx_1d)), maxval(aimag(Gsx_1d)), minval(aimag(Gsx_1d))

		print *, 'Gsy_1d'
		error = DftiCreateDescriptor(fft_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc)
		error = DftiComputeBackward(fft_desc, Gsy_1d)
		error = DftiFreeDescriptor(fft_desc)
		print '(4(1x,es12.5))', maxval(real(Gsy_1d)), minval(real(Gsy_1d)), maxval(aimag(Gsy_1d)), minval(aimag(Gsy_1d))

		print *, 'Gsz_1d'
		error = DftiCreateDescriptor(fft_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc)
		error = DftiComputeBackward(fft_desc, Gsz_1d)
		error = DftiFreeDescriptor(fft_desc)
		print '(4(1x,es12.5))', maxval(real(Gsz_1d)), minval(real(Gsz_1d)), maxval(aimag(Gsz_1d)), minval(aimag(Gsz_1d))

		Gcx_1d = Gcx_1d / sqrt(real(Nxyz))
		Gcy_1d = Gcy_1d / sqrt(real(Nxyz))
		Gcz_1d = Gcz_1d / sqrt(real(Nxyz))
		Gsx_1d = Gsx_1d / sqrt(real(Nxyz))
		Gsy_1d = Gsy_1d / sqrt(real(Nxyz))
		Gsz_1d = Gsz_1d / sqrt(real(Nxyz))

		allocate( sigcx(Nxf,Nyf,Nzf), sigcy(Nxf,Nyf,Nzf), sigcz(Nxf,Nyf,Nzf), sigsx(Nxf,Nyf,Nzf), sigsy(Nxf,Nyf,Nzf), sigsz(Nxf,Nyf,Nzf) )
		do m=1,Nzf
			do j=1,Nyf
				do i=1,Nxf
					sigcx(i,j,m) = real(Gcx_1d( (m-1)*Nxy + (j-1)*Nxf + i ))
					sigcy(i,j,m) = real(Gcy_1d( (m-1)*Nxy + (j-1)*Nxf + i ))
					sigcz(i,j,m) = real(Gcz_1d( (m-1)*Nxy + (j-1)*Nxf + i ))
					sigsx(i,j,m) = real(Gsx_1d( (m-1)*Nxy + (j-1)*Nxf + i ))
					sigsy(i,j,m) = real(Gsy_1d( (m-1)*Nxy + (j-1)*Nxf + i ))
					sigsz(i,j,m) = real(Gsz_1d( (m-1)*Nxy + (j-1)*Nxf + i ))
				end do
			end do
		end do
		deallocate( Gcx_1d, Gcy_1d, Gcz_1d, Gsx_1d, Gsy_1d, Gsz_1d )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		print *, '!!!!!'
		print *, 'TESTS'
		allocate( test(Nxf-2, Nyf-2, Nzf-2) )
		m = (Nxf-2)*(Nyf-2)*(Nzf-2)
		print *, 'divergence of solonoidal mode'
		test = 0.5_8 * (  sigsx(3:Nxf,     2:(Nyf-1), 2:(Nzf-1)) - sigsx(1:(Nxf-2), 2:(Nyf-1), 2:(Nzf-1)) + &
				& sigsy(2:(Nxf-1), 3:Nyf,     2:(Nzf-1)) - sigsy(2:(Nxf-1), 1:(Nyf-2), 2:(Nzf-1)) + &
				& sigsz(2:(Nxf-1), 2:(Nyf-1), 3:Nzf)     - sigsz(2:(Nxf-1), 2:(Nyf-1), 1:(Nzf-2)) )
		print '(3(1x,es12.5))', maxval(abs(test)), minval(abs(test)), sum(abs(test))/m
		print *, 'divergence of compressive mode'
		test = 0.5_8*(    sigcx(3:Nxf,     2:(Nyf-1), 2:(Nzf-1)) - sigcx(1:(Nxf-2), 2:(Nyf-1), 2:(Nzf-1)) + &
				& sigcy(2:(Nxf-1), 3:Nyf,     2:(Nzf-1)) - sigcy(2:(Nxf-1), 1:(Nyf-2), 2:(Nzf-1)) + &
				& sigcz(2:(Nxf-1), 2:(Nyf-1), 3:Nzf)     - sigcz(2:(Nxf-1), 2:(Nyf-1), 1:(Nzf-2)) )
		print '(3(1x,es12.5))', maxval(abs(test)), minval(abs(test)), sum(abs(test))/m
		print *, 'rotor of solonoidal mode'
		test = 0.5_8 * sqrt( &
				&  ( (sigsz(2:(Nxf-1), 3:Nyf, 2:(Nzf-1)) - sigsz(2:(Nxf-1), 1:(Nyf-2), 2:(Nzf-1))) - & 
				&    (sigsy(2:(Nxf-1), 2:(Nyf-1), 3:Nzf) - sigsy(2:(Nxf-1), 2:(Nyf-1), 1:(Nzf-2))) )**2 + &
				&  ( (sigsx(2:(Nxf-1), 2:(Nyf-1), 3:Nzf) - sigsx(2:(Nxf-1), 2:(Nyf-1), 1:(Nzf-2))) - & 
				&    (sigsz(3:Nxf, 2:(Nyf-1), 2:(Nzf-1)) - sigsz(1:(Nxf-2), 2:(Nyf-1), 2:(Nzf-1))) )**2 + &
				&  ( (sigsy(3:Nxf, 2:(Nyf-1), 2:(Nzf-1)) - sigsy(1:(Nxf-2), 2:(Nyf-1), 2:(Nzf-1))) - & 
				&    (sigsx(2:(Nxf-1), 3:Nyf, 2:(Nzf-1)) - sigsx(2:(Nxf-1), 1:(Nyf-2), 2:(Nzf-1))) )**2   &
				&  )
		print '(3(1x,es12.5))', maxval(abs(test)), minval(abs(test)), sum(abs(test))/m
		print *, 'rotor of compressive mode'
		test = 0.5_8 * sqrt( &
				&  ( (sigcz(2:(Nxf-1), 3:Nyf, 2:(Nzf-1)) - sigcz(2:(Nxf-1), 1:(Nyf-2), 2:(Nzf-1))) - & 
				&    (sigcy(2:(Nxf-1), 2:(Nyf-1), 3:Nzf) - sigcy(2:(Nxf-1), 2:(Nyf-1), 1:(Nzf-2))) )**2 + &
				&  ( (sigcx(2:(Nxf-1), 2:(Nyf-1), 3:Nzf) - sigcx(2:(Nxf-1), 2:(Nyf-1), 1:(Nzf-2))) - & 
				&    (sigcz(3:Nxf, 2:(Nyf-1), 2:(Nzf-1)) - sigcz(1:(Nxf-2), 2:(Nyf-1), 2:(Nzf-1))) )**2 + &
				&  ( (sigcy(3:Nxf, 2:(Nyf-1), 2:(Nzf-1)) - sigcy(1:(Nxf-2), 2:(Nyf-1), 2:(Nzf-1))) - & 
				&    (sigcx(2:(Nxf-1), 3:Nyf, 2:(Nzf-1)) - sigcx(2:(Nxf-1), 1:(Nyf-2), 2:(Nzf-1))) )**2   &
				&  )
		print '(3(1x,es12.5))', maxval(abs(test)), minval(abs(test)), sum(abs(test))/m
		deallocate( test )

		allocate( test(Nxf, Nyf, Nzf) )
		print *, 'difference magnitude from full field'
		test = sqrt( (pertx-sigsx-sigcx)**2 + (perty-sigsy-sigcy)**2 + (pertz-sigsz-sigcz)**2 )
		print '(3(1x,es12.5))', maxval(abs(test)), minval(abs(test)), sum(abs(test))/Nxyz

		print *, 'dot product'
		test = 2.0_8*( sigcx*sigsx + sigcy*sigsy + sigcz*sigsz )
		print '(4(1x,es12.5))', maxval(abs(test)), minval(abs(test)), sum(abs(test))/Nxyz, sum(test)/Nxyz

		print *, 'total power fourier space: t, s/t c/t'
		print '(3(1x,es12.5))', norms**2+normc**2, norms**2 / (norms**2+normc**2), normc**2 / (norms**2+normc**2)
		print *, 'total power real space: t, s/t c/t, dot/t'
		norms = sum(sigsx**2 + sigsy**2 + sigsz**2)
		normc = sum(sigcx**2 + sigcy**2 + sigcz**2)
		print '(4(1x,es12.5))', norms+normc+sum(test), norms/(norms+normc+sum(test)), normc/(norms+normc+sum(test)), sum(test)/(norms+normc+sum(test))
		print *, 'average power ratios real space: s/t c/t, dot/t'
		test = test + sigsx**2 + sigsy**2 + sigsz**2 + sigcx**2 + sigcy**2 + sigcz**2
		m = count( mask = test .gt. 1.d-6 )
		print '(3(1x,es12.5))',   sum( (sigsx**2 + sigsy**2 + sigsz**2)/test, mask = test .gt. 1d-6 )/m, sum( (sigcx**2 + sigcy**2 + sigcz**2)/test, mask = test .gt. 1d-6 )/m, & 
					& sum( 2.0_8*(sigcx*sigsx + sigcy*sigsy + sigcz*sigsz)/test, mask = test .gt. 1d-6 )/m


		deallocate( test )
	end subroutine fft_decompose
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine local_decompose
	! decomposes turbulence into solenoidal and compressive modes locally in real space (using strain rate tensor)

		implicit none
		real(8) :: local_sigcx(3,3,3),local_sigcy(3,3,3),local_sigcz(3,3,3),local_sigsx(3,3,3),local_sigsy(3,3,3),local_sigsz(3,3,3)
		real(8) :: divx, divy, divz
		integer :: i, j, m

		print *, 'Decomposing using strain rate tensor'
		allocate( sigcx(Nxf, Nyf, Nzf), sigcy(Nxf, Nyf, Nzf), sigcz(Nxf, Nyf, Nzf), sigsx(Nxf, Nyf, Nzf), sigsy(Nxf, Nyf, Nzf), sigsz(Nxf, Nyf, Nzf) )
		sigcx = 0.0_8
		sigcy = 0.0_8
		sigcz = 0.0_8
		sigsx = 0.0_8
		sigsy = 0.0_8
		sigsz = 0.0_8

		local_sigcx(:,:,:) = 0.0_8
		local_sigcy(:,:,:) = 0.0_8
		local_sigcz(:,:,:) = 0.0_8
		local_sigsx(:,:,:) = 0.0_8
		local_sigsy(:,:,:) = 0.0_8
		local_sigsz(:,:,:) = 0.0_8
		if( local_version .eq. 1 ) then
			do m=2,Nzf-1
				do j=2,Nyf-1
					do i=2,Nxf-1
						divx = (pertx(i+1,j,m) - pertx(i-1,j,m)) / 2.0_8
						divy = (perty(i,j+1,m) - perty(i,j-1,m)) / 2.0_8
						divz = (pertz(i,j,m+1) - pertz(i,j,m-1)) / 2.0_8
						local_sigcx(1, 1:3, 1:3) = -1.0_8 * divx
						local_sigcx(3, 1:3, 1:3) = divx
						local_sigcy(1:3, 1, 1:3) = -1.0_8 * divy
						local_sigcy(1:3, 3, 1:3) = divy
						local_sigcz(1:3, 1:3, 1) = -1.0_8 * divz
						local_sigcz(1:3, 1:3, 3) = divz

						local_sigsx(1:3, 1:3, 1:3) = pertx(i-1:i+1, j-1:j+1, m-1:m+1) - pertx(i, j, m) - local_sigcx(1:3, 1:3, 1:3)
						local_sigsy(1:3, 1:3, 1:3) = perty(i-1:i+1, j-1:j+1, m-1:m+1) - perty(i, j, m) - local_sigcy(1:3, 1:3, 1:3)
						local_sigsz(1:3, 1:3, 1:3) = pertz(i-1:i+1, j-1:j+1, m-1:m+1) - pertz(i, j, m) - local_sigcz(1:3, 1:3, 1:3)

						if( divx + divy + divz .ge. 0.0_8 ) then
							sigcx(i,j,m) = sum( local_sigcx**2 + local_sigcy**2 + local_sigcz**2 ) / 27.0_8
						else
							sigcy(i,j,m) = sum( local_sigcx**2 + local_sigcy**2 + local_sigcz**2 ) / 27.0_8
						end if
						sigsx(i,j,m) = sum( local_sigsx**2 + local_sigsy**2 + local_sigsz**2 ) / 27.0_8
						sigcz(i,j,m) = 2.0_8 * sum( local_sigcx*local_sigsx + local_sigcy*local_sigsy + local_sigcz*local_sigsz ) / 27.0_8
					end do
				end do
			end do
		elseif( local_version .eq. 2 ) then
			do m=2,Nzf-1
				do j=2,Nyf-1
					do i=2,Nxf-1
						divx = (pertx(i+1,j,m)-pertx(i-1,j,m) + perty(i,j+1,m)-perty(i,j-1,m) + pertz(i,j,m+1)-pertz(i,j,m-1)) / 6.0_8
						local_sigcx(1, 1:3, 1:3) = -1.0_8 * divx
						local_sigcx(3, 1:3, 1:3) = divx
						local_sigcy(1:3, 1, 1:3) = -1.0_8 * divx
						local_sigcy(1:3, 3, 1:3) = divx
						local_sigcz(1:3, 1:3, 1) = -1.0_8 * divx
						local_sigcz(1:3, 1:3, 3) = divx

						local_sigsx(1:3, 1:3, 1:3) = pertx(i-1:i+1, j-1:j+1, m-1:m+1) - pertx(i, j, m) - local_sigcx(1:3, 1:3, 1:3)
						local_sigsy(1:3, 1:3, 1:3) = perty(i-1:i+1, j-1:j+1, m-1:m+1) - perty(i, j, m) - local_sigcy(1:3, 1:3, 1:3)
						local_sigsz(1:3, 1:3, 1:3) = pertz(i-1:i+1, j-1:j+1, m-1:m+1) - pertz(i, j, m) - local_sigcz(1:3, 1:3, 1:3)

						if( divx .ge. 0.0_8 ) then
							sigcx(i,j,m) = sum( local_sigcx**2 + local_sigcy**2 + local_sigcz**2 ) / 27.0_8
						else
							sigcy(i,j,m) = sum( local_sigcx**2 + local_sigcy**2 + local_sigcz**2 ) / 27.0_8
						end if
						sigsx(i,j,m) = sum( local_sigsx**2 + local_sigsy**2 + local_sigsz**2 ) / 27.0_8
						sigcz(i,j,m) = 2.0_8 * sum( local_sigcx*local_sigsx + local_sigcy*local_sigsy + local_sigcz*local_sigsz ) / 27.0_8
					end do
				end do
			end do
		end if
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
		integer,intent(in) :: Gwidth, nx, smooth, window
		real(8),intent(inout) :: grid(nx)
		complex(8),allocatable :: fft_in1d(:), gauss_1d(:)
		real(8) :: Hwidth, Hwidth2, rx
		integer :: i, error, gw, fw, thx(2), nxhalf
		type(DFTI_DESCRIPTOR), POINTER :: fft_desc1

		print *, 'entering 1d smoothing'
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

		print *, 'defining fft_in1d matrix'
		allocate( fft_in1d(nx), stat=i )                      !!! deallocated at the end of this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation fft_in1d. stat= ', i
			stop
		end if
		print *, 'defining gaussian_1d matrix'
		allocate( gauss_1d(nx), stat=i )                      !!! deallocated at the end of this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation gauss_1d. stat= ', i
			stop
		end if

		rx = real(nx+1)/2.0_8
		fft_in1d(:) = grid(:)
		if( window .eq. 2 ) then
			Hwidth2 = (real(Gwidth)/2.0_8)**2
			do i=1,nx
				gauss_1d(i) = 0.5_8**( (real(rx-i)**2) / Hwidth2 ) 
			end do
		elseif( window .eq. 1 ) then
			Hwidth = real(Gwidth)/2.0_8
			gauss_1d(:) = 0.0_8
			thx(1) = ceiling( rx - Hwidth )
			thx(2) = ceiling( rx + Hwidth ) - 1
			do i=thx(1),thx(2)
				gauss_1d(i) = 1.0_8
			end do
		end if
		gauss_1d(:) = gauss_1d(:) / sum(gauss_1d(:))
!		print *, 'defined 1d vectors'

		error = DftiCreateDescriptor(fft_desc1, DFTI_DOUBLE, DFTI_COMPLEX, 1, nx )
		error = DftiCommitDescriptor(fft_desc1)
		error = DftiComputeForward(fft_desc1, fft_in1d)
		error = DftiFreeDescriptor(fft_desc1)

		error = DftiCreateDescriptor(fft_desc1, DFTI_DOUBLE, DFTI_COMPLEX, 1, nx )
		error = DftiCommitDescriptor(fft_desc1)
		error = DftiComputeForward(fft_desc1, gauss_1d)
		error = DftiFreeDescriptor(fft_desc1)

		fft_in1d(:) = fft_in1d(:)*gauss_1d(:)
!		print *, 'convolved'

		deallocate(gauss_1d)

		error = DftiCreateDescriptor(fft_desc1, DFTI_DOUBLE, DFTI_COMPLEX, 1, nx )
		error = DftiCommitDescriptor(fft_desc1)
		error = DftiComputeBackward(fft_desc1, fft_in1d)
		error = DftiFreeDescriptor(fft_desc1)

		print *, 'making output smoothed matrix'
		gw = 1-smooth
		fw = 1-2*smooth
		fft_in1d = fft_in1d / nx
		print *, 'max and min, real and imag parts of inverse FFT'
		print '(4(1x,es12.5))', maxval(real(fft_in1d)), minval(real(fft_in1d)), maxval(aimag(fft_in1d)), minval(aimag(fft_in1d))
		nxhalf = nx/2
		do i=1,nxhalf
			grid(i)        = gw*grid(i)        - fw*fft_in1d(nxhalf+i)
			grid(nxhalf+i) = gw*grid(nxhalf+i) - fw*fft_in1d(i)
		end do
		deallocate(fft_in1d)
		print *, 'done with 1d smoothing'
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
		integer,intent(in) :: Gwidth, nx, ny, smooth, window
		real(8),intent(inout) :: grid(nx, ny)
		complex(8),allocatable :: fft_in1d(:), gauss_1d(:)
		real(8) :: Hwidth, Hwidth2, rx, ry
		integer :: i, j, error, length(2), nxy, gw, fw, thx(2), thy(2), nxhalf, nyhalf
		type(DFTI_DESCRIPTOR), POINTER :: fft_desc2

		print *, 'entering 2d smoothing'
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

!		print *, 'defining fft_in1d matrix'
		allocate( fft_in1d(nxy), stat=i )                      !!! deallocated at the end of this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation fft_in1d. stat= ', i
			stop
		end if
!		print *, 'defining gaussian_1d matrix'
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
			Hwidth2 = (real(Gwidth)/2.0_8)**2
			do j=1,ny
				do i=1,nx
					gauss_1d( (j-1)*nx + i ) = 0.5_8**( (real(rx-i)**2 + real(ry-j)**2) / Hwidth2 ) 
				end do
			end do
		elseif( window .eq. 1 ) then
			Hwidth = real(Gwidth)/2.0_8
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
!		print *, 'defined 1d vectors'

		error = DftiCreateDescriptor(fft_desc2, DFTI_DOUBLE, DFTI_COMPLEX, 2, length )
		error = DftiCommitDescriptor(fft_desc2)
		error = DftiComputeForward(fft_desc2, fft_in1d)
		error = DftiFreeDescriptor(fft_desc2)

		error = DftiCreateDescriptor(fft_desc2, DFTI_DOUBLE, DFTI_COMPLEX, 2, length )
		error = DftiCommitDescriptor(fft_desc2)
		error = DftiComputeForward(fft_desc2, gauss_1d)
		error = DftiFreeDescriptor(fft_desc2)

		fft_in1d(:) = fft_in1d(:)*gauss_1d(:)
!		print *, 'convolved'

		deallocate(gauss_1d)

		error = DftiCreateDescriptor(fft_desc2, DFTI_DOUBLE, DFTI_COMPLEX, 2, length )
		error = DftiCommitDescriptor(fft_desc2)
		error = DftiComputeBackward(fft_desc2, fft_in1d)
		error = DftiFreeDescriptor(fft_desc2)

		gw = 1-smooth
		fw = 1-2*smooth
		fft_in1d = fft_in1d / nxy
		print *, 'max and min, real and imag parts of inverse FFT'
		print '(4(1x,es12.5))', maxval(real(fft_in1d)), minval(real(fft_in1d)), maxval(aimag(fft_in1d)), minval(aimag(fft_in1d))
		nxhalf = nx/2
		nyhalf = ny/2
		do j=1,nyhalf
			do i=1,nxhalf
				grid(i, j)               = gw*grid(i, j)               - fw*fft_in1d( (nyhalf+j-1)*nx + nxhalf+i )
				grid(nxhalf+i, nyhalf+j) = gw*grid(nxhalf+i, nyhalf+j) - fw*fft_in1d( (j-1)*nx + i )
				grid(i, nyhalf+j)        = gw*grid(i, nyhalf+j)        - fw*fft_in1d( (j-1)*nx + nxhalf+i )
				grid(nxhalf+i, j)        = gw*grid(nxhalf+i, j)        - fw*fft_in1d( (nyhalf+j-1)*nx + i )
			end do
		end do
		deallocate(fft_in1d)
		print *, 'done with 2d smoothing'
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
		integer,intent(in) :: Gwidth, nx, ny, nz, smooth, window
		real(8),intent(inout) :: grid(nx, ny, nz)
		complex(8),allocatable :: fft_in1d(:), gauss_1d(:)
		real(8) :: Hwidth, Hwidth2, rx, ry, rz
		integer :: i, j, k, error, length(3), nxy, nxyz, gw, fw, thx(2), thy(2), thz(2), nxhalf, nyhalf, nzhalf
		type(DFTI_DESCRIPTOR), POINTER :: fft_desc3

		print *, 'entering 3d smoothing'
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

		print *, 'defining fft_in1d matrix'
		allocate( fft_in1d(nxyz), stat=i )				!!! deallocated in this subroutine !!!
		if(i.ne.0) then
			print *, 'error in allocation fft_in1d. stat= ', i
			stop
		end if
		print *, 'defining gaussian_1d matrix'
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
			Hwidth2  = (real(Gwidth)/2.0_8)**2
			do k=1,nz
				do j=1,ny
					do i=1,nx
						gauss_1d( (k-1)*nxy + (j-1)*nx + i ) = 0.5_4**( (real(rx-i)**2 + real(ry-j)**2 + real(rz-k)**2)/Hwidth2 )
					end do
				end do
			end do
		elseif( window .eq. 1 ) then
			Hwidth = real(Gwidth)/2.0_8
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
!		print *, 'defined 1d vectors'

		error = DftiCreateDescriptor(fft_desc3, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc3)
		error = DftiComputeForward(fft_desc3, fft_in1d)
		error = DftiFreeDescriptor(fft_desc3)

		error = DftiCreateDescriptor(fft_desc3, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc3)
		error = DftiComputeForward(fft_desc3, gauss_1d)
		error = DftiFreeDescriptor(fft_desc3)

		fft_in1d(:) = fft_in1d(:)*gauss_1d(:)
!		print *, 'convolved'

		deallocate(gauss_1d)

		error = DftiCreateDescriptor(fft_desc3, DFTI_DOUBLE, DFTI_COMPLEX, 3, length )
		error = DftiCommitDescriptor(fft_desc3)
		error = DftiComputeBackward(fft_desc3, fft_in1d)
		error = DftiFreeDescriptor(fft_desc3)

		print *, 'making output smoothed matrix'
		gw = 1-smooth
		fw = 1-2*smooth
		fft_in1d = fft_in1d / nxyz
		print *, 'max and min, real and imag parts of inverse FFT'
		print '(4(1x,es12.5))', maxval(real(fft_in1d)), minval(real(fft_in1d)), maxval(aimag(fft_in1d)), minval(aimag(fft_in1d))
		nxhalf = nx/2
		nyhalf = ny/2
		nzhalf = nz/2
		do k=1,nzhalf
			do j=1,nyhalf
				do i=1,nxhalf
					grid(i, j, k)                      = gw*grid(i, j, k)                      - fw*fft_in1d( (nzhalf+k-1)*nxy + (nyhalf+j-1)*nx + nxhalf+i )
					grid(nxhalf+i, nyhalf+j, nzhalf+k) = gw*grid(nxhalf+i, nyhalf+j, nzhalf+k) - fw*fft_in1d( (k-1)*nxy + (j-1)*nx + i )
					grid(i, j, nzhalf+k)               = gw*grid(i, j, nzhalf+k)               - fw*fft_in1d( (k-1)*nxy + (nyhalf+j-1)*nx + nxhalf+i )
					grid(nxhalf+i, nyhalf+j, k)        = gw*grid(nxhalf+i, nyhalf+j, k)        - fw*fft_in1d( (nzhalf+k-1)*nxy + (j-1)*nx + i )
					grid(i, nyhalf+j, k)               = gw*grid(i, nyhalf+j, k)               - fw*fft_in1d( (nzhalf+k-1)*nxy + (j-1)*nx + nxhalf+i )
					grid(nxhalf+i, j, nzhalf+k)        = gw*grid(nxhalf+i, j, nzhalf+k)        - fw*fft_in1d( (k-1)*nxy + (nyhalf+j-1)*nx + i )
					grid(nxhalf+i, j, k)               = gw*grid(nxhalf+i, j, k)               - fw*fft_in1d( (nzhalf+k-1)*nxy + (nyhalf+j-1)*nx + i )
					grid(i, nyhalf+j, nzhalf+k)        = gw*grid(i, nyhalf+j, nzhalf+k)        - fw*fft_in1d( (k-1)*nxy + (j-1)*nx + nxhalf+i )
				end do
			end do
		end do
		deallocate(fft_in1d)
		print *, 'done with 3D smoothing'
	end subroutine fft_smooth_3D
end module decomp
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program main
use parameters
use generate
use decomp

	implicit none
	character(len=256) :: filename

	print *, 'reading parameters'
	call read_param()
	print *, 'done'
	call open_files()
	call make_perturbations()
	call turb_decomp()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contains
	subroutine read_param()
!!!!!!!!!! Read in parameter file and check for inconsistencies !!!!!!!!!!
		implicit none
		character(len=24) :: s
		integer :: j

		write(filename,'(a)') './parameter_input.dat'
		open(unit=14,file=filename)

		read(14,'(a24,F11.9)') s,pi
		print '(a24,F11.9)', s,pi
		read(14,'(a24,I4)') s,Nseed
		print '(a24,I4)', s,Nseed
		read(14,'(a24,I4)') s,kmin
		print '(a24,I4)', s,kmin
		read(14,'(a24,I4)') s,kmax
		print '(a24,I4)', s,kmax
		read(14,'(a24,I4)') s,Nx0
		print '(a24,I4)', s,Nx0
		read(14,'(a24,I4)') s,Ny0
		print '(a24,I4)', s,Ny0
		read(14,'(a24,I4)') s,Nz0
		print '(a24,I4)', s,Nz0
		read(14,'(a24,F4.2)') s,alpha
		print '(a24,F4.2)', s,alpha
		read(14,'(a24,F4.2)') s,f_solenoidal
		print '(a24,F4.2)', s,f_solenoidal
		read(14,'(a24,I1)') s,j
		print '(a24,I1)', s,j
		if(j .eq. 1) then
			use_sine = .true.
		elseif(j .eq. 0) then
			use_sine = .false.
		else
			print *, 'Problem with use_sine in parameter input file'
			print *, 'It should be equal to 0 or 1, but instead it is:'
			print *, j
			stop
		end if
		read(14,'(a24,F5.2)') s,Vrot
		print '(a24,F5.2)', s,Vrot
		read(14,'(a24,F5.2)') s,beta
		print '(a24,F5.2)', s,beta
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
		read(14,'(a24,I4)') s,pre_smooth_FWHM
		print '(a24,I4)', s,pre_smooth_FWHM
		read(14,'(a24,I1)') s,decomp_type
		print '(a24,I1)', s,decomp_type
		read(14,'(a24,I4)') s,type_3_FWHM
		print '(a24,I4)', s,type_3_FWHM
		read(14,'(a24,I1)') s,j
		print '(a24,I1)', s,j
		if(j .eq. 1) then
			fourier_decomp = .true.
		elseif(j .eq. 0) then
			fourier_decomp = .false.
		else
			print *, 'Problem with fourier_decomp in parameter input file'
			print *, 'It should be equal to 0 or 1, but instead it is:'
			print *, j
			stop
		end if
		read(14,'(a24,I1)') s,local_version
		print '(a24,I1)', s,local_version
		read(14,'(a24,F4.2)') s,sphererad
		print '(a24,F4.2)', s,sphererad
		read(14,'(a24,F4.2)') s,cylheigh
		print '(a24,F4.2)', s,cylheigh
		read(14,'(a24,I1)') s,j
		print '(a24,I1)', s,j
		if(j .eq. 1) then
			Ebox = .true.
		elseif(j .eq. 0) then
			Ebox = .false.
		else
			print *, 'Problem with Ebox in parameter input file'
			print *, 'It should be equal to 0 or 1, but instead it is:'
			print *, j
			stop
		end if

		!!! checks
		if( MOD(Nx0,2) .ne. 0 ) then
			print *, 'Nx0 must be even, but it is odd'
			print *, 'Nx0=',Nx0
			print *, 'adding 1 to Nx0'
			Nx0 = Nx0 + 1
		end if
		if( MOD(Ny0,2) .ne. 0 ) then
			print *, 'Ny0 must be even, but it is odd'
			print *, 'Ny0=',Ny0
			print *, 'adding 1 to Ny0'
			Ny0 = Ny0 + 1
		end if
		if( MOD(Nz0,2) .ne. 0 ) then
			print *, 'Nz0 must be even, but it is odd'
			print *, 'Nz0=',Nz0
			print *, 'adding 1 to Nz0'
			Nz0 = Nz0 + 1
		end if
		if(Ny0 .gt. Nx0) then
			j = Ny0
			Ny0 = Nx0
			Nx0 = j
		end if
		if(Nz0 .gt. Nx0) then
			j = Nz0
			Nz0 = Ny0
			Ny0 = Nx0
			Nx0 = j
		elseif(Nz0 .gt. Ny0) then
			j = Nz0
			Nz0 = Ny0
			Ny0 = j
		end if

		if( kmax .le. kmin ) then
			print *, 'ERROR'
			print *, 'must have Lmin<Lmax'
			print *, 'kmax, kmin, Nx0, Ny0, Nz0'
			print *, kmax, kmin, Nx0, Ny0, Nz0
			stop
		end if
		if( f_solenoidal .lt. 0.0_8 .or. f_solenoidal .gt. 1.0_8 ) then
			print *, 'NOTE'
			print *, 'must have 0<f_solenoidal<1 to set ratio'
			print *, f_solenoidal
			print *, 'fourier field will remain purely Gaussian'
		end if
		if( Vrot .le. 0.0_8 ) then
			print *, 'NOTE'
			print *, 'no rotation will be added, since Vrot is not positive'
			print *, Vrot
			Vrot = 0.0_8
		end if
		if(pre_smooth_dim .ne. 2 .and. pre_smooth_dim .ne. 3) then
			print *, 'Problem with pre_smooth_dim in parameter input file'
			print *, 'It should be equal to 2 or 3, but instead it is:'
			print *, pre_smooth_dim
			stop
		end if
		if( decomp_type .ne. 1 .and. decomp_type .ne. 2 .and. decomp_type .ne. 3 ) then
			print *, 'Problem with decomp_type in parameter input file'
			print *, 'It should be equal to 1 or 2 or 3, but instead it is:'
			print *, decomp_type
			stop
		end if
		if( local_version .ne. 1 .and. local_version .ne. 2 ) then
			print *, 'Problem with local_version in parameter input file'
			print *, 'It should be equal to 1 or 2, but instead it is:'
			print *, local_version
			stop
		end if
		sphererad = min(sphererad, 1.0_8)
		if( sphererad .le. 0.0_8 ) then
			print *, 'NOTE'
			print *, 'must have 0<sphererad to extract region'
			print *, sphererad
			print *, 'no sub-volume will be extracted'
		end if
		cylheigh = min(cylheigh, 1.0_8)
		if( cylheigh .le. 0.0_8 ) then
			print *, 'NOTE'
			print *, 'must have 0<cylheigh to extract cylinder'
			print *, cylheigh
			print *, 'no cylinder will be extracted'
		end if

	end subroutine read_param
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine open_files()
		implicit none

		write(output_dirname,'(a)') './outputs'
		write(filename,'(a,a)') 'mkdir -p ',trim(output_dirname)
		call system(trim(filename))

		write(output_dirname,'(a,a,i4.4,a,i4.4,a,i4.4,a,i4.4,a,i4.4)') trim(output_dirname),'/Nx0_',Nx0,'_Ny0_',Ny0,'_Nz0_',Nz0,'_kmin_',kmin,'_kmax_',kmax
		write(filename,'(a,a)') 'mkdir -p ',trim(output_dirname)
		call system(trim(filename))

		if( f_solenoidal .ge. 0.0_8 .and. f_solenoidal .le. 1.0_8 ) then
			write(output_dirname,'(a,a,f4.2,a,f4.2)') trim(output_dirname),'/fs_',f_solenoidal,'_alpha_',alpha
		else
			write(output_dirname,'(a,a,f4.2)') trim(output_dirname),'/fs_Gauss_alpha_',alpha
		end if

		if( vrot .gt. 0.0_8 ) then
			if( beta .ge. 0.0_8 ) then
				write(output_dirname,'(a,a,f4.2,a,f4.2)') trim(output_dirname),'_Vrot_',Vrot,'_beta_',beta
			else
				write(output_dirname,'(a,a,f4.2,a,f4.2)') trim(output_dirname),'_Vrot_',Vrot,'_beta_m',abs(beta)
			end if
		else
			write(output_dirname,'(a,a)') trim(output_dirname),'_NoRot'
		end if
		write(filename,'(a,a)') 'mkdir -p ',trim(output_dirname)
		call system(trim(filename))

		if( fourier_decomp ) then
			write(output_dirname,'(a,a,i1)') trim(output_dirname),'/fourier_decomp_',decomp_type
		else
			write(output_dirname,'(a,a,i1,a,i1)') trim(output_dirname),'/local_v',local_version,'_decomp_',decomp_type
		end if
		if( decomp_type .eq. 3 ) then
			write(output_dirname,'(a,a,i4.4)') trim(output_dirname),'_',type_3_FWHM
		end if
		write(filename,'(a,a)') 'mkdir -p ',trim(output_dirname)
		call system(trim(filename))

		if( pre_smooth ) then
			write(output_dirname,'(a,a,i1,a,i4.4)') trim(output_dirname),'/',pre_smooth_dim,'D_smooth_',pre_smooth_FWHM
		else
			write(output_dirname,'(a,a)') trim(output_dirname),'/no_pre_smooth'
		end if
		write(filename,'(a,a)') 'mkdir -p ',trim(output_dirname)
		call system(trim(filename))

		if( Ebox ) then
			write(output_dirname,'(a,a,f4.2,a,f4.2,a)') trim(output_dirname),'/R_',sphererad,'_H_',cylheigh,'Ebox'
		else
			write(output_dirname,'(a,a,f4.2,a,f4.2)') trim(output_dirname),'/R_',sphererad,'_H_',cylheigh
		end if
		write(filename,'(a,a)') 'mkdir -p ',trim(output_dirname)
		call system(trim(filename))

	end subroutine open_files
end program main

