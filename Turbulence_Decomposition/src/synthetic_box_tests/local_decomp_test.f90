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
	integer :: type_2_FWHM			!!! FWHM of Gaussian used to smooth the rotation curve before subtracting it (decomp_type=2)
	integer :: Nvert			!!! Number of vertical bins to use when subtracting rotation curve (decomp_type=2)
	integer :: type_3_FWHM			!!! FWHM of Gaussian used for defining field to be decomposed if decomp_type=3, in number of cells
	logical :: fourier_decomp		!!! perform Helmholtz decomposition in Fourier space. Otherwise use strain rate tensor in real space
	integer :: local_version		!!! 1 for diagonal of strain rate tensor, 2 for full divergence times {\vec {r}}

	!!! Variables for volume of region within which turbulence is defined and decomposed
	real(8) :: sphererad	! If >0, then perturbations are set to zero outside spherical region, and the perturbation field is shifted and renormalized to keep the center of mass velocity at zero and the variance at unity; the spherical region cut out is centered at the center of the perturbation cube, and has a radius given by the value of this parameter, with sphererad = 1 corresponding to the spherical region going all the way to the edge of the perturbation cube
	real(8) :: cylheigh	! If >0, then perturbations are set to zero outside cylindrical region, and the perturbation field is shifted and renormalized to keep the center of mass velocity at zero and the variance at unity; the cylindrical region cut out is centered at the center of the perturbation cube, and has a height given by the value of this parameter, with cylheigh = 1 corresponding to the height of the region going all the way to the edge of the perturbation cube. The radius of the cylinder is given by sphererad
	logical :: Ebox		! If extracting sphere or cylinder, rather than set cells outside geometric region within the large grid to 0, create smaller grid which won't have any 0s.
	integer :: rebin	! For local decomposition, rebin data so that each new cell is 'rebin>1' old cells.

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
		complex(8),allocatable :: fxs1d(:), fys1d(:), fzs1d(:), fxc1d(:), fyc1d(:), fzc1d(:)
		real(8) :: norms, normc, rescale
		real(8),allocatable :: pertsx(:,:,:), pertsy(:,:,:), pertsz(:,:,:), pertcx(:,:,:), pertcy(:,:,:), pertcz(:,:,:)
		real(8),allocatable :: test(:,:,:), Skmag(:,:,:)
		real(8) :: rad
		logical,allocatable :: kmask(:,:,:)
		integer :: i, j, m, n1, length(3), error
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

		print *, '!!! TEST BRUNT AND FEDERRATH !!!'
		allocate( test(Nx0,Ny0,Nz0) )
		rad = real(Nx0)*real(Ny0)*real(Nz0)
		print *, 'x comp: max, min, mean'
		test = abs(fxc)**2 - ( abs(fxc)**2 + abs(fyc)**2 + abs(fzc)**2 ) * ( kx**2 / kmag**2 )
		print '(3(1x,es12.5))', maxval(test, mask=kmag.ne.0.0_8), minval(test, mask=kmag.ne.0.0_8), sum(test, mask=kmag.ne.0.0_8) / rad
		print *, 'y comp: max, min, mean'
		test = abs(fyc)**2 - ( abs(fxc)**2 + abs(fyc)**2 + abs(fzc)**2 ) * ( ky**2 / kmag**2 )
		print '(3(1x,es12.5))', maxval(test, mask=kmag.ne.0.0_8), minval(test, mask=kmag.ne.0.0_8), sum(test, mask=kmag.ne.0.0_8) / rad
		print *, 'z comp: max, min, mean'
		test = abs(fzc)**2 - ( abs(fxc)**2 + abs(fyc)**2 + abs(fzc)**2 ) * ( kz**2 / kmag**2 )
		print '(3(1x,es12.5))', maxval(test, mask=kmag.ne.0.0_8), minval(test, mask=kmag.ne.0.0_8), sum(test, mask=kmag.ne.0.0_8) / rad
		print *, 'x sol: max, min, mean'
		test = abs(fxs)**2 - ( abs(fxs)**2 + abs(fys)**2 + abs(fzs)**2 ) * (ky**2+kz**2) / (2.0_8*kmag**2)
		print '(3(1x,es12.5))', maxval(test, mask=kmag.ne.0.0_8), minval(test, mask=kmag.ne.0.0_8), sum(test, mask=kmag.ne.0.0_8) / rad
		print *, 'y sol: max, min, mean'
		test = abs(fys)**2 - ( abs(fxs)**2 + abs(fys)**2 + abs(fzs)**2 ) * (kx**2+kz**2) / (2.0_8*kmag**2)
		print '(3(1x,es12.5))', maxval(test, mask=kmag.ne.0.0_8), minval(test, mask=kmag.ne.0.0_8), sum(test, mask=kmag.ne.0.0_8) / rad
		print *, 'z sol: max, min, mean'
		test = abs(fzs)**2 - ( abs(fxs)**2 + abs(fys)**2 + abs(fzs)**2 ) * (kx**2+ky**2) / (2.0_8*kmag**2)
		print '(3(1x,es12.5))', maxval(test, mask=kmag.ne.0.0_8), minval(test, mask=kmag.ne.0.0_8), sum(test, mask=kmag.ne.0.0_8) / rad
		print *, 'z sol z=0 plane: max, min, mean'
		rad = real(Nx0)*real(Ny0)
		print '(3(1x,es12.5))', maxval(test, mask=kmag.ne.0.0_8.and.abs(kz).lt.1.d-9), minval(test, mask=kmag.ne.0.0_8.and.abs(kz).lt.1.d-9), sum(test, mask=kmag.ne.0.0_8.and.abs(kz).lt.1.d-9) / rad
		deallocate(test)
		print *, '!!! DONE !!!'

		deallocate( kx, ky, kz, kmag )
		if( use_sine ) then
			deallocate( Skmag )
		end if

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
!					rad2 = abs(real(m)-rz) / (0.5_8 * real(Nzf))
!					pertx(i,j,m) = pertx(i,j,m) + Vrot * (rad**(beta-1)) * real(ry-j) / ( cosh(rad2) )**2
!					perty(i,j,m) = perty(i,j,m) + Vrot * (rad**(beta-1)) * real(i-rx) / ( cosh(rad2) )**2
					pertx(i,j,m) = pertx(i,j,m) + Vrot * (rad**(beta-1)) * real(ry-j)
					perty(i,j,m) = perty(i,j,m) + Vrot * (rad**(beta-1)) * real(i-rx)
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
	real(8),allocatable :: sigcxx(:,:,:,:), sigcyy(:,:,:,:), sigczz(:,:,:,:), sigsxx(:,:,:,:), sigsyy(:,:,:,:), sigszz(:,:,:,:)
	real(8) :: sigma_ratio(8)

contains
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine turb_def()
		implicit none
		integer :: i, j, k, m, n, n1, Nrot
		real(8) :: rad, rx, ry, vrot_interp, rm, vxt, vyt, Hw, norm
		real(8),allocatable :: Vx(:,:), Vy(:,:), vrot_curve(:,:), rrot(:), test(:,:)
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

			if( type_2_FWHM .gt. 0.0_8 ) then
				Hw = 0.5_8 * type_2_FWHM
				allocate( test(Nrot,Nvert), stat=i )
				if(i.ne.0) then
					print *, 'error in allocation of test. stat= ', i
					stop
				end if
				print *, 'smoothing rotation curve'
				test(:,:) = vrot_curve(:,:)
				do j=1,Nrot
					norm = sum( 0.5_8**((abs(rrot(:)-rrot(j))/Hw)**2), mask = abs(rrot(:)-rrot(j)) .le. 5.0_8*Hw )
					do i=1,Nvert
						vrot_curve(j,i) = sum( test(:,i) * 0.5_8**((abs(rrot(:)-rrot(j))/Hw)**2), mask = abs(rrot(:)-rrot(j)) .le. 5.0_8*Hw ) / norm
					end do
				end do
				deallocate( test )
			end if

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
!								vrot_interp = log10( vrot_curve(m-1,n) ) + log10( rad/rrot(m-1) ) * log10( vrot_curve(m,n)/vrot_curve(m-1,n) ) / &
!									&     log10( rrot(m)/rrot(m-1) )
								vrot_interp = vrot_curve(m-1,n) + ( rad - rrot(m-1) ) * ( vrot_curve(m,n) - vrot_curve(m-1,n) ) / ( rrot(m) - rrot(m-1) )
							else
!								vrot_interp = log10( vrot_curve(m,n) ) + log10( rad/rrot(m) ) * log10( vrot_curve(m+1,n)/vrot_curve(m,n) ) / &
!									&     log10( rrot(m+1)/rrot(m) )
								vrot_interp = vrot_curve(m,n) + ( rad - rrot(m) ) * ( vrot_curve(m+1,n) - vrot_curve(m,n) ) / ( rrot(m+1) - rrot(m) )
							end if
!							vrot_interp = 10**(vrot_interp)
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

	end subroutine turb_def
!_____________________________________________________________________________________________________________________________________________________________________

	subroutine turb_decomp()
	! Final part
	! decomposes turbulence into compreesive and solenoidal components, calculates the global ratios, and creates projected 2D plots of them
	! snapshot is snapshot number from main loop, needed for io purposes
		implicit none
		integer :: i, j, k, m
		real(8) :: Ei, fsol, rad
		real(8),allocatable :: test(:,:,:)
		character(len=256) :: filename

		! Define the field that will be decomposed
		call turb_def()
		print *, 'Starting decomposition'
		Ei = sum(pertx**2 + perty**2 + pertz**2) / (Nxf*Nyf*Nzf)
		print '(a3,1x,es12.5)', 'RMS', sqrt(Ei)

		call local_decompose
		deallocate( pertx, perty, pertz )
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		print *, 'total power real space'
		Ei = sum(sigcx + sigcy + sigcz + sigsx + sigsy + sigsz) + sum(sigcxx + sigcyy + sigczz) + sum( sigsxx(:,:,:,1:6) + sigsyy(:,:,:,1:6) + sigszz(:,:,:,1:6) )
		if( f_solenoidal .ge. 0.0_8 .and. f_solenoidal .le. 1.0_8 ) then
			fsol = (2.0_8*f_solenoidal**2) / ( 2.0_8*f_solenoidal**2 + (1.0_8-f_solenoidal)**2 )
		else
			fsol = 2.0_8 / 3.0_8
		end if
		write(20,'(46(1x,es12.5))') fsol, sum(sigcx)/Ei, sum(sigcy)/Ei, sum(sigcz)/Ei, sum(sigsx)/Ei, sum(sigsy)/Ei, sum(sigsz)/Ei, sum(sigcxx(:,:,:,1))/Ei, & 
			& sum(sigcyy(:,:,:,1))/Ei, sum(sigczz(:,:,:,1))/Ei, sum(sigsxx(:,:,:,1))/Ei, sum(sigsyy(:,:,:,1))/Ei, sum(sigszz(:,:,:,1))/Ei, sum(sigcxx(:,:,:,2))/Ei, & 
			& sum(sigcyy(:,:,:,2))/Ei, sum(sigczz(:,:,:,2))/Ei, sum(sigsxx(:,:,:,2))/Ei, sum(sigsyy(:,:,:,2))/Ei, sum(sigszz(:,:,:,2))/Ei, sum(sigcxx(:,:,:,3))/Ei, & 
			& sum(sigcyy(:,:,:,3))/Ei, sum(sigczz(:,:,:,3))/Ei, sum(sigsxx(:,:,:,3))/Ei, sum(sigsyy(:,:,:,3))/Ei, sum(sigszz(:,:,:,3))/Ei, sum(sigcxx(:,:,:,4))/Ei, & 
			& sum(sigcyy(:,:,:,4))/Ei, sum(sigczz(:,:,:,4))/Ei, sum(sigsxx(:,:,:,4))/Ei, sum(sigsyy(:,:,:,4))/Ei, sum(sigszz(:,:,:,4))/Ei, sum(sigcxx(:,:,:,5))/Ei, & 
			& sum(sigcyy(:,:,:,5))/Ei, sum(sigczz(:,:,:,5))/Ei, sum(sigsxx(:,:,:,5))/Ei, sum(sigsyy(:,:,:,5))/Ei, sum(sigszz(:,:,:,5))/Ei, sum(sigcxx(:,:,:,6))/Ei, & 
			& sum(sigcyy(:,:,:,6))/Ei, sum(sigczz(:,:,:,6))/Ei, sum(sigsxx(:,:,:,6))/Ei, sum(sigsyy(:,:,:,6))/Ei, sum(sigszz(:,:,:,6))/Ei, sum(sigcxx(:,:,:,7))/Ei, & 
			& sum(sigcyy(:,:,:,7))/Ei, sum(sigczz(:,:,:,7))/Ei
		print *, 'ii*ii'
		print '(3(1x,es12.5))', sum(sigcx)/Ei, sum(sigcy)/Ei, sum(sigcz)/Ei
		print *, 'ii*jj'
		print '(3(1x,es12.5))', sum(sigsx)/Ei, sum(sigsy)/Ei, sum(sigsz)/Ei
		print *, 'ij*ij'
		print '(6(1x,es12.5))', sum(sigcxx(:,:,:,1))/Ei,sum(sigcyy(:,:,:,1))/Ei,sum(sigczz(:,:,:,1))/Ei,sum(sigsxx(:,:,:,1))/Ei,sum(sigsyy(:,:,:,1))/Ei,sum(sigszz(:,:,:,1))/Ei
		print *, 'ij*ji'
		print '(3(1x,es12.5))', sum(sigcxx(:,:,:,2))/Ei,sum(sigcyy(:,:,:,2))/Ei,sum(sigczz(:,:,:,2))/Ei
		print *, 'xx*ij'
		print '(6(1x,es12.5))', sum(sigsxx(:,:,:,2))/Ei,sum(sigsyy(:,:,:,2))/Ei,sum(sigszz(:,:,:,2))/Ei,sum(sigcxx(:,:,:,3))/Ei,sum(sigcyy(:,:,:,3))/Ei,sum(sigczz(:,:,:,3))/Ei
		print *, 'yy*ij'
		print '(6(1x,es12.5))', sum(sigsxx(:,:,:,3))/Ei,sum(sigsyy(:,:,:,3))/Ei,sum(sigszz(:,:,:,3))/Ei,sum(sigcxx(:,:,:,4))/Ei,sum(sigcyy(:,:,:,4))/Ei,sum(sigczz(:,:,:,4))/Ei
		print *, 'zz*ij'
		print '(6(1x,es12.5))', sum(sigsxx(:,:,:,4))/Ei,sum(sigsyy(:,:,:,4))/Ei,sum(sigszz(:,:,:,4))/Ei,sum(sigcxx(:,:,:,5))/Ei,sum(sigcyy(:,:,:,5))/Ei,sum(sigczz(:,:,:,5))/Ei
		print *, 'xy*{xz,yz,zx,zy}; yz*zx; zx*zy'
		print '(6(1x,es12.5))', sum(sigsxx(:,:,:,5))/Ei,sum(sigsyy(:,:,:,5))/Ei,sum(sigszz(:,:,:,5))/Ei,sum(sigcxx(:,:,:,6))/Ei,sum(sigcyy(:,:,:,7))/Ei,sum(sigczz(:,:,:,7))/Ei
		print *, 'xz*{yx,yz,zy}; yx*{yz,zx,zy}'
		print '(6(1x,es12.5))', sum(sigcyy(:,:,:,6))/Ei,sum(sigczz(:,:,:,6))/Ei,sum(sigsxx(:,:,:,6))/Ei,sum(sigsyy(:,:,:,6))/Ei,sum(sigszz(:,:,:,6))/Ei,sum(sigcxx(:,:,:,7))/Ei

		Ei = sum( sigcx+sigcy+sigcz + sigcxx(:,:,:,1)+sigcyy(:,:,:,1)+sigczz(:,:,:,1) + sigsxx(:,:,:,1)+sigsyy(:,:,:,1)+sigszz(:,:,:,1) + 2.0_8*(sigsx+sigsy+sigsz - sigcxx(:,:,:,2)-sigcyy(:,:,:,2)-sigczz(:,:,:,2)) )
		print '(3(1x,es12.5))', sum( sigcx + sigcy + sigcz + 2.0_8*(sigsx + sigsy + sigsz) )/Ei, &
		& sum( sigcxx(:,:,:,1)+sigcyy(:,:,:,1)+sigczz(:,:,:,1)+sigsxx(:,:,:,1)+sigsyy(:,:,:,1)+sigszz(:,:,:,1) - 2.0_8*(sigcxx(:,:,:,2)+sigcyy(:,:,:,2)+sigczz(:,:,:,2)) )/Ei, &
		& fsol

		deallocate( sigcx, sigcy, sigcz, sigsx, sigsy, sigsz )
		deallocate( sigcxx, sigcyy, sigczz, sigsxx, sigsyy, sigszz )

	end subroutine turb_decomp

!_____________________________________________________________________________________________________________________________________________________________________

	subroutine local_decompose
	! decomposes turbulence into solenoidal and compressive modes locally in real space (using strain rate tensor)

		implicit none
!		real(8) :: local_sigcx(9,9,9),local_sigcy(9,9,9),local_sigcz(9,9,9),local_sigsx(9,9,9),local_sigsy(9,9,9),local_sigsz(9,9,9)
!		real(8) :: a, b, c, d
		real(8) :: local_sigcx(3,3,3),local_sigcy(3,3,3),local_sigcz(3,3,3),local_sigsx(3,3,3),local_sigsy(3,3,3),local_sigsz(3,3,3)
		real(8) :: divxx, divxy, divxz, divyx, divyy, divyz, divzx, divzy, divzz
		real(8),allocatable :: tempx(:,:,:), tempy(:,:,:), tempz(:,:,:)
		integer :: i, j, m, ii, jj, mm, ni, nj, nm, n1i, n1j, n1m

		print *, 'Decomposing using strain rate tensor'
		if( rebin .gt. 1 .and. rebin .le. Nxf .and. rebin .le. Nyf .and. rebin .le. Nzf ) then
			print *, 'rebinning data'
			print *, 'rebin=',rebin
			ii = floor( real(Nxf)/real(rebin) )
			jj = floor( real(Nyf)/real(rebin) )
			mm = floor( real(Nzf)/real(rebin) )
			allocate( tempx(ii,jj,mm), tempy(ii,jj,mm), tempz(ii,jj,mm) )
			do m=1,mm
				nm  = rebin*(m-1) + 1
				n1m = rebin*m
				do j=1,jj
					nj  = rebin*(j-1) + 1
					n1j = rebin*j
					do i=1,ii
						ni  = rebin*(i-1) + 1
						n1i = rebin*i
						tempx(i,j,m) = sum( pertx(ni:n1i, nj:n1j, nm:n1m) ) / (real(rebin)**3)
						tempy(i,j,m) = sum( perty(ni:n1i, nj:n1j, nm:n1m) ) / (real(rebin)**3)
						tempz(i,j,m) = sum( pertz(ni:n1i, nj:n1j, nm:n1m) ) / (real(rebin)**3)
					end do
				end do
			end do
			deallocate( pertx, perty, pertz )
			Nxf = ii
			Nyf = jj
			Nzf = mm
			allocate( pertx(Nxf,Nyf,Nzf), perty(Nxf,Nyf,Nzf), pertz(Nxf,Nyf,Nzf) )
			pertx = tempx
			perty = tempy
			pertz = tempz
			deallocate( tempx, tempy, tempz )
			print *, 'done'
		end if
		allocate( sigcx(Nxf, Nyf, Nzf), sigcy(Nxf, Nyf, Nzf), sigcz(Nxf, Nyf, Nzf), sigsx(Nxf, Nyf, Nzf), sigsy(Nxf, Nyf, Nzf), sigsz(Nxf, Nyf, Nzf) )
		allocate( sigcxx(Nxf, Nyf, Nzf, 7), sigcyy(Nxf, Nyf, Nzf, 7), sigczz(Nxf, Nyf, Nzf, 7), sigsxx(Nxf, Nyf, Nzf, 6), sigsyy(Nxf, Nyf, Nzf, 6), sigszz(Nxf, Nyf, Nzf, 6) )
		sigcx = 0.0_8
		sigcy = 0.0_8
		sigcz = 0.0_8
		sigsx = 0.0_8
		sigsy = 0.0_8
		sigsz = 0.0_8

		sigcxx = 0.0_8
		sigcyy = 0.0_8
		sigczz = 0.0_8
		sigsxx = 0.0_8
		sigsyy = 0.0_8
		sigszz = 0.0_8

		local_sigcx(:,:,:) = 0.0_8
		local_sigcy(:,:,:) = 0.0_8
		local_sigcz(:,:,:) = 0.0_8
		local_sigsx(:,:,:) = 0.0_8
		local_sigsy(:,:,:) = 0.0_8
		local_sigsz(:,:,:) = 0.0_8
		do m=2,Nzf-1
			do j=2,Nyf-1
				do i=2,Nxf-1
					divxx = (pertx(i+1,j,m) - pertx(i-1,j,m)) / 2.0_8
					divxy = (pertx(i,j+1,m) - pertx(i,j-1,m)) / 2.0_8
					divxz = (pertx(i,j,m+1) - pertx(i,j,m-1)) / 2.0_8
					divyx = (perty(i+1,j,m) - perty(i-1,j,m)) / 2.0_8
					divyy = (perty(i,j+1,m) - perty(i,j-1,m)) / 2.0_8
					divyz = (perty(i,j,m+1) - perty(i,j,m-1)) / 2.0_8
					divzx = (pertz(i+1,j,m) - pertz(i-1,j,m)) / 2.0_8
					divzy = (pertz(i,j+1,m) - pertz(i,j-1,m)) / 2.0_8
					divzz = (pertz(i,j,m+1) - pertz(i,j,m-1)) / 2.0_8

					sigcx(i,j,m) = divxx**2
					sigcy(i,j,m) = divyy**2
					sigcz(i,j,m) = divzz**2
					sigsx(i,j,m) = divxx*divyy
					sigsy(i,j,m) = divxx*divzz
					sigsz(i,j,m) = divyy*divzz

					sigcxx(i,j,m,1) = divxy**2
					sigcyy(i,j,m,1) = divyx**2
					sigczz(i,j,m,1) = divxz**2
					sigsxx(i,j,m,1) = divzx**2
					sigsyy(i,j,m,1) = divyz**2
					sigszz(i,j,m,1) = divzy**2
					sigcxx(i,j,m,2) = divxy*divyx
					sigcyy(i,j,m,2) = divxz*divzx
					sigczz(i,j,m,2) = divyz*divzy

					sigsxx(i,j,m,2) = divxx*divxy
					sigsyy(i,j,m,2) = divxx*divxz
					sigszz(i,j,m,2) = divxx*divyx
					sigcxx(i,j,m,3) = divxx*divyz
					sigcyy(i,j,m,3) = divxx*divzx
					sigczz(i,j,m,3) = divxx*divzy

					sigsxx(i,j,m,3) = divyy*divxy
					sigsyy(i,j,m,3) = divyy*divxz
					sigszz(i,j,m,3) = divyy*divyx
					sigcxx(i,j,m,4) = divyy*divyz
					sigcyy(i,j,m,4) = divyy*divzx
					sigczz(i,j,m,4) = divyy*divzy
					sigsxx(i,j,m,4) = divzz*divxy
					sigsyy(i,j,m,4) = divzz*divxz
					sigszz(i,j,m,4) = divzz*divyx
					sigcxx(i,j,m,5) = divzz*divyz
					sigcyy(i,j,m,5) = divzz*divzx
					sigczz(i,j,m,5) = divzz*divzy

					sigsxx(i,j,m,5) = divxy*divxz
					sigsyy(i,j,m,5) = divxy*divyz
					sigszz(i,j,m,5) = divxy*divzx
					sigcxx(i,j,m,6) = divxy*divzy

					sigcyy(i,j,m,6) = divxz*divyx
					sigczz(i,j,m,6) = divxz*divyz
					sigsxx(i,j,m,6) = divxz*divzy

					sigsyy(i,j,m,6) = divyx*divyz
					sigszz(i,j,m,6) = divyx*divzx
					sigcxx(i,j,m,7) = divyx*divzy

					sigcyy(i,j,m,7) = divyz*divzx

					sigczz(i,j,m,7) = divzx*divzy
				end do
			end do
		end do

!			a = 1.0_8 / 280.0_8
!			b = 4.0_8 / 105.0_8
!			c = 1.0_8 / 5.0_8
!			d = 4.0_8 / 5.0_8
!			do m=5,Nzf-4
!				do j=5,Nyf-4
!					do i=5,Nxf-4
!						divxx =   a*pertx(i-4,j,m) - b*pertx(i-3,j,m) + c*pertx(i-2,j,m) - d*pertx(i-1,j,m) + & 
!							& d*pertx(i+1,j,m) - c*pertx(i+2,j,m) + b*pertx(i+3,j,m) - a*pertx(i+4,j,m)
!						divxy =   a*pertx(i,j-4,m) - b*pertx(i,j-3,m) + c*pertx(i,j-2,m) - d*pertx(i,j-1,m) + & 
!							& d*pertx(i,j+1,m) - c*pertx(i,j+2,m) + b*pertx(i,j+3,m) - a*pertx(i,j+4,m)
!						divxz =   a*pertx(i,j,m-4) - b*pertx(i,j,m-3) + c*pertx(i,j,m-2) - d*pertx(i,j,m-1) + & 
!							& d*pertx(i,j,m+1) - c*pertx(i,j,m+2) + b*pertx(i,j,m+3) - a*pertx(i,j,m+4)
!
!						divyx =   a*perty(i-4,j,m) - b*perty(i-3,j,m) + c*perty(i-2,j,m) - d*perty(i-1,j,m) + & 
!							& d*perty(i+1,j,m) - c*perty(i+2,j,m) + b*perty(i+3,j,m) - a*perty(i+4,j,m)
!						divyy =   a*perty(i,j-4,m) - b*perty(i,j-3,m) + c*perty(i,j-2,m) - d*perty(i,j-1,m) + & 
!							& d*perty(i,j+1,m) - c*perty(i,j+2,m) + b*perty(i,j+3,m) - a*perty(i,j+4,m)
!						divyz =   a*perty(i,j,m-4) - b*perty(i,j,m-3) + c*perty(i,j,m-2) - d*perty(i,j,m-1) + & 
!							& d*perty(i,j,m+1) - c*perty(i,j,m+2) + b*perty(i,j,m+3) - a*perty(i,j,m+4)
!
!						divzx =   a*pertz(i-4,j,m) - b*pertz(i-3,j,m) + c*pertz(i-2,j,m) - d*pertz(i-1,j,m) + & 
!							& d*pertz(i+1,j,m) - c*pertz(i+2,j,m) + b*pertz(i+3,j,m) - a*pertz(i+4,j,m)
!						divzy =   a*pertz(i,j-4,m) - b*pertz(i,j-3,m) + c*pertz(i,j-2,m) - d*pertz(i,j-1,m) + & 
!							& d*pertz(i,j+1,m) - c*pertz(i,j+2,m) + b*pertz(i,j+3,m) - a*pertz(i,j+4,m)
!						divzz =   a*pertz(i,j,m-4) - b*pertz(i,j,m-3) + c*pertz(i,j,m-2) - d*pertz(i,j,m-1) + & 
!							& d*pertz(i,j,m+1) - c*pertz(i,j,m+2) + b*pertz(i,j,m+3) - a*pertz(i,j,m+4)
!
!						do mm=1,9
!							do jj=1,9
!								do ii=1,9
!									local_sigcx(ii, jj, mm) = real(ii-5) * divxx
!									local_sigcy(ii, jj, mm) = real(jj-5) * divyy
!									local_sigcz(ii, jj, mm) = real(mm-5) * divzz
!
!									local_sigsx(ii, jj, mm) = real(jj-5) * divxy + real(mm-5) * divxz
!									local_sigsy(ii, jj, mm) = real(ii-5) * divyx + real(mm-5) * divyz
!									local_sigsz(ii, jj, mm) = real(ii-5) * divzx + real(jj-5) * divzy
!								end do
!							end do
!						end do
!!!						local_sigsx(1:9, 1:9, 1:9) = pertx(i-4:i+4, j-4:j+4, m-4:m+4) - pertx(i, j, m) - local_sigcx(1:9, 1:9, 1:9)
!!!						local_sigsy(1:9, 1:9, 1:9) = perty(i-4:i+4, j-4:j+4, m-4:m+4) - perty(i, j, m) - local_sigcy(1:9, 1:9, 1:9)
!!!						local_sigsz(1:9, 1:9, 1:9) = pertz(i-4:i+4, j-4:j+4, m-4:m+4) - pertz(i, j, m) - local_sigcz(1:9, 1:9, 1:9)
!
!						if( divxx + divyy + divzz .ge. 0.0_8 ) then
!							sigcx(i,j,m) = sum( local_sigcx**2 + local_sigcy**2 + local_sigcz**2 ) / 729.0_8
!						else
!							sigcy(i,j,m) = sum( local_sigcx**2 + local_sigcy**2 + local_sigcz**2 ) / 729.0_8
!						end if
!						sigsx(i,j,m) = sum( local_sigsx**2 + local_sigsy**2 + local_sigsz**2 ) / 729.0_8
!						sigcz(i,j,m) = 2.0_8 * sum( local_sigcx*local_sigsx + local_sigcy*local_sigsy + local_sigcz*local_sigsz ) / 729.0_8
!					end do
!				end do
!			end do

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
	integer :: i

	print *, 'reading parameters'
	call read_param()
	print *, 'done'
	call open_files()
	write(filename,'(a,a)') trim(output_dirname),'/elements.out'
	open(unit=20,file=filename)
	do i=0,100
		f_solenoidal = real(i) * 0.01_8
		call make_perturbations()
		call turb_decomp()
	end do
	close(unit=20)

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
		read(14,'(a24,I4)') s,type_2_FWHM
		print '(a24,I4)', s,type_2_FWHM
		read(14,'(a24,I4)') s,Nvert
		print '(a24,I4)', s,Nvert
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
		read(14,'(a24,I3)') s,rebin
		print '(a24,I3)', s,rebin

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

		write(output_dirname,'(a,a,i4.4,i4.4,i4.4,a,i4.4,i4.4)') trim(output_dirname),'/Nxyz0_',Nx0,Ny0,Nz0,'_kmin_max',kmin,kmax
		write(filename,'(a,a)') 'mkdir -p ',trim(output_dirname)
		call system(trim(filename))

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

		write(output_dirname,'(a,a,i1)') trim(output_dirname),'/decomp_',decomp_type
		if( decomp_type .eq. 2 ) then
			write(output_dirname,'(a,a,i4.4)') trim(output_dirname),'_',type_2_FWHM
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

