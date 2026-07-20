!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
use amr_parameters
use hydro_parameters
	implicit none
	integer ::nn                            ! Number of cells
	real(dp)::dx                            ! Cell size
	real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
	real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
	integer  :: id, iu, iv, ip, ic
	real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
	real(dp),dimension(1:nvector,1:MAXPERTURB),save::Amp_up_tot,Amp_down_tot
	integer  :: i, j, ivar
	real(dp) :: dx_levelmax, sigma, gauss_sigma, ramp, ramp_up, ramp_down
	real(dp) :: lambda, kx, Amp_down_over_Amp_up, xx, zz, zz_up, zz_down, non_eigen_pert
	real(dp) :: rad, z_up, z_down, V_rel, D_b, D_s, P_b, P_s, Vx_b, Vx_s, Vz_b, Vz_s, C_b, C_s, var_b, var_s	!!! Unperturbed variables
	complex(dp) :: imag=(0,1), omega, Q_b, Q_s, expon1, expon2
	complex(dp) :: P_pert_up, P_pert_mid, P_pert_down, P_pert, D_pert_up, D_pert_mid, D_pert_down, D_pert, Vx_pert_up, Vx_pert_mid, Vx_pert_down, Vx_pert, Vz_pert_up, Vz_pert_mid, Vz_pert_down, Vz_pert, Amp_up

!	omega = (20.2371, 19.1188)			!!! Complex frequency of the eigenmode --> calculated numerically GIVEN the wavelength ( lambda_over_rad )

	id=1; iu=2; iv=3; ip=4; ic=5			!!! density, long vel, trans vel, pres, color
	!!! Unperturbed initialization !!!
	! Make sure slab is centerd around z=0.5, and only use relative velocity
	rad    = 0.5 * abs(x_center(2) - x_center(1))	!!! Radius of unperturbed slab
	z_up   = 0.5 + rad				!!! upper edge of unperturbed CENTERED slab
	z_down = 0.5 - rad				!!! lower edge of unperturbed CENTERED slab

	D_b    = d_region(1)				!!! unperturbed background density
	D_s    = d_region(2)				!!! unperturbed slab density
	P_b    = p_region(1)				!!! unperturbed background pressure
	P_s    = p_region(2)				!!! unperturbed slab pressure
	C_b    = sqrt(gamma * P_b / D_b)		!!! Speed of sound in background
	C_s    = sqrt(gamma * P_s / D_s)		!!! Speed of sound in slab

	V_rel  = abs(u_region(2) - u_region(1))*C_b	!!! Relative velocity of slab compared to background <--> NOTE: u_region is M_b, the Mach number with respect to the background
	Vx_b   = 0.0					!!! unperturbed background longitudinal velocity is ZERO in the analytical solution
	Vx_s   = V_rel					!!! unperturbed slab longitudinal velocity is thus the relative velocity
	Vz_b   = 0.0					!!! unperturbed transverse velocity is ZERO
	Vz_s   = 0.0					!!! unperturbed transverse velocity is ZERO
	var_b  = var_region(1,1)			!!! color (passive scalar) for background
	var_s  = var_region(2,1)			!!! color (passive scalar) for slab

	!!! For smoothing
	dx_levelmax = 0.5d0**nlevelmax
	sigma       = sigma_smooth * dx_levelmax
	gauss_sigma = gauss_sigma_smooth * dx_levelmax

	if( num_khi_eigenmodes .gt. 0 ) then
		if( .not. flat ) then
			!!! Calculate total amplitude of surface perturbation
			do j=1,num_khi_eigenmodes
				lambda = pert_lambda_over_diam(j) * 2.0 * rad
				kx     = 2.0 * acos(-1.d0) / lambda
				omega  = pert_omega(j)
				P_pert = pert_normalized_amp(j) * P_b		!!! A in my analytic language

				Q_b    = sqrt( kx**2 - ( (omega - kx*Vx_b) / C_b )**2 )
				Amp_up = -P_pert * Q_b / (D_b * omega**2)	!!! h in my analytic language - from A

				do i=1,nn
					Amp_up_tot(i,j) = real( Amp_up * exp( imag * kx * x(i,1) ) )
				end do
				if( pert_mode(j) .eq. 1 ) then		!!! S modes (fiducial)
					Amp_down_over_Amp_up = 1.0	!!! perturbation in slab upper vs lower boundaries. Should be +/-1.0 for Sinusoidal/Pinch modes respectively
				else if( pert_mode(j) .eq. 2 ) then	!!! P modes
					Amp_down_over_Amp_up = -1.0
				end if
				Amp_down_tot(:,j) = Amp_down_over_Amp_up * Amp_up_tot(:,j)		!!! Amplitude of perturbation in slab boundary, on the lower face	
			end do
		else
			Amp_up_tot(:,:) = 0.0
			Amp_down_tot(:,:) = 0.0
		end if
	else
		if( perturbed_var .eq. 0 ) then
			lambda   = lambda_over_diam * 2.0 * rad
			kx       = 2.0 * acos(-1.d0) / lambda
			do i=1,nn
				Amp_up_tot(i,1) = gauss_sigma * cos( kx * x(i,1) )
			end do
			Amp_down_tot(:,1) = -1.0 * Amp_up_tot(:,1)		!!! Symmetric P mode
		else
			Amp_up_tot(:,:) = 0.0
			Amp_down_tot(:,:) = 0.0
		end if
	end if

!	write(*,*) 'max x, min x'
!	write(*,*) maxval(x(1:nn,1)), minval(x(1:nn,1))
!	write(*,*) 'max H+, min H+, max H-, min H-'
!	write(*,*) maxval(Amp_up_tot(1:nn,1:num_khi_eigenmodes)), minval(Amp_up_tot(1:nn,1:num_khi_eigenmodes)), maxval(Amp_down_tot(1:nn,1:num_khi_eigenmodes)), minval(Amp_down_tot(1:nn,1:num_khi_eigenmodes))
	!!! Add unperturbed
	j = max(1, num_khi_eigenmodes)
	if(smooth .eq. 0) then
		do i=1,nn
			zz_up    = z_up   + sum(Amp_up_tot(i,1:j))
			zz_down  = z_down + sum(Amp_down_tot(i,1:j))
			if( x(i,2) .gt. zz_up ) then
				q(i,id) = D_b
				q(i,iu) = Vx_b
				q(i,iv) = Vz_b
				q(i,ip) = P_b
				q(i,ic) = var_b
			else if( x(i,2) .lt. zz_down ) then
				q(i,id) = D_b
				q(i,iu) = Vx_b
				q(i,iv) = Vz_b
				q(i,ip) = P_b
				q(i,ic) = var_b
			else
				q(i,id) = D_s
				q(i,iu) = Vx_s
				q(i,iv) = Vz_s
				q(i,ip) = P_s
				q(i,ic) = var_s
			end if
		end do
	elseif(smooth .eq. 1 .or. smooth .eq. 2) then
		do i=1,nn
			zz_up    = z_up   + sum(Amp_up_tot(i,1:j))
			zz_down  = z_down + sum(Amp_down_tot(i,1:j))

			ramp = 0.25 * ( 1 + tanh((x(i,2)-zz_down)/sigma) ) * ( 1 + tanh((zz_up-x(i,2))/sigma) )
			q(i,id) = D_b   + ramp * (D_s   - D_b)
			q(i,iu) = Vx_b  + ramp * (Vx_s  - Vx_b)
			q(i,iv) = Vz_b  + ramp * (Vz_s  - Vz_b)
			q(i,ip) = P_b   + ramp * (P_s   - P_b)
			q(i,ic) = var_b + ramp * (var_s - var_b)
		end do
	end if
	if(num_khi_eigenmodes .gt. 0) then
		!!! Calculate total perturbations
		do j=1,num_khi_eigenmodes
			if( pert_mode(j) .eq. 1 ) then		!!! S modes (fiducial)
				Amp_down_over_Amp_up = 1.0	!!! perturbation in slab upper vs lower boundaries. Should be +/-1.0 for Sinusoidal/Pinch modes respectively
			elseif( pert_mode(j) .eq. 2 ) then	!!! P modes
				Amp_down_over_Amp_up = -1.0
			end if

			lambda   = pert_lambda_over_diam(j) * 2.0 * rad
			kx       = 2.0 * acos(-1.d0) / lambda
			omega    = pert_omega(j)

			!!! complex wavenumbers and perturbation amplitudes !!!
			Q_b = sqrt( kx**2 - ( (omega - kx*Vx_b) / C_b )**2 )
			Q_s = sqrt( kx**2 - ( (omega - kx*Vx_s) / C_s )**2 )

			P_pert_up   =  pert_normalized_amp(j) * P_b		!!! A in my analytic language
			P_pert_down = -Amp_down_over_Amp_up * P_pert_up		!!! D in my analytic language
			P_pert_mid  =  P_pert_up

			D_pert_up   = -P_pert_up   * (Q_b**2 - kx**2) / ( (omega - kx*Vx_b)**2 )
			D_pert_down = -P_pert_down * (Q_b**2 - kx**2) / ( (omega - kx*Vx_b)**2 )
			D_pert_mid  = -P_pert_mid  * (Q_s**2 - kx**2) / ( (omega - kx*Vx_s)**2 )

			Vx_pert_up   = -P_pert_up   / ( D_b * (Vx_b - omega/kx) )
			Vx_pert_down = -P_pert_down / ( D_b * (Vx_b - omega/kx) )
			Vx_pert_mid  = -P_pert_mid  / ( D_s * (Vx_s - omega/kx) )

			Vz_pert_up   = -imag * Q_b * P_pert_up   / ( D_b * (kx*Vx_b - omega) )
			Vz_pert_down =  imag * Q_b * P_pert_down / ( D_b * (kx*Vx_b - omega) )
                        if( pert_mode(j) .eq. 1 ) then          !!! S modes (fiducial)
				Vz_pert_mid  = ( imag * Q_s * P_pert_mid / tanh(Q_s*rad) ) / ( D_s * (kx*Vx_s - omega) )
                        elseif( pert_mode(j) .eq. 2 ) then          !!! P modes
				Vz_pert_mid  = ( imag * Q_s * P_pert_mid * tanh(Q_s*rad) ) / ( D_s * (kx*Vx_s - omega) )
                        end if
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			do i=1,nn
				xx = x(i,1)
				zz = x(i,2) - 0.5

				zz_up   = z_up   + sum(Amp_up_tot(i,1:num_khi_eigenmodes))	! z_up = 0.5 + rad
				zz_down = z_down + sum(Amp_down_tot(i,1:num_khi_eigenmodes))	! z_down = 0.5 - rad

				ramp_up   = 0.5 * ( 1 + tanh((x(i,2)-zz_up)  /sigma) )	! 0 below, 0 in, 1 above the slab
				ramp_down = 0.5 * ( 1 + tanh((x(i,2)-zz_down)/sigma) ) 	! 0 below, 1 in, 1 above the slab

				!!! Add perturbations
				if( x(i,2) .gt. zz_up ) then
					expon1 = exp( imag * kx * xx - Q_b * (x(i,2)-zz_up) ) ! zz-rad or x(i,2)-zz_up=zz-(rad+Amp_up) --> Both include exp(Qb*rad)
					expon2 = expon1
					if(smooth .eq. 0 .or. smooth .eq. 1) then
						P_pert  = P_pert_up
						D_pert  = D_pert_up
						Vx_pert = Vx_pert_up
						Vz_pert = Vz_pert_up
					elseif(smooth .eq. 2) then
						P_pert  = P_pert_mid  + ramp_up * (P_pert_up  - P_pert_mid)
						D_pert  = D_pert_mid  + ramp_up * (D_pert_up  - D_pert_mid)
						Vx_pert = Vx_pert_mid + ramp_up * (Vx_pert_up - Vx_pert_mid)
						Vz_pert = Vz_pert_mid + ramp_up * (Vz_pert_up - Vz_pert_mid)
					end if
				else if( x(i,2) .lt. zz_down ) then
					expon1 = exp( imag * kx * xx + Q_b * (x(i,2)-zz_down) ) ! zz+rad or x(i,2)-zz_down=zz+(rad-Amp_down) --> include exp(Qb*rad)
					expon2 = expon1
					if(smooth .eq. 0 .or. smooth .eq. 1) then
						P_pert  = P_pert_down
						D_pert  = D_pert_down
						Vx_pert = Vx_pert_down
						Vz_pert = Vz_pert_down
					elseif(smooth .eq. 2) then
						if( pert_mode(j) .eq. 1 ) then		!!! S modes (fiducial)
							P_pert  = P_pert_down  + ramp_down * (-P_pert_mid  - P_pert_down )
							D_pert  = D_pert_down  + ramp_down * (-D_pert_mid  - D_pert_down )
							Vx_pert = Vx_pert_down + ramp_down * (-Vx_pert_mid - Vx_pert_down)
							Vz_pert = Vz_pert_down + ramp_down * ( Vz_pert_mid - Vz_pert_down)
						elseif( pert_mode(j) .eq. 2 ) then	!!! P modes
							P_pert  = P_pert_down  + ramp_down * ( P_pert_mid  - P_pert_down )
							D_pert  = D_pert_down  + ramp_down * ( D_pert_mid  - D_pert_down )
							Vx_pert = Vx_pert_down + ramp_down * ( Vx_pert_mid - Vx_pert_down)
							Vz_pert = Vz_pert_down + ramp_down * (-Vz_pert_mid - Vz_pert_down)
						end if
					end if
				else
					if( pert_mode(j) .eq. 1 ) then		!!! S modes (fiducial)
						expon1 = exp( imag * kx * xx ) * sinh( Q_s * abs(zz) ) / sinh(Q_s * rad)
						expon2 = exp( imag * kx * xx ) * cosh( Q_s * abs(zz) ) / cosh(Q_s * rad)
					elseif( pert_mode(j) .eq. 2 ) then	!!! P modes
						expon1 = exp( imag * kx * xx ) * cosh( Q_s * abs(zz) ) / cosh(Q_s * rad)
						expon2 = exp( imag * kx * xx ) * sinh( Q_s * abs(zz) ) / sinh(Q_s * rad)
					end if
					if(smooth .eq. 0 .or. smooth .eq. 1) then
						if( zz .ge. 0 ) then
							P_pert  = P_pert_mid
							D_pert  = D_pert_mid
							Vx_pert = Vx_pert_mid
							Vz_pert = Vz_pert_mid
						else
							if( pert_mode(j) .eq. 1 ) then		!!! S modes (fiducial)
								P_pert  = -P_pert_mid
								D_pert  = -D_pert_mid
								Vx_pert = -Vx_pert_mid
								Vz_pert =  Vz_pert_mid
							elseif( pert_mode(j) .eq. 2 ) then	!!! P modes
								P_pert  =  P_pert_mid
								D_pert  =  D_pert_mid
								Vx_pert =  Vx_pert_mid
								Vz_pert = -Vz_pert_mid
							end if
						end if
					elseif(smooth .eq. 2) then
						if( zz .ge. 0 ) then
							P_pert  = P_pert_mid  + ramp_up * (P_pert_up  - P_pert_mid)
							D_pert  = D_pert_mid  + ramp_up * (D_pert_up  - D_pert_mid)
							Vx_pert = Vx_pert_mid + ramp_up * (Vx_pert_up - Vx_pert_mid)
							Vz_pert = Vz_pert_mid + ramp_up * (Vz_pert_up - Vz_pert_mid)
						else
							if( pert_mode(j) .eq. 1 ) then		!!! S modes (fiducial)
								P_pert  = P_pert_down  + ramp_down * (-P_pert_mid  - P_pert_down )
								D_pert  = D_pert_down  + ramp_down * (-D_pert_mid  - D_pert_down )
								Vx_pert = Vx_pert_down + ramp_down * (-Vx_pert_mid - Vx_pert_down)
								Vz_pert = Vz_pert_down + ramp_down * ( Vz_pert_mid - Vz_pert_down)
							elseif( pert_mode(j) .eq. 2 ) then	!!! P modes
								P_pert  = P_pert_down  + ramp_down * ( P_pert_mid  - P_pert_down )
								D_pert  = D_pert_down  + ramp_down * ( D_pert_mid  - D_pert_down )
								Vx_pert = Vx_pert_down + ramp_down * ( Vx_pert_mid - Vx_pert_down)
								Vz_pert = Vz_pert_down + ramp_down * (-Vz_pert_mid - Vz_pert_down)
							end if
						end if
					end if
				end if
				q(i,ip) = q(i,ip) + real(P_pert  * expon1)
				q(i,id) = q(i,id) + real(D_pert  * expon1)
				q(i,iu) = q(i,iu) + real(Vx_pert * expon1)
				q(i,iv) = q(i,iv) + real(Vz_pert * expon2)
			end do
		end do
	elseif( num_khi_eigenmodes .eq. 0 .and. perturbed_var .ne. 0 ) then
		if( perturbed_var .eq. 1 ) then
			non_eigen_pert = normalized_amp * D_b
		elseif( perturbed_var .eq. 2 .or. perturbed_var .eq. 3 ) then
			non_eigen_pert = normalized_amp * V_rel
		elseif( perturbed_var .eq. 4 ) then
			non_eigen_pert = normalized_amp * P_b
		else
			non_eigen_pert = 0.0
		end if
		lambda   = lambda_over_diam * 2.0 * rad
		kx       = 2.0 * acos(-1.d0) / lambda
		do i=1,nn
			ramp = exp( -1.d0 * (x(i,2)-z_down)**2 / (2*gauss_sigma**2) ) + exp( -1.d0 * (x(i,2)-z_up)**2 / (2*gauss_sigma**2) )
			q(i,perturbed_var) = q(i,perturbed_var) + ramp * non_eigen_pert * cos(kx*x(i,1))
		end do
	end if

	! Convert primitive to conservative variables
	! density -> density
	u(1:nn,1)=q(1:nn,1)
	! velocity -> momentum
	u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
	u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
	! kinetic energy
	u(1:nn,ndim+2)=0.0d0
	u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,2)**2	!!! 0.5 * rho * Vx^2
	u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,3)**2	!!! 0.5 * rho * Vy^2
	! pressure -> total fluid energy
	u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)	!!! P / (gamma - 1)
	! passive scalars
	do ivar=ndim+3,nvar
		u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
	end do

end subroutine condinit
