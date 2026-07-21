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
	integer  :: id, iu, iv, iw, ip, iz, ic
	real(dp),dimension(1:nvector,1:nvar),save::q          ! Primitive variables
	real(dp),dimension(1:nvector),save::xx,yy,zz,rr,phi   ! Centred coordinates
	real(dp),dimension(1:nvector,1:MAXPERTURB),save::Radius_pert
	integer  :: i, j, k, ivar, nmode
	real(dp) :: Pi
	real(dp) :: dx_levelmax, sigma_smooth, sigma_gauss, ramp
	real(dp) :: klong, rad_up, pert_amp, pert_amp0
	real(dp) :: V_rel, D_b, D_s, P_b, P_s, Vlong_b, Vlong_s, Vtrans_b, Vtrans_s, C_b, C_s, col_b, col_s, met_b, met_s	!!! Unperturbed variables
	complex(dp) :: imag=(0,1), omega, Q_b, Q_s, expon1, expon2, Amp_up
	complex(dp) :: P_pert_out, P_pert_in, P_pert, D_pert_out, D_pert_in, D_pert
	complex(dp) :: Vlong_pert_out, Vlong_pert_in, Vlong_pert	! Vz
	complex(dp) :: Vtrans_pert_out, Vtrans_pert_in, Vtrans_pert	! 2d: Vx, 3d: Vr

	Pi = acos(-1.d0)
	id = 1
	iu = 2
	iv = 3
	iw = 4
	ip = ndim+2
	iz = ndim+3
	ic = ndim+4
	!!! Centre coordinates
#if NDIM==2
	zz(1:nn) = x(1:nn,1)
	xx(1:nn) = x(1:nn,2) - 0.5
	rr(1:nn) = abs(xx(1:nn))
	do i=1,nn
		if( xx(i) .ge. 0.0 ) then
			phi(i) = 0.0
		else
			phi(i) = Pi
		end if
	end do
#endif
#if NDIM==3
	zz(1:nn) = x(1:nn,1)
	xx(1:nn) = x(1:nn,2) - 0.5
	yy(1:nn) = x(1:nn,3) - 0.5
	rr(1:nn) = sqrt( xx(1:nn)**2 + yy(1:nn)**2 )
	do i=1,nn
		if( xx(i) .eq. 0.0 .and. yy(i) .eq. 0.0 ) then
			phi(i) = 0.0
		elseif( xx(i) .ge. 0.0 .and. yy(i) .ge. 0.0 ) then
			phi(i) = atan( yy(i) / xx(i) )
		elseif( xx(i) .lt. 0.0 ) then
			phi(i) = Pi + atan( yy(i) / xx(i) )
		else
			phi(i) = 2.0*Pi + atan( yy(i) / xx(i) )
		end if
	end do
#endif
	!!! Unperturbed variables
	D_b    = d_region(1)				!!! unperturbed background density
	D_s    = d_region(2)				!!! unperturbed slab density
	P_b    = p_region(1)				!!! unperturbed background pressure
	P_s    = p_region(2)				!!! unperturbed slab pressure
	C_b    = sqrt(gamma * P_b / D_b)		!!! Speed of sound in background
	C_s    = sqrt(gamma * P_s / D_s)		!!! Speed of sound in slab
	V_rel  = abs(u_region(2) - u_region(1))*C_b	!!! Relative velocity of slab compared to background <--> NOTE: u_region is M_b, the Mach number with respect to the background
	Vlong_b   = 0.0					!!! unperturbed background longitudinal velocity is ZERO in the analytical solution
	Vlong_s   = V_rel				!!! unperturbed slab longitudinal velocity is thus the relative velocity
	Vtrans_b   = 0.0				!!! unperturbed transverse velocity is ZERO
	Vtrans_s   = 0.0				!!! unperturbed transverse velocity is ZERO
	col_b  = col_region(1)				!!! color (passive scalar) for background
	col_s  = col_region(2)				!!! color (passive scalar) for slab
	met_b  = met_region(1)				!!! metalicity (passive scalar) for background
	met_s  = met_region(2)				!!! metalicity (passive scalar) for slab

	!!! For smoothing
	dx_levelmax  = 0.5d0**nlevelmax
	sigma_smooth = smooth_sigma * Rstream
	sigma_gauss  = gauss_sigma  * Rstream

	!!! Set up interface perturbations
	pert_amp0 = 0.0
	do k=1,num_symmetry_modes
		do j=1,Nwavelength
			pert_amp0 = pert_amp0 + ( ( 1.0 + (j-1)*wavenumber_range(2)/wavenumber_range(1) )**(PS_slope) )**2
		end do
	end do
	pert_amp0 = sqrt(pert_amp0)
	if( perturbed_var .eq. 5 ) then			! Eigenmodes - currently supports 2d only
		if( .not. flat ) then
			!!! Calculate total amplitude of surface perturbation
			do k=1,num_symmetry_modes
				do j=1,Nwavelength
					nmode = j + (k-1)*Nwavelength
					pert_amp = (normalized_amp/pert_amp0) * ( 1.0 + (j-1)*wavenumber_range(2)/wavenumber_range(1) )**(PS_slope)
					klong  = 2.0 * Pi * (wavenumber_range(1) + (j-1)*wavenumber_range(2))
					omega  = real_omega(nmode) + imag*imag_omega(nmode)
					P_pert = pert_amp * P_b							!!! A in Mandelker+2016 language
					Q_b    = sqrt( klong**2 - ( (omega - klong*Vlong_b) / C_b )**2 )
					Amp_up = -P_pert * Q_b / (D_b * omega**2)				!!! h in Mandelker+2016 language - from A
					Radius_pert(1:nn,nmode) = real( Amp_up * exp( imag * (klong*zz(1:nn) + symmetry_modes(k)*phi(1:nn) + pert_phase(nmode)) ) )
				end do
			end do
		else
			Radius_pert(1:nn,1:Nmode_tot) = 0.0
		end if
	elseif( perturbed_var .eq. 0 ) then
		pert_amp = (normalized_amp/pert_amp0) * ( 1.0 + (j-1)*wavenumber_range(2)/wavenumber_range(1) )**(PS_slope) 
		pert_amp = pert_amp * Rstream
		do k=1,num_symmetry_modes
			do j=1,Nwavelength
				nmode = j + (k-1)*Nwavelength
				klong  = 2.0 * Pi * (wavenumber_range(1) + (j-1)*Wavenumber_range(2))
				Radius_pert(1:nn,nmode) = pert_amp * cos( klong*zz(1:nn) + symmetry_modes(k)*phi(1:nn) + pert_phase(nmode) )
			end do
		end do
	else
		Radius_pert(1:nn,1:Nmode_tot) = 0.0
	end if

	!!! Add unperturbed values
	if(smooth .eq. 0) then
		do i=1,nn
			rad_up = Rstream + sum(Radius_pert(i,1:Nmode_tot))
			if( rr(i) .gt. rad_up ) then
				q(i,id) = D_b
				q(i,iu) = Vlong_b
				q(i,iv) = Vtrans_b
				q(i,iw) = Vtrans_b
				q(i,ip) = P_b
				q(i,iz) = met_b
				q(i,ic) = col_b
			else
				q(i,id) = D_s
				q(i,iu) = Vlong_s
				q(i,iv) = Vtrans_s
				q(i,iw) = Vtrans_s
				q(i,ip) = P_s
				q(i,iz) = met_s
				q(i,ic) = col_s
			end if
		end do
	elseif(smooth .eq. 1 .or. smooth .eq. 2) then
		do i=1,nn
			rad_up = Rstream + sum(Radius_pert(i,1:Nmode_tot))
			ramp = 0.5 * ( 1 + tanh( (rad_up-rr(i)) / sigma_smooth ) )	! 0 outside r>Rs, 1 inside r<Rs
			q(i,id) = D_b      + ramp * (D_s      - D_b)
			q(i,iu) = Vlong_b  + ramp * (Vlong_s  - Vlong_b)
			q(i,iv) = Vtrans_b + ramp * (Vtrans_s - Vtrans_b)
			q(i,iw) = Vtrans_b + ramp * (Vtrans_s - Vtrans_b)
			q(i,ip) = P_b      + ramp * (P_s      - P_b)
			q(i,iz) = met_b    + ramp * (met_s    - met_b)
			q(i,ic) = col_b    + ramp * (col_s    - col_b)
		end do
	end if

	!!! Calculate total perturbations
	if( perturbed_var .eq. 5 ) then			! Eigenmodes - currently supports 2d only
		do k=1,num_symmetry_modes
			do j=1,Nwavelength
				nmode = j + (k-1)*Nwavelength
				pert_amp = (normalized_amp/pert_amp0) * ( 1.0 + (j-1)*wavenumber_range(2)/wavenumber_range(1) )**(PS_slope)
				klong  = 2.0 * Pi * (wavenumber_range(1) + (j-1)*wavenumber_range(2))
				omega  = real_omega(nmode) + imag*imag_omega(nmode)

				!!! complex wavenumbers and perturbation amplitudes !!!
				Q_b = sqrt( klong**2 - ( (omega - klong*Vlong_b) / C_b )**2 )
				Q_s = sqrt( klong**2 - ( (omega - klong*Vlong_s) / C_s )**2 )

				P_pert_out = pert_amp * P_b				!!! A in my analytic language
				P_pert_in  = P_pert_out

				D_pert_out = -P_pert_out * (Q_b**2 - klong**2) / ( (omega - klong*Vlong_b)**2 )
				D_pert_in  = -P_pert_in  * (Q_s**2 - klong**2) / ( (omega - klong*Vlong_s)**2 )

				Vlong_pert_out = -P_pert_out / ( D_b * (Vlong_b - omega/klong) )
				Vlong_pert_in  = -P_pert_in  / ( D_s * (Vlong_s - omega/klong) )

				Vtrans_pert_out = -imag * Q_b * P_pert_out / ( D_b * (klong*Vlong_b - omega) )
				Vtrans_pert_in  =  imag * Q_s * P_pert_in  / ( D_s * (klong*Vlong_s - omega) )
	                        if( mod(symmetry_modes(k),2) .eq. 1 ) then		! S modes
					Vtrans_pert_in  = Vtrans_pert_in * ( exp(2.0*Q_s*Rstream) + 1.0 ) / ( exp(2.0*Q_s*Rstream) - 1.0 )
	                        else							! P modes
					Vtrans_pert_in  = Vtrans_pert_in * ( exp(2.0*Q_s*Rstream) - 1.0 ) / ( exp(2.0*Q_s*Rstream) + 1.0 )
	                        end if
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				do i=1,nn
					rad_up = Rstream + sum(Radius_pert(i,1:Nmode_tot))
					ramp = 0.5 * ( 1 + tanh( (rad_up-rr(i)) / sigma_smooth ) )	! 0 outside r>Rs, 1 inside r<Rs
					if( rr(i) .gt. rad_up ) then
						expon1 = exp( imag * (klong*zz(i) + symmetry_modes(k)*phi(i)     + pert_phase(nmode)) - Q_b * (rr(i)-rad_up) )
						expon2 = exp( imag * (klong*zz(i) + (symmetry_modes(k)+1)*phi(i) + pert_phase(nmode)) - Q_b * (rr(i)-rad_up) )
						if(smooth .eq. 0 .or. smooth .eq. 1) then
							P_pert      = P_pert_out
							D_pert      = D_pert_out
							Vlong_pert  = Vlong_pert_out
							Vtrans_pert = Vtrans_pert_out
						elseif(smooth .eq. 2) then
							P_pert      = P_pert_out      + ramp * (P_pert_in      - P_pert_out)
							D_pert      = D_pert_out      + ramp * (D_pert_in      - D_pert_out)
							Vlong_pert  = Vlong_pert_out  + ramp * (Vlong_pert_in  - Vlong_pert_out)
							Vtrans_pert = Vtrans_pert_out + ramp * (Vtrans_pert_in - Vtrans_pert_out)
						end if
					else
						if( mod(symmetry_modes(k),2) .eq. 1 ) then		! S modes
							expon1 = exp( imag * (klong*zz(i) + pert_phase(nmode)) ) * & 
								& ( exp(Q_s*xx(i)) - exp(-Q_s*xx(i)) ) / ( exp(Q_s*Rstream) - exp(-Q_s*Rstream) )
							expon2 = exp( imag * (klong*zz(i) + pert_phase(nmode)) ) * & 
								& ( exp(Q_s*xx(i)) + exp(-Q_s*xx(i)) ) / ( exp(Q_s*Rstream) + exp(-Q_s*Rstream) )
						else							! P modes
							expon1 = exp( imag * (klong*zz(i) + pert_phase(nmode)) ) * & 
								& ( exp(Q_s*xx(i)) + exp(-Q_s*xx(i)) ) / ( exp(Q_s*Rstream) + exp(-Q_s*Rstream) )
							expon2 = exp( imag * (klong*zz(i) + pert_phase(nmode)) ) * & 
								& ( exp(Q_s*xx(i)) - exp(-Q_s*xx(i)) ) / ( exp(Q_s*Rstream) - exp(-Q_s*Rstream) )
						end if
						if(smooth .eq. 0 .or. smooth .eq. 1) then
							P_pert  = P_pert_in
							D_pert  = D_pert_in
							Vlong_pert = Vlong_pert_in
							Vtrans_pert = Vtrans_pert_in
						elseif(smooth .eq. 2) then
							P_pert      = P_pert_out      + ramp * (P_pert_in      - P_pert_out)
							D_pert      = D_pert_out      + ramp * (D_pert_in      - D_pert_out)
							Vlong_pert  = Vlong_pert_out  + ramp * (Vlong_pert_in  - Vlong_pert_out)
							Vtrans_pert = Vtrans_pert_out + ramp * (Vtrans_pert_in - Vtrans_pert_out)
						end if
					end if
					q(i,ip) = q(i,ip) + real(P_pert  * expon1)
					q(i,id) = q(i,id) + real(D_pert  * expon1)
					q(i,iu) = q(i,iu) + real(Vlong_pert * expon1)
					q(i,iv) = q(i,iv) + real(Vtrans_pert * expon2)
				end do
			end do
		end do
	elseif( perturbed_var .ne. 0 ) then
		do k=1,num_symmetry_modes
			do j=1,Nwavelength
				nmode    = j + (k-1)*Nwavelength
				pert_amp = (normalized_amp/pert_amp0) * ( 1.0 + (j-1)*wavenumber_range(2)/wavenumber_range(1) )**(PS_slope)
				klong    = 2.0 * Pi * (wavenumber_range(1) + (j-1)*wavenumber_range(2))
				if( perturbed_var .eq. 1 ) then
					pert_amp = pert_amp * D_b
					do i=1,nn
						ramp = exp( -1.d0 * (rr(i)-Rstream)**2 / (2*sigma_gauss**2) )
						q(i,id) = q(i,id) + ramp * pert_amp * cos( klong*zz(i) + symmetry_modes(k)*phi(i) + pert_phase(nmode) )
					end do
				elseif( perturbed_var .eq. 2 ) then
					pert_amp = pert_amp * V_rel
					do i=1,nn
						ramp = exp( -1.d0 * (rr(i)-Rstream)**2 / (2*sigma_gauss**2) )
						q(i,iu) = q(i,iu) + ramp * pert_amp * cos( klong*zz(i) + symmetry_modes(k)*phi(i) + pert_phase(nmode) )
					end do
				elseif( perturbed_var .eq. 3 ) then
					pert_amp = pert_amp * C_s
					do i=1,nn
						ramp = exp( -1.d0 * (rr(i)-Rstream)**2 / (2*sigma_gauss**2) )
						q(i,iv) = q(i,iv) + ramp * pert_amp * cos( klong*zz(i) + symmetry_modes(k)*phi(i) + pert_phase(nmode) ) * cos(phi(i))
						q(i,iw) = q(i,iw) + ramp * pert_amp * cos( klong*zz(i) + symmetry_modes(k)*phi(i) + pert_phase(nmode) ) * sin(phi(i))
					end do
				elseif( perturbed_var .eq. 4 ) then
					pert_amp = pert_amp * P_b
					do i=1,nn
						ramp = exp( -1.d0 * (rr(i)-Rstream)**2 / (2*sigma_gauss**2) )
						q(i,ip) = q(i,ip) + ramp * pert_amp * cos( klong*zz(i) + symmetry_modes(k)*phi(i) + pert_phase(nmode) )
					end do
				else
					pert_amp = 0.0
				end if
			end do
		end do
	end if

	! Convert primitive to conservative variables
	! density -> density
	u(1:nn,1)=q(1:nn,1)
	! velocity -> momentum
	u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
	u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#if NDIM>2
	u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
	! kinetic energy
	u(1:nn,ndim+2)=0.0d0
	u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,2)**2	!!! 0.5 * rho * Vlong^2
	u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,3)**2	!!! 0.5 * rho * Vtrans^2
#if NDIM>2
	u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,4)**2
#endif
	! pressure -> total fluid energy
	u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)	!!! P / (gamma - 1)
	! passive scalars
	do ivar=ndim+3,nvar
		u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
	end do

end subroutine condinit
