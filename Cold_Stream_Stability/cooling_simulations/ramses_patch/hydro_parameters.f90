module hydro_parameters
  use amr_parameters

  ! Number of independant variables
#ifndef NENER
  integer,parameter::nener=0
#else
  integer,parameter::nener=NENER
#endif
#ifndef NVAR
  integer,parameter::nvar=ndim+2+nener
#else
  integer,parameter::nvar=NVAR
#endif
  ! Size of hydro kernel
  integer,parameter::iu1=-1
  integer,parameter::iu2=+4
  integer,parameter::ju1=(1-ndim/2)-1*(ndim/2)
  integer,parameter::ju2=(1-ndim/2)+4*(ndim/2)
  integer,parameter::ku1=(1-ndim/3)-1*(ndim/3)
  integer,parameter::ku2=(1-ndim/3)+4*(ndim/3)
  integer,parameter::if1=1
  integer,parameter::if2=3
  integer,parameter::jf1=1
  integer,parameter::jf2=(1-ndim/2)+3*(ndim/2)
  integer,parameter::kf1=1
  integer,parameter::kf2=(1-ndim/3)+3*(ndim/3)

  ! Imposed boundary condition variables
  real(dp),dimension(1:MAXBOUND,1:nvar)::boundary_var
  real(dp),dimension(1:MAXBOUND)::d_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::p_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::u_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::v_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::w_bound=0.0d0
#if NENER>0
  real(dp),dimension(1:MAXBOUND,1:NENER)::prad_bound=0.0
#endif
#if NVAR>NDIM+2+NENER
  real(dp),dimension(1:MAXBOUND,1:NVAR-NDIM-2-NENER)::var_bound=0.0
#endif
  ! Refinement parameters for hydro
  real(dp)::err_grad_d=-1.0  ! Density gradient
  real(dp)::err_grad_u=-1.0  ! Velocity gradient
  real(dp)::err_grad_p=-1.0  ! Pressure gradient
  real(dp)::floor_d=1.d-10   ! Density floor
  real(dp)::floor_u=1.d-10   ! Velocity floor
  real(dp)::floor_p=1.d-10   ! Pressure floor
  real(dp)::mass_sph=0.0D0   ! mass_sph
#if NENER>0
  real(dp),dimension(1:NENER)::err_grad_prad=-1.0
#endif
#if NVAR>NDIM+2+NENER
  real(dp),dimension(1:NVAR-NDIM-2)::err_grad_var=-1.0
#endif
  real(dp),dimension(1:MAXLEVEL)::jeans_refine=-1.0

  ! Initial conditions hydro variables
  real(dp),dimension(1:MAXREGION)::d_region=0.
  real(dp),dimension(1:MAXREGION)::u_region=0.
  real(dp),dimension(1:MAXREGION)::v_region=0.
  real(dp),dimension(1:MAXREGION)::w_region=0.
  real(dp),dimension(1:MAXREGION)::p_region=0.
  real(dp),dimension(1:MAXREGION)::met_region=0.
  real(dp),dimension(1:MAXREGION)::col_region=0.
#if NENER>0
  real(dp),dimension(1:MAXREGION,1:NENER)::prad_region=0.0
#endif
#if NVAR>NDIM+2+NENER
  real(dp),dimension(1:MAXREGION,1:NVAR-NDIM-2-NENER)::var_region=0.0
#endif

  ! Kelvin Helmholtz
  real(dp)::Rstream=0.0 				! Stream radius / Slab half-thickness
  integer::Nmode_tot=0					! Total number of modes. THIS IS NOT READ FROM THE NAMELIST, BUT RATHER COMPUTED BASED ON OTHER ENTRIES
  integer::Nwavelength=0				! Number of wavelengths per symmetry mode. THIS IS NOT READ FROM THE NAMELIST, BUT RATHER COMPUTED BASED ON OTHER ENTRIES
  integer,parameter::MAXPERTURB=1000			! Including all wavelengths and all symmetry modes
  integer,parameter::MAXSYMMETRY=10  			! Number of allowed different symmetry modes.
  integer::num_symmetry_modes=0				! Number of different symmetry modes
  integer,dimension(1:MAXSYMMETRY)::symmetry_modes=-1	! 0-Pinch (P), 1-Helical (S), 2-Elliptical, 3-Triangular, ...
  integer,dimension(1:3)::wavenumber_range=0		! (jmin, Delta_j, jmax) --> k=2*Pi*j --> lambda=1/j --> Range of perturbed wavenumbers for each symmetry mode
  integer::perturbed_var=-1				! Same for all perturbations in a run: 0-Interface, 1-Density, 2-V_long, 3-V_trans, 4-Pressure, 5-Eigenmode
  real(dp)::normalized_amp=0.0				! Normalized amplitude of the perturbation with jmin (largest wavelength). Larger j perturbations determined by PS_slope
  real(dp)::PS_slope=0.0				! Slope of power law determining perturbation amplitude as a function of wavenumber
  integer::phase_seed=0					! seed for random phases. 0 is default for random.f90
  real(dp),dimension(1:MAXPERTURB)::pert_phase=0.0	! Phases for different perturbations, in the range [0,2*pi). Different for each perturbation (wavenumbers and symmetry modes)
  real(dp),dimension(1:MAXPERTURB)::real_omega=0.0	! For Eigenmode perturbations, solutions for complex frequency
  real(dp),dimension(1:MAXPERTURB)::imag_omega=0.0	! For Eigenmode perturbations, solutions for complex frequency
  integer::smooth=-1					! 0 for no smoothing, 1 for smoothing unperturbed profile only, 2 for smoothing eigenmode perturbations as well
  logical::flat=.false.					! For Eigenmodes only, whether to perturb the interface as well.
  real(dp)::smooth_sigma=0				! Width of smoothing function for unperturbed setup and/or eigenmode perturbations
  real(dp)::gauss_sigma=0				! Width of Gaussian defining perturbation region for non-eigenmode perturbations
  real(dp)::Tmax_cool=1.d5				! Shut off cooling above this temperature
  real(dp)::z_UV=2.0					! Redshift for the Haardt and Madau UV background

  ! Hydro solver parameters
  integer ::niter_riemann=10
  integer ::slope_type=1
  real(dp)::gamma=1.4d0
  real(dp),dimension(1:512)::gamma_rad=1.33333333334d0
  real(dp)::courant_factor=0.5d0
  real(dp)::difmag=0.0d0
  real(dp)::smallc=1.d-10
  real(dp)::smallr=1.d-10
  character(LEN=10)::scheme='muscl'
  character(LEN=10)::riemann='llf'

  ! Interpolation parameters
  integer ::interpol_var=0
  integer ::interpol_type=1

  ! Passive variables index
  integer::imetal=6
  integer::idelay=6
  integer::ixion=6
  integer::ichem=6

end module hydro_parameters
