subroutine read_hydro_params(nml_ok)
  use amr_commons
  use hydro_commons
  use amr_parameters
  use hydro_parameters
  use pm_commons
  use pm_parameters
  Use random
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical::nml_ok
  !--------------------------------------------------
  ! Local variables  
  !--------------------------------------------------
  integer::i,idim,nboundary_true=0
  integer ,dimension(1:MAXBOUND)::bound_type
  real(dp)::scale,ek_bound
  integer,dimension( ncpu, IRandNumSize ) :: allseed

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/init_params/filetype,initfile,multiple,nregion,region_type &
       & ,x_center,y_center,z_center,aexp_ini &
       & ,length_x,length_y,length_z,exp_region &
#if NVAR>NDIM+2
       & ,var_region &
#endif
       & ,d_region,u_region,v_region,w_region,p_region,met_region,col_region &
       & ,Rstream,num_symmetry_modes,symmetry_modes &
       & ,wavenumber_range,perturbed_var,normalized_amp,PS_slope &
       & ,phase_seed,pert_phase,real_omega,imag_omega &
       & ,smooth,flat,smooth_sigma,gauss_sigma       
  namelist/hydro_params/gamma,courant_factor,smallr,smallc &
       & ,niter_riemann,slope_type,difmag &
       & ,pressure_fix,beta_fix,scheme,riemann
  namelist/refine_params/x_refine,y_refine,z_refine,r_refine &
       & ,a_refine,b_refine,exp_refine,jeans_refine,mass_cut_refine &
       & ,m_refine,mass_sph,err_grad_d,err_grad_p,err_grad_u &
       & ,floor_d,floor_u,floor_p,ivar_refine,var_cut_refine &
       & ,interpol_var,interpol_type
  namelist/boundary_params/nboundary,bound_type &
       & ,ibound_min,ibound_max,jbound_min,jbound_max &
       & ,kbound_min,kbound_max &
       & ,d_bound,u_bound,v_bound,w_bound,p_bound
  namelist/physics_params/cooling,haardt_madau,metal,isothermal,bondi &
       & ,m_star,t_star,n_star,T2_star,g_star,del_star,eps_star,jeans_ncells &
       & ,eta_sn,yield,rbubble,f_ek,ndebris,f_w,mass_gmc &
       & ,J21,a_spec,z_ave,z_reion,n_sink,bondi,delayed_cooling &
       & ,self_shielding,smbh,agn,rsink_max,msink_max &
       & ,units_density,units_time,units_length,Tmax_cool,z_UV

  ! Read namelist file
  rewind(1)
  read(1,NML=init_params,END=101)
  goto 102
101 write(*,*)' You need to set up namelist &INIT_PARAMS in parameter file'
  call clean_stop
102 rewind(1)
  if(nlevelmax>levelmin)read(1,NML=refine_params)
  rewind(1)
  if(hydro)read(1,NML=hydro_params)
  rewind(1)
  read(1,NML=boundary_params,END=103)
  simple_boundary=.true.
  goto 104
103 simple_boundary=.false.
104 if(nboundary>MAXBOUND)then
    write(*,*) 'Error: nboundary>MAXBOUND'
    call clean_stop
  end if
  rewind(1)
  read(1,NML=physics_params,END=105)
105 continue
#ifdef ATON
  if(aton)call read_radiation_params(1)
#endif

  !--------------------------------------------------
  ! Check for consistancy in KHI parameters and set phases
  !--------------------------------------------------
#if NDIM==1
	write(*,*)' ERROR: KHI requires Ndim=2 or Ndim=3'
	call clean_stop
#endif
	if( Rstream .lt. 4.0 / (2.0**nlevelmax) ) then
		write(*,*)' ERROR: The stream radius must be at least 4 high resolution cells'
		call clean_stop
	end if
	if( num_symmetry_modes .le. 0 ) then	
		write(*,*)' ERROR: You need to select number of symmetry modes in namelist file'
		call clean_stop
	end if
	if( num_symmetry_modes .gt. MAXSYMMETRY ) then
		write(*,*)' ERROR: Too many symmetry modes. Increase MAXSYMMETRY'
		call clean_stop
	end if
	if( minval(symmetry_modes(1:num_symmetry_modes)) .eq. -1 ) then
		write(*,*)' ERROR: You need to define which symmetry modes to excite in namelist file'
		call clean_stop
	end if
	if( wavenumber_range(1) .lt. 1 ) then
		write(*,*)' ERROR: Bad wavenumber range in namelist file. jmin<1'
		call clean_stop
	end if
	if( wavenumber_range(2) .lt. 1 ) then
		write(*,*)' ERROR: Bad wavenumber range in namelist file. Delta_j<1'
		call clean_stop
	end if
	if( wavenumber_range(3) .gt. floor(2.0**(nlevelmax-1)) ) then
		write(*,*)' ERROR: Bad wavenumber range in namelist file. jmax>2^(Lmax-1)'
		call clean_stop
	end if
	if( wavenumber_range(3) .lt. wavenumber_range(1) ) then
		write(*,*)' ERROR: Bad wavenumber range in namelist file. jmax<jmin'
		call clean_stop
	end if
	if( mod( wavenumber_range(3)-wavenumber_range(1), wavenumber_range(2) ) .ne. 0 ) then
		write(*,*)' ERROR: Bad wavenumber range in namelist file. Make sure that jmax-jmin+1 divides by Delta_j'
		call clean_stop
	end if
	Nwavelength = ( (wavenumber_range(3)-wavenumber_range(1)) / wavenumber_range(2) ) + 1
	Nmode_tot   = num_symmetry_modes * Nwavelength
	if( Nmode_tot .gt. MAXPERTURB ) then
		write(*,*)' ERROR: Too many perturbations. Increase MAXPERTURB'
		call clean_stop
	end if
	if( perturbed_var .lt. 0 .or. perturbed_var .gt. 5 ) then
		write(*,*)' ERROR: You need to select perturbed variable from 0-5 in namelist file'
		call clean_stop
	end if
	if( normalized_amp .le. 0.0 ) then
		write(*,*)' ERROR: You need to define perturbation amplitude in namelist file'
		call clean_stop
	end if
	if( minval(abs(pert_phase(1:Nmode_tot))) .lt. 1.d-4 .and. maxval(abs(pert_phase(1:Nmode_tot))) .lt. 1.d-4 ) then
		if(myid==1) write(*,*)' Generating random phases'
		call rans(ncpu,phase_seed,allseed)
		localseed = allseed(myid,1:IRandNumSize)
		do i=1,Nmode_tot
			call Ranf(localseed, pert_phase(i))
		end do
	else
		write(*,*)' Using pre-defined phases phases'
		if( minval(pert_phase(1:Nmode_tot)) .lt. 0.0 .or. maxval(pert_phase(1:Nmode_tot)) .ge. 1.0 ) then
			write(*,*)' ERROR: Phases must be in range [0, 1)*2*Pi'
			call clean_stop
		end if
	endif
	if(myid==1) write(*,*)' Nmodes= ',Nmode_tot
	pert_phase(1:Nmode_tot) = ( pert_phase(1:Nmode_tot) - minval(pert_phase(1:Nmode_tot)) ) * 2.0*acos(-1.d0)
	do i=1,Nmode_tot
		if(myid==1) write(*,*) pert_phase(i)
	end do
#if NDIM>2
	if( perturbed_var .eq. 5 ) then
		write(*,*)' ERROR: Sorry, the eigenmode option is currently only available for 2d runs. We are working to fix the issue. Good day.'
		call clean_stop
	end if
#endif
	if( perturbed_var .eq. 5 .and. minval(abs(imag_omega(1:Nmode_tot))) .eq. 0.0 .and. minval(abs(real_omega(1:Nmode_tot))) .eq. 0.0 ) then
		write(*,*)' ERROR: You selected eigenmodes but did not provide enough eigenfrequencies'
		call clean_stop
	end if
	if( smooth .lt. 0 .or. smooth .gt. 2 ) then
		write(*,*)' ERROR: You need to select smoothing from 0-2 in namelist file'
		call clean_stop
	end if
	if( smooth_sigma .le. 0 ) then
		write(*,*)' ERROR: Bad smoothing width in namelist file: <0'
		call clean_stop
	end if
	if( gauss_sigma .le. 0 ) then
		write(*,*)' ERROR: Bad perturbation width in namelist file: <0'
		call clean_stop
	end if

  !--------------------------------------------------
  ! Check for star formation
  !--------------------------------------------------
  if(t_star>0)then
     star=.true.
     pic=.true.
  else if(eps_star>0)then
     t_star=0.1635449*(n_star/0.1)**(-0.5)/eps_star
     star=.true.
     pic=.true.
  endif

  !--------------------------------------------------
  ! Check for metal
  !--------------------------------------------------
  if(metal.and.nvar<(ndim+3))then
     if(myid==1)write(*,*)'Error: metals need nvar >= ndim+3'
     if(myid==1)write(*,*)'Modify hydro_parameters.f90 and recompile'
     nml_ok=.false.
  endif

  !-------------------------------------------------
  ! This section deals with hydro boundary conditions
  !-------------------------------------------------
  if(simple_boundary.and.nboundary==0)then
     simple_boundary=.false.
  endif

  if (simple_boundary)then

     ! Compute new coarse grid boundaries
     do i=1,nboundary
        if(ibound_min(i)*ibound_max(i)==1.and.ndim>0.and.bound_type(i)>0)then
           nx=nx+1
           if(ibound_min(i)==-1)then
              icoarse_min=icoarse_min+1
              icoarse_max=icoarse_max+1
           end if
           nboundary_true=nboundary_true+1
        end if
     end do
     do i=1,nboundary
        if(jbound_min(i)*jbound_max(i)==1.and.ndim>1.and.bound_type(i)>0)then
           ny=ny+1
           if(jbound_min(i)==-1)then
              jcoarse_min=jcoarse_min+1
              jcoarse_max=jcoarse_max+1
           end if
           nboundary_true=nboundary_true+1
        end if
     end do
     do i=1,nboundary
        if(kbound_min(i)*kbound_max(i)==1.and.ndim>2.and.bound_type(i)>0)then
           nz=nz+1
           if(kbound_min(i)==-1)then
              kcoarse_min=kcoarse_min+1
              kcoarse_max=kcoarse_max+1
           end if
           nboundary_true=nboundary_true+1
        end if
     end do

     ! Compute boundary geometry
     do i=1,nboundary
        if(ibound_min(i)*ibound_max(i)==1.and.ndim>0.and.bound_type(i)>0)then
           if(ibound_min(i)==-1)then
              ibound_min(i)=icoarse_min+ibound_min(i)
              ibound_max(i)=icoarse_min+ibound_max(i)
              if(bound_type(i)==1)boundary_type(i)=1
              if(bound_type(i)==2)boundary_type(i)=11
              if(bound_type(i)==3)boundary_type(i)=21
           else
              ibound_min(i)=icoarse_max+ibound_min(i)
              ibound_max(i)=icoarse_max+ibound_max(i)
              if(bound_type(i)==1)boundary_type(i)=2
              if(bound_type(i)==2)boundary_type(i)=12
              if(bound_type(i)==3)boundary_type(i)=22
           end if
           if(ndim>1)jbound_min(i)=jcoarse_min+jbound_min(i)
           if(ndim>1)jbound_max(i)=jcoarse_max+jbound_max(i)
           if(ndim>2)kbound_min(i)=kcoarse_min+kbound_min(i)
           if(ndim>2)kbound_max(i)=kcoarse_max+kbound_max(i)
        else if(jbound_min(i)*jbound_max(i)==1.and.ndim>1.and.bound_type(i)>0)then
           ibound_min(i)=icoarse_min+ibound_min(i)
           ibound_max(i)=icoarse_max+ibound_max(i)
           if(jbound_min(i)==-1)then
              jbound_min(i)=jcoarse_min+jbound_min(i)
              jbound_max(i)=jcoarse_min+jbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=3
              if(bound_type(i)==2)boundary_type(i)=13
              if(bound_type(i)==3)boundary_type(i)=23
           else
              jbound_min(i)=jcoarse_max+jbound_min(i)
              jbound_max(i)=jcoarse_max+jbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=4
              if(bound_type(i)==2)boundary_type(i)=14
              if(bound_type(i)==3)boundary_type(i)=24
           end if
           if(ndim>2)kbound_min(i)=kcoarse_min+kbound_min(i)
           if(ndim>2)kbound_max(i)=kcoarse_max+kbound_max(i)
        else if(kbound_min(i)*kbound_max(i)==1.and.ndim>2.and.bound_type(i)>0)then
           ibound_min(i)=icoarse_min+ibound_min(i)
           ibound_max(i)=icoarse_max+ibound_max(i)
           jbound_min(i)=jcoarse_min+jbound_min(i)
           jbound_max(i)=jcoarse_max+jbound_max(i)
           if(kbound_min(i)==-1)then
              kbound_min(i)=kcoarse_min+kbound_min(i)
              kbound_max(i)=kcoarse_min+kbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=5
              if(bound_type(i)==2)boundary_type(i)=15
              if(bound_type(i)==3)boundary_type(i)=25
           else
              kbound_min(i)=kcoarse_max+kbound_min(i)
              kbound_max(i)=kcoarse_max+kbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=6
              if(bound_type(i)==2)boundary_type(i)=16
              if(bound_type(i)==3)boundary_type(i)=26
           end if
        end if
     end do
     do i=1,nboundary
        ! Check for errors
        if( (ibound_min(i)<0.or.ibound_max(i)>(nx-1)) .and. (ndim>0) .and.bound_type(i)>0 )then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along X direction',i
           nml_ok=.false.
        end if
        if( (jbound_min(i)<0.or.jbound_max(i)>(ny-1)) .and. (ndim>1) .and.bound_type(i)>0)then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along Y direction',i
           nml_ok=.false.
        end if
        if( (kbound_min(i)<0.or.kbound_max(i)>(nz-1)) .and. (ndim>2) .and.bound_type(i)>0)then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along Z direction',i
           nml_ok=.false.
        end if
     end do
  end if
  nboundary=nboundary_true
  if(simple_boundary.and.nboundary==0)then
     simple_boundary=.false.
  endif

  !--------------------------------------------------
  ! Compute boundary conservative variables
  !--------------------------------------------------
  do i=1,nboundary
     boundary_var(i,1)=MAX(d_bound(i),smallr)
     boundary_var(i,2)=d_bound(i)*u_bound(i)
#if NDIM>1
     boundary_var(i,3)=d_bound(i)*v_bound(i)
#endif
#if NDIM>2
     boundary_var(i,4)=d_bound(i)*w_bound(i)
#endif
     ek_bound=0.0d0
     do idim=1,ndim
        ek_bound=ek_bound+0.5d0*boundary_var(i,idim+1)**2/boundary_var(i,1)
     end do
     boundary_var(i,ndim+2)=ek_bound+P_bound(i)/(gamma-1.0d0)
  end do

  !-----------------------------------
  ! Rearrange level dependent arrays
  !-----------------------------------
  do i=nlevelmax,levelmin,-1
     jeans_refine(i)=jeans_refine(i-levelmin+1)
  end do
  do i=1,levelmin-1
     jeans_refine(i)=-1.0
  end do

  !-----------------------------------
  ! Sort out passive variable indices
  !-----------------------------------
  imetal=ndim+3
  idelay=imetal
  if(metal)idelay=imetal+1
  ixion=idelay
  if(delayed_cooling)ixion=idelay+1
  ichem=ixion
  if(aton)ichem=ixion+1

end subroutine read_hydro_params

