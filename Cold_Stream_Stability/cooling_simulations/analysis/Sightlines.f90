module parameters
    implicit none
    integer,parameter :: Nmax = 118000000
    real(8),parameter :: pi=3.141592654_8, pi2=2.0_8*pi
    real(8),parameter :: gamma_s = 5.0_8/3.0_8
    real(8),parameter :: KB = 1.38d-16
    real(8),parameter :: mu = 0.59
    real(8),parameter :: mp = 1.67d-24
    real(8),parameter :: kpc = 3.0856d21
    real(8),parameter :: Unit_length = 2.9568d23
    real(8),parameter :: Unit_time   = 2.083720d16
    real(8),parameter :: Unit_dens   = 2.184210526d-28
    real(8),parameter :: Zstream     = 0.0006
    real(8),parameter :: Zback       = 0.002

!!! Arrays for gas
!!! ---------------------------------
    character(len=256),allocatable :: ART_file_name(:)
    real(4),allocatable :: cell_size_gas(:), xgas(:), ygas(:), zgas(:), vxgas(:), vygas(:), vzgas(:), density_gas(:), pressure_gas(:), color_gas(:)
    real(4) :: res
    integer :: Ngas

!!! Stream data
!!! ---------------------------------
    real(4),allocatable :: tsnap(:)

end module parameters
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module read_binary
use parameters
    implicit none
contains
    subroutine data_gas(filename)
        implicit none
        character (*),intent (in) :: filename
        integer :: Ntrack

        open ( 12 , file = filename, form = 'unformatted', convert = 'big_endian')
        cell_size_gas(:) = 0.0_4
        xgas(:) = 0.0_4
        ygas(:) = 0.0_4
        zgas(:) = 0.0_4
        vxgas(:) = 0.0_4
        vygas(:) = 0.0_4
        vzgas(:) = 0.0_4
        density_gas(:) = 0.0_4
        pressure_gas(:) = 0.0_4
        color_gas(:) = 0.0_4
        Ngas=1

        Ntrack = (Nmax - mod(Nmax,10)) / 10
        DO WHILE (Ngas.lt.Nmax)
            if ( Ngas .eq. 1 .or.  mod(Ngas,Ntrack) .eq. 0 ) then
                print*, 'i of Nmax',Ngas,'of',Nmax
            end if
            read (12,end=6) cell_size_gas(Ngas), xgas(Ngas), ygas(Ngas), zgas(Ngas), vxgas(Ngas), vygas(Ngas), vzgas(Ngas),&
            & density_gas(Ngas), pressure_gas(Ngas), color_gas(Ngas)
            Ngas=Ngas+1
        end do
 6   	continue
        close (12)
        if( Ngas.eq.Nmax ) then
            print *, 'DATA GAS - error reading data: too many cells. Change Nmax.'
            print *, 'Ngas=',Ngas
            stop
        end if
        Ngas=Ngas-1	
        print *, 'Ngas', Ngas

        print *, 'xlimits:', minval(xgas(1:Ngas)), maxval(xgas(1:Ngas))
        print *, 'ylimits:', minval(ygas(1:Ngas)), maxval(ygas(1:Ngas))
        print *, 'zlimits:', minval(zgas(1:Ngas)), maxval(zgas(1:Ngas))

    end subroutine data_gas
end module read_binary
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module splitter
!!! Splits and refines cells for interpolation using oct-tree method.
!!! Also contains routines for rotating data to arbitrary frame
use parameters
    implicit none

    real(4) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, vec(3), vec_vel(3)
    real(4),allocatable :: xprime(:), yprime(:), zprime(:), vxprime(:), vyprime(:), vzprime(:), rprime(:)
    real(4) :: saxis1(3), saxis2(3), saxis3(3), theta, phi
    real(4),parameter :: ocx(8) = (/ 1.0_4,  1.0_4, -1.0_4, -1.0_4,  1.0_4,  1.0_4, -1.0_4, -1.0_4 /)
    real(4),parameter :: ocy(8) = (/ 1.0_4, -1.0_4,  1.0_4, -1.0_4,  1.0_4, -1.0_4,  1.0_4, -1.0_4 /)
    real(4),parameter :: ocz(8) = (/ 1.0_4,  1.0_4,  1.0_4,  1.0_4, -1.0_4, -1.0_4, -1.0_4, -1.0_4 /)
contains
!!! CREATE ORTHONORMAL CO-ORDINATE AXES
!!! ----------------------------------------
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

!!! DEALLOCATE PRIMED VALUES AND RE-ALLOCATE ACCORDINGLY
!!! ----------------------------------------
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

!!! DEALLOCATE PRIMED VALUES
!!! ----------------------------------------
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

!!! ROTATE PARTICLE DATA - NO CELL SPLITTING
!!! ----------------------------------------
    subroutine split0(x, y, z, rcm, vx, vy, vz, vcm)
        implicit none
        real(4), intent(in) :: x, y, z, rcm(3)
        real(4), optional, intent(in) :: vx, vy, vz, vcm(3)

        call allocation(1)

        vec(:) = (/ x-rcm(1),y-rcm(2),z-rcm(3) /)               !!! kpc
        xprime(1) = dot_product(vec,saxis1)                     !!! kpc
        yprime(1) = dot_product(vec,saxis2)                     !!! kpc
        zprime(1) = dot_product(vec,saxis3)                     !!! kpc
        rprime(1) = sqrt( xprime(1)**2 + yprime(1)**2 )         !!! kpc

        if(present(vx) .and. present(vy) .and. present(vz) .and. present(vcm) ) then
            vec_vel(:) = (/ vx-vcm(1),vy-vcm(2),vz-vcm(3) /)    !!! km/s
            vxprime(1) = dot_product(vec_vel,saxis1)            !!! km/s
            vyprime(1) = dot_product(vec_vel,saxis2)            !!! km/s
            vzprime(1) = dot_product(vec_vel,saxis3)            !!! km/s
        end if
    end subroutine split0

!!! ROTATE HIGHEST DENSITY GAS DATA - SPLIT CELLS ONCE
!!! --------------------------------------------------
    subroutine split1(x, y, z, rcm, cell_size, vx, vy, vz, vcm)

        implicit none
        real(4), intent(in) :: x, y, z, rcm(3), cell_size
        real(4), optional, intent(in) :: vx, vy, vz, vcm(3)
        integer :: k1

        call allocation(8)

        if(present(vx) .and. present(vy) .and. present(vz) .and. present(vcm) ) then
            do k1=1,8
                x1 = x + ocx(k1)*cell_size/4.0_4                    !!! kpc
                y1 = y + ocy(k1)*cell_size/4.0_4                    !!! kpc
                z1 = z + ocz(k1)*cell_size/4.0_4                    !!! kpc

                vec(:) = (/ x1-rcm(1),y1-rcm(2),z1-rcm(3) /)        !!! kpc
                xprime(k1) = dot_product(vec,saxis1)                !!! kpc
                yprime(k1) = dot_product(vec,saxis2)                !!! kpc
                zprime(k1) = dot_product(vec,saxis3)                !!! kpc
                rprime(k1) = sqrt( xprime(k1)**2 + yprime(k1)**2 )  !!! kpc

                vec_vel(:) = (/ vx-vcm(1),vy-vcm(2),vz-vcm(3) /)    !!! km/s
                vxprime(k1) = dot_product(vec_vel,saxis1)           !!! km/s
                vyprime(k1) = dot_product(vec_vel,saxis2)           !!! km/s
                vzprime(k1) = dot_product(vec_vel,saxis3)           !!! km/s
            end do
        else
            do k1=1,8
                x1 = x + ocx(k1)*cell_size/4.0_4                    !!! kpc
                y1 = y + ocy(k1)*cell_size/4.0_4                    !!! kpc
                z1 = z + ocz(k1)*cell_size/4.0_4                    !!! kpc

                vec(:) = (/ x1-rcm(1),y1-rcm(2),z1-rcm(3) /)        !!! kpc
                xprime(k1) = dot_product(vec,saxis1)                !!! kpc
                yprime(k1) = dot_product(vec,saxis2)                !!! kpc
                zprime(k1) = dot_product(vec,saxis3)                !!! kpc
                rprime(k1) = sqrt( xprime(k1)**2 + yprime(k1)**2 )  !!! kpc
            end do
        end if
    end subroutine split1

!!! ROTATE LEVEL 2 GAS DATA - SPLIT CELLS TWICE
!!! -------------------------------------------
    subroutine split2(x, y, z, rcm, cell_size, vx, vy, vz, vcm)
        implicit none
        real(4), intent(in) :: x, y, z, rcm(3), cell_size
        real(4), optional, intent(in) :: vx, vy, vz, vcm(3)
        integer :: k1,k2

        call allocation(64)

        if(present(vx) .and. present(vy) .and. present(vz) .and. present(vcm) ) then
            do k1=1,8
                x1 = x + ocx(k1)*cell_size/4.0_4                            !!! kpc
                y1 = y + ocy(k1)*cell_size/4.0_4                            !!! kpc
                z1 = z + ocz(k1)*cell_size/4.0_4                            !!! kpc
                do k2=1,8
                    x2 = x1 + ocx(k2)*cell_size/8.0_4                       !!! kpc
                    y2 = y1 + ocy(k2)*cell_size/8.0_4                       !!! kpc
                    z2 = z1 + ocz(k2)*cell_size/8.0_4                       !!! kpc

                    vec(:) = (/ x2-rcm(1),y2-rcm(2),z2-rcm(3) /)            !!! kpc
                    xprime( k2+8*(k1-1) ) = dot_product(vec,saxis1)         !!! kpc
                    yprime( k2+8*(k1-1) ) = dot_product(vec,saxis2)         !!! kpc
                    zprime( k2+8*(k1-1) ) = dot_product(vec,saxis3)         !!! kpc
                    rprime( k2+8*(k1-1) ) = sqrt( xprime( k2+8*(k1-1) )**2 + yprime( k2+8*(k1-1) )**2 )     !!! kpc

                    vec_vel(:) = (/ vx-vcm(1),vy-vcm(2),vz-vcm(3) /)        !!! km/s
                    vxprime( k2+8*(k1-1) ) = dot_product(vec_vel,saxis1)    !!! km/s
                    vyprime( k2+8*(k1-1) ) = dot_product(vec_vel,saxis2)    !!! km/s
                    vzprime( k2+8*(k1-1) ) = dot_product(vec_vel,saxis3)    !!! km/s
                end do
            end do
        else
            do k1=1,8
                x1 = x + ocx(k1)*cell_size/4.0_4                            !!! kpc
                y1 = y + ocy(k1)*cell_size/4.0_4                            !!! kpc
                z1 = z + ocz(k1)*cell_size/4.0_4                            !!! kpc
                do k2=1,8
                    x2 = x1 + ocx(k2)*cell_size/8.0_4                       !!! kpc
                    y2 = y1 + ocy(k2)*cell_size/8.0_4                       !!! kpc
                    z2 = z1 + ocz(k2)*cell_size/8.0_4                       !!! kpc

                    vec(:) = (/ x2-rcm(1),y2-rcm(2),z2-rcm(3) /)            !!! kpc
                    xprime( k2+8*(k1-1) ) = dot_product(vec,saxis1)         !!! kpc
                    yprime( k2+8*(k1-1) ) = dot_product(vec,saxis2)         !!! kpc
                    zprime( k2+8*(k1-1) ) = dot_product(vec,saxis3)         !!! kpc
                    rprime( k2+8*(k1-1) ) = sqrt( xprime( k2+8*(k1-1) )**2 + yprime( k2+8*(k1-1) )**2 )     !!! kpc
                end do
            end do
        end if
    end subroutine split2

!!! ROTATE LEVEL 3 GAS DATA - SPLIT CELLS 3 TIMES
!!! ---------------------------------------------
    subroutine split3(x, y, z, rcm, cell_size, vx, vy, vz, vcm)
        implicit none
        real(4), intent(in) :: x, y, z, rcm(3), cell_size
        real(4), optional, intent(in) :: vx, vy, vz, vcm(3)
        integer :: k1,k2,k3

        call allocation(512)

        if(present(vx) .and. present(vy) .and. present(vz) .and. present(vcm) ) then
            do k1=1,8
                x1 = x + ocx(k1)*cell_size/4.0_4                                            !!! kpc
                y1 = y + ocy(k1)*cell_size/4.0_4                                            !!! kpc
                z1 = z + ocz(k1)*cell_size/4.0_4                                            !!! kpc
                do k2=1,8
                    x2 = x1 + ocx(k2)*cell_size/8.0_4                                       !!! kpc
                    y2 = y1 + ocy(k2)*cell_size/8.0_4                                       !!! kpc
                    z2 = z1 + ocz(k2)*cell_size/8.0_4                                       !!! kpc
                    do k3=1,8
                        x3 = x2 + ocx(k3)*cell_size/16.0_4                                  !!! kpc
                        y3 = y2 + ocy(k3)*cell_size/16.0_4                                  !!! kpc
                        z3 = z2 + ocz(k3)*cell_size/16.0_4                                  !!! kpc

                        vec(:) = (/ x3-rcm(1),y3-rcm(2),z3-rcm(3) /)                        !!! kpc
                        xprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = dot_product(vec,saxis1)       !!! kpc
                        yprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = dot_product(vec,saxis2)       !!! kpc
                        zprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = dot_product(vec,saxis3)       !!! kpc
                        rprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = sqrt( xprime(k3+8*((k2-1)+ 8*(k1-1)))**2 + &
                        & yprime(k3+8*((k2-1)+ 8*(k1-1)))**2 )                              !!! kpc

                        vec_vel(:) = (/ vx-vcm(1),vy-vcm(2),vz-vcm(3) /)                    !!! km/s
                        vxprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = dot_product(vec_vel,saxis1)  !!! km/s
                        vyprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = dot_product(vec_vel,saxis2)  !!! km/s
                        vzprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = dot_product(vec_vel,saxis3)  !!! km/s
                    end do
                end do
            end do
        else
            do k1=1,8
                x1 = x + ocx(k1)*cell_size/4.0_4                                            !!! kpc
                y1 = y + ocy(k1)*cell_size/4.0_4                                            !!! kpc
                z1 = z + ocz(k1)*cell_size/4.0_4                                            !!! kpc
                do k2=1,8
                    x2 = x1 + ocx(k2)*cell_size/8.0_4                                       !!! kpc
                    y2 = y1 + ocy(k2)*cell_size/8.0_4                                       !!! kpc
                    z2 = z1 + ocz(k2)*cell_size/8.0_4                                       !!! kpc
                    do k3=1,8
                        x3 = x2 + ocx(k3)*cell_size/16.0_4                                  !!! kpc
                        y3 = y2 + ocy(k3)*cell_size/16.0_4                                  !!! kpc
                        z3 = z2 + ocz(k3)*cell_size/16.0_4                                  !!! kpc

                        vec(:) = (/ x3-rcm(1),y3-rcm(2),z3-rcm(3) /)                        !!! kpc
                        xprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = dot_product(vec,saxis1)       !!! kpc
                        yprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = dot_product(vec,saxis2)       !!! kpc
                        zprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = dot_product(vec,saxis3)       !!! kpc
                        rprime( k3+8*( (k2-1)+ 8*(k1-1) ) ) = sqrt( xprime(k3+8*((k2-1)+ 8*(k1-1)))**2 + &
                        & yprime(k3+8*((k2-1)+ 8*(k1-1)))**2 )                              !!! kpc
                    end do
                end do
            end do
        end if
    end subroutine split3

!!! ROTATE LEVEL 4 GAS DATA - SPLIT CELLS 4 TIMES
!!! ---------------------------------------------
    subroutine split4(x, y, z, rcm, cell_size, vx, vy, vz, vcm)
        implicit none
        real(4), intent(in) :: x, y, z, rcm(3), cell_size
        real(4), optional, intent(in) :: vx, vy, vz, vcm(3)
        integer :: k1,k2,k3,k4

        call allocation(4096)

        if(present(vx) .and. present(vy) .and. present(vz) .and. present(vcm) ) then
            do k1=1,8
                x1 = x + ocx(k1)*cell_size/4.0_4                                                            !!! kpc
                y1 = y + ocy(k1)*cell_size/4.0_4                                                            !!! kpc
                z1 = z + ocz(k1)*cell_size/4.0_4                                                            !!! kpc
                do k2=1,8
                    x2 = x1 + ocx(k2)*cell_size/8.0_4                                                       !!! kpc
                    y2 = y1 + ocy(k2)*cell_size/8.0_4                                                       !!! kpc
                    z2 = z1 + ocz(k2)*cell_size/8.0_4                                                       !!! kpc
                    do k3=1,8
                        x3 = x2 + ocx(k3)*cell_size/16.0_4                                                  !!! kpc
                        y3 = y2 + ocy(k3)*cell_size/16.0_4                                                  !!! kpc
                        z3 = z2 + ocz(k3)*cell_size/16.0_4                                                  !!! kpc
                        do k4=1,8
                            x4 = x3 + ocx(k4)*cell_size/32.0_4                                              !!! kpc
                            y4 = y3 + ocy(k4)*cell_size/32.0_4                                              !!! kpc
                            z4 = z3 + ocz(k4)*cell_size/32.0_4                                              !!! kpc

                            vec(:) = (/ x4-rcm(1),y4-rcm(2),z4-rcm(3) /)                                    !!! kpc
                            xprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = dot_product(vec,saxis1)     !!! kpc
                            yprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = dot_product(vec,saxis2)     !!! kpc
                            zprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = dot_product(vec,saxis3)     !!! kpc
                            rprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = &
                            & sqrt( xprime(k4+8*((k3-1)+8*((k2-1)+8*(k1-1))))**2 + &
                            & yprime(k4+8*((k3-1)+8*((k2-1)+8*(k1-1))))**2 )                                !!! kpc

                            vec_vel(:) = (/ vx-vcm(1),vy-vcm(2),vz-vcm(3) /)                                !!! km/s
                            vxprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = dot_product(vec_vel,saxis1)!!! km/s
                            vyprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = dot_product(vec_vel,saxis2)!!! km/s
                            vzprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = dot_product(vec_vel,saxis3)!!! km/s
                        end do
                    end do
                end do
            end do
        else
            do k1=1,8
                x1 = x + ocx(k1)*cell_size/4.0_4                                                            !!! kpc
                y1 = y + ocy(k1)*cell_size/4.0_4                                                            !!! kpc
                z1 = z + ocz(k1)*cell_size/4.0_4                                                            !!! kpc
                do k2=1,8
                    x2 = x1 + ocx(k2)*cell_size/8.0_4                                                       !!! kpc
                    y2 = y1 + ocy(k2)*cell_size/8.0_4                                                       !!! kpc
                    z2 = z1 + ocz(k2)*cell_size/8.0_4                                                       !!! kpc
                    do k3=1,8
                        x3 = x2 + ocx(k3)*cell_size/16.0_4                                                  !!! kpc
                        y3 = y2 + ocy(k3)*cell_size/16.0_4                                                  !!! kpc
                        z3 = z2 + ocz(k3)*cell_size/16.0_4                                                  !!! kpc
                        do k4=1,8
                            x4 = x3 + ocx(k4)*cell_size/32.0_4                                              !!! kpc
                            y4 = y3 + ocy(k4)*cell_size/32.0_4                                              !!! kpc
                            z4 = z3 + ocz(k4)*cell_size/32.0_4                                              !!! kpc

                            vec(:) = (/ x4-rcm(1),y4-rcm(2),z4-rcm(3) /)                                    !!! kpc
                            xprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = dot_product(vec,saxis1)     !!! kpc
                            yprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = dot_product(vec,saxis2)     !!! kpc
                            zprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = dot_product(vec,saxis3)     !!! kpc
                            rprime( k4+8*( (k3-1)+ 8*( (k2-1)+ 8*(k1-1) ) ) ) = &
                            & sqrt( xprime(k4+8*((k3-1)+8*((k2-1)+8*(k1-1))))**2 + &
                            & yprime(k4+8*((k3-1)+8*((k2-1)+8*(k1-1))))**2 )                                !!! kpc
                        end do
                    end do
                end do
            end do
        end if
    end subroutine split4
!!! ROTATE LEVEL 5 GAS DATA - SPLIT CELLS 5 TIMES
!!! ---------------------------------------------
    subroutine split5(x, y, z, rcm, cell_size, vx, vy, vz, vcm)
        implicit none
        real(4), intent(in) :: x, y, z, rcm(3), cell_size
        real(4), optional, intent(in) :: vx, vy, vz, vcm(3)
        integer :: k1,k2,k3,k4,k5

        call allocation(32768)

        if(present(vx) .and. present(vy) .and. present(vz) .and. present(vcm) ) then
            do k1=1,8
                x1 = x + ocx(k1)*cell_size/4.0_4                                                                            !!! kpc
                y1 = y + ocy(k1)*cell_size/4.0_4                                                                            !!! kpc
                z1 = z + ocz(k1)*cell_size/4.0_4                                                                            !!! kpc
                do k2=1,8
                    x2 = x1 + ocx(k2)*cell_size/8.0_4                                                                       !!! kpc
                    y2 = y1 + ocy(k2)*cell_size/8.0_4                                                                       !!! kpc
                    z2 = z1 + ocz(k2)*cell_size/8.0_4                                                                       !!! kpc
                    do k3=1,8
                        x3 = x2 + ocx(k3)*cell_size/16.0_4                                                                  !!! kpc
                        y3 = y2 + ocy(k3)*cell_size/16.0_4                                                                  !!! kpc
                        z3 = z2 + ocz(k3)*cell_size/16.0_4                                                                  !!! kpc
                        do k4=1,8
                            x4 = x3 + ocx(k4)*cell_size/32.0_4                                                              !!! kpc
                            y4 = y3 + ocy(k4)*cell_size/32.0_4                                                              !!! kpc
                            z4 = z3 + ocz(k4)*cell_size/32.0_4                                                              !!! kpc
                            do k5=1,8
                                x5 = x4 + ocx(k5)*cell_size/64.0_4                                                          !!! kpc
                                y5 = y4 + ocy(k5)*cell_size/64.0_4                                                          !!! kpc
                                z5 = z4 + ocz(k5)*cell_size/64.0_4                                                          !!! kpc

                                vec(:) = (/ x5-rcm(1),y5-rcm(2),z5-rcm(3) /)                                                !!! kpc
                                xprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec,saxis1)     !!! kpc
                                yprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec,saxis2)     !!! kpc
                                zprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec,saxis3)     !!! kpc
                                rprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = &
                                & sqrt( xprime( k5+ 8*((k4-1)+ 8*((k3-1)+8*((k2-1)+8*(k1-1)))) )**2 + &
                                & yprime( k5+ 8*((k4-1)+ 8*((k3-1)+8*((k2-1)+8*(k1-1)))) )**2 )                             !!! kpc

                                vec_vel(:) = (/ vx-vcm(1),vy-vcm(2),vz-vcm(3) /)                                            !!! km/s
                                vxprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec_vel,saxis1)!!! km/s
                                vyprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec_vel,saxis2)!!! km/s
                                vzprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec_vel,saxis3)!!! km/s
                            end do
                        end do
                    end do
                end do
            end do
        else
            do k1=1,8
                x1 = x + ocx(k1)*cell_size/4.0_4                                                                            !!! kpc
                y1 = y + ocy(k1)*cell_size/4.0_4                                                                            !!! kpc
                z1 = z + ocz(k1)*cell_size/4.0_4                                                                            !!! kpc
                do k2=1,8
                    x2 = x1 + ocx(k2)*cell_size/8.0_4                                                                       !!! kpc
                    y2 = y1 + ocy(k2)*cell_size/8.0_4                                                                       !!! kpc
                    z2 = z1 + ocz(k2)*cell_size/8.0_4                                                                       !!! kpc
                    do k3=1,8
                        x3 = x2 + ocx(k3)*cell_size/16.0_4                                                                  !!! kpc
                        y3 = y2 + ocy(k3)*cell_size/16.0_4                                                                  !!! kpc
                        z3 = z2 + ocz(k3)*cell_size/16.0_4                                                                  !!! kpc
                        do k4=1,8
                            x4 = x3 + ocx(k4)*cell_size/32.0_4                                                              !!! kpc
                            y4 = y3 + ocy(k4)*cell_size/32.0_4                                                              !!! kpc
                            z4 = z3 + ocz(k4)*cell_size/32.0_4                                                              !!! kpc
                            do k5=1,8
                                x5 = x4 + ocx(k5)*cell_size/64.0_4                                                          !!! kpc
                                y5 = y4 + ocy(k5)*cell_size/64.0_4                                                          !!! kpc
                                z5 = z4 + ocz(k5)*cell_size/64.0_4                                                          !!! kpc

                                vec(:) = (/ x5-rcm(1),y5-rcm(2),z5-rcm(3) /)                                                !!! kpc
                                xprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec,saxis1)     !!! kpc
                                yprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec,saxis2)     !!! kpc
                                zprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = dot_product(vec,saxis3)     !!! kpc
                                rprime( k5+ 8*((k4-1)+ 8*((k3-1)+ 8*((k2-1)+ 8*(k1-1) ) ) ) ) = &
                                & sqrt( xprime( k5+ 8*((k4-1)+ 8*((k3-1)+8*((k2-1)+8*(k1-1)))) )**2 + &
                                & yprime( k5+ 8*((k4-1)+ 8*((k3-1)+8*((k2-1)+8*(k1-1)))) )**2 )                             !!! kpc
                            end do
                        end do
                    end do
                end do
            end do
        end if
    end subroutine split5
end module splitter

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module sightlines
use parameters
use read_binary
use splitter
    implicit none
contains
    subroutine get_sightline(snap_num, res_multiply, S1, S2)

        implicit none
        integer,intent(in) :: snap_num, res_multiply
        real(4),intent(in) :: S1(3), S2(3)
        character(len=256) :: filename, format_string
        integer :: i, j, k, m, kp(2), Ntrack, ngrid
        real(4),allocatable :: s_prof(:), vol_prof(:), density_prof(:),  color_prof(:), temp_prof(:), vlos_prof(:)
        real(4) :: Slen, Svec(3), Shat(3)
        real(4) :: rcen(3), vcen(3), split, xmax, xmin, ymax, ymin, zmax, zmin, xgrid, ygrid, zgrid

        print *, 'MAIN ROUTINE'
        print *, 'SNAP NUM=', snap_num

        Svec = S2-S1
        Slen = Sqrt(dot_product(Svec,Svec))
        Shat = Svec / Slen
        res = res_multiply * minval( cell_size_gas(1:Ngas) )
        ngrid = int(Slen/res)
        print *, 'res=',  res,' ngrid=', ngrid

        if( S1(1) .ge. S2(1) ) then
            xmax = S1(1) + res 
            xmin = S2(1) - res
        else
            xmax = S2(1) + res 
            xmin = S1(1) - res
        end if
        if( S1(2) .ge. S2(2) ) then
            ymax = S1(2) + res 
            ymin = S2(2) - res
        else
            ymax = S2(2) + res 
            ymin = S1(2) - res
        end if
        if( S1(3) .ge. S2(3) ) then
            zmax = S1(3) + res 
            zmin = S2(3) - res
        else
            zmax = S2(3) + res 
            zmin = S1(3) - res
        end if

        allocate( s_prof(ngrid), vol_prof(ngrid), density_prof(ngrid),  color_prof(ngrid), temp_prof(ngrid), vlos_prof(ngrid) )
        s_prof       = -0.5_4*Slen + (/ (i,i=1,ngrid) /)*res + 0.5_4*res
        s_prof       = s_prof * Unit_length / kpc
        vol_prof     = 0.0_4
        density_prof = 0.0_4
        color_prof   = 0.0_4
        temp_prof    = 0.0_4
        vlos_prof    = 0.0_4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        rcen(:) = S1
        vcen(:) = 0.0_4
        saxis3(:) = Shat
        call axes(saxis1,saxis2,saxis3)

        Ntrack = (Ngas - mod(Ngas,10)) / 10
        do i=1,Ngas
            if( i .eq. 1 .or.  mod(i,Ntrack) .eq. 0 ) then
                print*, 'i of Ngas',i,'of',Ngas
            end if
            if( ygas(i) .le. ymax + cell_size_gas(i) .and. ygas(i) .ge. ymin - cell_size_gas(i) ) then  ! Start with y because I think we can always choose sightline with y=const because of rotational symmetry in YZ plane
                if( xgas(i) .le. xmax + cell_size_gas(i) .and. xgas(i) .ge. xmin - cell_size_gas(i) ) then  ! Delta x will likely usually be smaller than Delta z
                    if( zgas(i) .le. zmax + cell_size_gas(i) .and. zgas(i) .ge. zmin - cell_size_gas(i) ) then
                        if(cell_size_gas(i)<1.5_4*res) then
                            call split0(xgas(i), ygas(i), zgas(i), rcen(:), vxgas(i), vygas(i), vzgas(i), vcen(:))
                        else if(cell_size_gas(i)<3.0_4*res) then
                            call split1(xgas(i), ygas(i), zgas(i), rcen(:), cell_size_gas(i), vxgas(i), vygas(i), vzgas(i), vcen(:))
                        else if(cell_size_gas(i)<5.0_4*res) then
                            call split2(xgas(i), ygas(i), zgas(i), rcen(:), cell_size_gas(i), vxgas(i), vygas(i), vzgas(i), vcen(:))
                        else if(cell_size_gas(i)<9.0_4*res) then
                            call split3(xgas(i), ygas(i), zgas(i), rcen(:), cell_size_gas(i), vxgas(i), vygas(i), vzgas(i), vcen(:))
                        else if(cell_size_gas(i)<17.0_4*res) then
                            call split4(xgas(i), ygas(i), zgas(i), rcen(:), cell_size_gas(i), vxgas(i), vygas(i), vzgas(i), vcen(:))
                        else
                            call split5(xgas(i), ygas(i), zgas(i), rcen(:), cell_size_gas(i), vxgas(i), vygas(i), vzgas(i), vcen(:))
                        end if
                        k = size(xprime(:))
                        split = log(sngl(k)) / log(8.0_4)
                        do j=1,k
                            if( rprime(j) .le. res ) then   !!! This is binary whether the cell center is within the cylinder or not. Can try to smooth this radially
                                zgrid = zprime(j)/res
                                kp(1) = floor( zgrid )
                                kp(2) = kp(1)+1
                                do m=1,2
                                    if(kp(m) .ge. 1 .and. kp(m) .le. ngrid) then
                                        vol_prof(kp(m))     = vol_prof(kp(m)) + &
                                        & ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( zgrid-kp(3-m) )

                                        density_prof(kp(m)) = density_prof(kp(m)) + density_gas(i) * &
                                        & ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( zgrid-kp(3-m) )

                                        color_prof(kp(m))   = color_prof(kp(m))   + color_gas(i) * density_gas(i) * &
                                        & ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( zgrid-kp(3-m) )

                                        temp_prof(kp(m))    = temp_prof(kp(m))    + pressure_gas(i) * &
                                        & ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( (zgrid-kp(3-m)) )

                                        vlos_prof(kp(m))    = vlos_prof(kp(m))    + vzprime(j) * density_gas(i) * &
                                        & ( ( cell_size_gas(i) / (2.0_4**split) )**3 ) * abs( zgrid-kp(3-m) )
                                    end if
                                end do
                            end if
                        end do
                        call deallocate_primes()
                    end if
                end if
            end if
        end do

        vlos_prof    = ( vlos_prof    / density_prof ) * ( (Unit_length / Unit_time) / 1.d5 )   ! km/s
        color_prof   = ( color_prof   / density_prof ) * ( Zstream - Zback ) + Zback            ! Z absolute units
        temp_prof    = ( temp_prof    / density_prof ) * ( (Unit_length / Unit_time)**2 ) * ( mu * mp / KB )    ! K
        density_prof = ( density_prof / vol_prof )     * Unit_dens                                              ! gr/cm^3

        do i=1,ngrid
            write(30,'(5(1x,ES12.5))') s_prof(i), density_prof(i), temp_prof(i), color_prof(i), vlos_prof(i)
        end do
        call deallocate_primes()
        deallocate( s_prof, vol_prof, density_prof, color_prof, temp_prof, vlos_prof )

    end subroutine get_sightline
end module sightlines

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program main
    use parameters
    use read_binary
    use splitter
    use sightlines
    implicit none
    integer :: i, j, k, Nsnapshot, Nres, Nsnap, NSL
    real(4) :: s_start(3), s_end(3)
    real(4),allocatable :: xstart(:), ystart(:), zstart(:), xend(:), yend(:), zend(:)
    character(len=256) :: filename, format_string, input_arg

    if(iargc().ge.1) then
        call getarg(1,input_arg)
        read(input_arg,*) Nsnap
        if(iargc().ge.2) then
            call getarg(2,input_arg)
            read(input_arg,*) Nres
        else
            Nres = 1
        end if
    else
        Nsnap = 1
        Nres = 1
    end if
    print *, 'Nsnap, Nres'
    print *, Nsnap, Nres

    open(unit=16,file='./output/time.txt')
    read(16,*) Nsnapshot
    print *, Nsnapshot
    allocate( ART_file_name(Nsnapshot), tsnap(Nsnapshot) )
    do k=1,Nsnapshot
        read(16,'(F6.4)') tsnap(k)
        write(ART_file_name(k),'(a,F6.4,a)') './output/ART_format_t',tsnap(k),'.dat'
    end do
    close(unit=16)
    print *, trim(ART_file_name(1))
    print *, trim(ART_file_name(Nsnapshot))

    print *, 'alocating gas arrays'
    allocate( cell_size_gas(Nmax), xgas(Nmax), ygas(Nmax), zgas(Nmax), vxgas(Nmax), vygas(Nmax), vzgas(Nmax), density_gas(Nmax), pressure_gas(Nmax), color_gas(Nmax) )

    print *, ''
    print *, 't=',tsnap(Nsnap)
    call data_gas(trim(ART_file_name(Nsnap)))

    ! s_start = (/ 0.5_4, 0.5_4 + 2.0_4/32.0_4, 0.5_4 - 6.0_4/32.0_4 /) ! Starting position of sightline
    ! s_end   = (/ 0.5_4, 0.5_4 + 2.0_4/32.0_4, 0.5_4 + 6.0_4/32.0_4 /) ! Ending position of sightline

    NSL = 20
    !NSL = 30
    allocate( xstart(NSL), ystart(NSL), zstart(NSL), xend(NSL), yend(NSL), zend(NSL) )

    xstart = (/ 0.20, 0.40, 0.60, 0.80, 0.20, 0.40, 0.60, 0.80, 0.20, 0.40, 0.60, 0.80, 0.20, 0.40, 0.60, 0.80, 0.20, 0.40, 0.60, 0.80 /)
    xend   = (/ 0.20, 0.40, 0.60, 0.80, 0.20, 0.40, 0.60, 0.80, 0.20, 0.40, 0.60, 0.80, 0.20, 0.40, 0.60, 0.80, 0.20, 0.40, 0.60, 0.80 /)
    ystart = (/ 0.40, 0.40, 0.40, 0.40, 0.45, 0.45, 0.45, 0.45, 0.50, 0.50, 0.50, 0.50, 0.55, 0.55, 0.55, 0.55, 0.60, 0.60, 0.60, 0.60 /)
    yend   = (/ 0.40, 0.40, 0.40, 0.40, 0.45, 0.45, 0.45, 0.45, 0.50, 0.50, 0.50, 0.50, 0.55, 0.55, 0.55, 0.55, 0.60, 0.60, 0.60, 0.60 /)
    
    ! xstart = (/ 0.7513, 0.2551, 0.5060, 0.6991, 0.8909, 0.9593, 0.5472, 0.1386, 0.1493, 0.2575, 0.8407, 0.2543, 0.8143, 0.2435, 0.9293, &   ! Matlab generated random numbers from 0 to 1
    !        &   0.3500, 0.1966, 0.2511, 0.6160, 0.4733, 0.3517, 0.8308, 0.5853, 0.5497, 0.9172, 0.2858, 0.7572, 0.7537, 0.3804, 0.5678 /)
    ! xend   = (/ 0.0759, 0.0540, 0.5308, 0.7792, 0.9340, 0.1299, 0.5688, 0.4694, 0.0119, 0.3371, 0.1622, 0.7943, 0.3112, 0.5285, 0.1656, &   ! Matlab generated random numbers from 0 to 1
    !        &   0.6020, 0.2630, 0.6541, 0.6892, 0.7482, 0.4505, 0.0838, 0.2290, 0.9133, 0.1524, 0.8258, 0.5383, 0.9961, 0.0782, 0.4427 /)
    ! ystart = (/ 0.3525, 0.6732, 0.3142, 0.6031, 0.6190, 0.6383, 0.3442, 0.4624, 0.4100, 0.6125, 0.4743, 0.6540, 0.3807, 0.4114, 0.3671, &   ! Matlab generated random numbers from 0.5-6/32 to 0.5+6/32
    !        &   0.3635, 0.6385, 0.5299, 0.5187, 0.3669, 0.6324, 0.5458, 0.4441, 0.5050, 0.4632, 0.3410, 0.4025, 0.3587, 0.3815, 0.4025 /)
    ! yend   = (/ 0.4690, 0.3311, 0.6510, 0.6668, 0.4966, 0.4960, 0.4391, 0.6500, 0.4510, 0.3542, 0.6051, 0.4587, 0.4031, 0.4640, 0.3487, &   ! Matlab generated random numbers from 0.5-6/32 to 0.5+6/32
    !        &   0.3620, 0.6658, 0.6711, 0.5282, 0.3349, 0.4005, 0.4449, 0.6204, 0.3183, 0.3286, 0.3759, 0.5559, 0.5869, 0.5554, 0.4816/)

    zstart(1:NSL) = 0.5_4 - 8.0_4/32.0_4
    zend(1:NSL)   = 0.5_4 + 8.0_4/32.0_4
    
    do i=1,NSL
        s_start = (/ xstart(i), ystart(i), zstart(i) /) ! Starting position of sightline
        s_end   = (/ xend(i),   yend(i),   zend(i)   /) ! Ending position of sightline

        call open_files( Nsnap, Nres, s_start, s_end )
        call get_sightline( Nsnap, Nres, s_start, s_end )
        close(unit=30)
    end do
    deallocate( xstart, ystart, zstart, xend, yend, zend )
    deallocate( cell_size_gas, xgas, ygas, zgas, vxgas, vygas, vzgas, density_gas, pressure_gas, color_gas )
    deallocate( ART_file_name, tsnap )
    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contains
    subroutine open_files(snap_num, Sres, S1, S2)
        implicit none
        integer,intent(in) :: snap_num, Sres
        real(4),intent(in) :: S1(3), S2(3)
        character(len=256) :: dirname

        print *, 'enter open files'

        write(filename,'(a)') 'mkdir -p ./sightlines'
        call system(filename)
        filename = ''

        write(filename,'(a,F6.4,a,F6.4,a,F6.4,a,F6.4,a,F6.4,a,F6.4)') 'mkdir -p ./sightlines/r1_',S1(1),'_',S1(2),'_',S1(3),'_r2_',S2(1),'_',S2(2),'_',S2(3)
        call system(filename)
        filename = ''

        write(dirname,'(a,F6.4,a,F6.4,a,F6.4,a,F6.4,a,F6.4,a,F6.4)') './sightlines/r1_',S1(1),'_',S1(2),'_',S1(3),'_r2_',S2(1),'_',S2(2),'_',S2(3)

        write(filename,'(a,a,I5.5,a,I1,a)') trim(dirname),'/Snap_',snap_num,'_Nres_',Sres,'.txt'
        open(unit=30,file=filename,form='formatted')
        filename = ''

        print *, 'exit open files'
    end subroutine open_files

end program main

