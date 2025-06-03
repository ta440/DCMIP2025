module ic_horiz_mount_flow
    !-----------------------------------------------------------------------
    !
    ! Purpose: Set idealized initial conditions for the horizontal
    !          mountain flow test cases. These both lead to interesting
    !          orographically-driven dynamics.
    !          There are two cases that can be set:
    !            i) A gap flow.
    !            ii) A vortex street
    !   
    !       
    !    
    !        
    !
    !-----------------------------------------------------------------------
    use cam_logfile,         only: iulog
    use shr_kind_mod,        only: r8 => shr_kind_r8
    use cam_abortutils,      only: endrun
    use spmd_utils,          only: masterproc
    use shr_sys_mod,         only: shr_sys_flush
  
    use physconst, only : cpair, gravit, rearth, pi, omega, rair, cappa, &
                          latvap, zvir, tmelt, rh2o
    use hycoef, only : hyai, hybi, hyam, hybm, ps0
  
    implicit none
    private
  
    !=======================================================================
    !    Mountain-induced Rossby wave test case parameters
    !=======================================================================
    real(r8), parameter, private ::             &
         u0                     = 10._r8,       &  ! maximum initial velocity (m/s)
         T0                     = 288._r8,      &  ! Isothermal temperature (K)
         h0                     = 1500._r8,     &  ! peak height of the mountain (m)
         small_fact             = 20._r8,       &  ! Small earth factor
         mount_longitude        = 180._r8,      &  ! Longitude of the mountain (deg)
         p0                     = 1.e5_r8,      &  ! Reference surface pressure (Pa)
         p_sp                   = 1.e5_r8,    &  ! Surface pressure at the south pole (Pa)
         e0_sat                 = 610.78_r8        ! Saturation vapour pressure at the freezing/melting point tmelt
         
    real(r8), parameter :: deg2rad = pi/180._r8    ! conversion to radians
  
    ! Public interface
    public :: horiz_mount_flow_set_ic
  
  contains
  
    subroutine horiz_mount_flow_set_ic(vcoord, latvals, lonvals, U, V, T, PS, PHIS, &
                                       Q, m_cnst, mask, verbose, test_enum)
      use dyn_tests_utils,     only: vc_moist_pressure, vc_dry_pressure, vc_height
      use constituents,        only: cnst_name
      use const_init,          only: cnst_init_default
      use inic_analytic_utils, only: analytic_ic_is_moist
  
      !-----------------------------------------------------------------------
      !
      ! Purpose: Set mountain-induced Rossby wave initial values for dynamics state variables
      !
      !-----------------------------------------------------------------------
  
      ! Dummy arguments
      integer, intent(in)               :: vcoord
      real(r8),           intent(in)    :: latvals(:) ! lat in degrees (ncol)
      real(r8),           intent(in)    :: lonvals(:) ! lon in degrees (ncol)
                                                      ! z_k for vccord 1)
      real(r8), optional, intent(inout) :: U(:,:)     ! zonal velocity
      real(r8), optional, intent(inout) :: V(:,:)     ! meridional velocity
      real(r8), optional, intent(inout) :: T(:,:)     ! temperature
      real(r8), optional, intent(inout) :: PS(:)      ! surface pressure
      real(r8), optional, intent(out)   :: PHIS(:)    ! surface geopotential
      real(r8), optional, intent(inout) :: Q(:,:,:)   ! tracer (ncol, lev, m)
      integer,  optional, intent(in)    :: m_cnst(:)  ! tracer indices (reqd. if Q)
      logical,  optional, intent(in)    :: mask(:)    ! Only init where .true.
      logical,  optional, intent(in)    :: verbose    ! For internal use
      integer, optional, intent(in)     :: test_enum  ! Choice of test case
      ! Local variables
      logical, allocatable              :: mask_use(:)
      logical                           :: verbose_use
      logical                           :: lu,lv,lt,lq,l3d_vars
      integer                           :: i, k, m
      integer                           :: ncol
      integer                           :: nlev
      integer                           :: ncnst
      character(len=*), parameter       :: subname = 'HORIZ_MOUNT_FLOW_SET_IC'
      real(r8)                          :: r(size(latvals))
      real(r8)                          :: mount_lon, mount_lat
      real(r8)                          :: N2, epsilon
      real(r8)                          :: tmp1, tmp2
      real(r8)                          :: p_k(size(latvals))
      real(r8)                          :: tmp(size(latvals))
      real(r8)                          :: surf_h(size(latvals))
      real(r8)                          :: surf_p(size(latvals))
      real(r8)                          :: mount(size(latvals))
      real(r8)                          :: gap(size(latvals))
      real(r8)                          :: qmax
      real(r8)                          :: d, lon_d, lat_d, gap_d
      real(r8)                          :: lon_d_prop, lat_d_prop, gap_d_prop
      real(r8)                          :: lon_expon, lat_expon, gap_expon
      real(r8)                          :: lon_scale, lat_scale, gap_scale
      real(r8)                          :: mount_latitude

      !---------------------------------------
      ! Parameters pertaining to the choice of test
      !---------------------------------------

      if (test_enum == 1) then
        ! Gap flow parameters
        mount_latitude = 0._r8

        ! Length scales of the topography
        lon_d = 400.e3_r8/small_fact
        lat_d = 3000.e3_r8/small_fact
        gap_d = 250.e3_r8/small_fact

        ! Proportion of h0 after moving
        ! lon_d/lat_d/gap_d away from the centre
        ! of the topography
        lon_d_prop = 0.1_r8
        lat_d_prop = 0.1_r8
        gap_d_prop = 0.1_r8

        ! Topographical exponent
        lon_expon = 10._r8
        lat_expon = 10._r8
        gap_expon = 10._r8
      else
        ! Vortex street parameters
        mount_latitude = 20._r8
        d = 250.e3_r8/small_fact
      end if

  
      !--------------------------------------
      ! Helpful parameters
      !-------------------------------------
  
      mount_lon = mount_longitude * deg2rad  ! Longitude of the mountain in rad
      mount_lat = mount_latitude * deg2rad   ! Latitiude of the mountain in rad
      N2        = gravit**2/(cpair*T0)       ! Squared Brunt-Vaisalla frequency (1/s)
      epsilon   = rair/rh2o                  ! Ratio of specific gas constants

      allocate(mask_use(size(latvals)))
      if (present(mask)) then
        if (size(mask_use) /= size(mask)) then
          call endrun(subname//': input, mask, is wrong size')
        end if
        mask_use = mask
      else
        mask_use = .true.
      end if
  
      if (present(verbose)) then
        verbose_use = verbose
      else
        verbose_use = .true.
      end if
  
      ncol = size(latvals, 1)
      nlev = -1
  
      !
      ! We do not yet handle height-based vertical coordinates
      if (vcoord == vc_height) then
        call endrun(subname//':  height-based vertical coordinate not currently supported')
      end if
  
      !
      !******************************
      ! 
      ! Surface geopotential, surface pressure, and specific humidity
      !
      !*****************************

      if (test_enum == 1) then
        ! Compute the scales for gap topography
        lon_scale = lon_d/rearth * (log(1._r8/lon_d_prop) ** (-1._r8/lon_expon))
        lat_scale = lat_d/rearth * (log(1._r8/lat_d_prop) ** (-1._r8/lat_expon))
        gap_scale = gap_d/rearth * (log(1._r8/gap_d_prop) ** (-1._r8/gap_expon))

        ! Construct the gap topography
        where(mask_use)
            mount(:) = exp( - ((lonvals(:) - mount_lon)/lon_scale)**lon_expon - ((latvals(:) - mount_lat)/lat_scale)**lat_expon )
            gap(:) = 1._r8 - exp(-((latvals(:) - mount_lat)/gap_scale)**gap_expon)
            surf_h(:) = h0 * mount(:) * gap(:)
        end where
      else
        ! Construct the Gaussian mountain
        where(mask_use)
            r(:) = rearth * acos(sin(mount_lat)*sin(latvals(:))+cos(mount_lat)*cos(latvals(:))*cos(lonvals(:)-mount_lon))
            surf_h(:) = h0 * exp(-(r(:)/d)**2.0_r8)
        end where
      end if
  
      tmp1 = (rearth * N2 * u0)/(2.0_r8*cappa*(gravit**2.0_r8))*(u0/rearth + 2.0_r8*omega)
      tmp2 = N2/((gravit**2)*cappa)

      where(mask_use)
        ! Hydrostatically-balanced surface pressure
        surf_p(:) = p_sp*exp(-tmp1*((sin(latvals(:)))**2.0_r8 - 1.0_r8) - tmp2*gravit*surf_h(:))
      end where

      if (analytic_ic_is_moist()) then
        qmax = 0.8_r8*epsilon*e0_sat/p0*exp(-latvap/rh2o*(1/T0 - 1/tmelt))
      else
        qmax = 0
      end if
      !
      !*******************************
      !
      ! Initialize PHIS
      !
      !*******************************
      !
      if (present(PHIS)) then
        where(mask_use)
          PHIS(:) = gravit * surf_h(:)
        end where
        if(masterproc .and. verbose_use) then
          write(iulog,*) '          PHIS initialized by "',subname,'"'
        end if
      end if
      !
      !*******************************
      !
      ! Initialize surface pressure
      !
      !*******************************
      !
      if (present(PS)) then
        where(mask_use)
          PS(:) = surf_p(:)
        end where
  
        if(masterproc .and. verbose_use) then
          write(iulog,*) '          PS initialized by "',subname,'"'
        end if
      end if
      !
      !*******************************
      !
      ! Initialize 3D vars
      !
      !*******************************
      !
      lu = present(U)
      lv = present(V)
      lT = present(T)
      lq = present(Q)
      l3d_vars = lu .or. lv .or. lt .or.lq
      nlev = -1
      if (l3d_vars) then
        if (lu) nlev = size(U, 2)
        if (lv) nlev = size(V, 2)
        if (lt) nlev = size(T, 2)
        if (lq) nlev = size(Q, 2)
  
        if (lu) then
          do k = 1, nlev
            where(mask_use)
              U(:,k) = u0 * cos(latvals(:)) 
            end where
          end do
          if(masterproc.and. verbose_use) then
            write(iulog,*) '          U initialized by "',subname,'"'
          end if
        end if
        if (lv) then
          do k = 1, nlev
            where(mask_use)
              V(:,k) = 0.0_r8
            end where
          end do
          if(masterproc.and. verbose_use) then
            write(iulog,*) '          V initialized by "',subname,'"'
          end if
        end if
        if (lt) then
         do k = 1, nlev
            where(mask_use)
              ! Isothermal temperature profile
              p_k(:) = hyam(k)*p0 + hybm(k)*surf_p(:)
              T(:,k) = T0/(1.0_r8 + zvir*qmax*p_k(:))
            end where
         end do
          if(masterproc.and. verbose_use) then
            write(iulog,*) '          T initialized by "',subname,'"'
          end if
        end if
        if (lq) then
          do k = 1, nlev
            where(mask_use)
              Q(:,k,1) = qmax*(hyam(k)*p0 + hybm(k)*surf_p(:))/surf_p(:)
            end where
          end do
          if(masterproc.and. verbose_use) then
            write(iulog,*) '         ', trim(cnst_name(m_cnst(1))), ' initialized by "',subname,'"'
          end if
        end if
      end if
  
      if (lq) then
        ncnst = size(m_cnst, 1)
        if ((vcoord == vc_moist_pressure) .or. (vcoord == vc_dry_pressure)) then
          do m = 2, ncnst
            call cnst_init_default(m_cnst(m), latvals, lonvals, Q(:,:,m_cnst(m)),&
                 mask=mask_use, verbose=verbose_use, notfound=.false.)
          end do
        end if
      end if
  
      deallocate(mask_use)
  
    end subroutine horiz_mount_flow_set_ic
  
  end module ic_horiz_mount_flow
