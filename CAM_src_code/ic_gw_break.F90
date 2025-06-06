
module ic_gw_break
  !-----------------------------------------------------------------------
  !
  ! Purpose: Set idealized initial conditions for the Ullrich, Melvin,
  !          Jablonowski and Staniforth (QJRMS, 2014) baroclinic
  !          instability test.
  !
  !-----------------------------------------------------------------------
  use cam_logfile,         only: iulog
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use cam_abortutils,      only: endrun
  use spmd_utils,          only: masterproc

  use physconst, only : rair, gravit, rearth, pi, omega, epsilo
  use hycoef,    only : hyai, hybi, hyam, hybm, ps0

  implicit none
  private

  real(r8), parameter :: deg2rad = pi/180.0_r8

  !=======================================================================
  !    Baroclinic wave test case parameters
  !=======================================================================
  real(r8), parameter, private :: Mvap = 0.608_r8           ! Ratio of molar mass dry air/water vapor
  real(r8), parameter, private :: psurf_moist = 100000.0_r8 ! Moist surface pressure

  real(r8), parameter, private ::     &
       T0E        = 310.0_r8,         & ! Temperature at equatorial surface (K)
       T0P        = 240.0_r8,         & ! Temperature at polar surface (K)
       B          = 2.0_r8,           & ! Jet half-width parameter
       KK         = 3.0_r8,           & ! Jet width parameter
       lapse      = 0.005_r8,          & ! Lapse rate parameter
       lapse_strat =-0.005d0  ,      &
       strat_start = 20d3      ,      &
       T_weight = 0.5d0  

  real(r8), parameter, private ::     &
       pertu0     = 0.0_r8,           & ! SF Perturbation wind velocity (m/s)
       pertr      = 1.0_r8/6.0_r8,    & ! SF Perturbation radius (Earth radii)
       pertup     = 0.0_r8,           & ! Exp. perturbation wind velocity (m/s)
       pertexpr   = 0.1_r8,           & ! Exp. perturbation radius (Earth radii)
       pertlon    = pi/9.0_r8,        & ! Perturbation longitude
       pertlat    = 2.0_r8*pi/9.0_r8, & ! Perturbation latitude
       pertz      = 15000.0_r8,       & ! Perturbation height cap
       dxepsilon  = 1.0e-5_r8           ! Small value for numerical derivatives

  ! =========================================================================
  ! Mountain parameters
  ! =========================================================================
  real(r8), parameter, private ::  &
       mountain_amplitude = 4000,                &  ! Peak amplitude of the mountain (m)
       mountain_lon_0 = 0.7777777777_r8 * pi,    &  ! Longitudinal center of mountain(rad) 
       mountain_lon_1 = 0.4_r8 * pi,   &            ! in case of two mountains: longitudinal center of second mountain
       mountain_lon_width =  14.0 * deg2rad, &       ! Distance between 10th percentile, ridge mountain
       mountain_lat_0 = pi / 4.0_r8,           &    ! Latitidinal center of mountain (rad) 
       mountain_lat_width = 40.0 * deg2rad,    &    ! for a ridge mountain
       mountain_lat_scale = mountain_lat_width / (2.0_r8*(-log(0.1_r8))**(1.0/6.0) ),  & ! for a ridge  mountain
       mountain_lon_scale = mountain_lon_width / (2.0_r8*(-log(0.1_r8))**(1.0/2.0) ),  & ! for a ridge  mountain
       mountain_halfwidth = 1._r8/8._r8             ! halfwidth of the Gaussian mountain (if used) without Earth's radius=1 (a/8 with a=1) 


  real(r8), parameter, private ::     &
       moistqlat  = 2.0_r8*pi/9.0_r8, & ! Humidity latitudinal width
       moistqp    = 34000.0_r8,       & ! Humidity vertical pressure width
       moistq0    = 0.018_r8            ! Maximum specific humidity

  real(r8), parameter, private ::     &
       eps        = 1.0e-13_r8,       & ! Iteration threshold
       qv_min     = 1.0e-12_r8          ! Min specific humidity value


  integer,  parameter :: deep  = 0   ! Deep (1) or Shallow (0) test case
  integer,  parameter :: pertt = 0   ! 0: exponential, 1: streamfunction
  real(r8), parameter :: bigx  = 1.0 ! Factor for a reduced size earth

  !
  ! Gauss nodes and weights
  !
  integer , parameter :: num_gauss = 10
  real(r8), parameter, dimension(num_gauss), private :: gaussx =(/&
       -0.97390652851717_r8,-0.865063366689_r8,-0.67940956829902_r8,-0.4333953941292_r8,-0.14887433898163_r8,&
       0.14887433898163_r8,0.4333953941292_r8,0.679409568299_r8,0.86506336668898_r8,0.97390652851717_r8/)

  real(r8), parameter, dimension(num_gauss), private :: gaussw =(/&
       0.06667134430869_r8,0.1494513491506_r8,0.219086362516_r8,0.26926671931_r8,0.29552422471475_r8,        &
       0.2955242247148_r8,0.26926671931_r8,0.21908636251598_r8,0.1494513491506_r8,0.0666713443087_r8/)

  ! Public interface
  public :: gw_break_set_ic

contains

  subroutine gw_break_set_ic(vcoord,latvals, lonvals, zint, U, V, W, T, PS, PHIS, &
       Q, m_cnst, mask, verbose)
    use dyn_tests_utils,     only: vc_moist_pressure, vc_dry_pressure, vc_height
    use constituents,        only: cnst_name
    use const_init,          only: cnst_init_default
    use inic_analytic_utils, only: analytic_ic_is_moist

    !-----------------------------------------------------------------------
    !
    ! Purpose: Set baroclinic wave initial values for dynamics state variables
    !
    !-----------------------------------------------------------------------

    ! Dummy arguments
    integer,            intent(in)    :: vcoord     ! vertical coordinate type
    real(r8),           intent(in)    :: latvals(:) ! lat in degrees (ncol)
    real(r8),           intent(in)    :: lonvals(:) ! lon in degrees (ncol)
    real(r8), optional, intent(in)    :: zint(:,:)  ! interface height (ncol,ilev), ordered top to bottom
    real(r8), optional, intent(inout) :: U(:,:)     ! zonal velocity
    real(r8), optional, intent(inout) :: V(:,:)     ! meridional velocity
    real(r8), optional, intent(inout) :: W(:,:)     ! vertical velocity
    real(r8), optional, intent(inout) :: T(:,:)     ! temperature
    real(r8), optional, intent(inout) :: PS(:)      ! surface pressure
    real(r8), optional, intent(out)   :: PHIS(:)    ! surface geopotential
    real(r8), optional, intent(inout) :: Q(:,:,:)   ! tracer (ncol, lev, m)
    integer,  optional, intent(in)    :: m_cnst(:)  ! tracer indices (reqd. if Q)
    logical,  optional, intent(in)    :: mask(:)    ! only init where .true.
    logical,  optional, intent(in)    :: verbose    ! for internal use

    ! Local variables
    logical, allocatable              :: mask_use(:)
    logical                           :: verbose_use
    integer                           :: i, k, m
    integer                           :: ncol
    integer                           :: nlev
    integer                           :: ncnst
    character(len=*), parameter       :: subname = 'BC_WAV_SET_IC'
    real(r8)                          :: ztop,ptop
    real(r8)                          :: uk,vk,Tvk,qk,pk !mid-level state
    real(r8)                          :: psurface
    real(r8), allocatable             :: psurface_temp(:)
    real(r8), allocatable             :: zsurface_temp(:)
    real(r8)                          :: wvp,qdry
    logical                           :: lU, lV, lW, lT, lQ, l3d_vars
    logical                           :: cnst1_is_moisture
    real(r8), allocatable             :: pdry_half(:), pwet_half(:),zdry_half(:),zk(:)
    real(r8), allocatable             :: zmid(:,:)          ! layer midpoint heights for test tracer initialization
    real(r8)                          :: r(size(latvals))   ! great circle distance (unit circle, radius = 1.)
    logical                           :: mountain_gaussian  ! flag for mountain shape

    ! flag that distinguishes between two mountain shapes
    !mountain_gaussian          = .true.        ! flag that either selects a Gaussian mountain (.true.) or the oscillatory mountain (.false.)
    mountain_gaussian          = .false.     ! flag that either selects a Gaussian mountain (.true.) or the oscillatory mountain (.false.)

    if ((vcoord == vc_moist_pressure) .or. (vcoord == vc_dry_pressure)) then
      !
      ! pressure-based vertical coordinate
      !
      ptop = hyai(1) * ps0
      if (ptop > 1.0e5_r8) then
        call endrun(subname//' ERROR: For iterate_z_given_pressure to work ptop must be less than 100hPa')
      end if
      ztop      = iterate_z_given_pressure(ptop,.false.,ptop,0.0_r8,-1000._r8) !Find height of top pressure surface

    else if (vcoord == vc_height) then
       !
       ! height-based vertical coordinate
       !
       if (present(zint)) then
          ztop = zint(1,1)
       else
          call endrun(subname//' ERROR: z-based vertical coordinate requires using optional arg zint')
       end if
    else
      call endrun(subname//' ERROR: vcoord value out of range')
    end if

    allocate(mask_use(size(latvals)))
    if (present(mask)) then
      if (size(mask_use) /= size(mask)) then
        call endrun(subname//': input, mask, is wrong size')
      end if
      mask_use = mask
    else
      mask_use = .true.
    end if


    !==================================
    ! Define the shape of the mountain
    !==================================
    allocate(zsurface_temp(size(latvals, 1)))
    do i=1,size(latvals,1)
        if (mask_use(i)) then
          if (mountain_gaussian) then
    !       Gaussian mountain 
    !       great circle distance without the Earth's radius (unit circle)
            r(i) = acos( sin(mountain_lat_0)*sin(latvals(i)) + cos(mountain_lat_0)*cos(latvals(i))*cos(lonvals(i)-mountain_lon_1))
            zsurface_temp(i) = exp(- (r(i)/mountain_halfwidth)**2._r8 )
          else
    !       single ridge mountain as defined by Hughes and Jablonowski (2022)
            zsurface_temp(i) = exp(-1.0_r8 * (((latvals(i) - mountain_lat_0)/ mountain_lat_scale)**6 + ((lonvals(i) - mountain_lon_1)/ mountain_lon_scale)**2))
    !       Second ridge mountain, commented out
            zsurface_temp(i) = zsurface_temp(i) + exp(-1.0_r8 * (((latvals(i) - mountain_lat_0)/ mountain_lat_scale)**6 + ((lonvals(i) - (mountain_lon_0))/mountain_lon_scale)**2))
          endif
    !     multiply the mountain shape with the peak amplitude
          zsurface_temp(i) = mountain_amplitude * zsurface_temp(i) 
        end if
    end do

    allocate(psurface_temp(size(latvals, 1)))
    do i=1,size(latvals, 1)
        if (mask_use(i)) then
        psurface_temp(i) = moist_pressure_given_z(zsurface_temp(i), &
                                                        latvals(i))
        end if
    end do
    ! ==================================

    if (present(verbose)) then
      verbose_use = verbose
    else
      verbose_use = .true.
    end if

    if(masterproc .and. verbose .and. present(PS)) then
      write(iulog,*) subname, ': Model top (in km) is at z= ',ztop/1000.0_r8
    end if

    ncol = size(latvals, 1)
    nlev = -1
    !
    !*******************************
    !
    ! initialize surface pressure
    !
    !*******************************
    !
    if (present(PS)) then
      if (vcoord == vc_moist_pressure .or. vcoord == vc_height) then
        where(mask_use)
          PS = psurface_temp
        end where
      else if(vcoord == vc_dry_pressure) then
        !
        ! compute dry surface pressure (subtract water vapor in coloumn)
        !
        do i=1,ncol
          if (mask_use(i)) then
            wvp = weight_of_water_vapor_given_z(zsurface_temp(i),latvals(i),ztop)
            ps(i) = psurface_temp(i)-wvp
          end if
        end do
      endif

      if(masterproc .and. verbose_use) then
        write(iulog,*) '          PS initialized by "',subname,'"'
      end if
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
        PHIS = gravit * zsurface_temp
      end where
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          PHIS initialized by "',subname,'"'
      end if
    end if
    !
    !*******************************
    !
    ! Initialize 3D vars
    !
    !
    !*******************************
    !
    lu = present(U)
    lv = present(V)
    lw = present(W)
    lT = present(T)
    lq = present(Q)
    l3d_vars = lu .or. lv .or. lw .or. lt .or. lq
    nlev = -1
    if (l3d_vars) then
      if (lu) nlev = size(U, 2)
      if (lv) nlev = size(V, 2)
      if (lw) nlev = size(W, 2)
      if (lt) nlev = size(T, 2)

      if (lq) then
         nlev = size(Q, 2)
         ! check whether first constituent in Q is water vapor.
         cnst1_is_moisture = m_cnst(1) == 1
         allocate(zmid(size(Q, 1),nlev))         
      end if

      if (lw) then
        do k = 1, nlev
          where(mask_use)
            W(:,k) = 0.0_r8
          end where
        end do
        if(masterproc.and. verbose_use) then
          write(iulog,*) '          W (nonhydrostatic) initialized by "',subname,'"'
        end if
      end if

      allocate(zk(nlev))
      if ((lq.or.lt) .and. (vcoord == vc_dry_pressure)) then
        allocate(pdry_half(nlev+1))
        allocate(pwet_half(nlev+1))
        allocate(zdry_half(nlev+1))
      end if
      do i=1,ncol
        if (mask_use(i)) then
          if (vcoord == vc_moist_pressure) then
            psurface = psurface_temp(i)
            wvp = -99
          else if (vcoord == vc_dry_pressure) then
            !
            ! convert surface pressure to dry
            !
            wvp = weight_of_water_vapor_given_z(zsurface_temp(i),latvals(i),ztop)
            psurface = psurface_temp(i)-wvp
          end if

          if (vcoord == vc_moist_pressure .or. vcoord == vc_dry_pressure) then
             do k=1,nlev
                ! compute pressure levels
                pk = hyam(k)*ps0 + hybm(k)*psurface
                ! find height of pressure surface
                zk(k) = iterate_z_given_pressure(pk,(vcoord == vc_dry_pressure),ptop,latvals(i),ztop)
             end do
          else if (vcoord == vc_height) then
             zk = 0.5_r8*(zint(i,1:nlev) + zint(i,2:nlev+1))
          end if

          if (lq) then
             zmid(i,:) = zk(:)
          end if

          do k=1,nlev
            !
            ! wind components
            !
            if (lu.or.lv) call uv_given_z(zk(k),uk,vk,latvals(i),lonvals(i))
            if (lu) U(i,k)   = uk
            if (lv) V(i,k)   = vk
            !
            ! temperature and moisture for moist vertical coordinates
            !
            if ( (lq .or. lt) .and. &
                 (vcoord==vc_moist_pressure .or. vcoord==vc_height) ) then
              if (analytic_ic_is_moist()) then
                pk = moist_pressure_given_z(zk(k),latvals(i))
                qk = qv_given_moist_pressure(pk,latvals(i))
              else
                qk = 0.d0
              end if
              if (lq .and. cnst1_is_moisture) Q(i,k,1) = qk
              if (lt) then
                tvk    = Tv_given_z(zk(k),latvals(i))
                T(i,k) = tvk / (1.d0 + Mvap * qk)
              end if
            end if
          end do
          !
          ! temperature and moisture for dry-mass vertical coordinates
          !
          if ((lq.or.lt).and. (vcoord==vc_dry_pressure)) then
            !
            ! compute dry pressure vertical coordinate
            !
            pdry_half(1) = hyai(1)*ps0 + hybi(1)*psurface
            pwet_half(1) = pdry_half(1)
            zdry_half(1) = ztop
            do k=2,nlev+1
              pdry_half(k) =  hyai(k)*ps0 + hybi(k)*psurface
              ! find height of pressure surfaces corresponding moist pressure
              zdry_half(k) = iterate_z_given_pressure(pdry_half(k),.true.,ptop,latvals(i),ztop)
              pwet_half(k) = pdry_half(k)+weight_of_water_vapor_given_z(zdry_half(k),latvals(i),ztop)
            end do

            do k=1,nlev
              if (analytic_ic_is_moist()) then
                qdry =((pwet_half(k+1)-pwet_half(k))/(pdry_half(k+1)-pdry_half(k)))-1.0_r8
                qdry = MAX(qdry,qv_min/(1.0_r8-qv_min))
              else
                qdry = 0.0_r8
              end if
              if (lq .and. cnst1_is_moisture) then
                Q(i,k,1) = qdry
              end if
              if (lt) then
                !
                ! convert virtual temperature to temperature
                !
                tvk    = Tv_given_z(zk(k),latvals(i))
                T(i,k) = tvk*(1.0_r8+qdry)/(1.0_r8+(1.0_r8/epsilo)*qdry)
              end if
            end do
          end if
        end if
      end do
      if(lu .and. masterproc.and. verbose_use)  write(iulog,*) '          U initialized by "',subname,'"'
      if(lv .and. masterproc.and. verbose_use)  write(iulog,*) '          V initialized by "',subname,'"'
      if(lt .and. masterproc.and. verbose_use)  write(iulog,*) '          T initialized by "',subname,'"'
      if(lq .and. cnst1_is_moisture .and. masterproc.and. verbose_use)  write(iulog,*) &
           '          ', trim(cnst_name(m_cnst(1))), ' initialized by "',subname,'"'
    end if

    if (lq) then

       ncnst = size(m_cnst, 1)

       do m = 1, ncnst

          ! water vapor already done above
          if (m_cnst(m) == 1) cycle

          call cnst_init_default(m_cnst(m), latvals, lonvals, Q(:,:,m),&
               mask=mask_use, verbose=verbose_use, notfound=.false.,&
               z=zmid)               
          
       end do
    end if   ! lq

    deallocate(mask_use)

    deallocate(zsurface_temp)

    deallocate(psurface_temp)

    if (l3d_vars) then
      deallocate(zk)
      if ((lq.or.lt) .and. (vcoord == vc_dry_pressure)) then
        deallocate(pdry_half)
        deallocate(pwet_half)
        deallocate(zdry_half)
      end if
    end if
  end subroutine gw_break_set_ic

  real(r8) FUNCTION iterate_z_given_pressure(p,ldry_mass_vertical_coordinates,ptop,lat,ztop)

    real(r8), INTENT(IN)  :: &
         p,              &! Pressure (Pa)
         ptop,&! Pressure (Pa)
         lat,&! latitude
         ztop

    logical, INTENT(IN)  :: ldry_mass_vertical_coordinates

    integer :: ix

    real(r8) :: z0, z1, z2
    real(r8) :: p0, p1, p2
    z0 = 0.0_r8
    z1 = 10000.0_r8
    if (ldry_mass_vertical_coordinates) then
       p0 = weight_of_dry_air_given_z(z0,ptop,lat,ztop)
       p1 = weight_of_dry_air_given_z(z1,ptop,lat,ztop)
    else
       p0 =  moist_pressure_given_z(z0,lat)
       p1 =  moist_pressure_given_z(z1,lat)
    endif

    DO ix = 1, 1000
       z2 = z1 - (p1 - p) * (z1 - z0) / (p1 - p0)
       if (ldry_mass_vertical_coordinates) then
          p2 = weight_of_dry_air_given_z(z2,ptop,lat,ztop)
       else
          p2 = moist_pressure_given_z(z2,lat)
       end if

       IF (ABS(p2 - p)/p < eps.or.ABS(z1-z2)<eps.or.ABS(p1-p2)<eps) THEN
          EXIT
       END IF

       z0 = z1
       p0 = p1

       z1 = z2
       p1 = p2
    END DO
    if (ix==1001) then
      write(iulog,*) "p,p1,z1",p,p1,z1
      call endrun('iteration did not converge in iterate_z_given_pressure')
    end if
    iterate_z_given_pressure = z2
  END FUNCTION iterate_z_given_pressure

  real(r8) FUNCTION moist_pressure_given_z(z,lat)

    real(r8), INTENT(IN) :: z,lat
    real(r8)             :: aref
    real(r8)             :: T0, constA, constB, constC, constH, scaledZ, constAA
    real(r8)             :: inttau1, inttau2
    real(r8)             :: rratio, inttermT
    !--------------------------------------------
    ! Constants
    !--------------------------------------------
    aref = rearth / bigX

    T0 = 0.5_r8 * (T0E + T0P)
    constA = 1.0_r8 / lapse
    constAA = 1.d0 / lapse_strat
    constB = (T0 - T0P) / (T0 * T0P)
    constC = 0.5_r8 * (KK + 2.0_r8) * (T0E - T0P) / (T0E * T0P)
    constH = Rair * T0 / gravit

    scaledZ = z / (B * constH)

    !--------------------------------------------
    !    tau values
    !--------------------------------------------
    inttau1 = constB * z * exp(- scaledZ**2)
    inttau2 = constC * z * exp(- scaledZ**2)

    if (z > strat_start) then
        inttau1 = inttau1 + constAA * (exp(lapse_strat * (z - strat_start)/T0)-1.0d0 - (lapse_strat/T0* (z - strat_start)))
        inttau1 = inttau1 + constA * (exp(lapse * strat_start / T0) - 1.0d0) 
        inttau1 = inttau1 + constA * lapse/T0 * exp(lapse * strat_start / T0) * (z - strat_start)
    else
        inttau1 = inttau1 + constA * (exp(lapse * z / T0) - 1.0d0)
    end if

    !--------------------------------------------
    !    radius ratio
    !--------------------------------------------
    if (deep .eq. 0) then
       rratio = 1.0_r8
    else
       rratio = (z + aref) / aref;
    end if

    !--------------------------------------------
    !    interior term on temperature expression
    !--------------------------------------------
    inttermT = (rratio * cos(lat))**KK &
         - KK / (KK + 2.0_r8) * (rratio * cos(lat))**(KK + 2.0_r8)

    !--------------------------------------------
    !    hydrostatic pressure
    !--------------------------------------------
    moist_pressure_given_z = psurf_moist * exp(- gravit / Rair * (inttau1 - inttau2 * inttermT))
  END FUNCTION moist_pressure_given_z

  real(r8) FUNCTION Tv_given_z(z,lat)

    real(r8), INTENT(IN) :: z, lat
    real(r8)             :: aref,T0, constA, constB, constC, constH, scaledZ, constAA
    real(r8)             :: tau1, tau2, rratio, inttermT
    !--------------------------------------------
    ! Constants
    !--------------------------------------------
    aref = rearth / bigX

    T0 = 0.5_r8 * (T0E + T0P)
    constA = 1.0_r8 / lapse
    constAA = 1.d0 / lapse_strat
    constB = (T0 - T0P) / (T0 * T0P)
    constC = 0.5_r8 * (KK + 2.0_r8) * (T0E - T0P) / (T0E * T0P)
    constH = Rair * T0 / gravit

    scaledZ = z / (B * constH)

    !--------------------------------------------
    !    tau values
    !--------------------------------------------
    tau1 = constB * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)
    tau2 = constC * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)

    if (z > strat_start) then
        tau1 = tau1 + constAA * lapse_strat/T0 * (exp(lapse_strat * ((z - strat_start) / T0))-1.d0)
        tau1 = tau1 + constA * lapse / T0 * exp(lapse * strat_start / T0)
    else
        tau1 = tau1 + constA * lapse / T0 * exp(lapse * z / T0) 
    end if

    !--------------------------------------------
    !    radius ratio
    !--------------------------------------------
    if (deep .eq. 0) then
       rratio = 1.0_r8
    else
       rratio = (z + aref) / aref;
    end if

    !--------------------------------------------
    !    interior term on temperature expression
    !--------------------------------------------
    inttermT = (rratio * cos(lat))**KK &
         - KK / (KK + 2.0_r8) * (rratio * cos(lat))**(KK + 2.0_r8)

    !--------------------------------------------
    !    temperature
    !--------------------------------------------
    Tv_given_z = 1.0_r8 / (rratio**2 * (tau1 - tau2 * inttermT))
  END FUNCTION Tv_given_z

  SUBROUTINE uv_given_z(z,u,v,lat,lon)

    real(r8), INTENT(IN)  :: z, lat, lon
    real(r8), INTENT(OUT) :: u,v
    real(r8)              :: aref, omegaref
    real(r8)              :: T0, constH, constC, scaledZ, inttau2, rratio
    real(r8)              :: inttermU, bigU, rcoslat, omegarcoslat
    !------------------------------------------------
    !   Compute test case constants
    !------------------------------------------------
    aref = rearth / bigx
    omegaref = omega * bigx

    T0 = 0.5_r8 * (T0E + T0P)

    constH = Rair * T0 / gravit

    constC = 0.5_r8 * (KK + 2.0_r8) * (T0E - T0P) / (T0E * T0P)

    scaledZ = z / (B * constH)

    inttau2 = constC * z * exp(- scaledZ**2)

    ! radius ratio
    if (deep .eq. 0) then
       rratio = 1.0_r8
    else
       rratio = (z + aref) / aref;
    end if
    !-----------------------------------------------------
    !   Initialize velocity field
    !-----------------------------------------------------
    inttermU = (rratio * cos(lat))**(KK - 1.0_r8) - (rratio * cos(lat))**(KK + 1.0_r8)
    bigU = gravit / aref * KK * inttau2 * inttermU * Tv_given_z(z,lat)
    if (deep .eq. 0) then
       rcoslat = aref * cos(lat)
    else
       rcoslat = (z + aref) * cos(lat)
    end if

    omegarcoslat = omegaref * rcoslat

    u = - omegarcoslat + sqrt(omegarcoslat**2 + rcoslat * bigU)
    v = 0.0_r8

    !-----------------------------------------------------
    !   Add perturbation to the velocity field
    !-----------------------------------------------------

    ! Exponential type
    if (pertt .eq. 0) then
       u = u + evaluate_exponential(z,lat,lon)

       ! Stream function type
    elseif (pertt .eq. 1) then
       u = u - 1.0_r8 / (2.0_r8 * dxepsilon) *                       &
            ( evaluate_streamfunction(lon, lat + dxepsilon, z)    &
            - evaluate_streamfunction(lon, lat - dxepsilon, z))

       v = v + 1.0_r8 / (2.0_r8 * dxepsilon * cos(lat)) *            &
            ( evaluate_streamfunction(lon + dxepsilon, lat, z)    &
            - evaluate_streamfunction(lon - dxepsilon, lat, z))
    end if
  END SUBROUTINE uv_given_z

  !-----------------------------------------------------------------------
  !    Exponential perturbation function
  !-----------------------------------------------------------------------
  real(r8) FUNCTION evaluate_exponential(z,lat,lon)

    real(r8), INTENT(IN)  :: z ! Altitude (meters)
    real(r8), INTENT(IN)  :: lat
    real(r8), INTENT(IN)  :: lon

    real(r8) :: greatcircler, perttaper

    ! Great circle distance
    greatcircler = 1.0_r8 / pertexpr &
         * acos(sin(pertlat) * sin(lat) + cos(pertlat) * cos(lat) * cos(lon - pertlon))

    ! Vertical tapering of stream function
    if (z < pertz) then
       perttaper = 1.0_r8 - 3.0_r8 * z**2 / pertz**2 + 2.0_r8 * z**3 / pertz**3
    else
       perttaper = 0.0_r8
    end if

    ! Zonal velocity perturbation
    if (greatcircler < 1.0_r8) then
       evaluate_exponential = pertup * perttaper * exp(- greatcircler**2)
    else
       evaluate_exponential = 0.0_r8
    end if

  END FUNCTION evaluate_exponential

  !-----------------------------------------------------------------------
  !    Stream function perturbation function
  !-----------------------------------------------------------------------
  real(r8) FUNCTION evaluate_streamfunction(z,lon_local,lat_local)

    real(r8), INTENT(IN)  :: lon_local
    real(r8), INTENT(IN)  :: lat_local
    real(r8), INTENT(IN)  :: z             ! Altitude (meters)

    real(r8) :: greatcircler, perttaper, cospert

    ! Great circle distance
    greatcircler = 1.0_r8 / pertr &
         * acos(sin(pertlat) * sin(lat_local) + cos(pertlat) * cos(lat_local) * cos(lon_local - pertlon))

    ! Vertical tapering of stream function
    if (z < pertz) then
       perttaper = 1.0_r8 - 3.0_r8 * z**2 / pertz**2 + 2.0_r8 * z**3 / pertz**3
    else
       perttaper = 0.0_r8
    end if

    ! Horizontal tapering of stream function
    if (greatcircler < 1.0_r8) then
       cospert = cos(0.5_r8 * pi * greatcircler)
    else
       cospert = 0.0_r8
    end if

    evaluate_streamfunction = &
         (- pertu0 * pertr * perttaper * cospert**4)

  END FUNCTION evaluate_streamfunction

  real(r8) FUNCTION qv_given_moist_pressure(pwet,lat)
    use inic_analytic_utils, only: analytic_ic_is_moist

    real(r8), INTENT(IN)  :: pwet, lat

    real(r8)  :: eta
    if (.not. analytic_ic_is_moist()) then
      qv_given_moist_pressure = 0.0_r8
    else
      eta = pwet/psurf_moist
      if (eta > 0.1_r8) then  ! intialize q if p > 100 hPa
        qv_given_moist_pressure = moistq0 * exp(- (lat/moistqlat)**4)          &
             * exp(- ((eta-1.0_r8)*psurf_moist/moistqp)**2)
      else
        qv_given_moist_pressure = qv_min ! above 100 hPa set q to 1e-12 to avoid supersaturation
      endif
    end if
  END FUNCTION qv_given_moist_pressure

  real(r8) FUNCTION weight_of_water_vapor_given_z(z,lat, ztop)
    use inic_analytic_utils, only: analytic_ic_is_moist

    real(r8), INTENT(IN) :: z,lat, ztop
    real (r8)            :: xm,xr,integral
    real(r8)             :: qv, z1, z2, Tv,pwet, ztmp
    integer              :: jgw

    if (.not. analytic_ic_is_moist()) then
      !
      ! dry case
      !
      weight_of_water_vapor_given_z = 0.0_r8
    else
      z1=z
      z2=ztop
      xm=0.5_r8*(z1+z2)
      xr=0.5_r8*(z2-z1)
      integral=0
      do jgw=1,num_gauss
        ztmp=xm+gaussx(jgw)*xr
        pwet = moist_pressure_given_z(ztmp,lat); qv= qv_given_moist_pressure(pwet,lat);Tv= Tv_given_z(ztmp,lat)
        integral=integral+gaussw(jgw)*gravit*pwet*qv/(Rair*Tv)
      enddo
      integral=0.5_r8*(z2-z1)*integral    ! Scale the answer to the range of integration.
      weight_of_water_vapor_given_z = integral
    end if
  end FUNCTION weight_of_water_vapor_given_z


  real(r8) FUNCTION weight_of_dry_air_given_z(z,ptop,lat,ztop)

    real (r8), INTENT(IN) :: z,ptop, lat, ztop
    real (r8)             :: xm,xr,integral
    real(r8)              :: qv, z1, z2, Tv,pwet, ztmp
    integer               :: jgw

    z1=z
    z2=ztop
    xm=0.5*(z1+z2)
    xr=0.5*(z2-z1)
    integral=0
    do jgw=1,num_gauss
       ztmp=xm+gaussx(jgw)*xr
       pwet = moist_pressure_given_z(ztmp,lat); qv= qv_given_moist_pressure(pwet,lat);Tv= Tv_given_z(ztmp,lat)
       integral=integral+gaussw(jgw)*gravit*pwet*(1-qv)/(Rair*Tv)
    enddo
    integral=0.5_r8*(z2-z1)*integral    ! Scale the answer to the range of integration.
    weight_of_dry_air_given_z = integral+ptop
  end FUNCTION weight_of_dry_air_given_z

end module
