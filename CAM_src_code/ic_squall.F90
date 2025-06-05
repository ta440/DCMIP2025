module ic_squall
  !-----------------------------------------------------------------------
  !
  ! Purpose: Set idealized initial conditions for a squall line test case
  !          based on Klemp et al. JAMES 2015 and Zarzycki et al GMD 2016
  !          supercell test case.
  !
  !-----------------------------------------------------------------------
  use cam_logfile,         only: iulog
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use cam_abortutils,      only: endrun
  use spmd_utils,          only: masterproc

  use physconst, only : rair, gravit, rearth, pi, omega, epsilo, cpair, latvap, rh2o, cappa
  use hycoef,    only : hyai, hybi, hyam, hybm, ps0

  implicit none
  private

  real(r8), parameter :: deg2rad = pi/180.0_r8

  !=======================================================================
  !    Test case parameters
  !=======================================================================
  real(r8), parameter, private :: Mvap = 0.608_r8           ! Ratio of molar mass dry air/water vapor
  real(r8), private            :: psurf_moist  ! Moist surface pressure temp variable


  INTEGER(4), PARAMETER ::            &
       nz         = 40         ,      & ! number of vertical levels in init
       nphi       = 16         ,      & ! number of meridional points in init
       pert_num   =  9                  ! number of thermal bubble perturbations
    
    
  REAL(8), PARAMETER ::               &
       z1         = 0.0d0      ,      & ! lower sample altitude
       z2         = 50000.0d0           ! upper sample altitude
          
  !REAL(8), PARAMETER :: a          = 6371220.0d0 ! reference X = 1 Earth radius
    
  REAL(8), PARAMETER ::               &
       theta0     = 300.d0     ,      & ! theta at the equatorial surface
       theta_tr   = 343.d0     ,      & ! theta at the tropopause
       z_tr       = 12000.d0   ,      & ! altitude at the tropopause
       T_tr       = 213.d0     ,      & ! temperature at the tropopause
       pseq       = 100000.0d0 ,      & ! surface pressure at equator (Pa)
       p0         = 100000.0d0          ! reference pressure for pot. temp and exner
    
  REAL(8), PARAMETER ::               &
       us         = 12.d0      ,      & ! maximum zonal wind velocity
       uc         =  5.d0      ,       & ! coordinate reference velocity
       zs         = 2500.d0    ,      & ! lower altitude of maximum velocity
       zt         = 1000.d0             ! transition distance of velocity
 
  REAL(8), PARAMETER ::               &
       pert_dtheta = 3.d0         ,   & ! perturbation magnitude
       pert_lonc   = 120.d0 * deg2rad,& ! perturbation longitude
       pert_latc   = 0.d0         ,   & ! perturbation latitude
       pert_rh     = 5000.d0     ,    & ! perturbation horiz. halfwidth in m
       pert_zc     = 1500.d0      ,   & ! perturbation center altitude
       pert_rz     = 1500.d0      ,   & ! perturbation vert. halfwidth
       pert_spacing = 10000.d0          ! spacing between perturbations in m
  integer,  parameter :: pert = 1       ! 0: no bubble, 1: bubble included
  
  
  

  !-----------------------------------------------------------------------
  !    Coefficients computed from initialization
  !----------------------------------------------------------------------- 
  INTEGER(4)                  :: initialized = 0

  REAL(8), DIMENSION(nphi)    :: phicoord
  REAL(8), DIMENSION(nz)      :: zcoord
  REAL(8), DIMENSION(nphi,nz) :: thetavyz
  REAL(8), DIMENSION(nphi,nz) :: exneryz
  REAL(8), DIMENSION(nz)      :: qveq


  !=======================================================================
  !    Baroclinic wave test case parameters
  !=======================================================================
  real(r8), parameter, private ::     &
       T0E        = 310.0_r8,         & ! Temperature at equatorial surface (K)
       T0P        = 240.0_r8,         & ! Temperature at polar surface (K)
       B          = 2.0_r8,           & ! Jet half-width parameter
       KK         = 3.0_r8,           & ! Jet width parameter
       lapse      = 0.005_r8            ! Lapse rate parameter

  real(r8), parameter, private ::     &
       pertu0     = 0.5_r8,           & ! SF Perturbation wind velocity (m/s)
       pertr      = 1.0_r8/6.0_r8,    & ! SF Perturbation radius (Earth radii)
       pertup     = 1.0_r8,           & ! Exp. perturbation wind velocity (m/s)
       pertexpr   = 0.1_r8,           & ! Exp. perturbation radius (Earth radii)
       pertlon    = pi/9.0_r8,        & ! Perturbation longitude
       pertlat    = 2.0_r8*pi/9.0_r8, & ! Perturbation latitude
       pertz      = 15000.0_r8,       & ! Perturbation height cap
       dxepsilon  = 1.0e-5_r8           ! Small value for numerical derivatives

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
  public :: squall_set_ic

contains

  subroutine squall_set_ic(vcoord,latvals, lonvals, zint, U, V, W, T, PS, PHIS, &
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

    ! Northern hemisphere latitude
    REAL(8) :: nh_lat

    ! Dummy arguments
    integer,            intent(in)    :: vcoord     ! vertical coordinate type
    real(r8),           intent(in)    :: latvals(:) ! lat in degrees (ncol)
    real(r8),           intent(in)    :: lonvals(:) ! lon in degrees (ncol)
    real(r8), optional, intent(in)    :: zint(:,:)  ! interface height (ncol,ilev), ordered top to bottom
    real(r8), optional, intent(inout) :: U(:,:)     ! zonal velocity
    real(r8), optional, intent(inout) :: V(:,:)     ! meridional velocity
    real(r8), optional, intent(inout) :: W(:,:)     ! vertical velocity   (nonhydrostatic)  
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
    character(len=*), parameter       :: subname = 'SQUALL_WAV_SET_IC'
    real(r8)                          :: ztop,ptop
    real(r8)                          :: zk,uk,vk,Tvk,qk,pk,rhok !mid-level state
    real(r8)                          :: psurface
    real(r8)                          :: wvp,qdry
    logical                           :: lu, lv, lw, lt, lq, l3d_vars
    logical                           :: cnst1_is_moisture
    real(r8), allocatable             :: pdry_half(:), pwet_half(:),zdry_half(:)
    real(r8), allocatable             :: zmid(:,:) ! layer midpoint heights for test tracer initialization
    real(r8)                          :: temp1,temp2,temp3 ! temp variables for initial ztop calculation


    !
    !*******************************
    !
    ! generate initial conditions on Chebyshev nodes
    !
    !*******************************
    !

    call supercell_init()

    
    if ((vcoord == vc_moist_pressure) .or. (vcoord == vc_dry_pressure)) then
      !
      ! pressure-based vertical coordinate
      !
      ptop = hyai(1) * ps0
      if (ptop > 1.0e5_r8) then
        call endrun(subname//' ERROR: For iterate_z_given_pressure to work ptop must be less than 100hPa')
      end if
      ! calculate ztop
      call supercell_p(0.0_r8, 0.0_r8, ptop, ztop,temp1,temp2,temp3,pert)
      
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
      ! sample PS
      
      if (vcoord == vc_moist_pressure .or. vcoord == vc_height) then
        do i=1,ncol
          if (mask_use(i)) then
            CALL supercell_z(lonvals(i), latvals(i), 0.0_r8, psurf_moist, temp1, temp2, temp3, pert)
            PS(i) = psurf_moist
          end if
        end do
      else if(vcoord == vc_dry_pressure) then
        !
        ! compute dry surface pressure (subtract water vapor in coloumn)
        !
        do i=1,ncol
          if (mask_use(i)) then
            CALL supercell_z(lonvals(i), latvals(i), 0.0_r8, psurf_moist, temp1, temp2, temp3, pert)
            wvp = weight_of_water_vapor_given_z(0.0_r8,latvals(i),ztop)
            PS(i) = psurf_moist-wvp
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
      PHIS = 0.0_r8
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
    if (l3d_vars) then ! if at least one of U, V, T, or Q is present
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

      if ((lq.or.lt) .and. (vcoord == vc_dry_pressure)) then
        allocate(pdry_half(nlev+1))
        allocate(pwet_half(nlev+1))
        allocate(zdry_half(nlev+1))
      end if
      do i=1,ncol ! loop over columns (i.e. horizontal)
        if (mask_use(i)) then

          ! calculate surface pressure at col i, which we need to calculate p
          if (vcoord == vc_moist_pressure) then
            CALL supercell_z(lonvals(i), latvals(i), 0.0_r8, psurface, temp1, temp2, temp3, pert)
            wvp = -99
          else if (vcoord == vc_dry_pressure) then
            !
            ! convert surface pressure to dry
            !
            wvp = weight_of_water_vapor_given_z(0.0_r8,latvals(i),ztop)
            CALL supercell_z(lonvals(i), latvals(i), 0.0_r8, psurf_moist, temp1, temp2, temp3, pert)
            psurface = psurf_moist-wvp
          end if

          do k=1,nlev ! loop over levels

            ! compute z or p at midpoints and sample variables
            if (vcoord == vc_moist_pressure .or. vcoord == vc_dry_pressure) then
              ! compute pressure levels
              pk = hyam(k)*ps0 + hybm(k)*psurface

              ! sample z at pk and sample thetav, rho, qv
              CALL supercell_p(lonvals(i), latvals(i), pk, zk, tvk, rhok, qk, pert)                
                
            else if (vcoord == vc_height) then
               zk = 0.5_r8*(zint(i,k) + zint(i,k+1))
               ! sample thetav, rho, and qv at zk
               
               CALL supercell_z(lonvals(i), latvals(i), zk, pk,     &
                               tvk, rhok, qk, pert)
            end if

          
            if (lq) then
               zmid(i,k) = zk
            end if
            !
            ! wind components
            !
            if (lu) U(i,k) = zonal_velocity(zk, latvals(i))
            if (lv) V(i,k) = 0.0_r8
            !
            ! temperature and moisture for moist vertical coordinates
            !
            if ( (lq .or. lt) .and. &
                 (vcoord==vc_moist_pressure .or. vcoord==vc_height) ) then
                 
              if (lq .and. cnst1_is_moisture) Q(i,k,1) = qk
              
              if (lt) then
                T(i,k) = tvk / (1.d0 + Mvap * qk) * (pk / p0)**(rair/cpair)
              end if
            end if
          end do
          !
          ! temperature and moisture for dry-mass vertical coordinates
          !
!          if ((lq.or.lt).and. (vcoord==vc_dry_pressure)) then
!            !
!            ! compute dry pressure vertical coordinate
!            !
!            pdry_half(1) = hyai(1)*ps0 + hybi(1)*psurface
!            pwet_half(1) = pdry_half(1)
!            zdry_half(1) = ztop
!            do k=2,nlev+1
!              pdry_half(k) =  hyai(k)*ps0 + hybi(k)*psurface
!              ! find height of pressure surfaces corresponding moist pressure
!              zdry_half(k) = iterate_z_given_pressure(pdry_half(k),.true.,ptop,latvals(i),ztop)
!              pwet_half(k) = pdry_half(k)+weight_of_water_vapor_given_z(zdry_half(k),latvals(i),ztop)
!            end do
!
!            do k=1,nlev
!              if (analytic_ic_is_moist()) then
!                qdry =((pwet_half(k+1)-pwet_half(k))/(pdry_half(k+1)-pdry_half(k)))-1.0_r8
!                qdry = MAX(qdry,qv_min/(1.0_r8-qv_min))
!              else
!                qdry = 0.0_r8
!              end if
!              if (lq .and. cnst1_is_moisture) then
!                Q(i,k,1) = qdry
!              end if
!              if (lt) then
!                !
!                ! convert virtual temperature to temperature
!                !
!                tvk    = Tv_given_z(zk(k),latvals(i))
!                T(i,k) = tvk*(1.0_r8+qdry)/(1.0_r8+(1.0_r8/epsilo)*qdry)
!              end if
!            end do
!          end if
        end if
      end do
      if(lu .and. masterproc.and. verbose_use)  write(iulog,*) '          U initialized by "',subname,'"'
      if(lv .and. masterproc.and. verbose_use)  write(iulog,*) '          V initialized by "',subname,'"'
      if(lt .and. masterproc.and. verbose_use)  write(iulog,*) '          T initialized by "',subname,'"'
      if(lq .and. cnst1_is_moisture .and. masterproc.and. verbose_use)  write(iulog,*) &
           '          ', trim(cnst_name(m_cnst(1))), ' initialized by "',subname,'"'
    end if
    !
    !*******************************
    !
    ! Initialize W 
    !
    !*******************************
    !    
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


    if (lq) then ! set other consitutients to constant defaults

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
    if (l3d_vars) then
      if ((lq.or.lt) .and. (vcoord == vc_dry_pressure)) then
        deallocate(pdry_half)
        deallocate(pwet_half)
        deallocate(zdry_half)
      end if
    end if
  end subroutine squall_set_ic

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
    real(r8)             :: T0, constA, constB, constC, constH, scaledZ
    real(r8)             :: inttau1, inttau2
    real(r8)             :: rratio, inttermT
    !--------------------------------------------
    ! Constants
    !--------------------------------------------
    aref = rearth / bigX

    T0 = 0.5_r8 * (T0E + T0P)
    constA = 1.0_r8 / lapse
    constB = (T0 - T0P) / (T0 * T0P)
    constC = 0.5_r8 * (KK + 2.0_r8) * (T0E - T0P) / (T0E * T0P)
    constH = rair * T0 / gravit

    scaledZ = z / (B * constH)

    !--------------------------------------------
    !    tau values
    !--------------------------------------------

    inttau1 = constA * (exp(lapse * z / T0) - 1.0_r8) &
         + constB * z * exp(- scaledZ**2)
    inttau2 = constC * z * exp(- scaledZ**2)
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
    moist_pressure_given_z = psurf_moist * exp(- gravit / rair * (inttau1 - inttau2 * inttermT))
  END FUNCTION moist_pressure_given_z

  real(r8) FUNCTION Tv_given_z(z,lat)

    real(r8), INTENT(IN) :: z, lat
    real(r8)             :: aref,T0, constA, constB, constC, constH, scaledZ
    real(r8)             :: tau1, tau2, rratio, inttermT
    !--------------------------------------------
    ! Constants
    !--------------------------------------------
    aref = rearth / bigX

    T0 = 0.5_r8 * (T0E + T0P)
    constA = 1.0_r8 / lapse
    constB = (T0 - T0P) / (T0 * T0P)
    constC = 0.5_r8 * (KK + 2.0_r8) * (T0E - T0P) / (T0E * T0P)
    constH = rair * T0 / gravit

    scaledZ = z / (B * constH)

    !--------------------------------------------
    !    tau values
    !--------------------------------------------
    tau1 = constA * lapse / T0 * exp(lapse * z / T0) &
         + constB * (1.0_r8 - 2.0_r8 * scaledZ**2) * exp(- scaledZ**2)
    tau2 = constC * (1.0_r8 - 2.0_r8 * scaledZ**2) * exp(- scaledZ**2)

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
        integral=integral+gaussw(jgw)*gravit*pwet*qv/(rair*Tv)
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
       pwet = moist_pressure_given_z(ztmp,lat)
       qv= qv_given_moist_pressure(pwet,lat)
       Tv= Tv_given_z(ztmp,lat)
       integral=integral+gaussw(jgw)*gravit*pwet*(1-qv)/(rair*Tv)
    enddo
    integral=0.5_r8*(z2-z1)*integral    ! Scale the answer to the range of integration.
    weight_of_dry_air_given_z = integral+ptop
  end FUNCTION weight_of_dry_air_given_z



!=======================================================================
!    Generate the supercell initial conditions
!=======================================================================
  SUBROUTINE supercell_init() &
    BIND(c, name = "supercell_init")

    IMPLICIT NONE

    ! d/dphi and int(dphi) operators
    REAL(8), DIMENSION(nphi,nphi) :: ddphi, intphi

    ! d/dz and int(dz) operators
    REAL(8), DIMENSION(nz, nz) :: ddz, intz

    ! Buffer matrices for computing SVD of d/dphi operator
    REAL(8), DIMENSION(nphi,nphi) :: ddphibak
    REAL(8), DIMENSION(nphi,nphi) :: svdpu, svdpvt
    REAL(8), DIMENSION(nphi)      :: svdps
    REAL(8), DIMENSION(5*nphi)    :: pwork

    ! Buffer matrices for computing SVD of d/dz operator
    REAL(8), DIMENSION(nz, nz) :: ddzbak
    REAL(8), DIMENSION(nz, nz) :: svdzu, svdzvt
    REAL(8), DIMENSION(nz)     :: svdzs
    REAL(8), DIMENSION(5*nz)   :: zwork

    ! Buffer data for calculation of SVD
    INTEGER(4) :: lwork, info

    ! Sampled values of ueq**2 and d/dz(ueq**2)
    REAL(8), DIMENSION(nphi, nz) :: ueq2, dueq2

    ! Buffer matrices for iteration
    REAL(8), DIMENSION(nphi, nz) :: phicoordmat, dztheta, rhs, irhs
  
    ! Buffer for sampled potential temperature at equator
    REAL(8), DIMENSION(nz) :: thetaeq

    ! Buffer for computed equatorial Exner pressure and relative humidity
    REAL(8), DIMENSION(nz) :: exnereq, H

    ! Variables for calculation of equatorial profile
    REAL(8) :: exnereqs, p, T, qvs, qv

    ! Error metric
    REAL(8) :: err

    ! Loop indices
    INTEGER(4) :: i, k, iter

    ! Chebyshev nodes in the phi direction
    do i = 1, nphi
      phicoord(i) = - cos(dble(i-1) * pi / dble(nphi-1))
      phicoord(i) = 0.25d0 * pi * (phicoord(i) + 1.0d0)
    end do

    ! Matrix of phis
    do k = 1, nz
      phicoordmat(:,k) = phicoord
    end do

    ! Chebyshev nodes in the z direction
    do k = 1, nz
      zcoord(k) = - cos(dble(k-1) * pi / dble(nz-1))
      zcoord(k) = z1 + 0.5d0*(z2-z1)*(zcoord(k)+1.0d0)
    end do

    ! Compute the d/dphi operator
    do i = 1, nphi
      call diff_lagrangian_polynomial_coeffs( &
        nphi, phicoord, ddphi(:,i), phicoord(i))
    end do

    ! Zero derivative at pole
    ddphi(:,nphi) = 0.0d0

    ! Compute the d/dz operator
    do k = 1, nz
      call diff_lagrangian_polynomial_coeffs( &
        nz, zcoord, ddz(:,k), zcoord(k))
    end do

    ! Compute the int(dphi) operator via pseudoinverse
    lwork = 5*nphi

    ddphibak = ddphi
    call DGESVD('A', 'A', &
       nphi, nphi, ddphibak, nphi, &
       svdps, svdpu, nphi, svdpvt, nphi, &
       pwork, lwork, info)

    if (info .ne. 0) then
      write(iulog,*) 'Unable to compute SVD of d/dphi matrix'
      stop
    end if

    do i = 1, nphi
      if (abs(svdps(i)) .le. 1.0d-12) then
        call DSCAL(nphi, 0.0d0, svdpu(1,i), 1)
      else
        call DSCAL(nphi, 1.0d0 / svdps(i), svdpu(1,i), 1)
      end if
    end do
    call DGEMM('T', 'T', &
      nphi, nphi, nphi, 1.0d0, svdpvt, nphi, svdpu, nphi, 0.0d0, &
      intphi, nphi)

    ! Compute the int(dz) operator via pseudoinverse
    lwork = 5*nz

    ddzbak = ddz
    call DGESVD('A', 'A', &
       nz, nz, ddzbak, nz, &
       svdzs, svdzu, nz, svdzvt, nz, &
      zwork, lwork, info)

    if (info .ne. 0) then
      write(iulog,*) 'Unable to compute SVD of d/dz matrix'
      stop
    end if

    do i = 1, nz
      if (abs(svdzs(i)) .le. 1.0d-12) then
        call DSCAL(nz, 0.0d0, svdzu(1,i), 1)
      else
        call DSCAL(nz, 1.0d0 / svdzs(i), svdzu(1,i), 1)
      end if
    end do
    call DGEMM('T', 'T', &
      nz, nz, nz, 1.0d0, svdzvt, nz, svdzu, nz, 0.0d0, &
      intz, nz)

    ! Sample the equatorial velocity field and its derivative
    do k = 1, nz
      ueq2(1,k) = zonal_velocity(zcoord(k), 0.0d0)
      ueq2(1,k) = ueq2(1,k)**2
    end do
    !write(iulog,*) "klev, Ueq2",1,ueq2(1,1)
    do k = 1, nz
      dueq2(1,k) = dot_product(ddz(:,k), ueq2(1,:))
    end do
    !write(iulog,*) "klev, dUeq2",1,dueq2(1,1)
    do i = 2, nphi
      ueq2(i,:) = ueq2(1,:)
      dueq2(i,:) = dueq2(1,:)
    end do

    ! Initialize potential temperature at equator
    do k = 1, nz
      thetaeq(k) = equator_theta(zcoord(k))
      H(k) = equator_relative_humidity(zcoord(k))
    end do
    thetavyz(1,:) = thetaeq
    !write(iulog,*) "klev, thetavyz",1,thetavyz(1,1)


    ! Exner pressure at the equatorial surface
    exnereqs = (pseq / p0)**(rair/cpair)
    !write(iulog,*) "exnereqs",exnereqs
    

    ! Iterate on equatorial profile
    do iter = 1, 12

      ! Calculate Exner pressure in equatorial column (p0 at surface)
      rhs(1,:) = - gravit / cpair / thetavyz(1,:)
      do k = 1, nz
        exnereq(k) = dot_product(intz(:,k), rhs(1,:))
      end do
      do k = 2, nz
        exnereq(k) = exnereq(k) + (exnereqs - exnereq(1))
      end do
      exnereq(1) = exnereqs

      ! Calculate new pressure and temperature
      do k = 1, nz
        p = p0 * exnereq(k)**(cpair/rair)
        T = thetaeq(k) * exnereq(k)

        qvs = saturation_mixing_ratio(p, T)
        qveq(k) = qvs * H(k)

        thetavyz(1,k) = thetaeq(k) * (1.d0 + Mvap * qveq(k))
      end do
    end do

    !do k = 1, nz
    !  write(iulog,*) "exnereq * thetaeq", exnereq(k) * thetaeq(k)
    !end do

    ! Iterate on remainder of domain
    do iter = 1, 12

      ! Compute d/dz(theta)
      do i = 1, nphi
        do k = 1, nz
          dztheta(i,k) = dot_product(ddz(:,k), thetavyz(i,:))
        end do
      end do

      ! Compute rhs
      rhs = sin(2.0d0*phicoordmat)/(2.0d0*gravit) &
            * (ueq2 * dztheta - thetavyz * dueq2)

      ! Integrate
      do k = 1, nz
        do i = 1, nphi
          irhs(i,k) = dot_product(intphi(:,i), rhs(:,k))
        end do
      end do

      ! Apply boundary conditions (fixed Dirichlet condition at equator)
      do i = 2, nphi
        irhs(i,:) = irhs(i,:) + (thetavyz(1,:) - irhs(1,:))
      end do
      irhs(1,:) = thetavyz(1,:)

      ! Compute difference after iteration
      !err = sum(irhs - thetavyz)
      !write(iulog,*) "Balance iter error",iter, err

      ! Update iteration
      thetavyz = irhs
    end do

    

    ! Calculate pressure through remainder of domain
    rhs = - ueq2 * sin(phicoordmat) * cos(phicoordmat) / cpair / thetavyz

    do k = 1, nz
      do i = 1, nphi
        exneryz(i,k) = dot_product(intphi(:,i), rhs(:,k))
      end do
      do i = 2, nphi
        exneryz(i,k) = exneryz(i,k) + (exnereq(k) - exneryz(1,k))
      end do

      exneryz(1,k) = exnereq(k)
    end do

    !do k = 1, nz
    !  write(iulog,*) "exneryz at nphi=1,k=",k,exneryz(1,k)
    !end do

    ! Initialization successful
    initialized = 1
    !write(iulog, *) "Intialization of squall line conditions (1 True, 0 False):", initialized

  END SUBROUTINE supercell_init

!-----------------------------------------------------------------------
!    Calculate pointwise pressure and temperature
!-----------------------------------------------------------------------
  SUBROUTINE supercell_z(lon, lat, z, p, thetav, rho, q, pert)

    REAL(8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (m)

    INTEGER, INTENT(IN) :: pert  ! 1 if perturbation should be included
                                 ! 0 if no perturbation should be included

    ! Evaluated variables
    REAL(8), INTENT(OUT) :: p, thetav, rho, q

    ! Northern hemisphere latitude
    REAL(8) :: nh_lat

    ! Angular measure of pert_spacing
    REAL(8) :: angular_pert_spacing

    ! Pointwise Exner pressure
    REAL(8) :: exner

    ! Assembled variable values in a column
    REAL(8), DIMENSION(nz) :: varcol

    ! Coefficients for computing a polynomial fit in each coordinate
    REAL(8), DIMENSION(nphi) :: fitphi
    REAL(8), DIMENSION(nz)   :: fitz

    ! Loop indices
    INTEGER(4) :: k

    ! Northern hemisphere latitude
    if (lat .le. 0.0d0) then
      nh_lat = -lat
    else
      nh_lat = lat
    end if

    ! Perform fit
    CALL lagrangian_polynomial_coeffs(nz, zcoord, fitz, z)
    CALL lagrangian_polynomial_coeffs(nphi, phicoord, fitphi, nh_lat)

    ! Obtain exner pressure of background state
    do k = 1, nz
      varcol(k) = dot_product(fitphi, exneryz(:,k))
    end do
    exner = dot_product(fitz, varcol)
    p = p0 * exner**(cpair/rair)

    ! Sample the initialized fit at this point for theta_v
    do k = 1, nz
      varcol(k) = dot_product(fitphi, thetavyz(:,k))
    end do
    thetav = dot_product(fitz, varcol)

    ! Sample water vapor mixing ratio
    q = dot_product(fitz, qveq)

    ! Fixed density
    rho = p / (rair * exner * thetav)
    
    if (pert .ne. 0) then
        ! place pert_num thermal perturbations centered at pert_lonc, pert_latc
        
        ! spaced apart by pert_spacing 
        angular_pert_spacing = pert_spacing / rearth
        
        do k = 1, pert_num, 1
            thetav = thetav &
               + thermal_perturbation(lon, lat, z,pert_lonc,pert_latc+(k & 
               - REAL(pert_num+1,8)/2.d0)*angular_pert_spacing) &
               * (1.d0 + Mvap * q)
        end do

    end if

    !write(iulog, *) "p before pert update: ", p
    ! Updated pressure
    p = p0 * (rho * rair * thetav / p0)**(cpair/(cpair-rair))
    !write(iulog, *) "p after pert update: ", p

  END SUBROUTINE supercell_z

!-----------------------------------------------------------------------
!    Calculate pointwise z and temperature given pressure
!-----------------------------------------------------------------------
  SUBROUTINE supercell_p(lon, lat, p, z, thetav, rho, q, pert)

    REAL(8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                p             ! Pressure (Pa)

    INTEGER, INTENT(IN) :: pert  ! 1 if perturbation should be included
                                 ! 0 if no perturbation should be included

    ! Evaluated variables
    REAL(8), INTENT(OUT) :: z, thetav, rho, q

    ! Bounding interval and sampled values
    REAL(8) :: za, zb, zc, pa, pb, pc

    ! Iterate
    INTEGER(4) :: iter

    za = z1  !sfc
    zb = z2  !aloft

    ! works only for moist pressure coordinats
    CALL supercell_z(lon, lat, za, pa, thetav, rho, q, pert)
    CALL supercell_z(lon, lat, zb, pb, thetav, rho, q, pert)

    !write(iulog,*), pa
    !write(iulog,*), pb

    if (pa .lt. p) then
      write(iulog,*) pa, p
      call endrun('Requested pressure out of range on top, adjust sample interval')
    end if
    if (pb .gt. p) then
      write(iulog,*) pb, p
      call endrun('Requested pressure out of range on top, adjust sample interval')
    end if

    ! Iterate using fixed point method
    do iter = 1, 10000

      zc = (za * (pb - p) - zb * (pa - p)) / (pb - pa)
      
      !if (iter .lt. 10) then
      !    write(iulog,*) "za,zb,pa,pb",za,zb,pa,pb
      !end if

      CALL supercell_z(lon, lat, zc, pc, thetav, rho, q, pert)

      !write(*,*) pc
      !write(*,*) p
      !write(*,*) abs((pc - p) / p)

!      if (abs((pc - p) / p) .lt. 1.d-12) then
!        !write(*,*), iter
!        exit
!      end if

      IF (ABS(pc - p)/p < eps.or.ABS(zb-zc)<eps.or.ABS(pb-pc)<eps) THEN
          EXIT
      END IF

      if (pc .gt. p) then
        za = zc
        pa = pc
      else
        zb = zc
        pb = pc
      end if
    end do

    if (iter .eq. 10001) then
      write(iulog,*) "p,pc,zc",p,pc,zc
      write(iulog,*) "za,zb,pa,pb",za,zb,pa,pb
      call endrun('Iteration failed to converge in supercell_p')
    end if

    z = zc

  END SUBROUTINE supercell_p

!-----------------------------------------------------------------------
!    Calculate pointwise z and temperature given pressure
!-----------------------------------------------------------------------
  REAL(8) FUNCTION thermal_perturbation(lon, lat, z, pert_lon, pert_lat)

    REAL(8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z,          & ! Altitude (m)
                pert_lon,   & ! Perturbation Longitude (radians)
                pert_lat      ! Perturbation Latitude (radians)
                
    
    ! Great circle radius from the perturbation centerpoint
    REAL(8) :: gr

    ! Approximately spherical radius from the perturbation centerpoint
    REAL(8) :: Rtheta

    gr = rearth*acos(sin(pert_lat)*sin(lat) + &
         (cos(pert_lat)*cos(lat)*cos(lon-pert_lon)))

    Rtheta = sqrt((gr/pert_rh)**2 + ((z - pert_zc) / pert_rz)**2)

    if (Rtheta .le. 1.d0) then
      thermal_perturbation = pert_dtheta * (cos(0.5d0 * pi * Rtheta))**2
    else
      thermal_perturbation = 0.0d0
    end if

  END FUNCTION thermal_perturbation

!-----------------------------------------------------------------------
!    Calculate the reference zonal velocity
!-----------------------------------------------------------------------
  REAL(8) FUNCTION zonal_velocity(z, lat)

    IMPLICIT NONE

    REAL(8), INTENT(IN) :: z, lat 

    ! coefficients for smooth transition polynomial
    REAL(8) ::  A,  &  ! z^2 coefficient
                B,  &  ! z coefficient
                C      ! constant coefficient
    
    ! Define parameters such that the linear shear layer
    ! and constant wind layer are connected by a smooth,
    ! continuous quadtratic.
    A = -zs / 4.0d0 / zt
    B = (1 + zs / zt) / 2.0d0
    C = 1 / 2.0d0 - zs / zt / 4.0d0 - zt / zs / 4.0d0
    
    
    if (z .le. zs - zt) then
      zonal_velocity = us * (z / zs) - uc
    elseif (abs(z - zs) .le. zt) then
      zonal_velocity = &
        (C + B * z/zs + A * (z**2)/(zs**2)) * us - uc
    else
      zonal_velocity = us - uc
    end if

    zonal_velocity = zonal_velocity * cos(lat)
      
  END FUNCTION zonal_velocity

!-----------------------------------------------------------------------
!    Calculate pointwise theta at the equator at the given altitude
!-----------------------------------------------------------------------
  REAL(8) FUNCTION equator_theta(z)

    IMPLICIT NONE

    REAL(8), INTENT(IN) :: z

    if (z .le. z_tr) then
      equator_theta = &
        theta0 + (theta_tr - theta0) * (z / z_tr)**(1.25d0)
    else
      equator_theta = &
        theta_tr * exp(gravit/cpair/T_tr * (z - z_tr))
    end if

  END FUNCTION equator_theta

!-----------------------------------------------------------------------
!    Calculate pointwise relative humidity (in %) at the equator at the
!    given altitude
!-----------------------------------------------------------------------
  REAL(8) FUNCTION equator_relative_humidity(z)

    IMPLICIT NONE

    REAL(8), INTENT(IN) :: z

    if (z .le. z_tr) then
      equator_relative_humidity = 1.0d0 - 0.75d0 * (z / z_tr)**(1.25d0)
    else
      equator_relative_humidity = 0.25d0
    end if

  END FUNCTION equator_relative_humidity

!-----------------------------------------------------------------------
!    Calculate saturation mixing ratio (in kg/kg) in terms of pressure
!    (in Pa) and temperature (in K)
!-----------------------------------------------------------------------
  REAL(8) FUNCTION saturation_mixing_ratio(p, T)

    IMPLICIT NONE

    REAL(8), INTENT(IN)  :: &
                p,        & ! Pressure in Pa
                T           ! Temperature

    saturation_mixing_ratio = &
      380.d0 / p * exp(17.27d0 * (T - 273.d0) / (T - 36.d0))

    if (saturation_mixing_ratio > 0.014) then
      saturation_mixing_ratio = 0.014
    end if

  END FUNCTION saturation_mixing_ratio

!-----------------------------------------------------------------------
!    Calculate coefficients for a Lagrangian polynomial
!-----------------------------------------------------------------------
  SUBROUTINE lagrangian_polynomial_coeffs(npts, x, coeffs, xs)

    IMPLICIT NONE

    ! Number of points to fit
    INTEGER(4), INTENT(IN) :: npts

    ! Sample points to fit
    REAL(8), DIMENSION(npts), INTENT(IN) :: x

    ! Computed coefficients
    REAL(8), DIMENSION(npts), INTENT(OUT) :: coeffs

    ! Point at which sample is taken
    REAL(8), INTENT(IN) :: xs

    ! Loop indices
    INTEGER(4) :: i, j
    
    ! Compute the Lagrangian polynomial coefficients
    do i = 1, npts
      coeffs(i) = 1.0d0
      do j = 1, npts
        if (i .eq. j) then
          cycle
        end if
        coeffs(i) = coeffs(i) * (xs - x(j)) / (x(i) - x(j))
      end do
    end do

  END SUBROUTINE lagrangian_polynomial_coeffs

!-----------------------------------------------------------------------
!    Calculate coefficients of the derivative of a Lagrangian polynomial
!-----------------------------------------------------------------------
  SUBROUTINE diff_lagrangian_polynomial_coeffs(npts, x, coeffs, xs)

    IMPLICIT NONE

    ! Number of points to fit
    INTEGER(4), INTENT(IN) :: npts

    ! Sample points to fit
    REAL(8), DIMENSION(npts), INTENT(IN) :: x

    ! Computed coefficients
    REAL(8), DIMENSION(npts), INTENT(OUT) :: coeffs

    ! Point at which sample is taken
    REAL(8), INTENT(IN) :: xs

    ! Loop indices
    INTEGER(4) :: i, j, imatch

    ! Buffer sum
    REAL(8) :: coeffsum, differential

    ! Check if xs is equivalent to one of the values of x
    imatch = (-1)
    do i = 1, npts
      if (abs(xs - x(i)) < 1.0d-14) then
        imatch = i
        exit
      end if
    end do

    ! Equivalence detected; special treatment required
    if (imatch .ne. (-1)) then
      do i = 1, npts
        coeffs(i) = 1.0d0
        coeffsum = 0.0d0

        do j = 1, npts
          if ((j .eq. i) .or. (j .eq. imatch)) then
            cycle
          end if

          coeffs(i) = coeffs(i) * (xs - x(j)) / (x(i) - x(j))
          coeffsum = coeffsum + 1.0 / (xs - x(j))
        end do

        if (i .ne. imatch) then
          coeffs(i) = coeffs(i)                   &
            * (1.0 + (xs - x(imatch)) * coeffsum) &
            / (x(i) - x(imatch))
        else
          coeffs(i) = coeffs(i) * coeffsum
        end if
      end do

    ! No equivalence; simply differentiate Lagrangian fit
    else
      call lagrangian_polynomial_coeffs(npts, x, coeffs, xs)

      do i = 1, npts
        differential = 0.0d0
        do j = 1, npts
          if (i .eq. j) then
            cycle
          end if
          differential = differential + 1.0 / (xs - x(j))
        end do
        coeffs(i) = coeffs(i) * differential
      end do
    end if

  END SUBROUTINE


end module ic_squall
