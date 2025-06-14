!< \section arg_table_held_suarez_1994
!! \htmlinclude held_suarez_1994.html
module held_suarez_1994
  !-----------------------------------------------------------------------
  !
  ! Purpose: Implement idealized Held-Suarez forcings
  !    Held, I. M., and M. J. Suarez, 1994: 'A proposal for the
  !    intercomparison of the dynamical cores of atmospheric general
  !    circulation models.'
  !    Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
  !
  !-----------------------------------------------------------------------

  use ccpp_kinds, only: kind_phys

  ! For debugging output
  use cam_logfile,         only: iulog
  use spmd_utils,          only: masterproc

  ! For constructing pressure levels
  use hycoef,    only : hyai, hybi, hyam, hybm


  implicit none
  private
  save

  public :: held_suarez_1994_init
  public :: held_suarez_1994_run

  !!
  !! Model constants, reset in init call
  !!
  real(kind_phys)            :: pref    = 0.0_kind_phys  ! Surface pressure


  !!!!!!!!!!!!!!!!!!!!!!!!
  ! Constants need for the Rayleigh friction routine
  real(kind_phys), parameter     :: sday = 86400._kind_phys
  real(kind_phys), parameter     :: Rd = 287.04_kind_phys              ! gas constant of dry air 
  real(kind_phys), parameter     :: T0 = 288._kind_phys             ! Isothermal temperature
  real(kind_phys), parameter     :: gravit = 9.80616_kind_phys    ! Earth's gravitational accelaration



!=======================================================================
contains
!=======================================================================

!> \section arg_table_held_suarez_1994_init Argument Table
!! \htmlinclude held_suarez_1994_init.html
  subroutine held_suarez_1994_init(pref_in, errmsg, errflg)

    !! Dummy arguments
    real(kind_phys),    intent(in)  :: pref_in

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ' '
    errflg = 0

    pref   = pref_in

  end subroutine held_suarez_1994_init

!> \section arg_table_held_suarez_1994_run Argument Table
!! \htmlinclude held_suarez_1994_run.html
  subroutine held_suarez_1994_run(pver, ncol, pref_mid_norm, clat, cappa, &
       cpair, pmid, ps, uwnd, vwnd, u_reference, v_reference, temp, du, dv, ds, ztodt, scheme_name, errmsg, errflg)

    !
    ! Input arguments
    !
    integer,  intent(in)  :: pver                    ! Num vertical levels
    integer,  intent(in)  :: ncol                    ! Num active columns
    real(kind_phys), intent(in)  :: pref_mid_norm(:) ! reference pressure normalized by surface pressure
    real(kind_phys), intent(in)  :: clat(:)          ! latitudes(radians) for columns
    real(kind_phys), intent(in)  :: cappa(:,:)       ! ratio of dry air gas constant to specific heat at constant pressure
    real(kind_phys), intent(in)  :: cpair(:,:)       ! specific heat of dry air at constant pressure
    real(kind_phys), intent(in)  :: pmid(:,:)        ! mid-point pressure
    real(kind_phys), intent(in)  :: ps(:)            ! Surface pressure
    real(kind_phys), intent(in)  :: uwnd(:,:)        ! zonal wind (m/s)
    real(kind_phys), intent(in)  :: vwnd(:,:)        ! meridional wind (m/s)
    real(kind_phys), intent(in)  :: u_reference(:,:)        ! zonal wind (m/s)
    real(kind_phys), intent(in)  :: v_reference(:,:)        ! meridional wind (m/s)
 
    real(kind_phys), intent(in)  :: temp(:,:)        ! Temperature (K)
    !
    ! Output arguments
    !
    real(kind_phys),   intent(out) :: du(:,:)   ! Zonal wind tend
    real(kind_phys),   intent(out) :: dv(:,:)   ! Meridional wind tend
    real(kind_phys),   intent(out) :: ds(:,:)   ! Heating rate
    real(kind_phys),   intent(in)  :: ztodt     ! Physics timestep
    character(len=64), intent(out) :: scheme_name
    character(len=512),intent(out):: errmsg
    integer,           intent(out):: errflg
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i, k          ! Longitude, level indices

    real(kind_phys)            :: half_pi  ! 0.5*pi
    real(kind_phys)            :: kr, kr_imp, kr_base ! RF friction coefficients
    real(kind_phys)            :: tau      ! Time scale of Rayleigh friction
    real(kind_phys)            :: tau_inv  ! Inverse time scale of Rayleigh friction
    real(kind_phys)            :: H        ! Isothermal scale height
    real(kind_phys)            :: z_i_top(ncol) ! Model interface top, and top two midpoint values in m
    real(kind_phys)            :: p_c(ncol)    ! Cut-off pressure, in Pa
    real(kind_phys)            :: p_i_top(ncol)      ! Model top pressure, in Pa
    real(kind_phys)            :: RF_layer_depth 
    real(kind_phys)            :: z_mid1(ncol), z_mid2(ncol), dz(ncol)
    


    ! User defined parameters

    ! Choose the time scale
    !tau = sday ! 1 day
    tau = sday*5.0_kind_phys

    ! Depth of the RF layer
    RF_layer_depth = 15000._kind_phys


    !
    !-----------------------------------------------------------------------
    !

    errmsg = ' '
    errflg = 0
    scheme_name = "HELD_SUAREZ"

    ! Compute pi/2, noting that arctan(1) = pi/4
    half_pi = 2._kind_phys*atan(1._kind_phys)

    ! Inverse time scale
    tau_inv = 1._kind_phys/tau ! Units of s^-1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute the cut-off pressure at approximately
    ! 10 km below the model top. This assumes an isothermal 
    ! atmosphere with scale height H.
    
    H = Rd*T0/gravit ! Scale height

    z_mid1(:) = H*log(pref/pmid(:,1))
    z_mid2(:) = H*log(pref/pmid(:,2))

    dz(:) = z_mid1(:) - z_mid2(:)
    z_i_top(:) = z_mid1(:) + dz(:)/2

    ! Approximate the model top from Taylor series:
    p_i_top(:) = pmid(:,1) - (dz(:)*pref/(2*H))*exp(-z_mid1(:)/H)

    ! Compute the cut-off pressure as 10 km
    ! below z top. This value can vary due to the
    ! levels deforming with topography
    p_c(:) = pref * exp(-(z_i_top(:) - RF_layer_depth)/H)  

    ! For sanity, check what upper pressure and z values we get.
    if(masterproc) then
      write(iulog,*) 'z top is', z_i_top
      write(iulog,*) 'p top is', p_i_top
      write(iulog,*) 'pref is', pref
      write(iulog,*) 'cut-off pressure is', p_c
    end if


    !
    !-----------------------------------------------------------------------
    !
    ! Hacked Held/Suarez routine
    !
    !----------------------------------------------------------------------
    !
    ! Apply Rayleigh friction in the upper 10 km
    !
    du(:,:) = 0._kind_phys
    dv(:,:) = 0._kind_phys

    do k = 1, pver
      do i = 1, ncol
        if (pmid(i,k) < p_c(i)) then
          kr = tau_inv*( sin(half_pi * log(p_c(i)/pmid(i,k))/log(p_c(i)/p_i_top(i)))**2._kind_phys )
          kr_imp = 1._kind_phys/(1._kind_phys + kr*ztodt) - 1._kind_phys
          kr_base = kr*ztodt/(1._kind_phys + kr*ztodt)

          du(i,k) = uwnd(i,k)*kr_imp + kr_base*u_reference(i,k)
          dv(i,k) = vwnd(i,k)*kr_imp + kr_base*v_reference(i,k) 
        end if
      end do
    end do

  end subroutine held_suarez_1994_run

end module held_suarez_1994
