
 A summary of the key changes to enable the nonhydrostatic version of FV3 in CAM.

!!!!!!!!!!!!
1. Modify [CAM_clone]/src/dynamics/fv3/dyn_comp.F90

  a: Declare three new variables, lines 133-144, below 
  ! fv_core namelist variables - these namelist variables defined in fv3 library without fv3_

  Add two variables in the real(r8) group:
  fv3_a_imp, fv3_p_fac,

  and one in the logical group:
  fv3_use_logp.

  b: Lines 175-186, starting with 'namelist /fv_core_nml/ ' add the three new variables: 

  fv3_a_imp, fv3_p_fac, fv3_use_logp.

  c: Switch off the command to break the code if using the nonhydrostatic version.
  At line 240: comment out the following code:
  ! Non-hydrostatic runs not currently supported
  if (.not.fv3_hydrostatic) &
       call endrun('dyn_readnl: ERROR FV3 Non-hydrostatic option is not supported, set namelist fv3_hydrostatic = .true.')
    
  d: Lines 246-308, with comment '! write fv3 dycore namelist options to log', add lines
  for the new variables:
  
    write (iulog,*) '  fv3_a_imp                 = ',fv3_a_imp
    write (iulog,*) '  fv3_p_fac                 = ',fv3_p_fac
    write (iulog,*) '  fv3_use_logp              = ',fv3_use_logp 

 e: Add a statement to initialise W for a nonhydrostatic model. After line 1089, where the 
 V initial state has been set, add the following lines:
 
     dbuf3=0._r8
      if (.not. Atm(mytile)%flagstruct%hydrostatic) then
          call analytic_ic_set_ic(vcoord, latvals_rad, lonvals_rad, glob_ind,            &
                W=dbuf3(:,:,:))
          n=0
          do j = js, je
             do i = is, ie
                ! W. This can be nonzero for nonhydrostatic models
                n=n+1
                atm(mytile)%w(i,j,:) = dbuf3(n, :, 1)
             end do
          end do
      end if

  f: Add a statement for W if an initial condition is read in? 
    First, check if things run without this.

  g: Call a halo exchange for W.
  In the if statement code starting at line 1387:
  !!  Initialize non hydrostatic variables if needed
  if (.not. Atm(mytile)%flagstruct%hydrostatic) then
  
  After the three 'end do' statements and just before the 'end if',
  insert the following line:
  call mpp_update_domains( atm(mytile)%w,    Atm(mytile)%domain )

!!!!!!!!!!!!!!!!!!!!
2. Add information on the new variables into
  [CAM_clone]/bld/namelist_files/namelist_definition.xml

  Somewhere after line 8890
  <!-- for fv3 dynamical core -->
  
  Add the following code sections:

  a)

<entry id="fv3_a_imp" type="real" category="dyn_fv3"
       group="fv_core_nml" valid_values="" >
Real: Controls behavior of the non-hydrostatic solver. Values greater
than 0.5 enable the semi-implicit solver, in which the value of a_imp
controls the time-off-centering: use a_imp = 1.0 for a fully backward
time-stepping. For consistency, the sum of beta and a_imp should
be 1 when the semi-implicit solver is used. The semi-implicit algo-
rithm is substantially more efficient except at very high (km-scale)
resolutions with an acoustic time step of a few seconds or less. 0.75
by default. Proper values are 0, or between 0.5 and 1. This variable is
only used if hydrostatic =.false.
Default: 1.
</entry>

  b)

<entry id="fv3_p_fac" type="real" category="dyn_fv3"
       group="fv_core_nml" valid_values="" >
Real: Safety factor for minimum nonhydrostatic pressures, which
will be limited so the full pressure is no less than p_fac times the
hydrostatic pressure. This is only of concern in mid-top or high-
top models with very low pressures near the model top, and has
no effect in most simulations. The pressure limiting activates only
when model is in danger of blowup due to unphysical negative
total pressures. 0.05 by default. Only used if hydrostatic =.false.
and the semi-implicit solver is used. Proper range is 0 to 0.25.
Default: 0.05
</entry>

  c) 

<entry id="fv3_use_logp" type="real" category="dyn_fv3"
     group="fv_core_nml" valid_values="" >
Logical: Enables a variant of the Lin pressure-gradient force 
algorithm, which uses the logarithm of pressure instead of the Exner
function (as in Lin 1997). This yields more accurate results for re-
gions that are nearly isothermal. Ignored if hydrostatic = true.
Default: FALSE
</entry>

!!!!!!!!!!!!!!!!!
3. In [CAM_clone]/bld/namelist_files/namelist_defaults_cam.xml, add
defaults for the three new variables:

After line 3213
<!-- ================================================================== -->
<!-- Defaults for FV3                                                   -->
<!-- ================================================================== -->

Add new lines lines before/amidst
<fv3_adjust_dry_mass                                                        > .false.   </fv3_adjust_dry_mass>

New lines:
<fv3_a_imp                                                                  > 1.0D0     </fv3_a_imp>
<fv3_p_fac                                                                  > 0.05D0    </fv3_p_fac>
<fv3_use_logp                                                               > .false.   </fv3_use_logp>

!!!!!!!!!!!!!!!!!!
4. Modify [CAM_clone]/bld/build-namelist

Around line 3981:
# FV dycore
if ($dyn eq 'fv') {

Add three new lines:

add_default($nl, 'fv3_a_imp');
add_default($nl, 'fv3_p_fac');
add_default($nl, 'fv3_use_logp');

!!!!!!!!!!!!!!!!
5. Fix a bug in a nonhydrostatic FV3 routine.

Copy the file into the source override directory:

cp [CAM_clone]/src/dynamics/fv3/atmos_cubed_sphere/model/nh_utils.F90 [CAM_clone]/src/dynamics/fv3/src_override/nh_utils.F90

Edit the new nh_utils.F90 file. 

Around line 1258
! Start the w-solver
    do k=2, km
       do i=is, ie

Replace the current double-do loop with:

! Start the w-solver
    do k=2, km
       do i=is, ie
#ifdef MOIST_CAPPA
          aa(i,k) = t1g*0.5*(gm2(i,k-1)+gm2(i,k))/(dz2(i,k-1)+dz2(i,k)) * pem(i,k)
#else
          aa(i,k) = t1g/(dz2(i,k-1)+dz2(i,k)) * pem(i,k)
#endif
       enddo

!!!!!!!!!!!!!!!!!!!!
6. Add W as a variable for the initial conditions in 
[CAM_clone]/src/dynamics/tests/inic_analytic.F90

  a) At line 41, add W as a variable in subroutine dyn_set_inic_col()
and in the ! Dummy arguments list, add:

real(r8), optional, intent(inout) :: W(:,:)      ! vertical velocity (nonhydrostatic)

b) Add W as a variable to subroutine dyn_set_inic_cblock(),
and in the ! Dummy arguments list, add:

real(r8), optional, intent(inout) :: W(:,:)      ! vertical velocity (nonhydrostatic)

c) Add W as an argument to all the initial condition calls,
starting at line 161:
select case(trim(analytic_ic_type))

For example:
case('held_suarez_1994')
  call hs94_set_ic(latvals, lonvals, U=U, V=V, W=W, T=T, PS=PS, PHIS=PHIS_OUT,  &
       Q=Q, m_cnst=m_cnst, mask=mask_use, verbose=verbose_use)


d) Make calls for W as per the other prognostics.
This requires a lot of changes (8, by my count).
Anywhere there is a code with if (present(U)) etc. 
a corresponding block for W needs to be added.
For example,

if(present(W)) then
  call get_input_shape(W, 'W', mname, size1, size2, size3, subname)
end if

!!!!!!!!!!!!!!!!!!!!!!!
7. For each initial condition in the [CAM_clone]/src/dynamics/tests/initial_conditions
directory, modify these scripts to use W. This requires three steps:
 a) Add W as an argument to the [IC_name]_set_ic() routine.
 b) In the ! Dummy arguments list, add:

   real(r8), optional, intent(inout) :: W(:,:)     ! vertical velocity (nonhydrostatic)

 c) POSSIBLE EXTRA:
 If there are variables lu, lv ... (e.g in ic_baro_dry_jw06.F90)
 modify the corresponding !Local variables line to add lw:

 logical                           :: lu,lv,lw,lt,lq,l3d_vars

and when initialising the 3D vars, modify the code to:

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
l3d_vars = lu .or. lv .or. lw .or. lt .or.lq
nlev = -1
  if (l3d_vars) then
  if (lu) nlev = size(U, 2)
  if (lv) nlev = size(V, 2)
  if (lw) nlev = size(W, 2)
  if (lt) nlev = size(T, 2)
  if (lq) nlev = size(Q, 2)

c) Add an initialisation of W, 
e.g. set to zero for held_suarez/us_standard_atmosphere by checking if W is present:

if (present(W)) then
  nlev = size(W, 2)
  do k = 1, nlev
    where(mask_use)
      W(:,k) = 0.0_r8
    end where
  end do
  if(masterproc .and. verbose_use) then
    write(iulog,*) '          W (nonhydrostatic) initialized by "',subname,'"'
  end if
end if

set to zero for baro_dry_jw06/baroclinic by checking if lw is present:

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


!!!!!!!!!!!!!!!!!
8. Additional steps for small Earth tests.

a) Comment out a line in [CAM_clone]/cime_config/config_component.xml that 
forces the radius to the large Earth value:
Line 362:

<value grid="a%C[0-9]"> rearth = 6.37122D6 </value>

b)
The small Earth radius and omega need to be set in fv_arrays.F90 instead.
First, copy this into the src_override:
cp [CAM_clone]/src/dynamics/fv3/atmos_cubed_sphere/model/fv_arrays.F90 [CAM_clone]/src/dynamics/fv3/src_override/fv_arrays.F90
