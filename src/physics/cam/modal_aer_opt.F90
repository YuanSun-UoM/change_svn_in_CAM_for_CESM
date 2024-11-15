module modal_aer_opt

! parameterizes aerosol coefficients using chebychev polynomial
! parameterize aerosol radiative properties in terms of
! surface mode wet radius and wet refractive index

! Ghan and Zaveri, JGR 2007.

! uses Wiscombe's (1979) mie scattering code


use shr_kind_mod,      only: r8 => shr_kind_r8, shr_kind_cl
use ppgrid,            only: pcols, pver, pverp
use constituents,      only: pcnst
use spmd_utils,        only: masterproc
use ref_pres,          only: top_lev => clim_modal_aero_top_lev
use physconst,         only: rhoh2o, rga, rair
use radconstants,      only: nswbands, nlwbands, idx_sw_diag, idx_uv_diag, idx_nir_diag
use rad_constituents,  only: n_diag, rad_cnst_get_call_list, rad_cnst_get_info, rad_cnst_get_aer_mmr, &
                             rad_cnst_get_aer_props, rad_cnst_get_mode_props, rad_cnst_get_mode_num
use physics_types,     only: physics_state

use physics_buffer, only : pbuf_get_index,physics_buffer_desc, pbuf_get_field
use pio,               only: file_desc_t, var_desc_t, pio_inq_dimlen, pio_inq_dimid, pio_inq_varid, &
                             pio_get_var, pio_nowrite, pio_closefile
use cam_pio_utils,     only: cam_pio_openfile
use cam_history,       only: addfld, add_default, outfld, horiz_only
use cam_history_support, only: fillvalue
use cam_logfile,       only: iulog
use perf_mod,          only: t_startf, t_stopf
use cam_abortutils,    only: endrun

use modal_aero_wateruptake, only: modal_aero_wateruptake_dr
use modal_aero_calcsize,    only: modal_aero_calcsize_diag
use wv_saturation, only: qsat_water

implicit none
private
save

public :: modal_aer_opt_readnl, modal_aer_opt_init, modal_aero_sw, modal_aero_lw


character(len=*), parameter :: unset_str = 'UNSET'

! Namelist variables:
character(shr_kind_cl)      :: modal_optics_file = unset_str   ! full pathname for modal optics dataset
character(shr_kind_cl)      :: water_refindex_file = unset_str ! full pathname for water refractive index dataset

! Dimension sizes in coefficient arrays used to parameterize aerosol radiative properties
! in terms of refractive index and wet radius
integer, parameter :: ncoef=5, prefr=7, prefi=10

real(r8) :: xrmin, xrmax

! refractive index for water read in read_water_refindex
complex(r8) :: crefwsw(nswbands) ! complex refractive index for water visible
complex(r8) :: crefwlw(nlwbands) ! complex refractive index for water infrared

! physics buffer indices
integer :: dgnumwet_idx = -1
integer :: qaerwat_idx  = -1

character(len=4) :: diag(0:n_diag) = (/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ', &
                                       '_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)
    integer,parameter :: nx = 5     !! number of points in x
    integer,parameter :: ny = 5     !! number of points in y
    integer,parameter :: nz = 50     !! number of points in z
   integer, parameter :: bins=30
    real(r8) :: cabs_table(nx,ny,nz,bins),cext_table(nx,ny,nz,bins)
    real(r8) :: cabs_nobc_table(nx,ny,nz,bins),cext_nobc_table(nx,ny,nz,bins)
    real(r8) :: asy_table(nx,ny,nz,bins)
!===============================================================================
CONTAINS
!===============================================================================

subroutine modal_aer_opt_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'modal_aer_opt_readnl'

   namelist /modal_aer_opt_nl/ water_refindex_file
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'modal_aer_opt_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, modal_aer_opt_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   call mpibcast(water_refindex_file, len(water_refindex_file), mpichar, 0, mpicom)
#endif


end subroutine modal_aer_opt_readnl

!===============================================================================

subroutine modal_aer_opt_init()

   use ioFileMod,        only: getfil
   use phys_control,     only: phys_getopts

   ! Local variables

   integer  :: i, m
   real(r8) :: rmmin, rmmax       ! min, max aerosol surface mode radius treated (m)
   character(len=256) :: locfile
   
   logical           :: history_amwg            ! output the variables used by the AMWG diag package
   logical           :: history_aero_optics     ! output aerosol optics diagnostics
   logical           :: history_dust            ! output dust diagnostics

   logical :: call_list(0:n_diag)
   integer :: ilist, nmodes, m_ncoef, m_prefr, m_prefi
   integer :: errcode

   character(len=*), parameter :: routine='modal_aer_opt_init'
   character(len=10) :: fldname
   character(len=128) :: lngname

    integer :: ii,jj,kk,idi
   !----------------------------------------------------------------------------
 open (70, file='/nuist/scratch/liuc/ChenGZ/my_cesm_sandbox/cime/scripts/mietable/cabs_table.txt')
 open (75, file='/nuist/scratch/liuc/ChenGZ/my_cesm_sandbox/cime/scripts/mietable/cext_table.txt')
 open (80, file='/nuist/scratch/liuc/ChenGZ/my_cesm_sandbox/cime/scripts/mietable/cabs_nobc_table.txt')
 open (85, file='/nuist/scratch/liuc/ChenGZ/my_cesm_sandbox/cime/scripts/mietable/cext_nobc_table.txt')
 open (90, file='/nuist/scratch/liuc/ChenGZ/my_cesm_sandbox/cime/scripts/mietable/asy_table.txt')
  do idi=1,bins ! refr_shell(10)
    do kk=1,nz  ! refi_shell(5)
      do jj=1,ny  ! Rp_bc
        do ii=1,nx ! Rp
    
    read(70,*)cabs_table(ii,jj,kk,idi)
    read(75,*)cext_table(ii,jj,kk,idi)
    read(80,*)cabs_nobc_table(ii,jj,kk,idi)
    read(85,*)cext_nobc_table(ii,jj,kk,idi)
    read(90,*)asy_table(ii,jj,kk,idi)

	enddo
      enddo
    enddo
  enddo
close(70)
close(75)
close(80)
close(85)
close(90)

   rmmin = 0.01e-6_r8
   rmmax = 25.e-6_r8
   xrmin = log(rmmin)
   xrmax = log(rmmax)

   ! Check that dimension sizes in the coefficient arrays used to
   ! parameterize aerosol radiative properties are consistent between this
   ! module and the mode physprop files.
   call rad_cnst_get_call_list(call_list)
   do ilist = 0, n_diag
      if (call_list(ilist)) then
         call rad_cnst_get_info(ilist, nmodes=nmodes)
         do m = 1, nmodes
            call rad_cnst_get_mode_props(ilist, m, ncoef=m_ncoef, prefr=m_prefr, prefi=m_prefi)
            if (m_ncoef /= ncoef .or. m_prefr /= prefr .or. m_prefi /= prefi) then
               write(iulog,*) routine//': ERROR - file and module values do not match:'
               write(iulog,*) '   ncoef:', ncoef, m_ncoef
               write(iulog,*) '   prefr:', prefr, m_prefr
               write(iulog,*) '   prefi:', prefi, m_prefi
               call endrun(routine//': ERROR - file and module values do not match')
            end if
         end do
      end if
   end do

   ! Initialize physics buffer indices for dgnumwet and qaerwat.  Note the implicit assumption
   ! that the loops over modes in the optics calculations will use the values for dgnumwet and qaerwat
   ! that are set in the aerosol_wet_intr code.
   dgnumwet_idx = pbuf_get_index('DGNUMWET',errcode)
   if (errcode < 0) then
      call endrun(routine//' ERROR: cannot find physics buffer field DGNUMWET')
   end if
   qaerwat_idx  = pbuf_get_index('QAERWAT',errcode)
   if (errcode < 0) then
      call endrun(routine//' ERROR: cannot find physics buffer field QAERWAT')
   end if

   call getfil(water_refindex_file, locfile)
   call read_water_refindex(locfile)
   if (masterproc) write(iulog,*) "modal_aer_opt_init: read water refractive index file:", trim(locfile)

   call phys_getopts(history_amwg_out        = history_amwg, &
                     history_aero_optics_out = history_aero_optics, &
                     history_dust_out        = history_dust )

   ! Add diagnostic fields to history output.

   call addfld ('EXTINCT',    (/ 'lev' /), 'A','/m','Aerosol extinction 550 nm, day only',                   &
                flag_xyfill=.true.)
   call addfld ('EXTINCTUV',  (/ 'lev' /), 'A','/m','Aerosol extinction 350 nm, day only',                   &
                flag_xyfill=.true.)
   call addfld ('EXTINCTNIR', (/ 'lev' /), 'A','/m','Aerosol extinction 1020 nm, day only',                  &
                flag_xyfill=.true.)
   call addfld ('ABSORB',     (/ 'lev' /), 'A','/m','Aerosol absorption, day only',                          &
                flag_xyfill=.true.)
   call addfld ('AODVIS',     horiz_only,  'A','  ','Aerosol optical depth 550 nm, day only',                &
                flag_xyfill=.true.)
   call addfld ('AODVISst',   horiz_only,  'A','  ','Stratospheric aerosol optical depth 550 nm, day only',  &
                flag_xyfill=.true.)
   call addfld ('AODUV',      horiz_only,  'A','  ','Aerosol optical depth 350 nm, day only',                &
                flag_xyfill=.true.)
   call addfld ('AODUVst',    horiz_only,  'A','  ','Stratospheric aerosol optical depth 350 nm, day only',  &
                flag_xyfill=.true.)
   call addfld ('AODNIR',     horiz_only,  'A','  ','Aerosol optical depth 1020 nm, day only',               &
                flag_xyfill=.true.)
   call addfld ('AODNIRst',   horiz_only,  'A','  ','Stratospheric aerosol optical depth 1020 nm, day only', &
                flag_xyfill=.true.)
   call addfld ('AODABS',     horiz_only,  'A','  ','Aerosol absorption optical depth 550 nm, day only',     &
                flag_xyfill=.true.)
   call addfld ('AODxASYM',   horiz_only,  'A','  ','Aerosol optical depth 550 * asymmetry factor, day only',&
                flag_xyfill=.true.)
   call addfld ('EXTxASYM',   (/ 'lev' /), 'A','  ','extinction 550 nm * asymmetry factor, day only',        &
                flag_xyfill=.true.)

   call addfld ('EXTINCTdn',    (/ 'lev' /), 'A','/m','Aerosol extinction 550 nm, day night',                 &
                  flag_xyfill=.true.)
   call addfld ('EXTINCTUVdn',  (/ 'lev' /), 'A','/m','Aerosol extinction 350 nm, day night',                 &
                  flag_xyfill=.true.)
   call addfld ('EXTINCTNIRdn', (/ 'lev' /), 'A','/m','Aerosol extinction 1020 nm, day night',                &
                  flag_xyfill=.true.)
   call addfld ('ABSORBdn',     (/ 'lev' /), 'A','/m','Aerosol absorption, day night',                        &
                  flag_xyfill=.true.)
   call addfld ('AODVISdn',     horiz_only,  'A','  ','Aerosol optical depth 550 nm, day night',              &
                  flag_xyfill=.true.)
   call addfld ('AODVISstdn',   horiz_only,  'A','  ','Stratospheric aerosol optical depth 550 nm, day night',&
                  flag_xyfill=.true.)
   call addfld ('AODUVdn',      horiz_only,  'A','  ','Aerosol optical depth 350 nm, day night',              &
                  flag_xyfill=.true.)
   call addfld ('AODUVstdn',    horiz_only,  'A','  ','Stratospheric aerosol optical depth 350 nm, day night',&
                  flag_xyfill=.true.)
   call addfld ('AODNIRdn',     horiz_only,  'A','  ','Aerosol optical depth 1020 nm, day night',             &
                  flag_xyfill=.true.)
   call addfld ('AODNIRstdn',   horiz_only,  'A','  ','Stratospheric aerosol optical depth 1020 nm, day night',&
                 flag_xyfill=.true.)
   call addfld ('AODABSdn',     horiz_only,  'A','  ','Aerosol absorption optical depth 550 nm, day night',   &
                  flag_xyfill=.true.)
   call addfld ('AODxASYMdn',   horiz_only,  'A','  ','Aerosol optical depth 550 * asymmetry factor, day night',&
                flag_xyfill=.true.)
   call addfld ('EXTxASYMdn',   (/ 'lev' /), 'A','  ','extinction 550 * asymmetry factor, day night',        &
                flag_xyfill=.true.)
   call addfld ('refr', (/ 'lev' /), 'A','  ','refr', flag_xyfill=.true.)
   call addfld ('refi', (/ 'lev' /), 'A','  ','refi', flag_xyfill=.true.)
   call addfld ('num', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.) 
   call addfld ('specdens_m', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('mass', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('air_density', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('vol_fra', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('refr_shell', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('refi_shell', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('vol_fradust', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('vol_frass', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('vol_frasoa', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('vol_frapoa', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('vol_frasul', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('vol_fram4',(/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('vol_frawater',(/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('refr_core', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('refi_core', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('aodabsmode1', horiz_only, 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('specdens_mwet', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('specdens_nobc', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('mass_water', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('absorbBC', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('vol_all', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('mass_bc', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('absorbBC4', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('pabs_bc', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
   call addfld ('mass_bc4', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
         call addfld ('palb_layer', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
         call addfld ('tau_layer', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
         call addfld ('asym_layer', (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)


   call rad_cnst_get_info(0, nmodes=nmodes)

   do m = 1, nmodes

      write(fldname,'(a,i1)') 'BURDEN', m
      write(lngname,'(a,i1)') 'Aerosol burden, day only, mode ', m
      call addfld (fldname, horiz_only, 'A', 'kg/m2', lngname, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

      write(fldname,'(a,i1)') 'AODMODE', m
      write(lngname,'(a,i1)') 'Aerosol optical depth, day only, 550 nm mode ', m
      call addfld (fldname, horiz_only, 'A', '  ', lngname, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

      write(fldname,'(a,i1)') 'AODDUST', m
      write(lngname,'(a,i1,a)') 'Aerosol optical depth, day only, 550 nm mode ',m,' from dust'
      call addfld (fldname, horiz_only, 'A', '  ', lngname, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

      write(fldname,'(a,i1)') 'BURDENdn', m
      write(lngname,'(a,i1)') 'Aerosol burden, day night, mode ', m
      call addfld (fldname, horiz_only, 'A', 'kg/m2', lngname, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

      write(fldname,'(a,i1)') 'AODdnMODE', m
      write(lngname,'(a,i1)') 'Aerosol optical depth 550 nm, day night, mode ', m
      call addfld (fldname, horiz_only, 'A', '  ', lngname, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

      write(fldname,'(a,i1)') 'AODdnDUST', m
      write(lngname,'(a,i1,a)') 'Aerosol optical depth 550 nm, day night, mode ',m,' from dust'
      call addfld (fldname, horiz_only, 'A', '  ', lngname, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

      write(fldname,'(a,i1)') 'radsurf', m
      write(lngname,'(a,i1,a)') 'radsurf',m
      call addfld (fldname,(/ 'lev' /),'A', 'm', lngname, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

      write(fldname,'(a,i1)') 'specpext', m
      write(lngname,'(a,i1,a)') 'specpext',m
      call addfld (fldname,(/ 'lev' /),'A', ' ', lngname, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

      write(fldname,'(a,i1)') 'specpabs', m
      write(lngname,'(a,i1,a)') 'specpabs',m
      call addfld (fldname,(/ 'lev' /),'A', ' ', lngname, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

      write(fldname,'(a,i1)') 'wetvols',m
      write(lngname,'(a,i1,a)') 'wetvols',m
      call addfld (fldname,(/ 'lev' /),'A', ' ', lngname, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

   enddo

   call addfld ('AODDUST',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust, day only',         &
        flag_xyfill=.true.)
   call addfld ('AODSO4',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from SO4, day only',          &
        flag_xyfill=.true.)
   call addfld ('AODPOM',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from POM, day only',          &
        flag_xyfill=.true.)
   call addfld ('AODSOA',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from SOA, day only',          &
        flag_xyfill=.true.)
   call addfld ('AODBC',         horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from BC, day only',           &
        flag_xyfill=.true.)
   call addfld ('AODSS',         horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from seasalt, day only',      &
        flag_xyfill=.true.)
   call addfld ('AODABSBC',      horiz_only, 'A','  ',    'Aerosol absorption optical depth 550 nm from BC, day only',&
        flag_xyfill=.true.)
   call addfld ('BURDENDUST',    horiz_only, 'A','kg/m2', 'Dust aerosol burden, day only'        ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSO4',     horiz_only, 'A','kg/m2', 'Sulfate aerosol burden, day only'     ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENPOM',     horiz_only, 'A','kg/m2', 'POM aerosol burden, day only'         ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSOA',     horiz_only, 'A','kg/m2', 'SOA aerosol burden, day only'         ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENBC',      horiz_only, 'A','kg/m2', 'Black carbon aerosol burden, day only',                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSEASALT', horiz_only, 'A','kg/m2', 'Seasalt aerosol burden, day only'     ,                    &
        flag_xyfill=.true.)
   call addfld ('SSAVIS',        horiz_only, 'A','  ',    'Aerosol single-scatter albedo, day only',                  &
        flag_xyfill=.true.)

   call addfld ('AODDUSTdn',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust, day night',         &
        flag_xyfill=.true.)
   call addfld ('AODSO4dn',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from SO4, day night',          &
        flag_xyfill=.true.)
   call addfld ('AODPOMdn',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from POM, day night',          &
        flag_xyfill=.true.)
   call addfld ('AODSOAdn',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from SOA, day night',          &
        flag_xyfill=.true.)
   call addfld ('AODBCdn',         horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from BC, day night',           &
        flag_xyfill=.true.)
   call addfld ('AODSSdn',         horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from seasalt, day night',      &
        flag_xyfill=.true.)
   call addfld ('AODABSBCdn',      horiz_only, 'A','  ',    'Aerosol absorption optical depth 550 nm from BC, day night',&
        flag_xyfill=.true.)
   call addfld ('BURDENDUSTdn',    horiz_only, 'A','kg/m2', 'Dust aerosol burden, day night'        ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSO4dn',     horiz_only, 'A','kg/m2', 'Sulfate aerosol burden, day night'     ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENPOMdn',     horiz_only, 'A','kg/m2', 'POM aerosol burden, day night'         ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSOAdn',     horiz_only, 'A','kg/m2', 'SOA aerosol burden, day night'         ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENBCdn',      horiz_only, 'A','kg/m2', 'Black carbon aerosol burden, day night',                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSEASALTdn', horiz_only, 'A','kg/m2', 'Seasalt aerosol burden, day night'     ,                    &
        flag_xyfill=.true.)
   call addfld ('SSAVISdn',        horiz_only, 'A','  ',    'Aerosol single-scatter albedo, day night',                  &
        flag_xyfill=.true.)

 
   if (history_amwg) then 
      call add_default ('AODDUST1'     , 1, ' ')
      call add_default ('AODDUST3'     , 1, ' ')
      call add_default ('AODDUST'      , 1, ' ')
      call add_default ('AODVIS'       , 1, ' ')
   end if

   if (history_dust) then 
      call add_default ('AODDUST1'     , 1, ' ')
      call add_default ('AODDUST2'     , 1, ' ')
      call add_default ('AODDUST3'     , 1, ' ')
   end if

   if (history_aero_optics) then 
      call add_default ('AODDUST1'     , 1, ' ')
      call add_default ('AODDUST3'     , 1, ' ')
      call add_default ('ABSORB'       , 1, ' ')
      call add_default ('AODMODE1'     , 1, ' ')
      call add_default ('AODMODE2'     , 1, ' ')
      call add_default ('AODMODE3'     , 1, ' ')
      call add_default ('AODVIS'       , 1, ' ')
      call add_default ('AODUV'        , 1, ' ')
      call add_default ('AODNIR'       , 1, ' ')
      call add_default ('AODABS'       , 1, ' ')
      call add_default ('AODABSBC'     , 1, ' ')
      call add_default ('AODDUST'      , 1, ' ')
      call add_default ('AODSO4'       , 1, ' ')
      call add_default ('AODPOM'       , 1, ' ')
      call add_default ('AODSOA'       , 1, ' ')
      call add_default ('AODBC'        , 1, ' ')
      call add_default ('AODSS'        , 1, ' ')
      call add_default ('BURDEN1'      , 1, ' ')
      call add_default ('BURDEN2'      , 1, ' ')
      call add_default ('BURDEN3'      , 1, ' ')
      call add_default ('BURDENDUST'   , 1, ' ')
      call add_default ('BURDENSO4'    , 1, ' ')
      call add_default ('BURDENPOM'    , 1, ' ')
      call add_default ('BURDENSOA'    , 1, ' ')
      call add_default ('BURDENBC'     , 1, ' ')
      call add_default ('BURDENSEASALT', 1, ' ')
      call add_default ('SSAVIS'       , 1, ' ')
      call add_default ('EXTINCT'      , 1, ' ')
      call add_default ('AODxASYM'     , 1, ' ')
      call add_default ('EXTxASYM'     , 1, ' ')
     
      call add_default ('AODdnDUST1'     , 1, ' ')
      call add_default ('AODdnDUST3'     , 1, ' ')
      call add_default ('ABSORBdn'       , 1, ' ')
      call add_default ('AODdnMODE1'     , 1, ' ')
      call add_default ('AODdnMODE2'     , 1, ' ')
      call add_default ('AODdnMODE3'     , 1, ' ')
      call add_default ('AODVISdn'       , 1, ' ')
      call add_default ('AODUVdn'        , 1, ' ')
      call add_default ('AODNIRdn'       , 1, ' ')
      call add_default ('AODABSdn'       , 1, ' ')
      call add_default ('AODABSBCdn'     , 1, ' ')
      call add_default ('AODDUSTdn'      , 1, ' ')
      call add_default ('AODSO4dn'       , 1, ' ')
      call add_default ('AODPOMdn'       , 1, ' ')
      call add_default ('AODSOAdn'       , 1, ' ')
      call add_default ('AODBCdn'        , 1, ' ')
      call add_default ('AODSSdn'        , 1, ' ')
      call add_default ('BURDENdn1'      , 1, ' ')
      call add_default ('BURDENdn2'      , 1, ' ')
      call add_default ('BURDENdn3'      , 1, ' ')
      call add_default ('BURDENDUSTdn'   , 1, ' ')
      call add_default ('BURDENSO4dn'    , 1, ' ')
      call add_default ('BURDENPOMdn'    , 1, ' ')
      call add_default ('BURDENSOAdn'    , 1, ' ')
      call add_default ('BURDENBCdn'     , 1, ' ')
      call add_default ('BURDENSEASALTdn', 1, ' ')
      call add_default ('SSAVISdn'       , 1, ' ')
      call add_default ('EXTINCTdn'      , 1, ' ')
      call add_default ('AODxASYMdn'     , 1, ' ')
      call add_default ('EXTxASYMdn'     , 1, ' ')
  end if

   do ilist = 1, n_diag
      if (call_list(ilist)) then
         
         call addfld ('EXTINCT'//diag(ilist),  (/ 'lev' /), 'A','/m', &
              'Aerosol extinction', flag_xyfill=.true.)
         call addfld ('ABSORB'//diag(ilist),   (/ 'lev' /), 'A','/m', &
              'Aerosol absorption', flag_xyfill=.true.)
         call addfld ('AODVIS'//diag(ilist),   horiz_only,  'A','  ', &
              'Aerosol optical depth 550 nm', flag_xyfill=.true.)
         call addfld ('AODVISst'//diag(ilist), horiz_only,  'A','  ', &
              'Stratospheric aerosol optical depth 550 nm', flag_xyfill=.true.)
         call addfld ('AODABS'//diag(ilist),   horiz_only,  'A','  ', &
              'Aerosol absorption optical depth 550 nm', flag_xyfill=.true.)

         call addfld ('EXTINCTdn'//diag(ilist),    (/ 'lev' /), 'A','/m',&
              'Aerosol extinction 550 nm, day night', flag_xyfill=.true.)
         call addfld ('ABSORBdn'//diag(ilist),     (/ 'lev' /), 'A','/m',&
              'Aerosol absorption, day night',  flag_xyfill=.true.)
         call addfld ('AODVISdn'//diag(ilist),     horiz_only,  'A','  ',&
              'Aerosol optical depth 550 nm, day night',  flag_xyfill=.true.)
         call addfld ('AODVISstdn'//diag(ilist),   horiz_only,  'A','  ',&
              'Stratospheric aerosol optical depth 550 nm, day night', flag_xyfill=.true.)
         call addfld ('AODABSdn'//diag(ilist),     horiz_only,  'A','  ',&
              'Aerosol absorption optical depth 550 nm, day night',  flag_xyfill=.true.)
         call addfld ('EXTxASYMdn'//diag(ilist),   (/ 'lev' /), 'A','  ',&
              'extinction 550 * asymmetry factor, day night',  flag_xyfill=.true.)
         call addfld ('EXTxASYM'//diag(ilist),   (/ 'lev' /), 'A','  ',&
              'extinction 550 nm * asymmetry factor, day only',   flag_xyfill=.true.)
         call addfld ('palb_layer'//diag(ilist), (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
         call addfld ('tau_layer'//diag(ilist), (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)
         call addfld ('asym_layer'//diag(ilist), (/ 'lev' /), 'A','  ',' ', flag_xyfill=.true.)

         if (history_aero_optics) then
            call add_default ('EXTINCT'//diag(ilist), 1, ' ')
            call add_default ('ABSORB'//diag(ilist),  1, ' ')
            call add_default ('AODVIS'//diag(ilist),  1, ' ')
            call add_default ('AODVISst'//diag(ilist),  1, ' ')
            call add_default ('AODABS'//diag(ilist),  1, ' ')
         end if

      end if
   end do

end subroutine modal_aer_opt_init

!===============================================================================

subroutine modal_aero_sw(list_idx, state, pbuf, nnite, idxnite, &
                         tauxar, wa, ga, fa)

   ! calculates aerosol sw radiative properties
   
   use tropopause, only : tropopause_findChemTrop

   integer,             intent(in) :: list_idx       ! index of the climate or a diagnostic list
   type(physics_state), intent(in), target :: state          ! state variables
   
   type(physics_buffer_desc), pointer :: pbuf(:)
   integer,             intent(in) :: nnite          ! number of night columns
   integer,             intent(in) :: idxnite(nnite) ! local column indices of night columns

   real(r8), intent(out) :: tauxar(pcols,0:pver,nswbands) ! layer extinction optical depth
   real(r8), intent(out) :: wa(pcols,0:pver,nswbands)     ! layer single-scatter albedo
   real(r8), intent(out) :: ga(pcols,0:pver,nswbands)     ! asymmetry factor
   real(r8), intent(out) :: fa(pcols,0:pver,nswbands)     ! forward scattered fraction

   ! Local variables
   integer :: i, ifld, isw, k, l, m, nc, ns
   integer :: lchnk                    ! chunk id
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nmodes
   integer :: nspec
   integer :: troplevchem(pcols)       ! Chemical tropopause level
   integer :: istat

   real(r8) :: mass(pcols,pver)        ! layer mass
   real(r8) :: air_density(pcols,pver) ! (kg/m3)

   real(r8),    pointer :: specmmr(:,:)        ! species mass mixing ratio
   real(r8)             :: specdens            ! species density (kg/m3)
   complex(r8), pointer :: specrefindex(:)     ! species refractive index
   character*32         :: spectype            ! species type
   real(r8)             :: hygro_aer           ! 

   real(r8), pointer :: dgnumwet(:,:)     ! number mode wet diameter
   real(r8), pointer :: qaerwat(:,:)      ! aerosol water (g/g)

   real(r8), pointer :: dgnumdry_m(:,:,:) ! number mode dry diameter for all modes
   real(r8), pointer :: dgnumwet_m(:,:,:) ! number mode wet diameter for all modes
   real(r8), pointer :: qaerwat_m(:,:,:)  ! aerosol water (g/g) for all modes
   real(r8), pointer :: wetdens_m(:,:,:)  ! 
   real(r8), pointer :: hygro_m(:,:,:)  ! 
   real(r8), pointer :: dryvol_m(:,:,:)  ! 
   real(r8), pointer :: dryrad_m(:,:,:)  ! 
   real(r8), pointer :: drymass_m(:,:,:)  ! 
   real(r8), pointer :: so4dryvol_m(:,:,:)  ! 
   real(r8), pointer :: naer_m(:,:,:)  ! 

   real(r8) :: sigma_logr_aer         ! geometric standard deviation of number distribution
   real(r8) :: radsurf(pcols,pver)    ! aerosol surface mode radius
   real(r8) :: logradsurf(pcols,pver) ! log(aerosol surface mode radius)
   real(r8) :: cheb(ncoef,pcols,pver)

   real(r8)    :: refr(pcols)     ! real part of refractive index
   real(r8)    :: refi(pcols)     ! imaginary part of refractive index
   complex(r8) :: crefin(pcols)   ! complex refractive index
   real(r8), pointer :: refrtabsw(:,:) ! table of real refractive indices for aerosols
   real(r8), pointer :: refitabsw(:,:) ! table of imag refractive indices for aerosols
   real(r8), pointer :: extpsw(:,:,:,:) ! specific extinction
   real(r8), pointer :: abspsw(:,:,:,:) ! specific absorption
   real(r8), pointer :: asmpsw(:,:,:,:) ! asymmetry factor

   real(r8) :: vol(pcols)      ! volume concentration of aerosol specie (m3/kg)
   real(r8) :: dryvol(pcols)   ! volume concentration of aerosol mode (m3/kg)
   real(r8) :: watervol(pcols) ! volume concentration of water in each mode (m3/kg)
   real(r8) :: wetvol(pcols)   ! volume concentration of wet mode (m3/kg)

   integer  :: itab(pcols), jtab(pcols)
   real(r8) :: ttab(pcols), utab(pcols)
   real(r8) :: cext(pcols,ncoef), cabs(pcols,ncoef), casm(pcols,ncoef)
   real(r8) :: pext(pcols)     ! parameterized specific extinction (m2/kg)
   real(r8) :: specpext(pcols) ! specific extinction (m2/kg)
   real(r8) :: dopaer(pcols)   ! aerosol optical depth in layer
   real(r8) :: pabs(pcols)     ! parameterized specific absorption (m2/kg)
   real(r8) :: pasm(pcols)     ! parameterized asymmetry factor
   real(r8) :: palb(pcols)     ! parameterized single scattering albedo

   ! Diagnostics
   real(r8) :: extinct(pcols,pver)
   real(r8) :: extinctnir(pcols,pver)
   real(r8) :: extinctuv(pcols,pver)
   real(r8) :: absorb(pcols,pver)
   real(r8) :: aodvis(pcols)               ! extinction optical depth
   real(r8) :: aodvisst(pcols)             ! stratospheric extinction optical depth
   real(r8) :: aodabs(pcols)               ! absorption optical depth
   real(r8) :: asymvis(pcols)              ! asymmetry factor * optical depth
   real(r8) :: asymext(pcols,pver)         ! asymmetry factor * extinction

   real(r8) :: aodabsbc(pcols)             ! absorption optical depth of BC

   real(r8) :: ssavis(pcols)
   real(r8) :: dustvol(pcols)              ! volume concentration of dust in aerosol mode (m3/kg)

   real(r8) :: burden(pcols)
   real(r8) :: burdendust(pcols), burdenso4(pcols), burdenbc(pcols), &
               burdenpom(pcols), burdensoa(pcols), burdenseasalt(pcols)

   real(r8) :: aodmode(pcols)
   real(r8) :: dustaodmode(pcols)          ! dust aod in aerosol mode

   real(r8) :: specrefr, specrefi
   real(r8) :: scatdust(pcols), scatso4(pcols), scatbc(pcols), &
               scatpom(pcols), scatsoa(pcols), scatseasalt(pcols)
   real(r8) :: absdust(pcols), absso4(pcols), absbc(pcols), &
               abspom(pcols), abssoa(pcols), absseasalt(pcols)
   real(r8) :: hygrodust(pcols), hygroso4(pcols), hygrobc(pcols), &
               hygropom(pcols), hygrosoa(pcols), hygroseasalt(pcols)

   real(r8) :: scath2o, absh2o, sumscat, sumabs, sumhygro
   real(r8) :: aodc                        ! aod of component

   ! total species AOD
   real(r8) :: dustaod(pcols), so4aod(pcols), bcaod(pcols), &
               pomaod(pcols), soaaod(pcols), seasaltaod(pcols)




   logical :: savaervis ! true if visible wavelength (0.55 micron)
   logical :: savaernir ! true if near ir wavelength (~0.88 micron)
   logical :: savaeruv  ! true if uv wavelength (~0.35 micron)

   real(r8) :: aoduv(pcols)               ! extinction optical depth in uv
   real(r8) :: aoduvst(pcols)             ! stratospheric extinction optical depth in uv
   real(r8) :: aodnir(pcols)              ! extinction optical depth in nir
   real(r8) :: aodnirst(pcols)            ! stratospheric extinction optical depth in nir


   character(len=32) :: outname

   ! debug output
   integer, parameter :: nerrmax_dopaer=1000
   integer  :: nerr_dopaer = 0
   real(r8) :: volf            ! volume fraction of insoluble aerosol
   character(len=*), parameter :: subname = 'modal_aero_sw'

   real(r8)    :: specmmr_m(pcols)
   real(r8)    :: specmmr_m1(pcols,pver)
   real(r8)             :: specdens_m1
   real(r8)             :: specdens_m(pcols,pver)
   real(r8) :: specdens_mwet(pcols,pver)
   real(r8) :: specdens_nobc(pcols,pver)
   real(r8)    :: refr1(pcols,pver)
   real(r8)    :: refr_shell(pcols,pver)
   real(r8)    :: refi_shell(pcols,pver)
   real(r8)    :: refr_core(pcols,pver)
   real(r8)    :: refi_core(pcols,pver)
   real(r8)    :: refi1(pcols,pver)
   real(r8) :: vol_m(pcols)
   real(r8) :: dryvol_mode(pcols)
   real(r8) :: wetvols(pcols,pver)
   real(r8) :: vol_sol(pcols)
   real(r8) :: dryvol_sol(pcols)
   real(r8) :: vol_insol(pcols)
   real(r8) :: vol_dust(pcols)
   real(r8) :: vol_poa(pcols)
   real(r8) :: vol_soa(pcols)
   real(r8) :: vol_sul(pcols)
   real(r8) :: vol_ss(pcols)
   real(r8) :: dryvol_insol(pcols)
   real(r8) :: dryvol_dust(pcols)
   real(r8) :: dryvol_poa(pcols)
   real(r8) :: dryvol_soa(pcols)
   real(r8) :: dryvol_sul(pcols)
   real(r8) :: dryvol_ss(pcols)
   complex(r8) :: crefin_sol(pcols)
   complex(r8) :: crefin_insol(pcols)
   real(r8) :: vol_fra(pcols,pver)
   real(r8) :: vol_frasul(pcols,pver)
   real(r8) :: vol_fradust(pcols,pver)
   real(r8) :: vol_frapoa(pcols,pver)
   real(r8) :: vol_frasoa(pcols,pver)
   real(r8) :: vol_frass(pcols,pver)
   real(r8) :: vol_fram4(pcols,pver)
   complex(r8) :: A(pcols)  ! used in maxwell-Garnett mixing
   complex(r8) :: B(pcols)
   real(r8) :: specpabs(pcols) ! specific absorption (m2/kg)
   real(r8) :: specpext_vis(pcols,pver)
   real(r8) :: specpabs_vis(pcols,pver)
   real(r8) :: aodabsmode1(pcols)
   real(r8), pointer :: mode_num(:,:)
   real(r8) :: mass_water(pcols,pver)
   real(r8) :: specmmr_nobc(pcols)
   real(r8) :: absorbBC(pcols,pver)
   real(r8) :: absorbBC4(pcols,pver)
   real(r8) :: vol_all(pcols,pver)
   real(r8) :: vol_frawater(pcols,pver)
   real(r8) :: mass_bc(pcols,pver)
   real(r8) :: mass_bc4(pcols,pver)
   real(r8) :: palb_layer(pcols,pver)
   real(r8) :: tau_layer(pcols,pver)
   real(r8) :: asym_layer(pcols,pver)

   ! new-coreshell
   real(r8),parameter :: pi=3.14159
   real(r8) :: c1, c2, c, Rp(bins), Rm(pcols,pver), Rm_bc(pcols,pver)
   real(r8) :: volume_bc(bins), n_bc, n_all(bins)
   real(r8) :: N0_bc(pcols,pver), xman,rzc,asy1,asy2,asy3
   real(r8) :: out, out1, out2, out3, out4, out5, out6, out7
   real(r8) :: yy(bins), yy1(bins), yy2(bins),yy3(bins),yy4(bins)
   real(r8) :: yy5(bins), yy6(bins), yy7(bins)
   integer :: k1, j, isize
   REAL(r8) :: GSCA,QEXT,QSCA,qback,qabs,xcor
   real(r8) :: resu,imsu,rebc,imbc
   COMPLEX(r8) :: REFREL,m1,m2,s1,s2
   Integer,  parameter :: nang= 50
   real(r8) :: theta(nang), n1(bins), n2(bins), pabs_bc(pcols,pver), pext_bc(pcols,pver)
   real(r8) :: xcor_nobc, xman_nobc, qabs_nobc, qext_nobc, qsca_nobc
   real(r8) :: n_bc1, rad, dens_wet(pcols,pver), Reff_bc(pcols,pver)
   real(r8) :: rad_shell(pcols,pver),Rp_bc(bins)
   real(r8) :: qabs_lens,qext_lens,qsca_lens
   COMPLEX(r8) :: m2_lens
   real(r8) :: pabs_para(pcols)
   real(r8) :: pext_para(pcols)
   real(r8) :: sigma_m          ! geometric standard deviation of number distribution for accumulation mode
   real(r8) :: Reff_nobc(pcols,pver)
   real(r8) :: N0_nobc(pcols,pver)

   ! mie lookup table
    real(r8) :: ccabs(bins),ccext(bins),cabs_nobc(bins),cext_nobc(bins),asy(bins)
   real(r8) :: vol_fra_input(pcols),refr_shell_input(pcols),refi_shell_input(pcols)
   real(r8) :: dens_wet_input(pcols)

   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   ! initialize output variables
   tauxar(:ncol,:,:) = 0._r8
   wa(:ncol,:,:)     = 0._r8
   ga(:ncol,:,:)     = 0._r8
   fa(:ncol,:,:)     = 0._r8

   ! zero'th layer does not contain aerosol
   tauxar(1:ncol,0,:)  = 0._r8
   wa(1:ncol,0,:)      = 0.925_r8
   ga(1:ncol,0,:)      = 0.850_r8
   fa(1:ncol,0,:)      = 0.7225_r8

   mass(:ncol,:)        = state%pdeldry(:ncol,:)*rga
   air_density(:ncol,:) = state%pmid(:ncol,:)/(rair*state%t(:ncol,:))

   ! diagnostics for visible band summed over modes
   extinct(1:ncol,:)     = 0.0_r8
   absorb(1:ncol,:)      = 0.0_r8
   absorbBC(1:ncol,:)      = 0.0_r8
   absorbBC4(1:ncol,:)      = 0.0_r8
   aodvis(1:ncol)        = 0.0_r8
   aodvisst(1:ncol)      = 0.0_r8
   aodabs(1:ncol)        = 0.0_r8
   burdendust(:ncol)     = 0.0_r8
   burdenso4(:ncol)      = 0.0_r8
   burdenpom(:ncol)      = 0.0_r8
   burdensoa(:ncol)      = 0.0_r8
   burdenbc(:ncol)       = 0.0_r8
   burdenseasalt(:ncol)  = 0.0_r8
   ssavis(1:ncol)        = 0.0_r8
   asymvis(1:ncol)       = 0.0_r8
   asymext(1:ncol,:)     = 0.0_r8

   aodabsbc(:ncol)       = 0.0_r8 
   dustaod(:ncol)        = 0.0_r8
   so4aod(:ncol)         = 0.0_r8
   pomaod(:ncol)         = 0.0_r8
   soaaod(:ncol)         = 0.0_r8
   bcaod(:ncol)          = 0.0_r8
   seasaltaod(:ncol)     = 0.0_r8

   ! diags for other bands
   extinctuv(1:ncol,:)   = 0.0_r8
   extinctnir(1:ncol,:)  = 0.0_r8
   aoduv(:ncol)          = 0.0_r8
   aodnir(:ncol)         = 0.0_r8
   aoduvst(:ncol)        = 0.0_r8
   aodnirst(:ncol)       = 0.0_r8
   aodabsmode1(1:ncol)        = 0.0_r8
   call tropopause_findChemTrop(state, troplevchem)

   ! loop over all aerosol modes
   call rad_cnst_get_info(list_idx, nmodes=nmodes)
   if (list_idx == 0) then
      ! water uptake and wet radius for the climate list has already been calculated
      call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet_m)
      call pbuf_get_field(pbuf, qaerwat_idx,  qaerwat_m)
   else
      ! If doing a diagnostic calculation then need to calculate the wet radius
      ! and water uptake for the diagnostic modes
      allocate(dgnumdry_m(pcols,pver,nmodes),  dgnumwet_m(pcols,pver,nmodes), &
               qaerwat_m(pcols,pver,nmodes),   wetdens_m(pcols,pver,nmodes), &
               hygro_m(pcols,pver,nmodes),     dryvol_m(pcols,pver,nmodes), &
               dryrad_m(pcols,pver,nmodes),    drymass_m(pcols,pver,nmodes),  &
               so4dryvol_m(pcols,pver,nmodes), naer_m(pcols,pver,nmodes),     stat=istat)
      if (istat > 0) then
         call endrun('modal_aero_sw: allocation FAILURE: arrays for diagnostic calcs')
      end if
  ! remove BC component (keep Dp constant)
      call modal_aero_calcsize_diag(state, pbuf, list_idx, dgnumdry_m, hygro_m, &
                                    dryvol_m, dryrad_m, drymass_m, so4dryvol_m, naer_m)  
      call modal_aero_wateruptake_dr(state, pbuf, list_idx, dgnumdry_m, dgnumwet_m, &
                                     qaerwat_m, wetdens_m,  hygro_m, dryvol_m, dryrad_m, &
                                     drymass_m, so4dryvol_m, naer_m)
  !    call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet_m)
  !    call pbuf_get_field(pbuf, qaerwat_idx,  qaerwat_m)
  ! remove BC component (keep Dp constant)
   endif

   do m = 1, nmodes

      ! diagnostics for visible band for each mode
      burden(:ncol)       = 0._r8
      aodmode(1:ncol)     = 0.0_r8
      dustaodmode(1:ncol) = 0.0_r8

      dgnumwet => dgnumwet_m(:,:,m)
      qaerwat  => qaerwat_m(:,:,m)

      ! get mode properties
      call rad_cnst_get_mode_props(list_idx, m, sigmag=sigma_logr_aer, refrtabsw=refrtabsw , &
         refitabsw=refitabsw, extpsw=extpsw, abspsw=abspsw, asmpsw=asmpsw)

      ! get mode info
      call rad_cnst_get_info(list_idx, m, nspec=nspec)

      ! calc size parameter for all columns
      call modal_size_parameters(ncol, sigma_logr_aer, dgnumwet, radsurf, logradsurf, cheb)

      ! get mode number mixing ratio
      call rad_cnst_get_mode_num(list_idx, 1, 'a', state, pbuf, mode_num)

      do isw = 1, nswbands
         savaervis = (isw .eq. idx_sw_diag)
         savaeruv  = (isw .eq. idx_uv_diag)
         savaernir = (isw .eq. idx_nir_diag)

         do k = top_lev, pver

            ! form bulk refractive index
            crefin(:ncol) = (0._r8, 0._r8)
            dryvol(:ncol) = 0._r8
            dustvol(:ncol) = 0._r8

            scatdust(:ncol)     = 0._r8
            absdust(:ncol)      = 0._r8
            hygrodust(:ncol)    = 0._r8
            scatso4(:ncol)      = 0._r8
            absso4(:ncol)       = 0._r8
            hygroso4(:ncol)     = 0._r8
            scatbc(:ncol)       = 0._r8
            absbc(:ncol)        = 0._r8
            hygrobc(:ncol)      = 0._r8
            scatpom(:ncol)      = 0._r8
            abspom(:ncol)       = 0._r8
            hygropom(:ncol)     = 0._r8
            scatsoa(:ncol)      = 0._r8
            abssoa(:ncol)       = 0._r8
            hygrosoa(:ncol)     = 0._r8
            scatseasalt(:ncol)  = 0._r8
            absseasalt(:ncol)   = 0._r8
            hygroseasalt(:ncol) = 0._r8

            crefin_sol(:ncol) = (0._r8, 0._r8)
            dryvol_sol(:ncol) = 0._r8
            crefin_insol(:ncol) = (0._r8, 0._r8)
            dryvol_insol(:ncol) = 0._r8
            dryvol_ss(:ncol) = 0._r8
            dryvol_sul(:ncol) = 0._r8
            dryvol_dust(:ncol) = 0._r8
            dryvol_poa(:ncol) = 0._r8
            dryvol_soa(:ncol) = 0._r8
            dryvol_mode(:ncol) = 0._r8
            specmmr_m(:ncol) = 0._r8
            specdens_m1         = 0._r8
            specmmr_nobc(:ncol) =0._r8

            ! aerosol species loop
            do l = 1, nspec
               call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, specmmr)
               call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                           refindex_aer_sw=specrefindex, spectype=spectype, &
                                           hygro_aer=hygro_aer)

               do i = 1, ncol
                  vol(i)      = specmmr(i,k)/specdens
                  dryvol(i)   = dryvol(i) + vol(i)
                  crefin(i)   = crefin(i) + vol(i)*specrefindex(isw)
               end do

               if (m==1 .and. savaervis .and. trim(spectype) == 'black-c') then
                 do i = 1, ncol
                   mass_bc(i,k) = specmmr(i,k)
                 enddo
               endif
               if (m==4 .and. savaervis .and. trim(spectype) == 'black-c') then
                 do i = 1, ncol
                   mass_bc4(i,k) = specmmr(i,k)
                 enddo
               endif

               if (trim(spectype) == 'black-c') then
                 do i = 1, ncol
                   vol_insol(i) = specmmr(i,k)/specdens
                   dryvol_insol(i)   = dryvol_insol(i) + vol_insol(i)
                   crefin_insol(i)   = crefin_insol(i) + vol_insol(i)*specrefindex(isw)
                 enddo
               else
                 do i = 1, ncol
                   vol_sol(i) = specmmr(i,k)/specdens
                   dryvol_sol(i)   = dryvol_sol(i) + vol_sol(i)
                   crefin_sol(i)   = crefin_sol(i) + vol_sol(i)*specrefindex(isw)
                   specmmr_nobc(i) = specmmr_nobc(i) + specmmr(i,k)
                 enddo
               endif

               if(m==1 .and. savaervis) then
                 do i=1, ncol
                  vol_m(i)      = specmmr(i,k)/specdens
                  dryvol_mode(i) = dryvol_mode(i) + vol_m(i)
                  specmmr_m(i) = specmmr_m(i)+vol_m(i)*specdens
                 enddo
               endif


               ! compute some diagnostics for visible band only
               if (savaervis) then

                  specrefr = real(specrefindex(isw))
                  specrefi = aimag(specrefindex(isw))

                  do i = 1, ncol
                     burden(i) = burden(i) + specmmr(i,k)*mass(i,k)
                  end do

                  if (trim(spectype) == 'dust') then
                     do i = 1, ncol
                        vol_dust(i) = specmmr(i,k)/specdens
                        dryvol_dust(i) = dryvol_dust(i) + vol_dust(i)
                        burdendust(i) = burdendust(i) + specmmr(i,k)*mass(i,k)
                        dustvol(i)    = vol(i)
                        scatdust(i)   = vol(i)*specrefr
                        absdust(i)    = -vol(i)*specrefi
                        hygrodust(i)  = vol(i)*hygro_aer
                     end do
                  end if

                  if (trim(spectype) == 'sulfate') then
                     do i = 1, ncol
                        vol_sul(i) = specmmr(i,k)/specdens
                        dryvol_sul(i) = dryvol_sul(i) + vol_sul(i)
                        burdenso4(i) = burdenso4(i) + specmmr(i,k)*mass(i,k)
                        scatso4(i)   = vol(i)*specrefr
                        absso4(i)    = -vol(i)*specrefi
                        hygroso4(i)  = vol(i)*hygro_aer
                     end do
                  end if
                  if (trim(spectype) == 'black-c') then
                     do i = 1, ncol
                        burdenbc(i) = burdenbc(i) + specmmr(i,k)*mass(i,k)
                        scatbc(i)   = vol(i)*specrefr
                        absbc(i)    = -vol(i)*specrefi
                        hygrobc(i)  = vol(i)*hygro_aer
                   end do
                  end if
                  if (trim(spectype) == 'p-organic') then
                     do i = 1, ncol
                        vol_poa(i) = specmmr(i,k)/specdens
                        dryvol_poa(i) = dryvol_poa(i) + vol_poa(i)
                        burdenpom(i) = burdenpom(i) + specmmr(i,k)*mass(i,k)
                        scatpom(i)   = vol(i)*specrefr
                        abspom(i)    = -vol(i)*specrefi
                        hygropom(i)  = vol(i)*hygro_aer
                      end do
                  end if
                  if (trim(spectype) == 's-organic') then
                     do i = 1, ncol
                        vol_soa(i) = specmmr(i,k)/specdens
                        dryvol_soa(i) = dryvol_soa(i) + vol_soa(i)
                        burdensoa(i) = burdensoa(i) + specmmr(i,k)*mass(i,k)
                        scatsoa(i)   = vol(i)*specrefr
                        abssoa(i)    = -vol(i)*specrefi
                        hygrosoa(i)  = vol(i)*hygro_aer
                     end do
                  end if
                  if (trim(spectype) == 'seasalt') then
                     do i = 1, ncol
                        vol_ss(i) = specmmr(i,k)/specdens
                        dryvol_ss(i) = dryvol_ss(i) + vol_ss(i)
                        burdenseasalt(i) = burdenseasalt(i) + specmmr(i,k)*mass(i,k)
                        scatseasalt(i)   = vol(i)*specrefr
                        absseasalt(i)    = -vol(i)*specrefi
                        hygroseasalt(i)  = vol(i)*hygro_aer
                      end do
                  end if

               end if
            end do ! species loop

            if (m==1 .and. savaervis) then
              do i=1,ncol
                specdens_m(i,k)=specmmr_m(i)/dryvol_mode(i)
                specdens_mwet(i,k)=(specmmr_m(i)+qaerwat(i,k))/(dryvol_mode(i)+qaerwat(i,k)/rhoh2o)
                specdens_nobc(i,k) = (specmmr_nobc(i) +qaerwat(i,k))/(dryvol_sol(i)+qaerwat(i,k)/rhoh2o) 
                mass_water(i,k) = qaerwat(i,k)
              enddo
            endif

            do i = 1, ncol
               watervol(i) = qaerwat(i,k)/rhoh2o
               wetvol(i) = watervol(i) + dryvol(i)
               if (watervol(i) < 0._r8) then
                  if (abs(watervol(i)) .gt. 1.e-1_r8*wetvol(i)) then
                     write(iulog,'(a,2e10.2,a)') 'watervol,wetvol=', &
                        watervol(i), wetvol(i), ' in '//subname
                  end if
                  watervol(i) = 0._r8
                  wetvol(i) = dryvol(i)
               end if
               !-----------------
               crefin_sol(i) = crefin_sol(i) + watervol(i)*crefwsw(isw)
               crefin_sol(i) = crefin_sol(i)/max(watervol(i)+ dryvol_sol(i),1.e-60_r8)
               if (list_idx==0)then
                crefin_insol(i)= crefin_insol(i)/dryvol_insol(i)
               endif
               if (list_idx>0)then
                crefin_insol(i)=(1._r8, 0._r8)
               endif
               if ( (m==1).and. savaervis) then
                 vol_frass(i,k)=dryvol_ss(i)
                 vol_fradust(i,k)=dryvol_dust(i)
                 vol_frapoa(i,k)=dryvol_poa(i)
                 vol_frasoa(i,k)=dryvol_soa(i)
                 vol_frasul(i,k)=dryvol_sul(i)
                 vol_fra(i,k)=dryvol_insol(i)
                 vol_frawater(i,k) = watervol(i)
                 vol_all(i,k) = dryvol(i)+watervol(i)
                 dens_wet(i,k) = (vol_frass(i,k)*1900 + vol_fradust(i,k)*2600 + vol_frapoa(i,k)*1000 + &
                                vol_frasoa(i,k)*1000 + vol_frasul(i,k)*1770 + &
                                vol_fra(i,k)*1700 + vol_frawater(i,k)*1000)/ vol_all(i,k)
                 refr_shell(i,k)=real(crefin_sol(i))
                 refi_shell(i,k)=abs(aimag(crefin_sol(i)))
                 refr_core(i,k)=real(crefin_insol(i))
                 refi_core(i,k)=abs(aimag(crefin_insol(i)))
                 ! optical module input 
                 vol_fra_input(i)=dryvol_insol(i)
                 refr_shell_input(i)=real(crefin_sol(i))
                 refi_shell_input(i)=abs(aimag(crefin_sol(i)))
                 dens_wet_input(i) = (vol_frass(i,k)*1900 + vol_fradust(i,k)*2600 + vol_frapoa(i,k)*1000 + &
                                vol_frasoa(i,k)*1000 + vol_frasul(i,k)*1770 + &
                                vol_fra(i,k)*1700 + vol_frawater(i,k)*1000)/ vol_all(i,k)
               endif
               if (m==4.and. savaervis) then
                  vol_fram4(i,k)=dryvol_insol(i)/(dryvol(i)+watervol(i))
               endif
               !-----

               ! volume mixing
               crefin(i) = crefin(i) + watervol(i)*crefwsw(isw)
               crefin(i) = crefin(i)/max(wetvol(i),1.e-60_r8)
               refr(i)   = real(crefin(i))
               refi(i)   = abs(aimag(crefin(i)))
            end do

            if (m==1 .and. savaervis) then
              do i=1,ncol
                refr1(i,k)=refr(i)
                refi1(i,k)=refi(i)
              enddo
            endif
            if (savaervis) then
              do i=1,ncol
                wetvols(i,k)=wetvol(i)
              enddo
            endif
          !------------ core-shell optical calculation ----------------------

             itab(:ncol) = 0
             call binterp(extpsw(:,:,:,isw), ncol, ncoef, prefr, prefi, &
                          refr, refi, refrtabsw(:,isw), refitabsw(:,isw), &
                           itab, jtab, ttab, utab, cext)
             call binterp(abspsw(:,:,:,isw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtabsw(:,isw), refitabsw(:,isw), &
                         itab, jtab, ttab, utab, cabs)
            call binterp(asmpsw(:,:,:,isw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtabsw(:,isw), refitabsw(:,isw), &
                         itab, jtab, ttab, utab, casm)

          !--------------------------
          if (savaervis .and. (m==1)  ) then ! remove BC 
          !-------------------------
          do i=1,ncol 
            if (logradsurf(i,k) .le. xrmax) then
              call opt(refr_shell_input(i),refi_shell_input(i),dens_wet_input(i),Radsurf(i,k),vol_fra_input(i),&
                                  mode_num(i,k),pabs(i),pext(i),pasm(i))
            else
               pext(i) = 1.5_r8/(radsurf(i,k)*rhoh2o) ! geometric optics
            endif

!            pasm(i) = 0.5_r8*casm(i,1)
!            do nc = 2, ncoef
!               pasm(i) = pasm(i) + cheb(nc,i,k)*casm(i,nc) 
!            enddo
            specpext(i)=pext(i)
            pext(i) = pext(i)*wetvol(i)*rhoh2o
            ! make sure that pext >=0
            pext(i) = max(0._r8,pext(i))
            specpabs(i) = pabs(i)
            pabs(i) = pabs(i)*wetvol(i)*rhoh2o
            pabs(i) = max(0._r8,pabs(i))
            pabs(i) = min(pext(i),pabs(i))
            palb(i) = 1._r8-pabs(i)/max(pext(i),1.e-40_r8)
            palb(i) = 1._r8-pabs(i)/max(pext(i),1.e-40_r8)
            dopaer(i) = pext(i)*mass(i,k)
            if (m==1 .and. savaervis) then
               palb_layer(i,k)=palb(i)
               tau_layer(i,k)=dopaer(i)
               asym_layer(i,k)=pasm(i)
            endif

           enddo ! do i=1,ncol
        !--------------------------------- 
          else  ! if m is not equal  1
        !----------------------------------


            ! call t_stopf('binterp')

            ! parameterized optical properties
            do i=1,ncol

               if (logradsurf(i,k) .le. xrmax) then
                  pext(i) = 0.5_r8*cext(i,1)
                  do nc = 2, ncoef
                     pext(i) = pext(i) + cheb(nc,i,k)*cext(i,nc)
                  enddo
                  pext(i) = exp(pext(i))
               else
                  pext(i) = 1.5_r8/(radsurf(i,k)*rhoh2o) ! geometric optics
               endif

               ! convert from m2/kg water to m2/kg aerosol
               specpext(i) = pext(i)
               pext(i) = pext(i)*wetvol(i)*rhoh2o
               pabs(i) = 0.5_r8*cabs(i,1)
               pasm(i) = 0.5_r8*casm(i,1)
               do nc = 2, ncoef
                  pabs(i) = pabs(i) + cheb(nc,i,k)*cabs(i,nc)
                  pasm(i) = pasm(i) + cheb(nc,i,k)*casm(i,nc)
               enddo
               specpabs(i) = pabs(i)
               pabs(i) = pabs(i)*wetvol(i)*rhoh2o
               pabs(i) = max(0._r8,pabs(i))
               pabs(i) = min(pext(i),pabs(i))

               palb(i) = 1._r8-pabs(i)/max(pext(i),1.e-40_r8)
               palb(i) = 1._r8-pabs(i)/max(pext(i),1.e-40_r8)

               dopaer(i) = pext(i)*mass(i,k)
            end do
        !-----------------------  
        endif ! if m is not equal 1
        !-------------------------


            if (savaervis) then
              do i=1,ncol
                specpext_vis(i,k)= specpext(i)
                specpabs_vis(i,k)= specpabs(i)
              enddo
            endif

            if (savaeruv) then
               do i = 1, ncol
                 extinctuv(i,k) = extinctuv(i,k) + dopaer(i)*air_density(i,k)/mass(i,k)
                 aoduv(i) = aoduv(i) + dopaer(i)
                  if (k.le.troplevchem(i)) then
                    aoduvst(i) = aoduvst(i) + dopaer(i)
                  end if
               end do
            end if

            if (savaernir) then
               do i = 1, ncol
                  extinctnir(i,k) = extinctnir(i,k) + dopaer(i)*air_density(i,k)/mass(i,k)
                  aodnir(i) = aodnir(i) + dopaer(i)
                  if (k.le.troplevchem(i)) then
                    aodnirst(i) = aodnirst(i) + dopaer(i)
                  end if
               end do
            endif

            ! Save aerosol optical depth at longest visible wavelength
            ! sum over layers
            if (savaervis) then
               ! aerosol extinction (/m)
               do i = 1, ncol
                  extinct(i,k) = extinct(i,k) + dopaer(i)*air_density(i,k)/mass(i,k)
                  absorb(i,k)  = absorb(i,k) + pabs(i)*air_density(i,k)
                  aodvis(i)    = aodvis(i) + dopaer(i)
                  aodabs(i)    = aodabs(i) + pabs(i)*mass(i,k)
                  aodmode(i)   = aodmode(i) + dopaer(i)
                  ssavis(i)    = ssavis(i) + dopaer(i)*palb(i)
                  asymvis(i)    = asymvis(i) + dopaer(i)*pasm(i)
                  asymext(i,k)  = asymext(i,k) + dopaer(i)*pasm(i)*air_density(i,k)/mass(i,k)
                  if (k.le.troplevchem(i)) then
                    aodvisst(i) = aodvisst(i) + dopaer(i)
                  end if

                  if (m==1) then
                     aodabsmode1(i) = aodabsmode1(i) +pabs(i)*mass(i,k)
                  endif

                  if (wetvol(i) > 1.e-40_r8) then

                     dustaodmode(i) = dustaodmode(i) + dopaer(i)*dustvol(i)/wetvol(i)

                     ! partition optical depth into contributions from each constituent
                     ! assume contribution is proportional to refractive index X volume

                     scath2o        = watervol(i)*real(crefwsw(isw))
		     absh2o         = -watervol(i)*aimag(crefwsw(isw))
		     sumscat        = scatso4(i) + scatpom(i) + scatsoa(i) + scatbc(i) + &
                                      scatdust(i) + scatseasalt(i) + scath2o
		     sumabs         = absso4(i) + abspom(i) + abssoa(i) + absbc(i) + &
                                      absdust(i) + absseasalt(i) + absh2o
                     sumhygro       = hygroso4(i) + hygropom(i) + hygrosoa(i) + hygrobc(i) + &
                                      hygrodust(i) + hygroseasalt(i)

                     scatdust(i)    = (scatdust(i) + scath2o*hygrodust(i)/sumhygro)/sumscat
                     absdust(i)     = (absdust(i) + absh2o*hygrodust(i)/sumhygro)/sumabs

                     scatso4(i)     = (scatso4(i) + scath2o*hygroso4(i)/sumhygro)/sumscat
                     absso4(i)      = (absso4(i) + absh2o*hygroso4(i)/sumhygro)/sumabs

                     scatpom(i)     = (scatpom(i) + scath2o*hygropom(i)/sumhygro)/sumscat
                     abspom(i)      = (abspom(i) + absh2o*hygropom(i)/sumhygro)/sumabs

                     scatsoa(i)     = (scatsoa(i) + scath2o*hygrosoa(i)/sumhygro)/sumscat
                     abssoa(i)      = (abssoa(i) + absh2o*hygrosoa(i)/sumhygro)/sumabs

                     scatbc(i)      = (scatbc(i) + scath2o*hygrobc(i)/sumhygro)/sumscat
                     absbc(i)       = (absbc(i) + absh2o*hygrobc(i)/sumhygro)/sumabs

                     scatseasalt(i) = (scatseasalt(i) + scath2o*hygroseasalt(i)/sumhygro)/sumscat
                     absseasalt(i)  = (absseasalt(i) + absh2o*hygroseasalt(i)/sumhygro)/sumabs
                     
                     aodabsbc(i)    = aodabsbc(i) + absbc(i)*dopaer(i)*(1.0_r8-palb(i))
                     if (m==1) then
                       absorbBC(i,k) = absorbBC(i,k) + absbc(i)*pabs(i)*air_density(i,k)
                     endif
                     if (m==4) then
                       absorbBC4(i,k) = absorbBC4(i,k) + absbc(i)*pabs(i)*air_density(i,k)
                     endif

                     aodc           = (absdust(i)*(1.0_r8 - palb(i)) + palb(i)*scatdust(i))*dopaer(i)
                     dustaod(i)     = dustaod(i) + aodc

                     aodc           = (absso4(i)*(1.0_r8 - palb(i)) + palb(i)*scatso4(i))*dopaer(i)
                     so4aod(i)      = so4aod(i) + aodc

                     aodc           = (abspom(i)*(1.0_r8 - palb(i)) + palb(i)*scatpom(i))*dopaer(i)
                     pomaod(i)      = pomaod(i) + aodc

                     aodc           = (abssoa(i)*(1.0_r8 - palb(i)) + palb(i)*scatsoa(i))*dopaer(i)
                     soaaod(i)      = soaaod(i) + aodc

                     aodc           = (absbc(i)*(1.0_r8 - palb(i)) + palb(i)*scatbc(i))*dopaer(i)
                     bcaod(i)       = bcaod(i) + aodc

                     aodc           = (absseasalt(i)*(1.0_r8 - palb(i)) + palb(i)*scatseasalt(i))*dopaer(i)
                     seasaltaod(i)  = seasaltaod(i) + aodc

                  endif

               end do
            endif

            do i = 1, ncol

               if ((dopaer(i) <= -1.e-10_r8) .or. (dopaer(i) >= 30._r8)) then

                  if (dopaer(i) <= -1.e-10_r8) then
                     write(iulog,*) "ERROR: Negative aerosol optical depth &
                          &in this layer."
                  else
                     write(iulog,*) "WARNING: Aerosol optical depth is &
                          &unreasonably high in this layer."
                  end if

                  write(iulog,*) 'dopaer(', i, ',', k, ',', m, ',', lchnk, ')=', dopaer(i)
                  ! write(iulog,*) 'itab,jtab,ttab,utab=',itab(i),jtab(i),ttab(i),utab(i)
                  write(iulog,*) 'k=', k, ' pext=', pext(i), ' specext=', specpext(i)
                  write(iulog,*) 'wetvol=', wetvol(i), ' dryvol=', dryvol(i), ' watervol=', watervol(i)
                  ! write(iulog,*) 'cext=',(cext(i,l),l=1,ncoef)
                  ! write(iulog,*) 'crefin=',crefin(i)
                  write(iulog,*) 'nspec=', nspec
                  ! write(iulog,*) 'cheb=', (cheb(nc,m,i,k),nc=2,ncoef)
                  do l = 1, nspec
                     call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, specmmr)
                     call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                                 refindex_aer_sw=specrefindex)
                     volf = specmmr(i,k)/specdens
                     write(iulog,*) 'l=', l, 'vol(l)=', volf
                     write(iulog,*) 'isw=', isw, 'specrefindex(isw)=', specrefindex(isw)
                     write(iulog,*) 'specdens=', specdens
                  end do

                  nerr_dopaer = nerr_dopaer + 1
!                  if (nerr_dopaer >= nerrmax_dopaer) then
                  if (dopaer(i) < -1.e-10_r8) then
                     write(iulog,*) '*** halting in '//subname//' after nerr_dopaer =', nerr_dopaer
                     call endrun('exit from '//subname)
                  end if

               end if
            end do

            do i=1,ncol
                
               tauxar(i,k,isw) = tauxar(i,k,isw) + dopaer(i)
               wa(i,k,isw)     = wa(i,k,isw)     + dopaer(i)*palb(i)
               ga(i,k,isw)     = ga(i,k,isw)     + dopaer(i)*palb(i)*pasm(i)
               fa(i,k,isw)     = fa(i,k,isw)     + dopaer(i)*palb(i)*pasm(i)*pasm(i)
! test sensitivity of AOD, SSA and asymmetry on BC DRF
!               if (list_idx==0) then
!                tauxar(i,k,isw)=0.1
!                wa(i,k,isw)=0.1*0.95
!                ga(i,k,isw)=0.1*0.95*0.78
!                fa(i,k,isw)=0.1*0.95*0.78*0.78
!               else
!                tauxar(i,k,isw)=0.1*0.5
!                wa(i,k,isw)=0.1*0.95*0.5
!                ga(i,k,isw)=0.1*0.95*0.78*0.5
!                fa(i,k,isw)=0.1*0.95*0.78*0.78*0.5
!               endif 
            end do

         end do ! pver

!                 if (list_idx > 0 .and. m==1 .and. savaervis) then
!                    open(40,file='N0_BC.txt')
!                    open(50,file='refr_core.txt')
!                  do k = top_lev, pver
!                    do i =1, ncol
!                      write(40,*)N0_BC(i,k)
!                      write(50,*)refr_core(i,k)
!                    enddo
!                  enddo
!                 endif
!                 close(40)
!                 close(50)

      end do ! sw bands
         ! output without BC in diagnostic calculation

      ! mode diagnostics
      ! The diagnostics are currently only output for the climate list.  Code mods will
      ! be necessary to provide output for the rad_diag lists.
      if (list_idx == 0) then

         write(outname,'(a,i1)') 'BURDENdn', m
         call outfld(trim(outname), burden, pcols, lchnk)

         write(outname,'(a,i1)') 'AODdnMODE', m
         call outfld(trim(outname), aodmode, pcols, lchnk)

         write(outname,'(a,i1)') 'AODdnDUST', m
         call outfld(trim(outname), dustaodmode, pcols, lchnk)
         
         do i = 1, nnite
            burden(idxnite(i))  = fillvalue
            aodmode(idxnite(i)) = fillvalue
            dustaodmode(idxnite(i)) = fillvalue
         end do

         write(outname,'(a,i1)') 'BURDEN', m
         call outfld(trim(outname), burden, pcols, lchnk)

         write(outname,'(a,i1)') 'AODMODE', m
         call outfld(trim(outname), aodmode, pcols, lchnk)

         write(outname,'(a,i1)') 'AODDUST', m
         call outfld(trim(outname), dustaodmode, pcols, lchnk)

         write(outname,'(a,i1)') 'radsurf', m
         call outfld(trim(outname), radsurf, pcols, lchnk)

         write(outname,'(a,i1)') 'specpabs', m
         call outfld(trim(outname), specpabs_vis, pcols, lchnk)

         write(outname,'(a,i1)') 'specpext', m
         call outfld(trim(outname), specpext_vis, pcols, lchnk)

         write(outname,'(a,i1)') 'wetvols',m
         call outfld(trim(outname),wetvols, pcols, lchnk)

      end if

   end do ! nmodes

  ! remove BC component (keep Dp constant)
   if (list_idx > 0) then
      deallocate(dgnumdry_m)
      deallocate(dgnumwet_m)
      deallocate(qaerwat_m)
      deallocate(wetdens_m)
      deallocate(hygro_m)
      deallocate(dryvol_m)
      deallocate(dryrad_m)
      deallocate(drymass_m)
      deallocate(so4dryvol_m)
      deallocate(naer_m)
   end if

   ! Output visible band diagnostics for quantities summed over the modes
   ! These fields are put out for diagnostic lists as well as the climate list.

   call outfld('EXTINCTdn'//diag(list_idx),  extinct, pcols, lchnk)
   call outfld('ABSORBdn'//diag(list_idx),   absorb,  pcols, lchnk)
   call outfld('AODVISdn'//diag(list_idx),   aodvis,  pcols, lchnk)
   call outfld('AODABSdn'//diag(list_idx),   aodabs,  pcols, lchnk)
   call outfld('AODVISstdn'//diag(list_idx), aodvisst,pcols, lchnk)
   call outfld('EXTxASYMdn'//diag(list_idx), asymext, pcols, lchnk)
   
   do i = 1, nnite
      extinct(idxnite(i),:) = fillvalue
      absorb(idxnite(i),:)  = fillvalue
      aodvis(idxnite(i))    = fillvalue
      aodabs(idxnite(i))    = fillvalue
      aodvisst(idxnite(i))  = fillvalue
      asymext(idxnite(i),:) = fillvalue
      refr1(idxnite(i),:)     = fillvalue
      refi1(idxnite(i),:)     = fillvalue
      specdens_m(idxnite(i),:) = fillvalue
      mass(idxnite(i),:)  = fillvalue
      air_density(idxnite(i),:)  = fillvalue
      vol_fra(idxnite(i),:) = fillvalue
      refr_shell(idxnite(i),:) = fillvalue
      refi_shell(idxnite(i),:) = fillvalue
      vol_fradust(idxnite(i),:) = fillvalue
      vol_frasoa(idxnite(i),:) = fillvalue
      vol_frapoa(idxnite(i),:) = fillvalue
      vol_frasul(idxnite(i),:) = fillvalue
      vol_frass(idxnite(i),:) = fillvalue
      vol_fram4(idxnite(i),:) = fillvalue
      refr_core(idxnite(i),:) = fillvalue
      refi_core(idxnite(i),:) = fillvalue
      aodabsmode1(idxnite(i))    = fillvalue
      specdens_mwet(idxnite(i),:) = fillvalue
      specdens_nobc(idxnite(i),:) = fillvalue
      mass_water(idxnite(i),:) = fillvalue
      absorbBC(idxnite(i),:)  = fillvalue
      vol_all(idxnite(i),:) = fillvalue
      vol_frawater(idxnite(i),:) = fillvalue
      mass_bc(idxnite(i),:) = fillvalue
      absorbBC4(idxnite(i),:)  = fillvalue
      mass_bc4(idxnite(i),:) = fillvalue
!      pabs_bc(idxnite(i),:)  = fillvalue
        palb_layer(idxnite(i),:) = fillvalue
        tau_layer(idxnite(i),:) = fillvalue
        asym_layer(idxnite(i),:) = fillvalue

   end do

   call outfld('EXTINCT'//diag(list_idx),  extinct, pcols, lchnk)
   call outfld('ABSORB'//diag(list_idx),   absorb,  pcols, lchnk)
   call outfld('AODVIS'//diag(list_idx),   aodvis,  pcols, lchnk)
   call outfld('AODABS'//diag(list_idx),   aodabs,  pcols, lchnk)
!   call outfld('AODVISst'//diag(list_idx), aodvisst,pcols, lchnk)
   call outfld('EXTxASYM'//diag(list_idx), asymext, pcols, lchnk)
      call outfld('palb_layer'//diag(list_idx), palb_layer,    pcols, lchnk)
      call outfld('tau_layer'//diag(list_idx),         tau_layer,    pcols, lchnk)
      call outfld('asym_layer'//diag(list_idx),         asym_layer,    pcols, lchnk)
!   call outfld('refr'//diag(list_idx),   refr1,  pcols, lchnk)
!   call outfld('refi'//diag(list_idx),   refi1,  pcols, lchnk)
!   call outfld('num'//diag(list_idx),   mode_num,  pcols, lchnk)
!   call outfld('specdens_m'//diag(list_idx),  specdens_m, pcols, lchnk)
!   call outfld('mass'//diag(list_idx),   mass,  pcols, lchnk)
!   call outfld('air_density'//diag(list_idx),   air_density,  pcols, lchnk)
!   call outfld('vol_fra'//diag(list_idx),  vol_fra, pcols, lchnk)
!   call outfld('refr_shell'//diag(list_idx),  refr_shell, pcols, lchnk)
!   call outfld('refi_shell'//diag(list_idx),  refi_shell, pcols, lchnk)
!   call outfld('vol_frass'//diag(list_idx),  vol_frass, pcols, lchnk)
!   call outfld('vol_fradust'//diag(list_idx),  vol_fradust, pcols, lchnk)
!   call outfld('vol_frasoa'//diag(list_idx),  vol_frasoa, pcols, lchnk)
!   call outfld('vol_frapoa'//diag(list_idx),  vol_frapoa, pcols, lchnk)
!   call outfld('vol_frasul'//diag(list_idx),  vol_frasul, pcols, lchnk)
!   call outfld('vol_fram4'//diag(list_idx),  vol_fram4, pcols, lchnk)
!   call outfld('refr_core'//diag(list_idx),  refr_core, pcols, lchnk)
!   call outfld('refi_core'//diag(list_idx),  refi_core, pcols, lchnk)
!   call outfld('aodabsmode1'//diag(list_idx),   aodabsmode1,  pcols, lchnk)
!   call outfld('specdens_mwet'//diag(list_idx),  specdens_mwet, pcols, lchnk)
!   call outfld('specdens_nobc'//diag(list_idx), specdens_nobc, pcols, lchnk)
!   call outfld('mass_water'//diag(list_idx),  mass_water, pcols, lchnk)
!   call outfld('absorbBC'//diag(list_idx),   absorbBC,  pcols, lchnk)
!   call outfld('absorbBC4'//diag(list_idx),   absorbBC4,  pcols, lchnk)
!   call outfld('vol_all'//diag(list_idx),  vol_all, pcols, lchnk)
!   call outfld('vol_frawater'//diag(list_idx),  vol_frawater, pcols, lchnk)
!   call outfld('mass_bc'//diag(list_idx),  mass_bc, pcols, lchnk)
!   call outfld('pabs_bc'//diag(list_idx),   pabs_bc,  pcols, lchnk)
!   call outfld('mass_bc4'//diag(list_idx),  mass_bc4, pcols, lchnk)

   ! These diagnostics are output only for climate list
   if (list_idx == 0) then
      do i = 1, ncol
         if (aodvis(i) > 1.e-10_r8) then
            ssavis(i) = ssavis(i)/aodvis(i)
         else
            ssavis(i) = 0.925_r8
         endif
      end do
      
      call outfld('SSAVISdn',        ssavis,        pcols, lchnk)
      call outfld('AODxASYMdn',      asymvis,       pcols, lchnk)

      call outfld('EXTINCTUVdn',     extinctuv,     pcols, lchnk)
      call outfld('EXTINCTNIRdn',    extinctnir,    pcols, lchnk)
      call outfld('AODUVdn',         aoduv,         pcols, lchnk)
      call outfld('AODNIRdn',        aodnir,        pcols, lchnk)
      call outfld('AODUVstdn',       aoduvst,       pcols, lchnk)
      call outfld('AODNIRstdn',      aodnirst,      pcols, lchnk)

      call outfld('BURDENDUSTdn',    burdendust,    pcols, lchnk)
      call outfld('BURDENSO4dn' ,    burdenso4,     pcols, lchnk)
      call outfld('BURDENPOMdn' ,    burdenpom,     pcols, lchnk)
      call outfld('BURDENSOAdn' ,    burdensoa,     pcols, lchnk)
      call outfld('BURDENBCdn'  ,    burdenbc,      pcols, lchnk)
      call outfld('BURDENSEASALTdn', burdenseasalt, pcols, lchnk)

      call outfld('AODABSBCdn',      aodabsbc,      pcols, lchnk)

      call outfld('AODDUSTdn',       dustaod,       pcols, lchnk)
      call outfld('AODSO4dn',        so4aod,        pcols, lchnk)
      call outfld('AODPOMdn',        pomaod,        pcols, lchnk)
      call outfld('AODSOAdn',        soaaod,        pcols, lchnk)
      call outfld('AODBCdn',         bcaod,         pcols, lchnk)
      call outfld('AODSSdn',         seasaltaod,    pcols, lchnk)


      do i = 1, nnite
         ssavis(idxnite(i))     = fillvalue
         asymvis(idxnite(i))    = fillvalue

         aoduv(idxnite(i))      = fillvalue
         aodnir(idxnite(i))     = fillvalue
         aoduvst(idxnite(i))    = fillvalue
         aodnirst(idxnite(i))   = fillvalue
         extinctuv(idxnite(i),:)  = fillvalue
         extinctnir(idxnite(i),:) = fillvalue

         burdendust(idxnite(i)) = fillvalue
         burdenso4(idxnite(i))  = fillvalue
         burdenpom(idxnite(i))  = fillvalue
         burdensoa(idxnite(i))  = fillvalue
         burdenbc(idxnite(i))   = fillvalue
         burdenseasalt(idxnite(i)) = fillvalue

         aodabsbc(idxnite(i))   = fillvalue

         dustaod(idxnite(i))    = fillvalue
         so4aod(idxnite(i))     = fillvalue
         pomaod(idxnite(i))     = fillvalue
         soaaod(idxnite(i))     = fillvalue
         bcaod(idxnite(i))      = fillvalue
         seasaltaod(idxnite(i)) = fillvalue
       end do

      call outfld('SSAVIS',        ssavis,        pcols, lchnk)
      call outfld('AODxASYM',      asymvis,       pcols, lchnk)

      call outfld('EXTINCTUV',     extinctuv,     pcols, lchnk)
      call outfld('EXTINCTNIR',    extinctnir,    pcols, lchnk)
      call outfld('AODUV',         aoduv,         pcols, lchnk)
      call outfld('AODNIR',        aodnir,        pcols, lchnk)
      call outfld('AODUVst',       aoduvst,       pcols, lchnk)
      call outfld('AODNIRst',      aodnirst,      pcols, lchnk)

      call outfld('BURDENDUST',    burdendust,    pcols, lchnk)
      call outfld('BURDENSO4' ,    burdenso4,     pcols, lchnk)
      call outfld('BURDENPOM' ,    burdenpom,     pcols, lchnk)
      call outfld('BURDENSOA' ,    burdensoa,     pcols, lchnk)
      call outfld('BURDENBC'  ,    burdenbc,      pcols, lchnk)
      call outfld('BURDENSEASALT', burdenseasalt, pcols, lchnk)

      call outfld('AODABSBC',      aodabsbc,      pcols, lchnk)

      call outfld('AODDUST',       dustaod,       pcols, lchnk)
      call outfld('AODSO4',        so4aod,        pcols, lchnk)
      call outfld('AODPOM',        pomaod,        pcols, lchnk)
      call outfld('AODSOA',        soaaod,        pcols, lchnk)
      call outfld('AODBC',         bcaod,         pcols, lchnk)
      call outfld('AODSS',         seasaltaod,    pcols, lchnk)
   end if

end subroutine modal_aero_sw

!===============================================================================

subroutine modal_aero_lw(list_idx, state, pbuf, tauxar)

   ! calculates aerosol lw radiative properties

   integer,             intent(in)  :: list_idx ! index of the climate or a diagnostic list
   type(physics_state), intent(in), target :: state    ! state variables
   
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(out) :: tauxar(pcols,pver,nlwbands) ! layer absorption optical depth

   ! Local variables
   integer :: i, ifld, ilw, k, l, m, nc, ns
   integer :: lchnk                    ! chunk id
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nmodes
   integer :: nspec
   integer :: istat

   real(r8), pointer :: dgnumwet(:,:)  ! wet number mode diameter (m)
   real(r8), pointer :: qaerwat(:,:)   ! aerosol water (g/g)

   real(r8), pointer :: dgnumdry_m(:,:,:) ! number mode dry diameter for all modes
   real(r8), pointer :: dgnumwet_m(:,:,:) ! number mode wet diameter for all modes
   real(r8), pointer :: qaerwat_m(:,:,:)  ! aerosol water (g/g) for all modes
   real(r8), pointer :: wetdens_m(:,:,:)  ! 
   real(r8), pointer :: hygro_m(:,:,:)  !
   real(r8), pointer :: dryvol_m(:,:,:)  !
   real(r8), pointer :: dryrad_m(:,:,:)  !
   real(r8), pointer :: drymass_m(:,:,:)  !
   real(r8), pointer :: so4dryvol_m(:,:,:)  !
   real(r8), pointer :: naer_m(:,:,:)  !

   real(r8) :: sigma_logr_aer          ! geometric standard deviation of number distribution
   real(r8) :: alnsg_amode             ! log of geometric standard deviation of number distribution
   real(r8) :: xrad(pcols)
   real(r8) :: cheby(ncoef,pcols,pver)  ! chebychef polynomials

   real(r8) :: mass(pcols,pver) ! layer mass

   real(r8),    pointer :: specmmr(:,:)        ! species mass mixing ratio
   real(r8)             :: specdens            ! species density (kg/m3)
   complex(r8), pointer :: specrefindex(:)     ! species refractive index

   real(r8) :: vol(pcols)       ! volume concentration of aerosol specie (m3/kg)
   real(r8) :: dryvol(pcols)    ! volume concentration of aerosol mode (m3/kg)
   real(r8) :: wetvol(pcols)    ! volume concentration of wet mode (m3/kg)
   real(r8) :: watervol(pcols)  ! volume concentration of water in each mode (m3/kg)
   real(r8) :: refr(pcols)      ! real part of refractive index
   real(r8) :: refi(pcols)      ! imaginary part of refractive index
   complex(r8) :: crefin(pcols) ! complex refractive index
   real(r8), pointer :: refrtablw(:,:) ! table of real refractive indices for aerosols
   real(r8), pointer :: refitablw(:,:) ! table of imag refractive indices for aerosols
   real(r8), pointer :: absplw(:,:,:,:) ! specific absorption

   integer  :: itab(pcols), jtab(pcols)
   real(r8) :: ttab(pcols), utab(pcols)
   real(r8) :: cabs(pcols,ncoef)
   real(r8) :: pabs(pcols)      ! parameterized specific absorption (m2/kg)
   real(r8) :: dopaer(pcols)    ! aerosol optical depth in layer

   integer, parameter :: nerrmax_dopaer=1000
   integer  :: nerr_dopaer = 0
   real(r8) :: volf             ! volume fraction of insoluble aerosol

   character(len=*), parameter :: subname = 'modal_aero_lw'
   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   ! initialize output variables
   tauxar(:ncol,:,:) = 0._r8

   ! dry mass in each cell
   mass(:ncol,:) = state%pdeldry(:ncol,:)*rga

   ! loop over all aerosol modes
   call rad_cnst_get_info(list_idx, nmodes=nmodes)

   if (list_idx == 0) then
      ! water uptake and wet radius for the climate list has already been calculated
      call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet_m)
      call pbuf_get_field(pbuf, qaerwat_idx,  qaerwat_m)
   else
      ! If doing a diagnostic calculation then need to calculate the wet radius
      ! and water uptake for the diagnostic modes
      allocate(dgnumdry_m(pcols,pver,nmodes),  dgnumwet_m(pcols,pver,nmodes), &
               qaerwat_m(pcols,pver,nmodes),   wetdens_m(pcols,pver,nmodes), &
               hygro_m(pcols,pver,nmodes),     dryvol_m(pcols,pver,nmodes), &
               dryrad_m(pcols,pver,nmodes),    drymass_m(pcols,pver,nmodes),  &
               so4dryvol_m(pcols,pver,nmodes), naer_m(pcols,pver,nmodes),     stat=istat)

      if (istat > 0) then
         call endrun('modal_aero_lw: allocation FAILURE: arrays for diagnostic calcs')
      end if
      call modal_aero_calcsize_diag(state, pbuf, list_idx, dgnumdry_m, hygro_m, &
                                    dryvol_m, dryrad_m, drymass_m, so4dryvol_m, naer_m)  
      call modal_aero_wateruptake_dr(state, pbuf, list_idx, dgnumdry_m, dgnumwet_m, &
                                     qaerwat_m, wetdens_m,  hygro_m, dryvol_m, dryrad_m, &
                                     drymass_m, so4dryvol_m, naer_m)
   endif

   do m = 1, nmodes

      dgnumwet => dgnumwet_m(:,:,m)
      qaerwat  => qaerwat_m(:,:,m)

      ! get mode properties
      call rad_cnst_get_mode_props(list_idx, m, sigmag=sigma_logr_aer, refrtablw=refrtablw , &
         refitablw=refitablw, absplw=absplw)

      ! get mode info
      call rad_cnst_get_info(list_idx, m, nspec=nspec)

      ! calc size parameter for all columns
      ! this is the same calculation that's done in modal_size_parameters, but there
      ! some intermediate results are saved and the chebyshev polynomials are stored
      ! in a array with different index order.  Could be unified.
      do k = top_lev, pver
         do i = 1, ncol
            alnsg_amode = log( sigma_logr_aer )
            ! convert from number diameter to surface area
            xrad(i) = log(0.5_r8*dgnumwet(i,k)) + 2.0_r8*alnsg_amode*alnsg_amode
            ! normalize size parameter
            xrad(i) = max(xrad(i), xrmin)
            xrad(i) = min(xrad(i), xrmax)
            xrad(i) = (2*xrad(i)-xrmax-xrmin)/(xrmax-xrmin)
            ! chebyshev polynomials
            cheby(1,i,k) = 1.0_r8
            cheby(2,i,k) = xrad(i)
            do nc = 3, ncoef
               cheby(nc,i,k) = 2.0_r8*xrad(i)*cheby(nc-1,i,k)-cheby(nc-2,i,k)
            end do
         end do
      end do

      do ilw = 1, nlwbands

         do k = top_lev, pver

            ! form bulk refractive index. Use volume mixing for infrared
            crefin(:ncol) = (0._r8, 0._r8)
            dryvol(:ncol) = 0._r8

            ! aerosol species loop
            do l = 1, nspec
               call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, specmmr)
               call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                           refindex_aer_lw=specrefindex)

               do i = 1, ncol
                  vol(i)    = specmmr(i,k)/specdens
                  dryvol(i) = dryvol(i) + vol(i)
                  crefin(i) = crefin(i) + vol(i)*specrefindex(ilw)
               end do
            end do

            do i = 1, ncol
               watervol(i) = qaerwat(i,k)/rhoh2o
               wetvol(i)   = watervol(i) + dryvol(i)
               if (watervol(i) < 0.0_r8) then
                  if (abs(watervol(i)) .gt. 1.e-1_r8*wetvol(i)) then
                     write(iulog,*) 'watervol,wetvol,dryvol=',watervol(i),wetvol(i),dryvol(i),' in '//subname
                  end if
                  watervol(i) = 0._r8
                  wetvol(i)   = dryvol(i)
               end if

               crefin(i) = crefin(i) + watervol(i)*crefwlw(ilw)
               if (wetvol(i) > 1.e-40_r8) crefin(i) = crefin(i)/wetvol(i)
               refr(i) = real(crefin(i))
               refi(i) = aimag(crefin(i))
            end do

            ! interpolate coefficients linear in refractive index
            ! first call calcs itab,jtab,ttab,utab
            itab(:ncol) = 0
            call binterp(absplw(:,:,:,ilw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtablw(:,ilw), refitablw(:,ilw), &
                         itab, jtab, ttab, utab, cabs)

            ! parameterized optical properties
            do i = 1, ncol
               pabs(i) = 0.5_r8*cabs(i,1)
               do nc = 2, ncoef
                  pabs(i) = pabs(i) + cheby(nc,i,k)*cabs(i,nc)
               end do
               pabs(i)   = pabs(i)*wetvol(i)*rhoh2o
               pabs(i)   = max(0._r8,pabs(i))
               dopaer(i) = pabs(i)*mass(i,k)
            end do

            do i = 1, ncol

               if ((dopaer(i) <= -1.e-10_r8) .or. (dopaer(i) >= 20._r8)) then

                  if (dopaer(i) <= -1.e-10_r8) then
                     write(iulog,*) "ERROR: Negative aerosol optical depth &
                          &in this layer."
                  else
                     write(iulog,*) "WARNING: Aerosol optical depth is &
                          &unreasonably high in this layer."
                  end if

                  write(iulog,*) 'dopaer(',i,',',k,',',m,',',lchnk,')=', dopaer(i)
                  write(iulog,*) 'k=',k,' pabs=', pabs(i)
                  write(iulog,*) 'wetvol=',wetvol(i),' dryvol=',dryvol(i),     &
                     ' watervol=',watervol(i)
                  write(iulog,*) 'cabs=', (cabs(i,l),l=1,ncoef)
                  write(iulog,*) 'crefin=', crefin(i)
                  write(iulog,*) 'nspec=', nspec
                  do l = 1,nspec
                     call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, specmmr)
                     call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                                 refindex_aer_lw=specrefindex)
                     volf = specmmr(i,k)/specdens
                     write(iulog,*) 'l=',l,'vol(l)=',volf
                     write(iulog,*) 'ilw=',ilw,' specrefindex(ilw)=',specrefindex(ilw)
                     write(iulog,*) 'specdens=',specdens
                  end do

                  nerr_dopaer = nerr_dopaer + 1
                  if (nerr_dopaer >= nerrmax_dopaer .or. dopaer(i) < -1.e-10_r8) then
                     write(iulog,*) '*** halting in '//subname//' after nerr_dopaer =', nerr_dopaer
                     call endrun()
                  end if

               end if
            end do

            do i = 1, ncol
               tauxar(i,k,ilw) = tauxar(i,k,ilw) + dopaer(i)
            end do

         end do ! k = top_lev, pver

      end do  ! nlwbands

   end do ! m = 1, nmodes

   if (list_idx > 0) then
      deallocate(dgnumdry_m)
      deallocate(dgnumwet_m)
      deallocate(qaerwat_m)
      deallocate(wetdens_m)
      deallocate(hygro_m)
      deallocate(dryvol_m)
      deallocate(dryrad_m)
      deallocate(drymass_m)
      deallocate(so4dryvol_m)
      deallocate(naer_m)
   end if

end subroutine modal_aero_lw

!===============================================================================
! Private routines
!===============================================================================

subroutine read_water_refindex(infilename)

   ! read water refractive index file and set module data

   character*(*), intent(in) :: infilename   ! modal optics filename

   ! Local variables

   integer            :: i, ierr
   type(file_desc_t)  :: ncid              ! pio file handle
   integer            :: did               ! dimension ids
   integer            :: dimlen            ! dimension lengths
   type(var_desc_t)   :: vid               ! variable ids
   real(r8) :: refrwsw(nswbands), refiwsw(nswbands) ! real, imaginary ref index for water visible
   real(r8) :: refrwlw(nlwbands), refiwlw(nlwbands) ! real, imaginary ref index for water infrared
   !----------------------------------------------------------------------------

   ! open file
   call cam_pio_openfile(ncid, infilename, PIO_NOWRITE)

   ! inquire dimensions.  Check that file values match parameter values.

   ierr = pio_inq_dimid(ncid, 'lw_band', did)
   ierr = pio_inq_dimlen(ncid, did, dimlen)
   if (dimlen .ne. nlwbands) then
      write(iulog,*) 'lw_band len=', dimlen, ' from ', infilename, ' ne nlwbands=', nlwbands
      call endrun('read_modal_optics: bad lw_band value')
   endif

   ierr = pio_inq_dimid(ncid, 'sw_band', did)
   ierr = pio_inq_dimlen(ncid, did, dimlen)
   if (dimlen .ne. nswbands) then
      write(iulog,*) 'sw_band len=', dimlen, ' from ', infilename, ' ne nswbands=', nswbands
      call endrun('read_modal_optics: bad sw_band value')
   endif

   ! read variables
   ierr = pio_inq_varid(ncid, 'refindex_real_water_sw', vid)
   ierr = pio_get_var(ncid, vid, refrwsw)

   ierr = pio_inq_varid(ncid, 'refindex_im_water_sw', vid)
   ierr = pio_get_var(ncid, vid, refiwsw)

   ierr = pio_inq_varid(ncid, 'refindex_real_water_lw', vid)
   ierr = pio_get_var(ncid, vid, refrwlw)

   ierr = pio_inq_varid(ncid, 'refindex_im_water_lw', vid)
   ierr = pio_get_var(ncid, vid, refiwlw)

   ! set complex representation of refractive indices as module data
   do i = 1, nswbands
      crefwsw(i)  = cmplx(refrwsw(i), abs(refiwsw(i)),kind=r8)
   end do
   do i = 1, nlwbands
      crefwlw(i)  = cmplx(refrwlw(i), abs(refiwlw(i)),kind=r8)
   end do

   call pio_closefile(ncid)

end subroutine read_water_refindex

!===============================================================================

subroutine modal_size_parameters(ncol, sigma_logr_aer, dgnumwet, radsurf, logradsurf, cheb)

   integer,  intent(in)  :: ncol
   real(r8), intent(in)  :: sigma_logr_aer  ! geometric standard deviation of number distribution
   real(r8), intent(in)  :: dgnumwet(:,:)   ! aerosol wet number mode diameter (m)
   real(r8), intent(out) :: radsurf(:,:)    ! aerosol surface mode radius
   real(r8), intent(out) :: logradsurf(:,:) ! log(aerosol surface mode radius)
   real(r8), intent(out) :: cheb(:,:,:)

   integer  :: i, k, nc
   real(r8) :: alnsg_amode
   real(r8) :: explnsigma
   real(r8) :: xrad(pcols) ! normalized aerosol radius
   !-------------------------------------------------------------------------------

   alnsg_amode = log(sigma_logr_aer)
   explnsigma = exp(2.0_r8*alnsg_amode*alnsg_amode)

   do k = top_lev, pver
      do i = 1, ncol
         ! convert from number mode diameter to surface area
         radsurf(i,k) = 0.5_r8*dgnumwet(i,k)*explnsigma
         logradsurf(i,k) = log(radsurf(i,k))
         ! normalize size parameter
         xrad(i) = max(logradsurf(i,k),xrmin)
         xrad(i) = min(xrad(i),xrmax)
         xrad(i) = (2._r8*xrad(i)-xrmax-xrmin)/(xrmax-xrmin)
         ! chebyshev polynomials
         cheb(1,i,k) = 1._r8
         cheb(2,i,k) = xrad(i)
         do nc = 3, ncoef
            cheb(nc,i,k) = 2._r8*xrad(i)*cheb(nc-1,i,k)-cheb(nc-2,i,k)
         end do
      end do
   end do

end subroutine modal_size_parameters

!===============================================================================

    subroutine binterp(table,ncol,km,im,jm,x,y,xtab,ytab,ix,jy,t,u,out)

        !     bilinear interpolation of table
        !
        implicit none
        integer im,jm,km,ncol
        real(r8) table(km,im,jm),xtab(im),ytab(jm),out(pcols,km)
        integer i,ix(pcols),ip1,j,jy(pcols),jp1,k,ic,ip1m(pcols),jp1m(pcols),ixc,jyc
        real(r8) x(pcols),dx,t(pcols),y(pcols),dy,u(pcols),tu(pcols),tuc(pcols),tcu(pcols),tcuc(pcols)
        
        if(ix(1).gt.0) go to 30
        if(im.gt.1)then
            do ic=1,ncol
                do i=1,im
                    if(x(ic).lt.xtab(i))go to 10
                enddo
                10 ix(ic)=max0(i-1,1)
                ip1=min(ix(ic)+1,im)
                dx=(xtab(ip1)-xtab(ix(ic)))
                if(abs(dx).gt.1.e-20_r8)then
                    t(ic)=(x(ic)-xtab(ix(ic)))/dx
                else
                    t(ic)=0._r8
                endif
            end do
        else
            ix(:ncol)=1
            t(:ncol)=0._r8
        endif
        if(jm.gt.1)then
            do ic=1,ncol
                do j=1,jm
                    if(y(ic).lt.ytab(j))go to 20
                enddo
                20 jy(ic)=max0(j-1,1)
                jp1=min(jy(ic)+1,jm)
                dy=(ytab(jp1)-ytab(jy(ic)))
                if(abs(dy).gt.1.e-20_r8)then
                    u(ic)=(y(ic)-ytab(jy(ic)))/dy
                else
                    u(ic)=0._r8
                endif
            end do
        else
            jy(:ncol)=1
            u(:ncol)=0._r8
        endif
        30 continue
        do ic=1,ncol
            tu(ic)=t(ic)*u(ic)
            tuc(ic)=t(ic)-tu(ic)
            tcuc(ic)=1._r8-tuc(ic)-u(ic)
            tcu(ic)=u(ic)-tu(ic)
            jp1m(ic)=min(jy(ic)+1,jm)
            ip1m(ic)=min(ix(ic)+1,im)
        enddo
        do ic=1,ncol
            jyc=jy(ic)
            ixc=ix(ic)
            jp1=jp1m(ic)
            ip1=ip1m(ic)
            do k=1,km
                out(ic,k) = tcuc(ic) * table(k,ixc,jyc) + tuc(ic) * table(k,ip1,jyc) + &
                              tu(ic) * table(k,ip1,jp1) + tcu(ic) * table(k,ixc,jp1)
            end do
        end do
        return
    end subroutine binterp

subroutine lookuptable(refr_shell,refi_shell,Rp_BC,&
		cabs,cext,cabs_nobc,cext_nobc,asy,n)
  implicit none
  Integer inpind
  Complex(r8) epscor, epsman
  Complex(r8) m1, s1, m2, s2, x1, x2
! Variables passed to subroutine:
  Integer  i, j, k, kk
  Integer k1,k2,isize
  integer, intent(in) :: n
  real(r8),intent(in) :: refr_shell,refi_shell,Rp_bc(n)
  Real(r8)  c1, c2, a(n)
  real(r8),intent(out) ::  cabs(n),cext(n),cabs_nobc(n),cext_nobc(n),asy(n)
  real(r8) c,qabs
  integer row,time,id
!***********************************************************************
! parameters for optical interpolation  

    real(r8) :: refr_coating(nx),refi_coating(ny),dens_wet_table(nz)
    integer :: jj
    real(r8) :: Rp_bc_table(nz)
    real(r8),dimension(n,3,max(nx,ny,nz)) :: t_pabs,t_pext,t_pabs_bc,t_pext_bc,t_pasm
    real(r8),dimension(n,3,max(nx,ny,nz)) :: t_pabs_nobc,t_pext_nobc
    integer :: z
!***********************************************************************

!  call linspace(minval(refr_shell),maxval(refr_shell),nx,refr_coating)
!  call linspace(minval(refi_shell),maxval(refi_shell),ny,refi_coating)
  call linspace(1.3_r8,1.6_r8,nx,refr_coating)
  call linspace(0._r8,6.0e-3_r8,ny,refi_coating)
  call linspace(0._r8,1200e-9_r8,nz,Rp_BC_table)

!Rp=100.0e-9
!Rp_bc=50.0e-9
!refr_shell=1.5
!refi_shell=0.0012
do k=1,nz
 do j=1,ny 
  do i=1,nx 
    ! interpolate for Rp
    t_pabs(1:n,1,i)=cabs_table(i,j,k,1:n)
    t_pext(1:n,1,i)=cext_table(i,j,k,1:n)
    t_pabs_nobc(1:n,1,i)=cabs_nobc_table(i,j,k,1:n)
    t_pext_nobc(1:n,1,i)=cext_nobc_table(i,j,k,1:n)
    t_pasm(1:n,1,i)=asy_table(i,j,k,1:n)
  enddo
   !!! interpolate for refr_shell
   do jj=1,n
    call interp1(refr_coating,t_pabs(jj,1,1:nx),nx,refr_shell,t_pabs(jj,2,j))
    call interp1(refr_coating,t_pext(jj,1,1:nx),nx,refr_shell,t_pext(jj,2,j))
    call interp1(refr_coating,t_pabs_nobc(jj,1,1:nx),nx,refr_shell,t_pabs_nobc(jj,2,j))
    call interp1(refr_coating,t_pext_nobc(jj,1,1:nx),nx,refr_shell,t_pext_nobc(jj,2,j))
    call interp1(refr_coating,t_pasm(jj,1,1:nx),nx,refr_shell,t_pasm(jj,2,j))
   enddo
 enddo
   !!! interpolate for refi_shell
   do jj=1,n
    call interp1(refi_coating,t_pabs(jj,2,1:ny),ny,refi_shell,t_pabs(jj,3,k))
    call interp1(refi_coating,t_pext(jj,2,1:ny),ny,refi_shell,t_pext(jj,3,k))
    call interp1(refi_coating,t_pabs_nobc(jj,2,1:ny),ny,refi_shell,t_pabs_nobc(jj,3,k))
    call interp1(refi_coating,t_pext_nobc(jj,2,1:ny),ny,refi_shell,t_pext_nobc(jj,3,k))
    call interp1(refi_coating,t_pasm(jj,2,1:ny),ny,refi_shell,t_pasm(jj,3,k))
   enddo
enddo
 !!! interpolate for Rp_BC
do jj=1,n
 call interp1(Rp_BC_table,t_pabs(jj,3,1:nz),nz,Rp_bc(jj),cabs(jj))
 call interp1(Rp_BC_table,t_pext(jj,3,1:nz),nz,Rp_bc(jj),cext(jj))
 call interp1(Rp_BC_table,t_pabs_nobc(jj,3,1:nz),nz,Rp_bc(jj),cabs_nobc(jj))
 call interp1(Rp_BC_table,t_pext_nobc(jj,3,1:nz),nz,Rp_bc(jj),cext_nobc(jj))
 call interp1(Rp_BC_table,t_pasm(jj,3,1:nz),nz,Rp_bc(jj),asy(jj))
enddo

End subroutine lookuptable

subroutine  linspace(start,finish,num,output)
 implicit none
 integer :: k
 integer,intent(in):: num
 real(r8),intent(in) :: start,finish
 real(r8),intent(out) :: output(num)
 real(r8) :: i 
  k = 1
  Do i = start, finish, (finish-start)/(num-1)
    output(k) =i
    k = k + 1
  End Do
  output(num) = finish

end
subroutine  logspace(start,finish,num,output)
 implicit none
 integer :: k
 integer,intent(in):: num
 real(r8),intent(in) :: start,finish
 real(r8),intent(out) :: output(num)
 real(r8) :: i 
  k = 1
  Do i = start, finish, (finish-start)/(num-1)
    output(k) =10**i
    k = k + 1
  End Do
  output(num) = 10**finish

end

subroutine interp1( x, y, n, t, tmp )
        implicit none
        integer :: i, k1, k2
        integer :: n
        real(r8)  :: x(n), y(n), t, tmp, k
 
        !// 如果插值节点与被插值节点相同，则被插值节点处的值等于原节点处的值
        do i = 1, n
                if ( abs( t - x(i) ) < 1d-15 ) then
                        tmp = y(i)
                        return
                end if
        end do
        
        !// 寻找被插值节点属于的区间[X_i, X_i+1]
        i = 1
        do
                if ( t < x(i) ) then
                        exit
                else
                        i = i + 1
                end if
        end do
        
        k2 = i; k1 = k2 - 1
 
        k = ( y(k2) - y(k1) ) / ( x(k2) - x(k1) )
        tmp = y(k1) + k*(t-x(k1))
 
End subroutine interp1

subroutine opt(refr_shell,refi_shell,dens_wet,Radsurf,vol_fra,mode_num,&
               pabs,pext,pasm)
        implicit none
   real(r8),parameter :: pi=3.14159
   real(r8) :: c1, c2, c, Rp(bins), Rm, Rm_bc
   real(r8) ::  n_bc, n_all(bins)
   real(r8) :: N0_bc, xman,rzc,asy1,asy2,asy3
   real(r8) :: out, out1, out2, out3, out4, out5, out6, out7
   real(r8) :: yy(bins), yy1(bins), yy2(bins),yy3(bins),yy4(bins)
   real(r8) :: yy5(bins), yy6(bins), yy7(bins)
   integer :: k1, j, isize
   real(r8) ::  n1(bins), n2(bins), pabs_bc(pcols,pver), pext_bc(pcols,pver)
   real(r8) :: n_bc1, rad, Reff_bc
   real(r8) :: rad_shell,Rp_bc(bins)
   real(r8) :: sigma_m          ! geometric standard deviation of number distribution for accumulation mode
   real(r8) :: vol_fra,refr_shell,refi_shell,dens_wet,Radsurf,mode_num,pabs,pext,pasm
    real(r8) :: ccabs(bins),ccext(bins),cabs_nobc(bins),cext_nobc(bins),asy(bins)
           sigma_m=1.8
           c1=1.0
           c2=1200.0
           k1=1
           do c=c1,c2,(c2-c1)/(bins-1)
             Rp(k1)=c*1e-9
             k1=k1+1
           enddo
           Rp(bins)=c2*1e-9
           Rm=radsurf/exp(2*log(sigma_m)**2)
           Rm_bc=min(Rm,35.e-9_r8)
           ! CT varies with sizes in the accumulation mode
           rad_shell=Rm-Rm_bc
           ! CT is constantly 35nm
!               rad_shell(i,k)=35e-9
           
            Reff_bc=Rm_bc*exp(1.5*(log(1.8))**2)
            N0_BC=vol_fra/(4.0/3.0*pi*Reff_bc**3)


           ! use default parameterization to keep scattering part constant
!               pabs_para(i) = 0.5_r8*cabs(i,1)
!               do nc = 2, ncoef
!                  pabs_para(i) = pabs_para(i) + cheb(nc,i,k)*cabs(i,nc)
!               enddo
           ! ----------------------------------------------------
              yy=0.0
              yy1=0.0
              yy2=0.0
              yy3=0.0
              yy4=0.0
              yy5=0.0
              yy6=0.0
              yy7=0.0
              out1=0.0
              out2=0.0
              out3=0.0
              out4=0.0
              out5=0.0
              out6=0.0
              out7=0.0

              ccabs(1:bins)=0.0
              ccext(1:bins)=0.0
              cabs_nobc(1:bins)=0.0
              cext_nobc(1:bins)=0.0
              asy(1:bins)=0.0
              Rp_BC=Rp-rad_shell
              Rp_BC=max(Rp_BC,0.0)
              call lookuptable(refr_shell,refi_shell,Rp_BC,&
		               ccabs,ccext,cabs_nobc,cext_nobc,asy,bins)

!              ccabs=1.e-12_r8
!              ccext=2.e-12_r8
!              cabs_nobc=0_r8
!              cext_nobc=1.e-12_r8
!              asy=0.65

!              m2_lens=cmplx(resu, 0.0)
              
              n_all=mode_num/(sqrt(2*3.14159)*Rp*log(sigma_m))*exp(-0.5*((log(Rp)-log(Rm))/log(sigma_m))**2)
              n1=N0_BC/(sqrt(2*pi)*(Rp-rad_shell)*log(1.8))*exp(-0.5*((log(Rp-rad_shell)-log(Rm_bc))/log(1.8))**2)
              do j=1,bins
                 if ( Rp(j)<rad_shell )then
                   n1(j)=0._r8
                 endif
              enddo
              n1=min(n_all,n1)
              n2=n_all-n1
              yy=n1*ccabs + n2*cabs_nobc
              yy2=n1*ccext + n2*cext_nobc
              yy1=n1*pi*4.0/3.0*Rp**3*dens_wet + n2*pi*4.0/3.0*Rp**3*dens_wet
              yy7=asy*(yy2-yy)
              out1=sum( 0.5* (yy(2:bins)+yy(1:bins-1)) *(Rp(2:bins)-Rp(1:bins-1)) ) ! bulk absorption
              out2=sum( 0.5* (yy1(2:bins)+yy1(1:bins-1)) *(Rp(2:bins)-Rp(1:bins-1)) ) ! mass of bulk aerosol
              out3=sum( 0.5* (yy2(2:bins)+yy2(1:bins-1)) *(Rp(2:bins)-Rp(1:bins-1)) ) ! bulk extinction
              out7=sum( 0.5* (yy7(2:bins)+yy7(1:bins-1)) *(Rp(2:bins)-Rp(1:bins-1)) ) ! asy * sca
             
!              do j=1,bins
!                 n_all=mode_num(i,k)/(sqrt(2*3.14159)*Rp(j)*log(sigma_m))*exp(-0.5*((log(Rp(j))-log(Rm(i,k)))/log(sigma_m))**2)
!                 n1=N0_BC(i,k)/(sqrt(2*pi)*(Rp(j)-rad_shell(i,k))*log(1.8))*exp(-0.5*((log(Rp(j)-rad_shell(i,k))-log(Rm_bc(i,k)))/log(1.8))**2)
!                 if ( Rp(j)<rad_shell(i,k) )then
!                   n1=0
!                 endif
!                 n1=min(n_all,n1)
!                 n2=n_all-n1
!                 yy(j)=n1*ccabs(j) + n2*cabs_nobc(j)
!                 yy2(j)=n1*ccext(j) + n2*cext_nobc(j)
!                 yy1(j)=n1*pi*4.0/3.0*Rp(j)**3*dens_wet(i,k) + n2*pi*4.0/3.0*Rp(j)**3*dens_wet(i,k)

!                 n_bc1=N0_BC(i,k)/(sqrt(2*pi)*Rp(j)*log(1.8))*exp(-0.5*((log(Rp(j))-log(Rm_bc(i,k)))/log(1.8))**2)
!                 yy3(j)=n_bc1*pi*4.0/3.0*Rp(j)**3*1700
!                 yy4(j)=n1*qabs_lens*pi*Rp(j)**2
!                 yy5(j)=n1*qext_lens*pi*Rp(j)**2
!                 yy7(j)=asy(j)*(yy2(j)-yy(j))

!                 if (j>1) then
!                  out1=out1 + (yy(j)+yy(j-1))*0.5 *(Rp(j)-Rp(j-1)) ! bulk absorption 
!                  out2=out2 + (yy1(j)+yy1(j-1))*0.5 *(Rp(j)-Rp(j-1)) ! mass of bulk aerosol
!                  out3=out3 + (yy2(j)+yy2(j-1))*0.5 *(Rp(j)-Rp(j-1)) ! bulk extinction
!                  out4=out4 + (yy3(j)+yy3(j-1))*0.5 *(Rp(j)-Rp(j-1)) ! mass of BC
!                  out5=out5 + (yy4(j)+yy4(j-1))*0.5 *(Rp(j)-Rp(j-1)) ! BC-containing absorption
!                  out6=out6 + (yy5(j)+yy5(j-1))*0.5 *(Rp(j)-Rp(j-1)) ! BC-containing extinction
!                  out7=out7 + (yy7(j)+yy7(j-1))*0.5 *(Rp(j)-Rp(j-1)) ! asy * sca
!                 endif
!               enddo ! j=1,bins

                  pabs=out1/out2
                ! pext calculated by core-shell Mie theory
                  pext=out3/out2

                  pasm=out7/(out3-out1)
!                   pabs_bc(i,k)=out5/out4
!                  pext_bc(i,k)=out6/out4

            ! parameterization by steve Ghan
!                  pext_para(i) = 0.5_r8*cext(i,1)
!                  do nc = 2, ncoef
!                     pext_para(i) = pext_para(i) + cheb(nc,i,k)*cext(i,nc)
!                  enddo
!                  pext_para(i) = exp(pext_para(i))
!                  pext(i)=pabs(i)+(pext_para(i)-pabs_para(i))
!                  pext(i)=pext_para(i)

      end subroutine opt

end module modal_aer_opt
