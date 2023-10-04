MODULE mo_prp

  ! 
  ! Module to perform online Partial Radiation Perturbations
  ! (PRP) in ECHAM. 
  !
  ! Partial radiation perturbations (dR) and heating rates
  ! (dQ) are calculated by averaging the forward and back-
  ! ward calls:
  ! 
  ! dR_x = ((a-c) - (b-d))/2
  !
  ! where a is the flux with the current state, b with the 
  ! stored state, c is with the current state but state 
  ! variable x from the stored state, and d is the stored 
  ! state with state variable x from the current state.
  !
  ! The module employs its own time management, such that
  ! calls can be made every prp_dt hours, e.g. 2, 10 or 22, 
  ! which will all yield eventually covering the full diurnal
  ! cycle. The intervals have to be at radiation calls, 
  ! because the subroutine prp is called from within 
  ! radiation.
  !
  ! It is recommended to set prp_dt to 2 hours, while 
  ! writing the reference state (prprecord = .true.), even 
  ! if less frequent PRP calls are to be used.
  !
  ! Online averaging is done by default to monthly means.
  ! The interval can be set at runtime by prp_avg.
  ! 


USE mo_linked_list,    ONLY: t_stream 
USE mo_kind,           ONLY: wp
USE mo_time_event,     ONLY: io_time_event, time_event

IMPLICIT NONE

PUBLIC :: init_prp, construct_prp_stream, destruct_prp_stream, prp, prp_read_state

! From the prpctl namelist:

CHARACTER(len=256), PUBLIC, SAVE    :: prppath   = '/default_prppath/'  
CHARACTER(len=64),  PUBLIC, SAVE    :: prp_ref_expname = 'expid'
INTEGER,            PUBLIC, SAVE    :: prp_ref_year = 1979
LOGICAL,            PUBLIC, SAVE    :: prprecord = .true.
LOGICAL,            PUBLIC, SAVE    :: prp3doutput = .false. ! if .true. then 3D output is generated
LOGICAL,            PUBLIC, SAVE    :: prp_mpiesm = .true. ! if .true. then echam is executed within mpiesm: folders then change
INTEGER,            PUBLIC, SAVE    :: prp_dt = 2
TYPE(io_time_event),PUBLIC, SAVE    :: prp_write
TYPE(io_time_event),PUBLIC, SAVE    :: prp_read
TYPE(io_time_event),PUBLIC, SAVE    :: prp_avg   = io_time_event(1,'months','last',0)
TYPE(time_event),   PUBLIC, SAVE    :: ev_trigprp
LOGICAL,            PUBLIC, SAVE    :: l_trigprp  = .true.
LOGICAL,            PUBLIC, SAVE    :: prp_sp  = .false. ! set on .true. to diagnose total/direct/indirect effects of aerosols and 'direct' effect of clouds (without indirect effect)
LOGICAL,            PUBLIC, SAVE    :: prp_single_sp  = .false. ! diagnose direct/indirect effects of each individual plumes
LOGICAL,            PUBLIC, SAVE    :: prp_feedbacks  = .true. ! added to only diagnose simple plumes direct/indorect effects in historical simulations (avoid to many calls when not interested in climate feedbacks)
INTEGER,            PUBLIC, SAVE    :: iaero_ref = 3      ! must be set to the same case as the reference simulation (default is 3 for piControl)
LOGICAL,            PUBLIC, SAVE    :: prprads = .false. ! if .true. then surface generated
LOGICAL,            PUBLIC, SAVE    :: prpraf = .false. ! if .true. then clear-sky output is generated
LOGICAL,            PUBLIC, SAVE    :: prp3drad = .false. ! if .true. then 3D radiation output is generated


! Module-wide variables:

INTEGER,            PUBLIC ,SAVE    :: prp_j = 0
INTEGER,            PUBLIC ,SAVE    :: read_istep = -1
REAL(wp),           PUBLIC, POINTER :: filetime(:)
REAL(wp),           PUBLIC, POINTER :: filetod(:)
REAL(wp),           PUBLIC, POINTER :: filedom(:)
REAL(wp),           PUBLIC, POINTER :: timediff(:)
INTEGER,            PUBLIC, SAVE    :: oldmonth = -1
INTEGER, PARAMETER :: nplumes   = 9  !< Number of plumes

! Output PRP stream:

  TYPE (t_stream),    PUBLIC, POINTER :: prp_stream
  
  ! Time variables for output stream:

  REAL(wp),           PUBLIC, POINTER :: ptr_tod(:,:)            ! Time of day, between zero and one
  REAL(wp),           PUBLIC, POINTER :: ptr_dom(:,:)            ! Day of month

  ! State variables:
  
 ! REAL(wp),           PUBLIC, POINTER :: ptr_iaero(:)
 ! REAL(wp),           PUBLIC, POINTER :: ptr_cemiss(:)
  REAL(wp),           PUBLIC, POINTER :: ptr_loglac(:,:)
  REAL(wp),           PUBLIC, POINTER :: ptr_loland(:,:) 
  REAL(wp),           PUBLIC, POINTER :: ptr_ktype(:,:)          !< type of convection
  REAL(wp),           PUBLIC, POINTER :: ptr_pgeom1(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: ptr_pp_fl(:,:,:) 
  REAL(wp),           PUBLIC, POINTER :: ptr_pp_hl(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: ptr_pp_sfc(:,:)  
  REAL(wp),           PUBLIC, POINTER :: ptr_cos_mu0(:,:)        ! solar zenith angle

  REAL(wp),           PUBLIC, POINTER :: ptr_q_vap(:,:,:)        ! water vapor field
  REAL(wp),           PUBLIC, POINTER :: ptr_tk_fl(:,:,:)        ! atmospheric temperature at full levels
  REAL(wp),           PUBLIC, POINTER :: ptr_tk_hl(:,:,:)        ! atmospheric temperature at half levels
  REAL(wp),           PUBLIC, POINTER :: ptr_tk_sfc(:,:)         ! surface temperature field
  REAL(wp),           PUBLIC, POINTER :: ptr_q_liq(:,:,:)        ! cloud liquid content
  REAL(wp),           PUBLIC, POINTER :: ptr_q_ice(:,:,:)        ! cloud ice content
  REAL(wp),           PUBLIC, POINTER :: ptr_cdnc(:,:,:)         ! cloud droplet number concentration
  REAL(wp),           PUBLIC, POINTER :: ptr_cld_frc(:,:,:)      ! cloud fraction
  REAL(wp),           PUBLIC, POINTER :: ptr_cld_cvr(:,:)        ! projected cloud cover
  REAL(wp),           PUBLIC, POINTER :: ptr_alb_vis_dir(:,:)    ! box surface albedo visible range direct
  REAL(wp),           PUBLIC, POINTER :: ptr_alb_vis_dif(:,:)    ! box surface albedo visible range diffuse
  REAL(wp),           PUBLIC, POINTER :: ptr_alb_nir_dir(:,:)    ! box surface albedo NIR range direct
  REAL(wp),           PUBLIC, POINTER :: ptr_alb_nir_dif(:,:)    ! box surface albedo NIR range diffuse

  REAL(wp),           PUBLIC, POINTER :: tropo_p(:,:)            ! tropopause pressure

  ! Composition:

  REAL(wp),           PUBLIC, POINTER :: ptr_m_co2(:,:,:)        ! CO2 as 3D field
  REAL(wp),           PUBLIC, POINTER :: ptr_m_o3(:,:,:)         
  REAL(wp),           PUBLIC, POINTER :: ptr_m_ch4(:,:,:) 
  REAL(wp),           PUBLIC, POINTER :: ptr_m_n2o(:,:,:) 
  REAL(wp),           PUBLIC, POINTER :: ptr_m_cfc1(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: ptr_m_cfc2(:,:,:) 
  REAL(wp),           PUBLIC, POINTER :: ptr_m_o2(:,:,:) 

  ! PRP diagnostics:

  REAL(wp),           PUBLIC, POINTER :: dR_trad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_srad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_trads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_srads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_traf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_sraf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_trafs(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_srafs(:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_tmp_trad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_tmp_srad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_tmp_trads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_tmp_srads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_tmp_traf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_tmp_sraf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_tmp_trafs(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_tmp_srafs(:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_plk_trad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_plk_srad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_plk_trads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_plk_srads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_plk_traf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_plk_sraf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_plk_trafs(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_plk_srafs(:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_lr_trad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_lr_srad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_lr_trads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_lr_srads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_lr_traf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_lr_sraf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_lr_trafs(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_lr_srafs(:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_str_trad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_str_srad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_str_trads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_str_srads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_str_traf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_str_sraf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_str_trafs(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_str_srafs(:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_vap_trad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_vap_srad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_vap_trads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_vap_srads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_vap_traf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_vap_sraf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_vap_trafs(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_vap_srafs(:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_cld_trad0(:,:)   ! total effect of clouds. 
  REAL(wp),           PUBLIC, POINTER :: dR_cld_srad0(:,:)   ! If anthropogenic aerosols (simple-plumes) are present,
  REAL(wp),           PUBLIC, POINTER :: dR_cld_trads(:,:)   ! this includes the indirect effect of aerosols on clouds
  REAL(wp),           PUBLIC, POINTER :: dR_cld_srads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_cld_traf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_cld_sraf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_cld_trafs(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_cld_srafs(:,:)
  ! -------------------------------------------------------  ! Pointers added to the stream only if l_spprp is .true.
  REAL(wp),           PUBLIC, POINTER :: dR_cdd_trad0(:,:)   ! direct effect of clouds only, without indirect effect of simple-plumes (Twomey effect)
  REAL(wp),           PUBLIC, POINTER :: dR_cdd_srad0(:,:)   ! If the simple-plumes are not present, this should be equal to the total effect of clouds above.
  REAL(wp),           PUBLIC, POINTER :: dR_cdd_trads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_cdd_srads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_cdd_traf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_cdd_sraf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_cdd_trafs(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_cdd_srafs(:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_aer_trad0(:,:)   ! all aerosols and all effects (volcanoes, simple-plumes direct/indirect ...), depending on the scenario
  REAL(wp),           PUBLIC, POINTER :: dR_aer_srad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_aer_trads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_aer_srads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_aer_traf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_aer_sraf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_aer_trafs(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_aer_srafs(:,:) 

  REAL(wp),           PUBLIC, POINTER :: dR_sp_trad0(:,:)    ! direct + indirect effects of simple-plumes
  REAL(wp),           PUBLIC, POINTER :: dR_sp_srad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_sp_trads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_sp_srads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_sp_traf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_sp_sraf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_sp_trafs(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_sp_srafs(:,:) 

  REAL(wp),           PUBLIC, POINTER :: dR_spd_trad0(:,:)   ! direct effect of simple-plumes
  REAL(wp),           PUBLIC, POINTER :: dR_spd_srad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_spd_trads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_spd_srads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_spd_traf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_spd_sraf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_spd_trafs(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_spd_srafs(:,:) 

  REAL(wp),           PUBLIC, POINTER :: dR_spi_trad0(:,:)   ! indirect effect of simple-plumes
  REAL(wp),           PUBLIC, POINTER :: dR_spi_srad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_spi_trads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_spi_srads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_spi_traf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_spi_sraf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_spi_trafs(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_spi_srafs(:,:) 
! -------------------------------------------------------------------------------------------
  REAL(wp),           PUBLIC, POINTER :: dR_alb_trad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_alb_srad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_alb_trads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_alb_srads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_alb_traf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_alb_sraf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_alb_trafs(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_alb_srafs(:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_co2_trad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_co2_srad0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_co2_trads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_co2_srads(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_co2_traf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_co2_sraf0(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_co2_trafs(:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_co2_srafs(:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_tmp_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_tmp_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_tmp_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_tmp_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_plk_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_plk_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_plk_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_plk_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_lr_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_lr_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_lr_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_lr_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_str_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_str_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_str_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_str_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_vap_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_vap_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_vap_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_vap_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_cld_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_cld_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_cld_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_cld_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_alb_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_alb_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_alb_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_alb_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_co2_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_co2_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_co2_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_co2_sraf(:,:,:)
! ---------------------------------------------------------------------------------------
  REAL(wp),           PUBLIC, POINTER :: dR_cdd_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_cdd_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_cdd_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_cdd_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_aer_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_aer_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_aer_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_aer_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_sp_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_sp_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_sp_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_sp_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_spd_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_spd_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_spd_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_spd_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dR_spi_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_spi_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_spi_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dR_spi_sraf(:,:,:)
!-----------------------------------------------------------------------------------------
  REAL(wp),           PUBLIC, POINTER :: dQ_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dQ_tmp_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_tmp_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_tmp_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_tmp_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dQ_plk_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_plk_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_plk_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_plk_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dQ_lr_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_lr_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_lr_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_lr_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dQ_str_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_str_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_str_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_str_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dQ_vap_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_vap_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_vap_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_vap_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dQ_cld_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_cld_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_cld_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_cld_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dQ_alb_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_alb_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_alb_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_alb_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dQ_co2_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_co2_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_co2_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_co2_sraf(:,:,:)
! ----------------------------------------------------------------------------------------
  REAL(wp),           PUBLIC, POINTER :: dQ_cdd_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_cdd_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_cdd_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_cdd_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dQ_aer_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_aer_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_aer_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_aer_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dQ_sp_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_sp_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_sp_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_sp_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dQ_spd_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_spd_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_spd_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_spd_sraf(:,:,:)

  REAL(wp),           PUBLIC, POINTER :: dQ_spi_trad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_spi_srad(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_spi_traf(:,:,:)
  REAL(wp),           PUBLIC, POINTER :: dQ_spi_sraf(:,:,:)

  type real_ptr
    REAL(wp), PUBLIC, POINTER :: p(:,:)
  end type real_ptr 

  TYPE(real_ptr), PUBLIC, dimension(nplumes) :: dR_single_sp_srad0     ! TOA single simple-plumes 
  TYPE(real_ptr), PUBLIC, dimension(nplumes) :: dR_single_sp_trad0
  TYPE(real_ptr), PUBLIC, dimension(nplumes) :: dR_single_spd_srad0
  TYPE(real_ptr), PUBLIC, dimension(nplumes) :: dR_single_spd_trad0
  TYPE(real_ptr), PUBLIC, dimension(nplumes) :: dR_single_spi_srad0
  TYPE(real_ptr), PUBLIC, dimension(nplumes) :: dR_single_spi_trad0

  TYPE(real_ptr), PUBLIC, dimension(nplumes) :: dR_single_sp_sraf0     ! clear-sky TOA single simple-plumes 
  TYPE(real_ptr), PUBLIC, dimension(nplumes) :: dR_single_sp_traf0
  TYPE(real_ptr), PUBLIC, dimension(nplumes) :: dR_single_spd_sraf0
  TYPE(real_ptr), PUBLIC, dimension(nplumes) :: dR_single_spd_traf0
  TYPE(real_ptr), PUBLIC, dimension(nplumes) :: dR_single_spi_sraf0
  TYPE(real_ptr), PUBLIC, dimension(nplumes) :: dR_single_spi_traf0

  type real_ptr_3d
  REAL(wp), PUBLIC, POINTER :: p(:,:,:)
  end type real_ptr_3d

  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dQ_single_sp_srad    ! heating-rates single simple-plumes
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dQ_single_sp_trad
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dQ_single_spd_srad
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dQ_single_spd_trad
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dQ_single_spi_srad
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dQ_single_spi_trad

  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dQ_single_sp_sraf    ! clear-sky heating-rates single simple-plumes
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dQ_single_sp_traf
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dQ_single_spd_sraf
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dQ_single_spd_traf
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dQ_single_spi_sraf
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dQ_single_spi_traf

  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dR_single_sp_srad    ! 3D single simple-plumes
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dR_single_sp_trad
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dR_single_spd_srad
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dR_single_spd_trad
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dR_single_spi_srad
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dR_single_spi_trad

  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dR_single_sp_sraf    ! 3D single simple-plumes clear-sky
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dR_single_sp_traf
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dR_single_spd_sraf
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dR_single_spd_traf
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dR_single_spi_sraf
  TYPE(real_ptr_3d), PUBLIC, dimension(nplumes) :: dR_single_spi_traf

  TYPE(real_ptr), PUBLIC, dimension(nplumes+1) :: ptr_aod_sp  ! AOD field of individual plumes (1,2,3...) and all plumes (10)
  TYPE(real_ptr), PUBLIC, dimension(nplumes+1) :: ptr_aod_bg  ! AOD field of individual plumes (1,2,3...) and all plumes (10)
! ----------------------------------------------------------------------------------------
  CONTAINS

!-----------------------------------------------
!-----------------------------------------------

  SUBROUTINE init_prp
  !
  ! This subroutine reads the PRP namelist.
  ! 
  ! It is called from initialize.f90
  !

    USE mo_namelist,       ONLY: open_nml, position_nml, close_nml, &
                                 POSITIONED, MISSING
    USE mo_mpi,            ONLY: p_parallel_io, p_io, p_bcast
    USE mo_exception,      ONLY: finish, message
    USE mo_submodel,       ONLY: print_value
    USE mo_filename,       ONLY: out_filetype, trac_filetype
    USE	mo_time_control,   ONLY: p_bcast_event, echam_ev_init


    IMPLICIT NONE

    integer:: inml, iunit, ierr

    namelist /prpctl/    prppath, prp_ref_expname, prp_ref_year, &
                         prprecord, prp_dt, prp_avg, prp3doutput, prp_mpiesm, &
                         prp_sp, iaero_ref, prp_single_sp, prprads, prpraf, prp3drad, prp_feedbacks

    call message('','Initializing mo_prp')

    if (p_parallel_io) then
       
     inml = open_nml ('namelist.echam')
     iunit = position_nml ('PRPCTL', inml, status=ierr)
     SELECT CASE (ierr)
     CASE (POSITIONED)
       trac_filetype = 0
       READ (iunit, prpctl)
       IF(trac_filetype == 0) trac_filetype = out_filetype
     CASE(MISSING)
       CALL message('init_prp', 'Cannot find namelist PRPCTL in namelist.echam!')
     END SELECT
     CALL close_nml(inml)

    end if

    call p_bcast(prppath, p_io)
    call p_bcast(prp_ref_expname, p_io)
    call p_bcast(prp_ref_year, p_io)
    call p_bcast(prprecord, p_io)
    call p_bcast(prp3doutput, p_io)
    call p_bcast(prp_mpiesm, p_io)
    call p_bcast(prp_sp, p_io)
    call p_bcast(iaero_ref, p_io)
    call p_bcast(prp_single_sp, p_io)
    call p_bcast(prp_feedbacks, p_io)
    call p_bcast(prprads, p_io)
    call p_bcast(prp3drad, p_io)
    call p_bcast(prpraf, p_io)
    call p_bcast(prp_dt, p_io)
    call p_bcast_event(prp_avg, p_io)

    call message('','----------------------------')
    call message('','prppath: '// trim(prppath))
    call message('','prp_ref_expname: '// trim(prp_ref_expname))
    call print_value('prp_ref_year', prp_ref_year)
    call print_value('prprecord', prprecord)
    call print_value('prp3doutput', prp3doutput)
    call print_value('prp_mpiesm: echam is run within mpiesm', prp_mpiesm)
    call print_value('prp_dt', prp_dt)
    call print_value('iaero_ref', iaero_ref)
    !call print_value('prp anthropogenic aerosols (Simple-plumes)', prp_sp)
    !call print_value('prp individual anthropogenic aerosol sources', prp_single_sp)
    !call print_value('prp climate feedbacks', prp_feedbacks)
    call print_value('prprads', prprads)
    call print_value('prp3drad', prp3drad)
    call print_value('prpraf', prpraf)
    call print_value('prp_dt', prp_dt)
    call message('','----------------------------')

    prp_write = io_time_event(prp_dt,'hours','first',1)
    prp_read  = io_time_event(prp_dt,'hours','first',0)
    CALL echam_ev_init(ev_trigprp,prp_read,'PRP events','present')

  END SUBROUTINE init_prp

!-----------------------------------------------
!-----------------------------------------------
!-----------------------------------------------
!-----------------------------------------------

  SUBROUTINE construct_prp_stream

    !
    ! This subroutine creates an I/O stream for PRP. Depending
    ! on whether prprecord is true or false, the stream will be
    ! defined differently.
    !

    USE mo_exception,      ONLY: message
    USE mo_memory_base,    ONLY: new_stream, add_stream_element
    USE mo_linked_list,    ONLY: NETCDF
    USE mo_time_event,     ONLY: io_time_event
    USE mo_control,        ONLY: nlev, nlevp1

    IMPLICIT NONE
!-------------------------------------------------------------------------------
    INTEGER :: iplume !< loop over all plumes
    character(len=12), dimension(nplumes) :: str_dR_single_sp_srad0
    character(len=12), dimension(nplumes) :: str_dR_single_sp_trad0
    character(len=13), dimension(nplumes) :: str_dR_single_spd_srad0
    character(len=13), dimension(nplumes) :: str_dR_single_spd_trad0
    character(len=13), dimension(nplumes) :: str_dR_single_spi_srad0
    character(len=13), dimension(nplumes) :: str_dR_single_spi_trad0

    character(len=12), dimension(nplumes) :: str_dR_single_sp_sraf0
    character(len=12), dimension(nplumes) :: str_dR_single_sp_traf0
    character(len=13), dimension(nplumes) :: str_dR_single_spd_sraf0
    character(len=13), dimension(nplumes) :: str_dR_single_spd_traf0
    character(len=13), dimension(nplumes) :: str_dR_single_spi_sraf0
    character(len=13), dimension(nplumes) :: str_dR_single_spi_traf0

    character(len=11), dimension(nplumes) :: str_dR_single_sp_trad
    character(len=11), dimension(nplumes) :: str_dR_single_sp_srad
    character(len=12), dimension(nplumes) :: str_dR_single_spd_trad
    character(len=12), dimension(nplumes) :: str_dR_single_spd_srad
    character(len=12), dimension(nplumes) :: str_dR_single_spi_trad
    character(len=12), dimension(nplumes) :: str_dR_single_spi_srad

    character(len=11), dimension(nplumes) :: str_dR_single_sp_traf
    character(len=11), dimension(nplumes) :: str_dR_single_sp_sraf
    character(len=12), dimension(nplumes) :: str_dR_single_spd_traf
    character(len=12), dimension(nplumes) :: str_dR_single_spd_sraf
    character(len=12), dimension(nplumes) :: str_dR_single_spi_traf
    character(len=12), dimension(nplumes) :: str_dR_single_spi_sraf

    character(len=11), dimension(nplumes) :: str_dQ_single_sp_trad
    character(len=11), dimension(nplumes) :: str_dQ_single_sp_srad
    character(len=12), dimension(nplumes) :: str_dQ_single_spd_trad
    character(len=12), dimension(nplumes) :: str_dQ_single_spd_srad
    character(len=12), dimension(nplumes) :: str_dQ_single_spi_trad
    character(len=12), dimension(nplumes) :: str_dQ_single_spi_srad

    character(len=11), dimension(nplumes) :: str_dQ_single_sp_traf
    character(len=11), dimension(nplumes) :: str_dQ_single_sp_sraf
    character(len=12), dimension(nplumes) :: str_dQ_single_spd_traf
    character(len=12), dimension(nplumes) :: str_dQ_single_spd_sraf
    character(len=12), dimension(nplumes) :: str_dQ_single_spi_traf
    character(len=12), dimension(nplumes) :: str_dQ_single_spi_sraf

    character(len=7), dimension(nplumes) :: str_aod_single_sp
    character(len=7), dimension(nplumes) :: str_aod_single_bg
!--------------------------------------------------------------------------------
    call message('','Create: construct_prp_stream')
    
    IF (prprecord) THEN
       CALL new_stream(prp_stream, 'prp',filetype=NETCDF, &
                       lpost=.TRUE.,lpout=.TRUE.,lrerun=.FALSE.,linit=.FALSE., &
                       interval=prp_write)
    ELSE
       CALL new_stream(prp_stream, 'prp',filetype=NETCDF, &
                       lpost=.TRUE.,lpout=.TRUE.,lrerun=.FALSE.,linit=.FALSE., &
                       interval=prp_avg)
    END IF

    CALL add_stream_element(prp_stream, 'q_vap', ptr_q_vap, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'tk_fl', ptr_tk_fl, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'tk_hl', ptr_tk_hl, lpost=prprecord, lrerun=.FALSE., klev=nlevp1)
    CALL add_stream_element(prp_stream, 'tk_sfc', ptr_tk_sfc, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'q_liq', ptr_q_liq, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'q_ice', ptr_q_ice, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'cdnc', ptr_cdnc, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'cld_frc', ptr_cld_frc, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'cld_cvr', ptr_cld_cvr, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'alb_vis_dir', ptr_alb_vis_dir, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'alb_vis_dif', ptr_alb_vis_dif, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'alb_nir_dir', ptr_alb_nir_dir, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'alb_nir_dif', ptr_alb_nir_dif, lpost=prprecord, lrerun=.FALSE.)

 !   CALL add_stream_element(prp_stream, 'iaero', ptr_iaero, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'loland', ptr_loland, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'loglac', ptr_loglac, lpost=prprecord, lrerun=.FALSE.)
  !  CALL add_stream_element(prp_stream, 'cemiss', ptr_cemiss, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'cos_mu0', ptr_cos_mu0, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'm_co2', ptr_m_co2, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'm_o3', ptr_m_o3, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'm_ch4', ptr_m_ch4, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'm_n2o', ptr_m_n2o, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'm_cfc1', ptr_m_cfc1, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'm_cfc2', ptr_m_cfc2, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'm_o2', ptr_m_o2, lpost=prprecord, lrerun=.FALSE.)    
    CALL add_stream_element(prp_stream, 'pgeom1', ptr_pgeom1, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'pp_fl', ptr_pp_fl, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'pp_hl', ptr_pp_hl, lpost=prprecord, lrerun=.FALSE., klev=nlevp1)
    CALL add_stream_element(prp_stream, 'pp_sfc', ptr_pp_sfc, lpost=prprecord, lrerun=.FALSE.)
    CALL add_stream_element(prp_stream, 'ktype', ptr_ktype, lpost=prprecord, lrerun=.FALSE.)

    IF (prprecord) THEN
       
       CALL add_stream_element(prp_stream, 'tod', ptr_tod, lpost=.TRUE., lrerun=.FALSE.)
       CALL add_stream_element(prp_stream, 'dom', ptr_dom, lpost=.TRUE., lrerun=.FALSE.)

    END IF

    IF (.not.prprecord) THEN

       CALL add_stream_element(prp_stream, 'tropo_p', tropo_p, lpost=.FALSE., lrerun=.FALSE., laccu=.FALSE.)

       CALL add_stream_element(prp_stream, 'dR_trad0', dR_trad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_srad0', dR_srad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_trads', dR_trads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_srads', dR_srads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_traf0', dR_traf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_sraf0', dR_sraf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_trafs', dR_trafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_srafs', dR_srafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_trad', dR_trad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_srad', dR_srad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_traf', dR_traf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_sraf', dR_sraf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)

       IF (prp_feedbacks) THEN
       CALL add_stream_element(prp_stream, 'dR_tmp_trad0', dR_tmp_trad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_tmp_srad0', dR_tmp_srad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_tmp_trads', dR_tmp_trads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_tmp_srads', dR_tmp_srads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_tmp_traf0', dR_tmp_traf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_tmp_sraf0', dR_tmp_sraf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_tmp_trafs', dR_tmp_trafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_tmp_srafs', dR_tmp_srafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

       CALL add_stream_element(prp_stream, 'dR_plk_trad0', dR_plk_trad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_plk_srad0', dR_plk_srad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_plk_trads', dR_plk_trads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_plk_srads', dR_plk_srads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_plk_traf0', dR_plk_traf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_plk_sraf0', dR_plk_sraf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_plk_trafs', dR_plk_trafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_plk_srafs', dR_plk_srafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

       CALL add_stream_element(prp_stream, 'dR_lr_trad0', dR_lr_trad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_lr_srad0', dR_lr_srad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_lr_trads', dR_lr_trads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_lr_srads', dR_lr_srads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_lr_traf0', dR_lr_traf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_lr_sraf0', dR_lr_sraf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_lr_trafs', dR_lr_trafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_lr_srafs', dR_lr_srafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

       CALL add_stream_element(prp_stream, 'dR_str_trad0', dR_str_trad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_str_srad0', dR_str_srad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_str_trads', dR_str_trads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_str_srads', dR_str_srads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_str_traf0', dR_str_traf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_str_sraf0', dR_str_sraf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_str_trafs', dR_str_trafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_str_srafs', dR_str_srafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

       CALL add_stream_element(prp_stream, 'dR_vap_trad0', dR_vap_trad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_vap_srad0', dR_vap_srad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_vap_trads', dR_vap_trads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_vap_srads', dR_vap_srads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_vap_traf0', dR_vap_traf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_vap_sraf0', dR_vap_sraf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_vap_trafs', dR_vap_trafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_vap_srafs', dR_vap_srafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

       CALL add_stream_element(prp_stream, 'dR_cld_trad0', dR_cld_trad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_cld_srad0', dR_cld_srad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_cld_trads', dR_cld_trads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_cld_srads', dR_cld_srads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_cld_traf0', dR_cld_traf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_cld_sraf0', dR_cld_sraf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_cld_trafs', dR_cld_trafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_cld_srafs', dR_cld_srafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

       CALL add_stream_element(prp_stream, 'dR_aer_trad0', dR_aer_trad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_aer_srad0', dR_aer_srad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_aer_trads', dR_aer_trads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_aer_srads', dR_aer_srads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_aer_traf0', dR_aer_traf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_aer_sraf0', dR_aer_sraf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_aer_trafs', dR_aer_trafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_aer_srafs', dR_aer_srafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       
       CALL add_stream_element(prp_stream, 'dR_alb_trad0', dR_alb_trad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_alb_srad0', dR_alb_srad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_alb_trads', dR_alb_trads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_alb_srads', dR_alb_srads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_alb_traf0', dR_alb_traf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_alb_sraf0', dR_alb_sraf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_alb_trafs', dR_alb_trafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_alb_srafs', dR_alb_srafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

       CALL add_stream_element(prp_stream, 'dR_co2_trad0', dR_co2_trad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_co2_srad0', dR_co2_srad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_co2_trads', dR_co2_trads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_co2_srads', dR_co2_srads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_co2_traf0', dR_co2_traf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_co2_sraf0', dR_co2_sraf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_co2_trafs', dR_co2_trafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
       CALL add_stream_element(prp_stream, 'dR_co2_srafs', dR_co2_srafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
             
       CALL add_stream_element(prp_stream, 'dR_tmp_trad', dR_tmp_trad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_tmp_srad', dR_tmp_srad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_tmp_traf', dR_tmp_traf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_tmp_sraf', dR_tmp_sraf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)

       CALL add_stream_element(prp_stream, 'dR_plk_trad', dR_plk_trad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_plk_srad', dR_plk_srad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_plk_traf', dR_plk_traf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_plk_sraf', dR_plk_sraf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)

       CALL add_stream_element(prp_stream, 'dR_lr_trad', dR_lr_trad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_lr_srad', dR_lr_srad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_lr_traf', dR_lr_traf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_lr_sraf', dR_lr_sraf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)

       CALL add_stream_element(prp_stream, 'dR_str_trad', dR_str_trad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_str_srad', dR_str_srad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_str_traf', dR_str_traf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_str_sraf', dR_str_sraf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)

       CALL add_stream_element(prp_stream, 'dR_vap_trad', dR_vap_trad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_vap_srad', dR_vap_srad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_vap_traf', dR_vap_traf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_vap_sraf', dR_vap_sraf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)

       CALL add_stream_element(prp_stream, 'dR_cld_trad', dR_cld_trad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_cld_srad', dR_cld_srad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_cld_traf', dR_cld_traf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_cld_sraf', dR_cld_sraf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)

       CALL add_stream_element(prp_stream, 'dR_aer_trad', dR_aer_trad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_aer_srad', dR_aer_srad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_aer_traf', dR_aer_traf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_aer_sraf', dR_aer_sraf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)

       CALL add_stream_element(prp_stream, 'dR_alb_trad', dR_alb_trad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_alb_srad', dR_alb_srad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_alb_traf', dR_alb_traf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_alb_sraf', dR_alb_sraf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)

       CALL add_stream_element(prp_stream, 'dR_co2_trad', dR_co2_trad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_co2_srad', dR_co2_srad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_co2_traf', dR_co2_traf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       CALL add_stream_element(prp_stream, 'dR_co2_sraf', dR_co2_sraf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       END IF

      IF (prp_sp) THEN
        CALL add_stream_element(prp_stream, 'dR_cdd_trad0', dR_cdd_trad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_cdd_srad0', dR_cdd_srad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
 
        CALL add_stream_element(prp_stream, 'dR_sp_trad0', dR_sp_trad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_sp_srad0', dR_sp_srad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
 
        CALL add_stream_element(prp_stream, 'dR_spd_trad0', dR_spd_trad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_spd_srad0', dR_spd_srad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
 
        CALL add_stream_element(prp_stream, 'dR_spi_trad0', dR_spi_trad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_spi_srad0', dR_spi_srad0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

        CALL add_stream_element(prp_stream, 'aod_sp', ptr_aod_sp(nplumes+1)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'aod_bg', ptr_aod_bg(nplumes+1)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

        IF (prp_single_sp) THEN
         str_dR_single_sp_trad0 = (/'dR_sp1_trad0', 'dR_sp2_trad0', 'dR_sp3_trad0', 'dR_sp4_trad0', 'dR_sp5_trad0', 'dR_sp6_trad0', 'dR_sp7_trad0', 'dR_sp8_trad0', 'dR_sp9_trad0'/)
         str_dR_single_sp_srad0 = (/'dR_sp1_srad0', 'dR_sp2_srad0', 'dR_sp3_srad0', 'dR_sp4_srad0', 'dR_sp5_srad0', 'dR_sp6_srad0', 'dR_sp7_srad0', 'dR_sp8_srad0', 'dR_sp9_srad0'/)
         str_dR_single_spd_trad0 = (/'dR_spd1_trad0', 'dR_spd2_trad0', 'dR_spd3_trad0', 'dR_spd4_trad0', 'dR_spd5_trad0', 'dR_spd6_trad0', 'dR_spd7_trad0', 'dR_spd8_trad0', 'dR_spd9_trad0'/)
         str_dR_single_spd_srad0 = (/'dR_spd1_srad0', 'dR_spd2_srad0', 'dR_spd3_srad0', 'dR_spd4_srad0', 'dR_spd5_srad0', 'dR_spd6_srad0', 'dR_spd7_srad0', 'dR_spd8_srad0', 'dR_spd9_srad0'/)
         str_dR_single_spi_trad0 = (/'dR_spi1_trad0', 'dR_spi2_trad0', 'dR_spi3_trad0', 'dR_spi4_trad0', 'dR_spi5_trad0', 'dR_spi6_trad0', 'dR_spi7_trad0', 'dR_spi8_trad0', 'dR_spi9_trad0'/)
         str_dR_single_spi_srad0 = (/'dR_spi1_srad0', 'dR_spi2_srad0', 'dR_spi3_srad0', 'dR_spi4_srad0', 'dR_spi5_srad0', 'dR_spi6_srad0', 'dR_spi7_srad0', 'dR_spi8_srad0', 'dR_spi9_srad0'/)
         
         str_aod_single_sp = (/'aod_sp1', 'aod_sp2', 'aod_sp3', 'aod_sp4', 'aod_sp5', 'aod_sp6', 'aod_sp7', 'aod_sp8', 'aod_sp9'/)
         str_aod_single_bg = (/'aod_bg1', 'aod_bg2', 'aod_bg3', 'aod_bg4', 'aod_bg5', 'aod_bg6', 'aod_bg7', 'aod_bg8', 'aod_bg9'/)
        
         DO iplume=1,nplumes
           CALL add_stream_element(prp_stream, str_dR_single_sp_trad0(iplume), dR_single_sp_trad0(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
           CALL add_stream_element(prp_stream, str_dR_single_sp_srad0(iplume), dR_single_sp_srad0(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
           CALL add_stream_element(prp_stream, str_dR_single_spd_trad0(iplume), dR_single_spd_trad0(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
           CALL add_stream_element(prp_stream, str_dR_single_spd_srad0(iplume), dR_single_spd_srad0(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
           CALL add_stream_element(prp_stream, str_dR_single_spi_trad0(iplume), dR_single_spi_trad0(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
           CALL add_stream_element(prp_stream, str_dR_single_spi_srad0(iplume), dR_single_spi_srad0(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

           CALL add_stream_element(prp_stream, str_aod_single_sp(iplume), ptr_aod_sp(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
           CALL add_stream_element(prp_stream, str_aod_single_bg(iplume), ptr_aod_bg(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

         END DO
        END IF
      IF (prprads) THEN
        CALL add_stream_element(prp_stream, 'dR_cdd_trads', dR_cdd_trads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_cdd_srads', dR_cdd_srads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_sp_trads', dR_sp_trads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_sp_srads', dR_sp_srads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_spd_trads', dR_spd_trads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_spd_srads', dR_spd_srads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_spi_trads', dR_spi_trads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_spi_srads', dR_spi_srads, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
      END IF
      IF (prpraf) THEN
        CALL add_stream_element(prp_stream, 'dR_cdd_traf0', dR_cdd_traf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_cdd_sraf0', dR_cdd_sraf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_cdd_trafs', dR_cdd_trafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_cdd_srafs', dR_cdd_srafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_sp_traf0', dR_sp_traf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_sp_sraf0', dR_sp_sraf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_sp_trafs', dR_sp_trafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_sp_srafs', dR_sp_srafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_spd_traf0', dR_spd_traf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_spd_sraf0', dR_spd_sraf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_spd_trafs', dR_spd_trafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_spd_srafs', dR_spd_srafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_spi_traf0', dR_spi_traf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_spi_sraf0', dR_spi_sraf0, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_spi_trafs', dR_spi_trafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
        CALL add_stream_element(prp_stream, 'dR_spi_srafs', dR_spi_srafs, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

        IF (prp_single_sp) THEN
          str_dR_single_sp_traf0 = (/'dR_sp1_traf0', 'dR_sp2_traf0', 'dR_sp3_traf0', 'dR_sp4_traf0', 'dR_sp5_traf0', 'dR_sp6_traf0', 'dR_sp7_traf0', 'dR_sp8_traf0', 'dR_sp9_traf0'/)
          str_dR_single_sp_sraf0 = (/'dR_sp1_sraf0', 'dR_sp2_sraf0', 'dR_sp3_sraf0', 'dR_sp4_sraf0', 'dR_sp5_sraf0', 'dR_sp6_sraf0', 'dR_sp7_sraf0', 'dR_sp8_sraf0', 'dR_sp9_sraf0'/)
          str_dR_single_spd_traf0 = (/'dR_spd1_traf0', 'dR_spd2_traf0', 'dR_spd3_traf0', 'dR_spd4_traf0', 'dR_spd5_traf0', 'dR_spd6_traf0', 'dR_spd7_traf0', 'dR_spd8_traf0', 'dR_spd9_traf0'/)
          str_dR_single_spd_sraf0 = (/'dR_spd1_sraf0', 'dR_spd2_sraf0', 'dR_spd3_sraf0', 'dR_spd4_sraf0', 'dR_spd5_sraf0', 'dR_spd6_sraf0', 'dR_spd7_sraf0', 'dR_spd8_sraf0', 'dR_spd9_sraf0'/)
          str_dR_single_spi_traf0 = (/'dR_spi1_traf0', 'dR_spi2_traf0', 'dR_spi3_traf0', 'dR_spi4_traf0', 'dR_spi5_traf0', 'dR_spi6_traf0', 'dR_spi7_traf0', 'dR_spi8_traf0', 'dR_spi9_traf0'/)
          str_dR_single_spi_sraf0 = (/'dR_spi1_sraf0', 'dR_spi2_sraf0', 'dR_spi3_sraf0', 'dR_spi4_sraf0', 'dR_spi5_sraf0', 'dR_spi6_sraf0', 'dR_spi7_sraf0', 'dR_spi8_sraf0', 'dR_spi9_sraf0'/)
          
          DO iplume=1,nplumes
            CALL add_stream_element(prp_stream, str_dR_single_sp_traf0(iplume), dR_single_sp_traf0(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
            CALL add_stream_element(prp_stream, str_dR_single_sp_sraf0(iplume), dR_single_sp_sraf0(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
            CALL add_stream_element(prp_stream, str_dR_single_spd_traf0(iplume), dR_single_spd_traf0(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
            CALL add_stream_element(prp_stream, str_dR_single_spd_sraf0(iplume), dR_single_spd_sraf0(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
            CALL add_stream_element(prp_stream, str_dR_single_spi_traf0(iplume), dR_single_spi_traf0(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
            CALL add_stream_element(prp_stream, str_dR_single_spi_sraf0(iplume), dR_single_spi_sraf0(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          END DO
        END IF
      END IF
      END IF
       
      IF (prp3doutput) THEN
          
          CALL add_stream_element(prp_stream, 'dQ_trad', dQ_trad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_srad', dQ_srad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_traf', dQ_traf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_sraf', dQ_sraf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          IF (prp_feedbacks) THEN
          CALL add_stream_element(prp_stream, 'dQ_tmp_trad', dQ_tmp_trad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_tmp_srad', dQ_tmp_srad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_tmp_traf', dQ_tmp_traf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_tmp_sraf', dQ_tmp_sraf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

          CALL add_stream_element(prp_stream, 'dQ_plk_trad', dQ_plk_trad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_plk_srad', dQ_plk_srad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_plk_traf', dQ_plk_traf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_plk_sraf', dQ_plk_sraf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

          CALL add_stream_element(prp_stream, 'dQ_lr_trad', dQ_lr_trad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_lr_srad', dQ_lr_srad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_lr_traf', dQ_lr_traf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_lr_sraf', dQ_lr_sraf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

          CALL add_stream_element(prp_stream, 'dQ_str_trad', dQ_str_trad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_str_srad', dQ_str_srad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_str_traf', dQ_str_traf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_str_sraf', dQ_str_sraf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

          CALL add_stream_element(prp_stream, 'dQ_vap_trad', dQ_vap_trad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_vap_srad', dQ_vap_srad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_vap_traf', dQ_vap_traf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_vap_sraf', dQ_vap_sraf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

          CALL add_stream_element(prp_stream, 'dQ_cld_trad', dQ_cld_trad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_cld_srad', dQ_cld_srad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_cld_traf', dQ_cld_traf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_cld_sraf', dQ_cld_sraf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

          CALL add_stream_element(prp_stream, 'dQ_aer_trad', dQ_aer_trad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_aer_srad', dQ_aer_srad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_aer_traf', dQ_aer_traf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_aer_sraf', dQ_aer_sraf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

          CALL add_stream_element(prp_stream, 'dQ_alb_trad', dQ_alb_trad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_alb_srad', dQ_alb_srad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_alb_traf', dQ_alb_traf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_alb_sraf', dQ_alb_sraf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)

          CALL add_stream_element(prp_stream, 'dQ_co2_trad', dQ_co2_trad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_co2_srad', dQ_co2_srad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_co2_traf', dQ_co2_traf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_co2_sraf', dQ_co2_sraf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          END IF
          IF (prp_sp) THEN
          CALL add_stream_element(prp_stream, 'dQ_cdd_trad', dQ_cdd_trad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_cdd_srad', dQ_cdd_srad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_sp_trad', dQ_sp_trad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_sp_srad', dQ_sp_srad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_spd_trad', dQ_spd_trad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_spd_srad', dQ_spd_srad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_spi_trad', dQ_spi_trad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_spi_srad', dQ_spi_srad, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          IF (prpraf) THEN
          CALL add_stream_element(prp_stream, 'dQ_cdd_traf', dQ_cdd_traf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_cdd_sraf', dQ_cdd_sraf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_sp_traf', dQ_sp_traf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_sp_sraf', dQ_sp_sraf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_spd_traf', dQ_spd_traf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_spd_sraf', dQ_spd_sraf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_spi_traf', dQ_spi_traf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          CALL add_stream_element(prp_stream, 'dQ_spi_sraf', dQ_spi_sraf, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
          END IF
          IF (prp3drad) THEN
            CALL add_stream_element(prp_stream, 'dR_cdd_trad', dR_cdd_trad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
            CALL add_stream_element(prp_stream, 'dR_cdd_srad', dR_cdd_srad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
            CALL add_stream_element(prp_stream, 'dR_cdd_traf', dR_cdd_traf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
            CALL add_stream_element(prp_stream, 'dR_cdd_sraf', dR_cdd_sraf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       
            CALL add_stream_element(prp_stream, 'dR_sp_trad', dR_sp_trad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
            CALL add_stream_element(prp_stream, 'dR_sp_srad', dR_sp_srad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
            CALL add_stream_element(prp_stream, 'dR_sp_traf', dR_sp_traf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
            CALL add_stream_element(prp_stream, 'dR_sp_sraf', dR_sp_sraf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       
            CALL add_stream_element(prp_stream, 'dR_spd_trad', dR_spd_trad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
            CALL add_stream_element(prp_stream, 'dR_spd_srad', dR_spd_srad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
            CALL add_stream_element(prp_stream, 'dR_spd_traf', dR_spd_traf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
            CALL add_stream_element(prp_stream, 'dR_spd_sraf', dR_spd_sraf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
       
            CALL add_stream_element(prp_stream, 'dR_spi_trad', dR_spi_trad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
            CALL add_stream_element(prp_stream, 'dR_spi_srad', dR_spi_srad, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
            CALL add_stream_element(prp_stream, 'dR_spi_traf', dR_spi_traf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)
            CALL add_stream_element(prp_stream, 'dR_spi_sraf', dR_spi_sraf, lpost=.FALSE., lrerun=.FALSE., laccu=.TRUE., klev=nlevp1)

            IF (prp_single_sp) THEN
              str_dR_single_sp_trad = (/'dR_sp1_trad', 'dR_sp2_trad', 'dR_sp3_trad', 'dR_sp4_trad', 'dR_sp5_trad', 'dR_sp6_trad', 'dR_sp7_trad', 'dR_sp8_trad', 'dR_sp9_trad'/)
              str_dR_single_sp_srad = (/'dR_sp1_srad', 'dR_sp2_srad', 'dR_sp3_srad', 'dR_sp4_srad', 'dR_sp5_srad', 'dR_sp6_srad', 'dR_sp7_srad', 'dR_sp8_srad', 'dR_sp9_srad'/)
              str_dR_single_spd_trad = (/'dR_spd1_trad', 'dR_spd2_trad', 'dR_spd3_trad', 'dR_spd4_trad', 'dR_spd5_trad', 'dR_spd6_trad', 'dR_spd7_trad', 'dR_spd8_trad', 'dR_spd9_trad'/)
              str_dR_single_spd_srad = (/'dR_spd1_srad', 'dR_spd2_srad', 'dR_spd3_srad', 'dR_spd4_srad', 'dR_spd5_srad', 'dR_spd6_srad', 'dR_spd7_srad', 'dR_spd8_srad', 'dR_spd9_srad'/)
              str_dR_single_spi_trad = (/'dR_spi1_trad', 'dR_spi2_trad', 'dR_spi3_trad', 'dR_spi4_trad', 'dR_spi5_trad', 'dR_spi6_trad', 'dR_spi7_trad', 'dR_spi8_trad', 'dR_spi9_trad'/)
              str_dR_single_spi_srad = (/'dR_spi1_srad', 'dR_spi2_srad', 'dR_spi3_srad', 'dR_spi4_srad', 'dR_spi5_srad', 'dR_spi6_srad', 'dR_spi7_srad', 'dR_spi8_srad', 'dR_spi9_srad'/)
  
              str_dR_single_sp_traf = (/'dR_sp1_traf', 'dR_sp2_traf', 'dR_sp3_traf', 'dR_sp4_traf', 'dR_sp5_traf', 'dR_sp6_traf', 'dR_sp7_traf', 'dR_sp8_traf', 'dR_sp9_traf'/)
              str_dR_single_sp_sraf = (/'dR_sp1_sraf', 'dR_sp2_sraf', 'dR_sp3_sraf', 'dR_sp4_sraf', 'dR_sp5_sraf', 'dR_sp6_sraf', 'dR_sp7_sraf', 'dR_sp8_sraf', 'dR_sp9_sraf'/)
              str_dR_single_spd_traf = (/'dR_spd1_traf', 'dR_spd2_traf', 'dR_spd3_traf', 'dR_spd4_traf', 'dR_spd5_traf', 'dR_spd6_traf', 'dR_spd7_traf', 'dR_spd8_traf', 'dR_spd9_traf'/)
              str_dR_single_spd_sraf = (/'dR_spd1_sraf', 'dR_spd2_sraf', 'dR_spd3_sraf', 'dR_spd4_sraf', 'dR_spd5_sraf', 'dR_spd6_sraf', 'dR_spd7_sraf', 'dR_spd8_sraf', 'dR_spd9_sraf'/)
              str_dR_single_spi_traf = (/'dR_spi1_traf', 'dR_spi2_traf', 'dR_spi3_traf', 'dR_spi4_traf', 'dR_spi5_traf', 'dR_spi6_traf', 'dR_spi7_traf', 'dR_spi8_traf', 'dR_spi9_traf'/)
              str_dR_single_spi_sraf = (/'dR_spi1_sraf', 'dR_spi2_sraf', 'dR_spi3_sraf', 'dR_spi4_sraf', 'dR_spi5_sraf', 'dR_spi6_sraf', 'dR_spi7_sraf', 'dR_spi8_sraf', 'dR_spi9_sraf'/)
  
              str_dQ_single_sp_trad = (/'dQ_sp1_trad', 'dQ_sp2_trad', 'dQ_sp3_trad', 'dQ_sp4_trad', 'dQ_sp5_trad', 'dQ_sp6_trad', 'dQ_sp7_trad', 'dQ_sp8_trad', 'dQ_sp9_trad'/)
              str_dQ_single_sp_srad = (/'dQ_sp1_srad', 'dQ_sp2_srad', 'dQ_sp3_srad', 'dQ_sp4_srad', 'dQ_sp5_srad', 'dQ_sp6_srad', 'dQ_sp7_srad', 'dQ_sp8_srad', 'dQ_sp9_srad'/)
              str_dQ_single_spd_trad = (/'dQ_spd1_trad', 'dQ_spd2_trad', 'dQ_spd3_trad', 'dQ_spd4_trad', 'dQ_spd5_trad', 'dQ_spd6_trad', 'dQ_spd7_trad', 'dQ_spd8_trad', 'dQ_spd9_trad'/)
              str_dQ_single_spd_srad = (/'dQ_spd1_srad', 'dQ_spd2_srad', 'dQ_spd3_srad', 'dQ_spd4_srad', 'dQ_spd5_srad', 'dQ_spd6_srad', 'dQ_spd7_srad', 'dQ_spd8_srad', 'dQ_spd9_srad'/)
              str_dQ_single_spi_trad = (/'dQ_spi1_trad', 'dQ_spi2_trad', 'dQ_spi3_trad', 'dQ_spi4_trad', 'dQ_spi5_trad', 'dQ_spi6_trad', 'dQ_spi7_trad', 'dQ_spi8_trad', 'dQ_spi9_trad'/)
              str_dQ_single_spi_srad = (/'dQ_spi1_srad', 'dQ_spi2_srad', 'dQ_spi3_srad', 'dQ_spi4_srad', 'dQ_spi5_srad', 'dQ_spi6_srad', 'dQ_spi7_srad', 'dQ_spi8_srad', 'dQ_spi9_srad'/)
  
              str_dQ_single_sp_traf = (/'dQ_sp1_traf', 'dQ_sp2_traf', 'dQ_sp3_traf', 'dQ_sp4_traf', 'dQ_sp5_traf', 'dQ_sp6_traf', 'dQ_sp7_traf', 'dQ_sp8_traf', 'dQ_sp9_traf'/)
              str_dQ_single_sp_sraf = (/'dQ_sp1_sraf', 'dQ_sp2_sraf', 'dQ_sp3_sraf', 'dQ_sp4_sraf', 'dQ_sp5_sraf', 'dQ_sp6_sraf', 'dQ_sp7_sraf', 'dQ_sp8_sraf', 'dQ_sp9_sraf'/)
              str_dQ_single_spd_traf = (/'dQ_spd1_traf', 'dQ_spd2_traf', 'dQ_spd3_traf', 'dQ_spd4_traf', 'dQ_spd5_traf', 'dQ_spd6_traf', 'dQ_spd7_traf', 'dQ_spd8_traf', 'dQ_spd9_traf'/)
              str_dQ_single_spd_sraf = (/'dQ_spd1_sraf', 'dQ_spd2_sraf', 'dQ_spd3_sraf', 'dQ_spd4_sraf', 'dQ_spd5_sraf', 'dQ_spd6_sraf', 'dQ_spd7_sraf', 'dQ_spd8_sraf', 'dQ_spd9_sraf'/)
              str_dQ_single_spi_traf = (/'dQ_spi1_traf', 'dQ_spi2_traf', 'dQ_spi3_traf', 'dQ_spi4_traf', 'dQ_spi5_traf', 'dQ_spi6_traf', 'dQ_spi7_traf', 'dQ_spi8_traf', 'dQ_spi9_traf'/)
              str_dQ_single_spi_sraf = (/'dQ_spi1_sraf', 'dQ_spi2_sraf', 'dQ_spi3_sraf', 'dQ_spi4_sraf', 'dQ_spi5_sraf', 'dQ_spi6_sraf', 'dQ_spi7_sraf', 'dQ_spi8_sraf', 'dQ_spi9_sraf'/)
              
              DO iplume=1,nplumes
                CALL add_stream_element(prp_stream, str_dR_single_sp_trad(iplume), dR_single_sp_trad(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dR_single_sp_srad(iplume), dR_single_sp_srad(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dR_single_spd_trad(iplume), dR_single_spd_trad(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dR_single_spd_srad(iplume), dR_single_spd_srad(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dR_single_spi_trad(iplume), dR_single_spi_trad(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dR_single_spi_srad(iplume), dR_single_spi_srad(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
  
                CALL add_stream_element(prp_stream, str_dQ_single_sp_trad(iplume), dQ_single_sp_trad(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dQ_single_sp_srad(iplume), dQ_single_sp_srad(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dQ_single_spd_trad(iplume), dQ_single_spd_trad(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dQ_single_spd_srad(iplume), dQ_single_spd_srad(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dQ_single_spi_trad(iplume), dQ_single_spi_trad(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dQ_single_spi_srad(iplume), dQ_single_spi_srad(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                IF (prpraf) THEN
                CALL add_stream_element(prp_stream, str_dR_single_sp_traf(iplume), dR_single_sp_traf(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dR_single_sp_sraf(iplume), dR_single_sp_sraf(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dR_single_spd_traf(iplume), dR_single_spd_traf(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dR_single_spd_sraf(iplume), dR_single_spd_sraf(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dR_single_spi_traf(iplume), dR_single_spi_traf(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dR_single_spi_sraf(iplume), dR_single_spi_sraf(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                
                CALL add_stream_element(prp_stream, str_dQ_single_sp_traf(iplume), dQ_single_sp_traf(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dQ_single_sp_sraf(iplume), dQ_single_sp_sraf(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dQ_single_spd_traf(iplume), dQ_single_spd_traf(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dQ_single_spd_sraf(iplume), dQ_single_spd_sraf(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dQ_single_spi_traf(iplume), dQ_single_spi_traf(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                CALL add_stream_element(prp_stream, str_dQ_single_spi_sraf(iplume), dQ_single_spi_sraf(iplume)%p, lpost=.TRUE., lrerun=.FALSE., laccu=.TRUE.)
                END IF
              END DO
            END IF
          END IF
        END IF    
      END IF
    END IF
  END SUBROUTINE construct_prp_stream


!-----------------------------------------------
!-----------------------------------------------

  SUBROUTINE destruct_prp_stream

    USE mo_memory_base,    ONLY: delete_stream

    IMPLICIT NONE
    
    CALL delete_stream(prp_stream)

  END SUBROUTINE destruct_prp_stream

!-----------------------------------------------
!-----------------------------------------------

  SUBROUTINE prp( &
         & iaero           ,kproma          ,kbdim           ,klev            ,& 
         & krow            ,ktrac           ,ktype           ,nb_sw           ,&
         & loland          ,loglac          ,cemiss          ,cos_mu0         ,&
         & pgeom1          ,alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,&
         & alb_nir_dif     ,pp_fl           ,pp_hl           ,pp_sfc          ,&
         & tk_fl           ,tk_hl           ,tk_sfc          ,xq_vap          ,&
         & xq_liq          ,xq_ice          ,cdnc            ,cld_frc         ,&
         & cld_cvr         ,xm_o3           ,xm_co2          ,xm_ch4          ,&
         & xm_n2o          ,xm_cfc          ,xm_o2           ,pxtm1             )

    !
    ! This subroutine is called from radiation. Depending on 
    ! prprecord, it will either (true) fetch the state at the
    ! radiation call, or (false) perform PRP diagnostics.
    !

    USE mo_decomposition,  ONLY: ldc => local_decomposition
    USE mo_mpi,            ONLY: p_parallel_io
    USE mo_time_conversion,ONLY: TC_convert, time_native, print_date, tc_get
    USE mo_time_control,   ONLY: current_date, lstart
    USE mo_exception,      ONLY: message
    USE mo_submodel,       ONLY: print_value
  !  USE mo_geoloc,         ONLY: amu0_x  !< actual zenith angle
    USE mo_physical_constants,      ONLY: grav, cpd
    USE mo_psrad_interface


    IMPLICIT NONE

    INTEGER,INTENT(IN)  ::             &
         kproma,                       & !< number of longitudes
         kbdim,                        & !< first dimension of 2-d arrays
         krow,                         & !< first dimension of 2-d arrays
         klev,                         & !< number of levels
         ktrac,                        & !< number of tracers
         ktype(kbdim),                 & !< type of convection
         nb_sw                           !< number of shortwave bands
    
    INTEGER,INTENT(IN) :: iaero(2)
    
    LOGICAL,INTENT(IN) ::              &
         loland(kbdim),                & !< land sea mask, land=.true.
         loglac(kbdim)                   !< glacier mask, glacier=.true.

    REAL(WP),INTENT(IN)  ::            &
         cemiss,                       & !< surface emissivity
         cos_mu0(kbdim),               & !< mu0 for solar zenith angle
         pgeom1(kbdim,klev),           & !< geopotential above ground
         alb_vis_dir(kbdim),           & !< surface albedo for vis range and dir light
         alb_nir_dir(kbdim),           & !< surface albedo for NIR range and dir light
         alb_vis_dif(kbdim),           & !< surface albedo for vis range and dif light
         alb_nir_dif(kbdim),           & !< surface albedo for NIR range and dif light
         pp_fl(kbdim,klev),            & !< full level pressure in Pa
         pp_hl(kbdim,klev+1),          & !< half level pressure in Pa
         pp_sfc(kbdim),                & !< surface pressure in Pa
         tk_fl(kbdim,klev),            & !< full level temperature in K
         tk_hl(kbdim,klev+1),          & !< half level temperature in K
         tk_sfc(kbdim),                & !< surface temperature in K
         xq_vap(kbdim,klev),           & !< specific humidity in g/g
         xq_liq(kbdim,klev),           & !< specific liquid water content
         xq_ice(kbdim,klev),           & !< specific ice content in g/g
         cdnc(kbdim,klev),             & !< cloud nuclei concentration
         cld_frc(kbdim,klev),          & !< fractional cloud cover
         cld_cvr(kbdim),               & !< total cloud cover in m2/m2
         xm_o3(kbdim,klev),            & !< o3 mass mixing ratio
         xm_co2(kbdim,klev),           & !< co2 mass mixing ratio
         xm_ch4(kbdim,klev),           & !< ch4 mass mixing ratio
         xm_n2o(kbdim,klev),           & !< n2o mass mixing ratio
         xm_cfc(kbdim,klev,2),         & !< cfc volume mixing ratio
         xm_o2(kbdim,klev),            & !< o2  mass mixing ratio
         pxtm1(kbdim,klev,ktrac)         !< tracer mass mixing ratios

    REAL(WP) :: loland_real(kbdim),          &
                loglac_real(kbdim)

    ! Variable to control the different scenarios to get the different effect of aerosols
    INTEGER :: plume_number, zplume_number, iplume
    INTEGER, DIMENSION(2) ::                        &
                xiaero=(/ 9, 10 /),                 &
                ziaero=(/ 3, 10 /),                 &
                iaero_sp_ref=(/ 5, 10 /),           &
                iaero_sp=(/ 9, 10 /),               &
                iaero_direct=(/ 10, 10 /),          &
                iaero_indirect=(/ 11, 10 /)

    !  Local variables for the stored state:
    INTEGER :: jl,                      & !< loop index
               zktype(kbdim)              !< type of convection
    LOGICAL :: zloland(kbdim),          & !< land sea mask, land=.true.
               zloglac(kbdim)             !< glacier mask, glacier=.true.

    REAL(WP)::                         &
        ! ziaero_real,                  &
         zcos_mu0(kbdim),              & !< mu0 for solar zenith angle
         zcemiss,                      & !< surface emissivity
         zpgeom1(kbdim,klev),          & !< geopotential above ground
         zpp_fl(kbdim,klev),           & !< full level pressure in Pa
         zpp_hl(kbdim,klev+1),         & !< half level pressure in Pa
         zpp_sfc(kbdim),               & !< surface pressure in Pa
         zktype_real(kbdim),           & !< type of convection in real type to be saved to the stream
         zloland_real(kbdim),          & !< land sea mask in real type to be saved to the stream
         zloglac_real(kbdim),          & !< glacier mask in real type to be saved to the stream

         zalb_vis_dir(kbdim),          & !< surface albedo for vis range and dir light
         zalb_nir_dir(kbdim),          & !< surface albedo for NIR range and dir light
         zalb_vis_dif(kbdim),          & !< surface albedo for vis range and dif light
         zalb_nir_dif(kbdim),          & !< surface albedo for NIR range and dif light
         ztk_fl(kbdim,klev),           & !< full level temperature in K
         ztk_hl(kbdim,klev+1),         & !< half level temperature in K
         y_tk_fl(kbdim,klev),        & !< full level temperature in K used in Planck, Lapse-rate and Stratosphere calculation
         y_tk_hl(kbdim,klev+1),      & !< half level temperature in K used in Planck, Lapse-rate and Stratosphere calculation
         ztk_sfc(kbdim),               & !< surface temperature in K
         zq_vap(kbdim,klev),           & !< specific humidity in g/g
         zq_liq(kbdim,klev),           & !< specific liquid water content
         zq_ice(kbdim,klev),           & !< specific ice content in g/g
         zcdnc(kbdim,klev),            & !< cloud nuclei concentration
         zcld_frc(kbdim,klev),         & !< fractional cloud cover
         zcld_cvr(kbdim),              & !< total cloud cover in m2/m2
         zm_co2(kbdim,klev),           & !< co2 mass mixing ratio
         zm_o3(kbdim,klev),            & !< o3 mass mixing ratio
         zm_ch4(kbdim,klev),           & !< ch4 mass mixing ratio
         zm_n2o(kbdim,klev),           & !< n2o mass mixing ratio
         zm_cfc(kbdim,klev,2),         & !< cfc volume mixing ratio
         zm_cfc1(kbdim,klev),          & ! divide cfc 4d into two var for stream
         zm_cfc2(kbdim,klev),          &
         zm_o2(kbdim,klev)               !< o2  mass mixing ratio

    integer                        :: year, month, day, hour, minute, second, k, i, j
    type(time_native)              :: my_date

    REAL(WP)::                         & 
         aflx_dnlw(kbdim,klev+1),      &
         aflx_dnsw(kbdim,klev+1),      &
         aflx_dnlw_clr(kbdim,klev+1),  &
         aflx_dnsw_clr(kbdim,klev+1),  &
         aflx_uplw(kbdim,klev+1),      &
         aflx_upsw(kbdim,klev+1),      &
         aflx_uplw_clr(kbdim,klev+1),  &
         aflx_upsw_clr(kbdim,klev+1),  &
         bflx_dnlw(kbdim,klev+1),      &
         bflx_dnsw(kbdim,klev+1),      &
         bflx_dnlw_clr(kbdim,klev+1),  &
         bflx_dnsw_clr(kbdim,klev+1),  &
         bflx_uplw(kbdim,klev+1),      &
         bflx_upsw(kbdim,klev+1),      & 
         bflx_uplw_clr(kbdim,klev+1),  &
         bflx_upsw_clr(kbdim,klev+1),  &
         cflx_dnlw(kbdim,klev+1),      &
         cflx_dnsw(kbdim,klev+1),      &
         cflx_dnlw_clr(kbdim,klev+1),  &
         cflx_dnsw_clr(kbdim,klev+1),  &
         cflx_uplw(kbdim,klev+1),      &
         cflx_upsw(kbdim,klev+1),      &
         cflx_uplw_clr(kbdim,klev+1),  &
         cflx_upsw_clr(kbdim,klev+1),  &
         dflx_dnlw(kbdim,klev+1),      &
         dflx_dnsw(kbdim,klev+1),      &
         dflx_dnlw_clr(kbdim,klev+1),  &
         dflx_dnsw_clr(kbdim,klev+1),  &
         dflx_uplw(kbdim,klev+1),      &
         dflx_upsw(kbdim,klev+1),      &
         dflx_uplw_clr(kbdim,klev+1),  &
         dflx_upsw_clr(kbdim,klev+1),  &
         toa_solar_irr(kbdim),         &
         nir_sfc(kbdim),               &
         nir_dff_sfc(kbdim),           &
         vis_sfc(kbdim),               &
         vis_dff_sfc(kbdim),           &
         dvis_sfc(kbdim),              &
         dpar_sfc(kbdim),              &
         par_dff_sfc(kbdim),           &
         dcm(kbdim,klev),              &
        ! zenith(kbdim),                &
         dtprps

!-----------------------------------------------
    IF (iaero(1)==0) THEN   ! not sure, but try to take into account the condition of multiple radiation calls in radiation code
      xiaero(1) = 0
      ziaero(1) = 0
      iaero_direct(1) = 0
      iaero_indirect(1) = 0
      iaero_sp(1) = 0
      iaero_sp_ref(1) = 0
    ELSE
      xiaero(1) = iaero(1)
      ziaero(1) = iaero_ref
    END IF

    IF (l_trigprp) THEN
    IF (prprecord) THEN
       
       DO jl = 1, kproma !< convert logical masks in binary real type to be saved to the stream
         loland_real(jl)            = merge(1.d0, 0.d0, loland(jl))    
         loglac_real(jl)            = merge(1.d0, 0.d0, loglac(jl))
       END DO

      ! ptr_iaero(krow)                               = real(iaero)    
       ptr_loland(1:ldc%nproma,krow)                 = loland_real
       ptr_loglac(1:ldc%nproma,krow)                 = loglac_real
     !  ptr_cemiss(krow)                              = cemiss
       ptr_pgeom1(1:ldc%nproma,1:ldc%nlev,krow)      = pgeom1
       ptr_pp_fl(1:ldc%nproma,1:ldc%nlev,krow)       = pp_fl
       ptr_pp_hl(1:ldc%nproma,1:ldc%nlev+1,krow)     = pp_hl
       ptr_pp_sfc(1:ldc%nproma,krow)                 = pp_sfc
       ptr_cos_mu0(1:ldc%nproma,krow)                = cos_mu0
       ptr_ktype(1:ldc%nproma,krow)                  = real(ktype) !< pointers must always be real?

       ptr_q_vap(1:ldc%nproma,1:ldc%nlev,krow)       = xq_vap
       ptr_tk_fl(1:ldc%nproma,1:ldc%nlev,krow)       = tk_fl
       ptr_tk_hl(1:ldc%nproma,1:ldc%nlev+1,krow)     = tk_hl
       ptr_tk_sfc(1:ldc%nproma,krow)                 = tk_sfc
       ptr_q_liq(1:ldc%nproma,1:ldc%nlev,krow)       = xq_liq
       ptr_q_ice(1:ldc%nproma,1:ldc%nlev,krow)       = xq_ice
       ptr_cdnc(1:ldc%nproma,1:ldc%nlev,krow)        = cdnc
       ptr_cld_frc(1:ldc%nproma,1:ldc%nlev,krow)     = cld_frc
       ptr_cld_cvr(1:ldc%nproma,krow)                = cld_cvr
       ptr_alb_vis_dir(1:ldc%nproma,krow)            = alb_vis_dir
       ptr_alb_nir_dir(1:ldc%nproma,krow)            = alb_nir_dir
       ptr_alb_vis_dif(1:ldc%nproma,krow)            = alb_vis_dif
       ptr_alb_nir_dif(1:ldc%nproma,krow)            = alb_nir_dif
       ptr_m_co2(1:ldc%nproma,1:ldc%nlev,krow)       = xm_co2
       ptr_m_o3(1:ldc%nproma,1:ldc%nlev,krow)        = xm_o3
       ptr_m_ch4(1:ldc%nproma,1:ldc%nlev,krow)       = xm_ch4
       ptr_m_n2o(1:ldc%nproma,1:ldc%nlev,krow)       = xm_n2o
       ptr_m_cfc1(1:ldc%nproma,1:ldc%nlev,krow)      = xm_cfc(:,:,1)
       ptr_m_cfc2(1:ldc%nproma,1:ldc%nlev,krow)      = xm_cfc(:,:,2)
       ptr_m_o2(1:ldc%nproma,1:ldc%nlev,krow)        = xm_o2      

       ! Create a local date:
       CALL TC_convert(current_date,my_date)
       CALL print_date(my_date,'Date')
       CALL tc_get(my_date, year, month, day, hour, minute, second)

       ptr_tod(1:ldc%nproma,krow)  = hour/24. + minute/(24.*60.) + second/(24.*3600.)
       ptr_dom(1:ldc%nproma,krow)  = day

    ELSEIF (.not.lstart) THEN

      ! ziaero_real = ptr_iaero(krow)                
       zloland_real = ptr_loland(1:ldc%nproma,krow) 
       zloglac_real = ptr_loglac(1:ldc%nproma,krow)

       DO jl = 1, kproma                     !< convert masks in real type back to logical (must exist a more sophisticated way, similar to merge)
         IF (abs(zloland_real(jl))==1) THEN
            zloland(jl) = .true.
         ELSE 
            zloland(jl) = .false.
         END IF
         IF (abs(zloglac_real(jl))==1) THEN
            zloglac(jl) = .true.
         ELSE 
            zloglac(jl) = .false.
         END IF
       END DO
      
      ! zcemiss      = ptr_cemiss(krow)
       zcemiss      = cemiss
       zpgeom1      = ptr_pgeom1(1:ldc%nproma,1:ldc%nlev,krow)
       zpp_fl       = ptr_pp_fl(1:ldc%nproma,1:ldc%nlev,krow)
       zpp_hl       = ptr_pp_hl(1:ldc%nproma,1:ldc%nlev+1,krow)
       zpp_sfc      = ptr_pp_sfc(1:ldc%nproma,krow)
       zcos_mu0     = ptr_cos_mu0(1:ldc%nproma,krow)
       zktype_real  = ptr_ktype(1:ldc%nproma,krow)
       zktype       = int(zktype_real)              !< pointer is real, but ktype is integer

       zalb_vis_dir = ptr_alb_vis_dir(1:ldc%nproma,krow)
       zalb_nir_dir = ptr_alb_nir_dir(1:ldc%nproma,krow)
       zalb_vis_dif = ptr_alb_vis_dif(1:ldc%nproma,krow)
       zalb_nir_dif = ptr_alb_nir_dif(1:ldc%nproma,krow)
       ztk_fl       = ptr_tk_fl(1:ldc%nproma,1:ldc%nlev,krow)
       ztk_hl       = ptr_tk_hl(1:ldc%nproma,1:ldc%nlev+1,krow)
       ztk_sfc      = ptr_tk_sfc(1:ldc%nproma,krow)
       zq_vap       = ptr_q_vap(1:ldc%nproma,1:ldc%nlev,krow)
       zq_liq       = ptr_q_liq(1:ldc%nproma,1:ldc%nlev,krow)
       zq_ice       = ptr_q_ice(1:ldc%nproma,1:ldc%nlev,krow)
       zcdnc        = ptr_cdnc(1:ldc%nproma,1:ldc%nlev,krow)
       zcld_frc     = ptr_cld_frc(1:ldc%nproma,1:ldc%nlev,krow)
       zcld_cvr     = ptr_cld_cvr(1:ldc%nproma,krow)
       zm_co2       = ptr_m_co2(1:ldc%nproma,1:ldc%nlev,krow)
       zm_o3        = ptr_m_o3(1:ldc%nproma,1:ldc%nlev,krow)
       zm_ch4       = ptr_m_ch4(1:ldc%nproma,1:ldc%nlev,krow)
       zm_n2o       = ptr_m_n2o(1:ldc%nproma,1:ldc%nlev,krow)
       zm_cfc1      = ptr_m_cfc1(1:ldc%nproma,1:ldc%nlev,krow)
       zm_cfc2      = ptr_m_cfc2(1:ldc%nproma,1:ldc%nlev,krow)
       zm_cfc(:,:,1) = zm_cfc1
       zm_cfc(:,:,2) = zm_cfc2
       zm_o2         = ptr_m_o2(1:ldc%nproma,1:ldc%nlev,krow) 


       CALL message('','Call psrad_interface with current state')

      ! zenith(:)                = amu0_x(1:ldc%nproma,krow)
       dtprps                   = prp_dt*3600._wp
       dcm(1:ldc%nproma,1:klev) = grav/cpd/(pp_hl(1:ldc%nproma,2:klev+1) - &
                                  pp_hl(1:ldc%nproma,1:klev))*24._wp*3600._wp*dtprps

       CALL psrad_interface( &
         & iaero           ,kproma          ,kbdim            ,klev             ,& 
         & krow            ,ktrac           ,ktype            ,nb_sw            ,& 
         & loland          ,loglac          ,cemiss           ,cos_mu0          ,& 
         & pgeom1          ,alb_vis_dir     ,alb_nir_dir      ,alb_vis_dif      ,&
         & alb_nir_dif     ,pp_fl           ,pp_hl            ,pp_sfc           ,&
         & tk_fl           ,tk_hl           ,tk_sfc           ,xq_vap           ,&
         & xq_liq          ,xq_ice          ,cdnc             ,cld_frc          ,&
         & cld_cvr         ,xm_o3           ,xm_co2           ,xm_ch4           ,&
         & xm_n2o          ,xm_cfc          ,xm_o2            ,pxtm1            ,&
         & aflx_uplw       ,aflx_uplw_clr   ,aflx_dnlw        ,aflx_dnlw_clr    ,&
         & aflx_upsw       ,aflx_upsw_clr   ,aflx_dnsw        ,aflx_dnsw_clr    ,& 
         & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
         & par_dff_sfc                                                          )

     !  CALL print_value('Temperature A: ', tk_fl(1,:))

       CALL message('','Call psrad_interface with saved state')

       CALL psrad_interface( &
         & ziaero          ,kproma          ,kbdim            ,klev             ,& 
         & krow            ,ktrac           ,zktype           ,nb_sw            ,&
         & zloland         ,zloglac         ,zcemiss          ,zcos_mu0         ,& 
         & zpgeom1         ,zalb_vis_dir    ,zalb_nir_dir     ,zalb_vis_dif     ,&
         & zalb_nir_dif    ,zpp_fl          ,zpp_hl           ,zpp_sfc          ,& 
         & ztk_fl          ,ztk_hl          ,ztk_sfc          ,zq_vap           ,&
         & zq_liq          ,zq_ice          ,zcdnc            ,zcld_frc         ,&
         & zcld_cvr        ,zm_o3           ,zm_co2           ,zm_ch4           ,&
         & zm_n2o          ,zm_cfc          ,zm_o2            ,pxtm1            ,&
         & bflx_uplw       ,bflx_uplw_clr   ,bflx_dnlw        ,bflx_dnlw_clr    ,&
         & bflx_upsw       ,bflx_upsw_clr   ,bflx_dnsw        ,bflx_dnsw_clr    ,&
         & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
         & par_dff_sfc                                                          )

       dR_trad0(1:ldc%nproma,krow)  = dR_trad0(1:ldc%nproma,krow) + &
                                      dtprps*((aflx_dnlw(:,1) - aflx_uplw(:,1)) - &
                                              (bflx_dnlw(:,1) - bflx_uplw(:,1)))
       dR_srad0(1:ldc%nproma,krow)  = dR_srad0(1:ldc%nproma,krow) + &
                                      dtprps*((aflx_dnsw(:,1) - aflx_upsw(:,1)) - &
                                              (bflx_dnsw(:,1) - bflx_upsw(:,1)))
       dR_trads(1:ldc%nproma,krow)  = dR_trads(1:ldc%nproma,krow) + &
                                      dtprps*((aflx_dnlw(:,klev+1) - aflx_uplw(:,klev+1)) - & 
                                              (bflx_dnlw(:,klev+1) - bflx_uplw(:,klev+1)))
       dR_srads(1:ldc%nproma,krow)  = dR_srads(1:ldc%nproma,krow) + &
                                      dtprps*((aflx_dnsw(:,klev+1) - aflx_upsw(:,klev+1)) - &
                                              (bflx_dnsw(:,klev+1) - bflx_upsw(:,klev+1)))
       dR_traf0(1:ldc%nproma,krow)  = dR_traf0(1:ldc%nproma,krow) + &
                                      dtprps*((aflx_dnlw_clr(:,1) - aflx_uplw_clr(:,1)) - &
                                              (bflx_dnlw_clr(:,1) - bflx_uplw_clr(:,1)))
       dR_sraf0(1:ldc%nproma,krow)  = dR_sraf0(1:ldc%nproma,krow) + &
                                      dtprps*((aflx_dnsw_clr(:,1) - aflx_upsw_clr(:,1)) - &
                                              (bflx_dnsw_clr(:,1) - bflx_upsw_clr(:,1)))
       dR_trafs(1:ldc%nproma,krow)  = dR_trafs(1:ldc%nproma,krow) + &
                                      dtprps*((aflx_dnlw_clr(:,klev+1) - aflx_uplw_clr(:,klev+1)) - &
                                              (bflx_dnlw_clr(:,klev+1) - bflx_uplw_clr(:,klev+1)))
       dR_srafs(1:ldc%nproma,krow)  = dR_srafs(1:ldc%nproma,krow) + &
                                      dtprps*((aflx_dnsw_clr(:,klev+1) - aflx_upsw_clr(:,klev+1)) - &
                                              (bflx_dnsw_clr(:,klev+1) - bflx_upsw_clr(:,klev+1)))

       dR_trad(1:ldc%nproma,1:klev+1,krow)  = dR_trad(1:ldc%nproma,1:klev+1,krow) + &
                                              dtprps*((aflx_dnlw(:,:) - aflx_uplw(:,:)) - &
						      (bflx_dnlw(:,:) - bflx_uplw(:,:)))
       dR_srad(1:ldc%nproma,1:klev+1,krow)  = dR_srad(1:ldc%nproma,1:klev+1,krow) + &
                                              dtprps*((aflx_dnsw(:,:) - aflx_upsw(:,:)) - &
						      (bflx_dnsw(:,:) - bflx_upsw(:,:)))
       dR_traf(1:ldc%nproma,1:klev+1,krow)  = dR_traf(1:ldc%nproma,1:klev+1,krow) + &
                                              dtprps*((aflx_dnlw_clr(:,:) - aflx_uplw_clr(:,:)) - &
						      (bflx_dnlw_clr(:,:) - bflx_uplw_clr(:,:)))
       dR_sraf(1:ldc%nproma,1:klev+1,krow)  = dR_sraf(1:ldc%nproma,1:klev+1,krow) + &
                                              dtprps*((aflx_dnsw_clr(:,:) - aflx_upsw_clr(:,:)) - &
						      (bflx_dnsw_clr(:,:) - bflx_upsw_clr(:,:)))

     IF (prp3doutput) THEN
       dQ_trad(1:ldc%nproma,1:klev,krow)  = dQ_trad(1:ldc%nproma,1:klev,krow) + &
                                             dcm(:,1:klev)*(((aflx_dnlw(:,1:klev) - aflx_uplw(:,1:klev)) - &
                                                             (aflx_dnlw(:,2:klev+1) - aflx_uplw(:,2:klev+1)) - &
                                                            ((bflx_dnlw(:,1:klev) - bflx_uplw(:,1:klev)) - &
                                                             (bflx_dnlw(:,2:klev+1) - bflx_dnlw(:,2:klev+1)))))
       dQ_srad(1:ldc%nproma,1:klev,krow)  = dQ_srad(1:ldc%nproma,1:klev,krow) + &
                                             dcm(:,1:klev)*(((aflx_dnsw(:,1:klev) - aflx_upsw(:,1:klev)) - &
                                                             (aflx_dnsw(:,2:klev+1) - aflx_upsw(:,2:klev+1)) - &
                                                            ((bflx_dnsw(:,1:klev) - bflx_upsw(:,1:klev)) - &
                                                             (bflx_dnsw(:,2:klev+1) - bflx_dnsw(:,2:klev+1)))))
       dQ_traf(1:ldc%nproma,1:klev,krow)  = dQ_traf(1:ldc%nproma,1:klev,krow) + &
                                             dcm(:,1:klev)*(((aflx_dnlw_clr(:,1:klev) - aflx_uplw_clr(:,1:klev)) - &
                                                             (aflx_dnlw_clr(:,2:klev+1) - aflx_uplw_clr(:,2:klev+1)) - &
                                                            ((bflx_dnlw_clr(:,1:klev) - bflx_uplw_clr(:,1:klev)) - &
                                                             (bflx_dnlw_clr(:,2:klev+1) - bflx_dnlw_clr(:,2:klev+1)))))
       dQ_sraf(1:ldc%nproma,1:klev,krow)  = dQ_sraf(1:ldc%nproma,1:klev,krow) + &
                                             dcm(:,1:klev)*(((aflx_dnsw_clr(:,1:klev) - aflx_upsw_clr(:,1:klev)) - &
                                                             (aflx_dnsw_clr(:,2:klev+1) - aflx_upsw_clr(:,2:klev+1)) - &
                                                            ((bflx_dnsw_clr(:,1:klev) - bflx_upsw_clr(:,1:klev)) - &
                                                             (bflx_dnsw_clr(:,2:klev+1) - bflx_dnsw_clr(:,2:klev+1)))))
     END IF

     IF (prp_feedbacks) THEN
       ! Calculate temperature feedback by averaging forward and backward calls:
      CALL message('','temperature feedback')
       CALL psrad_interface( &
         & iaero           ,kproma          ,kbdim            ,klev             ,& 
         & krow            ,ktrac           ,ktype            ,nb_sw            ,&
         & loland          ,loglac          ,cemiss           ,cos_mu0          ,&
         & pgeom1          ,alb_vis_dir     ,alb_nir_dir      ,alb_vis_dif      ,&
         & alb_nir_dif     ,pp_fl           ,pp_hl            ,pp_sfc           ,&
         & ztk_fl          ,ztk_hl          ,ztk_sfc          ,xq_vap           ,&
         & xq_liq          ,xq_ice          ,cdnc             ,cld_frc          ,&
         & cld_cvr         ,xm_o3           ,xm_co2           ,xm_ch4           ,&
         & xm_n2o          ,xm_cfc          ,xm_o2            ,pxtm1            ,&
         & cflx_uplw       ,cflx_uplw_clr   ,cflx_dnlw        ,cflx_dnlw_clr    ,&
         & cflx_upsw       ,cflx_upsw_clr   ,cflx_dnsw        ,cflx_dnsw_clr    ,&
         & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
         & par_dff_sfc                                                          ) 

       CALL psrad_interface( &
         & ziaero          ,kproma          ,kbdim            ,klev             ,& 
         & krow            ,ktrac           ,zktype           ,nb_sw            ,&
         & zloland         ,zloglac         ,zcemiss          ,zcos_mu0         ,&
         & zpgeom1         ,zalb_vis_dir    ,zalb_nir_dir     ,zalb_vis_dif     ,&
         & zalb_nir_dif    ,zpp_fl          ,zpp_hl           ,zpp_sfc          ,&
         & tk_fl           ,tk_hl           ,tk_sfc           ,zq_vap           ,&
         & zq_liq          ,zq_ice          ,zcdnc            ,zcld_frc         ,&
         & zcld_cvr        ,zm_o3           ,zm_co2           ,zm_ch4           ,&
         & zm_n2o          ,zm_cfc          ,zm_o2            ,pxtm1            ,&
         & dflx_uplw       ,dflx_uplw_clr   ,dflx_dnlw        ,dflx_dnlw_clr    ,&
         & dflx_upsw       ,dflx_upsw_clr   ,dflx_dnsw        ,dflx_dnsw_clr    ,&
         & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
         & par_dff_sfc                                                          ) 

       dR_tmp_trad0(1:ldc%nproma,krow)  = dR_tmp_trad0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw(:,1) - aflx_uplw(:,1)) - & 
							 (cflx_dnlw(:,1) - cflx_uplw(:,1)) + &
                                                         (dflx_dnlw(:,1) - dflx_uplw(:,1)) - &
                                                         (bflx_dnlw(:,1) - bflx_uplw(:,1))) 
       dR_tmp_srad0(1:ldc%nproma,krow)  = dR_tmp_srad0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw(:,1) - aflx_upsw(:,1)) - & 
							 (cflx_dnsw(:,1) - cflx_upsw(:,1)) + &
                                                         (dflx_dnsw(:,1) - dflx_upsw(:,1)) - &
                                                         (bflx_dnsw(:,1) - bflx_upsw(:,1))) 
       dR_tmp_trads(1:ldc%nproma,krow)  = dR_tmp_trads(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw(:,klev+1) - aflx_uplw(:,klev+1)) - & 
							 (cflx_dnlw(:,klev+1) - cflx_uplw(:,klev+1)) + &
                                                         (dflx_dnlw(:,klev+1) - dflx_uplw(:,klev+1)) - &
                                                         (bflx_dnlw(:,klev+1) - bflx_uplw(:,klev+1))) 
       dR_tmp_srads(1:ldc%nproma,krow)  = dR_tmp_srads(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw(:,klev+1) - aflx_upsw(:,klev+1)) - & 
							 (cflx_dnsw(:,klev+1) - cflx_upsw(:,klev+1)) + &
                                                         (dflx_dnsw(:,klev+1) - dflx_upsw(:,klev+1)) - &
                                                         (bflx_dnsw(:,klev+1) - bflx_upsw(:,klev+1))) 
       dR_tmp_traf0(1:ldc%nproma,krow)  = dR_tmp_traf0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw_clr(:,1) - aflx_uplw_clr(:,1)) - & 
							 (cflx_dnlw_clr(:,1) - cflx_uplw_clr(:,1)) + &
                                                         (dflx_dnlw_clr(:,1) - dflx_uplw_clr(:,1)) - &
                                                         (bflx_dnlw_clr(:,1) - bflx_uplw_clr(:,1))) 
       dR_tmp_sraf0(1:ldc%nproma,krow)  = dR_tmp_sraf0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw_clr(:,1) - aflx_upsw_clr(:,1)) - & 
							 (cflx_dnsw_clr(:,1) - cflx_upsw_clr(:,1)) + &
                                                         (dflx_dnsw_clr(:,1) - dflx_upsw_clr(:,1)) - &
                                                         (bflx_dnsw_clr(:,1) - bflx_upsw_clr(:,1))) 
       dR_tmp_trafs(1:ldc%nproma,krow)  = dR_tmp_trafs(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw_clr(:,klev+1) - aflx_uplw_clr(:,klev+1)) - & 
							 (cflx_dnlw_clr(:,klev+1) - cflx_uplw_clr(:,klev+1)) + &
                                                         (dflx_dnlw_clr(:,klev+1) - dflx_uplw_clr(:,klev+1)) - &
                                                         (bflx_dnlw_clr(:,klev+1) - bflx_uplw_clr(:,klev+1))) 
       dR_tmp_srafs(1:ldc%nproma,krow)  = dR_tmp_srafs(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw_clr(:,klev+1) - aflx_upsw_clr(:,klev+1)) - & 
							 (cflx_dnsw_clr(:,klev+1) - cflx_upsw_clr(:,klev+1)) + &
                                                         (dflx_dnsw_clr(:,klev+1) - dflx_upsw_clr(:,klev+1)) - &
                                                         (bflx_dnsw_clr(:,klev+1) - bflx_upsw_clr(:,klev+1))) 

       dR_tmp_trad(1:ldc%nproma,1:klev+1,krow)  = dR_tmp_trad(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnlw(:,:) - aflx_uplw(:,:)) - & 
								 (cflx_dnlw(:,:) - cflx_uplw(:,:)) + &
								 (dflx_dnlw(:,:) - dflx_uplw(:,:)) - &
								 (bflx_dnlw(:,:) - bflx_uplw(:,:))) 
       dR_tmp_srad(1:ldc%nproma,1:klev+1,krow)  = dR_tmp_srad(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnsw(:,:) - aflx_upsw(:,:)) - & 
								 (cflx_dnsw(:,:) - cflx_upsw(:,:)) + &
								 (dflx_dnsw(:,:) - dflx_upsw(:,:)) - &
								 (bflx_dnsw(:,:) - bflx_upsw(:,:))) 
       dR_tmp_traf(1:ldc%nproma,1:klev+1,krow)  = dR_tmp_traf(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnlw_clr(:,:) - aflx_uplw_clr(:,:)) - & 
								 (cflx_dnlw_clr(:,:) - cflx_uplw_clr(:,:)) + &
								 (dflx_dnlw_clr(:,:) - dflx_uplw_clr(:,:)) - &
								 (bflx_dnlw_clr(:,:) - bflx_uplw_clr(:,:))) 
       dR_tmp_sraf(1:ldc%nproma,1:klev+1,krow)  = dR_tmp_sraf(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnsw_clr(:,:) - aflx_upsw_clr(:,:)) - & 
								 (cflx_dnsw_clr(:,:) - cflx_upsw_clr(:,:)) + &
								 (dflx_dnsw_clr(:,:) - dflx_upsw_clr(:,:)) - &
								 (bflx_dnsw_clr(:,:) - bflx_upsw_clr(:,:)))

     IF (prp3doutput) THEN  
       dQ_tmp_trad(1:ldc%nproma,1:klev,krow)  = dQ_tmp_trad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnlw(:,1:klev) - aflx_uplw(:,1:klev)) - &
                                                   (aflx_dnlw(:,2:klev+1) - aflx_uplw(:,2:klev+1))) - &
                                                  ((cflx_dnlw(:,1:klev) - cflx_uplw(:,1:klev)) - & 
                                                   (cflx_dnlw(:,2:klev+1) - cflx_uplw(:,2:klev+1)))) - &
                                                 (((bflx_dnlw(:,1:klev) - bflx_uplw(:,1:klev)) - & 
                                                   (bflx_dnlw(:,2:klev+1) - bflx_uplw(:,2:klev+1))) - &
                                                  ((dflx_dnlw(:,1:klev) - dflx_uplw(:,1:klev)) - &
                                                   (dflx_dnlw(:,2:klev+1) - dflx_uplw(:,2:klev+1)))))
       dQ_tmp_srad(1:ldc%nproma,1:klev,krow)  = dQ_tmp_srad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnsw(:,1:klev) - aflx_upsw(:,1:klev)) - &
                                                   (aflx_dnsw(:,2:klev+1) - aflx_upsw(:,2:klev+1))) - &
                                                  ((cflx_dnsw(:,1:klev) - cflx_upsw(:,1:klev)) - & 
                                                   (cflx_dnsw(:,2:klev+1) - cflx_upsw(:,2:klev+1)))) - &
                                                 (((bflx_dnsw(:,1:klev) - bflx_upsw(:,1:klev)) - & 
                                                   (bflx_dnsw(:,2:klev+1) - bflx_upsw(:,2:klev+1))) - &
                                                  ((dflx_dnsw(:,1:klev) - dflx_upsw(:,1:klev)) - &
                                                   (dflx_dnsw(:,2:klev+1) - dflx_upsw(:,2:klev+1)))))
       dQ_tmp_traf(1:ldc%nproma,1:klev,krow)  = dQ_tmp_traf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnlw_clr(:,1:klev) - aflx_uplw_clr(:,1:klev)) - &
                                                   (aflx_dnlw_clr(:,2:klev+1) - aflx_uplw_clr(:,2:klev+1))) - &
                                                  ((cflx_dnlw_clr(:,1:klev) - cflx_uplw_clr(:,1:klev)) - & 
                                                   (cflx_dnlw_clr(:,2:klev+1) - cflx_uplw_clr(:,2:klev+1)))) - &
                                                 (((bflx_dnlw_clr(:,1:klev) - bflx_uplw_clr(:,1:klev)) - & 
                                                   (bflx_dnlw_clr(:,2:klev+1) - bflx_uplw_clr(:,2:klev+1))) - &
                                                  ((dflx_dnlw_clr(:,1:klev) - dflx_uplw_clr(:,1:klev)) - &
                                                   (dflx_dnlw_clr(:,2:klev+1) - dflx_uplw_clr(:,2:klev+1)))))
       dQ_tmp_sraf(1:ldc%nproma,1:klev,krow)  = dQ_tmp_sraf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnsw_clr(:,1:klev) - aflx_upsw_clr(:,1:klev)) - &
                                                   (aflx_dnsw_clr(:,2:klev+1) - aflx_upsw_clr(:,2:klev+1))) - &
                                                  ((cflx_dnsw_clr(:,1:klev) - cflx_upsw_clr(:,1:klev)) - & 
                                                   (cflx_dnsw_clr(:,2:klev+1) - cflx_upsw_clr(:,2:klev+1)))) - &
                                                 (((bflx_dnsw_clr(:,1:klev) - bflx_upsw_clr(:,1:klev)) - & 
                                                   (bflx_dnsw_clr(:,2:klev+1) - bflx_upsw_clr(:,2:klev+1))) - &
                                                  ((dflx_dnsw_clr(:,1:klev) - dflx_upsw_clr(:,1:klev)) - &
                                                   (dflx_dnsw_clr(:,2:klev+1) - dflx_upsw_clr(:,2:klev+1)))))
     END IF


       ! Calculate Planck feedback by averaging forward and backward calls:
      CALL message('','Planck feedback')
       DO k=1,kproma
          WHERE (tropo_p(k,krow).lt.pp_fl(k,:))
             y_tk_fl(k,:) = tk_fl(k,:) + (ztk_sfc(k)-tk_sfc(k))
          ELSEWHERE
             y_tk_fl(k,:) = tk_fl(k,:)
          END WHERE

          WHERE (tropo_p(k,krow).lt.pp_hl(k,:))
             y_tk_hl(k,:) = tk_hl(k,:) + (ztk_sfc(k)-tk_sfc(k))
          ELSEWHERE
             y_tk_hl(k,:) = tk_hl(k,:)
          END WHERE
       END DO

       CALL psrad_interface( &
         & iaero           ,kproma          ,kbdim            ,klev             ,& 
         & krow            ,ktrac           ,ktype            ,nb_sw            ,&
         & loland          ,loglac          ,cemiss           ,cos_mu0          ,&
         & pgeom1          ,alb_vis_dir     ,alb_nir_dir      ,alb_vis_dif      ,&
         & alb_nir_dif     ,pp_fl           ,pp_hl            ,pp_sfc           ,&
         & y_tk_fl         ,y_tk_hl         ,ztk_sfc          ,xq_vap           ,&
         & xq_liq          ,xq_ice          ,cdnc             ,cld_frc          ,&
         & cld_cvr         ,xm_o3           ,xm_co2           ,xm_ch4           ,&
         & xm_n2o          ,xm_cfc          ,xm_o2            ,pxtm1            ,&
         & cflx_uplw       ,cflx_uplw_clr   ,cflx_dnlw        ,cflx_dnlw_clr    ,&
         & cflx_upsw       ,cflx_upsw_clr   ,cflx_dnsw        ,cflx_dnsw_clr    ,&
         & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
         & par_dff_sfc                                                          ) 


       DO k=1,kproma
          WHERE (tropo_p(k,krow).lt.pp_fl(k,:))
             y_tk_fl(k,:) = ztk_fl(k,:) + (tk_sfc(k)-ztk_sfc(k))
          ELSEWHERE
             y_tk_fl(k,:) = ztk_fl(k,:)
          END WHERE

          WHERE (tropo_p(k,krow).lt.pp_hl(k,:))
             y_tk_hl(k,:) = ztk_hl(k,:) + (tk_sfc(k)-ztk_sfc(k))
          ELSEWHERE
             y_tk_hl(k,:) = ztk_hl(k,:)
          END WHERE
       END DO

       CALL psrad_interface( &
         & ziaero          ,kproma          ,kbdim            ,klev             ,& 
         & krow            ,ktrac           ,zktype           ,nb_sw            ,&
         & zloland         ,zloglac         ,zcemiss          ,zcos_mu0         ,& 
         & zpgeom1         ,zalb_vis_dir    ,zalb_nir_dir     ,zalb_vis_dif     ,&
         & zalb_nir_dif    ,zpp_fl          ,zpp_hl           ,zpp_sfc          ,&
         & y_tk_fl         ,y_tk_hl         ,tk_sfc           ,zq_vap           ,&
         & zq_liq          ,zq_ice          ,zcdnc            ,zcld_frc         ,&
         & zcld_cvr        ,zm_o3           ,zm_co2           ,zm_ch4           ,&
         & zm_n2o          ,zm_cfc          ,zm_o2            ,pxtm1            ,&
         & dflx_uplw       ,dflx_uplw_clr   ,dflx_dnlw        ,dflx_dnlw_clr    ,&
         & dflx_upsw       ,dflx_upsw_clr   ,dflx_dnsw        ,dflx_dnsw_clr    ,&
         & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
         & par_dff_sfc                                                          )

       dR_plk_trad0(1:ldc%nproma,krow)  = dR_plk_trad0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw(:,1) - aflx_uplw(:,1)) - & 
							 (cflx_dnlw(:,1) - cflx_uplw(:,1)) + &
                                                         (dflx_dnlw(:,1) - dflx_uplw(:,1)) - &
                                                         (bflx_dnlw(:,1) - bflx_uplw(:,1))) 
       dR_plk_srad0(1:ldc%nproma,krow)  = dR_plk_srad0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw(:,1) - aflx_upsw(:,1)) - & 
							 (cflx_dnsw(:,1) - cflx_upsw(:,1)) + &
                                                         (dflx_dnsw(:,1) - dflx_upsw(:,1)) - &
                                                         (bflx_dnsw(:,1) - bflx_upsw(:,1))) 
       dR_plk_trads(1:ldc%nproma,krow)  = dR_plk_trads(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw(:,klev+1) - aflx_uplw(:,klev+1)) - & 
							 (cflx_dnlw(:,klev+1) - cflx_uplw(:,klev+1)) + &
                                                         (dflx_dnlw(:,klev+1) - dflx_uplw(:,klev+1)) - &
                                                         (bflx_dnlw(:,klev+1) - bflx_uplw(:,klev+1))) 
       dR_plk_srads(1:ldc%nproma,krow)  = dR_plk_srads(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw(:,klev+1) - aflx_upsw(:,klev+1)) - & 
							 (cflx_dnsw(:,klev+1) - cflx_upsw(:,klev+1)) + &
                                                         (dflx_dnsw(:,klev+1) - dflx_upsw(:,klev+1)) - &
                                                         (bflx_dnsw(:,klev+1) - bflx_upsw(:,klev+1))) 
       dR_plk_traf0(1:ldc%nproma,krow)  = dR_plk_traf0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw_clr(:,1) - aflx_uplw_clr(:,1)) - & 
							 (cflx_dnlw_clr(:,1) - cflx_uplw_clr(:,1)) + &
                                                         (dflx_dnlw_clr(:,1) - dflx_uplw_clr(:,1)) - &
                                                         (bflx_dnlw_clr(:,1) - bflx_uplw_clr(:,1))) 
       dR_plk_sraf0(1:ldc%nproma,krow)  = dR_plk_sraf0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw_clr(:,1) - aflx_upsw_clr(:,1)) - & 
							 (cflx_dnsw_clr(:,1) - cflx_upsw_clr(:,1)) + &
                                                         (dflx_dnsw_clr(:,1) - dflx_upsw_clr(:,1)) - &
                                                         (bflx_dnsw_clr(:,1) - bflx_upsw_clr(:,1))) 
       dR_plk_trafs(1:ldc%nproma,krow)  = dR_plk_trafs(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw_clr(:,klev+1) - aflx_uplw_clr(:,klev+1)) - & 
							 (cflx_dnlw_clr(:,klev+1) - cflx_uplw_clr(:,klev+1)) + &
                                                         (dflx_dnlw_clr(:,klev+1) - dflx_uplw_clr(:,klev+1)) - &
                                                         (bflx_dnlw_clr(:,klev+1) - bflx_uplw_clr(:,klev+1))) 
       dR_plk_srafs(1:ldc%nproma,krow)  = dR_plk_srafs(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw_clr(:,klev+1) - aflx_upsw_clr(:,klev+1)) - & 
							 (cflx_dnsw_clr(:,klev+1) - cflx_upsw_clr(:,klev+1)) + &
                                                         (dflx_dnsw_clr(:,klev+1) - dflx_upsw_clr(:,klev+1)) - &
                                                         (bflx_dnsw_clr(:,klev+1) - bflx_upsw_clr(:,klev+1))) 

       dR_plk_trad(1:ldc%nproma,1:klev+1,krow)  = dR_plk_trad(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnlw(:,:) - aflx_uplw(:,:)) - & 
								 (cflx_dnlw(:,:) - cflx_uplw(:,:)) + &
								 (dflx_dnlw(:,:) - dflx_uplw(:,:)) - &
								 (bflx_dnlw(:,:) - bflx_uplw(:,:))) 
       dR_plk_srad(1:ldc%nproma,1:klev+1,krow)  = dR_plk_srad(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnsw(:,:) - aflx_upsw(:,:)) - & 
								 (cflx_dnsw(:,:) - cflx_upsw(:,:)) + &
								 (dflx_dnsw(:,:) - dflx_upsw(:,:)) - &
								 (bflx_dnsw(:,:) - bflx_upsw(:,:))) 
       dR_plk_traf(1:ldc%nproma,1:klev+1,krow)  = dR_plk_traf(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnlw_clr(:,:) - aflx_uplw_clr(:,:)) - & 
								 (cflx_dnlw_clr(:,:) - cflx_uplw_clr(:,:)) + &
								 (dflx_dnlw_clr(:,:) - dflx_uplw_clr(:,:)) - &
								 (bflx_dnlw_clr(:,:) - bflx_uplw_clr(:,:))) 
       dR_plk_sraf(1:ldc%nproma,1:klev+1,krow)  = dR_plk_sraf(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnsw_clr(:,:) - aflx_upsw_clr(:,:)) - & 
								 (cflx_dnsw_clr(:,:) - cflx_upsw_clr(:,:)) + &
								 (dflx_dnsw_clr(:,:) - dflx_upsw_clr(:,:)) - &
								 (bflx_dnsw_clr(:,:) - bflx_upsw_clr(:,:)))

     IF (prp3doutput) THEN
       dQ_plk_trad(1:ldc%nproma,1:klev,krow)  = dQ_plk_trad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnlw(:,1:klev) - aflx_uplw(:,1:klev)) - &
                                                   (aflx_dnlw(:,2:klev+1) - aflx_uplw(:,2:klev+1))) - &
                                                  ((cflx_dnlw(:,1:klev) - cflx_uplw(:,1:klev)) - & 
                                                   (cflx_dnlw(:,2:klev+1) - cflx_uplw(:,2:klev+1)))) - &
                                                 (((bflx_dnlw(:,1:klev) - bflx_uplw(:,1:klev)) - & 
                                                   (bflx_dnlw(:,2:klev+1) - bflx_uplw(:,2:klev+1))) - &
                                                  ((dflx_dnlw(:,1:klev) - dflx_uplw(:,1:klev)) - &
                                                   (dflx_dnlw(:,2:klev+1) - dflx_uplw(:,2:klev+1)))))
       dQ_plk_srad(1:ldc%nproma,1:klev,krow)  = dQ_plk_srad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnsw(:,1:klev) - aflx_upsw(:,1:klev)) - &
                                                   (aflx_dnsw(:,2:klev+1) - aflx_upsw(:,2:klev+1))) - &
                                                  ((cflx_dnsw(:,1:klev) - cflx_upsw(:,1:klev)) - & 
                                                   (cflx_dnsw(:,2:klev+1) - cflx_upsw(:,2:klev+1)))) - &
                                                 (((bflx_dnsw(:,1:klev) - bflx_upsw(:,1:klev)) - & 
                                                   (bflx_dnsw(:,2:klev+1) - bflx_upsw(:,2:klev+1))) - &
                                                  ((dflx_dnsw(:,1:klev) - dflx_upsw(:,1:klev)) - &
                                                   (dflx_dnsw(:,2:klev+1) - dflx_upsw(:,2:klev+1)))))
       dQ_plk_traf(1:ldc%nproma,1:klev,krow)  = dQ_plk_traf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnlw_clr(:,1:klev) - aflx_uplw_clr(:,1:klev)) - &
                                                   (aflx_dnlw_clr(:,2:klev+1) - aflx_uplw_clr(:,2:klev+1))) - &
                                                  ((cflx_dnlw_clr(:,1:klev) - cflx_uplw_clr(:,1:klev)) - & 
                                                   (cflx_dnlw_clr(:,2:klev+1) - cflx_uplw_clr(:,2:klev+1)))) - &
                                                 (((bflx_dnlw_clr(:,1:klev) - bflx_uplw_clr(:,1:klev)) - & 
                                                   (bflx_dnlw_clr(:,2:klev+1) - bflx_uplw_clr(:,2:klev+1))) - &
                                                  ((dflx_dnlw_clr(:,1:klev) - dflx_uplw_clr(:,1:klev)) - &
                                                   (dflx_dnlw_clr(:,2:klev+1) - dflx_uplw_clr(:,2:klev+1)))))
       dQ_plk_sraf(1:ldc%nproma,1:klev,krow)  = dQ_plk_sraf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnsw_clr(:,1:klev) - aflx_upsw_clr(:,1:klev)) - &
                                                   (aflx_dnsw_clr(:,2:klev+1) - aflx_upsw_clr(:,2:klev+1))) - &
                                                  ((cflx_dnsw_clr(:,1:klev) - cflx_upsw_clr(:,1:klev)) - & 
                                                   (cflx_dnsw_clr(:,2:klev+1) - cflx_upsw_clr(:,2:klev+1)))) - &
                                                 (((bflx_dnsw_clr(:,1:klev) - bflx_upsw_clr(:,1:klev)) - & 
                                                   (bflx_dnsw_clr(:,2:klev+1) - bflx_upsw_clr(:,2:klev+1))) - &
                                                  ((dflx_dnsw_clr(:,1:klev) - dflx_upsw_clr(:,1:klev)) - &
                                                   (dflx_dnsw_clr(:,2:klev+1) - dflx_upsw_clr(:,2:klev+1)))))
     END IF



       ! Calculate Lapse-rate feedback by averaging forward and backward calls:
       CALL message('','laspe-rate feedback')
       DO k=1,kproma
          WHERE (tropo_p(k,krow).lt.pp_fl(k,:))
             y_tk_fl(k,:) = ztk_fl(k,:) - (ztk_sfc(k)-tk_sfc(k))
          ELSEWHERE
             y_tk_fl(k,:) = tk_fl(k,:)
          END WHERE

          WHERE (tropo_p(k,krow).lt.pp_hl(k,:))
             y_tk_hl(k,:) = ztk_hl(k,:) - (ztk_sfc(k)-tk_sfc(k))
          ELSEWHERE
             y_tk_hl(k,:) = tk_hl(k,:)
          END WHERE
       END DO

       CALL psrad_interface( &
         & iaero           ,kproma          ,kbdim            ,klev             ,& 
         & krow            ,ktrac           ,ktype            ,nb_sw            ,&
         & loland          ,loglac          ,cemiss           ,cos_mu0          ,& 
         & pgeom1          ,alb_vis_dir     ,alb_nir_dir      ,alb_vis_dif      ,&
         & alb_nir_dif     ,pp_fl           ,pp_hl            ,pp_sfc           ,&
         & y_tk_fl         ,y_tk_hl         ,tk_sfc           ,xq_vap           ,&
         & xq_liq          ,xq_ice          ,cdnc             ,cld_frc          ,&
         & cld_cvr         ,xm_o3           ,xm_co2           ,xm_ch4           ,&
         & xm_n2o          ,xm_cfc          ,xm_o2            ,pxtm1            ,&
         & cflx_uplw       ,cflx_uplw_clr   ,cflx_dnlw        ,cflx_dnlw_clr    ,&
         & cflx_upsw       ,cflx_upsw_clr   ,cflx_dnsw        ,cflx_dnsw_clr    ,&
         & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
         & par_dff_sfc                                                          )

       DO k=1,kproma
          WHERE (tropo_p(k,krow).lt.pp_fl(k,:))
             y_tk_fl(k,:) = tk_fl(k,:) - (tk_sfc(k)-ztk_sfc(k))
          ELSEWHERE
             y_tk_fl(k,:) = ztk_fl(k,:)
          END WHERE

          WHERE (tropo_p(k,krow).lt.pp_hl(k,:))
             y_tk_hl(k,:) = tk_hl(k,:) - (tk_sfc(k)-ztk_sfc(k))
          ELSEWHERE
             y_tk_hl(k,:) = ztk_hl(k,:)
          END WHERE
       END DO

       CALL psrad_interface( &
         & ziaero          ,kproma          ,kbdim            ,klev             ,& 
         & krow            ,ktrac           ,zktype           ,nb_sw            ,&
         & zloland         ,zloglac         ,zcemiss          ,zcos_mu0         ,& 
         & zpgeom1         ,zalb_vis_dir    ,zalb_nir_dir     ,zalb_vis_dif     ,&
         & zalb_nir_dif    ,zpp_fl          ,zpp_hl           ,zpp_sfc          ,&
         & y_tk_fl         ,y_tk_hl         ,ztk_sfc          ,zq_vap           ,&
         & zq_liq          ,zq_ice          ,zcdnc            ,zcld_frc         ,&
         & zcld_cvr        ,zm_o3           ,zm_co2           ,zm_ch4           ,&
         & zm_n2o          ,zm_cfc          ,zm_o2            ,pxtm1            ,&
         & dflx_uplw       ,dflx_uplw_clr   ,dflx_dnlw        ,dflx_dnlw_clr    ,&
         & dflx_upsw       ,dflx_upsw_clr   ,dflx_dnsw        ,dflx_dnsw_clr    ,&
         & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
         & par_dff_sfc                                                          ) 

       dR_lr_trad0(1:ldc%nproma,krow)  = dR_lr_trad0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw(:,1) - aflx_uplw(:,1)) - & 
							 (cflx_dnlw(:,1) - cflx_uplw(:,1)) + &
                                                         (dflx_dnlw(:,1) - dflx_uplw(:,1)) - &
                                                         (bflx_dnlw(:,1) - bflx_uplw(:,1))) 
       dR_lr_srad0(1:ldc%nproma,krow)  = dR_lr_srad0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw(:,1) - aflx_upsw(:,1)) - & 
							 (cflx_dnsw(:,1) - cflx_upsw(:,1)) + &
                                                         (dflx_dnsw(:,1) - dflx_upsw(:,1)) - &
                                                         (bflx_dnsw(:,1) - bflx_upsw(:,1))) 
       dR_lr_trads(1:ldc%nproma,krow)  = dR_lr_trads(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw(:,klev+1) - aflx_uplw(:,klev+1)) - & 
							 (cflx_dnlw(:,klev+1) - cflx_uplw(:,klev+1)) + &
                                                         (dflx_dnlw(:,klev+1) - dflx_uplw(:,klev+1)) - &
                                                         (bflx_dnlw(:,klev+1) - bflx_uplw(:,klev+1))) 
       dR_lr_srads(1:ldc%nproma,krow)  = dR_lr_srads(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw(:,klev+1) - aflx_upsw(:,klev+1)) - & 
							 (cflx_dnsw(:,klev+1) - cflx_upsw(:,klev+1)) + &
                                                         (dflx_dnsw(:,klev+1) - dflx_upsw(:,klev+1)) - &
                                                         (bflx_dnsw(:,klev+1) - bflx_upsw(:,klev+1))) 
       dR_lr_traf0(1:ldc%nproma,krow)  = dR_lr_traf0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw_clr(:,1) - aflx_uplw_clr(:,1)) - & 
							 (cflx_dnlw_clr(:,1) - cflx_uplw_clr(:,1)) + &
                                                         (dflx_dnlw_clr(:,1) - dflx_uplw_clr(:,1)) - &
                                                         (bflx_dnlw_clr(:,1) - bflx_uplw_clr(:,1))) 
       dR_lr_sraf0(1:ldc%nproma,krow)  = dR_lr_sraf0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw_clr(:,1) - aflx_upsw_clr(:,1)) - & 
							 (cflx_dnsw_clr(:,1) - cflx_upsw_clr(:,1)) + &
                                                         (dflx_dnsw_clr(:,1) - dflx_upsw_clr(:,1)) - &
                                                         (bflx_dnsw_clr(:,1) - bflx_upsw_clr(:,1))) 
       dR_lr_trafs(1:ldc%nproma,krow)  = dR_lr_trafs(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw_clr(:,klev+1) - aflx_uplw_clr(:,klev+1)) - & 
							 (cflx_dnlw_clr(:,klev+1) - cflx_uplw_clr(:,klev+1)) + &
                                                         (dflx_dnlw_clr(:,klev+1) - dflx_uplw_clr(:,klev+1)) - &
                                                         (bflx_dnlw_clr(:,klev+1) - bflx_uplw_clr(:,klev+1))) 
       dR_lr_srafs(1:ldc%nproma,krow)  = dR_lr_srafs(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw_clr(:,klev+1) - aflx_upsw_clr(:,klev+1)) - & 
							 (cflx_dnsw_clr(:,klev+1) - cflx_upsw_clr(:,klev+1)) + &
                                                         (dflx_dnsw_clr(:,klev+1) - dflx_upsw_clr(:,klev+1)) - &
                                                         (bflx_dnsw_clr(:,klev+1) - bflx_upsw_clr(:,klev+1))) 

       dR_lr_trad(1:ldc%nproma,1:klev+1,krow)  = dR_lr_trad(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnlw(:,:) - aflx_uplw(:,:)) - & 
								 (cflx_dnlw(:,:) - cflx_uplw(:,:)) + &
								 (dflx_dnlw(:,:) - dflx_uplw(:,:)) - &
								 (bflx_dnlw(:,:) - bflx_uplw(:,:))) 
       dR_lr_srad(1:ldc%nproma,1:klev+1,krow)  = dR_lr_srad(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnsw(:,:) - aflx_upsw(:,:)) - & 
								 (cflx_dnsw(:,:) - cflx_upsw(:,:)) + &
								 (dflx_dnsw(:,:) - dflx_upsw(:,:)) - &
								 (bflx_dnsw(:,:) - bflx_upsw(:,:))) 
       dR_lr_traf(1:ldc%nproma,1:klev+1,krow)  = dR_lr_traf(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnlw_clr(:,:) - aflx_uplw_clr(:,:)) - & 
								 (cflx_dnlw_clr(:,:) - cflx_uplw_clr(:,:)) + &
								 (dflx_dnlw_clr(:,:) - dflx_uplw_clr(:,:)) - &
								 (bflx_dnlw_clr(:,:) - bflx_uplw_clr(:,:))) 
       dR_lr_sraf(1:ldc%nproma,1:klev+1,krow)  = dR_lr_sraf(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnsw_clr(:,:) - aflx_upsw_clr(:,:)) - & 
								 (cflx_dnsw_clr(:,:) - cflx_upsw_clr(:,:)) + &
								 (dflx_dnsw_clr(:,:) - dflx_upsw_clr(:,:)) - &
								 (bflx_dnsw_clr(:,:) - bflx_upsw_clr(:,:)))

     IF (prp3doutput) THEN
       dQ_lr_trad(1:ldc%nproma,1:klev,krow)  = dQ_lr_trad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnlw(:,1:klev) - aflx_uplw(:,1:klev)) - &
                                                   (aflx_dnlw(:,2:klev+1) - aflx_uplw(:,2:klev+1))) - &
                                                  ((cflx_dnlw(:,1:klev) - cflx_uplw(:,1:klev)) - & 
                                                   (cflx_dnlw(:,2:klev+1) - cflx_uplw(:,2:klev+1)))) - &
                                                 (((bflx_dnlw(:,1:klev) - bflx_uplw(:,1:klev)) - & 
                                                   (bflx_dnlw(:,2:klev+1) - bflx_uplw(:,2:klev+1))) - &
                                                  ((dflx_dnlw(:,1:klev) - dflx_uplw(:,1:klev)) - &
                                                   (dflx_dnlw(:,2:klev+1) - dflx_uplw(:,2:klev+1)))))
       dQ_lr_srad(1:ldc%nproma,1:klev,krow)  = dQ_lr_srad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnsw(:,1:klev) - aflx_upsw(:,1:klev)) - &
                                                   (aflx_dnsw(:,2:klev+1) - aflx_upsw(:,2:klev+1))) - &
                                                  ((cflx_dnsw(:,1:klev) - cflx_upsw(:,1:klev)) - & 
                                                   (cflx_dnsw(:,2:klev+1) - cflx_upsw(:,2:klev+1)))) - &
                                                 (((bflx_dnsw(:,1:klev) - bflx_upsw(:,1:klev)) - & 
                                                   (bflx_dnsw(:,2:klev+1) - bflx_upsw(:,2:klev+1))) - &
                                                  ((dflx_dnsw(:,1:klev) - dflx_upsw(:,1:klev)) - &
                                                   (dflx_dnsw(:,2:klev+1) - dflx_upsw(:,2:klev+1)))))
       dQ_lr_traf(1:ldc%nproma,1:klev,krow)  = dQ_lr_traf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnlw_clr(:,1:klev) - aflx_uplw_clr(:,1:klev)) - &
                                                   (aflx_dnlw_clr(:,2:klev+1) - aflx_uplw_clr(:,2:klev+1))) - &
                                                  ((cflx_dnlw_clr(:,1:klev) - cflx_uplw_clr(:,1:klev)) - & 
                                                   (cflx_dnlw_clr(:,2:klev+1) - cflx_uplw_clr(:,2:klev+1)))) - &
                                                 (((bflx_dnlw_clr(:,1:klev) - bflx_uplw_clr(:,1:klev)) - & 
                                                   (bflx_dnlw_clr(:,2:klev+1) - bflx_uplw_clr(:,2:klev+1))) - &
                                                  ((dflx_dnlw_clr(:,1:klev) - dflx_uplw_clr(:,1:klev)) - &
                                                   (dflx_dnlw_clr(:,2:klev+1) - dflx_uplw_clr(:,2:klev+1)))))
       dQ_lr_sraf(1:ldc%nproma,1:klev,krow)  = dQ_lr_sraf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnsw_clr(:,1:klev) - aflx_upsw_clr(:,1:klev)) - &
                                                   (aflx_dnsw_clr(:,2:klev+1) - aflx_upsw_clr(:,2:klev+1))) - &
                                                  ((cflx_dnsw_clr(:,1:klev) - cflx_upsw_clr(:,1:klev)) - & 
                                                   (cflx_dnsw_clr(:,2:klev+1) - cflx_upsw_clr(:,2:klev+1)))) - &
                                                 (((bflx_dnsw_clr(:,1:klev) - bflx_upsw_clr(:,1:klev)) - & 
                                                   (bflx_dnsw_clr(:,2:klev+1) - bflx_upsw_clr(:,2:klev+1))) - &
                                                  ((dflx_dnsw_clr(:,1:klev) - dflx_upsw_clr(:,1:klev)) - &
                                                   (dflx_dnsw_clr(:,2:klev+1) - dflx_upsw_clr(:,2:klev+1)))))
      END IF


       ! Calculate Stratospheric temperature feedback by averaging forward and backward calls:
      CALL message('','stratospheric temperature feedback')
       DO k=1,kproma
          WHERE (tropo_p(k,krow).lt.pp_fl(k,:))
             y_tk_fl(k,:) = tk_fl(k,:)
          ELSEWHERE
             y_tk_fl(k,:) = ztk_fl(k,:)
          END WHERE

          WHERE (tropo_p(k,krow).lt.pp_hl(k,:))
             y_tk_hl(k,:) = tk_hl(k,:)
          ELSEWHERE
             y_tk_hl(k,:) = ztk_hl(k,:)
          END WHERE
       END DO

       CALL psrad_interface( &
         & iaero           ,kproma          ,kbdim            ,klev             ,& 
         & krow            ,ktrac           ,ktype            ,nb_sw            ,&
         & loland          ,loglac          ,cemiss           ,cos_mu0          ,& 
         & pgeom1          ,alb_vis_dir     ,alb_nir_dir      ,alb_vis_dif      ,&
         & alb_nir_dif     ,pp_fl           ,pp_hl            ,pp_sfc           ,&
         & y_tk_fl         ,y_tk_hl         ,tk_sfc           ,xq_vap           ,&
         & xq_liq          ,xq_ice          ,cdnc             ,cld_frc          ,&
         & cld_cvr         ,xm_o3           ,xm_co2           ,xm_ch4           ,&
         & xm_n2o          ,xm_cfc          ,xm_o2            ,pxtm1            ,&
         & cflx_uplw       ,cflx_uplw_clr   ,cflx_dnlw        ,cflx_dnlw_clr    ,&
         & cflx_upsw       ,cflx_upsw_clr   ,cflx_dnsw        ,cflx_dnsw_clr    ,&
         & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
         & par_dff_sfc                                                          )

       DO k=1,kproma
          WHERE (tropo_p(k,krow).lt.pp_fl(k,:))
             y_tk_fl(k,:) = ztk_fl(k,:)
          ELSEWHERE
             y_tk_fl(k,:) = tk_fl(k,:)
          END WHERE

          WHERE (tropo_p(k,krow).lt.pp_hl(k,:))
             y_tk_hl(k,:) = ztk_hl(k,:)
          ELSEWHERE
             y_tk_hl(k,:) = tk_hl(k,:)
          END WHERE
       END DO

       CALL psrad_interface( &
         & ziaero          ,kproma          ,kbdim            ,klev             ,& 
         & krow            ,ktrac           ,zktype           ,nb_sw            ,&
         & zloland         ,zloglac         ,zcemiss          ,zcos_mu0         ,& 
         & zpgeom1         ,zalb_vis_dir    ,zalb_nir_dir     ,zalb_vis_dif     ,&
         & zalb_nir_dif    ,zpp_fl          ,zpp_hl           ,zpp_sfc          ,& 
         & y_tk_fl         ,y_tk_hl         ,ztk_sfc          ,zq_vap           ,&
         & zq_liq          ,zq_ice          ,zcdnc            ,zcld_frc         ,&
         & zcld_cvr        ,zm_o3           ,zm_co2           ,zm_ch4           ,&
         & zm_n2o          ,zm_cfc          ,zm_o2            ,pxtm1            ,&
         & dflx_uplw       ,dflx_uplw_clr   ,dflx_dnlw        ,dflx_dnlw_clr    ,&
         & dflx_upsw       ,dflx_upsw_clr   ,dflx_dnsw        ,dflx_dnsw_clr    ,&
         & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
         & par_dff_sfc                                                          )

       dR_str_trad0(1:ldc%nproma,krow)  = dR_str_trad0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw(:,1) - aflx_uplw(:,1)) - & 
							 (cflx_dnlw(:,1) - cflx_uplw(:,1)) + &
                                                         (dflx_dnlw(:,1) - dflx_uplw(:,1)) - &
                                                         (bflx_dnlw(:,1) - bflx_uplw(:,1))) 
       dR_str_srad0(1:ldc%nproma,krow)  = dR_str_srad0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw(:,1) - aflx_upsw(:,1)) - & 
							 (cflx_dnsw(:,1) - cflx_upsw(:,1)) + &
                                                         (dflx_dnsw(:,1) - dflx_upsw(:,1)) - &
                                                         (bflx_dnsw(:,1) - bflx_upsw(:,1))) 
       dR_str_trads(1:ldc%nproma,krow)  = dR_str_trads(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw(:,klev+1) - aflx_uplw(:,klev+1)) - & 
							 (cflx_dnlw(:,klev+1) - cflx_uplw(:,klev+1)) + &
                                                         (dflx_dnlw(:,klev+1) - dflx_uplw(:,klev+1)) - &
                                                         (bflx_dnlw(:,klev+1) - bflx_uplw(:,klev+1))) 
       dR_str_srads(1:ldc%nproma,krow)  = dR_str_srads(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw(:,klev+1) - aflx_upsw(:,klev+1)) - & 
							 (cflx_dnsw(:,klev+1) - cflx_upsw(:,klev+1)) + &
                                                         (dflx_dnsw(:,klev+1) - dflx_upsw(:,klev+1)) - &
                                                         (bflx_dnsw(:,klev+1) - bflx_upsw(:,klev+1))) 
       dR_str_traf0(1:ldc%nproma,krow)  = dR_str_traf0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw_clr(:,1) - aflx_uplw_clr(:,1)) - & 
							 (cflx_dnlw_clr(:,1) - cflx_uplw_clr(:,1)) + &
                                                         (dflx_dnlw_clr(:,1) - dflx_uplw_clr(:,1)) - &
                                                         (bflx_dnlw_clr(:,1) - bflx_uplw_clr(:,1))) 
       dR_str_sraf0(1:ldc%nproma,krow)  = dR_str_sraf0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw_clr(:,1) - aflx_upsw_clr(:,1)) - & 
							 (cflx_dnsw_clr(:,1) - cflx_upsw_clr(:,1)) + &
                                                         (dflx_dnsw_clr(:,1) - dflx_upsw_clr(:,1)) - &
                                                         (bflx_dnsw_clr(:,1) - bflx_upsw_clr(:,1))) 
       dR_str_trafs(1:ldc%nproma,krow)  = dR_str_trafs(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw_clr(:,klev+1) - aflx_uplw_clr(:,klev+1)) - & 
							 (cflx_dnlw_clr(:,klev+1) - cflx_uplw_clr(:,klev+1)) + &
                                                         (dflx_dnlw_clr(:,klev+1) - dflx_uplw_clr(:,klev+1)) - &
                                                         (bflx_dnlw_clr(:,klev+1) - bflx_uplw_clr(:,klev+1))) 
       dR_str_srafs(1:ldc%nproma,krow)  = dR_str_srafs(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw_clr(:,klev+1) - aflx_upsw_clr(:,klev+1)) - & 
							 (cflx_dnsw_clr(:,klev+1) - cflx_upsw_clr(:,klev+1)) + &
                                                         (dflx_dnsw_clr(:,klev+1) - dflx_upsw_clr(:,klev+1)) - &
                                                         (bflx_dnsw_clr(:,klev+1) - bflx_upsw_clr(:,klev+1))) 

       dR_str_trad(1:ldc%nproma,1:klev+1,krow)  = dR_str_trad(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnlw(:,:) - aflx_uplw(:,:)) - & 
								 (cflx_dnlw(:,:) - cflx_uplw(:,:)) + &
								 (dflx_dnlw(:,:) - dflx_uplw(:,:)) - &
								 (bflx_dnlw(:,:) - bflx_uplw(:,:))) 
       dR_str_srad(1:ldc%nproma,1:klev+1,krow)  = dR_str_srad(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnsw(:,:) - aflx_upsw(:,:)) - & 
								 (cflx_dnsw(:,:) - cflx_upsw(:,:)) + &
								 (dflx_dnsw(:,:) - dflx_upsw(:,:)) - &
								 (bflx_dnsw(:,:) - bflx_upsw(:,:))) 
       dR_str_traf(1:ldc%nproma,1:klev+1,krow)  = dR_str_traf(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnlw_clr(:,:) - aflx_uplw_clr(:,:)) - & 
								 (cflx_dnlw_clr(:,:) - cflx_uplw_clr(:,:)) + &
								 (dflx_dnlw_clr(:,:) - dflx_uplw_clr(:,:)) - &
								 (bflx_dnlw_clr(:,:) - bflx_uplw_clr(:,:))) 
       dR_str_sraf(1:ldc%nproma,1:klev+1,krow)  = dR_str_sraf(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnsw_clr(:,:) - aflx_upsw_clr(:,:)) - & 
								 (cflx_dnsw_clr(:,:) - cflx_upsw_clr(:,:)) + &
								 (dflx_dnsw_clr(:,:) - dflx_upsw_clr(:,:)) - &
								 (bflx_dnsw_clr(:,:) - bflx_upsw_clr(:,:)))

     IF (prp3doutput) THEN
       dQ_str_trad(1:ldc%nproma,1:klev,krow)  = dQ_str_trad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnlw(:,1:klev) - aflx_uplw(:,1:klev)) - &
                                                   (aflx_dnlw(:,2:klev+1) - aflx_uplw(:,2:klev+1))) - &
                                                  ((cflx_dnlw(:,1:klev) - cflx_uplw(:,1:klev)) - & 
                                                   (cflx_dnlw(:,2:klev+1) - cflx_uplw(:,2:klev+1)))) - &
                                                 (((bflx_dnlw(:,1:klev) - bflx_uplw(:,1:klev)) - & 
                                                   (bflx_dnlw(:,2:klev+1) - bflx_uplw(:,2:klev+1))) - &
                                                  ((dflx_dnlw(:,1:klev) - dflx_uplw(:,1:klev)) - &
                                                   (dflx_dnlw(:,2:klev+1) - dflx_uplw(:,2:klev+1)))))
       dQ_str_srad(1:ldc%nproma,1:klev,krow)  = dQ_str_srad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnsw(:,1:klev) - aflx_upsw(:,1:klev)) - &
                                                   (aflx_dnsw(:,2:klev+1) - aflx_upsw(:,2:klev+1))) - &
                                                  ((cflx_dnsw(:,1:klev) - cflx_upsw(:,1:klev)) - & 
                                                   (cflx_dnsw(:,2:klev+1) - cflx_upsw(:,2:klev+1)))) - &
                                                 (((bflx_dnsw(:,1:klev) - bflx_upsw(:,1:klev)) - & 
                                                   (bflx_dnsw(:,2:klev+1) - bflx_upsw(:,2:klev+1))) - &
                                                  ((dflx_dnsw(:,1:klev) - dflx_upsw(:,1:klev)) - &
                                                   (dflx_dnsw(:,2:klev+1) - dflx_upsw(:,2:klev+1)))))
       dQ_str_traf(1:ldc%nproma,1:klev,krow)  = dQ_str_traf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnlw_clr(:,1:klev) - aflx_uplw_clr(:,1:klev)) - &
                                                   (aflx_dnlw_clr(:,2:klev+1) - aflx_uplw_clr(:,2:klev+1))) - &
                                                  ((cflx_dnlw_clr(:,1:klev) - cflx_uplw_clr(:,1:klev)) - & 
                                                   (cflx_dnlw_clr(:,2:klev+1) - cflx_uplw_clr(:,2:klev+1)))) - &
                                                 (((bflx_dnlw_clr(:,1:klev) - bflx_uplw_clr(:,1:klev)) - & 
                                                   (bflx_dnlw_clr(:,2:klev+1) - bflx_uplw_clr(:,2:klev+1))) - &
                                                  ((dflx_dnlw_clr(:,1:klev) - dflx_uplw_clr(:,1:klev)) - &
                                                   (dflx_dnlw_clr(:,2:klev+1) - dflx_uplw_clr(:,2:klev+1)))))
       dQ_str_sraf(1:ldc%nproma,1:klev,krow)  = dQ_str_sraf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnsw_clr(:,1:klev) - aflx_upsw_clr(:,1:klev)) - &
                                                   (aflx_dnsw_clr(:,2:klev+1) - aflx_upsw_clr(:,2:klev+1))) - &
                                                  ((cflx_dnsw_clr(:,1:klev) - cflx_upsw_clr(:,1:klev)) - & 
                                                   (cflx_dnsw_clr(:,2:klev+1) - cflx_upsw_clr(:,2:klev+1)))) - &
                                                 (((bflx_dnsw_clr(:,1:klev) - bflx_upsw_clr(:,1:klev)) - & 
                                                   (bflx_dnsw_clr(:,2:klev+1) - bflx_upsw_clr(:,2:klev+1))) - &
                                                  ((dflx_dnsw_clr(:,1:klev) - dflx_upsw_clr(:,1:klev)) - &
                                                   (dflx_dnsw_clr(:,2:klev+1) - dflx_upsw_clr(:,2:klev+1)))))
     END IF

       ! Calculate water vapor feedback:
        CALL message('','water vapor feedback')
       CALL psrad_interface( &
         & iaero           ,kproma          ,kbdim            ,klev             ,& 
         & krow            ,ktrac           ,ktype            ,nb_sw            ,& 
         & loland          ,loglac          ,cemiss           ,cos_mu0          ,& 
         & pgeom1          ,alb_vis_dir     ,alb_nir_dir      ,alb_vis_dif      ,&
         & alb_nir_dif     ,pp_fl           ,pp_hl            ,pp_sfc           ,&
         & tk_fl           ,tk_hl           ,tk_sfc           ,zq_vap           ,&
         & xq_liq          ,xq_ice          ,cdnc             ,cld_frc          ,&
         & cld_cvr         ,xm_o3           ,xm_co2           ,xm_ch4           ,&
         & xm_n2o          ,xm_cfc          ,xm_o2            ,pxtm1            ,&
         & cflx_uplw   ,cflx_uplw_clr,cflx_dnlw       ,cflx_dnlw_clr    ,&
         & cflx_upsw   ,cflx_upsw_clr,cflx_dnsw       ,cflx_dnsw_clr    ,&
         & vis_dff_sfc     ,dpar_sfc         ,nir_dff_sfc     ,dvis_sfc         ,&
         & par_dff_sfc                                                          )


       CALL psrad_interface( &
         & ziaero          ,kproma          ,kbdim            ,klev             ,& 
         & krow            ,ktrac           ,zktype           ,nb_sw            ,&
         & zloland         ,zloglac         ,zcemiss          ,zcos_mu0         ,& 
         & zpgeom1         ,zalb_vis_dir    ,zalb_nir_dir     ,zalb_vis_dif     ,&
         & zalb_nir_dif    ,zpp_fl          ,zpp_hl           ,zpp_sfc          ,& 
         & ztk_fl          ,ztk_hl          ,ztk_sfc          ,xq_vap           ,&
         & zq_liq          ,zq_ice          ,zcdnc            ,zcld_frc         ,&
         & zcld_cvr        ,zm_o3           ,zm_co2           ,zm_ch4           ,&
         & zm_n2o          ,zm_cfc          ,zm_o2            ,pxtm1            ,&
         & dflx_uplw   ,dflx_uplw_clr,dflx_dnlw       ,dflx_dnlw_clr    ,&
         & dflx_upsw   ,dflx_upsw_clr,dflx_dnsw       ,dflx_dnsw_clr    ,&
         & vis_dff_sfc     ,dpar_sfc         ,nir_dff_sfc     ,dvis_sfc         ,&
         & par_dff_sfc                                                          )

       dR_vap_trad0(1:ldc%nproma,krow)  = dR_vap_trad0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw(:,1) - aflx_uplw(:,1)) - & 
							 (cflx_dnlw(:,1) - cflx_uplw(:,1)) + &
                                                         (dflx_dnlw(:,1) - dflx_uplw(:,1)) - &
                                                         (bflx_dnlw(:,1) - bflx_uplw(:,1))) 
       dR_vap_srad0(1:ldc%nproma,krow)  = dR_vap_srad0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw(:,1) - aflx_upsw(:,1)) - & 
							 (cflx_dnsw(:,1) - cflx_upsw(:,1)) + &
                                                         (dflx_dnsw(:,1) - dflx_upsw(:,1)) - &
                                                         (bflx_dnsw(:,1) - bflx_upsw(:,1))) 
       dR_vap_trads(1:ldc%nproma,krow)  = dR_vap_trads(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw(:,klev+1) - aflx_uplw(:,klev+1)) - & 
							 (cflx_dnlw(:,klev+1) - cflx_uplw(:,klev+1)) + &
                                                         (dflx_dnlw(:,klev+1) - dflx_uplw(:,klev+1)) - &
                                                         (bflx_dnlw(:,klev+1) - bflx_uplw(:,klev+1))) 
       dR_vap_srads(1:ldc%nproma,krow)  = dR_vap_srads(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw(:,klev+1) - aflx_upsw(:,klev+1)) - & 
							 (cflx_dnsw(:,klev+1) - cflx_upsw(:,klev+1)) + &
                                                         (dflx_dnsw(:,klev+1) - dflx_upsw(:,klev+1)) - &
                                                         (bflx_dnsw(:,klev+1) - bflx_upsw(:,klev+1))) 
       dR_vap_traf0(1:ldc%nproma,krow)  = dR_vap_traf0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw_clr(:,1) - aflx_uplw_clr(:,1)) - & 
							 (cflx_dnlw_clr(:,1) - cflx_uplw_clr(:,1)) + &
                                                         (dflx_dnlw_clr(:,1) - dflx_uplw_clr(:,1)) - &
                                                         (bflx_dnlw_clr(:,1) - bflx_uplw_clr(:,1))) 
       dR_vap_sraf0(1:ldc%nproma,krow)  = dR_vap_sraf0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw_clr(:,1) - aflx_upsw_clr(:,1)) - & 
							 (cflx_dnsw_clr(:,1) - cflx_upsw_clr(:,1)) + &
                                                         (dflx_dnsw_clr(:,1) - dflx_upsw_clr(:,1)) - &
                                                         (bflx_dnsw_clr(:,1) - bflx_upsw_clr(:,1))) 
       dR_vap_trafs(1:ldc%nproma,krow)  = dR_vap_trafs(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw_clr(:,klev+1) - aflx_uplw_clr(:,klev+1)) - & 
							 (cflx_dnlw_clr(:,klev+1) - cflx_uplw_clr(:,klev+1)) + &
                                                         (dflx_dnlw_clr(:,klev+1) - dflx_uplw_clr(:,klev+1)) - &
                                                         (bflx_dnlw_clr(:,klev+1) - bflx_uplw_clr(:,klev+1))) 
       dR_vap_srafs(1:ldc%nproma,krow)  = dR_vap_srafs(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw_clr(:,klev+1) - aflx_upsw_clr(:,klev+1)) - & 
							 (cflx_dnsw_clr(:,klev+1) - cflx_upsw_clr(:,klev+1)) + &
                                                         (dflx_dnsw_clr(:,klev+1) - dflx_upsw_clr(:,klev+1)) - &
                                                         (bflx_dnsw_clr(:,klev+1) - bflx_upsw_clr(:,klev+1))) 

       dR_vap_trad(1:ldc%nproma,1:klev+1,krow)  = dR_vap_trad(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnlw(:,:) - aflx_uplw(:,:)) - & 
								 (cflx_dnlw(:,:) - cflx_uplw(:,:)) + &
								 (dflx_dnlw(:,:) - dflx_uplw(:,:)) - &
								 (bflx_dnlw(:,:) - bflx_uplw(:,:))) 
       dR_vap_srad(1:ldc%nproma,1:klev+1,krow)  = dR_vap_srad(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnsw(:,:) - aflx_upsw(:,:)) - & 
								 (cflx_dnsw(:,:) - cflx_upsw(:,:)) + &
								 (dflx_dnsw(:,:) - dflx_upsw(:,:)) - &
								 (bflx_dnsw(:,:) - bflx_upsw(:,:))) 
       dR_vap_traf(1:ldc%nproma,1:klev+1,krow)  = dR_vap_traf(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnlw_clr(:,:) - aflx_uplw_clr(:,:)) - & 
								 (cflx_dnlw_clr(:,:) - cflx_uplw_clr(:,:)) + &
								 (dflx_dnlw_clr(:,:) - dflx_uplw_clr(:,:)) - &
								 (bflx_dnlw_clr(:,:) - bflx_uplw_clr(:,:))) 
       dR_vap_sraf(1:ldc%nproma,1:klev+1,krow)  = dR_vap_sraf(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnsw_clr(:,:) - aflx_upsw_clr(:,:)) - & 
								 (cflx_dnsw_clr(:,:) - cflx_upsw_clr(:,:)) + &
								 (dflx_dnsw_clr(:,:) - dflx_upsw_clr(:,:)) - &
								 (bflx_dnsw_clr(:,:) - bflx_upsw_clr(:,:)))

     IF (prp3doutput) THEN
       dQ_vap_trad(1:ldc%nproma,1:klev,krow)  = dQ_vap_trad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnlw(:,1:klev) - aflx_uplw(:,1:klev)) - &
                                                   (aflx_dnlw(:,2:klev+1) - aflx_uplw(:,2:klev+1))) - &
                                                  ((cflx_dnlw(:,1:klev) - cflx_uplw(:,1:klev)) - & 
                                                   (cflx_dnlw(:,2:klev+1) - cflx_uplw(:,2:klev+1)))) - &
                                                 (((bflx_dnlw(:,1:klev) - bflx_uplw(:,1:klev)) - & 
                                                   (bflx_dnlw(:,2:klev+1) - bflx_uplw(:,2:klev+1))) - &
                                                  ((dflx_dnlw(:,1:klev) - dflx_uplw(:,1:klev)) - &
                                                   (dflx_dnlw(:,2:klev+1) - dflx_uplw(:,2:klev+1)))))
       dQ_vap_srad(1:ldc%nproma,1:klev,krow)  = dQ_vap_srad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnsw(:,1:klev) - aflx_upsw(:,1:klev)) - &
                                                   (aflx_dnsw(:,2:klev+1) - aflx_upsw(:,2:klev+1))) - &
                                                  ((cflx_dnsw(:,1:klev) - cflx_upsw(:,1:klev)) - & 
                                                   (cflx_dnsw(:,2:klev+1) - cflx_upsw(:,2:klev+1)))) - &
                                                 (((bflx_dnsw(:,1:klev) - bflx_upsw(:,1:klev)) - & 
                                                   (bflx_dnsw(:,2:klev+1) - bflx_upsw(:,2:klev+1))) - &
                                                  ((dflx_dnsw(:,1:klev) - dflx_upsw(:,1:klev)) - &
                                                   (dflx_dnsw(:,2:klev+1) - dflx_upsw(:,2:klev+1)))))
       dQ_vap_traf(1:ldc%nproma,1:klev,krow)  = dQ_vap_traf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnlw_clr(:,1:klev) - aflx_uplw_clr(:,1:klev)) - &
                                                   (aflx_dnlw_clr(:,2:klev+1) - aflx_uplw_clr(:,2:klev+1))) - &
                                                  ((cflx_dnlw_clr(:,1:klev) - cflx_uplw_clr(:,1:klev)) - & 
                                                   (cflx_dnlw_clr(:,2:klev+1) - cflx_uplw_clr(:,2:klev+1)))) - &
                                                 (((bflx_dnlw_clr(:,1:klev) - bflx_uplw_clr(:,1:klev)) - & 
                                                   (bflx_dnlw_clr(:,2:klev+1) - bflx_uplw_clr(:,2:klev+1))) - &
                                                  ((dflx_dnlw_clr(:,1:klev) - dflx_uplw_clr(:,1:klev)) - &
                                                   (dflx_dnlw_clr(:,2:klev+1) - dflx_uplw_clr(:,2:klev+1)))))
       dQ_vap_sraf(1:ldc%nproma,1:klev,krow)  = dQ_vap_sraf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnsw_clr(:,1:klev) - aflx_upsw_clr(:,1:klev)) - &
                                                   (aflx_dnsw_clr(:,2:klev+1) - aflx_upsw_clr(:,2:klev+1))) - &
                                                  ((cflx_dnsw_clr(:,1:klev) - cflx_upsw_clr(:,1:klev)) - & 
                                                   (cflx_dnsw_clr(:,2:klev+1) - cflx_upsw_clr(:,2:klev+1)))) - &
                                                 (((bflx_dnsw_clr(:,1:klev) - bflx_upsw_clr(:,1:klev)) - & 
                                                   (bflx_dnsw_clr(:,2:klev+1) - bflx_upsw_clr(:,2:klev+1))) - &
                                                  ((dflx_dnsw_clr(:,1:klev) - dflx_upsw_clr(:,1:klev)) - &
                                                   (dflx_dnsw_clr(:,2:klev+1) - dflx_upsw_clr(:,2:klev+1)))))
     END IF


       ! Calculate total cloud feedback (with Twomey):

     SELECT CASE (iaero(1))
     CASE DEFAULT
      iaero_indirect(1) = iaero_ref            ! cloud droplet number in state d calculated from simple-plumes if they are present
      iaero_direct(1) = iaero(1)
     CASE (0)      
      iaero_direct(1) = 0 
      iaero_indirect(1) = 0 
     CASE (8)                              
       iaero_direct(1) = 12 
       iaero_indirect(1) = 11   
     CASE (9)                              ! Depending whether or not volcanoes in current state
       iaero_direct(1) = 10                   ! the direct effect in state c changes
       iaero_indirect(1) = 11
     END SELECT                     
     CALL message('','total cloud feedback')
      CALL psrad_interface( &
        & iaero_direct    ,kproma          ,kbdim            ,klev             ,& 
        & krow            ,ktrac           ,ktype            ,nb_sw            ,& 
        & loland          ,loglac          ,cemiss           ,cos_mu0          ,& 
        & pgeom1          ,alb_vis_dir     ,alb_nir_dir      ,alb_vis_dif      ,&
        & alb_nir_dif     ,pp_fl           ,pp_hl            ,pp_sfc           ,&
        & tk_fl           ,tk_hl           ,tk_sfc           ,xq_vap           ,&
        & zq_liq          ,zq_ice          ,zcdnc            ,zcld_frc         ,&
        & zcld_cvr        ,xm_o3           ,xm_co2           ,xm_ch4           ,&
        & xm_n2o          ,xm_cfc          ,xm_o2            ,pxtm1            ,&
        & cflx_uplw       ,cflx_uplw_clr,cflx_dnlw           ,cflx_dnlw_clr    ,&
        & cflx_upsw       ,cflx_upsw_clr,cflx_dnsw           ,cflx_dnsw_clr    ,&
        & vis_dff_sfc     ,dpar_sfc         ,nir_dff_sfc     ,dvis_sfc         ,&
        & par_dff_sfc                                                          )

      CALL psrad_interface( &
        & iaero_indirect  ,kproma          ,kbdim            ,klev             ,& 
        & krow            ,ktrac           ,zktype           ,nb_sw            ,&
        & zloland         ,zloglac         ,zcemiss          ,zcos_mu0         ,& 
        & zpgeom1         ,zalb_vis_dir    ,zalb_nir_dir     ,zalb_vis_dif     ,&
        & zalb_nir_dif    ,zpp_fl          ,zpp_hl           ,zpp_sfc          ,& 
        & ztk_fl          ,ztk_hl          ,ztk_sfc          ,zq_vap           ,&
        & xq_liq          ,xq_ice          ,cdnc             ,cld_frc          ,&
        & cld_cvr         ,zm_o3           ,zm_co2           ,zm_ch4           ,&
        & zm_n2o          ,zm_cfc          ,zm_o2            ,pxtm1            ,&
        & dflx_uplw       ,dflx_uplw_clr   ,dflx_dnlw        ,dflx_dnlw_clr    ,&
        & dflx_upsw       ,dflx_upsw_clr   ,dflx_dnsw        ,dflx_dnsw_clr    ,&
        & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
        & par_dff_sfc                                                          )

      dR_cld_trad0(1:ldc%nproma,krow)  = dR_cld_trad0(1:ldc%nproma,krow) + &
                                         dtprps*0.5_wp*((aflx_dnlw(:,1) - aflx_uplw(:,1)) - & 
              (cflx_dnlw(:,1) - cflx_uplw(:,1)) + &
                                                        (dflx_dnlw(:,1) - dflx_uplw(:,1)) - &
                                                        (bflx_dnlw(:,1) - bflx_uplw(:,1))) 
      dR_cld_srad0(1:ldc%nproma,krow)  = dR_cld_srad0(1:ldc%nproma,krow) + &
                                         dtprps*0.5_wp*((aflx_dnsw(:,1) - aflx_upsw(:,1)) - & 
              (cflx_dnsw(:,1) - cflx_upsw(:,1)) + &
                                                        (dflx_dnsw(:,1) - dflx_upsw(:,1)) - &
                                                        (bflx_dnsw(:,1) - bflx_upsw(:,1))) 
      dR_cld_trads(1:ldc%nproma,krow)  = dR_cld_trads(1:ldc%nproma,krow) + &
                                         dtprps*0.5_wp*((aflx_dnlw(:,klev+1) - aflx_uplw(:,klev+1)) - & 
              (cflx_dnlw(:,klev+1) - cflx_uplw(:,klev+1)) + &
                                                        (dflx_dnlw(:,klev+1) - dflx_uplw(:,klev+1)) - &
                                                        (bflx_dnlw(:,klev+1) - bflx_uplw(:,klev+1))) 
      dR_cld_srads(1:ldc%nproma,krow)  = dR_cld_srads(1:ldc%nproma,krow) + &
                                         dtprps*0.5_wp*((aflx_dnsw(:,klev+1) - aflx_upsw(:,klev+1)) - & 
              (cflx_dnsw(:,klev+1) - cflx_upsw(:,klev+1)) + &
                                                        (dflx_dnsw(:,klev+1) - dflx_upsw(:,klev+1)) - &
                                                        (bflx_dnsw(:,klev+1) - bflx_upsw(:,klev+1))) 
      dR_cld_traf0(1:ldc%nproma,krow)  = dR_cld_traf0(1:ldc%nproma,krow) + &
                                         dtprps*0.5_wp*((aflx_dnlw_clr(:,1) - aflx_uplw_clr(:,1)) - & 
              (cflx_dnlw_clr(:,1) - cflx_uplw_clr(:,1)) + &
                                                        (dflx_dnlw_clr(:,1) - dflx_uplw_clr(:,1)) - &
                                                        (bflx_dnlw_clr(:,1) - bflx_uplw_clr(:,1))) 
      dR_cld_sraf0(1:ldc%nproma,krow)  = dR_cld_sraf0(1:ldc%nproma,krow) + &
                                         dtprps*0.5_wp*((aflx_dnsw_clr(:,1) - aflx_upsw_clr(:,1)) - & 
              (cflx_dnsw_clr(:,1) - cflx_upsw_clr(:,1)) + &
                                                        (dflx_dnsw_clr(:,1) - dflx_upsw_clr(:,1)) - &
                                                        (bflx_dnsw_clr(:,1) - bflx_upsw_clr(:,1))) 
      dR_cld_trafs(1:ldc%nproma,krow)  = dR_cld_trafs(1:ldc%nproma,krow) + &
                                         dtprps*0.5_wp*((aflx_dnlw_clr(:,klev+1) - aflx_uplw_clr(:,klev+1)) - & 
              (cflx_dnlw_clr(:,klev+1) - cflx_uplw_clr(:,klev+1)) + &
                                                        (dflx_dnlw_clr(:,klev+1) - dflx_uplw_clr(:,klev+1)) - &
                                                        (bflx_dnlw_clr(:,klev+1) - bflx_uplw_clr(:,klev+1))) 
      dR_cld_srafs(1:ldc%nproma,krow)  = dR_cld_srafs(1:ldc%nproma,krow) + &
                                         dtprps*0.5_wp*((aflx_dnsw_clr(:,klev+1) - aflx_upsw_clr(:,klev+1)) - & 
              (cflx_dnsw_clr(:,klev+1) - cflx_upsw_clr(:,klev+1)) + &
                                                        (dflx_dnsw_clr(:,klev+1) - dflx_upsw_clr(:,klev+1)) - &
                                                        (bflx_dnsw_clr(:,klev+1) - bflx_upsw_clr(:,klev+1))) 

      dR_cld_trad(1:ldc%nproma,1:klev+1,krow)  = dR_cld_trad(1:ldc%nproma,1:klev+1,krow) + &
                                                 dtprps*0.5_wp*((aflx_dnlw(:,:) - aflx_uplw(:,:)) - & 
                (cflx_dnlw(:,:) - cflx_uplw(:,:)) + &
                (dflx_dnlw(:,:) - dflx_uplw(:,:)) - &
                (bflx_dnlw(:,:) - bflx_uplw(:,:))) 
      dR_cld_srad(1:ldc%nproma,1:klev+1,krow)  = dR_cld_srad(1:ldc%nproma,1:klev+1,krow) + &
                                                 dtprps*0.5_wp*((aflx_dnsw(:,:) - aflx_upsw(:,:)) - & 
                (cflx_dnsw(:,:) - cflx_upsw(:,:)) + &
                (dflx_dnsw(:,:) - dflx_upsw(:,:)) - &
                (bflx_dnsw(:,:) - bflx_upsw(:,:))) 
      dR_cld_traf(1:ldc%nproma,1:klev+1,krow)  = dR_cld_traf(1:ldc%nproma,1:klev+1,krow) + &
                                                 dtprps*0.5_wp*((aflx_dnlw_clr(:,:) - aflx_uplw_clr(:,:)) - & 
                (cflx_dnlw_clr(:,:) - cflx_uplw_clr(:,:)) + &
                (dflx_dnlw_clr(:,:) - dflx_uplw_clr(:,:)) - &
                (bflx_dnlw_clr(:,:) - bflx_uplw_clr(:,:))) 
      dR_cld_sraf(1:ldc%nproma,1:klev+1,krow)  = dR_cld_sraf(1:ldc%nproma,1:klev+1,krow) + &
                                                 dtprps*0.5_wp*((aflx_dnsw_clr(:,:) - aflx_upsw_clr(:,:)) - & 
                (cflx_dnsw_clr(:,:) - cflx_upsw_clr(:,:)) + &
                (dflx_dnsw_clr(:,:) - dflx_upsw_clr(:,:)) - &
                (bflx_dnsw_clr(:,:) - bflx_upsw_clr(:,:)))

    IF (prp3doutput) THEN
      dQ_cld_trad(1:ldc%nproma,1:klev,krow)  = dQ_cld_trad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                               ((((aflx_dnlw(:,1:klev) - aflx_uplw(:,1:klev)) - &
                                                  (aflx_dnlw(:,2:klev+1) - aflx_uplw(:,2:klev+1))) - &
                                                 ((cflx_dnlw(:,1:klev) - cflx_uplw(:,1:klev)) - & 
                                                  (cflx_dnlw(:,2:klev+1) - cflx_uplw(:,2:klev+1)))) - &
                                                (((bflx_dnlw(:,1:klev) - bflx_uplw(:,1:klev)) - & 
                                                  (bflx_dnlw(:,2:klev+1) - bflx_uplw(:,2:klev+1))) - &
                                                 ((dflx_dnlw(:,1:klev) - dflx_uplw(:,1:klev)) - &
                                                  (dflx_dnlw(:,2:klev+1) - dflx_uplw(:,2:klev+1)))))
      dQ_cld_srad(1:ldc%nproma,1:klev,krow)  = dQ_cld_srad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                               ((((aflx_dnsw(:,1:klev) - aflx_upsw(:,1:klev)) - &
                                                  (aflx_dnsw(:,2:klev+1) - aflx_upsw(:,2:klev+1))) - &
                                                 ((cflx_dnsw(:,1:klev) - cflx_upsw(:,1:klev)) - & 
                                                  (cflx_dnsw(:,2:klev+1) - cflx_upsw(:,2:klev+1)))) - &
                                                (((bflx_dnsw(:,1:klev) - bflx_upsw(:,1:klev)) - & 
                                                  (bflx_dnsw(:,2:klev+1) - bflx_upsw(:,2:klev+1))) - &
                                                 ((dflx_dnsw(:,1:klev) - dflx_upsw(:,1:klev)) - &
                                                  (dflx_dnsw(:,2:klev+1) - dflx_upsw(:,2:klev+1)))))
      dQ_cld_traf(1:ldc%nproma,1:klev,krow)  = dQ_cld_traf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                               ((((aflx_dnlw_clr(:,1:klev) - aflx_uplw_clr(:,1:klev)) - &
                                                  (aflx_dnlw_clr(:,2:klev+1) - aflx_uplw_clr(:,2:klev+1))) - &
                                                 ((cflx_dnlw_clr(:,1:klev) - cflx_uplw_clr(:,1:klev)) - & 
                                                  (cflx_dnlw_clr(:,2:klev+1) - cflx_uplw_clr(:,2:klev+1)))) - &
                                                (((bflx_dnlw_clr(:,1:klev) - bflx_uplw_clr(:,1:klev)) - & 
                                                  (bflx_dnlw_clr(:,2:klev+1) - bflx_uplw_clr(:,2:klev+1))) - &
                                                 ((dflx_dnlw_clr(:,1:klev) - dflx_uplw_clr(:,1:klev)) - &
                                                  (dflx_dnlw_clr(:,2:klev+1) - dflx_uplw_clr(:,2:klev+1)))))
      dQ_cld_sraf(1:ldc%nproma,1:klev,krow)  = dQ_cld_sraf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                               ((((aflx_dnsw_clr(:,1:klev) - aflx_upsw_clr(:,1:klev)) - &
                                                  (aflx_dnsw_clr(:,2:klev+1) - aflx_upsw_clr(:,2:klev+1))) - &
                                                 ((cflx_dnsw_clr(:,1:klev) - cflx_upsw_clr(:,1:klev)) - & 
                                                  (cflx_dnsw_clr(:,2:klev+1) - cflx_upsw_clr(:,2:klev+1)))) - &
                                                (((bflx_dnsw_clr(:,1:klev) - bflx_upsw_clr(:,1:klev)) - & 
                                                  (bflx_dnsw_clr(:,2:klev+1) - bflx_upsw_clr(:,2:klev+1))) - &
                                                 ((dflx_dnsw_clr(:,1:klev) - dflx_upsw_clr(:,1:klev)) - &
                                                  (dflx_dnsw_clr(:,2:klev+1) - dflx_upsw_clr(:,2:klev+1)))))
    END IF

          ! Calculate surface albedo feedback:
    CALL psrad_interface( &
          & iaero           ,kproma          ,kbdim            ,klev             ,& 
          & krow            ,ktrac           ,ktype            ,nb_sw            ,&
          & loland          ,loglac          ,cemiss           ,cos_mu0          ,& 
          & pgeom1          ,zalb_vis_dir    ,zalb_nir_dir     ,zalb_vis_dif     ,&
          & zalb_nir_dif    ,pp_fl           ,pp_hl            ,pp_sfc           ,& 
          & tk_fl           ,tk_hl           ,tk_sfc           ,xq_vap           ,&
          & xq_liq          ,xq_ice          ,cdnc             ,cld_frc          ,&
          & cld_cvr         ,xm_o3           ,xm_co2           ,xm_ch4           ,&
          & xm_n2o          ,xm_cfc          ,xm_o2            ,pxtm1            ,&
          & cflx_uplw       ,cflx_uplw_clr   ,cflx_dnlw        ,cflx_dnlw_clr    ,&
          & cflx_upsw       ,cflx_upsw_clr   ,cflx_dnsw        ,cflx_dnsw_clr    ,&
          & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
          & par_dff_sfc                                                          )

    CALL psrad_interface( &
          & ziaero          ,kproma          ,kbdim            ,klev             ,& 
          & krow            ,ktrac           ,zktype           ,nb_sw            ,& 
          & zloland         ,zloglac         ,zcemiss          ,zcos_mu0         ,& 
          & zpgeom1         ,alb_vis_dir     ,alb_nir_dir      ,alb_vis_dif      ,&
          & alb_nir_dif     ,zpp_fl          ,zpp_hl           ,zpp_sfc          ,&
          & ztk_fl          ,ztk_hl          ,ztk_sfc          ,zq_vap           ,&
          & zq_liq          ,zq_ice          ,zcdnc            ,zcld_frc         ,&
          & zcld_cvr        ,zm_o3           ,zm_co2           ,zm_ch4           ,&
          & zm_n2o          ,zm_cfc          ,zm_o2            ,pxtm1            ,&
          & dflx_uplw       ,dflx_uplw_clr   ,dflx_dnlw        ,dflx_dnlw_clr    ,&
          & dflx_upsw       ,dflx_upsw_clr   ,dflx_dnsw        ,dflx_dnsw_clr    ,&
          & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
          & par_dff_sfc                                                          )

        dR_alb_trad0(1:ldc%nproma,krow)  = dR_alb_trad0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw(:,1) - aflx_uplw(:,1)) - & 
                (cflx_dnlw(:,1) - cflx_uplw(:,1)) + &
                                                          (dflx_dnlw(:,1) - dflx_uplw(:,1)) - &
                                                          (bflx_dnlw(:,1) - bflx_uplw(:,1))) 
        dR_alb_srad0(1:ldc%nproma,krow)  = dR_alb_srad0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw(:,1) - aflx_upsw(:,1)) - & 
                (cflx_dnsw(:,1) - cflx_upsw(:,1)) + &
                                                          (dflx_dnsw(:,1) - dflx_upsw(:,1)) - &
                                                          (bflx_dnsw(:,1) - bflx_upsw(:,1))) 
        dR_alb_trads(1:ldc%nproma,krow)  = dR_alb_trads(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw(:,klev+1) - aflx_uplw(:,klev+1)) - & 
                (cflx_dnlw(:,klev+1) - cflx_uplw(:,klev+1)) + &
                                                          (dflx_dnlw(:,klev+1) - dflx_uplw(:,klev+1)) - &
                                                          (bflx_dnlw(:,klev+1) - bflx_uplw(:,klev+1))) 
        dR_alb_srads(1:ldc%nproma,krow)  = dR_alb_srads(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw(:,klev+1) - aflx_upsw(:,klev+1)) - & 
                (cflx_dnsw(:,klev+1) - cflx_upsw(:,klev+1)) + &
                                                          (dflx_dnsw(:,klev+1) - dflx_upsw(:,klev+1)) - &
                                                          (bflx_dnsw(:,klev+1) - bflx_upsw(:,klev+1))) 
        dR_alb_traf0(1:ldc%nproma,krow)  = dR_alb_traf0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw_clr(:,1) - aflx_uplw_clr(:,1)) - & 
                (cflx_dnlw_clr(:,1) - cflx_uplw_clr(:,1)) + &
                                                          (dflx_dnlw_clr(:,1) - dflx_uplw_clr(:,1)) - &
                                                          (bflx_dnlw_clr(:,1) - bflx_uplw_clr(:,1))) 
        dR_alb_sraf0(1:ldc%nproma,krow)  = dR_alb_sraf0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw_clr(:,1) - aflx_upsw_clr(:,1)) - & 
                (cflx_dnsw_clr(:,1) - cflx_upsw_clr(:,1)) + &
                                                          (dflx_dnsw_clr(:,1) - dflx_upsw_clr(:,1)) - &
                                                          (bflx_dnsw_clr(:,1) - bflx_upsw_clr(:,1))) 
        dR_alb_trafs(1:ldc%nproma,krow)  = dR_alb_trafs(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw_clr(:,klev+1) - aflx_uplw_clr(:,klev+1)) - & 
                (cflx_dnlw_clr(:,klev+1) - cflx_uplw_clr(:,klev+1)) + &
                                                          (dflx_dnlw_clr(:,klev+1) - dflx_uplw_clr(:,klev+1)) - &
                                                          (bflx_dnlw_clr(:,klev+1) - bflx_uplw_clr(:,klev+1))) 
        dR_alb_srafs(1:ldc%nproma,krow)  = dR_alb_srafs(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw_clr(:,klev+1) - aflx_upsw_clr(:,klev+1)) - & 
                (cflx_dnsw_clr(:,klev+1) - cflx_upsw_clr(:,klev+1)) + &
                                                          (dflx_dnsw_clr(:,klev+1) - dflx_upsw_clr(:,klev+1)) - &
                                                          (bflx_dnsw_clr(:,klev+1) - bflx_upsw_clr(:,klev+1))) 

        dR_alb_trad(1:ldc%nproma,1:klev+1,krow)  = dR_alb_trad(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnlw(:,:) - aflx_uplw(:,:)) - & 
                  (cflx_dnlw(:,:) - cflx_uplw(:,:)) + &
                  (dflx_dnlw(:,:) - dflx_uplw(:,:)) - &
                  (bflx_dnlw(:,:) - bflx_uplw(:,:))) 
        dR_alb_srad(1:ldc%nproma,1:klev+1,krow)  = dR_alb_srad(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnsw(:,:) - aflx_upsw(:,:)) - & 
                  (cflx_dnsw(:,:) - cflx_upsw(:,:)) + &
                  (dflx_dnsw(:,:) - dflx_upsw(:,:)) - &
                  (bflx_dnsw(:,:) - bflx_upsw(:,:))) 
        dR_alb_traf(1:ldc%nproma,1:klev+1,krow)  = dR_alb_traf(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnlw_clr(:,:) - aflx_uplw_clr(:,:)) - & 
                  (cflx_dnlw_clr(:,:) - cflx_uplw_clr(:,:)) + &
                  (dflx_dnlw_clr(:,:) - dflx_uplw_clr(:,:)) - &
                  (bflx_dnlw_clr(:,:) - bflx_uplw_clr(:,:))) 
        dR_alb_sraf(1:ldc%nproma,1:klev+1,krow)  = dR_alb_sraf(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnsw_clr(:,:) - aflx_upsw_clr(:,:)) - & 
                  (cflx_dnsw_clr(:,:) - cflx_upsw_clr(:,:)) + &
                  (dflx_dnsw_clr(:,:) - dflx_upsw_clr(:,:)) - &
                  (bflx_dnsw_clr(:,:) - bflx_upsw_clr(:,:)))

      IF (prp3doutput) THEN
        dQ_alb_trad(1:ldc%nproma,1:klev,krow)  = dQ_alb_trad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnlw(:,1:klev) - aflx_uplw(:,1:klev)) - &
                                                    (aflx_dnlw(:,2:klev+1) - aflx_uplw(:,2:klev+1))) - &
                                                  ((cflx_dnlw(:,1:klev) - cflx_uplw(:,1:klev)) - & 
                                                    (cflx_dnlw(:,2:klev+1) - cflx_uplw(:,2:klev+1)))) - &
                                                  (((bflx_dnlw(:,1:klev) - bflx_uplw(:,1:klev)) - & 
                                                    (bflx_dnlw(:,2:klev+1) - bflx_uplw(:,2:klev+1))) - &
                                                  ((dflx_dnlw(:,1:klev) - dflx_uplw(:,1:klev)) - &
                                                    (dflx_dnlw(:,2:klev+1) - dflx_uplw(:,2:klev+1)))))
        dQ_alb_srad(1:ldc%nproma,1:klev,krow)  = dQ_alb_srad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnsw(:,1:klev) - aflx_upsw(:,1:klev)) - &
                                                    (aflx_dnsw(:,2:klev+1) - aflx_upsw(:,2:klev+1))) - &
                                                  ((cflx_dnsw(:,1:klev) - cflx_upsw(:,1:klev)) - & 
                                                    (cflx_dnsw(:,2:klev+1) - cflx_upsw(:,2:klev+1)))) - &
                                                  (((bflx_dnsw(:,1:klev) - bflx_upsw(:,1:klev)) - & 
                                                    (bflx_dnsw(:,2:klev+1) - bflx_upsw(:,2:klev+1))) - &
                                                  ((dflx_dnsw(:,1:klev) - dflx_upsw(:,1:klev)) - &
                                                    (dflx_dnsw(:,2:klev+1) - dflx_upsw(:,2:klev+1)))))
        dQ_alb_traf(1:ldc%nproma,1:klev,krow)  = dQ_alb_traf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnlw_clr(:,1:klev) - aflx_uplw_clr(:,1:klev)) - &
                                                    (aflx_dnlw_clr(:,2:klev+1) - aflx_uplw_clr(:,2:klev+1))) - &
                                                  ((cflx_dnlw_clr(:,1:klev) - cflx_uplw_clr(:,1:klev)) - & 
                                                    (cflx_dnlw_clr(:,2:klev+1) - cflx_uplw_clr(:,2:klev+1)))) - &
                                                  (((bflx_dnlw_clr(:,1:klev) - bflx_uplw_clr(:,1:klev)) - & 
                                                    (bflx_dnlw_clr(:,2:klev+1) - bflx_uplw_clr(:,2:klev+1))) - &
                                                  ((dflx_dnlw_clr(:,1:klev) - dflx_uplw_clr(:,1:klev)) - &
                                                    (dflx_dnlw_clr(:,2:klev+1) - dflx_uplw_clr(:,2:klev+1)))))
        dQ_alb_sraf(1:ldc%nproma,1:klev,krow)  = dQ_alb_sraf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnsw_clr(:,1:klev) - aflx_upsw_clr(:,1:klev)) - &
                                                    (aflx_dnsw_clr(:,2:klev+1) - aflx_upsw_clr(:,2:klev+1))) - &
                                                  ((cflx_dnsw_clr(:,1:klev) - cflx_upsw_clr(:,1:klev)) - & 
                                                    (cflx_dnsw_clr(:,2:klev+1) - cflx_upsw_clr(:,2:klev+1)))) - &
                                                  (((bflx_dnsw_clr(:,1:klev) - bflx_upsw_clr(:,1:klev)) - & 
                                                    (bflx_dnsw_clr(:,2:klev+1) - bflx_upsw_clr(:,2:klev+1))) - &
                                                  ((dflx_dnsw_clr(:,1:klev) - dflx_upsw_clr(:,1:klev)) - &
                                                    (dflx_dnsw_clr(:,2:klev+1) - dflx_upsw_clr(:,2:klev+1)))))
      END IF


        ! Calculate direct CO2 forcing:

    CALL psrad_interface( &
          & iaero           ,kproma          ,kbdim            ,klev             ,& 
          & krow            ,ktrac           ,ktype            ,nb_sw            ,& 
          & loland          ,loglac          ,cemiss           ,cos_mu0          ,& 
          & pgeom1          ,alb_vis_dir     ,alb_nir_dir      ,alb_vis_dif      ,&
          & alb_nir_dif     ,pp_fl           ,pp_hl            ,pp_sfc           ,&
          & tk_fl           ,tk_hl           ,tk_sfc           ,xq_vap           ,&
          & xq_liq          ,xq_ice          ,cdnc             ,cld_frc          ,&
          & cld_cvr         ,xm_o3           ,zm_co2           ,xm_ch4           ,&
          & xm_n2o          ,xm_cfc          ,xm_o2            ,pxtm1            ,&
          & cflx_uplw       ,cflx_uplw_clr   ,cflx_dnlw        ,cflx_dnlw_clr    ,&
          & cflx_upsw       ,cflx_upsw_clr   ,cflx_dnsw        ,cflx_dnsw_clr    ,& 
          & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
          & par_dff_sfc                                                          )


    CALL psrad_interface( &
          & ziaero          ,kproma          ,kbdim            ,klev             ,& 
          & krow            ,ktrac           ,zktype           ,nb_sw            ,&
          & zloland         ,zloglac         ,zcemiss          ,zcos_mu0         ,& 
          & zpgeom1         ,zalb_vis_dir    ,zalb_nir_dir     ,zalb_vis_dif     ,&
          & zalb_nir_dif    ,zpp_fl          ,zpp_hl           ,zpp_sfc          ,& 
          & ztk_fl          ,ztk_hl          ,ztk_sfc          ,zq_vap           ,&
          & zq_liq          ,zq_ice          ,zcdnc            ,zcld_frc         ,&
          & zcld_cvr        ,zm_o3           ,xm_co2           ,zm_ch4           ,&
          & zm_n2o          ,zm_cfc          ,zm_o2            ,pxtm1            ,&
          & dflx_uplw       ,dflx_uplw_clr   ,dflx_dnlw        ,dflx_dnlw_clr    ,&
          & dflx_upsw       ,dflx_upsw_clr   ,dflx_dnsw        ,dflx_dnsw_clr    ,&
          & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
          & par_dff_sfc                                                          )

        dR_co2_trad0(1:ldc%nproma,krow)  = dR_co2_trad0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw(:,1) - aflx_uplw(:,1)) - & 
                (cflx_dnlw(:,1) - cflx_uplw(:,1)) + &
                                                          (dflx_dnlw(:,1) - dflx_uplw(:,1)) - &
                                                          (bflx_dnlw(:,1) - bflx_uplw(:,1))) 
        dR_co2_srad0(1:ldc%nproma,krow)  = dR_co2_srad0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw(:,1) - aflx_upsw(:,1)) - & 
                (cflx_dnsw(:,1) - cflx_upsw(:,1)) + &
                                                          (dflx_dnsw(:,1) - dflx_upsw(:,1)) - &
                                                          (bflx_dnsw(:,1) - bflx_upsw(:,1))) 
        dR_co2_trads(1:ldc%nproma,krow)  = dR_co2_trads(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw(:,klev+1) - aflx_uplw(:,klev+1)) - & 
                (cflx_dnlw(:,klev+1) - cflx_uplw(:,klev+1)) + &
                                                          (dflx_dnlw(:,klev+1) - dflx_uplw(:,klev+1)) - &
                                                          (bflx_dnlw(:,klev+1) - bflx_uplw(:,klev+1))) 
        dR_co2_srads(1:ldc%nproma,krow)  = dR_co2_srads(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw(:,klev+1) - aflx_upsw(:,klev+1)) - & 
                (cflx_dnsw(:,klev+1) - cflx_upsw(:,klev+1)) + &
                                                          (dflx_dnsw(:,klev+1) - dflx_upsw(:,klev+1)) - &
                                                          (bflx_dnsw(:,klev+1) - bflx_upsw(:,klev+1))) 
        dR_co2_traf0(1:ldc%nproma,krow)  = dR_co2_traf0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw_clr(:,1) - aflx_uplw_clr(:,1)) - & 
                (cflx_dnlw_clr(:,1) - cflx_uplw_clr(:,1)) + &
                                                          (dflx_dnlw_clr(:,1) - dflx_uplw_clr(:,1)) - &
                                                          (bflx_dnlw_clr(:,1) - bflx_uplw_clr(:,1))) 
        dR_co2_sraf0(1:ldc%nproma,krow)  = dR_co2_sraf0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw_clr(:,1) - aflx_upsw_clr(:,1)) - & 
                (cflx_dnsw_clr(:,1) - cflx_upsw_clr(:,1)) + &
                                                          (dflx_dnsw_clr(:,1) - dflx_upsw_clr(:,1)) - &
                                                          (bflx_dnsw_clr(:,1) - bflx_upsw_clr(:,1))) 
        dR_co2_trafs(1:ldc%nproma,krow)  = dR_co2_trafs(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw_clr(:,klev+1) - aflx_uplw_clr(:,klev+1)) - & 
                (cflx_dnlw_clr(:,klev+1) - cflx_uplw_clr(:,klev+1)) + &
                                                          (dflx_dnlw_clr(:,klev+1) - dflx_uplw_clr(:,klev+1)) - &
                                                          (bflx_dnlw_clr(:,klev+1) - bflx_uplw_clr(:,klev+1))) 
        dR_co2_srafs(1:ldc%nproma,krow)  = dR_co2_srafs(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw_clr(:,klev+1) - aflx_upsw_clr(:,klev+1)) - & 
                (cflx_dnsw_clr(:,klev+1) - cflx_upsw_clr(:,klev+1)) + &
                                                          (dflx_dnsw_clr(:,klev+1) - dflx_upsw_clr(:,klev+1)) - &
                                                          (bflx_dnsw_clr(:,klev+1) - bflx_upsw_clr(:,klev+1))) 

        dR_co2_trad(1:ldc%nproma,1:klev+1,krow)  = dR_co2_trad(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnlw(:,:) - aflx_uplw(:,:)) - & 
                  (cflx_dnlw(:,:) - cflx_uplw(:,:)) + &
                  (dflx_dnlw(:,:) - dflx_uplw(:,:)) - &
                  (bflx_dnlw(:,:) - bflx_uplw(:,:))) 
        dR_co2_srad(1:ldc%nproma,1:klev+1,krow)  = dR_co2_srad(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnsw(:,:) - aflx_upsw(:,:)) - & 
                  (cflx_dnsw(:,:) - cflx_upsw(:,:)) + &
                  (dflx_dnsw(:,:) - dflx_upsw(:,:)) - &
                  (bflx_dnsw(:,:) - bflx_upsw(:,:))) 
        dR_co2_traf(1:ldc%nproma,1:klev+1,krow)  = dR_co2_traf(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnlw_clr(:,:) - aflx_uplw_clr(:,:)) - & 
                  (cflx_dnlw_clr(:,:) - cflx_uplw_clr(:,:)) + &
                  (dflx_dnlw_clr(:,:) - dflx_uplw_clr(:,:)) - &
                  (bflx_dnlw_clr(:,:) - bflx_uplw_clr(:,:))) 
        dR_co2_sraf(1:ldc%nproma,1:klev+1,krow)  = dR_co2_sraf(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnsw_clr(:,:) - aflx_upsw_clr(:,:)) - & 
                  (cflx_dnsw_clr(:,:) - cflx_upsw_clr(:,:)) + &
                  (dflx_dnsw_clr(:,:) - dflx_upsw_clr(:,:)) - &
                  (bflx_dnsw_clr(:,:) - bflx_upsw_clr(:,:)))

      IF (prp3doutput) THEN
        dQ_co2_trad(1:ldc%nproma,1:klev,krow)  = dQ_co2_trad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnlw(:,1:klev) - aflx_uplw(:,1:klev)) - &
                                                    (aflx_dnlw(:,2:klev+1) - aflx_uplw(:,2:klev+1))) - &
                                                  ((cflx_dnlw(:,1:klev) - cflx_uplw(:,1:klev)) - & 
                                                    (cflx_dnlw(:,2:klev+1) - cflx_uplw(:,2:klev+1)))) - &
                                                  (((bflx_dnlw(:,1:klev) - bflx_uplw(:,1:klev)) - & 
                                                    (bflx_dnlw(:,2:klev+1) - bflx_uplw(:,2:klev+1))) - &
                                                  ((dflx_dnlw(:,1:klev) - dflx_uplw(:,1:klev)) - &
                                                    (dflx_dnlw(:,2:klev+1) - dflx_uplw(:,2:klev+1)))))
        dQ_co2_srad(1:ldc%nproma,1:klev,krow)  = dQ_co2_srad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnsw(:,1:klev) - aflx_upsw(:,1:klev)) - &
                                                    (aflx_dnsw(:,2:klev+1) - aflx_upsw(:,2:klev+1))) - &
                                                  ((cflx_dnsw(:,1:klev) - cflx_upsw(:,1:klev)) - & 
                                                    (cflx_dnsw(:,2:klev+1) - cflx_upsw(:,2:klev+1)))) - &
                                                  (((bflx_dnsw(:,1:klev) - bflx_upsw(:,1:klev)) - & 
                                                    (bflx_dnsw(:,2:klev+1) - bflx_upsw(:,2:klev+1))) - &
                                                  ((dflx_dnsw(:,1:klev) - dflx_upsw(:,1:klev)) - &
                                                    (dflx_dnsw(:,2:klev+1) - dflx_upsw(:,2:klev+1)))))
        dQ_co2_traf(1:ldc%nproma,1:klev,krow)  = dQ_co2_traf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnlw_clr(:,1:klev) - aflx_uplw_clr(:,1:klev)) - &
                                                    (aflx_dnlw_clr(:,2:klev+1) - aflx_uplw_clr(:,2:klev+1))) - &
                                                  ((cflx_dnlw_clr(:,1:klev) - cflx_uplw_clr(:,1:klev)) - & 
                                                    (cflx_dnlw_clr(:,2:klev+1) - cflx_uplw_clr(:,2:klev+1)))) - &
                                                  (((bflx_dnlw_clr(:,1:klev) - bflx_uplw_clr(:,1:klev)) - & 
                                                    (bflx_dnlw_clr(:,2:klev+1) - bflx_uplw_clr(:,2:klev+1))) - &
                                                  ((dflx_dnlw_clr(:,1:klev) - dflx_uplw_clr(:,1:klev)) - &
                                                    (dflx_dnlw_clr(:,2:klev+1) - dflx_uplw_clr(:,2:klev+1)))))
        dQ_co2_sraf(1:ldc%nproma,1:klev,krow)  = dQ_co2_sraf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnsw_clr(:,1:klev) - aflx_upsw_clr(:,1:klev)) - &
                                                    (aflx_dnsw_clr(:,2:klev+1) - aflx_upsw_clr(:,2:klev+1))) - &
                                                  ((cflx_dnsw_clr(:,1:klev) - cflx_upsw_clr(:,1:klev)) - & 
                                                    (cflx_dnsw_clr(:,2:klev+1) - cflx_upsw_clr(:,2:klev+1)))) - &
                                                  (((bflx_dnsw_clr(:,1:klev) - bflx_upsw_clr(:,1:klev)) - & 
                                                    (bflx_dnsw_clr(:,2:klev+1) - bflx_upsw_clr(:,2:klev+1))) - &
                                                  ((dflx_dnsw_clr(:,1:klev) - dflx_upsw_clr(:,1:klev)) - &
                                                    (dflx_dnsw_clr(:,2:klev+1) - dflx_upsw_clr(:,2:klev+1)))))
      END IF

        ! Calculates total aerosol forcing:
        CALL message('','total aerosol with current state')

        CALL psrad_interface( &
            & ziaero          ,kproma          ,kbdim            ,klev             ,& 
            & krow            ,ktrac           ,ktype            ,nb_sw            ,& 
            & loland          ,loglac          ,cemiss           ,cos_mu0          ,& 
            & pgeom1          ,alb_vis_dir     ,alb_nir_dir      ,alb_vis_dif      ,&
            & alb_nir_dif     ,pp_fl           ,pp_hl            ,pp_sfc           ,&
            & tk_fl           ,tk_hl           ,tk_sfc           ,xq_vap           ,&
            & xq_liq          ,xq_ice          ,zcdnc            ,cld_frc          ,&
            & cld_cvr         ,xm_o3           ,xm_co2           ,xm_ch4           ,&
            & xm_n2o          ,xm_cfc          ,xm_o2            ,pxtm1            ,&
            & cflx_uplw       ,cflx_uplw_clr   ,cflx_dnlw        ,cflx_dnlw_clr    ,&
            & cflx_upsw       ,cflx_upsw_clr   ,cflx_dnsw        ,cflx_dnsw_clr    ,&
            & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
            & par_dff_sfc                                                          )
        
        CALL message('','total aerosol with saved state')
        CALL psrad_interface( &
            & iaero           ,kproma          ,kbdim            ,klev             ,& 
            & krow            ,ktrac           ,zktype           ,nb_sw            ,&
            & zloland         ,zloglac         ,zcemiss          ,zcos_mu0         ,& 
            & zpgeom1         ,zalb_vis_dir    ,zalb_nir_dir     ,zalb_vis_dif     ,&
            & zalb_nir_dif    ,zpp_fl          ,zpp_hl           ,zpp_sfc          ,& 
            & ztk_fl          ,ztk_hl          ,ztk_sfc          ,zq_vap           ,&
            & zq_liq          ,zq_ice          ,cdnc             ,zcld_frc         ,&
            & zcld_cvr        ,zm_o3           ,zm_co2           ,zm_ch4           ,&
            & zm_n2o          ,zm_cfc          ,zm_o2            ,pxtm1            ,&
            & dflx_uplw       ,dflx_uplw_clr   ,dflx_dnlw        ,dflx_dnlw_clr    ,&
            & dflx_upsw       ,dflx_upsw_clr   ,dflx_dnsw        ,dflx_dnsw_clr    ,&
            & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
            & par_dff_sfc                                                          )

        dR_aer_trad0(1:ldc%nproma,krow)  = dR_aer_trad0(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw(:,1) - aflx_uplw(:,1)) - & 
                (cflx_dnlw(:,1) - cflx_uplw(:,1)) + &
                                                          (dflx_dnlw(:,1) - dflx_uplw(:,1)) - &
                                                          (bflx_dnlw(:,1) - bflx_uplw(:,1))) 
        dR_aer_srad0(1:ldc%nproma,krow)  = dR_aer_srad0(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw(:,1) - aflx_upsw(:,1)) - & 
                (cflx_dnsw(:,1) - cflx_upsw(:,1)) + &
                                                          (dflx_dnsw(:,1) - dflx_upsw(:,1)) - &
                                                          (bflx_dnsw(:,1) - bflx_upsw(:,1))) 
        dR_aer_trads(1:ldc%nproma,krow)  = dR_aer_trads(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw(:,klev+1) - aflx_uplw(:,klev+1)) - & 
                (cflx_dnlw(:,klev+1) - cflx_uplw(:,klev+1)) + &
                                                          (dflx_dnlw(:,klev+1) - dflx_uplw(:,klev+1)) - &
                                                          (bflx_dnlw(:,klev+1) - bflx_uplw(:,klev+1))) 
        dR_aer_srads(1:ldc%nproma,krow)  = dR_aer_srads(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw(:,klev+1) - aflx_upsw(:,klev+1)) - & 
                (cflx_dnsw(:,klev+1) - cflx_upsw(:,klev+1)) + &
                                                          (dflx_dnsw(:,klev+1) - dflx_upsw(:,klev+1)) - &
                                                          (bflx_dnsw(:,klev+1) - bflx_upsw(:,klev+1))) 
        dR_aer_traf0(1:ldc%nproma,krow)  = dR_aer_traf0(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw_clr(:,1) - aflx_uplw_clr(:,1)) - & 
                (cflx_dnlw_clr(:,1) - cflx_uplw_clr(:,1)) + &
                                                          (dflx_dnlw_clr(:,1) - dflx_uplw_clr(:,1)) - &
                                                          (bflx_dnlw_clr(:,1) - bflx_uplw_clr(:,1))) 
        dR_aer_sraf0(1:ldc%nproma,krow)  = dR_aer_sraf0(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw_clr(:,1) - aflx_upsw_clr(:,1)) - & 
                (cflx_dnsw_clr(:,1) - cflx_upsw_clr(:,1)) + &
                                                          (dflx_dnsw_clr(:,1) - dflx_upsw_clr(:,1)) - &
                                                          (bflx_dnsw_clr(:,1) - bflx_upsw_clr(:,1))) 
        dR_aer_trafs(1:ldc%nproma,krow)  = dR_aer_trafs(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw_clr(:,klev+1) - aflx_uplw_clr(:,klev+1)) - & 
                (cflx_dnlw_clr(:,klev+1) - cflx_uplw_clr(:,klev+1)) + &
                                                          (dflx_dnlw_clr(:,klev+1) - dflx_uplw_clr(:,klev+1)) - &
                                                          (bflx_dnlw_clr(:,klev+1) - bflx_uplw_clr(:,klev+1))) 
        dR_aer_srafs(1:ldc%nproma,krow)  = dR_aer_srafs(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw_clr(:,klev+1) - aflx_upsw_clr(:,klev+1)) - & 
                (cflx_dnsw_clr(:,klev+1) - cflx_upsw_clr(:,klev+1)) + &
                                                          (dflx_dnsw_clr(:,klev+1) - dflx_upsw_clr(:,klev+1)) - &
                                                          (bflx_dnsw_clr(:,klev+1) - bflx_upsw_clr(:,klev+1))) 

        dR_aer_trad(1:ldc%nproma,1:klev+1,krow)  = dR_aer_trad(1:ldc%nproma,1:klev+1,krow) + &
                                                    dtprps*0.5_wp*((aflx_dnlw(:,:) - aflx_uplw(:,:)) - & 
                  (cflx_dnlw(:,:) - cflx_uplw(:,:)) + &
                  (dflx_dnlw(:,:) - dflx_uplw(:,:)) - &
                  (bflx_dnlw(:,:) - bflx_uplw(:,:))) 
        dR_aer_srad(1:ldc%nproma,1:klev+1,krow)  = dR_aer_srad(1:ldc%nproma,1:klev+1,krow) + &
                                                    dtprps*0.5_wp*((aflx_dnsw(:,:) - aflx_upsw(:,:)) - & 
                  (cflx_dnsw(:,:) - cflx_upsw(:,:)) + &
                  (dflx_dnsw(:,:) - dflx_upsw(:,:)) - &
                  (bflx_dnsw(:,:) - bflx_upsw(:,:))) 
        dR_aer_traf(1:ldc%nproma,1:klev+1,krow)  = dR_aer_traf(1:ldc%nproma,1:klev+1,krow) + &
                                                    dtprps*0.5_wp*((aflx_dnlw_clr(:,:) - aflx_uplw_clr(:,:)) - & 
                  (cflx_dnlw_clr(:,:) - cflx_uplw_clr(:,:)) + &
                  (dflx_dnlw_clr(:,:) - dflx_uplw_clr(:,:)) - &
                  (bflx_dnlw_clr(:,:) - bflx_uplw_clr(:,:))) 
        dR_aer_sraf(1:ldc%nproma,1:klev+1,krow)  = dR_aer_sraf(1:ldc%nproma,1:klev+1,krow) + &
                                                    dtprps*0.5_wp*((aflx_dnsw_clr(:,:) - aflx_upsw_clr(:,:)) - & 
                  (cflx_dnsw_clr(:,:) - cflx_upsw_clr(:,:)) + &
                  (dflx_dnsw_clr(:,:) - dflx_upsw_clr(:,:)) - &
                  (bflx_dnsw_clr(:,:) - bflx_upsw_clr(:,:)))

      IF (prp3doutput) THEN
        dQ_aer_trad(1:ldc%nproma,1:klev,krow)  = dQ_aer_trad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                  ((((aflx_dnlw(:,1:klev) - aflx_uplw(:,1:klev)) - &
                                                    (aflx_dnlw(:,2:klev+1) - aflx_uplw(:,2:klev+1))) - &
                                                    ((cflx_dnlw(:,1:klev) - cflx_uplw(:,1:klev)) - & 
                                                    (cflx_dnlw(:,2:klev+1) - cflx_uplw(:,2:klev+1)))) - &
                                                  (((bflx_dnlw(:,1:klev) - bflx_uplw(:,1:klev)) - & 
                                                    (bflx_dnlw(:,2:klev+1) - bflx_uplw(:,2:klev+1))) - &
                                                    ((dflx_dnlw(:,1:klev) - dflx_uplw(:,1:klev)) - &
                                                    (dflx_dnlw(:,2:klev+1) - dflx_uplw(:,2:klev+1)))))
        dQ_aer_srad(1:ldc%nproma,1:klev,krow)  = dQ_aer_srad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                  ((((aflx_dnsw(:,1:klev) - aflx_upsw(:,1:klev)) - &
                                                    (aflx_dnsw(:,2:klev+1) - aflx_upsw(:,2:klev+1))) - &
                                                    ((cflx_dnsw(:,1:klev) - cflx_upsw(:,1:klev)) - & 
                                                    (cflx_dnsw(:,2:klev+1) - cflx_upsw(:,2:klev+1)))) - &
                                                  (((bflx_dnsw(:,1:klev) - bflx_upsw(:,1:klev)) - & 
                                                    (bflx_dnsw(:,2:klev+1) - bflx_upsw(:,2:klev+1))) - &
                                                    ((dflx_dnsw(:,1:klev) - dflx_upsw(:,1:klev)) - &
                                                    (dflx_dnsw(:,2:klev+1) - dflx_upsw(:,2:klev+1)))))
        IF (prpraf) THEN
        dQ_aer_traf(1:ldc%nproma,1:klev,krow)  = dQ_aer_traf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                  ((((aflx_dnlw_clr(:,1:klev) - aflx_uplw_clr(:,1:klev)) - &
                                                    (aflx_dnlw_clr(:,2:klev+1) - aflx_uplw_clr(:,2:klev+1))) - &
                                                    ((cflx_dnlw_clr(:,1:klev) - cflx_uplw_clr(:,1:klev)) - & 
                                                    (cflx_dnlw_clr(:,2:klev+1) - cflx_uplw_clr(:,2:klev+1)))) - &
                                                  (((bflx_dnlw_clr(:,1:klev) - bflx_uplw_clr(:,1:klev)) - & 
                                                    (bflx_dnlw_clr(:,2:klev+1) - bflx_uplw_clr(:,2:klev+1))) - &
                                                    ((dflx_dnlw_clr(:,1:klev) - dflx_uplw_clr(:,1:klev)) - &
                                                    (dflx_dnlw_clr(:,2:klev+1) - dflx_uplw_clr(:,2:klev+1)))))
        dQ_aer_sraf(1:ldc%nproma,1:klev,krow)  = dQ_aer_sraf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                  ((((aflx_dnsw_clr(:,1:klev) - aflx_upsw_clr(:,1:klev)) - &
                                                    (aflx_dnsw_clr(:,2:klev+1) - aflx_upsw_clr(:,2:klev+1))) - &
                                                    ((cflx_dnsw_clr(:,1:klev) - cflx_upsw_clr(:,1:klev)) - & 
                                                    (cflx_dnsw_clr(:,2:klev+1) - cflx_upsw_clr(:,2:klev+1)))) - &
                                                  (((bflx_dnsw_clr(:,1:klev) - bflx_upsw_clr(:,1:klev)) - & 
                                                    (bflx_dnsw_clr(:,2:klev+1) - bflx_upsw_clr(:,2:klev+1))) - &
                                                    ((dflx_dnsw_clr(:,1:klev) - dflx_upsw_clr(:,1:klev)) - &
                                                    (dflx_dnsw_clr(:,2:klev+1) - dflx_upsw_clr(:,2:klev+1)))))
        END IF
      END IF
    END IF  
  
  IF (prp_sp) THEN
      IF (iaero(1)==8 .or. iaero(1)==9) THEN   ! If scenario contains simple-plumes aerosols, calculates cloud feedback without Twomey effect of simple-plumes, 
                                         ! and calculates direct, indirect and total effets of simple-plumes

          ! Calculate cloud feedback only (without Twomey):
      CALL message('','cloud feedback')
      CALL psrad_interface( &
          & iaero           ,kproma          ,kbdim            ,klev             ,& 
          & krow            ,ktrac           ,ktype            ,nb_sw            ,& 
          & loland          ,loglac          ,cemiss           ,cos_mu0          ,& 
          & pgeom1          ,alb_vis_dir     ,alb_nir_dir      ,alb_vis_dif      ,&
          & alb_nir_dif     ,pp_fl           ,pp_hl            ,pp_sfc           ,&
          & tk_fl           ,tk_hl           ,tk_sfc           ,xq_vap           ,&
          & zq_liq          ,zq_ice          ,cdnc             ,zcld_frc         ,&
          & zcld_cvr        ,xm_o3           ,xm_co2           ,xm_ch4           ,&
          & xm_n2o          ,xm_cfc          ,xm_o2            ,pxtm1            ,&
          & cflx_uplw       ,cflx_uplw_clr   ,cflx_dnlw        ,cflx_dnlw_clr    ,&
          & cflx_upsw       ,cflx_upsw_clr   ,cflx_dnsw        ,cflx_dnsw_clr    ,&
          & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
          & par_dff_sfc                                                          )

      CALL psrad_interface( &
          & ziaero          ,kproma          ,kbdim            ,klev             ,& 
          & krow            ,ktrac           ,zktype           ,nb_sw            ,&
          & zloland         ,zloglac         ,zcemiss          ,zcos_mu0         ,& 
          & zpgeom1         ,zalb_vis_dir    ,zalb_nir_dir     ,zalb_vis_dif     ,&
          & zalb_nir_dif    ,zpp_fl          ,zpp_hl           ,zpp_sfc          ,& 
          & ztk_fl          ,ztk_hl          ,ztk_sfc          ,zq_vap           ,&
          & xq_liq          ,xq_ice          ,zcdnc            ,cld_frc          ,&
          & cld_cvr         ,zm_o3           ,zm_co2           ,zm_ch4           ,&
          & zm_n2o          ,zm_cfc          ,zm_o2            ,pxtm1            ,&
          & dflx_uplw       ,dflx_uplw_clr   ,dflx_dnlw        ,dflx_dnlw_clr    ,&
          & dflx_upsw       ,dflx_upsw_clr   ,dflx_dnsw        ,dflx_dnsw_clr    ,&
          & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
          & par_dff_sfc                                                          )

      dR_cdd_trad0(1:ldc%nproma,krow)  = dR_cdd_trad0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw(:,1) - aflx_uplw(:,1)) - & 
              (cflx_dnlw(:,1) - cflx_uplw(:,1)) + &
                                                        (dflx_dnlw(:,1) - dflx_uplw(:,1)) - &
                                                        (bflx_dnlw(:,1) - bflx_uplw(:,1))) 
      dR_cdd_srad0(1:ldc%nproma,krow)  = dR_cdd_srad0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw(:,1) - aflx_upsw(:,1)) - & 
              (cflx_dnsw(:,1) - cflx_upsw(:,1)) + &
                                                        (dflx_dnsw(:,1) - dflx_upsw(:,1)) - &
                                                        (bflx_dnsw(:,1) - bflx_upsw(:,1))) 
    IF (prprads) THEN
      dR_cdd_trads(1:ldc%nproma,krow)  = dR_cdd_trads(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw(:,klev+1) - aflx_uplw(:,klev+1)) - & 
              (cflx_dnlw(:,klev+1) - cflx_uplw(:,klev+1)) + &
                                                        (dflx_dnlw(:,klev+1) - dflx_uplw(:,klev+1)) - &
                                                        (bflx_dnlw(:,klev+1) - bflx_uplw(:,klev+1))) 
      dR_cdd_srads(1:ldc%nproma,krow)  = dR_cdd_srads(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw(:,klev+1) - aflx_upsw(:,klev+1)) - & 
              (cflx_dnsw(:,klev+1) - cflx_upsw(:,klev+1)) + &
                                                        (dflx_dnsw(:,klev+1) - dflx_upsw(:,klev+1)) - &
                                                        (bflx_dnsw(:,klev+1) - bflx_upsw(:,klev+1))) 
    END IF
    IF (prpraf) THEN
      dR_cdd_traf0(1:ldc%nproma,krow)  = dR_cdd_traf0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw_clr(:,1) - aflx_uplw_clr(:,1)) - & 
              (cflx_dnlw_clr(:,1) - cflx_uplw_clr(:,1)) + &
                                                        (dflx_dnlw_clr(:,1) - dflx_uplw_clr(:,1)) - &
                                                        (bflx_dnlw_clr(:,1) - bflx_uplw_clr(:,1))) 
      dR_cdd_sraf0(1:ldc%nproma,krow)  = dR_cdd_sraf0(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw_clr(:,1) - aflx_upsw_clr(:,1)) - & 
              (cflx_dnsw_clr(:,1) - cflx_upsw_clr(:,1)) + &
                                                        (dflx_dnsw_clr(:,1) - dflx_upsw_clr(:,1)) - &
                                                        (bflx_dnsw_clr(:,1) - bflx_upsw_clr(:,1))) 
      dR_cdd_trafs(1:ldc%nproma,krow)  = dR_cdd_trafs(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnlw_clr(:,klev+1) - aflx_uplw_clr(:,klev+1)) - & 
              (cflx_dnlw_clr(:,klev+1) - cflx_uplw_clr(:,klev+1)) + &
                                                        (dflx_dnlw_clr(:,klev+1) - dflx_uplw_clr(:,klev+1)) - &
                                                        (bflx_dnlw_clr(:,klev+1) - bflx_uplw_clr(:,klev+1))) 
      dR_cdd_srafs(1:ldc%nproma,krow)  = dR_cdd_srafs(1:ldc%nproma,krow) + &
                                          dtprps*0.5_wp*((aflx_dnsw_clr(:,klev+1) - aflx_upsw_clr(:,klev+1)) - & 
              (cflx_dnsw_clr(:,klev+1) - cflx_upsw_clr(:,klev+1)) + &
                                                        (dflx_dnsw_clr(:,klev+1) - dflx_upsw_clr(:,klev+1)) - &
                                                        (bflx_dnsw_clr(:,klev+1) - bflx_upsw_clr(:,klev+1))) 
    END IF
    IF (prp3doutput) THEN
      IF (prp3drad) THEN
      dR_cdd_trad(1:ldc%nproma,1:klev+1,krow)  = dR_cdd_trad(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnlw(:,:) - aflx_uplw(:,:)) - & 
                (cflx_dnlw(:,:) - cflx_uplw(:,:)) + &
                (dflx_dnlw(:,:) - dflx_uplw(:,:)) - &
                (bflx_dnlw(:,:) - bflx_uplw(:,:))) 
      dR_cdd_srad(1:ldc%nproma,1:klev+1,krow)  = dR_cdd_srad(1:ldc%nproma,1:klev+1,krow) + &
                                                  dtprps*0.5_wp*((aflx_dnsw(:,:) - aflx_upsw(:,:)) - & 
                (cflx_dnsw(:,:) - cflx_upsw(:,:)) + &
                (dflx_dnsw(:,:) - dflx_upsw(:,:)) - &
                (bflx_dnsw(:,:) - bflx_upsw(:,:))) 
        IF (prpraf) THEN
        dR_cdd_traf(1:ldc%nproma,1:klev+1,krow)  = dR_cdd_traf(1:ldc%nproma,1:klev+1,krow) + &
                                                    dtprps*0.5_wp*((aflx_dnlw_clr(:,:) - aflx_uplw_clr(:,:)) - & 
                  (cflx_dnlw_clr(:,:) - cflx_uplw_clr(:,:)) + &
                  (dflx_dnlw_clr(:,:) - dflx_uplw_clr(:,:)) - &
                  (bflx_dnlw_clr(:,:) - bflx_uplw_clr(:,:))) 
        dR_cdd_sraf(1:ldc%nproma,1:klev+1,krow)  = dR_cdd_sraf(1:ldc%nproma,1:klev+1,krow) + &
                                                    dtprps*0.5_wp*((aflx_dnsw_clr(:,:) - aflx_upsw_clr(:,:)) - & 
                  (cflx_dnsw_clr(:,:) - cflx_upsw_clr(:,:)) + &
                  (dflx_dnsw_clr(:,:) - dflx_upsw_clr(:,:)) - &
                  (bflx_dnsw_clr(:,:) - bflx_upsw_clr(:,:)))
        END IF
      END IF
      dQ_cdd_trad(1:ldc%nproma,1:klev,krow)  = dQ_cdd_trad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnlw(:,1:klev) - aflx_uplw(:,1:klev)) - &
                                                  (aflx_dnlw(:,2:klev+1) - aflx_uplw(:,2:klev+1))) - &
                                                  ((cflx_dnlw(:,1:klev) - cflx_uplw(:,1:klev)) - & 
                                                  (cflx_dnlw(:,2:klev+1) - cflx_uplw(:,2:klev+1)))) - &
                                                (((bflx_dnlw(:,1:klev) - bflx_uplw(:,1:klev)) - & 
                                                  (bflx_dnlw(:,2:klev+1) - bflx_uplw(:,2:klev+1))) - &
                                                  ((dflx_dnlw(:,1:klev) - dflx_uplw(:,1:klev)) - &
                                                  (dflx_dnlw(:,2:klev+1) - dflx_uplw(:,2:klev+1)))))
      dQ_cdd_srad(1:ldc%nproma,1:klev,krow)  = dQ_cdd_srad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnsw(:,1:klev) - aflx_upsw(:,1:klev)) - &
                                                  (aflx_dnsw(:,2:klev+1) - aflx_upsw(:,2:klev+1))) - &
                                                  ((cflx_dnsw(:,1:klev) - cflx_upsw(:,1:klev)) - & 
                                                  (cflx_dnsw(:,2:klev+1) - cflx_upsw(:,2:klev+1)))) - &
                                                (((bflx_dnsw(:,1:klev) - bflx_upsw(:,1:klev)) - & 
                                                  (bflx_dnsw(:,2:klev+1) - bflx_upsw(:,2:klev+1))) - &
                                                  ((dflx_dnsw(:,1:klev) - dflx_upsw(:,1:klev)) - &
                                                  (dflx_dnsw(:,2:klev+1) - dflx_upsw(:,2:klev+1)))))
      IF (prpraf) THEN
      dQ_cdd_traf(1:ldc%nproma,1:klev,krow)  = dQ_cdd_traf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnlw_clr(:,1:klev) - aflx_uplw_clr(:,1:klev)) - &
                                                  (aflx_dnlw_clr(:,2:klev+1) - aflx_uplw_clr(:,2:klev+1))) - &
                                                  ((cflx_dnlw_clr(:,1:klev) - cflx_uplw_clr(:,1:klev)) - & 
                                                  (cflx_dnlw_clr(:,2:klev+1) - cflx_uplw_clr(:,2:klev+1)))) - &
                                                (((bflx_dnlw_clr(:,1:klev) - bflx_uplw_clr(:,1:klev)) - & 
                                                  (bflx_dnlw_clr(:,2:klev+1) - bflx_uplw_clr(:,2:klev+1))) - &
                                                  ((dflx_dnlw_clr(:,1:klev) - dflx_uplw_clr(:,1:klev)) - &
                                                  (dflx_dnlw_clr(:,2:klev+1) - dflx_uplw_clr(:,2:klev+1)))))
      dQ_cdd_sraf(1:ldc%nproma,1:klev,krow)  = dQ_cdd_sraf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                ((((aflx_dnsw_clr(:,1:klev) - aflx_upsw_clr(:,1:klev)) - &
                                                  (aflx_dnsw_clr(:,2:klev+1) - aflx_upsw_clr(:,2:klev+1))) - &
                                                  ((cflx_dnsw_clr(:,1:klev) - cflx_upsw_clr(:,1:klev)) - & 
                                                  (cflx_dnsw_clr(:,2:klev+1) - cflx_upsw_clr(:,2:klev+1)))) - &
                                                (((bflx_dnsw_clr(:,1:klev) - bflx_upsw_clr(:,1:klev)) - & 
                                                  (bflx_dnsw_clr(:,2:klev+1) - bflx_upsw_clr(:,2:klev+1))) - &
                                                  ((dflx_dnsw_clr(:,1:klev) - dflx_upsw_clr(:,1:klev)) - &
                                                  (dflx_dnsw_clr(:,2:klev+1) - dflx_upsw_clr(:,2:klev+1)))))
      END IF
      END IF

          ! Calculates total simple-plumes forcing:
        iaero_direct(2) = 10
        iaero_indirect(2) = 10
        iaero_sp(2) = 10
        SELECT CASE (iaero(1))
        CASE (0)
          iaero_sp_ref(1) = 0
        CASE (8)
          iaero_sp_ref(1) = 5    
        CASE (9)
          iaero_sp_ref(1) = 3 
        CASE DEFAULT
        END SELECT

        CALL message('','total anthropogenic aerosol with current state')
        CALL psrad_interface( &
            & iaero_sp_ref    ,kproma          ,kbdim            ,klev             ,& 
            & krow            ,ktrac           ,ktype            ,nb_sw            ,& 
            & loland          ,loglac          ,cemiss           ,cos_mu0          ,& 
            & pgeom1          ,alb_vis_dir     ,alb_nir_dir      ,alb_vis_dif      ,&
            & alb_nir_dif     ,pp_fl           ,pp_hl            ,pp_sfc           ,&
            & tk_fl           ,tk_hl           ,tk_sfc           ,xq_vap           ,&
            & xq_liq          ,xq_ice          ,zcdnc            ,cld_frc          ,&
            & cld_cvr         ,xm_o3           ,xm_co2           ,xm_ch4           ,&
            & xm_n2o          ,xm_cfc          ,xm_o2            ,pxtm1            ,&
            & cflx_uplw       ,cflx_uplw_clr   ,cflx_dnlw        ,cflx_dnlw_clr    ,&
            & cflx_upsw       ,cflx_upsw_clr   ,cflx_dnsw        ,cflx_dnsw_clr    ,&
            & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
            & par_dff_sfc                                                          )
        
        CALL message('','total aerosol with saved state')
        CALL psrad_interface( &
            & iaero_sp        ,kproma          ,kbdim            ,klev             ,& 
            & krow            ,ktrac           ,zktype           ,nb_sw            ,&
            & zloland         ,zloglac         ,zcemiss          ,zcos_mu0         ,& 
            & zpgeom1         ,zalb_vis_dir    ,zalb_nir_dir     ,zalb_vis_dif     ,&
            & zalb_nir_dif    ,zpp_fl          ,zpp_hl           ,zpp_sfc          ,& 
            & ztk_fl          ,ztk_hl          ,ztk_sfc          ,zq_vap           ,&
            & zq_liq          ,zq_ice          ,cdnc             ,zcld_frc         ,&
            & zcld_cvr        ,zm_o3           ,zm_co2           ,zm_ch4           ,&
            & zm_n2o          ,zm_cfc          ,zm_o2            ,pxtm1            ,&
            & dflx_uplw       ,dflx_uplw_clr   ,dflx_dnlw        ,dflx_dnlw_clr    ,&
            & dflx_upsw       ,dflx_upsw_clr   ,dflx_dnsw        ,dflx_dnsw_clr    ,&
            & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
            & par_dff_sfc                                                          )

        dR_sp_trad0(1:ldc%nproma,krow)  = dR_sp_trad0(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw(:,1) - aflx_uplw(:,1)) - & 
                (cflx_dnlw(:,1) - cflx_uplw(:,1)) + &
                                                          (dflx_dnlw(:,1) - dflx_uplw(:,1)) - &
                                                          (bflx_dnlw(:,1) - bflx_uplw(:,1))) 
        dR_sp_srad0(1:ldc%nproma,krow)  = dR_sp_srad0(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw(:,1) - aflx_upsw(:,1)) - & 
                (cflx_dnsw(:,1) - cflx_upsw(:,1)) + &
                                                          (dflx_dnsw(:,1) - dflx_upsw(:,1)) - &
                                                          (bflx_dnsw(:,1) - bflx_upsw(:,1))) 
      IF (prprads) THEN
        dR_sp_trads(1:ldc%nproma,krow)  = dR_sp_trads(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw(:,klev+1) - aflx_uplw(:,klev+1)) - & 
                (cflx_dnlw(:,klev+1) - cflx_uplw(:,klev+1)) + &
                                                          (dflx_dnlw(:,klev+1) - dflx_uplw(:,klev+1)) - &
                                                          (bflx_dnlw(:,klev+1) - bflx_uplw(:,klev+1))) 
        dR_sp_srads(1:ldc%nproma,krow)  = dR_sp_srads(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw(:,klev+1) - aflx_upsw(:,klev+1)) - & 
                (cflx_dnsw(:,klev+1) - cflx_upsw(:,klev+1)) + &
                                                          (dflx_dnsw(:,klev+1) - dflx_upsw(:,klev+1)) - &
                                                          (bflx_dnsw(:,klev+1) - bflx_upsw(:,klev+1))) 
      END IF
      IF (prpraf) THEN
        dR_sp_traf0(1:ldc%nproma,krow)  = dR_sp_traf0(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw_clr(:,1) - aflx_uplw_clr(:,1)) - & 
                (cflx_dnlw_clr(:,1) - cflx_uplw_clr(:,1)) + &
                                                          (dflx_dnlw_clr(:,1) - dflx_uplw_clr(:,1)) - &
                                                          (bflx_dnlw_clr(:,1) - bflx_uplw_clr(:,1))) 
        dR_sp_sraf0(1:ldc%nproma,krow)  = dR_sp_sraf0(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw_clr(:,1) - aflx_upsw_clr(:,1)) - & 
                (cflx_dnsw_clr(:,1) - cflx_upsw_clr(:,1)) + &
                                                         (dflx_dnsw_clr(:,1) - dflx_upsw_clr(:,1)) - &
                                                          (bflx_dnsw_clr(:,1) - bflx_upsw_clr(:,1))) 
        dR_sp_trafs(1:ldc%nproma,krow)  = dR_sp_trafs(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw_clr(:,klev+1) - aflx_uplw_clr(:,klev+1)) - & 
                (cflx_dnlw_clr(:,klev+1) - cflx_uplw_clr(:,klev+1)) + &
                                                          (dflx_dnlw_clr(:,klev+1) - dflx_uplw_clr(:,klev+1)) - &
                                                          (bflx_dnlw_clr(:,klev+1) - bflx_uplw_clr(:,klev+1))) 
        dR_sp_srafs(1:ldc%nproma,krow)  = dR_sp_srafs(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw_clr(:,klev+1) - aflx_upsw_clr(:,klev+1)) - & 
                (cflx_dnsw_clr(:,klev+1) - cflx_upsw_clr(:,klev+1)) + &
                                                         (dflx_dnsw_clr(:,klev+1) - dflx_upsw_clr(:,klev+1)) - &
                                                          (bflx_dnsw_clr(:,klev+1) - bflx_upsw_clr(:,klev+1))) 
      END IF
      IF (prp3doutput) THEN
        IF (prp3drad) THEN
        dR_sp_trad(1:ldc%nproma,1:klev+1,krow)  = dR_sp_trad(1:ldc%nproma,1:klev+1,krow) + &
                                                    dtprps*0.5_wp*((aflx_dnlw(:,:) - aflx_uplw(:,:)) - & 
                  (cflx_dnlw(:,:) - cflx_uplw(:,:)) + &
                  (dflx_dnlw(:,:) - dflx_uplw(:,:)) - &
                  (bflx_dnlw(:,:) - bflx_uplw(:,:))) 
        dR_sp_srad(1:ldc%nproma,1:klev+1,krow)  = dR_sp_srad(1:ldc%nproma,1:klev+1,krow) + &
                                                    dtprps*0.5_wp*((aflx_dnsw(:,:) - aflx_upsw(:,:)) - & 
                  (cflx_dnsw(:,:) - cflx_upsw(:,:)) + &
                  (dflx_dnsw(:,:) - dflx_upsw(:,:)) - &
                  (bflx_dnsw(:,:) - bflx_upsw(:,:))) 
          IF (prpraf) THEN
            dR_sp_traf(1:ldc%nproma,1:klev+1,krow)  = dR_sp_traf(1:ldc%nproma,1:klev+1,krow) + &
                                                                dtprps*0.5_wp*((aflx_dnlw_clr(:,:) - aflx_uplw_clr(:,:)) - & 
                                                        (cflx_dnlw_clr(:,:) - cflx_uplw_clr(:,:)) + &
                                                        (dflx_dnlw_clr(:,:) - dflx_uplw_clr(:,:)) - &
                                                        (bflx_dnlw_clr(:,:) - bflx_uplw_clr(:,:))) 
            dR_sp_sraf(1:ldc%nproma,1:klev+1,krow)  = dR_sp_sraf(1:ldc%nproma,1:klev+1,krow) + &
                                                                dtprps*0.5_wp*((aflx_dnsw_clr(:,:) - aflx_upsw_clr(:,:)) - & 
                                                        (cflx_dnsw_clr(:,:) - cflx_upsw_clr(:,:)) + &
                                                        (dflx_dnsw_clr(:,:) - dflx_upsw_clr(:,:)) - &
                                                                  (bflx_dnsw_clr(:,:) - bflx_upsw_clr(:,:)))
          END IF
        END IF
        dQ_sp_trad(1:ldc%nproma,1:klev,krow)  = dQ_sp_trad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                  ((((aflx_dnlw(:,1:klev) - aflx_uplw(:,1:klev)) - &
                                                    (aflx_dnlw(:,2:klev+1) - aflx_uplw(:,2:klev+1))) - &
                                                    ((cflx_dnlw(:,1:klev) - cflx_uplw(:,1:klev)) - & 
                                                    (cflx_dnlw(:,2:klev+1) - cflx_uplw(:,2:klev+1)))) - &
                                                  (((bflx_dnlw(:,1:klev) - bflx_uplw(:,1:klev)) - & 
                                                    (bflx_dnlw(:,2:klev+1) - bflx_uplw(:,2:klev+1))) - &
                                                    ((dflx_dnlw(:,1:klev) - dflx_uplw(:,1:klev)) - &
                                                    (dflx_dnlw(:,2:klev+1) - dflx_uplw(:,2:klev+1)))))
        dQ_sp_srad(1:ldc%nproma,1:klev,krow)  = dQ_sp_srad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                  ((((aflx_dnsw(:,1:klev) - aflx_upsw(:,1:klev)) - &
                                                    (aflx_dnsw(:,2:klev+1) - aflx_upsw(:,2:klev+1))) - &
                                                    ((cflx_dnsw(:,1:klev) - cflx_upsw(:,1:klev)) - & 
                                                    (cflx_dnsw(:,2:klev+1) - cflx_upsw(:,2:klev+1)))) - &
                                                  (((bflx_dnsw(:,1:klev) - bflx_upsw(:,1:klev)) - & 
                                                    (bflx_dnsw(:,2:klev+1) - bflx_upsw(:,2:klev+1))) - &
                                                    ((dflx_dnsw(:,1:klev) - dflx_upsw(:,1:klev)) - &
                                                    (dflx_dnsw(:,2:klev+1) - dflx_upsw(:,2:klev+1)))))
        IF (prpraf) THEN
        dQ_sp_traf(1:ldc%nproma,1:klev,krow)  = dQ_sp_traf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                  ((((aflx_dnlw_clr(:,1:klev) - aflx_uplw_clr(:,1:klev)) - &
                                                    (aflx_dnlw_clr(:,2:klev+1) - aflx_uplw_clr(:,2:klev+1))) - &
                                                    ((cflx_dnlw_clr(:,1:klev) - cflx_uplw_clr(:,1:klev)) - & 
                                                    (cflx_dnlw_clr(:,2:klev+1) - cflx_uplw_clr(:,2:klev+1)))) - &
                                                  (((bflx_dnlw_clr(:,1:klev) - bflx_uplw_clr(:,1:klev)) - & 
                                                    (bflx_dnlw_clr(:,2:klev+1) - bflx_uplw_clr(:,2:klev+1))) - &
                                                    ((dflx_dnlw_clr(:,1:klev) - dflx_uplw_clr(:,1:klev)) - &
                                                    (dflx_dnlw_clr(:,2:klev+1) - dflx_uplw_clr(:,2:klev+1)))))
        dQ_sp_sraf(1:ldc%nproma,1:klev,krow)  = dQ_sp_sraf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                  ((((aflx_dnsw_clr(:,1:klev) - aflx_upsw_clr(:,1:klev)) - &
                                                    (aflx_dnsw_clr(:,2:klev+1) - aflx_upsw_clr(:,2:klev+1))) - &
                                                    ((cflx_dnsw_clr(:,1:klev) - cflx_upsw_clr(:,1:klev)) - & 
                                                    (cflx_dnsw_clr(:,2:klev+1) - cflx_upsw_clr(:,2:klev+1)))) - &
                                                  (((bflx_dnsw_clr(:,1:klev) - bflx_upsw_clr(:,1:klev)) - & 
                                                    (bflx_dnsw_clr(:,2:klev+1) - bflx_upsw_clr(:,2:klev+1))) - &
                                                    ((dflx_dnsw_clr(:,1:klev) - dflx_upsw_clr(:,1:klev)) - &
                                                    (dflx_dnsw_clr(:,2:klev+1) - dflx_upsw_clr(:,2:klev+1)))))
        END IF
      END IF

      ! Calculates direct-effect of simple-plumes:
        SELECT CASE (iaero(1))
        CASE (0)
          iaero_direct(1) = 0
          iaero_indirect(1) = 0
        CASE (8)
          iaero_direct(1) = 10              ! Only direct-effect of SP in state d (no volcano)
          iaero_indirect(1) = 13            ! Only indirect-effect of SP in state c, but if volcanoes -> added to aod
        CASE (9)
          iaero_direct(1) = 10              ! If experiment with SP only
          iaero_indirect(1) = 11
        CASE DEFAULT
        END SELECT

        CALL psrad_interface( &
            & iaero_indirect  ,kproma          ,kbdim            ,klev             ,& 
            & krow            ,ktrac           ,ktype            ,nb_sw            ,& 
            & loland          ,loglac          ,cemiss           ,cos_mu0          ,& 
            & pgeom1          ,alb_vis_dir     ,alb_nir_dir      ,alb_vis_dif      ,&
            & alb_nir_dif     ,pp_fl           ,pp_hl            ,pp_sfc           ,&
            & tk_fl           ,tk_hl           ,tk_sfc           ,xq_vap           ,&
            & xq_liq          ,xq_ice          ,cdnc             ,cld_frc          ,&
            & cld_cvr         ,xm_o3           ,xm_co2           ,xm_ch4           ,&
            & xm_n2o          ,xm_cfc          ,xm_o2            ,pxtm1            ,&
            & cflx_uplw       ,cflx_uplw_clr   ,cflx_dnlw        ,cflx_dnlw_clr    ,&
            & cflx_upsw       ,cflx_upsw_clr   ,cflx_dnsw        ,cflx_dnsw_clr    ,&
            & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
            & par_dff_sfc                                                          )
       ! CALL message('','direct-effect with saved state')
        CALL psrad_interface( &
            & iaero_direct    ,kproma          ,kbdim            ,klev             ,& 
            & krow            ,ktrac           ,zktype           ,nb_sw            ,&
            & zloland         ,zloglac         ,zcemiss          ,zcos_mu0         ,& 
            & zpgeom1         ,zalb_vis_dir    ,zalb_nir_dir     ,zalb_vis_dif     ,&
            & zalb_nir_dif    ,zpp_fl          ,zpp_hl           ,zpp_sfc          ,& 
            & ztk_fl          ,ztk_hl          ,ztk_sfc          ,zq_vap           ,&
            & zq_liq          ,zq_ice          ,zcdnc            ,zcld_frc         ,&
            & zcld_cvr        ,zm_o3           ,zm_co2           ,zm_ch4           ,&
            & zm_n2o          ,zm_cfc          ,zm_o2            ,pxtm1            ,&
            & dflx_uplw       ,dflx_uplw_clr   ,dflx_dnlw        ,dflx_dnlw_clr    ,&
            & dflx_upsw       ,dflx_upsw_clr   ,dflx_dnsw        ,dflx_dnsw_clr    ,&
            & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
            & par_dff_sfc                                                          )

        dR_spd_trad0(1:ldc%nproma,krow)  = dR_spd_trad0(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw(:,1) - aflx_uplw(:,1)) - & 
                (cflx_dnlw(:,1) - cflx_uplw(:,1)) + &
                                                          (dflx_dnlw(:,1) - dflx_uplw(:,1)) - &
                                                          (bflx_dnlw(:,1) - bflx_uplw(:,1))) 
        dR_spd_srad0(1:ldc%nproma,krow)  = dR_spd_srad0(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw(:,1) - aflx_upsw(:,1)) - & 
                (cflx_dnsw(:,1) - cflx_upsw(:,1)) + &
                                                          (dflx_dnsw(:,1) - dflx_upsw(:,1)) - &
                                                          (bflx_dnsw(:,1) - bflx_upsw(:,1))) 
      IF (prprads) THEN
        dR_spd_trads(1:ldc%nproma,krow)  = dR_spd_trads(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw(:,klev+1) - aflx_uplw(:,klev+1)) - & 
                (cflx_dnlw(:,klev+1) - cflx_uplw(:,klev+1)) + &
                                                          (dflx_dnlw(:,klev+1) - dflx_uplw(:,klev+1)) - &
                                                          (bflx_dnlw(:,klev+1) - bflx_uplw(:,klev+1))) 
        dR_spd_srads(1:ldc%nproma,krow)  = dR_spd_srads(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw(:,klev+1) - aflx_upsw(:,klev+1)) - & 
                (cflx_dnsw(:,klev+1) - cflx_upsw(:,klev+1)) + &
                                                          (dflx_dnsw(:,klev+1) - dflx_upsw(:,klev+1)) - &
                                                          (bflx_dnsw(:,klev+1) - bflx_upsw(:,klev+1))) 
      END IF
      IF (prpraf) THEN
        dR_spd_traf0(1:ldc%nproma,krow)  = dR_spd_traf0(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw_clr(:,1) - aflx_uplw_clr(:,1)) - & 
                (cflx_dnlw_clr(:,1) - cflx_uplw_clr(:,1)) + &
                                                          (dflx_dnlw_clr(:,1) - dflx_uplw_clr(:,1)) - &
                                                          (bflx_dnlw_clr(:,1) - bflx_uplw_clr(:,1))) 
        dR_spd_sraf0(1:ldc%nproma,krow)  = dR_spd_sraf0(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw_clr(:,1) - aflx_upsw_clr(:,1)) - & 
                (cflx_dnsw_clr(:,1) - cflx_upsw_clr(:,1)) + &
                                                          (dflx_dnsw_clr(:,1) - dflx_upsw_clr(:,1)) - &
                                                          (bflx_dnsw_clr(:,1) - bflx_upsw_clr(:,1))) 
        dR_spd_trafs(1:ldc%nproma,krow)  = dR_spd_trafs(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw_clr(:,klev+1) - aflx_uplw_clr(:,klev+1)) - & 
                (cflx_dnlw_clr(:,klev+1) - cflx_uplw_clr(:,klev+1)) + &
                                                          (dflx_dnlw_clr(:,klev+1) - dflx_uplw_clr(:,klev+1)) - &
                                                          (bflx_dnlw_clr(:,klev+1) - bflx_uplw_clr(:,klev+1))) 
        dR_spd_srafs(1:ldc%nproma,krow)  = dR_spd_srafs(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw_clr(:,klev+1) - aflx_upsw_clr(:,klev+1)) - & 
                (cflx_dnsw_clr(:,klev+1) - cflx_upsw_clr(:,klev+1)) + &
                                                          (dflx_dnsw_clr(:,klev+1) - dflx_upsw_clr(:,klev+1)) - &
                                                          (bflx_dnsw_clr(:,klev+1) - bflx_upsw_clr(:,klev+1))) 
      END IF
      IF (prp3doutput) THEN
        IF (prp3drad) THEN
        dR_spd_trad(1:ldc%nproma,1:klev+1,krow)  = dR_spd_trad(1:ldc%nproma,1:klev+1,krow) + &
                                                    dtprps*0.5_wp*((aflx_dnlw(:,:) - aflx_uplw(:,:)) - & 
                  (cflx_dnlw(:,:) - cflx_uplw(:,:)) + &
                  (dflx_dnlw(:,:) - dflx_uplw(:,:)) - &
                  (bflx_dnlw(:,:) - bflx_uplw(:,:))) 
        dR_spd_srad(1:ldc%nproma,1:klev+1,krow)  = dR_spd_srad(1:ldc%nproma,1:klev+1,krow) + &
                                                    dtprps*0.5_wp*((aflx_dnsw(:,:) - aflx_upsw(:,:)) - & 
                  (cflx_dnsw(:,:) - cflx_upsw(:,:)) + &
                  (dflx_dnsw(:,:) - dflx_upsw(:,:)) - &
                  (bflx_dnsw(:,:) - bflx_upsw(:,:))) 
          IF (prpraf) THEN
            dR_spd_traf(1:ldc%nproma,1:klev+1,krow)  = dR_spd_traf(1:ldc%nproma,1:klev+1,krow) + &
                                                      dtprps*0.5_wp*((aflx_dnlw_clr(:,:) - aflx_uplw_clr(:,:)) - & 
                    (cflx_dnlw_clr(:,:) - cflx_uplw_clr(:,:)) + &
                    (dflx_dnlw_clr(:,:) - dflx_uplw_clr(:,:)) - &
                    (bflx_dnlw_clr(:,:) - bflx_uplw_clr(:,:))) 
            dR_spd_sraf(1:ldc%nproma,1:klev+1,krow)  = dR_spd_sraf(1:ldc%nproma,1:klev+1,krow) + &
                                                      dtprps*0.5_wp*((aflx_dnsw_clr(:,:) - aflx_upsw_clr(:,:)) - & 
                    (cflx_dnsw_clr(:,:) - cflx_upsw_clr(:,:)) + &
                    (dflx_dnsw_clr(:,:) - dflx_upsw_clr(:,:)) - &
                    (bflx_dnsw_clr(:,:) - bflx_upsw_clr(:,:)))
          END IF
        END IF

        dQ_spd_trad(1:ldc%nproma,1:klev,krow)  = dQ_spd_trad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                  ((((aflx_dnlw(:,1:klev) - aflx_uplw(:,1:klev)) - &
                                                    (aflx_dnlw(:,2:klev+1) - aflx_uplw(:,2:klev+1))) - &
                                                    ((cflx_dnlw(:,1:klev) - cflx_uplw(:,1:klev)) - & 
                                                    (cflx_dnlw(:,2:klev+1) - cflx_uplw(:,2:klev+1)))) - &
                                                  (((bflx_dnlw(:,1:klev) - bflx_uplw(:,1:klev)) - & 
                                                    (bflx_dnlw(:,2:klev+1) - bflx_uplw(:,2:klev+1))) - &
                                                    ((dflx_dnlw(:,1:klev) - dflx_uplw(:,1:klev)) - &
                                                    (dflx_dnlw(:,2:klev+1) - dflx_uplw(:,2:klev+1)))))
        dQ_spd_srad(1:ldc%nproma,1:klev,krow)  = dQ_spd_srad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                  ((((aflx_dnsw(:,1:klev) - aflx_upsw(:,1:klev)) - &
                                                    (aflx_dnsw(:,2:klev+1) - aflx_upsw(:,2:klev+1))) - &
                                                    ((cflx_dnsw(:,1:klev) - cflx_upsw(:,1:klev)) - & 
                                                    (cflx_dnsw(:,2:klev+1) - cflx_upsw(:,2:klev+1)))) - &
                                                  (((bflx_dnsw(:,1:klev) - bflx_upsw(:,1:klev)) - & 
                                                    (bflx_dnsw(:,2:klev+1) - bflx_upsw(:,2:klev+1))) - &
                                                    ((dflx_dnsw(:,1:klev) - dflx_upsw(:,1:klev)) - &
                                                    (dflx_dnsw(:,2:klev+1) - dflx_upsw(:,2:klev+1)))))
        IF (prpraf) THEN
        dQ_spd_traf(1:ldc%nproma,1:klev,krow)  = dQ_spd_traf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                  ((((aflx_dnlw_clr(:,1:klev) - aflx_uplw_clr(:,1:klev)) - &
                                                    (aflx_dnlw_clr(:,2:klev+1) - aflx_uplw_clr(:,2:klev+1))) - &
                                                    ((cflx_dnlw_clr(:,1:klev) - cflx_uplw_clr(:,1:klev)) - & 
                                                    (cflx_dnlw_clr(:,2:klev+1) - cflx_uplw_clr(:,2:klev+1)))) - &
                                                  (((bflx_dnlw_clr(:,1:klev) - bflx_uplw_clr(:,1:klev)) - & 
                                                    (bflx_dnlw_clr(:,2:klev+1) - bflx_uplw_clr(:,2:klev+1))) - &
                                                    ((dflx_dnlw_clr(:,1:klev) - dflx_uplw_clr(:,1:klev)) - &
                                                    (dflx_dnlw_clr(:,2:klev+1) - dflx_uplw_clr(:,2:klev+1)))))
        dQ_spd_sraf(1:ldc%nproma,1:klev,krow)  = dQ_spd_sraf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                  ((((aflx_dnsw_clr(:,1:klev) - aflx_upsw_clr(:,1:klev)) - &
                                                    (aflx_dnsw_clr(:,2:klev+1) - aflx_upsw_clr(:,2:klev+1))) - &
                                                    ((cflx_dnsw_clr(:,1:klev) - cflx_upsw_clr(:,1:klev)) - & 
                                                    (cflx_dnsw_clr(:,2:klev+1) - cflx_upsw_clr(:,2:klev+1)))) - &
                                                  (((bflx_dnsw_clr(:,1:klev) - bflx_upsw_clr(:,1:klev)) - & 
                                                    (bflx_dnsw_clr(:,2:klev+1) - bflx_upsw_clr(:,2:klev+1))) - &
                                                    ((dflx_dnsw_clr(:,1:klev) - dflx_upsw_clr(:,1:klev)) - &
                                                    (dflx_dnsw_clr(:,2:klev+1) - dflx_upsw_clr(:,2:klev+1)))))
        END IF
      END IF

      SELECT CASE (iaero(1))
      CASE (0)
        iaero_direct(1) = 0
        iaero_indirect(1) = 0
      CASE (8)
        iaero_direct(1) = 12              ! Only direct-effect of SP in state c, but if volcanoes -> added to aod
        iaero_indirect(1) = 11            ! Only indirect-effect of SP in state d (no volcano)
      CASE (9)
        iaero_direct(1) = 10              ! If experiment with SP only
        iaero_indirect(1) = 11
      CASE DEFAULT
      END SELECT

      CALL psrad_interface( &
            & iaero_direct    ,kproma          ,kbdim            ,klev             ,& 
            & krow            ,ktrac           ,ktype            ,nb_sw            ,& 
            & loland          ,loglac          ,cemiss           ,cos_mu0          ,& 
            & pgeom1          ,alb_vis_dir     ,alb_nir_dir      ,alb_vis_dif      ,&
            & alb_nir_dif     ,pp_fl           ,pp_hl            ,pp_sfc           ,&
            & tk_fl           ,tk_hl           ,tk_sfc           ,xq_vap           ,&
            & xq_liq          ,xq_ice          ,zcdnc            ,cld_frc          ,&
            & cld_cvr         ,xm_o3           ,xm_co2           ,xm_ch4           ,&
            & xm_n2o          ,xm_cfc          ,xm_o2            ,pxtm1            ,&
            & cflx_uplw       ,cflx_uplw_clr   ,cflx_dnlw        ,cflx_dnlw_clr    ,&
            & cflx_upsw       ,cflx_upsw_clr   ,cflx_dnsw        ,cflx_dnsw_clr    ,&
            & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
            & par_dff_sfc                                                          )
     ! CALL message('','indirect-effect with saved state')
      CALL psrad_interface( &
            & iaero_indirect  ,kproma          ,kbdim            ,klev             ,& 
            & krow            ,ktrac           ,zktype           ,nb_sw            ,&
            & zloland         ,zloglac         ,zcemiss          ,zcos_mu0         ,& 
            & zpgeom1         ,zalb_vis_dir    ,zalb_nir_dir     ,zalb_vis_dif     ,&
            & zalb_nir_dif    ,zpp_fl          ,zpp_hl           ,zpp_sfc          ,& 
            & ztk_fl          ,ztk_hl          ,ztk_sfc          ,zq_vap           ,&
            & zq_liq          ,zq_ice          ,cdnc             ,zcld_frc         ,&
            & zcld_cvr        ,zm_o3           ,zm_co2           ,zm_ch4           ,&
            & zm_n2o          ,zm_cfc          ,zm_o2            ,pxtm1            ,&
            & dflx_uplw       ,dflx_uplw_clr   ,dflx_dnlw        ,dflx_dnlw_clr    ,&
            & dflx_upsw       ,dflx_upsw_clr   ,dflx_dnsw        ,dflx_dnsw_clr    ,&
            & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
            & par_dff_sfc                                                          )

        dR_spi_trad0(1:ldc%nproma,krow)  = dR_spi_trad0(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw(:,1) - aflx_uplw(:,1)) - & 
                (cflx_dnlw(:,1) - cflx_uplw(:,1)) + &
                                                          (dflx_dnlw(:,1) - dflx_uplw(:,1)) - &
                                                          (bflx_dnlw(:,1) - bflx_uplw(:,1))) 
        dR_spi_srad0(1:ldc%nproma,krow)  = dR_spi_srad0(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw(:,1) - aflx_upsw(:,1)) - & 
                (cflx_dnsw(:,1) - cflx_upsw(:,1)) + &
                                                          (dflx_dnsw(:,1) - dflx_upsw(:,1)) - &
                                                          (bflx_dnsw(:,1) - bflx_upsw(:,1))) 
      IF (prprads) THEN
        dR_spi_trads(1:ldc%nproma,krow)  = dR_spi_trads(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw(:,klev+1) - aflx_uplw(:,klev+1)) - & 
                (cflx_dnlw(:,klev+1) - cflx_uplw(:,klev+1)) + &
                                                          (dflx_dnlw(:,klev+1) - dflx_uplw(:,klev+1)) - &
                                                          (bflx_dnlw(:,klev+1) - bflx_uplw(:,klev+1))) 
        dR_spi_srads(1:ldc%nproma,krow)  = dR_spi_srads(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw(:,klev+1) - aflx_upsw(:,klev+1)) - & 
                (cflx_dnsw(:,klev+1) - cflx_upsw(:,klev+1)) + &
                                                          (dflx_dnsw(:,klev+1) - dflx_upsw(:,klev+1)) - &
                                                          (bflx_dnsw(:,klev+1) - bflx_upsw(:,klev+1))) 
      END IF
      IF (prpraf) THEN
        dR_spi_traf0(1:ldc%nproma,krow)  = dR_spi_traf0(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw_clr(:,1) - aflx_uplw_clr(:,1)) - & 
                (cflx_dnlw_clr(:,1) - cflx_uplw_clr(:,1)) + &
                                                          (dflx_dnlw_clr(:,1) - dflx_uplw_clr(:,1)) - &
                                                          (bflx_dnlw_clr(:,1) - bflx_uplw_clr(:,1))) 
        dR_spi_sraf0(1:ldc%nproma,krow)  = dR_spi_sraf0(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw_clr(:,1) - aflx_upsw_clr(:,1)) - & 
                (cflx_dnsw_clr(:,1) - cflx_upsw_clr(:,1)) + &
                                                          (dflx_dnsw_clr(:,1) - dflx_upsw_clr(:,1)) - &
                                                          (bflx_dnsw_clr(:,1) - bflx_upsw_clr(:,1))) 
        dR_spi_trafs(1:ldc%nproma,krow)  = dR_spi_trafs(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw_clr(:,klev+1) - aflx_uplw_clr(:,klev+1)) - & 
                (cflx_dnlw_clr(:,klev+1) - cflx_uplw_clr(:,klev+1)) + &
                                                          (dflx_dnlw_clr(:,klev+1) - dflx_uplw_clr(:,klev+1)) - &
                                                          (bflx_dnlw_clr(:,klev+1) - bflx_uplw_clr(:,klev+1))) 
        dR_spi_srafs(1:ldc%nproma,krow)  = dR_spi_srafs(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw_clr(:,klev+1) - aflx_upsw_clr(:,klev+1)) - & 
                (cflx_dnsw_clr(:,klev+1) - cflx_upsw_clr(:,klev+1)) + &
                                                          (dflx_dnsw_clr(:,klev+1) - dflx_upsw_clr(:,klev+1)) - &
                                                          (bflx_dnsw_clr(:,klev+1) - bflx_upsw_clr(:,klev+1))) 
      END IF
      IF (prp3doutput) THEN
        IF (prp3drad) THEN
        dR_spi_trad(1:ldc%nproma,1:klev+1,krow)  = dR_spi_trad(1:ldc%nproma,1:klev+1,krow) + &
                                                    dtprps*0.5_wp*((aflx_dnlw(:,:) - aflx_uplw(:,:)) - & 
                  (cflx_dnlw(:,:) - cflx_uplw(:,:)) + &
                  (dflx_dnlw(:,:) - dflx_uplw(:,:)) - &
                  (bflx_dnlw(:,:) - bflx_uplw(:,:))) 
        dR_spi_srad(1:ldc%nproma,1:klev+1,krow)  = dR_spi_srad(1:ldc%nproma,1:klev+1,krow) + &
                                                    dtprps*0.5_wp*((aflx_dnsw(:,:) - aflx_upsw(:,:)) - & 
                  (cflx_dnsw(:,:) - cflx_upsw(:,:)) + &
                  (dflx_dnsw(:,:) - dflx_upsw(:,:)) - &
                  (bflx_dnsw(:,:) - bflx_upsw(:,:))) 
          IF (prpraf) THEN
            dR_spi_traf(1:ldc%nproma,1:klev+1,krow)  = dR_spi_traf(1:ldc%nproma,1:klev+1,krow) + &
                                                      dtprps*0.5_wp*((aflx_dnlw_clr(:,:) - aflx_uplw_clr(:,:)) - & 
                                              (cflx_dnlw_clr(:,:) - cflx_uplw_clr(:,:)) + &
                                              (dflx_dnlw_clr(:,:) - dflx_uplw_clr(:,:)) - &
                                              (bflx_dnlw_clr(:,:) - bflx_uplw_clr(:,:))) 
            dR_spi_sraf(1:ldc%nproma,1:klev+1,krow)  = dR_spi_sraf(1:ldc%nproma,1:klev+1,krow) + &
                                                      dtprps*0.5_wp*((aflx_dnsw_clr(:,:) - aflx_upsw_clr(:,:)) - & 
                                              (cflx_dnsw_clr(:,:) - cflx_upsw_clr(:,:)) + &
                                              (dflx_dnsw_clr(:,:) - dflx_upsw_clr(:,:)) - &
                                              (bflx_dnsw_clr(:,:) - bflx_upsw_clr(:,:)))
          END IF
        END IF
        dQ_spi_trad(1:ldc%nproma,1:klev,krow)  = dQ_spi_trad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                  ((((aflx_dnlw(:,1:klev) - aflx_uplw(:,1:klev)) - &
                                                    (aflx_dnlw(:,2:klev+1) - aflx_uplw(:,2:klev+1))) - &
                                                    ((cflx_dnlw(:,1:klev) - cflx_uplw(:,1:klev)) - & 
                                                    (cflx_dnlw(:,2:klev+1) - cflx_uplw(:,2:klev+1)))) - &
                                                  (((bflx_dnlw(:,1:klev) - bflx_uplw(:,1:klev)) - & 
                                                    (bflx_dnlw(:,2:klev+1) - bflx_uplw(:,2:klev+1))) - &
                                                    ((dflx_dnlw(:,1:klev) - dflx_uplw(:,1:klev)) - &
                                                    (dflx_dnlw(:,2:klev+1) - dflx_uplw(:,2:klev+1)))))
        dQ_spi_srad(1:ldc%nproma,1:klev,krow)  = dQ_spi_srad(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                  ((((aflx_dnsw(:,1:klev) - aflx_upsw(:,1:klev)) - &
                                                    (aflx_dnsw(:,2:klev+1) - aflx_upsw(:,2:klev+1))) - &
                                                    ((cflx_dnsw(:,1:klev) - cflx_upsw(:,1:klev)) - & 
                                                    (cflx_dnsw(:,2:klev+1) - cflx_upsw(:,2:klev+1)))) - &
                                                  (((bflx_dnsw(:,1:klev) - bflx_upsw(:,1:klev)) - & 
                                                    (bflx_dnsw(:,2:klev+1) - bflx_upsw(:,2:klev+1))) - &
                                                    ((dflx_dnsw(:,1:klev) - dflx_upsw(:,1:klev)) - &
                                                    (dflx_dnsw(:,2:klev+1) - dflx_upsw(:,2:klev+1)))))
        IF (prpraf) THEN
        dQ_spi_traf(1:ldc%nproma,1:klev,krow)  = dQ_spi_traf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                  ((((aflx_dnlw_clr(:,1:klev) - aflx_uplw_clr(:,1:klev)) - &
                                                    (aflx_dnlw_clr(:,2:klev+1) - aflx_uplw_clr(:,2:klev+1))) - &
                                                    ((cflx_dnlw_clr(:,1:klev) - cflx_uplw_clr(:,1:klev)) - & 
                                                    (cflx_dnlw_clr(:,2:klev+1) - cflx_uplw_clr(:,2:klev+1)))) - &
                                                  (((bflx_dnlw_clr(:,1:klev) - bflx_uplw_clr(:,1:klev)) - & 
                                                    (bflx_dnlw_clr(:,2:klev+1) - bflx_uplw_clr(:,2:klev+1))) - &
                                                    ((dflx_dnlw_clr(:,1:klev) - dflx_uplw_clr(:,1:klev)) - &
                                                    (dflx_dnlw_clr(:,2:klev+1) - dflx_uplw_clr(:,2:klev+1)))))
        dQ_spi_sraf(1:ldc%nproma,1:klev,krow)  = dQ_spi_sraf(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                  ((((aflx_dnsw_clr(:,1:klev) - aflx_upsw_clr(:,1:klev)) - &
                                                    (aflx_dnsw_clr(:,2:klev+1) - aflx_upsw_clr(:,2:klev+1))) - &
                                                    ((cflx_dnsw_clr(:,1:klev) - cflx_upsw_clr(:,1:klev)) - & 
                                                    (cflx_dnsw_clr(:,2:klev+1) - cflx_upsw_clr(:,2:klev+1)))) - &
                                                  (((bflx_dnsw_clr(:,1:klev) - bflx_upsw_clr(:,1:klev)) - & 
                                                    (bflx_dnsw_clr(:,2:klev+1) - bflx_upsw_clr(:,2:klev+1))) - &
                                                    ((dflx_dnsw_clr(:,1:klev) - dflx_upsw_clr(:,1:klev)) - &
                                                    (dflx_dnsw_clr(:,2:klev+1) - dflx_upsw_clr(:,2:klev+1)))))
        END IF
      END IF
    
    ! single plumes. Only TOA calculated (most relevant)

    IF (prp_single_sp) THEN

      DO iplume=1, nplumes

        ! Calculates total forcing of individual plumes:
        iaero_sp(2)=iplume           ! single plume only in saved state
        xiaero(2)=10+iplume           ! all plumes but single plume in current state

       ! CALL message('','total single plume in current state')
        CALL psrad_interface( &  ! iaero includes volcanoes in historical simulation
            & xiaero          ,kproma          ,kbdim            ,klev             ,& 
            & krow            ,ktrac           ,ktype            ,nb_sw            ,& 
            & loland          ,loglac          ,cemiss           ,cos_mu0          ,& 
            & pgeom1          ,alb_vis_dir     ,alb_nir_dir      ,alb_vis_dif      ,&
            & alb_nir_dif     ,pp_fl           ,pp_hl            ,pp_sfc           ,&
            & tk_fl           ,tk_hl           ,tk_sfc           ,xq_vap           ,&
            & xq_liq          ,xq_ice          ,cdnc             ,cld_frc          ,&
            & cld_cvr         ,xm_o3           ,xm_co2           ,xm_ch4           ,&
            & xm_n2o          ,xm_cfc          ,xm_o2            ,pxtm1            ,&
            & cflx_uplw       ,cflx_uplw_clr   ,cflx_dnlw        ,cflx_dnlw_clr    ,&
            & cflx_upsw       ,cflx_upsw_clr   ,cflx_dnsw        ,cflx_dnsw_clr    ,&
            & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
            & par_dff_sfc                                                          )
        
       ! CALL message('','total single plume with saved state')
        CALL psrad_interface( &
            & iaero_sp        ,kproma          ,kbdim            ,klev             ,& 
            & krow            ,ktrac           ,zktype           ,nb_sw            ,&
            & zloland         ,zloglac         ,zcemiss          ,zcos_mu0         ,& 
            & zpgeom1         ,zalb_vis_dir    ,zalb_nir_dir     ,zalb_vis_dif     ,&
            & zalb_nir_dif    ,zpp_fl          ,zpp_hl           ,zpp_sfc          ,& 
            & ztk_fl          ,ztk_hl          ,ztk_sfc          ,zq_vap           ,&
            & zq_liq          ,zq_ice          ,cdnc             ,zcld_frc         ,&
            & zcld_cvr        ,zm_o3           ,zm_co2           ,zm_ch4           ,&
            & zm_n2o          ,zm_cfc          ,zm_o2            ,pxtm1            ,&
            & dflx_uplw       ,dflx_uplw_clr   ,dflx_dnlw        ,dflx_dnlw_clr    ,&
            & dflx_upsw       ,dflx_upsw_clr   ,dflx_dnsw        ,dflx_dnsw_clr    ,&
            & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
            & par_dff_sfc                                                          )

        dR_single_sp_trad0(iplume)%p(1:ldc%nproma,krow)  = dR_single_sp_trad0(iplume)%p(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw(:,1) - aflx_uplw(:,1)) - & 
                (cflx_dnlw(:,1) - cflx_uplw(:,1)) + &
                                                          (dflx_dnlw(:,1) - dflx_uplw(:,1)) - &
                                                          (bflx_dnlw(:,1) - bflx_uplw(:,1))) 
        dR_single_sp_srad0(iplume)%p(1:ldc%nproma,krow)  = dR_single_sp_srad0(iplume)%p(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw(:,1) - aflx_upsw(:,1)) - & 
                (cflx_dnsw(:,1) - cflx_upsw(:,1)) + &
                                                          (dflx_dnsw(:,1) - dflx_upsw(:,1)) - &
                                                          (bflx_dnsw(:,1) - bflx_upsw(:,1))) 
        IF (prpraf) THEN
          dR_single_sp_traf0(iplume)%p(1:ldc%nproma,krow)  = dR_single_sp_traf0(iplume)%p(1:ldc%nproma,krow) + &
                                              dtprps*0.5_wp*((aflx_dnlw_clr(:,1) - aflx_uplw_clr(:,1)) - & 
                  (cflx_dnlw_clr(:,1) - cflx_uplw_clr(:,1)) + &
                                                            (dflx_dnlw_clr(:,1) - dflx_uplw_clr(:,1)) - &
                                                            (bflx_dnlw_clr(:,1) - bflx_uplw_clr(:,1))) 
          dR_single_sp_sraf0(iplume)%p(1:ldc%nproma,krow)  = dR_single_sp_sraf0(iplume)%p(1:ldc%nproma,krow) + &
                                              dtprps*0.5_wp*((aflx_dnsw_clr(:,1) - aflx_upsw_clr(:,1)) - & 
                  (cflx_dnsw_clr(:,1) - cflx_upsw_clr(:,1)) + &
                                                            (dflx_dnsw_clr(:,1) - dflx_upsw_clr(:,1)) - &
                                                            (bflx_dnsw_clr(:,1) - bflx_upsw_clr(:,1))) 
        END IF
        IF (prp3doutput) THEN
          IF (prp3drad) THEN
          dR_single_sp_trad(iplume)%p(1:ldc%nproma,1:klev+1,krow)  = dR_single_sp_trad(iplume)%p(1:ldc%nproma,1:klev+1,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw(:,:) - aflx_uplw(:,:)) - & 
                                  (cflx_dnlw(:,:) - cflx_uplw(:,:)) + &
                                  (dflx_dnlw(:,:) - dflx_uplw(:,:)) - &
                                  (bflx_dnlw(:,:) - bflx_uplw(:,:))) 
          dR_single_sp_srad(iplume)%p(1:ldc%nproma,1:klev+1,krow)  = dR_single_sp_srad(iplume)%p(1:ldc%nproma,1:klev+1,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw(:,:) - aflx_upsw(:,:)) - & 
                                  (cflx_dnsw(:,:) - cflx_upsw(:,:)) + &
                                  (dflx_dnsw(:,:) - dflx_upsw(:,:)) - &
                                  (bflx_dnsw(:,:) - bflx_upsw(:,:))) 
            IF (prpraf) THEN
              dR_single_sp_traf(iplume)%p(1:ldc%nproma,1:klev+1,krow)  = dR_single_sp_traf(iplume)%p(1:ldc%nproma,1:klev+1,krow) + &
                                                                dtprps*0.5_wp*((aflx_dnlw_clr(:,:) - aflx_uplw_clr(:,:)) - & 
                                                      (cflx_dnlw_clr(:,:) - cflx_uplw_clr(:,:)) + &
                                                      (dflx_dnlw_clr(:,:) - dflx_uplw_clr(:,:)) - &
                                                      (bflx_dnlw_clr(:,:) - bflx_uplw_clr(:,:))) 
              dR_single_sp_sraf(iplume)%p(1:ldc%nproma,1:klev+1,krow)  = dR_single_sp_sraf(iplume)%p(1:ldc%nproma,1:klev+1,krow) + &
                                                                dtprps*0.5_wp*((aflx_dnsw_clr(:,:) - aflx_upsw_clr(:,:)) - & 
                                                      (cflx_dnsw_clr(:,:) - cflx_upsw_clr(:,:)) + &
                                                      (dflx_dnsw_clr(:,:) - dflx_upsw_clr(:,:)) - &
                                                      (bflx_dnsw_clr(:,:) - bflx_upsw_clr(:,:)))
            END IF
          END IF 
          dQ_single_sp_trad(iplume)%p(1:ldc%nproma,1:klev,krow)  = dQ_single_sp_trad(iplume)%p(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                    ((((aflx_dnlw(:,1:klev) - aflx_uplw(:,1:klev)) - &
                                                      (aflx_dnlw(:,2:klev+1) - aflx_uplw(:,2:klev+1))) - &
                                                      ((cflx_dnlw(:,1:klev) - cflx_uplw(:,1:klev)) - & 
                                                      (cflx_dnlw(:,2:klev+1) - cflx_uplw(:,2:klev+1)))) - &
                                                    (((bflx_dnlw(:,1:klev) - bflx_uplw(:,1:klev)) - & 
                                                      (bflx_dnlw(:,2:klev+1) - bflx_uplw(:,2:klev+1))) - &
                                                      ((dflx_dnlw(:,1:klev) - dflx_uplw(:,1:klev)) - &
                                                      (dflx_dnlw(:,2:klev+1) - dflx_uplw(:,2:klev+1)))))
          dQ_single_sp_srad(iplume)%p(1:ldc%nproma,1:klev,krow)  = dQ_single_sp_srad(iplume)%p(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                    ((((aflx_dnsw(:,1:klev) - aflx_upsw(:,1:klev)) - &
                                                      (aflx_dnsw(:,2:klev+1) - aflx_upsw(:,2:klev+1))) - &
                                                      ((cflx_dnsw(:,1:klev) - cflx_upsw(:,1:klev)) - & 
                                                      (cflx_dnsw(:,2:klev+1) - cflx_upsw(:,2:klev+1)))) - &
                                                    (((bflx_dnsw(:,1:klev) - bflx_upsw(:,1:klev)) - & 
                                                      (bflx_dnsw(:,2:klev+1) - bflx_upsw(:,2:klev+1))) - &
                                                      ((dflx_dnsw(:,1:klev) - dflx_upsw(:,1:klev)) - &
                                                      (dflx_dnsw(:,2:klev+1) - dflx_upsw(:,2:klev+1)))))
          IF (prpraf) THEN  
          dQ_single_sp_traf(iplume)%p(1:ldc%nproma,1:klev,krow)  = dQ_single_sp_traf(iplume)%p(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                    ((((aflx_dnlw_clr(:,1:klev) - aflx_uplw_clr(:,1:klev)) - &
                                                      (aflx_dnlw_clr(:,2:klev+1) - aflx_uplw_clr(:,2:klev+1))) - &
                                                      ((cflx_dnlw_clr(:,1:klev) - cflx_uplw_clr(:,1:klev)) - & 
                                                      (cflx_dnlw_clr(:,2:klev+1) - cflx_uplw_clr(:,2:klev+1)))) - &
                                                    (((bflx_dnlw_clr(:,1:klev) - bflx_uplw_clr(:,1:klev)) - & 
                                                      (bflx_dnlw_clr(:,2:klev+1) - bflx_uplw_clr(:,2:klev+1))) - &
                                                      ((dflx_dnlw_clr(:,1:klev) - dflx_uplw_clr(:,1:klev)) - &
                                                      (dflx_dnlw_clr(:,2:klev+1) - dflx_uplw_clr(:,2:klev+1)))))
          dQ_single_sp_sraf(iplume)%p(1:ldc%nproma,1:klev,krow)  = dQ_single_sp_sraf(iplume)%p(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                    ((((aflx_dnsw_clr(:,1:klev) - aflx_upsw_clr(:,1:klev)) - &
                                                      (aflx_dnsw_clr(:,2:klev+1) - aflx_upsw_clr(:,2:klev+1))) - &
                                                      ((cflx_dnsw_clr(:,1:klev) - cflx_upsw_clr(:,1:klev)) - & 
                                                      (cflx_dnsw_clr(:,2:klev+1) - cflx_upsw_clr(:,2:klev+1)))) - &
                                                    (((bflx_dnsw_clr(:,1:klev) - bflx_upsw_clr(:,1:klev)) - & 
                                                      (bflx_dnsw_clr(:,2:klev+1) - bflx_upsw_clr(:,2:klev+1))) - &
                                                      ((dflx_dnsw_clr(:,1:klev) - dflx_upsw_clr(:,1:klev)) - &
                                                      (dflx_dnsw_clr(:,2:klev+1) - dflx_upsw_clr(:,2:klev+1)))))
          END IF
        END IF

      ! single plumes direct effect. Only TOA calculated (most relevant)
        SELECT CASE (iaero(1))
        CASE (0)
          iaero_direct(1) = 0
          iaero_indirect(1) = 0
        CASE (8)
          iaero_direct(1) = 10              ! Only direct-effect of SP in state d
          iaero_indirect(1) = 13            ! Only indirect-effect of SP in state c, but if volcanoes -> added to aod
        CASE (9)
          iaero_direct(1) = 10              ! If experiment with SP only
          iaero_indirect(1) = 11
        CASE DEFAULT
        END SELECT
        iaero_direct(2)=iplume
        iaero_indirect(2)=10+iplume

       ! CALL message('','direct-effect single plume with current state')
        CALL psrad_interface( &
            & iaero_indirect  ,kproma          ,kbdim            ,klev             ,& 
            & krow            ,ktrac           ,ktype            ,nb_sw            ,& 
            & loland          ,loglac          ,cemiss           ,cos_mu0          ,& 
            & pgeom1          ,alb_vis_dir     ,alb_nir_dir      ,alb_vis_dif      ,&
            & alb_nir_dif     ,pp_fl           ,pp_hl            ,pp_sfc           ,&
            & tk_fl           ,tk_hl           ,tk_sfc           ,xq_vap           ,&
            & xq_liq          ,xq_ice          ,cdnc             ,cld_frc          ,&
            & cld_cvr         ,xm_o3           ,xm_co2           ,xm_ch4           ,&
            & xm_n2o          ,xm_cfc          ,xm_o2            ,pxtm1            ,&
            & cflx_uplw       ,cflx_uplw_clr   ,cflx_dnlw        ,cflx_dnlw_clr    ,&
            & cflx_upsw       ,cflx_upsw_clr   ,cflx_dnsw        ,cflx_dnsw_clr    ,&
            & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
            & par_dff_sfc                                                          )
       ! CALL message('','direct-effect single plume with saved state')
        CALL psrad_interface( &
            & iaero_direct    ,kproma          ,kbdim            ,klev             ,& 
            & krow            ,ktrac           ,zktype           ,nb_sw            ,&
            & zloland         ,zloglac         ,zcemiss          ,zcos_mu0         ,& 
            & zpgeom1         ,zalb_vis_dir    ,zalb_nir_dir     ,zalb_vis_dif     ,&
            & zalb_nir_dif    ,zpp_fl          ,zpp_hl           ,zpp_sfc          ,& 
            & ztk_fl          ,ztk_hl          ,ztk_sfc          ,zq_vap           ,&
            & zq_liq          ,zq_ice          ,zcdnc            ,zcld_frc         ,&
            & zcld_cvr        ,zm_o3           ,zm_co2           ,zm_ch4           ,&
            & zm_n2o          ,zm_cfc          ,zm_o2            ,pxtm1            ,&
            & dflx_uplw       ,dflx_uplw_clr   ,dflx_dnlw        ,dflx_dnlw_clr    ,&
            & dflx_upsw       ,dflx_upsw_clr   ,dflx_dnsw        ,dflx_dnsw_clr    ,&
            & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
            & par_dff_sfc                                                          )

        dR_single_spd_trad0(iplume)%p(1:ldc%nproma,krow)  = dR_single_spd_trad0(iplume)%p(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw(:,1) - aflx_uplw(:,1)) - & 
                (cflx_dnlw(:,1) - cflx_uplw(:,1)) + &
                                                          (dflx_dnlw(:,1) - dflx_uplw(:,1)) - &
                                                          (bflx_dnlw(:,1) - bflx_uplw(:,1))) 
        dR_single_spd_srad0(iplume)%p(1:ldc%nproma,krow)  = dR_single_spd_srad0(iplume)%p(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw(:,1) - aflx_upsw(:,1)) - & 
                (cflx_dnsw(:,1) - cflx_upsw(:,1)) + &
                                                          (dflx_dnsw(:,1) - dflx_upsw(:,1)) - &
                                                          (bflx_dnsw(:,1) - bflx_upsw(:,1))) 
        IF (prpraf) THEN
          dR_single_spd_traf0(iplume)%p(1:ldc%nproma,krow)  = dR_single_spd_traf0(iplume)%p(1:ldc%nproma,krow) + &
                                              dtprps*0.5_wp*((aflx_dnlw_clr(:,1) - aflx_uplw_clr(:,1)) - & 
                  (cflx_dnlw_clr(:,1) - cflx_uplw_clr(:,1)) + &
                                                            (dflx_dnlw_clr(:,1) - dflx_uplw_clr(:,1)) - &
                                                            (bflx_dnlw_clr(:,1) - bflx_uplw_clr(:,1))) 
          dR_single_spd_sraf0(iplume)%p(1:ldc%nproma,krow)  = dR_single_spd_sraf0(iplume)%p(1:ldc%nproma,krow) + &
                                              dtprps*0.5_wp*((aflx_dnsw_clr(:,1) - aflx_upsw_clr(:,1)) - & 
                  (cflx_dnsw_clr(:,1) - cflx_upsw_clr(:,1)) + &
                                                            (dflx_dnsw_clr(:,1) - dflx_upsw_clr(:,1)) - &
                                                            (bflx_dnsw_clr(:,1) - bflx_upsw_clr(:,1))) 
        END IF
        IF (prp3doutput) THEN
          IF (prp3drad) THEN
          dR_single_spd_trad(iplume)%p(1:ldc%nproma,1:klev+1,krow)  = dR_single_spd_trad(iplume)%p(1:ldc%nproma,1:klev+1,krow) + &
                                                      dtprps*0.5_wp*((aflx_dnlw(:,:) - aflx_uplw(:,:)) - & 
                                            (cflx_dnlw(:,:) - cflx_uplw(:,:)) + &
                                            (dflx_dnlw(:,:) - dflx_uplw(:,:)) - &
                                            (bflx_dnlw(:,:) - bflx_uplw(:,:))) 
          dR_single_spd_srad(iplume)%p(1:ldc%nproma,1:klev+1,krow)  = dR_single_spd_srad(iplume)%p(1:ldc%nproma,1:klev+1,krow) + &
                                                      dtprps*0.5_wp*((aflx_dnsw(:,:) - aflx_upsw(:,:)) - & 
                                            (cflx_dnsw(:,:) - cflx_upsw(:,:)) + &
                                            (dflx_dnsw(:,:) - dflx_upsw(:,:)) - &
                                            (bflx_dnsw(:,:) - bflx_upsw(:,:))) 
            IF (prpraf) THEN
              dR_single_spd_traf(iplume)%p(1:ldc%nproma,1:klev+1,krow)  = dR_single_spd_traf(iplume)%p(1:ldc%nproma,1:klev+1,krow) + &
                                                                  dtprps*0.5_wp*((aflx_dnlw_clr(:,:) - aflx_uplw_clr(:,:)) - & 
                                                        (cflx_dnlw_clr(:,:) - cflx_uplw_clr(:,:)) + &
                                                        (dflx_dnlw_clr(:,:) - dflx_uplw_clr(:,:)) - &
                                                        (bflx_dnlw_clr(:,:) - bflx_uplw_clr(:,:))) 
              dR_single_spd_sraf(iplume)%p(1:ldc%nproma,1:klev+1,krow)  = dR_single_spd_sraf(iplume)%p(1:ldc%nproma,1:klev+1,krow) + &
                                                                  dtprps*0.5_wp*((aflx_dnsw_clr(:,:) - aflx_upsw_clr(:,:)) - & 
                                                        (cflx_dnsw_clr(:,:) - cflx_upsw_clr(:,:)) + &
                                                        (dflx_dnsw_clr(:,:) - dflx_upsw_clr(:,:)) - &
                                                        (bflx_dnsw_clr(:,:) - bflx_upsw_clr(:,:)))
            END IF
          END IF
          dQ_single_spd_trad(iplume)%p(1:ldc%nproma,1:klev,krow)  = dQ_single_spd_trad(iplume)%p(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                    ((((aflx_dnlw(:,1:klev) - aflx_uplw(:,1:klev)) - &
                                                      (aflx_dnlw(:,2:klev+1) - aflx_uplw(:,2:klev+1))) - &
                                                      ((cflx_dnlw(:,1:klev) - cflx_uplw(:,1:klev)) - & 
                                                      (cflx_dnlw(:,2:klev+1) - cflx_uplw(:,2:klev+1)))) - &
                                                    (((bflx_dnlw(:,1:klev) - bflx_uplw(:,1:klev)) - & 
                                                      (bflx_dnlw(:,2:klev+1) - bflx_uplw(:,2:klev+1))) - &
                                                      ((dflx_dnlw(:,1:klev) - dflx_uplw(:,1:klev)) - &
                                                      (dflx_dnlw(:,2:klev+1) - dflx_uplw(:,2:klev+1)))))
          dQ_single_spd_srad(iplume)%p(1:ldc%nproma,1:klev,krow)  = dQ_single_spd_srad(iplume)%p(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                    ((((aflx_dnsw(:,1:klev) - aflx_upsw(:,1:klev)) - &
                                                      (aflx_dnsw(:,2:klev+1) - aflx_upsw(:,2:klev+1))) - &
                                                      ((cflx_dnsw(:,1:klev) - cflx_upsw(:,1:klev)) - & 
                                                      (cflx_dnsw(:,2:klev+1) - cflx_upsw(:,2:klev+1)))) - &
                                                    (((bflx_dnsw(:,1:klev) - bflx_upsw(:,1:klev)) - & 
                                                      (bflx_dnsw(:,2:klev+1) - bflx_upsw(:,2:klev+1))) - &
                                                      ((dflx_dnsw(:,1:klev) - dflx_upsw(:,1:klev)) - &
                                                      (dflx_dnsw(:,2:klev+1) - dflx_upsw(:,2:klev+1)))))
          IF (prpraf) THEN
          dQ_single_spd_traf(iplume)%p(1:ldc%nproma,1:klev,krow)  = dQ_single_spd_traf(iplume)%p(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                    ((((aflx_dnlw_clr(:,1:klev) - aflx_uplw_clr(:,1:klev)) - &
                                                      (aflx_dnlw_clr(:,2:klev+1) - aflx_uplw_clr(:,2:klev+1))) - &
                                                      ((cflx_dnlw_clr(:,1:klev) - cflx_uplw_clr(:,1:klev)) - & 
                                                      (cflx_dnlw_clr(:,2:klev+1) - cflx_uplw_clr(:,2:klev+1)))) - &
                                                    (((bflx_dnlw_clr(:,1:klev) - bflx_uplw_clr(:,1:klev)) - & 
                                                      (bflx_dnlw_clr(:,2:klev+1) - bflx_uplw_clr(:,2:klev+1))) - &
                                                      ((dflx_dnlw_clr(:,1:klev) - dflx_uplw_clr(:,1:klev)) - &
                                                      (dflx_dnlw_clr(:,2:klev+1) - dflx_uplw_clr(:,2:klev+1)))))
          dQ_single_spd_sraf(iplume)%p(1:ldc%nproma,1:klev,krow)  = dQ_single_spd_sraf(iplume)%p(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                    ((((aflx_dnsw_clr(:,1:klev) - aflx_upsw_clr(:,1:klev)) - &
                                                      (aflx_dnsw_clr(:,2:klev+1) - aflx_upsw_clr(:,2:klev+1))) - &
                                                      ((cflx_dnsw_clr(:,1:klev) - cflx_upsw_clr(:,1:klev)) - & 
                                                      (cflx_dnsw_clr(:,2:klev+1) - cflx_upsw_clr(:,2:klev+1)))) - &
                                                    (((bflx_dnsw_clr(:,1:klev) - bflx_upsw_clr(:,1:klev)) - & 
                                                      (bflx_dnsw_clr(:,2:klev+1) - bflx_upsw_clr(:,2:klev+1))) - &
                                                      ((dflx_dnsw_clr(:,1:klev) - dflx_upsw_clr(:,1:klev)) - &
                                                      (dflx_dnsw_clr(:,2:klev+1) - dflx_upsw_clr(:,2:klev+1)))))
          END IF
        END IF

        SELECT CASE (iaero(1))
        CASE (0)
          iaero_direct(1) = 0
          iaero_indirect(1) = 0
        CASE (8)
          iaero_direct(1) = 12              ! Only direct-effect of SP in state c, but if volcanoes -> added to aod
          iaero_indirect(1) = 11            ! Only indirect-effect of SP in state d (no volcano)
        CASE (9)
          iaero_direct(1) = 10              ! If experiment with SP only
          iaero_indirect(1) = 11
        CASE DEFAULT
        END SELECT
        iaero_direct(2)=10+iplume
        iaero_indirect(2)=iplume

        ! Calculates indirect-effect forcing:
       ! CALL message('','indirect-effect single plume with current state')
        CALL psrad_interface( &
              & iaero_direct    ,kproma          ,kbdim            ,klev             ,& 
              & krow            ,ktrac           ,ktype            ,nb_sw            ,& 
              & loland          ,loglac          ,cemiss           ,cos_mu0          ,& 
              & pgeom1          ,alb_vis_dir     ,alb_nir_dir      ,alb_vis_dif      ,&
              & alb_nir_dif     ,pp_fl           ,pp_hl            ,pp_sfc           ,&
              & tk_fl           ,tk_hl           ,tk_sfc           ,xq_vap           ,&
              & xq_liq          ,xq_ice          ,cdnc             ,cld_frc          ,&
              & cld_cvr         ,xm_o3           ,xm_co2           ,xm_ch4           ,&
              & xm_n2o          ,xm_cfc          ,xm_o2            ,pxtm1            ,&
              & cflx_uplw       ,cflx_uplw_clr   ,cflx_dnlw        ,cflx_dnlw_clr    ,&
              & cflx_upsw       ,cflx_upsw_clr   ,cflx_dnsw        ,cflx_dnsw_clr    ,&
              & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
              & par_dff_sfc                                                          )
       ! CALL message('','indirect-effect single plume with saved state')
        CALL psrad_interface( &
              & iaero_indirect  ,kproma          ,kbdim            ,klev             ,& 
              & krow            ,ktrac           ,zktype           ,nb_sw            ,&
              & zloland         ,zloglac         ,zcemiss          ,zcos_mu0         ,& 
              & zpgeom1         ,zalb_vis_dir    ,zalb_nir_dir     ,zalb_vis_dif     ,&
              & zalb_nir_dif    ,zpp_fl          ,zpp_hl           ,zpp_sfc          ,& 
              & ztk_fl          ,ztk_hl          ,ztk_sfc          ,zq_vap           ,&
              & zq_liq          ,zq_ice          ,cdnc             ,zcld_frc         ,&
              & zcld_cvr        ,zm_o3           ,zm_co2           ,zm_ch4           ,&
              & zm_n2o          ,zm_cfc          ,zm_o2            ,pxtm1            ,&
              & dflx_uplw       ,dflx_uplw_clr   ,dflx_dnlw        ,dflx_dnlw_clr    ,&
              & dflx_upsw       ,dflx_upsw_clr   ,dflx_dnsw        ,dflx_dnsw_clr    ,&
              & vis_dff_sfc     ,dpar_sfc        ,nir_dff_sfc      ,dvis_sfc         ,&
              & par_dff_sfc                                                          )
  
        dR_single_spi_trad0(iplume)%p(1:ldc%nproma,krow)  = dR_single_spi_trad0(iplume)%p(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnlw(:,1) - aflx_uplw(:,1)) - & 
                (cflx_dnlw(:,1) - cflx_uplw(:,1)) + &
                                                          (dflx_dnlw(:,1) - dflx_uplw(:,1)) - &
                                                          (bflx_dnlw(:,1) - bflx_uplw(:,1))) 
        dR_single_spi_srad0(iplume)%p(1:ldc%nproma,krow)  = dR_single_spi_srad0(iplume)%p(1:ldc%nproma,krow) + &
                                            dtprps*0.5_wp*((aflx_dnsw(:,1) - aflx_upsw(:,1)) - & 
                (cflx_dnsw(:,1) - cflx_upsw(:,1)) + &
                                                          (dflx_dnsw(:,1) - dflx_upsw(:,1)) - &
                                                          (bflx_dnsw(:,1) - bflx_upsw(:,1))) 
        IF (prpraf) THEN
          dR_single_spi_traf0(iplume)%p(1:ldc%nproma,krow)  = dR_single_spi_traf0(iplume)%p(1:ldc%nproma,krow) + &
                                              dtprps*0.5_wp*((aflx_dnlw_clr(:,1) - aflx_uplw_clr(:,1)) - & 
                  (cflx_dnlw_clr(:,1) - cflx_uplw_clr(:,1)) + &
                                                            (dflx_dnlw_clr(:,1) - dflx_uplw_clr(:,1)) - &
                                                            (bflx_dnlw_clr(:,1) - bflx_uplw_clr(:,1))) 
          dR_single_spi_sraf0(iplume)%p(1:ldc%nproma,krow)  = dR_single_spi_sraf0(iplume)%p(1:ldc%nproma,krow) + &
                                              dtprps*0.5_wp*((aflx_dnsw_clr(:,1) - aflx_upsw_clr(:,1)) - & 
                  (cflx_dnsw_clr(:,1) - cflx_upsw_clr(:,1)) + &
                                                            (dflx_dnsw_clr(:,1) - dflx_upsw_clr(:,1)) - &
                                                            (bflx_dnsw_clr(:,1) - bflx_upsw_clr(:,1))) 
        END IF
        IF (prp3doutput) THEN
          IF (prp3drad) THEN
          dR_single_spi_trad(iplume)%p(1:ldc%nproma,1:klev+1,krow)  = dR_single_spi_trad(iplume)%p(1:ldc%nproma,1:klev+1,krow) + &
                                                        dtprps*0.5_wp*((aflx_dnlw(:,:) - aflx_uplw(:,:)) - & 
                                              (cflx_dnlw(:,:) - cflx_uplw(:,:)) + &
                                              (dflx_dnlw(:,:) - dflx_uplw(:,:)) - &
                                              (bflx_dnlw(:,:) - bflx_uplw(:,:))) 
          dR_single_spi_srad(iplume)%p(1:ldc%nproma,1:klev+1,krow)  = dR_single_spi_srad(iplume)%p(1:ldc%nproma,1:klev+1,krow) + &
                                                        dtprps*0.5_wp*((aflx_dnsw(:,:) - aflx_upsw(:,:)) - & 
                                              (cflx_dnsw(:,:) - cflx_upsw(:,:)) + &
                                              (dflx_dnsw(:,:) - dflx_upsw(:,:)) - &
                                              (bflx_dnsw(:,:) - bflx_upsw(:,:))) 
            IF (prpraf) THEN
              dR_single_spi_traf(iplume)%p(1:ldc%nproma,1:klev+1,krow)  = dR_single_spi_traf(iplume)%p(1:ldc%nproma,1:klev+1,krow) + &
                                                              dtprps*0.5_wp*((aflx_dnlw_clr(:,:) - aflx_uplw_clr(:,:)) - & 
                                                    (cflx_dnlw_clr(:,:) - cflx_uplw_clr(:,:)) + &
                                                    (dflx_dnlw_clr(:,:) - dflx_uplw_clr(:,:)) - &
                                                    (bflx_dnlw_clr(:,:) - bflx_uplw_clr(:,:))) 
              dR_single_spi_sraf(iplume)%p(1:ldc%nproma,1:klev+1,krow)  = dR_single_spi_sraf(iplume)%p(1:ldc%nproma,1:klev+1,krow) + &
                                                              dtprps*0.5_wp*((aflx_dnsw_clr(:,:) - aflx_upsw_clr(:,:)) - & 
                                                    (cflx_dnsw_clr(:,:) - cflx_upsw_clr(:,:)) + &
                                                    (dflx_dnsw_clr(:,:) - dflx_upsw_clr(:,:)) - &
                                                    (bflx_dnsw_clr(:,:) - bflx_upsw_clr(:,:)))
            END IF
          END IF
          dQ_single_spi_trad(iplume)%p(1:ldc%nproma,1:klev,krow)  = dQ_single_spi_trad(iplume)%p(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                    ((((aflx_dnlw(:,1:klev) - aflx_uplw(:,1:klev)) - &
                                                      (aflx_dnlw(:,2:klev+1) - aflx_uplw(:,2:klev+1))) - &
                                                      ((cflx_dnlw(:,1:klev) - cflx_uplw(:,1:klev)) - & 
                                                      (cflx_dnlw(:,2:klev+1) - cflx_uplw(:,2:klev+1)))) - &
                                                    (((bflx_dnlw(:,1:klev) - bflx_uplw(:,1:klev)) - & 
                                                      (bflx_dnlw(:,2:klev+1) - bflx_uplw(:,2:klev+1))) - &
                                                      ((dflx_dnlw(:,1:klev) - dflx_uplw(:,1:klev)) - &
                                                      (dflx_dnlw(:,2:klev+1) - dflx_uplw(:,2:klev+1)))))
          dQ_single_spi_srad(iplume)%p(1:ldc%nproma,1:klev,krow)  = dQ_single_spi_srad(iplume)%p(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                    ((((aflx_dnsw(:,1:klev) - aflx_upsw(:,1:klev)) - &
                                                      (aflx_dnsw(:,2:klev+1) - aflx_upsw(:,2:klev+1))) - &
                                                      ((cflx_dnsw(:,1:klev) - cflx_upsw(:,1:klev)) - & 
                                                      (cflx_dnsw(:,2:klev+1) - cflx_upsw(:,2:klev+1)))) - &
                                                    (((bflx_dnsw(:,1:klev) - bflx_upsw(:,1:klev)) - & 
                                                      (bflx_dnsw(:,2:klev+1) - bflx_upsw(:,2:klev+1))) - &
                                                      ((dflx_dnsw(:,1:klev) - dflx_upsw(:,1:klev)) - &
                                                      (dflx_dnsw(:,2:klev+1) - dflx_upsw(:,2:klev+1)))))
          IF (prpraf) THEN
          dQ_single_spi_traf(iplume)%p(1:ldc%nproma,1:klev,krow)  = dQ_single_spi_traf(iplume)%p(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                    ((((aflx_dnlw_clr(:,1:klev) - aflx_uplw_clr(:,1:klev)) - &
                                                      (aflx_dnlw_clr(:,2:klev+1) - aflx_uplw_clr(:,2:klev+1))) - &
                                                      ((cflx_dnlw_clr(:,1:klev) - cflx_uplw_clr(:,1:klev)) - & 
                                                      (cflx_dnlw_clr(:,2:klev+1) - cflx_uplw_clr(:,2:klev+1)))) - &
                                                    (((bflx_dnlw_clr(:,1:klev) - bflx_uplw_clr(:,1:klev)) - & 
                                                      (bflx_dnlw_clr(:,2:klev+1) - bflx_uplw_clr(:,2:klev+1))) - &
                                                      ((dflx_dnlw_clr(:,1:klev) - dflx_uplw_clr(:,1:klev)) - &
                                                      (dflx_dnlw_clr(:,2:klev+1) - dflx_uplw_clr(:,2:klev+1)))))
          dQ_single_spi_sraf(iplume)%p(1:ldc%nproma,1:klev,krow)  = dQ_single_spi_sraf(iplume)%p(1:ldc%nproma,1:klev,krow) + 0.5_wp*dcm(:,1:klev)* &
                                                    ((((aflx_dnsw_clr(:,1:klev) - aflx_upsw_clr(:,1:klev)) - &
                                                      (aflx_dnsw_clr(:,2:klev+1) - aflx_upsw_clr(:,2:klev+1))) - &
                                                      ((cflx_dnsw_clr(:,1:klev) - cflx_upsw_clr(:,1:klev)) - & 
                                                      (cflx_dnsw_clr(:,2:klev+1) - cflx_upsw_clr(:,2:klev+1)))) - &
                                                    (((bflx_dnsw_clr(:,1:klev) - bflx_upsw_clr(:,1:klev)) - & 
                                                      (bflx_dnsw_clr(:,2:klev+1) - bflx_upsw_clr(:,2:klev+1))) - &
                                                      ((dflx_dnsw_clr(:,1:klev) - dflx_upsw_clr(:,1:klev)) - &
                                                      (dflx_dnsw_clr(:,2:klev+1) - dflx_upsw_clr(:,2:klev+1)))))
          END IF
        END IF
      END DO
      END IF
    END IF
  END IF

       CALL message('','Done with PRP')

    END IF
    END IF

  END SUBROUTINE prp

!-----------------------------------------------
!-----------------------------------------------

  SUBROUTINE prp_read_state

    !
    ! This subroutine reads stored data from files. It
    ! is called from stepon, to make sure that all data
    ! is read before PRP calculations are being done.
    !

    USE mo_read_netcdf77,  ONLY: read_var_hs_nf77_3d, read_var_nf77_1d, &
                                 read_diml_nf77, read_var_hs_nf77_2d,   &
                                 read_var_nf77_3d, read_var_hs_nf77_0d
    USE mo_exception,      ONLY: message
    USE mo_submodel,       ONLY: print_value
    USE mo_mpi,            ONLY: p_parallel_io, p_io
    USE mo_decomposition,  ONLY: ldc=>local_decomposition, global_decomposition
    USE mo_transpose,      ONLY: scatter_gp
    USE mo_kind,           ONLY: dp
    USE mo_time_control,   ONLY: get_time_step, current_date, event_state, lstart
    USE mo_time_conversion,ONLY: TC_convert, time_native, print_date, tc_get
    USE mo_time_event,     ONLY: time_event

    IMPLICIT NONE

    integer                        :: ierr, istep, prp_j
    integer                        :: year, month, day, hour, minute, second
    integer                        :: nfiletime, nlat, nlon
    character(len=256)             :: filename
    type(time_native)              :: my_date

    real(dp), POINTER              :: dummy(:,:,:)
    
    real(dp), POINTER              :: zcos_mu0(:,:)
  !  real(dp), POINTER              :: ziaero_real
    real(dp), POINTER              :: zloland_real(:,:)
    real(dp), POINTER              :: zloglac_real(:,:)
   ! real(dp), POINTER              :: zcemiss
    real(dp), POINTER              :: zpgeom1(:,:,:)
    real(dp), POINTER              :: zpp_fl(:,:,:)
    real(dp), POINTER              :: zpp_hl(:,:,:)
    real(dp), POINTER              :: zpp_sfc(:,:)
    real(dp), POINTER              :: zktype_real(:,:)
    
    real(dp), POINTER              :: zq_vap(:,:,:)
    real(dp), POINTER              :: ztk_fl(:,:,:)
    real(dp), POINTER              :: ztk_hl(:,:,:)
    real(dp), POINTER              :: ztk_sfc(:,:)
    real(dp), POINTER              :: zq_liq(:,:,:)
    real(dp), POINTER              :: zq_ice(:,:,:)
    real(dp), POINTER              :: zcdnc(:,:,:)
    real(dp), POINTER              :: zcld_frc(:,:,:)
    real(dp), POINTER              :: zcld_cvr(:,:)    
    real(dp), POINTER              :: zalb_vis_dir(:,:)
    real(dp), POINTER              :: zalb_vis_dif(:,:)
    real(dp), POINTER              :: zalb_nir_dir(:,:)
    real(dp), POINTER              :: zalb_nir_dif(:,:)
    real(dp), POINTER              :: zm_co2(:,:,:)    
    real(dp), POINTER              :: zm_o3(:,:,:)   
    real(dp), POINTER              :: zm_ch4(:,:,:)   
    real(dp), POINTER              :: zm_n2o(:,:,:)   
    real(dp), POINTER              :: zm_cfc1(:,:,:)
    real(dp), POINTER              :: zm_cfc2(:,:,:)    
    real(dp), POINTER              :: zm_o2(:,:,:)     
    real(dp), POINTER              :: ztropo_p(:,:)


    ! Determine if it is time to trigger PRP:
    l_trigprp = event_state(ev_trigprp, current_date) .OR. lstart

    IF (.NOT.prprecord.AND.l_trigprp) THEN

    ! Get and print model time step:
    istep=get_time_step()

    ! Create a local date:
    CALL TC_convert(current_date,my_date)
    CALL print_date(my_date,'Date')
    CALL tc_get(my_date, year, month, day, hour, minute, second)

    ! ------------------

    IF (p_parallel_io) THEN

      ! Read data:
      CALL message('','Reading fields for PRP from netcdf file:')

      ! Determine filename:
      CALL prp_filename(filename,prppath,prp_ref_expname,prp_ref_year,month)
      CALL message('',trim(filename))      

      ! First check for a new month, and if needed read the time from the file:
      IF (oldmonth.ne.month) THEN
        oldmonth =month
        nfiletime=read_diml_nf77(trim(filename),'time')
        nlat     =read_diml_nf77(trim(filename),'lat')
        nlon     =read_diml_nf77(trim(filename),'lon')
        IF (ASSOCIATED(filetime)) DEALLOCATE(filetime,timediff,filedom,filetod)
        ALLOCATE(filetime(nfiletime),timediff(nfiletime))
        ALLOCATE(filedom(nfiletime),filetod(nfiletime))
        ALLOCATE(dummy(nlon,nlat,nfiletime))
        CALL read_var_nf77_3d (trim(filename),'lon','lat','time','tod',dummy,ierr)
        filetod = dummy(1,1,:)
        CALL read_var_nf77_3d (trim(filename),'lon','lat','time','dom',dummy,ierr)
        filedom = dummy(1,1,:)
        DEALLOCATE(dummy)
        filetime = filedom + filetod

      END IF

      ! Determine timestep in file to read: 
      timediff = abs( day + hour/24. + minute/(24.*60.) + second/(24.*60.*60.) &
                     - filetime )
      prp_j    = minloc(timediff,1)

      call print_value('timediff',minval(timediff))
      call print_value('prp_j', prp_j)
      call print_value('filetime(prp_j)', filetime(prp_j))

      ! Allocate temporary variables and read in data:
      ALLOCATE(ztropo_p(ldc%nlon,ldc%nlat))
      ! path to tropopause is different in mpiesm. can not be the work-directory because it is cleared.
      CALL read_var_hs_nf77_2d(trim(prppath)//trim(prp_ref_expname)//'/outdata/echam6/tropopause.nc','lon','lat','time',month,'tropo',ztropo_p,ierr)

    !  ALLOCATE(ziaero_real(0))
    !  CALL read_var_hs_nf77_0d(trim(filename),'time',prp_j,'iaero',ziaero_real,ierr)
      ALLOCATE(zloland_real(ldc%nlon,ldc%nlat))
      CALL read_var_hs_nf77_2d(trim(filename),'lon','lat','time',prp_j,'loland',zloland_real,ierr)
      ALLOCATE(zloglac_real(ldc%nlon,ldc%nlat))
      CALL read_var_hs_nf77_2d(trim(filename),'lon','lat','time',prp_j,'loglac',zloglac_real,ierr)
  !    ALLOCATE(zcemiss(0))
  !    CALL read_var_hs_nf77_0d(trim(filename),'time',prp_j,'cemiss',zcemiss,ierr)
      ALLOCATE(zcos_mu0(ldc%nlon,ldc%nlat))
      CALL read_var_hs_nf77_2d(trim(filename),'lon','lat','time',prp_j,'cos_mu0',zcos_mu0,ierr)
      ALLOCATE(zpgeom1(ldc%nlon,ldc%nlev,ldc%nlat))
      CALL read_var_hs_nf77_3d(trim(filename),'lon','lev','lat','time',prp_j,'pgeom1',zpgeom1,ierr)
      ALLOCATE(zpp_fl(ldc%nlon,ldc%nlev,ldc%nlat))
      CALL read_var_hs_nf77_3d(trim(filename),'lon','lev','lat','time',prp_j,'pp_fl',zpp_fl,ierr)
      ALLOCATE(zpp_hl(ldc%nlon,ldc%nlev+1,ldc%nlat))
      CALL read_var_hs_nf77_3d(trim(filename),'lon','lev_2','lat','time',prp_j,'pp_hl',zpp_hl,ierr)
      ALLOCATE(zpp_sfc(ldc%nlon,ldc%nlat))
      CALL read_var_hs_nf77_2d(trim(filename),'lon','lat','time',prp_j,'pp_sfc',zpp_sfc,ierr)
      ALLOCATE(zktype_real(ldc%nlon,ldc%nlat))
      CALL read_var_hs_nf77_2d(trim(filename),'lon','lat','time',prp_j,'ktype',zktype_real,ierr)

      ALLOCATE(zq_vap(ldc%nlon,ldc%nlev,ldc%nlat))
      CALL read_var_hs_nf77_3d(trim(filename),'lon','lev','lat','time',prp_j,'q_vap',zq_vap,ierr)
      ALLOCATE(ztk_fl(ldc%nlon,ldc%nlev,ldc%nlat))
      CALL read_var_hs_nf77_3d(trim(filename),'lon','lev','lat','time',prp_j,'tk_fl',ztk_fl,ierr)
      ALLOCATE(ztk_hl(ldc%nlon,ldc%nlev+1,ldc%nlat))
      CALL read_var_hs_nf77_3d(trim(filename),'lon','lev_2','lat','time',prp_j,'tk_hl',ztk_hl,ierr)
      ALLOCATE(ztk_sfc(ldc%nlon,ldc%nlat))
      CALL read_var_hs_nf77_2d(trim(filename),'lon','lat','time',prp_j,'tk_sfc',ztk_sfc,ierr)
      ALLOCATE(zq_liq(ldc%nlon,ldc%nlev,ldc%nlat))
      CALL read_var_hs_nf77_3d(trim(filename),'lon','lev','lat','time',prp_j,'q_liq',zq_liq,ierr)
      ALLOCATE(zq_ice(ldc%nlon,ldc%nlev,ldc%nlat))
      CALL read_var_hs_nf77_3d(trim(filename),'lon','lev','lat','time',prp_j,'q_ice',zq_ice,ierr)
      ALLOCATE(zcdnc(ldc%nlon,ldc%nlev,ldc%nlat))
      CALL read_var_hs_nf77_3d(trim(filename),'lon','lev','lat','time',prp_j,'cdnc',zcdnc,ierr)
      ALLOCATE(zcld_frc(ldc%nlon,ldc%nlev,ldc%nlat))
      CALL read_var_hs_nf77_3d(trim(filename),'lon','lev','lat','time',prp_j,'cld_frc',zcld_frc,ierr)
      ALLOCATE(zcld_cvr(ldc%nlon,ldc%nlat))
      CALL read_var_hs_nf77_2d(trim(filename),'lon','lat','time',prp_j,'cld_cvr',zcld_cvr,ierr)
      ALLOCATE(zalb_vis_dir(ldc%nlon,ldc%nlat))
      CALL read_var_hs_nf77_2d(trim(filename),'lon','lat','time',prp_j,'alb_vis_dir',zalb_vis_dir,ierr)
      ALLOCATE(zalb_vis_dif(ldc%nlon,ldc%nlat))
      CALL read_var_hs_nf77_2d(trim(filename),'lon','lat','time',prp_j,'alb_vis_dif',zalb_vis_dif,ierr)
      ALLOCATE(zalb_nir_dir(ldc%nlon,ldc%nlat))
      CALL read_var_hs_nf77_2d(trim(filename),'lon','lat','time',prp_j,'alb_nir_dir',zalb_nir_dir,ierr)
      ALLOCATE(zalb_nir_dif(ldc%nlon,ldc%nlat))
      CALL read_var_hs_nf77_2d(trim(filename),'lon','lat','time',prp_j,'alb_nir_dif',zalb_nir_dif,ierr)
      ALLOCATE(zm_co2(ldc%nlon,ldc%nlev,ldc%nlat))
      CALL read_var_hs_nf77_3d(trim(filename),'lon','lev','lat','time',prp_j,'m_co2',zm_co2,ierr)
      ALLOCATE(zm_o3(ldc%nlon,ldc%nlev,ldc%nlat))
      CALL read_var_hs_nf77_3d(trim(filename),'lon','lev','lat','time',prp_j,'m_o3',zm_o3,ierr)
      ALLOCATE(zm_n2o(ldc%nlon,ldc%nlev,ldc%nlat))
      CALL read_var_hs_nf77_3d(trim(filename),'lon','lev','lat','time',prp_j,'m_n2o',zm_n2o,ierr)
      ALLOCATE(zm_ch4(ldc%nlon,ldc%nlev,ldc%nlat))
      CALL read_var_hs_nf77_3d(trim(filename),'lon','lev','lat','time',prp_j,'m_ch4',zm_ch4,ierr)
      ALLOCATE(zm_cfc1(ldc%nlon,ldc%nlev,ldc%nlat))
      CALL read_var_hs_nf77_3d(trim(filename),'lon','lev','lat','time',prp_j,'m_cfc1',zm_cfc1,ierr)
      ALLOCATE(zm_cfc2(ldc%nlon,ldc%nlev,ldc%nlat))
      CALL read_var_hs_nf77_3d(trim(filename),'lon','lev','lat','time',prp_j,'m_cfc2',zm_cfc2,ierr)
      ALLOCATE(zm_o2(ldc%nlon,ldc%nlev,ldc%nlat))
      CALL read_var_hs_nf77_3d(trim(filename),'lon','lev','lat','time',prp_j,'m_o2',zm_o2,ierr)

    END IF

    CALL scatter_gp(ztropo_p,tropo_p,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(ztropo_p)
    
    ! Distribute water vapor and deallocate:
  !  CALL scatter_gp(ziaero_real,ptr_iaero,global_decomposition)
  !  IF (p_parallel_io) DEALLOCATE(ziaero_real)
    CALL scatter_gp(zloland_real,ptr_loland,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zloland_real)
    CALL scatter_gp(zloglac_real,ptr_loglac,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zloglac_real)
  !  CALL scatter_gp(zcemiss,ptr_cemiss,global_decomposition)
  !  IF (p_parallel_io) DEALLOCATE(zcemiss)
    CALL scatter_gp(zcos_mu0,ptr_cos_mu0,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zcos_mu0)
    CALL scatter_gp(zpgeom1,ptr_pgeom1,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zpgeom1)
    CALL scatter_gp(zpp_fl,ptr_pp_fl,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zpp_fl)
    CALL scatter_gp(zpp_hl,ptr_pp_hl,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zpp_hl)
    CALL scatter_gp(zpp_sfc,ptr_pp_sfc,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zpp_sfc)
    CALL scatter_gp(zktype_real,ptr_ktype,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zktype_real)

    CALL scatter_gp(zq_vap,ptr_q_vap,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zq_vap)
    CALL scatter_gp(ztk_fl,ptr_tk_fl,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(ztk_fl)
    CALL scatter_gp(ztk_hl,ptr_tk_hl,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(ztk_hl)
    CALL scatter_gp(ztk_sfc,ptr_tk_sfc,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(ztk_sfc)
    CALL scatter_gp(zq_liq,ptr_q_liq,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zq_liq)
    CALL scatter_gp(zq_ice,ptr_q_ice,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zq_ice)
    CALL scatter_gp(zcdnc,ptr_cdnc,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zcdnc)
    CALL scatter_gp(zcld_frc,ptr_cld_frc,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zcld_frc)
    CALL scatter_gp(zcld_cvr,ptr_cld_cvr,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zcld_cvr)
    CALL scatter_gp(zalb_vis_dir,ptr_alb_vis_dir,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zalb_vis_dir)
    CALL scatter_gp(zalb_vis_dif,ptr_alb_vis_dif,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zalb_vis_dif)
    CALL scatter_gp(zalb_nir_dir,ptr_alb_nir_dir,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zalb_nir_dir)
    CALL scatter_gp(zalb_nir_dif,ptr_alb_nir_dif,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zalb_nir_dif)
    CALL scatter_gp(zm_co2,ptr_m_co2,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zm_co2)

    CALL scatter_gp(zm_o3,ptr_m_o3,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zm_o3)
    CALL scatter_gp(zm_ch4,ptr_m_ch4,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zm_ch4)
    CALL scatter_gp(zm_n2o,ptr_m_n2o,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zm_n2o)
    CALL scatter_gp(zm_cfc1,ptr_m_cfc1,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zm_cfc1)
    CALL scatter_gp(zm_cfc2,ptr_m_cfc2,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zm_cfc2)
    CALL scatter_gp(zm_o2,ptr_m_o2,global_decomposition)
    IF (p_parallel_io) DEALLOCATE(zm_o2)

    CALL message('','Done reading PRP data')

    END IF

  END SUBROUTINE prp_read_state
!-----------------------------------------------
!-----------------------------------------------

  SUBROUTINE prp_filename(filename,xpath,xexp,year,month)

    IMPLICIT NONE

    character(*)             :: filename
    integer                  :: year,month
    character(len=4)         :: yearstring
    character(len=2)         :: monthstring
    character(len=256)       :: xpath
    character(len=64)        :: xexp

    ! Create substrings for year and month to insert in filename:

    write(yearstring,'(i4)') year 
    IF (year.lt.1000) THEN
      yearstring    = adjustr(yearstring)
      yearstring(1:1) = '0'
    END IF
    write(monthstring,'(i2)') month
    IF (month.lt.10) THEN
      monthstring    = adjustr(monthstring)
      monthstring(1:1) = '0'
    END IF

    ! Concatenate the filename and path:
    
    ! file path when used in echam standalone                                                                                                                                                                                                                       
    filename=trim(xpath)//trim(xexp)//'/'//trim(xexp)//'_'//trim(yearstring) &
            //trim(monthstring)//'.01_prp.nc'

    ! file path when used in mpiesm                                                                                                                   
    IF (prp_mpiesm) THEN
        filename=trim(xpath)//trim(xexp)//'/work/run_'//trim(yearstring)// &
               '0101-'//trim(yearstring)//'1231/'//trim(xexp)//'_'//trim(yearstring)// &
       trim(monthstring)//'.01_prp.nc'
    END IF
    
  END SUBROUTINE prp_filename

!-----------------------------------------------
!-----------------------------------------------

END MODULE mo_prp
