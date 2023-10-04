!>
!!
!! @brief Module MO_SIMPLE_PLUMES: provides anthropogenic aerorol optial properties as a function of lat, lon
!!   height, time, and wavelength
!!
!! @remarks
!!
!! @author Bjorn Stevens & Karsten Peters MPI-M, Hamburg (v1-beta release 2015-12-05)
!!
!!         Stephanie Fiedler, MPI-M, Hamburg (bug fixes for annual cycle, Twomey effect, 2016-09-22)
!!         Stephanie Fiedler, MPI-M, Hamburg (revised vertical distribution to meter above sea level, 2016-09-29)
!!         Stephanie Fiedler, MPI-M, Hamburg (corrected artifical gradients, 2016-10-11)
!! 
!! $ID: n/a$
!!
!! @par Origin
!!   Based on code originally developed at the MPI by Karsten Peters, Bjorn Stevens, Stephanie Fiedler
!!   and Stefan Kinne with input from Thorsten Mauritsen and Robert Pincus
!!
!! @par Copyright
!! 
!
MODULE MO_SIMPLE_PLUMES

  USE mo_kind,             ONLY: dp
  USE mo_mpi,              ONLY: p_parallel_io, p_io, p_bcast, p_parallel
  USE mo_read_netcdf77,    ONLY: read_diml_nf77, read_var_nf77_1d, &
                                 read_var_nf77_2d, read_var_nf77_3d
  USE mo_exception,        ONLY: finish, message, message_text
  USE mo_time_conversion,  ONLY: &
         day_in_year,            &     !< get the present day in year
         time_days,              &     !< type for time variables
         time_native,            &     !< type for time variables
         TC_get, TC_set, TC_convert    !< get time components, set and convert dates
  USE mo_time_base,        ONLY: &     
         Get_JulianMonLen              !< get Julian month length
  USE mo_vphysc,           ONLY: &
         vphysc                        !< derived type containing geopotential at interaces
  USE mo_memory_g3b,       ONLY: & 
         geosp                         !< surface geopotential
  USE mo_geoloc,           ONLY: &
         philat_2d              ,&     !< 2d array of latitudes
         philon_2d                     !< 2d array of longitudes
  USE mo_physical_constants,ONLY: &
         grav                          !< gravity
  USE mo_srtm_setup,       ONLY: &
         sw_wv1 => wavenum1     ,&     !< smallest wave number in each of the sw bands
         sw_wv2 => wavenum2            !< largest wave number in each of the sw bands
  USE mo_rrtm_params,      ONLY: &
         jpb1                   ,&     !< index for lower sw band
         jpb2                          ! index for upper sw band
  USE mo_math_constants,   ONLY: &
         pi                     ,&     ! pi
         deg2rad                       ! pi/180
  USE mo_namelist,         ONLY: open_nml, position_nml, POSITIONED, &
                                 READ_ERROR, LENGTH_ERROR
  USE mo_submodel,       ONLY: print_value
  !USE mo_aod_sp_pointer, ONLY: ptr_aod_sp, ptr_aod_bg
  USE mo_prp,            ONLY: ptr_aod_sp, ptr_aod_bg

  IMPLICIT NONE

  INTEGER, SAVE :: fix_yr = -9999   !< Fix plume parameters to a given year, unless fix_yr==-9999
  REAL(dp), SAVE:: bg_multiplier = 1.0_dp ! Multiplier of background aerosol for calculating Twomey effect
  REAL(dp), SAVE:: spd_multiplier = 1.0_dp ! Multiplier of simple plumes direct aerosol effect

  PUBLIC            :: add_aop_sp

  CHARACTER(LEN=256), SAVE  :: cfname='MAC-SP.nc'

  INTEGER, PARAMETER ::                        &
       nplumes   = 9                          ,& !< Number of plumes
       nfeatures = 2                          ,& !< Number of features per plume
       ntimes    = 52                         ,& !< Number of times resolved per year (52 => weekly resolution)
       nyears    = 251                           !< Number of years of available forcing

  LOGICAL, SAVE ::                             &
       sp_initialized = .FALSE.                  !< parameter determining whether input needs to be read

  REAL(dp) ::                                  &
       plume_lat   (nplumes)                  ,& !< latitude where plume maximizes
       plume_lon   (nplumes)                  ,& !< longitude where plume maximizes
       beta_a      (nplumes)                  ,& !< parameter a for beta function vertical profile
       beta_b      (nplumes)                  ,& !< parameter b for beta function vertical profile
       aod_spmx    (nplumes)                  ,& !< aod at 550 for simple plume (maximum)
       aod_fmbg    (nplumes)                  ,& !< aod at 550 for fine mode natural background (for twomey effect)
       asy550      (nplumes)                  ,& !< asymmetry parameter for plume at 550nm
       ssa550      (nplumes)                  ,& !< single scattering albedo for plume at 550nm
       angstrom    (nplumes)                  ,& !< angstrom parameter for plume 
       sig_lon_E   (nfeatures,nplumes)        ,& !< Eastward extent of plume feature
       sig_lon_W   (nfeatures,nplumes)        ,& !< Westward extent of plume feature
       sig_lat_E   (nfeatures,nplumes)        ,& !< Southward extent of plume feature
       sig_lat_W   (nfeatures,nplumes)        ,& !< Northward extent of plume feature
       theta       (nfeatures,nplumes)        ,& !< Rotation angle of feature
       ftr_weight  (nfeatures,nplumes)        ,& !< Feature weights = (nfeatures + 1) to account for BB background
       time_weight (nfeatures,nplumes)        ,& !< Time-weights = (nfeatures +1) to account for BB background
       time_weight_bg (nfeatures,nplumes)     ,& !< as time_wight but for natural background in Twomey effect
       year_weight (nyears,nplumes)           ,& !< Yearly weight for plume
       ann_cycle   (nfeatures,ntimes,nplumes)    !< annual cycle for feature

  PUBLIC sp_aop_profile

CONTAINS
  !
  ! ------------------------------------------------------------------------------------------------------------------------
  ! SP_SETUP:  This subroutine should be called at initialization to read the netcdf data that describes the simple plume
  ! climatology.  The information needs to be either read by each processor or distributed to processors.
  !
  SUBROUTINE sp_setup
    !
    ! ---------- 
    !
    LOGICAL :: lex
    INTEGER :: xdmy,ierr, inml, iunit

    INCLUDE 'simpleplumesctl.inc'
    
    !
    ! ---------- 
    !    
    IF (p_parallel_io) THEN

      inml = open_nml ('namelist.echam')
      iunit = position_nml ('SIMPLEPLUMESCTL', inml, status=ierr)
      SELECT CASE (ierr)
      CASE (POSITIONED)
        READ(iunit, simpleplumesctl)
      CASE (LENGTH_ERROR)
          CALL finish ('mo_simple-plumes','internal error: namelist name is too long')
      CASE (READ_ERROR)
          CALL finish ('mo_simple-plumes','read error while seeking namelist simpleplumesctl')
      END SELECT
 
      IF (fix_yr == -9999) THEN
        WRITE (message_text,*) 'using transient anthropogenic aerosols'
        CALL message('mo_simple-plumes', message_text) 
      ELSE IF ((1850 <= fix_yr) .AND. (fix_yr <= 2100)) THEN
        WRITE (message_text,*) 'using constant ', fix_yr, 'anthropogenic aerosols'
        CALL message('mo_simple-plumes', message_text) 
      ELSE 
        CALL finish ('mo_simple-plumes','selected year out of 1850-2100 range')
      END IF 

      INQUIRE (file=TRIM(cfname), exist=lex)
      IF (.NOT. lex) THEN
        CALL finish('sp_setup, mo_simple-plumes', &
                    'file '//TRIM(cfname)//' does not exist')
      END IF
      xdmy=read_diml_nf77(TRIM(cfname),'plume_number')
      IF (xdmy /= nplumes) THEN
        CALL finish('sp_setup, mo_simple-plumes', &
                    'NetCDF improperly dimensioned -- plume_number' &
                    //TRIM(cfname))
      END IF
      xdmy=read_diml_nf77(TRIM(cfname),'plume_feature')
      IF (xdmy /= nfeatures) THEN
        CALL finish('sp_setup, mo_simple-plumes', &
                    'NetCDF improperly dimensioned -- plume_feature' &
                    //TRIM(cfname))
      END IF
      xdmy=read_diml_nf77(TRIM(cfname),'year_fr')
      IF (xdmy /= ntimes) THEN
        CALL finish('sp_setup, mo_simple-plumes', &
                    'NetCDF improperly dimensioned -- year_fr' &
                    //TRIM(cfname))
      END IF
      xdmy=read_diml_nf77(TRIM(cfname),'years')
      IF (xdmy /= nyears) THEN
        CALL finish('sp_setup, mo_simple-plumes', &
                    'NetCDF improperly dimensioned -- years' &
                    //TRIM(cfname))
      END IF
      !
      ! read variables that define the simple plume climatology
      !
      CALL read_var_nf77_1d(TRIM(cfname),'plume_number','plume_lat',plume_lat,ierr)
      CALL read_var_nf77_1d(TRIM(cfname),'plume_number','plume_lon',plume_lon,ierr)
      CALL read_var_nf77_1d(TRIM(cfname),'plume_number','beta_a',beta_a,ierr)
      CALL read_var_nf77_1d(TRIM(cfname),'plume_number','beta_b',beta_b,ierr)
      CALL read_var_nf77_1d(TRIM(cfname),'plume_number','aod_spmx',aod_spmx,ierr)
      CALL read_var_nf77_1d(TRIM(cfname),'plume_number','aod_fmbg',aod_fmbg,ierr)
      CALL read_var_nf77_1d(TRIM(cfname),'plume_number','ssa550',ssa550,ierr)
      CALL read_var_nf77_1d(TRIM(cfname),'plume_number','asy550',asy550,ierr)
      CALL read_var_nf77_1d(TRIM(cfname),'plume_number','angstrom',angstrom,ierr)
      CALL read_var_nf77_2d(TRIM(cfname),'plume_feature','plume_number', &
           'sig_lat_W',sig_lat_W,ierr)
      CALL read_var_nf77_2d(TRIM(cfname),'plume_feature','plume_number', &
           'sig_lat_E',sig_lat_E,ierr)
      CALL read_var_nf77_2d(TRIM(cfname),'plume_feature','plume_number', &
           'sig_lon_W',sig_lon_W,ierr)
      CALL read_var_nf77_2d(TRIM(cfname),'plume_feature','plume_number', &
           'sig_lon_E',sig_lon_E,ierr)
      CALL read_var_nf77_2d(TRIM(cfname),'plume_feature','plume_number', &
           'theta',theta,ierr)
      CALL read_var_nf77_2d(TRIM(cfname),'plume_feature','plume_number', &
           'ftr_weight',ftr_weight,ierr)
      CALL read_var_nf77_2d(TRIM(cfname),'years','plume_number', &
           'year_weight',year_weight,ierr)
      CALL read_var_nf77_3d(TRIM(cfname),'plume_feature','year_fr','plume_number', &
           'ann_cycle',ann_cycle,ierr)
    END IF

    IF (p_parallel) THEN
      CALL p_bcast(fix_yr,p_io)
      CALL p_bcast(bg_multiplier,p_io)
      CALL p_bcast(spd_multiplier,p_io)
      CALL p_bcast(plume_lat,p_io)
      CALL p_bcast(plume_lon,p_io)
      CALL p_bcast(beta_a,p_io)
      CALL p_bcast(beta_b,p_io)
      CALL p_bcast(aod_spmx,p_io)
      CALL p_bcast(aod_fmbg,p_io)
      CALL p_bcast(ssa550,p_io)
      CALL p_bcast(asy550,p_io)
      CALL p_bcast(angstrom,p_io)
      CALL p_bcast(sig_lat_W,p_io)
      CALL p_bcast(sig_lat_E,p_io)
      CALL p_bcast(sig_lon_W,p_io)
      CALL p_bcast(sig_lon_E,p_io)
      CALL p_bcast(theta,p_io)
      CALL p_bcast(ftr_weight,p_io)
      CALL p_bcast(year_weight,p_io)
      CALL p_bcast(ann_cycle,p_io)
    END IF

    sp_initialized = .TRUE.

    RETURN
  END SUBROUTINE sp_setup

  !
  ! ------------------------------------------------------------------------------------------------------------------------
  ! SET_TIME_WEIGHT:  The simple plume model assumes that meteorology constrains plume shape and that only source strength
  ! influences the amplitude of a plume associated with a given source region.   This routine retrieves the temporal weights
  ! for the plumes.  Each plume feature has its own temporal weights which varies yearly.  The annual cycle is indexed by
  ! week in the year and superimposed on the yearly mean value of the weight. 
  !
  SUBROUTINE set_time_weight(year_fr)
    !
    ! ---------- 
    !
    REAL(dp), INTENT(IN) ::  &
         year_fr           !< Fractional Year (1850.0 - 2100.99)

    INTEGER          ::  &
         iyear          ,& !< Integer year values between 1 and 156 (1850-2100) 
         iweek          ,& !< Integer index (between 1 and ntimes); for ntimes=52 this corresponds to weeks (roughly)
         iplume            ! plume number
    !
    ! ---------- 
    !
    iyear = MERGE(FLOOR(year_fr), fix_yr, fix_yr == -9999) - 1849
    iweek = FLOOR((year_fr - FLOOR(year_fr)) * ntimes) + 1

    IF ((iweek > ntimes) .OR. (iweek < 1) .OR. (iyear > nyears) .OR. (iyear < 1)) STOP 'Time out of bounds in set_time_weight'
    DO iplume=1,nplumes
      time_weight(1,iplume) = year_weight(iyear,iplume) * ann_cycle(1,iweek,iplume)
      time_weight(2,iplume) = year_weight(iyear,iplume) * ann_cycle(2,iweek,iplume)
      time_weight_bg(1,iplume) = ann_cycle(1,iweek,iplume)
      time_weight_bg(2,iplume) = ann_cycle(2,iweek,iplume)  
    END DO
    
    RETURN
  END SUBROUTINE set_time_weight
  !
  ! ------------------------------------------------------------------------------------------------------------------------
  ! SP_AOP_PROFILE:  This subroutine calculates the simple plume aerosol and cloud active optical properites based on the
  ! the simple plume fit to the MPI Aerosol Climatology (Version 2).  It sums over nplumes to provide a profile of aerosol
  ! optical properties on a host models vertical grid. 
  !
  SUBROUTINE sp_aop_profile                                                                           ( &
       nlevels        ,ncol           ,lambda         ,oro            ,lon            ,lat            , &
       year_fr        ,z              ,dz             ,dNovrN         ,aod_prof       ,ssa_prof       , &
       asy_prof       ,plume_number   ,krow)
    !
    ! ---------- 
    !
    USE mo_decomposition,  ONLY: ldc => local_decomposition
    !
    ! ---------- 
    !
    INTEGER, INTENT(IN)        :: &
         nlevels,                 & !< number of levels
         krow,                    & !< index for current block
         ncol,                    & !< number of columns
         plume_number

    REAL(dp), INTENT(IN)       :: &
         lambda,                  & !< wavelength
         year_fr,                 & !< Fractional Year (1903.0 is the 0Z on the first of January 1903, Gregorian)
         oro(ncol),               & !< orographic height (m)
         lon(ncol),               & !< longitude 
         lat(ncol),               & !< latitude
         z (ncol,nlevels),        & !< height above sea-level (m)
         dz(ncol,nlevels)           !< level thickness (difference between half levels)

    REAL(dp), INTENT(OUT)      :: &
         dNovrN(ncol)           , & !< anthropogenic increment to cloud drop number concentration 
         aod_prof(ncol,nlevels) , & !< profile of aerosol optical depth
         ssa_prof(ncol,nlevels) , & !< profile of single scattering albedo
         asy_prof(ncol,nlevels)     !< profile of asymmetry parameter

    INTEGER                    :: iplume, icol, k, inplume, jnplume

    REAL(dp)                   ::  &
         eta(ncol,nlevels),        & !< normalized height (by 15 km)
         z_beta(ncol,nlevels),     & !< profile for scaling column optical depth
         prof(ncol,nlevels),       & !< scaled profile (by beta function)
         beta_sum(ncol),           & !< vertical sum of beta function
         ssa(ncol),                & !< aerosol optical depth 
         asy(ncol),                & !< aerosol optical depth 
         cw_an(ncol),              & !< column weight for simple plume (anthropogenic) aod at 550 nm
         cw_bg(ncol),              & !< column weight for fine-mode industrial background aod at 550 nm
         caod_sp(ncol),            & !< column simple plume (anthropogenic) aod at 550 nm
         caod_bg(ncol),            & !< column fine-mode natural background aod at 550 nm
         a_plume1,                 & !< gaussian longitude factor for feature 1
         a_plume2,                 & !< gaussian longitude factor for feature 2
         b_plume1,                 & !< gaussian latitude factor for feature 1
         b_plume2,                 & !< gaussian latitude factor for feature 2
         delta_lat,                & !< latitude offset
         delta_lon,                & !< longitude offset
         delta_lon_t,              & !< threshold for maximum longitudinal plume extent used in transition from 360 to 0 degrees
         lon1,                     & !< rotated longitude for feature 1
         lat1,                     & !< rotated latitude for feature 2
         lon2,                     & !< rotated longitude for feature 1
         lat2,                     & !< rotated latitude for feature 2
         f1,                       & !< contribution from feature 1
         f2,                       & !< contribution from feature 2
         f3,                       & !< contribution from feature 1 in natural background of Twomey effect
         f4,                       & !< contribution from feature 2 in natural background of Twomey effect
         aod_550,                  & !< aerosol optical depth at 550nm
         aod_lmd,                  & !< aerosol optical depth at input wavelength
         lfactor                     !< factor to compute wavelength dependence of optical properties
    !
    ! ---------- 
    !
    ! initialize input data (by calling setup at first instance) 
    !
    IF (.NOT.sp_initialized) CALL sp_setup
    !
    ! get time weights
    !
    CALL set_time_weight(year_fr)
    !
    ! initialize variables, including output
    !
    DO k=1,nlevels
      DO icol=1,ncol
        aod_prof(icol,k) = 0.0
        ssa_prof(icol,k) = 0.0
        asy_prof(icol,k) = 0.0
        z_beta(icol,k)   = MERGE(1.0, 0.0, z(icol,k) >= oro(icol))
        eta(icol,k)      = MAX(0.0,MIN(1.0,z(icol,k)/15000.))
      END DO
    END DO
    DO icol=1,ncol
      dNovrN(icol)   = 1.0
      caod_sp(icol)  = 0.00
      caod_bg(icol)  = 0.02*bg_multiplier
    END DO
    !
    ! sum contribution from plumes to construct composite profiles of aerosol otpical properties
    !
    ! conditions to calculate single plume only
    SELECT CASE (plume_number)
    CASE (:9)                     ! total effect single plume
      inplume=plume_number
      jnplume=inplume
    CASE (21:29)                  ! direct effect single plume
      inplume=(plume_number-20)
      jnplume=inplume
    CASE (41:49)                  ! indirect effect single plume
      inplume=(plume_number-40)
      jnplume=inplume
    CASE DEFAULT                  ! loop over all plumes
      inplume=1  
      jnplume=nplumes
    END SELECT
    DO iplume=inplume,jnplume
      IF ((iplume+10)==plume_number) cycle   ! skip single plume for backward prp in total aerosol effect
        !
        ! calculate vertical distribution function from parameters of beta distribution
        !
        DO icol=1,ncol
          beta_sum(icol) = 0.
        END DO
        DO k=1,nlevels
          DO icol=1,ncol
            prof(icol,k)   = (eta(icol,k)**(beta_a(iplume)-1.) * (1.-eta(icol,k))**(beta_b(iplume)-1.))*dz(icol,k)
            beta_sum(icol) = beta_sum(icol) + prof(icol,k)
          END DO
        END DO
        DO k=1,nlevels
          DO icol=1,ncol
            prof(icol,k)   = prof(icol,k) / beta_sum(icol) * z_beta(icol,k)
          END DO
        END DO
        !
        ! calculate plume weights
        !
        DO icol=1,ncol
          !
          ! get plume-center relative spatial parameters for specifying amplitude of plume at given lat and lon
          !
          delta_lat = lat(icol) - plume_lat(iplume)
          delta_lon = lon(icol) - plume_lon(iplume)
          delta_lon_t = MERGE (260._dp, 180._dp, iplume == 1)
          delta_lon = MERGE ( delta_lon-SIGN(360._dp,delta_lon) , delta_lon , ABS(delta_lon) > delta_lon_t)

          a_plume1  = 0.5_dp / (MERGE(sig_lon_E(1,iplume), sig_lon_W(1,iplume), delta_lon > 0)**2)
          b_plume1  = 0.5_dp / (MERGE(sig_lat_E(1,iplume), sig_lat_W(1,iplume), delta_lon > 0)**2)
          a_plume2  = 0.5_dp / (MERGE(sig_lon_E(2,iplume), sig_lon_W(2,iplume), delta_lon > 0)**2)
          b_plume2  = 0.5_dp / (MERGE(sig_lat_E(2,iplume), sig_lat_W(2,iplume), delta_lon > 0)**2)
          !
          ! adjust for a plume specific rotation which helps match plume state to climatology.
          !
          lon1 =   COS(theta(1,iplume))*(delta_lon) + SIN(theta(1,iplume))*(delta_lat)
          lat1 = - SIN(theta(1,iplume))*(delta_lon) + COS(theta(1,iplume))*(delta_lat)
          lon2 =   COS(theta(2,iplume))*(delta_lon) + SIN(theta(2,iplume))*(delta_lat)
          lat2 = - SIN(theta(2,iplume))*(delta_lon) + COS(theta(2,iplume))*(delta_lat)
          !
          ! calculate contribution to plume from its different features, to get a column weight for the anthropogenic
          ! (cw_an) and the fine-mode background aerosol (cw_bg)
          !
          f1 = time_weight(1,iplume) * ftr_weight(1,iplume) * EXP(-1._dp* (a_plume1 * ((lon1)**2) + (b_plume1 * ((lat1)**2)))) 
          f2 = time_weight(2,iplume) * ftr_weight(2,iplume) * EXP(-1._dp* (a_plume2 * ((lon2)**2) + (b_plume2 * ((lat2)**2)))) 
          f3 = time_weight_bg(1,iplume) * ftr_weight(1,iplume) * EXP(-1.* (a_plume1 * ((lon1)**2) + (b_plume1 * ((lat1)**2)))) 
          f4 = time_weight_bg(2,iplume) * ftr_weight(2,iplume) * EXP(-1.* (a_plume2 * ((lon2)**2) + (b_plume2 * ((lat2)**2))))

          cw_an(icol) = f1 * aod_spmx(iplume) + f2 * aod_spmx(iplume)  
          cw_bg(icol) = f3 * aod_fmbg(iplume) + f4 * aod_fmbg(iplume) 
          !
          ! calculate wavelength-dependent scattering properties
          !
          lfactor   = MIN(1.0_dp,700.0_dp/lambda)
          ssa(icol) = (ssa550(iplume) * lfactor**4) / ((ssa550(iplume) * lfactor**4) + ((1-ssa550(iplume)) * lfactor))
          asy(icol) =  asy550(iplume) * SQRT(lfactor)
        END DO
        !
        ! distribute plume optical properties across its vertical profile weighting by optical depth and scaling for
        ! wavelength using the anstrom parameter. 
        !      
        lfactor = EXP(-angstrom(iplume) * LOG(lambda/550.0_dp))
        DO k=1,nlevels
          DO icol = 1,ncol
            SELECT CASE (plume_number)
            CASE (:19)  ! for 1 to 9: total effect individual plumes, for 10: total effect all plumes, for 11 to 19: total effect all plumes but individual plume
                aod_550          = prof(icol,k)     * cw_an(icol)
                aod_lmd          = aod_550          * lfactor
                caod_sp(icol)    = caod_sp(icol)    + aod_550
                caod_bg(icol)    = caod_bg(icol)    + prof(icol,k) * cw_bg(icol)
                asy_prof(icol,k) = asy_prof(icol,k) + aod_lmd * ssa(icol) * asy(icol)
                ssa_prof(icol,k) = ssa_prof(icol,k) + aod_lmd * ssa(icol)
                aod_prof(icol,k) = aod_prof(icol,k) + aod_lmd
                
            CASE (21:29)    ! direct effect single plume
                aod_550          = prof(icol,k)     * cw_an(icol)
                aod_lmd          = aod_550          * lfactor
                asy_prof(icol,k) = asy_prof(icol,k) + aod_lmd * ssa(icol) * asy(icol)
                ssa_prof(icol,k) = ssa_prof(icol,k) + aod_lmd * ssa(icol)
                aod_prof(icol,k) = aod_prof(icol,k) + aod_lmd

            CASE (31:39)    ! total effect all plumes but indirect effect of single plume
                aod_550          = prof(icol,k)     * cw_an(icol)
                aod_lmd          = aod_550          * lfactor
                asy_prof(icol,k) = asy_prof(icol,k) + aod_lmd * ssa(icol) * asy(icol)
                ssa_prof(icol,k) = ssa_prof(icol,k) + aod_lmd * ssa(icol)
                aod_prof(icol,k) = aod_prof(icol,k) + aod_lmd

                IF ((iplume+30)/=plume_number) THEN ! if not single plume, add indirect effect of other plumes                 
                  caod_sp(icol)    = caod_sp(icol)    + aod_550
                  caod_bg(icol)    = caod_bg(icol)    + prof(icol,k) * cw_bg(icol)  
                ENDIF
            
            CASE (41:49)     ! indirect effect single plume
                aod_550          = prof(icol,k)     * cw_an(icol)
                caod_sp(icol)    = caod_sp(icol)    + aod_550
                caod_bg(icol)    = caod_bg(icol)    + prof(icol,k) * cw_bg(icol)

            CASE (51:59)     ! total effect all plumes but direct effect of single plume
                aod_550          = prof(icol,k)     * cw_an(icol)
                caod_sp(icol)    = caod_sp(icol)    + aod_550
                caod_bg(icol)    = caod_bg(icol)    + prof(icol,k) * cw_bg(icol)

                IF ((iplume+50)/=plume_number) THEN ! if not single plume, add direct effect of other plumes
                  aod_lmd          = aod_550          * lfactor
                  asy_prof(icol,k) = asy_prof(icol,k) + aod_lmd * ssa(icol) * asy(icol)
                  ssa_prof(icol,k) = ssa_prof(icol,k) + aod_lmd * ssa(icol)
                  aod_prof(icol,k) = aod_prof(icol,k) + aod_lmd
                ENDIF
                  
            CASE DEFAULT
            END SELECT
          END DO
        END DO
    END DO
    !
    ! complete optical depth weighting
    !
    DO k=1,nlevels
      DO icol = 1,ncol
        asy_prof(icol,k) = MERGE(asy_prof(icol,k)/ssa_prof(icol,k), 0.0_dp, ssa_prof(icol,k) > TINY(1._dp))
        ssa_prof(icol,k) = MERGE(ssa_prof(icol,k)/aod_prof(icol,k), 1.0_dp, aod_prof(icol,k) > TINY(1._dp))
      END DO
    END DO
    !
    ! calcuate effective radius normalization (divisor) factor
    !
    DO icol=1,ncol
      dNovrN(icol) = LOG((1000.0_dp * (caod_sp(icol) + caod_bg(icol))) + 1.0_dp)/LOG((1000.0_dp * caod_bg(icol)) + 1.0_dp)
    END DO

    IF (plume_number<11) THEN
      ptr_aod_sp(plume_number)%p(1:ldc%nproma,krow) = caod_sp(:)
      ptr_aod_bg(plume_number)%p(1:ldc%nproma,krow) = caod_bg(:)
    END IF

    RETURN
  END SUBROUTINE sp_aop_profile

  SUBROUTINE sp_reference_time(prev_rad_date,            rad_date,            &
                               prev_radiation_date_1850, radiation_date_1850)

    TYPE (time_days), INTENT(in) :: rad_date, prev_rad_date
    TYPE (time_days), INTENT(out):: radiation_date_1850,prev_radiation_date_1850
    TYPE (time_native) :: xdate_1850
    INTEGER            :: xdate_yy,xdate_mm,xdate_dd,xdate_hr,xdate_mi,xdate_se
         call TC_convert(rad_date,xdate_1850)
         call TC_get(xdate_1850,xdate_yy,xdate_mm,xdate_dd,xdate_hr,xdate_mi,xdate_se)
         IF (xdate_mm == 2) xdate_dd = min(28,xdate_dd)
         call TC_set(1850,xdate_mm,xdate_dd,xdate_hr,xdate_mi,xdate_se,xdate_1850)
         call TC_convert(xdate_1850,radiation_date_1850)

         call TC_convert(prev_rad_date,xdate_1850)
         call TC_get(xdate_1850,xdate_yy,xdate_mm,xdate_dd,xdate_hr,xdate_mi,xdate_se)
         IF (xdate_mm == 2) xdate_dd = min(28,xdate_dd)
         call TC_set(1850,xdate_mm,xdate_dd,xdate_hr,xdate_mi,xdate_se,xdate_1850)
         call TC_convert(xdate_1850,prev_radiation_date_1850)
    
  END SUBROUTINE sp_reference_time
  ! ------------------------------------------------------------------------------------------------------------------------
  ! ADD_AOP_SP:  This subroutine provides the interface to simple plume (sp) fit to the MPI Aerosol Climatology (Version 2).
  ! It does so by collecting or deriving spatio-temporal information and calling the simple plume aerosol subroutine and
  ! incrementing the background aerosol properties (and effective radisu) with the anthropogenic plumes.
  !
  SUBROUTINE add_aop_sp                                                           ( &
       kproma         ,kbdim          ,klev           ,krow           ,nb_sw      , &
       aod_sw_vr      ,ssa_sw_vr      ,asy_sw_vr      ,x_cdnc         ,plume_number )
    !
    ! --- 0.0 Use-associated variables
    !
    USE mo_time_control,       ONLY: &
         next_date,                  & !< date at the next timestep (derived type variable)
         get_date_components           !< function returning date components (years, months, etc) optional arguments
    !
    ! --- 0.1 Variables passed through argument list
    INTEGER, INTENT(IN) ::            &
         kproma                      ,& !< number of elements in current block
         kbdim                       ,& !< block dimension (greater than or equal to kproma)
         klev                        ,& !< number of full levels
         krow                        ,& !< index for current block
         nb_sw                       ,& !< number of bands in short wave
         plume_number                   !< plume number to use: for all plumes -> 10

    REAL(dp), INTENT (INOUT) ::       &
         aod_sw_vr(kbdim,klev,nb_sw) ,& !< Aerosol shortwave optical depth
         ssa_sw_vr(kbdim,klev,nb_sw) ,& !< Aerosol single scattering albedo
         asy_sw_vr(kbdim,klev,nb_sw) ,& !< Aerosol asymmetry parameter
         x_cdnc(kbdim)                  !< Scale factor for Cloud Droplet Number Concentration
    !
    ! --- 0.2 Dummy variables
    !
    INTEGER ::                        &
         idy                         ,& !< Day of present year
         iyr                         ,& !< year
         jk                          ,& !< index for looping over vertical dimension
         jki                         ,& !< index for looping over vertical dimension for reversing
         jl                          ,& !< index for looping over block
         jwl                         ,& !< index for looping over wavelengths
         j_sw                           !< index for looping over wavelengths
    
    REAL(dp) ::                       &
         days_in_year                ,& !< Total number of days in present year
         year_fr                     ,& !< time in year fraction (1989.0 is 0Z on Jan 1 1989)
         lambda                      ,& !< wavelength at central band wavenumber [nm]
         z_sfc(kproma)               ,& !< surface height [m]
         lon_sp(kproma)              ,& !< longitude passed to sp
         lat_sp(kproma)              ,& !< latitude passed to sp
         z_fl_vr(kproma,klev)        ,& !< level height [m], vertically reversed indexing (1=lowest level)
         dz_vr(kproma,klev)          ,& !< level thickness [m], vertically reversed 
         sp_aod_vr(kproma,klev)      ,& !< simple plume aerosol optical depth, vertically reversed 
         sp_ssa_vr(kproma,klev)      ,& !< simple plume single scattering albedo, vertically reversed
         sp_asy_vr(kproma,klev)      ,& !< simple plume asymmetry factor, vertically reversed indexing
         sp_xcdnc(kproma)               !< drop number scale factor
    !
    ! --- 1.0 Calculate year fraction and proceed if in era covered by Simple Plume Model
    ! 
    CALL get_date_components(next_date, year=iyr)
    idy          = day_in_year(next_date)
    days_in_year = 365.0_dp - 28.0_dp + Get_JulianMonLen(iyr,2)
    year_fr      = iyr + REAL(idy,dp)/days_in_year
    
    IF (iyr > 1850.0_dp) THEN
      ! 
      ! --- 1.1 geographic information
      !
      DO jk=1,klev
        jki=klev-jk+1
        DO jl=1,kproma
          dz_vr  (jl,jk) = (vphysc%geohm1(jl,jki,krow)-vphysc%geohm1(jl,jki+1,krow))/grav
          z_fl_vr(jl,jk) = vphysc%geom1(jl,jki,krow)/grav
        END DO
      END DO
      z_sfc(:)  = geosp(1:kproma,krow)/grav
      lon_sp(:) = philon_2d(1:kproma,krow)
      lat_sp(:) = philat_2d(1:kproma,krow)
      ! 
      ! --- 1.2 Aerosol Shortwave properties
      !
      ! get aerosol optical properties in each band, and adjust effective radius
      !
      
      DO jwl = 1,nb_sw
        j_sw   = jpb1 + jwl - 1
        lambda = 1.e7_dp/ (0.5_dp * (sw_wv1(j_sw) + sw_wv2(j_sw)))
        CALL sp_aop_profile(                                                                  &
             klev               ,kproma             ,lambda              ,z_sfc(:)            , &
             lon_sp(:)          ,lat_sp(:)          ,year_fr             ,z_fl_vr(:,:)        , &
             dz_vr(:,:)         ,sp_xcdnc(:)        ,sp_aod_vr(:,:)      ,sp_ssa_vr(:,:)      , &
             sp_asy_vr(:,:)     ,plume_number       ,krow                                      )

        DO jk=1,klev
          DO jl=1,kproma
            asy_sw_vr(jl,jk,jwl) = asy_sw_vr(jl,jk,jwl) * ssa_sw_vr(jl,jk,jwl) * aod_sw_vr(jl,jk,jwl)    &
                 + sp_asy_vr(jl,jk)   * sp_ssa_vr(jl,jk)    * sp_aod_vr(jl,jk)
            ssa_sw_vr(jl,jk,jwl) = ssa_sw_vr(jl,jk,jwl) * aod_sw_vr(jl,jk,jwl)                           &
                 + sp_ssa_vr(jl,jk)   * sp_aod_vr(jl,jk)
            aod_sw_vr(jl,jk,jwl) = aod_sw_vr(jl,jk,jwl) + sp_aod_vr(jl,jk)
            asy_sw_vr(jl,jk,jwl) = MERGE(asy_sw_vr(jl,jk,jwl)/ssa_sw_vr(jl,jk,jwl),asy_sw_vr(jl,jk,jwl), &
                 ssa_sw_vr(jl,jk,jwl) > TINY(1.0_dp))
            ssa_sw_vr(jl,jk,jwl) = MERGE(ssa_sw_vr(jl,jk,jwl)/aod_sw_vr(jl,jk,jwl),ssa_sw_vr(jl,jk,jwl), &
                 aod_sw_vr(jl,jk,jwl) > TINY(1.0_dp))
          END DO
        END DO
      END DO
     
      DO jl=1,kproma
        x_cdnc(jl) = sp_xcdnc(jl)
      END DO
    END IF
    ! 
    ! --- 1.3 Aerosol longwave properties (presently blanck)
    !
  END SUBROUTINE add_aop_sp
END MODULE MO_SIMPLE_PLUMES
