      MODULE MAIN1
C***********************************************************************
C     MAIN1
C     AERMOD Model Data - Parameter, Variable and Array Declarations
C                         Global Data for All Modules
C
C***********************************************************************

      IMPLICIT NONE

c ----------------------------------------------------------------------
c --- AERMOD
c
c     Modified to accept additional building data required by
c     the PRIME building downwash modules
c
c     Changes are denoted in comment fields
c
c     Prepared by    Earth Tech, Inc
c     Prepared for   EPRI under contract WO3527-01
c ----------------------------------------------------------------------

C***********************************************************************
C     User Specified Model Parameters for Array Dimensions
C***********************************************************************

C --- Most array limits for data storage are now allocated at runtime.
      INTEGER, PARAMETER :: NWSCAT= 6, NKST= 6, NHR= 24,
     &                      NPAIR= 100, NHIANN= 10, 
     &                      NMXPM= 10, MXPLVL=50, MXGLVL=87

      INTEGER :: NYEARS

C**   NWSCAT = Number of Wind Speed Categories
C**   NKST   = Number of Stability Categories
C**   NHR    = Number of Hours in Met Data Loop
C**   NPAIR  = Number of Pairs of TXCONC and IDCONC for TOXXFILE Output
C**   NHIANN = Number of high period/annual averages to include in the
C**                   summary page of the output file (formerly controlled
C**                   by NVAL)
C**   NMXPM  = Number of high average High-N-High 24-hour PM2.5 averages or
C**                   High-N-High 1-hour NO2 or SO2 averages to include in the
C**                   summary table for PM-2.5/NO2/SO2 processing across years
C**   MXPLVL = Maximum number of levels in the observed profile file
C**   MXGLVL = Maximum number of levels in the gridded profiles (0 - 4000 m)
C**   NYEARS = Number of Years allowed for multi-year analyses for PM2.5,
C**                   1-hour NO2 and 1-hour SO2 standards, which are averaged
C**                   over the number of years modeled; a default value of 5
C**                   years is assumed, but the user can override the default
C**                   using the ME NUMYEARS keyword; e.g., setting NYEARS = 1
C**                   on the ME NUMYEARS keyword for 1 year of site-specific
C**                   met data will significantly reduce the memory requirements
C**                   for the MAXDCONT option.

C** The following array limits are set dynamically at runtime:
C**   NSRC   = Max Number of Sources
C**   NREC   = Max Number of Receptors
C**   NGRP   = Max Number of Source Groups
C**   NAVE   = Max Number of Short Term Averaging Periods
C**   NVAL   = Max Number of High Values by Receptor (RECTABLE Keyword)
C**   NTYP   = Max Number of Output Types per Run (CONC, DEPOS, DDEP and WDEP)
C**   NMAX   = Max Number of Overall Maximum Values (MAXTABLE Keyword)
C**   NSEC   = Number of Sectors for Building Dimensions
C**   NQF    = Number of Variable Emission Rate Factors Per Source
C**   NBF    = Number of Temporally Varying Background Concentrations
C**   NO3F   = Number of Temporally Varying Ozone Concentrations
C**   NPDMAX = Max Number of Particle Diameter Categories Per Source
C**   NVMAX  = Max Number of Vertices for AREA/AREACIRC/AREAPOLY and/or 
C**            OPENPIT Sources
C**   IXM    = Max Number of X-coord (Distance) Values Per Receptor Network
C**   IYM    = Max Number of Y-coord (Direction) Values Per Receptor Network
C**   NNET   = Max Number of Cartesian and/or Polar Receptor Networks
C**   NARC   = Maximum number of Receptor Groupings ('ARCs') for EVALCART
C**   NURB   = Maximum number of Urban Areas


C***********************************************************************
C     Programmer Specified Model Parameters
C***********************************************************************

      INTEGER, PARAMETER :: IFMAX=150, IKN=96, ISTRG=512, ILEN_FLD=200,
     &                      IERRN=289
C*#

C**   IFMAX  = Max Number of Fields Per Runstream Record
C**   IKN    = Number of Keywords
C**   ISTRG  = Max Length of Runstream Image Record
C**   ILEN_FLD = Max Length of Runstream Input Fields.  Also used 
C**              to specify length of input filenames and formats.
C**   IERRN  = Number of Error/Warning/Informational Messages


C***********************************************************************
C     Model Constants Specified as Parameters
C***********************************************************************

      DOUBLE PRECISION, PARAMETER :: 
     &                    G = 9.80616D0,       VONKAR = 0.4D0,
     &                    GOVRCP = 0.00977D0,  DCTODK = 273.16D0,
     &                    BETA1  = 0.6D0,      BETA2  = 0.4D0,
     &                    AT1PT2 = 1.2D0,      UMINGR = 0.01D0, 
     &                    GSIGV  = 0.073864D0, EFOLDH = 0.44D0,
     &                    SVUMIN = 0.05D0,     
     &                    SWMIN  = 0.02D0,     XVAL   = 0.0D0, 
     &                    SPTGMN = 0.002D0,    BSUBC  = 0.5D0,
     &                    SZCOEF = 2.15D0,     ALPHAR = 1.4D0,
     &                    LAMDAY = 2.3D0,      ASUBE  = 0.1D0,
     &                    REFPOP = 2.0D+6,     DELTRUR= 12.0D0,
     &                    RGAS   = 8.3145D0

C**   Declare SVMIN variable for minimum sigma-v, formerly specified as 
C**   PARAMETER SVMIN.  Default value of 0.2D0 is initialized in SUBROUTINE
C**   MODOPT; with SVMIN increased to 0.5D0 under the LOWWIND1 BETA option,
C**   and increased to 0.3D0 under the LOWWIND2 and LOWWIND3 options, with
C**   optional user inputs under the CO LOW_WIND keyword.
      DOUBLE PRECISION :: SVMIN

C**   Declare WSMIN variable for minimum wind speed, formerly set to 0.2828
C**   based on SQRT(2*SVmin*SVmin), where SVmin was 0.2. The default value 
C**   of WSMIN is also set to 0.2828, but may be adjusted under the LOWWIND
C**   BETA options using the CO LOW_WIND keyword.
      DOUBLE PRECISION :: WSMIN

C**   Declare FRANMAX variable for maximum meander factor, FRAN, used in
C**   the LowWind2 BETA option.  The "default" value for the LowWind2
C**   option is set at 0.95, but can be modified by the user under the
C**   optional LOW_WIND keyword on the CO pathway, within a range of 
C**   0.50 to 1.0, inclusive.
      DOUBLE PRECISION :: FRANMAX

C**   Set concentration unit conversion factors for NO2, SO2, CO, and 03
C**   for use with OZONEVAL, OZONEFIL, O3VALUES, and BACKGRND keywords;
C**   factors defined by pollutant for PPB-to-UG/M3 and PPM-to-UG/M3,
C**   based on reference temperature (25 C) and pressure (1013.25 mb).
C**   Note that factors for NO2 and SO2 are PPB/(UG/M3) and PPM/(UG/M3),
C**   and factors for CO and O3 are for (UG/M3)/PPB and (UG/M3)/PPM.
      DOUBLE PRECISION, PARAMETER ::  
     &                    NO2_PPB = 0.5319D0, NO2_PPM = 0.5319D-3,
     &                    SO2_PPB = 0.3823D0, SO2_PPM = 0.3823D-3,
     &                    CO_PPB  = 1.144D0,  CO_PPM  = 1.144D3,
     &                    O3_PPB  = 1.960D0,  O3_PPM  = 1.960D3

           
      DOUBLE PRECISION :: PI, TWOPI, RTOFPI, SRT2PI, RTOF2, RTPIBY2, 
     &                    DTORAD, RTODEG, 
     &                    THIRD

C**   PI     = PI               ! Initialized in sub.VARINI as 4.0D0*DATAN(1.0D0)
C**   TWOPI  = 2.*PI            ! Initialized in sub.VARINI
C**   RTOFPI = SQRT(PI)         ! Initialized in sub.VARINI
C**   SRT2PI = SQRT(2.*PI)      ! Initialized in sub.VARINI
C**   RTOF2  = SQRT(2.)         ! Initialized in sub.VARINI
C**   RTPIBY2= SQRT(PI/2.)      ! Initialized in sub.VARINI
C**   DTORAD = Degrees to Radians Conversion   ! Initialized in sub.VARINI as PI/180.0D0
C**   RTODEG = Radians to Degrees Conversion   ! Initialized in sub.VARINI as 180.0D0/PI
C**   THIRD  = 1.0/3.0          ! Initialized in sub.VARINI

C**   G      = Acceleration Due to Gravity (m/s**2)
C**   VONKAR = von Karman constant
C**   GOVRCP = Gravity divided by specific heat at constant pressure
C**   DCTODK = Degrees Celsius to kelvin conversion factor
C**   BETA1  = Coeff. in the calculation of 'direct' plume rise
C**   BETA2  = Coeff. in the calculation of buoyancy-induced dispersion
C**   AT1PT2 = The fraction of the mixed layer and above the mixed
C              layer through which a variable changes its value
C**   UMINGR = Minimum value for a gridded wind speed
C**   GSIGV  = Constant used in converting sigma_A to sigma_V
C**   EFOLDH = Constant in computation of dTHETA/dZ in stable atmosphere
C**   SVUMIN = Minimum value applied to Sigma_V / U when calculating
C**            Sigma_Y.
C**   SVMIN  = Minimum value applied to measured Sigma_V values
C**   SWMIN  = Minimum value applied to measured Sigma_W values
C**   XVAL   = Gradient in the miXing layer (unstable atmosphere)
C**   SPTGMN = Minimum vert. potential temp. gradient for stable atmosphere
C**   BSUBC  = Constant used in computing sigma_Z for surface layer releases
C**   SZCOEF = Coefficient of sigma-z to define the plume half-width
C**   ALPHAR = Parameter used in January 18, 1995 indirect source MCA
C**   LAMDAY = Parameter used in January 18, 1995 indirect source MCA
C**   ASUBE  = Parameter used in January 18, 1995 indirect source MCA
C**   REFPOP = Reference population for urban option (4/1/96 MCA)
C**   DELTRUR= Surface cooling in the rural area (4/1/96 MCA)
C**   RGAS   = ideal gas law constant = 8.3145 Pa-m^3/mol-K
C**

C**   NUMSYEFF = Number of effective sigma-y's to use with the 
C**              FASTALL non-DFAULT option; 
C**              receptors more than NUMSYEFF*SYEFF off the plume
C**              centerline are skipped to optimize runtime
      DOUBLE PRECISION, PARAMETER :: NUMSYEFF  = 4.0D0

C**   NUMSYEFF3= Number of effective sigma-y's to use with the 
C**              LowWind3 non-DFAULT/BETA option; 
C**              receptors more than NUMSYEFF*SYEFF off the plume
C**              centerline are skipped to optimize runtime
      DOUBLE PRECISION, PARAMETER :: NUMSYEFF3 = 6.0D0

C**   MAXDIST = Maximum transport distance for calculation;
C**             set to 80km for FASTALL or FASTAREA options;
C**             this was formerly associated with the TOXICS
C**             option, which is now obsolete.
C**             Set to 1.0D20 for applications w/o FASTALL or FASTAREA
      DOUBLE PRECISION :: MAXDIST


C***********************************************************************
C     Common Block for Input/Output File Units (Initialized in BLOCK DATA)
C***********************************************************************

      INTEGER :: INUNIT, IOUNIT, MFUNIT, MPUNIT, IERUNT, IERWRT,
     &           IDPUNT, IDPUN2, IRSUNT, IEVUNT, ITEVUT, IHREMI,
     &           IBGUNT(6), IO3UNT(6), INCUNT, ISUMUNT, DBGUNT, DBMUNT, 
     &           AREADBUNT, GDEPDBG, PDEPDBG, PRMDBUNT, PVMDBG, OLMDBG, 
     &           ARMDBG, RDISPUNT

C**   These input/output file units are initialized below in a DATA statement
C**   INUNIT = Input Runstream File Unit (Initialized to 7)
C**   IOUNIT = Main Printed Output File Unit (Initialized to 8)
C**   PVMDBG = PVMRM debug option file unit (initialized to 9)
C**   OLMDBG = OLM debug option file unit (initialized to 9)
C**   ARMDBG = ARM/ARM2 debug option file unit (initialized to 9)
C**   IERUNT = Temporary Error/Message File Unit (Initialized to 10)
C**   IERWRT = Permanent Detailed Error/Message File Unit (Init. to 11)
C**   IDPUNT = Main SAVEFILE Unit for Re-start Option (Init. to 12)
C**   IDPUN2 = Secondary SAVEFILE Unit for Re-start Option (Init. to 14)
C**   IRSUNT = INITFILE Unit for Re-start Option (Initialized to 15)
C**   IHREMI = Hourly Emission Parameters File Unit (Init. to 16)
C**   IEVUNT = Event File Unit for Use With EVENT Model Option (Init. to 17)
C**   ITEVUT = Temporary Event File Used to Store High Value Events for
C**            Summary Tables and for EVENTFIL Option (Initialized to 18)
C**   MFUNIT = Input Surface Met Data File Unit (Initialized to 19)
C**   INCUNT = INCLUDED File Unit (Initialized to 20)
C**   MPUNIT = Input Profile Met Data File Unit (Initialized to 21)
C**   ISUNIT = Surface Meteorology File for SCIM'd data (Initialized to 22)
C**   IPUNIT = Profile Meteorology File for SCIM'd data (Initialized to 23)
C**   DBGUNT = Debug Output File for Calculations (Init. to 24)
C**   DBMUNT = Debug Output File for Meteorology Profiles (Init. to 25)
C**   AREADBUNT = Debug Output File for AREA/LINE/OPENPIT sources (Init. to 26)
C**   PRMDBUNT  = Debug Output File for PRIME related Debug Infor (Init. to 27)
C**   ISUMUNT = Summary File under SUMMFILE Option (Init. to 28)
C**   GDEPDBG = Debug Output File for Gas Deposition Velocities (Init. to 29)
C**   PDEPDBG = Debug Output File for Particle Deposition Velocities (Init. to. 30)
C**   IBGUNT  = Hourly Background Data File(s) by Sector (= 2000 + sector number)
C**   IO3UNT  = Hourly Ozone Data File(s) for PVMRM and OLM Options (= 1000 + sector number)

C***********************************************************************
C     This is The Global Variable Definition Block for Runstream Data
C***********************************************************************

      LOGICAL BLINE, INFLD, MARK, ECHO

      CHARACTER PATH*2, PPATH*2, KEYWRD*8, PKEYWD*8, KEYWD*8, KTYPE*5,
     &          RUNST*1

      CHARACTER (LEN=ILEN_FLD) :: FIELD, INPFIL, OUTFIL, INCFIL
      CHARACTER (LEN=ISTRG)    :: RUNST1

      INTEGER ::  LOCB(IFMAX), LOCE(IFMAX), IFC, IDC1, IPNUM, IPPNUM
      DIMENSION   FIELD(IFMAX), KEYWD(IKN), RUNST(ISTRG)


C***********************************************************************
C     This is The Global Variable Definition Block for Error Handling
C***********************************************************************

      LOGICAL FATAL, ISTART, IFINIS, RECERR, ERRLST, EOF, ALLOC_ERR
      LOGICAL L_SkipMessages
      
      REAL    :: STORE                ! Estimate of memory storage requirement

      CHARACTER ERRMSG*50, ERRCOD*3, VERSN*6
      CHARACTER (LEN=6) :: C_METVER         ! Character string for met version
      CHARACTER (LEN=ILEN_FLD) :: MSGFIL

      DIMENSION  ERRMSG(IERRN), ERRCOD(IERRN)
      INTEGER :: ILINE, IQLINE, IBLINE, IOLINE, IERROR, IFTL, IWRN,
     &           INFO, ICLM, IMSG, NFATAL, NWARN, IPAGE, IPGSUM
C --- Met data array indices for use with MAXDCONT option
      INTEGER :: IHR_NDX, IYR_NDX
      DOUBLE PRECISION :: EXPLIM
      INTEGER :: ICSTAT(50), ISSTAT(50), IRSTAT(50), IMSTAT(50),
     &           IOSTAT(50), IESTAT(50)
      INTEGER :: INCSET, IXYSET, IEVSET, IHLSET, IFGSET


C***********************************************************************
C     This is The Global Variable Definition Block for COntrol Pathway
C***********************************************************************

      LOGICAL DFAULT, CONC, DEPOS, DDEP, WDEP, RURAL, URBAN, GRDRIS,
     &        NOSTD, NOBID, CLMPRO, MSGPRO, PERIOD, ANNUAL, MONTH,
     &        FLAT, ELEV, FLATSRCS, FLGPOL, RUN, EVENTS, RSTSAV, 
     &        RSTINP, DAYTAB, MXFILE, PPFILE, PLFILE, ANPOST, ANPLOT, 
     &        STATOK, MULTYR, TXFILE, RKFILE, SEASONHR, 
     &        MXDAILY, MXDAILY_BYYR, L_MAXDCONT,
     &        DDPLETE, WDPLETE, DRYDPLT, WETDPLT, NODRYDPLT, NOWETDPLT,
     &        FSTCMP, EVONLY, SOCONT, DETAIL, NEWMET, ARDPLETE, 
     &        PM25AVE, NO2AVE, SO2AVE, L_NO_PM25AVE, L_NO_NO2AVE, 
     &        L_NO_SO2AVE, NOCHKD, NOWARN, 
     &        DEBUG, METEORDBG, AREADBG, PRIMEDBG, PVMRMDBG, OLMDEBUG, 
     &        ARMDEBUG, ARM2DEBUG, DEPOSDBG,  
     &        L_WARNCHKD, SCIM, SCIMHR, 
     &        FASTAREA, FASTALL, L_NonDFAULT, 
     &        SCREEN, URBSTAB, PRM_FSTREC, ROMBERG,
     &        PVMRM, PVMRM2, PSDCREDIT, OLM, L_MULTURB, 
     &        L_PRESET_URBAN, L_UrbanTransition, L_URBAN_ALL,
     &        L_Urban, L_Rural,
     &        BETA, ARM, ARM2
C*#
C --- Logical variables used to track inconsistency between 'original'
C     results and results calculated during the MAXDCONT option 
C     internal "post-processing" and between 'original' results and
C     results calculated during the EVENT post-processing 
      LOGICAL L_EVENT_OrigConc_Warning, L_MAXDCONT_OrigConc_Warning

      LOGICAL :: L_EFFSIGY,   L_VECTORWS, 
     &           L_AdjUstar,  L_BULKRN,
     &           L_MMIF_Data, L_MMIF_Profile,
     &           L_LowWind1,  L_LowWind2, L_LowWind3,
     &           L_UserSVmin, L_UserWSmin, L_UserFRANmax,
     &           L_CCVR_Sub,  L_TEMP_Sub
      
      CHARACTER (LEN=ILEN_FLD) :: TITLE1, TITLE2
      CHARACTER RUNDAT*8, RUNTIM*8
      CHARACTER EVPARM*6, CHRAVE*5, CHIDEP*4, SOELEV*6, REELEV*6, 
     &          TGELEV*6, OUTTYP*5, NO2_FIELD4*3, SO2_FIELD4*3,
     &          PM25_FIELD4*3

      CHARACTER (LEN=ILEN_FLD) :: SAVFIL, SAVFL2, INIFIL, EVFILE,
     &                            DBGFIL, DBMFIL, DBAREAFIL, DBPVFIL, 
     &                                                       RDISPFIL,
     &                            DBOLMFIL, DBPRMFIL, DBARMFIL, 
     &                            DBARM2FIL, OZONFL(6), O3FILUNITS, 
     &                            O3VALUNITS, O3FORM(6), OzoneUnits, 
     &                            URBNAM

      DOUBLE PRECISION ::  O3CONC, O3BACK(6), NO2Equil, NO2Stack,
     &                     EV_O3CONC(NHR), O3SECT(6), O3_Max24hr(NHR,6),
     &                     ARM_1hr, ARM_Ann, ARM2_Min, ARM2_Max,
     &                     RatioARM2

      LOGICAL ::  O3MISS, L_O3Sector, L_O3Hourly
      LOGICAL ::  L_O3File(6), L_O3VAL(6), L_O3VALUES(6)
      LOGICAL ::  L_AO3MISS(24)

      INTEGER ::  NUMO3Sects, IO3SECT 
      INTEGER ::  NHIVAL, NMXVAL, NDUMP, NHIMXDLY

C --- Declare array for O3SECTs by hour for EVENT processing
      INTEGER ::  AO3SECT(NHR)

      INTEGER ::  NSRC, NREC, NGRP, NAVE, NVAL, NTYP, NMAX,
     &            NSEC, NQF, NBF, NO3F, NPDMAX, NNET, IXM, IYM,
     &            NEVE, NUMEVE, IEVENT, NARC, NOLM, NURB, NPSD,
     &            IO3MAX(6), IBGMAX(6) 

      INTEGER ::  NUMCONT     ! Number of contributing sources for PVMRM/PVMRM2
      
      INTEGER, ALLOCATABLE :: KAVE(:)

      LOGICAL, ALLOCATABLE :: EVAL(:)

      ALLOCATABLE ::  CHRAVE(:), CHIDEP(:,:), OUTTYP(:), URBNAM(:)
      CHARACTER (LEN=9)   :: MODOPS(24)
      CHARACTER (LEN=240) :: MODOPS_String


C***********************************************************************
C     This is The Global Variable Definition Block for SOurce Pathway
C***********************************************************************

      CHARACTER SRCID*12, SRCTYP*8, SOPCRD*1, SOGAS*1, URBSRC*1, 
     &          GRPID*8, EMILBL*40, OUTLBL*40, POLLUT*8, 
     &          QFLAG*8, BFLAG(6)*8, O3FLAG(6)*8, PERLBL*40, OLMID*8, 
     &          URBID*8, PSDID*8 

      CHARACTER (LEN=ILEN_FLD) :: HRFILE, BKGRND_File(6), BGFORM(6), 
     &                            BackUnits
C*#
      CHARACTER PREVSRCID*12
      CHARACTER PREVGRPID*8

C --- Declare logicals related to deposition options applicable to output types;
C     LDPART indicates that particle dry deposition is used
C     LWPART indicates that particle wet deposition is used
C     LDGAS  indicates that gaseous dry deposition is used
C     LWGAS  indicates that gaseous wet deposition is used
      LOGICAL LDPART, LWPART, LDGAS, LWGAS

C --- Declare logicals related to background concentration options:
C     L_BACKGRND indicates generally that background concentration options are used
C     L_BGHourly indicates generally that hourly background concentrations are used
C     L_BGSector indicates generally that sector-varying background concs are used
      LOGICAL L_BACKGRND, L_BGHourly, L_BGSector
C --- Declare logicals related to sector-varying background concentration options:
C     L_BGFile indicates that hourly background concentrations are used for given sector
C     L_BGValues indicates that non-hourly background concs are available for given sector
      LOGICAL L_BGFile(6), L_BGValues(6)
      
      LOGICAL, ALLOCATABLE :: GRP_BACK(:)

      DOUBLE PRECISION :: BGCONC, EV_BGCONC(NHR)
      DOUBLE PRECISION :: BGBACK(6), BGSECT(6)

      INTEGER ::  NUMBGSects, IBGSECT
      INTEGER ::  NSubBGHOUR                      ! No. of BGHOUR subs
      INTEGER ::  INDX_GRPALL                     ! Index for SrcGroup ALL for ARM/ARM2 options 

C --- Declare array for BGSECTs by hour for EVENT processing
      INTEGER ::  ABGSECT(NHR)

      DOUBLE PRECISION, ALLOCATABLE ::  AXS(:), AYS(:), AZS(:), AQS(:),
     &                        AHS(:), ATS(:), AVS(:), ADS(:),
     &                        ASYINI(:), ASZINI(:), ANO2_RATIO(:),
     &                        ADSFACT(:)
      DOUBLE PRECISION, ALLOCATABLE :: AAQS(:,:,:), AAHS(:,:,:), 
     &                        AATS(:,:,:), AAVS(:,:,:), 
     &                        AASYINI(:,:,:), AASZINI(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE ::  ADSBH(:,:), ADSBW(:,:),
     &                        ADSBL(:,:), ADSXADJ(:,:), ADSYADJ(:,:)

C*---	Added for LINE source
      DOUBLE PRECISION, ALLOCATABLE ::  AWIDTH(:), AXS1(:), AYS1(:),
     &                                             AXS2(:), AYS2(:) 
C*---

      INTEGER, ALLOCATABLE :: INPD(:), NDXSTK(:)
      
      DOUBLE PRECISION, ALLOCATABLE ::  QFACT(:,:), 
     &                                  O3VARY(:,:),
     &                                  BACKGRND(:,:)
      DOUBLE PRECISION :: EMICON, HAFLIF, DECOEF, VAIRMS, ZRDEP, VDPHOR
      DOUBLE PRECISION, ALLOCATABLE :: EMIFAC(:), APDIAM(:,:), 
     &                                 APHI(:,:), APDENS(:,:),
     &                                 AVGRAV(:,:), ATSTOP(:,:)

C*--- Variables for hourly emissions 
      DOUBLE PRECISION ::  HRQS, HRTS, HRVS, HRHS, HRSY, HRSZ
      DOUBLE PRECISION, ALLOCATABLE ::  EV_HRQS(:,:), EV_HRTS(:,:), 
     &                                  EV_HRVS(:,:), EV_HRHS(:,:),
     &                                  EV_HRSY(:,:), EV_HRSZ(:,:)

      INTEGER ::  FULLHRQ
C*----
C*#
      CHARACTER :: PSDSRCTYP*2
      INTEGER, ALLOCATABLE :: IGROUP(:,:), IGRP_OLM(:,:), IGRP_PSD(:,:),
     &                        IURBGRP(:,:)

      ALLOCATABLE ::  SRCID(:), SRCTYP(:), SOPCRD(:),
     &                SOGAS(:), URBSRC(:),
     &                GRPID(:), QFLAG(:), EMILBL(:),
     &                OUTLBL(:),
     &                PERLBL(:), OLMID(:), PSDID(:), URBID(:)
      LOGICAL, ALLOCATABLE :: L_OLMGRP(:), L_PSDGRP(:)
      LOGICAL, ALLOCATABLE :: L_HRLYSIG(:), L_FLATSRC(:),
     &                        L_WakeMessage(:)

      ALLOCATABLE :: PSDSRCTYP(:)


C**   NDXSTK   Index of the gridded height immediately below stack top
C**   FULLHRQ  Date/Time Stamp for Hourly Emissions Data
C**   HRFILE   Hourly Emissions Data File Name
C*#

C***********************************************************************
C     This is The Global Variable Definition Block for the New Area
C     Source Algorithm - 7/7/93
C
C*    Added XVERT,YVERT - Jayant Hardikar, PES, 7/20/94
C***********************************************************************

C**   NVMAX and NVMAX2 are now dynamically allocated at runtime for
C**   AREAPOLY and AREACIRC sources.  If AREACIRC sources are included, 
C**   then the maximum number of vertices needed for AREACIRC sources
C**   is based on the maximum number specified by the user using the
C**   optional NVERTS parameter + 4.  The default number of vertices 
C**   for an AREACIRC source is still set at 20 vertices, but that 
C**   value is only used to determine array storage if the the input 
C**   file includes an AREACIRC source without the NVERTS parameter.
C**   If AREAPOLY sources are included, then NVMAX is set to the maximum 
C**   number vertices specified for an AREAPOLY source + 8 (but not 
C**   less than the number needed for AREACIRC sources).  The '+ 8' 
C**   allows for the additional number of sides on an AREAPOLY source that
C**   may be needed to define the portion of the source being integrated,
C**   depending on the source/receptor/wind direction geometry.  This 
C**   allows for some complexity in the shape of an AREAPOLY source, but 
C**   very complex shapes may result in runtime errors (E406) and should 
C**   be avoided. 
C**   If no AREACIRC or AREAPOLY sources are included, but rectangular
C**   AREA sources and/or OPENPIT sources are included, then the value
C**   of NVMAX is set to 8 (4 sides plus 4).

C**   NVMAX2= NVMAX * 2

C**   Assign NVPOLY parameter for the number of additional vertices to 
C**   include for AREAPOLY sources, beyond the number of vertices used
C**   to define the source.  This is needed to account for the integrated
C**   portion of the area intersecting several sides of the polygon,
C**   depending on the geometry of the source, receptor, and wind direcion.
C**   The NVPOLY parameter is assigned a value of 12, which should work 
C**   in most cases.  NVPOLY is added to NVMAX in subroutine SRCSIZ in 
C**   the 'aermod.f' source file.
      INTEGER, PARAMETER :: NVPOLY = 12

      INTEGER :: IVERT, NVERT, NSEGS,
     &           NVMAX, NVMAX2, NPIT, NPNT, NVOL, NLINE, NBLINE,
     &           NAREA, NPOLY, NVTEMP, NCIRC, 
     &           NPTEMP
      DOUBLE PRECISION, ALLOCATABLE :: UVERT(:), VVERT(:), VNVERT(:),
     &                                 WVERT(:), UASEGS(:), UBSEGS(:),
     &                                 XVERT(:), YVERT(:)
      DOUBLE PRECISION, ALLOCATABLE :: SPA(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: AXINIT(:), AYINIT(:), AANGLE(:),
     &                                 AXVERT(:,:), AYVERT(:,:),
     &                                 RADIUS(:), AXCNTR(:), AYCNTR(:)
      INTEGER, ALLOCATABLE :: NVERTS(:)

      LOGICAL LSEG


C***********************************************************************
C     This is The Global Variable Definition Block for the New OPENPIT
C     Source Algorithm - 7/19/94
C***********************************************************************

      DOUBLE PRECISION, PARAMETER   :: ALPHA = 0.029D0
      DOUBLE PRECISION, ALLOCATABLE :: AALPHA(:), APDEFF(:), AVOLUM(:)
      DOUBLE PRECISION, ALLOCATABLE :: EFRAC(:), QPART(:)
      DOUBLE PRECISION :: PALPHA, THETA, PDEFF, PDREL, PITFRA, QEFF
      DOUBLE PRECISION :: PITLEN, PITWID, PITL, EMIHGT, XEFF, YEFF


C*    ALPHA     = Proportionality Constant for determining Escape Fraction
C*    AALPHA    = Array of Angles of Long Pit Dimension w.r.t.
C*                North for Each Source
C*    PALPHA    = Angle of Long Pit Dimension w.r.t. North for
C*                the Current Source
C*    THETA     = Wind Direction angle w.r.t Long Axis
C*                of the Pit
C*    APDEFF    = Array of Effective Pit Depths for Each Source
C*    PDEFF     = Effective Pit Depths for Current Source
C*    PDREL     = Relative Pit Depth
C*    AVOLUM    = Array of Volume of the OPENPIT Sources
C*    EFRAC     = Array of Escape Fractions
C*    QPART     = Array of Adjusted Emission Rates
C*    PITFRA    = Fractional Size of the Effective Pit Area
C*    PITLEN    = Length of the Pit
C*    PITWID    = Width of the Pit
C*    PITL      = Along-Wind Length of the Pit
C*    EMIHGT    = Height of Emissions Above Base of Pit
C*    XEFF      = X-dimension of Effective Pit
C*    YEFF      = Y-dimension of Effective Pit



c----------------------------------------------------------------------
c --- COMMON BLOCK /DRYGAS/ -- Dry deposition parameters        CALPUFF
c----------------------------------------------------------------------
c
      DOUBLE PRECISION, ALLOCATABLE :: pdiff(:),pdiffw(:),rmolwt(:),
     &                                 alphas(:),react(:),henry(:),
     &                                 rcli(:),finemass(:), scf(:)
     
      LOGICAL, ALLOCATABLE :: L_METHOD2(:)
      
      INTEGER :: ISEAS_GD(12), ILAND_GD(36), NCLOUD
      
      DOUBLE PRECISION :: rm, rcut, qsw, xlai, vdepg, uservd, zsubp, 
     &                    delta_z, FO, FSEAS2, FSEAS5, fracsat, 
     &                    liqcont, denom, xnu

      DOUBLE PRECISION :: Wold, Wnew, f2, EsTa

      CHARACTER (LEN = 40) ::  REFSPE

      LOGICAL  LUSERVD
C
C     REFSPE      - Reference Species (Default is SO2)
c
c --- COMMON BLOCK /DRYGAS/ Variables:
c       PDIFF(NSRC) - real    - Molecular diffusivity (m**2/s)
c                               of each pollutant.
c                               SEE NOTE #1
c      PDIFFW(NSRC) - real    - Molecular diffusivity in water (m**2/s)
c                               of each pollutant.
c                               SEE NOTE #1
c      RMOLWT(NSRC) - real    - Molecular weight of pollutant (g/mol)
c                               of each pollutant.
c      ALPHAS(NSRC) - real    - Solubility enhancement factor due
c                               to the aqueous phase reactivity of
c                               the pollutant.
c       REACT(NSRC) - real    - Reactivity factor for each
c                               pollutant.
c                RM - real    - Mesophyll resistance (s/m)
c                               SEE NOTE #2
c       HENRY(NSRC) - real    - Henry's law constant (ratio of
c                               gas to aqueous phase concentration
c                               of the pollutant).
c              RCUT - real    - Cuticle resistance (s/m).
c
c  NOTE #1: Input units of this variable are cm**2/s.  Conversion to m**2/s
c           is made internally in the SETUP phase.
c
c  NOTE #2: Input units of s/cm are converted to s/m in the SETUP phase.



C***********************************************************************
C     This is The Global Variable Definition Block for REceptor Pathway
C***********************************************************************

      LOGICAL ISTA, IEND, NEWID

      CHARACTER NETID*8, NETIDT*8, PNETID*8, NTID*8, NTTYP*8,
     &          RECTYP*2, PXSOID*12, PESOID*12, ARCID*8

      DOUBLE PRECISION, ALLOCATABLE ::  AXR(:), AYR(:), AZELEV(:),
     &                                  AZFLAG(:), AZHILL(:)
      INTEGER, ALLOCATABLE :: IREF(:), NDXARC(:)
      ALLOCATABLE ::  NETID(:), RECTYP(:), NTID(:),
     &                NTTYP(:), ARCID(:)
      INTEGER ::      ICOUNT, JCOUNT, IZE, IZH, IZF, IRZE, IRZH, IRZF,
     &                IRXR, IRYR, IRHZ, INNET
      DOUBLE PRECISION ::  XINT, YINT
      DOUBLE PRECISION, ALLOCATABLE ::  XCOORD(:,:), YCOORD(:,:),
     &                                  XORIG(:), YORIG(:)
      INTEGER, ALLOCATABLE :: NETSTA(:), NETEND(:),
     &                        NUMXPT(:), NUMYPT(:)

C**  AZHILL Hill Height Associated with the Receptor
C**  HCRIT  Critical dividing streamline associated with the receptor


C***********************************************************************
C     This is The Global Variable Definition Block for MEteorology Pathway
C***********************************************************************

      CHARACTER SFNAME*40, UANAME*40, ONNAME*40, ALAT*10, ALON*10

      CHARACTER (LEN=ILEN_FLD) :: METINP, SCIM_SFCFIL, SCIM_PROFIL,
     &                            PROINP
      CHARACTER (LEN=ILEN_FLD) :: METFRM, PROFRM
      CHARACTER :: MMIF_Version*27
      
      LOGICAL SCIMOUT
      
C**** Logical flags for met data version, L_OldMetVer is used to flag
C     an outdated met version date in the surface file header record;
C     L_NAD_ADJ_Flags is used to flag cases where surface file header shows
C     current met version date, but surface file lacks additional fields
C     introduced with version 11059 for the wind data source/adjustment.
      LOGICAL :: L_OldMetVer, L_NAD_ADJ_Flags

C**** Include logical variable to track whether current year being
C     processed is a leap year or not, in order to properly handle
C     DAYRANGE inputs based on MN/DY for multiple-year met data
      LOGICAL :: L_LeapYear

      INTEGER :: ISDATE, IEDATE, ISYR, ISMN, ISDY, ISHR, IEYR,
     &           IEMN, IEDY, IEHR, IPROC(366), IPROCL(366), 
     &           ISYEAR, IUYEAR, IOYEAR,
     &           IDSURF, IDUAIR, IDSITE, ISJDAY, IEJDAY,
     &           NDAYS, INCRST,
     &           ISTRT_CENT, ISTRT_WIND,
C        RWB/MJ - allow for SCIM option - May, 1998.
     &           NREGSTART, NREGINT, IFIRSTHR, ISUNIT, IPUNIT,
     &           NSKIPTOT, IMETMSG

      INTEGER :: FULL_YYMMDD, IEDATE_YYMMDD

      DOUBLE PRECISION :: UCAT(5), ROTANG,
     &                    VIRTPNT_URB(NKST), VIRTPNT_RUR(NKST), VP_FACT
      DOUBLE PRECISION :: SFX, SFY, UAX, UAY, ONX, ONY
   

C***********************************************************************
C     This is The Global Variable Definition Block for METEXT
C***********************************************************************

      LOGICAL CLMHR, MSGHR, UNSTAB, NEUTRL, STABLE,
     &        RUNERR, PFLERR, ENDMON, METHDR,
     &        HOURLY, L_DayOfWeekOpts

      LOGICAL, ALLOCATABLE :: L_MorningTrans(:), AL_MorningTrans(:,:,:),
     &                        ACLMHR(:,:), AMSGHR(:,:),
     &                        ASTABLE(:,:), AUNSTAB(:,:),
     &                        AURBSTAB(:,:)
      
      INTEGER ::  KSTMSG
      INTEGER ::  IHOUR, IYEAR, IMONTH, IDAY, KURDAT, JDAY, ISEAS,
     &            KHOUR, KYEAR, KMONTH, KDAY, KURPFL, NTOTHRS,
     &            IPHOUR, IPYEAR, IPDATE, IPCODE, KST,
     &            IYR, IDAY_OF_WEEK, IDAY_OF_WEEK7, NPLVLS, NTGLVL,
     &            IFLAG(MXPLVL)
      INTEGER ::  JDAY_PREV
      INTEGER ::  FULLDATE
      DOUBLE PRECISION ::  SFCHF, USTAR, WSTAR, VPTGZI, ZICONV, ZIMECH,
     &                     OBULEN, SFCZ0, BOWEN, ALBEDO, UREF, WDREF,
     &                     UREFHT, TA, TREFHT, ZI, AFV,
     &                     BVF, BVPRIM, XLAT, TSIGN, ZIRUR,
     &                     PRATE, PREC1, PREC2, TOTAL_PRECIP,
     &                     UREF10, RURUSTR, RUROBULEN, RH, SFCP

      DOUBLE PRECISION, ALLOCATABLE :: URBPOP(:), URBZ0(:), ZIURB(:), 
     &                     URBWSTR(:),
     &                     URBUSTR(:), URBOBULEN(:)

      INTEGER, ALLOCATABLE :: IKST(:,:), IAPCODE(:,:), NACLOUD(:,:)
      
      DOUBLE PRECISION, ALLOCATABLE :: APRATE(:,:), AQSW(:,:),ARH(:,:), 
     &                                 ASFCP(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: ASFCHF(:,:), AUREF(:,:), 
     &                                 AUREFHT(:,:), ATA(:,:), 
     &                                 ATREFHT(:,:), AWDREF(:,:), 
     &                                 AUSTAR(:,:), AWSTAR(:,:),
     &                                 AZICONV(:,:), AZIMECH(:,:), 
     &                                 AOBULEN(:,:), AVPTGZI(:,:),
     &                                 ASFCZ0(:,:), ABOWEN(:,:),
     &                                 AALBEDO(:,:), AWNEW(:,:), 
     &                                 AWOLD(:,:), AESTA(:,:), 
     &                                 AF2(:,:), APREC1(:,:), 
     &                                 APREC2(:,:),
     &                                 AO3CONC(:,:), ABGCONC(:,:),
     &                                 AKST(:,:),ABLTA(:,:)

      INTEGER :: ISTRHOUR, ISTRDY, ISTRMN, 
     &           IENDHOUR, IENDDY, IENDMN, NUMYRS, NREMAIN, NDX4ZI

      INTEGER, ALLOCATABLE :: ANDX4ZI(:,:)

      DOUBLE PRECISION :: PFLHT(MXPLVL), PFLWD(MXPLVL), PFLWS(MXPLVL),
     &                    PFLTA(MXPLVL), PFLSA(MXPLVL), PFLSW(MXPLVL),
     &                    PFLSV(MXPLVL), PFLTG(MXPLVL), PFLTGZ(MXPLVL)
     
      DOUBLE PRECISION, ALLOCATABLE :: APFLHT(:,:,:), APFLWD(:,:,:),
     &                                 APFLWS(:,:,:), APFLTA(:,:,:),
     &                                 APFLSA(:,:,:), APFLSW(:,:,:),
     &                                 APFLSV(:,:,:), APFLTG(:,:,:),
     &                                 APFLTGZ(:,:,:)

      INTEGER, ALLOCATABLE :: AIFLAG(:,:,:)
      
      INTEGER, ALLOCATABLE :: ANPLVLS(:,:), ANTGLVL(:,:)
      
      DOUBLE PRECISION :: GRIDHT(MXGLVL), GRIDWD(MXGLVL), 
     &                    GRIDWS(MXGLVL), GRIDSW(MXGLVL), 
     &                    GRIDSV(MXGLVL), GRIDTG(MXGLVL),
     &                    GRIDPT(MXGLVL),
C---  Add density profile for PRIME
     &                    GRIDRHO(MXGLVL),
C---  Add tubulence dissipation rate (epsilon) profile for PVMRM/PVMRM2
     &                    GRIDEPS(MXGLVL)
      DOUBLE PRECISION :: GRDSWR(MXGLVL), GRDSVR(MXGLVL), 
     &                    GRDTGR(MXGLVL),
     &                    GRDPTR(MXGLVL)

      DOUBLE PRECISION, ALLOCATABLE :: GRDSWU(:,:), GRDSVU(:,:), 
     &                                 GRDTGU(:,:), GRDPTU(:,:)

      DOUBLE PRECISION, ALLOCATABLE :: 
     &                    AGRIDHT(:,:,:), AGRIDWD(:,:,:), 
     &                    AGRIDWS(:,:,:), AGRIDSW(:,:,:), 
     &                    AGRIDSV(:,:,:), AGRIDTG(:,:,:),
     &                    AGRIDPT(:,:,:),
C---  Add density profile for PRIME
     &                    AGRIDRHO(:,:,:),
C---  Add tubulence dissipation rate (epsilon) profile for PVMRM/PVMRM2
     &                    AGRIDEPS(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: AGRDSWR(:,:,:), AGRDSVR(:,:,:), 
     &                                 AGRDTGR(:,:,:), AGRDPTR(:,:,:),
     &                                 AUATZI(:,:),  ASVATZI(:,:),
     &                                 ASWATZI(:,:), APTATZI(:,:),
     &                                 AUAVG(:,:),  ASVAVG(:,:),
     &                                 ASWAVG(:,:), APTAVG(:,:)

      DOUBLE PRECISION, ALLOCATABLE :: AGRDSWU(:,:,:,:), 
     &                                 AGRDSVU(:,:,:,:),
     &                                 AGRDTGU(:,:,:,:), 
     &                                 AGRDPTU(:,:,:,:),
     &                                 AZIURB(:,:,:),
     &                                 AURBWSTR(:,:,:),
     &                                 AURBUSTR(:,:,:),
     &                                 AURBOBULEN(:,:,:),
     &                                 ARURUSTR(:,:), 
     &                                 ARUROBULEN(:,:)

      DOUBLE PRECISION :: TG4PFL, TG4XTR,
     &                    THSTAR, SVAVG, SWAVG, UAVG,
     &                    SVATZI, SWATZI, UATZI,
     &                    PTATZI, UATHE, SVATHE, SWATHE,
     &                    UAVH3, SVAVH3, SWAVH3, SWRMAX

C**   BVF    = Brunt-Vaisala frequency
C**   BVPRIM = 0.7*BVF
C**   SFCHF  = Surface heat flux (W/sq m)
C**   USTAR  = Surface friction velocity (m/s)
C**   WSTAR  = Convective scaling velocity (m/s)
C**   VPTGZI = Vertical potential temperature gradient from ZI to ZI+500
C**            (degrees/m)
C**   ZICONV = Hourly convective mixing height estimated by AERMET (m)
C**   ZIMECH = Hourly mechanical mixing height estimated by AERMET (m)
C**   OBULEN = Monin-Obukhov length (m)
C**   SFCZ0  = Surface roughness length (m)
C**   BOWEN  = Bowen ratio = sensible heat flux/latent heat flux
C**   ALBEDO = Albedo at the earth's surface (nondimensional)
C**   UREF   = Reference height wind speed (m/s)
C**   WDREF  = Reference height wind direction (degrees from north)
C**   UREFHT = Reference height for winds (m) (first nonmissing level
C**            of wind speed AND direction above 7.0*SFCZ0)
C**   TA     = Ambient temperature at a reference height (kelvin)
C**   TREFHT = Reference height for temperature (m) (first nonmissing
C**            level of temperature)
C**   ZI     = The mixing height used by AERMOD after any manipulation
C**            and massaging (m)
C**   NPLVLS = Number of levels in the observed hourly profile data
C**   NTGLVL = Number of levels of observed potential temperature gradient
C**   IFLAG  = Top of profile flag: 1 = top level, 0 = level below top
C**   PFLHT  = Profile height above local ground level (m)
C**   PFLWD  = Profile wind direction (degrees from north)
C**   PFLWS  = Profile wind speed (m/s)
C**   PFLTA  = Profile ambient temperature (kelvins)
C**   PFLSA  = Profile sigma_A (degrees)
C**   PFLSW  = Profile sigma_W (m/s)
C**   PFLSV  = Profile sigma_V (m/s), computed from sigma_A and wind speed
C**   PFLTG  = Profile of Vertical Potential Temperature Gradient (kelvin/m)
C**   PFLTGZ = Profile of VPTG heights (midpoint of interval) (m)
C**   GRIDHT = Gridded height (m)
C**   GRIDWD = Gridded wind direction (degrees from north)
C**   GRIDWS = Gridded wind speed (m/s)
C**   GRIDSW = Gridded sigma_W (m/s)
C**   GRIDSV = Gridded sigma_V (m/s)
C**   GRIDTG = Gridded vertical potential temperature gradient (deg/m)
C**   GRIDPT = Gridded potential temperature profile
C**   GRIDRHO= Gridded density profile
C**   GRIDEPS= Gridded tubulence dissipation rate (epsilon) profile for PVMRM/PVMRM2
C**   TG4PFL = Potential temperature gradient at 2.0 meters
C**   TG4XTR = Potential temperature gradient at 100.0 meters
C**   XLAT   = Station latitude, decimal degrees
C**   TSIGN  = Sign used for turning of wind: 1.0 for northern hemis.
C**                                          -1.0 for southern hemis.
C**   NDX4ZI = Index of gridded height immediately below ZI
C**   SVAVG  = Average sigma_V from the surface to ZI (m/s)
C**   SWAVG  = Average sigma_W from the surface to ZI (m/s)
C**   UAVG   = Average wind speed from the surface to ZI (m/s)
C**   SVATZI = sigma_V at ZI (m/s)
C**   SWATZI = sigma_W at ZI (m/s)
C**   UATZI  = Wind speed at ZI (m)
C**   PTATZI = Potential temperature at ZI (kelvin)
C**   SVATHE = Average sigma_V from the surface to HS for HS > ZI (m/s)
C**   SWATHE = Average sigma_W from the surface to HS for HS > ZI (m/s)
C**   UATHE  = Average wind speed from the surface to HS for HS > ZI (m/s)
C**   SVAVH3 = Average sigma_V from the surface to HE3 for penetrated plume
C**   SWAVH3 = Average sigma_W from the surface to HE3 for penetrated plume
C**   UAVH3  = Average wind speed from the surface to HE3 for penetrated plume
C**   SWRMAX = Residual vertical turbulence, average of measured sigma-w
C**            above ZI or 0.02 * UATZI


C***********************************************************************
C     This is The Global Variable Definition Block for Calculation
C***********************************************************************

      LOGICAL CALCS, WAKE
      LOGICAL SURFAC

      DOUBLE PRECISION PHID1, PHID2, PHIN1, PHIN2

      INTEGER :: IREC,   ISRC,   IGRP,   IAVE,   ITYP,  ISET,
     &           NUMREC, NUMSRC, NUMGRP, NUMAVE, NUMARC, NUMTYP,
     &           NUMCAP, NUMHOR, 
C ---            Include counters for all source types
     &           NUMPNT, NUMVOL, NUMAREA, NUMLINE, NUMPIT,
     &           NUMFLAT, IBKGRD(6), IO3SET(6),
     &           ICYEAR, NURBSRC, NUMURB, NPD, IFVSEC,
     &           IUCAT, IOLM, NUMOLM, IPSD, NUMPSD, IURB
      DOUBLE PRECISION :: XS, YS, ZS, QS, HS, DS, VS, TS, SYINIT,
     &                    SZINIT, XINIT, YINIT, ANGLE, XCNTR, YCNTR,
     &                    DSFACT, DSBH, DSBW,
c --- PRIME Modification -------------------------------------------
     &                    DSBL, XADJ, YADJ, B_SUBS, B_SUBL, RSCALE,
c ------------------------------------------------------------------
     &                    D, VD, E, WDRAD, WDSIN, WDCOS, ZBASE

      DOUBLE PRECISION, ALLOCATABLE :: PDIAM(:), PHI(:), PDENS(:),
     &                                 VGRAV(:), TSTOP(:), SCHMIDT(:),
     &                                 VDEP(:), WQCOR(:), DQCOR(:),
     &                                 PSCVRT(:), WASHOUT(:), ECOLL(:),
     &                                 AWDSIN(:), AWDCOS(:),
     &                                 AAFV(:)
      DOUBLE PRECISION :: WQCORG, GSCVRT, DQCORG, WASHOUTG, VSETL
      DOUBLE PRECISION :: XR, YR, X, Y, ZELEV, ZFLAG, ZR, ZEFF, DISTR,
     &                    ZHILL, HCRIT, ZRT, XDIST
      DOUBLE PRECISION :: HE, HSP, HEFLAT, HTER, HEMWAK, HEDHH, ZB, ZM,
     &                    HED1, HED2, HEN1, HEN2, HE3, HPEN, HED1M,
     &                    HED2M, HEN1M, HEN2M, HE3M, HSBL, QSUBN, QSUB3,
     &                    XY, XZ, SBID, FM, FB, DTDZ, DHF, DHFAER, DHP,
     &                    DHP1,DHP2, DHP3, DELT, DHPB, DHPM, XF, XMAX, 
     &                    XFM, XFB, XRAD, WPB, DHCRIT, HTEFF, CENTER,
     &                    Z4GAMMA, XTR4GAMMA
      DOUBLE PRECISION :: HESETL, HE3SETL, HV
      DOUBLE PRECISION :: US, SVS, SWS, TGS, TYS, PTS, UP, WDIR, DA,
     &                    ZLY, ZLB, RINIT, CB, CM, QTK, PPF, PSUBS, FHC,
     &                    SY, SYB, SYN, SY3, SZ, SZUPR, SYAMB, SZAMB,
     &                    SZAS, SZAD1, SZAD2, SZAN1, SZAN2, SYAN, SZA3,
     &                    SZB, SZBD, SZBN, SZ3, SZD1, SZD2, SZN1, SZN2,
     &                    SZEFF, SZSURF, SYA3, SYB3, SZB3, VSY3, VSIGY,
     &                    VSIGZ, VSYN, VSZD1,VSZD2, VSZN1, VSZN2, VSZ3,
     &                    SZD1M, SZD2M, SZN1M, SZN2M, SZ3M, U3, SV3, 
     &                    SW3, TGP
      DOUBLE PRECISION :: FSUBY, FSUBYD, FSUBYN, FSUBY3
      DOUBLE PRECISION :: FSUBZ, FSUBZD, FSUBZN, FSUBZ3, 
     &                    PHEE, FOPT, CWRAP, CLIFT, XMDBG,
     &                    CWRAPC, CLIFTC, FSUBYC, FSBY3C
      DOUBLE PRECISION :: UEFF, SVEFF, SWEFF, TGEFF,
     &                    UEFFD, SVEFFD, SWEFFD,
     &                    UEFFN, SVEFFN, SWEFFN,
     &                    UEFF3, SVEFF3, SWEFF3, TGEFF3,
     &                    EPSEFF, EPSEFFD, EPSEFF3,
     &                    XMIXED, XFINAL, ZMIDMX,
     &                    SIGVEFF, SIGVEFFD, SIGVEFF3,
     &                    SYEFF
      DOUBLE PRECISION :: SKEW, R, ALPHPD, BETAPD, ASUB1, ASUB2,
     &                    BSUB1, BSUB2, LAMDA1, LAMDA2
      DOUBLE PRECISION :: CHIW, CHIDW, CHINW, CHI3W,
     &                    CHIL, CHIDL, CHINL, CHI3L
      DOUBLE PRECISION :: GAMFACT

C**   AZSAVG = Average stack base elevation (m)
C**   ZBASE  = Base elevation used for potential temperature profile (m MSL)
C**   US     = Wind speed at stack height (m/s)
C**   UP     = Stack top wind speed for plume rise computations
C**   WDIR   = Stack top wind direction used for plume transport
C**   SVS    = sigma_V at stack height (m/s)
C**   SWS    = sigma_W at stack height (m/s)
C**   TGS    = Potential temperature gradient at stack top
C**   PTS    = Stack top potential temperature for plume rise
C**   xxEFF  = "effective" value for parameter xx
C**   FHC    = Function of Plume material above HCRIT
c**   PHEE   = "PHI" Term : Fraction of Plume Below Hcrit
C**   FSUBY  = Fy Term (Horizontal Gaussian term)
C**   FSUBYN = Fy Term (Horizontal Gaussian term) for the
C**            Indirect Source
C**   FSUBY3 = Fy Term (Horizontal Gaussian term) for the
C**            Penetrated Source
C**   HEDx   = Effective Source Heights for Direct Plume,
C**            x corresponding to each of the 2 distributions
C**   HENx   = Effective Source Heights for Indirect Plume,
C**            x corresponding to each of the 2 distributions
C**   HE3    = Effective Source Height for Penetrated Plume
C**   HEDxM  = Effective Source Heights for Direct Plume at Xm,
C**            x corresponding to each of the 2 distributions
C**   HENxM  = Effective Source Heights for Indirect Plume at Xm,
C**            x corresponding to each of the 2 distributions
C**   HE3M   = Effective Source Height for Penetrated Plume at Xm

C**   QSUBN  = Source Term for Indirect Source
C**   QSUB3  = Source Term for Penetrated Source
C**   SKEW   = Skewness of the Vertical Velocity
C**   R      = Lagrangian Correlation Coefficient
C**   ALPHPD = ALPHA Coefficient for the CBL PDF
C**   BETAPD = BETA  Coefficient for the CBL PDF
C**   ASUB1  =
C**   ASUB2  =
C**   BSUB1  =
C**   BSUB2  =
C**   LAMDA1 = Relative Frequencies of Updrafts
C**   LAMDA2 = Relative Frequencies of Downdrafts

C**   QSUM        = Sum of emissions for merged plumes in PVMRM/PVMRM2 option
C**   SUM_NO2RAT  = Sum of NO2/NOx ratios for merged plumes in PVMRM/PVMRM2

      DOUBLE PRECISION :: QSUM,  SUM_NO2RAT
      DOUBLE PRECISION :: QSUM3, SUM3_NO2RAT
      
      DOUBLE PRECISION :: CWMAX, CWMIN, CWMAX3, CWMIN3, 
     &                    DWMAX, DWMIN, DWMAX3, DWMIN3
      DOUBLE PRECISION :: HMNH,  HMXH,  HMNT,   HMXT, 
     &                    HMNH3, HMXH3, HMNT3,  HMXT3

C***********************************************************************
C     This is The Global Variable Definition Block for EVent Pathway
C***********************************************************************

      CHARACTER EVNAME*10, EVGRP*8
      INTEGER, ALLOCATABLE ::  EVAPER(:), EVDATE(:), EVJDAY(:),
     &                         IDXEV(:)

      ALLOCATABLE ::  EVNAME(:), EVGRP(:)




C***********************************************************************
C     This is The Global Variable Definition Block for OUtput Pathway
C***********************************************************************

      LOGICAL OUTPART, SUMMFILE, L_NoHeader(8), EVALFIL, TOXXFIL

      LOGICAL, ALLOCATABLE :: ANPART(:), ALLPARTS(:), ALLPARTG(:)

      CHARACTER (LEN=ILEN_FLD) :: THRFIL, PSTFIL, PLTFIL, ANNPST,
     &                            ANNPLT, THRFRM, PSTFRM, PLTFRM,
     &                            TOXFIL, SEAHRS, RNKFIL, RNKFRM,
     &                            EVLFIL, SUMFIL, MXDFRM,
     &                            MAXDLY, MAXDLY_BYYR, MAXDCONT_FILE

C --- Variable for specifying format for file outputs (default = 'FIX')
      CHARACTER (LEN = 3) :: FILE_FORMAT
      
      INTEGER, ALLOCATABLE :: NHIAVE(:,:), MAXAVE(:), IMXVAL(:),
     &                        IDYTAB(:), MAXFLE(:,:), IPSTFL(:,:), 
     &                        IPLTFL(:,:,:), IANPST(:), IANPLT(:), 
     &                        INHI(:), ITOXFL(:), ISEAHR(:), 
     &                        IMXDLY(:), IMXDLY_BYYR(:), MAXDCONT(:),
     &                        IRNKFL(:), IRKVAL(:), IANPART(:)
     
      DOUBLE PRECISION, ALLOCATABLE :: THRESH(:,:), TOXTHR(:),
     &                                 MAXD_THRESH(:)
     
      DOUBLE PRECISION, ALLOCATABLE :: AXR_SAV(:), AYR_SAV(:), 
     &                                 AZELEV_SAV(:), AZFLAG_SAV(:), 
     &                                 AZHILL_SAV(:)

      INTEGER, ALLOCATABLE :: IMXUNT(:,:), IPSUNT(:,:), IPSFRM(:,:), 
     &                        IPLUNT(:,:,:), IAPUNT(:), 
     &                        IANFRM(:), IPPUNT(:), ITXUNT(:), 
     &                        IRKUNT(:), IELUNT(:), IUPART(:),
     &                        ISHUNT(:), IMDUNT(:), IMDUNT_BYYR(:),
     &                        IMXDCUNT(:), MXD_RANK(:,:)

      ALLOCATABLE ::  THRFIL(:,:), PSTFIL(:,:), PLTFIL(:,:,:), 
     &                ANNPST(:), ANNPLT(:), TOXFIL(:), SEAHRS(:),
     &                RNKFIL(:), EVLFIL(:), MAXDLY(:),
     &                MAXDLY_BYYR(:), MAXDCONT_FILE(:)

      INTEGER, ALLOCATABLE :: IDCONC(:,:)

      INTEGER :: ITAB, NXTOX, NYTOX, NHOURS, IPAIR

      DOUBLE PRECISION, ALLOCATABLE :: TXCONC(:,:)



C***********************************************************************
C     This is The Global Variable Definition Block for Working Space
C***********************************************************************

      CHARACTER WORKID*12, DUMMY*12

      INTEGER :: IMIT, INUM, IDUM, INDAVE, INDGRP, INDVAL,
     &           ISC, IOERRN, NCPP, NRPP, NGPP, NPPX, NPPY
      REAL :: FNUM
      DOUBLE PRECISION :: DNUM

      ALLOCATABLE ::          WORKID(:)
      INTEGER, ALLOCATABLE :: IWRK2(:,:)

C     Declare Temporary Work Arrays for ZELEV and ZFLAG Receptor Data
      DOUBLE PRECISION, ALLOCATABLE :: ZETMP1(:), ZETMP2(:)
      DOUBLE PRECISION, ALLOCATABLE :: ZFTMP1(:), ZFTMP2(:)
      DOUBLE PRECISION, ALLOCATABLE :: ZHTMP1(:), ZHTMP2(:)


C***********************************************************************
C     AERMOD Model Data - Array Names, Array Limits, Named Common Blocks
C                         Necessary for Model Results
C     MODIFIED - 4/17/95   Output CONC/DEPOS in same model run
C***********************************************************************


C***********************************************************************
C     This is The Global Variable Definition Block For The Maximum
C     Value, Highest Value, Average Value, Annual Average Value and
C     Model Result Arrays.  Also Included are Calm/Missing Flag Arrays.
C***********************************************************************


      CHARACTER HCLMSG, MCLMSG, HMCLM

      DOUBLE PRECISION, ALLOCATABLE ::  HRVAL(:), AVEVAL(:,:,:,:),
     &                                  AERVAL(:), PRMVAL(:)
      DOUBLE PRECISION, ALLOCATABLE ::  HIVALU(:,:,:,:,:),
     &                                  HMAX(:,:,:,:)
      INTEGER, ALLOCATABLE ::  HMLOC(:,:,:,:),
     &                         HMDATE(:,:,:,:),
     &                         NHIDAT(:,:,:,:,:),
     &                         NHIDATMXD(:,:,:),
     &                         NHIDATMXD_BYYR(:,:,:,:)

      DOUBLE PRECISION, ALLOCATABLE ::  ANNVAL(:,:,:), AMXVAL(:,:,:),
     &                                  SHVALS(:,:,:,:,:), MXDVAL(:,:), 
     &                                  HIMXDLY(:,:,:),
     &                                  HIMXDLY_BYYR(:,:,:,:)
      INTEGER, ALLOCATABLE ::  IMXLOC(:,:,:), IMXDHR(:,:)
      INTEGER              ::  IANHRS, IANCLM, IANMSG,
     &                         NSEAHR(4,24), NSEACM(4,24)
      DOUBLE PRECISION, ALLOCATABLE ::  RMXVAL(:,:,:,:)
      INTEGER, ALLOCATABLE ::  MXDATE(:,:,:,:),
     &                         MXLOCA(:,:,:,:)
      INTEGER, ALLOCATABLE ::  NUMHRS(:), NUMCLM(:), NUMMSG(:)
      ALLOCATABLE ::           HCLMSG(:,:,:,:,:),
     &                         MCLMSG(:,:,:,:),
     &                         HMCLM(:,:,:,:)

      DOUBLE PRECISION, ALLOCATABLE ::  SUMANN(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE ::  SUMHNH(:,:,:), MXPMVAL(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE ::  SUMVAL_MAXD(:,:,:,:)

      INTEGER, ALLOCATABLE ::  MXPMLOC(:,:,:)

      DOUBLE PRECISION, ALLOCATABLE ::  CHI(:,:,:), HECNTR(:,:), 
     &                                  HECNTR3(:,:), PPFACT(:,:),
     &                                  UEFFS(:,:), UEFF3S(:,:),
     &                                  EPSEF(:,:), EPSEF3(:,:),
     &                                  FOPTS(:,:), 
     &                                  ABVAL(:,:), BCVAL(:,:)

      DOUBLE PRECISION, ALLOCATABLE ::  ARCMAX(:), QMAX(:), DXMAX(:),
     &                                  UMAX(:),
     &                        SVMAX(:), SWMAX(:), SYMAX(:), SY3MX(:),
     &                        U3MAX(:), HEMAX(:), ARCCL(:), SZMAX(:),
     &                        CHIDMW(:), CHINMW(:), CHI3MW(:),
     &                        CHIDML(:), CHINML(:), CHI3ML(:),
     &                        HSBLMX(:)

      LOGICAL, ALLOCATABLE :: CHIMASK(:,:,:)

C***********************************************************************
C     This is The Global Variable Definition Block For The
C     EVENT Model Result Arrays
C***********************************************************************

      DOUBLE PRECISION, ALLOCATABLE ::  EV_AVEVAL(:), HRVALS(:,:), 
     &                                  GRPVAL(:,:), BACKHR(:,:),
     &                                  GRPAVE(:), BACKAVE(:),
     &                                  BACKANN(:), BACKSEASHR(:,:,:)

C --- Declare allocatable array for Original Event Concentrations;
C     to be compared to GRPAVE calculated value for QA purpose
      DOUBLE PRECISION, ALLOCATABLE :: EV_OrigConc(:)

      INTEGER ::  EV_NUMHRS, EV_NUMCLM, EV_NUMMSG, ISTAHR, IENDHR


C***********************************************************************
C    
C     BLOCK DATA area for initializing global data
C
C***********************************************************************
C***********************************************************************
C     Initialize Model Version Number, VERSN (Year, Julian Day), as a
C     Character Variable
C***********************************************************************

C---- VERSN is now a 6-character variable to accomodate leading qualifier
C     character, such as 'B' for Beta version or 'D' for Draft version.
      DATA VERSN /' 15181'/   ! June 30, 2015

C     Initialize C_METVER to blanks in case there is an error opening 
C     the surface met file or with reading the version date, otherwise
C     C_METVER will be undefined when writing the page headers.
      DATA C_METVER /'      '/

C --- Initialize array of hourly O3 values to 40 ppb (78.4 ug/m3)
      DATA O3_Max24hr /144*78.4D0/   ! 144 = 24hr/day*6sectors

C***********************************************************************
C     Input/Output File Units and Input/Output File Names
C***********************************************************************

      DATA INUNIT/ 7/, IOUNIT/ 8/, PVMDBG/ 9/, OLMDBG/ 9/, ARMDBG/ 9/,
     &     IERUNT/10/, IERWRT/11/, IDPUNT/12/, IDPUN2/14/, IRSUNT/15/, 
     &     IHREMI/16/, IEVUNT/17/, ITEVUT/18/, MFUNIT/19/, INCUNT/20/,  ! IBGUNT is now an array for BGSECTOR
     &     MPUNIT/21/, ISUNIT/22/, IPUNIT/23/, DBGUNT/24/, DBMUNT/25/, 
     &     AREADBUNT/26/, PRMDBUNT/27/, ISUMUNT/28/, GDEPDBG/29/, 
     &     PDEPDBG/30/,                                                 ! IO3UNT is now an array for O3SECTOR
     &     RDISPUNT/ 3/                                                 ! RDISPUNT is for RELDISP debug file for PVMRM
C*#
      DATA INPFIL/' '/, OUTFIL/' '/


C***********************************************************************
C     Initialize Keyword Array
C***********************************************************************

      INTEGER, PRIVATE :: I

      DATA (KEYWD(I),I=1,IKN) /
     &   'STARTING','FINISHED','TITLEONE','TITLETWO','MODELOPT',
     &   'AVERTIME','POLLUTID','HALFLIFE','DCAYCOEF','DEBUGOPT',
     &   'ELEVUNIT','FLAGPOLE','RUNORNOT','EVENTFIL','SAVEFILE',
     &   'INITFILE','MULTYEAR','ERRORFIL','GASDEPDF','GDSEASON',
     &   'GDLANUSE','GASDEPVD','URBANOPT','NO2EQUIL','NO2STACK',
     &   'OZONEVAL','OZONEFIL','O3VALUES','OZONUNIT','O3SECTOR',
     7   'LOW_WIND','ARMRATIO',
     &   'LOCATION','SRCPARAM','BUILDHGT','BUILDWID','BUILDLEN',
     &   'XBADJ   ','YBADJ   ','EMISFACT','EMISUNIT','PARTDIAM',
     &   'MASSFRAX','PARTDENS','METHOD_2','CONCUNIT','DEPOUNIT',
     &   'GASDEPOS','HOUREMIS','NO2RATIO','AREAVERT','URBANSRC',
     &   'SRCGROUP','OLMGROUP','PSDGROUP','BACKGRND','BACKUNIT',
     &   'BGSECTOR','BLPINPUT',
     &   'INCLUDED','EVENTPER','EVENTLOC','GRIDCART','GRIDPOLR',
     &   'DISCCART','DISCPOLR','EVALCART','SURFFILE','PROFFILE',
     &   'PROFBASE','SURFDATA','UAIRDATA','SITEDATA','STARTEND',
     &   'DAYRANGE','SCIMBYHR','WDROTATE','WINDCATS','NUMYEARS',
     &   'RECTABLE','MAXTABLE','DAYTABLE','SUMMFILE','MAXIFILE',
     &   'POSTFILE','PLOTFILE','TOXXFILE','SEASONHR','EVENTOUT',
     &   'RANKFILE','EVALFILE','FILEFORM','MAXDAILY','MXDYBYYR',
     &   'MAXDCONT','NOHEADER'/


C***********************************************************************
C     Initialize Miscellaneous Variables
C***********************************************************************

      DATA IPROC /366*1/, IPROCL/366*1/, EXPLIM /-50.0D0/
      DATA UCAT /1.54D0, 3.09D0, 5.14D0, 8.23D0, 10.8D0/
      DATA MODOPS /24*'         '/

C***********************************************************************
C     Initialize distance factors used in determining when to switch
C     to point source approximation for area sources under the FASTAREA
C     option (formerly the TOXICS option).
C***********************************************************************

C     STAB. CLASS         A      B      C       D       E       F
C                        ***    ***    ***     ***     ***     ***
      DATA VIRTPNT_URB /3.5D0, 3.5D0, 5.5D0, 10.5D0, 15.5D0, 15.5D0/,
     &     VIRTPNT_RUR /3.5D0, 5.5D0, 7.5D0, 12.5D0, 15.5D0, 25.5D0/


C***********************************************************************
C     Initialize Setup Status Arrays
C***********************************************************************

      DATA ICSTAT/50*0/, ISSTAT/50*0/, IRSTAT/50*0/, IMSTAT/50*0/,
     &     IOSTAT/50*0/, IESTAT/50*0/

C***********************************************************************
C     Cross-reference table of Keywords by I?STAT array index:
C     
C     Array    Index  Keyword            Array    Index  Keyword
C     -----    -----  -------            -----    -----  -------    
C     ICSTAT:    1 = starting            ISSTAT:    1 = starting
C                2 = titleone                       2 = location
C                3 = titletwo                       3 = srcparam
C                4 = modelopt                       4 = buildhgt
C                5 = avertime                       5 = buildwid
C                6 = pollutid                       7 = emisfact
C                7 = halflife                       8 = emisunit
C                8 = dcaycoef                       9 = partdiam
C               11 = flagpole                      10 = massfrax
C               12 = runornot                      11 = partdens
C               13 = eventfil                      15 = elevunit
C               14 = savefile                      16 = houremis
C               15 = initfile                      17 = concunit
C               16 = multyear                      18 = depounit
C               17 = errorfil                      19 = areavert
C               18 = gdseason                      20 = included
C               19 = gasdepdf                      21 = buildlen
C               20 = gdlanuse                      22 = xbadj   
C               21 = gasdepvd                      23 = ybadj   
C               22 = debugopt                      24 = srcgroup
C               23 = urbanopt                      26 = gasdepos
C               24 = ozoneval                      27 = method_2
C               25 = O3VALUES                      28 = urbansrc
C               26 = ozonefil                      29 = no2ratio
C               27 = ozonunit                      30 = olmgroup
C               28 = no2stack                      34 = psdgroup
C               29 = no2equil                      40 = backgrnd
C               30 = low_wind                      41 = backunit
C               31 = o3sector                      42 = bgsector
C               50 = finished                      43 = blavgval
C                                                  50 = finished

C     IRSTAT:    1 = starting            IESTAT:    1 = starting
C                2 = gridcart            (EVENT)    2 = eventper
C                3 = gridpolr                       3 = eventloc
C                4 = disccart                      10 = included
C                5 = discpolr                      50 = finished
C                8 = evalcart    
C                9 = elevunit    
C               11 = included    
C               50 = finished    
C    
C     IMSTAT:    1 = starting
C                2 = surffile
C                3 = proffile
C                4 = surfdata
C                5 = uairdata
C                6 = startend
C                7 = dayrange
C                8 = wdrotate
C                9 = sitedata
C               10 = profbase
C               11 = windcats
C               12 = scimbyhr
C               13 = numyears
C               50 = finished
C    
C     IOSTAT:    1 = starting            IOSTAT:    1 = starting
C (non-EVENT)    2 = rectable            (EVENT)    2 = eventout
C                3 = maxtable                      13 = fileform
C                4 = daytable                      25 = finished
C                5 = maxifile
C                6 = postfile
C                7 = plotfile
C                8 = toxxfile
C                9 = seasonhr
C               10 = rankfile
C               11 = evalfile
C               12 = summfile
C               13 = fileform
C               14 = maxdaily
C               15 = mxdybyyr
C               16 = maxdcont
C               18 = noheader
C               50 = finished
C                                            
C***********************************************************************


C***********************************************************************
C     Initialize Gridded Profile Height Array
C***********************************************************************

      DATA GRIDHT / 
     &   0.0D0, 0.5D0, 1.0D0, 2.0D0, 4.0D0,  8.0D0, 14.0D0, 20.0D0,
     &  30.0D0,  40.0D0,  50.0D0,  60.0D0,  70.0D0,  80.0D0,  90.0D0,
     &  100.0D0, 120.0D0, 140.0D0, 160.0D0, 180.0D0, 200.0D0, 250.0D0,
     &  300.0D0, 350.0D0, 400.0D0, 450.0D0, 500.0D0, 550.0D0, 600.0D0,
     &  650.0D0, 700.0D0, 750.0D0, 800.0D0, 850.0D0, 900.0D0, 950.0D0,
     & 1000.0D0, 1050.0D0, 1100.0D0, 1150.0D0, 1200.0D0, 1250.0D0,
     & 1300.0D0, 1350.0D0, 1400.0D0, 1450.0D0, 1500.0D0, 1550.0D0,
     & 1600.0D0, 1650.0D0, 1700.0D0, 1750.0D0, 1800.0D0, 1850.0D0,
     & 1900.0D0, 1950.0D0, 2000.0D0, 2100.0D0, 2200.0D0, 2300.0D0,
     & 2400.0D0, 2500.0D0, 2600.0D0, 2700.0D0, 2800.0D0, 2900.0D0,
     & 3000.0D0, 3100.0D0, 3200.0D0, 3300.0D0, 3400.0D0, 3500.0D0,
     & 3600.0D0, 3700.0D0, 3800.0D0, 3900.0D0, 4000.0D0, 4100.0D0,
     & 4200.0D0, 4300.0D0, 4400.0D0, 4500.0D0, 4600.0D0, 4700.0D0,
     & 4800.0D0, 4900.0D0, 5000.0D0/


C***********************************************************************
C     Initialize Error Code and Message Arrays
C***********************************************************************

      DATA ERRCOD(  1)/'100'/,
     & ERRMSG(  1)/'Invalid Pathway Specified. The Troubled Pathway is'/
      DATA ERRCOD(  2)/'105'/,
     & ERRMSG(  2)/'Invalid Keyword Specified. The Troubled Keyword is'/
      DATA ERRCOD(  3)/'109'/,
     & ERRMSG(  3)/'Too many fields specified on runstream image; MAX='/
      DATA ERRCOD(  4)/'110'/,
     & ERRMSG(  4)/'Keyword is Not Valid for This Pathway.  Keyword is'/
      DATA ERRCOD(  5)/'111'/,
     & ERRMSG(  5)/'LOW_WIND Keyword specified without LOWWIND Options'/
      DATA ERRCOD(  6)/'112'/,
     & ERRMSG(  6)/'User-specified minimum Sigma-V on LOW_WIND Keyword'/
      DATA ERRCOD(  7)/'113'/,
     & ERRMSG(  7)/'User-specified minimum WindSpeed on LOW_WIND Keywd'/
      DATA ERRCOD(  8)/'114'/,
     & ERRMSG(  8)/'User-specified maximum FRAN on the LOW_WIND Keywrd'/
      DATA ERRCOD(  9)/'115'/,
     & ERRMSG(  9)/'STARTING or FINISHED Out of Sequence:  Pathway =  '/
      DATA ERRCOD( 10)/'116'/,
     & ERRMSG( 10)/'Vector Wind Speeds specified on MODELOPT Keyword  '/
      DATA ERRCOD( 11)/'120'/,
     & ERRMSG( 11)/'Pathway is Out of Sequence:  Pathway =            '/
      DATA ERRCOD( 12)/'121'/,
     & ERRMSG( 12)/'LowWind1 Beta Option specified on MODELOPT Keyword'/
      DATA ERRCOD( 13)/'122'/,
     & ERRMSG( 13)/'LowWind2 Beta Option specified on MODELOPT Keyword'/
      DATA ERRCOD( 14)/'125'/,
     & ERRMSG( 14)/'Missing FINISHED-Runstream File Incomplete: ISTAT='/
      DATA ERRCOD( 15)/'130'/,
     & ERRMSG( 15)/'Missing Mandatory Keyword.  The Missing Keyword is'/
      DATA ERRCOD( 16)/'131'/,
     & ERRMSG( 16)/'Minimum sigmav value (SVmin) for LowWind1 Beta Opt'/
      DATA ERRCOD( 17)/'132'/,
     & ERRMSG( 17)/'Minimum sigmav value (SVmin) for LW2/LW3 Beta Opt '/
      DATA ERRCOD( 18)/'133'/,
     & ERRMSG( 18)/'Maximum FRAN value (FRANmax) for LW2/LW3 Beta Opt '/
      DATA ERRCOD( 19)/'135'/,
     & ERRMSG( 19)/'Nonrepeatable Keyword or Recursed INCLUDED: Keywrd'/
      DATA ERRCOD( 20)/'136'/,
     & ERRMSG( 20)/'Option Conflict - Multiple LowWind Opts specified:'/
      DATA ERRCOD( 21)/'140'/,
     & ERRMSG( 21)/'Invalid Order of Keyword.  The Troubled Keyword is'/
      DATA ERRCOD( 22)/'141'/,
     & ERRMSG( 22)/'Conflicting Options for NO2 conversion specified: '/
      DATA ERRCOD( 23)/'142'/,
     & ERRMSG( 23)/'Following Keyword Invalid Without PVMRM or OLM:   '/
      DATA ERRCOD( 24)/'143'/,
     & ERRMSG( 24)/'Following Keyword Invalid Without PVMRM Option:   '/
      DATA ERRCOD( 25)/'144'/,
     & ERRMSG( 25)/'Following Keyword Invalid Without OLM Option:     '/
      DATA ERRCOD( 26)/'145'/,
     & ERRMSG( 26)/'Following Keyword Invalid Without ARM or ARM2:    '/
      DATA ERRCOD( 27)/'146'/,
     & ERRMSG( 27)/'PSDGROUP Keyword Specified without PSDCREDIT Opt. '/
      DATA ERRCOD( 28)/'147'/,
     & ERRMSG( 28)/'Following Option is Invalid with PSDCREDIT Option:'/
      DATA ERRCOD( 29)/'148'/,
     & ERRMSG( 29)/'Both OZONEVAL and O3VALUES keywords are specified '/
      DATA ERRCOD( 30)/'149'/,
     & ERRMSG( 30)/'Conflicting options specified on MODELOPT keyword:'/
      DATA ERRCOD( 31)/'150'/,
     & ERRMSG( 31)/'Conflicting Options: MULTYEAR Option with         '/
      DATA ERRCOD( 32)/'151'/,
     & ERRMSG( 32)/'Non-DFAULT NoUrbTran option selected on MODELOPT  '/
      DATA ERRCOD( 33)/'152'/,
     & ERRMSG( 33)/'ELEVUNIT card must be first for this Pathway:     '/
      DATA ERRCOD( 34)/'153'/,
     & ERRMSG( 34)/'Conflicting Opts: MAXDCONT with Re-Start or MULTYR'/
      DATA ERRCOD( 35)/'154'/,
     & ERRMSG( 35)/'Conflicting options:  SCIM cannot be used with    '/
      DATA ERRCOD( 36)/'155'/,
     & ERRMSG( 36)/'Conflicting Decay Keyword. Inputs Ignored for     '/
      DATA ERRCOD( 37)/'156'/,
     & ERRMSG( 37)/'Option ignored - not valid with SCIM.  Option =   '/
      DATA ERRCOD( 38)/'157'/,
     & ERRMSG( 38)/'Wet SCIM Not Supported - Wet SCIM Inputs Ignored  '/
      DATA ERRCOD( 39)/'158'/,
     & ERRMSG( 39)/'EMISUNIT Keyword Used With More Than 1 Output Type'/
      DATA ERRCOD( 40)/'159'/,
     & ERRMSG( 40)/'EMISUNIT Keyword Used With the Following Keyword: '/
      DATA ERRCOD( 41)/'160'/,
     & ERRMSG( 41)/'Duplicate ORIG Secondary Keyword for GRIDPOLR:    '/
      DATA ERRCOD( 42)/'161'/,
     & ERRMSG( 42)/'MAXDCONT option already defined for source group: '/
      DATA ERRCOD( 43)/'162'/,
     & ERRMSG( 43)/'Option only applies to 1-hr NO2 or 1-hr SO2 NAAQS:'/
      DATA ERRCOD( 44)/'163'/,
     & ERRMSG( 44)/'Option only applies to 24h PM25, 1h NO2 or 1h SO2:'/
      DATA ERRCOD( 45)/'164'/,
     & ERRMSG( 45)/'NOHEADER selected for non-specified output option:'/
      DATA ERRCOD( 46)/'165'/,
     & ERRMSG( 46)/'Inconsistent temporally-varying BACKGRND options: '/
      DATA ERRCOD( 47)/'166'/,
     & ERRMSG( 47)/'BGSECTOR/O3SECTOR option invalid w/o BG/O3 Inputs:'/
      DATA ERRCOD( 48)/'167'/,
     & ERRMSG( 48)/'Inconsistent temporally-varying O3VALUES options: '/
      DATA ERRCOD( 49)/'168'/,
     & ERRMSG( 49)/'Hourly BACKGRND already specified for this sector:'/
      DATA ERRCOD( 50)/'170'/,
     & ERRMSG( 50)/'Invalid Secondary Keyword for Receptor Grid:      '/
      DATA ERRCOD( 51)/'171'/,
     & ERRMSG( 51)/'Sector ID specified without Sector-varying Option:'/
      DATA ERRCOD( 52)/'175'/,
     & ERRMSG( 52)/'Missing Secondary Keyword END for Receptor Grid:  '/
      DATA ERRCOD( 53)/'180'/,
     & ERRMSG( 53)/'Conflicting Secondary Keyword for Receptor Grid:  '/
      DATA ERRCOD( 54)/'181'/,
     & ERRMSG( 54)/'BULKRN Delta-T & SolarRad option for SBL was used '/
      DATA ERRCOD( 55)/'182'/,
     & ERRMSG( 55)/'MMIF-generated meteorological inputs were used    '/
      DATA ERRCOD( 56)/'183'/,
     & ERRMSG( 56)/'Non-DFAULT option for MMIF-generated data without '/
      DATA ERRCOD( 57)/'184'/,
     & ERRMSG( 57)/'PROFFILE heights > 999m; inputs could be from MMIF'/
      DATA ERRCOD( 58)/'185'/,
     & ERRMSG( 58)/'Either No Sources or No Receptors are specified!!!'/
      DATA ERRCOD( 59)/'186'/,
     & ERRMSG( 59)/'THRESH_1MIN 1-min ASOS wind speed threshold used  '/
      DATA ERRCOD( 60)/'187'/,
     & ERRMSG( 60)/'ADJ_U* Beta Option for Low Winds used in AERMET   '/
      DATA ERRCOD( 61)/'188'/,
     & ERRMSG( 61)/'Non-Default ADJ_U* Option used in AERMET without  '/
      DATA ERRCOD( 62)/'189'/,
     & ERRMSG( 62)/'No Keywords for OU Path and No PERIOD/ANNUAL Aves.'/
      DATA ERRCOD( 63)/'190'/,
     & ERRMSG( 63)/'Incompatible Option Used With SAVEFILE or INITFILE'/
      DATA ERRCOD( 64)/'191'/,
     & ERRMSG( 64)/'PM25, 1h NO2 or SO2 w/o MAXIFILE incompatible with'/
      DATA ERRCOD( 65)/'192'/,
     & ERRMSG( 65)/'FASTALL option also implies use of FASTAREA option'/
      DATA ERRCOD( 66)/'193'/,
     & ERRMSG( 66)/'Units keyword specified without appropriate option'/
      DATA ERRCOD( 67)/'194'/,
     & ERRMSG( 67)/'DEBUGOPT input option is invalid or not applicable'/
      DATA ERRCOD( 68)/'195'/,
     & ERRMSG( 68)/'Incompatible Keyword used with GASDEPVD option    '/
      DATA ERRCOD( 69)/'196'/,
     & ERRMSG( 69)/'Gas deposition algorithms are non-DFAULT options  '/
      DATA ERRCOD( 70)/'197'/,
     & ERRMSG( 70)/'METHOD_2 for particulates is a non-DFAULT option  '/
      DATA ERRCOD( 71)/'198'/,
     & ERRMSG( 71)/'TOXICS Option obsolete; see Users Guide Addendum  '/
      DATA ERRCOD( 72)/'199'/,
     & ERRMSG( 72)/'Non-DFAULT BETA Option Required for               '/

      DATA ERRCOD( 73)/'200'/,
     & ERRMSG( 73)/'Missing Parameter(s). No Options Specified For    '/
      DATA ERRCOD( 74)/'201'/,
     & ERRMSG( 74)/'Not Enough Parameters Specified For the Keyword of'/
      DATA ERRCOD( 75)/'202'/,
     & ERRMSG( 75)/'Too Many Parameters Specified For the Keyword of  '/
      DATA ERRCOD( 76)/'203'/,
     & ERRMSG( 76)/'Invalid Parameter Specified.  Troubled Parameter: '/
      DATA ERRCOD( 77)/'204'/,
     & ERRMSG( 77)/'Regulatory DFAULT Conflicts with Non-DFAULT Option'/
      DATA ERRCOD( 78)/'205'/,
     & ERRMSG( 78)/'No Option Parameter Setting.  Forced by Default to'/
      DATA ERRCOD( 79)/'206'/,
     & ERRMSG( 79)/'Regulatory DFAULT Overrides Non-DFAULT Option For '/
      DATA ERRCOD( 80)/'207'/,
     & ERRMSG( 80)/'No Parameters Specified. Default Values Will Used.'/
      DATA ERRCOD( 81)/'208'/,
     & ERRMSG( 81)/'Illegal Numerical Field Encountered in            '/
      DATA ERRCOD( 82)/'209'/,
     & ERRMSG( 82)/'Negative Value Appears For Non-negative Variable. '/
      DATA ERRCOD( 83)/'210'/,
     & ERRMSG( 83)/'Num Ranked values on RANKFILE > MAXTABLE value for'/

      DATA ERRCOD( 84)/'211'/,
     & ERRMSG( 84)/'Duplicate Averaging Period Specified for Keyword  '/
      DATA ERRCOD( 85)/'212'/,
     & ERRMSG( 85)/'END Encountered Without (X,Y) Points Properly Set '/
      DATA ERRCOD( 86)/'213'/,
     & ERRMSG( 86)/'ELEV Input Inconsistent With Option: Input Ignored'/
      DATA ERRCOD( 87)/'214'/,
     & ERRMSG( 87)/'ELEV Input Inconsistent With Option: Defaults Used'/
      DATA ERRCOD( 88)/'215'/,
     & ERRMSG( 88)/'FLAG Input Inconsistent With Option: Input Ignored'/
      DATA ERRCOD( 89)/'216'/,
     & ERRMSG( 89)/'FLAG Input Inconsistent With Option: Defaults Used'/
      DATA ERRCOD( 90)/'217'/,
     & ERRMSG( 90)/'More Than One Delimiter In A Field for Keyword    '/
      DATA ERRCOD( 91)/'218'/,
     & ERRMSG( 91)/'Number of (X,Y) Points Does Not Match Number of   '/
      DATA ERRCOD( 92)/'219'/,
     & ERRMSG( 92)/'Urban ID field is too long (>8); first 12 char:   '/
      DATA ERRCOD( 93)/'220'/,
     & ERRMSG( 93)/'Missing Origin (Use Default = 0,0) In GRIDPOLR    '/
      DATA ERRCOD( 94)/'221'/,
     & ERRMSG( 94)/'Missing Dist or Direction Setting In Polar Network'/
      DATA ERRCOD( 95)/'222'/,
     & ERRMSG( 95)/'O3SECTOR or BGSECTOR Value is out of order:       '/
      DATA ERRCOD( 96)/'223'/,
     & ERRMSG( 96)/'Missing Distance or Degree Field in               '/
C --- New messages '224' '225' and '226' added for undefined SrcID on 
C     the SRCGROUP, OLMGROUP, or PSDGROUP keywords
      DATA ERRCOD( 97)/'224'/,
     & ERRMSG( 97)/'SrcID specified on SRCGROUP keyword not defined:  '/
      DATA ERRCOD( 98)/'225'/,
     & ERRMSG( 98)/'SrcID specified on OLMGROUP keyword not defined:  '/
      DATA ERRCOD( 99)/'226'/,
     & ERRMSG( 99)/'SrcID specified on PSDGROUP keyword not defined:  '/
      DATA ERRCOD(100)/'227'/,
     & ERRMSG(100)/'O3SECTOR or BGSECTOR Width is out of range:       '/
      DATA ERRCOD(101)/'228'/,
     & ERRMSG(101)/'Default(s) Used for Missing Parameters on Keyword '/
      DATA ERRCOD(102)/'229'/,
     & ERRMSG(102)/'Too Many Parameters - Inputs Ignored on Keyword   '/
      DATA ERRCOD(103)/'230'/,
     & ERRMSG(103)/'Source ID field is too long (>12); first 12 chars:'/
      DATA ERRCOD(104)/'231'/,
     & ERRMSG(104)/'Too Many Numerical Values Specified for           '/
      DATA ERRCOD(105)/'232'/,
     & ERRMSG(105)/'OLMGroup ID field is too long (>8); first 12 char:'/
      DATA ERRCOD(106)/'233'/,
     & ERRMSG(106)/'Building Dimensions Specified for Non-POINT Source'/
      DATA ERRCOD(107)/'234'/,
     & ERRMSG(107)/'Too Many Sectors Input for                        '/
      DATA ERRCOD(108)/'235'/,
     & ERRMSG(108)/'Num of SRCGRPs exceeds limit for EVT name; Set=999'/
      DATA ERRCOD(109)/'236'/,
     & ERRMSG(109)/'Not Enough BUILDHGTs Specified for SourceID       '/
      DATA ERRCOD(110)/'237'/,
     & ERRMSG(110)/'Not Enough BUILDWIDs Specified for SourceID       '/
      DATA ERRCOD(111)/'238'/,
     & ERRMSG(111)/'Not Enough BACKGRND Concentration Values Specified'/
      DATA ERRCOD(112)/'239'/,
     & ERRMSG(112)/'Not Enough QFACTs Specified for SourceID          '/
      DATA ERRCOD(113)/'240'/,
     & ERRMSG(113)/'Inconsistent Number of Particle Categories for    '/
      DATA ERRCOD(114)/'241'/,
     & ERRMSG(114)/'Not Enough BUILDLENs Specified for SourceID       '/
      DATA ERRCOD(115)/'242'/,
     & ERRMSG(115)/'No Particle Cat. or Gas Depos. Specified for SRCID'/
      DATA ERRCOD(116)/'243'/,
     & ERRMSG(116)/'Wet depos (DEPOS, WDEP, WETDPLT) incompatible with'/
      DATA ERRCOD(117)/'244'/,
     & ERRMSG(117)/'Source parameters are missing or incomplete for   '/
      DATA ERRCOD(118)/'245'/,
     & ERRMSG(118)/'SrcGroup ID field is too long (>8); first 12 char:'/
      DATA ERRCOD(119)/'246'/,
     & ERRMSG(119)/'Not Enough XBADJs Specified for SourceID          '/
      DATA ERRCOD(120)/'247'/,
     & ERRMSG(120)/'Not Enough YBADJs Specified for SourceID          '/
      DATA ERRCOD(121)/'248'/,
     & ERRMSG(121)/'Either BGVALs or BGFILE missing for this sector:  '/
      DATA ERRCOD(122)/'249'/,
     & ERRMSG(122)/'Source elevation is missing (-9999.0); SRCID =    '/
      DATA ERRCOD(123)/'250'/,
     & ERRMSG(123)/'Duplicate XPNT/DIST or YPNT/DIR Specified for GRID'/
      DATA ERRCOD(124)/'252'/,
     & ERRMSG(124)/'Duplicate Receptor Network ID Specified.  NETID = '/
      DATA ERRCOD(125)/'253'/,
     & ERRMSG(125)/'PSDGROUP ID field is too long (>8); first 12 char:'/
      DATA ERRCOD(126)/'256'/,
     & ERRMSG(126)/'EVALFILE Option Used Without EVALCART Receptors   '/
      DATA ERRCOD(127)/'259'/,
     & ERRMSG(127)/'Receptor elevation is missing (-9999.0); IREC =   '/
      DATA ERRCOD(128)/'260'/,
     & ERRMSG(128)/'Number of EMISFACT/O3VALUES/BACKGRND values > max:'/
      DATA ERRCOD(129)/'261'/,
     & ERRMSG(129)/'Not Enough O3VALUES Ozone Concentrations Specified'/
      DATA ERRCOD(130)/'262'/,
     & ERRMSG(130)/'First Vertex Does Not Match LOCATION for AREAPOLY '/
      DATA ERRCOD(131)/'264'/,
     & ERRMSG(131)/'Too Many Vertices Specified for AREAPOLY Source   '/
      DATA ERRCOD(132)/'265'/,
     & ERRMSG(132)/'Not Enough Vertices Specified for AREAPOLY Source '/
      DATA ERRCOD(133)/'266'/,
     & ERRMSG(133)/'Invalid shape defined (area=0) for AREAPOLY source'/
      DATA ERRCOD(134)/'   '/,
     & ERRMSG(134)/'                                                  '/
      DATA ERRCOD(135)/'271'/,
     & ERRMSG(135)/'O3FILE w/o O3VALs; full conv for hrs with miss O3 '/
      DATA ERRCOD(136)/'272'/,
     & ERRMSG(136)/'Upper bound rank > Lower bound rank for MAXDCONT: '/
      DATA ERRCOD(137)/'273'/,
     & ERRMSG(137)/'Range of ranks for MAXDCONT THRESH Opt is limited:'/
      DATA ERRCOD(138)/'   '/,
     & ERRMSG(138)/'                                                  '/
      DATA ERRCOD(139)/'279'/,
     & ERRMSG(139)/'Multiple URBANOPT/URBANSRC inputs not allowed for:'/
      DATA ERRCOD(140)/'280'/,
     & ERRMSG(140)/'Number of Output Types Specified Exceeds Max:NTYP='/
      DATA ERRCOD(141)/'   '/,
     & ERRMSG(141)/'                                                  '/
      DATA ERRCOD(142)/'282'/,
     & ERRMSG(142)/'Following SRCID Included in Multiple OLMGROUPs:   '/
      DATA ERRCOD(143)/'283'/,
     & ERRMSG(143)/'OZONEVAL, O3VALUES or OZONEFIL Keyword Needed for '/
      DATA ERRCOD(144)/'284'/,
     & ERRMSG(144)/'Invalid POLLUTID Specified for PVMRM/OLM; Must Use'/
      DATA ERRCOD(145)/'285'/,
     & ERRMSG(145)/'BACKGROUND and BACKGRND are invalid as Source IDs '/
      DATA ERRCOD(146)/'286'/,
     & ERRMSG(146)/'Following SRCID Included in Multiple PSDGROUPs:   '/
      DATA ERRCOD(147)/'287'/,
     & ERRMSG(147)/'PSDGROUP ID Must be INCRCONS, RETRBASE or NONRBASE'/
      DATA ERRCOD(148)/'288'/,
     & ERRMSG(148)/'Use of "*" for repeated values not meaningful for '/
      DATA ERRCOD(149)/'289'/,
     & ERRMSG(149)/'Source defined as both particulate and gaseous    '/
      DATA ERRCOD(150)/'290'/,
     & ERRMSG(150)/'This array limit exceeded; possible coding error: '/
      DATA ERRCOD(151)/'291'/,
     & ERRMSG(151)/'Filename specified is too long. Maximum length =  '/
      DATA ERRCOD(152)/'292'/,
     & ERRMSG(152)/'Potential problem with Fortran format specifier:  '/
      DATA ERRCOD(153)/'293'/,
     & ERRMSG(153)/'User-specified met data format not used;  use FREE'/
      DATA ERRCOD(154)/'294'/,
     & ERRMSG(154)/'PERIOD and ANNUAL averages are both selected for  '/
      DATA ERRCOD(155)/'295'/,
     & ERRMSG(155)/'Invalid Averaging Period Specified for SCREEN Mode'/
      DATA ERRCOD(156)/'296'/,
     & ERRMSG(156)/'Averaging Period .NE. 1-Hr for TOXXFILE Option    '/
      DATA ERRCOD(157)/'297'/,
     & ERRMSG(157)/'Aver. Period must be .LE. 24 for EVENT Processing '/
      DATA ERRCOD(158)/'298'/,
     & ERRMSG(158)/'Results reported for source group ALL include     '/
      DATA ERRCOD(159)/'299'/,
     & ERRMSG(159)/'SRCGROUP ALL is missing, but is NOT required for  '/

      DATA ERRCOD(160)/'300'/,
     & ERRMSG(160)/'Specified SRCID Has Not Been Defined Yet: KEYWORD='/
      DATA ERRCOD(161)/'301'/,
     & ERRMSG(161)/'Urban Area ID Has Not Been Defined.  URBID =      '/
      DATA ERRCOD(162)/'302'/,
     & ERRMSG(162)/'Following SRCID Included in Multiple Urban Areas: '/
      DATA ERRCOD(163)/'303'/,
     & ERRMSG(163)/'Urban ID has already been defined.  URBID =       '/
      DATA ERRCOD(164)/'305'/,
     & ERRMSG(164)/'Stack height > or = EPA formula height for SRCID: '/
      DATA ERRCOD(165)/'310'/,
     & ERRMSG(165)/'Attempt to Define Duplicate LOCATION Card for SRC:'/
      DATA ERRCOD(166)/'313'/,
     & ERRMSG(166)/'Attempt to Define Duplicate EVENTPER card for     '/
      DATA ERRCOD(167)/'315'/,
     & ERRMSG(167)/'Attempt to Define Duplicate SRCPARAM Card for SRC:'/
      DATA ERRCOD(168)/'317'/,
     & ERRMSG(168)/'Specified SRCID is not included in any PSDGROUP:  '/
      DATA ERRCOD(169)/'318'/,
     & ERRMSG(169)/'No Sources Defined for Urban Area.  URBID =       '/
      DATA ERRCOD(170)/'319'/,
     & ERRMSG(170)/'No Sources Included in Specified Source Group:    '/
      DATA ERRCOD(171)/'320'/,
     & ERRMSG(171)/'Input Parameter May Be Out-of-Range for Parameter '/
      DATA ERRCOD(172)/'321'/,
     & ERRMSG(172)/'BACKGROUND concs are NOT included in any SRCGROUP!'/
      DATA ERRCOD(173)/'322'/,
     & ERRMSG(173)/'Release Height Exceeds Effective Depth for OPENPIT'/
      DATA ERRCOD(174)/'323'/,
     & ERRMSG(174)/'BACKGRND included w/o BACKGRND keyword for SrcGrp:'/
      DATA ERRCOD(175)/'324'/,
     & ERRMSG(175)/'Release Height Exceeds 3000 Meters for SRCID:     '/
      DATA ERRCOD(176)/'325'/,
     & ERRMSG(176)/'Negative Exit Velocity (Set=1.0E-5) for SRCID:    '/
      DATA ERRCOD(177)/'330'/,
     & ERRMSG(177)/'Mass Fraction Parameters Do Not Sum to 1. for Src '/
      DATA ERRCOD(178)/'332'/,
     & ERRMSG(178)/'Mass Fraction Parameter Out-of-Range for Source   '/
      DATA ERRCOD(179)/'334'/,
     & ERRMSG(179)/'Particle Density Out-of-Range for Source          '/
      DATA ERRCOD(180)/'335'/,
     & ERRMSG(180)/'Particle Diameter Out-of-Range for Source         '/
      DATA ERRCOD(181)/'336'/,
     & ERRMSG(181)/'NO2RATIO Missing or Invalid for OLM/PVMRM - SrcID:'/
      DATA ERRCOD(182)/'338'/,
     & ERRMSG(182)/'Neg Emis Rate Cannot be Used with OLM/PVMRM. Src: '/
      DATA ERRCOD(183)/'340'/,
     & ERRMSG(183)/'Possible Error in PROFBASE Input:  Value is < 0   '/
      DATA ERRCOD(184)/'341'/,
     & ERRMSG(184)/'Emissions in HOUREMIS file < -90; set to 0.0 for  '/
      DATA ERRCOD(185)/'342'/,
     & ERRMSG(185)/'Src ID Mismatch in Hourly Emissions File for ID = '/
      DATA ERRCOD(186)/'344'/,
     & ERRMSG(186)/'Missing HOUREMIS fields; EmisRate set = 0. KURDAT='/
      DATA ERRCOD(187)/'345'/,
     & ERRMSG(187)/'Problem processing the HOUREMIS file.   KURDAT =  '/
      DATA ERRCOD(188)/'346'/,
     & ERRMSG(188)/'Too many fields for HOUREMIS file.     KURDAT =   '/
      DATA ERRCOD(189)/'350'/,
     & ERRMSG(189)/'Julian Day Out Of Range at                        '/
      DATA ERRCOD(190)/'352'/,
     & ERRMSG(190)/'The "H6H" field is no longer required for MULTYEAR'/
      DATA ERRCOD(191)/'353'/,
     & ERRMSG(191)/'Urban Roughness Length (m) May Be Out-of-Range:   '/
      DATA ERRCOD(192)/'360'/,
     & ERRMSG(192)/'2-Digit Year Specified: Valid for Range 1950-2049 '/
      DATA ERRCOD(193)/'361'/,
     & ERRMSG(193)/'Multiyear PERIOD/ANNUAL values for NO2/SO2 require'/
      DATA ERRCOD(194)/'362'/,
     & ERRMSG(194)/'Multiyear 1h NO2/SO2 processing not applicable for'/
      DATA ERRCOD(195)/'363'/,
     & ERRMSG(195)/'Multiyr 24h/Ann PM25 processing not applicable for'/
      DATA ERRCOD(196)/'365'/,
     & ERRMSG(196)/'Year Input is Greater Than 2147                   '/
      DATA ERRCOD(197)/'370'/,
     & ERRMSG(197)/'Invalid Date: 2/29 In a Non-leap Year.            '/
      DATA ERRCOD(198)/'380'/,
     & ERRMSG(198)/'This Input Variable is Out-of-Range:              '/
      DATA ERRCOD(199)/'381'/,
     & ERRMSG(199)/'Latitude in Surface File Is Not Valid:            '/
      DATA ERRCOD(200)/'382'/,
     & ERRMSG(200)/'Error Decoding Latitude:                          '/
      DATA ERRCOD(201)/'384'/,
     & ERRMSG(201)/'Not enough fields specified for HOUREMIS; KURDAT ='/
      DATA ERRCOD(202)/'386'/,
     & ERRMSG(202)/'PARTDIAM and METHOD_2 specified for same SRCID:   '/
      DATA ERRCOD(203)/'387'/,
     & ERRMSG(203)/'METHOD_2 option already specified for this SRCID: '/
      DATA ERRCOD(204)/'390'/,
     & ERRMSG(204)/'Aspect ratio (L/W) of LINE source greater than 100'/
      DATA ERRCOD(205)/'391'/,
     & ERRMSG(205)/'Aspect ratio (L/W) of AREA source greater than 100'/
      DATA ERRCOD(206)/'392'/,
     & ERRMSG(206)/'Aspect ratio (L/W) of OPENPIT is greater than 10  '/
      DATA ERRCOD(207)/'394'/,
     & ERRMSG(207)/'Met data may be from outdated version of AERMET:  '/
      DATA ERRCOD(208)/'395'/,
     & ERRMSG(208)/'Met. Data Error; Incompatible Version of AERMET:  '/
      DATA ERRCOD(209)/'396'/,
     & ERRMSG(209)/'AERMET Version Out-dated or Non-standard; Version:'/
      DATA ERRCOD(210)/'397'/,
     & ERRMSG(210)/'SCREEN option used without use of SCREEN Met Data '/
      DATA ERRCOD(211)/'398'/,
     & ERRMSG(211)/'SCREEN met used without specifying SCREEN option  '/
      DATA ERRCOD(212)/'399'/,
     & ERRMSG(212)/'EXP format specified with no applicable file types'/

      DATA ERRCOD(213)/'400'/,
     & ERRMSG(213)/'Output values exceed format limit; use OU FILEFORM'/
      DATA ERRCOD(214)/'405'/,
     & ERRMSG(214)/'Value of PHEE Exceeds 1.0 on KURDAT =             '/
      DATA ERRCOD(215)/'406'/,
     & ERRMSG(215)/'Number of Vertices Exceeds Max (NVMAX) for SRCID: '/
      DATA ERRCOD(216)/'409'/,
     & ERRMSG(216)/'Error Allocating Storage for Setup/Result Arrays! '/
      DATA ERRCOD(217)/'410'/,
     & ERRMSG(217)/'Wind Direction Out-of-Range.  KURDAT =            '/
      DATA ERRCOD(218)/'413'/,
     & ERRMSG(218)/'Number of Threshold Events > 999999 for Ave Period'/
      DATA ERRCOD(219)/'415'/,
     & ERRMSG(219)/'MAXDCONT THRESH not reached within range of ranks '/
      DATA ERRCOD(220)/'420'/,
     & ERRMSG(220)/'Wind Speed Out-of-Range.   KURDAT =               '/
      DATA ERRCOD(221)/'430'/,
     & ERRMSG(221)/'Ambient Temperature Data Out-of-Range.  KURDAT =  '/
      DATA ERRCOD(222)/'432'/,
     & ERRMSG(222)/'Friction Velocity Out-of-Range.   KURDAT =        '/
      DATA ERRCOD(223)/'435'/,
     & ERRMSG(223)/'Surface Roughness Length Out-of-Range.  KURDAT =  '/
      DATA ERRCOD(224)/'438'/,
     & ERRMSG(224)/'Convective Velocity Data Out-of-Range.  KURDAT =  '/
      DATA ERRCOD(225)/'439'/,
     & ERRMSG(225)/'Monin-Obukhov Length Out-of-Range.  KURDAT =      '/
      DATA ERRCOD(226)/'440'/,
     & ERRMSG(226)/'Calm Hour Identified in Meteorology Data File at  '/
      DATA ERRCOD(227)/'441'/,
     & ERRMSG(227)/'Vert Pot Temp Grad abv ZI set to min .005, KURDAT='/
      DATA ERRCOD(228)/'442'/,
     & ERRMSG(228)/'Vert Pot Temp Grad abv ZI exceeds 0.1 K/m, KURDAT='/
      DATA ERRCOD(229)/'450'/,
     & ERRMSG(229)/'Record Out of Sequence in Meteorological File at: '/
      DATA ERRCOD(230)/'452'/,
     & ERRMSG(230)/'Missing hourly BACKGRND w/o BGSUB, KURDAT/Sector ='/
      DATA ERRCOD(231)/'453'/,
     & ERRMSG(231)/'BGSUB for missing hourly BACKGRND, KURDAT/Sector ='/
      DATA ERRCOD(232)/'454'/,
     & ERRMSG(232)/'Date/time Mismatch: BACKGRND File, KURDAT/Sector ='/
      DATA ERRCOD(233)/'455'/,
     & ERRMSG(233)/'Date/time Mismatch: Hourly Emission File, KURDAT ='/
      DATA ERRCOD(234)/'456'/,
     & ERRMSG(234)/'Date/time Mismatch on Surface & Profile.  KURDAT ='/
      DATA ERRCOD(235)/'457'/,
     & ERRMSG(235)/'Date/time Mismatch: OZONEFIL File, KURDAT/Sector ='/
      DATA ERRCOD(236)/'458'/,
     & ERRMSG(236)/'O3SUB for missing hourly O3 value, KURDAT/Sector ='/
      DATA ERRCOD(237)/'459'/,
     & ERRMSG(237)/'No Hrly O3 & No Sub; Use Full Conversion, KURDAT ='/
      DATA ERRCOD(238)/'460'/,
     & ERRMSG(238)/'Missing Hour Identified in Meteor. Data File at   '/
      DATA ERRCOD(239)/'465'/,
     & ERRMSG(239)/'Number of Profile Levels Exceeds Max:   MXPLVL =  '/
      DATA ERRCOD(240)/'470'/,
     & ERRMSG(240)/'Mixing Height Value is < or = 0.0.   KURDAT =     '/
      DATA ERRCOD(241)/'474'/,
     & ERRMSG(241)/'WS RefHt invalid (<0.001); Not msg or clm: KURDAT='/
      DATA ERRCOD(242)/'475'/,
     & ERRMSG(242)/'WS reference height is higher than 100m.  KURDAT ='/
      DATA ERRCOD(243)/'480'/,
     & ERRMSG(243)/'Less than 1yr for MULTYEAR, MAXDCONT or ANNUAL Ave'/
      DATA ERRCOD(244)/'481'/,
     & ERRMSG(244)/'Data Remaining After End of Year. Number of Hours='/
      DATA ERRCOD(245)/'482'/,
     & ERRMSG(245)/'Too many years modeled for 24h-PM25 1h-NO2 1h-SO2:'/
      DATA ERRCOD(246)/'483'/,
     & ERRMSG(246)/'User Start Date is Earlier Than Start of Met File '/
      DATA ERRCOD(247)/'484'/,
     & ERRMSG(247)/'Restart Date < STARTEND date or start of Met File '/
      DATA ERRCOD(248)/'485'/,
     & ERRMSG(248)/'MULTYR DataGap; Restart Date < STARTEND or MetFile'/
      DATA ERRCOD(249)/'486'/,
     & ERRMSG(249)/'MULTYR Date Overlap; STARTEND Date < Restart Date '/
      DATA ERRCOD(250)/'487'/,
     & ERRMSG(250)/'MULTYR Date Overlap; MetFile Start < Restart Date '/
      DATA ERRCOD(251)/'488'/,
     & ERRMSG(251)/'First met HR.ne.1; ST results may not be valid    '/
      DATA ERRCOD(252)/'489'/,
     & ERRMSG(252)/'First met HR.ne.1; EV results may not be valid for'/
      DATA ERRCOD(253)/'490'/,
     & ERRMSG(253)/'Problem reading SURFFILE date for EVENTS; MNDYHR ='/
      DATA ERRCOD(254)/'491'/,
     & ERRMSG(254)/'MAXDCONT option requires 1st Hr of met data = 01; '/
      DATA ERRCOD(255)/'492'/,
     & ERRMSG(255)/'SURFDATA YR .NE. 1st YR of file, adj to match file'/
      DATA ERRCOD(256)/'493'/,
     & ERRMSG(256)/'SURFDATA YR must match 1st YR of file for DAYRANGE'/
      DATA ERRCOD(257)/'495'/,
     & ERRMSG(257)/'Surface met file does not include enough variables'/
      DATA ERRCOD(258)/'496'/,
     & ERRMSG(258)/'Total precipitation in SURFFILE is zero (0.0) with'/
      DATA ERRCOD(259)/'497'/,
     & ERRMSG(259)/'Possible code ERROR!!! EVENT mismatch for EVENTID:'/
      DATA ERRCOD(260)/'498'/,
     & ERRMSG(260)/'Possible code ERROR! MAXDCONT mismatch GRP/RNK/REC'/
      DATA ERRCOD(261)/'499'/,
     & ERRMSG(261)/'PRIME plume rise error; check stack parameters for'/

      DATA ERRCOD(262)/'500'/,
     & ERRMSG(262)/'Fatal Error Occurs Opening the Data File of       '/
      DATA ERRCOD(263)/'501'/,
     & ERRMSG(263)/'Dup Filename! Fatal Error Opening the Data File of'/
      DATA ERRCOD(264)/'510'/,
     & ERRMSG(264)/'Fatal Error Occurs During Reading of the File of  '/
      DATA ERRCOD(265)/'520'/,
     & ERRMSG(265)/'Fatal Error Occurs During Writing to the File of  '/
      DATA ERRCOD(266)/'530'/,
     & ERRMSG(266)/'CAUTION! Met Station ID Mismatch with SURFFILE for'/
      DATA ERRCOD(267)/'531'/,
     & ERRMSG(267)/'CAUTION! Met Station ID Missing from SURFFILE for '/
      DATA ERRCOD(268)/'540'/,
     & ERRMSG(268)/'No RECTABLE/MAXTABLE/DAYTABLE for Average Period  '/
      DATA ERRCOD(269)/'550'/,
     & ERRMSG(269)/'File Unit/Name Conflict for the Output Option:    '/
      DATA ERRCOD(270)/'555'/,
     & ERRMSG(270)/'File Unit/Name conflict across options: GRP# AVE  '/
      DATA ERRCOD(271)/'560'/,
     & ERRMSG(271)/'User Specified File Unit .LE. 30 for OU Keyword:  '/
      DATA ERRCOD(272)/'565'/,
     & ERRMSG(272)/'Possible Conflict With Dynamically Allocated FUNIT'/
      DATA ERRCOD(273)/'570'/,
     & ERRMSG(273)/'Problem Reading Temporary Event File for Event:   '/
      DATA ERRCOD(274)/'580'/,
     & ERRMSG(274)/'End of File Reached Trying to Read the File of    '/
      DATA ERRCOD(275)/'585'/,
     & ERRMSG(275)/'Output data file for INITFILE option was not found'/
      DATA ERRCOD(276)/'590'/,
     & ERRMSG(276)/'The INITFILE filename matches a SAVEFILE filename '/
      DATA ERRCOD(277)/'592'/,
     & ERRMSG(277)/'MAXIFILE includes data past start of MULTYEAR run '/
      DATA ERRCOD(278)/'593'/,
     & ERRMSG(278)/'POSTFILE includes data past start of MULTYEAR run '/
C --- Included new messages regarding special processing for 1hr NO2/SO2
C     and 24hr PM25
      DATA ERRCOD(279)/'276'/,
     & ERRMSG(279)/'Special proc for 1h-NO2/SO2 24hPM25 NAAQS disabled'/
      DATA ERRCOD(280)/'277'/,
     & ERRMSG(280)/'Specified option not applicable for this pollutant'/
C --- Included new message regarding QSUM = 0.0 runtime error in PVMRM; 
C     this should never occur, but could indicate a programming error
      DATA ERRCOD(281)/'411'/,
     & ERRMSG(281)/'Possible ERROR in PVMRM_CALC! QSUM=0.0 @ Rec# Date'/
      DATA ERRCOD(282)/'412'/,
     & ERRMSG(282)/'Possible ERROR in PVMRM_CALC! QSUM=0.0 @ Evt# Date'/
C --- Additional messages 
      DATA ERRCOD(283)/'316'/,
     & ERRMSG(283)/'Specified SRCID is not included in any SRCGROUP:  '/
      DATA ERRCOD(284)/'314'/,
     & ERRMSG(284)/'Specified GRP index and SRC index is duplicated:  '/
      DATA ERRCOD(285)/'123'/,
     & ERRMSG(285)/'LowWind3 Beta Option specified on MODELOPT Keyword'/
      DATA ERRCOD(286)/'479'/,
     & ERRMSG(286)/'Potential temperature gradient is out-of-range:   '/
      DATA ERRCOD(287)/'494'/,
     & ERRMSG(287)/'Meander factor (FRAN) exceeds max @ YR MN DY ISRC:'/

C --- New messages for buoyant line processing
      DATA ERRCOD(288)/'388'/,
     & ERRMSG(288)/'Input buoyant line sources not in correct order:  '/
      DATA ERRCOD(289)/'389'/,
     & ERRMSG(289)/'Rotated buoyant line sources not in correct order:'/

      END MODULE MAIN1


      MODULE BUOYANT_LINE
C***********************************************************************
C     This is The Global Variable Definition Block for the New BUOYLINE
C     Source Algorithm - January 2015
C***********************************************************************

      TYPE BLINEDATA
         INTEGER   (kind=4)   :: ISRCNUM
         CHARACTER (len=12)   :: SRCID
         DOUBLE PRECISION     :: XBEG, YBEG, XEND, YEND                 ! untranslated, unrotated
         DOUBLE PRECISION     :: XBEG_TR1, YBEG_TR1, XEND_TR1, YEND_TR1 ! translated, rotate #1 w/TCOR
         DOUBLE PRECISION     :: ELEV, BLQS, BLHS
      END TYPE BLINEDATA
      
C     ISRCNUM = source number in the list of of sources for the model run
C     SRCID   = source ID in the control file for this line number
C     XBEG    = x-coordinate of beginning of line (as entered by user)
C     YBEG    = y-coordinate of beginning of line (as entered by user)
C     XEND    = x-coordinate of end of line (as entered by user)
C     YEND    = y-coordinate of end of line (as entered by user)
C     XBEG_TR1 = x-coordinate of beginning of translated, rotated line with TCOR
C     YBEG_TR1 = y-coordinate of beginning of translated, rotated line with TCOR
C     XEND_TR1 = x-coordinate of end of translated, rotated line with TCOR
C     YEND_TR1 = y-coordinate of end of translated, rotated line with TCOR
C     ELEV    = elevation of line
C     BLQS    = emission rate of line source in g/s
C     BLHS    = release height

      TYPE (BLINEDATA), ALLOCATABLE :: BLINEPARMS (:)

C     Whereas the source coordinates in BLINEPARMS are for the entire
C      line, the following are for the segments of the line.
      DOUBLE PRECISION, ALLOCATABLE :: XS_SCS(:,:), YS_SCS(:)
      DOUBLE PRECISION, ALLOCATABLE :: XS_RCS(:,:), YS_RCS(:,:)

C      _SCS = source coordinate system, initial translation/rotation and
C             performed only once
C      _RCS = rotated coordinate system, translated/rotated for wind dir
C
      LOGICAL          :: L_BLSOURCE
      INTEGER          :: NBLINES
      DOUBLE PRECISION :: BLAVGBLEN, BLAVGBHGT, BLAVGBWID, BLAVGLWID,
     &                    BLAVGBSEP, BLAVGFPRM, BLTA

C     NBLINES   = number of lines in the buoyant line source

C     The following are input via the source keyword BLAVGVAL:
C     BLAVGBLEN  = average building length (m)
C     BLAVGBHGT  = average building height (m)
C     BLAVGBWID  = average building width (m)
C     BLAVGBSEP  = average building separation (m)
C     BLAVGFPRM  = average buoyancy parameter (m^4/s^3)
C     BLTA       = ambient temperature read from AERMET's surface file
C                  (saved because TA is modified if there are point sources)

      DOUBLE PRECISION, dimension (7)  :: BL_XDIST = (0.0D0)
      DOUBLE PRECISION, dimension (7)  :: DH = (0.0D0)
      DOUBLE PRECISION, dimension (10) :: DEL
      DOUBLE PRECISION :: BL_XFB, LEFF, LD, R0, BL_XFINAL, BL_XFS
      DOUBLE PRECISION :: XOR, YOR, TCOR, SINTCOR, COSTCOR
      DOUBLE PRECISION :: FPRMXNB, XMATCH
      
C     LEFF     = effective building length
C     R0       = edge radius
C     LD       = LEFF * sin (theta), where theta = angle between flow
C                vector and the orientation of the line source
C     BL_XDIST = array of intermediate distances for downswsh calculations
C     FPRMXNB  = fprime * number of lines
C     BL_XFB   = distance to full buoyancy
C     BL_XFS   = distance to final rise for stable conditions
C     XMATCH   =
C     AVFACT   =
C     TCOR     = rotation angle to align the x-axis with the long side 
C                of the buoyant line source

C     For the translated/rotated receptors
C      to be allocated with the number of receptors
      DOUBLE PRECISION, ALLOCATABLE :: XR_SCS(:), YR_SCS(:)
      DOUBLE PRECISION, ALLOCATABLE :: XR_RCS(:), YR_RCS(:)
C
C      _SCS = source coordinate system, initial translation/rotation and
C             performed only once
C      _RCS = rotated coordinate system, translated/rotated for wind dir
C
C     The following arrays are used to save the coordinates from the 
C       first tranlastion/rotation of the receptors when MAXDCONT is used
      DOUBLE PRECISION, ALLOCATABLE :: XR_SCS_SAV(:), YR_SCS_SAV(:)


C     Arrays for total and partial concentration from lines, 
C      to be allocated with the number of receptors
      DOUBLE PRECISION, ALLOCATABLE :: CHIBL(:), PARTCH(:)

      END MODULE BUOYANT_LINE
      
C-----------------------------------------------------------------------
c --- The following MODULE subprograms replace the *.pri "INCLUDE"
c     files formerly used for global data storage for PRIME, and the
c     /PLU/-named COMMON block used in a few subroutines.
c --- R.W. Brode, U.S. EPA/OAQPS/AQMG, 10/19/2009


      MODULE PRIME_PARAMS

c --- Formerly part of PARAMS.PRI "INCLUDE" File:
c
c----------------------------------------------------------------------
c --- PARAMETER statements                                        PRIME
c----------------------------------------------------------------------

      INTEGER, PARAMETER :: io5=7,io6=8

c --- FORTRAN I/O unit numbers:
c           IO5 - Control file                  - input  - formatted
c           IO6 - List file                     - output - formatted
c

      END MODULE PRIME_PARAMS


      MODULE PRIME_NUMPARM

c --- Formerly NUMPARM.PRI "INCLUDE" File:
c
c----------------------------------------------------------------------
c --- COMMON BLOCK /NUMPARM/ -- Parameters used in the            PRIME
c                               numerical plume rise algorithm
c----------------------------------------------------------------------
c
      INTEGER, PARAMETER :: mxnw=5000
      INTEGER, PARAMETER :: mxent=10
      INTEGER, PARAMETER :: mxentp1=mxent+1
      INTEGER :: nstep, nent
      DOUBLE PRECISION :: gravi,rgas,zmin,ds,slast,rp,
     1           alphap(mxent),betap(mxent),xcat(mxentp1)

c
c --- GENERAL PARAMETER definitions:
c          MXNW - Maximum number of downwind distances for numerical
c                 plume rise integration (should be set equal to
c                 SLAST/DS)
c         MXENT - Maximum number of perturbed entrainment coefficients
c                 entered
c
c --- FORTRAN I/O unit numbers:
c           IO5 - Control file                  - input  - formatted
c           IO6 - List file                     - output - formatted
c
c --- NUMPARM Global Variables:
c
c         GRAVI - real    - Acceleration due to gravity (m/s**2)
c          RGAS - real    - Gas constant (m**2/s**2/deg. K)
c          ZMIN - real    - Minimum plume centerline height (m)
c            DS - real    - Step size (m) in the numerical plume
c                           rise algorithm
c         NSTEP - integer - Internal save frequency of plume rise
c                           calculations (i.e., every DS*NSTEP meters)
c                           (NOTE: this the frequency with which the
c                           results are saved internally -- not that
c                           passed back from the NUMRISE routine)
c         SLAST - real    - Termination distance (m) of the plume rise
c                           calculation
c            RP - real    - Radiation coefficient (kg/m**2/deg. K**3/s)
c   ALPHAP(mxent) - real array - Perturbed entrainment coefficients
c                                (parallel)
c    BETAP(mxent) - real array - Perturbed entrainment coefficients
c                                (normal)
c   XCAT(mxentp1) - real array - Downwind distances (m) for which each
c                                perturbed entrainment coefficient
c                                (ALPHAP, BETAP) is valid (NENT+1 values
c                                for NENT entrainment coefficients).
c            NENT - integer    - Number of perturbed entrainment
c                                coefficients entered

      END MODULE PRIME_NUMPARM


      MODULE PRIME_DFSN

c --- Formerly DFSN.PRI "INCLUDE" File:
c
c----------------------------------------------------------------------
c --- COMMON BLOCK /DFSN/ -- Parameters used in the            PRIME
c                            PRIME turbulence and diffusion
c                            subroutines
c----------------------------------------------------------------------
c
      DOUBLE PRECISION :: afac,xbyrmax,wiz0,wiy0,wfz,wfy,
     &        dua_ua,xdecay,xdecayi
c
c --- DFSN Global Variables:
c
c          AFAC - real    - Diffusion transitions to ambient (with
c                           virtual source) when wake turbulence decays
c                           to AFAC*(ambient turbulence intensity) for
c                           PG classes 4, 5, and 6
c       XBYRMAX - real    - Upper limit on distance from upwind face
c                           of bldg to transition point for ambient
c                           diffusion
c       WIZ,WIY - real    - Base Turbulence intensities in wake
c       WFZ,WFY - real    - Scaling factors for sigmaz and sigmay
c        DUA_UA - real    - [Ua-U]/Ua in wake at downwind face of bldg
c                                U: average speed in wake
c                               Ua: ambient speed
c         DECAY - real    - Exponent for turbulence intensity change
c                           with distance from downwind face of bldg
c        DECAYI - real    - 1/DECAY
c     RURLIZ(6) - real    - Rural turbulence intensities in z
c     RURLIY(6) - real    - Rural turbulence intensities in y
c     URBNIZ(6) - real    - Urban turbulence intensities in z
c     URBNIY(6) - real    - Urban turbulence intensities in y
c --- Ambient turbulence intensities are inferred from Briggs (1973)
c --- "Diffusion estimation for small emissions", ATDL-106;

      END MODULE PRIME_DFSN


      MODULE PRIME_WAKEDAT

c --- Formerly WAKEDAT.PRI "INCLUDE" File:
c
c----------------------------------------------------------------------
c --- COMMON BLOCK /WAKEDAT/ -- Parameters used in the            PRIME
c                               PRIME wake and streamline
c                               subroutines
c----------------------------------------------------------------------
c
      logical lrurl
      INTEGER, PARAMETER :: mxntr=50
      INTEGER :: nwak,ncav
      DOUBLE PRECISION :: Hb,Wb,xLb,Rb,HR,xLR,xLC,
     &                    xbadj,ybadj,Ub,Urh,
     &                    xwak(mxntr),szwak(mxntr),sywak(mxntr),
     &                    drwak(mxntr),
     &                    xcav(mxntr),szcav(mxntr),sycav(mxntr),
     &                    fqcav,
     &                    vsigy, vsigz, vsigyc, vsigzc, zint

      DOUBLE PRECISION :: third  ! constant = 1/3 used in various places
                                 ! initialized for PRIME in sub. WAKINI

c --- GENERAL PARAMETER definitions:
c         MXNTR - Maximum number of downwind distances for which
c                 numerical plume rise will be reported
c
c --- WAKEDAT Global Variables:
c
c            HB - real    - Building height (m)
c            WB - real    - Building width (crosswind) - (m)
c           XLB - real    - Building length (alongwind) - (m)
c            RB - real    - Scale length (m)
c            HR - real    - Maximum cavity height (m) above ground
c           XLR - real    - Length of downwind cavity (m) from
c                           downwind face of building
c           XLC - real    - Length of roof cavity (m)
c         XBADJ - real    - Distance along the wind from the stack to
c                           the origin of the building (upwind center
c                           of effective building)
c         YBADJ - real    - Distance crosswind from the stack to
c                           the origin of the building (upwind center
c                           of effective building)
c            Ub - real    - Wind speed (m/s) at the height of bldg
c           Urh - real    - Wind speed (m/s) at release height
c
c          NWAK - integer - Number of downwind distances at which
c                           wake properties are tabulated (LE mxntr)
c   XWAK(mxntr) - real    - Downwind distance (m) from source
c  SZWAK(mxntr) - real    - Sigma-z (m) at position XWAK
c  SYWAK(mxntr) - real    - Sigma-y (m) at position XWAK
c  DRWAK(mxntr) - real    - Plume growth rate at position XWAK expressed
c                           as d/dx(plume radius) for equivalent top-hat
c          NCAV - integer - Number of downwind distances at which
c                           wake properties of cavity source are
c                           tabulated (LE mxntr)
c   XCAV(mxntr) - real    - Downwind distance (m) from primary source
c  SZCAV(mxntr) - real    - Sigma-z (m) for cavity source
c  SYCAV(mxntr) - real    - Sigma-y (m) for cavity source
c         FQCAV - real    - Fraction of plume mass captured by cavity
c         ISTAB - integer - PG stability class
c         LRURL - logical - Rural dispersion when .TRUE.
c         VSIGZ - real    - Virtual source sigma (m) for sigma-z beyond wake
c         VSIGY - real    - Virtual source sigma (m) for sigma-y beyond wake
c        VSIGZC - real    - Virtual source sigma (m) for sigma-z beyond wake
c                           for cavity source
c        VSIGYC - real    - Virtual source sigma (m) for sigma-y beyond wake
c                           for cavity source

      END MODULE PRIME_WAKEDAT


      MODULE PRIME_AMBIENT

c --- Formerly AMBIENT.PRI "INCLUDE" File:
c
c----------------------------------------------------------------------
c --- COMMON BLOCK /AMBIENT/ -- Selected met. data at one         PRIME
c                               grid cell;  used in numerical
c                               plume rise computation
c----------------------------------------------------------------------
c
      INTEGER, PARAMETER :: mxnz=100
      INTEGER, PARAMETER :: mxnzp1=mxnz+1
      INTEGER :: NZA
      DOUBLE PRECISION :: uamb(mxnz),ramb(mxnz),dedz(mxnzp1),tamb(mxnz),
     1           zfacea(mxnzp1),zgpta(mxnz),tamb0,ramb0,adia,ptgrad0

c --- GENERAL PARAMETER definitions:
c          MXNZ - Maximum number of vertical layers in
c                 the meteorological data
c
c --- COMMON BLOCK /AMBIENT/ Variables:
c
c                    NZA - integer - Number of layers
c             UAMB(mxnz) - real    - Wind speed profile (m/s) - winds
c                                    defined at cell CENTERS
c             RAMB(mxnz) - real    - Ambient air density profile
c                                    (kg/m**3) - defined at cell CENTERS
c           DEDZ(mxnzp1) - real    - Pot. temperature gradient profile
c                                    (deg. K/m) - defined at cell FACES
c             TAMB(mxnz) - real    - Temperature profile (deg .K) -
c                                    defined at cell CENTERS
c         ZFACEA(mxnzp1) - real    - Heights of layer faces (m)
c            ZGPTA(mxnz) - real    - Heights of layer centers (m)
c                  TAMB0 - real    - Surface air temperature (deg. K)
c                  RAMB0 - real    - Surface air density (kg/m**3)
c                   ADIA - real    - Dry adiabatic lapse rate (deg. K/m)
c                PTGRAD0 - real    - Minimum potential temperature lapse
c                                    rate (deg. K/m)

      END MODULE PRIME_AMBIENT


      MODULE PRIME_PLU

c --- Formerly COMMON /PLU/ in selected PRIME subroutines:
c
c----------------------------------------------------------------------
c --- Notation --- in (KG,M,S) units
c               S:      LENGTH ALONG PLUME CENTERLINE
c               X:      PLUME LOCATION (downwind from source)
c               Y:      PLUME LOCATION (crosswind from source)
c               Z:      PLUME HEIGHT
c               R:      PLUME RADIUS
c               U:      PLUME HORIZONTAL (ALONGWIND) VELOCITY COMPONENT
c               V:      PLUME CROSSWIND VELOCITY COMPONENT
c               W:      PLUME VERTICAL VELOCITY COMPONENT
c               USC:    VELOCITY ALONG PLUME CENTERLINE
c               PHI:    ANGLE BETWEEN PLUME TRAJECTORY AND GROUND
c               DEN:    PLUME DENSITY
c               TP:     PLUME TEMPERATURE
c----------------------------------------------------------------------

      DOUBLE PRECISION :: S,X,Y,Z,R,U,V,W,USC,PHI,DEN,TP

      END MODULE PRIME_PLU
