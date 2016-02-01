      PROGRAM AERMOD
C=======================================================================
C            MAIN Module of the AMS/EPA Regulatory Model - AERMOD
C                           Version Dated 15181
C 
C                               June 30, 2015
C
C               *** SEE AERMOD MODEL CHANGE BULLETIN MCB#11 ***
C
C       ON THE SUPPORT CENTER FOR REGULATORY AIR MODELS (SCRAM) WEBSITE
C
C                      http://www.epa.gov/scram001/
C
C=======================================================================
C
C       This revised version of AERMOD (dated 15181) includes several
C       modifications relative to the previous version (dated 14134); 
C       see AERMOD Model Change Bullentin (MCB) #11 and AERMOD User's 
C       Guide Addendum.
C
C       http://www.epa.gov/ttn/scram/dispersion_prefrec.htm#aermod
C
C-----  MODIFIED BY:    Roger W. Brode
C                       U.S. EPA, OAQPS/AQAD
C                       Air Quality Modeling Group
C
C                       June 30, 2015
C
C-----  MODIFIED FROM:        AERMOD
C                       (Version Dated 14134)
C
C
C=======================================================================
C
C       This revised version of AERMOD (dated 14134) includes several
C       modifications relative to the previous version (dated 13350); 
C       see AERMOD Model Change Bullentin (MCB) #10 and AERMOD User's 
C       Guide Addendum.
C
C       http://www.epa.gov/ttn/scram/dispersion_prefrec.htm#aermod
C
C-----  MODIFIED BY:    Roger W. Brode
C                       U.S. EPA, OAQPS/AQAD
C                       Air Quality Modeling Group
C
C                       May 14, 2014
C
C-----  MODIFIED FROM:        AERMOD
C                       (Version Dated 13350)
C
C
C=======================================================================
C
C       This revised version of AERMOD (dated 13350) includes the
C       following modifications relative to the previous version
C       (dated 12345); see MCB#9 and AERMOD User's Guide Addendum.
C
C-----  Bug Fixes:
C
C       1.  Modified subroutine HRQEXT to remove the FATAL = .TRUE. 
C           statement for cases when the exit temp < 200K (about -100F) 
C           with non-zero emissions in the HOUREMIS file.  This case was 
C           intended to be treated as a non-fatal warning condition, but 
C           was inadvertently treated as a fatal error in version 12345 
C           by including the FATAL = .TRUE. statement. A fatal error 
C           message is still generated if the hourly exit temperature is 
C           less than 200K (about -100F) and the exit velocity is greater 
C           than 200m/s, unless the hourly emissions are zero (0) for 
C           that hour, since this would likely indicate that exit velocity
C           and exit temperature were swithed in the HOUREMIS file.
C
C       2.  Modified subroutine HRLOOP to correct a bug associated with 
C           use of the DAYRANGE keyword for multiple years of meteorological 
c           data, where the YR/MN/DY conversion to Julian may be incorrect in
C           relation to leap-years vs. non-leap-years. This problem was 
c           partially corrected in version 12345.
c
C       3.  Modified subroutine BGREAD to move the unit conversion for 
C           hourly background concentrations to follow the READ statements 
C           to avoid "double counting" unit conversion for non-hourly 
C           background since unit conversion for BGFILL has already been
C           applied in sub_BGVAL.
C
C       4.  Moved calculation of the center of effective area source for 
C           OPENPIT sources from subroutine ARDIST to subroutine PITEFF. 
C           Previous versions would have skipped calculation of center 
C           coordinates if the first receptor was located inside the 
C           actual OPENPIT source.
C
C       5.  Modified subroutine AREAIN to include an additional check 
C           for the number of sides exceeding the maximum number for an 
C           AREAPOLY source and to issue an AERMOD error message to avoid 
C           a compiler-generated runtime error.
C
C       6.  Modified subroutines wake_dfsn and wake_dfsn2 to correct a
C           problem with virtual sigma-z values for the cavity source,
C           which could result in anomalously high concentrations in
C           some cases due to the initial dispersion coefficient for 
C           the "cavity source" being set to zero (0).
C           
C       7.  Modified subroutine OLM_CALC to include a check for the 
C           OLM value (variable OLMVAL) being zero(0) to avoid a
C           zero-devide in calculating PercentNO2.
C
C       8.  Modified subroutine HEFF to include plume rise for penetrated
C           source, DHP3, for purposes of calculating plume height at 20km 
C           for use in wet deposition/depletion calculations.
C           
C       9.  Modified subroutine IBLVAL to use a minimum value of 5m for 
C           ZHI in the calculation of effective parameters.
C           
C       10. Modified subroutine O3READ to avoid including a substituted 
C           ozone (O3) value for hours with missing O3 data in the calculation 
C           of the maximum ozone value within the past 24 hours, which may be
C           used in the adjustment of the hourly O3 value under stable 
C           conditions.
C
C
C-----  Enhancements:
C 
C       1.  Incorporated new options for modeling NO2, including a new 
C           Default option for the Ambient Ratio Method (ARM) and a 
C           non-Default BETA option for the Ambient Ratio Method - 2 
C           (ARM2).
C       
C       2.  Incorporated options to vary background ozone (O3) data by 
C           wind sector (new CO O3SECTOR keyword), and options to vary 
C           background concentrations for the pollutant being modeled by 
C           wind sector (new SO BGSECTOR keyword).  Any of the existing 
C           options for specifying background data can be used with the
C           new sector-varying options, and will continue to work as 
C           before if no sectors are defined.  The applicable sector 
C           is determined by the flow vector (downwind) based on the 
C           wind direction in the surface meteorological data file.
C           
C       3.  Added new "debug" output file options for the OLM option
C           and for the new ARM and ARM2 options for NO2. Also added
C           a new PRIME debug option to separate the debug information
C           associated with the PRIME downwash algorithm from the 
C           non-PRIME related information provided under the MODEL
C           debug option.  Also removed all debug inforation from
C           the main 'aermod.out' file.
C
C
C-----  Miscellaneous:
C
C       1.  Modified module MAIN1 and subroutine SRCSIZ to include a
C           new NVPOLY parameter for the number of additional vertices 
C           to include for AREAPOLY sources, beyond the number of 
C           vertices used to define the source.  This is needed to 
C           account for the integrated portion of the area intersecting 
C           several sides of the polygon, depending on the geometry of 
C           the source, receptor, and wind direcion. The NVPOLY parameter 
C           is assigned a value of 12, which should work in most cases, 
C           but users are cautioned to avoid overly complex shapes for 
C           AREAPOLY sources to avoid model runtime errors. Reducing the 
C           complexity of AREAPOLY sources will also reduce model runime 
C           in most cases. The value of NVPOLY is added to NVMAX (determined 
C           based on the maximum number of area sources vertices defined 
C           by the user) in subroutine SRCSIZ for purposes of allocating 
C           arrays for AREAPOLY sources. The previous version of AERMOD 
C           added a value of 8 to NVMAX in subroutine SRCSIZ.
C           
C       2.  Incorporated checks on the consistency of concentrations that
C           are calculated during the EVENT processing mode with the original
C           value calculated during normal processing, which is included on 
C           the EVENTPER keyword. Warning messages are generated if the 
C           results are inconsistent, which may indicate a potential error
C           in the AERMOD code.
C
C       3.  Modified several subroutines to include the AERMET version date
C           associated with the input meteorological data in the header records
C           of the 'aermod.out' file and all other output file options. Also
C           modified the earliest acceptable version date for AERMET inputs
C           from 11059 to 12345. However, users are encouraged to use input
C           meteorological data based on the latest version of the AERMET 
C           processor (dated 13350).
C           
C       4.  Modified subroutines MEREAD, MEOPEN, and METEXT to extract met
C           data station IDs from the header record of the surface met data
C           file based on the index of keywords within the header record 
C           to provide more flexibility in the positioning of station IDs 
C           within the surface header record.
C
C       5.  Modified some error/warning messages to improve the clarity of the
C           message and to provide more flexibility in the error reporting.
C           For example:
C           
C           a. Errors associated with allocating arrays to have been modified 
C              to identify the type of array(s) being allocated when the
C              error occurs;
C           
C           b. Additional error handling and reporting has been incorporated 
C              for problems associated with opening, reading, and/or writing 
C              to external data files;
C           
C           c. Subroutine JULIAN was modified to provide additional information 
C              errors occur to faciliate diagnosing and correcting problems;
C           
C           d. Subroutines ERRHDL, TERRST, and other subroutines that call these
C              routines, were modified to include up to 12 characters in the module 
C              names included in error/warning/informational messages. The previous 
C              version limited the module name to 6 characters.
C
C
C-----  MODIFIED BY:    Roger W. Brode
C                       U.S. EPA, OAQPS/AQAD
C                       Air Quality Modeling Group
C
C                       December 16, 2013
C
C-----  MODIFIED FROM:        AERMOD
C                       (Version Dated 12345)
C
C
C
C=======================================================================
C
C----- NOTE: Summaries of AERMOD model changes prior to version 13350 
C            have been moved to the end of the 'aermod.f' file.
C
C=======================================================================
C
C     Variable Declarations
      USE MAIN1
      USE BUOYANT_LINE
      
      IMPLICIT NONE
      CHARACTER MODNAM*12

C --- Declare character strings for use in the header records of 
C     debug files for NO2 options
      CHARACTER (LEN=18) :: NO2DebugString
      CHARACTER (LEN=50) :: NO2DebugStrNonEvt, NO2DebugStrEvent

      INTEGER :: I, IDSTAT, IASTAT, IBSTAT
      LOGICAL :: L_OPENED, HIT_THRESH
      
C     Variable Initializations
      MODNAM = 'MAIN'
      NO2DebugString    = ''
      NO2DebugStrNonEvt = ''
      NO2DebugStrEvent  = ''
      
      NTOTHRS = 0
      NYEARS  = 0
      IBSTAT  = 0
      
      FATAL  = .FALSE.
      RUNERR = .FALSE.
      L_OPENED = .FALSE.
      HIT_THRESH = .FALSE.
      L_NoHeader(:) = .FALSE.

C     Open the Temporary File for Error Messages Generated from the Program
      OPEN(UNIT=IERUNT,FILE='ERRMSG.TMP',STATUS='REPLACE')

CCLC     Command line arguments removed.  Use AERMOD.INP and AERMOD.OUT.
CCLC     Retrieve Input and Output File Names From Command Line,
CCLC     ---   CALL GETCOM
CCL      CALL GETCOM (' AERMOD ',ILEN_FLD,INPFIL,OUTFIL)

C     Initialize TITLE1 and TITLE2 in case of fatal errors during PRESET phase
      TITLE1 = ' '
      TITLE2 = ' '
C     Get Run Date and Time using Fortran 90 functions      ---   CALL DATIME
      CALL DATIME (RUNDAT, RUNTIM)

C     Open Input and Output Files                           ---   CALL FILOPN
      CALL FILOPN

C     Preprocess Setup Information to Determine Data Storage Needs
      CALL PRESET

C --- Check for FATAL error occurring during PRESET processing, and
C     bypass further processing
      IF (FATAL) GOTO 12345
      
      IF (NYEARS .EQ. 0) THEN
C ---    Number of years for MAXDCONT arrays has not be specified,
C        set to default value of 5 years
         NYEARS = 5
      END IF
      
      IF (.NOT. EVONLY) THEN
C        OPEN The Temporary File to Store Events for EVENT File;
         OPEN(UNIT=ITEVUT,FILE='EVENT.TMP',STATUS='REPLACE')
C        Initialize the Event Counter
         IEVENT = 0
      END IF

C     Allocate SETUP Array Storage
      CALL ALLSETUP

      IF (ALLOC_ERR) THEN
C        Error occurred during allocation of Setup arrays.
C        Issue error messages and ABORT model run.
         WRITE(IOUNIT,*) ' '
         WRITE(IOUNIT,*) '  ERROR OCCURRED DURING ALLOCATION OF SETUP ',
     &                   'ARRAYS! ABORTING MODEL EXECUTION!'
         WRITE(IOUNIT,10901) NSRC,NGRP,NREC,NSEC,NQF,NBF,NPDMAX,NVMAX,
     &                       NURB,NOLM,NPSD,NNET,IXM,IYM,NAVE,NTYP,
     &                       nval,nhiann,nmax
10901    FORMAT(/'   ARRAY PARAMETER SETTINGS: ',/
     &           '         NSRC   = ', I8,/
     &           '         NGRP   = ', I8,/
     &           '         NREC   = ', I8,/
     &           '         NSEC   = ', I8,/
     &           '         NQF    = ', I8,/
     &           '         NBF    = ', I8,/
     &           '         NPDMAX = ', I8,/
     &           '         NVMAX  = ', I8,/
     &           '         NURB   = ', I8,/
     &           '         NOLM   = ', I8,/
     &           '         NPSD   = ', I8,/
     &           '         NNET   = ', I8,/
     &           '         IXM    = ', I8,/
     &           '         IYM    = ', I8,/
     &           '         NAVE   = ', I8,/
     &           '         NTYP   = ', I8,/
     &           '         NVAL   = ', I8,/
     &           '         NHIANN = ', I8,/
     &           '         NMAX   = ', I8)

         WRITE(IOUNIT,*)
         WRITE(IOUNIT,9057) STORE
 9057    FORMAT(/'   Estimated Storage Requirements of Model = ',
     &          F9.1,' MB of RAM.'/)

C        Write error message to terminal
         WRITE(*,*) ' '
         WRITE(*,*)'  ERROR OCCURRED DURING ALLOCATION OF SETUP ',
     &             'ARRAYS! ABORTING MODEL EXECUTION!'
         WRITE(*,9057) STORE

         GO TO 9999
      END IF

C     Variable Initializations                              ---   CALL VARINI
      CALL VARINI

C     Process The Model Setup Information                   ---   CALL SETUP
      IF (EVONLY) THEN
         CALL EV_SETUP
      ELSE
         CALL SETUP
      END IF

C --- Open files for NO2 Debug Options, including PVMRM, OLM, 
C     ARM, and ARM2, and write file header information. First
C     define char string explaining order of results for both
C     Non-EVENT and EVENT processing.
      NO2DebugStrNonEvt = 
     &   'Hrly non-EVT results sorted by hour/rec/srcgrp/src'
      NO2DebugStrEvent  =
     &   'Hrly & Ave EVT results sorted by day/event/grp/src'

      IF (PVMRM) THEN
         NO2DebugString = 'PVMRM Debug File: '
      ELSE IF (PVMRM2) THEN
         NO2DebugString = 'PVMRM2 Debug File:'
      ELSE IF (OLM) THEN
         NO2DebugString = 'OLM Debug File:   '
      ELSE IF (ARM) THEN
         NO2DebugString = 'ARM Debug File:   '
      ELSE IF (ARM2) THEN
         NO2DebugString = 'ARM2 Debug File:  '
      END IF

c --- Text included in second header record of NO2 option debug files:
c Non-Event runs:
c OLM Debug File:   "Hrly non-EVT results sorted by hour/rec/srcgrp/src"
c ARM Debug File:   "Hrly non-EVT results sorted by hour/rec/srcgrp/src"
c ARM2 Debug File:  "Hrly non-EVT results sorted by hour/rec/srcgrp/src"
c PVMRM Debug File: "Hrly non-EVT results sorted by hour/rec/srcgrp/src"
c Event runs:
c OLM Debug File:   "Hrly & Ave EVT results sorted by day/event/grp/src"
c ARM Debug File:   "Hrly & Ave EVT results sorted by day/event/grp/src"
c ARM2 Debug File:  "Hrly & Ave EVT results sorted by day/event/grp/src"
c PVMRM Debug File: "Hrly & Ave EVT results sorted by day/event/grp/src"

C --- Open file with PVMRM debugging output, but first check 
C     for potential file unit conflict
      IF (PVMRMDBG) THEN
         L_OPENED = .FALSE.
         INQUIRE (FILE=DBPVFIL,OPENED=L_OPENED)
         IF( .NOT.L_OPENED )THEN
            INQUIRE (UNIT=PVMDBG,OPENED=L_OPENED)
            IF (.NOT.L_OPENED) THEN
              OPEN(UNIT=PVMDBG,FILE=DBPVFIL,STATUS='REPLACE')
C             Write the standard header information with AERMOD/AERMET 
C             versions, 1st title, and rundat/runtim to the debug output file
              WRITE(PVMDBG,9028) VERSN, TITLE1(1:68), RUNDAT
 9028         FORMAT('*** AERMOD - VERSION ',A6,' ***',3X,'*** ',A68,
     &                                                  ' ***',8X,A8)
              IF( .NOT.EVONLY )THEN
C ---          Write header record with column labels for Non-EVENT results
               WRITE(PVMDBG,9029) C_METVER, NO2DebugString, 
     &                            NO2DebugStrNonEvt, RUNTIM
 9029          FORMAT('*** AERMET - VERSION ',A6,' ***',3X,'*** ',
     &                                  A18, A50,' ***',8X,A8/)
C              Write the model options (MODOPS) to the debug output file
               WRITE ( PVMDBG, 200 ) 
     &                        MODOPS_String(1:LEN_TRIM(MODOPS_String))
  200          FORMAT ( ' OPTIONS: ', A /)
              
               WRITE(PVMDBG,9001)
 9001          FORMAT(1x,'TYPE    DOM_SRCID     DATE    IREC    GRPID',
     &                '   ISRC    SRCID       NUMCONT     DISTDOM     ',
     &                'MAXCONC_NOx     O3CONC        O3MOLES      ',
     &                'NOxMOLES       BHORIZ        BVERT         ',
     &                'BVERT3       ', 
     &                'PLUMEVOL      OrigConc  X  PercentNO2 =    ',
     &                'HRVAL         AVEVAL')
              ELSE
C ---          Write header record with column labels for EVENT results
               WRITE(PVMDBG,9029) C_METVER, NO2DebugString, 
     &                            NO2DebugStrEvent, RUNTIM
C              Write the model options (MODOPS) to the debug output file
               WRITE ( PVMDBG, 200 ) 
     &                        MODOPS_String(1:LEN_TRIM(MODOPS_String))
               WRITE(PVMDBG,90011)
90011          FORMAT(1x,'TYPE    DOM_SRCID     DATE    IEVE     ',
     &                'EVENTID    IAVE    GRPID   ISRC    SRCID     ',
     &                'NUMCONT   DISTDOM     MAXCONC_NOx     O3CONC   ',
     &                '     O3MOLES      NOxMOLES       BHORIZ        ',
     &                'BVERT         BVERT3       PLUMEVOL      ',
     &                'OrigConc  X  PercentNO2 =    ',
     &                'HRVAL         AVEVAL')
              ENDIF

            ELSE
C ---          Unit is already opened, issue error message
               CALL ERRHDL(PATH,MODNAM,'E','501','PVMRMDBG')
            END IF
         ELSE
C ---       Unit is already opened, issue error message
            CALL ERRHDL(PATH,MODNAM,'E','501','PVMRMDBG')
         END IF
      END IF

C --- Add separate debug file for relative dispersion coefficients for PVMRM/PVMRM2 options;
C     the RELDISP debug file will automatically be named RELDISP.DBG
      IF (PVMRMDBG) THEN
         L_OPENED = .FALSE.
         INQUIRE (FILE=RDISPFIL,OPENED=L_OPENED)
         IF( .NOT.L_OPENED )THEN
            INQUIRE (UNIT=RDISPUNT,OPENED=L_OPENED)
            IF (.NOT.L_OPENED) THEN
              OPEN(UNIT=RDISPUNT,FILE=RDISPFIL,STATUS='REPLACE')
C             Write the standard header information with AERMOD/AERMET 
C             versions, 1st title, and rundat/runtim to the debug output file
              WRITE(RDISPUNT,9028) VERSN, TITLE1(1:68), RUNDAT
              IF( .NOT.EVONLY )THEN
C ---          Write header record with column labels for Non-EVENT results
               WRITE(RDISPUNT,9029) C_METVER, NO2DebugString, 
     &                            NO2DebugStrNonEvt, RUNTIM
C              Write the model options (MODOPS) to the debug output file
               WRITE ( RDISPUNT, 200 ) 
     &                        MODOPS_String(1:LEN_TRIM(MODOPS_String))

C ---          PVMRM and MODEL DEBUG options have been selected; include the
C              RELDISP Calcs in the RELDISP DEBUG file

              ELSE
C ---          Write header record with column labels for EVENT results
               WRITE(RDISPUNT,9029) C_METVER, NO2DebugString, 
     &                            NO2DebugStrEvent, RUNTIM
C              Write the model options (MODOPS) to the debug output file
               WRITE ( RDISPUNT, 200 ) 
     &                        MODOPS_String(1:LEN_TRIM(MODOPS_String))

C ---          PVMRM and MODEL DEBUG options have been selected; include the
C              RELDISP Calcs in the RELDISP DEBUG file

              END IF
            ELSE
C ---          Unit is already opened, issue error message
               CALL ERRHDL(PATH,MODNAM,'E','501','PVMRMDBG')
            END IF
         ELSE
C ---       Unit is already opened, issue error message
            CALL ERRHDL(PATH,MODNAM,'E','501','PVMRMDBG')
         END IF
      END IF


C --- Open file with OLM debugging output, but first check 
C     for potential file unit conflict
      IF (OLMDEBUG) THEN
         L_OPENED = .FALSE.
         INQUIRE (FILE=DBOLMFIL,OPENED=L_OPENED)
         IF (.NOT.L_OPENED) THEN
            INQUIRE (UNIT=OLMDBG,OPENED=L_OPENED)
            IF (.NOT.L_OPENED) THEN
              OPEN(UNIT=OLMDBG,FILE=DBOLMFIL,STATUS='REPLACE')
C             Write the standard header information with AERMOD/AERMET 
C             versions, 1st title, and rundat/runtim to the debug output file
              WRITE(OLMDBG,9028) VERSN, TITLE1(1:68), RUNDAT
              IF( .NOT.EVONLY )THEN
                WRITE(OLMDBG,9029) C_METVER, NO2DebugString, 
     &                             NO2DebugStrNonEvt, RUNTIM
C               Write the model options (MODOPS) to the debug output file
                WRITE ( OLMDBG, 200 ) 
     &                        MODOPS_String(1:LEN_TRIM(MODOPS_String))
C ---           Write header record with column labels for Non-EVENT results
                WRITE(OLMDBG,9002)
9002            FORMAT(3X,'DATE    IREC    GRPID   ISRC    SRCID',7X,
     &          'IOLM    OLMID     O3CONC',7X,'OLMVAL',3X,
     &          'X   NO2Ratio  =  NO2VAL',7X,'NO_VAL',6X,'OrigConc  ',
     &          'X  PercentNO2 =  HRVAL',8X,'AVEVAL')

              ELSE
                WRITE(OLMDBG,9029) C_METVER, NO2DebugString, 
     &                             NO2DebugStrEvent, RUNTIM
C ---            Write header record with column labels for EVENT results
C                Write the model options (MODOPS) to the debug output file
                 WRITE ( OLMDBG, 200 ) 
     &                        MODOPS_String(1:LEN_TRIM(MODOPS_String))
                 WRITE(OLMDBG,90021)
90021            FORMAT(3X,'DATE    IEVE   EVENTID   IAVE    GRPID',
     &           3x,'ISRC',4x,'SRCID       IOLM	 OLMID       O3CONC',7x,
     &           'OLMVAL   X   NO2Ratio  =  NO2VAL       NO_VAL',6x,
     &           'OrigConc  X  PercentNO2 =  HRVAL       EV_AVEVAL')

             END IF
            ELSE
C ---          Unit is already opened, issue error message
               CALL ERRHDL(PATH,MODNAM,'E','501','OLMDEBUG')
            END IF
         ELSE
C ---       Unit is already opened, issue error message
            CALL ERRHDL(PATH,MODNAM,'E','501','OLMDEBUG')
         END IF
      END IF

C --- Open file with ARM debugging output, but first check 
C     for potential file unit conflict
      IF (ARMDEBUG) THEN
         L_OPENED = .FALSE.
         INQUIRE (FILE=DBARMFIL,OPENED=L_OPENED)
         IF (.NOT.L_OPENED) THEN
            INQUIRE (UNIT=ARMDBG,OPENED=L_OPENED)
            IF (.NOT.L_OPENED) THEN
              OPEN(UNIT=ARMDBG,FILE=DBARMFIL,STATUS='REPLACE')
C             Write the standard header information with AERMOD/AERMET 
C             versions, 1st title, and rundat/runtim to the debug output file
              WRITE(ARMDBG,9028) VERSN, TITLE1(1:68), RUNDAT

              IF( .NOT.EVONLY )THEN
                WRITE(ARMDBG,9029) C_METVER, NO2DebugString, 
     &                             NO2DebugStrNonEvt, RUNTIM
C               Write the model options (MODOPS) to the debug output file
                WRITE ( ARMDBG, 200 ) 
     &                        MODOPS_String(1:LEN_TRIM(MODOPS_String))
C ---           Write header record with column labels for Non-EVENT results
                WRITE(ARMDBG,9003)
9003            FORMAT(3X,'DATE',4x,'IREC    GRPID   ISRC    SRCID',11x,
     &                   'OrigConc X  ARM_Ratio =    HRVAL',7X,'AVEVAL')
              ELSE
                WRITE(ARMDBG,9029) C_METVER, NO2DebugString, 
     &                             NO2DebugStrEvent, RUNTIM
C               Write the model options (MODOPS) to the debug output file
                WRITE ( ARMDBG, 200 ) 
     &                        MODOPS_String(1:LEN_TRIM(MODOPS_String))
C ---           Write header record with column labels for EVENT results
                WRITE(ARMDBG,90031)
90031           FORMAT(3X,'DATE',4X,'IEVE   EVENTID   IAVE    GRPID',3X,
     &               'ISRC    SRCID',11X,'OrigConc X  ARM_Ratio =',4X,
     &               'HRVAL',6X,'EV_AVEVAL')
              END IF

            ELSE
C ---          Unit is already opened, issue error message
               CALL ERRHDL(PATH,MODNAM,'E','501','ARMDEBUG')
            END IF
         ELSE
C ---       Unit is already opened, issue error message
            CALL ERRHDL(PATH,MODNAM,'E','501','ARMDEBUG')
         END IF
      END IF

C --- Open file with ARM2 debugging output, but first check 
C     for potential file unit conflict
      IF (ARM2DEBUG) THEN
         L_OPENED = .FALSE.
         INQUIRE (FILE=DBARM2FIL,OPENED=L_OPENED)
         IF (.NOT.L_OPENED) THEN
           INQUIRE (UNIT=ARMDBG,OPENED=L_OPENED)
           IF (.NOT.L_OPENED) THEN
             OPEN(UNIT=ARMDBG,FILE=DBARM2FIL,STATUS='REPLACE')
C            Write the standard header information with AERMOD/AERMET 
C            versions, 1st title, and rundat/runtim to the debug output file
             WRITE(ARMDBG,9028) VERSN, TITLE1(1:68), RUNDAT

             IF( .NOT.EVONLY )THEN
               WRITE(ARMDBG,9029) C_METVER, NO2DebugString, 
     &                            NO2DebugStrNonEvt, RUNTIM
C              Write the model options (MODOPS) to the debug output file
               WRITE ( ARMDBG, 200 ) 
     &                        MODOPS_String(1:LEN_TRIM(MODOPS_String))
C ---          Write header record with column labels for Non-EVENT results
               WRITE(ARMDBG,9004)
9004           FORMAT(3X,'DATE',4x,'IREC    GRPID   ISRC    SRCID',11x,
     &                   'NOXCONC',5x,'OrigConc  X  ARM2Ratio =',3x,
     &                   'HRVAL',8x,'AVEVAL')
             ELSE
               WRITE(ARMDBG,9029) C_METVER, NO2DebugString, 
     &                            NO2DebugStrEvent, RUNTIM
C              Write the model options (MODOPS) to the debug output file
               WRITE ( ARMDBG, 200 ) 
     &                        MODOPS_String(1:LEN_TRIM(MODOPS_String))
C ---          Write header record with column labels for EVENT results
               WRITE(ARMDBG,90041)
90041          FORMAT(3X,'DATE',4x,'IEVE   EVENTID   IAVE    GRPID',3x,
     &              'ISRC',4x,'SRCID',11x,'NOXCONC',6x,'OrigConc X  ',
     &              'ARM2Ratio =',4x,'HRVAL',6x,'EV_AVEVAL')
             END IF
           ELSE
C ---         Unit is already opened, issue error message
              CALL ERRHDL(PATH,MODNAM,'E','501','ARM2DEBUG')
           END IF
         ELSE
C ---       Unit is already opened, issue error message
            CALL ERRHDL(PATH,MODNAM,'E','501','ARM2DEBUG')
         END IF
      END IF

C --- Write the model options and debug data template to the
C     debug file if MODEL is specified on the DEBUGOPT keyword
      IF( DEBUG )THEN
C        Write the standard header information with AERMOD/AERMET 
C        versions, 1st title, and rundat/runtim to the debug output file
         WRITE(DBGUNT,9028) VERSN, TITLE1(1:68), RUNDAT
         IF( EVONLY )THEN
C ---       Write header information for EVENT results
            WRITE(DBGUNT,90291) C_METVER,
     &                  'MODEL Debug File: EVENT Processing   ', RUNTIM
         ELSE
C ---       Write header information for Non-EVENT results
            WRITE(DBGUNT,90291) C_METVER,
     &                  'MODEL Debug File: NonEVENT Processing', RUNTIM
         ENDIF
90291    FORMAT('*** AERMET - VERSION ',A6,' ***',3X,'*** ',
     &                             A37,31X,' ***',8X,A8/)
C ---    Include string of MODEL Options in DEBUG file header
         WRITE (DBGUNT, 200)  MODOPS_String(1:LEN_TRIM(MODOPS_String))
      END IF

C --- Write the model options and debug data template to the
C     debug file if METEOR is specified on the DEBUGOPT keyword
      IF( METEORDBG )THEN 
         WRITE(DBMUNT,9028) VERSN, TITLE1(1:68), RUNDAT
         IF( EVONLY )THEN
C ---       Write header information for EVENT results
            WRITE(DBMUNT,90292) C_METVER,
     &              'METEOR Debug File: EVENT Processing   ', RUNTIM
         ELSE
C ---       Write header information for Non-EVENT results
            WRITE(DBMUNT,90292) C_METVER,
     &              'METEOR Debug File: NonEVENT Processing', RUNTIM
         ENDIF
90292    FORMAT('*** AERMET - VERSION ',A6,' ***',3X,'*** ',
     &                                 A38,30X,' ***',8X,A8/)
C ---    Include string of MODEL Options in DEBUG file header
         WRITE (DBMUNT, 200)  MODOPS_String(1:LEN_TRIM(MODOPS_String))
      END IF

C --- Write the model options and debug data template to the
C     debug file if PRIME is specified on the DEBUGOPT keyword
C     and model run includes sources with building downwash (NSEC > 0)
      IF (PRIMEDBG .AND. NSEC .GT. 0) THEN
C ---    Write main header records for PRIME debug file
         WRITE(PRMDBUNT,9028) VERSN, TITLE1(1:68), RUNDAT
         IF( EVONLY )THEN
C ---       Write header information for EVENT results
            WRITE(PRMDBUNT,90293) C_METVER, 
     &                  'PRIME Debug File: EVENT Processing   ', RUNTIM
         ELSE
C ---       Write header information for NonEVENT results
            WRITE(PRMDBUNT,90293) C_METVER, 
     &                  'PRIME Debug File: NonEVENT Processing', RUNTIM
         ENDIF
90293    FORMAT('*** AERMET - VERSION ',A6,' ***',3X,'*** ',
     &                              A37,31X,' ***',8X,A8/)
C        Write the model options (MODOPS) to the debug output file
         WRITE ( PRMDBUNT,200 ) MODOPS_String(1:LEN_TRIM(MODOPS_String))
      END IF

C --- Write the model options and debug data template to the
C     debug file if AREA/LINE is specified on the DEBUGOPT keyword
C     and model run includes AREA, LINE, or OPENPIT sources
      IF( AREADBG .AND. (NAREA.GT.0 .OR. NCIRC.GT.0 .OR. 
     &                   NLINE.GT.0 .OR. NPIT.GT.0) )THEN
C ---    Write main header records for METEOR debug file, then reset METDBGHDR = .F.
         WRITE(AREADBUNT,9028) VERSN, TITLE1(1:68), RUNDAT
         IF( EVONLY )THEN
C ---       Write header information for EVENT results
            WRITE(AREADBUNT,90294) C_METVER, 
     &                    'AREA Debug File: EVENT Processing   ', RUNTIM
         ELSE
C ---       Write header information for NonEVENT results
            WRITE(AREADBUNT,90294) C_METVER, 
     &                    'AREA Debug File: NonEVENT Processing', RUNTIM
         ENDIF
90294    FORMAT('*** AERMET - VERSION ',A6,' ***',3X,'*** ',
     &                                 A36,32X,' ***',8X,A8/)
C        Write the model options (MODOPS) to the debug output file
         WRITE( AREADBUNT,200 ) MODOPS_String(1:LEN_TRIM(MODOPS_String))
      END IF

C     Open file for GDEP output from gas dry deposition algorithms,
C     but first check for potential file unit conflict
      IF ((DEBUG .OR. DEPOSDBG) .AND. LDGAS) THEN
         L_OPENED = .FALSE.
         INQUIRE (FILE='GDEP.DAT',OPENED=L_OPENED)
         IF (.NOT.L_OPENED) THEN
            INQUIRE (UNIT=GDEPDBG,OPENED=L_OPENED)
            IF (.NOT.L_OPENED) THEN
               OPEN(UNIT=GDEPDBG,FILE='GDEP.DAT',STATUS='REPLACE')
            ELSE
C ---          Unit is already opened, issue error message
               CALL ERRHDL(PATH,MODNAM,'E','501','GDEP.DAT')
               RUNERR = .TRUE.
            END IF
         ELSE
C ---       Unit is already opened, issue error message
            CALL ERRHDL(PATH,MODNAM,'E','501','GDEP.DAT')
            RUNERR = .TRUE.
         END IF
      END IF

C     Open file for PDEP output from particle dry deposition algorithms,
C     but first check for potential file unit conflict
      IF ((DEBUG .OR. DEPOSDBG) .AND. LDPART) THEN
         L_OPENED = .FALSE.
         INQUIRE (FILE='PDEP.DAT',OPENED=L_OPENED)
         IF (.NOT.L_OPENED) THEN
            INQUIRE (UNIT=PDEPDBG,OPENED=L_OPENED)
            IF (.NOT.L_OPENED) THEN
               OPEN(UNIT=PDEPDBG,FILE='PDEP.DAT',STATUS='REPLACE')
            ELSE
C ---          Unit is already opened, issue error message
               CALL ERRHDL(PATH,MODNAM,'E','501','PDEP.DAT')
               RUNERR = .TRUE.
            END IF
         ELSE
C ---       Unit is already opened, issue error message
            CALL ERRHDL(PATH,MODNAM,'E','501','PDEP.DAT')
            RUNERR = .TRUE.
         END IF
      END IF

C     Deallocate Temporary Storage
      DEALLOCATE  (IWRK2, STAT=IDSTAT)
      IF (IDSTAT .NE. 0) THEN
         CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
      END IF
      IF (.NOT. EVONLY) THEN
         DEALLOCATE  (ZETMP1,ZETMP2,ZHTMP1,ZHTMP2,ZFTMP1,ZFTMP2,
     &                STAT=IDSTAT)
         IF (IDSTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
         END IF
      END IF

C     Allocate Array Storage for Results                    ---   CALL ALLRESULT
      CALL ALLRESULT

      IF (ALLOC_ERR) THEN
C        Error occurred during allocation of Results arrays.
C        Issue error message and skip initialization of results arrays.
         WRITE(IOUNIT,*) ' '
         WRITE(IOUNIT,*) '  ERROR OCCURRED DURING ALLOCATION OF RESULT',
     &                   ' ARRAYS!'
         WRITE(IOUNIT,10902) NSRC,NGRP,NREC,NSEC,NQF,NBF,NPDMAX,NVMAX,
     &                       NURB,NOLM,NPSD,NNET,IXM,IYM,NAVE,NTYP,
     &                       NHIVAL,NHIANN,NMXVAL,NYEARS
10902    FORMAT(/'   ARRAY PARAMETER SETTINGS: ',/
     &           '         NSRC   = ', I8,/
     &           '         NGRP   = ', I8,/
     &           '         NREC   = ', I8,/
     &           '         NSEC   = ', I8,/
     &           '         NQF    = ', I8,/
     &           '         NBF    = ', I8,/
     &           '         NPDMAX = ', I8,/
     &           '         NVMAX  = ', I8,/
     &           '         NURB   = ', I8,/
     &           '         NOLM   = ', I8,/
     &           '         NPSD   = ', I8,/
     &           '         NNET   = ', I8,/
     &           '         IXM    = ', I8,/
     &           '         IYM    = ', I8,/
     &           '         NAVE   = ', I8,/
     &           '         NTYP   = ', I8,/
     &           '         NHIVAL = ', I8,/
     &           '         NHIANN = ', I8,/
     &           '         NMXVAL = ', I8,/
     &           '         NYEARS = ', I8)

         WRITE(IOUNIT,*)
         WRITE(IOUNIT,9057) STORE

C        Write error message to terminal
         WRITE(*,*) ' '
         WRITE(*,*) '  ERROR OCCURRED DURING ALLOCATION OF RESULT',
     &              ' ARRAYS!'
         WRITE(*,9057) STORE

         GO TO 9999

      ELSE IF (.NOT. EVONLY) THEN
C        No Errors During Allocation of Results Arrays
C        Initialize Results Arrays With Zeroes              ---   CALL RESINI
         CALL RESINI
      END IF

C     Determine Number of Setup Messages by Message Type    ---   CALL TERRST
      CALL TERRST

c --- Set up common for PRIME numerical rise algorithm      ---   CALL NUMPR1
      CALL NUMPR1

c --- Set up common for PRIME building cavity model         ---   CALL PRIME1
      CALL PRIME1

      IF (.NOT.RUN .OR. FATAL .OR. IWRN .GT. 0) THEN
C        Write Out Summary Of Setup Error/Message Stats     ---   CALL SUMTBL
         WRITE(IOUNIT,9111)
 9111    FORMAT(//2X,'*** Message Summary For AERMOD Model Setup ***'/)
         CALL SUMTBL(IOUNIT)
      END IF

12345 continue

      IF (FATAL) THEN
         WRITE(*,99111)
99111    FORMAT('+','Fatal Error Occurred During Setup Phase!')
         WRITE(IOUNIT,9112)
 9112    FORMAT(/4X,'**************************************',
     &          /4X,'*** SETUP Finishes UN-successfully ***',
     &          /4X,'**************************************'/)
      ELSE
         WRITE(IOUNIT,9113)
 9113    FORMAT(/1X,'***********************************',
     &          /1X,'*** SETUP Finishes Successfully ***',
     &          /1X,'***********************************'/)

C        Print Summary of the Input Data                       ---   CALL INPSUM
         CALL INPSUM

C        Write Headers to GDEP.DAT and PDEP.DAT Files for new deposition algorithms
         IF ((DEBUG .OR. DEPOSDBG) .AND. LDGAS) THEN
C           Write the model options (MODOPS) to the debug output file
            WRITE ( GDEPDBG, 200 ) 
     &                      MODOPS_String(1:LEN_TRIM(MODOPS_String))
            WRITE(GDEPDBG,9901)
 9901       FORMAT(1X,'YYMMDDHH',3X,'ISRC',4X,'Ra',12X,'Rb',12X,'Rc',
     &             12X,'Vdepg')
         END IF
         IF ((DEBUG .OR. DEPOSDBG) .AND. LDPART) THEN
            WRITE(PDEPDBG,9902)
C           Write the model options (MODOPS) to the debug output file
            WRITE ( PDEPDBG, 200 ) 
     &                      MODOPS_String(1:LEN_TRIM(MODOPS_String))
 9902       FORMAT(1X,'YYMMDDHH',3X,'ISRC',1X,'ICAT',2X,'Method No.',
     &             3X,'Ra',12X,'Rp',12X,'Vg(i)',9x,'Vdep(i)')
         END IF

      END IF

      IF (.NOT.FATAL .AND. RUN .AND. EVONLY) THEN
C        No Fatal Errors in Setup and RUN Option Selected and EVENT Processing

C        Process The Data For Each Event                    ---   CALL EVLOOP
         CALL EVLOOP

      ELSE IF (.NOT.FATAL .AND. RUN .AND. .NOT.EVONLY) THEN
C        No Fatal Errors in Setup and RUN Option Selected and Normal Processing

C        Reinitialize Results Arrays With Zeroes            ---   CALL RESINI
         CALL RESINI

         IF (RSTINP) THEN
C           Initialize Results Arrays from Re-start File    ---   CALL RSINIT
            CALL RSINIT
         END IF

C        Process The Hourly Meteorological Data             ---   CALL HRLOOP
         CALL HRLOOP

C ---    Check total precipitation if wet deposition is being used
         IF ((WDPLETE .OR. DEPOS .OR. WDEP) .AND. 
     &                                TOTAL_PRECIP .LT. 0.0001D0) THEN
C ---       Write warning message for no precip with wet deposition
            CALL ERRHDL(PATH,MODNAM,'W','496','WetDepos')
         END IF
         
         IF ((PM25AVE .OR. NO2AVE .OR. SO2AVE .OR. ANNUAL)
     &                                           .AND. MULTYR
     &                                           .AND. .NOT.RUNERR) THEN
C ---       Results arrays for MULTYEAR applications WITH ANNUAL average,
C           or other outputs averaged across years, need to be "dumped" to 
C           SAVEFILE BEFORE calculating averages
C ---                                                       ---   CALL RSDUMP
            CALL RSDUMP

            IF (SEASONHR .AND. .NOT.RUNERR) THEN
C ---          Calculate averages for season by hour-of-day results
               IF (CONC) THEN
                  CALL SHAVE
C ---             Check for values exceeding fixed-format field width (F13.8)
C                 without FILE_FORMAT = 'EXP'
                  IF (FILE_FORMAT .NE. 'EXP' .AND. 
     &                         MAXVAL(SHVALS) .GT. 9999.99999999D0) THEN
                     CALL ERRHDL(PATH,MODNAM,'W','400','= EXP')
                  END IF
               END IF
            END IF

         END IF

         IF ( (PM25AVE .OR. NO2AVE .OR. SO2AVE .OR. ANNUAL)
     &                                           .AND. .NOT.RUNERR) THEN
C ---       Compute averages of the High-N-High 24-hr PM25, 1-hr NO2, 
C           1-hr SO2, and annual values
            IF (NUMYRS .GT. 0) THEN
               DO IGRP = 1, NUMGRP
                  DO IREC = 1, NUMREC
                     IF (PM25AVE .OR. NO2AVE .OR. SO2AVE) THEN
                        SUMHNH(IREC,IGRP,1:NHIVAL) = 
     &                  SUMHNH(IREC,IGRP,1:NHIVAL) / DBLE(NUMYRS)
                     END IF
                     IF (ANNUAL) THEN
                        DO ITYP = 1, NUMTYP
                           ANNVAL(IREC,IGRP,ITYP) =
     &                                 SUMANN(IREC,IGRP,ITYP) / 
     &                                                      DBLE(NUMYRS)
                        END DO
                     END IF
                  END DO
               END DO
            ELSE
               IF( ANNUAL )THEN
C                 Write Error Message: Number of Years = 0 for ANNUAL
                  CALL ERRHDL(PATH,MODNAM,'E','480','NUMYRS=0')
                  RUNERR = .TRUE.
               ELSEIF( NO2AVE )THEN
C                 Write Error Message: Need complete years for 1-hr SO2/NO2 and 24-hr PM25
                  CALL ERRHDL(PATH,MODNAM,'E','480','1hr NO2AVE')
                  RUNERR = .TRUE.
               ELSEIF( SO2AVE )THEN
C                 Write Error Message: Need complete years for 1-hr SO2/NO2 and 24-hr PM25
                  CALL ERRHDL(PATH,MODNAM,'E','480','1hr SO2AVE')
                  RUNERR = .TRUE.
               ELSEIF( PM25AVE )THEN
C                 Write Error Message: Need complete years for 1-hr SO2/NO2 and 24-hr PM25
                  CALL ERRHDL(PATH,MODNAM,'E','480','24hr PM25AVE')
                  RUNERR = .TRUE.
               END IF
            END IF
            IF (NREMAIN .NE. 0) THEN
C              Write Warning Message: Met Data Remains After End of Last Year
               IF (.NOT. L_SkipMessages) THEN
                  WRITE(DUMMY,'(I8)') NREMAIN
                  CALL ERRHDL(PATH,MODNAM,'W','481',DUMMY)
               END IF
            END IF
         END IF

         IF ((PERIOD.OR.ANNUAL) .AND. (.NOT. RUNERR) .AND.
     &                                                NTOTHRS.GT.0) THEN
C ---       PERIOD Average Selected and No Runtime/Meteorology Errors
            IF (CONC .AND. PERIOD) THEN
C              Calculate Period Average Concentrations      ---   CALL PERAVE
               CALL PERAVE
            END IF
C ---       Check for values exceeding fixed-format field width (F13.5)
C           without FILE_FORMAT = 'EXP'
            IF (FILE_FORMAT .NE. 'EXP' .AND. 
     &                        MAXVAL(ANNVAL) .GT. 9999999.99999D0) THEN
               CALL ERRHDL(PATH,MODNAM,'W','400','= EXP')
            END IF
            DO ITYP = 1, NUMTYP
C              Select Highest PERIOD Values by Source Group ---   CALL HIPER
               CALL HIPER
            END DO
            IF (ANPOST) THEN
C              Write PERIOD/ANNUAL Results to Post File     ---   CALL PSTANN
               CALL PSTANN
            END IF
            IF (ANPLOT) THEN
C              Write PERIOD/ANNUAL Results to Plot File     ---   CALL PLTANN
               CALL PLTANN
            END IF
         END IF

         IF (MULTYR .AND. .NOT.RUNERR .AND. 
     &                    .NOT.(ANNUAL .OR. PM25AVE .OR. NO2AVE .OR. 
     &                                       SO2AVE)) THEN
C ---       Results arrays for MULTYEAR applications WITHOUT ANNUAL average,
C           or other outputs averaged across years, need to be "dumped" to 
C           SAVEFILE AFTER calculating averages
C ---                                                       ---   CALL RSDUMP
            CALL RSDUMP
         END IF

         IF (.NOT.(MULTYR .AND. (ANNUAL .OR. PM25AVE .OR. 
     &                                        NO2AVE .OR. 
     &                                        SO2AVE) ) .AND.
     &                         SEASONHR .AND. .NOT.RUNERR) THEN
            IF (CONC) THEN
               CALL SHAVE
C ---          Check for values exceeding fixed-format field width (F13.8)
C              without FILE_FORMAT = 'EXP'
               IF (FILE_FORMAT .NE. 'EXP' .AND. 
     &                         MAXVAL(SHVALS) .GT. 9999.99999999D0) THEN
                  CALL ERRHDL(PATH,MODNAM,'W','400','= EXP')
               END IF
            END IF
         END IF

         IF (PLFILE .AND. (.NOT. RUNERR)) THEN
C           Write Short Term High Values to Plot File       ---   CALL PLOTFL
C ---       Check for values exceeding fixed-format field width (F13.5)
C           without FILE_FORMAT = 'EXP'
            IF (FILE_FORMAT .NE. 'EXP') THEN
               IF (PM25AVE .OR. NO2AVE .OR. SO2AVE) THEN
                  IF (MAXVAL(SUMHNH) .GT. 9999999.99999D0) THEN
                     CALL ERRHDL(PATH,MODNAM,'W','400','= EXP')
                  END IF
               ELSE IF (.NOT.PM25AVE .AND. .NOT.NO2AVE .AND. .NOT.SO2AVE
     &                      .AND.MAXVAL(HIVALU).GT.9999999.99999D0) THEN
                  CALL ERRHDL(PATH,MODNAM,'W','400','= EXP')
               END IF
            END IF
C ---       Call plotfile routine
            CALL PLOTFL
         END IF

         IF (.NOT. RUNERR) THEN
C ---       Check for values exceeding fixed-format field width (F13.5)
C           without FILE_FORMAT = 'EXP'
            IF (.NOT.PLFILE .AND. FILE_FORMAT .NE. 'EXP') THEN
               IF (PM25AVE .OR. NO2AVE .OR. SO2AVE) THEN
                  IF (MAXVAL(SUMHNH) .GT. 9999999.99999D0) THEN
                     CALL ERRHDL(PATH,MODNAM,'W','400','= EXP')
                  END IF
               ELSE IF (.NOT.PM25AVE .AND. .NOT.NO2AVE .AND. .NOT.SO2AVE
     &                      .AND.MAXVAL(HIVALU).GT.9999999.99999D0) THEN
                  CALL ERRHDL(PATH,MODNAM,'W','400','= EXP')
               END IF
            END IF
C           Print Out Model Results                         ---   CALL OUTPUT
            CALL OUTPUT
         END IF

C ---    Check for MAXDCONT options to evaluate source group contributions 
C        based on rank for PM2.5 24hr, NO2 1hr or SO2 1hr NAAQS
         IF (.NOT. RUNERR .AND. 
     &         L_MAXDCONT .AND. 
     &            (PM25AVE .OR. NO2AVE .OR. SO2AVE) ) THEN

            IF (PVMRMDBG) THEN
C ---          PVMRM Debug option selected; print header record to delimit
C              debug information related to MAXDCONT processing
               WRITE(PVMDBG,9001)
            END IF

C ---       Allocate arrays to save receptor data;
C           also allocate array to store summed
C           contributions for max daily 1-hour averages
            ALLOCATE  (AXR_SAV(NREC), AYR_SAV(NREC), 
     &           AZELEV_SAV(NREC), AZFLAG_SAV(NREC), 
     &           AZHILL_SAV(NREC), 
     &           SUMVAL_MAXD(NVAL,NGRP,NGRP,NREC), 
     &           STAT=IASTAT)
            IF (IASTAT .NE. 0) THEN
               CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
               ALLOC_ERR = .TRUE.
               WRITE(IOUNIT,*) '  Error Occurred During Allocation ',
     &                         'Rec Arrays for MAXDCONT!'
            END IF

            IF (L_BLSOURCE) THEN
              ALLOCATE (XR_SCS_SAV(NREC), YR_SCS_SAV(NREC), STAT=IBSTAT)
            ENDIF
            IF (IBSTAT .NE. 0) THEN
               CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
               ALLOC_ERR = .TRUE.
               WRITE(IOUNIT,*) '  Error Occurred During Allocation ',
     &                         'BuoyLine Rec Arrays for MAXDCONT!'
            END IF

            IF( IASTAT .NE. 0 .OR. IBSTAT .NE. 0) THEN
               GO TO 9999
            END IF
            
            CALL MAXDCONT_LOOP

         END IF

      ELSE
C ---    FATAL error occurred during PRESETUP phase; initialize MODOPS string 
C        before call to HEADER to avoid runtime error for undefined string
         MODOPS_String = ''
            
      END IF

C     Determine Number of Errors/Messages by Message Type   ---   CALL TERRST
      CALL TERRST
      
C     Write Summary of Message Stats for Model Execution    ---   CALL SUMTBL
      CALL HEADER(IOUNIT)
      WRITE(IOUNIT,9114)
 9114 FORMAT(/1X,'*** Message Summary : AERMOD Model Execution ***'/)
 
      CALL SUMTBL(IOUNIT)
      
      IF (SUMMFILE) THEN
C        Write Summary of Messages to optional SUMMFILE      
         CALL HEADER(ISUMUNT)
         WRITE(ISUMUNT,9114)
         CALL SUMTBL(ISUMUNT)
      END IF

C     Skip to here if error occurs during allocation of arrays
 9999 CONTINUE

      IF (FATAL .OR. RUNERR) THEN
         IF (RUNERR) THEN
            WRITE(*,99112)
99112       FORMAT('+','Fatal Error Occurred During Runtime Phase!')
         END IF
         WRITE(IOUNIT,9115)
 9115    FORMAT(/4X,'***************************************',
     &          /4X,'*** AERMOD Finishes UN-successfully ***',
     &          /4X,'***************************************'/)
      ELSE
C ---    AERMOD finished without any "fatal" errors, but issue
C        warning to IOUNIT if MAXDCONT or EVENT processing 
C        inconsistencies occurred, before message that AERMOD 
C        finished successfully
         IF (L_MAXDCONT_OrigConc_Warning) THEN
C ---       MAXDCONT processing inconsistency warning
            WRITE(*,91161)
            WRITE(IOUNIT,91161)
91161       FORMAT(/4X,'NOTE: MAXDCONT Inconsistency Warning Issued!'/)
         ELSE IF (L_EVENT_OrigConc_Warning) THEN
C ---       EVENT processing inconsistency warning
            WRITE(*,91162)
            WRITE(IOUNIT,91162)
91162       FORMAT(/4X,'NOTE: EVENT Inconsistency Warning(s) Issued!')
         END IF
         WRITE(IOUNIT,9116)
9116     FORMAT(/4X,'************************************',
     &          /4X,'*** AERMOD Finishes Successfully ***',
     &          /4X,'************************************'/)
      END IF

      IF (ERRLST) THEN
C        OPEN and Write Out Permanent Error Message File    ---   CALL MSGWRT
         OPEN(UNIT=IERWRT,FILE=MSGFIL,STATUS='REPLACE',
     &        FORM='FORMATTED')
         CALL MSGWRT
         CLOSE(IERWRT)
      END IF

C     Close and Delete The Error Message And EVENT Temporary Files
      CLOSE(IERUNT,STATUS='DELETE')
      CLOSE(ITEVUT,STATUS='DELETE')

      STOP
      END

      SUBROUTINE HRLOOP
C***********************************************************************
C                 HRLOOP Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Controls Main Calculation Loop Through
C                 Hourly Meteorological Data
C
C        PROGRAMMER: Roger Brode, Jeff Wang
C
C        DATE:    March 2, 1992
C
C        MODIFIED:   Incorporated non-DFAULT/BETA ARM2 option for NO2
C                    Mark Podrez, RTP Environmental Associates, Inc.
C                    R. Brode, US EPA, OAQPS, AQMG, December 16, 2013
C
C        MODIFIED:  To include error handling for data mismatches between
C                   the hourly emissions and meteorological data files.
C                   Included error hanlding for end-of-file (EOF) for the
C                   meteorological data files occurring before the user-
C                   specified end-date (STARTEND keyword).  Also removed
C                   code related to "wet scimming" option, which is not
C                   supported in AERMOD.
C                   R.W. Brode, U.S. EPA/OAQPS/AQMG, 10/19/2009
C
C        MODIFIED:  To include the PVMRM and OLM options for
C                   modeling conversion of NOx to NO2.
C                   R. W. Brode, MACTEC (f/k/a PES), Inc., 07/27/04
C
C        MODIFIED:  To incorporate modifications to date processing
C                   for Y2K compliance, including use of date window
C                   variables (ISTRT_WIND and ISTRT_CENT) and calculation
C                   of 10-digit date variable (FULLDATE) with 4-digit
C                   year for date comparisons.
C                   Also modified to include SCIM option.
C                   R.W. Brode, PES, Inc., 5/12/99
C
C        MODIFIED:  To correct problems with the post-1997 PM10
C                   calculations involving leap years, and to
C                   add the year to the status message.
C                   R.W. Brode, PES, Inc. - 12/2/98
C
C        MODIFIED:  Changes to accommodate the post-1997 PM10
C                   calculations for average H4H 24-hour averages
C                   and ANNUAL averages.
C                   R.W. Brode, PES, Inc. - 8/14/98
C
C        MODIFIED:  Minor change to logic of IF block to correct
C                   potential problem with STARTEND keyword for
C                   non-sequential meteorological data sets.
C                   R.W. Brode, PES, Inc. - 4/22/96
C
C        MODIFIED:  To Include TOXXFILE Option - 9/29/92
C
C        INPUTS:  Source, Receptor and Setup Options
C
C        OUTPUTS: Update Hourly Results
C
C        CALLED FROM:   MAIN
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12
      CHARACTER BGReadErr*5, BGEndErr*5
      CHARACTER O3ReadErr*5, O3EndErr*5

      INTEGER :: I, ILSAVE

      DOUBLE PRECISION :: RDUM

C     Logical variable to identify READ errors in hourly BACKGROUND files
      LOGICAL :: L_BGReadErr, L_O3ReadErr

C     Variable Initializations
      MODNAM = 'HRLOOP'
      EOF = .FALSE.
      L_BGReadErr = .FALSE.
      L_O3ReadErr = .FALSE.
      IF (ALLOCATED(L_MorningTrans)) L_MorningTrans(:) = .FALSE.
      KURDAT   = 0
      FULLDATE = 0

C     Begin Hourly LOOP
      HOUR_LOOP: DO WHILE (FULLDATE.LT.IEDATE .AND. .NOT.EOF)
C        Retrieve One Hour of Meteorology                   ---   CALL METEXT
C        Call to METEXT also determines the sector IDs for 
C        BACKGRND and/or OZONE data, IBGSECT and IO3SECT
         CALL METEXT

C ---    Check for runtime error generated in call to METEXT;
C        Exit HOUR_LOOP if runtime error found
         IF (RUNERR) THEN
            EXIT HOUR_LOOP
         END IF
         
         IF (FULLDATE.GE.ISDATE .AND. FULLDATE.LE.IEDATE .AND. 
     &       ( (.NOT.L_LeapYear.AND.IPROC (JDAY).EQ.1) .OR.
     &              (L_LeapYear.AND.IPROCL(JDAY).EQ.1) ) .AND.
     &                                                  .NOT.EOF) THEN
C           Increment counter for total number of hours processed
            IF (.NOT.L_SkipMessages) NTOTHRS = NTOTHRS + 1
         ELSE IF (FULLDATE.LT.IEDATE .AND. IEDATE.LT.2147123124 .AND.
     &                                                         EOF) THEN
C ---       End of met data file(s) reached before user-specified End Date
C           Issue fatal error message
            CALL ERRHDL(PATH,MODNAM,'E','580','MET-DATA')
C           Exit hourly loop
            EXIT HOUR_LOOP
         ELSE IF (EOF .OR. FULLDATE .GT. IEDATE) THEN
C ---       End of File or data period has been reached; EXIT hour loop
            EXIT HOUR_LOOP
         END IF

C        Save ILINE as ILSAVE and Initialize ILINE
         ILSAVE = ILINE

         IF (HOURLY) THEN
C           Process Hourly Emissions from File
C           Begin Source Loop
            DO ISRC = 1, NUMSRC
               IF (QFLAG(ISRC) .EQ. 'HOURLY') THEN
C*                Increment IQLINE counter to reflect line number of HOUREMIS file
                  IQLINE = IQLINE + 1
C*                Retrieve Source Parameters for This Hour     ---   CALL HRQREAD
                  CALL HRQREAD(ISRC)
C*                Check for Date and Time Consistency with Met Data; 
C*                If Failed, Issue Fatal Error
                  IF (EOF) THEN
C*                   Write Error Message - EOF reached in hourly emission file
                     CALL ERRHDL(PATH,MODNAM,'E','580','HOUREMIS')
                     RUNERR = .TRUE.
                  ELSE IF (FULLDATE .NE. FULLHRQ) THEN
C*                   WRITE Error Message - Date mismatch
                     WRITE(DUMMY,'(I10.10)') FULLDATE
                     CALL ERRHDL(PATH,MODNAM,'E','455',DUMMY)
                     RUNERR = .TRUE.
                  END IF
C*                Extract source parameters to standard arrays, if not RUNERR
                  IF (.NOT. RUNERR) THEN
                     CALL HRQEXT(ISRC)
                  END IF

                  IF (.NOT.RSTINP .AND. L_MAXDCONT .AND. 
     &                                  FULLDATE.GE.ISDATE) THEN
C ---                Save hourly emissions for MAXDCONT option
                     AAQS(IHR_NDX,IYR_NDX,ISRC) = AQS(ISRC)
        
                     IF (SRCTYP(ISRC)(1:5) .EQ. 'POINT') THEN
                        AATS(IHR_NDX,IYR_NDX,ISRC) = ATS(ISRC)
                        AAVS(IHR_NDX,IYR_NDX,ISRC) = AVS(ISRC)
                     ELSE IF (SRCTYP(ISRC) .EQ. 'VOLUME' .AND. 
     &                                           L_HRLYSIG(ISRC)) THEN
                        AAHS(IHR_NDX,IYR_NDX,ISRC)    = AHS(ISRC)
                        AASYINI(IHR_NDX,IYR_NDX,ISRC) = ASYINI(ISRC)
                        AASZINI(IHR_NDX,IYR_NDX,ISRC) = ASZINI(ISRC)
                     ELSE IF (SRCTYP(ISRC)(1:4) .EQ. 'AREA' .AND. 
     &                                           L_HRLYSIG(ISRC)) THEN
                        AAHS(IHR_NDX,IYR_NDX,ISRC)    = AHS(ISRC)
                        AASZINI(IHR_NDX,IYR_NDX,ISRC) = ASZINI(ISRC)
                     ELSE IF (SRCTYP(ISRC) .EQ. 'LINE' .AND. 
     &                                           L_HRLYSIG(ISRC)) THEN
                        AAHS(IHR_NDX,IYR_NDX,ISRC)    = AHS(ISRC)
                        AASZINI(IHR_NDX,IYR_NDX,ISRC) = ASZINI(ISRC)
                     END IF
                  END IF

               END IF
            END DO
C*          End Source Loop
         END IF

C ---    Check for runtime error generated in call to HRQREAD;
C        Exit HOUR_LOOP if runtime error found
         IF (RUNERR) EXIT HOUR_LOOP
         
C*----
C        Save ILINE as ILSAVE and Initialize ILINE
         ILSAVE = ILINE

         IF (L_BACKGRND) THEN
C ---       Process Background Concentration inputs, starting with Hourly BG

            IF (L_BGHourly) THEN
C*             Increment IBLINE counter to reflect line number of hourly BACKGRND file
               IBLINE = IBLINE + 1

C*             Retrieve hourly background concentrations      ---   CALL BGEXT
               CALL BGEXT(L_BGReadErr,BGReadErr,BGEndErr)

C*             Check for issues reading hourly BG file
               IF (EOF) THEN
C*                Write Error Message - EOF reached in hourly background file
                  WRITE(DUMMY,'(''BGFILE '',A5)') BGEndErr
                  CALL ERRHDL(PATH,MODNAM,'E','580',DUMMY)
                  RUNERR = .TRUE.
               ELSE IF (L_BGReadErr) THEN
C*                Write Error Message - READ error in hourly background file
                  WRITE(DUMMY,'(''BGFILE '',A5)') BGReadErr
                  CALL ERRHDL(PATH,MODNAM,'E','510',DUMMY)
                  RUNERR = .TRUE.
               END IF

            ELSE

C ---          Check for temporally-varying background to substitute for missing hours
               IF (IBGSECT .GT. 0) THEN
C                 Valid IBGSECT value
                  IF (L_BGValues(IBGSECT)) THEN
                     CALL BGVAL(IBGSECT,BGCONC)
                  ELSE
                     BGCONC = 0.0D0
                  END IF
               ELSE
C ---             IBGSECT is missing (calm or missing hour); set BGCONC = 0.0
C                 Note: this should result in a fatal error (452)
                  BGCONC = 0.0D0
               END IF
            END IF

         END IF
            
C ---    Check for runtime error generated in call to BGEXT;
C        Exit HOUR_LOOP if runtime error found
         IF (RUNERR) EXIT HOUR_LOOP
         
         IF (L_BACKGRND .AND. .NOT.RSTINP .AND. L_MAXDCONT .AND. 
     &                         FULLDATE.GE.ISDATE) THEN
C ---       Save hourly background concentration for MAXDCONT option
            ABGCONC(IHR_NDX,IYR_NDX) = BGCONC
         END IF
            
C*----
C        Retrive ILINE From ILSAVE
         ILINE = ILSAVE

         IF (PVMRM .OR. PVMRM2 .OR. OLM) THEN
C-----      Read Ozone Data File if available
            IF (L_O3Hourly) THEN
C*             Increment IOLINE counter to reflect line number of Hourly O3 file
               IOLINE = IOLINE + 1

C ---          Extract O3 value from hourly data file; O3EXT also reads a record 
C              for O3FILEs available for other sectors to keep files synchronized, 
C              so premature EOF for one file should result in date mismatch error
               IF (.NOT. EOF) THEN
                  CALL O3EXT(L_O3ReadErr,O3ReadErr,O3EndErr)
               END IF

               IF (EOF) THEN
C*                Write Error Message - EOF reached in hourly O3 file
                  WRITE(DUMMY,'(''O3FILE '',A5)') O3EndErr
                  CALL ERRHDL(PATH,MODNAM,'E','580',DUMMY)
                  RUNERR = .TRUE.
               ELSE IF (L_O3ReadErr) THEN
C*                Write Error Message - EOF reached in hourly O3 file
                  WRITE(DUMMY,'(''O3FILE '',A5)') O3ReadErr
                  CALL ERRHDL(PATH,MODNAM,'E','510',DUMMY)
                  RUNERR = .TRUE.
               END IF

            ELSE IF (IO3SECT .GT. 0) THEN
               IF (L_O3VALUES(IO3SECT)) THEN
C ---             Use ozone concentration from O3VALUES keyword
                  CALL OZONVALS(IO3SECT,O3CONC)
               ELSE IF (L_O3VAL(IO3SECT)) THEN
C ---             Use single "background" O3 value from OZONEVAL keyword
                  O3CONC = O3BACK(IO3SECT)
               ELSE
C ---             Set O3CONC to 0.0 for full conversion (subject to 
C                 equilibrium ratio)
                  O3CONC = 0.0D0
               END IF
            ELSE
C ----         IO3SECT is 0 due to calm/missing hour; set O3CONC to 0.0
               O3CONC = 0.0D0
            END IF

            IF (.NOT.RSTINP .AND. L_MAXDCONT .AND. 
     &                            FULLDATE.GE.ISDATE) THEN
C ---          Save hourly ozone concentration for MAXDCONT option
               AO3CONC(IHR_NDX,IYR_NDX) = O3CONC
            END IF
         END IF

C ---    Check for runtime error generated in call to O3EXT;
C        Exit HOUR_LOOP if runtime error found
         IF (RUNERR) EXIT HOUR_LOOP
         
C*----
C        Retrive ILINE From ILSAVE
         ILINE = ILSAVE

C*       Check for IHOUR = 1 and Write Update to the Screen For PC Version
         IF ((IHOUR.EQ.1 .OR. ILINE.EQ.1) .AND. .NOT.NOCHKD) THEN
C*          Write Out Update to the Screen by Julian Day
            WRITE(*,909) JDAY, IYR
 909        FORMAT('+','Now Processing Data For Day No. ',I4,' of ',I4)
         ELSE IF (NOCHKD) THEN
C*          Write Out Update to the Screen by Hour
            WRITE(*,910) KURDAT
 910        FORMAT('+','Now Processing Data For     ',I8.8)
         END IF
C*----
C*#
         IF (SCIM .AND. .NOT.EOF) THEN
            SCIMHR = .FALSE.

C           User has specified SCIM option.  Check for whether current
C           hour is to be sampled, and whether to write sampled met
C           data to output file.

C           Keep track of total no. of hours.
C           Also, keep track of dry & wet, and calm & missing hours
C           Note:  Under SCIM option, IANHRS/IANCLM/IANMSG (see below) pertain
C                  to no. of hours sampled.
            NSKIPTOT = NSKIPTOT + 1

            IF( ILINE .LE. 24 .AND. IHOUR .EQ. NREGSTART )THEN
C              Current hour is to be sampled - first SCIM'd hour.
               IFIRSTHR = ILINE
               SCIMHR   = .TRUE.
            ELSE IF( ILINE .GT. NREGSTART .AND.
     &               MOD( ILINE-IFIRSTHR, NREGINT ) .EQ. 0 )THEN
C              Current hour is to be sampled - SCIM'd hour
               SCIMHR   = .TRUE.
            ELSE
C              Current hour is NOT to be sampled. Check for end of year first.
               CALL CHK_ENDYR
               CYCLE HOUR_LOOP
            END IF

            IF (SCIMOUT) THEN
C              Write sampled meteorology to SCIM'd met data file
               CALL METSUM
            END IF
         END IF

         IF (FULLDATE.GE.ISDATE .AND. FULLDATE.LE.IEDATE .AND.
     &       ( (.NOT.L_LeapYear.AND.IPROC (JDAY).EQ.1) .OR.
     &              (L_LeapYear.AND.IPROCL(JDAY).EQ.1) ) .AND.
     &                .NOT.EOF .AND. .NOT.RUNERR) THEN

C ---       Check for calm winds or missing met data, for which model
C           calculations cannot be made; increment counters for number
C           of hours, but do not include background concentrations, if
C           specified through the BACKGRND keyword.
            IF (CLMHR .AND. CLMPRO) THEN
C              Check for Calm Hr & Processing and Increment Counters
               DO IAVE = 1, NUMAVE
                  NUMHRS(IAVE) = NUMHRS(IAVE) + 1
                  NUMCLM(IAVE) = NUMCLM(IAVE) + 1
               END DO
               IF (PERIOD .OR. ANNUAL) THEN
                  IF (.NOT.SCIM .OR. (SCIM.AND.SCIMHR)) THEN
                     IANHRS = IANHRS + 1
                     IANCLM = IANCLM + 1
                  END IF
               END IF
               IF (SEASONHR) THEN
                  NSEAHR(ISEAS,IHOUR) = NSEAHR(ISEAS,IHOUR) + 1
                  NSEACM(ISEAS,IHOUR) = NSEACM(ISEAS,IHOUR) + 1
               END IF
            ELSE IF (MSGHR .AND. MSGPRO) THEN
C              Check for Missing Hour & Processing and Increment Counters
               DO IAVE = 1, NUMAVE
                  NUMHRS(IAVE) = NUMHRS(IAVE) + 1
                  NUMMSG(IAVE) = NUMMSG(IAVE) + 1
               END DO
               IF (PERIOD .OR. ANNUAL) THEN
                  IF (.NOT.SCIM .OR. (SCIM.AND.SCIMHR)) THEN
                     IANHRS = IANHRS + 1
                     IANMSG = IANMSG + 1
                  END IF
               END IF
               IF (SEASONHR) THEN
                  NSEAHR(ISEAS,IHOUR) = NSEAHR(ISEAS,IHOUR) + 1
                  NSEACM(ISEAS,IHOUR) = NSEACM(ISEAS,IHOUR) + 1
               END IF
            ELSE IF (ZI .LE. 0.0D0) THEN
C              Write Out The Informational Message & Increment Counters
               IF (.NOT. L_SkipMessages) THEN
                  WRITE(DUMMY,'(I8.8)') KURDAT
                  CALL ERRHDL(PATH,MODNAM,'I','470',DUMMY)
               END IF
               DO IAVE = 1, NUMAVE
                  NUMHRS(IAVE) = NUMHRS(IAVE) + 1
               END DO
               IF (PERIOD .OR. ANNUAL) THEN
                  IF (.NOT.SCIM .OR. (SCIM.AND.SCIMHR)) THEN
                     IANHRS = IANHRS + 1
                  END IF
               END IF
               IF (SEASONHR) THEN
                  NSEAHR(ISEAS,IHOUR) = NSEAHR(ISEAS,IHOUR) + 1
               END IF
            ELSE
C              Set CALCS Flag, Increment Counters & Calculate HRVAL
               CALCS = .TRUE.
               DO IAVE = 1, NUMAVE
                  NUMHRS(IAVE) = NUMHRS(IAVE) + 1
               END DO
               IF (PERIOD .OR. ANNUAL) THEN
                  IF (.NOT.SCIM .OR. (SCIM.AND.SCIMHR)) THEN
                     IANHRS = IANHRS + 1
                  END IF
               END IF
               IF (SEASONHR) THEN
                  NSEAHR(ISEAS,IHOUR) = NSEAHR(ISEAS,IHOUR) + 1
               END IF

C              Time/Date Marker for DEBUG Output
               IF (DEBUG) THEN
                  WRITE(DBGUNT,*)
                  WRITE(DBGUNT,*) '--------------------------------',
     &                            '--------------------'
                  WRITE(DBGUNT,*) '---  JDAY, IHOUR =  ',JDAY,IHOUR
                  WRITE(DBGUNT,*) '--------------------------------',
     &                            '--------------------'
               END IF

C ---          Calculate CONC or DEPOS Values               ---   CALL CALC
               CALL CALC
            END IF

            IF (.NOT.CLMHR .AND. .NOT.MSGHR) THEN
C ---          Non-calm, non-missing hour; apply NO2 options as appropriate

               IF ((PVMRM .OR. PVMRM2) .AND. .NOT.PSDCREDIT) THEN
C ---             Process Hourly Values for PVMRM Option
                  CALL PVMRM_CALC('ALLSRCS')
                  
               ELSE IF ((PVMRM .OR. PVMRM2) .AND. PSDCREDIT) THEN
C ---             Process Hourly Values for PVMRM Option and PSD credits
C ---             Need to process two separate sets of sources - the
C                 increment consumption sources ('NAAQSRC') and the 
C                 increment expanding sources ('ALLBASE')
                  CALL PVMRM_CALC('NAAQSRC')
                  CALL PVMRM_CALC('ALLBASE')

               ELSE IF (OLM) THEN
C ---             Process Hourly Values for OLM Option
                  CALL OLM_CALC

               ELSE IF (ARM) THEN
C ---             Process Hourly Values for ARM Option
                  CALL ARM_CALC

               ELSE IF (ARM2) THEN
C ---             Process Hourly Values for ARM2 Option
                  CALL ARM2_CALC

               END IF

            END IF

C           Begin Averaging Period LOOP
            DO IAVE = 1, NUMAVE
C              Check for End of Averaging Period
               IF (MOD(IHOUR,KAVE(IAVE)).EQ.0 .OR.
     &            (KAVE(IAVE).EQ.720 .AND. ENDMON)) THEN
                  IF (CONC) THEN
C                    Calculate Applicable Averages          ---   CALL AVER
                     CALL AVER
                  END IF
C                 Update High Value Arrays                  ---   CALL HIVALS
                  CALL HIVALS

                  IF( (NO2AVE .OR. SO2AVE) .AND. KAVE(IAVE).EQ.1 )THEN
C ---                Loop through SRCGRPs again to get max daily 1hr cumulative value
                     DO IGRP = 1, NUMGRP
                        DO IREC = 1, NUMREC
                           IF (AVEVAL(IREC,IGRP,IAVE,1) .GT. 
     &                                          MXDVAL(IREC,IGRP)) THEN
                            MXDVAL(IREC,IGRP) = AVEVAL(IREC,IGRP,IAVE,1)
                            IMXDHR(IREC,IGRP) = IHOUR
                           END IF
                        END DO
                     END DO
                  END IF
                  
                  IF( PM25AVE .AND. MOD(IHOUR,24).EQ.0 .AND. 
     &                                           KAVE(IAVE).EQ.24 )THEN
C ---                Loop through source groups again to get max daily 24-hr cumulative value
                     DO IGRP = 1, NUMGRP
                        DO IREC = 1, NUMREC
                           IF (AVEVAL(IREC,IGRP,IAVE,1) .GT. 
     &                                          MXDVAL(IREC,IGRP)) THEN
                            MXDVAL(IREC,IGRP) = AVEVAL(IREC,IGRP,IAVE,1)
                            IMXDHR(IREC,IGRP) = IHOUR
                           END IF
                        END DO
                     END DO
                  END IF
            
                  IF (DAYTAB .AND. IDYTAB(IAVE).EQ.1) THEN
                     DO ITYP = 1, NUMTYP
C                       Print Out Daily Value Tables        ---   CALL PRTDAY
                        CALL PRTDAY
                     END DO
                  END IF
                  IF (MXFILE) THEN
C                    Write Max Values (>Thresh) to File     ---   CALL MAXFIL
                     CALL MAXFIL
                  END IF
                  IF (PPFILE) THEN
C                    Write Values to Postprocessor File     ---   CALL POSTFL
                     CALL POSTFL
                  END IF
                  IF (TXFILE) THEN
C                    Write Values to TOXXFILE File (9/29/92) ---  CALL TOXXFL
                     CALL TOXXFL
                  END IF
C                 Flush Block Average Values in AVEVAL Array for This IAVE
                  AVEVAL(1:NUMREC,1:NUMGRP,IAVE,1:NUMTYP) = 0.0D0
               END IF
            END DO
C           End Averaging Period LOOP

C ---       Check for PM25AVE, NO2AVE or SO2AVE to update daily
C           maximum value arrays; also output to MAXDAILY file, 
C           if requested
            IF (PM25AVE .OR. NO2AVE .OR. SO2AVE) THEN
               IF (MOD(IHOUR,24).EQ.0) THEN
C ---             End of day reached, call MXDLYFL
                  CALL MXDLYFL
               END IF
            END IF
            
            IF (RSTSAV .AND. IHOUR.EQ.24) THEN
               NDAYS = NDAYS + 1
               IF (NDAYS .EQ. INCRST) THEN
C                 Save Results to File for Later Re-start   ---   CALL RSDUMP
                  CALL RSDUMP
                  NDAYS = 0
               END IF
            END IF

C           Flush HRVAL Arrays (1:NUMTYP)
            HRVAL(:)   = 0.0D0
            AERVAL(:)  = 0.0D0
            PRMVAL(:)  = 0.0D0
            IF (ALLOCATED(BACKAVE)) BACKAVE(:)  = 0.0D0
            IF (ALLOCATED(BACKHR))  BACKHR(:,:) = 0.0D0

            IF (PVMRM .OR. PVMRM2 .OR. OLM .OR. ARM2 .OR. ARM) THEN
C              Flush CHI(NUMREC,NUMSRC,NUMTYP) Array 
               CHI(:,:,:) = 0.0D0
               IF (PSDCREDIT) THEN
C                 Flush ABVAL(NUMREC,NUMTYP) and BCVAL(NUMREC,NUMTYP) Arrays
                  ABVAL(:,:) = 0.0D0
                  BCVAL(:,:) = 0.0D0
               END IF
            END IF

         END IF

C        Check for end of year of data for PM25, NO2, SO2, or MULTYR processing;
C        but skip if NOCHKD option or WARNCHKD option is used (this also includes
C        SCREEN option since SCREEN ==> NOCHKD)
         IF (FULLDATE.GT.ISDATE .AND. .NOT.EOF .AND. .NOT.NOCHKD .AND.
     &                                           .NOT.L_WARNCHKD .AND.
     &                           (PM25AVE .OR. NO2AVE .OR. SO2AVE .OR.
     &                                         ANNUAL .OR. MULTYR)) THEN

            CALL CHK_ENDYR

         ELSEIF(FULLDATE.EQ.ISDATE .AND..NOT.EOF .AND..NOT.NOCHKD .AND.
     &                                            .NOT.L_WARNCHKD .AND.
     &                           (PM25AVE .OR. NO2AVE .OR. SO2AVE .OR.
     &                                         ANNUAL .OR. MULTYR)) THEN

            NREMAIN = NREMAIN + 1

         END IF

C        Reset CALCS and ENDMON Flags
         CALCS  = .FALSE.
         ENDMON = .FALSE.

C        Save precipitation rates for two previous hours
         prec2 = prec1
         prec1 = Prate

      END DO HOUR_LOOP
C     End Hourly LOOP

C     Check for TOXXFILE Option, Fill Buffer and Dump to File - 9/29/92
      IF (TXFILE) THEN
         IDUM = 0
         RDUM = 0.0D0
         DO IAVE = 1, NUMAVE
            IF (ITOXFL(IAVE) .EQ. 1) THEN
C              Fill Rest of Buffer With Zeroes and Write to TOXXFILE
               DO I = IPAIR+1, NPAIR
                  IDCONC(IAVE,I) = IDUM
                  TXCONC(IAVE,I) = RDUM
               END DO
               WRITE(ITXUNT(IAVE)) (IDCONC(IAVE,I),I=1,NPAIR)
               WRITE(ITXUNT(IAVE)) (TXCONC(IAVE,I),I=1,NPAIR)
               CLOSE(ITXUNT(IAVE))
            END IF
         END DO
      END IF

C     Write Out Update to the Screen for PC Version
      WRITE(*,919)
 919  FORMAT('+','Now Processing Output Options               ')

      RETURN
      END

      SUBROUTINE JULIAN(INYR,INMN,INDY,JDY)
C***********************************************************************
C                 JULIAN Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:    CONVERT YR/MN/DY DATE TO JULIAN DAY (1-366),
C                    INCLUDES TEST FOR 100 AND 400 YEAR CORRECTIONS TO
C                    HANDLE 4 DIGIT YEARS BEYOND 2099 AND BEFORE 1901
C                    (WILL WORK WITH 2 DIGIT YR FOR PERIOD 1901-2099)
C
C        PROGRAMMER: Roger Brode, Jeff Wang
C
C        DATE:    March 2, 1992
C
C        INPUTS:     YEAR,  INYR (2 OR 4 DIGIT)
C                    MONTH, INMN
C                    DAY,   INDY
C
C        OUTPUT:     JULIAN DAY,  JDY (1-366)
C
C        CALLED FROM:   DAYRNG
C
C        ERROR HANDLING:   Checks for Invalid Month or Day
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: NDAY(12), IDYMAX(12)
      INTEGER :: INYR, INMN, INDY, JDY

C     Variable Initializations
      DATA NDAY/0,31,59,90,120,151,181,212,243,273,304,334/
      DATA IDYMAX/31,29,31,30,31,30,31,31,30,31,30,31/
      MODNAM = 'JULIAN'
      JDY = 0

C     Check for Invalid Month or Day
      IF (INMN.LT.1 .OR. INMN.GT.12) THEN
C        WRITE Error Message    ! Invalid Month
         WRITE(DUMMY,'(''MONTH = '',I2)') INMN
         CALL ERRHDL(PATH,MODNAM,'E','203',DUMMY)
         RUNERR = .TRUE.
         GO TO 999
      ELSE IF (INDY .GT. IDYMAX(INMN)) THEN
C        WRITE Error Message    ! Invalid Day
         WRITE(DUMMY,'(''DAY='',I2,'' MO='',I2)') INDY,INMN
         CALL ERRHDL(PATH,MODNAM,'E','203',DUMMY)
         RUNERR = .TRUE.
         GO TO 999
      END IF

C     Determine JULIAN Day Number; For Non-Leap Year First
      IF ((MOD(INYR,4) .NE. 0) .OR.
     &    (MOD(INYR,100) .EQ. 0 .AND. MOD(INYR,400) .NE. 0)) THEN
C        Not a Leap Year
         IF (INMN.NE.2 .OR. (INMN.EQ.2 .AND. INDY.LE.28)) THEN
            JDY = INDY + NDAY(INMN)
         ELSE
C           WRITE Error Message    ! Invalid Date; 2/29 in a Non-Leap Year
            WRITE(DUMMY,'("YR= ",I4)') INYR
            CALL ERRHDL(PATH,MODNAM,'E','370',DUMMY)
            JDY = 60
            RUNERR = .TRUE.
         END IF
      ELSE
C        Leap Year
         JDY = INDY + NDAY(INMN)
         IF (INMN .GT. 2)  JDY = JDY + 1
      END IF

 999  CONTINUE

      RETURN
      END

      SUBROUTINE GREGOR(INYR,INMN,JDY,IDY)
C***********************************************************************
C                 GREGOR Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:    CONVERT JULIAN DAY (1-366) TO DAY OF MONTH,
C                    INCLUDES TEST FOR 100 AND 400 YEAR CORRECTIONS TO
C                    HANDLE 4 DIGIT YEARS BEYOND 2099 AND BEFORE 1901
C                    (WILL WORK WITH 2 DIGIT YR FOR PERIOD 1901-2099)
C
C        PROGRAMMER: Roger Brode, Jeff Wang
C
C        DATE:    March 2, 1992
C
C        INPUTS:     YEAR,       INYR (2 OR 4 DIGIT)
C                    MONTH,      INMN
C                    JULIAN DAY, JDY (1-366)
C
C        OUTPUT:     DAY OF MONTH, IDY
C
C        CALLED FROM:   METEXT
C
C        ERROR HANDLING:   Checks for Invalid Month or Day
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: NDAY(12)
      INTEGER :: INYR, INMN, IDY, JDY

C     Variable Initializations
      DATA NDAY/0,31,59,90,120,151,181,212,243,273,304,334/
      MODNAM = 'GREGOR'

C     Check for Invalid Month or Julian Day
      IF (INMN.LT.1 .OR. INMN.GT.12) THEN
C        WRITE Error Message    ! Invalid Month
         CALL ERRHDL(PATH,MODNAM,'E','203','MONTH')
         GO TO 999
      ELSE IF (JDY.LT.1 .OR. JDY.GT.366) THEN
C        WRITE Error Message    ! Invalid Julian Day
         CALL ERRHDL(PATH,MODNAM,'E','203','Juli Day')
         GO TO 999
      END IF

C     Determine Day-of-Month Number; For Non-Leap Year First
      IF ((MOD(INYR,4) .NE. 0) .OR.
     &    (MOD(INYR,100).EQ.0 .AND. MOD(INYR,400).NE.0)) THEN
C        Not a Leap Year
         IDY = JDY - NDAY(INMN)
      ELSE
C        Leap Year
         IDY = JDY - NDAY(INMN)
         IF (INMN .GT. 2)  IDY = IDY - 1
      END IF

 999  CONTINUE

      RETURN
      END

      SUBROUTINE HRQREAD (IS)
C***********************************************************************
C*                  HRQREAD Module of AERMOD
C* 
C*         PURPOSE: To Assign Hourly Source Parameters
C* 
C*         PROGRAMMER:  Jayant Hardikar, Roger Brode
C* 
C*         DATE:    September 15, 1993
C* 
C*         INPUTS:  Current Source Number Being Processed
C* 
C*         OUTPUTS: Source Arrays
C*
C*         Revision History:
C*
C*         MODIFIED:  Added processing of hourly emissions file for 
C*                    a buoyant line source.
C*                    Amec Foster Wheeler, 6/30/2015
C*
C*         MODIFIED:  Included check on length of FIELD(7) before assigning
C*                    to the HRSOID variable to avoid runtime error if
C*                    SRCIDs in HOUREMIS file exceed 12 character limit.
C*                    R.W. Brode, U.S. EPA/OAQPS/AQMG, 03/19/2014
C*
C*         MODIFIED:  Check for use of 4-digit year in HOUREMIS file, and
C*                    adjust if needed for comparison to KURDAT from the
C*                    met data file.
C*                    Incorporated options to specify hourly-varying
C*                    release heights and initial dispersion coefficients
C*                    for VOLUME and AREA sources.
C*                    R.W. Brode, U.S. EPA/OAQPS/AQMG, 10/19/2009
C*         
C*         MODIFIED:  Corrected processing of missing parameters for 
C*                    point sources to assign all parameters to 0.0 if
C*                    any of the parameters are missing, in conformance
C*                    with Section 3.3.9 of the AERMOD User's Guide.
C*                    R.W. Brode, U.S. EPA/OAQPS/AQMG, 01/24/2007
C*         
C*         MODIFIED:  REMOVED THE 'POINT' SOURCE CONDITION, SO IT APPLIES 
C*                    TO ALL SOURCE TYPES, EXCEPT SAVING THE TEMP & VEL
C* 
C*         CALLED FROM:  HRLOOP
C************************************************************************
C*
C*    Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: I, IS
      INTEGER :: IHYEAR, IHMON, IHDAY, IHHOUR, IHYEAR2
      INTEGER :: ILSAVE
      CHARACTER (LEN=20) :: RDFRM

      CHARACTER (LEN=12) :: HRSOID
      
C*    Variable Initializations
      MODNAM = 'HRQREAD'
      
C*    Assign IQLINE counter to ILINE for passing to ERRHDL if needed, save as ILSAVE first
      ILSAVE = ILINE
      ILINE  = IQLINE

C*    READ Record to Buffers, A'num' and 'num'A1, where num=ISTRG
C*    Length of ISTRG is Set in PARAMETER Statement in MAIN1
C     Setup READ format and ECHO format for runstream record,
C     based on the ISTRG PARAMETER (set in MAIN1)
      WRITE(RDFRM,9100) ISTRG, ISTRG
 9100 FORMAT('(A',I4.4,',T1,',I4.4,'A1)')
      READ (IHREMI,RDFRM,END=888,ERR=99) RUNST1, (RUNST(I), I=1, ISTRG)
C*
C*    Convert Lower Case to Upper Case Letters              ---   CALL LWRUPR
      CALL LWRUPR
C*      
C*    Define Fields on Card                                 ---   CALL DEFINE
      CALL DEFINE
C*
C*    Get the Contents of the Fields                        ---   CALL GETFLD
      CALL GETFLD
C*
C*    Check for number of fields - error if less than 7.
      IF (IFC .LT. 7) THEN
         WRITE(DUMMY,'(I8)') KURDAT
         CALL ERRHDL(PATH,MODNAM,'E','384',DUMMY)
         RUNERR = .TRUE.
         GO TO 999
      END IF
C*         
C*    Assign the Fields to Local Varables and Check The Numerical Field
C*
      CALL STONUM(FIELD(3), ILEN_FLD, FNUM, IMIT)
      IHYEAR = NINT(FNUM)
      IF (IMIT .NE. 1) THEN
         CALL ERRHDL(PATH,MODNAM,'E','208','HOUREMIS')
         RUNERR = .TRUE.
         GO TO 999
      END IF

      CALL STONUM(FIELD(4), ILEN_FLD, FNUM, IMIT)
      IHMON = NINT(FNUM)
      IF (IMIT .NE. 1) THEN
         CALL ERRHDL(PATH,MODNAM,'E','208','HOUREMIS')
         RUNERR = .TRUE.
         GO TO 999
      END IF

      CALL STONUM(FIELD(5), ILEN_FLD, FNUM, IMIT)
      IHDAY = NINT(FNUM)
      IF (IMIT .NE. 1) THEN
         CALL ERRHDL(PATH,MODNAM,'E','208','HOUREMIS')
         RUNERR = .TRUE.
         GO TO 999
      END IF

      CALL STONUM(FIELD(6), ILEN_FLD, FNUM, IMIT)
      IHHOUR = NINT(FNUM)
      IF (IMIT .NE. 1) THEN
         CALL ERRHDL(PATH,MODNAM,'E','208','HOUREMIS')
         RUNERR = .TRUE.
         GO TO 999
      END IF

C --- Check for use of 2-digit year in HOUREMIS file, adjust to 4-digit
C     year for comparison with FULLDATE based on met data file
      IF (IHYEAR .LE. 99) THEN
         IHYEAR2 = IHYEAR
         IF (IHYEAR2 .GE. ISTRT_WIND .AND. 
     &                        IHYEAR2 .LE. 99) THEN
            IHYEAR = ISTRT_CENT*100 + IHYEAR2
         ELSE IF (IHYEAR2 .LT. ISTRT_WIND) THEN
            IHYEAR = (ISTRT_CENT+1)*100 + IHYEAR2
         END IF
      END IF

C --- Calculate current date (YYYYMMDDHH) from HOUREMIS file record, FULLHRQ
      FULLHRQ = IHYEAR*1000000 + IHMON*10000 + IHDAY*100 + IHHOUR

C --- Assign source ID but check for field length > 12 first
      IF( LEN_TRIM(FIELD(7)) .LE. 12 ) THEN
         HRSOID = FIELD(7)
      ELSE
         HRSOID = FIELD(7)(1:12)
      END IF

C*    Check for Source ID Consistency ; If Failed Issue Error
      IF ( HRSOID .NE. SRCID(IS) ) THEN
         WRITE(DUMMY,'(A12)') SRCID(IS)
         CALL ERRHDL(PATH,MODNAM,'E','342',SRCID(IS))
         RUNERR = .TRUE.
         GO TO 999
      END IF

      IF (IFC .EQ. 7) THEN
C*       All parameters missing for this hour/source - WRITE Warning Message
C*       Assign zeros to all parameters
         IF (.NOT. L_SkipMessages) THEN
            WRITE(DUMMY,'(I10.10)') FULLHRQ
            CALL ERRHDL(PATH,MODNAM,'W','344',DUMMY)
         END IF
         HRQS = 0.0D0
         HRTS = 0.0D0
         HRVS = 0.0D0
         HRHS = 0.0D0
         HRSY = 0.0D0
         HRSZ = 0.0D0         
      
      ELSE IF (SRCTYP(IS)(1:5) .EQ. 'POINT' .AND. IFC.EQ.10) THEN
C*       Assign emission rate, exit temperature and exit velocity
C*       for POINT sources

         CALL STODBL(FIELD(8), ILEN_FLD, HRQS, IMIT)
         IF (IMIT .NE. 1) THEN
            CALL ERRHDL(PATH,MODNAM,'E','208','HOUREMIS')
            RUNERR = .TRUE.
            GO TO 999
C ---    Check for large negative values, could be a missing indicator
         ELSE IF ( HRQS .LE. -90.0D0 ) THEN
C*          Assume emissions are missing; assign value of 0.0 and issue Warning;
C*          Note that AERMOD User's Guide (p. 3-38, 3-39) implies that blanks 
C*          should be used for missing values.
            HRQS = 0.0D0  
            WRITE(DUMMY,'(I10.10)') FULLHRQ
            CALL ERRHDL(PATH,MODNAM,'W','341',DUMMY)
         END IF
        
         CALL STODBL(FIELD(9), ILEN_FLD, HRTS, IMIT)
         IF (IMIT .NE. 1) THEN
            CALL ERRHDL(PATH,MODNAM,'E','208','HOUREMIS')
            RUNERR = .TRUE.
            GO TO 999
         END IF

         CALL STODBL(FIELD(10), ILEN_FLD, HRVS, IMIT)
         IF (IMIT .NE. 1) THEN
            CALL ERRHDL(PATH,MODNAM,'E','208','HOUREMIS')
            RUNERR = .TRUE.
            GO TO 999
         END IF

      ELSE IF (SRCTYP(IS) .EQ. 'VOLUME' .AND. IFC.EQ.11) THEN
C*       Assign emission rate, release height and initial sigmas
C*       for VOLUME source.
C*       Assign logical variable indicating hourly sigmas, L_HRLYSIG
         IF (ILSAVE .EQ. 1) THEN
            L_HRLYSIG(IS) = .TRUE.
         ELSE IF (ILSAVE .GT. 1 .AND. .NOT. L_HRLYSIG(IS)) THEN
C*          This volume source should not include hourly sigmas;
C*          issue error message
            WRITE(DUMMY,'(I10.10)') FULLHRQ
            CALL ERRHDL(PATH,MODNAM,'E','345',DUMMY)
            HRQS = 0.0D0
            RUNERR = .TRUE.
            GO TO 999
         END IF

         CALL STODBL(FIELD(8), ILEN_FLD, HRQS, IMIT)
         IF (IMIT .NE. 1) THEN
            CALL ERRHDL(PATH,MODNAM,'E','208','HOUREMIS')
            RUNERR = .TRUE.
            GO TO 999
C ---    Check for large negative values, could be a missing indicator
         ELSE IF ( HRQS .LE. -90.0D0 ) THEN
C*          Assume emissions are missing; assign value of 0.0 and issue Warning;
C*          Note that AERMOD User's Guide (p. 3-38, 3-39) implies that blanks 
C*          should be used for missing values.
            HRQS = 0.0D0  
            WRITE(DUMMY,'(I10.10)') FULLHRQ
            CALL ERRHDL(PATH,MODNAM,'W','341',DUMMY)
         END IF
        
         CALL STODBL(FIELD(9), ILEN_FLD, HRHS, IMIT)
         IF (IMIT .NE. 1) THEN
            CALL ERRHDL(PATH,MODNAM,'E','208','HOUREMIS')
            RUNERR = .TRUE.
            GO TO 999
         END IF

         CALL STODBL(FIELD(10), ILEN_FLD, HRSY, IMIT)
         IF (IMIT .NE. 1) THEN
            CALL ERRHDL(PATH,MODNAM,'E','208','HOUREMIS')
            RUNERR = .TRUE.
            GO TO 999
         END IF

         CALL STODBL(FIELD(11), ILEN_FLD, HRSZ, IMIT)
         IF (IMIT .NE. 1) THEN
            CALL ERRHDL(PATH,MODNAM,'E','208','HOUREMIS')
            RUNERR = .TRUE.
            GO TO 999
         END IF

      ELSE IF (SRCTYP(IS) .EQ. 'VOLUME' .AND. IFC.EQ.8) THEN
C*       Assign emission rate for volume sources
C*       Check logical variable indicating hourly sigmas, L_HRLYSIG
         IF (L_HRLYSIG(IS)) THEN
C*          WRITE Error Message; Hourly Sigmas must be used for all hours         
            WRITE(DUMMY,'(I10.10)') FULLHRQ
            CALL ERRHDL(PATH,MODNAM,'E','345',DUMMY)
            HRQS = 0.0D0
            RUNERR = .TRUE.
            GO TO 999
         END IF

         CALL STODBL(FIELD(8), ILEN_FLD, HRQS, IMIT)
         IF (IMIT .NE. 1) THEN
            CALL ERRHDL(PATH,MODNAM,'E','208','HOUREMIS')
            RUNERR = .TRUE.
            GO TO 999
C ---    Check for large negative values, could be a missing indicator
         ELSE IF ( HRQS .LE. -90.0D0 ) THEN
C*          Assume emissions are missing; assign value of 0.0 and issue Warning;
C*          Note that AERMOD User's Guide (p. 3-38, 3-39) implies that blanks 
C*          should be used for missing values.
            HRQS = 0.0D0  
            WRITE(DUMMY,'(I10.10)') FULLHRQ
            CALL ERRHDL(PATH,MODNAM,'W','341',DUMMY)
         END IF

      ELSE IF ((SRCTYP(IS)(1:4).EQ.'AREA' .OR. SRCTYP(IS).EQ.'LINE')
     &                                    .AND. IFC.EQ.10) THEN
C*       Assign emission rate for AREA and LINE sources
C*       Assign logical variable indicating hourly sigmas, L_HRLYSIG
         IF (ILSAVE .EQ. 1) THEN
            L_HRLYSIG(IS) = .TRUE.
         ELSE IF (ILSAVE .GT. 1 .AND. .NOT. L_HRLYSIG(IS)) THEN
            WRITE(DUMMY,'(I10.10)') FULLHRQ
            CALL ERRHDL(PATH,MODNAM,'E','345',DUMMY)
            HRQS = 0.0D0
            RUNERR = .TRUE.
            GO TO 999
         END IF

         CALL STODBL(FIELD(8), ILEN_FLD, HRQS, IMIT)
         IF (IMIT .NE. 1) THEN
            CALL ERRHDL(PATH,MODNAM,'E','208','HOUREMIS')
            RUNERR = .TRUE.
            GO TO 999
C ---    Check for large negative values, could be a missing indicator
         ELSE IF ( HRQS .LE. -90.0D0 ) THEN
C*          Assume emissions are missing; assign value of 0.0 and issue Warning;
C*          Note that AERMOD User's Guide (p. 3-38, 3-39) implies that blanks 
C*          should be used for missing values.
            HRQS = 0.0D0  
            WRITE(DUMMY,'(I10.10)') FULLHRQ
            CALL ERRHDL(PATH,MODNAM,'W','341',DUMMY)
         END IF
        
         CALL STODBL(FIELD(9), ILEN_FLD, HRHS, IMIT)
         IF (IMIT .NE. 1) THEN
            CALL ERRHDL(PATH,MODNAM,'E','208','HOUREMIS')
            RUNERR = .TRUE.
            GO TO 999
         END IF

         CALL STODBL(FIELD(10), ILEN_FLD, HRSZ, IMIT)
         IF (IMIT .NE. 1) THEN
            CALL ERRHDL(PATH,MODNAM,'E','208','HOUREMIS')
            RUNERR = .TRUE.
            GO TO 999
         END IF

      ELSE IF ((SRCTYP(IS)(1:4).EQ.'AREA' .OR. SRCTYP(IS).EQ.'LINE')
     &                                    .AND. IFC.EQ.8) THEN
C*       Assign emission rate for AREA and LINE sources
C*       Check logical variable indicating hourly sigmas, L_HRLYSIG
         IF (L_HRLYSIG(IS)) THEN
C*          WRITE Error Message; Hourly Sigmas must be used for all hours
            WRITE(DUMMY,'(I10.10)') FULLHRQ
            CALL ERRHDL(PATH,MODNAM,'E','345',DUMMY)
            RUNERR = .TRUE.
            HRQS = 0.0D0
            GO TO 999
         END IF
         
         CALL STODBL(FIELD(8), ILEN_FLD, HRQS, IMIT)
         IF (IMIT .NE. 1) THEN
            CALL ERRHDL(PATH,MODNAM,'E','208','HOUREMIS')
            RUNERR = .TRUE.
            GO TO 999
C ---    Check for large negative values, could be a missing indicator
         ELSE IF ( HRQS .LE. -90.0D0 ) THEN
C*          Assume emissions are missing; assign value of 0.0 and issue Warning;
C*          Note that AERMOD User's Guide (p. 3-38, 3-39) implies that blanks 
C*          should be used for missing values.
            HRQS = 0.0D0  
            WRITE(DUMMY,'(I10.10)') FULLHRQ
            CALL ERRHDL(PATH,MODNAM,'W','341',DUMMY)
         END IF

      ELSE IF (SRCTYP(IS) .EQ. 'OPENPIT' .AND. IFC.EQ.8) THEN
C*       Assign emission rate for OPENPIT sources
         CALL STODBL(FIELD(8), ILEN_FLD, HRQS, IMIT)
         IF (IMIT .NE. 1) THEN
            CALL ERRHDL(PATH,MODNAM,'E','208','HOUREMIS')
            RUNERR = .TRUE.
            GO TO 999
C ---    Check for large negative values, could be a missing indicator
         ELSE IF ( HRQS .LE. -90.0D0 ) THEN
C*          Assume emissions are missing; assign value of 0.0 and issue Warning;
C*          Note that AERMOD User's Guide (p. 3-38, 3-39) implies that blanks 
C*          should be used for missing values.
            HRQS = 0.0D0  
            WRITE(DUMMY,'(I10.10)') FULLHRQ
            CALL ERRHDL(PATH,MODNAM,'W','341',DUMMY)
         END IF
  
      ELSE IF (SRCTYP(IS) .EQ. 'BUOYLINE' .AND. IFC.EQ.8) THEN
C*       Assign emission rate for BUOYANT LINE sources

         CALL STODBL(FIELD(8), ILEN_FLD, HRQS, IMIT)
         IF (IMIT .NE. 1) THEN
            CALL ERRHDL(PATH,MODNAM,'E','208','HOUREMIS')
            RUNERR = .TRUE.
            GO TO 999
C ---    Check for large negative values, could be a missing indicator
         ELSE IF ( HRQS .LE. -90.0D0 ) THEN
C*          Assume emissions are missing; assign value of 0.0 and issue Warning;
C*          Note that AERMOD User's Guide (p. 3-38, 3-39) implies that blanks
C*          should be used for missing values.
            HRQS = 0.0D0
            WRITE(DUMMY,'(I10.10)') FULLHRQ
            CALL ERRHDL(PATH,MODNAM,'W','341',DUMMY)
         END IF

      ELSE IF (SRCTYP(IS)(1:5).EQ.'POINT' .AND. IFC.GT.10) THEN
C*       Too many parameters - WRITE Error Message
C*       Assign zeros to all parameters
         WRITE(DUMMY,'(I10.10)') FULLHRQ
         CALL ERRHDL(PATH,MODNAM,'E','346',DUMMY)
         HRQS = 0.0D0
         HRTS = 0.0D0
         HRVS = 0.0D0
         RUNERR = .TRUE.
                       
      ELSE IF (SRCTYP(IS).EQ.'VOLUME' .AND.
     &        ((L_HRLYSIG(IS) .AND. IFC.GT.11) .OR.
     &    (.NOT.L_HRLYSIG(IS) .AND. IFC.GT.8))) THEN
C*       Too many parameters - WRITE Error Message
C*       Assign zeros to all parameters
         WRITE(DUMMY,'(I10.10)') FULLHRQ
         CALL ERRHDL(PATH,MODNAM,'E','346',DUMMY)
         HRQS = 0.0D0
         HRHS = 0.0D0
         HRSY = 0.0D0
         HRSZ = 0.0D0
         RUNERR = .TRUE.
                              
      ELSE IF ((SRCTYP(IS)(1:4).EQ.'AREA' .OR. SRCTYP(IS).EQ.'LINE')
     &         .AND. ((L_HRLYSIG(IS) .AND. IFC.GT.10) .OR.
     &           (.NOT.L_HRLYSIG(IS) .AND. IFC.GT.8))) THEN
C*       Too many parameters - WRITE Error Message
C*       Assign zeros to all parameters
         WRITE(DUMMY,'(I10.10)') FULLHRQ
         CALL ERRHDL(PATH,MODNAM,'E','346',DUMMY)
         HRQS = 0.0D0
         HRHS = 0.0D0
         HRSZ = 0.0D0
         RUNERR = .TRUE.
                              
      ELSE IF (SRCTYP(IS).EQ.'OPENPIT' .AND. IFC.GT.8) THEN
C*       Too many parameters - WRITE Error Message
C*       Assign zeros to all parameters
         WRITE(DUMMY,'(I10.10)') FULLHRQ
         CALL ERRHDL(PATH,MODNAM,'E','346',DUMMY)
         HRQS = 0.0D0
         RUNERR = .TRUE.
                              
      ELSE IF (SRCTYP(IS) .EQ. 'BUOYLINE' .AND. IFC.GT.8) THEN
C*       Too many parameters - WRITE Error Message
C*       Assign zeros to all parameters
         WRITE(DUMMY,'(I10.10)') FULLHRQ
         CALL ERRHDL(PATH,MODNAM,'E','346',DUMMY)
         HRQS = 0.0D0
         RUNERR = .TRUE.

      ELSE IF (SRCTYP(IS)(1:5) .EQ. 'POINT') THEN
C*       Some missing parameters - WRITE Error Message
C*       Assign zeros to all parameters
         WRITE(DUMMY,'(I10.10)') FULLHRQ
         CALL ERRHDL(PATH,MODNAM,'E','384',DUMMY)
         HRQS = 0.0D0
         HRTS = 0.0D0
         HRVS = 0.0D0
         RUNERR = .TRUE.

      ELSE IF (SRCTYP(IS) .EQ. 'VOLUME' .AND.
     &        ((L_HRLYSIG(IS) .AND. IFC.LT.11) .OR.
     &    (.NOT.L_HRLYSIG(IS) .AND. IFC.LT.8))) THEN
C*       Some missing parameters - WRITE Error Message
C*       Assign zeros to all parameters
         WRITE(DUMMY,'(I10.10)') FULLHRQ
         CALL ERRHDL(PATH,MODNAM,'E','384',DUMMY)
         HRQS = 0.0D0
         HRHS = 0.0D0
         HRSY = 0.0D0
         HRSZ = 0.0D0
         RUNERR = .TRUE.

      ELSE IF ((SRCTYP(IS)(1:4).EQ.'AREA' .OR. SRCTYP(IS).EQ.'LINE')
     &        .AND. ((L_HRLYSIG(IS) .AND. IFC.LT.10) .OR.
     &    (.NOT.L_HRLYSIG(IS) .AND. IFC.LT.8))) THEN
C*       Some missing parameters - WRITE Error Message
C*       Assign zeros to all parameters
         WRITE(DUMMY,'(I10.10)') FULLHRQ
         CALL ERRHDL(PATH,MODNAM,'E','384',DUMMY)
         HRQS = 0.0D0
         HRHS = 0.0D0
         HRSZ = 0.0D0
         RUNERR = .TRUE.

      ELSE IF (SRCTYP(IS) .EQ. 'BUOYLINE') THEN
C*       Some missing parameters - WRITE Error Message
C*       Assign zeros to all parameters
         WRITE(DUMMY,'(I10.10)') FULLHRQ
         CALL ERRHDL(PATH,MODNAM,'E','384',DUMMY)
         HRQS = 0.0D0
         RUNERR = .TRUE.

      ELSE
C*       Problem processing HOUREMIS record - WRITE Error Message
C*       Assign zeros to emission rate
         WRITE(DUMMY,'(I10.10)') FULLHRQ
         CALL ERRHDL(PATH,MODNAM,'E','345',DUMMY)
         HRQS = 0.0D0
         RUNERR = .TRUE.

      END IF

      GO TO 999

C*    Write Error Message for Error Reading Hourly Emissions File
 99   CALL ERRHDL(PATH,MODNAM,'E','510','HOUREMIS')
      RUNERR = .TRUE.
      GO TO 999

888   CONTINUE

      EOF = .TRUE.

999   RETURN
      END
      
      SUBROUTINE HRQEXT (IS)
C***********************************************************************
C*                  HRQEXT Module of AERMOD
C* 
C*         PURPOSE: To Assign Hourly Source Parameters
C* 
C*         PROGRAMMER:  Jayant Hardikar, Roger Brode
C* 
C*         DATE:    September 15, 1993
C* 
C*         INPUTS:  Current Source Number Being Processed
C* 
C*         OUTPUTS: Source Arrays
C*
C*         Revision History:
C*
C*         MODIFIED:  Incorporated options to specify hourly-varying
C*                    release heights and initial dispersion coefficients
C*                    for VOLUME and AREA sources.
C*                    R.W. Brode, U.S. EPA/OAQPS/AQMG, 10/19/2009
C*         
C*         MODIFIED:  Corrected processing of missing parameters for 
C*                    point sources to assign all parameters to 0.0 if
C*                    any of the parameters are missing, in conformance
C*                    with Section 3.3.9 of the AERMOD User's Guide.
C*                    R.W. Brode, U.S. EPA/OAQPS/AQMG, 01/24/2007
C*         
C*         MODIFIED:  REMOVED THE 'POINT' SOURCE CONDITION, SO IT APPLIES 
C*                    TO ALL SOURCE TYPES, EXCEPT SAVING THE TEMP & VEL
C* 
C*         CALLED FROM:  HRLOOP
C************************************************************************
C*
C*    Variable Declarations
      USE MAIN1
      USE BUOYANT_LINE
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: IS

C*    Variable Initializations
      MODNAM = 'HRQEXT'
      
C*    Assign the Hourly Emission Parameters to the appropriate arrays
      IF (EVONLY) THEN

         AQS(IS) = EV_HRQS(IS,IHOUR)

         IF (SRCTYP(IS)(1:5) .EQ. 'POINT') THEN
            ATS(IS) = EV_HRTS(IS,IHOUR)
            AVS(IS) = EV_HRVS(IS,IHOUR)
         ELSE IF (SRCTYP(IS) .EQ. 'VOLUME' .AND. L_HRLYSIG(IS)) THEN
            AHS(IS)    = EV_HRHS(IS,IHOUR)
            ASYINI(IS) = EV_HRSY(IS,IHOUR)
            ASZINI(IS) = EV_HRSZ(IS,IHOUR)
         ELSE IF (SRCTYP(IS)(1:4) .EQ. 'AREA' .AND. L_HRLYSIG(IS)) THEN
            AHS(IS)    = EV_HRHS(IS,IHOUR)
            ASZINI(IS) = EV_HRSZ(IS,IHOUR)
         ELSE IF (SRCTYP(IS) .EQ. 'LINE' .AND. L_HRLYSIG(IS)) THEN
            AHS(IS)    = EV_HRHS(IS,IHOUR)
            ASZINI(IS) = EV_HRSZ(IS,IHOUR)
         END IF

      ELSE
      
         AQS(IS) = HRQS

         IF (SRCTYP(IS)(1:5) .EQ. 'POINT') THEN
            ATS(IS) = HRTS
            AVS(IS) = HRVS
         ELSE IF (SRCTYP(IS) .EQ. 'VOLUME' .AND. L_HRLYSIG(IS)) THEN
            AHS(IS)    = HRHS
            ASYINI(IS) = HRSY
            ASZINI(IS) = HRSZ
         ELSE IF (SRCTYP(IS)(1:4) .EQ. 'AREA' .AND. L_HRLYSIG(IS)) THEN
            AHS(IS)    = HRHS
            ASZINI(IS) = HRSZ
         ELSE IF (SRCTYP(IS) .EQ. 'LINE' .AND. L_HRLYSIG(IS)) THEN
            AHS(IS)    = HRHS
            ASZINI(IS) = HRSZ
         END IF

      END IF
      
C*    Perform QA Error Checking on Source Parameters
C       Since the only buoyant line source parameter that can vary hourly 
C       is the emission rate, there is nothing to check here.
C*

      IF (SRCTYP(IS)(1:5) .EQ. 'POINT') THEN
         IF (ATS(IS) .EQ. 0.0D0) THEN
C*          Set Temperature to Small Negative Value for Ambient Releases
            ATS(IS) = -1.0D-5
         ELSE IF (ATS(IS) .GT. 2000.0D0) THEN
C*          WRITE Warning Message:  Exit Temp. > 2000K
            IF (.NOT. L_SkipMessages) THEN
               CALL ERRHDL(PATH,MODNAM,'W','320','HRTS')
            END IF
         ELSE IF ( (DABS(AQS(IS)).GT.0.0D0) .AND. 
     &             (ATS(IS).GT.0.0D0) .AND. (ATS(IS).LT.200.0D0) .AND.
     &                                      (AVS(IS).GT.200.0D0) ) THEN
C*          Exit temp < 200K (about -100F) and exit velocity > 200m/s
C*          with non-zero emissions; Incorrect units may have been 
C*          used or ATS and AVS may have been switched;  
C*          WRITE Fatal Error Message
            IF (.NOT. L_SkipMessages) THEN
               CALL ERRHDL(PATH,MODNAM,'E','320','HRTS')
               RUNERR = .TRUE.
            END IF
         ELSE IF ( (DABS(AQS(IS)).GT.0.0D0) .AND. 
     &             (ATS(IS).GT.0.0D0) .AND. (ATS(IS).LT.200.0D0) ) THEN
C*          Exit temp < 200K (about -100F) with non-zero emissions; 
C*          Incorrect units may have been used or ATS and AVS may 
C*          have been switched;  
C*          WRITE Warnign Message
            IF (.NOT. L_SkipMessages) THEN
               CALL ERRHDL(PATH,MODNAM,'W','320','HRTS')
            END IF
         END IF

         IF (AVS(IS) .LT. 0.0D0) THEN
C*          WRITE Warning Message:  Negative or Zero Exit Velocity
            IF (.NOT. L_SkipMessages) THEN
               CALL ERRHDL(PATH,MODNAM,'W','325',SRCID(IS))
            END IF
C*          Set to Small Value to Avoid Zero-divide and Underflow
            AVS(IS) = 1.0D-5
         ELSE IF (AVS(IS) .LT. 1.0D-5) THEN
C*          Set to Small Value to Avoid Zero-divide and Underflow
            AVS(IS) = 1.0D-5
         ELSE IF (AVS(IS) .GT. 50.0D0) THEN
C*          WRITE Informational Message:  Exit Velocity > 50.0 m/s
            IF (.NOT. L_SkipMessages) THEN
               CALL ERRHDL(PATH,MODNAM,'I','320','HRVS')
            END IF
         END IF
         
      ELSE IF (SRCTYP(IS) .EQ. 'VOLUME') THEN        
         IF (AHS(IS) .LT. 0.0D0) THEN
C           WRITE Error Message:  Negative Release Height
            CALL ERRHDL(PATH,MODNAM,'E','209','HRHS')
         ELSE IF (AHS(IS) .GT. 100.0D0) THEN
C           WRITE Warning Message:  Large Release Height (> 100M)
            IF (.NOT. L_SkipMessages) THEN
               CALL ERRHDL(PATH,MODNAM,'W','320','HRHS')
            END IF
         ELSE IF (AHS(IS) .GT. 3000.0D0) THEN
C           WRITE Error Message:  Large Release Height (> 3000M)
            CALL ERRHDL(PATH,MODNAM,'E','324',SRCID(IS))
            RUNERR = .TRUE.
         END IF
 
         IF (ASYINI(IS) .LT. 0.0D0) THEN
C           WRITE Warning Message:  Negative Initial Lateral Parameter
            CALL ERRHDL(PATH,MODNAM,'E','209','HRSY')
C           Set to Small Value to Avoid Zero-divide and Underflow
            ASYINI(IS) = 1.0D-5
         ELSE IF (ASYINI(IS) .LT. 1.0D-5) THEN
C           WRITE Warning Message:  Small Initial Lateral Parameter
            IF (.NOT. L_SkipMessages) THEN
               CALL ERRHDL(PATH,MODNAM,'W','320','HRSY')
            END IF
C           Set to Small Value to Avoid Zero-divide and Underflow
            ASYINI(IS) = 1.0D-5
         ELSE IF (ASYINI(IS) .GT. 200.0D0) THEN
C           WRITE Warning Message:  Large Initial Lateral Parameter (> 200m)
            IF (.NOT. L_SkipMessages) THEN
               CALL ERRHDL(PATH,MODNAM,'W','320','HRSY')
            END IF
         END IF
 
         IF (ASZINI(IS) .LT. 0.0D0) THEN
C           WRITE Warning Message:  Negative Initial Vertical Parameter
            CALL ERRHDL(PATH,MODNAM,'E','209','HRSZ')
C           Set to Small Value to Avoid Zero-divide and Underflow
            ASZINI(IS) = 1.0D-5
         ELSE IF (ASZINI(IS) .LT. 1.0D-5) THEN
C           WRITE Warning Message:  Small Initial Lateral Parameter
            IF (.NOT. L_SkipMessages) THEN
               CALL ERRHDL(PATH,MODNAM,'W','320','HRSZ')
            END IF
C           Set to Small Value to Avoid Zero-divide and Underflow
            ASZINI(IS) = 1.0D-5
         ELSE IF (ASZINI(IS) .GT. 200.0D0) THEN
C           WRITE Warning Message:  Large Initial Vertical Parameter (> 200m)
            IF (.NOT. L_SkipMessages) THEN
               CALL ERRHDL(PATH,MODNAM,'W','320','HRSZ')
            END IF
         END IF
         
      ELSE IF (SRCTYP(IS)(1:4).EQ.'AREA' .OR. SRCTYP(IS).EQ.'LINE') THEN
         IF (AHS(IS) .LT. 0.0D0) THEN
C           WRITE Error Message:  Negative Release Height
            CALL ERRHDL(PATH,MODNAM,'E','209','HRHS')
         ELSE IF (AHS(IS) .GT. 100.0D0) THEN
C           WRITE Warning Message:  Large Release Height (> 100M)
            IF (.NOT. L_SkipMessages) THEN
               CALL ERRHDL(PATH,MODNAM,'W','320','HRHS')
            END IF
         ELSE IF (AHS(IS) .GT. 3000.0D0) THEN
C           WRITE Error Message:  Large Release Height (> 3000M)
            CALL ERRHDL(PATH,MODNAM,'E','324',SRCID(IS))
            RUNERR = .TRUE.
         END IF
 
         IF (ASZINI(IS) .LT. 0.0D0) THEN
C           WRITE Warning Message:  Negative Initial Vertical Parameter
            CALL ERRHDL(PATH,MODNAM,'E','209','HRSZ')
C           Set to Small Value to Avoid Zero-divide and Underflow
            ASZINI(IS) = 1.0D-5
         ELSE IF (ASZINI(IS) .LT. 1.0D-5) THEN
C           WRITE Warning Message:  Small Initial Lateral Parameter
            IF (.NOT. L_SkipMessages) THEN
               CALL ERRHDL(PATH,MODNAM,'W','320','HRSZ')
            END IF
C           Set to Small Value to Avoid Zero-divide and Underflow
            ASZINI(IS) = 1.0D-5
         ELSE IF (ASZINI(IS) .GT. 200.0D0) THEN
C           WRITE Warning Message:  Large Initial Vertical Parameter (> 200m)
            IF (.NOT. L_SkipMessages) THEN
               CALL ERRHDL(PATH,MODNAM,'W','320','HRSZ')
            END IF
         END IF
         
      END IF

999   RETURN
      END

      SUBROUTINE O3EXT(L_ReadErr,ReadErr,EndErr)
C***********************************************************************
C*                  O3EXT Module of AERMOD
C*
C*         PURPOSE: To extract hourly ozone data for PVMRM and OLM options
C*
C*         PROGRAMMER:  Roger W. Brode, PES, Inc.
C*
C*         DATE:    May 6, 2002
C*
C*         INPUTS:
C*
C*         OUTPUTS:
C*
C*         CALLED FROM:  HRLOOP
C************************************************************************
C*
C*    Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      CHARACTER ReadErr*5, EndErr*5

      DOUBLE PRECISION :: O3SUB(6), O3TEMP, O3MAX24, O3MIN

      INTEGER :: IO3YR, IO3MN, IO3DY, IO3HR, IO3YR2, I
      INTEGER :: FULLO3HR(6)
      LOGICAL :: L_ReadErr

C*    Variable Initializations
      MODNAM  = 'O3EXT'
      ReadErr = ''
      EndErr  = ''
      L_ReadErr = .FALSE.
      FULLO3HR(:) = 0

C --- Initialize logical for missing data to FALSE
      O3MISS = .FALSE.

C*    Assign IOLINE counter to ILINE for passing to ERRHDL if needed
      ILINE = IOLINE
C*
C --- Read a record in all hourly O3 files, but only read concentation 
C     based on the applicable file for the current sector
C     First initialize O3CONC to 0.0, and O3TEMP, O3SUB and O3MIN to -99.
      O3CONC   = 0.0D0
      O3TEMP   = -99.0D0
      O3SUB(:) = -99.0D0
      O3MIN    = -99.0D0

      DO I = 1, NUMO3Sects
C ---    Loop through O3SECTORs

C ---    Reinitialize O3SUB for this sector
         O3SUB(I) = -99.0D0

C ---    Check for non-hourly O3 values to substitute for missing data
         IF (L_O3VALUES(I)) THEN
            CALL OZONVALS(I,O3SUB(I))
         ELSE IF (L_O3VAL(I)) THEN
            O3SUB(I) = O3BACK(I)
         ELSE
            O3SUB(I) = 0.0D0
         END IF

C ---    Check for hour O3FILE for this sector
         IF (L_O3File(I)) THEN

C ---       Hourly O3 file available for current sector

            IF (I .EQ. IO3SECT) THEN
C ---          This is the applicable sector for this hour; read next hour of O3 data

C ---          Initialize ReadErr and EndErr for this sector
               WRITE(ReadErr,'(''SECT'',I1)') I

               IF (O3FORM(I) .EQ. 'FREE') THEN
                  READ(IO3UNT(I),*,ERR=99,END=9991) IO3YR, IO3MN, 
     &                                              IO3DY, IO3HR, 
     &                                              O3CONC
               ELSE
                  READ(IO3UNT(I),O3FORM(I),ERR=99,END=9991) 
     &                                              IO3YR, IO3MN, 
     &                                              IO3DY, IO3HR,
     &                                              O3CONC
               END IF

C ---          Check for use of 2-digit year in OZONEFIL file, adjust to 4-digit
C              year for comparison with FULLDATE based on met data file
               IF (IO3YR .LE. 99) THEN
                  IO3YR2 = IO3YR
                  IF (IO3YR2 .GE. ISTRT_WIND .AND. 
     &                                 IO3YR2 .LE. 99) THEN
                     IO3YR  = ISTRT_CENT*100 + IO3YR2
                  ELSE IF (IO3YR2 .LT. ISTRT_WIND) THEN
                     IO3YR  = (ISTRT_CENT+1)*100 + IO3YR2
                  END IF
               END IF
           
C ---          Calculate full date for this hour of O3 data
               FULLO3HR(I) = IO3YR*1000000 + IO3MN*10000 + IO3DY*100
     &                                                   + IO3HR

               IF (O3CONC .GE. 0.0D0 .AND. O3CONC .LT. 900.0D0) THEN
C ---             Valid hourly value; convert to ug/m3 if needed
                  IF (O3FILUNITS .EQ. 'PPB') THEN
                     O3CONC = O3CONC * O3_PPB
                  ELSE IF (O3FILUNITS .EQ. 'PPM') then
                     O3CONC = O3CONC * O3_PPM
                  END IF
C ---             Valid hourly O3 value; check for application of O3MIN value
                  IF (STABLE) THEN
C                    Use min of 40 ppb (78.4ug/m3) and max from previous 24 hrs
                     O3MAX24 = MIN ( 78.40D0, 
     &                         MAXVAL( O3_Max24hr(:,IO3SECT) ) )
C                    Adjust minimum O3 value based on OBULEN
                     IF (OBULEN .GT.  0.0D0 .AND. 
     &                   OBULEN .LE. 50.0D0) THEN
                        O3MIN = O3MAX24
                     ELSE IF (OBULEN .GT. 250.0D0) THEN
                        O3MIN = 0.0D0
                     ELSE
                        O3MIN = O3MAX24 * (250.D0 - OBULEN) /200.D0
                     END IF
                  ELSE
                     O3MIN = -9.0D0
                  END IF
C ---             Save this hour's O3CONC (in ug/m3) to array of previous 
C                 24 values, before applying minimum value
                  O3_Max24hr(IO3HR,IO3SECT) = O3CONC
                  O3CONC = MAX( O3CONC, O3MIN )
               ELSE IF (L_O3VALUES(IO3SECT) .OR.
     &                     L_O3VAL(IO3SECT)) THEN
C ---             Hourly O3 value is missing; assign O3SUB value based on
C                 O3VALUES or OZONEVAL inputs
                  O3CONC = O3SUB(IO3SECT)
C ---             Assign 0.0 to O3_Max24hr array for this hour
                  O3_Max24hr(IO3HR,IO3SECT) = 0.0D0
               ELSE
C ---             Assign O3MISS logical to TRUE
                  O3MISS = .TRUE.
C ---             Assign 0.0 to O3_Max24hr array for this sector
                  O3_Max24hr(IO3HR,IO3SECT) = 0.0D0
               END IF

               GO TO 9992

9991           CONTINUE
C              End-of-file reached, set logical flag
               EOF = .TRUE.

               WRITE(EndErr, '(''SECT'',I1)') I           

9992           CONTINUE

            ELSE
C ---          This is not applicable sector for this hour; read record with temp data

C ---          Initialize ReadErr and EndErr for this sector
               WRITE(ReadErr,'(''SECT'',I1)') I

               IF (O3FORM(I) .EQ. 'FREE') THEN
                  READ(IO3UNT(I),*,ERR=99,END=9993) IO3YR, IO3MN, 
     &                                              IO3DY, IO3HR,
     &                                              O3TEMP
               ELSE
                  READ(IO3UNT(I),O3FORM(I),ERR=99,END=9993) 
     &                                              IO3YR, IO3MN, 
     &                                              IO3DY, IO3HR,
     &                                              O3TEMP
               END IF

C ---          Check for use of 2-digit year in OZONEFIL file, adjust to 4-digit
C              year for comparison with FULLDATE based on met data file
               IF (IO3YR .LE. 99) THEN
                  IO3YR2 = IO3YR
                  IF (IO3YR2 .GE. ISTRT_WIND .AND. 
     &                                 IO3YR2 .LE. 99) THEN
                     IO3YR  = ISTRT_CENT*100 + IO3YR2
                  ELSE IF (IO3YR2 .LT. ISTRT_WIND) THEN
                     IO3YR  = (ISTRT_CENT+1)*100 + IO3YR2
                  END IF
               END IF
           
C ---          Calculate full date for this hour of O3 data
               FULLO3HR(I) = IO3YR*1000000 + IO3MN*10000 + IO3DY*100
     &                                                   + IO3HR

               IF (O3TEMP .GE. 0.0D0 .AND. O3TEMP .LT. 900.0D0) THEN
C ---             Valid hourly value; convert to ug/m3 if needed
                  IF (O3FILUNITS .EQ. 'PPB') THEN
                     O3TEMP = O3TEMP * O3_PPB
                  ELSE IF (O3FILUNITS .EQ. 'PPM') then
                     O3TEMP = O3TEMP * O3_PPM
                  END IF
C ---             Save this hour's O3CONC (in ug/m3) to array of previous 
C                 24 values for this sector
                  O3_Max24hr(IO3HR,I) = O3TEMP
               ELSE IF (L_O3VALUES(I) .OR.
     &                     L_O3VAL(I)) THEN
C ---             Hourly O3 value is missing; assign O3SUB value;
C                 these have already been converted to ug/m3
                  O3TEMP = O3SUB(I)
C ---             Assign 0.0 to O3_Max24hr array so that substituted value will
C                 not be used in determining max value from previous 24 hours
                  O3_Max24hr(IO3HR,I) = 0.0D0
               ELSE
C ---             Assign 0.0 to O3_Max24hr array
                  O3_Max24hr(IO3HR,I) = 0.0D0
               END IF

               GO TO 9994

9993           CONTINUE
C              End-of-file reached, set logical flag
               EOF = .TRUE.

               WRITE(EndErr, '(''SECT'',I1)') I             

9994           CONTINUE

            END IF

         END IF

      END DO   ! END of O3Sector Loop

      IF (O3MISS) THEN
C ---    No O3 value available for this hour; assign 0.0 to O3CONC 
C        and issue informational message
         O3CONC = 0.0D0
         IF (.NOT. L_SkipMessages) THEN
            WRITE(DUMMY,'(I10.10)') FULLDATE
            CALL ERRHDL(PATH,MODNAM,'I','459',DUMMY)
         END IF
      END IF

C*    Check for Date and Time Consistency Across all Sectors; If Failed, Issue Fatal Error
      DO I = 1, NUMO3Sects
         IF (L_O3File(I)) THEN
            IF (FULLDATE .NE. FULLO3HR(I)) THEN
C*             WRITE Error Message - Date mismatch
               WRITE(DUMMY,'(I10.10,''S'',I1)') FULLDATE, I
               CALL ERRHDL(PATH,MODNAM,'E','457',DUMMY)
C*             Set RUNERR logical and skip to end
               RUNERR = .TRUE.
               GO TO 1000
            END IF
         END IF
      END DO

C     Date/Time consistency checks were ok; skip to end
      GO TO 1000

C*    Write Error Message for Error Reading Hourly Ozone File
 99   CONTINUE
      L_ReadErr = .TRUE.
      WRITE(ReadErr,'(''SECT'',I1)') I
      CALL ERRHDL(PATH,MODNAM,'E','510',ReadErr)
      RUNERR = .TRUE.
      
      GO TO 1000

 999  CONTINUE
C     End-of-file reached, set logical flag
      EOF = .TRUE.
 
C --- End of file reached on O3 file; check FULLDATE vs. FULLO3HR before returning
C*    for Date and Time Consistency Across all Sectors; If Failed, Issue Fatal Error
      DO I = 1, NUMO3Sects
         IF (L_O3File(I)) THEN
            IF (FULLDATE .NE. FULLO3HR(I)) THEN
C*             WRITE Error Message - Date mismatch
               WRITE(DUMMY,'(I10.10,''S'',I1)') FULLDATE, I
               CALL ERRHDL(PATH,MODNAM,'E','457',DUMMY)
               RUNERR = .TRUE.
            END IF
         END IF
      END DO

1000  RETURN
      END

      SUBROUTINE BGEXT(L_ReadErr,ReadErr,EndErr)
C***********************************************************************
C*                  BGEXT Module of AERMOD
C*
C*         PURPOSE: To extract hourly background concentrations
C*
C*         PROGRAMMER:  Roger W. Brode
C*
C*         DATE:        February 28, 2011
C*
C*         INPUTS:
C*
C*         OUTPUTS:
C*
C*         CALLED FROM:  HRLOOP
C************************************************************************
C*
C*    Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      CHARACTER ReadErr*5, EndErr*5

      INTEGER :: IBGYR, IBGMN, IBGDY, IBGHR, IBGYR2, I
      INTEGER :: FULLBGHR(6)
      DOUBLE PRECISION :: BGHrVal, BGSUB(6)
      LOGICAL :: L_ReadErr

C*    Variable Initializations
      MODNAM  = 'BGEXT'
      ReadErr = ''
      EndErr  = ''
      L_ReadErr = .FALSE.

C*    Assign IBLINE counter to ILINE for passing to ERRHDL if needed
      ILINE = IBLINE
C*
C --- Read a record in all hourly BG files, but only read concentation 
C     based on the applicable file for the current sector
C     First initialize BGHrVal to -99.0
      BGHrVal = -99.0D0

C --- Loop through all BACKGRND Sectors
      DO I = 1, NUMBGSects

C ---    Reinitialize BGSUB for this sector
         BGSUB(I) = 0.0D0

C ---    Check for temporally-varying background to substitute for missing hours
         IF (L_BGValues(I)) THEN
            CALL BGVAL(I,BGSUB(I))
         ELSE
C           No temporally-varying background for this sector; set BGSUB = 0.0
            BGSUB(I) = 0.0D0
         END IF

C ---    Check for HOURLY BACKGRDN data for this sector and read the data
         IF (L_BGFile(I)) THEN

C ---       Hourly BG file available for current sector

C ---       Check for whether this is the applicable sector for this hour
            IF (I .EQ. IBGSECT) THEN

C ---          Initialize ReadErr and EndErr for this sector
               WRITE(ReadErr,'(''SECT'',I1)') I

               IF (BGFORM(I) .EQ. 'FREE') THEN
                  READ(IBGUNT(I),*,ERR=99,END=9991) IBGYR,IBGMN,
     &                                              IBGDY,IBGHR,
     &                                              BGHrVal
               ELSE
                  READ(IBGUNT(I),BGFORM(I),ERR=99,END=9991) 
     &                                              IBGYR,IBGMN, 
     &                                              IBGDY,IBGHR, 
     &                                              BGHrVal
               END IF

               IF (BGHrVal .GT. 0.0D0) THEN
C ---             Valid hourly value; convert to ug/m3 if needed
                  IF (POLLUT .EQ. 'NO2') THEN
                     IF (BackUnits .EQ. 'PPB') THEN
                        BGHrVal = BGHrVal / NO2_PPB
                     ELSE IF (BackUnits .EQ. 'PPM') THEN
                        BGHrVal = BGHrVal / NO2_PPM
                     END IF
                  ELSE IF (POLLUT .EQ. 'SO2') THEN
                     IF (BackUnits .EQ. 'PPB') THEN
                        BGHrVal = BGHrVal / SO2_PPB
                     ELSE IF (BackUnits .EQ. 'PPM') THEN
                        BGHrVal = BGHrVal / SO2_PPM
                     END IF
                  ELSE IF (POLLUT .EQ. 'CO') THEN
                     IF (BackUnits .EQ. 'PPB') THEN
                        BGHrVal = BGHrVal * CO_PPB
                     ELSE IF (BackUnits .EQ. 'PPM') THEN
                        BGHrVal = BGHrVal * CO_PPM
                     END IF
                  END IF
               END IF

C ---          Check for use of 2-digit year in BACKGRND file, adjust to 4-digit
C              year for comparison with FULLDATE based on met data file
               IF (IBGYR .LE. 99) THEN
                  IBGYR2 = IBGYR
                  IF (IBGYR2 .GE. ISTRT_WIND .AND. 
     &                                 IBGYR2 .LE. 99) THEN
                     IBGYR  = ISTRT_CENT*100 + IBGYR2
                  ELSE IF (IBGYR2 .LT. ISTRT_WIND) THEN
                     IBGYR  = (ISTRT_CENT+1)*100 + IBGYR2
                  END IF
               END IF

C*             Assign full date for this HOURLY BACKGRND file
               FULLBGHR(I) = IBGYR*1000000 + IBGMN*10000 + IBGDY*100 +
     &                                                     IBGHR

               GO TO 9992

9991           CONTINUE
C              End-of-file reached, set logical flag
               EOF = .TRUE.

               WRITE(EndErr, '(''SECT'',I1)') I              

9992           CONTINUE

            ELSE
C ---          This is not applicable sector for this hour; read record without data
C
C ---          Initialize ReadErr and EndErr for this sector
               WRITE(ReadErr,'(''SECT'',I1)') I

               IF (BGFORM(I) .EQ. 'FREE') THEN
                  READ(IBGUNT(I),*,ERR=99,END=9993) IBGYR, IBGMN,
     &                                              IBGDY, IBGHR
               ELSE
                  READ(IBGUNT(I),BGFORM(I),ERR=99,END=9993) 
     &                                              IBGYR, IBGMN, 
     &                                              IBGDY, IBGHR
               END IF

C ---          Check for use of 2-digit year in BACKGRND file, adjust to 4-digit
C              year for comparison with FULLDATE based on met data file
               IF (IBGYR .LE. 99) THEN
                  IBGYR2 = IBGYR
                  IF (IBGYR2 .GE. ISTRT_WIND .AND. 
     &                                 IBGYR2 .LE. 99) THEN
                     IBGYR  = ISTRT_CENT*100 + IBGYR2
                  ELSE IF (IBGYR2 .LT. ISTRT_WIND) THEN
                     IBGYR  = (ISTRT_CENT+1)*100 + IBGYR2
                  END IF
               END IF

C*             Assign full date for this HOURLY BACKGRND file
               FULLBGHR(I) = IBGYR*1000000 + IBGMN*10000 + IBGDY*100 +
     &                                                     IBGHR

               GO TO 9994

9993           CONTINUE
C              End-of-file reached, set logical flag
               EOF = .TRUE.

               WRITE(EndErr, '(''SECT'',I1)') I         

9994           CONTINUE

            END IF

         END IF   ! End of L_BGFile IF Block

      END DO      ! End of BGSector Loop

C*    Check for Date and Time Consistency Across all Sectors; If Failed, Issue Fatal Error
      DO I = 1, NUMBGSects
         IF (L_BGFile(I)) THEN
            IF (FULLDATE .NE. FULLBGHR(I)) THEN
C*             WRITE Error Message - Date mismatch
               WRITE(DUMMY,'(I10.10,''s'',I1)') FULLDATE, I
               CALL ERRHDL(PATH,MODNAM,'E','454',DUMMY)
               RUNERR = .TRUE.
            END IF
         END IF
      END DO

      IF (RUNERR) RETURN

      IF (BGHrVal .LT. 0.0D0) THEN
C ---    Hourly BGCONC is missing; look for substitution values
         IF (IBGSECT .GT. 0) THEN
C ---       Valid BGSECT defined, check for hourly values for this 
C           sector, and then for non-hourly values to substitute
            IF (L_BGFile(IBGSECT)) THEN
               IF (L_BGValues(IBGSECT)) THEN
C                 Hourly background value is missing but non-hourly 
C                 values have been specified for substitution, 
C                 which were processed in subroutine BGVAL;
                  BGCONC = BGSUB(IBGSECT)
C                 Write informational message 
                  WRITE(DUMMY,'(I10.10,''s'',I1)') FULLDATE, IBGSECT
                  CALL ERRHDL(PATH,MODNAM,'I','453',DUMMY)
C                 Increment counter for number of missing BGval substitutions
                  NSubBGHOUR = NSubBGHOUR + 1
               ELSE
C                 Hourly background value is missing for this sector and no 
C                 non-hourly values specified for substitution;
C                 Write Error message 
                  WRITE(DUMMY,'(I10.10,''s'',I1)') FULLDATE, IBGSECT
                  CALL ERRHDL(PATH,MODNAM,'E','452',DUMMY)
                  RUNERR = .TRUE.
                  GO TO 1000
               END IF
            ELSE
               IF (L_BGValues(IBGSECT)) THEN
C                 Hourly background value is missing but non-hourly 
C                 values have been specified for substitution, 
C                 which were processed in subroutine BGVAL;
                  BGCONC = BGSUB(IBGSECT)
               END IF
            END IF
         ELSE
C ---       IBGSECT = 0 due to calm or msg hr; BGSUB = 0.0D0; exit
            BGCONC = 0.0D0
         END IF
      ELSE
         BGCONC = BGHrVal
      END IF

      RETURN

C*    Write Error Message for Error Reading Hourly Background File
 99   CONTINUE
      L_ReadErr = .TRUE.
      
      RETURN

 999  CONTINUE
C     End-of-file reached, set logical flag
      EOF = .TRUE.
 
C --- End of file reached on hourly BACKGRND file; check FULLDATE vs. FULLBGHR
C     before returning
      DO I = 1, NUMBGSects
         IF (L_BGFile(I)) THEN
            IF (FULLDATE .NE. FULLBGHR(I)) THEN
C*             WRITE Error Message - Date mismatch
               WRITE(DUMMY,'(I10.10,''S'',I1)') FULLDATE, I
               CALL ERRHDL(PATH,MODNAM,'E','454',DUMMY)
               RUNERR = .TRUE.
            END IF
         END IF
      END DO

1000  RETURN
      END

      SUBROUTINE ERRHDL(PATHWY,MODNAM,INERTP,INERCD,INPMSG)
C***********************************************************************
C                 ERRHDL Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: A General Error Handling Procedure
C
C        PROGRAMMER: Jeff Wang
C
C        DATE:    March 2, 1992
C
C        MODIFIED:  Sets upper limit on line number included in error
C                   message to avoid overflowing the field; also increased
C                   field length for last message field from 8 to 12 to
C                   accommodate 12 character source IDs.
C                   R.W. Brode, U.S. EPA/OAQPS/AQMG, 02/28/2011
C         
C        INPUTS:  Error Code, Occur Locations
C
C        OUTPUTS: Error Message, Error Statistics..etc.
C
C        CALLED FROM:  (This Is An Utility Programm)
C***********************************************************************
C
C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE

      INTEGER :: I, ILINE_PRT
      CHARACTER ERRMG1*50, PATHWY*2, INERTP*1, INERCD*3, ICODE*3,
     &          INPMSG*(*), MODNAM*(*), TMPMOD*6, TMPMSG*12
      LOGICAL FOUND

C     Variable Initializations
      IERROR = IERROR + 1
      FOUND = .FALSE.
      I = 1

C     Check for Occurrence of 'E' Error Type, and Set FATAL Switch
      IF (INERTP .EQ. 'E') THEN
         FATAL = .TRUE.
         NFATAL = NFATAL + 1
         IF (NFATAL .EQ. 999) THEN
C           Number Of Fatal Errors Has Reached Limit of 999
            ERRMG1 = 'Number of Fatal Errors Has Reached Limit of 999'
            TMPMOD = 'ERRHDL'
            ICODE  = '999'
            TMPMSG = ' '
            ILINE_PRT = MIN(ILINE,99999999)
            WRITE(IERUNT,1111) PATHWY,INERTP,ICODE,ILINE_PRT,TMPMOD,
     &                         ERRMG1,TMPMSG
            GO TO 999
         ELSE IF (NFATAL .GT. 999) THEN
C           Skip Any More Error WRITEs
            GO TO 999
         END IF
      END IF

C     Go To Match The Error Massage
      DO WHILE (.NOT.FOUND .AND. I.LE.IERRN)
         IF (INERCD .EQ. ERRCOD(I)) THEN
            ERRMG1 = ERRMSG(I)
            FOUND = .TRUE.
         END IF
         I = I + 1
      END DO

      IF (.NOT. FOUND) THEN
         WRITE(ERRMG1,1001)
 1001    FORMAT('SYSTEM ERROR: MESSAGE NOT FOUND FOR THIS NUMBER!')
      END IF

C --- Set upper limit on ILINE to avoid write error
      ILINE_PRT = MIN(ILINE,99999999)
C     Write Out The Error Message
      WRITE(IERUNT,1111) PATHWY,INERTP,INERCD,ILINE_PRT,
     &                   MODNAM(1:MIN(LEN_TRIM(MODNAM),12)),ERRMG1,
     &                   INPMSG(1:MIN(LEN_TRIM(INPMSG),12))
 1111 FORMAT(A2,1X,A1,A3,I8,1X,A12,': ',A50,1X,A12)

 999  RETURN
      END

      SUBROUTINE TERRST
C***********************************************************************
C                 TERRST Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: To Determine Total Error/Message Statistics
C
C        PROGRAMMER:  Jeff Wang, Roger Brode
C
C        DATE:    March 2, 1992
C
C        MODIFIED:  Corrected issues with determining number of calm 
C                   and/or missing hours only for the data period that
C                   is processed.  Also increased field length for last 
C                   message field from 8 to 12 to accommodate 12 character 
C                   source IDs.
C                   R.W. Brode, U.S. EPA/OAQPS/AQMG, 02/28/2011
C         
C        MODIFIED:  To remove reference to legacy ISCST3 option (HE>ZI)
C                   for plume height above mixing height.
C                   Determine number of calm and/or missing hours only
C                   for the data period processed, and determine number
C                   of hours processed.
C                   R.W. Brode, U.S. EPA/OAQPS/AQMG, 10/19/2009
C
C        INPUTS:  Error Message Temporary File
C
C        OUTPUTS: Total Number of Messages by Message Type
C
C        CALLED FROM:  This is A Utility Program
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: IERRLN
      INTEGER :: ICYR, ICMN, ICDY, ICHR, ICDAT, ICDAT8, ICJDY
      CHARACTER ERRTP*1, ERRCD*3, ERRMG1*50, ERRMG2*12, INPFLD3*3
      CHARACTER INPFLD12*12

C     Variable Initialization
      MODNAM = 'TERRST'
      IFTL = 0
      IWRN = 0
      INFO = 0
      ICLM = 0
      IMSG = 0
      ICYR = 0
      ICMN = 0
      ICDY = 0
      ICHR = 0
      ICDAT = 0
      ICDAT8 = 0
      ICJDY  = 0
      DNUM   = 0.0D0
      EOF = .FALSE.

C     Rewind the Temporary Error/Message File
      REWIND IERUNT

      DO WHILE (.NOT. EOF)
         READ(IERUNT,1116,END=99,ERR=9999) PATH,ERRTP,ERRCD,IERRLN,
     &                                     MODNAM,ERRMG1,ERRMG2

C        Convert error code from character string to number
         INPFLD3 = ERRCD
         CALL STONUM(INPFLD3,3,FNUM,IMIT)

         IF (ERRTP .EQ. 'E') THEN
            IFTL = IFTL + 1
         ELSE IF (ERRTP .EQ. 'W') THEN
            IWRN = IWRN + 1
         ELSE IF (ERRTP .EQ. 'I') THEN
            INFO = INFO + 1
            IF (NINT(FNUM) .EQ. 440) THEN
C ---          Determine if this calm hour is during period 
C              of data processed; convert date field from 
C              character string to number (using double precision)
               INPFLD12 = ERRMG2
               CALL STODBL(INPFLD12,12,DNUM,IMIT)
               ICDAT = IDNINT(DNUM)
               ICYR  = ICDAT/1000000
               IF (RSTINP .OR. IMSTAT(6).EQ.0 .OR.
     &            (ICDAT.GE.ISDATE .AND. ICDAT.LE.IEDATE) ) THEN
C ---             This hour is between start and end dates, 
C                 or this is a restarted model run, now
C                 determine Julian day and check IPROC array 
C                 for DAYRANGE.
                  ICMN = (ICDAT/10000) - (ICDAT/1000000)*100
                  ICDY = (ICDAT/100) - (ICDAT/10000)*100
                  IF (ICMN.GT.0 .AND. ICDY.GT.0) THEN
                     CALL JULIAN(ICYR,ICMN,ICDY,ICJDY)
                  ELSE
                     ICJDY = 0
                     CYCLE
                  END IF
                  IF (IPROC(ICJDY) .EQ. 1) THEN
C ---                Message for Calm Hour, Increment Calm Counter
                     ICLM = ICLM + 1
                  END IF
               END IF
            ELSE IF (NINT(FNUM) .EQ. 460) THEN
C ---          Determine if this missing hour is during period 
C              of data processed;  convert date field from
C              character string to number (using double precision)
               INPFLD12 = ERRMG2
               CALL STODBL(INPFLD12,12,DNUM,IMIT)
               ICDAT = IDNINT(DNUM)
               ICYR   = ICDAT/1000000
               IF (RSTINP .OR. IMSTAT(6).EQ.0 .OR.
     &            (ICDAT.GE.ISDATE .AND. ICDAT.LE.IEDATE) ) THEN
C ---             This hour is between start and end dates, 
C                 or this is a restarted model run, now
C                 determine Julian day and check IPROC array 
C                 for DAYRANGE.
                  ICMN = (ICDAT/10000) - (ICDAT/1000000)*100
                  ICDY = (ICDAT/100) - (ICDAT/10000)*100
                  IF (ICMN.GT.0 .AND. ICDY.GT.0) THEN
                     CALL JULIAN(ICYR,ICMN,ICDY,ICJDY)
                  ELSE
                     ICJDY = 0
                     CYCLE
                  END IF
                  IF (IPROC(ICJDY) .EQ. 1) THEN
C ---                Message for Missing Hour, Increment Missing Counter
                     IMSG = IMSG + 1
                  END IF
               END IF
            END IF
         END IF

         GO TO 11
 99      EOF = .TRUE.
 11      CONTINUE
      END DO

 1116 FORMAT(A2,1X,A1,A3,I8,1X,A12,2X,A50,1X,A12)

C     Use BACKSPACE To Reposition Temporary Error Message File Ahead of EOF;
C     This Is Needed in Order To Allow For Additional Message Writes
      BACKSPACE IERUNT

      GO TO 1000

C     WRITE Error Message: Error Reading Temp Error Message File
 9999 CALL ERRHDL(PATH,MODNAM,'E','510','ERRORMSG')

 1000 RETURN
      END

      SUBROUTINE SUMTBL(IOUNT)
C***********************************************************************
C                 SUMTBL Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: To Print Out The Error Summary Table
C
C        PROGRAMMER:  Jeff Wang, Roger Brode
C
C        DATE:    March 2, 1992
C
C        MODIFIED:  Increased field length for last message field from 
C                   8 to 12 to accommodate 12 character source IDs.
C                   R.W. Brode, U.S. EPA/OAQPS/AQMG, 02/28/2011
C         
C        MODIFIED:  To remove reference to legacy ISCST3 option (HE>ZI)
C                   for plume height above mixing height.
C                   Include the number of hours processed from the
C                   meteorological data file.
C                   R.W. Brode, U.S. EPA/OAQPS/AQMG, 10/19/2009
C
C        INPUTS:  Error Message Temporary File
C
C        OUTPUTS: Summary Of Errors
C
C        CALLED FROM:  This is A Utility Program
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      REAL    :: PERCENT
      INTEGER :: J, IERRLN, IOUNT
      CHARACTER ERRTP*1, ERRCD*3, ERRMG1*50, ERRMG2*12

C     Variable Initialization
      MODNAM = 'SUMTBL'

C     Write Out The Total Error Statistics
      WRITE(IOUNT,*) ' --------- Summary of Total Messages --------'
      WRITE(IOUNT,*) ' '
      WRITE(IOUNT,9014) IFTL
 9014 FORMAT(' A Total of   ',I10,' Fatal Error Message(s)')
      WRITE(IOUNT,9015) IWRN
 9015 FORMAT(' A Total of   ',I10,' Warning Message(s)')
      WRITE(IOUNT,9016) INFO
 9016 FORMAT(' A Total of   ',I10,' Informational Message(s)')
      IF (NTOTHRS .GT. 0) THEN
         WRITE(IOUNT,90171) NTOTHRS
90171    FORMAT(/,' A Total of   ',I10,' Hours Were Processed')
         WRITE(IOUNT,9017) ICLM
 9017    FORMAT(/,' A Total of   ',I10,' Calm Hours Identified')
C        Calculate percentage of missing hours, and check for > 10 percent.
         PERCENT = 100. * (FLOAT(IMSG)/FLOAT(NTOTHRS))
         WRITE(IOUNT,9018) IMSG, PERCENT
 9018    FORMAT(/,' A Total of   ',I10,' Missing Hours Identified (',
     &          F6.2,' Percent)')
         IF (NSubBGHOUR .GT. 0) THEN
            WRITE(IOUNT,99018) NSubBGHOUR
99018       FORMAT(/,' A Total of   ',I10,' Missing Hourly BACKGRND ',
     &            'Values Substituted')
         END IF
         IF (PERCENT .GT. 10.0) THEN
            WRITE(IOUNT,9019)
 9019       FORMAT(/,' CAUTION!:  Number of Missing Hours Exceeds 10 ',
     &             'Percent of Total!',/,12X,'Data May Not Be ',
     &             'Acceptable for Regulatory Applications.',/,12X,
     &             'See Section 5.3.2 of "Meteorological Monitoring ',
     &             'Guidance',/,12X,'for Regulatory Modeling ',
     &             'Applications" (EPA-454/R-99-005).')
         END IF
C ---    Output total precipipation if wet deposition algorithms are used
         IF (WDPLETE .OR. DEPOS .OR. WDEP) THEN
            WRITE(IOUNT,9020) TOTAL_PRECIP, TOTAL_PRECIP/25.4D0
 9020       FORMAT(/,' Met Data File Includes ',F10.2,' Millimeters (',
     &                 F10.3,' Inches) of Precipitation')
         END IF
      END IF
      WRITE(IOUNT,*) ' '

C     Write Out All The Fatal Error Messages
      WRITE(IOUNT,*) ' '
      WRITE(IOUNT,*) '   ******** FATAL ERROR MESSAGES ******** '
      REWIND IERUNT
      EOF = .FALSE.
      J = 0
      DO WHILE (.NOT. EOF)
         READ(IERUNT,1116,END=99,ERR=9999) PATH,ERRTP,ERRCD,IERRLN,
     &                                     MODNAM,ERRMG1,ERRMG2
         IF (ERRTP .EQ. 'E') THEN
            J = J + 1
            WRITE(IOUNT,1117) PATH,ERRTP,ERRCD,IERRLN,MODNAM(1:12),
     &                        ERRMG1,ERRMG2
         END IF
         GO TO 11
 99      EOF = .TRUE.
 11      CONTINUE
      END DO

C     If No Fatal Error Messages, Then Write 'NONE'
      IF (J .EQ. 0) THEN
         WRITE(IOUNT,*) '              ***  NONE  ***         '
         WRITE(IOUNT,*) ' '
      END IF

C     Write Out All The Warning Messages
      WRITE(IOUNT,*) ' '
      WRITE(IOUNT,*) '   ********   WARNING MESSAGES   ******** '
      REWIND IERUNT
      EOF = .FALSE.
      J = 0
      DO WHILE (.NOT. EOF)
         READ(IERUNT,1116,END=999,ERR=9999) PATH,ERRTP,ERRCD,IERRLN,
     &                                      MODNAM,ERRMG1,ERRMG2
         IF (ERRTP .EQ. 'W') THEN
            J = J + 1
            IF (.NOT. NOWARN) THEN
               IF (J .LE. 999) THEN
                  WRITE(IOUNT,1117) PATH,ERRTP,ERRCD,IERRLN,
     &                              MODNAM(1:12),ERRMG1,ERRMG2
               ELSE
                  WRITE(IOUNT,*) 'More Than 999 Warning Messages ',
     &                           'Found.  See ERRORFIL Output for',
     &                           ' the Remainder.'
                  EOF = .TRUE.
               END IF
            END IF
         END IF
         GO TO 111
 999     EOF = .TRUE.
 111     CONTINUE
      END DO

C     If No Warning Messages, Then Write 'NONE'
      IF (J .EQ. 0) THEN
         WRITE(IOUNT,*) '              ***  NONE  ***        '
         WRITE(IOUNT,*) ' '
      ELSE IF (NOWARN) THEN
         WRITE(IOUNT,*) ' ** WARNINGS SUPPRESSED BY NOWARN OPTION **'
         WRITE(IOUNT,*) ' '
      END IF

 1116 FORMAT(A2,1X,A1,A3,I8,1X,A12,2X,A50,1X,A12)
 1117 FORMAT(1X,A2,1X,A1,A3,I8,1X,A12,': ',A50,1X,A12)

C     Use BACKSPACE To Reposition Temporary Error Message File Ahead of EOF;
C     This Is Needed in Order To Allow For Additional Message Writes
      BACKSPACE IERUNT

      GO TO 1000

C     WRITE Error Message: Error Reading Temp Error Message File
 9999 CALL ERRHDL(PATH,MODNAM,'E','510','ERRORMSG')

 1000 RETURN
      END

      SUBROUTINE MSGWRT
C***********************************************************************
C                 MSGWRT Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: To Print Out The Error Summary Table
C
C        PROGRAMMER: Jeff Wang
C
C        DATE:    March 2, 1992
C
C        INPUTS:  Input Error Message File
C
C        OUTPUTS: The Error Message File
C
C        CALLED FROM:  This is A Utility Program
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: IERRLN
      CHARACTER ERRTP*1, ERRCD*3, ERRMG1*50, ERRMG2*12

C     Variable Initialization
      MODNAM = 'MSGWRT'

C     Write Out The Header Of The Message File
      WRITE(IERWRT,*) ' '
      WRITE(IERWRT,*) '   ************ Error Message List *************'
      WRITE(IERWRT,*) ' '
      WRITE(IERWRT,*) '   PW     --- Pathway                           '
      WRITE(IERWRT,*) '   Code   --- Error Type + Error Code           '
      WRITE(IERWRT,*) '   L#     --- The Line Number Where Error Occurs'
      WRITE(IERWRT,*) '   ModNam --- Module Name In Which Error Occurs '
      WRITE(IERWRT,*) '   Hints  --- Hints For The Possible Solution   '
      WRITE(IERWRT,*) '   *********************************************'
      WRITE(IERWRT,*) ' '
      WRITE(IERWRT,1114)
      WRITE(IERWRT,1115)
 1114 FORMAT('PW CODE    L#      MODNAM     ',18X,'ERROR MESSAGES',22X,
     &                                                          'HINTS')
 1115 FORMAT('-- ---- ------- ------------  ',50('-'),' ------------')
      WRITE(IERWRT,*) ' '
      REWIND IERUNT
      EOF = .FALSE.

      DO WHILE (.NOT. EOF)
         READ(IERUNT,1116,END=99,ERR=999) PATH,ERRTP,ERRCD,IERRLN,
     &                                    MODNAM,ERRMG1,ERRMG2
         WRITE(IERWRT,1117) PATH,ERRTP,ERRCD,IERRLN,
     &                      MODNAM(1:12),ERRMG1,ERRMG2
         GO TO 11
 99      EOF = .TRUE.
 11      CONTINUE
      END DO

 1116 FORMAT(A2,1X,A1,A3,I8,1X,A12,2X,A50,1X,A12)
 1117 FORMAT(A2,1X,A1,A3,I8,1X,A12,': ',A50,1X,A12)

      GO TO 1000

C     WRITE Error Message: Error Reading Temp Error Message File
 999  CALL ERRHDL(PATH,MODNAM,'E','510','ERRORMSG')

 1000 RETURN
      END

      SUBROUTINE PNPOLY (PX,PY,X,Y,N,INOUT)
C----------------------------------------------------------------------
C     Courtesy: Jay Sandhu
C               email: jsandhu@esri.com
C
C
C Please cite David H. Douglas, COLLECTED ALGORITHMS, Cambridge MA:
C Harvard Laboratory for Computer Graphics, 1974
C
C This is my reinvention buster.
C 1974 1974 1974 1974 1974 1974 1974 1974 1974 1974 1974 1974
C
C>>>PNPY
C     .................................................................
C
C        SUBROUTINE PNPOLY
C
C        PURPOSE
C           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON
C
C        USAGE
C           CALL PNPOLY (PX, PY, X, Y, N, INOUT )
C
C        DESCRIPTION OF THE PARAMETERS
C           PX      - X-COORDINATE OF POINT IN QUESTION.
C           PY      - Y-COORDINATE OF POINT IN QUESTION.
C           X       - N LONG VECTOR CONTAINING X-COORDINATES OF
C                     VERTICES OF POLYGON.
C           Y       - N LONG VECTOR CONTAINING Y-COORDINATES OF
C                     VERTICES OF POLYGON.
C           N       - NUMBER OF VERTICES IN THE POLYGON.
C           INOUT   - THE SIGNAL RETURNED:
C                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,
C                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,
C                      1 IF THE POINT IS INSIDE OF THE POLYGON.
C
C        REMARKS
C           THE VERTICES MAY BE LISTED IN CLOCKWISE OR ANTICLOCKWISE
C           ORDER.  FOR THIS SUBROUTINE A POINT IS CONSIDERED INSIDE
C           THE POLYGON IF IT IS LOCATED IN THE ENCLOSED AREA DEFINED
C           BY THE LINE FORMING THE POLYGON.
C           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING
C           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX
C           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING
C           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.
C           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.
C           PNPOLY CAN HANDLE ANY NUMBER OF VERTICES IN THE POLYGON.
C           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 6/72.
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           A VERTICAL SEMI-INFINITE LINE IS DRAWN UP FROM THE POINT
C           IN QUESTION. IF IT CROSSES THE POLYGON AN ODD NUMBER OF
C           TIMES, THE POINT IS INSIDE THE POLYGON.
C
C --- Modified to use an Internal Function for EOR rather than a
C     Statement Function, which has been identified as obsolescent 
C --- in Fortran 95.  R.W. Brode, EPA/OAQPS/AQMG, 10/19/2009
C     .................................................................
C

      IMPLICIT NONE

      INTEGER I, J, N, INOUT
      DOUBLE PRECISION X(N), Y(N), XI, YI, XJ, YJ, PX, PY
      LOGICAL IX, IY
      LOGICAL JX, JY
      LOGICAL L_EOR
      
      L_EOR = .FALSE.

      INOUT=-1

      DO I=1,N
         XI=X(I)-PX
         YI=Y(I)-PY
C        CHECK WHETHER THE POINT IN QUESTION IS AT THIS VERTEX.
         IF (XI.EQ.0.0D0 .AND. YI.EQ.0.0D0) THEN
            INOUT = 0
            EXIT
         END IF
C        J IS NEXT VERTEX NUMBER OF POLYGON.
         J=1+MOD(I,N)
         XJ=X(J)-PX
         YJ=Y(J)-PY
C        IS THIS LINE OF 0 LENGTH ?
         IF (XI.EQ.XJ .AND. YI.EQ.YJ) CYCLE
         IX=XI.GE.0.0D0
         IY=YI.GE.0.0D0
         JX=XJ.GE.0.0D0
         JY=YJ.GE.0.0D0
C        CHECK WHETHER (PX,PY) IS ON VERTICAL SIDE OF POLYGON.
         L_EOR = EOR(IY,JY)
         IF (XI.EQ.0.0D0 .AND. XJ.EQ.0.0D0 .AND. L_EOR) THEN
            INOUT = 0
            EXIT
         END IF
C        CHECK WHETHER (PX,PY) IS ON HORIZONTAL SIDE OF POLYGON.
         L_EOR = EOR(IX,JX)
         IF (YI.EQ.0.0D0 .AND. YJ.EQ.0.0D0 .AND. L_EOR) THEN
            INOUT = 0
            EXIT
         END IF
C        CHECK WHETHER BOTH ENDS OF THIS SIDE ARE COMPLETELY 1) TO RIGHT
C        OF, 2) TO LEFT OF, OR 3) BELOW (PX,PY).
         L_EOR = EOR(IX,JX)
         IF (.NOT.((IY.OR.JY).AND.L_EOR)) CYCLE
C        DOES THIS SIDE OBVIOUSLY CROSS LINE RISING VERTICALLY FROM (PX,PY)
         L_EOR = EOR(IX,JX)
         IF (.NOT.(IY.AND.JY.AND.L_EOR)) THEN
            IF ((YI*XJ-XI*YJ)/(XJ-XI) .LT. 0.0D0) THEN
               CYCLE
            ELSE IF ((YI*XJ-XI*YJ)/(XJ-XI) .EQ. 0.0D0) THEN
               INOUT = 0
               EXIT
            ELSE
               INOUT = -INOUT
            END IF
         ELSE
            INOUT = -INOUT
         END IF

      END DO

C     "EXCLUSIVE OR" Internal FUNCTION, EOR:
      CONTAINS
        LOGICAL FUNCTION EOR(IX,IY)
          LOGICAL IX, IY
          EOR = (IX.OR.IY) .AND. .NOT.(IX.AND.IY)
        END FUNCTION

      END

      SUBROUTINE ALLSETUP
C***********************************************************************
C                 ALLSETUP Module
C
C        PURPOSE: Allocate Array Storage for SETUP
C
C        PROGRAMMER: Roger Brode, PES, Inc.
C
C        DATE:    September 21, 1996
C
C        MODIFIED:  To include CHIBL, PARTCH, buoyant line parameters,
C                   and source and receptor coordinates in a rotated 
C                   system for a buoyant line source
C                   Amec Foster Wheeler, 06/30/15
C
C        MODIFIED:  To include ADSFACT, AWDSIN, and AWDCOS.
C                   R. W. Brode, MACTEC (f/k/a PES), Inc., 08/02/05
C
C        INPUTS:
C
C
C        OUTPUTS:
C
C        CALLED FROM:  MAIN
C
C        ERROR HANDLING:   Checks for error allocating arrays
C***********************************************************************
C
C     Variable Declarations
      USE MAIN1
      USE BUOYANT_LINE
      
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: IASTAT

C     Variable Initializations
      MODNAM = 'ALLSET'
      ALLOC_ERR = .FALSE.

      IF (NUMBGsects .LT. 1) NUMBGsects = 1
      IF (NUMO3sects .LT. 1) NUMO3sects = 1

      ALLOCATE  (KAVE(NAVE), CHRAVE(NAVE), CHIDEP(6,NTYP),
     &           OUTTYP(NTYP),STAT=IASTAT)
      IF (IASTAT .NE. 0) THEN
         CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
         ALLOC_ERR = .TRUE.
      END IF

      ALLOCATE  (AXS(NSRC), AYS(NSRC), AZS(NSRC), AQS(NSRC),
     &           AHS(NSRC), ATS(NSRC), AVS(NSRC), ADS(NSRC),
     &           ASYINI(NSRC), ASZINI(NSRC),
     &           AXS1(NSRC), AYS1(NSRC), AXS2(NSRC), AYS2(NSRC), 
     &           AWIDTH(NSRC), ADSFACT(NSRC), NDXSTK(NSRC),
     &           AWDSIN(NSRC), AWDCOS(NSRC), AAFV(NSRC), 
     &           INPD(NSRC), EVAL(NSRC), URBSRC(NSRC), 
     &           L_HRLYSIG(NSRC), L_FLATSRC(NSRC), 
     &           L_METHOD2(NSRC), L_WakeMessage(NSRC),
     &           IGROUP(NSRC,NGRP), SRCID(NSRC), SRCTYP(NSRC),
     &           SOPCRD(NSRC), SOGAS(NSRC), O3VARY(NO3F,NUMO3sects),
     &           GRPID(NGRP), QFLAG(NSRC), EMILBL(NTYP),
     &           OUTLBL(NTYP), PERLBL(NTYP), BACKGRND(NBF,NUMBGsects),
     &           EMIFAC(NTYP), GRP_BACK(NGRP),
     &           STAT=IASTAT)
      IF (IASTAT .NE. 0) THEN
         CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
         ALLOC_ERR = .TRUE.
         WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                   'Basic Source Arrays!'

      END IF

      IF (NSEC .GT. 0) THEN
         ALLOCATE  (ADSBH(NSEC,NSRC), ADSBW(NSEC,NSRC),
     &              ADSBL(NSEC,NSRC), ADSXADJ(NSEC,NSRC),
     &              ADSYADJ(NSEC,NSRC),
     &              STAT=IASTAT)
         IF (IASTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
            ALLOC_ERR = .TRUE.
            WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                      'Source Downwash Arrays!'

         END IF
      END IF

      IF (NQF .GT. 0) THEN
         ALLOCATE  (QFACT(NQF,NSRC), STAT=IASTAT)
         IF (IASTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
            ALLOC_ERR = .TRUE.
            WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                      'Source Emission Factor Arrays!'
         END IF
         
      END IF

      IF (NPDMAX .GT. 0) THEN
         ALLOCATE  (APDIAM(NPDMAX,NSRC), APHI(NPDMAX,NSRC),
     &              APDENS(NPDMAX,NSRC), AVGRAV(NPDMAX,NSRC), 
     &              ATSTOP(NPDMAX,NSRC),
     &              EFRAC(NPDMAX), QPART(NPDMAX),
     &              PDIAM(NPDMAX), PHI(NPDMAX), PDENS(NPDMAX),
     &              VGRAV(NPDMAX), TSTOP(NPDMAX), SCHMIDT(NPDMAX),
     &              VDEP(NPDMAX), SCF(NPDMAX), WQCOR(NPDMAX),
     &              DQCOR(NPDMAX), PSCVRT(NPDMAX), WASHOUT(NPDMAX),
     &              ECOLL(NPDMAX), finemass(NSRC),
     &              STAT=IASTAT)
         IF (IASTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
            ALLOC_ERR = .TRUE.
            WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                      'Source Particle Deposition Arrays!'
         END IF

      END IF

      ALLOCATE  (pdiff(NSRC), pdiffw(NSRC), rmolwt(NSRC), alphas(NSRC),
     &           react(NSRC), henry(NSRC), rcli(NSRC),  
     &           STAT=IASTAT)
      IF (IASTAT .NE. 0) THEN
         CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
         ALLOC_ERR = .TRUE.
         WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                   'Gas Deposition Arrays!'
      END IF

      IF (NVMAX .GT. 0) THEN
         ALLOCATE  (AXINIT(NSRC), AYINIT(NSRC), AANGLE(NSRC),
     &              AXVERT(NVMAX,NSRC), AYVERT(NVMAX,NSRC),
     &              UVERT(NVMAX), VVERT(NVMAX), VNVERT(NVMAX),
     &              WVERT(NVMAX), UASEGS(NVMAX), UBSEGS(NVMAX),
     &              XVERT(NVMAX), YVERT(NVMAX),
     &              SPA(NVMAX,2),
     &              AALPHA(NSRC), APDEFF(NSRC), AVOLUM(NSRC),
     &              RADIUS(NSRC), NVERTS(NSRC), AXCNTR(NSRC),
     &              AYCNTR(NSRC),
     &              STAT=IASTAT)
         IF (IASTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
            ALLOC_ERR = .TRUE.
            WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                      'Area Source Arrays!'
         END IF
      END IF

      IF (NBLINES .GT. 0) THEN
         ALLOCATE (BLINEPARMS(NBLINES), XS_SCS(NBLINES,129), 
     &             YS_SCS(NBLINES), XS_RCS(NBLINES,129), 
     &             YS_RCS(NBLINES,129), STAT=IASTAT)
         IF (IASTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
            ALLOC_ERR = .TRUE.
            WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                      'Buoyant Line Arrays!'
         END IF
      ENDIF         

      IF (NURB .GT. 0) THEN
         ALLOCATE  (IURBGRP(NSRC,NURB), URBID(NURB), URBNAM(NURB),
     &              URBPOP(NURB), URBZ0(NURB),
     &              ZIURB(NURB), URBWSTR(NURB), URBUSTR(NURB),
     &              URBOBULEN(NURB),
     &              GRDSWU(MXGLVL,NURB), GRDSVU(MXGLVL,NURB),
     &              GRDTGU(MXGLVL,NURB), GRDPTU(MXGLVL,NURB),
     &              L_MorningTrans(NURB),
     &              STAT=IASTAT)
         IF (IASTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
            WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                      'Urban Arrays!'
         END IF
      END IF

      IF (EVONLY) THEN
         ALLOCATE  (EV_HRQS(NSRC,NHR), EV_HRTS(NSRC,NHR),
     &              EV_HRVS(NSRC,NHR), EV_HRHS(NSRC,NHR),
     &              EV_HRSY(NSRC,NHR), EV_HRSZ(NSRC,NHR),
     &              EVAPER(NEVE), EVDATE(NEVE), EVJDAY(NEVE),
     &              IDXEV(NEVE), AXR(NEVE), AYR(NEVE), AZELEV(NEVE),
     &              AZFLAG(NEVE), AZHILL(NEVE), EVNAME(NEVE),
     &              EVGRP(NEVE), EV_OrigConc(NEVE), STAT=IASTAT)
         IF (IASTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
            ALLOC_ERR = .TRUE.
            WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                      'EVENT Processing Arrays!'
         END IF
      END IF

      IF (.NOT. EVONLY) THEN
         ALLOCATE  (AXR(NREC), AYR(NREC), AZELEV(NREC),
     &              AZFLAG(NREC), AZHILL(NREC), IREF(NREC),
     &              NETID(NREC), RECTYP(NREC),
     &              NDXARC(NREC), ARCID(NARC),
     &              NTID(NNET), NTTYP(NNET),
     &              XCOORD(IXM,NNET), YCOORD(IYM,NNET),
     &              XORIG(NNET), YORIG(NNET),
     &              NETSTA(NNET), NETEND(NNET),
     &              NUMXPT(NNET), NUMYPT(NNET), STAT=IASTAT)
         IF (IASTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
            ALLOC_ERR = .TRUE.
            WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                      'Receptor Arrays!'
         END IF
      END IF

      IF (NBLINES .GT. 0) THEN
         ALLOCATE (XR_SCS(NREC), YR_SCS(NREC), 
     &             XR_RCS(NREC), YR_RCS(NREC), STAT=IASTAT)
         IF (IASTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
            ALLOC_ERR = .TRUE.
            WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                      'Receptor Arrays for Buoyant Line!'
         END IF
      ENDIF         

      ALLOCATE  (NHIAVE(NVAL,NAVE), MAXAVE(NAVE), IMXVAL(NAVE),
     &           IDYTAB(NAVE), MAXFLE(NGRP,NAVE),
     &           IPSTFL(NGRP,NAVE), IPLTFL(NVAL,NGRP,NAVE),
     &           IANPST(NGRP), IANPLT(NGRP), INHI(NAVE),
     &           ITOXFL(NAVE), IRNKFL(NAVE), IRKVAL(NAVE),
     &           THRESH(NGRP,NAVE), TOXTHR(NAVE),
     &           IMXUNT(NGRP,NAVE), IPSUNT(NGRP,NAVE),
     &           IPSFRM(NGRP,NAVE), IPLUNT(NVAL,NGRP,NAVE),
     &           IAPUNT(NGRP), IANFRM(NGRP), IPPUNT(NGRP),
     &           ITXUNT(NAVE), IRKUNT(NAVE), IELUNT(NSRC),
     &           THRFIL(NGRP,NAVE), PSTFIL(NGRP,NAVE),
     &           PLTFIL(NVAL,NGRP,NAVE), ANNPST(NGRP),
     &           ANNPLT(NGRP), TOXFIL(NAVE), RNKFIL(NAVE),
     &           EVLFIL(NSRC), ISEAHR(NGRP), SEAHRS(NGRP),
     &           ISHUNT(NGRP), IMXDLY(NGRP), IMDUNT(NGRP), 
     &           MAXDLY(NGRP), IMXDLY_BYYR(NGRP), 
     &           IMDUNT_BYYR(NGRP), MAXDLY_BYYR(NGRP),  
     &           MAXDCONT(NGRP), IMXDCUNT(NGRP),
     &           MAXDCONT_FILE(NGRP), MXD_RANK(NGRP,2),
     &           MAXD_THRESH(NGRP), 
     &           STAT=IASTAT)
      IF (IASTAT .NE. 0) THEN
         CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
         ALLOC_ERR = .TRUE.
         WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                   'Output Option Arrays!'
      END IF

      IF (NBLINES .GT. 0) THEN
         ALLOCATE (CHIBL(NREC), PARTCH(NREC), STAT=IASTAT)
         IF (IASTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
            ALLOC_ERR = .TRUE.
            WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                      'Concentration Arrays for Buoyant Line!'
         END IF
      ENDIF         

      ALLOCATE  (IDCONC(NAVE,NPAIR), TXCONC(NAVE,NPAIR), STAT=IASTAT)
      IF (IASTAT .NE. 0) THEN
         CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
         ALLOC_ERR = .TRUE.
         WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                   'TOXXFILE Arrays!'
      END IF

      ALLOCATE  (WORKID(NSRC+1), IWRK2(NSRC,13), STAT=IASTAT)
      IF (IASTAT .NE. 0) THEN
         CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
         ALLOC_ERR = .TRUE.
         WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                   'Temporary Source Arrays!'
      END IF

      IF (.NOT. EVONLY) THEN
         ALLOCATE  (ZETMP1(NREC), ZETMP2(NREC),
     &              ZHTMP1(NREC), ZHTMP2(NREC),
     &              ZFTMP1(NREC), ZFTMP2(NREC), STAT=IASTAT)
         IF (IASTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
            ALLOC_ERR = .TRUE.
            WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                      'Temporary Receptor Arrays!'
         END IF
      END IF

      IF (PVMRM .OR. PVMRM2 .OR. OLM .OR. ARM2 .OR. ARM) THEN
         ALLOCATE (ANO2_RATIO(NSRC), CHI(NREC,NSRC,NTYP), STAT=IASTAT)
         IF (IASTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
            ALLOC_ERR = .TRUE.
            WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                      'PVMRM/OLM/ARM2/ARM CHI Array!'
         END IF
         
         IF (PVMRM .OR. PVMRM2) THEN
            ALLOCATE (HECNTR(NREC,NSRC), HECNTR3(NREC,NSRC),
     &                UEFFS(NREC,NSRC), UEFF3S(NREC,NSRC),
     &                EPSEF(NREC,NSRC), EPSEF3(NREC,NSRC),
     &                FOPTS(NREC,NSRC), PPFACT(NREC,NSRC), STAT=IASTAT)
            IF (IASTAT .NE. 0) THEN
               CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
               ALLOC_ERR = .TRUE.
               WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                         'PVMRM Source Data Arrays!'
            END IF
         END IF
         
         IF (OLM) THEN
            ALLOCATE (OLMID(NOLM), L_OLMGRP(NSRC),
     &                IGRP_OLM(NSRC,NOLM), STAT=IASTAT)
            IF (IASTAT .NE. 0) THEN
               CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
               ALLOC_ERR = .TRUE.
               WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                         'OLMGROUP Source Data Arrays!'
            END IF
         END IF
         
         IF (PSDCREDIT) THEN
            ALLOCATE (PSDSRCTYP(NSRC), PSDID(NPSD), L_PSDGRP(NSRC),
     &                IGRP_PSD(NSRC,NPSD),
     &                ABVAL(NREC,NTYP), BCVAL(NREC,NTYP), STAT=IASTAT)
            IF (IASTAT .NE. 0) THEN
               CALL ERRHDL(PATH,MODNAM,'E','409','Setup Arrays')
               ALLOC_ERR = .TRUE.
               WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                         'PSDCREDIT Source Data Arrays!'
            END IF
         END IF
      END IF

      RETURN
      END

      SUBROUTINE ALLRESULT
C***********************************************************************
C                 ALLRESULT Module
C
C        PURPOSE: Allocate Array Storage for Results
C
C        PROGRAMMER: Roger Brode, PES, Inc.
C
C        DATE:       September 21, 1996
C 
C        MODIFIED:   Modified subroutine ALLRESULT to eliminate the arrays of 
C                    profile met data (the observed data in the PROFFILE input) 
C                    by hour-of-year, level, and year for use with the MAXDCONT 
C                    option. The profile met data arrays are not needed for the 
C                    MAXDCONT option and their removal reduces the memory 
C                    requirements for that option. Additional adjustments to the
C                    array allocations for the MAXDCONT option were made based
C                    on the options used for a specific model run in order to 
C                    minimize the memory requirements for the MAXDCONT option.
C                    Subroutine ALLRESULT was also modified to improve the 
C                    accuracy of the memory storage estimates included on the
C                    first page of the AERMOD.OUT file.
C                    R.W. Brode, U.S. EPA/OAQPS/AQMG, 02/29/2012
C   
C        MODIFIED:   Added calculation of STORE, estimated memory
C                    storage requirements, to report if allocation
C                    errors occur.
C                    R.W. Brode, U.S. EPA/OAQPS/AQMG, 10/19/2009
C
C                    Changed parameter for allocating the number of
C                    high annual/period averages from NHIVAL to NHIANN.
C                    R.W. Brode, PES, Inc.,  4/3/98
C
C        INPUTS:
C
C
C        OUTPUTS:
C
C        CALLED FROM:  MAIN
C
C        ERROR HANDLING:   Checks for error allocating arrays
C***********************************************************************
C
C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: IASTAT
C     Declare Real Variables used to estimate memory requirements
      REAL    :: RSRC, RSEC, RGRP, RREC, RURB, RARC, RAVE,
     &           RHIVAL, RTYP, RMAXVAL, RNET, RXM , RYM , REVE, ROLM,
     &           RPSD, RQF, RPDMAX, RVMAX, RPAIR, RHIANN, RHIMXDLY,
     &           RYEARS

C     Variable Initializations
      MODNAM = 'ALLRESULT'
      ALLOC_ERR = .FALSE.

C     NARC was initially set to NREC prior to SETUP, now set NARC = NUMARC
      NARC = NUMARC

C --- Assign maximum value from IRKVAL array for the OU RANKFILE option
C     and from the IMXVAL array for the OU MAXTABLE option to NMXVAL
      NMXVAL = MAX( MAXVAL(IRKVAL), MAXVAL(IMXVAL) )

C     Assign array limits to REAL for calculation of STORE
      RSRC   = REAL(NSRC)
      RSEC   = REAL(NSEC)
      RGRP   = REAL(NGRP)
      RREC   = REAL(NREC)
      RURB   = REAL(NURB)
      RARC   = REAL(NARC)
      RAVE   = REAL(NAVE)
      RHIVAL = REAL(NHIVAL)
      RTYP   = REAL(NTYP)
      RMAXVAL= REAL(NMXVAL)
      RNET   = REAL(NNET)
      RXM    = REAL(IXM )
      RYM    = REAL(IYM )
      REVE   = REAL(NEVE)
      ROLM   = REAL(NOLM)
      RPSD   = REAL(NPSD)
      RQF    = REAL(NQF)
      RPDMAX = REAL(NPDMAX)
      RVMAX  = REAL(NVMAX)
      RPAIR  = REAL(NPAIR)
      RHIANN = REAL(NHIANN)
      RYEARS = REAL(NYEARS)
      RHIMXDLY = REAL(NHIMXDLY)

      STORE = 0.0
      IF (.NOT. EVONLY) THEN
C        Calculate Approximate Allocated Storage Requirements
         STORE = RSRC*(54.+RQF+5.*RSEC+5.*RPDMAX+2.*RVMAX+
     &           0.5*(RGRP+RURB))+
     &           RPDMAX*14. +
     &           RREC*(9.+RHIVAL*RGRP*RAVE*RTYP*1.75+RGRP*RAVE*RTYP+
     &                 2.*RGRP*RTYP+RGRP) +
     &           RARC*20. + RNET*(9.+RXM+RYM) +
     &           RHIVAL*(RGRP*RAVE*RTYP*3.)+
     &           RMAXVAL*(RGRP*RAVE*RTYP*2.25) +
     &           RHIANN*1.5*RGRP*RTYP +
     &           RAVE*(12.+2.*RPAIR+3.*RHIVAL*RGRP+8.*RGRP) +
     &           RGRP*11. + RTYP*38. + RVMAX*20.
         IF (SEASONHR) THEN
            STORE = STORE + ( 4.*24.*RREC*RGRP*RTYP )
         END IF
         IF (PVMRM .OR. PVMRM2 .OR. OLM) THEN
            STORE = STORE + ( RSRC + RSRC*RREC*RTYP )
            IF (PVMRM .OR. PVMRM2) THEN
               STORE = STORE + ( 5.*RSRC*RREC + RSRC )
               IF (PSDCREDIT) THEN
                  STORE = STORE + ( RSRC*(2.+RPSD) + RREC*2.*RTYP )
               END IF
            ELSE IF (OLM) THEN
               STORE = STORE + ( 0.5*ROLM*RSRC + ROLM + RSRC )
            END IF
         END IF
         IF (PM25AVE .OR. NO2AVE .OR. SO2AVE) THEN
            STORE = STORE + RREC*RGRP*(1.5+1.5*RHIMXDLY+
     &                                     1.5*RHIMXDLY*RYEARS)
         END IF
         IF (L_MAXDCONT) THEN
            STORE = STORE + RYEARS*(19.*8784.+10.*8784.+
     &                               7.*8784.*REAL(MXGLVL))
            IF (LDGAS.OR.LDPART.OR.LWGAS.OR.LWPART) THEN
               STORE = STORE + RYEARS*14.*8784.
            END IF
            IF (NSEC .GT. 0) THEN
               STORE = STORE + RYEARS*8784.*REAL(MXGLVL)
            END IF
            IF (PVMRM .OR. PVMRM2) THEN
               STORE = STORE + RYEARS*8784.*REAL(MXGLVL)
            END IF
            IF (PVMRM .OR. PVMRM2 .OR. OLM) THEN
               STORE = STORE + RYEARS*8784.
            END IF
            IF (L_BACKGRND) THEN
               STORE = STORE + RYEARS*8784.
            END IF
            IF (HOURLY) THEN
               STORE = STORE + RYEARS*REAL(2*8784)*RSRC
               IF (NPNT .GT. 0) THEN
                  STORE = STORE + RYEARS*REAL(2*8784)*RSRC
               END IF
               IF (NVOL .GT. 0 .OR. NVMAX .GT. 0) THEN
                  STORE = STORE + RYEARS*REAL(2*8784)*RSRC
               END IF
            END IF            
            IF (NURB .GT. 0) THEN
               STORE = STORE + REAL(4*8784*MXGLVL*NYEARS*NURB) +
     &                         REAL(4*8784*MXGLVL*NYEARS) +
     &                         REAL(5*8784*NYEARS*NURB)
            END IF
         END IF
         STORE = STORE*8./1.048576E6 + 3.5
      END IF

      ALLOCATE  (HRVAL(NUMTYP), AERVAL(NUMTYP), PRMVAL(NUMTYP),
     &           STAT=IASTAT)

C --- Allocate BACKAVE and BACKANN arrays to store total background contribution to GRP value
      ALLOCATE  (BACKAVE(NUMGRP), BACKANN(NUMGRP), STAT=IASTAT)

      IF (NARC .GT. 0) THEN
         ALLOCATE  (ARCMAX(NARC), QMAX(NARC), DXMAX(NARC), UMAX(NARC),
     &              SVMAX(NARC), SWMAX(NARC), SYMAX(NARC), SY3MX(NARC),
     &              U3MAX(NARC), HEMAX(NARC), ARCCL(NARC), SZMAX(NARC),
     &              CHIDMW(NARC), CHINMW(NARC), CHI3MW(NARC),
     &              CHIDML(NARC), CHINML(NARC), CHI3ML(NARC),
     &              HSBLMX(NARC),
     &              STAT=IASTAT)
      END IF

      IF (.NOT. EVONLY) THEN
         ALLOCATE  (AVEVAL(NUMREC,NUMGRP,NUMAVE,NUMTYP),
     &              HIVALU(NUMREC,NHIVAL,NUMGRP,NUMAVE,NUMTYP),
     &              HMAX(NHIVAL,NUMGRP,NUMAVE,NUMTYP),
     &              HMLOC(NHIVAL,NUMGRP,NUMAVE,NUMTYP),
     &              HMDATE(NHIVAL,NUMGRP,NUMAVE,NUMTYP),
     &              NHIDAT(NUMREC,NHIVAL,NUMGRP,NUMAVE,NUMTYP),
     &              STAT=IASTAT)
         IF (IASTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Result Array')
            ALLOC_ERR = .TRUE.
            WRITE(IOUNIT,*) ' '
            WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                      'Short Term Average Results Arrays!'
         END IF

         IF (PERIOD .OR. ANNUAL) THEN
            ALLOCATE  (ANNVAL(NUMREC,NUMGRP,NUMTYP),
     &                 AMXVAL(NHIANN,NUMGRP,NUMTYP),
     &                 IMXLOC(NHIANN,NUMGRP,NUMTYP),
     &                 STAT=IASTAT)
            IF (IASTAT .NE. 0) THEN
               CALL ERRHDL(PATH,MODNAM,'E','409','Result Array')
               ALLOC_ERR = .TRUE.
               WRITE(IOUNIT,*) ' '
               WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                         'Long Term Average Results Arrays!'
            END IF
         END IF

         ALLOCATE  (RMXVAL(NMXVAL,NUMGRP,NUMAVE,NUMTYP),
     &              MXDATE(NMXVAL,NUMGRP,NUMAVE,NUMTYP),
     &              MXLOCA(NMXVAL,NUMGRP,NUMAVE,NUMTYP),
     &              NUMHRS(NUMAVE), NUMCLM(NUMAVE), NUMMSG(NUMAVE),
     &              STAT=IASTAT)
         IF (IASTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Result Array')
            ALLOC_ERR = .TRUE.
            WRITE(IOUNIT,*) ' '
            WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                      'Overall Maximum Results Arrays!'
         END IF

         IF (SEASONHR) THEN
            ALLOCATE (SHVALS(NUMREC,NUMGRP,4,24,NUMTYP),
     &                STAT=IASTAT)
            IF (IASTAT .NE. 0) THEN
               CALL ERRHDL(PATH,MODNAM,'E','409','Result Array')
               ALLOC_ERR = .TRUE.
               WRITE(IOUNIT,*) ' '
               WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                         'SEASONHR Results Array!'
            END IF

            IF (L_BACKGRND) THEN
               ALLOCATE (BACKSEASHR(NUMGRP,4,24),STAT=IASTAT)
               IF (IASTAT .NE. 0) THEN
                  CALL ERRHDL(PATH,MODNAM,'E','409','Result Array')
                  ALLOC_ERR = .TRUE.
                  WRITE(IOUNIT,*) ' '
                  WRITE(IOUNIT,*) '  Error Occurred During Allocation ',
     &                            'of SEASONHR BACKGRND Array!'
               END IF
            END IF
         END IF

         IF (PM25AVE .OR. NO2AVE .OR. SO2AVE) THEN
            NHIMXDLY = MAX( NHIMXDLY, NHIVAL )
            ALLOCATE (MXDVAL(NUMREC,NUMGRP),
     &                HIMXDLY(NUMREC,NUMGRP,NHIMXDLY),
     &                HIMXDLY_BYYR(NUMREC,NUMGRP,NHIMXDLY,NYEARS),
     &                IMXDHR(NUMREC,NUMGRP),
     &                NHIDATMXD(NUMREC,NUMGRP,NHIMXDLY),
     &                NHIDATMXD_BYYR(NUMREC,NUMGRP,NHIMXDLY,NYEARS),
     &                STAT=IASTAT)
            IF (IASTAT .NE. 0) THEN
               CALL ERRHDL(PATH,MODNAM,'E','409','Result Array')
               ALLOC_ERR = .TRUE.
               WRITE(IOUNIT,*) ' '
               WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                         'MAXDAILY Results Arrays!'
           END IF
         END IF

         ALLOCATE  (HCLMSG(NUMREC,NHIVAL,NUMGRP,NUMAVE,NUMTYP),
     &              MCLMSG(NMXVAL,NUMGRP,NUMAVE,NUMTYP),
     &              HMCLM(NHIVAL,NUMGRP,NUMAVE,NUMTYP),STAT=IASTAT)
         IF (IASTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Result Array')
            ALLOC_ERR = .TRUE.
            WRITE(IOUNIT,*) ' '
            WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                      'High Value Result Flag Arrays!'
         END IF

         IF (ANNUAL) THEN
            ALLOCATE  (SUMANN(NUMREC,NUMGRP,NUMTYP),
     &                 STAT=IASTAT)
            IF (IASTAT .NE. 0) THEN
               CALL ERRHDL(PATH,MODNAM,'E','409','Result Array')
               ALLOC_ERR = .TRUE.
               WRITE(IOUNIT,*) ' '
               WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                         'ANNUAL Results Arrays!'
            END IF
         END IF

         IF (PM25AVE .OR. NO2AVE .OR. SO2AVE) THEN
            ALLOCATE  (SUMHNH(NUMREC,NUMGRP,NHIVAL),
     &                 MXPMVAL(NMXPM,NUMGRP,NHIVAL), 
     &                 MXPMLOC(NMXPM,NUMGRP,NHIVAL),
     &                 STAT=IASTAT)
            IF (IASTAT .NE. 0) THEN
               CALL ERRHDL(PATH,MODNAM,'E','409','Result Array')
               ALLOC_ERR = .TRUE.
               WRITE(IOUNIT,*) ' '
               IF (PM25AVE) THEN
                  WRITE(IOUNIT,*) '  Error Occurred During Allocation',
     &                            ' of PM-2.5 24-hr Results Arrays!'
               ELSE IF (NO2AVE) THEN
                  WRITE(IOUNIT,*) '  Error Occurred During Allocation',
     &                            ' of NO2 1-hr Results Arrays!'
               ELSE IF (SO2AVE) THEN
                  WRITE(IOUNIT,*) '  Error Occurred During Allocation',
     &                            ' of SO2 1-hr Results Arrays!'
               END IF
            END IF
         END IF

      END IF

      IF (EVONLY) THEN
         ALLOCATE  (ASFCHF(NHR,1), AUREF(NHR,1), 
     &              AUREFHT(NHR,1), ATA(NHR,1), 
     &              ATREFHT(NHR,1), AWDREF(NHR,1), 
     &              AUSTAR(NHR,1), AWSTAR(NHR,1),
     &              AZICONV(NHR,1), AZIMECH(NHR,1), 
     &              AOBULEN(NHR,1), AVPTGZI(NHR,1),
     &              ASFCZ0(NHR,1), ABOWEN(NHR,1),
     &              AALBEDO(NHR,1), AWNEW(NHR,1), 
     &              AWOLD(NHR,1), AESTA(NHR,1), 
     &              AKST(NHR,1), ABLTA(NHR,1),
     &              AF2(NHR,1), APREC1(NHR,1), AQSW(NHR,1),
     &              APRATE(NHR,1), ARH(NHR,1), ASFCP(NHR,1),
     &              APREC2(NHR,1), IAPCODE(NHR,1), NACLOUD(NHR,1),
     &              ACLMHR(NHR,1), AMSGHR(NHR,1),
     &              AUNSTAB(NHR,1), ASTABLE(NHR,1),
     &              AURBSTAB(NHR,1),
     &              ANPLVLS(NHR,1), ANTGLVL(NHR,1), 
     &              AO3CONC(NHR,1), ABGCONC(NHR,1),
     &              AAQS(NHR,1,NSRC), AAHS(NHR,1,NSRC),
     &              AAVS(NHR,1,NSRC), AATS(NHR,1,NSRC),
     &              AASYINI(NHR,1,NSRC), AASZINI(NHR,1,NSRC),
     &              AIFLAG(NHR,MXPLVL,1),
     &              APFLHT(NHR,MXPLVL,1), APFLWD(NHR,MXPLVL,1),
     &              APFLWS(NHR,MXPLVL,1), APFLTA(NHR,MXPLVL,1),
     &              APFLSA(NHR,MXPLVL,1), APFLSW(NHR,MXPLVL,1),
     &              APFLSV(NHR,MXPLVL,1), APFLTG(NHR,MXPLVL,1),
     &              APFLTGZ(NHR,MXPLVL,1),
     &              EV_AVEVAL(NSRC), HRVALS(NHR,NSRC), 
     &              GRPVAL(NGRP,NHR), BACKHR(NGRP,NHR),
     &              GRPAVE(NGRP),
     &              STAT=IASTAT)
         IF (IASTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Result Array')
            ALLOC_ERR = .TRUE.
            WRITE(IOUNIT,*) ' '
            WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                      'Event Arrays!'
         END IF
      ELSE IF (L_MAXDCONT) THEN
         ALLOCATE (ASFCHF(8784,NYEARS), AUREF(8784,NYEARS), 
     &             AUREFHT(8784,NYEARS), ATA(8784,NYEARS), 
     &             ATREFHT(8784,NYEARS), AWDREF(8784,NYEARS), 
     &             AUSTAR(8784,NYEARS), AWSTAR(8784,NYEARS),
     &             AZICONV(8784,NYEARS), AZIMECH(8784,NYEARS), 
     &             AOBULEN(8784,NYEARS), AVPTGZI(8784,NYEARS),
     &             ASFCZ0(8784,NYEARS), 
     &             AKST(8784,NYEARS), ABLTA(8784,NYEARS),
     &             ACLMHR(8784,NYEARS), AMSGHR(8784,NYEARS),
     &             AUNSTAB(8784,NYEARS), ASTABLE(8784,NYEARS),
     &             AURBSTAB(8784,NYEARS), 
     &             ANTGLVL(8784,NYEARS), 
     &             AGRIDHT(8784,MXGLVL,NYEARS), 
     &             AGRIDWD(8784,MXGLVL,NYEARS),
     &             AGRIDWS(8784,MXGLVL,NYEARS), 
     &             AGRIDSW(8784,MXGLVL,NYEARS),
     &             AGRIDSV(8784,MXGLVL,NYEARS),
     &             AGRIDTG(8784,MXGLVL,NYEARS),
     &             AGRIDPT(8784,MXGLVL,NYEARS),
     &             AUATZI(8784,NYEARS),
     &             ASVATZI(8784,NYEARS),
     &             ASWATZI(8784,NYEARS),
     &             AUAVG(8784,NYEARS),
     &             ASVAVG(8784,NYEARS),
     &             ASWAVG(8784,NYEARS),
     &             APTATZI(8784,NYEARS),
     &             ANDX4ZI(8784,NYEARS),
     &             ARURUSTR(8784,NYEARS),
     &             ARUROBULEN(8784,NYEARS),
     &             STAT=IASTAT)
         IF (IASTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Result Array')
            ALLOC_ERR = .TRUE.
            WRITE(IOUNIT,*) ' '
            WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                      'MAXDCONT Arrays!'
         END IF

C ---    Allocate arrays for NO2 options, PVMRM or OLM
         IF (PVMRM .OR. PVMRM2) THEN
            ALLOCATE( AGRIDEPS(8784,MXGLVL,NYEARS),
     &                AO3CONC(8784,NYEARS), STAT=IASTAT )
         END IF
         IF (OLM) THEN
            ALLOCATE( AO3CONC(8784,NYEARS), STAT=IASTAT )
         END IF

         IF (L_BACKGRND) THEN
            ALLOCATE( ABGCONC(8784,NYEARS),STAT=IASTAT )
         END IF
         IF (NSEC .GT. 0) THEN
            ALLOCATE( AGRIDRHO(8784,MXGLVL,NYEARS),STAT=IASTAT )
         END IF
         IF (LDGAS .OR. LDPART .OR. LWPART .OR. LWGAS) THEN
            ALLOCATE( ABOWEN(8784,NYEARS),
     &           AALBEDO(8784,NYEARS), AWNEW(8784,NYEARS), 
     &           AWOLD(8784,NYEARS), AESTA(8784,NYEARS), 
     &           AF2(8784,NYEARS), APREC1(8784,NYEARS), 
     &           APREC2(8784,NYEARS), APRATE(8784,NYEARS), 
     &           ARH(8784,NYEARS), ASFCP(8784,NYEARS),AQSW(8784,NYEARS),
     &           IAPCODE(8784,NYEARS), NACLOUD(8784,NYEARS),STAT=IASTAT)
         END IF

         IF (IASTAT .NE. 0) THEN
            CALL ERRHDL(PATH,MODNAM,'E','409','Result Array')
            ALLOC_ERR = .TRUE.
            WRITE(IOUNIT,*) ' '
            WRITE(IOUNIT,*) '  Error Occurred During Allocation of ',
     &                      'MAXDCONT Arrays!'
         END IF

         IF (HOURLY) THEN
            ALLOCATE(
     &           AAQS(8784,NYEARS,NSRC), AAHS(8784,NYEARS,NSRC),
     &           STAT=IASTAT)
            IF (NPNT .GT. 0) THEN
               ALLOCATE(
     &           AAVS(8784,NYEARS,NSRC), AATS(8784,NYEARS,NSRC),
     &           STAT=IASTAT)
            END IF
            IF (NVOL .GT. 0 .OR. NVMAX .GT. 0) THEN
               ALLOCATE(
     &           AASYINI(8784,NYEARS,NSRC), AASZINI(8784,NYEARS,NSRC),
     &           STAT=IASTAT)
            END IF
            IF (IASTAT .NE. 0) THEN
               CALL ERRHDL(PATH,MODNAM,'E','409','Result Array')
               ALLOC_ERR = .TRUE.
               WRITE(IOUNIT,*) ' '
               WRITE(IOUNIT,*) '  Error Occurred During Allocation',
     &                         ' of MAXDCONT Arrays for HOUREMIS!'
            END IF
         END IF

         IF (NURB .GT. 0) THEN
            ALLOCATE(
     &           AGRDSWR(8784,MXGLVL,NYEARS), 
     &           AGRDSVR(8784,MXGLVL,NYEARS),
     &           AGRDTGR(8784,MXGLVL,NYEARS),
     &           AGRDPTR(8784,MXGLVL,NYEARS),
     &           AGRDSWU(8784,MXGLVL,NYEARS,NURB), 
     &           AGRDSVU(8784,MXGLVL,NYEARS,NURB),
     &           AGRDTGU(8784,MXGLVL,NYEARS,NURB), 
     &           AGRDPTU(8784,MXGLVL,NYEARS,NURB),
     &           AZIURB(8784,NYEARS,NURB),
     &           AURBWSTR(8784,NYEARS,NURB),
     &           AURBUSTR(8784,NYEARS,NURB),
     &           AURBOBULEN(8784,NYEARS,NURB),
     &           AL_MorningTrans(8784,NYEARS,NURB),
     &           STAT=IASTAT)
            IF (IASTAT .NE. 0) THEN
               CALL ERRHDL(PATH,MODNAM,'E','409','Result Array')
               ALLOC_ERR = .TRUE.
               WRITE(IOUNIT,*) ' '
               WRITE(IOUNIT,*) '  Error Occurred During Allocation',
     &                         ' of MAXDCONT Arrays for URBANOPT!'
            END IF

         END IF
      END IF

      RETURN
      END

CCLC  Note that SUBROUTINE GETCOM is no longer called to retrieve
CCLC  the command line arguments.  The input and output filenames
CCLC  hardwired by the model as AERMOD.INP and AERMOD.OUT, respectively.
CCL      SUBROUTINE GETCOM (MODEL,LENGTH,INPFIL,OUTFIL)
CCLC***********************************************************************
CCLC
CCLC        GETCOM Module of the AMS/EPA Regulatory Model - AERMOD
CCLC
CCLC        PURPOSE: Controls Retrieving Input and Output File Names From
CCLC                 the Command Line for PCs
CCLC
CCLC        PROGRAMMER: Roger Brode
CCLC
CCLC        DATE:    March 2, 1992
CCLC
CCLC        MODIFIED:   To use ILEN_FLD (passed in as LENGTH) to define
CCLC                    the length of the INPFIL and OUTFIL variables,
CCLC                    and to specify length of the command line as
CCLC                    a PARAMETER, initially set to 150.  Also set up
CCLC                    conditional compilation statements (commented out)
CCLC                    to facilitate compilation by Compaq Visual Fortran.
CCLC                    R.W. Brode, PES, Inc. - 12/2/98
CCLC
CCLC        MODIFIED:   Jayant Hardikar, PES, Inc.
CCLC                    - Length of command line for Lahey version changed
CCLC                      from 80 to 120 characters - 4/19/93
CCLC                    - Adapted for DEPMET/PMERGE - 7/29/94
CCLC
CCLC        INPUTS:  Command Line
CCLC
CCLC        OUTPUTS: Input Runstream File Name
CCLC                 Output Print File Name
CCLC
CCLC        CALLED FROM:   MAIN
CCLC***********************************************************************
CCLC
CCLC     Variable Declarations
CCLC     For compilation with Compaq Visual Fortran Compiler, delete the string
CCLC     'CCVF' from columns 1-4 in this subroutine (using a null replacement).
CCLC     This will allow the Compaq compiler to conditionally compile the
CCLC     appropriate code for retrieving the command line arguments.
CCL!DEC$ DEFINE CVF
CCL!DEC$ IF DEFINED (CVF)
CCL      USE DFLIB
CCL!DEC$ ENDIF
CCL      IMPLICIT NONE
CCL
CCL      INTEGER LENGTH
CCL      CHARACTER (LEN=LENGTH) :: INPFIL, OUTFIL
CCL      CHARACTER (LEN=8)      :: MODEL
CCL!DEC$ IF DEFINED (CVF)
CCLC     Declare 2-Byte Integer for Field Number of Command Line Argument
CCL      INTEGER*2 IARG, IFCNT, ISTAT
CCL!DEC$ ELSEIF DEFINED (LAHEY)
CCLC     Declare the COMLIN Variable to Hold Contents of Command Line for Lahey
CCL      INTEGER , PARAMETER :: LENCL = 150
CCL      CHARACTER (LEN=LENCL) :: COMLIN
CCL      INTEGER LOCB(LENCL), LOCE(LENCL), I, IFCNT
CCL      LOGICAL INFLD
CCL
CCL      COMLIN = ' '
CCL!DEC$ ENDIF
CCL
CCL!DEC$ IF DEFINED (CVF)
CCLC************************************************************CVF START
CCLC     Use Microsoft/DEC Functions NARGS and GETARG To Retrieve
CCLC     Contents of Command Line
CCL      IFCNT = NARGS()
CCLC     IFCNT Is The Number Of Arguments on Command Line Including Program
CCL      IF (IFCNT .NE. 3) THEN
CCLC        Error on Command Line.  Write Error Message and STOP
CCL         WRITE(*,660) MODEL
CCL         STOP
CCL      ELSE
CCLC        Retrieve First Argument as Input File Name
CCL         IARG = 1
CCL         CALL GETARG(IARG,INPFIL,ISTAT)
CCLC        Retrieve Second Argument as Output File Name
CCL         IARG = 2
CCL         CALL GETARG(IARG,OUTFIL,ISTAT)
CCL      END IF
CCLC************************************************************CVF STOP
CCL
CCL!DEC$ ELSEIF DEFINED (LAHEY)
CCLC************************************************************LAHEY START
CCLC     Use Lahey Function GETCL To Retrieve Contents of Command Line.
CCLC     Retrieve Input and Output File Names From the COMLIN Variable.
CCL      CALL GETCL(COMLIN)
CCL      INFLD = .FALSE.
CCL      IFCNT = 0
CCL      DO I = 1, LENCL
CCL         IF (.NOT.INFLD .AND. COMLIN(I:I) .NE. ' ') THEN
CCL            INFLD = .TRUE.
CCL            IFCNT = IFCNT + 1
CCL            LOCB(IFCNT) = I
CCL         ELSE IF (INFLD .AND. COMLIN(I:I) .EQ. ' ') THEN
CCL            INFLD = .FALSE.
CCL            LOCE(IFCNT) = I - 1
CCL         END IF
CCL      END DO
CCL      IF (IFCNT .NE. 2) THEN
CCLC        Error on Command Line.  Write Error Message and STOP
CCL         WRITE(*,660) MODEL
CCL         STOP
CCL      END IF
CCL      INPFIL = COMLIN(LOCB(1):LOCE(1))
CCL      OUTFIL = COMLIN(LOCB(2):LOCE(2))
CCLC************************************************************LAHEY STOP
CCL
CCL!DEC$ ENDIF
CCL
CCL  660 FORMAT (' COMMAND LINE ERROR: ',A8,' input_file output_file')
CCL
CCL      RETURN
CCL      END


      SUBROUTINE DATIME ( DCALL, TCALL )
C***********************************************************************
C                 DATIME Module
C
C        PURPOSE: Obtain the system date and time
C
C        PROGRAMMER: Jim Paumier, PES, Inc.
C
C        DATE:    April 15, 1994
C
C        MODIFIED:   Uses Fortran 90 DATE_AND_TIME routine.
C                    R.W. Brode, PES, 8/14/98
C
C        INPUTS:  none
C
C        OUTPUTS: Date and time in character format
C
C        CALLED FROM:  RUNTIME
C***********************************************************************
C
C     Variable Declarations
      IMPLICIT NONE

      CHARACTER DCALL*8, TCALL*8
      CHARACTER CDATE*8, CTIME*10, CZONE*5
      INTEGER :: IDATETIME(8)
      INTEGER :: IPTYR, IPTMON, IPTDAY, IPTHR, IPTMIN, IPTSEC

      DCALL = ' '
      TCALL = ' '

C     Call Fortran 90 date and time routine
      CALL DATE_AND_TIME (CDATE, CTIME, CZONE, IDATETIME)

C     Convert year to two digits and store array variables
      IPTYR  = IDATETIME(1) - 100 * INT(IDATETIME(1)/100)
      IPTMON = IDATETIME(2)
      IPTDAY = IDATETIME(3)
      IPTHR  = IDATETIME(5)
      IPTMIN = IDATETIME(6)
      IPTSEC = IDATETIME(7)

C     Write Date and Time to Character Variables, DCALL & TCALL
      WRITE(DCALL, '(2(I2.2,"/"),I2.2)' ) IPTMON, IPTDAY, IPTYR
      WRITE(TCALL, '(2(I2.2,":"),I2.2)' ) IPTHR, IPTMIN, IPTSEC

      RETURN
      END

      SUBROUTINE FILOPN
C***********************************************************************
C                 FILOPN Module
C
C        PURPOSE: Obtain the system date and time
C
C        PROGRAMMER: Roger Brode, PES, Inc.
C
C        DATE:    December 6, 1994
C
C        MODIFIED:   Remove non-standard option for 
C                    CARRIAGECONTROL='Fortran' to control
C                    page feed in aermod.out file.  ASCII form
C                    feed character is used in subroutine HEADER
C                    to insert page feed instead of using Fortan
C                    carriage control.
C                    R.W. Brode, USEPA/OAQPS/AQMG, October 19, 2009
C
C        INPUTS:  Input filename, INPFIL
C                 Output filename, OUTFIL
C
C        OUTPUTS: Openned files
C
C        CALLED FROM:  HEADER
C
C        ERROR HANDLING:   Checks errors openning files
C***********************************************************************
C
C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      MODNAM = 'FILOPN'
      
C     OPEN Input Runstream File, Unit INUNIT=7
      DUMMY = 'RUN-STRM'
      OPEN (UNIT=INUNIT,FILE='aermod.inp',ACTION='READ',ERR=99,
     &      STATUS='OLD')

C     OPEN Print Output File, Unit IOUNIT=8
      DUMMY = 'OUTPUT'
      OPEN (UNIT=IOUNIT,FILE='aermod.out',
     &      ERR=99,STATUS='REPLACE')

C     Write Out Update to the Screen
      WRITE(*,909)
 909  FORMAT('+','Now Processing SETUP Information')

      GO TO 1000

C     WRITE Error Message:  Error Opening File
 99   CALL ERRHDL('  ',MODNAM,'E','500',DUMMY)

C     Check for Error Opening Runstream File and STOP
      IF (DUMMY .EQ. 'RUN-STRM') THEN
         WRITE(*,919)
 919     FORMAT('+','Error Opening Runstream Input File!  Aborting.')
         STOP
      END IF

 1000 CONTINUE

      RETURN
      END

      SUBROUTINE HEADER(IOUNT)
C***********************************************************************
C                 HEADER Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Control Page Feed and Header Information for
C                 Printed File Output
C
C        PROGRAMMER: Roger Brode, Jeff Wang
C
C        DATE:    September 28, 1993
C
C        MODIFIED:   Use ASCII form feed character [ACHAR(12)] for 
C                    page feed in 'aermap.out' file rather then
C                    CARRIAGECONTROL='Fortran', which is not a
C                    standard Fortran option.
C                    Include adjustments to header format for large
C                    page numbers.
C                    Include output file unit argument to support
C                    output to main 'aermod.out' file and to the
C                    optional SUMMFILE.
C                    R.W. Brode, USEPA/OAQPS/AQMG, October 19, 2009
C
C        MODIFIED:   Replace DEPLETE parameter for plume depletion option
C                    with DDPLETE and WDPLETE in the list of model options
C                    for Wet & Dry depletion.
C                    D. Strimaitis, SRC - 11/8/93
C
C        MODIFIED:   Header modified for draft version of model with new
C                    area source and deposition algorithms - 9/28/93
C
C        MODIFIED:   To add DEPLETE parameter for plume depletion option
C                    to the list of model options
C                    D. Strimaitis, SRC - 2/15/93
C
C        INPUTS:  Page Number from COMMON
C
C        OUTPUTS: Page Feed and Header
C
C        CALLED FROM:  (This Is An Utility Program)
C***********************************************************************
C
C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: I, J, IOUNT, ILEN
      CHARACTER FFEED*1

C     Variable Initializations
      MODNAM = 'HEADER'

C*    FFEED is ASCII form-feed character
      J = 12
      FFEED  = ACHAR(J)

C     Increment Page Number Counter
      IF (IOUNT .EQ. IOUNIT) THEN
         IPAGE = IPAGE + 1
      ELSE IF (IOUNT .EQ. ISUMUNT) THEN
         IPGSUM = IPGSUM + 1
      END IF

C     Write Header to Printed Output File
      WRITE(IOUNT,9028) FFEED, VERSN, TITLE1(1:68), RUNDAT
      IF (IOUNT .EQ. IOUNIT) THEN
C        Adjust format statement based on page number
         IF (IPAGE .LE. 999) THEN
            WRITE(IOUNT,9029) C_METVER, TITLE2(1:68), RUNTIM, IPAGE
         ELSE IF (IPAGE .LE. 99999) THEN
            WRITE(IOUNT,90291) C_METVER, TITLE2(1:68), RUNTIM, IPAGE
         ELSE IF (IPAGE .LE. 9999999) THEN
            WRITE(IOUNT,90292) C_METVER, TITLE2(1:68), RUNTIM, IPAGE
         ELSE
            WRITE(IOUNT,90292) C_METVER, TITLE2(1:68), RUNTIM, 
     &                                              MIN(IPAGE,99999999)
         END IF
      ELSE IF (IOUNT .EQ. ISUMUNT) THEN
         WRITE(IOUNT,9029) C_METVER, TITLE2(1:68), RUNTIM, IPGSUM
      END IF
      ILEN = LEN_TRIM( MODOPS_String )
      IF (ILEN .LE. 110) THEN
         WRITE(IOUNT,9030) MODOPS_String(1:LEN_TRIM(MODOPS_String))
      ELSE
         WRITE(IOUNT,9030) MODOPS_String(1:LEN_TRIM(MODOPS_String))
         WRITE(IOUNT,9040) MODOPS_String(LEN_TRIM(MODOPS_String)+1:)
      END IF

 9028 FORMAT(A1,1X,'*** AERMOD - VERSION ',A6,' ***',3X,'*** ',A68,
     &       ' ***',8X,A8)
 9029 FORMAT(1X,'*** AERMET - VERSION ',A6,' ***',3X,'*** ',A68,' ***',
     &       8X,A8,/T120,'PAGE',I4)
90291 FORMAT(1X,'*** AERMET - VERSION ',A6,' ***',3X,'*** ',A68,' ***',
     &       8X,A8,/T118,'PAGE',I6)
90292 FORMAT(1X,'*** AERMET - VERSION ',A6,' ***',3X,'*** ',A68,' ***',
     &       8X,A8,/T116,'PAGE',I8)
 9030 FORMAT(1X,'**MODELOPTs: ',1X,A:)
 9040 FORMAT(4X,A:)

      RETURN
      END

      SUBROUTINE DCDLAT ()
C***********************************************************************
C            DCDLAT Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: To decode the hemisphere and latitude from
C                 the character variable ALAT (record 1 in scalar file)
C
C        PROGRAMMER: Jim Paumier, PES, Inc.
C
C        DATE:       September 30, 1993
C
C        INPUTS:  ALAT, the character variable latitude from AERMET
C
C        ASSUMPTIONS:  The first field in the first record of the
C                      scalar input file contains the latitude
C
C        OUTPUTS: Hemisphere (NORTH or SOUTH), latitude and sign (TSIGN)
C                 for turning of wind with height
C
C        CALLED FROM:  HRLOOP
C***********************************************************************

C---- Variable declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER      NORS, SORN

C---- Data initialization
      MODNAM = 'DCDLAT'
      PATH   = 'ME'

C---- Determine if the letter 'N' or 'n' is in the latitude field
      NORS = INDEX(ALAT,'N') + INDEX(ALAT,'n')

      IF( NORS .NE. 0 )THEN

C        The latitude is in the northern hemisphere; decode the latitude

         TSIGN = 1.0D0
         READ( ALAT, '(F9.1)',ERR=1000 ) XLAT

C        Write a message if the latitude is too far north

         IF( XLAT .GT. 90.0D0  .OR.  XLAT .LT. 0.0D0 )THEN
C           Write a warning to the user - latitude out-of-range
            CALL ERRHDL(PATH,MODNAM,'W','381',ALAT)
         END IF

      ELSE

C        The latitude may be in the southern hemisphere

         SORN = INDEX(ALAT,'S') + INDEX(ALAT,'s')
         IF( SORN .NE. 0 )THEN
            TSIGN = -1.0D0
            READ( ALAT, '(F9.1)',ERR=1000 ) XLAT

            IF( XLAT .GT. 90.0D0  .OR.  XLAT .LT. 0.0D0 )THEN
C              Write a warning to the user - latitude out-of-range
               CALL ERRHDL(PATH,MODNAM,'W','381',ALAT)
            END IF


         ELSE
C           Write a warning to the user - error decoding the latitude
            CALL ERRHDL(PATH,MODNAM,'W','382',ALAT)

         END IF

      END IF

      GO TO 999

1000  CONTINUE
C     Write a warning to the user - error decoding the latitude
      CALL ERRHDL(PATH,MODNAM,'W','382',ALAT)

999   RETURN
      END

      SUBROUTINE PRESET
C***********************************************************************
C                 PRESET Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Preprocesses SETUP Information to Determine Data
C                 Storage Requirements
C
C        PROGRAMMER: Roger Brode
C
C        DATE:    September 24, 1996
C
C        MODIFIED:   Added new NUMYEARS option to specify the number
C                    of years of met data for allocating arrays for 
C                    the MAXDCONT option. The default number of years
C                    is still five (5) years.
C                    Also modified to allow for the use of URBANSRC ALL 
C                    on the SO pathway to indicate that all sources are 
C                    to be treated as URBAN sources. This option assumes 
C                    that only one (1) urban area has been defined using 
C                    the CO URBANOPT keyword.
C                    R.W. Brode, U.S. EPA/OAQPS/AQMG, 02/29/2012
C
C        MODIFIED:   Added calculation of STORE, estimated memory
C                    storage requirements, to report if allocation
C                    errors occur.
C                    R.W. Brode, U.S. EPA/OAQPS/AQMG, 10/19/2009
C
C        MODIFIED:   To include new options incorporated in version
C                    dated 06341.
C                    R.W. Brode, U.S. EPA/OAQPS/AQMG - 12/07/2006
C
C        MODIFIED:   To check for NO ECHO in the input file.
C                    R.W. Brode, PES, Inc. - 12/2/98
C
C        INPUTS:  Input Runstream File
C
C        OUTPUTS: Array Sizes
C
C        CALLED FROM:   MAIN
C***********************************************************************
C
C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: I, J, K, ISPRD, IEPRD
C     Declare Real Variables used to estimate memory requirements
      REAL    :: RSRC, RSEC, RGRP, RREC, RURB, RARC, RAVE,
     &           RVAL, RTYP, RMAX, RNET, RXM , RYM , REVE, ROLM, RPSD,
     &           RQF, RPDMAX, RVMAX, RPAIR, RHIANN

      LOGICAL NOPATH, NOKEY, L_O3data
      CHARACTER RDFRM*20
      CHARACTER LPRD*8, HPRD*8, NCHR1(10)*8, NCHR2(10)*5
      LOGICAL RMARK
      LOGICAL L_LW1_tmp, L_LW2_tmp, L_LW3_tmp
      CHARACTER INPFLD*2, PATHWY(7)*2
      INTERFACE
         SUBROUTINE EXPATH(INPFLD,PATHWY,IPN,NOPATH)
            CHARACTER (LEN=2), INTENT(IN) :: INPFLD
            CHARACTER (LEN=2), INTENT(IN), DIMENSION(:) :: PATHWY
            INTEGER, INTENT(IN) :: IPN
            LOGICAL, INTENT(OUT) :: NOPATH
         END SUBROUTINE EXPATH
      END INTERFACE

C     Variable Initializations
      DATA (NCHR1(I),I=1,10) /'FIRST','SECOND','THIRD','FOURTH',
     &                        'FIFTH','SIXTH','SEVENTH','EIGHTH',
     &                        'NINTH','TENTH'/
      DATA (NCHR2(I),I=1,10) /'1ST','2ND','3RD','4TH','5TH',
     &                        '6TH','7TH','8TH','9TH','10TH'/

C     Variable Initializations
      MODNAM = 'PRESET'
      PREVSRCID = '        '
      PREVGRPID = '        '
      PATH  = '  '
      PPATH = '  '
      EOF = .FALSE.
      NSEC   = 0
      NPDMAX = 0
      NQF    = 0
      NBF    = 0
      NURB   = 0
      NVMAX  = 0
      ILINE  = 0

C     Initialize PATHWY array
      PATHWY(1) = 'CO'
      PATHWY(2) = 'SO'
      PATHWY(3) = 'RE'
      PATHWY(4) = 'ME'
      PATHWY(5) = 'OU'
      PATHWY(6) = '**'
      PATHWY(7) = 'EV'

      IPNUM  = 0
      IPPNUM = 0
C     Counters for the Receptor Groups
      IREC = 0
      ISTA = .FALSE.
      IEND = .FALSE.
      NEWID = .TRUE.
C     Initialize logical for urban option and multiple urban areas
      L_PRESET_URBAN = .FALSE.
      L_MULTURB      = .FALSE.
C     Initialize logical for the 'URBANSRC ALL' option
      L_URBAN_ALL = .FALSE.
C     Initialize logical for sector-varying O3 data
      L_O3Sector = .FALSE.
C     Initialize logical for O3 data inputs
      L_O3data = .FALSE.
C     Initialize logical for sector-varying Background data
      L_BGSector = .FALSE.
C     Initialize logical for BG data inputs
      L_BACKGRND = .FALSE.
C     Initialize temporary logicals for LowWind1, LowWind2, and LowWind3 BETA options
      L_LW1_tmp = .FALSE.
      L_LW2_tmp = .FALSE.
      L_LW3_tmp = .FALSE.
C     Initialize file format to 'FIX'; will be overridden if
C     user specified 'EXP' format on OU FILEFORM keyword
      FILE_FORMAT = 'FIX'

C --- Initialize all NO2 options to .FALSE.
      ARM    = .FALSE.
      ARM2   = .FALSE.
      OLM    = .FALSE.
      PVMRM  = .FALSE.
      PVMRM2 = .FALSE.

C --- Initialize MODOPTS_String character variable
      MODOPS_String = ''

C     Setup READ format and ECHO format for runstream record,
C     based on the ISTRG PARAMETER (set in MAIN1)
      WRITE(RDFRM,9100) ISTRG, ISTRG
 9100 FORMAT('(A',I4.4,',T1,',I4.4,'A1)')

C --- First loop through Runstream Records to identify
C     O3SECTOR option and BGSECTOR/BACKGRND options
C     since they are interdependant;
C     This avoids imposing a requirement on the
C     order of these keywords:
C     LOOP Through Input Runstream Records
      DO WHILE (.NOT. EOF)

C        Increment the Line Counter
         ILINE = ILINE + 1

C        READ Record to Buffers, as A'num' and 'num'A1 for 'num' = ISTRG.
C        Length of ISTRG is Set in PARAMETER Statement in MAIN1
         READ (INUNIT,RDFRM,END=998) RUNST1, (RUNST(I), I = 1, ISTRG)

C        Check for blank input record and cycle
         IF (LEN_TRIM(RUNST1) .EQ. 0) CYCLE

C        Convert Lower Case to Upper Case Letters           ---   CALL LWRUPR
         CALL LWRUPR

C        Define Fields on Card                              ---   CALL DEFINE
         CALL DEFINE

C        Get the Contents of the Fields                     ---   CALL GETFLD
         CALL GETFLD

C        If Blank Line, Then CYCLE to Next Card
         IF (BLINE) GO TO 10

C        Check for 'NO ECHO' In First Two Fields
         IF (FIELD(1) .EQ. 'NO' .AND. FIELD(2) .EQ. 'ECHO') THEN
C           Skip record with NO ECHO during PRESET stage of processing
            GO TO 10
         END IF

C        Extract Pathway ID From Field 1                    ---   CALL EXPATH
         CALL EXPATH(FIELD(1),PATHWY,7,NOPATH)

C        For Invalid Pathway and Comment Lines Skip to Next Record
         IF (NOPATH) THEN
C           Skip Error Message for PRESET stage of processing
            PATH = PPATH
            GO TO 10
         ELSE IF (PATH .EQ. '**') THEN
            GO TO 10
         END IF

C        Extract Keyword From Field 2                       ---   CALL EXKEY
         CALL EXKEY(FIELD(2),NOKEY)

         IF (NOKEY) THEN
C           No Keyword - Skip Error Message for PRESET stage
               PKEYWD = KEYWRD
            GO TO 10
         END IF

C        Save Current Path and Path Number as Previous Path and Number
         PPATH = PATH
         IPPNUM = IPNUM

C        First process cards to determine whether O3SECTOR and/or BGSECTOR keywords are used
         IF (PATH .EQ. 'CO') THEN

            IF (KEYWRD .EQ. 'O3SECTOR') THEN
C ---          Assign logical variable for O3SECTORs
               L_O3Sector = .TRUE.
               CYCLE
            
            ELSEIF (KEYWRD .EQ. 'OZONEFIL' .OR. 
     &              KEYWRD .EQ. 'O3VALUES' .OR. 
     &              KEYWRD .EQ. 'OZONEVAL') THEN
               L_O3Data = .TRUE.
               CYCLE

            ELSEIF (KEYWRD .EQ. 'EVENTFIL') THEN
               CYCLE

            END IF

         ELSE IF (PATH .EQ. 'SO') THEN

            IF (KEYWRD .EQ. 'BGSECTOR') THEN
C ---          Assign logical variable for BGSECTORs
               L_BGSector = .TRUE.
               CYCLE
            ELSE IF (KEYWRD .EQ. 'BACKGRND') THEN
C ---          Assign logical variable for BACKGRND
               L_BACKGRND = .TRUE.
               CYCLE
            ELSE IF (KEYWRD .EQ. 'INCLUDED' .AND.
     &               (.NOT.L_BGSector .OR. .NOT.L_BACKGRND)) THEN
C ---          Call PREINCLUD since BGSECTOR and/or BACKGRND may
C              be in an INCLUDED file
               CALL PREINCLUD
               CYCLE
            ELSE IF (KEYWRD .EQ. 'URBANSRC') THEN
               L_Urban = .TRUE.
               IF (FIELD(3) .EQ. 'ALL') THEN
                  L_Rural = .FALSE.
               END IF
            END IF

C ---       Exit pre-PRESET loop if both BGSECTOR and BACKGRND 
C           have been specified
            IF (L_BGSector .AND. L_BACKGRND) THEN
               EXIT
            END IF

         ELSE IF (PATH .EQ. 'RE') THEN
C ---       CO and SO pathways already processed:
C           Exit pre-PRESET loop
            EXIT

         END IF

         GO TO 10
 998     EOF = .TRUE.
 10      CONTINUE
      END DO

C --- Check for O3SECTOR keyword without O3 data inputs
      IF (L_O3Sector .AND. .NOT.L_O3Data) THEN
C ----   Issue warning message for O3SECTOR without O3 data
         CALL ERRHDL(PATH,MODNAM,'W','166','O3SECTOR')
      END IF
      IF (L_BGSector .AND. .NOT.L_BACKGRND) THEN
C ----   Issue warning message for O3SECTOR without O3 data
         CALL ERRHDL(PATH,MODNAM,'W','166','BGSECTOR')
      END IF

C --- Rewind Input Runstream File for complete pre-setup
      REWIND INUNIT

      PATH  = '  '
      PPATH = '  '
      IPNUM  = 0
      IPPNUM = 0
      EOF = .FALSE.

C     LOOP Through Input Runstream Records
      DO WHILE (.NOT. EOF)

C        Increment the Line Counter
         ILINE = ILINE + 1

C        READ Record to Buffers, as A'num' and 'num'A1 for 'num' = ISTRG.
C        Length of ISTRG is Set in PARAMETER Statement in MAIN1
         READ (INUNIT,RDFRM,END=999) RUNST1, (RUNST(I), I = 1, ISTRG)

C        Check for blank input record and cycle
         IF (LEN_TRIM(RUNST1) .EQ. 0) CYCLE

C        Convert Lower Case to Upper Case Letters           ---   CALL LWRUPR
         CALL LWRUPR

C        Define Fields on Card                              ---   CALL DEFINE
         CALL DEFINE

C        Get the Contents of the Fields                     ---   CALL GETFLD
         CALL GETFLD

C        Check for 'NO ECHO' In First Two Fields
         IF (FIELD(1) .EQ. 'NO' .AND. FIELD(2) .EQ. 'ECHO') THEN
C           Skip record with NO ECHO during PRESET stage of processing
            CYCLE
         END IF

C        Extract Pathway ID From Field 1                    ---   CALL EXPATH
         CALL EXPATH(FIELD(1),PATHWY,7,NOPATH)

C        For Invalid Pathway and Comment Lines Skip to Next Record
         IF (NOPATH) THEN
C           Skip Error Message for PRESET stage of processing
            PATH = PPATH
            CYCLE
         ELSE IF (PATH .EQ. '**') THEN
            CYCLE
         END IF

C        Extract Keyword From Field 2                       ---   CALL EXKEY
         CALL EXKEY(FIELD(2),NOKEY)

         IF (NOKEY) THEN
C           Invalid Keyword - Skip Error Message for PRESET stage
            PKEYWD = KEYWRD
            CYCLE
         END IF

C        Save Current Path and Path Number as Previous Path and Number
         PPATH = PATH
         IPPNUM = IPNUM

C        Process Cards to Determine Storage Requirements
         IF (PATH .EQ. 'CO') THEN
            IF (KEYWRD .EQ. 'MODELOPT') THEN
               DO I = 3, IFC
                  IF (FIELD(I) .EQ. 'CONC'  .OR.
     &                FIELD(I) .EQ. 'DEPOS' .OR.
     &                FIELD(I) .EQ. 'DDEP'  .OR.
     &                FIELD(I) .EQ. 'WDEP') THEN
                     NTYP = NTYP + 1
                  ELSE IF ((FIELD(I).EQ.'LOWWIND') .OR. 
     &                     (FIELD(I).EQ.'LOWWIND1')) THEN
                     L_LW1_tmp = .TRUE.
                  ELSE IF (FIELD(I) .EQ. 'LOWWIND2') THEN
                     L_LW2_tmp = .TRUE.
                  ELSE IF (FIELD(I) .EQ. 'LOWWIND3') THEN
                     L_LW3_tmp = .TRUE.
                  END IF
C                 Set PVMRM and OLM logicals for use in ALLSETUP
                  IF (FIELD(I) .EQ. 'PVMRM') THEN
                     PVMRM = .TRUE.
                  ELSE IF (FIELD(I) .EQ. 'PVMRM2') THEN
                     PVMRM2 = .TRUE.
                  ELSE IF (FIELD(I) .EQ. 'PSDCREDIT' )THEN
                     PSDCREDIT = .TRUE.
C----                Number of "source groups" will be set to 2 below for
C                    PSDCREDIT applications, to account for hardwired
C                    'NAAQS' and 'PSDINC' source groups, otherwise it 
C                    would be overwritten in SRCSIZ
                  ELSE IF (FIELD(I) .EQ. 'OLM') THEN
                     OLM = .TRUE.
                  ELSE IF (FIELD(I) .EQ. 'ARM2') THEN
                     ARM2 = .TRUE.
                  ELSE IF (FIELD(I) .EQ. 'ARM') THEN
                     ARM = .TRUE.
                  END IF
               END DO

            ELSE IF (KEYWRD .EQ. 'AVERTIME') THEN
               DO I = 3, IFC
                  IF (FIELD(I).NE.'PERIOD' .AND. 
     &                FIELD(I).NE.'ANNUAL') THEN
                     NAVE = NAVE + 1
                  END IF
               END DO

            ELSE IF (KEYWRD .EQ. 'URBANOPT') THEN
               NURB = NURB + 1
C----          Set preliminary flag for URBAN option, used to allow flexibility in 
C              order of CO pathway keywords for URBANOPT
               L_PRESET_URBAN = .TRUE.
               IF (NURB .GT. 1) THEN
                  L_MULTURB = .TRUE.
               END IF

            ELSE IF (KEYWRD .EQ. 'O3SECTOR') THEN
C ---          Assign logical variable for O3SECTORs
               L_O3Sector = .TRUE.
C ---          Set maximum array limit for number of ozone sectors
               NUMO3Sects = IFC - 2

            ELSE IF (KEYWRD .EQ. 'O3VALUES') THEN
C ---          Set maximum array limit for temporally-varying
C              ozone concentrations for O3VALUES keyword
C              Assign field index to 4 if O3SECTOR is used, otherwise 3
               IF (L_O3Sector) THEN
                  K = 4
               ELSE
                  K = 3
               END IF
               IF (FIELD(K) .EQ. 'ANNUAL') THEN
                  NO3F = MAX( NO3F, 1)
               ELSE IF (FIELD(K) .EQ. 'SEASON') THEN
                  NO3F = MAX( NO3F, 4)
               ELSE IF (FIELD(K) .EQ. 'MONTH') THEN
                  NO3F = MAX( NO3F, 12)
               ELSE IF (FIELD(K) .EQ. 'HROFDY') THEN
                  NO3F = MAX( NO3F, 24)
               ELSE IF (FIELD(K) .EQ. 'WSPEED') THEN
                  NO3F = MAX( NO3F, 6)
               ELSE IF (FIELD(K) .EQ. 'SEASHR') THEN
                  NO3F = MAX( NO3F, 96)
               ELSE IF (FIELD(K) .EQ. 'HRDOW') THEN
                  NO3F = MAX( NO3F, 72)
               ELSE IF (FIELD(K) .EQ. 'HRDOW7') THEN
                  NO3F = MAX( NO3F, 168)
               ELSE IF (FIELD(K) .EQ. 'SHRDOW') THEN
                  NO3F = MAX( NO3F, 288)
               ELSE IF (FIELD(K) .EQ. 'SHRDOW7') THEN
                  NO3F = MAX( NO3F, 672)
               ELSE IF (FIELD(K) .EQ. 'MHRDOW') THEN
                  NO3F = MAX( NO3F, 864)
               ELSE IF (FIELD(K) .EQ. 'MHRDOW7') THEN
                  NO3F = MAX( NO3F, 2016)
               END IF

            ELSE IF (KEYWRD .EQ. 'LOW_WIND') THEN
C----          Check for LOW_WIND keyword for users to adjust parameters
C              used with the LowWind1 or LowWind2 BETA options, to avoid
C              duplicate conflicting messages regarding these parameters
C----          First check for LowWind1 option, then for LowWind2 option
               IF (L_LW1_tmp) THEN
                  IF (IFC.GT.3) THEN
C----                User specified both SVmin and WSmin
                     L_UserSVmin = .TRUE.
                     L_UserWSmin = .TRUE.
                  ELSE
C----                User specified only SVmin
                     L_UserSVmin = .TRUE.
                  END IF

               ELSE IF (L_LW2_tmp) THEN
                  IF (L_LW2_tmp .AND. IFC.GT.4) THEN
C----                User specified SVmin, WSmin, and FRANmax
                     L_UserSVmin   = .TRUE.
                     L_UserWSmin   = .TRUE.
                     L_UserFRANmax = .TRUE.
                  ELSE IF (L_LW2_tmp .AND. IFC.GT.3) THEN
C----                User specified both SVmin and WSmin
                     L_UserSVmin   = .TRUE.
                     L_UserWSmin   = .TRUE.
                  ELSE
C----                User specified only SVmin
                     L_UserSVmin   = .TRUE.
                  END IF
 
               ELSE IF (L_LW3_tmp) THEN
                  IF (L_LW3_tmp .AND. IFC.GT.4) THEN
C----                User specified SVmin, WSmin, and FRANmax
                     L_UserSVmin   = .TRUE.
                     L_UserWSmin   = .TRUE.
                     L_UserFRANmax = .TRUE.
                  ELSE IF (L_LW3_tmp .AND. IFC.GT.3) THEN
C----                User specified both SVmin and WSmin
                     L_UserSVmin   = .TRUE.
                     L_UserWSmin   = .TRUE.
                  ELSE
C----                User specified only SVmin
                     L_UserSVmin   = .TRUE.
                  END IF
               END IF

            END IF

         ELSE IF (PATH .EQ. 'SO') THEN
            CALL SRCSIZ

         ELSE IF (PATH .EQ. 'RE') THEN
            EVONLY = .FALSE.
            CALL RECSIZ

         ELSE IF (PATH .EQ. 'EV') THEN
            EVONLY = .TRUE.
            IF (KEYWRD .EQ. 'EVENTPER') THEN
               NEVE = NEVE + 1
            ELSE IF (KEYWRD .EQ. 'INCLUDED') THEN
               CALL PREINCLUD
            END IF

         ELSE IF (PATH .EQ. 'ME' .AND. KEYWRD .EQ. 'SURFDATA') THEN
C           Read start year from SURFDATA card to establish date window
            CALL SET_WINDOW

         ELSE IF (PATH .EQ. 'ME' .AND. KEYWRD .EQ. 'NUMYEARS') THEN
C ---       Set number of years for allocating the MAXDCONT arrays
            IF (IFC .EQ. 3) THEN
               CALL STONUM(FIELD(3),ILEN_FLD,FNUM,IMIT)
               IF (IMIT .EQ. 1) THEN
                  NYEARS = NINT(FNUM)
               END IF
            END IF

         ELSE IF (PATH .EQ. 'OU') THEN
            IF(KEYWRD .EQ. 'RECTABLE') THEN
C              Begin LOOP Through Fields
               DO I = 4, IFC
C ---             Skip processing of fields if IFC > IFMAX
                  IF (I .GT. IFMAX) EXIT
C                 Retrieve The High Value
                  CALL FSPLIT(PATH,KEYWRD,FIELD(I),ILEN_FLD,'-',RMARK,
     &                        LPRD,HPRD)
                  ISPRD = 0
                  IEPRD = 0
C ---             First check for simple numeric value
                  CALL STONUM(LPRD,ILEN_FLD,FNUM,IMIT)
                  IF (IMIT .EQ. 1) THEN
                     ISPRD = INT(FNUM)
                  END IF   
                  CALL STONUM(HPRD,ILEN_FLD,FNUM,IMIT)
                  IF (IMIT .EQ. 1) THEN
                     IEPRD = INT(FNUM)
                  END IF
C ---             Now check for character strings NCHR1 or NCHR2
                  DO J = 1, 10
                     IF (LPRD.EQ.NCHR1(J) .OR.
     &                   LPRD.EQ.NCHR2(J)) ISPRD = J
                     IF (HPRD.EQ.NCHR1(J) .OR.
     &                   HPRD.EQ.NCHR2(J)) IEPRD = J
                  END DO
                  IF (ISPRD .GT. 999 .OR. IEPRD .GT. 999) THEN
C                    Write Error Message:Illegal Parameter Field
                     CALL ERRHDL(PATH,MODNAM,'E','203','HIVALU')
                     CYCLE
                  END IF
                  IF (ISPRD .GT. NVAL) THEN
                     NVAL = ISPRD
                  END IF
                  IF (IEPRD .GT. NVAL) THEN
                     NVAL = IEPRD
                  END IF
C              End LOOP Through Fields
               END DO

            ELSE IF (KEYWRD .EQ. 'MAXTABLE' .OR. 
     &               KEYWRD .EQ. 'RANKFILE') THEN
C              Set Number of Maximum Values to Sort
               CALL STONUM(FIELD(4),ILEN_FLD,FNUM,IMIT)
               IF (IMIT .NE. 1) THEN
C                 Invalid Numerical Field
                  GO TO 999
               END IF
               INUM = NINT(FNUM)
               IF (INUM .GT. NMXVAL) THEN
                  NMXVAL = INUM
               END IF

            ELSE IF (KEYWRD .EQ. 'SEASONHR') THEN
C              Set SEASONHR logical flag to account for SHVALS array needs
               SEASONHR = .TRUE.
            
            ELSE IF (KEYWRD .EQ. 'MAXDCONT') THEN
C              Set MAXDCONT logical flag to account for MAXDCONT array needs
               L_MAXDCONT = .TRUE.
            
            ELSE IF (KEYWRD .EQ. 'FILEFORM' .AND.
     &               FIELD(3)(1:3) .EQ. 'EXP') THEN
C ---          Check for FILEFORM keyword with 'EXP' format in order to 
C              include correct format in output file headers
               FILE_FORMAT = 'EXP'
            
            ELSE IF (KEYWRD .EQ. 'NOHEADER') THEN
C              Set NOHEADER logical flag to suppress output file headers
               DO I = 3, IFC
                  IF (FIELD(I) .EQ. 'ALL') THEN
C                    No headers for any ouput file type
                     L_NoHeader(:) = .TRUE.
                     EXIT
                  ELSE IF (FIELD(I) .EQ. 'MAXIFILE') THEN
C                    No headers for MAXIFILE
                     L_NoHeader(1) = .TRUE.
                  ELSE IF (FIELD(I) .EQ. 'POSTFILE') THEN
C                    No headers for POSTFILE
                     L_NoHeader(2) = .TRUE.
                  ELSE IF (FIELD(I) .EQ. 'PLOTFILE') THEN
C                    No headers for PLOTFILE
                     L_NoHeader(3) = .TRUE.
                  ELSE IF (FIELD(I) .EQ. 'SEASONHR') THEN
C                    No headers for SEASONHR
                     L_NoHeader(4) = .TRUE.
                  ELSE IF (FIELD(I) .EQ. 'RANKFILE') THEN
C                    No headers for RANKFILE
                     L_NoHeader(5) = .TRUE.
                  ELSE IF (FIELD(I) .EQ. 'MAXDAILY') THEN
C                    No headers for MAXDAILY
                     L_NoHeader(6) = .TRUE.
                  ELSE IF (FIELD(I) .EQ. 'MXDYBYYR') THEN
C                    No headers for MXDYBYYR
                     L_NoHeader(7) = .TRUE.
                  ELSE IF (FIELD(I) .EQ. 'MAXDCONT') THEN
C                    No headers for MAXDCONT
                     L_NoHeader(8) = .TRUE.
                  END IF
               END DO
      
            ELSE IF (KEYWRD .EQ. 'FINISHED') THEN
C              Check for 'OU FINISHED' Card.  Exit DO WHILE Loop By Branching
C              to Statement 999 in Order to Avoid Reading a ^Z "End of File"
C              Marker That May Be Present For Some Editors.
               GO TO 999

            END IF
            
         END IF

C        Store the Current Keyword as the Previous Keyword
         PKEYWD = KEYWRD

         GO TO 11
 999     EOF = .TRUE.
 11      CONTINUE
      END DO

C---- Check for PSDCREDIT option and set number of source groups
      IF (PSDCREDIT) THEN
C----    Set number of "source groups" to 2 for PSDCREDIT applications,
C        to account for hardwired 'NAAQS' and 'PSDINC' source groups
         NGRP = 2
      END IF

C --- Determine maximum number of vertices for AREA sources, including 
C     AREAPOLY and AREACIRC source types, LINE sources, and OPENPIT sources.  
      IF (NAREA .EQ. 0 .AND. NLINE .EQ. 0 .AND. NPIT .EQ. 0) THEN
C        No area, line, or openpit sources, set NVMAX to 0
         NVMAX = 0
      ELSE
         IF (NPOLY .EQ. 0 .AND. NCIRC .EQ. 0) THEN
C           No AREAPOLY or AREACIRC sources; initialize NVMAX to 8 
C           for rectangular AREA, LINE, and/or OPENPIT sources
            NVMAX = 8
         ELSE
C           Assign NVMAX to at least 8 to handle rectangular AREA,
C           LINE, and/or OPENPIT sources
            NVMAX = MAX( NVMAX, 8 )
         END IF
      END IF
C     Calculate value of NVMAX2
      NVMAX2 = NVMAX*2

C     Rewind File and Reinitialize Line Number Counter for SETUP
      REWIND INUNIT
      ILINE = 0
      PNETID = '        '

C     Ensure that certain array limits are not < 1.
      NSRC = MAX( NSRC, 1)
      NGRP = MAX( NGRP, 1)
      NREC = MAX( NREC, 1)
C     Set NARC = NREC temporarily for allocating setup arrays
      NARC = NREC
      NAVE = MAX( NAVE, 1)
      NTYP = MAX( NTYP, 1)
      NNET = MAX( NNET, 1)
      IXM  = MAX( IXM , 1)
      IYM  = MAX( IYM , 1)

C     Assign array limits to REAL for calculation of STORE
      RSRC   = REAL(NSRC)
      RSEC   = REAL(NSEC)
      RGRP   = REAL(NGRP)
      RREC   = REAL(NREC)
      RURB   = REAL(NURB)
      RARC   = REAL(NARC)
      RAVE   = REAL(NAVE)
      RVAL   = REAL(NVAL)
      RTYP   = REAL(NTYP)
      RMAX   = REAL(NMAX)
      RNET   = REAL(NNET)
      RXM    = REAL(IXM )
      RYM    = REAL(IYM )
      REVE   = REAL(NEVE)
      ROLM   = REAL(NOLM)
      RPSD   = REAL(NPSD)
      RQF    = REAL(NQF)
      RPDMAX = REAL(NPDMAX)
      RVMAX  = REAL(NVMAX)
      RPAIR  = REAL(NPAIR)
      RHIANN = REAL(NHIANN)

      RETURN
      END

      SUBROUTINE PREINCLUD
C***********************************************************************
C*                PREINCLUD Module of the AMS/EPA Regulatory Model - AERMOD
C*
C*       PURPOSE: To read an external receptor/source file using the
C*                INCLUDED keyword.
C*
C*       PROGRAMMER: Roger Brode
C*
C*       DATE:    September 24, 1996
C*
C*       MODIFIED:   
C*                   
C*       INPUTS: 
C*
C*       OUTPUTS:
C*               
C*
C*       CALLED FROM:   PRESET, SRCSIZ, RECSIZ
C***********************************************************************
        
C*    Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: I, ITEMPL
      LOGICAL NOPATH, NOKEY      
      CHARACTER RDFRM*20
      CHARACTER INPFLD*2, PATHWY(7)*2
      INTERFACE
         SUBROUTINE EXPATH(INPFLD,PATHWY,IPN,NOPATH)
            CHARACTER (LEN=2), INTENT(IN) :: INPFLD
            CHARACTER (LEN=2), INTENT(IN), DIMENSION(:) :: PATHWY
            INTEGER, INTENT(IN) :: IPN
            LOGICAL, INTENT(OUT) :: NOPATH
         END SUBROUTINE EXPATH
      END INTERFACE

C*    Variable Initializations
      MODNAM = 'PREINCLUD'
      EOF = .FALSE.
      ILINE  = 0
      ITEMPL = 0
      
C     Initialize PATHWY array
      PATHWY(1) = 'CO'
      PATHWY(2) = 'SO'
      PATHWY(3) = 'RE'
      PATHWY(4) = 'ME'
      PATHWY(5) = 'OU'
      PATHWY(6) = '**'
      PATHWY(7) = 'EV'

C     Setup READ format and ECHO format for runstream record,
C     based on the ISTRG PARAMETER (set in MAIN1)
      WRITE(RDFRM,9100) ISTRG, ISTRG
 9100 FORMAT('(A',I4.4,',T1,',I4.4,'A1)')
      
      IF (IFC .EQ. 3) THEN
C        Retrieve Included Filename as Character Substring to Maintain Case
         IF ((LOCE(3)-LOCB(3)) .LE. (ILEN_FLD - 1) ) THEN
C           Retrieve Filename as Character Substring to Maintain Original Case
C           Also Check for Filename Larger Than ILEN_FLD Characters
            INCFIL = RUNST1(LOCB(3):LOCE(3))
            OPEN (INCUNT,FILE=INCFIL,ACTION='READ',STATUS='OLD',ERR=99)
         ELSE
C           WRITE Error Message:  INCFIL Field is Too Long
            WRITE(DUMMY,'(I8)') ILEN_FLD
            CALL ERRHDL(PATH,MODNAM,'E','291',DUMMY)
            RETURN
         END IF

      ELSE IF (IFC .GE. 4) THEN
C        Too Many Parameters
         RETURN
      ELSE
C        No Parameters Specified
         RETURN
      END IF

      GO TO 1001

C     Write Out Error Message for File OPEN Error
99    CALL ERRHDL(PATH,MODNAM,'E','500','INCFILE ')
      RETURN

1001  CONTINUE

C     LOOP Through Input Runstream Records
      DO WHILE (.NOT. EOF)

C        Increment the Line Counter.  It was Initially Set to 1, to Handle
C        the Code in Subroutine DEFINE
         ILINE = ILINE + 1
         
C        READ Record to Buffers, as A'num' and 'num'A1 where 'num' = ISTRG.
C        Length of ISTRG is Set in PARAMETER Statement in MAIN1
         READ (INCUNT,RDFRM,END=999,ERR=888) RUNST1, 
     &                                      (RUNST(I), I = 1, ISTRG)

C        If ILINE=1, reset ILINE temporarily to avoid the
C        check for column shift in subroutine DEFINE
         IF (ILINE .EQ. 1) THEN
            ILINE  = 2
            ITEMPL = 1
         END IF

C        Check for blank input record and cycle
         IF (LEN_TRIM(RUNST1) .EQ. 0) CYCLE

C        Convert Lower Case to Upper Case Letters           ---   CALL LWRUPR
         CALL LWRUPR

C        Define Fields on Card                              ---   CALL DEFINE
         CALL DEFINE

C        Reset ILINE if needed
         IF (ITEMPL .EQ. 1) THEN
            ILINE  = 1
            ITEMPL = 0
         END IF

C        Get the Contents of the Fields                     ---   CALL GETFLD
         CALL GETFLD

C        Check for 'NO ECHO' In First Two Fields
         IF (FIELD(1) .EQ. 'NO' .AND. FIELD(2) .EQ. 'ECHO') THEN
C           Skip record with NO ECHO during PREINCLUD stage of processing
            CYCLE
         END IF

C        Extract Pathway ID From Field 1                    ---   CALL EXPATH
         CALL EXPATH(FIELD(1),PATHWY,7,NOPATH)

C        For Invalid Pathway and Comment Lines Skip to Next Record
         IF (NOPATH) THEN
C           Skip Error Message for PREINCLUD stage of processing
            PATH = PPATH
            CYCLE
         ELSE IF (PATH .EQ. '**') THEN
            CYCLE
         END IF

C        Extract Keyword From Field 2                       ---   CALL EXKEY
         CALL EXKEY(FIELD(2),NOKEY)

         IF (NOKEY) THEN
C           Invalid Keyword - Skip Error Message for PREINCLUD stage
            PKEYWD = KEYWRD
            CYCLE
         END IF

C        Save Current Path and Path Number as Previous Path and Number
         PPATH = PATH
         IPPNUM = IPNUM

C        Process Input Card Based on Pathway
         IF (PATH .EQ. 'SO') THEN
C           Process SOurce Pathway Cards                    ---   CALL SRCSIZ
            CALL SRCSIZ
         ELSE IF (PATH .EQ. 'RE') THEN
C           Process REceptor Pathway Cards                  ---   CALL RECSIZ
            CALL RECSIZ
         ELSE IF (PATH .EQ. 'EV') THEN
            IF (KEYWRD .EQ. 'EVENTPER') THEN
               NEVE = NEVE + 1
            END IF
         END IF

C        Store the Current Keyword as the Previous Keyword
         PKEYWD = KEYWRD

         GO TO 11
 999     EOF = .TRUE.
         REWIND INCUNT
 11      CONTINUE

      END DO
      EOF = .FALSE.

      GO TO 1002
      
 888  CONTINUE
C --- Error occurred reading the included file, issue error message
      CALL ERRHDL(PATH,MODNAM,'E','510','INCLUDED')
      RUNERR = .TRUE.
      
1002  CONTINUE
      
C     Close the INCLUDED File
      CLOSE (INCUNT)
      
      RETURN
      END

      SUBROUTINE SRCSIZ
C***********************************************************************
C                 SRCSIZ Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: To preprocess source inputs to determine
C                 storage requirements
C
C        PROGRAMMER: Roger Brode
C
C        DATE:    September 24, 1996
C
C        MODIFIED:   Modified to add counter for number of lines in a 
C                    buoyant line source
C                    Amec Foster Wheeler, 06/30/2015
C
C        MODIFIED:   Modified to allow for the use of URBANSRC ALL on 
C                    the SO pathway to indicate that all sources are 
C                    to be treated as URBAN sources. This option assumes 
C                    that only one (1) urban area has been defined using 
C                    the CO URBANOPT keyword.
C                    R.W. Brode, U.S. EPA/OAQPS/AQMG, 02/29/2012
C
C        MODIFIED:   To include options to vary emissions by
C                    hour-of-day, and day-of-week (HRDOW and HRDOW7).
C                    Modified method for determining maximum number
C                    of vertices for AREAPOLY sources for more precise
C                    and efficient memory allocation.  Also included
C                    allocation of arrays for building downwash data.
C                    R.W. Brode, U.S. EPA/OAQPS/AQMG, 10/19/2009
C
C        MODIFIED:   To include new options incorporated in version
C                    dated 06341.
C                    R.W. Brode, U.S. EPA/OAQPS/AQMG - 12/07/2006
C
C        MODIFIED:   To include options to vary emissions by month,
C                    hour-of-day, and day-of-week (MHRDOW and MHRDOW7).
C                    R.W. Brode, MACTEC (f/k/a PES), Inc., 06/22/05
C
C        MODIFIED:   To include an option to vary emissions by season,
C                    hour-of-day, and day-of-week (SHRDOW).
C                    R.W. Brode, PES, 4/10/2000
C
C        INPUTS:  Pathway (RE) and Keyword
C
C        OUTPUTS: Receptor Arrays
C                 Receptor Setup Status Switches
C
C        CALLED FROM:   PRESET
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      USE BUOYANT_LINE
      
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: I, K, NUM_SRCIDS
      CHARACTER (LEN=12) :: TMPSRCID, SAVESRCID
      ALLOCATABLE :: TMPSRCID(:), SAVESRCID(:)

      SAVE TMPSRCID, SAVESRCID, NUM_SRCIDS
      
C     Variable Initializations
      MODNAM = 'SRCSIZ'

      IF (KEYWRD .EQ. 'STARTING') THEN
C        Initialize Counters and Set Status Switch
         NSRC = 0
         NGRP = 0
         NOLM = 0
         NPSD = 0
         NQF  = 0
         NBF  = 0
         NSEC = 0
         NPIT = 0
         NPNT = 0
         NVOL = 0
         NAREA = 0
         NPOLY = 0
         NCIRC = 0
         NLINE = 0
C ---    Add variable for number of lines in buoyant source
         NBLINES = 0
C ---    Add variable for number of urban sources
         NURBSRC = 0
         NVMAX  = 0
         NVTEMP = 0
         NPTEMP = 0
         NPDMAX = 0
         PREVSRCID = '        '
         PREVGRPID = '        '

      ELSE IF (KEYWRD .EQ. 'LOCATION') THEN
         NSRC = NSRC + 1
         IF (FIELD(4) .EQ. 'POINT') THEN
            NPNT = NPNT + 1
         ELSE IF (FIELD(4) .EQ. 'VOLUME') THEN
            NVOL = NVOL + 1
         ELSE IF (FIELD(4) .EQ. 'LINE') THEN
            NLINE = NLINE + 1
         ELSE IF (FIELD(4) .EQ. 'BUOYLINE') THEN
            NBLINES = NBLINES + 1
         ELSE IF (FIELD(4)(1:4) .EQ. 'AREA') THEN
            NAREA = NAREA + 1
            IF (FIELD(4) .EQ. 'AREAPOLY') THEN
               NPOLY = NPOLY + 1
            ELSE IF (FIELD(4) .EQ. 'AREACIRC') THEN
C ---          Increment counter for number of AREACIRC sources
               NCIRC = NCIRC + 1
C ---          Save AREACIRC source IDs in temporary arrays in order
C              to check SRCPARAM keyword inputs for NVERTS parameter.
C              First check allocation status of TMPSRCID array; 
C              if allocated, then save TMPSRCID array, deallocate,
C              and reallocate based on current number of sources.
C              The end result of TMPSRCID will be an array with 
C ---          only the source IDs for AREACIRC sources.
               IF (ALLOCATED(TMPSRCID)) THEN
                  IF (ALLOCATED(SAVESRCID)) THEN
                     SAVESRCID = TMPSRCID
                     NUM_SRCIDS = SIZE(TMPSRCID)
                     DEALLOCATE (TMPSRCID)
                     ALLOCATE (TMPSRCID(NSRC))
                     TMPSRCID(1:NUM_SRCIDS)  = SAVESRCID
                     TMPSRCID(NUM_SRCIDS+1:) = ' '
                  END IF
               ELSE
                  ALLOCATE (TMPSRCID(NSRC))
                  ALLOCATE (SAVESRCID(NSRC))
                  NUM_SRCIDS = NSRC
               END IF
               TMPSRCID(NSRC) = FIELD(3)
               NUM_SRCIDS = NSRC
            END IF
         ELSE IF (FIELD(4) .EQ. 'OPENPIT') THEN
            NPIT = NPIT + 1
         END IF

      ELSE IF (KEYWRD .EQ. 'SRCPARAM') THEN
C ---    Check for AREACIRC sources with user-specified NVERTS
         DO I = 1, NSRC
C ---       Exit loop if number of temporary AREACIRC source IDs 
C           is less than current loop index (I)
            IF (NUM_SRCIDS .LT. I) EXIT
            
            IF (FIELD(3) .EQ. TMPSRCID(I)) THEN
C ---          This is an AREACIRC source: check for NVERTS input
               IF (IFC .GE. 7) THEN
C ---             Set maximum number of vertices based on 
C                 user-specified number for AREACIRC source
                  CALL STODBL(FIELD(7),ILEN_FLD,DNUM,IMIT)
C                 Check The Numerical Field
                  IF (IMIT .NE. 1) THEN
                     CALL ERRHDL(PATH,MODNAM,'E','208',KEYWRD)
                  ELSE
C ---                Adjust NVMAX if needed based on number of
C                    vertices specified for this AREACIRC source +4
                     NVMAX = MAX( NVMAX, IDNINT(DNUM)+4 )
                  END IF
               ELSE
C ---             User did not specify number of vertices for
C                 this AREACIRC source; adjust NVMAX if needed
C                 based on default value of 20 for NVERTS (+4)
                  NVMAX = MAX( NVMAX, 24 )
               END IF
C ---          AREACIRC source ID was found, EXIT loop
               EXIT
            END IF
            
         END DO
         
      ELSE IF (KEYWRD .EQ. 'AREAVERT') THEN
         IF (FIELD(3) .EQ. PREVSRCID) THEN
            NVTEMP = NVTEMP + IFC - 3
C ---       Set NVMAX based on current number of 
C           vertices for this source (NVTEMP/2) + NVPOLY
C           (where the NVPOLY PARAMETER is assigned a value of 12)
C           to account for maximum number of sides 
C           for transect through source
            NVMAX  = MAX( NVMAX, NVPOLY+NINT(FLOAT(NVTEMP/2)) )
         ELSE
C ---       This is first AREAVERT keyword for this AREAPOLY source.
C           Assign NVTEMP based on number of data fields specified.
            NVTEMP = IFC - 3
            NVMAX  = MAX( NVMAX, NVPOLY+NINT(FLOAT(NVTEMP/2)) )
            PREVSRCID = FIELD(3)
         END IF

      ELSE IF (KEYWRD.EQ.'PARTDIAM') THEN
         IF (FIELD(3) .EQ. PREVSRCID) THEN
            NPTEMP = NPTEMP + IFC - 3
C ---       Set NPDMAX based on current number of 
C           particle size categories for this source 
            NPDMAX = MAX( NPDMAX, NPTEMP )
         ELSE
            NPTEMP = IFC - 3
            NPDMAX = MAX( NPDMAX, NPTEMP )
            PREVSRCID = FIELD(3)
         END IF

      ELSE IF ((KEYWRD.EQ.'BUILDHGT' .OR.
     &          KEYWRD.EQ.'BUILDWID' .OR.
     &          KEYWRD.EQ.'BUILDLEN')) THEN
         NSEC = 36

      ELSE IF (KEYWRD .EQ. 'METHOD_2') THEN
         NPDMAX = MAX( NPDMAX, 1 )

      ELSE IF (KEYWRD .EQ. 'EMISFACT') THEN
         IF (FIELD(4) .EQ. 'SEASON') THEN
            NQF = MAX( NQF, 4)
         ELSE IF (FIELD(4) .EQ. 'MONTH') THEN
            NQF = MAX( NQF, 12)
         ELSE IF (FIELD(4) .EQ. 'HROFDY') THEN
            NQF = MAX( NQF, 24)
         ELSE IF (FIELD(4) .EQ. 'WSPEED') THEN
            NQF = MAX( NQF, 6)
         ELSE IF (FIELD(4) .EQ. 'SEASHR') THEN
            NQF = MAX( NQF, 96)
         ELSE IF (FIELD(4) .EQ. 'HRDOW') THEN
            NQF = MAX( NQF, 72)
         ELSE IF (FIELD(4) .EQ. 'HRDOW7') THEN
            NQF = MAX( NQF, 168)
         ELSE IF (FIELD(4) .EQ. 'SHRDOW') THEN
            NQF = MAX( NQF, 288)
         ELSE IF (FIELD(4) .EQ. 'SHRDOW7') THEN
            NQF = MAX( NQF, 672)
         ELSE IF (FIELD(4) .EQ. 'MHRDOW') THEN
            NQF = MAX( NQF, 864)
         ELSE IF (FIELD(4) .EQ. 'MHRDOW7') THEN
            NQF = MAX( NQF, 2016)
         END IF

      ELSE IF (KEYWRD .EQ. 'BGSECTOR') THEN
C ---    Assign logical variable for BGSECTORs
         L_BGSector = .TRUE.
C ---    Set maximum array limit for number of BACKGRND sectors 
         NUMBGSects = IFC - 2

      ELSE IF (KEYWRD .EQ. 'BACKGRND') THEN
         L_BACKGRND = .TRUE.
C ---    Assign field index to 4 if BGSECTOR is used, otherwise 3
         IF (L_BGSector) THEN
            K = 4
         ELSE
            K = 3
         END IF
         IF (FIELD(K) .EQ. 'ANNUAL') THEN
            NBF = MAX( NBF, 1)
         ELSE IF (FIELD(K) .EQ. 'SEASON') THEN
            NBF = MAX( NBF, 4)
         ELSE IF (FIELD(K) .EQ. 'MONTH') THEN
            NBF = MAX( NBF, 12)
         ELSE IF (FIELD(K) .EQ. 'HROFDY') THEN
            NBF = MAX( NBF, 24)
         ELSE IF (FIELD(K) .EQ. 'WSPEED') THEN
            NBF = MAX( NBF, 6)
         ELSE IF (FIELD(K) .EQ. 'SEASHR') THEN
            NBF = MAX( NBF, 96)
         ELSE IF (FIELD(K) .EQ. 'HRDOW') THEN
            NBF = MAX( NBF, 72)
         ELSE IF (FIELD(K) .EQ. 'HRDOW7') THEN
            NBF = MAX( NBF, 168)
         ELSE IF (FIELD(K) .EQ. 'SHRDOW') THEN
            NBF = MAX( NBF, 288)
         ELSE IF (FIELD(K) .EQ. 'SHRDOW7') THEN
            NBF = MAX( NBF, 672)
         ELSE IF (FIELD(K) .EQ. 'MHRDOW') THEN
            NBF = MAX( NBF, 864)
         ELSE IF (FIELD(K) .EQ. 'MHRDOW7') THEN
            NBF = MAX( NBF, 2016)
         END IF

      ELSE IF (KEYWRD .EQ. 'URBANSRC' .AND. 
     &                     IFC .EQ. 3 .AND. FIELD(3) .EQ. 'ALL'
     &                                .AND. L_PRESET_URBAN) THEN
         IF (.NOT. L_MULTURB) THEN
C           Set logical for URBANSRC ALL option (not applicable for
C           multiple urban areas)
            L_URBAN_ALL = .TRUE.
C ---       Assign number of urban sources
            NURBSRC = NSRC
         ELSE 
C           Issue ERROR message for URBANSRC ALL option with multiple 
C           urban areas
            CALL ERRHDL(PATH,MODNAM,'E','279','URBANSRC ALL')
         END IF

C --- Check for number of urbansrc
      ELSE IF (KEYWRD .EQ. 'URBANSRC' .AND. .NOT.L_MULTURB) THEN
         NURBSRC = NURBSRC + IFC - 2
      
      ELSE IF (KEYWRD .EQ. 'URBANSRC' .AND. L_MULTURB) THEN
         NURBSRC = NURBSRC + IFC - 3
      
      ELSE IF (KEYWRD .EQ. 'OLMGROUP') THEN
         IF (NOLM .EQ. 0) PREVGRPID = '        '
         IF (FIELD(3) .NE. PREVGRPID) THEN
            NOLM = NOLM + 1
            PREVGRPID = FIELD(3)
         END IF

      ELSE IF (KEYWRD .EQ. 'PSDGROUP') THEN
         IF (NPSD .EQ. 0) PREVGRPID = '        '
         IF (FIELD(3) .NE. PREVGRPID) THEN
            NPSD = NPSD + 1
            PREVGRPID = FIELD(3)
         END IF

      ELSE IF (KEYWRD .EQ. 'SRCGROUP') THEN
         IF (NGRP .EQ. 0) PREVGRPID = '        '
         IF (FIELD(3) .NE. PREVGRPID) THEN
            NGRP = NGRP + 1
            PREVGRPID = FIELD(3)
         END IF

      ELSE IF (KEYWRD .EQ. 'INCLUDED') THEN
         CALL PREINCLUD

      ELSE IF (KEYWRD .EQ. 'FINISHED') THEN
         RETURN

      END IF

 999  RETURN
      END

      SUBROUTINE RECSIZ
C***********************************************************************
C                 RECSIZ Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: To preprocess receptor inputs to determine
C                 storage requirements
C
C        PROGRAMMER: Roger Brode
C
C        DATE:    September 24, 1996
C
C        INPUTS:  Pathway (RE) and Keyword
C
C        OUTPUTS: Receptor Arrays
C                 Receptor Setup Status Switches
C
C        CALLED FROM:   PRESET
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

C     Variable Initializations
      MODNAM = 'RECSIZ'

      IF (KEYWRD .EQ. 'STARTING') THEN
C        Initialize Counters and Set Status Switch
         NREC = 0
         NNET = 0
         IXM  = 0
         IYM  = 0
         PXSOID = ' '
         ISTA = .FALSE.
      ELSE IF (KEYWRD .EQ. 'GRIDCART') THEN
C        Process Cartesian Grid Receptor Network            ---   CALL PRECART
         CALL PRECART
      ELSE IF (KEYWRD .EQ. 'GRIDPOLR') THEN
C        Process Polar Receptor Network                     ---   CALL PREPOLR
         CALL PREPOLR
      ELSE IF (KEYWRD .EQ. 'DISCCART') THEN
         NREC = NREC + 1
      ELSE IF (KEYWRD .EQ. 'EVALCART') THEN
         NREC = NREC + 1
      ELSE IF (KEYWRD .EQ. 'DISCPOLR') THEN
         NREC = NREC + 1
      ELSE IF (KEYWRD .EQ. 'INCLUDED') THEN
         CALL PREINCLUD
      ELSE IF (KEYWRD .EQ. 'FINISHED') THEN
         RETURN
      END IF

 999  RETURN
      END

      SUBROUTINE PRECART
C***********************************************************************
C                 PRECART Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Processes Cartesian Grid Receptor Network Inputs
C
C        PROGRAMMER:  Roger Brode
C
C        DATE:    September 24, 1996
C
C        INPUTS:  Input Runstream Image Parameters
C
C        OUTPUTS: Cartesian Grid Receptor Network Inputs
C
C        CALLED FROM:   RECSIZ
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

C     Variable Initializations
      MODNAM = 'PRECART'

C     READ in the Netid and Nettype
      IF (IFC .LT. 3) THEN
C        Missing Data Field
         GO TO 999
      END IF
      NETIDT = FIELD(3)
      IF (.NOT.NEWID .AND. (NETIDT.EQ.'    ' .OR.
     &    NETIDT.EQ.'XYINC' .OR. NETIDT.EQ.'XPNTS' .OR.
     &    NETIDT.EQ.'YPNTS' .OR. NETIDT.EQ.'ELEV' .OR.
     &    NETIDT.EQ.'HILL'  .OR.
     &    NETIDT.EQ.'FLAG'  .OR. NETIDT.EQ.'END')) THEN
         NETIDT = PNETID
         KTYPE = FIELD(3)
      ELSE IF (.NOT.NEWID .AND. NETIDT.EQ.PNETID) THEN
         KTYPE = FIELD(4)
      ELSE IF (NEWID .AND. NETIDT.NE.' ') THEN
         NEWID = .FALSE.
         KTYPE = FIELD(4)
C        The Keyword Counter
         NNET = NNET + 1
      ELSE
C        Invalid Secondary Keyword
         GO TO 999
      END IF

C     Start to Set Up the Network
      IF (KTYPE .EQ. 'STA') THEN
C        Initialize Logical Control Variables
         ISTA = .TRUE.
         IEND = .FALSE.
         NEWID = .FALSE.
         RECERR = .FALSE.
C        Set Counters of Calculation Field
         ICOUNT = 0
         JCOUNT = 0
      ELSE IF (KTYPE .EQ. 'XYINC') THEN
C        Set the Uniform Spacing Receptor Network           ---   CALL PREGENCAR
         CALL PREGENCAR
      ELSE IF (KTYPE.EQ.'XPNTS' .OR. KTYPE.EQ.'YPNTS') THEN
C        Set the Non-uniform Spacing Receptor Network       ---   CALL PREXYPNTS
         CALL PREXYPNTS
      ELSE IF (KTYPE .EQ. 'END') THEN
         IEND = .TRUE.
         IF (.NOT. RECERR) THEN
            NREC = NREC + ICOUNT*JCOUNT
         END IF
         ISTA = .FALSE.
         NEWID = .TRUE.

      ELSE IF (KTYPE.NE.'ELEV' .AND. KTYPE.NE.'FLAG' .AND.
     &         KTYPE.NE.'HILL') THEN
C        Invalid Secondary Keyword
         RECERR = .TRUE.
         GO TO 999

      END IF

      PNETID = NETIDT

 999  RETURN
      END

      SUBROUTINE PREGENCAR
C***********************************************************************
C                 PREGENCAR Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Generates Cartesian Grid Receptor Network With
C                 Uniform Spacing
C
C        PROGRAMMER: Roger Brode
C
C        DATE:    September 24, 1996
C
C        INPUTS:  Input Runstream Image Parameters
C
C        OUTPUTS: Cartesian Grid Receptor Network With Uniform
C                 Spacing
C
C        CALLED FROM:   PRECART
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: I, K
      REAL    :: TEMPR(6)
      DOUBLE PRECISION :: TEMPD(6), XDELTA, YDELTA
      LOGICAL ERROR

C     Variable Initializations
      MODNAM = 'PREGENCAR'
      ERROR = .FALSE.

C     Check for Location of Secondary Keyword, XYINC
      DO I = 1, IFC
         IF (FIELD(I) .EQ. 'XYINC') THEN
            ISC = I + 1
         END IF
      END DO

C     Determine Whether There Are Enough Parameter Fields
      IF (IFC .EQ. ISC-1) THEN
C        Missing Parameter
         RECERR = .TRUE.
         GO TO 999
      ELSE IF (IFC .GT. ISC+5) THEN
C        Too Many Parameters
         RECERR = .TRUE.
         GO TO 999
      ELSE IF (IFC .LT. ISC+5) THEN
C        Too Few Parameters
         RECERR = .TRUE.
         GO TO 999
      END IF

C     Input The Numerical Values
      DO K = 1,6
         IF (K .EQ. 2 .OR. K .EQ. 5) THEN
            CALL STONUM(FIELD(ISC + K-1),ILEN_FLD,TEMPR(K),IMIT)
C           Check The Numerical Field
            IF (IMIT .EQ. -1) THEN
               ERROR = .TRUE.
               RECERR = .TRUE.
            END IF
         ELSE
            CALL STODBL(FIELD(ISC + K-1),ILEN_FLD,TEMPD(K),IMIT)
C           Check The Numerical Field
            IF (IMIT .EQ. -1) THEN
               ERROR = .TRUE.
               RECERR = .TRUE.
            END IF
         END IF
      END DO

      IF (ERROR) THEN
         ERROR = .FALSE.
         GO TO 999
      END IF

C     Assign Values to Appropriate Variables for Generated Network
      XINT   = TEMPD(1)
      ICOUNT = NINT(TEMPR(2))
      XDELTA = TEMPD(3)
      YINT   = TEMPD(4)
      JCOUNT = NINT(TEMPR(5))
      YDELTA = TEMPD(6)

C     Assign Them to the Coordinate Arrays
      IF (ICOUNT .GT. IXM) THEN
         IXM = ICOUNT
      END IF
      IF (JCOUNT .GT. IYM) THEN
         IYM = JCOUNT
      END IF

 999  RETURN
      END

      SUBROUTINE PREXYPNTS
C***********************************************************************
C                 PREXYPNTS Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Processes Cartesian Grid x,y Input Value
C
C        PROGRAMMER: Roger Brode
C
C        DATE:    September 24, 1996
C
C        INPUTS:  Input Runstream Image Parameters
C
C        OUTPUTS: Cartesian Grid x,y Input Value
C
C        CALLED FROM:   PRECART
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: I, JSET

C     Variable Initializations
      MODNAM = 'PREXYPNTS'

      IF (KTYPE .EQ. 'XPNTS') THEN
C        Check for Location of Secondary Keyword, XPNTS
         DO I = 1, IFC
            IF (FIELD(I) .EQ. 'XPNTS') THEN
               ISC = I + 1
            END IF
         END DO

C        Determine Whether There Are Enough Parameter Fields
         IF (IFC .EQ. ISC-1) THEN
C           Missing Parameter
            RECERR = .TRUE.
            GO TO 999
         END IF

         ISET = ICOUNT
         DO I = ISC, IFC
            CALL STONUM(FIELD(I),ILEN_FLD,FNUM,IMIT)
C           Check The Numerical Field
            IF (IMIT .EQ. -1) THEN
               RECERR = .TRUE.
            END IF
            ISET = ISET + 1
            IF (ISET .GT. IXM) THEN
               IXM = ISET
            END IF
         END DO
         ICOUNT = ISET

      ELSE IF (KTYPE .EQ. 'YPNTS') THEN
C        Check for Location of Secondary Keyword, YPNTS
         DO I = 1, IFC
            IF (FIELD(I) .EQ. 'YPNTS') THEN
               ISC = I + 1
            END IF
         END DO

C        Determine Whether There Are Enough Parameter Fields
         IF (IFC .EQ. ISC-1) THEN
C           Missing Parameter
            RECERR = .TRUE.
            GO TO 999
         END IF

         JSET = JCOUNT
         DO I = ISC, IFC
            CALL STONUM(FIELD(I),ILEN_FLD,FNUM,IMIT)
C           Check The Numerical Field
            IF (IMIT .EQ. -1) THEN
               RECERR = .TRUE.
            END IF
            JSET = JSET + 1
            IF (JSET .GT. IYM) THEN
               IYM = JSET
            END IF
         END DO
         JCOUNT = JSET

      END IF

 999  RETURN
      END

      SUBROUTINE PREPOLR
C***********************************************************************
C                 PREPOLR Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Processes Polar Grid Receptor Network Inputs
C
C        PROGRAMMER:  Roger Brode
C
C        DATE:    September 24, 1996
C
C        INPUTS:  Input Runstream Image Parameters
C
C        OUTPUTS: Polar Receptor Network Inputs
C
C        CALLED FROM:   RECSIZ
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

C     Variable Initializations
      MODNAM = 'PREPOLR'

      IF (IFC .LT. 3) THEN
C        Missing Data Field
         GO TO 999
      END IF

C     READ in the Netid and Nettype
      NETIDT = FIELD(3)
      IF (.NOT.NEWID .AND. (NETIDT.EQ.'    ' .OR.
     &    NETIDT.EQ.'ORIG' .OR. NETIDT.EQ.'DIST' .OR.
     &    NETIDT.EQ.'DDIR' .OR. NETIDT.EQ.'ELEV' .OR.
     &    NETIDT.EQ.'HILL' .OR.
     &    NETIDT.EQ.'FLAG' .OR. NETIDT.EQ.'GDIR' .OR.
     &    NETIDT.EQ.'END')) THEN
         NETIDT = PNETID
         KTYPE = FIELD(3)
      ELSE IF (.NOT.NEWID .AND. NETIDT.EQ.PNETID) THEN
         KTYPE = FIELD(4)
      ELSE IF (NEWID .AND. NETIDT.NE.'    ') THEN
         NEWID = .FALSE.
         KTYPE = FIELD(4)
C        The Keyword Counter
         NNET = NNET + 1
      ELSE
C        Invalid Secondary Keyword
         RECERR = .TRUE.
         GO TO 999
      END IF

C     Start to Set Up the Network
      IF (KTYPE .EQ. 'STA') THEN
         ISTA = .TRUE.
         IEND = .FALSE.
         NEWID = .FALSE.
         RECERR = .FALSE.
         ICOUNT = 0
         JCOUNT = 0
      ELSE IF (KTYPE .EQ. 'DIST') THEN
C        Read in the Distance Set                           ---   CALL PREPOLDST
         CALL PREPOLDST
      ELSE IF (KTYPE .EQ. 'GDIR') THEN
         CALL PREGENPOL
      ELSE IF (KTYPE .EQ. 'DDIR') THEN
         CALL PRERADRNG
      ELSE IF (KTYPE .EQ. 'END') THEN
         IEND = .TRUE.
C        Get the Final Result
         IF (.NOT. RECERR) THEN
            NREC = NREC + ICOUNT*JCOUNT
         END IF
         ISTA = .FALSE.
         NEWID = .TRUE.

      ELSE IF (KTYPE.NE.'ELEV' .AND. KTYPE.NE.'FLAG' .AND.
     &         KTYPE.NE.'HILL' .AND. KTYPE.NE.'ORIG') THEN
C        Invalid Secondary Keyword
         RECERR = .TRUE.
         GO TO 999

      END IF

      PNETID = NETIDT

 999  RETURN
      END

      SUBROUTINE PREPOLDST
C***********************************************************************
C                 PREPOLDST Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Gets Distances for the Polar Network
C
C        PROGRAMMER: Roger Brode
C
C        DATE:    September 24, 1996
C
C        INPUTS:  Input Runstream Image Parameters
C
C        OUTPUTS: Polar Network Distance Input Value
C
C        CALLED FROM:   PREPOLR
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: I

C     Variable Initializations
      MODNAM = 'PREPOLDST'

C     Skip the Unrelated Fields
      DO I = 1, IFC
         IF (FIELD(I) .EQ. 'DIST') THEN
            ISC = I + 1
         END IF
      END DO

C     Determine Whether There Are Enough Parameter Fields
      IF (IFC .EQ. ISC-1) THEN
C        Missing Parameter
         RECERR = .TRUE.
         GO TO 999
      END IF

      ISET = ICOUNT

      DO I = ISC, IFC
         CALL STONUM(FIELD(I),ILEN_FLD,FNUM,IMIT)
C        Check The Numerical Field
         IF (IMIT .EQ. -1) THEN
            RECERR = .TRUE.
         END IF
         ISET = ISET + 1
         IF (ISET .GT. IXM) THEN
            IXM = ISET
         END IF
      END DO

      ICOUNT = ISET

 999  RETURN
      END

      SUBROUTINE PREGENPOL
C***********************************************************************
C                 PREGENPOL Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Generates Polar Receptor Network With
C                 Uniform Spacing
C
C        PROGRAMMER: Roger Brode
C
C        DATE:    September 24, 1996
C
C        INPUTS:  Input Runstream Image Parameters
C
C        OUTPUTS: Polar Receptor Network With Uniform Direction Spacing
C
C        CALLED FROM:   PREPOLR
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: I, K
      DOUBLE PRECISION :: TEMPP(3), DIRINI, DIRINC
      LOGICAL ERROR

C     Variable Initializations
      MODNAM = 'PREGENPOL'
      ERROR = .FALSE.

C     Check for the Location of the Secondary Keyword, GDIR
      DO I = 1, IFC
         IF (FIELD(I) .EQ. 'GDIR') THEN
            ISC = I + 1
         END IF
      END DO

C     Determine Whether There Are Enough Parameter Fields
      IF (IFC .EQ. ISC-1) THEN
C        Missing Parameter
         RECERR = .TRUE.
         GO TO 999
      ELSE IF (IFC .LT. ISC+2) THEN
C        Not Enough Parameters
         RECERR = .TRUE.
         GO TO 999
      ELSE IF (IFC .GT. ISC+2) THEN
C        Too Many Parameters
         RECERR = .TRUE.
         GO TO 999
      END IF

C     Input Numerical Values
      DO K = 1, 3
         CALL STODBL(FIELD(ISC + K-1),ILEN_FLD,TEMPP(K),IMIT)
C        Check The Numerical Field
         IF (IMIT .EQ. -1) THEN
            RECERR = .TRUE.
            ERROR = .TRUE.
         END IF
      END DO

      IF (ERROR) THEN
         ERROR = .FALSE.
         GO TO 999
      END IF

      JCOUNT = IDNINT(TEMPP(1))
      DIRINI = TEMPP(2)
      DIRINC = TEMPP(3)

C     Assign Them to the Coordinate Arrays
      IF (JCOUNT .GT. IYM) THEN
         IYM = JCOUNT
      END IF

 999  RETURN
      END

      SUBROUTINE PRERADRNG
C***********************************************************************
C                 PRERADRNG Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Processes Non-Uniform Polar Network Value
C
C        PROGRAMMER: Roger Brode
C
C        DATE:    September 24, 1996
C
C        INPUTS:  Input Runstream Image Parameters
C
C        OUTPUTS: Polar Network Directions in Non-Uniform Spacing
C
C        CALLED FROM:   PREPOLR
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: I

C     Variable Initializations
      MODNAM = 'PRERADRNG'

C     Skip the non-useful Fields
      DO I = 1, IFC
         IF (FIELD(I) .EQ. 'DDIR') THEN
            ISC = I + 1
         END IF
      END DO

C     Determine Whether There Are Enough Parameter Fields
      IF (IFC .EQ. ISC-1) THEN
C        Error Message: Missing Parameter
         RECERR = .TRUE.
         GO TO 999
      END IF

      ISET = JCOUNT

      DO I = ISC, IFC
         CALL STONUM(FIELD(I),ILEN_FLD,FNUM,IMIT)
C        Check The Numerical Field
         IF (IMIT .EQ. -1) THEN
            RECERR = .TRUE.
         END IF
         ISET = ISET + 1
         IF (ISET .GT. IYM) THEN
            IYM = ISET
         END IF
      END DO

      JCOUNT = ISET

 999  RETURN
      END

      SUBROUTINE SET_WINDOW
C***********************************************************************
C                 SET_WINDOW Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Preprocess Meteorology Surface Data Card (SURFDATA)
C                 to Set Date Window for Y2K Fixes
C
C        PROGRAMMER: Roger Brode, PES, Inc.
C
C        DATE:    April 29, 1999
C
C        MODIFICATIONS:
C
C                    To subtract 1 from ISTRT_WIND in case data file
C                    contains data from end of previous year.
C                    R.W. Brode, PES, Inc.  8/28/01
C
C        INPUTS:  Input Runstream Image Parameters
C
C        OUTPUTS: Starting Century, ISTRT_CENT                    [I4]
C                 Starting Year for 2-digit Window, ISTRT_WIND    [I4]
C
C        ERROR HANDLING:   Checks for Too Few Parameters;
C                          Checks for Invalid Numeric Fields;
C                          Checks for Too Many Parameters
C
C        CALLED FROM:   PRESET
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

C     Variable Initializations
      MODNAM = 'SET_WINDOW'

      IF (IFC .LT. 4) THEN
         GO TO 999
      ELSE IF (IFC .GT. 7) THEN
         GO TO 999
      END IF

      CALL STONUM(FIELD(4),ILEN_FLD,FNUM,IMIT)
C     Check The Numerical Field
      IF (IMIT .NE. 1) THEN
         GO TO 999
      END IF
      ISYEAR = NINT(FNUM)
      IF (ISYEAR .LT. 100) THEN
C        Write warning message for 2-digit year, and set default "windowing"
C        variables, ISTRT_CENT (=19) and ISTRT_WIND (=50).
         IF (.NOT. L_SkipMessages) THEN
            CALL ERRHDL(PATH,MODNAM,'W','360',KEYWRD)
         END IF
         ISTRT_CENT = 19
         ISTRT_WIND = 50
      ELSE
C        Determine starting century (ISTRT_CENT) and starting year for
C        window (ISTRT_WIND) from 4-digit input
         ISTRT_CENT = ISYEAR/100
         ISTRT_WIND = ISYEAR - ISTRT_CENT*100
C        Subtract 1 from ISTRT_WIND in case data file contains data
C        from end of previous year
         ISTRT_WIND = ISTRT_WIND - 1
         IF (ISTRT_WIND .LT. 0) THEN
            ISTRT_WIND = 0
         END IF
C        Check for year .ge. 2148 to avoid integer overflow on FULLDATE
         IF (ISTRT_CENT .GE. 21 .AND. ISTRT_WIND .GE. 48) THEN
            CALL ERRHDL(PATH,MODNAM,'E','365',KEYWRD)
            ISTRT_CENT = 21
            ISTRT_WIND = 47
         END IF
      END IF

      GO TO 1000

 999  CONTINUE
C     For error in processing assume 1900 for start century and 50 for window
      ISTRT_CENT = 19
      ISTRT_WIND = 50

 1000 RETURN
      END

      SUBROUTINE CHK_ENDYR
C***********************************************************************
C                 CHK_ENDYR Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Checks date for "end-of-year" for use in ANNUAL
C                 averages and PM-2.5 processing.
C
C        PROGRAMMER: Roger Brode
C
C        DATE:
C
C        MODIFIED:   To allow user-specified rank for PM2.5 processing
C                    to accommodate latest guidance for PM2.5 modeling.
C                    R.W. Brode, U.S. EPA/OAQPS/AQMG, 02/28/2011
C               
C                    To check for allocation status prior to
C                    initializing allocatable arrays, and include
C                    maximum annual value arrays.
C                    R.W. Brode, U.S. EPA/OAQPS/AQMG, 10/19/2009
C
C        INPUTS:  Input Runstream Image Parameters
C
C        OUTPUTS: Plant Boundary Receptor Location Inputs
C
C        CALLED FROM:   HRLOOP
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      IMPLICIT NONE
      CHARACTER MODNAM*12

      INTEGER :: IEND_DAY, N

C     Variable Initializations
      MODNAM = 'CHK_ENDYR'

      IF( (IENDMN.EQ.2.AND.IENDDY.EQ.29.AND.IMONTH.EQ.2) .AND.
     &    (MOD(IYR,4).NE.0) .OR.
     &    (MOD(IYR,100).EQ.0 .AND. MOD(IYR,400).NE.0)) THEN
C        Set End Day to 28 for non-leap year February
         IEND_DAY = 28
      ELSE
         IEND_DAY = IENDDY
      END IF

      IF (IMONTH.EQ.IENDMN .AND. IDAY.EQ.IEND_DAY .AND.
     &    IHOUR.EQ.IENDHOUR) THEN
C        End of year reached, increment counter and store High-N-High (HNH) values
         NUMYRS = NUMYRS + 1
C        Reset hour counter for MAXDCONT met data arrays
         IHR_NDX = 0
         IF (ANNUAL) THEN
C ---       Calculate ANNUAL averages
            CALL PERAVE
C ---       Check for ANNUAL POSTFILE
            IF( ANPOST )THEN
               CALL PSTANN
            ENDIF
C ---       Sum the annual averages
            SUMANN(:,:,:) = SUMANN(:,:,:) + ANNVAL(:,:,:)
C           Re-initialize the annual counters and array
            IANHRS  = 0
            IANCLM  = 0
            IANMSG  = 0
            NSKIPTOT = 0
            IF (ALLOCATED(ANNVAL))  ANNVAL  = 0.0D0
            IF (ALLOCATED(AMXVAL))  AMXVAL  = 0.0D0
            IF (ALLOCATED(IMXLOC))  IMXLOC  = 0
         END IF
         IF ((PM25AVE .OR. NO2AVE .OR. SO2AVE) .AND. NUMAVE.GE.1) THEN
C ---       Sum the High-N-High 24-hour values for PM-2.5, 
C           or High-N-High 1-hour values for NO2
            DO N = 1, NVAL
               SUMHNH(1:NUMREC,1:NUMGRP,N) = 
     &         SUMHNH(1:NUMREC,1:NUMGRP,N) + 
     &            HIMXDLY(1:NUMREC,1:NUMGRP,N)
               IF (NUMYRS .LE. NYEARS) THEN
                  HIMXDLY_BYYR(1:NUMREC,1:NUMGRP,N,NUMYRS) =
     &                 HIMXDLY(1:NUMREC,1:NUMGRP,N)
                  NHIDATMXD_BYYR(1:NUMREC,1:NUMGRP,N,NUMYRS) =
     &                 NHIDATMXD(1:NUMREC,1:NUMGRP,N)
               ELSE
C ---             Write Error Message        ! Too many years
                  WRITE(DUMMY,'(''NYR='',I4)') NYEARS
                  CALL ERRHDL(PATH,MODNAM,'E','482',DUMMY)
                  RUNERR = .TRUE.
               END IF
               
               IF (MXDAILY_BYYR .AND. NHIAVE(N,1) .EQ. 1) THEN
                  CALL MXDYBYYR(N)
               ENDIF
            END DO
C           Re-initialize the MAXDAILY Value Arrays used for
C           PM25/NO2/SO2 Processing averaged across years
            IF (ALLOCATED(HIMXDLY))   HIMXDLY   = 0.0D0
            IF (ALLOCATED(NHIDATMXD)) NHIDATMXD = 0
         END IF
         NREMAIN = 0
      ELSE
C        Increment counter for number of hours remaining after
C        the end of the last year
         NREMAIN = NREMAIN + 1
      END IF

      RETURN
      END

      SUBROUTINE MAXDCONT_LOOP
C***********************************************************************
C                MAXDCONT_LOOP Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:    Control "post-processing" for OU MAXDCONT option
C
C        PROGRAMMER: Roger W. Brode
C
C        DATE:       February 28, 2011
C
C        MODIFIED:   Modified to add 'saved' arrays for buoyant line source
C                    Amec Foster Wheeler, 06/30/2015
C
C        MODIFIED:   Modified subroutine MAXDCONT_LOOP to correct problems 
C                    with the MAXDCONT option for applications that vary 
C                    emissions (EMISFACT), background ozone data (O3VALUES), 
C                    or background concentrations (BACKGRND) with a day-of-week 
C                    component, e.g., SHRDOW or SHRDOW7, etc. Also modified 
C                    subroutine MAXDCONT_LOOP to include checks on the 
C                    consistency between results in the SUMVAL_MAXD and 
C                    SUMHNH arrays for the "target" source group under 
C                    the MAXDCONT option. 
C                    R.W. Brode, U.S. EPA/OAQPS/AQMG, 02/29/2012
C
C        INPUTS:  
C
C        OUTPUTS: 
C
C        CALLED FROM:   MAIN
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      USE BUOYANT_LINE
      IMPLICIT NONE
      CHARACTER MODNAM*12
      INTEGER I

      INTEGER :: NDAY(12)
      INTEGER :: INYR, INMN, INDY, IJDY, JDY
      INTEGER :: ISDATE_SAV, IEDATE_SAV
      INTEGER :: JGRP, IVAL, IR, IG
      INTEGER :: ICYR, ICYR2, IMN, IDY, IHR, ICJDAY, IPJDAY
      LOGICAL :: HIT_THRESH

C     Variable Initializations
      DATA NDAY/0,31,59,90,120,151,181,212,243,273,304,334/

C     Variable Initializations

      MODNAM = 'MAXDCONT_LOOP'
      PATH   = 'CN'

      HIT_THRESH = .FALSE.

C --- Reinitialize NUMHRS, NUMCLM and NUMMSG
      NUMHRS(:) = 0
      NUMCLM(:) = 0
      NUMMSG(:) = 0
C     Initialize __VAL arrays (1:NUMTYP)
      HRVAL   = 0.0D0
      AERVAL  = 0.0D0
      AVEVAL  = 0.0D0
      IF (ALLOCATED(BACKAVE)) BACKAVE(:) = 0.0D0
      IF (ALLOCATED(CHI)) CHI(:,:,:) = 0.0D0

C --- Copy standard receptor arrays to "saved" arrays
      AXR_SAV = AXR
      AYR_SAV = AYR
      AZELEV_SAV = AZELEV
      AZHILL_SAV = AZHILL
      AZFLAG_SAV = AZFLAG
C     For buoyant line sources - results of first tranlation/rotation
      XR_SCS_SAV = XR_SCS
      YR_SCS_SAV = YR_SCS
      
C --- Save original start and end dates (ISDATE and IEDATE)
C     for processing messages in TERRST; these variables
C     include 4-digit year (for comparisons to FULLDATE)
      ISDATE_SAV = ISDATE
      IEDATE_SAV = IEDATE
      
C --- Reset number of receptors (NUMREC) to 1
C     for use in max daily contribution analyses
      NUMREC = 1
C --- Initialized SUMVAL_MAXD array for max daily contributions
      SUMVAL_MAXD = 0.0D0

C --- Set logical flag to skip messages while re-processing
C     the meteorological, BACKGRND, and ozone data
      L_SkipMessages = .TRUE.

C --- Loop through "target" source groups for 
C     max daily contributions option (OU MAXDCONT)
      DO JGRP = 1, NUMGRP
      
         IF (MAXDCONT(JGRP) .GT. 0) THEN
C ---       Max daily contribution results selected
C           for this source group

C ---       Loop through user-specified ranks for max daily
C           contribution analysis for this source group
            DO IVAL = MXD_RANK(JGRP,1), MXD_RANK(JGRP,2)

C ---          Write message to screen indicating current Rank
               WRITE(*,909) GRPID(JGRP), IVAL
 909           FORMAT
     &('+','Now Processing MAXDCONT for Group ',A8,' and Rank No. ',I4)
                  
               DO IR = 1, NREC
C ---             Loop through receptors; but skip receptor if
C                 value is below user-specified threshold
                  IF (SUMHNH(IR,JGRP,IVAL) .LE. 0.0D0 .OR.
     &                SUMHNH(IR,JGRP,IVAL) .LT. MAXD_THRESH(JGRP)) CYCLE
                  
C ---             Assign data from the "saved" arrays 
C                 to array index 1 for standard arrays
                  AXR(1) = AXR_SAV(IR) 
                  AYR(1) = AYR_SAV(IR) 
                  AZELEV(1) = AZELEV_SAV(IR) 
                  AZHILL(1) = AZHILL_SAV(IR) 
                  AZFLAG(1) = AZFLAG_SAV(IR) 
                                 
                  IF (L_BLSOURCE) THEN
                     XR_SCS(1) = XR_SCS_SAV(IR)
                     YR_SCS(1) = YR_SCS_SAV(IR)
                  END IF
C ---             Loop through number of years processed for
C                 max daily contribution analysis
                  DO I = 1, NUMYRS
C ---                Assign "start date" and "end date" for
C                    max daily 1-hr value associated with
C                    this rank
                     IEDATE = NHIDATMXD_BYYR(IR,JGRP,IVAL,I)
C                    Convert start date from 8-digits to 10-digits
                     ICYR2  = IEDATE/1000000
                     IF (ICYR2 .GE. ISTRT_WIND .AND. 
     &                                    ICYR2 .LE. 99) THEN
                        ICYR   = ISTRT_CENT*100 + ICYR2
                        IEDATE = ISTRT_CENT*100000000 + 
     &                           NHIDATMXD_BYYR(IR,JGRP,IVAL,I)
                     ELSE IF (ICYR2 .LT. ISTRT_WIND) THEN
                        ICYR   = (ISTRT_CENT+1)*100 + ICYR2
                        IEDATE = (ISTRT_CENT+1)*100000000 +
     &                           NHIDATMXD_BYYR(IR,JGRP,IVAL,I)
                     END IF
                     
                     IF (NO2AVE .OR. SO2AVE) THEN
C ---                   Assign start date to end date for NO2 or SO2
C                       since we're only processing 1 hour at a time
                        ISDATE = IEDATE
                     ELSE IF (PM25AVE) THEN
C ---                   Subtract 23 from end date for PM2.5 since
C                       these will always be 24-hour averages
                        ISDATE = IEDATE - 23
                     END IF

                     DO FULLDATE = ISDATE, IEDATE
C ---                   Calculate IHR_NDX, hour-of-year array index 
C                       to extract met data for MAXDCONT option,
C                       based on FULLDATE.
                        IMN = (FULLDATE/10000) - 
     &                        (FULLDATE/1000000)*100
                        IDY = (FULLDATE/100) - 
     &                        (FULLDATE/10000)*100
                        IHR =  FULLDATE - (FULLDATE/100)*100
                        CALL JULIAN(ICYR,IMN,IDY,IJDY)

C ---                   Determine julian day for current year based on 
C                       orignial start month/day
                        CALL JULIAN(ICYR,ISMN,ISDY,ICJDAY)
                        
C ---                   Check for invalid Julian days (=0), indicating a
C                       potential problem with the NHIDATMXD_BYYR array.
C                       An error message will already have been issued in
C                       subroutine Julian, but exit loop to avoid Fortran
C                       runtime errors since IHR_NDX may not be valid.
                        IF (IJDY .EQ. 0 .OR. ICJDAY .EQ. 0) THEN
                           EXIT
                        END IF

                        IF (ISJDAY .EQ. 1) THEN
C ---                      Data starts on Jan. 1
                           IHR_NDX = 24*(IJDY-ISJDAY)+(IHR-ISHR)+1
                        ELSE IF (IJDY .GE. ICJDAY) THEN
C ---                      Data does not start on Jan. 1, but "event"
C                          jday is .ge. start jday
                           IHR_NDX = 24*(IJDY-ICJDAY)+(IHR-ISHR)+1
                        ELSE
C ---                      Data does not start on Jan. 1, and "event"
C                          jday is .lt. start jday; calculation of the
C                          "hour-of-year" index must account for potential
C                          influence of leap year.
C ---                      Determine julian day for previous year based on 
C                          start month/day
                           CALL JULIAN(ICYR-1,ISMN,ISDY,IPJDAY)
                           IF (IPJDAY .GT. ICJDAY) THEN
C ---                        Account for leap year
                             IHR_NDX = 24*(IJDY-IPJDAY+366)+(IHR-ISHR)+1
                           ELSE
C ---                        No leap year adjustment needed
                             IHR_NDX = 24*(IJDY-IPJDAY+365)+(IHR-ISHR)+1
                           END IF
                        END IF
C ---                   Assign year arrary index for extracting met data 
                        IYR_NDX = I

C ---                   Set other global date variables for this hour                        
                        IMONTH = IMN
                        IDAY   = IDY
                        IHOUR  = IHR
                        IYEAR  = ICYR2    ! 2-digit year
                        IYR    = ICYR     ! 4-digit year
                        KURDAT = ICYR2*1000000 + IMONTH*10000 + 
     &                           IDAY*100 + IHOUR
                        KURPFL = KURDAT
                     
C ---                   Call MAXDCALC subroutine to calculate
C                       max daily contributions for this hour
                        CALL MAXDCALC
                     END DO

                     DO IG = 1, NUMGRP
C ---                   Loop through source groups to determine
C                       contributions to target source group
                        SUMVAL_MAXD(IVAL,IG,JGRP,IR) = 
     &                  SUMVAL_MAXD(IVAL,IG,JGRP,IR) + 
     &                              AVEVAL(1,IG,1,1)
                     END DO
                     
C ---                Reinitialize AVEVAL array
                     AVEVAL(:,:,:,:) = 0.0D0
C ---                Reinitialize NUMHRS, NUMCLM and NUMMSG
                     NUMHRS(:) = 0
                     NUMCLM(:) = 0
                     NUMMSG(:) = 0
C                    Initialize __VAL arrays (1:NUMTYP)
                     HRVAL   = 0.0D0
                     AERVAL  = 0.0D0
C                             
                     IF (ALLOCATED(CHI)) CHI(:,:,:) = 0.0D0

                  END DO  ! end loop on years
               
                  DO IG = 1, NUMGRP
C ---                Divide sumval_maxd results by number of years
C                    to get averages across number of years modeled
                     SUMVAL_MAXD(IVAL,IG,JGRP,IR) = 
     &               SUMVAL_MAXD(IVAL,IG,JGRP,IR)/
     &                               DBLE(NUMYRS)
                  END DO

C ---             Check for consistency of SUMVAL_MAXD and SUMHNH arrays for
C                 "target" source group under MAXDCONT option
                  IF (SUMHNH(IR,JGRP,IVAL) .GT. 0.0D0) THEN
                     IF( DABS( (SUMVAL_MAXD(IVAL,JGRP,JGRP,IR) -
     &                          SUMHNH(IR,JGRP,IVAL)) )/
     &                          SUMHNH(IR,JGRP,IVAL) .GT. 5.0D-6) THEN
C ---                   Arrays don't match; issue warning message indicating
C                       potential coding error. A warning to the default 
C                       output unit will also be issued at the end of the run.
                        WRITE(DUMMY,'("G",I3.3,"R",I2.2,"R",I4.4)')
     &                        MIN(JGRP,999), MIN(IVAL,99), MIN(IR,9999)
                        CALL ERRHDL(PATH,MODNAM,'W','498',DUMMY)
                        L_MAXDCONT_OrigConc_Warning = .TRUE.
                     END IF
                  END IF
C              
               END DO   ! end loop over receptors

C ---          Call subroutine to write MAXDCONT output file
C              for this source group and rank
               CALL MAXDCNT_FILE(JGRP,IVAL)

C ---          Check for value below threshold for MAXDCONT
               IF (MAXVAL(SUMHNH(1:NREC,JGRP,IVAL)) .LT. 
     &                                  MAXD_THRESH(JGRP)) THEN
C ---             All values for this rank are below threshold;
C                 reset upper bound of ranks to process, set flag
C                 to indicate that threshold was reached, and 
C                 EXIT the loop over ranks
                  MXD_RANK(JGRP,2) = IVAL
                  HIT_THRESH = .TRUE.
                  EXIT
               END IF

            END DO   ! end loop over ranks

C ---       Check for whether MAXD_THRESH specified by user was
C           not reached
            IF (MAXD_THRESH(JGRP) .GT. 0.0D0 .AND. 
     &                                    .NOT.HIT_THRESH) THEN
C ---          User-specified threshold was not reached within the
C              range of ranks analyzed, based on the RECTABLE keyword;
C              issue warning message
               WRITE(DUMMY,'(''GRP '',A8)') GRPID(JGRP)
               CALL ERRHDL(PATH,MODNAM,'W','415',DUMMY)
            END IF
C ---       Reset HIT_THRESH flag for next source group
            HIT_THRESH = .FALSE.
               
         END IF
               
      END DO    ! end loop over source groups

C --- Copy saved receptor arrays to standard arrays
      AXR = AXR_SAV 
      AYR = AYR_SAV 
      AZELEV = AZELEV_SAV 
      AZHILL = AZHILL_SAV 
      AZFLAG = AZFLAG_SAV 

      IF (L_BLSOURCE) THEN
         XR_SCS = XR_SCS_SAV
         YR_SCS = YR_SCS_SAV
      END IF

C --- Reset number of receptors
      NUMREC = NREC
C --- Reset original start and end dates
      ISDATE = ISDATE_SAV
      IEDATE = IEDATE_SAV

      RETURN
      END

      SUBROUTINE MAXD_METEXT
C***********************************************************************
C                MAXD_METEXT Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE: Controls extraction of meteorological data from 
C                 arrays for use with OU MAXDCONT option
C
C        PROGRAMMER: Roger W. Brode
C
C        DATE:       February 28, 2011
C
C        MODIFIED:   Modified to PG stability parameter for the buoyant 
C                    line source
C                    Amec Foster Wheeler, 06/30/2015
C
C        INPUTS:  
C
C        OUTPUTS: 
C
C        CALLED FROM:   MAXDCALC
C***********************************************************************

C     Variable Declarations
      USE MAIN1
      USE BUOYANT_LINE
      IMPLICIT NONE
      CHARACTER MODNAM*12
      INTEGER I

C     Variable Initializations
      MODNAM = 'MAXD_METEXT'
      PATH   = 'MX'

C     Save Value of Last YR/MN/DY/HR and Previous Hour
      IPDATE = KURDAT
      IPYEAR = IYR
      IPHOUR = IHOUR

C     Set Meteorological Variables for This Hour
      SFCHF  = ASFCHF(IHR_NDX,IYR_NDX)
      UREF   = AUREF(IHR_NDX,IYR_NDX)
      UREFHT = AUREFHT(IHR_NDX,IYR_NDX)
      TA     = ATA(IHR_NDX,IYR_NDX)
      TREFHT = ATREFHT(IHR_NDX,IYR_NDX)
      WDREF  = AWDREF(IHR_NDX,IYR_NDX)
      USTAR  = AUSTAR(IHR_NDX,IYR_NDX)
      WSTAR  = AWSTAR(IHR_NDX,IYR_NDX)
      ZICONV = AZICONV(IHR_NDX,IYR_NDX)
      ZIMECH = AZIMECH(IHR_NDX,IYR_NDX)
      OBULEN = AOBULEN(IHR_NDX,IYR_NDX)
      VPTGZI = AVPTGZI(IHR_NDX,IYR_NDX)
      SFCZ0  = ASFCZ0(IHR_NDX,IYR_NDX)
      KST    = NINT( AKST(IHR_NDX,IYR_NDX) )
      BLTA   = ABLTA(IHR_NDX,IYR_NDX)
      
      IF (LDGAS .OR. LDPART .OR. LWPART .OR. LWGAS) THEN
         BOWEN  = ABOWEN(IHR_NDX,IYR_NDX)
         ALBEDO = AALBEDO(IHR_NDX,IYR_NDX)
         IPCODE = IAPCODE(IHR_NDX,IYR_NDX)
         PRATE  = APRATE(IHR_NDX,IYR_NDX)
         RH     = ARH(IHR_NDX,IYR_NDX)
         SFCP   = ASFCP(IHR_NDX,IYR_NDX)
         NCLOUD = NACLOUD(IHR_NDX,IYR_NDX)
         QSW    = AQSW(IHR_NDX,IYR_NDX)
         Wnew   = AWnew(IHR_NDX,IYR_NDX)
         f2     = Af2(IHR_NDX,IYR_NDX)
         EsTa   = AEsTa(IHR_NDX,IYR_NDX)
         Prec1  = APrec1(IHR_NDX,IYR_NDX)
         Prec2  = APrec2(IHR_NDX,IYR_NDX)
      END IF
      
      RURUSTR   = ARURUSTR(IHR_NDX,IYR_NDX)  
      RUROBULEN = ARUROBULEN(IHR_NDX,IYR_NDX)

      CLMHR = ACLMHR(IHR_NDX,IYR_NDX)
      MSGHR = AMSGHR(IHR_NDX,IYR_NDX)
      
      UNSTAB  = AUNSTAB(IHR_NDX,IYR_NDX)
      STABLE  = ASTABLE(IHR_NDX,IYR_NDX)
      URBSTAB = AURBSTAB(IHR_NDX,IYR_NDX)

      NDX4ZI = ANDX4ZI(IHR_NDX,IYR_NDX)
      UATZI  = AUATZI(IHR_NDX,IYR_NDX)  
      SVATZI = ASVATZI(IHR_NDX,IYR_NDX) 
      SWATZI = ASWATZI(IHR_NDX,IYR_NDX) 
      UAVG   = AUAVG(IHR_NDX,IYR_NDX)  
      SVAVG  = ASVAVG(IHR_NDX,IYR_NDX) 
      SWAVG  = ASWAVG(IHR_NDX,IYR_NDX) 
      PTATZI = APTATZI(IHR_NDX,IYR_NDX) 

      GRIDWD(1:MXGLVL)  = AGRIDWD(IHR_NDX,1:MXGLVL,IYR_NDX)
      GRIDWS(1:MXGLVL)  = AGRIDWS(IHR_NDX,1:MXGLVL,IYR_NDX) 
      GRIDSW(1:MXGLVL)  = AGRIDSW(IHR_NDX,1:MXGLVL,IYR_NDX)
      GRIDSV(1:MXGLVL)  = AGRIDSV(IHR_NDX,1:MXGLVL,IYR_NDX)
      GRIDTG(1:MXGLVL)  = AGRIDTG(IHR_NDX,1:MXGLVL,IYR_NDX)
      GRIDPT(1:MXGLVL)  = AGRIDPT(IHR_NDX,1:MXGLVL,IYR_NDX)
      IF (NSEC .GT. 0) THEN
         GRIDRHO(1:MXGLVL) = AGRIDRHO(IHR_NDX,1:MXGLVL,IYR_NDX)
      END IF
      IF (PVMRM .OR. PVMRM2) THEN
         GRIDEPS(1:MXGLVL) = AGRIDEPS(IHR_NDX,1:MXGLVL,IYR_NDX)
      END IF
      
      IF (NURB .GT. 0) THEN
         GRDSWR(1:MXGLVL) = AGRDSWR(IHR_NDX,1:MXGLVL,IYR_NDX) 
         GRDSVR(1:MXGLVL) = AGRDSVR(IHR_NDX,1:MXGLVL,IYR_NDX)
         GRDTGR(1:MXGLVL) = AGRDTGR(IHR_NDX,1:MXGLVL,IYR_NDX)
         GRDPTR(1:MXGLVL) = AGRDPTR(IHR_NDX,1:MXGLVL,IYR_NDX)
      
         DO I = 1, NURB
            GRDSWU(1:MXGLVL,I) = AGRDSWU(IHR_NDX,1:MXGLVL,IYR_NDX,I) 
            GRDSVU(1:MXGLVL,I) = AGRDSVU(IHR_NDX,1:MXGLVL,IYR_NDX,I)
            GRDTGU(1:MXGLVL,I) = AGRDTGU(IHR_NDX,1:MXGLVL,IYR_NDX,I) 
            GRDPTU(1:MXGLVL,I) = AGRDPTU(IHR_NDX,1:MXGLVL,IYR_NDX,I)
            ZIURB(I)           = AZIURB(IHR_NDX,IYR_NDX,I)    
            URBWSTR(I)         = AURBWSTR(IHR_NDX,IYR_NDX,I)  
            URBUSTR(I)         = AURBUSTR(IHR_NDX,IYR_NDX,I)  
            URBOBULEN(I)       = AURBOBULEN(IHR_NDX,IYR_NDX,I)
            L_MorningTrans(I)  = AL_MorningTrans(IHR_NDX,IYR_NDX,I)
         END DO
      END IF
      
      IF (.NOT.CLMHR .AND. .NOT.MSGHR) THEN
C ---    Set Meteorological Variables for Current Hour
         CALL SET_METDATA
      END IF

 999  RETURN
      END

C=======================================================================
C
C-----  Summaries of AERMOD model changes prior to version 13350 are 
C       provided below:
C       
C=======================================================================
C
C       This revised version of AERMOD (dated 12345) includes the
C       following modifications relative to the previous version
C       (dated 12060); see MCB#8 and AERMOD User's Guide Addendum.
C
C-----  Bug Fixes:
C
C       1.  Adjustments to wind speeds based on the assumption that input
C           wind speeds are vector (or resultant) mean winds have been 
C           removed (see Eq. 112 on page 79 of the AERMOD Model Formulation 
C           Document, EPA-454/R-03-004). This is considered a formulation bug
C           fix since current EPA guidance for site-specific meteorological 
C           monitoring (EPA-454/R-99-005) recommends that scalar mean wind 
C           speeds be used in steady-state Gaussian dispersion models. 
C           Furthermore, all wind speeds derived from NWS or FAA airport 
C           data represent scalar mean wind speeds.  
C       
C           An option has also been included on the MODELOPT keyword 
C           on the CO pathway (VECTORWS) that allows users to specify that 
C           input wind speeds are vector means, in which case the previous 
c           adjustments will be included. The new VECTORWS option is not 
C           linked with the DFAULT option, but users should be able to 
c           confirm that input speeds are vector means in order to justify 
c           use of the option.  However, scalar mean wind speeds, if 
C           available, should be used based on the current guidance.  
C           This change affected several subroutines.
C    
C       2.  Modified subroutine HRLOOP to correct a problem associated
C           with missing hourly ozone data during MAXDCONT post-processing
C           when no O3 values have been specified for substitution through
C           the OZONEVAL or O3VALUES keywords on the CO pathway. The 
C           MAXDCONT results in those cases were invalid and may have 
C           shown up as negative concentrations.
C    
C       3.  Modified subroutines DAYRNG, METEXT, and SET_DATES to correct a 
C           bug associated with use of the DAYRANGE keyword for multiple years
C           of meteorological data, where the YR/MN/DY conversion to Julian
C           may be incorrect.  Also included checks of the year specified on
C           the SURFDATA keyword with the first year included in the surface
C           meteorological data file. The ISYEAR variable is now adjusted to 
C           match the first year of the data file, unless the DAYRANGE keyword 
C           is being used. In that case, a fatal error is generated since the 
C           ISYEAR variable is used to determine leap vs. non-leap years for 
C           the YR/MN/DY conversion to Julian day for the DAYRANGE option.
C           
C       4.  Modified subroutine HRQREAD to check for large negative hourly 
C           emissions (< -90), which may be used as missing indicators.
C           Since AERMOD allows inputs of negative emissions for use in 
C           emission credit calculations, negative values used as missing 
C           indicators in the HOUREMIS file resulted in negative hourly 
C           concentrations in the previous versions.  Warning messages are 
C           generated and the emission rate is set to zero (0) for these cases.
C           
C       5.  Modified subroutine EVALFL to address a potential problem 
C           with URBAN applications where the L_MorningTrans logical
C           variable was not defined.
C
C       6.  Modified MAIN program and subroutine PRESET to check for 
C           duplicate STARTING keywords on the SO or RE pathways,since
C           that would reset the array limits for setup arrays to zero
C           during the PRESET phase, resulting in array subscript 
C           out-of-bounds runtime errors.
C
C
C-----  Enhancements:
C    
C       1.  Incorporated two new BETA (non-Default) options to address
C           concerns regarding model performance under low wind speed
C           conditions.  The LOWWIND1 option increases the minimum value
C           of sigma-v from 0.2 to 0.5 m/s and "turns off" the horizontal
C           meander component. The LOWWIND2 option increases the minimum
C           value of sigma-v from 0.2 to 0.3 m/s, and incorporates the 
C           meander component, with some adjustments to the algorithm, 
C           including an upper limit on the meander factor (FRAN) of 0.95.
C           A new LOW_WIND keyword has been added to the CO pathway that
C           allows users to adjust the minimum sigma-v value (within a
C           range of 0.01 to 1.0 m/s), and the minimum wind speed value
C           (within a range from 0.01 to 1.0 m/s), with a default value
C           of 0.2828 m/s, consistent with the default applied in previous
C           versions based on SQRT(2*SVmin*SVmin) with SVmin=0.2. The new
C           LOW_WIND keyword also allows users to adjust the maximum value
C           for the meander factor (FRAN) within a range of 0.50 to 1.0,
C           inclusive, when the LOWWIND2 option is used.  Both of the new 
C           LowWind BETA options also modify the the adjustment of vector 
C           mean wind speeds (based on Eq. 112, p. 79, of the AERMOD Model 
C           Formulation Document) to use the original values of sigma-v
C           before they are adjusted based on SVmin. The new LOWWIND1
C           and LOWWIND2 options are mutually exclusive and the model will
C           issue a fatal warning message if both options are specified.
C   
C           These new LowWind BETA options can also be used in conjunction 
C           with the new option in AERMET (v12345) to adjust the surface 
C           friction velocity (U*) under low-wind/stable conditions (ADJ_U*), 
C           based on Qian, W., and A. Venkatram, 2011: "Performance of Steady-
C           State Dispersion Models Under Low Wind-Speed Conditions", Boundary-
C           Layer Meteorology, 138, 475-491. The ADJ_U* option in AERMET is 
C           also considered a non-Default option, and users must also specify 
C           the BETA option on the MODELOPT keyword in AERMOD in order to 
C           process meteorological inputs derived using the new ADJ_U* option. 
C           See the AERMET Model Change Bulletin (MCB) #3 and AERMET User's 
C           Guide Addendum for additional details regarding the ADJ_U* option 
C           in AERMET.
C           
C           The new LowWind BETA options in AERMOD and the new ADJ_U*
C           option in AERMET are considered to be non-Default options
C           and are therefore subject to the alternative model provisions
C           in Section 3.2 of Appendix W (40 CFR Part 51). Users should
C           coordinate with the appropriate reviewing authority regarding
C           the procedures and requirements for approval of these BETA
C           options for regulatory modeling applications.
C          
C       2.  A new LINE source type has been included that allows users
C           to specify line-type sources based on a start-point and 
C           end-point of the line and the width of the line, as an 
C           alternative to the current AREA source type for rectangular
C           sources. The LINE source type utilizes the same routines as
C           the AREA source type, and will give identical results for 
C           equivalent source inputs.  The LINE source type also includes
C           an optional initial sigma-z parameter to account for initial
C           dilution of the emissions. As with the AREA source type, the
C           LINE source type does not include the horizontal meander 
C           component in AERMOD. This change affected several subroutines
C           and also includes a new LPARM subroutine.
C
C
C-----  Miscellaneous:
C
C       1.  Modified subroutine HRQEXT to include additional range checks
C           on hourly stack exit velocities and exit temperatures input
C           through the HOUREMIS option. A fatal error is generated if
C           the hourly exit temperature is less than 200K (about -100F)
C           and the exit velocity is greater than 200m/s, unless the hourly 
C           emissions are zero (0) for that hour. This may indicate that 
C           incorrect units for exit temperature have been used, or that 
C           the order of exit temperature and exit velocity may have been 
C           reversed in the HOUREMIS file. Also added a new warning message 
C           for exit velocities larger than 250 m/s. Comparable changes 
C           were incorporated also in subroutine PPARM for inputs on the 
C           SRCPARAM keyword.
C
C       2.  Modified subroutine MEOPEN to adjust the acceptable AERMET
C           version date from 06341 to 11059. AERMOD will no longer run
C           with met data based on version 06341 of AERMET.  AERMOD will 
C           run using met data based on version 11059. However, a warning 
C           message will be generated, and users are strongly encouraged 
C           to update their meteorological data to version 12345 of AERMET
C           due to the scope of changes included in that update.
C
C           
C-----  MODIFIED BY:    Roger W. Brode
C                       U.S. EPA, OAQPS/AQAD
C                       Air Quality Modeling Group
C
C                       December 10, 2012
C
C-----  MODIFIED FROM:        AERMOD
C                       (Version Dated 12060)
C
C
C=======================================================================
C
C       This revised version of AERMOD (dated 12060) includes the
C       following modifications relative to the previous version
C       (dated 11353); see MCB#7.
C
C-----  Bug Fixes:
C
C       1.  Modified subroutine MAXDCONT_LOOP to correct problems with the
C           MAXDCONT option for applications that vary emissions (EMISFACT), 
C           background ozone data (O3VALUES), or background concentrations 
C           (BACKGRND) with a day-of-week component, e.g., SHRDOW or SHRDOW7, 
C           etc. Also modified subroutine MAXDCALC to correct the MAXDCONT 
C           option for applications that include hourly emissions (HOUREMIS) 
C           for at least one source, but not all sources. Both of these bugs 
C           associated with the MAXDCONT option may have caused significant 
C           errors in the source group contributions in the MAXDCONT output 
C           file. 
C   
C           Subroutine MAXDCONT_LOOP was also modified to include checks on 
C           the consistency between results in the SUMVAL_MAXD and SUMHNH 
C           arrays for the "target" source group under the MAXDCONT option. 
C           The SUMVAL_MAXD array values are recalulated during the MAXDCONT 
C           processing, whereas the SUMHNH array values are calculated during 
C           normal processing. Results in these two arrays should match for 
C           the "target" source group. These consistency checks will serve
C           to flag any potential problems with the MAXDCONT option.
C
C       2.  Modified subroutine O3READ to correct a problem with the OLM and
C           PVMRM options under the EVENT processing mode when only the 
C           CO OZONEVAL keyword is used to specify a background ozone value,
C           without an hourly ozone file through the CO OZONEFIL keyword or
C           varying ozone values through the CO O3VALUES keyword.  Previous
C           versions effectively assumed a background ozone value of 0 (zero)
C           in these cases, resulting in underestimated source contributions 
C           in the EVENT output file.
C
C       3.  Modified subroutine LTOPG to use code from CTDMPLUS (also used
C           in CALPUFF) for estimating PG class from Monin-Obukhov length 
C           and surface roughness based on Golder (1972). The CTDMPLUS code
C           more closely matches the PG-class curves in Figure 4 of the Golder 
C           (1972) paper. The LTOPG routine is only used to determine the
C           stability-dependent distance for transitioning to a virtual point 
C           source approximiation for area sources under the FASTAREA
C           or FASTALL options.
C           
C-----  Enhancements:
C   
C       1.  Modified subroutines PRESET, SRCSIZ, and URBANS to allow for 
C           the use of URBANSRC ALL on the SO pathway to indicate that 
C           all sources are to be treated as URBAN sources. This option 
C           assumes that only one (1) urban area has been defined using 
C           the CO URBANOPT keyword.
C   
C       2.  Modified subroutines PRESET and MECARD and added subroutine 
C           NUMYR to allow the user to specify the number of years of 
C           meteorological data that are being processed for a particular 
C           run.  The option is exercised with the new NUMYEARS keyword on 
C           the ME pathway. The value specified on the NUMYEARS keyword 
C           is used to allocate storage for the arrays that are used in 
C           the MAXDCONT option, and allows the user to reduce the memory 
C           storage requirements under the MAXDCONT option when less than 
C           five (5) years of met data are being used. The default number 
C           of years used for the MAXDCONT array allocation without the 
C           NUMYEARS keyword is still 5 years (formerly specified by 
C           the NYEARS PARAMETER).
C   
C-----  Miscellaneous:
C   
C       1.  Modified subroutine ALLRESULT to eliminate the arrays of 
C           profile met data (the observed data in the PROFFILE input) 
C           by hour-of-year, level, and year for use with the MAXDCONT 
C           option. The profile met data arrays are not needed for the 
C           MAXDCONT option and their removal reduces the memory 
C           requirements for that option. Additional adjustments to the
C           array allocations for the MAXDCONT option were made based
C           on the options used for a specific model run in order to 
C           minimize the memory requirements for the MAXDCONT option.
C           Subroutine ALLRESULT was also modified to improve the 
C           accuracy of the memory storage estimates included on the
C           first page of the AERMOD.OUT file.
C   
C       2.  Modified subroutine OUMAXD_CONT to include checks on the use 
C           of a limited range of ranks (specified on the OU RECTABLE 
C           keyword) with the THRESH option on the OU MAXDCONT keyword. 
C           A fatal error message will be generated if the range of ranks 
C           specified is less than or equal to the design value rank for 
C           the specified pollutant plus 4, i.e., a fatal error will be 
C           generated if the range of ranks is less than or equal to 8 
C           for 1-hr SO2, or less than or equal to 12 for 1-hr NO2 or 
C           24-hr PM2.5. A non-fatal warning message is also generated 
C           if the range of ranks is less than or equal to the design 
C           value rank plus 20, i.e., if the range of ranks is less than 
C           or equal to 24 for 1-hr SO2, or less than or equal to 28 for 
C           1-hr NO2 or 24-hr PM2.5.
C
C
C-----  MODIFIED BY:    Roger W. Brode
C                       U.S. EPA, OAQPS/AQAD
C                       Air Quality Modeling Group
C
C                       February 29, 2012
C
C-----  MODIFIED FROM:          AERMOD
C                       (Version Dated 11353)
C
C
C=======================================================================
C
C       This revised version of AERMOD (dated 11353) includes the
C       following modifications relative to the previous version
C       (dated 11103); see MCB#6.
C
C-----  Bug Fixes:
C
C       1.  Modified subroutines METEXT, MAXD_METEXT, ALLRESULT to fix the
C           MAXDCONT option for URBAN applications. Source contributions 
C           reported in MAXDCONT files generated by the previous version
C           could be erroneous for applications involving URBAN sources.
C
C       2.  Modified subroutine OUTQA to include logical variable PLFILE 
C           for PLOTFILEs in the check for output file options that use
C           the FILEFORM keyword specifying whether fixed-format or 
C           exp-format is used. Previous versions erroneously issued 
C           warning message W399 if PLOTFILE was the only relevant output 
C           option used with the FILEFORM keyword.
C
C       3.  Modified subroutine O3VALS to correct the test for number of 
C           parameters to be .GE. 4 to allow for the ANNUAL option.
C
C       4.  Modified subroutine METEXT to only increment the hour index and 
C           year index if MAXDCONT option is being used. This change avoids 
C           unnecessary runtime error if the met data file includes a single 
C           hour (or more) beyond the maximum number of years specified by the 
C           NYEARS PARAMETER, which is initialized to 5 in the MAIN1 module.
C
C       5.  Modified subroutine PRTO3VALS to use the O3FLAG variable which
C           identifies the user-specified option for defining temporally-
C           varying background ozone concentrations. The previous version 
C           erroneously referenced the BFLAG variable associated with 
C           user-specified BACKGROUND concentrations. This error would not
C           have affected modeled concentrations based on the previous 
C           version of the model, but could result in a runtime error due 
C           to an array subscript out-of-bounds for some applications, or 
C           could result in skipping the summary of user-specified ozone
C           values in the output file for some applications.
C
C
C-----  MODIFIED BY:    Roger W. Brode
C                       U.S. EPA, OAQPS/AQAD
C                       Air Quality Modeling Group
C
C                       December 19, 2011
C
C-----  MODIFIED FROM:          AERMOD
C                       (Version Dated 11103)
C
C
C=======================================================================
C
C       This revised version of AERMOD (dated 11103) includes the
C       following modifications relative to the previous version
C       (dated 11059); see MCB#5 and updated User's Guide Addendum.
C
C-----  Bug Fixes:
C
C       1.  Modified subroutine MXDLYFL to include an IF-THEN block 
C           to account for cases with NHIMXDLY = 1, i.e., only the 
C           1st-highest rank was selected on the RECTABLE keyword.  
C           The previous version (11059) resulted in all short-term 
C           values being 0.0 if only the 1st-highest rank was selected 
C           for applications involving the special processing for 
C           daily maximum values (24hr PM25, 1hr NO2 and 1hr SO2).
C           This bug did not affect the results for applications
C           that included other ranks on the RECTABLE keyword,
C           including single ranks other than 1st-highest or any
C           range of ranks, or applications that do not utilize
C           the special processing for the 24hr PM25, 1hr NO2 and 
C           1hr SO2 NAAQS.
C
C-----  Miscellaneous:
C
C       1.  Modified subroutines O3FIL and BACK_GRND to include checks
C           for potential problems with the Fortran format specifier 
C           for hourly ozone files and hourly background files.  The 
C           hourly ozone and/or background concentrations may have been 
C           assigned values of zero (0) in previous versions of AERMOD 
C           if the user-specified Fortran format includes an integer 
C           (I) format to read the concentration values.  This could
C           significantly affect modeled concentrations for the OLM
C           and PVMRM options based on hourly ozone data, without any
C           clear indication of a problem.  The user-specified format
C           must use an integer format to read the date variables and
C           a real format (F, E, or D) to read the concentration 
C           variable.  Warning messages will be generated if the format
C           does not meet these requirements.  In addition, the compiler
C           options for the AERMOD executable distributed on SCRAM have 
C           been modified to include an option to check for consistency 
C           between variable types and format specifiers, which will 
C           cause AERMOD to issue a fatal error when reading the file 
C           in cases where real variables are read with an integer
C           format, or vice versa.  The AERMOD User's Guide Addendum 
C           has also been modified to clarify the requirements for 
C           user-specified Fortran formats with these options.
C
C
C-----  MODIFIED BY:    Roger W. Brode
C                       U.S. EPA, OAQPS/AQAD
C                       Air Quality Modeling Group
C
C                       April 13, 2011
C
C-----  MODIFIED FROM:          AERMOD
C                       (Version Dated 11059)
C
C
C=======================================================================
C
C       This revised version of AERMOD (dated 11059) includes the
C       following modifications relative to the previous version
C       (dated 09292); see MCB#4 and updated User's Guide Addendum:
C
C-----  Bug Fixes:
C
C       1.  Modified subroutines PVMRM_CALC, MAJOR_CONT, and 
C           MOLES_NOX to include calls to subroutine SETSRC for 
C           all source types, and modified subroutine SETSRC to 
C           initialize arrays for area source dimensions and 
C           initial sigmas to correct initialization problems 
C           with the PVMRM option.  Also modified subroutine 
C           RELDISP to more appropriately account for initial
C           sigmas for volume and area sources in the calculation 
C           of relative dispersion coefficients for determining 
C           plume volume.
C
C       2.  Modified subroutines PVMRM_CALC and MOLES_NOX to
C           include calls to subroutine EMFACT to apply emission
C           factors, if appropriate, in order to use the 
C           EMISFACT-adjusted emission rates in the calculations
C           of the moles of NOx. Previous versions used the 
C           emission rate specified on the SRCPARAM keyword
C           to calculate the models of NOx.
C
C       NOTE:  These bugs related to the PVMRM option could have 
C              significant impacts on modeled concentrations using 
C              the PVMRM option with the EMISFACT option and/or with
C              mixed source types.  The magnitude and bias of the
C              differences associated with these bugs will depend on 
C              the specifics of the application, but there may
C              be a tendency to overestimate NO2 concentrations 
C              for POINT sources due to the first item when the 
C              application also included VOLUME sources with large
C              initial sigma values.
C
C       3.  Modified subroutine HRLOOP to skip call to CHK_ENDYR
C           to determine whether the end of the year has been 
C           reached if the NOCHKD or WARNCHKD options are specified, 
C           since the end-of-year processing is only applicable to
C           sequential meteorological data.  Since the NOCHKD option
C           is invoked if the SCREEN option is specified, the call
C           to CHK_ENDYR will also be skipped for screening met data
C           with AERSCREEN.
C           
C       4.  Modified subroutine CHKDAT to correct problems with 
C           date sequence checks involving gaps of complete days 
C           or complete years.  A gap of complete calendar years 
C           is allowed even without the NOCHKD or WARNCHKD options
C           to account for missing years due to data completeness
C           issues for applications involving design values that
C           are based on multi-year averages (e.g., 1hr NO2/SO2
C           and 24hr PM2.5).
C        
C       5.  Modified subroutine TERRST to correct problems with 
C           determining the number of calm and/or missing hours 
C           only for the data period that is processed.   
C
C       6.  Modified subroutine EVLINI to include initializations
C           for SZMAX(:) and HSBLMX(:) arrays used for EVALFILEs.
C
C       7.  Corrected subroutines OUMXFL, OUPOST, OUPLOT, PERPST, 
C           PERPLT, OUTOXX, OUSEAS, OURANK, and OUEVAL to set lower 
C           limit for acceptable user-specified file units to 31, 
C           in order to avoid potential conflicts with hardcoded 
C           file units. Previous version incorrectly used '.LT. 30' 
C           rather than '.LE. 30' for checks on user-specified file 
C           units, allowing for potential file units conflicts with 
C           the debug file for particle deposition velocities.
C
C-----  Enhancements:
C
C       1.  A number of enhancements have been incorporated to more 
C           fully support the form of the new 1-hour NO2 and SO2 
C           NAAQS, as well as the 24-hour PM2.5 standard.  The 
C           form of these standards are similar in that they are 
C           based on a ranked percentile value averaged over the 
C           number of years processed.  To more fully support 
C           implementation of recent guidance on these NAAQS, 
C           the RECTABLE keyword has been modified to allow 
C           user-specified ranks of short-term averages (for all 
C           pollutants) up to the 999th highest value.  The 
C           previous version of AERMOD was limited to the 
C           10th-highest value and also restricted the rank 
C           for the 24-hour PM2.5 NAAQS to the 8th-highest value
C           (corresponding to the 98th percentile of daily values
C           during a year).  Note that the range of ranks specified
C           on the RECTABLE keyword (not the individual ranks) also
C           determines the range of ranks that may be considered 
C           with the new MAXDCONT option, described below.
C
C       2.  Added new MAXDAILY option on the OU pathway to output a 
C           summary of daily maximum 1-hour values for each day
C           processed.  These files provide an interim output product
C           that may be used to analyze new 1-hour NO2 and SO2 NAAQS
C           based on a specified percentile rank of daily maximum 
C           1-hour values.
C
C       3.  Added new MXDYBYYR option on the OU pathway to output a
C           summary of daily maximum 1-hour values by year and rank.  
C           These files provide an interim output product that may
C           be used to analyze the new 1-hour NO2 and SO2 NAAQS based 
C           on a specified percentile rank of daily maximum 1-hour 
C           values.
C
C       4.  Added new MAXDCONT option on the OU pathway to output a
C           summary of source group contributions to high ranked values
C           for a target group, averaged across the number of years 
C           processed and paired in time and space.  The new MAXDCONT 
C           option is applicable to daily maximum values for the 24-hour 
C           PM2.5 NAAQS and the new 1-hour NO2 and SO2 standards, and
C           can be used to determine whether a source or a group of 
C           sources contributes significantly to modeled violations of
C           the NAAQS, paired in time and space.
C
C       5.  For applications addressing the 24-hour PM2.5 standard 
C           or the 1-hour NO2 and SO2 standards, which are based on 
C           ranked values averaged across the number of years modeled, 
C           the PLOTFILE option has been enhanced to include values 
C           for each of the years processed based on the specified 
C           rank, in addition to the multi-year average.
C
C       6.  Added new BACKGRND option on the SO pathway to allow 
C           users to specify background concentrations, which can
C           be added to impacts from modeled emission sources to 
C           determine cumulative impacts.  Background concentrations
C           can be varied temporally using options similar to the
C           EMISFACT keyword for temporally-varying source emissions.
C           The new BACKGRND keyword also allows an option to use
C           a separate file of background concentrations on an hourly
C           basis. Applications with hourly background concentrations
C           can also include temporally-varying background values 
C           based on the EMISFACT options, such as SEASHR for season
C           by hour-of-day, which are used to substitute for missing
C           hourly values.
C
C       7.  For applications using the OLM or PVMRM options for NO2, 
C           a new option for specifying background ozone concentrations 
C           has been incorporated.  Similar to the new BACKGRND keyword, 
C           the new O3VALUES keyword on the CO pathway allows the user 
C           to specify temporally-varying background O3 concentrations 
C           using many of the same options available on the EMISFACT 
C           keyword for source emission factors.  The O3VALUES keyword 
C           can be used by itself or in conjunction with an hourly 
C           ozone file.  In the latter case, the O3VALUES are used to 
C           fill in for missing values in the hourly ozone file.
C           A separate OZONUNIT keyword is also available on the CO
C           pathway to specify units for the concentrations input 
C           through the O3VALUES keyword.
C
C       8.  Incorporated the equilibrium NO2/NOx ratio component of the
C           PVMRM option into the OLM option for estimating conversion 
C           from NOx emissions to ambient NO2 concentrations. The same
C           NO2EQUIL keyword on the CO pathway can be used to specify
C           the equilibrium ratio for either option, and a default 
C           ratio of 0.90 is assumed for both options if the NO2EQUIL
C           option is omitted.
C 
C       9.  Modified subroutine DEBOPT to allow user to specify
C           debug output only for PVMRM or deposition options
C           on the DEBUGOPT keyword, avoiding large output files 
C           under the MODEL debug option. Debug output for deposition 
C           options will still be generated if the MODEL debug option 
C           is selected. See AERMOD User's Guide Addendum for details 
C           on the DEBUGOPT keyword.  Also assigned file unit 9 to 
C           variable PVMDBG for the PVMRM debug file, and adjusted 
C           the PVMRM debug output to report total PercentNO2, 
C           including in-stack NO2/NOx contribution.
C
C      10.  A modification to the urban option has been implemented 
C           to address issues with the transition from the nighttime 
C           urban boundary layer to the daytime convective boundary 
C           layer. Under the new default urban option, the model will
C           continue to apply the urban boundary layer approach for
C           urban sources until the daytime convective boundary layer
C           exceeds the population-dependent urban boundary layer 
C           height. This enhancement is desribed in more detail in 
C           Appendix E of the updated AERMOD User's Guide Addendum.  
C           A non-DFAULT option has also been included to allow
C           users to revert to the original urban implementation.
C
C      11.  Increased the maximum length of source IDs from 8 to 12
C           characters, and increased the length of EVENT names from
C           8 to 10 characters, involving modifications to several 
C           subroutines.
C
C      12.  Included a new NOHEADER keyword on the OU pathway to allow
C           users to suppress the file header records for formatted 
C           output files.
C
C
C-----  Miscellaneous:
C
C       1.  Modified subroutines CHK_ENDYR, PRTPM25, PRTPM25SUM,
C           OUHIGH, and PRTOPT to allow for user-specified rank for
C           processing PM2.5 24-hour averages to accommodate current
C           recommendations for PM2.5 modeling. Also changed the
C           name of array used to store these values from SUMH8H 
C           to SUMHNH.
C
C       2.  The table of distances used in calculating the dominant 
C           plume volume for the PVMRM option (in subroutine PLUME_VOL)
C           was adjusted to use a more logical progression of distance 
C           intervals and to reduce the total number of intervals. This 
C           change may affect results slightly, but will also reduce 
C           model runtime. 
C
C       3.  Subroutine WAKFLG was modified to no longer ignore potential 
C           downwash effects for stack heights that equal or exceed the 
C           EPA formula height.  The determination of whether building
C           downwash effects apply is based on the criterion implemented
C           within the PRIME downwash algorithm.
C
C       4.  Modified subroutine URBPOP to adjust the limit for issuing 
C           a warning for urban population out-of-range from 10,000 to 
C           21,206, which corresponds to a population density of 
C           750/sq-km for an area within a 3km radius, consistent with 
C           the Appendix W criterion for urban/rural determination based 
C           on the population density.
C
C       5.  Several miscellaneous changes to address output formatting 
C           issues, replace DO loops with array assignments for array
C           initializations, and other minor code cleanup.
C
C       6.  Moved setup-related subroutines for EVENT processing option
C           from the 'evcalc.f' source file to the 'evset.f' source file.
C
C
C-----  MODIFIED BY:    Roger W. Brode
C                       U.S. EPA, OAQPS/AQAD
C                       Air Quality Modeling Group
C
C                       February 28, 2011
C
C-----  MODIFIED FROM:          AERMOD
C                       (Version Dated 09292)
C
C=======================================================================
C
C       This revised version of AERMOD (dated 09292) includes the
C       following modifications relative to the previous version
C       (dated 07026); see MCB#3:
C
C-----  Bug Fixes:
C
C       1.  Modified subroutine OLM_CALC to correct initialization 
C           problem with OLMGROUP keyword option.  The NO2VAL and NO_VAL 
C           arrays need to be reinitialized for each receptor when the 
C           OLMGROUP keyword is used.  This bug may result in significant
C           errors in concentration estimates for applications of the
C           OLM option with OLMGROUPs.  More details regarding this bug
C           are provided in the Addendum to MCB#3.
C
C       2.  Modified subroutine CHK_ENDYR to check for allocation
C           status of allocatable arrays during array initializations;
C           this has caused runtime errors on some compilers related to 
C           initializing unallocated ALLOCATABLE arrays when ANNUAL or 
C           PERIOD averages were calculated.         
C
C       3.  Modified ISSTAT() array indices in subroutine SOCARD to
C           eliminate potential conflicts among different options
C           referencing the same index.
C
C       4.  Modified subroutine MEOPEN to correct problems with
C           processing of the 'SCREEN' option in the AERMET version
C           date field of the header record of the surface file to 
C           support the use of screening meteorological data with 
C           AERSCREEN.
C
C       5.  Modified subroutine SETIDG to use 'I8' format for an internal
C           read of the source group ID from NUMID to avoid runtime errors 
C           with some compilers.
C
C       6.  Modified subroutines METEXT and MEREAD to include a more 
C           robust check for stable or missing hours in calculation
C           of solar irradiance (QSW) for use in dry deposition
C           calculations.
C
C       7.  Modified subroutine OUTQA to correct write statement for
C           message '540' to accommodate MONTH average option, labeled
C           as 720-hr option. Also modified subroutine OUTQA to include 
C           error checking for file name and file unit conflicts across 
C           output file options.
C
C       8.  Modified subroutines OUMXFL, OUPOST, OUPLOT, PERPST, PERPLT,
C           OUTOXX, OUSEAS, OURANK, and OUEVAL to set lower limit for
C           acceptable user-specified file units to 31, in order to avoid 
C           potential conflicts with hardcoded file units.  Also modified 
C           file units for deposition debug files to be within the range
C           established for hardcoded file units.
C
C       9.  Modified subroutine METH_2 to check for out-of-range inputs
C           for fine mass fraction (finemass) for the Method 2 particle
C           deposition option.
C
C      10.  Modified subroutine PLUME_VOL to include a call to subroutine
C           PENFCT to calculate plume penetration factor.  This corrects
C           potential errors with the PVMRM option for penetrated plumes.
C    
C      11.  Modified subroutine SIGZ to include a check on the sign of 
C           TGEFF before calculation of BVFRQ to avoid potential runtime 
C           errors.
C
C      12.  Modified subroutine SRCQA to check for source being defined
C           as both particulate and gaseous.  Also modified SRCQA to check
C           for invalid shape of AREAPOLY source, identified with
C           an area = 0.
C
C      13.  Modified subroutine VDP1 to output ZRDEP instead of undefined
C           variable ZREF.
C
C      14.  Modified subroutine METQA to include error checking for a 
C           value of wind speed reference height (UREFHT) of less than 
C           0.001m for non-missing and non-calm hours. Fatal error issued
C           to avoid runtime errors that may occur for values entered as 
C           zero.
C
C      15.  Modified several subroutines to address potential problems
C           related to options that assume full years of meteorological
C           data.  This includes modifications to subroutines STAEND and 
C           METEXT for setting the month, day, and hour for the "end of 
C           the year", based on the first hour of the meteorological data 
C           file or on the STARTEND keyword, to resolve potential problems 
C           for PM-2.5 and ANNUAL average applications in which the 
C           "starting hour" is not 01.  Modified subroutine METEXT to 
C           determine "end of the year" for the MULTYEAR option, and 
C           modified subroutine HRLOOP to call CHK_ENDYR for the MULTYEAR 
C           option to allow error checking for use of incomplete years or 
C           periods longer than a year with the MULTYEAR option.  Modified 
C           subroutine RSINIT to adjust the "starting hour" read from a 
C           re-start file (based on the INITFILE or MULTYEAR keyword) by 
C           adding 1 hour for the re-started model run.  Also improved the 
C           error handling in subroutine STAEND for user-specified start and
C           end dates, including a check for the user-specified STARTEND
C           period being less than 1 complete year when the ANNUAL average
C           option and/or the MULTYEAR option are specified, resulting
C           in a setup error.
C
C      16.  Modified code in subroutines RSINIT and METEXT to include
C           checks for potential conflicts between the first date of met
C           data file or the user-specified start date on STARTEND keyword
C           being later than the "start date" from a re-start file.  This 
C           data gap results in a fatal error message for INITFILE re-starts
C           and a warning message for MULTYEAR re-starts.  Also checks for
C           overlapping periods for MULTYEAR applications with start date
C           from STARTEND keyword being earlier than start date from the 
C           MULTYEAR re-start file.  This condition results in a fatal error
C           message.
C
C      17.  Modified subruotines METEXT and MEREAD to check for the
C           meteorological data file starting hour, and issue appropriate
C           warnings regarding short-term averages for the first day
C           of data if the starting hour is not 01, since they may 
C           represent partial data periods.  Modifications were also 
C           made to subroutine MEREAD to cycle through the meteorological 
C           data file until the hour matches the hour loop index for 
C           EVENT processing.
C
C      18.  Modified MAIN program to adjust the order of calls to 
C           subroutine RSDUMP to write results arrays to the SAVEFILE
C           for ANNUAL averages and for PERIOD averages with the MULTYEAR
C           option.  These adjustments ensure that the summary of overall
C           maximum PERIOD averages with the MULTYEAR option are based on
C           the maximum individual PERIOD values across the years processed, 
C           while the maximum ANNUAL averages reflect the multi-year 
C           average of the ANNUAL values across the years processed.  
C           The order of calculating averages by season and hour-of-day 
C           with the SEASONHR option was also adjusted to ensure that the 
C           SEASONHR averages are correct across the years of data processed.
C           Subroutine PRTSUM was modified to include a message in the
C           summary tables of maximum PERIOD and high ranked short-term
C           averages to indicate the number of years represented by the 
C           results when the MULTYEAR option is used.  Subroutine PRTOPT 
C           was also modified to clarify what PERIOD averages represent 
C           when the MULTYEAR option is used.
C
C      19.  Modified several subrountines, including HRLOOP, METEXT,
C           HRQREAD, O3READ, and O3EXT, to provide additional error
C           handling for premature "end-of-file" (EOF) condition for
C           input meteorological data, hourly emissions data (HOUREMIS
C           keyword), and ozone data files.  A fatal error is generated 
C           if the data files end before the "end date" specified by the 
C           user on the STARTEND keyword.
C
C      20.  Modified subroutines OUMXFL and OUPOST to correct problems
C           with MAXIFILE and POSTFILE outputs for re-started model 
C           runs using the SAVEFILE/INITFILE option.  The previous 
C           version of AERMOD used the 8-digit date from the MAXIFILE
C           and POSTFILE files to compare with the "full" 10-digit
C           date read from the INITFILE, and also read the 8-digit 
C           date from the wrong columns for POSTFILEs.  Also added
C           error handling for missing MAXIFILE and/or POSTFILE
C           files with INITFILE option, and for MAXIFILE and/or
C           POSTFILE files with data past the start date for the
C           MULTYEAR option.
C
C      21.  Modified subroutine INCLUD to include the EV pathway to
C           ensure proper handling of INCLUDED files for events on the
C           EV pathway.  Also incorporated additional error handling
C           for INCLUDED files.
C
C      22.  Modified assignment of ILAND_NDX based on flow vector (AFV)
C           in subroutine VDP to be consistent with assignment of the
C           flow vector sector for downwash, IFVSEC.
C
C      23.  Modified subroutine OCALC to correct the calculation of the
C           adjusted emission rate for the point source approximation 
C           for OPENPIT sources under the FASTAREA option to use the 
C           length and width of the "effective" area source rather than 
C           length and width of the original openpit source.  Also 
C           modified subroutine ARDIST to recalculate the coordinates
C           for the center of the "effective" area source used for the
C           OPENPIT source algorithm, rather than using the center
C           coordinates of the original OPENPIT source.
C
C      24.  Modified several subroutines to correct and clarify aspects
C           of the logic related to use of the deposition algorithms
C           for both gas and particle deposition, especially in relation
C           to the depletion options for applications including both
C           gaseous and particulate emissions.  Previous versions of
C           AERMOD would allow a particulate source and a gaseous source
C           in the same model run for calculating concentrations only,
C           but would indicate that dry and wet depletion were being 
C           applied, even if no deposition parameters were specified 
C           for the gaseous source.  The updated version will only 
C           allow both gaseous and particulate sources in the same model
C           run for concentrations if deposition parameters are also 
C           specified for the gaseous source, or if dry and wet depletion 
C           are turned off by specifying the NODRYDPLT and NOWETDPLT 
C           options on the MODELOPT keyword.  Specifying NODRYDPLT and
C           NOWETDPLT removes any requirement for deposition parameters 
C           for the gaseous source. See the AERMOD User's Guide Addendum 
C           for more details regarding use of deposition algorithms. 
C
C-----  Enhancements:
C
C       1.  Incorporated additional options under EMISFACT keyword to 
C           vary emissions by hour-of-day and day-of-week (HRDOW and 
C           HRDOW7).
C
C       2.  Improved efficiency of allocating array storage, including
C           more precise allocation of array limits for AREAPOLY and
C           AREACIRC sources, more precise allocation for number of 
C           particle size categories for particulate sources, and 
C           allocating the building downwash arrays only when needed.  
C           Also included additional informational messages to identify 
C           where allocation errors occur with a list of array limits 
C           and preliminary memory storage estimate to facilitate 
C           diagnosis of problems.  Modifications made to MAIN program, 
C           MODULE MAIN1, and subroutines ALLSETUP, ALLRESULT, PRESET, 
C           and SRCSIZ.
C
C       3.  Modified subroutines HRQREAD (renamed from HQREAD), HRQEXT,
C           and EV_HRQEXT to allow option of specifying hourly varying 
C           release height and initial dispersion coefficients for VOLUME 
C           and AREA sources (including AREAPOLY and AREACIRC) through 
C           the optional HOUREMIS file.  Also eliminated redundant code 
C           for processing hourly emission files for 'normal' vs. 'event' 
C           processing to ensure consistency and simplify code maintenance.
C
C       4.  Increased maximum length for filenames to 200 (controlled by
C           the ILEN_FLD parameter in modules.f), and the maximum input
C           string length to 512 (controlled by the ISTRG parameter in
C           modules.f).  Also modified subroutine DEFINE to allow double 
C           quotes (") as field delimiters in the 'aermod.inp' file to 
C           support filenames with embedded spaces.
C
C       5.  Modified subroutines SRCQA and OCALC to allow for use of the 
C           non-DFAULT METHOD_2 option for particulate emissions and for
C           non-particulate (gaseous) emissions for OPENPIT sources.
C
C       6.  Modified several subroutines to allow for the non-DFAULT 
C           option of FLAT terrain to be specified for individual sources,
C           allowing both FLAT and ELEV terrain treatments within the 
C           same model run (see Section 4.1 of the AERMOD Implementation
C           Guide regarding modeling of sources with terrain-following 
C           plumes in sloped terrain).  This non-DFAULT option is 
C           activated by specifying both 'ELEV' and 'FLAT' on the 
C           CO MODELOPT keyword.  The user identifies which sources to
C           model with the FLAT terrain option by specifying 'FLAT' (not
C           case-sensitive) on the SO LOCATION keyword in place of the 
C           field for source elevation.  Sources treated as FLAT terrain 
C           will be identified in the 'aermod.out' file in the summary 
C           tables of source input parameters.  The 'aermod.out' page 
C           header of model options, and other output file headers, 
C           will include the field 'FLAT and ELEV' to identify use 
C           of this feature.                    
C
C       7.  Incorporated a non-DFAULT option to optimize model runtime
C           for POINT and VOLUME sources through an alternative 
C           implementation of the horizontal meander algorithm,
C           which preserves the centerline concentration based on 
C           the DFAULT meander algorithm, but uses of an effective 
C           sigma-y for the lateral spread of the plume. This eliminates 
C           the upwind component of dispersion that occurs under the 
C           DFAULT meander algorithm for POINT and VOLUME sources, 
C           which should significantly reduce model runtime.  This 
C           option is selected by including 'FASTALL' on the 
C           CO MODELOPT keyword.  The FASTALL option also activates 
C           the optimized treatment for AREA sources (including AREAPOLY, 
C           AREACIRC, and OPENPIT sources) formerly associated with the 
C           TOXICS option, which is now obsolete (see Item #22 below 
C           under 'Miscellaneous' and the AERMOD User's Guide Addendum).
C
C       8.  Incorporated non-DFAULT option for user-specified dry 
C           deposition velocity for gaseous emissions, using the GASDEPVD 
C           keyword on the CO pathway.  The user-specified dry deposition 
C           velocity will be used for all gaseous sources.  Particulate 
C           sources can also be included in the same model run, but no 
C           wet deposition calculations (WDEP, DEPOS or WETDPLT) will
C           be allowed.
C
C       9.  Modified OPEN statements for input data files to use the
C           ACTION='READ' specifier to allow file sharing across 
C           multiple model runs.  This includes the 'aermod.inp' input
C           file, INCLUDED files referenced from 'aermod.inp', surface 
C           and profile meteorological inputs, hourly emission files, 
C           and hourly ozone data files.
C
C      10.  Added new SUMMFILE option on the OU pathway to output the 
C           summary of high ranked values to a separate file.  The new
C           SUMMFILE includes the "MODEL SETUP OPTIONS SUMMARY" page, 
C           the summary of meteorological data periods processed and
C           summary of first 24 hours of meteorological data, the
C           summaries of high ranked values provided at the end of 
C           the standard output file, and the summary of error and 
C           warning messages.
C
C      11.  Added new FILEFORM option on the OU pathway to allow
C           the user to specify exponential-formatted rather than
C           fixed-format output of model results (CONC, DEPOS, DDEP, 
C           and/or WDEP) for external results files.  Output file 
C           options affected by this feature are MAXIFILE, PLOTFILE,
C           POSTFILE (using PLOT format), RANKFILE, and the SEASONHR
C           file.  The default will continue to be fixed format for
C           these output files.  Results included in the 'aermod.out' 
C           file and optional SUMMFILE output file are not affected
C           by this new keyword.  Also included checks on values that
C           may exceed the output format limit, and a warning message
C           is generated if the FILEFORM = EXP option is not used.
C
C      12.  Modified several subroutines to include model run date
C           and run time in the header records for output files,
C           including the MAXIFILE, PLOTFILE, POSTFILE (using PLOT 
C           format), RANKFILE, and SEASONHR file options.
C
C      13.  Modified several subroutines to include a new option 
C           on the MODELOPT keyword to issue warnings rather than
C           fatal errors for records out of date sequence in the
C           meteorological data files.  The new option is selected
C           with the WARNCHKD parameter on the MODELOPT keyword.  
C           This option is primarily intended for multi-year 
C           meteorological data files that include gaps between 
C           years of meteorological data.  The WARNCHKD option is 
C           allowed under the regulatory DFAULT option, but should 
C           only be used when processing files with data gaps is 
C           clearly documented and justified.
C
C-----  Miscellaneous:
C
C       1.  Modified subroutine URBOPT to prohibit use of urban roughness
C           length not equal to 1.0m for regulatory DFAULT applications, 
C           as discussed in Section 5.3 of the AERMOD Implementation Guide.  
C           Also modified limits on urban roughness length used to generate 
C           warning messages.
C
C       2.  Modified code (most modules and subroutines) to explicitly
C           declare most non-integer variables as DOUBLE PRECISION.
C           This addresses some long-standing concerns regarding the 
C           potential sensitivity of the model to precision involving 
C           UTM coordinates, with the full Northing coordinate near 
C           the limit of single precision.  More serious concerns have 
C           arisen recently with a case showing a consistent negative 
C           bias of about 3 percent for modeled concentrations based on 
C           source group ALL as compared to the sum of concentrations 
C           from the same sources grouped separately.  These differences 
C           were attributable to precision errors, with the group ALL 
C           results biased low due to the impacts from many sources being 
C           truncated as a result of the wide range of impacts across  
C           sources.  This change should also result generally in more 
C           consistent concentration estimates from AERMOD across different
C           compilers, compiler options, and computing platforms.
C
C       3.  Modified code in subroutines METEXT and MEREAD for processing
C           multi-year meteorological data files to determine whether an 
C           embedded header record is included between individual calendar
C           years of data in the surface file, which may occur if the files
C           have been concatenated.  An error while attempting to read 
C           the data as a header record will be interpreted to mean that 
C           there is no embedded header record, while the presence of a 
C           colon (':') will be interpreted to mean that there is an embedded 
C           header record. AERMOD should accept data with or without the
C           embedded header records.  However, if multi-year surface files
C           based on non-calendar year periods are concatenated, then the 
C           embedded header records must be removed before input to AERMOD.
C           The revised code also issues warning messages if UAIR and 
C           SURF IDs do not match inputs in the runstream file for multiple 
C           years since AERMOD allows mismatch (with warning) for single 
C           year files.
C
C       4.  Modified subroutines PRTANN, SPRTHT, and PRTPM25 to include
C           EVALCART receptors with DISCCART receptors for output tables.
C
C       5.  Expanded array sizes (in 'modules.f') to 50 for ICSTAT,
C           ISSTAT, IRSTAT, IMSTAT, and IOSTAT arrays.  Also changed the
C           index used for the FINISHED keyword on each pathway to 50
C           instead of 25.
C
C       6.  Removed obsolete code inherited from ISCST3 code, including
C           TGSET.FOR and other references to the 'TG' (terrain grid)
C           pathway; reference to IRSTAT(7) in subroutine RECARD and
C           other references to BOUNDARY receptors in subroutine PRTDAY;
C           'HE>ZI' option; and NWET parameter.  Also removed obsolete
C           MODULE DEPVAR from the 'modules.f' source file.
C
C       7.  Replaced *.pri "INCLUDE" files used for global data storage
C           in PRIME subroutines with MODULE subprograms. The new MODULE
C           subprograms are contained in the 'modules.f' source file.
C           Also imposed explicit variable type declarations through use
C           of IMPLICIT NONE within the PRIME portions of the AERMOD code.
C
C       8.  Incorporated several modifications to subroutine PRTOPT, 
C           including a more "refined" estimate of memory storage
C           requirements, including the DFAULT urban roughness 
C           length and a more complete summary of options when DFAULT
C           option is not specified, correcting cosmetic problems with
C           the output option summary for EVENT vs. standard processing,
C           clarification of deposition option logic, inclusion of 
C           emission and output units for both concentrations and
C           deposition when needed, and more complete explanations 
C           for some options.
C
C       9.  Modified subroutine SRCQA to issue warning messages for 
C           sources that are not included in any SRCGROUPs, modified
C           checks for urban areas without urban sources, and changed
C           warning to fatal error for urban areas with no urban sources.
C           Additional error-checking for OLMGROUPs and PSDGROUPSs has
C           also been included.
C
C      10.  Modified criterion for issuing a warning message regarding
C           aspect ratio for rectangular AREA sources being out of range,
C           from an aspect ratio of 10:1 to an aspect ratio of 100:1. 
C           The upper limit of aspect ratio for stable performance of 
C           the numerical integration algorithm for area sources has not 
C           been fully tested and documented, and may vary depending on 
C           the specifics of the application.  A ratio of 10:1 is probably 
C           too strict and may unnecessarily lead to a large number of 
C           warning messages in some cases.  Users should always carefully 
C           quality assure the source inputs provided to AERMOD for 
C           accuracy and appropriateness for the application.
C
C      11.  Added error message (number 499) to handle PRIME plume rise
C           error when maximum number of downwind distances for numerical
C           plume rise (MXNW) is exceeded. This error condition has
C           occurred with application of Model Clearinghouse procedure
C           for capped stacks (using large effective stack diameter),
C           which is not appropriate with PRIME algorithms. Also
C           modified subroutine NUMRISE to include additional information
C           regarding this error message. 
C
C      12.  Modified subroutine PRTSRC to remove reference to "STABILITY
C           CATEGORY" and correct format statement 9024 for 'WSPEED' 
C           EMISFACT option (inherited from legacy ISCST3 code for
C           'STAR' option).
C
C      13.  Modified subroutines PPARM, VPARM, APARM, APPARM, ACPARM, and 
C           OPARM to issue a fatal error message for source release heights 
C           that exceed 3,000 meters.  Such abnormally high release heights 
C           have resulted in runtime errors and are considered beyond the 
C           range of applicability for AERMOD.
C
C      14.  Modified subroutines RANKFL and VARINI to increase number of
C           high ranked values allowed in the RANKFILE option without
C           overflowing the field, from 999 to 999,999.
C
C      15.  Modified subroutine HEADER in aermod.f to use ASCII
C           form feed character (ACHAR(12)) in 'aermod.out' file to 
C           eliminate need for the non-standard CARRIAGECONTROL='FORTRAN'
C           option in the OPEN statement for the 'aermod.out' file.  
C           Modified subroutine FILOPN to eliminate the use of the 
C           CARRIAGECONTROL='FORTRAN' option when opening the
C           'aermod.out' file.  Also modified subroutine SETUP 
C           to remove '1X' from format statements for echoing 
C           the runstream inputs to the output file since Fortran 
C           carriage-control is no longer applied.  These modifications 
C           will improve the portability of the AERMOD code to other
C           compilers and platforms.
C
C      16.  Modified several subroutines to reflect updates to the
C           listing of MODELOPT options based on other modifications
C           to the model.  The summary of model options in the page
C           header of output files has been modified to more clearly
C           indicate when the regulatory default mode (DFAULT) is in
C           effect and whether non-DFAULT options are being used for
C           applications where the DFAULT option has not been specified.
C           Also included error checking for conflicting options being
C           included on the MODELOPT keyword, such as the DRYDPLT and
C           NODRYDPLT options.
C
C      17.  Modified subroutine PNPOLY to use an Internal Function 
C           for EOR rather than a Statement Function, which has been 
C           identified as obsolescent in Fortran 95.  
C
C      18.  Modified subroutine METEXT to check for the presence of 
C           additional variables in the surface meteorological file 
C           needed for use of the deposition options when required.
C           A fatal error is issued if less than the minimum number 
C           of variables is included in the surface file.
C
C      19.  Modified subroutines SURFIL and PROFIL in meset.f to remove 
C           the optional user-specified READ format for the surface and 
C           profile meteorological input files.  AERMOD assumes FREE 
C           format for all input meteorological data files, i.e., the data 
C           fields must be space or comma-delimited.  Note that the 
C           "header" record for the surface file must conform to the 
C           format used by AERMET.
C
C      20.  Modified subroutines SOLOCA and RECARD to check for missing
C           source or receptor elevations, coded as -9999.0 by AERMAP.
C           Missing elevations will generate a fatal error, requiring the
C           user to resolve the issue(s) associated with the missing 
C           elevations.
C
C      21.  Modified subroutine FILOPN in aermod.f to use lower case for
C           'aermod.inp' and 'aermod.out' files.  Also modified Fortran
C           source code filenames to use lower case with '.f' file
C           extensions.  These changes may be important for operating
C           systems that recognize case-sensitive file names.
C
C      22.  Modified several subroutines to remove the non-DFAULT 
C           TOXICS option, inherited by AERMOD from the ISCST3 model 
C           code.  The area source optimizations previously activated 
C           with the TOXICS option are now associated with the new
C           non-DFAULT 'FASTAREA' option, and are also included under
C           the new non-DFAULT 'FASTALL' option (see Item #7 above
C           under 'Enhancements' and the AERMOD User's Guide Addendum).
C           The TOXICS option is no longer required to utilize the 
C           gas deposition options or the METHOD_2 option for particle
C           deposition.  However, these options are still considered
C           to be non-DFAULT options within AERMOD.  Also modified
C           several subroutines to clarify reporting of deposition
C           options in the 'aermod.out' file, including previously 
C           undocumented options to turn off dry depletion (NODRYDPLT)
C           and wet depletion (NOWETDPLT).
C
C      23.  Modified subroutine MODOPT in coset.f to use fatal error
C           (message number 204) for several options that conflict 
C           with the regulatory DFAULT option, including the BETA,
C           OLM, PVMRM, PSDCREDIT, SCIM, FASTALL, FASTAREA, and 
C           AREADPLT options.
C
C      24.  Modified subroutines COCARD and MECARD to skip keywords
C           that are not valid for EVENT processing runs, including
C           EVENTFIL, SAVEFILE, INITFILE, and MULTYEAR on the CO 
C           pathway, and STARTEND and DAYRANGE on the ME pathway.  
C           These keywords are not echoed to the input file created 
C           with the EVENTFIL keyword, so this change would only 
C           impact cases where these keywords were added to the
C           EVENTFIL after the file was created by the standard 
C           run, or EVENT input files created through other means.
C
C      25.  Modified subroutines HRQREAD and O3EXT to allow for 
C           4-digit years in the HOUREMIS and OZONEFIL data files.
C
C      26.  Modified subroutine CHKREC to include a check for 
C           receptors beyond MAXDIST from sources, using the center 
C           coordinates for AREA/AREACIRC/AREAPOLY and OPENPIT 
C           sources.  MAXDIST is set to 80km under the non-DFAULT 
C           FASTALL and FASTAREA options.
C
C      27.  Modified subroutine MYEAR to treat the 'H6H' field 
C           for the MULTYEAR keyword as optional, with a warning 
C           indicating that it is no longer required.
C
C      28.  Modified subroutine METQA to calculate the total 
C           precipitation amount from the surface met data file.
C           Also modified subroutine SUMTBL to include the total
C           precipitation amount with the message summary in the
C           'aermod.out' file, and included a warning message for
C           applications using wet deposition algorithms with a 
C           total precipitation amount of zero (0).
C
C      29.  Modified subroutine TERRST to determine the number of
C           calm and/or missing hours only for the meteorological
C           data period being processed when the STARTEND and/or 
C           DAYRANGE keywords are used.  Also modified subroutine 
C           SUMTBL to include the number of hours processed from 
C           the meteorological data file in the message summary.
C
C      30.  Modified subroutines PRTSUM and PRTPM25SUM to adjust 
C           format of column headers and other write statements in
C           the summary tables of high ranked values, including removal 
C           of '1X' used to skip the Fortran carriage-control character, 
C           which is no longer needed.
C
C      31.  Miscellaneous code clean-up, including the removal of 
C           extraneous Fortran 'SAVE' statements and specifying the
C           SAVE attribute only where needed for local variable type
C           declaration statements, removal of subroutines PRESOINC
C           and PREREINC to eliminate unnecessary code redundancy, 
C           removal of unused subroutine EV_CHKDAT, removal of obsolete
C           PRIME subroutine NUMMET, removal of unused PARAMETERS in 
C           MODULE MAIN1, replacement of DO loops for variable 
C           initialization with array assignments, improved consistency 
C           of error handling for processing of numeric inputs from 
C           the runstream file, adjusting header format in subroutine
C           SHOUT for SEASONHR files with multiple output typse, and 
C           removal of code associated with "wet SCIM'ing" option 
C           which is not supported in AERMOD.  The input format for
C           the SCIMBYHR keyword has also been adjusted to remove 
C           the parameters associated with wet SCIM'ing, but the
C           model should still process input files based on the 
C           previous format while issuing a warning message.
C
C
C-----  MODIFIED BY:    Roger W. Brode
C                       U.S. EPA, OAQPS/AQAD
C                       Air Quality Modeling Group
C
C                       October 19, 2009
C
C-----  MODIFIED FROM:          AERMOD
C                       (Version Dated 07026)
C
C=======================================================================
C
C       This revised version of AERMOD (dated 07026) includes the
C       following modifications relative to the previous version
C       (dated 06341); see MCB#2:
C
C-----  Bug Fixes:
C
C       1.  Corrected array indexing problem for POINT, POINTCAP, and
C           POINTHOR sources in subroutine SOPARM that could result 
C           in incorrect processing of SRCPARAM data for some sources 
C           and/or incorrect issuance of fatal runtime errors for number
C           of parameters specified (error code 'E201').
C
C       2.  Modified subroutine SRCQA to include QA checks on the
C           number of emission factors input for MHRDOW and MHRDOW7
C           options.
C
C       3.  Modified subroutines SRCSIZ and PRESOINC to correct
C           potential problem with assignment of array size for
C           MHRDOW option.
C
C       4.  Modified subroutine SRCQA to use DOUBLE PRECISION in the
C           calculation of area and centroid coordinates for AREAPOLY
C           sources. This change avoids problems encountered with
C           the Compaq Visual Fortran compiler producing erroneous
C           results for some compiler options.
C
C       5.  Modified subroutine MEREAD to assign non-array logical
C           variables STABLE and UNSTAB for use in subroutine COMPTG
C           for EVENT processing mode.
C
C       6.  Modified subroutines PCALC, VCALC, ACALC and OCALC to
C           include calls to subroutine HEFF prior to calculation of
C           zsubp for deposition applications.
C
C       7.  Modified subroutine HRQEXT for the HOUREMIS option to
C           correct processing of missing parameters for point sources
C           to assign all parameters to 0.0 if any of the parameters 
C           are missing, in conformance with Section 3.3.9 of the 
C           AERMOD User's Guide.
C
C       8.  Modified subroutine VDP to prevent potential zero-divide
C           condition for cases with zero specific humidity deficit.
C
C       9.  Modified subroutines VARINI and RESINI to check for 
C           allocation status of allocatable arrays during array 
C           initializations.
C
C-----  Miscellaneous:
C
C       1.  Included additional error-checking for the non-DFAULT,
C           BETA-test PSDCREDIT option.
C
C       2.  Modified subroutine PRTSRC to correct format for printing
C           QFLAG from A6 to A7 to accommodate SHRDOW7 and MHRDOW7,
C           and other minor adjustments to formatting.
C
C       3.  Modified code to read the header record of the surface file
C           in subroutines MEOPEN, METEXT, and MEREAD to include a 
C           separate test on the AERMET version date field under the 
C           SCREEN option, to allow for the future use of screening 
C           meteorology that is not directly linked to a specific
C           version of the AERMET processor.  Also modified subroutine
C           METDAT to include meteorological data version date in the
C           summary of met data inputs.
C
C
C-----  MODIFIED BY:    Roger W. Brode
C                       U.S. EPA, OAQPS/AQAD
C                       Air Quality Modeling Group
C
C                       January 24, 2007
C
C       MODIFIED FROM:          AERMOD
C                       (Version Dated 06341)
C
C=======================================================================
C
C       This revised version of AERMOD (dated 06341) includes the
C       following modifications relative to the previous version
C       (dated 04300); see MCB#1:
C
C-----  Bug Fixes:
C
C       1.  Initialize HE = HS prior to first call to ADISY in
C           subroutines ACALC and OCALC.
C
C       2.  Replaced undefined variable lbd with .TRUE. in calls to
C           WAKE_XSIG to ignore BID for "outside" cavity source
C           in subroutine CAV_SRC.
C
C       3.  Modified subroutine PRMDELH to save L_INWAKE to local
C           variable.
C
C       4.  Modified subroutines WAKE_DFSN and WAKE_DFSN2 to initialize
C           dummy variables zkdum and ykdum to 0.0.  Otherwise these
C           variables may be undefined in subrountine WAKE_SIG under
C           some circumstances.
C
C       5.  Moved call to RSDUMP for MULTYR option from subroutine
C           HIPER to MAIN to resolve minor discrepancy in ANNUAL
C           averages for MULTYR option.
C
C       6.  Added arrays to store WDSIN and WDCOS by source for use
C           in determining the wind direction for the dominant source
C           in the PVMRM option.
C
C       7.  Deleted code in subroutine PCALC that adjusts DHCRIT based
C           on distance to well-mixed layer for stable conditions for
C           consistency with Model Formulation Document.
C
C       8.  Added error checking for missing GASDEPOS inputs.
C
C       9.  Corrected variable type from INTEGER to REAL for variable
C           NEWRAD (effective radius) in subroutine GENCIR for AREACIRC
C           sources.
C
C-----  Enhancements:
C
C       1.  A new 'BETA' option has been added to the CO MODELOPT card
C           to identify and allow use of new features added to the model
C           that are still in a draft BETA-test status.  The  BETA option
C           is a non-DFAULT option, and will be overridden if the DFAULT
C           option is specified.
C
C           The following two draft enhancements are included in
C           this update to AERMOD under the BETA option:
C
C           a) options for capped stacks (source type = POINTCAP) and for
C              horizontal releases (source type = POINTHOR); and
C           b) the PSDCREDIT option for PVMRM to account for NO/NO2 plume
C              chemistry of combined plumes in the computation PSD credits.
C
C           Inclusion of these draft BETA-test options does not imply
C           any endorsement of their use for regulatory or non-regulatory
C           applications of the model. In addition, the designation of
C           BETA-test to these draft enhancements does not imply that
C           these options have completed rigorous internal (Alpha) testing
C           prior to being included in a public release of the model.
C
C       2.  Additional options under EMISFACT keyword to vary emissions
C           by month, hour-of-day, and day-of-week (MHRDOW and MHRDOW7).
C
C       3.  Expanded urban option to allow multiple urban areas to
C           be defined through multiple URBANOPT cards.  A new
C           Urban ID parameter has been added to link sources with
C           specific urban areas.  Existing input files with a single
C           urban area do NOT need to be modified to run with this
C           version of the model.
C
C       4.  The following changes have been made to ensure consistency
C           with the current PM NAAQS:
C
C           a) Added special processing for PM-2.5 to calculate design
C              values in accordance with the PM NAAQS.  The design value
C              for 24-hour averages is based on the high-eighth-high (H8H)
C              averaged over N years, as an unbiased surrogate for the 98th
C              percentile.  The long-term design value for PM-2.5 is based
C              on the highest annual average concentration averaged over
C              N years using the ANNUAL keyword on the AVERTIME card.
C
C           b) The "post-1997" PM-10 processing based on H4H averaged over
C              N years has been removed since that standard was vacated.
C              The PM-10 design value for 24-hour averages is based on the
C              high-sixth-high (H6H) over five years (for NWS data), or more
C              generally by the high-N+1-high value over N years. This can
C              be accomplished using the existing CO MULTYEAR option with
C              multiple 1-year input met data files, or using single
C              five-year data files without the MULTYEAR option.
C
C       5.  Added option to specify initial in-stack NO2 ratio for
C           PVMRM and OLM options with CO NO2STACK card.  Default ratio
C           without the NO2STACK card for OLM is 0.10.  No default ratio
C           has been determined for PVMRM, so user must either use the
C           CO NO2STACK card to initialize ratio for all sources, or
C           specify ratio for each source using the SO NO2RATIO card.
C           The SO NO2RATIO card can be used to override the value
C           specified on the CO NO2STACK card for specific sources.
C
C       6.  The maximum number of vertices for an AREAPOLY source is now
C           allocated dynamically at runtime.  The previous version of
C           AERMOD allowed a maximum of 20 vertices.  The maximum number
C           of vertices allowed for AREACIRC sources is not allocated
C           dynamically.  However, the maximum number of vertices allowed
C           for AREACIRC sources will always be at least 50, and if both
C           AREAPOLY and AREACIRC sources are included in the same model
C           run, then the maximum number of vertices identified for the
C           AREAPOLY sources will also define the maximum number allowed
C           for the AREACIRC sources if it is larger than 50.
C
C-----  Miscellaneous:
C
C       1.  Added range check on value of vertical potential temperature
C           gradient above ZI (variable VPTGZI), to avoid problems with
C           data provided from sources other than AERMET.  A minimum
C           value of 0.005 K/m is applied for consistency with AERMET,
C           and a warning message is also generated for values larger
C           than 0.10 K/m.
C
C       2.  Refined the process of dynamically allocating array storage
C           to skip allocation of arrays that are not needed based on the
C           model options selected, in order to reduce unnecessary memory
C           usage.
C
C       3.  Removed calls to subroutines VDP and SCAVRAT from subroutine
C           PLUME_VOL.
C
C       4.  Removed references to obsolete BOUNDARY keyword for receptors
C           in RECSIZ and PREREINC subroutines to allocate receptor arrays,
C           and in subroutines PRTANN, SPRTHT, and PRTPM25 of OUTPUT.FOR.
C
C       5.  Changed AMAX1/AMIN1 intrinsics to MAX/MIN.
C
C
C-----  MODIFIED BY:    Roger W. Brode
C                       U.S. EPA, OAQPS/AQAD
C                       Air Quality Modeling Group
C
C                       James O. Paumier (PSDCREDIT option)
C                       MACTEC Federal Programs, Inc.
C
C                       December 7, 2006
C
C       MODIFIED FROM:          AERMOD
C                       (Version Dated 04300)
C
C=======================================================================
C
C       This DRAFT version (dated 04300) includes the Plume Volume Molar
C       Ratio Method (PVMRM) and the Ozone Limiting Method (OLM) for
C       modeling conversion of NOx to NO2.  This work was supported by
C       BP Exploration (Alaska), Inc., Phillips Exploration, Inc.,
C       and the Alaska Department of Environmental Conservation.
C
C       This DRAFT version (dated 04300) also includes the following
C       modifications:
C
C       1.  Dry depletion (DRYDPLT) and wet depletion (WETDPLT) are no
C           longer optional for deposition applications.  These options
C           for removal of mass from the plume due to dry and/or wet
C           deposition processes will automatically be invoked for
C           applications in which dry and/or wet deposition are
C           considered.  The DRYDPLT and WETDPLT options on the
C           MODELOPT card will be ignored, and need not be removed
C           from the model input file for the model to run.
C
C       2.  Correction made to area source algorithm, subroutine PLUMEF,
C           to include a call to CRITDS to calculate the critical
C           dividing streamline height for gaseous pollutants.  Also
C           modified PLUMEF to correct a problem with the AREADPLT option.
C
C       3.  Corrections made to area source and openpit algorithms,
C           in subroutines ACALC and OCALC, to include tilted plume
C           for point source approximation of particle emissions, and
C           to include reinitialization of __VAL arrays at end of
C           receptor loop (reinitializations also included in PCALC and
C           VCALC for point and volume sources for consistency).  The
C           latter correction fixes a potential problem with particle
C           emissions for area sources when the point source
C           approximation is used under the TOXICS option.
C
C       4.  Corrected calling arguments for call to WAKE_SIG from
C           subroutine WAKE_DFSN2, to use wakiz and wakiy instead of
C           turbz and turby.
C
C       5.  Minor correction made to wet deposition calculations to
C           include lateral term (FSUBY) in weighting of direct
C           and penetrated source contributions for WETFLUX.
C
C       6.  Modified suroutine PRMCALC to place receptor on centerline
C           of cavity plumes by setting Y2 = 0.0 for SCREEN option.
C
C       7.  Modified subroutine SRCQA to calculate equivalent XINIT
C           and YINIT values for AREAPOLY sources to allow for
C           calculation of area of source under TOXICS option and
C           for PVMRM option.  Also modified SRCQA to include a more
C           refined computation of centroid for AREAPOLY sources.
C
C       8.  Included check in subroutine METQA for absolute values of
C           Monin-Obukhov length (OBULEN) less than 1.0.  Adjustment
C           of OBULEN is made to limit ABS(OBULEN) .GE. 1.0.  The
C           sign of OBULEN is assigned the opposite of the sign of the
C           heat flux if OBULEN is 0.0.  This limit on OBULEN is
C           already applied in AERMET, so this change in AERMOD will
C           only affect input data generated by other means.
C
C       9.  Moved call to SUB. METDAT ahead of call to SUB. SET_METDATA
C           to avoid potential problem with negative (missing) 
C           precipitation for first hour.
C
C      10.  Added range check on gas deposition parameters to trap
C           on input of zero (0.0) values.
C
C      11.  Modified subroutine METQA to reduce number of extraneous
C           warning messages, especially for hours with missing
C           meteorological data.  Also modified range check for missing
C           wind direction in subroutine CHKMSG.
C
C      12.  Modified PLOTFILE output to include date field.
C
C      13.  Modifications to some debug output statements based on
C           code provided by ENSR.
C
C       MODIFIED BY:    Roger W. Brode
C                       MACTEC Federal Programs, Inc.
C                       (formerly known as PES, Inc.)
C                       October 26, 2004
C
C       MODIFIED FROM:          AERMOD
C                       (Version Dated 04079)
C
C=======================================================================
C
C       This revised DRAFT version (dated 04079) incorporates
C       modifications to the wet deposition algorithms for both
C       gaseous and particle emissions.  For both gaseous and particle
C       wet deposition, the wet fluxes have been corrected to include
C       a factor of 3600.*SQRT(2*PI) in the denominator.  The factor
C       of 3600 was needed to correct a unit conversion error between
C       seconds and hours in the final calculation of the flux.
C       The factor of SQRT(2*PI) is needed to complete the integrated
C       vertical term.  A problem causing potential runtime errors for
C       volume and area sources with dry depletion was also corrected.
C
C       In addition to the corrections identified above, the particle
C       wet deposition algorithms were also modified to include an
C       algorithm for calculating the collision efficiency as a function
C       of particle size and raindrop size.  The previous version of
C       the model included a fixed value of 4.0e-4 for the collision
C       efficiency.
C
C
C       MODIFIED BY:    Roger W. Brode
C                       MACTEC Federal Programs, Inc.
C                       (formerly known as PES, Inc.)
C                       March 19, 2004
C
C       MODIFIED FROM:           AERMOD
C                         (Version Dated 03273)
C
C=======================================================================
C
C       This DRAFT version (dated 03273) incorporates wet and dry
C       deposition algorithms based on the draft ANL report (Wesely,
C       et. al, 2001), with modifications to the wet deposition
C       algorithms based on peer review comments.  The dry deposition
C       algorithms include dry depletion based on the simple source
C       depletion method.
C
C       NOTE:  The wet SCIM'ing and output by particle size options 
C       from the ISCST3 model have not been implemented yet in AERMOD.
C
C       This version includes the following modifications relative to
C       the previous draft (dated 03213):
C
C       1.  Removed depletion for the "inside cavity source" from the
C           PRIME calculations.
C
C       2.  Moved the code to adjust for TS < TA (used to model a fixed
C           delta TS-TA) from SUBROUTINE SETSRC back to SUBROUTINE
C           FLUXES.  This corrects some minor discrepancies between
C           the consequence analysis results for non-buoyant sources
C           relative to version 02222.
C
C       3.  Removed the dry particle deposition code associated with
C           DFAULT mode in ISCST3.  The only dry particle deposition
C           in AERMOD is based on the ANL report for Methods 1 and 2.
C           This also corrects a logic problem if neither TOXICS nor
C           DFAULT options are specified.
C
C       4.  Modified FUNCTION F2INT to only call DELTAH for point
C           sources during plume depletion calculation.
C
C       5.  Added identification of urban sources and Method 2 sources
C           in summary of source inputs.
C
C       6.  Included optional dry depletion option for Method 2.
C
C       7.  Additional code cleanup and documentatino, including removal
C           of unused data arrays associated with ISCST3 depletion code
C           and moving the call to SUBROUTINE METINI to follow the call
C           to SETSRC in SUBROUTINE OCALC for open pit sources.
C
C
C       Version 03213 (August 1, 2003) included the following
C       modifications relative to the previous draft (dated 03171):
C
C       1.  Corrects problem with calculation of f2 term used in gas
C           dry deposition.  The calculation of Wnew and f2 had to be
C           moved outside the source loop to properly account for
C           accumulated precipitation over three hour period.
C
C       2.  Changed the definition of the "top of the plume" to be based
C           on the plume centerline height plus 2.15*sigma-z, evaluated
C           at 20km downwind.  The previous definition was based on
C           3.9*sigma-z.
C
C       3.  Terrain effects are now incorporated in the dry depletion
C           calculation.  The terrain elevation is linearly interpolated
C           between the source base elevation and the terrain elevation
C           at the receptor.  The hill height scale is linearly
C           interpolated between the stack release height at the
C           source and the hill height scale from AERMAP at the receptor
C           location.
C
C       4.  Modified the short-term EVENT processing option to be
C           compatible with the deposition algorithms.  If more than
C           one output type is selected in the normal model run that
C           generates the EVENT input file, the events will be defined
C           based on the first output type, in the order of CONC, DEPOS,
C           DDEP, WDEP.
C
C       5.  Modified summary of first 24-hours of met data to include
C           additional parameters for deposition applications.
C
C       6.  Modified the program to use free-formatted READ for the
C           surface meteorological input file for all cases.  Decision
C           on whether to read additional parameters needed for deposition
C           is based on logical variables associated with deposition
C           calculations.
C
C       7.  Incorporated a patch in SUBROUTINE UNLUMP of PRIME.FOR to
C           avoid potential math error for downwash calculations.
C           The plume temperature calculated by NUMRISE is limited to
C           be greater than or equal to the ambient temperature minus
C           10 K.  This avoids a potential SQRT of a negative number.
C
C       MODIFIED BY:    Roger W. Brode
C                       MACTEC Federal Programs, Inc.
C                       (formerly known as PES, Inc.)
C                       September 30, 2003
C
C       MODIFIED FROM:           AERMOD
C                         (Version Dated 02222)
C
C=======================================================================
C
C       This DRAFT version (dated 02222) includes fixes to the following
C       bugs:
C
C       1) modification to ACALC to avoid potential math errors
C          for AREAPOLY sources;
C       2) correction to METSUM to output missing temperatures correctly
C          for the SCIM option;
C       3) correction to EMVARY to replace 'STAR' option with
C          'WSPEED' option;
C       4) modified CAV_SRC to keep "outside" cavity source in array
C          element 3 for cases when no "inside" cavity source contribution
C          occurs;
C       5) corrected meander algorithm to combine "plume" and
C          "pancake" components of concentrations rather than just
C          blending the lateral dispersion term, removed limit on pancake
C          term to be smaller than plume lateral term, and removed meander
C          from the PRIME component for sources subject to building
C          downwash;
C       6) corrected problem in NUMRISE to avoid referencing an
C          undefined variable (xbi);
C       7) corrected calling arguments for call to WAKE_SIG from
C          subroutine WAKE_DFSN, to use wakiz and wakiy instead of
C          turbz, and turby; and
C       8) correction of Z_iuo in subroutine URBCALC from 500 to 400
C          meters.
C
C       Modified to remove command line arguments for specifying
C       input and output file names, and use hardwired names of
C       AERMOD.INP and AERMOD.OUT.
C
C       This version includes an adjustment to ustar and L for urban stable
C       cases, by equating the "convective" sigma-w based on the urban
C       "convective" w* with the mechanical sigma-w based on u*
C       evaluated at a height of 7 times the urban roughness length.
C       The URBANOPT keyword was modified to allow the user to input
C       the urban roughness length as an optional parameter following
C       the optional city name.  If no urban roughness length is input,
C       then the model assumes an urban roughness length of 1.0 meter.
C
C       This version also includes a modification to the minimum layer
C       depth near the ground used to calculate effective parameters.
C       A minimum layer depth of 5 meters is used instead of 2 meters
C       if the plume centroid height and receptor height are both below
C       5 meters.
C
C       MODIFIED BY:    Roger W. Brode
C                       PES, Inc.
C                       September 10, 2002
C
C       MODIFIED FROM:           AERMOD
C                         (Version Dated 01247)
C
C=======================================================================
C
C       This draft version (dated 01247) includes the PRIME building
C       downwash algorithms based on the ISC-PRIME model (dated 99207).
C       For cases involving building downwash, the model calculates a
C       non-wake contribution using the AERMOD algorithms, a wake
C       contribution using the PRIME algorithms, and blends the two
C       results using a factor called GAMFACT that varies based on the
C       location of the receptor relative to the wake.  For receptors
C       within the wake region, where the lateral and vertical boundaries
C       are defined by the wake half-width and height, respectively, and
C       the longitudinal boundary is defined by a distance equal to 15R
C       or the point where wake turbulence intensity decays to ambient
C       turbulence intensity, whichever is greater, measured from the
C       upwind edge of the building, GAMFACT is set equal to 1.0 (i.e.,
C       uses the PRIME result only).  The PRIME algorithm has been
C       modified to use the AERMOD meteorological profiles and
C       definitions of ambient turbulence.
C
C       This version also includes a modification to the DTHETA/DZ 
C       profile (TEMPGRID.FOR) for extrapolating above the highest
C       measurement height for cases with observed temperature profiles,
C       a modification to the upper limit on the integration for
C       HCRIT (CALC2.FOR), a correction to the calculation of FYPAN
C       for meander (CALC2.FOR), and an adjustment to ISTRT_WIND used
C       for Y2K compliance to subtract 1 from ISTRT_WIND in case the
C       meteorological data file contains data from the end of the
C       previous year.
C
C       MODIFIED BY:    Roger W. Brode
C                       PES, Inc.
C                       September 4, 2001
C
C       MODIFIED FROM:           AERMOD
C                         (Version Dated 00357)
C
C=======================================================================
C
C       This version (dated 00357) includes enhancements based on the
C       current ISCST3 model (dated 00101).  These include:  1) the use of
C       globally allocatable arrays for data storage; 2) expanded data
C       structures to allow for output of concentration and deposition in
C       a single model run (for use when deposition algorithms are added
C       to AERMOD); 3) EVENT processing for short-term culpability analyses;
C       4) post-1997 PM10 processing; 5) TOXICS option enhancements such as
C       optimizations for area sources, the SCIM option, and SEASONHR output
C       file option; 6) explicit treatment of multiple-year meteorological
C       data files and ANNUAL averages; 7) the SHRDOW and SHRDOW7 options
C       for specifying emissions that vary by season, hour-of-day, and
C       day-of-week; and 8) improved data structures for field length and
C       filename lengths.
C
C       The following modifications have also been made to correct errors:
C       1) the PARAMETER LAMDAY was implicitly typed as integer in previous
C       versions, and is now explicitly typed as real, which has an impact
C       on the plume height calculations for the indirect source; 2) the
C       calculation of the plume centroid height for unstable conditions
C       for area sources was moved from SUB. ACALC to be included inside
C       the area source integration in SUB. PLUMEF and SUB. PWIDTH, which
C       affects results for area sources during unstable conditions,
C       especially for receptors located inside an elongated area source
C       with the wind blowing along the long dimension of the source; and
C       3) the STABLE and UNSTAB logical variables are now assigned prior
C       to the call to COMPTG in SUB. METEXT, potentially affecting the
C       observed DTHETA/DZ profile for the first stable hour in a day;
C       4) the calculation of wind direction (WDIR) at stack top in SUB.
C       METINI and at stack height plus 0.5*deltaH in SUB. PCALC was
C       corrected to account for possible 0-360 crossover in the profile;
C       5) corrections were made to the vertical terms to correct the
C       reflection component for receptors located below stack base elevation;
C       6) corrected calculation of FYPAN term for meander in SUB. FYTERM
C       to use radial distance (DISTR) instead of downwind distance (X); and
C       7) application of meander to both stable and unstable conditions.
C       Additional modifications were made to consolidate redundant code
C       and simplify future maintenance activities.
C
C       MODIFIED BY:    Roger W. Brode
C                       PES, Inc.
C                       December 22, 2000
C
C       MODIFIED FROM:           AERMOD
C                         (Version Dated 99351)
C
C========================================================================
C
C                              (Version Dated 99351)
C                                December 17, 1999
C
C        This version (99351) includes the following corrections to the
C        implementation of the Schulman-Scire downwash algorithm:  1) added
C        call to DHPSS in calculation of plume centroid height (CENTER) in
C        subroutine PCALC; and 2) modified to use SZ3LB and SY3LB based on
C        building enhanced dispersion curves only in subroutine DHPSS.
C        Additional modifications were made to improve consistency with
C        ISCST3 implementation of Schulman-Scire downwash algorithm.
C        Also includes changes to subroutine PCCODE to facilitate
C        compilation of the model using the DEC Visual Fortran compiler.
C        The output file unit number, IOUNIT, was also changed from 6 to 9
C        in order for runtime status update to appear on the screen for
C        DEC-compiled executables.  Minor, inconsequential changes were
C        also made to comment headers and variable declarations in
C        SIGGRID.FOR.
C
C
C        MODIFIED FROM:
C                              (Version Dated 99211)
C                                  July 30, 1999
C
C========================================================================
C
C        This version (99211) incorporates modifications for Y2K compliance.
C        Uses a window of 1950 to 2049 for 2-digit years.  Will utilize
C        4-digit year if input for surface and profile files using FREE
C        format (the default read format still reads a 2-digit year).
C        Changes also include calculation of a 10-digit date variable
C        (FULLDATE) with 4-digit year for date comparisons, and changes of
C        the output formats for the 8-digit variable, KURDAT, to I8.8 to
C        include leading zeros.  The date and time routines used for the
C        page headers have been modified to use the standard Fortran 90
C        routines, and minor changes have been made to remove obsolescent
C        features from the code.  Changes also include a correction to
C        a variable name in SUBROUTINE DELTAH.
C
C
C        MODIFIED FROM:
C                              (Version Dated 98314)
C                                November 10, 1998
C
C========================================================================
C
C        This version (dated 98314) incorporates modifications for the final
C        draft AERMOD Model. This draft represents the final version of the
C        model prior to the Notice of Proposed Rulemaking for including
C        AERMOD in the modeling Guideline and subsequent public comment period.
C        Changes to the model are too numerous to list here in detail.
C        Changes to the interface are incorporated in the revised
C        AERMOD user's guide, and pertain primarily to the regulatory
C        default option on the CO MODELOPT card, and removal of the
C        developmental options for terrain affects.  The CO TERRHGTS
C        keyword is now obsolete, and a new mandatory ME PROFBASE
C        keyword has been added for inputting the base elevation above
C        MSL for the gridded potential temperature profile.
C
C
C        MODIFIED FROM:
C                               (Version Dated 98022)
C                                  January 22, 1998
C
C========================================================================
C
C        This version (dated 98022) incorporates modifications for the revised
C        draft AERMOD Model. This draft represents the final Phase I version
C        of the model.  This version includes the use a minimum sigma-v of
C        0.2 m/s. It also incorporates corrections to the downwash algorithms.
C        It limits the number of iterations on inhomogeneity to 1,
C        i.e., effective parameters are calculated based on an average
C        for the layer from plume centerline to 2.15 sigma-z, where
C        sigma-z is based on parameters at plume centerline height.
C        The transport wind direction is based on the modpoint between
C        stack height and "final" plume height.  Modifications have also
C        been made to the dtheta/dz profile, the Tly and Tlz used for stable
C        plumes above the CBL, and the height of the effective reflecting
C        surface for stable plumes.
C
C
C        MODIFIED FROM:
C                               (Version Dated 97350)
C                                 December 16, 1997
C
C========================================================================
C
C        This draft includes an option for specifying the transport wind
C        direction.  The option is specified on the MODELOPT card, where
C        WDOPT1 is for wind direction taken at stack height, and
C        WDOPT2 is for wind direction at the midpoint between stack height
C        and "final" plume height.  The default option is WDOPT1.
C        This version does not include conditional compilation code to
C        support the Microsoft Fortran compiler.  This version of the
C        code is compatible with the Lahey F77L3 and Lahey LF/90 compilers.
C        The deposition algorithm for AERMOD is still under development,
C        and is not operational in this version.
C                   R.W. Brode, PES, Inc., December 16, 1997
C
C
C        MODIFIED FROM:
C                               (Version Dated 97064)
C                                   March 5, 1997
C
C========================================================================
C
C        Incorporates modifications for second round of Beta Testing,
C        including code clean-up and removing obsolete options.  Since
C        the urban stable boundary layer algorithm in AERMOD is still
C        under development, the keywords for implementing the urban
C        option in AERMOD have been disabled for this revised draft
C        Beta release of the model.
C                   R.W. Brode, PES, Inc., March 5, 1997
C
C        MODIFIED FROM:
C                               (Version Dated 96239)
C                                   August 26, 1996
C
C========================================================================
C
C        Includes modification to wind direction gridded profiling in
C        SUB. GRDWD, correcting a problem with wind directions backing
C        through 360 degrees.  This could result in a level of observed
C        wind direction being erroneously counted as missing.
C                   R.W. Brode, PES, Inc., August 26, 1996
C
C        MODIFIED FROM:
C                               (Version Dated 96228)
C                                   August 15, 1996
C
C========================================================================
C
C        Includes modifications for low wind/low turbulence cases.  Also
C        changed definition of "valid lower bound" to be 7zo instead of
C        20zo for consistency with recent changes to AERMET.  The terrain
C        options have also been reinstated on the CO MODELOPT card as
C        TERROPT1 and TERROPT2.  The command line developmental option
C        switch has been enabled in this version.
C                   R.W. Brode, PES, Inc., August 15, 1996
C
C        MODIFIED FROM:
C
C                               (Version Dated 96198)
C                                    July 16, 1996
C
C========================================================================
C
C        Includes Urban Boundary Layer Option based on the Model Coding
C        Abstract by Akula Ventkatram dated 4/1/96.  Modified inputs
C        include two new keywords:
C
C                 CO URBANOPT  Urbpop  (Urbnam)
C                    where Urbpop is the population of the urban area, and
C                          Urbnam is an optional character field
C                          for the name of the urban area.
C
C                 SO URBANSRC  Srcid's  and/or  Srcrng's
C                    where Srcid's identifies individual sources to be
C                          modeled as urban sources, and Srcrng's
C                          identifies a range of sources to be modeled
C                          as urban.
C
C
C        MODIFIED FROM:
C                               (Version Dated 96131)
C                                    May 10, 1996
C
C        Includes OPTG3 and OPTG4 for stable plume reflection options.
C
C        MODIFIED FROM:
C                               (Version Dated 96053)
C                                 February 22, 1996
C
C        Includes flow vector, AFV, in the EVALFILE output.
C
C        MODIFIED FROM:
C                               (Version Dated 96046)
C                                 February 15, 1996
C
C        Includes new sigma-v and sigma-w profiles coded by Bob Paine,
C        some modifications to TEMPGRID.FOR to avoid discontinuities in
C        the VPTG profile, reinstates command-line input for developmental
C        options for the Lahey version, modifies the default option settings
C        to remove stable plume reflections, reinstates the original stable
C        profile for TLz, increases the maximum number of iterations in
C        IBLVAL from 5 to 20, and includes patches for smoothed h < 0 and
C        for mixing heights at (or near) 5,000m.
C                   R.W. Brode, PES, Inc., February 15, 1996
C
C        MODIFIED FROM:
C                               (Version Dated 95272)
C
C        Includes a SCREEN mode option on CO MODELOPT card,
C        addition of AREAPOLY and AREACIRC source types,
C        and INCLUDED keyword option for including data from
C        an external file for the SO and RE pathways.  The
C        INCLUDED option is intended as a link to AERMAP and
C        for use with the screening version of AERMOD.
C
C        MODIFIED FROM:
C                               (Version Dated 95188)
C
C        Hardcoded option settings for the Beta release of AERMOD.
C
C        MODIFIED FROM:
C                               (Version Dated 95066)
C
C        "Optionized" version of AERMOD for Developmental Evaluation.
C        Developmental options are selected by use of an additional
C        10-character command line argument, e.g. '1213121111'.
C        If the additional command line argument is not present, then
C        the model will default to the Base Model.  This is equivalent
C        to using an additional command line argument of '1111111111'.
C        Note that the single quotes are not included on the command line.
C                   R.W. Brode, PES, Inc., January 27, 1995
C
C        Version 95066 includes a few minor fixes, some modifications to
C        the EVALFILE output, and the volume source option (mostly untested).
C        The fixes affect the SBL plume reflection option (OPTG2), and two of
C        the inhomogeneity options (OPTD2 and OPTD3).  In the latter cases,
C        the effective wind speed was allowed to be less than the effective
C        sigma-w.  Also, included EVALCART receptors with DISCCART receptors
C        for output purposes (RECTABLE, DAYTABLE, etc., and INPSUM.FOR).
C                   R.W. Brode, PES, Inc., February 16, 1995
C
C        Base Case Model for AERMOD Developmental Evaluation - 12/8/94
C
C        MODIFIED FROM:    ISC2 Short Term Model - ISCST2
C                               (Version Dated 93109)
C
C        MODIFIED FROM:         (Version Dated 92273)
C
C        MODIFIED FROM:         (Version Dated 92062)
C
C        PURPOSE: Controls Overall Flow and Processing of ISCST2 Model
C
C        PROGRAMMED BY: Roger W. Brode
C                       James O. Paumier
C                       Jayant A. Hardikar
C                       Pacific Environmental Services, Inc.
C                       P.O. Box 12077
C                       Research Triangle Park, North Carolina  27709
C
C        DATE:    November 9, 1993
C
C        INPUTS:  Command Line Options
C
C        OUTPUTS: Model Results
C
C========================================================================