#include "DARWIN_OPTIONS.h"

CBOP
C !ROUTINE: DARWIN_GENERATE_INDICES

C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_GENERATE_INDICES( myThid )

C !DESCRIPTION:
C     Generate parameters for plankton types using an allometric approach

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_INDICES.h"
#include "DARWIN_RADTRANS.h"
#include "DARWIN_PARAMS.h"
#include "DARWIN_TRAITPARAMS.h"
#include "DARWIN_TRAITS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid
CEOP

#ifdef ALLOW_DARWIN

C !LOCAL VARIABLES: ====================================================
C     msgBuf   :: Informational/error meesage buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER jp,jp2,jp3,jp4,jpm,jps,jz,g,nmom,ntrt
#ifndef SPEAD_CONTINUOUS_TRAIT
      _RL logvol
#endif

C In a continuous model, this subroutine only sets the traits that do not
C vary within a group (Le Gland, 28/01/2021)
C Is this cumbersome separation really necessary ?

C ======================================================================
C compute traits from trait parameters

C     WRITE(msgBuf,'(A)') 'DARWIN_GENERATE_CONSTANT: entering:'
C     CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
C    &                    SQUEEZE_RIGHT, myThid )
C     print*,'DARWIN_GENERATE_CONSTANT: entering:'

#ifdef SPEAD_CONTINUOUS_TRAIT

      jpm = 1
      jps = 1
      print*,'--------------------------------------'
      print*,'DARWIN_GENERATE_INDICES:'
      DO g=1,nGroup
C The number of traits can be set for each species
        ntrt = grp_traitbvol(g)+grp_traittopt(g)+grp_traitparopt(g)
        nmom = nint(1 + 1.5*ntrt + 0.5*ntrt**2)
        print*,'group',g,'ntrt',ntrt,'nmom',nmom
        DO jp=1,grp_nplank(g)
C First moment index of species jps
          iplank(jps) = jpm
          num_trait(jps) = ntrt
#ifdef SPEAD_CONTINUOUS_COVARIANCE
          num_cov(jps)   = nint(0.5*ntrt*(ntrt-1))
#else
          num_cov(jps)   = 0
#endif
          DO jp2=1,nmom
          print*,'group',g,'species',jps,'moment',jpm
            group(jpm) = g
C Species index "plank" is also necessary in darwin_spead_rates.F
            plank(jpm) = jps
            jpm = jpm+1
            biovol(jpm) = biovol0(g)
            phytoTempOptimum(jpm) = Topt0(g)
            PARopt(jpm) = PARopt0(g)
          END DO
          jps = jps+1
        END DO
      END DO
      print*,'--------------------------------------'

#else /* NOT spead_continuous_trait */

C compute cell volumes in micrometer^3
C
C in decreasing precedence (if bold quantity is set):
C
C   V = GRP_BIOVOL(jp,g)
C   V = 10**(logvolbase+(GRP_BIOVOLIND(jp,g)-1)*logvolinc)
C   V = 10**(logvolbase+(LOGVOL0IND(g)+jp-2)*logvolinc)
C   V = BIOVOL0(g)*biovolfac(g)**(jp-1)
C
C     if logvol0ind is set, use it to compute biovol0
      DO g=1,nGroup
        IF (logvol0ind(g) .GT. 0) THEN
          IF (biovol0(g) .NE. 0.0 _d 0) THEN
            WRITE(msgBuf,'(2A)') 'DARWIN_GENERATE_INDICES: ',
     &      'cannot set both biovol0 and logvol0ind'
            CALL PRINT_ERROR( msgBuf, myThid )
            STOP 'ABNORMAL END: S/R DARWIN_GENERATE_INDICES'
          ENDIF
          logvol = logvolbase + (logvol0ind(g)-1)*logvolinc
          biovol0(g) = 10 _d 0 ** logvol
          biovolfac(g) = 10 _d 0 ** logvolinc
        ENDIF
      ENDDO

      DO g=1,nGroup
       DO jp=1,grp_nbiovol(g)
        IF (grp_logvolind(jp,g) .GT. 0 .AND.
     &      grp_biovol(jp,g) .GT. 0.0 _d 0) THEN
          WRITE(msgBuf,'(2A)') 'DARWIN_GENERATE_INDICES: ',
     &    'cannot set both grp_biovol and grp_logvolind'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R DARWIN_GENERATE_INDICES'
        ELSEIF (grp_logvolind(jp,g) .GT. 0) THEN
C         pick from logrange
          logvol = logvolbase + (grp_logvolind(jp,g)-1)*logvolinc
          grp_biovol(jp,g) = 10 _d 0 ** logvol
        ELSEIF (grp_biovol(jp,g) .EQ. 0) THEN
          IF (biovol0(g) .GT. 0.0 _d 0) THEN
           grp_biovol(jp,g) = biovol0(g) * biovolfac(g)**(jp-1)
          ELSE
           WRITE(msgBuf,'(3A)') 'DARWIN_GENERATE_INDICES: ',
     &     'Need to set one of grp_biovol, grp_logvolind, ',
     &     'biovol0, logvol0ind'
           CALL PRINT_ERROR( msgBuf, myThid )
           STOP 'ABNORMAL END: S/R DARWIN_GENERATE_INDICES'
          ENDIF
        ENDIF
        biovol_bygroup(jp,g) = grp_biovol(jp,g)
       ENDDO
       DO jp=grp_nbiovol(g)+1,nplank
        IF (grp_biovol(jp,g) .NE. 0.0 _d 0) THEN
          WRITE(msgBuf,'(2A,I5,A,I5,A)') 'DARWIN_GENERATE_INDICES: ',
     &    'index',jp,'group',g,'out of range for grp_biovol'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R DARWIN_GENERATE_INDICES'
        ENDIF
        IF (grp_logvolind(jp,g) .NE. 0) THEN
          WRITE(msgBuf,'(2A,I5,A,I5,A)') 'DARWIN_GENERATE_INDICES: ',
     &    'index',jp,'group',g,'out of range for grp_logvolind'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R DARWIN_GENERATE_INDICES'
        ENDIF
        biovol_bygroup(jp,g) = 0.0 _d 0
       ENDDO
      ENDDO

C Compute optimal temperature in degree Celsius
C
C in decreasing precedence (if bold quantity is set):
C
C   Topt = GRP_TOPT(jp,g)
C   Topt = Topt0(g)+Toptstep(g)*(jp-1)

      DO g=1,nGroup
        DO jp=1,grp_ntopt(g)
          IF (grp_Topt(jp,g) .EQ. 0) THEN
            IF (Topt0(g) .GT. 0.0 _d 0) THEN
              grp_Topt(jp,g) = Topt0(g) + Toptstep(g)*(jp-1)
            ELSE
              WRITE(msgBuf,'(3A)') 'DARWIN_GENERATE_INDICES: ',
     &        'Need to set one of grp_Topt, Topt0, '
              CALL PRINT_ERROR( msgBuf, myThid )
              STOP 'ABNORMAL END: S/R DARWIN_GENERATE_INDICES'
            ENDIF
          ENDIF
        Topt_bygroup(jp,g) = grp_Topt(jp,g)
        ENDDO
        DO jp=grp_ntopt(g)+1,nplank
          IF (grp_Topt(jp,g) .NE. 0.0 _d 0) THEN
            WRITE(msgBuf,'(2A,I5,A,I5,A)') 'DARWIN_GENERATE_INDICES',
     &      ': index',jp,'group',g,'out of range for grp_Topt'
            CALL PRINT_ERROR( msgBuf, myThid )
            STOP 'ABNORMAL END: S/R DARWIN_GENERATE_INDICES'
          ENDIF
          Topt_bygroup(jp,g) = 0.0 _d 0
        ENDDO
      ENDDO

C Compute optimal irradiance in W/m2
C
C in decreasing precedence (if bold quantity is set):
C
C   PARopt = GRP_PAROPT(j,g)
C   PARopt = PARopt0(g)+PARoptfac(g)*(j-1)

      DO g=1,nGroup
        DO jp=1,grp_nparopt(g)
          IF (grp_PARopt(jp,g) .EQ. 0) THEN
            IF (PARopt0(g) .GE. 0.0 _d 0) THEN
              grp_PARopt(jp,g) = PARopt0(g) * PARoptfac(g)**(jp-1)
            ELSE
              WRITE(msgBuf,'(3A)') 'DARWIN_GENERATE_INDICES: ',
     &        'Need to set one of grp_PARopt, PARopt0, '
              CALL PRINT_ERROR( msgBuf, myThid )
              STOP 'ABNORMAL END: S/R DARWIN_GENERATE_INDICES'
            ENDIF
          ENDIF
        PARopt_bygroup(jp,g) = grp_PARopt(jp,g)
        ENDDO
        DO jp=grp_nparopt(g)+1,nplank
          IF (grp_PARopt(jp,g) .NE. 0.0 _d 0) THEN
            WRITE(msgBuf,'(2A,I5,A,I5,A)') 'DARWIN_GENERATE_INDICES',
     &      ': index',jp,'group',g,'out of range for grp_PARopt'
            CALL PRINT_ERROR( msgBuf, myThid )
            STOP 'ABNORMAL END: S/R DARWIN_GENERATE_INDICES'
          ENDIF
          PARopt_bygroup(jp,g) = 0.0 _d 0
        ENDDO
      ENDDO

C Three trait sorting (Le Gland, 07/04/2021)
C First change PARopt, then Topt, then biovol, then group
      jp = 1
      DO g=1,ngroup
        DO jp2=1,grp_nbiovol(g)
          DO jp3=1,grp_ntopt(g)
            DO jp4=1,grp_nparopt(g)
              IF (jp .GT. nPlank) THEN
                WRITE(msgBuf,'(2A)') 'DARWIN_GENERATE_INDICES: ',
     &          'need SUM(grp_nbiovol*grp_ntopt) = nPlank,
     &           nPlank too small'
                CALL PRINT_ERROR( msgBuf, myThid )
                STOP 'ABNORMAL END: S/R DARWIN_GENERATE_INDICES'
              ENDIF
              group(jp)  = g
              igroup(jp) = jp2
              jgroup(jp) = jp3
              kgroup(jp) = jp4
              biovol(jp) = grp_biovol(jp2,g)
              phytoTempOptimum(jp) = grp_Topt(jp3,g)
              PARopt(jp) = grp_PARopt(jp4,g)
C             ksatPAR(jp) = LOG(1+PARchi(jp)) 
C    &                    / (11.574*PARopt(jp)/2.5)
C             kinhPAR(jp) = ksatPAR(jp) / PARchi(jp) 
C             ksatPAR(jp) = LOG(1+12.0) 
C    &                    / (11.574*PARopt(jp)/2.5)
C             kinhPAR(jp) = ksatPAR(jp) / 12.0 
C Necessary for compatibility with continuous model
              plank(jp) = jp
              iplank(jp) = jp
              jp = jp+1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF (jp .NE. nPlank + 1) THEN
        WRITE(msgBuf,'(2A,2I4)') 'DARWIN_GENERATE_INDICES: ',
     &    'need SUM(grp_nbiovol*grp_ntopt) = nPlank, not ',jp-1,nPlank
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R DARWIN_GENERATE_INDICES'
      ENDIF

#endif /* spead_continuous_trait */

#endif  /*ALLOW_DARWIN*/

      RETURN
      END
