#include "DARWIN_OPTIONS.h"

CBOP
C !ROUTINE: DARWIN_SPEAD_WEIGHT
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_SPEAD_WEIGHT(
     I     myTime,myIter,myThid,
     I     mn_tr,sd_tr,
#ifdef SPEAD_CONTINUOUS_COVARIANCE
     I     cv_tr,
#endif
#ifdef DARWIN_DEBUG
     I     iG, jG, k, dT,
#endif
     O     step_tr,
     O     wght)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_INDICES.h"
#include "DARWIN_TRAITS.h"

C !INPUT PARAMETERS: ===================================================
C  myTime :: current time
C  myIter :: current iteration number
C  mn_tr  :: mean trait values
C  sd_tr  :: trait standard deviations
C  cv_tr  :: trait covariances
      INTEGER myThid, myIter
      _RL myTime
      _RL mn_tr(nPhyp,nTrait), sd_tr(nPhyp,nTrait)
#ifdef SPEAD_CONTINUOUS_COVARIANCE
      _RL cv_tr(nPhyp,nCov)
#endif
#ifdef DARWIN_DEBUG
      INTEGER iG, jG, k
      _RL dT
#endif

C !OUTPUT PARAMETERS: ==================================================
C  wght :: weight of each phytoplankton tracer
C          Used in the derivatives of non-living tracers
C          Avoid explicitly manipulating derivatives of tracer fluxes
C          with respect to trait
C          Also secure mass conservation
      _RL step_tr(nPhyp,nTrait)
      _RL wght(nTrac)

c !LOCAL VARIABLES: ====================================================
      INTEGER jp, jt, jt1, jt2, jt3, jn, jn1, jn2, jc
C Error messages for debugging (Le Gland, 16/03/2021)
      CHARACTER*(MAX_LEN_MBUF) msgBuf

C Trait steps to estimate derivatives (Le Gland, 10/03/2021)
      _RL dbv, dto, dpo

C Minimum value for optimal irradiance at a given cell volume
C     _RL PARoptmin

C Trait value for each tracer and each variable trait (Le Gland, 31/03/2021)
      _RL vl_tr(nTrac,nTrait)

      _RL tmp

C Derivatives computed numerically to reduce code complexity
C The number of call to the biogeochemical model is equal to the number of tracers

C     WRITE(msgBuf,'(A)') 'DARWIN_SPEAD_RATES: variables declared:'
C     CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
C    &                    SQUEEZE_RIGHT, myThid )
C     print*,'DARWIN_SPEAD_RATES: variables declared'


C-----------------------------------------------------------------------
C Step 1: Set weight and master trait values for each tracer
C-----------------------------------------------------------------------

C Steps to compute numerical trait derivatives
      dbv = 0.1
      dto = 0.1
      dpo = 0.01

C Traits used to compute trait derivatives
      DO jp = 1, nPhyp

        jt = 1
        IF (traitbvol(jp) .NE. 0) THEN
          vl_tr(iplank(jp),jt) = mn_tr(jp,jt)
#ifdef SPEAD_DERIVATIVES_LOCAL
          step_tr(jp,jt) = dbv
#else
          CALL RANDOM_NUMBER(tmp)
          step_tr(jp,jt) = SIGN(MAX(dbv,sd_tr(jp,jt)),tmp-0.5)
#endif
          jt = jt + 1
        ENDIF
        IF (traittopt(jp) .NE. 0) THEN
          vl_tr(iplank(jp),jt) = mn_tr(jp,jt)
#ifdef SPEAD_DERIVATIVES_LOCAL
          step_tr(jp,jt) = dto
#else
          CALL RANDOM_NUMBER(tmp)
          step_tr(jp,jt) = SIGN(MAX(dto,sd_tr(jp,jt)),tmp-0.5)
#endif
          jt = jt + 1
        ENDIF
        IF (traitparopt(jp) .NE. 0) THEN
          vl_tr(iplank(jp),jt) = mn_tr(jp,jt)
#ifdef SPEAD_DERIVATIVES_LOCAL
          step_tr(jp,jt) = dpo
#else
          CALL RANDOM_NUMBER(tmp)
          step_tr(jp,jt) = SIGN(MAX(dpo,sd_tr(jp,jt)),tmp-0.5)
#endif
        ENDIF

        jn = iplank(jp) + 1
        DO jt = 1, num_trait(jp)
          vl_tr(jn,jt)   = mn_tr(jp,jt) - step_tr(jp,jt)
          vl_tr(jn+1,jt) = mn_tr(jp,jt) + step_tr(jp,jt)
          DO jt2 = 1, num_trait(jp)
            IF (jt2 .NE. jt) THEN
              vl_tr(jn,jt2)   = mn_tr(jp,jt2)
              vl_tr(jn+1,jt2) = mn_tr(jp,jt2)
            ENDIF
          END DO
          wght(jn) = (1./2.) * (sd_tr(jp,jt)/step_tr(jp,jt))**2
          wght(jn+1) = wght(jn)
          jn = jn + 2
        ENDDO

#ifdef SPEAD_CONTINUOUS_COVARIANCE
        DO jt1 = 1, num_trait(jp)-1
          DO jt2 = jt1+1, num_trait(jp)
            vl_tr(jn,jt1) = mn_tr(jp,jt1) + step_tr(jp,jt1)
            vl_tr(jn,jt2) = mn_tr(jp,jt2) + step_tr(jp,jt2)
            DO jt3 = 1, num_trait(jp)
              IF ((jt3 .NE. jt1) .AND. (jt3 .NE. jt2)) THEN
                vl_tr(jn,jt3) = mn_tr(jp,jt3)
              ENDIF
            ENDDO
            jc = NINT((jt1-1)*(num_trait(jp)-1-jt1/2.) + jt2 - 1)
            wght(jn) = cv_tr(jp,jc) / (step_tr(jp,jt1)*step_tr(jp,jt2))
            jn1 = iplank(jp) + 2*jt1
            jn2 = iplank(jp) + 2*jt2
            wght(jn1) = wght(jn1) - wght(jn)
            wght(jn2) = wght(jn2) - wght(jn)
C           print*, 'Covariance loop  ', 'num_trait:', num_trait(jp),
C    &              'jc:', jc, 'jn:', jn, 'vr_tr:', vl_tr(jn,:)
            jn = jn + 1
          ENDDO
        ENDDO
#endif

        wght(iplank(jp)) = 1.0 _d 0 - SUM(wght(iplank(jp)+1:jn-1))
      ENDDO

C-----------------------------------------------------------------------
C Step 2: Set the values of all traits for all tracers
C-----------------------------------------------------------------------

      DO jp = 1, nPhyp
        DO jn = iplank(jp), iplank(jp)+2*num_trait(jp)+num_cov(jp)
          jt = 1
          IF (traitbvol(jp) .NE. 0) THEN
            biovol(jn) = EXP(vl_tr(jn,jt))
            jt = jt + 1
          ENDIF
          IF (traittopt(jp) .NE. 0) THEN
            phytoTempOptimum(jn) = vl_tr(jn,jt)
            jt = jt + 1
          ENDIF
          IF (traitparopt(jp) .NE. 0) THEN
            PARopt(jn) = EXP(vl_tr(jn,jt))
C           ksatPAR(jn) = LOG(1+PARchi(jp))
C    &                  / (11.574*PARopt(jp)/2.5)
CC   &                  / (darwin_inscal_PAR*PARopt(jp)/2.5)
C           kinhPAR(jp) = ksatPAR(jp) / PARchi(jp)
C convert to µEin/m2/s (Le Gland, 27/06/2024)
C           PARopt(jp) = PARopt(jp) * 11.574 / 2.5
C Maximum growth rate depends on optimal irradiance
C           g = group(jp)
C           PARoptmin = a_PARoptmin(g) * biovol(jp)**b_PARoptmin(g)
C           PARoptmin = PARoptmin * 11.574 / 2.5
C           PCmax(jp) = MAX(0., a_PCmax(g) * biovol(jp)**b_PCmax(g)
C    &                * (1.0 - PARoptmin/PARopt(jp))
C    &                * (1.0 - (a_biovolmin(g)/biovol(jp)) ) )
C         ELSE
C           PCmax(jp) = MAX(0., a_PCmax(g) * biovol(jp)**b_PCmax(g) )
          ENDIF
#ifdef DARWIN_DEBUG
C         IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG) THEN
            print*, 'species', jp, 'tracer', jn, 'trait', jt,
     &              'biovol', biovol(jn), 'topt', phytoTempOptimum(jn),
     &              'PARopt', PARopt(jn)
C         ENDIF
#endif
        END DO
      END DO

      DO jp = nPhoto+1, nTrac
        wght(jp) = 1.0 _d 0
      ENDDO

#ifdef DARWIN_DEBUG
C     IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG) THEN
        print*,'Weight_derivatives', wght, SUM(wght(1:nTrac))
C     ENDIF
#endif

#ifdef SPEAD_CONTINUOUS_TRAIT
      CALL DARWIN_GENERATE_ALLOMETRIC(myThid)
#endif

      RETURN
      END SUBROUTINE
