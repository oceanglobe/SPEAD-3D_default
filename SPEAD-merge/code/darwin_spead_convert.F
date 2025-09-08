#include "DARWIN_OPTIONS.h"

CBOP
C !ROUTINE: DARWIN_SPEAD_CONVERT
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_SPEAD_CONVERT(
     U     Ptr,
     I     mode,myTime,myIter,myThid)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_INDICES.h"
#include "DARWIN_TRAITS.h"

#include "DARWIN_TRAITPARAMS.h"

C !INPUT PARAMETERS: ===================================================
C  mode   :: convert to trait (-1) or to conserved (1) moments
C  myTime :: current time
C  myIter :: current iteration number
C  myThid :: thread number
      INTEGER mode, myThid, myIter
      _RL myTime

C !INPUT/OUTPUT PARAMETERS: ============================================
C  pTr    :: darwin plankton model tracers
      _RL Ptr(nTrac)

c !LOCAL VARIABLES: ====================================================
      _RL X(nPlank)
      _RL mn(nTrait)
      _RL sd(nTrait)
      _RL maxcov
      INTEGER jp, jt, jt1, jt2, jc

#ifdef DARWIN_DEBUG
      print*, '----------  DARWIN_SPEAD_CONVERT:  ----------'
      print*, 'mode', mode, 'Ptr input', Ptr
C    &        'min_tr',min_tr, 'max_tr',max_tr, 'max_vr_tr',max_vr_tr,
C    &        'min_bvol',min_bvol
#endif

C Mode "-1" -> Convert conserved moments to trait moments
      IF (mode .EQ. -1) THEN
       DO jp = 1, nPlank
        X(jp) = MAX(0., Ptr(iplank(jp)))
        IF (X(jp) .GT. 0) THEN
         DO jt = 1, num_trait(jp)
          Ptr(iplank(jp)+jt) = MIN(max_tr(jp,jt),MAX(min_tr(jp,jt),
     &                         Ptr(iplank(jp)+jt)/X(jp)))
          Ptr(iplank(jp)+num_trait(jp)+jt) = MIN(max_vr_tr(jp,jt),
     &                      MAX(0.0,Ptr(iplank(jp)+num_trait(jp)+jt)
     &                      /X(jp) - Ptr(iplank(jp)+jt)**2))
          sd(jt) = sqrt(Ptr(iplank(jp)+num_trait(jp)+jt))
         ENDDO
#ifdef SPEAD_CONTINUOUS_COVARIANCE
         DO jt1 = 1, num_trait(jp)-1
          DO jt2 = jt1+1, num_trait(jp)
           jc = NINT((jt1-1)*(num_trait(jp)-1-jt1/2.)) + jt2 - 1
           maxcov = 0.99*sd(jt1)*sd(jt2)
           Ptr(iplank(jp)+2*num_trait(jp)+jc) = MIN(maxcov,MAX(-maxcov,
     &                         Ptr(iplank(jp)+2*num_trait(jp)+jc)/X(jp)
     &                       - Ptr(iplank(jp)+jt1)*Ptr(iplank(jp)+jt2)))
          ENDDO
         ENDDO
#endif
        ELSE
         DO jt = 1, 2*num_trait(jp)
          Ptr(jp+jt) = 0.0 _d 0
         ENDDO
#ifdef SPEAD_ALLOW_COV
         DO jc = 1, num_cov(jp)
          Ptr(jp+2*num_trait(jp)+jc) = 0.0 _d 0
         ENDDO
#endif
        ENDIF
       ENDDO

C Mode "1" -> Convert trait moments to conserved moments
      ELSEIF (mode .EQ. 1) THEN
       DO jp = 1, nPlank
        X(jp) = MAX(0., Ptr(iplank(jp)))
        DO jt = 1, num_trait(jp)
         mn(jt) = Ptr(iplank(jp)+jt)
         Ptr(iplank(jp)+jt) = mn(jt)*X(jp)
         Ptr(iplank(jp)+num_trait(jp)+jt) = X(jp)*
     &         (Ptr(iplank(jp)+num_trait(jp)+jt) + mn(jt)**2)
        ENDDO
#ifdef SPEAD_CONTINUOUS_COVARIANCE
        DO jt1 = 1, num_trait(jp)-1
         DO jt2 = jt1+1, num_trait(jp)
          jc = NINT((jt1-1)*(num_trait(jp)-1-jt1/2.) + jt2 - 1)
          Ptr(iplank(jp)+2*num_trait(jp)+jc) = X(jp)*
     &     (Ptr(iplank(jp)+2*num_trait(jp)+jc) + mn(jt1)*mn(jt2))
         ENDDO
        ENDDO
#endif
       ENDDO
      ENDIF

#ifdef DARWIN_DEBUG
      print*, 'Ptr output', Ptr
      print*, '---------------------------------------------'
#endif

      RETURN
      END SUBROUTINE
