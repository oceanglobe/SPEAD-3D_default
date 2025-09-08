#include "DARWIN_OPTIONS.h"

CBOP
C !ROUTINE: DARWIN_SPEAD_DERIVATIVES
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_SPEAD_DERIVATIVES(
     I     X,
     I     acom,gcom,
     I     a_1t,
     I     step_tr,
     I     vr_tr,sd_tr,
#ifdef SPEAD_CONTINUOUS_COVARIANCE
     I     cv_tr,cr_tr,
#endif
#ifdef DARWIN_GRAZING_SWITCH
     I     coeff_KTW,
#endif
#ifdef DARWIN_DEBUG
     I     iG, jG, k, dT,
#endif
     U     diags,
     U     gTr,
     I     myTime,myIter,myThid)

C !DESCRIPTION: CONVERT NET GROWTH RATES OF PLANKTON TRACERS INTO
C               DERIVATIVES WITH RESPECT TO TRAITS
C               THEN COMPUTE TENDENCIES OF TRAIT MOMENTS

C !USES: ===============================================================
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_INDICES.h"
#include "DARWIN_DIAGS.h"
#include "DARWIN_PARAMS.h"
#include "DARWIN_TRAITS.h"
#include "DARWIN_TRAITPARAMS.h"

C !INPUT PARAMETERS: ===================================================
C  X         :: biomass (mmolC/m3) of phytoplankton species
C  acom   :: net growth rate (gross growth - death) of each species
C  gcom   :: gross growth rate of each species
C  a_1t   :: net growth rate of each tracer (continuous-trait)
C  vr_tr  :: trait variances
C  sd_tr  :: trait standard deviations
C  cv_tr  :: trait covariances
C  cr_tr  :: trait correlations
C  coeff_KTW :: multiplicator of diversity under Kill-The-Winner grazing
C  myTime    :: current time
C  myIter    :: current iteration number
C  myThid    :: thread number
      _RL X(nPlank)
      _RL acom(nTrac), gcom(nTrac)
      _RL a_1t(nTrac)
      _RL step_tr(nPhyp,nTrait)
      _RL vr_tr(nPhyp,nTrait),sd_tr(nPhyp,nTrait)
#ifdef SPEAD_CONTINUOUS_COVARIANCE
      _RL cv_tr(nPhyp,nCov),cr_tr(nPhyp,nCov)
#endif
#ifdef DARWIN_GRAZING_SWITCH
      _RL coeff_KTW(nPlank)
#endif
#ifdef DARWIN_DEBUG
      INTEGER iG, jG, k
      _RL dT
#endif
      INTEGER myThid, myIter
      _RL myTime

C !INPUT/OUTPUT PARAMETERS: ============================================
C  gTr    :: accumulates computed tendencies
      _RL diags(darwin_nDiag)
      _RL gTr(nDarwin)

C !OUTPUT PARAMETERS: ==================================================
CEOP

#ifdef SPEAD_CONTINUOUS_TRAIT

c !LOCAL VARIABLES: ====================================================
      INTEGER j, l
      INTEGER jz, jp, jt, jt1, jt2, jt3, jn, jn1, jn2
      INTEGER jc, jc12, jc13, jc23
      INTEGER g
C  a_0  :: net growth rate of dominant ecotype of each species
C  a_1  :: first derivatives of net growth rate with respect to traits
C  a_2  :: second derivative
C  a_11 :: cross-derivatives
#ifdef SPEAD_CONTINUOUS_TRAIT
      _RL a_0(nPhyp)
      _RL a_d1(nPhyp,nTrait), a_d2(nPhyp,nTrait)
      _RL dmndt(nPhyp,nTrait), dvrdt(nPhyp,nTrait)
#ifdef SPEAD_CONTINUOUS_COVARIANCE
      _RL a_d11(nPhyp,nCov)
      _RL dcvdt(nPhyp,nCov)
#endif
#endif
      CHARACTER*(MAX_LEN_MBUF) msgBuf

C======================================================================
C derivatives with respect to trait

      DO jp = 1, nTrac
        IF ( X(plank(jp)) .GT. 0 _d 0 ) THEN
          a_1t(jp) = a_1t(jp) / X(plank(jp))
        ENDIF
      ENDDO

      DO jp = 1, nPlank
        IF ( X(jp) .GT. 0 _d 0 ) THEN
          gcom(jp) = gcom(jp) / X(jp)
          acom(jp) = acom(jp) / X(jp)
        ENDIF
      ENDDO

      DO jp = 1, nPhyp

        a_0(jp) = a_1t(iplank(jp))

        jn = iplank(jp) + 1
        DO jt = 1, num_trait(jp)
          a_d1(jp,jt) = (a_1t(jn+1) - a_1t(jn)) / (2*step_tr(jp,jt))
          a_d2(jp,jt) = (a_1t(jn+1) + a_1t(jn) - 2*a_1t(iplank(jp)))
     &                  / (step_tr(jp,jt)**2)
          a_d2(jp,jt) = MIN(0., a_d2(jp,jt))
          jn = jn + 2
        ENDDO
        DO jt = num_trait(jp)+1, nTrait
          a_d1(jp,jt) = 0.0 _d 0
          a_d2(jp,jt) = 0.0 _d 0
        ENDDO
#ifdef SPEAD_CONTINUOUS_COVARIANCE
        DO jt1 = 1, num_trait(jp)-1
         DO jt2 = jt1+1, num_trait(jp)
          jc = NINT((jt1-1)*(num_trait(jp)-1-jt1/2.0)) + jt2 - 1
          jn1 = iplank(jp) + 2*jt1
          jn2 = iplank(jp) + 2*jt2
          a_d11(jp,jc) = ( a_1t(jn) + a_1t(iplank(jp))
     &                   - a_1t(jn1) - a_1t(jn2) )
     &                   / (step_tr(jp,jt1)*step_tr(jp,jt2))
          jn = jn + 1
         ENDDO
        ENDDO
        DO jc = num_cov(jp)+1, nCov
          a_d11(jp,jc) = 0.0 _d 0
        ENDDO
#endif

      END DO

C======================================================================
C tendencies of trait moments

C tendencies directly computed in the biogeochemical functions
C (this can be changed if necessary
      DO jp = 1, nPlank
        gTr(ic+iplank(jp)-1) = gTr(ic+iplank(jp)-1) + acom(jp)*X(jp)
      ENDDO

      DO jp = 1, nPhyp

       DO jt = 1, num_trait(jp)
C       dmndt(jp,jt) = dmndt(jp,jt) + a_d1(jp,jt)*vr_tr(jp,jt)
C       dvrdt(jp,jt) = dvrdt(jp,jt) + a_d2(jp,jt)*(vr_tr(jp,jt)**2)
C    &                + 2*gcom(jp)*numut_tr(jp,jt)     
        dmndt(jp,jt) = a_d1(jp,jt)*vr_tr(jp,jt)
        dvrdt(jp,jt) = a_d2(jp,jt)*(vr_tr(jp,jt)**2) 
     &                 + 2*gcom(jp)*numut_tr(jp,jt)     
       END DO
C       print*,'-------------------------------------------------------'
C       print*, 'dvrdt step 1:', dvrdt(1,:)

#ifdef SPEAD_CONTINUOUS_COVARIANCE
C Perhaps the matrix form could help us in finding the right equation
C Terms occurring when num_trait(jp) > 2
       DO jt1 = 1, num_trait(jp)-1
        DO jt2 = jt1+1, num_trait(jp)
         jc = NINT((jt1-1)*(num_trait(jp)-1-jt1/2.0)) + jt2 - 1
         dmndt(jp,jt1) = dmndt(jp,jt1) + a_d1(jp,jt2)*cv_tr(jp,jc)
         dmndt(jp,jt2) = dmndt(jp,jt2) + a_d1(jp,jt1)*cv_tr(jp,jc)
         dvrdt(jp,jt1) = dvrdt(jp,jt1)
     &                 + 2*a_d11(jp,jc)*vr_tr(jp,jt1)*cv_tr(jp,jc)
     &                 + a_d2(jp,jt2)*(cv_tr(jp,jc)**2)
         dvrdt(jp,jt2) = dvrdt(jp,jt2)
     &                 + 2*a_d11(jp,jc)*vr_tr(jp,jt2)*cv_tr(jp,jc)
     &                 + a_d2(jp,jt1)*(cv_tr(jp,jc)**2)
C        dcvdt(jp,jc)  = dcvdt(jp,jc)
C    &                 + a_d2(jp,jt1)*vr_tr(jp,jt1)*cv_tr(jp,jc)
C    &                 + a_d2(jp,jt2)*vr_tr(jp,jt2)*cv_tr(jp,jc)
C    &   + a_d11(jp,jc)*(vr_tr(jp,jt1)*vr_tr(jp,jt2)+(cv_tr(jp,jc)**2))
         dcvdt(jp,jc)  = a_d2(jp,jt1)*vr_tr(jp,jt1)*cv_tr(jp,jc)
     &                 + a_d2(jp,jt2)*vr_tr(jp,jt2)*cv_tr(jp,jc)
     &   + a_d11(jp,jc)*(vr_tr(jp,jt1)*vr_tr(jp,jt2)+(cv_tr(jp,jc)**2))
        ENDDO
       ENDDO
C       print*, 'dvrdt step 2:', dvrdt(1,:)
C       print*, 'dcvdt step 2:', dcvdt(1,:)

       DO jt1 = 1, num_trait(jp)-1
        DO jt2 = jt1+1, num_trait(jp)
         jc = NINT((jt1-1)*(num_trait(jp)-1-jt1/2.0)) + jt2 - 1
C Terms occurring when num_trait(jp) > 3
         DO jt3 = jt2+1, num_trait(jp)
              jc13 = jc + jt3 - jt2
              jc23 = jc13
     &             + NINT((jt2-jt1)*(num_trait(jp)-(jt1+jt2+1)/2.0))
C             print*, 'jt1', jt1, 'jt2', jt2, 'jt3', jt3, 'jc', jc,
C    &                'jc13', jc13, 'jc23', jc23
              dvrdt(jp,jt1) = dvrdt(jp,jt1)
     &                   + 2*a_d11(jp,jc23)*cv_tr(jp,jc)*cv_tr(jp,jc13)
              dvrdt(jp,jt2) = dvrdt(jp,jt2)
     &                   + 2*a_d11(jp,jc13)*cv_tr(jp,jc)*cv_tr(jp,jc23)
              dvrdt(jp,jt3) = dvrdt(jp,jt3)
     &                   + 2*a_d11(jp,jc)*cv_tr(jp,jc13)*cv_tr(jp,jc23)
              dcvdt(jp,jc)  = dcvdt(jp,jc)
     &                + a_d11(jp,jc13) * ( vr_tr(jp,jt1)*cv_tr(jp,jc23)
     &                                  + cv_tr(jp,jc)*cv_tr(jp,jc13) )
     &                + a_d11(jp,jc23) * ( vr_tr(jp,jt2)*cv_tr(jp,jc13)
     &                                  + cv_tr(jp,jc)*cv_tr(jp,jc23) )
     &                + a_d2(jp,jt3) * cv_tr(jp,jc13) * cv_tr(jp,jc23)

              dcvdt(jp,jc13)  = dcvdt(jp,jc13)
     &               + a_d11(jp,jc)   * ( vr_tr(jp,jt1)*cv_tr(jp,jc23)
     &                                 + cv_tr(jp,jc13)*cv_tr(jp,jc) )
     &               + a_d11(jp,jc23) * ( vr_tr(jp,jt3)*cv_tr(jp,jc)
     &                               + cv_tr(jp,jc13)*cv_tr(jp,jc23) )
     &               + a_d2(jp,jt2) * cv_tr(jp,jc) * cv_tr(jp,jc23)

              dcvdt(jp,jc23)  = dcvdt(jp,jc23)
     &               + a_d11(jp,jc13) * ( vr_tr(jp,jt3)*cv_tr(jp,jc)
     &                               + cv_tr(jp,jc13)*cv_tr(jp,jc23) )
     &               + a_d11(jp,jc) * ( vr_tr(jp,jt2)*cv_tr(jp,jc13)
     &                                 + cv_tr(jp,jc)*cv_tr(jp,jc23) )
     &               + a_d2(jp,jt1) * cv_tr(jp,jc13) * cv_tr(jp,jc)

C Terms occurring when num_trait(jp) > 4 (unchecked)
C             DO jt4 = jt3+1, num_trait(j)
C               jc14 = jc   + jt4 - jt1
C               jc24 = jc23 + jt4 - jt3
C               jc34 = jc24
C    &               + nint((jt3-jt2)*(num_trait(j)-(jt2+jt3+1)/2.))
C               dcvdt(j,jc) = dcvdt(j,jc)
C    &               + a_d11(j,jc34) * ( cv_tr(j,jc13)*cv_tr(j,jc24)
C    &                                 + cv_tr(j,jc14)*cv_tr(j,jc23) )
C               dcvdt(j,jc13) = dcvdt(j,jc13)
C    &               + a_d11(j,jc24) * ( cv_tr(j,jc)  *cv_tr(j,jc34)
C    &                                 + cv_tr(j,jc14)*cv_tr(j,jc23) )
C               dcvdt(j,jc14) = dcvdt(j,jc14)
C    &               + a_d11(j,jc23) * ( cv_tr(j,jc)  *cv_tr(j,jc34)
C    &                                 + cv_tr(j,jc13)*cv_tr(j,jc24) )
C               dcvdt(j,jc23) = dcvdt(j,jc14)
C    &               + a_d11(j,jc14) * ( cv_tr(j,jc)  *cv_tr(j,jc34)
C    &                                 + cv_tr(j,jc13)*cv_tr(j,jc24) )
C               dcvdt(j,jc24) = dcvdt(j,jc24)
C    &               + a_d11(j,jc13) * ( cv_tr(j,jc)  *cv_tr(j,jc34)
C    &                                 + cv_tr(j,jc14)*cv_tr(j,jc23) )
C               dcvdt(j,jc34) = dcvdt(j,jc34)
C    &               + a_d11(j,jc)   * ( cv_tr(j,jc13)*cv_tr(j,jc24)
C    &                                 + cv_tr(j,jc14)*cv_tr(j,jc23) )
C             ENDDO
            ENDDO
          ENDDO
        ENDDO
C       print*, 'dvrdt step 3:', dvrdt(1,:)
C       print*, 'dcvdt step 3:', dcvdt(1,:)
#endif
      
C       print*, 'coeff_KTW:', coeff_KTW  
C       print*, 'vr_tr:', vr_tr(jp,:)
C       print*, 'sd_tr:', sd_tr(jp,:)
C       print*, 'cv_tr:', cv_tr(jp,:)
C       print*, 'cr_tr:', cr_tr(jp,:)


C Parameterization of Kill The Winner (KTW) preferential grazing 
#ifdef DARWIN_GRAZING_SWITCH 
#ifdef DARWIN_GRAZING_SWITCH_1trait
        dvrdt(jp,1) = dvrdt(jp,1) + coeff_KTW(jp)*vr_tr(jp,1)
        DO jt = 2, num_trait(jp)
          dvrdt(jp,jt) = dvrdt(jp,jt) 
     &            + coeff_KTW(jp)*vr_tr(jp,jt)*(cr_tr(jp,jt-1)**2)
        ENDDO
#ifdef SPEAD_CONTINUOUS_COVARIANCE
        DO jt1 = 1, num_trait(jp)-1
          DO jt2 = jt1+1, num_trait(jp)
            jc = NINT((jt1-1)*(num_trait(jp)-jt1/2.0)) + jt2 - jt1
            IF (jt1 .EQ. 1) THEN
              dcvdt(jp,jc) = dcvdt(jp,jc) + coeff_KTW(jp)*cv_tr(jp,jc)
            ELSE
             dcvdt(jp,jc) = dcvdt(jp,jc) + coeff_KTW(jp)*
     &       sd_tr(jp,jt1)*sd_tr(jp,jt2)*cr_tr(jp,jt1-1)*cr_tr(jp,jt2-1)
            ENDIF
          ENDDO
        ENDDO
#endif
#else /* NOT darwin_grazing_switch_1trait */
        DO jt = 1, num_trait(jp)
          dvrdt(jp,jt) = dvrdt(jp,jt) + coeff_KTW(jp)*vr_tr(jp,jt)
        ENDDO
        DO jc = 1, num_cov(jp)
          dcvdt(jp,jc) = dcvdt(jp,jc) + coeff_KTW(jp)*cv_tr(jp,jc)
        ENDDO
#endif /* darwin_grazing_switch_1trait */
#endif /* darwin_grazing_switch */

C Trends of conserved tracers
        jn = ic + iplank(jp)
        DO jt = 1, num_trait(jp)
          gTr(jn) = gTr(jn) + dmndt(jp,jt)
          jn = jn + 1
        ENDDO
        DO jt = 1, num_trait(jp)
          gTr(jn) = gTr(jn) + dvrdt(jp,jt)
          jn = jn + 1
        ENDDO
        DO jt1 = 1,num_trait(jp)-1
          DO jt2 = jt1+1,num_trait(jp)
            jc = NINT((jt1-1)*(num_trait(jp)-jt1/2.0)) + jt2 - jt1
            gTr(jn) = gTr(jn) + dcvdt(jp,jc)
            jn = jn + 1
          ENDDO
        ENDDO

      ENDDO
#ifdef DARWIN_DEBUG
C     IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG) THEN
        print*, 'dX:', acom(1)*X(1)
        print*, 'dmndt', dmndt(1,:)
        print*, 'dvrdt', dvrdt(1,:)
        print*, 'dcvdt:', dcvdt(1,:)
CC      print*, 'Ptr:', Ptr(ic:ic+nTrac-2)
        print*, 'gTr:', gTr(ic:ic+nTrac-2)
        print*,'-------------------------------------------------------'
C     ENDIF
#endif

C     STOP 'error'

C Check that no NaN is present (Le Gland, 13/05/2021)
      DO jp=1,nPhyp
        IF (acom(jp) /= acom(jp) .AND. myIter > 1) THEN
#ifdef DARWIN_DEBUG
          WRITE(msgBuf,'(2A,2I6,A,3I6,A,2F20.6,4(A,3F20.6))')
     &    'DARWIN_PLANKTON: ',
     &    'acom(jp) has NaN value', jp, myIter, ' Location:', iG,jG,k,
#else
          WRITE(msgBuf,'(2A,2I6,A,2F20.6,4(A,3F20.6))')
     &    'DARWIN_PLANKTON: ',
     &    'acom(jp) has NaN value', j, myIter,
#endif
     &    '  Concentrations:', X(jp), X(jp+1),
CC   &    '  Mean trait values:', mn_tr(j,1), mn_tr(j,2), mn_tr(j,3),
     &    '  Reference trait values:', ref_tr(jp,1), ref_tr(jp,2),
     &    ref_tr(j,3),
     &    '  Trait variances:', vr_tr(jp,1), vr_tr(jp,2), vr_tr(jp,3),
     &    '  Trait covariances:', cv_tr(jp,1), cv_tr(jp,2), cv_tr(jp,3)
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R DARWIN_PLANKTON'
        ENDIF
      ENDDO

C     STOP 'error'

C Moment derivatives
#ifdef SPEAD_DIAG_MD
      DO jp = 1,nPhyp
        jc = 1
        DO jt1 = 1,nTrait
          diags(iDmnt+(jp-1)*nTrait+jt1-1) = dmndt(jp,jt1)
          diags(iDvrt+(jp-1)*nTrait+jt1-1) = dvrdt(jp,jt1)
          diags(iDvrs+(jp-1)*nTrait+jt1-1) = 2*gcom(jp)*numut_tr(jp,jt1)
#ifdef DARWIN_GRAZING_SWITCH
#ifdef DARWIN_GRAZING_SWITCH_1trait
          IF(jt1 .EQ. 1) THEN
            diags(iDvrs+(jp-1)*nTrait+jt1-1) =
     &      diags(iDvrs+(jp-1)*nTrait+jt1-1) + coeff_KTW(jp)*vr_tr(jp,1)
          ELSE
            diags(iDvrs+(jp-1)*nTrait+jt1-1) = diags(iDvrs+(jp-1)*nTrait
     &      +jt1-1) + coeff_KTW(jp)*vr_tr(jp,jt1)*(cr_tr(jp,jt1-1)**2)
          ENDIF
#else
     &                                     + coeff_KTW(jp)*vr_tr(jp,jt1)
#endif
#endif
          DO jt2 = jt1+1,nTrait
            diags(iDcvt+(jp-1)*nCov+jc-1) = dcvdt(jp,jc)
#ifdef DARWIN_GRAZING_SWITCH
#ifdef DARWIN_GRAZING_SWITCH_1trait
          IF(jt1 .EQ. 1) THEN
            diags(iDcvs+(jp-1)*nCov+jc-1) = coeff_KTW(jp)*cv_tr(jp,jc)
          ELSE
            diags(iDcvs+(jp-1)*nCov+jc-1) = coeff_KTW(jp)*sd_tr(jp,jt1)*
     &      sd_tr(jp,jt2)*cr_tr(jp,jt1-1)*cr_tr(jp,jt2-1)
          ENDIF
#else
            diags(iDcvs+(jp-1)*nCov+jc-1) = coeff_KTW(jp)*cv_tr(jp,jc)
#endif
#else
            diags(iDcvs+(jp-1)*nCov+jc-1) = 0
#endif
            jc = jc + 1
          ENDDO
        ENDDO
      ENDDO
#endif

C     STOP 'error'

#ifdef DARWIN_DEBUG
C     IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG) THEN
        print*,'myTime',myTime,'myIter',myIter
        print*,'gcom:', gcom(1), 'acom:', acom(1),'a_0:',a_0(1)
        print*,'a_1t:', a_1t
        print*,'a_d1:', a_d1(1,:)
        print*,'a_d2:', a_d2(1,:)
        print*,'a_d11:', a_d11(1,:)
        print*,'alternative acom:', a_0(1)+(1./2.)*vr_tr(1,1)*a_d2(1,1)
     &         + (1./2.)*vr_tr(1,2)*a_d2(1,2) + cv_tr(1,1)*a_d11(1,1)
     &         + (1./2.)*vr_tr(1,3)*a_d2(1,3) + cv_tr(1,2)*a_d11(1,2)
     &         + cv_tr(1,3)*a_d11(1,3)
#ifdef DARWIN_GRAZING_SWITCH
        print*,'a_KTW', a_KTW, 'coeff_KTW', coeff_KTW
#endif
C     ENDIF
#endif
C       STOP 'error'

#endif /* spead_continuous_trait */

      RETURN
      END SUBROUTINE
