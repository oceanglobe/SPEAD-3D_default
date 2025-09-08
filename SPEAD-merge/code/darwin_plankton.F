#include "DARWIN_OPTIONS.h"

CBOP
C !ROUTINE: DARWIN_PLANKTON
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_PLANKTON(
     I     Ptr,
     O     gTr,
     O     chlout, diags,
     I     PAR, temp,
#ifdef DARWIN_DEBUG
     I     iG, jG, k, dT,
#endif
     I     myTime,myIter,myThid)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_INDICES.h"
#include "DARWIN_TRAITS.h"
#include "DARWIN_DIAGS.h"
#ifdef DARWIN_DEBUG
#include "DARWIN_PARAMS.h"
#include "DARWIN_TRAITPARAMS.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  Ptr    :: darwin model tracers
C  PAR    :: PAR in uEin/s/m2
C         :: either non-spectral (tlam=1) or waveband total
C  myTime :: current time
C  myIter :: current iteration number
C  myThid :: thread number
      _RL Ptr(nDarwin)
      _RL PAR(nlam)
      _RL temp
#ifdef DARWIN_DEBUG
      INTEGER iG, jG, k
      _RL dT
#endif
      INTEGER myThid, myIter
      _RL myTime

C !INPUT/OUTPUT PARAMETERS: ============================================

C !OUTPUT PARAMETERS: ==================================================
C  gTr    :: accumulates computed tendencies
C  chlout :: computed acclimated chlorophyll if not dynamic
C  diags  :: Darwin diagnostics
      _RL gTr(nDarwin)
      _RL chlout(nPhyp)
      _RL diags(darwin_nDiag)
CEOP

#ifdef ALLOW_DARWIN

c !LOCAL VARIABLES: ====================================================
      INTEGER j, l
      INTEGER jz, jp
C Le Gland (18/12/2020)
      INTEGER g
C Indices for traits and covariances (Le Gland, 29/03/2021)
      INTEGER jt, jt1, jt2, jt3, jt4
      INTEGER jc, jc13, jc23, jc14, jc24, jc34
C Error messages for debugging (Le Gland, 16/03/2021)
      CHARACTER*(MAX_LEN_MBUF) msgBuf

      _RL photoTempFunc(nTrac), grazTempFunc(nTrac)
      _RL reminTempFunc, uptakeTempFunc, mortTempFunc, mort2TempFunc

#ifdef SPEAD_DIC
      _RL  DIC,  PIC
#endif
      _RL  NH4,  NO2,  NO3,  PO4,  FeT
      _RL  DOC,  POC
#ifdef SPEAD_SILICA
      _RL  SiO2
#endif
#ifdef SPEAD_STOICHIOMETRY_OM
      _RL  DON,  DOP,  DOFe,  PON,  POP,  POFe
#ifdef SPEAD_SILICA
      _RL  POSi
#endif
#endif
#ifdef DARWIN_ALLOW_CARBON
      _RL  O2,  ALK
#endif

      _RL X(nPlank), acom(nPlank), gcom(nPlank), dXdt(nPlank)
      _RL step_tr(nPhyp,nTrait)
      _RL wght(nTrac)
C I removed the quotas for legibility 
#ifdef SPEAD_CONTINUOUS_TRAIT
C Means, variances and standard deviations
      _RL mn_tr(nPhyp,nTrait), vr_tr(nPhyp,nTrait)
      _RL sd_tr(nPhyp,nTrait)
C Growth and net growth rates (Le Gland, 10/03/2021)
      _RL a_1t(nTrac)
C     _RL a_0(nPhyp), a_d1(nPhyp,nTrait), a_d2(nPhyp,nTrait)
#ifdef SPEAD_CONTINUOUS_COVARIANCE
      _RL cv_tr(nPhyp,nCov)
      _RL cr_tr(nPhyp,nCov)
C     _RL a_d11(nPhyp,nCov)
#endif
#endif /* spead_continuous_trait */

#ifdef DARWIN_GRAZING_SWITCH
      _RL coeff_KTW(nPlank)
#endif

C Failed attempt to cap competitive exclusion in the continuous model
C     _RL mcex

C-----------------------------------------------------------------------
C Step 1: define non-living tracers
C-----------------------------------------------------------------------
C==== make all bio fields non-negative and compute quotas ==============

C     print*,'DARWIN_PLANKTON: step 1 starting:'

#ifdef SPEAD_DIC
      DIC  = MAX(0., Ptr(iDIC))
      PIC  = MAX(0., Ptr(iPIC))
#endif
      NH4  = MAX(0., Ptr(iNH4))
      NO2  = MAX(0., Ptr(iNO2))
      NO3  = MAX(0., Ptr(iNO3))
      PO4  = MAX(0., Ptr(iPO4))
#ifdef SPEAD_SILICA
      SiO2 = MAX(0., Ptr(iSiO2))
#endif
      FeT  = MAX(0., Ptr(iFeT))

      DOC  = MAX(0., Ptr(iDOC))
      POC  = MAX(0., Ptr(iPOC))
#ifdef SPEAD_STOICHIOMETRY_OM
      DON  = MAX(0., Ptr(iDON))
      DOP  = MAX(0., Ptr(iDOP))
      DOFe = MAX(0., Ptr(iDOFe))
      PON  = MAX(0., Ptr(iPON))
      POP  = MAX(0., Ptr(iPOP))
      POFe = MAX(0., Ptr(iPOFe))
#ifdef SPEAD_SILICA
      POSi = MAX(0., Ptr(iPOSi))
#endif
#endif

#ifdef DARWIN_ALLOW_CARBON
      O2   = MAX(0., Ptr(iO2))
      ALK  = MAX(0., Ptr(iALK))
#endif
#ifdef DARWIN_ALLOW_CDOM
      CDOM = MAX(0., Ptr(iCDOM))
#endif

#ifdef DARWIN_DEBUG
CC    IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG) THEN
CC       print*,'Concentrations',myiter,k,j,
        print*,'nTrac:  ',  nTrac,
     &         'nPhoto:  ', nPhoto,
     &         'nPlank:  ', nPlank,
     &         'nPhyp:  ',  nPhyp
        print*,'Concentrations:  ',
     &     'NH4', iNH4, NH4, Ptr(iNH4),
     &     'NO2', iNO2, NO2, Ptr(iNO2),
     &     'NO3', iNO3, NO3, Ptr(iNO3),
     &     'PO4', iPO4, PO4, Ptr(iPO4),
     &     'FeT', iFeT, FeT, Ptr(iFeT),
     &     'DOC', iDOC, DOC, Ptr(iDOC),
     &     'POC', iPOC, POC, Ptr(iPOC),
     &     'DIN', NH4+NO2+NO3,
     &     'Ntot', NH4+NO2+NO3+(DOC+POC+SUM(X))*R_NC(1),
#ifdef DARWIN_ALLOW_CARBON
     &     'O2' , iO2 , O2 , Ptr(iO2),
#endif
     &     ''
CC    ENDIF
#endif

C-----------------------------------------------------------------------
C-------------------------- End Step 1 ---------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C Step 2: define plankton biomass (discrete and continuous)
C         and trait moments (continuous only)
C-----------------------------------------------------------------------

      DO jp = 1, nPlank
       X(jp) = MAX(0., Ptr(ic+iplank(jp)-1))
#ifdef SPEAD_CONTINUOUS_TRAIT
       IF (jp .LE. nPhyp) THEN
        DO jt = 1, num_trait(jp)
         mn_tr(jp,jt) = Ptr(ic+iplank(jp)-1+jt)
         vr_tr(jp,jt) = Ptr(ic+iplank(jp)-1+num_trait(jp)+jt)
C        IF(vr_tr(jp,jt) /= vr_tr(jp,jt)) THEN 
C         vr_tr(jp,jt) = 0 _d 0
C        ENDIF
         sd_tr(jp,jt) = sqrt(vr_tr(jp,jt))
        END DO
#ifdef SPEAD_CONTINUOUS_COVARIANCE
        DO jt1 = 1, num_trait(jp) - 1
         DO jt2 = jt1+1, num_trait(jp)
          jc = NINT((jt1-1)*(num_trait(jp)-1-jt1/2.0)) + jt2 - 1
          cv_tr(jp,jc) =  Ptr(ic+iplank(jp)-1+2*num_trait(jp)+jc)
C         IF(cv_tr(jp,jt) /= cv_tr(jp,jt)) THEN 
C          cv_tr(jp,jt) = 0 _d 0
C         ENDIF
          IF (sd_tr(jp,jt1) > 0 .AND. sd_tr(jp,jt2) > 0) THEN
           cr_tr(jp,jc) = cv_tr(jp,jc) / (sd_tr(jp,jt1)*sd_tr(jp,jt2))
          ELSE
           cr_tr(jp,jc) = 0.0 _d 0
          ENDIF
C         IF(cr_tr(jp,jt) /= cr_tr(jp,jt)) THEN 
C           cr_tr(jp,jt) = 0 _d 0
C         ELSE
C           cr_tr(jp,jt) = MIN(0.99,MAX(-0.99,cr_tr(jp,jt)))
C         ENDIF
         END DO
        END DO
#endif
       ENDIF
#endif /* spead_continuous_trait */

      ENDDO

#ifdef DARWIN_DEBUG
C     IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG) THEN
        print*,'Plankton concentrations', X(:)
#ifdef SPEAD_CONTINUOUS_TRAIT
        print*,'Mean traits:', mn_tr(1,:)
        print*,'Trait variances:', vr_tr(1,:)
        print*,'Trait standard deviations:', sd_tr(1,:)
#ifdef SPEAD_CONTINUOUS_COVARIANCE
        print*,'Trait covariances:', cv_tr(1,:)
        print*,'Trait correlations:', cr_tr(1,:)
#endif
#endif /* spead_continuous_trait */
C     ENDIF
#endif /* DARWIN_DEBUG */

C-----------------------------------------------------------------------
C-------------------------- End Step 2 ---------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C Step 3: Compute uptake, growth, death and their derivatives with
C respect to trait
C-----------------------------------------------------------------------

C Check that no NaN is present
      DO jp=1,nPlank
        IF (X(jp) /= X(jp)) THEN
#ifdef DARWIN_DEBUG
C         WRITE(msgBuf,'(2A,2I6,A,3I6)') 'DARWIN_PLANKTON: ',
C    &    'X(jp) has NaN value', j, myIter, ' Location:', iG, jG, k
#else
          WRITE(msgBuf,'(2A,2I6)') 'DARWIN_PLANKTON: ',
     &    'X(jp) has NaN value', j, myIter
#endif
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R DARWIN_PLANKTON'
        ENDIF
      ENDDO

      DO jp = 1, nDarwin
        gTr(jp) = 0.0 _d 0
      ENDDO
      DO jp = 1, darwin_nDiag
        diags(jp) = 0.0 _d 0
      ENDDO
      DO jp = 1, nPlank
        acom(jp) = 0.0 _d 0
        gcom(jp) = 0.0 _d 0
      ENDDO
#ifdef SPEAD_CONTINUOUS_TRAIT
      DO jp = 1, nTrac
        a_1t(jp) = 0.0 _d 0
      ENDDO
#endif

C Set the traits and weight of each plankton species
#ifdef SPEAD_CONTINUOUS_TRAIT
      CALL DARWIN_SPEAD_WEIGHT(
     I           myTime,myIter,myThid,
     I           mn_tr,sd_tr,
#ifdef SPEAD_CONTINUOUS_COVARIANCE
     I           cv_tr,
#endif
#ifdef DARWIN_DEBUG
     I     iG, jG, k, dT,
#endif
     O           step_tr,
     O           wght)
#else /* NOT spead_continuous_trait */
      DO jp = 1, nTrac
        wght(jp) = 1.0 _d 0
      ENDDO
#endif /* spead_continuous_trait */

C temperature-dependent multiplicators
      CALL DARWIN_TEMPFUNC(
     I         temp,
     O         photoTempFunc, grazTempFunc, reminTempFunc,
     O         mortTempFunc, mort2TempFunc, uptakeTempFunc,
     I         myThid)

C Photosynthetic growth and nutrient uptake
      CALL DARWIN_SPEAD_PHYTO(
     I     PAR,photoTempFunc,uptakeTempFunc,
     I     NH4,NO2,NO3,PO4,FeT,
#ifdef SPEAD_DIC
     I     DIC,
#endif
#ifdef SPEAD_SILICA
     I     SiO2,
#endif
     O     chlout,
     U     diags,
     U     gTr,
#ifdef DARWIN_DEBUG
     I     iG, jG, k, dT,
#endif
     I     X,wght,
     U     acom, gcom,
#ifdef SPEAD_CONTINUOUS_TRAIT
     U     a_1t,
#endif
     I     myTime,myIter,myThid)

C Bacterial growth and remineralization of organic matter
C Place to add biogeochemical processes from PISCES or Eco-3M ?
C (Or perhaps in another function where X is not an input)
      CALL DARWIN_SPEAD_REMIN(
     I     PAR,reminTempFunc,
     I     NH4,NO2,NO3,
#ifdef SPEAD_DIC
     I     PIC,
#endif
     I     DOC,POC,
#ifdef SPEAD_STOICHIOMETRY_OM
     I     DON,DOP,DOFe,PON,POP,POFe,
#ifdef SPEAD_SILICA
     I     POSi,
#endif
#endif /* spead_stoichiometry_om */
#ifdef DARWIN_ALLOW_CARBON
     I     O2,
#endif
#ifdef DARWIN_ALLOW_CDOM
     I     CDOM,
#endif
     U     diags,
     U     gTr,
#ifdef DARWIN_DEBUG
     I     iG, jG, k, dT,
#endif
     I     X,wght,
     U     acom,gcom,
#ifdef SPEAD_CONTINUOUS_TRAIT
     U     a_1t,
#endif
     I     myTime,myIter,myThid)

C Grazing, mortality and exudation
      CALL DARWIN_SPEAD_GRAZING(
     I     reminTempFunc,grazTempFunc,mortTempFunc,mort2TempFunc,
     U     diags,
     U     gTr,
#ifdef DARWIN_DEBUG
     I     iG, jG, k, dT,
#endif
     I     X,wght,
     U     acom,gcom,
#ifdef SPEAD_CONTINUOUS_TRAIT
     U     a_1t,
#ifdef DARWIN_GRAZING_SWITCH
     O     coeff_KTW,
#endif
#endif /* spead_continuous_trait */
     I     myTime,myIter,myThid)

C Derivatives with respect to traits and trait tendencies
#ifdef SPEAD_CONTINUOUS_TRAIT
      CALL DARWIN_SPEAD_DERIVATIVES(
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
#else /* NOT spead_continuous_trait */
      DO jp = 1, nPlank
        gTr(ic+jp-1) = gTr(ic+jp-1) + acom(jp)
      ENDDO
#endif

C     STOP 'error in DARWIN_PLANKTON'

C acom and a0 are redundant but I keep them for checking

#ifdef DARWIN_DEBUG
CC    IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG) THEN
        print*,'Concentration derivatives',
     &     gTr(iNO3:iDOC), gTr(iPOC),
#ifdef DARWIN_ALLOW_CARBON
     &     gTr(iO2:iALK),
#endif
     &     ''
CC    ENDIF
#endif

C-----------------------------------------------------------------------
C-------------------------- End Step 3 ---------------------------------
C-----------------------------------------------------------------------

#ifdef DARWIN_DEBUG
C     IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG) THEN
        print*, 'Ptr:', Ptr(ic:ic+nTrac-2)
#ifndef SPEAD_CONTINUOUS_TRAIT
        print*, 'sum:', SUM(Ptr(ic:ic+nTrac-2)), SUM(X)
#endif
        print*, 'gTr:', gTr(ic:ic+nTrac-2)
#ifndef SPEAD_CONTINUOUS_TRAIT
        print*, 'sum:', SUM(gTr(ic:ic+nTrac-2))
#endif
        print*,'-------------------------------------------------------'
C     ENDIF
#endif

C Check that no NaN is present (Le Gland, 13/05/2021)
C     STOP 'error in DARWIN_PLANKTON'
      DO jp =1,nPhyp
        IF (acom(jp) /= acom(jp) .AND. myIter > 1) THEN
#ifdef DARWIN_DEBUG
          WRITE(msgBuf,'(2A,2I6,A,3I6,A,2F20.6,4(A,3F20.6))')
     &    'DARWIN_PLANKTON: ',
     &    'acom(j) has NaN value', jp, myIter, ' Location:', iG,jG,k,
#else
          WRITE(msgBuf,'(2A,2I6,A,2F20.6,4(A,3F20.6))')
     &    'DARWIN_PLANKTON: ',
     &    'acom(j) has NaN value', jp, myIter,
#endif
#ifdef SPEAD_CONTINUOUS_TRAIT
     &    '  Mean trait values:', mn_tr(jp,1), mn_tr(jp,2), mn_tr(jp,3),
     &    '  Reference trait values:', ref_tr(jp,1), ref_tr(jp,2),
     &    ref_tr(j,3),
     &    '  Trait variances:', vr_tr(jp,1), vr_tr(jp,2), vr_tr(jp,3),
     &    '  Trait covariances:', cv_tr(jp,1), cv_tr(jp,2), cv_tr(jp,3),
#endif
     &    '  Concentrations:', X(jp)
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R DARWIN_PLANKTON'
        ENDIF
      ENDDO
C     IF (myTime .GT. 4000) STOP 'error in DARWIN_PLANKTON'

C discrete trait diffusion (continuous TD is managed in DERIVATIVES)
#ifdef SPEAD_TRAIT_DIFFUSION
#ifndef SPEAD_CONTINUOUS_TRAIT
      CALL DARWIN_SPEAD_TD_DISCRETE(
     I     X,
     I     gcom,
     U     gTr,
     I     myTime,myIter,myThid)
#endif
#endif /* spead_trait_diffusion */

C-----------------------------------------------------------------------
C-------------------------- End Step 4 ---------------------------------
C-----------------------------------------------------------------------

#endif /* ALLOW_DARWIN */


      RETURN
      END SUBROUTINE
