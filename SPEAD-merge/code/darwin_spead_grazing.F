#include "DARWIN_OPTIONS.h"

CBOP
C !ROUTINE: DARWIN_SPEAD_GRAZING
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_SPEAD_GRAZING(
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

C !DESCRIPTION: PLANKTON GRAZING, MORTALITY and EXUDATION

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
C  grazTempFunc,reminTempFunc, :: temperature response functions
C  mortTempFunc,mort2TempFunc  :: 
C  X      :: biomass (mmolC/m3) of phytoplankton species
C  wght   :: weight of each phytoplankton tracer (continuous-trait)
C  myTime :: current time
C  myIter :: current iteration number
C  myThid :: thread number
      _RL grazTempFunc(nTrac)
      _RL reminTempFunc, mortTempFunc, mort2TempFunc
#ifdef DARWIN_DEBUG
      INTEGER iG, jG, k
      _RL dT
#endif
      _RL X(nPlank)
      _RL wght(nTrac)
      INTEGER myThid, myIter
      _RL myTime

C !INPUT/OUTPUT PARAMETERS: ============================================
C  gTr    :: accumulates computed tendencies
C  diags  :: darwin diagnostics
C  acom   :: net growth rate (gross growth - death) of each species
C  gcom   :: gross growth rate of each species
C  a_1t   :: net growth rate of each tracer (continuous-trait)
      _RL gTr(nDarwin)
      _RL diags(darwin_nDiag)
      _RL gcom(nplank), acom(nplank)
#ifdef SPEAD_CONTINUOUS_TRAIT
      _RL a_1t(nTrac)
#endif

C !OUTPUT PARAMETERS: ==================================================
C coeff_KTW :: multiplicator of diversity under Kill-The-Winner grazing
#ifdef SPEAD_CONTINUOUS_TRAIT
#ifdef DARWIN_GRAZING_SWITCH
      _RL coeff_KTW(nplank)
#endif
#endif /* spead_continuous_trait */
CEOP

#ifdef ALLOW_DARWIN

c !LOCAL VARIABLES: ====================================================
      INTEGER jz, jp, jp2
C Error messages for debugging (Le Gland, 16/03/2021)
      CHARACTER*(MAX_LEN_MBUF) msgBuf

      _RL Qn(nTrac), Qp(nTrac), Qsi(nTrac), Qfe(nTrac)

C for grazing

      _RL regQc, regQn, regQp, regQfe
      _RL sumprey, sumpref, grazphy

      _RL preygraz(nTrac), preygrazexp(nTrac), predgrazc(nTrac)

      _RL Xi(nTrac,nTrac)

      _RL predexpc, predexpn, predexpp, predexpfe
      _RL graz2OC, graz2ON, graz2OP, graz2OFe
      _RL graz2POC, graz2PON, graz2POP, graz2POSi, graz2POFe
      _RL graz2PIC

      _RL expfrac

      _RL Xe
      _RL mortX, mortX2

      _RL exude_DOC
      _RL exude_DON
      _RL exude_DOP
      _RL exude_DOFe

      _RL exude_PIC
      _RL exude_POC
      _RL exude_PON
      _RL exude_POP
      _RL exude_POSi
      _RL exude_POFe

      _RL mort_c(nTrac)

      _RL tmp

#ifdef SPEAD_DIC
      _RL respir
#endif
      _RL respir_c

#ifdef DARWIN_ALLOW_CDOM
      _RL graz2CDOM, exude_CDOM
#endif

C Assign quotas when they are constant 
      DO jp = 1, nTrac
#ifndef DARWIN_ALLOW_NQUOTA
        Qn(jp) = R_NC(jp)
#endif
#ifndef DARWIN_ALLOW_PQUOTA
        Qp(jp) = R_PC(jp)
#endif
#ifndef DARWIN_ALLOW_NQUOTA
        Qsi(jp) = R_SiC(jp)
#endif
#ifndef DARWIN_ALLOW_NQUOTA
        Qfe(jp) = R_FeC(jp)
#endif
      ENDDO

#ifdef DARWIN_GRAZING_SWITCH
#ifdef SPEAD_CONTINUOUS_TRAIT
      DO jp = 1, nplank
        coeff_KTW(jp) = 0.0 _d 0
      ENDDO
#endif
      DO jz=1,nTrac
      IF (isPred(jz).NE.0) THEN
       DO jp=1,nTrac
        Xi(jp,jz) = palat(jp,jz)*X(plank(jp))
#ifndef SPEAD_CONTINUOUS_TRAIT
#ifdef DARWIN_GRAZING_SWITCH_1trait
        DO jp2=1,nTrac
         IF ((jp.NE.jp2) .AND. ( group(jp) .EQ.  group(jp2)) 
     &                   .AND. (igroup(jp) .EQ. igroup(jp2))) THEN
          Xi(jp,jz) = Xi(jp,jz) + palat(jp2,jz)*X(plank(jp2))
         ENDIF
        ENDDO
#endif
#endif /* spead_continuous_trait */
       ENDDO
      ENDIF
      ENDDO
#endif /* darwin_grazing_switch */

C==== grazing ==========================================================

      DO jp = 1, nTrac
       preygraz(jp)    = 0.0 _d 0
       preygrazexp(jp) = 0.0 _d 0
       predgrazc(jp)   = 0.0 _d 0
#ifdef DARWIN_ALLOW_NQUOTA
       predgrazn(jp)  = 0.0 _d 0
#endif
#ifdef DARWIN_ALLOW_PQUOTA
       predgrazp(jp)  = 0.0 _d 0
#endif
#ifdef DARWIN_ALLOW_FEQUOTA
       predgrazfe(jp) = 0.0 _d 0
#endif
      ENDDO
      graz2OC   = 0.0
      graz2POC  = 0.0
#ifdef SPEAD_DIC
      graz2PIC  = 0.0
#endif
#ifdef SPEAD_STOICHIOMETRY_OM
      graz2ON   = 0.0
      graz2OP   = 0.0
      graz2OFE  = 0.0
      graz2PON  = 0.0
      graz2POP  = 0.0
      graz2POFE = 0.0
#ifdef SPEAD_SILICA
      graz2POSI = 0.0
#endif
#endif /* spead_stoichiometry_om */

      regQn  = 1.0
      regQp  = 1.0
      regQfe = 1.0
      regQc  = 1.0

C=======================================================================
      DO jz = 1, nTrac
       IF (isPred(jz).NE.0) THEN

C       regulate grazing near full quota
        regQc = 1.0 _d 0
#ifdef DARWIN_ALLOW_NQUOTA
        regQn = MAX(0., MIN(1., (Qnmax(jz)-Qn(jz))/
     &                          (Qnmax(jz)-Qnmin(jz)) ))
        regQc = MIN(regQc, 1.0 _d 0 - regQn)
        regQn = regQn**hillnumGraz
#endif
#ifdef DARWIN_ALLOW_PQUOTA
        regQp = MAX(0., MIN(1., (Qpmax(jz)-Qp(jz))/
     &                          (Qpmax(jz)-Qpmin(jz)) ))
        regQc = MIN(regQc, 1.0 _d 0 - regQp)
        regQp = regQp**hillnumGraz
#endif
#ifdef DARWIN_ALLOW_FEQUOTA
        regQfe= MAX(0., MIN(1., (Qfemax(jz)-Qfe(jz))/
     &                          (Qfemax(jz)-Qfemin(jz)) ))
        regQc = MIN(regQc, 1.0 _d 0 - regQfe)
        regQfe=regQfe**hillnumGraz
#endif
        regQc = regQc**hillnumGraz

        sumprey = 0.0
        sumpref = 0.0
        DO jp = 1, nTrac
         IF (palat(jp,jz).NE.0 _d 0) THEN
          sumprey = sumprey + palat(jp,jz)*X(plank(jp))*wght(jp)
C#ifdef DARWIN_GRAZING_SWITCH
C          sumpref = sumpref + (palat(jp,jz)**a_KTW)
C     &              * (X(plank(jp))**a_KTW) * wght(jp)
C#else
C          sumpref = sumpref + palat(jp,jz)*X(plank(jp))*wght(jp)
C#endif
          sumpref = sumpref + palat(jp,jz)*X(plank(jp))*wght(jp)
#ifdef DARWIN_GRAZING_SWITCH
     &                       *Xi(jp,jz)**(a_KTW-1)
#endif
         ENDIF
        ENDDO
        sumprey = MAX(0.0, sumprey - phygrazmin)
        sumpref = MAX(phygrazmin, sumpref)
        tmp = grazemax(jz)*grazTempFunc(jz)**tempGraz(jz)*X(plank(jz))*
     &    (sumprey**hollexp/(sumprey**hollexp+kgrazesat(jz)**hollexp))*
     &    (1.0 - EXP(-inhib_graz*sumprey))**inhib_graz_exp

        predexpc  = 0.0 _d 0
        predexpn  = 0.0 _d 0
        predexpp  = 0.0 _d 0
        predexpfe = 0.0 _d 0
        DO jp = 1, nTrac
         IF (palat(jp,jz).NE.0 _d 0) THEN
C#ifdef DARWIN_GRAZING_SWITCH
C          grazphy = tmp*((palat(jp,jz)*X(plank(jp)))**a_KTW)/sumpref
C#else
C          grazphy = tmp*palat(jp,jz)*X(plank(jp))/sumpref
C#endif
          grazphy = tmp*palat(jp,jz)*X(plank(jp))/sumpref
#ifdef DARWIN_GRAZING_SWITCH
     &                 *Xi(jp,jz)**(a_KTW-1)
#endif

          expFrac = ExportFracPreyPred(jp,jz)

          preygraz(jp) = preygraz(jp) + grazphy*wght(jz)
          preygrazexp(jp) = preygrazexp(jp) + expFrac*grazphy*wght(jz)

          predgrazc(jz) = predgrazc(jz) + grazphy*asseff(jp,jz)*wght(jp)
          predexpc = predexpc + expFrac*grazphy*asseff(jp,jz)
     &                          *wght(jp)*wght(jz)
#ifdef DARWIN_ALLOW_NQUOTA
          predgrazn(jz) = predgrazn(jz) + grazphy*asseff(jp,jz)*
     &                                    regQn*Qn(jp)*wght(jp)
          predexpn = predexpn + expFrac*grazphy*asseff(jp,jz)*
     &                          regQn*Qn(jp)*wght(jp)*wght(jz)
#endif
#ifdef DARWIN_ALLOW_PQUOTA
          predgrazp(jz) = predgrazp(jz) + grazphy*asseff(jp,jz)*
     &                                 regQp*Qp(jp)*wght(jp)
          predexpp = predexpp + expFrac*grazphy*asseff(jp,jz)*
     &                          regQp*Qp(jp)*wght(jp)*wght(jz)
#endif
#ifdef DARWIN_ALLOW_FEQUOTA
          predgrazfe(jz) = predgrazfe(jz) + grazphy*asseff(jp,jz)*
     &                                      regQfe*Qfe(jp)*wght(jp)
          predexpfe = predexpfe + expFrac*grazphy*asseff(jp,jz)*
     &                            regQfe*Qfe(jp)*wght(jp)*wght(jz)
#endif
#ifdef DARWIN_DEBUG
C         IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG) THEN
            print*, 'predator', jz, 'prey', jp, 'palat', palat(jp,jz),
     &              'sumprey', sumprey, 'grazemax(jz)', grazemax(jz),
     &              'grazTempFunc(jz)', grazTempFunc(jz), 'tmp', tmp,
     &              'grazphy', grazphy, 'preygraz(jp)', preygraz(jp),
     &              'predgrazc(jz)', predgrazc(jz), 
CC   &              'Temperature', temp,
     &              'hollexp', hollexp, 'kgrazesat(jz)', kgrazesat(jz),
     &              'inhib_graz', inhib_graz, 'inhib_graz_exp',
     &              inhib_graz_exp, 'pred. conc.', X(plank(jz))
C         ENDIF
#endif
         ENDIF
        ENDDO

C organic-matter gain will be total preygraz - predgraz
        graz2OC   = graz2OC  - predgrazc(jz)*wght(jz)
        graz2POC  = graz2POC - predexpc

#ifdef SPEAD_STOICHIOMETRY_OM
#ifdef DARWIN_ALLOW_NQUOTA
        graz2ON   = graz2ON  - predgrazn(jz)*wght(jz)
        graz2PON  = graz2PON - predexpn
#else
        graz2ON   = graz2ON  - predgrazc(jz)*Qn(jz)*wght(jz)
        graz2PON  = graz2PON - predexpc     *Qn(jz)
#endif

#ifdef DARWIN_ALLOW_PQUOTA
        graz2OP   = graz2OP  - predgrazp(jz)*wght(jz)
        graz2POP  = graz2POP - predexpp
#else
        graz2OP   = graz2OP  - predgrazc(jz)*Qp(jz)*wght(jz)
        graz2POP  = graz2POP - predexpc     *Qp(jz)
#endif

#ifdef DARWIN_ALLOW_FEQUOTA
        graz2OFe  = graz2OFe  - predgrazfe(jz)*wght(jz)
        graz2POFe = graz2POFe - predexpfe
#else
        graz2OFe  = graz2OFe  - predgrazc(jz)*Qfe(jz)*wght(jz)
        graz2POFe = graz2POFe - predexpc     *Qfe(jz)
#endif
#endif /* spead_stoichiometry_om */

       ENDIF
C     end predator loop
      ENDDO

      DO jp = 1, nTrac
       IF (isPrey(jp).NE.0) THEN
         graz2OC   = graz2OC   + preygraz(jp)*wght(jp)
         graz2POC  = graz2POC  + preygrazexp(jp)*wght(jp)
#ifdef DARWIN_ALLOW_CARBON
         graz2PIC = graz2PIC + preygraz(jp)*R_PICPOC(jp)*wght(jp)
#endif
#ifdef SPEAD_STOICHIOMETRY_OM
         graz2ON   = graz2ON   + preygraz(jp)*Qn(jp)*wght(jp)
         graz2OP   = graz2OP   + preygraz(jp)*Qp(jp)*wght(jp)
         graz2OFe  = graz2OFe  + preygraz(jp)*Qfe(jp)*wght(jp)
         graz2PON  = graz2PON  + preygrazexp(jp)*Qn(jp)*wght(jp)
         graz2POP  = graz2POP  + preygrazexp(jp)*Qp(jp)*wght(jp)
         graz2POFe = graz2POFe + preygrazexp(jp)*Qfe(jp)*wght(jp)
#ifdef SPEAD_SILICA
         graz2POSi = graz2POSi + preygraz(jp)*Qsi(jp)*wght(jp)
#endif
#endif /* spead_stoichiometry_om */
       ENDIF
      ENDDO

C==== tendencies =======================================================

      gTr(iDOC )=gTr(iDOC ) + graz2OC  - graz2POC
      gTr(iPOC )=gTr(iPOC ) + graz2POC

#ifdef SPEAD_STOICHIOMETRY_OM
      gTr(iDON )=gTr(iDON ) + graz2ON  - graz2PON
      gTr(iDOP )=gTr(iDOP ) + graz2OP  - graz2POP
      gTr(iDOFe)=gTr(iDOFe) + graz2OFe - graz2POFe
      gTr(iPON )=gTr(iPON ) + graz2PON
      gTr(iPOP )=gTr(iPOP ) + graz2POP
      gTr(iPOFe)=gTr(iPOFe) + graz2POFe
#ifdef SPEAD_SILICA
      gTr(iPOSi)=gTr(iPOSi) + graz2POSi
#endif
#endif /* spead_stoichiometry_om */

#ifdef SPEAD_DIC
      gTr(iPIC )=gTr(iPIC ) + graz2PIC
#endif

#ifdef DARWIN_ALLOW_CDOM
      graz2CDOM = fracCDOM*(graz2OP - graz2POP)
      gTr(iCDOM)=gTr(iCDOM) + graz2CDOM
      gTr(iDOC )=gTr(iDOC )             - R_CP_CDOM*graz2CDOM
#ifdef SPEAD_STOICHIOMETRY_OM
      gTr(iDON )=gTr(iDON )             - R_NP_CDOM*graz2CDOM
      gTr(iDOP )=gTr(iDOP ) - graz2CDOM
      gTr(iDOFe)=gTr(iDOFe)             - R_FeP_CDOM*graz2CDOM
#endif
#endif /* darwin_allow_cdom */

      DO jp = 1, nTrac
      IF (isPrey(jp).NE.0) THEN
C      gTr(iplank(plank(jp)))= gTr(iplank(plank(jp))) 
C    &                        - preygraz(jp)*wght(jp)
       acom(plank(jp)) = acom(plank(jp)) - preygraz(jp)*wght(jp)
#ifdef SPEAD_CONTINUOUS_TRAIT
       a_1t(jp) = a_1t(jp) - preygraz(jp)
#ifdef DARWIN_GRAZING_SWITCH
C I Added the conditions on X>0 to avoid Nans (Boris 07/01/24)
        IF (X(plank(jp)) .GT. 0) THEN
           coeff_KTW(plank(jp)) = coeff_KTW(plank(jp)) + preygraz(jp)
     &             *(1 _d 0 - (1 _d 0)/a_KTW)*wght(jp)/X(plank(jp))
        ENDIF
        
#ifdef DARWIN_DEBUG
C       IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG) THEN
          print*, 'coeff_KTW', coeff_KTW(plank(jp)), 'preygraz',
     &            preygraz(jp), 'a_KTW', a_KTW, '1-1/a',
     &            (1 _d 0 - (1 _d 0)/a_KTW), 'wght', wght(jp),
     &            'X', X(plank(jp)), 'added coeff_KTW', preygraz(jp)*
     &            (1 _d 0 - (1 _d 0)/a_KTW)*wght(jp)/X(plank(jp))
C       ENDIF
#endif
#endif /* darwin_grazing_switch */
#endif /* spead_continuous_trait */
      ENDIF
      ENDDO

C tendencies with quotas would be here

C     DO jz = 1, nplank
      DO jz = 1, nTrac
      IF (isPred(jz).NE.0) THEN
C      gTr(iplank(plank(jp))) = gTr(iplank(plank(jp))) 
C    &                         + preygraz(jp)*wght(jp)
       acom(plank(jz)) = acom(plank(jz)) + predgrazc(jz)*wght(jz)
       gcom(plank(jz)) = gcom(plank(jz)) + predgrazc(jz)*wght(jz)
#ifdef SPEAD_CONTINUOUS_TRAIT
       a_1t(jz) = a_1t(jz) + predgrazc(jz)
#endif

C tendencies with quotas would be here

      ENDIF
      ENDDO

#ifdef DARWIN_DIAG_PERTYPE
      DO jp = 1, nTrac
        diags(iGRplank+plank(jp)-1) = diags(iGRplank+plank(jp)-1)
     &                              + preygraz(jp)*wght(jp)
      ENDDO
      DO jz = 1, nTrac
        diags(iGrGn+plank(jz)-1) = diags(iGrGn+plank(jz)-1)
     &                           + predgrazc(jz)*wght(jz)
      ENDDO
#endif

C==== mortality ========================================================
      exude_DOC  = 0.0 _d 0
      exude_POC  = 0.0 _d 0
#ifdef SPEAD_DIC
      exude_PIC  = 0.0 _d 0
      respir     = 0.0 _d 0
#endif
#ifdef SPEAD_STOICHIOMETRY_OM
      exude_DON  = 0.0 _d 0
      exude_DOP  = 0.0 _d 0
      exude_DOFe = 0.0 _d 0
      exude_PON  = 0.0 _d 0
      exude_POP  = 0.0 _d 0
      exude_POFe = 0.0 _d 0
#ifdef SPEAD_SILICA
      exude_POSi = 0.0 _d 0
#endif
#endif /* spead_stoichiometry_om */

      DO jp = 1, nTrac
        Xe = MAX(0 _d 0, X(plank(jp)) - Xmin(jp))
        mortX = mort(jp)*Xe*mortTempFunc**tempMort(jp)
        mortX2= mort2(jp)*Xe*Xe*mort2TempFunc**tempMort2(jp)

        mort_c(jp) = mortX + mortX2

        exude_DOC = exude_DOC
     &            + (1.-ExportFracMort(jp)) *mortX *wght(jp)
     &            + (1.-ExportFracMort2(jp))*mortX2*wght(jp)
        exude_POC = exude_POC
     &            + ExportFracMort(jp) *mortX *wght(jp)
     &            + ExportFracMort2(jp)*mortX2*wght(jp)

#ifdef SPEAD_STOICHIOMETRY_OM
        exude_DON = exude_DON
     &            + (1.-ExportFracMort(jp)) *mortX *Qn(jp)*wght(jp)
     &            + (1.-ExportFracMort2(jp))*mortX2*Qn(jp)*wght(jp)
        exude_PON = exude_PON
     &            + ExportFracMort(jp) *mortX *Qn(jp)*wght(jp)
     &            + ExportFracMort2(jp)*mortX2*Qn(jp)*wght(jp)

        exude_DOP = exude_DOP
     &            + (1.-ExportFracMort(jp)) *mortX *Qp(jp)*wght(jp)
     &            + (1.-ExportFracMort2(jp))*mortX2*Qp(jp)*wght(jp)
        exude_POP = exude_POP
     &            + ExportFracMort(jp) *mortX *Qp(jp)*wght(jp)
     &            + ExportFracMort2(jp)*mortX2*Qp(jp)*wght(jp)

        exude_DOFe = exude_DOFe
     &           + (1.-ExportFracMort(jp)) *mortX *Qfe(jp)*wght(jp)
     &           + (1.-ExportFracMort2(jp))*mortX2*Qfe(jp)*wght(jp)
        exude_POFe = exude_POFe
     &           + ExportFracMort(jp) *mortX *Qfe(jp)*wght(jp)
     &           + ExportFracMort2(jp)*mortX2*Qfe(jp)*wght(jp)

#ifdef SPEAD_SILICA
        exude_POSi = exude_POSi + mort_c(jp)*Qsi(jp)*wght(jp)
#endif
#endif /* spead_stoichiometry_om */

        respir_c = respRate(jp)*Xe*reminTempFunc
#ifdef SPEAD_DIC
        respir = respir + respir_c*wght(jp)
        exude_PIC = exude_PIC + mort_c(jp)*R_PICPOC(jp)*wght(jp)
#endif

C       gTr(iplank(plank(jp))) = gTr(iplank(plank(jp))) 
C    &                           - (mort_c(jp) + respir_c) * wght(jp)
        acom(plank(jp)) = acom(plank(jp))
     &                    - (mort_c(jp) + respir_c) * wght(jp)
#ifdef SPEAD_CONTINUOUS_TRAIT
        a_1t(jp) = a_1t(jp) - mort_c(jp) - respir_c
#endif

C tendencies with quotas would be here

#ifdef DARWIN_ALLOW_EXUDE
        exude_DOC = exude_DOC
     &  + (1.-ExportFracExude(jp))*kexcc(jp)*Xe*wght(jp)
        exude_POC = exude_POC
     &  +     ExportFracExude(jp) *kexcc(jp)*Xe*wght(jp)
#ifdef SPEAD_STOICHIOMETRY_OM
        exude_DON = exude_DON
     &  + (1.-ExportFracExude(jp))*kexcn(jp)*Xe*Qn(jp)*wght(jp)
        exude_PON = exude_PON
     &  +     ExportFracExude(jp) *kexcn(jp)*Xe*Qn(jp)*wght(jp)
        exude_DOP = exude_DOP
     &  + (1.-ExportFracExude(jp))*kexcp(jp)*Xe*Qp(jp)*wght(jp)
        exude_POP = exude_POP
     &  +     ExportFracExude(jp) *kexcp(jp)*Xe*Qp(jp)*wght(jp)
        exude_DOFe = exude_DOFe
     &  + (1.-ExportFracExude(jp))*kexcfe(jp)*Xe*Qfe(jp)*wght(jp)
        exude_POFe = exude_POFe
     &  +     ExportFracExude(jp) *kexcfe(jp)*Xe*Qfe(jp)*wght(jp)
#ifdef SPEAD_SILICA
        exude_POSi = exude_POSi
     &             + kexcsi(jp)*Xe*Qsi(jp)*wght(jp)
#endif
#endif /* spead_stoichiometry_om */

C       gTr(iplank(plank(jp))) = gTr(iplank(plank(jp))) 
C    &                           - kexcc(jp)*Xe*wght(jp)
        acom(plank(jp)) = acom(plank(jp)) - kexcc(jp)*Xe*wght(jp)
#ifdef SPEAD_CONTINUOUS_TRAIT
        a_1t(jp) = a_1t(jp) - kexcc(jp)*Xe
#endif

C tendencies with quotas would be here
#endif /* darwin_allow_exude */

#ifdef DARWIN_DEBUG
C       IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG) THEN
          print*, 'tracer', jp, 'species', plank(jp), 
#ifdef SPEAD_CONTINUOUS_TRAIT
     &            'a_1t', a_1t(jp),
#endif
     &            'acom', acom(plank(jp)), 'mort_c', mort_c(jp),
     &            'mortX', mortX, 'mortX2', mortX2, 'mort', mort(jp),
     &            'mort2', mort2(jp), 'mortTempFunc', mortTempFunc,
     &            'mort2TempFunc', mort2TempFunc, 'respir_c',
     &            respir_c, 'kexcc', kexcc(jp), 'Xe', Xe,
     &            'respRate', respRate(jp)
C       ENDIF
#endif
      ENDDO

C tendencies with quotas (Chl) would be here

#ifdef SPEAD_DIC
      gTr(iDIC )=gTr(iDIC ) + respir
      gTr(iPIC )=gTr(iPIC ) + exude_PIC
#endif
      gTr(iDOC )=gTr(iDOC ) + exude_DOC
      gTr(iPOC )=gTr(iPOC ) + exude_POC
#ifdef SPEAD_STOICHIOMETRY_OM
      gTr(iDON )=gTr(iDON ) + exude_DON
      gTr(iDOP )=gTr(iDOP ) + exude_DOP
      gTr(iDOFe)=gTr(iDOFe) + exude_DOFe
      gTr(iPON )=gTr(iPON ) + exude_PON
      gTr(iPOP )=gTr(iPOP ) + exude_POP
      gTr(iPOFe)=gTr(iPOFe) + exude_POFe
#ifdef SPEAD_SILICA
      gTr(iPOSi)=gTr(iPOSi) + exude_POSi
#endif
#endif /* spead_stoichiometry_om */

#ifdef DARWIN_ALLOW_CDOM
      exude_CDOM = fracCDOM*exude_DOP
      gTr(iCDOM)=gTr(iCDOM) + exude_CDOM
      gTr(iDOC )=gTr(iDOC )              - R_CP_CDOM*exude_CDOM
#ifdef SPEAD_STOICHIOMETRY_OM
      gTr(iDON )=gTr(iDON )              - R_NP_CDOM*exude_CDOM
      gTr(iDOP )=gTr(iDOP ) - exude_CDOM
      gTr(iDOFe)=gTr(iDOFe)              - R_FeP_CDOM*exude_CDOM
#endif
#endif /* darwin_allow_cdom */

#endif /* ALLOW_DARWIN */

C     STOP 'error'      

      RETURN
      END SUBROUTINE
