#include "DARWIN_OPTIONS.h"

CBOP
C !ROUTINE: DARWIN_SPEAD_REMIN
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_SPEAD_REMIN(
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

C !DESCRIPTION: REMINERALIZATION OF ORGANIC MATTER (DOM+POM) BY BACTERIA
C               EXPLICIT OR PARAMETERIZED

C !USES: ===============================================================
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_INDICES.h"
#include "DARWIN_DIAGS.h"
#include "DARWIN_PARAMS.h"
#include "DARWIN_TRAITS.h"

C !INPUT PARAMETERS: ===================================================
C  DOC,POC,DON,DOP,DOFe,PON,POP,POFe, :: darwin biogeochemical tracers
C  POSi,CDOM                          :: organic matter 
C  reminTempFunc :: temperature response functions
C  X      :: biomass (mmolC/m3) of phytoplankton species
C  wght   :: weight of each phytoplankton tracer (continuous-trait)
C  myTime :: current time
C  myIter :: current iteration number
C  myThid :: thread number
      _RL PAR(nlam)
      _RL reminTempFunc
      _RL X(nPlank)
      _RL wght(nTrac)
      _RL NH4, NO2, NO3
#ifdef SPEAD_DIC
      _RL PIC
#endif
      _RL DOC, POC
#ifdef SPEAD_STOICHIOMETRY_OM
      _RL DON, DOP, DOFe, PON, POP, POFe
#ifdef SPEAD_SILICA
      _RL POSi
#endif
#endif /* spead_stoichiometry_om */
#ifdef DARWIN_DEBUG
      INTEGER iG, jG, k
      _RL dT
#endif
#ifdef DARWIN_ALLOW_CARBON
      _RL O2
#endif
#ifdef DARWIN_ALLOW_CDOM
      _RL CDOM
#endif
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
CEOP

c !LOCAL VARIABLES: ====================================================
      INTEGER jp
C Error messages for debugging (Le Gland, 16/03/2021)
      CHARACTER*(MAX_LEN_MBUF) msgBuf

      _RL muPOC, muDOC 
#ifdef SPEAD_STOICHIOMETRY_OM
      _RL muPON, muPOP, muPOFe, muDON, muDOP, muDOFe
#endif
      _RL muO, mu

      _RL uptakePOC, uptakeDOC, uptakePOP, uptakePOFe
#ifdef SPEAD_STOICHIOMETRY_OM
      _RL uptakePON, uptakePOP, uptakePOFe
      _RL uptakeDON, uptakeDOP, uptakeDOFe
#endif
      _RL uptakeO2
      _RL uptakeNO3

      _RL respPOC, respPON, respPOP, respPOFe, respPOSi
      _RL respDOC, respDON, respDOP, respDOFe

      _RL hydrolPOC, solubilPOC
#ifdef SPEAD_STOICHIOMETRY_OM 
      _RL hydrolPON, hydrolPOP, hydrolPOFe
      _RL solubilPON, solubilPOP, solubilPOFe
#endif

      _RL consumPOC, consumDOC
#ifdef SPEAD_STOICHIOMETRY_OM
      _RL consumPON, consumPOP, consumPOFe, consumPOSi
      _RL consumDON, consumDOP, consumDOFe
#endif
      _RL consumNO3
      _RL consumO2

      _RL PCm
      _RL growth

      _RL reminDOC, reminDON, reminDOP, reminDOFe
      _RL reminPOC, reminPON, reminPOP, reminPOSi, reminPOFe
#ifdef SPEAD_DIC
      _RL disscPIC
#endif

      _RL prodNO2, prodNO3

#ifdef DARWIN_ALLOW_CDOM
      _RL CDOM
      _RL reminPOC_CDOM, reminPON_CDOM, reminPOP_CDOM, reminPOFe_CDOM
      _RL degrCDOM_DOC,  degrCDOM_DON,  degrCDOM_DOP,  degrCDOM_DOFe
#endif

#ifdef DARWIN_ALLOW_DENIT
      _RL denit, denitNH4
#endif

      _RL PARtot

C     WRITE(msgBuf,'(A)') 'DARWIN_SPEAD_RATES: step 1 completed:'
C     CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
C    &                    SQUEEZE_RIGHT, myThid )

      consumNO3     = 0.0 _d 0
      solubilPOC    = 0.0 _d 0
      consumDOC     = 0.0 _d 0
      consumPOC     = 0.0 _d 0
#ifdef SPEAD_STOICHIOMETRY_OM
      solubilPON    = 0.0 _d 0
      solubilPOP    = 0.0 _d 0
      solubilPOFe   = 0.0 _d 0
      consumDON     = 0.0 _d 0
      consumDOP     = 0.0 _d 0
      consumDOFe    = 0.0 _d 0
      consumPON     = 0.0 _d 0
      consumPOP     = 0.0 _d 0
      consumPOFe    = 0.0 _d 0
#ifdef SPEAD_SILICA
      consumPOSi    = 0.0 _d 0
#endif
#endif /* spead_stoichiometry_om */
      consumNO3     = 0.0 _d 0
#ifdef DARWIN_ALLOW_CARBON
      consumO2      = 0.0 _d 0
#endif
      reminPON      = 0.0 _d 0
      reminPOP      = 0.0 _d 0
      reminPOC      = 0.0 _d 0
      reminPOFe     = 0.0 _d 0
      reminPOSi     = 0.0 _d 0
      reminDON      = 0.0 _d 0
      reminDOP      = 0.0 _d 0
      reminDOC      = 0.0 _d 0
      reminDOFe     = 0.0 _d 0
      prodNO2       = 0.0 _d 0
      prodNO3       = 0.0 _d 0

      PARtot = SUM(PAR)


C=======================================================================
C==== bacteria =========================================================

      DO jp = 1, nPlank
       IF (bactType(jp) .NE. 0) THEN

#ifdef DARWIN_ALLOW_CARBON
        uptakeO2   = 0.0 _d 0
#endif
        uptakeNO3  = 0.0 _d 0
        uptakeDOC  = 0.0 _d 0
        uptakePOC  = 0.0 _d 0
        hydrolPOC  = 0.0 _d 0
#ifdef SPEAD_STOICHIOMETRY_OM
        uptakeDON  = 0.0 _d 0
        uptakeDOP  = 0.0 _d 0
        uptakeDOFe = 0.0 _d 0
        uptakePON  = 0.0 _d 0
        uptakePOP  = 0.0 _d 0
        uptakePOFe = 0.0 _d 0
        hydrolPON  = 0.0 _d 0
        hydrolPOP  = 0.0 _d 0
        hydrolPOFe = 0.0 _d 0
#endif
        respPOC    = 0.0 _d 0
        respPON    = 0.0 _d 0
        respPOP    = 0.0 _d 0
        respPOFe   = 0.0 _d 0
        respDOC    = 0.0 _d 0
        respDON    = 0.0 _d 0
        respDOP    = 0.0 _d 0
        respDOFe   = 0.0 _d 0
        growth     = 0.0 _d 0

        IF (isAerobic(jp) .NE. 0) THEN
#ifdef DARWIN_ALLOW_CARBON
          muO = yieldO2(jp)*pcoefO2*O2
#else
C No limitation by oxygen
          muO = yieldO2(jp)*pcoefO2*(10.0**6)
#endif
        ELSEIF (isDenit(jp) .NE. 0) THEN
          muO = yieldNO3(jp)*pmaxDIN*NO3/(NO3 + ksatDIN)*reminTempFunc
        ENDIF

C       POM-consuming (particle-associated)
        IF (bactType(jp) .EQ. 1) THEN

          PCm = yield(jp)*PCmax(jp)*reminTempFunc
          muPOC  = PCm*POC/(POC + ksatPOC(jp))
          mu = muPOC
#ifdef SPEAD_STOICHIOMETRY_OM
          muPON  = PCm*PON/(PON + ksatPON(jp))
          muPOP  = PCm*POP/(POP + ksatPOP(jp))
          muPOFe = PCm*POFe/(POFe + ksatPOFe(jp))
          mu = MIN(mu, muPON, muPOP, muPOFe, muO)
#endif

          growth = mu*X(plank(jp))

          uptakePOC = alpha_hydrol*growth/yield(jp)
#ifdef SPEAD_STOICHIOMETRY_OM
          uptakePON  = uptakePOC*R_NC(jp)
          uptakePOP  = uptakePOC*R_PC(jp)
          uptakePOFe = uptakePOC*R_FeC(jp)
#endif
C         O2/NO3 is only used for the part of POC that is metabolized:
#ifdef DARWIN_ALLOW_CARBON
          uptakeO2  = isAerobic(jp)*growth/yieldO2(jp)
#endif
          uptakeNO3 = isDenit(jp)  *growth/yieldNO3(jp)

C         This is the part of POM that is hydrolized into DOM:
          hydrolPOC = (alpha_hydrol-1)*growth/yield(jp)
#ifdef SPEAD_STOICHIOMETRY_OM
          hydrolPON  = hydrolPOC*R_NC(jp)
          hydrolPOP  = hydrolPOC*R_PC(jp)
          hydrolPOFe = hydrolPOC*R_FeC(jp)
#endif

C         These are the bacteria products for remineralization of POM:
          respPOC = growth*(1/yield(jp)-1)
          respPON  = respPOC*R_NC(jp)
          respPOP  = respPOC*R_PC(jp)
          respPOFe = respPOC*R_FeC(jp)

C       DOM-consuming (free-living):
        ELSEIF (bactType(jp) .EQ. 2) THEN

          PCm = yield(jp)*PCmax(jp)*reminTempFunc
          muDOC  = PCm*DOC/(DOC + ksatDOC(jp))
          mu = muDOC
#ifdef SPEAD_STOICHIOMETRY_OM
          muDON  = PCm*DON/(DON + ksatDON(jp))
          muDOP  = PCm*DOP/(DOP + ksatDOP(jp))
          muDOFe = PCm*DOFe/(DOFe + ksatDOFe(jp))
          mu = MIN(mu, muDON, muDOP, muDOFe, muO)
#endif

          growth = mu*X(plank(jp))

          uptakeDOC = growth/yield(jp)
#ifdef SPEAD_STOICHIOMETRY_OM
          uptakeDON  = uptakeDOC*R_NC(jp)
          uptakeDOP  = uptakeDOC*R_PC(jp)
          uptakeDOFe = uptakeDOC*R_FeC(jp)
#endif
#ifdef DARWIN_ALLOW_CARBON
          uptakeO2 = isAerobic(jp)*growth/yieldO2(jp)
#endif
          uptakeNO3 = isDenit(jp)*growth/yieldNO3(jp)

C         DOC respired to DIC
          respDOC = growth*(1/yield(jp)-1)
          respDON  = respDOC*R_NC(jp)
          respDOP  = respDOC*R_PC(jp)
          respDOFe = respDOC*R_FeC(jp)

        ENDIF

#ifdef DARWIN_DIAG_PERTYPE
        diags(iPPplank+plank(jp)-1) = diags(iPPplank+plank(jp)-1)
     &                             + growth*wght(jp)
#endif

C       gTr(iplank(plank(jp)))=gTr(iplank(plank(jp)))+growth*wght(jp)
        acom(plank(jp)) = acom(plank(jp)) + growth*wght(jp)
        gcom(plank(jp)) = gcom(plank(jp)) + growth*wght(jp)
#ifdef SPEAD_CONTINUOUS_TRAIT
        a_1t(jp) = a_1t(jp) + growth
#endif

C==== Cumulative consum, remin, and prod ===============================
        consumNO3  = consumNO3  + uptakeNO3*wght(jp)

C       add B consum and accumulating remin, and prod:
        consumO2 = consumO2 + uptakeO2*wght(jp)

        consumDOC  = consumDOC  + uptakeDOC *wght(jp)
#ifdef SPEAD_STOICHIOMETRY_OM
        consumDON  = consumDON  + uptakeDON *wght(jp)
        consumDOP  = consumDOP  + uptakeDOP *wght(jp)
        consumDOFe = consumDOFe + uptakeDOFe*wght(jp)
#endif

        consumPOC  = consumPOC  + uptakePOC *wght(jp)
#ifdef SPEAD_STOICHIOMETRY_OM
        consumPON  = consumPON  + uptakePON *wght(jp)
        consumPOP  = consumPOP  + uptakePOP *wght(jp)
        consumPOFe = consumPOFe + uptakePOFe*wght(jp)
#endif

        reminPOC  = reminPOC  + respPOC *wght(jp)
        reminPON  = reminPON  + respPON *wght(jp)
        reminPOP  = reminPOP  + respPOP *wght(jp)
        reminPOFe = reminPOFe + respPOFe*wght(jp)

        solubilPOC  = solubilPOC  + hydrolPOC *wght(jp)
#ifdef SPEAD_STOICHIOMETRY_OM
        solubilPON  = solubilPON  + hydrolPON *wght(jp)
        solubilPOP  = solubilPOP  + hydrolPOP *wght(jp)
        solubilPOFe = solubilPOFe + hydrolPOFe*wght(jp)
#endif

        reminDOC  = reminDOC  + respDOC *wght(jp)
        reminDON  = reminDON  + respDON *wght(jp)
        reminDOP  = reminDOP  + respDOP *wght(jp)
        reminDOFe = reminDOFe + respDOFe*wght(jp)

       ENDIF
C     j loop end
      ENDDO

C=======================================================================
C=======================================================================

C parameterized remineralization; want to set all K except KPOSi to zero
C if running with bacteria
      respDOC  = reminTempFunc*KDOC *DOC
      respPOC  = reminTempFunc*KPOC *POC
#ifdef SPEAD_STOICHIOMETRY_OM
      respDON  = reminTempFunc*KDON *DON
      respDOP  = reminTempFunc*KDOP *DOP
      respDOFe = reminTempFunc*KDOFe*DOFe
      respPON  = reminTempFunc*KPON *PON
      respPOP  = reminTempFunc*KPOP *POP
      respPOFe = reminTempFunc*KPOFe*POFe
#ifdef SPEAD_SILICA
      respPOSi = reminTempFunc*KPOSi*POSi
#endif
#else /* NOT spead_stoichiometry_om */
C The ratios must be the same for all species
      respDON  = respDOC*R_NC(1)
      respDOP  = respDOC*R_PC(1)
      respDOFe = respDOC*R_FeC(1)
      respPON  = respPOC*R_NC(1)
      respPOP  = respPOC*R_PC(1)
      respPOFe = respPOC*R_FeC(1)
#ifdef SPEAD_SILICA
      respPOSi = respPOC*R_SiC(1)
#endif
#endif /* spead_stoichiometry_om */

      consumDOC  = consumDOC  + respDOC
      consumPOC  = consumPOC  + respPOC
#ifdef SPEAD_STOICHIOMETRY_OM
      consumDON  = consumDON  + respDON
      consumDOP  = consumDOP  + respDOP
      consumDOFe = consumDOFe + respDOFe
      consumPON  = consumPON  + respPON
      consumPOP  = consumPOP  + respPOP
      consumPOFe = consumPOFe + respPOFe
#ifdef SPEAD_SILICA
      consumPOSi = consumPOSi + respPOSi
#endif
#endif /* spead_stoichiometry_om */

      reminDOC  = reminDOC  + respDOC
      reminDON  = reminDON  + respDON
      reminDOP  = reminDOP  + respDOP
      reminDOFe = reminDOFe + respDOFe
      reminPOC  = reminPOC  + respPOC
      reminPON  = reminPON  + respPON
      reminPOP  = reminPOP  + respPOP
      reminPOFe = reminPOFe + respPOFe
#ifdef SPEAD_SILICA
      reminPOSi = reminPOSi + respPOSi
#endif

#ifdef DARWIN_ALLOW_CARBON
      consumO2  = consumO2  + respDOP*R_OP
#ifndef DARWIN_ALLOW_CDOM
      consumO2  = consumO2  + respPOP*R_OP
#endif
#endif

#ifdef SPEAD_DIC
      disscPIC = Kdissc*PIC
#endif

c nitrogen chemistry
c NH4 -> NO2 -> NO3 by bacterial action, parameterized
      prodNO2 = knita*NH4
      prodNO3 = knitb*NO2
      IF (PAR_oxi .NE. 0.0 _d 0) THEN
        prodNO2 = prodNO2*MAX(0.0, 1.0 - PARtot/PAR_oxi)
        prodNO3 = prodNO3*MAX(0.0, 1.0 - PARtot/PAR_oxi)
      ENDIF

#ifdef DARWIN_ALLOW_CDOM
      reminPOP_CDOM = fracCDOM*reminPOP
      reminPOC_CDOM = R_CP_CDOM*reminPOP_CDOM
      reminPON_CDOM = R_NP_CDOM*reminPOP_CDOM
      reminPOFe_CDOM = R_FeP_CDOM*reminPOP_CDOM
c degradation of  CDOM - high when bleached by light
      degrCDOM_DOP = reminTempFunc*CDOM*
     &     (CDOMdegrd+CDOMbleach*MIN(1.0 _d 0, PARtot/PARCDOM))
      degrCDOM_DOC  = R_CP_CDOM  * degrCDOM_DOP
      degrCDOM_DON  = R_NP_CDOM  * degrCDOM_DOP
      degrCDOM_DOFe = R_FeP_CDOM * degrCDOM_DOP
#endif

#ifdef DARWIN_ALLOW_DENIT
      IF (O2 .LT. O2crit .AND. NO3 .LT. NO3crit) THEN
        consumDOC = 0. _d 0
        consumDOP  = 0. _d 0
        consumDON  = 0. _d 0
        consumDOFe = 0. _d 0
        consumPOC = 0. _d 0
        consumPOP  = 0. _d 0
        consumPON  = 0. _d 0
        consumPOFe = 0. _d 0
        reminDOC = 0. _d 0
        reminDOP  = 0. _d 0
        reminDON  = 0. _d 0
        reminDOFe = 0. _d 0
        reminPOC = 0. _d 0
        reminPOP  = 0. _d 0
        reminPON  = 0. _d 0
        reminPOFe = 0. _d 0
#ifdef DARWIN_ALLOW_CDOM
        reminPOC_cdom = 0. _d 0
        reminPOP_cdom  = 0. _d 0
        reminPON_cdom  = 0. _d 0
        reminPOFe_cdom = 0. _d 0
        degrCDOM_DOP = reminTempFunc*CDOM*
     &                 CDOMbleach*MIN(1.0 _d 0, PARtot/PARCDOM)
        degrCDOM_DON  = R_NP_CDOM  * degrCDOM_DOP
        degrCDOM_DOFe = R_FeP_CDOM * degrCDOM_DOP
        degrCDOM_DOC  = R_CP_CDOM  * degrCDOM_DOP
#endif
      ENDIF
#endif /* darwin_allow_denit */

C==== apply tendencies =================================================

#ifdef DARWIN_ALLOW_CARBON
c loss of O2 by remineralization
      IF (O2 .GT. O2crit) THEN
        gTr(iO2)=gTr(iO2) - consumO2
      ENDIF
      gTr(iALK)=gTr(iALK) - prodNO3 + 2.0 _d 0*disscPIC
#endif /* DARWIN_ALLOW_CARBON */

#ifdef SPEAD_DIC
      gTr(iDIC )=gTr(iDIC ) + reminDOC + disscPIC
#endif
      gTr(iNH4 )=gTr(iNH4 ) + reminDON - prodNO2
      gTr(iNO2 )=gTr(iNO2 ) + prodNO2 - prodNO3
      gTr(iNO3 )=gTr(iNO3 ) + prodNO3
      diags(iDenitN) = 0.0
#ifdef DARWIN_ALLOW_DENIT
      IF (O2 .LT. O2crit) THEN
        denitNH4 = reminDON
        denit = denit_NP*reminDOP
#ifndef DARWIN_ALLOW_CDOM
        denitNH4 = denitNH4 + reminPON
        denit = denit + denit_NP*reminPOP
#endif
        diags(iDenit) = denit
        gTr(iNH4)=gTr(iNH4) - denitNH4
        gTr(iNO3)=gTr(iNO3) - denit_NO3/denit_np*denit
        gTr(iALK)=gTr(iALK) + denit_NO3/denit_np*denit
        diags(iDenitN) = denitNH4 + denit_NO3/denit_np*denit
      ENDIF
#endif /* darwin_allow_denit */

      gTr(iPO4 )=gTr(iPO4 ) + reminDOP
      gTr(iFeT )=gTr(iFeT ) + reminDOFe
#ifdef SPEAD_SILICA
      gTr(iSiO2)=gTr(iSiO2) + reminPOSi
#endif

C     DOC is created by #4 PA-assoc solubilization and consumed by #5
      gTr(iDOC )=gTr(iDOC ) + solubilPOC - consumDOC
#ifdef SPEAD_STOICHIOMETRY_OM
      gTr(iDON )=gTr(iDON ) + solubilPON - consumDON
      gTr(iDOP )=gTr(iDOP ) + solubilPOP - consumDOP
      gTr(iDOFe)=gTr(iDOFe) + solubilPOFe - consumDOFe
#endif

#ifdef SPEAD_DIC
      gTr(iPIC )=gTr(iPIC ) - disscPIC
#endif
      gTr(iPOC )=gTr(iPOC ) - consumPOC
#ifdef SPEAD_STOICHIOMETRY_OM
      gTr(iPON )=gTr(iPON ) - consumPON
      gTr(iPOP )=gTr(iPOP ) - consumPOP
      gTr(iPOFe)=gTr(iPOFe) - consumPOFe
#ifdef SPEAD_SILICA
      gTr(iPOSi)=gTr(iPOSi) - consumPOSi
#endif
#endif /* spead_stoichiometry_om */

#ifdef DARWIN_ALLOW_CDOM
      gTr(iDOC )=gTr(iDOC ) + reminPOC  - reminPOC_CDOM  + degrCDOM_DOC
#ifdef SPEAD_STOICHIOMETRY_OM
      gTr(iDON )=gTr(iDON ) + reminPON  - reminPON_CDOM  + degrCDOM_DON
      gTr(iDOP )=gTr(iDOP ) + reminPOP  - reminPOP_CDOM  + degrCDOM_DOP
      gTr(iDOFe)=gTr(iDOFe) + reminPOFe - reminPOFe_CDOM + degrCDOM_DOFe
#endif
      gTr(iCDOM)=gTr(iCDOM) + reminPOP_CDOM - degrCDOM_DOP
#else /* NOT darwin_allow_cdom */
#ifdef SPEAD_DIC
      gTr(iDIC )=gTr(iDIC ) + reminPOC
#endif
      gTr(iNH4 )=gTr(iNH4 ) + reminPON
      gTr(iPO4 )=gTr(iPO4 ) + reminPOP
      gTr(iFeT )=gTr(iFeT ) + reminPOFe
#endif /* darwin_allow_cdom */

      diags(iConsDIN) = diags(iConsDIN) + consumNO3

C     print*, '----------  DARWIN_SPEAD_REMIN:  ----------'
C     print*, 'reminPON' , reminPON , 'reminDON' , reminDON
C     print*, 'consumPOC', consumPOC, 'consumDOC', consumDOC
C     print*, 'solubilPOC', solubilPOC, 'growth', growth
      
C     STOP 'error'

      RETURN
      END SUBROUTINE
