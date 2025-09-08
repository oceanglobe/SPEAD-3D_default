#include "DARWIN_OPTIONS.h"

CBOP
C !ROUTINE: DARWIN_SPEAD_PHYTO
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_SPEAD_PHYTO(
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
     U     acom,gcom,
#ifdef SPEAD_CONTINUOUS_TRAIT
     U     a_1t,
#endif
     I     myTime,myIter,myThid)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE
#include "EEPARAMS.h"
#ifdef ALLOW_RADTRANS
#include "RADTRANS_SIZE.h"
#endif
#include "DARWIN_SIZE.h"
#include "DARWIN_INDICES.h"
#include "DARWIN_DIAGS.h"
#include "DARWIN_RADTRANS.h"
#include "DARWIN_PARAMS.h"
#include "DARWIN_TRAITS.h"
#include "DARWIN_TRAITPARAMS.h"

C !INPUT PARAMETERS: ===================================================
C  NH4,NO2,NO3,PO4,FeT,DIC,SiO2 :: darwin biogeochemical tracers
C                               :: inorganic nutrients
C  photoTempFunc,uptakeTempFunc :: temperature response functions
C  PAR    :: PAR in uEin/s/m2
C         :: either non-spectral (tlam=1) or waveband total
C  X      :: biomass (mmolC/m3) of phytoplankton species
C  wght   :: weight of each phytoplankton tracer (continuous-trait)
C  myTime :: current time
C  myIter :: current iteration number
C  myThid :: thread number
C Continuous model: temperature dependance is managed within darwin_plankton (Le Gland, 14/01/2021)
      _RL PAR(nlam)
      _RL photoTempFunc(nplank)
      _RL uptakeTempFunc
#ifdef SPEAD_DIC
      _RL DIC
#endif
      _RL NH4, NO2, NO3, PO4, FeT
#ifdef SPEAD_SILICA
      _RL SiO2
#endif
#ifdef DARWIN_DEBUG
      INTEGER iG, jG, k
      _RL dT
#endif
#ifdef DARWIN_ALLOW_CARBON
      _RL O2, ALK
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
      _RL acom(nPlank), gcom(nPlank)
#ifdef SPEAD_CONTINUOUS_TRAIT
      _RL a_1t(nTrac)
#endif

C !OUTPUT PARAMETERS: ==================================================
C  chlout :: computed acclimated chlorophyll if not dynamic
      _RL chlout(nPhyp)
CEOP

c !LOCAL VARIABLES: ====================================================
      INTEGER jp, jt, l, g
C Error messages for debugging (Le Gland, 16/03/2021)
      CHARACTER*(MAX_LEN_MBUF) msgBuf

C Minimum value for optimal irradiance at a given cell volume
      _RL PARoptmin

C Parameter of the steepness of the limitation by light
C following Edward et al 2015 (Boris Sauterey 22/08/2023)
      _RL alpha

C Continuous model: deriv. from one tracer only (Le Gland, 11/03/2021)
      _RL Qn(nTrac),Qp(nTrac),Qsi(nTrac),Qfe(nTrac)

      _RL regQ

      _RL limitpCO2
      _RL limitNH4, limitNO, limitNO2, limitNO3
      _RL fracNH4, fracNO2, fracNO3
      _RL ntot
      _RL limitn, limitp, limitsi, limitfe, limitnut
      _RL limitI
      _RL ngrow

#ifdef SPEAD_DIC
      _RL uptakeDIC
      _RL consumDIC, consumDIC_PIC
#endif
      _RL uptakeNH4, uptakeNO2, uptakeNO3, uptakeN
      _RL uptakePO4, uptakeSiO2, uptakeFeT
      _RL consumNH4, consumNO2, consumNO3
      _RL consumPO4, consumSiO2, consumFeT

      _RL alpha_I, alpha_I_growth
      _RL PCm, PC
      _RL acclim, chl2c, rhochl, synthChl
      _RL growth
      _RL Ek, EkoverE

      _RL PARtot, PAReff

      _RL tmp

#ifdef DARWIN_ALLOW_DENIT
      _RL denit, denitNH4
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

C Continuous-trait: traits must be allocated at each time and location
C     print*, 'a_PCmax:', a_PCmax, 'b_PCmax:', b_PCmax
C     print*, 'PCmax:', PCmax
C     print*, 'biovol:', biovol
#ifdef SPEAD_CONTINUOUS_TRAIT
C     print*, "TEST: ksatPAR", ksatPAR, kinhPAR
C     CALL DARWIN_GENERATE_ALLOMETRIC(myThid)
C     print*, "TEST: ksatPAR", ksatPAR, kinhPAR
#endif
C     print*, 'a_PCmax:', a_PCmax, 'b_PCmax:', b_PCmax
C     print*, 'PCmax:', PCmax
C     print*, 'biovol:', biovol

C     WRITE(msgBuf,'(A)') 'DARWIN_SPEAD_RATES: step 1 completed:'
C     CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
C    &                    SQUEEZE_RIGHT, myThid )

#ifdef SPEAD_DIC
      consumDIC  = 0.0
      consumDIC_PIC = 0.0
#endif
      consumNH4  = 0.0
      consumNO2  = 0.0
      consumNO3  = 0.0
      consumPO4  = 0.0
      consumSiO2 = 0.0
      consumFeT  = 0.0

#ifdef SPEAD_CONTINUOUS_TRAIT
      DO jp = 1, nTrac
        a_1t(jp) = 0.0 _d 0
      ENDDO
#endif

      DO jp = 1, nPlank
        gcom(jp) = 0.0 _d 0
        acom(jp) = 0.0 _d 0
      ENDDO

      DO jp = 1, nPhyp
        chlout(jp) = 0.0 _d 0
      ENDDO

      DO l=1,darwin_nDiag
        diags(l) = 0.0
      ENDDO

C     PARtot = SUM(PAR)
      PARtot = SUM(PAR) * .2174

C=======================================================================
C==== phytoplankton ====================================================

      DO jp = 1, nPhoto
       IF (isPhoto(jp) .NE. 0) THEN

C==== uptake and nutrient limitation ===================================
C       for quota elements, growth is limiteed by available quota,
C       for non-quota elements, by available nutrients in medium

C       to not use PO4, ..., set ksatPO4=0 and vmaxPO4=0 (if DARWIN_ALLOW_PQUOTA)
C       or R_PC=0 (if not)
C       the result will be limitp = 1, uptakePO4 = 0

c phosphorus
        IF (ksatPO4(jp) .GT. 0 _d 0) THEN
         limitp = PO4/(PO4 + ksatPO4(jp))
        ELSE
         limitp = 1 _d 0
        ENDIF
#ifdef DARWIN_ALLOW_PQUOTA
        regQ = MAX(0., MIN(1., (Qpmax(jp)-Qp(jp))/
     &                         (Qpmax(jp)-Qpmin(jp)) ))
        regQ = regQ**hillnumUptake
        uptakePO4 = vmaxPO4(jp)*limitp*regQ*uptakeTempFunc*X(plank(jp))
c       normalized Droop limitation
        limitp = MIN(1., (1.0-Qpmin(jp)/MAX(Qpmin(jp), Qp(jp)))/
     &                   (1.0-Qpmin(jp)/Qpmax(jp)))
#endif

#ifdef SPEAD_SILICA
c silica
        IF (ksatSiO2(jp) .GT. 0 _d 0) THEN
         limitsi = SiO2/(SiO2 + ksatSiO2(jp))
        ELSE
         limitsi = 1 _d 0
        ENDIF
#ifdef DARWIN_ALLOW_SIQUOTA
        regQ = MAX(0., MIN(1., (Qsimax(jp) - Qsi(jp))/
     &                         (Qsimax(jp) - Qsimin(jp)) ))
        regQ = regQ**hillnumUptake
        uptakeSiO2 = vmaxSiO2(jp)*limitsi*regQ*uptakeTempFunc
     &                                          *X(plank(jp))

c       linear limitation
        limitsi = MAX(0., MIN(1., (Qsi(jp) - Qsimin(jp))/
     &                         (Qsimax(jp) - Qsimin(jp)) ))
#endif
#else /* spead_silica */
        limitsi = 1.0 _d 0
#endif /* spead_silica */

c iron
        IF (ksatFeT(jp) .GT. 0 _d 0) THEN
         limitfe = FeT/(FeT + ksatFeT(jp))
        ELSE
         limitfe = 1 _d 0
        ENDIF
#ifdef DARWIN_ALLOW_FEQUOTA
        regQ = MAX(0., MIN(1., (Qfemax(jp)-Qfe(jp))/
     &                         (Qfemax(jp)-Qfemin(jp)) ))
        regQ = regQ**hillnumUptake
        uptakeFeT = vmaxFeT(jp)*limitfe*regQ*uptakeTempFunc*X(plank(jp))

c       normalized Droop limitation
        limitfe = MIN(1., (1.0-Qfemin(jp)/MAX(Qfemin(jp), Qfe(jp)))/
     &                    (1.0-Qfemin(jp)/Qfemax(jp)))
#endif

c nitrogen
        limitNH4 = useNH4(jp)*NH4/(NH4 + ksatNH4(jp))
        ntot = useNH4(jp) * NH4 / ksatNH4(jp)
     &       + useNO2(jp) * NO2 / ksatNO2(jp)
     &       + useNO3(jp) * NO3 / ksatNO3(jp)
        limitn = ntot / (1 + ntot)
        limitNO = limitn - limitNH4
        IF (NO2 .GT. 0 _d 0) THEN
          limitNO2 = limitNO * useNO2(jp)*NO2
     &             / (NO2 + useNO3(jp)*NO3*ksatNO2(jp)/ksatNO3(jp))
        ELSE
          limitNO2 = 0 _d 0
        ENDIF
        limitNO3 = limitNO - limitNO2
C       normalize to sum 1
        IF (limitn .GT. 0 _d 0) THEN
          fracNH4 = limitNH4/limitn
          fracNO2 = limitNO2/limitn
          fracNO3 = limitNO3/limitn
        ELSE
          fracNH4 = 0 _d 0
          fracNO2 = 0 _d 0
          fracNO3 = 0 _d 0
        ENDIF
C if diazo, all fracN* == 0 but want no N limitation
C Perhaps it should be MIN(1.0, MAX(limitn,diazo(jp)))
        limitn = MIN(1.0, limitn + diazo(jp))
#ifdef DARWIN_ALLOW_NQUOTA
c       have nitrogen quota
        regQ = MAX(0., MIN(1., (Qnmax(jp)-Qn(jp))/
     &                         (Qnmax(jp)-Qnmin(jp)) ))
        regQ = regQ**hillnumUptake
        uptakeNH4=vmaxNH4(jp)*limitNH4*regQ*uptakeTempFunc*X(plank(jp))
        uptakeNO2=vmaxNO2(jp)*limitNO2*regQ*uptakeTempFunc*X(plank(jp))
        uptakeNO3=vmaxNO3(jp)*limitNO3*regQ*uptakeTempFunc*X(plank(jp))
#ifdef DARWIN_ALLOW_FEQUOTA
#ifdef DARWIN_NITRATE_FELIMIT
        uptakeNO3 = uptakeNO3 * limitfe
#endif
#endif
        uptakeN = MAX(uptakeNH4 + uptakeNO2 + uptakeNO3,
     &            vmaxN(jp)*regQ*uptakeTempFunc*X(plank(jp))*diazo(jp))

c       linear limitation
        limitn = MAX(0., MIN(1., (Qn(jp) - Qnmin(jp))/
     &                        (Qnmax(jp) - Qnmin(jp)) ))
#endif /* darwin_allow_nquota */
        ngrow = 1.0 _d 0

        limitnut = MIN(limitn, limitp, limitsi)
#ifndef DARWIN_ALLOW_FEQUOTA
        limitnut = MIN(limitnut, limitfe)
#endif

        limitpCO2 = 1.0 _d 0

C==== growth ===========================================================
#ifdef DARWIN_ALLOW_GEIDER

        alpha_I = 0 _d 0
        DO l = 1, nlam
          alpha_I = alpha_I + alphachl(jp,l)*PAR(l)
        ENDDO
C       NB: for quota, PCmax(j) = Vmax_c(j)
        PCm = PCmax(jp)*limitnut*photoTempFunc(jp)*limitpCO2

        IF (PCm .GT. 0.0) THEN
          acclim = MAX(chl2cmin(jp), MIN(chl2cmax(jp),
     &             chl2cmax(jp)/(1+(chl2cmax(jp)*alpha_I)/(2*PCm)) ))
        ELSE
          acclim = chl2cmin(jp)
        ENDIF

#ifdef DARWIN_ALLOW_CHLQUOTA
C       quotas are already relative to carbon
        chl2c = QChl(jp)
#else
        chl2c = acclim
#endif

        alpha_I_growth = alpha_I
C a la quota
#ifdef DARWIN_ALLOW_FEQUOTA
        alpha_I_growth = alpha_I_growth*limitfe
#endif

C       carbon-specific growth rate
C       PC = PCm*(1-EXP(-alpha_I_growth*chl2c/MAX(EPS, PCm)))
        IF (PCm .GT. 0.0 .AND. PARtot .GT. PARmin) THEN
          PC = PCm*(1-EXP(-alpha_I_growth*chl2c/PCm))
        ELSE
          PC = 0.0 _d 0
        ENDIF

        IF (inhibGeider(jp) .GT. 0.0) THEN
C         "total" PAR:
          tmp = alpha_I/alpha_mean(jp)
          Ek = PCm/(chl2c*alpha_mean(jp))
          EkoverE = Ek / tmp
          IF (tmp .GE. Ek) THEN
            PC = PC*EkoverE*inhibGeider(jp)
          ENDIF
        ENDIF

#else /* NOT darwin_allow_geider */

        IF (PARtot .GT. 0.2174*PARmin*
     &                  (biovol(jp)**b_PARoptmin(group(jp)))) THEN
C       IF (PARtot .GT. PARmin) THEN
CC        only 1 waveband without DARWIN_ALLOW_GEIDER
C         limitI = (1.0 _d 0 - EXP(-PARtot*ksatPAR(jp)))*
C    &             EXP(-PARtot*kinhPAR(jp)) * normI(jp)
C From Edwards et al. (2015)
          alpha  = 2.0/PARopt(jp)
CC        tmp = MIN(1.0,MAX(0.05,limitnut*photoTempFunc(j)*limitpCO2*
CC   &          (biovol(jp)**b_PCmax(group(jp)))))
CC        tmp = MIN(1.0,MAX(0.1,0.5*limitnut*photoTempFunc(jp)*limitpCO2*
CC   &          (biovol(jp)**b_PCmax(group(jp)))))
CC        alpha  = 1.0/(tmp*PARopt(j)) ! Le Gland (2024-07-02)
CC        PAReff = sum(PAR) - PARmin*(biovol(jp)**b_PARoptmin(group(j)))
CC        PAReff = MAX(0.0,PAReff)
          limitI = PARtot/(PARtot**2/(alpha*PARopt(jp)**2)+(1-2/
     &             (alpha*PARopt(jp)))*PARtot+1/alpha)
CC        limitI = PAReff/(PAReff**2/(alpha*PARopt(jp)**2)+(1-2/
CC   &             (alpha*PARopt(jp)))*PAReff+1/alpha)
CC        print*, 'limitI:', limitI, 'alpha:', alpha
          PAReff = 0.2174*PARmin*(biovol(jp)**b_PARoptmin(group(jp)))
          limitI = MAX(0.0, (limitI-alpha*PAReff) / 
     &             MAX(0.8,(1-alpha*PAReff)))
CC        print*, 'limitI:', limitI
          PC = PCmax(jp)*limitnut*limitI*photoTempFunc(jp)*limitpCO2
        ELSE
          PC = 0.0 _d 0
        ENDIF
        synthChl = 0.0

#endif /* darwin_allow_geider */

        growth = PC*ngrow*X(plank(jp))
#ifdef SPEAD_DIC
        uptakeDIC = growth
#endif

C non-quota elements are taken up with fixed stoichiometry
#ifndef DARWIN_ALLOW_NQUOTA
        uptakeN = growth*R_NC(jp)
        uptakeNH4 = uptakeN*fracNH4
        uptakeNO2 = uptakeN*fracNO2
        uptakeNO3 = uptakeN*fracNO3
#endif
#ifndef DARWIN_ALLOW_PQUOTA
        uptakePO4 = growth*R_PC(jp)
#endif
#ifndef DARWIN_ALLOW_SIQUOTA
        uptakeSiO2 = growth*R_SiC(jp)
#endif
#ifndef DARWIN_ALLOW_FEQUOTA
        uptakeFeT = growth*R_FeC(jp)
#endif

C==== chlorophyll ======================================================
#ifdef DARWIN_ALLOW_GEIDER
# ifdef DARWIN_ALLOW_CHLQUOTA
#  ifdef DARWIN_ALLOW_NQUOTA
C       Geider 1998
        IF (PARtot .GT. PARmin) THEN
         IF (alpha_I*chl2c .GT. 0.0 _d 0) THEN
C          rhochl = Chl2Nmax/(alpha_I*chl2c)*ngrow ???
          rhochl = Chl2Nmax*PC*ngrow/(alpha_I*chl2c)
         ELSE
          rhochl = Chl2Nmax
         ENDIF
         synthChl = rhochl*uptakeN
        ELSE
         synthChl = 0 _d 0
        ENDIF
CTODO: this seems wrong!  Should only happen when PARtot.GT.PARmin
        uptakeDIC = uptakeDIC - synthcost*uptakeN
#  else /* not DARWIN_ALLOW_NQUOTA */
#   ifdef DARWIN_GEIDER_RHO_SYNTH
        IF (alpha_I .GT. 0.0 _d 0 .AND. acclim .GT. 0.0 _d 0) THEN
          rhochl = Chl2Cmax(jp)*PC*ngrow/(alpha_I*acclim)
        ELSE
          rhochl = 0.0 _d 0    ! should be Chl2Cmax(j) ?????
        ENDIF
        synthChl = rhochl*growth +
     &                    acclimtimescl(jp)*(acclim-chl2c)*X(plank(jp))
#   else
        synthChl = acclim*growth +
     &                    acclimtimescl(jp)*(acclim-chl2c)*X(plank(jp))
#   endif
#  endif /* not DARWIN_ALLOW_NQUOTA */
# else /* not DARWIN_ALLOW_CHLQUOTA */
C       chlout(j) = X(j)*Qc(j)*chl2c
C Could be done out of the loop if chl2c is a constant
        chlout(plank(jp)) = chlout(plank(jp))
     &                   + X(plank(jp))*chl2c*wght(jp)
        synthChl = 0.0
# endif /* DARWIN_ALLOW_CHLQUOTA */
#endif /* DARWIN_ALLOW_GEIDER */
C=======================================================================
#ifdef SPEAD_DIC
        consumDIC_PIC = consumDIC_PIC + uptakeDIC*R_PICPOC(jp)*wght(jp)
        consumDIC  = consumDIC  + uptakeDIC *wght(jp)
#endif
        consumNH4  = consumNH4  + uptakeNH4 *wght(jp)
        consumNO2  = consumNO2  + uptakeNO2 *wght(jp)
        consumNO3  = consumNO3  + uptakeNO3 *wght(jp)
        consumPO4  = consumPO4  + uptakePO4 *wght(jp)
        consumSiO2 = consumSiO2 + uptakeSiO2*wght(jp)
        consumFeT  = consumFeT  + uptakeFeT *wght(jp)

C       gTr(iplank(plank(jp)))=gTr(iplank(plank(jp)))+growth*wght(jp)
        acom(plank(jp)) = acom(plank(jp)) + growth*wght(jp)
        gcom(plank(jp)) = gcom(plank(jp)) + growth*wght(jp)
#ifdef SPEAD_CONTINUOUS_TRAIT
        a_1t(jp) = a_1t(jp) + growth
#endif

C gTr and a_1t should be set around here

        diags(iPP) = diags(iPP) + growth*wght(jp)
#ifdef DARWIN_DIAG_PERTYPE
        diags(iPPplank+plank(jp)-1) = diags(iPPplank+plank(jp)-1)
     &                             + growth*wght(jp)
#endif
        IF (diazo(jp) .GT. 0.0 _d 0) THEN
         diags(iNfix) = diags(iNfix)
     &   + (uptakeN-uptakeNH4-uptakeNO2-uptakeNO3)*wght(jp)
        ENDIF

C=======================================================================

#ifdef DARWIN_DEBUG
C       IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG) THEN
CC        print*,'j', j, 'plank(j)', plank(j), 'uptakeDIC', uptakeDIC,
          print*,'jp', jp, 'plank(jp)', plank(jp),
     &           'wght(jp)', wght(jp), 'X(jp)', X(jp),
#ifdef SPEAD_CONTINUOUS_TRAIT
     &   'a_1t(jp)', a_1t(jp),
#endif
     &   'acom(plank(jp))', acom(plank(jp)),
     &   'gcom(plank(jp))', gcom(plank(jp))

          print*,'PC', PC, 'PCmax(jp)', PCmax(jp), 'limitnut', limitnut,
     &   'limitI', limitI, 'PAReff', PAReff,
     &   'photoTempFunc(jp)', photoTempFunc(jp),
     &   'limitpCO2', limitpCO2, 'PARtot', PARtot, 'PARopt', PARopt(jp)

          print*,'limitn', limitn, 'limitNH4', limitNH4, 'limitNO',
     &    limitNO, 'limitNO2', limitNO2, 'limitNO3', limitNO3,
     &   'limitp', limitp, 'limitsi', limitsi,
     &   'limitfe', limitfe, 'normI(jp)', normI(jp)

          print*, 'ksatNH4(jp)',ksatNH4(jp),'ksatNO2(jp)',ksatNO2(jp),
     &    'ksatNO3(jp)',ksatNO3(jp),'ksatPO4',ksatPO4(jp),
     &    'ksatSiO2(jp)',ksatSiO2(jp),'ksatPAR(jp)',ksatPAR(jp),
     &    'ksatFeT(jp)',ksatFeT(jp)
C       ENDIF
#endif

C tendencies with quotas would be here

#ifdef DARWIN_DEBUG
CC      IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG) THEN
         print*,'uptake',myiter,k,jp,
C    &     uptakeDIC,
     &     uptakeNH4,
     &     uptakeNO2,
     &     uptakeNO3,
     &     uptakeN,
     &     uptakePO4,
C    &     uptakeSiO2,
     &     uptakeFeT
CC      ENDIF
#endif

C      isPhoto(j)
       ENDIF

C     j
      ENDDO

#ifdef SPEAD_DIC
      gTr(iDIC)  = gTr(iDIC)  - consumDIC - consumDIC_PIC
#endif
      gTr(iNH4)  = gTr(iNH4)  - consumNH4
      gTr(iNO2)  = gTr(iNO2)  - consumNO2
      gTr(iNO3)  = gTr(iNO3)  - consumNO3 
      gTr(iPO4)  = gTr(iPO4)  - consumPO4 
      gTr(iFeT)  = gTr(iFeT)  - consumFeT 
#ifdef SPEAD_SILICA
      gTr(iSiO2) = gTr(iSiO2) - consumSiO2
#endif
#ifdef DARWIN_ALLOW_CARBON
      gTr(iO2)   = gTr(iO2)   + R_OP*consumPO4
      gTr(iALK)  = gTr(iALK)  + consumNO3 - 2.0*consumDIC_PIC
#endif
C     DO jp = 1, nPhyp
C       gTr(ic+iplank(jp)-1) = gTr(ic+iplank(jp)) + acom(jp)
C     ENDDO

      diags(iConsDIN) = consumNH4 + consumNO2 + consumNO3
      diags(iConsPO4) = consumPO4
      diags(iConsFe)  = consumFeT
#ifdef SPEAD_SILICA
      diags(iConsSi)  = consumSiO2
#endif

C If extreme mean traits occur, growth rate should be zero. If NaNs
C occur, all derivatives should be forced to be zero to prevent failure
C If this is all created by transport, future transport may solve it
C (Le Gland, 28/05/2021)
C Change position in the code to still allow for mortality, grazing,
C and chemical reactions (Le Gland, 01/06/2021)
      DO jp = 1, nplank
        IF (gcom(jp) /= gcom(jp) .OR. gcom(jp) .LT. -1.0 _d -6) THEN
#ifdef DARWIN_DEBUG
CC        IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG) THEN
            print*,'SPEAD is buggy (NaN): all derivatives are set to',
     &             ' zero to prevent failure', 'coordinates', iG, jG, k,
     &             ' gcom: ', gcom(jp)
CC        ENDIF
#endif
C         chlout = 0 _d 0
C         diags = 0 _d 0
C         gTr  = 0 _d 0
C         acom = 0 _d 0
C         gcom = 0 _d 0
C         a_1t = 0 _d 0
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE
