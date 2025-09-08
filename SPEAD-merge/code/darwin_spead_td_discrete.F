#include "DARWIN_OPTIONS.h"

CBOP
C !ROUTINE: DARWIN_SPEAD_TD_DISCRETE
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_SPEAD_TD_DISCRETE(
     I     X,
     I     gcom,
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
C  gcom   :: gross growth rate of each species
C  myTime    :: current time
C  myIter    :: current iteration number
C  myThid    :: thread number
      _RL X(nPlank)
      _RL gcom(nTrac)
      INTEGER myThid, myIter
      _RL myTime

C !INPUT/OUTPUT PARAMETERS: ============================================
C  gTr    :: accumulates computed tendencies
      _RL gTr(nDarwin)

C !OUTPUT PARAMETERS: ==================================================
CEOP

c !LOCAL VARIABLES: ====================================================
      INTEGER g, jp
      _RL tdout_bvol(nPlank), tdnet_bvol(nPlank)
      _RL tdout_topt(nPlank), tdnet_topt(nPlank)
      _RL tdout_paropt(nPlank), tdnet_paropt(nPlank)

C==== trait diffusion ==================================================

      DO jp = 1, nPlank
        g = group(jp)
        IF (grp_nbiovol(g) .GT. 1) THEN
          tdout_bvol(jp) = numut_bvol(g)*gcom(jp) /
     &    (LOG(biovolfac(g)))**2
        ENDIF
        IF (grp_ntopt(g) .GT. 1) THEN
          tdout_topt(jp) = numut_topt(g)*gcom(jp) /
     &    (Toptstep(g))**2
        ENDIF
        IF (grp_nparopt(g) .GT. 1) THEN
          tdout_paropt(jp) = numut_paropt(g)*gcom(jp)
     &    / (LOG(PARoptfac(g)))**2
        ENDIF
      ENDDO
      
      DO jp = 1, nplank
        g = group(jp)
        IF (grp_nbiovol(g) .GT. 1) THEN
          IF (igroup(jp) .EQ. 1) THEN
            tdnet_bvol(jp) = tdout_bvol(jp+grp_ntopt(g)*grp_nparopt(g))
     &                     - tdout_bvol(jp)
          ELSEIF (igroup(jp) .EQ. grp_nbiovol(g)) THEN
            tdnet_bvol(jp) = tdout_bvol(jp-grp_ntopt(g)*grp_nparopt(g))
     &                     - tdout_bvol(jp)
          ELSE
            tdnet_bvol(jp) = tdout_bvol(jp-grp_ntopt(g)*grp_nparopt(g))
     &                     + tdout_bvol(jp+grp_ntopt(g)*grp_nparopt(g))
     &                     - 2*tdout_bvol(jp)
          ENDIF
        ELSE
          tdnet_bvol(jp) = 0.0
        ENDIF
        IF (grp_ntopt(g) .GT. 1) THEN
          IF (jgroup(jp) .EQ. 1) THEN
            tdnet_topt(jp) = tdout_topt(jp+grp_nparopt(g))
     &                     - tdout_topt(jp)
          ELSEIF (jgroup(jp) .EQ. grp_ntopt(g)) THEN
            tdnet_topt(jp) = tdout_topt(jp-grp_nparopt(g)) 
     &                     - tdout_topt(jp)
          ELSE
            tdnet_topt(jp) = tdout_topt(jp-grp_nparopt(g)) 
     &      + tdout_topt(jp+grp_nparopt(g)) - 2*tdout_topt(jp)
          ENDIF
        ELSE
          tdnet_topt(jp) = 0.0
        ENDIF
        IF (grp_nparopt(g) .GT. 1) THEN
          IF (kgroup(jp) .EQ. 1) THEN
            tdnet_paropt(jp) = tdout_paropt(jp+1) - tdout_paropt(jp)
          ELSEIF (kgroup(jp) .EQ. grp_nparopt(g)) THEN
            tdnet_paropt(jp) = tdout_paropt(jp-1) - tdout_paropt(jp)
          ELSE
            tdnet_paropt(jp) = tdout_paropt(jp-1) + tdout_paropt(jp+1)
     &                       - 2*tdout_paropt(jp)
          ENDIF
        ELSE
          tdnet_paropt(jp) = 0.0
        ENDIF
        gTr(ic+jp-1)=gTr(ic+jp-1) + tdnet_bvol(jp) + tdnet_topt(jp)
     &                            + tdnet_paropt(jp)
      ENDDO

      RETURN
      END SUBROUTINE
