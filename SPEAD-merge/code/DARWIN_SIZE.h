#ifdef ALLOW_DARWIN

CBOP
C    !ROUTINE: DARWIN_SIZE.h
C    !INTERFACE:
C #include DARWIN_SIZE.h
C Necessary to know if the discrete or the continuous version is activated ?
#include "DARWIN_OPTIONS.h"

C    !DESCRIPTION:
C Contains dimensions and index ranges for cell model.

      INTEGER nTrac, nPlank, nPhyp, nGroup, nPhoto, nopt
      INTEGER nTrait

      parameter(nGroup=2)
      parameter(nopt=1)
      parameter(nTrait=3)
#ifdef SPEAD_CONTINUOUS_TRAIT
      parameter(nPlank=2)
      parameter(nPhyp=1)
#ifdef SPEAD_CONTINUOUS_COVARIANCE
      INTEGER nCov
      parameter(nCov=(nTrait-1)*nTrait/2)
#endif
      parameter(nPhoto=nPhyp*(nTrait+1)*(nTrait+2)/2)
      parameter(nTrac=nPhoto+nPlank-nPhyp)
#else /* NOT spead_continuous_trait */
      parameter(nTrac=344) 
      parameter(nPlank=344) 
      parameter(nPhyp=343) 
      parameter(nPhoto=343) 
#endif /* spead_continuous_trait */

#ifndef ALLOW_RADTRANS
      INTEGER nlam
      parameter(nlam=1)
#endif

CEOP
#endif /* ALLOW_DARWIN */
