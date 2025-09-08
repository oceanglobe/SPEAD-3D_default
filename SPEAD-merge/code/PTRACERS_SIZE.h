#ifdef ALLOW_PTRACERS

CBOP
C    !ROUTINE: PTRACERS_SIZE.h
C    !INTERFACE:
C #include PTRACERS_SIZE.h
C Necessary to know if the discrete or the continuous version is activated ?
#include "DARWIN_OPTIONS.h"
 
C    !DESCRIPTION:
C Contains passive tracer array size (number of tracers).

C PTRACERS_num defines how many passive tracers are allocated/exist.
C  and is set here (default 1)
C
C     Number of tracers
      INTEGER PTRACERS_num
#ifdef SPEAD_CONTINUOUS_TRAIT
      PARAMETER(PTRACERS_num = 18 )
#else /* NOT spead_continuous_trait */
      PARAMETER(PTRACERS_num = 351) 
#endif


#ifdef ALLOW_AUTODIFF_TAMC
      INTEGER    iptrkey
      INTEGER    maxpass
      PARAMETER( maxpass     = PTRACERS_num + 2 )
#endif

CEOP
#endif /* ALLOW_PTRACERS */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
