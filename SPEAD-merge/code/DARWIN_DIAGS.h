#ifdef ALLOW_DARWIN

CBOP
C    !ROUTINE: DARWIN_DIAGS.h
C    !INTERFACE:
C #include DARWIN_DIAGS.h

C    !DESCRIPTION:
C Contains indices into diagnostics array

      integer iPP
      integer iNfix
      integer iDenit
      integer iDenitN
      integer iPPplank
      integer iGRplank
      integer iGrGn
      integer iConsDIN
      integer iConsPO4
      integer iConsSi
      integer iConsFe
C To output terms of the moment derivatives (Le Gland, 16/06/2022)
#if defined SPEAD_DIAG_MD && defined SPEAD_CONTINUOUS_TRAIT
      integer iDmnt
      integer iDvrt
      integer iDcvt
      integer iDvrs
      integer iDcvs
#endif
      integer darwin_nDiag
      PARAMETER(iPP=     1)
      PARAMETER(iNfix=   2)
      PARAMETER(iDenit=  3)
      PARAMETER(iDenitN= 4)
      PARAMETER(iConsPO4=5)
      PARAMETER(iConsSi= 6)
      PARAMETER(iConsFe= 7)
      PARAMETER(iConsDIN=8)
      PARAMETER(iPPplank=9)
#ifdef DARWIN_DIAG_PERTYPE
      PARAMETER(iGRplank=iPPplank+nplank)
      PARAMETER(iGrGn=iGRplank+nplank)
C To output terms of the moment derivatives (Le Gland, 16/06/2022)
#if defined SPEAD_DIAG_MD && defined SPEAD_CONTINUOUS_TRAIT
      PARAMETER(iDmnt=iGrGn+nphyp*nTrait)
      PARAMETER(iDvrt=iDmnt+nphyp*nTrait)
      PARAMETER(iDcvt=iDvrt+nphyp*nTrait)
      PARAMETER(iDvrs=iDcvt+nphyp*nCov)
      PARAMETER(iDcvs=iDvrs+nphyp*nTrait)
      PARAMETER(darwin_nDiag=iDcvs+nphyp*nCov-1)
#else
      PARAMETER(darwin_nDiag=iGrGn+nplank-1)
#endif
#else
      PARAMETER(iGRplank=iPPplank)
      PARAMETER(iGrGn=iGRplank)
      PARAMETER(darwin_nDiag=iGrGn-1)
#endif

CEOP
#endif /* ALLOW_DARWIN */
