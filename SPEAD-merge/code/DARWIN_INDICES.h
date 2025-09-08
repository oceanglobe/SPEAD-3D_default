#ifdef ALLOW_DARWIN

CBOP
C    !ROUTINE: DARWIN_INDICES.h
C    !INTERFACE:
C #include DARWIN_INDICES.h

C    !DESCRIPTION:
C Contains indices into ptracer array

      INTEGER eDIC
      INTEGER iNO3
      INTEGER iNO2
      INTEGER iNH4
      INTEGER iPO4
      INTEGER iFeT
      INTEGER eSiO2
      INTEGER iDOC
      INTEGER eDOFe
      INTEGER iPOC
      INTEGER ePOSi
      INTEGER ePIC
      INTEGER ic
      INTEGER eCARBON
      INTEGER eCDOM
      INTEGER ec
      INTEGER en
      INTEGER ep
      INTEGER efe
      INTEGER esi
      INTEGER eChl
      INTEGER nDarwin
#ifdef SPEAD_DIC
       INTEGER iDIC
       PARAMETER (iDIC   =1)
       PARAMETER (eDIC   =iDIC)
#else
       PARAMETER (eDIC   =0)
#endif
      PARAMETER (iNO3   =eDIC+1)
      PARAMETER (iNO2   =iNO3 +1)
      PARAMETER (iNH4   =iNO2 +1)
      PARAMETER (iPO4   =iNH4 +1)
      PARAMETER (iFeT   =iPO4 +1)
#ifdef SPEAD_SILICA
       INTEGER iSiO2
       PARAMETER (iSiO2  =iFeT +1)
       PARAMETER (eSiO2  =iSiO2)
#else
       PARAMETER (eSiO2  =iFeT)
#endif
      PARAMETER (iDOC   =eSiO2+1)
#ifdef SPEAD_STOICHIOMETRY_OM
       INTEGER iDON, iDOP, iDOFe
       PARAMETER (iDON   =iDOC +1)
       PARAMETER (iDOP   =iDON +1)
       PARAMETER (iDOFe  =iDOP +1)
       PARAMETER (eDOFe  =iDOFe)
#else
       PARAMETER (eDOFe  =iDOC)
#endif
      PARAMETER (iPOC   =eDOFe+1)
#ifdef SPEAD_STOICHIOMETRY_OM
       INTEGER iPON, iPOP, iPOFe
       PARAMETER (iPON   =iPOC +1)
       PARAMETER (iPOP   =iPON +1)
       PARAMETER (iPOFe  =iPOP +1)
#ifdef SPEAD_SILICA
        INTEGER iPOSi
        PARAMETER (iPOSi  =iPOFe+1)
        PARAMETER (ePOSi  =iPOSi)
#else
        PARAMETER (ePOSi  =iPOFe)
#endif /* SPEAD_SILICA */
#else
       PARAMETER (ePOSi  =iPOC)
#endif /* SPEAD_STOICHIOMETRY_OM */
#ifdef SPEAD_DIC
       INTEGER iPIC
       PARAMETER (iPIC   =ePOSi+1)
       PARAMETER (ePIC   =iPIC)
#else
       PARAMETER(ePIC = ePOSi)
#endif
#ifdef DARWIN_ALLOW_CARBON
       INTEGER iALK
       INTEGER iO2
       PARAMETER (iALK   =ePIC +1)
       PARAMETER (iO2    =iALK +1)
       PARAMETER (eCARBON=iO2)
#else
       PARAMETER (eCARBON=ePIC)
#endif
#ifdef DARWIN_ALLOW_CDOM
       INTEGER iCDOM
       PARAMETER (iCDOM  =eCARBON+1)
       PARAMETER (eCDOM  =iCDOM)
#else
       PARAMETER (eCDOM  =eCARBON)
#endif
       PARAMETER (ic     =eCDOM+1)
       PARAMETER (ec     =ic   +nTrac-1)
#ifdef DARWIN_ALLOW_NQUOTA
       INTEGER in
       PARAMETER (in     =ec  +1)
       PARAMETER (en     =in  +nTrac-1)
#else
       PARAMETER (en     =ec)
#endif
#ifdef DARWIN_ALLOW_PQUOTA
       INTEGER ip
       PARAMETER (ip     =en  +1)
       PARAMETER (ep     =ip  +nTrac-1)
#else
       PARAMETER (ep     =en)
#endif
#ifdef DARWIN_ALLOW_FEQUOTA
       INTEGER ife
       PARAMETER (ife    =ep  +1)
       PARAMETER (efe    =ife +nTrac-1)
#else
       PARAMETER (efe    =ep)
#endif
#ifdef DARWIN_ALLOW_SIQUOTA
       INTEGER isi
       PARAMETER (isi    =efe +1)
       PARAMETER (esi    =isi +nTrac-1)
#else
       PARAMETER (esi    =efe)
#endif
#ifdef DARWIN_ALLOW_CHLQUOTA
       INTEGER iChl
       PARAMETER (iChl   =esi +1)
       PARAMETER (eChl   =iChl+nPhoto-1)
#else
       PARAMETER (eChl   =efe)
#endif
       PARAMETER (nDarwin=eChl)

CEOP
#endif /* ALLOW_DARWIN */
