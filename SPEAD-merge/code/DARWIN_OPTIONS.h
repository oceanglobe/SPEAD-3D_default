#ifndef DARWIN_OPTIONS_H
#define DARWIN_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#ifdef ALLOW_DARWIN

#include "CPP_OPTIONS.h"

CBOP
C    !ROUTINE: DARWIN_OPTIONS.h
C    !INTERFACE:

C    !DESCRIPTION:
C options for darwin package
CEOP

C tracer selection

C enable (or disable) nitrogen quotas for all plankton
#undef  DARWIN_ALLOW_NQUOTA

C enable (or disable) phosphorus quotas for all plankton
#undef  DARWIN_ALLOW_PQUOTA

C enable (or disable) iron quotas for all plankton
#undef  DARWIN_ALLOW_FEQUOTA

C enable (or disable) silica quotas for all plankton
#undef  DARWIN_ALLOW_SIQUOTA

C enable (or disable) chlorophyll quotas for all phototrophs
#undef  DARWIN_ALLOW_CHLQUOTA

C enable (or disable) a dynamic CDOM tracer
#undef  DARWIN_ALLOW_CDOM

C enable (or disable) air-sea carbon exchange and Alk and O2 tracers
#undef  DARWIN_ALLOW_CARBON

C optional bits

C enable (or disable) denitrification code
#undef  DARWIN_ALLOW_DENIT

C enable (or disable) separate exudation of individual elements
#undef  DARWIN_ALLOW_EXUDE

C enable (or disable) old virtualflux code for DIC and Alk
#undef  ALLOW_OLD_VIRTUALFLUX

C reduce nitrate uptake by iron limitation factor
#undef DARWIN_NITRATE_FELIMIT

C allow organic matter to sink into bottom (sedimentize)
C#define DARWIN_BOTTOM_SINK
#undef DARWIN_BOTTOM_SINK


C light

C compute average PAR in layer, assuming exponential decay
C (ignored when radtrans package is used)
#undef  DARWIN_AVPAR

C enable (or disable) GEIDER light code
#undef  DARWIN_ALLOW_GEIDER

C use rho instead of acclimated Chl:C for chlorophyll synthesis
#undef  DARWIN_GEIDER_RHO_SYNTH

C initialize chl as in darwin2 (with radtrans package)
#undef  DARWIN_CHL_INIT_LEGACY

C scattering coefficients are per Chlorophyll (with radtrans package)
#undef  DARWIN_SCATTER_CHL

C make diagnostics for instrinsic optical properties available
#define DARWIN_DIAG_IOP


C grazing

C for quadratic grazing as in darwin2+quota
C 1 represents selective grazing on first trait only (Le Gland et al, 19/04/2024)
#define DARWIN_GRAZING_SWITCH
#define DARWIN_GRAZING_SWITCH_1trait

C compute palat from size ratios
#undef  DARWIN_ALLOMETRIC_PALAT

C turn off grazing temperature dependence
#undef  DARWIN_NOZOOTEMP
#undef  DARWIN_TIME_GRAZING

C temperature

C turn off all temperature dependence
#undef  DARWIN_NOTEMP

C select temperature version: 1, 2 or 3
#define DARWIN_TEMP_VERSION 6

C restrict phytoplankton growth to a temperature range
#define DARWIN_TEMP_RANGE


C iron

C restrict maximum free iron
#define DARWIN_MINFE

C enable particle scavenging code
#undef  DARWIN_PART_SCAV

C enable variable iron sediment source
#undef  DARWIN_IRON_SED_SOURCE_VARIABLE

C revert to old variable iron sediment source in terms of POP
#undef  DARWIN_IRON_SED_SOURCE_POP

C diagnostics

C include code for per-type diagnostics
#define  DARWIN_DIAG_PERTYPE


C debugging

C turn on debugging output
#undef  DARWIN_DEBUG

C compute and print global element totals
#define DARWIN_ALLOW_CONS

C value for unused traits
#define DARWIN_UNUSED 0


C deprecated

C base particle scavenging on POP as in darwin2
#undef  DARWIN_PART_SCAV_POP


C random trait generation
C these are for darwin_generate_random
C assign traits based on random numbers as in darwin2
#undef  DARWIN_RANDOM_TRAITS

C set traits for darwin2 2-species setup (requires DARWIN_RANDOM_TRAITS)
#undef  DARWIN_TWO_SPECIES_SETUP

C set traits for darwin2 9-species setup (requires DARWIN_RANDOM_TRAITS)
#undef  DARWIN_NINE_SPECIES_SETUP

C enable diazotrophy when using (requires DARWIN_RANDOM_TRAITS)
#undef  DARWIN_ALLOW_DIAZ

C If defined, the continuous approach to trait diversity is used
#define SPEAD_CONTINUOUS_TRAIT
C If defined, ntrait must be superior or equal to 2
C If undefined, traits are assumed to be independent (cov = 0)
#define SPEAD_CONTINUOUS_COVARIANCE

C If defined, output moment derivatives
#define SPEAD_DIAG_MD

C If defined, derivatives with respect to traits are computed locally
C If undefined, differences over the whole trait distribution are used
#undef  SPEAD_DERIVATIVES_LOCAL

C If undefined, DIC and PIC is not modeled
#undef  SPEAD_DIC

C If undefined, silica and diatoms are not modeled
#undef  SPEAD_SILICA

C If undefined, all DOM and POM have the same stoichiometry
C and only DOC and POC are simulated
#undef  SPEAD_STOICHIOMETRY_OM 
  
C If defined the phytoplankton is a Darwian Demon
#undef  DARWINIAN_DEMON_SIZE
#undef  DARWINIAN_DEMON_TEMP
#undef  DARWINIAN_DEMON_LIGHT

C Better to switch TD off in discrete mode if not used
#define SPEAD_TRAIT_DIFFUSION

#endif /* ALLOW_DARWIN */
#endif /* DARWIN_OPTIONS_H */
