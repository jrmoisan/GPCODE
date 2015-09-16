!> @brief
!>  This module declares variables and parameters used in the Fasham model and 
!!  sets a value in some of the parameters.
!>
!> @details
!>  This module declares variables and parameters used in the Fasham model and 
!!  sets a value in some of the parameters.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

MODULE Fasham_Variables_module

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!

!---------------------------------------------------------------------------  

USE kinds_mod 

IMPLICIT none


REAL (KIND=r8b) :: aNO3, aNH4, DON, DET, bact
REAL (KIND=r8b) :: Q1, Q2, Q

REAL (KIND=r8b) :: nitro
REAL (KIND=r8b) :: phyto
REAL (KIND=r8b) :: zoo

REAL (KIND=r8b) :: par0, cloudy, cloud, z, delz, oktas, tday, ratio, trans, albedo, saltd, cosz
REAL (KIND=r8b) :: zenith, alatd, alatr, dec, th, E_0, value, sigma, parz, tmp, beta, tau
REAL (KIND=r8b) :: U2, U1, S, enot, date, temp
REAL (KIND=r8b) :: thour,  dayn


REAL (KIND=r8b) :: prey_growth_rate
REAL (KIND=r8b) :: predator_feeding_rate
REAL (KIND=r8b) :: predator_assim
REAL (KIND=r8b) :: predator_biomass_feeding_rate
REAL (KIND=r8b) :: predator_mortality_rate

!-----------------------------------------------------------------------------------------

REAL (KIND=r8b),PARAMETER ::     alpha   =  0.025D+0    ! initial slope of the P-I curve [(W m-2)-1 d-1]
REAL (KIND=r8b),PARAMETER ::     aK1     =  0.5D+0      ! half-saturation for phytoplankton NO3 uptake [mMol N m-3]
REAL (KIND=r8b),PARAMETER ::     aK2     =  0.5D+0      ! half-saturation for phytoplankton NH4 uptake [mMol N m-3]
REAL (KIND=r8b),PARAMETER ::     amu1    =  0.045D+0    ! phytoplankton specific mortality rate [d-1]
REAL (KIND=r8b),PARAMETER ::     akc     =  0.03D+0     ! light attenuation by phytoplankton [m^2 mMol N)-1]
REAL (KIND=r8b),PARAMETER ::     gamma1  =  0.05D+0     ! fraction of total primary production that is exuded [n.d.]
REAL (KIND=r8b),PARAMETER ::     phi     =  1.5D+0      ! phytoplankton ammonium inhibition parameter [(mMol N)-1]
REAL (KIND=r8b),PARAMETER ::     g       =  1.0D+0      ! maximum zooplankton growth rate [d-1]
REAL (KIND=r8b),PARAMETER ::     beta1   =  0.75D+0     ! zooplankton assimilation efficiency of zooplankton [n.d.]
REAL (KIND=r8b),PARAMETER ::     beta2   =  0.75D+0     ! zooplankton assimilation efficiency of phytoplankton [n.d.]
REAL (KIND=r8b),PARAMETER ::     beta3   =  0.75D+0     ! zooplankton assimilation efficiency of bacteria [n.d.]
REAL (KIND=r8b),PARAMETER ::     amu2    =  0.1D+0      ! zooplankton specific excretion rate [d-1]
REAL (KIND=r8b),PARAMETER ::     amu5    =  0.05D+0     ! zooplankton specific mortality rate [d-1]
REAL (KIND=r8b),PARAMETER ::     aK3     =  1.0D+0      ! zooplankton half-saturation conts. for ingestion [d-1]
REAL (KIND=r8b),PARAMETER ::     omega   =  0.33D+0     ! detrital fraction of zooplankton mortality [n.d.]
REAL (KIND=r8b),PARAMETER ::     epsilon =  0.75D+0     ! ammonium fraction of zooplankton excretion [n.d.]
REAL (KIND=r8b),PARAMETER ::     Vb      =  2.0D+0      ! bacteria maximum growth rate [d-1]
REAL (KIND=r8b),PARAMETER ::     Vp      =  2.9D+0      ! phyto maximum growth rate [d-1]
REAL (KIND=r8b),PARAMETER ::     amu3    =  0.05D+0     ! bacteria specific excretion rate [d-1]
REAL (KIND=r8b),PARAMETER ::     aK4     =  0.5D+0      ! bacteria half-saturation rate for uptake [(mMol N) m-3]
REAL (KIND=r8b),PARAMETER ::     eta     =  0.6D+0      ! ammonium/DON uptake ratio [n.d.]
REAL (KIND=r8b),PARAMETER ::     amu4    =  0.05D+0     ! detrital breakdown rate [d-1]
REAL (KIND=r8b),PARAMETER ::     V       =  1.0D+0      ! detrital sinking rate [m d-1]
REAL (KIND=r8b),PARAMETER ::     p1      =  1.0D+0      ! zooplankton preference for phytoplankton [n.d.]
REAL (KIND=r8b),PARAMETER ::     p2      =  1.0D+0      ! zooplankton preference for bacteria [n.d.]
REAL (KIND=r8b),PARAMETER ::     p3      =  1.0D+0      ! zooplankton preference for detritus [n.d.]
REAL (KIND=r8b),PARAMETER ::     aN0     =  2.0D+0      ! concentration of NO3 below the mixed-layer [(mMol N) m-3]


!-----------------------------------------------------------------------------------------

! parameters as in Table 1; Fasham et al. [JMR, 48, 591-639, 1990]

REAL (KIND=r8b),PARAMETER ::     akw = 0.04D+0     ! light attenuation due to sea water [m-1]
REAL (KIND=r8b),PARAMETER ::     am  = 0.1D+0      ! cross-thermocline mixing rate [m d-1]


!-----------------------------------------------------------------------------------------

! Variables to define indexes of species in b_tmp and in Numerical Code Solutions

INTEGER (KIND=i4b),PARAMETER :: SPECIES_NITRATE                      = -1
INTEGER (KIND=i4b),PARAMETER :: SPECIES_AMMONIUM                     = -2
INTEGER (KIND=i4b),PARAMETER :: SPECIES_DISSOLVED_ORGANIC_NITROGEN   = -3
INTEGER (KIND=i4b),PARAMETER :: SPECIES_DETRITUS                     = -4
INTEGER (KIND=i4b),PARAMETER :: SPECIES_BACTERIA                     = -5
INTEGER (KIND=i4b),PARAMETER :: SPECIES_PHYTOPLANKTON                = -6
INTEGER (KIND=i4b),PARAMETER :: SPECIES_ZOOPLANKTON                  = -7

INTEGER (KIND=i4b) :: SPECIES_Zoo
INTEGER (KIND=i4b) :: SPECIES_Phyto
INTEGER (KIND=i4b) :: SPECIES_Nitro



! Variables to define indexes of forcing functions in Numerical_CODE_Forcing_Functions

INTEGER (KIND=i4b),PARAMETER :: FORCING_MIXED_LAYER_DEPTH         = -5001
INTEGER (KIND=i4b),PARAMETER :: FORCING_MLD_CHANGE_NON_MOTILE     = -5002
INTEGER (KIND=i4b),PARAMETER :: FORCING_MLD_CHANGE_MOTILE         = -5003
INTEGER (KIND=i4b),PARAMETER :: FORCING_LIGHT_LIMITED_GROWTH_RATE = -5004



END MODULE Fasham_Variables_module
