module Fasham_Variables_module

use kinds_mod 

implicit none


real(kind=r8b) :: aNO3, aNH4, DON, DET, bact
real(kind=r8b) :: Q1, Q2, Q

real(kind=r8b) :: nitro
real(kind=r8b) :: phyto
real(kind=r8b) :: zoo

real(kind=r8b) :: par0, cloudy, cloud, z, delz, oktas, tday, ratio, trans, albedo, saltd, cosz
real(kind=r8b) :: zenith, alatd, alatr, dec, th, E_0, value, sigma, parz, tmp, beta, tau
real(kind=r8b) :: U2, U1, S, enot, date, temp
real(kind=r8b) :: thour,  dayn


real(kind=r8b) :: prey_growth_rate
real(kind=r8b) :: predator_feeding_rate
real(kind=r8b) :: predator_assim
real(kind=r8b) :: predator_biomass_feeding_rate
real(kind=r8b) :: predator_mortality_rate

!-----------------------------------------------------------------------------------------

real(kind=r8b),parameter ::     alpha   =  0.025D+0    ! initial slope of the P-I curve [(W m-2)-1 d-1]
real(kind=r8b),parameter ::     aK1     =  0.5D+0      ! half-saturation for phytoplankton NO3 uptake [mMol N m-3]
real(kind=r8b),parameter ::     aK2     =  0.5D+0      ! half-saturation for phytoplankton NH4 uptake [mMol N m-3]
real(kind=r8b),parameter ::     amu1    =  0.045D+0    ! phytoplankton specific mortality rate [d-1]
real(kind=r8b),parameter ::     akc     =  0.03D+0     ! light attenuation by phytoplankton [m^2 mMol N)-1]
real(kind=r8b),parameter ::     gamma1  =  0.05D+0     ! fraction of total primary production that is exuded [n.d.]
real(kind=r8b),parameter ::     phi     =  1.5D+0      ! phytoplankton ammonium inhibition parameter [(mMol N)-1]
real(kind=r8b),parameter ::     g       =  1.0D+0      ! maximum zooplankton growth rate [d-1]
real(kind=r8b),parameter ::     beta1   =  0.75D+0     ! zooplankton assimilation efficiency of zooplankton [n.d.]
real(kind=r8b),parameter ::     beta2   =  0.75D+0     ! zooplankton assimilation efficiency of phytoplankton [n.d.]
real(kind=r8b),parameter ::     beta3   =  0.75D+0     ! zooplankton assimilation efficiency of bacteria [n.d.]
real(kind=r8b),parameter ::     amu2    =  0.1D+0      ! zooplankton specific excretion rate [d-1]
real(kind=r8b),parameter ::     amu5    =  0.05D+0     ! zooplankton specific mortality rate [d-1]
real(kind=r8b),parameter ::     aK3     =  1.0D+0      ! zooplankton half-saturation conts. for ingestion [d-1]
real(kind=r8b),parameter ::     omega   =  0.33D+0     ! detrital fraction of zooplankton mortality [n.d.]
real(kind=r8b),parameter ::     epsilon =  0.75D+0     ! ammonium fraction of zooplankton excretion [n.d.]
real(kind=r8b),parameter ::     Vb      =  2.0D+0      ! bacteria maximum growth rate [d-1]
real(kind=r8b),parameter ::     Vp      =  2.9D+0      ! phyto maximum growth rate [d-1]
real(kind=r8b),parameter ::     amu3    =  0.05D+0     ! bacteria specific excretion rate [d-1]
real(kind=r8b),parameter ::     aK4     =  0.5D+0      ! bacteria half-saturation rate for uptake [(mMol N) m-3]
real(kind=r8b),parameter ::     eta     =  0.6D+0      ! ammonium/DON uptake ratio [n.d.]
real(kind=r8b),parameter ::     amu4    =  0.05D+0     ! detrital breakdown rate [d-1]
real(kind=r8b),parameter ::     V       =  1.0D+0      ! detrital sinking rate [m d-1]
real(kind=r8b),parameter ::     p1      =  1.0D+0      ! zooplankton preference for phytoplankton [n.d.]
real(kind=r8b),parameter ::     p2      =  1.0D+0      ! zooplankton preference for bacteria [n.d.]
real(kind=r8b),parameter ::     p3      =  1.0D+0      ! zooplankton preference for detritus [n.d.]
real(kind=r8b),parameter ::     aN0     =  2.0D+0      ! concentration of NO3 below the mixed-layer [(mMol N) m-3]


!-----------------------------------------------------------------------------------------

! parameters as in Table 1; Fasham et al. [JMR, 48, 591-639, 1990]

real(kind=r8b),parameter ::     akw = 0.04D+0     ! light attenuation due to sea water [m-1]
real(kind=r8b),parameter ::     am  = 0.1D+0      ! cross-thermocline mixing rate [m d-1]


!-----------------------------------------------------------------------------------------

! Variables to define indexes of species in b_tmp and in Numerical Code Solutions

integer(kind=i4b),parameter :: SPECIES_NITRATE                      = -1
integer(kind=i4b),parameter :: SPECIES_AMMONIUM                     = -2
integer(kind=i4b),parameter :: SPECIES_DISSOLVED_ORGANIC_NITROGEN   = -3
integer(kind=i4b),parameter :: SPECIES_DETRITUS                     = -4
integer(kind=i4b),parameter :: SPECIES_BACTERIA                     = -5
integer(kind=i4b),parameter :: SPECIES_PHYTOPLANKTON                = -6
integer(kind=i4b),parameter :: SPECIES_ZOOPLANKTON                  = -7

integer(kind=i4b) :: SPECIES_Zoo
integer(kind=i4b) :: SPECIES_Phyto
integer(kind=i4b) :: SPECIES_Nitro



! Variables to define indexes of forcing functions in Numerical_CODE_Forcing_Functions

integer(kind=i4b),parameter :: FORCING_MIXED_LAYER_DEPTH         = -5001
integer(kind=i4b),parameter :: FORCING_MLD_CHANGE_NON_MOTILE     = -5002
integer(kind=i4b),parameter :: FORCING_MLD_CHANGE_MOTILE         = -5003
integer(kind=i4b),parameter :: FORCING_LIGHT_LIMITED_GROWTH_RATE = -5004



end module Fasham_Variables_module
