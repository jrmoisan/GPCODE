!> @brief
!>  Computes some forcing functions for the Fasham model.

!> @details
!>  Computes some forcing functions for the Fasham model.

!> @author  Dave Coulter
!> @date June 26, 2013     Dave Coulter        

!> @param[in] species - current values of the Fasham variables
!> @param[in] day
!> @param[in] aMLD
!> @param[out] aJ
!> @param[out] L_bad  - true if error computing aJ

SUBROUTINE JQforce(species, day, aMLD, aJ, L_bad)
!---------------------------------------------------------------------------  
!
! File:   Fasham_Forcing.f90
! Author: Dave Coulter
!
! Created on June 26, 2013, 1:37 PM
!
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

USE kinds_mod

USE fasham_variables_module
USE GP_variables_module

    IMPLICIT none

    REAL (KIND=r8b), INTENT(IN) :: species(n_CODE_Equations), day, aMLD
    REAL (KIND=r8b), INTENT(OUT) :: aJ
    REAL (KIND=r8b) :: daym(14),coktas(14)

    REAL (KIND=r8b),PARAMETER :: solar=1353.D+0 ! the solar max is from Brock, 1981
    REAL (KIND=r8b),PARAMETER :: pi = 3.141592653589793D0

    INTEGER (KIND=i4b) :: iz,i

    LOGICAL :: L_bad

    DATA daym /0.D+0,16.D+0,46.D+0,75.D+0,105.D+0,136.D+0,166.D+0,&
               197.D+0,228.D+0,258.D+0,289.D+0, 319.D+0, 350.D+0, 365.D+0/

    ! OWS Papa CLOUD CLIMATOLOGY
    !offdata coktas /6.29D+0,6.26D+0,6.31D+0,6.31D+0,6.32D+0,6.70D+0,&
    !           7.12D+0,7.26D+0,6.93D+0,6.25D+0,6.19D+0,6.23D+0,6.31D+0,6.29D+0/

    ! BATS CLOUD CLIMATOLOGY
    DATA coktas /4.00D+0,4.00D+0,4.00D+0,4.00D+0,4.00D+0,4.00D+0,4.00D+0,&
                 4.00D+0,4.00D+0,4.00D+0,4.00D+0,4.00D+0,4.00D+0,4.00D+0/

    ! BATS COADS CLOUD CLIMATOLOGY
    !off data coktas /5.25D+0,5.36D+0,5.38D+0,5.18D+0,4.84D+0,4.65D+0,&
    !                 4.58D+0,3.87D+0,3.81D+0,4.15D+0,4.70D+0,4.81D+0,5.14D+0,5.25D+0/


    !---------------------------------------------------------------------------------
 
    L_bad = .FALSE.

    ! Copy phytoplankton to phyto

    phyto = species(MAX (1,ABS (SPECIES_PHYTOPLANKTON)))


    dec=-0.406D+0*COS (2.D+0*pi*day/365.D+0) ! in radians and from Evans and Parlsow, 1985
    alatr=alatd*2.D+0*pi/360.D+0


!   correction of solar constant due to the ellipticity of the earth orbit

    enot=1.D+0+(0.033D+0*COS (2.D+0*pi*day/365.D+0))   ! E_0 from [Duffie and Beckman, 1980]


    th=(ABS (thour-12.D+0)*15.D+0)*(2.D+0*pi/360.D+0)     ! hour angle [15.D+0 ==> degrees per hour]
    cosz=((SIN (dec)*SIN (alatr))+(COS (dec)*COS (alatr)*COS (th)))
    zenith=ACOS (cosz)  ! zenith ==> solar zenith angle in radians
    saltd=90.D+0-((360.D+0*zenith)/(2.D+0*pi))  ! solar altitude in degrees
    albedo=0.04D+0
    trans=1.0D+0

    DO  i=1,13
        IF ( day .ge. daym(i) .and. day .le. daym(i+1) ) THEN
            ratio=(day-daym(i))/(daym(i+1)-daym(i))
            oktas=coktas(i)+(ratio*(coktas(i+1)-coktas(i)))
        END IF
    END DO ! i

    cloud=oktas/8.D+0


!   cloud corr from Smith and Dobson, [Ecological Modeling, 14, 1-19, 1981]


    IF ( ABS (zenith) .le. 1.55D+0 ) THEN
        cloudy=0.0375D+0+(cosz*EXP (-0.24D+0/cosz)*((cloud*EXP (-0.07D+0/cosz))+1.D+0-cloud))
    ELSE
        cloudy=0.0375D+0
    END IF



!   calculate the light value only during the day [night = zenith .ge. pi/2.D+0]

    IF ( ABS (zenith) .lt. pi/2.D+0 ) THEN
        par0=solar*trans*enot*cosz*0.43D+0*(1.D+0-albedo)*cloudy
    ELSE
        par0=0.D+0
    END IF


!   calculate aJ by integrating (simple trapezoidal) PvsI terms over the aMLD


    aJ=0.D+0
    IF ( par0 .gt. 0.D+0) THEN

        delz=aMLD/100.D+0
        DO  iz=0,100
            z=FLOAT (iz)*delz



            IF ( ABS ( z * ( akw +  akc*ABS (phyto) ) ) < 150.0d0 ) THEN 
                parz = par0 * EXP ( -z * ( akw +  akc*ABS (phyto) ) )
            ELSE
                parz = 1.0D-65
            END IF 


            IF ( parz > 1.0D+100 ) THEN
                L_bad = .true. 
                RETURN
            END IF ! parz > 1.0D+100


            IF ( Vp**2 + (alpha * parz)**2 > 0.0d0 ) THEN

                tmp = (Vp*alpha*parz) / ( SQRT ( Vp**2 + (alpha * parz)**2 ) )
            ELSE

                tmp = 0.0d0

            END IF ! Vp**2 + (alpha * parz)**2 > 0.0d0



            IF ( ISNAN ( tmp ) ) THEN
                L_bad = .true. 
                RETURN
            END IF ! ISNAN ( tmp ) 


            IF ( iz .eq. 0 .or. iz .eq. 100 ) THEN
                aJ=aJ+(0.5D+0*tmp)
            ELSE IF ( iz .gt. 0 .and. iz .lt. 100 ) THEN
                aJ=aJ+tmp
            END IF
        END DO


        aJ=aJ*delz/aMLD

    END IF !  par0 .gt. 0.D+0


END SUBROUTINE JQforce



SUBROUTINE mldforce(day, h, aMLD, L_bad )

 
!---------------------------------------------------------------------------  
!> @author  Dave Coulter
!> @date June 26, 2013     Dave  Coulter       

!> @brief
!>  Computes the aMLD and h(t) terms of forcing functions for the Fasham model.

!> @details
!>  Computes the aMLD and h(t) terms of forcing functions for the Fasham model.

!> @param[in] day
!>
!> @param[out] h 
!> @param[out] aMLD
!> @param[out] L_bad = true if error in computation

!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 

!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! subroutine to determine the aMLD and h(t) terms
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

USE fasham_variables_module
USE GP_variables_module

    IMPLICIT none

    REAL (KIND=r8b), INTENT(IN) :: day
    REAL (KIND=r8b), INTENT(OUT) :: h, aMLD
    INTEGER (KIND=i4b), parameter :: n=14
    REAL (KIND=r8b) :: daym(14),cmld(14)
    INTEGER (KIND=i4b) nloop,iloop,i

    DATA daym /0.D+0,16.D+0,46.D+0,75.D+0,105.D+0,136.D+0,&
               166.D+0,197.D+0,228.D+0,258.D+0,289.D+0,319.D+0,350.D+0,365.D+0/

    ! forcing data for SUPER [a.k.a. OWS Papa] site
    ! forcing data for OWS Papa site [data digitized from Matear, 1995]

    DATA cmld /87.0D+0,92.5D+0,100.0D+0,93.5D+0,67.0D+0,45.0D+0,&
               30.0D+0,19.0D+0,22.5D+0,30.0D+0,45.0D+0,65.0D+0,82.0D+0,87.0D+0/

    LOGICAL :: L_bad 

!----------------------------------------------------------------------------------

    L_bad = .FALSE.

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    ! calculate the aMLD and h values
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    ! find a [nloop+1] day moving average for h

    nloop=30     ! make the number EVEN to keep it symmetrical
    th=0.D+0


    DO  iloop=1,nloop+1

        tday=day-(nloop/2.D+0)+iloop

        IF ( tday .lt.   0.D+0 ) tday=tday+365.D+0
        IF ( tday .ge. 365.D+0 ) tday=tday-365.D+0

        DO  i=1,n-1
            IF ( tday .ge. daym(i) .and. tday .lt. daym(i+1) ) THEN
                th=th+((cmld(i+1)-cmld(i))/(daym(i+1)-daym(i)))
            END IF
        END DO ! i 

    END DO ! iloop

    h=th/(nloop+1)


    DO  i=1,n-1
        IF ( day .ge. daym(i) .and. day .le. daym(i+1) ) THEN
            ratio=(day-daym(i))/(daym(i+1)-daym(i))
            aMLD=cmld(i)+(ratio*(cmld(i+1)-cmld(i)))
        END IF
    END DO ! i

    IF ( ISNAN ( h ) .or. ABS (h) > 1.0D100 ) THEN
        L_bad = .true. 
        RETURN
    END IF ! ISNAN ( h ) 

    IF ( ISNAN ( aMLD ) .or. ABS (aMLD) > 1.0D100 ) THEN
        L_bad = .true. 
        RETURN
    END IF ! ISNAN ( aMLD ) 


    RETURN

END SUBROUTINE mldforce
