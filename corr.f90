SUBROUTINE CORR( X,Y,N,IWRITE,C , &
                 dt, sse_min_time, sse_max_time, sse_low_wt  ) 
!                                                                       
!     PURPOSE--THIS SUBROUTINE COMPUTES THE                             
!              SAMPLE CORRELATION COEFFICIENT                           
!              BETWEEN THE 2 SETS OF DATA IN THE INPUT VECTORS X AND Y. 
!              THE SAMPLE CORRELATION COEFFICIENT WILL BE A SINGLE      
!              PRECISION VALUE BETWEEN -1.0 AND 1.0 (INCLUSIVELY).      
!     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF         
!                                (UNSORTED) OBSERVATIONS                
!                                WHICH CONSTITUTE THE FIRST SET         
!                                OF DATA.                               
!                     --Y      = THE SINGLE PRECISION VECTOR OF         
!                                (UNSORTED) OBSERVATIONS                
!                                WHICH CONSTITUTE THE SECOND SET        
!                                OF DATA.                               
!                     --N      = THE INTEGER NUMBER OF OBSERVATIONS     
!                                IN THE VECTOR X, OR EQUIVALENTLY,      
!                                THE INTEGER NUMBER OF OBSERVATIONS     
!                                IN THE VECTOR Y.                       
!                     --IWRITE = AN INTEGER FLAG CODE WHICH             
!                                (IF SET TO 0) WILL SUPPRESS            
!                                THE PRINTING OF THE                    
!                                SAMPLE CORRELATION COEFFICIENT         
!                                AS IT IS COMPUTED;                     
!                                OR (IF SET TO SOME INTEGER             
!                                VALUE NOT EQUAL TO 0),                 
!                                LIKE, SAY, 1) WILL CAUSE               
!                                THE PRINTING OF THE                    
!                                SAMPLE CORRELATION COEFFICIENT         
!                                AT THE TIME IT IS COMPUTED.            
!     OUTPUT ARGUMENTS--C      = THE SINGLE PRECISION VALUE OF THE      
!                                COMPUTED SAMPLE CORRELATION COEFFICIENT
!                                BETWEEN THE 2 SETS OF DATA             
!                                IN THE INPUT VECTORS X AND Y.          
!                                THIS SINGLE PRECISION VALUE            
!                                WILL BE BETWEEN -1.0 AND 1.0           
!                                (INCLUSIVELY).                         
!     OUTPUT--THE COMPUTED SINGLE PRECISION VALUE OF THE                
!             SAMPLE CORRELATION COEFFICIENT BETWEEN THE 2 SETS         
!             OF DATA IN THE INPUT VECTORS X AND Y.                     
!     PRINTING--NONE, UNLESS IWRITE HAS BEEN SET TO A NON-ZERO          
!               INTEGER, OR UNLESS AN INPUT ARGUMENT ERROR              
!               CONDITION EXISTS.                                       
!     RESTRICTIONS--THERE IS NO RESTRICTION ON THE MAXIMUM VALUE        
!                   OF N FOR THIS SUBROUTINE.                           
!     OTHER DATAPAC   SUBROUTINES NEEDED--NONE.                         
!     FORTRAN LIBRARY SUBROUTINES NEEDED--SQRT.                         
!     MODE OF INTERNAL OPERATIONS--SINGLE PRECISION.                    
!     LANGUAGE--ANSI FORTRAN.                                           
!     REFERENCES--KENDALL AND STUART, THE ADVANCED THEORY OF            
!                 STATISTICS, VOLUME 1, EDITION 2, 1963, PAGES 235-236. 
!               --KENDALL AND STUART, THE ADVANCED THEORY OF            
!                 STATISTICS, VOLUME 2, EDITION 1, 1961, PAGES 292-293. 
!               --SNEDECOR AND COCHRAN, STATISTICAL METHODS,            
!                 EDITION 6, 1967, PAGES 172-198.                       
!     WRITTEN BY--JAMES J. FILLIBEN                                     
!                 STATISTICAL ENGINEERING LABORATORY (205.03)           
!                 NATIONAL BUREAU OF STANDARDS                          
!                 WASHINGTON, D. C. 20234                               
!                 PHONE:  301-921-2315                                  
!     ORIGINAL VERSION--JUNE      1972.                                 
!     UPDATED         --SEPTEMBER 1975.                                 
!     UPDATED         --NOVEMBER  1975.                                 
!                                                                       
!---------------------------------------------------------------------  

use kinds_mod 

implicit none

integer, intent(in) :: n
integer, intent(in) :: iwrite

real(kind=r8b),intent(in)    ::  dt
real(kind=r8b),intent(in)    ::  sse_min_time
real(kind=r8b),intent(in)    ::  sse_max_time
real(kind=r8b),intent(in)    ::  sse_low_wt

real(kind=r8b)  ::  xi  
real(kind=r8b)  ::  sse_wt


real(kind=r8b), dimension(n) ::     X
real(kind=r8b), dimension(n) ::     Y

real(kind=r8b)               ::     C

real(kind=r8b)  :: AN 
real(kind=r8b)  ::  HOLD
real(kind=r8b)  ::  XBAR
real(kind=r8b)  ::  YBAR
!                                                                       
real(kind=r8b)  ::  SUM1
real(kind=r8b)  ::  SUM2
real(kind=r8b)  ::  SUM3
integer(kind=i4b) :: ipr    
integer(kind=i4b) :: IFLAG
integer(kind=i4b) :: I

!---------------------------------------------------------------------  

      IPR=6 

!     CHECK THE INPUT ARGUMENTS FOR ERRORS                              

      AN = N 

      C=0.0d0 

      IFLAG=0 
      IF( N .LT. 1 )then
          WRITE(IPR,25) 
          write(ipr,'(A)') 'The third argument is the number of points = dimension of the input array'
          WRITE(IPR,47)N 
          RETURN 
      endif !(N.LT.1)

      if( N .EQ. 1 )then
          WRITE(IPR,28) 
          RETURN 
      endif ! N.EQ.1

      HOLD=X(1) 
      DO  60 I=2,N 
          IF(X(I).NE.HOLD)GO TO 65 
   60 ENDDO 

      WRITE(IPR, 9)HOLD 

      IFLAG=1 
   65 HOLD=Y(1) 
      DO  70 I=2,N 
      IF(Y(I).NE.HOLD)GO TO 80 
   70 END DO 

      WRITE(IPR,19)HOLD 

      IFLAG=1 
   80 IF(IFLAG.EQ.1)RETURN 


   90 CONTINUE 

    9 FORMAT(/1x ,'>>>>> NON-FATAL DIAGNOSTIC--THE FIRST  INPUT ARGUMENT&
         & (A VECTOR) TO THE CORR   SUBROUTINE HAS ALL ELEMENTS = ', &
         E15.8,' <<<<<'/)


   19 FORMAT(/1x ,'>>>>> NON-FATAL DIAGNOSTIC--THE SECOND INPUT ARGUMENT&
             & (A VECTOR) TO THE CORR   SUBROUTINE HAS ALL ELEMENTS = ', &
             E15.8,' <<<<<'/)

   25 FORMAT(/1x, '>>>>> FATAL ERROR--THE THIRD  INPUT ARGUMENT TO THE  CORR&
              &   SUBROUTINE IS NON-POSITIVE <<<<<')

   28 FORMAT(/1x ,'>>>>> NON-FATAL DIAGNOSTIC--THE THIRD  INPUT ARGUMENT &
                  &TO THE CORR   SUBROUTINE HAS THE VALUE 1 <<<<<'/)

   47 FORMAT(1x, '>>>>> THE VALUE OF THE ARGUMENT IS ',I8   ,' <<<<<'/) 
!                                                                       
!-----START POINT-----------------------------------------------------  
!                                                                       
      XBAR=0.0d0
      YBAR=0.0d0 

      DO 100 I=1,N 
          XBAR=XBAR+X(I)
          YBAR=YBAR+Y(I)
  100 ENDDO 

      XBAR=XBAR/AN 
      YBAR=YBAR/AN 
!                                                                       
      SUM1=0.0d0 
      SUM2=0.0d0 
      SUM3=0.0d0 

      DO 200  I=1,N 

      SUM1=SUM1+(X(I)-XBAR)*(Y(I)-YBAR) !orig 
      SUM2=SUM2+(X(I)-XBAR)**2 !orig 
      SUM3=SUM3+(Y(I)-YBAR)**2 !orig 


  200 ENDDO 

      SUM2=SQRT(SUM2) 
      SUM3=SQRT(SUM3) 
      C   =SUM1/(SUM2*SUM3) 

      IF(IWRITE.NE.0)WRITE(IPR,205)N,C 

  205 FORMAT(1x, &
       'THE LINEAR        CORRELATION COEFFICIENT OF THE 2 SETS OF ', &
             I6,' OBSERVATIONS IS ',F14.5)                            

      RETURN 

      END                                           
