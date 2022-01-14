#include "Preprocessor"

submodule ( Gradient_Form ) Gradient_Kernel

  use Basics
  
  implicit none
  
contains


  module procedure Compute_CGS_Kernel

    integer ( KDI ) :: &
      iS, &
      iF, &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVP, iaVM, &
      lV, uV
    real ( KDR ) :: &
       fM,  fP, &  !-- f_Minus, f_Plus
      xAM, xAP
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption

    lV  =  1
    where ( shape ( XA )  >  1 )
      lV  =  oV
    end where
    
    uV  =  1
    where ( shape ( XA )  >  1 )
      uV  =  shape ( XA )  -  oV
    end where
    uV ( iD )  =  size ( XA, dim = iD )  -  oV  +  1 
      
    iaS  =  0
    iaS ( iD )  =  1
    
    if ( UseDevice ) then
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 4 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP private ( iF, iaVP, iaVM, fM, fP, xAM, xAP )
      do iS  =  1,  size ( iaSlctd )
        do kV  =  lV ( 3 ),  uV ( 3 ) 
          do jV  =  lV ( 2 ),  uV ( 2 )
            do iV  =  lV ( 1 ),  uV ( 1 )

              iF  =  iaSlctd ( iS )

              iaVM  =  [ iV, jV, kV ]  -  iaS
              iaVP  =  [ iV, jV, kV ]  +  iaS

              fM  =  F ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ), iF )
              fP  =  F ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF )

              xAM  =  XA ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
              xAP  =  XA ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

              dFdX ( iV, jV, kV, iS )  =  ( fP - fM ) / ( xAP - xAM )

            end do !-- iV
          end do !-- jV
        end do !-- kV
      end do !-- iS
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else !-- use host
              
      !$OMP parallel do collapse ( 4 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP private ( iF, iaVP, iaVM, fM, fP, xAM, xAP )
      do iS  =  1,  size ( iaSlctd )
        do kV  =  lV ( 3 ),  uV ( 3 ) 
          do jV  =  lV ( 2 ),  uV ( 2 )
            do iV  =  lV ( 1 ),  uV ( 1 )

              iF  =  iaSlctd ( iS )

              iaVM  =  [ iV, jV, kV ]  -  iaS
              iaVP  =  [ iV, jV, kV ]  +  iaS

              fM  =  F ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ), iF )
              fP  =  F ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF )

              xAM  =  XA ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
              xAP  =  XA ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

              dFdX ( iV, jV, kV, iS )  =  ( fP - fM ) / ( xAP - xAM )

            end do !-- iV
          end do !-- jV
        end do !-- kV
      end do !-- iS
      !$OMP end parallel do
      
    end if !-- UseDevice
        
  end procedure Compute_CGS_Kernel


end submodule Gradient_Kernel
