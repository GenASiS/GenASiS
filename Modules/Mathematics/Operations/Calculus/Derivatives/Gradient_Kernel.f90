#include "Preprocessor"

submodule ( Gradient_Form ) Gradient_Kernel

  use Basics 
  
  implicit none
  
contains 

  module procedure ComputeChart_SL_G_Kernel 

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &
      lV, uV
    real ( KDR ) :: &
      dX_L, dX_R, &
      dV_L, dV_R
    logical ( KDL ) :: &
      UseDevice
      
    real ( KDR ), pointer :: &
      Theta
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      UseLimiter
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    lV = 1
    where ( shape ( dV_I ) > 1 )
      lV = oV + 1
    end where
    lV ( iD ) = oV
    
    uV = 1
    where ( shape ( dV_I ) > 1 )
      uV = shape ( dV_I ) - oV
    end where
    uV ( iD ) = size ( dV_I, dim = iD ) - oV + 1 
      
    iaS = 0
    iaS ( iD ) = +1

    if ( present ( ThetaOption ) .and. present ( UseLimiterOption ) ) then
      !-- FIXME: XL refuses to compile association with optional variable
      !associate &
      !  ( Theta => ThetaOption, &
      !    UseLimiter => UseLimiterOption )
      
      Theta => ThetaOption
      UseLimiter => UseLimiterOption
      
      if ( UseDevice ) then

        !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
        !$OMP& private ( iV, jV, kV, iaVS, dV_L, dV_R, dX_L, dX_R )
        do kV = lV ( 3 ), uV ( 3 ) 
          do jV = lV ( 2 ), uV ( 2 )
            do iV = lV ( 1 ), uV ( 1 )

              iaVS = [ iV, jV, kV ] + iaS

              dV_L = dV_I ( iV, jV, kV )
              dV_R = dV_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
              dX_L = dX_I ( iV, jV, kV )
              dX_R = dX_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
                
              if ( UseLimiter ( iV, jV, kV ) > 0.0_KDR ) then
                dVdX ( iV, jV, kV ) &
                  = ( sign ( 0.5_KDR, dV_L ) + sign ( 0.5_KDR, dV_R ) ) &
                    * min ( abs ( Theta * dV_L / dX_L ), &
                            abs ( Theta * dV_R / dX_R ), &
                            abs ( ( dX_R ** 2  *  dV_L  +  dX_L ** 2  *  dV_R ) &
                                  / ( dX_L * dX_R * ( dX_L + dX_R ) ) ) )
              else
                dVdX ( iV, jV, kV )&
                  =  ( dX_R ** 2  *  dV_L  +  dX_L ** 2  *  dV_R ) &
                     / ( dX_L * dX_R * ( dX_L + dX_R ) )
              end if

            end do !-- iV
          end do !-- jV
        end do !-- kV
        !$OMP end OMP_TARGET_DIRECTIVE parallel do
      
      
      else

        !$OMP  parallel do collapse ( 3 ) &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
        !$OMP& private ( iV, jV, kV, iaVS, dV_L, dV_R, dX_L, dX_R )
        do kV = lV ( 3 ), uV ( 3 ) 
          do jV = lV ( 2 ), uV ( 2 )
            do iV = lV ( 1 ), uV ( 1 )

              iaVS = [ iV, jV, kV ] + iaS

              dV_L = dV_I ( iV, jV, kV )
              dV_R = dV_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
              dX_L = dX_I ( iV, jV, kV )
              dX_R = dX_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
                
              if ( UseLimiter ( iV, jV, kV ) > 0.0_KDR ) then
                dVdX ( iV, jV, kV ) &
                  = ( sign ( 0.5_KDR, dV_L ) + sign ( 0.5_KDR, dV_R ) ) &
                    * min ( abs ( Theta * dV_L / dX_L ), &
                            abs ( Theta * dV_R / dX_R ), &
                            abs ( ( dX_R ** 2  *  dV_L  +  dX_L ** 2  *  dV_R ) &
                                  / ( dX_L * dX_R * ( dX_L + dX_R ) ) ) )
              else
                dVdX ( iV, jV, kV )&
                  =  ( dX_R ** 2  *  dV_L  +  dX_L ** 2  *  dV_R ) &
                     / ( dX_L * dX_R * ( dX_L + dX_R ) )
              end if

            end do !-- iV
          end do !-- jV
        end do !-- kV
        !$OMP end parallel do

      end if
      
      !-- FIXME: see above
      !end associate !-- Theta
    else if ( present ( ThetaOption ) ) then
      !-- FIXME: see above
      !associate ( Theta => ThetaOption )
      Theta => ThetaOption
      
      if ( UseDevice ) then
      
        !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
        !$OMP& private ( iV, jV, kV, iaVS, dV_L, dV_R, dX_L, dX_R )
        do kV = lV ( 3 ), uV ( 3 ) 
          do jV = lV ( 2 ), uV ( 2 )
            do iV = lV ( 1 ), uV ( 1 )

              iaVS = [ iV, jV, kV ] + iaS

              dV_L = dV_I ( iV, jV, kV )
              dV_R = dV_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
              dX_L = dX_I ( iV, jV, kV )
              dX_R = dX_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
                
              dVdX ( iV, jV, kV ) &
                = ( sign ( 0.5_KDR, dV_L ) + sign ( 0.5_KDR, dV_R ) ) &
                  * min ( abs ( Theta * dV_L / dX_L ), &
                          abs ( Theta * dV_R / dX_R ), &
                          abs ( ( dX_R ** 2  *  dV_L  +  dX_L ** 2  *  dV_R ) &
                                / ( dX_L * dX_R * ( dX_L + dX_R ) ) ) )

            end do !-- iV
          end do !-- jV
        end do !-- kV
        !$OMP end OMP_TARGET_DIRECTIVE parallel do
      
      else

        !$OMP parallel do collapse ( 3 ) &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
        !$OMP& private ( iV, jV, kV, iaVS, dV_L, dV_R, dX_L, dX_R )
        do kV = lV ( 3 ), uV ( 3 ) 
          do jV = lV ( 2 ), uV ( 2 )
            do iV = lV ( 1 ), uV ( 1 )

              iaVS = [ iV, jV, kV ] + iaS

              dV_L = dV_I ( iV, jV, kV )
              dV_R = dV_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
              dX_L = dX_I ( iV, jV, kV )
              dX_R = dX_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
                
              dVdX ( iV, jV, kV ) &
                = ( sign ( 0.5_KDR, dV_L ) + sign ( 0.5_KDR, dV_R ) ) &
                  * min ( abs ( Theta * dV_L / dX_L ), &
                          abs ( Theta * dV_R / dX_R ), &
                          abs ( ( dX_R ** 2  *  dV_L  +  dX_L ** 2  *  dV_R ) &
                                / ( dX_L * dX_R * ( dX_L + dX_R ) ) ) )

            end do !-- iV
          end do !-- jV
        end do !-- kV
        !$OMP end parallel do

      end if

      !-- FIXME: see above
      !end associate !-- Theta
    else
    
      if ( UseDevice ) then

        !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
        !$OMP& private ( iV, jV, kV, iaVS, dV_L, dV_R, dX_L, dX_R )
        do kV = lV ( 3 ), uV ( 3 ) 
          do jV = lV ( 2 ), uV ( 2 )
            do iV = lV ( 1 ), uV ( 1 )

              iaVS = [ iV, jV, kV ] + iaS

              dV_L = dV_I ( iV, jV, kV )
              dV_R = dV_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
              dX_L = dX_I ( iV, jV, kV )
              dX_R = dX_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
                
              dVdX ( iV, jV, kV )&
                =  ( dX_R ** 2  *  dV_L  +  dX_L ** 2  *  dV_R ) &
                   / ( dX_L * dX_R * ( dX_L + dX_R ) )
              
            end do !-- iV
          end do !-- jV
        end do !-- kV
        !$OMP end OMP_TARGET_DIRECTIVE parallel do

      else

        !$OMP parallel do collapse ( 3 ) &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
        !$OMP& private ( iV, jV, kV, iaVS, dV_L, dV_R, dX_L, dX_R )
        do kV = lV ( 3 ), uV ( 3 ) 
          do jV = lV ( 2 ), uV ( 2 )
            do iV = lV ( 1 ), uV ( 1 )

              iaVS = [ iV, jV, kV ] + iaS

              dV_L = dV_I ( iV, jV, kV )
              dV_R = dV_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
              dX_L = dX_I ( iV, jV, kV )
              dX_R = dX_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
                
              dVdX ( iV, jV, kV )&
                =  ( dX_R ** 2  *  dV_L  +  dX_L ** 2  *  dV_R ) &
                   / ( dX_L * dX_R * ( dX_L + dX_R ) )
              
            end do !-- iV
          end do !-- jV
        end do !-- kV
        !$OMP end parallel do

      end if

    end if
    
  end procedure ComputeChart_SL_G_Kernel


end submodule Gradient_Kernel
