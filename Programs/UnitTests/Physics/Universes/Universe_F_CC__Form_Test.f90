program Universe_F_CC__Form_Test

  !-- Universe_Fluid_Box__Form_Test

  use Basics
  use Mathematics
  use Fluids
  use Universe_F_CC__Form

  implicit none

  type ( TimerForm ), pointer :: &
    T_A, &
    T_W
  type ( Universe_F_CC_Form ), allocatable, target :: &
    U

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Universe_F_CC__Form_Test', DimensionalityOption = '2D' )

  allocate ( U )
  call U % Initialize &
         ( FluidType = 'DUST', &
           GravitationType = 'NEWTON_SG' )
  call U % Show ( )

  associate ( I  =>  U % Integrator )

  I % System  =>  U

  call SetFluid ( )

  T_A  =>  I % Timer_A ( )
  call T_A % Start ( )
  call I % Analyze ( T_A )
  call T_A % Stop ( )

  T_W  =>  I % Timer_W ( )
  call T_W % Start ( )
  call I % Write ( T_W )
  call T_W % Stop ( )

  end associate !-- I

  deallocate ( U )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine SetFluid ( )

    integer ( KDI ) :: &
      iC, jC, kC  !-- iCell, etc.
    real ( KDR ) :: &
      dS
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      S_1, S_2, S_3

    call Show ( 'Setting Fluid' )

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_D_Form )
    select type ( A  =>  I % X )
      class is ( Atlas_SCG_Form )
    select type ( C  =>  A % Chart_GS )
      class is ( Chart_GS_C_Form )
    associate &
      (  FV  =>  F % Storage_GS % Value, &
        nCB  =>  C % nCellsBrick, &
        iaB  =>  C % iaBrick )

    call C % SetFieldPointer ( FV ( :, F % MOMENTUM_DENSITY_D_1 ), S_1 )
    call C % SetFieldPointer ( FV ( :, F % MOMENTUM_DENSITY_D_2 ), S_2 )
    call C % SetFieldPointer ( FV ( :, F % MOMENTUM_DENSITY_D_3 ), S_3 )

    do kC = 1, nCB ( 3 )
      do jC = 1, nCB ( 2 )
        do iC = 1, nCB ( 1 )
          S_1 ( iC, jC, kC )  =  1.0_KDR
          S_2 ( iC, jC, kC )  =  2.0_KDR
          S_3 ( iC, jC, kC )  =  3.0_KDR
        end do
      end do
    end do

    if ( C % nDimensions < 2 ) &
      return

    dS  =  CONSTANT % PI  /  C % nCells ( 2 )

    do kC = 1, nCB ( 3 )
      do jC = 1, nCB ( 2 )
        do iC = 1, nCB ( 1 )
          S_2 ( iC, jC, kC )  &
            =  S_2 ( iC, jC, kC )  &
               +  sin ( ( ( iaB ( 2 ) - 1 ) * nCB ( 2 )  +  jC - 0.5_KDR )  &
                        *  dS )
        end do
      end do
    end do

    if ( C % nDimensions < 3 ) &
      return

    dS  =  2.0_KDR * CONSTANT % PI  /  C % nCells ( 3 )

    do kC = 1, nCB ( 3 )
      do jC = 1, nCB ( 2 )
        do iC = 1, nCB ( 1 )
          S_3 ( iC, jC, kC )  &
            =  S_3 ( iC, jC, kC )  &
               +  sin ( ( ( iaB ( 3 ) - 1 ) * nCB ( 3 )  +  kC - 0.5_KDR )  &
                        *  dS )
        end do
      end do
    end do

    end associate !-- FV, etc.
    end select !-- C
    end select !-- A
    end select !-- F
    end select !-- I

  end subroutine SetFluid


end program Universe_F_CC__Form_Test
