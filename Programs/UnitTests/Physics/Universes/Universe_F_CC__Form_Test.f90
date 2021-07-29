program Universe_F_CC__Form_Test

  !-- Universe_Fluid_Box__Form_Test

  use Basics
  use Mathematics
  use Fluids
  use Universe_F_CC__Form

  implicit none

  type ( Universe_F_CC_Form ), allocatable :: &
    U

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Universe_F_CC__Form_Test', DimensionalityOption = '2D' )

  allocate ( U )
  call U % Initialize &
         ( FluidType = 'DUST', &
           GravitationType = 'NEWTON_SG' )
  call U % Show ( )
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

    end associate !-- FV, etc.
    end select !-- C
    end select !-- A
    end select !-- F
    end select !-- I

  end subroutine SetFluid


end program Universe_F_CC__Form_Test
