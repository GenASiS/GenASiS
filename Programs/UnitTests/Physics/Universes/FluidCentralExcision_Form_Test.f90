program FluidCentralExcision_Form_Test

  use Basics
  use FluidCentralExcision_Form

  implicit none

  type ( FluidCentralExcisionForm ), allocatable :: &
    FCE

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FluidCentralExcision_Form_Test', DimensionalityOption = '2D' )

  allocate ( FCE )
  call FCE % Initialize &
         ( PROGRAM_HEADER % Name, &
           FluidType = 'IDEAL', &
           GeometryType = 'NEWTONIAN', &
           DimensionlessOption = .true. )
  call FCE % OpenManifoldStreams ( VerboseStreamOption = .true. )

  call SetFluid ( FCE )
  call FCE % Write ( )

  call CoarsenFluid ( FCE )
  call FCE % Write ( )

  deallocate ( FCE )

  deallocate ( PROGRAM_HEADER )

end program FluidCentralExcision_Form_Test


subroutine SetFluid ( FCE )

  use Basics
  use Mathematics
  use StressEnergies
  use FluidCentralExcision_Form

  type ( FluidCentralExcisionForm ), intent ( inout ) :: &
    FCE

  integer ( KDI ) :: &
    iC, jC, kC  !-- iCell, etc.
  real ( KDR ) :: &
    dS
  real ( KDR ), dimension ( :, :, : ), pointer :: &
    S_2, S_3
  class ( Fluid_D_Form ), pointer :: &
    F

  select type ( FA => FCE % Current_ASC )
  class is ( Fluid_ASC_Form )
  F => FA % Fluid_D ( )

  select type ( PS => FCE % PositionSpace )
  class is ( Atlas_SC_Form )

  select type ( C => PS % Chart )
  class is ( Chart_SLD_C_Template )

  associate &
    ( iaB => C % iaBrick, &
      nCB => C % nCellsBrick )

  call C % SetVariablePointer &
         ( F % Value ( :, F % MOMENTUM_DENSITY_D ( 2 ) ), S_2 )
  call C % SetVariablePointer &
         ( F % Value ( :, F % MOMENTUM_DENSITY_D ( 3 ) ), S_3 )

  if ( C % nDimensions < 2 ) &
    return

  dS  =  CONSTANT % PI  /  C % nCells ( 2 )

  do kC = 1, nCB ( 3 )
    do jC = 1, nCB ( 2 )
      do iC = 1, nCB ( 1 )
        S_2 ( iC, jC, kC )  &
          =  sin ( ( ( iaB ( 2 ) - 1 ) * nCB ( 2 )  +  jC - 0.5_KDR )  *  dS )
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
          =  sin ( ( ( iaB ( 3 ) - 1 ) * nCB ( 3 )  +  kC - 0.5_KDR )  *  dS )
      end do
    end do
  end do

  end associate !-- iaB, etc.
  end select !-- C
  end select !-- PS
  end select !-- FA
  nullify ( F, S_2, S_3 )

end subroutine SetFluid


subroutine CoarsenFluid ( FCE )

  use Basics
  use StressEnergies
  use FluidCentralExcision_Form

  type ( FluidCentralExcisionForm ), intent ( inout ) :: &
    FCE

  type ( StorageForm ) :: &
    S
  class ( Fluid_D_Form ), pointer :: &
    F

  select type ( FA => FCE % Current_ASC )
  class is ( Fluid_ASC_Form )
  F => FA % Fluid_D ( )

  call S % Initialize ( F, iaSelectedOption = F % iaConserved )
  call FCE % CoarsenSingularities ( S )

  end select !-- FA
  nullify ( F )

end subroutine CoarsenFluid
