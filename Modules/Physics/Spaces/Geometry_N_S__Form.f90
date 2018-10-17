!-- Geometry_N_S represents Newtonian geometry, plus an overloading of
!   ComputeReconstruction to add reconstruction of gradient of Phi for use
!   in the divergence of a gravitational stress as an alternative to a 
!   gravitational force source term.

module Geometry_N_S__Form

  !-- Geometry_Newtonian_Stress__Form

  use Basics
  use Mathematics
  use Geometry_N__Form

  implicit none
  private

  type, public, extends ( Geometry_N_Form ) :: Geometry_N_S_Form
  contains
    procedure, public, pass :: &
      InitializeAllocate_G
    procedure, private, pass ( G ) :: &
      ComputeReconstruction_CSL
    final :: &
      Finalize
  end type Geometry_N_S_Form

    private :: &
      ComputeReconstruction_CSL_Kernel

contains


  subroutine InitializeAllocate_G &
               ( G, CoordinateSystem, CoordinateUnit, nValues, &
                 VariableOption, VectorOption, NameOption, ClearOption, &
                 UnitOption, VectorIndicesOption )

    class ( Geometry_N_S_Form ), intent ( inout ) :: &
      G
    character ( * ), intent ( in ) :: &
      CoordinateSystem
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      CoordinateUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption

    if ( G % Type == '' ) &
      G % Type = 'a Geometry_N_S'

    call G % Geometry_N_Form % Initialize &
           ( CoordinateSystem, CoordinateUnit, nValues, &
             VariableOption, VectorOption, NameOption, ClearOption, &
             UnitOption, VectorIndicesOption )

  end subroutine InitializeAllocate_G


  subroutine ComputeReconstruction_CSL ( G_I, CSL, G, nDimensions, iDimension )

    type ( StorageForm ), intent ( inout ) :: &
      G_I
    class ( ChartHeader_SL_Form ), intent ( in ) :: &
      CSL
    class ( Geometry_N_S_Form ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      iDimension

    integer ( KDI ) :: &
      iD
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      dX_L, dX_R, &
      GradPhi, &
      GradPhi_I

    call G % GeometryFlatForm % ComputeReconstruction &
           ( G_I, CSL, nDimensions, iDimension )

    call CSL % SetVariablePointer &
           ( G % Value ( :, G % WIDTH_LEFT_U ( iDimension ) ), dX_L )
    call CSL % SetVariablePointer &
           ( G % Value ( :, G % WIDTH_RIGHT_U ( iDimension ) ), dX_R )

    do iD = 1, nDimensions
      call CSL % SetVariablePointer &
             ( G % Value ( :, G % POTENTIAL_GRADIENT_D ( iD ) ), GradPhi )
      call CSL % SetVariablePointer &
             ( G_I % Value ( :, G % POTENTIAL_GRADIENT_D ( iD ) ), GradPhi_I )
      call ComputeReconstruction_CSL_Kernel &
             ( GradPhi, dX_L, dX_R, iDimension, &
               CSL % nGhostLayers ( iDimension ), GradPhi_I )
    end do

    nullify ( dX_L, dX_R, GradPhi, GradPhi_I )

  end subroutine ComputeReconstruction_CSL


  impure elemental subroutine Finalize ( G )

    type ( Geometry_N_S_Form ), intent ( inout ) :: &
      G

    !-- Empty routine to trigger finalization of parent type

  end subroutine Finalize


  subroutine ComputeReconstruction_CSL_Kernel &
               ( GradPhi, dX_L, dX_R, iD, oV, GradPhi_I )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      GradPhi, &
      dX_L, dX_R
    integer ( KDI ), intent ( in ) :: &
      iD, &
      oV   
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      GradPhi_I

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &
      lV, uV

    lV = 1
    where ( shape ( GradPhi ) > 1 )
      lV = oV + 1
    end where
    
    uV = 1
    where ( shape ( GradPhi ) > 1 )
      uV = shape ( GradPhi ) - oV
    end where
    uV ( iD ) = size ( GradPhi, dim = iD ) - oV + 1 
      
    iaS = 0
    iaS ( iD ) = -1
      
    !$OMP parallel do private ( iV, jV, kV, iaVS )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )

          iaVS = [ iV, jV, kV ] + iaS

          GradPhi_I ( iV, jV, kV )  &
            =  (    GradPhi ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )  &
                      *  dX_L ( iV, jV, kV )  &
                 +  GradPhi ( iV, jV, kV )  &
                      *  dX_R ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) )  &
               /  ( dX_R ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )  &
                    +  dX_L ( iV, jV, kV ) )

        end do !-- iV
      end do !-- jV
    end do !-- kV
    !$OMP end parallel do
      
  end subroutine ComputeReconstruction_CSL_Kernel


end module Geometry_N_S__Form
