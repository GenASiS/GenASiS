module Tally_F_P__Form

  !-- Tally_Fluid_Perfect__Form

  use Basics
  use Mathematics
  use Fluid_P__Template
  use Tally_F_D__Form
  
  implicit none
  private
  
    integer ( KDI ), private, parameter :: &
      N_INTEGRALS_PERFECT = 2
  
  type, public, extends ( Tally_F_D_Form ) :: Tally_F_P_Form
    integer ( KDI ) :: &
      N_INTEGRALS_PERFECT = N_INTEGRALS_PERFECT, &
      FLUID_ENERGY        = 0, &
      INTERNAL_ENERGY     = 0
  contains
    procedure, public, pass :: &
      InitializeFluid
    procedure, public, pass :: &
      SelectVariables
    final :: &
      Finalize
    procedure, public, pass :: &
      ComputeInteriorIntegrand
    procedure, public, pass :: &
      ComputeBoundaryIntegrand_CSL
  end type Tally_F_P_Form


contains


  subroutine InitializeFluid &
               ( T, A, MassUnitOption, EnergyUnitOption, MomentumUnitOption, &
                 AngularMomentumUnitOption )
    
    class ( Tally_F_P_Form ), intent ( inout ) :: &
      T
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      MassUnitOption, &
      EnergyUnitOption, &
      MomentumUnitOption, &
      AngularMomentumUnitOption

    integer ( KDI ) :: &
      oI     !-- oIntegral

    oI = T % N_INTEGRALS_DUST
    if ( T % N_INTEGRALS == 0 ) &
      T % N_INTEGRALS = oI + T % N_INTEGRALS_PERFECT
    
    call T % Tally_F_D_Form % Initialize &
           ( A, MassUnitOption, EnergyUnitOption, MomentumUnitOption, &
             AngularMomentumUnitOption )
    
    T % FLUID_ENERGY     = oI + 1
    T % INTERNAL_ENERGY  = oI + 2
    
    T % Variable ( oI + 1 : oI + T % N_INTEGRALS_PERFECT ) &
      = [ 'FluidEnergy                    ', &
          'InternalEnergy                 ' ]

    if ( present ( EnergyUnitOption ) ) then
      T % Unit ( oI + 1 : oI + T % N_INTEGRALS_PERFECT ) &
        = [ EnergyUnitOption, EnergyUnitOption ]
    else
      T % Unit ( oI + 1 : oI + T % N_INTEGRALS_PERFECT ) &
        = [ UNIT % IDENTITY, UNIT % IDENTITY ]
    end if

    call T % SelectVariables ( A )

    T % Atlas => A

  end subroutine InitializeFluid

  
  subroutine SelectVariables ( T, A ) 
    
    class ( Tally_F_P_Form ), intent ( inout ) :: &
      T
    class ( AtlasHeaderForm ), intent ( in ) :: &
      A

    class ( GeometryFlatForm ), pointer :: &
      G

    if ( allocated ( T % iaSelected ) ) &
      deallocate ( T % iaSelected )

    select type ( A )
    class is ( Atlas_SC_Form )

    G => A % Geometry ( )
    
    select type ( G )
    type is ( GeometryFlatForm )
      T % nSelected = 10
      allocate ( T % iaSelected ( T % nSelected ) )
      T % iaSelected &
        = [ T % BARYON_NUMBER, T % MOMENTUM, T % FLUID_ENERGY, &
            T % INTERNAL_ENERGY, T % KINETIC_ENERGY, T % ANGULAR_MOMENTUM ]
    
    ! type is ( GeometryNewtonianForm )
    !   allocate ( T % iaSelected ( 8 ) )
    !   T % iaSelected &
    !     = [ T % BARYON_NUMBER, T % MOMENTUM, T % KINETIC_ENERGY, &
    !         T % FLUID_ENERGY, T % GRAVITATIONAL_ENERGY, &
    !         T % TOTAL_ENERGY ]
    
    class default 
      call Show ( 'This type is not implemented yet', CONSOLE % WARNING )
      call Show ( 'Tally_F_P_Form % SelectVariables', &
                  'Subroutine', CONSOLE % WARNING )
            
    end select !-- G

    class default 
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Tally_F_P__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SelectVariables', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    nullify ( G )

  end subroutine SelectVariables
  

  impure elemental subroutine Finalize ( T )
  
    type ( Tally_F_P_Form ), intent ( inout ) :: &
      T

    !-- Trigger finalization of parent

  end subroutine Finalize
  

  subroutine ComputeInteriorIntegrand &
               ( T, Integrand, C, G, nDimensions )

    class ( Tally_F_P_Form ), intent ( inout ) :: &
      T
    type ( Real_1D_Form ), dimension ( : ), intent ( inout ) :: &
      Integrand
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    type ( GeometryFlatForm ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      nDimensions

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iI     !-- iIntegral

    select type ( C )
    class is ( Fluid_P_Template )

    call T % Tally_F_D_Form &
           % ComputeInteriorIntegrand ( Integrand, C, G, nDimensions )

    associate &
      ( CE => C % Value ( :, C % CONSERVED_ENERGY ), &
        IE => C % Value ( :, C % INTERNAL_ENERGY ) )

    do iS = 1, T % nSelected
      iI = T % iaSelected ( iS )
      if ( iI == T % FLUID_ENERGY ) then
        call Copy ( CE, Integrand ( iS ) % Value )
      else if ( iI == T % INTERNAL_ENERGY ) then
        call Copy ( IE, Integrand ( iS ) % Value )
      end if !-- iI
    end do !-- iS

    end associate !-- CE, etc.
    end select !-- F

  end subroutine ComputeInteriorIntegrand


  subroutine ComputeBoundaryIntegrand_CSL &
               ( T, Integrand, C, CSL, G, BoundaryFluence )

    class ( Tally_F_P_Form ), intent ( inout ) :: &
      T
    type ( Real_3D_Form ), dimension ( :, : ), intent ( inout ) :: &
      Integrand
    class ( CurrentTemplate ), intent ( in ) :: & 
      C
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    type ( GeometryFlatForm ), intent ( in ) :: &
      G
    type ( Real_3D_Form ), dimension ( :, : ), intent ( in ) :: &
      BoundaryFluence

    integer ( KDI ) :: &
      iF, &   !-- iFace
      iS, &   !-- iSelected
      iI, &   !-- iIntegral
      iFluence, &
      iEnergy

    select type ( C )
    class is ( Fluid_P_Template )

    call T % Tally_F_D_Form % ComputeBoundaryIntegrand_CSL &
           ( Integrand, C, CSL, G, BoundaryFluence )

    do iFluence = 1, C % N_CONSERVED
      if ( C % iaConserved ( iFluence ) == C % CONSERVED_ENERGY ) &
        iEnergy = iFluence
    end do !-- iFluence

    associate ( Cnnct => CSL % Atlas % Connectivity )
    do iF = 1, Cnnct % nFaces
      associate ( CE => BoundaryFluence ( iEnergy, iF ) % Value )
      do iS = 1, T % nSelected
        iI = T % iaSelected ( iS )
        if ( iI == T % FLUID_ENERGY ) then
          call CopyCollapse ( CE, Integrand ( iS, iF ) % Value )
        end if !-- iI
      end do !-- iS
      end associate !-- CE
    end do !-- iF
    end associate !-- Cnnct

    end select !-- C

  end subroutine ComputeBoundaryIntegrand_CSL


end module Tally_F_P__Form
