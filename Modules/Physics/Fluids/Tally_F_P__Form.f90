module Tally_F_P__Form

  !-- Tally_Fluid_Perfect__Form

  use Basics
  use Mathematics
  use Gravitations
  use Units_F__Form
  use Fluid_P__Form
  use Tally_F_D__Form

  implicit none
  private
  
    integer ( KDI ), private, parameter :: &
      N_INTEGRALS_P = 3
  
  type, public, extends ( Tally_F_D_Form ) :: Tally_F_P_Form
    integer ( KDI ) :: &
      N_INTEGRALS_P   = N_INTEGRALS_P, &
      FLUID_ENERGY    = 0, &
      INTERNAL_ENERGY = 0, &
      ENTROPY         = 0
  contains
    procedure, public, pass :: &
      InitializeFluid
    procedure, public, pass :: &
      SelectVariables
    final :: &
      Finalize
    procedure, public, pass :: &
      ComputeInteriorIntegrand_G
    procedure, public, pass :: &
      ComputeBoundaryIntegrand_G
  end type Tally_F_P_Form

    private :: &
      ComputeDensity_S


contains


  subroutine InitializeFluid ( T, G, Units )
    
    class ( Tally_F_P_Form ), intent ( inout ) :: &
      T
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G
    class ( Units_F_Form ), intent ( in ) :: &
      Units

    integer ( KDI ) :: &
      oI     !-- oIntegral

    oI  =  T % N_INTEGRALS_D
    if ( T % nIntegrals  ==  0 ) &
      T % nIntegrals  =  oI  +  T % N_INTEGRALS_P
    
    call T % Tally_F_D_Form % Initialize ( G, Units )
    
    T % FLUID_ENERGY     =  oI + 1
    T % INTERNAL_ENERGY  =  oI + 2
    T % ENTROPY          =  oI + 3
    
    T % Variable ( oI + 1 : oI + T % N_INTEGRALS_P ) &
      = [ 'FluidEnergy   ', &
          'InternalEnergy', &
          'Entropy       ' ]

    T % Unit ( oI + 1 : oI + T % N_INTEGRALS_P ) &
      = [ Units % Energy, &
          Units % Energy, &
          Units % Energy  /  Units % Temperature ]

    call T % SelectVariables ( )

  end subroutine InitializeFluid

  
  subroutine SelectVariables ( T ) 
    
    class ( Tally_F_P_Form ), intent ( inout ) :: &
      T

    if ( allocated ( T % iaSelected ) ) &
      deallocate ( T % iaSelected )

    select type ( G  =>  T % Geometry )
    type is ( Gravitation_G_Form )
      T % nSelected  =  11
      allocate ( T % iaSelected ( T % nSelected ) )
      T % iaSelected &
        = [ T % BARYON_NUMBER, &
            T % MOMENTUM, &
            T % FLUID_ENERGY, &
            T % INTERNAL_ENERGY, &
            T % KINETIC_ENERGY, &
            T % ANGULAR_MOMENTUM, &
            T % ENTROPY ]
    class is ( Gravitation_N_H_Form )
      T % nSelected  =  13
      allocate ( T % iaSelected ( T % nSelected ) )
      T % iaSelected &
        = [ T % BARYON_NUMBER, &
            T % MOMENTUM, &
            T % FLUID_ENERGY, &
            T % INTERNAL_ENERGY, &
            T % KINETIC_ENERGY, &
            T % ANGULAR_MOMENTUM, &
            T % ENTROPY, &
            T % GRAVITATIONAL_ENERGY, &
            T % TOTAL_ENERGY ]
    class default 
      call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
      call Show ( 'Tally_F_P_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SelectVariables', 'subroutine', CONSOLE % ERROR )            
      call PROGRAM_HEADER % Abort ( )
    end select !-- G

  end subroutine SelectVariables
  

  impure elemental subroutine Finalize ( T )
  
    type ( Tally_F_P_Form ), intent ( inout ) :: &
      T

  end subroutine Finalize
  

  subroutine ComputeInteriorIntegrand_G ( T, CS )

    class ( Tally_F_P_Form ), intent ( inout ) :: &
      T
    class ( FieldSet_BM_Form ), intent ( in ) :: &
      CS

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iI     !-- iIntegral

    call T % Tally_F_D_Form % ComputeInteriorIntegrand_G ( CS )

    select type ( CS )
      class is ( Fluid_P_Form )
    associate &
      ( CSV  =>  CS % Storage_GS % Value, &
        IV  =>  T % InteriorIntegral % Integrand % Storage_GS % Value )
    associate &
      ( BE  =>  CSV ( :, CS % ENERGY_DENSITY_B ), &
        IE  =>  CSV ( :, CS % ENERGY_DENSITY_C ), &
        SB  =>  CSV ( :, CS % ENTROPY_PER_BARYON ), &
        N   =>  CSV ( :, CS % BARYON_DENSITY_C ) )

    do iS  =  1, T % nSelected
      iI  =  T % iaSelected ( iS )
      if ( iI  ==  T % FLUID_ENERGY ) then
        call Copy ( BE, IV ( :, iS ) )
      else if ( iI  ==  T % TOTAL_ENERGY ) then
        call Copy ( BE, IV ( :, iS ) )
      else if ( iI  ==  T % INTERNAL_ENERGY ) then
        call Copy ( IE, IV ( :, iS ) )
      else if ( iI  ==  T % ENTROPY ) then
        call ComputeDensity_S ( SB, N, IV ( :, iS ) )
      end if !-- iI
    end do !-- iS

    end associate !-- CE, etc.
    end associate !-- CSV, etc.
    end select !-- CS

  end subroutine ComputeInteriorIntegrand_G


  subroutine ComputeBoundaryIntegrand_G ( T, CS, C, BF )

    class ( Tally_F_P_Form ), intent ( inout ) :: &
      T
    class ( FieldSet_BM_Form ), intent ( in ) :: &
      CS
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    type ( Real_3D_Form ), dimension ( :, : ), intent ( in ) :: &
      BF

    integer ( KDI ) :: &
      iF, &   !-- iFace
      iS, &   !-- iSelected
      iI, &   !-- iIntegral
      iFluence, &
      iEnergy

    call T % Tally_F_D_Form % ComputeBoundaryIntegrand_G ( CS, C, BF )

    select type ( CS )
      class is ( Fluid_P_Form )
    associate &
      ( I   =>  T % BoundaryIntegral % Integrand, &
        Cy  =>  C % Connectivity )
       
    do iFluence  =  1,  CS % nBalanced
      if ( CS % iaBalanced ( iFluence )  ==  CS % ENERGY_DENSITY_B ) &
        iEnergy = iFluence
    end do !-- iFluence

    do iF  =  1, Cy % nFaces
      associate  ( CE  =>  BF ( iEnergy, iF ) % Value )
      do iS  =  1, T % nSelected
        iI  =  T % iaSelected ( iS )
        if ( iI  ==  T % FLUID_ENERGY ) then
          call Copy ( CE, I ( iS, iF ) % Value )
        else if ( iI  ==  T % TOTAL_ENERGY ) then
          call Copy ( CE, I ( iS, iF ) % Value )
        end if !-- iI
      end do !-- iS
      end associate !-- CE
    end do !-- iF

    end associate !-- I, etc.
    end select !-- CS

  end subroutine ComputeBoundaryIntegrand_G


  subroutine ComputeDensity_S ( SB, N, I )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      SB, &
      N
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      I

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( I )

    !$OMP parallel do
    do iV  =  1,  nV
      I ( iV )  =  SB ( iV )  *  N ( iV )
    end do
    !$OMP end parallel do

  end subroutine ComputeDensity_S


end module Tally_F_P__Form
