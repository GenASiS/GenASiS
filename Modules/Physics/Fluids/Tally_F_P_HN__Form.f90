module Tally_F_P_HN__Form

  !-- Tally_Fluid_Perfect_HeavyNucleus_Form

  use Basics
  use Mathematics
  use Gravitations
  use Units_F__Form
  use Fluid_P_HN__Form
  use Tally_F_P__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_INTEGRALS_HN = 1
  
  type, public, extends ( Tally_F_P_Form ) :: Tally_F_P_HN_Form
    integer ( KDI ) :: &
      N_INTEGRALS_HN = N_INTEGRALS_HN, &
      ELECTRON_NUMBER = 0
  contains
    procedure, private, pass :: &
      InitializeFluid
    procedure, public, pass :: &
      SelectVariables
    final :: &
      Finalize
    procedure, public, pass :: &
      ComputeInteriorIntegrand_G
    procedure, public, pass :: &
      ComputeBoundaryIntegrand_G
  end type Tally_F_P_HN_Form


contains


  subroutine InitializeFluid ( T, G, Units )
    
    class ( Tally_F_P_HN_Form ), intent ( inout ) :: &
      T
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G
    class ( Units_F_Form ), intent ( in ) :: &
      Units

    integer ( KDI ) :: &
      oI     !-- oIntegral

    oI  =  T % N_INTEGRALS_D  +  T % N_INTEGRALS_P
    if ( T % nIntegrals  ==  0 ) &
      T % nIntegrals  =  oI  +  T % N_INTEGRALS_HN
    
    call T % Tally_F_P_Form % Initialize ( G, Units )
    
    T % ELECTRON_NUMBER  =  oI + 1
    
    T % Variable ( oI + 1 : oI + T % N_INTEGRALS_HN ) &
      = [ 'ElectronNumber' ]

    T % Unit ( oI + 1 : oI + T % N_INTEGRALS_HN ) &
      = [ Units % Number ]

    call T % SelectVariables ( )

  end subroutine InitializeFluid

  
  subroutine SelectVariables ( T ) 
    
    class ( Tally_F_P_HN_Form ), intent ( inout ) :: &
      T

    if ( allocated ( T % iaSelected ) ) &
      deallocate ( T % iaSelected )

    select type ( G  =>  T % Geometry )
    type is ( Gravitation_G_Form )
      T % nSelected = 12
      allocate ( T % iaSelected ( T % nSelected ) )
      T % iaSelected &
        = [ T % BARYON_NUMBER, &
            T % ELECTRON_NUMBER, &
            T % MOMENTUM, &
            T % FLUID_ENERGY, &
            T % INTERNAL_ENERGY, &
            T % KINETIC_ENERGY, &
            T % ANGULAR_MOMENTUM, &
            T % ENTROPY ]    
    class is ( Gravitation_N_H_Form )
      T % nSelected = 14
      allocate ( T % iaSelected ( T % nSelected ) )
      T % iaSelected &
        = [ T % BARYON_NUMBER, &
            T % ELECTRON_NUMBER, &
            T % MOMENTUM, &
            T % FLUID_ENERGY, &
            T % INTERNAL_ENERGY, &
            T % KINETIC_ENERGY, &
            T % ANGULAR_MOMENTUM, &
            T % ENTROPY, &
            T % GRAVITATIONAL_ENERGY, &
            T % TOTAL_ENERGY ]    
    class default 
      call Show ( 'This type is not implemented yet', CONSOLE % WARNING )
      call Show ( 'Tally_F_P_HN__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SelectVariables', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- G

  end subroutine SelectVariables
  

  impure elemental subroutine Finalize ( T )
  
    type ( Tally_F_P_HN_Form ), intent ( inout ) :: &
      T

  end subroutine Finalize
  

  subroutine ComputeInteriorIntegrand_G ( T, CS )

    class ( Tally_F_P_HN_Form ), intent ( inout ) :: &
      T
    class ( FieldSet_BM_Form ), intent ( in ) :: &
      CS

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iI     !-- iIntegral

    call T % Tally_F_P_Form % ComputeInteriorIntegrand_G ( CS )

    select type ( CS )
      class is ( Fluid_P_HN_Form )
    associate &
      ( CSV  =>  CS % Storage_GS % Value, &
        IV  =>  T % InteriorIntegral % Integrand % Storage_GS % Value )
    associate &
      ( DE  =>  CSV ( :, CS % ELECTRON_DENSITY_B ) )

    do iS  =  1,  T % nSelected
      iI  =  T % iaSelected ( iS )
      if ( iI  ==  T % ELECTRON_NUMBER ) then
        call Copy ( DE, IV ( :, iS ) )
      end if !-- iI
    end do !-- iS

    end associate !-- DE
    end associate !-- CSV, etc.
    end select !-- CS

  end subroutine ComputeInteriorIntegrand_G


  subroutine ComputeBoundaryIntegrand_G ( T, CS, C, BF )

    class ( Tally_F_P_HN_Form ), intent ( inout ) :: &
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
      iElectron

    call T % Tally_F_P_Form % ComputeBoundaryIntegrand_G ( CS, C, BF )

    select type ( CS )
      class is ( Fluid_P_HN_Form )
    associate &
      ( I   =>  T % BoundaryIntegral % Integrand, &
        Cy  =>  C % Connectivity )
       
    do iFluence = 1, CS % nBalanced
      if ( CS % iaBalanced ( iFluence ) == CS % ELECTRON_DENSITY_B ) &
        iElectron = iFluence
    end do !-- iFluence

    do iF = 1, Cy % nFaces
      associate ( CE  =>  BF ( iElectron, iF ) % Value )
      do iS  =  1,  T % nSelected
        iI  =  T % iaSelected ( iS )
        if ( iI  ==  T % ELECTRON_NUMBER ) then
          call Copy ( CE, I ( iS, iF ) % Value )
        end if !-- iI
      end do !-- iS
      end associate !-- CE
    end do !-- iF

    end associate !-- I, etc.
    end select !-- CS

  end subroutine ComputeBoundaryIntegrand_G


end module Tally_F_P_HN__Form
