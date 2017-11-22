module Tally_RT__Form

  !-- Tally_RayleighTaylor_Form

  use Basics
  use Mathematics
  use Fluid_P_P__Form
  use Tally_F_P__Form
  
  implicit none
  private

  type, public, extends ( Tally_F_P_Form ) :: Tally_RT_Form
     real ( KDR ) :: &
       Acceleration
  contains
    procedure, public, pass :: &
      SetAcceleration
    procedure, public, pass :: &
      SelectVariables
    procedure, public, pass :: &
      ComputeInteriorIntegrand
  end type Tally_RT_Form

contains


  subroutine SetAcceleration ( T, Acceleration )

    class ( Tally_RT_Form ), intent ( inout ) :: &
      T
    real ( KDR ), intent ( in ) :: &
      Acceleration

    T % Acceleration = Acceleration
    
  end subroutine SetAcceleration


  subroutine SelectVariables ( T, A ) 
    
    class ( Tally_RT_Form ), intent ( inout ) :: &
      T
    class ( AtlasHeaderForm ), intent ( in ) :: &
      A

    if ( allocated ( T % iaSelected ) ) &
      deallocate ( T % iaSelected )

    T % nSelected = 12

    allocate ( T % iaSelected ( T % nSelected ) )

    T % iaSelected &
      = [ T % BARYON_NUMBER, T % MOMENTUM, T % FLUID_ENERGY, &
          T % INTERNAL_ENERGY, T % KINETIC_ENERGY, T % ANGULAR_MOMENTUM, &
          T % GRAVITATIONAL_ENERGY, T % TOTAL_ENERGY ]
    
  end subroutine SelectVariables
  

  subroutine ComputeInteriorIntegrand &
               ( T, Integrand, C, G, nDimensions )

    class ( Tally_RT_Form ), intent ( inout ) :: &
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
    class is ( Fluid_P_P_Form )

    call T % Tally_F_P_Form &
           % ComputeInteriorIntegrand ( Integrand, C, G, nDimensions )

    associate &
      ( CE => C % Value ( :, C % CONSERVED_ENERGY ), &
        N  => C % Value ( :, C % COMOVING_DENSITY ), &
        Y  => G % Value ( :, G % CENTER ( 2 ) ), &
        A  => T % Acceleration )

    do iS = 1, T % nSelected
      iI = T % iaSelected ( iS )
      if ( iI == T % GRAVITATIONAL_ENERGY ) then
        Integrand ( iS ) % Value = N * A * Y
      else if ( iI == T % TOTAL_ENERGY ) then
        Integrand ( iS ) % Value = CE  +  N * A * Y
      end if !-- iI
    end do !-- iS

    end associate !-- CE, etc.
    end select !-- F

  end subroutine ComputeInteriorIntegrand


end module Tally_RT__Form
