module ComputeGravity_Command

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies

  implicit none
  private

  public :: &
    ComputeGravity

contains


  subroutine ComputeGravity ( S )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S

    type ( StorageForm ) :: &
      GradPhi
    class ( Geometry_N_Form ), pointer :: &
      G
    class ( Fluid_D_Form ), pointer :: &
      F

    select type ( PS => S % Current_ASC % Atlas_SC )
    class is ( Atlas_SC_Form )

    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )
    G => GA % Geometry_N ( )

    select type ( FA => S % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_D ( )

    call GA % ComputeGravity &
           ( S % Current_ASC, &
             iBaryonMass = F % BARYON_MASS, &
             iBaryonDensity = F % COMOVING_BARYON_DENSITY )

    ! select type ( G )
    ! class is ( Geometry_N_S_Form )
    ! select type ( C => PS % Chart )
    ! class is ( Chart_SLD_Form )
    !   call GradPhi % Initialize &
    !          ( G, iaSelectedOption = [ G % POTENTIAL_GRADIENT_D ] )
    !   call C % ExchangeGhostData ( GradPhi )
    !   call PS % ApplyBoundaryConditionsFaces ( GradPhi )  
    ! end select !-- C
    ! end select !-- G

    end select !-- FA
    end select !-- GA
    end select !-- PS
    nullify ( F, G )

  end subroutine ComputeGravity


end module ComputeGravity_Command
