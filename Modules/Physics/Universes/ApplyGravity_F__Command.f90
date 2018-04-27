module ApplyGravity_F__Command

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies

  implicit none
  private

  public :: &
    ApplyGravity_F

    private :: &
      ApplyGravityMomentum

contains


  subroutine ApplyGravity_F &
               ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_F
    type ( VariableGroupForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      iMomentum
    type ( TimerForm ), pointer :: &
      Timer
    class ( Geometry_N_Form ), pointer :: &
      G

    Timer => PROGRAM_HEADER % TimerPointer ( S % iTimerSources )
    if ( associated ( Timer ) ) call Timer % Start ( )

    select type ( PS => S % Current_ASC % Atlas_SC )
    class is ( Atlas_SC_Form )

    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )
    G  =>  GA % Geometry_N ( )

    select type ( F => Fluid )
    class is ( Fluid_D_Form )

    call GA % ComputeGravity &
           ( S % Current_ASC, &
             iBaryonMass = F % BARYON_MASS, &
             iBaryonDensity = F % COMOVING_BARYON_DENSITY )

    select type ( FS => Sources_F )
    class is ( Sources_F_Form )

    select type ( C => PS % Chart )
    class is ( Chart_SLD_Form )
      do iD = 1, C % nDimensions
        call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( iD ), &
                      iMomentum )
        call ApplyGravityMomentum &
               ( F % Value ( :, F % BARYON_MASS ), &
                 F % Value ( :, F % CONSERVED_BARYON_DENSITY ), &
                 G % Value ( :, G % GRAVITATIONAL_ACCELERATION_D ( iD ) ), &
                 TimeStep, S % B ( iStage ), &
                 Increment % Value ( :, iMomentum ), & 
                 FS % Value ( :, FS % GRAVITATIONAL_S_D ( iD ) ) )
      end do !-- iD
    end select !-- C

    end select !-- FS
    end select !-- F
    end select !-- GA
    end select !-- PS
    nullify ( G )

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ApplyGravity_F


  subroutine ApplyGravityMomentum ( M, N, GradPhi, dt, Weight_RK, K, S )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      GradPhi
    real ( KDR ) :: &
      dt, &
      Weight_RK
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      K, &
      S

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues

    nValues = size ( K )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      K ( iV )  =  K ( iV )  -  M ( iV )  *  N ( iV )  *  GradPhi ( iV )  &
                                *  dt
      S ( iV )  =  S ( iV )  -  M ( iV )  *  N ( iV )  *  GradPhi ( iV )  &
                                *  Weight_RK
    end do
    !$OMP end parallel do

  end subroutine ApplyGravityMomentum


end module ApplyGravity_F__Command
