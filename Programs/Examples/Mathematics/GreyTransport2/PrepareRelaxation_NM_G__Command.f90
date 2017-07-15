module PrepareRelaxation_NM_G__Command

  use Basics
  use Mathematics
!   use RadiationMoments_Form
!   use Sources_RM__Form
  use NeutrinoMoments_G__Form
  use PrepareRelaxation_RM__Command

  implicit none
  private

  public :: &
    PrepareRelaxation_NM_G

    private :: &
      PrepareRelaxation_NM_G_Kernel

contains


  subroutine PrepareRelaxation_NM_G &
               ( S, Sources_NM_G, IncrementExplicit, DampingCoefficient, &
                 NeutrinoMoments_G, Chart, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_NM_G
    type ( VariableGroupForm ), intent ( inout ) :: &
      IncrementExplicit, &
      DampingCoefficient
    class ( CurrentTemplate ), intent ( in ) :: &
      NeutrinoMoments_G
    class ( ChartTemplate ), intent ( in ) :: &
      Chart
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    integer ( KDI ) :: &
      iNumber

    call PrepareRelaxation_RM &
           ( S, Sources_NM_G, IncrementExplicit, DampingCoefficient, &
             NeutrinoMoments_G, Chart, TimeStep, iStage)

    select type ( NM => NeutrinoMoments_G )
    class is ( NeutrinoMoments_G_Form )

    call Search ( NM % iaConserved, NM % CONSERVED_NUMBER, &
                  iNumber )

    select type ( Chart )
    class is ( Chart_SL_Template )

    associate ( I => NM % Interactions )

!     select type ( SNM => Sources_NM )
!     class is ( Sources_NM_Form )
      call PrepareRelaxation_NM_G_Kernel &
             ( IncrementExplicit % Value ( :, iNumber ), &
               DampingCoefficient % Value ( :, iNumber ), &
!                SNM % Value ( :, SNM % NET_EMISSION_E ), &
!                SNM % Value ( :, SNM % NET_EMISSION_S_D ( 1 ) ), &
!                SNM % Value ( :, SNM % NET_EMISSION_S_D ( 2 ) ), &
!                SNM % Value ( :, SNM % NET_EMISSION_S_D ( 3 ) ), &
               Chart % IsProperCell, &
               I % Value ( :, I % EMISSIVITY_N ), &
               I % Value ( :, I % OPACITY_N ), &
               NM % Value ( :, NM % COMOVING_NUMBER ), &
               NM % Value ( :, NM % CONSERVED_NUMBER ), &
               TimeStep, S % B ( iStage ), CONSTANT % SPEED_OF_LIGHT )
!     end select !-- SNM

    end associate !-- I

    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'PrepareRelaxation_NM_G__Command', 'module', CONSOLE % ERROR )
      call Show ( 'PrepareRelaxation_NM_G', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Chart

    end select !-- NM
    
  end subroutine PrepareRelaxation_NM_G
    

  subroutine PrepareRelaxation_NM_G_Kernel &
               ( KV_N, DCV_N, &
!                 SVNE_E, SVNE_S_1, SVNE_S_2, SVNE_S_3, &
                 IsProperCell, Xi_N, Chi_N, N, D, dT, Weight_RK, c )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      KV_N, &
      DCV_N
!      SVNE_E, SVNE_S_1, SVNE_S_2, SVNE_S_3
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Xi_N, Chi_N, &
      N, D
    real ( KDR ), intent ( in ) :: &
      dT, &
      Weight_RK, &
      c

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( KV_N )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( .not. IsProperCell ( iV ) ) &
        cycle

      KV_N ( iV )  &
        =  KV_N ( iV )  &
           +  c * dT  &
              * ( Xi_N ( iV )  &
                  -  Chi_N ( iV )  *  ( N ( iV )  -  D ( iV ) ) )

      DCV_N ( iV )  =  c * Chi_N ( iV )

    end do
    !$OMP end parallel do

!     !$OMP parallel do private ( iV )
!     do iV = 1, nV
!       if ( .not. IsProperCell ( iV ) ) &
!         cycle

!       !-- Approximate since E, S_1, etc. do not include implicit update

!       SVNE_E ( iV )  &
!         =  SVNE_E ( iV )  &
!            +  Weight_RK * c * ( Xi_J ( iV )  -  Chi_J ( iV ) * E ( iV ) )

!       SVNE_S_1 ( iV )  &
!         =  SVNE_S_1 ( iV )  -  Weight_RK * c * Chi_H ( iV ) * S_1 ( iV )
!       SVNE_S_2 ( iV )  &
!         =  SVNE_S_2 ( iV )  -  Weight_RK * c * Chi_H ( iV ) * S_2 ( iV )
!       SVNE_S_2 ( iV )  &
!         =  SVNE_S_2 ( iV )  -  Weight_RK * c * Chi_H ( iV ) * S_3 ( iV )

!    end do
!    !$OMP end parallel do

  end subroutine PrepareRelaxation_NM_G_Kernel
  

end module PrepareRelaxation_NM_G__Command
