module Slope_DFV_F_F_P_HN__Form

  !-- Slope_DivergenceFiniteVolume_Flat_Fluid_Perfect_HeavyNucleus__Form

  use Basics
  use Mathematics
  use Fluid_P_HN__Form

  implicit none
  private

  type, public, extends ( Slope_DFV_F_Form ) :: Slope_DFV_F_F_P_HN_Form
    real ( KDR ) :: &
      GammaHigh, &
      GammaLow
  contains
    procedure, private, pass :: &
      InitializeAllocate_F
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_DFV_F_F_P_HN_Form


contains


  subroutine InitializeAllocate_F &
               ( S, RS, DP_1D, SuffixOption, IgnorabilityOption )

    class ( Slope_DFV_F_F_P_HN_Form ), intent ( inout ) :: &
      S
    class ( RiemannSolver_HLL_Form ), intent ( in ), target :: &
      RS
    type ( DivergencePartElement ), dimension ( : ), intent ( in ) :: &
      DP_1D
    character ( * ), intent ( in ), optional :: &
      SuffixOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_DFV_F_F_P_HN'

    call S % Slope_DFV_F_Form % Initialize &
           ( RS, DP_1D, SuffixOption, IgnorabilityOption )

    S % GammaHigh  =  2.2_KDR
    S % GammaLow   =  1.3_KDR
    call PROGRAM_HEADER % GetParameter ( S % GammaHigh, 'GammaHigh' )
    call PROGRAM_HEADER % GetParameter ( S % GammaLow,  'GammaLow' )

  end subroutine InitializeAllocate_F


  subroutine Compute ( S, T_Option, iS_Option )

    class ( Slope_DFV_F_F_P_HN_Form ), intent ( inout ) :: &
      S
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option
    integer ( KDI ), intent ( in ), optional :: &
      iS_Option

    integer ( KDI ) :: &
      iR_1, iR_2, iR_3, &
      iR, &
      iTh, &
      iPh, &
      iMomentum_1
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      Gamma, &
      P, &
      PF, &  !-- Pressure force, -dPdR
      R

    call S % Slope_DFV_F_Form % Compute ( T_Option, iS_Option )

    select type ( A  =>  S % Atlas )
    class is ( Atlas_SCG_CC_Form )

    select type ( F  =>  S % RiemannSolver % CurrentSet )
      class is ( Fluid_P_HN_Form )
    associate &
      ( SP  =>  S % Component ( 2 ) % Element, &  !-- Pressure gradient 
        G   =>  F % Geometry )
    associate &
      (  C   =>  A % Chart_GS_CC, &
         FV  =>  F % Storage_GS % Value, &
        SPV  =>  SP % Storage_GS % Value, &
         GV  =>  G % Storage_GS % Value )
    
    call Search ( F % iaBalanced, F % MOMENTUM_DENSITY_D_1, iMomentum_1 )

    call C % SetFieldPointer (  FV ( :, F % ADIABATIC_INDEX ), Gamma )
    call C % SetFieldPointer (  FV ( :, F % PRESSURE ), P )
    call C % SetFieldPointer ( SPV ( :, iMomentum_1 ), PF )
    call C % SetFieldPointer (  GV ( :, G % CENTER_U_1 ), R )

    do iPh  =  1,  C % nCells ( 3 )
      do iTh  =  1,  C % nCells ( 2 )
        if ( Gamma ( 2, iTh, iPh )  >  S % GammaHigh ) then

          iR_1  =  0
          iR_2  =  0
          iR_3  =  0
          do iR  =  2,  C % nCellsBrick ( 1 )
            if ( iR_1  ==  0 ) then
              if ( Gamma ( iR, iTh, iPh )  <  S % GammaHigh ) &
                iR_1  =  iR
            else if ( iR_2  ==  0 ) then
              if ( Gamma ( iR, iTh, iPh )  <  S % GammaLow ) &
                iR_2  =  iR
            else
              if ( Gamma ( iR, iTh, iPh )  >  S % GammaLow ) then
                iR_3  =  iR
                exit
              end if
            end if
          end do !-- iR

!call Show ( [ iR_1, iR_2, iR_3 ], '>>> iR_1,2,3' )

          do iR  =  iR_1, iR_3
!call Show ( iR, '>>> iR' )
!call Show ( PF ( iR, iTh, iPh ), '>>> PF old' )
            PF ( iR, iTh, iPh )  &
              =  -    ( P ( iR + 1, iTh, iPh )  -  P ( iR - 1, iTh, iPh ) )  &
                   /  ( R ( iR + 1, iTh, iPh )  -  R ( iR - 1, iTh, iPh ) )
!call Show ( PF ( iR, iTh, iPh ), '>>> PF new' )
          end do !--iR

        end if !-- Bulk nuclear matter
      end do !-- iTh
    end do !-- iPh

    end associate !-- C, etc.
    end associate !-- G, etc.
    end select !--F
    
    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Slope_DFV_F_F_P_HN__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    call S % Slope_DFV_F_Form % AddComponents ( T_Option = T_Option )

  end subroutine Compute


  impure elemental subroutine Finalize ( S )

    type ( Slope_DFV_F_F_P_HN_Form ), intent ( inout ) :: &
      S

  end subroutine Finalize


end module Slope_DFV_F_F_P_HN__Form
