module LinearAdvection_Form

  use Basics
  use Mathematics

  implicit none
  private

  type, public :: LinearAdvectionForm
    real ( KDR ), private :: &
      Speed, &
      Length, &
      Density, &
      DimensionFactor
    character ( LDL ), private :: &
      CoordinateSystem, &
      AdvectionType
    type ( Integrator_CSA_Form ), allocatable :: &
      Integrator
    type ( CurrentSet_A_Form ), allocatable :: &
      Reference_A, &
      Difference_A
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Evolve
    final :: &
      Finalize
  end type LinearAdvectionForm

  class ( LinearAdvectionForm ), public, pointer :: &
    LINEAR_ADVECTION => null ( )  !-- Makes instance of LinearAdvection
                                  !   accessible to SetInitial and SetReference

    private :: &
      SetParameters, &
      SetInitial, &
      SetReference

      private :: &
        ComputeError


contains


  subroutine Initialize ( LA, CoordinateSystem, AdvectionType )

    class ( LinearAdvectionForm ), intent ( inout ), target :: &
      LA
    character ( * ), intent ( in ) :: &
      CoordinateSystem, &
      AdvectionType

    call Show ( 'Initializing a LinearAdvection' )

    LINEAR_ADVECTION  =>  LA

    LA % CoordinateSystem  =  CoordinateSystem
    LA % AdvectionType     =  AdvectionType

    allocate ( LA % Integrator )
    associate ( I  =>  LA % Integrator )

    call I % Initialize ( )
    I % SetInitial    =>  SetInitial
!    I % SetReference  =>  SetReference

    associate ( CSA  =>  I % CurrentSet_X_A )
    call CSA % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'OUTFLOW   ' ], iDimension = 1 )
    select case ( trim ( CoordinateSystem ) )
    case ( 'SPHERICAL' )
      call CSA % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 2 )
    end select !-- CoordinateSystem
    end associate !-- CSA

    call SetParameters ( LA )

    end associate !-- Integrator

  end subroutine Initialize


  subroutine Evolve ( LA )

    class ( LinearAdvectionForm ), intent ( inout ) :: &
      LA

    associate ( I  =>  LA % Integrator )
    call I % Evolve ( )
    end associate !-- I
   
  end subroutine Evolve


  impure elemental subroutine Finalize ( LA )

    type ( LinearAdvectionForm ), intent ( inout ) :: &
      LA

    if ( allocated ( LA % Difference_A ) ) &
      deallocate ( LA % Difference_A )
    if ( allocated ( LA % Reference_A ) ) &
      deallocate ( LA % Reference_A )
    if ( allocated ( LA % Integrator ) ) &
      deallocate ( LA % Integrator )

    call Show ( 'Finalizing a LinearAdvection' )

  end subroutine Finalize


  subroutine SetParameters ( LA )

    class ( LinearAdvectionForm ), intent ( inout ) :: &
      LA

    associate &
      ( I  =>  LA % Integrator )
    select type ( A  =>  I % X_A )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )
    associate &
      ( V_0  =>  LINEAR_ADVECTION % Speed, &
        L_0  =>  LINEAR_ADVECTION % Length, &
        N_0  =>  LINEAR_ADVECTION % Density, &
        D    =>  LINEAR_ADVECTION % DimensionFactor )

    select case ( trim ( LINEAR_ADVECTION % CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      D  =  1.0_KDR
    case ( 'SPHERICAL' )
      D  =  3.0_KDR
    end select

    select case ( trim ( LINEAR_ADVECTION % AdvectionType ) )
    case ( 'CONTRACTION' )

      V_0  =  -1.0_KDR
      L_0  =   C % MaxCoordinate ( 1 )  *  0.8_KDR
      call PROGRAM_HEADER % GetParameter ( V_0, 'Speed'  )
      call PROGRAM_HEADER % GetParameter ( L_0, 'Length' )

      I % T_Finish  =  log ( 0.25_KDR ** D )  *  L_0 / ( D * V_0 )

    case ( 'EXPANSION' )

      V_0  =  +1.0_KDR
      L_0  =   C % MaxCoordinate ( 1 )  *  0.2_KDR
      call PROGRAM_HEADER % GetParameter ( V_0, 'Speed'  )
      call PROGRAM_HEADER % GetParameter ( L_0, 'Length' )

      I % T_Finish  =  log ( 4.0_KDR ** D )  *  L_0 / ( D * V_0 )

   end select

    N_0  =  1.0_KDR
    call PROGRAM_HEADER % GetParameter ( N_0, 'Density' )

    call I % Show ( )
    call Show ( 'LinearAdvection Parameters' )
    call Show ( LA % CoordinateSystem, 'CoordinateSystem' )
    call Show ( LA % AdvectionType,    'AdvectionType' )
    call Show ( V_0, 'Speed' )
    call Show ( L_0, 'Length' )
    call Show ( N_0, 'Density' )

    end associate !-- V_0, etc.
    end associate !-- C
    end select !-- A
    end associate !-- I

  end subroutine SetParameters


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( I )
      class is ( Integrator_CSA_Form )
    associate &
      ( CSA  =>  I % CurrentSet_X_A )

    call SetDensity ( CSA, T = 0.0_KDR )

    end associate !-- CSA
    end select !-- I

  end subroutine SetInitial


  subroutine SetReference ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( I )
      class is ( Integrator_CSA_Form )
    associate &
      ( CSA    =>  I % CurrentSet_X_A, &
        CSA_R  =>  LINEAR_ADVECTION % Reference_A, &
        CSA_D  =>  LINEAR_ADVECTION % Difference_A )
    associate &
      ( CSC    =>  CSA   % FieldSet_C ( 1 ) % Element, &
        CSC_R  =>  CSA_R % FieldSet_C ( 1 ) % Element, &
        CSC_D  =>  CSA_D % FieldSet_C ( 1 ) % Element )
    associate &
      ( CSV    =>  CSC   % Storage_FSC % Storage % Value, &
        CSV_R  =>  CSC_R % Storage_FSC % Storage % Value, &
        CSV_D  =>  CSC_D % Storage_FSC % Storage % Value )
    
    call SetDensity ( CSA_R, I % T )

    CSV_D  =  CSV  -  CSV_R

    call ComputeError ( CSC, CSC_R )

    end associate !-- CSV, etc.
    end associate !-- CSC, etc.
    end associate !-- CSA, etc.
    end select !-- I

  end subroutine SetReference


  subroutine SetDensity ( CSA, T )

    class ( CurrentSet_A_Form ), intent ( inout ) :: &
      CSA
    real ( KDR ), intent ( in ) :: &
      T

    real ( KDR ) :: &
      L_T, &
      N_T, &
      E_T, &
      Dim   

    associate &
      ( GA  =>  CSA % Geometry_A )
    select type ( GC  =>  GA % FieldSet_C ( 1 ) % Element )
      class is ( Geometry_F_C_Form )
    select type ( CSC  =>  CSA % FieldSet_C ( 1 ) % Element )
      class is ( CurrentSet_C_Form )
    associate &
      (  GV  =>   GC % Storage_FSC % Storage % Value, &
        CSV  =>  CSC % Storage_FSC % Storage % Value )
    associate &
      ( X_1  =>   GV ( :,  GC % CENTER_U_1 ), &
        N    =>  CSV ( :, CSC % DENSITY_DEFAULT ) )
    associate &
      ( V_0  =>  LINEAR_ADVECTION % Speed, &
        L_0  =>  LINEAR_ADVECTION % Length, &
        N_0  =>  LINEAR_ADVECTION % Density, &
        D    =>  LINEAR_ADVECTION % DimensionFactor )

    call CSC % SetVelocityLinear ( V_0, L_0 )

    E_T  =  exp ( - V_0 / L_0  *  T )

    N_T  =  N_0  *  E_T ** D
    L_T  =  L_0  *  E_T ** ( - 1.0_KDR / D )

    where ( X_1  <=  L_T )
      N  =  N_T
    elsewhere
      N  =  0.0_KDR
    end where

    end associate !-- V_0, etc.
    end associate !-- X_1, etc.
    end associate !-- GV, etc.
    end select !-- CSC
    end select !-- GC
    end associate !-- GA

  end subroutine SetDensity


  subroutine ComputeError ( FSC, FSC_R )

    class ( FieldSet_C_Form ), intent ( in ) :: &
      FSC, &
      FSC_R

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF     !-- iField
    type ( CollectiveOperation_R_Form ) :: &
      CO 

    call Show ( 'Computing error' )
    call Show ( FSC % Name, 'FieldSet' )

    associate ( nF  =>  FSC % nFields )

    select type ( C  =>  FSC % Chart )
    class is ( Chart_GS_Form )

    associate ( Cmm  =>  C % Communicator )
    call CO % Initialize &
           ( Cmm, nOutgoing = [ 2 * nF ], nIncoming = [ 2 * nF ] )
    end associate !-- Cmm

    do iS  =  1, nF
      iF  =  FSC % iaSelected ( iS )
      associate &
        ( F_R  =>  FSC_R % Storage_FSC % Storage % Value ( :, iF ), &
          F    =>  FSC   % Storage_FSC % Storage % Value ( :, iF ) )

      !-- proper cells only
      CO % Outgoing % Value ( iS )  &
        =  sum ( pack ( abs ( F  -  F_R ), mask = C % ProperCell ) )
      CO % Outgoing % Value ( nF + iS )  &
        =  sum ( pack ( abs ( F_R ), mask = C % ProperCell ) )

!      !-- with ghost cells
!      CO % Outgoing % Value ( iF )  &
!        =  sum ( abs ( F  -  F_R ) )

      end associate !-- F_R, etc.
    end do !-- iF
    end select !-- G

    call CO % Reduce ( REDUCTION % SUM )
    associate &
      ( Norm_D  =>        CO % Incoming % Value (      1 :      nF ), &
        Norm_R  =>  max ( CO % Incoming % Value ( nF + 1 : nF + nF ), &
                          sqrt ( tiny ( 0.0_KDR ) ) ) )
    call Show ( Norm_D  /  Norm_R , 'L1 Error' )
    end associate !-- Norm_D, etc.

    end associate !-- nF

  end subroutine ComputeError


end module LinearAdvection_Form
