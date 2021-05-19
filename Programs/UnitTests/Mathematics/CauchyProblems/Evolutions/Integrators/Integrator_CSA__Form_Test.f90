program Integrator_CSA__Form_Test

  !-- Integrator_CurrentStreamAtlas__Form_Test

  use Basics
  use Manifolds
  use Fields
  use Integrators

  implicit none

  integer ( KDI ) :: &
    nPeriods
  integer ( KDI ), dimension ( 3 ) :: &
    nWavelengths
  real ( KDR ) :: &
    Offset, &
    Amplitude, &
    Speed, &
    Period
  real ( KDR ), dimension ( 3 ) :: &
    Wavenumber
  type ( CurrentSet_A_Form ), allocatable :: &
    CSA_R, &  !-- Reference
    CSA_D     !-- Difference
  type ( Integrator_CSA_Form ), allocatable :: &
    I

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Integrator_CSA__Form_Test', DimensionalityOption = '2D' )

  allocate ( I )
  call I % Initialize ( )
  I % SetInitial    =>  SetInitial
  I % SetReference  =>  SetReference

  call SetParameters ( I )

  allocate ( CSA_R, CSA_D )
  associate &
    ( GA  =>  I % Geometry_X_A, &
      SA  =>  I % Checkpoint_X_A )
  call CSA_R % Initialize ( GA, NameOption = 'Reference' )
  call CSA_D % Initialize ( GA, NameOption = 'Difference' )
  call CSA_R % SetStream ( SA )
  call CSA_D % SetStream ( SA )
  end associate !-- GA, etc.

  call I % Show ( )
  call CSA_R % Show ( )
  call CSA_D % Show ( )

  call I % Evolve ( )

  deallocate ( CSA_D, CSA_R )
  deallocate ( I )
  deallocate ( PROGRAM_HEADER )

contains


  subroutine SetParameters ( I )

    class ( Integrator_CSA_Form ), intent ( inout ) :: &
      I

    select type ( C  =>  I % X_A % Chart ( 1 ) % Element )
      class is ( Chart_GS_Form )

    nWavelengths  =  0
    nWavelengths ( 1 : C % nDimensions )  =  1
    call PROGRAM_HEADER % GetParameter ( nWavelengths, 'nWavelengths' )

    Offset     =  2.0_KDR
    Amplitude  =  1.0_KDR
    Speed      =  1.0_KDR
    call PROGRAM_HEADER % GetParameter ( Offset, 'Offset' )
    call PROGRAM_HEADER % GetParameter ( Amplitude, 'Amplitude' )
    call PROGRAM_HEADER % GetParameter ( Speed, 'Speed' )

    associate ( BoxSize  =>  C % MaxCoordinate  -  C % MinCoordinate )
    where ( BoxSize  >  0.0_KDR )
      Wavenumber  =  nWavelengths / BoxSize
    elsewhere
      Wavenumber  =  0.0_KDR
    end where
    end associate !-- BoxSize

    associate &
      ( K      =>  Wavenumber, &
        Abs_K  =>  sqrt ( dot_product ( Wavenumber, Wavenumber ) ), &
        V      =>  Speed )
    Period  =  1.0_KDR / ( Abs_K * V )
    call Show ( Period, 'Period' )
    end associate !-- K, etc.

    nPeriods   =  1
    call PROGRAM_HEADER % GetParameter ( nPeriods, 'nPeriods' )

    I % T_Finish  =  nPeriods * Period

    call Show ( 'Test Parameters' )
    call Show ( Offset, 'Offset' )
    call Show ( Amplitude, 'Amplitude' )
    call Show ( Speed, 'Speed' )
    call Show ( Period, 'Period' )
    call Show ( nPeriods, 'nPeriods' )

    end select !--   C

  end subroutine SetParameters


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( I )
      class is ( Integrator_CSA_Form )
    associate ( CSA  =>  I % CurrentSet_X_A )

    call SetWave ( CSA, T = 0.0_KDR )

    end associate !-- CSA
    end select !-- I

  end subroutine SetInitial


  subroutine SetReference ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    call SetWave ( CSA_R, I % T )

    select type ( I )
      class is ( Integrator_CSA_Form )
    associate &
      ( CSA  =>  I % CurrentSet_X_A )
    associate &
      ( CSC    =>  CSA   % FieldSet_C ( 1 ) % Element, &
        CSC_R  =>  CSA_R % FieldSet_C ( 1 ) % Element, &
        CSC_D  =>  CSA_D % FieldSet_C ( 1 ) % Element )
    associate &
      ( CSV    =>  CSC   % Storage_FSC % Storage % Value, &
        CSV_R  =>  CSC_R % Storage_FSC % Storage % Value, &
        CSV_D  =>  CSC_D % Storage_FSC % Storage % Value )
    
    CSV_D  =  CSV  -  CSV_R

    call ComputeError ( CSC, CSC_R )

    end associate !-- CSV, etc.
    end associate !-- CSC, etc.
    end associate !-- CSA
    end select !-- I

  end subroutine SetReference


  subroutine SetWave ( CSA, T )

    class ( CurrentSet_A_Form ), intent ( inout ) :: &
      CSA
    real ( KDR ), intent ( in ) :: &
      T

    associate ( GA  =>  CSA % Geometry_A )
    select type ( CSC  =>  CSA % FieldSet_C ( 1 ) % Element )
      class is ( CurrentSet_C_Form )
    select type ( GC  =>  GA % FieldSet_C ( 1 ) % Element )
      class is ( Geometry_F_C_Form )

    associate &
      (     X  =>   GC % Storage_FSC % Storage &
                       % Value ( :, GC % CENTER_U_1 ), &
            Y  =>   GC % Storage_FSC % Storage &
                       % Value ( :, GC % CENTER_U_2 ), &
            Z  =>   GC % Storage_FSC % Storage &
                       % Value ( :, GC % CENTER_U_3 ), &
          Rho  =>  CSC % Storage_FSC % Storage &
                       % Value ( :, CSC % DENSITY_DEFAULT ), &
            V  =>  CSC % VelocityDefault_U, &
            K  =>  Wavenumber, &
        Abs_K  =>  sqrt ( dot_product ( Wavenumber, Wavenumber ) ), &
        TwoPi  =>  2.0_KDR  *  CONSTANT % PI )

    V ( 1 )  =  Speed  *  K ( 1 )  /  Abs_K
    V ( 2 )  =  Speed  *  K ( 2 )  /  Abs_K
    V ( 3 )  =  Speed  *  K ( 3 )  /  Abs_K
    
    Rho  =  Offset  &
            +  Amplitude  &
               *  sin ( TwoPi * (    K ( 1 )  *  ( X  -  V ( 1 )  *  T ) &
                                  +  K ( 2 )  *  ( Y  -  V ( 2 )  *  T ) &
                                  +  K ( 3 )  *  ( Z  -  V ( 3 )  *  T ) ) )

    end associate !-- Rho, etc.
    end select !-- GC
    end select !-- CSC
    end associate !-- GA

  end subroutine SetWave


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
    call Show (    CO % Incoming % Value (      1 :      nF )  &
                /  CO % Incoming % Value ( nF + 1 : nF + nF ), 'L1 Error' )

    end associate !-- nF

  end subroutine ComputeError


end program Integrator_CSA__Form_Test
