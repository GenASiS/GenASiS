program Integrator_CS__Form_Test

  !-- Integrator_CurrentStream__Form_Test

  use Basics
  use Manifolds
  use Fields
  use Integrators

  implicit none

  integer ( KDI ) :: &
    iD, &
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
  type ( CurrentSetForm ), allocatable :: &
    CS_R, &  !-- Reference
    CS_D     !-- Difference
  type ( Integrator_CS_Form ), allocatable :: &
    I

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Integrator_CS__Form_Test', DimensionalityOption = '2D' )

  allocate ( I )
  call I % Initialize ( )
  I % SetInitial    =>  SetInitial
  I % SetReference  =>  SetReference

  associate ( CS  =>  I % CurrentSet_X )
  do iD  =  1, 3
    call CS % SetBoundaryConditionsFace &
           ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = iD )
  end do !-- iD
  end associate !-- CS

  call SetParameters ( I )

  allocate ( CS_R, CS_D )
  associate &
    ( G  =>  I % Geometry_X, &
      S  =>  I % Checkpoint_X )
  call CS_R % Initialize ( G, NameOption = 'Reference' )
  call CS_D % Initialize ( G, NameOption = 'Difference' )
  call CS_R % SetStream ( S )
  call CS_D % SetStream ( S )
  end associate !-- G, etc.

  call CS_R % Show ( )
  call CS_D % Show ( )

  call I % Evolve ( )

  deallocate ( CS_D, CS_R )
  deallocate ( I )
  deallocate ( PROGRAM_HEADER )

contains


  subroutine SetParameters ( I )

    class ( Integrator_CS_Form ), intent ( inout ) :: &
      I

    select type ( A  =>  I % X )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

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

    end associate !-- C
    end select    !-- A

  end subroutine SetParameters


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( I )
      class is ( Integrator_CS_Form )
    associate ( CS  =>  I % CurrentSet_X )

    call SetWave ( CS, T = 0.0_KDR )

    end associate !-- CS
    end select !-- I

  end subroutine SetInitial


  subroutine SetReference ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    call SetWave ( CS_R, I % T )

    select type ( I )
      class is ( Integrator_CS_Form )
    associate &
      ( CS  =>  I % CurrentSet_X )
    associate &
      ( CSV    =>  CS   % Storage_GS % Value, &
        CSV_R  =>  CS_R % Storage_GS % Value, &
        CSV_D  =>  CS_D % Storage_GS % Value )
    
    CSV_D  =  CSV  -  CSV_R

    call ComputeError ( CS, CS_R )

    end associate !-- CSV, etc.
    end associate !-- CS
    end select !-- I

  end subroutine SetReference


  subroutine SetWave ( CS, T )

    class ( CurrentSetForm ), intent ( inout ) :: &
      CS
    real ( KDR ), intent ( in ) :: &
      T

    associate ( G  =>  CS % Geometry )

    associate &
      (  GV  =>   G % Storage_GS % Value, &
        CSV  =>  CS % Storage_GS % Value )
    associate &
      (     X  =>   GV ( :,  G % CENTER_U_1 ), &
            Y  =>   GV ( :,  G % CENTER_U_2 ), &
            Z  =>   GV ( :,  G % CENTER_U_3 ), &
          Rho  =>  CSV ( :, CS % DENSITY_CS ), &
          V_1  =>  CSV ( :, CS % VELOCITY_CS_U_1 ), &
          V_2  =>  CSV ( :, CS % VELOCITY_CS_U_2 ), &
          V_3  =>  CSV ( :, CS % VELOCITY_CS_U_3 ), &
            K  =>  Wavenumber, &
        Abs_K  =>  sqrt ( dot_product ( Wavenumber, Wavenumber ) ), &
        TwoPi  =>  2.0_KDR  *  CONSTANT % PI )

    Rho  =  Offset  &
            +  Amplitude  &
               *  sin ( TwoPi * (    K ( 1 )  *  ( X  -  V_1  *  T ) &
                                  +  K ( 2 )  *  ( Y  -  V_2  *  T ) &
                                  +  K ( 3 )  *  ( Z  -  V_3  *  T ) ) )

    V_1  =  Speed  *  K ( 1 )  /  Abs_K
    V_2  =  Speed  *  K ( 2 )  /  Abs_K
    V_3  =  Speed  *  K ( 3 )  /  Abs_K
    
    end associate !-- X, etc.
    end associate !-- GV, etc.
    end associate !-- G

  end subroutine SetWave


  subroutine ComputeError ( FS, FS_R )

    class ( FieldSetForm ), intent ( in ) :: &
      FS, &
      FS_R

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF     !-- iField
    type ( CollectiveOperation_R_Form ) :: &
      CO 

    call Show ( 'Computing error' )
    call Show ( FS % Name, 'FieldSet' )

    select type ( A  =>  FS % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

    associate ( nF  =>  FS % nFields )

    associate ( Cr  =>  C % Communicator )
    call CO % Initialize &
           ( Cr, nOutgoing = [ 2 * nF ], nIncoming = [ 2 * nF ] )
    end associate !-- Cr

    do iS  =  1, nF
      iF  =  FS % iaSelected ( iS )
      associate &
        ( F_R  =>  FS_R % Storage_GS % Value ( :, iF ), &
          F    =>  FS   % Storage_GS % Value ( :, iF ) )

      !-- proper cells only
      CO % Outgoing % Value ( iS )  &
        =  sum ( pack ( abs ( F  -  F_R ), mask = C % ProperCell ) )
      CO % Outgoing % Value ( nF + iS )  &
        =  sum ( pack ( abs ( F_R ), mask = C % ProperCell ) )

!      !-- with ghost cells
!      CO % Outgoing % Value ( iF )  &
!        =  sum ( abs ( F  -  F_R ) )

      end associate !-- F_R, etc.
    end do !-- iS

    call CO % Reduce ( REDUCTION % SUM )
    associate &
      ( Norm_D  =>        CO % Incoming % Value (      1 :      nF ), &
        Norm_R  =>  max ( CO % Incoming % Value ( nF + 1 : nF + nF ), &
                          sqrt ( tiny ( 0.0_KDR ) ) ) )
    call Show ( Norm_D  /  Norm_R , 'L1 Error' )
    end associate !-- Norm_D, etc.

    end associate !-- nF
    end associate !-- C
    end select    !-- A

  end subroutine ComputeError


end program Integrator_CS__Form_Test
