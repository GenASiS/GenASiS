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
  type ( Integrator_CSA_Form ), allocatable :: &
    I

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Integrator_CSA__Form_Test', DimensionalityOption = '2D' )

  allocate ( I )
  call I % Initialize ( )

  call SetParameters ( I )

  call I % Show ( )
  I % SetInitial  =>  SetInitial

  call I % Evolve ( )

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


end program Integrator_CSA__Form_Test
