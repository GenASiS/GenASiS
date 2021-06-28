program RiemannSolver_HLL_A__Form_Test

  !-- RiemannSolver_HartenLaxVanLeer_Atlas__Form_Test

  use Basics
  use Manifolds
  use Fields
  use Slopes

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_Form ), allocatable :: &
    A
  type ( Stream_A_Form ), allocatable :: &
    SA
  type ( Geometry_F_A_Form ), allocatable :: &
    GA
  type ( CurrentSet_A_Form ), allocatable :: &
    CSA
  type ( RiemannSolver_HLL_A_Form ), allocatable :: &
    RSA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'RiemannSolver_HLL_A__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           PeriodicOption = [ .true., .true., .true. ] )

  allocate ( SA )
  call SA % Initialize ( A, GIS )

  allocate ( GA )
  call GA % Initialize ( A )

  allocate ( CSA )
  call CSA % Initialize( GA )
  call CSA % SetStream ( SA )

  allocate ( RSA )
  call RSA % Initialize ( CSA )
  call SA % AddFieldSet ( RSA )

  call   A % Show ( )
  call CSA % Show ( )
  call RSA % Show ( )
  call  SA % Show ( )

  call SetWave ( CSA, GA )
  associate ( C  =>  A % Chart ( 1 ) % Element )
  call TestRiemannSolver ( RSA, iD = 1 )
  if ( C % nDimensions  >  1 ) &
    call TestRiemannSolver ( RSA, iD = 2 )
  if ( C % nDimensions  >  2 ) &
    call TestRiemannSolver ( RSA, iD = 3 )
  end associate !-- C

  deallocate ( RSA )
  deallocate ( CSA )
  deallocate ( GA )
  deallocate ( SA )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine SetWave ( CSA, GA )

    class ( CurrentSet_A_Form ), intent ( inout ) :: &
      CSA
    class ( Geometry_F_A_Form ), intent ( in ), target :: &
      GA

    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Offset, &
      Amplitude, &
      Speed
    real ( KDR ), dimension ( 3 ) :: &
      Wavenumber

    select type ( CSC  =>  CSA % FieldSet_C ( 1 ) % Element )
    class is ( CurrentSet_C_Form )

    select type ( GC  =>  GA % FieldSet_C ( 1 ) % Element )
    class is ( Geometry_F_C_Form )

    select type ( C  =>  CSC % Chart )
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
      (     X  =>   GC % Storage_FSC % Storage &
                       % Value ( :, GC % CENTER_U_1 ), &
            Y  =>   GC % Storage_FSC % Storage &
                       % Value ( :, GC % CENTER_U_2 ), &
            Z  =>   GC % Storage_FSC % Storage &
                       % Value ( :, GC % CENTER_U_3 ), &
          Rho  =>  CSC % Storage_FSC % Storage &
                       % Value ( :, CSC % DENSITY_DEFAULT ), &
            K  =>  Wavenumber, &
        Abs_K  =>  sqrt ( dot_product ( Wavenumber, Wavenumber ) ), &
        TwoPi  =>  2.0_KDR  *  CONSTANT % PI )

    Rho  =  Offset  &
            +  Amplitude  &
               *  sin ( TwoPi * (    K ( 1 ) * X  &
                                  +  K ( 2 ) * Y  &
                                  +  K ( 3 ) * Z  ) )

    call CSC % SetVelocityConstant ( Wavenumber, Speed )
    
    end associate !-- Rho, etc.
    end select !-- C
    end select !-- GC
    end select !-- CSC

  end subroutine SetWave


  subroutine TestRiemannSolver ( RSA, iD )

    class ( RiemannSolver_HLL_A_Form ), intent ( inout ) :: &
      RSA
    integer ( KDI ), intent ( in ) :: &
      iD

    call RSA % Compute ( iD )

    call GIS % Open ( GIS % ACCESS_CREATE )
    call SA % Write ( )
    call GIS % Close ( )

  end subroutine TestRiemannSolver


end program RiemannSolver_HLL_A__Form_Test
