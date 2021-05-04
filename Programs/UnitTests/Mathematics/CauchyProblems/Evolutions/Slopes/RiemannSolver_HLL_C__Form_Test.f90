program RiemannSolver_HLL_C__Form_Test

  !-- RiemannSolver_HartenLaxVanLeer_Chart__Form_Test

  use Basics
  use Manifolds
  use Fields
  use Slopes

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Chart_GS_Form ), allocatable :: &
    C
  type ( FieldSet_C_Form ), allocatable :: &
    BC, &
    OB_IL_C, OB_IR_C, &
    OF_IL_C, OF_IR_C, &
    OE_IL_C, OE_IR_C
  type ( Stream_C_Form ), allocatable :: &
    SC
  type ( Geometry_F_C_Form ), allocatable :: &
    GC
  type ( CurrentSet_C_Form ), allocatable :: &
    CSC
  type ( FluxSet_C_Form ), allocatable :: &
    FSC
  type ( Eigenspeeds_F_C_Form ), allocatable :: &
    EC
  type ( Reconstruction_C_Form ), allocatable :: &
    RBC, RFC, REC
  type ( RiemannSolver_HLL_C_Form ), allocatable :: &
    RSC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'RiemannSolver_HLL_C__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( C )
  call C % Initialize &
         ( PeriodicOption = [ .true., .true., .true. ] )

  allocate ( SC )
  call SC % Initialize ( C, GIS )

  allocate ( GC )
  call GC % Initialize ( C )

  allocate ( CSC )
  call CSC % Initialize ( GC ) 
  call CSC % SetStream ( SC )

  allocate ( BC )
  call BC % Initialize ( CSC, CSC % iaBalanced, NameOption = 'Balanced' )

  allocate ( FSC )
  call FSC % Initialize ( CSC )

  allocate ( EC )
  call EC % Initialize ( CSC )

  allocate ( OB_IL_C, OB_IR_C )
  call OB_IL_C % Initialize &
         ( C, &
           NameOption = 'R_' // trim ( BC % Name ) // '_IL', &
           nFieldsOption = BC % nFields )
  call OB_IR_C % Initialize &
         ( C, &
           NameOption = 'R_' // trim ( BC % Name ) // '_IR', &
           nFieldsOption = BC % nFields )

  allocate ( OF_IL_C, OF_IR_C )
  call OF_IL_C % Initialize &
         ( C, &
           NameOption = 'R_' // trim ( FSC % Name ) // '_IL', &
           nFieldsOption = FSC % nFields )
  call OF_IR_C % Initialize &
         ( C, &
           NameOption = 'R_' // trim ( FSC % Name ) // '_IR', &
           nFieldsOption = FSC % nFields )

  allocate ( OE_IL_C, OE_IR_C )
  call OE_IL_C % Initialize &
         ( C, &
           NameOption = 'R_' // trim ( EC % Name ) // '_IL', &
           nFieldsOption = EC % nFields )
  call OE_IR_C % Initialize &
         ( C, &
           NameOption = 'R_' // trim ( EC % Name ) // '_IR', &
           nFieldsOption = EC % nFields )

  allocate ( RBC, RFC, REC )
  call RBC % Initialize ( GC,  BC, OB_IL_C, OB_IR_C )
  call RFC % Initialize ( GC, FSC, OF_IL_C, OF_IR_C )
  call REC % Initialize ( GC,  EC, OE_IL_C, OE_IR_C )

  allocate ( RSC )
  call RSC % Initialize ( RBC, RFC, REC, EC, FSC, CSC ) 
  call SC % AddFieldSet ( RSC )

  call   C % Show ( )
  call CSC % Show ( )
  call RSC % Show ( )
  call  SC % Show ( )

  call SetWave ( CSC, GC )
  call TestRiemannSolver ( RSC, iD = 1 )
  if ( C % nDimensions  >  1 ) &
    call TestRiemannSolver ( RSC, iD = 2 )
  if ( C % nDimensions  >  2 ) &
    call TestRiemannSolver ( RSC, iD = 3 )

  deallocate ( RSC )
  deallocate ( REC, RFC, RBC )
  deallocate ( OE_IR_C, OE_IL_C )
  deallocate ( OF_IR_C, OF_IL_C )
  deallocate ( OB_IR_C, OB_IL_C )
  deallocate ( EC )
  deallocate ( FSC )
  deallocate ( BC )
  deallocate ( CSC )
  deallocate ( GC )
  deallocate ( SC )
  deallocate ( C )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine SetWave ( CSC, GC )

    class ( CurrentSet_C_Form ), intent ( inout ) :: &
      CSC
    class ( Geometry_F_C_Form ), intent ( in ), target :: &
      GC

    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Offset, &
      Amplitude, &
      Speed
    real ( KDR ), dimension ( 3 ) :: &
      Wavenumber

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
            V  =>  CSC % VelocityDefault_U, &
            K  =>  Wavenumber, &
        Abs_K  =>  sqrt ( dot_product ( Wavenumber, Wavenumber ) ), &
        TwoPi  =>  2.0_KDR  *  CONSTANT % PI )

    Rho  =  Offset  &
            +  Amplitude  &
               *  sin ( TwoPi * (    K ( 1 ) * X  &
                                  +  K ( 2 ) * Y  &
                                  +  K ( 3 ) * Z  ) )

    V ( 1 )  =  Speed  *  K ( 1 )  /  Abs_K
    V ( 2 )  =  Speed  *  K ( 2 )  /  Abs_K
    V ( 3 )  =  Speed  *  K ( 3 )  /  Abs_K
    
    end associate !-- Rho, etc.
    end select !-- C

  end subroutine SetWave


  subroutine TestRiemannSolver ( RSC, iD )

    class ( RiemannSolver_HLL_C_Form ), intent ( inout ) :: &
      RSC
    integer ( KDI ), intent ( in ) :: &
      iD

    call RSC % Compute ( iD )

    call GIS % Open ( GIS % ACCESS_CREATE )
    call SC % Write ( )
    call GIS % Close ( )

  end subroutine TestRiemannSolver


end program RiemannSolver_HLL_C__Form_Test
