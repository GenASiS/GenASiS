program Reconstruction_C__Form_Test

  !-- Reconstruction_Chart__Form_Test

  use Basics
  use Manifolds
  use Fields
  use Slopes

  implicit none

  integer ( KDI ) :: &
    iD
  character ( 1 ), dimension ( 3 ) :: &
    D  =  [ 'X', 'Y', 'Z' ]
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Chart_GS_Form ), allocatable :: &
    C
  type ( FieldSet_C_Form ), allocatable :: &
    FSC
  type ( FieldSet_C_Form ), dimension ( : ), allocatable :: &
    FSC_IL, FSC_IR, &
     DC_IL,  DC_IR
  type ( Stream_C_Form ), allocatable :: &
    SC
  type ( Geometry_F_C_Form ), allocatable :: &
    GC
  type ( Reconstruction_C_Form ), allocatable :: &
    RC_0, RC_1, RC_2

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Reconstruction_C__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( C )
  call C % Initialize &
         ( PeriodicOption = [ .true., .true., .true. ] )

  allocate ( FSC )
  call FSC % Initialize ( C )

  associate ( nD  =>  C % nDimensions )

  allocate ( FSC_IL ( nD ), FSC_IR ( nD ) )
  allocate (  DC_IL ( nD ),  DC_IR ( nD ) )
  do iD = 1, nD
    call FSC_IL ( iD ) % Initialize &
           ( C, NameOption = 'Fields_IL_' // D ( iD ) )
    call FSC_IR ( iD ) % Initialize &
           ( C, NameOption = 'Fields_IR_' // D ( iD ) )
    call DC_IL ( iD ) % Initialize &
          ( C, NameOption = 'Difference_IL_' // D ( iD ) )
    call DC_IR ( iD ) % Initialize &
          ( C, NameOption = 'Difference_IR_' // D ( iD ) )
  end do !-- iD

  allocate ( GC )
  call GC % Initialize ( C )

  allocate ( RC_0 )
  allocate ( RC_1 )
  allocate ( RC_2 )
  call RC_0 % Initialize &
         ( GC, FSC, &
           NameOption = 'Reconstruction_0', &
           StreamedOption = .true., &
           OrderOption = 0 )
  call RC_1 % Initialize &
         ( GC, FSC, &
           NameOption = 'Reconstruction_1', &
           StreamedOption = .true., &
           OrderOption = 1 )
  call RC_2 % Initialize &
         ( GC, FSC, &
           NameOption = 'Reconstruction_2', &
           StreamedOption = .true., &
           OrderOption = 2 )

  allocate ( SC )
  call SC % Initialize ( C, GIS )
  call SC % AddFieldSet ( FSC )
  do iD = 1, nD
    call SC % AddFieldSet ( FSC_IL ( iD ) )
    call SC % AddFieldSet ( FSC_IR ( iD ) )
    call SC % AddFieldSet (  DC_IL ( iD ) )
    call SC % AddFieldSet (  DC_IR ( iD ) )
  end do !-- iD

  call   C   % Show ( )
  call FSC   % Show ( )
  call  GC   % Show ( )
  call  RC_0 % Show ( )
  call  RC_1 % Show ( )
  call  RC_2 % Show ( )
  call  SC   % Show ( )

  call SetWave ( FSC, GC )
  do iD = 1, nD
    call SetReference ( FSC_IL ( iD ), FSC_IR ( iD ), GC, iD )
  end do !-- iD

  call TestReconstruction ( RC_0, SC, DC_IL, DC_IR, FSC_IL, FSC_IR )

  end associate !-- nD

  deallocate ( RC_2 )
  deallocate ( RC_1 )
  deallocate ( RC_0 )
  deallocate ( GC )
  deallocate ( DC_IR, DC_IL )
  deallocate ( FSC_IR, FSC_IL )
  deallocate ( FSC )
  deallocate ( C )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine SetWave ( FSC, GC, X_Option, Y_Option, Z_Option )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC
    class ( Geometry_F_C_Form ), intent ( in ), target :: &
      GC
    real ( KDR ), dimension ( : ), intent ( in ), target, optional :: &
      X_Option, Y_Option, Z_Option

    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Offset, &
      Amplitude
    real ( KDR ), dimension ( 3 ) :: &
      Wavenumber
    real ( KDR ), dimension ( : ), pointer :: &
      X, Y, Z

    select type ( C  =>  FSC % Chart )
    class is ( Chart_GS_Form )

    nWavelengths  =  0
    nWavelengths ( 1 : C % nDimensions )  =  1
    call PROGRAM_HEADER % GetParameter ( nWavelengths, 'nWavelengths' )

    Offset     =  2.0_KDR
    Amplitude  =  1.0_KDR
    call PROGRAM_HEADER % GetParameter ( Offset, 'Offset' )
    call PROGRAM_HEADER % GetParameter ( Amplitude, 'Amplitude' )

    associate ( BoxSize => C % MaxCoordinate - C % MinCoordinate )
    where ( BoxSize > 0.0_KDR )
      Wavenumber = nWavelengths / BoxSize
    elsewhere
      Wavenumber = 0.0_KDR
    end where
    end associate !-- BoxSize

    if ( present ( X_Option ) ) then
      X  =>  X_Option
    else
      X  =>  GC % Storage_FSC % Storage % Value ( :, GC % CENTER_U_1 )
    end if

    if ( present ( Y_Option ) ) then
      Y  =>  Y_Option
    else
      Y  =>  GC % Storage_FSC % Storage % Value ( :, GC % CENTER_U_2 )
    end if

    if ( present ( Z_Option ) ) then
      Z  =>  Z_Option
    else
      Z  =>  GC % Storage_FSC % Storage % Value ( :, GC % CENTER_U_3 )
    end if

    associate &
      (     N  =>  FSC % Storage_FSC % Storage % Value ( :, 1 ), &
            K  =>  Wavenumber, &
        TwoPi  =>  2.0_KDR  *  CONSTANT % PI )

    N  =  Offset  &
          +  Amplitude  &
             *  sin ( TwoPi * ( K ( 1 ) * X  +  K ( 2 ) * Y  +  K ( 3 ) * Z  ) )

    end associate !-- N, etc.
    end select !-- C

  end subroutine SetWave


  subroutine SetReference ( FSC_IL, FSC_IR, GC, iD )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC_IL, FSC_IR
    class ( Geometry_F_C_Form ), intent ( in ) :: &
      GC
    integer ( KDI ), intent ( in ) :: &
      iD

    associate &
      ( C  =>  GC % Storage_FSC % Storage % Value ( :, GC % CENTER_U ( iD ) ), &
        W  =>  GC % Storage_FSC % Storage % Value ( :, GC % WIDTH_U ( iD ) ) )

    select case ( iD )
    case ( 1 )
      call SetWave ( FSC_IL, GC, X_Option  =  C  -  0.5 * W )
      call SetWave ( FSC_IR, GC, X_Option  =  C  +  0.5 * W )
    case ( 2 )
      call SetWave ( FSC_IL, GC, Y_Option  =  C  -  0.5 * W )
      call SetWave ( FSC_IR, GC, Y_Option  =  C  +  0.5 * W )
    case ( 3 )
      call SetWave ( FSC_IL, GC, Z_Option  =  C  -  0.5 * W )
      call SetWave ( FSC_IR, GC, Z_Option  =  C  +  0.5 * W )
    end select !-- iD

    end associate !-- C, etc.

  end subroutine SetReference


  subroutine TestReconstruction ( RC, SC, DC_IL, DC_IR, FSC_IL, FSC_IR )

    class ( Reconstruction_C_Form ), intent ( inout ) :: &
      RC
    type ( Stream_C_Form ), intent ( inout ) :: &
      SC
    type ( FieldSet_C_Form ), dimension ( : ), intent ( inout ) :: &
      DC_IL, DC_IR
    type ( FieldSet_C_Form ), dimension ( : ), intent ( in ) :: &
      FSC_IL, FSC_IR

    associate ( nD  =>  RC % FieldSet_C % Chart % nDimensions )

    do iD  =  1, nD

      call RC % Compute ( iD )

      associate & 
        ( OV_IL  =>  RC % Output_IL_C % Storage_FSC % Storage % Value, &
          OV_IR  =>  RC % Output_IR_C % Storage_FSC % Storage % Value, &
          FV_IL  =>  FSC_IL ( iD ) % Storage_FSC % Storage % Value, &
          FV_IR  =>  FSC_IR ( iD ) % Storage_FSC % Storage % Value, &
          DV_IL  =>   DC_IL ( iD ) % Storage_FSC % Storage % Value, &
          DV_IR  =>   DC_IR ( iD ) % Storage_FSC % Storage % Value )

      DV_IL  =  OV_IL  -  FV_IL
      DV_IR  =  OV_IR  -  FV_IR

      end associate !-- OV_IL, etc.

    end do !-- iD

    call GIS % Open ( GIS % ACCESS_CREATE )
    call SC % Write ( )
    call GIS % Close ( )

    end associate !-- nD

  end subroutine TestReconstruction


end program Reconstruction_C__Form_Test
