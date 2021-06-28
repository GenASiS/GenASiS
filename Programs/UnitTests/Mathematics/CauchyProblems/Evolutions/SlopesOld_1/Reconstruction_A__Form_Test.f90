program Reconstruction_A__Form_Test

  !-- Reconstruction_Atlas__Form_Test

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
  type ( Atlas_SCG_Form ), allocatable :: &
    A
  type ( FieldSet_A_Form ), allocatable :: &
    FSA
  type ( FieldSet_A_Form ), dimension ( : ), allocatable :: &
    FS_I_A, &
     D_IL_A,  D_IR_A
  type ( Stream_A_Form ), allocatable :: &
    SA
  type ( Geometry_F_A_Form ), allocatable :: &
    GA
  type ( Reconstruction_A_Form ), allocatable :: &
    RA_0, RA_1, RA_2

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Reconstruction_A__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           PeriodicOption = [ .true., .true., .true. ] )

  allocate ( FSA )
  call FSA % Initialize ( A )

  associate ( nD  =>  A % Chart ( 1 ) % Element % nDimensions )

  allocate ( FS_I_A ( nD ) )
  allocate ( D_IL_A ( nD ), D_IR_A ( nD ) )
  do iD = 1, nD
    call FS_I_A ( iD ) % Initialize &
           ( A, NameOption = 'Fields_I_' // D ( iD ) )
    call D_IL_A ( iD ) % Initialize &
          ( A, NameOption = 'Difference_IL_' // D ( iD ) )
    call D_IR_A ( iD ) % Initialize &
          ( A, NameOption = 'Difference_IR_' // D ( iD ) )
  end do !-- iD

  allocate ( SA )
  call SA % Initialize ( A, GIS )
  call SA % AddFieldSet ( FSA )
  do iD = 1, nD
    call SA % AddFieldSet ( FS_I_A  ( iD ) )
    call SA % AddFieldSet (  D_IL_A ( iD ) )
    call SA % AddFieldSet (  D_IR_A ( iD ) )
  end do !-- iD

  allocate ( GA )
  call GA % Initialize ( A )

  allocate ( RA_0 )
  allocate ( RA_1 )
  allocate ( RA_2 )
  call RA_0 % Initialize &
         ( GA, FSA, &
           NameOption = 'Reconstruction_0', &
           OrderOption = 0 )
  call RA_1 % Initialize &
         ( GA, FSA, &
           NameOption = 'Reconstruction_1', &
           OrderOption = 1 )
  call RA_2 % Initialize &
         ( GA, FSA, &
           NameOption = 'Reconstruction_2', &
           OrderOption = 2 )

  call   A   % Show ( )
  call FSA   % Show ( )
  call  GA   % Show ( )
  call  SA   % Show ( )
  call  RA_0 % Show ( )
  call  RA_1 % Show ( )
  call  RA_2 % Show ( )

  call SetWave ( FSA, GA )
  do iD = 1, nD
    call SetReference ( FS_I_A ( iD ), GA, iD )
  end do !-- iD

  call TestReconstruction ( RA_0, SA, D_IL_A, D_IR_A, FS_I_A )
  call TestReconstruction ( RA_1, SA, D_IL_A, D_IR_A, FS_I_A )
  call TestReconstruction ( RA_2, SA, D_IL_A, D_IR_A, FS_I_A )

  end associate !-- nD

  deallocate ( RA_2 )
  deallocate ( RA_1 )
  deallocate ( RA_0 )
  deallocate ( GA )
  deallocate ( SA )
  deallocate ( D_IR_A, D_IL_A )
  deallocate ( FS_I_A )
  deallocate ( FSA )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine SetWave ( FSA, GA, X_Option, Y_Option, Z_Option )

    class ( FieldSet_A_Form ), intent ( inout ) :: &
      FSA
    class ( Geometry_F_A_Form ), intent ( in ), target :: &
      GA
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

    associate ( FSC  =>  FSA % FieldSet_C ( 1 ) % Element )

    select type ( GC  =>  GA % FieldSet_C ( 1 ) % Element )
    class is ( Geometry_F_C_Form )

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
    end select !-- GC
    end associate !-- FSC

  end subroutine SetWave


  subroutine SetReference ( FS_I_A, GA, iD )

    class ( FieldSet_A_Form ), intent ( inout ) :: &
      FS_I_A
    class ( Geometry_F_A_Form ), intent ( in ) :: &
      GA
    integer ( KDI ), intent ( in ) :: &
      iD

    select type ( GC  =>  GA % FieldSet_C ( 1 ) % Element )
    class is ( Geometry_F_C_Form )

    associate &
      ( C  =>  GC % Storage_FSC % Storage % Value ( :, GC % CENTER_U ( iD ) ), &
        W  =>  GC % Storage_FSC % Storage % Value ( :, GC % WIDTH_U ( iD ) ) )

    select case ( iD )
    case ( 1 )
      call SetWave ( FS_I_A, GA, X_Option  =  C  -  0.5 * W )
    case ( 2 )
      call SetWave ( FS_I_A, GA, Y_Option  =  C  -  0.5 * W )
    case ( 3 )
      call SetWave ( FS_I_A, GA, Z_Option  =  C  -  0.5 * W )
    end select !-- iD

    end associate !-- C, etc.
    end select !-- GC

  end subroutine SetReference


  subroutine TestReconstruction ( RA, SA, D_IL_A, D_IR_A, FS_I_A )

    class ( Reconstruction_A_Form ), intent ( inout ) :: &
      RA
    type ( Stream_A_Form ), intent ( inout ) :: &
      SA
    type ( FieldSet_A_Form ), dimension ( : ), intent ( inout ) :: &
      D_IL_A, D_IR_A
    type ( FieldSet_A_Form ), dimension ( : ), intent ( in ) :: &
      FS_I_A

    associate ( RC  =>  RA % Reconstruction_C ( 1 ) % Element )
        
    associate ( nD  =>  RC % FieldSet_C % Chart % nDimensions )

    do iD  =  1, nD

      associate &
        ( D_IL_C  =>  D_IL_A ( iD ) % FieldSet_C ( 1 ) % Element, &
          D_IR_C  =>  D_IR_A ( iD ) % FieldSet_C ( 1 ) % Element, &
          FS_I_C  =>  FS_I_A ( iD ) % FieldSet_C ( 1 ) % Element )

      call RC % Compute ( iD )
      call CompareFieldSets ( RC % Output_IL_C, FS_I_C, iD )
      call CompareFieldSets ( RC % Output_IR_C, FS_I_C, iD )

      associate & 
        ( FV_I   =>  FS_I_C % Storage_FSC % Storage % Value, &
          OV_IL  =>  RC % Output_IL_C % Storage_FSC % Storage % Value, &
          OV_IR  =>  RC % Output_IR_C % Storage_FSC % Storage % Value, &
          DV_IL  =>  D_IL_C % Storage_FSC % Storage % Value, &
          DV_IR  =>  D_IR_C % Storage_FSC % Storage % Value )

      DV_IL  =  OV_IL  -  FV_I
      DV_IR  =  OV_IR  -  FV_I

      end associate !-- OV_IL, etc.
      end associate !-- D_IL_C, etc.

    end do !-- iD

    call GIS % Open ( GIS % ACCESS_CREATE )
    call SA % Write ( )
    call GIS % Close ( )

    end associate !-- nD
    end associate !-- RC, etc.

  end subroutine TestReconstruction


  subroutine CompareFieldSets ( FSC, FSC_R, iD )

    class ( FieldSet_C_Form ), intent ( in ) :: &
      FSC, &
      FSC_R
    integer ( KDI ), intent ( in ) :: &
      iD

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF     !-- iField
    type ( CollectiveOperation_R_Form ) :: &
      CO 

    call Show ( 'Comparing FieldSet with reference' )
    call Show ( FSC % Name, 'FieldSet' )
    call Show ( iD, 'iDimension' )

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

  end subroutine CompareFieldSets


end program Reconstruction_A__Form_Test
