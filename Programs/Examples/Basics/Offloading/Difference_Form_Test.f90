module Difference_Form_Test__Form

  use Basics
  use Difference_Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'Difference_Form_Test'
  
  type, public :: Difference_Form_Test_Form
    !type ( GridImageStreamForm ), allocatable :: &
    !  GridImageStream
    type ( DifferenceForm ), allocatable :: &
      Difference
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type Difference_Form_Test_Form

contains


  subroutine Initialize ( DFT, Name )

    class ( Difference_Form_Test_Form ), intent ( inout ), target :: &
      DFT
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      iC, &  !-- iCompute
      nDimensions, &
      nCompute
    integer ( KDI ), dimension ( 3 ) :: &
      nValues
    character ( 1 ) :: &
      iD_String
    type ( VariableGroupForm ) :: &
      Gaussian
    type ( VariableGroupForm ), dimension ( 1 ) :: &
      d_Gaussian_Inner
    
!;    !-- Stream
!;
!;    allocate ( DFT % GridImageStream )
!;    associate ( GIS => DFT % GridImageStream )
!;    call GIS % Initialize &
!;           ( Name, CommunicatorOption = PROGRAM_HEADER % Communicator )
    
    !-- Values to be differenced

    nValues = [ 32, 32, 32 ]
    call PROGRAM_HEADER % GetParameter ( nValues, 'nValues' )
    
    call Gaussian % Initialize &
           ( [ product( nValues ), 1 ], NameOption = 'TestField', &
             VariableOption = [ 'Gaussian                       ' ] )

!;    associate &
!;      ( R0 => ( C % MaxCoordinate + C % MinCoordinate ) / 2.0_KDR, &
!;        L  => ( C % MaxCoordinate - C % MinCoordinate ), &
!;        X  => G % Value ( :, G % CENTER_U ( 1 ) ), &
!;        Y  => G % Value ( :, G % CENTER_U ( 2 ) ), &
!;        Z  => G % Value ( :, G % CENTER_U ( 3 ) ), &
!;        GaussianValue => Gaussian % Value ( :, 1 ) )
!;    associate &
!;      ( SigmaSq  => dot_product ( L, L ) / 64.0_KDR, &
!;        RadiusSq => ( X - R0 ( 1 ) ) ** 2  +  ( Y - R0 ( 2 ) ) ** 2 &
!;                    +  ( Z - R0 ( 3 ) ) ** 2 )
!; 
!;    GaussianValue = exp ( - RadiusSq / ( 2.0_KDR * SigmaSq ) )
    
      call random_number ( Gaussian % Value  ( :, 1 ) )
!;    end associate !-- Sigma, etc.
!;    end associate !-- R0, etc.

    !-- Difference

    allocate ( DFT % Difference )
    associate ( D => DFT % Difference )
    call D % Initialize ( 'GaussianDifference', shape ( Gaussian % Value ) )
    associate ( OI => D % OutputInner )

    nCompute = 100
    call PROGRAM_HEADER % GetParameter ( nCompute, 'nCompute' )

    call Show ( 'Iterating D % Compute' )
    
    nDimensions = count ( nValues > 1 )
    do iC = 1, nCompute
      if ( mod ( iC, 10 ) == 0 ) call Show ( iC, 'iC' )
      do iD = 1, nDimensions
        call D % Compute ( Gaussian, nValues, iD )
      end do !-- iD
    end do !-- iC

!;    !-- An extra iteration to output to disk
!;    do iD = 1, C % nDimensions
!;
!;      call D % Compute ( C, Gaussian, iD )
!;
!;      write ( iD_String, fmt = ' ( i1 ) ' ) iD
!;
!;      associate ( dGI => d_Gaussian_Inner ( iD ) )
!;      call dGI % Initialize &
!;             ( [ OI % nValues, OI % nVariables ], &
!;               VariableOption = Gaussian % Variable, &
!;               NameOption = 'd_' // trim ( Gaussian % Name ) // '_' &
!;                            // iD_String )
!;      call Copy ( OI % Value, dGI % Value )
!;      end associate !-- dGI
!;
!;    end do !-- iD

    end associate !-- OI
    end associate !-- D

!;    !-- Write
!;
!;    call C % AddFieldImage ( Gaussian, iStream = 1 )
!;    do iD = 1, C % nDimensions
!;      call C % AddFieldImage ( d_Gaussian_Inner ( iD ), iStream = 1 )
!;    end do
!;
!;    call GIS % Open ( GIS % ACCESS_CREATE )
!;    call A % Write ( iStream = 1 )
!;    call GIS % Close ( )

    !-- Cleanup

!;    end associate !-- GIS
    
  end subroutine Initialize


  subroutine Finalize ( DFT )

    type ( Difference_Form_Test_Form ), intent ( inout ) :: &
      DFT

    deallocate ( DFT % Difference )
!;    deallocate ( DFT % GridImageStream )

  end subroutine Finalize


end module Difference_Form_Test__Form



program Difference_Form_Test

  use Basics
  use Difference_Form_Test__Form

  implicit none
  
  type ( Difference_Form_Test_Form ), allocatable :: &
    DFT
    
  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )
    
  allocate ( DFT )
  call DFT % Initialize ( PROGRAM_HEADER % Name )
  deallocate ( DFT )

  deallocate ( PROGRAM_HEADER )

end program Difference_Form_Test
