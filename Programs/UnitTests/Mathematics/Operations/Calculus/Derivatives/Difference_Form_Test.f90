module Difference_Form_Test__Form

  use Basics
  use Manifolds
  use Difference_Form
  use Difference_GPU__Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'Difference_Form_Test'
  
  type, public :: Difference_Form_Test_Form
    type ( GridImageStreamForm ), allocatable :: &
      GridImageStream
    type ( Atlas_SC_Form ), allocatable :: &
      Atlas
    type ( DifferenceForm ), allocatable :: &
      Difference
    type ( Difference_GPU_Form ), allocatable :: &
      Difference_GPU
    type ( TimerForm ) :: &
      T_CPU, &
      T_GPU
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
      nCompute, &
      nVariables
    character ( 1 ) :: &
      iD_String
    type ( StorageForm ) :: &
      Gaussian
    type ( StorageForm ), dimension ( ATLAS % MAX_DIMENSIONS ) :: &
      d_Gaussian_Inner
    class ( GeometryFlatForm ), pointer :: &
      G
    
    !-- Stream

    allocate ( DFT % GridImageStream )
    associate ( GIS => DFT % GridImageStream )
    call GIS % Initialize &
           ( Name, CommunicatorOption = PROGRAM_HEADER % Communicator )
    
    !-- Atlas, etc.

    allocate ( DFT % Atlas )
    associate ( A => DFT % Atlas )
    call A % Initialize ( Name, PROGRAM_HEADER % Communicator )
    call A % CreateChart ( )
    call A % SetGeometry ( )
    
    call A % OpenStream ( GIS, '1', iStream = 1 )

    G => A % Geometry ( )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
       
    !-- Values to be differenced
    
    nVariables = 1 
    call PROGRAM_HEADER % GetParameter ( nVariables, 'nVariables' )

    call Gaussian % Initialize &
           ( [ G % nValues, nVariables ], NameOption = 'TestField', &
             VariableOption &
               = spread ( 'Gaussian                       ', &
                          dim = 1, ncopies = nVariables ) )
             

    associate &
      ( R0 => ( C % MaxCoordinate + C % MinCoordinate ) / 2.0_KDR, &
        L  => ( C % MaxCoordinate - C % MinCoordinate ), &
        X  => G % Value ( :, G % CENTER_U ( 1 ) ), &
        Y  => G % Value ( :, G % CENTER_U ( 2 ) ), &
        Z  => G % Value ( :, G % CENTER_U ( 3 ) ), &
        GaussianValue => Gaussian % Value ( :, 1 ) )
    associate &
      ( SigmaSq  => dot_product ( L, L ) / 64.0_KDR, &
        RadiusSq => ( X - R0 ( 1 ) ) ** 2  +  ( Y - R0 ( 2 ) ) ** 2 &
                    +  ( Z - R0 ( 3 ) ) ** 2 )
 
    GaussianValue = exp ( - RadiusSq / ( 2.0_KDR * SigmaSq ) )

    end associate !-- Sigma, etc.
    end associate !-- R0, etc.
    
    !-- Timers
    
    associate &
      ( T_CPU => DFT % T_CPU, &
        T_GPU => DFT % T_GPU )
    
    call T_CPU % Initialize ( 'Difference_CPU', Level = 1 )
    call T_GPU % Initialize ( 'Difference_GPU', Level = 1 )

    !-- Difference

    allocate ( DFT % Difference )
    allocate ( DFT % Difference_GPU )
    associate &
      ( D_CPU => DFT % Difference, &
        D_GPU => DFT % Difference_GPU )
    
    call D_CPU % Initialize ( 'GaussianDifference', shape ( Gaussian % Value ) )
    call D_GPU % Initialize ( 'GaussianDifference', shape ( Gaussian % Value ) )
    

    nCompute = 100
    call PROGRAM_HEADER % GetParameter ( nCompute, 'nCompute' )

    call Show ( 'Iterating D_CPU % Compute' )
    call T_CPU % Start ( )
    
    do iC = 1, nCompute
      if ( mod ( iC, 10 ) == 0 ) call Show ( iC, 'iC' )
      do iD = 1, C % nDimensions
        call D_CPU % Compute ( C, Gaussian, iD )
      end do !-- iD
    end do !-- iC
    
    call T_CPU % Stop ( )
    call T_CPU % ShowTotal ( CONSOLE % INFO_1 )
    
    call Show ( 'Iterating D_GPU % Compute' )
    
    call T_GPU % Start ( )
    
    do iC = 1, nCompute
      if ( mod ( iC, 10 ) == 0 ) call Show ( iC, 'iC' )
      call D_GPU % Compute ( C, Gaussian )
    end do !-- iC
    
    call T_GPU % Stop ( )
    call T_GPU % ShowTotal ( CONSOLE % INFO_1 )
    
    call Show ( T_CPU % TotalTime / T_GPU % TotalTime, 'GPU SpeedUp Factor' )
    call Show ( T_GPU % TotalTime / T_CPU % TotalTime, 'CPU SpeedUp Factor' )
    
    
    associate ( OI => D_CPU % OutputInner )

    !-- An extra iteration to output to disk
    do iD = 1, C % nDimensions

      call D_CPU % Compute ( C, Gaussian, iD )

      write ( iD_String, fmt = ' ( i1 ) ' ) iD

      associate ( dGI => d_Gaussian_Inner ( iD ) )
      call dGI % Initialize &
             ( [ OI % nValues, OI % nVariables ], &
               VariableOption = Gaussian % Variable, &
               NameOption = 'd_' // trim ( Gaussian % Name ) // '_' &
                            // iD_String )
      call Copy ( OI % Value, dGI % Value )
      end associate !-- dGI

    end do !-- iD

    end associate !-- OI

    end associate !-- Timers
    
    end associate !-- D

    !-- Write

    call C % AddFieldImage ( Gaussian, iStream = 1 )
    do iD = 1, C % nDimensions
      call C % AddFieldImage ( d_Gaussian_Inner ( iD ), iStream = 1 )
    end do

    call GIS % Open ( GIS % ACCESS_CREATE )
    call A % Write ( iStream = 1 )
    call GIS % Close ( )

    !-- Cleanup

    end select !-- C
    end associate !-- A
    end associate !-- GIS
    nullify ( G )
    
  end subroutine Initialize


  subroutine Finalize ( DFT )

    type ( Difference_Form_Test_Form ), intent ( inout ) :: &
      DFT

    deallocate ( DFT % Difference )
    deallocate ( DFT % Atlas )
    deallocate ( DFT % GridImageStream )

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
