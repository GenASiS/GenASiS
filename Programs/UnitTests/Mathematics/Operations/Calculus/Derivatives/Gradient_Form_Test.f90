module Gradient_Form_Test__Form

  use Basics
  use Manifolds
  use Gradient_Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'Gradient_Form_Test'
  
  type, public :: Gradient_Form_Test_Form
    type ( GridImageStreamForm ), allocatable :: &
      GridImageStream
    type ( Atlas_SC_Form ), allocatable :: &
      Atlas
    type ( GradientForm ), allocatable :: &
      Gradient
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type Gradient_Form_Test_Form

    integer, private, parameter :: &
      GAUSSIAN = 1, &
      TOP_HAT  = 2, &
      N_FIELDS = 2

contains


  subroutine Initialize ( GFT, Name )

    class ( Gradient_Form_Test_Form ), intent ( inout ), target :: &
      GFT
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      iC, &  !-- iCompute
      nCompute
    character ( 1 ) :: &
      iD_String
    type ( VariableGroupForm ) :: &
      TestFields
    type ( VariableGroupForm ), dimension ( ATLAS % MAX_DIMENSIONS ) :: &
      Grad_TestFields, &
      GradLimited_TestFields
    class ( GeometryFlatForm ), pointer :: &
      G

    !-- Stream

    allocate ( GFT % GridImageStream )
    associate ( GIS => GFT % GridImageStream )
    call GIS % Initialize &
           ( Name, CommunicatorOption = PROGRAM_HEADER % Communicator )
    
    !-- Atlas, etc.

    allocate ( GFT % Atlas )
    associate ( A => GFT % Atlas )
    call A % Initialize ( Name, PROGRAM_HEADER % Communicator )
    call A % CreateChart ( )
    call A % SetGeometry ( )
    
    call A % OpenStream ( GIS, '1', iStream = 1 )

    G => A % Geometry ( )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
       
    !-- Test fields

    call TestFields % Initialize &
           ( [ G % nValues, N_FIELDS ], NameOption = 'TestFields', &
             VariableOption = [ 'Gaussian                       ', &
                                'TopHat                         ' ] )

    associate &
      ( R0 => ( C % MaxCoordinate + C % MinCoordinate ) / 2.0_KDR, &
        L  => ( C % MaxCoordinate - C % MinCoordinate ), &
        X  => G % Value ( :, G % CENTER_U ( 1 ) ), &
        Y  => G % Value ( :, G % CENTER_U ( 2 ) ), &
        Z  => G % Value ( :, G % CENTER_U ( 3 ) ), &
        GaussianValue => TestFields % Value ( :, GAUSSIAN ), &
        TopHatValue   => TestFields % Value ( :, TOP_HAT ) )
    associate &
      ( SigmaSq  => dot_product ( L, L ) / 64.0_KDR, &
        RadiusSq => ( X - R0 ( 1 ) ) ** 2  +  ( Y - R0 ( 2 ) ) ** 2 &
                    +  ( Z - R0 ( 3 ) ) ** 2 )
    associate &
       ( TopHatRadius => sqrt ( SigmaSq ), &
         Radius => sqrt ( RadiusSq ) )
 
    GaussianValue = exp ( - RadiusSq / ( 2.0_KDR * SigmaSq ) )

    where ( Radius <= TopHatRadius )
      TopHatValue = 1.0_KDR
    elsewhere
      TopHatValue = 0.0_KDR
    end where

    end associate !-- TopHatRadius, etc.
    end associate !-- Sigma, etc.
    end associate !-- R0, etc.

    !-- Gradient

    allocate ( GFT % Gradient )
    associate ( Grad => GFT % Gradient )
    call Grad % Initialize &
           ( 'TestFieldsGradient', shape ( TestFields % Value ) )
    associate ( GO => Grad % Output )

    nCompute = 100
    call PROGRAM_HEADER % GetParameter ( nCompute, 'nCompute' )

    call Show ( 'Iterating Grad % Compute' )
    
    do iC = 1, nCompute
      call Show ( iC, 'iC' )
      do iD = 1, C % nDimensions
        call Grad % Compute ( C, TestFields, iD )
        call Grad % Compute ( C, TestFields, iD, &
                              LimiterParameterOption = 1.4_KDR )
      end do !-- iD
    end do !-- iC

    !-- An extra iteration to output to disk
    do iD = 1, C % nDimensions

      write ( iD_String, fmt = ' ( i1 ) ' ) iD

      !-- Without limiter

      call Grad % Compute ( C, TestFields, iD )

      associate ( Grad_TF => Grad_TestFields ( iD ) )
      call Grad_TF % Initialize &
             ( [ GO % nValues, GO % nVariables ], &
               VariableOption = TestFields % Variable, &
               NameOption = 'Grad_' // trim ( TestFields % Name ) // '_' &
                            // iD_String )
      call Copy ( GO % Value, Grad_TF % Value )
      end associate !-- Grad_TF

      !-- With limiter

      call Grad % Compute ( C, TestFields, iD, &
                            LimiterParameterOption = 1.4_KDR )

      associate ( GradLimited_TF => GradLimited_TestFields ( iD ) )
      call GradLimited_TF % Initialize &
             ( [ GO % nValues, GO % nVariables ], &
               VariableOption = TestFields % Variable, &
               NameOption = 'Grad_' // trim ( TestFields % Name ) // &
                            '_Limited_' // iD_String )
      call Copy ( GO % Value, GradLimited_TF % Value )
      end associate !-- Grad_TF

    end do !-- iD

    end associate !-- GO
    end associate !-- Grad

    !-- Write

    call C % AddFieldImage ( TestFields, iStream = 1 )
    do iD = 1, C % nDimensions
      call C % AddFieldImage ( Grad_TestFields ( iD ), iStream = 1 )
      call C % AddFieldImage ( GradLimited_TestFields ( iD ), iStream = 1 )
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


  subroutine Finalize ( GFT )

    type ( Gradient_Form_Test_Form ), intent ( inout ) :: &
      GFT

    deallocate ( GFT % Gradient )
    deallocate ( GFT % Atlas )
    deallocate ( GFT % GridImageStream )

  end subroutine Finalize


end module Gradient_Form_Test__Form



program Gradient_Form_Test

  use Basics
  use Gradient_Form_Test__Form

  implicit none
  
  type ( Gradient_Form_Test_Form ), allocatable :: &
    GFT
    
  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )
    
  allocate ( GFT )
  call GFT % Initialize ( PROGRAM_HEADER % Name )
  deallocate ( GFT )

  deallocate ( PROGRAM_HEADER )

end program Gradient_Form_Test
