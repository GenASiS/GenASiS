module Laplacian_M_H__Form

  !-- Laplacian_Multipole_Header_Form

  use Basics
  use Calculus
  use Manifolds
  use Fields

  implicit none
  private
  
    type, private :: Parameters_P_Form
      integer ( KDI ) :: &
        iDegree, &  !-- iL
        iOrder      !-- iM
      class ( Laplacian_M_H_Form ), pointer :: &
        Laplacian => null ( )
    end type Parameters_P_Form

  type, public :: Laplacian_M_H_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nRadialCells = 0, &
      nAngularMoments = 0, &
      nEquations = 0, &
      MaxDegree = 0, &  !-- Max L
      MaxOrder  = 0     !-- Max M
    integer ( KDI ) :: &
      iTimer = 0, &
      iTimer_CM = 0, &  !-- ClearMoments
      iTimer_LM = 0, &  !-- LocalMoments
      iTimer_RM = 0, &  !-- ReduceMoments
      iTimer_AM = 0     !-- AddMoments
    logical ( KDL ) :: &
      DeviceMemory = .false., &
      PinnedMemory = .false., &
      DevicesCommunicate = .false.
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      MyAngularMoment_3D => null ( ), &
        AngularMoment_3D => null ( ), &
       RadialMoment_R_3D => null ( ), &
       RadialMoment_I_3D => null ( )
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    character ( LDL ), dimension ( : ), allocatable :: &
      AngularFunctionName, &
      MomentName
    type ( StorageForm ), allocatable :: &
      d_Radius_3_3, &
      RadialFunctions_R, &
      RadialFunctions_I, &
      MyAngularMoments, &
        AngularMoments, &
      RadialMoments_R, &
      RadialMoments_I, &
      DeltaFactor
    type ( CollectiveOperation_R_Form ), allocatable :: &
      CO_AngularMoments
    type ( IntegralForm ), allocatable :: &
      Integral_P
    type ( Parameters_P_Form ), allocatable :: &
      Parameters_P
  contains
    procedure, public, pass :: &
      Initialize_H
    procedure, public, pass :: &
      Show => Show_L
    procedure, public, pass :: &
      Timer
    procedure, public, pass :: &
      ComputeMoments
    final :: &
      Finalize
    procedure, private, pass :: &
      SetParameters
    procedure, private, pass :: &
      SetParameters_A
    procedure, private, pass :: &
      SetKernelFunctions
    procedure, private, pass :: &
      AllocateMoments
    procedure, private, pass :: &
      ComputeAngularMomentsLocal
    procedure, private, pass :: &
      ComputeRadialMoments
    procedure, public, pass :: &
      ShowMoments
    procedure, public, nopass :: &
      AssociatedLegendre
  end type Laplacian_M_H_Form

    private :: &
      AllocateReduction, &
      AssignMomentPointers

    private :: &
      Integrand_P   

    private :: &
      ComputeRadialMomentsKernel

    interface

      module subroutine ComputeRadialMomentsKernel &
                          ( RM_R, RM_I, AM, RF_R, RF_I, dR33, nE, nAM, nR, &
                            UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
          RM_R, RM_I  !-- RadialMoment_Regular, _Irregular
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          AM  !-- AngularMoment
        real ( KDR ), dimension ( :, : ), intent ( in ) :: &
          RF_R, RF_I  !-- RadialFunction_Regular, _Irregular
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          dR33  !-- dRadius_3_3
        integer ( KDI ), intent ( in ) :: &
          nE, &   !-- nEquations
          nAM, &  !-- nAngularMoments
          nR      !-- nRadial
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine

    end interface


contains


  subroutine Initialize_H ( L, G, MaxDegree, nEquations )

    class ( Laplacian_M_H_Form ), intent ( inout ), target :: &
      L
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      MaxDegree, &
      nEquations

    L % IGNORABILITY  =  G % IGNORABILITY

    if ( L % Type  ==  '' ) &
      L % Type  =  'a Laplacian_M' 

    L % Name  =  'Laplacian'

    call Show ( 'Initializing ' // trim ( L % Type ), L % IGNORABILITY )
    call Show ( L % Name, 'Name', L % IGNORABILITY )

    allocate ( L % Parameters_P )
    allocate ( L % Integral_P )
    associate &
      ( IP  =>  L % Integral_P, &
        PP  =>  L % Parameters_P )

      call IP % Initialize ( PP )

      PP % Laplacian  =>  L
      IP % Integrand  =>  Integrand_P

    end associate !-- IP, etc.

    call L % SetParameters ( G, MaxDegree, nEquations )
    call L % SetKernelFunctions ( )
    call L % AllocateMoments ( )

  end subroutine Initialize_H


  subroutine Show_L ( L )

    class ( Laplacian_M_H_Form ), intent ( in ) :: &
      L

   character ( LDL ), dimension ( : ), allocatable :: &
     TypeWord

    call Split ( L % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', L % IGNORABILITY )
    call Show ( L % Name, 'Name', L % IGNORABILITY )

    call Show ( L % MaxDegree, 'MaxDegree (l)', L % IGNORABILITY )
    call Show ( L % MaxOrder, 'MaxOrder (m)', L % IGNORABILITY )
    call Show ( L % nRadialCells, 'nRadialCells', L % IGNORABILITY )
    call Show ( L % nAngularMoments, 'nAngularMoments', L % IGNORABILITY )
    call Show ( L % nEquations, 'nEquations', L % IGNORABILITY )
    call Show ( L % DeviceMemory, 'DeviceMemory', L % IGNORABILITY )
    call Show ( L % PinnedMemory, 'PinnedMemory', L % IGNORABILITY )
    call Show ( L % DevicesCommunicate, 'DevicesCommunicate', L % IGNORABILITY )

  end subroutine Show_L


  function Timer ( L, Level ) result ( T )

    class ( Laplacian_M_H_Form ), intent ( inout ) :: &
      L
    integer ( KDI ), intent ( in ) :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = L % iTimer, &
               Name = trim ( L % Name ) // '_CmptMmnts', &
               Level = Level )

  end function Timer


  subroutine ComputeMoments ( L, Source, T_Option )

    class ( Laplacian_M_H_Form ), intent ( inout ) :: &
      L
    class ( FieldSetForm ), intent ( inout ) :: &
      Source
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    type ( TimerForm ), pointer :: &
      T_CM, &
      T_LM, &
      T_RM, &
      T_AM

    if ( present ( T_Option ) ) then
      T_CM  =>  PROGRAM_HEADER % Timer &
                  ( Handle = L % iTimer_CM, &
                    Name = trim ( L % Name ) // '_ClrMmnts', &
                    Level = T_Option % Level + 1 )
      T_LM  =>  PROGRAM_HEADER % Timer &
                  ( Handle = L % iTimer_LM, &
                    Name = trim ( L % Name ) // '_LclMmnts', &
                    Level = T_Option % Level + 1 )
      T_RM  =>  PROGRAM_HEADER % Timer &
                  ( Handle = L % iTimer_RM, &
                    Name = trim ( L % Name ) // '_RdcMmnts', &
                    Level = T_Option % Level + 1 )
      T_AM  =>  PROGRAM_HEADER % Timer &
                  ( Handle = L % iTimer_AM, &
                    Name = trim ( L % Name ) // '_AddMmnts', &
                    Level = T_Option % Level + 1 )
    else
      T_CM  =>  null ( )
      T_LM  =>  null ( )
      T_RM  =>  null ( )
      T_AM  =>  null ( )
    end if

    call Show ( 'Computing Moments', L % IGNORABILITY + 4 )

    associate &
      (   AM  =>  L %   AngularMoments, &
        MyAM  =>  L % MyAngularMoments )

    if ( associated ( T_CM ) ) call T_CM % Start ( )
    call MyAM % Clear ( )
    if ( associated ( T_CM ) ) call T_CM % Stop ( )

    if ( associated ( T_LM ) ) call T_LM % Start ( )
    call L % ComputeAngularMomentsLocal ( Source )
    if ( associated ( T_LM ) ) call T_LM % Stop ( )

    if ( associated ( T_RM ) ) call T_RM % Start ( )
    if ( .not. L % DevicesCommunicate ) call MyAM % UpdateHost ( ) 
    call L % CO_AngularMoments % Reduce ( REDUCTION % SUM )
    if ( .not. L % DevicesCommunicate ) call AM % UpdateDevice ( ) 
    if ( associated ( T_RM ) ) call T_RM % Stop ( )

    if ( associated ( T_AM ) ) call T_AM % Start ( )
    call L % ComputeRadialMoments ( )
    if ( associated ( T_AM ) ) call T_AM % Stop ( )
    
    end associate !-- M, etc.

  end subroutine ComputeMoments


  impure elemental subroutine Finalize ( L )

    type ( Laplacian_M_H_Form ), intent ( inout ) :: &
      L

    if ( allocated ( L % CO_AngularMoments ) ) &
      deallocate ( L % CO_AngularMoments )

    if ( allocated ( L % DeltaFactor ) ) &
      deallocate ( L % DeltaFactor )
    if ( allocated ( L % RadialMoments_I ) ) &
      deallocate ( L % RadialMoments_I )
    if ( allocated ( L % RadialMoments_R ) ) &
      deallocate ( L % RadialMoments_R )
    if ( allocated ( L % AngularMoments ) ) &
      deallocate ( L % AngularMoments )
    if ( allocated ( L % MyAngularMoments ) ) &
      deallocate ( L % MyAngularMoments )
    if ( allocated ( L % RadialFunctions_I ) ) &
      deallocate ( L % RadialFunctions_I )
    if ( allocated ( L % RadialFunctions_R ) ) &
      deallocate ( L % RadialFunctions_R )
    if ( allocated ( L % d_Radius_3_3 ) ) &
      deallocate ( L % d_Radius_3_3 )

    nullify ( L % RadialMoment_I_3D )
    nullify ( L % RadialMoment_R_3D )
    nullify ( L % AngularMoment_3D )
    nullify ( L % MyAngularMoment_3D )

    if ( allocated ( L % MomentName ) ) &
      deallocate ( L % MomentName )
    if ( allocated ( L % AngularFunctionName ) ) &
      deallocate ( L % AngularFunctionName )

    if ( L % Name == '' ) &
      return

    call Show ( 'Finalizing ' // trim ( L % Type ), L % IGNORABILITY )
    call Show ( L % Name, 'Name', L % IGNORABILITY )
    
  end subroutine Finalize


  subroutine SetParameters ( L, G, MaxDegree, nEquations )

    class ( Laplacian_M_H_Form ), intent ( inout ) :: &
      L
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G
    integer ( KDI ), intent ( in ) :: &
      MaxDegree, &
      nEquations

    integer ( KDI ) :: &
      iL, &  !-- iDegree
      iM, &  !-- iOrder
      iE, &  !-- iEquation
      iA, &  !-- iAngular
      iAE    !-- iAngularEquations  
    character ( 1 ) :: &
      iE_Label
    character ( 2 ) :: &
      iL_Label, &
      iM_Label

    L % MaxDegree  =  MaxDegree
    L % MaxOrder   =  MaxDegree

    call L % SetParameters_A ( G )

    associate &
      (  L_Max  =>  L % MaxDegree, &
         M_Max  =>  L % MaxOrder, &
            nA  =>  L % nAngularMoments, &
            nE  =>  L % nEquations )

    nA = 0
    do iM  =  0, M_Max
      do iL  =  iM, L_Max
        nA  =  nA + 1
      end do
    end do
    nA  =  2 * nA  !-- Sine and Cosine parts 

    allocate ( L % AngularFunctionName ( nA ) )
    iA = 1
    do iM  =  0, M_Max
      do iL  =  iM, L_Max
        write ( iL_Label, fmt = '(i2.2)' ) iL
        write ( iM_Label, fmt = '(i2.2)' ) iM
        L % AngularFunctionName ( iA )  &
          =  'AngularFunction_' // iL_Label // '_' // iM_Label // '_Cos'
        L % AngularFunctionName ( iA + 1 )  &
          =  'AngularFunction_' // iL_Label // '_' // iM_Label // '_Sin'
        iA  =  iA + 2
      end do  !-- iL
    end do  !-- iM

    nE  =  nEquations

    allocate ( L % MomentName ( nA * nE ) )
    iAE = 1
    do iE  =  1, nE
      do iM  =  0, M_Max
        do iL  =  iM, L_Max
          write ( iL_Label, fmt = '(i2.2)' ) iL
          write ( iM_Label, fmt = '(i2.2)' ) iM
          write ( iE_Label, fmt = '(i1.1)' ) iE
          L % MomentName ( iAE )  &
            =  'Moment_' // iL_Label // '_' // iM_Label // '_Cos_Eq_' &
               // iE_Label
          L % MomentName ( iAE + 1 )  &
            =  'Moment_' // iL_Label // '_' // iM_Label // '_Sin_Eq_' &
               // iE_Label
          iAE  =  iAE + 2
        end do  !-- iL
      end do  !-- iM
    end do  !-- iE

    allocate ( L % DeltaFactor )
    associate ( DF => L % DeltaFactor )
    call DF % Initialize ( [ nA, 1 ] )
    iA = 1
    do iM  =  0, M_Max
      do iL  =  iM, L_Max
        if ( iM == 0 ) then
          DF % Value ( iA : iA + 1, 1 )  &
            =  1.0_KDR / ( 2.0_KDR * iL  +  1.0_KDR )
        else
          DF % Value ( iA : iA + 1, 1 )  &
            =  2.0_KDR / ( 2.0_KDR * iL  +  1.0_KDR )
        end if
        iA  =  iA + 2
      end do  !-- iL
    end do  !-- iM
    if ( L % DeviceMemory ) then
      call DF % AllocateDevice ( AssociateVariablesOption = .false. )
      call DF % UpdateDevice ( )
    end if
    end associate !-- DF

    end associate !-- L_Max, etc.

  end subroutine SetParameters


  subroutine SetParameters_A ( L, G )

    class ( Laplacian_M_H_Form ), intent ( inout ) :: &
      L
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G

    call Show ( 'Subroutine should be overidden', CONSOLE % ERROR )
    call Show ( 'Laplacian_M_H__Form', 'module', CONSOLE % ERROR )
    call Show ( 'SetParameters_A', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine SetParameters_A


  subroutine SetKernelFunctions ( L )

    class ( Laplacian_M_H_Form ), intent ( inout ) :: &
      L

    call Show ( 'Subroutine should be overidden', CONSOLE % ERROR )
    call Show ( 'Laplacian_M_H__Form', 'module', CONSOLE % ERROR )
    call Show ( 'SetKernelFunctions', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine SetKernelFunctions


  subroutine AllocateMoments ( L )

    class ( Laplacian_M_H_Form ), intent ( inout ) :: &
      L

!     if ( allocated ( L % AngularMoments ) ) &
!       deallocate ( L % AngularMoments )
!     if ( allocated ( L % MyAngularMoments ) ) &
!       deallocate ( L % MyAngularMoments )

    allocate ( L % AngularMoments )
    allocate ( L % MyAngularMoments )
    allocate ( L % RadialMoments_R )
    allocate ( L % RadialMoments_I )
    associate &
      (         AM  =>  L % AngularMoments, &
              MyAM  =>  L % MyAngularMoments, &
              RM_R  =>  L % RadialMoments_R, &
              RM_I  =>  L % RadialMoments_I, &
               nAM  =>  L % nAngularMoments, &
                nR  =>  L % nRadialCells, &
                nE  =>  L % nEquations )

    call   AM % Initialize &
             ( [ nR, nAM * nE ], &
               NameOption = 'AngularMoments', &
               VariableOption = L % MomentName, &
               PinnedOption = L % PinnedMemory )
    call MyAM % Initialize &
             ( [ nR, nAM * nE ], &
               NameOption = 'MyAngularMoments', &
               VariableOption = L % MomentName, &
               PinnedOption = L % PinnedMemory )
    call RM_R % Initialize &
             ( [ nR + 1, nAM * nE ], &
               NameOption = 'RadialMoments_R', &
               VariableOption = L % MomentName, &
               PinnedOption = L % PinnedMemory )
    call RM_I % Initialize &
             ( [ nR + 1, nAM * nE ], &
               NameOption = 'RadialMoments_I', &
               VariableOption = L % MomentName, &
               PinnedOption = L % PinnedMemory )
    if ( L % DeviceMemory ) then
      call   AM % AllocateDevice ( AssociateVariablesOption = .false. )
      call MyAM % AllocateDevice ( AssociateVariablesOption = .false. )
      call RM_R % AllocateDevice ( AssociateVariablesOption = .false. )
      call RM_I % AllocateDevice ( AssociateVariablesOption = .false. )
    end if

    call Clear ( MyAM % Value, UseDeviceOption = L % DeviceMemory )

    call AllocateReduction ( L, AM % Value, MyAM % Value )

    call AssignMomentPointers &
           ( L, L %  AngularMoments % Value, L % MyAngularMoments % Value, &
                L % RadialMoments_R % Value, L %  RadialMoments_I % Value, &
                L % AngularMoment_3D,        L % MyAngularMoment_3D, &
                L % RadialMoment_R_3D,         L % RadialMoment_I_3D )

    end associate !-- M, etc.

  end subroutine AllocateMoments


  subroutine ComputeRadialMoments ( L )

    class ( Laplacian_M_H_Form ), intent ( inout ) :: &
      L

    call ComputeRadialMomentsKernel &
           ( L % RadialMoment_R_3D, L % RadialMoment_I_3D, &
             L % AngularMoment_3D, L % RadialFunctions_R % Value, &
             L % RadialFunctions_I % Value, L % d_Radius_3_3 % Value ( :, 1 ), &
             L % nEquations, L % nAngularMoments, L % nRadialCells, &
             UseDeviceOption = L % DeviceMemory )
    
  end subroutine ComputeRadialMoments


  subroutine ShowMoments ( L )

    class ( Laplacian_M_H_Form ), intent ( inout ) :: &
      L

    integer ( KDI ) :: &
      iL, &  !-- iDegree
      iM, &  !-- iOrder
      iE, &  !-- iEquation
      iA, &  !-- iAngular
      iAE    !-- iAngularEquations  

    associate &
      (  L_Max  =>  L % MaxDegree, &
         M_Max  =>  L % MaxOrder, &
            nE  =>  L % nEquations )
    
    if ( L % DeviceMemory ) then
      call L %  AngularMoments % UpdateHost ( )
      call L % RadialMoments_R % UpdateHost ( )
      call L % RadialMoments_I % UpdateHost ( )
    end if

    call Show ( 'Displaying Angular Moments' )

    iAE = 1
    do iE  =  1, nE
      iA = 1
      do iM  =  0, M_Max
        do iL  =  iM, L_Max
          !-- Cos and Sin
          call Show ( L % AngularMoment_3D ( :, iA, iE ), &
                      L % MomentName ( iAE ) )
          call Show ( L % AngularMoment_3D ( :, iA + 1, iE ), &
                      L % MomentName ( iAE + 1) )
           iA  =   iA + 2
          iAE  =  iAE + 2
        end do  !-- iL
      end do  !-- iM
    end do  !-- iE

    call Show ( 'Displaying Regular Radial Moments' )

    iAE = 1
    do iE  =  1, nE
      iA = 1
      do iM  =  0, M_Max
        do iL  =  iM, L_Max
          !-- Cos and Sin
          call Show ( L % RadialMoment_R_3D ( :, iA, iE ), &
                      L % MomentName ( iAE ) )
          call Show ( L % RadialMoment_R_3D ( :, iA + 1, iE ), &
                      L % MomentName ( iAE + 1) )
           iA  =   iA + 2
          iAE  =  iAE + 2
        end do  !-- iL
      end do  !-- iM
    end do  !-- iE

    call Show ( 'Displaying Irregular Moments' )

    iAE = 1
    do iE  =  1, nE
      iA = 1
      do iM  =  0, M_Max
        do iL  =  iM, L_Max
          !-- Cos and Sin
          call Show ( L % RadialMoment_I_3D ( :, iA, iE ), &
                      L % MomentName ( iAE ) )
          call Show ( L % RadialMoment_I_3D ( :, iA + 1, iE ), &
                      L % MomentName ( iAE + 1) )
           iA  =   iA + 2
          iAE  =  iAE + 2
        end do  !-- iL
      end do  !-- iM
    end do  !-- iE

    end associate !-- L_Max, etc.

  end subroutine ShowMoments


  function AssociatedLegendre ( X, L, M ) result ( P_LM )
  
    !-- Normalized, see Numerical Recipes Third Edition, Section 6.7

    real ( KDR ), intent ( in ) :: &
      X
    integer ( KDI ), intent ( in ) :: &
      L, &
      M
    real ( KDR ) :: &
      P_LM

    integer ( KDI ) :: &
      iM, &
      iL
    real ( KDR ) :: &
      P_MM, &
      P_MP1_M, &  !-- P_M+1_M
      P_LL, &
      O_M_X2, &  !-- 1 - x ** 2
      Factor, &
      FactorOld, &
      FourPi

    if ( M < 0 .or. M > L .or. abs ( X ) > 1.0_KDR ) then
      call Show ( 'Arguments out of range', CONSOLE % ERROR )
      call Show ( M, 'M', CONSOLE % ERROR )
      call Show ( L, 'L', CONSOLE % ERROR )
      call Show ( X, 'X', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_Template', 'module', CONSOLE % ERROR )
      call Show ( 'AssociatedLegendre', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    FourPi  =  4.0_KDR  *  CONSTANT % PI

    !-- Compute P_mm

    P_MM  =  1.0_KDR

    if ( M > 0 ) then
      O_M_X2  =  ( 1.0_KDR - X ) * ( 1.0_KDR + X )
      Factor  =  1.0_KDR
      do iM  =  1, M
          P_MM  =  P_MM  *  Factor / ( Factor + 1.0_KDR )  *  O_M_X2  
        Factor  =  Factor + 2.0_KDR
      end do !-- iM
    end if !-- iM > 0

    P_MM  =  sqrt ( ( 2 * M + 1 ) * P_MM / FourPi )
    if ( mod ( M, 2 ) == 1 )  &
      P_MM  =  - P_MM

    if ( L == M ) then

      P_LM  =  P_MM
      return

    else

      !-- Compute P_lm

      P_MP1_M  =  sqrt ( 2.0_KDR * M  +  3.0_KDR )  *  X  *  P_MM

      if ( L == M + 1 ) then
        P_LM  =  P_MP1_M
        return

      else 

        !-- Compute P_lm, l > m + 1
        FactorOld  =  sqrt ( 2.0_KDR * M  +  3.0_KDR )
        do iL  =  M + 2, L
              Factor  =  sqrt ( ( 4.0_KDR * iL * iL  - 1.0_KDR )  &
                                /  ( iL * iL  -  M * M ) )
                P_LL  =  Factor * ( X * P_MP1_M  -  P_MM / FactorOld )
           FactorOld  =  Factor
                P_MM  =  P_MP1_M
             P_MP1_M  =  P_LL
        end do !-- iL

        P_LM  =  P_LL
        return

      end if  !-- L == M + 1
    end if !-- L == M

  end function AssociatedLegendre


  subroutine ComputeAngularMomentsLocal ( L, Source )

    class ( Laplacian_M_H_Form ), intent ( inout ) :: &
      L
    class ( FieldSetForm ), intent ( inout ) :: &
      Source

    call Show ( 'Subroutine should be overidden', CONSOLE % ERROR )
    call Show ( 'Laplacian_M_H__Form', 'module', CONSOLE % ERROR )
    call Show ( 'ComputeAngularMomentsLocal', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine ComputeAngularMomentsLocal


  subroutine AllocateReduction ( L, AM_Value, MyAM_Value )

    class ( Laplacian_M_H_Form ), intent ( inout ) :: &
      L
    real ( KDR ), dimension ( :, : ), intent ( in ), target, contiguous :: &
        AM_Value, &
      MyAM_Value

    real ( KDR ), dimension ( : ), pointer :: &
        AngularMoment_1D, &
      MyAngularMoment_1D

      AngularMoment_1D ( 1 : size (   AM_Value ) )  =>    AM_Value
    MyAngularMoment_1D ( 1 : size ( MyAM_Value ) )  =>  MyAM_Value

    if ( allocated ( L % CO_AngularMoments ) ) &
      deallocate ( L % CO_AngularMoments )
    allocate ( L % CO_AngularMoments )
    associate &
      (  CO  =>  L % CO_AngularMoments, &
        PHC  =>  PROGRAM_HEADER % Communicator )
      call CO % Initialize &
             ( PHC, OutgoingValue = MyAngularMoment_1D, &
               IncomingValue = AngularMoment_1D )
      if ( L % DevicesCommunicate ) &
        call CO % AllocateDevice ( )
    end associate !-- RM, etc.

    nullify ( AngularMoment_1D, MyAngularMoment_1D )

  end subroutine AllocateReduction


  subroutine AssignMomentPointers &
               ( L, AM_2D, MyAM_2D, RM_R_2D, RM_I_2D, &
                    AM_3D, MyAM_3D, RM_R_3D, RM_I_3D )

    class ( Laplacian_M_H_Form ), intent ( in ) :: &
      L
    real ( KDR ), dimension ( :, : ), intent ( in ), target, contiguous :: &
        AM_2D, &
      MyAM_2D, &
      RM_R_2D, &
      RM_I_2D
    real ( KDR ), dimension ( :, :, : ), intent ( out ), pointer :: &
        AM_3D, &
      MyAM_3D, &
      RM_R_3D, &
      RM_I_3D

    associate &
      (  nE => L % nEquations, &
        nAM => L % nAngularMoments, &
         nR => L % nRadialCells )

      AM_3D ( 1 : nR,     1 : nAM, 1 : nE )  =>    AM_2D
    MyAM_3D ( 1 : nR,     1 : nAM, 1 : nE )  =>  MyAM_2D
    RM_R_3D ( 1 : nR + 1, 1 : nAM, 1 : nE )  =>  RM_R_2D
    RM_I_3D ( 1 : nR + 1, 1 : nAM, 1 : nE )  =>  RM_I_2D

    end associate  !-- nE, nA, nR

  end subroutine AssignMomentPointers


  function Integrand_P ( Parameters, X ) result ( IP )

    class ( * ), intent ( in ) :: &
      Parameters
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      IP

    select type ( P => Parameters )
      type is ( Parameters_P_Form )
    associate &
      ( L   =>  P % Laplacian, &
        iL  =>  P % iDegree, &
        iM  =>  P % iOrder )

    IP  =  L % AssociatedLegendre ( X, iL, iM )

    end associate !-- L, etc.
    end select !-- Parameters

  end function Integrand_P


end module Laplacian_M_H__Form
