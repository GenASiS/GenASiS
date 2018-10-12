!-- MeasuredValue_Form handles numbers with labels (typically units) providing
!   means of dealing with units, for example. 

module MeasuredValue_Form
  
  use KIND_DEFAULT_Singleton
  use KIND_BIG_Singleton
  use LEN_DEFAULT_Singleton
  use Split_Command
  use Join_Command
  
  implicit none
  private
  
  type, public :: MeasuredValueForm
    real ( KDR ) :: &
      Number = 1.0_KDR
    character ( LDL, KDCH ) :: &
      Unit = ''
    character ( LDL, KDCH ) :: &
      Label = '' !-- Actual unit name if being used as a unit
    character ( LDL, KBCH ) :: &
      Unit_UCS = ''
    character ( LDL, KBCH ) :: &
      Label_UCS = '' !-- Actual unit name if being used as a unit
  contains
    procedure, private, pass :: &
      Initialize_MV
    procedure, private, pass :: &
      Initialize_MV_Label
    procedure, private, pass :: &
      Initialize_MV_From_MV
    procedure, private, pass :: &
      Initialize_MV_From_MV_Label
    generic :: &
      Initialize &
        => Initialize_MV, Initialize_MV_Label, &
           Initialize_MV_From_MV, Initialize_MV_From_MV_Label
    procedure, private, pass :: &
      Initialize_UCS_MV
    procedure, private, pass :: &
      Initialize_UCS_MV_Label
    procedure, private, pass :: &
      Initialize_UCS_MV_From_MV_Label
    generic :: &
      Initialize_UCS &
        => Initialize_UCS_MV, Initialize_UCS_MV_Label, &
           Initialize_UCS_MV_From_MV_Label
    procedure, private, pass :: &
      Addition_MV_MV
    procedure, private, pass :: &
      Addition_MV_Real
    procedure, private, pass ( Addend ) :: &
      AdditionReal_MV
    procedure, private, pass :: &
      Addition_MV_Integer
    procedure, private, pass ( Addend ) :: &
      AdditionInteger_MV
    generic :: &
      operator ( + ) &
        => Addition_MV_MV, Addition_MV_Real, AdditionReal_MV, &
           Addition_MV_Integer, AdditionInteger_MV
    procedure, private, pass :: &
      Subtraction_MV_MV
    procedure, private, pass :: &
      Subtraction_MV_Real
    procedure, private, pass ( Subtrahend ) :: &
      SubtractionReal_MV
    procedure, private, pass :: &
      Subtraction_MV_Integer
    procedure, private, pass ( Subtrahend ) :: &
      SubtractionInteger_MV
    generic :: &
      operator ( - ) &
        => Subtraction_MV_MV, Subtraction_MV_Real, SubtractionReal_MV, &
           Subtraction_MV_Integer, SubtractionInteger_MV
    procedure, private, pass :: &
      Product_MV_MV
    procedure, private, pass :: &
      Product_MV_Real
    procedure, private, pass ( Multiplicand ) :: &
      ProductReal_MV
    procedure, private, pass :: &
      Product_MV_Integer
    procedure, private, pass ( Multiplicand ) :: &
      ProductInteger_MV
    generic :: &
      operator ( * ) &
        => Product_MV_MV, Product_MV_Real, ProductReal_MV, &
           Product_MV_Integer, ProductInteger_MV
    procedure, private, pass :: &
      Power_MV_Integer
!    procedure, private, pass :: &
!      Power_MV_Real
    generic :: &
      operator ( ** ) => Power_MV_Integer!, Power_MV_Real
    procedure, private, pass :: &
      Quotient_MV_MV
    procedure, private, pass :: &
      Quotient_MV_Real
    procedure, private, pass ( Divisor ) :: &
      QuotientReal_MV
    procedure, private, pass :: &
      Quotient_MV_Integer
    procedure, private, pass ( Divisor ) :: &
      QuotientInteger_MV
    generic :: &
      operator ( / ) &
        => Quotient_MV_MV, Quotient_MV_Real, QuotientReal_MV, &
           Quotient_MV_Integer, QuotientInteger_MV
    procedure, private, pass ( MV ) :: &
      AssignReal_MV
    generic :: &
      assignment ( = ) => AssignReal_MV
    procedure, private, pass :: &
      EqualTo_MV_MV
    procedure, private, pass :: &
      EqualTo_MV_Real
    procedure, private, pass ( MV ) :: &
      EqualToReal_MV
    procedure, private, pass :: &
      EqualTo_MV_Integer
    procedure, private, pass ( MV ) :: &
      EqualToInteger_MV
    generic :: &
      operator ( == ) &
        => EqualTo_MV_MV, EqualTo_MV_Real, EqualToReal_MV, &
           EqualTo_MV_Integer, EqualToInteger_MV
    procedure, private, pass :: &
      NotEqualTo_MV_MV
    procedure, private, pass :: &
      NotEqualTo_MV_Real
    procedure, private, pass ( MV ) :: &
      NotEqualToReal_MV
    procedure, private, pass :: &
      NotEqualTo_MV_Integer
    procedure, private, pass ( MV ) :: &
      NotEqualToInteger_MV
    generic :: &
      operator ( /= ) &
        => NotEqualTo_MV_MV, NotEqualTo_MV_Real, NotEqualToReal_MV, &
           NotEqualTo_MV_Integer, NotEqualToInteger_MV
    procedure, private, pass :: &
      GreaterThan_MV_MV
    procedure, private, pass :: &
      GreaterThan_MV_Real
    procedure, private, pass ( MV ) :: &
      GreaterThanReal_MV
    procedure, private, pass :: &
      GreaterThan_MV_Integer
    procedure, private, pass ( MV ) :: &
      GreaterThanInteger_MV
    generic :: &
      operator ( > ) &
        => GreaterThan_MV_MV, GreaterThan_MV_Real, GreaterThanReal_MV, &
           GreaterThan_MV_Integer, GreaterThanInteger_MV
    procedure, private, pass :: &
      LessThan_MV_MV
    procedure, private, pass :: &
      LessThan_MV_Real
    procedure, private, pass ( MV ) :: &
      LessThanReal_MV
    procedure, private, pass :: &
      LessThan_MV_Integer
    procedure, private, pass ( MV ) :: &
      LessThanInteger_MV
    generic :: &
      operator ( < ) &
        => LessThan_MV_MV, LessThan_MV_Real, LessThanReal_MV, &
           LessThan_MV_Integer, LessThanInteger_MV
    procedure, private, pass :: &
      GreaterThanEqualTo_MV_MV
    procedure, private, pass :: &
      GreaterThanEqualTo_MV_Real
    procedure, private, pass ( MV ) :: &
      GreaterThanEqualToReal_MV
    procedure, private, pass :: &
      GreaterThanEqualTo_MV_Integer
    procedure, private, pass ( MV ) :: &
      GreaterThanEqualToInteger_MV
    generic :: &
      operator ( >= ) &
        => GreaterThanEqualTo_MV_MV, GreaterThanEqualTo_MV_Real, &
           GreaterThanEqualToReal_MV, GreaterThanEqualTo_MV_Integer, &
           GreaterThanEqualToInteger_MV
    procedure, private, pass :: &
      LessThanEqualTo_MV_MV
    procedure, private, pass :: &
      LessThanEqualTo_MV_Real
    procedure, private, pass ( MV ) :: &
      LessThanEqualToReal_MV
    procedure, private, pass :: &
      LessThanEqualTo_MV_Integer
    procedure, private, pass ( MV ) :: &
      LessThanEqualToInteger_MV
    generic :: &
      operator ( <= ) &
        => LessThanEqualTo_MV_MV, LessThanEqualTo_MV_Real, &
           LessThanEqualToReal_MV, LessThanEqualTo_MV_Integer, &
           LessThanEqualToInteger_MV
    procedure, private, nopass :: &
      ProductUnit
    procedure, private, nopass :: &
      ProductUnit_UCS
    procedure, private, nopass :: &
      PowerUnitInteger
    procedure, private, nopass :: &
      PowerUnitInteger_UCS
!    procedure, private, nopass :: &
!      PowerUnitReal
!     generic :: &
!       PowerUnit => PowerUnitInteger!, PowerUnitReal
  end type MeasuredValueForm

    private :: &
      AnalyzeUnit

    integer ( KDI ), private, parameter :: &
      ExponentMinus = KIND_BIG % SUPERSCRIPT_MINUS
    integer ( KDI ), dimension ( 9 ), private, parameter :: &
      Exponent_I &
        = [ KIND_BIG % SUPERSCRIPT_1, KIND_BIG % SUPERSCRIPT_2, &
            KIND_BIG % SUPERSCRIPT_3, KIND_BIG % SUPERSCRIPT_4, &
            KIND_BIG % SUPERSCRIPT_5, KIND_BIG % SUPERSCRIPT_6, &
            KIND_BIG % SUPERSCRIPT_7, KIND_BIG % SUPERSCRIPT_8, &
            KIND_BIG % SUPERSCRIPT_9 ]
    character ( 1, KBCH ), private, parameter :: &
      Exponent_1_Char   = char ( KIND_BIG % SUPERSCRIPT_1, KBCH ), &
      Exponent_2_Char   = char ( KIND_BIG % SUPERSCRIPT_2, KBCH ), &
      Exponent_3_Char   = char ( KIND_BIG % SUPERSCRIPT_3, KBCH ), &
      Exponent_4_Char   = char ( KIND_BIG % SUPERSCRIPT_4, KBCH ), &
      Exponent_5_Char   = char ( KIND_BIG % SUPERSCRIPT_5, KBCH ), &
      Exponent_6_Char   = char ( KIND_BIG % SUPERSCRIPT_6, KBCH ), &
      Exponent_7_Char   = char ( KIND_BIG % SUPERSCRIPT_7, KBCH ), &
      Exponent_8_Char   = char ( KIND_BIG % SUPERSCRIPT_8, KBCH ), &
      Exponent_9_Char   = char ( KIND_BIG % SUPERSCRIPT_9, KBCH ), &
      ExponentMinusChar = char ( KIND_BIG % SUPERSCRIPT_MINUS, KBCH )
    character ( 1, KBCH ), dimension ( 9 ), private, parameter :: &
      Exponent_I_Char &
        = [ Exponent_1_Char, Exponent_2_Char, Exponent_3_Char, &
            Exponent_4_Char, Exponent_5_Char, Exponent_6_Char, &
            Exponent_7_Char, Exponent_8_Char, Exponent_9_Char ]
    character ( 9, KBCH ), private, parameter :: &
      ExponentCharacters &
        =    Exponent_1_Char // Exponent_2_Char // Exponent_3_Char &
          // Exponent_4_Char // Exponent_5_Char // Exponent_6_Char &
          // Exponent_7_Char // Exponent_8_Char // Exponent_9_Char

contains


  subroutine Initialize_MV ( MV, Unit, Number )
    
    class ( MeasuredValueForm ), intent ( inout ) :: &
      MV
    character ( * ), intent ( in ) :: &
      Unit
    real ( KDR ), intent ( in ) :: &
      Number
      
    MV % Number   = Number
    MV % Unit     = Unit
    MV % Unit_UCS = Unit

  end subroutine Initialize_MV
  
  
  subroutine Initialize_MV_Label ( MV, Label, Unit, Number )
    
    class ( MeasuredValueForm ), intent ( inout ) :: &
      MV
    character ( * ), intent ( in ) :: &
      Label, &
      Unit
    real ( KDR ), intent ( in ) :: &
      Number
      
    MV % Number    = Number
    MV % Unit      = Unit
    MV % Label     = Label
    MV % Unit_UCS  = Unit
    MV % Label_UCS = Label

  end subroutine Initialize_MV_Label
  
  
  subroutine Initialize_MV_From_MV ( MV_Target, MV_Source )
    
    class ( MeasuredValueForm ), intent ( inout ) :: &
      MV_Target
    type ( MeasuredValueForm ), intent ( in ) :: &
      MV_Source
            
    MV_Target % Number    = MV_Source % Number
    MV_Target % Unit      = MV_Source % Unit
    MV_Target % Label     = MV_Source % Label
    MV_Target % Unit_UCS  = MV_Source % Unit_UCS
    MV_Target % Label_UCS = MV_Source % Label_UCS
  
  end subroutine Initialize_MV_From_MV
  
  
  subroutine Initialize_MV_From_MV_Label ( MV_Target, MV_Source, Label )
    
    class ( MeasuredValueForm ), intent ( inout ) :: &
      MV_Target
    type ( MeasuredValueForm ), intent ( in ) :: &
      MV_Source
    character ( * ), intent ( in ) :: &
      Label
            
    MV_Target % Number    = MV_Source % Number
    MV_Target % Unit      = MV_Source % Unit
    MV_Target % Label     = Label
    MV_Target % Unit_UCS  = MV_Source % Unit_UCS
    MV_Target % Label_UCS = Label
  
  end subroutine Initialize_MV_From_MV_Label
  
  
  subroutine Initialize_UCS_MV ( MV, Unit_UCS, Unit, Number )
    
    class ( MeasuredValueForm ), intent ( inout ) :: &
      MV
    character ( *, KBCH ), intent ( in ) :: &
      Unit_UCS
    character ( * ), intent ( in ) :: &
      Unit
    real ( KDR ), intent ( in ) :: &
      Number
      
    MV % Number   = Number
    MV % Unit     = Unit
    MV % Unit_UCS = Unit_UCS

  end subroutine Initialize_UCS_MV
  
  
  subroutine Initialize_UCS_MV_Label &
               ( MV, Label_UCS, Unit_UCS, Label, Unit, Number )
    
    class ( MeasuredValueForm ), intent ( inout ) :: &
      MV
    character ( *, KBCH ), intent ( in ) :: &
      Label_UCS
    character ( *, KBCH ), intent ( in ) :: &
      Unit_UCS
    character ( * ), intent ( in ) :: &
      Label
    character ( * ), intent ( in ) :: &
      Unit
    real ( KDR ), intent ( in ) :: &
      Number
      
    MV % Number    = Number
    MV % Unit      = Unit
    MV % Label     = Label
    MV % Unit_UCS  = Unit_UCS
    MV % Label_UCS = Label_UCS

  end subroutine Initialize_UCS_MV_Label
  
  
  subroutine Initialize_UCS_MV_From_MV_Label &
               ( MV_Target, MV_Source, Label_UCS, Label )
    
    class ( MeasuredValueForm ), intent ( inout ) :: &
      MV_Target
    type ( MeasuredValueForm ), intent ( in ) :: &
      MV_Source
    character ( *, KBCH ), intent ( in ) :: &
      Label_UCS
    character ( * ), intent ( in ) :: &
      Label
            
    MV_Target % Number    = MV_Source % Number
    MV_Target % Unit      = MV_Source % Unit
    MV_Target % Label     = Label
    MV_Target % Unit_UCS  = MV_Source % Unit_UCS
    MV_Target % Label_UCS = Label_UCS
  
  end subroutine Initialize_UCS_MV_From_MV_Label
  
  
  elemental function Addition_MV_MV ( Augend, Addend ) result ( A_MV_MV )

    class ( MeasuredValueForm ), intent ( in ) :: &
      Augend
    type ( MeasuredValueForm ), intent ( in ) :: &
      Addend
    type ( MeasuredValueForm ) :: &
      A_MV_MV
    
      A_MV_MV % Number = Augend % Number + Addend % Number
    
      if ( trim ( Augend % Unit ) == trim ( Addend % Unit ) ) then
        A_MV_MV % Unit     = Augend % Unit
        A_MV_MV % Unit_UCS = Augend % Unit_UCS
      else
        A_MV_MV % Unit     = 'Undefined'
        A_MV_MV % Unit_UCS = 'Undefined'
      end if
  
      if ( trim ( Augend % Label ) == trim ( Addend % Label ) ) then
        A_MV_MV % Label     = Augend % Label
        A_MV_MV % Label_UCS = Augend % Label_UCS
      else
        A_MV_MV % Label     = 'Undefined'
        A_MV_MV % Label_UCS = 'Undefined'
      end if
  
  end function Addition_MV_MV
  
  
  elemental function Addition_MV_Real &
                       ( Augend, Addend ) result ( A_MV_R )
  
    class ( MeasuredValueForm ), intent ( in ) :: &
      Augend
    real ( KDR ), intent ( in ) :: &
      Addend
    real ( KDR ) :: &
      A_MV_R
      
    A_MV_R = Augend % Number + Addend
      
  end function Addition_MV_Real
   
  
  elemental function AdditionReal_MV &
                       ( Augend, Addend ) result ( A_R_MV )
  
    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "Augend + Addend" 

    real ( KDR ), intent ( in ) :: &
      Augend
    class ( MeasuredValueForm ), intent ( in ) :: &
      Addend
    real ( KDR ) :: &
      A_R_MV
      
    A_R_MV = Augend + Addend % Number
      
  end function AdditionReal_MV
   
  
  elemental function Addition_MV_Integer &
                       ( Augend, Addend ) result ( A_MV_I )
  
    class ( MeasuredValueForm ), intent ( in ) :: &
      Augend
    integer ( KDI ), intent ( in ) :: &
      Addend
    real ( KDR ) :: &
      A_MV_I
      
    A_MV_I = Augend % Number + real ( Addend, KDR )
      
  end function Addition_MV_Integer
   
  
  elemental function AdditionInteger_MV &
                       ( Augend, Addend ) result ( A_I_MV )
  
    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "Augend + Addend" 

    integer ( KDI ), intent ( in ) :: &
      Augend
    class ( MeasuredValueForm ), intent ( in ) :: &
      Addend
    real ( KDR ) :: &
      A_I_MV
      
    A_I_MV = real ( Augend, KDR ) + Addend % Number
      
  end function AdditionInteger_MV
   
  
  elemental function Subtraction_MV_MV &
                       ( Minuend, Subtrahend ) result ( S_MV_MV )
  
    class ( MeasuredValueForm ), intent ( in ) :: &
      Minuend
    class ( MeasuredValueForm ), intent ( in ) :: &
      Subtrahend
    type ( MeasuredValueForm ) :: &
      S_MV_MV
      
    S_MV_MV % Number = Minuend % Number - Subtrahend % Number
      
    if ( trim ( Minuend % Unit ) == trim ( Subtrahend % Unit ) ) then
      S_MV_MV % Unit     = Minuend % Unit
      S_MV_MV % Unit_UCS = Minuend % Unit_UCS
    else
      S_MV_MV % Unit     = 'Undefined'
      S_MV_MV % Unit_UCS = 'Undefined'
    end if
    
    if ( trim ( Minuend % Label ) == trim ( Subtrahend % Label ) ) then
      S_MV_MV % Label     = Minuend % Label
      S_MV_MV % Label_UCS = Minuend % Label_UCS
    else
      S_MV_MV % Label     = 'Undefined'
      S_MV_MV % Label_UCS = 'Undefined'
    end if
    
  end function Subtraction_MV_MV
  
  
  elemental function Subtraction_MV_Real &
                       ( Minuend, Subtrahend ) result ( S_MV_R )
  
    class ( MeasuredValueForm ), intent ( in ) :: &
      Minuend
    real ( KDR ), intent ( in ) :: &
      Subtrahend
    real ( KDR ) :: &
      S_MV_R
      
    S_MV_R = Minuend % Number - Subtrahend
      
  end function Subtraction_MV_Real
   
  
  elemental function SubtractionReal_MV &
                       ( Minuend, Subtrahend ) result ( S_R_MV )
  
    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "Minuend - Subtrahend" 

    real ( KDR ), intent ( in ) :: &
      Minuend
    class ( MeasuredValueForm ), intent ( in ) :: &
      Subtrahend
    real ( KDR ) :: &
      S_R_MV
      
    S_R_MV = Minuend - Subtrahend % Number
      
  end function SubtractionReal_MV
   
  
  elemental function Subtraction_MV_Integer &
                       ( Minuend, Subtrahend ) result ( S_MV_I )
  
    class ( MeasuredValueForm ), intent ( in ) :: &
      Minuend
    integer ( KDI ), intent ( in ) :: &
      Subtrahend
    real ( KDR ) :: &
      S_MV_I
      
    S_MV_I = Minuend % Number - real ( Subtrahend, KDR )
      
  end function Subtraction_MV_Integer
   
  
  elemental function SubtractionInteger_MV &
                       ( Minuend, Subtrahend ) result ( S_I_MV )
  
    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "Minuend - Subtrahend" 

    integer ( KDI ), intent ( in ) :: &
      Minuend
    class ( MeasuredValueForm ), intent ( in ) :: &
      Subtrahend
    real ( KDR ) :: &
      S_I_MV
      
    S_I_MV = real ( Minuend, KDR ) - Subtrahend % Number
      
  end function SubtractionInteger_MV
   
  
  impure elemental function Product_MV_MV &
                     ( Multiplier, Multiplicand ) result ( P_MV_MV )
  
    class ( MeasuredValueForm ), intent ( in ) :: &
      Multiplier, &
      Multiplicand
    type ( MeasuredValueForm ) :: &
      P_MV_MV

    P_MV_MV % Number = Multiplier % Number * Multiplicand % Number

    P_MV_MV % Unit &
      = P_MV_MV % ProductUnit ( Multiplier % Unit, Multiplicand % Unit )
    P_MV_MV % Label &
      = P_MV_MV % ProductUnit ( Multiplier % Label, Multiplicand % Label )

    if ( KBCH > KDCH ) then
      P_MV_MV % Unit_UCS &
        = P_MV_MV % ProductUnit_UCS &
            ( Multiplier % Unit_UCS, Multiplicand % Unit_UCS )
      P_MV_MV % Label_UCS &
        = P_MV_MV % ProductUnit_UCS &
            ( Multiplier % Label_UCS, Multiplicand % Label_UCS )
    else
      P_MV_MV % Unit_UCS &
        = P_MV_MV % ProductUnit &
            ( Multiplier % Unit, Multiplicand % Unit )
      P_MV_MV % Label_UCS &
        = P_MV_MV % ProductUnit &
            ( Multiplier % Label, Multiplicand % Label )
    end if

  end function Product_MV_MV
  
  
  elemental function Product_MV_Real &
                       ( Multiplier, Multiplicand ) result ( P_MV_R )
    
    class ( MeasuredValueForm ), intent ( in ) :: &
      Multiplier
    real ( KDR ), intent ( in ) :: &
      Multiplicand
    type ( MeasuredValueForm ) :: &
      P_MV_R
    
    P_MV_R % Number    = Multiplier % Number * Multiplicand
    P_MV_R % Unit      = Multiplier % Unit
    P_MV_R % Label     = Multiplier % Label
    P_MV_R % Unit_UCS  = Multiplier % Unit_UCS
    P_MV_R % Label_UCS = Multiplier % Label_UCS
  
  end function Product_MV_Real
  

  elemental function ProductReal_MV &
                        ( Multiplier, Multiplicand ) result ( P_R_MV )
    
    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "Multiplier * Multiplicand" 

    real ( KDR ), intent ( in ) :: &
      Multiplier
    class ( MeasuredValueForm ), intent ( in ) :: &
      Multiplicand
    type ( MeasuredValueForm ) :: &
      P_R_MV
    
    P_R_MV % Number    = Multiplier * Multiplicand % Number
    P_R_MV % Unit      = Multiplicand % Unit
    P_R_MV % Label     = Multiplicand % Label
    P_R_MV % Unit_UCS  = Multiplicand % Unit_UCS
    P_R_MV % Label_UCS = Multiplicand % Label_UCS
    
  end function ProductReal_MV
  
  
  elemental function Product_MV_Integer &
                       ( Multiplier, Multiplicand ) result ( P_MV_I )
    
    class ( MeasuredValueForm ), intent ( in ) :: &
      Multiplier
    integer ( KDI ), intent ( in ) :: &
      Multiplicand
    type ( MeasuredValueForm ) :: &
      P_MV_I
    
    P_MV_I % Number    = Multiplier % Number * real ( Multiplicand, KDR )
    P_MV_I % Unit      = Multiplier % Unit
    P_MV_I % Label     = Multiplier % Label
    P_MV_I % Unit_UCS  = Multiplier % Unit_UCS
    P_MV_I % Label_UCS = Multiplier % Label_UCS
  
  end function Product_MV_Integer
  

  elemental function ProductInteger_MV &
                        ( Multiplier, Multiplicand ) result ( P_I_MV )
    
    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "Multiplier * Multiplicand" 

    integer ( KDI ), intent ( in ) :: &
      Multiplier
    class ( MeasuredValueForm ), intent ( in ) :: &
      Multiplicand
    type ( MeasuredValueForm ) :: &
      P_I_MV
    
    P_I_MV % Number    = real ( Multiplier, KDR ) * Multiplicand % Number
    P_I_MV % Unit      = Multiplicand % Unit
    P_I_MV % Label     = Multiplicand % Label
    P_I_MV % Unit_UCS  = Multiplicand % Unit_UCS
    P_I_MV % Label_UCS = Multiplicand % Label_UCS
    
  end function ProductInteger_MV
  
  
  elemental function Power_MV_Integer ( Base, Exponent ) result ( P_MV_I )
    
    class ( MeasuredValueForm ), intent ( in ) :: &
      Base
    integer ( KDI ), intent ( in ) :: &
      Exponent
    type ( MeasuredValueForm ) :: &
      P_MV_I
    
    P_MV_I % Number = Base % Number ** Exponent

    P_MV_I % Unit  = P_MV_I % PowerUnitInteger ( Base % Unit, Exponent )
    P_MV_I % Label = P_MV_I % PowerUnitInteger ( Base % Label, Exponent )

    if ( KBCH > KDCH ) then
      P_MV_I % Unit_UCS &
        = P_MV_I % PowerUnitInteger_UCS ( Base % Unit_UCS, Exponent )
      P_MV_I % Label_UCS &
        = P_MV_I % PowerUnitInteger_UCS ( Base % Label_UCS, Exponent )
    else
      P_MV_I % Unit_UCS &
        = P_MV_I % PowerUnitInteger ( Base % Unit, Exponent )
      P_MV_I % Label_UCS &
        = P_MV_I % PowerUnitInteger ( Base % Label, Exponent )
    end if

  end function Power_MV_Integer
  

  ! elemental function Power_MV_Real ( Base, Exponent ) result ( P_MV_R )
    
  !   class ( MeasuredValueForm ), intent ( in ) :: &
  !     Base
  !   real ( KDR ), intent ( in ) :: &
  !     Exponent
  !   type ( MeasuredValueForm ) :: &
  !     P_MV_R
    
  !   P_MV_R % Number = Base % Number ** Exponent
    
  !   P_MV_R % Unit = P_MV_R % PowerUnit ( Base % Unit, Exponent )
    
  !   P_MV_R % Label = P_MV_R % PowerUnit ( Base % Label, Exponent )
    
  ! end function Power_MV_Real
  

  impure elemental function Quotient_MV_MV ( Dividend, Divisor ) &
                              result ( Q_MV_MV )
    
    class ( MeasuredValueForm ), intent ( in ) :: &
      Dividend, &
      Divisor
    type ( MeasuredValueForm ) :: &
      P_MV_I, &
      Q_MV_MV
    
    !-- FIXME: Workaround for GCC 5.1 bug 66257
    P_MV_I  = Divisor ** ( -1 )
    
    Q_MV_MV = Dividend * P_MV_I

  end function Quotient_MV_MV
  
  
  elemental function Quotient_MV_Real ( Dividend, Divisor ) result ( Q_MV_R )
    
    class ( MeasuredValueForm ), intent ( in ) :: &
      Dividend
    real ( KDR ), intent ( in ) :: &
      Divisor
    type ( MeasuredValueForm ) :: &
      Q_MV_R
    
    Q_MV_R = Dividend * Divisor ** ( -1 )

  end function Quotient_MV_Real
  
  
  elemental function QuotientReal_MV ( Dividend, Divisor ) result ( Q_R_MV )
    
    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "Dividend / Divisor" 

    real ( KDR ), intent ( in ) :: &
      Dividend
    class ( MeasuredValueForm ), intent ( in ) :: &
      Divisor
    type ( MeasuredValueForm ) :: &
      P_MV_I, &
      Q_R_MV
    
    !-- FIXME: Workaround for GCC 5.1 bug 66257
    P_MV_I = Divisor ** ( -1 )
    
    Q_R_MV = Dividend * P_MV_I

  end function QuotientReal_MV
  
  
  elemental function Quotient_MV_Integer ( Dividend, Divisor ) result ( Q_MV_I )
    
    class ( MeasuredValueForm ), intent ( in ) :: &
      Dividend
    integer ( KDI ), intent ( in ) :: &
      Divisor
    type ( MeasuredValueForm ) :: &
      Q_MV_I
    
    Q_MV_I = Dividend * real ( Divisor, KDR ) ** ( -1 )

  end function Quotient_MV_Integer
  
  
  elemental function QuotientInteger_MV ( Dividend, Divisor ) result ( Q_I_MV )
    
    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "Dividend / Divisor" 

    integer ( KDI ), intent ( in ) :: &
      Dividend
    class ( MeasuredValueForm ), intent ( in ) :: &
      Divisor
    type ( MeasuredValueForm ) :: &
      P_MV_I, &
      Q_I_MV
    
    !-- FIXME: Workaround for GCC 5.1 bug 66257
    P_MV_I = Divisor ** ( -1 )
    
    Q_I_MV = real ( Dividend, KDR ) * P_MV_I

  end function QuotientInteger_MV
  
  
  elemental subroutine AssignReal_MV ( A, MV )

    real ( KDR ), intent ( inout ) :: &
      A
    class ( MeasuredValueForm ), intent ( in ) :: &
      MV

    A = MV % Number

  end subroutine AssignReal_MV


  elemental function EqualTo_MV_MV ( MV_1, MV_2 ) result ( ET )

    class ( MeasuredValueForm ), intent ( in ) :: &
      MV_1, MV_2
    logical ( KDL ) :: &
      ET

    ET = .false.

    if ( MV_1 % Unit == MV_2 % Unit &
         .and. MV_1 % Number == MV_2 % Number ) &
      ET = .true.

  end function EqualTo_MV_MV


  elemental function EqualTo_MV_Real ( MV, A ) result ( ET )

    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    real ( KDR ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      ET

    ET = .false.

    if ( MV % Number == A ) ET = .true.

  end function EqualTo_MV_Real


  elemental function EqualToReal_MV ( A, MV ) result ( ET )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A == MV" 

    real ( KDR ), intent ( in ) :: &
      A
    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    logical ( KDL ) :: &
      ET

    ET = .false.

    if ( A == MV % Number ) ET = .true.

  end function EqualToReal_MV


  elemental function EqualTo_MV_Integer ( MV, A ) result ( ET )

    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    integer ( KDI ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      ET

    ET = .false.

    if ( MV % Number == real ( A, KDR ) ) ET = .true.

  end function EqualTo_MV_Integer


  elemental function EqualToInteger_MV ( A, MV ) result ( ET )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A == MV" 

    integer ( KDI ), intent ( in ) :: &
      A
    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    logical ( KDL ) :: &
      ET

    ET = .false.

    if ( real ( A, KDR ) == MV % Number ) ET = .true.

  end function EqualToInteger_MV


  elemental function NotEqualTo_MV_MV ( MV_1, MV_2 ) result ( NET )

    class ( MeasuredValueForm ), intent ( in ) :: &
      MV_1, MV_2
    logical ( KDL ) :: &
      NET

    NET = .false.

    if ( MV_1 % Unit /= MV_2 % Unit &
         .or. MV_1 % Number /= MV_2 % Number ) &
      NET = .true.

  end function NotEqualTo_MV_MV


  elemental function NotEqualTo_MV_Real ( MV, A ) result ( NET )

    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    real ( KDR ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      NET

    NET = .false.

    if ( MV % Number /= A ) NET = .true.

  end function NotEqualTo_MV_Real


  elemental function NotEqualToReal_MV ( A, MV ) result ( NET )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A /= MV" 

    real ( KDR ), intent ( in ) :: &
      A
    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    logical ( KDL ) :: &
      NET

    NET = .false.

    if ( A /= MV % Number ) NET = .true.

  end function NotEqualToReal_MV


  elemental function NotEqualTo_MV_Integer ( MV, A ) result ( NET )

    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    integer ( KDI ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      NET

    NET = .false.

    if ( MV % Number /= real ( A, KDR ) ) NET = .true.

  end function NotEqualTo_MV_Integer


  elemental function NotEqualToInteger_MV ( A, MV ) result ( NET )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A /= MV" 

    integer ( KDI ), intent ( in ) :: &
      A
    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    logical ( KDL ) :: &
      NET

    NET = .false.

    if ( real ( A, KDR ) /= MV % Number ) NET = .true.

  end function NotEqualToInteger_MV


  elemental function GreaterThan_MV_MV ( MV_1, MV_2 ) result ( GT )

    class ( MeasuredValueForm ), intent ( in ) :: &
      MV_1, MV_2
    logical ( KDL ) :: &
      GT

    GT = .false.

    if ( MV_1 % Unit == MV_2 % Unit &
         .and. MV_1 % Number > MV_2 % Number ) &
      GT = .true.

  end function GreaterThan_MV_MV


  elemental function GreaterThan_MV_Real ( MV, A ) result ( GT )

    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    real ( KDR ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      GT

    GT = .false.

    if ( MV % Number > A ) GT = .true.

  end function GreaterThan_MV_Real


  elemental function GreaterThanReal_MV ( A, MV ) result ( GT )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A > MV" 

    real ( KDR ), intent ( in ) :: &
      A
    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    logical ( KDL ) :: &
      GT

    GT = .false.

    if ( A > MV % Number ) GT = .true.

  end function GreaterThanReal_MV


  elemental function GreaterThan_MV_Integer ( MV, A ) result ( GT )

    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    integer ( KDI ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      GT

    GT = .false.

    if ( MV % Number > real ( A, KDR ) ) GT = .true.

  end function GreaterThan_MV_Integer


  elemental function GreaterThanInteger_MV ( A, MV ) result ( GT )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A > MV" 

    integer ( KDI ), intent ( in ) :: &
      A
    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    logical ( KDL ) :: &
      GT

    GT = .false.

    if ( real ( A, KDR ) > MV % Number ) GT = .true.

  end function GreaterThanInteger_MV


  elemental function LessThan_MV_MV ( MV_1, MV_2 ) result ( LT )

    class ( MeasuredValueForm ), intent ( in ) :: &
      MV_1, MV_2
    logical ( KDL ) :: &
      LT

    LT = .false.

    if ( MV_1 % Unit == MV_2 % Unit &
         .and. MV_1 % Number < MV_2 % Number ) &
      LT = .true.

  end function LessThan_MV_MV


  elemental function LessThan_MV_Real ( MV, A ) result ( LT )

    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    real ( KDR ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      LT

    LT = .false.

    if ( MV % Number < A ) LT = .true.

  end function LessThan_MV_Real


  elemental function LessThanReal_MV ( A, MV ) result ( LT )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A < MV" 

    real ( KDR ), intent ( in ) :: &
      A
    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    logical ( KDL ) :: &
      LT

    LT = .false.

    if ( A < MV % Number ) LT = .true.

  end function LessThanReal_MV


  elemental function LessThan_MV_Integer ( MV, A ) result ( LT )

    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    integer ( KDI ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      LT

    LT = .false.

    if ( MV % Number < real ( A, KDR ) ) LT = .true.

  end function LessThan_MV_Integer


  elemental function LessThanInteger_MV ( A, MV ) result ( LT )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A < MV" 

    integer ( KDI ), intent ( in ) :: &
      A
    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    logical ( KDL ) :: &
      LT

    LT = .false.

    if ( real ( A, KDR ) < MV % Number ) LT = .true.

  end function LessThanInteger_MV


  elemental function GreaterThanEqualTo_MV_MV ( MV_1, MV_2 ) result ( GTE )

    class ( MeasuredValueForm ), intent ( in ) :: &
      MV_1, MV_2
    logical ( KDL ) :: &
      GTE

    GTE = .false.

    if ( MV_1 % Unit == MV_2 % Unit &
         .and. MV_1 % Number >= MV_2 % Number ) &
      GTE = .true.

  end function GreaterThanEqualTo_MV_MV


  elemental function GreaterThanEqualTo_MV_Real ( MV, A ) result ( GTE )

    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    real ( KDR ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      GTE

    GTE = .false.

    if ( MV % Number >= A ) GTE = .true.

  end function GreaterThanEqualTo_MV_Real


  elemental function GreaterThanEqualToReal_MV ( A, MV ) result ( GTE )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A >= MV" 

    real ( KDR ), intent ( in ) :: &
      A
    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    logical ( KDL ) :: &
      GTE

    GTE = .false.

    if ( A >= MV % Number ) GTE = .true.

  end function GreaterThanEqualToReal_MV


  elemental function GreaterThanEqualTo_MV_Integer ( MV, A ) result ( GTE )

    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    integer ( KDI ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      GTE

    GTE = .false.

    if ( MV % Number >= real ( A, KDR ) ) GTE = .true.

  end function GreaterThanEqualTo_MV_Integer


  elemental function GreaterThanEqualToInteger_MV ( A, MV ) result ( GTE )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A >= MV" 

    integer ( KDI ), intent ( in ) :: &
      A
    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    logical ( KDL ) :: &
      GTE

    GTE = .false.

    if ( real ( A, KDR ) >= MV % Number ) GTE = .true.

  end function GreaterThanEqualToInteger_MV


  elemental function LessThanEqualTo_MV_MV ( MV_1, MV_2 ) result ( LTE )

    class ( MeasuredValueForm ), intent ( in ) :: &
      MV_1, MV_2
    logical ( KDL ) :: &
      LTE

    LTE = .false.

    if ( MV_1 % Unit == MV_2 % Unit &
         .and. MV_1 % Number <= MV_2 % Number ) &
      LTE = .true.

  end function LessThanEqualTo_MV_MV


  elemental function LessThanEqualTo_MV_Real ( MV, A ) result ( LTE )

    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    real ( KDR ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      LTE

    LTE = .false.

    if ( MV % Number <= A ) LTE = .true.

  end function LessThanEqualTo_MV_Real


  elemental function LessThanEqualToReal_MV ( A, MV ) result ( LTE )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A <= MV" 

    real ( KDR ), intent ( in ) :: &
      A
    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    logical ( KDL ) :: &
      LTE

    LTE = .false.

    if ( A <= MV % Number ) LTE = .true.

  end function LessThanEqualToReal_MV


  elemental function LessThanEqualTo_MV_Integer ( MV, A ) result ( LTE )

    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    integer ( KDI ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      LTE

    LTE = .false.

    if ( MV % Number <= real ( A, KDR ) ) LTE = .true.

  end function LessThanEqualTo_MV_Integer


  elemental function LessThanEqualToInteger_MV ( A, MV ) result ( LTE )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A <= MV" 

    integer ( KDI ), intent ( in ) :: &
      A
    class ( MeasuredValueForm ), intent ( in ) :: &
      MV
    logical ( KDL ) :: &
      LTE

    LTE = .false.

    if ( real ( A, KDR ) <= MV % Number ) LTE = .true.

  end function LessThanEqualToInteger_MV


  recursive function ProductUnit ( Unit_1, Unit_2 ) result ( P_U )

    character ( * ), intent ( in ) :: &
      Unit_1, &
      Unit_2
    character ( LDL ) :: &
      P_U

    integer ( KDI ) :: &
      iU_1, iU_2, &
      Caret, &
      Dot, &
      Exponent_1, &
      Exponent_2, &
      Exponent
    character ( LDL ) :: &
      U_1, &
      U_2, &
      UnitBase_1, &
      UnitBase_2, &
      ExponentString
    character ( LDL ), dimension ( : ), allocatable :: &
      UnitPiece_1, &
      UnitPiece_2
  
    P_U = ''

    if ( trim ( Unit_1 ) == '' .and. trim ( Unit_2 ) == '' ) &
      return

    call Split ( Unit_1, ' ', UnitPiece_1 )
    call Split ( Unit_2, ' ', UnitPiece_2 )
    if ( size ( UnitPiece_1 ) > 1 .or. size ( UnitPiece_2 ) > 1 ) then
      do iU_1 = 1, size ( UnitPiece_1 )
        do iU_2 = 1, size ( UnitPiece_2 )
          Caret = index ( UnitPiece_1 ( iU_1 ), '^' )
          if ( Caret > 0 ) then
            UnitBase_1 = UnitPiece_1 ( iU_1 ) ( : Caret - 1 )
          else
            UnitBase_1 = UnitPiece_1 ( iU_1 )
          end if
          if ( index ( UnitPiece_2 ( iU_2 ), trim ( UnitBase_1 ) ) > 0 ) then
            UnitPiece_1 ( iU_1 ) &
              = ProductUnit ( UnitPiece_1 ( iU_1 ), UnitPiece_2 ( iU_2 ) )
            UnitPiece_2 ( iU_2 ) = ''
          end if
        end do !-- iU_2
      end do !-- iU_1
      call Join ( UnitPiece_1, ' ', U_1 )
      call Join ( UnitPiece_2, ' ', U_2 )
      P_U = adjustl ( trim ( U_1 ) // ' ' // trim ( U_2 ) )
      return
    end if

    !-- Analyze Unit_1
    Caret = index ( Unit_1, '^' )     
    if ( Caret > 0 ) then
      UnitBase_1 = Unit_1 ( : Caret - 1 )
      Dot = index ( Unit_1 ( Caret + 1 : ), '.' )
      if ( Dot > 0 ) then
        print *, 'ERROR: Real exponents not implemented in Product_MV_MV'
        stop
      else
        read ( Unit_1 ( Caret + 1 : ), fmt = '(i7)') Exponent_1
      end if
    else
      UnitBase_1 = Unit_1
      Exponent_1 = 1
    end if

    !-- Analyze Unit_2
    Caret = index ( Unit_2, '^' )     
    if ( Caret > 0 ) then
      UnitBase_2 = Unit_2 ( : Caret - 1 )
      Dot = index ( Unit_2 ( Caret + 1 : ), '.' )
      if ( Dot > 0 ) then
        print *, 'ERROR: Real exponents not implemented in Product_MV_MV'
        stop
      else
        read ( Unit_2 ( Caret + 1 : ), fmt = '(i7)') Exponent_2
      end if
    else
      UnitBase_2 = Unit_2
      Exponent_2 = 1
    end if

    if ( trim ( UnitBase_1 ) == trim ( UnitBase_2 ) ) then
      Exponent = Exponent_1 + Exponent_2
      if ( Exponent == 0 ) then
        P_U = ''
      else if ( Exponent == 1 ) then
        P_U = trim ( UnitBase_1 )
      else 
        write ( ExponentString, fmt ='(i7)' ) Exponent
        P_U = trim ( UnitBase_1 ) // '^' &
              // trim ( adjustl ( ExponentString ) )
      end if
    else
      P_U = adjustl ( trim ( Unit_1 ) // ' ' // trim ( Unit_2 ) )
    end if

! print*, 'UnitBase_1 ', trim ( UnitBase_1 )
! print*, 'UnitBase_2 ', trim ( UnitBase_2 )
! print*, 'Exponent_1 ', Exponent_1
! print*, 'Exponent_2 ', Exponent_2
! print*, 'Exponent', Exponent

  end function ProductUnit


  recursive function ProductUnit_UCS ( Unit_1, Unit_2 ) result ( P_U )

    character ( *, KBCH ), intent ( in ) :: &
      Unit_1, &
      Unit_2
    character ( LDL, KBCH ) :: &
      P_U

    integer ( KDI ) :: &
      iU_1, iU_2, &
      Exponent_1, &
      Exponent_2, &
      Exponent
    character ( LDL, KBCH ) :: &
      U_1, &
      U_2, &
      UnitBase_1, &
      UnitBase_2!, &
!       ExponentString
    character ( LDL, KBCH ), dimension ( : ), allocatable :: &
      UnitPiece_1, &
      UnitPiece_2
  
    P_U = KBCH_''

    if ( trim ( Unit_1 ) == KBCH_'' .and. trim ( Unit_2 ) == KBCH_'' ) &
      return

    call Split_KBCH ( Unit_1, KBCH_' ', UnitPiece_1 )
    call Split_KBCH ( Unit_2, KBCH_' ', UnitPiece_2 )
    if ( size ( UnitPiece_1 ) > 1 .or. size ( UnitPiece_2 ) > 1 ) then
      do iU_1 = 1, size ( UnitPiece_1 )
        do iU_2 = 1, size ( UnitPiece_2 )
          call AnalyzeUnit &
                 ( UnitPiece_1 ( iU_1 ), UnitBase_1, Exponent_1 )
          if ( index ( UnitPiece_2 ( iU_2 ), trim ( UnitBase_1 ) ) > 0 ) then
            UnitPiece_1 ( iU_1 ) &
              = ProductUnit_UCS ( UnitPiece_1 ( iU_1 ), UnitPiece_2 ( iU_2 ) )
            UnitPiece_2 ( iU_2 ) = KBCH_''
          end if
        end do !-- iU_2
      end do !-- iU_1
      call Join_KBCH ( UnitPiece_1, KBCH_' ', U_1 )
      call Join_KBCH ( UnitPiece_2, KBCH_' ', U_2 )
      P_U = adjustl ( trim ( U_1 ) // KBCH_' ' // trim ( U_2 ) )
      return
    end if

    call AnalyzeUnit ( Unit_1, UnitBase_1, Exponent_1 )
    call AnalyzeUnit ( Unit_2, UnitBase_2, Exponent_2 )

    if ( trim ( UnitBase_1 ) == trim ( UnitBase_2 ) ) then
      Exponent = Exponent_1 + Exponent_2
      if ( Exponent == 0 ) then
        P_U = ''
      else if ( Exponent == 1 ) then
        P_U = trim ( UnitBase_1 )
      else 
        if ( Exponent < 0 ) &
          UnitBase_1 = trim ( UnitBase_1 ) // ExponentMinusChar
        P_U = trim ( UnitBase_1 ) // Exponent_I_Char ( abs ( Exponent ) )
      end if
    else !-- different units
      P_U = adjustl ( trim ( Unit_1 ) // KBCH_' ' // trim ( Unit_2 ) )
    end if

! print*, 'UnitBase_1 ', trim ( UnitBase_1 )
! print*, 'UnitBase_2 ', trim ( UnitBase_2 )
! print*, 'Exponent_1 ', Exponent_1
! print*, 'Exponent_2 ', Exponent_2
! print*, 'Exponent', Exponent

  end function ProductUnit_UCS


  elemental function PowerUnitInteger ( Unit, Exponent ) result ( P_U )

    character ( * ), intent ( in ) :: &
      Unit
    integer ( KDI ), intent ( in ) :: &
      Exponent
    character ( LDL ) :: &
      P_U
    
    integer ( KDI ) :: &
      iU, &  !-- iUnit
      Caret, &
      Dot, &
      OldIntegerExponent
    real ( KDR ) :: &
      OldRealExponent
    character ( LDL ) :: &
      ExponentString, &
      Scratch
    character ( LDL ), dimension ( : ), allocatable :: &
      UnitPiece
      
    P_U = ''

    if ( Unit == '' ) &
      return

    call Split ( Unit, ' ', UnitPiece )
    
    do iU = 1, size ( UnitPiece )

      Scratch = UnitPiece ( iU )
      Caret = index ( Scratch, '^' )
      
      if ( Caret > 0 ) then
        Dot = index ( Scratch ( Caret + 1 : ), '.' )
        if ( Dot > 0 ) then
          read ( Scratch ( Caret + 1 : ), fmt = '(f7.4)' ) OldRealExponent
          write &
            ( ExponentString, fmt = '(f7.4)' )  ( OldRealExponent * Exponent )
        else
          read ( Scratch ( Caret + 1 : ), fmt = '(i7)') &
               OldIntegerExponent
          write &
            ( ExponentString, fmt ='(i7)' ) ( OldIntegerExponent *  Exponent )
        end if
        UnitPiece ( iU ) &
          = trim ( Scratch ( : Caret ) ) // trim ( adjustl ( ExponentString ) )
      else
        if ( Exponent == 1 ) exit
        write ( ExponentString, fmt = '(i7)' ) Exponent
        UnitPiece ( iU ) = &
          trim ( Scratch ) // '^' // trim ( adjustl ( ExponentString ) )
      end if

    end do 

    call Join ( UnitPiece, ' ', P_U )
      
  end function PowerUnitInteger


  elemental function PowerUnitInteger_UCS ( Unit, Exponent ) result ( P_U )

    character ( *, KBCH ), intent ( in ) :: &
      Unit
    integer ( KDI ), intent ( in ) :: &
      Exponent
    character ( LDL, KBCH ) :: &
      P_U
    
    integer ( KDI ) :: &
      iU, &   !-- iUnit
      OldExponent, &
      NewExponent
    character ( LDL, KBCH ) :: &
      UnitBase
    character ( LDL, KBCH ), dimension ( : ), allocatable :: &
      UnitPiece
      
    P_U = KBCH_''

    if ( Unit == KBCH_'' ) &
      return

    call Split_KBCH ( Unit, KBCH_' ', UnitPiece )
    
    do iU = 1, size ( UnitPiece )

      call AnalyzeUnit ( UnitPiece ( iU ), UnitBase, OldExponent )

      NewExponent = OldExponent * Exponent
  
      if ( NewExponent /= 1 ) then
        if ( NewExponent < 0 ) &
          UnitBase = trim ( UnitBase ) // ExponentMinusChar
        UnitPiece ( iU ) &
          = trim ( UnitBase ) &
            // Exponent_I_Char ( abs ( NewExponent ) )
      end if

    end do 

    call Join_KBCH ( UnitPiece, KBCH_' ', P_U )
      
  end function PowerUnitInteger_UCS


  ! elemental function PowerUnitReal ( Unit, Exponent ) result ( P_U )

  !   character ( *, KBCH ), intent ( in ) :: &
  !     Unit
  !   real ( KDR ), intent ( in ) :: &
  !     Exponent
  !   character ( LDL, KBCH ) :: &
  !     P_U
    
  !   integer ( KDI ) :: &
  !     iU, &  !-- iUnit
  !     Caret, &
  !     Dot, &
  !     OldIntegerExponent
  !   real ( KDR ) :: &
  !     OldRealExponent
  !   character ( LDL, KBCH ) :: &
  !     ExponentString, &
  !     Scratch
  !   character ( LDL, KBCH ), dimension ( : ), allocatable :: &
  !     UnitPiece
      
  !   P_U = KBCH_''

  !   if ( Unit == KBCH_'' ) return

  !   call Split ( Unit, KBCH_' ', UnitPiece )
    
  !   do iU = 1, size ( UnitPiece )

  !     Scratch = UnitPiece ( iU )
  !     Caret = index ( Scratch, KBCH_'^' )
     
  !     if ( Caret > 0 ) then
  !       Dot = index ( Scratch ( Caret + 1 : ), KBCH_'.' )
  !       if ( Dot > 0 ) then
  !         read ( Scratch ( Caret + 1 : ), fmt = '(f7.4)' ) OldRealExponent
  !         write &
  !           ( ExponentString, fmt = '(f7.4)' )  ( OldRealExponent * Exponent )
  !       else
  !         read ( Scratch ( Caret + 1 : ), fmt = '(i7)') &
  !              OldIntegerExponent
  !         write &
  !           ( ExponentString, fmt ='(f7.4)' ) ( OldIntegerExponent * Exponent )
  !       end if
  !       if ( trim ( adjustl ( ExponentString ) ) == KBCH_'1.0000' ) then
  !         UnitPiece ( iU ) = trim ( Scratch ( : Caret - 1 ) )
  !       else
  !         UnitPiece ( iU ) &
  !           = trim ( Scratch ( : Caret ) ) &
  !             // trim ( adjustl ( ExponentString ) )
  !       end if
  !     else
  !       write ( ExponentString, fmt = '(f7.4)' ) Exponent
  !       if ( trim ( adjustl ( ExponentString ) ) == KBCH_'1.0000' ) then
  !         UnitPiece ( iU ) = trim ( Scratch )
  !       else
  !         UnitPiece ( iU ) = &
  !           trim ( Scratch ) // KBCH_'^' // trim ( adjustl ( ExponentString ) )
  !       end if
  !     end if

  !   end do
    
  !   call Join ( UnitPiece, KBCH_' ', P_U )
      
  ! end function PowerUnitReal


  pure subroutine AnalyzeUnit ( Unit, UnitBase, Exponent )

    character ( *, KBCH ), intent ( in ) :: &
      Unit
    character ( *, KBCH ), intent ( out ) :: &
      UnitBase
    integer ( KDI ), intent ( out ) :: &
      Exponent

    integer ( KDI ) :: &
      ltU  !-- ltUnit
    logical ( KDL ) :: &
      HasExponent, &
      NegativeExponent

    ltU = len_trim ( Unit )

    if ( ltU == 0 ) then
      UnitBase = ''
      Exponent = 0
      return
    end if

    HasExponent = any ( ichar ( Unit ( ltU : ), KBCH ) == Exponent_I )

    if ( HasExponent ) then
      Exponent = index ( ExponentCharacters, Unit ( ltU : ltU ) )
      NegativeExponent &
        = ichar ( Unit ( ltU - 1 : ), KBCH ) == ExponentMinus
      if ( NegativeExponent ) then
        Exponent  =  - Exponent
        UnitBase = Unit ( 1 : ltU - 2 )
      else
        UnitBase = Unit ( 1 : ltU - 1 )
      end if
    else !-- no exponent
      UnitBase = Unit
      Exponent = 1
    end if

  end subroutine AnalyzeUnit


end module MeasuredValue_Form
