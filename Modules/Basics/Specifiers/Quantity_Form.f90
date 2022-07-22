!-- Quantity_Form handles numbers with labels (typically units) providing
!   means of dealing with units, for example. 

module Quantity_Form
  
  use KIND_DEFAULT_Singleton
  use KIND_BIG_Singleton
  use LEN_DEFAULT_Singleton
  use Split_Command
  use Join_Command
  
  implicit none
  private
  
  type, public :: QuantityForm
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
      Initialize_Q
    procedure, private, pass :: &
      Initialize_Q_Label
    procedure, private, pass :: &
      Initialize_Q_From_Q
    procedure, private, pass :: &
      Initialize_Q_From_Q_Label
    generic :: &
      Initialize &
        => Initialize_Q, Initialize_Q_Label, &
           Initialize_Q_From_Q, Initialize_Q_From_Q_Label
    procedure, private, pass :: &
      Initialize_UCS_Q
    procedure, private, pass :: &
      Initialize_UCS_Q_Label
    procedure, private, pass :: &
      Initialize_UCS_Q_From_Q_Label
    generic :: &
      Initialize_UCS &
        => Initialize_UCS_Q, Initialize_UCS_Q_Label, &
           Initialize_UCS_Q_From_Q_Label
    procedure, private, pass :: &
      Addition_Q_Q
    procedure, private, pass :: &
      Addition_Q_Real
    procedure, private, pass ( Addend ) :: &
      AdditionReal_Q
    procedure, private, pass :: &
      Addition_Q_Integer
    procedure, private, pass ( Addend ) :: &
      AdditionInteger_Q
    generic :: &
      operator ( + ) &
        => Addition_Q_Q, Addition_Q_Real, AdditionReal_Q, &
           Addition_Q_Integer, AdditionInteger_Q
    procedure, private, pass :: &
      Subtraction_Q_Q
    procedure, private, pass :: &
      Subtraction_Q_Real
    procedure, private, pass ( Subtrahend ) :: &
      SubtractionReal_Q
    procedure, private, pass :: &
      Subtraction_Q_Integer
    procedure, private, pass ( Subtrahend ) :: &
      SubtractionInteger_Q
    generic :: &
      operator ( - ) &
        => Subtraction_Q_Q, Subtraction_Q_Real, SubtractionReal_Q, &
           Subtraction_Q_Integer, SubtractionInteger_Q
    procedure, private, pass :: &
      Product_Q_Q
    procedure, private, pass :: &
      Product_Q_Real
    procedure, private, pass ( Multiplicand ) :: &
      ProductReal_Q
    procedure, private, pass :: &
      Product_Q_Integer
    procedure, private, pass ( Multiplicand ) :: &
      ProductInteger_Q
    generic :: &
      operator ( * ) &
        => Product_Q_Q, Product_Q_Real, ProductReal_Q, &
           Product_Q_Integer, ProductInteger_Q
    procedure, private, pass :: &
      Power_Q_Integer
!    procedure, private, pass :: &
!      Power_Q_Real
    generic :: &
      operator ( ** ) => Power_Q_Integer!, Power_Q_Real
    procedure, private, pass :: &
      Quotient_Q_Q
    procedure, private, pass :: &
      Quotient_Q_Real
    procedure, private, pass ( Divisor ) :: &
      QuotientReal_Q
    procedure, private, pass :: &
      Quotient_Q_Integer
    procedure, private, pass ( Divisor ) :: &
      QuotientInteger_Q
    generic :: &
      operator ( / ) &
        => Quotient_Q_Q, Quotient_Q_Real, QuotientReal_Q, &
           Quotient_Q_Integer, QuotientInteger_Q
    procedure, private, pass ( Q ) :: &
      AssignReal_Q
    generic :: &
      assignment ( = ) => AssignReal_Q
    procedure, private, pass :: &
      EqualTo_Q_Q
    procedure, private, pass :: &
      EqualTo_Q_Real
    procedure, private, pass ( Q ) :: &
      EqualToReal_Q
    procedure, private, pass :: &
      EqualTo_Q_Integer
    procedure, private, pass ( Q ) :: &
      EqualToInteger_Q
    generic :: &
      operator ( == ) &
        => EqualTo_Q_Q, EqualTo_Q_Real, EqualToReal_Q, &
           EqualTo_Q_Integer, EqualToInteger_Q
    procedure, private, pass :: &
      NotEqualTo_Q_Q
    procedure, private, pass :: &
      NotEqualTo_Q_Real
    procedure, private, pass ( Q ) :: &
      NotEqualToReal_Q
    procedure, private, pass :: &
      NotEqualTo_Q_Integer
    procedure, private, pass ( Q ) :: &
      NotEqualToInteger_Q
    generic :: &
      operator ( /= ) &
        => NotEqualTo_Q_Q, NotEqualTo_Q_Real, NotEqualToReal_Q, &
           NotEqualTo_Q_Integer, NotEqualToInteger_Q
    procedure, private, pass :: &
      GreaterThan_Q_Q
    procedure, private, pass :: &
      GreaterThan_Q_Real
    procedure, private, pass ( Q ) :: &
      GreaterThanReal_Q
    procedure, private, pass :: &
      GreaterThan_Q_Integer
    procedure, private, pass ( Q ) :: &
      GreaterThanInteger_Q
    generic :: &
      operator ( > ) &
        => GreaterThan_Q_Q, GreaterThan_Q_Real, GreaterThanReal_Q, &
           GreaterThan_Q_Integer, GreaterThanInteger_Q
    procedure, private, pass :: &
      LessThan_Q_Q
    procedure, private, pass :: &
      LessThan_Q_Real
    procedure, private, pass ( Q ) :: &
      LessThanReal_Q
    procedure, private, pass :: &
      LessThan_Q_Integer
    procedure, private, pass ( Q ) :: &
      LessThanInteger_Q
    generic :: &
      operator ( < ) &
        => LessThan_Q_Q, LessThan_Q_Real, LessThanReal_Q, &
           LessThan_Q_Integer, LessThanInteger_Q
    procedure, private, pass :: &
      GreaterThanEqualTo_Q_Q
    procedure, private, pass :: &
      GreaterThanEqualTo_Q_Real
    procedure, private, pass ( Q ) :: &
      GreaterThanEqualToReal_Q
    procedure, private, pass :: &
      GreaterThanEqualTo_Q_Integer
    procedure, private, pass ( Q ) :: &
      GreaterThanEqualToInteger_Q
    generic :: &
      operator ( >= ) &
        => GreaterThanEqualTo_Q_Q, GreaterThanEqualTo_Q_Real, &
           GreaterThanEqualToReal_Q, GreaterThanEqualTo_Q_Integer, &
           GreaterThanEqualToInteger_Q
    procedure, private, pass :: &
      LessThanEqualTo_Q_Q
    procedure, private, pass :: &
      LessThanEqualTo_Q_Real
    procedure, private, pass ( Q ) :: &
      LessThanEqualToReal_Q
    procedure, private, pass :: &
      LessThanEqualTo_Q_Integer
    procedure, private, pass ( Q ) :: &
      LessThanEqualToInteger_Q
    generic :: &
      operator ( <= ) &
        => LessThanEqualTo_Q_Q, LessThanEqualTo_Q_Real, &
           LessThanEqualToReal_Q, LessThanEqualTo_Q_Integer, &
           LessThanEqualToInteger_Q
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
  end type QuantityForm

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


  subroutine Initialize_Q ( Q, Unit, Number )
    
    class ( QuantityForm ), intent ( inout ) :: &
      Q
    character ( * ), intent ( in ) :: &
      Unit
    real ( KDR ), intent ( in ) :: &
      Number
      
    Q % Number   = Number
    Q % Unit     = Unit
    Q % Unit_UCS = Unit

  end subroutine Initialize_Q
  
  
  subroutine Initialize_Q_Label ( Q, Label, Unit, Number )
    
    class ( QuantityForm ), intent ( inout ) :: &
      Q
    character ( * ), intent ( in ) :: &
      Label, &
      Unit
    real ( KDR ), intent ( in ) :: &
      Number
      
    Q % Number    = Number
    Q % Unit      = Unit
    Q % Label     = Label
    Q % Unit_UCS  = Unit
    Q % Label_UCS = Label

  end subroutine Initialize_Q_Label
  
  
  subroutine Initialize_Q_From_Q ( Q_Target, Q_Source )
    
    class ( QuantityForm ), intent ( inout ) :: &
      Q_Target
    type ( QuantityForm ), intent ( in ) :: &
      Q_Source
            
    Q_Target % Number    = Q_Source % Number
    Q_Target % Unit      = Q_Source % Unit
    Q_Target % Label     = Q_Source % Label
    Q_Target % Unit_UCS  = Q_Source % Unit_UCS
    Q_Target % Label_UCS = Q_Source % Label_UCS
  
  end subroutine Initialize_Q_From_Q
  
  
  subroutine Initialize_Q_From_Q_Label ( Q_Target, Q_Source, Label )
    
    class ( QuantityForm ), intent ( inout ) :: &
      Q_Target
    type ( QuantityForm ), intent ( in ) :: &
      Q_Source
    character ( * ), intent ( in ) :: &
      Label
            
    Q_Target % Number    = Q_Source % Number
    Q_Target % Unit      = Q_Source % Unit
    Q_Target % Label     = Label
    Q_Target % Unit_UCS  = Q_Source % Unit_UCS
    Q_Target % Label_UCS = Label
  
  end subroutine Initialize_Q_From_Q_Label
  
  
  subroutine Initialize_UCS_Q ( Q, Unit_UCS, Unit, Number )
    
    class ( QuantityForm ), intent ( inout ) :: &
      Q
    character ( *, KBCH ), intent ( in ) :: &
      Unit_UCS
    character ( * ), intent ( in ) :: &
      Unit
    real ( KDR ), intent ( in ) :: &
      Number
      
    Q % Number   = Number
    Q % Unit     = Unit
    Q % Unit_UCS = Unit_UCS

  end subroutine Initialize_UCS_Q
  
  
  subroutine Initialize_UCS_Q_Label &
               ( Q, Label_UCS, Unit_UCS, Label, Unit, Number )
    
    class ( QuantityForm ), intent ( inout ) :: &
      Q
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
      
    Q % Number    = Number
    Q % Unit      = Unit
    Q % Label     = Label
    Q % Unit_UCS  = Unit_UCS
    Q % Label_UCS = Label_UCS

  end subroutine Initialize_UCS_Q_Label
  
  
  subroutine Initialize_UCS_Q_From_Q_Label &
               ( Q_Target, Q_Source, Label_UCS, Label )
    
    class ( QuantityForm ), intent ( inout ) :: &
      Q_Target
    type ( QuantityForm ), intent ( in ) :: &
      Q_Source
    character ( *, KBCH ), intent ( in ) :: &
      Label_UCS
    character ( * ), intent ( in ) :: &
      Label
            
    Q_Target % Number    = Q_Source % Number
    Q_Target % Unit      = Q_Source % Unit
    Q_Target % Label     = Label
    Q_Target % Unit_UCS  = Q_Source % Unit_UCS
    Q_Target % Label_UCS = Label_UCS
  
  end subroutine Initialize_UCS_Q_From_Q_Label
  
  
  elemental function Addition_Q_Q ( Augend, Addend ) result ( A_Q_Q )

    class ( QuantityForm ), intent ( in ) :: &
      Augend
    type ( QuantityForm ), intent ( in ) :: &
      Addend
    type ( QuantityForm ) :: &
      A_Q_Q
    
      A_Q_Q % Number = Augend % Number + Addend % Number
    
      if ( trim ( Augend % Unit ) == trim ( Addend % Unit ) ) then
        A_Q_Q % Unit     = Augend % Unit
        A_Q_Q % Unit_UCS = Augend % Unit_UCS
      else
        A_Q_Q % Unit     = 'Undefined'
        A_Q_Q % Unit_UCS = 'Undefined'
      end if
  
      if ( trim ( Augend % Label ) == trim ( Addend % Label ) ) then
        A_Q_Q % Label     = Augend % Label
        A_Q_Q % Label_UCS = Augend % Label_UCS
      else
        A_Q_Q % Label     = 'Undefined'
        A_Q_Q % Label_UCS = 'Undefined'
      end if
  
  end function Addition_Q_Q
  
  
  elemental function Addition_Q_Real &
                       ( Augend, Addend ) result ( A_Q_R )
  
    class ( QuantityForm ), intent ( in ) :: &
      Augend
    real ( KDR ), intent ( in ) :: &
      Addend
    real ( KDR ) :: &
      A_Q_R
      
    A_Q_R = Augend % Number + Addend
      
  end function Addition_Q_Real
   
  
  elemental function AdditionReal_Q &
                       ( Augend, Addend ) result ( A_R_Q )
  
    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "Augend + Addend" 

    real ( KDR ), intent ( in ) :: &
      Augend
    class ( QuantityForm ), intent ( in ) :: &
      Addend
    real ( KDR ) :: &
      A_R_Q
      
    A_R_Q = Augend + Addend % Number
      
  end function AdditionReal_Q
   
  
  elemental function Addition_Q_Integer &
                       ( Augend, Addend ) result ( A_Q_I )
  
    class ( QuantityForm ), intent ( in ) :: &
      Augend
    integer ( KDI ), intent ( in ) :: &
      Addend
    real ( KDR ) :: &
      A_Q_I
      
    A_Q_I = Augend % Number + real ( Addend, KDR )
      
  end function Addition_Q_Integer
   
  
  elemental function AdditionInteger_Q &
                       ( Augend, Addend ) result ( A_I_Q )
  
    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "Augend + Addend" 

    integer ( KDI ), intent ( in ) :: &
      Augend
    class ( QuantityForm ), intent ( in ) :: &
      Addend
    real ( KDR ) :: &
      A_I_Q
      
    A_I_Q = real ( Augend, KDR ) + Addend % Number
      
  end function AdditionInteger_Q
   
  
  elemental function Subtraction_Q_Q &
                       ( Minuend, Subtrahend ) result ( S_Q_Q )
  
    class ( QuantityForm ), intent ( in ) :: &
      Minuend
    class ( QuantityForm ), intent ( in ) :: &
      Subtrahend
    type ( QuantityForm ) :: &
      S_Q_Q
      
    S_Q_Q % Number = Minuend % Number - Subtrahend % Number
      
    if ( trim ( Minuend % Unit ) == trim ( Subtrahend % Unit ) ) then
      S_Q_Q % Unit     = Minuend % Unit
      S_Q_Q % Unit_UCS = Minuend % Unit_UCS
    else
      S_Q_Q % Unit     = 'Undefined'
      S_Q_Q % Unit_UCS = 'Undefined'
    end if
    
    if ( trim ( Minuend % Label ) == trim ( Subtrahend % Label ) ) then
      S_Q_Q % Label     = Minuend % Label
      S_Q_Q % Label_UCS = Minuend % Label_UCS
    else
      S_Q_Q % Label     = 'Undefined'
      S_Q_Q % Label_UCS = 'Undefined'
    end if
    
  end function Subtraction_Q_Q
  
  
  elemental function Subtraction_Q_Real &
                       ( Minuend, Subtrahend ) result ( S_Q_R )
  
    class ( QuantityForm ), intent ( in ) :: &
      Minuend
    real ( KDR ), intent ( in ) :: &
      Subtrahend
    real ( KDR ) :: &
      S_Q_R
      
    S_Q_R = Minuend % Number - Subtrahend
      
  end function Subtraction_Q_Real
   
  
  elemental function SubtractionReal_Q &
                       ( Minuend, Subtrahend ) result ( S_R_Q )
  
    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "Minuend - Subtrahend" 

    real ( KDR ), intent ( in ) :: &
      Minuend
    class ( QuantityForm ), intent ( in ) :: &
      Subtrahend
    real ( KDR ) :: &
      S_R_Q
      
    S_R_Q = Minuend - Subtrahend % Number
      
  end function SubtractionReal_Q
   
  
  elemental function Subtraction_Q_Integer &
                       ( Minuend, Subtrahend ) result ( S_Q_I )
  
    class ( QuantityForm ), intent ( in ) :: &
      Minuend
    integer ( KDI ), intent ( in ) :: &
      Subtrahend
    real ( KDR ) :: &
      S_Q_I
      
    S_Q_I = Minuend % Number - real ( Subtrahend, KDR )
      
  end function Subtraction_Q_Integer
   
  
  elemental function SubtractionInteger_Q &
                       ( Minuend, Subtrahend ) result ( S_I_Q )
  
    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "Minuend - Subtrahend" 

    integer ( KDI ), intent ( in ) :: &
      Minuend
    class ( QuantityForm ), intent ( in ) :: &
      Subtrahend
    real ( KDR ) :: &
      S_I_Q
      
    S_I_Q = real ( Minuend, KDR ) - Subtrahend % Number
      
  end function SubtractionInteger_Q
   
  
  impure elemental function Product_Q_Q &
                     ( Multiplier, Multiplicand ) result ( P_Q_Q )
  
    class ( QuantityForm ), intent ( in ) :: &
      Multiplier, &
      Multiplicand
    type ( QuantityForm ) :: &
      P_Q_Q

    P_Q_Q % Number = Multiplier % Number * Multiplicand % Number

    P_Q_Q % Unit &
      = P_Q_Q % ProductUnit ( Multiplier % Unit, Multiplicand % Unit )
    P_Q_Q % Label &
      = P_Q_Q % ProductUnit ( Multiplier % Label, Multiplicand % Label )

    if ( KBCH > KDCH ) then
      P_Q_Q % Unit_UCS &
        = P_Q_Q % ProductUnit_UCS &
            ( Multiplier % Unit_UCS, Multiplicand % Unit_UCS )
      P_Q_Q % Label_UCS &
        = P_Q_Q % ProductUnit_UCS &
            ( Multiplier % Label_UCS, Multiplicand % Label_UCS )
    else
      P_Q_Q % Unit_UCS &
        = P_Q_Q % ProductUnit &
            ( Multiplier % Unit, Multiplicand % Unit )
      P_Q_Q % Label_UCS &
        = P_Q_Q % ProductUnit &
            ( Multiplier % Label, Multiplicand % Label )
    end if

  end function Product_Q_Q
  
  
  elemental function Product_Q_Real &
                       ( Multiplier, Multiplicand ) result ( P_Q_R )
    
    class ( QuantityForm ), intent ( in ) :: &
      Multiplier
    real ( KDR ), intent ( in ) :: &
      Multiplicand
    type ( QuantityForm ) :: &
      P_Q_R
    
    P_Q_R % Number    = Multiplier % Number * Multiplicand
    P_Q_R % Unit      = Multiplier % Unit
    P_Q_R % Label     = Multiplier % Label
    P_Q_R % Unit_UCS  = Multiplier % Unit_UCS
    P_Q_R % Label_UCS = Multiplier % Label_UCS
  
  end function Product_Q_Real
  

  elemental function ProductReal_Q &
                        ( Multiplier, Multiplicand ) result ( P_R_Q )
    
    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "Multiplier * Multiplicand" 

    real ( KDR ), intent ( in ) :: &
      Multiplier
    class ( QuantityForm ), intent ( in ) :: &
      Multiplicand
    type ( QuantityForm ) :: &
      P_R_Q
    
    P_R_Q % Number    = Multiplier * Multiplicand % Number
    P_R_Q % Unit      = Multiplicand % Unit
    P_R_Q % Label     = Multiplicand % Label
    P_R_Q % Unit_UCS  = Multiplicand % Unit_UCS
    P_R_Q % Label_UCS = Multiplicand % Label_UCS
    
  end function ProductReal_Q
  
  
  elemental function Product_Q_Integer &
                       ( Multiplier, Multiplicand ) result ( P_Q_I )
    
    class ( QuantityForm ), intent ( in ) :: &
      Multiplier
    integer ( KDI ), intent ( in ) :: &
      Multiplicand
    type ( QuantityForm ) :: &
      P_Q_I
    
    P_Q_I % Number    = Multiplier % Number * real ( Multiplicand, KDR )
    P_Q_I % Unit      = Multiplier % Unit
    P_Q_I % Label     = Multiplier % Label
    P_Q_I % Unit_UCS  = Multiplier % Unit_UCS
    P_Q_I % Label_UCS = Multiplier % Label_UCS
  
  end function Product_Q_Integer
  

  elemental function ProductInteger_Q &
                        ( Multiplier, Multiplicand ) result ( P_I_Q )
    
    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "Multiplier * Multiplicand" 

    integer ( KDI ), intent ( in ) :: &
      Multiplier
    class ( QuantityForm ), intent ( in ) :: &
      Multiplicand
    type ( QuantityForm ) :: &
      P_I_Q
    
    P_I_Q % Number    = real ( Multiplier, KDR ) * Multiplicand % Number
    P_I_Q % Unit      = Multiplicand % Unit
    P_I_Q % Label     = Multiplicand % Label
    P_I_Q % Unit_UCS  = Multiplicand % Unit_UCS
    P_I_Q % Label_UCS = Multiplicand % Label_UCS
    
  end function ProductInteger_Q
  
  
  elemental function Power_Q_Integer ( Base, Exponent ) result ( P_Q_I )
    
    class ( QuantityForm ), intent ( in ) :: &
      Base
    integer ( KDI ), intent ( in ) :: &
      Exponent
    type ( QuantityForm ) :: &
      P_Q_I
    
    P_Q_I % Number = Base % Number ** Exponent

    P_Q_I % Unit  = P_Q_I % PowerUnitInteger ( Base % Unit, Exponent )
    P_Q_I % Label = P_Q_I % PowerUnitInteger ( Base % Label, Exponent )

    if ( KBCH > KDCH ) then
      P_Q_I % Unit_UCS &
        = P_Q_I % PowerUnitInteger_UCS ( Base % Unit_UCS, Exponent )
      P_Q_I % Label_UCS &
        = P_Q_I % PowerUnitInteger_UCS ( Base % Label_UCS, Exponent )
    else
      P_Q_I % Unit_UCS &
        = P_Q_I % PowerUnitInteger ( Base % Unit, Exponent )
      P_Q_I % Label_UCS &
        = P_Q_I % PowerUnitInteger ( Base % Label, Exponent )
    end if

  end function Power_Q_Integer
  

  ! elemental function Power_Q_Real ( Base, Exponent ) result ( P_Q_R )
    
  !   class ( QuantityForm ), intent ( in ) :: &
  !     Base
  !   real ( KDR ), intent ( in ) :: &
  !     Exponent
  !   type ( QuantityForm ) :: &
  !     P_Q_R
    
  !   P_Q_R % Number = Base % Number ** Exponent
    
  !   P_Q_R % Unit = P_Q_R % PowerUnit ( Base % Unit, Exponent )
    
  !   P_Q_R % Label = P_Q_R % PowerUnit ( Base % Label, Exponent )
    
  ! end function Power_Q_Real
  

  impure elemental function Quotient_Q_Q ( Dividend, Divisor ) &
                              result ( Q_Q_Q )
    
    class ( QuantityForm ), intent ( in ) :: &
      Dividend, &
      Divisor
    type ( QuantityForm ) :: &
      P_Q_I, &
      Q_Q_Q
    
    !-- FIXME: Workaround for GCC 5.1 bug 66257
    !P_Q_I  = Divisor ** ( -1 )
    
    !Q_Q_Q = Dividend * P_Q_I
    Q_Q_Q = Dividend * Divisor ** ( -1 )

  end function Quotient_Q_Q
  
  
  elemental function Quotient_Q_Real ( Dividend, Divisor ) result ( Q_Q_R )
    
    class ( QuantityForm ), intent ( in ) :: &
      Dividend
    real ( KDR ), intent ( in ) :: &
      Divisor
    type ( QuantityForm ) :: &
      Q_Q_R
    
    Q_Q_R = Dividend * Divisor ** ( -1 )

  end function Quotient_Q_Real
  
  
  elemental function QuotientReal_Q ( Dividend, Divisor ) result ( Q_R_Q )
    
    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "Dividend / Divisor" 

    real ( KDR ), intent ( in ) :: &
      Dividend
    class ( QuantityForm ), intent ( in ) :: &
      Divisor
    type ( QuantityForm ) :: &
      P_Q_I, &
      Q_R_Q
    
    !-- FIXME: Workaround for GCC 5.1 bug 66257
    P_Q_I = Divisor ** ( -1 )
    
    Q_R_Q = Dividend * P_Q_I

  end function QuotientReal_Q
  
  
  elemental function Quotient_Q_Integer ( Dividend, Divisor ) result ( Q_Q_I )
    
    class ( QuantityForm ), intent ( in ) :: &
      Dividend
    integer ( KDI ), intent ( in ) :: &
      Divisor
    type ( QuantityForm ) :: &
      Q_Q_I
    
    Q_Q_I = Dividend * real ( Divisor, KDR ) ** ( -1 )

  end function Quotient_Q_Integer
  
  
  elemental function QuotientInteger_Q ( Dividend, Divisor ) result ( Q_I_Q )
    
    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "Dividend / Divisor" 

    integer ( KDI ), intent ( in ) :: &
      Dividend
    class ( QuantityForm ), intent ( in ) :: &
      Divisor
    type ( QuantityForm ) :: &
      P_Q_I, &
      Q_I_Q
    
    !-- FIXME: Workaround for GCC 5.1 bug 66257
    P_Q_I = Divisor ** ( -1 )
    
    Q_I_Q = real ( Dividend, KDR ) * P_Q_I

  end function QuotientInteger_Q
  
  
  elemental subroutine AssignReal_Q ( A, Q )

    real ( KDR ), intent ( inout ) :: &
      A
    class ( QuantityForm ), intent ( in ) :: &
      Q

    A = Q % Number

  end subroutine AssignReal_Q


  elemental function EqualTo_Q_Q ( Q_1, Q_2 ) result ( ET )

    class ( QuantityForm ), intent ( in ) :: &
      Q_1, Q_2
    logical ( KDL ) :: &
      ET

    ET = .false.

    if ( Q_1 % Unit == Q_2 % Unit &
         .and. Q_1 % Number == Q_2 % Number ) &
      ET = .true.

  end function EqualTo_Q_Q


  elemental function EqualTo_Q_Real ( Q, A ) result ( ET )

    class ( QuantityForm ), intent ( in ) :: &
      Q
    real ( KDR ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      ET

    ET = .false.

    if ( Q % Number == A ) ET = .true.

  end function EqualTo_Q_Real


  elemental function EqualToReal_Q ( A, Q ) result ( ET )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A == Q" 

    real ( KDR ), intent ( in ) :: &
      A
    class ( QuantityForm ), intent ( in ) :: &
      Q
    logical ( KDL ) :: &
      ET

    ET = .false.

    if ( A == Q % Number ) ET = .true.

  end function EqualToReal_Q


  elemental function EqualTo_Q_Integer ( Q, A ) result ( ET )

    class ( QuantityForm ), intent ( in ) :: &
      Q
    integer ( KDI ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      ET

    ET = .false.

    if ( Q % Number == real ( A, KDR ) ) ET = .true.

  end function EqualTo_Q_Integer


  elemental function EqualToInteger_Q ( A, Q ) result ( ET )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A == Q" 

    integer ( KDI ), intent ( in ) :: &
      A
    class ( QuantityForm ), intent ( in ) :: &
      Q
    logical ( KDL ) :: &
      ET

    ET = .false.

    if ( real ( A, KDR ) == Q % Number ) ET = .true.

  end function EqualToInteger_Q


  elemental function NotEqualTo_Q_Q ( Q_1, Q_2 ) result ( NET )

    class ( QuantityForm ), intent ( in ) :: &
      Q_1, Q_2
    logical ( KDL ) :: &
      NET

    NET = .false.

    if ( Q_1 % Unit /= Q_2 % Unit &
         .or. Q_1 % Number /= Q_2 % Number ) &
      NET = .true.

  end function NotEqualTo_Q_Q


  elemental function NotEqualTo_Q_Real ( Q, A ) result ( NET )

    class ( QuantityForm ), intent ( in ) :: &
      Q
    real ( KDR ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      NET

    NET = .false.

    if ( Q % Number /= A ) NET = .true.

  end function NotEqualTo_Q_Real


  elemental function NotEqualToReal_Q ( A, Q ) result ( NET )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A /= Q" 

    real ( KDR ), intent ( in ) :: &
      A
    class ( QuantityForm ), intent ( in ) :: &
      Q
    logical ( KDL ) :: &
      NET

    NET = .false.

    if ( A /= Q % Number ) NET = .true.

  end function NotEqualToReal_Q


  elemental function NotEqualTo_Q_Integer ( Q, A ) result ( NET )

    class ( QuantityForm ), intent ( in ) :: &
      Q
    integer ( KDI ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      NET

    NET = .false.

    if ( Q % Number /= real ( A, KDR ) ) NET = .true.

  end function NotEqualTo_Q_Integer


  elemental function NotEqualToInteger_Q ( A, Q ) result ( NET )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A /= Q" 

    integer ( KDI ), intent ( in ) :: &
      A
    class ( QuantityForm ), intent ( in ) :: &
      Q
    logical ( KDL ) :: &
      NET

    NET = .false.

    if ( real ( A, KDR ) /= Q % Number ) NET = .true.

  end function NotEqualToInteger_Q


  elemental function GreaterThan_Q_Q ( Q_1, Q_2 ) result ( GT )

    class ( QuantityForm ), intent ( in ) :: &
      Q_1, Q_2
    logical ( KDL ) :: &
      GT

    GT = .false.

    if ( Q_1 % Unit == Q_2 % Unit &
         .and. Q_1 % Number > Q_2 % Number ) &
      GT = .true.

  end function GreaterThan_Q_Q


  elemental function GreaterThan_Q_Real ( Q, A ) result ( GT )

    class ( QuantityForm ), intent ( in ) :: &
      Q
    real ( KDR ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      GT

    GT = .false.

    if ( Q % Number > A ) GT = .true.

  end function GreaterThan_Q_Real


  elemental function GreaterThanReal_Q ( A, Q ) result ( GT )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A > Q" 

    real ( KDR ), intent ( in ) :: &
      A
    class ( QuantityForm ), intent ( in ) :: &
      Q
    logical ( KDL ) :: &
      GT

    GT = .false.

    if ( A > Q % Number ) GT = .true.

  end function GreaterThanReal_Q


  elemental function GreaterThan_Q_Integer ( Q, A ) result ( GT )

    class ( QuantityForm ), intent ( in ) :: &
      Q
    integer ( KDI ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      GT

    GT = .false.

    if ( Q % Number > real ( A, KDR ) ) GT = .true.

  end function GreaterThan_Q_Integer


  elemental function GreaterThanInteger_Q ( A, Q ) result ( GT )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A > Q" 

    integer ( KDI ), intent ( in ) :: &
      A
    class ( QuantityForm ), intent ( in ) :: &
      Q
    logical ( KDL ) :: &
      GT

    GT = .false.

    if ( real ( A, KDR ) > Q % Number ) GT = .true.

  end function GreaterThanInteger_Q


  elemental function LessThan_Q_Q ( Q_1, Q_2 ) result ( LT )

    class ( QuantityForm ), intent ( in ) :: &
      Q_1, Q_2
    logical ( KDL ) :: &
      LT

    LT = .false.

    if ( Q_1 % Unit == Q_2 % Unit &
         .and. Q_1 % Number < Q_2 % Number ) &
      LT = .true.

  end function LessThan_Q_Q


  elemental function LessThan_Q_Real ( Q, A ) result ( LT )

    class ( QuantityForm ), intent ( in ) :: &
      Q
    real ( KDR ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      LT

    LT = .false.

    if ( Q % Number < A ) LT = .true.

  end function LessThan_Q_Real


  elemental function LessThanReal_Q ( A, Q ) result ( LT )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A < Q" 

    real ( KDR ), intent ( in ) :: &
      A
    class ( QuantityForm ), intent ( in ) :: &
      Q
    logical ( KDL ) :: &
      LT

    LT = .false.

    if ( A < Q % Number ) LT = .true.

  end function LessThanReal_Q


  elemental function LessThan_Q_Integer ( Q, A ) result ( LT )

    class ( QuantityForm ), intent ( in ) :: &
      Q
    integer ( KDI ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      LT

    LT = .false.

    if ( Q % Number < real ( A, KDR ) ) LT = .true.

  end function LessThan_Q_Integer


  elemental function LessThanInteger_Q ( A, Q ) result ( LT )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A < Q" 

    integer ( KDI ), intent ( in ) :: &
      A
    class ( QuantityForm ), intent ( in ) :: &
      Q
    logical ( KDL ) :: &
      LT

    LT = .false.

    if ( real ( A, KDR ) < Q % Number ) LT = .true.

  end function LessThanInteger_Q


  elemental function GreaterThanEqualTo_Q_Q ( Q_1, Q_2 ) result ( GTE )

    class ( QuantityForm ), intent ( in ) :: &
      Q_1, Q_2
    logical ( KDL ) :: &
      GTE

    GTE = .false.

    if ( Q_1 % Unit == Q_2 % Unit &
         .and. Q_1 % Number >= Q_2 % Number ) &
      GTE = .true.

  end function GreaterThanEqualTo_Q_Q


  elemental function GreaterThanEqualTo_Q_Real ( Q, A ) result ( GTE )

    class ( QuantityForm ), intent ( in ) :: &
      Q
    real ( KDR ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      GTE

    GTE = .false.

    if ( Q % Number >= A ) GTE = .true.

  end function GreaterThanEqualTo_Q_Real


  elemental function GreaterThanEqualToReal_Q ( A, Q ) result ( GTE )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A >= Q" 

    real ( KDR ), intent ( in ) :: &
      A
    class ( QuantityForm ), intent ( in ) :: &
      Q
    logical ( KDL ) :: &
      GTE

    GTE = .false.

    if ( A >= Q % Number ) GTE = .true.

  end function GreaterThanEqualToReal_Q


  elemental function GreaterThanEqualTo_Q_Integer ( Q, A ) result ( GTE )

    class ( QuantityForm ), intent ( in ) :: &
      Q
    integer ( KDI ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      GTE

    GTE = .false.

    if ( Q % Number >= real ( A, KDR ) ) GTE = .true.

  end function GreaterThanEqualTo_Q_Integer


  elemental function GreaterThanEqualToInteger_Q ( A, Q ) result ( GTE )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A >= Q" 

    integer ( KDI ), intent ( in ) :: &
      A
    class ( QuantityForm ), intent ( in ) :: &
      Q
    logical ( KDL ) :: &
      GTE

    GTE = .false.

    if ( real ( A, KDR ) >= Q % Number ) GTE = .true.

  end function GreaterThanEqualToInteger_Q


  elemental function LessThanEqualTo_Q_Q ( Q_1, Q_2 ) result ( LTE )

    class ( QuantityForm ), intent ( in ) :: &
      Q_1, Q_2
    logical ( KDL ) :: &
      LTE

    LTE = .false.

    if ( Q_1 % Unit == Q_2 % Unit &
         .and. Q_1 % Number <= Q_2 % Number ) &
      LTE = .true.

  end function LessThanEqualTo_Q_Q


  elemental function LessThanEqualTo_Q_Real ( Q, A ) result ( LTE )

    class ( QuantityForm ), intent ( in ) :: &
      Q
    real ( KDR ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      LTE

    LTE = .false.

    if ( Q % Number <= A ) LTE = .true.

  end function LessThanEqualTo_Q_Real


  elemental function LessThanEqualToReal_Q ( A, Q ) result ( LTE )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A <= Q" 

    real ( KDR ), intent ( in ) :: &
      A
    class ( QuantityForm ), intent ( in ) :: &
      Q
    logical ( KDL ) :: &
      LTE

    LTE = .false.

    if ( A <= Q % Number ) LTE = .true.

  end function LessThanEqualToReal_Q


  elemental function LessThanEqualTo_Q_Integer ( Q, A ) result ( LTE )

    class ( QuantityForm ), intent ( in ) :: &
      Q
    integer ( KDI ), intent ( in ) :: &
      A
    logical ( KDL ) :: &
      LTE

    LTE = .false.

    if ( Q % Number <= real ( A, KDR ) ) LTE = .true.

  end function LessThanEqualTo_Q_Integer


  elemental function LessThanEqualToInteger_Q ( A, Q ) result ( LTE )

    !-- Convention on argument ordering is broken here to maintain 
    !   correspondence with order of entities in the expression
    !   "A <= Q" 

    integer ( KDI ), intent ( in ) :: &
      A
    class ( QuantityForm ), intent ( in ) :: &
      Q
    logical ( KDL ) :: &
      LTE

    LTE = .false.

    if ( real ( A, KDR ) <= Q % Number ) LTE = .true.

  end function LessThanEqualToInteger_Q


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
        print *, 'ERROR: Real exponents not implemented in Product_Q_Q'
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
        print *, 'ERROR: Real exponents not implemented in Product_Q_Q'
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


end module Quantity_Form
