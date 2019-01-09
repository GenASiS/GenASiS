!-- This wrapper exists for use by CONSOLE_Singleton.

subroutine ShowCharacter &
             ( Character, Description, IgnorabilityOption, &
               DisplayRankOption, nLeadingLinesOption, nTrailingLinesOption )

    use Specifiers
    use Show_Command

    implicit none

    character ( * ), intent ( in ) :: &
      Character, &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    call Show ( Character, Description, IgnorabilityOption, &
                DisplayRankOption, nLeadingLinesOption, nTrailingLinesOption )

end subroutine ShowCharacter
