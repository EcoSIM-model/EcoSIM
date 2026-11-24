module StringToolsMod
implicit none

  private
  save
  character(len=*), parameter :: mod_filename = &
  __FILE__
  
  public :: parse_var_val_string
  public :: are_strings_equal_icase
  public :: extract_number_and_unit
  public :: is_substring_present
  public :: to_lower_string
  public :: count_delimited_items
  public :: string_compaction
contains

  SUBROUTINE parse_var_val_string(input_string, var_array, val_array, num_pairs_found)
  !
  !Description
  !> @brief Parses a key-value string separated by semicolons.
  !> @details
  !>   Takes an input string of the format "key1[val1];key2[val2]"
  !>   and separates it into two parallel arrays, one for keys and one for values.
  !>
  !> @param[in] input_string The string to be parsed.
  !> @param[out] var_array A character array to store the extracted 'var' strings.
  !> @param[out] val_array A character array to store the extracted 'val' strings.
  !> @param[out] num_pairs_found The total number of valid pairs found.

    IMPLICIT NONE

    ! --- Arguments ---
    CHARACTER(LEN=*), INTENT(IN) :: input_string
    ! Output arrays: Assumed-size. The caller's array must be large enough.
    ! The string length (e.g., 100) is defined by the caller's declaration.
    CHARACTER(LEN=*), DIMENSION(:), INTENT(OUT) :: var_array, val_array
    INTEGER, INTENT(OUT) :: num_pairs_found

    ! --- Local Variables ---
    CHARACTER(LEN=LEN(input_string)) :: remaining_string
    CHARACTER(LEN=LEN(input_string)) :: current_pair
    INTEGER :: pair_index
    INTEGER :: delim_pos
    INTEGER :: bracket_open_pos, bracket_close_pos
    INTEGER :: max_pairs

    ! --- Implementation ---

    ! Initialize outputs
    num_pairs_found = 0
    var_array = ""  ! Clear the output array
    val_array = ""
    
    ! Get the max number of pairs we can store from the provided array size
    max_pairs = SIZE(var_array)
    
    ! Use a temporary, adjustable string
    remaining_string = TRIM(input_string)

    pair_index = 0
    
    DO WHILE (LEN_TRIM(remaining_string) > 0)
      pair_index = pair_index + 1
      
      ! Check if we've exceeded the array bounds
      IF (pair_index > max_pairs) THEN
          WRITE(*, '(A)') "Warning: More pairs found than array size. Stopping."
          pair_index = pair_index - 1 ! Don't count this last one
          EXIT
      END IF
      
      ! Find the next pair delimiter (;)
      delim_pos = INDEX(remaining_string, ';')
      
      IF (delim_pos == 0) THEN
        ! This is the last pair (or the only pair)
        current_pair = TRIM(remaining_string)
        remaining_string = "" ! Signal to exit the loop
      ELSE
        ! This is not the last pair
        current_pair = TRIM(remaining_string(:delim_pos-1))
        remaining_string = remaining_string(delim_pos+1:)
      END IF
      
      ! Now, parse the current_pair (e.g., "var1[val1]")
      bracket_open_pos = INDEX(current_pair, '[')
      bracket_close_pos = INDEX(current_pair, ']')
      
      ! Basic validation:
      ! 1. Must have an open bracket (not at the start)
      ! 2. Must have a close bracket after the open bracket
      IF (bracket_open_pos > 1 .AND. &
          bracket_close_pos > (bracket_open_pos + 1)) THEN
          
        ! Extract variable
        var_array(pair_index) = TRIM(current_pair(:bracket_open_pos-1))
        
        ! Extract value
        val_array(pair_index) = TRIM(current_pair(bracket_open_pos+1 : bracket_close_pos-1))
        
      ELSE
        ! Malformed pair, skip it
        IF (LEN_TRIM(current_pair) > 0) THEN
          WRITE(*, '(A, A)') "Warning: Skipping malformed pair: ", TRIM(current_pair)
        END IF
        pair_index = pair_index - 1 ! Don't increment the count
      END IF
      
    END DO
    
    num_pairs_found = pair_index

  END SUBROUTINE parse_var_val_string

!------------------------------------------------------------------------------------------


  LOGICAL FUNCTION are_strings_equal_icase(str1, str2)
  !
  !Description
  !> @brief Compares two strings, ignoring case.
  !> @param[in] str1 The first input string.
  !> @param[in] str2 The second input string.
  !> @return .TRUE. if the strings are equal (ignoring case), .FALSE. otherwise.
    IMPLICIT NONE
    
    ! --- Arguments ---
    CHARACTER(LEN=*), INTENT(IN) :: str1
    CHARACTER(LEN=*), INTENT(IN) :: str2
    
    ! --- Local Variables ---
    INTEGER :: len1, len2
    INTEGER :: i
    CHARACTER(LEN=1) :: c1, c2
    
    ! --- Implementation ---
    
    ! Get the significant length of each string, ignoring trailing spaces
    len1 = LEN_TRIM(str1)
    len2 = LEN_TRIM(str2)
    
    ! If the lengths are different, they can't be equal.
    IF (len1 /= len2) THEN
      are_strings_equal_icase = .FALSE.
      RETURN
    END IF
    
    ! If lengths are the same, compare character by character
    DO i = 1, len1
      c1 = to_lower_char(str1(i:i))
      c2 = to_lower_char(str2(i:i))
      
      IF (c1 /= c2) THEN
        are_strings_equal_icase = .FALSE.
        RETURN
      END IF
    END DO
    
    ! All characters matched
    are_strings_equal_icase = .TRUE.
    
  END FUNCTION are_strings_equal_icase

!------------------------------------------------------------------------------------------


  CHARACTER(LEN=1) FUNCTION to_lower_char(c)
  !
  !Description
  !> @brief Converts a single character to lowercase.
  !> @param[in] c The character to convert.
  !> @return The lowercase version of the character.
  IMPLICIT NONE

  ! --- Arguments ---
  CHARACTER(LEN=1), INTENT(IN) :: c

  ! --- Local Variables ---
  INTEGER :: char_val
  INTEGER, PARAMETER :: offset = IACHAR('a') - IACHAR('A')

  ! --- Implementation ---
  char_val = IACHAR(c)

  ! Check if the character is an uppercase letter (A-Z)
  IF (char_val >= IACHAR('A') .AND. char_val <= IACHAR('Z')) THEN
    ! Convert to lowercase by adding the ASCII offset
    to_lower_char = ACHAR(char_val + offset)
  ELSE
    ! Return the original character
    to_lower_char = c
  END IF
    
  END FUNCTION to_lower_char
!------------------------------------------------------------------------------------------

  
  FUNCTION extract_number_and_unit(input_string, number_out, unit_out)result(success_out)
  !
  !Description
  !> @brief Extracts a number and a multi-character unit from a string.
  !> @details
  !>   Parses strings like "4K", "23.5cm", "100Pa", etc.
  !>   The string must start with a number and end with a unit.
  !>   The unit is assumed to start at the first alphabetic character.
  !>
  !> @param[in] input_string The string to be parsed.
  !> @param[out] number_out The extracted REAL number.
  !> @param[out] unit_out The extracted multi-char unit.
  !> @param[out] success_out .TRUE. if parsing was successful, .FALSE. otherwise.

    IMPLICIT NONE
    
    ! --- Arguments ---
    CHARACTER(LEN=*), INTENT(IN) :: input_string
    REAL, INTENT(OUT) :: number_out
    CHARACTER(LEN=*), INTENT(OUT) :: unit_out ! Changed to multi-char
    LOGICAL :: success_out
    
    ! --- Local Variables ---
    CHARACTER(LEN=LEN(input_string)) :: trimmed_str
    CHARACTER(LEN=LEN(input_string)) :: number_str ! Must be large enough
    INTEGER :: L
    INTEGER :: i
    INTEGER :: first_letter_pos
    INTEGER :: unit_val
    LOGICAL :: is_letter
    INTEGER :: read_status
    
    ! --- Implementation ---
    
    ! Initialize outputs to a failed state
    success_out = .FALSE.
    number_out = 0.0
    unit_out = ' '
    
    trimmed_str = TRIM(input_string)
    L = LEN_TRIM(trimmed_str)
    
    ! Find the position of the *first* character that is a letter
    first_letter_pos = 0
    DO i = 1, L
      unit_val = IACHAR(trimmed_str(i:i))
      is_letter = (unit_val >= IACHAR('A') .AND. unit_val <= IACHAR('Z')) .OR. &
                  (unit_val >= IACHAR('a') .AND. unit_val <= IACHAR('z')) .or. unit_val==iachar('%')
      IF (is_letter) THEN
        first_letter_pos = i
        EXIT
      END IF
    END DO
    
    ! Check for valid parse:
    ! 1. A letter must have been found (first_letter_pos > 0)
    ! 2. The letter must not be the very first character (first_letter_pos > 1)
    !    (This ensures there is a number part, e.g. "10cm" not "cm")
    IF (first_letter_pos > 1) THEN
      
      ! Everything before the first letter is the number
      number_str = trimmed_str(1 : first_letter_pos - 1)
      
      ! Everything from the first letter to the end is the unit
      unit_out = trimmed_str(first_letter_pos : L)
      
      ! Try to read the number part from the string
      ! IOSTAT is a non-zero value if the read fails
      READ(number_str, *, IOSTAT=read_status) number_out
      
      IF (read_status == 0) THEN
        ! The read was successful
        success_out = .TRUE.
      ELSE
        ! Read failed (e.g. "A 5 G" -> number_str = "A 5 ")
        ! Reset outputs to be safe
        success_out = .FALSE.
        number_out = 0.0
        unit_out = ' '
      END IF
      
    ELSE
      ! No letter found (e.g. "100") or
      ! letter is at the start (e.g. "K" or "cm")
      ! This is a parse failure according to the requirements.
      success_out = .FALSE.
      RETURN
    END IF
    
  END function extract_number_and_unit
!------------------------------------------------------------------------------------------
  FUNCTION to_lower_string(input_str)
    IMPLICIT NONE
    
    ! --- Arguments ---
    CHARACTER(LEN=*), INTENT(IN) :: input_str
    
    ! 2. Now define the function return type using the input's length
    CHARACTER(LEN=LEN(input_str)) :: to_lower_string    
    ! --- Local Variables ---
    INTEGER :: i
    
    ! --- Implementation ---
    to_lower_string = input_str ! Start by copying the string
    
    DO i = 1, LEN_TRIM(input_str)
      to_lower_string(i:i) = to_lower_char(input_str(i:i))
    END DO
    
  END FUNCTION to_lower_string
!------------------------------------------------------------------------------------------
  LOGICAL FUNCTION is_substring_present(x, y, ignore_case)
    IMPLICIT NONE
    
    ! --- Arguments ---
    CHARACTER(LEN=*), INTENT(IN) :: x
    CHARACTER(LEN=*), INTENT(IN) :: y
    LOGICAL, INTENT(IN), OPTIONAL :: ignore_case
    
    ! --- Local Variables ---
    INTEGER :: pos
    CHARACTER(LEN=LEN(x)) :: search_x
    CHARACTER(LEN=LEN(y)) :: search_y
    LOGICAL :: case_sensitive_search
    
    ! --- Implementation ---
    
    ! Determine if case-insensitive search is requested
    IF (PRESENT(ignore_case)) THEN
      case_sensitive_search = .NOT. ignore_case
    ELSE
      case_sensitive_search = .TRUE. ! Default to case-sensitive
    END IF
    
    IF (case_sensitive_search) THEN
      ! Perform a case-sensitive search using INDEX
      pos = INDEX(y, x)
    ELSE
      ! Convert both strings to lowercase for case-insensitive comparison
      ! (Need a helper function for this)
      search_x = to_lower_string(x)
      search_y = to_lower_string(y)
      pos = INDEX(search_y, search_x)
    END IF
    
    ! INDEX returns 0 if not found, non-zero if found.
    is_substring_present = (pos /= 0)
    
  END FUNCTION is_substring_present

!------------------------------------------------------------------------------------------

  FUNCTION string_compaction(input_string) RESULT(compacted_string)
  !> Normalize commas, NUL padding and line endings to spaces before field parsing.
  !> Preserve string length and other characters; repeated spaces are handled by
  !> count_delimited_items or list-directed input.
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: input_string
    CHARACTER(LEN=LEN(input_string)) :: compacted_string
    INTEGER :: i

    compacted_string = input_string
    DO i = 1, LEN_TRIM(compacted_string)
      SELECT CASE(compacted_string(i:i))
      CASE(',', ACHAR(0), ACHAR(10), ACHAR(13))
        compacted_string(i:i) = ' '
      END SELECT
    END DO
  END FUNCTION string_compaction

!------------------------------------------------------------------------------------------

  INTEGER FUNCTION count_delimited_items(input_string, delimiter)
  !
  !Description
  !> @brief Counts the number of items separated by a specified delimiter.
  !> @details
  !>   Parses a delimited string (e.g., "item1,item2,item3" with delimiter ',')
  !>   and counts the number of items. An empty string returns 0.
  !>   A space or tab delimiter counts whitespace-separated fields, ignoring
  !>   leading/trailing whitespace and treating runs of spaces/tabs as one separator.
  !>   Other consecutive delimiters are treated as separating empty items.
  !>
  !> @param[in] input_string The string to be parsed.
  !> @param[in] delimiter The delimiter character (e.g., ',', ';', ' ').
  !> @return The number of items found, or 0 if input_string is empty.

    IMPLICIT NONE

    ! --- Arguments ---
    CHARACTER(LEN=*), INTENT(IN) :: input_string
    CHARACTER(LEN=1), INTENT(IN) :: delimiter

    ! --- Local Variables ---
    INTEGER :: i
    INTEGER :: item_count
    INTEGER :: string_len
    LOGICAL :: in_item

    ! --- Implementation ---

    ! Initialize count
    item_count = 0

    ! Get trimmed length of input string
    string_len = LEN_TRIM(input_string)

    IF (delimiter == ' ' .OR. delimiter == ACHAR(9)) THEN
      in_item = .FALSE.
      DO i = 1, string_len
        IF (input_string(i:i) == ' ' .OR. input_string(i:i) == ACHAR(9)) THEN
          in_item = .FALSE.
        ELSE IF (.NOT. in_item) THEN
          item_count = item_count + 1
          in_item = .TRUE.
        END IF
      END DO
      count_delimited_items = item_count
      RETURN
    END IF

    ! If string is empty, return 0
    IF (string_len == 0) THEN
      count_delimited_items = 0
      RETURN
    END IF

    ! Count items: at least 1 item in non-empty string, plus 1 for each delimiter
    item_count = 1

    DO i = 1, string_len
      IF (input_string(i:i) == delimiter) THEN
        item_count = item_count + 1
      END IF
    END DO

    count_delimited_items = item_count

  END FUNCTION count_delimited_items

end module StringToolsMod
