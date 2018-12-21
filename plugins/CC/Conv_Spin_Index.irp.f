integer function Conv_Spin_Index(p)

  implicit none

  BEGIN_DOC
!  convert Spinorbital index to spacial orbital integral index
  END_DOC

  integer,intent(in):: p

! Divide integer by integer return the result of the eclidean division as integer

  Conv_Spin_Index = (p+1)/2 
end function Conv_Spin_Index
