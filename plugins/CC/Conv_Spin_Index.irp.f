integer function Conv_Spin_Index(p)

  implicit none

  integer,intent(in):: p

  Conv_Spin_Index = floor((p+1)/2)
end function Conv_Spin_Index
