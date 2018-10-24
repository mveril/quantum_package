subroutine getCe(Ce)
  BEGIN_DOC
! Core energy 
  END_DOC
  implicit none
  double precision::Ve,Te
  double precision,intent(out)::Ce
  call getTe(Te)
  call getVe(Ve)
  Ce=Te+Ve
  print *, '******'
  print *, 'Core energy', Ce
  print *, '******'
end
