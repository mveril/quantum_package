subroutine getGe(Ge)
  BEGIN_DOC
! Bielectronic energy 
  END_DOC
  implicit none
  double precision::Je,Ke
  double precision,intent(out)::Ge
  call getJe(Je)
  call getKe(Ke)
  Ge=2d0*Je-Ke
  print *, '******'
  print *, 'Bielectronoic energy', Ge
  print *, '******'
end
