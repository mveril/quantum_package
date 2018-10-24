subroutine getHFe(HFe)
  BEGIN_DOC
! Hartree Fock energy 
  END_DOC
  implicit none
  double precision::Ce,Ge
  double precision,intent(out)::HFe
  call getCe(Ce)
  call getGe(Ge)
  HFe=nuclear_repulsion+Ce+Ge
  print *, '******'
  print *, 'HFEnergy', HFe
  print *, '******'
end
