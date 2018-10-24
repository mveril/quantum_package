program CC
  implicit none
  double precision::HFe
  BEGIN_DOC
! Coupled Cluster module
  END_DOC
  call getHFe(HFe)
  print *, '******'
  print *, "ΔHFe", HFE-hf_energy
  print *, '******'
end
