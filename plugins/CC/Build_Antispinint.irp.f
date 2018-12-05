subroutine Build_ijkl_antispinint(anti)
  
  BEGIN_DOC
  ! ijkl integrals
  END_DOC
  
  implicit none
  double precision,intent(out) :: anti(n_spin_occ,n_spin_occ,n_spin_occ,n_spin_occ)
  double precision :: spinint
  integer :: i,j,k,l
  anti(:,:,:,:) = 0d0

  do i=1, n_spin_occ
    do j=1,n_spin_occ
      do k=1,n_spin_occ
        do l=1,n_spin_occ
           anti(i,j,k,l) = spinint(i,j,k,l)-spinint(i,j,l,k)
        enddo
      enddo 
    enddo
  enddo
end subroutine

subroutine Build_ijab_antispinint(anti)

  BEGIN_DOC
  ! ijab integrals
  END_DOC

  implicit none

  DOUBLE PRECISION,intent(out) :: anti(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  double precision :: ijab_spinint
  integer :: i,j,a,b
  anti(:,:,:,:) = 0d0

  do i=1, n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          anti(i,j,a,b) = ijab_spinint(i,j,a,b)-ijab_spinint(i,j,b,a)
        enddo
      enddo 
    enddo
  enddo
end subroutine

subroutine Build_abcd_antispinint(anti)
  BEGIN_DOC
  ! abcd integrals
  END_DOC

  implicit none

  double precision,intent(out) :: anti(n_spin_virt,n_spin_virt,n_spin_virt,n_spin_virt)
  double precision :: abcd_spinint
  integer :: a,b,c,d
  anti(:,:,:,:) = 0d0

  do a=1, N_SPIN_Virt
    do b=1,n_spin_virt
      do c=1,n_spin_virt
        do d=1,n_spin_virt
          anti(a,b,c,d) = abcd_spinint(a,b,c,d)-abcd_spinint(a,b,d,c)
        enddo
      enddo 
    enddo
  enddo
end subroutine

subroutine Build_iajb_antispinint(anti)
  BEGIN_DOC
  ! iajb integrals
  END_DOC

  implicit none
  
  double precision,intent(out) :: anti(n_spin_occ,n_spin_virt,n_spin_occ,n_spin_virt)
  double precision :: iajb_spinint,iabj_spinint
  integer :: i,a,j,b
  anti(:,:,:,:) = 0d0
 
  do i=1, n_spin_occ
    do a=1,n_spin_virt
      do j=1,n_spin_occ
        do b=1,n_spin_virt
          anti(i,a,j,b) = iajb_spinint(i,a,j,b)-iabj_spinint(i,a,b,j)
        enddo
      enddo 
    enddo
  enddo
end subroutine


subroutine Build_abij_antispinint(anti)
  BEGIN_DOC
  ! abij integrals
  END_DOC
  implicit none
  
  double precision,intent(out)::anti(n_spin_virt,n_spin_virt,n_spin_occ,n_spin_occ)
  double precision::abij_spinint
  integer :: a,b,i,j
  anti(:,:,:,:) = 0d0

  do a=1,n_spin_virt
    do b=1,n_spin_virt
      do i=1, n_spin_occ
        do j=1,n_spin_occ
          anti(a,b,i,j) = abij_spinint(a,b,i,j)-abij_spinint(a,b,j,i)
        enddo
      enddo 
    enddo
  enddo
end subroutine
