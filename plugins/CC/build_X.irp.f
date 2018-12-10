subroutine build_X1(X1,ijab_antispinint,t2)
  
  BEGIN_DOC
  ! X1 matrix 
  END_DOC
  
  implicit none
  
  
  double precision,intent(out) :: X1(n_spin_occ,n_spin_occ,n_spin_occ,n_spin_occ)
  double precision,intent(in),dimension(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt) :: ijab_antispinint,t2
  integer :: k,l,i,j
  integer :: c,d

  X1(:,:,:,:) = 0d0

  do k=1,n_spin_occ
    do l=1,n_spin_occ
      do i=1,n_spin_occ
        do j=1,n_spin_occ
          do c=1,n_spin_virt
            do d=1,n_spin_virt
              X1(k,l,i,j) += ijab_antispinint(k,l,c,d)*t2(i,j,c,d)
            end do
          end do
        end do
      end do
    end do
  end do

end subroutine build_X1

subroutine build_X2(X2,ijab_antispinint,t2)
  
  BEGIN_DOC
  ! X2 matrix
  END_DOC

  implicit none
  
  double precision,intent(out) :: X2(n_spin_virt,n_spin_virt)
  double precision,intent(in),dimension(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt) :: ijab_antispinint,t2
  integer :: b,c
  integer::k,l,d

  X2(:,:) = 0d0

  do b=1,n_spin_virt
    do c=1,n_spin_virt
      do k=1,n_spin_occ
        do l=1,n_spin_occ
          do d=1,n_spin_virt
            X2(b,c) += ijab_antispinint(k,l,c,d)*t2(k,l,b,d)
          enddo
        enddo
    enddo
  enddo
enddo

end subroutine build_X2

subroutine build_X3(X3,ijab_antispinint,t2)

  BEGIN_DOC
  ! X3 matrix
  END_DOC

  implicit none

  double precision,intent(out) :: X3(n_spin_occ,n_spin_occ)
  double precision,intent(in),dimension(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt) :: ijab_antispinint,t2
  integer :: k,j
  integer :: l,c,d
  
  X3(:,:) = 0d0

  do k=1,n_spin_occ
    do j=1,n_spin_occ
      do d=1,n_spin_virt 
        do c=1,n_spin_virt
          do l=1,n_spin_occ
            X3(k,j) += ijab_antispinint(k,l,c,d)*t2(j,l,c,d)
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine build_X3

subroutine build_X4(X4,ijab_antispinint,t2)
  
  BEGIN_DOC
  ! X4 matrix
  END_DOC

  implicit none
  
  double precision,intent(out) :: X4(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  double precision,intent(in),dimension(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt) :: ijab_antispinint,t2
  integer::i,l,a,d
  integer::k,c

  X4(:,:,:,:) = 0d0

  do i=1,n_spin_occ
    do l=1,n_spin_occ
      do a=1,n_spin_virt
        do d=1,n_spin_virt
          do k=1,n_spin_occ
            do c=1,n_spin_virt
              X4(i,l,a,d) += ijab_antispinint(k,l,c,d)*t2(i,k,a,c)
            end do
          end do
        end do
      end do
    end do
  end do

end subroutine build_X4
