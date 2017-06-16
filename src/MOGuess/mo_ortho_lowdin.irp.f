BEGIN_PROVIDER [double precision, ao_ortho_lowdin_coef, (ao_num_align,ao_num)]
  implicit none
  BEGIN_DOC
! matrix of the coefficients of the mos generated by the 
! orthonormalization by the S^{-1/2} canonical transformation of the aos
! ao_ortho_lowdin_coef(i,j) = coefficient of the ith ao on the jth ao_ortho_lowdin orbital
  END_DOC
  integer                        :: i,j,k,l
  double precision               :: accu
  double precision, allocatable  :: tmp_matrix(:,:)
  allocate (tmp_matrix(ao_num_align,ao_num))
  tmp_matrix(:,:) = 0.d0
  do j=1, ao_num
    tmp_matrix(j,j) = 1.d0
  enddo
  call ortho_lowdin(ao_overlap,ao_num_align,ao_num,tmp_matrix,ao_num_align,ao_num)
  do i=1, ao_num
    do j=1, ao_num
      ao_ortho_lowdin_coef(j,i) = tmp_matrix(i,j)
    enddo
  enddo
  deallocate(tmp_matrix)
END_PROVIDER

BEGIN_PROVIDER [double precision, ao_ortho_lowdin_overlap, (ao_num_align,ao_num)]
  implicit none
  BEGIN_DOC
! overlap matrix of the ao_ortho_lowdin
! supposed to be the Identity
  END_DOC
  integer                        :: i,j,k,l
  double precision               :: c
  do j=1, ao_num
    do i=1, ao_num
      ao_ortho_lowdin_overlap(i,j) = 0.d0
    enddo
  enddo
  do k=1, ao_num
    do j=1, ao_num
      c = 0.d0
      do l=1, ao_num
        c +=  ao_ortho_lowdin_coef(j,l) * ao_overlap(k,l)
      enddo
      do i=1, ao_num
        ao_ortho_lowdin_overlap(i,j) += ao_ortho_lowdin_coef(i,k) * c
      enddo
    enddo
  enddo
END_PROVIDER
