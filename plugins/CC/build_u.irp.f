subroutine build_u(u,oooo_db_spin_int,ovov_db_spin_int,abcd_antispinint,t2)
  
  BEGIN_DOC
  ! U matrix (defined on equation 13 of pople paper 
  END_DOC
  
  implicit none
  
  
  double precision,intent(out) :: u(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  
  double precision,intent(in) :: oooo_db_spin_int(n_spin_occ,n_spin_occ,n_spin_occ,n_spin_occ)
  double precision,intent(in) :: ovov_db_spin_int(n_spin_occ,n_spin_virt,n_spin_occ,n_spin_virt)
  double precision,intent(in) :: abcd_antispinint(n_spin_virt,n_spin_virt,n_spin_virt,n_spin_virt)
  double precision,intent(in) :: t2(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  
  double precision :: pcd, pkl, pkc 
  integer :: i,j,a,b
  integer :: k,l,c,d

  u(:,:,:,:) = 0d0

  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt

          pcd = 0d0
          do c=1,n_spin_virt
            do d=1,n_spin_virt
              pcd += abcd_antispinint(a,b,c,d)*t2(i,j,c,d)
            end do
          end do

          pkl = 0d0
          pkc = 0d0
          do k=1,n_spin_occ

            do l=1 ,n_spin_occ
              pkl += oooo_db_spin_int(k,l,i,j)*t2(k,l,a,b)
            end do

            do c=1,n_spin_virt
              pkc = pkc                                   &
                  - ovov_db_spin_int(k,b,j,c)*t2(i,k,a,c) &
                  + ovov_db_spin_int(k,a,j,c)*t2(i,k,b,c) &
                  - ovov_db_spin_int(k,a,i,c)*t2(j,k,b,c) &
                  + ovov_db_spin_int(k,b,i,c)*t2(j,k,a,c)
            end do

          end do

          u(i,j,a,b) = 0.5d0*pcd + 0.5d0*pkl + pkc

        end do
      end do
    end do
  end do

end subroutine build_u
