program check_ortho
 implicit none

 call do_print
 print *,  '--'
 print *,  '--'
 call orthonormalize_mos
 call do_print
end

subroutine do_print
 implicit none
 integer :: i,j
 double precision :: off_diag, diag

 off_diag = 0.d0
 diag = 0.d0
 do j=1,mo_tot_num
   do i=1,mo_tot_num
     off_diag += abs(mo_overlap(i,j))
   enddo
   diag += abs(mo_overlap(j,j))
   off_diag -= abs(mo_overlap(j,j))
 enddo
 print *,  'Diag      = ', abs(1.d0-diag/dble(mo_tot_num))
 print *,  'Off-Diag  = ', off_diag/(dble(mo_tot_num)**2-dble(mo_tot_num))
 print *,  '--'

end


