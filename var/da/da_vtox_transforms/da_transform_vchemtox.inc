subroutine da_transform_vchemtox(grid, vchem, be)

   !-----------------------------------------------------------------------
   ! Purpose: Transform chemical CV to model space
   !-----------------------------------------------------------------------

   implicit none

   type(domain),  intent(inout) :: grid
   type(xchem_type), intent(in) :: vchem   ! Grid point/EOF equivalent.
   type(be_type), intent(in)    :: be   ! Background error structure.

   integer :: ij, ic, i, j, k, m  ! Loop counters
   integer :: mz
   real :: temp

   if (trace_use) call da_trace_entry("da_transform_vchemtox")

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, ic, k, m, j, i, temp )

   do ij = 1 , grid%num_tiles

      grid%xachem%chem_ic(:,grid%j_start(ij):grid%j_end(ij),:,:) = 0.0

      do ic = PARAM_FIRST_SCALAR, num_chem
         mz = be % v12(ic-1) % mz     !EIGENMODES

         do k = kts, kte  !levels
            do m = 1, mz
               do j = grid%j_start(ij), grid%j_end(ij)
                  temp = be % v12(ic-1) % evec (j,k,m) * be % v12(ic-1) % val (j,m)

                  do i = its, ite
                      grid%xachem%chem_ic (i,j,k,ic) = grid%xachem%chem_ic (i,j,k,ic) + &
                                                       temp * vchem%chem_ic(i,j,m,ic)
                  end do

               end do
            end do
         end do
      end do
   end do
   !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
#include "HALO_CHEM_XA.inc"     
#endif

   if (trace_use) call da_trace_exit("da_transform_vchemtox")

end subroutine da_transform_vchemtox
