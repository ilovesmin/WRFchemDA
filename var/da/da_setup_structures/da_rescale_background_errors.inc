subroutine da_rescale_background_errors (var_scaling, len_scaling, &
                                          ds, s, be_sub)

   !---------------------------------------------------------------------------
   ! Purpose: Rescale wrfvar background errors.
   !
   ! Method:  Empirical scaling and inclusion of recursive filter rescaling.
   !---------------------------------------------------------------------------

   implicit none

   real, intent(in)                 :: var_scaling       ! Variance factor.
   real, intent(in)                 :: len_scaling       ! Lengthscale factor.
   real, intent(in)                 :: ds                ! Resolution (m)
   real*8, intent(inout)            :: s(:)              ! RF lengthscale.
   type (be_subtype), intent(inout) :: be_sub            ! Backgrd error struct.
    
   integer                          :: mz                ! Vertical truncation.
   integer                          :: m
   real*8, allocatable              :: rf_scale_factor(:)! RF rescaling.

   if (trace_use_dull) call da_trace_entry("da_rescale_background_errors")

   write(unit=stdout,fmt='(3x,"Scaling: var, len, ds:",3f9.2)') &
                                     var_scaling, len_scaling, ds

   !--------------------------------------------------------------------------
   ! [1.0] Initialise:
   !--------------------------------------------------------------------------

   mz = be_sub % mz

   !--------------------------------------------------------------------------
   ! [2.0] Perform various rescalings:
   !--------------------------------------------------------------------------

   if (mz > 0) then

      ! [2.1] Empirical rescaling of lengthscales:
      s(1:mz) = len_scaling * s(1:mz)
   
      if (print_detail_be) then
         !write(unit=stdout,fmt='(a,a)')trim(be_sub % name), ' Error Lengthscales (m=1):'
         !do m = 1, mz
         !   write(unit=stdout,fmt='(a,i4,1pe13.5)')be_sub % name, m, s(m)
         !end do
      end if
      
      ! [2.2] Make lengthscale nondimensional:
      s(1:mz) = s(1:mz) / ds

      ! [2.3] Empirical rescaling of variances:
      be_sub % val(:,:) = var_scaling * be_sub % val(:,:)

      ! Calculate recursive filter rescaling:

      allocate(rf_scale_factor(1:mz))

      call da_calculate_rf_factors(s(:), be_sub % rf_alpha(:), &
                                    rf_scale_factor(:))

      do m = 1, mz
         be_sub % val(:,m) = rf_scale_factor(m) * be_sub % val(:,m)
      end do
      !m = 1
      !write(unit=stdout,fmt='(a,a,3f12.5,3f7.1)') 'rescale_be (j=89,m=1):',&
      !      be_sub%name, s(m),be_sub%val(89,m), var_scaling, len_scaling, ds
                                       
      deallocate (rf_scale_factor)   

   end if

   if (trace_use_dull) call da_trace_exit("da_rescale_background_errors")

end subroutine da_rescale_background_errors


