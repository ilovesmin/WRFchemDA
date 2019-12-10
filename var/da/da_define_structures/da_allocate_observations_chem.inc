subroutine da_allocate_observations_chem (iv)

   !---------------------------------------------------------------------------
   ! Purpose: Allocate components of chem observation structure.
   !---------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(inout) :: iv     ! Observation structure.

   integer :: i

   if (trace_use) call da_trace_entry("da_allocate_observations_chem")

   if (iv%info(chem_surf)%nlocal    > 0) allocate(iv%chem_surf   (1:iv%info(chem_surf)%nlocal))
   if (iv%info(chem_acft)%nlocal    > 0) allocate(iv%chem_acft   (1:iv%info(chem_acft)%nlocal))

   i = chem_surf
      if (iv%info(i)%nlocal > 0) then
         allocate (iv%info(i)%name(iv%info(i)%nlocal))
         allocate (iv%info(i)%platform(iv%info(i)%nlocal))
         allocate (iv%info(i)%id(iv%info(i)%nlocal))
         allocate (iv%info(i)%date_char(iv%info(i)%nlocal))
         allocate (iv%info(i)%levels(iv%info(i)%nlocal))
         allocate (iv%info(i)%lat(iv%info(i)%max_lev,iv%info(i)%nlocal))
         allocate (iv%info(i)%lon(iv%info(i)%max_lev,iv%info(i)%nlocal))
         allocate (iv%info(i)%elv(iv%info(i)%nlocal))
         allocate (iv%info(i)%pstar(iv%info(i)%nlocal))

         allocate (iv%info(i)%slp(iv%info(i)%nlocal))
         allocate (iv%info(i)%pw(iv%info(i)%nlocal))

         allocate (iv%info(i)%x  (kms:kme,iv%info(i)%nlocal))
         allocate (iv%info(i)%y  (kms:kme,iv%info(i)%nlocal))
         allocate (iv%info(i)%i  (kms:kme,iv%info(i)%nlocal))
         allocate (iv%info(i)%j  (kms:kme,iv%info(i)%nlocal))
         allocate (iv%info(i)%dx (kms:kme,iv%info(i)%nlocal))
         allocate (iv%info(i)%dxm(kms:kme,iv%info(i)%nlocal))
         allocate (iv%info(i)%dy (kms:kme,iv%info(i)%nlocal))
         allocate (iv%info(i)%dym(kms:kme,iv%info(i)%nlocal))
         allocate (iv%info(i)%k  (iv%info(i)%max_lev,iv%info(i)%nlocal))
         allocate (iv%info(i)%dz (iv%info(i)%max_lev,iv%info(i)%nlocal))
         allocate (iv%info(i)%dzm(iv%info(i)%max_lev,iv%info(i)%nlocal))
         allocate (iv%info(i)%zk (iv%info(i)%max_lev,iv%info(i)%nlocal))
         allocate (iv%info(i)%proc_domain(iv%info(i)%max_lev,iv%info(i)%nlocal))
         allocate (iv%info(i)%thinned(iv%info(i)%max_lev,iv%info(i)%nlocal))
         allocate (iv%info(i)%obs_global_index(iv%info(i)%nlocal))

         iv%info(i)%proc_domain(:,:)  = .false.
         iv%info(i)%thinned(:,:)      = .false.
         iv%info(i)%zk(:,:)           = missing_r
      end if

   i = chem_acft
      if (iv%info(i)%nlocal > 0) then
         allocate (iv%info(i)%name(iv%info(i)%nlocal))
         allocate (iv%info(i)%platform(iv%info(i)%nlocal))
         allocate (iv%info(i)%id(iv%info(i)%nlocal))
         allocate (iv%info(i)%date_char(iv%info(i)%nlocal))
         allocate (iv%info(i)%levels(iv%info(i)%nlocal))
         allocate (iv%info(i)%lat(1:num_platform,iv%info(i)%nlocal))
         allocate (iv%info(i)%lon(1:num_platform,iv%info(i)%nlocal))
         allocate (iv%info(i)%elv(iv%info(i)%nlocal))
         allocate (iv%info(i)%pstar(iv%info(i)%nlocal))

         allocate (iv%info(i)%slp(iv%info(i)%nlocal))
         allocate (iv%info(i)%pw(iv%info(i)%nlocal))

         allocate (iv%info(i)%x  (1:num_platform,iv%info(i)%nlocal))
         allocate (iv%info(i)%y  (1:num_platform,iv%info(i)%nlocal))
         allocate (iv%info(i)%i  (1:num_platform,iv%info(i)%nlocal))
         allocate (iv%info(i)%j  (1:num_platform,iv%info(i)%nlocal))
         allocate (iv%info(i)%dx (1:num_platform,iv%info(i)%nlocal))
         allocate (iv%info(i)%dxm(1:num_platform,iv%info(i)%nlocal))
         allocate (iv%info(i)%dy (1:num_platform,iv%info(i)%nlocal))
         allocate (iv%info(i)%dym(1:num_platform,iv%info(i)%nlocal))
         allocate (iv%info(i)%k  (1:num_platform,iv%info(i)%nlocal))
         allocate (iv%info(i)%dz (1:num_platform,iv%info(i)%nlocal))
         allocate (iv%info(i)%dzm(1:num_platform,iv%info(i)%nlocal))
         allocate (iv%info(i)%zk (1:num_platform,iv%info(i)%nlocal))
         allocate (iv%info(i)%proc_domain(1:num_platform,iv%info(i)%nlocal))
         allocate (iv%info(i)%thinned(1:num_platform,iv%info(i)%nlocal))
         allocate (iv%info(i)%obs_global_index(iv%info(i)%nlocal))

         iv%info(i)%proc_domain(:,:)  = .false.
         iv%info(i)%thinned(:,:)      = .false.
         iv%info(i)%zk(:,:)           = missing_r
      end if

   if (trace_use) call da_trace_exit("da_allocate_observations_chem")

end subroutine da_allocate_observations_chem

