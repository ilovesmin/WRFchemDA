subroutine da_transform_xtoy_chem_sfc (grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   type (domain),  intent(inout) :: grid
   type (iv_type), intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa) (linear)

   integer :: n,ichem       ! Loop counter.

   real, allocatable :: model_chemic(:,:)
   real, allocatable :: model_chemic_surf(:,:)
   real, allocatable :: model_rho(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_chem_sfc")

   if (sfc_assi_options == sfc_assi_options_1) then
      allocate (model_rho        (iv%info(chemic_surf)%n1:iv%info(chemic_surf)%n2,num_chem))
      allocate (model_chemic     (iv%info(chemic_surf)%n1:iv%info(chemic_surf)%n2,num_chem))
      allocate (model_chemic_surf(iv%info(chemic_surf)%n1:iv%info(chemic_surf)%n2,num_chemic_surf))

      ! Initialize first
      model_rho = 0.0
      model_chemic = 0.0
      model_chemic_surf = 0.0

      ! Interpolate horizontally:
      do ichem = PARAM_FIRST_SCALAR ,num_chem
         call da_interp_lin_2d (grid%xachem%chem_ic(:,:,1,ichem),  iv%info(chemic_surf), 1, model_chemic(:,ichem))
         call da_interp_lin_2d (grid%xb%rho(:,:,1),                iv%info(chemic_surf), 1, model_rho(:,ichem))
      end do

      if (chem_cv_options == 108) then      ! racm_soa_vbs_da

       if ( chemicda_opt == 1 ) then        ! pm2.5 only

           do n = p_chem_ic_so4aj, p_chem_ic_p25i
              model_chemic_surf(:,p_chemsi_pm25)=model_chemic_surf(:,p_chemsi_pm25)+model_rho(:,n)*model_chemic(:,n)
           enddo

       else if ( chemicda_opt == 2 ) then   ! pm10 only

           do n = p_chem_ic_so4aj, p_chem_ic_soila
              model_chemic_surf(:,p_chemsi_pm10)=model_chemic_surf(:,p_chemsi_pm10)+model_rho(:,n)*model_chemic(:,n)
           enddo

       else if ( chemicda_opt == 3 .or. chemicda_opt == 5) then

           do n = p_chem_ic_so4aj, p_chem_ic_p25i
              model_chemic_surf(:,p_chemsi_pm25)=model_chemic_surf(:,p_chemsi_pm25)+model_rho(:,n)*model_chemic(:,n)
           enddo

           do n = p_chem_ic_so4aj, p_chem_ic_soila
              model_chemic_surf(:,p_chemsi_pm10)=model_chemic_surf(:,p_chemsi_pm10)+model_rho(:,n)*model_chemic(:,n)
           enddo

       endif

       if ( chemicda_opt == 4 .or. chemicda_opt == 5 ) then

          model_chemic_surf(:,p_chemsi_so2)=model_chemic(:,p_chem_ic_so2)
          model_chemic_surf(:,p_chemsi_no2)=model_chemic(:,p_chem_ic_no2)
          model_chemic_surf(:,p_chemsi_o3) =model_chemic(:,p_chem_ic_o3 )
          model_chemic_surf(:,p_chemsi_co) =model_chemic(:,p_chem_ic_co )

       endif

      endif !(chem_cv_options == 108)

      do ichem = PARAM_FIRST_SCALAR, num_chemic_surf 
          do n=iv%info(chemic_surf)%n1,iv%info(chemic_surf)%n2
             if(iv%chemic_surf(n)%chem(ichem)%qc >= obs_qc_pointer) then
                 y%chemic_surf(n)%chem(ichem) = model_chemic_surf(n,ichem)
             else
                 y%chemic_surf(n)%chem(ichem) = 0.0
             endif
          end do
      end do

      deallocate (model_rho)
      deallocate (model_chemic)
      deallocate (model_chemic_surf)

   end if !(sfc_assi_options == sfc_assi_options_1) then

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_chem_sfc")

end subroutine da_transform_xtoy_chem_sfc
