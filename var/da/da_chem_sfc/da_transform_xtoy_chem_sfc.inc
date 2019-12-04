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
   real :: psum, bsum, pm2

   real, allocatable :: model_rho(:,:)          ! model rho for each species at obs site [kg/m^3]
   real, allocatable :: model_chemic(:,:)       ! model incr. in aerosol/gas species at obs site
   real, allocatable :: model_chemic_surf(:,:)  ! adjustment in aerosol/gas species at obs site

   ! Needed for met and chem DA at the same time as rho will be changed by p and t increments.
!   real, allocatable :: model_drho(:,:)         ! model incr. in rho for each species at obs site [kg/m^3]
!   real, allocatable :: model_chemic_b(:,:)      ! model aerosol/gas species at obs site (background)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_chem_sfc")

   if (sfc_assi_options == sfc_assi_options_1) then
      allocate (model_rho        (iv%info(chemic_surf)%n1:iv%info(chemic_surf)%n2,num_chem))
      allocate (model_chemic     (iv%info(chemic_surf)%n1:iv%info(chemic_surf)%n2,num_chem))
!      allocate (model_chemic_b   (iv%info(chemic_surf)%n1:iv%info(chemic_surf)%n2,num_chem))
      allocate (model_chemic_surf(iv%info(chemic_surf)%n1:iv%info(chemic_surf)%n2,num_chemic_surf))

      ! Interpolate horizontally:
      do ichem = PARAM_FIRST_SCALAR ,num_chem
!         call da_interp_lin_2d (grid%xbchem%chem_ic(:,:,1,ichem), iv%info(chemic_surf), 1, model_chemic_b(:,ichem))
         call da_interp_lin_2d (grid%xachem%chem_ic(:,:,1,ichem), iv%info(chemic_surf), 1, model_chemic(:,ichem))
         call da_interp_lin_2d (grid%xb%rho(:,:,1),               iv%info(chemic_surf), 1, model_rho   (:,ichem))
      end do

      !do n = 1, 10
      !psum = sum(model_chemic   (n,p_chem_ic_so4aj:p_chem_ic_p25i))   ! xachem
      !bsum = sum(model_chemic_b (n,p_chem_ic_so4aj:p_chem_ic_p25i))   ! xbchem
      ! pm2 = model_rho(n,6)*psum
      !print '(A,I2,7F9.3)',"da_transform_xtoy_chem_sfc: n,sum(xachem),rho,pm2,sum(xbchem),pm2_b=",&
      !       n,psum,model_rho(n,6),pm2,bsum,model_rho(n,6)*bsum
      !end do

      if (chem_cv_options == 108) then      ! racm_soa_vbs_da

         if (chemicda_opt == 1 .or. chemicda_opt == 3 .or. chemicda_opt == 4) then    ! pm2.5
             model_chemic_surf(:,p_chemsi_pm25) = 0.0
             do ichem = p_chem_ic_so4aj, p_chem_ic_p25i
                model_chemic_surf(:,p_chemsi_pm25) = model_chemic_surf(:,p_chemsi_pm25) + &
                                                     model_rho(:,ichem)  * model_chemic(:,ichem)
             enddo
         endif  !(chemicda_opt == 1 .or. chemicda_opt == 3 .or. chemicda_opt == 4) 

         if (chemicda_opt == 2) then        ! pm10 only
             model_chemic_surf(:,p_chemsi_pm10) = 0.0
             do ichem = p_chem_ic_so4aj, p_chem_ic_soila
                model_chemic_surf(:,p_chemsi_pm10) = model_chemic_surf(:,p_chemsi_pm10) + &
                                                     model_rho(:,ichem) * model_chemic(:,ichem)
             enddo
         endif !(chemicda_opt == 2)

         if (chemicda_opt == 3 .or. chemicda_opt == 4) then    ! pm10 after pm2.5
             ! pm10 - pm2.5 residual
             model_chemic_surf(:,p_chemsi_pm10) = 0.0
             do ichem = p_chem_ic_antha, p_chem_ic_soila
                model_chemic_surf(:,p_chemsi_pm10) = model_chemic_surf(:,p_chemsi_pm10) + &
                                                     model_rho(:,ichem) * model_chemic(:,ichem)
             enddo
         endif !(chemicda_opt == 3 .or. chemicda_opt == 4)

         if ( chemicda_opt == 4 .or. chemicda_opt == 5 ) then
              model_chemic_surf(:,p_chemsi_so2)=model_chemic(:,p_chem_ic_so2)
              model_chemic_surf(:,p_chemsi_no2)=model_chemic(:,p_chem_ic_no2)
              model_chemic_surf(:,p_chemsi_o3) =model_chemic(:,p_chem_ic_o3 )
              model_chemic_surf(:,p_chemsi_co) =model_chemic(:,p_chem_ic_co )
         endif

         do ichem = PARAM_FIRST_SCALAR, num_chemic_surf
          do n=iv%info(chemic_surf)%n1,iv%info(chemic_surf)%n2
             if(iv%chemic_surf(n)%chem(ichem)%qc >= obs_qc_pointer) then
                 y%chemic_surf(n)%chem(ichem) = model_chemic_surf(n,ichem)
             else
                 y%chemic_surf(n)%chem(ichem) = 0.0
             endif
             !print '(A,2I2,7F9.3)',"da_transform_xtoy_chem_sfc: n,ichemda,y=",n,ichem,y%chemic_surf(n)%chem(ichem)
          end do
         end do

      endif !(chem_cv_options == 108)

      deallocate (model_rho)
      deallocate (model_chemic)
      deallocate (model_chemic_surf)

   end if !(sfc_assi_options == sfc_assi_options_1) then

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_chem_sfc")

end subroutine da_transform_xtoy_chem_sfc
