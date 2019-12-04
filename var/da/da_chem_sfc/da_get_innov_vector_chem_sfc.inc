subroutine da_get_innov_vector_chem_sfc( it, num_qcstat_conv, grid, ob, iv)
   !-----------------------------------------------------------------------
   ! Purpose: TBD    
   !-----------------------------------------------------------------------

   implicit none

   integer,          intent(in)    :: it      ! External iteration.
   type(domain),     intent(in)    :: grid    ! first guess state.
   type(y_type),     intent(inout) :: ob      ! Observation structure.
   type(iv_type),    intent(inout) :: iv      ! O-B structure.
   integer,          intent(inout) :: num_qcstat_conv(:,:,:,:)

   integer :: n,ichem        ! Loop counters
   integer :: imax

   real, allocatable :: model_chemic_surf(:,:)
   real, allocatable :: model_chemic(:,:)
   real, allocatable :: model_rho(:,:)

   real, parameter  :: scaleL = 3000.         ! scaling factor for obs error [meters]
   real, parameter  :: maxomb = 100.          ! maximum threshold for abs(o-b)
   real             :: err_o, err_r, maxv, maxe, maxb, maxo

   if (trace_use) call da_trace_entry("da_get_innov_vector_chem_sfc")
   print*,'da_get_innov_vector_chem_sfc: chem_cv_options, chemicda_opt =', chem_cv_options,chemicda_opt

   if ( it > 1 ) then
      do n=iv%info(chemic_surf)%n1,iv%info(chemic_surf)%n2
         do ichem = PARAM_FIRST_SCALAR, num_chemic_surf
            if (iv%chemic_surf(n)%chem(ichem)%qc == fails_error_max) iv%chemic_surf(n)%chem(ichem)%qc = 0
         end do
      end do
   end if

   ! [0.0]  Get chem surf conc.

   allocate (model_chemic_surf(iv%info(chemic_surf)%n1:iv%info(chemic_surf)%n2,num_chemic_surf))
   allocate (model_chemic     (iv%info(chemic_surf)%n1:iv%info(chemic_surf)%n2,num_chem))
   allocate (model_rho        (iv%info(chemic_surf)%n1:iv%info(chemic_surf)%n2,num_chem))

   do ichem = PARAM_FIRST_SCALAR ,num_chem
      ! [1.0] Get horizontal interpolation weights: 
      call da_interp_lin_2d (grid%xbchem%chem_ic(:,:,1,ichem),  iv%info(chemic_surf), 1, model_chemic(:,ichem))
      call da_interp_lin_2d (grid%xb    %rho    (:,:,1),        iv%info(chemic_surf), 1, model_rho   (:,ichem))
   end do

   if (chem_cv_options == 108) then  ! racm_soa_vbs_da

      if (chemicda_opt == 1 .or. chemicda_opt == 3 .or. chemicda_opt == 4) then

         ! pm2.5
         model_chemic_surf(:,p_chemsi_pm25)=model_rho(:,p_chem_ic_so4aj) * &
                                           (model_chemic(:,p_chem_ic_so4aj) +model_chemic(:,p_chem_ic_so4ai) + &
                                            model_chemic(:,p_chem_ic_nh4aj) +model_chemic(:,p_chem_ic_nh4ai) + &
                                            model_chemic(:,p_chem_ic_no3aj) +model_chemic(:,p_chem_ic_no3ai) + &
                                            model_chemic(:,p_chem_ic_naaj ) +model_chemic(:,p_chem_ic_naai ) + &
                                            model_chemic(:,p_chem_ic_claj ) +model_chemic(:,p_chem_ic_clai ) + &
                                            model_chemic(:,p_chem_ic_asoa1j)+model_chemic(:,p_chem_ic_asoa1i)+ &
                                            model_chemic(:,p_chem_ic_asoa2j)+model_chemic(:,p_chem_ic_asoa2i)+ &
                                            model_chemic(:,p_chem_ic_asoa3j)+model_chemic(:,p_chem_ic_asoa3i)+ &
                                            model_chemic(:,p_chem_ic_asoa4j)+model_chemic(:,p_chem_ic_asoa4i)+ &
                                            model_chemic(:,p_chem_ic_bsoa1j)+model_chemic(:,p_chem_ic_bsoa1i)+ &
                                            model_chemic(:,p_chem_ic_bsoa2j)+model_chemic(:,p_chem_ic_bsoa2i)+ &
                                            model_chemic(:,p_chem_ic_bsoa3j)+model_chemic(:,p_chem_ic_bsoa3i)+ &
                                            model_chemic(:,p_chem_ic_bsoa4j)+model_chemic(:,p_chem_ic_bsoa4i)+ &
                                            model_chemic(:,p_chem_ic_orgpaj)+model_chemic(:,p_chem_ic_orgpai)+ &
                                            model_chemic(:,p_chem_ic_ecj   )+model_chemic(:,p_chem_ic_eci   )+ &
                                            model_chemic(:,p_chem_ic_p25j  )+model_chemic(:,p_chem_ic_p25i  ))
      endif

      if (chemicda_opt == 2) then

         ! pm10
         model_chemic_surf(:,p_chemsi_pm10)=model_rho(:,p_chem_ic_so4aj) * &
                                           (model_chemic(:,p_chem_ic_so4aj) +model_chemic(:,p_chem_ic_so4ai) + &
                                            model_chemic(:,p_chem_ic_nh4aj) +model_chemic(:,p_chem_ic_nh4ai) + &
                                            model_chemic(:,p_chem_ic_no3aj) +model_chemic(:,p_chem_ic_no3ai) + &
                                            model_chemic(:,p_chem_ic_naaj ) +model_chemic(:,p_chem_ic_naai ) + &
                                            model_chemic(:,p_chem_ic_claj ) +model_chemic(:,p_chem_ic_clai ) + &
                                            model_chemic(:,p_chem_ic_asoa1j)+model_chemic(:,p_chem_ic_asoa1i)+ &
                                            model_chemic(:,p_chem_ic_asoa2j)+model_chemic(:,p_chem_ic_asoa2i)+ &
                                            model_chemic(:,p_chem_ic_asoa3j)+model_chemic(:,p_chem_ic_asoa3i)+ &
                                            model_chemic(:,p_chem_ic_asoa4j)+model_chemic(:,p_chem_ic_asoa4i)+ &
                                            model_chemic(:,p_chem_ic_bsoa1j)+model_chemic(:,p_chem_ic_bsoa1i)+ &
                                            model_chemic(:,p_chem_ic_bsoa2j)+model_chemic(:,p_chem_ic_bsoa2i)+ &
                                            model_chemic(:,p_chem_ic_bsoa3j)+model_chemic(:,p_chem_ic_bsoa3i)+ &
                                            model_chemic(:,p_chem_ic_bsoa4j)+model_chemic(:,p_chem_ic_bsoa4i)+ &
                                            model_chemic(:,p_chem_ic_orgpaj)+model_chemic(:,p_chem_ic_orgpai)+ &
                                            model_chemic(:,p_chem_ic_ecj   )+model_chemic(:,p_chem_ic_eci   )+ &
                                            model_chemic(:,p_chem_ic_p25j  )+model_chemic(:,p_chem_ic_p25i  )+ &
                                            model_chemic(:,p_chem_ic_antha )+model_chemic(:,p_chem_ic_seas ) + &
                                            model_chemic(:,p_chem_ic_soila ))

      endif !(chemicda_opt == 2) then

      if (chemicda_opt == 3 .or. chemicda_opt == 4) then   ! pm2.5 and pm10

         ! pm10 - pm2.5
         model_chemic_surf(:,p_chemsi_pm10)=model_rho(:,p_chem_ic_antha) * &
                                           (model_chemic(:,p_chem_ic_antha )+model_chemic(:,p_chem_ic_seas ) + &
                                            model_chemic(:,p_chem_ic_soila ))

      endif !(chemicda_opt == 3 .or. chemicda_opt == 4)

      if (chemicda_opt == 4 .or. chemicda_opt == 5) then

          model_chemic_surf(:,p_chemsi_so2)=model_chemic(:,p_chem_ic_so2)
          model_chemic_surf(:,p_chemsi_no2)=model_chemic(:,p_chem_ic_no2)
          model_chemic_surf(:,p_chemsi_o3) =model_chemic(:,p_chem_ic_o3 )
          model_chemic_surf(:,p_chemsi_co) =model_chemic(:,p_chem_ic_co )

      endif !(chemicda_opt == 4 ) then

   end if   !(chem_cv_options == 108) then  ! racm_soa_vbs_da

   do ichem = PARAM_FIRST_SCALAR, num_chemic_surf
 
      do n=iv%info(chemic_surf)%n1,iv%info(chemic_surf)%n2

         ! [2.0] Initialise components of innovation vector: 

         iv % chemic_surf(n) % chem(ichem) % inv   = 0.0
         iv % chemic_surf(n) % chem(ichem) % error = missing_r

         ! [3.0] To compute the 'inv':
         if (ob % chemic_surf(n) % chem(ichem) > missing_r .and. iv % chemic_surf(n) % chem(ichem) % qc >=  obs_qc_pointer) then

            ! compute innovation (= o - f)
            iv % chemic_surf(n) % chem(ichem) % inv = ob % chemic_surf(n) % chem(ichem) - model_chemic_surf(n,ichem)
            err_o = 1.5 + 0.0075 * ob % chemic_surf(n) % chem(ichem)
            err_r = 0.5 * err_o * sqrt(grid % dx / scaleL)
            iv % chemic_surf(n) % chem(ichem) % error = sqrt(err_o**2 + err_r**2)

            ! observation error assignment (into iv % error) - use the same formula for all species
  	    !  iv % chemic_surf(n) % chem(ichem) % error = sqrt((1.0+0.0075*ob % chemic_surf(n) % chem(ichem))**2+ &
            !      0.5*(1.0+0.0075*ob % chemic_surf(n) % chem(ichem))*sqrt(grid%dx/scaleL))

            ! Wei's error formula below
            ! iv % chemic_surf(n) % chem(ichem) % error = 1.5*sqrt((1.0+0.0075*ob % chemic_surf(n) % chem(ichem))**2+ &
            !      0.5*(1.0+0.0075*ob % chemic_surf(n) % chem(ichem))*sqrt(grid%dx/2000.0))

            ! Quality Control
            ! ----------------
            if (abs(iv % chemic_surf(n) % chem(ichem) % inv) .ge. maxomb) then
                iv % chemic_surf(n) % chem(ichem) % inv   = 0.0   !missing_r
                iv % chemic_surf(n) % chem(ichem) % qc    = missing
                iv % chemic_surf(n) % chem(ichem) % error = missing_r
            end if

         end if	   ! good obs      
      end do
   end do

   deallocate(model_rho)
   deallocate(model_chemic)
   deallocate(model_chemic_surf)
 
   ! -----------------------------------------------------------------------
   ! [4.0] Perform optional maximum error check:
   !-----------------------------------------------------------------------

   if ( check_max_iv ) &
      call da_check_max_iv_chem_sfc(iv,ob, it, num_qcstat_conv)

   if (trace_use) call da_trace_exit("da_get_innov_vector_chem_sfc")

end subroutine da_get_innov_vector_chem_sfc
