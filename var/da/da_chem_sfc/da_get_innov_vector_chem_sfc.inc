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

   real, allocatable :: model_chemic_surf(:,:)
   real, allocatable :: model_chemic(:,:)
   real, allocatable :: model_rho(:,:)

   real, parameter  :: scaleL = 3000.         ! scaling factor for obs error [meters]
   real, parameter  :: maxomb = 100.          ! maximum threshold for abs(o-b)
   real             :: err_o, err_r

   if (trace_use) call da_trace_entry("da_get_innov_vector_chem_sfc")

   write(*,*) &
   'da_get_innov_vector_chem_sfc: PARAM_FIRST_SCALAR, num_chemic_surf, num_chem, scaleL, maxomb =', &
                                  PARAM_FIRST_SCALAR, num_chemic_surf, num_chem, scaleL, maxomb

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

   if (chemicda_opt == 4 .or. chemicda_opt == 5) then     ! unit conversion - common for all chem_cv_options
       model_chemic_surf(:,p_chemsi_so2)=model_rho(:,p_chem_ic_so2)*(64.06/28.964*1000)*(model_chemic(:,p_chem_ic_so2))
       model_chemic_surf(:,p_chemsi_no2)=model_rho(:,p_chem_ic_no2)*(46.01/28.964*1000)*(model_chemic(:,p_chem_ic_no2))
       model_chemic_surf(:,p_chemsi_o3) =model_rho(:,p_chem_ic_o3 )*(48.00/28.964*1000)*(model_chemic(:,p_chem_ic_o3 ))
       model_chemic_surf(:,p_chemsi_co) =model_rho(:,p_chem_ic_co )*(28.01/28.964*1000)*(model_chemic(:,p_chem_ic_co ))
   endif !(chemicda_opt == 4 ) then

   if (chem_cv_options == 108) then  ! racm_soa_vbs_da

      if (chemicda_opt <= 4 ) then   ! pm assimilation

      ! Initialize first
      model_chemic_surf(:,p_chemsi_pm25) = 0.

      ! pm2.5
      do n = p_chem_ic_so4aj, p_chem_ic_p25i 
      model_chemic_surf(:,p_chemsi_pm25)=model_chemic_surf(:,p_chemsi_pm25)+model_chemic(:,n)*model_rho(:,n)
      enddo

      ! pm10
      if (chemicda_opt > 1 .and. chemicda_opt <= 4 ) then
          model_chemic_surf(:,p_chemsi_pm10)=model_chemic_surf(:,p_chemsi_pm25) + &
          model_rho(:,p_chem_ic_so4aj)*(model_chemic(:,p_chem_ic_antha)+model_chemic(:,p_chem_ic_seas)+model_chemic(:,p_chem_ic_soila))
      endif !(chemicda_opt > 1)

      endif !(chemicda_opt <= 4 )

   end if   !(chem_cv_options == 108) then  ! racm_soa_vbs_da

   do ichem = PARAM_FIRST_SCALAR, num_chemic_surf
 
      do n=iv%info(chemic_surf)%n1,iv%info(chemic_surf)%n2

         ! [2.0] Initialise components of innovation vector: 

         iv % chemic_surf(n) % chem(ichem) % inv   = missing_r
         iv % chemic_surf(n) % chem(ichem) % qc    = missing
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
                iv % chemic_surf(n) % chem(ichem) % inv   = missing_r
                iv % chemic_surf(n) % chem(ichem) % qc    = missing
                iv % chemic_surf(n) % chem(ichem) % error = missing_r
            end if
            if (chemicda_opt == 4 .or. chemicda_opt == 5 .or. chemicda_opt == 9) then   ! gas assimilation
              if (ichem  == num_chemic_surf ) then    ! co assimilation
                if (abs(iv % chemic_surf(n) % chem(ichem) % inv) .ge. maxomb*1000) then
                    iv % chemic_surf(n) % chem(ichem) % inv   = missing_r
                    iv % chemic_surf(n) % chem(ichem) % qc    = missing
		    iv % chemic_surf(n) % chem(ichem) % error = missing_r
                end if
              end if
            end if
            ! ----------------

         end if	   ! good obs      
      end do
   end do

   deallocate(model_chemic)
   deallocate(model_chemic_surf)
   deallocate(model_rho)
 
   ! -----------------------------------------------------------------------
   ! [4.0] Perform optional maximum error check:
   !-----------------------------------------------------------------------

   if ( check_max_iv ) &
      call da_check_max_iv_chem_sfc(iv,ob, it, num_qcstat_conv)

   if (trace_use) call da_trace_exit("da_get_innov_vector_chem_sfc")

end subroutine da_get_innov_vector_chem_sfc


