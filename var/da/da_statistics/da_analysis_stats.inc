subroutine da_analysis_stats  (grid, stats_unit, stats_unit2)
   
   !------------------------------------------------------------------------
   ! Purpose: Calculate min, max, mean and RMS of input 1d field.
   !------------------------------------------------------------------------
#if (WRF_CHEM == 1)
use module_state_description, only : num_chem, PARAM_FIRST_SCALAR                  &
!                     ,p_chem_ic_p25,  p_chem_ic_sulf                                &
!                     ,p_chem_ic_bc1, p_chem_ic_bc2, p_chem_ic_oc1, p_chem_ic_oc2    &
!                     ,p_chem_ic_dust_1, p_chem_ic_dust_2, p_chem_ic_seas_1, p_chem_ic_seas_2  &
!                     ,p_chem_ic_bc_a01, p_chem_ic_bc_a02, p_chem_ic_bc_a03    &
!                     ,p_chem_ic_oc_a01, p_chem_ic_oc_a02, p_chem_ic_oc_a03    &
!                     ,p_chem_ic_so4_a01, p_chem_ic_so4_a02, p_chem_ic_so4_a03    &
!                     ,p_chem_ic_no3_a01, p_chem_ic_no3_a02, p_chem_ic_no3_a03    &
!                     ,p_chem_ic_nh4_a01, p_chem_ic_nh4_a02, p_chem_ic_nh4_a03    &
!                     ,p_chem_ic_cl_a01, p_chem_ic_cl_a02, p_chem_ic_cl_a03    &
!                     ,p_chem_ic_na_a01, p_chem_ic_na_a02, p_chem_ic_na_a03    &
!                     ,p_chem_ic_oin_a01, p_chem_ic_oin_a02, p_chem_ic_oin_a03 &
                     ,p_chem_ic_so4aj,  p_chem_ic_so4ai,  p_chem_ic_nh4aj,  p_chem_ic_nh4ai,  &
                      p_chem_ic_no3aj,  p_chem_ic_no3ai,  p_chem_ic_naaj,   p_chem_ic_naai,   &
                      p_chem_ic_claj,   p_chem_ic_clai, &
                      p_chem_ic_asoa1j, p_chem_ic_asoa1i, p_chem_ic_asoa2j, p_chem_ic_asoa2i, &
                      p_chem_ic_asoa3j, p_chem_ic_asoa3i, p_chem_ic_asoa4j, p_chem_ic_asoa4i, &
                      p_chem_ic_bsoa1j, p_chem_ic_bsoa1i, p_chem_ic_bsoa2j, p_chem_ic_bsoa2i, &
                      p_chem_ic_bsoa3j, p_chem_ic_bsoa3i, p_chem_ic_bsoa4j, p_chem_ic_bsoa4i, &
                      p_chem_ic_orgpaj, p_chem_ic_orgpai, p_chem_ic_ecj,    p_chem_ic_eci,    &
                      p_chem_ic_p25j,   p_chem_ic_p25i,   p_chem_ic_antha,  p_chem_ic_seas,   &
                      p_chem_ic_soila,  p_chem_ic_so2,    p_chem_ic_no2,    p_chem_ic_o3,     p_chem_ic_co
#endif

   implicit none

   type (domain), intent (in) :: grid
   integer,       intent (in) :: stats_unit ! Output unit for stats.
   integer,       intent (in) :: stats_unit2 ! Output unit for chem stats.
   
   integer :: i, j, k
   integer :: ij_g, ijk_g   ! ij, ijk for global domain. 
   integer :: kdim  ! k range 

   real    :: um, vm, tm, pm, qm , qcwm, qrnm ! On local domain.
   real    :: qcim, qsnm, qgrm
   real    :: rij_g, rijk_g      ! On global domain.

   type (maxmin_field_type) :: max_u(kts:kte), max_v(kts:kte),     &
                               max_t(kts:kte), max_p(kts:kte),     &
                               max_q(kts:kte),                     &
                               min_u(kts:kte), min_v(kts:kte),     &
                               min_t(kts:kte), min_p(kts:kte),     &
                               min_q(kts:kte)

   type (maxmin_field_type) :: max_qcw(kts:kte), max_qrn(kts:kte), &
                               max_qci(kts:kte), max_qsn(kts:kte), &
                               max_qgr(kts:kte),                   &
                               min_qcw(kts:kte), min_qrn(kts:kte), &
                               min_qci(kts:kte), min_qsn(kts:kte), &
                               min_qgr(kts:kte)
#if (WRF_CHEM == 1)
   type (maxmin_field_type) :: max_chem(kts:kte,num_chem), min_chem(kts:kte,num_chem)
   real    :: chemm(num_chem), chemv(kts:kte,num_chem)
   integer :: ic
#endif

   real                     :: uv(kts:kte), vv(kts:kte),     &
                               tv(kts:kte), pv(kts:kte),     &
                               qv(kts:kte)
   real                     :: qcwv(kts:kte), qrnv(kts:kte), &
                               qciv(kts:kte), qsnv(kts:kte), &
                               qgrv(kts:kte)

   call da_trace_entry("da_analysis_stats")

   kdim = kte-kts+1

   ij_g = mix * mjy
   ijk_g = ij_g * mkz

   rij_g  = 1.0/real(ij_g)
   rijk_g = 1.0/real(ijk_g)

   if (rootproc) then
      write(unit=stats_unit, fmt='(/a/)') ' Minimum of gridded analysis increments'
      select case ( cloud_cv_options )
      case ( 0 )
         write(unit=stats_unit, fmt='(6a/)') &
              ' Lvl         ', &
              'u     i    j          ', &
              'v     i    j          ', &
              't     i    j          ', &
              'p     i    j          ', &
              'q     i    j'
      case ( 1 )
         write(unit=stats_unit, fmt='(8a/)') &
              ' Lvl         ', &
              'u     i    j          ', &
              'v     i    j          ', &
              't     i    j          ', &
              'p     i    j          ', &
              'q     i    j          ', &
              'qcw   i    j          ', &
              'qrn   i    j'
      case ( 2, 3 )
         write(unit=stats_unit, fmt='(11a/)') &
              ' Lvl         ', &
              'u     i    j          ', &
              'v     i    j          ', &
              't     i    j          ', &
              'p     i    j          ', &
              'q     i    j          ', &
              'qcw   i    j          ', &
              'qrn   i    j          ', &
              'qci   i    j          ', &
              'qsn   i    j          ', &
              'qgr   i    j'
      end select

#if (WRF_CHEM == 1)
      write(unit=stats_unit2, fmt='(/a/)') ' Minimum of gridded analysis increments'
      select case ( chem_cv_options )
       case ( 112 )
          write(unit=stats_unit2, fmt='(11a/)') &
               ' Lvl         ', &
               'p25      i    j          ', &
               'sulf     i    j          ', &
               'bc1      i    j          ', &
               'bc2      i    j          ', &
               'oc1      i    j          ', &
               'oc2      i    j          ', &
               'dust_1   i    j          ', &
               'dust_2   i    j          ', &
               'seas_1   i    j          ', &
               'seas_2   i    j'
       case ( 108 )
          write(unit=stats_unit2, fmt='(9a/)') &
               ' Lvl, so4aj i j, so4ai i j, nh4aj i j, nh4ai i j, no3aj i j,', &
               'no3ai  i j, naaj i j, naai i j, claj i j, clai i j,', &
               'asoa1j i j, asoa1i i j, asoa2j i j, asoa2i i j, ', &
               'asoa3j i j, asoa3i i j, asoa4j i j, asoa4i i j,', &
               'bsoa1j i j, bsoa1i i j, bsoa2j i j, bsoa2i i j, ', &
               'bsoa3j i j, bsoa3i i j, bsoa4j i j, bsoa4i i j,', &
               'orgpaj i j, orgpai i j, ecj i j,    eci i j,', &
               'p25j   i j,   p25i i j, antha  i j, seas   i j, soila i j,', &
               'so2 i j,    no2 i j,    o3 i j,     co i j'
      end select
#endif
   end if

   call da_maxmin_in_field(grid%xa%u(its:ite,jts:jte,kts:kte), max_u, min_u)
   call da_proc_maxmin_combine(kdim, max_u, min_u)
   call da_maxmin_in_field(grid%xa%v(its:ite,jts:jte,kts:kte), max_v, min_v)
   call da_proc_maxmin_combine(kdim, max_v, min_v)
   call da_maxmin_in_field(grid%xa%t(its:ite,jts:jte,kts:kte), max_t, min_t)
   call da_proc_maxmin_combine(kdim, max_t, min_t)
   call da_maxmin_in_field(grid%xa%p(its:ite,jts:jte,kts:kte), max_p, min_p)
   call da_proc_maxmin_combine(kdim, max_p, min_p)
   call da_maxmin_in_field(grid%xa%q(its:ite,jts:jte,kts:kte), max_q, min_q)
   call da_proc_maxmin_combine(kdim, max_q, min_q)

   if ( cloud_cv_options >= 1 ) then
      call da_maxmin_in_field(grid%xa%qcw(its:ite,jts:jte,kts:kte), max_qcw, min_qcw)
      call da_proc_maxmin_combine(kdim, max_qcw, min_qcw)
      call da_maxmin_in_field(grid%xa%qrn(its:ite,jts:jte,kts:kte), max_qrn, min_qrn)
      call da_proc_maxmin_combine(kdim, max_qrn, min_qrn)
   end if

   if ( cloud_cv_options >= 2 ) then
      call da_maxmin_in_field(grid%xa%qci(its:ite,jts:jte,kts:kte), max_qci, min_qci)
      call da_proc_maxmin_combine(kdim, max_qci, min_qci)
      call da_maxmin_in_field(grid%xa%qsn(its:ite,jts:jte,kts:kte), max_qsn, min_qsn)
      call da_proc_maxmin_combine(kdim, max_qsn, min_qsn)
      call da_maxmin_in_field(grid%xa%qgr(its:ite,jts:jte,kts:kte), max_qgr, min_qgr)
      call da_proc_maxmin_combine(kdim, max_qgr, min_qgr)
   end if

#if (WRF_CHEM == 1)
   if ( chem_cv_options >= 1 ) then
   do ic=PARAM_FIRST_SCALAR, num_chem
     call da_maxmin_in_field(grid%xachem%chem_ic(its:ite,jts:jte,kts:kte, ic), max_chem(:,ic), min_chem(:,ic))
     call da_proc_maxmin_combine(kdim, max_chem(:,ic), min_chem(:,ic))
   end do
   end if
#endif

   um = 999999.0
   vm = 999999.0
   tm = 999999.0
   pm = 999999.0
   qm = 999999.0

   qcwm = 999999.0
   qrnm = 999999.0
   qcim = 999999.0
   qsnm = 999999.0
   qgrm = 999999.0

#if (WRF_CHEM == 1)
   chemm = 999999.0
#endif

   do k = kts, kte
      if (rootproc) then
         if ( abs(min_q(k)%value) < 1.e-30 ) min_q(k)%value = 0.0
         select case ( cloud_cv_options )
         case ( 0 )
            write(unit=stats_unit, fmt='(i4,4(f12.4,2i5),e12.4,2i5)')    k, &
               min_u(k), min_v(k), min_t(k), min_p(k), min_q(k)
         case ( 1 )
            write(unit=stats_unit, fmt='(i4,4(f12.4,2i5),3(e12.4,2i5))') k, &
               min_u(k), min_v(k), min_t(k), min_p(k), min_q(k),            &
               min_qcw(k), min_qrn(k)
         case ( 2, 3 )
            write(unit=stats_unit, fmt='(i4,4(f12.4,2i5),6(e12.4,2i5))') k, &
               min_u(k), min_v(k), min_t(k), min_p(k), min_q(k),            &
               min_qcw(k), min_qrn(k), min_qci(k), min_qsn(k), min_qgr(k)
         end select

#if (WRF_CHEM == 1)
         select case ( chem_cv_options )
         !case ( 112 )   ! mozcart_da
         !   write(unit=stats_unit2, fmt='(i4,4(f12.4,2i5),6(e12.4,2i5))') k, &
         !      min_chem(k, p_chem_ic_p25), min_chem(k, p_chem_ic_sulf), &
         !      min_chem(k, p_chem_ic_bc1), min_chem(k, p_chem_ic_bc2), min_chem(k, p_chem_ic_oc1), min_chem(k, p_chem_ic_oc2), &
         !      min_chem(k, p_chem_ic_dust_1), min_chem(k, p_chem_ic_dust_2), min_chem(k, p_chem_ic_seas_1), min_chem(k, p_chem_ic_seas_1)
         case ( 108 )     ! racm_soa_vbs_da
            write(unit=stats_unit2, fmt='(i4,39(f12.4,2i5))') k, &
               min_chem(k, p_chem_ic_so4aj),  min_chem(k, p_chem_ic_so4ai),  min_chem(k, p_chem_ic_nh4aj),  min_chem(k, p_chem_ic_nh4ai), &
               min_chem(k, p_chem_ic_no3aj),  min_chem(k, p_chem_ic_no3ai),  min_chem(k, p_chem_ic_naaj),   min_chem(k, p_chem_ic_naai),  &
               min_chem(k, p_chem_ic_asoa1j), min_chem(k, p_chem_ic_asoa1i), min_chem(k, p_chem_ic_asoa2j), min_chem(k, p_chem_ic_asoa2i), &
               min_chem(k, p_chem_ic_asoa3j), min_chem(k, p_chem_ic_asoa3i), min_chem(k, p_chem_ic_asoa4j), min_chem(k, p_chem_ic_asoa4i), &
               min_chem(k, p_chem_ic_bsoa1j), min_chem(k, p_chem_ic_bsoa1i), min_chem(k, p_chem_ic_bsoa2j), min_chem(k, p_chem_ic_bsoa2i), &
               min_chem(k, p_chem_ic_bsoa3j), min_chem(k, p_chem_ic_bsoa3i), min_chem(k, p_chem_ic_bsoa4j), min_chem(k, p_chem_ic_bsoa4i), &
               min_chem(k, p_chem_ic_orgpaj), min_chem(k, p_chem_ic_orgpai), min_chem(k, p_chem_ic_ecj),    min_chem(k, p_chem_ic_eci),    &
               min_chem(k, p_chem_ic_p25j),   min_chem(k, p_chem_ic_p25i),   min_chem(k, p_chem_ic_antha),  min_chem(k, p_chem_ic_seas),   &
               min_chem(k, p_chem_ic_soila),  min_chem(k, p_chem_ic_claj),   min_chem(k, p_chem_ic_clai),  &
               min_chem(k, p_chem_ic_so2),    min_chem(k, p_chem_ic_no2),    min_chem(k, p_chem_ic_o3),     min_chem(k, p_chem_ic_co)
         !case ( 7 )     ! cbmz_mosaic_4bin_da
         !   write(unit=stats_unit2, fmt='(i4,24f12.4)') k, &
         !      min_chem(k, p_chem_ic_bc_a01), min_chem(k, p_chem_ic_bc_a02), min_chem(k, p_chem_ic_bc_a03), &
         !      min_chem(k, p_chem_ic_oc_a01), min_chem(k, p_chem_ic_oc_a02), min_chem(k, p_chem_ic_oc_a03), &
         !      min_chem(k, p_chem_ic_so4_a01), min_chem(k, p_chem_ic_so4_a02), min_chem(k, p_chem_ic_so4_a03), &
         !      min_chem(k, p_chem_ic_no3_a01), min_chem(k, p_chem_ic_no3_a02), min_chem(k, p_chem_ic_no3_a03), &
         !      min_chem(k, p_chem_ic_nh4_a01), min_chem(k, p_chem_ic_nh4_a02), min_chem(k, p_chem_ic_nh4_a03), &
         !      min_chem(k, p_chem_ic_cl_a01), min_chem(k, p_chem_ic_cl_a02), min_chem(k, p_chem_ic_cl_a03), &
         !      min_chem(k, p_chem_ic_na_a01), min_chem(k, p_chem_ic_na_a02), min_chem(k, p_chem_ic_na_a03), &
         !      min_chem(k, p_chem_ic_oin_a01), min_chem(k, p_chem_ic_oin_a02), min_chem(k, p_chem_ic_oin_a03)
         end select
#endif
      end if ! rootproc

      um=minval(min_u(:)%value)
      vm=minval(min_v(:)%value)
      tm=minval(min_t(:)%value)
      pm=minval(min_p(:)%value)
      qm=minval(min_q(:)%value)
      if ( cloud_cv_options >= 1 ) then
         qcwm=minval(min_qcw(:)%value)
         qrnm=minval(min_qrn(:)%value)
      end if
      if ( cloud_cv_options >= 2 ) then
         qcim=minval(min_qci(:)%value)
         qsnm=minval(min_qsn(:)%value)
         qgrm=minval(min_qgr(:)%value)
      end if
   end do

#if (WRF_CHEM == 1)
   do ic=PARAM_FIRST_SCALAR, num_chem
      chemm(ic)=minval(min_chem(:,ic)%value)
   end do
#endif

   if (rootproc) then
      select case ( cloud_cv_options )
      case ( 0 )
         write(unit=stats_unit, fmt='(a,4(f12.4,10x),e12.4)') ' ALL', &
            um, vm, tm, pm, qm
      case ( 1 )
         write(unit=stats_unit, fmt='(a,4(f12.4,10x),3(e12.4,10x))') ' ALL', &
            um, vm, tm, pm, qm, qcwm, qrnm
      case ( 2, 3 )
         write(unit=stats_unit, fmt='(a,4(f12.4,10x),6(e12.4,10x))') ' ALL', &
            um, vm, tm, pm, qm, qcwm, qrnm, qcim, qsnm, qgrm
      end select

#if (WRF_CHEM == 1)
      select case ( chem_cv_options )
      !case ( 112 )
      !    write(unit=stats_unit2, fmt='(a,4(f12.4,10x),6(e12.4,10x))') ' ALL', &
      !       chemm(p_chem_ic_p25), chemm(p_chem_ic_sulf), chemm(p_chem_ic_bc1), chemm(p_chem_ic_bc2), chemm(p_chem_ic_oc1), &
      !       chemm(p_chem_ic_oc2), chemm(p_chem_ic_dust_1), chemm(p_chem_ic_dust_2), chemm(p_chem_ic_seas_1), chemm(p_chem_ic_seas_2)
      case ( 108 )     ! racm_soa_vbs_da
          write(unit=stats_unit2, fmt='(a,39(f12.4,5x))') ' ALL', &
               chemm(p_chem_ic_so4aj),  chemm(p_chem_ic_so4ai),  chemm(p_chem_ic_nh4aj),  chemm(p_chem_ic_nh4ai), &
               chemm(p_chem_ic_no3aj),  chemm(p_chem_ic_no3ai),  chemm(p_chem_ic_naaj),   chemm(p_chem_ic_naai),  &
               chemm(p_chem_ic_asoa1j), chemm(p_chem_ic_asoa1i), chemm(p_chem_ic_asoa2j), chemm(p_chem_ic_asoa2i), &
               chemm(p_chem_ic_asoa3j), chemm(p_chem_ic_asoa3i), chemm(p_chem_ic_asoa4j), chemm(p_chem_ic_asoa4i), &
               chemm(p_chem_ic_bsoa1j), chemm(p_chem_ic_bsoa1i), chemm(p_chem_ic_bsoa2j), chemm(p_chem_ic_bsoa2i), &
               chemm(p_chem_ic_bsoa3j), chemm(p_chem_ic_bsoa3i), chemm(p_chem_ic_bsoa4j), chemm(p_chem_ic_bsoa4i), &
               chemm(p_chem_ic_orgpaj), chemm(p_chem_ic_orgpai), chemm(p_chem_ic_ecj),    chemm(p_chem_ic_eci),    &
               chemm(p_chem_ic_p25j),   chemm(p_chem_ic_p25i),   chemm(p_chem_ic_antha),  chemm(p_chem_ic_seas),   &
               chemm(p_chem_ic_soila),  chemm(p_chem_ic_claj),   chemm(p_chem_ic_clai),  &
               chemm(p_chem_ic_so2),    chemm(p_chem_ic_no2),    chemm(p_chem_ic_o3),     chemm(p_chem_ic_co)
      !case ( 7 )
      !    write(unit=stats_unit2, fmt='(a,24(f12.4,10x))') ' ALL', &
      !       chemm(p_chem_ic_bc_a01), chemm(p_chem_ic_bc_a02), chemm(p_chem_ic_bc_a03), &
      !       chemm(p_chem_ic_oc_a01), chemm(p_chem_ic_oc_a02), chemm(p_chem_ic_oc_a03), &
      !       chemm(p_chem_ic_so4_a01), chemm(p_chem_ic_so4_a02), chemm(p_chem_ic_so4_a03), &
      !       chemm(p_chem_ic_no3_a01), chemm(p_chem_ic_no3_a02), chemm(p_chem_ic_no3_a03), &
      !       chemm(p_chem_ic_nh4_a01), chemm(p_chem_ic_nh4_a02), chemm(p_chem_ic_nh4_a03), &
      !       chemm(p_chem_ic_cl_a01), chemm(p_chem_ic_cl_a02), chemm(p_chem_ic_cl_a03), &
      !       chemm(p_chem_ic_na_a01), chemm(p_chem_ic_na_a02), chemm(p_chem_ic_na_a03), &
      !       chemm(p_chem_ic_oin_a01), chemm(p_chem_ic_oin_a02), chemm(p_chem_ic_oin_a03)
      end select
#endif

      write(unit=stats_unit, fmt='(/a/)') &
         ' Maximum of gridded analysis increments'

#if (WRF_CHEM == 1)
      write(unit=stats_unit2, fmt='(/a/)') ' Maximum of gridded analysis increments'
      select case ( chem_cv_options )
       case ( 112 )
          write(unit=stats_unit2, fmt='(11a/)') &
               ' Lvl         ', &
               'so4aj    i    j          ', &
               'so4ai    i    j          ', &
               'nh4aj    i    j          ', &
               'nh4ai    i    j          ', &
               'no3aj    i    j          ', &
               'no3ai    i    j          ', &
               'naaj     i    j          ', &
               'naai     i    j          ', &
               'claj     i    j          ', &
               'clai     i    j          ', &
               'asoa1j   i    j          ', &
               'asoa1i   i    j          ', &
               'asoa2j   i    j          ', &
               'asoa2i   i    j          ', &
               'asoa3j   i    j          ', &
               'asoa3i   i    j          ', &
               'asoa4j   i    j          ', &
               'asoa4i   i    j          ', &
               'bsoa1j   i    j          ', &
               'bsoa1i   i    j          ', &
               'bsoa2j   i    j          ', &
               'bsoa2i   i    j          ', &
               'bsoa3j   i    j          ', &
               'bsoa3i   i    j          ', &
               'bsoa4j   i    j          ', &
               'bsoa4i   i    j          ', &
               'orgpaj   i    j          ', &
               'orgpai   i    j          ', &
               'ecj      i    j          ', &
               'eci      i    j          ', &
               'p25j     i    j          ', &
               'p25i     i    j          ', &
               'antha    i    j          ', &
               'seas     i    j          ', &
               'soila    i    j          ', &
               'so2      i    j          ', &
               'no2      i    j          ', &
               'o3       i    j          ', &
               'co       i    j'
      end select
#endif

      select case ( cloud_cv_options )
      case ( 0 )
         write(unit=stats_unit, fmt='(6a/)') &
            ' Lvl         ', &
            'u     i    j          ', &
            'v     i    j          ', &
            't     i    j          ', &
            'p     i    j          ', &
            'q     i    j'
      case ( 1 )
         write(unit=stats_unit, fmt='(8a/)') &
            ' Lvl         ', &
            'u     i    j          ', &
            'v     i    j          ', &
            't     i    j          ', &
            'p     i    j          ', &
            'q     i    j          ', &
            'qcw   i    j          ', &
            'qrn   i    j'
      case ( 2, 3 )
         write(unit=stats_unit, fmt='(11a/)') &
            ' Lvl         ', &
            'u     i    j          ', &
            'v     i    j          ', &
            't     i    j          ', &
            'p     i    j          ', &
            'q     i    j          ', &
            'qcw   i    j          ', &
            'qrn   i    j          ', &
            'qci   i    j          ', &
            'qsn   i    j          ', &
            'qgr   i    j'
      end select
   end if !rootproc

   do k = kts, kte
      if (rootproc) then
         if ( abs(max_q(k)%value) < 1.e-30 ) max_q(k)%value = 0.0
         select case ( cloud_cv_options )
         case ( 0 )
            write(unit=stats_unit, fmt='(i4,4(f12.4,2i5),e12.4,2i5)') k,    &
               max_u(k), max_v(k), max_t(k), max_p(k), max_q(k)
         case ( 1 )
            write(unit=stats_unit, fmt='(i4,4(f12.4,2i5),3(e12.4,2i5))') k, &
               max_u(k), max_v(k), max_t(k), max_p(k), max_q(k),            &
               max_qcw(k), max_qrn(k)
         case ( 2, 3 )
            write(unit=stats_unit, fmt='(i4,4(f12.4,2i5),6(e12.4,2i5))') k, &
               max_u(k), max_v(k), max_t(k), max_p(k), max_q(k),            &
               max_qcw(k), max_qrn(k), max_qci(k), max_qsn(k), max_qgr(k)
         end select

#if (WRF_CHEM == 1)
         select case ( chem_cv_options )
         !case ( 112 )
         !   write(unit=stats_unit2, fmt='(i4,4(f12.4,2i5),6(e12.4,2i5))') k, &
         !      max_chem(k, p_chem_ic_p25), max_chem(k, p_chem_ic_sulf), &
         !      max_chem(k, p_chem_ic_bc1), max_chem(k, p_chem_ic_bc2), max_chem(k, p_chem_ic_oc1), max_chem(k, p_chem_ic_oc2), &
         !      max_chem(k, p_chem_ic_dust_1), max_chem(k, p_chem_ic_dust_2), max_chem(k, p_chem_ic_seas_1), max_chem(k, p_chem_ic_seas_1)
         case ( 108 )     ! racm_soa_vbs_da
            write(unit=stats_unit2, fmt='(i4,39(f12.4,2i5))') k, &
               max_chem(k, p_chem_ic_so4aj),  max_chem(k, p_chem_ic_so4ai),  max_chem(k, p_chem_ic_nh4aj),  max_chem(k, p_chem_ic_nh4ai), &
               max_chem(k, p_chem_ic_no3aj),  max_chem(k, p_chem_ic_no3ai),  max_chem(k, p_chem_ic_naaj),   max_chem(k, p_chem_ic_naai),  &
               max_chem(k, p_chem_ic_asoa1j), max_chem(k, p_chem_ic_asoa1i), max_chem(k, p_chem_ic_asoa2j), max_chem(k, p_chem_ic_asoa2i), &
               max_chem(k, p_chem_ic_asoa3j), max_chem(k, p_chem_ic_asoa3i), max_chem(k, p_chem_ic_asoa4j), max_chem(k, p_chem_ic_asoa4i), &
               max_chem(k, p_chem_ic_bsoa1j), max_chem(k, p_chem_ic_bsoa1i), max_chem(k, p_chem_ic_bsoa2j), max_chem(k, p_chem_ic_bsoa2i), &
               max_chem(k, p_chem_ic_bsoa3j), max_chem(k, p_chem_ic_bsoa3i), max_chem(k, p_chem_ic_bsoa4j), max_chem(k, p_chem_ic_bsoa4i), &
               max_chem(k, p_chem_ic_orgpaj), max_chem(k, p_chem_ic_orgpai), max_chem(k, p_chem_ic_ecj),    max_chem(k, p_chem_ic_eci),    &
               max_chem(k, p_chem_ic_p25j),   max_chem(k, p_chem_ic_p25i),   max_chem(k, p_chem_ic_antha),  max_chem(k, p_chem_ic_seas),   &
               max_chem(k, p_chem_ic_soila),  max_chem(k, p_chem_ic_claj),   max_chem(k, p_chem_ic_clai),  &
               max_chem(k, p_chem_ic_so2),    max_chem(k, p_chem_ic_no2),    max_chem(k, p_chem_ic_o3),     max_chem(k, p_chem_ic_co)
         !case ( 7 )
         !   write(unit=stats_unit2, fmt='(i4,24(f12.4))') k, &
         !      max_chem(k, p_chem_ic_bc_a01), max_chem(k, p_chem_ic_bc_a02), max_chem(k, p_chem_ic_bc_a03), &
         !      max_chem(k, p_chem_ic_oc_a01), max_chem(k, p_chem_ic_oc_a02), max_chem(k, p_chem_ic_oc_a03), &
         !      max_chem(k, p_chem_ic_so4_a01), max_chem(k, p_chem_ic_so4_a02), max_chem(k, p_chem_ic_so4_a03), &
         !      max_chem(k, p_chem_ic_no3_a01), max_chem(k, p_chem_ic_no3_a02), max_chem(k, p_chem_ic_no3_a03), &
         !      max_chem(k, p_chem_ic_nh4_a01), max_chem(k, p_chem_ic_nh4_a02), max_chem(k, p_chem_ic_nh4_a03), &
         !      max_chem(k, p_chem_ic_cl_a01), max_chem(k, p_chem_ic_cl_a02), max_chem(k, p_chem_ic_cl_a03), &
         !      max_chem(k, p_chem_ic_na_a01), max_chem(k, p_chem_ic_na_a02), max_chem(k, p_chem_ic_na_a03), &
         !      max_chem(k, p_chem_ic_oin_a01), max_chem(k, p_chem_ic_oin_a02), max_chem(k, p_chem_ic_oin_a03)
         end select
#endif
      end if

      um=maxval(max_u(:)%value)
      vm=maxval(max_v(:)%value)
      tm=maxval(max_t(:)%value)
      pm=maxval(max_p(:)%value)
      qm=maxval(max_q(:)%value)
      if ( cloud_cv_options >= 1 ) then
         qcwm=maxval(max_qcw(:)%value)
         qrnm=maxval(max_qrn(:)%value)
      end if
      if ( cloud_cv_options >= 2 ) then
         qcim=maxval(max_qci(:)%value)
         qsnm=maxval(max_qsn(:)%value)
         qgrm=maxval(max_qgr(:)%value)
      end if
   end do

#if (WRF_CHEM == 1)
   do ic=PARAM_FIRST_SCALAR, num_chem
      chemm(ic)=maxval(max_chem(:,ic)%value)
   end do
#endif

   if (rootproc) then
      select case ( cloud_cv_options )
      case ( 0 )
        write(unit=stats_unit, fmt='(a,4(f12.4,10x),e12.4)') ' ALL', &
           um, vm, tm, pm, qm
      case ( 1 )
        write(unit=stats_unit, fmt='(a,4(f12.4,10x),3(e12.4,10x))') ' ALL', &
           um, vm, tm, pm, qm, qcwm, qrnm
      case ( 2, 3 )
        write(unit=stats_unit, fmt='(a,4(f12.4,10x),6(e12.4,10x))') ' ALL', &
           um, vm, tm, pm, qm, qcwm, qrnm, qcim, qsnm, qgrm
      end select

#if (WRF_CHEM == 1)
      select case ( chem_cv_options )
      !case ( 112 )
      !    write(unit=stats_unit2, fmt='(a,4(f12.4,10x),6(e12.4,10x))') ' ALL', &
      !       chemm(p_chem_ic_p25), chemm(p_chem_ic_sulf), chemm(p_chem_ic_bc1), chemm(p_chem_ic_bc2), chemm(p_chem_ic_oc1), &
      !       chemm(p_chem_ic_oc2), chemm(p_chem_ic_dust_1), chemm(p_chem_ic_dust_2), chemm(p_chem_ic_seas_1), chemm(p_chem_ic_seas_2)
      case ( 108 )     ! racm_soa_vbs_da
          write(unit=stats_unit2, fmt='(a,39(f12.4,10x))') ' ALL', &
               chemm(p_chem_ic_so4aj),  chemm(p_chem_ic_so4ai),  chemm(p_chem_ic_nh4aj),  chemm(p_chem_ic_nh4ai), &
               chemm(p_chem_ic_no3aj),  chemm(p_chem_ic_no3ai),  chemm(p_chem_ic_naaj),   chemm(p_chem_ic_naai),  &
               chemm(p_chem_ic_asoa1j), chemm(p_chem_ic_asoa1i), chemm(p_chem_ic_asoa2j), chemm(p_chem_ic_asoa2i), &
               chemm(p_chem_ic_asoa3j), chemm(p_chem_ic_asoa3i), chemm(p_chem_ic_asoa4j), chemm(p_chem_ic_asoa4i), &
               chemm(p_chem_ic_bsoa1j), chemm(p_chem_ic_bsoa1i), chemm(p_chem_ic_bsoa2j), chemm(p_chem_ic_bsoa2i), &
               chemm(p_chem_ic_bsoa3j), chemm(p_chem_ic_bsoa3i), chemm(p_chem_ic_bsoa4j), chemm(p_chem_ic_bsoa4i), &
               chemm(p_chem_ic_orgpaj), chemm(p_chem_ic_orgpai), chemm(p_chem_ic_ecj),    chemm(p_chem_ic_eci),    &
               chemm(p_chem_ic_p25j),   chemm(p_chem_ic_p25i),   chemm(p_chem_ic_antha),  chemm(p_chem_ic_seas),   &
               chemm(p_chem_ic_soila),  chemm(p_chem_ic_claj),   chemm(p_chem_ic_clai),  &
               chemm(p_chem_ic_so2),    chemm(p_chem_ic_no2),    chemm(p_chem_ic_o3),     chemm(p_chem_ic_co)
      !case ( 7 )
      !    write(unit=stats_unit2, fmt='(a,24(f12.4,10x))') ' ALL', &
      !       chemm(p_chem_ic_bc_a01), chemm(p_chem_ic_bc_a02), chemm(p_chem_ic_bc_a03),  &
      !       chemm(p_chem_ic_oc_a01), chemm(p_chem_ic_oc_a02), chemm(p_chem_ic_oc_a03),  &
      !       chemm(p_chem_ic_so4_a01), chemm(p_chem_ic_so4_a02), chemm(p_chem_ic_so4_a03),  &
      !       chemm(p_chem_ic_no3_a01), chemm(p_chem_ic_no3_a02), chemm(p_chem_ic_no3_a03),  &
      !       chemm(p_chem_ic_nh4_a01), chemm(p_chem_ic_nh4_a02), chemm(p_chem_ic_nh4_a03),  &
      !       chemm(p_chem_ic_cl_a01), chemm(p_chem_ic_cl_a02), chemm(p_chem_ic_cl_a03),  &
      !       chemm(p_chem_ic_na_a01), chemm(p_chem_ic_na_a02), chemm(p_chem_ic_na_a03),  &
      !       chemm(p_chem_ic_oin_a01), chemm(p_chem_ic_oin_a02), chemm(p_chem_ic_oin_a03)
      end select
#endif

      write(unit=stats_unit, fmt='(/a/)') ' Mean of gridded analysis increments'

#if (WRF_CHEM == 1)
      write(unit=stats_unit2, fmt='(/a/)') ' Mean of gridded analysis increments'
      select case ( chem_cv_options )
       case ( 112 )
          write(unit=stats_unit2, fmt='(a/)') &
            ' Lvl        p25           sulf           bc1           bc2           oc1           oc2         dust_1         dust_2         seas_1         seas_2'
      end select
#endif

      select case ( cloud_cv_options )
      case ( 0 )
         write(unit=stats_unit, fmt='(a/)') &
            ' Lvl        u           v           t           p           q'
      case ( 1 )
         write(unit=stats_unit, fmt='(a/)') &
            ' Lvl        u           v           t           p           q           qcw         qrn'
      case ( 2, 3 )
         write(unit=stats_unit, fmt='(a/)') &
            ' Lvl        u           v           t           p           q           qcw         qrn         qci         qsn         qgr'
      end select
   end if !rootproc

   um = 0.0
   vm = 0.0
   tm = 0.0
   pm = 0.0
   qm = 0.0

#if (WRF_CHEM == 1)
   chemm = 0.0
   if ( chem_cv_options >= 1 ) then
   do ic=PARAM_FIRST_SCALAR, num_chem
     do k = kts, kte
        chemv(k,ic) = sum(grid%xachem%chem_ic(its:ite,jts:jte,k,ic))
     end do
   call da_proc_sum_real (chemv(:,ic))
   end do
   end if
#endif

   do k = kts, kte
      uv(k) = sum(grid%xa%u(its:ite,jts:jte,k))
      vv(k) = sum(grid%xa%v(its:ite,jts:jte,k))
      tv(k) = sum(grid%xa%t(its:ite,jts:jte,k))
      pv(k) = sum(grid%xa%p(its:ite,jts:jte,k))
      qv(k) = sum(grid%xa%q(its:ite,jts:jte,k))
   end do
   call da_proc_sum_real (uv)
   call da_proc_sum_real (vv)
   call da_proc_sum_real (tv)
   call da_proc_sum_real (pv)
   call da_proc_sum_real (qv)

   if ( cloud_cv_options >= 1 ) then
      qcwm = 0.0
      qrnm = 0.0
       do k = kts, kte
          qcwv(k) = sum(grid%xa%qcw(its:ite,jts:jte,k))
          qrnv(k) = sum(grid%xa%qrn(its:ite,jts:jte,k))
       end do
      call da_proc_sum_real (qcwv)
      call da_proc_sum_real (qrnv)
   end if
   if ( cloud_cv_options >= 2 ) then
      qcim = 0.0
      qsnm = 0.0
      qgrm = 0.0
       do k = kts, kte
          qciv(k) = sum(grid%xa%qci(its:ite,jts:jte,k))
          qsnv(k) = sum(grid%xa%qsn(its:ite,jts:jte,k))
          qgrv(k) = sum(grid%xa%qgr(its:ite,jts:jte,k))
       end do
      call da_proc_sum_real (qciv)
      call da_proc_sum_real (qsnv)
      call da_proc_sum_real (qgrv)
   end if

   if (rootproc) then
      do k = kts, kte
         select case ( cloud_cv_options )
         case ( 0 )
            write(unit=stats_unit, fmt='(i4,4f12.4,e12.4)')  k, &
               uv(k)*rij_g, vv(k)*rij_g, tv(k)*rij_g, &
               pv(k)*rij_g, qv(k)*rij_g
         case ( 1 )
            write(unit=stats_unit, fmt='(i4,4f12.4,3e12.4)') k, &
               uv(k)*rij_g, vv(k)*rij_g, tv(k)*rij_g, &
               pv(k)*rij_g, qv(k)*rij_g,                        &
               qcwv(k)*rij_g, qrnv(k)*rij_g
         case ( 2, 3 )
            write(unit=stats_unit, fmt='(i4,4f12.4,6e12.4)') k, &
               uv(k)*rij_g, vv(k)*rij_g, tv(k)*rij_g, &
               pv(k)*rij_g, qv(k)*rij_g,                        &
               qcwv(k)*rij_g, qrnv(k)*rij_g, qciv(k)*rij_g,     &
               qsnv(k)*rij_g, qgrv(k)*rij_g
         end select

#if (WRF_CHEM == 1)
         select case ( chem_cv_options )
         !case ( 112 )
         !   write(unit=stats_unit2, fmt='(i4,10e14.4)') k, &
         !      chemv(k, p_chem_ic_p25)*rij_g, chemv(k, p_chem_ic_sulf)*rij_g, &
         !      chemv(k, p_chem_ic_bc1)*rij_g, chemv(k, p_chem_ic_bc2)*rij_g, chemv(k, p_chem_ic_oc1)*rij_g, chemv(k, p_chem_ic_oc2)*rij_g, &
         !      chemv(k, p_chem_ic_dust_1)*rij_g, chemv(k, p_chem_ic_dust_2)*rij_g, chemv(k, p_chem_ic_seas_2)*rij_g, chemv(k, p_chem_ic_seas_2)*rij_g
         case ( 108 )     ! racm_soa_vbs_da
            write(unit=stats_unit2, fmt='(i4,39f12.5)') k, &
               chemv(k, p_chem_ic_so4aj)*rij_g,  chemv(k, p_chem_ic_so4ai)*rij_g,  chemv(k, p_chem_ic_nh4aj)*rij_g,  chemv(k, p_chem_ic_nh4ai)*rij_g, &
               chemv(k, p_chem_ic_no3aj)*rij_g,  chemv(k, p_chem_ic_no3ai)*rij_g,  chemv(k, p_chem_ic_naaj)*rij_g,   chemv(k, p_chem_ic_naai)*rij_g,  &
               chemv(k, p_chem_ic_asoa1j)*rij_g, chemv(k, p_chem_ic_asoa1i)*rij_g, chemv(k, p_chem_ic_asoa2j)*rij_g, chemv(k, p_chem_ic_asoa2i)*rij_g, &
               chemv(k, p_chem_ic_asoa3j)*rij_g, chemv(k, p_chem_ic_asoa3i)*rij_g, chemv(k, p_chem_ic_asoa4j)*rij_g, chemv(k, p_chem_ic_asoa4i)*rij_g, &
               chemv(k, p_chem_ic_bsoa1j)*rij_g, chemv(k, p_chem_ic_bsoa1i)*rij_g, chemv(k, p_chem_ic_bsoa2j)*rij_g, chemv(k, p_chem_ic_bsoa2i)*rij_g, &
               chemv(k, p_chem_ic_bsoa3j)*rij_g, chemv(k, p_chem_ic_bsoa3i)*rij_g, chemv(k, p_chem_ic_bsoa4j)*rij_g, chemv(k, p_chem_ic_bsoa4i)*rij_g, &
               chemv(k, p_chem_ic_orgpaj)*rij_g, chemv(k, p_chem_ic_orgpai)*rij_g, chemv(k, p_chem_ic_ecj)*rij_g,    chemv(k, p_chem_ic_eci)*rij_g,    &
               chemv(k, p_chem_ic_p25j)*rij_g,   chemv(k, p_chem_ic_p25i)*rij_g,   chemv(k, p_chem_ic_antha)*rij_g,  chemv(k, p_chem_ic_seas)*rij_g,   &
               chemv(k, p_chem_ic_soila)*rij_g,  chemv(k, p_chem_ic_claj)*rij_g,   chemv(k, p_chem_ic_clai)*rij_g,  &
               chemv(k, p_chem_ic_so2)*rij_g,    chemv(k, p_chem_ic_no2)*rij_g,    chemv(k, p_chem_ic_o3)*rij_g,     chemv(k, p_chem_ic_co)*rij_g
         !case ( 7 )
         !   write(unit=stats_unit2, fmt='(i4,24e14.4)') k, &
         !      chemv(k, p_chem_ic_bc_a01)*rij_g, chemv(k, p_chem_ic_bc_a02)*rij_g, chemv(k, p_chem_ic_bc_a03)*rij_g,  &
         !      chemv(k, p_chem_ic_oc_a01)*rij_g, chemv(k, p_chem_ic_oc_a02)*rij_g, chemv(k, p_chem_ic_oc_a03)*rij_g,  &
         !      chemv(k, p_chem_ic_so4_a01)*rij_g, chemv(k, p_chem_ic_so4_a02)*rij_g, chemv(k, p_chem_ic_so4_a03)*rij_g,  &
         !      chemv(k, p_chem_ic_no3_a01)*rij_g, chemv(k, p_chem_ic_no3_a02)*rij_g, chemv(k, p_chem_ic_no3_a03)*rij_g,  &
         !      chemv(k, p_chem_ic_nh4_a01)*rij_g, chemv(k, p_chem_ic_nh4_a02)*rij_g, chemv(k, p_chem_ic_nh4_a03)*rij_g,  &
         !      chemv(k, p_chem_ic_cl_a01)*rij_g, chemv(k, p_chem_ic_cl_a02)*rij_g, chemv(k, p_chem_ic_cl_a03)*rij_g,  &
         !      chemv(k, p_chem_ic_na_a01)*rij_g, chemv(k, p_chem_ic_na_a02)*rij_g, chemv(k, p_chem_ic_na_a03)*rij_g,  &
         !      chemv(k, p_chem_ic_oin_a01)*rij_g, chemv(k, p_chem_ic_oin_a02)*rij_g, chemv(k, p_chem_ic_oin_a03)*rij_g
         end select
         do ic=PARAM_FIRST_SCALAR, num_chem
           chemm(ic)=chemm(ic)+chemv(k,ic)
         end do
#endif

         um=um+uv(k)
         vm=vm+vv(k)
         tm=tm+tv(k)
         pm=pm+pv(k)
         qm=qm+qv(k)
         if ( cloud_cv_options >= 1 ) then
            qcwm = qcwm + qcwv(k)
            qrnm = qrnm + qrnv(k)
         end if
         if ( cloud_cv_options >= 2 ) then
            qcim = qcim + qciv(k)
            qsnm = qsnm + qsnv(k)
            qgrm = qgrm + qgrv(k)
         end if
      end do !k loop
   end if !rootproc

   if (rootproc) then
      select case ( cloud_cv_options )
      case ( 0 )
         write(unit=stats_unit, fmt='(a,4f12.4,e12.4)')  ' ALL',   &
            um*rijk_g, vm*rijk_g, tm*rijk_g, pm*rijk_g, qm*rijk_g
      case ( 1 )
         write(unit=stats_unit, fmt='(a,4f12.4,3e12.4)') ' ALL',   &
            um*rijk_g, vm*rijk_g, tm*rijk_g, pm*rijk_g, qm*rijk_g, &
            qcwm*rijk_g, qrnm*rijk_g
      case ( 2, 3 )
         write(unit=stats_unit, fmt='(a,4f12.4,6e12.4)') ' ALL',   &
            um*rijk_g, vm*rijk_g, tm*rijk_g, pm*rijk_g, qm*rijk_g, &
            qcwm*rijk_g, qrnm*rijk_g, qcim*rijk_g, qsnm*rijk_g, qgrm*rijk_g
      end select

#if (WRF_CHEM == 1)
         select case ( chem_cv_options )
         !case ( 112 )
         !write(unit=stats_unit2, fmt='(a,10e14.4)') ' ALL',   &                                                                                      
         !   chemv(k, p_chem_ic_p25)*rijk_g, chemv(k, p_chem_ic_sulf)*rijk_g, &
         !   chemv(k, p_chem_ic_bc1)*rijk_g, chemv(k, p_chem_ic_bc2)*rijk_g, chemv(k, p_chem_ic_oc1)*rijk_g, chemv(k, p_chem_ic_oc2)*rijk_g, &
         !   chemv(k, p_chem_ic_dust_1)*rijk_g, chemv(k, p_chem_ic_dust_2)*rijk_g, chemv(k, p_chem_ic_seas_2)*rijk_g, chemv(k, p_chem_ic_seas_2)*rijk_g
         case ( 108 )     ! racm_soa_vbs_da
         write(unit=stats_unit2, fmt='(a,39f12.5)') ' ALL',   &                                                                                      
            chemv(k, p_chem_ic_so4aj)*rijk_g,  chemv(k, p_chem_ic_so4ai)*rijk_g,  chemv(k, p_chem_ic_nh4aj)*rijk_g,  chemv(k, p_chem_ic_nh4ai)*rijk_g, &
            chemv(k, p_chem_ic_no3aj)*rijk_g,  chemv(k, p_chem_ic_no3ai)*rijk_g,  chemv(k, p_chem_ic_naaj)*rijk_g,   chemv(k, p_chem_ic_naai)*rijk_g,  &
            chemv(k, p_chem_ic_asoa1j)*rijk_g, chemv(k, p_chem_ic_asoa1i)*rijk_g, chemv(k, p_chem_ic_asoa2j)*rijk_g, chemv(k, p_chem_ic_asoa2i)*rijk_g, &
            chemv(k, p_chem_ic_asoa3j)*rijk_g, chemv(k, p_chem_ic_asoa3i)*rijk_g, chemv(k, p_chem_ic_asoa4j)*rijk_g, chemv(k, p_chem_ic_asoa4i)*rijk_g, &
            chemv(k, p_chem_ic_bsoa1j)*rijk_g, chemv(k, p_chem_ic_bsoa1i)*rijk_g, chemv(k, p_chem_ic_bsoa2j)*rijk_g, chemv(k, p_chem_ic_bsoa2i)*rijk_g, &
            chemv(k, p_chem_ic_bsoa3j)*rijk_g, chemv(k, p_chem_ic_bsoa3i)*rijk_g, chemv(k, p_chem_ic_bsoa4j)*rijk_g, chemv(k, p_chem_ic_bsoa4i)*rijk_g, &
            chemv(k, p_chem_ic_orgpaj)*rijk_g, chemv(k, p_chem_ic_orgpai)*rijk_g, chemv(k, p_chem_ic_ecj)*rijk_g,    chemv(k, p_chem_ic_eci)*rijk_g,    &
            chemv(k, p_chem_ic_p25j)*rijk_g,   chemv(k, p_chem_ic_p25i)*rijk_g,   chemv(k, p_chem_ic_antha)*rijk_g,  chemv(k, p_chem_ic_seas)*rijk_g,   &
            chemv(k, p_chem_ic_soila)*rijk_g,  chemv(k, p_chem_ic_claj)*rijk_g,   chemv(k, p_chem_ic_clai)*rijk_g,  &
            chemv(k, p_chem_ic_so2)*rijk_g,    chemv(k, p_chem_ic_no2)*rijk_g,    chemv(k, p_chem_ic_o3)*rijk_g,     chemv(k, p_chem_ic_co)*rijk_g
         !case ( 7 )
         !write(unit=stats_unit2, fmt='(a,24e14.4)') ' ALL',   &
         !   chemv(k, p_chem_ic_bc_a01)*rijk_g, chemv(k, p_chem_ic_bc_a02)*rijk_g, chemv(k, p_chem_ic_bc_a03)*rijk_g, &
         !   chemv(k, p_chem_ic_oc_a01)*rijk_g, chemv(k, p_chem_ic_oc_a02)*rijk_g, chemv(k, p_chem_ic_oc_a03)*rijk_g, &
         !   chemv(k, p_chem_ic_so4_a01)*rijk_g, chemv(k, p_chem_ic_so4_a02)*rijk_g, chemv(k, p_chem_ic_so4_a03)*rijk_g, &
         !   chemv(k, p_chem_ic_no3_a01)*rijk_g, chemv(k, p_chem_ic_no3_a02)*rijk_g, chemv(k, p_chem_ic_no3_a03)*rijk_g, &
         !   chemv(k, p_chem_ic_nh4_a01)*rijk_g, chemv(k, p_chem_ic_nh4_a02)*rijk_g, chemv(k, p_chem_ic_nh4_a03)*rijk_g, &
         !   chemv(k, p_chem_ic_cl_a01)*rijk_g, chemv(k, p_chem_ic_cl_a02)*rijk_g, chemv(k, p_chem_ic_cl_a03)*rijk_g, &
         !   chemv(k, p_chem_ic_na_a01)*rijk_g, chemv(k, p_chem_ic_na_a02)*rijk_g, chemv(k, p_chem_ic_na_a03)*rijk_g, &
         !   chemv(k, p_chem_ic_oin_a01)*rijk_g, chemv(k, p_chem_ic_oin_a02)*rijk_g, chemv(k, p_chem_ic_oin_a03)*rijk_g
         end select
#endif


      write(unit=stats_unit, fmt='(/a/)') ' RMSE of gridded analysis increments'

#if (WRF_CHEM == 1)
      write(unit=stats_unit2, fmt='(/a/)') ' RMSE of gridded analysis increments'
      select case ( chem_cv_options )
       case ( 1 )
          write(unit=stats_unit2, fmt='(a/)') &
            ' Lvl        p25           sulf           bc1           bc2           oc1           oc2         dust_1         dust_2         seas_1         seas_2'
      end select
#endif

      select case ( cloud_cv_options )
      case ( 0 )
         write(unit=stats_unit, fmt='(a/)') &
            ' Lvl        u           v           t           p           q'
      case ( 1 )
         write(unit=stats_unit, fmt='(a/)') &
            ' Lvl        u           v           t           p           q           qcw         qrn'
      case ( 2, 3 )
         write(unit=stats_unit, fmt='(a/)') &
            ' Lvl        u           v           t           p           q           qcw         qrn         qci          qsn         qgr'
      end select
   end if !rootproc

   um = 0.0
   vm = 0.0
   tm = 0.0
   pm = 0.0
   qm = 0.0
   uv = 0.0
   vv = 0.0
   tv = 0.0
   pv = 0.0
   qv = 0.0

#if (WRF_CHEM == 1)
   chemm = 0.0
   chemv = 0.0
#endif

   qcwv = 0.0
   qrnv = 0.0
   qciv = 0.0
   qsnv = 0.0
   qgrv = 0.0
   qcwm = 0.0
   qrnm = 0.0
   qcim = 0.0
   qsnm = 0.0
   qgrm = 0.0

   do k = kts, kte
      do j=jts,jte
         do i=its,ite
            uv(k) = uv(k) + grid%xa%u(i,j,k) * grid%xa%u(i,j,k)
            vv(k) = vv(k) + grid%xa%v(i,j,k) * grid%xa%v(i,j,k)
            tv(k) = tv(k) + grid%xa%t(i,j,k) * grid%xa%t(i,j,k)
            pv(k) = pv(k) + grid%xa%p(i,j,k) * grid%xa%p(i,j,k)
            qv(k) = qv(k) + grid%xa%q(i,j,k) * grid%xa%q(i,j,k)
         end do
      end do
   end do
   call da_proc_sum_real (uv)
   call da_proc_sum_real (vv)
   call da_proc_sum_real (tv)
   call da_proc_sum_real (pv)
   call da_proc_sum_real (qv)

#if (WRF_CHEM == 1)
   do ic=PARAM_FIRST_SCALAR, num_chem
   do k = kts, kte
      do j=jts,jte
         do i=its,ite
            chemv(k,ic) = chemv(k,ic) + grid%xachem%chem_ic(i,j,k,ic) * grid%xachem%chem_ic(i,j,k,ic)
         end do
      end do
   call da_proc_sum_real (chemv(:,ic))
   end do
   end do
#endif

   if ( cloud_cv_options >= 1 ) then
      do k = kts, kte
         do j=jts,jte
            do i=its,ite
               qcwv(k) = qcwv(k) + grid%xa%qcw(i,j,k) * grid%xa%qcw(i,j,k)
               qrnv(k) = qrnv(k) + grid%xa%qrn(i,j,k) * grid%xa%qrn(i,j,k)
            end do
         end do
      end do
      call da_proc_sum_real (qcwv)
      call da_proc_sum_real (qrnv)
   end if

   if ( cloud_cv_options >= 2 ) then
      do k = kts, kte
         do j=jts,jte
            do i=its,ite
               qciv(k) = qciv(k) + grid%xa%qci(i,j,k) * grid%xa%qci(i,j,k)
               qsnv(k) = qsnv(k) + grid%xa%qsn(i,j,k) * grid%xa%qsn(i,j,k)
               qgrv(k) = qgrv(k) + grid%xa%qgr(i,j,k) * grid%xa%qgr(i,j,k)
            end do
         end do
      end do
      call da_proc_sum_real (qciv)
      call da_proc_sum_real (qsnv)
      call da_proc_sum_real (qgrv)
   end if

   if (rootproc) then
      do k = kts, kte
         select case ( cloud_cv_options )
         case ( 0 )
            write(unit=stats_unit, fmt='(i4,4f12.4,e12.4)')  k, &
               sqrt(uv(k)*rij_g),   &
               sqrt(vv(k)*rij_g),   &
               sqrt(tv(k)*rij_g),   &
               sqrt(pv(k)*rij_g),   &
               sqrt(qv(k)*rij_g)
         case ( 1 )
            write(unit=stats_unit, fmt='(i4,4f12.4,3e12.4)') k, &
               sqrt(uv(k)*rij_g),   &
               sqrt(vv(k)*rij_g),   &
               sqrt(tv(k)*rij_g),   &
               sqrt(pv(k)*rij_g),   &
               sqrt(qv(k)*rij_g),   &
               sqrt(qcwv(k)*rij_g), &
               sqrt(qrnv(k)*rij_g)
         case ( 2, 3 )
            write(unit=stats_unit, fmt='(i4,4f12.4,6e12.4)') k, &
               sqrt(uv(k)*rij_g),   &
               sqrt(vv(k)*rij_g),   &
               sqrt(tv(k)*rij_g),   &
               sqrt(pv(k)*rij_g),   &
               sqrt(qv(k)*rij_g),   &
               sqrt(qcwv(k)*rij_g), &
               sqrt(qrnv(k)*rij_g), &
               sqrt(qciv(k)*rij_g), &
               sqrt(qsnv(k)*rij_g), &
               sqrt(qgrv(k)*rij_g)
         end select

#if (WRF_CHEM == 1)
         select case ( chem_cv_options )
         !case ( 112 )
         !   write(unit=stats_unit2, fmt='(i4,10e14.4)') k, &
         !      sqrt(chemv(k, p_chem_ic_p25)*rij_g),   &
         !      sqrt(chemv(k, p_chem_ic_sulf)*rij_g),   &
         !      sqrt(chemv(k, p_chem_ic_bc1)*rij_g),   &
         !      sqrt(chemv(k, p_chem_ic_bc2)*rij_g),   &
         !      sqrt(chemv(k, p_chem_ic_oc1)*rij_g),   &
         !      sqrt(chemv(k, p_chem_ic_oc2)*rij_g), &
         !      sqrt(chemv(k, p_chem_ic_dust_1)*rij_g), &
         !      sqrt(chemv(k, p_chem_ic_dust_2)*rij_g), &
         !      sqrt(chemv(k, p_chem_ic_seas_1)*rij_g), &
         !      sqrt(chemv(k, p_chem_ic_seas_2)*rij_g)
         case ( 108 )
            write(unit=stats_unit2, fmt='(i4,39f12.3)') k, &
               sqrt(chemv(k, p_chem_ic_so4aj)*rij_g), &
               sqrt(chemv(k, p_chem_ic_so4ai)*rij_g), &
               sqrt(chemv(k, p_chem_ic_nh4aj)*rij_g), &
               sqrt(chemv(k, p_chem_ic_nh4ai)*rij_g), &
               sqrt(chemv(k, p_chem_ic_no3aj)*rij_g), &
               sqrt(chemv(k, p_chem_ic_no3ai)*rij_g), &
               sqrt(chemv(k, p_chem_ic_naaj)*rij_g), &
               sqrt(chemv(k, p_chem_ic_naai)*rij_g), &
               sqrt(chemv(k, p_chem_ic_claj)*rij_g), &
               sqrt(chemv(k, p_chem_ic_clai)*rij_g), &
               sqrt(chemv(k, p_chem_ic_asoa1j)*rij_g), &
               sqrt(chemv(k, p_chem_ic_asoa1i)*rij_g), &
               sqrt(chemv(k, p_chem_ic_asoa2j)*rij_g), &
               sqrt(chemv(k, p_chem_ic_asoa2i)*rij_g), &
               sqrt(chemv(k, p_chem_ic_asoa3j)*rij_g), &
               sqrt(chemv(k, p_chem_ic_asoa3i)*rij_g), &
               sqrt(chemv(k, p_chem_ic_asoa4j)*rij_g), &
               sqrt(chemv(k, p_chem_ic_asoa4i)*rij_g), &
               sqrt(chemv(k, p_chem_ic_bsoa1j)*rij_g), &
               sqrt(chemv(k, p_chem_ic_bsoa1i)*rij_g), &
               sqrt(chemv(k, p_chem_ic_bsoa2j)*rij_g), &
               sqrt(chemv(k, p_chem_ic_bsoa2i)*rij_g), &
               sqrt(chemv(k, p_chem_ic_bsoa3j)*rij_g), &
               sqrt(chemv(k, p_chem_ic_bsoa3i)*rij_g), &
               sqrt(chemv(k, p_chem_ic_bsoa4j)*rij_g), &
               sqrt(chemv(k, p_chem_ic_bsoa4i)*rij_g), &
               sqrt(chemv(k, p_chem_ic_orgpaj)*rij_g), &
               sqrt(chemv(k, p_chem_ic_orgpai)*rij_g), &
               sqrt(chemv(k, p_chem_ic_ecj)*rij_g), &
               sqrt(chemv(k, p_chem_ic_eci)*rij_g), &
               sqrt(chemv(k, p_chem_ic_p25j)*rij_g), &
               sqrt(chemv(k, p_chem_ic_p25i)*rij_g), &
               sqrt(chemv(k, p_chem_ic_antha)*rij_g), &
               sqrt(chemv(k, p_chem_ic_seas)*rij_g), &
               sqrt(chemv(k, p_chem_ic_soila)*rij_g), &
               sqrt(chemv(k, p_chem_ic_so2)*rij_g), &
               sqrt(chemv(k, p_chem_ic_no2)*rij_g), &
               sqrt(chemv(k, p_chem_ic_o3)*rij_g), &
               sqrt(chemv(k, p_chem_ic_co)*rij_g)
         !case ( 7 )
         !   write(unit=stats_unit2, fmt='(i4,24e14.4)') k, &
         !      sqrt(chemv(k, p_chem_ic_bc_a01)*rij_g), sqrt(chemv(k, p_chem_ic_bc_a02)*rij_g), sqrt(chemv(k, p_chem_ic_bc_a03)*rij_g),  &
         !      sqrt(chemv(k, p_chem_ic_oc_a01)*rij_g), sqrt(chemv(k, p_chem_ic_oc_a02)*rij_g), sqrt(chemv(k, p_chem_ic_oc_a03)*rij_g),  &
         !      sqrt(chemv(k, p_chem_ic_so4_a01)*rij_g), sqrt(chemv(k, p_chem_ic_so4_a02)*rij_g), sqrt(chemv(k, p_chem_ic_so4_a03)*rij_g),  &
         !      sqrt(chemv(k, p_chem_ic_no3_a01)*rij_g), sqrt(chemv(k, p_chem_ic_no3_a02)*rij_g), sqrt(chemv(k, p_chem_ic_no3_a03)*rij_g),  &
         !      sqrt(chemv(k, p_chem_ic_nh4_a01)*rij_g), sqrt(chemv(k, p_chem_ic_nh4_a02)*rij_g), sqrt(chemv(k, p_chem_ic_nh4_a03)*rij_g),  &
         !      sqrt(chemv(k, p_chem_ic_cl_a01)*rij_g), sqrt(chemv(k, p_chem_ic_cl_a02)*rij_g), sqrt(chemv(k, p_chem_ic_cl_a03)*rij_g),  &
         !      sqrt(chemv(k, p_chem_ic_na_a01)*rij_g), sqrt(chemv(k, p_chem_ic_na_a02)*rij_g), sqrt(chemv(k, p_chem_ic_na_a03)*rij_g),  &
         !      sqrt(chemv(k, p_chem_ic_oin_a01)*rij_g), sqrt(chemv(k, p_chem_ic_oin_a02)*rij_g), sqrt(chemv(k, p_chem_ic_oin_a03)*rij_g)
         end select
         do ic=PARAM_FIRST_SCALAR, num_chem
           chemm(ic)=chemm(ic)+chemv(k,ic)
         end do
#endif

         um=um+uv(k)
         vm=vm+vv(k)
         tm=tm+tv(k)
         pm=pm+pv(k)
         qm=qm+qv(k)
         if ( cloud_cv_options >= 1 ) then
            qcwm=qcwm+qcwv(k)
            qrnm=qrnm+qrnv(k)
         end if
         if ( cloud_cv_options >= 2 ) then
            qcim=qcim+qciv(k)
            qsnm=qsnm+qsnv(k)
         end if
      end do !k loop
   end if !rootproc

   if (rootproc) then
      select case ( cloud_cv_options )
      case ( 0 )
         write(unit=stats_unit, fmt='(a,4f12.4,e12.4)')  ' ALL', &
            sqrt(um*rijk_g), sqrt(vm*rijk_g), sqrt(tm*rijk_g),   &
            sqrt(pm*rijk_g), sqrt(qm*rijk_g)
      case ( 1 )
         write(unit=stats_unit, fmt='(a,4f12.4,3e12.4)') ' ALL', &
            sqrt(um*rijk_g), sqrt(vm*rijk_g), sqrt(tm*rijk_g),   &
            sqrt(pm*rijk_g), sqrt(qm*rijk_g),                    &
            sqrt(qcwm*rijk_g), sqrt(qrnm*rijk_g)
      case ( 2, 3 )
         write(unit=stats_unit, fmt='(a,4f12.4,6e12.4)') ' ALL', &
            sqrt(um*rijk_g), sqrt(vm*rijk_g), sqrt(tm*rijk_g),   &
            sqrt(pm*rijk_g), sqrt(qm*rijk_g),                    &
            sqrt(qcwm*rijk_g), sqrt(qrnm*rijk_g), sqrt(qcim*rijk_g), &
            sqrt(qsnm*rijk_g), sqrt(qgrm*rijk_g)
      end select

#if (WRF_CHEM == 1)
         select case ( chem_cv_options )
         !case ( 112 )
         !write(unit=stats_unit2, fmt='(a,10e14.4)') ' ALL', &
         !sqrt(chemm(p_chem_ic_p25)*rijk_g), sqrt(chemm(p_chem_ic_sulf)*rijk_g),   &
         !sqrt(chemm(p_chem_ic_bc1)*rijk_g), sqrt(chemm(p_chem_ic_bc2)*rijk_g), sqrt(chemm(p_chem_ic_oc1)*rijk_g), sqrt(chemm(p_chem_ic_oc2)*rijk_g), &
         !sqrt(chemm(p_chem_ic_dust_1)*rijk_g), sqrt(chemm(p_chem_ic_dust_2)*rijk_g), &
         !sqrt(chemm(p_chem_ic_seas_1)*rijk_g), sqrt(chemm(p_chem_ic_seas_2)*rijk_g)
         case ( 108 )
         write(unit=stats_unit2, fmt='(a,39f12.5)') ' ALL', &
         sqrt(chemm(p_chem_ic_so4aj)*rijk_g),  sqrt(chemm(p_chem_ic_so4ai)*rijk_g),  &
         sqrt(chemm(p_chem_ic_nh4aj)*rijk_g),  sqrt(chemm(p_chem_ic_nh4ai)*rijk_g),  &
         sqrt(chemm(p_chem_ic_no3aj)*rijk_g),  sqrt(chemm(p_chem_ic_no3ai)*rijk_g),  &
         sqrt(chemm(p_chem_ic_naaj)*rijk_g),   sqrt(chemm(p_chem_ic_naai)*rijk_g),   &
         sqrt(chemm(p_chem_ic_claj)*rijk_g),   sqrt(chemm(p_chem_ic_clai)*rijk_g),   &
         sqrt(chemm(p_chem_ic_asoa1j)*rijk_g), sqrt(chemm(p_chem_ic_asoa1i)*rijk_g), &
         sqrt(chemm(p_chem_ic_asoa2j)*rijk_g), sqrt(chemm(p_chem_ic_asoa2i)*rijk_g), &
         sqrt(chemm(p_chem_ic_asoa3j)*rijk_g), sqrt(chemm(p_chem_ic_asoa3i)*rijk_g), &
         sqrt(chemm(p_chem_ic_asoa4j)*rijk_g), sqrt(chemm(p_chem_ic_asoa4i)*rijk_g), &
         sqrt(chemm(p_chem_ic_bsoa1j)*rijk_g), sqrt(chemm(p_chem_ic_bsoa1i)*rijk_g), &
         sqrt(chemm(p_chem_ic_bsoa2j)*rijk_g), sqrt(chemm(p_chem_ic_bsoa2i)*rijk_g), &
         sqrt(chemm(p_chem_ic_bsoa3j)*rijk_g), sqrt(chemm(p_chem_ic_bsoa3i)*rijk_g), &
         sqrt(chemm(p_chem_ic_bsoa4j)*rijk_g), sqrt(chemm(p_chem_ic_bsoa4i)*rijk_g), &
         sqrt(chemm(p_chem_ic_orgpaj)*rijk_g), sqrt(chemm(p_chem_ic_orgpai)*rijk_g), &
         sqrt(chemm(p_chem_ic_ecj)*rijk_g),    sqrt(chemm(p_chem_ic_eci)*rijk_g),    &
         sqrt(chemm(p_chem_ic_p25j)*rijk_g),   sqrt(chemm(p_chem_ic_p25i)*rijk_g),   &
         sqrt(chemm(p_chem_ic_antha)*rijk_g),  sqrt(chemm(p_chem_ic_seas)*rijk_g),   &
         sqrt(chemm(p_chem_ic_soila)*rijk_g),  sqrt(chemm(p_chem_ic_so2)*rijk_g),    &
         sqrt(chemm(p_chem_ic_no2)*rijk_g),    sqrt(chemm(p_chem_ic_o3)*rijk_g),     &
         sqrt(chemm(p_chem_ic_co)*rijk_g)
         !case ( 7 )
         !write(unit=stats_unit2, fmt='(a,24e14.4)') ' ALL', &
         !sqrt(chemm(p_chem_ic_bc_a01)*rijk_g), sqrt(chemm(p_chem_ic_bc_a02)*rijk_g), sqrt(chemm(p_chem_ic_bc_a03)*rijk_g), &
         !sqrt(chemm(p_chem_ic_oc_a01)*rijk_g), sqrt(chemm(p_chem_ic_oc_a02)*rijk_g), sqrt(chemm(p_chem_ic_oc_a03)*rijk_g), &
         !sqrt(chemm(p_chem_ic_so4_a01)*rijk_g), sqrt(chemm(p_chem_ic_so4_a02)*rijk_g), sqrt(chemm(p_chem_ic_so4_a03)*rijk_g), &
         !sqrt(chemm(p_chem_ic_no3_a01)*rijk_g), sqrt(chemm(p_chem_ic_no3_a02)*rijk_g), sqrt(chemm(p_chem_ic_no3_a03)*rijk_g), &
         !sqrt(chemm(p_chem_ic_nh4_a01)*rijk_g), sqrt(chemm(p_chem_ic_nh4_a02)*rijk_g), sqrt(chemm(p_chem_ic_nh4_a03)*rijk_g), &
         !sqrt(chemm(p_chem_ic_cl_a01)*rijk_g), sqrt(chemm(p_chem_ic_cl_a02)*rijk_g), sqrt(chemm(p_chem_ic_cl_a03)*rijk_g), &
         !sqrt(chemm(p_chem_ic_na_a01)*rijk_g), sqrt(chemm(p_chem_ic_na_a02)*rijk_g), sqrt(chemm(p_chem_ic_na_a03)*rijk_g), &
         !sqrt(chemm(p_chem_ic_oin_a01)*rijk_g), sqrt(chemm(p_chem_ic_oin_a02)*rijk_g), sqrt(chemm(p_chem_ic_oin_a03)*rijk_g)
         end select
#endif
   end if

   call da_trace_exit("da_analysis_stats")

contains

#include "da_maxmin_in_field.inc"

end subroutine da_analysis_stats


