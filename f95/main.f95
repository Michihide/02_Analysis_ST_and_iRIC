program main
  use mndr

  case  = '20170906_40_5x5'          ! <- folder name

  call set_parameter
  call read_fln
  call read_arry
  call set_pss
  call read_roughness
  call cal_critical_shear_stress
    
  ! ---- read vtk for dz/dt & dudt(obs) ----
  do n = 1, nend
     call set_arry
     call read_schalar_vtk(dev_bl, 4, 3, pss_st)   ; bl2(:,:,n)    = dev_bl(:,:)
     call read_schalar_vtk(u_cal , 6, 5, pss_iric) ; u_cal2(:,:,n) = u_cal(:,:)
     call deallocate
  enddo
  
  do n = 1, nend, intbl
     call set_arry
     
     ! ---- read vtk ----
     call read_schalar_vtk(dep_obs, 4, 1, pss_st)    ;  call read_schalar_vtk(wl_obs,  4, 2, pss_st)
     call read_schalar_vtk(dev_bl,  4, 3, pss_st)    ;  call read_schalar_vtk(bl    ,  4, 4, pss_st)

!!$     call read_schalar_vtk(dep_obs, 4, 2, pss_st)    ;  call read_schalar_vtk(wl_obs,  4, 1, pss_st)
!!$     call read_schalar_vtk(dev_bl,  4, 4, pss_st)    ;  call read_schalar_vtk(bl    ,  4, 3, pss_st)

!!$     call read_schalar_vtk(dwl_obs, 5, 1, pss_st)    ;  call read_schalar_vtk(wl_obs,  5, 2, pss_st)
!!$     call read_schalar_vtk(dev_bl,  5, 3, pss_st)    ;  call read_schalar_vtk(bl    ,  5, 4, pss_st)
!!$     call read_schalar_vtk(dep_obs, 5, 5, pss_st)

     call read_schalar_vtk(dep_cal, 6, 1, pss_iric)  ;  call read_schalar_vtk(wl_cal,  6, 2, pss_iric)
     call read_schalar_vtk(u_cal,   6, 5, pss_iric)  ;  call read_schalar_vtk(v_cal ,  6, 6, pss_iric)

     dwl_obs(1,:) = wl_obs(1,:)
     do i = 2, iend
        dwl_obs(i,:) = wl_obs(i,:) + 0.05*i*(1./200)
     enddo
     
     ! ---- cal slope & curvature ----
     call calculate_slope(wl_obs, iwi_obs, iwj_obs, iw_obs)
     call calculate_slope(wl_cal, iwi_cal, iwj_cal, iw_cal)
     call calculate_slope(bl    , ibi    , ibj    , ib)
     call calculate_curvature(bl, cbj    , cbi    , cb)
     call calculate_energy_slope(nc, dep_cal, u_cal, v_cal, iei, iej, ie)
     vel_cal(:,:)    = sqrt(u_cal(:,:)**2 + v_cal(:,:)**2)
     fr(:,:)         = u_cal(:,:) / sqrt(g * dep_cal(:,:))

     
     ! ---- cal shields_number &  M_equation ----
     call calculate_shields_number(dep_cal, iei, iej, taue_i, taue_j, taue) ; tautau(:,:) = taue_i(:,:) / tauc
     call calculate_Exner_equation_mpm(nc, dep_cal, u_cal, qb_mpm, dz_mpm)
     call calculate_M_equation_mpm_s (nc, dep_cal, iei, ibi, u_cal, Pe_s, Ms, dz_Ms)
     call calculate_M_equation_mpm_u (nc, dep_cal, iei, ibi, u_cal, Pe_u, Mu, dz_Mu)
     call calculate_M_equation_mpm_df(nc, dep_cal, iei, ibi, cbi, u_cal, Md, Df, Pe_d, dz_Md)

     do j = 1, jend
        ave_crs_dbl(j) = sum(dev_bl(:,j)) / iend
        ave_crs_dwl(j) = sum(dwl_obs(:,j)) / iend
     enddo
     
!!$     do i = 1, iend
!!$        do j = 1, jend
!!$           write(*,*)(dz_Mu(i,j)/(1)), (Mu(i,j)*ibi(i,j)), Mu(i,j)*iei(i,j), & 
!!$                Mu(i,j) / 2 * cbi(i,j) * ( ( x(i,j) - x(i+1,j) ) - Mu(i,j)*(1) )
!!$           write(*,*)abs(dz_Mu(i,j)/((time2(n+1)-time2(n)) * 60)), abs(Mu(i,j)*ibi(i,j)), mu(i,df)*iei(i,j), & 
!!$                abs(Mu(i,j) / 2 * cbi(i,j) * (abs(x(i,j)-x(i+1,j)) - Mu(i,j)*((time2(n+1)-time2(n))*60)))
!!$        enddo
!!$     enddo
     
     do i = 1, iend
        do j = 1, jend
           lcl(i,j) = dz_Md(i,j) / (10**3 * 600)
           adv(i,j) = Md(i,j) * ibi(i,j)
           dff(i,j) = Df(i,j) * cbi(i,j)
           frc(i,j) = Md(i,j) * iei(i,j)
        enddo
     enddo
     
!!$     Md(:,:) = Md(:,:) * 10**3
     
!!$     call calculate_M_equation_mpm_2d(nc, dep_cal, iei, iej, ibi, ibj, u_cal, v_cal, dhdx_mpm, dhdy_mpm, Mx, My, Mxy, dz_mpm)
!!$     call cal_nondim_celerity(taue_i, fr, iei, Mu, tau_nd, fr_nd, Mu_nd)
          
     call calclate_and_compare_dzdt_dhdx     
     call cmp_each_item

     count = 0
!!$     call filter(dif_dep, Mu)
     call filter(Fr, Mu)
     ave_Mu(n) = sum(Mu(:,:)) / (iend*jend-count)

! ---- output_vtk ----
     call output_schalar_vtk
     call out_vector_vtk

! ---- drw_cntr ----
     call drw_cntr
     
     call deallocate
  end do
  
end program main




subroutine set_arry
  use mndr
  
  ! ---- read_schalar_vtk & read_vector_vtk ----
  allocate(x(iend,jend), y(iend,jend))
  allocate(dep_obs(iend,jend), wl_obs(iend,jend), dev_bl(iend,jend), bl(iend,jend), dwl_obs(iend,jend))
  allocate(dep_cal(iend,jend), wl_cal(iend,jend), u_cal(iend,jend) , v_cal(iend,jend))

  ! ---- cal_slope & curvature ----  
  allocate(iwi_obs(iend,jend), iwj_obs(iend,jend), iw_obs(iend,jend))
  allocate(iwi_cal(iend,jend), iwj_cal(iend,jend), iw_cal(iend,jend))
  allocate(cbj(iend,jend)    , cbi(iend,jend)    , cb(iend,jend))
  allocate(ibi(iend,jend)    , ibj(iend,jend)    , ib(iend,jend))
  allocate(iei(iend,jend)    , iej(iend,jend)    , ie(iend,jend))

  ! ---- cal_tau, as, dz ----
  allocate(vel_cal(iend,jend), taue_i(iend,jend), taue_j(iend,jend), taue(iend,jend), tautau(iend,jend), Fr(iend,jend))
  allocate(qb_mpm(iend,jend), dz_mpm(iend,jend))
  allocate(Pe_s(iend,jend), Ms(iend,jend), dz_Ms(iend,jend))
  allocate(Pe_u(iend,jend), Mu(iend,jend), dz_Mu(iend,jend))
  allocate(Pe_d(iend,jend), Md(iend,jend), Df(iend,jend), dz_Md(iend,jend))
  allocate(lcl(iend,jend) , adv(iend,jend), dff(iend,jend), frc(iend,jend))
  allocate(Mx(iend,jend)  , My(iend,jend), Mxy(iend,jend))

  ! ---- cal and cmp dzdt & dhdx ----
  allocate(dz_obs(iend,jend)   , dz_amf(iend,jend))
  allocate(dif_dep(iend,jend)  , dif_dzs(iend,jend)  , dif_dzu(iend,jend), dif_dzd(iend,jend), dif_dzmpm(iend,jend))
  allocate(dhdx_obs(iend,jend) , dhdy_obs(iend,jend))
  allocate(dhdx_Mu(iend,jend)  , dhdy_Mu(iend,jend))
  allocate(dhdx_Ms(iend,jend)  , dhdx_Md(iend,jend))
  allocate(dif_dhdxs(iend,jend), dif_dhdys(iend,jend))
  allocate(dif_dhdxu(iend,jend), dif_dhdyu(iend,jend))
  
  allocate(ave_crs_dbl(jend), ave_crs_dwl(jend))
  
  ! ---- cal_nondim ----
  allocate(tau_nd(iend,jend), fr_nd(iend,jend), Mu_nd(iend,jend))  
 
  ! ---- cmp_each_item ---- 
  allocate(dudt(iend,jend), dudx(iend,jend), asu(iend,jend))

end subroutine set_arry



subroutine deallocate
  use mndr
  
  ! ---- read_schalar_vtk & read_vector_vtk ----
  deallocate(x, y)
  deallocate(dep_obs, wl_obs, dev_bl, bl, dwl_obs)
  deallocate(dep_cal, wl_cal, u_cal, v_cal)
    
  ! ---- cal_slope ----
  deallocate(iwi_obs, iwj_obs, iw_obs)
  deallocate(iwi_cal, iwj_cal, iw_cal)
  deallocate(cbj, cbi, cb)
  deallocate(iei, iej, ie)
  deallocate(ibi, ibj, ib)
  deallocate(Fr , vel_cal)    
  
  ! ---- cal_tau, pe, dz ----
  deallocate(taue_i, taue_j, taue, tautau)
  deallocate(qb_mpm, dz_mpm)
  deallocate(Pe_s, Ms, dz_Ms)
  deallocate(Pe_u, Mu, dz_Mu)
  deallocate(Pe_d, Md, Df, dz_Md)
  deallocate(lcl , adv, dff, frc)
  deallocate(Mx  , My, Mxy)

  ! ---- cal_nondim ----
  deallocate(tau_nd,fr_nd,Mu_nd)
  
  ! ---- cal and cmp dzdt & dhdx ----  
  deallocate(dz_obs   , dz_amf)
  deallocate(dif_dep  , dif_dzs  , dif_dzu, dif_dzd, dif_dzmpm)
  deallocate(dhdx_obs , dhdy_obs)
  deallocate(dhdx_Mu  , dhdy_Mu)
  deallocate(dhdx_Ms  , dhdx_Md)
  deallocate(dif_dhdxs, dif_dhdys)
  deallocate(dif_dhdxu, dif_dhdyu)

  deallocate(ave_crs_dbl, ave_crs_dwl)
  
 ! ---- cmp_each_item ----   
  deallocate(dudt, dudx, asu) 
  
end subroutine deallocate



subroutine set_parameter
  use mndr

  pi    = acos(-1.0d0)
  g     = 9.8
  intbl = 1
  s     = 1.65
  dm    = 0.00076  

end subroutine set_parameter



subroutine read_fln
  use mndr  
  
  call system('ls -1 ../00_input_vtk/'//trim(case)//'/vtk_schalar_gaussian/ > ./file.txt')

  open(100,file = 'file.txt')  

  n = 0
  do
     n = n + 1
     read(100, *, end = 998)
  enddo
! 998 ! write(*,*) 'Number of file =', n - 1
998  nend = n - 1
  close(100)
  
  allocate(fln(nend),fln2(nend))
  allocate(time1(nend))
  allocate(time2(nend))
  
  open(100,file = 'file.txt')  
  read(100,*)(fln(n), n = 1, nend)
  close(100)

  open(100,file = 'file.txt')  
  read(100,*)(time1(n), n = 1, nend)
  close(100)

  do n = 1, nend
     read(time1(n),'(f4.0)')time2(n)
  enddo
  
  do n = 1, nend
     read(fln(n),'(a)')fln2(n)
  enddo
  
end subroutine read_fln



subroutine read_arry
  use mndr  
  open(100,file = '../00_input_vtk/'//trim(case)//'/vtk_schalar_gaussian/'//trim(fln(1)))
  read(100,'(a)')(chr, i = 1, 4)
  read(100,*)chr, iend, jend
  close(100)
  allocate(bl2(iend, jend, nend), u_cal2(iend,jend,nend), mut(iend,jend,nend), ave_Mu(nend))
end subroutine read_arry



subroutine read_schalar_vtk(z, kend, knm, pss)
  use mndr
  character*100, intent(in) :: pss
  integer,intent(in)        :: kend, knm
  real*8,intent(out)        :: z(iend,jend)

  if(pss .eq. pss_st)then
     open(100,file = '../'//trim(pss)//'/'//trim(case)//'/vtk_schalar_gaussian/'//trim(fln(n)))
  elseif(pss .eq. pss_iric)then
     open(100,file = '../'//trim(pss)//'/'//trim(case)//'/vtk_schalar_gaussian/'//trim(fln(n)(1:4)//'.vtk'))
  endif
     
  read(100,'(a)')(chr, i = 1, 6)

  do j = 1, jend
     do i = 1, iend
        read(100,*)x(i,j), y(i,j)
     enddo
  enddo
  y(:,:) = y(:,:) - minval(y)
  
  read(100,'(a)')(chr, i = 1, 3)  
  do k = 1, kend
     read(100,*)   
     do j = 1, jend
        do i = 1, iend
           if(k.eq.knm)then
              read(100,*)z(i,j)
           elseif(k.ne.knm)then              
              read(100,*)
           endif
        end do
     end do
  end do  
  close(100)
  
end subroutine read_schalar_vtk



subroutine set_pss
  use mndr
  pss_st   = '00_input_vtk'
  pss_iric = '01_input_result_of_iric'  
end subroutine set_pss



subroutine read_roughness
  use mndr
  real*8 :: k
  allocate(nc(nend))
  
  open(100, file = '01_input_manning/'//trim(case)//'_Manning.txt', status = 'old')
  do n = 1, nend
     read(100,*)k, nc(n)
  enddo
  close(100)
  
end subroutine read_roughness



subroutine read_vector_vtk(z1, z2, kend, knm, pss)
  use mndr
  character*100, intent(in) :: pss
  integer,intent(in)        :: kend, knm
  real*8,intent(out)        :: z1(iend,jend),z2(iend,jend)

  open(100,file = '../'//trim(pss)//'/'//trim(case)//'/vtk_vector_gaussian/'//trim(fln(n)))
  read(100,'(a)')(chr, i = 1, 6)

  do j = 1, jend
     do i = 1, iend
        read(100,'()')
     enddo
  enddo
  y(:,:) = y(:,:) - minval(y)
  
  read(100,'(a)')(chr, i = 1, 3)  
  do k = 1, kend
     read(100,*)   
     do j = 1, jend
        do i = 1, iend
           if(k.eq.knm)then
              read(100,*)z1(i,j), z2(i,j)
           elseif(k.ne.knm)then              
              read(100,*)
           endif
        end do
     end do
  end do  
  close(100)
  
end subroutine read_vector_vtk



subroutine moving_average(z)
  use mndr
  real*8, intent(inout)  :: z(iend,jend)
  mend = 6
  do m = 1, mend
     do i = 1, iend
        do j = 1, jend
           if(j.eq.1 .and. i.ge.2 .and. i.le.iend-1)then
              z(i,j) = sum(z(i-1:i+1,j:j+1)) / 6
           elseif(j.eq.jend .and. i.ge.2 .and. i.le.iend-1)then
              z(i,j) = sum(z(i-1:i+1,j-1:j)) / 6
           elseif(i.eq.1 .and. j.ge.2 .and. j.le.jend-1)then
              z(i,j) = sum(z(i:i+1,j-1:j+1)) / 6              
           elseif(i.eq.iend .and. j.ge.2 .and. j.le.jend-1)then
              z(i,j) = sum(z(i-1:i,j-1:j+1)) / 6                            
           elseif(i.eq.1 .and. j.eq.1)then 
              z(i,j) = sum(z(i:i+1,j:j+1)) / 4
           elseif(i.eq.1 .and. j.eq.jend)then 
              z(i,j) = sum(z(i:i+1,j-1:j)) / 4
           elseif(i.eq.iend .and. j.eq.1)then              
              z(i,j) = sum(z(i-1:i,j:j+1)) / 4
           elseif(i.eq.iend .and. j.eq.jend)then
              z(i,j) = sum(z(i-1:i,j-1:j)) / 4
           else
              z(i,j) = sum(z(i-1:i+1,j-1:j+1)) / 9
           endif
        enddo
     enddo
  enddo
end subroutine moving_average



subroutine cal_critical_shear_stress
  use mndr
  real*8 :: ustrc
  
  if(dm*100 .ge. 0.303)then
     ustrc = sqrt(80.9*(dm*100)) / 100
     tauc  = ustrc**2 / (s * g * dm)
  else if(dm*100 .ge. 0.118 .and. dm*100 .lt. 0.303)then
     ustrc = sqrt(134.6 * (dm*100)**(31./32)) / 100
     tauc  = ustrc**2 / (s * g * dm)
  else if(dm*100 .ge. 0.0565 .and. dm*100 .lt. 0.118)then
     ustrc = sqrt(55.0 * (dm*100)) / 100
     tauc  = ustrc**2 / (s * g * dm)
  else if(dm*100 .ge. 0.0065 .and. dm*100 .lt. 0.0565)then
     ustrc = sqrt(8.41 * (dm*100)**(11./32)) / 100
     tauc  = ustrc**2 / (s * g * dm)
  else if(dm*100 .lt. 0.0065)then
     ustrc = sqrt(226 * (dm*100)) / 100
     tauc  = ustrc**2 / (s * g * dm)
  endif

end subroutine cal_critical_shear_stress



subroutine calculate_velocity(rough, dep, slpi, slpj, veli, velj)
  use mndr
  real*8              :: vel(iend,jend), slp(iend,jend)
  real*8, intent(out) :: veli(iend,jend), velj(iend,jend)
  real*8, intent(in)  :: rough(nend), dep(iend,jend), slpi(iend,jend), slpj(iend,jend)

  do i = 1, iend
     do j = 1, jend
        veli(i,j) = 1./ rough(n) * dep(i,j)**(2./3.) * abs(slpi(i,j))**(1./2.)
        velj(i,j) = 1./ rough(n) * dep(i,j)**(2./3.) * abs(slpj(i,j))**(1./2.)
        vel(i,j) = sqrt(veli(i,j)**2 + veli(i,j)**2)
        veli(i, j) = vel(i,j) * cos( atan(slpj(i, j) / slpi(i, j)) )
        velj(i, j) = vel(i,j) * sin( atan(slpj(i, j) / slpi(i, j)) )
     enddo
  enddo

  if(n.eq.1)then
     velo_obs = sum(vel(:,:)) / (iend*jend)
  endif
  
end subroutine calculate_velocity



subroutine calculate_slope(wl, slpi, slpj, slp)
  use mndr
  real*8,intent(in)  :: wl(iend,jend)
  real*8,intent(out) :: slpi(iend,jend), slpj(iend,jend), slp(iend,jend)

  do i = 1, iend
     do j = 1, jend
        if(i.eq.1)then
           slpi(i,j) = ( wl(i,j) - wl(i+1,j) ) / ( x(i,j) - x(i+1,j) )
        else if(i.eq.iend)then
           slpi(i,j) = ( wl(i,j) - wl(i-1,j) ) / ( x(i,j) - x(i-1,j) )
        else
           slpi(i,j) = ( wl(i-1,j) - wl(i+1,j) ) / ( x(i-1,j) - x(i+1,j) )
        endif
     enddo
  enddo

  do i = 1, iend
     do j = 1, jend
        if(j.eq.1)then
           slpj(i,j) = ( wl(i,j) - wl(i,j+1) ) / ( y(i,j) - y(i,j+1) )
        else if(j.eq.jend)then
           slpj(i,j) = ( wl(i,j) - wl(i,j-1) ) / ( y(i,j) - y(i,j-1) )
        else
           slpj(i,j) = ( wl(i,j-1) - wl(i,j+1) ) / ( y(i,j-1) - y(i,j+1) )
        endif
     enddo
  enddo

  do i = 1, iend
     do j = 1, jend
        slp(i,j) = sqrt(slpi(i,j)**2 + slpj(i,j)**2)
     enddo
  enddo
  
end subroutine calculate_slope



subroutine calculate_curvature(wl, curi, curj, cur)
  use mndr
  real*8,intent(in)  :: wl(iend,jend)
  real*8,intent(out) :: curi(iend,jend), curj(iend,jend), cur(iend,jend)

  do i = 1, iend
     do j = 1, jend
        if(i.eq.1)then
           curi(i,j) = ( (wl(i+2,j) - wl(i+1,j)) / (x(i+2,j) - x(i+1,j)) - (wl(i+1,j) - wl(i,j)) / (x(i+1,j) - x(i,j)) ) &
                * ( (x(i+1,j) - x(i,j)) / 2.0 )
        else if(i.eq.iend)then
           curi(i,j) = ( (wl(i,j) - wl(i-1,j)) / (x(i,j) - x(i-1,j)) - (wl(i-1,j) - wl(i-2,j)) / (x(i-1,j) - x(i-2,j)) ) &
                * ( (x(i,j) - x(i-2,j)) / 2.0 )
        else
           curi(i,j) = ( (wl(i+1,j) - wl(i,j)) / (x(i+1,j) - x(i,j)) - (wl(i,j) - wl(i-1,j)) / (x(i,j) - x(i-1,j)) ) &
                * ( (x(i+1,j) - x(i-1,j)) / 2.0 )
        endif
     enddo
  enddo
  
  do i = 1, iend
     do j = 1, jend
        if(j.eq.1)then
           curj(i,j) = ( (wl(i,j+2) - wl(i,j+1)) / (y(i,j+2) - y(i,j+1)) - (wl(i,j+1) - wl(i,j)) / (y(i,j+1) - y(i,j)) ) &
                * ( (y(i,j+1) - y(i,j)) / 2.0 )
        else if(j.eq.jend)then
           curj(i,j) = ( (wl(i,j) - wl(i,j-1)) / (y(i,j) - y(i,j-1)) - (wl(i,j-1) - wl(i,j-2)) / (y(i,j-1) - y(i,j-2)) ) &
                * ( (y(i,j) - y(i,j-2)) / 2.0 )
        else
           curj(i,j) = ( (wl(i,j+1) - wl(i,j)) / (y(i,j+1) - y(i,j)) - (wl(i,j) - wl(i,j-1)) / (y(i,j) - y(i,j-1)) ) &
                * ( (y(i,j+1) - y(i,j-1)) / 2.0 )
        endif
     enddo
  enddo
  
  do i = 1, iend
     do j = 1, jend
        cur(i,j) = sqrt(curi(i,j)**2 + curj(i,j)**2)
     enddo
  enddo
  
end subroutine calculate_curvature



subroutine calculate_energy_slope(rough, dep, veli, velj, slpi, slpj, slp)
  use mndr
  real*8              :: vel(iend,jend)
  real*8, intent(in)  :: rough(nend), dep(iend,jend), veli(iend,jend), velj(iend,jend)
  real*8, intent(out) :: slpi(iend,jend), slpj(iend,jend), slp(iend,jend)

  do i = 1, iend
     do j = 1, jend
        vel(i, j) = sqrt(veli(i,j)**2 + velj(i,j)**2)
        slp(i, j) = ( rough(n) * vel(i,j) / dep(i,j)**(2./3.) )**2
        slpi(i, j) = slp(i, j) * cos( atan(velj(i, j) / veli(i, j)) )
        slpj(i, j) = slp(i, j) * sin( atan(velj(i, j) / veli(i, j)) )
        if(isnan(slpi(i,j)) .or. isnan(slpi(i,j)) .or. isnan(slpi(i,j)))then
           slpi(i,j) = 0.0
           slpj(i,j) = 0.0
           slp(i,j)  = 0.0
        endif
     enddo
  enddo

  if(n.eq.1)then
     velo_cal = sum(vel(:,:)) / (iend*jend)
  endif
  
end subroutine calculate_energy_slope



subroutine calculate_shields_number(dep, slpi, slpj, taui, tauj, tau)
  use mndr
  real*8, intent(in)  :: dep(iend, jend), slpi(iend, jend), slpj(iend, jend)
  real*8, intent(out) :: taui(iend, jend), tauj(iend, jend), tau(iend, jend)
    
  do j = 1, jend
     do i = 1, iend
        taui(i, j) = ( dep(i, j) * slpi(i, j) / (s * dm) )
        tauj(i, j) = ( dep(i, j) * slpj(i, j) / (s * dm) )
        tau(i, j)  = ( sqrt( taui(i, j)**2 + tauj(i, j)**2 ) ) / tauc
     enddo
  enddo
  
end subroutine calculate_shields_number



subroutine calculate_Exner_equation_mpm(nb1, h, v, qbo, dzdt)
  use mndr
  real*8, intent(in)  :: nb1(nend), h(iend,jend), v(iend,jend)
  real*8, intent(out) :: qbo(iend,jend), dzdt(iend,jend)
  real*8              :: tau

  do i = 1, iend
     do j = 1,jend
        tau      = nb1(n)**2 * v(i,j)**2 / (s * dm * h(i,j)**(1./3.))
        qbo(i,j) = 8 * (tau - tauc)**(3./2.) * sqrt(s * g * dm**3)
        if((tau - tauc).lt.0)then
           qbo(i,j) = 0.0
        endif
     enddo
  enddo

  do i = 1, iend
     do j = 1,jend
        if(i.eq.1)then
           dzdt(i,j) = - 1. / (1.-0.4) * - ( ( qbo(i,j) - qbo(i+1,j) )   / (x(i,j) - x(i+1,j)) ) &
                * 10**3 * (time2(n+1)-time2(n)) * 60 * intbl
        elseif(i.eq.iend)then
           dzdt(i,j) = - 1. / (1.-0.4) * - ( ( qbo(i-1,j) - qbo(i,j) )   / (x(i-1,j) - x(i,j)) ) &
                * 10**3 * (time2(n+1)-time2(n)) * 60 * intbl
        else
           dzdt(i,j) = - 1. / (1.-0.4) * - ( ( qbo(i-1,j) - qbo(i+1,j) ) / (x(i-1,j) - x(i+1,j)) ) &
                * 10**3 * (time2(n+1)-time2(n)) * 60 * intbl
        endif
     enddo
  enddo
        
end subroutine calculate_Exner_equation_mpm



subroutine calculate_M_equation_mpm_s(nb1, h, slp1, slp2, v, Pec, Aso, dzdt)
  use mndr
  real*8,intent(in)  :: nb1(nend),h(iend,jend), slp1(iend,jend), slp2(iend,jend),v(iend,jend)
  real*8,intent(out) :: Pec(iend,jend), Aso(iend,jend), dzdt(iend,jend)
  real*8             :: Fro, tau

  do i = 1, iend
     do j = 1,jend
        tau = nb1(n)**2 * v(i,j)**2 / (s * dm * h(i,j)**(1./3.))
        Fro = v(i,j) / sqrt(g * h(i,j))
        Aso(i,j) = 4*(tau - tauc)**(1./2.)*sqrt(s * g * dm**3)*slp1(i,j) / (s*dm * (1-0.4) * ((1 + Fro)/2.))
        if((tau - tauc).lt.0)then
           Aso(i,j) = 0.0
        endif
        Aso(i,j) = Aso(i,j) * 10**3
        Pec(i,j) = slp1(i,j) / slp2(i,j)
        dzdt(i,j) = (-Aso(i,j) * (slp2(i,j) + slp1(i,j))) * (time2(n+1)-time2(n)) * 60 * intbl
     enddo
  enddo
  
end subroutine calculate_M_equation_mpm_s



subroutine calculate_M_equation_mpm_u(nb1, h, slp1, slp2, v, Pec, Aso, dzdt)
  use mndr
  real*8,intent(in)  :: nb1(nend),h(iend,jend), slp1(iend,jend), slp2(iend,jend),v(iend,jend)
  real*8,intent(out) :: Pec(iend,jend), Aso(iend,jend), dzdt(iend,jend)
  real*8             :: Fro, tau

  do i = 1, iend
     do j = 1,jend
        tau = nb1(n)**2 * v(i,j)**2 / (s * dm * h(i,j)**(1./3.))
        Fro = v(i,j) / sqrt(g * h(i,j))
        Aso(i,j) = 4.*(tau - tauc)**(1./2.)*sqrt(s * g * dm**3)*slp1(i,j) / (s*dm * (1-0.4) * (1 - 4./9.*Fro**2))
        if((tau - tauc).le.0)then
           Aso(i,j) = 0.0
        endif
        if(fro.eq.1.5)then
           write(*,*)'Fr=1.5'
        endif
        Aso(i,j) = Aso(i,j) * 10**3
        if(isnan(Aso(i,j)))then
           Aso(i,j) = 0.0
        endif
        Pec(i,j) = slp1(i,j) / slp2(i,j)
        dzdt(i,j) = (-Aso(i,j) * (slp2(i,j) + slp1(i,j))) * (time2(n+1)-time2(n)) * 60 * intbl
     enddo
  enddo
  
end subroutine calculate_M_equation_mpm_u



subroutine calculate_M_equation_mpm_df(nb1, h, slp1, slp2, cur, v, Aso, Dfo, Pec, dzdt)
  use mndr
  real*8,intent(in)  :: nb1(nend),h(iend,jend), slp1(iend,jend), slp2(iend,jend), cur(iend,jend), v(iend,jend)
  real*8,intent(out) :: Aso(iend,jend), Dfo(iend,jend), Pec(iend,jend), dzdt(iend,jend)
  real*8             :: Fro, tau, muc
  muc = 0.47
  
  do i = 1, iend
     do j = 1,jend
        tau = nb1(n)**2 * v(i,j)**2 / (s * dm * h(i,j)**(1./3.))
        Fro = v(i,j) / sqrt(g * h(i,j))
        Aso(i,j) = 4*(tau - tauc - tauc/muc*slp2(i,j))**(1./2.)*sqrt(s*g*dm**3)*slp1(i,j) / (s*dm * (1-0.4) * (1 - 4./9.*Fro**2))
        Dfo(i,j) = 12*sqrt(s*g*dm**3)*tauc/(muc*(1-0.4)) * (tau - tauc - tauc/muc*slp2(i,j))**(1./2.)
        Pec(i,j)  = (Aso(i,j)*dm) / Dfo(i,j)
        if((tau - tauc - tauc/muc*slp2(i,j)).lt.0)then
           Aso(i,j) = 0.0
           Dfo(i,j) = 0.0
           Pec(i,j) = 0.0
        endif
        dzdt(i,j) = (-Aso(i,j) * (slp2(i,j) + slp1(i,j)) + Dfo(i,j)*cur(i,j)) * 10**3 * (time2(n+1)-time2(n)) * 60 * intbl
     enddo
  enddo
  
end subroutine calculate_M_equation_mpm_df



subroutine calculate_M_equation_mpm_2d(nb1, h, slp1_i, slp1_j, slp2_i, slp2_j, u, v, dhdx, dhdy, Mox, Moy, Moxy, dzdt)
  use mndr
  real*8,intent(in)  :: nb1(nend), h(iend,jend), slp1_i(iend,jend), slp1_j(iend,jend), slp2_i(iend,jend), slp2_j(iend,jend)
  real*8,intent(in)  :: u(iend,jend), v(iend,jend)
  real*8,intent(out) :: Mox(iend,jend), Moy(iend,jend), Moxy(iend,jend), dzdt(iend,jend), dhdx(iend,jend), dhdy(iend,jend)
  real*8             :: Fro_i, Fro_j, tau_i, tau_j, M1, M2, M3, M4, ubyv, vbyu
  real*4             :: sig1, sig2
  
  do i = 1, iend
     do j = 1,jend
        tau_i = nb1(n)**2 * u(i,j)**2 / (s * dm * h(i,j)**(1./3.))
        tau_j = nb1(n)**2 * v(i,j)**2 / (s * dm * h(i,j)**(1./3.))
        Fro_i = u(i,j) / sqrt(g * h(i,j))
        Fro_j = v(i,j) / sqrt(g * h(i,j))

        ubyv = u(i,j)/v(i,j)
        vbyu = v(i,j)/u(i,j)

        if(h(i,j) .eq. 0)then
           Fro_i     = 0
           Fro_j     = 0
        endif
        
        if(u(i,j) .eq. 0 .or. v(i,j) .eq. 0)then
           ubyv = 0.0
           vbyu = 0.0
        endif
        
        sig1  = u(i,j)
        sig2  = v(i,j)

        
        M2    = 4 * sign(1.0, sig1) * (tau_i - tauc)**(1./2.) * sqrt(s*g*dm**3) * slp1_i(i,j) * (1 + 4./3. * Fro_j**2) &
             / (s * dm * (1-0.4) * (1 + 4./3.* (Fro_i**2 + Fro_j**2 + Fro_i**2 * Fro_j**2)) )

        M3    = 4 * sign(1.0, sig2) * (tau_j -tauc)**(1./2.) * sqrt(s*g*dm**3) * slp1_j(i,j) * (2./3. * vbyu * Fro_i**2) &
             / (s * dm * (1-0.4) * (1 + 4./3.* (Fro_i**2 + Fro_j**2 + Fro_i**2 * Fro_j**2)) )

        M4    = 4 * sign(1.0, sig2) * (tau_j - tauc)**(1./2.) * sqrt(s*g*dm**3) * slp1_j(i,j) * (1 + 4./3. * Fro_i**2) &
             / (s * dm * (1-0.4) * (1 + 4./3.* (Fro_i**2 + Fro_j**2 + Fro_i**2 * Fro_j**2)) )

        M1    = 4 * sign(1.0, sig1) * (tau_i - tauc)**(1./2.) * sqrt(s*g*dm**3) * slp1_i(i,j) * (2./3. * ubyv * Fro_j**2) & 
             / (s * dm * (1-0.4) * (1 + 4./3.* (Fro_i**2 + Fro_j**2 + Fro_i**2 * Fro_j**2)) )


        if((tau_i - tauc) .lt. 0 .or. h(i,j) .eq. 0.)then
           M2     = 0
           M1     = 0
        endif
        if((tau_j - tauc) .lt. 0 .or. h(i,j) .eq. 0.)then
           M3     = 0
           M4     = 0
        endif
        if(v(i,j) .eq. 0 .or. h(i,j) .eq. 0.)then
           M1     = 0
        endif
        if(u(i,j) .eq. 0 .or. h(i,j) .eq. 0.)then
           M3     = 0
        endif

        
        dhdx(i,j) = ((2./3. * ubyv * Fro_j**2) * (slp1_j(i,j) - slp2_j(i,j)) &
             + (1 + 4./3. * Fro_j**2) * (slp2_i(i,j) - slp1_i(i,j))) &
             / (1 + 4./3.* (Fro_i**2 + Fro_j**2 + Fro_i**2 * Fro_j**2))

        ! dhdx(i,j) = (slp2_i(i,j) - slp1_i(i,j)) / (1 + 4./3.* Fro_i**2 )
                
        dhdy(i,j) = ((2./3. * vbyu * Fro_i**2) * (slp1_i(i,j) - slp2_i(i,j)) &
             + (1 + 4./3. * Fro_i**2) * (slp2_j(i,j) - slp1_j(i,j))) &
             / (1 + 4./3.* (Fro_i**2 + Fro_j**2 + Fro_i**2 * Fro_j**2))

        ! dhdy(i,j) = (slp2_j(i,j) - slp1_j(i,j)) / (1 + 4./3.* Fro_j**2 )               
        
        Mox(i,j)  = (M2 - M3) * 10**3
        Moy(i,j)  = (M4 - M1) * 10**3
        
        Moxy(i,j) = sqrt( Mox(i,j)**2 + Moy(i,j)**2 )
        
        dzdt(i,j) = ( Mox(i,j) * (slp2_i(i,j) - slp1_i(i,j)) + Moy(i,j) * (slp2_j(i,j) - slp1_j(i,j)) ) &
             * (time2(n+1)-time2(n)) * 60 * intbl

     enddo
  enddo
  
end subroutine calculate_M_equation_mpm_2d




subroutine calculate_M_equation_amf(nb1, h, slp1, slp2, v, Pec, Aso, dzdt)
  use mndr
  real*8,intent(in)  :: nb1(nend),h(iend,jend), slp1(iend,jend), slp2(iend,jend),v(iend,jend)
  real*8,intent(out) :: Pec(iend,jend), Aso(iend,jend), dzdt(iend,jend)
  real*8             :: Fro, tau, ustr, ustrc

  ustrc = sqrt(0.034*s*g*dm)
  do i = 1, iend
     do j = 1,jend
        tau      = nb1(n)**2 * v(i,j)**2 / (s * dm * h(i,j)**(1./3.))
        ustr     = sqrt(tau * s * g * dm)
        Fro      = v(i,j) / sqrt(g * h(i,j))
        Aso(i,j)   = (119 * (1 - ustrc/ustr) * (-tauc+3*tau) * sqrt(s*g*dm**3) * slp1(i,j)) / (6*sqrt(tau)*s*dm*(1-Fro**2)*(1-0.4))
        if((tau - tauc).lt.0)then
           Aso(i,j) = 0.0
        endif
        Pec(i,j) = slp2(i,j) / slp1(i,j)
        ! dzdt(i,j) = (Aso(i,j) * (slp2(i,j) - slp1(i,j))) * 10**3 * 600  ! <--- st:10min
        dzdt(i,j) = (Aso(i,j) * (slp2(i,j) - slp1(i,j))) * 10**3 * 300    ! <--- st:5min
     enddo
  enddo
end subroutine calculate_M_equation_amf



subroutine filter(fro,aso)
  use mndr
  real*8,intent(in)    :: fro(iend,jend)
  real*8,intent(inout) :: Aso(iend,jend)
  count = 0
  do i = 1, iend
     do j = 1, jend
        if(fro(i,j).gt.1.0)then
           count = count + 1
           Aso(i,j) = 0.0
        endif
     enddo
  enddo  
end subroutine filter



subroutine cal_nondim_celerity(tau, fro, slp, M_in, tau_non, fr_non, M_non)
  use mndr
  real*8,intent(in)   :: tau(iend,jend), fro(iend,jend), slp(iend,jend), M_in(iend,jend)
  real*8,intent(out)  :: tau_non(iend,jend), fr_non(iend,jend), M_non(iend,jend)

  ! Aso(i,j) = 4 * (tau - tauc)**(1./2.) * sqrt(s * g * dm**3)*slp1(i,j) / (s*dm * (1-0.4) * (1 + 4./3.*Fro**2))
  
  do i = 1, iend
     do j = 1, jend
        if(tautau(i,j).lt.1.0)then
           tau_non(i,j) = 0.0
           fr_non(i,j)  = 0.0
           M_non(i,j)   = 0.0
        else
           tau_non(i,j) = (tau(i,j) - tauc)**(1./2.)
           fr_non(i,j)  = 1. / (1. + 4./3. * Fro(i,j)**2)
           M_non(i,j)   = tau_non(i,j) * fr_non(i,j) * slp(i,j)
        endif
     enddo
  enddo
    
end subroutine cal_nondim_celerity



subroutine calclate_and_compare_dzdt_dhdx
  use mndr
  ! ---- dzdt ----
  if(n.eq.nend)then
     dz_obs(:,:) = 0.0
  else
     dz_obs(:,:) = (bl2(:,:,n+intbl) - bl2(:,:,n)) * 10**3
  endif
  
  dif_dep(:,:)   = abs(( dep_obs(:,:) - dep_cal(:,:) ) / dep_obs(:,:) * 100)     
  dif_dzs (:,:)  = abs(dz_obs(:,:) - dz_Ms(:,:))  / 0.76 * 100
  dif_dzu(:,:)   = abs(dz_obs(:,:) - dz_Mu(:,:))  / 0.76 * 100
  dif_dzd(:,:)   = abs(dz_obs(:,:) - dz_Md(:,:))  / 0.76 * 100
  dif_dzmpm(:,:) = abs(dz_obs(:,:) - dz_mpm(:,:)) / 0.76 * 100

  do i = 1, iend
     do j = 1, jend
        dhdx_Ms(i,j) = (ibi(i,j) - iei(i,j)) / (1 - fr(i,j)**2)
        dhdx_Mu(i,j) = (ibi(i,j) - iei(i,j)) / (1 + 4./3. * fr(i,j)**2)
     enddo
  enddo

  ! ---- dhdx,dhdy ----
  do i = 1, iend
     do j = 1, jend
        if(i.eq.1)then
           dhdx_obs(i,j)      = (dep_obs(i+1,j) - dep_obs(i,j))   / (x(i+1,j) - x(i,j))
        elseif(i.eq.iend)then
           dhdx_obs(i,j)      = (dep_obs(i,j) - dep_obs(i-1,j))   / (x(i,j) - x(i-1,j))
        else
           dhdx_obs(i,j)      = (dep_obs(i+1,j) - dep_obs(i-1,j)) / (x(i+1,j) - x(i-1,j))
        endif
     enddo
  enddo

  do i = 1, iend
     do j = 1, jend
        if(j.eq.1)then
           dhdy_obs(i,j)      = (dep_obs(i,j+1) - dep_obs(i,j))   / (y(i,j+1) - y(i,j))
        elseif(j.eq.jend)then
           dhdy_obs(i,j)      = (dep_obs(i,j) - dep_obs(i,j-1))   / (y(i,j) - y(i,j-1))
        else
           dhdy_obs(i,j)      = (dep_obs(i,j+1) - dep_obs(i,j-1)) / (y(i,j+1) - y(i,j-1))
        endif
     enddo
  enddo

  dif_dhdxs(:,:) = abs((dhdx_Ms(:,:) - dhdx_obs(:,:)) / dhdx_obs(:,:)) * 100
  dif_dhdxu(:,:) = abs((dhdx_Mu(:,:) - dhdx_obs(:,:)) / dhdx_obs(:,:)) * 100  
  dif_dhdyu(:,:) = abs((dhdy_Mu(:,:) - dhdy_obs(:,:)) / dhdy_obs(:,:)) * 100

end subroutine calclate_and_compare_dzdt_dhdx



subroutine cmp_each_item
  use mndr
  do i = 1, iend
     do j = 1, jend
        if(i.eq.1)then
           dudt(i,j) = abs(u_cal2(i,j,n+intbl) - u_cal2(i,j,n)) / (600*intbl)
           dudx(i,j) = u_cal2(i,j,n) * (u_cal2(i+1,j,n) - u_cal2(i,j,n)) / (x(i+1,j) - x(i,j))
        elseif(i.eq.iend)then
           dudt(i,j) = abs(u_cal2(i,j,n+intbl) - u_cal2(i,j,n)) / (600*intbl)
           dudx(i,j) = u_cal2(i,j,n) * (u_cal2(i,j,n) - u_cal2(i-1,j,n)) / (x(i,j) - x(i-1,j))
        else
           dudt(i,j) = abs(u_cal2(i,j,n+intbl) - u_cal2(i,j,n)) / (600*intbl)
           dudx(i,j) = u_cal2(i,j,n) * (u_cal2(i+1,j,n) - u_cal2(i-1,j,n)) / (x(i+1,j) - x(i-1,j))
        endif
     enddo
  enddo
end subroutine cmp_each_item



subroutine output_schalar_vtk
  use mndr
  call system('mkdir ../03_output_vtk')
  call system('mkdir ../03_output_vtk/'//trim(case))
  call system('mkdir ../03_output_vtk/'//trim(case)//'/vtk_schalar_gaussian')    

  open(100,file='../03_output_vtk/'//trim(case)//'/vtk_schalar_gaussian/'//trim(fln(n)))  
  write(100,'(a)') '# vtk DataFile Version 3.1'
  write(100,'(a,i4.4,a)') 'out.vtk'
  write(100,'(a)')'ASCII'
  write(100,'(a)')'DATASET STRUCTURED_GRID'
  write(100,'(a,5x,i4,5x,i4,5x,i4)')'DIMENSIONS', iend, jend, 1
  write(100,'(a,5x,i8,5x,a)')'POINTS', iend * jend,'float'

  do j = 1, jend
     do i = 1, iend
        write(100,*)x(i,j), y(i,j), 0.d0
     enddo
  enddo

  write(100,*)
  write(100,'(a,5x,i8)')'POINT_DATA', iend * jend
  write(100,'(a,5x,i8)')'FIELD FieldData', 11
 
  call out_scalar(dev_bl,          'Deviation_of_bedlevel')
  call out_scalar(dep_cal,         'depth')    
  call out_scalar(wl_obs,          'Water_elevation')
  call out_scalar(Mu,              'celerity_unsteady')
  call out_scalar(tautau,          'tautau')
  call out_scalar(u_cal,           'velocity_i')
  call out_scalar(dz_Mu,         'dz_Mu')
  call out_scalar(dz_obs,         'dz_obs')
  call out_scalar(dz_mpm,         'dz_mpm')
  call out_scalar(dif_dzu,         'dif_dz_u')
  call out_scalar(dif_dzmpm,       'dif_d_mpm')
  
  ! call out_scalar(iei,               'Energy_slope')
  ! call out_scalar(tautau,            'effective_sheilds_number')
  ! call out_scalar(flux,      'flux')
  ! call out_scalar(dev_wl_cal(:,:)-0.057, 'dif_wl')
  ! call out_scalar((dev_wl_cal(:,:)-0.057)/0.057*100, 'dif_wl_persent')
  ! call out_scalar(abs(dz_mpm),       'dz_mpm')
  ! call out_scalar(dz_obs,            'dz_obs')
  ! call out_scalar(dif_dep,           'dif_dep')
  
  close(100)
end subroutine output_schalar_vtk



subroutine out_scalar(f, name)
  use mndr
  real*8,intent(in)       :: f(iend, jend)
  character(*),intent(in) :: name

  write(100,'(a,5x,i4,5x,i8,5x,a)')name, 1, iend*jend, 'double'
  do j = 1, jend
     do i = 1, iend
        write(100,*)f(i,j)
     end do
  end do

end subroutine out_scalar



subroutine out_vector_vtk
  use mndr

  intrvl_i = 1
  intrvl_j = 1
  
  ic = 0
  do i = 1, iend, intrvl_i
     ic = ic + 1
  enddo

  jc = 0
  do j = 1, jend, intrvl_j
     jc = jc + 1
  enddo

  ijc = ic * jc

  call system('mkdir ../03_output_vtk')
  call system('mkdir ../03_output_vtk/'//trim(case))
  call system('mkdir ../03_output_vtk/'//trim(case)//'/vtk_vector_gaussian')    
  
  open( 100,file='../03_output_vtk/'//trim(case)//'/vtk_vector_gaussian/'//trim(fln(n)))  
  write(100,'(a)') '# vtk DataFile Version 3.1'
  write(100,'(a,i4.4,a)') 'out.vtk'
  write(100,'(a)')'ASCII'
  write(100,'(a)')'DATASET STRUCTURED_GRID'
  write(100,'(a,5x,i4,5x,i4,5x,i4)')'DIMENSIONS', ic, jc, 1
  write(100,'(a,5x,i8,5x,a)')'POINTS', ijc,'float'

  do j = 1, jend, intrvl_j
     do i = 1, iend, intrvl_i
        write(100,*)x(i,j), y(i,j), 0.d0
     enddo
  enddo

  write(100,*)
  write(100,'(a,5x,i8)')'POINT_DATA', ijc
  write(100,'(a,5x,i8)')'FIELD FieldData', 3
 
!!$  call out_vector(iei,     iej,     'energy_slope')
  call out_vector(taue_i,  taue_j,  'tau')
  call out_vector(u_cal,   v_cal,   'velocity')
!!$  call out_vector(Mx,      My,      'Celerity')
  call out_vector(ibi,     ibj,     'bed_slope')
!!$  call out_vector(iwi_cal, iwj_cal, 'Water_slope_cal')
  
  close(100)

end subroutine out_vector_vtk



subroutine out_vector(fx, fy, name)
  use mndr
  real*8,intent(in)       :: fx(iend,jend), fy(iend,jend)
  character(*),intent(in) :: name

  write(100,'(a,5x,i4,5x,i8,5x,a)')name, 3, ijc, 'double'
  do j = 1, jend, intrvl_j
     do i = 1, iend, intrvl_i
        if( fx(i,j) .lt. 0 )then
           write(100,*)0, fy(i,j), 0
        else
           write(100,*)fx(i,j), fy(i,j), 0
        endif
     end do
  end do
end subroutine out_vector
