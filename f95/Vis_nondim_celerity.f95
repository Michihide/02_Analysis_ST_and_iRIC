! **** Draw 2d contour ****
subroutine drw_cntr
  use mndr
  use plplot

  if(n .eq. 1)then
     plparseopts_rc = plparseopts(PL_PARSE_FULL)
     if(plparseopts_rc .ne. 0) stop "plparseopt error"
     call plscol0( 0,255,255,255)
     call plscol0(15,  0,  0,  0)
     call plscmap0n(0)
     call plinit
  endif

  call pladv(0)
  call plschr(1.5d0, 2.0d0)  

  pdx   = 0.32d0
  pdy   = 0.20d0

  
! ! ---- cntr of D.B.L ---------------
!   barsname = 'D.B.L [m]'
!   axis_l   = 'Width(m)'
!   ! axis_t = ''
!   axis_t   = trim(fln2(n))
!   axis_b   = ''
!   ! axis_b   = 'Distance from upstream [m]'
!   clrnm    = 3

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!   ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!   zmin = 0.04d0 ; zmax = 0.06d0

!   pxmin = 0.05d0 ; pxmax = 0.45d0
!   pymin = 0.85d0 ; pymax = 0.95d0  
    
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call set_clr_cntr(clrnm, zmin, zmax)
!   call plot_cntr(1, iend, 1, jend, x, y, dev_bl)
!   call plcol0(9)
!   call pljoin(x(1,jend-2),y(1,jend-2),x(iend,jend-2),y(iend,jend-2))
!   ! do i = 1, iend
!   !    do j = 1, jend
!   !       if(dev_bl(i,j) .le. 0.060)then
!   !          call plot_point(1, x(i,j), y(i,j), 3, 1)
!   !       endif
!   !    enddo
!   ! enddo
!   ! call plcol0(9)
!   ! call pljoin(x(1,1), y(1,1), x(1,jend), y(1,jend))
!   ! call pljoin(x(450,1), y(450,1), x(450,jend), y(450,jend))  ! <- 1A
!   ! call pljoin(x(620,1), y(620,1), x(620,jend), y(620,jend))  ! <- 1B

!   ! call pljoin(x(70,1), y(70,1), x(70,jend), y(70,jend))
!   ! call pljoin(x(400,1), y(400,1), x(400,jend), y(400,jend))  ! <- 1C, 1D
  
!   ! call pljoin(x(1,jend-1), y(1,jend/2), x(iend,jend/2), y(iend,jend/2))
!   call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))  
! !------------------------------------------

  

! ! ---- cntr of D.B.L ---------------
!   barsname = 'Ms [mm/s] '
!   axis_l   = 'Width(m)'
!   axis_t   = ''
!   axis_b   = ''
!   clrnm    = 5

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!   ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!   zmin = 0.0d0 ; zmax = 1.0d0

!   pymin = pymin - pdy ; pymax = pymax - pdy
    
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call set_clr_cntr(clrnm, zmin, zmax)
!   call plot_cntr(1, iend, 1, jend, x, y, Ms)
!   call plcol0(9)
!   call pljoin(x(1,jend-2),y(1,jend-2),x(iend,jend-2),y(iend,jend-2))
!   call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
! !------------------------------------------


! ! ---- cntr of D.B.L ---------------
!   barsname = 'Mu [mm/s]'
!   axis_l   = 'Width(m)'
!   axis_t   = ''
!   axis_b   = ''
!   clrnm    = 5

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!   ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!   zmin = 0.0d0 ; zmax = 1.00d0

!   pymin = pymin - pdy ; pymax = pymax - pdy
    
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call set_clr_cntr(clrnm, zmin, zmax)
!   call plot_cntr(1, iend, 1, jend, x, y, Mu)
!   call plcol0(9)
!   call pljoin(x(1,jend-2),y(1,jend-2),x(iend,jend-2),y(iend,jend-2))
!   call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
! !------------------------------------------

! ! ---- cntr of D.B.L ---------------
!   barsname = 'M(unsteady)/u0 '
!   axis_l   = 'Width(m)'
!   axis_t   = ''
!   axis_b   = ''
!   clrnm    = 5

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!   ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!   zmin = minval(asu) - (maxval(asu)-minval(asu))*0.1d0
!   zmax = maxval(asu) + (maxval(asu)-minval(asu))*0.1d0  
!   zmin = 0.0d0 ; zmax = 0.20d0

!   pymin = pymin - pdy ; pymax = pymax - pdy
    
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call set_clr_cntr(clrnm, zmin, zmax)
!   call plot_cntr(1, iend, 1, jend, x, y, asu)
!   call plcol0(9)
!   call pljoin(x(1,jend-2),y(1,jend-2),x(iend,jend-2),y(iend,jend-2))
!   call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
! !------------------------------------------



! ! ---- cntr of D.B.L ---------------
!   barsname   = 'dif_dz(obs - steady)[%]'
!   axis_l   = 'Width(m)'
!   axis_t = ''
!   ! axis_t   = trim(fln2(n))
!   axis_b   = ''
!   ! axis_b   = 'Distance from upstream [m]'
!   clrnm    = 2

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!   ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!   zmin =  0.0d0
!   zmax =  1000.0d0

!   ! pxmin = 0.1d0 ; pxmax = 0.9d0
!   pymin = pymin - pdy ; pymax = pymax - pdy
    
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call set_clr_cntr(clrnm, zmin, zmax)
!   call plot_cntr(1, iend, 1, jend, x, y, dif_dz2)
!   call plcol0(9)
!   call pljoin(x(1,jend-2),y(1,jend-2),x(iend,jend-2),y(iend,jend-2))
!   call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
! !------------------------------------------


  
! ! ---- cntr of D.B.L ---------------
!   barsname = 'Mx [mm/s]'
!   axis_l   = 'Width(m)'
!   axis_t   = ''
!   axis_b   = ''
!   clrnm    = 5

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!   ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!   ! zmin = minval(Mx) - (maxval(Mx)-minval(Mx))*0.1d0
!   ! zmax = maxval(Mx) + (maxval(Mx)-minval(Mx))*0.1d0  
!   zmin = 0.0d0
!   zmax =  0.30d0

!   pymin = pymin - pdy ; pymax = pymax - pdy
    
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call set_clr_cntr(clrnm, zmin, zmax)
!   call plot_cntr(1, iend, 1, jend, x, y, Mx)
!   call plcol0(9)
!   call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
! !------------------------------------------



  
! ! ---- cntr of D.B.L ---------------
!   barsname = 'My [mm/s]'
!   axis_l   = 'Width(m)'
!   axis_t   = ''
!   axis_b   = ''
!   clrnm    = 3

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!   ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!   ! zmin = minval(My) - (maxval(My)-minval(My))*0.1d0
!   ! zmax = maxval(My) + (maxval(My)-minval(My))*0.1d0  
!   zmin = -0.10d0
!   zmax =  0.10d0

!   pymin = pymin - pdy ; pymax = pymax - pdy
    
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call set_clr_cntr(clrnm, zmin, zmax)
!   call plot_cntr(1, iend, 1, jend, x, y, My)
!   call plcol0(9)
!   call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
! !------------------------------------------



! ! ---- cntr of D.B.L ---------------
!   barsname = 'ie y'
!   axis_l   = 'Width(m)'
!   axis_t   = ''
!   axis_b   = ''
!   clrnm    = 3

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!   ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!   ! zmin = minval(My) - (maxval(My)-minval(My))*0.1d0
!   ! zmax = maxval(My) + (maxval(My)-minval(My))*0.1d0  
!   zmin = -0.0010d0
!   zmax =  0.0010d0

!   pymin = pymin - pdy ; pymax = pymax - pdy
    
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call set_clr_cntr(clrnm, zmin, zmax)
!   call plot_cntr(1, iend, 1, jend, x, y, iej)
!   call plcol0(9)
!   call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
! !------------------------------------------


  
! ! ---- cntr of D.B.L ---------------
!   barsname = 'Fr '
!   axis_l   = 'Width(m)'
!   axis_t   = ''
!   axis_b   = ''
!   clrnm    = 3

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!   ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!   zmin = 0.5d0 ; zmax = 1.5d0

!   pymin = pymin - pdy ; pymax = pymax - pdy
    
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call set_clr_cntr(clrnm, zmin, zmax)
!   call plot_cntr(1, iend, 1, jend, x, y, fr)
!   call plcol0(9)
!   call pljoin(x(1,jend-2),y(1,jend-2),x(iend,jend-2),y(iend,jend-2))
!   call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
! !------------------------------------------


  
!---- lng of dev_bl ---------------
  axis_l   = 'D.B.L [m]'
  axis_t   = trim(fln2(n))
  axis_b = ''

  ! pdy = 0.25d0  
  pxmin = 0.10d0 ; pxmax = 0.550d0
  pymin = 0.8d0  ; pymax = 0.95d0  

  xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
  xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
  ! ymin = minval(dev_bl(:,jend/2)) - (maxval(dev_bl(:,jend/2))-minval(dev_bl(:,jend/2)))*0.1d0
  ! ymax = maxval(dev_bl(:,jend/2)) + (maxval(dev_bl(:,jend/2))-minval(dev_bl(:,jend/2)))*0.1d0  
  ymin = 0.04d0 ; ymax = 0.06d0
  
  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
  call plot_line(iend, x(:,jend/2), dev_bl(:,jend-2), 1)
  call plcol0(7)
  call pljoin(0.0d0, 0.048d0, 10.0d0, 0.048d0)

  do i = 1, iend
     if(Mu_nd(i,jend-2) .eq. 0)then
        call plot_point(1, x(i,jend-2), dev_bl(i,jend-2), 3, 3)
     endif
  enddo
!------------------------------------------  

  

! !---- lng of M/v ---------------
!   ! axis_l   = 'Ms/sqrt(gh) [R]'
!   axis_l   = 'Mu/u [R] & Mu/sqrt(gh) [B]'
!   axis_t = ''
!   axis_b = ''

!   pdy = 0.22d0  
!   pymin = pymin - pdy ; pymax = pymax - pdy

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
  
!   ! ymin = minval(ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2))) - (maxval(ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2))) &
!   !      -minval(ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2))))*0.1d0
!   ! ymax = maxval(ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2))) + (maxval(ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2))) &
!   !      -minval(ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2))))*0.1d0
  
!   ! ymin = 0.0d0 ; ymax = 0.20d0
  
!   ! call set_window_first_axis(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   ! call plot_line(iend, x(:,jend/2), ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2)), 1)
!   ! ! call plcol0(7)
!   ! ! call pljoin(0.0d0, 1.0d0, 10.0d0, 0.0d0)

!   ! axis_l   = 'Mu/sqrt(gh) [B]'
  
!   ! ymin = minval(mu(:,jend-2)/10**(3)/sqrt(g*dep_cal(:,jend-2))) - (maxval(mu(:,jend-2)/10**(3)/sqrt(g*dep_cal(:,jend-2))) &
!   !      -minval(mu(:,jend-2)/10**(3)/sqrt(g*dep_cal(:,jend-2))))*0.1d0
!   ! ymax = maxval(mu(:,jend-2)/10**(3)/sqrt(g*dep_cal(:,jend-2))) + (maxval(mu(:,jend-2)/10**(3)/sqrt(g*dep_cal(:,jend-2))) &
!   !      -minval(mu(:,jend-2)/10**(3)/sqrt(g*dep_cal(:,jend-2))))*0.1d0
  
!   ymin = 0.0d0 ; ymax = 0.3d0
  
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call plot_line(iend, x(:,jend/2), (mu(:,jend-2)/10**(3)/sqrt(g*dep_cal(:,jend-2)))*100, 9)
!   call plot_line(iend, x(:,jend/2), (mu(:,jend-2)/10**(3)/u_cal(:,jend-2))*100, 1)
!   call plcol0(7)
!   call pljoin(0.0d0, 0.045d0, 10.0d0, 0.045d0)

!   do i = 1, iend
!      if(Mu_nd(i,jend-2) .eq. 0)then
!         call plot_point(1, x(i,jend-2), (mu(i,jend-2)/10**(3)/u_cal(i,jend-2))*100, 3, 3)
!         call plot_point(1, x(i,jend-2), (mu(i,jend-2)/10**(3)/sqrt(g*dep_cal(i,jend-2)))*100, 3, 3)
!      endif
!   enddo

! !------------------------------------------  

  

! !---- lng of dif_dz ---------------
!   axis_l   = 'dif_dz(s)=[R] & (u)=[B]'
!   Axis_t = ''
!   axis_b = ''

!   pymin = pymin - pdy ; pymax = pymax - pdy

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ! ymin = minval(ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2))) - (maxval(ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2))) &
!   !      -minval(ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2))))*0.1d0
!   ! ymax = maxval(ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2))) + (maxval(ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2))) &
!   !      -minval(ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2))))*0.1d0
!   ymin = 0.0d0 ; ymax = 1000.0d0
  
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call plot_line(iend, x(:,jend/2), dif_dz2(:,jend-2), 1)
!   call plot_line(iend, x(:,jend/2), dif_dz(:,jend-2), 9)

    
!   call plcol0(7)
!   call pljoin(0.0d0, 200.0d0, 10.0d0, 200.0d0)

!   do i = 1, iend
!      if(Mu_nd(i,jend-2) .eq. 0)then
!         call plot_point(1, x(i,jend-2), dif_dz2(i,jend-2), 3, 3)
!         call plot_point(1, x(i,jend-2), dif_dz(i,jend-2), 3, 3)
!      endif
!   enddo
! !------------------------------------------  



!---- lng of dif_dz ---------------
  axis_l  = 'dif_dz (u)=[R], (s)=[B]'
  Axis_t  = ''
  axis_b  = 'fr'

  pxmin = 0.10d0 ; pxmax = 0.30d0
  pymin = 0.45d0  ; pymax = 0.75d0  

  xmin = minval(fr) - (maxval(fr)-minval(fr))*0.02d0
  xmax = maxval(fr) + (maxval(fr)-minval(fr))*0.02d0
  ymin = minval(dif_dz) - (maxval(dif_dz) - minval(dif_dz))*0.1d0
  ymax = maxval(dif_dz) + (maxval(dif_dz) - minval(dif_dz))*0.1d0
  xmin = 0.0d0 ; xmax = 2.0d0
  ymin = 0.0d0 ; ymax = 1000.0d0
  
  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
  call plot_point(iend*jend, fr, dif_dz,  1, 1)
  call plot_point(iend*jend, fr, dif_dz2, 9, 1)
  call plcol0(7)
  call pljoin(0.8d0, 0.0d0, 0.8d0, 1000.0d0)  
!------------------------------------------  



!---- lng of dif_dz ---------------
  axis_l  = ''
  Axis_t  = ''
  axis_b  = 'Mu/sqrt(gh)'

  pxmin = 0.350d0 ; pxmax = 0.55d0


  xmin = minval(asu) - (maxval(asu)-minval(asu))*0.02d0
  xmax = maxval(asu) + (maxval(asu)-minval(asu))*0.02d0
  ymin = minval(dif_dz) - (maxval(dif_dz) - minval(dif_dz))*0.1d0
  ymax = maxval(dif_dz) + (maxval(dif_dz) - minval(dif_dz))*0.1d0
  xmin = 0.0d0 ; xmax = 0.2d0
  ymin = 0.0d0 ; ymax = 1000.0d0
  
  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
  call plot_point(iend*jend, asu, dif_dz,  1, 1)
  call plot_point(iend*jend, asu, dif_dz2, 9, 1)
  call plcol0(7)
  call pljoin(0.05d0, 0.0d0, 0.05d0, 1000.0d0)    
!------------------------------------------  


  
! !---- lng of dif_dz ---------------
!   axis_l   = 'Fr'
!   Axis_t = ''
!   axis_b = ''

!   pymin = pymin - pdy ; pymax = pymax - pdy

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ! ymin = minval(ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2))) - (maxval(ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2))) &
!   !      -minval(ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2))))*0.1d0
!   ! ymax = maxval(ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2))) + (maxval(ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2))) &
!   !      -minval(ms(:,jend-2)/sqrt(g*dep_cal(:,jend-2))))*0.1d0
!   ymin = 0.5d0 ; ymax = 1.5d0
  
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call plot_line(iend, x(:,jend/2), fr(:,jend-2), 1)

!   do i = 1, iend
!      if(Mu_nd(i,jend-2) .eq. 0)then
!         call plot_point(1, x(i,jend-2), fr(i,jend-2), 3, 3) 
!      endif
!   enddo
  
!   call plcol0(7)
!   call pljoin(0.0d0, 0.80d0, 10.0d0, 0.80d0)
! !------------------------------------------  


  
! !---- lng of M-M0 ---------------
!   axis_l   = 'Mu_nd-Mu_nd0'
!   axis_t = ''
!   axis_b = ''

!   pdy = 0.18d0  
!   pymin = pymin - pdy ; pymax = pymax - pdy

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(mu_nd(:,jend-2)-mu_nd0) - (maxval(mu_nd(:,jend-2)-mu_nd0)-minval(mu_nd(:,jend-2)-mu_nd0))*0.1d0
!   ymax = maxval(mu_nd(:,jend-2)-mu_nd0) + (maxval(mu_nd(:,jend-2)-mu_nd0)-minval(mu_nd(:,jend-2)-mu_nd0))*0.1d0  
!   ! ymin = 0.0d0 ; ymax = 0.20d0
  
!   call set_window_first_axis(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call plot_line(iend, x(:,jend/2), mu_nd(:,jend-2)-mu_nd0, 1)
!   call plcol0(7)
!   call pljoin(0.0d0, 0.0d0, 10.0d0, 0.0d0)

!   axis_l   = 'Mu_nd'
!   ymin = minval(mu_nd(:,jend-2)) - (maxval(mu_nd(:,jend-2))-minval(mu_nd(:,jend-2)))*0.1d0
!   ymax = maxval(mu_nd(:,jend-2)) + (maxval(mu_nd(:,jend-2))-minval(mu_nd(:,jend-2)))*0.1d0    
!   call set_window_second_axis(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call plot_line(iend, x(:,jend/2), mu_nd(:,jend-2), 9)

!   do i = 1, iend
!      if(Mu_nd(i,jend-2) .eq. 0)then
!         call plot_point(1, x(i,jend-2), mu_nd(i,jend-2), 3, 3)
!      endif
!   enddo
! !------------------------------------------  

  

! !---- lng of iei-iei0 ---------------
!   axis_l = 'iei-iei0'
!   axis_t = ''
!   axis_b = ''

!   pymin = pymin - pdy ; pymax = pymax - pdy

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(iei(:,jend-2)-iei0) - (maxval(iei(:,jend-2)-iei0)-minval(iei(:,jend-2)-iei0))*0.1d0
!   ymax = maxval(iei(:,jend-2)-iei0) + (maxval(iei(:,jend-2)-iei0)-minval(iei(:,jend-2)-iei0))*0.1d0  
!   ! ymin =  0.0d0 ; ymax =  1000.0d0
  
!   call set_window_first_axis(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call plot_line(iend, x(:,jend/2), iei(:,jend-2)-iei0, 1)
!   call plcol0(7)
!   call pljoin(0.0d0, 0.0d0, 10.0d0, 0.0d0)

!   axis_l = 'iei'
!   ymin = minval(iei(:,jend-2)) - (maxval(iei(:,jend-2))-minval(iei(:,jend-2)))*0.1d0
!   ymax = maxval(iei(:,jend-2)) + (maxval(iei(:,jend-2))-minval(iei(:,jend-2)))*0.1d0    
!   call set_window_second_axis(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call plot_line(iend, x(:,jend/2), iei(:,jend-2), 9)

!   do i = 1, iend
!      if(Mu_nd(i,jend-2) .eq. 0)then
!         call plot_point(1, x(i,jend-2), iei(i,jend-2), 3, 3)
!      endif
!   enddo

! !------------------------------------------  



! !---- lng of tau_nd-tau_nd0 ---------------
!   axis_l = 'tau_nd-tau_nd0'
!   axis_t = ''
!   axis_b = ''

!   pymin = pymin - pdy ; pymax = pymax - pdy

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(tau_nd(:,jend-2)-tau_nd0) - (maxval(tau_nd(:,jend-2)-tau_nd0)-minval(tau_nd(:,jend-2)-tau_nd0))*0.1d0
!   ymax = maxval(tau_nd(:,jend-2)-tau_nd0) + (maxval(tau_nd(:,jend-2)-tau_nd0)-minval(tau_nd(:,jend-2)-tau_nd0))*0.1d0  
!   ! ymin =  0.0d0 ; ymax =  1000.0d0
  
!   call set_window_first_axis(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call plot_line(iend, x(:,jend/2), tau_nd(:,jend-2)-tau_nd0, 1)
!   call plcol0(7)
!   call pljoin(0.0d0, 0.0d0, 10.0d0, 0.0d0)

!   axis_l = 'tau_nd'
!   ymin = minval(tau_nd(:,jend-2)) - (maxval(tau_nd(:,jend-2))-minval(tau_nd(:,jend-2)))*0.1d0
!   ymax = maxval(tau_nd(:,jend-2)) + (maxval(tau_nd(:,jend-2))-minval(tau_nd(:,jend-2)))*0.1d0    
!   call set_window_second_axis(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call plot_line(iend, x(:,jend/2), tau_nd(:,jend-2), 9)

!   do i = 1, iend
!      if(Mu_nd(i,jend-2) .eq. 0)then
!         call plot_point(1, x(i,jend-2), tau_nd(i,jend-2), 3, 3)
!      endif
!   enddo

! !------------------------------------------  



! !---- lng of tau_nd-tau_nd0 ---------------
!   axis_l = 'fr_nd-fr_nd0'
!   axis_t = ''
!   axis_b = ''

!   pymin = pymin - pdy ; pymax = pymax - pdy

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(fr_nd(:,jend-2)-fr_nd0) - (maxval(fr_nd(:,jend-2)-fr_nd0)-minval(fr_nd(:,jend-2)-fr_nd0))*0.1d0
!   ymax = maxval(fr_nd(:,jend-2)-fr_nd0) + (maxval(fr_nd(:,jend-2)-fr_nd0)-minval(fr_nd(:,jend-2)-fr_nd0))*0.1d0  
!   ! ymin =  0.0d0 ; ymax =  1000.0d0
  
!   call set_window_first_axis(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call plot_line(iend, x(:,jend/2), fr_nd(:,jend-2)-fr_nd0, 9)
!   call plcol0(7)
!   call pljoin(0.0d0, 0.0d0, 10.0d0, 0.0d0)

!   do i = 1, iend
!      if(Mu_nd(i,jend-2) .eq. 0)then
!         call plot_point(1, x(i,jend-2), fr_nd(i,jend-2), 3, 3)
!      endif
!   enddo


!   axis_l = 'fr'
!   ! ymin = minval(fr(:,jend-2)) - (maxval(fr(:,jend-2))-minval(fr(:,jend-2)))*0.1d0
!   ! ymax = maxval(fr(:,jend-2)) + (maxval(fr(:,jend-2))-minval(fr(:,jend-2)))*0.1d0  
!   ymin = 0.5d0 ; ymax =1.5d0
  
!   call set_window_second_axis(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call plot_line(iend, x(:,jend/2), fr(:,jend-2), 1)
!   call plcol0(7)
!   call pljoin(0.0d0, 0.75d0, 10.0d0, 0.75d0)
  
!   do i = 1, iend
!      if(Mu_nd(i,jend-2) .eq. 0)then
!         call plot_point(1, x(i,jend-2), fr(i,jend-2), 3, 3)
!      endif
!   enddo
! !------------------------------------------  


  
! !---- cross of dev_bl ---------------
!   axis_l   = 'cross_ave_dev_bl [m]'
!   axis_t = ''
!   axis_b = ''

!   pxmin = 0.10d0      ; pxmax = 0.40d0
!   pymin = 0.65d0 ; pymax = 0.9d0  
!   ! pymin = pymin - pdy ; pymax = pymax - pdy

!   xmin = minval(y) - (maxval(y)-minval(y))*0.02d0
!   xmax = maxval(y) + (maxval(y)-minval(y))*0.02d0
!   ! ymin = minval(ave_crs_dev_bl) - (maxval(ave_crs_dev_bl)-minval(ave_crs_dev_bl))*0.1d0
!   ! ymax = maxval(ave_crs_dev_bl) + (maxval(ave_crs_dev_bl)-minval(ave_crs_dev_bl))*0.1d0  
!   ymin = 0.043d0 ; ymax = 0.058d0
  
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call plot_line(jend, y(1,:), ave_crs_dev_bl(:), 9)
! !------------------------------------------  


  
! ! ---- Histgram of dzdt(obs) - dzdt(mpm) ------------
!   axis_l   = 'number of points'
!   axis_t   = '' 
!   axis_b   = 'M [mm/s]'

!   pxmin = 0.65d0  ; pxmax = 0.85d0
!   ! pymin = 0.75d0  ; pymax = 0.95d0
!   pymin = pymin - pdy ; pymax = pymax - pdy

!   xmin =  0.0d0  ; xmax = 500.0d0
!   ! ymin =  0.0d0  ; ymax = 1000.0d0
!   xmin = minval(as_mpm) - (maxval(as_mpm)-minval(as_mpm))*0.1d0
!   xmax = maxval(as_mpm) + (maxval(as_mpm)-minval(as_mpm))*0.1d0  

!   wdth= 20.0
  
!   pend = (xmax - xmin)/wdth + 2
!   print*, 'Number of section:', pend
!   allocate(hst_x(pend),hst_y(pend))

!   do p = 1, pend
!      con = 0
!      do j = 1, jend
!         do i = 1, iend
!            if(((p - 1)*wdth + xmin < as_mpm(i,j)) .and. (as_mpm(i,j) <= p*wdth + xmin))then
!               con = con + 1
!            endif
!         enddo
!      enddo
!      hst_x(p) = (p - 1)*wdth + xmin
!      hst_y(p) = con / (iend*jend) * 100
!   enddo
!   call set_window(xmin - wdth/2., xmax+wdth/2, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax) 
!   call plot_hstgrm(pend, hst_x, hst_y, wdth, 3)
!   deallocate(hst_x, hst_y)
! ! ------------------------------------------


  
! ! ---- cntr of D.B.L ---------------
!   barsname = 'dz(cal) [mm]'
!   axis_l   = 'Width(m)'
!   axis_t   = ''
!   axis_b   = ''
!   clrnm    = 5

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!   ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!   zmin = 0.0d0 ; zmax = 2000.0d0

!   pymin = pymin - pdy ; pymax = pymax - pdy
    
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call set_clr_cntr(clrnm, zmin, zmax)
!   call plot_cntr(1, iend, 1, jend, x, y, abs(dz_mpm / 0.76 * 100))
!   call plcol0(9)
!   call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
! !------------------------------------------

! ! ---- cntr of D.B.L ---------------
!   barsname = 'dif_dz [%]'
!   axis_l   = 'Width(m)'
!   axis_t   = ''
!   axis_b   = ''
!   clrnm    = 2

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!   ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!   zmin = 0.0d0 ; zmax = 500.0d0

!   pymin = pymin - pdy ; pymax = pymax - pdy
    
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call set_clr_cntr(clrnm, zmin, zmax)
!   call plot_cntr(1, iend, 1, jend, x, y, dif_dz)
!   call plcol0(9)
!   call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
! !------------------------------------------


! ! ---- cntr of D.B.L ---------------
!   barsname = 'dhdx(obs)'
!   axis_l   = 'Width(m)'
!   axis_t   = ''
!   axis_b   = ''
!   clrnm    = 3

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!   ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!   zmin = -0.01d0 ; zmax = 0.01d0

!   pymin = pymin - pdy ; pymax = pymax - pdy
    
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call set_clr_cntr(clrnm, zmin, zmax)
!   call plot_cntr(1, iend, 1, jend, x, y, dhdx)
!   call plcol0(9)
!   call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
! !------------------------------------------

! ! ---- cntr of D.B.L ---------------
!   barsname = 'dhdx(Unst)'
!   axis_l   = 'Width(m)'
!   axis_t   = ''
!   axis_b   = ''
!   clrnm    = 5

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!   ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!   ! zmin = -0.01d0 ; zmax = 0.01d0
!   zmin = 0.0d0 ; zmax = 1000.0d0

!   pymin = pymin - pdy ; pymax = pymax - pdy
    
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call set_clr_cntr(clrnm, zmin, zmax)
!   ! call plot_cntr(1, iend, 1, jend, x, y, dhdx_mpm)
!   call plot_cntr(1, iend, 1, jend, x, y, difdh_mpm)
!   call plcol0(9)
!   call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
! !------------------------------------------

! ! ---- cntr of D.B.L ---------------
!   barsname = 'dhdx(Nonuni)'
!   axis_l   = 'Width(m)'
!   axis_t   = ''
!   axis_b   = ''
!   clrnm    = 5

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!   ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!   ! zmin = -0.01d0 ; zmax = 0.01d0
!   zmin = 0.0d0 ; zmax = 1000.0d0

!   pymin = pymin - pdy ; pymax = pymax - pdy
    
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call set_clr_cntr(clrnm, zmin, zmax)
!   ! call plot_cntr(1, iend, 1, jend, x, y, dhdx_mpm2)
!   call plot_cntr(1, iend, 1, jend, x, y, difdh_mpm2)  
!   call plcol0(9)
!   call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
! !------------------------------------------

! ! ---- cntr of celerity ---------------  
!   barsname = 'Fr'
!   axis_l = 'Width(m)' ; axis_t = '' ; axis_b = '' ;  clrnm    = 3

!   xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!   xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!   ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!   ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!   zmin = 0.0d0
!   zmax = 2.0d0
  
!   pdy = 0.30d0
!   pymin = pymin - pdy
!   pymax = pymax - pdy
    
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call set_clr_cntr(clrnm, zmin, zmax)
!   call plot_cntr(1, iend, 1, jend, x, y, Fr)
!   call plcol0(0)
!   call plot_point(1, x1(n), y1(n), 15, 17)
!   call plot_point(1, x2(n), y2(n), 15, 17)
!   call plot_point(1, x3(n), y3(n), 15, 17)
!   call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
! !------------------------------------------

! !---- time change ---------------
!   axis_l   = 'As [mm/s]'
!   axis_t   = ''
!   axis_b   = ''
!   pxmin = 0.1d0 ; pxmax = 0.4d0
!   pymin = 0.8d0 ; pymax = 0.95d0
  
!   xmin = minval(time2(7:nend)) - (maxval(time2(7:nend))-minval(time2(7:nend)))*0.02d0
!   xmax = maxval(time2(7:nend)) + (maxval(time2(7:nend))-minval(time2(7:nend)))*0.02d0
!   ymin = 0.0d0
!   ymax = 0.6d0
  
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call plot_point(n-6, time2(7:n), as_cal1(7:n), 9, 21)
!   call plot_point(n-6, time2(7:n), as_cal4(7:n), 3, 0)
!   write(*,*)as_cal4(n)
!   if(n.ge.7)then
!      call plot_point(n-6, time2(7:n), as_obs1(7:n), 1, 18)
!   endif
! !------------------------------------------  
! !---- time change ---------------
!   axis_l   = 'As [mm/s]'
!   axis_t   = ''
!   axis_b   = ''
!   pdy = 0.2d0
!   pymin = pymin - pdy ; pymax = pymax - pdy

!   ! xmin = minval(time2(7:nend)) - (maxval(time2(7:nend))-minval(time2(7:nend)))*0.02d0
!   ! xmax = maxval(time2(7:nend)) + (maxval(time2(7:nend))-minval(time2(7:nend)))*0.02d0
!   ! ymin = 0.0d0
!   ! ymax = 0.8d0
  
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call plot_point(n-6, time2(7:n), as_cal2(7:n), 9, 21)
!   call plot_point(n-6, time2(7:n), as_cal5(7:n), 3, 0)
!   if(n.ge.7)then
!      call plot_point(n-6, time2(7:n), as_obs2(7:n), 1, 18)
!   endif
! !------------------------------------------  

  ! !---- time change ---------------
!   axis_l   = 'As [mm/s]'
!   axis_b   = 'Time [min]'

!   ! pxmin = 0.65d0 ; pxmax = 0.95d0
!   pymin = pymin - pdy ; pymax = pymax - pdy
!   ! pxmin = 0.7d0 ; pxmax = 0.9d0
  
!   ! xmin = minval(time2(7:nend)) - (maxval(time2(7:nend))-minval(time2(7:nend)))*0.02d0
!   ! xmax = maxval(time2(7:nend)) + (maxval(time2(7:nend))-minval(time2(7:nend)))*0.02d0
!   ! ymin = 0.0d0
!   ! ymax = 0.6d0
  
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call plot_point(n-6, time2(7:n), as_cal3(7:n), 9, 21)
!   call plot_point(n-6, time2(7:n), as_cal6(7:n), 3, 0)
!   if(n.ge.7)then
!      call plot_point(n-6, time2(7:n), As_obs3(7:n), 1, 18)
!   endif
! !------------------------------------------  

  
! !---- time change ---------------
!   axis_l   = 'As(cal)/As(obs)'
!   axis_t   = ''

!   pxmin = 0.6d0 ; pxmax = 0.95d0
!   pymin = 0.75d0 ; pymax = 0.95d0

!   xmin = minval(time2(1:nend)) - (maxval(time2(1:nend))-minval(time2(1:nend)))*0.02d0
!   xmax = maxval(time2(1:nend)) + (maxval(time2(1:nend))-minval(time2(1:nend)))*0.02d0
!   ymin = 0.0d0
!   ymax = 500.0d0
  
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call plot_point(n, time2(1:n), abs(as_cal1(1:n)-As_obs1(1:n))/As_obs1(1:n)*100, 9, 21)
! !------------------------------------------  
! !---- time change ---------------
!   axis_l   = 'As(cal)/As(obs)'
!   axis_t   = ''
!   pdy = 0.30d0
!   pymin = pymin - pdy ; pymax = pymax - pdy

!   ! xmin = minval(time2(7:nend)) - (maxval(time2(7:nend))-minval(time2(7:nend)))*0.02d0
!   ! xmax = maxval(time2(7:nend)) + (maxval(time2(7:nend))-minval(time2(7:nend)))*0.02d0
!   ymin = 0.0d0
!   ymax = 500.0d0
  
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call plot_point(n, time2(1:n), abs(as_cal2(1:n)-As_obs2(1:n))/As_obs2(1:n)*100, 9, 21)
! !------------------------------------------  
! !---- time change ---------------
!   axis_l   = 'As(cal)/As(obs)'
!   axis_t   = ''

!   ! pxmin = 0.65d0 ; pxmax = 0.95d0
!   pymin = pymin - pdy ; pymax = pymax - pdy
!   ! pxmin = 0.7d0 ; pxmax = 0.9d0
  
!   ! xmin = minval(time2(7:nend)) - (maxval(time2(7:nend))-minval(time2(7:nend)))*0.02d0
!   ! xmax = maxval(time2(7:nend)) + (maxval(time2(7:nend))-minval(time2(7:nend)))*0.02d0
!   ymin = 0.0d0
!   ymax = 500.0d0
  
!   call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!   call plot_point(n, time2(1:n), abs(as_cal3(1:n)-As_obs3(1:n))/As_obs3(1:n)*100, 9, 21)
! !------------------------------------------  
  

! ! ---- cal for cntr of Pe,ib,As ---------------
!   icount(:,:) = 0.0  ; icount2(:,:) = 0.0 ; z4(:,:) = 0.0
!   x4_min  =  0.0d0   ; x4_max  = 0.05d0
!   y4_min  = -0.05d0  ; y4_max  = 0.05d0
!   x4(1,:) =  x4_min  ; y4(:,1) = y4_min
!   do m = 2, mend
!      do l = 2, lend
!         x4(m,:) = x4(m-1,:) + (x4_max-x4_min)/(mend-1)
!         y4(:,l) = y4(:,l-1) + (y4_max-y4_min)/(lend-1)
!      enddo
!   enddo  
!   do m = 1, mend-1
!      do l = 1, lend-1
!         do i = 1, iend
!            do j = 1, jend
!               ! if( (Pe(i,j).gt.x4(m,l) .and. Pe(i,j).lt.x4(m+1,l)) .and. (ibi(i,j).gt.y4(m,l) .and. ibi(i,j).lt.y4(m,l+1)) )then
!               if( (iei(i,j).gt.x4(m,l) .and. iei(i,j).lt.x4(m+1,l)) .and. (ibi(i,j).gt.y4(m,l) .and. ibi(i,j).lt.y4(m,l+1)) )then
!                  z4(m,l) = z4(m,l) + As_mpm(i,j)
!                  icount(m,l) = icount(m,l) + 1
!                  ! icount2(m,l) = icount2(m,l) + 1
!               endif
!            enddo
!         enddo
!      enddo
!   enddo
!   allocate(z5(icount(mm,ll))) ; p = 1
!   ave_z4(:,:) = 0.0
!   do m = 1, mend-1
!      do l = 1, lend-1
!         do i = 1, iend
!            do j = 1, jend
!               ! if( (Pe(i,j).gt.x4(m,l) .and. Pe(i,j).lt.x4(m+1,l)) .and. (ibi(i,j).gt.y4(m,l) .and. ibi(i,j).lt.y4(m,l+1)) )then
!               if( (iei(i,j).gt.x4(m,l) .and. iei(i,j).lt.x4(m+1,l)) .and. (ibi(i,j).gt.y4(m,l) .and. ibi(i,j).lt.y4(m,l+1)) )then
!                  ave_z4(m,l) = z4(m,l) / icount(m,l)
!                  ! write(*,*)icount(m,l), ave_z4(m,l)
!                  if(m.eq.mm .and. l.eq.ll)then
!                     z5(p) = As_mpm(i,j)
!                     p = p + 1                    
!                  endif
!               endif
!            enddo
!         enddo
!      enddo
!   enddo
!   ! call output_for_3d_figure
  if(n .eq. nend)then
     call plend
  endif
end subroutine drw_cntr
!===================================================================================================================================








! Plot line
!===================================================================================================================================
subroutine plot_line(iend, xp, yp, icol)
!===================================================================================================================================
  use plplot
  implicit none
  integer,intent(in) :: iend, icol
  real*8,intent(in)  :: xp(0:iend-1), yp(0:iend-1)

  call plcol0(icol)
  call plline(xp, yp)
end subroutine plot_line
!===================================================================================================================================


! Plot line reverse
!===================================================================================================================================
subroutine plot_line_rvrs(iend, xp, yp, icol)
!===================================================================================================================================
  use plplot
  implicit none
  integer            :: i
  integer,intent(in) :: iend, icol
  real*8,intent(in)  :: xp(0:iend-1), yp(0:iend-1)
  real*8             :: xpr(0:iend-1), ypr(0:iend-1)  

  do i = 0, iend-1
     ! xpr(i) = xp(iend-i-1)
     ypr(i) = yp(iend-i-1)
  enddo
  
  call plcol0(icol)
  call plline(xp, ypr)
end subroutine plot_line_rvrs
!===================================================================================================================================

! Plot point
!===================================================================================================================================
subroutine plot_point(iend, xp, yp, icol, symbl)
!===================================================================================================================================  
  use plplot
  implicit none
  integer,intent(in) :: iend, icol, symbl
  real*8,intent(in)  :: xp(0:iend-1), yp(0:iend-1)

  call plcol0(icol)
  call plpoin(xp, yp, symbl)
end subroutine plot_point
!===================================================================================================================================

! Plot point rvrs
!===================================================================================================================================
subroutine plot_point_rvrs(iend, xp, yp, icol, symbl)
!===================================================================================================================================  
  use plplot
  implicit none
  integer            :: i
  integer,intent(in) :: iend, icol, symbl
  real*8,intent(in)  :: xp(0:iend-1), yp(0:iend-1)
  real*8             :: xpr(0:iend-1), ypr(0:iend-1)  

  do i = 0, iend-1
     ! xpr(i) = xp(iend-i-1)
     ypr(i) = yp(iend-i-1)
  enddo

  call plcol0(icol)
  call plpoin(xp, ypr, symbl)
end subroutine plot_point_rvrs
!===================================================================================================================================


! Set window
!===================================================================================================================================
subroutine set_window(axmin, axmax, aymin, aymax, axis_b, axis_l, axis_t, pxmin, pxmax, pymin, pymax)
!===================================================================================================================================
  use plplot
  real*8,intent(in)       :: axmin, axmax, aymin, aymax, pxmin, pxmax, pymin, pymax
  character(*),intent(in) :: axis_b, axis_l, axis_t

  call plvpor(pxmin, pxmax, pymin, pymax)
  call plwind(axmin, axmax ,aymin, aymax)
  call plcol0(15) 
  call plbox('bcnst', 0.0d0, 5, 'bcnstv', 0.0d0, 0)
  ! call plbox('bcnstlv', 0.0d0, 5, 'bcnstv', 0.0d0, 0)
  call plmtex('b', 3.0d0, 0.5d0, 0.5d0, trim(axis_b))
  call plmtex('l', 4.5d0, 0.5d0, 0.5d0, trim(axis_l))
  call plmtex('t', 1.d0, 0.5d0, 0.5d0, trim(axis_t))
end subroutine set_window
!===================================================================================================================================

! ! Set window_3d
! !===================================================================================================================================
! subroutine set_window_3d(axmin, axmax, aymin, aymax, axis_b, axis_l, axis_t, pxmin, pxmax, pymin, pymax, pzmin, pzmax)
! !===================================================================================================================================
!   use plplot
!   real*8,intent(in)       :: axmin, axmax, aymin, aymax, pxmin, pxmax, pymin, pymax, pzmin, pzmax
!   character(*),intent(in) :: axis_b, axis_l, axis_t

!   call plcol0(15)
!   call plvpor(pxmin, pxmax, pymin, pymax)
!   call plwind(axmin, axmax ,aymin, aymax)
!   call plw3d(0.5d0, 0.8d0, 1.3d0, -3.0d0, 3.0d0, -3.0d0, 3.0d0, pzmin, pzmax, 40.0d0, 20.0d0)
!   call plbox3('bnstu', axis_b, 0.0d0, 0, 'bnstu', axis_l, 0.0d0, 0,'bcdmnstuv', axis_t, 0.0d0, 0)
!   ! call plmesh(x(:xpts), y(:ypts), z(:xpts,:ypts), 1)
  
! end subroutine set_window_3d
! !===================================================================================================================================

! ! 
! !===================================================================================================================================
! subroutine plot_3d(axmin, axmax, aymin, aymax, axis_b, axis_l, axis_t, pxmin, pxmax, pymin, pymax, pzmin, pzmax)
! !===================================================================================================================================
!   use plplot
!   implicit none
!   integer,intent(in) :: iend, icol, symbl
!   real*8,intent(in)  :: xp(0:iend-1), yp(0:iend-1),z(mend,lend)

!   call plmesh(x(:xpts), y(:ypts), z(:xpts,:ypts), 0)
! end subroutine plot_3d
! !===================================================================================================================================

  
! Set window first_axis
!===================================================================================================================================
subroutine set_window_first_axis(axmin, axmax, aymin, aymax, axis_b, axis_l, axis_t, pxmin, pxmax, pymin, pymax)
!===================================================================================================================================
  use plplot
  real*8,intent(in)       :: axmin, axmax, aymin, aymax, pxmin, pxmax, pymin, pymax
  character(*),intent(in) :: axis_b, axis_l, axis_t

  call plvpor(pxmin, pxmax, pymin, pymax)
  call plwind(axmin, axmax ,aymin, aymax)
  call plcol0(15) 
  call plbox('bcnst', 0.0d0, 5, 'bnstv', 0.0d0, 0)
  call plmtex('b', 3.0d0, 0.5d0, 0.5d0, trim(axis_b))
  call plmtex('l', 5.0d0, 0.5d0, 0.5d0, trim(axis_l))
  call plmtex('t', 1.d0, 0.5d0, 0.5d0, trim(axis_t))
  
end subroutine set_window_first_axis
!===================================================================================================================================

! Set window second axis
!===================================================================================================================================
subroutine set_window_second_axis(axmin, axmax, aymin, aymax, axis_b, axis_r, axis_t, pxmin, pxmax, pymin, pymax)
!===================================================================================================================================
  use plplot
  real*8,intent(in)       :: axmin, axmax, aymin, aymax, pxmin, pxmax, pymin, pymax
  character(*),intent(in) :: axis_b, axis_r, axis_t

  call plvpor(pxmin, pxmax, pymin, pymax)
  call plwind(axmin, axmax ,aymin, aymax)
  call plcol0(15) 
  call plbox('bcnst', 0.0d0, 5, 'cmstv', 0.0d0, 0)
  call plmtex('b', 3.0d0, 0.5d0, 0.5d0, trim(axis_b))
  call plmtex('r', 5.0d0, 0.5d0, 0.5d0, trim(axis_r))
  call plmtex('t', 1.d0, 0.5d0, 0.5d0, trim(axis_t))
  
end subroutine set_window_second_axis
!===================================================================================================================================


! Set window title
!===================================================================================================================================
subroutine set_window_title(axmin, axmax, aymin, aymax, axis_b, axis_l, axis_t, pxmin, pxmax, pymin, pymax, title)
!===================================================================================================================================
  use plplot
  real*8,intent(in)       :: axmin, axmax, aymin, aymax, pxmin, pxmax, pymin, pymax
  character(*),intent(in) :: axis_b, axis_l, axis_t, title

  call plvpor(pxmin, pxmax, pymin, pymax)
  call plwind(axmin, axmax ,aymin, aymax)
  call plcol0(15) 
  call plbox('bcnst', 0.0d0, 5, 'bcnstv', 0.0d0, 0)
  call plmtex('b', 3.0d0, 0.5d0, 0.5d0, trim(axis_b))
  call plmtex('l', 4.5d0, 0.5d0, 0.5d0, trim(axis_l))
  call plmtex('t', 1.d0, 0.5d0, 0.5d0, trim(axis_t))
  call plmtex('t', 1.d0, 0.4d0, 0.1d0, trim(title))
end subroutine set_window_title
!===================================================================================================================================

!##### Drawing a histgram ###########################
!===================================================================================================================================
subroutine plot_hstgrm(iend, x, y, wdth, icol)
!===================================================================================================================================
  use plplot
  implicit none
  integer            :: i
  integer,intent(in) :: iend, icol
  real*8,intent(in)  :: x(iend), y(iend), wdth
  real*8             :: xp(0:4), yp(0:4)

  do i = 1, iend
     xp(0) = x(i) - (wdth/2.d0)
     xp(1) = x(i) + (wdth/2.d0)
     xp(2) = xp(1)
     xp(3) = xp(0)
     xp(4) = xp(0)
     
     yp(0) = 0.d0
     yp(1) = yp(0)
     yp(2) = y(i)
     yp(3) = yp(2)
     yp(4) = yp(0)

     call plcol0(icol)
     call plfill(xp, yp)
     call plcol0(15)     
     call plline(xp, yp)     
  enddo
  
end subroutine plot_hstgrm
!===================================================================================================================================

! Set color of contour
!===================================================================================================================================
subroutine set_clr_cntr(Clr_nm, zmin, zmax)
!===================================================================================================================================  
  use Val_Cmap
  integer,intent(in) :: Clr_nm
  real*8,intent(in)  :: zmin, zmax

  iclrDst = 16
  cqmin = zmin
  cqmax = zmax
  
  if(Clr_nm.eq.3)then
     iru(1) = 0
     igu(1) = 0
     ibu(1) = 255
     iru(2) = 255
     igu(2) = 255
     ibu(2) = 255
     iru(3) = 255
     igu(3) = 0
     ibu(3) = 0
  elseif(Clr_nm.eq.2)then
     iru(1) = 255
     igu(1) = 255
     ibu(1) = 255
     iru(2) = 0
     igu(2) = 0
     ibu(2) = 0
  elseif(Clr_nm.eq.4)then
     iru(1) = 0
     igu(1) = 0
     ibu(1) = 255
     iru(2) = 0
     igu(2) = 255
     ibu(2) = 255
     iru(3) = 255
     igu(3) = 255
     ibu(3) = 0
     iru(4) = 255
     igu(4) = 0
     ibu(4) = 0
  elseif(Clr_nm.eq.5)then
     iru(1) = 0
     igu(1) = 0
     ibu(1) = 255
     ! iru(1) = 255
     ! igu(1) = 255
     ! ibu(1) = 255
     iru(2) = 0
     igu(2) = 255
     ibu(2) = 255
     iru(3) = 0
     igu(3) = 255
     ibu(3) = 20
     iru(4) = 255
     igu(4) = 255
     ibu(4) = 0
     iru(5) = 255
     igu(5) = 0
     ibu(5) = 0
     ! iru(5) = 255
     ! igu(5) = 255
     ! ibu(5) = 255
  endif

  if(Clr_nm .eq. 2)call remaprgb_TwoClr
  if(Clr_nm .eq. 3)call remaprgb_ThreClr
  if(Clr_nm .eq. 4)call remaprgb_FourClr
  if(Clr_nm .eq. 5)call remaprgb_FiveClr

end subroutine set_clr_cntr
!===================================================================================================================================


! Draw contour
!===================================================================================================================================
subroutine plot_cntr(imin, imax, jmin, jmax, x, y, z)
!===================================================================================================================================  
  use Val_Cmap
  integer               :: i, j
  integer,intent(in)    :: imin, imax, jmin, jmax
  real*8                :: zave
  real*8,dimension(0:4) :: xp, yp
  real*8,dimension(imin:imax,jmin:jmax),intent(in) :: x, y, z   

  do j = jmin, jmax-1
     do i = imin, imax-1
        xp(0) = x(i,j)
        yp(0) = y(i,j)
        xp(1) = x(i+1,j)
        yp(1) = y(i+1,j)
        xp(2) = x(i+1,j+1)
        yp(2) = y(i+1,j+1)
        xp(3) = x(i,j+1)
        yp(3) = y(i,j+1)
        xp(4) = x(i,j)
        yp(4) = y(i,j)
        zave = sum(z(i:i+1,j:j+1)) / 4.d0
        call clrdqnty(zave)        
        call plfill(xp, yp)
     enddo
  enddo

end subroutine plot_cntr
!===================================================================================================================================


! Draw color bars
!===================================================================================================================================
subroutine colbar(axmin, axmax, aymin, aymax, pxmin, pxmax, pymin, pymax, barname)
!===================================================================================================================================  
  use Val_Cmap
  integer                 :: i, j 
  real*8                  :: cqt, xp(0:4), yp(0:4)
  real*8,intent(in)       :: axmin, axmax, aymin, aymax, pxmin, pxmax, pymin, pymax
  character(*),intent(in) :: barname

  call plvpor(pxmin, pxmax, pymin, pymax)
  call plwind(axmin, axmax, aymin, aymax)
  call plcol0(15)
  call plbox('bc', 0.0d0, 0, 'bcmst', 0.0d0, 0)      
  call plmtex('l', 1.d0, 0.5d0, 0.5d0, barname)   

  dpRng   = aymax - aymin
  RngDlt  = dpRng / iclrDst
  
  do i = 0, iclrDst-1
     xp(0) = axmin
     yp(0) = aymin + RngDlt * (i)
     xp(1) = axmax
     yp(1) = yp(0)
     xp(2) = xp(1)
     yp(2) = aymax + RngDlt * (i)
     xp(3) = xp(0)
     yp(3) = yp(2)
     xp(4) = xp(0)
     yp(4) = yp(0)
     cqt   = rng(i)
     call clrdqnty(cqt)
     call plfill(xp, yp)
  enddo

end subroutine colbar
!===================================================================================================================================


! ---- Output for_3d_figure ----
!===================================================================================================================================
subroutine output_for_3d_figure
  use mndr
  call system('mkdir ../07_3D_index/01_input_number')
  call system('mkdir ../07_3D_index/01_input_number/'//trim(case))    

  open(100,file='../07_3D_index/01_input_number/'//trim(case)//'/'//trim(fln(n)(1:8)//'.txt'))  
  do l = 1, lend
     if(l.eq.1)then
        write(100,*)' ',x4(1:mend,1)
     else
        write(100, *)y4(1,l), ave_z4(1:mend,l)
     endif
  enddo
  close(100)
end subroutine output_for_3d_figure
!===================================================================================================================================
