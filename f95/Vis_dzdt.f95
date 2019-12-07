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

  pdx   = 0.32d0 ; pdy   = 0.15d0

  
  ! ---- cntr of D.B.L ---------------
  barsname = 'D.B.L [m]'
  axis_l   = 'Width(m)'
  axis_t   = trim(fln2(n))
  axis_b   = ''
  clrnm    = 3

  xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
  xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
  ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
  ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!!$  zmin = minval(dev_bl) - (maxval(dev_bl)-minval(dev_bl))*0.1d0
!!$  zmax = maxval(dev_bl) + (maxval(dev_bl)-minval(dev_bl))*0.1d0
  zmin = 0.04d0 ; zmax = 0.06d0
!!$  zmin = 0.03d0 ; zmax = 0.05d0

  pxmin = 0.1d0  ; pxmax = 0.90d0
  pymin = 0.85d0 ; pymax = 0.95d0  
    
  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
  call set_clr_cntr(clrnm, zmin, zmax)
  call plot_cntr(1, iend, 1, jend, x, y, dev_bl)
  do i = 1, iend
     do j = 1, jend
        if(tautau(i,j).lt.1)then
           call plot_point(1,x(i,j), y(i,j), 7, 16)
        endif
     enddo
  enddo
  call plcol0(9)
  call pljoin(x(1,jend-2),y(1,jend-2),x(iend,jend-2),y(iend,jend-2))
  call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
  !------------------------------------------

  !---- cntr of dz_obs ---------------
  barsname = 'dz_obs[mm]'
  axis_l   = 'Width(m)'
  axis_t   = ''
  axis_b   = ''
  clrnm    = 3

  xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
  xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
  ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
  ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
  zmin =  -10.0d0 ; zmax =   10.0d0
  
  pymin = pymin - pdy ; pymax = pymax - pdy
    
  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
  call set_clr_cntr(clrnm, zmin, zmax)
  call plot_cntr(1, iend, 1, jend, x, y, dz_obs)
  do i = 1, iend
     do j = 1, jend
        if(tautau(i,j).lt.1)then
           call plot_point(1,x(i,j), y(i,j), 7, 16)
        endif
     enddo
  enddo
  call plcol0(9)
  call pljoin(x(1,jend-2),y(1,jend-2),x(iend,jend-2),y(iend,jend-2))
  call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
  !------------------------------------------
  
!!$  !---- cntr of dz_Mu ---------------
!!$  barsname = 'dz_mpm[mm]'
!!$  axis_l   = 'Width(m)'
!!$  axis_t   = ''
!!$  axis_t   = trim(fln2(n))
!!$  axis_b   = ''
!!$  axis_b   = ''
!!$  clrnm    = 3
!!$
!!$  xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!!$  xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!!$  ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!!$  ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!!$  zmin =  -10.0d0 ; zmax =   10.0d0
!!$
!!$  pymin = pymin - pdy ; pymax = pymax - pdy
!!$
!!$  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!!$  call set_clr_cntr(clrnm, zmin, zmax)
!!$  call plot_cntr(1, iend, 1, jend, x, y, dz_mpm)
!!$  call plcol0(9)
!!$  call pljoin(x(1,jend-2),y(1,jend-2),x(iend,jend-2),y(iend,jend-2))
!!$  call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
!!$  !------------------------------------------

!---- cntr of dz_Mu ---------------
  barsname = 'dz_Mu[mm]'
  axis_l   = 'Width(m)'
  axis_t   = ''
  axis_t   = trim(fln2(n))
  axis_b   = ''
  axis_b   = ''
  clrnm    = 3

  xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
  xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
  ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
  ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
  zmin =  -10.0d0 ; zmax =   10.0d0

  pymin = pymin - pdy ; pymax = pymax - pdy
    
  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
  call set_clr_cntr(clrnm, zmin, zmax)
  call plot_cntr(1, iend, 1, jend, x, y, dz_Mu)
  do i = 1, iend
     do j = 1, jend
        if(tautau(i,j).lt.1)then
           call plot_point(1,x(i,j), y(i,j), 7, 16)
        endif
     enddo
  enddo
  call plcol0(9)
  call pljoin(x(1,jend-2),y(1,jend-2),x(iend,jend-2),y(iend,jend-2))
  call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
  !------------------------------------------

!!$  !---- cntr of dz_Mu ---------------
!!$  barsname = 'dz_Md[mm]'
!!$  axis_l   = 'Width(m)'
!!$  axis_t   = ''
!!$  axis_t   = trim(fln2(n))
!!$  axis_b   = ''
!!$  axis_b   = ''
!!$  clrnm    = 3
!!$
!!$  xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!!$  xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!!$  ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!!$  ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!!$  zmin =  -10.0d0 ; zmax =   10.0d0
!!$
!!$  pymin = pymin - pdy ; pymax = pymax - pdy
!!$    
!!$  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!!$  call set_clr_cntr(clrnm, zmin, zmax)
!!$  call plot_cntr(1, iend, 1, jend, x, y, dz_Md)
!!$  call plcol0(9)
!!$  call pljoin(x(1,jend-2),y(1,jend-2),x(iend,jend-2),y(iend,jend-2))
!!$  call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
!!$!------------------------------------------ 

!!$  !---- cntr of D.B.L ---------------
!!$  barsname = 'dif_dzmpm [%]'
!!$  axis_l   = 'Width(m)'
!!$  axis_t   = ''
!!$  axis_t   = ''
!!$  axis_b   = ''
!!$  clrnm    = 2
!!$
!!$  xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!!$  xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!!$  ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!!$  ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!!$  zmin = 0.0d0
!!$  zmax = 500.0d0
!!$
!!$  pymin = pymin - pdy ; pymax = pymax - pdy
!!$    
!!$  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!!$  call set_clr_cntr(clrnm, zmin, zmax)
!!$  call plot_cntr(1, iend, 1, jend, x, y, dif_dzmpm)
!!$  call plcol0(9)
!!$  call pljoin(x(1,jend-2),y(1,jend-2),x(iend,jend-2),y(iend,jend-2))
!!$  call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
!!$!------------------------------------------

!---- cntr of D.B.L ---------------
  barsname   = 'dif_dzu [%]'
  axis_l   = 'Width(m)'
  axis_t = ''
  axis_t   = trim(fln2(n))
  axis_b   = ''
  axis_b   = 'Distance from upstream [m]'
  clrnm    = 2

  xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
  xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
  ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
  ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
  zmin = 0.0d0
  zmax = 2000.0d0

  pymin = pymin - pdy ; pymax = pymax - pdy
    
  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
  call set_clr_cntr(clrnm, zmin, zmax)
  call plot_cntr(1, iend, 1, jend, x, y, dif_dzu)
  do i = 1, iend
     do j = 1, jend
        if(tautau(i,j).lt.1)then
           call plot_point(1,x(i,j), y(i,j), 7, 16)
        endif
     enddo
  enddo
  call plcol0(9)
  call pljoin(x(1,jend-2),y(1,jend-2),x(iend,jend-2),y(iend,jend-2))
  call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
!------------------------------------------


!!$!---- cntr of D.B.L ---------------
!!$  barsname = 'dif_dzd [%]'
!!$  axis_l   = 'Width(m)'
!!$  axis_t   = ''
!!$  axis_t   = trim(fln2(n))
!!$  axis_b   = ''
!!$  axis_b   = 'Distance from upstream [m]'
!!$  clrnm    = 2
!!$
!!$  xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!!$  xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!!$  ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!!$  ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!!$  zmin = 0.0d0
!!$  zmax = 500.0d0
!!$
!!$  pymin = pymin - pdy ; pymax = pymax - pdy
!!$    
!!$  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!!$  call set_clr_cntr(clrnm, zmin, zmax)
!!$  call plot_cntr(1, iend, 1, jend, x, y, dif_dzd)
!!$  call plcol0(9)
!!$  call pljoin(x(1,jend-2),y(1,jend-2),x(iend,jend-2),y(iend,jend-2))
!!$  call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
!!$!------------------------------------------

!!$  !---- cntr of D.B.L ---------------
!!$  barsname = 'dif_dzu [%]'
!!$  axis_l   = 'Width(m)'
!!$  axis_t   = ''
!!$  axis_t   = ''
!!$  axis_b   = ''
!!$  axis_b   = ''
!!$  clrnm    = 2
!!$
!!$  xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!!$  xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!!$  ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
!!$  ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
!!$  zmin = 0.0d0
!!$  zmax = 500.0d0
!!$
!!$  pymin = pymin - pdy ; pymax = pymax - pdy
!!$    
!!$  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!!$  call set_clr_cntr(clrnm, zmin, zmax)
!!$  call plot_cntr(1, iend, 1, jend, x, y, dif_dzu)
!!$  call plcol0(9)
!!$  call pljoin(x(1,jend-2),y(1,jend-2),x(iend,jend-2),y(iend,jend-2))
!!$  call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
!!$  !------------------------------------------



!!$  !---- lng of D.B.L ---------------
!!$  axis_l   = 'D.B.L [m]'
!!$  axis_t   = ''
!!$  axis_t   = trim(fln2(n))
!!$  axis_b   = ''
!!$  axis_b   = ''
!!$
!!$  xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!!$  xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!!$  ymin = 0.04d0
!!$  ymax = 0.06d0
!!$
!!$  pxmin = 0.55d0  ; pxmax = 0.9d0
!!$  pymin = 0.85d0 ; pymax = 0.95d0  
!!$    
!!$  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!!$  call plot_line(iend, x(:,jend-2), dev_bl(:,jend-2), 1)
!!$  !------------------------------------------
!!$
!!$  !---- lng of dif_dz ---------------
!!$  axis_l   = 'dif_dz [%] '
!!$  axis_t   = ''
!!$  axis_t   = ''
!!$  axis_b   = ''
!!$  axis_b   = 'Distance from upstream [m]'
!!$
!!$  xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
!!$  xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
!!$  ymin = 0.0d0
!!$  ymax = 800.0d0
!!$
!!$  pymin = pymin - pdy ; pymax = pymax - pdy
!!$    
!!$  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!!$  call plot_line(iend, x(:,jend-2), dif_dzu(:,jend-2), 9)
!!$  call plot_line(iend, x(:,jend-2), dif_dzmpm(:,jend-2), 1)
!!$  !------------------------------------------

!!$  !---- time change ---------------
!!$  axis_l   = 'time [min]'
!!$  axis_t   = 'z [m]'
!!$  axis_b   = ''
!!$  
!!$  pxmin = 0.75d0 ; pxmax = 0.9d0
!!$  pymin = 0.1d0 ; pymax = 0.95d0
!!$
!!$  xmin = 0.048d0 ; xmax = 0.050d0
!!$  ymax = minval(time2(1:10)) - (maxval(time2(1:10))-minval(time2(1:10)))*0.02d0
!!$  ymin = maxval(time2(1:10)) + (maxval(time2(1:10))-minval(time2(1:10)))*0.02d0
!!$  
!!$  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
!!$  call plot_point(n, bl2(iend/2,jend/2,1:n), time2(1:n), 9, 21)
!!$  !------------------------------------------  

!!$! ---- Histgram of dzdt(obs) - dzdt(mpm) ------------
!!$  axis_l   = 'number of points'
!!$  axis_t   = '' 
!!$  axis_b   = 'Pe'
!!$  pdy      = 0.30d0
!!$  
!!$  pxmin = 0.65d0 ; pxmax = 0.85d0
!!$  pymin = pymin - pdy ; pymax = pymax - pdy
!!$
!!$  xmin = 0.0d0 ; xmax = 0.02d0
!!$  ymin = 0.0d0 ; ymax = 50.0d0
!!$
!!$  wdth = (xmax - xmin) / 40
!!$  pend = (xmax - xmin)/wdth + 2
!!$  allocate(hst_x(pend),hst_y(pend))
!!$
!!$  do p = 1, pend
!!$     con = 0
!!$     do j = 1, jend
!!$        do i = 1, iend
!!$           if(tautau(i,j) .ge. 1.0 )then
!!$              if(((p - 1)*wdth + xmin < Pe_d(i,j)) .and. (Pe_d(i,j) <= p*wdth + xmin))then
!!$                 con = con + 1
!!$              endif
!!$           endif
!!$        enddo
!!$     enddo
!!$     hst_x(p) = (p - 1)*wdth + xmin
!!$     hst_y(p) = con / (iend*jend) * 100
!!$  enddo
!!$  
!!$  call set_window(xmin - wdth/2., xmax+wdth/2, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax) 
!!$  call plot_hstgrm(pend, hst_x, hst_y, wdth, 3)
!!$  
!!$  deallocate(hst_x, hst_y)
!!$! ------------------------------------------

!!$! ---- Histgram of dzdt(obs) - dzdt(mpm) ------------
!!$  axis_l   = 'number of points'
!!$  axis_t   = '' 
!!$  axis_b   = 'Pe'
!!$  pdy      = 0.30d0
!!$  
!!$  pxmin = 0.65d0 ; pxmax = 0.85d0
!!$  pymin = pymin - pdy ; pymax = pymax - pdy
!!$
!!$  xmin = 0.0d0 ; xmax = 0.02d0
!!$  ymin = 0.0d0 ; ymax = 50.0d0
!!$
!!$  wdth = (xmax - xmin) / 40
!!$  pend = (xmax - xmin)/wdth + 2
!!$  allocate(hst_x(pend),hst_y(pend))
!!$
!!$  do p = 1, pend
!!$     con = 0
!!$     do j = 1, jend
!!$        do i = 1, iend
!!$           if(((p - 1)*wdth + xmin < Pe_d(i,j)) .and. (Pe_d(i,j) <= p*wdth + xmin))then
!!$              con = con + 1
!!$           endif
!!$        enddo
!!$     enddo
!!$     hst_x(p) = (p - 1)*wdth + xmin
!!$     hst_y(p) = con / (iend*jend) * 100
!!$  enddo
!!$  
!!$  call set_window(xmin - wdth/2., xmax+wdth/2, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax) 
!!$  call plot_hstgrm(pend, hst_x, hst_y, wdth, 3)
!!$  
!!$  deallocate(hst_x, hst_y)
!!$! ------------------------------------------

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
  call plmtex('b', 3.0d0, 0.5d0, 0.5d0, trim(axis_b))
  call plmtex('l', 4.5d0, 0.5d0, 0.5d0, trim(axis_l))
  call plmtex('t', 1.d0, 0.5d0, 0.5d0, trim(axis_t))
  
end subroutine set_window
!===================================================================================================================================

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
  call plmtex('l', 4.5d0, 0.5d0, 0.5d0, trim(axis_l))
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
  call plmtex('r', 2.0d0, 0.5d0, 0.5d0, trim(axis_r))
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
  elseif(Clr_nm.eq.5)then
     iru(1) = 0
     igu(1) = 0
     ibu(1) = 255
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
