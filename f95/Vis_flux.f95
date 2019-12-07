! **** Draw 2d contour ****
subroutine drw_cntr
  use mndr
  use plplot

  if(n.eq.1 .and. i.eq.1)then
     plparseopts_rc = plparseopts(PL_PARSE_FULL)
     if(plparseopts_rc .ne. 0) stop "plparseopt error"
     call plscol0( 0,255,255,255)
     call plscol0(15,  0,  0,  0)
     call plscmap0n(0)
     call plinit
  endif

  call pladv(0)
  call plschr(1.5d0, 2.0d0)  

  pdx   = 0.32d0 ; pdy   = 0.18d0

  
! ---- cntr of D.B.L ---------------
  barsname = 'D.B.L [m]'
  axis_l   = 'Width(m)'
  ! axis_t = ''
  axis_t   = trim(fln2(n))
  axis_b   = ''
  ! axis_b   = 'Distance from upstream [m]'
  clrnm    = 3

  xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
  xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
  ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
  ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
  zmin = 0.04d0 ; zmax = 0.06d0

  pxmin = 0.1d0  ; pxmax = 0.90d0
  pymin = 0.82d0 ; pymax = 0.95d0  
    
  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
  call set_clr_cntr(clrnm, zmin, zmax)
  call plot_cntr(1, iend, 1, jend, x, y, dev_bl)
  call plcol0(9)
  call pljoin(x(i,1),y(i,1),x(i,jend),y(i,jend))
  call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
!------------------------------------------


  
! ---- cntr of tau ---------------
  barsname = 'tau'
  axis_l   = 'Width(m)'
  axis_t   = trim(fln2(n))
  axis_b   = ''
  clrnm    = 5

  xmin = minval(x) - (maxval(x)-minval(x))*0.02d0
  xmax = maxval(x) + (maxval(x)-minval(x))*0.02d0
  ymin = minval(y) - (maxval(y)-minval(y))*0.1d0
  ymax = maxval(y) + (maxval(y)-minval(y))*0.1d0
  zmin = minval(taue_i) - (maxval(taue_i)-minval(taue_i))*0.1d0
  zmax = maxval(taue_i) + (maxval(taue_i)-minval(taue_i))*0.1d0
  ! zmin = 0.04d0 ; zmax = 0.06d0

  pymin = pymin - pdy ; pymax = pymax - pdy
  
  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)
  call set_clr_cntr(clrnm, zmin, zmax)
  call plot_cntr(1, iend, 1, jend, x, y, taue_i)
  call plcol0(9)
  call pljoin(x(i,1),y(i,1),x(i,jend),y(i,jend))
  call colbar(0.d0, 1.d0, zmin, zmax, pxmax+0.02d0, pxmax+0.03d0, pymin, pymax, trim(barsname))
!------------------------------------------


  
!---- cross section of tau & flux ---------------
  axis_l   = 'tau'
  axis_t   = trim(fln2(n))
  axis_b   = ''

  pxmin = 0.1d0  ; pxmax = 0.40d0
  pymin = pymin - pdy ; pymax = pymax - pdy

  xmin = minval(y) - (maxval(y)-minval(y))*0.02d0
  xmax = maxval(y) + (maxval(y)-minval(y))*0.02d0
  ymin = minval(taue_i(i,:)) - (maxval(taue_i(i,:))-minval(taue_i(i,:)))*0.1d0
  ymax = maxval(taue_i(i,:)) + (maxval(taue_i(i,:))-minval(taue_i(i,:)))*0.1d0
  
  call set_window_first_axis(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)  
  call plot_line(jend, y(i,:), taue_i(i,:), 3)

  ! axis_r   = 'flux'  
  ! ymin = minval(taue_i) - (maxval(taue_i)-minval(taue_i))*0.1d0
  ! ymax = maxval(taue_i) + (maxval(taue_i)-minval(taue_i))*0.1d0
  ! call set_window_second_axis(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_r), trim(axis_t), pxmin, pxmax, pymin, pymax)  
  ! call plot_line(jend, x(i,:), dz_obs(i,:), 1)
!------------------------------------------  



!---- cross section of bed ---------------
  axis_l   = 'D.B.L [m]'
  axis_t   = ''
  axis_b   = 'Distance from left bank [m]'

  pymin = pymin - pdy ; pymax = pymax - pdy

  xmin = minval(y) - (maxval(y)-minval(y))*0.02d0
  xmax = maxval(y) + (maxval(y)-minval(y))*0.02d0
  ymin = minval(dev_bl(i,:)) - (maxval(dev_bl(i,:))-minval(dev_bl(i,:)))*0.1d0
  ymax = maxval(dev_bl(i,:)) + (maxval(dev_bl(i,:))-minval(dev_bl(i,:)))*0.1d0
  
  call set_window(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_l), trim(axis_t), pxmin, pxmax, pymin, pymax)  
  call plot_line(jend, y(i,:), dev_bl(i,:), 1)

  ! axis_r   = 'flux'  
  ! ymin = minval(taue_i) - (maxval(taue_i)-minval(taue_i))*0.1d0
  ! ymax = maxval(taue_i) + (maxval(taue_i)-minval(taue_i))*0.1d0
  ! call set_window_second_axis(xmin, xmax, ymin, ymax, trim(axis_b), trim(axis_r), trim(axis_t), pxmin, pxmax, pymin, pymax)  
  ! call plot_line(jend, x(i,:), dz_obs(i,:), 1)
!------------------------------------------  

  
  if(n .eq. nend .and. i .eq. iend)then
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
  ! elseif(Clr_nm.eq.4)then
  !    iru(1) = 0
  !    igu(1) = 255
  !    ibu(1) = 255
  !    iru(2) = 255
  !    igu(2) = 255
  !    ibu(2) = 255
  !    iru(3) = 255
  !    igu(3) = 153
  !    ibu(3) = 0     
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
