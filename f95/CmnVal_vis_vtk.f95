Module mndr
  implicit none
  
  integer                                :: m,    p, count
  integer                                :: mend, pend
  integer                                :: clrnm, intbl
 

  ! ---- set_parameter ----
  real*8                                 :: pi, q, b, g, s, dm

  ! ---- read_fln ----
  integer                                :: n, nend
  character*100,dimension(:),allocatable :: fln
  character*4,dimension(:),allocatable   :: fln2
  character*4,dimension(:),allocatable   :: time1
  real*8,dimension(:),allocatable        :: time2

  ! ---- read_array ----
  integer                                :: i, iend, j, jend
 
  ! ---- set_pss ----
  character*100                          :: pss_st, pss_iric

  ! ---- read_roughness ----
  real*8,dimension(:),allocatable        :: nc

  ! ---- cal_critical_shear_stress ----
  real*8                                 :: tauc

  ! ==== set_array at bellow ========  ! You could have to refer to delete if you don't wanna use them
  ! ---- read vtk for dz/dt & du/dt (obs) ----
  real*8,dimension(:,:,:),allocatable    :: bl2, u_cal2, mut  
  
  ! ---- read vtk of ST (updated 2019/12/05 for data of 2019) ----
  real*8,dimension(:,:),allocatable      :: x, y
  real*8,dimension(:,:),allocatable      :: dep_obs, wl_obs, dev_bl, bl, dwl_obs

  ! ---- read vtk of iRIC ----
  real*8,dimension(:,:),allocatable      :: dep_cal, wl_cal, u_cal, v_cal
  
  ! ---- cal_slope & curvature & velocity ----
  real*8,dimension(:,:),allocatable      :: iwj_obs, iwi_obs, iw_obs
  real*8,dimension(:,:),allocatable      :: iwj_cal, iwi_cal, iw_cal
  real*8,dimension(:,:),allocatable      :: ibi    , ibj    , ib  
  real*8,dimension(:,:),allocatable      :: cbj    , cbi    , cb
  real*8,dimension(:,:),allocatable      :: iej    , iei    , ie
  real*8,dimension(:,:),allocatable      :: Fr     , vel_cal

  ! ---- calculate shields_number  &  M equation  &  Exner equation ----
  real*8,dimension(:,:),allocatable      :: taue_i, taue_j, taue , tautau
  real*8,dimension(:,:),allocatable      :: qb_mpm, dz_mpm
  real*8,dimension(:,:),allocatable      :: Pe_s  , Ms    , dz_Ms
  real*8,dimension(:,:),allocatable      :: Pe_u  , Mu    , dz_Mu
  real*8,dimension(:,:),allocatable      :: Pe_d  , Md    , Df   , dz_Md, dz_amf
  real*8,dimension(:,:),allocatable      :: lcl   , adv   , dff  , frc
  real*8,dimension(:,:),allocatable      :: Mx    , My    , Mxy
  real*8,dimension(:),allocatable        :: ave_Mu

  ! ---- cal perstent difference between 2 physical quantity ----
  character*100                          :: diano

  ! ---- calculate temporal slope ----
  real*8,dimension(:,:),allocatable      :: dz_obs  

  ! ---- calculate perstent difference between 2 physical quantity ----
  real*8,dimension(:,:),allocatable      :: dif_dep  , dif_dzs , dif_dzd, dif_dzu, dif_dzmpm
  
  ! ---- I should make subroutine to calculate below parameter ----  
  real*8,dimension(:,:),allocatable      :: dhdx_obs , dhdy_obs
  real*8,dimension(:,:),allocatable      :: dhdx_Mu  , dhdy_Mu
  real*8,dimension(:,:),allocatable      :: dhdx_Ms  , dhdx_Md
  real*8,dimension(:,:),allocatable      :: dif_dhdxs, dif_dhdys
  real*8,dimension(:,:),allocatable      :: dif_dhdxu, dif_dhdyu

  ! ---- cal average crs ans lng (unuse now(2019/12/05), It's just stay in order to use if I wanna use) ----
  real*8,dimension(:),allocatable        :: ave_crs_dbl, ave_crs_dwl
  
  ! ---- cmp_each_item  (unuse now(2019/12/05), It's just stay in order to use if I wanna use) ---- 
  real*8,dimension(:,:),allocatable      :: dudt, dudx, asu
  
  ! ----  cal_nondim_celerity  (unuse now(2019/12/05), It's just stay in order to use if I wanna use) ----
  real*8,dimension(:,:),allocatable      :: tau_nd, fr_nd, Mu_nd
  real*8                                 :: Mu_nd0, tau_nd0, fr_nd0, iei0
  ! ===============================================================================
  
  ! ---- output vtk ----  
  integer                                :: ic, jc, ijc, intrvl_i, intrvl_j
  
  ! ---- drw_cntr ----
  integer :: plparseopts_rc
  character*100                          :: chr     , case  , stage , nchr
  character*100                          :: barsname, axis_b, axis_l, axis_t, axis_r
  real*8                                 :: pxmin   , pxmax , pymin , pymax , pdy   , pzmin, pzmax
  real*8                                 :: xmin    , xmax  , ymin  , ymax  , zmin  , zmax
  real*8                                 :: con     , con1  , wdth
  real*8,dimension(:),allocatable        :: hst_x   , hst_y
  real*8                                 :: blim

End Module mndr
