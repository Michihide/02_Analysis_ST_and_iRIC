Module mndr
  implicit none
  
  integer                             :: i,    j,    n,    m,    p, count
  integer                             :: iend, jend, nend, mend, pend
  integer                             :: clrnm, intbl
  real*8                              :: pi, q, b, g, dx, tauc, s, dm

  ! ---- set_pss ----
  character*100 :: pss_st, pss_iric

    
  ! ---- read_roughness ----
  real*8,dimension(:),allocatable     :: nc

  
  ! ---- read vtk for dz/dt & du/dt (obs) ----
  real*8,dimension(:,:,:),allocatable :: bl2, u_cal2, mut
  
  
  ! ---- read_schalar_vtk & read_vector_vtk ----
  real*8,dimension(:,:),allocatable   :: x, y
  real*8,dimension(:,:),allocatable   :: dep_obs, wl_obs, dev_bl, bl, dwl_obs
  real*8,dimension(:,:),allocatable   :: dep_cal, wl_cal, u_cal, v_cal

  
  ! ---- cal_slope & curvature & velocity ----
  real*8,dimension(:,:),allocatable   :: iwj_obs, iwi_obs, iw_obs
  real*8,dimension(:,:),allocatable   :: iwj_cal, iwi_cal, iw_cal
  real*8,dimension(:,:),allocatable   :: ibi    , ibj    , ib  
  real*8,dimension(:,:),allocatable   :: cbj    , cbi    , cb
  real*8,dimension(:,:),allocatable   :: iej    , iei    , ie
  real*8,dimension(:,:),allocatable   :: Fr     , vel_cal


  ! ---- cal shields_number &  M_equation ----
  real*8,dimension(:,:),allocatable   :: taue_i, taue_j, taue , tautau
  real*8,dimension(:,:),allocatable   :: qb_mpm, dz_mpm
  real*8,dimension(:,:),allocatable   :: Pe_s  , Ms    , dz_Ms
  real*8,dimension(:,:),allocatable   :: Pe_u  , Mu    , dz_Mu
  real*8,dimension(:,:),allocatable   :: Pe_d  , Md    , Df   , dz_Md, dz_amf
  real*8,dimension(:,:),allocatable   :: lcl   , adv   , dff  , frc
  real*8,dimension(:,:),allocatable   :: Mx    , My    , Mxy
  real*8,dimension(:),allocatable     :: ave_Mu
  
  ! ---- cal and cmp dzdt & dhdx ----
  real*8,dimension(:,:),allocatable   :: dz_obs
  real*8,dimension(:,:),allocatable   :: dif_dep  , dif_dzs , dif_dzd, dif_dzu, dif_dzmpm
  real*8,dimension(:,:),allocatable   :: dhdx_obs , dhdy_obs
  real*8,dimension(:,:),allocatable   :: dhdx_Mu  , dhdy_Mu
  real*8,dimension(:,:),allocatable   :: dhdx_Ms  , dhdx_Md
  real*8,dimension(:,:),allocatable   :: dif_dhdxs, dif_dhdys
  real*8,dimension(:,:),allocatable   :: dif_dhdxu, dif_dhdyu

  ! ---- cal average crs ans lng ----
  real*8,dimension(:),allocatable     :: ave_crs_dbl, ave_crs_dwl
  
  ! ---- cmp_each_item ---- 
  real*8,dimension(:,:),allocatable   :: dudt, dudx, asu

  
  ! ----  cal_nondim_celerity ----
  real*8,dimension(:,:),allocatable   :: tau_nd, fr_nd, Mu_nd
  real*8                              :: Mu_nd0, tau_nd0, fr_nd0, iei0
  
  ! ---- output vtk ----  
  integer                             :: ic, jc, ijc, intrvl_i, intrvl_j

  
  ! ---- drw_cntr ----
  integer :: plparseopts_rc
  character*100                          :: chr, case, stage, nchr
  character*4,dimension(:),allocatable   :: fln2
  character*100                          :: barsname, axis_b, axis_l, axis_t, axis_r
  real*8                                 :: pxmin, pxmax, pymin, pymax, pdy, pzmin, pzmax
  real*8                                 :: xmin, xmax, ymin, ymax, zmin, zmax
  real*8                                 :: lim
  real*8                                 :: con, con1, wdth
  real*8,dimension(:),allocatable        :: hst_x, hst_y
  character*100,dimension(:),allocatable :: fln
  character*4,dimension(:),allocatable   :: time1
  real*8,dimension(:),allocatable        :: time2

End Module mndr
