module para_fige
  implicit none
  integer          ::   ista,    lg,lgcmdx,  lsta,    lt
  integer          ::     lz,   mtb,   mtt,  ndir,  neqt
  integer          ::   nind,   nmx,  nobj,  nsta
  parameter(  ndir=3     )
  parameter(  nind=3     )
  parameter(  neqt=7     )
  parameter(  nsta=50    )
  parameter(  lsta=7     )
  parameter(  ista=2     )
  parameter(  nobj=600   )
  parameter(   nmx=500   )
  parameter(lgcmdx=1316  )
end module para_fige
!
module boundary
  use para_fige
implicit none
  integer          ::  kexl,  nbd,nbdko
  integer         ,allocatable :: crotr(:),imaxb(:),iminb(:),jmaxb(:),jminb(:)
  integer         ,allocatable :: kmaxb(:),kminb(:),  lbd(:),lbdko(:), mdnc(:)
  integer         ,allocatable ::   mmb(:),  mpb(:),  mpc(:), mper(:),  mpn(:)
  integer         ,allocatable ::   mpr(:),  nba(:), nbdc(:), ndcc(:), ndlb(:)
  integer         ,allocatable ::  ndrr(:), nfba(:), nfbc(:), nfbn(:), nfbr(:)
  integer         ,allocatable ::  nfei(:),srotr(:),new2old_f(:),tab_raccord(:)
  double precision,allocatable :: bc(:,:)
  character(len=4),allocatable :: cl(:)
  character(len=2),allocatable :: indfl(:)
end module boundary
!
module maillage
  use para_fige
implicit none
  integer          ::    kcaldis,   kecrdis,   klecdis,       kvn,       lgx
  integer          ::        lzx,   mdimtbx,   mdimtcx,   mdimtnx,   mdimtrx
  integer          ::    mdimubx,      mtax,      mtbx,      mtcx,      mtnx
  integer          ::       mtrx,  ndimctbx,  ndimntbx,   ndimubx,save_mtbx
  integer          ::      neqtx,     nptot
  integer         ,allocatable ::    id1(:),   id2(:),   ii1(:),   ii2(:),   jd1(:)
  integer         ,allocatable ::    jd2(:),   jj1(:),   jj2(:),   kd1(:),   kd2(:)
  integer         ,allocatable ::    kk1(:),   kk2(:),lbdrat(:),   nnc(:),  nnfb(:)
  integer         ,allocatable ::    nnn(:),   npc(:),  npfb(:),   npn(:),npbrat(:),nbdrat(:)
end module maillage
!
module para_var
    use maillage
implicit none
  integer          ::    ccg2,   cfg2,   cng2,   ip00,   ip11
  integer          ::    ip12,   ip13,   ip21,   ip31,   ip40
  integer          ::    ip41,   ip42,   ip43,   ip44,   ip60
  integer          ::   kdimg,  kdimk,  kdimv,mdimtbf,mdimtcf
  integer          :: mdimtnf,mdimtrf, mdimub,ndimctf,ndimnts
  integer          :: ndimntu, ndimub,   nvar
  double precision :: ccg,cfg,cng
  parameter(ndimub =1200000)
  parameter(ndimctf=1000000)
  parameter(ndimnts=1000000)
  parameter(ndimntu=1)
  parameter(kdimg  =1)
  parameter(kdimv  =1)
  parameter(kdimk  =1)
  parameter(mdimub =4000)
  parameter(mdimtbf=9000)
  parameter(mdimtnf=9000)
  parameter(mdimtcf=16000)
  parameter(mdimtrf=1)
  parameter(nvar   =7)
!  parameter(    ccg=1./3.) !3D
!  parameter(    cng=1./3.) !3D
  parameter(    ccg=1./2.)
  parameter(    cng=1./2.)
  parameter(    cfg=1.)
!parameter(     ccg2=3) !3D
!parameter(     cng2=3) !3D
  parameter(   ccg2=2)
  parameter(   cng2=2)
  parameter(   cfg2=1)
end module para_var
!
module definition
  use para_fige
implicit none
  integer          :: klomg
  double precision ::              aa1,             anz,             dnz,             ha1,             omg
  double precision ::              pa1,           perio,             pnz,          protat,          ptrans
  double precision ::              rnz,            roa1,            ronz,             ta1,             tnz
  double precision :: varst(nsta,lsta)
end module definition
!
module chainecarac
implicit none
  character(len=24) :: c0,c1,c2,c3
  character(len=32) :: cb,cc,cd,ch,ci,cf,cm,cr,cs
  data c0/'VALEUR NON SIGNIFICATIVE'/
  data c1/'valeur par defaut       '/
  data c2/'                        '/
  data c3/'valeur anterieure       '/
  data cb/'                                '/
  data cc/' character trop long!           '/
  data cd/' pas de domaine deja cree!      '/
  data ch/' character attendu!             '/
  data ci/' entier attendu!                '/
  data cf/' pas de frontiere deja creee!   '/
  data cm/' liste d''entiers attendue!     '/
  data cr/' reel attendu!                  '/
  data cs/' suite de commande attendue!    '/
  character(len=7)  :: equat,equatt
  character(len=80) :: titrt1
  character(len=4)  :: config
end module chainecarac
!
module constantes
implicit none
    integer          :: intmx, linx
    double precision :: degrad,  pis2,raddeg,reelmn,reelmx
  data linx/132/
!  data intmx/999999/
!  data reelmx/999999999./
!  data reelmn/1.e-30/
!  data pis2/1.570796327/
!  data raddeg/57.29577951/
!  data degrad/0.01745329252/
end module constantes
!
module proprieteflu
implicit none
    double precision ::    cp,   cv,  gam, gam1, gam2
    double precision ::  gam3, gam4, gam5,pinfl,   pr
    double precision ::   prt,   ql,   rd,reynz,  rgp
end module proprieteflu
!
module kcle
  use para_fige
implicit none
    integer          ::              kanz,          kconfig,             kcte,          kdiscsv,             kdnz
    integer          ::           kdt1min,           kequat,             kgam,         kicychr0,         kicytur0
    integer          ::             kilim,         kischema,             kkcg,            kkdtl,         kkdualns
    integer          ::             kkexl,            kkfmg,           kklomg,           kkprec,           kkvisq
    integer          ::              kkvn,             klgx,             klzx,         kmdimtbx,         kmdimtcx
    integer          ::          kmdimtnx,         kmdimtrx,         kmdimubx,            kmtax,            kmtbx
    integer          ::             kmtcx,            kmtnx,            kmtrx,           kmuscl,             knba
    integer          ::          kncychro,         kncyexpl,         kncyresi,         kncysave
    integer          ::          kncyturb,        kndimctbx,        kndimntbx,         kndimubx,           kneqtx
    integer          ::              knfi,           kniter,           knitur,            knumt,             komg
    integer          ::            kperio,           kpinfl,             kpnz,              kpr,             kprt
    integer          ::               kql,              krd,           kreynz,             krnz,            kronz
    integer          ::           ktitrt1,             ktnz,             ktol,           ktolke,kvarst(nsta,lsta)
    integer          ::               kxk
    integer         ,allocatable ::     keta(:),    kki2(:),    kki4(:),    kkmf(:),   klmax(:),      kncycle(:)
end module kcle
!
module schemanum
  use para_fige
implicit none
    integer          ::    icychr0,      ilim,   ischema,       kcg,      kdtl
    integer          ::    kdualns,   kdualto,      kfmg,     klroe,     kprec
    integer          ::     ktrans,     kvisq,     lacou,       mgl,     muscl
    integer          ::    ncychro,   ncyexpl,   ncyresi,   ncysave
    integer          ::        nfi,     niter,     nitur,      numt
    double precision ::      al0,     be0,     cga,     cte,  dt1min
    double precision ::   epsroe,    freq,  resite,  reske1,reskeite
    double precision ::   resno1,   rki2t,   rki4t,     rm0,     tol
    double precision ::    tolke,  vrtalp,   vrtcz,  vrtlre,  vrtmac
    double precision ::   vrtxre,  vrtzre,      x0,      xk,      y0
    double precision ::       z0
    integer         ,allocatable ::  kmf(:),lmax(:),ncycle(:)
    double precision,allocatable :: eta(:),ki2(:),ki4(:)
   character(len=4) :: discsv
end module schemanum
!
module modeleturb
  use para_fige
implicit none
    integer          ::  icytur0, imxclko,  kclkep,   kcmut,  kcutke
    integer          ::   kditur,   keasm,  keinit,   kesst, kfludis
    integer          ::  kfracom,   kinke,  komsst,  kparoi,  ksecmb
    integer          ::  ktransi,   kutau,   kwsst,  lparoi,  naprng
    integer          ::   ncycke, ncycrac,ncytuke0, ncyturb
    double precision ::   alfae,  alfak, allfa0,allfae0,   beta
    double precision ::   beta1,  beta2,  betae,  betas, bethae
    double precision ::   bheta,    cb1,    cb2,    cc1,    cc2
    double precision ::    ccmu,    ce2,   ceta,    cgl,   cke1
    double precision ::    cke2,  cklb1,  ckle2,  cllog,    cmu
    double precision ::   cmukl,    ct1,    ct2,    ct3,    ct4
    double precision ::     cv1,    cw1,    cw2,    cw3, drtrac
    double precision ::    epse,   epsk, epspid,epstaud,epsvord
    double precision ::   kappa, okappa,rapvisq,  rbeta, rkplus
    double precision ::  roeinf, rokinf,  romeg,    rrk,  rtrac
    double precision ::    sige,   sigk,  sigma, sigma1, sigma2
    double precision ::  sigmab, sigmak, sigmal,  sigmd, sigme1
    double precision ::  sigme2, sigmeb,  sigmk,  sigmw,utaumin
    double precision ::    vkar,  wsig1,  wsig2, xkappa,    yp0
    double precision,allocatable :: pctvort(:)
end module modeleturb
!
module sortiefichier
  use para_fige
implicit none
    integer          ::    don1,    imp,  inia1,  inig1,   kdac
    integer          ::   kdacf,   kdav,   kdgc,  kdgcf,   kdgv
    integer          ::     kfa,    kfb,    kfc,    kfi,    kfn
    integer          ::     kfr,   kimp,   kres,  kvglo,    lec
    integer          :: lsortie,  nbfll,  nfreq,    out,    sec
    integer          ::    sor1,   sor2,   sor3,  sorf1,  sorf2
    double precision ::             alpha0,             beta0,            p0spi0,            q0spi0,              sref
    double precision :: temp_array(neqt,3),                v0,             xlref,              xref,              yref
    double precision ::               zref
    integer         ,allocatable :: nmfint(:)
  data kimp/3/
  data lec,imp,out,sec,sor1,sor2,kfa,sor3/11,12,13,14,15,16,17,18/
  data kdgv,kdav,kdgc,kdac/21,22,23,24/
!  data inig1,kfi,kfb,kfn,kfc,kfr,kdgcf,kdacf,kres/31,32,33,34,35,36,37,38,39/
  data don1,inia1,sorf1,sorf2/41,42,43,44/
end module sortiefichier
