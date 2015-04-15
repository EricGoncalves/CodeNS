module para_fige
  implicit none
  integer          ::   ista,    lg,lgcmdx,  lsta,    lt
  integer          ::     lz,   mtb,   mtt,  ndir,  neqt
  integer          ::   nind,   nmx,  nobj,  nsta
  parameter(  ndir=3     )
  parameter(  nind=3     )
  parameter(    lz=50    )
  parameter(    lg=6     )
  parameter(    lt=lz*lg )
  parameter(   mtb=600   )
  parameter(   mtt=mtb*lg)
  parameter(  neqt=7     )
  parameter(  nsta=50    )
  parameter(  lsta=7     )
  parameter(  ista=2     )
  parameter(  nobj=600   )
  parameter(   nmx=500   )
  parameter(lgcmdx=1316  )
end module para_fige
!
module para_var
implicit none
  integer          ::    ccg2,   cfg2,   cng2,   ip00,   ip11
  integer          ::    ip12,   ip13,   ip21,   ip31,   ip40
  integer          ::    ip41,   ip42,   ip43,   ip44,   ip60
  integer          ::   kdimg,  kdimk,  kdimv,mdimtbf,mdimtcf
  integer          :: mdimtnf,mdimtrf, mdimub,ndimctf,ndimnts
  integer          :: ndimntu, ndimub,   nvar
  double precision :: ccg,cfg,cng
  parameter(ndimub =120000)
  parameter(ndimctf=100000)
  parameter(ndimnts=100000)
  parameter(ndimntu=1)
  parameter(kdimg  =1)
  parameter(kdimv  =1)
  parameter(kdimk  =1)
  parameter(mdimub =400)
  parameter(mdimtbf=900)
  parameter(mdimtnf=900)
  parameter(mdimtcf=160)
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
  parameter(ip00=ndimub)
  parameter(ip11=ndimctf+kdimg*ndimctf/ccg2)
  parameter(ip12=kdimv*(ip11-1)+1)
  parameter(ip13=kdimk*(ip11-1)+1)
  parameter(ip21=ndimnts+kdimg*ndimnts/cng2+ndimntu)
  parameter(ip31=1+3*(ndimnts+kdimg*ndimnts/cng2))
  parameter(ip40=mdimub)
  parameter(ip41=mdimtbf+kdimg*mdimtbf/cfg2)
  parameter(ip42=mdimtbf+kdimg*mdimtbf/cfg2)
  parameter(ip43=mdimtbf+kdimg*mdimtbf/cfg2)
  parameter(ip44=mdimtbf+kdimg*mdimtbf/cfg2)
  parameter(ip60=nvar)
end module para_var
!
module boundary
  use para_fige
implicit none
  integer          :: crotr(mtb),imaxb(mtt),iminb(mtt),jmaxb(mtt),jminb(mtt)
  integer          ::       kexl,kmaxb(mtt),kminb(mtt),  lbd(mtt),lbdko(mtt)
  integer          ::  mdnc(mtt),  mmb(mtt),  mpb(mtt),  mpc(mtt), mper(mtt)
  integer          ::   mpn(mtt),  mpr(mtt),  nba(mtb),       nbd, nbdc(mtb)
  integer          ::      nbdko, ndcc(mtb), ndlb(mtb), ndrr(mtb), nfba(mtb)
  integer          ::  nfbc(mtb), nfbn(mtb), nfbr(mtb), nfei(mtb),srotr(mtb)
  double precision :: bc(mtb,ista*lsta)
  character(len=4) :: cl(mtb)
  character(len=2) :: indfl(mtb)
end module boundary
!
module maillage
  use para_fige
implicit none
  integer          ::     id1(lt),    id2(lt),    ii1(lt),    ii2(lt),    jd1(lt)
  integer          ::     jd2(lt),    jj1(lt),    jj2(lt),    kcaldis,    kd1(lt)
  integer          ::     kd2(lt),    kecrdis,    kk1(lt),    kk2(lt),    klecdis
  integer          ::         kvn,lbdrat(mtb),        lgx,        lzx,    mdimtbx
  integer          ::     mdimtcx,    mdimtnx,    mdimtrx,    mdimubx,       mtax
  integer          ::        mtbx,       mtcx,       mtnx,       mtrx, nbdrat(lz)
  integer          ::    ndimctbx,   ndimntbx,    ndimubx,      neqtx,    nnc(lt)
  integer          ::    nnfb(lt),    nnn(lt), npbrat(lz),    npc(lt),   npfb(lt)
  integer          ::     npn(lt),      nptot
end module maillage
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
  integer          ::           kdt1min,           kequat,         keta(lt),             kgam,         kicychr0
  integer          ::          kicytur0,            kilim,         kischema,             kkcg,            kkdtl
  integer          ::          kkdualns,            kkexl,            kkfmg,         kki2(lt),         kki4(lt)
  integer          ::            kklomg,         kkmf(lt),           kkprec,           kkvisq,             kkvn
  integer          ::              klgx,        klmax(lt),             klzx,         kmdimtbx,         kmdimtcx
  integer          ::          kmdimtnx,         kmdimtrx,         kmdimubx,            kmtax,            kmtbx
  integer          ::             kmtcx,            kmtnx,            kmtrx,           kmuscl,             knba
  integer          ::          kncychro,      kncycle(lg),         kncyexpl,         kncyresi,         kncysave
  integer          ::          kncyturb,        kndimctbx,        kndimntbx,         kndimubx,           kneqtx
  integer          ::              knfi,           kniter,           knitur,            knumt,             komg
  integer          ::      kpctvort(lt),           kperio,           kpinfl,             kpnz,              kpr
  integer          ::              kprt,              kql,              krd,           kreynz,             krnz
  integer          ::             kronz,          ktitrt1,             ktnz,             ktol,           ktolke
  integer          :: kvarst(nsta,lsta),              kxk
end module kcle
!
module schemanum
  use para_fige
implicit none
  integer          ::    icychr0,      ilim,   ischema,       kcg,      kdtl
  integer          ::    kdualns,   kdualto,      kfmg,     klroe,   kmf(lt)
  integer          ::      kprec,    ktrans,     kvisq,     lacou,  lmax(lt)
  integer          ::        mgl,     muscl,   ncychro,ncycle(lg),   ncyexpl
  integer          ::    ncyresi,   ncysave,       nfi,     niter,     nitur
  integer          ::       numt
  double precision ::      al0,     be0,     cga,     cte,  dt1min
  double precision ::   epsroe, eta(lt),    freq, ki2(lt), ki4(lt)
  double precision ::   resite,  reske1,reskeite,  resno1,   rki2t
  double precision ::    rki4t,     rm0,     tol,   tolke,  vrtalp
  double precision ::    vrtcz,  vrtlre,  vrtmac,  vrtxre,  vrtzre
  double precision ::       x0,      xk,      y0,      z0
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
  double precision ::       alfae,      alfak,     allfa0,    allfae0,       beta
  double precision ::       beta1,      beta2,      betae,      betas,     bethae
  double precision ::       bheta,        cb1,        cb2,        cc1,        cc2
  double precision ::        ccmu,        ce2,       ceta,        cgl,       cke1
  double precision ::        cke2,      cklb1,      ckle2,      cllog,        cmu
  double precision ::       cmukl,        ct1,        ct2,        ct3,        ct4
  double precision ::         cv1,        cw1,        cw2,        cw3,     drtrac
  double precision ::        epse,       epsk,     epspid,    epstaud,    epsvord
  double precision ::       kappa,     okappa,pctvort(lt),    rapvisq,      rbeta
  double precision ::      rkplus,     roeinf,     rokinf,      romeg,        rrk
  double precision ::       rtrac,       sige,       sigk,      sigma,     sigma1
  double precision ::      sigma2,     sigmab,     sigmak,     sigmal,      sigmd
  double precision ::      sigme1,     sigme2,     sigmeb,      sigmk,      sigmw
  double precision ::     utaumin,       vkar,      wsig1,      wsig2,     xkappa
  double precision ::         yp0
end module modeleturb
!
module sortiefichier
  use para_fige
implicit none
  integer          ::        don1,        imp,      inia1,      inig1,       kdac
  integer          ::       kdacf,       kdav,       kdgc,      kdgcf,       kdgv
  integer          ::         kfa,        kfb,        kfc,        kfi,        kfn
  integer          ::         kfr,       kimp,       kres,      kvglo,        lec
  integer          ::     lsortie,      nbfll,      nfreq,nmfint(mtb),        out
  integer          ::         sec,       sor1,       sor2,       sor3,      sorf1
  integer          ::       sorf2
  double precision :: alpha0, beta0,p0spi0,q0spi0,  sref,temp_array(neqt,3)
  double precision ::     v0, xlref,  xref,  yref,  zref
  data kimp/3/
  data lec,imp,out,sec,sor1,sor2,kfa,sor3/11,12,13,14,15,16,17,18/
  data kdgv,kdav,kdgc,kdac/21,22,23,24/
!  data inig1,kfi,kfb,kfn,kfc,kfr,kdgcf,kdacf,kres/31,32,33,34,35,36,37,38,39/
  data don1,inia1,sorf1,sorf2/41,42,43,44/
end module sortiefichier
