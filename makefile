SHELL=/bin/sh

SRCS    =                modules.f90             at_dlist.f90         \
at_ecrdist.f90           at_fidist.f90           at_grdist.f90        \
at_lecdist.f90           at_lecopt.f90           atccc.f90            \
atccfp.f90               atctranske.f90          atdist_1.f90         \
atdist_2.f90             atdist_3.f90            atecrfp.f90          \
atindnor.f90             atcaldis.f90            atcaldis3.f90        \
atintrans.f90            atlecdon.f90            atparoi.f90          \
convich.f90              b1_cpbd.f90             b1_cpfw.f90          \
b1_crbds.f90             b1_crdms.f90            b1_dffw.f90          \
b1_dfgm.f90              b1_dfnm.f90             b1_dfnzst.f90        \
b1_dfph.f90              b1_dfpmcfg.f90          b1_dfpmdsd.f90       \
b1_dfpmdtd.f90           b1_dfpmdtg.f90          b1_dfpmimd.f90       \
b1_dfpmtbkeg.f90         b1_dfpmtbn.f90          b1_dfst.f90          \
b1_dfst0.f90             b1_dftl1.f90            b1_dpbd.f90          \
b1_end.f90               b1_inbdb.f90            b1_inbdc.f90         \
b1_inbdn.f90             b1_infw.f90             b1_ingr.f90          \
b1_intn.f90              b1_secpfw.f90           b1_svfw.f90          \
b1_svgr.f90              b2_crbds.f90            b2_crdms.f90         \
b2_dpdim.f90             b3_crdms.f90            cctcmd.f90           \
synterr.f90              entier.f90              valenti.f90          \
reel.f90                 valreel.f90             lentier.f90          \
vallent.f90              tcmd_cpbd.f90           tcmd_cpfw.f90        \
tcmd_crbds.f90           tcmd_crdms.f90          tcmd_dffw.f90        \
tcmd_dfgm.f90            tcmd_dfnm.f90           tcmd_dfnzst.f90      \
tcmd_dfph.f90            tcmd_dfpmcfg.f90        tcmd_dfpmdsd.f90     \
tcmd_dfpmdtd.f90         tcmd_dfpmdtg.f90        tcmd_dfpmimd.f90     \
tcmd_dfpmtbkeg.f90       tcmd_dfpmtbn.f90        tcmd_dfst.f90        \
tcmd_dftl1.f90           tcmd_dpbd.f90           tcmd_inbdb.f90       \
tcmd_inbdc.f90           tcmd_inbdn.f90          tcmd_infw.f90        \
tcmd_ingr.f90            tcmd_intn.f90           tcmd_secpfw.f90      \
tcmd_svfw.f90            tcmd_svgr.f90           rfve.f90             \
rfvc.f90                 rfvr.f90                utidd.f90            \
rbse.f90                 rbtc.f90                rbte.f90             \
rbtr.f90                 rbvc.f90                rbve.f90             \
clidd.f90                clchoc.f90              cldebit.f90          \
cldebit_prcd.f90         clextr.f90              clgli2.f90           \
idirch.f90               clglis.f90              clidd2_prcd.f90      \
clidi.f90                clidi0.f90              lparoi1.f90          \
lparoi2.f90              lparoi4.f90             lparoi3d.f90         \
cllparoi1.f90            cllparoi2.f90           clnrd.f90            \
clnrd_prcd.f90           clpara.f90              clpari.f90           \
clprd.f90                clprd_prcd.f90          clsym.f90            \
clvrt.f90                utidi.f90               utitfr.f90           \
utinia.f90               utnrd.f90               utpari.f90           \
utprd.f90                utreadav.f90            utsor.f90            \
utsorfr.f90              utvrt.f90               utdebit.f90          \
utdon_gen.f90            rbvr.f90                lecture_acou.f90     \
rbord.f90                gtcmd.f90               splcmd.f90           \
rdcmd.f90                rfsc.f90                rfspstc.f90          \
rfspstf.f90              rfsr.f90                residu.f90           \
writda.f90               writdg.f90              cvccg.f90            \
rscpsv.f90               rbc.f90                 cpbd.f90             \
c_cpbd.f90               extmhg.f90              norm.f90             \
snorm.f90                svol.f90                vervol.f90           \
metrics.f90              metric.f90              metric2.f90          \
metric3.f90              teq_grads.f90           teq_gradv.f90        \
met_brkl.f90             met_brko.f90            met_chmut.f90        \
met_cutke.f90            met_cutked.f90          met_cutsa.f90        \
met_dist.f90             met_difsa.f90           met_dual.f90         \
met_dual2.f90            met_fludc.f90           met_fludcsa.f90      \
met_fludko.f90           met_fludmt.f90          met_gradtr.f90       \
met_klmut.f90            met_inikl.f90           met_komut.f90        \
met_kocmut.f90           met_iniko.f90           met_inisa.f90        \
zvismo.f90               zgrad.f90               zgrad2.f90           \
zfluto.f90               met_iniuttau.f90        met_intep3.f90       \
met_bare.f90             met_bark.f90            met_rbve.f90         \
met_rbvr.f90             met_rbsr.f90            met_rbvc.f90         \
met_rfvc.f90             met_rfve.f90            met_bord.f90         \
met_brad.f90             met_cut.f90             met_exgr.f90         \
met_ini.f90              met_inmut.f90           met_prod.f90         \
met_laplaciens.f90       met_rbsc.f90            met_rbse.f90         \
met_roe.f90              met_roe2o.f90           met_roe2oh.f90       \
met_smkes.f90            met_smker.f90           met_smkec.f90        \
met_smkesas.f90          met_klnmut.f90          met_komutr.f90       \
met_kobmut.f90           met_klsmut.f90          met_kemut.f90        \
met_ke2mut.f90           met_kemutm.f90          met_kemutr.f90       \
met_pardis.f90           met_parko.f90           met_mtcorf1.f90      \
met_samut.f90            met_smch.f90         \
met_smdes.f90            met_smke.f90            met_smkl.f90         \
met_smklsas.f90          met_smko.f90            met_smkor.f90        \
met_smmt.f90             met_smmtr.f90           met_smsa.f90         \
met_uttau.f90            pgrad.f90               lpke1.f90            \
lpke2.f90                lpke.f90                lpkl1.f90            \
lpkl2.f90                lpkl.f90                lpkomega1.f90        \
lpkomega2.f90            lpkomega.f90            lpsa1.f90            \
lpsa2.f90                lpsa.f90                lpker1.f90           \
lpker.f90                lpkomegar1.f90          lpkomegar.f90        \
lp2kl1.f90               lp2kl3d.f90             lp2kl.f90            \
lp2ke1.f90               lp2ke.f90               lp2kw1.f90           \
lp2kw.f90                lp2sa1.f90              lp2sa.f90            \
met_klrmut.f90           met_smsasas.f90         utinig.f90           \
met_yplus.f90            at_cutke.f90            sch_acou.f90         \
sch_ausmp.f90            sch_ausmp_prcd.f90      sch_dual.f90         \
sch_dual2.f90            sch_duin.f90            sch_duup.f90         \
sch_duup2.f90            sch_jameson.f90         sch_jameson_pond.f90 \
sch_jameson3.f90         sch_jameson3pond.f90    sch_jam3_turb.f90    \
sch_roe.f90              sch_roe_prcd.f90        sch_roe_pond.f90     \
sch_roe_pond_prcd.f90    sch_roe_euler.f90       sch_ausmp_pond.f90   \
sch_turb.f90             sch_turb_pond.f90       sch_hllc.f90         \
sch_hllc_prcd.f90        sch_hllc_euler.f90      sch_rusanov.f90      \
sch_rusanov_prcd.f90     sch_weno3.f90           sch_weno3split.f90   \
sch_weno3split2.f90      sch_weno3_3d.f90        sch_weno3pond.f90    \
sch_weno3pond_split2.f90 sch_weno3pond_3d.f90    sch_weno5.f90        \
sch_weno5_3d.f90         sch_weno5pond.f90       smg_cf_2d.f90        \
smg_cf_3d.f90            smg_cf.f90              smg_cn.f90           \
smg_fcm.f90              smg_fcr.f90             smg_fcs.f90          \
smg_fcv.f90              zpres.f90               dissip_jameson.f90   \
dissip_jameson_prcd2.f90 dissip_jameson_turb.f90 zvisqc.f90           \
smg_res.f90              smg_flu.f90             smg_upc.f90          \
chronos.f90              chronos_prcd.f90        chrono.f90           \
impli2_eqt.f90           impli2_eqt_3d.f90       implimf_2d.f90       \
implimf_prcd2.f90        implimf_3d.f90          implimf_eu.f90       \
met_num.f90              prcd_turkel.f90         sch_expli.f90        \
atsch_num.f90            dua_resro.f90           sortieplot.f90       \
sortietest.f90           sortieplot2.f90         sortieplot3.f90      \
cccca.f90                cccva.f90               chdacc.f90           \
svfw.f90                 smg_num.f90             cpfw.f90             \
c_cpfw.f90               defcpbd.f90             defdffw.f90          \
defdfgm.f90              defdfnm.f90             defdfnzst.f90        \
defdfph.f90              defdfpmcfg.f90          defdfpmdsd.f90       \
defdfpmdtd.f90           defdfpmdtg.f90          defdfpmimd.f90       \
defdfpmtbn.f90           defdfst.f90             defdftl1.f90         \
defintn.f90              defsecpfw.f90           def.f90              \
inivec.f90               inimem.f90              initbs.f90           \
initcs.f90               initis.f90              initns.f90           \
crbds.f90                cvcvg.f90               chdgcv.f90           \
crdms.f90                dffw.f90                dfgm.f90             \
dfnzst.f90               dfph.f90                dfpmcfg.f90          \
dpbdi.f90                dpbdb.f90               dpbdc.f90            \
dpbdn.f90                dpbd.f90                eend.f90             \
readda.f90               readdg.f90              inbdbad.f90          \
inbdbdf.f90              inbdbfl.f90             inbdbst.f90          \
inbdb.f90                c_inbdb.f90             inbdc.f90            \
inbdn.f90                infw.f90                ingr.f90             \
svdual.f90               svgr.f90                c_crbds.f90          \
c_crdms.f90              c_dffw.f90              c_dfgm.f90           \
c_dfnm.f90               c_dfnzst.f90            c_dfph.f90           \
c_dfpmcfg.f90            c_dfpmdsd.f90           c_dfpmdtd.f90        \
c_dfpmdtg.f90            c_dfpmimd.f90           c_dfpmtbkeg.f90      \
c_dfpmtbn.f90            c_dfst.f90              c_dfst0.f90          \
c_dftl1.f90              c_dpbd.f90              c_dpdim.f90          \
c_end.f90                c_inbdc.f90             c_inbdn.f90          \
c_infw.f90               c_ingr.f90              c_intn.f90           \
c_nzst.f90               c_secpfw.f90            c_svbd.f90           \
c_svbdb.f90              c_svbdc.f90             c_svbdn.f90          \
c_svfw.f90               c_svgr.f90              partition.f90        \
solve.f90 

OBJS =  ${SRCS:.f90=.o}

# Tunable parameters
#
# CF		Name of the fortran compiling system to use
# LDFLAGS	Flags to the loader
# LIBS		List of libraries
# CMD		Name of the executable
# PROFLIB	Library needed for profiling
#
#FC =  gfortran
FC =  ifort
#CMD =	 solver_air
CMD =	 solver_tic

# To perform the default compilation, use the first line
#FFLAGS =  -O2  -ffree-line-length-none
#LDFLAGS = -O2 -fdefault-double-8 -fdefault-real-8
#FFLAGS =  -O2   -fdefault-double-8 -fdefault-real-8 -ffree-line-length-none
FFLAGS =  -O2 -openmp -parallel -threads
LDFLAGS = -O2 -openmp -parallel -threads


# Lines from here on down should not need to be changed.
#
%.o: %.f90
	$(FC) $(FFLAGS) -c $<
all:		$(CMD)

$(CMD):		$(OBJS)
	  $(FC) $(FFLAGS) -o $(@) $(OBJS) $(LIBS)

clean:	
	rm -f *.o *.mod *.optrpt *__genmod.f90
