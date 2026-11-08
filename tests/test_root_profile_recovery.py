"""Exercise production root-domain and recovery-demand procedures in compiled Fortran."""
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
BUILD=Path(os.environ.get('ECOSIM_TEST_BUILD',ROOT/'build/Darwin-arm64-double-Release'))

class RootProfileRecoveryTest(unittest.TestCase):
    def test_domain_cap_and_living_root_recovery(self):
        source=(ROOT/'f90src/Plant_bgc/RootMod.F90').read_text()
        names=['LimitPrimaryElongationResp','RootRecoverySinks','GetInLayerRootDepth','StructuralRootElongation','PrimeRootsExtension','RelocateRootBoundaryReserves']
        owners={'plt_site':'PlantSite','plt_morph':'PlantMorphology','plt_biom':'PlantBiomass','plt_allom':'PlantAllometry','plt_pheno':'PlantPhenology','plt_ew':'PlantEnergyWater','plt_rbgc':'PlantRootBGC','plt_bgcr':'PlantBGCRates','plt_photo':'PlantPhotosynthesis','plt_soilchem':'PlantSoilChemistry'}
        fields=set()
        for name in names:
            start=source.index('  subroutine '+name+'(');end=source.index('  end subroutine '+name,start)
            fields.update(re.findall(r'(plt_\w+)%([\w]+)',source[start:end]))
        allocation=[]
        for owner,field in sorted(fields):
            api=(ROOT/'f90src/APIData'/(owners[owner]+'APIData.F90')).read_text()
            decl=re.search(r'::\s*'+field+r'\s*(\([^\n]*?\))?',api,re.I);self.assertIsNotNone(decl,field)
            rank=(decl[1] or '').count(':')
            if not rank:continue
            shape={1:'1',2:'3,1',3:'3,1,1',4:'3,3,1,1',5:'3,2,3,1,1'}[rank]
            if field=='RootMycoNonstElms_rpvr':shape='3,2,3,1'
            if field in ['Root1stSpecLen_pft']:shape='2,1'
            if field in ['Root2ndRadius_rpvr','Root2ndEffLen4uptk_rpvr']:shape='2,3,1'
            if field=='FracRootElmAllocm':shape='3,2'
            if field=='Root1stRadius_pvr':shape='2,3,1'
            if owner in ['plt_rbgc','plt_ew','plt_bgcr'] and rank==3:shape='2,3,1'
            if field=='RootMyco1stSinkC_rpvr':shape='3,1,1'
            if owner=='plt_photo' and rank==2:shape='1,1'
            if owner in ['plt_site','plt_soilchem'] and rank==1:shape='3'
            if field=='CumSoilThickness_vr':shape='0:3'
            allocation += [f'allocate({owner}%{field}({shape}))',f'{owner}%{field}=0']
        driver='''program test_profile_recovery
use data_kind_mod,only:r8=>DAT_KIND_R8,yearIJ_type
use RootMod,only:LimitPrimaryElongationResp,RootRecoverySinks,GetInLayerRootDepth,StructuralRootElongation,RelocateRootBoundaryReserves
use PlantAPICommonData
use PlantBGCPars
use ElmIDMod
use EcoSIMCtrlMod,only:lcoarseroot
'''+ '\n'.join('use '+v+'APIData,only:'+k for k,v in owners.items())+'''
implicit none
real(r8)::u,v,fine,support,saved(3),loc(1,3),effective(1,3),expected
logical::changed
type(yearIJ_type)::yearIJ
real(r8)::before(3),after(3),litter(3),resp,sinks(2,3),tips(1),temps(3)
allocate(JZ1,MaxNumRootAxes);JZ1=3;MaxNumRootAxes=1
pltpar%jroots=2;pltpar%MaxNumRootAxes=1;pltpar%k_fine_comp=1
'''+ '\n'.join(allocation)+'''
plt_site%NK=3;plt_site%NU=1;plt_site%MaxNumRootLays=2
plt_site%CumSoilThickness_vr=[0._r8,.1_r8,.3_r8,.6_r8]
plt_site%DLYR3=[.1_r8,.2_r8,.3_r8]
plt_morph%NGTopRootLayer_pft=1;plt_morph%SeedDepth_pft=.01_r8
plt_morph%NumStructuralRootAxes_pft=1;plt_morph%Myco_pft=2
plt_morph%MaxSoilLays4Root_pft=3
plt_morph%Root1stDepz_raxes=.35_r8;plt_morph%Root1stLenPP_rpvr=.25_r8
plt_biom%RootMycoNonstElms_rpvr(:,1,3,1)=[1._r8,.1_r8,.01_r8]
plt_biom%RootMycoNonstElms_rpvr(:,2,3,1)=[2._r8,.2_r8,.02_r8]
saved=sum(sum(plt_biom%RootMycoNonstElms_rpvr(:,:,:,1),dim=3),dim=2)
RTDPX=0._r8
call RelocateRootBoundaryReserves(1,changed)
if(.not.changed)stop 9
call GetInLayerRootDepth(1,loc,effective)
if(plt_morph%Root1stDepz_raxes(1,1)/=.3_r8)stop 10
if(plt_morph%NRoot1stTipLay_raxes(1,1)/=2.or.plt_morph%MaxSoilLays4Root_pft(1)/=2)stop 11
if(plt_morph%Root1stLenPP_rpvr(2,1,1)/=.2_r8)stop 12
if(any(plt_biom%RootMycoNonstElms_rpvr(:,:,3,1)/=0._r8))stop 13
if(any(abs(sum(sum(plt_biom%RootMycoNonstElms_rpvr(:,:,:,1),dim=3),dim=2)-saved)>1.e-14_r8))stop 14
if(loc(1,3)/=0._r8)stop 15
! Remaining length bounds growth respiration before any C/N/P debit.
plt_morph%Root1stSpecLen_pft=10._r8
plt_allom%FracRootElmAllocm=.5_r8
plt_morph%Num1stAxesPerStructRootXPOP_pft=.2_r8
plt_morph%Root1stDepz_raxes=.299_r8
u=1._r8;v=.5_r8
call LimitPrimaryElongationResp(1,1,2,1,.2_r8,.8_r8,u,v)
expected=(.3_r8-.299_r8)*.2_r8/(10._r8*.5_r8)*.2_r8/.8_r8
if(abs(u-expected)>1.e-18_r8.or.abs(v-expected)>1.e-18_r8)stop 20
! At boundary no elongation carbon is charged; zero resource availability stays zero.
plt_morph%Root1stDepz_raxes=.3_r8;u=1._r8;v=0._r8
call LimitPrimaryElongationResp(1,1,2,1,.2_r8,.8_r8,u,v)
if(u/=0._r8.or.v/=0._r8)stop 21
! No artificial axes or tissue: only living primary segments acquire basal demands.
plt_pheno%iPlantRootProfile_pft=3;plt_pheno%iPlant2ndGrothPattern_pft=1
plt_biom%Root1stActStructElms_rpvr=0._r8
plt_morph%Root1stLenPP_rpvr=.1_r8
plt_morph%FineRootBranchFreq_pft=100._r8
plt_morph%Root2ndRadius_rpvr=.00025_r8
plt_morph%Root2ndEffLen4uptk_rpvr=0._r8;Root2ndTipLen4uptk=.002_r8
lcoarseroot=.true.;fine=0._r8
call RootRecoverySinks(1,1,1,1.e-7_r8,fine,support)
if(fine/=0._r8.or.support/=0._r8)stop 30
plt_biom%Root1stActStructElms_rpvr(1,1,1,1)=1.e-12_r8
call RootRecoverySinks(1,1,1,1.e-7_r8,fine,support)
if(fine<=0._r8.or.support<=0._r8)stop 31
! Existing conductance is retained if larger than the bud demand.
fine=1._r8
call RootRecoverySinks(1,1,1,1.e-7_r8,fine,support)
if(fine/=1._r8)stop 32
lcoarseroot=.false.;fine=0._r8
call RootRecoverySinks(1,1,1,1.e-7_r8,fine,support)
if(fine/=0._r8.or.support/=0._r8)stop 33
! Execute the actual positive-growth procedure through the reserve debit and
! geometry update. The cap must stop at the boundary without losing C/N/P.
lcoarseroot=.true.
plt_morph%Root1stDepz_raxes=.299_r8
plt_morph%Root1stAxesTipDepz2Surf_pft=.299_r8
plt_morph%Root1stLenPP_rpvr=.1_r8
plt_morph%Root1stLenPP_rpvr(2,1,1)=.199_r8
plt_morph%NRoot1stTipLay_raxes=2
plt_morph%Root1stMaxRadius1_pft=.001_r8
plt_morph%Root1stMaxRadius_pft=.001_r8
plt_morph%MainBranchNum_pft=1
plt_morph%Root1stSpecLen_pft=10._r8
plt_morph%Num1stAxesPerStructRootXPOP_pft=.2_r8
plt_biom%ZERO4Groth_pft=1.e-30_r8
plt_biom%RootMycoNonstElms_rpvr=0._r8
plt_biom%RootMycoNonstElms_rpvr(:,1,2,1)=[1._r8,.1_r8,.01_r8]
plt_biom%Root1stActStructElms_rpvr=0._r8
plt_biom%Root1stActStructElms_rpvr(:,2,1,1)=[1.e-6_r8,5.e-9_r8,5.e-10_r8]
plt_biom%RootMyco1stStrutElms_rpvr=plt_biom%Root1stActStructElms_rpvr
plt_biom%RootMyco1stElm_raxs(:,1,1)=plt_biom%Root1stActStructElms_rpvr(:,2,1,1)
plt_allom%RootBiomGrosYld_pft=.8_r8
plt_allom%CNRTS_pft=.01_r8;plt_allom%CPRTS_pft=.001_r8
plt_allom%FracRootElmAllocm=.5_r8
plt_rbgc%RAutoRootO2Limter_rpvr=1._r8
plt_site%ZERO=1.e-30_r8
RmSpecPlant=.01_r8;EMODR=50._r8
plt_ew%PSIRoot_pvr=-.1_r8
sinks=1._r8;tips=1._r8;temps=1._r8;litter=0._r8;resp=0._r8
before=plt_biom%RootMycoNonstElms_rpvr(:,1,2,1)+plt_biom%RootMyco1stStrutElms_rpvr(:,2,1,1)
call StructuralRootElongation(yearIJ,1,1,2,2,1,.5_r8,tips,sinks,1._r8,temps,1._r8,1._r8,.2_r8,.005_r8,.0005_r8,litter,resp)
after=plt_biom%RootMycoNonstElms_rpvr(:,1,2,1)+plt_biom%RootMyco1stStrutElms_rpvr(:,2,1,1)+litter
after(1)=after(1)-resp
if(any(abs(after-before)>1.e-12_r8))stop 40
if(abs(plt_morph%Root1stDepz_raxes(1,1)-.3_r8)>1.e-14_r8)stop 41
if(any(plt_biom%RootMycoNonstElms_rpvr(:,1,2,1)<0._r8))stop 42
print *,'PASS: root domain, conservative reserve relocation, pre-debit elongation cap, living-only recovery demand'
end program
'''
        source=source.replace('  public :: RootBGCModel','  public :: RootBGCModel, '+', '.join(names),1)
        flags=(BUILD/'f90src/Plant_bgc/CMakeFiles/plant_bgc.dir/flags.make').read_text()
        inc=shlex.split(re.search(r'^Fortran_INCLUDES = (.*)$',flags,re.M)[1])+['-I'+str(BUILD/'f90src/Plant_bgc')]
        opts=shlex.split(re.search(r'^Fortran_FLAGS = (.*)$',flags,re.M)[1])
        link=shlex.split((BUILD/'drivers/ecosim/CMakeFiles/ecosim.f90.x.dir/link.txt').read_text())
        if '--' in link:link=link[link.index('--')+1:]
        cleaned=[];i=0
        while i<len(link):
            if link[i]=='-o':i+=2;continue
            if not link[i].endswith('.o'):cleaned.append(link[i])
            i+=1
        with tempfile.TemporaryDirectory(prefix='ecosim_profile_') as work:
            w=Path(work);(w/'RootMod.F90').write_text(source);(w/'test.F90').write_text(driver)
            for cmd in [cleaned[:1]+inc+opts+['-J'+work,'-c',str(w/'RootMod.F90'),'-o',str(w/'RootMod.o')],cleaned[:1]+['-I'+work]+inc+opts+[str(w/'test.F90'),str(w/'RootMod.o'),'-o',str(w/'test.x')]+cleaned[1:],[str(w/'test.x')]]:
                p=subprocess.run(cmd,cwd=BUILD/'drivers/ecosim',capture_output=True,text=True)
                self.assertEqual(p.returncode,0,p.stdout+p.stderr)
            self.assertIn('PASS:',p.stdout)

if __name__=='__main__':unittest.main()
