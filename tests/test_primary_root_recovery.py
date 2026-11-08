"""Exercise the production woody-root recovery and maintenance routines.

Private routines are exposed only in a temporary copy of RootMod. Existing build
compiler flags and libraries are reused; no production procedure is substituted.
"""
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BUILD = Path(os.environ.get('ECOSIM_TEST_BUILD', ROOT / 'build/Darwin-arm64-double-Release'))


class PrimaryRootRecoveryTest(unittest.TestCase):
    def test_surviving_segments_resources_and_maintenance(self):
        source = (ROOT / 'f90src/Plant_bgc/RootMod.F90').read_text()
        owners = {'plt_site':'PlantSite', 'plt_morph':'PlantMorphology',
                  'plt_biom':'PlantBiomass', 'plt_allom':'PlantAllometry',
                  'plt_pheno':'PlantPhenology', 'plt_photo':'PlantPhotosynthesis',
                  'plt_ew':'PlantEnergyWater', 'plt_rbgc':'PlantRootBGC',
                  'plt_bgcr':'PlantBGCRates', 'plt_soilchem':'PlantSoilChemistry'}
        associations = {}
        for procedure in ['StructuralRootThickening','ReserveMediumRootMaintenance',
                          'RootThickenWaterFunc','CRootsLive2DeadTransformation']:
            start = source.lower().index('  subroutine '+procedure.lower()+'(')
            end = source.lower().index('  end subroutine '+procedure.lower(), start)
            for alias,owner,field in re.findall(r'(\w+)\s*=>\s*(plt_\w+)%([\w]+)',source[start:end]):
                if not any(o==owner and f.lower()==field.lower() for o,f in associations):
                    associations[(owner,field)] = alias
        allocation=[]
        for (owner,field),alias in associations.items():
            api=(ROOT/'f90src/APIData'/ (owners[owner]+'APIData.F90')).read_text()
            decl=re.search(r'::\s*'+field+r'\s*(\([^\n]*?\))?',api,re.I)
            self.assertIsNotNone(decl,field)
            rank=(decl[1] or '').count(':')
            if rank==0: continue
            if rank==1: shape='2' if owner in ['plt_site','plt_soilchem'] else '1'
            elif owner=='plt_biom':
                shape={2:'3,1',3:'3,1,1',4:'3,2,1,1',5:'3,2,2,1,1'}[rank]
                if field=='RootMycoNonstElms_rpvr':shape='3,2,2,1'
            elif owner=='plt_allom':shape='3,1'
            elif owner in ['plt_ew','plt_rbgc','plt_bgcr']:shape={2:'2,1',3:'2,2,1'}[rank]
            else:shape={2:'2,1',3:'2,1,1'}[rank]
            if field=='Root1stRadius_pvr':shape='2,2,1'
            if field=='CumSoilThickness_vr':shape='0:2'
            allocation += [f'allocate({owner}%{field}({shape}))', f'{owner}%{field}=0']
        aliases=', &\n'.join(f'{a}=>{o}%{f}' for (o,f),a in associations.items())
        driver='''program recovery_test
use data_kind_mod,only:r8=>DAT_KIND_R8,yearIJ_type
use RootMod,only:StructuralRootThickening,ReserveMediumRootMaintenance
use PlantAPICommonData
use PlantBGCPars
use EcoSIMCtrlMod,only:lcoarseroot
use ElmIDMod
use EcosimConst
'''+ '\n'.join('use '+v+'APIData,only:'+k for k,v in owners.items())+'''
implicit none
type(yearIJ_type)::yearIJ
real(r8)::sink(2,2),primarysink(2,1),temp(2),litter(3),resp,before(3),after(3),oldc,growth
real(r8)::tips(1),unlim(1),limited(1),expected
logical::found(2,1),done(1)
integer::testcase
allocate(JZ1);JZ1=2
pltpar%jroots=2;pltpar%MaxNumRootAxes=1
'''+ '\n'.join(allocation)+'''
associate(&
'''+aliases+''')
lcoarseroot=.true.
RmSpecPlant=.01_r8;EMODR=50._r8;TurgPSIMin4OrganExtens=.1_r8
StalkBiomGrowthYld_pft=.8_r8
rECLiveCRoot_pft(:,1)=[1._r8,.005_r8,.0005_r8]
rECDeadCRoot_pft(:,1)=rECLiveCRoot_pft(:,1)
CRootActVolPerMassC_pft=5.e-6_r8
Root1stLenPP_rpvr=.1_r8;Num1stAxesPerStructRootXPOP_pft=.01_r8
Root1stRadius_rpvr=.0001_r8
PSIRootTurg_vr=.5_r8;SoilModulus4RootRadialexp_vr=0._r8
RAutoRootO2Limter_rpvr=1._r8;ZERO4Groth_pft=1.e-30_r8
sink=1._r8;primarysink=.5_r8;temp=1._r8
! Existing structural N is below the former cutoff. Recovery must use reserves.
do testcase=1,9
 Root1stActStructElms_rpvr(:,1,1,1)=[1.e-6_r8,5.e-9_r8,5.e-10_r8]
 Root1stLigStructElms_rpvr=0._r8
 RootMyco1stStrutElms_rpvr=Root1stActStructElms_rpvr
 RootMyco1stElm_raxs(:,1,1)=Root1stActStructElms_rpvr(:,1,1,1)
 RootMycoNonstElms_rpvr(:,1,1,1)=[1._r8,.1_r8,.01_r8]
 Root1stLenPP_rpvr=.1_r8;Num1stAxesPerStructRootXPOP_pft=.01_r8
 Root1stRadius_rpvr=.0001_r8
 if(testcase==2)RootMycoNonstElms_rpvr(ielmn,1,1,1)=0._r8
 if(testcase==3)RootMycoNonstElms_rpvr(ielmp,1,1,1)=0._r8
 if(testcase==4)then
  Root1stActStructElms_rpvr(:,1,1,1)=0._r8
  Root1stLigStructElms_rpvr(ielmc,1,1,1)=1._r8
  RootMyco1stStrutElms_rpvr=Root1stLigStructElms_rpvr
 endif
 if(testcase==5)Root1stLenPP_rpvr=0._r8
 if(testcase==6)Num1stAxesPerStructRootXPOP_pft=0._r8
 if(testcase==7)RootMycoNonstElms_rpvr(ielmc,1,1,1)=0._r8
 RAutoRootO2Limter_rpvr=1._r8
 if(testcase==8)RAutoRootO2Limter_rpvr=0._r8
 if(testcase==9)then
  Root1stActStructElms_rpvr=Root1stActStructElms_rpvr*.001_r8
  RootMyco1stStrutElms_rpvr=Root1stActStructElms_rpvr
  RootMyco1stElm_raxs(:,1,1)=Root1stActStructElms_rpvr(:,1,1,1)
  RootMycoNonstElms_rpvr=RootMycoNonstElms_rpvr*.001_r8
  Num1stAxesPerStructRootXPOP_pft=Num1stAxesPerStructRootXPOP_pft*.001_r8
 endif
 before=RootMycoNonstElms_rpvr(:,1,1,1)+RootMyco1stStrutElms_rpvr(:,1,1,1)
 oldc=RootMyco1stStrutElms_rpvr(ielmc,1,1,1);litter=0._r8;resp=0._r8
 call StructuralRootThickening(yearIJ,1,1,1,1,.1_r8,sink,primarysink,1._r8,temp,litter,resp,testcase==7.or.testcase==8)
 after=RootMycoNonstElms_rpvr(:,1,1,1)+RootMyco1stStrutElms_rpvr(:,1,1,1)+litter
 after(ielmc)=after(ielmc)-resp
 if(any(abs(after-before)>1.e-12_r8))stop 10
 if(any(RootMycoNonstElms_rpvr(:,1,1,1)<0._r8))stop 11
 if(testcase==1)then
  if(RootMyco1stStrutElms_rpvr(ielmc,1,1,1)<=oldc)stop 12
  growth=RootMyco1stStrutElms_rpvr(ielmc,1,1,1)-oldc
 elseif(testcase==9)then
  if(abs((RootMyco1stStrutElms_rpvr(ielmc,1,1,1)-oldc)/.001_r8-growth)>1.e-16_r8)stop 14
 else
  if(RootMyco1stStrutElms_rpvr(ielmc,1,1,1)/=oldc)stop 13
 endif
enddo
! Reservation must include primary maintenance below the former N threshold.
! Medium and primary N are equal, so medium gets half the limited capacity.
RAutoRootO2Limter_rpvr=1._r8
Root1stLenPP_rpvr=.1_r8;Num1stAxesPerStructRootXPOP_pft=.01_r8
Root1stActStructElms_rpvr(:,1,1,1)=[1.e-6_r8,5.e-9_r8,5.e-10_r8]
RootMediumStructElms_rpvr(:,1,1,1)=Root1stActStructElms_rpvr(:,1,1,1)
iPlantRootProfile_pft=3;iPlant2ndGrothPattern_pft=1
NumStructuralRootAxes_pft=1;NRoot1stTipLay_raxes=2
MaxNumRootLays=2;Root1stAxesTipDepz2Surf_pft=1._r8
CumSoilThickness_vr=0._r8;flag2ndGrowth_pvr=.true.
fTGrowRootP_vr=1._r8;VMXC=1._r8
found=.false.;done=.false.;tips=0._r8
RootMycoNonstElms_rpvr(ielmc,1,1,1)=5.e-11_r8
call ReserveMediumRootMaintenance(1,1,1,temp,1._r8,tips,found,done,unlim,limited)
expected=2.5e-11_r8
if(abs(limited(1)-expected)>1.e-20_r8)stop 20
if(abs(RootMycoNonstElms_rpvr(ielmc,1,1,1)-expected)>1.e-20_r8)stop 21
end associate
print *,'PASS: low-N recovery, C/N/P/O2 limits, population scaling, inactive/empty geometry, CNP balance, maintenance reservation'
end program
'''
        source=source.replace('  public :: RootBGCModel','  public :: RootBGCModel, StructuralRootThickening, ReserveMediumRootMaintenance',1)
        flags=(BUILD/'f90src/Plant_bgc/CMakeFiles/plant_bgc.dir/flags.make').read_text()
        includes=shlex.split(re.search(r'^Fortran_INCLUDES = (.*)$',flags,re.M)[1])
        includes+=['-I'+str(BUILD/'f90src/Plant_bgc')]
        opts=shlex.split(re.search(r'^Fortran_FLAGS = (.*)$',flags,re.M)[1])
        link=shlex.split((BUILD/'drivers/ecosim/CMakeFiles/ecosim.f90.x.dir/link.txt').read_text())
        if '--' in link:link=link[link.index('--')+1:]
        cleaned=[];i=0
        while i<len(link):
            if link[i]=='-o':i+=2;continue
            if not link[i].endswith('.o'):cleaned.append(link[i])
            i+=1
        with tempfile.TemporaryDirectory(prefix='ecosim_recovery_') as work:
            w=Path(work);(w/'RootMod.F90').write_text(source);(w/'test.F90').write_text(driver)
            commands=[cleaned[:1]+includes+opts+['-J'+work,'-c',str(w/'RootMod.F90'),'-o',str(w/'RootMod.o')],
                      cleaned[:1]+['-I'+work]+includes+opts+[str(w/'test.F90'),str(w/'RootMod.o'),'-o',str(w/'test.x')]+cleaned[1:],
                      [str(w/'test.x')]]
            for command in commands:
                p=subprocess.run(command,cwd=BUILD/'drivers/ecosim',capture_output=True,text=True)
                self.assertEqual(p.returncode,0,p.stdout+p.stderr)
            self.assertIn('PASS:',p.stdout)


if __name__=='__main__':
    unittest.main()
