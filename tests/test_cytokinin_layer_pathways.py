"""Compile and exercise the production layer-pathway cytokinin routines.

Expose private procedures only in a temporary module; link the actual EcoSIM
libraries. Tests check inventories including production, decay, and shoot export.
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


class CytokininPathwaysTest(unittest.TestCase):
    def test_conservation_topology_and_geometry(self):
        source = (ROOT / 'f90src/Plant_bgc/RootMod.F90').read_text()
        owners = {'plt_site': 'PlantSite', 'plt_morph': 'PlantMorphology',
                  'plt_biom': 'PlantBiomass', 'plt_pheno': 'PlantPhenology',
                  'plt_ew': 'PlantEnergyWater', 'plt_rbgc': 'PlantRootBGC'}
        fields = set()
        for name in ['CytoKininDynamics', 'RootCytokininInventories', 'MediumTransportGeometry', 'CalcPotentialSapArea']:
            start = source.lower().index(('  function ' if name.startswith('Calc') else '  subroutine ') + name.lower()+'(')
            end = source.lower().index('  end '+('function ' if name.startswith('Calc') else 'subroutine ')+name.lower(), start)
            fields.update(re.findall(r'(plt_\w+)%([\w]+)', source[start:end]))
        allocation = []
        for owner, field in sorted(fields):
            api = (ROOT/'f90src/APIData'/(owners[owner]+'APIData.F90')).read_text()
            decl = re.search(r'::\s*'+field+r'\s*(\([^\n]*?\))?', api, re.I)
            self.assertIsNotNone(decl, field)
            rank = (decl[1] or '').count(':')
            if not rank:
                continue
            shape = {1:'1', 2:'3,1', 3:'3,2,1', 4:'2,3,2,1'}[rank]
            if owner == 'plt_site': shape = '3'
            if field in ['Root2ndXSecArea_pft', 'RootPorosity_pft']: shape = '2,1'
            if field in ['RPlantRootH2OUptk_pvr', 'Root2ndXNumL_rpvr']: shape = '2,3,1'
            if field in ['RootMediumStructElms_rpvr', 'RootMediumYoungElms_rpvr']: shape = '3,3,2,1'
            allocation += [f'allocate({owner}%{field}({shape}))',f'{owner}%{field}=0']
        driver = '''program pathway_test
use data_kind_mod,only:r8=>DAT_KIND_R8
use RootMod,only:CytoKininDynamics,RootCytokininInventories
use PlantAPICommonData
use PlantBGCPars
use EcoSIMCtrlMod,only:lcoarseroot
use ElmIDMod
'''+ '\n'.join('use '+v+'APIData,only:'+k for k,v in owners.items())+'''
implicit none
real(r8)::fm(2,3,2),mm(3,2),pm(3,2),af(2,3,2),am(3,2),ap(3,2),sink(2,3,2)
real(r8)::before,after,production,err,upper,expected
integer::scenario
allocate(JZ1);JZ1=3
pltpar%jroots=2;pltpar%MaxNumRootAxes=2
'''+ '\n'.join(allocation)+'''
plt_site%NU=1;plt_site%NK=3
plt_site%DLYR3=.2_r8;plt_site%PlantPopuLive_pft=.1_r8
plt_morph%Myco_pft=2;plt_morph%NumStructuralRootAxes_pft=2
plt_morph%MaxSoilLays4Root_pft=2
plt_pheno%iPlantRootProfile_pft=3;plt_pheno%iPlant2ndGrothPattern_pft=1
plt_morph%Radius95pctMature_pft=.05_r8
plt_morph%xylemPhi_min_pft=.4_r8;plt_morph%xylemPhi_max_pft=.7_r8
plt_morph%StalkAveRadius_pft=.02_r8
plt_morph%Root2ndXSecArea_pft=1.e-7_r8
plt_morph%RootPorosity_pft=.3_r8
plt_morph%CRootActVolPerMassC_pft=acos(-1._r8)*(.0005_r8)**2*.4_r8
lcoarseroot=.true.;sink=1._r8
! Every scenario checks the whole-root signal inventory, including both categories.
do scenario=1,10
 plt_morph%RootMediumYoungXNum_rpvr=0._r8
 plt_morph%RootMediumYoungLength_rpvr=0._r8
 plt_morph%RootMediumYoungElms_rpvr=0._r8
 plt_morph%Root1stRadius_rpvr=.001_r8
 plt_morph%RootMediumRadius_rpvr=.0005_r8
 plt_morph%RootMediumLength_rpvr=.4_r8
 plt_morph%RootMediumXNum_rpvr=2._r8
 plt_biom%RootMediumStructElms_rpvr=1._r8
 plt_morph%Num1stAxesPerStructRootXPOP_pft=.1_r8
 plt_morph%Root2ndXNum_rpvr=1._r8
 plt_morph%Root2ndXNumL_rpvr=2._r8
 plt_morph%Root2ndLen_rpvr=1._r8
 plt_morph%Root2ndVH2O_rpvr=7.e-8_r8
 plt_morph%CRootLumenArea_rpvr=1.e-6_r8
 plt_rbgc%Cytokinin2ndConc_rpvr=.1_r8
 plt_rbgc%CytokininMRConc_rpvr=.2_r8
 plt_rbgc%Cytokinin1stConc_rpvr=.3_r8
 plt_rbgc%Root2ndProdCytok_rpvr=1.e-9_r8
 plt_rbgc%RootMRProdCytok_rpvr=2.e-9_r8
 plt_ew%RPlantRootH2OUptk_pvr=-1.e-6_r8
 if(scenario==1)plt_ew%RPlantRootH2OUptk_pvr=-1.e3_r8 ! huge Courant number
 if(scenario==2.or.scenario==5.or.scenario==6)plt_ew%RPlantRootH2OUptk_pvr=0._r8
 if(scenario==7)plt_ew%RPlantRootH2OUptk_pvr=1.e-6_r8 ! reverse flow: no acropetal export
 if(scenario==4)then
  plt_biom%RootMediumStructElms_rpvr=0._r8
  plt_morph%RootMediumLength_rpvr=0._r8
  plt_morph%RootMediumXNum_rpvr=0._r8
  plt_rbgc%RootMRProdCytok_rpvr=0._r8
 endif
 if(scenario==10)then
  ! The young and established lumen volumes must share one conservative signal inventory.
  plt_morph%RootMediumYoungXNum_rpvr=.5_r8
  plt_morph%RootMediumYoungLength_rpvr=.005_r8
  plt_morph%RootMediumYoungElms_rpvr=.0125_r8
 endif
 call RootCytokininInventories(1,fm,mm,pm)
 before=sum(fm)+sum(mm)+sum(pm)
 production=sum(plt_rbgc%Root2ndProdCytok_rpvr(:,1:2,:,1))+sum(plt_rbgc%RootMRProdCytok_rpvr(1:2,:,1))
 if(scenario==5)then
  ! Geometry grows after the snapshot: inventory must be unchanged except explicit fluxes.
  plt_morph%Root1stRadius_rpvr=.002_r8
  plt_morph%RootMediumRadius_rpvr=.001_r8
  plt_morph%RootMediumLength_rpvr=.8_r8
  plt_biom%RootMediumStructElms_rpvr=8._r8
  plt_morph%Root2ndLen_rpvr=2._r8
  plt_morph%Root2ndVH2O_rpvr=14.e-8_r8
 endif
 if(scenario==9)plt_morph%Root1stRadius_rpvr=0._r8
 if(scenario==6)then
  ! A disappearing compartment has an explicit loss, not silent signal creation.
  plt_morph%Root2ndLen_rpvr=0._r8;plt_morph%Root2ndVH2O_rpvr=0._r8
  plt_biom%RootMediumStructElms_rpvr=0._r8
  plt_morph%RootMediumLength_rpvr=0._r8;plt_morph%RootMediumXNum_rpvr=0._r8
 endif
 if(scenario==3.or.scenario==8)then
  ! No local upper-layer uptake: deeper uptake must not flush upper primary roots.
  plt_ew%RPlantRootH2OUptk_pvr(:,1,1)=0._r8
  if(scenario==8)plt_ew%RPlantRootH2OUptk_pvr(:,2,1)=-1.e3_r8
 endif
 call CytoKininDynamics(1,1,1,sink,fm,mm,pm)
 call RootCytokininInventories(1,af,am,ap)
 after=sum(af)+sum(am)+sum(ap)
 err=after+plt_rbgc%CytokininShootExport_pft(1)+plt_rbgc%CytokininRootLoss_pft(1)-before-production
 if(abs(err)>1.e-12_r8*(before+production))then
  print *,scenario,err,before,after
  stop 10
 endif
 if(minval(af)<0.or.minval(am)<0.or.minval(ap)<0)stop 11
 if(scenario==2.or.scenario==5.or.scenario==6.or.scenario==7)then
  if(plt_rbgc%CytokininShootExport_pft(1)/=0._r8)stop 12
 endif
 if(scenario==5)then
  expected=(fm(1,1,1)+1.e-9_r8)*exp(-kDCytof(1))
  if(abs(af(1,1,1)-expected)>1.e-20_r8)stop 16
 endif
 if(scenario==3)upper=plt_rbgc%Cytokinin1stConc_rpvr(1,1,1)
 if(scenario==8)then
  if(plt_rbgc%Cytokinin1stConc_rpvr(1,1,1)/=upper)stop 13
 endif
 if(scenario==4)then
  if(any(plt_rbgc%CytokininMRConc_rpvr(1:2,:,1)/=0._r8))stop 14
 endif
enddo
! No structural axes: diagnostics reset, with no stale export from previous PFT/call.
plt_morph%NumStructuralRootAxes_pft=0
call CytoKininDynamics(1,1,1,sink,fm,mm,pm)
if(plt_rbgc%CytokininShootExport_pft(1)/=0._r8)stop 15
print *,'PASS: conservation, zero/high/reverse flow, layer independence, direct routing, geometry, disappearance, reset'
end program
'''
        source=source.replace('  public :: RootBGCModel','  public :: RootBGCModel, CytoKininDynamics, RootCytokininInventories',1)
        flags=(BUILD/'f90src/Plant_bgc/CMakeFiles/plant_bgc.dir/flags.make').read_text()
        includes=shlex.split(re.search(r'^Fortran_INCLUDES = (.*)$',flags,re.M)[1])+['-I'+str(BUILD/'f90src/Plant_bgc')]
        opts=shlex.split(re.search(r'^Fortran_FLAGS = (.*)$',flags,re.M)[1])
        link=shlex.split((BUILD/'drivers/ecosim/CMakeFiles/ecosim.f90.x.dir/link.txt').read_text())
        if '--' in link: link=link[link.index('--')+1:]
        cleaned=[];i=0
        while i<len(link):
            if link[i]=='-o':i+=2;continue
            if not link[i].endswith('.o'):cleaned.append(link[i])
            i+=1
        with tempfile.TemporaryDirectory(prefix='ecosim_cytokinin_') as work:
            w=Path(work);(w/'RootMod.F90').write_text(source);(w/'test.F90').write_text(driver)
            commands=[cleaned[:1]+includes+opts+['-J'+work,'-c',str(w/'RootMod.F90'),'-o',str(w/'RootMod.o')],
                      cleaned[:1]+['-I'+work]+includes+opts+[str(w/'test.F90'),str(w/'RootMod.o'),'-o',str(w/'test.x')]+cleaned[1:],
                      [str(w/'test.x')]]
            for command in commands:
                p=subprocess.run(command,cwd=BUILD/'drivers/ecosim',capture_output=True,text=True)
                self.assertEqual(p.returncode,0,p.stdout+p.stderr)
            self.assertIn('PASS:',p.stdout)


if __name__=='__main__': unittest.main()
