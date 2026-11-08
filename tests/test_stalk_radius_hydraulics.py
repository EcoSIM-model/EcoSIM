"""Compile the production hydraulic procedure and exercise its transport scaling.

The procedure is private, so compile its unchanged source body in a test host
against the built EcoSIM API modules. No alternative hydraulic implementation
is substituted. Set ECOSIM_TEST_BUILD to select a configured build.
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


class StalkRadiusHydraulicsTest(unittest.TestCase):
    def test_geometry_population_paths_and_legacy_switch(self):
        source = (ROOT / 'f90src/Plant_bgc/UptakesMod.F90').read_text()
        start = source.index('  subroutine CalcPlantHydroResistance(')
        end = source.index('  end subroutine CalcPlantHydroResistance', start)
        procedure = source[start:end] + '  end subroutine CalcPlantHydroResistance\n'
        associations = re.findall(r'(\w+)\s*=>\s*(plt_\w+)%([\w]+)', procedure)
        # Allocate every pointer selector in the production associate block,
        # including selectors not used by the tested pathways.
        allocation = []
        for alias, owner, field in associations:
            api = next((ROOT / 'f90src/APIData').glob({'plt_site':'PlantSiteAPIData.F90',
                'plt_morph':'PlantMorphologyAPIData.F90','plt_pheno':'PlantPhenologyAPIData.F90',
                'plt_soilchem':'PlantSoilChemistryAPIData.F90','plt_ew':'PlantEnergyWaterAPIData.F90',
                'plt_biom':'PlantBiomassAPIData.F90'}[owner])).read_text()
            declaration = re.search(r'::\s*'+field+r'\s*(\([^\n]*?\))?', api, re.I)
            self.assertIsNotNone(declaration, field)
            rank = (declaration[1] or '').count(':')
            if rank == 0:
                continue  # The site scalars have ordinary storage.
            shape = '1' if rank == 1 and owner in ['plt_morph','plt_pheno','plt_biom','plt_ew'] else '2'
            if rank == 2:
                shape = '3,2' if field == 'VLWatMicPM_vr' else '2,1'
            if rank == 3:
                shape = '2,2,1'
            allocation.append(f'  allocate({owner}%{field}({shape}))')
        aliases = ', &\n    '.join(f'{a} => {o}%{f}' for a,o,f in associations)
        harness = '''program test_stalk_radius
  use data_kind_mod, only: r8 => DAT_KIND_R8, yearIJ_type
  use EcoSIMCtrlMod, only: lcoarseroot
  use RootMod, only: MediumTransportGeometry
  use PlantAPICommonData
  use PlantSiteAPIData, only: plt_site
  use PlantMorphologyAPIData, only: plt_morph
  use PlantPhenologyAPIData, only: plt_pheno
  use PlantSoilChemistryAPIData, only: plt_soilchem
  use PlantEnergyWaterAPIData, only: plt_ew
  use PlantBiomassAPIData, only: plt_biom
  use PlantBGCPars
  use UptakePars
  use EcoSIMSolverPar, only: NPH
  use minimathmod
  use EcosimConst
  use DebugToolMod
  use PlantMathFuncMod
  implicit none
  character(len=*), parameter :: mod_filename='hydraulic test'
  type(yearIJ_type) :: yearIJ
  real(r8) :: radial(2,2),radius(2,2),efflen(2,2),soil(2,2),total(2,2),cd,grav
  real(r8) :: base,wide,tall,twopop,paths,expected,area,depth,phi
  logical :: hasroot(2,2)
  allocate(JZ1);JZ1=2
  pltpar%jroots=2
''' + '\n'.join(allocation) + '\n  associate( &\n    ' + aliases + ''')
    NU=1;ZERO=1.e-30_r8;ZEROS2=1.e-30_r8;ZERO4Groth_pft=1.e-30_r8
    MaxSoilLays4Root_pft=1;Myco_pft=1
    DLYR3=1._r8;CumSoilThickMidL_vr=0.5_r8
    VLSoilPoreMicP_vr=1._r8;VLWatMicPM_vr=1._r8;VLMicP_vr=1._r8
    HYCDMicP4RootUptake_vr=1.e30_r8;THETW_vr=1._r8
    CanopyHeightLive_pft=10._r8;StalkAveRadius_pft=0.01_r8
    PlantPopuLive_pft=1._r8;PSICanopy_pft=0._r8
    iPlantRootProfile_pft=3;iPlant2ndGrothPattern_pft=1
    xylemPhi_min_pft=0.4_r8;xylemPhi_max_pft=0.7_r8;Radius95pctMature_pft=0.05_r8
    RootSingleVesselArea_pft=PiCON*(15.e-6_r8)**2
    RootSingleVesselRstaxial_pft=1._r8
    Root1stXNumL_pvr=10._r8;Root2ndXNumL_rpvr=1._r8
    Root2ndRadius_rpvr=0.00025_r8;Root2ndMaxRadius_pft=0.00025_r8
    RootAbsorbLenPerPlant_pvr=1._r8;Root2ndEffLen4uptk_rpvr=1._r8
    CRootLumenArea_pvr=1.e20_r8;MRootLumenArea_pvr=1.e20_r8
    RootMediumLength_pvr=0._r8;RootFineFrac2Med_rpvr=0._r8
    RootRadialResist_pft=1.e-30_r8;Root2ndAxialResist_pft=1.e-30_r8
    StalkAxialResist_pft=4.e9_r8
    EMODW=50._r8;NPH=1
    radial=1._r8;radius=0.00025_r8;efflen=1._r8
    lcoarseroot=.true.
    call evaluate(base)
    depth=dmax*tanh(0.01_r8/dmax)
    area=PiCON*depth*(0.02_r8-depth)
    phi=0.4_r8+0.3_r8*(1._r8-exp(-0.6_r8))
    expected=area*phi/RootSingleVesselArea_pft(1)/8._r8
    call close(base,expected,1)
    StalkAveRadius_pft=0.02_r8;call evaluate(wide)
    if(wide<=base)stop 2
    StalkAveRadius_pft=0.01_r8;CanopyHeightLive_pft=20._r8;call evaluate(tall)
    call close(tall,base/2._r8,3)
    CanopyHeightLive_pft=10._r8;PlantPopuLive_pft=2._r8;call evaluate(twopop)
    call close(twopop,2._r8*base,4)
    PlantPopuLive_pft=1._r8;MaxSoilLays4Root_pft=2;Myco_pft=2
    call evaluate(paths);call close(paths,base,5)
    ! Redistribution among unequal parallel paths must also conserve capacity.
    RootAbsorbLenPerPlant_pvr(1,1,1)=7._r8
    call evaluate(paths);call close(paths,base,6)
    RootAbsorbLenPerPlant_pvr=1._r8;MaxSoilLays4Root_pft=1;Myco_pft=1
    Root1stXNumL_pvr=5000._r8;call evaluate(paths);call close(paths,base,7)
    StalkAveRadius_pft=0._r8;call evaluate(paths)
    if(any(hasroot).or.paths>1.e-20_r8)stop 8
    StalkAveRadius_pft=0.01_r8;PlantPopuLive_pft=0._r8;call evaluate(paths)
    if(any(hasroot).or.paths>1.e-20_r8)stop 9
    PlantPopuLive_pft=1._r8;lcoarseroot=.false.;Root1stXNumL_pvr=10._r8
    call evaluate(paths)
    expected=1.e4_r8*10._r8/(4.e9_r8*8._r8)
    call close(paths,expected,10)
    StalkAveRadius_pft=0.02_r8;call evaluate(paths);call close(paths,expected,11)
    lcoarseroot=.true.;iPlant2ndGrothPattern_pft=0
    call evaluate(paths);call close(paths,expected,12)
  end associate
  print *, 'PASS: stalk geometry, height, population, shared paths, empty stem, and legacy switch'
contains
  subroutine evaluate(value)
    real(r8), intent(out) :: value
    call CalcPlantHydroResistance(yearIJ,1,radial,radius,efflen,soil,total,cd,grav,hasroot)
    value=cd
  end subroutine evaluate
  subroutine close(value,reference,code)
    real(r8),intent(in)::value,reference
    integer,intent(in)::code
    if(abs(value-reference)>1.e-10_r8*abs(reference))then
      print *,code,value,reference
      stop 99
    endif
  end subroutine close
''' + procedure + '\nend program test_stalk_radius\n'
        command = shlex.split((BUILD / 'drivers/ecosim/CMakeFiles/ecosim.f90.x.dir/link.txt').read_text())
        if '--' in command:
            command = command[command.index('--')+1:]
        cleaned=[];i=0
        while i<len(command):
            if command[i]=='-o':i+=2;continue
            if not command[i].endswith('.o'):cleaned.append(command[i])
            i+=1
        includes=sorted({'-I'+str(p.parent.resolve()) for p in BUILD.rglob('*.mod')})
        with tempfile.TemporaryDirectory(prefix='ecosim_stalk_test_') as work:
            program=Path(work)/'test.F90';exe=Path(work)/'test.x'
            program.write_text(harness)
            p=subprocess.run(cleaned[:1]+includes+['-ffree-line-length-none','-fcheck=all',str(program),'-o',str(exe)]+cleaned[1:],
                             cwd=BUILD/'drivers/ecosim',capture_output=True,text=True)
            self.assertEqual(p.returncode,0,p.stdout+p.stderr)
            p=subprocess.run([str(exe)],cwd=work,capture_output=True,text=True)
            self.assertEqual(p.returncode,0,p.stdout+p.stderr)
            self.assertIn('PASS:',p.stdout)


if __name__=='__main__':
    unittest.main()
