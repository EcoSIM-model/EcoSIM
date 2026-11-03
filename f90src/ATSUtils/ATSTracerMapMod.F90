module ATSTracerMapMod
!
! Maps ATS transport component indices to EcoSIM TracerIDMod species
! indices. The canonical order defined here MUST match the
! "component names" array in the ATS transport XML.
!
! Also provides bidirectional conversion subroutines:
!   UnpackMoleFracToEcoSIM  — ATS mole fractions  -> EcoSIM concentrations
!   PackSourceFromEcoSIM    — EcoSIM source/sinks  -> ATS source terms
!
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use TracerIDMod
  implicit none
  private

  character(len=*), parameter :: mod_filename = __FILE__

  integer, parameter, public :: PHASE_GAS = 1, PHASE_SOLUTE = 2

  type, public :: tracer_map_entry
    integer  :: ats_index       ! 0-based index in a_MFrac
    integer  :: ecosim_index    ! TracerIDMod index (idg_* or ids_*)
    integer  :: phase           ! PHASE_GAS or PHASE_SOLUTE
    real(r8) :: mol_weight      ! g/mol
    character(len=16) :: name
  end type

  integer, public :: n_mapped_species = 0
  type(tracer_map_entry), allocatable, public :: species_map(:)

  public :: InitSpeciesMap
  public :: UnpackMoleFracToEcoSIM
  public :: PackSourceFromEcoSIM

contains

  subroutine InitSpeciesMap(num_components)
    implicit none
    integer, intent(in) :: num_components
    integer :: k

    if (allocated(species_map)) deallocate(species_map)
    allocate(species_map(num_components))
    k = 0

    ! Canonical order — must match XML "component names"
    ! Gases first (8 species)
    if (num_components > k) then; k=k+1
      species_map(k) = tracer_map_entry(k-1, idg_CO2, PHASE_GAS, 44.01_r8,  'CO2')
    endif
    if (num_components > k) then; k=k+1
      species_map(k) = tracer_map_entry(k-1, idg_CH4, PHASE_GAS, 16.04_r8,  'CH4')
    endif
    if (num_components > k) then; k=k+1
      species_map(k) = tracer_map_entry(k-1, idg_O2,  PHASE_GAS, 32.00_r8,  'O2')
    endif
    if (num_components > k) then; k=k+1
      species_map(k) = tracer_map_entry(k-1, idg_N2,  PHASE_GAS, 28.01_r8,  'N2')
    endif
    if (num_components > k) then; k=k+1
      species_map(k) = tracer_map_entry(k-1, idg_N2O, PHASE_GAS, 44.01_r8,  'N2O')
    endif
    if (num_components > k) then; k=k+1
      species_map(k) = tracer_map_entry(k-1, idg_H2,  PHASE_GAS,  2.016_r8, 'H2')
    endif
    if (num_components > k) then; k=k+1
      species_map(k) = tracer_map_entry(k-1, idg_AR,  PHASE_GAS, 39.95_r8,  'Ar')
    endif
    if (num_components > k) then; k=k+1
      species_map(k) = tracer_map_entry(k-1, idg_NH3, PHASE_GAS, 17.03_r8,  'NH3')
    endif

    ! Nutrients (if num_components > 8)
    if (num_components > k) then; k=k+1
      species_map(k) = tracer_map_entry(k-1, ids_NH4,  PHASE_SOLUTE, 18.04_r8, 'NH4')
    endif
    if (num_components > k) then; k=k+1
      species_map(k) = tracer_map_entry(k-1, ids_NO3,  PHASE_SOLUTE, 62.00_r8, 'NO3')
    endif
    if (num_components > k) then; k=k+1
      species_map(k) = tracer_map_entry(k-1, ids_NO2,  PHASE_SOLUTE, 46.01_r8, 'NO2')
    endif
    if (num_components > k) then; k=k+1
      species_map(k) = tracer_map_entry(k-1, ids_H1PO4, PHASE_SOLUTE, 96.99_r8, 'H1PO4')
    endif
    if (num_components > k) then; k=k+1
      species_map(k) = tracer_map_entry(k-1, ids_H2PO4, PHASE_SOLUTE, 97.99_r8, 'H2PO4')
    endif

    n_mapped_species = k

    write(*,*) '(ATSTracerMapMod) Initialized species map with ', n_mapped_species, ' species'
    do k = 1, n_mapped_species
      write(*,'(A,I3,A,A,A,I3,A,I2)') '  k=', k, ' ', trim(species_map(k)%name), &
        ' ecosim_id=', species_map(k)%ecosim_index, ' phase=', species_map(k)%phase
    enddo
  end subroutine InitSpeciesMap


  subroutine UnpackMoleFracToEcoSIM(NYS)
  !
  ! Convert ATS mole fractions [mol/mol] in a_MFrac to EcoSIM
  ! concentration arrays [g m⁻³] and mass arrays [g d⁻²].
  !
  ! a_MFrac(L,NY,k) * a_LDENS(L,NY) gives mol/m³ of liquid;
  ! multiply by mol_weight to get g/m³.
  !
    use SharedDataMod,        only : a_MFrac, a_LDENS, a_NU, a_NL
    use SoilBGCDataType,      only : trcg_gascl_vr, trcg_gasml_vr, &
                                     trc_solcl_vr, trcs_solml_vr
    use SoilWaterDataType,    only : VLWatMicP_vr
    use SoilPropertyDataType, only : VGeomLayer_vr
    implicit none
    integer, intent(in) :: NYS

    integer  :: NY, L, k, eid, NX
    real(r8) :: conc_g_m3

    NX = 1
    do NY = 1, NYS
      do L = a_NU(NY), a_NL(NY)
        do k = 1, n_mapped_species
          eid = species_map(k)%ecosim_index
          ! mol/mol * mol/m³ * g/mol = g/m³
          conc_g_m3 = a_MFrac(L,NY,species_map(k)%ats_index+1) &
                    * a_LDENS(L,NY) * species_map(k)%mol_weight
          
          if (species_map(k)%phase == PHASE_GAS) then
            trcg_gascl_vr(eid, L, NY, NX) = conc_g_m3
            trcg_gasml_vr(eid, L, NY, NX) = conc_g_m3 * VGeomLayer_vr(L,NY,NX)
            ! ATS passes dissolved concentrations; populate solml so GetChemicalConcsInSoil
            ! can compute trc_solcl_vr correctly (gas mass alone is zeroed in saturated soil)
            trcs_solml_vr(eid, L, NY, NX) = conc_g_m3 * VLWatMicP_vr(L, NY, NX)
          else
            trc_solcl_vr(eid, L, NY, NX)  = conc_g_m3
            ! mass = concentration * liquid water volume in layer
            trcs_solml_vr(eid, L, NY, NX) = conc_g_m3 * VLWatMicP_vr(L,NY,NX)
          endif
        enddo
      enddo
    enddo
  end subroutine UnpackMoleFracToEcoSIM


  subroutine PackSourceFromEcoSIM(NYS)
  !
  ! Convert EcoSIM microbial source/sink rates to ATS transport
  ! source terms in a_MFrac_ss_source [mol m⁻³ s⁻¹].
  !
  ! Gas uptake (trcs_RMicbUptake_vr): positive = consumption.
  !   ATS convention: positive = production. So negate.
  !   Units: [g d⁻² h⁻¹] → [mol m⁻³ s⁻¹]
  !   = (-rate / mol_weight) / Volume / 3600
  !
  ! Nutrient release (RNut_MicbRelease_vr): positive = release.
  !   ATS convention: positive = production. Same sign.
  !   Units: same conversion.
  !
    use SharedDataMod,    only : a_MFrac_ss_source, a_Volume, a_NU, a_NL
    use SoilBGCDataType,  only : trcs_RMicbUptake_vr, RNut_MicbRelease_vr
    implicit none
    integer, intent(in) :: NYS

    integer  :: NY, L, k, eid, NX
    real(r8) :: rate_g, rate_mol_m3_s

    NX = 1
    do NY = 1, NYS
      do L = a_NU(NY), a_NL(NY)
        do k = 1, n_mapped_species
          eid = species_map(k)%ecosim_index

          if (species_map(k)%phase == PHASE_GAS) then
            ! idg_AR (=7) is inert: trcs_RMicbUptake_vr upper bound is idg_NH3-1=7
            if (eid <= idg_NH3-1) then
              rate_g = trcs_RMicbUptake_vr(eid, L, NY, NX)
            else
              rate_g = 0._r8
            endif
            ! Negate: EcoSIM positive=consumption, ATS positive=production
            rate_mol_m3_s = (-rate_g / species_map(k)%mol_weight) &
                          / a_Volume(L,NY) / 3600.0_r8

          else ! PHASE_SOLUTE
            rate_g = RNut_MicbRelease_vr(eid, L, NY, NX)
            ! Same sign: EcoSIM positive=release = ATS positive=production
            rate_mol_m3_s = (rate_g / species_map(k)%mol_weight) &
                          / a_Volume(L,NY) / 3600.0_r8
          endif

          a_MFrac_ss_source(L, NY, species_map(k)%ats_index+1) = rate_mol_m3_s
        enddo
      enddo
    enddo
  end subroutine PackSourceFromEcoSIM

end module ATSTracerMapMod
