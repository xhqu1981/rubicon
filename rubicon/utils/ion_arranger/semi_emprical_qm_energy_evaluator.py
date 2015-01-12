import copy
import os
import shutil

from custodian import Custodian
from monty.tempfile import ScratchDir
from pymatgen.core.structure import Molecule
from pymatgen.core.units import Energy

from rubicon.io.mopacio.custodian.mopac_error_handlers import MopacErrorHandler
from rubicon.io.mopacio.custodian.mopacjob import MopacJob
from rubicon.io.mopacio.mopacio import MopOutput, MopTask
from rubicon.utils.ion_arranger.energy_evaluator import EnergyEvaluator
from rubicon.utils.ion_arranger.hard_sphere_energy_evaluators import HardSphereEnergyEvaluator, AtomicRadiusUtils, \
    ContactDetector, LargestContactGapEnergyEvaluator


__author__ = 'xiaohuiqu'


class SemiEmpricalQuatumMechanicalEnergyEvaluator(EnergyEvaluator):

    def __init__(self, ob_mol, ob_cation, ob_anion, total_charge, num_cation, num_anion,
                 lower_covalent_radius_scale=2.0, lower_metal_radius_scale=0.8,
                 upper_covalent_radius_scale=3.0, upper_metal_radius_scale=1.5,):
        from rubicon.utils.ion_arranger.ion_arranger import IonPlacer
        mol_coords = IonPlacer.normalize_molecule(ob_mol)
        super(SemiEmpricalQuatumMechanicalEnergyEvaluator, self).__init__(mol_coords)
        self.total_charge = total_charge
        self.lower_sphere = self._construct_hardsphere_energy_evaluator(
            lower_covalent_radius_scale, lower_metal_radius_scale,
            mol_coords, ob_mol, ob_cation, ob_anion)
        self.contact_detector = self._construct_contact_detector(
            upper_covalent_radius_scale, upper_metal_radius_scale,
            mol_coords, ob_mol, ob_cation, ob_anion)
        self.gravitation = self._construct_largest_cap_energy_evaluator(
            upper_covalent_radius_scale, upper_metal_radius_scale,
            num_cation, num_anion, mol_coords, ob_mol, ob_cation, ob_anion)
        self.mol_species = IonPlacer.get_mol_species(ob_mol)
        self.cation_species = IonPlacer.get_mol_species(ob_cation)
        self.anion_species = IonPlacer.get_mol_species(ob_anion)
        self.run_number = 1
        self.best_energy = 0.0

    def calc_energy(self, cation_coords, anion_coords):
        energy = self.lower_sphere.calc_energy(cation_coords, anion_coords)
        if energy > HardSphereEnergyEvaluator.overlap_energy * 0.9:
            return energy
        if not self.contact_detector.is_contact(cation_coords, anion_coords):
            return self.gravitation.calc_energy(cation_coords, anion_coords)
        mol = self._get_super_molecule(cation_coords, anion_coords)
        energy = self.run_mopac(mol)
        energy = round(energy, 3)
        # coarse grained energy,
        # make potential energy surface simpler
        return energy

    def run_mopac(self, mol):
        cur_dir = os.getcwd()
        all_errors = set()
        energy = 0.0
        with ScratchDir(rootpath=cur_dir, copy_from_current_on_enter=False, copy_to_current_on_exit=True):
            order_text = ["st", "nd", "th"]
            title = "Salt Alignment {}{} Calculation".format(
                self.run_number, order_text[self.run_number-1 if self.run_number < 3 else 2])
            self.run_number += 1
            mop = MopTask(mol, self.total_charge, "opt", title, "PM7", {"CYCLES": 1000})
            mop.write_file("mol.mop")
            job = MopacJob()
            handler = MopacErrorHandler()
            c = Custodian(handlers=[handler], jobs=[job], max_errors=50)
            custodian_out = c.run()
            mopout = MopOutput("mol.out")
            energy = Energy(mopout.data["energies"][-1][-1], "eV").to("Ha")
            if energy < self.best_energy:
                self.best_energy = energy
                shutil.copy("mol.out", os.path.join(cur_dir, "best_mol.out"))
            for run in custodian_out:
                for correction in run['corrections']:
                    all_errors.update(correction['errors'])
            if len(all_errors) > 0:
                out_error_logs_file = os.path.join(cur_dir, "out_error_logs.txt")
                bak_cmd = "cat {} >> {}".format("mol.out", out_error_logs_file)
                os.system(bak_cmd)
        return energy

    def _get_super_molecule(self, cation_coords, anion_coords):
        super_mol_species = []
        super_mol_coords_au = []
        super_mol_species.extend(copy.deepcopy(self.mol_species))
        super_mol_coords_au.extend(copy.deepcopy(self.mol_coords))
        for cc in cation_coords:
            super_mol_species.extend(copy.deepcopy(self.cation_species))
            super_mol_coords_au.extend(copy.deepcopy(cc))
        for ac in anion_coords:
            super_mol_species.extend(copy.deepcopy(self.anion_species))
            super_mol_coords_au.extend(copy.deepcopy(ac))
        super_mol_coords_ang = [[x/AtomicRadiusUtils.angstrom2au for x in coord] for coord in super_mol_coords_au]
        return Molecule(super_mol_species, super_mol_coords_ang)

    @staticmethod
    def _construct_hardsphere_energy_evaluator(covalent_radius_scale, metal_radius_scale,
                                               mol_coords, ob_mol, ob_cation, ob_anion):
        rad_util = AtomicRadiusUtils(covalent_radius_scale, metal_radius_scale)
        mol_radius = rad_util.get_radius(ob_mol)
        cation_radius = rad_util.get_radius(ob_cation)
        anion_radius = rad_util.get_radius(ob_anion)
        sphere = HardSphereEnergyEvaluator(
            mol_coords, mol_radius, cation_radius, anion_radius)
        return sphere

    @staticmethod
    def _construct_contact_detector(covalent_radius_scale, metal_radius_scale,
                                    mol_coords, ob_mol, ob_cation, ob_anion):
        rad_util = AtomicRadiusUtils(covalent_radius_scale, metal_radius_scale)
        mol_radius = rad_util.get_radius(ob_mol)
        cation_radius = rad_util.get_radius(ob_cation)
        anion_radius = rad_util.get_radius(ob_anion)
        detector = ContactDetector(mol_coords, mol_radius, cation_radius, anion_radius)
        return detector

    @staticmethod
    def _construct_largest_cap_energy_evaluator(covalent_radius_scale, metal_radius_scale,
                                                num_cation, num_anion,
                                                mol_coords, ob_mol, ob_cation, ob_anion):
        rad_util = AtomicRadiusUtils(covalent_radius_scale, metal_radius_scale)
        mol_radius = rad_util.get_radius(ob_mol)
        cation_radius = rad_util.get_radius(ob_cation)
        anion_radius = rad_util.get_radius(ob_anion)
        from rubicon.utils.ion_arranger.ion_arranger import IonPlacer
        bounder = IonPlacer.get_bounder(mol_coords, ob_cation, ob_anion, num_cation, num_anion)
        max_cap = max(bounder.upper_bound) * 2.0 / AtomicRadiusUtils.angstrom2au
        evaluator = LargestContactGapEnergyEvaluator(
            mol_coords, mol_radius, cation_radius, anion_radius, max_cap, threshold=0.01)
        return evaluator