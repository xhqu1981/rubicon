import shlex
import subprocess
from rubicon.gff.gaussian_input_old import GaussianInput

__author__ = 'navnidhirajput'

from fireworks import FireTaskBase, explicit_serialize, FWAction
from pymatgen import write_mol, Molecule
import glob
from pymatgen.io import gaussian


@explicit_serialize
class WritegaussianGeoTask(FireTaskBase):
    """
    Writes Gaussian Input files for geometry optimization.

    """

    _fw_name = "Gaussian Geometry Writer"

    def run_task(self, mol, charge=None, spin_multiplicity=None):
        moleculelist = glob.glob(
            "/Users/navnidhirajput/Dropbox/solvent_molecules/*")
        for filename in moleculelist:
            mol = Molecule.from_file(filename)
            gaus_lines = gaussian.GaussianInput(mol, charge=charge,
                                                spin_multiplicity = spin_multiplicity,
                                                title='created by gaussian_geo_task from' + ' ' + filename[
                                                                                                    48:],
                                                functional="b3lyp",
                                                basis_set="aug-cc-pvdz",
                                                route_parameters={
                                                    'opt': "(calcfc,tight)",
                                                    'int': "ultrafine",
                                                    "\n# SCF": "tight  nosymm test"},
                                                input_parameters=None,
                                                link0_parameters={
                                                    "%mem": "256MW",
                                                    "%NProcShared": 4,
                                                    "%LindaWorker": "localhost",
                                                    "%chk": filename[
                                                            48:-4] + ".chk"},
                                                dieze_tag="#",
                                                gen_basis=None)

        with open('mol_geo.gau', 'w') as f:
            gaus_lines.write_file('mol_geo.gau', cart_coords=True)

        subprocess.check_call(shlex.split("g09 "), stdin=f)

        prev_gaussian_geo = shlex.os.path.join(shlex.os.getcwd(), filename[
                                                            48:-4]+'.out')

        update_spec = {'prev_gaussian_geo': prev_gaussian_geo}

        return FWAction(update_spec=update_spec)




class WritegaussianFreqTask(FireTaskBase):
    """
    Writes Gaussian Input files for frequency and charge calculation.

    """

    _fw_name = "Gaussian frequency and charge Writer"

    def run_task(self, mol, charge=None, spin_multiplicity=None):
        moleculelist = glob.glob(
            "/Users/navnidhirajput/Dropbox/solvent_molecules/*")
        for filename in moleculelist:
            mol = Molecule.from_file(filename)

            gaus_freq_charge = gaussian.GaussianInput(mol, charge = charge,
                                                      spin_multiplicity=spin_multiplicity,
                                                      title='created by gaussian_frq_task from' + ' ' + filename[
                                                                                                          48:],
                                                      functional="b3lyp",
                                                      basis_set="aug-cc-pvdz  freq",
                                                      route_parameters={
                                                          '\n# geom': "allcheck",
                                                          'guess': "read",
                                                          "SCF": "tight  nosymm test",
                                                      "pop":"MK iop(6/33=2,6/41=10,6/42=10,7/33=1)"},
                                                      input_parameters=None,
                                                      link0_parameters={
                                                          "%mem": "256MW",
                                                          "%NProcShared": 4,
                                                          "%LindaWorker": "localhost",
                                                          "%chk": filename[
                                                                  48:-4] + ".chk"},
                                                      dieze_tag="#",
                                                      gen_basis=None)
        with open('mol_freq.gau', 'w') as f:
            gaus_freq_charge.write_file('mol_freq.gau', cart_coords=True)

        #with open("mol_freq.gau") as f:
        subprocess.check_call(shlex.split("g09 "), stdin=f)

        return FWAction()


# moleculelist = glob.glob("/Users/navnidhirajput/Dropbox/solvent_molecules/*")
# for filename in moleculelist:
#     mol = Molecule.from_file(filename)
#     task_geo = WritegaussianGeoTask()
#     task_geo.run_task(mol, -1,1)
#     task_freq = WritegaussianFreqTask()
#     task_freq.run_task(mol, -1, 1)
