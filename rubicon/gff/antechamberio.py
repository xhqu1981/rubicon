from collections import namedtuple
from rubicon.gff.ffmol import FFmol

__author__ = 'navnidhirajput'

""" A wrapper for AntechamberRunner which generates force field files
    for a specified molecule using pdb file as input
"""

import subprocess
from pymatgen import write_mol
import shlex
from gff import Gff
from monty.io import ScratchDir
import tempfile
from rubicon.gff.topology import TopMol




class AntechamberRunner():
    """
    A wrapper for AntechamberRunner software

    """
    def __init__(self, mols):
        """
        Args:
            mols: List of molecules
        """

        self.mols = mols


    def _convert_to_pdb(self, molecule, filename=None):
        """
        generate pdb file for a given molecule
        """
        write_mol(molecule, filename)


    def _run_parmchk(self, filename='ANTECHAMBER_AC.AC'):
        """
        run parmchk using ANTECHAMBER_AC.AC file

        Args:
            filename: pdb file of the molecule
        """
        command_parmchk = (
        'parmchk -i ' + filename + ' -f ac -o mol.frcmod -a Y')
        return_cmd = subprocess.call(shlex.split(command_parmchk))
        return return_cmd


    def get_ff_top_mol(self,mol,filename=None):
        """
        run antechamber for a specified pdb file then run parmchk
        to generate missing force field parameters. Store and return
        the force field and topology information in ff_mol.

        Args:
            filename: pdb file of the molecule
            mol: molecule
        Returns:
            ff_mol : Object of FFmol which contains information of
                     force field and topology
        """
        scratch = tempfile.gettempdir()

        with ScratchDir(scratch,
                        copy_to_current_on_exit=False) as d:

            self._convert_to_pdb(mol, 'mol.pdb')
            command = (
            'antechamber -i ' + filename + ' -fi pdb -o ' + filename[:-4] +
            " -fo charmm")
            return_cmd = subprocess.call(shlex.split(command))
            self.molname = filename.split('.')[0]
            self._run_parmchk()
            top = TopMol.from_file('mol.rtf')
            gff = Gff.from_forcefield_para('ANTECHAMBER.FRCMOD')
            gff.read_atom_index(mol,'ANTECHAMBER_AC.AC')
            gff.read_charges()
            mol.add_site_property("atomname",(gff.atom_index.values()))
        ffmol= FFmol(gff,top)
        return ffmol














