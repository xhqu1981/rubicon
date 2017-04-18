# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import argparse

from pymatgen.core.structure import Molecule
from pymatgen.io.qchem import QcInput, QcNucVeloc
from pymatgen.io.xyz import MXYZ

__author__ = 'xiaohuiqu'


def main():
    parser = argparse.ArgumentParser(
        description="Replace the atom coordinates of the first job in QChem input file with the coordinates from"
                    "an XYZ file")
    parser.add_argument("-i", "--input", dest="input", type=str,
                        required=True,
                        help="the QChem input filename")
    parser.add_argument("-c", "--coords", dest="coords", type=str,
                        required=True,
                        help="The XYZ file contains the new coords")
    parser.add_argument("-v", "--velocity", dest="velocity", type=str,
                        default=None,
                        help="The AIMD velocity file")
    parser.add_argument("-o", "--output", dest="output", type=str,
                        required=True,
                        help="the QChem input filename with the coordinates from the XYZ file")
    options = parser.parse_args()
    qcinp = QcInput.from_file(options.input)
    charge, spin = qcinp.jobs[0].charge, qcinp.jobs[0].spin_multiplicity
    if options.velocity is None:
        new_mol = Molecule.from_file(options.coords)
    else:
        mxyz = MXYZ.from_file(options.coords)
        new_mol = mxyz.molecules[-1]
        qcinp.jobs[0].params["rem"].pop("aimd_init_veloc", None)
        qcnv = QcNucVeloc(options.velocity)
        qcinp.jobs[0].set_velocities(qcnv.velocities[-1])
    if charge is not None:
        new_mol.set_charge_and_spin(charge, spin)
    qcinp.jobs[0].mol = new_mol
    qcinp.write_file(options.output)
    print(
        "created new QChem input file {new_file} using {old_inp} as an template and filled with coordinates " \
        "from {coord_file}".format(old_inp=options.input,
                                   coord_file=options.coords,
                                   new_file=options.output))


if __name__ == "__main__":
    main()
