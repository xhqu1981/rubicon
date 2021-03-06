# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

from fireworks.core.firework import FireTaskBase, FWAction
from fireworks.utilities.fw_serializers import FWSerializable

from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.matproj.snl import StructureNL
from rubicon.utils.snl.egsnl_mongo import EGSNLMongoAdapter


class AddEGSNLTask(FireTaskBase, FWSerializable):
    """
    Add a new SNL into the SNL database, and build duplicate groups
    """

    _fw_name = "Add EG SNL Task"

    def run_task(self, fw_spec):

        sma = EGSNLMongoAdapter.auto_load()
        if isinstance(fw_spec['snl'], dict):
            snl = StructureNL.from_dict(fw_spec['snl'])
        else:
            snl = fw_spec['snl']
        egsnl, snlgroup_id = sma.add_snl(snl)

        mol = egsnl.structure
        bb = BabelMolAdaptor(mol)
        pbmol = bb.pybel_mol
        inchi_root = pbmol.write(str("inchi")).strip()

        return FWAction(update_spec={'egsnl': egsnl.as_dict(),
                                     'snlgroup_id': snlgroup_id,
                                     'inchi_root': inchi_root})
