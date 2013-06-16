#!/usr/bin/env python

"""
TODO: Modify module doc.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "6/15/13"


import os
import logging
import datetime

from pymongo import MongoClient

from pymatgen.apps.borg.hive import AbstractDrone
from pymatgen.io.nwchemio import NwOutput
from pymatgen.util.io_utils import clean_json

logger = logging.getLogger(__name__)


class DeltaSCFNwChemToDbTaskDrone(AbstractDrone):
    """
    Assimilates a delta scf nwchem run and inserts it into the database.
    """

    #Version of this db creator document.
    __version__ = "0.1.0"

    def __init__(self, host="127.0.0.1", port=27017, database="mg_core_dev",
                 user=None, password=None,  collection="mol_tasks",
                 simulate_mode=False, update_duplicates=True):
        """
        Args:
            host:
                Hostname of database machine. Defaults to 127.0.0.1 or
                localhost.
            port:
                Port for db access. Defaults to mongo's default of 27017.
            database:
                Actual database to access. Defaults to "mg_core_dev".
            user:
                User for db access. Requires write access. Defaults to None,
                which means no authentication.
            password:
                Password for db access. Requires write access. Defaults to
                None, which means no authentication.
            collection:
                Collection to query. Defaults to "mol_tasks".
            simulate_mode:
                Allows one to simulate db insertion without actually performing
                the insertion.
        """
        self.host = host
        self.database = database
        self.user = user
        self.password = password
        self.collection = collection
        self.port = port
        self.simulate = simulate_mode
        self.update_duplicates = update_duplicates
        if not simulate_mode:
            conn = MongoClient(self.host, self.port, j=True)
            db = conn[self.database]
            if self.user:
                db.authenticate(self.user, self.password)
            if db.counter.find({"_id": "mol_taskid"}).count() == 0:
                db.counter.insert({"_id": "mol_taskid", "c": 1})

    def assimilate(self, path):
        """
        Parses nwchem runs. Then insert the result into the db. and return the
        mol_id or doc of the insertion.

        Returns:
            If in simulate_mode, the entire doc is returned for debugging
            purposes. Else, only the task_id of the inserted doc is returned.
        """
        try:
            d = self.get_task_doc(path)
            tid = self._insert_doc(d)
            return tid
        except Exception as ex:
            import traceback
            print traceback.format_exc(ex)
            logger.error(traceback.format_exc(ex))
            return False

    @classmethod
    def get_task_doc(cls, path):
        """
        Get the entire task doc for a path, including any post-processing.
        """
        logger.info("Getting task doc for base dir :{}".format(path))
        nwo = NwOutput(path)
        data = nwo.data
        mol = data[0]["molecules"][0]
        d = {"path": os.path.abspath(path),
             "pretty_formula": mol.composition.reduced_formula,
             "calculations": dict(zip(["GeomOpt", "Freq", "SCF", "IE_SCF",
                                       "EA_SCF"], data)),
             "EA": data[2]["energies"][-1] - data[3]["energies"][-1],
             "IE": data[4]["energies"][-1] - data[2]["energies"][-1]}
        return clean_json(d)

    def _insert_doc(self, d):
        if not self.simulate:
            # Perform actual insertion into db. Because db connections cannot
            # be pickled, every insertion needs to create a new connection
            # to the db.
            conn = MongoClient(self.host, self.port)
            db = conn[self.database]
            if self.user:
                db.authenticate(self.user, self.password)
            coll = db[self.collection]

            # Insert dos data into gridfs and then remove it from the dict.
            # DOS data tends to be above the 4Mb limit for mongo docs. A ref
            # to the dos file is in the dos_fs_id.
            result = coll.find_one({"path": d["path"]},
                                   fields=["path", "task_id"])
            if result is None or self.update_duplicates:
                d["last_updated"] = datetime.datetime.today()
                if result is None:
                    if ("task_id" not in d) or (not d["task_id"]):
                        d["task_id"] = db.counter.find_and_modify(
                            query={"_id": "mol_taskid"},
                            update={"$inc": {"c": 1}}
                        )["c"]
                    logger.info("Inserting {} with taskid = {}"
                                .format(d["path"], d["task_id"]))
                elif self.update_duplicates:
                    d["task_id"] = result["task_id"]
                    logger.info("Updating {} with taskid = {}"
                                .format(d["path"], d["task_id"]))

                coll.update({"path": d["path"]}, {"$set": d}, upsert=True)
                return d["task_id"]
            else:
                logger.info("Skipping duplicate {}".format(d["dir_name"]))
        else:
            d["task_id"] = 0
            logger.info("Simulated insert into database for {} with task_id {}"
                        .format(d["path"], d["task_id"]))
            return d

    def get_valid_paths(self, path):
        """
        There are some restrictions on the valid directory structures:

        1. There can be only one vasp run in each directory. Nested directories
           are fine.
        2. Directories designated "relax1", "relax2" are considered to be 2
           parts of an aflow style run.
        3. Directories containing vasp output with ".relax1" and ".relax2" are
           also considered as 2 parts of an aflow style run.
        """
        (parent, subdirs, files) = path

        return [os.path.join(parent, f) for f in files
                if f.endswith(".nwout")]

    def convert(self, d):
        return d

    def __str__(self):
        return self.__class__.__name__

    @classmethod
    def from_dict(cls, d):
        return cls(**d["init_args"])

    @property
    def to_dict(self):
        init_args = {"host": self.host, "port": self.port,
                     "database": self.database, "user": self.user,
                     "password": self.password,
                     "collection": self.collection,
                     "simulate_mode": self.simulate}
        output = {"name": self.__class__.__name__,
                  "init_args": init_args, "version": __version__}
        return output



import unittest


from pymatgen.apps.borg.queen import BorgQueen


test_dir = os.path.join(os.path.dirname(__file__), "..", "..",
                        'test_files')


class DeltaSCFNwChemToDbTaskDroneTest(unittest.TestCase):

    def test_assimilate(self):
        drone = DeltaSCFNwChemToDbTaskDrone(
            collection="mol_task_test")
        q = BorgQueen(drone)
        q.serial_assimilate(test_dir)


if __name__ == "__main__":
    unittest.main()