# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import json
import logging
import os
import shutil

from monty.io import zopen
from monty.os.path import zpath
from rubicon.workflows.qchem.wf_settings import EG_RUN_LOCS

__author__ = 'xiaohuiqu'


def get_eg_block_part(m_dir):
    parent_dir = os.path.basename(m_dir)
    grandpa_dir = os.path.basename(os.path.dirname(m_dir))
    return grandpa_dir, parent_dir


def get_eg_dir_loc(m_dir):
    if os.path.exists(m_dir):
        return m_dir

    all_locs = []
    if "GARDEN_LOC" in os.environ:
        all_locs.append(os.environ["GARDEN_LOC"])
    for scr_key in ["GSCRATCH", "SCRATCH", "SCRATCH2"]:
        for loc in EG_RUN_LOCS:
            if scr_key in os.environ:
                scr = os.environ[scr_key]
                all_locs.append(os.path.join(scr, loc))

    if "NERSC_HOST" in os.environ and os.environ["NERSC_HOST"] == "edison":
        all_locs.append("/scratch3/scratchdirs/xhqu/calculations/prod")
        all_locs.append("/scratch3/scratchdirs/xhqu/calculations/test")

    grandpa_dir, parent_dir = get_eg_block_part(m_dir)

    for preamble in all_locs:
        new_loc = os.path.join(preamble, grandpa_dir, parent_dir)
        if os.path.exists(new_loc):
            return new_loc

    raise ValueError('get_loc() -- dir does not exist!! Make sure your base '
                     'directory is listed in RUN_LOCS of wf_settings.py\n'
                     'and your garden location is set in the environment'
                     ' variable GARDEN_LOC')


def get_eg_file_loc(m_file):
    m_dir = os.path.dirname(m_file)
    filename = os.path.basename(m_file)
    dir_loc = get_eg_dir_loc(m_dir)
    return os.path.join(dir_loc, filename)


def move_to_eg_garden(m_dir):
    logger = logging.getLogger('MV_TO_GARDEN')
    if "GARDEN_LOC" not in os.environ:
        logger.info("GARDEN_LOC not available, nothing will be moved")
        return m_dir
    garden_part = os.path.abspath(os.environ["GARDEN_LOC"])
    if os.path.exists(m_dir) or os.path.exists(m_dir + ".gz"):
        if not os.path.isdir(m_dir):
            m_dir = os.path.dirname(m_dir)
    else:
        raise ValueError("The folder \"{}\" doesn't exist".format(m_dir))
    grandpa_dir, parent_dir = get_eg_block_part(m_dir)
    dest_dir = os.path.join(garden_part, grandpa_dir, parent_dir)
    dest_parent_dir = os.path.join(garden_part, grandpa_dir)
    if not os.path.exists(dest_parent_dir):
        # noinspection PyBroadException
        try:
            os.makedirs(dest_parent_dir)
        except:
            if not os.path.exists(dest_parent_dir):
                raise ValueError("Couldn't create parent folder \"{}\" in "
                                 "garden".format(dest_parent_dir))
    if os.path.exists(dest_dir):
        raise ValueError("A same name folder \"{}\" already exists in garden"
                         .format(dest_dir))
    logger.info("Trying to move folder {} to {}".format(m_dir, dest_dir))
    shutil.move(m_dir, dest_dir)
    logger.info("move to garden successfully completed")
    return dest_dir


def get_defuse_causing_qchem_fwid(qcout_path):
    dirname = os.path.dirname(qcout_path)
    fw_spec_path = os.path.join(dirname, "FW.json")
    with zopen(zpath(fw_spec_path), 'rt') as f:
        fw_dict = json.load(f)
    fw_id = fw_dict["fw_id"]
    return fw_id
