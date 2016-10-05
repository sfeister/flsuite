#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
helloMPI.py: A "Hello, world" for MPI4PY module.

Created by Scott Feister on Wed Oct 05 11:05:35 2016
CONTENT COPIED DIRECTLY FROM: https://github.com/erdc-cm/mpi4py/blob/master/demo/helloworld.py

Example call:
mpirun -np 4 python helloMPI.py

"""

from mpi4py import MPI
import sys

size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name()

sys.stdout.write(
    "Hello, World! I am process %d of %d on %s.\n"
    % (rank, size, name))