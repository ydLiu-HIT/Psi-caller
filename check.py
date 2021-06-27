#!/usr/bin/env python
# coding=utf-8

import os
import sys
from subprocess import check_output, PIPE, Popen

def CheckFileExist(fn, sfx=""):
    if not os.path.isfile(fn+sfx):
        sys.exit("Error: %s not found" %(fn+sfx))
    return os.path.abspath(fn+sfx)


def CheckCmdExist(cmd):
    if not isinstance(cmd, str):
        return False

    try:
        check_output("which %s" % (cmd), shell=True)
    except:
        sys.exit("Error: %s executable not found" %(cmd))
    return cmd

def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=sys.stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)
