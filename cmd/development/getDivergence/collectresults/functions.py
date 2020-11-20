#!/usr/bin/env python3
import os
import tempfile

def mkdir_p(dir):
    'make a directory if doesnt exist'
    if not os.path.exists(dir):
        os.mkdir(dir)