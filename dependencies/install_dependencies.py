#!/usr/bin/env python3
'''
Install dependencies
'''

import os
import stat
import platform
import subprocess


if __name__ == '__main__':

  # Download DSSP

  dssp_bin = 'dssp-2.0.4-linux-i386'

  system = platform.system()

  if system == 'Windows':
    dssp_bin = 'dssp-2.0.4-win32.exe'

  if not os.path.exists(dssp_bin):
    subprocess.check_call(['wget', 'ftp://ftp.cmbi.ru.nl/pub/software/dssp/' + dssp_bin])
  
  st = os.stat(dssp_bin)
  os.chmod(dssp_bin, st.st_mode | stat.S_IEXEC)

