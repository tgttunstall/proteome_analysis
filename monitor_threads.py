#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 15:18:39 2024

@author: tunstall
"""
###############################################################################
#================
# without psutil
#================

import subprocess
import time

# path to the thread-spawning program
program = "python thread_test.py"

# start the subprocess
process = subprocess.Popen(program.split())

print("\nChecking no. of threads using linux proc")


# allow the subprocess to initialise
time.sleep(2)

try:
    # Read thread count from /proc/<pid>/status
    proc_status_path = f"/proc/{process.pid}/status"
    with open(proc_status_path, "r") as f:
        for line in f:
            if line.startswith("Threads:"):
                thread_count = int(line.split(":")[1].strip())
                break
    print(f"Number of threads used by the program: {thread_count}")
finally:
    # clean up the subprocess
    process.terminate()
    process.wait()
    print("Subprocess terminated.")
###############################################################################
#==============
# with psutil
#==============
import subprocess
import psutil
import time

# path to the thread-spawning program
program = "python thread_test.py"

# start the subprocess
process = subprocess.Popen(program.split())

print("\nChecking no. of threads using psutil")


# allow the subprocess to initialise
time.sleep(2)
try:
     # monitor thread count using psutil
     ps_process = psutil.Process(process.pid)
     thread_count = ps_process.num_threads()
     print(f"Number of threads used by the program: {thread_count}")
finally:
     # clean up the subprocess
     process.terminate()
     process.wait()
     print("Subprocess terminated.")

###############################################################################
