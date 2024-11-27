#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 15:20:47 2024

@author: tunstall
"""

import threading
import time

# Define a simple thread worker function
def worker():
    while True:
        time.sleep(1)

# Create and start 5 threads
for i in range(5):
    thread = threading.Thread(target=worker, name=f"Thread-{i+1}")
    thread.start()

# Keep the main thread alive
while True:
    time.sleep(1)
    

