#!/bin/bash

alias dcc_interactive="ssh -L 8090:localhost:8090 cjp47@dcc-login.oit.duke.edu | tee /dev/tty | python3 ~/.ssh/run_webbrowser.py"

echo "import sys
import webbrowser

for line in sys.stdin:
    if "OPEN_ON_LOCAL[" in line:
        line = line.replace("0.0.0.0", "127.0.0.1")
        webbrowser.open(line.split("OPEN_ON_LOCAL[")[1].split("]")[0])" > ~/.ssh/run_webbrowser.py