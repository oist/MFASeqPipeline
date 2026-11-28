# Mock of VeloGP's command line, see tests/mock/mocktool.py
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), *[".."] * 4))
import mocktool
mocktool.main("velogp")
