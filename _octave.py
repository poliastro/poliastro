import os.path

from oct2py import octave

OCTAVE_PATH = os.path.join(os.path.dirname(__file__), "octave")
octave.addpath(OCTAVE_PATH)
