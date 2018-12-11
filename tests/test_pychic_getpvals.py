"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""
import os
import pytest # pylint: disable=unused-import
from CHiC.tool.pyCHiC import pyCHiC
import pandas as pd
import numpy as np

def test_pychic_getpvals():
   """
   Function to test the getPvals function from pyCHiC
   """
   input_getpvals = pd.read_csv("./data/getpvals_input.csv", sep="\t")
   output_getpvals = pd.read_csv("./data/getpvals_output.csv", sep="\t")

   dispersion = 2.5563913

   pyCHiC_object = pyCHiC()

   current_output = pyCHiC_object.getPvals(input_getpvals, dispersion)

   pd.testing.assert_frame_equal(output_getpvals, current_output)
