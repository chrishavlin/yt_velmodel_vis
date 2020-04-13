"""
test script for plotting TF sequence

"""
from yt_velmodel_vis import animator as an
import numpy as np
import os
import matplotlib.pyplot as plt

# choose model and plot settings
a=an.TFsequence('cold_vol','cold_TF','./output/tf_tester/animate')
a.buildGif('./output/tf_tester/cold.gif')

a=an.TFsequence('hot_vol','hot_TF','./output/tf_tester/animate')
a.buildGif('./output/tf_tester/hot.gif')
