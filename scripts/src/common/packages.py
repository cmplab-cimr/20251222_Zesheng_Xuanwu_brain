import itertools as it
import warnings
from collections import defaultdict, abc
import gzip
import pickle
import pathlib
import os
import copy
import multiprocessing as mp
import argparse
import enum
import datetime


class ValueEnum(enum.Enum):
    def __str__(self):
        return self.value

    def startswith(self, substr):
        return self.value.startswith(substr)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import manifold
from sklearn import decomposition
import tqdm
import xlsxwriter
import openpyxl
from scipy.stats import linregress
