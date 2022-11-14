import sys
import os
import numpy as np
import copy as cp
import math as m
from alive_progress import alive_bar, config_handler
import random
import pandas as pd
import time

import matplotlib.pyplot as plt
from matplotlib import rc

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../../'))
sys.path.append(os.path.join(current_path, '../../../netgrid'))
sys.path.append(os.path.join(current_path, '../../../raplea'))

import hybridPNM as pnm

import raplea as rl

start = time.time()

from netgrid import Netgrid, saveFilesCollectionToFile

bar_type = random.choice(['smooth', 'classic', 'brackets', 'blocks', 'bubbles', 'solid', 'checks',
                          'circles', 'squares', 'ruler', 'ruler2', 'fish', 'scuba'])
spinner_type = random.choice(['classic', 'stars', 'twirl', 'twirls', 'horizontal', 'vertical',
                              'waves', 'waves2', 'waves3', 'dots', 'triangles', 'brackets',
                              'arrow', 'arrows', 'arrows2', 'arrows_in', 'arrows_out', 'radioactive'])
config_handler.set_global(bar=bar_type, spinner=spinner_type)


def smooth_bc(x_ini, x_bc, time, time_curr):
    if time_curr > time:
        return x_bc
    else:
        return x_ini + (x_bc - x_ini) * time_curr / time


rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"
# plt.style.use('dracula')

np.set_printoptions(threshold=np.inf, linewidth=300)

# ==========================================================================================
#                                       input                                              #
# ==========================================================================================

dt = 5.e-2
t_max = 100.
t_bc = t_max
t_step_out_coarse = t_max / 100
t_step_out_fine = t_max / 1000


visc0 = 1.e-3
visc1 = 2.e-5
gamma = 1
beta = 7.2e-2
theta = (m.pi / 180.) * 30.
xi = 0.

P_ini = 3.0e6
PBcTypes = {"inlet": "Dirichlet", "outlet": "Dirichlet"}
valuePInlet = 4.e6
valuePInletIni = P_ini
valuePOutlet = 4.e6
valuePOutletIni = P_ini

S0_ini = 1.
S0BcTypes = {"inlet": "Neumann", "outlet": "Neumann"}
valueS0Inlet = 1.
valueS0Outlet = 1.

scaleM = 0.5
scaleM *= 0.4
phiM = 0.087
kM = 4. * 1.e-17

S0RM = 0.0
S1RM = 0.0

P0M_ini = P_ini
S0M_ini = 0.

round_output_time = 4

# ==========================================================================================
#                                    netgrid for pnm                                       #
# ==========================================================================================

voxelSize = 0.000028

VTotalFractures = 3.3487776e-10

poresCoordinates = {0: [0.00253672, 0.], 1: [0.00209, 0.00164], 2: [0.00185, 0.00195], 3: [0.00272, 0.],
                    4: [0.00297467, 0.000560006], 5: [0.00532, 0.0010158], 6: [0.00563499, 0.000946919],
                    7: [0.0058464, 0.00111452],
                    8: [0.00728, 0.00117005], 9: [0.00518, 0.001925], 10: [0.00532061, 0.00230941],
                    11: [0.00653, 0.0025], 12: [0.00728, 0.00266],
                    13: [0.0017866, 0.00231], 14: [0.0015792, 0.00259], 15: [0.001538, 0.0028], 16: [0.00143, 0.003115],
                    17: [0.0016462, 0.0032008],
                    18: [0.001227, 0.0036233], 19: [0.00119, 0.00381], 20: [0.00112, 0.00392], 21: [0.00094, 0.00408],
                    22: [0.0007, 0.00408], 23: [0.00049, 0.00425],
                    24: [0, 0.00414], 25: [0.00102, 0.00427], 26: [0.00112, 0.00427], 27: [0.00119, 0.00407],
                    28: [0.0021, 0.0021], 29: [0.00214, 0.00228], 30: [0.0021, 0.00255], 31: [0.00208, 0.0027],
                    32: [0.00196, 0.0031], 33: [0.0018, 0.0035],
                    34: [0.0017, 0.0039],
                    35: [0.002, 0.00357], 36: [0.0021, 0.00365], 37: [0.0023, 0.00372], 38: [0.0024, 0.00375],
                    39: [0.00257, 0.00381], 40: [0.00285, 0.00387],
                    41: [0.0036, 0.004021], 42: [0.0044, 0.0041], 43: [0.005, 0.00418], 44: [0.0058, 0.0043],
                    45: [0.0065, 0.00448], 46: [0.0069, 0.00452], 47: [0.007025, 0.0048],
                    48: [0.00728, 0.0049], 49: [0.00224, 0.00187], 50: [0.0020, 0.00392], 51: [0.00215, 0.00425],
                    52: [0.00259, 0.004265], 53: [0.00276, 0.00436], 54: [0.00292, 0.00453],
                    55: [0.00329, 0.00453], 56: [0.00378, 0.00490], 57: [0.00392, 0.00515], 58: [0.0042, 0.00528],
                    59: [0.00493, 0.00539],
                    60: [0.00728, 0.00536], 61: [0.00399, 0.00555], 62: [0.00383, 0.00619], 63: [0.00381, 0.00630],
                    64: [0.00369, 0.00632], 65: [0.00319, 0.00634],
                    66: [0.00231, 0.00606], 67: [0.00247, 0.00538], 68: [0.00143, 0.00490], 69: [0.00144, 0.00455],
                    70: [0.00161, 0.00424],
                    71: [0.00128, 0.00529], 72: [0.00112, 0.00574], 73: [0.00119, 0.00658], 74: [0.00114, 0.00732],
                    75: [0.00105, 0.00826],
                    76: [0.00099, 0.00854], 77: [0.00091, 0.00917], 78: [0.00070, 0.00970], 79: [0.00055, 0.00985],
                    80: [0., 0.00971], 81: [0.00390, 0.00637],
                    82: [0.00408, 0.00643], 83: [0.00430, 0.00654], 84: [0.00446, 0.00667], 85: [0.004545, 0.00684],
                    86: [0.00446, 0.00701],
                    87: [0.00383, 0.00658], 88: [0.00421, 0.00721], 89: [0.00478, 0.00772], 90: [0.00529, 0.00778],
                    91: [0.0055, 0.00788],
                    92: [0.00602, 0.00789], 93: [0.00697, 0.00801], 94: [0.00728, 0.00836], 95: [0.00728, 0.00781],
                    96: [0.00499, 0.00830], 97: [0.00469, 0.00843],
                    98: [0.00455, 0.00884], 99: [0.00434, 0.00931], 100: [0.00413, 0.00963], 101: [0.00388, 0.00987],
                    102: [0.00396, 0.01008], 103: [0.00506, 0.01012], 104: [0.00585, 0.01022], 105: [0.00591, 0.01003],
                    106: [0.00686, 0.01012], 107: [0.00728, 0.01021],
                    108: [0.00332, 0.01005], 109: [0.00268, 0.00998], 110: [0.00234, 0.00991], 111: [0.00208, 0.0102],
                    112: [0.00122, 0.01012], 113: [0.00108, 0.00994], 114: [0.00109, 0.00973], 115: [0.00087, 0.00973],
                    116: [0.00117, 0.01045], 117: [0.00046, 0.01012],
                    118: [0.00028, 0.01050], 119: [0.00029, 0.01075], 120: [0.00092, 0.01158], 121: [0.00076, 0.0119],
                    122: [0.00119, 0.01365], 123: [0.00105, 0.01425], 124: [0.00133, 0.01477], 125: [0.00081, 0.01551],
                    126: [0.00098, 0.01605], 127: [0.00000, 0.01718],
                    128: [0.00387, 0.01027], 129: [0.00372, 0.01055], 130: [0.00375, 0.01096], 131: [0.00469, 0.01089],
                    132: [0.00473, 0.01078], 133: [0.00493, 0.01067], 134: [0.00584, 0.01071], 135: [0.00728, 0.01065],
                    136: [0.00172, 0.01484], 137: [0.00237, 0.01462], 138: [0.00287, 0.01455], 139: [0.00340, 0.01414],
                    140: [0.00525, 0.01302], 141: [0.0051, 0.01258], 142: [0.00517, 0.01172], 143: [0.00602, 0.01177],
                    144: [0.00703, 0.01183], 145: [0.00728, 0.01177], 146: [0.00047, 0.01596],
                    147: [0.00127, 0.01144], 148: [0.00131, 0.01116], 149: [0.00158, 0.01066], 150: [0.00020, 0.01092],
                    151: [0.00035, 0.01126],
                    152: [0.00277, 0.01421], 153: [0.00311, 0.01365], 154: [0.00360, 0.01293], 155: [0.00360, 0.01215],
                    156: [0.00388, 0.01183],
                    157: [0.00411, 0.01160], 158: [0.00589, 0.01134],

                    159: [0.00686, 0.01197], 160: [0.00600, 0.01267],
                    161: [0.00161, 0.01512], 162: [0.00119, 0.01571],
                    163: [0.00054, 0.01158]
                    }

throatsPores = {0: [0, 1], 1: [1, 2], 2: [3, 4], 3: [4, 5], 4: [5, 6], 5: [6, 7], 6: [7, 8], 7: [5, 9], 8: [9, 10],
                9: [10, 11], 10: [11, 12], 11: [2, 13], 12: [13, 14], 13: [14, 15], 14: [15, 16], 15: [16, 17],
                16: [16, 18], 17: [18, 19],
                18: [19, 20], 19: [20, 21], 20: [21, 22], 21: [22, 23], 22: [23, 24], 23: [21, 25], 24: [25, 26],
                25: [26, 27], 26: [27, 20],
                27: [2, 28], 28: [28, 29], 29: [29, 30], 30: [30, 31], 31: [31, 32], 32: [32, 33], 33: [33, 34],
                34: [34, 19],
                35: [33, 35], 36: [35, 36], 37: [36, 37], 38: [37, 38], 39: [38, 39], 40: [39, 40], 41: [40, 41],
                42: [41, 42], 43: [42, 43],
                44: [43, 44], 45: [44, 45], 46: [45, 46], 47: [46, 47], 48: [47, 48], 49: [28, 49], 50: [34, 50],
                51: [50, 51], 52: [51, 52],
                53: [52, 53], 54: [53, 54], 55: [54, 55], 56: [55, 56], 57: [56, 57], 58: [57, 58], 59: [58, 59],
                60: [59, 60],
                61: [58, 61], 62: [61, 62], 63: [62, 63], 64: [63, 64], 65: [64, 65], 66: [65, 66], 67: [66, 67],
                68: [67, 68], 69: [68, 69],
                70: [69, 70], 71: [70, 34], 72: [68, 71], 73: [71, 72], 74: [72, 73], 75: [73, 74], 76: [74, 75],
                77: [75, 76], 78: [76, 77], 79: [77, 78],
                80: [78, 79], 81: [79, 80], 82: [63, 81], 83: [81, 82], 84: [82, 83], 85: [83, 84], 86: [84, 85],
                87: [85, 86],
                88: [63, 87], 89: [87, 88], 90: [88, 89], 91: [89, 90], 92: [90, 91], 93: [91, 92], 94: [92, 93],
                95: [93, 94], 96: [93, 95],
                97: [91, 96], 98: [96, 97], 99: [97, 98], 100: [98, 99], 101: [99, 100], 102: [100, 101],
                103: [101, 102], 104: [102, 103],
                105: [103, 104], 106: [104, 105], 107: [105, 106], 108: [106, 107], 109: [102, 108], 110: [108, 109],
                111: [109, 110], 112: [110, 111],
                113: [111, 112], 114: [112, 113], 115: [113, 114], 116: [114, 115], 117: [115, 78], 118: [112, 116],
                119: [79, 117], 120: [117, 118],
                121: [118, 119], 122: [119, 120], 123: [120, 121], 124: [121, 122], 125: [122, 123], 126: [123, 124],
                127: [124, 125], 128: [125, 126], 129: [126, 127],
                130: [102, 128], 131: [128, 129], 132: [129, 130], 133: [130, 131], 134: [131, 132], 135: [132, 133],
                136: [133, 134],
                137: [134, 135], 138: [124, 136], 139: [136, 137],
                140: [137, 138], 141: [138, 139], 142: [139, 140], 143: [140, 141], 144: [141, 142], 145: [142, 143],
                146: [143, 144], 147: [144, 145], 148: [125, 146],
                149: [120, 147], 150: [147, 148], 151: [148, 149], 152: [119, 150], 153: [150, 151],
                154: [137, 152], 155: [152, 153], 156: [153, 154], 157: [154, 155], 158: [155, 156], 159: [156, 157],
                160: [143, 158], 161: [158, 134],

                162: [144, 159], 163: [159, 160], 164: [136, 161], 165: [161, 162], 166: [121, 163],
                167: [0, 3]
                }

throatsWidths = {0: 0.00028, 1: 0.00019, 2: 0.0002, 3: 0.000145, 4: 0.00014, 5: 0.00012, 6: 0.0001, 7: 0.00015,
                 8: 0.00021, 9: 0.00021, 10: 0.00017, 11: 0.00019, 12: 0.00021, 13: 0.00014, 14: 0.000185, 15: 0.00012,
                 16: 0.000165,
                 17: 0.0001, 18: 0.000152, 19: 0.000115, 20: 0.000105, 21: 0.000115, 22: 0.000125, 23: 0.0001,
                 24: 0.0002, 25: 0.00007, 26: 0.0001,
                 27: 0.0001, 28: 0.0001, 29: 0.0001, 30: 0.0001, 31: 0.0001, 32: 0.0001, 33: 0.0001, 34: 0.0001,
                 35: 0.0001, 36: 0.0001, 37: 0.0001, 38: 0.0001, 39: 0.0001, 40: 0.0001, 41: 0.0001, 42: 0.0001,
                 43: 0.0001,
                 44: 0.0001, 45: 0.0001, 46: 0.0001, 47: 0.0001, 48: 0.0001, 49: 0.0001, 50: 0.0001, 51: 0.0001,
                 52: 0.0001, 53: 0.0001, 54: 0.0001,
                 55: 0.0001, 56: 0.0001, 57: 0.0001, 58: 0.0001, 59: 0.0001, 60: 0.0001, 61: 0.0001, 62: 0.0001,
                 63: 0.0001,
                 64: 0.0001, 65: 0.0001, 66: 0.0001, 67: 0.0001, 68: 0.0001, 69: 0.0001, 70: 0.0001, 71: 0.0001,
                 72: 0.0001, 73: 0.0001,
                 74: 0.0001, 75: 0.0001, 76: 0.0001, 77: 0.0001, 78: 0.0001, 79: 0.0001, 80: 0.0001, 81: 0.0001,
                 82: 0.0001, 83: 0.0001,
                 84: 0.0001, 85: 0.0001, 86: 0.0001, 87: 0.0001, 88: 0.0001, 89: 0.0001, 90: 0.0001, 91: 0.0001,
                 92: 0.0001, 93: 0.0001,
                 94: 0.0001, 95: 0.0001, 96: 0.0001, 97: 0.0001, 98: 0.0001, 99: 0.0001, 100: 0.0001, 101: 0.0001,
                 102: 0.0001, 103: 0.0001,
                 104: 0.0001, 105: 0.0001, 106: 0.0001, 107: 0.0001, 108: 0.0001, 109: 0.0001, 110: 0.0001,
                 111: 0.0001, 112: 0.0001, 113: 0.0001, 114: 0.0001,
                 115: 0.0001, 116: 0.0001, 117: 0.0001, 118: 0.0001, 119: 0.0001, 120: 0.0001, 121: 0.0001,
                 122: 0.0001, 123: 0.0001, 124: 0.0001, 125: 0.0001,
                 126: 0.0001, 127: 0.0001, 128: 0.0001, 129: 0.0001, 130: 0.0001, 131: 0.0001, 132: 0.0001,
                 133: 0.0001, 134: 0.0001, 135: 0.0001, 136: 0.0001, 137: 0.0001, 138: 0.0001, 139: 0.0001, 140: 0.0001,
                 141: 0.0001,
                 142: 0.0001, 143: 0.0001, 144: 0.0001, 145: 0.0001,
                 146: 0.0001, 147: 0.0001, 148: 0.0001, 149: 0.0001, 150: 0.0001, 151: 0.0001, 152: 0.0001,
                 153: 0.0001, 154: 0.0001, 155: 0.0001, 156: 0.0001, 157: 0.0001, 158: 0.0001, 159: 0.0001, 160: 0.0001,
                 161: 0.0001, 162: 0.0001, 163: 0.0001,
                 164: 0.0001, 165: 0.0001, 166: 0.0001, 167: 0.0001, 168: 0.0001, 169: 0.0001, 170: 0.0001}

throatsDepths = {0: voxelSize, 1: voxelSize, 2: voxelSize, 3: voxelSize, 4: voxelSize, 5: voxelSize, 6: voxelSize, 7: voxelSize,
                 8: voxelSize, 9: voxelSize, 10: voxelSize, 11: voxelSize, 12: voxelSize, 13: voxelSize, 14: voxelSize, 15: voxelSize,
                 16: voxelSize,
                 17: voxelSize, 18: voxelSize, 19: voxelSize, 20: voxelSize, 21: voxelSize, 22: voxelSize, 23: voxelSize, 24: voxelSize,
                 25: voxelSize,
                 26: voxelSize, 27: voxelSize, 28: voxelSize, 29: voxelSize, 30: voxelSize, 31: voxelSize, 32: voxelSize, 33: voxelSize,
                 34: voxelSize,
                 35: voxelSize, 36: voxelSize, 37: voxelSize, 38: voxelSize, 39: voxelSize, 40: voxelSize, 41: voxelSize, 42: voxelSize,
                 43: voxelSize,
                 44: voxelSize, 45: voxelSize, 46: voxelSize, 47: voxelSize, 48: voxelSize, 49: voxelSize, 50: voxelSize, 51: voxelSize,
                 52: voxelSize, 53: voxelSize, 54: voxelSize,
                 55: voxelSize, 56: voxelSize, 57: voxelSize, 58: voxelSize, 59: voxelSize, 60: voxelSize, 61: voxelSize, 62: voxelSize,
                 63: voxelSize,
                 64: voxelSize, 65: voxelSize, 66: voxelSize, 67: voxelSize, 68: voxelSize, 69: voxelSize, 70: voxelSize, 71: voxelSize,
                 72: voxelSize, 73: voxelSize,
                 74: voxelSize, 75: voxelSize, 76: voxelSize, 77: voxelSize, 78: voxelSize, 79: voxelSize, 80: voxelSize, 81: voxelSize,
                 82: voxelSize, 83: voxelSize,
                 84: voxelSize, 85: voxelSize, 86: voxelSize, 87: voxelSize, 88: voxelSize, 89: voxelSize, 90: voxelSize, 91: voxelSize,
                 92: voxelSize, 93: voxelSize,
                 94: voxelSize, 95: voxelSize, 96: voxelSize, 97: voxelSize, 98: voxelSize, 99: voxelSize, 100: voxelSize,
                 101: voxelSize, 102: voxelSize, 103: voxelSize,
                 104: voxelSize, 105: voxelSize, 106: voxelSize, 107: voxelSize, 108: voxelSize, 109: voxelSize, 110: voxelSize,
                 111: voxelSize, 112: voxelSize, 113: voxelSize, 114: voxelSize,
                 115: voxelSize, 116: voxelSize, 117: voxelSize, 118: voxelSize, 119: voxelSize, 120: voxelSize, 121: voxelSize,
                 122: voxelSize, 123: voxelSize, 124: voxelSize, 125: voxelSize,
                 126: voxelSize, 127: voxelSize, 128: voxelSize, 129: voxelSize, 130: voxelSize, 131: voxelSize, 132: voxelSize,
                 133: voxelSize, 134: voxelSize, 135: voxelSize, 136: voxelSize, 137: voxelSize, 138: voxelSize, 139: voxelSize,
                 140: voxelSize,
                 141: voxelSize,
                 142: voxelSize, 143: voxelSize, 144: voxelSize, 145: voxelSize,
                 146: voxelSize, 147: voxelSize, 148: voxelSize, 149: voxelSize, 150: voxelSize, 151: voxelSize, 152: voxelSize,
                 153: voxelSize, 154: voxelSize, 155: voxelSize, 156: voxelSize, 157: voxelSize, 158: voxelSize, 159: voxelSize,
                 160: voxelSize,
                 161: voxelSize, 162: voxelSize, 163: voxelSize,
                 164: voxelSize, 165: voxelSize, 166: voxelSize, 167: voxelSize, 168: voxelSize, 169: voxelSize, 170: voxelSize
                 }

throatsWidthsChange1 = {0: 1., 1: 1., 2: 1., 3: 1., 4: 1., 5: 1., 6: 1., 7: 1., 8: 0.95, 9: 1., 10: 0.9, 11: 0.97,
                        12: 0.97,
                        13: 0.95, 14: 1., 15: 1., 16: 1.07, 17: 1., 18: 1., 19: 1., 20: 1., 21: 1., 22: 1., 23: 1.,
                        24: 1.,
                        25: 1., 26: 1., 27: 1.38, 28: 1.85, 29: 1.33, 30: 1.2, 31: 1.5, 32: 1.2, 33: 2.3, 34: 1.15,
                        35: 1.2, 36: 0.7,
                        37: 0.7, 38: 0.7, 39: 0.8, 40: 0.65, 41: 0.9, 42: 0.9, 43: 0.9, 44: 0.9, 45: 1., 46: 1., 47: 1.,
                        48: 1.,
                        49: 1.34, 50: 1.1, 51: 1.48, 52: 1.45, 53: 1., 54: 1.45, 55: 1., 56: 1.15, 57: 1.2, 58: 1.15,
                        59: 1.1,
                        60: 1.25,
                        61: 1.2, 62: 0.75, 63: 0.9, 64: 1., 65: 0.65, 66: 1., 67: 1., 68: 1., 69: 1.7, 70: 2., 71: 1.65,
                        72: 1.8,
                        73: 1.35, 74: 2.6, 75: 2., 76: 2.2, 77: 2.2, 78: 2.2, 79: 2.3, 80: 1.75, 81: 1., 82: 0.8,
                        83: 0.73,
                        84: 1.13,
                        85: 0.87, 86: 1.15, 87: 1., 88: 1., 89: 1.1, 90: 1.27, 91: 0.95, 92: 1., 93: 1.1, 94: 1.,
                        95: 0.8, 96: 0.6,
                        97: 1.17, 98: 1., 99: 1.2, 100: 1.15, 101: 1.1, 102: 1.1, 103: 1.11, 104: 0.7, 105: 0.85,
                        106: 1.14,
                        107: 0.8,
                        108: 0.9, 109: 0.7, 110: 0.8, 111: 0.8, 112: 0.8, 113: 0.87, 114: 1.15, 115: 1.15, 116: 0.9,
                        117: 2.5,
                        118: 0.9, 119: 1.5, 120: 1.48, 121: 1.55, 122: 1.1, 123: 1.65, 124: 1.9, 125: 2.15, 126: 1.65,
                        127: 2.0,
                        128: 1.2, 129: 2.2, 130: 1.3, 131: 0.9, 132: 1., 133: 1., 134: 1.2, 135: 0.7, 136: 1., 137: 1.2,
                        138: 1.15, 139: 2.5, 140: 1.3, 141: 1., 142: 1.3, 143: 1., 144: 1., 145: 1., 146: 0.75, 147: 1.,
                        148: 1.1, 149: 0.8, 150: 1., 151: 1., 152: 1., 153: 0.9, 154: 1.1, 155: 1.2, 156: 1.2, 157: 1.2,
                        158: 1.5, 159: 0.8, 160: 1.15, 161: 1.2, 162: 1.2, 163: 0.8, 164: 1.3, 165: 1., 166: 1., 167: 1.5,
                        168: 1., 169: 1., 170: 1.}

throatsWidthsChange2 = {0: 1., 1: 1., 2: 1., 3: 1., 4: 1., 5: 1., 6: 1., 7: 1., 8: 1., 9: 1., 10: 1., 11: 1., 12: 1.,
                        13: 1., 14: 1., 15: 1., 16: 1., 17: 1., 18: 1., 19: 1., 20: 1., 21: 1., 22: 1., 23: 1., 24: 1.,
                        25: 1., 26: 1., 27: 1., 28: 1., 29: 1., 30: 1., 31: 1., 32: 1., 33: 1., 34: 1., 35: 1., 36: 1.,
                        37: 1., 38: 1., 39: 1., 40: 1., 41: 1., 42: 1., 43: 1., 44: 1., 45: 1., 46: 1., 47: 1., 48: 1.,
                        49: 1., 50: 1., 51: 1., 52: 1., 53: 1., 54: 1., 55: 1., 56: 1., 57: 1., 58: 1., 59: 1., 60: 1.,
                        61: 1., 62: 1., 63: 1., 64: 1., 65: 1., 66: 1., 67: 1., 68: 1., 69: 1., 70: 1., 71: 1., 72: 1.,
                        73: 1., 74: 1., 75: 1., 76: 1., 77: 1., 78: 1., 79: 1., 80: 1., 81: 1., 82: 1., 83: 1., 84: 1.,
                        85: 1., 86: 1., 87: 1., 88: 1., 89: 1., 90: 1., 91: 1., 92: 1., 93: 1., 94: 1., 95: 1., 96: 1.,
                        97: 1., 98: 1., 99: 1., 100: 1., 101: 1., 102: 1., 103: 1., 104: 1., 105: 1., 106: 1., 107: 1.,
                        108: 1., 109: 1., 110: 1., 111: 1., 112: 1., 113: 1., 114: 1., 115: 1., 116: 1., 117: 1.,
                        118: 1., 119: 1., 120: 1., 121: 1., 122: 1., 123: 1., 124: 1., 125: 1., 126: 1., 127: 1.,
                        128: 1., 129: 1., 130: 1., 131: 1., 132: 1., 133: 1., 134: 1., 135: 1., 136: 1., 137: 1.,
                        138: 1., 139: 1., 140: 1., 141: 1., 142: 1., 143: 1., 144: 1., 145: 1., 146: 1., 147: 1.,
                        148: 1., 149: 1., 150: 1., 151: 1., 152: 1., 153: 1., 154: 1., 155: 1., 156: 1., 157: 1.,
                        158: 1., 159: 1., 160: 1., 161: 1., 162: 1., 163: 1., 164: 1., 165: 1., 166: 1., 167: 1.,
                        168: 1., 169: 1., 170: 1.}

# Throats alteration No1
for key, value in throatsWidths.items():
    throatsWidths[key] = value * throatsWidthsChange1[key]
# Throats alteration No2
for key, value in throatsWidths.items():
    throatsWidths[key] = value * throatsWidthsChange2[key]    


  
deltaV = 1e-13
minCellsN = 10
inletPores = {24, 80, 127}
outletPores = {145, 135, 107, 94, 95, 60, 48, 12, 8}


netgrid = Netgrid(poresCoordinates, throatsPores,
                  throatsWidths, throatsDepths, deltaV, minCellsN,
                  inletPores, outletPores)
                  
VTotalFracturesCurr = 0.
for throat in throatsPores:
    VTotalFracturesCurr += netgrid._throatsLs[throat] * netgrid._throatsAs[throat]
    
for throat in throatsPores:
    throatsWidths[throat] *= VTotalFractures / VTotalFracturesCurr
    
                  
netgrid = Netgrid(poresCoordinates, throatsPores,
                  throatsWidths, throatsDepths, deltaV, minCellsN,
                  inletPores, outletPores)
# print('PN cellsN:', netgrid._cellsN)

# ==========================================================================================
#                                    netgrid for matrix                                    #
# ==========================================================================================


poresCoordinatesM = {}
throatsPoresM = {}
throatsWidthsM = {}
throatsDepthsM = {}

deltaVM = scaleM * deltaV / phiM
deltaVM = scaleM * deltaV
deltaVM = deltaVM
minCellsNM = 20
inletPoresM = set()
outletPoresM = set()

VTotalMatrix = 3.2037846399999994e-09
VTotalMatrix *= 1.25 # temporary

fracturesLength = 0.
for throat in throatsPores:
    fracturesLength += netgrid._throatsLs[throat]

throatsVM = {}
for throat in throatsPores:
    throatsVM[throat] = netgrid._throatsLs[throat] * VTotalMatrix / fracturesLength

poreM = 0
for throat, pores in throatsPores.items():
    width = throatsWidths[throat]
    depth = throatsDepths[throat]
    L = netgrid._throatsLs[throat]
    LM = throatsVM[throat] / depth / L * scaleM
    pore_a = pores[0]
    pore_b = pores[1]
    coordinates_a = poresCoordinates[pore_a]
    coordinates_b = poresCoordinates[pore_b]
    alpha = -m.atan2(coordinates_b[0] - coordinates_a[0], coordinates_b[1] - coordinates_a[1])
    cellsN = netgrid._throatsCellsNs[throat]
    countCell = 0
    for cell in netgrid._throatsCells[throat]:
        point_a = [coordinates_a[0] + (coordinates_b[0] - coordinates_a[0]) * countCell / cellsN,
                   coordinates_a[1] + (coordinates_b[1] - coordinates_a[1]) * countCell / cellsN]
        point_b = [coordinates_a[0] + (coordinates_b[0] - coordinates_a[0]) * (countCell + 1) / cellsN,
                   coordinates_a[1] + (coordinates_b[1] - coordinates_a[1]) * (countCell + 1) / cellsN]
        point_c = [(point_a[0] + point_b[0]) / 2. + m.cos(alpha) * width / 2.,
                   (point_a[1] + point_b[1]) / 2. + m.sin(alpha) * width / 2.]
        point_pc = [point_c[0] + m.cos(alpha) * LM,
                    point_c[1] + m.sin(alpha) * LM]

        poresCoordinatesM[poreM] = point_pc
        poresCoordinatesM[poreM + 1] = point_c
        throatsPoresM[cell] = [poreM, poreM + 1]
        throatsWidthsM[cell] = L / cellsN
        throatsDepthsM[cell] = depth
        inletPoresM.add(poreM)
        outletPoresM.add(poreM + 1)

        poreM += 2
        countCell += 1

  
netgridM = Netgrid(poresCoordinatesM, throatsPoresM,
                   throatsWidthsM, throatsDepthsM, deltaVM, minCellsNM,
                   inletPoresM, outletPoresM)
# print('Matrix cellsN:', netgridM._cellsN)
# ==========================================================================================
#                                        raplea                                            #
# ==========================================================================================


alphaCoeffs = [0, 0, 0, 4.3549791928771896e-06, -10.161618116713441]

# temporary
# for i in range(len(alphaCoeffs)):
#    alphaCoeffs[i] *= 1.0122
    

pars = rl.Pars(S0R=S0RM, S1R=S1RM, M0=1., M1=1., MC=0., m0=1., m1=1., mC=1., mu0=visc0, mu1=visc1,
               alphaCoeffs=alphaCoeffs, sorptionType='kR1Linear', is0Wet=True, dt=dt)

alpha1M_ini = pars.alpha1(P0M_ini)

throatsNM = netgridM._throatsN
facesNM = netgridM._facesN
cellsNM = netgridM._cellsN
uM = [np.zeros(facesNM, dtype=np.float64),
      np.zeros(facesNM, dtype=np.float64)]
P0M = [np.full(cellsNM, P0M_ini, dtype=np.float64),
       np.full(cellsNM, P0M_ini, dtype=np.float64)]
S0M = [np.full(cellsNM, S0M_ini, dtype=np.float64),
       np.full(cellsNM, S0M_ini, dtype=np.float64)]
alpha1M = [np.full(cellsNM, alpha1M_ini, dtype=np.float64),
           np.full(cellsNM, alpha1M_ini, dtype=np.float64)]

data = {"uM": np.zeros(throatsNM, dtype=np.float64),
        "qM": np.zeros(throatsNM, dtype=np.float64),
        "S0M": np.zeros(throatsNM, dtype=np.float64),
        "fM": np.zeros(throatsNM, dtype=np.float64),
        "deltaAlpha1M": np.zeros(throatsNM, dtype=np.float64),
        "avAlpha1M": np.full(throatsNM, alpha1M_ini, dtype=np.float64),
        "deltaS0M": np.zeros(throatsNM, dtype=np.float64),
        "avS0M": np.full(throatsNM, S0M_ini, dtype=np.float64),
        "S0": np.full(throatsNM, S0_ini, dtype=np.float64),
        "P": np.full(throatsNM, P_ini, dtype=np.float64),
        "KFactor": np.ones(throatsNM, dtype=np.float64)}

matrix = rl.Matrix(phiM, kM, scaleM, netgridM, pars, data, uM, P0M, S0M, alpha1M)
matrix.setSolverPs('sparseLU')
matrix.setSolverSs('sparseLU')

iPrevM = matrix._iPrev
uM_array_times = [np.copy(uM[iPrevM])]
uM_av_times = [np.mean(uM[iPrevM])]
P0M_array_times = [np.copy(P0M[iPrevM])]
P0M_av_times = [np.mean(P0M[iPrevM])]
S0M_array_times = [np.copy(S0M[iPrevM])]
S0M_av_times = [np.mean(S0M[iPrevM])]
alpha1M_array_times = [np.copy(alpha1M[iPrevM])]
alpha1M_av_times = [np.mean(alpha1M[iPrevM])]

# ==========================================================================================
#                                     hybridPNM                                            #
# ==========================================================================================


props = pnm.Props({"dt": float(dt),
                   "visc0": float(visc0), "visc1": float(visc0),
                   "theta": float(theta), "gamma": float(gamma), "beta": float(beta),
                   "dim": "2D", "xi": float(xi)})

u = np.zeros(netgrid._facesN, dtype=np.float64)
dS0 = np.zeros(netgrid._facesN, dtype=np.float64)
P = [np.full(netgrid._cellsN, P_ini, dtype=np.float64), np.full(netgrid._cellsN, P_ini, dtype=np.float64)]
S0 = [np.full(netgrid._cellsN, S0_ini, dtype=np.float64), np.full(netgrid._cellsN, S0_ini, dtype=np.float64)]

locP = pnm.LocP(netgrid, props, u, dS0, P, S0, data)
conP = pnm.ConP(netgrid, props, u, dS0, P, S0, data)
equationP = pnm.Equation(locP, conP, P)
equationP._bcTypes = PBcTypes
equationP._bc = {"inlet": valuePInlet, "outlet": valuePOutlet}
equationP._solver = 'sparseLU'

locS = pnm.LocS(netgrid, props, u, dS0, P, S0, data)
conS = pnm.ConS(netgrid, props, u, dS0, P, S0, data)
equationS = pnm.Equation(locS, conS, S0)
equationS._bcTypes = S0BcTypes
equationS._bc = {"inlet": valueS0Inlet, "outlet": valueS0Outlet}
equationS._solver = 'sparseLU'


def cfd_procedure():
    equationP.calculateCoefficients()
    equationP.fillMatrix()
    equationP.solve()

    equationS.calculateCoefficients()
    equationS.fillMatrix()
    equationS.solve()
    equationS.calculateData()

    equationP.iterateIndices()
    equationS.iterateIndices()


iPrev = equationP._iPrev
P_array_times = [np.copy(P[iPrev])]
P_av_times = [np.mean(P[iPrev])]
u_array_times = [np.copy(u)]
u_av_times = [np.mean(np.fabs(u))]
dS0_array_times = [np.copy(dS0)]
dS0_av_times = [np.mean(dS0)]
S0_array_times = [np.copy(S0[iPrev])]
S0_av_times = [np.mean(S0[iPrev])]
S0S1 = np.copy(S0[iPrev]) * (1. - np.copy(S0[iPrev]))
S0S1_array_times = [np.copy(S0S1)]
S0S1_av_times = [np.mean(S0S1)]

# ==========================================================================================
#                                     time loop                                            #
# ==========================================================================================

if t_step_out_coarse < dt:
    t_step_out_coarse = dt

if t_step_out_fine < dt:
    t_step_out_fine = dt

t_stepsN = round(t_max / dt)
step_out_coarse = round(t_step_out_coarse / dt)
print(step_out_coarse)
step_out_fine = round(t_step_out_fine / dt)
print(step_out_fine)
ts_out_coarse = [0]
ts_out_fine = [0]
t_cur = 0
with alive_bar(t_stepsN, title='time loop') as bar:
    for i in range(t_stepsN):
        t_cur += dt

        # pnm
        equationP._bc = {"inlet": smooth_bc(valuePInletIni, valuePInlet, t_bc, t_cur),
                         "outlet": smooth_bc(valuePOutletIni, valuePOutlet, t_bc, t_cur)}
        cfd_procedure()
        iPrev = equationP._iPrev
        S0S1 = np.copy(S0[iPrev]) * (1. - np.copy(S0[iPrev]))
        iPrev = equationP._iPrev

        # matrix
        matrix.imposeBc()
        matrix.solve()
        matrix.calculateData()
        matrix.iterateIndices()
        iPrevM = matrix._iPrev

        if (i + 1) % step_out_coarse == 0:
            ts_out_coarse.append(t_cur)
            # pnm
            P_array_times.append(np.copy(P[iPrev]))
            S0_array_times.append(np.copy(S0[iPrev]))
            u_array_times.append(np.copy(u))
            dS0_array_times.append(np.copy(dS0))
            S0S1_array_times.append(np.copy(S0S1))
            # matrix
            uM_array_times.append(np.copy(uM[iPrevM]))
            P0M_array_times.append(np.copy(P0M[iPrevM]))
            S0M_array_times.append(np.copy(S0M[iPrevM]))
            alpha1M_array_times.append(np.copy(alpha1M[iPrevM]))

        if (i + 1) % step_out_fine == 0:
            ts_out_fine.append(t_cur)
            # pnm
            P_av_times.append(np.mean(P[iPrev]))
            u_av_times.append(np.mean(np.fabs(u)))
            dS0_av_times.append(np.mean(dS0))
            S0_av_times.append(np.mean(S0[iPrev]))
            S0S1_av_times.append(np.mean(S0S1))
            # matrix
            uM_av_times.append(np.mean(uM[iPrevM]))
            P0M_av_times.append(np.mean(P0M[iPrevM]))
            S0M_av_times.append(np.mean(S0M[iPrevM]))
            alpha1M_av_times.append(np.mean(alpha1M[iPrevM]))

        bar()

end = time.time()
total_time = end - start
print("\n"+ str(round(total_time, 2)) +' s')
# ==========================================================================================
#                                     pnm VTK output                                       #
# ==========================================================================================


cells_file_name = 'pnm_real.pvd'
cells_files_names = list()
cells_files_descriptions = list()
for i in range(len(ts_out_coarse)):
    time = ts_out_coarse[i]
    P = P_array_times[i]
    S0 = S0_array_times[i]
    S0S1 = S0S1_array_times[i]
    u = u_array_times[i]
    dS0 = dS0_array_times[i]
    netgrid.cellsArrays = {'P': P, 'S0': S0, 'S0S1': S0S1}
    netgrid.facesArrays = {'u': u, 'dS0': dS0}
    cells_files_names.append(str(round(time, round_output_time)) + '.vtu')
    cells_files_descriptions.append(str(round(time, round_output_time)))
    netgrid.saveCells(cells_files_names[-1])
    saveFilesCollectionToFile(cells_file_name, cells_files_names, cells_files_descriptions)

# ==========================================================================================
#                                     matrix VTK output                                    #
# ==========================================================================================


cells_file_name = 'pnm_real_M.pvd'
cells_files_names = list()
cells_files_descriptions = list()
for i in range(len(ts_out_coarse)):
    time = ts_out_coarse[i]
    netgridM.cellsArrays = {'P0M': P0M_array_times[i],
                            'S0M': S0M_array_times[i],
                            'alpha1M': alpha1M_array_times[i]}
    netgridM.facesArrays = {'uM': uM_array_times[i]}
    cells_files_names.append(str(round(time, round_output_time)) + '_M.vtu')
    cells_files_descriptions.append(str(round(time, round_output_time)))
    netgridM.saveCells(cells_files_names[-1])
    saveFilesCollectionToFile(cells_file_name, cells_files_names, cells_files_descriptions)

# ==========================================================================================
#                                       csv output                                         #
# ==========================================================================================
d = {'ts_out_fine': ts_out_fine, 'S0M_av_times': S0M_av_times, 'alpha1M_av_times': alpha1M_av_times,
'P0M_av_times': P0M_av_times}
df = pd.DataFrame(data=d)
df.to_csv('data_pnm.csv')
# ==========================================================================================
#                             pnm and matrix matplotlib output                             #
# ==========================================================================================
fig_width = 12
fig, axs = plt.subplots(4, 2, figsize=(fig_width, fig_width))

# 1th set
ax0 = axs[0, 0]
ax0.plot(ts_out_fine, S0_av_times, label='$S_{\\mathit{0}av}$', marker='o', markersize=2, color='tab:green')
ax0.set_ylabel('$S_{\\mathit{0}}$')
ax0.set_xlabel('$t, s$')
ax0.legend(loc=0)

# 2th set
ax1 = axs[1, 0]
ax1.plot(ts_out_fine, dS0_av_times, label='$dS_{\\mathit{0}av}$', marker='o', markersize=2, color='tab:blue')
ax1.set_ylabel('$dS_{\\mathit{0}av}$')
ax1.set_xlabel('$t, s$')
ax1.legend(loc=0)

# 3th set
ax2 = axs[2, 0]
ax2.plot(ts_out_fine, P_av_times, label='$P_{av}$', marker='o', markersize=2, color='tab:cyan')
ax2.set_ylabel('$P, Pa$')
ax2.set_xlabel('$t, s$')
ax2.legend(loc=0)

# 4th set
ax3 = axs[3, 0]
ax3.plot(ts_out_fine, u_av_times, label='$u_{av}$', marker='o', markersize=2, color='tab:orange')
ax3.set_ylabel('$u, m/s$')
ax3.set_xlabel('$t, s$')
ax3.legend(loc=0)

# 5th set
ax4 = axs[0, 1]
ax4.plot(ts_out_fine, S0M_av_times, label='$S_{\\mathit{0}avM}$', marker='o', markersize=2, color='tab:green')
ax4.set_ylabel('$S_{\\mathit{0}}$')
ax4.set_xlabel('$t, s$')
ax4.legend(loc=0)

# 6th set
ax5 = axs[1, 1]
ax5.plot(ts_out_fine, alpha1M_av_times, label='$\\alpha_{\\mathit{1}avM}$', marker='o', markersize=2, color='tab:red')
ax5.set_ylabel('$\\alpha_{\\mathit{1}}$')
ax5.set_xlabel('$t, s$')
ax5.legend(loc=0)

# 7th set
ax6 = axs[2, 1]
ax6.plot(ts_out_fine, P0M_av_times, label='$P_{\\mathit{0}avM}$', marker='o', markersize=2, color='tab:cyan')
ax6.set_ylabel('$P_{\\mathit{0}}, Pa$')
ax6.set_xlabel('$t, s$')
ax6.legend(loc=0)

# 8th set
ax7 = axs[3, 1]
ax7.plot(ts_out_fine, uM_av_times, label='$u_{avM}$', marker='o', markersize=2, color='tab:orange')
ax7.set_ylabel('$u, m/s$')
ax7.set_xlabel('$t, s$')
ax7.legend(loc=0)

plt.tight_layout()
plt.show()
