# configuration.py

import numpy as np
import numerov_qex as nv

"""
information of the coupled-channel system:
we can just add channels
for gaussian+well in 0->0 potential we can do stuff like...
nv.PotentialInfo("gaussian", [-1.0, 1.0], 0, 0),
nv.PotentialInfo("well",     [-3.0, 2.0], 0, 0),
"""

# coupled-channels config
# nv.ChannelInfo(channel ID, name, m1, m2, ell, z1, z2)
system = nv.ChannelSystem([
    nv.ChannelInfo(0, "channel-0 (l=0)", 1.0, 1.0, 0, +1.0, -1.0),
    nv.ChannelInfo(1, "channel-1 (l=0)", 1.0, 1.5, 0, +1.0, -1.0),
    #nv.ChannelInfo(2, "channel-2 (l=0)", 1.0, 1.7, 0, +0.0, +0.0),
])

# potentials
# nv.PotentialInfo(name, pars, channel ID, channel ID)
# support: "gaussian", "well"
# append -> auto sum
potential = nv.PotentialMatrix([
    nv.PotentialInfo("gaussian", [-5.0, 1.0], 0, 0),
    nv.PotentialInfo("gaussian", [-4.0, 1.0], 1, 1),
    nv.PotentialInfo("gaussian", [-1.0, 1.0], 0, 1),
    nv.PotentialInfo("gaussian", [-1.0, 1.0], 1, 0),
], system)