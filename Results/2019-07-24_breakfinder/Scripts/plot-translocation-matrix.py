'''
plot-translocation-matrix
==========

Plot the new contact matrix for the T2E translocation
'''

from __future__ import division, absolute_import, print_function
import os.path as path
import numpy as np
import cooler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ==============================================================================
# Functions
# ==============================================================================


# def scale_matrix(m, lower, upper):
#     '''
#     Scale a contact matrix to new lower and upper bounds

#     Parameters
#     ----------
#     m : numpy.array
#         Matrix to be scaled
#     lower : float
#         New lower bound
#     upper : float
#         New upper bound
#     '''
#     # lower and upper bounds of m
#     finite_vals = m[np.isfinite(np.log10(m))]
#     min_val = np.min(finite_vals)
#     max_val = np.max(finite_vals)
#     print(min_val, max_val)
#     # return m * upper / max_val
#     return lower + (m - min_val) / (max_val - min_val) * (upper - lower)


# ==============================================================================
# Main
# ==============================================================================
# contact matrix for PCa13848
c = cooler.Cooler(path.join("..", "..", "Data", "Processed", "2019-06-18_PCa-LowC-sequencing", "Contacts", "PCa13848.mcool::/resolutions/10000"))
loci = [
    "chr14:34000000-35719999",
    "chr21:38500000-41500000",
    "chr14:35720000-37000000",
]
m_00 = c.matrix().fetch(loci[0], loci[0])
m_01 = c.matrix().fetch(loci[0], loci[1])[:, ::-1] # flip columns corresponding to insertion orientation
m_02 = c.matrix().fetch(loci[0], loci[2])
m_11 = c.matrix().fetch(loci[1], loci[1])[::-1, ::-1]  # flip rows and columns corresponding to insertion orientation
m_12 = c.matrix().fetch(loci[1], loci[2])[::-1, :] # flip rows corresponding to insertion orientation
m_22 = c.matrix().fetch(loci[2], loci[2])

# combine into new block matrix
new_m = np.block([[m_00, m_01, m_02], [m_01.transpose(), m_11, m_12], [m_02.transpose(), m_12.transpose(), m_22]])

# # scale inter-chromosomal matrices so that they have the same range of values as the intra-chromosomal ones
# scale = [1, 0]
# for m in [m_00, m_02, m_11, m_22]:
#     finite_vals = m[np.isfinite(np.log10(m))]
#     min_val = np.min(finite_vals)
#     max_val = np.max(finite_vals)
#     if min_val < scale[0]:
#         scale[0] = min_val
#     if max_val > scale[1]:
#         scale[1] = max_val

# find the range of values of the inter-chromosomal matrices
# scaled_m_01 = scale_matrix(m_01, scale[0], scale[1])

# scaled_m = np.block([
#     [m_00, scaled_m_01, m_02],
#     [m_01.transpose(), m_11, m_12],
#     [m_02.transpose(), m_12.transpose(), m_22]
# ])

# plot as heatmap
plt.imshow(np.log10(new_m), cmap="YlOrRd", interpolation=None)
plt.axvline(x=m_00.shape[1])
plt.axvline(x=m_00.shape[1] + m_11.shape[1])
# plt.axhline(y=m_00.shape[0])
# plt.axhline(y=m_00.shape[0] + m_11.shape[0])
plt.savefig("Plots/translocation-matrix.png")
plt.close()

# # plot scaled matrix as heatmap
# plt.imshow(np.log10(scaled_m), cmap="YlOrRd", interpolation=None)
# plt.axvline(x=m_00.shape[1])
# plt.axvline(x=m_00.shape[1] + m_11.shape[1])
# # plt.axhline(y=m_00.shape[0])
# # plt.axhline(y=m_00.shape[0] + m_11.shape[0])
# plt.savefig("Plots/translocation-matrix.scaled.png")
# plt.close()
