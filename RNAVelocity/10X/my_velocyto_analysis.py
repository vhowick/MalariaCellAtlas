import sys
import numpy as np
import pandas as pd 
import matplotlib
import matplotlib.pyplot as plt
import loompy
import velocyto as vcy
import logging
from sklearn.svm import SVR
from sklearn.linear_model import LinearRegression
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.interpolate import interp1d

# plotting utility functions
def despline():
    ax1 = plt.gca()
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')


def minimal_xticks(start, end):
    end_ = np.around(end, -int(np.log10(end))+1)
    xlims = np.linspace(start, end_, 5)
    xlims_tx = [""]*len(xlims)
    xlims_tx[0], xlims_tx[-1] = f"{xlims[0]:.0f}", f"{xlims[-1]:.02f}"
    plt.xticks(xlims, xlims_tx)


def minimal_yticks(start, end):
    end_ = np.around(end, -int(np.log10(end))+1)
    ylims = np.linspace(start, end_, 5)
    ylims_tx = [""]*len(ylims)
    ylims_tx[0], ylims_tx[-1] = f"{ylims[0]:.0f}", f"{ylims[-1]:.02f}"
    plt.yticks(ylims, ylims_tx)


# save raw output 
vlm = vcy.VelocytoLoom("/lustre/scratch117/cellgen/team218/TA/Malaria_RNAVelocity/10XOutput/velocyto_output/possorted_genome_bam_UVJ4K.loom")

np.savetxt("out_S_matrix.csv", vlm.S, delimiter=",")
np.savetxt("out_U_matrix.csv", vlm.U, delimiter=",")
np.savetxt("out_A_matrix.csv", vlm.A, delimiter=",")

df = pd.DataFrame(vlm.ca)
df.to_csv("out_ca_matrix.csv")

df = pd.DataFrame(vlm.ra)
df.to_csv("out_ra_matrix.csv")

# filtering & normalization
vlm = vcy.VelocytoLoom("/lustre/scratch117/cellgen/team218/TA/Malaria_RNAVelocity/10XOutput/velocyto_output/possorted_genome_bam_UVJ4K.loom")

vlm.filter_cells(bool_array=vlm.initial_Ucell_size > np.percentile(vlm.initial_Ucell_size, 0.4))
vlm.score_detection_levels(min_expr_counts=40, min_cells_express=30)
vlm.filter_genes(by_detection_levels=True)
vlm.score_detection_levels(min_expr_counts=0, min_cells_express=0, min_expr_counts_U=25, min_cells_express_U=20)
vlm.filter_genes(by_detection_levels=True)
# best with sample and expression scaling
vlm._normalize_S(relative_size=vlm.initial_cell_size,
                 target_size=np.mean(vlm.initial_cell_size))
vlm._normalize_U(relative_size=vlm.initial_Ucell_size,
                 target_size=np.mean(vlm.initial_Ucell_size))
# PCA
vlm.perform_PCA()
plt.plot(np.cumsum(vlm.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(vlm.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
n_comps

#knn smoothing
k = 500
vlm.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8, b_maxl=k*4, n_jobs=16)

# Fit model
vlm.fit_gammas(limit_gamma=False, fit_offset=False)
vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift(assumption="constant_velocity")
vlm.extrapolate_cell_at_t(delta_t=1.)

vlm.estimate_transition_prob(hidim="Sx_sz", embed="pca", transform="sqrt", psc=1, n_neighbors=2000, knn_random=True, sampled_fraction=0.5)
vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=False)
vlm.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=300)

plt.figure(None,(14,14))
quiver_scale = 60

plt.scatter(vlm.embedding.components_[:, 0], vlm.embedding.components_[:, 1], c="0.8", alpha=0.2, s=10, edgecolor="")





