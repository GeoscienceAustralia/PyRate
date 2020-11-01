import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

PI = np.pi

cmap = mpl.cm.Spectral
# norm = mpl.colors.Normalize(vmin=-PI/2, vmax=PI/2)


def plot_closure(closure: np.ndarray):
    nrows, ncols, n_loops = closure.shape
    fig = plt.figure()

    plt_rows = np.int(np.sqrt(n_loops))
    plt_cols = n_loops//plt_rows
    if n_loops % plt_rows:
        plt_cols += 1

    tot_plots = 1
    for p_r in range(plt_rows):
        for p_c in range(plt_cols):
            ax = fig.add_subplot(plt_rows, plt_cols, tot_plots)
            data = closure[:, :, p_r * (p_c-1) + p_c]
            im = ax.imshow(data, vmin=-PI/2, vmax=PI/2, cmap=cmap)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(im, cax=cax)
            if tot_plots == n_loops:
                break
            tot_plots += 1

    # ax = fig.add_subplot(plt_rows, plt_cols, tot_plots+1)
    # fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap), cax=ax, orientation='horizontal', label='radians')

    plt.show()