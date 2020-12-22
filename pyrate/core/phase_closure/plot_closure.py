from typing import List
import numpy as np
from pyrate.core.phase_closure.mst_closure import WeightedLoop

PI = np.pi
# norm = mpl.colors.Normalize(vmin=-PI/2, vmax=PI/2)


def plot_closure(closure: np.ndarray, loops: List[WeightedLoop]):
    try:
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        cmap = mpl.cm.Spectral
    except ImportError as e:
        raise ImportError(e)

    nrows, ncols, n_loops = closure.shape

    plt_rows = np.int(np.sqrt(n_loops))
    plt_cols = n_loops//plt_rows
    if n_loops % plt_rows:
        plt_cols += 1

    fig = plt.figure(figsize=(6*plt_rows, 4*plt_cols))

    tot_plots = 1
    for p_r in range(plt_rows):
        for p_c in range(plt_cols):
            ax = fig.add_subplot(plt_rows, plt_cols, tot_plots)
            data = closure[:, :, p_r * (p_c-1) + p_c]
            loop = loops[p_r * (p_c-1) + p_c]
            leg = ',\n'.join([swe.SignedEdge.edge.first.isoformat() + '-' + swe.SignedEdge.edge.second.isoformat() for swe in loop.loop])
            im = ax.imshow(data, vmin=-PI/2, vmax=PI/2, cmap=cmap)
            text = ax.text(10, 10, leg, bbox={'facecolor': 'white', 'pad': 5})
            text.set_fontsize(min(10, int(n_loops/10)))

            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(im, cax=cax)
            if tot_plots == n_loops:
                break
            tot_plots += 1

    # ax = fig.add_subplot(plt_rows, plt_cols, tot_plots+1)
    # fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap), cax=ax, orientation='horizontal', label='radians')

    plt.savefig(f'Closure-{len(loops)}.png')
