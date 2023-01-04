import logging
from glycogen_module.core import GlycogenStructure
import numpy as np

from pylab import pi, plt
logger = logging.getLogger(__name__)


def plot_structure(gs: GlycogenStructure, colors: list[str] = ['black', 'limegreen'], figsize: tuple[int] = (14, 10), size: int = 30):
    plot_nre_in_different_color = True if len(colors) > 1 else False
    markeredgecolor = colors[0]
    if plot_nre_in_different_color:
        nre_color = colors[1]

    XLIST = []
    YLIST = []
    ZLIST = []

    XNRE = []
    YNRE = []
    ZNRE = []

    for chain in gs.chains:
        XNRE.append(GlycogenStructure.L*chain.glucose_positions[-1][0])
        YNRE.append(GlycogenStructure.L*chain.glucose_positions[-1][1])
        ZNRE.append(GlycogenStructure.L*chain.glucose_positions[-1][2])
        for pos in chain.glucose_positions:
            #logger.debug(f"chain {chain.id}, glucose_positions: {chain.glucose_positions}.Pos: {pos}")
            XLIST.append(pos[0])
            YLIST.append(pos[1])
            ZLIST.append(pos[2])

    XLIST, YLIST, ZLIST = GlycogenStructure.L * \
        np.asarray(XLIST),  GlycogenStructure.L * \
        np.asarray(YLIST), GlycogenStructure.L * np.asarray(ZLIST)

    ax = plt.figure(figsize=figsize)
    ax = ax.add_subplot(projection='3d')

    ax.plot(XLIST, YLIST, ZLIST, 'o', markersize=7, markeredgewidth=0.2,
            markeredgecolor=markeredgecolor, label='glycogen 3d')

    if plot_nre_in_different_color:
        ax.plot(XNRE, YNRE, ZNRE, 'o', markersize=7, markeredgewidth=0.2,
                color=nre_color, markeredgecolor=markeredgecolor, label='glycogen 3d')

    ax.view_init(elev=15, azim=60)
    plt.xlabel('x [nm]', fontsize='15')

    plt.ylabel('y [nm]', fontsize='15')
    ax.set(xlim=(-size, size), ylim=(-size, size), zlim=(-size, size))

    ax.set_zlabel(r'z [nm]', fontsize='25', rotation=200, labelpad=20)

    plt.show()
    ax.legend()
