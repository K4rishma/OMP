import numpy as np
from matplotlib import pyplot as plt
from tqdm import trange


def vectors2d_over_bmode_pdf(x, y, u, v, bmode, savefile, bmode_drange=(-40, 0), vrange=None, imsize=None, tlims=None,
                             scl=0.1, vcmap='inferno'): ## show vector maps
    import matplotlib.animation as manimation
    from matplotlib.backends.backend_pdf import PdfPages
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    mag = np.sqrt(u ** 2 + v ** 2)
    if tlims is None:
        tlims = [0, u.shape[-1]]
    if vrange is None:
        vrange = [np.nanmin(mag[..., tlims[0]:tlims[1]]), np.nanmax(mag[..., tlims[0]:tlims[1]])]
    if imsize is None:
        imsize = [x.min(), x.max(), y.min(), y.max()]
    # Plot vector maps (multipage pdf)
    fig, axs = plt.subplots()
    piv = axs.quiver(x[..., 0], y[..., 0], u[..., 0], v[..., 0], mag[..., 0],
                     cmap=vcmap, angles='xy', scale_units='xy', units='xy', width=0.25,
                     clim=vrange, scale=scl)
    if bmode is not None:
        bg = axs.imshow(bmode[..., 0], cmap='gray', vmin=bmode_drange[0], vmax=bmode_drange[1], aspect=1.0,
                        interpolation='hamming', extent=imsize)
    divider = make_axes_locatable(axs)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(piv, cax=cax)
    cbar.set_label('Vector Magnitude [m/s]')
    plt.setp(axs.get_xticklabels(), rotation=30, fontsize=10)
    axs.set_xlim([imsize[0], imsize[1]])
    axs.set_ylim([imsize[2], imsize[3]])
    axs.set_xlabel('Azimuth [mm]')
    axs.set_ylabel('Depth [mm]')
    # plt.tight_layout()
    #
    with PdfPages(savefile) as pdf:
        for i in trange(tlims[0], tlims[1]):
            try:
                piv.set_UVC(u[..., i], v[..., i], mag[..., i])
                if x.shape == u.shape:
                    piv.set_offsets(np.vstack((x[..., i], y[..., i])).T)
                axs.set_title(str(i) + 'ms')
                if bmode is not None:
                    bg.set_data(bmode[..., i])
                plt.draw()
                pdf.savefig()
                # time.sleep(0.05)
            except KeyboardInterrupt:
                break


def stream_lines_over_bmode_pdf(x, y, u, v, bmode, savefile, bmode_drange=(-40, 0), vrange=None, imsize=None,
                                tlims=None, lw=1.0, vcmap='jet'): ## show vector maps
    import matplotlib.animation as manimation
    from matplotlib.backends.backend_pdf import PdfPages
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.colors import Normalize

    mag = np.sqrt(u**2 + v**2)
    if tlims is None:
        tlims = [0, u.shape[-1]]
    if vrange is None:
        vrange = [np.nanmin(mag[..., tlims[0]:tlims[1]]), np.nanmax(mag[..., tlims[0]:tlims[1]])]
    if imsize is None:
        imsize = [x.min(), x.max(), y.min(), y.max()]

    fig, ax = plt.subplots()
    plt.setp(ax.get_xticklabels(), rotation=30, fontsize=10)
    plt.xlabel('Azimuth [mm]')
    plt.ylabel('Depth [mm]')
    ax.set_xlim([imsize[0], imsize[1]])
    ax.set_ylim([imsize[2], imsize[3]])
    # lw = 5*velmag/velmag.max()
    bim = ax.imshow(bmode[..., 0], cmap='gray', vmin=bmode_drange[0], vmax=bmode_drange[1], aspect=1.0,
                    interpolation='hamming', extent=imsize)
    stream = ax.streamplot(x[:, 0], y[0, :], u[..., 0].T, v[..., 0].T, density=1, color=mag[..., 0].T,
                           cmap=vcmap, linewidth=lw, norm=Normalize(vmin=vrange[0], vmax=vrange[1], clip=False))
    cbar = plt.colorbar(stream.lines, orientation='vertical')
    cbar.set_label('Vector Magnitude [m/s]')
    with PdfPages(savefile) as pdf:
        for i in trange(tlims[0], tlims[1]):
            try:
                # tic = time()
                bim.set_data(bmode[..., i])
                stream.lines.remove()
                ax.patches = []
                stream = ax.streamplot(x[:, 0], y[0, :], u[..., i].T, v[..., i].T, density=1,
                                       color=mag[..., i].T, cmap=vcmap, linewidth=lw,
                                       norm=Normalize(vmin=vrange[0], vmax=vrange[1], clip=True))
                ax.set_title(str(i) + ' ms')
                plt.draw()
                pdf.savefig()
                # print(time() - tic)
            except KeyboardInterrupt:
                break