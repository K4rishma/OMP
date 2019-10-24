import sys
import os
import matplotlib.animation as manimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import smoothing
import makeVectorStream as vstream
import PlotFuncs
import PlotOpts as opts
from tqdm import trange


def saveQuiverPlots(filename, ftitle=None):
    # determine save directory
    savefilename = filename + '_vec.mp4'
    print(savefilename)
    savedir, savename = os.path.split(savefilename)
    os.makedirs(savedir, exist_ok=True)
    if ftitle is None:
        ftitle = savename
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    matfile = PlotFuncs.loadmat(filename + '_piv.mat')  # load matfile containing synthetic data

    u = np.transpose(np.real(matfile['ufilt']), (1, 0, 2)).astype(np.float32)
    v = np.transpose(np.real(matfile['vfilt']), (1, 0, 2)).astype(np.float32)
    x = matfile['xfilt'].T.astype(np.float32)
    y = matfile['zfilt'].T.astype(np.float32)
    vmask = np.transpose(matfile['vmaskfilt'], (1, 0, 2))
    vmask[np.isnan(vmask)] = 0
    vmask = vmask.astype(bool)
    wcorr = np.transpose(np.real(matfile['Wcorrfilt']), (1, 0, 2)).astype(np.float32)

    binfile = matfile['bmode_path']
    datatype = matfile['BmodeType']
    datashape = matfile['BmodeShape']
    mask = matfile['mask']
    mask[np.isnan(mask)] = 0
    mpoly = matfile['mPoly']
    roi = matfile['roi'].ravel().astype(int)
    tRes = float(matfile['tRes'])
    nAve = float(matfile['n_ave'])
    xRes = float(matfile['xRes'])
    yRes = float(matfile['zRes'])
    first = int(matfile['piv_start'])
    last = int(matfile['piv_stop'])
    origin = matfile['origin']
    unit = matfile['unit']
    if unit == 'm/s':
        unit_factor = 1000.0
    elif unit == 'mm/s':
        unit_factor = 1.0
    elif unit == 'cm/s':
        unit_factor = 10.0
    else:
        exit('unit not recognised. Only m/s, cm/s and mm/s are currently recognised')
    del matfile
    roi = np.array([roi[0], roi[1], roi[0]+roi[2], roi[1]+roi[3]])  # [xmin, ymin. xmax, ymax]
    bmode = np.fromfile(binfile, datatype)
    bmode = bmode.reshape(datashape[::-1])
    bmode = np.transpose(bmode, (2, 1, 0))  # [t, y, x] -> [y, x, t]

    # smoothing
    sm_str = opts.smooth
    if opts.smooth == 'fourier':
        u[~vmask] = 0
        v[~vmask] = 0
        u, v = smoothing.Fourier_smooth(u, v, prop=(0.75, 0.75, 0.30))
        u[~vmask] = np.nan
        v[~vmask] = np.nan

    if opts.smooth == 'gauss':
        u[~vmask] = 0
        v[~vmask] = 0
        u, v = smoothing.gauss_smooth(u, v, std=[0.5, 0.5, 0.1], trunc=[1.0, 1.0, 1.0])
        u[~vmask] = np.nan
        v[~vmask] = np.nan

    if opts.smooth == 'svd':
        ni, nj, ns = u.shape
        u[np.isnan(u)] = 0
        v[np.isnan(v)] = 0
        S = np.concatenate((u.reshape((ni*nj, ns), order='F'), v.reshape((ni*nj, ns), order='F')), axis=0)
        U, s, V = sp.linalg.svd(S, lapack_driver='gesvd')
        s[40:-1] = 0
        # s[2:] = 0
        Sigma = sp.linalg.diagsvd(s, U.shape[0], V.shape[0])
        S1 = np.dot(U, np.dot(Sigma, V))
        S1 = S1.reshape((ni, nj, 2, ns), order='F')
        u = np.squeeze(S1[:, :, 0, :])
        v = np.squeeze(S1[:, :, 1, :])

    if opts.temp_filt:
        sm_str += '_tm'
        u, v = smoothing.temp_mov_ave_gaussian(u, v, sigma=2.0, trunc=3.0, axes=2, plotwindow=True)

    imsize = origin

    # cleanup and calibrate vector maps
    u[np.isnan(u)] = 0
    v[np.isnan(v)] = 0

    # find magnitudes for colormaps
    velmag = np.sqrt(u ** 2 + v ** 2)
    # velmag_signed = velmag*np.sign(v)
    velmax = np.nanmax(velmag)  # 99th percentile velocity magnitude in the data set
    # velmin = np.nanmin(velmag)
    # velmax_signed = np.nanmax(velmag_signed)
    # velmin_signed = np.nanmin(velmag_signed)
    # velprc99 = np.nanpercentile(velmag, 99)
    # velmean = np.nanmean(velmag)
    print('Max velocity in field: {}'.format(velmax))
    maxind = np.asarray(np.unravel_index(np.nanargmax(velmag), velmag.shape))

    velprc95 = np.nanpercentile(velmag, 95)
    scl = (velprc95/np.sqrt(np.diff(x[:, 0]).mean()**2 + np.diff(y[0, :]).mean()**2))/opts.sf
    u_vis = u  # mag_norm*np.cos(theta)
    v_vis = v  # mag_norm*np.sin(theta)
    masked = np.where(~vmask)
    u_vis[masked] = np.nan
    v_vis[masked] = np.nan
    velmag_vis = velmag
    velmag_vis[masked] = np.nan

    ## Make movie of images
    # logo = image.imread('K:/emclogo.png')
    # logo_asp = logo.shape[0] / logo.shape[1]
    # logo_width = 30
    if opts.cmapmax == 'prc99':
        vcmax = np.nanpercentile(velmag, 99)
    elif isinstance(opts.cmapmax, float):
        vcmax = opts.cmapmax
    else:
        vcmax = velmax


    imsize_rev = [imsize[0], imsize[1], imsize[3], imsize[2]]
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='Movie Test', artist='',
                    comment='')
    writer = FFMpegWriter(fps=30, codec='h264', metadata=metadata)
    figsize_inches = (5, 4)
    fig, ax = plt.subplots(figsize=figsize_inches)
    i = maxind[2]
    bim = ax.imshow(bmode[..., i], cmap='gray', vmin=opts.drange[0], vmax=opts.drange[1], aspect=1.0,
                    interpolation='hamming', extent=imsize, origin='lower', alpha=0.8)
    quiv = ax.quiver(x, y, u_vis[..., i], v_vis[..., i], velmag_vis[..., i],
                     cmap=opts.vcmap, angles='xy', scale_units='xy', scale=scl, units='xy', width=0.4,
                     clim=[0, vcmax])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "5%", pad="3%")
    cbar = fig.colorbar(mappable=quiv, cax=cax)
    cbar.set_label('Vector Magnitude [m/s]')
    # plt.setp(ax.get_xticklabels(), rotation=30, fontsize=10)
    ax.set_xlabel('Azimuth [mm]')
    ax.set_ylabel('Depth [mm]')
    ax.set_facecolor((1, 1, 1))
    ax.set_xlim([imsize_rev[0], imsize_rev[1]])
    ax.set_ylim([imsize_rev[2], imsize_rev[3]])
    ax.set_title("{:5.1f}".format((0 + first) * tRes * 1000) + ' ms')
    # ax.imshow(logo, extent=(33, 33 + logo_width, 14 + logo_asp * logo_width, 14), aspect='equal')
    # plt.ioff()
    fig.tight_layout(pad=0.1)
    with writer.saving(fig, savefilename, dpi=opts.dpi):
        for i in trange(u.shape[2]):
            try:
                bim.set_data(bmode[..., i])
                quiv.set_UVC(u_vis[..., i], v_vis[..., i], velmag_vis[..., i])
                ax.set_title(ftitle + f'-> {(i + first) * tRes * 1000:5.1f} ms ')
                writer.grab_frame()
            except KeyboardInterrupt:
                break


def test_command_line(arg1, arg2):
    print(arg1)
    print(arg2)
    print('Done')


if __name__ == '__main__':
    saveQuiverPlots(sys.argv[1], sys.argv[2])
    # test_command_line(sys.argv[1], sys.argv[2])
