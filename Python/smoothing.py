import numpy as np

def POD(U):
    from scipy.linalg import eigh
    C = np.dot(U.conj().T, U)  # create covariance matrix
    lmb, V = eigh(C)  # eigenvalue decomp
    lmb[lmb <= 0] = np.spacing(1)
    lmb = np.sqrt(lmb[::-1])  # sort eigvals descending, take sqrt
    V = V[:, ::-1]  # sort eigvecs descending
    Phi = np.dot(U, V)  # project original dataset onto eig vecs
    Phi /= lmb  # normalized Phi
    A = np.dot(Phi.conj().T,  U)  # determine coefficients
    return [Phi, A, lmb]

def POD_smooth(u, v, ens=15, trunc=-1):
    from scipy.ndimage.filters import uniform_filter1d
    dim0 = u.shape
    U = np.concatenate((u.reshape((-1, u.shape[-1])), v.reshape((-1, v.shape[-1]))), axis=0)
    nanidx = np.where(np.isnan(U))
    U[nanidx] = 0
    Phi, A, lmb = POD(U)
    A = uniform_filter1d(A, size=ens, axis=1)
    A[trunc:] = 0
    U = np.dot(Phi, A)
    U[nanidx] = np.nan
    u1 = U[:np.prod(dim0[:-1])].reshape(dim0)
    v1 = U[np.prod(dim0[:-1]):].reshape(dim0)
    return [u1, v1]

def DMD(X, r=-1, dt=1):
    from scipy.linalg import svd, eig, inv
    if r == -1:
        r = X.shape[1]-1
    X1 = X[:, :-1]
    X2 = X[:, 1:]

    U, s, V = svd(X1, full_matrices=False, lapack_driver='gesvd')
    V = V.T
    r = np.min((r, U.shape[1]))

    U_r = U[:, :r]
    S_r = np.diag(s[:r])
    V_r = V[:, :r]
    XVE = np.dot(np.dot(X2, V_r), inv(S_r))
    A_tilde = np.dot(U_r.conj().T, XVE)
    lmb, W_r = eig(A_tilde)
    # Phi = np.dot(XVE, W_r)
    Phi = np.dot(U_r, W_r)

    omega = np.log(lmb)/dt

    b = sp.linalg.lstsq(Phi, X1[:, 0])[0]

    mm1 = X1.shape[1]
    t = np.r_[0:mm1+1]*dt
    td = b[:, np.newaxis]*np.exp(np.outer(omega, t))
    Xdmd = np.real(np.dot(Phi, td))
    return [Phi, omega, lmb, b, td, Xdmd]

def DMD_smooth(u, v, trunc=-1):
    dim0 = u.shape
    X = np.concatenate((u.reshape((-1, u.shape[-1])), v.reshape((-1, v.shape[-1]))), axis=0).astype(np.double)
    nanidx = np.where(np.isnan(X))
    X[nanidx] = 0
    Phi, omega, lmb, b, td, X = DMD(X, trunc, 1)
    X[nanidx] = np.nan
    u1 = X[:np.prod(dim0[:-1])].reshape(dim0)
    v1 = X[np.prod(dim0[:-1]):].reshape(dim0)
    # Phi_u = Phi[:np.prod(dim0[:-1])].reshape((dim0[0], dim0[1], -1))
    # Phi_v = Phi[np.prod(dim0[:-1]):].reshape((dim0[0], dim0[1], -1))
    return [u1, v1]

def Fourier_smooth(u, v, prop=(0.5, 0.5, 0.2), offset=(0, 0, 0)):
    from scipy.signal import tukey
    dim0 = u.shape
    U = u + 1j*v
    nanidx = np.where(np.isnan(U))
    U[nanidx] = 0
    Ufft = np.fft.fftshift(np.fft.fftn(U, axes=(0, 1, 2)), axes=(0, 1, 2))
    wsize = np.rint(np.asarray(dim0)*np.asarray(prop)).astype(int)
    winx = tukey(wsize[0], 0.2, sym=True)
    winy = tukey(wsize[1], 0.2, sym=True)
    winz = tukey(wsize[2], 0.2, sym=True)
    win = winx[:, np.newaxis, np.newaxis]*winy[:, np.newaxis]*winz
    mctr = np.rint(np.asarray(dim0)/2).astype(int)
    maskarr = np.zeros(dim0)
    id0 = np.rint(mctr+offset-wsize/2).astype(int)
    id1 = id0 + win.shape
    maskarr[id0[0]:id1[0], id0[1]:id1[1], id0[2]:id1[2]] = win
    U = np.fft.ifftn(np.fft.ifftshift(Ufft*maskarr, axes=(0, 1, 2)), axes=(0, 1, 2))
    u1 = U.real
    v1 = U.imag
    return [u1, v1]

def gauss_smooth(u, v, std=[0.5, 0.5, 0.1], trunc=[3.0, 3.0, 1.0]):
    from scipy.ndimage.filters import gaussian_filter1d
    u = gaussian_filter1d(u, std[0], axis=0, order=0, mode='nearest', truncate=trunc[0])  # x ax
    u = gaussian_filter1d(u, std[1], axis=1, order=0, mode='nearest', truncate=trunc[1])  # y ax
    u = gaussian_filter1d(u, std[2], axis=2, order=0, mode='nearest', truncate=trunc[2])  # t ax
    v = gaussian_filter1d(v, std[0], axis=0, order=0, mode='nearest', truncate=trunc[0])  # x ax
    v = gaussian_filter1d(v, std[1], axis=1, order=0, mode='nearest', truncate=trunc[1])  # y ax
    v = gaussian_filter1d(v, std[2], axis=2, order=0, mode='nearest', truncate=trunc[2])  # t ax
    return (u, v)

def temp_mov_ave(u, v, ens, axes=None):
    from scipy.ndimage.filters import uniform_filter1d
    if axes is None:
        axes = u.ndim-1
    u = uniform_filter1d(u, ens, axis=axes)
    v = uniform_filter1d(v, ens, axis=axes)
    return (u, v)

def temp_mov_ave_gaussian(u, v, sigma=0.6, trunc=4.0, axes=None, plotwindow=False):
    from scipy.ndimage.filters import gaussian_filter1d
    if axes is None:
        axes = u.ndim-1
    radius = int(trunc * sigma + 0.5)
    print(f'gaussian window length = {2*radius} frames')
    if plotwindow:
        import matplotlib.pyplot as plt
        p = np.polynomial.Polynomial([0, 0, -0.5 / (sigma * sigma)])
        x = np.arange(-radius, radius + 1)
        phi_x = np.exp(p(x), dtype=np.double)
        phi_x /= phi_x.sum()
        plt.plot(x, phi_x)
        plt.xlabel('Pixels')
        plt.ylabel('Normalized Magnitude')
        plt.title('Convolution Kernel')
        plt.ylim([0, phi_x.max()])
    u = gaussian_filter1d(u, sigma=sigma, truncate=trunc, axis=axes, mode='nearest')
    v = gaussian_filter1d(v, sigma=sigma, truncate=trunc, axis=axes, mode='nearest')

    return (u, v)