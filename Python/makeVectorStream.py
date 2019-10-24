from tqdm import trange
import numpy as np

def in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    from scipy.spatial import Delaunay
    if not isinstance(hull, Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p) >= 0


def plot_in_hull(p, hull):
    """
    plot relative to `in_hull` for 2d data
    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection, LineCollection

    from scipy.spatial import Delaunay
    if not isinstance(hull, Delaunay):
        hull = Delaunay(hull)

    # plot triangulation
    poly = PolyCollection(hull.points[hull.vertices], facecolors='w', edgecolors='b')
    plt.clf()
    plt.title('in hull')
    plt.gca().add_collection(poly)
    plt.plot(hull.points[:, 0], hull.points[:, 1], 'o', hold=1)

    # plot the convex hull
    edges = set()
    edge_points = []

    def add_edge(i, j):
        """Add a line between the i-th and j-th points, if not in the list already"""
        if (i, j) in edges or (j, i) in edges:
            # already added
            return
        edges.add((i, j))
        edge_points.append(hull.points[[i, j]])

    for ia, ib in hull.convex_hull:
        add_edge(ia, ib)

    lines = LineCollection(edge_points, color='g')
    plt.gca().add_collection(lines)
    plt.show()

    # plot tested points `p` - black are inside hull, red outside
    inside = in_hull(p, hull)
    plt.plot(p[inside, 0], p[inside, 1], '.g')
    plt.plot(p[-inside, 0], p[-inside, 1], '.r')


def make_roi_from_mask(mask):
    temp = np.argwhere(mask)
    (y0, x0), (y1, x1) = temp.min(0), temp.max(0) + 1
    return np.asarray([y0, x0, y1, x1], dtype=int)

def convert_to_delauny(p):
    from scipy.spatial import Delaunay
    return Delaunay(p)

def makeMaskDelauny(mask, roi, origin, xRes, yRes):
    from skimage.feature import canny
    # convert mask to convex hull in spatial coordinates
    x_cntr = origin[0]
    y_cntr = origin[2]
    mask_edge = canny(mask.astype(np.float32), sigma=0.01)
    mask_y, mask_x = np.where(mask_edge)
    mask_edge = mask_edge[roi[1]:roi[3], roi[0]:roi[2]]
    mask_x, mask_y = mask_x.astype(np.float64), mask_y.astype(np.float64)
    mask_y = (mask_y + y_cntr) * yRes
    mask_x = (mask_x + x_cntr) * xRes
    imsize = [mask_x.min(), mask_x.max(), mask_y.max(), mask_y.min()]
    imsize_rev = [mask_x.max(), mask_x.min(), mask_y.min(), mask_y.max()]
    mask_pnts = np.vstack((mask_x*0.95, mask_y*0.95)).T
    mask_delauny = convert_to_delauny(mask_pnts)
    return (mask_delauny, imsize, mask_edge)

def makeMaskDelaunyScaled(mask, origin, xRes, yRes):
    from skimage.feature import canny
    # convert mask to convex hull in spatial coordinates
    x_cntr = origin[0]  # + roi[0]*xRes
    y_cntr = origin[2]  # + roi[1]*yRes
    mask_edge = canny(mask.astype(np.float32), sigma=0.01)
    mask_y, mask_x = np.where(mask_edge)
    mask_x, mask_y = mask_x.astype(np.float64), mask_y.astype(np.float64)
    mask_y = mask_y * yRes + y_cntr
    mask_x = mask_x * xRes + x_cntr
    imsize = [mask_x.min(), mask_x.max(), mask_y.max(), mask_y.min()]
    mask_pnts = np.vstack((mask_x, mask_y)).T
    mask_delauny = convert_to_delauny(mask_pnts)
    return (mask_delauny, mask_edge)


def vectorstream(x, y, u, v, mask_delauny, exts, tRes, randomize=True, move_along=True, vel_weighting=True, seed_cnt=1000,
                 rnd_val=0.2, unit_factor=1000.0):
    from scipy.interpolate import RectBivariateSpline
    # generate vector-streams
    #unit_factor = 1000.0  # for converting between mm and m. If u/v in mm/s then 1, if m/s then 1000.
    seed_x = np.empty((seed_cnt, u.shape[2] + 1))
    seed_x.fill(np.nan)
    seed_y = seed_x.copy()
    seed_u = np.zeros_like(seed_x)
    seed_v = seed_u.copy()
    seed_mag = seed_u.copy()
    vel_weight = [0.0, 0.0]  # weighting factor for seed placement based on velocity


    #place seeds randomly
    seed_x[..., 0] = np.random.random(seed_cnt) * (exts[1] - exts[0]) + exts[0]
    seed_y[..., 0] = np.random.random(seed_cnt) * (exts[3] - exts[2]) + exts[2]
    seed_x[..., 0], seed_y[..., 0] = ensure_in_poly(seed_x[..., 0], seed_y[..., 0], exts, mask_delauny,
                                                    weighting=vel_weight)

    for i in trange(u.shape[2]):
        if randomize:  # and (i % 10 == 0):
            rnd_inds = np.random.randint(0, seed_cnt, int(seed_cnt * rnd_val)).astype(np.uint)
            seed_x[rnd_inds, i] = np.random.random(int(seed_cnt * rnd_val)) * (exts[1] - exts[0]) + exts[0]
            seed_y[rnd_inds, i] = np.random.random(int(seed_cnt * rnd_val)) * (exts[3] - exts[2]) + exts[2]
            seed_x[rnd_inds, i], seed_y[rnd_inds, i] = ensure_in_poly(seed_x[rnd_inds, i], seed_y[rnd_inds, i], exts,
                                                                      mask_delauny, weighting=vel_weight)

        if move_along:
            if i > 10:
                pos_change_x = np.mean(np.diff(seed_x[:, i-10:i], n=1, axis=1), axis=1)
                pos_change_y = np.mean(np.diff(seed_y[:, i - 10:i], n=1, axis=1), axis=1)
                pos_change = np.sqrt(pos_change_x**2 + pos_change_y**2)
                stationary, = np.where(pos_change < 0.01)
                seed_x[stationary, i] = np.random.random(stationary.size) * (exts[1] - exts[0]) + exts[0]
                seed_y[stationary, i] = np.random.random(stationary.size) * (exts[3] - exts[2]) + exts[2]
                seed_x[stationary, i], seed_y[stationary, i] = ensure_in_poly(seed_x[stationary, i],
                                                                              seed_y[stationary, i], exts, mask_delauny,
                                                                              weighting=vel_weight)

        intF_u = RectBivariateSpline(x[:, 0], y[0, :], u[..., i])
        intF_v = RectBivariateSpline(x[:, 0], y[0, :], v[..., i])

        seed_u[..., i] = intF_u(seed_x[..., i], seed_y[..., i], grid=False)
        seed_v[..., i] = intF_v(seed_x[..., i], seed_y[..., i], grid=False)
        seed_mag[..., i] = np.sqrt(seed_u[..., i] ** 2 + seed_v[..., i] ** 2)

        seed_x[..., i + 1] = seed_x[..., i] + seed_u[..., i] * tRes * unit_factor
        seed_y[..., i + 1] = seed_y[..., i] + seed_v[..., i] * tRes * unit_factor

        if vel_weighting:
            if i > 10:
                u_mean = seed_u[:, i - 10:i].mean()
                vel_weight[0] = -u_mean * tRes * unit_factor
                v_mean = seed_v[:, i - 10:i].mean()
                vel_weight[1] = -v_mean * tRes * unit_factor

        seed_x[..., i+1], seed_y[..., i+1] = ensure_in_poly(seed_x[..., i+1], seed_y[..., i+1], exts, mask_delauny,
                                                            weighting=vel_weight)


    return (seed_x, seed_y, seed_u, seed_v)


def ensure_in_poly(x, y, exts, mask, weighting=(0, 0)):
    all_in_poly = False
    while not all_in_poly:
        invalid, = np.where(np.invert(in_hull(np.vstack((x, y)).T, mask)))
        x[invalid] = np.random.random(invalid.size) * (exts[1] - exts[0]) + exts[0] + weighting[0]
        y[invalid] = np.random.random(invalid.size) * (exts[3] - exts[2]) + exts[2] + weighting[1]
        all_in_poly = invalid.size == 0  # loop until only valid seeds are present

    return (x, y)
