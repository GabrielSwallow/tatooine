import pluto
import numpy as np

colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff00', #'Gold', 
          '#a65628', '#f781bf', '#00b38f', '#6600cc']

# Physical constants in cgs units
constants = {
    'G': 6.6726e-8, # Gravitational Constant in cm^3 / (g s^2)
    'MSUN': 2.e33, # Solar Mass in g
    'RSUN': 6.96e10, # Solar Radius in cm
    'MEARTH': 5.9736e27, # Earth Mass in g
    'REARTH': 6.378136e8, # Earth Radius in cm
    'MJUPITER': 1.899e30, # Jupiter Mass in g
    'RJUPIRER': 6.9911e9, # Jupiter Radius in cm
    'AU': 1.496e13, # Astronomical unit in cm
    'LY': 9.4607e17, # Light-year in cm
    'PC': 3.0857e18, # Parsec in cm
    'MH': 1.6733e-24, # Hydrogen atom mass in g
    'KB': 1.3806505e-16, #Boltzmann constant in erg / K         
}

Kepler_47_constants = {
    'a_b_AU': 0.08145,
    'a_b_cm': 0.08145 * 1.496e+13,
    'period_days': 7.4483648,
    'period_s': 7.4483648 * 86400,
}

# binary = {'a': 1.0, 'e': 0.5, 'q': 0.1, 'm1': 2.4, 'm2': 0.24}
binary = {
    'a': 1.0,
    'e': 0.2,
    'q': 0.1,
    'm1': 0.957,
    'm2': 0.342,
}

def moving_average(array, n=7):
    weights = np.ones(n)/n
    return np.convolve(weights, array)[n-1:-n+1]


def get_t0_in_years(a, M):
    t0_in_seconds = np.sqrt( (a*constants['AU'])**3 
                            /(constants['G']*M*constants['MSUN']))
    year_in_seconds = 3.1556e7
    return t0_in_seconds / year_in_seconds


def get_t0_in_orbits(a):
    return 1.0/(2*np.pi*np.sqrt(a**3))


def orbits_to_years(orbits, a, M):
    year_in_seconds = 3.1556e7
    orbit_in_seconds = 2*np.pi*np.sqrt( (a*constants['AU'])**3 
                                       /(constants['G']*M*constants['MSUN']))
    return orbits*orbit_in_seconds / year_in_seconds


def years_to_orbits(years, a, M):
    year_in_seconds = 3.1556e7
    orbit_in_seconds = 2*np.pi*np.sqrt( (a*constants['AU'])**3
                                       /(constants['G']*M*constants['MSUN']))
    return years*year_in_seconds / orbit_in_seconds


def solve_kepler_equation(e, a, t, max_iter=200, epsilon=1.0e-9):
    mean_anomaly = np.sqrt(1.0/(a**3))*t
    ecc_anomaly = mean_anomaly

    for i in range(1, max_iter+1):
        ff = ecc_anomaly - e*np.sin(ecc_anomaly) - mean_anomaly
        fs = 1.0 - e*np.cos(ecc_anomaly)
        new_ecc_anomaly = ecc_anomaly - ff/fs
        err = ecc_anomaly - new_ecc_anomaly

        if (i > 2 and np.absolute(err) < epsilon):
            break

        ecc_anomaly = new_ecc_anomaly

    return ecc_anomaly


def kepler_position(e, a, t):
    ecc_anomaly = solve_kepler_equation(e, a, t)
    phi = (2.0*np.arctan2(np.sqrt(1.0 + e) * np.sin(0.5*ecc_anomaly),
                          np.sqrt(1.0 - e) * np.cos(0.5*ecc_anomaly)))
    if phi < 0.0:
        phi += 2*np.pi

    r = (a * (1.0 -e **2)) / (1.0 + e * np.cos(phi))

    return r, phi


def kepler_orbit(binary):
    phi = np.linspace(0.0, 2.0*np.pi, 200)
    r = binary['a'] * (1.0 - binary['e']**2) / (1.0 + binary['e']*np.cos(phi))

    return r, phi

def binary_position(binary, r, phi):
    q = binary['q']

    x1 = q/(1.0+q) * r * np.cos(phi)
    y1 = q/(1.0+q) * r * np.sin(phi)

    x2 = -1.0/(1.0+q) * r * np.cos(phi)
    y2 = -1.0/(1.0+q) * r * np.sin(phi)

    return x1, y1, x2, y2


def azimuthal_average(var, dx2):
    avg = var * dx2
    return np.sum(avg, axis=0) / (2*np.pi)


def fargo_azimuthal_avg(var, Nr, Nphi):
    phi = np.linspace(0.0, 2*np.pi, Nphi+1, endpoint=True)

    avg = np.zeros(Nr)

    for i in range(Nr):
        for j in range(Nphi):
            avg[i] += var[j,i] * (phi[j+1] - phi[j])
        avg[i] /= 2*np.pi
    
    return avg


def orbital_elements(R, phi, vR, vphi):
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)

    factor1 = -R*vR*vphi
    factor2 = R*vphi*vphi - 1.0

    e_x = -factor1*sin_phi + factor2*cos_phi
    e_y =  factor1*cos_phi + factor2*sin_phi

    e = np.sqrt(e_x*e_x + e_y*e_y)
    pericenter = np.arctan2(e_y, e_x)
    pericenter = np.where(pericenter < 0.0, pericenter + 2*np.pi, pericenter)

    a = 1.0 / (2.0/R - (vR*vR + vphi*vphi))

    return a, e, pericenter


def radial_orbital_elements(sigma, vR, vphi, grid):
    R, phi = np.meshgrid(grid['centers'][0], grid['centers'][1])
    dVx1, dVx2 = np.meshgrid(grid['dV'][0], grid['dV'][1])
    dV = dVx1*dVx2

    a, e, pericenter = orbital_elements(R, phi, vR, vphi)

    dm = sigma*dV

    mass_ring = np.sum(dm, axis=0)
    ecc_ring = np.sum(e*dm, axis=0)/mass_ring
    peri_ring = np.sum(pericenter*dm, axis=0)/mass_ring

    return ecc_ring, peri_ring

def mass_weighted_orbital_elements(sigma, vR, vphi, grid, Rmax):
    R, phi = np.meshgrid(grid['centers'][0], grid['centers'][1])
    dVx1, dVx2 = np.meshgrid(grid['dV'][0], grid['dV'][1])
    dV = dVx1*dVx2

    a, e, pericenter = orbital_elements(R, phi, vR, vphi)

    dm = np.where(R <= Rmax, sigma*dV, 0.0)
    disc_mass = np.sum(dm)
    disc_ecc = np.sum(e*dm)/disc_mass
    disc_pericenter = np.sum(pericenter*dm)/disc_mass

    return disc_ecc, disc_pericenter


def gap_size(radius, avg_sigma):
    max_index = np.argmax(avg_sigma)

    sigma_peak = avg_sigma[max_index]
    sigma_05 = 0.5*sigma_peak
    sigma_01 = 0.1*sigma_peak

    R_peak = radius[max_index]
    R_05 = 0.0
    R_01 = 0.0

    for R, sigma in zip(radius, avg_sigma):
        if sigma < sigma_05:
            continue
        else:
            R_05 = R
            break

    for R, sigma in zip(radius, avg_sigma):
        if sigma < sigma_01:
            continue
        else:
            R_01 = R
            break

    return (R_peak, R_05, R_01)


def get_tgap_in_orbits(t, w, start_orbit=6000, offset=1000, min_diff=200):
    indices_forward = []
    indices_backward = []
    for k in range(offset, len(w)-offset):
        if t[k] < start_orbit:
            continue

        if (w[k]   <= np.pi and 
            w[k+1] >  np.pi and
            w[k-offset] <= np.pi and
            w[k+offset] >  np.pi):
            indices_forward.append(k)

        if (w[k]   >= np.pi and 
            w[k+1] <  np.pi and
            w[k-offset] >= np.pi and
            w[k+offset] <  np.pi):
            indices_backward.append(k)

    n, s = 0, 0
    for k in range(len(indices_forward)-1):
        diff = t[indices_forward[k+1]] - t[indices_forward[k]]
        if diff > min_diff:
            s += diff
            n += 1
    tgap_forward = s / n
    print('Forward periods {0:d}'.format(n))

    n, s  = 0, 0
    for k in range(len(indices_backward)-1):
        diff = t[indices_backward[k+1]] - t[indices_backward[k]]
        if diff > min_diff:
            s += diff
            n += 1
    tgap_backward = s / n
    print('Backward periods {0:d}'.format(n))

    tgap = 0.5*(tgap_forward + tgap_backward)
    print('T_gap: {0:.2f} orbits\n'.format(tgap))
    return tgap


def get_tgap_convergence(t, w, start_orbit=0, offset=1000, min_diff=200):
    indices_forward = []
    indices_backward = []
    for k in range(offset, len(w)-offset):
        if t[k] < start_orbit:
            continue

        if (w[k]   <= np.pi and 
            w[k+1] >  np.pi and
            w[k-offset] <= np.pi and
            w[k+offset] >  np.pi):
            indices_forward.append(k)

        if (w[k]   >= np.pi and 
            w[k+1] <  np.pi and
            w[k-offset] >= np.pi and
            w[k+offset] <  np.pi):
            indices_backward.append(k)

    n_f = []
    diff_f = []
    count = 0
    for k in range(len(indices_forward)-1):
        diff = t[indices_forward[k+1]] - t[indices_forward[k]]
        if diff > min_diff:
            diff_f.append(diff)
            n_f.append(count)
            count += 1

    n_b = []
    diff_b = []
    count = 0
    for k in range(len(indices_backward)-1):
        diff = t[indices_backward[k+1]] - t[indices_backward[k]]
        if diff > min_diff:
            diff_b.append(diff)
            n_b.append(count)
            count += 1

    return (np.array(n_f), np.array(diff_f), np.array(n_b), np.array(diff_b))


def calc_gap_ellipse(sigma, radius, phi, fraction=0.1):
    NX2 = phi.size
    jmax, imax = np.unravel_index(sigma.argmax(), sigma.shape)
    sigma_max = sigma[jmax, imax]
    
    R_perihel = 0.0
    for R, sig, in zip(radius, sigma[jmax, :]):
        if sig < fraction * sigma_max:
            continue
        else:
            R_aphel = R
            break

    if jmax >= NX2//2:
        jperi = jmax - NX2//2
    else:
        jperi = jmax + NX2//2

    for R, sig in zip(radius, sigma[jperi, :]):
        if sig < fraction * sigma_max:
            continue
        else:
            R_perihel = R
            break

    if R_perihel == 0.0:
        for R, sig in zip(radius, sigma[jperi, :]):
            if sig < 0.05*sigma_max:
                continue
            else:
                R_perihel = R
                break

    phi_0 = phi[jmax]
    a = 0.5*(R_aphel + R_perihel)
    e = R_aphel - a
    b = np.sqrt(a**2 - e**2)
    x_0 = e*np.cos(phi_0)
    y_0 = e*np.sin(phi_0)

    return x_0, y_0, phi_0, a, b, e/a
    

def get_ellipse_points(x0, y0, phi0, a, b):
    phi = np.arange(0.0, 2.0*np.pi, 0.01)
    x_ell = x0 + a*np.cos(phi)*np.cos(phi0) - b*np.sin(phi)*np.sin(phi0)
    y_ell = y0 + a*np.cos(phi)*np.sin(phi0) + b*np.sin(phi)*np.cos(phi0)

    return x_ell, y_ell

def kepler_velocity(radius: float) -> float: 
    # cm_per_a_b = Kepler_47_constants['a_b_cm']
    # s_per_period = Kepler_47_constants['period_s']
    # length_scale = data.read_units()['length']
    # mass_scale = data.read_units()['length']

    central_mass = binary['m1'] + binary['m2']
    # speed_in_cm_per_s = np.sqrt(constants['G']*(central_mass*mass_scale)/ (radius*length_scale))
    speed = np.sqrt(central_mass/ radius)
    return speed
    # speed = speed_in_a_b_per_period * (cm_per_a_b) * (1/s_per_period)
    # angular_velocity = speed / radius
    # return angular_velocity

def x_coord(r: float, phi: float) -> float:
    return r*np.cos(phi)

def y_coord(r: float, phi: float) -> float:
    return r*np.sin(phi)

def r_coord(x: float, y: float) -> float:
    return np.sqrt(x**2 + y**2)

def time_split(time: list, t_min: int, t_max: int) -> list:
    for i in range(len(time)):
        i_min = i
        if time[i] >= t_min: break
        else: continue
    
    for i_inv in range(len(time)):
        i = len(time) - 1 - i_inv
        i_max = i
        if time[i] <= t_max: break
        else:continue
    
    return i_min, i_max
