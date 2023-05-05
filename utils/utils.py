import numpy as np

def dipole_basis_vectors(lat,lon=None,cartesian=False):
    """
    qhat, phat = dipole_basis_vectors(lat)

    INPUT
    =====
    lat : (N,) Array of magnetic latitudes in degrees

    OUTPUT
    ======
    qhat : (2,N) Field-aligned basis vector in spherical coordinates (r component, then theta component)
    phat : (2,N) Field-perp basis vector in spherical coordinates (r component, then theta component)
    """

    latshape = lat.shape
    lat = lat.ravel()

    if cartesian:
        assert lon is not None,"Must provide longitude to get Cartesian components"
        assert lon.shape == latshape

        lon = lon.ravel()

    latr = np.radians(lat)
    denom = np.sqrt(1.+3.*np.sin(latr)**2)

    qhatr = -2. * np.sin(latr)/denom
    qhatt = -np.cos(latr)/denom

    phatr = np.cos(latr)/denom
    phatt = -2 * np.sin(latr)/denom

    if cartesian:
        thetar = np.pi/2 - latr
        phir = np.radians(lon)
        qhatx = qhatr*np.sin(thetar)*np.cos(phir) + qhatt * np.cos(thetar)*np.cos(phir)
        qhaty = qhatr*np.sin(thetar)*np.sin(phir) + qhatt * np.cos(thetar)*np.sin(phir)
        qhatz = qhatr*np.cos(thetar) - qhatt*np.sin(thetar)

        phatx = phatr*np.sin(thetar)*np.cos(phir) + phatt * np.cos(thetar)*np.cos(phir)
        phaty = phatr*np.sin(thetar)*np.sin(phir) + phatt * np.cos(thetar)*np.sin(phir)
        phatz = phatr*np.cos(thetar) - phatt*np.sin(thetar)

        qhat = np.vstack([qhatx,qhaty,qhatz]).T
        phat = np.vstack([phatx,phaty,phatz]).T

        qhat = qhat.reshape((*latshape,3))
        phat = phat.reshape((*latshape,3))

        return (qhat,phat)

    else:
        qhat = np.vstack([qhatr,qhatt]).T
        phat = np.vstack([phatr,phatt]).T
        qhat = qhat.reshape((*latshape,2))
        phat = phat.reshape((*latshape,2))

        return (qhat,phat)


def ECEF2geodetic(x, y, z, degrees = True):
    """ 
    convert from x, y, z ECEF to geodetic coordinates, h, lat, lon.
    returns h, lat, lon

    input should be given in km

    Using the Zhu algorithm described in:
    J. Zhu, "Conversion of Earth-centered Earth-fixed coordinates 
    to geodetic coordinates," IEEE Transactions on Aerospace and 
    Electronic Systems, vol. 30, pp. 957-961, 1994.
    """

    d2r = np.pi/180
    r2d = 180  /np.pi
    WGS84_e2 = 0.00669437999014
    WGS84_a  = 6378.137

    conv = r2d if degrees else 1.

    # compute the lon straighforwardly
    lon = ((np.arctan2(y, x)*180/np.pi) % 360)/180*np.pi * conv

    # input params
    e = np.sqrt(WGS84_e2)
    a = WGS84_a
    b = a*np.sqrt(1 - e**2)
    r = np.sqrt(x**2 + y**2)

    # and enter the algorithm...
    l = e**2 / 2.
    m = (r / a)**2
    n = ((1 - e**2) * z / b)**2
    i = -(2 * l**2 + m + n) / 2.
    k = l**2 * (l**2 - m - n)
    q = (m + n - 4 * l**2)**3/216. + m * n * l**2
    D = np.sqrt((2 * q - m * n * l**2) * m * n * l**2 + 0j)
    beta = i/3. - (q + D)**(1/3.) - (q - D)**(1/3.) 
    t = np.sqrt(np.sqrt(beta**2 - k) - (beta + i)/2) - np.sign(m - n) * np.sqrt((beta - i) / 2)
    r0 = r / (t + l)
    z0 = (1 - e**2) * z / (t - l)
    lat = np.arctan(z0 / ((1 - e**2) * r0)) * conv
    h = np.sign(t - 1 + l) * np.sqrt((r - r0)**2 + (z - z0)**2)

    return h.real, lat.real, lon



def spherical_to_dipole(r,theta):
    """
    r normalized: r = R/R_E
    theta in radians
    """
    q = np.cos(theta) / r**2
    p = r/np.sin(theta)**2

    return q,p

# transformation equations for dipole coordinates from Kageyama et al (2005)
def dipole_to_spherical(q,p):
        """
        q = (R_E/r)^2 * cos(theta)
        p = (r/R_E)/sin^2(theta)

        NOTE: These differ slightly from definitions of dipole coordinates given by Kageyama et al (2005), who define μ = - cos(theta)/r, and Χ = sin^2(theta)/r, such that

        q = -μ;
        p = 1/Χ
        """

        c1 = 2**(7/3)/3**(1/3)
        c2 = 2**(1/3)*3**(2/3)

        zeta = (q*p**2)**2

        gamma = lambda zeta: (9*zeta + np.sqrt(3)*np.sqrt(27*zeta**2+256*zeta**3))**(1/3)
        wfunc = lambda zeta: -c1/gamma(zeta)+gamma(zeta)/c2/zeta
        ufunc = lambda zeta: -0.5*np.sqrt(wfunc(zeta))+0.5*np.sqrt(-wfunc(zeta)+2/zeta/np.sqrt(wfunc(zeta)))
        
        u = ufunc(zeta)
        r = u*p
        theta = np.arcsin(np.sqrt(u))

        return r,theta
