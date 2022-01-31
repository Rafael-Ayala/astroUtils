import numpy as np

def gravityGradientSphericalCoords(Pnm, dPnm, Cnm, Snm, lat, lon, r, R, GM, n, m):
    dUr = 0
    dUlat = 0
    dUlon = 0
    g1 = -GM/(r**2)
    g2 = GM/r
    for i in range(0, n+1):
        x1 = g1 * (i+1) * (R/r)**i
        x2 = g2 * (R/r)**i
        for j in range(0, i+1):
            dUr += x1*Pnm[i,j]*(Cnm[i,j]*np.cos(j*lon) + Snm[i,j]*np.sin(j*lon))
            dUlat += x2*dPnm[i,j]*(Cnm[i,j]*np.cos(j*lon) + Snm[i,j]*np.sin(j*lon))
            dUlon += x2*j*Pnm[i,j]*(Snm[i,j]*np.cos(j*lon) - Cnm[i,j]*np.sin(j*lon))
    return dUr, dUlat, dUlon

