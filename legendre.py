import numpy as np

def legendre(n, m, lat):
    cosLat = np.cos(lat)
    sinLat = np.sin(lat)
    pnm = np.zeros((n+1, m+1))
    dpnm = np.zeros((n+1, m+1))
    pnm[0,0] = 1
    pnm[1,1] = np.sqrt(3)*cosLat
    dpnm[1,1] = -np.sqrt(3)*sinLat
    # Diagonal coefficients
    for i in range(2, n+1):
        pnm[i,i] = np.sqrt((2*i+1)/(2*i))*cosLat*pnm[i-1,i-1]
        dpnm[i,i] = np.sqrt((2*i+1)/(2*i))*(cosLat*dpnm[i-1,i-1] - sinLat*pnm[i-1,i-1])
    # Horizontal step 1 (coefficients adjacent to diagonals)
    for i in range(1, n+1):
        pnm[i,i-1] = np.sqrt(2*i+1)*sinLat*pnm[i-1,i-1]
        dpnm[i,i-1] = np.sqrt(2*i+1)*(cosLat*pnm[i-1,i-1] + sinLat*dpnm[i-1,i-1])
    # All other row-wise recursions:
    for j in range(0, n-2):
        for i in range(2+j, n):
            pnm[i,j] = (np.sqrt((2*i+1)/((i-j)*(i+j)))*
                ((np.sqrt(2*i-1)*sinLat*pnm[i-1,j]) -
                (np.sqrt(((i+j-1)*(i-j-1))/(2*i-3))*pnm[i-2,j])))
            dpnm[i,j] = (np.sqrt((2*i+1)/((i-j)*(i+j)))*
                ((np.sqrt(2*i-1)*sinLat*dpnm[i-1,j]) +
                (np.sqrt(2*i-1)*cosLat*pnm[i-1,j]) -
                (np.sqrt(((i+j-1)*(i-j-1))/(2*i-3))*dpnm[i-2,j])))
    return pnm, dpnm

