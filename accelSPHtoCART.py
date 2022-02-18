import numpy as np

def accelSPHtoCART(position_bodyFixed, dUdr, dUdlat, dUdlon):
    d = np.sqrt(position_bodyFixed[0]**2 + position_bodyFixed[1]**2 + position_bodyFixed[2]**2)
    r2xy = position_bodyFixed[0]**2 + position_bodyFixed[1]**2
    ax = (1/d*dUdr-position_bodyFixed[2]/(d**2*np.sqrt(r2xy))*dUdlat)*position_bodyFixed[0]-(1/r2xy*dUdlon)*position_bodyFixed[1]
    ay = (1/d*dUdr-position_bodyFixed[2]/(d**2*np.sqrt(r2xy))*dUdlat)*position_bodyFixed[1]+(1/r2xy*dUdlon)*position_bodyFixed[0]
    az = 1/d*dUdr*position_bodyFixed[2]+np.sqrt(r2xy)/d**2*dUdlat
    return ax,ay,az
