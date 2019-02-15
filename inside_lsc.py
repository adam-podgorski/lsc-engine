def enterLSC(launch_parameters,i_vec,n1=1,n2=1.5):
    R = sn.reflRatio(n1,n2,theta_i)
    n_vec = (0,0,1)
    if R == 1:
        return False
    else:
        return sn.refract(n1,n2,i_vec,n_vec,theta_i)

def checkAbsorption():

def checkNextBoundary(point,vector):