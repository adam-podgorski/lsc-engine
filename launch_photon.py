import numpy as np
import snells as sn

def LSCDimensions(lsc_dimensions):
    return float(lsc_dimensions[0]),float(lsc_dimensions[1]),float(lsc_dimensions[2])

def launchVariables(lsc_dimensions,lambda_cdf):
    # LSC Dimensions
    l_x,l_y,l_z = LSCDimensions(lsc_dimensions)
    # Wavelengths
    wavelengths = list(lambda_cdf[0])
    cdf = list(lambda_cdf[1])
    # Random launch variables, cartesian coordinates of the photon, 
    # plus azimuthal and incident angle 
    x_i = np.random.random_sample()*l_x
    y_i = np.random.random_sample()*l_y
    z_i = l_z
    gamma_i = np.random.normal()*2*np.pi # azimuthal angle
    theta_i = 0.501*np.pi+np.random.normal()*np.pi # incident angle
    lambda_i = np.interp(np.random.random(),cdf,wavelengths)
    return np.array([x_i,y_i,z_i,gamma_i,theta_i,lambda_i])

def maxDomain(lsc_dimensions,max_offsets=[1.2,1.2,2.0]):
    return np.array([max_offsets[i]*lsc_dimensions[i] for i in range(3)])

def vectorIntercepts(launch_parameters,bc_domain):
    x_i,y_i,z_i,gamma_i,theta_i = launch_parameters[0:5]
    x_max,y_max,z_max = bc_domain

    t_s = [0]*5 # 5 Values, three cases
    t_s[0] = x_i/(np.cos(gamma_i)*np.sin(theta_i))
    t_s[1] = (x_max-x_i)/(np.cos(gamma_i)*np.sin(theta_i))
    t_s[2] = y_i/(np.sin(gamma_i)*np.sin(theta_i))
    t_s[3] = (y_max-y_i)/(np.sin(gamma_i)*np.sin(theta_i))
    t_s[4] = (z_max-z_i)/(np.cos(theta_i))

    # Get smallest negative value of t
    array = np.array(t_s)
    t_intersect_idx = list(array).index(max(array[array<0]))

    # Compute boundary points using solar geometry
    x_0 = x_i + np.cos(gamma_i)*np.sin(theta_i)*t_s[t_intersect_idx]
    y_0 = y_i + np.sin(gamma_i)*np.sin(theta_i)*t_s[t_intersect_idx]
    z_0 = z_i + np.cos(theta_i)*t_s[t_intersect_idx]
    return np.array([x_0,y_0,z_0])

def launchPhoton(lsc_dimensions,lambda_cdf):
    # LSC Dimensions
    l_x,l_y,l_z = LSCDimensions(lsc_dimensions)
    # Only launch if photon is launched from certain range
    launch = True
    while launch:
        launch_parameters = launchVariables(lsc_dimensions,lambda_cdf)
        bc_domain = maxDomain(lsc_dimensions)
        boundary_pt = vectorIntercepts(launch_parameters,bc_domain)
        if not(boundary_pt[0] < -0.2*bc_domain[0]):
            launch = False
        if not(boundary_pt[1] < -0.2*bc_domain[1]):
            launch = False
    # Return point on LSC, vector and incident/azimuthal angles for later
    return launch_parameters, boundary_pt

def createIncidentVector(launch_pt,boundary_pt):
    # Requires launch_photon
    # Determine incident ray as a vector, normalised
    v_xi = launch_pt[0] - boundary_pt[0]
    v_yi = launch_pt[1] - boundary_pt[1]
    v_zi = launch_pt[2] - boundary_pt[2]
    incident_ray_vector = np.array([v_xi,v_yi,v_zi])/np.linalg.norm([v_xi,v_yi,v_zi])
    return incident_ray_vector

    
def plotIncidentRay(lsc_dimensions,launch_pt,boundary_pt,ax):  
    # Requires launch_photon
    # Plot line between these two points
    x_s = [boundary_pt[0], launch_pt[0]]
    y_s = [boundary_pt[1], launch_pt[1]]
    z_s = [boundary_pt[2], launch_pt[2]]
    ax.plot(x_s,y_s,z_s,color='red') 

def printLaunchParameters(launch_parameters):
    print(launch_parameters)
    print("Launch parameters: (x="+str(launch_parameters[0])
        +",y="+str(launch_parameters[1])
        +",z="+str(launch_parameters[2])+")\n")
    print("Azimuthal angle = "+str(launch_parameters[3]))
    print("Incidence angle = "+str(launch_parameters[4]))
    print("Wavelength = "+str(launch_parameters[5]))

