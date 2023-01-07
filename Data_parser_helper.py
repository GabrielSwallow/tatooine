import re
import numpy as np
import UI_helper
import Navigation_helper

def get_planet_masses(data_name: str):
    directories = Navigation_helper.Directories(data_name)
    with open(directories.planet_ini) as f:
        lines = f.readlines()
        nbodies = lines[21:]
        masses = [1e-4]
        # masses = [float(nbody_data.split()[0]) for nbody_data in nbodies]
    return masses

def getAnalysisOutMetaInfo(data_name: str):
    '''
    returns 
    (
        dbl_num_orbits_per_out,
        dbl_num_azimuthal_angle_per_out,
        analysis_num_azimuthal_angle_per_out,
    )
    '''
    directories = Navigation_helper.Directories(data_name)
    with open(directories.pluto_log_filename) as dbl_txt:
        lines = dbl_txt.readlines()

        dbl_str = lines[79]
        dbl_list = list(map(int, re.findall('\d+', dbl_str)))
        dbl_num_orbits_per_out = dbl_list[3]
        dbl_num_azimuthal_angle_per_out = dbl_list[2]

        analysis_str = lines[86]
        analysis_list = list(map(int, re.findall('\d+', analysis_str)))
        analysis_num_azimuthal_angle_per_out = analysis_list[2]
    
    return (
        dbl_num_orbits_per_out,
        dbl_num_azimuthal_angle_per_out,
        analysis_num_azimuthal_angle_per_out,
    )


    
def getNbodyInformation_out(data_name: str, obj: int):
    '''
    returns 
    (
        time,
        a,
        e,
        omega,
        anomoly,
    )
    '''
    directories = Navigation_helper.Directories(data_name)

    dbl_num_orbits_per_out, _, _ = getAnalysisOutMetaInfo(data_name)

    (
        unfiltered_time,
        object_id,
        unfiltered_a,
        unfiltered_e,
        unfiltered_omega,
        unfiltered_anomaly,
    ) = np.loadtxt(directories.nbody_elements_filename, usecols=(0,1,2,3,6,7), unpack = True)

    time = dbl_num_orbits_per_out * unfiltered_time[object_id == obj]

    if len(time) == 0: 
        max_object_id = findNumBodies(directories.out_dir) - 1
        print('\nMax object_id = {}. You selected {}. Try again'.format(max_object_id, obj))
        obj, obj_des = UI_helper.selectObjectToPlot()
        getNbodyInformation_out(data_name, obj)
        return
    
    return (
        time,
        unfiltered_a[object_id == obj],
        unfiltered_e[object_id == obj],
        unfiltered_omega[object_id == obj],
        unfiltered_anomaly[object_id == obj]
    )

def getNbodyInformation_dat(data_name: str, obj: int):
    '''
    returns 
    (
        time (T_bin),
        a (a_bin),
        e,
        anomoly,
    )
    '''
    directories = Navigation_helper.Directories(data_name)

    (
        object_id,
        unfiltered_time,
        unfiltered_a,
        unfiltered_e,
        unfiltered_anomoly,
    ) = np.loadtxt(directories.nbody_elements_data_filename, skiprows=12, usecols=(0,1,2,3,7), unpack = True)

    # TODO: make this 10 be the number of out measurements per orbit
    time = unfiltered_time[object_id == obj] / (2*np.pi)

    if len(time) == 0: 
        max_object_id = findNumBodies(directories.out_dir)
        print('\nMax object_id = {}. You selected {}. Try again'.format(max_object_id, obj))
        getNbodyInformation_dat(directories.data_name, obj)
        return
    
    return (
        time,
        unfiltered_a[object_id == obj],
        unfiltered_e[object_id == obj],
        unfiltered_anomoly[object_id == obj],
    )

def getNbodyCoordinates(data_name: str, obj: int):
    '''
    returns 
    (
        time,
        x,
        y,
        z,
        vx,
        vy,
        vz,
    )
    '''
    directories = Navigation_helper.Directories(data_name)

    (
        object_id,
        unfiltered_time,
        unfiltered_x,
        unfiltered_y,
        unfiltered_z,
        unfiltered_vx,
        unfiltered_vy,
        unfiltered_vz,
    ) = np.loadtxt(directories.nbody_elements_coordinates_filename, skiprows=9, unpack = True)

    time = unfiltered_time[object_id == obj] / (2*np.pi)

    if len(time) == 0: 
        max_object_id = findNumBodies(directories.out_dir)
        print('\nMax object_id = {}. You selected {}. Try again'.format(max_object_id, obj))
        getNbodyCoordinates(directories.data_name, obj)
        return
    
    return (
        time,
        unfiltered_x[object_id == obj],
        unfiltered_y[object_id == obj],
        unfiltered_z[object_id == obj],
        unfiltered_vx[object_id == obj],
        unfiltered_vy[object_id == obj],
        unfiltered_vz[object_id == obj],
    )

def findNumBodies(out_dir: str) -> int:
    # directories = Navigation_helper.Directories(data_name)
    _, pid, _, _ = np.loadtxt('{}/nbody.out'.format(out_dir), usecols=(0,1,3,4), unpack=True)
    return len(np.unique(pid))

def get_averages_data(data_name: str):
    '''
    returns
    (
        time,
        disk_mass,
        disk_eccentricity,
        disk_periapsis,
        inner_disk_eccentricity,
        inner_disk_periapsis,
        sigma_min,
        iter,
        accretion
    )
    '''
    # TODO: once we add accretion of multiple planets, the last few columns will be optional 
    # depending on how many objects we have. Make this general
    directories = Navigation_helper.Directories(data_name)
    return np.loadtxt(directories.averages, skiprows=9, unpack = True)
