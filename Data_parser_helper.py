import re
import numpy as np
import UI_helper

def getNbodyInformation_out(out_dir: str, obj: int):
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
    pluto_log_filename = out_dir + '/pluto.log'
    nbody_elements_filename = out_dir + '/nbody_orbital_elements.out'
    
    dbl_txt = open(pluto_log_filename)
    dbl_strs = dbl_txt.readlines()
    dbl_str = dbl_strs[79]
    dbl_list = list(map(int, re.findall('\d+', dbl_str)))
    dbl_time = dbl_list[3]
    dbl_txt.close()

    (
        unfiltered_time,
        object_id,
        unfiltered_a,
        unfiltered_e,
        unfiltered_omega,
        unfiltered_anomaly,
    ) = np.loadtxt(nbody_elements_filename, usecols=(0,1,2,3,6,7), unpack = True)

    time = dbl_time * unfiltered_time[object_id == obj]

    if len(time) == 0: 
        max_object_id = findNumBodies(out_dir) - 1
        print('\nMax object_id = {}. You selected {}. Try again'.format(max_object_id, obj))
        obj, obj_des = UI_helper.selectObjectToPlot()
        getNbodyInformation_out(out_dir, obj)
        return
    
    return (
        time,
        unfiltered_a[object_id == obj],
        unfiltered_e[object_id == obj],
        unfiltered_omega[object_id == obj],
        unfiltered_anomaly[object_id == obj]
    )

def getNbodyInformation_dat(out_dir: str, obj: int):
    '''
    returns 
    (
        time (T_bin),
        a (a_bin),
        e,
        period
    )
    '''
    nbody_elements_data_filename = out_dir + '/nbody_orbital_elements.dat'

    (
        object_id,
        unfiltered_time,
        unfiltered_a,
        unfiltered_e,
        unfiltered_period,
    ) = np.loadtxt(nbody_elements_data_filename, skiprows=12, usecols=(0,1,2,3,8), unpack = True)

    # TODO: make this 10 be the number of out measurements per orbit
    time = unfiltered_time[object_id == obj] / (2*np.pi)

    if len(time) == 0: 
        max_object_id = findNumBodies(out_dir)
        print('\nMax object_id = {}. You selected {}. Try again'.format(max_object_id, obj))
        getNbodyInformation_dat(out_dir, obj)
        return
    
    return (
        time,
        unfiltered_a[object_id == obj],
        unfiltered_e[object_id == obj],
        unfiltered_period[object_id == obj],
    )

def findNumBodies(out_dir: str) -> int:
    _, pid, _, _ = np.loadtxt('{}/nbody.out'.format(out_dir), usecols=(0,1,3,4), unpack=True)
    return len(np.unique(pid))