"""
Module to read PLUTO output data.
"""
import os
import sys
import re
import collections
import numpy as np

class Pluto:
    """
    Class for reading the various output files of the PLUTO code.

    Class variables:

    constants  -- dictionary containing physical constants in cgs units

    length                  cm
    mass                    g
    time                    s
    velocity                cm/s
    force                   dyne = g cm / s^2
    energy                  erg = g cm^2 /s^2
    pressure                barye = g / (cm s^2)
    dynamic viscosity       g / (cm s)
    kinematic viscosity     cm^2 / s

    Atomic mass unit:       AMU    = 1.66053886e-24     [g]
    Astronomical unit:      AU     = 1.49597892e13      [cm]
    Speed of light:         C      = 2.99792458e10      [cm/s]
    Electron volt:          Ev     = 1.602176463158e-12 [erg]
    Gravitational constant: G      = 6.6726e-8          [cm^3/(g s^2)]
    Planck constant:        H      = 6.62606876e-27     [erg s]
    Boltzmann constant:     Kb     = 1.3806505e-16      [erg / K]
    Light year              LY     = 0.9461e18          [cm]
    Proton mass             MP     = 1.67262171e-24     [g]
    Neutron mass            MN     = 1.67492728e-24     [g]
    Electron mass           ME     = 9.1093826e-28      [g]
    Hydrogen atom mass      MH     = 1.6733e-24         [g]
    Solar mass              MSUN   = 2.0e33             [g]
    Earth mass              MEARTH = 5.9736e27          [g]
    Jupiter mass            MJUP   = 1.899e30           [g]
    Avogadro constant       NA     = 6.0221367e23
    Parsec                  PC     = 3.0856775807e18    [cm]
    Earth radius            REARTH = 6.378136e8         [cm]
    Jupiter radius          RJUP   = 6.9911e9           [cm]
    Stephan-Boltzmann const SIGMA  = 5.67051e-5         [erg/(cm^2 K^4 s)]

    Instance variables:

    path                -- String containing the path to the read data.
    grid                -- Dictionary containing information about the grid.
    time                -- Dictionary containing information about the output times.
    units               -- Dictionary containing information about the code units.
    user_def_parameters -- Dictionary containing the user defined parameters
                           from pluto.ini
    last                -- Integer containing the number of the last output file.
    """
    constants = {'AMU'    : 1.66053886e-24,
                 'AU'     : 1.49597892e13,
                 'C'      : 2.99792458e10,
                 'Ev'     : 1.602176463158e-12,
                 'G'      : 6.6726e-8,
                 'H'      : 6.62606876e-27,
                 'KB'     : 1.3806505e-16,
                 'LY'     : 0.9461e18,
                 'MP'     : 1.67262171e-24,
                 'MN'     : 1.67492728e-24,
                 'ME'     : 9.1093826e-28,
                 'MH'     : 1.6733e-24,
                 'MSUN'   : 2.0e33,
                 'MEARTH' : 5.9736e27,
                 'MJUP'   : 1.899e30,
                 'NA'     : 6.0221367e23,
                 'PC'     : 3.0856775807e18,
                 'REARTH' : 6.378136e8,
                 'RSUN'   : 6.96e10,
                 'RJUP'   : 6.9911e9,
                 'SIGMA'  : 5.67051e-5}

    def __init__(self, path=os.getcwd()):
        self.path = path
        self.outdir = self.get_outdir()

        print("Reading from directory: {0}".format(self.path))
        self.grid = self.read_grid()
        self.time = self.read_time()
        self.units = self.read_units()
        self.user_def_parameters = self.get_user_def_parameters()

        self.last = int(self.time['t'].size-1)
        

    def get_outdir(self):
        """
        Reads the output directory from pluto.ini

        Returns
        String contatining the output directory
        """
        outdir = ''
        with open('{0}/../pluto.ini'.format(self.path), 'r') as f:
            for line in f:
                if 'output_dir' in line:
                    outdir = line.split()[1]
                    break

        return os.path.basename(outdir)


    def get_cell_volume(self, direction, geometry, xl, xr):
        """
        Calculates the one dimensional cell volume.

        Parameters:
        direction -- Integer for the direction (x1 == 0, x2 == 1, x3 == 3).
        geometry  -- String containing the geometry of the grid (CARTESIAN,
                     POLAR, CYLYNDRICAL, SPHERICAL).
        xl        -- Float containing the left interface of the cell.
        xr        -- Float containing the right interface of the cell.
        """
        if direction == 0:
            if geometry in ["POLAR", "CYLINDRYCAL"]:
                dV = 0.5*(xr**2 - xl**2)
            elif geometry == "SPHERICAL":
                dV = (xr**3 - xl**3)/3.0
            else:
                dV = xr-xl
        elif direction == 1:
            if geometry == "SPHERICAL":
                dV = np.cos(xl) - np.cos(xr)
            else:
                dV = xr-xl
        else:
            dV = xr-xl

        return dV

    def read_grid(self):
        """
        Reads information stored in grid.out.

        Returns:
        Dictionary containing:
            total   -- Integer containing total number of cells.
            cells   -- Named tupel containing number of cells in each
                       direction.
            faces   -- Named tupel containing 1D-numpy-arrays of cell faces in
                       each direction.
            centers -- Named tupel containing 1D-numpy-arrays of cell centers in
                       each direction.
            dx      -- Named tupel containing 1D-numpy-arrays of cell lengths.
            dV      -- Named tupel containing 1D-numpy-arrays of cell volumes.
        """

        gridfile_path = "{0}/grid.out".format(self.path)

        try:
            gridfile = open(gridfile_path, 'r')
        except:
            print("Could not open {0}".format(gridfile_path))
            sys.exit(1)
        else:
            print("Reading grid.out")
            dimensions = 0
            direction = 0
            number_of_cells = np.ones(3, dtype=np.int)
            Cells = collections.namedtuple('Cells', ['X1', 'X2', 'X3'])

            # Parse header
            for line in gridfile:
                if not line.startswith('#'):
                    break

                if "DIMENSIONS:" in line:
                    dimensions = int(line.split()[2])

                if "GEOMETRY:" in line:
                    geometry = line.split()[2]

                for d in range(dimensions):
                    if "X{0:d}:".format(d+1) in line:
                        number_of_cells[d] = int(re.search(', (\d+) point\(s\),', line).group(1))

            number_of_cells = Cells._make(number_of_cells)

            centers = Cells(np.zeros(number_of_cells.X1, dtype=np.float_),
                            np.zeros(number_of_cells.X2, dtype=np.float_),
                            np.zeros(number_of_cells.X3, dtype=np.float_))
            dx = Cells(np.zeros(number_of_cells.X1, dtype=np.float_),
                       np.zeros(number_of_cells.X2, dtype=np.float_),
                       np.zeros(number_of_cells.X3, dtype=np.float_))
            dV = Cells(np.zeros(number_of_cells.X1, dtype=np.float_),
                       np.zeros(number_of_cells.X2, dtype=np.float_),
                       np.zeros(number_of_cells.X3, dtype=np.float_))
            interfaces = Cells(np.zeros(number_of_cells.X1+1, dtype=np.float_),
                               np.zeros(number_of_cells.X2+1, dtype=np.float_),
                               np.zeros(number_of_cells.X3+1, dtype=np.float_))

            for line in gridfile:
                # Don't read comments
                if line.startswith('#'):
                    continue

                # Read interfaces
                try:
                    no, xl, xr = line.split()
                    index, xl, xr = int(no)-1, float(xl), float(xr)
                except:
                    direction += 1
                    continue

                # Calculate different cell values
                interfaces[direction][index] = xl
                centers[direction][index] = (xl + xr)*0.5
                dx[direction][index] = xr-xl
                dV[direction][index] = self.get_cell_volume(direction, 
                                                            geometry, xr, xl)

                # Append last interface
                if index+1 == number_of_cells[direction]:
                    interfaces[direction][index+1] = xr

            gridfile.close()

            return {'total'   : np.prod(number_of_cells),
                    'cells'   : number_of_cells,
                    'faces'   : interfaces,
                    'centers' : centers,
                    'dx'      : dx,
                    'dV'      : dV}

    def read_time(self):
        """
        Reads information stored in dlb.out.

        Returns:
        Dictionary containing:
            t        -- 1D-numpy-array containing the output times.
            dt       -- 1D-numpy-array containing the time steps at the output
                        times.
            nstep    -- 1D-numpy-array containing the number of steps at the
                        output times.
            format   -- String storing the file format (single_file or
                        multiple_files)
            variable -- List of strings storing the written variables
                        (e.g. ['rho', 'vx1', 'vx2', 'prs'])
        """
        timefile_path = "{0}/dbl.out".format(self.path)

        try:
            timefile = open(timefile_path, 'r')
        except:
            print("Could not open {0}".format(timefile_path))
            sys.exit(1)
        else:
            print("Reading dbl.out")

            tdata = [[], [], []]
            file_format = None
            variables = {}

            for i, line in enumerate(timefile):
                values = line.split()
                if i == 0:
                    file_format = values[4]

                    for j in range(6, len(values)):
                        variables[values[j]] = j-6

                tdata[0].append(float(values[1]))
                tdata[1].append(float(values[2]))
                tdata[2].append(float(values[3]))

            timefile.close

            return {'t'         : np.array(tdata[0]),
                    'dt'        : np.array(tdata[1]),
                    'nstep'     : np.array(tdata[2]),
                    'format'    : file_format,
                    'variables' : variables}

    def read_units(self):
        """
        Reads units in cgs from pluto.log

        Returns:
        code_units -- Dictionary containing units in cgs of the following
                      quantities:
                        - density                 [g/cm^3]
                        - pressure                [dyne/cm^2 = g/(cm * s^2)]
                        - length                  [cm]
                        - velocity                [cm/s]
                        - temperature * p/rho*mu  [K]
                        - time                    [s]
        """
        # Default units
        code_units = {}
        code_units['density'] = 1.67262171e-24
        code_units['length'] = 1.49597892e13
        code_units['velocity'] = 1.0e5

        unitfile_path = "{0}/pluto.log".format(self.path)

        try:
            unitfile = open(unitfile_path, 'r')
        except:
            print("Could not open {0}".format(unitfile_path))
            print("Using default units")
            code_units['pressure'] = (code_units['density']
                                      *code_units['velocity']**2)
            code_units['temperature'] = (code_units['velocity']**2
                                         *Pluto.constants['AMU']
                                         /Pluto.constants['KB'])
            code_units['time'] = code_units['length']/code_units['velocity']
            return code_units
        else:
            print("Reading pluto.log")

            for line in unitfile:
                if "[Density]:" in line:
                    code_units['density'] = (float(line.split()[1]))
                elif "[Pressure]:" in line:
                    code_units['pressure'] = (float(line.split()[1]))
                elif "[Velocity]:" in line:
                    code_units['velocity'] = (float(line.split()[1]))
                elif "[Length]:" in line:
                    code_units['length'] = (float(line.split()[1]))
                elif "[Temperature]:" in line:
                    code_units['temperature'] = (float(line.split()[1]))
                elif "[Time]:" in line:
                    code_units['time'] = (float(line.split()[1]))

            code_units['pressure'] = (code_units['density']
                                      *code_units['velocity']**2)
            code_units['time'] = code_units['length'] / code_units['velocity']
            code_units['mass'] = code_units['density']*code_units['length']**3

            return code_units

    def primitive_variable(self, variable, file_number=-1):
        """
        Reads cell-centered values of a primitive variable stored in
        variable.number.dbl.

        Parameters:   
        variable    -- String containing the name of the primitive variable
                       (e.g. 'rho', 'vx1').
        file_number -- Integer containing the file number. Default is the last
                       output file (file_number = -1)

        Returns:
        data -- 3D-Numpy-Array containing the cell-centered values of the
                primitive variable. Indexing data[k, j, i] with
                    i -- x1-direction
                    j -- x2-direction
                    k -- x3-direction
        """
        if file_number == -1:
            number = self.last
        else:
            number = file_number

        if self.time['format'] == 'single_file':
            binaryfile = "{0}/data.{1:04d}.dbl".format(self.path, number)
        else:
            binaryfile = "{0}/{1}.{2:04d}.dbl".format(self.path, variable, number)

        try:
            linear_data = np.fromfile(binaryfile, dtype=np.float_)
        except:
            print("Could not open {0}".format(binaryfile))
            raise
            #sys.exit(1)
        else:

            if self.time['format'] == 'single_file':
                print("Reading {0} from data.{1:04d}.dbl".format(variable, number))
                start = self.time['variables'][variable]*self.grid['total']
                end = (self.time['variables'][variable]+1)*self.grid['total']
                linear_data = linear_data[start:end]
            else:
                print("Reading {0}.{1:04d}.dbl".format(variable, number))

            data = linear_data.reshape(tuple(reversed(self.grid['cells'])))

            return data
            
    def all_variables(self, number=-1):
        """
        Reads cell-centered values of a primitive variable stored in
        variable.number.dbl.

        Parameters:
        file_number -- Integer containing the file number. Default is the last
                   output file (file_number = -1)

        Returns:
        a dict of {varname: data}
        where
        data -- 3D-Numpy-Array containing the cell-centered values of the
            primitive variable. Indexing data[k, j, i] with
                i -- x1-direction
                j -- x2-direction
                k -- x3-direction
        """

        if self.time['format'] == 'single_file':
            binaryfile = "{0}/data.{1:04d}.dbl".format(self.path, number)
        else:
            raise NotImplementedError

        if number % 10 == 0:
            print("Reading {}".format(binaryfile))

        vars = self.time['variables']
        shape = (len(vars),) + tuple(reversed(self.grid['cells']))
        linear_data = np.fromfile(binaryfile, dtype=np.float_)
        array_data = linear_data.reshape(shape)
        return array_data

    def get_user_def_parameters(self):
        """
        Reads the user-defined parameters from pluto.ini.
        
        Parameters:
        
        Returns:
        Dictionary containing the parameters
        """


        ini_file = '{0}/../pluto.ini'.format(self.path)
        sections = {}
        sections['[undefined]'] = current_section = {}
        with open(ini_file , 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('['):
                    sections[line] = current_section = {}
                elif line:
                    key, value = line.split(None, 1)
                    current_section[key] = value

        return {
            key: float(value)
            for key, value in sections['[Parameters]'].items()
        }
                

    def quantity(self, q):
        """
        Reads a user-defined quantity stored in pluto.ini.

        Parameters:
        q    -- String containing the name of the quantity as it appears in
                pluto.ini.

        Returns:
        Float containing the value of q
        """
        inifile_path = "{0}/../pluto.ini".format(self.path)

        try:
            inifile = open(inifile_path, 'r')
        except:
            print("Could not open {0}".format(inifile_path))
            sys.exit(1)
        else:
            print("Reading {0} from pluto.ini".format(q))

            for line in inifile:
                if q in line:
                    return float(line.split()[1])
