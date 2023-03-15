all_data_dir = '../Data/'
global_plots_dir = '../Plots/'
cart  = False
logsc = True
nbody = True
var   = "rho" # rho, vx1, vx2
nts   = 100 # output time step
a_bin = 1
size = 20
Rmax = 70
fs = 14 # font size
C0, C1, C2 = 'k', 'b', 'y'
size = size*a_bin

name = 'Kep-47' #name of the simulation

class astrophysical_object():
    def __init__(self, id: int, name: str, shorthand_name: str) -> None:
        self.id = id
        self.name = name
        self.shorthand_name = shorthand_name

class data_id():
    def __init__(self, name: str, out_file: int, legend_name: str) -> None:
        self.name = name
        self.out_file = out_file
        self.legend_name = legend_name

objects_in_kep47 = [
    astrophysical_object(0, 'Star 1', '1'),
    astrophysical_object(1, 'Star 2', '2'),
    astrophysical_object(2, 'kep47b', '47b'),
    astrophysical_object(3, 'kep47c', '47d'),
    astrophysical_object(4, 'kep47d', '47c'),
]