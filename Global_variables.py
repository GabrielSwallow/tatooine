all_data_dir = '../Data/'
cart  = False
logsc = True
nbody = True
var   = "rho" # rho, vx1, vx2
nts   = 100 # output time step
a_bin = 1
size = 12
Rmax = 70
fs = 14 # font size
C0, C1, C2 = 'k', 'b', 'y'
size = size*a_bin

name = 'Kep-47' #name of the simulation

class astrophysical_object():
    def __init__(self, id: int, name: str) -> None:
        self.id = id
        self.name = name

objects_in_kep47 = [
    astrophysical_object(0, 'Star 1'),
    astrophysical_object(1, 'Star 2'),
    astrophysical_object(2, 'kep47b'),
    astrophysical_object(3, 'kep47c'),
    astrophysical_object(4, 'kep47d'),
]