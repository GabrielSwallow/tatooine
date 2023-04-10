from tools import Unit_conv

all_data_dir = '../Data/'
global_plots_dir = '../Plots/'
cart  = False
logsc = True
nbody = True
var   = "rho" # rho, vx1, vx2
nts   = 100 # output time step
a_bin = 1
size = 15 # in a_b
Rmax = 70
fs = 14 # font size
C0, C1, C2 = 'k', 'b', 'y'
overwrite_plots = True

size = size * Unit_conv.distance(a_bin)




name = 'Kep-47' #name of the simulation

class data_id():
    def __init__(self, name: str, out_file: int, legend_name: str) -> None:
        self.name = name
        self.out_file = out_file
        self.legend_name = legend_name

class astrophysical_object():
    def __init__(self, id: int, name: str, shorthand_name: str, radius: float) -> None:
        self.id = id
        self.name = name
        self.shorthand_name = shorthand_name
        self.radius = radius

# class astrophysical_objects():
#     def __init__(self) -> None:
#         self.Star1 = astrophysical_object(0, 'Star 1', '1', 1.),
#         self.Star2 = astrophysical_object(1, 'Star 2', '2', 1.),
#         self.kep47b = astrophysical_object(2, 'kep47b', '47b', 3.5322),
#         self.kep47c = astrophysical_object(3, 'kep47c', '47d', 8.5844),
#         self.kep47d = astrophysical_object(4, 'kep47d', '47c', 11.8330),
#         self.cavity = astrophysical_object(5, 'Cavity', 'gap', 0)

KepStar1_astrophysical_object = astrophysical_object(0, 'Star 1', '1', 1.)
KepStar2_astrophysical_object = astrophysical_object(1, 'Star 2', '2', 1.)
Kep47b_astrophysical_object = astrophysical_object(2, 'kep47b', '47b', 3.53)
Kep47d_astrophysical_object = astrophysical_object(3, 'kep47d', '47d', 8.5844)
Kep47c_astrophysical_object = astrophysical_object(4, 'kep47c', '47c', 11.8330)
cavity_astrophysical_object = astrophysical_object(5, 'Cavity', 'gap', 0.)
instability_limit_astrophysical_object = astrophysical_object(6, 'instability limit', 'crit', 2.31)
Jupiter_astrophysical_object = astrophysical_object(7, 'Jupiter', 'J', 0.)

objects_in_kep47 = [
    KepStar1_astrophysical_object,
    KepStar2_astrophysical_object,
    Kep47b_astrophysical_object,
    Kep47c_astrophysical_object,
    Kep47d_astrophysical_object,
    cavity_astrophysical_object,
    instability_limit_astrophysical_object,
]

# 3
# 2
# 0
# small
# y
# 2
# 3
# medium
# y
# 2
# 5
# large
# y
# 9
# 6
# none
# asdf
