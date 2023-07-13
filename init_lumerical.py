import os, sys
import numpy as np
import importlib.util


#default path for current release 
sys.path.append('C:\\Program Files\\Lumerical\\FDTD\\api\\python') 
sys.path.append(os.path.dirname(__file__)) #Current directory
os.add_dll_directory('C:/Program Files/Lumerical/FDTD/api/python')


#default path for current release 
spec_win = importlib.util.spec_from_file_location('lumapi', 'C:\\Program Files\\Lumerical\\FDTD\\api\\python\\lumapi.py')

#Functions that perform the actual loading
lumapi = importlib.util.module_from_spec (spec_win) #windows
spec_win.loader.exec_module(lumapi)



if __name__=='__main__':
    fdtd = lumapi.FDTD()
    # nm = 1e-9
    # um = 1e-9
    # lamb = np.linspace(-450*nm, 450*nm, 5*nm)
    # freq_point = len(lamb)

    # height = 600*nm
    # unitcell = 200*nm
    # dia = np.linspace(50*nm, 190*nm, 10*nm)

    # mname = 'Transmission' # monitor name
    # total = len(height) * len(unitcell) * len(dia)

    # print(fdtd.materialexists('TiO2'))