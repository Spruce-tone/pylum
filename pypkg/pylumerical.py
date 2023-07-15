import os, sys
from os.path import join as opj
import numpy as np
import importlib.util

#default path for current release 
PYAPI_PATH = 'C:/Program Files/Lumerical/v221/api/python'
sys.path.append('C:/Program Files/Lumerical/v221/api/python') 
sys.path.append(os.path.abspath('')) #Current directory
# os.add_dll_directory('C:/Program Files/Lumerical/FDTD/api/python')

#default path for current release 
spec_win = importlib.util.spec_from_file_location('lumapi', os.path.join(PYAPI_PATH, 'lumapi.py'))

#Functions that perform the actual loading
lumapi = importlib.util.module_from_spec(spec_win) #windows
spec_win.loader.exec_module(lumapi)



class Pylum(lumapi.FDTD):
    def __init__(self):
        super().__init__()

        self.define_vars()

        self.init_fdtd()

    def init_fdtd(self):
        self.addfdtd()
        self.switchtolayout()
        self.deleteall()
        self.update_material_lib()
        self.set_FDTD()
        self.add_planewave()

    def define_vars(self):
        self.m = 1
        self.nm = 1e-9 * self.m
        self.um = 1e-6 * self.m

        self.MATLIB_PATH = './material_lib'
        self.light_vel = self.c() # light velocity, 299792458 m/s

    def set_wavelength(self):
        pass

    def add_rect(self, name: str='substrate',
                 material: str='SiO2 (Glass) - Palik', 
                 x: float=0, y: float=0):
        self.addrect()
        self.set('name', name)
        self.set('x', x)
        self.set('y', y)
        self.set('material', material)

    def add_circle(self, name: str='pillar', 
                   x: float=0, y: float=0, 
                   material: str='SiO2 (Glass) - Palik'):
        self.addcircle()
        self.set('name', name)
        self.set('x', x)
        self.set('y', y)
        self.set('material', material)

    def update_material_lib(self):
        material_list = os.listdir(self.MATLIB_PATH)

        for fname in material_list:
            if fname.endswith('.txt'):
                self.load_material(fname)

    def load_material(self, fname: str):
        mat_name = fname.split('.')[0]
        mat_path = opj(self.MATLIB_PATH, fname)
        sampled = self.readdata(mat_path)

        lamb = sampled[:, 0] # wavelength
        n = sampled[:, 1] # refractive index
        k = sampled[:, 2] # extinction coefficient
        complex_sampled = np.zeros((len(sampled), 2), dtype=np.complex128)
        
        # wavelength to frequency conversion
        # c = lambda * nu (light velocity = wavelength * frequency)
        complex_sampled[:, 0] = self.light_vel / (lamb * self.nm)
        complex_sampled[:, 1] = (n**2 - k**2) + 1j*(2*n*k) # real and imaginary part of permittivity
        
        if not self.materialexists(mat_name):
            Mat_add = self.addmaterial('Sampled 3D data')
            self.setmaterial(Mat_add, 'name', mat_name)
        self.setmaterial(mat_name, 'sampled data', complex_sampled)
            

    def add_planewave(self, sname: str='source',
                      inj_axis: str='z-axis',
                      direc: str='Forward',
                      pol_angle: float=0,
                      phase: float=0,
                      x: float=0, y: float=0, z: float=-200e-9,
                      lambda_st: float=400e-9, 
                      lambda_ed: float=800e-9):
        self.addplane
        self.set('name', sname) 
        self.set('injection axis', inj_axis)
        self.set('direction', direc) 
        self.set('polarization angle', pol_angle)
        self.set('phase', phase)
        self.set('x', x) 
        self.set('y', y) 
        self.set('z', z) 
        self.set('wavelength start', lambda_st)
        self.set('wavelength stop', lambda_ed)

    def set_FDTD(self,
                dim: str='3D',
                xlb='Periodic', xub='Periodic',
                ylb='Periodic', yub='Periodic',
                zlb='PML', zub='PML',
                x=0, y=0,
                mesh_type='auto non-uniform',
                mesh_acc=2,
                dt_stb_factor=0.95,
                sim_time=1000e-14,
                auto_shut_min=0.5e-4):
        self.addfdtd()
        self.select('FDTD')
        self.set('dimension', dim)
        self.set('x min bc', xlb)
        self.set('x max bc', xub)
        self.set('y min bc', ylb)
        self.set('y max bc', yub)
        self.set('z min bc', zlb)
        self.set('z max bc', zub)
        self.set('x', x) 
        self.set('y', y) 
        self.set('mesh type', mesh_type)
        self.set('mesh accuracy', mesh_acc)
        self.set('dt stability factor', dt_stb_factor)
        self.set('simulation time', sim_time)
        self.set('auto shutoff min', auto_shut_min)

# function add_monitor(mname, mtype,
#     ov_glob_m, iswlinspace, frq_list,
#     issrclim, x, y) {
#     addprofile
#     set('name', mname) 
#     set('monitor type', mtype)
#     set('override global monitor settings', ov_glob_m) 
#     set('use wavelength spacing', iswlinspace)
#     set('frequency points', frq_list) 
#     set('use source limits', issrclim)
#     set('x', x) 
#     set('y', y) 
# }



# function add_rect(name, x, y, material) {
#     addrect
#     set('name', name)
#     set('x', x) 
#     set('y', y) 
#     set('material', material)   
# }

# function add_circ(name, x, y, material){
#     addcircle
#     set('name', name)
#     set('x', x) 
#     set('y', y) 
#     set('material', material)      
# }


# function unitcell_dependency(unitcell, src_name){
#     # unit cell size dependent model setting
#     setnamed(src_name, 'x span', unitcell*2)
#     setnamed(src_name, 'y span', unitcell*2)
#     setnamed('transmission', 'x span', unitcell*2)
#     setnamed('transmission', 'y span', unitcell*2)
#     setnamed('FDTD', 'x span', unitcell)
#     setnamed('FDTD', 'y span', unitcell)
#     setnamed('substrate', 'x span', unitcell*2)
#     setnamed('substrate', 'y span', unitcell*2)
# }


# function height_dependency(height){
#     # height size dependent model setting
#     um = 1e-6
#     nm = 1e-9
    
#     setnamed('FDTD', 'z max', height + 1*um)
#     setnamed('FDTD', 'z min', -1*um)
#     setnamed('transmission', 'z max', height + 200*nm)
    
#     setnamed('pillar', 'z min', 0) 
#     setnamed('pillar', 'z max', height)
    
#     setnamed('substrate', 'z max', 0)
#     setnamed('substrate', 'z min', -2*um)

# }



# function screening(unit_cell, height, src_name, 
# dia, frq_point, lambda, total) {
#     ucell_len = length(unit_cell)
#     hlen = length(height)
#     nm = 1e-9
#     um = 1e-6
    
#     for (idx=1:ucell_len){
#         U = unit_cell(idx)
#         unitcell_dependency(unitcell=U, src_name=src_name)
        
#         for (jdx=1:hlen) {
#             H = height(jdx)
#             height_dependency(H)
            
#             ?'Running reference'
#             select('pillar')
#             set('enabled', 0)
#             run
           
#             T = transmission('transmission') # get transmission spectra
#             REF_phase_Ex = matrix(frq_point,1)
#             REF_phase_Ey = matrix(frq_point,1)
#             REF_Ex = matrix(frq_point,1)
#             REF_Ey = matrix(frq_point,1)
#             REF_Ez = matrix(frq_point,1)
            
#             for (i = 1:frq_point) {
#                 E = farfieldexact('transmission',0,0,1e-4,i) # For normalization, record farfield eletric field in the case of no structure.
#                 REF_Ex(i) = E(1) REF_Ey(i) = E(2)  
#             }
#             matlabsave("REF_phase"+"_U="+num2str(unitcell(idx)/nm),lambda,REF_Ex,REF_Ey,T,H,U)	    
            
#             switchtolayout    
            
#             select('pillar')
#             set('enabled', 1)
            
#             # PARAMETERS SWEEP
#             # ================================================================================
#             counter = 1
            
#             for (w = 1:length(dia)) { 
#                 wid = dia(w)
                        
#                     setnamed('pillar', "radius", wid/2)            
                    
#                     ?"Running "+num2str(counter)+" out of "+num2str(total)  
#                     ?"Diameter: "+num2str(wid/nm)                    
#                     run
                    
                    
#                     T = transmission('transmission') # get transmission spectra
#                     nEx= matrix(frq_point,1)
#                     nEy = matrix(frq_point,1)
#                     Ex = matrix(frq_point,1)
#                     Ey = matrix(frq_point,1)
#                     Ez = matrix(frq_point,1)
                    
#                     for (i = 1:frq_point) {
#                         E = farfieldexact('transmission',0,0,1e-4,i)
#                         Ex(i) = E(1) Ey(i) = E(2)
#                         nEx(i) = Ex(i)/REF_Ex(i)  # normalized (to the case of only substrate) electric field
#                         nEy(i) = Ey(i)/REF_Ex(i)
#                     }
                    
#                     matlabsave("phases_"+"H="+num2str(H/nm)+"_U="+num2str(U/nm)+"_"+num2str(counter),lambda, Ex, Ey,nEx, nEy, T,H,U,wid)	    
#                     switchtolayout
#                     counter = counter + 1
#             }
#         }
#     }

# }



# # call matlab function to count required time for this simulation
# #matlab('simulation_time=toc')
# #matlabget(simulation_time)
# #?'Total elapsed time: ' + num2str(simulation_time) + 's'











    # lamb = np.linspace(-450*nm, 450*nm, 100)
    # freq_point = len(lamb)

    # height = 600*nm
    # unitcell = 200*nm
    # dia = np.linspace(50*nm, 190*nm, 5)

    # mname = 'Transmission' # monitor name
    # total = len(height) * len(unitcell) * len(dia)

    # print(fdtd.materialexists('TiO2'))