import os, sys
from os.path import join as opj
import numpy as np
import importlib.util
from pypkg.utils import *

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
    def __init__(self, fname: str):
        super().__init__()
        self.fname = fname

        self.define_vars()
        self.init_fdtd()

    def init_fdtd(self):
        self.addfdtd()
        self.switchtolayout()
        self.deleteall()

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
            if (fname.endswith('.txt')) and (not self.materialexists(fname.split('.')[0])):
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
        self.addplane()
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
                xlb: str='Periodic', xub: str='Periodic',
                ylb: str='Periodic', yub: str='Periodic',
                zlb: str='PML', zub: str='PML',
                x: float=0, y: float=0,
                mesh_type: str='auto non-uniform',
                mesh_acc: int=2,
                dt_stb_factor: float=0.95,
                sim_time: float=1000e-14,
                auto_shut_min: float=0.5e-4):
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

    def add_monitor(self,
                    mname: str, mtype: str,
                    ov_glob_m: int=1, 
                    iswlinspace: int=1, 
                    num_frq_points: int=100,
                    issrclim: int=1, 
                    x: float=0, y: float=0):
        self.addprofile()
        self.set('name', mname) 
        self.set('monitor type', mtype)
        self.set('override global monitor settings', ov_glob_m) 
        self.set('use wavelength spacing', iswlinspace)
        self.set('frequency points', num_frq_points) 
        self.set('use source limits', issrclim)
        self.set('x', x) 
        self.set('y', y) 

    def screening(self, 
                  unit_cell, height, src_name, 
                dia, frq_point, lamb, total, mname):
        ucell_len = len(unit_cell)
        hlen = len(height)

        for idx, ucell in enumerate(unit_cell):
            self.unicell_dependency(ucell, src_name, mname)
            for jdx, H in enumerate(height):
                self.height_dependency(H, mname)

                self.reference_field(struct_name='pillar',
                                     mname=mname,
                                     lamb=lamb,
                                     height=H,
                                     unitcell=ucell)



    def reference_field(self, struct_name: str,
                        mname: str, lamb: np.ndarray,
                        height: float, unitcell: float):
        data = {}

        self.select(struct_name)
        self.set('enabled', False)
        self.run()

        T = self.transmission(mname)
        REF_phase_Ex = np.zeros((len(lamb), 1))
        REF_phase_Ey = np.zeros((len(lamb), 1))
        REF_Ex = np.zeros((len(lamb), 1))
        REF_Ey = np.zeros((len(lamb), 1))
        REF_Ez = np.zeros((len(lamb), 1))

        for ldx, w in enumerate(lamb):
             # For normalization, record farfield eletric field in the case of no structure.
            # E = self.farfieldexact(mname=mname, x=0, y=0, z=1e-4, opt={'f': ldx})
            E = self.farfieldexact(mname, 0, 0, 1e-4, ldx+1)
            REF_Ex[ldx] = E.squeeze()[0]
            REF_Ey[ldx] = E.squeeze()[1]

        data['wavelength'] = lamb
        data['Ex_ref'] = REF_Ex
        data['Ey_ref'] = REF_Ey
        data['transmission'] = T
        data['height'] = height
        data['unit_cell'] = unitcell
         
        SavePickle(data, f'{self.fname}_ref_phase_[unit_cell]_{unitcell/self.nm}_[height]_{height/self.nm}')

        # matlabsave("REF_phase"+"_U="+num2str(unitcell(idx)/nm),lambda,REF_Ex,REF_Ey,T,H,U);	    
            
        self.switchtolayout()    
        self.select(struct_name)
        self.set('enabled', True)
        
    def unicell_dependency(self, unitcell, src_name, mname):
        # unit cell size dependent model setting
        self.setnamed(src_name, 'x span', unitcell*2)
        self.setnamed(src_name, 'y span', unitcell*2)
        # self.setnamed('transmission', 'x span', unitcell*2)
        # self.setnamed('transmission', 'y span', unitcell*2)
        self.setnamed(mname, 'x span', unitcell)
        self.setnamed(mname, 'y span', unitcell)
        self.setnamed('FDTD', 'x span', unitcell)
        self.setnamed('FDTD', 'y span', unitcell)
        self.setnamed('substrate', 'x span', unitcell*2)
        self.setnamed('substrate', 'y span', unitcell*2)

    def height_dependency(self, height, mname):
        # height size dependent model setting
        self.setnamed('FDTD', 'z max', height + 1*self.um)
        self.setnamed('FDTD', 'z min', -1*self.um)
        # self.setnamed('transmission', 'z', height + 200*self.nm)
        self.setnamed(mname, 'z', height + 200*self.nm)
        self.setnamed('pillar', 'z min', 0) 
        self.setnamed('pillar', 'z max', height)
        self.setnamed('substrate', 'z max', 0)
        self.setnamed('substrate', 'z min', -2*self.um)







