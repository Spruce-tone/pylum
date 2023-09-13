import os, sys
from os.path import join as opj
import numpy as np
import pandas as pd
import importlib.util
from pypkg.utils import *
from typing import Union, List
from tqdm.autonotebook import tqdm

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

M = 1
MM =  1e-3 * M
UM = 1e-6 * M
NM = 1e-9 * M


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
        self.MATLIB_PATH = './material_lib'
        self.light_vel = self.c() # light velocity, 299792458 m/s
        self.reflection_monitor_name = 'reflection'
        self.transmission_monitor_name = 'transmission'
        self.source_name = 'source'
        self.substrate_name = 'substrate'
        self.structure_name = 'structure'

    def set_wavelength(self, lamb_st: float, lamb_ed: float, num_points: int, endpoint: bool=True):
        self.lamb = np.linspace(lamb_st, lamb_ed, num_points, endpoint=endpoint)
        self.num_frq = len(self.lamb)
        self.lamb_idx = np.array(range(self.num_frq))

    def set_height(self, height_list: Union[List, np.ndarray]):
        if type(height_list)==list:  
            self.height = np.array(height_list)
        elif type(height_list)==np.ndarray:
            assert height_list.ndim==1, 'height list should be 1D array'
            self.height = height_list
        self.height_idx = np.array(range(len(self.height)))

    def set_period(self, period_list: Union[List, np.ndarray]):
        if (type(period_list)==list):  
            self.period = np.array(period_list)
        elif type(period_list)==np.ndarray:
            assert period_list.ndim==1, 'height list should be 1D array'
            self.period = period_list
        self.period_idx = np.array(range(len(self.period)))

    def add_substrate(self, 
                 material: str='SiO2 (Glass) - Palik', 
                 x: float=0, y: float=0):
        self.addrect()
        self.set('name', self.substrate_name)
        self.set('x', x)
        self.set('y', y)
        self.set('material', material)

    def add_struct(self, struct_type: str='circle',
                   lattice: str='square', 
                   x: float=0, y: float=0, 
                   material: str='SiO2 (Glass) - Palik'):
        
        if struct_type=='circle':
            self.struct_type = struct_type

            if lattice=='square':
                self.lattice = lattice
                self.addcircle()
                self.set('name', self.structure_name)
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
        complex_sampled[:, 0] = self.light_vel / (lamb * NM)
        complex_sampled[:, 1] = (n**2 - k**2) + 1j*(2*n*k) # real and imaginary part of permittivity
        
        if not self.materialexists(mat_name):
            Mat_add = self.addmaterial('Sampled 3D data')
            self.setmaterial(Mat_add, 'name', mat_name)
        self.setmaterial(mat_name, 'sampled data', complex_sampled)
            
    def add_planewave(self, 
                      inj_axis: str='z-axis',
                      direc: str='Forward',
                      pol_angle: float=0,
                      phase: float=0,
                      x: float=0, y: float=0, z: float=-200*NM,
                      lambda_st: float=400e-9, 
                      lambda_ed: float=800e-9):
        self.addplane()
        self.set('name', self.source_name) 
        self.set('injection axis', inj_axis)
        self.set('direction', direc) 
        self.set('polarization angle', pol_angle)
        self.set('phase', phase)
        self.set('x', x) 
        self.set('y', y) 
        self.set('z', z) 
        self.set('wavelength start', self.lamb[0])
        self.set('wavelength stop', self.lamb[-1])

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

    def add_transmission_monitor(self,
                    mtype: str='2D Z-normal',
                    ov_glob_m: int=1, 
                    iswlinspace: int=1, 
                    issrclim: int=1, 
                    x: float=0, y: float=0):
        self.addprofile()
        self.set('name', self.transmission_monitor_name) 
        self.set('monitor type', mtype)
        self.set('override global monitor settings', ov_glob_m) 
        self.set('use wavelength spacing', iswlinspace)
        self.set('frequency points', self.num_frq) 
        self.set('use source limits', issrclim)
        self.set('x', x) 
        self.set('y', y) 

    def add_reflection_monitor(self,
                    mtype: str='2D Z-normal',
                    ov_glob_m: int=1, 
                    iswlinspace: int=1, 
                    issrclim: int=1, 
                    x: float=0, y: float=0):
        self.addprofile()
        self.set('name', self.reflection_monitor_name) 
        self.set('monitor type', mtype)
        self.set('override global monitor settings', ov_glob_m) 
        self.set('use wavelength spacing', iswlinspace)
        self.set('frequency points', self.num_frq) 
        self.set('use source limits', issrclim)
        self.set('x', x) 
        self.set('y', y) 

    def simul_span_params(self, FDTD_height_offset: float=None,
                           source_depth: float=None,
                           substrate_depth: float=None,
                           T_monitor_height: float=None,
                           R_monitor_depth: float=None,
                           far_field_height: float=100*UM):
        
        if FDTD_height_offset==None:
            self.FDTD_height_offset = 1*UM
        else:
            self.FDTD_height_offset = FDTD_height_offset

        if source_depth==None:
            self.source_depth = -0.2 * self.FDTD_height_offset
        else:
            self.source_depth = source_depth

        if substrate_depth==None:
            self.substrate_depth = -2 * self.FDTD_height_offset
        else:
            self.substrate_depth = substrate_depth

        if T_monitor_height==None:
            self.T_monitor_height = 0.2 * self.FDTD_height_offset
        else:
            self.T_monitor_height = T_monitor_height

        if R_monitor_depth==None:
            self.R_monitor_depth = -0.4 * self.FDTD_height_offset
        else:
            self.R_monitor_depth = R_monitor_depth

        self.far_field_height = far_field_height

    def simul_span_setting(self, period: float, height: float):
        # Source
        self.setnamed(self.source_name, 'x span', period*2)
        self.setnamed(self.source_name, 'y span', period*2)
        self.setnamed(self.source_name, 'z', self.source_depth)

        # Transmission monitor
        self.setnamed(self.transmission_monitor_name, 'x span', period) # CLEO used period*2
        self.setnamed(self.transmission_monitor_name, 'y span', period)
        self.setnamed(self.transmission_monitor_name, 'z', height + self.T_monitor_height)

        # Reflection monitor
        self.setnamed(self.reflection_monitor_name, 'x span', period) # CLEO used period*2
        self.setnamed(self.reflection_monitor_name, 'y span', period)
        self.setnamed(self.reflection_monitor_name, 'z', self.R_monitor_depth)

        # FDTD grid setting
        self.setnamed('FDTD', 'x span', period)
        self.setnamed('FDTD', 'y span', period)
        self.setnamed('FDTD', 'z max', height + self.FDTD_height_offset)
        self.setnamed('FDTD', 'z min', -self.FDTD_height_offset)


    def set_sweep_params(self, data: dict):
        if self.struct_type=='circle':
            if self.lattice=='square':
                data['diameter'] = [(x, idx) for idx, x in enumerate(data['diameter'])]
                data['period'] = [(x, idx) for idx, x in enumerate(data['period'])]
                data['height'] = [(x, idx) for idx, x in enumerate(data['height'])]

                df = pd.DataFrame(dict_product(data))

                df['diameter_idx'] = df['diameter'].map(lambda x: x[1]) 
                df['diameter'] = df['diameter'].map(lambda x: x[0]) 

                df['period_idx'] = df['period'].map(lambda x: x[1]) 
                df['period'] = df['period'].map(lambda x: x[0]) 

                df['height_idx'] = df['height'].map(lambda x: x[1]) 
                df['height'] = df['height'].map(lambda x: x[0]) 
                
                # physical size limit
                df = df.loc[df.diameter <= df.period, :]
                df.drop(['material'], axis=1, inplace=True)

                self.sweep_df = df
                
        self.structure_material = [struct_mat[0] for struct_mat in data['material']]
        self.substrate_material = [struct_mat[1] for struct_mat in data['material']]
        self.material_idx = list(range(len(data['material'])))

    def set_sweep_structure(self, struct_mat, sub_mat, sweep_params):
        if self.struct_type=='circle':
            if self.lattice=='square':
                self.setnamed(self.structure_name, 'material', struct_mat)
                self.setnamed(self.structure_name, 'radius', sweep_params['diameter']/2)
                self.setnamed(self.structure_name, 'z min', 0) 
                self.setnamed(self.structure_name, 'z max', sweep_params['height'])

                # Substrate setting 
                self.setnamed(self.substrate_name, 'material', sub_mat)
                self.setnamed(self.substrate_name, 'x span', sweep_params['period']*2)
                self.setnamed(self.substrate_name, 'y span', sweep_params['period']*2)
                self.setnamed(self.substrate_name, 'z max', 0)
                self.setnamed(self.substrate_name, 'z min', self.substrate_depth)

    def screening(self, save_path: str=None):
        if save_path==None:
            save_path = f'./{self.fname}.fsp'
        self.save(save_path)

        for mdx, struct_mat, sub_mat in tqdm(zip(self.material_idx, self.structure_material, self.substrate_material)):
            ref_data = []
            sim_data = []

            for params in tqdm(self.sweep_df.to_dict('records')):
                self.simul_span_setting(period=params['period'], height=params['height'])
                self.set_sweep_structure(struct_mat=struct_mat,
                                         sub_mat=sub_mat,
                                         sweep_params=params)
                
                ref_data, REF_Ex, REF_Ey = self.reference_field(ref_data=ref_data,
                                                                sweep_params=params)
                sim_data = self.field_simulation(sim_data=sim_data, sweep_params=params,
                                                 REF_Ex=REF_Ex, REF_Ey=REF_Ey)
                
                SavePickle(pd.DataFrame(ref_data), f'{self.fname}_ref_phase_material_{mdx:02d}_{struct_mat}_{sub_mat}.pickle')
                SavePickle(pd.DataFrame(sim_data), f'{self.fname}_phase_material_{mdx:02d}_{struct_mat}_{sub_mat}.pickle')


    def field_simulation(self, REF_Ex: np.ndarray, REF_Ey: np.ndarray,
                         sweep_params: dict, sim_data: List):
        self.run()

        Transmission = self.transmission(self.transmission_monitor_name)
        Reflection = self.transmission(self.reflection_monitor_name)

        if self.num_frq==1:
            Transmission = np.array([Transmission])
            Reflection = np.array([Reflection])
        elif self.num_frq > 1:
            Transmission = Transmission.squeeze()
            Reflection = Transmission.squeeze()  

        nEx = np.zeros((self.num_frq,), dtype=np.complex128)
        nEy = np.zeros((self.num_frq,), dtype=np.complex128)
        Ex = np.zeros((self.num_frq,), dtype=np.complex128)
        Ey = np.zeros((self.num_frq,), dtype=np.complex128)
        Ez = np.zeros((self.num_frq,), dtype=np.complex128)
    
        for ldx, W, T, R in zip(self.lamb_idx, self.lamb, Transmission, Reflection):
            # For normalization, record farfield eletric field in the case of no structure.
            E = self.farfieldexact(self.transmission_monitor_name, 0, 0, self.far_field_height, ldx+1)
            Ex[ldx] = E.squeeze()[0]
            Ey[ldx] = E.squeeze()[1]
            nEx[ldx] = Ex[ldx] / REF_Ex[ldx] # normalized (to the case of only substrate) electric field
            nEy[ldx] = Ey[ldx] / REF_Ey[ldx]

            update_params = sweep_params.copy()

            update_params.update(
                        {'wavelength_idx' : ldx,
                        'wavelength' : W,
                        'reflection' : -R,
                        'transmission' : T,
                        'Ex' : E.squeeze()[0],
                        'nEx' : nEx[ldx],
                        'Ey' : E.squeeze()[1],
                        'nEy' : nEy[ldx]})
            sim_data.append(update_params)

        self.switchtolayout()
        return sim_data

    def reference_field(self, ref_data: List,
                        sweep_params: dict):
        self.select(self.structure_name)
        self.set('enabled', False)
        self.run()

        Transmission = self.transmission(self.transmission_monitor_name)
        Reflection = self.transmission(self.reflection_monitor_name)

        if self.num_frq==1:
            Transmission = np.array([Transmission])
            Reflection = np.array([Reflection])
        elif self.num_frq > 1:
            Transmission = Transmission.squeeze()
            Reflection = Transmission.squeeze()          

        REF_phase_Ex = np.zeros((self.num_frq,), dtype=np.complex128)
        REF_phase_Ey = np.zeros((self.num_frq,), dtype=np.complex128)
        REF_Ex = np.zeros((self.num_frq,), dtype=np.complex128)
        REF_Ey = np.zeros((self.num_frq,), dtype=np.complex128)
        REF_Ez = np.zeros((self.num_frq,), dtype=np.complex128)

        for ldx, W, T, R in zip(self.lamb_idx, self.lamb, Transmission, Reflection):
            # For normalization, record farfield eletric field in the case of no structure.
            E = self.farfieldexact(self.transmission_monitor_name, 0, 0, self.far_field_height, ldx+1)
            REF_Ex[ldx] = E.squeeze()[0]
            REF_Ey[ldx] = E.squeeze()[1]

            update_params = sweep_params.copy()

            update_params.update(
                        {'wavelength_idx' : ldx,
                        'wavelength' : W,
                        'reflection' : -R,
                        'transmission' : T,
                        'Ex_ref' : E.squeeze()[0],
                        'Ey_ref' : E.squeeze()[1]})
            ref_data.append(update_params)

        self.switchtolayout()    
        self.select(self.structure_name)
        self.set('enabled', True)
        return ref_data, REF_Ex, REF_Ey



