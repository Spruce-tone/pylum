from ctypes import Union
from re import L
import numpy as np
from typing import Union, Tuple, List
from scipy.special import jv, j1
import matplotlib.pyplot as plt

def ucomb(x: Union[np.ndarray, int]) -> Union[np.ndarray, int]:
    '''
    unit sample comb function
    sequence of unit values for x=integer value
    round is used to truncate roundoff error

    Parameters
    ----------
    x : Union[np.ndarray, int]
        grid or single value of coordinates

    Returns
    -------
    Union[np.ndarray, int]
        grid or single value of comb function
    '''
    x = np.round(x * 10**6) / 10**6
    x = ((x%1)==0)*1
    return x        

def rect(x: Union[np.ndarray, int]) -> Union[np.ndarray, int]:
    '''
    make rect function
    if input is x/(2w), w is the half-width of the rect function

    Parameters
    ----------
    x : Union[np.ndarray, int]
        grid of coordinates

    Returns
    -------
    Union[np.ndarray, int]
        grid of rect function
    '''    
    return (np.abs(x) < 1/2)*1

def triangle(x: Union[np.ndarray, int]) -> Union[np.ndarray, int]:
    return (np.abs(x) <= 1 ) * (1 - np.abs(x))

def circle_func(x_array: np.ndarray, y_array: np.ndarray, r: Union[int, float]) -> np.ndarray:
    return (np.sqrt(x_array**2 + y_array**2) <= r)*1

def jinc(x: Union[np.ndarray, int], a: float=1.0) -> Union[np.ndarray, int]:
    mask = (x!=0)
    out = np.pi*np.ones_like(x)
    out[mask] = a**2*j1(2*np.pi*a*x[mask])/x[mask]/a
    return out

def FresnelNumber(w: float, wavelength: float, z: float) -> float:
    '''
    Fresnel number (N_F)
    N_F < 1 
        commonly accepted that the observation plane is in Fresnel region, 
        where the Fresnel approximation typically lead to useful results
        For relatively smooth fields over the source aperture, 
        Fresnel approximation can be applicable up to even 20 or 30


    Parameters
    ----------
    w : float
        half-width of a squre aperture or radius of a circular aperture
    wavelength : float
        wavelength
    z : float
        propagation distance

    Returns
    -------
    float
        Fresnel number
    '''    
    return w**2 / (wavelength * z)

def propa_criterion(dx: float, wavelength: float, z: float, L: float):
    '''
    propagation design criterion

    Parameters
    ----------
    dx : float
        sampling interval at source plane
    wavelength : float
        wavelength
    z : float
        propagation distance
    L : float
        array side length
    '''    
    # determine sampling regime
    criterion = np.round(wavelength * z / L * 1e6)/1e6
    norm_cirterion = np.round(criterion / dx * 1e6)/1e6
    if dx > criterion:
        methods = 'TF'
        sampling = 'TF (H) oversampling, IR (h) undersampling'
        B1 = 1 / (2*dx) # source bandwidth limit
    elif dx == criterion:
        methods = 'TF'
        sampling = 'critical sampling'
        B1 = 1 / (2*dx) # source bandwidth limit 
    elif dx < criterion:
        methods = 'IR'
        sampling = 'TF (H) undersampling, IR (h) oversampling'
        B1 = 1 / (2 * criterion) # source bandwidth limit 

    print(f'dx : {dx} | wavelength : {wavelength} | propagation distance : {z} | sidelength : {L}')
    print(f'criterion (norm.) : {norm_cirterion} | cirterion : {criterion} \nmethod : {methods} | {sampling} | B1 : {B1}\n')
    return methods

def GaussianField(X: np.ndarray, Y: np.ndarray, A0: float, w0: float) -> np.ndarray:
    '''
    make Gaussian source field

    Parameters
    ----------
    X : np.ndarray
        X coordinates array
    Y : np.ndarray
        Y coordinates array
    A0 : float
        amplitude
    w0 : float
        source beam e^-2 radius (at z=0)

    Returns
    -------
    np.ndarray
        Gaussian field
    '''    
    U0 = A0 * np.exp(-(X**2 + Y**2) / w0**2)
    return U0

def RayleighRange(w0: float, wavelength: float) -> float:
    '''
    Rayleigh Range

    Parameters
    ----------
    w0 : float
        source beam e^-2 radius (at z=0)
    wavelength : float
        wavelength

    Returns
    -------
    float
        Rayleigh Range (zR)
    '''    
    zR = np.pi * w0**2 / wavelength
    return zR

def GaussianPropRadius(w0: float, z: float, zR: float) -> float:
    '''
    Beam radius of propagating Gaussian beam at distance z 

    Parameters
    ----------
    w0 : float
        source beam e^-2 radius (at z=0)
    z : float
        propagation distance
    zR : float
        Rayleigh range (zR = np.pi * w0**2 / wavelength)

    Returns
    -------
    float
        Beam radius of propagating Gaussian beam at distance z 
    '''    
    wz = w0 * np.sqrt(1 + (z / zR)**2)
    return wz

def LaserPropIntensity(X: np.ndarray, Y: np.ndarray, I0: float, w0: float, wz: float) -> np.ndarray:
    '''
    The irradiance (intensity) distribution of a Gaussian laser beam (TEM 00 mode)
    propagating in the z direction

    Parameters
    ----------
    X : np.ndarray
        X coordinates array
    Y : np.ndarray
        Y coordinates array
    I0 : float
        source irradiance at beam center (x, y = 0)
    w0 : float
        source beam e^-2 radius (at z=0)
    wz : float
        beam radius at distance z

    Returns
    -------
    np.ndarray
        Intensity distribution of Gaussian beam at z
        I(x, y, z)
    '''    
    I = I0 * (w0 / wz)**2 * np.exp(-2 * (X**2 + Y**2) / wz**2)
    return I

def coh_cutoff_frq(wxp: float, wavelength: float, zxp: float) -> float:
    '''
    coherent cutoff frequency

    Parameters
    ----------
    wxp : float
        raidius of exit pupil
    wavelength : float
        wavelength
    zxp : float
        exit pupil distance

    Returns
    -------
    float
        coherent cutoff frequency
    '''    
    return wxp/(wavelength * zxp)

def spatial_grid(L: Union[Tuple[float, float], float], M: Union[Tuple[int, int], int]) -> Tuple[np.array, np.array, np.array, np.array]:
    '''
    make spatial grid

    Parameters
    ----------
    L : Union[Tuple[float, float], float]
        side length 
        (Lx, Ly) or Lx
    M : Union[Tuple[int, int], int]
        No. sample of each side
        (Mx, My) or Mx

    Returns
    -------
    Tuple[np.array, np.array, np.array, np.array]
        length of x_array (y_array) is Mx (My) 
        shape of X_grid and Y_grid is (Mx, My) 
        [x_array, y_array,  X_grid, Y_grid]
    '''    
    assert type(M) in [tuple, int], 'length of M (No. of sample) should be tuple or int'
    
    if type(M)==tuple:
        assert type(L)==tuple, 'type of L (side length) should be tuple'
        assert len(L)==2, 'length of L (side length) should be 2'
        assert len(M)==2, 'length of M (No. of sample) should be 2'
        Lx, Ly = L
        Mx, My = M
        dx = Lx / Mx
        dy = Ly / My
        x = np.linspace(-Lx/2, Lx/2 - dx, Mx, endpoint=True)
        y = np.linspace(-Ly/2, Ly/2 - dy, My, endpoint=True)

    else:
        dx = L / M
        x = np.linspace(-L/2, L/2 - dx, M, endpoint=True)
        y = x

    X, Y = np.meshgrid(x, y)
    return (x, y, X, Y)

def freq_grid(L: Union[Tuple[float, float], float], M: Union[Tuple[int, int], int]) -> Tuple[np.array, np.array, np.array, np.array]:
    '''
    make frequency grid

    Parameters
    ----------
    L : Union[Tuple[float, float], float]
        side length 
        (Lx, Ly) or Lx
    M : Union[Tuple[int, int], int]
        No. sample of each side
        (Mx, My) or Mx

    Returns
    -------
    Tuple[np.array, np.array, np.array, np.array]
        length of fx_array (fy_array) is Mx (My) 
        shape of Fx_grid and Fx_grid is (Mx, My) 
        [fx_array, fy_array,  Fx_grid, Fy_grid]
    '''    
    assert type(M) in [tuple, int], 'length of M (No. of sample) should be tuple or int'
    
    if type(M)==tuple:
        assert type(L)==tuple, 'type of L (side length) should be tuple'
        assert len(L)==2, 'length of L (side length) should be 2'
        assert len(M)==2, 'length of M (No. of sample) should be 2'
        Lx, Ly = L
        Mx, My = M
        dx = Lx / Mx
        dy = Ly / My
        fx = np.linspace(-1/(2*dx), 1/(2*dx) - 1/L, Mx, endpoint=True)
        fy = np.linspace(-1/(2*dy), 1/(2*dy) - 1/L, My, endpoint=True)

    else:
        dx = L / M
        fx = np.linspace(-1/(2*dx), 1/(2*dx) - 1/L, M, endpoint=True)
        fy = fx

    Fx, Fy = np.meshgrid(fx, fy)
    return (fx, fy, Fx, Fy)

def fig_gen(nrows: int=1, ncols: int=1, scale_factor: int=5, **kwargs):
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*scale_factor, nrows*scale_factor))
    
    if nrows==1 and ncols==1:
        return fig, ax
    else:
        ax = ax.ravel()
        return fig, ax
    
def axis_label(ax, **kwargs):
    if kwargs.get('set_xlabel', None) is not None:
        ax.set_xlabel(kwargs.get('set_xlabel', None))
        
    if kwargs.get('set_ylabel', None) is not None:
        ax.set_ylabel(kwargs.get('set_ylabel', None))
        
    if kwargs.get('set_xticks', None) is not None:
        ax.set_xticks(kwargs.get('set_xticks', None))
        
    if kwargs.get('set_yticks', None) is not None:
        ax.set_yticks(kwargs.get('set_yticks', None))
        
    if kwargs.get('set_title', None) is not None:
        ax.set_title(kwargs.get('set_title', None))    
    return ax

