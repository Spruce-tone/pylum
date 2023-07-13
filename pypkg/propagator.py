from codecs import latin_1_decode
from typing import List, Tuple
import numpy as np

def RSpropTF(u1: np.ndarray, L: float, wavelength: float, z: float) -> np.ndarray:
    '''
    Rayleigh-Sommerfeld propagation - transfer function approach
    assume samse x, y side lengths and uniform sampling

    Parameters
    ----------
    u1 : np.ndarray
        source field plane
    L : float
        source and observation plane side length
    wavelength : float
        wavelength
    z : float
        propagation distance

    Returns
    -------
    np.ndarray
        observation plane field
    '''    
    M, N = u1.shape                             # get input field array size
    dx = L/M                                    # sample interval
    k = 2*np.pi / wavelength                    # wavenumber

    fx = np.linspace(-1/(2*dx), 1/(2*dx) - 1/L, M, endpoint=True) # frequency coordinates
    Fx, Fy = np.meshgrid(fx, fx)

    H = np.exp(1j*k*z*np.sqrt(1-(wavelength*Fx)**2 - (wavelength*Fy)**2))
    H = np.fft.fftshift(H)
    U1 = np.fft.fft2(np.fft.fftshift(u1))
    U2 = U1 * H
    u2 = np.fft.ifftshift(np.fft.ifft2(U2))
    return u2    

def RSpropIR(u1: np.ndarray, L: float, wavelength: float, z: float) -> np.ndarray:
    '''
    Rayleigh-Sommerfeld propagation - impulse response approach
    assume samse x, y side lengths and uniform sampling

    Parameters
    ----------
    u1 : np.ndarray
        source field plane
    L : float
        source and observation plane side length
    wavelength : float
        wavelength
    z : float
        propagation distance

    Returns
    -------
    np.ndarray
        observation plane field
    '''    
    M, N = u1.shape                             # get input field array size
    dx = L/M                                    # sample interval
    k = 2*np.pi / wavelength                    # wavenumber

    x = np.linspace(-L/2, L/2-dx, M, endpoint=True) # spatial coordinates
    X, Y = np.meshgrid(x, x)

    r = np.sqrt(z**2 + X**2 + Y**2)
    h = z/(1j*wavelength) * np.exp(1j*k*r) / r**2 # impulse response
    H = np.fft.fft2(np.fft.fftshift(h)) * dx**2
    U1 = np.fft.fft2(np.fft.fftshift(u1))
    U2 = U1 * H
    u2 = np.fft.ifftshift(np.fft.ifft2(U2))
    return u2

def FresnelpropTF(u1: np.ndarray, L: float, wavelength: float, z: float) -> np.ndarray:
    '''
    propagation - transfer function approach
    assume samse x, y side lengths and uniform sampling

    Parameters
    ----------
    u1 : np.ndarray
        source field plane
    L : float
        source and observation plane side length
    wavelength : float
        wavelength
    z : float
        propagation distance

    Returns
    -------
    np.ndarray
        observation plane field
    '''    
    M, N = u1.shape                             # get input field array size
    dx = L/M                                    # sample interval
    k = 2*np.pi / wavelength                    # wavenumber

    fx = np.linspace(-1/(2*dx), 1/(2*dx) - 1/L, M, endpoint=True) # frequency coordinates
    Fx, Fy = np.meshgrid(fx, fx)

    H = np.exp(-1j*np.pi*wavelength*z*(Fx**2 + Fy**2))*np.exp(1j*k*z) # transfer function
    H = np.fft.fftshift(H)
    U1 = np.fft.fft2(np.fft.fftshift(u1))
    U2 = U1 * H
    u2 = np.fft.ifftshift(np.fft.ifft2(U2))
    return u2

def FresnelpropIR(u1: np.ndarray, L: float, wavelength: float, z: float) -> np.ndarray:
    '''
    propagation - impulse response approach
    assume samse x, y side lengths and uniform sampling

    Parameters
    ----------
    u1 : np.ndarray
        source field plane
    L : float
        source and observation plane side length
    wavelength : float
        wavelength
    z : float
        propagation distance

    Returns
    -------
    np.ndarray
        observation plane field
    '''    
    M, N = u1.shape                             # get input field array size
    dx = L/M                                    # sample interval
    k = 2*np.pi / wavelength                    # wavenumber

    x = np.linspace(-L/2, L/2-dx, M, endpoint=True) # spatial coordinates
    X, Y = np.meshgrid(x, x)
    
    h = 1 / (1j*wavelength*z) * np.exp(1j*k/(2*z) * (X**2 + Y**2)) * np.exp(1j*k*z) # impulse response
    H = np.fft.fft2(np.fft.fftshift(h)) * dx**2
    U1 = np.fft.fft2(np.fft.fftshift(u1))
    U2 = U1 * H
    u2 = np.fft.ifftshift(np.fft.ifft2(U2))
    return u2

def FraunhoferProp(u1: np.ndarray, L1: float, wavelength: float, z: float) -> Tuple[np.ndarray, float]:
    '''
    Fraunhofer propagation
    assume uniform sampling

    Parameters
    ----------
    u1 : np.ndarray
        source plane field
    L1 : float
        source plane side length
    wavelength : float
        wavelength
    z : float
        propagation distance

    Returns
    -------
    List[np.ndarray, float]
        u2 : np.ndarray
            observation plane field
        L2 : float
            observation plane side length

    '''    
    M, N = u1.shape                             # get input field array size
    dx1 = L1/M                                   # sample interval
    k = 2*np.pi / wavelength                    # wavenumber

    L2 = wavelength*z / dx1                     # observation plane sidelength
    dx2 = wavelength*z / L1         # observation plane sample interval
    x2 = np.linspace(-L2/2, L2/2-dx2, num=M, endpoint=True) # observation plane coordinates

    X2, Y2 = np.meshgrid(x2, x2)

    c = np.exp(1j*k*z) / (1j*wavelength*z) * np.exp(1j*k/(2*z) * (X2**2 + Y2**2))
    u2 = c * np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(u1))) * dx1**2

    return (u2, L2)

def FraunhoferProp1D(u1: np.ndarray, L1: float, wavelength: float, z: float) -> Tuple[np.ndarray, float]:
    '''
    Fraunhofer propagation
    assume uniform sampling

    Parameters
    ----------
    u1 : np.ndarray
        source plane field [1 x N] 1D
    L1 : float
        source plane side length [1 x N] 1D
    wavelength : float
        wavelength
    z : float
        propagation distance

    Returns
    -------
    List[np.ndarray, float]
        u2 : np.ndarray
            observation plane field [1 x N] 1D
        L2 : float
            observation plane side length

    '''    
    M = u1.shape[0]                             # get input field array size
    dx1 = L1/M                                   # sample interval
    k = 2*np.pi / wavelength                    # wavenumber

    L2 = wavelength*z / dx1                     # observation plane sidelength
    dx2 = wavelength*z / L1         # observation plane sample interval
    x2 = np.linspace(-L2/2, L2/2-dx2, num=M, endpoint=True) # observation plane coordinates

    # X2, Y2 = np.meshgrid(x2, x2)

    c = np.exp(1j*k*z) / (1j*wavelength*z) * np.exp(1j*k/(2*z) * (x2**2))
    u2 =  np.sqrt(c)* np.fft.ifftshift(np.fft.fft(np.fft.fftshift(u1))) * dx1
    return (u2, L2, dx2, x2)


def prop2step(u1: np.ndarray, L1: float, L2: float, wavelength: float, z: float) -> np.ndarray:
    '''
    2 strep Fresnel diffraction method
    assume uniform sampling and square array

    Parameters
    ----------
    u1 : np.ndarray
        complex field at source plane
    L1 : float
        source plane side-length
    L2 : float
        observation plane side-length
    wavelength : float
        wavelength
    z : float
        propagation distance

    Returns
    -------
    np.ndarray
        output field at observation plane
    '''    
    M, N = u1.shape                             # get input field array size
    k = 2*np.pi / wavelength                    # wavenumber 

    # source plane
    dx1 = L1/M                                   # sample interval
    x1 = np.linspace(-L1/2, L1/2-dx1, M, endpoint=True) # spatial coordinates
    X, Y = np.meshgrid(x1, x1)
    u = u1 * np.exp(1j*k / (2*z*L1) * (L1-L2) * (X**2 + Y**2))
    u = np.fft.fft2(np.fft.fftshift(u))

    # dummy (frequency) plane
    fx1 = np.linspace(-1/(2*dx1), 1/(2*dx1) - 1/L1, num=M, endpoint=True)
    fx1 = np.fft.fftshift(fx1)
    FX1, FY1 = np.meshgrid(fx1, fx1)
    u = np.exp(-1j*np.pi*wavelength*z * L1/L2 * (FX1**2 + FY1**2)) * u
    u = np.fft.ifftshift(np.fft.ifft2(u))

    # observation plane
    dx2 = L2/M
    x2 = np.linspace(-L2/2, L2/2 - dx2, num=M, endpoint=True)
    X, Y = np.meshgrid(x2, x2)

    u2 = L2/L1 * np.exp(1j*k*z) * np.exp(-1j*k/(2*z*L2) * (L1-L2) * (X**2 + Y**2)) * u
    u2 = u2 * (dx1**2/dx2**2) # x1 to x2 scale adjustment
    return u2, x2, dx2
