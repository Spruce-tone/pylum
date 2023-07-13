from ctypes import Union
from re import L
import numpy as np
from typing import Union
from scipy.special import jv, j1
from pypkg.utils import *

def tilt(uin: np.ndarray, L: float, wavelength: float, alpha: float, theta: float) -> np.ndarray:
    '''
    tilt function

    Parameters
    ----------
    uin : np.ndarray
        input filed
    L : float
        side length
    wavelength : float
        wavelength
    alpha : float
        tilt angle (in degree)
    theta : float
        rotation angle (in degree) (x axis 0)

    Returns
    -------
    np.ndarray
        out field
    '''    
    M, _ = uin.shape # input field array size
    dx = L/M # sample interval
    k = 2 * np.pi / wavelength # wavenumber

    x = np.linspace(-L/2, L/2-dx, M, endpoint=True) # coordinates
    X, Y = np.meshgrid(x, x)

    alpha_deg = np.deg2rad(alpha)
    theta_deg = np.deg2rad(theta)
    uout = uin * np.exp(1j * k * (X*np.cos(theta_deg) + Y*np.sin(theta_deg)) * np.tan(alpha_deg))
    return uout

def focus(uin: np.ndarray, L: float, wavelength: float, zf: float) -> np.ndarray:
    '''
    transmittance function for focus

    Parameters
    ----------
    uin : np.ndarray
        input filed
    L : float
        side length
    wavelength : float
        wavelength
    zf : float
        focal distance (+ converge, - diverge)

    Returns
    -------
    np.ndarray
        out field
    '''  
    M, _ = uin.shape # get input array size
    dx = L/M    # sample interval
    k = 2 * np.pi / wavelength # wavenumber

    x = np.linspace(-L/2, L/2-dx, M, endpoint=True) # coordinates
    X, Y = np.meshgrid(x, x)

    uout = uin * np.exp(-1j * k / (2*zf) * (X**2 + Y**2))
    return uout

def zoneplate(uin: np.ndarray, L: float, wavelength: float, zf: float, w: float) -> np.ndarray:
    '''
    transmittance function for focus

    Parameters
    ----------
    uin : np.ndarray
        input filed
    L : float
        side length
    wavelength : float
        wavelength
    zf : float
        focal distance (+ converge, - diverge)
    w : float
        radius of zone plate

    Returns
    -------
    np.ndarray
        out field
    '''  
    M, _ = uin.shape # get input array size
    dx = L/M    # sample interval
    k = 2 * np.pi / wavelength # wavenumber

    x = np.linspace(-L/2, L/2-dx, M, endpoint=True) # coordinates
    X, Y = np.meshgrid(x, x)

    uout = uin * 1/2 * (1 + np.cos(k/(2*zf) * (X**2 + Y**2))) * circle_func(X, Y, w)
    return uout

def signzoneplate(uin: np.ndarray, L: float, wavelength: float, zf: float, w: float) -> np.ndarray:
    '''
    transmittance function for focus

    Parameters
    ----------
    uin : np.ndarray
        input filed
    L : float
        side length
    wavelength : float
        wavelength
    zf : float
        focal distance (+ converge, - diverge)
    w : float
        radius of zone plate

    Returns
    -------
    np.ndarray
        out field
    '''  
    M, _ = uin.shape # get input array size
    dx = L/M    # sample interval
    k = 2 * np.pi / wavelength # wavenumber

    x = np.linspace(-L/2, L/2-dx, M, endpoint=True) # coordinates
    X, Y = np.meshgrid(x, x)

    uout = uin * 1/2 * (1 + np.sign(np.cos(k/(2*zf) * (X**2 + Y**2)))) * circle_func(X, Y, w)
    return uout

def cosGrating(X: np.ndarray, Y: np.ndarray, P: float, D: float) -> np.ndarray:
    '''
    transmittance function for cosine magnitude grating

    Parameters
    ----------
    X, Y : np.ndarray
        input coordinates array
    P : float
        period of grating
    D : float
        square aperture of grating 

    Returns
    -------
    np.ndarray
        graing intensity profile
    '''    
    u = 1/2 * (1 - np.cos(2*np.pi*X/P)) * rect(X/D) * rect(Y/D)
    return u

def squareGrating(X: np.ndarray, Y: np.ndarray, P: float, D: float) -> np.ndarray:
    '''
    transmittance function for squre magnitude grating

    Parameters
    ----------
    X, Y : np.ndarray
        input coordinates array
    P : float
        period of grating
    D : float
        square aperture of grating 

    Returns
    -------
    np.ndarray
        graing intensity profile
    '''
    x = X[0, :]
    fc = np.fft.fft(np.fft.fftshift(ucomb(x/P)))
    fr = np.fft.fft(np.fft.fftshift(rect(x/(P/2))))
    ux = np.fft.ifftshift(np.fft.ifft(fc*fr)) # 1D convolution rect and comb function
    u1 = np.tile(ux, (len(X[:, 0]), 1)) # replicate to 2D
    u1 = u1 * rect(X/D) * rect(Y/D) # set aperture size
    return u1

def cosPhaseGrating(X: np.ndarray, Y: np.ndarray, P: float, D: float, m: float) -> np.ndarray:
    '''
    transmittance function for cosine magnitude grating

    Parameters
    ----------
    X, Y : np.ndarray
        input coordinates array
    P : float
        period of grating
    D : float
        square aperture of grating 
    m : float
        period factor

    Returns
    -------
    np.ndarray
        graing intensity profile
    '''    
    u = np.exp(1j * np.pi/m * np.cos(2*np.pi*X/P)) * rect(X/D) * rect(Y/D)
    u = np.where(rect(X/D) * rect(Y/D) > 0, u, 0)
    return u

def trianglePhaseGrating(X: np.ndarray, Y: np.ndarray, P: float, D: float, m: float) -> np.ndarray:
    '''
    transmittance function for squre magnitude grating

    Parameters
    ----------
    X, Y : np.ndarray
        input coordinates array
    P : float
        period of grating
    D : float
        square aperture of grating 
    m : float
        period factor

    Returns
    -------
    np.ndarray
        graing intensity profile
    '''
    x = X[0, :]
    fc = np.fft.fft(np.fft.fftshift(ucomb(x/P)))
    fr = np.fft.fft(np.fft.fftshift(1j * np.pi / m * triangle(x/(P/2))))
    ux = np.fft.ifftshift(np.fft.ifft(fc*fr)) # 1D convolution rect and comb function
    u1 = np.tile(ux, (len(x), 1)) # replicate to 2D
    u1 = np.exp(u1)
    u1 = u1 * rect(X/D) * rect(Y/D) # set aperture size
    u1 = np.where(rect(X/D) * rect(Y/D) > 0, u1, 0)
    return u1

def squarePhaseGrating(X: np.ndarray, Y: np.ndarray, P: float, D: float, m: float) -> np.ndarray:
    '''
    transmittance function for squre magnitude grating

    Parameters
    ----------
    X, Y : np.ndarray
        input coordinates array
    P : float
        period of grating
    D : float
        square aperture of grating 
    m : float
        period factor

    Returns
    -------
    np.ndarray
        graing intensity profile
    '''
    x = X[0, :]
    fc = np.fft.fft(np.fft.fftshift(ucomb(x/P)))
    fr = np.fft.fft(np.fft.fftshift(1j * np.pi / m * rect(x/(P/2))))
    ux = np.fft.ifftshift(np.fft.ifft(fc*fr)) # 1D convolution rect and comb function
    u1 = np.tile(ux, (len(X[:, 0]), 1)) # replicate to 2D
    u1 = np.exp(u1)
    u1 = u1 * rect(X/D) * rect(Y/D) # set aperture size
    u1 = np.where(rect(X/D) * rect(Y/D) > 0, u1, 0)
    return u1

def squareGrating1D(X: np.ndarray, P: float, D: float) -> np.ndarray:
    '''
    transmittance function for squre magnitude grating

    Parameters
    ----------
    X, Y : np.ndarray
        input coordinates array
    P : float
        period of grating
    D : float
        square aperture of grating 

    Returns
    -------
    np.ndarray
        graing intensity profile
    '''
    x = X[0, :]
    fc = np.fft.fft(np.fft.fftshift(ucomb(x/P)))
    fr = np.fft.fft(np.fft.fftshift(rect(x/(P/2))))
    ux = np.fft.ifftshift(np.fft.ifft(fc*fr)) # 1D convolution rect and comb function
    u1 = ux * rect(X/D) # set aperture size
    return u1

def triangleGrating1D(X: np.ndarray, P: float, D: float) -> np.ndarray:
    '''
    transmittance function for squre magnitude grating

    Parameters
    ----------
    X, Y : np.ndarray
        input coordinates array
    P : float
        period of grating
    D : float
        square aperture of grating 

    Returns
    -------
    np.ndarray
        graing intensity profile
    '''
    x = X[0, :]
    fc = np.fft.fft(np.fft.fftshift(ucomb(x/P)))
    fr = np.fft.fft(np.fft.fftshift(triangle(x/(P/2))))
    ux = np.fft.ifftshift(np.fft.ifft(fc*fr)) # 1D convolution rect and comb function
    u1 = ux * rect(x/D) # set aperture size
    return u1

def cylindricalfocus(uin: np.ndarray, L: float, wavelength: float, 
                    zf: float, axis: str='X') -> np.ndarray:      
    '''
    transmittance function for focus

    Parameters
    ----------
    uin : np.ndarray
        input filed
    L : float
        side length
    wavelength : float
        wavelength
    zf : float
        focal distance (+ converge, - diverge)
    axis : str
        determine axis where is focused, by default 'X'
        'X' or 'Y' is possible

    Returns
    -------
    np.ndarray
        out field
    '''  
    M, _ = uin.shape # get input array size
    dx = L/M    # sample interval
    k = 2 * np.pi / wavelength # wavenumber

    x = np.linspace(-L/2, L/2-dx, M, endpoint=True) # coordinates
    X, Y = np.meshgrid(x, x)
    
    if axis=='X':
        uout = uin * np.exp(-1j * k / (2*zf) * (X**2))
    else:
        uout = uin * np.exp(-1j * k / (2*zf) * (Y**2))
    return uout