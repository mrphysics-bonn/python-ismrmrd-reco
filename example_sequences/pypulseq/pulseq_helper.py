"""
helper functions for spiral pulseq sequence    
"""   

import numpy as np
import matplotlib.pyplot as plt
import math

dt_grad  = 10e-6     # gradient raster [s]
fw_shift = 3.3e-6 # unsigned fat water shift [ppm]

#############
# FFTs
#############

def ifft(sig, dim=None):
    """ Computes the Fourier transform from k-space to image space 
    along a given or all dimensions

    :param img: image space data
    :param dim: vector of dimensions to transform
    :returns: data in k-space (along transformed dimensions)
    """
    import collections.abc
    if dim is None:
        dim = range(sig.ndim)
    elif not isinstance(dim, collections.abc.Iterable):
        dim = [dim]

    sig = np.fft.ifftshift(sig, axes=dim)
    sig = np.fft.ifftn(sig, axes=dim)
    sig = np.fft.fftshift(sig, axes=dim)

    return sig

def fft(sig, dim=None):
    """ Computes the Fourier transform from image space to k-space
    along a given or all dimensions

    :param img: image space data
    :param dim: vector of dimensions to transform
    :returns: data in k-space (along transformed dimensions)
    """
    import collections.abc
    if dim is None:
        dim = range(sig.ndim)
    elif not isinstance(dim, collections.abc.Iterable):
        dim = [dim]

    sig = np.fft.ifftshift(sig, axes=dim)
    sig = np.fft.fftn(sig, axes=dim)
    sig = np.fft.fftshift(sig, axes=dim)

    return sig

################
# Functions for gradient calculations
################

def merge_ramps(grads: list, system=None):
    """
    Merge the ramps of trapezoidal gradients

    grads: List of trapezoidal gradients in merging order. Ramps should have equal length.
    """

    delay = grads[0].delay # initial delay of first gradient
    for grad in grads:
        grad.delay = delay
        delay += grad.rise_time + grad.flat_time 
        
    return add_gradients(grads, system)

def add_gradients(grads: list, system=None):
    """
    Adds up gradient events
    """
    from pypulseq.make_arbitrary_grad import make_arbitrary_grad

    if system==None:
        raise ValueError('Provide the MR System limits.')

    # First gradient defines channel
    channel = grads[0].channel

    # read in gradient waveforms
    grad_length = []
    grad_list = []
    for grad in grads:
        w = waveform_from_seqblock(grad, system)
        grad_list.append(w)
        grad_length.append(len(w))
    
    # prolong waveforms to maximum length and add them up
    length=max(grad_length)
    added_grad = np.array(np.zeros(length))
    for grad in grad_list:
        grad = np.append(grad, np.zeros(length-len(grad)))
        added_grad += grad

    return make_arbitrary_grad(channel=channel, waveform=added_grad, system=system)

def trapezoid(amplitude, rise_time, flat_time, fall_time, dt=1e-5):
    # Time segments
    t_rise = np.arange(dt/2, rise_time, dt)
    t_flat = np.arange(dt/2, flat_time, dt)
    t_fall = np.arange(dt/2, fall_time, dt)

    # Signal segments
    rise = (amplitude / rise_time) * t_rise
    plateau = np.ones_like(t_flat) * amplitude
    fall = amplitude - (amplitude / fall_time) * t_fall

    # Concatenate full waveform
    waveform = np.concatenate((rise, plateau, fall))

    return waveform

def waveform_from_seqblock(grad, system):
    """
    extracts gradient waveform from Pypulseq sequence block
    """

    if grad.type == 'trap':
        waveform = trapezoid(grad.amplitude, grad.rise_time, grad.flat_time, grad.fall_time, dt=system.grad_raster_time)
    else:
        waveform = grad.waveform

    waveform = np.concatenate((np.zeros(round(grad.delay / system.grad_raster_time)), waveform))

    return waveform

def trap_from_area(area, system, slewrate = None, max_grad = None):
    """
    Calculate minimum time trapezoidal gradient from gradient area/moment
    In the last step the amplitude is recalculated to get the right moment and sign of the gradient

    area: Gradient area [1/m]
    system: system configuration from Pypulseq
    slewrate: (Optional) set maximum slewrate [T/m/s]
    max_grad: (Optional) set maximum gradient
    """
    if slewrate is not None:
        tmp_slew = system.max_slew
        system.max_slew = slewrate*system.gamma

    if max_grad is not None:
        tmp_grad = system.max_grad
        system.max_grad = max_grad * system.gamma

    if abs(area) < system.max_grad*round_up_to_raster(system.max_grad/system.max_slew, decimals=5):
        ftop = 0
        amp = np.sqrt(abs(area)*system.max_slew)
        ramp = round_up_to_raster(amp/system.max_slew, decimals=5) 
        amp = area/ramp

    else:
        amp = system.max_grad
        ramp = round_up_to_raster(amp/system.max_slew, decimals=5)
        ftop = round_up_to_raster(abs(area)/amp - ramp, decimals=5)
        amp = area/(ftop+ramp)

    if slewrate is not None:
        system.max_slew = tmp_slew

    if max_grad is not None:
        system.max_grad = tmp_grad

    return amp, ftop, ramp

def calc_triang_wf(amp, ftop, ramp):
    """ Calculate triangular waveform from:
    amp: amplitude [Hz]
    ftop: flat top [s]
    ramp: ramptime [s]
    """

    ramp_wf = amp * np.arange(0.5, int(ramp/10e-6+0.5)) / int(ramp/10e-6+0.5)
    wf = np.concatenate((ramp_wf, amp*np.ones(int(ftop/10e-6+0.5)), ramp_wf[::-1]))
    return wf

def rot_grad(gx, gy, phi):
    """
    rotate gradient with 2D rotation matrix
    This rotation direction matches the direction in the IDEA sequence
    """
    gx_rot = np.cos(phi)*gx+np.sin(phi)*gy
    gy_rot = -np.sin(phi)*gx+np.cos(phi)*gy    
    
    return gx_rot, gy_rot

#############
# Time raster functions
#############

def round_up_to_raster(number, decimals=5, tol=1e-10):
    """
    Round number up to a specific number of decimal places.
    Rounds up only if the digit beyond the desired precision exceeds a tolerance.
    This avoids rounding up for tiny floating point errors.
    
    Parameters:
    - number: float, the value to round
    - decimals: int, number of decimal places
    - tol: float, minimum excess to consider as real (not floating point noise)
    """
    multiplier = 10 ** decimals
    scaled = number * multiplier
    rounded = math.floor(scaled)

    if scaled - rounded > tol:
        rounded += 1

    return rounded / multiplier

def trunc_to_raster(number, decimals=0):
    """
    Returns a value truncated to a specific number of decimal places.
    """
    if not isinstance(decimals, int):
        raise TypeError("decimal places must be an integer.")
    elif decimals < 0:
        raise ValueError("decimal places has to be 0 or more.")
    elif decimals == 0:
        return math.trunc(number)

    factor = 10.0 ** decimals
    return math.trunc(number * factor) / factor

#############
# Gradient resonance check
#############

def check_resonances(grads, scanner, seq=None):
    """ Checks is the maximum frequency of the gradient waveform is in the resonance range.

    grads: list of gradient waveform on gradient raster (10us), e.g. [grad_x, grad_y]
    resonances: list of tuples describing the resonance bands in [Hz], e.g. [(100,200), (1000,1200)]
    seq: Provide sequence object to plot resonances of whole sequence (optional)
    """

    if scanner == '7tplus':
        # 7T Plus resonances
        resonances = [(500,600), (930, 1280)]
    elif scanner == 'terra':
        # 7T Terra resonances
        resonances = [(312,422), (450, 650), (900,1250)]
    elif scanner == 'skyra':
        # 3T Skyra resonances
        resonances = [(535, 635), (1010,1230)]
    elif scanner == 'connectom':
        # 3T Connectom resonances
        resonances = [(280,340), (546, 646), (1000,1500)]
    else:
        print(f"\033[93mWARNING:\033[0m No resonance frequencies defined for scanner {scanner}. Fallback to 7T Plus resonances.")
        resonances = [(500,600), (930, 1280)]

    freq_max = []
    grad_ft_list = []
    for key, grad in enumerate(grads):
        if len(grad)<20000:
            grad = np.concatenate((grad,np.zeros(20000-len(grad)))) # add zeros for higher freq resolution
        grad_ft = fft(grad, dim=-1)
        freq = np.arange(-1/(2*dt_grad), 1/(2*dt_grad), 1/(dt_grad*len(grad)))
        argmax = np.argmax(abs(grad_ft[len(grad)//2:]), axis=-1)
        freq_max.append(freq[len(grad)//2 + argmax]) # peak frequencies from gradient waveform
        grad_ft_list.append(grad_ft)
        for res in resonances:
            if (res[0] <= freq_max[key] <= res[1]):
                print(f"\033[91mWARNING:\033[0m Frequency peak {freq_max[key]:.0f} Hz of axis {key} is in the forbidden range of {resonances}. "
                      "Plot gradient spectrum and check strength of peak.")
    
    grad_ft_arr = np.array(grad_ft_list)
    plt.figure()
    for key, grad_ft in enumerate(grad_ft_arr):
        plt.plot(freq, abs(grad_ft), label=f'Grad axis {key}')
        plt.xlim(0, 2000)
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Magnitude')
    plt.legend()
    for res in resonances:
        plt.axvspan(res[0], res[1], color='red', alpha=0.3)

    if seq is not None:
        resonances_for_spectrum = [{'frequency': res[0]+(res[1]-res[0])/2, 'bandwidth': (res[1]-res[0])} for res in resonances]
        seq.calculate_gradient_spectrum(acoustic_resonances=resonances_for_spectrum)

    return freq_max
    