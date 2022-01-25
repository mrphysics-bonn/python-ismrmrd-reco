
from inspect import signature
import math
import numpy as np
import datetime
import os

import ismrmrd

from pypulseq.Sequence.sequence import Sequence
from pypulseq.calc_duration import calc_duration
from pypulseq.make_adc import make_adc
from pypulseq.make_delay import make_delay
from pypulseq.make_sinc_pulse import make_sinc_pulse
from pypulseq.make_trap_pulse import make_trapezoid
from pypulseq.opts import Opts

import pulseq_helper as ph
from prot import create_hdr

# Parameters
seq_name = 'gre_b0mapping_hashed' # sequence/protocol filename
TE = [2.04e-3, 4.08e-3] # [TE1, TE2] atm only 2 echo times supported
fov = 220e-3
res = 2e-3
flip_angle = 15
rf_dur = 0.8e-3
tbp = 2
slices = 1
slice_res = 2e-3
readout_bw = 1200
max_grad = 35
max_slew = 150
prepscans = 40 # number of dummy preparation scans

# System limits
seq = Sequence()
rf_dead_time = 100e-6 # lead time before rf can be applied
rf_ringdown_time = 30e-6 # coil hold time (20e-6) + frequency reset time (10e-6)
system = Opts(max_grad=max_grad, grad_unit='mT/m', max_slew=max_slew, slew_unit='T/m/s', rf_dead_time=rf_dead_time, rf_ringdown_time=rf_ringdown_time)

# RF
rf, gz, gz_reph, rf_del = make_sinc_pulse(flip_angle=flip_angle * math.pi / 180, duration=rf_dur, slice_thickness=slice_res,
                            apodization=0.5, time_bw_product=tbp, system=system, return_gz=True, return_delay=True)

# Calculate readout gradient and ADC parameters
delta_k = 1 / fov
Nx = Ny = int(fov/res+0.5)
samples = 2*Nx
gx_flat_time_us = int(1e6/readout_bw) # readout_bw is in Hz/Px
dwelltime_us = gx_flat_time_us / samples
gx_flat_time = round(1e-6*dwelltime_us*samples, 5)

# Gradients
gx = make_trapezoid(channel='x', flat_area=Nx * delta_k, flat_time=gx_flat_time, system=system)
gx_mid = make_trapezoid(channel='x', area=-gx.area, system=system)
gx_pre = make_trapezoid(channel='x', area=-gx.area / 2, system=system)
gx_pre.delay = calc_duration(gz_reph) - calc_duration(gx_pre)
phase_areas = (np.arange(Ny) - Ny / 2) * delta_k

# spoilers
gx_spoil = make_trapezoid(channel='x', area=2 * Nx * delta_k, system=system)
gz_spoil = make_trapezoid(channel='z', area=4 / slice_res, system=system)

# calculate TE delays
max_gy_pre = make_trapezoid(channel='y', area=max(abs(phase_areas)), system=system)
gy_pre_dur = calc_duration(max_gy_pre)
min_TE = np.ceil((gz.fall_time + gz.flat_time / 2 + calc_duration(gx_pre, max_gy_pre, gz_reph) + calc_duration(gx) / 2) / seq.grad_raster_time) * seq.grad_raster_time
delay_TE1 = TE[0] - min_TE
delay_TE2 = TE[1] - TE[0] - calc_duration(gx) - calc_duration(gx_mid)
if delay_TE1 < 0:
    raise ValueError(f"TE 1 too small by {1e3*abs(delay_TE1)} ms. Increase readout bandwidth.")
if delay_TE2 < 0:
    raise ValueError(f"TE 2 too small by {1e3*abs(delay_TE2)} ms. Increase readout bandwidth.")

# ADC with 2x oversampling
adc = make_adc(num_samples=samples, dwell=1e-6*dwelltime_us, delay=gx.rise_time, system=system)

# RF spoiling
rf_spoiling_inc = 117
rf_phase = 0
rf_inc = 0

#%% Set up protocol

if os.path.isfile(seq_name+'.h5'):
    os.remove(seq_name+'.h5')
prot = ismrmrd.Dataset(seq_name+'.h5')
hdr = ismrmrd.xsd.ismrmrdHeader()
params_hdr = {"trajtype": "cartesian", "fov": fov*1e3, "res": res*1e3, "slices": slices, "slice_res": slice_res, "nintl": Ny, "ncontrast": len(TE)}
create_hdr(hdr, params_hdr)

#%% Set up sequence

seq.set_definition("Name", seq_name) # protocol name is saved in Siemens header for FIRE reco
seq.set_definition("FOV", [fov, fov, slice_res])
seq.set_definition("Slice_Thickness", "%f" %(slice_res*slices))

# Noise scans
noisescans = 16
noise_samples = 256
noise_adc = make_adc(system=system, num_samples=256, dwell=1e-6*dwelltime_us)
noise_delay = make_delay(d=ph.round_up_to_raster(noise_adc.duration+1e-3,decimals=5))
for k in range(noisescans):
    seq.add_block(noise_adc, noise_delay)
    acq = ismrmrd.Acquisition()
    acq.setFlag(ismrmrd.ACQ_IS_NOISE_MEASUREMENT)
    prot.append_acquisition(acq)

if slices%2 == 1:
    slc = 0
else:
    slc = 1
for s in range(slices):
    if s==int(slices/2+0.5): 
        if slices%2 == 1:
            slc = 1
        else:
            slc = 0
    rf.freq_offset = gz.amplitude * slice_res * (slc - (slices - 1) / 2)

    # prepscans
    for d in range(prepscans):
        rf.phase_offset = rf_phase / 180 * np.pi
        adc.phase_offset = rf_phase / 180 * np.pi
        rf_inc = divmod(rf_inc + rf_spoiling_inc, 360.0)[1]
        rf_phase = divmod(rf_phase + rf_inc, 360.0)[1]

        seq.add_block(rf, gz, rf_del)
        gy_pre = make_trapezoid(channel='y', area=0, duration=gy_pre_dur, system=system)
        seq.add_block(gx_pre, gy_pre, gz_reph)
        seq.add_block(make_delay(delay_TE1))
        seq.add_block(gx)
        seq.add_block(gx_mid)
        seq.add_block(make_delay(delay_TE2))
        seq.add_block(gx)
        gy_pre.amplitude = -gy_pre.amplitude
        seq.add_block(gx_spoil, gy_pre, gz_spoil)

    # imaging scans
    for i in range(Ny):
        rf.phase_offset = rf_phase / 180 * np.pi
        adc.phase_offset = rf_phase / 180 * np.pi
        rf_inc = divmod(rf_inc + rf_spoiling_inc, 360.0)[1]
        rf_phase = divmod(rf_phase + rf_inc, 360.0)[1]

        seq.add_block(rf, gz, rf_del)
        gy_pre = make_trapezoid(channel='y', area=phase_areas[i], duration=gy_pre_dur, system=system)
        seq.add_block(gx_pre, gy_pre, gz_reph)
        seq.add_block(make_delay(delay_TE1))
        seq.add_block(gx, adc)
        seq.add_block(gx_mid)
        seq.add_block(make_delay(delay_TE2))
        seq.add_block(gx, adc)
        gy_pre.amplitude = -gy_pre.amplitude
        seq.add_block(gx_spoil, gy_pre, gz_spoil)

        if prot is not None:
            for k in range(len(TE)):
                acq = ismrmrd.Acquisition()
                acq.idx.kspace_encode_step_1 = i
                acq.idx.kspace_encode_step_2 = 0 # only 2D atm
                acq.idx.slice = slc
                acq.idx.contrast = k
                if i == Ny-1:
                    acq.setFlag(ismrmrd.ACQ_LAST_IN_SLICE)
                prot.append_acquisition(acq)
            
    slc += 2 # acquire every 2nd slice, afterwards fill slices inbetween
    
# use arrays to save b-values and directions in protocol
prot.append_array("echo_times", np.asarray(TE))

# write sequence and add hash to protocol
seq.write(seq_name+'.seq')
seq_hash = seq.get_hash()
signature = ismrmrd.xsd.userParameterStringType()
signature.name = 'seq_signature'
signature.value = seq_hash
hdr.userParameters.userParameterString.append(signature)

prot.write_xml_header(hdr.toXML('utf-8'))
prot.close()

# Optional: Add first chars of hash to sequence name
os.rename(seq_name+'.seq', f'{seq_name}_{seq_hash[:5]}.seq')
os.rename(seq_name+'.h5', f'{seq_name}_{seq_hash[:5]}.h5')
