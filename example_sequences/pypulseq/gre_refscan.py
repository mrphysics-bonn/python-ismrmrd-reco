""" This function adds a GRE reference scan to a sequence for calculation of sensitivity maps

WIP: - make this faster (e.g. with partial Fourier?)
     - combine this with a B0 field map
     - maybe 3D acquisition?

"""

import math
import numpy as np

import ismrmrd

from pypulseq.Sequence.sequence import Sequence
from pypulseq.calc_duration import calc_duration
from pypulseq.make_adc import make_adc
from pypulseq.make_delay import make_delay
from pypulseq.make_sinc_pulse import make_sinc_pulse
from pypulseq.make_trap_pulse import make_trapezoid
from pypulseq.opts import Opts

import pulseq_helper as ph

def gre_refscan(seq, prot=None, system=Opts(), params=None):

    # decrease slew rate a bit
    save_slew = system.max_slew
    system.max_slew = 100 * system.gamma
    if params is None:
        params = {"fov":210e-3, "res":3e-3, "flip_angle":12, "rf_dur":1e-3, "tbp": 2, "slices":1, "slice_res":2e-3, "dist_fac":0, "readout_bw": 600}

    # RF
    rf, gz, gz_reph, rf_del = make_sinc_pulse(flip_angle=params["flip_angle"] * math.pi / 180, duration=params["rf_dur"], slice_thickness=params["slice_res"],
                                apodization=0.5, time_bw_product=params["tbp"], system=system, return_gz=True, return_delay=True)

    # Calculate readout gradient and ADC parameters
    delta_k = 1 / params["fov"]
    Nx = Ny = int(params["fov"]/params["res"]+0.5)
    samples = 2*Nx # 2x oversampling
    gx_flat_time_us = int(1e6/params["readout_bw"]) # readout_bw is in Hz/Px
    dwelltime = ph.trunc_to_raster(1e-6*gx_flat_time_us / samples, decimals=7)
    gx_flat_time = round(dwelltime*samples, 5)
    if (1e5*gx_flat_time %2 == 1):
        gx_flat_time += 10e-6 # even flat time
    diff_flat_adc = gx_flat_time - (dwelltime*samples)

    # Gradients
    gx_flat_area = Nx * delta_k * (gx_flat_time / (dwelltime*samples)) # compensate for longer flat time than ADC
    gx = make_trapezoid(channel='x', flat_area=gx_flat_area, flat_time=gx_flat_time, system=system)
    gx_pre = make_trapezoid(channel='x', area=-gx.area / 2, duration=1.4e-3, system=system)
    phase_areas = (np.arange(Ny) - Ny / 2) * delta_k

    # reduce slew rate of spoilers to avoid stimulation
    gx_spoil = make_trapezoid(channel='x', area=2 * Nx * delta_k, system=system, max_slew=120*system.gamma)
    gz_spoil = make_trapezoid(channel='z', area=4 / params["slice_res"], system=system, max_slew=120*system.gamma)

    # take minimum TE rounded up to .1 ms
    min_TE = np.ceil((gz.fall_time + gz.flat_time / 2 + calc_duration(gx_pre) + calc_duration(gx) / 2) / seq.grad_raster_time) * seq.grad_raster_time
    TE = ph.round_up_to_raster(min_TE, decimals=4)
    delay_TE = TE-min_TE

    # take minimum TR rounded up to .1 ms
    min_TR = calc_duration(gx_pre) + calc_duration(gz) + calc_duration(gx) + delay_TE + calc_duration(gx_spoil, gz_spoil)
    TR = ph.round_up_to_raster(min_TR, decimals=4)
    delay_TR = TR - min_TR

    # ADC with 2x oversampling
    adc = make_adc(num_samples=samples, dwell=dwelltime, delay=gx.rise_time+diff_flat_adc/2, system=system)
    
    # RF spoiling
    rf_spoiling_inc = 117
    rf_phase = 0
    rf_inc = 0

    # build sequence
    prepscans = 40 # number of dummy preparation scans

    if params["slices"]%2 == 1:
        slc = 0
    else:
        slc = 1
    for s in range(params["slices"]):
        if s==int(params["slices"]/2+0.5): 
            if params["slices"]%2 == 1:
                slc = 1
            else:
                slc = 0
        rf.freq_offset = gz.amplitude * params["slice_res"] * (slc - (params["slices"] - 1) / 2) * (1+params["dist_fac"]*1e-2)

        # prepscans
        for d in range(prepscans):
            rf.phase_offset = rf_phase / 180 * np.pi
            adc.phase_offset = rf_phase / 180 * np.pi
            rf_inc = divmod(rf_inc + rf_spoiling_inc, 360.0)[1]
            rf_phase = divmod(rf_phase + rf_inc, 360.0)[1]

            seq.add_block(rf, gz, rf_del)
            gy_pre = make_trapezoid(channel='y', area=phase_areas[0], duration=1.4e-3, system=system)
            seq.add_block(gx_pre, gy_pre, gz_reph)
            seq.add_block(make_delay(delay_TE))
            seq.add_block(gx)
            gy_pre.amplitude = -gy_pre.amplitude
            seq.add_block(make_delay(delay_TR), gx_spoil, gy_pre, gz_spoil)

        # imaging scans
        for i in range(Ny):
            rf.phase_offset = rf_phase / 180 * np.pi
            adc.phase_offset = rf_phase / 180 * np.pi
            rf_inc = divmod(rf_inc + rf_spoiling_inc, 360.0)[1]
            rf_phase = divmod(rf_phase + rf_inc, 360.0)[1]

            seq.add_block(rf, gz, rf_del)
            gy_pre = make_trapezoid(channel='y', area=phase_areas[i], duration=1.4e-3, system=system)
            seq.add_block(gx_pre, gy_pre, gz_reph)
            seq.add_block(make_delay(delay_TE))
            seq.add_block(gx, adc)
            gy_pre.amplitude = -gy_pre.amplitude
            seq.add_block(make_delay(delay_TR), gx_spoil, gy_pre, gz_spoil)

            if prot is not None:
                acq = ismrmrd.Acquisition()
                acq.idx.kspace_encode_step_1 = i
                acq.idx.kspace_encode_step_2 = 0 # only 2D atm
                acq.idx.slice = slc
                # acq.idx.average = avg
                acq.setFlag(ismrmrd.ACQ_IS_PARALLEL_CALIBRATION)
                if i == Ny-1:
                    acq.setFlag(ismrmrd.ACQ_LAST_IN_SLICE)
                prot.append_acquisition(acq)
                
        slc += 2 # acquire every 2nd slice, afterwards fill slices inbetween

    delay_end = make_delay(d=2) # 5s delay after reference scan to allow for relaxation
    seq.add_block(delay_end)
    system.max_slew = save_slew
    