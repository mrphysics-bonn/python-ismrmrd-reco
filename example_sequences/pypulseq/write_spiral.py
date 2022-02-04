# Spiral Pulseq Sequence

#%%

import numpy as np
import matplotlib.pyplot as plt
import h5py
import ismrmrd
import os
import datetime
from pathlib import Path # create directories
import shutil # copy files

from pypulseq.make_arbitrary_grad import make_arbitrary_grad
from pypulseq.Sequence.sequence import Sequence
from pypulseq.make_adc import make_adc
from pypulseq.make_sinc_pulse import make_sinc_pulse
from pypulseq.make_gauss_pulse import make_gauss_pulse
from pypulseq.make_trap_pulse import make_trapezoid
from pypulseq.make_delay import make_delay
from pypulseq.make_digital_output_pulse import make_digital_output_pulse
from pypulseq.opts import Opts
from pypulseq.calc_duration import calc_duration

import spiraltraj
import pulseq_helper as ph
from prot import create_hdr
from gre_refscan import gre_refscan

#%% Parameters 
"""
PyPulseq units (SI): 
time:       [s] (not [ms] as in documentation)
spatial:    [m]
gradients:  [Hz/m] (gamma*T/m)
grad area:  [1/m]
flip angle: [rad]

Some units get converted below, others have to stay in non-SI units as spiral calculation needs different units
"""

# General
seq_name        = 'spiralout_gre_fatsat_3T' # sequence/protocol filename

# Sequence - Contrast and Geometry
fov             = 220       # field of view [mm]
TR              = 30       # repetition time [ms]
TE              = 5        # echo time [ms]
res             = 1       # in plane resolution [mm]                   
slice_res       = 1       # slice thickness [mm]
slices          = 1         # number of slices
dist_fac        = 0          # distance factor for slices [%]
averages        = 1         # number of averages

refscan         = True      # Cartesian reference scan for sensmaps
prepscans       = 64        # number of preparation/dummy scans before GRE
noisescans      = 16        # number of noise scans

# ADC
os_factor = 2               # oversampling factor (automatic 2x os from Siemens is not applied)

# RF
flip_angle      = 10        # flip angle of excitation pulse [°]
rf_dur          = 3       # RF duration [ms]
tbp_exc         = 5         # time bandwidth product excitation pulse
rf_spoiling     = False     # RF spoiling

fatsat          = True      # Fat saturation pulse
fatsat_dur      = 2.1       # duration of fatsat pulse [ms] (increase at low TR to avoid SAR problems)

# Gradients
max_slew        = 120       # maximum slewrate [T/m/s] (system limit)
spiral_slew     = 100       # maximum slew rate of spiral gradients - for pre Emph: set lower than max_slew
max_grad        = 40        # maximum gradient amplitude [mT/m] (system limit)
max_grad_sp     = 30        # maximum gradient amplitude of spiral gradients - for pre_emph: set lower than max_grad

Nintl           = 15         # spiral interleaves
redfac          = 1         # reduction/acceleration factor
spiraltype      = 1         # 1: Spiral Out, 4: ROI, WIP: other spiral waveforms
spiral_os       = 1         # spiral oversampling in center

#%% Limits, checks and preparations

# Set System limits
rf_dead_time = 100e-6 # lead time before rf can be applied
rf_ringdown_time = 30e-6 # coil hold time (20e-6) + frequency reset time (10e-6)
system = Opts(max_grad=max_grad, grad_unit='mT/m', max_slew=max_slew, slew_unit='T/m/s', rf_dead_time=rf_dead_time, rf_ringdown_time=rf_ringdown_time)

# convert parameters to Pulseq units
TR          *= 1e-3 # [s]
TE          *= 1e-3 # [s]
rf_dur      *= 1e-3 # [s]
slice_res   *= 1e-3 # [m]
fatsat_dur  *= 1e-3 # [s]

if Nintl/redfac%1 != 0:
    raise ValueError('Number of interleaves is not multiple of reduction factor')

if spiraltype!=1 and spiraltype!=4:
    ValueError('Right now only spiraltype 1 (spiral out) and 4 (ROI) possible.')
    
#%% RF Pulse and slab/slice selection gradient

# make rf pulse and calculate duration of excitation and rewinding
rf, gz, gz_rew, rf_del = make_sinc_pulse(flip_angle=flip_angle*np.pi/180, system=system, duration=rf_dur, slice_thickness=slice_res,
                            apodization=0.5, time_bw_product=tbp_exc, use='excitation', return_gz=True, return_delay=True)
exc_to_rew = rf_del.delay - rf_dur/2 - rf.delay # time from middle of rf pulse to rewinder, rf_del.delay equals the block length
rew_dur = calc_duration(gz_rew)

# RF spoiling parameters
rf_spoiling_inc = 50 # increment of RF spoiling [°]
rf_phase        = 0 
rf_inc          = 0

# Fat saturation
if fatsat:
    fatsat_bw = 1000 # bandwidth [Hz] (1000 Hz is approx. used by Siemens)
    fatsat_fa = 110 # flip angle [°]

    rf_fatsat, fatsat_del = make_gauss_pulse(flip_angle=fatsat_fa*np.pi/180, duration=fatsat_dur, bandwidth=fatsat_bw, freq_offset=ph.fw_shift, system=system, return_delay=True)

# echo time delay
min_te = exc_to_rew + rew_dur
if min_te > TE:
    raise ValueError('Minimum TE is {} ms.'.format(min_te*1e3))
te_delay = make_delay(d = TE - min_te)

#%% Spiral Readout Gradients

# Parameters spiral trajectory:

# parameter         description               default value
# ---------        -------------              --------------

# nitlv:      number of spiral interleaves        15
# res:        resolution                          1 mm
# fov:        target field of view                192 mm
# max_amp:    maximum gradient amplitude          42 mT/m
# min_rise:   minimum gradient risetime           5 us/(mT/m)
# spiraltype: 1: spiral out                   
#             2: spiral in                        
#             3: double spiral                    x
#             4: ROI
#             5: RIO
# spiral_os:  spiral oversampling in center       1

# Rotation of Spirals
if spiraltype==3:
    max_rot     = np.pi    # maximum rotation angle of spirals [rad]
else:
    max_rot     = 2*np.pi  

# read in Spirals [T/m]
min_rise_sp = 1/spiral_slew * 1e3
spiral_calc = spiraltraj.calc_traj(nitlv=Nintl, fov=fov, res=res, spiraltype=spiraltype, min_rise=min_rise_sp, max_amp=max_grad_sp, spiral_os=spiral_os)
spiral_calc = np.asarray(spiral_calc)
spiral_x = 1e-3*spiral_calc[:,0]
spiral_y = 1e-3*spiral_calc[:,1]

N_spiral = len(spiral_x)
readout_dur = N_spiral*system.grad_raster_time # readout duration [s]

# write spiral readout blocks to list
spirals = [{'deph': [None, None], 'spiral': [None, None], 'reph': [None, None]} for k in range(Nintl)]
reph_dur = []
save_sp = np.zeros((Nintl, 2, N_spiral)) # save gradients for FIRE reco
rot_angle = np.linspace(0, max_rot, Nintl, endpoint=False)
for k in range(Nintl):
    # rotate spiral gradients for shot selection
    sp_x, sp_y = ph.rot_grad(spiral_x, spiral_y, rot_angle[k])

    save_sp[k,0,:] = sp_x
    save_sp[k,1,:] = sp_y

    # unit to [Hz/m], make spiral gradients
    sp_x *= system.gamma
    sp_y *= system.gamma
    if spiraltype==1:
        spiral_delay = 2e-5 # delay the spirals to have some ADC samples before the start of the spiral
    else:
        spiral_delay = 0 # no spiral delay as sampling points should be mirrored in ROI
    spirals[k]['spiral'][0] = make_arbitrary_grad(channel='x', waveform=sp_x, delay=spiral_delay, system=system)
    spirals[k]['spiral'][1] = make_arbitrary_grad(channel='y', waveform=sp_y, delay=spiral_delay, system=system)

    if spiraltype==1:
        # calculate rephaser area
        area_x = sp_x.sum()*system.grad_raster_time
        area_y = sp_y.sum()*system.grad_raster_time

        # calculate rephasers and make gradients
        amp_x, ftop_x, ramp_x = ph.trap_from_area(-area_x, system, slewrate = 100) # reduce slew rate to 100 T/m/s to avoid stimulation
        amp_y, ftop_y, ramp_y = ph.trap_from_area(-area_y, system, slewrate = 100)
        spirals[k]['reph'][0] = make_trapezoid(channel='x', system=system, amplitude=amp_x, flat_time=ftop_x, rise_time=ramp_x)
        spirals[k]['reph'][1] = make_trapezoid(channel='y', system=system, amplitude=amp_y, flat_time=ftop_y, rise_time=ramp_y)
        reph_dur.append(max(ftop_x+2*ramp_x, ftop_y+2*ramp_y))


# check for acoustic resonances (checks only spirals)
freq_max = ph.check_resonances([spiral_x,spiral_y])

#%% Gradient Spoiler on slice axis

spoiler_area = 2*gz.flat_area - gz.area/2 # 2x moment under excitation pulse
amp_spoil, ftop_spoil, ramp_spoil = ph.trap_from_area(spoiler_area, system, slewrate=100) # reduce slew rate to 100 T/m/s to avoid stimulation
spoiler_z = make_trapezoid(channel='z',system=system, amplitude=amp_spoil, flat_time=ftop_spoil, rise_time=ramp_spoil, delay=100e-6)
spoiler_dur = calc_duration(spoiler_z)

#%% ADC

max_grad_sp_cmb = 1e3*np.max(np.sqrt(abs(spiral_x)**2+abs(spiral_y)**2))
dwelltime = 1/(system.gamma*max_grad_sp_cmb*fov*os_factor)*1e6 # ADC dwelltime [s]
dwelltime = ph.trunc_to_raster(dwelltime, decimals=7) # truncate dwelltime to 100 nanoseconds (scanner limit)
min_dwelltime = 1e-6
if dwelltime < min_dwelltime:
    dwelltime = min_dwelltime
print("ADC dwelltime: {}".format(dwelltime))

num_samples = round((readout_dur+spiral_delay)/dwelltime)
if num_samples%2==1:
    num_samples += 1 # even number of samples

if num_samples <= 8192:
    num_segments = 1
    print('Number of ADCs: {}.'.format(num_samples))
else:
    # the segment duration has to be on the gradient raster
    # increase number of segments or samples/segments to achieve this
    # number of samples and number of samples per segment should always be an even number
    num_segments = 2
    if (num_samples/num_segments % 2 != 0):
        num_samples += 2
    segm_dur = 1e5 * dwelltime * num_samples/num_segments # segment duration [10us - gradient raster]
    while (not round(segm_dur,ndigits=5).is_integer() or num_samples/num_segments > 8192):
        if num_samples/num_segments > 8192:
            num_segments += 1
            while (num_samples/num_segments % 2 != 0):
                num_samples += 2
        else:
            num_samples += 2*num_segments
        segm_dur = 1e5 * dwelltime * num_samples/num_segments 
    print('ADC has to be segmented!! Number of ADCs: {}. Per segment: {}. Segments: {}.'.format(num_samples,num_samples/num_segments,num_segments))

    # self check
    if (num_samples/num_segments % 2 != 0 or num_samples % 2 != 0 or not round(segm_dur,ndigits=5).is_integer()):
        raise ValueError("Check if number of samples and number of samples per segment are even. Check if segment duration is on gradient raster time.")

if num_samples > 65535: # max of uint16 used by ISMRMRD
    raise ValueError("Too many samples for ISMRMRD format - lower the oversampling factor or take more interleaves")

adc = make_adc(system=system, num_samples=num_samples, dwell=dwelltime)
adc_dur = num_samples * dwelltime
adc_delay = ph.round_up_to_raster(adc_dur+200e-6, decimals=5) # add small delay after readout for ADC frequency reset event and to avoid stimulation by rephaser
adc_delay = make_delay(d=adc_delay)

#%% Set up protocol for FIRE reco and write header

if os.path.isfile(seq_name+'.h5'):
    os.remove(seq_name+'.h5')
prot = ismrmrd.Dataset(seq_name+'.h5')
hdr = ismrmrd.xsd.ismrmrdHeader()
t_min = TE + dwelltime/2 # save trajectory starting point for B0-correction
params_hdr = {"fov": fov, "res": res, "slices": slices, "slice_res": slice_res, "nintl":int(Nintl/redfac), "avg": averages,
                "nsegments": num_segments, "dwelltime": dwelltime, "traj_delay": spiral_delay, "t_min": t_min, "trajtype": "spiral"}
create_hdr(hdr, params_hdr)
prot.write_xml_header(hdr.toXML('utf-8'))

#%% Add sequence blocks to sequence & write acquisitions to protocol

# Set up the sequence
seq = Sequence()

# Definitions section in seq file
seq.set_definition("Name", seq_name) # protocol name is saved in Siemens header for FIRE reco
seq.set_definition("FOV", [1e-3*fov, 1e-3*fov, slice_res]) # for FOV positioning
seq.set_definition("Slice_Thickness", "%f" % (slice_res*(1+dist_fac*1e-2)*(slices-1)+slice_res)) # we misuse this to show the total covered head area in the GUI
if num_segments > 1:
    seq.set_definition("MaxAdcSegmentLength", "%d" % int(num_samples/num_segments+0.5)) # for automatic ADC segment length setting

# Noise scans
noise_samples = 256
noise_dwelltime = 2e-6
noise_adc = make_adc(system=system, num_samples=256, dwell=noise_dwelltime)
noise_delay = make_delay(d=ph.round_up_to_raster(noise_adc.duration+1e-3,decimals=5)) # add some more time to the ADC delay to be safe
for k in range(noisescans):
    seq.add_block(noise_adc, noise_delay)
    acq = ismrmrd.Acquisition()
    acq.setFlag(ismrmrd.ACQ_IS_NOISE_MEASUREMENT)
    prot.append_acquisition(acq)

# Perform cartesian reference scan: if selected / for accelerated spirals / for long readouts
if refscan == False:
    if redfac > 1:
        refscan = True
        print("Accelerated scan: Activate Cartesian reference scan.")

if refscan:
    res_refscan = res*1e-3 * 2
    flip_refscan = 15
    bw_refscan = 800
    params_ref = {"fov":fov*1e-3, "res":res_refscan, "slices":slices, "slice_res":slice_res, "dist_fac": dist_fac, "flip_angle":flip_refscan,
     "rf_dur":rf_dur, "tbp": tbp_exc, "readout_bw": bw_refscan}
    gre_refscan(seq, prot=prot, system=system, params=params_ref)


""" GRE

The following code generates a spoiled gradient echo (GRE) spiral sequence.
"""

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
    rf.freq_offset = gz.amplitude * slice_res * (slc - (slices - 1) / 2) * (1+dist_fac*1e-2)

    # preparation scans without adc
    for k in range(prepscans): 
        if rf_spoiling:
            rf.phase_offset  = rf_phase / 180 * np.pi
            adc.phase_offset = rf_phase / 180 * np.pi
        if fatsat:
            rf_fatsat.phase_offset = rf_phase / 180 * np.pi
            seq.add_block(rf_fatsat, fatsat_del)
            seq.add_block(spoiler_z)
        rf_inc = divmod(rf_inc + rf_spoiling_inc, 360.0)[1]
        rf_phase = divmod(rf_phase + rf_inc, 360.0)[1]
        
        seq.add_block(rf,gz,rf_del)
        seq.add_block(gz_rew)
        seq.add_block(te_delay)
        seq.add_block(spirals[0]['spiral'][0], spirals[0]['spiral'][1], adc_delay)
        if spiraltype==1:
            seq.add_block(spirals[0]['reph'][0], spirals[0]['reph'][1], spoiler_z)
            min_tr = exc_to_rew + TE + adc_delay.delay + max(reph_dur[0], spoiler_dur)
        else:
            seq.add_block(spoiler_z)
            min_tr = exc_to_rew + TE + adc_delay.delay + spoiler_dur  

        if fatsat:
            min_tr += fatsat_del.delay + spoiler_dur
        if TR < min_tr:
            raise ValueError('Minimum TR is {} ms.'.format(min_tr*1e3))
        tr_delay = make_delay(d=TR-min_tr)
        seq.add_block(tr_delay)

    # imaging scans
    for avg in range(averages):
        for n in range(int(Nintl/redfac)):
            if rf_spoiling:
                rf.phase_offset  = rf_phase / 180 * np.pi
                adc.phase_offset = rf_phase / 180 * np.pi
            if fatsat:
                rf_fatsat.phase_offset = rf_phase / 180 * np.pi # always use RF spoiling for fat sat pulse
                seq.add_block(rf_fatsat, fatsat_del)
                seq.add_block(spoiler_z)
            rf_inc = divmod(rf_inc + rf_spoiling_inc, 360.0)[1]
            rf_phase = divmod(rf_phase + rf_inc, 360.0)[1]

            # excitation
            seq.add_block(rf,gz,rf_del)
            seq.add_block(gz_rew)

            # spiral readout block with spoiler gradient
            seq.add_block(te_delay)
            seq.add_block(spirals[n*redfac]['spiral'][0], spirals[n*redfac]['spiral'][1], adc, adc_delay)
            if spiraltype==1:
                seq.add_block(spirals[n*redfac]['reph'][0], spirals[n*redfac]['reph'][1], spoiler_z)
                min_tr = exc_to_rew + TE + adc_delay.delay + max(reph_dur[n*redfac], spoiler_dur)
            else:
                seq.add_block(spoiler_z)
                min_tr = exc_to_rew + TE + adc_delay.delay + spoiler_dur

            # delay at end of one TR
            if fatsat:
                min_tr += fatsat_del.delay + spoiler_dur
            if TR < min_tr:
                raise ValueError('Minimum TR is {} ms.'.format(min_tr*1e3))
            tr_delay = make_delay(d=TR-min_tr)
            seq.add_block(tr_delay)

            # add protocol information
            for seg in range(num_segments):
                acq = ismrmrd.Acquisition()
                if (n == int(Nintl/redfac) - 1) and (seg == num_segments - 1):
                    acq.setFlag(ismrmrd.ACQ_LAST_IN_SLICE)
                acq.idx.kspace_encode_step_1 = n
                acq.idx.slice = slc
                acq.idx.average = avg
                acq.idx.segment = seg

                # save gradient only in first segment to save space
                if seg == 0:
                    # we misuse the trajectory field for the gradient array - channels are hardcoded
                    acq.resize(trajectory_dimensions = save_sp.shape[1], number_of_samples=save_sp.shape[2], active_channels=32)
                    acq.traj[:] = np.swapaxes(save_sp[n*redfac],0,1) # [samples, dims]
                prot.append_acquisition(acq)
    
    slc += 2 # acquire every 2nd slice, afterwards fill slices inbetween

        # intl
    # avg
# slices

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
