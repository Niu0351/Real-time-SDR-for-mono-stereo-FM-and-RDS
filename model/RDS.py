import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math
import fmPllBlockRDS
import fmRRC

from fmSupportLib import fmDemodArctan
from fmMonoBasic import lp_impulse, conv

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

demod_fs = 2.4e6//10

audio_Fs = 48e3
audio_Fc = 16e3
resampling_taps = 151*19
RDS_taps = 151
audio_decim = 5

def bandpass(fb, fe, fs, N):
    coeff = np.zeros(N)
    norm_center=((fe+fb)/2)/(fs/2)
    norm_pass=(fe-fb)/(fs/2)
    for i in range(N):
        if (i==(N-1)/2):
            coeff[i] = norm_pass
        else:
            coeff[i] = norm_pass*math.sin(np.pi*(norm_pass/2)*(i-(N-1)/2))/(np.pi*(norm_pass/2)*(i-(N-1)/2))
        coeff[i]=coeff[i]*math.cos(i*np.pi*norm_center)
        coeff[i]=coeff[i]*math.sin(i*np.pi/N)*math.sin(i*np.pi/N)
    return coeff

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is normalized between -1 and +1 and interleaved
	in_fname = "../data/my_samples_u8.raw"
	iq_data = np.fromfile(in_fname, dtype='uint8')
	iq_data = np.subtract((iq_data.astype('float32')), 128)
	iq_data = np.divide(iq_data, 128)
	# iq_data = iq_data[:len(iq_data)]
	print("Read raw RF data from \"" + in_fname + "\" in float32 format")
	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, \
							 rf_Fc/(rf_Fs/2), \
							 window=('hann'))
	# rf_coeff = lp_impulse(rf_Fc, rf_Fs, rf_taps)

	# coefficients for the filters to extract audio
	RDS_channel_coeff = bandpass(54e3, 60e3, demod_fs, RDS_taps) # signal.firwin(RDS_taps, [54e3, 60e3], window=('hann'), fs=demod_Fs)
	# RDS_extra_bp_coeff = bandpass(55e3, 59e3, demod_fs, RDS_taps) # signal.firwin(stereo_taps, [113.5e3, 114.5e3], window=('hann'), fs=demod_Fs) # correct the sampling rate for filters

	RDS_carrier_coeff = bandpass(113.5e3, 114.5e3, demod_fs, RDS_taps) # signal.firwin(stereo_taps, [113.5e3, 114.5e3], window=('hann'), fs=demod_Fs) # correct the sampling rate for filters

	RDS_demod_LPF_coeff = signal.firwin(RDS_taps, 3e3, window=('hann'), fs=demod_fs)
	RDS_sampling_coeff = signal.firwin(resampling_taps, 16e3, window=('hann'), fs=demod_fs*19)
	RDS_RRC_coeff = fmRRC.impulseResponseRootRaisedCosine(demod_fs/80*19, 151)

	# set up drawing
	# fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)
	# fig.subplots_adjust(hspace = 1.0)

	# select a block_size that is in KB and
	# a multiple of decimation factors
	block_size = 1024 * rf_decim * audio_decim * 2
	block_count = 0

	# states needed for continuity in block processing
	state_i_lpf_100k = signal.lfilter_zi(rf_coeff, 1)
	state_q_lpf_100k = signal.lfilter_zi(rf_coeff, 1)
	state_pll = {
		"integrator" : 0.0,
		"phaseEst" : 0.0,
		"feedbackI" : 1.0,
		"feedbackQ" : 0.0,
		"ncoOut_I_buf" : 1.0,
		"ncoOut_Q_buf" : 0.0,
		"trigOffset": 0
	}
	state_phase = 0
	buf_I = 0.0
	buf_Q = 0.0

	# add state as needed for the mono channel filter
	RDS_channel_state = signal.lfilter_zi(RDS_channel_coeff, 1.0)
	RDS_carrier_state = signal.lfilter_zi(RDS_carrier_coeff, 1.0)
	RDS_demod_LPF_state = []
	RDS_demod_LPF_state.append(signal.lfilter_zi(RDS_demod_LPF_coeff, 1.0))
	RDS_demod_LPF_state.append(signal.lfilter_zi(RDS_demod_LPF_coeff, 1.0))
	RDS_sampling_state = []
	RDS_sampling_state.append(signal.lfilter_zi(RDS_sampling_coeff, 1.0))
	RDS_sampling_state.append(signal.lfilter_zi(RDS_sampling_coeff, 1.0))
	RDS_RRC_state = [signal.lfilter_zi(RDS_RRC_coeff, 1.0), signal.lfilter_zi(RDS_RRC_coeff, 1.0)]
	# audio buffer that stores all the audio blocks
	stereo_out_right = np.array([]) # to be updated by you during in-lab
	stereo_out_left = np.array([])
	mono_out = np.array([]) # to be updated by you during in-lab
	# if the number of samples in the last block is less than the block size
	# it is fine to ignore the last few samples from the raw IQ file
	while (block_count+1)*block_size < len(iq_data):
		# if you wish to have shorter runtimes while troubleshooting
		# you can control the above loop exit condition as you see fit

		print('Processing block ' + str(block_count))
		# Extract I and Q signals
		i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
												  iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
												  zi=state_i_lpf_100k)
		#
		q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
												  iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
												  zi=state_q_lpf_100k)

		# i_filt, state_i_lpf_100k = conv(iq_data[(block_count)*block_size:(block_count+1)*block_size:2], rf_coeff, state_i_lpf_100k)
		# q_filt, state_q_lpf_100k = conv(iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2], rf_coeff, state_q_lpf_100k)

		# downsample the FM channel
		i_ds = i_filt[::rf_decim]
		q_ds = q_filt[::rf_decim]

		# FM demodulator
		# fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)
		fm_demod, buf_I, buf_Q = fmDemodArctan(i_ds, q_ds, block_count, buf_I, buf_Q)

		# Extract RDS channel signal
		RDS_channel, RDS_channel_state = signal.lfilter(RDS_channel_coeff, 1.0, fm_demod, zi=RDS_channel_state)
		
		# Square non-linearity
		RDS_carrier_squared = RDS_channel*RDS_channel
		
		# Carrier Recovery
		RDS_carrier_filtered, RDS_carrier_state = signal.lfilter(RDS_carrier_coeff, 1.0, RDS_carrier_squared, zi=RDS_carrier_state)
		
		# Generate PLL 57k tone from filtered carrier signal
		RDS_I, RDS_Q, state_pll = fmPllBlockRDS.fmPll(RDS_carrier_filtered, 114e3, demod_fs, state_pll, ncoScale=0.5, normBandwidth=0.005)

		RDS_demod_inputs = [RDS_I, RDS_Q]
		demod_input_size = len(RDS_demod_inputs[0])-1

		#Set up arrays for demodulation process, two channels in each, one for i one for q
		mixed_RDS = [np.zeros(demod_input_size), np.zeros(demod_input_size)]
		filtered_mixed_RDS = [np.zeros(demod_input_size), np.zeros(demod_input_size)]
		resampled_RDS_A = [np.zeros(demod_input_size*19), np.zeros(demod_input_size*19)]
		resampled_RDS_B = [np.zeros(demod_input_size*19),np.zeros(demod_input_size*19)]
		resampled_RDS_C = [np.zeros(demod_input_size*19//80),np.zeros(demod_input_size*19//80)]
		RRC_RDS_outputs = [np.zeros(demod_input_size*19//80), np.zeros(demod_input_size*19//80)]
		
		# Mixing
		mixed_RDS[0] = RDS_demod_inputs[0][:-1] * RDS_channel
		mixed_RDS[1] = RDS_demod_inputs[1][:-1] * RDS_channel
		

		# Filter pre-rational resampler
		filtered_mixed_RDS[0], RDS_demod_LPF_state[0] = signal.lfilter(RDS_demod_LPF_coeff, 1.0, mixed_RDS[0], zi=RDS_demod_LPF_state[0])
		filtered_mixed_RDS[1], RDS_demod_LPF_state[1] = signal.lfilter(RDS_demod_LPF_coeff, 1.0, mixed_RDS[1], zi=RDS_demod_LPF_state[1])
				
		# Rational resampler
		# Upsample by 19
		for i in range(0, len(filtered_mixed_RDS[0])-1):
			resampled_RDS_A[0][i*19] = filtered_mixed_RDS[0][i]
			resampled_RDS_A[1][i*19] = filtered_mixed_RDS[1][i]

		# Rational resamlping filter
		resampled_RDS_B[0], RDS_sampling_state[0]= signal.lfilter(RDS_sampling_coeff, 1.0, resampled_RDS_A[0], zi=RDS_sampling_state[0])
		resampled_RDS_B[1], RDS_sampling_state[1]= signal.lfilter(RDS_sampling_coeff, 1.0, resampled_RDS_A[1], zi=RDS_sampling_state[1])
		
		resampled_RDS_B[0] = resampled_RDS_B[0]
		resampled_RDS_B[1] = resampled_RDS_B[1]

		#Downsample by 80
		resampled_RDS_C[0] = resampled_RDS_B[0][::80]
		resampled_RDS_C[1] = resampled_RDS_B[1][::80]

		# Apply RRC
		RRC_RDS_outputs[0], RDS_RRC_state[0] = signal.lfilter(RDS_RRC_coeff, 1.0, resampled_RDS_C[0], zi=RDS_RRC_state[0])
		RRC_RDS_outputs[1], RDS_RRC_state[1] = signal.lfilter(RDS_RRC_coeff, 1.0, resampled_RDS_C[1], zi=RDS_RRC_state[1])

		# Data Recovery
		chunk0 = RRC_RDS_outputs[0][:24]
		offset = np.argmax(np.absolute(chunk0))
		recovered_data = [[],[]]
		for i in range(offset, len(RRC_RDS_outputs[0]), 24):
			recovered_data[0].append(RRC_RDS_outputs[0][i])
			recovered_data[1].append(RRC_RDS_outputs[1][i])


		# Debug graph section
		if block_count == 1:
			# # Plot PLL output and input
			# plt.plot(range(len(RDS_demod_inputs[0])), RDS_demod_inputs[0]/max(RDS_demod_inputs[0][500:])*max(RDS_carrier_filtered))
			# plt.plot(range(len(RDS_carrier_filtered)), RDS_carrier_filtered)
			# plt.show()

			# # Plot the mixed signal
			# plt.plot(range(len(mixed_RDS[0])), mixed_RDS[0])
			# plt.show()
			
			# # Plot the PSD of the mixed signal
			# plt.psd(mixed_RDS[0], NFFT=512, Fs=demod_fs)
			# plt.show()

			# # Plot the filtered mixed signal
			# plt.plot(range(len(filtered_RDS_A[0])), filtered_RDS_A[0])
			# plt.show()
			
			# # Plot the upsampled signal
			# plt.plot(range(len(filtered_RDS_B[0])), filtered_RDS_B[0])
			# plt.show()

			# # Plot the output of the resampling filter
			# plt.plot(range(len(filtered_RDS_C[0])), filtered_RDS_C[0])
			# plt.show()

			# Plot the downsampled signal
			plt.plot(range(len(resampled_RDS_C[0])), resampled_RDS_C[0], range(len(RRC_RDS_outputs[0])), np.zeros((len(RRC_RDS_outputs[0]))))
			plt.show()

			# Plot RRC output
			plt.plot(range(len(RRC_RDS_outputs[0])), np.zeros((len(RRC_RDS_outputs[0]))), range(len(RRC_RDS_outputs[0])), RRC_RDS_outputs[0])
			plt.stem(range(offset, len(RRC_RDS_outputs[0]), 24), recovered_data[0])
			plt.show()

		block_count += 1

	print('Finished processing the raw I/Q samples')
	left_channel = np.int16(stereo_out_left*32767)
	right_channel = np.int16(stereo_out_right*32767)
	# write audio data to a .wav file (assumes audio_data samples are -1 to +1)
	wavfile.write("../data/fmStereoBlockRight.wav", int(audio_Fs), right_channel)
	wavfile.write("../data/fmStereoBlockLeft.wav", int(audio_Fs), left_channel)
	wavfile.write("../data/fmMonoBlock-Stereo.wav", int(audio_Fs), np.int16((mono_out/2)*32767))	

	# uncomment assuming you wish to show some plots
	plt.show()