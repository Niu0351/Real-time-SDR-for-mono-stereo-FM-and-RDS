import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math
import fmPllBlock

from fmSupportLib import fmDemodArctan
from fmMonoBasic import lp_impulse, conv

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

demod_Fs = rf_Fs//10

audio_Fs = 48e3
audio_Fc = 16e3
audio_taps = 81
stereo_taps = 181
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
	iq_data = np.subtract((iq_data.astype('float32')), 127.5)
	iq_data = np.divide(iq_data, 127.5)
	# iq_data = iq_data[:len(iq_data)]
	print("Read raw RF data from \"" + in_fname + "\" in float32 format")
	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, \
							 rf_Fc, \
							 window=('hann'), fs=rf_Fs)
	# rf_coeff = lp_impulse(rf_Fc, rf_Fs, rf_taps)

	# coefficients for the filters to extract audio
	audio_coeff = signal.firwin(audio_taps, audio_Fc, window=('hann'), fs=demod_Fs)
	stereo_channel_coeff = bandpass(22e3, 54e3, demod_Fs, stereo_taps)#signal.firwin(stereo_taps, [22e3, 54e3], window=('hann'),  fs=demod_Fs)
	stereo_carrier_coeff = bandpass(18500, 19500, demod_Fs, stereo_taps)#signal.firwin(stereo_taps, 18500, window=('hann'),  fs=demod_Fs)

	# # set up drawing
	# fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)
	# fig.subplots_adjust(hspace = 1.0)

	# select a block_size that is in KB and
	# a multiple of decimation factors
	block_size = 1024 * rf_decim * audio_decim
	block_count = 0

	# states needed for continuity in block processing
	state_i_lpf_state = signal.lfilter_zi(rf_coeff, 1)
	state_q_lpf_state = signal.lfilter_zi(rf_coeff, 1)
	state_pll = {
		"integrator" : 0.0,
		"phaseEst" : 0.0,
		"feedbackI" : 1.0,
		"feedbackQ" : 0.0,
		"ncoOut_buf" : 1.0,
		"trigOffset": 0
	}
	state_phase = 0
	buf_I = 0.0
	buf_Q = 0.0

	# add state as needed for the mono channel filter
	mono_lpf_state = signal.lfilter_zi(audio_coeff, 1.0)
	stereo_channel_state = signal.lfilter_zi(stereo_channel_coeff, 1.0)
	stereo_carrier_state = signal.lfilter_zi(stereo_carrier_coeff, 1.0)
	stereo_processor_state = signal.lfilter_zi(audio_coeff, 1.0)
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
		# print(np.isnan((iq_data[(block_count)*block_size:(block_count+1)*block_size:2])).sum())
		# print("%NaN = " + str(100*np.isnan((iq_data[(block_count)*block_size:(block_count+1)*block_size:2])).sum()/block_size))
		# filter to extract the FM channel (I samples are even, Q samples are odd)
		# filter to extract the FM channel (I samples are even, Q samples are odd)
		i_filt, state_i_lpf_state = signal.lfilter(rf_coeff, 1.0, \
												  iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
												  zi=state_i_lpf_state)
		#
		q_filt, state_q_lpf_state = signal.lfilter(rf_coeff, 1.0, \
												  iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
												  zi=state_q_lpf_state)

		# i_filt, state_i_lpf_100k = conv(iq_data[(block_count)*block_size:(block_count+1)*block_size:2], rf_coeff, state_i_lpf_100k)
		# q_filt, state_q_lpf_100k = conv(iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2], rf_coeff, state_q_lpf_100k)

		# downsample the FM channel
		i_ds = i_filt[::rf_decim]
		q_ds = q_filt[::rf_decim]

		# FM demodulator
		# fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)
		fm_demod, buf_I, buf_Q = fmDemodArctan(i_ds, q_ds, block_count, buf_I, buf_Q)

		# extract the mono audio data through filtering
		# audio_filt = ... change as needed
		mono_filt, mono_lpf_state = signal.lfilter(audio_coeff, 1.0, \
												   fm_demod,
	     											   zi=mono_lpf_state)
		# audio_filt, mono_lpf_100k = conv(fm_demod, audio_coeff, mono_lpf_100k)              

		# downsample audio data
		# audio_block = ... change as needed
		mono_block = mono_filt[::audio_decim]

		stereo_channel, stereo_channel_state = signal.lfilter(stereo_channel_coeff, 1.0, fm_demod, zi=stereo_channel_state)		

		stereo_carrier_filtered, stereo_carrier_state = signal.lfilter(stereo_carrier_coeff, 1.0, fm_demod, zi=stereo_carrier_state)

		stereo_carrier, state_pll = fmPllBlock.fmPll(stereo_carrier_filtered, 19e3, demod_Fs, state_pll, ncoScale=2.0, normBandwidth=0.005)
		
		if block_count == 1:
			# plt.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim))
			# plt.show()
			print("Size of Stereo Carrier After PLL: " + str(np.size(stereo_carrier)))
			print("Size of Stereo Carrier Before PLL: " + str(np.size(stereo_carrier_filtered)))
			pllout = plt.plot(range(0, 300), stereo_carrier[500:800], range(0, 300), stereo_carrier_filtered[500:800])#, range(0, len(fm_demod)), fm_demod)
			plt.show()
			# scatter = plt.scatter(RDS_I, RDS_Q, s=10)
			# plt.show()

		stereo_carrier = stereo_carrier[:-1]

		# mixer
		stereo_mixed = np.zeros(len(stereo_carrier)-1)
		for i in range(0, len(stereo_mixed)-1):
			stereo_mixed[i] = stereo_carrier[i] * stereo_channel[i]

		stereo_processed, stereo_processor_state = signal.lfilter(audio_coeff, 1.0, \
												   stereo_mixed,
												   zi=stereo_processor_state)
		stereo_processed = stereo_processed[::audio_decim]
		stereo_block = [np.zeros(len(stereo_processed)-1), np.zeros(len(stereo_processed)-1)]
		
		# combiner
		for i in range(0, len(stereo_block[0])):
			stereo_block[0][i] = (mono_block[i] + stereo_processed[i])/2
			stereo_block[1][i] = (mono_block[i] - stereo_processed[i])/2

		# concatenate most recently processed audio_block
		# to the previous blocks stored in audio_data
		
		mono_out = np.concatenate((mono_out, mono_block))
		stereo_out_left = np.concatenate((stereo_out_left, stereo_block[0]))
		stereo_out_right = np.concatenate((stereo_out_right, stereo_block[1]))

		# to save runtime select the range of blocks to log iq_data
		# this includes both saving binary files as well plotting PSD
		# below we assume we want to plot for graphs for blocks 10 and 11
		# if block_count >= 10 and block_count < 12:
		# 	# PSD after FM demodulation
		# 	ax0.clear()
		# 	ax0.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
		# 	ax0.set_ylabel('PSD (dB/Hz)')
		# 	ax0.set_xlabel('Freq (kHz)')
		# 	ax0.set_title('Demodulated FM (block ' + str(block_count) + ')')
		# 	# output binary file name (where samples are written from Python)
		# 	fm_demod_fname = "../data/fm_demod_Stereo" + str(block_count) + ".bin"
		# 	# create binary file where each sample is a 32-bit float
		# 	fm_demod.astype('float32').tofile(fm_demod_fname)

		# 	# PSD after extracting mono audio
		# 	ax1.clear()
		# 	ax1.psd(mono_block, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
		# 	ax1.set_ylabel('PSD (dB/Hz)')
		# 	ax1.set_xlabel('Freq (kHz)')
		# 	ax1.set_title('Mono Audio')

		# 	# PSD after decimating mono audio
		# 	ax2.clear()
		# 	ax2.psd(stereo_block[0], NFFT=512, Fs=audio_Fs/1e3)
		# 	ax2.set_ylabel('PSD (dB/Hz)')
		# 	ax2.set_xlabel('Freq (kHz)')
		# 	ax2.set_title('Stereo Audio left')

		# 	ax2.psd(stereo_block[1], NFFT=512, Fs=audio_Fs/1e3)
		# 	ax2.set_ylabel('PSD (dB/Hz)')
		# 	ax2.set_xlabel('Freq (kHz)')
		# 	ax2.set_title('Stereo Audio right')

		# 	# save figure to file
		# 	fig.savefig("../data/fmStereoBlock" + str(block_count) + ".png")

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