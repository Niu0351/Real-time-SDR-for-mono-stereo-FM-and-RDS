import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math
import fmPll

from fmSupportLib import fmDemodArctan
from fmMonoBasic import lp_impulse, conv

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

audio_Fs = 48e3
audio_Fc = 16e3
audio_taps = 151
stereo_taps = 151
audio_decim = 5

def stereo_extraction(fm_demod):
	bandpass_coeff = signal.firwin(stereo_taps, [22e3/(audio_Fs/2), 54e3/(audio_Fs/2)], window=('hann'), fs=audio_Fs)
	print("size of message bandpass coeff: " + str(len(bandpass_coeff)))
	filtered_stero = signal.lfilter(fm_demod, 1, bandpass_coeff)
	print("size of filtered message: " + str(len(filtered_stero)))

	return filtered_stero


def stereo_carrier_recovery(fm_demod):
	bandpass_coeff = signal.firwin(stereo_taps, [18500/(audio_Fs/2), 19500/(audio_Fs/2)], window=('hann'), fs=audio_Fs)
	print("size of bandpass carrier coeff: " + str(len(bandpass_coeff)))
	filtered_demod = signal.lfilter(fm_demod, 1, bandpass_coeff)
	print("size of filtered carrier: " + str(len(filtered_demod)))
	clocked_output = fmPll.fmPll(filtered_demod, 19000, rf_Fs)
	print("size of pll out: " + str(len(clocked_output)))

	return clocked_output[:-1]


def stereo_mix(stereo_carrier, stereo_channel):
	return np.multiply(stereo_carrier, stereo_channel)


def digital_filter_and_downsample(full_scale_signal):
	lp_coeff = signal.firwin(audio_taps, 16000/(audio_Fs/2), window=('hann'))
	
	# In c++ the filter and downsample will be combined to reduce comuptations required
	filtered_signal = signal.lfilter(full_scale_signal, 1, lp_coeff)

	return filtered_signal[::audio_decim]


def build_stereo(stereo, mono):
	return [np.add(stereo, mono), np.subtract(mono, stereo)]


def rf_frontend():
	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))
	#  rf_coeff = lp_impulse(rf_Fc, rf_Fs, rf_taps)

	# filter to extract the FM channel (I samples are even, Q samples are odd)
	i_filt = signal.lfilter(rf_coeff, 1.0, iq_data[0::2])
	q_filt = signal.lfilter(rf_coeff, 1.0, iq_data[1::2])
	
	block_size = int(len(iq_data[0::2]))
	# print(block_size)
	# i_filt, buf_wave = conv(iq_data[0::2], rf_coeff, np.zeros(len(rf_coeff)-1))
	# q_filt, buf_wave = conv(iq_data[1::2], rf_coeff, np.zeros(len(rf_coeff)-1))

	# downsample the FM channel
	i_ds = i_filt[::rf_decim]
	q_ds = q_filt[::rf_decim]
	print("ARCTAN DEMOD BEGIN (LEN I = " + str(len(i_ds)) + ")")
	# FM demodulator (check the library)
	fm_demod = fmDemodArctan(i_ds, q_ds)
	print("ARCTAN DEMOD EMD")
	return fm_demod
	# we use a dummy because there is no state for this single-pass model



if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is normalized between -1 and +1 and interleaved
	in_fname = "../data/test2.raw" #in_fname = "../data/test2.raw"
	iq_data = np.fromfile(in_fname, dtype='float32')
	iq_data = iq_data[:len(iq_data)//3]
	print("Read raw RF data from \"" + in_fname + "\" in float32 format")
	(fm_demod, prev_phase) = rf_frontend()

	print("RF FRONTEND COMPLETE")

	# set up drawing
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)
	fig.subplots_adjust(hspace = 1.0)

	# PSD after FM demodulation
	ax0.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
	ax0.set_ylabel('PSD (db/Hz)')
	ax0.set_title('Demodulated FM')

#----------------------------------------------------------------------------------------------#

	mono_final = digital_filter_and_downsample(fm_demod)
	
	print("MONO COMPLETE")

	stereo_channel = stereo_extraction(fm_demod)
	print("MESSAGE EXTRACTION")
	stereo_carrier = stereo_carrier_recovery(fm_demod)
	print("CARRIER EXTRACTION COMPLETE")
	mixed_stereo = stereo_mix(stereo_carrier, stereo_channel)
	print("STEREO MIXING COMPLETE")
	
	stereo_filtered = digital_filter_and_downsample(mixed_stereo)
	
	print("FINAL STEREO FILTERING COMPLETE")

	stereo_final = build_stereo(stereo_filtered, mono_final)
	
	print("FINAL AUDIO COMPLETE COMPLETE")



	ax1.psd(mono_final, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
	ax1.set_ylabel('PSD (db/Hz)')
	ax1.set_title('Mono Audio')

	# PSD after decimating mono audio
	ax2.psd(stereo_final[0], NFFT=512, Fs=audio_Fs/1e3)
	ax2.set_ylabel('PSD (db/Hz)')
	ax2.set_title('Stereo Audio 0')

	ax2.psd(stereo_final[1], NFFT=512, Fs=audio_Fs/1e3)
	ax2.set_ylabel('PSD (db/Hz)')
	ax2.set_title('Stereo Audio 1')

	# save PSD plots
	fig.savefig("../data/fmStereoBasic.png") #fig.savefig("../data/fmStereoBasic.png")
	plt.show()

	# write audio data to file (assumes audio_data samples are -1 to +1)
	wavfile.write("../data/fmStereoBasic.wav", int(audio_Fs), np.int16((stereo_final)*32767)) #wavfile.write("../data/fmStereoBasic.wav", int(audio_Fs), np.int16((stereo_final/2)*32767))
# during FM transmission audio samples in the mono channel will contain
# the sum of the left and right audio channels; hence, we first
# divide by two the audio sample value and then we rescale to fit
# in the range offered by 16-bit signed int representation
