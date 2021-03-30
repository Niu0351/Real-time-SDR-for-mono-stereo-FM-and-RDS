#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy import signal
import math

# use "custom" fmDemodArctan
from fmSupportLib import fmDemodArctan

# the radio-frequency (RF) sampling rate
# this sampling rate is either configured on RF hardware
# or documented when a raw file with IQ samples is provided
rf_Fs = 2.4e6

# the cutoff frequency to extract the FM channel from raw IQ data
rf_Fc = 100e3

# the number of taps for the low-pass filter to extract the FM channel
# this default value for the width of the impulse response should be changed
# depending on some target objectives, like the width of the transition band
# and/or the minimum expected attenuation from the pass to the stop band
rf_taps = 151

# the decimation rate when reducing the front end sampling rate (i.e., RF)
# to a smaller samping rate at the intermediate frequency (IF) where
# the demodulated data will be split into the mono/stereo/radio data channels
rf_decim = 10

# audio sampling rate (we assume audio will be at 48 KSamples/sec)
audio_Fs = 48e3

# complete your own settings for the mono channel
# (cutoff freq, audio taps, decimation rate, ...)
audio_Fc = 16e3
audio_taps = 151
audio_decim = 5

# own impulse function
def lp_impulse(Fc,Fs,N_taps):

	h = np.zeros(N_taps)
	norm_cut = Fc/(Fs/2)

	for i in range(0,N_taps):
		if i == (N_taps-1)/2:
			h[i] = norm_cut
		else:
			h[i] = norm_cut*((np.sin(math.pi*norm_cut*(i-(N_taps-1)/2)))/
							 (math.pi*norm_cut*(i-(N_taps-1)/2)))

		h[i] = h[i] * (np.sin(i*math.pi/N_taps))**2

	return h

# Own convolution
def conv(x,h,buf_wave):

	#create a new array for storing the previous blocks info
	# buf_wave = np.zeros(int(len(h)-1))
	filtered_data = np.zeros(shape=x.shape)

	# while True:
	cycle = 0

	for n in range(0,len(x)):
		buf = 0

		for k in range(0,int(len(h))):

			#to check if the sample need to be computed within previous block info
			if (cycle+buf) < (int(len(h)-1)):
				# filtered_data[position+n] = filtered_data[position+n] + buf_wave[cycle+buf]*h[int(len(h))-1-k]
				filtered_data[n] += buf_wave[cycle+buf]*h[int(len(h))-1-k]
				buf = buf + 1

			else:
				# filtered_data[position+n] = filtered_data[position+n] + x[position+n+k-int(len(h))+1]*h[int(len(h))-1-k]
				filtered_data[n] += x[n+k-int(len(h))+1]*h[int(len(h))-1-k]

		cycle += 1

	#store the new info from this block and would be used in the next block
	buf_wave = x[(int(len(x)) - int(len(h)) +1):]

	return filtered_data, buf_wave

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is normalized between -1 and +1 and interleaved
	in_fname = "../data/test2.raw"
	iq_data = np.fromfile(in_fname, dtype='float32')
	print("Read raw RF data from \"" + in_fname + "\" in float32 format")

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

	# FM demodulator (check the library)
	fm_demod, buf_I, buf_Q = fmDemodArctan(i_ds, q_ds)
	# we use a dummy because there is no state for this single-pass model

	# set up drawing
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)
	fig.subplots_adjust(hspace = 1.0)

	# PSD after FM demodulation
	ax0.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
	ax0.set_ylabel('PSD (db/Hz)')
	ax0.set_title('Demodulated FM')

	# coefficients for the filter to extract mono audio
	# audio_coeff = np.array([]) # to be updated by you during in-lab
	audio_coeff = signal.firwin(audio_taps, audio_Fc/(24e4/2), window=('hann'))
	# audio_coeff = lp_impulse(audio_Fc, 24e4, audio_taps)

	# extract the mono audio data through filtering
	# audio_filt = np.array([]) # to be updated by you during in-lab
	audio_filt = signal.lfilter(audio_coeff, 1.0, fm_demod)
	# audio_filt = conv(fm_demod,audio_coeff,block_size)

	# you should uncomment the plots below once you have processed the data

	# PSD after extracting mono audio
	ax1.psd(audio_filt, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
	ax1.set_ylabel('PSD (db/Hz)')
	ax1.set_title('Extracted Mono')

	# downsample audio data
	# audio_data = np.array([]) # to be updated by you during in-lab
	audio_data = audio_filt[::audio_decim]

	# PSD after decimating mono audio
	ax2.psd(audio_data, NFFT=512, Fs=audio_Fs/1e3)
	ax2.set_ylabel('PSD (db/Hz)')
	ax2.set_title('Mono Audio')

	# save PSD plots
	fig.savefig("../data/fmMonoBasic.png")
	plt.show()

	# write audio data to file (assumes audio_data samples are -1 to +1)
	wavfile.write("../data/fmMonoBasic.wav", int(audio_Fs), np.int16((audio_data/2)*32767))
# during FM transmission audio samples in the mono channel will contain
# the sum of the left and right audio channels; hence, we first
# divide by two the audio sample value and then we rescale to fit
# in the range offered by 16-bit signed int representation
