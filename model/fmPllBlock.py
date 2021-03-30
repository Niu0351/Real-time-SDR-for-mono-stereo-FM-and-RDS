#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import numpy as np
import math

def fmPll(pllIn, freq, Fs, state,\
		ncoScale = 1.0, phaseAdjust = 0.0, normBandwidth = 0.01):

	"""
	pllIn 	 		array of floats
					input signal to the PLL (assume known frequency)

	freq 			float
					reference frequency to which the PLL locks

	Fs  			float
					sampling rate for the input/output signals

	ncoScale		float
					frequency scale factor for the NCO output

	phaseAdjust		float
					phase adjust to be added to the NCO only

	normBandwidth	float
					normalized bandwidth for the loop filter
					(relative to the sampling rate)

	state 			to be added

	"""
	
	# note: state saving will be needed for block processing


	# scale factors for proportional/integrator terms
	# these scale factors were derived assuming the following:
	# damping factor of 0.707 (1 over square root of 2)
	# there is no oscillator gain and no phase detector gain
	Cp = 2.666
	Ci = 3.555

	# gain for the proportional term
	Kp = (normBandwidth)*Cp
	# gain for the integrator term
	Ki = (normBandwidth*normBandwidth)*Ci

	ncoOut = np.zeros(len(pllIn)+1)
	ncoOut[0] = state["ncoOut_buf"]
	for k in range(len(pllIn)):

		# phase detector
		errorI = pllIn[k] * (+state["feedbackI"])  # complex conjugate of the
		errorQ = pllIn[k] * (-state["feedbackQ"])  # feedback complex exponential

		# four-quadrant arctangent discriminator for phase error detection
		errorD = math.atan2(errorQ, errorI)
		# if (k < 100):
		# 	print("errorD: " + str(errorD))

		# loop filter
		state["integrator"] = state["integrator"] + Ki*errorD

		# update phase estimate
		state["phaseEst"] = state["phaseEst"] + Kp*errorD + state["integrator"]

		# internal oscillator
		trigArg = 2*math.pi*(freq/Fs)*(state["trigOffset"]+k+1) + state["phaseEst"]
		state["feedbackI"] = math.cos(trigArg)
		state["feedbackQ"] = math.sin(trigArg)
		ncoOut[k+1] = math.cos(trigArg*ncoScale + phaseAdjust)

	state["trigOffset"] += len(pllIn)

	state["ncoOut_buf"] = ncoOut[-1]
	# for stereo only the in-phase NCO component should be returned
	return ncoOut, state
	# for RDS add also the quadrature NCO component to the output

if __name__ == "__main__":

	pass
