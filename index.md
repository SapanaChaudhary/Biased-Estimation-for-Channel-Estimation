### Biased Estimator for OFDMA

With Vaishnavi Adella (now at Qualcomm, India), Sai Charan Thoutam (now at Qualcomm, India)

Implemented the [paper](http://ieeexplore.ieee.org/document/6471293/) : Biased estimators with adaptive shrinkage targets for orthogonal frequency division multiple access channel estimation

### Channel estimation : What and Why of it 

* Channel Estimation is the means of characterising channel effects-scattering, fading and power decay
* Can be done using either decision feedback scheme or the known pilot symbols
* Here : Using user specific pilots
* Why is it crucial ?
    Huge transmit power spent on it 
    Incorrect estimates lead to residual cancellation errors
    Necessary for high data rates

### The unbiased estimators and the scope for biased estimator

* Channel statistics – known 
* Estimation methods such as modified least squares require wide band pilots
* Finite or large number of pilots / wide band pilots
* The two-dimensional (2D)-minimum mean square error (2D-MMSE) methods [1, 2] can be applied using the pilots in the RB. However, optimal MMSE estimator requires the knowledge of the channel statistics which are seldom known accurately at the receiver.
* MVUE-Zero bias as constraint and acheive CRLB
* ML estimate is biased in order to reduce the MSE.

### Motivation for JS estimator

* 2D-MMSE can be used when channel statistics is unknown
* Assumption: Ideally band limited and time limited uniform scattering function
Is that even possible ? Google says 'the only time and band limited signal is zero'. A little more of research and we find: A time limited signal with most of its energy contained in the band is a good approximation for both time and band limited signal. 
* 2D-MMSE : Performance degrades if the robust filter has finite number of taps [2]
* JS-estimator used to bridge the gap between robust 2D MMSE and optimal MMSE

### Main contribution of the work

* Applying biased estimation techniques for localised CFR estimation
* Adaptively choosing a vector shrinkage target for the biased estimator that best reflects the time–frequency selectivity of the CFR using the aid of multiple hypothesis tests 
* Choosing the thresholds for the hypothesis tests such that the probability of incorrect choice of shrinkage targets is bounded.

