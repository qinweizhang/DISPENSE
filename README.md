# SoSAIDPRESS
Stack of Spiral AIded Diffusion PREpared Stimulated-echo Sequence

**This repo describes all the details for the Stack of Spiral AIded Diffusion PREpared Stimulated-echo Sequence (SoSAIDPRESS) implementation.**

---

* Add SQDPNav and AQDPNav to test. Not finish yet. GR object is not defined yet.

* Set AQDPNav:dest to DAS only. AQDPNav is just a copy of AQbase. Problem not solved yet.

* Fixed the AQDPNav crush problem. Reason: AQDPNav enabled during receiver gain. Unexpected data passed for receiver gain calculation and crushes the scan. Now AQDPNav is disabled during receiver gain phase.
AQDPNav and GRDPNav need to be further calculated.

* Add DP-GRdr mechanism to vary dephasing and rephasing gradient for different shots. Trying to eliminate the influence from previous shots. V4.01

* Change to another scheme for flipping DP-GRr-d. Now the gradient rotate in the P-S plane with a cycle of 4 shots.  V4.02 Kerry

* Declared new GR object for DP Spiral Navigator.
SQDPNav is defined in using existing function.
The existing function is modified for DP Navigator purpose.
Remove AQDPNav and SQDPNav code in PDF. They are handled in the mpirfe_sq__g.c now.
To do: only one spiral is defined now. GR waveform needs to be extended for 8 spiral acquisition (continuously).
Gradient in S direction needs to be added. V4.1

* Add slice encoding gradients. And introduced 3 parameter to EC.  

* Set AQDPNav:dest = RECON_DAS instead of DAS. DAS data will not end up in reconstructor, therefore not be exported.
Also, proper RC parameters need to be set according to the samples from AQDPNav. Not solve yet.

* Add spoiler in the SQDPNav; Fixed spiral in gradient error.

* Fixed the AQDPNav sample size problem. Maarten modified the recon parameter and generated two exe files. Now  RC_ad_hoc_int_array[ 0 ] and [1] are used to control the parameter and pass the correct sample size to Recon. V4.2

* Add a new MP parameter to control phase cycling pattern for RFt2prep[last]. It's possible to do phase less cycling now. But need lots of SDE editing.

* Also synced the GR repetition option for spiral gradient. Too long single gradient waveform will cause crush. V4.4.6

* Set AQ::DPNav->dest to RECON. Sending data to DAS will crush the scan.
Set slab selective RF::t2prep->freq_step1 = GR::s->str for accurate slice off-center. V4.4.7

* Test push from Atom
