% Delta-Sigma Toolbox
% Version 2020.0
% R. Schreier, now retired from Analog Devices
%
% Consult the DSToolbox.pdf and OnePageStory files for more complete documentation.
%
% Modulator synthesis and simulation
%   synthesizeNTF - Noise transfer function (NTF) synthesis.
%   clans         - "Closed-loop analysis of noise shapers"
%		    (NTF synthesis for multi-bit modulators).
%		    Requires the Optimization Toolbox.
%   simulateDSM   - Simulate a delta-sigma modulator using either 
%		    its NTF or its state-space description.
%   simulateSNR   - Simulate a DSM with sine wave inputs of varying amplitudes 
%		    using simulateDSM and record the resulting SNR.
%   simulateMS    - Simulate vector-based mismatch-shaping.
%
% Modulator realization
%   realizeNTF	  - Compute coefficients for one of the supported modulator topologies.
%   stuffABCD	  - Create the state-space description of a modulator given
%                   the coefficients for a particular topology.
%   scaleABCD     - Perform dynamic-range scaling.
%   mapABCD       - Convert a state-space description back to coefficients.
%   calculateTF   - Calculate the NTF and STF from the ABCD matrix.
%   realizeNTF_ct - Realize an NTF with a continuous-time loop filter.
%   mapCtoD	  - map a continuous-time system to discrete-time.
%
% Quadrature functions
%   synthesizeQNTF- Quadrature noise transfer function (NTF) synthesis.
%   simulateQDSM  - Simulate a quadrature delta-sigma modulator using either 
%		    its NTF or its state-space description.
%   realizeQNTF	  - Compute coefficients for one of the supported quadrature modulator topologies.
%   calculateQTF  - Calculate the NTF, STF, INTF and ISTF from the ABCDr matrix.
%   simulateQESL  - Simulate the element selection logic in a quadrature differential mismatch-shaping DAC.
%
% Other functions related to delta-sigma 
%   infnorm       - Calculate the infinity norm (maximum gain) of a DT TF.
%   predictSNR    - SNR prediction for binary modulators 
%		    (Uses the describing function method of Ardalan & Paulos)
%   designHBF     - Design a multiplierless half-band filter which uses
%                   Saramaki's hierarchical filter structure.
%   designPBF     - Design a symmetric polynomial-based filter using
%                   Hunter's method
%   findPIS       - Compute a positively-invariant set for a DSM. (The 
%                   PosInvSet sub-directory will need to be added to your PATH)
%
% Demonstrations and Examples
%   dsdemo1       - Synthesize a 5th-order lowpass and an 8th-order bandpass NTF.
%   dsdemo2       - Time-domain simulation and SNR calculation.
%   dsdemo3       - Modulator realization and dynamic range scaling.
%   dsdemo4       - Audio demonstration of MOD1 and MOD2
%   dsdemo5	  - Simulate the element selection logic of a mismatch-shaping DAC.
%   dsdemo6	  - Design a hardware-efficient halfband filter.
%   dsdemo7	  - Find positively-invariant sets for second-order (and third-order) modulators.
%   dsdemo8       - Continuous-time bandpass modulator design using LC tanks.
%   dsexample1	  - Discrete-time lowpass/bandpass/quadrature modulator.
%   dsexample2	  - Continuous-time lowpass modulator.
%
% Copyright (c) 1993-2020 R. Schreier
