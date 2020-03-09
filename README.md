# muon_hits

Calculate total number of hits per MAPMT for the SciCRT detector. For this purpose we need to fit a parabola to the
pedestal data (per file) and then calculate the number of counts per scintillator bar. Because we use 4-fold coincidence
for detecting muon in the telescope (2 MAPMTs X and Y top and 2 MAPMTs X and Y bottom), the sum of of two MAPMTs hits
should be ~ Total events.

The threshold is fixed to 290 ADC.

The softwares requires numpy and scipy to work.
