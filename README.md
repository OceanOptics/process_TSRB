# process_TSRB
Code to read the files output from SatCon and produce Rrs spectra and chlorophyll concentrations

This code is a set of functions to process `.dat` files produced by the SatCon software, which processes `.raw` files recorded by a Satlantic HyperPro in buoy mode, a.k.a. a tethered surface radiometer buoy (TSRB)

The code outputs remote-sensing reflectance (Rrs), and chlorophyll concentration calculated with the NASA OC algorithm.

## To get started
* Download the folder of functions, six total:
`calculate_Rrs.m`
`get_ap_oc_from_Rrs.m`
`get_water_iops.m`
`read_datfiles.m`
`read_tiltfile.m`
`run_process_TSRB.m`
  
* The driver code is `run_process_TSRB.m` and the complete input and output details are included in the header of the `.m` file. 
    
Please see the header of `run_process_TSRB.m` for further information on inputs and outputs. 
