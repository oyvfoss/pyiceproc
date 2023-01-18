TO DO
------

Sea ice draft
''''''''''''''

- Function to calculate draft from depth/sound speed/tilt/altimeter dist/ow ss correction
  - Return data, not an updated array (need to do this recursively for OW correction..) 
  - Call the raw data (*surface_position* or something (not *draft*))
  - Select only the *is_ice* portion :math:`\rightarrow` this is the draft  
     - Need to think about this! Do we actually only want draft measurements from full (100% SIC)
       ensembles? (I *think* we want to only use all-4-beam FOM samples - then do ensemble statistics 
       that dont include not-ice samples). Then out draft **only contains measurements from times where we
       are sure that we are measuring sea ice** - this seems like the right way to go.  

  - Function for cleaning draft data (crazy outliers, etc) - do this before ensemble averaging.
   - Include rejection based on quality criterion!

  - Assign the resulting distance and surface position to the dataset (maybe as ``DRAFT_INITIAL`` or similar?)
    - Want to run OW script below and recompute

- Function to compute a long-term OW mean (to be used for OW correction)
  - Want to use ``scipy.uniform_filter1d`` or similar for running stats (avoid the more uncommon 
    libraries I normally use..) 
  - Need visualization/evaluation and relatively easy customization.
    - Tink of which parameter is useful here :math:`\rightarrow` Window size of runnig mean?    
  - Export the appropriate correction factor..

- Wrapper function to apply OW correction! (resulting in "DRAFT", "distance" fields..)

Ice velocity
''''''''''''
- Mask each SAMPLE by FOM criterion
- Basic outlier and parameter editing on SAMPLE data.
- Ensemble statistics (median?) - rejecting invalid ensembles.
- Magdec correction (can be applied after ensemble averaging).
- Consider a sound speed correction (but drop if it is comlpetely negligible anyways..)
- Save variables to xr dataset

Ocean velocity
''''''''''''''
- Compute bin depths from (updated) depth.
  - I don't think it's worth doing sound speed/density corrections to bin depths/velocity magnitude.
 
- Mask by sidelobe interference and above water measurements.
 - Account for ice draft here! (Probably the easiest approach.) 
 - Reject bins with nothing (<0.1%) in them.
  
- Mask SAMPLES based on the normal, successive threshold criteria
- Ensemble statistics (median?) - rejecting invalid ensembles.
- Magdec correction (can be applied after ensemble averaging).
- Save variables to xr dataset.

Other functionality
'''''''''''''''''''

- Chopping start/end (I think I have good starting point code for this somewhere)
- Basic visualisation/statistics (no need to go overboard)
- Export to reduced xr Dataset (drop non-useful stuff)
- Export to nc
 - For analysis - but maybe eventually also for publication with appropriate conventions.. 

- Interpolation of ocean velocity onto fixed depth (do this at the ensemble stage, 
  cumbersome and not too useful to do it to SAMPLEs)

Documentation
''''''''''''''

- Look over README.md.
- Make working example in README.md
- Make more detailed/realistic example (notebook) 
- Look into sphinx automatic documentation from functions..

Maybe
'''''

- Look into burst stuff
- Time means etc
- Packaging for install