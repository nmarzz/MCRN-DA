# Models README
Please update this read me as goals evolve

Goals:
 - Get shallow water model working with LETKF

Notes:
- 2019 11 13: copied create_analysis.m, letkf.m, letkfstep.m, mgs.m, winobs.m from the LETKF files from August 1st and wrote them into this directory.
 Added drivL96.m as drivedit.m and initialize.m as initedit.m
- 2019 11 25: 
   Copied init_swe.m (a SWE file) and added it to this directory.
   Copied formod.m (a SWE file) and added it to this directory.
   Copied lax_wendroff.m (a SWE file) and added it to this directory.
   Moved the creation of the "truth" state to the driver file (it happens when init_swe is called).
   Replaced initialized.m with createbackground.m where createbackground turns the "truth" into an ensemble by adding noise.
   In letkfstep.m, replaced rkfixed with a loop calling formod for each member of the ensemble
   Calling letkfstep.m is super slow...
   In runnew.m replaced rkfixed call with a formod call
