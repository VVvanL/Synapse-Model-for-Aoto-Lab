mex -O synapse_C.c

mex -O -DABSORB_AT_GLIAL_SHEATH   -output synapse_absorb_at_glia_C    synapse_C.c

mex -O -DABSORB_AT_CLEFT_BOUNDARY -output synapse_absorb_at_cleftbd_C synapse_C.c
