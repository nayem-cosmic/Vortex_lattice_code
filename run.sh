#!/bin/sh
echo compiling...

gfortran vlm_v2.F90
./a.out

python2 vlm_mesh.py &
python2 dp_coeff.py &
python2 d_lift_span.py &
python2 d_lift_chord.py &
python2 d_drag_span.py &
python2 d_drag_chord.py &

echo done!
