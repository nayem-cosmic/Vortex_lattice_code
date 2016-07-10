#!/bin/bash
echo Compiling...
gfortran -o vlm_part1.out vlm_part1.F90
gfortran -o vlm_part2.out vlm_part2.F90

echo Part 1
./vlm_part1.out
echo Do you want to continue with part 2? '(y/n)'
echo N.B. Part 2 contains angles vs. coefficients graph calculations etc.
read  ans
if [ $ans == 'y' ];then
    echo Part 2
    ./vlm_part2.out
fi
    
echo 
echo Graph Generating...
python2 vlm_mesh.py &
python2 dp_coeff.py &
python2 d_lift_span.py &
python2 d_lift_chord.py &
python2 d_drag_span.py &
python2 d_drag_chord.py &
if [ $ans == 'y' ];then
    python2 alphacl.py &
    python2 alphacm.py &
    python2 alphacd.py &
    python2 lambdacl.py &
    python2 lambdacm.py &
    python2 lambdacd.py &
    python2 phicl.py &
    python2 phicm.py &
    python2 phicd.py &
    python2 chcl.py &
    python2 chcd.py &
fi

echo All done!
echo Thank you...
