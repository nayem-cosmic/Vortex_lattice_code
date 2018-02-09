#!/bin/bash

rm figures/*
echo All previous figures deleted.

echo Compiling...
gfortran -o vlm_part1.out vlm_part1.F90
gfortran -o vlm_part2.out vlm_part2.F90

echo Part 1
./vlm_part1.out
echo Part 2
./vlm_part2.out

while true;do
    echo Graph will be generated of:
    echo -e '\t''(1)'  Mesh
    echo -e '\t''(2)'  Pressure Coefficients at Different Places
    echo -e '\t''(3)'  Lift vs. Span
    echo -e '\t''(4)'  Lift vs. Chord
    echo -e '\t''(5)'  Drag vs. Span
    echo -e '\t''(6)'  Drag vs. Chord
    echo -e '\t''(7)'  CL, CM, CD vs. Alpha
    echo -e '\t''(8)'  CL, CM, CD vs. Lambda
    echo -e '\t''(9)'  CL, CM, CD vs. Phi
    echo -e '\t''(10)' CL, CD vs . Height Above Ground
    echo -e '\t''(11)' Generate All Graphs
    echo -e '\t''(0)' Exit
    read  ans
    if [ $ans == '1' ];then
        python2 vlm_mesh.py &
    elif [ $ans == '2' ];then
        python2 dp_coeff.py &
    elif [ $ans == '3' ];then
        python2 d_lift_span.py &
    elif [ $ans == '4' ];then
        python2 d_lift_chord.py &
    elif [ $ans == '5' ];then
        python2 d_drag_span.py &
    elif [ $ans == '6' ];then
        python2 d_drag_chord.py &
    elif [ $ans == '7' ];then
        python2 alphacl.py &
        python2 alphacm.py &
        python2 alphacd.py &
    elif [ $ans == '8' ];then
        python2 lambdacl.py &
        python2 lambdacm.py &
        python2 lambdacd.py &
    elif [ $ans == '9' ];then
        python2 phicl.py &
        python2 phicm.py &
        python2 phicd.py &
    elif [ $ans == '10' ];then
        python2 chcl.py &
        python2 chcd.py &
    elif [ $ans == '11' ];then
        python2 vlm_mesh.py &
        python2 dp_coeff.py &
        python2 d_lift_span.py &
        python2 d_lift_chord.py &
        python2 d_drag_span.py &
        python2 d_drag_chord.py &
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
    elif [ $ans == '0' ];then
        break
    else
        echo Please enter correct option!
        echo 
        continue
    fi
done

echo All done!
echo Thank you...
