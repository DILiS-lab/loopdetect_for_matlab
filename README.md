# LoopDetect_for_Matlab

This repository contains a suite of functions to perform feedback loop detection in ordinary differential equation (ODE) models in Matlab. Feedback loops (cycles, circuits) are paths from one node (a variable of the ODE) to itself without visiting any other node twice, and they have important regulatory functions. Together with the loop length it is also reported whether the loop is a positive or a negative feedback loop. An upper limit of the returned number of feedback loops can be entered to limit the runtime (which scales with feedback loop count).

## Quickstart
Download and unzip the contents of this folder. Open Matlab and navigate to the dowloaded folder or add its location to the Matlab path. The following call reports (up to 10) feedback loops for an ODE system determined by a function, here *func_POSm4*, at variable values *s_star* (here these are all equal to 1, column vector). The function values are only allowed to depend on the variable values, other parameters and time have to be fixed.
```
s_star=[1,1,1,1]';
klin=ones(1,8); knonlin=[2.5,3]; t=0;
loop_list=find_loops_vset(@(x)func_POSm4(t,x,klin,knonlin),s_star,10);
%this returns the full list of feedback loops (up to 10)
loop_list{1}
%this returns only the first loop
first_loop=loop_list{1}.loop{1}
```

## Workflow and documentation
A possible workflow and useful command are described in detail in the file *workflow_LoopDetect_Matlab.html* (also as pdf), that relies on the live script workflow_LoopDetect_Matlab.m.

Each function file contains a description with examples of usage which can be called using Matlab's *help* function. The m-files are furthermore documented in in the folder *doc_from_m2html*, documentation generated from m2html (https://www.artefact.tk/software/matlab/m2html/).

## License

All code is licensed under the GNU GPLv3, LoopDetect, LoopDetect_for_Matlab, LoopDetect_for_R, LoopDetect_for_Python, Copyright (C) 2020  Katharina Baum.





