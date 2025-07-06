% Template of the matlab script file for collecting some of the 
% results of your written exam (Modern Design of Control Systems)
% 
s=tf('s');
%
% Below type your Surname (family name) and your ID number 
%
nome=['Ayanmanesh Motlaghmofrad-StudentID'];
%
% Below type the chosen nominal value of the uncertain variable k
%
k_n=  ;
%
% Below type the chosen nominal value of the uncertain variable p1
%
p1_n=    ;
%
% Below type the chosen nominal value of the uncertain variable p2
%
p2_n=    ;
%
% In the following lines, type the requested weighting functions, denoting 
% 's' as the independent complex variable. A unitary default value, 
% e.g. WS=tf(1,1), is typed to indicate that you do not have a result for
% that item. Below type the weighting function WS
%
WS=tf(1,1);
%
% Below type the weighting function WT
%
WT=tf(1,1);
%
% Below type the weighting function Wu
%
Wu=tf(1,1);
%
% Below type the weighting function W1
%
W1=tf(1,1);
%
% Below type the weighting function W2
%
W2=tf(1,1);
%
% Below type the weighting function W1mod
%
W1mod=tf(1,1);
%
% Below type the weighting function W2mod
%
W2mod=tf(1,1);
%
% Below type other possible weighting functions. 
%
WSmu=tf(1,1);
%
WTmu=tf(1,1);
%
Wunew=tf(1,1);
%
% Below type the transfer function Gcmod of the designed controller as 
% obtained from the optimizer
%
Gcmod=tf(1,1);
%
% Below type the transfer function Gc of the designed controller in its 
% final form
%
Gc=tf(1,1);
%
% In the lines below, write obtained nominal performance in the time domain
% "0" (zero) is the default value if you do not have a result
% Below type the obtained rise time
%
r_tr=0             ;
%
% Below type the obtained settling time
%
r_ts=0             ;
%
% Below type the obtained percentage overshoot
%
r_sovr=0             ;
%
% Below type the maximum value of the absolute output error in the presence
% of da
%
r_yda=0             ;
%
% Below type the maximum value of the absolute output error in the presence
% of dp
%
r_ydp=0             ;
%
% Below type the amplitude of the steady-state output error in the presence
% of dp
%
r_ydp_ss=0             ;
%
% Below type the amplitude of the steady-state output error in the presence
% of ds
%
r_yds_ss=0             ;
%
% Below type the maximum value of the Lower Bound of mu for Robust
% Performance:
%
muRPLBmax=0;
%
%
% Below type the maximum value of the Upper Bound of mu for Robust
% Performance:
%
muRPUBmax=0;
%
%
% Below type the maximum value of the Lower Bound of mu for Robust
% Stability:
%
muRSLBmax=0;
%
%
% Below type the maximum value of the Upper Bound of mu for Robust
% Stability:
%
muRSUBmax=0;
%
% Do not modify the lines below
%
save Results_MDCS.mat  k_n p1_n p2_n WS WT Wu W1 W2 W1mod W2mod WSmu WTmu Wunew ...
     Gcmod Gc r_tr r_ts r_sovr r_yda r_ydp r_yds_ss r_ydp_ss ...
     muRPLBmax muRSLBmax muRPUBmax muRSUBmax
%
% Fill in and run this script file, in order to get the related .mat file
% 
% This script file and the related .mat file must be placed in
% sub-directory MCDS\Results
%
%