
%100 MHz 30 kHz SCS, 300ns delay-spread
%fs=122.88e6;
%SCS=30e3;
%NPRB=273;
%DS=300e-9;

function [H] = get_tdl(fs,SCS,PRBS,DS,chan)

  if (strcmp(chan,"tdlc")==1)
  coef_start=1;
  tdl=[1 0 -4.4 
       2 0.2099 -1.2 
       3 0.2219 -3.5
4 0.2329 -5.2 
5 0.2176 -2.5 
6 0.6366 0 
7 0.6448 -2.2 
8 0.6560 -3.9 
9 0.6584 -7.4 
10 0.7935 -7.1 
11 0.8213 -10.7 
12 0.9336 -11.1 
13 1.2285 -5.1 
14 1.3083 -6.8 
15 2.1704 -8.7 
16 2.7105 -13.2 
17 4.2589 -13.9 
18 4.6003 -13.9 
19 5.4902 -15.8 
20 5.6077 -17.1 
21 6.3065 -16 
22 6.6374 -15.7 
23 7.0427 -21.6 
24 8.6523 -22.8];

elseif (strcmp(chan,"tdld")==1)
  coef_start=2;
  tdl=[1 0 -.00147
       2 0.035 -18.8
       3 0.612 -21
       4 1.363 -22.8
       5 1.405 -17.9
       6 1.804 -20.1
       7 2.596 -21.9
       8 1.775 -22.9
       9 4.042 -27.8
       10 7.937 -23.6
       11 9.424 -24.8
       12 9.708 -30.0
	 13 12.525 -27.7];

else
  fprintf("unknown channel %s\n",chan);
  error("exiting");
  end

tdl_norm = sum(10.^(.1*tdl(:,3)));

PRB=1:12;
fgrid=[];
for (PRBs=PRBS)  
  fgrid= [fgrid SCS*((PRBs*12)+PRB)];
end
Htdl=[];

for (f=fgrid),
  Htdl=[Htdl ; (exp(-sqrt(-1)*2*pi*f*tdl(coef_start:end,2)*DS).*sqrt(10.^(.1*tdl(coef_start:end,3)))).'];
end

gv=(sqrt(.5)*randn(1,length(tdl(coef_start:end,2))) + sqrt(-.5)*randn(1,length(tdl(coef_start:end,2))))';

H = Htdl * gv / sqrt(tdl_norm);
if (coef_start == 2) % Ricean channel, so add the mean
  H=H+exp(-sqrt(-1)*2*pi*f*tdl(1,2)*DS).*sqrt(10.^(.1*tdl(1,3))).' / sqrt(tdl_norm);
end

