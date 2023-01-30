x0=zeros(1,127);
x1=x0;
x0(1:7)=[1 0 0 0 0 0 0];
x1(1:7)=[1 0 0 0 0 0 0];

for (n=1:120)
  x0(n+7)=rem(x0(n+4)+x0(n),2);
  x1(n+7)=rem(x1(n+1)+x1(n),2);
end

d = zeros(336*3,127);
for (nid2=0:2)
  for (nid1=0:335)
    m0=15*floor(nid1/112) + (5*nid2);
    m1=rem(nid1,112);
    for (n=0:126)
      d(1+nid2+(nid1*3),1+n) = (1-2*x0(1+rem(n+m0,127)))*(1-2*x1(1+rem(n+m1,127))); 
    end
  end
end
