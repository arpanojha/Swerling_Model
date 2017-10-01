clc;
%gaussian rcs
x=0:0.5:180;
n=361;
for i=1:n
    x1=x(i)*x(i);
end
x0=sqrt(x1/n);%rms
me=mean(x);   %mean
obj = gmdistribution(0.1,1,0.1);  %gaussian noise
k=random(obj,n);      
k=abs(k/max(k));
k=3.1.*k;
% netrcs=netrcs/max(netrcs);
% k=netrcs+k';

rcs1(1)=mean(k(1:60));
rcs1(2)=mean(k(61:120));
rcs1(3)=mean(k(121:180));
rcs1(4)=mean(k(181:240));
rcs1(5)=mean(k(241:300));
rcs1(6)=mean(k(301:361));


%% 
%swerling 3
nfa=1e9;
np=100;
snrbar=rcs1;
%snrbar = 10.0.^(snrbar/10.);
eps = 0.00000001;
delmax = .00001;
delta =10000.;
% Calculate the threshold Vt
pfa = np * log(2) / nfa;
sqrtpfa = sqrt(-log10(pfa));
sqrtnp = sqrt(np);
vt0 = np - sqrtnp + 2.3 * sqrtpfa * (sqrtpfa + sqrtnp - 1.0);
vt = vt0;
Vprev=vt0;
er=0.000000001;
counter=0;
while er < abs(Vprev/10000)
    igf = gammainc(Vprev,np);
    g=(0.5^(np/nfa))-igf;
    gprime=exp(-Vprev)*Vprev^(np-1);
    gprime= gprime/factorial(np-1);
    Vt=Vprev-(g/gprime);
    er = abs(Vt-Vprev);
    Vprev=Vt;
    counter=counter+1;
end
vt=Vt;
temp1 = vt ./ gadd(1.0 , 0.5 .* np .*snrbar);
temp2 =gadd( 1.0 , 2.0 ./ (np .* snrbar));
temp3 = 2.0 .* (np - 2.0) ./ (np .* snrbar);
ko = exp(-temp1) .* temp2.^(np-2.) .* (1.0 + temp1 - temp3);
if (np <= 2)
pd3 = ko;
return
else
temp4 = vt.^(np-1.) .* exp(-vt) ./ (temp1 .* factorial(np-2.));
temp5 = vt ./ (1.0 + 2.0 ./ (np .*snrbar));
pd3 = temp4 + 1.0 - gammainc(vt,np-1.) + ko .*gammainc(temp5,np-1.);
end
pd3=max(pd3,0.);
% plot(angle,pd3);
% grid;
% xlabel ('Aspect angle-degrees');
%  ylabel ('Pd');
 %%
 %swerling 2 
SNR=snrbar;

vt0 = np - sqrt(np)+(2.3*sqrt(-log10(pfa))*(sqrt(-log10(pfa))+sqrt(np)-1));
Vprev=vt0;
er=0.000000001;
counter=0;
while er < abs(Vprev/10000)
    igf = gammainc(Vprev,np);
    g=(0.5^(np/nfa))-igf;
    gprime=exp(-Vprev)*Vprev^(np-1);
    gprime= gprime/factorial(np-1);
    Vt=Vprev-(g/gprime);
    er = abs(Vt-Vprev);
    Vprev=Vt;
    counter=counter+1;
end

vt=Vt;

if np <=50
     a2=(vt/1+SNR) ;
     gamma1=gammainc(a2,np);
     pd2= 1-gamma1;
else
    a1=gadd(vt,-np.*(gadd(1,SNR)));
    omega = sqrt(np).*gadd(1,SNR);
    V=gdivide(a1,omega);
    b1=erfc(V/sqrt(2))/2;
    b2=exp(-(V.^2)/2)/sqrt(2*pi);
    c3= -1/(3*sqrt(np));
    c4= 1/(4*np);
   
    c6=(c3^2)/2;
    b3=c3.*gadd(V.^2,-1);
    b4=c4.*gadd((3.*V),-(V.^3));
    b5=-c6.*(gadd(V.^5,-10).*gadd((V.^3),(15*V)));
    bnew= b3+b4+b5;
    bnew=bnew.*b2;
    pd2=gadd(b1,-bnew);
end
pd2=min(pd2,1.);
% plot(angle,pd2);
%  grid;
%  xlabel ('Aspect angle-degrees');
%  ylabel ('Pd');

%%
%swerling1
vt0 = np- sqrtnp+2.3*sqrtpfa*(sqrtpfa+sqrtnp-1);
vt=vt0;

Vprev=vt0;
er=0.000000001;
counter=0;
while er < abs(Vprev/10000)
    igf = gammainc(Vprev,np);
    g=(0.5^(np/nfa))-igf;
    gprime=exp(-Vprev)*Vprev^(np-1);
    gprime= gprime/factorial(np-1);
    Vt=Vprev-(g/gprime);
    er = abs(Vt-Vprev);
    Vprev=Vt;
    counter=counter+1;
end
vt=Vt;
if np ==1
    temp= gdivide(-vt,gadd(1,snrbar));
    pd=exp(temp);
    
end
temp1=1+np*snrbar;
temp2=gdivide(1,(np*snrbar));
temp = 1+temp2;
val1=temp.^(np-1);
igf1=gammainc(vt,np-1);
v0=vt(1);
igf2=gammainc(gdivide(vt(1),temp),np-1);
pd1 = 1 + gadd(-igf1(1),gmultiply(val1,gmultiply(igf2,exp(-gdivide(vt(1),temp1)))));
%%
%swerling4

vt0 = np - sqrtnp + 2.3 * sqrtpfa * (sqrtpfa + sqrtnp - 1.0);
vt = vt0;
Vprev=vt0;
er=0.000000001;
counter=0;
while er < abs(Vprev/10000)
    igf = gammainc(Vprev,np);
    g=(0.5^(np/nfa))-igf;
    gprime=exp(-Vprev)*Vprev^(np-1);
    gprime= gprime/factorial(np-1);
    Vt=Vprev-(g/gprime);
    er = abs(Vt-Vprev);
    Vprev=Vt;
    counter=counter+1;
end
vt=Vt;
h8 = snrbar /2.0;
beta = 1.0 + h8;
beta2 = (2.0 .* beta.^2) - 1.0;
beta3 = 2.0 .* beta.^3;
if (np >= 50)
temp1 = (2.0 * beta) -1;
omegabar = sqrt(np * temp1);
c3 = (beta3 - 1.) ./ 3.0 ./ beta2 ./ omegabar;
c4 = (beta3 .* beta3 - 1.0) ./ 4. ./ np ./beta2 ./beta2;

c6 = c3 .* c3 /2.0;
V = (vt - np * (1.0 + snrbar)) ./ omegabar;
Vsqr = V .*V;
val1 = exp(-Vsqr / 2.0) / sqrt( 2.0 * pi);
val2 = (c3 .* (V.^2 -1.0)) + (c4 .* V .* (3.0 - V.^2)) - (c6 .* V .* ((V.^4) - (10. * V.^2) + 15.0));
q = 0.5 * erfc (V/sqrt(2.0));
pd4 = q - val1 .* val2;

else
snr = 1.0;
gamma0 = gammainc(vt./beta,np);
a1 = (vt ./ beta).^np ./ (factorial(np) .* exp(vt./beta));
sum = gamma0;
for i = 1:1:np
temp1 = 1;
if (i == 1)
ai = a1;
else
ai = (vt ./ beta) .* a1 ./ (np + i -1);
end
a1 = ai;
gammai = gamma0 - ai;
gamma0 = gammai;
a1 = ai;
for ii = 1:1:i
temp1 = temp1 * (np + 1 - ii);
end
term = (snrbar /2.0).^i .* gammai .* temp1 ./ factorial(i);
sum = sum + term;
end
pd4 = 1.0 - sum ./ beta.^np;
end
pd4 = max(pd4,0.);

%%
%find plane aspect angle

array_size=length(pdtgt);

%%%%%%%%%%%
%%%%%%%%%%%

reqangles=120

%%%%%%%%%%%
%%%%%%%%%%%
%added 1 airplane
reqangle=ceil(reqangles*1793/180);
avgangle=ceil((reqangles*6/180));
pd4(avgangle)=0.05*pd4(avgangle)+(0.95*pdtgt(reqangle));

for i=1:6
    for m=1:array_size
        error1(i,m)=abs(pd4(i)-pdtgt(m));
    end
      err1(i)=min(error1(i,:));
end

 angles=15:30:178;
err1
plot(angle,pdtgt,angle,error1,angles,pd4);
aspct=min(err1);
for i=1:6
if err1(i)==aspct
    p=i
end
end


% for k=1:array_size
%        temp=error1(p,k);
%        temp1=err1(p);
%      if temp==temp1
%          l= k;
%      end
% end

if p~=1
    for k=(p-1)*298:p*298
        temp=error1(p,k);
        temp1=err1(p);
        if abs(temp-temp1)<=0.0001
            l= k;
        end
        
    end
else
    for k=1:p*298
        temp=error1(p,k);
        temp1=err1(p);
        if abs(temp-temp1)<=0.0001
            l= k;
        end
        
    end
end

aspct_angle=l*180/array_size;
expct_aspct=p*30;
if p==1
    'plane between 0-30'
end
if p==2
    'plane between 30-60'
end
if p==3
    'plane between 60-90'
end
if p==4
    'plane between 90-120'
end
if p==5
    'plane between 120-150'
end
if p==6
    'plane between 150-180'
end
%%
% percentage sureity of airplane
percentage_sureity_of_approx=(pd2(p)-err1(p))*100/pd2(p)
pd1=pd1';
pd2=pd2';
pd3=pd3';
pd4=pd4';
p1=[pd1 pd2 pd3 pd4];
save('swerlingdata.dat', 'p1', '-ASCII');
type swerlingdata.dat;
save('planedata.dat', 'p','ans',  '-ASCII');
type planedata.dat;
%b=load('planedata.dat');
a=load('swerlingdata.dat');