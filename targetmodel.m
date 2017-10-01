%%
% %body
freq=9e9;
r=4;
h=50;
index = 0;
eps =0.00001;

lambda = 3.0e+8 / freq;
% Compute RCS from zero aspect to broadside
for theta = 0.0:.1:90-.5
index = index +1;
theta = theta * pi /180.;
rcsbody(index) = (lambda * r * sin(theta) /(8. * pi * (cos(theta))^2)) + eps;
end
% Compute RCS for broadside specular
theta = pi/2;
index = index +1;
rcsbody(index) = (2. * pi * h^2 * r / lambda )+ eps;
% Compute RCS from 90 to 180 degrees
for theta = 90+.5:.1:180.
index = index + 1;
theta = theta * pi / 180.;
rcsbody(index) = ( lambda * r * sin(theta) / (8. * pi * (cos(theta))^2)) + eps;
end
% Plot results
delta= 180/(index-1);
angle = 0:delta:180;
figure(2);
%plot(angle,10*log10(rcsbody),'k');
plot(angle,rcsbody,'k');
grid;
xlabel ('Aspect angle - degrees');
ylabel ('RCS - dBsm');
title ('body');
%%
% % %wings
a=35;
b=10;
freq=9e9;
count1=0;
count2=0;count3=0;
A = a * b / 2.;
lambda = 3.e+8 / 9.5e+8;
ka = 2. * pi / lambda;
kb = 2. *pi / lambda;
% Compute theta vector
%  = 0.01:.05:179;
theta_deg=-90:delta:90;
phi = (pi /180.) .* theta_deg;
theta=(pi /180.) .*75;

alpha = ka * cos(phi) .* sin(theta);
beta = kb * sin(phi) .* sin(theta);

sigmao1 = 0.25 *sin(phi).^2 .* ((2. * a / b) .* cos(phi) .* sin(beta) - sin(phi) .* sin(2. .* alpha)).^2;
fact1 = (alpha).^2 - (.5 .* beta).^2;
fact2 = (sin(alpha).^2 - sin(.5 .* beta).^2).^2;
fact1=abs(fact1(:)');
sigmao = gdivide(gadd(fact2 , sigmao1) , (fact1));
rcswings = ((4. * pi * A^2 / lambda^2) .* cos(theta).^2 .* sigmao) + eps;

rcsdb = 10. *log10(rcswings);
figure(1);
%plot(angle, 10. *log10(rcswings));
plot(angle,rcswings,'k');
grid;
xlabel ('Apsect angle - degrees');
ylabel ('RCS - dBsm');
title ('wings');
%%
%exhaust
r1=3;
r2=4;
h=1;
indicator=1;
freq=9e9;
index = 0;
eps = 0.000001;
%lambda = 3.0e+8 / freq;
% Comput half cone angle, alpha
alpha = atan(( r2 - r1)/h);
% Compute z1 and z2
z2 = r2 / tan(alpha);
z1 = r1 / tan(alpha);
delta = (z2^1.5 - z1^1.5)^2;
factor = (8. * pi * delta) / (9. * lambda);
large_small_end = indicator;
if (large_small_end == 1)
% Compute normal incidence, large end

normal_incidence = (180./pi) * ((pi /2) + alpha);
% Compute RCS from zero aspect to normal incidence
for theta = 0.001:.1:normal_incidence-.5
index = index +1;
theta = theta * pi /180.;
rcsexhaust(index) = (lambda * z1 * tan(alpha) *(tan(theta - alpha))^2) / ...
(8. * pi *sin(theta)) + eps;
end
%Compute broadside RCS
index = index +1;
rcs_normal = factor * sin(alpha) / ((cos(alpha))^4) + eps;
rcs(index) = rcs_normal;
% Compute RCS from broad side to 180 degrees
for theta = normal_incidence+.5:.1:180
index = index + 1;
theta = theta * pi / 180. ;
rcsexhaust(index) = (lambda * z2 * tan(alpha) *(tan(theta - alpha))^2) / ...
(8. * pi *sin(theta)) + eps;
end
else
% Compute normal incidence, small end
normal_incidence = (180./pi) * ((pi /2) - alpha);
% Compute RCS from zero aspect to normal incidence (large end)
for theta = 0.001:.1:normal_incidence-.5
index = index +1;
theta = theta * pi /180.;
rcsexhaust(index) = (lambda * z1 * tan(alpha) *(tan(theta + alpha))^2) / ...
(8. * pi *sin(theta)) + eps;
end
%Compute broadside RCS
index = index +1;
rcs_normal = factor * sin(alpha) / ((cos(alpha))^4) + eps;
rcsexhaust(index) = rcs_normal;
% Compute RCS from broad side to 180 degrees (small end of frustum)
for theta = normal_incidence+.5:.1:180
index = index + 1;
theta = theta * pi / 180. ;
rcsexhaust(index) = (lambda * z2 * tan(alpha) *(tan(theta + alpha))^2) / ...
(8. * pi *sin(theta)) + eps;
end
end
rcsexhaust(1793)=rcsexhaust(1792);
rcsexhaust(1792)=rcsexhaust(1791);
rcsexhaust(1793)=rcsexhaust(1791);
rcsexhaust(1)=rcsexhaust(2);
% Plot RCS versus aspect angle
delta = 180 /index;
angle = 0:delta:180;
figure(3);
%plot (angle,10*log10(rcsexhaust),'k');
plot(angle,rcsexhaust,'k');
grid;
xlabel ('Apsect angle - degrees');
ylabel ('RCS - dBsm');
title ('exhaust');
%%
%nose
r1=0.0001;
r2=4;
h=1;
indicator=2;
freq=9e9;
index = 0;
eps = 0.000001;
%lambda = 3.0e+8 / freq;
% Comput half cone angle, alpha
alpha = atan(( r2 - r1)/h);
% Compute z1 and z2
z2 = r2 / tan(alpha);
z1 = r1 / tan(alpha);
delta = (z2^1.5 - z1^1.5)^2;
factor = (8. * pi * delta) / (9. * lambda);
large_small_end = indicator;
if (large_small_end == 1)
% Compute normal incidence, large end

normal_incidence = (180./pi) * ((pi /2) + alpha);
% Compute RCS from zero aspect to normal incidence
for theta = 0.001:.1:normal_incidence-.5
index = index +1;
theta = theta * pi /180.;
rcsnose(index) = (lambda * z1 * tan(alpha) *(tan(theta - alpha))^2) / ...
(8. * pi *sin(theta)) + eps;
end
%Compute broadside RCS
index = index +1;
rcs_normal = factor * sin(alpha) / ((cos(alpha))^4) + eps;
rcs(index) = rcs_normal;
% Compute RCS from broad side to 180 degrees
for theta = normal_incidence+.5:.1:180
index = index + 1;
theta = theta * pi / 180. ;
rcsnose(index) = (lambda * z2 * tan(alpha) *(tan(theta - alpha))^2) / ...
(8. * pi *sin(theta)) + eps;
end
else
% Compute normal incidence, small end
normal_incidence = (180./pi) * ((pi /2) - alpha);
% Compute RCS from zero aspect to normal incidence (large end)
for theta = 0.001:.1:normal_incidence-.5
index = index +1;
theta = theta * pi /180.;
rcsnose(index) = (lambda * z1 * tan(alpha) *(tan(theta + alpha))^2) / ...
(8. * pi *sin(theta)) + eps;
end
%Compute broadside RCS
index = index +1;
rcs_normal = factor * sin(alpha) / ((cos(alpha))^4) + eps;
rcsexhaust(index) = rcs_normal;
% Compute RCS from broad side to 180 degrees (small end of frustum)
for theta = normal_incidence+.5:.1:180
index = index + 1;
theta = theta * pi / 180. ;
rcsnose(index) = (lambda * z2 * tan(alpha) *(tan(theta + alpha))^2) / ...
(8. * pi *sin(theta)) + eps;
end
end
% Plot RCS versus aspect angle
delta = 180 /index;
angle = 0:delta:180;
rcsnose(1793)=rcsnose(1792);
figure(4);
% plot (angle,10*log10(rcsnose),'k');
plot(angle,rcsnose,'k');
grid;
xlabel ('Apsect angle - degrees');
ylabel ('RCS - dBsm');
title ('nose');
%%
%tail wing
A = a * b / 2.;
lambda = 3.e+8 / 9.5e+8;
phi = pi / 2.;
ka = 2. * pi / lambda;
kb = 2. *pi / lambda;
% Compute theta vector
theta_deg = -90:delta:90;
theta = (pi /180.) .* theta_deg;
alpha = ka * cos(phi) .* sin(theta);
beta = kb * sin(phi) .* sin(theta);

rcstail = (4. * pi * A^2 / lambda^2) .* cos(theta).^2 .* (sin(beta ./ 2)).^4./ (beta./2).^4 + eps;
rcstail=max(rcstail,0.001);
figure(5);
%plot(angle,10*log10(rcstail),'k');
plot(angle,rcstail,'k');
grid;
title ('tail wing');
%%
%net rcs;

netrcs= rcsnose+rcstail+rcsexhaust+rcsbody+rcswings;
netrcs(1792)=netrcs(1791);
netrcs(1793)=netrcs(1791);
netrcs(1)=netrcs(2);
figure(6);
normrcs=netrcs/max(netrcs);
plot(angle,10*log10(netrcs),'k');
%plot(angle,netrcs,'k');
grid;
 xlabel ('Apsect angle - degrees');
 ylabel ('RCS - dBsm');
 title ('plane rcs at aspect angles');
 %%
 %swerling3
nfa=1e9;
np=100;
snrbar=abs(normrcs);
snrbar = 10.0.^(snrbar/10.);
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
pdtgt=pd3;
figure(7);
plot(angle,pd3,angle,mean(pd3));
grid;
xlabel ('Aspect angle-degrees');
ylabel ('Pd');
