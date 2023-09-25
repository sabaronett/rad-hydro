data=importdata('Rad.radhst',' ',4);
data=data.data;
dimdata=size(data);
dataline=dimdata(1,1);
angles=importdata('Rad_angles.txt');
dim=size(angles);
Nang=dim(1,1);
Ix=zeros(Nang,1);
Iy=zeros(Nang,1);
for i=1:Nang
   nx=angles(i,2);
   ny=angles(i,3);
   intensity=data(dataline,22+i);
   Ix(i,1)=intensity*nx;
   Iy(i,1)=intensity*ny;    
end

% now create analytic solution
vel=2.95619;
lorz=(1/(1-vel^2/100))^0.5;
I0=data(dataline,13)/(4*pi);
nz1=8.819171e-01;
nz2=3.333333e-01;
npoint=1000;
dtheta=2*pi/npoint;
nxsolution1=zeros(npoint,1);
nysolution1=zeros(npoint,1);
nxsolution2=zeros(npoint,1);
nysolution2=zeros(npoint,1);
norm1=(1-nz1^2)^0.5;
norm2=(1-nz2^2)^0.5;
for i=1:npoint
    theta=(i-0.5)*dtheta;
    tran1=1/(lorz*(1-vel*cos(theta)*norm1/10))^4;
    nxsolution1(i,1)=I0*tran1*cos(theta)*norm1;
    nysolution1(i,1)=I0*tran1*sin(theta)*norm1;
    tran2=1/(lorz*(1-vel*cos(theta)*norm2/10))^4;
    nxsolution2(i,1)=I0*tran2*cos(theta)*norm2;
    nysolution2(i,1)=I0*tran2*sin(theta)*norm2;    
end

output=zeros(npoint,4);
output(:,1)=nxsolution1;
output(:,2)=nysolution1;
output(:,3)=nxsolution2;
output(:,4)=nysolution2;

save('anglesolution.txt','output','-ascii');