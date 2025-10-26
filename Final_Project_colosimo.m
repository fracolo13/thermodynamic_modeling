%% Thermomechanical code with variable conductivity and thermal Advection for Slab break off scenario

% 1) Clear memory and figures
clear all;
clf;

%2) Define numerical model
xsize = 500000; %horizontal model size, m
ysize =400000 ; %vertical model size, m

Nx = 101; % horizontal resolution of Eulerian nodes
Ny = 81; % vertical resolution of Eulerian nodes

Nxm = (Nx-1)*5; % horizontal resolution of Lagrangian markers (5 * 5 markers per cell)
Nym = (Ny-1)*5; % vertical resolution of Lagrangian markers
Nm=Nxm*Nym; % total number of Lagrangian markers

dx = xsize/(Nx-1);% horizontal grid step, m
dy = ysize/(Ny-1); %vertical grid step, m

x = 0:dx:(xsize+dx); %horizontal coordinates of points, m
y = 0:dy:(ysize+dy); % vertical coordinates of points, m

x_vx=0:dx:(xsize+dx); %horizontal coordinates of points for vx, m
y_vx=-dy/2:dy:(ysize+dy/2); %vertical coordinates of points for vx, m

x_vy=-dx/2:dx:(xsize+dx/2); %horizontal coordinates of points for vy, m
y_vy=0:dy:(ysize+dy); %vertical coordinates ofpoints for vy, m

x_p=(-dx/2):dx:(xsize+dx/2); %horizontal coordinates of points for p and T, m 
y_p=(-dy/2):dy:(ysize+dy/2);%vertical coordinates of points for p and T, m

dxm = xsize/Nxm; % horizontal grid step of moving markers
dym = ysize/Nym; %vertical grid step of moving markers (1 step more necessary than for eularian nodes)

%Arrays for Eulerian grid nodes
ETAB = zeros(Ny+1,Nx+1) ; % Matrix for etaB
ETAP = zeros(Ny+1,Nx+1) ; % Matrix for etaP
RHOvy = zeros(Ny+1,Nx+1); %Matrix for density on vy points
Kx=zeros(Ny+1,Nx+1); %Matrix for qx thermal conductivity
Ky=zeros(Ny+1,Nx+1); %Matrix for qy thermal conductivity
T0 = zeros(Ny+1,Nx+1);%Matrix for initial temperature
Tdt=zeros(Ny+1,Nx+1);%Matrix for temperature after DT
Hrad=zeros(Ny+1,Nx+1);%Matrix for radioactive heating
Ha=zeros(Ny+1,Nx+1);%Matrix for adiabatic heating
Hs=zeros(Ny+1,Nx+1);%Matrix for shear heating
alpha = zeros(Ny+1,Nx+1);%Matrix for thermal exp. coeff
beta = zeros(Ny+1,Nx+1);
RHOp = zeros(Ny+1,Nx+1);%Matrix for density on pressure points
Sigmay=zeros(Ny+1,Nx+1);%Matrix for yield stress 

SXX= zeros(Ny+1,Nx+1);% Matrices for deviatoric stresses
SYY= zeros(Ny+1,Nx+1);
SXY= zeros(Ny+1,Nx+1);
EXX= zeros(Ny+1,Nx+1);% Matrices for deviatoris deformation rates
EYY= zeros(Ny+1,Nx+1);
EXY= zeros(Ny+1,Nx+1);

%Arrays for Lagrangian grid nodes
xm = zeros(Nm,1); % horizontal coordinates for moving markers
ym = zeros(Nm,1); % vertical coordinates for moving markers
vxm = zeros(Nm,1); % horizontal velocity of markers
vym = zeros(Nm,1); % vertical velocity of markers
RHOm = zeros(Nm,1);%Density value for markers, kg/m^3
CP=zeros(Nm,1); %heat capacity
ETAm = zeros(Nm,1); %Viscosity for markers, Pa*s
Tm = zeros(Nm,1); % Temperature for markers
RHOCPm = zeros(Nm,1); % isobaric volumetric heatcapacity on markers
Km = zeros(Nm,1); %conductivity on markers
alpham = zeros(Nm,1); %thermal expansion coefficient (1/K)
betam=zeros(Nm,1);
Hradm=zeros(Nm,1); %radioactive Heating W/m^3
Sigmay=zeros(Nm,1);%yield stress
Pm=zeros(Nm,1); %pressure
EXXm=zeros(Nm,1);
EYYm=zeros(Nm,1);
EXYm=zeros(Nm,1);
EII=zeros(Nm,1); %deviatoric strain rate
SII=zeros(Nm,1); %deviatoric stress
%Types for distinction dry and wet olivine and slab, plate, neck, mantle for viscosity calculation
%0 = air
%1 = dry olivine in lithosphere and slab 
%2 = dry olivine in necking area
%3 = wet olivine in mantle
type=zeros(Nm,1); 


%Define BC-type: 
%Free slip: 
%vx=0, dvy/dx=0
%vy=0, dvx/dy=0
%BCtype=-1
%
%No slip: 
%vx=0
%vy=0
%BCtype=1
%
BCtype=-1;%insert +-1 
%
%
%Timestep calculations
%Possible between 0.1 and 1; max displacement
dcoef = 0.5; %0.1 bc of drunken sailor instability--> needs more time iterations
Tchangemax=50; %max temperature change per timestep

%% Constants
gy = 10; % gravitational acceleration, m/s^2
%Temperature dependent thermal conductivity
K0=0.73; % Initial thermal conductivity
a=1293; %Coefficient for T dependent thermal conductivity
b=77; %Coefficient for T dependent thermal conductivity

%non-linear viscosity

%boundaries for viscosity
ETAmin = 1e+18; 
ETAmax=1e+24; 

%dry olivine
ADd=2.5e-17;
nd=3.5;
Vad=8e-6;
Ead=532000;
RHO0d=3400; 
%wet olivine
ADw = 2e-21;
nw=4;
Vaw=4e-6;
Eaw=471000;
RHO0w=3300;



%isobaric heat capacties J/(kgK)
CP_rock= 1000;
CP_air=3.3e+6 ;


%% Variable density, viscosity, thermal conductivity and heat capacity
m = 1; %starting point of iteration

T_BCi = 1573; %Temperature for horizontal boundary condition, lower boundary
T_BCair = 273; %Temperature for horizontal bc, upper boundary
for jm=1:1:Nxm
    for im = 1:1:Nym        
        xm(m) = dxm/2+(jm-1)*dxm ;%+ (rand-0.5)*dxm; %remove rand term for reference values
        ym(m) = dym/2+(im-1)*dym ;%+ (rand-0.5)*dym; 
        %define RHOm, ETAm, Km and T0 on lagrangian markers
        if(ym(m)<ysize*0.125) %Sticky air layer --> 0-50 km depth (of 400)
            type(m)=0; %no type bc no changing viscosity
            RHOm(m) = 1; %Density
            ETAm(m) = 1e+18; %viscosity           
            RHOCPm(m)=RHOm(m)*CP_air; %volumetric heat capacity           
            Tm(m) = 273; %Temperature
            Hradm(m) = 0; %Radioactive heating 
            alpham(m) = 0; %Thermal expansion coeff
            betam(m)=0; %Coefficient for temperature dependent density
            Km(m)=3000;
        end
        if(ym(m)>=ysize*0.125 && ym(m)<=ysize*0.25) %Plate --> 50-100 km depth
            type(m)=1; %type dry olivin and lithosphere
            RHOm(m) = 3400; %Initial Density
            ETAm(m) = 1e+23; %viscosity           
            RHOCPm(m)=RHOm(m)*CP_rock; %volumetric heat capacity           
            Tm(m) = T_BCair + (ym(m)-ysize*0.125)*((T_BCi-T_BCair)/(ysize*0.25-ysize*0.125));  %Temperature gradient in plate
            Hradm(m) = 2e-8; %Radioactive heating 
            alpham(m) = 3e-5; %Thermal expansion coeff
            betam(m)=1e-11; %Coefficient for temperature dependent density
            Sigmay(m)=1e+8; %yield stress
            
            %add tilted Temperature gradient where slab goes down
            if(ym(m)>(-5*xm(m)+1200000) && ym(m)>(xm(m)-180000))
                Tm(m) =  T_BCair + (ym(m)+180000-xm(m))*((T_BCi-T_BCair)/(60000)) ;
            end
            
            Km(m)=K0+a/(Tm(m)+b); %T dependent thermal conductivity
        end
        if(ym(m)>ysize*0.25) %Mantle --> deeper than 100 km 
            type(m)=3; %type wet olivin and mantle
            RHOm(m) = 3250; %Initial Density
            ETAm(m) = 1e+20; %viscosity           
            RHOCPm(m)=RHOm(m)*CP_rock; %volumetric heat capacity           
            Tm(m) = 1573; %Temperature 
            Hradm(m) = 3e-8; %Radioactive heating 
            alpham(m) = 3e-5; %Thermal expansion coeff
            betam(m)=1e-11; %Coefficient for temperature dependent density
            Sigmay(m)=5e+7; %yield stress 
            Km(m)=K0+a/(Tm(m)+b); %T dependent thermal conductivity
        end
        if(ym(m)>ysize*0.25 && ym(m)<ysize*0.625 && xm(m)<(ym(m)+180000) && xm(m)>(ym(m)+120000)) %whole slab
            type(m)=1; %type dry olivin and lithosphere
            RHOm(m) = 3400; %Initial Density
            ETAm(m) = 1e+23; %viscosity           
            RHOCPm(m)=RHOm(m)*CP_rock; %volumetric heat capacity           
            Tm(m) =  T_BCair + (ym(m)+180000-xm(m))*((T_BCi-T_BCair)/(60000)) ; %Temperature gradient in plate
            Hradm(m) = 2e-8; %Radioactive heating 
            alpham(m) = 3e-5; %Thermal expansion coeff
            betam(m)=1e-11; %Coefficient for temperature dependent density
            Sigmay(m)=1e+8; %yield stress
            Km(m)=K0+a/(Tm(m)+b); %T dependent thermal conductivity
        end
        if(ym(m)>ysize*0.25 && ym(m)<ysize*0.375 && xm(m)<(ym(m)+180000) && xm(m)>(ym(m)+120000)) %Necking area of slab
            type(m)=2; %type dry olivin and necking area
            RHOm(m) = 3400; %Initial Density
            ETAm(m) = 1e+23; %viscosity
            RHOCPm(m)=RHOm(m)*CP_rock; %volumetric heat capacity
            Tm(m) =  T_BCair + (ym(m)+180000-xm(m))*((T_BCi-T_BCair)/(60000)) ;%Temperature gradient in plate
            Hradm(m) = 2e-8; %Radioactive heating
            alpham(m) = 3e-5; %Thermal expansion coeff
            betam(m)=1e-11; %Coefficient for temperature dependent density
            Sigmay(m)=2e+7; %yield stress
            Km(m)=K0+a/(Tm(m)+b); %T dependent thermal conductivity
        end
        
        m=m+1;
    end
end

% Arrays for subgrid temperature diffusion

%Eulerian nodes:
TREij= zeros(Ny+1,Nx+1); %non-relaxation degree
DTsubij=zeros(Ny+1,Nx+1); %DT subgrid
DTremij=zeros(Ny+1,Nx+1); %DT remain

%Lagrangian markers
TREm = zeros(Nm,1);
T0m=zeros(Nm,1);
DTm0=zeros(Nm,1);
DTm=zeros(Nm,1);
DTsubm=zeros(Nm,1);
DTremm=zeros(Nm,1);


%%

%3)Define global matrixes L() and R()

uk=3;% number of unknowns in mechanical equation
ukt=1;% number of unknowns in thermal equation

N=(Nx+1)*(Ny+1)*uk; %Number of unknowns = number of equations --> + 1 ghost row and 1 ghost column, for 3 unknowns to solve for
L=sparse(N,N); %left hand side values
R=zeros(N,1); %right hand side values

NT=(Nx+1)*(Ny+1)*ukt; %Number of unknowns = number of equations --> for 1 unknown to solve for
LT=sparse(NT,NT); %left hand side values
RT=zeros(NT,1); %right hand side values
pscale=2*min(ETAm)/(dx+dy);
%%
%4) Timeloop
ntimestep = 50; 
dt=1e+11; %initial time step for drunken sailor instability
timesum=0; %initial timesum
for t = 1:1:ntimestep
    %
    %%4a) Interpolation of RHO, ETA, Km, RHOCPm and T0 from Lagrangian Markers to Eulerian nodes
    % Starting sums for weighted interpolations
    RHOvywtsum=zeros(Ny+1, Nx+1);% sums for RHO on vy points; numerator
    wtsumvy = zeros(Ny+1, Nx+1);% denominator
    
    ETABwtsum =zeros(Ny+1, Nx+1); %sums for ETA on basic nodal points
    wtsumEB = zeros(Ny+1, Nx+1);
    
    ETAPwtsum = zeros(Ny+1, Nx+1); %sums for ETA on P points
    wtsumEP = zeros(Ny+1, Nx+1);
    
    Kxwtsum=zeros(Ny+1, Nx+1);% sums for Kx on vx points; numerator
    wtsumkx = zeros(Ny+1, Nx+1);% denominator
    
    Kywtsum =zeros(Ny+1, Nx+1); %sums for Ky on vy points
    wtsumky = zeros(Ny+1, Nx+1);
    
    RHOCPwtsum = zeros(Ny+1, Nx+1); %sums for RHOCP on P/T points
    wtsumrhocp = zeros(Ny+1, Nx+1);
    
    T0wtsum = zeros(Ny+1,Nx+1); %sums for T0 on P/T points
    wtsumt = zeros(Ny+1,Nx+1);
    
    Hradwtsum=zeros(Ny+1, Nx+1); %sums for radioactive heating on P/T points
    wtsumhrad=zeros(Ny+1, Nx+1);
    
    alphawtsum=zeros(Ny+1, Nx+1); %thermal expansion coef. on P/T points
    wtsumalpha=zeros(Ny+1, Nx+1);
    
    rhopwtsum=zeros(Ny+1, Nx+1); %density on P/T points
    wtsumrhop=zeros(Ny+1, Nx+1);
 
    %Going through all Lagrangian markers
    for m = 1:1:Nm
        %% RHOvy interpolation
        %
        % Indexing:
        j = fix((xm(m)-x_vy(1))/dx)+1;
        i = fix((ym(m)-y_vy(1))/dy)+1;
   
        % Distances between Eulerian nodes and Lagrangian marker for vy
        % points
        dxmj= abs(xm(m)-x_vy(j)); 
        dymi = abs(ym(m)-y_vy(i));
        % x = nodal grid point, o = lagrangian marker
        %(i,j)             (i,j+1)
        %x------------------x
        %|         |        |
        %|         |dymi    |
        %|         o--------|
        %|            dxmi  |
        %|                  |
        %x------------------x
        %(i+1,j)           (i+1,j+1)
        %Definition of weights (for four nodal points)
        wtmij = (1-dxmj/dx)*(1-dymi/dy);%indexes i and j not necessary bc all markers are inside grid
        wtmi1j =(1-dxmj/dx)*(dymi/dy);
        wtmij1 = (dxmj/dx)*(1-dymi/dy);
        wtmi1j1 = (dxmj/dx)*(dymi/dy);
        
        %Define sums
        RHOvywtsum(i,j) = RHOvywtsum(i,j) + wtmij*RHOm(m) ;
        RHOvywtsum(i+1,j) = RHOvywtsum(i+1,j) + wtmi1j*RHOm(m) ;
        RHOvywtsum(i,j+1) = RHOvywtsum(i,j+1) + wtmij1*RHOm(m) ;
        RHOvywtsum(i+1,j+1) = RHOvywtsum(i+1,j+1) + wtmi1j1*RHOm(m) ;
        
        wtsumvy(i,j) = wtsumvy(i,j) + wtmij ;
        wtsumvy(i+1,j) = wtsumvy(i+1,j) + wtmi1j ;
        wtsumvy(i,j+1) = wtsumvy(i,j+1) + wtmij1 ;
        wtsumvy(i+1,j+1) = wtsumvy(i+1,j+1) + wtmi1j1 ;
        
        %% ETAB interpolation
        % Indexing:
        j = fix((xm(m)-x(1))/dx)+1;
        i = fix((ym(m)-y(1))/dy)+1;
        % 
      
        
        % Distances between Eulerian nodes and Lagrangian marker for basic
        % nodal points
        
        dxmj= abs(xm(m)-x(j)); 
        dymi = abs(ym(m)-y(i));
        %
        %Definition of weights (for four nodal points)
        wtmij = (1-dxmj/dx)*(1-dymi/dy);
        wtmi1j =(1-dxmj/dx)*(dymi/dy);
        wtmij1 = (dxmj/dx)*(1-dymi/dy);
        wtmi1j1 = (dxmj/dx)*(dymi/dy);
        
        %Define sums
        ETABwtsum(i,j) = ETABwtsum(i,j) + wtmij*ETAm(m) ;
        ETABwtsum(i+1,j) = ETABwtsum(i+1,j) + wtmi1j*ETAm(m) ;
        ETABwtsum(i,j+1) = ETABwtsum(i,j+1) + wtmij1*ETAm(m) ;
        ETABwtsum(i+1,j+1) = ETABwtsum(i+1,j+1) + wtmi1j1*ETAm(m) ;
        
        wtsumEB(i,j) = wtsumEB(i,j) + wtmij ;
        wtsumEB(i+1,j) = wtsumEB(i+1,j) + wtmi1j ;
        wtsumEB(i,j+1) = wtsumEB(i,j+1) + wtmij1 ;
        wtsumEB(i+1,j+1) = wtsumEB(i+1,j+1) + wtmi1j1 ;
        
        %% ETAP interpolation
        % Indexing:
        j = fix((xm(m)-x_p(1))/dx)+1;
        i = fix((ym(m)-y_p(1))/dy)+1;
      
        % Distances between Eulerian nodes and Lagrangian marker for
        % pressure
        % nodal points
        
        dxmj= abs(xm(m)-x_p(j)); 
        dymi = abs(ym(m)-y_p(i));
        %
        %Definition of weights (for four nodal points)
        wtmij = (1-dxmj/dx)*(1-dymi/dy);
        wtmi1j =(1-dxmj/dx)*(dymi/dy);
        wtmij1 = (dxmj/dx)*(1-dymi/dy);
        wtmi1j1 = (dxmj/dx)*(dymi/dy);
        
        %Define sums
        ETAPwtsum(i,j) = ETAPwtsum(i,j) + wtmij*ETAm(m) ;
        ETAPwtsum(i+1,j) = ETAPwtsum(i+1,j) + wtmi1j*ETAm(m) ;
        ETAPwtsum(i,j+1) = ETAPwtsum(i,j+1) + wtmij1*ETAm(m) ;
        ETAPwtsum(i+1,j+1) = ETAPwtsum(i+1,j+1) + wtmi1j1*ETAm(m) ;
        
        wtsumEP(i,j) = wtsumEP(i,j) + wtmij ;
        wtsumEP(i+1,j) = wtsumEP(i+1,j) + wtmi1j ;
        wtsumEP(i,j+1) = wtsumEP(i,j+1) + wtmij1 ;
        wtsumEP(i+1,j+1) = wtsumEP(i+1,j+1) + wtmi1j1 ;
        
        %% Kmx interpolation
        j = fix((xm(m)-x_vx(1))/dx)+1;
        i = fix((ym(m)-y_vx(1))/dy)+1;
        
        % Distances between Eulerian nodes and Lagrangian marker for vx
        % points
        dxmj= abs(xm(m)-x_vx(j));
        dymi = abs(ym(m)-y_vx(i));
        %x_vx = nodal points; o = lagrangian marker
        %(i,j)             (i,j+1)
        %x------------------x
        %|         |        |
        %|         |dymi    |
        %|         o--------|
        %|            dxmi  |
        %|                  |
        %x------------------x
        %(i+1,j)           (i+1,j+1)
        %Definition of weights (for four nodal points)
        wtmij = (1-dxmj/dx)*(1-dymi/dy);%indexes i and j not necessary bc all markers are inside grid
        wtmi1j =(1-dxmj/dx)*(dymi/dy);
        wtmij1 = (dxmj/dx)*(1-dymi/dy);
        wtmi1j1 = (dxmj/dx)*(dymi/dy);
        %Define sums
        Kxwtsum(i,j) = Kxwtsum(i,j) + wtmij*Km(m) ;
        Kxwtsum(i+1,j) = Kxwtsum(i+1,j) + wtmi1j*Km(m) ;
        Kxwtsum(i,j+1) = Kxwtsum(i,j+1) + wtmij1*Km(m) ;
        Kxwtsum(i+1,j+1) = Kxwtsum(i+1,j+1) + wtmi1j1*Km(m) ;
        
        wtsumkx(i,j) = wtsumkx(i,j) + wtmij ;
        wtsumkx(i+1,j) = wtsumkx(i+1,j) + wtmi1j ;
        wtsumkx(i,j+1) = wtsumkx(i,j+1) + wtmij1 ;
        wtsumkx(i+1,j+1) = wtsumkx(i+1,j+1) + wtmi1j1 ;
        
        %% Kmy interpolation
        j = fix((xm(m)-x_vy(1))/dx)+1;
        i = fix((ym(m)-y_vy(1))/dy)+1;
        
        % Distances between Eulerian nodes and Lagrangian marker for vx
        % points
        dxmj= abs(xm(m)-x_vy(j));
        dymi = abs(ym(m)-y_vy(i));
        %x_vy = nodal points; o = lagrangian marker
        %(i,j)             (i,j+1)
        %x------------------x
        %|         |        |
        %|         |dymi    |
        %|         o--------|
        %|            dxmi  |
        %|                  |
        %x------------------x
        %(i+1,j)           (i+1,j+1)
        %Definition of weights (for four nodal points)
        wtmij = (1-dxmj/dx)*(1-dymi/dy);%indexes i and j not necessary bc all markers are inside grid
        wtmi1j =(1-dxmj/dx)*(dymi/dy);
        wtmij1 = (dxmj/dx)*(1-dymi/dy);
        wtmi1j1 = (dxmj/dx)*(dymi/dy);
        %Define sums
        Kywtsum(i,j) = Kywtsum(i,j) + wtmij*Km(m) ;
        Kywtsum(i+1,j) = Kywtsum(i+1,j) + wtmi1j*Km(m) ;
        Kywtsum(i,j+1) = Kywtsum(i,j+1) + wtmij1*Km(m) ;
        Kywtsum(i+1,j+1) = Kywtsum(i+1,j+1) + wtmi1j1*Km(m) ;
        
        wtsumky(i,j) = wtsumky(i,j) + wtmij ;
        wtsumky(i+1,j) = wtsumky(i+1,j) + wtmi1j ;
        wtsumky(i,j+1) = wtsumky(i,j+1) + wtmij1 ;
        wtsumky(i+1,j+1) = wtsumky(i+1,j+1) + wtmi1j1 ;
        
        
        %% RHOCP interpolation
        %
        % Indexing:
        j = fix((xm(m)-x_p(1))/dx)+1;
        i = fix((ym(m)-y_p(1))/dy)+1;
        
        % Distances between Eulerian nodes and Lagrangian marker for vy
        % points
        dxmj= abs(xm(m)-x_p(j));
        dymi = abs(ym(m)-y_p(i));
        % x = nodal grid point, o = lagrangian marker
        %(i,j)             (i,j+1)
        %x------------------x
        %|         |        |
        %|         |dymi    |
        %|         o--------|
        %|            dxmi  |
        %|                  |
        %x------------------x
        %(i+1,j)           (i+1,j+1)
        %Definition of weights (for four nodal points)
        wtmij = (1-dxmj/dx)*(1-dymi/dy);%indexes i and j not necessary bc all markers are inside grid
        wtmi1j =(1-dxmj/dx)*(dymi/dy);
        wtmij1 = (dxmj/dx)*(1-dymi/dy);
        wtmi1j1 = (dxmj/dx)*(dymi/dy);
        
        %Define sums
        RHOCPwtsum(i,j) = RHOCPwtsum(i,j) + wtmij*RHOCPm(m) ;
        RHOCPwtsum(i+1,j) = RHOCPwtsum(i+1,j) + wtmi1j*RHOCPm(m) ;
        RHOCPwtsum(i,j+1) = RHOCPwtsum(i,j+1) + wtmij1*RHOCPm(m) ;
        RHOCPwtsum(i+1,j+1) = RHOCPwtsum(i+1,j+1) + wtmi1j1*RHOCPm(m) ;
        
        wtsumrhocp(i,j) = wtsumrhocp(i,j) + wtmij ;
        wtsumrhocp(i+1,j) = wtsumrhocp(i+1,j) + wtmi1j ;
        wtsumrhocp(i,j+1) = wtsumrhocp(i,j+1) + wtmij1 ;
        wtsumrhocp(i+1,j+1) = wtsumrhocp(i+1,j+1) + wtmi1j1 ;
        
        
        %% T0 interpolation
        %
        % Indexing:
        j = fix((xm(m)-x_p(1))/dx)+1;
        i = fix((ym(m)-y_p(1))/dy)+1;
        
        % Distances between Eulerian nodes and Lagrangian marker for vy
        % points
        dxmj= abs(xm(m)-x_p(j));
        dymi = abs(ym(m)-y_p(i));
        % x = nodal grid point, o = lagrangian marker
        %(i,j)             (i,j+1)
        %x------------------x
        %|         |        |
        %|         |dymi    |
        %|         o--------|
        %|            dxmi  |
        %|                  |
        %x------------------x
        %(i+1,j)           (i+1,j+1)
        %Definition of weights (for four nodal points)
        wtmij = (1-dxmj/dx)*(1-dymi/dy);%indexes i and j not necessary bc all markers are inside grid
        wtmi1j =(1-dxmj/dx)*(dymi/dy);
        wtmij1 = (dxmj/dx)*(1-dymi/dy);
        wtmi1j1 = (dxmj/dx)*(dymi/dy);
        
        %Define sums
        T0wtsum(i,j) = T0wtsum(i,j) + wtmij*Tm(m)*RHOCPm(m) ;%Heat conservative interpolation
        T0wtsum(i+1,j) = T0wtsum(i+1,j) + wtmi1j*Tm(m)*RHOCPm(m) ;
        T0wtsum(i,j+1) = T0wtsum(i,j+1) + wtmij1*Tm(m)*RHOCPm(m) ;
        T0wtsum(i+1,j+1) = T0wtsum(i+1,j+1) + wtmi1j1*Tm(m)*RHOCPm(m) ;
        
        wtsumt(i,j) = wtsumt(i,j) + wtmij*RHOCPm(m) ;
        wtsumt(i+1,j) = wtsumt(i+1,j) + wtmi1j*RHOCPm(m) ;
        wtsumt(i,j+1) = wtsumt(i,j+1) + wtmij1*RHOCPm(m);
        wtsumt(i+1,j+1) = wtsumt(i+1,j+1) + wtmi1j1*RHOCPm(m);
            
        %% Hrad interpolation
        % Indexing:
        j = fix((xm(m)-x_p(1))/dx)+1;
        i = fix((ym(m)-y_p(1))/dy)+1;
      
        % Distances between Eulerian nodes and Lagrangian marker for
        % pressure
        % nodal points
        
        dxmj= abs(xm(m)-x_p(j)); 
        dymi = abs(ym(m)-y_p(i));
        %
        %Definition of weights (for four nodal points)
        wtmij = (1-dxmj/dx)*(1-dymi/dy);
        wtmi1j =(1-dxmj/dx)*(dymi/dy);
        wtmij1 = (dxmj/dx)*(1-dymi/dy);
        wtmi1j1 = (dxmj/dx)*(dymi/dy);
        
        %Define sums
        Hradwtsum(i,j) = Hradwtsum(i,j) + wtmij*Hradm(m) ;
        Hradwtsum(i+1,j) = Hradwtsum(i+1,j) + wtmi1j*Hradm(m) ;
        Hradwtsum(i,j+1) = Hradwtsum(i,j+1) + wtmij1*Hradm(m) ;
        Hradwtsum(i+1,j+1) = Hradwtsum(i+1,j+1) + wtmi1j1*Hradm(m) ;
        
        wtsumhrad(i,j) = wtsumhrad(i,j) + wtmij ;
        wtsumhrad(i+1,j) = wtsumhrad(i+1,j) + wtmi1j ;
        wtsumhrad(i,j+1) = wtsumhrad(i,j+1) + wtmij1 ;
        wtsumhrad(i+1,j+1) = wtsumhrad(i+1,j+1) + wtmi1j1 ; 
        
        
          %% alpha interpolation
        % Indexing:
        j = fix((xm(m)-x_p(1))/dx)+1;
        i = fix((ym(m)-y_p(1))/dy)+1;
      
        % Distances between Eulerian nodes and Lagrangian marker for
        % pressure
        % nodal points
        
        dxmj= abs(xm(m)-x_p(j)); 
        dymi = abs(ym(m)-y_p(i));
        %
        %Definition of weights (for four nodal points)
        wtmij = (1-dxmj/dx)*(1-dymi/dy);
        wtmi1j =(1-dxmj/dx)*(dymi/dy);
        wtmij1 = (dxmj/dx)*(1-dymi/dy);
        wtmi1j1 = (dxmj/dx)*(dymi/dy);
        
        %Define sums
        alphawtsum(i,j) = alphawtsum(i,j) + wtmij*alpham(m) ;
        alphawtsum(i+1,j) = alphawtsum(i+1,j) + wtmi1j*alpham(m) ;
        alphawtsum(i,j+1) = alphawtsum(i,j+1) + wtmij1*alpham(m) ;
        alphawtsum(i+1,j+1) = alphawtsum(i+1,j+1) + wtmi1j1*alpham(m) ;
        
        wtsumalpha(i,j) = wtsumalpha(i,j) + wtmij ;
        wtsumalpha(i+1,j) = wtsumalpha(i+1,j) + wtmi1j ;
        wtsumalpha(i,j+1) = wtsumalpha(i,j+1) + wtmij1 ;
        wtsumalpha(i+1,j+1) = wtsumalpha(i+1,j+1) + wtmi1j1 ; 
        
           %% RHOp interpolation
        % Indexing:
        j = fix((xm(m)-x_p(1))/dx)+1;
        i = fix((ym(m)-y_p(1))/dy)+1;
      
        % Distances between Eulerian nodes and Lagrangian marker for
        % pressure
        % nodal points
        
        dxmj= abs(xm(m)-x_p(j)); 
        dymi = abs(ym(m)-y_p(i));
        %
        %Definition of weights (for four nodal points)
        wtmij = (1-dxmj/dx)*(1-dymi/dy);
        wtmi1j =(1-dxmj/dx)*(dymi/dy);
        wtmij1 = (dxmj/dx)*(1-dymi/dy);
        wtmi1j1 = (dxmj/dx)*(dymi/dy);
        
        %Define sums
        rhopwtsum(i,j) = rhopwtsum(i,j) + wtmij*RHOm(m) ;
        rhopwtsum(i+1,j) = rhopwtsum(i+1,j) + wtmi1j*RHOm(m) ;
        rhopwtsum(i,j+1) = rhopwtsum(i,j+1) + wtmij1*RHOm(m) ;
        rhopwtsum(i+1,j+1) = rhopwtsum(i+1,j+1) + wtmi1j1*RHOm(m) ;
        
        wtsumrhop(i,j) = wtsumrhop(i,j) + wtmij ;
        wtsumrhop(i+1,j) = wtsumrhop(i+1,j) + wtmi1j ;
        wtsumrhop(i,j+1) = wtsumrhop(i,j+1) + wtmij1 ;
        wtsumrhop(i+1,j+1) = wtsumrhop(i+1,j+1) + wtmi1j1 ; 
        
    end 
    %Calculate Eulerian RHO and ETAB/P, Kx, Ky, RHOCP, T0

    for j=1:1:Nx+1
        for i=1:1:Ny+1
            if(wtsumvy(i,j)>0)
                RHOvy(i,j) = RHOvywtsum(i,j)/wtsumvy(i,j);
            end
            if(wtsumEB(i,j)>0)
                ETAB(i,j) = ETABwtsum(i,j)/wtsumEB(i,j);
            end
            if(wtsumEP(i,j)>0)
                ETAP(i,j) = ETAPwtsum(i,j)/wtsumEP(i,j);
            end 
            if(wtsumkx(i,j)>0)
                Kx(i,j) = Kxwtsum(i,j)/wtsumkx(i,j);
            end
            if(wtsumky(i,j)>0)
                Ky(i,j) = Kywtsum(i,j)/wtsumky(i,j);
            end
            if(wtsumrhocp(i,j)>0)
                RHOCP(i,j) = RHOCPwtsum(i,j)/wtsumrhocp(i,j);
            end
            if(wtsumt(i,j)>0)
                T0(i,j) = T0wtsum(i,j)/wtsumt(i,j);
            end
            if(wtsumhrad(i,j)>0)
                Hrad(i,j)=Hradwtsum(i,j)/wtsumhrad(i,j);
            end
            if(wtsumalpha(i,j)>0)
                alpha(i,j)=alphawtsum(i,j)/wtsumalpha(i,j);
            end
            if(wtsumrhop(i,j)>0)
                RHOp(i,j)=rhopwtsum(i,j)/wtsumrhop(i,j);
            end
            
        end
    end
    
    %4b) Boundaries for T0   
    %horizontal boundary conditions: no slip --> (T1+T2)/2=T_BCi
    T0(1,:)=2*T_BCair - T0(2,:); %T1=2*T_BCi -T2   
    T0(Ny+1,:)=2*T_BCi - T0(Ny,:); %T2=2*T_BCi -T1    
    %vertical boundary conditions: free slip --> (T2-T1)/2=0   
    T0(:,1)=T0(:,2); %T1=T2   
    T0(:,Nx+1)=T0(:,Nx); %T2=T1
   
    %5) Introducing thermomechanical iteration
    nitermax=2;% at least 2
    %Probing increase of the timestep
    dt=dt*1.1; % the coefficient 1.1 implies gradual (by 10%max increase in the timestep size (eventhough we have deleted 1e30 --> so we dont have deadlock with small temp)
   
    for niter=1:1:nitermax %try to adjust timestep within iteration
    
        % 6) Composing of global matrixes L(), R()
        % going through all points of the grid (2D means two loops) 
        %First loop - horizontal index j
        for j=1:1:Nx+1 %+1 ghost column
            %Second loop - vertical index i
            for i=1:1:Ny+1 % +1 ghost row
                %
                %Define global indexes
                gvx=((j-1)*(Ny+1)+(i-1))*uk+1; %index for vx
                gvy=gvx+1;%index for vy
                gp=gvx+2; %index for p
                %
                %4a) vx through x-stokes
                %Decide which equation to solve: BC or Poisson equation
                if(j==Nx+1)%Ghost points-->0 (Fictious BC)
                    L(gvx,gvx)=1;
                    R(gvx,1)=0;
                elseif(j==1 || j==Nx) %real BC on j axis
                    %BC; vx(i,j)=0 --> 1*S(g)=0
                    L(gvx,gvx)=1; % left hand side for ; (g of eq left, and of unknown right)
                    R(gvx,1)= 0;% right hand side
                elseif(i==1) %real BC on i axis
                    %BC: (-+1)vx3+1vx4)=0
                    %left hand side
                    L(gvx,gvx+uk)=1;%vx3
                    L(gvx,gvx)=BCtype;%vx4
                    R(gvx,1)=0;%right hand side
                elseif(i==Ny+1)
                    L(gvx, gvx-uk)=BCtype;%vx3
                    L(gvx,gvx)=1;%vx4
                    R(gvx,1)=0;
                else
                    %
                    %                    S(gvx-3)
                    %                      vx_2
                    %          S(gvy-3)      |                S(gvy+(Ny+1)*3-3)
                    %           vy1      SIGMA'xy1            vy3
                    %                      ETAB(i-1,j)
                    %                        |
                    %                        |
                    %                        |
                    %S(gvx-(Ny+1)*3)      S(gvx)       S(gvx+(Ny+1)*3)
                    %vx_1--------P1--------vx_3------P2-------vx_5
                    %          S(gp)         |      S(gp+(Ny+1)*3
                    %           ETAP(i,j)    |       ETAP(i,j+1)
                    %                      ETAB(i,j)
                    %           vy2       SIGMA'xy2           vy4
                    %          S(gvy)        |               S(gvy+(Ny+1)*3)
                    %                     S(gvx+3)
                    %                      vx_4
                    %
                    %(SIGMA'xx2-SIGMA'xx1)/dx+(SIGMA'xy2-SIGMA'xy1/dy)-(Px2-Px1)/dx=0
                    %SIGMA'xx1=2ETAP(i,j)(vx3-vx1)/dx
                    %SIGMA'xx2=2ETAP(i,j+1)(vx5-vx3)/dx
                    %SIGMA'xy1=ETAB(i-1,j)((vx3-vx2)/dy+(vy3-vy1)/dx
                    %SIGMA'xy2=ETAB(i,j)((vx4-vx3)/dy+(vy4-vy2)/dx))
                    %Left hand side (11 unknowns)
                    
                    L(gvx,gvx-(Ny+1)*uk)=2*ETAP(i,j)/dx^2;%vx1
                    L(gvx,gvx-uk)=ETAB(i-1,j)/dy^2;%vx2
                    L(gvx,gvx)=-2*ETAP(i,j+1)/dx^2-2*ETAP(i,j)/dx^2-ETAB(i,j)/dy^2-ETAB(i-1,j)/dy^2;%vx3
                    L(gvx, gvx+uk)=ETAB(i,j)/dy^2;%vx4
                    L(gvx, gvx+(Ny+1)*uk)=2*ETAP(i,j+1)/dx^2;%vx5
                    
                    L(gvx, gvy-uk)=ETAB(i-1,j)/(dx*dy);%vy1
                    L(gvx, gvy)=-ETAB(i,j)/(dx*dy);%vy2
                    L(gvx, gvy+(Ny+1)*uk-uk)=-ETAB(i-1,j)/(dx*dy);%vy3
                    L(gvx, gvy+(Ny+1)*uk)=ETAB(i,j)/(dx*dy);%vy4
                    
                    L(gvx, gp)= pscale/dx;%P1
                    L(gvx,gp+(Ny+1)*uk)=-pscale/dx;%P2
                    
                    %Right hand side of x-stokes
                    R(gvx,1)=0;
                end
                
                %
                %vy through y-stokes
                if(i==Ny+1)%Ghost points --> 0 (Fictious BC)
                    L(gvy,gvy)=1;
                    R(gvy,1)=0;
                elseif(j==1)% real BC in j direction
                    %Left hand side
                    %BC: (+-1)vy3+1vy5=0
                    L(gvy,gvy)=BCtype;%vy3
                    L(gvy, gvy+(Ny+1)*uk) = 1;%vy5
                    R(gvy,1)= 0;% right hand side
                elseif(j==Nx+1)
                    L(gvy,gvy-(Ny+1)*uk)=BCtype;%vy3
                    L(gvy,gvy)=1;%vy5
                    R(gvy,1)=0;%right hand side
                elseif(i==1 || i==Ny)%real BC in i direction
                    %BC: 1vy=0
                    L(gvy,gvy)=1;%left hand side
                    R(gvy,1)=0;%right hand side
                else
                    %
                    %                            S(gvy-3)
                    %                              vy_2
                    %         S(gvx-(Ny+1)*3)        |               S(gvx)
                    %                             ETAP(i,j)
                    %           vx1              P1 --> S(gvy)          vx3
                    %                            SIGMA'YY1
                    %                                |
                    %         ETAB(i,j-1)      S(gvy)               ETAB(i,j)
                    %vy_1-----SIGMA'xy1------------vy_3-------------SIGMA'xy2-------------vy_5
                    %S(gvy-(Ny+1)*3)                 |                             S(gvy+(Ny+1)*3)
                    %           vx2            P2 --> S(gvy+3)           vx4
                    %           S(gvx-(Ny+1)*3+3)     SIGMA'YY2        S(gvx+3)
                    %                              ETAP(i+1,j)
                    %                                |
                    %                              S(gvy+3)
                    %                                vy_4
                    %
                    %(SIGMA'xy2-SIGMA'xy1)/dx+(SIGMA'yy2-SIGMA'yy1)/dy-(P2-P1)/dy=-RHOvy(i,j)*gy
                    %SIGMA'yy1 = 2ETAP(i,j)(vy3-vy2)/dy
                    %SIGMA'yy2 = 2ETAP(i+1,j)(vy4-vy3)/dy
                    %SIGMA'xy1=ETAB(i,j-1)((vy3-vy1)/dx+(vx2-vx1)/dy)
                    %SIGMA'xy2=ETAB(i,j)((vy5-vy3)7dx+(vx4-vx3)/dy)
                    %
                    %Corrections for drunken sailor instability:
                    %vy=vy3:
                    %corrvy3=-dt*gy*(RHOvy(i+1,j)-RHOvy(i-1,j))/(2*dy);
                    dRHOvydy=(RHOvy(i+1,j)-RHOvy(i-1,j))/(2*dy);%-dt*dRHOvydy*gy
                    %vx=(vx1+vx2+vx3+vx4)/4
                    %corrvx=-dt/4*gy*(RHOvy(i,j+1)-RHOvy(i,j-1))/(2*dx);
                    dRHOvydx=(RHOvy(i,j+1)-RHOvy(i,j-1))/(2*dx);%-1/4*dt*dRHOvydx*gy
                    %
                    %
                    %Left hand side (11 unknowns)
                    L(gvy,gvy-(Ny+1)*uk)=ETAB(i,j-1)/dx^2;%vy1
                    L(gvy,gvy-uk)=2*ETAP(i,j)/dy^2;%vy2
                    L(gvy,gvy)=-ETAB(i,j)/dx^2-ETAB(i,j-1)/dx^2-2*ETAP(i+1,j)/dy^2-2*ETAP(i,j)/dy^2-dt*dRHOvydy*gy;%vy3
                    L(gvy, gvy+uk)=2*ETAP(i+1,j)/dy^2;%vy4
                    L(gvy, gvy+(Ny+1)*uk)=ETAB(i,j)/dx^2;%vy5
                    
                    L(gvy,gvx-(Ny+1)*uk)=ETAB(i,j-1)/(dx*dy)-1/4*dt*dRHOvydx*gy; %vx1
                    L(gvy, gvx-(Ny+1)*uk+uk)=-ETAB(i,j-1)/(dx*dy)-1/4*dt*dRHOvydx*gy; %vx2
                    L(gvy, gvx)=-ETAB(i,j)/(dx*dy)-1/4*dt*dRHOvydx*gy;%vx3
                    L(gvy, gvx+uk)=ETAB(i,j)/(dx*dy)-1/4*dt*dRHOvydx*gy;%vx4
                    
                    L(gvy, gp)= pscale/dy;%P1
                    L(gvy,gp+uk)=-pscale/dy;%P2
                    
                    %Right hand side of y-stokes
                    R(gvy,1)=-RHOvy(i,j)*gy;
                end
                %
                %P with continuity equation
                if(j==1 || j==Nx+1 || i==1 || i==Ny+1) %ghost points (Fictious BC)
                    L(gp,gp)=pscale;
                    R(gp,1)=0;
                elseif(j==2 && i==2) %BC of real boundaries
                    % 1vp=10^9 Pa
                    %1vp=10^5 Pa for sticky air bc
                    L(gp,gp)=pscale;%left hand side
                    R(gp,1)=1e+5;%right hand side%lower bc of sticky air
                else
                    % Continuity equaiton: dvx/dx + dvy/dy=0
                    %                    S(gvy-3)
                    %                     vy1
                    %                      |
                    %                      |
                    %vx1-------------------P-------------------vx2
                    %S(gvx-(Ny+1)*3)      S(gp)               S(gvx)
                    %                      |
                    %                      |
                    %                     vy2
                    %                    S(gvy)
                    % (vx2-vx1)/dx +(vy2-vy1)/dy = 0
                    %
                    %Left hand side (4 unknowns)
                    L(gp, gvx-(Ny+1)*uk)=-pscale/dx;%vx1
                    L(gp, gvx) = pscale/dx;%vx2
                    L(gp, gvy-uk)=-pscale/dy;%vy1
                    L(gp, gvy)=pscale/dy;%vy2
                    %Right hand side
                    R(gp,1) = 0;
                end
                
            end
            
        end
        
        %7) Solve global matrix
        S=L\R; %contains a mixture of solutions
        
        % 8) Reload S (1D array to 2D)
        vx = zeros(Ny+1,Nx+1); %Create geometrical arrays
        vy=zeros(Ny+1,Nx+1);
        pr=zeros(Ny+1,Nx+1);
        
        % going through all points of the grid (2D means two loops) --> geometrical
        % view
        %First loop - horizontal index j
        for j=1:1:Nx+1
            %Second loop - vertical index i
            for i=1:1:Ny+1
                %Define global indexes
                gvx=((j-1)*(Ny+1)+(i-1))*uk+1;
                gvy=gvx+1;
                gp=gvx+2;
                %Reload S(gvx)--> vx(i,j)
                vx(i,j) = S(gvx);
                %Reload S(gvy) --> vy(i,j)
                vy(i,j)=S(gvy);
                %Reload S(gp)--> P(i,j)
                pr(i,j)=S(gp)*pscale;
            end
        end
        
        %9) adjust time step (not in final iteration)
        if(niter<nitermax)
            %vmax velocity:
            vxmax = max(max(abs(vx)));%2D --> 2times max
            vymax= max(max(abs(vy)));
            
            %starting dt
            
            if(dt > (dx*dcoef/vxmax))
                dt=dx*dcoef/vxmax;
               % dtmechanicalvx=dt
            end
            if(dt > (dy*dcoef/vymax))
                dt=dy*dcoef/vymax;
              %  dtmechanicalvy=dt
            end
          
        end
        
        
        %% -------------------------------------------------------------------------------------------
        %10) Thermal equation
        
        %% 10a) Adiabatic heating
        %averaged vy velocity
        %Ha=T*alpha*DP/DT=(T0+Tdt)/2*alpha*rhop*gy*vyaverage
        for j=2:1:Nx
            for i=2:1:Ny
                vyaverage=1/2*(vy(i,j)+vy(i-1,j));
                HaLH(i,j)=gy/2*alpha(i,j)*RHOp(i,j)*vyaverage;%left hand side of T3
                HaRH(i,j)=gy/2*alpha(i,j)*RHOp(i,j)*T0(i,j)*vyaverage;%right hand side
                %Ha(i,j)=HaLH(i,j)*Tdt(i,j)+HaRH(i,j);
            
            end
        end
        %% 10b) Shear heating
        %Hs=SXX*EXX + SYY*EYY + 2SXY*EXY
        %Sigma=2*ETA*Epsilon
        %SigmaXY&EpsilonXY
        for j=1:1:Nx %on basic nodes 
            for i=1:1:Ny
                EXY(i,j)=1/2*((vx(i+1,j)-vx(i,j))/dy+(vy(i,j+1)-vy(i,j))/dx);
                SXY(i,j)=2*ETAB(i,j)*EXY(i,j); 
            end
        end
        %EpsilonXX  &EpsilonYY: 
        for j=2:1:Nx  %on pressure nodes
            for i=2:1:Ny 
                EXX(i,j)=(vx(i,j)-vx(i,j-1))/dx;
              
                EYY(i,j)=(vy(i,j)-vy(i-1,j))/dy;
               
            end
        end
        
        %Symmetry BC for EXX, EYY, EXY and  isostatic pressure assumption for P
        %excluding corners in x direction:
        for j=2:1:Nx
            EXX(1,j)=EXX(2,j);
            EYY(1,j)=EYY(2,j);

            EXX(Ny+1,j)=EXX(Ny,j);
            EYY(Ny+1,j)=EYY(Ny,j);
            
            pr(1,j)=pr(2,j)-RHOvy(1,j)*gy*dy;
            pr(Ny+1,j)=pr(Ny,j)+RHOvy(Ny,j)*gy*dy;
        end
        %including corners in y direction 
        for i=1:1:Ny+1
            EXX(i,1)=EXX(i,2);
            EYY(i,1)=EYY(i,2);
            EXX(i,Nx+1)=EXX(i,Nx);
            EYY(i,Nx+1)=EYY(i,Nx);
            
            pr(i,1)=pr(i,2);
            pr(i,Nx+1)=pr(i,Nx);
        end 
        
       % SXX and SYY:
       
       for i=1:1:Ny
           for j =1:1:Nx
               SXX(i,j)=2*ETAP(i,j)*EXX(i,j);
               SYY(i,j)=2*ETAP(i,j)*EYY(i,j);
           end
       end
       
        %Calculate Shearheating
        for j=2:1:Nx
            for i=2:1:Ny
                Hs(i,j)=SXX(i,j)*EXX(i,j)+SYY(i,j)*EYY(i,j)+2*1/4*(SXY(i,j)*EXY(i,j)+SXY(i-1,j)*EXY(i-1,j)+SXY(i-1,j-1)*EXY(i-1,j-1)+SXY(i,j-1)*EXY(i,j-1));
            end
        end

        %11) Composing of global Matrixes LT() and RT()
        for j=1:1:Nx+1
            for i=1:1:Ny+1
                %Global index for 1 unknown Tdt
                gt=((j-1)*(Ny+1)+(i-1))*ukt+1;
                
                %horizontal boundaries: no slip BC:(T1+T2)/2=1573=T_BCi
                if(j>1 && j<Nx+1) %excludes corners
                    if(i==1)%upper boundary --> sticky air
                        %Left hand side
                        LT(gt,gt)=1;
                        LT(gt,gt+ukt)=1;
                        %Right hand side
                        RT(gt,1)=2*T_BCair;
                        
                    elseif(i==Ny+1)%lower boundary --> rest of mantle
                        %noslip
                        LT(gt,gt)=1;
                        LT(gt,gt-ukt)=1;
                        %Right hand side
                        RT(gt,1)=2*T_BCi;
                        
                    end
                end
                %vertical boundaries: free slip BC: (T1+T2)/2=0
                %includes corner
                
                if(j==1)
                    %Lefthandside
                    LT(gt,gt+(Ny+1)*ukt)=+1;%T2
                    LT(gt, gt)=-1;%T1
                    
                    %Righthandside
                    RT(gt,1)=0;
                    
                elseif(j==Nx+1)
                    %Lefthandside
                    LT(gt,gt)=1;%T2
                    LT(gt,gt-(Ny+1)*ukt)=-1;%T1
                    
                    %Righthandside
                    RT(gt,1)=0;
                    
                    
                end
                
                if(j>1 && i>1 && j <Nx+1 && i < Ny+1)
                    %inner points
                    % Temperature equaiton: RHOCP*DT/Dt=-dqx/dx-dqy/dy
                    %
                    %                    S(gt-ukt)
                    %                     Tdt2
                    %                      |
                    %                      o Ky(i-1,j)
                    %                      |
                    %         Kx(i,j-1)              Kx(i,j)
                    %Tdt1--------x----------Tdt3-------x---------Tdt5
                    %S(gt-(Ny+1)*ukt)      S(gt)               S(gt+(Ny+1)*ukt)
                    %                      |
                    %                      o Ky(i,j)
                    %                      |
                    %                     Tdt4
                    %                    S(gt+ukt)
                    %
                    %
                    %RHOCP*(T3-T0)/dt = (Kx(i,j)(T5-T3)/dx^2-Kx(i,j-1)(T3-T1)/dx^2 +
                    %Ky(i,j)(T4-T3)/dy^2 - Ky(i-1,j)(T3-T2)(dy^2)
                    %Left hand side (5 unknowns)
                    LT(gt, gt-(Ny+1)*ukt)=-Kx(i,j-1)/dx^2; %T1
                    LT(gt, gt-ukt)=-Ky(i-1,j)/dy^2;%T2
                    LT(gt, gt)= RHOCP(i,j)/dt+(Kx(i,j)+Kx(i,j-1))/dx^2+(Ky(i,j)+Ky(i-1,j))/dy^2-HaLH(i,j);%T3
                    LT(gt, gt+ukt)=-Ky(i,j)/dy^2;%T4
                    LT(gt, gt+(Ny+1)*ukt)=-Kx(i,j)/dx^2;%T5
                    %Right hand side               
                    RT(gt,1)=RHOCP(i,j)/dt*T0(i,j)+HaRH(i,j)+Hrad(i,j)+Hs(i,j); 
                end
            end
        end
        
        %12) Solving of global matrixes
        ST=LT\RT;
        
        %13) Reloading ST
        
        Tdt=zeros(Ny+1, Nx+1);
        
        for j=1:1:Nx+1
            for i=1:1:Ny+1
                %global index
                gt=((j-1)*(Ny+1)+(i-1))*ukt+1;
                
                %reload S(gt) into T(i,j)
                Tdt(i,j)=ST(gt);
            end
        end
        %Calculation of total adiabatic heating after solving for Tdt
        %(otherwise first Tdt=0)
        for j=2:1:Nx
            for i=2:1:Ny
                Ha(i,j)=HaLH(i,j)*Tdt(i,j)+HaRH(i,j);
            end
        end
        
        %14) Compute temperature change
        DT=Tdt-T0;
        
        %Do not change time step during the last iteration
        if(niter<nitermax)
            DTmax=max(max(abs(DT)));
            if(DTmax>Tchangemax) %if we have large temp change we need smaller time change
                dt=dt*Tchangemax/DTmax/1.1; %For safety we make dt 10% less than the allowed dt
                %dtthermal = dt
            end
            
        end
        
        
        
    end %Exit thermomechanical iteration --> makes code slower but more correct
   
    
    
        % 15)Visualise before moving markers
    figure(1)
    subplot(3,4,1); colormap('jet');
    pcolor(x,y,RHOvy) %Adjust y axis because its not a valid density anymore
   % hold on
    %[c,h]=contour(x_vy,y_vy(1:Ny),RHOvy(1:Ny,:),[-1 1650.5 5000],'m','LineWidth',2);
   % clabel(c,h,'LabelSpacing',1100,'FontSize',10,'Color','m');
   % hold off
    colorbar
    shading flat %flat because otherwise circle is asymmetrical with this low resolution
    axis ij image %Image size proportional to coordinates, vertical axis  upside down
    title("density, kg/m^3")

    subplot(3,4,2) ;colormap('jet');
    pcolor(x,y,log(ETAB))
    colorbar
    shading flat
    axis ij image %Image size proportional to coordinates, vertical axis  upside down
    title("Viscosity on normal grid points, Pa*s")

    subplot(3,4,3); colormap('jet');
    pcolor(x_p,y_p,T0)
    colorbar
    shading flat
    axis ij image %Image size proportional to coordinates, vertical axis  upside down
    title("Temperature, K")

    subplot(3,4,4);  colormap('jet');
    pcolor(x_p,y_p,pr)
    colorbar
    shading interp
    axis ij image %Image size proportional to coordinates, vertical axis  upside down
    title("Pressure, Pa")

    subplot(3,4,5);  colormap('jet');
    pcolor(x_vx,y_vx,vx)
    colorbar
    shading interp
    axis ij image %Image size proportional to coordinates, vertical axis  upside down
    title("vx-velocity, m/s")

    subplot(3,4,6);  colormap('jet');
    pcolor(x_vy,y_vy,vy)
    colorbar
    shading interp
    axis ij image %Image size proportional to coordinates, vertical axis  upside down
    title("vy-velocity, m/s")

    subplot(3,4,7);  colormap('jet');
    pcolor(x_p, y_p, RHOCP)
    colorbar
    shading flat
    axis ij image %Image size proportional to coordinates, vertical axis  upside down
    title("RHOCP, J/m^3/K")
    
    subplot(3,4,8);  colormap('jet');
    pcolor(x_vx,y_vx,log(Kx))
    colorbar
    shading flat
    axis ij image %Image size proportional to coordinates, vertical axis  upside down
    title("log(Kx),W/m/K")
    
    subplot(3,4,9);  colormap('jet');
    pcolor(x_p,y_p,alpha)
    colorbar
    shading flat
    axis ij image %Image size proportional to coordinates, vertical axis  upside down
    title("alpha, 1/K")
    
    subplot(3,4,10);  colormap('jet');
    pcolor(x_p,y_p,Hrad)
    colorbar
    shading flat
    axis ij image %Image size proportional to coordinates, vertical axis  upside down
    title("Hr, W/m^3")
    
    subplot(3,4,11);  colormap('jet');
    pcolor(x_p,y_p,Hs)
    colorbar
    shading flat
    axis ij image %Image size proportional to coordinates, vertical axis  upside down
    title("Hs, W/m^3")
    
    subplot(3,4,12);  colormap('jet');
    pcolor(x_p,y_p,Ha)
    colorbar
    shading flat
    axis ij image %Image size proportional to coordinates, vertical axis  upside down
    title("Ha, W/m^3")
    
    sgtitle(strcat('Timestep:', int2str(t)));
 
    %16) interpolation of eulerian temperature on to markers via DT
        for m=1:1:Nm
            %for vxm
            % Indexing:
            j = fix((xm(m)-x_p(1))/dx)+1;
            i = fix((ym(m)-y_p(1))/dy)+1;
           
            % Distances between Eulerian nodes and Lagrangian marker for
            % Tdt
            % nodal points
            
            dxmj= abs(xm(m)-x_p(j));
            dymi = abs(ym(m)-y_p(i));
            
            wtmij = (1-dxmj/dx)*(1-dymi/dy);
            wtmi1j =(1-dxmj/dx)*(dymi/dy);
            wtmij1 = (dxmj/dx)*(1-dymi/dy);
            wtmi1j1 = (dxmj/dx)*(dymi/dy);
            
            %Interpolate temperature
            if(t==1) %first temperature interpolation (other time steps in step 17))
                
                Tm(m) = Tdt(i,j)*wtmij + Tdt(i+1,j)*wtmi1j + Tdt(i,j+1)*wtmij1 + Tdt(i+1,j+1)*wtmi1j1;
         
            end 
           
        end 
    

    
    
    %17) Subgrid temperature diffusion
    
    %17.1) Compute Non-relaxation degree
    %for inner nodes
    for i=2:1:Ny
        for j=2:1:Nx
            TREij(i,j)=exp(-1/RHOCP(i,j)*((Kx(i,j)+Kx(i,j-1))/dx^2+(Ky(i,j)+Ky(i-1,j))/dy^2)*dt); 
        end
    end
    %External nodes - Symmetry condition 
     %excluding corners in x direction:
        for j=2:1:Nx
            TREij(1,j)=TREij(2,j);
            TREij(Ny+1,j)=TREij(Ny,j);
        end
        %including corners in y direction 
        for i=1:1:Ny+1
            TREij(i,1)=TREij(i,2);
            TREij(i,Nx+1)=TREij(i,Nx);
        end 
        
     %17.2) Calculate DT subgrid on marker
     
      %sums for DTm sub interpolation
      DTsubmwtsum = zeros(Ny+1,Nx+1);
      RHOCPsubwtsum = zeros(Ny+1,Nx+1);
   
     %Interpolate TRE on marker and T0 on marker 
     for m=1:1:Nm
         % Indexing for interpolations on pressure nodes
        j = fix((xm(m)-x_p(1))/dx)+1;
        i = fix((ym(m)-y_p(1))/dy)+1;
        
        % Distances between Eulerian nodes and Lagrangian marker for
        % Tdt
        % nodal points
        
        dxmj= abs(xm(m)-x_p(j));
        dymi = abs(ym(m)-y_p(i));
        
        wtmij = (1-dxmj/dx)*(1-dymi/dy);
        wtmi1j =(1-dxmj/dx)*(dymi/dy);
        wtmij1 = (dxmj/dx)*(1-dymi/dy);
        wtmi1j1 = (dxmj/dx)*(dymi/dy);
        
        %Interpolate TRE
        TREm(m)= TREij(i,j)*wtmij+TREij(i+1,j)*wtmi1j+ TREij(i,j+1)*wtmij1+TREij(i+1,j+1)*wtmi1j1;
        T0m(m)=T0(i,j)*wtmij+T0(i+1,j)*wtmi1j+ T0(i,j+1)*wtmij1+T0(i+1,j+1)*wtmi1j1;
    
        DTm0(m)=Tm(m)-T0m(m);
        
        DTm(m)=DTm0(m)*TREm(m);
        
        DTsubm(m)=DTm(m)-DTm0(m); %DTm subgrid
        
        
        %17.3) interpolate DTm subgrid onto nodal points ; weighted and heat balanced
        % Indexing for interpolations on pressure nodes
        j = fix((xm(m)-x_p(1))/dx)+1;
        i = fix((ym(m)-y_p(1))/dy)+1;
        
        % Distances between Eulerian nodes and Lagrangian marker for
        % Tdt
        % nodal points
        
        dxmj= abs(xm(m)-x_p(j));
        dymi = abs(ym(m)-y_p(i));
        
        wtmij = (1-dxmj/dx)*(1-dymi/dy);
        wtmi1j =(1-dxmj/dx)*(dymi/dy);
        wtmij1 = (dxmj/dx)*(1-dymi/dy);
        wtmi1j1 = (dxmj/dx)*(dymi/dy);
        %Define sums
        
        DTsubmwtsum(i,j) =DTsubmwtsum(i,j) + wtmij*DTsubm(m)*RHOCPm(m) ;%Heat conservative interpolation
        DTsubmwtsum(i+1,j) = DTsubmwtsum(i+1,j) + wtmi1j*DTsubm(m)*RHOCPm(m) ;
        DTsubmwtsum(i,j+1) =DTsubmwtsum(i,j+1) + wtmij1*DTsubm(m)*RHOCPm(m) ;
        DTsubmwtsum(i+1,j+1) = DTsubmwtsum(i+1,j+1) + wtmi1j1*DTsubm(m)*RHOCPm(m) ;
        
        RHOCPsubwtsum(i,j) = RHOCPsubwtsum(i,j) + wtmij*RHOCPm(m) ;
        RHOCPsubwtsum(i+1,j) = RHOCPsubwtsum(i+1,j) + wtmi1j*RHOCPm(m) ;
        RHOCPsubwtsum(i,j+1) = RHOCPsubwtsum(i,j+1) + wtmij1*RHOCPm(m);
        RHOCPsubwtsum(i+1,j+1) = RHOCPsubwtsum(i+1,j+1) + wtmi1j1*RHOCPm(m);
        
     end
     
     for j=1:1:Nx+1
         for i=1:1:Ny+1
             if(RHOCPsubwtsum(i,j)>0)
                 DTsubij(i,j)=DTsubmwtsum(i,j)/RHOCPsubwtsum(i,j);
             end
         end
     end
     
     %17.4) compute remaining DT
     DTremij=DT-DTsubij; 
     
     %17.5) interpolate DT remaining onto markers 
     for m=1:1:Nm
         % Indexing for interpolations on pressure nodes
         j = fix((xm(m)-x_p(1))/dx)+1;
         i = fix((ym(m)-y_p(1))/dy)+1;
         
         % Distances between Eulerian nodes and Lagrangian marker for
         % Tdt
         % nodal points
         
         dxmj= abs(xm(m)-x_p(j));
         dymi = abs(ym(m)-y_p(i));
         
         wtmij = (1-dxmj/dx)*(1-dymi/dy);
         wtmi1j =(1-dxmj/dx)*(dymi/dy);
         wtmij1 = (dxmj/dx)*(1-dymi/dy);
         wtmi1j1 = (dxmj/dx)*(dymi/dy);
         
         %Interpolate DT remaining
         DTremm(m)= DTremij(i,j)*wtmij+DTremij(i+1,j)*wtmi1j+ DTremij(i,j+1)*wtmij1+DTremij(i+1,j+1)*wtmi1j1;
   
        %17.6) Update Tm = Tm + DT subgrid + DT remaining 
        if(t>1)
            Tm(m)=Tm(m)+DTsubm(m)+DTremm(m);
        end
        
     end
     
      
        
    %18) Calculate temperature and pressure dependent viscosity and density
    for m=1:1:Nm
        % Indexing for interpolations on pressure nodes: Pm, EXXm, EYYm
        j = fix((xm(m)-x_p(1))/dx)+1;
        i = fix((ym(m)-y_p(1))/dy)+1;
        
        % Distances between Eulerian nodes and Lagrangian marker for
        % Tdt
        % nodal points
        
        dxmj= abs(xm(m)-x_p(j));
        dymi = abs(ym(m)-y_p(i));
        
        wtmij = (1-dxmj/dx)*(1-dymi/dy);
        wtmi1j =(1-dxmj/dx)*(dymi/dy);
        wtmij1 = (dxmj/dx)*(1-dymi/dy);
        wtmi1j1 = (dxmj/dx)*(dymi/dy);
        
        %Interpolate Pressure, EXX, EYY
        Pm(m)= pr(i,j)*wtmij+pr(i+1,j)*wtmi1j+ pr(i,j+1)*wtmij1+pr(i+1,j+1)*wtmi1j1;
        
        EXXm(m)= EXX(i,j)*wtmij+EXX(i+1,j)*wtmi1j+ EXX(i,j+1)*wtmij1+EXX(i+1,j+1)*wtmi1j1;
        EYYm(m)= EYY(i,j)*wtmij+EYY(i+1,j)*wtmi1j+ EYY(i,j+1)*wtmij1+EYY(i+1,j+1)*wtmi1j1;
        
        % Indexing for interpolations on basic nodes: EXYm
        j = fix((xm(m)-x(1))/dx)+1;
        i = fix((ym(m)-y(1))/dy)+1;
        
        % Distances between Eulerian nodes and Lagrangian marker for
        % Tdt
        % nodal points
        
        dxmj= abs(xm(m)-x(j));
        dymi = abs(ym(m)-y(i));
        
        wtmij = (1-dxmj/dx)*(1-dymi/dy);
        wtmi1j =(1-dxmj/dx)*(dymi/dy);
        wtmij1 = (dxmj/dx)*(1-dymi/dy);
        wtmi1j1 = (dxmj/dx)*(dymi/dy);
        
        %Interpolate Pressure, EXX, EYY
        EXYm(m)= EXY(i,j)*wtmij+EXY(i+1,j)*wtmi1j+ EXY(i,j+1)*wtmij1+EXY(i+1,j+1)*wtmi1j1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calculation of ETAm and RHOm for different types
        if(type(m)==1 || type(m)==2 ) %in plate and slab (incl. necking area)
            %ETA calculation
            EII(m) = ((EXXm(m)^2+EYYm(m)^2)/2+EXYm(m)^2)^0.5;
            ETAm(m) = 0.5/(ADd^(1/nd))*EII(m)^(1/nd-1)*exp((Ead+Pm(m)*Vad)/(8.314*Tm(m)*nd));
            %Strength test
            SII(m)= 2*EII(m)*ETAm(m);
            
            if(SII(m)>Sigmay(m)) %if deviatoric stress higher than yield strength reduce eta
                ETAm(m) = Sigmay(m)/(2*EII(m));
            end
            %ETA boundaries
            if(ETAm(m)>ETAmax)
                ETAm(m)=ETAmax;
            elseif(ETAm(m)<ETAmin)
                ETAm(m)=ETAmin;
            end
            
            %Density calculation
            RHOm(m) = RHO0d*(1+betam(m)*(Pm(m)-1e+5))/(1+alpham(m)*(Tm(m)-273));
            
            %update RHOCP
            RHOCPm(m) = RHOm(m)*CP_rock;
            
            %update K
            Km(m) = K0+a/(Tm(m)+b);
        end
  
        
        if(type(m)==3) %in mantle
            %ETA calculation
            EII(m) = ((EXXm(m)^2+EYYm(m)^2)/2+EXYm(m)^2)^0.5;
            ETAm(m) = 0.5/(ADw^(1/nw))*EII(m)^(1/nw-1)*exp((Eaw+Pm(m)*Vaw)/(8.314*Tm(m)*nw));
            %Strength test
            SII(m)= 2*EII(m)*ETAm(m);
            
            if(SII(m)>Sigmay(m)) %if deviatoric stress higher than yield strength reduce eta
                ETAm(m) = Sigmay(m)/(2*EII(m));
            end
            %ETA boundaries
            if(ETAm(m)>ETAmax)
                ETAm(m)=ETAmax;
            elseif(ETAm(m)<ETAmin)
                ETAm(m)=ETAmin;
            end
            
            %Density calculation
            RHOm(m) = RHO0w*(1+betam(m)*(Pm(m)-1e+5))/(1+alpham(m)*(Tm(m)-273));
            
            %update RHOCP
            RHOCPm(m) = RHOm(m)*CP_rock;
            
            %update K
            Km(m) = K0+a/(Tm(m)+b);
        end    
    end
    
    
    
    
    %  19) interpolation of eulerian velocity on to lagrangian marker
    for m=1:1:Nm
        
        %Classical 4.th Order Runge-Kutta advection Scheme
        xA=xm(m); %initial positions
        yA=ym(m);
        
        x_RK = xA; %points that change with every iteration to xA xB and xC
        y_RK= yA;
        
        for RK = 1:1:4 %four times for A B C D
            % Distances between Eulerian nodes and Lagrangian marker for
            % vx
            % nodal points
            
            %for vxm
            % Indexing:
            j = fix((x_RK-x_vx(1))/dx)+1;
            i = fix((y_RK-y_vx(1))/dy)+1;
            
            dxmj= abs(x_RK-x_vx(j));
            dymi = abs(y_RK-y_vx(i));
            
            wtmij = (1-dxmj/dx)*(1-dymi/dy);
            wtmi1j =(1-dxmj/dx)*(dymi/dy);
            wtmij1 = (dxmj/dx)*(1-dymi/dy);
            wtmi1j1 = (dxmj/dx)*(dymi/dy);
            
            vxm(RK)= vx(i,j)*wtmij + vx(i+1,j)*wtmi1j + vx(i,j+1)*wtmij1 + vx(i+1,j+1)*wtmi1j1;
            
            %Magical correction for continuity based velocity interpolation
            %
            
            %more explicit in two steps
            %if dxmj>dx/2
            %   if j<(Nx-1)
            %        vxcorr ...
            %    end
            %else
            %   if j>1
            %       vxcorr ...
            %   end
            %end
            if(dxmj>(dx/2) && j<(Nx-1)) %1. to the right 2. boundary at the right is not good
                corrvxX=0.5*(dxmj/dx-0.5)^2*((vx(i,j)-2*vx(i,j+1)+vx(i,j+2))*(1-dymi/dy)+(vx(i+1,j)-2*vx(i+1,j+1)+vx(i+1,j+2))*(dymi/dy));
                vxm(RK)=vxm(RK)+corrvxX;
                %
            elseif(dxmj <= (dx/2)&& j>1)%1. to the left 2. boundary at the left  is not good
                corrvxY=0.5*(dxmj/dx-0.5)^2*((vx(i,j-1)-2*vx(i,j)+vx(i,j+1))*(1-dymi/dy)+(vx(i+1,j-1)-2*vx(i+1,j)+vx(i+1,j+1))*(dymi/dy));
                vxm(RK)=vxm(RK)+corrvxY;
            end
            
            %for vym
            % Indexing:
            j = fix((x_RK-x_vy(1))/dx)+1;
            i = fix((y_RK-y_vy(1))/dy)+1;
            
            % Distances between Eulerian nodes and Lagrangian marker for
            % vy
            % nodal points
            
            dxmj= abs(x_RK-x_vy(j));
            dymi = abs(y_RK-y_vy(i));
            
            wtmij = (1-dxmj/dx)*(1-dymi/dy);
            wtmi1j =(1-dxmj/dx)*(dymi/dy);
            wtmij1 = (dxmj/dx)*(1-dymi/dy);
            wtmi1j1 = (dxmj/dx)*(dymi/dy);
            
            vym(RK)= vy(i,j)*wtmij + vy(i+1,j)*wtmi1j + vy(i,j+1)*wtmij1 + vy(i+1,j+1)*wtmi1j1;
            
            %Magical correction for continuity based velocity interpolation
            
            if(dymi>(dy/2) && i<(Ny-1)) %1. to the right 2. boundary at the right is not good
                corrvyX=0.5*(dymi/dy-0.5)^2*((vy(i,j)-2*vy(i+1,j)+vy(i+2,j))*(1-dxmj/dx)+(vy(i,j+1)-2*vy(i+1,j+1)+vy(i+2,j+1))*(dxmj/dx));
                vym(RK)=vym(RK)+corrvyX;
                
            elseif(dymi <= (dy/2)&& i>1) %1. to the left 2. boundary to the left is not good
                corrvyY=0.5*(dymi/dy-0.5)^2*((vy(i-1,j)-2*vy(i,j)+vy(i+1,j))*(1-dxmj/dx)+(vy(i-1,j+1)-2*vy(i,j+1)+vy(i+1,j+1))*(dxmj/dx));
                vym(RK)=vym(RK)+corrvyY;
            end
            
            %Runge Kutta Adveciton for velocity interpolation
            %Point A: --> pre-defined before RK loop
            %xA=xm(m)
            %yA=ym(m)
            %
            %Point B:
            %xB=xA+dt/2*vxmA
            %yB=yA+dt/2*vymA
            %
            %Point C:
            %xC=xA+dt/2*vxmB
            %yC=yA+dt/2*vymB
            %
            %Point D: 
            %xD=xA+ dt*vxmC
            %yD=yA+dt*vymC
            
            %points B and C --> B calc w/ vel A RK==1, C calc w/ vel B RK
            %==2
            if(RK<3)%==1 ||RK==2)
                x_RK=xA+dt/2*vxm(RK);
                y_RK=yA+dt/2*vym(RK);
                
            end
            %point D calculated with velocity C --> RK==3
            if(RK==3)
                x_RK=xA+dt*vxm(RK);
                y_RK=yA+dt*vym(RK);
            end 
             
        end
        %Calculate effective marker velocities after 4 iterations of interpolation
        %vxmeff = 1/6*(vxmA+2vxmB+2vxmC+vxmD)
        vxmeff=1/6*(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4));% or maybe other indexing w/m m+1, m+2, m+3, m+4
        %vymeff = 1/6*(vymA+2vymB+2vxmC+vxmD)
        vymeff=1/6*(vym(1)+2*vym(2)+2*vym(3)+vym(4)); 
      
        %18) moving marker 
        xm(m) = xm(m) + vxmeff*dt;
        ym(m) = ym(m) + vymeff*dt;
        
        
           
        
    end
   timesum = timesum+dt;
    
end
    %Controll values
format long




aaa(1,1)=pr(27,12);
aaa(2,1)=vx(27,12);
aaa(3,1)=vy(27,12);
aaa(4,1)=Tdt(27,12) ;
aaa(5,1)=dt; 
aaa(6,1)=timesum ;
aaa(7,1)=Hrad(27,12);
aaa(8,1)=Hs(27,12);
aaa(9,1)=Ha(27,12);

