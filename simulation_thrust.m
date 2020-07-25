clc;
    D = 0.037;           %stator width > taken from LIM drawing: 0.457 mm x 81 = 37.017 mm = 0.037 m
    p = 4;               %number of poles > taken from spreadsheet
    dr = 0.01049;        %rotor plate thickness > from spaceX drawing 0.413+-0.008 in = 0.01049 m
    sigma = 2.5*10^7;  %conductivity of the aluminum I-beam/ or rotor plate > taken from reference = 2.5x10^7 siemens/m (SI units)
    s = 0.1;             %slip > taken from spreadsheet
    vr = 40;             %mechanical speed of the rotor > SS desired velocity from spreadsheet
    vs = vr/(1-s);       %stator sync speed > SS desired velocity
    N = 54;               %number of turns of the windings > taken from spreadsheet (9 turns per slot * 6 parallel wires)
    I = 15;              %steady state current > taken from FDP
    g = (2*0.003)+dr;     %air gap distance = 2*air gap + thickness of the track = 2 * 0.003 + rotor plate thickness > 0.003 m taken from spreadsheet
    mu = (4*pi*10^-7)*(5000);           %mu0*muR=mu > muR is the permeability of the iron core > (4*pi*10^-7 H/m)*(5000 H/m)
    f = 120;             %source frequency > taken from spreadsheet
    L = 0.7366;              %length of stator > taken from LIM drawing (29 inches = 0.7366 m)
    tau = L/(2*p);       %pole pitch > taken from Donovan's drawing (24 slots / 3 phases | 3 coils | 3 slots | 3 phases make one pole | 3.625
    %w = 2*pi*f;          %angular frequency given by w=2*pi/period, w=2*pi*f
    B1 = (3*sqrt(2)*N*I)/sqrt((pi*g/mu)^2+((vs-vr)^2*(sigma*dr*tau)^2));      %amplitude of the magnetic flux density
    %phi = atan((pi*g)/(mu*sigma*dr*tau*(vs-vr)));                             %phase shift of the magnetic flux density
    %y                %position along the y-axis (direction of linear travel)
    %Jz = -sigma*(vs-vr)*B1*e^(1i*((pi*y/tau)-(w*t)+phi));                     %induced eddy currents
    a = sigma*(vs-vr)*B1;

Ft = D*dr*a*B1*L