clc
clear all
format short g
%% Start of Code   (L=400 , e1=20 , e2=5 , e3=10)
%Input Parameters
%mirror parameters
w=23.18;%width of mirrors
I_sun=518;%sun radiation=I_sun
T_amb=296;%ambient temperature
f=2;%focal distance
betta=0.07;%concentrate coefficient=1-(parabola specular reflectance)
%environment parameters (Air at 296~300 K)
k_e=26.3e-3;%thermal conductivity of environment=k_e(modifed from white-fluid mechanics)
betta_e=0.00337;%volumetric thermal expansion coefficient for environment
v_e=15.89e-6;%kinematic viscisity for environment
alpha_e=22.5e-6;%thermal diffusivity for environment
C_p_e=1007;%specific heat at constant pressure for environment
mu_e=184.6e-7;%viscosity for environment
ro_e=1.1614;%Desity of Air at Environment Temperature
%glass (pyrex) parameters
k_g=1.4;%thermal conductivity of glass=k_g
w_c=0.6;%width of glass
t_g=0.01;%thickness of glass
gamma=0.95;%coefficient of glass transient=gamma
alpha_g=0.05;%absorptivity coefficient of glass
epsilon_g=0.82;%emmisivity of glass
%imprisoned air parameters (Air at 400 K)
k_c=33.8e-3;%thermal conductivity of cavity=k_c
betta_c=0.0025;%volumetric thermal expansion coefficient for cavity
v_c=26.41e-6;%kinematic viscisity for cavity
alpha_c=38.3e-6;%thermal diffusivity for cavity
C_p_c=1014;%specific heat at constant pressure for cavity
mu_c=230.1e-7;%viscosity for cavity
w_s=0.49;%surface's width for bottom of pipes
D=0.1;%heigth of cavity
%absorber parameters
D_a_in=0.084;%internal diameter of absorber=D_a_in
D_a_ex=0.098;%external diameter of absorber=D_a_ex
r_a_ex=D_a_ex/2;
alpha_a=0.96;%absorptivity coefficient of absorber
epsilon_a=0.22;%emmisivity of absorber=epsilon_a
k_a=126.7;%thermal conductivity of absorber=k_a
N=input('number of pipes:');%number of pipes
%Nano Particle Parameter
fi=input('Volume fraction of nano particle:');%Volume Fraction of Nano Particle
ro_p=19300;%Density of Nano Particle
ro_f_ref=999;%Fluid Density in the 20 degree celcius
d_p=67.5e-9;%Diameter of Nano Particle
k_p=310;%Thermal Conductivity of Nano Particle
C_p_p=129;%Specific Heat for Nano Paricle
n=3;%empirical shape factor of Nano Particle
Av=6.022*(10e23);%Avogadro
M=0.0180153;%Mole Mass of base Fluid
Anp=-4.34023;
Bnp=-105.60799;
Cnp=527.45631;
Dnp=-151.74505;
Enp=-903.41949;
Fnp=2814.30560;
alphanp=-0.00355;
bettanp=0.03564;
gammanp=-0.10898;
sigmanp=0.08845;
epsilonnp=0.02737;
omeganp=-0.04025;
%general parameters
pi=3.14;%pi number=3.14
g=9.81;%gravitational acceleration=9.81
sigma=0.000000056;%Stefan–Boltzmann constant (sigma)=5.67*10^-8
Cr_g=(f/(D+t_g-r_a_ex));
Cr_a=(f/D_a_ex);
%velosity of wind=v
v_air=input('value of wind velosity:');
if v_air<0.1 % when value of wind velosity < 0.1
    %% Start of Calculations for Preheating Section
    %fluid parameters for moderate temperature of 315 (K)
    L=input('value of pipe length:');%length of absorber at Pre Heating Section
    n1=input('please input number of elements for absorber:');
    Length1=L/n1;% length of each absorber element at Pre Heating Section (meter)
    V_f1=0.1;%mass rate of fluid
    I_11=I_sun*betta*alpha_g*Cr_g;% eq 2
    I_21=(I_sun*betta*gamma*alpha_a*Cr_a)/N;% eq 5
    T_out1=zeros(n1,1);
    T_a_ex1=zeros(n1,1);
    T_a_in1=zeros(n1,1);
    T_g_in1=zeros(n1,1);
    T_g_ex1=zeros(n1,1);
    Nu_e1=zeros(n1,1);
    Nu_c1=zeros(n1,1);
    h_nf1=zeros(n1,1);
    %input temperature=T_i
    T_in1=input('Value of Fluid Temperature at Start Process:');%ok
    for i=1:n1 % quantity of elements
        % initial temperatures= Tstar_a_ex, Tstar_g_ex, Tstar_g_in, T_g_ex
        Tstar_a_ex1=360;
        Tstar_g_in1=345;
        Tstar_g_ex1=340;
        while (abs(T_a_ex1-Tstar_a_ex1))>0.01 & (abs(T_g_in1-Tstar_g_in1))>0.01 & (abs(T_g_ex1-Tstar_g_ex1))>0.01
            mu_f1=((2.1897e-11)*(T_in1.^4))-((3.055e-8)*(T_in1.^3))+((1.6028e-5)*(T_in1.^2))-(0.0037524*T_in1)+0.33158;
            C_p_f1=((1.1105e-5)*(T_in1.^3))-(0.0031078*(T_in1.^2))-(1.478*T_in1)+4631.9;
            k_f1=((1.5362e-8)*(T_in1.^3))-((2.261e-05)*(T_in1.^2))+(0.010879*T_in1)-1.0294;
            ro_f1=((-1.5629e-5)*(T_in1.^3))+(0.011778*(T_in1.^2))-(3.0726*T_in1)+1227.8;
            ro_nf1=(ro_p*fi)+(ro_f1*(1-fi));
            A_i1=pi*((D_a_in/2).^2);
            mdat1=V_f1*A_i1*ro_nf1;
            d_f=0.1*((6*M)/(Av*pi*ro_f_ref));
            mu_nf1=mu_f1/(1-(34.87*((d_p/d_f).^(-0.3))*(fi.^1.03)));
            k_nf1=k_f1*((k_p+((n-1)*k_f1)-((n-1)*fi*(k_f1-k_p)))/(k_p+((n-1)*k_f1)+(fi*(k_f1-k_p))));
            C_p_nf1=((fi*ro_p*C_p_p)+((1-fi)*ro_f1*C_p_f1))/ro_nf1;
            Pr_nf1=(C_p_nf1*mu_nf1)/k_nf1;
            Re_nf1=(ro_nf1*V_f1*D_a_in)/(mu_nf1);
            Nu_nf1=(0.023)*(Re_nf1.^(4/5))*(Pr_nf1.^0.4);
            A11=(w_c*k_g)/t_g;% From eq 3
            A21=(2*pi*k_a)/(log(D_a_ex/D_a_in));% From eq 6
            h_nf1=(Nu_nf1*k_nf1)/D_a_in;
            A31=h_nf1*pi*D_a_in;% From eq 7
            A41=(mdat1*C_p_f1)/Length1;% From eq 19
            Nu_c1=0.52*(((g*betta_c*(D_a_ex.^3))/(v_c*alpha_c)).^0.2)*((Tstar_a_ex1-Tstar_g_in1).^0.2);% eq 34
            h_c1=(Nu_c1*k_c)/D;% From eq 33
            A51=(h_c1*Length1)/N;% From eq 32
            Nu_e1=0.15*(((g*betta_e*(w_c.^3))/(v_e*alpha_e)).^(1/3))*((Tstar_g_ex1-T_amb).^(1/3));% eq 27
            h_e1=(k_e*Nu_e1)/w_c;% From eq 26
            A61=h_e1*w;% From eq 25
            b11=I_11*w;% From eq 1
            b21=0;% From eq 2
            b31=I_21*w;% From eq 3
            b41=0;% From eq 4
            b51=-(1/2)*A31*T_in1;% From eq 5
            b61=(-A41)*T_in1;% From eq 6
            b71=0;% From eq 7
            b81=(-A61)*T_amb;% From eq 8
            b91=0;% From eq 9
            b101=0;% From eq 10
            b111=0;% From eq 11
            b121=0;% From eq 12
            A1=[1,0,0,0,0,0,0,0,0,0,0,0% eq 1
                 0,1,0,0,0,0,0,-A11,A11,0,0,0% eq 2
                 0,0,1,0,0,0,0,0,0,0,0,0% eq 3
                 0,0,0,1,0,0,0,0,0,-A21,A21,0% eq 4
                 0,0,0,0,1,0,0,0,0,0,-A31,(A31)/2% eq 5
                 0,0,0,0,1,0,0,0,0,0,0,-A41% eq 6
                 0,0,0,0,0,1,0,0,A51,-A51,0,0% eq 7
                 0,0,0,0,0,0,1,-A61,0,0,0,0% eq 8
                 0,1,0,0,0,1,0,0,0,0,0,0% eq 9
                 1,0,0,0,0,N,-1,0,0,0,0,0% eq 10
                 0,0,1,0,-1,-1,0,0,0,0,0,0% eq 11
                 0,0,0,1,-1,0,0,0,0,0,0,0];% eq 12
                 B1=[b11;b21;b31;b41;b51;b61;b71;b81;b91;b101;b111;b121];
                 %X=[qdat_g_srad;qdat_g_cond;qdat_a_srad;qdat_a_cond;qdat_af_conv;qdat_ag_conv;q_ge_conv;T_g_ex;T_g_in;T_a_ex;T_a_in;T_out];
                 X1=A1\B1;
                 T_a_ex1=X1(10,:);
                 T_g_in1=X1(9,:);
                 T_g_ex1=X1(8,:);
                 Delta11=(T_a_ex1-Tstar_a_ex1);
                 Delta21=(T_g_in1-Tstar_g_in1);
                 Delta31=(T_g_ex1-Tstar_g_ex1);
                 Tstar_a_ex1=Tstar_a_ex1+0.5*(Delta11);
                 Tstar_g_in1=Tstar_g_in1+0.5*(Delta21);
                 Tstar_g_ex1=Tstar_g_ex1+0.5*(Delta31);
        end       
        T_out1(i)=X1(12,:);
        T_in1=T_out1(i);
        TEMP11(i)=X1(12,:);
        TF1=TEMP11';%TF : Temperature of Fluid
        T_a_in1=X1(11,:);
        TEMP12(i)=X1(11,:);
        TIS1=TEMP12';%TIS : Temperature of Internal Surface for pipe
        Flux1(i)=X1(5,:);
        HFaf1=Flux1';%HFaf : Convection Heat Flux from Absorber to Fluid
        Nuselt_c1(i)=Nu_c1;
        Nuselt_c1=Nuselt_c1';
        Nuselt_e1(i)=Nu_e1;
        Nuselt_e1=Nuselt_e1';
        HT1(i)=h_nf1; 
        HT11=HT1';
        if T_a_in1>373.15
            break
        end
    end
    disp 'Distribution of Internal Surface Temperature at Pre Heating Section is:';
    disp (TIS1);
    y11=TIS1;
    sizey11=size(y11);
    s=sizey11(1,1);
    a11=linspace(L/n1,((L/n1)*s),s)';% because i have not first point so start from L/n1; n1 = numbber of the distances
    c11=polyfit(a11,y11,s-1);%because i have not first point so degree of polynomial is n1-1
    disp 'The Internal Surface Temperature at Start of Preheating Section:'
    T_a_in01=polyval(c11,0)%Temperature of internal surface at start of preheating section
    sizec11=size(c11);
    C1=sizec11(1,2);%for determine number of members
    DELTA1=c11(1,C1)-373.15;%for determine equation of TIS1 at 373.15
    c11(:,C1)=[];%for determine equation of TIS1 at 373.15
    d1=[c11,DELTA1];%equation of TIS1 at 373.15
    q1=roots(d1);
    z1=size(q1);
    disp 'The Length of Preheating Section:'
    j=1:1:z1(1,1);
    x1=zeros(j,1);
    for j=1:1:z1(1,1)
        M1=q1(j);
        N1=imag(q1(j));
        if M1>0 & M1<L & N1==0
            x1(j)=M1;
        end
    end
    x11=nonzeros(x1);
    l1=min(x11)
    disp 'Distribution of Fluid Temperature at Pre Heating Section is:';
    disp (TF1);
    a12=linspace(0,(L/n1)*s,s+1)';
    y12=[303.15;TF1];
    c12=polyfit(a12,y12,s);
    disp 'The Fluid Temperature at End of Preheating Section:'
    T1=polyval(c12,l1)
    disp 'The Different Temperature of Internal Surface & Fluid at End of Preheating Section:'
    deltaT1=373.15-T1
    disp 'Distribution of Convection Heat Flux from Absorber to Fluid at Preheating Section is:';
    disp (HFaf1);
    a13=linspace(L/n1,L,s)';
    y13=HFaf1;
    c13=polyfit(a13,y13,s-1);
    disp 'The Heat Flux at Start of Preheating Section:';
    HF10=polyval(c13,0);
    disp 'The Heat Flux at End of Preheating Section:';
    HF11=polyval(c13,l1);
    %% Start of Calculation for Subcooled Section
    %fluid parameters for moderate temperature of 355 (K)
    hlv2=2317000;%input('value of enthalpy at subcooled section:');%Enthalpy of Fluid at Subcooled Section (j/kg)
    deltaPsat2=1;%input('value of saturated pressure at subcooled section:');%Vapor-Pressure difference corresponding to superheat temperature
    Psat2=0.4163;
    V_f2=input('Please input value of fluid velocity in sub-cooled section:');%fluid velosity=V
    L1=L-l1;%length of absorber at Subcoled Section
    n2=input('please input number of elements for absorber:');
    Length2=L1/n2;%length of each absorber element at Subcoled Section (meter)
    I_12=I_sun*betta*alpha_g*Cr_g;% eq 2
    I_22=(I_sun*betta*gamma*alpha_a*Cr_a)/N;% eq 5
    T_out2=zeros(n2,1);
    T_a_ex2=zeros(n2,1);
    T_a_in2=zeros(n2,1);
    T_g_ex2=zeros(n2,1);
    T_g_in2=zeros(n2,1);
    Nu_c2=zeros(n2,1);
    Nu_e2=zeros(n2,1);
    h_sc=zeros(n2,1);
    %input temperature=T_i
    T_in2=T1;% (Temperature for end of Pre Heating Section)
    for i=1:n2 % quantity of elements
        % initial temperatures= Tstar_a_ex, Tstar_g_ex, Tstar_g_in, T_g_ex
        Tstar_a_ex2=510;
        Tstar_a_in2=500;
        Tstar_g_in2=407.485;
        Tstar_g_ex2=420;
        Tstar_out2=473.15;
        while (abs(T_a_in2-Tstar_a_in2))>0.01 & (abs(T_out2-Tstar_out2))>0.01 & (abs(T_a_ex2-Tstar_a_ex2))>0.01 & (abs(T_g_in2-Tstar_g_in2))>0.01 & (abs(T_g_ex2-Tstar_g_ex2))>0.01
            mu_f2=((2.1897e-11)*(T_in2.^4))-((3.055e-8)*(T_in2.^3))+((1.6028e-5)*(T_in2.^2))-(0.0037524*T_in2)+0.33158;
            C_p_f2=((1.1105e-5)*(T_in2.^3))-(0.0031078*(T_in2.^2))-(1.478*T_in2)+4631.9;
            k_f2=((1.5362e-8)*(T_in2.^3))-((2.261e-05)*(T_in2.^2))+(0.010879*T_in2)-1.0294;
            ro_f2=((-1.5629e-5)*(T_in2.^3))+(0.011778*(T_in2.^2))-(3.0726*T_in2)+1227.8;
            ro_v2=0.260;
            ro_nf2=(ro_p*fi)+(ro_f2*(1-fi));
            A_i2=pi*((D_a_in/2).^2);
            mdat2=V_f2*A_i1*ro_nf2;
            d_f=0.1*((6*M)/(Av*pi*ro_f_ref));
            mu_nf2=mu_f2/(1-(34.87*((d_p/d_f).^(-0.3))*(fi.^1.03)));
            k_nf2=k_f2*((k_p+((n-1)*k_f2)-((n-1)*fi*(k_f2-k_p)))/(k_p+((n-1)*k_f2)+(fi*(k_f2-k_p))));
            C_p_nf2=((fi*ro_p*C_p_p)+((1-fi)*ro_f2*C_p_f2))/ro_nf2;
            A12=(w_c*k_g)/t_g;% From eq 3
            A22=(2*pi*k_a)/(log(D_a_ex/D_a_in));% From eq 6
            if fi==0
                hlv_nf2=2317000;
            else
                C02=(Anp*(fi.^5))+(Bnp*(fi.^4))+(Cnp*(fi.^3))+(Dnp*(fi.^2))+(Enp*fi)+Fnp;
                C12=(alphanp*(fi.^5))+(bettanp*(fi.^4))+(gammanp*(fi.^3))+(sigmanp*(fi.^2))+(epsilonnp*fi)+omeganp;
                hlv_nf2=C02*(deltaPsat2.^C12);
            end
            Pr_nf2=(C_p_nf2*mu_nf2)/k_nf2;% eq 10
            Re_nf2=(ro_f2*V_f2*D_a_in)/(mu_nf2);% eq 11
            Nu_nf2=(0.023)*(Re_nf2.^(4/5))*(Pr_nf2.^0.4);% eq 9
            h_nf2=(Nu_nf2*k_nf2)/D_a_in;% eq 8
            G=0.00122*(((k_nf2.^0.79)*(C_p_nf2.^0.45)*(ro_f2.^0.49))/((sigma.^0.5)*(mu_nf2.^0.29)*(hlv_nf2.^0.24)*(ro_v2.^0.24)))*(deltaPsat2.^0.75);% % From eq 14
            S_chen=1/(1+((2.53*(10.^-6))*(Re_nf2.^1.17)));
            h_nb=G*((Tstar_a_in2-Tstar_out2).^0.24);% eq 14
            A32=h_nf2*pi*D_a_in;% From eq 13
            A332=h_nb*S_chen*pi*D_a_in;
            A42=(mdat2*C_p_f2)/Length2;% From eq 19
            Nu_c2=0.52*(((g*betta_c*(D_a_ex.^3))/(v_c*alpha_c)).^0.2)*(Tstar_a_ex2-Tstar_g_in2).^0.2;% eq 34
            h_c2=(Nu_c2*k_c)/D;% eq 33
            A52=(h_c2*Length2)/N;% From eq 32
            Nu_e2=0.15*(((g*betta_e*(w_c.^3))/(v_e*alpha_e)).^(1/3))*((Tstar_g_ex2-T_amb).^(1/3));% eq 27
            h_e2=(k_e*Nu_e2)/w_c;% eq 26
            A62=h_e2*w;% From eq 25
            b12=((I_12)*w);% From eq 1
            b22=0;% From eq 2
            b32=((I_22)*w);% From eq 3
            b42=0;% From eq 4
            b52=-373.15*A332;% T,sat = T,out So : b5=0----% From eq 5
            b62=(-A42)*T_in2;% From eq 6
            b72=0;% From eq 7
            b82=(-A62)*T_amb;% From eq 8
            b92=0;% From eq 9
            b102=0;% From eq 10
            b112=0;% From eq 11
            b122=0;% From eq 12
            A2=[1,0,0,0,0,0,0,0,0,0,0,0% eq 1
                 0,1,0,0,0,0,0,-A12,A12,0,0,0% eq 2
                 0,0,1,0,0,0,0,0,0,0,0,0% eq 3
                 0,0,0,1,0,0,0,0,0,-A22,A22,0% eq 4
                 0,0,0,0,1,0,0,0,0,0,-(A32+A332),A32% eq 5
                 0,0,0,0,1,0,0,0,0,0,0,-A42% eq 6
                 0,0,0,0,0,1,0,0,A52,-A52,0,0% eq 7
                 0,0,0,0,0,0,1,-A62,0,0,0,0% eq 8
                 0,1,0,0,0,1,0,0,0,0,0,0% eq 9
                 1,0,0,0,0,N,-1,0,0,0,0,0% eq 10
                 0,0,1,0,-1,-1,0,0,0,0,0,0% eq 11
                 0,0,0,1,-1,0,0,0,0,0,0,0];% eq 12
                 B2=[b12;b22;b32;b42;b52;b62;b72;b82;b92;b102;b112;b122];
                 %X=[qdat_g_srad;qdat_g_cond;qdat_a_srad;qdat_a_cond;qdat_af_conv;qdat_ag_conv;q_ge_conv;T_g_ex;T_g_in;T_a_ex;T_a_in;T_out];
                 X2=A2\B2;
                 T_a_ex2=X2(10,:);
                 T_a_in2=X2(11,:);
                 T_g_in2=X2(9,:);
                 T_g_ex2=X2(8,:);
                 T_out2=X2(12,:);
                 Delta12=(T_g_ex2-Tstar_g_ex2);
                 Delta22=(T_a_ex2-Tstar_a_ex2);
                 Delta32=(T_g_in2-Tstar_g_in2);
                 Delta42=(T_a_in2-Tstar_a_in2);
                 Delta52=(T_out2-Tstar_out2);
                 Tstar_g_ex2=Tstar_g_ex2+0.5*(Delta12);
                 Tstar_a_ex2=Tstar_a_ex2+0.5*(Delta22);
                 Tstar_g_in2=Tstar_g_in2+0.5*(Delta32);
                 Tstar_a_in2=Tstar_a_in2+0.5*(Delta42);
                 Tstar_out2=Tstar_out2+0.5*(Delta52);
        end       
        T_out2(i)=X2(12,:);
        T_in2=T_out2(i);
        TEMP21(i)=X2(12,:);
        TF2=TEMP21';%TF : Temperature of Fluid
        TEMP22(i)=X2(11,:);
        TIS2=TEMP22';%TIS : Temperature of Internal Surface
        Flux2(i)=X2(5,:);
        HFaf2=Flux2';%HFaf : Convection Heat Flux from Absorber to Fluid
        Nuselt_c2(i)=Nu_c2;
        Nuselt_c2=Nuselt_c2';
        Nuselt_e2(i)=Nu_e2;
        Nuselt_e2=Nuselt_e2';
        HT2(i)=h_nb;
        HT22=HT2';
        if T_out2>373.15
            break
        end
    end
    disp 'Distribution of Fluid Temperature at Subcooled Section is:';
    disp (TF2);
    y21=[T1;TF2];
    sizey21=size(y21);
    s=sizey21(1,1);
    a21=linspace(0,(L1/n2)*s,s)';
    c21=polyfit(a21,y21,s);
    disp 'The Fluid Temperature at Start of Subcooled Section:'
    T_in02=c21(1,s+1)
    sizec21=size(c21);
    C2=sizec21(1,2);% for determine number of members
    DELTA2=c21(1,C2)-373.15;% for determine equation of TF2 at 373.15
    c21(:,C2)=[];% for determine equation of TF2 at 373.15
    d2=[c21,DELTA2];% equation of TF2 at 373.15
    q2=roots(d2);
    z2=size(q2);
    disp 'The Length of Subcooled Section:'
    j=1:1:z2(1,1);
    x2=zeros(j,1);
    for j=1:1:z2(1,1)
        M2=q2(j);
        N2=imag(q2(j));
        if M2>0 & M2<L1 & N2==0
            x2(j)=M2;
        end
    end
    x22=nonzeros(x2);
    l2=min(x22)
    disp 'Distribution of Internal Surface Temperature at Subcooled Section is:';
    disp (TIS2);
    a22=linspace(0,(L1/n2)*s,s)';
    y22=[373.15;TIS2];
    size22=size(y22);
    c22=polyfit(a22,y22,s);%change size... with n2
    disp 'The Internal Surface Temperature at Start of Subcooled Section:'
    T_a_in02=polyval(c22,0)
    disp 'The The Internal Surface Temperature at End of Subcooled Section:'
    T2=polyval(c22,l2)
    disp 'The Different Temperature of Internal Surface & Fluid at End of Subcooled Section:'
    deltaT2=T2-373.15
    disp 'Distribution of Convection Heat Flux from Absorber to Fluid at Subcooled Section is:';
    disp (HFaf2);
    a23=linspace(0,(L1/n2)*s,s)';
    y23=[HF11;HFaf2];
    size23=size(y23);
    c23=polyfit(a23,y23,s);
    disp 'The Heat Flux at Start of Preheating Section:';
    HF20=polyval(c23,0)
    disp 'The Heat Flux at End of Preheating Section:';
    HF21=polyval(c23,l2);
    %% Start of Calculation for Two-Phase Flow Section
    C_p_f3=4217;%specific heat at constant pressure=C_p
    mu_f3=279e-006;%viscosity=mu
    k_f3=680e-003;%thermal conductivity of fluid=k_f
    ro_l3=957.85;% water density
    ro_v3=0.595;% vapor density
    ro_lv3=0.598;%vapor-liquid density
    V_lv3=1.67185;
    V_f3=input('Please input value of fluid velocity in two-phase section:');%fluid velosity=V
    l3=L-(l1+l2);%length of absorber
    ro_nf3=(ro_p*fi)+(ro_l3*(1-fi));
    d_f=0.1*((6*M)/(Av*pi*ro_f_ref));
    mu_nf3=mu_f3/(1-(34.87*((d_p/d_f).^(-0.3))*(fi.^1.03)));
    k_nf3=k_f3*((k_p+((n-1)*k_f3)-((n-1)*fi*(k_f3-k_p)))/(k_p+((n-1)*k_f3)+(fi*(k_f3-k_p))));
    C_p_nf3=((fi*ro_p*C_p_p)+((1-fi)*ro_l3*C_p_f3))/ro_nf3;
    Pr_nf3=(C_p_nf3*mu_nf3)/k_nf3;
    Re_nf3=(ro_nf3*V_f3*D_a_in)/(mu_nf3);
    Nu_nf3=(0.023)*(Re_nf3.^(4/5))*(Pr_nf3.^0.4);
    h_l3=(Nu_nf3*k_nf3)/D_a_in;
    deltaPsat3=1;
    pr3=1;% eq 18
    n3=input('please input number of elements for absorber:');
    Length3=l3/n3;%length of each absorber element (meter)
    I_13=I_sun*betta*alpha_g*Cr_g;% eq 2
    I_23=(I_sun*betta*gamma*alpha_a*Cr_a)/N;% eq 5
    x_out3=zeros(n3,1);
    T_a_ex3=zeros(n3,1);
    T_g_ex3=zeros(n3,1);
    T_g_in3=zeros(n3,1);
    Nu_c3=zeros(n3,1);
    Nu_e3=zeros(n3,1);
    h_tp3=zeros(n3,1);
    %input temperature=T_i
    x_in3=0;% (vapor quality for end of Subcooled Section)
    for i=1:n3 % quantity of elements
        % initial temperatures= Tstar_a_ex, Tstar_g_ex, Tstar_g_in, T_g_ex
        Tstar_a_ex3=400;
        Tstar_g_in3=390;
        Tstar_g_ex3=395;
        xstar_out3=0.9;%25.12.2016;18.01.2017;delete of T_out and add the x=vapore quality
        while (abs(x_out3-xstar_out3))>0.0001 & (abs(T_a_ex3-Tstar_a_ex3))>0.01 & (abs(T_g_in3-Tstar_g_in3))>0.01 & (abs(T_g_ex3-Tstar_g_ex3))>0.01
            A13=(w_c*k_g)/t_g;% From eq 3
            A23=(2*pi*k_a)/(log(D_a_ex/D_a_in));% From eq 6
            h_tp3=((h_l3*(pr3.^0.38)*((1-xstar_out3).^0.8))+(3.8*(xstar_out3.^0.76)*((1-xstar_out3).^0.04)))*(pr3.^(-0.38));% eq 16
            A33=h_tp3*pi*D_a_in;% From eq 15
            if fi==0
                hlv_nf3=2087580;
                ro_nf3=ro_lv3;
            else
                C03=(Anp*(fi.^5))+(Bnp*(fi.^4))+(Cnp*(fi.^3))+(Dnp*(fi.^2))+(Enp*fi)+Fnp;
                C13=(alphanp*(fi.^5))+(bettanp*(fi.^4))+(gammanp*(fi.^3))+(sigmanp*(fi.^2))+(epsilonnp*fi)+omeganp;
                hlv_nf3=C03*(deltaPsat3.^C13);
            end
            A43=(ro_nf3*hlv_nf3*V_lv3)/Length3;% From eq 21
            Nu_c3=0.52*(((g*betta_c*(D_a_ex.^3))/(v_c*alpha_c)).^0.2)*(Tstar_a_ex3-Tstar_g_in3).^0.2;% eq 34
            h_c3=(Nu_c3*k_c)/D;% eq 33
            A53=(h_c3*Length3)/N;% From eq 32
            Nu_e3=0.15*(((g*betta_e*(w_c.^3))/(v_e*alpha_e)).^(1/3))*((Tstar_g_ex3-T_amb).^(1/3));% eq 27
            h_e3=(Nu_e3*k_e)/w_c;% eq 26
            A63=h_e3*w;% From eq 25
            b13=((I_13)*w);% From eq 1
            b23=0;% From eq 2
            b33=((I_23)*w);% From eq 3
            b43=0;% From eq 4
            b53=-373.15*A33;% From eq 5
            b63=(-A43)*x_in3;% From eq 6
            b73=0;% From eq 7
            b83=(-A63)*T_amb;% From eq 8
            b93=0;% From eq 9
            b103=0;% From eq 10
            b113=0;% From eq 11
            b123=0;% From eq 12
            A3=[1,0,0,0,0,0,0,0,0,0,0,0% eq 1
                 0,1,0,0,0,0,0,-A13,A13,0,0,0% eq 2
                 0,0,1,0,0,0,0,0,0,0,0,0% eq 3
                 0,0,0,1,0,0,0,0,0,-A23,A23,0% eq 4
                 0,0,0,0,1,0,0,0,0,0,-A33,0% eq 5
                 0,0,0,0,1,0,0,0,0,0,0,-A43% eq 6
                 0,0,0,0,0,1,0,0,A53,-A53,0,0% eq 7
                 0,0,0,0,0,0,1,-A63,0,0,0,0% eq 8
                 0,1,0,0,0,1,0,0,0,0,0,0% eq 9
                 1,0,0,0,0,N,-1,0,0,0,0,0% eq 10
                 0,0,1,0,-1,-1,0,0,0,0,0,0% eq 11
                 0,0,0,1,-1,0,0,0,0,0,0,0];% eq 12
                 B3=[b13;b23;b33;b43;b53;b63;b73;b83;b93;b103;b113;b123];
                 %X=[qdat_g_srad;qdat_g_cond;qdat_a_srad;qdat_a_cond;qdat_af_conv;qdat_ag_conv;q_ge_conv;T_g_ex;T_g_in;T_a_ex;T_a_in;x_out];
                 X3=A3\B3;
                 T_a_ex3=X3(10,:);
                 T_g_in3=X3(9,:);
                 T_g_ex3=X3(8,:);
                 x_out3=X3(12,:);
                 Delta13=(T_g_ex3-Tstar_g_ex3);
                 Delta23=(T_a_ex3-Tstar_a_ex3);
                 Delta33=(T_g_in3-Tstar_g_in3);
                 Delta43=(x_out3-xstar_out3);
                 Tstar_g_ex3=Tstar_g_ex3+0.5*(Delta13);
                 Tstar_a_ex3=Tstar_a_ex3+0.5*(Delta23);
                 Tstar_g_in3=Tstar_g_in3+0.5*(Delta33);
                 xstar_out3=xstar_out3+0.5*(Delta43);
        end
        x_out3(i)=X3(12,:);
        x_in3=x_out3(i);
        Quality31(i)=X3(12,:);
        vq3=Quality31';%vq : vapore quality
        TEMP32(i)=X3(11,:);
        TIS3=TEMP32';%TIS : Temperature of Internal Surface
        Flux3(i)=X3(5,:);
        HFaf3=Flux3';%HFaf : Convection Heat Flux from Absorber to Fluid
        Nuselt_c3(i)=Nu_c3;
        Nuselt_c3=Nuselt_c3';
        Nuselt_e3(i)=Nu_e3;
        Nuselt_e3=Nuselt_e3';
        HT3(i)=h_tp3;
        HT33=HT3';
    end
    disp 'Distribution of Internal Surface Temperature at Two-Phase Flow Section is:';
    disp (TIS3);
    a31=linspace(0,l3,n3+1)';%because i have not first point so start from L/n1; n1 = numbber of the distances
    y31=[T2;TIS3];
    c31=polyfit(a31,y31,n3);%because i have not first point so degree of polynomial is n1-1
    disp 'The Internal Surface Temperature at Start of Two-Phase Flow Section:'
    T_a_in03=c31(1,n3+1)
    disp 'Distribution of Vapore Quality at Two-Phase Flwo Section is:';
    disp (vq3);
    a32=linspace(0,l3,n3+1)';
    y32=[0;vq3];
    c32=polyfit(a32,y32,n3);
    disp 'Distribution of Convection Heat Flux from Absorber to Fluid at Two-Phase Flow Section is:';
    disp (HFaf3);
else
    %% when value of wind velosity > 0.1
    syms qdat_g_srad1 qdat_g_cond1 qdat_ag_conv1 qdat_ag_trad1 qdat_af_conv1 qdat_a_srad1 qdat_a_cond1 q_ge_conv1 T_g_ex1 T_g_in1 T_a_ex1 T_a_in1 T_out1
    %fluid parameters for moderate temperature of 335 (K)
    C_p_f1=input('value of specific heat for fluid at preheating section in the 335 (K):');%specific heat at constant pressure=C_p
    mu_f1=input('value of fluid viscosity at preheating section in the 335 (K):');%viscosity=mu
    k_f1=input('value of thermal conductivity for fluid at preheating section in the 335 (K):');%thermal conductivity of fluid=k_f
    ro1=input('value of fluid density at preheating section in the 335 (K):');%density=ro
    V_f1=0.07;%fluid velosity=V
    L=input('value of pipe length:');%length of absorber at Pre Heating Section
    n1=input('please input number of elements for absorber:');
    Length1=L/n1;% length of each absorber element at Pre Heating Section (meter)
    qdatsun_g1=I_sun*betta*w*Cr_g;
    qdatsun_a1=I_sun*betta*w*Cr_a;
    F_rad=1;% Based on Fig. 13.3 (p.865) in the Incropera
    Pr_f1=(C_p_f1*mu_f1)/k_f1;
    Re_f1=(ro1*V_f1*D_a_in)/(mu_f1);
    Nu_f1=(0.023)*(Re_f1.^(4/5))*(Pr_f1.^0.4);
    h_f1=(Nu_f1*k_f1)/D_a_in;
    A_i1=pi*((D_a_in/2).^2);
    mdat1=A_i1*ro1*V_f1;
    T_out1=zeros(n1,1);
    T_a_ex1=zeros(n1,1);
    T_a_in1=zeros(n1,1);
    T_g_ex1=zeros(n1,1);
    T_g_in1=zeros(n1,1);
    %input temperature=T_i
    T_in1=input('Value of Fluid Temperature at Start Process:');%ok
    for i=1:n1 % quantity of elements
        i=i;
        % initial temperatures= Tstar_a_ex, Tstar_g_ex, Tstar_g_in, T_g_ex
        Tstar_a_ex1=500;
        Tstar_g_in1=407.485;
        while (abs(T_a_ex1-Tstar_a_ex1))>0.01 & (abs(T_g_in1-Tstar_g_in1))>0.01
            A11=(2*pi*k_g)/(log(D_g_ex/D_g_in));%ok
            A21=(2*pi*k_a)/(log(D_a_ex/D_a_in));%ok
            A31=h_f1*pi*D_a_in;%ok
            A41=(mdat1*C_p_f1)/Length1;%ok
            Ra_L=((g*betta_c)/(v_c*alpha_c))*(Tstar_a_ex1-Tstar_g_in1)*(D_a_ex.^3);
            Ra_c1=(((log(D_g_in*D_a_ex)).^4)/((L.^3)*(((D_a_ex.^(-0.6))+(D_g_in.^(-0.6))).^5)))*Ra_L;
            if (Ra_c1>10.^2) & (Ra_c1<10.^7)
               Pr_c1=(C_p_c*mu_c)/k_c;
               k_eff1=k_c*0.386*((Pr_c1/(0.861+Pr_c1)).^0.25)*(Ra_c1.^0.25);
            elseif (Ra_c1<10.^2)
                k_eff1=k_c;
            end
            h_c1=(2*k_eff1)/(D_a_ex*(log(D_g_in/D_a_ex)));
            A51=h_c1*pi*D_a_ex;%ok
            A61=F_rad*sigma*pi*D_a_ex*(Tstar_a_ex1+Tstar_g_in1)*((Tstar_a_ex1.^2)+(Tstar_g_in1.^2));%ok
            Re_e1=(ro_e*v_air*D_g_ex)/(mu_e);
            Pr_e1=(C_p_e*mu_e)/k_e;
            Nu_e1=0.3+(((0.62*(Re_e1.^0.5)*(Pr_e1.^(1/3)))/((1+((0.4/Pr_e1).^(2/3))).^0.25))*((1+((Re_e1/282000).^(5/8))).^0.8));
            h_e1=(Nu_e1*k_e)/D_g_ex;
            A71=h_e1*pi*D_g_ex;%ok
            b11=((qdatsun_g1)*alpha_g);%ok
            b21=0;%ok
            b31=((qdatsun_a1)*gamma*alpha_a);%ok
            b41=0;%ok
            b51=-(1/2)*A31*T_in1;%ok
            b61=(-A41)*T_in1;%ok
            b71=0;%ok
            b81=0;%ok
            b91=(-A71)*T_amb;%ok
            b101=0;%ok
            b111=0;%ok
            b121=0;%ok
            b131=0;%ok
            A1=[1,0,0,0,0,0,0,0,0,0,0,0,0%ok 1
                 0,1,0,0,0,0,0,0,-A11,A11,0,0,0%ok 2
                 0,0,1,0,0,0,0,0,0,0,0,0,0%ok 3
                 0,0,0,1,0,0,0,0,0,0,-A21,A21,0%ok 4
                 0,0,0,0,1,0,0,0,0,0,0,-A31,(A31)/2%ok 5
                 0,0,0,0,1,0,0,0,0,0,0,0,-A41%ok 6
                 0,0,0,0,0,1,0,0,0,-A51,A51,0,0%ok 7
                 0,0,0,0,0,0,1,0,0,-A61,A61,0,0%ok 8
                 0,0,0,0,0,0,0,1,-A71,0,0,0,0%ok 9
                 1,0,0,0,0,1,1,-1,0,0,0,0,0%ok 10
                 0,0,1,0,-1,-1,-1,0,0,0,0,0,0%ok 11
                 1,-1,0,0,0,0,0,0,0,0,0,0,0%ok 12;(12.01.2017)
                 0,0,0,1,-1,0,0,0,0,0,0,0,0];%ok 13
                 B1=[b11;b21;b31;b41;b51;b61;b71;b81;b91;b101;b111;b121;b131];
                 %X=[qdat_g_srad;qdat_g_cond;qdat_a_srad;qdat_a_cond;qdat_af_conv;qdat_ag_conv;qdat_ag_trad;q_ge_conv;T_g_ex;T_g_in;T_a_ex;T_a_in;T_out];
                 X1=A1\B1;
                 T_a_ex1=X1(11,:);
                 T_g_in1=X1(10,:);
                 T_g_ex1=X1(9,:);
                 F11=(T_a_ex1-Tstar_a_ex1);
                 F21=(T_g_in1-Tstar_g_in1);
                 Tstar_a_ex1=Tstar_a_ex1+0.5*(F11);
                 Tstar_g_in1=Tstar_g_in1+0.5*(F21);
        end       
        T_out1(i)=X1(13,:);
        T_in1=T_out1(i);
        TEMP11(i)=X1(13,:);
        TF1=TEMP11';%TF : Temperature of Fluid
        TEMP12(i)=X1(12,:);
        TIS1=TEMP12';%TIS : Temperature of Internal Surface
        Flux1(i)=X1(5,:);
        HFaf1=Flux1';%HFaf : Convection Heat Flux from Absorber to Fluid
    end
    disp 'Distribution of Internal Surface Temperature at Pre Heating Section is:';
    disp (TIS1);
    a11=linspace(L/n1,L,n1)';%because i have not first point so start from L/n1; n1 = numbber of the distances
    y11=TIS1;
    c11=polyfit(a11,y11,n1-1);%because i have not first point so degree of polynomial is n1-1
    disp 'The Internal Surface Temperature at Start of Preheating Section:'
    T_a_in01=c11(1,n1);%Temperature of internal surface at start of preheating section
    sizec11=size(c11);
    C1=sizec11(1,2);%for determine number of members
    DELTA1=c11(1,C1)-373.15;%for determine equation of TIS1 at 373.15
    c11(:,C1)=[];%for determine equation of TIS1 at 373.15
    d1=[c11,DELTA1];%equation of TIS1 at 373.15
    q1=roots(d1);
    z1=size(q1);
    disp 'The Length of Preheating Section:'
    j=1:1:z1(1,1);
    x1=zeros(j,1);
    for j=1:1:z1(1,1)
        M1=q1(j);
        N1=imag(q1(j));
        if M1>0 & M1<L & N1==0
            x1(j)=M1;
        end
    end
    x11=nonzeros(x1);
    l1=min(x11)
    disp 'Distribution of Fluid Temperature at Pre Heating Section is:';
    disp (TF1);
    a12=linspace(0,L,n1+1)';
    y12=[293.15;TF1];
    c12=polyfit(a12,y12,n1);%change size... with n1
    disp 'The Fluid Temperature at End of Preheating Section:'%
    T1=polyval(c12,l1)%
    disp 'The Different Temperature of Internal Surface & Fluid at End of Preheating Section:'
    deltaT1=373.15-T1;
    disp 'Distribution of Convection Heat Flux from Absorber to Fluid at Preheating Section is:';
    disp (HFaf1);
    a13=linspace(L/n1,L,n1)';
    y13=HFaf1;
    c13=polyfit(a13,y13,n1-1);
    disp 'The Heat Flux at Start of Preheating Section:';
    HF10=c13(1,n1);
    disp 'The Heat Flux at End of Preheating Section:';
    HF11=polyval(c13,l1);
    %% Start of Calculation for Subcooled Section
    %fluid parameters for moderate temperature of 355 (K)
    C_p_f2=input('value of specific heat for fluid at subcooled section:');%specific heat at constant pressure=C_p
    mu_f2=input('value of fluid viscosity at subcooled section:');%viscosity=mu
    k_f2=input('value of thermal conductivity for fluid at subcooled section:');%thermal conductivity of fluid=k_f
    ro_Length2=input('value of liquid density at subcooled section:');% water density
    ro_v2=input('value of vapor density at subcooled section:');% vapor density
    hlv2=input('value of enthalpy at subcooled section:');%Enthalpy of Fluid at Subcooled Section (j/kg)
    Psat=input('value of saturated pressure at subcooled section:');%Saturate Pressure of Fluid at Subcooled Section
    V_f2=0.07;%fluid velosity=V
    L1=L-l1;%length of absorber at Subcoled Section
    n2=input('please input number of elements for absorber:');
    Length2=L1/n2;%length of each absorber element at Subcoled Section (meter)
    qdatsun_g2=I_sun*betta*w*Cr_g;
    qdatsun_a2=I_sun*betta*w*Cr_a;
    F_rad=1;% Based on Fig. 13.3 (p.865) in the Incropera
    A_i2=pi*((D_a_in/2).^2);
    mdat2=A_i2*ro_Length2*V_f2;
    T_out2=zeros(n2,1);
    T_a_ex2=zeros(n2,1);
    T_a_in2=zeros(n2,1);
    T_g_ex2=zeros(n2,1);
    T_g_in2=zeros(n2,1);
    %input temperature=T_i
    T_in2=T1;% (Temperature for end of Pre Heating Section)
    for i=1:n2 % quantity of elements
        i=i;
        % initial temperatures= Tstar_a_ex, Tstar_g_ex, Tstar_g_in, T_g_ex
        Tstar_a_ex2=880;
        Tstar_a_in2=800;
        Tstar_g_in2=507.485;
        Tstar_out2=473.15;
        while (abs(T_a_ex2-Tstar_a_ex2))>0.01 & (abs(T_g_in2-Tstar_g_in2))>0.01 & (abs(T_a_in2-Tstar_a_in2))>0.01 & (abs(T_out2-Tstar_out2))>0.01
            A12=(2*pi*k_g)/(log(D_g_ex/D_g_in));%ok
            A22=(2*pi*k_a)/(log(D_a_ex/D_a_in));%ok
            G=0.00122*(((k_f2.^0.79)*(C_p_f2.^0.45)*(ro_Length2.^0.49))/((sigma.^0.5)*(mu_f2.^0.29)*(hlv2.^0.24)*(ro_v2.^0.24)))*(Psat.^0.75);
            h_nb=G*((Tstar_a_in2-Tstar_out2).^0.99);%h_nb : Heat Transfer of Nucleate Boiling, rev03 : T,sat = T,out
            A32=h_nb*pi*D_a_in;
            A42=(mdat2*C_p_f2)/Length2;%ok
            Ra_L1=((g*betta_c)/(v_c*alpha_c))*(Tstar_a_ex2-Tstar_g_in2)*(D_a_ex.^3);
            Ra_c2=(((log(D_g_in*D_a_ex)).^4)/((L1.^3)*(((D_a_ex.^(-0.6))+(D_g_in.^(-0.6))).^5)))*Ra_L1;
            if (Ra_c2>10.^2) & (Ra_c2<10.^7)
               Pr_c2=(C_p_c*mu_c)/k_c;
               k_eff2=k_c*0.386*((Pr_c2/(0.861+Pr_c2)).^0.25)*(Ra_c2.^0.25);
            elseif (Ra_c2<10.^2)
                k_eff2=k_c;
            end
            h_c2=(2*k_eff2)/(D_a_ex*(log(D_g_in/D_a_ex)));
            A52=h_c2*pi*D_a_ex;%ok
            A62=F_rad*sigma*pi*D_a_ex*(Tstar_a_ex2+Tstar_g_in2)*((Tstar_a_ex2.^2)+(Tstar_g_in2.^2));%ok
            Re_e2=(ro_e*v_air*D_g_ex)/(mu_e);
            Pr_e2=(C_p_e*mu_e)/k_e;
            Nu_e2=0.3+(((0.62*(Re_e2.^0.5)*(Pr_e2.^(1/3)))/((1+((0.4/Pr_e2).^(2/3))).^0.25))*((1+((Re_e2/282000).^(5/8))).^0.8));
            h_e2=(Nu_e2*k_e)/D_g_ex;
            A72=h_e2*pi*D_g_ex;%ok
            b12=((qdatsun_g2)*alpha_g);%ok
            b22=0;%ok
            b32=((qdatsun_a2)*gamma*alpha_a);%ok
            b42=0;%ok
            b52=0;%T,sat = T,out So : b5=0
            b62=(-A42)*T_in2;%ok
            b72=0;%ok
            b82=0;%ok
            b92=(-A72)*T_amb;%ok
            b102=0;%ok
            b112=0;%ok
            b122=0;%ok
            b132=0;%ok
            A2=[1,0,0,0,0,0,0,0,0,0,0,0,0%ok 1
                 0,1,0,0,0,0,0,0,-A12,A12,0,0,0%ok 2
                 0,0,1,0,0,0,0,0,0,0,0,0,0%ok 3
                 0,0,0,1,0,0,0,0,0,0,-A22,A22,0%ok 4
                 0,0,0,0,1,0,0,0,0,0,0,-A32,A32%ok 5
                 0,0,0,0,1,0,0,0,0,0,0,0,-A42%ok 6
                 0,0,0,0,0,1,0,0,0,-A52,A52,0,0%ok 7
                 0,0,0,0,0,0,1,0,0,-A62,A62,0,0%ok 8
                 0,0,0,0,0,0,0,1,-A72,0,0,0,0%ok 9
                 1,0,0,0,0,1,1,-1,0,0,0,0,0%ok 10
                 0,0,1,0,-1,-1,-1,0,0,0,0,0,0%ok 11
                 1,-1,0,0,0,0,0,0,0,0,0,0,0%ok 12
                 0,0,0,1,-1,0,0,0,0,0,0,0,0];%ok 13
                 B2=[b12;b22;b32;b42;b52;b62;b72;b82;b92;b102;b112;b122;b132];
                 %X=[qdat_g_srad;qdat_g_cond;qdat_a_srad;qdat_a_cond;qdat_af_conv;qdat_ag_conv;qdat_ag_trad;q_ge_conv;T_g_ex;T_g_in;T_a_ex;T_a_in;T_out];
                 X2=A2\B2;
                 T_a_ex2=X2(11,:);
                 T_a_in2=X2(12,:);
                 T_g_in2=X2(10,:);
                 T_out2=X2(13,:);
                 F12=(T_a_ex2-Tstar_a_ex2);
                 F22=(T_g_in2-Tstar_g_in2);
                 F42=(T_a_in2-Tstar_a_in2);
                 F52=(T_out2-Tstar_out2);
                 Tstar_a_ex2=Tstar_a_ex2+0.5*(F12);
                 Tstar_g_in2=Tstar_g_in2+0.5*(F22);
                 Tstar_a_in2=Tstar_a_in2+0.5*(F42);
                 Tstar_out2=Tstar_out2+0.5*(F52);
        end       
        T_out2(i)=X2(13,:);
        T_in2=T_out2(i);
        TEMP21(i)=X2(13,:);
        TF2=TEMP21';%TF : Temperature of Fluid
        TEMP22(i)=X2(12,:);
        TIS2=TEMP22';%TIS : Temperature of Internal Surface
        Flux2(i)=X2(5,:);
        HFaf2=Flux2';%HFaf : Convection Heat Flux from Absorber to Fluid
    end
    disp 'Distribution of Fluid Temperature at Subcooled Section is:';
    disp (TF2);
    a21=linspace(0,L1,n2+1)';
    y21=[T1;TF2];
    c21=polyfit(a21,y21,n2);
    disp 'The Fluid Temperature at Start of Subcooled Section:'
    T_in02=c21(1,n2+1);
    sizec21=size(c21);
    C2=sizec21(1,2);%for determine number of members
    DELTA2=c21(1,C2)-373.15;%for determine equation of TF2 at 373.15
    c21(:,C2)=[];%for determine equation of TF2 at 373.15
    d2=[c21,DELTA2];%equation of TF2 at 373.15
    q2=roots(d2);
    z2=size(q2);
    disp 'The Length of Subcooled Section:'
    j=1:1:z2(1,1);
    x2=zeros(j,1);
    for j=1:1:z2(1,1)
        M2=q2(j);
        N2=imag(q2(j));
        if M2>0 & M2<L1 & N2==0
            x2(j)=M2;
        end
    end
    x22=nonzeros(x2);
    l2=min(x22)
    disp 'Distribution of Internal Surface Temperature at Subcooled Section is:';
    disp (TIS2);
    a22=linspace(0,L1,n2+1)';
    y22=[373.15;TIS2];
    size22=size(y22);
    c22=polyfit(a22,y22,n2);%change size... with n2
    disp 'The Internal Surface Temperature at Start of Subcooled Section:'
    T_a_in02=c22(1,n2+1);
    disp 'The The Internal Surface Temperature at End of Subcooled Section:'
    T2=polyval(c22,l2)
    disp 'The Different Temperature of Internal Surface & Fluid at End of Subcooled Section:'
    deltaT2=T2-373.15
    disp 'Distribution of Convection Heat Flux from Absorber to Fluid at Subcooled Section is:';
    disp (HFaf2);
    a23=linspace(0,L1,n2+1)';
    y23=[HF11;HFaf2];
    size23=size(y23);
    c23=polyfit(a23,y23,n2);
    disp 'The Heat Flux at Start of Preheating Section:';
    HF20=c23(1,n2+1);
    disp 'The Heat Flux at End of Preheating Section:';
    HF21=polyval(c23,l2);
    %% Start of Calculation for Two-Phase Flow Section
    C_p_f3=4217;%specific heat at constant pressure=C_p
    mu_f3=279e-006;%viscosity=mu
    k_f3=680e-003;%thermal conductivity of fluid=k_f
    ro_l3=957.85;% water density
    ro_v3=0.595;% vapor density
    ro_lv3=0.598;%vapor-liquid density
    V_lv3=1.67185;
    h_lv3=2087580;%Enthalpy of Fluid at Subcooled Section (j/kg)
    V_f3=0.07;%fluid velosity=V
    Pr_f3=(C_p_f3*mu_f3)/k_f3;
    Re_f3=(ro_l3*V_f3*D_a_in)/(mu_f3);
    Nu_f3=(0.023)*(Re_f3.^(4/5))*(Pr_f3.^0.4);
    h_l3=(Nu_f3*k_f3)/D_a_in;
    pr3=0.031;
    l3=L-(l1+l2);%length of absorber;............ok
    n3=input('please input number of elements for absorber:');.........ok
    Length3=l3/n3;%length of each absorber element (meter).......ok
    qdatsun_g3=I_sun*betta*w*Cr_g;
    qdatsun_a3=I_sun*betta*w*Cr_a;
    F_rad=1;% Based on Fig. 13.3 (p.865) in the Incropera
    x_out3=zeros(n3,1);
    T_a_ex3=zeros(n3,1);
    T_g_ex3=zeros(n3,1);
    T_g_in3=zeros(n3,1);
    %input temperature=T_i
    x_in3=0;% (vapor quality for end of Subcooled Section)
    for i=1:n3 % quantity of elements
        i=i;
        % initial temperatures= Tstar_a_ex, Tstar_g_ex, Tstar_g_in, T_g_ex
        Tstar_a_ex3=880;
        Tstar_g_in3=507.485;
        xstar_out3=0.9;%delete of T_out and add the x=vapore quality
        while (abs(T_a_ex3-Tstar_a_ex3))>0.01 & (abs(T_g_in3-Tstar_g_in3))>0.01 & (abs(x_out3-xstar_out3))>0.01
            A13=(2*pi*k_g)/(log(D_g_ex/D_g_in));%ok
            A23=(2*pi*k_a)/(log(D_a_ex/D_a_in));%ok
            h_tp3=((h_l3*(pr3.^0.38)*((1-xstar_out3).^0.8))+(3.8*(xstar_out3.^0.76)*((1-xstar_out3).^0.04)))*(pr3.^(-0.38));
            A33=h_tp3*pi*D_a_in;
            A43=(ro_lv3*h_lv3*V_lv3)/Length3;%ok
            Ra_l3=((g*betta_c)/(v_c*alpha_c))*(Tstar_a_ex3-Tstar_g_in3)*(D_a_ex.^3);
            Ra_c3=(((log(D_g_in*D_a_ex)).^4)/((l3.^3)*(((D_a_ex.^(-0.6))+(D_g_in.^(-0.6))).^5)))*Ra_l3;
            if (Ra_c3>10.^2) & (Ra_c3<10.^7)
               Pr_c3=(C_p_c*mu_c)/k_c;
               k_eff3=k_c*0.386*((Pr_c3/(0.861+Pr_c3)).^0.25)*(Ra_c3.^0.25);
            elseif (Ra_c3<10.^2)
                k_eff3=k_c;
            end
            h_c3=(2*k_eff3)/(D_a_ex*(log(D_g_in/D_a_ex)));
            A53=h_c3*pi*D_a_ex;%ok
            A63=F_rad*sigma*pi*D_a_ex*(Tstar_a_ex3+Tstar_g_in3)*((Tstar_a_ex3.^2)+(Tstar_g_in3.^2));%ok
            Re_e3=(ro_e*v_air*D_g_ex)/(mu_e);
            Pr_e3=(C_p_e*mu_e)/k_e;
            Nu_e3=0.3+(((0.62*(Re_e3.^0.5)*(Pr_e3.^(1/3)))/((1+((0.4/Pr_e3).^(2/3))).^0.25))*((1+((Re_e3/282000).^(5/8))).^0.8));
            h_e3=(Nu_e3*k_e)/D_g_ex;
            A73=h_e3*pi*D_g_ex;%ok
            b13=((qdatsun_g3)*alpha_g);%ok
            b23=0;%ok
            b33=((qdatsun_a3)*gamma*alpha_a);%ok
            b43=0;%ok
            b53=-373.15*A33;%ok
            b63=(-A43)*x_in3;%ok
            b73=0;%ok
            b83=0;%ok
            b93=(-A73)*T_amb;%ok
            b103=0;%ok
            b113=0;%ok
            b123=0;%ok
            b133=0;%ok
            A3=[1,0,0,0,0,0,0,0,0,0,0,0,0%ok 1
                 0,1,0,0,0,0,0,0,-A13,A13,0,0,0%ok 2
                 0,0,1,0,0,0,0,0,0,0,0,0,0%ok 3
                 0,0,0,1,0,0,0,0,0,0,-A23,A23,0%ok 4
                 0,0,0,0,1,0,0,0,0,0,0,-A33,0%ok 5
                 0,0,0,0,1,0,0,0,0,0,0,0,-A43%ok 6
                 0,0,0,0,0,1,0,0,0,-A53,A53,0,0%ok 7
                 0,0,0,0,0,0,1,0,0,-A63,A63,0,0%ok 8
                 0,0,0,0,0,0,0,1,-A73,0,0,0,0%ok 9
                 1,0,0,0,0,1,1,-1,0,0,0,0,0%ok 10
                 0,0,1,0,-1,-1,-1,0,0,0,0,0,0%ok 11
                 1,-1,0,0,0,0,0,0,0,0,0,0,0%ok 12
                 0,0,0,1,-1,0,0,0,0,0,0,0,0];%ok 13
                 B3=[b13;b23;b33;b43;b53;b63;b73;b83;b93;b103;b113;b123;b133];
                 %X=[qdat_g_srad;qdat_g_cond;qdat_a_srad;qdat_a_cond;qdat_af_conv;qdat_ag_conv;qdat_ag_trad;q_ge_conv;T_g_ex;T_g_in;T_a_ex;T_a_in;x_out];
                 X3=A3\B3;
                 T_a_ex3=X3(11,:);
                 T_g_in3=X3(10,:);
                 x_out3=X3(13,:);
                 F13=(T_a_ex3-Tstar_a_ex3);
                 F23=(T_g_in3-Tstar_g_in3);
                 F43=(x_out3-xstar_out3);
                 Tstar_a_ex3=Tstar_a_ex3+0.5*(F13);
                 Tstar_g_in3=Tstar_g_in3+0.5*(F23);
                 xstar_out3=xstar_out3+0.5*(F43);
        end
        x_out3(i)=X3(13,:);
        x_in3=x_out3(i);
        TEMP31(i)=X3(13,:);
        vq3=TEMP31';%vq : vapore quality
        TEMP32(i)=X3(12,:);
        TIS3=TEMP32';%TIS : Temperature of Internal Surface
        Flux3(i)=X3(5,:);
        HFaf3=Flux3';%HFaf : Convection Heat Flux from Absorber to Fluid
    end
    disp 'Distribution of Internal Surface Temperature at Two-Phase Flow Section is:';
    disp (TIS3);
    a31=linspace(0,l3,n3+1)';%because i have not first point so start from L/n1; n1 = numbber of the distances
    y31=[T2;TIS3];
    c31=polyfit(a31,y31,n3);%because i have not first point so degree of polynomial is n1-1
    disp 'The Internal Surface Temperature at Start of Two-Phase Flow Section:'
    T_a_in03=c31(1,n3+1);
    sizec33=size(c31);
    C3=sizec33(1,2);%for determine number of members
    deltaT=input('Value of Different Temperature at End of Two-Phase Flow Section:');
    T_a_in3=373.15+deltaT;
    B3=c31(1,C3)-T_a_in3;%for determine equation of TIS1 at 433.15
    c31(:,C3)=[];%for determine equation of TIS1 at 373.15
    d3=[c31,B3];%equation of TIS1 at 373.15
    q3=roots(d3);
    z3=size(q3);
    disp 'The Length of Preheating Section:'
    j=1:1:z3(1,1);
    x3=zeros(j,1);
    for j=1:1:z3(1,1)
        M3=q3(j);
        N3=imag(q3(j));
        if M3>0 & M3<l3 & N3==0
            x3(j)=M3;
        end
    end
    x33=max(x3)
    disp 'Distribution of Vapore Quality at Two-Phase Flwo Section is:';
    disp (vq3);
    a32=linspace(0,l3,n3+1)';
    y32=[0;vq3];
    c32=polyfit(a32,y32,n3);%change size... with n2
    disp 'The The Internal Surface Temperature at End of Subcooled Section:'
    x_out=polyval(c32,x33);
    disp 'Distribution of Convection Heat Flux from Absorber to Fluid at Two-Phase Flow Section is:';
    disp (HFaf3);
end