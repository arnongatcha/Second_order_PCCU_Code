%%                               NGATCHA NDENGNA ARNO ROLAND 
%%  NATIONAL ADVANCED SCHOOL OF MARITIME  AND OCEAN SCIENCE AND TECHNOLOGY (NASMOST) OF EBOLOWA UNIVERSITY
%%                              BOX P.O. 292 KRIBI, CAMEROON.
%% 
%% CODE  FOR SEDIMENT TRANSPORT WITH TURBULENCE. THE MODEL COMPUTED IS AVAILABLE IN THE PAPER
%% ENTITLED " A high order Path-Conservative Central-Upwind Arbitrary DERivative (PCCU-ADER)
%% method for a generalized high order sediment transport model" 
%% THE CODE IS WRITTEN IN 1D VERSION AND IS USED FOR ALL THE 1D TESTS PRESENT IN THIS PAPER. 
%%  



close all
clear ; clf

set(0,'defaultTextInterpreter','latex')
set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',16,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',16,...
'DefaultLineLineWidth',3.0,...
'DefaultAxesBox','on',...
'defaultAxesLineWidth',1.0,...
'DefaultFigureColor','w',...
'DefaultLineMarkerSize',1.75)


% gravitation g
g    = 9.81 ;
pore = 0.4 ; 

% nx = number of cells in x-direction
nx   = 100 ;
 
% initialisation 
x     =  zeros(1,nx) ;   
hue   =  zeros(1,nx) ; % Water discharge in x-direction
hve   =  zeros(1,nx) ; % Water discharge in y-direction
hce   =  zeros(1,nx) ; 
he    =  zeros(1,nx) ;  % Water depth
ce    =  zeros(1,nx) ;  % Averaged Sediment concentration
ue    =  zeros(1,nx) ;  %   Averaged Water velocity in x-direction
ve    =  zeros(1,nx) ;  % Averaged Water velocity in y-direction
P11e  =  zeros(1,nx) ;  % Shear stress component (in x-direction)
P12e  =  zeros(1,nx)  ;  % Shear stress component (in shear plane)
P22e  =   zeros(1,nx)  ; % Shear stress component (in y-direction)
be    =  zeros(1,nx) ;   % Bottom evolution
Ee    =  zeros(1,nx) ;   
qe    =  zeros(1,nx) ;   % Water discharge
hK11e =  zeros(1,nx) ;  
hK12e =  zeros(1,nx) ; 
hK22e =  zeros(1,nx) ;
K11e  =   zeros(1,nx)  ; % Kinetic energy in x-direction
K12e  =   zeros(1,nx)  ; % Kinetic energy  in shear plane
K22e  =   zeros(1,nx)  ; % Kinetic energy in y-direction
  
SSW   = 1 ;            % SSW = Shear Shallow Water = Shallow Water with turbulence
SW    = 0 ;            % SW = Shallow Water without turbulence
 
TEST_1 = 1;          % Perform dambreak with sediment transport see Ngatcha et al (2024)
TEST_2 = 0 ;        % Perform dambreak without sediment transport 
% Other tests can be easily implemented. Please add  TEST_3, TEST_4, etc. 

Model = SSW; 

        switch Model

            case SSW;
                a = 1 ;
               
            case SW
                a=0 ;
          
        end 
    
        
 Second_order = 1 ;
 First_order  = 0 ;
  
 
for i=1:nx
   
    % Dambreak with sediment transport following BENKHALDOUN ET AL (2012)
    
  %  Tests = dambreak ;  
   % Cla = 1 ;   
   
    % dambreak without sediment transport (active please if necessary)
    
     Tests = TEST_2 ;    
     Cla = 0 ;
   
 
       switch Tests
         
         case TEST_1
             
 % Domain length            
 lg = 50000;  

% dx = step size space
dx   = lg/nx ; 

cfl      =  0.1        ;   % courant-Fredirich-Lewy
tstop    = 300        ;
tscreen  = 0           ;
dtscreen = tstop/10.0  ;    
 x(i) = (i-1)*dx  ; 
    if x(i) <= lg/2   
     
 
     
     %% DAM BREAK TEST FAYSAL BENKHALDOUM 
     ue(i)    = 0.;  
     ve(i)    = 0.    ;
     he(i)    = 40 ; 
  
     be(i)    = 0.0  ;
     ce(i)    = 0.001 ; 
     hce(i)   = he(i)*ce(i);
     
     P11e(i)  =  a*1.00E-04;
     P12e(i)  =  0.00  ;  
     P22e(i)  =  0*a*1.00E-04 ; 
     

     hK11e(i) =  0.5*he(i)*P11e(i) + 0.5*he(i)*ue(i)*ue(i) ;    %  
     hK12e(i) =  0.5*he(i)*P12e(i) + 0.5*he(i)*ue(i)*ve(i) ;    %
     hK22e(i) =  0.5*he(i)*P22e(i) + 0.5*he(i)*ve(i)*ve(i) ;    %  
          
  else
     ue(i)    =  0;  
     ve(i)    =  0;
     he(i)    =  2; 
   
     be(i)    = 0.0  ;
     ce(i)    = 0.001  ; 
     hce(i)   = he(i)*ce(i);
     
     
     P11e(i)  =  a*1.00E-04 ;  
     P12e(i)  =  0.00  ;  
     P22e(i)  =  0*a*1.00E-04 ; 
     
     
     hK11e(i) =  0.5*he(i)*P11e(i) + 0.5*he(i)*ue(i)*ue(i) ;    %  
     hK12e(i) =  0.5*he(i)*P12e(i) + 0.5*he(i)*ue(i)*ve(i) ;    %
     hK22e(i) =  0.5*he(i)*P22e(i) + 0.5*he(i)*ve(i)*ve(i) ;    %  
     
     end  
  
       case TEST_2  % Modified dam break 
           
% The domain length
lg = 1 ;  
% The step space
dx        =  lg/nx       ; 

cfl       =  0.1         ;   % courant-Fredirich-Lewy
tstop     =  0.5         ; % final time
tscreen   =  0           ;
dtscreen  = tstop/10.0   ;
     x(i) = (i-1)*dx     ; % Points of discretization 
     
   
       if x(i) <= lg/2  
             Classical = 0 ; 
     ue(i) =  0.1;  
     ve(i)  = 0.2    ;
     he(i)  = 0.01 ; 
     hue(i) = he(i)*ue(i) ; 
     hve(i) = he(i)*ve(i) ;
     be(i)  = 0.  ;
     ce(i)  = 0.0 ; 
     
     P11e(i) = 4E-02 ;  
     P12e(i) = 1E-08 ;  
     P22e(i) = 4E-02 ; 
     
     hK11e(i) =   0.5*he(i)*P11e(i) + 0.5*he(i)*ue(i)*ue(i) ;     %  
     hK12e(i) =   0.5*he(i)*P12e(i) + 0.5*he(i)*ue(i)*ve(i) ;     % 
     hK22e(i) =   0.5*he(i)*P22e(i) + 0.5*he(i)*ve(i)*ve(i)  ;      % 
         
  else
     ue(i) =  0.1   ;  
     ve(i) = -0.2     ;
     he(i) =  0.02 ; 
     hue(i) = he(i)*ue(i);  
     hve(i) = he(i)*ve(i);
     be(i) = 0.  ;
     ce(i) = 0.0  ; 
     
     P11e(i) = 4E-02 ;  
     P12e(i) = 1E-08 ;  
     P22e(i) = 4E-02 ; 
     
  
     hK11e(i) =   0.5*he(i)*P11e(i) + 0.5*he(i)*ue(i)*ue(i) ;    %  
     hK12e(i) =   0.5*he(i)*P12e(i) + 0.5*he(i)*ue(i)*ve(i) ;     % 
     hK22e(i) =  0.5*he(i)*P22e(i)  + 0.5*he(i)*ve(i)*ve(i) ;    % 
       end
          
       end
  
end


% Initial condition for the unknown variables
c   = ce  ;
h   = he   ; 
u   = ue   ;
v   = ve   ;
hu  = hue ;
hv  = hve;
hc  = hce    ;
    
b   = be     ;
q   = qe     ;
 
% Initial condition for the shear stress tensor components

P11  = a*P11e  ; % a = 0 for shallow water model and a = 1 for distortion shallow water model
P12  = a*P12e  ;
P22  = a*P22e  ;

% Initial condition for the kinetic energy

hK11 = hK11e ;
hK12 = hK12e ;
hK22 = hK22e ;

% U is our unknown matrix. U(1,:,1)=h; U(1,:2)=hu; U(1,:,3)=hv; U(1,:,4)=b;

U(1,:,1) = h    ;
U(1,:,2) = hu   ; 
U(1,:,3) = hv   ; 
U(1,:,4) = hK11 ;
U(1,:,5) = hK12 ;
U(1,:,6) = hK22 ;   
U(1,:,7) = hc   ;
U(1,:,8) = b    ;
  
  
% vectors for circular shift
shiftp1 = circshift((1:length(x))',1);
shiftm1 = circshift((1:length(x))',-1);


% Initialisation parameters

t     = 0;
dt    = 0.001;
ii    = 1;   %temps initial et temps final
% nombre de round
numplots = 3;                          
tplot    = [0.06; tstop];

% lecture des solutionsg?n?r?es
%Uplot          = zeros([U; zeros(length(tplot)+1, nx)]);
Uplot(1,:,:,1) = U;
styles         = {'k:','k--','k-'};
k             = 0;

B1    = zeros(1, nx, 8); 
B2    = zeros(1, nx, 8); 
B3    = zeros(1, nx, 8);
B1_i    = zeros(1, nx, 8); 
B2_i    = zeros(1, nx, 8); 
B3_i    = zeros(1, nx, 8);
Fx    = zeros(1, nx, 8);
fluxx = zeros(1, nx, 8);

ajm_plus_x = zeros(1, nx);
ajm_moin_x = zeros(1, nx);

B1_fluct  = zeros(1, nx, 8); 
B2_fluct  = zeros(1, nx, 8); 
B3_fluct  = zeros(1, nx, 8); 

B1_fluct_i  = zeros(1, nx, 8); 
B2_fluct_i  = zeros(1, nx, 8); 
B3_fluct_i  = zeros(1, nx, 8); 

kmax = 10000;
     

while (t < tstop) && (k < kmax)
    
    % initialisation des inconnues
    
    told = t;
    t    = t + dt;

    % passage aux variables non conservatives 
   
    b   =   U(1,:,8) ;
    h   =   U(1,:,1) ;
    u   =   U(1,:,2)./h;
    v   =   U(1,:,3)./h;
    hc  =   U(1,:,7) ;
    c   =   hc./h ;
    huc =   h.*u.*c;
    hvc =   h.*v.*c;
    huv =   U(1,:,2).*U(1,:,3)./U(1,:,1);
    hvv =   U(1,:,3).*U(1,:,3)./U(1,:,1);
    huu =   U(1,:,2).*U(1,:,2)./U(1,:,1);
    hu  =   U(1,:,2)                   ;
    hv  =   U(1,:,3)                   ;
    
    hK11 = U(1,:,4)                   ; 
    hK12 = U(1,:,5)                   ;
    hK22 = U(1,:,6)                   ;
    
    K11 = U(1,:,4)./U(1,:,1);
    K12 = U(1,:,5)./U(1,:,1);
    K22 = U(1,:,6)./U(1,:,1);
    
    P11 =  a*(2*K11 - (U(1,:,2)./U(1,:,1)).^2) ;
    P12 =  a*(2*K12 - (U(1,:,3)./U(1,:,1)).*(U(1,:,2)./U(1,:,1)));
    P22 =  a*(2*K22 - (U(1,:,3)./U(1,:,1)).^2) ;
    
    hP11  =  a*h.*P11 ;
    hP22  =  a*h.*P22 ;
    hP12  =  a*h.*P12 ;

   % SECOND ORDER PROCEDURE 
   % We perform a reconstruction procedure to achieve second order of
   % accuracy. We can achieve first and second order accuracy using the
   % same procedure. However, several other procedures can be used. 
   
    Generalized_AENO = 0                                   ;
    classical_AENO = 1                                     ;
    Limiter =  Generalized_AENO ;
   % Limiter = classical_AENO ; 
     
     order = First_order; 
        switch order

            case Second_order;
                b = 1 ;
                
            case First_order
                b=0 ;
        end 
 
   dU     = (-U+U(1,shiftm1,:))./dx                               ;
   dUm    = dU(1,shiftm1,:)                                      ; 
                                    
   sdU    = sign(dU)                                             ;   
   sdUm   = sign(dUm)                                            ;  
   dUmp   = (U(1,shiftm1,:) - U(1,shiftp1,:))./2*dx                ;
  
   dUp    = (U-U(1,shiftp1,:))./dx                              ;
   dUpp   = dU(1,shiftp1,:)                                      ;
   sdUp   = sign(dUp)                                            ;   
   sdUpp  = sign(dUpp)                                           ;
   grad_p = min(sdUp, sdUpp)                                     ;
        
       switch Limiter 
         case classical_AENO   
           
   grad   = min(abs(dU), abs(dUm))                                       ;
   Um     = U(1,shiftm1,:)      +     0.5*b*dx.*grad.*(sdU + sdUm)/2     ; 
   Up     = U                   +     0.5*b*dx.*grad.*(sdU + sdUm)/2     ;
            
        case Generalized_AENO
   
   grad    = min(abs(dU), abs(dUm))                                       ;
   theta   = 1.;
   
   minU = zeros(1,nx) ;  maxU = zeros(1,nx) ;
   minmod = zeros(1,nx) ; minmod1 = zeros(1,nx) ;  minmod2 = zeros(1,nx) ;
   for i=1:nx   
       
   minU(i)    = min(theta*dU(i), dUmp(i)) ;    
   minmod1(i) = min(minU(i), theta*dUp(i)) ;
   
   maxU(i)    = max(theta*dU(i), dUmp(i)) ;    
   minmod2(i) = max(maxU(i), theta*dUp(i)) ;   
    if   dU(i)>0
       minmod(i) = minmod1(i) ; 
    end
       
      if dU(i)<0
       minmod(i) = minmod2(i) ;
      else
       minmod(i) =0 ;

      end
   end
     
   minmod_M = minmod(1,shiftm1,:) ;
    
   Um     = U(1,shiftm1,:)      -     0.5*b*dx.*minmod_M(i)                  ; 
   Up     = U(1,:,:)            +     0.5*b*dx.*minmod(i)                  ;
    
        end
           
    bm     =   U(1,shiftm1,8) ;
    hm     =   U(1,shiftm1,1) ;
    um     =   U(1,shiftm1,2)./hm; uoldm = um ;
    vm     =   U(1,shiftm1,3)./hm; voldm = vm ;
    hcm    =   U(1,shiftm1,7) ;
    cm     =   hcm./hm            ;
    hum    =   U(1,shiftm1,2) ;
    hvm    =   U(1,shiftm1,3) ; 
    hK11m  =   U(1,shiftm1,4)   ;
    hK12m  =   U(1,shiftm1,5) ; 
    hK22m  =   U(1,shiftm1,6) ; 
    
    bpp     =   U(1,shiftp1,8) ;
    hpp     =   U(1,shiftp1,1) ;
    upp     =   U(1,shiftp1,2)./hpp; uoldpp = upp ;
    vpp     =   U(1,shiftp1,3)./hpp; voldpp = vpp ;
    hcpp    =   U(1,shiftp1,7) ;
    cpp     =   hcpp./hpp      ;
    hupp    =   U(1,shiftp1,2) ;
    hvpp    =   U(1,shiftp1,3) ; 
    hK11pp  =   U(1,shiftp1,4)   ;
    hK12pp  =   U(1,shiftp1,5) ; 
    hK22pp  =   U(1,shiftp1,6) ; 
   
   
   
    %Uoldm = Uold(1,shiftm1,:);
    %Uoldp = Uold(1,shiftp1,:);
    
    b_m     =   Um(1,:,8) ;
    h_m     =   Um(1,:,1) ;
    u_m     =   Um(1,:,2)./h_m; uold_m = u ;
    v_m     =   Um(1,:,3)./h_m; vold_m = v ;
    hc_m    =   Um(1,:,7) ;
    c_m     =   hc_m./h_m          ;
    hu_m    =   Um(1,:,2) ;
    hv_m    =   Um(1,:,3) ; 
    hK11_m  =   Um(1,:,4) ;
    hK12_m  =   Um(1,:,5) ; 
    hK22_m  =   Um(1,:,6) ; 
    
    
    b_p     =   Up(1,:,8) ;
    h_p     =   Up(1,:,1) ;
    u_p     =   Up(1,:,2)./h; uoldp = u ;
    v_p     =   Up(1,:,3)./h; voldp = v ;
    hc_p    =   Up(1,:,7) ;
    c_p     =   hc_p./h_p ;
    hu_p    =   Up(1,:,2) ;
    hv_p    =   Up(1,:,3) ; 
    hK11_p  =   Up(1,:,4) ;
    hK12_p  =   Up(1,:,5) ; 
    hK22_p  =   Up(1,:,6) ; 
    
  % calculate lambda = u +/- sqrt(3P11 + gh) used for finding flux

  
    hs = 0.5*( h_p + h_m);
    cs = 0.5*( c_p + c_m);
     
      
    switch Model
        case SW
    us = 0.5*( u_p + u_m);
    vs = 0.5*( v_p + v_m);
    hu_moy = 0.5*( hu_m    +  hu_p) ;
    hv_moy = 0.5*( hv_m    +  hv_p) ;
    umoy   = 0.5*( u_m     +  u_p ) ;
    vmoy   = 0.5*( v_m     +  v_p ) ;
        case SSW
    us = 0.5*( u + um);
    vs = 0.5*( v + vm);       
    hu_moy = 0.5*( hu    +  hum) ;
    hv_moy = 0.5*( hv    +  hvm) ;
    umoy   = 0.5*( u     +  um ) ;
    vmoy   = 0.5*( v     +  vm ) ;
    
    end
    lambdaup_E =   us + sqrt(3*P11 + g*hs ) ;
    lambdaum_E =   us - sqrt(3*P11 + g*hs ) ;
    
    hp = 0.5*( h_p + h_m);
    
    switch Model
        case  SW
    up = 0.5*( u_p + u_m);
        case SSW
    up = 0.5*( u   + um);     
    end
    
    lambdaup_W = up + sqrt(3*P11 + g*hp ) ;
    lambdaum_W = up - sqrt(3*P11 + g*hp ) ;

    %       calcul des speeds  ai+1/2 and ai-1/2

    ajm_plusm_x = max([lambdaup_E(:); lambdaup_W(:)]);
    ajm_moinm_x = min([lambdaum_W(:); lambdaum_E(:)]);

    aa = max(ajm_plusm_x, -ajm_moinm_x) ;

    for i=1:nx
        ajm_plus_x(1,i) = max(lambdaup_E(i), lambdaup_W(i) );
        ajm_moin_x(1,i) = min(lambdaum_W(i), lambdaum_E(i) );
      %  ajm_moin_x(1,i) = min(ajm_moin_x(1,i), 0.0);
    end
    
    % Calcul des gradients
    
   switch Model
       case SW
    dh_dx    =  ( h_p  - h_m )./(dx) ; 
    db_dx    =  ( b_p  - b_m )./(dx) ;
    dc_dx    =  ( c_p  - c_m )./(dx) ;
       case SSW
    dh_dx    =  ( h_p  - h_m )./(dx) ; 
    db_dx    =  ( b_p  - b_m )./(dx) ;
    dc_dx    =  ( c_p  - c_m )./(dx) ;
   end

    % disretisation par les path valeur moyennes h,  hu, hv 
    % Discretization of the mean values of h, hu and hv  from nonconservative terms
    % in each edge  (i+1/2)^+  +  (i+1/2)^-
    
    h_moy  = 0.5*( h_m     +  h_p ) ;
    h2_moy = 0.5*((h_m).^2 +  h_p.^2) ;
  
    
    
   % Time step of the scheme
   
   dt     = cfl*dx/min(abs(aa)) ;              

   % adjust dt to produce plots at the right time

    if (ii<=length(tplot) && tplot(ii)>=told && tplot(ii)<=t+dt)
        dt = tplot(ii)-t;
        ii = ii + 1;
    end

    % CALCULATION  OF COMPONENTS OF PHYSICAL FLUX
    %--------------------------------------------- Friction source terms coefficients
    nn    = 0.03 ;
    Cf    = g*nn^2./(h.^(1/3)) ;  
    uee    = sqrt(Cf.*(u.^2 + v.^2))  ;
    d50   = 0.001;   %d50   =  0.008            ; 
    rhos  = 2650                   ;
    rhow  = 1000                   ;
    delta_rho = rhos - rhow        ;
    s     = rhos/rhow              ;   
    nu    = 1.2*10^(-6)            ;
    De    = d50*((s-1)*g/nu^2)^1/3 ; 
    
    %  Sediment fall Velocity 
    Ws    = sqrt((13.95*nu/d50)^2 + 1.09*g*s*d50)-13.95*nu/d50 ;
    
    Karman = 0.41;  % Von Karman number
    Theta = (uee)/((s-1)*d50);      %  Shields Parameter 
   
 
   
     % CALCULATION OF Erosion, Deposition, Sediment/water Entrainment source terms

     % pour E
     
     Theta_cr  = 0.045  ;  %critical Shields parameter; 
    %  E = zeros(nx) ;  
     ttheta = Theta  ;
     phi  =  0.015   ;     % Erosion force parameter
     
     if ttheta > Theta_cr
     E = phi.*(Theta-Theta_cr).*(h).^-1.*sqrt(u.^2 + v.^2).*d50.^(-0.2) ;
     else
     E=0.0 ;
     end
  
    % pour D
    beta = min(2,(1-pore)./(c))                ;
    D    =   beta.*(c).*Ws.*(1 - beta.*(c)).^(2)  ;  % m=2 hidding coefficient
  
    rhom = rhow*(1-c) + rhos*c ; 
    rho_0 = rhow*(pore) + rhos*(1-pore)  ;
    
    rhom_moy = (rhom) ;
    
    Ri   = rhom*g.*U(:,:,1)./(u.^2) ; 
    Ew   = 0.00153./(0.0204 + Ri);           

    dFu  = Ew.*sqrt((U(:,:,2).^2./U(:,:,1).^2 + U(:,:,3).^2./U(:,:,1).^2)) ; 
 
    %----------------------------------------------------------------------------

    % speed waves propagation ratio in sense of CU technique
    
    dif1_x    = ajm_plus_x./(ajm_plus_x-ajm_moin_x) ;
    dif1_x_m  = ajm_moin_x./(ajm_plus_x-ajm_moin_x) ;  


    %% CALCULATION OF NONCONSERVATIVE TERMS 
  
     h2_moy_i   = 0.5*( hm.^2  + hpp.^2 )                             ;
    
     switch Model
    
         case SW
      h_moy_i   = 0.25*(hm + hpp) ;  
     hu_moy_i   = 0.5*( hupp + hum)                                       ;
     hv_moy_i   = 0.5*( hvpp + hvm)                                   ;
     umoy_i     = 0.5*( upp + um )                                    ;
     vmoy_i     = 0.5*( vpp + vm )                                    ; 
         case SSW
     h_moy_i    = 0.25*(h_m + h_p) ;  % Well balanced strategy
     hu_moy_i   = 0.5*( hu + hum)                                     ;
     hv_moy_i   = 0.5*( hv + hvm)                                     ;
     umoy_i     = 0.5*( u + um )                                      ;
     vmoy_i     = 0.5*( v + vm )                                      ;      
     
     end 
     
     rhom_p     = rhow*(1-cpp) + rhos*cpp                              ; 
     rhom_m     = rhow*(1-cm)  + rhos*cm                               ; 
     rhomp     = rhow*(1-c_p) + rhos*c_p                               ; 
     rhomm     = rhow*(1-c_m)  + rhos*c_m                              ; 
     rho_0      = rhow*(pore)  + rhos*(1-pore)                         ;
     rhom_moy_i = 0.5*(rhom_m +rhom_p)                                 ;
     rhom_moy = 0.5*(rhomm +rhomp)                                     ;
     Ag = 0.009 ;
     ub   = -2*g*Ag*umoy/(1-pore) ;  % see Ngatcha et al (2024) in Taylor Francis
           
    for is = 1:nx
        B1_i(1,is,1:8) = 0.0;  B2_i(1,is,1:8) = 0.0;  B3_i(1,is,1:8) = 0.0; 
        
        % For water depth
        B1_i(1,is,2) = g*h_moy_i(1,is)    ;
        B1_i(1,is,4) = g*hu_moy_i(1,is)     ;
        B1_i(1,is,5) = 0.5*g*hv_moy_i(1,is) ;
      
        % for bed evolution
        B2_i(1,is,2) = g*h_moy_i(1,is)      ;
        B2_i(1,is,4) = g*hu_moy_i(1,is)     ;
        B2_i(1,is,5) = 0.5*g*hv_moy_i(1,is) ;
        B2_i(1,is,8) = ub(1,is);
       % for sediment concentration
   
        B3_i(1,is,2)  = (g*(h_moy_i(1,is))^2*delta_rho)/(2*rhom_moy_i(1,is))  ;
        B3_i(1,is,4)  = ((g*(h_moy_i(1,is))^2*delta_rho)/(2*rhom_moy_i(1,is)))*umoy_i(1,is) ;
        B3_i(1,is,5)  = 0.5*((g*(h_moy_i(1,is))^2*delta_rho)/(2*rhom_moy_i(1,is)))*vmoy_i(1,is) ;
    end
    
    
     switch  Model

         case SW;
    
        dh_dx_i    =  ( hpp  -  hm )./(dx) ; 
        db_dx_i    =  ( bpp  -  bm )./(dx) ;
        dc_dx_i    =  ( cpp  -  cm )./(dx) ; 
        
         case SSW
        dh_dx_i    =  ( hpp  -  hm )./(dx) ; 
        db_dx_i    =  ( bpp  -  bm )./(dx) ;
        dc_dx_i    =  ( cpp  -  cm )./(dx) ; 
     end
     
   for is=1:nx
        B1_fluct_i(1,is,:) = (B1_i(1,is,:))*dh_dx_i(1,is) ;
        B2_fluct_i(1,is,:) = (B2_i(1,is,:))*db_dx_i(1,is)  ;
        B3_fluct_i(1,is,:) = (B3_i(1,is,:))*dc_dx_i(1,is)  ;    
   end  
    
    
    % Calcul des termes nonconsertifs en h sur les ar?tes    
    for is = 1:nx
        B1(1,is,1:8) = 0.0;  B2(1,is,1:8) = 0.0;  B3(1,is,1:8) = 0.0; 
        
        % For water depth
        B1(1,is,2) = g*h_moy(1,is)  ;
        B1(1,is,4) = g*hu_moy(1,is) ;
        B1(1,is,5) = 0.5*g*hv_moy(1,is) ;
      
        % for bed evolution
        B2(1,is,2) = g*h_moy(1,is)  ;
        B2(1,is,4) = g*hu_moy(1,is) ;
        B2(1,is,5) = 0.5*g*hv_moy(1,is) ;
        B2(1,is,8) = ub(1,is);
       % for sediment concentration
   
        B3(1,is,2)  = (g*(h_moy(1,is))^2*delta_rho)/(2*rhom_moy(1,is))  ;
        B3(1,is,4)  = ((g*(h_moy(1,is))^2*delta_rho)/(2*rhom_moy(1,is)))*umoy(1,is) ;
        B3(1,is,5)  = 0.5*((g*(h_moy(1,is))^2*delta_rho)/(2*rhom_moy(1,is)))*vmoy(1,is) ;
    end
       
    B1m = B1(1,shiftp1,:) ;
    B2m = B2(1,shiftp1,:) ;
    B3m = B3(1,shiftp1,:) ;
   
    % Fluctuations in each meshes edge
    % H_psi = (a^+/a^+-a^-)*Hi+1/2 - (a^+/a^+-a^-)*Hi-1/2

   for is=1:nx
        B1_fluct(1,is,:) = (dif1_x(1,is)*(B1(1,is,:)) ...
            -  dif1_x_m(1,is)*(B1m(1,is,:)))*dh_dx(1,is) ;
        B2_fluct(1,is,:) = (dif1_x(1,is)*(B2(1,is,:)) ...
            -  dif1_x_m(1,is)*(B2m(1,is,:)))*db_dx(1,is) ;
        B3_fluct(1,is,:) = (dif1_x(1,is)*(B3(1,is,:)) ...
            -  dif1_x_m(1,is)*(B3m(1,is,:)))*dc_dx(1,is) ;    
   end
   
    B1m = B1(1,shiftp1,:) ;
    B2m = B2(1,shiftp1,:) ;
    B3m = B3(1,shiftp1,:) ;

    % Flux components    
    
    FK11_x =  K11.*hu + hu.*P11 ;           
    FK12_x =  2*K12.*hu  + hv.*K11 - hu.*u.*v ;     
    FK22_x =  K22.*hu + hv.*P12 ;    
    FK88   = 0.*hu  ;

    for is=1:nx
        Fx(1,is,1) = hu(1,is) ;
        Fx(1,is,2) = huu(1,is) + hP11(1,is) ;
        Fx(1,is,3) = huv(1,is) + hP12(1,is) ;
        Fx(1,is,4) = FK11_x(1,is) ;
        Fx(1,is,5) = FK12_x(1,is) ;
        Fx(1,is,6) = FK22_x(1,is) ;
        Fx(1,is,7) = huc(1,is) ;
        Fx(1,is,8) = 0*FK88(1,is) ;
    end

    Fxm = Fx(1,shiftm1,:) ; % F(W^+(i+1/2))
   % Um  = U(1,shiftm1,:) ;   % W^+(i+1/2)

    
     %Erosion/deposition source terms
     
     uv = huv./h ; vv = hvv./h ; uu = huu./h ;
     
       switch  Model

              case SSW;
     
     S_sed_echange = cat(3, (E-D)/(1-pore),...
      (E-D).*u/(1-pore),...
      (E-D).*v/(1-pore), ...
     ((E-D).*u.^2)/(1-pore),...
     (0.5*(E-D).*uv)/(1-pore),...
     (E-D).*vv/(1-pore),...
      (E-D) ,...
     -(E-D)/(1-pore)); 

 
     SU_friction = cat(3, zeros(size(x)),...
      -Cf.*(u).*(u.^2 +v.^2 + P11),...
      -Cf.*(U(:,:,3)./U(:,:,1)).*(U(:,:,2).^2./U(:,:,1).^2 +...
         U(:,:,3).^2./U(:,:,1).^2 ),...
      -Cf.*(U(:,:,2)./U(:,:,1)).^2.*(U(:,:,2).^2./U(:,:,1).^2 +...
         U(:,:,3).^2./U(:,:,1).^2 ) ,...
       (-0.5*Cf.*(U(:,:,2)./U(:,:,1)).*(U(:,:,2)./U(:,:,1))...
       -0.5*Cf.*(U(:,:,3)./U(:,:,1)).*(U(:,:,2)./U(:,:,1))).*(U(:,2).^2./U(:,:,1).^2 +...
         U(:,:,3).^2./U(:,:,1).^2) ,...
         -0.5*Cf.*(U(:,:,3)./U(:,:,1)).^2.*(U(:,2).^2./U(:,:,1).^2 +...
         U(:,:,3).^2./U(:,:,1).^2),...
         zeros(size(x)),...
         zeros(size(x))) ; 
 
 
            case SW       
                
     S_sed_echange = cat(3, (E-D)/(1-pore),...
      (E-D).*u./((1-pore)),...
     zeros(size(x)), ...
     zeros(size(x)),...
     zeros(size(x)),...
     zeros(size(x)),...
      (E-D) ,...
     -(E-D)/(1-pore)); 

     SU_friction = cat(3, zeros(size(x)),...
      -Cf.*(u).^2 ,...
      zeros(size(x)),...
      zeros(size(x)) ,...
       zeros(size(x)),...
       zeros(size(x)) ,...
         zeros(size(x)),...
         zeros(size(x))) ;
         
        end
     
    % calculation of CU fluxes at i+1/2
    for is = 1:nx
        Sr  = ajm_plus_x(is);
        Sl  = ajm_moin_x(is);
        
        fluxx(1,is,:) = (Sr*Fx(1,is,:) - Sl*Fxm(1,is,:) ...
                          - Sl*Sr*( Up(1,is,:) - Um(1,is,:)))/(Sr-Sl);
                      
                      
    end
    
      % STRESS GRADIENT DISCRETIZATION 
       k = 0.1 ; nu = 0.0000012 ; 
       Dgradx_P11  = k*(P11(:,shiftm1,1) - 2*P11 + P11(:,shiftp1,1))./(2*dx.^2) ;

       S_w_echange = cat(3, dFu, u.*dFu, u.*u.*dFu, ...
       zeros(size(x)), zeros(size(x)), zeros(size(x)), zeros(size(x)), zeros(size(x)));
   
       S_P_echange = cat(3, zeros(size(x)), zeros(size(x)),  Dgradx_P11 , ...
       zeros(size(x)), zeros(size(x)), zeros(size(x)), zeros(size(x)), 2*g*Ag/(1-pore)*Dgradx_P11);

       Diplus = (fluxx - fluxx(:,shiftp1,:))/dx    ;


    %Update solution with multiple source  terms 
    
    U =  U  - dt*(Diplus - B1_fluct - b*B1_fluct_i + Cla*(-B2_fluct - B3_fluct - b*B3_fluct_i- b*B2_fluct_i)) + ...
        Cla*(dt*S_sed_echange + dt*SU_friction + dt*S_P_echange + dt*S_w_echange) ;
            
    k = k+1;
    
    
    % impose boundary conditions on all the variables
    U(1,end,1) =  U(1,end-1,1);   U(1,1,1)   =  U(1,2,1);  
    % on hu
    U(:,end,2) =  U(:,end-1,2); U(:,1,2) =  U(:,2,2);
    % on hv
    U(:,end,3) = U(:,end-1,3);  U(:,1,3) = U(:,2,3);
    % on hK11
    U(:,end,4) =  U(:,end-1,4); U(:,1,4) =  U(:,2,4);
    % on hK12
    U(:,end,5) =  U(:,end-1,5); U(:,1,5) =  U(:,2,5);
    % on hK22
    U(:,end,6) =  U(:,end-1,6); U(:,1,6) =  U(:,2,6);
    % on c
    U(:,end,7) =  U(:,end-1,7); U(:,1,7) =  U(:,2,7);
    % on b
    U(:,end,8) =  U(:,end-1,8); U(:,1,8) =  U(:,2,8);

    
    % Updated solutions after applied boundary conditions
    
    % u = hu./h;             v = hv./h;
 
    h   = U(1,:,1)          ;
    b   = U(1,:,8)          ;
    u   = U(1,:,2)./U(1,:,1);
    v   = U(1,:,3)./U(1,:,1);
    c   = U(1,:,7)./U(1,:,1);
    
    K11 = 0.5*u.^2  + a*P11;
    K12 = 0.5*u.*v  + a*P12;
    K22 = 0.5*v.^2  + a*P22;
    
    P11 = a*(2*K11 - u.^2)      ;
    P12 = a*(2*K12 - u.*v)      ;
    P22 = a*(2*K22 - v.^2)      ;
  
    
    
   
      % % save the reference solution for SSW obtained by refining mesh
      % N=2500
      
% %href = h ; uref = u ; vref = v ;  P11ref = P11 ; P12ref = P12 ;  P22ref = P22 ; 
% load('solution_mesh_ref', 'href','uref', 'vref', 'P11ref', 'P12ref', 'P22ref') 
% 
% 
% href_1 = h ; uref_1 = u ; vref_1 = v ;  P11ref_1 = P11 ; P12ref_1 = P12 ;  P22ref_1 = P22 ; 
% load('solution_mesh_ref_111', 'href_1','uref_1', 'vref_1', 'P11ref_1', 'P12ref_1', 'P22ref_1') 
% 
% 
% href_2 = h ; uref_2 = u ; vref_2 = v ;  P11ref_2 = P11 ; P12ref_2 = P12 ;  P22ref_2 = P22 ; 
% load('solution_mesh_ref_2', 'href_2','uref_2', 'vref_2', 'P11ref_2', 'P12ref_2', 'P22ref_2') 
% 


disp(t)

   
   xref   = 0:lg/2500:lg-lg/2500 ; % refined solution (or reference solution) for N = 2500
   xref_1 = 0:lg/400:lg-lg/400 ;   % solution for N = 400
   xref_2 = 0:lg/200:lg-lg/200 ;  % Solution for  N = 200
   
  
            
figure(1)

hold on 
plot(x,  h ,'k' ,'linewidth',2 );
grid on 
title(' \bf Water depth level','interpreter','latex','FontSize',12);
xlabel('  $\mathbf{x[m]}$','interpreter','latex','FontSize',12);
ylabel({'$ \mathbf{h[m]}$'},'interpreter','latex','FontSize',12);


figure (2)
plot(x, u, 'k' ,'linewidth',2 );
grid on 
title(' \bf x-velocity','interpreter','latex','FontSize',12);
xlabel('  $\mathbf{x[m]}$','interpreter','latex','FontSize',12);
ylabel({'$ u[m/s]$'},'interpreter','latex','FontSize',12);

   
figure (3)
plot(x, v ,'k' ,'linewidth',2 );
grid on 
title(' \bf y-Velocity','interpreter','latex','FontSize',12);
xlabel('  $\mathbf{x[m]}$','interpreter','latex','FontSize',12);
ylabel({'$ v[m/s]$'},'interpreter','latex','FontSize',12);

    


figure(4)
plot(x, P11 ,'k' ,  'linewidth',2 );
grid on 
title(' \bf Shear stress component','interpreter','latex','FontSize',12);
xlabel('  $\mathbf{x[m]}$','interpreter','latex','FontSize',12);
ylabel({'$ \mathbf{P11[m^2/s^2]}$'},'interpreter','latex','FontSize',12);



figure(5)
plot(x, P12,'k' , 'linewidth',2 );
grid on 
title(' \bf Shear stress component','interpreter','latex','FontSize',12);
xlabel('  $\mathbf{x[m]}$','interpreter','latex','FontSize',12);
ylabel({'$ \mathbf{P12[m^2/s^2]}$'},'interpreter','latex','FontSize',12);


figure(6)
plot(x,  P22 ,'k', 'linewidth',2 );
grid on 
title(' \bf Shear stress component','interpreter','latex','FontSize',12);
xlabel('  $\mathbf{x[m]}$','interpreter','latex','FontSize',12);
ylabel({'$ \mathbf{P22[m^2/s^2]}$'},'interpreter','latex','FontSize',12);




    %if (ismember(t+dt,tplot))axis([-0.01 0.01 -0.01 0.01 0 0.02])
    if (any(tplot-t-dt==0))
        Uplot(:,:,:,ii) = U;      % store U for plotting
    end
    
    

end 
