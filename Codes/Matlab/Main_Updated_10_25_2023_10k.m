% Main_program.m
% Main program for SG's modeling
% Authors: Alex Stankovic and Alex Saric
%
  clear all

intervals = [1];
for i=1:length(intervals)
  clearvars -except i     %deletes all variables except X in workspace

  intervals = 1; % Inputs  
  dt=0.02%0.0578;              % Integration step
  opt1=1;                 % 1. Infinite bus, 2. With measurements in Point of Connection (ONLY WITH model==1!)
  n=6;                    % Model order
  model=1;                % 1. Large SG from Milano's book (Infinite bus (opt1=1), or With measurements in Point of Connection (opt1=2))
                          % 2. Diesel generator 5.5MVA (see ref. [2] in SG_Models_and_Parameters.doc)
                          % 3. Large Hydraulic unit (see ref. [4] in SG_Models_and_Parameters.doc)
                          % 4. Large Thermal unit (see ref. [4] in SG_Models_and_Parameters.doc)    
                          % 5. Small SG, 0.1MVA, 0.4kV (see ref. [3] in SG_Models_and_Parameters.doc)    
  if model==1,  
     D=0.5;H=2.53;ra=0;T1d0=4.75;T1q0=1.5;T2d0=0.06;T2q0=0.21;x1d=0.928;x1q=2.86;x2d=0.48;x2q=0.48;xd=5;xl=0.536;xq=4.88;
     %D=1.2;H=7.848/2;ra=0;T1d0=7.4;T1q0=0.3;T2d0=0.03;T2q0=0.033;x1d=0.049;x1q=0.073;x2d=0.037;x2q=0.065;xd=0.146;xl=0.536;xq=0.105; % Mark's input data
     if opt1==1,
        %x06=[1.1764;0.9873;1.2778;0.3430;0.9046;0.3430]; % Initial condition of state variables (delta,omega,e1q,e1d,e2q,e2d)
        %vf0=4.1468;Pm=0.7;V=1.09;theta=-0.0383;          % 
        x06=[0.5;0.98;2.13;0.02;1.93;0.02]; % Initial condition of state variables (delta,omega,e1q,e1d,e2q,e2d)
        vf0=4.2;Pm=0.7;V=1.09;theta=0;          % Initial condition of algebraic variables        
     elseif opt1==2,
        [x0,Y_PSAT,dt,time]=initial_condition(n,t_end);
     end
  elseif model==2,
     D=0;H=0.4;ra=0.01;T1d0=4.3;T1q0=0.85;T2d0=0.03;T2q0=0.05;x1d=0.169;x1q=0.228;x2d=0.135;x2q=0.2;xd=1.79;xl=0;xq=1.71;
     x06=[1.18;0.99;1.28;0.34;1.20;0.34];                % Initial condition of state variables (delta,omega,e1q,e1d,e2q,e2d)
     vf0=4.0;Pm=0.7;V=1.0;theta=0.0;                     % Initial condition of algebraic variables        
  elseif model==3,
     D=0;H=3.0;ra=0.01;T1d0=6;T1q0=0.85;T2d0=0.03;T2q0=0.05;x1d=0.35;x1q=0.228;x2d=0.135;x2q=0.3;xd=1.0;xl=0.15;xq=0.7;
     x06=[1.18;0.99;2.05;0.20;1.55;0.20];                % Initial condition of state variables (delta,omega,e1q,e1d,e2q,e2d)
     vf0=4.0;Pm=0.7;V=1.0;theta=0.0;                     % Initial condition of algebraic variables        
  elseif model==4,
     D=0;H=3.0;ra=0.0025;T1d0=7;T1q0=1;T2d0=0.035;T2q0=0.035;x1d=0.3;x1q=0.7;x2d=0.2;x2q=0.2;xd=1.8;xl=0.15;xq=1.8;
     x06=[1.18;0.99;2.05;0.20;1.55;0.20];                % Initial condition of state variables (delta,omega,e1q,e1d,e2q,e2d)
     vf0=4.0;Pm=0.7;V=1.0;theta=0.0;                     % Initial condition of algebraic variables        
  elseif model==5,
     D=0;H=0.5;ra=0.002;T1d0=0.2;T1q0=0.1;T2d0=0.001;T2q0=0.001;x1d=0.00546;x1q=0.329;x2d=0;x2q=0;xd=1.66;xl=0.14;xq=1.66;
     x06=[0.05;0.99;1.009;0.05;1.009;0.05];              % Initial condition of state variables (delta,omega,e1q,e1d,e2q,e2d)
     vf0=4.0;Pm=0.05;V=1.0;theta=0.0;                    % Initial condition of algebraic variables        
  end
  if opt1==1,
     x0(1:n,1)=x06(1:n,1);
  end
  
  freq=60;                % Frequency
  t_end=5;               % Time for transients 
  n_perturb_param=10000;     % Number of patterns for random perturbation of parameters
  %n_perturb_init_cond=10; % Number of patterns for random perturbation of initial conditions
  interval_param=10;      % Interval of random variation of parameters (in percents around the basic value)
  %interval_init_cond=5;   % Interval of random variation of initial condition (in percents around the basic value)
   
  tspan = 0:dt:t_end;    
  options = odeset('RelTol',1e-7,'AbsTol',1e-7*ones(1,n)); 
  
%Basic parameters and initial condition 
  if opt1==1,
     [t,XX]=ode45(@(t,x) SG_model_infinite_bus(t,x,freq,vf0,Pm,V,theta,xd,x1d,x2d,T1d0,T2d0,xq,x1q,x2q,T1q0,T2q0,H,D,n),tspan,x0,options);
  elseif opt1==2,
     [t,XX]=ode45(@(t,x) SG_model_with_measurements(t,x,freq,Y_PSAT,xd,x1d,x2d,T1d0,T2d0,xq,x1q,x2q,T1q0,T2q0,H,D,time,n,opt1),tspan,x0,options);
  end
  out_perturb_param(1,:,:)=XX(1:end,1:end);
  out_perturb_init_cond(1,:,:)=XX(1:end,1:end);
  if opt1==2,
     V=Y_PSAT(2,1);
     theta=Y_PSAT(3,1);
     vf0=Y_PSAT(1,1);
     Pm=Y_PSAT(4,1);
  end
% Perturbation of parameters
Parameter_values=zeros(n_perturb_param,12);
Behaviors =  zeros(n_perturb_param,length(tspan),6);
  for ii=1:n_perturb_param,
      disp ' PERTURBATION OF PARAMETERS' 
      ii/n_perturb_param *100
      [xdp,x1dp,x2dp,T1d0p,T2d0p,xqp,x1qp,x2qp,T1q0p,T2q0p,Hp,Dp] = perturb_param(interval_param,xd,x1d,x2d,T1d0,T2d0,xq,x1q,x2q,T1q0,T2q0,H,D); %Perturbation of parameters
  if opt1==1,
     [t,XX]=ode45(@(t,x) SG_model_infinite_bus(t,x,freq,vf0,Pm,V,theta,xdp,x1dp,x2dp,T1d0p,T2d0p,xqp,x1qp,x2qp,T1q0p,T2q0p,Hp,Dp,n),tspan,x0,options);
  elseif opt1==2,
     [t,XX]=ode45(@(t,x) SG_model_with_measurements(t,x,freq,Y_PSAT,xdp,x1dp,x2dp,T1d0p,T2d0p,xqp,x1qp,x2qp,T1q0p,T2q0p,Hp,Dp,time,n,opt1),tspan,x0,options);
  end
      out_perturb_param(ii+1,:,:)=XX(1:end,1:end);
      Behaviors(ii,:,:) = XX;
      Parameter_values(ii,:) = [xdp,x1dp,x2dp,T1d0p,T2d0p,xqp,x1qp,x2qp,T1q0p,T2q0p,Hp,Dp];
  end
save ('../../Data/Model_output.mat','Behaviors','Parameter_values');
plot_perturb_param_results(n_perturb_param,out_perturb_param,tspan,n);
%   disp ' PERTURBATION OF PARAMETERS'
%   disp 'Outputs are stored in out_perturb_param(I,J,K).mat' 
% %   disp 'I=basic+perturbations;   J=time points,    K=states' 
%   fprintf ('\n Press Enter to continue \n\n');
%    File_name_behaviors =  strcat('Behaviors_interval_crop_8_seconds_',num2str(interval_param),'.csv');
%    csvwrite(File_name_behaviors,Behaviors);
%    
%    File_name_parameters =  strcat('Parameters_for_crop_8_seconds_',num2str(interval_param),'.csv');
%    csvwrite(File_name_parameters,Parameter_values)
end
