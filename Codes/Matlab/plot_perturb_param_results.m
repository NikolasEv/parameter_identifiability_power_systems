% plot_perturb_param_results.m
% Plotting for parameters perturbation results
  function plot_perturb_param_results(n_perturb_param,out_perturb_param,tspan,n)
% Authors: Alex Stankovic and Alex Saric
% 
  figure
  subplot(3,2,1)
  plot(tspan,out_perturb_param(1,:,1),'LineStyle','-','LineWidth',1.5,'Color','red'); hold on, grid on; 
  plot(tspan,out_perturb_param(2:n_perturb_param,:,1),'LineStyle','-','LineWidth',1.0,'Color','blue'); hold on, grid on; 
  xlabel('{\fontname{Times}\it Time }{\fontname{Times} [s]}'); 
  ylabel('{\fontname{Times} \it \delta}{\fontname{Times} [pu]}');
  legend ('Basic','Perturbed')
  subplot(3,2,2)
  plot(tspan,out_perturb_param(1,:,2),'LineStyle','-','LineWidth',1.5,'Color','red'); hold on, grid on; 
  plot(tspan,out_perturb_param(2:n_perturb_param,:,2),'LineStyle','-','LineWidth',1.0,'Color','blue'); hold on, grid on; 
  xlabel('{\fontname{Times}\it Time }{\fontname{Times} [s]}'); 
  ylabel('{\fontname{Times} \it \omega}{\fontname{Times} [pu]}');
  legend ('Basic','Perturbed')
  subplot(3,2,3)
  if n>=3,
     plot(tspan,out_perturb_param(1,:,3),'LineStyle','-','LineWidth',1.5,'Color','red'); hold on, grid on; 
     plot(tspan,out_perturb_param(2:n_perturb_param,:,3),'LineStyle','-','LineWidth',1.0,'Color','blue'); hold on, grid on; 
  end
  xlabel('{\fontname{Times}\it Time }{\fontname{Times} [s]}'); 
  ylabel('{\fontname{Times} \it E''_q}{\fontname{Times} [pu]}');
  legend ('Basic','Perturbed')
  subplot(3,2,4)
  if n>=4,
     plot(tspan,out_perturb_param(1,:,4),'LineStyle','-','LineWidth',1.5,'Color','red'); hold on, grid on; 
     plot(tspan,out_perturb_param(2:n_perturb_param,:,4),'LineStyle','-','LineWidth',1.0,'Color','blue'); hold on, grid on; 
  end
  xlabel('{\fontname{Times}\it Time }{\fontname{Times} [s]}'); 
  ylabel('{\fontname{Times} \it E''_d}{\fontname{Times} [pu]}');
  legend ('Basic','Perturbed')
  subplot(3,2,5)
  if n>=5,
     plot(tspan,out_perturb_param(1,:,5),'LineStyle','-','LineWidth',1.5,'Color','red'); hold on, grid on; 
     plot(tspan,out_perturb_param(2:n_perturb_param,:,5),'LineStyle','-','LineWidth',1.0,'Color','blue'); hold on, grid on; 
     legend ('Basic','Perturbed')
  end
  xlabel('{\fontname{Times}\it Time }{\fontname{Times} [s]}'); 
  ylabel('{\fontname{Times} \it E''''_q}{\fontname{Times} [pu]}');
  subplot(3,2,6)
  if n>=6,
     plot(tspan,out_perturb_param(1,:,6),'LineStyle','-','LineWidth',1.5,'Color','red'); hold on, grid on; 
     plot(tspan,out_perturb_param(2:n_perturb_param,:,6),'LineStyle','-','LineWidth',1.0,'Color','blue'); hold on, grid on; 
     legend ('Basic','Perturbed')
  end
  xlabel('{\fontname{Times}\it Time }{\fontname{Times} [s]}'); 
  ylabel('{\fontname{Times} \it E''''_d}{\fontname{Times} [pu]}');
  %suptitle('Perturbation of parameters') 
