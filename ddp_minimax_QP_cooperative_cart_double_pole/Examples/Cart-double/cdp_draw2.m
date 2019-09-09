function cdp_draw2(xp, theta1, theta2, force,trial_num)

scale = 1;

l = 0.3*scale;
xmin = -2*scale;
xmax = 2*scale;
height = 0.07*scale;
width  = 0.25*scale;

font_size = 12;

% Compute positions
cart = [ xp + width,  height
  xp + width, -height
  xp - width, -height
  xp - width,  height
  xp + width,  height ];
pend2 = [xp, 0;
  xp-2*l*sin(theta1), cos(theta1)*2*l];
pend3 = [xp-2*l*sin(theta1), 2*l*cos(theta1);
  xp-2*l*sin(theta1)-2*l*sin(theta2), 2*l*cos(theta1)+2*l*cos(theta2)];

% plot cart double pendulum
clf
hold on

% plot(0,4*l,'k+','MarkerSize',2*font_size,'linewidth',2);
plot([xmin, xmax], [-height-0.03*scale, -height-0.03*scale], ...
  'Color','k','LineWidth',3);
% plot([0 force/20*xmax],[-0.3, -0.3].*scale, 'Color', 'g', 'LineWidth', font_size);


% Draw Cart
plot(cart(:,1), cart(:,2),'Color','k','LineWidth',1);           
fill(cart(:,1), cart(:,2),'g');
% Draw Pendulum2
plot(pend2(:,1), pend2(:,2),'Color','b','LineWidth', round(font_size/2)); 
 % Draw Pendulum3
plot(pend3(:,1), pend3(:,2),'Color','b','LineWidth', round(font_size/2));
% joint at cart
plot(xp,0,'o','MarkerSize', round((font_size+4)/2),'Color','y','markerface','y'); 
% 2nd joint
plot(pend3(1,1),pend3(1,2),'o','MarkerSize', ...
  round((font_size+4)/2),'Color','y','markerface','y');
% tip of 2nd joint
plot(pend3(2,1),pend3(2,2),'o','MarkerSize', ...
  round((font_size+4)/2),'Color','y','markerface','y'); 


set(gca,'DataAspectRatio',[1 1 1],'XLim',[xmin xmax],'YLim',[-1.4 1.4].*scale);

axis off
text = ['Trial#' num2str(trial_num)];
title(text);
drawnow;
