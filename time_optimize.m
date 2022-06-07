clear all;close all;
x0 = [0.05,0.025];
lb = [0,0];
ub = [0.15,0.1];
C1 = 0.001;
C2 = 0.01;
options = optimoptions('ga','PlotFcn', @gaplotbestf,'UseParallel',true);
sol = ga(@(x) swingTime(x,C1,C2),2,[],[],[],[],lb,ub,[],options);
%% test
max_time = 100;
dt = 1e-3;
max_cnt = max_time/dt;
l1 = sol(1);
l2 = sol(2);
lDE = 0.15;
lGF = 0.1;
h = 0.5;
m1 = 0.6;%kg
m2 = 0.5;%kg
J1 = 1/12*m1*lDE^2;%kgm^2
J2 = 1/12*m2*lGF^2;%kgm^2
g = 9.8;%m/s^2;
theta0 = [0,pi/2];
gif_delay = 0.05;

% initial condition
res = fkine2bar(theta0,l1,l2);
q0 = [res(1:2),theta0(1),res(3:end),theta0(2)];
dq0 = zeros(1,size(q0,2));

% acceleration
q = zeros(1,size(q0,2));
q(1,:) = q0;
dq = zeros(1,size(dq0,2));
dq(1,:) = dq0;
ddq = zeros(1,size(dq0,2));
M = diag([m1,m1,J1,m2,m2,J2]);
energy = zeros(1,1);
energy(1) = dq(1,:)*M*dq(1,:)'/2+q(1,2)*m1*g+q(1,5)*m2*g;
if 0.5*lDE-l1>=0 && 0.5*lGF-l2>=0
    min_energy(1) = (h+0.5*lDE-l1)*m1*g+(h-l1-0.5*lGF+l2)*m2*g;
    min_energy(2) = (h-0.5*lDE+l1)*m1*g+(h+l1-0.5*lGF+l2)*m2*g;
    min_energy = min(min_energy);
else
    if 0.5*lDE-l1>=0 && 0.5*lGF-l2<0
        min_energy(1) = (h+0.5*lDE-l1)*m1*g+(h-l1+0.5*lGF-l2)*m2*g;
        min_energy(2) = (h-0.5*lDE+l1)*m1*g+(h+l1+0.5*lGF-l2)*m2*g;
        min_energy = min(min_energy);
    else
        if 0.5*lDE-l1<0 && 0.5*lGF-l2>=0
            min_energy = (h+0.5*lDE-l1)*m1*g+(h-l1-0.5*lGF+l2)*m2*g;
        else
            min_energy = (h+0.5*lDE-l1)*m1*g+(h-l1+0.5*lGF-l2)*m2*g;
        end
    end
end

i = 1;
while true
    theta1  = q(i, 3);dtheta1  = dq(i, 3);
    theta2 = q(i, 6);dtheta2 = dq(i, 6);
    
    Phi_q = [ 1,  0, -sin(theta1)*(l1 - lDE/2), 0, 0,                         0;
        0,  1,  cos(theta1)*(l1 - lDE/2), 0, 0,                         0;
        -1,  0,      -(lDE*sin(theta1))/2, 1, 0, -sin(theta2)*(l2 - lGF/2);
        0, -1,       (lDE*cos(theta1))/2, 0, 1,  cos(theta2)*(l2 - lGF/2)];
    gamma = [dtheta1^2*cos(theta1)*(l1 - lDE/2);
        dtheta1^2*sin(theta1)*(l1 - lDE/2);
        (lDE*cos(theta1)*dtheta1^2)/2 + cos(theta2)*(l2 - lGF/2)*dtheta2^2;
        (lDE*sin(theta1)*dtheta1^2)/2 + sin(theta2)*(l2 - lGF/2)*dtheta2^2];
    %     hat_gamma = gamma - 2*alpha*(Phi_q*dq(i,:)')-beta^2*Phi([q(i,1:2),q(i,4:end)],q(i,3))';
    A = [M,Phi_q';Phi_q,zeros(size(Phi_q,1))];
    %     if i > max(length(t)/20,1)
    %         M_d = 0;
    %     end
    Q = [0;-m1*g;-C1*dtheta1;0;-m2*g;-C2*(dtheta2-dtheta1)];
    %     B = [Q;hat_gamma];
    B = [Q;gamma];
    ddq_temp = A\B;
    ddq(i,:) = ddq_temp(1:size(q,2))';
    if i < 2
        dq(i+1,:) = dq(i,:) + ddq(i,:)*dt;
        q(i+1,:) = q(i,:) + dq(i,:)*dt;
    else
        if i < 4
            dq(i+1,:) = dq(i,:) + (3/2*ddq(i,:)-1/2*ddq(i-1,:))*dt;
            q(i+1,:) = q(i,:) + (3/2*dq(i,:)-1/2*dq(i-1,:))*dt;
        else
            dq(i+1,:) = dq(i,:) + (55*ddq(i,:)-59*ddq(i-1,:)+37*ddq(i-2,:)-9*ddq(i-3,:))*dt/24;
            q(i+1,:) = q(i,:) + (55*dq(i,:)-59*dq(i-1,:)+37*dq(i-2,:)-9*dq(i-3,:))*dt/24;
        end
    end
    % choose theta as independent variable
    % position constraint
    res = fkine2bar([q(i+1,3),q(i+1,6)],l1,l2);
    q(i+1,:) = [res(1:2),q(i+1,3),res(3:end),q(i+1,6)];
    % velocity constraint
    Phi_v = [Phi_q(:,3),Phi_q(:,6)];
    Phi_u = [Phi_q(:,1:2),Phi_q(:,4:5)];
    b = -Phi_v*[dq(i+1,3);dq(i+1,6)];
    res = (Phi_u\b)';
    dq(i+1,:) = [res(1:2),dq(i+1,3),res(3:end),dq(i+1,6)];
    
    %     position_error(i+1,:) = norm(Phi([q(i+1,1:2),q(i+1,4:end)],q(i+1,3)));
    %     velocity_error(i+1,:) = norm(Phi_q*dq(i+1,:)');
    energy(i+1,:) = dq(i+1,:)*M*dq(i+1,:)'/2+q(i+1,2)*m1*g+q(i+1,5)*m2*g;
    if energy(i+1,:)<min_energy*1.001 || i> max_cnt
        swing_time = -i*dt;
        break
    end
    i = i + 1;
end

figure()
plot(energy)
hold on
plot(repmat(min_energy,size(energy)))
disp(['went on for ' num2str(-swing_time) ' sec(s)'])

%% position plotting
fig = figure;
counter = 1;
step_size = gif_delay/dt;
im = [];
for i = 1:step_size:size(q,1)
    plotMechanism(q(i,:),i*dt);
    drawnow
    F = getframe(fig);
    im{counter} = frame2im(F);
    counter = counter + 1;
end

% save file
filename = 'testAnimated.gif'; % Specify the output file name
for idx = 1:size(im,2)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',gif_delay);
    end
end