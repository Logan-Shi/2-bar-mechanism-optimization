function swing_time = swingTime(x,C1,C2)
max_time = 100;
dt = 1e-3;
max_cnt = max_time/dt;
l1 = x(1);
l2 = x(2);
lDE = 0.15;
lGF = 0.1;
h = 0.5;
m1 = 0.6;%kg
m2 = 0.5;%kg
J1 = 1/12*m1*lDE^2;%kgm^2
J2 = 1/12*m2*lGF^2;%kgm^2
g = 9.8;%m/s^2;
theta0 = [0,pi/2];

% initial condition
res = fkine2bar(theta0,l1,l2);
q = [res(1:2),theta0(1),res(3:end),theta0(2)];
dq = zeros(1,size(q,2));

% acceleration
M = diag([m1,m1,J1,m2,m2,J2]);
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

counter = 0;
while true
    theta1  = q(3);dtheta1  = dq(3);
    theta2 = q(6);dtheta2 = dq(6);
    
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
    ddq = ddq_temp(1:size(q,2))';
    
    q = q + dq*dt;
    dq = dq + ddq*dt;
    
    % choose theta as independent variable
    % position constraint
    res = fkine2bar([q(3),q(6)],l1,l2);
    q = [res(1:2),q(3),res(3:end),q(6)];
    % velocity constraint
    Phi_v = [Phi_q(:,3),Phi_q(:,6)];
    Phi_u = [Phi_q(:,1:2),Phi_q(:,4:5)];
    b = -Phi_v*[dq(3);dq(6)];
    res = (Phi_u\b)';
    dq = [res(1:2),dq(3),res(3:end),dq(6)];
    
    %     position_error(i+1,:) = norm(Phi([q(i+1,1:2),q(i+1,4:end)],q(i+1,3)));
    %     velocity_error(i+1,:) = norm(Phi_q*dq(i+1,:)');
    energy = dq*M*dq'/2+q(2)*m1*g+q(5)*m2*g;
    if energy<min_energy*1.001 || counter> max_cnt
        swing_time = -counter*dt;
        break
    end
    counter = counter + 1;
end

% plot(energy)
% hold on
% plot(repmat(min_energy,size(energy)))
% if finished
%     disp(['went on for ' num2str(t_end) ' sec(s)'])
% else
%     disp(['still going after ' num2str(t(end)) ' sec(s)'])
% end
end