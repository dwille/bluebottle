%clear all; close all; clc;
%% tetrad_analysis.m
% - for given simulation and parameters
% -- form all tetrads fitting requirements
% -- calculate
% --- tetrad shape characteristics over time
% --- tetrad size characteristics over time
% INPUTS:
% -- r0: 1D array of tetrad base lengths, as a multiple of the particle radius
% -- ts: Index of start time (TODO: change to actual time, then find index)
% -- te: Index of end time (TODO: change to actual time, then find index)
% -- tol: Position tolerance as a multiple of particle radius
function tetrad_analysis(r0, ts, te, tol)
load part_data.mat
load grid_data.mat

% Conver r0 and tol to multiples of radius
r0 = r0*dom.r;
tol = tol*dom.r; % TODO: tol should be a function of r0 -- maybe tol on angle?

% Sort out desired time
nInd = 1:length(time);
ind = find(time < ts | time > te);
nInd(ind) = [];
% Deal with incorrect time input
if (isempty(nInd) == 1)
  fprintf('ts = %f and te = %f\n', time(1), time(end));
  error('Desired time is not within the simulation time limits.');
end
time(ind) = [];
ts = nInd(1);
te = nInd(end);

% Track absolute position of particles
[X, Y, Z] = periodic_flip(Xp, Yp, Zp, dom.N, length(time), dom.xl, dom.yl, dom.zl);

fprintf('Looping... \n')
for rr = 1:length(r0)
  % find all tetrads that satisfy r0(rr) positioning within tol
  % TODO: Don't have explicit periodicity in funciton (i.e. z should be included if needed)
  T = form_tetrads(r0(rr), X(:,ts), Y(:,ts), Z(:,ts), dom, tol);
  % TODO: ensure that tetrads are unique
  if isequal(T, -ones(4,1))
    fprintf('\tNo tetrads found for r0 = %.2f\n', r0(rr))
    rcount(rr) = 0;
    r0(rr) = -1;
    continue;
  else
    fprintf('\tFound %d tetrads for r0 = %.2f\n', size(T,2), r0(rr))
    rcount(rr) = size(T,2);
  end

  % loop over all tetrads
  for tet = 1:size(T,2)
    % Pull part number and position vector and velocity vector
    p1.n = T(1,tet);
    p2.n = T(2,tet);
    p3.n = T(3,tet);
    p4.n = T(4,tet);
    p1.X = [X(p1.n, :); Y(p1.n, :); Z(p1.n, :)];
    p2.X = [X(p2.n, :); Y(p2.n, :); Z(p2.n, :)];
    p3.X = [X(p3.n, :); Y(p3.n, :); Z(p3.n, :)];
    p4.X = [X(p4.n, :); Y(p4.n, :); Z(p4.n, :)];
    p1.U = [Up(p1.n, :); Vp(p1.n, :); Wp(p1.n, :)];
    p2.U = [Up(p2.n, :); Vp(p2.n, :); Wp(p2.n, :)];
    p3.U = [Up(p3.n, :); Vp(p3.n, :); Wp(p3.n, :)];
    p4.U = [Up(p4.n, :); Vp(p4.n, :); Wp(p4.n, :)];

    % calculate reduced vectors
    rho_1 = (p2.X - p1.X)./sqrt(2);
    rho_2 = (2*p3.X - p2.X - p1.X)./sqrt(6);
    rho_3 = (3*p4.X - p3.X - p2.X - p1.X)./sqrt(12);

    W_1 = (p2.U - p1.U)./sqrt(2);
    W_2 = (2*p3.U - p2.U - p1.U)./sqrt(6);
    W_3 = (3*p4.U - p3.U - p2.U - p1.U)./sqrt(12);

    % loop over time
    for tt = 1:length(time)
      %plot3(p1.X(1,tt), p1.X(2,tt), p1.X(3,tt),'ok','MarkerSize', 5, 'MarkerFaceColor', 'k')
      %axis([dom.xs dom.xe dom.ys dom.ye dom.zs dom.ze])
      %xlabel('x')
      %ylabel('y')
      %zlabel('z')
      %hold on
      %plot3(p2.X(1,tt), p2.X(2,tt), p2.X(3,tt),'or','MarkerSize', 5, 'MarkerFaceColor', 'r')
      %plot3(p3.X(1,tt), p3.X(2,tt), p3.X(3,tt),'ob','MarkerSize', 5, 'MarkerFaceColor', 'b')
      %plot3(p4.X(1,tt), p4.X(2,tt), p4.X(3,tt),'og','MarkerSize', 5, 'MarkerFaceColor', 'g')
      %pts12 = [p1.X(:,tt)'; p2.X(:,tt)'];
      %pts13 = [p1.X(:,tt)'; p3.X(:,tt)'];
      %pts14 = [p1.X(:,tt)'; p4.X(:,tt)'];
      %pts23 = [p2.X(:,tt)'; p3.X(:,tt)'];
      %pts24 = [p2.X(:,tt)'; p4.X(:,tt)'];
      %pts34 = [p3.X(:,tt)'; p4.X(:,tt)'];
      %plot3(pts12(:,1), pts12(:,2), pts12(:,3), 'k-')
      %plot3(pts13(:,1), pts13(:,2), pts13(:,3), 'k-')
      %plot3(pts14(:,1), pts14(:,2), pts14(:,3), 'k-')
      %plot3(pts23(:,1), pts23(:,2), pts23(:,3), 'k-')
      %plot3(pts24(:,1), pts24(:,2), pts24(:,3), 'k-')
      %plot3(pts34(:,1), pts34(:,2), pts34(:,3), 'k-')
      %hold off
      %drawnow
      %pause(0.25)

      G = [rho_1(:,tt), rho_2(:,tt), rho_3(:,tt)];
      Gmom = G*transpose(G);
      eigVal = eig(Gmom);
      g1 = eigVal(1);
      g2 = eigVal(2);
      g3 = eigVal(3);
      %[eigVec, eigVal] = eigs(Gmom);
      %g1 = eigVal(1,1);
      %g2 = eigVal(2,2);
      %g3 = eigVal(3,3);
      %gv1 = eigVec(:,1);
      %gv2 = eigVec(:,2);
      %gv3 = eigVec(:,3);
      if abs(g1) < 1e-5
        g1 = 1e-5;
      end
      if (g1 < 0 | g2 < 0 | g3 < 0)
        G
        Gmom
        fprintf('g1 g2 g3 %f %f %f\n', g1, g2, g3);
        error('eigs less than zero')
      end

      W = [W_1(:,tt), W_2(:,tt), W_3(:,tt)];
      %K = 0.5*(W*transpose(G) + G*transpose(W));
      %[eigVec, eigVal] = eigs(K);
      %k1 = eigVal(1,1);
      %k2 = eigVal(2,2);
      %k3 = eigVal(3,3);
      %kv1 = eigVec(:,1);
      %kv2 = eigVec(:,2);
      %kv3 = eigVec(:,3);

      % Radius of gyration
      Rsq(tet,tt) = g1 + g2 + g3;

      % volume
      Vol(tet,tt) = (g1*g2*g3)^(1/2)/3;

      % shape factors
      I1(tet,tt) = g1/Rsq(tet,tt);
      I2(tet,tt) = g2/Rsq(tet,tt);
      I3(tet,tt) = g3/Rsq(tet,tt);

      % lambda factor
      Lambda(tet,tt) = ...
        Vol(tet,tt)^(2/3)./Rsq(tet,tt);

      % angle between g and k
      %theta1(rr,tet,tt) = acos(dot(g1, k1)/(norm(g1)*norm(k1)));
      %theta2(rr,tet,tt) = acos(dot(g2, k2)/(norm(g2)*norm(k2)));
      %theta3(rr,tet,tt) = acos(dot(g3, k3)/(norm(g3)*norm(k3)));

      % store kappa
      %kappa1(tet,tt) = k1;
      %kappa2(tet,tt) = k2;
      %kappa3(tet,tt) = k3;

      % Coarse grained velocity gradient tensor
      %M = G.^(-1) * W - eye(3)*trace(inv(G)*W)./3;
      %S = 0.5*(M + transpose(M));
      %O = 0.5*(M - transpose(M));
      % tensor dot product Oij Oij
      %s(tet,tt) = 2*sum(sum(O.*O));

      % Invariants
      %P(tet,tt) = -trace(M);
      %Q(tet,tt) = 0.5*(trace(M)^2 - trace(M*M));
      %R(tet,tt) = -det(M);

    end
  end

  % average over all tetrads
  avgRsq(rr,:) = mean(Rsq, 1);
  avgVol(rr,:) = mean(Vol, 1);
  avgI1(rr,:) = mean(I1, 1);
  avgI2(rr,:) = mean(I2, 1);
  avgI3(rr,:) = mean(I3, 1);
  avgLambda(rr,:) = mean(Lambda, 1);
  %avgTheta1(rr,:) = mean(theta1, 1);
  %avgTheta2(rr,:) = mean(theta2, 1);
  %avgTheta3(rr,:) = mean(theta3, 1);
  %avgK1(rr,:) = mean(kappa1, 1);
  %avgK2(rr,:) = mean(kappa2, 1);
  %avgK3(rr,:) = mean(kappa3, 1);
  % TODO: probability as a function of theta
  %theta1_percent(rr,:) = numel(theta1(theta1 < pi/6))/numel(theta1);
  %theta2_percent(rr,:) = numel(theta2(theta2 < pi/6))/numel(theta2);
  %theta3_percent(rr,:) = numel(theta3(theta3 < pi/6))/numel(theta3);
  % TODO: save M to get sym / anti sym parts
  %invariants.(['r0_' num2str(r0(rr))]).s = s;
  %invariants.(['r0_' num2str(r0(rr))]).P = P;
  %invariants.(['r0_' num2str(r0(rr))]).Q = Q;
  %invariants.(['r0_' num2str(r0(rr))]).R = R;
  %invariants.(['r0_' num2str(r0(rr))]).Rsq = Rsq;

  clearvars Rsq Vol I1 I2 I3 Lambda;
  %clearvars kappa1 kappa2 kappa3;
  %clearvars theta1 theta2 theta3; 
  %clearvars M P Q R s Rsq;
end
fprintf(' ... Done!\n')

% Save tetrad average values
save('tetrad_stats.mat', ... 
      'avgI1', 'avgI2', 'avgI3', 'avgLambda', 'avgRsq', 'avgVol', ...
      'r0', 'time', 'dom')
%      'avgTheta1', 'avgTheta2', 'avgTheta3', 'avgK1', 'avgK2', 'avgK3',...
%      'theta1_percent', 'theta2_percent', 'theta3_percent', ...
%      'P', 'Q', 'R',... 
%save('invariants.mat', 'invariants')
