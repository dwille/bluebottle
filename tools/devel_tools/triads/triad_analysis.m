%% triad_analysis.m
% - for given simulation and parameters
% -- form all triads fitting requirements
% -- calculate
% --- triad shape characteristics over time
% --- triad size characteristics over time
%
% USAGE:
% triad_analysis(r0, ts, dt, te, tol)
%
% INPUTS:
% -- r0: 1D array of triad base lengths, as a multiple of the particle radius
% -- ts: Start time
% -- dt: Sampling frequency of TIMESTEPS (in units of timesteps)
% -- te: End time
% -- tol: Position tolerance as a multiple of particle radius

function triad_analysis(r0, ts, dt, te, tol)
% TESTING:
%r0 = 4; ts = 500; dt = 10; te = 600; tol = 1;

load part_data.mat
load grid_data.mat

% Conver r0 and tol to multiples of radius
r0 = r0*dom.r;
tol = tol*dom.r;

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
[X, Y, Z] = periodic_flip(Xp, Yp, Zp, dom.N, length(time), ...
                          dom.xl, dom.yl, dom.zl);

fprintf('Looping... \n')
for rr = 1:length(r0)
  % find all triads that satisfy r0(rr) positioning within tol
  % TODO: Don't have explicit periodicity in function 
  % (i.e. z should be included if needed)
  T = form_triads(r0(rr), X(:,ts), Y(:,ts), Z(:,ts), dom, tol);
  if isequal(T, -ones(3,1))
    fprintf('\tNo triads found for r0 = %.2f\n', r0(rr))
    rcount(rr) = 0;
    r0(rr) = -1;
    continue;
  else
    fprintf('\tFound %d triads for r0 = %.2f\n', size(T,2), r0(rr))
    rcount(rr) = size(T,2);
  end

  % loop over all triads
  tetiter = 0;
  for tet = 1:size(T,2)
    tetiter = tetiter + 1;
    % Pull part number and position vector
    p1.n = T(1,tet);
    p2.n = T(2,tet);
    p3.n = T(3,tet);
    p1.X = [X(p1.n, :); Y(p1.n, :); Z(p1.n, :)];
    p2.X = [X(p2.n, :); Y(p2.n, :); Z(p2.n, :)];
    p3.X = [X(p3.n, :); Y(p3.n, :); Z(p3.n, :)];

    % calculate reduced vectors
    rho_1 = (p2.X - p1.X)./sqrt(2);
    rho_2 = (2*p3.X - p2.X - p1.X)./sqrt(6);
    rho_3 = [0;0;0;];

    % loop over time
    titer = 0;
    time = ts:10:te;
    for tt = ts:dt:te
      % Plot over time
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

      titer = titer + 1;
      G = [rho_1(:,tt), rho_2(:,tt), rho_3];
      Gmom = G*transpose(G);
      eigVal = eigs(Gmom);
      g1 = eigVal(1);
      g2 = eigVal(2);
      g3 = eigVal(3);
      if (g1 < 0 | g2 < 0)
        error('eigs less than zero')
      elseif (abs(g3) > 1e-10)
        error('g3 not zero')
      end

      % avg of squared lengths of each side
      Rsq(tetiter, titer) = g1 + g2;

      % area
      triA(tetiter, titer) = norm(cross(rho_1(:,tt), rho_2(:,tt)));

      % aspect ratio -- ratio of area to that of equilateral of same scale
      triW(tetiter, titer) = 4*triA(tetiter, titer)/(sqrt(3)*...
                               Rsq(tetiter, titer));

      % shape factors
      I1(tetiter, titer) = g1/Rsq(tetiter, titer);
      I2(tetiter, titer) = g2/Rsq(tetiter, titer);

      % chi factor
      chi(tetiter, titer) = 0.5*atan(2*dot(rho_1(:,tt), rho_2(:,tt))/...
                            (norm(rho_2(:,tt))^2 - norm(G*rho_1(:,tt))^2));
    end
  end

  avgRsq(rr,:) = mean(Rsq, 1);
  avgtriA(rr,:) = mean(triA, 1);
  avgtriW(rr,:) = mean(triW, 1);
  avgI1(rr,:) = mean(I1, 1);
  avgI2(rr,:) = mean(I2, 1);
  avgChi(rr,:) = mean(chi, 1);

  clearvars Rsq Vol I1 I2 Lambda;
end
fprintf(' ... Done!\n')

% Average values over all triads (dim=2)
time = ts:dt:te;
save('triad_stats.mat', 'avgI1', 'avgI2', 'avgChi', 'avgRsq', 'avgtriA',...
     'avgtriW', 'r0', 'time', 'dom')
