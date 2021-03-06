global model
% Plot stumbler outcome
figure(1)
set(1, 'units', 'normalized', 'position', [0.1 0.1 0.8 0.8])
h = axes;

% Lengths
Lstance     = model.Lstance;    % [m]
Lhip        = model.Lhip;       % [m]
Lthigh      = model.Lthigh;     % [m]
Lshank      = model.Lshank;     % [m]
Lfoot       = model.Lfoot;      % [m]

% Centres of gravity with respect to proximal joint
cgStance    = model.cgStance;   % [m]
cgThigh     = model.cgThigh;    % [m]
cgShank     = model.cgShank;    % [m]
cgFoot      = model.cgFoot;     % [m]

% % Part B: Brick coordinates
 bx1         = model.bx1;
 bx2         = model.bx2;
 by1         = model.by1;
 by2         = model.by2;
 bz1         = model.bz1;
 bz2         = model.bz2;

for jj = [3, 1]
    for kk = 1:jj:length(X_out.time)
        % Get angles from the correct index of X_out.signals.values
        % Example gamma1 = X_out.signals.values(kk,1);
        gamma1 = X_out.signals.values(kk,1);
        alpha2 =X_out.signals.values(kk,2);
        beta2= X_out.signals.values(kk,3);
        gamma2 =X_out.signals.values(kk,4);
        gamma3= X_out.signals.values(kk,5);
        gamma4=X_out.signals.values(kk,6);
        
        
        symb_Ti;
        % Mass locations (Hint; use values from Ti)
        %xmass = Ti([2]);
        %ymass = Ti([4]);
        %zmass = Ti([6]);
        xmass = [Ti(4);Ti(10);Ti(16);Ti(22)];
        ymass = [Ti(5);Ti(11);Ti(17);Ti(23)];
        zmass = [Ti(6);Ti(12);Ti(18);Ti(24)];
        
        
        
        
        
        % Joint locations (Hint; use values from Ti)
    
        xjoint = [0; [Ti(1);Ti(7);Ti(13);Ti(19);Ti(25)]];
        yjoint = [0; [Ti(2);Ti(8);Ti(14);Ti(20);Ti(26)]];
        zjoint = [0; [Ti(3);Ti(9);Ti(15);Ti(21);Ti(27)]];

        subplot(223) % Front view
        plot([0; Ti(3:3:end)], [0; Ti(2:3:end)], 'b-', 'linewidth', 2); hold on
        plot(zjoint, yjoint, 'bo', 'markerfacecolor', 'w');
        plot(zmass, ymass, 'ks', 'markerfacecolor', 'k');
        plot([-1 1], [0 0], 'k--'); hold off
        axis([-3 3 -3 3])
        title('Front view')

        subplot(221) % Top view
        plot([0; Ti(3:3:end)], [0; Ti(1:3:end)], 'b-', 'linewidth', 2); hold on
        plot(zjoint, xjoint, 'bo', 'markerfacecolor', 'w');
        plot(zmass, xmass, 'ks', 'markerfacecolor', 'k'); hold off
        axis([-3 3 -1 1]);
        title('Top view')

        subplot(224) % Right view
        plot([0; Ti(1:3:end)],  [0; Ti(2:3:end)], 'b-', 'linewidth', 2); hold on
        
        plot([bx1; bx2;bx2;bx1;bx1],  [by1; by1;by2;by2;by1], 'r-', 'linewidth', 2); hold on
        plot(xjoint, yjoint, 'bo', 'markerfacecolor', 'w');
        plot(xmass, ymass, 'ks', 'markerfacecolor', 'k');
        plot([-1 1], [0 0], 'k--'); hold off
        axis([-3 3 -3 3])
        text(-0.9, 0.9, ['TIME: ', num2str(X_out.time(kk)), ' s.'])
        title('Right view')

        subplot(2,2,2) % 3D view
        plot3([0; Ti(3:3:end)], [0; Ti(1:3:end)], [0; Ti(2:3:end)], 'b-', 'linewidth', 2); hold on
        plot3(zjoint, xjoint, yjoint, 'bo', 'markerfacecolor', 'w');
        plot3(zmass, xmass, ymass, 'ks', 'markerfacecolor', 'k'); hold off
        axis([-3 3 -3 3 -3 3])
        view([1 0.5 0.5])
        xlabel('z'); ylabel('x'); zlabel('y')
        grid on
        title('3D view')

        drawnow
    end
    pause(0.5)
end
disp('END of animation')