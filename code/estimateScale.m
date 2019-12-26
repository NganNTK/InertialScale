function [scale,gravity,bias] = ...
    estimateScale(A,b,scale0,gravity0,bias0,t)
% Estimation is performed in the frequency domain

fprintf('%s', repmat('-', 1, 60));
fprintf('\nFinal estimation in the frequency domain\n');
tic;

% Select valid range of frequencies [0Hz - 1.2Hz]
N = length(t);
fmax = 1.2;
fs = 1/mean(diff(t));
f = fs*(0:(N/2))/N;
freqRange = (f <= fmax);
fprintf('Upper limit for the frequencies = %.2f Hz\n',fmax);

% Optimize while enforcing gravity constraint: norm(g) = 9.82
% options = optimoptions(@fmincon, 'Display', 'off'); % original
options = optimoptions(@fmincon, 'Display', 'final-detailed');
% options = optimoptions(@fmincon, 'Display', 'iter'); % for debug
gConst = @gravityConstraint;

x0 = [scale0; gravity0; bias0];
% ----original----
% x = fmincon(@(x)minFunc(x, A, b, freqRange), ...
%     x0, [],[],[],[],[],[],gConst,options);
% ----original----
[x, fval, exitflag, output] = fmincon(@(x)minFunc(x, A, b, freqRange), ...
    	x0, [],[],[],[],[],[],gConst,options);
fprintf('Exit flag: %d\n', exitflag);   

% checking if fmincon is converged
if exitflag == 1
     disp('Optimization done!');
elseif exitflag == 2
    disp('Optimization not done! Do it again!');
    x02 = rand(size(x));
    [x2, fval2, exitflag2, output2] = fmincon(@(x)minFunc(x, A, b, freqRange), ...
    	x02, [],[],[],[],[],[],gConst,options);
    fprintf('Exit flag 2: %d\n', exitflag2);   
    if abs(x2(1) - x(1)) < 1e-10
        disp('Good');
    else
        fprintf('Error, diff_scale_x_x2 = %f\n', abs(x2(1) - x(1)));
    end
    fprintf('%s', repmat('-', 1, 60));
    fprintf('\nFinal estimates (second time)\n');
    fprintf('scale2 = %.4f ', x2(1));
    fprintf('\n');
end
% end


scale = x(1);
gravity = x(2:4);
bias = x(5:7);

fprintf('Finished in %.3f seconds\n', toc);

end


function [c,ceq] = gravityConstraint(x)

% c = [];
c = 1e-4 - x(1); % add constraint for scale, scale shouldn't be negative
ceq = norm([x(2) x(3) x(4)])-9.80;

end


function f = minFunc(x, A, b, freqRange)

Av = A*x; % Visual accelerations
Ai = b;    % Inertial accelerations

Av = [Av(1:3:end) Av(2:3:end) Av(3:3:end)];
Ai = [Ai(1:3:end) Ai(2:3:end) Ai(3:3:end)];

Fv = abs(fft(Av));
Fi = abs(fft(Ai));

Fv = Fv(freqRange,:);
Fi = Fi(freqRange,:);

f = (Fv - Fi).^2;
f = sum(f(:));

end