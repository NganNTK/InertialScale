% checkInput
function checkVisualInput(posV, qtV, tV)
% function checkVisualInput(posV, qtV, tV)
% sanity check for posV, qtV, and tV
% Inputs: posV: Nx3 matrix of x,y,z camera centres
%         qtV:  Nx4 matrix of quaternions of camera rotations
%         tV:   Nx1 vector of times in nanoseconds
%
% posV: check for smooth trajectory
% qtV: check if valid rotation, and if rotations are changing smoothly
% tV: check if time sampling is uniform

[N,junk] = size(posV);

% check pos
posDiff = diff(posV); % row difference
posDiffAbs = sum(abs(posDiff),2);
posDiffMean = mean(posDiffAbs);
[posDiffMax, iPosDiffMax] = max(posDiffAbs);
% fprintf('ID max of camera trajectory = %d\n', iPosDiffMax);
fprintf('Mean, max of camera trajectory = %f, %f\n', posDiffMean, posDiffMax);

% check rotations
qtRot = qt_dircos(qtV'); % convert to 3x3 rotation matrices
for ii = 1:N  
    qtRotDet(ii) = det(qtRot(:,:,ii)); %compute determinant
end
ind = find(qtRotDet < 1); % any det < 1?
if length(ind) == 0
    fprintf('All rotation matrices have det 1\n');
else
    fprintf('fraction of suspicious rotation matrices = %.2f\n', length(ind)/N);
end
qtDiff = diff(qtV);
qtDiffAbs = sum(abs(qtDiff),2);
qtDiffMean = mean(qtDiffAbs);
qtDiffMax = max(qtDiffAbs);
fprintf('Mean, max of quarternions = %f, %f\n', qtDiffMean, qtDiffMax);

% check time
tDiff = diff(tV);
tDiffMean = mean(tDiff);
tDiffStd = std(tDiff);
fprintf('Mean, std dev of timings = %f, %f\n', tDiffMean, tDiffStd);

% plots
subplot(3,1,1); plot(posDiffAbs, 'b'); title('Smoothness of camera trajectory');
subplot(3,1,2); plot(qtDiffAbs, 'r'); title('Smoothness of quaternions');
subplot(3,1,3); plot(tDiff, 'k'); title('Timing intervals');

% plot camera trajectory
disp('Plotting camera trajectory');
figure;
plot3(posV(:,1), posV(:,2), posV(:,3));
xlabel('x');
ylabel('y');
zlabel('z');
title('Camera trajectory');
hold on;
grid on;
plot3(posV(1,1),posV(1,2),posV(1,3),'ro');
plot3(posV(N,1),posV(N,2),posV(N,3),'g*');

end


