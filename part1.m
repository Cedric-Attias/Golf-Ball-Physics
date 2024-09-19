shot1 = load('shot1.mat');
shot2 = load('shot2.mat');
shot3 = load('shot3.mat');
shot4 = load('shot4.mat');
shot5 = load('shot5.mat');
shot6 = load('shot6.mat');
shot7 = load('shot7.mat');
shot8 = load('shot8.mat');
shot9 = load('shot9.mat');
shot10 = load('shot10.mat');

carry_real = [252 286 277 274 283 268 184 183 279 264];
carry_sim = [235.252 250.322 261.452 257.391 255.508 243.589 173.233 175.806 233.511 251.366];
offline_real = [45 -32 -11 -12 -1 -13 -17 -15 -53 -1];
offline_sim = [61.151 -50.543 -15.293 -34.855 1.0581 -44.548 -12.970 -24.331 -80.103 2.693];

sums = zeros([1 10]);


for i=1:length(carry_real)
    carry_diff = (carry_real(i) - carry_sim(i))^2;
    offline_diff = (offline_real(i) - offline_sim(i))^2;
    diff = sqrt(carry_diff + offline_diff);
    sums(i) = diff;
end

LPE = sum(sums)/10;



% figure
% axis equal
% grid on
% hold on
% plot3(shot1.X, shot1.Z, shot1.Y)
% plot3(shot2.X, shot2.Z, shot2.Y)
% plot3(shot3.X, shot3.Z, shot3.Y)
% plot3(shot4.X, shot4.Z, shot4.Y)
% plot3(shot5.X, shot5.Z, shot5.Y)
% plot3(shot6.X, shot6.Z, shot6.Y)
% plot3(shot7.X, shot7.Z, shot7.Y)
% plot3(shot8.X, shot8.Z, shot8.Y)
% plot3(shot9.X, shot9.Z, shot9.Y)
% plot3(shot10.X, shot10.Z, shot10.Y)
% set(gca, 'YDir','reverse')
% axis([0 280 -90 90 0 70])
% xlabel('X (yds)')
% ylabel('Z (yds)')
% zlabel('Y (yds)')
% title('3D Trajectory')
% view([-65, 8])
