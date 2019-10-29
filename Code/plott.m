code = train('train_new/', 4);
a = code{2};
b = code{1};
xa = a(1, :);
xb = b(1, :);
ya = a(5, :);
yb = b(5, :);

[s, fs] = audioread('train_new/s2.wav');
v1 = mfcc1(s, fs);

[s1, fs1] = audioread('train_new/s1.wav');
v2 = mfcc1(s1, fs1);

yv1 = v1(5, :);
xv1 = v1(1, :);

yv2 = v2(5, :);
xv2 = v2(1, :);

scatter(ya,xa,100,'r', 'filled')
hold on
scatter(yb, xb, 100, 'b', 'filled')
hold on 
scatter(yv1, xv1, 'r')
hold on 
scatter(yv2, xv2, 'b')

title('Clusters = 4')

