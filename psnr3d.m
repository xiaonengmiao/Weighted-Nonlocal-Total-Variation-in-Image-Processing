function y=psnr3d(f,g)
[m,n,k]=size(f);
x = 0;
for i = 1:k
    x = x + norm(f(:,:,i)-g(:,:,i),'fro');
end
y=-20*log10(x/(k*sqrt(m*n))/255);