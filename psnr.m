function y=psnr(f,g)
[m,n]=size(f);
y=-20*log10(norm(f-g,'fro')/sqrt(m*n)/255);