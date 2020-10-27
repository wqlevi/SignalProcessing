function tsnr=tSNR(x)
% [Revision] The revision was drafted by @Qi.Wang
tsnr = squeeze(mean(x,2) ./ std(x,0,2));
% [ABORTED] Following equation originated from @S.C. Choi
% ave = mean(x,2);
% n = size(x,2);
% d = sqrt(sum((x(1,:)- ave).^2)./n);
% tsnr = ave./d;
end