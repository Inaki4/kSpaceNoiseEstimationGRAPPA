function myfilt=WalshFilter(Rs,Rn)
    [U,S] = svd(inv(Rn)*Rs);
    myfilt = -U(:,1); 
end