function z = findValueMid(x1,x2,x3,v1,v2,v3,flux)

% initial guess

z1 = v2-.1;
aa = [1 x1 x1^2; 1 x2 x2^2 ; 1 x3 x3^2] \ [v1; z1; v3];
myfun3 = @(x) aa(1) + aa(2)*x + aa(3)*x.^2;
err1 = flux - quad(myfun3,x1,x3);

z2 = v2+.1;
aa = [1 x1 x1^2; 1 x2 x2^2 ; 1 x3 x3^2] \ [v1; z2; v3];
myfun3 = @(x) aa(1) + aa(2)*x + aa(3)*x.^2;
err2 = flux - quad(myfun3,x1,x3);

%if err1*err2<0
%    display('Found a sign change, applying bisection')
%end

while (z2-z1)>1e-15
    m = (z1+z2)/2;
    aa = [1 x1 x1^2; 1 x2 x2^2 ; 1 x3 x3^2] \ [v1; m; v3];
    myfun3 = @(x) aa(1) + aa(2)*x + aa(3)*x.^2;
    errm = flux - quad(myfun3,x1,x3);
    
    if err1*errm<0
        z2=m;
        err2=errm;
    else
        z1=m;
        err1=errm;
    end
    [z1 z2];
end
z = (z1+z2)/2;



