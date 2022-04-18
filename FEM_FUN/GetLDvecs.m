% create v_drag and v_lift

bdrydof2=[];
for ii=1:size(nodeco,1)
    x=nodeco(ii,1); y=nodeco(ii,2);
    if x< 0.000001 ||  y <0.000001 || y> .40999999 || abs( (x-.2)^2 + (y-.2)^2 ) < (0.05^2 + 0.0000001) || y>2.199999
        bdrydof2=[bdrydof2;2*ii-1;2*ii];
    end
end

% vdrag
% need to be 0 at all boundary, and <1,0> on cylinder
RHSvec = zeros(NVU+NPU,1) ;
Acoeff =  [StiffnessMatrix + 0*GradDivMatrix, PressureMatrix; PressureMatrix', sparse(NPU,NPU) ];

xypts2 = nodeco( bdrydof2(2:2:end) / 2 , :);
xpts = xypts2(:,1).';
ypts = xypts2(:,2).';
npts = size(xpts,2) ;

bdryvals = 0*bdrydof2;
for ii=1:npts
    x = xpts(ii);
    y = ypts(ii);
    if abs( (x-.2)^2 + (y-.2)^2 ) < (0.05^2 + 0.0000001)
        bdryvals(2*ii-1)=1;
        bdryvals(2*ii)=0;
    end
end

% Apply Dirichlet boundary conditions for velocity
bdrydirvals = bdryvals;
RHSvec = RHSvec - Acoeff(:,bdrydof2) * bdrydirvals;
Acoeff(:,bdrydof2) = 1e-30*Acoeff(:,bdrydof2);
Acoeff(bdrydof2,:) = 1e-30*Acoeff(bdrydof2,:);
for ii=1:size(bdrydof2,1)
        Acoeff(bdrydof2(ii),bdrydof2(ii))=1;
end
RHSvec(bdrydof2)=bdrydirvals;

x = Acoeff \ RHSvec;

CDvec = x(1:NVU);


%%%%%%%
% vlift
% need to be 0 at all boundary, and <0,1> on cylinder
RHSvec = zeros(NVU+NPU,1) ;
Acoeff =  [StiffnessMatrix + 0 * GradDivMatrix, PressureMatrix; PressureMatrix', sparse(NPU,NPU) ];


bdryvals = 0*bdrydof2;
for ii=1:npts
    x = xpts(ii);
    y = ypts(ii);
    if abs( (x-.2)^2 + (y-.2)^2 ) < (0.05^2 + 0.0000001)
        bdryvals(2*ii-1)=0;
        bdryvals(2*ii)=1;
    end
end

% Apply Dirichlet boundary conditions for velocity
bdrydirvals = bdryvals;
RHSvec = RHSvec - Acoeff(:,bdrydof2) * bdrydirvals;
Acoeff(:,bdrydof2) = 1e-30*Acoeff(:,bdrydof2);
Acoeff(bdrydof2,:) = 1e-30*Acoeff(bdrydof2,:);
for ii=1:size(bdrydof2,1)
        Acoeff(bdrydof2(ii),bdrydof2(ii))=1;
end
RHSvec(bdrydof2)=bdrydirvals;

x = Acoeff \ RHSvec;

CLvec = x(1:NVU);