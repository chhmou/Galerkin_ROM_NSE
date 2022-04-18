% This routine filters and deconvolves FilterAndDeconvolveThis and puts the
% solution back into it

% You can have either the Stokes or Helmholtz filters

if strcmp(filtertype,'Helmholtz')
    % (ubar,v) + a^2(grad ubar,grad v) = (u,v) forall v in Xh
    FilterMat = alpha^2 * StiffnessMatrix + MassMatrix + gd_parameter * GradDivMatrix;
    RHSvec = MassMatrix * FilterAndDeconvolveThis;
else
    % The Vh filter
    FilterMat = [alpha^2 * StiffnessMatrix + MassMatrix, PressureMatrix; -PressureMatrix', spalloc(NPU,NPU,10) ];
    RHSvec = [MassMatrix * FilterAndDeconvolveThis; zeros(NPU,1)];
end

Acoeff = FilterMat;

% Apply Dirichlet boundary conditions for velocity
% the dof numbers are in bdrydof, and we use the same bc for ubar as for u
bdrydirvals = FilterAndDeconvolveThis(bdrydof);
RHSvec = RHSvec - Acoeff(:,bdrydof) * bdrydirvals;
Acoeff(:,bdrydof) = 1e-16*Acoeff(:,bdrydof);
Acoeff(bdrydof,:) = 1e-16*Acoeff(bdrydof,:);
for ii=1:size(bdrydof,1)
        Acoeff(bdrydof(ii),bdrydof(ii))=1;
end
RHSvec(bdrydof)=bdrydirvals;

if strcmp(filtertype,'Stokes')
    Acoeff(:,end) = 1e-16 * Acoeff(:,end);
    Acoeff(end,:) = 1e-16 * Acoeff(end,:);
    Acoeff(end,end)=1;
    RHSvec(end)=0;
end

% Try to use a different solver here! cholinc + pcg?
xx = Acoeff \ RHSvec;

    

FilteredSoln = xx(1:NVU);

if N==0
    FilterAndDeconvolveThis = FilteredSoln;

else 
    % N==1
    
    % Filter the filtered velocity
    if strcmp(filtertype,'Helmholtz')
        RHSvec = MassMatrix * FilteredSoln;
    else
        RHSvec = [MassMatrix * FilteredSoln; zeros(NPU,1)];
    end

    Acoeff = FilterMat;

    bdrydirvals = getBdryVals(bdrydof,t);
    RHSvec = RHSvec - Acoeff(:,bdrydof) * bdrydirvals;
    Acoeff(:,bdrydof) = 1e-16*Acoeff(:,bdrydof);
    Acoeff(bdrydof,:) = 1e-16*Acoeff(bdrydof,:);
    for ii=1:size(bdrydof,1)
        Acoeff(bdrydof(ii),bdrydof(ii))=1;
    end
    RHSvec(bdrydof)=bdrydirvals;
    
    if strcmp(filtertype,'Stokes')
        Acoeff(:,end) = 1e-16 * Acoeff(:,end);
        Acoeff(end,:) = 1e-16 * Acoeff(end,:);
        Acoeff(end,end)=1;
        RHSvec(end)=0;
    end

    xx = Acoeff \ RHSvec;

    FilteredFilteredSoln = xx(1:NVU);
    FilterAndDeconvolveThis = 2*FilteredSoln - FilteredFilteredSoln;
end



clear Acoeff RHSvec FilterMat FilteredSoln FilteredFilteredSoln


