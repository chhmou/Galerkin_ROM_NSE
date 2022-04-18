% curl then filter and deconvolve!

% Still using ExtrapVel
ExtrapVel = 3/2*PastV(:,end) - 1/2*PastV(:,end-1); 

%% Filter and deconvolve the extrapolated curl 

% First we need to project the curl into Vh, so that we have a boundary
% condition for the filtering.
ProjMat = [MassMatrix, PressureMatrix; PressureMatrix', spalloc(NPU,NPU,10)];
RHSvec = zeros(NVU+NPU,1);
for iDom = 1:nDom
       itrg = DomPoint(iDom) ;
       while itrg ~= -1
          if strcmp(vel_bas_type, 'CtsQuad') == 1
                Vstart = [2*(elnode(itrg,1:3) - 1) + 1 , 2*(elnode(itrg,4:6) + nVert - 1) + 1];
                localUnk = reshape([Vstart ; Vstart+1],[12,1]) ;            
          end
          % (curl extrapvel,v)
          rhscurlextrapvel =  inner_prod_ten1_Vec(itrg,'quad_75','curlExtrap_2d', vel_bas_type) ;                    
        % Assemble the local matrices
        RHSvec(localUnk) = RHSvec(localUnk) + rhscurlextrapvel;
        itrg = TrgMeNxt(itrg) ;  
       end
end
% For the filter, we will use the same rhs vector
RHSFilter=RHSvec;
ProjMat(:,end)=1e-16*ProjMat(:,end);
ProjMat(end,:)=1e-16*ProjMat(end,:);
ProjMat(end,end)=1;
RHSvec(end)=0;
xx = ProjMat\RHSvec;
curlVelInVh = xx(1:NVU);

% Next we need to filter and deconvolve the curl of the extrapolated
% velocity, so filter and then filter again

% (phibar,v) + alpha^2 *(grad phibar,grad v) + (lambda,div v) = (curl extrapvel, v)
% (div phibar,q)=0
FilterMatrix = [MassMatrix + alpha^2 * StiffnessMatrix, PressureMatrix ; PressureMatrix', spalloc(NPU,NPU,100)];
% already have RHSFilter, so just need to set boundary condition
bdrydirvals = curlVelInVh(bdrydof);
RHSFilter = RHSFilter - FilterMatrix(:,bdrydof) * bdrydirvals;
FilterMatrix(:,bdrydof) = 1e-16*FilterMatrix(:,bdrydof);
FilterMatrix(bdrydof,:) = 1e-16*FilterMatrix(bdrydof,:);
for ii=1:size(bdrydof,1)
        FilterMatrix(bdrydof(ii),bdrydof(ii))=1;
end
RHSFilter(bdrydof)=bdrydirvals;
% Dirichlet "pressure" node
FilterMatrix(end,:)=1e-16*FilterMatrix(end,:);
FilterMatrix(:,end)=1e-16*FilterMatrix(:,end);
FilterMatrix(end,end)=1;
RHSFilter(end)=0;

xx = FilterMatrix \ RHSFilter;
FilteredCurl = xx(1:NVU);

% Now filter again
FilterMatrix = [MassMatrix + alpha^2 * StiffnessMatrix, PressureMatrix ; PressureMatrix', spalloc(NPU,NPU,100)];
RHSFilter = MassMatrix*FilteredCurl;
bdrydirvals = curlVelInVh(bdrydof);
RHSFilter = RHSFilter - FilterMatrix(:,bdrydof) * bdrydirvals;
FilterMatrix(:,bdrydof) = 1e-16*FilterMatrix(:,bdrydof);
FilterMatrix(bdrydof,:) = 1e-16*FilterMatrix(bdrydof,:);
for ii=1:size(bdrydof,1)
        FilterMatrix(bdrydof(ii),bdrydof(ii))=1;
end
RHSFilter(bdrydof)=bdrydirvals;
% Dirichlet "pressure" node
FilterMatrix(end,:)=1e-16*FilterMatrix(end,:);
FilterMatrix(:,end)=1e-16*FilterMatrix(:,end);
FilterMatrix(end,end)=1;
RHSFilter(end)=0;

xx = FilterMatrix \ RHSFilter;
FilteredCurl2 = xx(1:NVU);

omegadc_curl = 2*FilteredCurl - FilteredCurl2;