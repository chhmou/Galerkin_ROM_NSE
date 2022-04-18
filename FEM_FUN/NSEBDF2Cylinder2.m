%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear
format long

%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge GlobalVnew
global GlobalV  GlobalP  nu vel_bas_type  pre_bas_type problem sig
global GlobalPSV sumV t AllV dt




vel_bas_type = 'CtsQuad' ;  
pre_bas_type = 'DCtsLin' ;  

 problem = 'Cylinder'
 nu=1/1000;
 t=0;
 T=10;
 dt=0.002;
% meshname = 'cylmesh35.msh';
% meshname = '12Kmesh.msh';
 meshname = 'cylmesh8K.msh';
 %[nodeco, elnode, bdynde] = FFEM2ErvGrd2(meshname) ;
GenSVgrd(meshname);
 [nodeco, elnode, bdynde] = FFEM2ErvGrd2('SVgrid.msh') ;

gamma=0.0; 


nbdy = size(bdynde,1) ;
bdynde = [bdynde , zeros(nbdy , 1), -1*ones(nbdy , 2) ] ;
[nodeco, elnode, bdynde, bdyedge, nVert, nedge] =  midEDGEgen(nodeco, elnode, bdynde) ;
 
 GlobalV = zeros((nVert + nedge)*2, 1) ;
 
 bdynde=[];
 for i=1:size(nodeco,1)
     x=nodeco(i,1);
     y=nodeco(i,2);
    if x< 0.000001 ||  y <0.000001 || y> .40999999 || abs( (x-.2)^2 + (y-.2)^2 ) < (0.05^2 + 0.0000001) || x>2.1999
        bdynde=[bdynde;i];
    end
 end
bdrysize = 2*size(bdynde,1);
bdrydof=[];
for ii=1:bdrysize/2
    bdrydof = [bdrydof; 2*bdynde(ii,1)-1; 2*bdynde(ii,1)];
end

if strcmp(pre_bas_type, 'CtsLin') == 1 
   GlobalP = zeros(nVert,2) ;

elseif strcmp(pre_bas_type, 'DCtsLin') == 1 
   ntri = size(elnode,1) ;
   GlobalP = zeros(3*ntri,2) ;
   GlobalP(:,2) = -1 ;
   mapN2GP = reshape(elnode(:,1:3).',[3*ntri 1]) ;
end
NVU = size(GlobalV,1) 
NPU = size(GlobalP,1) 
Sdim = NVU+NPU;

% create the matrices that will be involved in each timestep
createConstantMatrices2;
%load mats55KSV
%save constmats16bary MassMatrix StiffnessMatrix PressureMassMatrix PressureMatrix GradDivMatrix


% initialize u0
GlobalV = initVel(nodeco);
AllV = [GlobalV,GlobalV];

%  solve Stokes with Cylinder for initial condition
RHSvec = zeros(NVU+NPU,1) ;
Acoeff =  [1/100*StiffnessMatrix + gamma * GradDivMatrix,  PressureMatrix; PressureMatrix', spalloc(NPU,NPU,100*NPU) ];
% Apply Dirichlet boundary conditions for velocity
bdrydirvals = GlobalV(bdrydof);
RHSvec = RHSvec - Acoeff(:,bdrydof) * bdrydirvals;
Acoeff(:,bdrydof) = 1e-30*Acoeff(:,bdrydof);
Acoeff(bdrydof,:) = 1e-30*Acoeff(bdrydof,:);
for ii=1:size(bdrydof,1)
        Acoeff(bdrydof(ii),bdrydof(ii))=1;
end
RHSvec(bdrydof)=bdrydirvals;
% % Dirichlet pressure node
Acoeff(:,NVU+1)=1e-30 * Acoeff(:,NVU+1);
Acoeff(NVU+1,:)=1e-30 * Acoeff(NVU+1,:);
Acoeff(NVU+1,NVU+1)=1;
RHSvec(NVU+1)=0;

Acoeff = Acoeff .* (abs(Acoeff)>1e-16);

soln = Acoeff \ RHSvec;
GlobalV = soln(1:NVU);
GlobalP = soln(NVU+1:end);

CalcErrDiv
divL2


% calculate basis function values and store it
basisValues = zeros(6,7,size(elnode,1));
gradBasisValues = zeros(6,2,7,size(elnode,1));
for triag_no=1:size(elnode,1)
    % Description of triangle.
    cotri(1:3,1) = nodeco(elnode(triag_no, 1:3), 1) ;
    cotri(1:3,2) = nodeco(elnode(triag_no, 1:3), 2) ;
    
    Jmat = [(cotri(2,1) - cotri(1,1)), (cotri(3,1) - cotri(1,1)) ; ...
        (cotri(2,2) - cotri(1,2)) , (cotri(3,2) - cotri(1,2)) ] ;
    detJ = abs(Jmat(1,1)*Jmat(2,2) - Jmat(1,2)*Jmat(2,1));
    JInv = inv(Jmat) ;
    
    detJglobal(triag_no)=detJ;

    % Evaluation of quadrature points and quadrature weights.
    [quad_pts, quad_wghts] = feval('quad_75') ;
    nqpts = size(quad_pts,1) ;

    % Adjust points and weights to account for size of true triangle.
    xy_pts = ( Jmat * quad_pts.' ).' ;
    xy_pts(:,1) = cotri(1,1) + xy_pts(:,1) ;
    xy_pts(:,2) = cotri(1,2) + xy_pts(:,2) ;
    quad_wghts = detJ * quad_wghts ;
    
    quad_wghtsglobal(triag_no,:)=quad_wghts;

    %  Evaluate Basis Functions and their Gradients at quad. points.
    [ten1a, Gradten1a] = feval(vel_bas_type, quad_pts) ;
    basisValues(:,:,triag_no)=ten1a;
    
    for iq = 1:nqpts
        Gradtrue1(:,:,iq) = Gradten1a(:,:,iq) * JInv ;
    end

    gradBasisValues(:,:,:,triag_no)=Gradtrue1;
    
end
    

% setup for lift/drag
GetLDvecs

BalanceTable=[];

Snapshots=zeros(NVU,5*round(1/dt));
snapcount=0;
 
 for jj=1:round(T/dt)
     
     t = jj*dt;
     NSEsolve_BDF2_noYosida


    if jj>2          
        CalcLiftDrag
    
        BalanceTable=[BalanceTable; t, energy,divL2, cd, cl]   
    end
    
     if t>5 
         snapcount=snapcount+1;
         Snapshots(:,snapcount)=GlobalV;
     end
    
      AllV = [AllV(:,2:end),GlobalV(:,1)];    
     
 end


save snapshotData25Kdt002SV_Re100


