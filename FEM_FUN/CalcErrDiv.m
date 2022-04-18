%%%%%%%%%%%%%%%%%%%
%
% This files contains the command for computing the L2 error
% for the pressure and velocity and the H1 error for the veolcity
%
% Note: If utrue, Gradutrue and ptrue are all set to zero this
% routine computes the appropriate norms for the pressure and velocity.

clear temp
GradvelErr = zeros(2,2) ;
velErr = zeros(2,1) ;
preErr = 0.0 ;
preErr2=0;
divL2 = 0.0 ;
divTotal=0.0;

errordata=[];

for itrg=1:size(elnode,1)
      
   % Description of triangle.
   triag_no = itrg ;
	cotri(1:3,1) = nodeco(elnode(triag_no, 1:3), 1) ;
	cotri(1:3,2) = nodeco(elnode(triag_no, 1:3), 2) ;
    
	Jmat = [(cotri(2,1) - cotri(1,1)), (cotri(3,1) - cotri(1,1)) ; ...
        	(cotri(2,2) - cotri(1,2)) , (cotri(3,2) - cotri(1,2)) ] ;
	detJ = abs(Jmat(1,1)*Jmat(2,2) - Jmat(1,2)*Jmat(2,1));
	JInv = inv(Jmat) ;


	% Evaluation of quadrature points and quadrature weights.
	[quad_pts, quad_wghts] = feval('quad_75') ;
	nqpts = size(quad_pts,1) ;

	% Adjust points and weights to account for size of true triangle.
	xy_pts = ( Jmat * quad_pts.' ).' ;
	xy_pts(:,1) = cotri(1,1) + xy_pts(:,1) ;
	xy_pts(:,2) = cotri(1,2) + xy_pts(:,2) ;
	quad_wghts = detJ * quad_wghts ;

	% Evaluate the velocity, Gradient of the velocity and pressure functions 
   % at the quadrature points. 
   Gradv_vals = GradvelFun(xy_pts, triag_no) ;
   vfun_vals = velFun(xy_pts, triag_no) ;
   pfun_vals = preFun(xy_pts, triag_no) ;
   
%    Gradvtru_vals = Gradvtrue(xy_pts, triag_no) ;
%    vtru_vals = vtrue(xy_pts, triag_no) ;
%    ptru_vals = ptrue(xy_pts, triag_no) ;
%    
%    preErr = preErr + quad_wghts * ( (ptru_vals - pfun_vals).^2 ).' ;
%    velErr(1) = velErr(1) + quad_wghts * ( (vtru_vals(1,:) - vfun_vals(1,:) ).^2 ).'  ;
%    velErr(2) = velErr(2) + quad_wghts * ( (vtru_vals(2,:) - vfun_vals(2,:) ).^2 ).'  ;
%    temp(1:nqpts) = (Gradvtru_vals(1,1,:) - Gradv_vals(1,1,:) ).^2 ;
%    GradvelErr(1,1) = GradvelErr(1,1) +  quad_wghts * temp.'  ;
%    temp(1:nqpts) = (Gradvtru_vals(1,2,:) - Gradv_vals(1,2,:) ).^2 ;
%    GradvelErr(1,2) = GradvelErr(1,2) + quad_wghts * temp.'  ;
%    temp(1:nqpts) = (Gradvtru_vals(2,1,:) - Gradv_vals(2,1,:) ).^2 ;
%    GradvelErr(2,1) = GradvelErr(2,1) + quad_wghts * temp.'  ;
%    temp(1:nqpts) = (Gradvtru_vals(2,2,:) - Gradv_vals(2,2,:) ).^2 ;
%    GradvelErr(2,2) = GradvelErr(2,2) + quad_wghts * temp.'  ;

   temp(1:nqpts) = (Gradv_vals(1,1,:) + Gradv_vals(2,2,:)).^2 ;
   divL2 = divL2 + quad_wghts * temp.' ;
   
   temp2(1:nqpts) = (Gradv_vals(1,1,:) + Gradv_vals(2,2,:)) ;
   divTotal = divTotal + quad_wghts * temp2.';
   
%   preErr2 = preErr2 + quad_wghts * ( (ptru_vals - pfun_vals ).^2 ).' ;
%   preErr2 = preErr2 + quad_wghts * ( (ptru_vals - pfun_vals + gamma*temp2).^2 ).' ;
   
%   itrg = TrgMeNxt(itrg) ;  
   
%    errordata=[errordata;xy_pts(:,1), xy_pts(:,2), vfun_vals(1,:)' - vtru_vals(1,:)',...
%         vfun_vals(2,:)' - vtru_vals(2,:)', abs(pfun_vals'- ptru_vals'), sqrt(temp)' ];
%    

end

%%%
VelL2Error = sqrt( velErr(1) + velErr(2) );
GradVelL2Error = sqrt( GradvelErr(1,1) + GradvelErr(1,2) + ...
   GradvelErr(2,1) + GradvelErr(2,2) );

VelH1Error = sqrt(GradVelL2Error^2 + VelL2Error^2) ;

PreL2Error = sqrt(preErr) ;
   
divL2 = sqrt(divL2);

ModPreL2Error = sqrt(preErr2);

divTotal;

%display([ num2str(VelL2Error) ', ' num2str(GradVelL2Error) ', ' num2str(divL2) ])

% [X,Y]=meshgrid(0:.01:1,0:.01:1);
% velerror = griddata(errordata(:,1),errordata(:,2),sqrt(errordata(:,3).^2 + errordata(:,4).^2 ),X,Y);
% contourf(X,Y,velerror)
% 
% figure
% preserror = griddata(errordata(:,1),errordata(:,2),errordata(:,5)  ,X,Y);
% contourf(X,Y,preserror)
% 
% figure
% diverror = griddata(errordata(:,1),errordata(:,2),errordata(:,6)  ,X,Y);
% contourf(X,Y,diverror)
