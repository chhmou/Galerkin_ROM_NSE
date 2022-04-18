%%%%%%%%%%%%%%%%%%%
%
% This files contains the command for computing the L2 error
% for the pressure and velocity and the H1 error for the veolcity
%
% Note: If utrue, Gradutrue and ptrue are all set to zero this
% routine computes the appropriate norms for the pressure and velocity.

GradvelErr = zeros(2,2) ;
velErr = zeros(2,1) ;
preErr = 0.0 ;
for iDom = 1:nDom
   itrg = DomPoint(iDom) ;
   while itrg ~= -1
      
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
   
   Gradvtru_vals = Gradvtrue(xy_pts, triag_no,t) ;
   vtru_vals = vtrue(xy_pts, triag_no,t) ;
   ptru_vals = ptrue(xy_pts, triag_no,t-dt/2) ;
   
   preErr = preErr + quad_wghts * ( (ptru_vals - pfun_vals).^2 ).' ;
   velErr(1) = velErr(1) + quad_wghts * ( (vtru_vals(1,:) - vfun_vals(1,:) ).^2 ).'  ;
   velErr(2) = velErr(2) + quad_wghts * ( (vtru_vals(2,:) - vfun_vals(2,:) ).^2 ).'  ;
   temp(1:nqpts) = (Gradvtru_vals(1,1,:) - Gradv_vals(1,1,:) ).^2 ;
   GradvelErr(1,1) = GradvelErr(1,1) +  quad_wghts * temp.'  ;
   temp(1:nqpts) = (Gradvtru_vals(1,2,:) - Gradv_vals(1,2,:) ).^2 ;
   GradvelErr(1,2) = GradvelErr(1,2) + quad_wghts * temp.'  ;
   temp(1:nqpts) = (Gradvtru_vals(2,1,:) - Gradv_vals(2,1,:) ).^2 ;
   GradvelErr(2,1) = GradvelErr(2,1) + quad_wghts * temp.'  ;
   temp(1:nqpts) = (Gradvtru_vals(2,2,:) - Gradv_vals(2,2,:) ).^2 ;
   GradvelErr(2,2) = GradvelErr(2,2) + quad_wghts * temp.'  ;
   
   itrg = TrgMeNxt(itrg) ;  
  end

end

%%%
VelL2Error = sqrt( velErr(1) + velErr(2) );
GradVelL2Error = sqrt( GradvelErr(1,1) + GradvelErr(1,2) + ...
   GradvelErr(2,1) + GradvelErr(2,2) );

VelH1Error = sqrt(GradVelL2Error^2 + VelL2Error^2) ;

PreL2Error = sqrt(preErr) ;
   

