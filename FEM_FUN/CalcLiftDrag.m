%%%%%%%%%%%%%%%%%%%
%
% This files contains the command for computing the L2 error
% for the pressure and velocity and the H1 error for the veolcity
%
% Note: If utrue, Gradutrue and ptrue are all set to zero this
% routine computes the appropriate norms for the pressure and velocity.

clear temp
cd = 0;
cl = 0;
divL2 = 0;
energy=0;
tic

for itrg=1:size(elnode,1);
    
    
   % Description of triangle.
   triag_no = itrg ;
% 	cotri(1:3,1) = nodeco(elnode(triag_no, 1:3), 1) ;
% 	cotri(1:3,2) = nodeco(elnode(triag_no, 1:3), 2) ;
%     
% 	Jmat = [(cotri(2,1) - cotri(1,1)), (cotri(3,1) - cotri(1,1)) ; ...
%         	(cotri(2,2) - cotri(1,2)) , (cotri(3,2) - cotri(1,2)) ] ;
%	detJ = abs(Jmat(1,1)*Jmat(2,2) - Jmat(1,2)*Jmat(2,1));
    detJ = detJglobal(itrg);

% 	% Evaluation of quadrature points and quadrature weights.
% 	[quad_pts, quad_wghts] = feval('quad_75') ;
% 	nqpts = size(quad_pts,1) ;
% 
% 	quad_wghts = detJ * quad_wghts ;

    quad_wghts = quad_wghtsglobal(triag_no,:);
    
    basvals = basisValues(:,:,triag_no);
    Gradtrue = gradBasisValues(:,:,:,triag_no);
    Gradx(:,:) = Gradtrue(:,1,:) ;
    Grady(:,:) = Gradtrue(:,2,:) ;


   
   
    Vstart = [2*(elnode(triag_no,1:3) - 1) + 1 , 2*(elnode(triag_no,4:6) + nVert - 1) + 1 ] ;
    Vel1 = GlobalV(Vstart, 1) ;
    VD1 = CDvec(Vstart, 1) ;
    VL1 = CLvec(Vstart, 1) ;
    VP1 = AllV(Vstart,end);    
    
    Vstart = Vstart + 1 ;
    Vel2 = GlobalV([Vstart], 1) ;
    VD2 = CDvec([Vstart], 1) ;
    VL2 = CLvec([Vstart], 1) ;
    VP2 = AllV(Vstart,end);
    
    
    Gradv_vals(1,1,:) = Vel1.' * Gradx ;
    Gradv_vals(1,2,:) = Vel2.' * Gradx ;
    Gradv_vals(2,1,:) = Vel1.' * Grady ;
    Gradv_vals(2,2,:) = Vel2.' * Grady ;
    
    vfun_vals(1,:) = Vel1.' * basvals ;
    vfun_vals(2,:) = Vel2.' * basvals ;
    
    vdrag(1,:) = VD1.' * basvals ;
    vdrag(2,:) = VD2.' * basvals ;
    
    vlift(1,:) = VL1.' * basvals ;
    vlift(2,:) = VL2.' * basvals ;
    
    vfun_valsP(1,:) = VP1.' * basvals ;
    vfun_valsP(2,:) = VP2.' * basvals ;   
    
    Gradvlift(1,1,:) = VL1.' * Gradx ;
    Gradvlift(1,2,:) = VL2.' * Gradx ;
    Gradvlift(2,1,:) = VL1.' * Grady ;
    Gradvlift(2,2,:) = VL2.' * Grady ;
    
    Gradvdrag(1,1,:) = VD1.' * Gradx ;
    Gradvdrag(1,2,:) = VD2.' * Gradx ;
    Gradvdrag(2,1,:) = VD1.' * Grady ;
    Gradvdrag(2,2,:) = VD2.' * Grady ;
    
    
    
    dvdt = (vfun_vals - vfun_valsP)/dt;
   
    avec1x(1:nqpts) = Gradv_vals(1,1,:) ;
    avec1y(1:nqpts) = Gradv_vals(2,1,:) ;
    avec2x(1:nqpts) = Gradv_vals(1,2,:) ;
    avec2y(1:nqpts) = Gradv_vals(2,2,:) ;
    
    NL_vals(1,:) = vfun_vals(1,:) .* avec1x + vfun_vals(2,:) .* avec1y ;
    NL_vals(2,:) = vfun_vals(1,:) .* avec2x + vfun_vals(2,:) .* avec2y ;
   
   term1D = quad_wghts * (dvdt(1,:) .* vdrag(1,:) + dvdt(2,:) .* vdrag(2,:) )' ; 
   term1L = quad_wghts * (dvdt(1,:) .* vlift(1,:) + dvdt(2,:) .* vlift(2,:) )' ;
   
   term2D = quad_wghts * (NL_vals(1,:).* vdrag(1,:) + NL_vals(2,:).* vdrag(2,:) )';
   term2L = quad_wghts * (NL_vals(1,:).* vlift(1,:) + NL_vals(2,:).* vlift(2,:) )';
   
   
   temp(1:nqpts) = ( Gradv_vals(1,1,:).*Gradvlift(1,1,:) + Gradv_vals(1,2,:).*Gradvlift(1,2,:) + Gradv_vals(2,1,:).*Gradvlift(2,1,:) + Gradv_vals(2,2,:).*Gradvlift(2,2,:) ) ;
   term3L = quad_wghts * temp.'  ;

   temp(1:nqpts) = ( Gradv_vals(1,1,:).*Gradvdrag(1,1,:) + Gradv_vals(1,2,:).*Gradvdrag(1,2,:) + Gradv_vals(2,1,:).*Gradvdrag(2,1,:) + Gradv_vals(2,2,:).*Gradvdrag(2,2,:)) ;
   term3D = quad_wghts * temp.'  ;


   cl = cl -20* (term1L + term2L +nu* term3L);
   cd = cd -20* (term1D + term2D +nu* term3D);
   
   temp(1:nqpts) = (Gradv_vals(1,1,:) + Gradv_vals(2,2,:)).^2 ;
   divL2 = divL2 + quad_wghts * temp.' ;
   
   energy = energy + quad_wghts * ( ( vfun_vals(1,:) ).^2 ).' + quad_wghts * ( ( vfun_vals(2,:) ).^2 ).'  ;
   
   
end
divL2 = sqrt(divL2);
energy = 1/2 * sqrt(energy);
toc
