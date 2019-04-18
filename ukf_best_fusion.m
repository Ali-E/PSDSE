close	all 
clear	all
	clc
	tic
	instance_number=1;
	temp_dir = sprintf('%s/ temporary_files',pwd);
	if	exist(temp_dir)
		rmdir(temp_dir , 's')
	end
	mkdir(temp_dir)% Open a fresh	temporary	directory
	if	1 % Define	fault
		fault(1).nodes = [5	6];
		fault(1).time = [6.1	6.16];
		fault(1).location = 0.5; % fault	distance	(from nodes(1)) 	relative	to the	line	length
	else
		fault = [ ] ;
    	end
    	[ Sl_star , Vl , genbus_id , loadbus_id , xd , V_0 , V_0_A , Sg_star,
	NOB, D, H, M, Gn, ns , T_end , Ts]= indi_conds;

%       NOB=9;
% 	T_end=10;
% 	Ts=0.03;
% 	genbus_id =[1;2;3];
% 	Gn=length ( genbus_id ) ;
% 	ns=2*Gn-2;
% 	%nlbus id =[1;	2;	3;	4;	] ;
% 	%loadbus id=setdiff ([1:NOB] , nlbus_id ) ;
% 	loadbus_id =[5;	6;	8];
% 	V_0=[1.04;	1.025;	1.025];
% 
%  	V_0_A=pi /180*[0.00;	9.3;	4.7];
% 
% 	i=sqrt (-1) ;
% 
% 	Sl_star =[1.25-0.5*1 i ;	0.9-0.3*1 i ; 1-0.35*1 i ] ;
% 
% 
%  	Vl = [0.996;	1.013;	1.016];
% 
% 	Sg_star=[ 0.716 - 1 i .*0.27;	1.63 - 1 i *0.067;	0.85 + 1 i.*0.109];
% 
% 	xd=[0.0608;	0.1198;	0.1813];
% 	H=[23.64;	6.4;	3.01];
% 
% 	M = H./( pi *60) ;
% 	D = 5*M; %	LET US ASSUME UNIFORM DAMPING FOR MODEL REDUCTION
	%=====================================
	% Describe Fault
	%=====================================

	n=2; %Ratio of smapling	rates	of Measurements
	%clear	latest_n
	latest_n=n;
    	save	latest_n ;
	ns=(length(xd))*2-2;
	%syms x
	%x=[x1 ; x2 ; x3 ; x4 ] ;
 	x_simu = State_Simulator_9(fault);
	[ Y_big , Yr, Rv, E0,	delta0 , Pm, RefNode , x0 ] = Calc_Y_bus_9(fault);

	%======================================
	% Initialization	before EKF
	%======================================

	k_end=floor(T_end/Ts) ;
	rz2 =0.02^2; % Fast PMU noise
	rx =0.003^2; % Process	noise
	rz1 =0.025^2; % Slow PMI noise
	ks=1;
	[a ,b] = size(x_simu);
	x_est=zeros(ns , k_end);
	x_sim=zeros(ns , k_end);
	x_21=zeros(ns , k_end);
	x_21=zeros(ns , ns , k_end);
	x_est1=zeros(ns , round( k_end/n)+1) ;
	%x0=[delta0(1:2) ;	0;0 ];
	x_est1(: ,1)=x0 ;
	x_est(: ,1)=x0 ;
	x_sim(: ,1)=x0 ;
	x_21(: ,1)=x0 ;

	%======================================
	% Now the FOR loop
	%======================================
	% These	select	the buses where PMU groups are	placed !	sel_va	-sel_vb	is	for
	% finest	group PMU	and sel_vc-sel_vd	is	for	coarsest	group PMUs
	sel_va =1; 
    sel_vb =5;
	sel_vc =3;
	sel_vd =9;

	if	~isempty(fault)
 Y_big0 = Y_big ;  Yr0 = Yr;  Rv0 = Rv;

	% NOTE: WE ASSUME THAT THERE IS ONLY ONE FAULT!

	for	iF = 1:1%length(fault)

	k_Fstart = floor(fault(iF).time(1)/Ts);
	k_Fend = floor(fault(iF).time(2)/Ts);

	%===========================================
	% Some additional	initialization %-------------------------------------
	%===========================================

 eval(sprintf('[f,x] = dx_for_9bus_%d_preF_sym' ,iF)) ; 
 Yr = Yr0(iF).preF ; 
 Rv = Rv0(iF).preF ;
	Y_big = Y_big0(iF).preF ;
	dx = f ;
	gx = dx.*Ts + x ;
	gx1 = dx.*(n*Ts) + x ;
	%dfdx = jacobian(dx , x) ;

	E1 = sym(E0(1:Gn,1));
	E2 = sym(0*E1) ;
	E1( setdiff([1:Gn] , RefNode) ,1) = E0( setdiff([1:Gn] ,RefNode) ,1) .* cos(x(1:Gn-1,1)) ;
	E2( setdiff([1:Gn] , RefNode) ,1) = E0( setdiff([1:Gn] ,RefNode) ,1) .* sin(x(1:Gn-1,1)) ;
	%
	%
%==============================================================
	% Construction	of Measurement System
	%
%==============================================================
	Yr1 = real(Yr) ; 
    Yr2 = imag(Yr) ;
	I1 =(Yr1*E1-Yr2*E2) ;
	I2 =(Yr1*E2+Yr2*E1) ;

 Rv1 = real(Rv) ; 
 Rv2 = imag(Rv) ;
	V1 =(Rv1*E1-Rv2*E2) ;
	V2 =(Rv1*E2+Rv2*E1) ;

	%=====================
	% 1) PMU-1 Measurement(COARSEST)
	%=====================

	mz1 =2*abs(( sel_vc - sel_vd))+2;
	% H1f = jacobian(z1f , x) ;
	z1 = zeros(mz1, k_end); % Creation	of	the array	for	storing	the simulated	values	of z1
	z1(: ,1) = compute_z2fsub(V1,V2,x , x0 , sel_vc ,sel_vd)+rz1 .* randn(mz1,1) ; % I n i t i a l	simulated	value	--- measurement equations: 2.20 and 2.21 in paper SCADA measurments
 
	%=====================
	% 2) PMU-2 Measurements(FINEST)
	%=====================
	%z2f=[ z2f_1(sel va : sel vb ,1) ; z2f_2(sel_va : sel_vb ,1) ] ;
	mz2 = 2*(abs( sel_vb - sel_va))+2;
	%	Af = dfdx*Ts + eye(length(x));%jacobian(gx , x) ;
	z2 = zeros(mz2, k_end); % Creation	of	the array	for storing	the simulated	values	of z2
	z2(: ,1)= compute_z2fsub(V1,V2,x , x0 , sel_va ,sel_vb)+rz2 .* randn(mz2,1) ; % I n i t i a l	simulated	value	of z2
	Q = rx .* eye(ns , ns) ; % Process	noise cov of	the dynamical system correspnding t to	finest	time scale
	Q1 =(rx*n) .* eye(ns , ns) ; % Process	noise cov of	the dynamical system correspnding to	coarsest	time scale
	R2 = rz2 .* eye(mz2,mz2) ; % Measurement noise cov 	corresponding to	finest PMU
	R1 = rz1 .* eye(mz1,mz1) ; % Measurement noise cov corresponding to	coarsest PMU

	P = zeros(ns , ns , k_end);
	P(: ,: ,1) = Q; % P corresponding to PMU group with	finest	time	scale
	P1 = zeros(ns , ns , round( k_end/n)+1) ;
	P1(: ,: ,1) = Q1; % P corresponding to PMU group with	coarsest	time	scale
	P_21(: ,: ,1) = P1(: ,: ,1) ; % P corresponding to	coarse	sensors with	finest	time	scale
	x_2=zeros(ns , k_end);
	x_2(: ,1) = x0 ;

	denom_a=(P(: ,: ,1) + P_21(: ,: ,1));
	%	a1=(P(: ,: ,1))*denom a ;
	%	a2=(P_21(: ,: ,1))*denom a ;
	x_2(: ,1)= P(: ,: ,1) *(denom_a\x_21(: ,1))+ P(: ,: ,21) *(denom_a\x_est(: ,1));

	sum_11=P_21(: ,: ,1) + P(: ,: ,1) ;
	L11 = chol(sum_11 , 'lower');
	IL11=inv(L11) ;
	P2(: ,: ,1)=(P(: ,: ,1))*IL11*(P_21(: ,: ,1));
	x_sc(: ,1)=x_est1(: ,1) ;
	Tr(1)=trace(P(: ,: ,1));
	Tr_21(1)=trace(P_21(: ,: ,1));
	Tr_2(1)=trace(P2(: ,: ,1));
	Ob=zeros(k_end ,1) ;
	Ob(1 ,1)=ns ;
	Ob1=zeros(k_end ,1) ;
	Ob1(1 ,1)=ns ;
	st='preF' ;
	for k = 2: k_end
	tic
	k
	if k==k_Fstart+1 && k <= k_Fend
	st='fault' ;
	eval(sprintf('[f,x] = dx_for_9bus_%d_fault_sym' , iF)) ;
	Yr = Yr0( iF).fault ; 
    Rv = Rv0( iF).fault ;
	Y_big = Y_big0(iF).fault ;

	dx = f ;
	gx = dx.*Ts + x ;
	gx1 = dx.*(n*Ts) + x ;
	dfdx = jacobian(dx , x) ;
%
%======================================================
	% Construction	of Measurement System
	%
%======================================================
	Yr1 = real(Yr) ; 
    Yr2 = imag(Yr) ;
	I1 =(Yr1*E1-Yr2*E2) ;
	I2 =(Yr1*E2+Yr2*E1) ;

	Rv1 = real(Rv) ; 
    Rv2 = imag(Rv) ;
	V1 =(Rv1*E1-Rv2*E2) ; 
    V2 =(Rv1*E2+Rv2*E1) ;


    elseif	k==k_Fend+1
	st='postF' ;
	eval(sprintf('[f,x] = dx_for_9bus_%d_postF_sym' , iF)) ;
	Yr = Yr0( iF).postF ; 
    Rv = Rv0( iF).postF ;
	Y_big = Y_big0(iF).postF ;

	dx = f ;
	gx = dx.*Ts + x ;
	gx1 = dx.*(n*Ts) + x ;
	%dfdx = jacobian(dx , x) ;

	%
%======================================================
	% Construction	of Measurement System
	%
%======================================================
	Yr1 = real(Yr) ; 
    Yr2 = imag(Yr) ;
	I1 =(Yr1*E1-Yr2*E2) ;
	I2 =(Yr1*E2+Yr2*E1) ;
	Rv1 = real(Rv) ; 
    Rv2 = imag(Rv) ;
	V1 =(Rv1*E1-Rv2*E2) ; 
    V2 =(Rv1*E2+Rv2*E1) ;
	end
	new_V =(k==2 || k==k_Fstart+1 || k==k_Fend+1);

	x_sim(: , k) =(x_simu(: , k))+ sqrt(rx).* randn(ns,1) ; % Simulation	of	noisy	states

	%==============================
	% Independent	finest PMU Estimation
	%==============================
	z2(: , k) = compute_z2fsub(V1,V2,x , x_sim(:,k) ,sel_va ,sel_vb)+ sqrt(rz2).* randn(mz2,1) ;

	% 1)	Prediction	Steps=========
	[ x_est(: , k) , P(: ,: , k) ] = UKF_pred1(gx , x ,	x_est(: , k-1) ,(P(: ,: , k-1)) ,Q, Ts ,st); % a-priori estimate

 [z2_UKF, Pxy, Py]=UKF_meas(z2(: , k) , x_est(: , k) , P(: ,: , k) ,x , R2, V1,V2, sel_va , sel_vb); % a- priori estimate --------------------------------- steps 3,4 of algorithm 3

	Lpy=chol(Py, 'lower');
	ILpy=inv(Lpy) ;
	K = Pxy*transpose(ILpy)*ILpy ; % Kalman Gaincalculation
	x_est(: , k) = x_est(: , k) + K*(z2(: , k) - z2_UKF) ; %State update based on measurement
	P(: ,: , k) = P(: ,: , k) - K*Py*transpose(K) ; %Covariance Update
	Tr(k) = trace(P(: ,: , k));	% Ob(k*Ts ,	1)=rank([H; H*A; H*(A^2) ;H*(A^3) ]) ;

	%==============================
	% Independent Coarse PMU Estimation
	%==============================
	if	(round((k-1)/n) -(k-1)/n)==0
	ks = ks+1;
	z1(: , ks)= compute_z2fsub(V1,V2,x , x_sim(: , k) ,	sel_vc ,sel_vd)+ sqrt(rz1).* randn(mz1,1) ;

	%Prediction

	% 1)	Prediction	Steps=========
	[ x_est1(: , ks), P1(: ,: , ks)] = UKF_pred1(gx1 ,x ,x_est1(: , ks-1) ,(P1(: ,: , ks-1)) ,Q1,n*Ts , st); % a-priori estimate


	% 1)	Filtering	Steps=========
	[z1_UKF, Pxy1 , Py1]=UKF_meas(z1(: , ks),x_est1(: , ks), P1(: ,: , ks), x , R1, V1,V2, sel_vc , sel_vd); % a-priori estimate
	Pxy1;
	Py1;
	Lpy1 = chol(Py1 , 'lower'); 
    ILpy1=inv(Lpy1) ;

	K1 = Pxy1*inv(Py1) ; % Kalman Gain calculation
	x_est1(: , ks)= x_est1(: , ks)+ K1*(z1(: , ks)-z1_UKF) ; % State update based on measurement
	P1(: ,: , ks)=	P1(: ,: , ks)- K1*Py1*transpose(K1) ; % Covariance Update
	% Ob1(k*Ts ,	1)=rank([H1; H1*A1; H1*(A1^2) ;H1*(A1^3) ]) ;

	%=====================

 x_21(: , k)=x_est1(: , ks); 
 P_21(: ,: , k)=P1(: ,: , ks);

	%
%==================================================
	%2b) Steps to	calculate	x_21 and P_21 when k=nl + p
	%
%=================================================
	else
	p = k-n*(ks-1)-1;
	gx2 = dx.*(p*Ts)+x ;
	[x_21(: , k),P_21(: ,: , k)]=UKF_pred1(dx ,x ,x_21(: , k-1) ,(x_21(: ,:,k-1)) ,Q,Ts,st);
	end
	Tr_21(k)=trace(P_21(: ,: , k));
	x_sc(: , k)=x_est1(: , ks);

	%==============================
	% Fusion of Fine and Coarse PMUs based Estimates
	%==============================
	denom_a=(P(: ,: , k) + P_21(: ,: , k)); 	%a1=(P(: ,: , k))*denom a ;
	%a2=(P_21(: ,: , k))*denom a ;
	%a2=eye(ns , ns)-a1 ;
	x_2(: , k)=(P(: ,: , k))*(denom_a\x_21(: , k))+(P_21(: ,: , k))*(denom_a\x_est(: , k));

	%add 2=P 21(: ,: , k) + P(: ,: , k) ;
	Lp2 = chol(denom_a , 'lower');
	ILp2=inv(Lp2) ;
	P2(: ,: , k)=(P(: ,: , k))*ILp2 *(( P_21(: ,: , k))) ;
	Tr_2(k)=trace(P2(: ,: , k));
	toc
	end
	end
	end

	for	ii =1:ns
	figure(ii)
	k=[1: k_end ] ;
	plot(k ,x_simu(ii , k) , '-')
	hold on
	plot(k ,x_est(ii , k) , '-.m')
   	hold on
	plot(k ,x_21(ii , k) , '-r ')
    hold on
	%plot(k ,x_2(ii , k) , '-.b ')
	hold on
	title('Comparison')
	legend('True State ' , 'Only fast PMU estimation' , ' Only coarse PMU Estimate' , 'Fused Estimate');
	hold on

	end



	xTr= x_sim - x_est ;
	xTr_21= x_sim - x_21 ;
	xTr_2= x_sim - x_2 ;

	figure(ns+1)
	k=[1: k_end ];
	plot(k*Ts ,	Tr(k) ,	' r ')
	hold on	
	plot(k*Ts ,	Tr_21(k) ,	'k ')
	hold on	
	plot(k*Ts ,	Tr_2(k) ,	'm')
	hold on
	%t i t l e(' Trace ')
	xlabel(' time(s)')
	ylabel('Trace of estimation	error covariances ');
	legend('PMU trace','SCADA trace' ,'Fusion');

	toc
