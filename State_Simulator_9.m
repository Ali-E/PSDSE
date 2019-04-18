	function	[ x_sim ,	x_sim_2]= State_Simulator_9(varargin)
	if ~isempty(varargin)
		fault = varargin {1};
	else
		fault = [];
	end
 [Y, Yr, Rv, E0 , delta0 , Pm, RefNode , x0]= Calc_Y_bus_9(fault) ;
[ Sl_star , Vl , genbus_id , loadbus_id , xd , V_0 , V_0_A , Sg_star ,NOB, D, H, M, Gn, ns , T_end , Ts]= indi_conds ;
 load latest_n ;
 n=latest_n ;
 
 	if	isempty( fault)

		ode_fname = write_dx(E0, Pm, M, Yr, D, RefNode) ; % ToDo 
		sym_fname = write_dx_sym(E0, Pm, M, Yr, D, RefNode ,ode_fname);

	else

	for	iF = 1: length(fault)
    ode_fname(iF).preF= sprintf ( 'dx_for_9bus_%d_preF' , iF) ;
write_dx(E0, Pm, M, Yr(iF).preF , D, RefNode, de_fname(iF).preF);
	sym_fname(iF).preF = write_dx_sym(E0, Pm, M, Yr(iF).preF , D, RefNode ,ode_fname(iF).preF);
    ode_fname(iF).fault= sprintf ( 'dx_for_9bus_%d_fault' , iF) ;
write_dx(E0, Pm, M, Yr(iF).fault , D, RefNode, de_fname(iF).fault);
	sym_fname(iF).fault = write_dx_sym(E0, Pm, M, Yr(iF).fault , D, RefNode ,ode_fname(iF).fault);
    ode_fname(iF).postF= sprintf ( 'dx_for_9bus_%d_postF' , iF) ;
write_dx(E0, Pm, M, Yr(iF).postF , D, RefNode, de_fname(iF).postF);
	sym_fname(iF).postF = write_dx_sym(E0, Pm, M, Yr(iF).postF , D, RefNode ,ode_fname(iF).postF);

	end

 end


	% STEP 3: Use ODE45 to	solve	for x the dx :
	%

	if	isempty(fault)
	[T, x]= ode45(ode_fname ,	[0 T_end],x0) ;
	else

	%NOW THIS PART TACITLY ASSUMES THAT THERE IS ONLY ONE FAULT!
	for	iF = 1:1%length(fault)
	Fstart = fault(iF).time(1) ;
	Fend = fault(iF).time(2) ;
	[T1, x1]= ode45(ode_fname(iF).preF ,[0	Fstart], x0) ;
    [T2, x2]= ode45( ode_fname(iF).fault ,	[ Fstart Fend],x1(end, :)) ;
    [T3, x3]= ode45( ode_fname(iF).postF ,	[ Fend T_end], x2(end, :)) ;
	T = [T1;T2(2:end) ;T3(2:end)];
	x = [ x1 ; x2(2:end, :); x3(2:end, :)] ;
    end
 end

% STEP 4: ODE45 returns x at arbitrary interpolate those states	times T. Now
	% at	arbitrary T values	corresponding


	if	isempty(fault)
	t = 0:Ts: T_end ;
	y = interp1(T,x , t);
	else
	t1 = 0:Ts: Fstart ;
	t2 = Fstart+Ts:Ts: Fend ;
	t3 = Fend+Ts:Ts: T_end ;

	y1 = interp1(T1, x1 , t1);
	y2 = interp1(T2, x2 , t2);
	y3 = interp1(T3, x3 , t3);

	y = [ y1 ; y2 ; y3];
	TS=0:Ts*n: floor(T_end/(n*Ts))*n*Ts;

	yy1=interp1([ t1 t2 t3],y ,TS) ;
	end

	x_sim = transpose(y) ;
    x_sim_2 = transpose(yy1) ;

	k_end = floor(T_end/Ts) ;

