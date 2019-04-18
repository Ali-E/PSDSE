	function	[Y, Yr, Rv, E0,	delta0 ,Pm,	RefNode ,x0] = Calc_Y_bus_9(varargin)
	NOB=9;
	T_end=10;
	Ts=0.03;
	genbus_id =[1;2;3];
	Gn=length (genbus_id) ;
	ns=2*Gn-2;
	%nlbus id =[1;	2;	3;	4;	] ;
	%loadbus id=setdiff ([1:NOB] , nlbus_id) ;
	loadbus_id =[5;	6;	8];
	V_0=[1.04;	1.025;	1.025];

 V_0_A=pi /180*[0.00;	9.3;	4.7];

	i=sqrt (-1) ;

	Sl_star =[1.25-0.5*1 i ;	0.9-0.3*1 i ; 1-0.35*1 i] ;


 Vl =[0.996;	1.013;	1.016];

	Sg_star=[0.716 - 1 i .*0.27;	1.63 - 1 i *0.067;	0.85 + 1 i.*0.109];

	xd=[0.0608;	0.1198;	0.1813];
	H=[23.64;	6.4;	3.01];

	M = H./(pi *60) ;
	D = 5*M;

	GG=Gn + NOB;
	if	~isempty(varargin)
	fault = varargin {1};
	else
	fault =[] ;
	end

	% Line data
	Line(1) = struct('b', 0.0 ,'bus1', 1 ,'bus2', 4 ,'r',0.0 ,'x',0.0576 ,'a',1) ;
	Line(2) = struct('b',	0.176 ,	'bus1', 4 ,	'bus2', 5 ,	'r',	0.010 ,	'x',	0.085 ,'a',1) ;
	Line(3) = struct('b',	0.158 ,	'bus1'	, 4 ,	'bus2'	, 6 ,	'r'	,	0.017 ,	'x'	,	0.092 ,'a',1) ;
	Line(4) = struct('b',	0.306 ,	'bus1'	, 5 ,	'bus2'	, 7 ,	'r'	,	0.032 ,	'x'	,	0.161 ,'a',1) ;
	Line(5) = struct('b',	0.358 ,	'bus1'	, 6 ,	'bus2'	, 9 ,	'r'	,	0.039 ,	'x'	,	0.170 ,'a',1) ;
	Line(6) = struct('b',	0.0 ,	'bus1'	, 7 ,	'bus2'	, 2 ,	'r'	,	0.0 ,	'x'	,	0.0625 ,'a',1) ;
	Line(7) = struct('b', 0.0 ,'bus1', 9 ,'bus2', 3 ,'r', 0.0 ,'x', 0.0586 ,'a',1) ;
	Line(8) = struct('b',	0.1490 ,	'bus1'	, 7 ,	'bus2'	, 8 ,	'r'	,	0.0085 ,	'x'	,	0.072 ,'a',1) ;
	Line(9) = struct('b',	0.2090 ,	'bus1'	, 9 ,	'bus2'	, 8 ,	'r'	,	0.0119 ,	'x'	,	0.1008 ,'a',1) ;

	%
%=====================================================================
	% This Y matrix is the network Y matrix for 9 Bus system from the Pete
	% Sauer and M. A.	Pai Book.
	%
%=====================================================================

	if	~isempty(fault)
	Y0 = zeros(NOB) ;
	for p = 1: length(Line)
	bus1 = Line(p).bus1 ;
	bus2 = Line(p).bus2 ;
	a=Line(p).a ;
	b = Line(p).b;
	z = Line(p).r + 1i*Line(p).x ;
	Y0(bus1 , bus2) = -1/(a*z) ;
	Y0(bus2 , bus1) = -1/(a*z) ;
	Y0(bus1 , bus1) = Y0(bus1 , bus1) + 1/((a^2)*z) + 1i*b;
	Y0(bus2 , bus2) = Y0(bus2 , bus2) + 1/((a^2)*z) + 1i*b;
	end
	for	iF = 1: length(fault)	% If	multiple	faults ,	then	this'iF'	indicate a	specific	fault
	% Find the	line	index 53	
    for	ln=1:length(Line)
	if	isempty(setdiff([Line(ln).bus1 Line(ln).bus2] ,fault(iF).nodes))
	break
	end
	end	% Hence ”ln” has the	faulted	line	index
	Y_fault = Y0;
    Y_fault(fault(iF).nodes(1) , fault(iF).nodes(1)) =Y_fault(fault(iF).nodes(1) , fault(iF).nodes(1)) -...
1/((Line(ln).a)^2*(Line(ln).r+1i *Line(ln).x)) +(1/ fault(iF).location) /((Line(ln).a)^2*(Line(ln).r+1i *Line(ln).x)) ;
Y_fault(fault(iF).nodes(2) , fault(iF).nodes(2)) =Y_fault(fault(iF).nodes(2) , fault(iF).nodes(2)) -...
1/((Line(ln).a)^2*(Line(ln).r+1i *Line(ln).x)) +(1/(1- fault(iF).location)) /((Line(ln).a)^2*(Line(ln).r+1i *Line(ln).x)) ;
	Y_fault(fault(iF).nodes(1) , fault(iF).nodes(2)) = 0;
    Y_fault(fault(iF).nodes(2) , fault(iF).nodes(1)) = 0;
	Y_postF = Y0;
	Y_postF(fault(iF).nodes(1) , fault(iF).nodes(1)) = Y_postF(fault(iF).nodes(1) , fault(iF).nodes(1)) - 1i*Line(ln).b - 1/((Line(ln).a)^2*(Line(ln).r+1i *Line(ln).x)) ;
	Y_postF(fault(iF).nodes(2) , fault(iF).nodes(2)) = Y_postF(fault(iF).nodes(2) , fault(iF).nodes(2)) - 1i*Line(ln).b - 1/((Line(ln).a)^2*(Line(ln).r+1i * Line(ln).x)) ;
	Y_postF(fault(iF).nodes(1) , fault(iF).nodes(2)) = 0; 
    Y_postF(fault(iF).nodes(2) , fault(iF).nodes(1)) = 0;

	Y(iF).preF = Y0;
	Y(iF).fault = Y_fault ; 
    Y(iF).postF = Y_postF ;
	end
    if	isempty(setdiff([Line(ln).bus1 Line(ln).bus2] ,fault(iF).nodes))
 	Y = zeros(NOB) ;
	for p = 1: length(Line)
	a=Line(p).a ;
	bus1 = Line(p).bus1 ;
	bus2 = Line(p).bus2 ;
	b = Line(p).b;
	z = Line(p).r + 1i*Line(p).x ;
	Y(bus1 , bus2) = -1/(a*z) ;
	Y(bus2 , bus1) = -1/(a*z) ;
	Y(bus1 , bus1) = Y(bus1 , bus1) + 1/(a^2*z) + 1i*b;
	Y(bus2 , bus2) = Y(bus2 , bus2) + 1/(a^2*z) + 1i*b;
	end

 end

	%
%=====================================================================
	% Calculate	equivalent impedance for	loads .
	%
%=====================================================================

 %Sl_star =[1.25 -0.5*1i ; 0.9-0.3*1i ; 1-0.35*1i];% CONJUGATE OF COMPLEX 99 % LOAD POWER.
 %Vl =[0.996;	1.013;	1.016];
	Yload= Sl_star./Vl.^2;
	%Yload=Yloas1 ;
	%loadbus_id =[5 6	8]';%[5;7;9];
	%
%====================================================================
	% Calculate	i n i t i a l	generator	internal	voltage Eg0 and i n i t i a l	rotor
	% angle	delta0 .	This	will	be used in	future .
	%
%=====================================================================

	%xd =[0.0608;	0.1198;	0.1813];
	%V_0 =[1.04;	1.025;	1.025]; % Generator bus voltage magnitude
	%V_0_A = pi /180*[0.00; 9.3; 4. 7]; % Generator bus volate angle(*IN RADIAN*)
	%Sg_star =[0.716 - 1i .*0.27; 1.63 - 1i *0.067; 0.85 + 1i .*0.109]; % The
	% generator i n i t i a l condition COMPLEX POWER CONJUGATE formed from power
	% flow	solution from the Pete Sauer/M.A.	Pai book

 Eg0 = V_0 +(1i *xd) .*(Sg_star ./V_0) ; % Generator internal voltage vector %NEED TO VERIFY
 % at	i n i t i a l	condition
				
	E0 = abs(Eg0) ; % The magnitude remains more or
	% less same	of	internal	voltage .	This
	derl0 = angle(Eg0) + V_0_A ;

	RefNode = 1;
	delta0 = derl0(setdiff(1:Gn, RefNode)) - derl0(RefNode) ;%[derl0(1)-derl0(3) ; derl0(2)-derl0(3)] ;
				
 % To find out the delta0 ,	i.e ,	the	i n i t i a l			
	% rotor	angle , we need to add the

	% NEEDS FURTHER INSPECTION.

 %	bus voltage	angle	to	it .
%====================================================================
	%Now we will	modify the Y_bus matrix ;
%1) by adding the	equivalent	admittance of	load value to the	load bus diagonal	elements
	%
%=====================================================================

	if	~isempty(fault)

	for	iF = 1: length(fault)
	for p = 1: length(loadbus_id)
	Y(iF).preF(loadbus_id(p),loadbus_id(p))=Y(iF).preF(loadbus_id(p),loadbus_id(p))+ Yload(p);
	Y(iF).fault(loadbus_id(p) , loadbus_id(p)) = Y(iF).fault(loadbus_id(p) , loadbus_id(p)) + Yload(p) ;
	Y(iF).postF(loadbus_id(p) , loadbus_id(p)) = Y(iF).postF(loadbus_id(p) , loadbus_id(p)) + Yload(p) ;
	end
	end

	else

	for p = 1: length(loadbus_id)
	Y(loadbus_id(p) , loadbus_id(p)) = Y(loadbus_id(p) ,loadbusid(p)) + Yload(p) ;
	end

 end

	%
%====================================================================
	%	Now Extended Y-bus to	include	the	internal Gen nodes
	%
%=====================================================================
	%genbus_id =[1 ; 2 ; 3] ;


	internal_id = NOB*(ones(length(genbus_id) ,1)) +[1: length(genbus_id)]';
	%[10:12]; % 10-->1, 11--> 2 , 12-->3

	if	~isempty(fault)

	for	iF = 1: length(fault)



	Yext(iF).preF =[Y(iF).preF zeros(NOB,Gn) ;	zeros(Gn,GG)] ;
	Yext(iF).fault =[Y(iF).fault	zeros(NOB,Gn) ;	zeros(Gn,GG)] ;
	Yext(iF).postF =[Y(iF).postF zeros(NOB,Gn) ;	zeros(Gn,GG)] ;

	for p=1:Gn
	Yext(iF).preF(genbus_id(p) , internal_id(p)) =	-1/(1i *xd(p,1)) ; % off-diagonal	entries
	Yext(iF).preF(internal_id(p) , genbus_id(p)) =...
-1/(1i *xd(p,1)) ;
	Yext(iF).preF(genbus_id(p) , genbus_id(p)) = Yext(iF).preF(genbus_id(p) , genbus_id(p)) + 1/(1i *xd(p,1)) ; % diagonal	entries
	Yext(iF).preF(internal_id(p) , internal_id(p)) =1/(1i *xd(p,1)) ;

	Yext(iF).fault(genbus_id(p) , internal_id(p)) =	-1/(1i *xd(p,1)) ; % off-diagonal	entries
	Yext(iF).fault(internal_id(p) , genbus_id(p)) =-1/(1i *xd(p,1)) ;
	Yext(iF).fault(genbus_id(p) , genbus_id(p)) = Yext(iF).fault(genbus_id(p) , genbus_id(p)) + 1/(1i * xd(p,1)) ; % diagonal entries
	Yext(iF).fault(internal_id(p) , internal_id(p)) =1/(1i *xd(p,1)) ;

	Yext(iF).postF(genbus_id(p) , internal_id(p)) =	-1/(1i *xd(p,1)) ; % off-diagonal	entries
	Yext(iF).postF(internal_id(p) , genbus_id(p)) =...
-1/(1i *xd(p,1)) ;
	Yext(iF).postF(genbus_id(p) , genbus_id(p)) = Yext(iF).postF(genbus_id(p) , genbus_id(p)) + 1/(1i * xd(p,1)) ; % diagonal entries
Yext(iF).postF(internal_id(p) , internal_id(p)) =1/(1i *xd(p,1)) ;
	end
	end

	else

	Yext =[Y zeros(NOB,Gn) ;	zeros(Gn,GG)] ;
	for p=1:Gn
	Yext1(genbus_id(p) , internal_id(p)) = -1/(1i *xd(p,1)) ;
% off-diagonal	entries
	Yext1(internal_id(p) , genbus_id(p)) = -1/(1i *xd(p,1)) ;
	Yext1(genbus_id(p) , genbus_id(p)) = Yext1(genbus_id(p), genbus_id(p)) + 1/(1i *xd(p,1)) ; % diagonal
	Yext1(internal_id(p) , internal_id(p)) = 1/(1i *xd(p,1));
	end

 end

	%
%=====================================================================
	% APPLY KRON REDUCTION: Let's split Yext into Yext =[Yll Ylg ; Ygl Ygg] ;
	%
%=====================================================================

	if	~isempty(fault)

	for	iF = 1: length(fault)

Yll = Yext(iF).preF(1:NOB,1:NOB) ;
	Ylg = Yext(iF).preF(1:NOB,(NOB+1) :(NOB + Gn)) ;
	Ygl = Yext(iF).preF((NOB+1) :(NOB + Gn) ,1:NOB) ;
	Ygg = Yext(iF).preF((NOB+1) :(NOB + Gn) ,(NOB+1) :(NOB +Gn)) ;

	Yr(iF).preF = Ygg - Ygl*inv(Yll)*Ylg
	Rv(iF).preF = - inv(Yll)*Ylg

	Yll = Yext(iF).fault(1:NOB,1:NOB) ;
	Ylg = Yext(iF).fault(1:NOB,(NOB+1) :(NOB + Gn)) ;
	Ygl = Yext(iF).fault((NOB+1) :(NOB + Gn) ,1:NOB) ;
	Ygg = Yext(iF).fault((NOB+1) :(NOB + Gn) ,(NOB+1) :(NOB + Gn)) ;

	Yr(iF).fault = Ygg - Ygl*inv(Yll)*Ylg
	Rv(iF).fault = - inv(Yll)*Ylg

	Yll = Yext(iF).postF(1:NOB,1:NOB) ;
	Ylg = Yext(iF).postF(1:NOB,(NOB+1) :(NOB + Gn)) ;
	Ygl = Yext(iF).postF((NOB+1) :(NOB + Gn) ,1:NOB) ;
	Ygg = Yext(iF).postF((NOB+1) :(NOB + Gn) ,(NOB+1) :(NOB + Gn)) ;

	Yr(iF).postF = Ygg - Ygl*inv(Yll)*Ylg
	Rv(iF).postF = - inv(Yll)*Ylg

	end


	else
	Yll = Yext1(1:NOB,1:NOB) ;
	Ylg = Yext1(1:NOB,(NOB+1) :GG) ;
	Ygl = Yext1((NOB+1):GG,1:NOB) ;
	Ygg = Yext1((NOB+1):GG,(NOB+1):GG) ;

	Yr = Ygg - Ygl*inv(Yll)*Ylg ;
	Rv = - inv(Yll)*Ylg ;
	end

	%
%====================================================================
	% ADDITONALY: We store	various model parameters and	i n i t i a l

	% to pass such as Pm, M, H for	future	use .
	%
%====================================================================
	Pe=zeros(Gn) ;
	if	~isempty(fault)
	% NOTE: WE ASSUME THERE IS ONLY ONE FAULT
	for	iF = 1:1%length(fault)
	for p = 1: length(genbus_id)
	Pe(p) = 0;
	for q = 1: length(genbus_id)
	Pe(p) = Pe(p) + E0(p)*E0(q) *(real(Yr(iF).preF(p, q))*cos(derl0(p)-derl0(q)) + imag(Yr(iF).preF...
(p, q))*sin(derl0(p)-derl0(q))) ;
	end
	end
	end
	else
	for p = 1: length(genbus_id)
	Pe(p) = 0;
	for q = 1: length(genbus_id)
	Pe(p) = Pe(p) + E0(p)*E0(q) *(real(Yr(p, q))*cos(derl0(p)-derl0(q)) + imag(Yr(p, q))*sin(derl0(p)-derl0(q))) ;
	end
	end
	end

 Pm = Pe; % Because , Pm = Pe in the	equilibrium/operating

H =[23.64;	6.4;	3.01];
M = H./(pi *60) ;
	D = 2*M; %	LET US ASSUME UNIFORM DAMPING FOR MODEL REDUCTION
	Ts = 0.01; % Chosen sampling time
	T_end = 10;
	x0 = zeros((Gn-1)*2 ,1) ;
% 	x0=[delta0(1:Gn-1,1) ;	zeros(Gn-1,1)] ;
	end
