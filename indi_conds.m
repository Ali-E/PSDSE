	function [ Sl_star , Vl , genbus_id ,loadbus_id , xd , V_0 , V_0_A , Sg_star , NOB, D, H, M, Gn, ns , T_end , Ts]=indi_conds
	NOB=9;
	T_end=10;
	Ts=0.03;
	genbus_id =[1;2;3];
	Gn=length ( genbus_id ) ;
	ns=2*Gn-2;
	%nlbus id =[1;	2;	3;	4;	] ;
	%loadbus id=setdiff ([1:NOB] , nlbus_id ) ;
	loadbus_id =[5;	6;	8];
	V_0=[1.04;	1.025;	1.025];

 	V_0_A=pi /180*[0.00;	9.3;	4.7];

	i=sqrt (-1) ;

	Sl_star =[1.25-0.5*1 i ;	0.9-0.3*1 i ; 1-0.35*1 i ] ;


 	Vl = [0.996;	1.013;	1.016];

	Sg_star=[ 0.716 - 1 i .*0.27;	1.63 - 1 i *0.067;	0.85 + 1 i.*0.109];

	xd=[0.0608;	0.1198;	0.1813];
	H=[23.64;	6.4;	3.01];

	M = H./( pi *60) ;
	D = 5*M; %	LET US ASSUME UNIFORM DAMPING FOR MODEL REDUCTION
    	end
