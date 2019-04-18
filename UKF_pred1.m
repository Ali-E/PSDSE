function[ x_est_new , P_new ] = UKF_pred1(gx , x ,	x_est , P,Q, Ts	,st	) % a-priori	estimate
global A	
[ sigma_pts , n]=UKF_calc_sigma( x_est , P ) ;
	g_sigma=zeros (n,2*n) ;
	ffnameis=(sprintf( ' dx_for_9bus_1_%s ' , st ) ) ;
	for	ii =1:2*n
	[ t , xx]=ode45(@(t,x) write_dx(t,x,A),	[0 , Ts] ,	sigma_pts (: , ii ) ) ;
	g_sigma (: , ii )=transpose (xx(end , : ) ) ;
	end
	x_est_new= (1/(2*n) )*sum( g_sigma ,2) ;

	sums=(g_sigma(: ,1)-x_est_new )*transpose( g_sigma(: ,1)- x_est_new ) ;
	for	i =2:2*n
	sums=sums + ( g_sigma(: , i )-x_est_new )*transpose( g_sigma (: , i )- x_est_new );
	end
	P_new=(1/(2*n) )*sums + Q;
