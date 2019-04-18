
	function	[ sigma_pts , n]=UKF_calc_sigma( x_est , P )
	n=length (P) ;
	sigma_pts=zeros (n,2*n) ;
	A=chol (n*P) ;
	x_tilda= [	transpose(A) , -transpose(A) ] ;
	sigma_pts=repmat( x_est ,1 , 2*n) + x_tilda ;
	end
