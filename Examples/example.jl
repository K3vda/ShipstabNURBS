#basis functions
begin
	P = [0 1 3 5 6 8 9;0 2 4 2 0 0 3]
	U = [0 0 0 0 2/5 3/5 3/5 1 1 1 1]
	n = size(P,2) - 1
	m = length(U) - 1
	p = m-n-1
	u = range(0, stop=1, length=100)
	funs = zeros(length(u),n+1)
	for i = 1:length(u)
	    funs[i,:] = BasisFuns(n,p,u[Int(i)],U)
	end

    
	plot(u,funs)
	plot!([2/5, 2/5],[0, 1], line=(:dash,2))
	plot!([3/5, 3/5],[0, 1], line=(:dash,2))
	png("Basisfuns")

    
	C = funs*P'
	plot(C[:,1],C[:,2], line=(:dot,2))
	plot!(P[1,:],P[2,:], line=(:dash,2))
	png("Example_curve")
end
#Bspline_example
begin
	P1 = [0 1 3 5 5 8 9;0 2 4 2 0 0 3]
	U1 = [0 0 0 0 2/5 3/5 3/5 1 1 1 1]
	n1 = size(P1,2) - 1
	m1 = length(U1) - 1
	p1 = m1-n1-1
	x=4
	u1 = range(0,1,length=1*10^x)
	C1 = zeros(size(P1,1),length(u1))
	for i = 1:length(u1)
    	C1[1,i] = bsplevalc(u1[i],U1,p1,P1[1,:])
    	C1[2,i] = bsplevalc(u1[i],U1,p1,P1[2,:])
	end

	plot(P1[1,:],P1[2,:],line=(:dot,2),show=false)
	plot!(C1[1,:],C1[2,:],line=(:dot,2),show=false)
	png("Bspline_example")
end
#Bspline_surface
begin
	Px2 = [0 3 6 9;0 3 6 9;0 3 6 9]'
	Py2 = [0 0 0 0;2 2 2 2;4 4 4 4]'
	Pz2 = [0 3 3 0;2 5 5 2;0 3 3 0]'
	U2 = [0 0 0 1/2 1 1 1]
	V2 = [0 0 0 1 1 1]
	n2 = size(Px2,2) - 1
	m2 = size(Px2,1) - 1
	h2 = length(U2) - 1
	k2 = length(V2) - 1
	
	q2 = k2-n2-1
	p2 = h2-m2-1
	
	x2 = 2
	u2 = range(0,1,length=30)
	v2 = range(0,1,length=30)
	
	lu2 = length(u2)
	lv2 = length(v2)
	Cx2 = zeros(lu2,lv2)
	Cy2 = zeros(lu2,lv2)
	Cz2 = zeros(lu2,lv2)
	
	for i = 1:length(u2)
	    for j = 1:length(v2)
	        Cx2[i,j] = bsplevals(u2[i],U2,v2[j],V2,p2,q2,Px2)
	        Cy2[i,j] = bsplevals(u2[i],U2,v2[j],V2,p2,q2,Py2)
	        Cz2[i,j] = bsplevals(u2[i],U2,v2[j],V2,p2,q2,Pz2)
	    end
	end
	plot(Cx2[1,:],Cy2[1,:],Cz2[1,:],show=false)
	for i = 1:length(u2)
	    plot!(Cx2[i,:],Cy2[i,:],Cz2[i,:])
	    plot!(Cx2[:,i],Cy2[:,i],Cz2[:,i])
	end
	for i = 1:size(Px2,1)
	   	plot!(Px2[i,:],Py2[i,:],Pz2[i,:])
	end
	for i = 1:size(Px2,2)
	    plot!(Px2[:,i],Py2[:,i],Pz2[:,i])
	end
	png("Bspline_surface")
end

#NURBS_spiral
begin
	P3 = [0            0
      0.55221      0.77547
     -0.62274       1.7993
      -2.7403      0.80463
      -2.9933      -2.3539
      0.22649      -4.7546
       4.8052      -3.0881
       6.1866       2.4768
       1.7955       7.4013
      -5.6108       6.4752
      -9.4769     -0.90493
       -5.236       -9.069
       4.7457      -10.392
       12.152      -2.3422
       9.6459       9.1974
      -2.0323       14.135
      -13.539       6.9797
      -14.385      -7.4159
      -2.4387      -16.962
       13.091      -12.482
       ]'
	U3 = [0 0 0 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 17 17 17]./17
	W3 = ones(1,size(P3,2))

	npts3 = 150
	nrbplot(P3,W3,U3,npts3)

	u3=0.3 # % for derivative
	n3 = size(P3,2) - 1
	m3 = length(U3) - 1
	p3 = m3-n3-1
	T3 = zeros(length(u3),2)
	C3 = T3
	for i = 1:length(u3)
		T3[i,:] = nrbder(P3,U3,W3,u3[i])
		C3[i,:],~ = nrbevalc(u3[i],U3,p3,P3,W3)
	end
	plot!(C3[:,1],C3[:,2],show=false)
	png("Image_spiral_nrbs")
end

#NURBSsurface
begin
	Px4 = [0 3 6 9;	0 3 6 9;0 3 6 9]'
	Py4 = [0 0 0 0;	2 2 2 2;4 4 4 4]'
	Pz4 = [0 3 3 0;2 5 5 2;0 3 3 0]'
	W4 = ones(4,3)
	U4 = [0 0 0 1/2 1 1 1]
	V4 = [0 0 0 1 1 1]
	
	nptsu4 = 20
	nptsv4 = 20
	nrbplot3(Px4,Py4,Pz4,W4,U4,V4,nptsu4,nptsv4,1)
	png("NURBSsurface")
end
