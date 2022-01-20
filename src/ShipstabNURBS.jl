module ShipstabNURBS

using Plots

export BasisFuns,Findspan,
nrbplot3,nrbplot,
bsplevalc,bsplevals,bsplder,
nrbevalc,nrbevals,nrbder,
nrbcircle,nrbcylinder,
revolX,revolX_half


function FindSpan(u,U)
	#
	# Short description:
	# 	computes the span index and the multiplicity of a a value u using a
	# 	binary search algorithm. Only for clamped b-splines and NURBS. 
	#
	# Input parameters: 
	# 	u: value between 0 and 1 for which we want to find the span index and
	#   multiplicity
	#
	#   U: the knot vector
	#
	# Output paramters:
	#   SpanIndex: the span index of u. if u belongs to the open interval
	#   [Uk,Uk+1[, we say that it belongs to the k-th knot index. In MatLab,
	#   indexes start from 1 and not 0. So to be mathematically correct,
	#   FindSpan computs the span index + 1. 
	#
	#   s: the multiplicity of u. If u is no knot, it has a multiplicity of 0
	#
	# ------------------------------------------------------------------------

	# % first, we delete the first and last elements of the knot vector to avoid
	# % problems while performing de Boor's algorithm. (later, we add 1 to the
	# % found span index because we deleted the first element)
	U = U[2:end-1]
	s=0
	m = length(U) - 1
	low = 1
	high = m+1
	mid = round(Int,(low+high)/2)
	if isempty(findall(isequal(u),U)) == 0 && u!=U[end]
    	s = length(findall(isequal(u),U))
    	SpanIndex = findlast(isequal(u),U) + 1 # (+1, because we set U(1) = [])
    	return Int(SpanIndex),s
	end

	if U[end] == u
    	SpanIndex = findfirst(isequal(u),U)
    	return Int(SpanIndex),s
	end
	
	SpanIndex = mid
	while u < U[Int(mid)] || u >= U[Int(mid+1)]
    	if u < U[Int(mid)] 
        	high = mid
    	else
        	low = mid
		end
    	mid = round((high+low)/2)
    	if mid == SpanIndex
        	mid = mid - 1
    	end
   		SpanIndex = mid
	end
	SpanIndex = SpanIndex + 1
	return Int(SpanIndex),s
end

function BasisFuns(n,p,u,U)
	#
	# Short description: 
	#
	# compute the coefficients N0,p(u), N1,p(u), ..., Nn,p(u) for u in [0,1]
	#
	# Input parameters:
	#
	# n: there are n+1 control points
	# p: degree of the curve
	# u: value of the parameter u
	# U: knot vector
	#
	# Output parameters:
	#
	# Array containing the n+1 values of basis functions for a given u
	#
	# (algorithm explained in notes)
	# ------------------------------------------------------------------------
	k = FindSpan(u,U)[1]
	N = zeros(1,n+1)
	N[k] = 1
	for d = 1:p
	    if k-d > 0
    	    N[k-d] = ((U[k+1] - u)/(U[k+1] - U[k-d+1]))*N[k-d+1]
	        if k-d == n+1
	            return N
	        end
	        for i = k-d+1:k-1
	            N[i] = ((u - U[i])/(U[i+d]-U[i]))*N[i] + ((U[i+d+1] - u)/(U[i+d+1] - U[i+1]))*N[i+1]
	        end
	        N[k] = ((u - U[k])/(U[k+d] - U[k]))*N[k]
	    else
	        for i = 1:k-1
	            N[i] = ((u - U[i])/(U[i+d]-U[i]))*N[i] + ((U[i+d+1] - u)/(U[i+d+1] - U[i+1]))*N[i+1]
	        end
	        N[k] = ((u - U[k])/(U[k+d] - U[k]))*N[k]
	    end
	end
	return N
end

function bsplevalc(u,U,p,P)
	# U [1 x m+1]
	#
	# Short description:
	#   evaluates a (clamped) B spline curve for a given value for u using De
	#   Boors's algorithm
	#
	# Input parameters:
	#   u: a value in the knot vector U between 0 and 1
	#
	#   p: the degree of the NURBS curve
	#
	#   P: a [1 X n+1] vector with the cartesian x OR y coördinates
	#   of the n+1 control points
	#
	#   W: matrix containing the weights corresponding to the control points
	#   (see function nrbplot)
	#
	#   U: Knot vector [1 X m+1]
	#
	# Output parameters: 
	#   C: the cartesian coördinates of the point on the B spline curve for a
	#   given value for u
	# ------------------------------------------------------------------------

	k,s = FindSpan(u,U)
	h = p-s
	if h <= 0
    	C = P[k-s]
    	return C
	end

	B = zeros(p-s+1,h+1)
	B[:,1] = P[k-p:k-s]

	for r = 1:h
    	for i = r+1:h+1
      		a = (u - U[i+k-p-1])/(U[i+k-r] - U[i+k-p-1])
      		B[i,r+1] = (1-a)*B[i-1,r]+a*B[i,r]
    	end
	end
	C = B[h+1,h+1]
	return C
end

function bsplevals(u,U,v,V,p,q,P)
	#
	# Short description:
	#   evaluates a (clamped) B-spline surface for given values for u and v using 
	#   De Boors's algorithm
	#
	# Input parameters:
	#   u,v: a value in the knot vectors U and V between 0 and 1
	#
	#   p,q: the degree in respectively the u and v direction
	#
	#   P: cartesian coordinates of the control points (x,y or z direction, see
	#   input parameters of the function nrbplot3)
	#
	#   W: matrix containing the weights corresponding to the control points
	#   (see function nrbplot3)
	#
	# Output parameters: 
	#   C: the cartesian coördinates of the point on the B spline surface for the
	#   given values for u and v
	#
	# ------------------------------------------------------------------------

	c,s = FindSpan(u,U)
	Q = zeros(1,p-s+length(P))
	for i = c-p:c-s
    	Q[i] = bsplevalc(v,V,q,P[i,:])
	end
	C = bsplevalc(u,U,p,Q)
	return C
end

function bsplder(P,U,u)
	Uprime = U[2:end-1]
	Q = zeros(size(P,1),size(P,2)-1)
	n = size(P,2) - 1
	m = length(U) - 1
	p = m-n-1

	for i = 1:size(Q,1)
    	for j = 1:size(Q,2)
	        Q[i,j] = p *(P[i,j+1]-P[i,j])/(U[j+p+1]-U[j+1])
	    end
	end

	T = zeros(1,size(Q,1))
	for i = 1:size(Q,1)
	    T[i] = bsplevalc(u,Uprime,p-1,Q[i,:])
	end
	return T
end

function nrbevalc(u,U,p,P,W)
	#
	# Short description:
	#   evaluates a (clamped) NURBS curve for a given value for u using De
	#   Boors's algorithm
	#
	# Input parameters:
	#   u: a value in the knot vector U between 0 and 1
	#
	#   p: the degree of the NURBS curve
	#
	#   P: a [2 X n+1] matrix with on the first row the cartesian x-coördinates
	#   of the n+1 control points and on the second row the associated
	#   y-coördinates.
	#
	#   W: matrix containing the weights corresponding to the control points
	#   (see function nrbplot)
	#
	#   U: Knot vector [1 X m+1]
	#
	# Output parameters:
	#   C: the cartesian coördinates of the point on the NURBS curve for a
	#   given value for u
	#
	#   Wn: the weight of the calculated point on the NURBS curve (needed for
	#   the calculation of NURBS surfaces)
	# ------------------------------------------------------------------------

	Pw = P.*W
	Pw = vcat(Pw, W)
	Qw = zeros(size(Pw,1),1)
	for i = 1:size(Pw,1)
    	Qw[i] = bsplevalc(u,U,p,Pw[i,:])
	end
	C = Qw[1:end-1]./Qw[end]
	Wn = Qw[end]
	return C,Wn
end

function nrbder(P,U,W,u)
	Pw = [P.*W;W]
	Uprime = U[2:end-1]
	Q = zeros(size(Pw,1),size(Pw,2)-1)
	n = size(Pw,2) - 1
	m = length(U) - 1
	p = m-n-1

	for i = 1:size(Q, 1)
	    for j = 1:size(Q,2)
	        Q[i,j] = p *(Pw[i,j+1]-Pw[i,j])/(U[j+p+1]-U[j+1])
	    end
	end

	Ader = zeros(1,size(Q,1))
	for i = 1:length(Ader)
	    Ader[i] = bsplevalc(u,Uprime,p-1,Q[i,:])
	end

	C,Wn = nrbevalc(u,U,p,P,W)
	T = zeros(1,length(Ader)-1)
	for i = 1:length(Ader)-1
	    T[i] = (Ader[i]-Ader[end]*C[i])/Wn
	end
	return T
end

function nrbevals(u,U,v,V,p,q,P,W)
	#
	# Short description:
	#   evaluates a (clamped) NURBS surface for given values for u and v using 
	#   De Boors's algorithm
	#
	# Input parameters:
	#   u,v: a value in the knot vectors U and V between 0 and 1
	#
	#   p,q: the degree in respectively the u and v direction
	#
	#   P: cartesian coordinates of the control points (x,y or z direction, see
	#   input parameters of the function nrbplot3)
	#
	#   W: matrix containing the weights corresponding to the control points
	#   (see function nrbplot3)
	#
	# Output parameters: 
	#   C: the cartesian coördinates of the point on the NURBS surface for the
	#   given values for u and v
	#
	# ------------------------------------------------------------------------

	c,s = FindSpan(u,U)
	Q = zeros(2,p-s+length(P))
	for i = c-p:c-s
    	Q[1,i] = nrbevalc(v,V,q,P[i,:]',W[i,:]')[1][1]
		Q[2,i] = nrbevalc(v,V,q,P[i,:]',W[i,:]')[2]
	end
	C, ~ = nrbevalc(u,U,p,Q[1,:]',Q[2,:]')
	return C
end

function nrbplot(P,W,U,npts)
	#
	# Short description:
	#
	#   Constructs a (clamped) NURBS curve
	#
	# Input parameters:
	#
	#   P: a [2 X n+1] matrix with on the first row the cartesian x-coördinates
	#   of the n+1 control points and on the second row the associated
	#   y-coördinates. 
	#
	#   W: a [n+1] line vector containing the weights for each control point
	#
	#   U: knot vector
	#
	#   npts: the number of points for which the NURBS curve will be calculated
	#
	# Output parameters:
	#
	#   C: a matrix with on the first row the calculated cartesian
	#   x-coördinates and on the second row the associated y-coördinates
	# ----------------------------------------------------------------------

	n = size(P,2) - 1
	m = length(U) - 1
	p = m-n-1
	u = range(0,1,length=npts)
	C = zeros(2,length(u))

	for i = 1:length(u)
    	C[:,i], ~ = nrbevalc(u[i],U,p,P,W)
	end
	
	plot(P[1,:],P[2,:], line=(:dash,2))
	plot!(C[1,:],C[2,:], line=(:dot,2))
	return C
end

function nrbcircle(npts)
	#
	# Short description:
	#   construcs a NURBS circle evaluated at npts points
	#
	# Output parameters:
	#   P: control points
	#   W: weight vector
	#   U: knot vector
	#   C: carthesian coördinates of the calculated points
	#
	# (see input parameters of nrbplot for a more detailed description)
	#
	#   --------------------------------------------------------------------

	P = [1 1 0 -1 -1 -1 0 1 1;0 1 1 1 0 -1 -1 -1 0]
	W = [1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1]
	U = [0 0 0 1/4 1/4 1/2 1/2 3/4 3/4 1 1 1]

	C = nrbplot(P,W,U,npts)
	return P,W,U,C
end

function nrbplot3(Px,Py,Pz,W,U,V,nptsu,nptsv,flag_plot)
	# Short description:
	#
	#   Constructs a (clamped) NURBS surface
	#
	# Input parameters:
	#
	#   Px, Py, Pz: [m+1 X n+1] matrices with the cartesian coordinates of the
	#   control points for respectively the x, y and z-axis.
	#
	#          ------------------- > U direction
	#          |p1,1        p1,n+1 |
	#          | .     .           |
	#          | .      .          |
	#          | .       .         |
	#          |pm+1,1     pm+1,n+1|
	#          ----------------
	#          |
	#          v V direction
	#
	#   W: the weight assigned to each control point, [m+1 X n+1] matrix.
	#
	#   U,V: knot vectors in respectively the u- and v-direction
	#
	#   nptsu, nptsv: the number of points that subdivide the parametric space
	#   in respectively the u- and v-direction. It determines the fineness of
	#   the mesh in the physical space.
	#
	# Output parameters:
	#
	#   Cx,Cy,Cz: cartesian coördinates of the NURBS surface in respectively
	#   the x,y and z direction
	#   --------------------------------------------------------------------

	n = size(Px,2) - 1
	m = size(Px,1) - 1
	h = length(U) - 1
	k = length(V) - 1

	q = k-n-1
	p = h-m-1

	u = range(0,1,length=nptsu)
	v = range(0,1,length=nptsv)

	lu = length(u)
	lv = length(v)
	Cx = zeros(lu,lv)
	Cy = zeros(lu,lv)
	Cz = zeros(lu,lv)

	for i = 1:length(u)
    	for j = 1:length(v)
        	Cx[i,j] = nrbevals(u[i],U,v[j],V,p,q,Px,W)[1]
        	Cy[i,j] = nrbevals(u[i],U,v[j],V,p,q,Py,W)[1]
        	Cz[i,j] = nrbevals(u[i],U,v[j],V,p,q,Pz,W)[1]
    	end
	end
	if flag_plot == 1
		plot(Cx[1,:],Cy[1,:],Cz[1,:])
    	for i = 2:length(u)
        	plot!(Cx[i,:],Cy[i,:],Cz[i,:])
    	end
    	for i =1:length(v)
        	plot!(Cx[:,i],Cy[:,i],Cz[:,i])
    	end
     	# for i = 1:size(Px,1)
      #    	plot3(Px(i,:),Py(i,:),Pz(i,:),'ko-','Linewidth',2)
     	# end
     	# for i = 1:size(Px,2)
      #    	plot3(Px(:,i),Py(:,i),Pz(:,i),'ko-','Linewidth',2)
     # 	# end
    	# xaxis!("x-axis")
    	# yaxis!("y-axis")
    	# zaxis!("z-axis")
	end
	return Cx, Cy, Cz
end

function nrbcylinder(nptsu,nptsv)
	#
	# Short description:
	#   constructs a NURBS cylinder evaluated at nptsu x nptsv points
	#
	# Output parameters:
	#
	#   see input parameters of the function nrbplot3 for a detailed
	#   description
	#
	#   --------------------------------------------------------------------

	Pz = [1 1 0 -1 -1 -1 0 1 1;1 1 0 -1 -1 -1 0 1 1]'
	Py = [0 1 1 1 0 -1 -1 -1 0;0 1 1 1 0 -1 -1 -1 0]'
	Px = [-1 -1 -1 -1 -1 -1 -1 -1 -1; 1 1 1 1 1 1 1 1 1]'
	
	W = [1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1; 1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1]'
	U = [0 0 0 1/4 1/4 1/2 1/2 3/4 3/4 1 1 1]
	V = [0 0 1 1]

	nrbplot3(Px,Py,Pz,W,U,V,nptsu,nptsv,1)
	return Px,Py,Pz,W,U,V
end

function revolX(P,W,U,nptsu,nptsv,flag_plot)
	# Constructs a surface of revolution by means of a b-spline curve in 2D.

	Pj = zeros(3,size(P,2))
	Pj[1,:] = P[1,:]
	Pj[3,:] = P[2,:]
	Wj = W

	Pi = [0 0 0 0 0 0 0 0 0;1 1 0 -1 -1 -1 0 1 1;0 1 1 1 0 -1 -1 -1 0]
	Wi = [1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1]
	V = [0 0 0 1/4 1/4 1/2 1/2 3/4 3/4 1 1 1]

	nPj = size(Pj,2)
	nPi = size(Pi,2)

	xj = Pj[1,:]
	zj = Pj[3,:]

	Px = zeros(nPj,nPi)
	Py = Px
	Pz = Px

	for i = 1:length(xj)
    	Px[i,:] = xj[i].* ones(1,nPi)
    	Py[i,:] = zj[i].* Pi[2,:]
    	Pz[i,:] = zj[i].* Pi[3,:]
	end

	W = zeros(size(Px))
	for i = 1:nPj
    	W[i,:] = Wj[i].*Wi
	end
	nrbplot3(Px,Py,Pz,W,U,V,nptsu,nptsv,flag_plot)
	return Px,Py,Pz,U,V,W
end

function revolX_half(P,W,U,nptsu,nptsv,flag_plot)
	# Constructs a surface of revolution by means of a b-spline curve in 2D.

	Pj = zeros(3,size(P,2))
	Pj[1,:] = P[1,:]
	Pj[3,:] = P[2,:]
	Wj = W

	Pi = [0 0 0 0 0 0 0 0 0;1 1 0 -1 -1 -1 0 1 1;0 1 1 1 0 -1 -1 -1 0]
	Wi = [1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1]
	Pi[:,1:4] = []
	Wi[1:4] = []
	V = [0 0 0 1/2 1/2 1 1 1]

	nPj = size(Pj,2)
	nPi = size(Pi,2)

	xj = Pj[1,:]
	zj = Pj[3,:]

	Px = zeros(nPj,nPi)
	Py = Px
	Pz = Px

	for i = 1:length(xj)
    	Px[i,:] = xj[i].* ones(1,nPi)
    	Py[i,:] = zj[i].* Pi[2,:]
    	Pz[i,:] = zj[i].* Pi[3,:]
	end

	W = zeros(size(Px))
	for i = 1:nPj
    	W[i,:] = Wj[i].*Wi
	end
	nrbplot3(Px,Py,Pz,W,U,V,nptsu,nptsv,flag_plot)
	return Px,Py,Pz,U,V,W
end

end # module
