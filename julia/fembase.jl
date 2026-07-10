using LinearAlgebra, SparseArrays

function elem1D2N( he, re, x, mat, src, gp, w )
# Element matrix for Laplacian, 1D, 2-node linear element - assumes mat & src constant over element
    # he     - element stiffness matrix
    # re     - element right hand side vector
    # x      - 2x1 nodal coordinates
    # mat(x) - material coefficient
    # src(x) - source
    # gp, w - not used
    
    N = [ x[1] 1
          x[2] 1 ] \ [ 1 0
                       0 1 ]
    center= [ 0.5*(x[1]+x[2]) ]
    he = mat(center)*abs(x[1]-x[2])*[ N[1,1]^2      N[1,1]*N[1,2] 
                             N[1,1]*N[1,2]  N[1,2]^2      ]
    re = abs(x[1]-x[2])*src(center)*ones(2,1)
    return nothing
end

function elem2D3N( he, re, x, mat, src, gp, w )
# Element matrix for Laplacian, 2D, 3-node linear element (triangle) - assumes mat & src constant over element
    # he     - element stiffness matrix
    # re     - elrmrnt right hand side vector
    # x      - 3x2 nodal coordinates
    # mat(x,y) - material coefficient
    # src(x,y) - source
    # gp, w - not used
    
	A = [x ones(3,1) ]
	S = det(A)/2
	N = A \ 1I
    for i=1:3
		for j=1:3
			he[i,j] = N[1,i]*N[1,j] + N[2,i]*N[2,j];
		end
	end
    center = [ sum(x[:,1])/3, sum(x[:,2])/3 ]
    he .*= mat(center)*S
    re .= S*src(center)/3
    return nothing
end

function elem3D4N( he, re, x, mat, src, gp, w )
# Element matrix for Laplacian, 3D, 4-node linear element (tetrahedra), constant material and source
    # he     - element stiffness matrix
    # re     - elrmrnt right hand side vector
    # x      - 4×3 nodal coordinates
    # epsfun(x,y,z) - material coefficient
    # rhofun(x,y,z) - source
    # gp, w - not used
    
	A = [ x ones(4,1) ]
	V = abs(det(A)/6)
	N = A \ 1I
	for i=1:4
		for j=1:4
			he[i,j] = N[1,i]*N[1,j] + N[2,i]*N[2,j] + N[3,i]*N[3,j]
		end
	end
    center = [ sum(x[:,1])/4, sum(x[:,2])/4, sum(x[:,3])/4 ]
    me = mat(center)
    se = src(center)
	he .*= me*V
    re .= se*V
    
	return nothing
end

function tetgauss(n::Int)
"""
Gaussin quadratures 
    
    gp, w = tetgauss(n)

Allowed n: 1,4,5,11,15,24

Returns:
    gp :: Matrix{Float64}   nGP × 4   (barycentric coordinates)
    w  :: Vector{Float64}   weights (sum(w)=1)
"""

    gp = Matrix{Float64}(undef, 0, 4)
    w  = Float64[]

    function add!(λs, weight)
        gp_new = reduce(vcat, [reshape(collect(λ), 1, 4) for λ in λs])
        gp = nothing
    end

    function append_points!(pts, ww)
        for p in pts
            gp = nothing
        end
    end

    function orbit4(a,b)
        return [
            a b b b
            b a b b
            b b a b
            b b b a
        ]
    end

    function orbit6(a,b)
        return [
            a a b b
            a b a b
            a b b a
            b a a b
            b a b a
            b b a a
        ]
    end

    function orbit12(a,b,c)
        return [
            a a b c
            a a c b
            a b a c
            a c a b
            b a a c
            c a a b
            a b c a
            a c b a
            b a c a
            c a b a
            b c a a
            c b a a
        ]
    end

    function push_rule!(pts, weight)
        for i in 1:size(pts,1)
            push!(w, 6.0*weight)
        end
        return pts
    end

    blocks = Matrix{Float64}[]

    if n == 1

        push!(blocks, push_rule!([0.25 0.25 0.25 0.25],
                                 0.166666666666666667))

    elseif n == 4

        push!(blocks, push_rule!(orbit4(0.5854101966249685,
                                        0.1381966011250105),
                                 0.0416666666666666667))

    elseif n == 5

        push!(blocks, push_rule!([0.25 0.25 0.25 0.25],
                                 -0.133333333333333333))

        push!(blocks, push_rule!(orbit4(0.5,
                                        0.166666666666666667),
                                 0.075000000000000000))

    elseif n == 11

        push!(blocks, push_rule!([0.25 0.25 0.25 0.25],
                                 -0.0131555555555555556))

        push!(blocks, push_rule!(orbit4(0.785714285714285714,
                                        0.0714285714285714285),
                                 0.00762222222222222222))

        push!(blocks, push_rule!(orbit6(0.399403576166799219,
                                        0.100596423833200785),
                                 0.0248888888888888889))

    elseif n == 15

        push!(blocks, push_rule!([0.25 0.25 0.25 0.25],
                                 0.0302836780970891856))

        push!(blocks, push_rule!(orbit4(0.0,
                                        0.333333333333333333),
                                 0.00602678571428571597))

        push!(blocks, push_rule!(orbit4(0.727272727272727273,
                                        0.0909090909090909091),
                                 0.0116452490860289742))

        push!(blocks, push_rule!(orbit6(0.0665501535736642813,
                                        0.433449846426335728),
                                 0.0109491415613864534))

    elseif n == 24

        push!(blocks, push_rule!(orbit4(0.356191386222544953,
                                        0.214602871259151684),
                                 0.00665379170969464506))

        push!(blocks, push_rule!(orbit4(0.877978124396165982,
                                        0.0406739585346113397),
                                 0.00167953517588677620))

        push!(blocks, push_rule!(orbit4(0.0329863295731730594,
                                        0.322337890142275646),
                                 0.00922619692394239843))

        push!(blocks, push_rule!(orbit12(0.0636610018750175299,
                                         0.269672331458315867,
                                         0.603005664791649076),
                                 0.00803571428571428248))

    else
        error("Unsupported tetrahedral quadrature rule: $n points. Use 1, 4, 5, 11, 15 or 24.")
    end

    gp = vcat(blocks...)

    return gp, w
end

function elem3D4N_num(he, re, x, mat, src, gp, w)
# Element matrix for Laplacian, 3D, 4-node linear element (tetrahedra) - numerical integration
    # he     - element stiffness matrix
    # re     - elrmrnt right hand side vector
    # x      - 4×3 nodal coordinates
    # mat(x,y,z) - material coefficient
    # src(x,y,z) - source
    # gp     - nGP × 4 barycentric coordinates of integration points
    # w      - nGP weights (sum(w) = 1)

    A = [x ones(4,1)]
    V = abs(det(A)/6)   # volume

    # shape functions
    C = A \ I(4)
    
    fill!(he, 0.0)
    fill!(re, 0.0)

    for k in eachindex(w)

        lambda = gp[k,:]

        # global coordinates of the k-th integration point
        xk = vec( lambda' * x )

        eps = mat(xk)
        rho = src(xk)

        # wartości funkcji kształtu
        N = lambda

        for i=1:4

            re[i] += rho*N[i]*w[k]*V

            for j=1:4

                grad =
                    C[1,i]*C[1,j] +
                    C[2,i]*C[2,j] +
                    C[3,i]*C[3,j]

                he[i,j] += eps*grad*w[k]*V

            end
        end
    end

    return nothing
end

function volume3D4N( x )
# 3D, 4-node linear element (tetrahedra) volume
	A = [ x ones(4,1) ]
	return abs(det(A)/6)
end

function mesh3DVolume( mesh )
    evol = [ volume3D4N( mesh["x"][mesh["v"][e,:],:] ) for e=1:size(mesh["v"],1) ]
    return sum(evol),evol
end

function localErr5pQuad3D4N( x, Fv )
# 5-point gaussian quadrature of the square of the linear vector function Fv over tetrahedra
# (x is 4 x 3 matrix of nodal coordinates, Fv is 4 x 3 matric of the values of Fv at the nodes)
    a=0.25
    b=1.0/6.0
    c=0.5
    p = [ [a,a,a], [b,b,b], [b,b,c], [b,c,b], [c,b,b] ]
    w = [ -2.0/15., 3.0/40., 3.0/40., 3.0/40., 3.0/40. ]
    N = [ -1. -1. -1. 1.
           1.  0.  0. 0.
           0.  1.  0. 0.
           0.  0.  1. 0.]
    dNl = [ -1 1 0 0
            -1 0 1 0
            -1 0 0 1 ]
    J = dNl*x
    sum = 0
    for k=1:5
        Fvp = Fv'*(N*[p[1]...;1])
        sum += dot(Fvp,Fvp)*w[k]
    end
    return abs(det(J))/6,sum*6
end

function localErr1pQuad3D4N( x, diff )
# Gradient-recovery error estimator for single tetrahedra
# 1-poin quadrature
	vol= det( [ x ones(4,1) ] )/6
	cdiff = vec(sum(diff,dims=1)/4) 
    # return element volume and squared err
    #println( (pcE-linE4) )
	return abs(vol),dot(cdiff,cdiff)
end

function aposteriori( mesh, V; E=nothing, linE=nothing, elemErr=localErr1pQuad3D4N, elemGrad=gradVIn3D4NElem )
# Gradient-recovery error estimator for tetrahedral mesh
	if E == nothing
		E = grad(mesh["x"],mesh["v"],V,elemGrad)
	end
	if linE == nothing
		linE= linear3Dgrad( mesh, E )
	end
    x= mesh["x"]
    nop= mesh["v"]
    localEE = zeros(size(nop,1),1)
    tot_vol = 0
	for e=1:size(nop,1)
		nds = nop[e,:]
		vol, err  = elemErr(x[nds,:], linE[nds,:].-E[e,:]')
        localEE[e] = vol*err
        tot_vol += vol
	end
	return vec(localEE),sum(localEE)/size(localEE,1),E,linE
end

function gradVIn3D4NElem( x, v )
# Value of the gradient of linear field v i element x
    A = [ x ones(4,1) ]
    N = A \ 1I
    return N[1:3,:]*v
end

function grad( x, nop, V, elemVals )
# Calculate grad of V over (x,nop) mesh
    grad = zeros(size(nop,1),size(x,2))
    for e=1:size(nop,1)
        nds = nop[e,:]
        grad[e,:] = elemVals( x[nds,:], V[nds] )
    end
    return grad
end

function linear3Dgrad( mesh, E)
# Recalculate piecewise-constant => linear gradient over 3D mesh
    x = mesh["x"]
    nop = mesh["v"]
    linG = zeros(size(x,1),size(E,2))
    if !( "iv" in keys(mesh)) || mesh["iv"] == nothing
        makeinop( mesh )
    end
    inop = mesh["iv"]
    for v=1:size(x,1)
        vol= 0
        for e in inop[v]
            nds = nop[e,:]
            evol = volume3D4N( x[nds,:] )
            linG[v,:] += E[e,:] * evol
            vol += evol
        end
        linG[v,:] /= vol
    end
    return linG
end

function apply_periodic_sets!(dofmap::Vector{Int},
                              periodicsets::Vector{<:AbstractSet{Int}})

    # Wszystkie węzły sz każdego zbioru są przenumerowane na minimum
    for s in periodicsets
        isempty(s) && continue

        master = minimum(s)

        for i in s
            dofmap[i] = master
        end
    end

    # 2. Kompresja
    renumber = Dict{Int,Int}()  # słownik do renumeracji
    used = 0

    for i in eachindex(dofmap)
        old = dofmap[i]

        if haskey(renumber, old)  # numer już jest w słowniku
            dofmap[i] = renumber[old]
        else                      # nie ma - dodajemy kolejny nr
            used += 1
            renumber[old] = used
            dofmap[i] = used
        end
    end

    return used    # liczba stopni swobody
end

function stiffnessPattern( ndof, nop, dofmap )
# Sparsity pattern in CSR/CSC format (the matrix is symmetric)
    ne,nen = size(nop)
    ias = [ Set{Int}() for i=1:ndof ]
    for e=1:ne
        for v=1:nen
            for w=v:nen
                push!(ias[dofmap[nop[e,v]]],dofmap[nop[e,w]])
                push!(ias[dofmap[nop[e,w]]],dofmap[nop[e,v]])
            end
        end
    end
    #for v=1:nv
    #    println( v, ":", ias[v] )
    #end
    ia = [ 1 for i=1:ndof+1]
    imx,mx = 0,0
    imi,mi = 0,100
    for v=1:ndof
        ia[v+1] = ia[v] + length(ias[v])
    end
    ja = [ 0 for i=1:(ia[ndof+1]-1) ]
    j= 1
    for v=1:ndof
        pja = sort( [ n for n in ias[v] ] )
        for k=1:size(pja,1)
            ja[j] = pja[k]
            j += 1
        end
    end
    return ia,ja
end

function sparseStructure( ndof, nop, dofmap )
# Generates pattern in CSC/CSR format, then converts it to row-indexes,col-indexes
    ip,ja = stiffnessPattern( ndof, nop, dofmap )
    ia = zeros(Int,size(ja))
    va = zeros(size(ja))
    l= 1
    for i=1:ndof
        for k=ip[i]:(ip[i+1]-1)
            ia[l] = i
            l += 1
        end
    end
    return sparse(ia,ja,va),ip,ja
end

function make_field(v, ne::Int)
# Helper for mes - allows use materials and sources constant over element or variable over domain
    if isa(v, AbstractVector)
        length(v) == ne || error("Vector length must be equal to number of elements")
        return e -> (x -> v[e])

    elseif isa(v, Function)
        return e -> (x -> v(x))

    else
        return e -> (x -> v)
    end
end

function mes( x, nop, mat, src, elem, gp, w, dofmap ) 
# Basic formulation of FEM equation system 
    n = maximum(dofmap)
    H,ip,ja = sparseStructure(n,nop,dofmap)
    R = zeros(n,1)
    ne, nve = size(nop)
    matfield = make_field(mat, ne)
    srcfield = make_field(src, ne)
    he = zeros(nve,nve)
	re = zeros(nve,1)
    nds = [0 for i=1:nve]
    ig,jg = 0,0
    for e=1:ne
		nds .= nop[e,:]
        me = matfield(e)
        se = srcfield(e)
        elem( he, re, x[nds,:], me, se, gp, w )
        for i=1:nve
            ig = dofmap[nds[i]]
            for j=1:nve
                jg = dofmap[nds[j]]
                H[ig,jg] += he[i,j]
            end
            R[ig] += re[i]
        end
    end
    return H,R,ip,ja
end

#=
LEGACY
function vmes( x, nop, mat, src, elem ) 
# Basic formulation of FEM equation system VECTORIZED version
    n = size(x,1)
    #H = spzeros(n,n)
    @time H,ip,ja = sparseStructure(n,nop)
    R = zeros(n,1)
    ne, nve = size(nop)
    A = ones(nve,nve)
    he = zeros(nve,nve)
	re = zeros(nve,1)
    nds = [0 for i=1:nve]
    for e=1:ne
        nds .= nop[e,:]
        elem( he, re, x[nds,:], mat[e], src[e] )
        H[nds[:],nds[:]] .+= he
        R[nds[:]] .+= re
    end
    return H,R,ip,ja
end
=#

function dbc( H, R, dofmap, nodes, vals, ip, ja )
# Introduce Dirichlet BC to (H,R)
    n = 0
    f = 0
	for i=1:length(nodes)
		n = dofmap[nodes[i]]
		for k=ip[n]:ip[n+1]-1
            if ja[k] != n
                H[n,ja[k]] = 0
                R[ja[k]] -= H[ja[k],n]*vals[i]
                H[ja[k],n] = 0
            else
                f = H[n,ja[k]]
            end
        end
		R[n] = f*vals[i]
	end
	return H,R
end

function solveBP( mesh, materials, sources, diribc; elem=elem3D4N, gp=nothing, w=nothing, solver="base", dofmap=nothing )
    nop = mesh["v"]
    x= mesh["x"]
    subdomain= mesh["s"]
    ne = size(nop,1)
    if dofmap == nothing
        dofmap = 1:size(x,1)
    end
    mat = materials
    if ! applicable( materials, x[1,:])
       mat= [materials[subdomain[e]] for e in 1:ne]
    end
    src = sources
    if ! applicable( sources, x[1,:] )
        src = [sources[subdomain[e]] for e in 1:ne]
    end
    #@time H,R = vmes(x,nop,eps,rho,elem) # seems to be many times slower
    @time H,R,ip,ja = mes(x,nop,mat,src,elem, gp, w, dofmap) 
    H,R = dbc(H,R,dofmap,diribc["bn"],diribc["bv"],ip,ja)
    for i = 1:size(H,1)
        if abs(H[i,i]) < 1e-15
            println( "Singular matrix for row", i, " |diagonal| < 1e-15" )
        end
    end
    sH = deepcopy(H)
    @time V = H \ R
    fullV = V[dofmap]
    return fullV,H,sH,R
end

#=
LEGACY
using AlgebraicMultigrid
function AMGsolveBP( mesh, materials, sources, diribc, elem=elem3D4N )
    nop = mesh["v"]
    x= mesh["x"]
    mat= mesh["s"]
    ne = size(nop,1)
    vn = size(x,1)
    rho = zeros(ne,1)
    eps = ones(ne,1)
    for e = 1:ne
        eps[e] = materials[mat[e]]
        rho[e] = sources[mat[e]]
    end
    @time H,R,ip,ja = mes(x,nop,eps,rho,elem) 
    H,R = dbc(H,R,diribc["bn"],diribc["bv"],ip,ja)
    @time V = solve(H, vec(R), RugeStubenAMG(), maxiter = 1, abstol = 1e-12)
    return V,H,R
end
=#