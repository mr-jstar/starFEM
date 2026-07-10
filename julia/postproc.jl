function in3DElem( x, p )
# Check if point p is inside of tetrahedra x
    for i=1:3
        if p[i] > maximum(x[:,i]) || p[i] < minimum(x[:,i])
            return false,0
        end
    end
    EPS = 6e-15
    A = [x ones(4,1) ]
    N = A \ 1I
    f=[p
    1]'*N
    low = minimum(f)
    high = maximum(f)
    err = maximum( [ 0-low, high-1 ] )
    return low >= -EPS && high <= 1.0+EPS, err
end

function locate3DElem( x, nop, p; subset=Set() )
# Locate element of tetrahedral mesh containing point p (test only subset of elelements)
    if isempty(subset)
        subset=1:size(nop,1)
    end
    for e=subset
        test,err = in3DElem( x[nop[e,:],:], p )
        if test
            return e
        end
    end   
    return 0     
end

using StaticArrays

function in3DElem_fast(X, T, e, p; EPS=6e-15, debug=false)
    i1 = T[e,1]
    i2 = T[e,2]
    i3 = T[e,3]
    i4 = T[e,4]

    # szybkie testy bbox 
    xmin1 = min(X[i1,1], X[i2,1], X[i3,1], X[i4,1])
    xmax1 = max(X[i1,1], X[i2,1], X[i3,1], X[i4,1])
    if p[1] < xmin1 || p[1] > xmax1
        return false, 0.0
    end

    xmin2 = min(X[i1,2], X[i2,2], X[i3,2], X[i4,2])
    xmax2 = max(X[i1,2], X[i2,2], X[i3,2], X[i4,2])
    if p[2] < xmin2 || p[2] > xmax2
        return false, 0.0
    end

    xmin3 = min(X[i1,3], X[i2,3], X[i3,3], X[i4,3])
    xmax3 = max(X[i1,3], X[i2,3], X[i3,3], X[i4,3])
    if p[3] < xmin3 || p[3] > xmax3
        return false, 0.0
    end

    x1 = X[i1,1]; y1 = X[i1,2]; z1 = X[i1,3]
    x2 = X[i2,1]; y2 = X[i2,2]; z2 = X[i2,3]
    x3 = X[i3,1]; y3 = X[i3,2]; z3 = X[i3,3]
    x4 = X[i4,1]; y4 = X[i4,2]; z4 = X[i4,3]

    # wyznaczamy wsp barycentryczne
    # p = x4 + l1*(x1-x4) + l2*(x2-x4) + l3*(x3-x4)
    A = @SMatrix [
        x1-x4  x2-x4  x3-x4
        y1-y4  y2-y4  y3-y4
        z1-z4  z2-z4  z3-z4
    ]

    b = @SVector [p[1]-x4, p[2]-y4, p[3]-z4]

    bc = A \ b

    l1 = bc[1]
    l2 = bc[2]
    l3 = bc[3]
    l4 = 1 - l1 - l2 - l3

    if debug
        @show e l1 l2 l3 l4
    end
    
    low = min(l1,l2,l3,l4)
    high = max(l1,l2,l3,l4)

    err = max(0.0 - low, high - 1.0)

    return low >= -EPS && high <= 1.0 + EPS, err
end

function locate3DElem_fast(X, T, p; subset=nothing)
    elems = subset === nothing ? axes(T,1) : subset

    for e in elems
        inside, err = in3DElem_fast(X, T, e, p)
        if inside
            return e
        end
    end

    return 0
end

function valueAt3DPoint( x, nop, V, p, elemVal )
# Calculate vaule of linear potential V at point p in 3D mesh
    #println(p)
    for e=1:size(nop,1)
        test,err = in3DElem_fast( x, nop, e, p )
        if test
            return elemVal(x[nop[e,:],:],V[nop[e,:]],p)
        end
    end   
    return NaN  
end

function cut3D(x, nop, V, d, dval, nn, elemVal )
# cuts potential V at dval of dimension d
    nx = ny  = 10
    if nn > 0
        nx = ny = nn
    end
    i= 0
    box= zeros(2,2)
    dms = [0,0,0]
    for dim=1:3
        if dim != d
            i+=1
            dms[i]= dim
            box[i,1] = minimum(x[:,dim])+0.0001
            box[i,2] = maximum(x[:,dim])-0.0001
        end
    end
    dms[3] = d
    x2= LinRange(box[1,1], box[1,2], nx )
    y2= LinRange(box[2,1], box[2,2], ny )
    v2 = zeros(nx,ny)
    p=zeros(3,1)
    p[dms[3]]= dval
    for i=1:nx
        p[dms[1]] = x2[i]
        for j=1:ny
            p[dms[2]] = y2[j]
            v2[i,j]= valueAt3DPoint( x, nop, V, p, elemVal)
        end
    end
    return x2,y2,v2
end

function scan3D(x, nop, V, nn, elemVal )
# 3D scan of linear potential V on rectangular grid nn x nn x nn
    nx = ny = nz = 10
    if nn > 0
        nx = ny = nz = nn
    end
    box= zeros(3,2)
    for dim=1:3
       box[dim,1] = minimum(x[:,dim])+0.001
       box[dim,2] = maximum(x[:,dim])-0.001
    end
    x2= LinRange(box[1,1], box[1,2], nx )
    y2= LinRange(box[2,1], box[2,2], ny )
    z2= LinRange(box[3,1], box[3,2], nz )
    v2 = zeros(nx,ny,nz)
    p=zeros(3,1)
    for i=1:nx
        p[1] = x2[i]
        for j=1:ny
            p[2] = y2[j]
            for k=1:nz
                p[3] = z2[k]
                v2[i,j,k]= valueAt3DPoint( x, nop, V, p, elemVal)
            end
        end
    end
    return x2,y2,z2,v2
end

function cut3DPiecewiseConst(x, nop, C, d, dval, nn)
# cuts piecewise constant field C at dval of dimension d
    nx = ny  = 10
    if nn > 0
        nx = ny = nn
    end
    i= 0
    box= zeros(2,2)
    dms = [0,0,0]
    for dim=1:3
        if dim != d
            i+=1
            dms[i]= dim
            box[i,1] = minimum(x[:,dim])+0.0001
            box[i,2] = maximum(x[:,dim])-0.0001
        end
    end
    dms[3] = d
    x2= LinRange(box[1,1], box[1,2], nx )
    y2= LinRange(box[2,1], box[2,2], ny )
    v2 = zeros(nx,ny)
    p=zeros(3,1)
    p[dms[3]]= dval
    for i=1:nx
        p[dms[1]] = x2[i]
        for j=1:ny
            p[dms[2]] = y2[j]
            e=locate3DElem( x, nop, p )
            if e > 0 
                v2[i,j]= C[e]
            end
        end
    end
    return x2,y2,v2
end

function scan3DPiecewiseConst(x, nop, C, nn, dxdydz=[1,1,1] )
# Fast 3D scan of piecewise constant field C on rectangular grid 
    box= zeros(3,2)
    for dim=1:3
       box[dim,1] = minimum(x[:,dim])+0.001
       box[dim,2] = maximum(x[:,dim])-0.001
    end
    nx = ny = nz = 10
    if nn > 0
        nx = ny = nz = nn
        x2= LinRange(box[1,1], box[1,2], nx )
        y2= LinRange(box[2,1], box[2,2], ny )
        z2= LinRange(box[3,1], box[3,2], nz )
        dxdydx = [ (box[1,2]-box[1,1])/(nx-1), (box[2,2]-box[2,1])/(ny-1), (box[3,2]-box[3,1])/(nz-1) ]
    else
        x2= box[1,1]:dxdydz[1]:box[1,2]
        nx= size(x2)[1]
        y2= box[2,1]:dxdydz[2]:box[2,2]
        ny= size(y2)[1]
        z2= box[3,1]:dxdydz[3]:box[3,2]
        nz= size(z2)[1]
    end
    
    ex = zeros(Int64,nx,ny,nz)
    errs = zeros(Float64,nx,ny,nz)
    for e=1:size(nop,1)
        v = x[nop[e,:],:]
        ix=1+ceil(Int64,minimum(v[:,1])/dxdydz[1])
        jx=minimum([nx,1+floor(Int64,maximum(v[:,1])/dxdydz[1])])
        iy=1+ceil(Int64,minimum(v[:,2])/dxdydz[2])
        jy=minimum([ny,1+floor(Int64,maximum(v[:,2])/dxdydz[2])])
        iz=1+ceil(Int64,minimum(v[:,3])/dxdydz[3])
        jz=minimum([nz,1+floor(Int64,maximum(v[:,3])/dxdydz[3])])
        if ix < 1 || iy < 1 || iz < 1
            println( e, " ", ix, " ", iy, " ", iz )
        end
        for i=ix:jx
            for j=iy:jy
                for k=iz:jz
                    test,err = in3DElem( v, [x2[i],y2[j],z2[k]] )
                    if test
                        if  ex[i,j,k] == 0 || ( ex[i,j,k] > 0 && errs[i,j,k] > err )
                            ex[i,j,k] = e
                            errs[i,j,k] = err
                        end
                    end
                end
            end
         end
    end

    v2 = zeros(nx,ny,nz)
    
    p=zeros(3,1)
    for i=1:nx
        p[1] = x2[i]
        for j=1:ny
            p[2] = y2[j]
            for k=1:nz
                p[3] = z2[k]
                if ex[i,j,k] > 0
                    v2[i,j,k] = C[ex[i,j,k]]
                else
                    #println( "no elem for ", i, ",", j, ",", k )
                    e=locate3DElem( x, nop, p )
                    if e > 0 
                        v2[i,j,k]= C[e]
                        ex[i,j,k]= e
                        errs[i,j,k]= 1.0
                    end
                end
            end
        end
    end
    return x2,y2,z2,v2,ex,errs
end

function valueIn3D4NElem( x, v, p )
# Compute value of linear field v at point p (p must lie within element x)
    A = [ x ones(4,1) ]
    N = A \ 1I
    return dot([p;1]'*N,v)
end

function flux( mesh, mats, E, s )
# Compute flux of mat*E through the surface s
    if !("faces" in keys(mesh)) || mesh["faces"] == nothing || findfirst(x->x==s,msh["faces"][:,1]) == nothing
        return nothing
    else
        if !("iv" in keys(mesh)) || mesh["iv"] == nothing
            makeinop( mesh )
        end
        flx=0.0
        tot_a=0.0
        for f in eachrow(mesh["faces"])
            if f[1] == s
                nds = f[2:4]
                es = [e for e in intersect(mesh["iv"][nds]...)]
                uv = mesh["x"][nds[1:2],:] .- mesh["x"][nds[3],:]'
                n=cross(uv[1,:],uv[2,:])
                a=0.5*norm(n)
                n /= norm(n)
                mE = E[es[1],:]*mats[mesh["s"][es[1]]]
                for i=2:length(es)
                    mE .+= E[es[i],:]*mats[mesh["s"][es[i]]]
                end
                mE /= length(es)
                #println( f, ":", nds, " e:", es, " n=", n, " mE=", mE, " a=", a)
                flx += abs(a*dot(n,mE))
                tot_a += a
            end
        end
        return flx,tot_a
    end
end

function cutTetra(x, V, dim, at; tol=1e-10)
# Cut tetrahedra x at x(:,dim)=at, V - linear field over tetrahedra
    a = x[1,dim]
    b = x[2,dim]
    c = x[3,dim]
    d = x[4,dim]

    if !(min(a,b,c,d) - tol <= at <= max(a,b,c,d) + tol)
        return 0, nothing
    end

    edges = ((1,2), (1,3), (1,4), (2,3), (2,4), (3,4))

    nv = size(V,2)
    cut = zeros(Float64, 4, 3 + nv)
    nc = 0

    function addpoint!(idx)
        px = x[idx,1]
        py = x[idx,2]
        pz = x[idx,3]

        for k in 1:nc # check for duplicates
            if abs(cut[k,1]-px) < tol &&
               abs(cut[k,2]-py) < tol &&
               abs(cut[k,3]-pz) < tol
                return
            end
        end

        nc += 1
        cut[nc,1] = px
        cut[nc,2] = py
        cut[nc,3] = pz

        for q in 1:nv
            cut[nc,3+q] = V[idx,q]
        end

        return
    end

    function addinterp!(i,j,t)
        px = x[i,1] + t*(x[j,1] - x[i,1])
        py = x[i,2] + t*(x[j,2] - x[i,2])
        pz = x[i,3] + t*(x[j,3] - x[i,3])

        for k in 1:nc  # check for duplicates
            if abs(cut[k,1]-px) < tol &&
               abs(cut[k,2]-py) < tol &&
               abs(cut[k,3]-pz) < tol
                return
            end
        end

        nc += 1
        cut[nc,1] = px
        cut[nc,2] = py
        cut[nc,3] = pz

        for q in 1:nv
            cut[nc,3+q] = V[i,q] + t*(V[j,q] - V[i,q])
        end

        return
    end

    for (i,j) in edges
        si = x[i,dim] - at
        sj = x[j,dim] - at

        if abs(si) <= tol && abs(sj) <= tol
            addpoint!(i)
            addpoint!(j)
        elseif abs(si) <= tol
            addpoint!(i)
        elseif abs(sj) <= tol
            addpoint!(j)
        elseif si * sj < 0
            t = -si / (sj - si)
            addinterp!(i,j,t)
        end
    end

    if nc < 3
        return 0, nothing
    end

    return nc - 2, cut
end

function multiCutMeshSLOW( msh, V, cuts; sub=nothing )
# Make n cuts of the field V spanned over msh
# Orthogonal cut in dimensions dims at values in levels
# If sub is not nothong, consider only subdomains listed in it.
    cx= [0.0 0.0 0.0]
    cv= [0 0 0 ]
    cV= zeros(1,size(V,2))
    nc= 0
    nn= 0
    if sub == nothing
        elems = 1:size(msh["v"],1)
    else
        elems = findall( x->x in sub, msh["s"] )
    end
    for i=1:length(cuts)
        dim,at = cuts[i]
        for e in elems
            #println(e)
            nds = vec(msh["v"][e,:])
            x= msh["x"][nds,:]
            if minimum(@view x[:, dim]) <= at <= maximum(@view x[:, dim])
                ntr,c= cutTetra( x, V[nds,:], dim, at )
                if ntr > 0
                    cx = vcat( cx, c[1:(ntr+2),1:3] )
                    cV = vcat( cV, c[1:(ntr+2),4:end] )
                    if ntr == 1
                        cv = vcat( cv, [ nn+1 nn+2 nn+3 ] )
                    end
                    if ntr == 2
                        cv = vcat( cv, [ nn+1 nn+2 nn+3; nn+2 nn+3 nn+4 ] )
                    end
                    if ntr == 3
                        cv = vcat( cv, [ nn+1 nn+2 nn+3; nn+2 nn+3 nn+4; nn+3 nn+4 nn+5 ] )
                    end
                    nc += 1
                    nn += ntr+2
                end
            end
            #println("---")
        end
    end
    return Dict( "x"=>cx[2:end,:], "v"=>cv[2:end,:]),cV[2:end,:]
end

function multiCutMesh(msh, V, cuts; sub=nothing)
# Make a n-cuts of the field V spanned over msh
# Orthogonal cut in dimensions dims at values in levels
# If sub is not nothing, consider only subdomains listed in it.
    elems = if sub === nothing
        1:size(msh["v"], 1)
    else
        findall(s -> s in sub, msh["s"])
    end
    X = msh["x"]
    T = msh["v"]
    cx_rows = Vector{Vector{Float64}}()
    ce_rows = Vector{NTuple{3,Int}}()
    cV_rows = Vector{Vector{Float64}}()
    result = Vector{Tuple{Dict{String,Any}, Matrix{Float64}}}()
    for cut in cuts
        dim, at = cut
        nn = 0
        empty!(cx_rows)
        empty!(ce_rows)
        empty!(cV_rows)
        for e in elems
            i1 = T[e,1]
            i2 = T[e,2]
            i3 = T[e,3]
            i4 = T[e,4]
            a = X[i1, dim]
            b = X[i2, dim]
            c = X[i3, dim]
            d = X[i4, dim]
            if min(a,b,c,d) <= at <= max(a,b,c,d)
                nds = [i1, i2, i3, i4]
                eX = X[nds, :]
                eV = V[nds, :]
                ntr, ctab = cutTetra(eX, eV, dim, at)
                if ntr > 0
                    np = ntr + 2
                    for j in 1:np
                        push!(cx_rows, collect(ctab[j, 1:3]))
                        push!(cV_rows, collect(ctab[j, 4:end]))
                    end
                    if ntr >= 1
                        push!(ce_rows, (nn+1, nn+2, nn+3))
                    end
                    if ntr >= 2
                        push!(ce_rows, (nn+2, nn+3, nn+4))
                    end
                    if ntr >= 3
                        push!(ce_rows, (nn+3, nn+4, nn+5))
                    end
                    nn += np
                end
            end
        end
        if isempty(cx_rows)
            cx = Matrix{Float64}(undef, 0, 3)
            ce = Matrix{Int}(undef, 0, 3)
            cV = Matrix{Float64}(undef, 0, size(V,2))
        else
            cx = reduce(vcat, permutedims.(cx_rows))
            ce = reduce(vcat, permutedims.(collect.(ce_rows)))
            cV = reduce(vcat, permutedims.(cV_rows))
        end
        push!(result, (Dict("x" => cx, "v" => ce),cV) )
    end

    return result
end

function linScan3D(msh, V, p1, p2, np=50)
    X = msh["x"]
    T = msh["v"]

    xmin = min.(p1, p2)
    xmax = max.(p1, p2)
    bbox = (xmin, xmax)

    elems = findall(e -> boBoxIntersect(elemBoBox3D4N(msh, e), bbox), axes(T,1))

    result = Matrix{Float64}(undef, np, 2)

    dl = (p2 - p1) / (np - 1)
    len_dl = norm(dl)

    p = similar(p1)

    last_e = 0

    for i in 1:np
        @. p = p1 + (i-1)*dl

        # najpierw spróbuj czy poprzedni element może być
        e = 0
        if last_e > 0 && in3DElem_fast( X, T, last_e, p )[1]
            e = last_e
        else # jeśli nie - szukaj
            e = locate3DElem(X, T, p, subset=elems)
        end

        result[i,1] = (i-1) * len_dl

        if e > 0
            result[i,2] = valueAt3DPoint(X, T, V, p, valueIn3D4NElem)
            last_e = e
        else
            result[i,2] = NaN
        end
    end

    return result
end

using PyPlot
using Plots

function plotMesh( x, nop )
# 2D mesh plot
	f=Plots.plot( x[ :,1], x[ :,2], seriestype=:scatter, label="" , aspect_ratio=:equal )

	nds = append!(nop[1,:], nop[1,1])
	xp = x[ nds, 1]
	yp = x[ nds, 2]
	for i=2:size(nop,1)
        nds = append!(nop[i,:],nop[i,1])
        xp = [xp x[nds,1]]
        yp = [yp x[nds,2]]
	end
 	Plots.plot!( xp, yp, label="", aspect_ratio=:equal )

	return f
end

using Makie

function mesh_contour!(
    ax,
    x::AbstractMatrix,
    tri::AbstractMatrix{<:Integer},
    u::AbstractMatrix;
    levels = range(extrema(u)..., length=11),
    at = (0.0, 0.0, 0.0),
    color = :black,
    linewidth = 1
)
    ax0, ay0, az0 = at

    segments = Vector{Point3f}()
    sizehint!(segments, 2 * size(tri, 1))

    for lev in levels
        empty!(segments)

        for e in axes(tri, 1)
            i1 = tri[e,1]
            i2 = tri[e,2]
            i3 = tri[e,3]

            x1 = x[i1,1]; y1 = x[i1,2]; z1 = x[i1,3]; u1 = u[i1]
            x2 = x[i2,1]; y2 = x[i2,2]; z2 = x[i2,3]; u2 = u[i2]
            x3 = x[i3,1]; y3 = x[i3,2]; z3 = x[i3,3]; u3 = u[i3]

            n = 0
            p1 = Point3f(0,0,0)
            p2 = Point3f(0,0,0)

            # edge 1-2
            s1 = lev - u1
            s2 = lev - u2
            if s1 * s2 < 0
                t = s1 / (u2 - u1)
                p = Point3f(
                    x1 + t*(x2-x1) + ax0,
                    y1 + t*(y2-y1) + ay0,
                    z1 + t*(z2-z1) + az0
                )
                n += 1
                if n == 1
                    p1 = p
                else
                    p2 = p
                end
            end

            # edge 2-3
            s2 = lev - u2
            s3 = lev - u3
            if s2 * s3 < 0
                t = s2 / (u3 - u2)
                p = Point3f(
                    x2 + t*(x3-x2) + ax0,
                    y2 + t*(y3-y2) + ay0,
                    z2 + t*(z3-z2) + az0
                )
                n += 1
                if n == 1
                    p1 = p
                else
                    p2 = p
                end
            end

            # edge 3-1
            s3 = lev - u3
            s1 = lev - u1
            if s3 * s1 < 0
                t = s3 / (u1 - u3)
                p = Point3f(
                    x3 + t*(x1-x3) + ax0,
                    y3 + t*(y1-y3) + ay0,
                    z3 + t*(z1-z3) + az0
                )
                n += 1
                if n == 1
                    p1 = p
                else
                    p2 = p
                end
            end

            if n == 2
                push!(segments, p1)
                push!(segments, p2)
            end
        end

        if !isempty(segments)
            linesegments!(
                ax,
                copy(segments);
                color = color,
                linewidth = linewidth
            )
        end
    end

    return nothing
end

############ LEGACY

function slow_linScan3D( msh, V, p1, p2, np=50 )
# Linear scan on segment <p1-p2> of the filed V spanned over mesh msh.
    m = [ p1 p2 ]
    bbox = (minimum(m,dims=2),maximum(m,dims=2))
    elems = findall( x->boBoxIntersect(elemBoBox3D4N(msh,x),bbox), 1:size(msh["v"],1) )
    #    println( size(elems) )
    result = zeros( np, 2 )
    dl = (p2-p1)/(np-1)
    for i=1:np
        p = p1 + (i-1)*dl
        e = locate3DElem( msh["x"], msh["v"], p, subset=elems )
        #=
        if e <= 0
            e = locate3DElem( msh["x"], msh["v"], p )
            println( "Why point ",p, " has not been found?" )
        end
        =#
        if e > 0
            result[i,:] = [ norm(p-p1) valueAt3DPoint( msh["x"], msh["v"], V, p, valueIn3D4NElem ) ]
        else
            result[i,:] = [ norm(p-p1) NaN ]
        end
    end
    return result   
end

function slow_mesh_contour!(
    ax,
    x::AbstractMatrix,
    tri::AbstractMatrix{<:Integer},
    u::AbstractMatrix;
    levels = range(minimum(u), maximum(u), length=11),
    at = (0.0, 0.0, 0.0),
    color = :black,
    linewidth = 1
)

    ax0, ay0, az0 = at

    for e in axes(tri,1)

        i1 = tri[e,1]
        i2 = tri[e,2]
        i3 = tri[e,3]

        x1 = x[i1,1]; y1 = x[i1,2]; z1 = x[i1,3]; u1 = u[i1]
        x2 = x[i2,1]; y2 = x[i2,2]; z2 = x[i2,3]; u2 = u[i2]
        x3 = x[i3,1]; y3 = x[i3,2]; z3 = x[i3,3]; u3 = u[i3]

        for lev in levels

            n = 0
            p1 = Point3f(0,0,0)
            p2 = Point3f(0,0,0)

            # ----- edge 1-2 -----
            if (lev-u1)*(lev-u2) < 0
                t = (lev-u1)/(u2-u1)
                p = Point3f(
                    x1 + t*(x2-x1) + ax0,
                    y1 + t*(y2-y1) + ay0,
                    z1 + t*(z2-z1) + az0
                )
                n += 1
                if n == 1
                    p1 = p
                else
                    p2 = p
                end
            end

            # ----- edge 2-3 -----
            if (lev-u2)*(lev-u3) < 0
                t = (lev-u2)/(u3-u2)
                p = Point3f(
                    x2 + t*(x3-x2) + ax0,
                    y2 + t*(y3-y2) + ay0,
                    z2 + t*(z3-z2) + az0
                )
                n += 1
                if n == 1
                    p1 = p
                else
                    p2 = p
                end
            end

            # ----- edge 3-1 -----
            if (lev-u3)*(lev-u1) < 0
                t = (lev-u3)/(u1-u3)
                p = Point3f(
                    x3 + t*(x1-x3) + ax0,
                    y3 + t*(y1-y3) + ay0,
                    z3 + t*(z1-z3) + az0
                )
                n += 1
                if n == 1
                    p1 = p
                else
                    p2 = p
                end
            end

            if n == 2
                lines!(ax, [p1,p2], color=color, linewidth=linewidth)
            end
        end
    end

    return nothing
end

function slow_Scan3DCnst(x, nop, C, nn, dxdydz=[1,1,1] )
# LEGACY!!!  
# 3D scan of piecewise constant field C on rectangular grid (legacy version, very slow)
    box= zeros(3,2)
    for dim=1:3
       box[dim,1] = minimum(x[:,dim])+0.001
       box[dim,2] = maximum(x[:,dim])-0.001
    end
    nx = ny = nz = 10
    if nn > 0
        nx = ny = nz = nn
        x2= LinRange(box[1,1], box[1,2], nx )
        y2= LinRange(box[2,1], box[2,2], ny )
        z2= LinRange(box[3,1], box[3,2], nz )
        dxdydx = [ (box[1,2]-box[1,1])/(nx-1), (box[2,2]-box[2,1])/(ny-1), (box[3,2]-box[3,1])/(nz-1) ]
    else
        x2= box[1,1]:dxdydz[1]:box[1,2]
        nx= size(x2)[1]
        y2= box[2,1]:dxdydz[2]:box[2,2]
        ny= size(y2)[1]
        z2= box[3,1]:dxdydz[3]:box[3,2]
        nz= size(z2)[1]
    end
    
    buckets=Dict(i => Set() for i = 1:nz)
    for e=1:size(nop,1)
        nds = nop[e,:]
        vz = x[nds,3]
        i=1+ceil(Int64,minimum(vz)/dxdydz[3])
        j=minimum([nz,1+floor(Int64,maximum(vz)/dxdydz[3])])
        for k=i:j
            push!(buckets[k], e)
        end
    end
    suma = 0
    for (key, val) in buckets
        #println(key, ": ", length(val))
        suma += length(val)
    end
    #println( size(nop,1), " -> ", suma )

    v2 = zeros(nx,ny,nz)
    i2 = zeros(nx,ny,nz)
    p=zeros(3,1)
    for i=1:nx
        p[1] = x2[i]
        for j=1:ny
            p[2] = y2[j]
            for k=1:nz
                p[3] = z2[k]
                e=locate3DElem( x, nop, p, subset=buckets[k] )
                if e > 0 
                    v2[i,j,k]= C[e]
                    i2[i,j,k]= e
                end
            end
        end
    end
    return x2,y2,z2,v2,i2
end
