function uziom(file)
    mesh=loadNG(file)
    V,rho,eps,bn,bv = solveBox( mesh )
    return mesh,V,rho,eps,bn,bv
end

function uziom_adaptive(file,maxit)
    mesh=loadNG(file)
    err = zeros(maxit,3)
    V = zeros(size(mesh["x"],1))
    bp = Dict()
    for it=1:maxit
        V,bp = solveBox( mesh )
        E = grad( mesh["x"], mesh["v"],V,gradVIn3D4NElem)
        localEE,meanEE,_,_= aposteriori( mesh, V, E=E)
        maxEE= maximum(localEE)
        err[it,:]=[size(mesh["x"],1),meanEE,maxEE]
        if it < maxit
            refineMesh( mesh, localEE )
        end
    end
    return mesh,V,err,bp
end

function solveBox( mesh ) 
    nop = mesh["v"]
    x= mesh["x"]
    mat= mesh["s"]
    ne = size(nop,1)
    vn = size(x,1)
    rho=zeros(ne,1)
    eps=ones(ne,1)
    for e=1:ne
        if mat[e] > 0
            eps[e] = 100
        end
    end
    H,R=mes(x,nop,eps,rho,elem3D4N)
    height = maximum(x[:,3])
    bn=[]
    bv=[]
    for v=1:size(x,1)
        if  x[v,3]==0
            push!(bn,v)
            push!(bv,0)
        end
        if  x[v,3] > 0.9999*height
            push!(bn,v)
            push!(bv,height)
        end
    end
    H,R=dbc(H,R,bn,bv)
    V= H\R
    return V,Dict("rho"=>rho,"eps"=>eps,"dbn"=>bn,"dbv"=>bv)
end

function hotspots(x,nop,V,threshold,gradf)
    ne = size(nop,1)
    E=zeros(ne,3)
    for e=1:ne
        nds=nop[e,:]
        E[e,:] = gradf(x[nds,:],V[nds])'
    end
    mE = sqrt.(sum(E.^2,dims=2))
    Emax = maximum(mE)
    #println( Emax )
    Et = threshold*Emax
    #println( Et )
    cld = []
    for e=1:ne
        nds=nop[e,:]
        xx = x[nds,1]
        yy = x[nds,2]
        if mE[e] >= Et
            push!(cld,[[minimum(xx),maximum(xx),maximum(xx),minimum(xx),minimum(xx)],[minimum(yy),minimum(yy),maximum(yy),maximum(yy),minimum(yy)]] )
        end
    end
    return cld,E,mE
end
