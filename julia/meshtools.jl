function quality3D4N( x )
# Assess quality of tetra x (4 by 3 matrix of vertex coordinates) as a ratio of the volume related to the volume of circumscribed sphere.
# The Vtet/Vsphere is rescaled, so a regular tetra shall get Q=1.
    d= x[2:4,:] .- x[1,:]'
    c= d \ sum(abs2.(d),dims=2)
    r= sqrt(dot(c,c))/2
    #println( c, " ", r, " ", det(d)/6)
    V = abs(det(d)/6)
    return (9*sqrt(3)/8*V)^(1/3)/r
end

function meshQuality( mesh )
# Assess quality of tetras in the mesh x (4 by 3 matrix of vertex coordinates) as a ratio of the volume related to the volume of 
    quals = zeros(size(mesh["v"],1),1)
    for e=1:size(mesh["v"],1)
        quals[e] = quality3D4N( mesh["x"][mesh["v"][e,:],:] )
		if quals[e] < 0 
			println( "Warning: volume of element # $e < 0 !!!" )
		end
    end
    return quals
end

function splitTetraCentral( mesh, e )
# Split e-th tetrahedra into four (add center of grawity as a new vertex)
    one = size(mesh["v"],1)
    onv = size(mesh["x"],1)
    nds = mesh["v"][e,:]
    mesh["x"]= vcat( mesh["x"], (sum(mesh["x"][nds,:],dims=1)/4) )
    if "b" in keys(mesh) && mesh["b"] != nothing
        mesh["b"] = vcat( mesh["b"], [0] )
    end
    nv = size(mesh["x"],1)
    mesh["v"][e,:] = [nds[1:3];nv]'
    mesh["v"]= vcat( mesh["v"], [nds[2];nds[4];nds[3];nv]')
    mesh["v"]= vcat( mesh["v"], [nds[3];nds[4];nds[1];nv]')
    mesh["v"]= vcat( mesh["v"], [nds[1];nds[4];nds[2];nv]')
    mesh["s"]= vcat( mesh["s"], [mesh["s"][e],mesh["s"][e],mesh["s"][e]] )
    if "iv" in keys(mesh)
         delete!(mesh,"iv")
    end
    #println( onv, "->", size(mesh["x"],1), "  ", one, "->", size(mesh["v"],1) )
end

function moveVertexToCentroid( mesh, v )
# Move vertex v to the conter of the surrounding hull
    neigh= Set()
    ovol = 0.0
    for e in mesh["iv"][v]
        union!(neigh,mesh["v"][e,:])
        ovol += volume3D4N( mesh["x"][mesh["v"][e,:],:] )
    end
    oldx = mesh["x"][v,:]
    mesh["x"][v,:] = vec(sum(mesh["x"][collect(neigh),:],dims=1))/length(neigh)
    nvol = 0.0
    negv = (0,0.0)
    for e in mesh["iv"][v]
        evol = volume3D4N( mesh["x"][mesh["v"][e,:],:] )
        if evol < 0
            negv = (e,evol)
        end
        nvol += evol
    end
    if abs(ovol-nvol) > 0.001*ovol || negv[1] != 0
        #println( "Vertex ", v, ": ", oldx, "=>", mesh["x"][v,:], " : volume changed: ", ovol, "->", nvol )
        #for e in mesh["iv"][v]
        #    println( "V(element ", e, ")= ", volume3D4N( mesh["x"][mesh["v"][e,:],:] ) )
        #end
        mesh["x"][v,:] = oldx;
        #println( "Move withdrawed" )
    #else
    #   println( v, " ok" )
    end
end

function moveToCentroids( mesh )
# Move all free vertices of the mesh towards the centroids
    if !("iv" in keys(mesh)) || mesh["iv"] == nothing
        makeinop( mesh )
    end
    if !("b" in keys(mesh)) || mesh["b"] == nothing
        mesh["b"] = zeros(size(mesh["x"],1),1)
    end
    newx = zeros(size(mesh["x"]))
    for v=1:size(mesh["x"],1)
        if mesh["b"][v] == 0
            moveVertexToCentroid( mesh, v )
        end
    end
end

function elemCG3D4N(msh, e )
# Center of gravity for element e
    return sum(msh["x"][msh["v"][e,:],:],dims=1)/size(msh["v"][e,:],1)
end

function elemBoBox3D4N(msh, e )
# Bounding Box for element e
    return minimum(msh["x"][msh["v"][e,:],:],dims=1),maximum(msh["x"][msh["v"][e,:],:],dims=1)
end

function boBoxIntersect( bb1, bb2 )
# Checks if bb1 and bb2 intersect
    for i=1:length(bb1)
        if bb2[2][i] < bb1[1][i] || bb1[2][i] < bb2[1][i]
            return false
        end
    end
    return true
end

function pointInBoBox( p, bb )
# Checks if p is within bb
    return dot((p .>= bb[1]),(p .<= bb[2])) == 3
end

function bbox4surf( msh, surfID )
# Calculates bounding box for the given surface
    min,max=[1e7 1e7 1e7],[-1e7 -1e7 -1e7]
    xd = zeros(3,3)
    for f in eachrow(msh["faces"])
        if f[1] == surfID
            xd = msh["x"][f[2:4],:]
            for d=1:3
                if minimum( xd[:,d] ) < min[d] 
                    min[d] = minimum( xd[:,d] )
                end
                if maximum( xd[:,d] ) > max[d] 
                    max[d] = maximum( xd[:,d] )
                end
            end
        end
    end
    return ( min, max )
end

function elemsInBBox( msh, bbox )
    elems = Set{Int}()
    for e=1:size(msh["v"],1)
        if boBoxIntersect( elemBoBox3D4N( msh, e ), bbox )
        #if pointInBoBox( elemCG3D4N( msh,e ), bbox )
            push!( elems, e )
        end
    end
    return elems
end

function pointLineDist( line, point )
# Calculates distance of point (x,y,z) from line defined by point (x0,y0,z0) and vector (dx,dy,dz)
# Returns: distance, direction with respect to (x0,y0,z0) and the nearst point on the line 
    u = line[2]/norm(line[2]) # (dx,dy,dz) normalized
    v = point - line[1] # vector: line-basis to the point
    d = dot(u,v) # |u|*|v|*sin(u,v), |u| == 1, so we get |v|*sin(u,v) == the distance from line[1] to the point of the line nearest to the point
    n = cross(u,v) # |u|*|v|*cos(u,v),  |u| == 1, so we get |v|*cos(u,v) == the distance (may be -distance as well)
    return norm(n),sign(d),line[1]+d*u # the sign shows on which side of the (x0,y0,z0) the point is
end

function longestEdge3D4N( mesh, e )
# Longest egde(s) of element e
    x = mesh["x"]
    nds = vec(mesh["v"][e,:])
    v= x[nds,:]
    v0=[ v[2:4,:].-v[1,:]'; v[3:4,:].-v[2,:]'; v[4,:]'-v[3,:]' ]
    lgt=map(norm,eachrow(v0))
    return findall(x->x==maximum(lgt),lgt)
end

function edgeNo( nds, v1, v2 )
# Number of tetrahedra edge v1-v2 in element defined by nds
    i1 = findfirst( ==(v1), nds )
    i2 = findfirst( ==(v2), nds )
    if i1 == nothing || i2 == nothing
        return 0
    end
    fi = minimum( [i1, i2 ] )
    se = maximum( [i1, i2 ] )
    if fi == 1
        return se-1
    else 
        if fi == 2
            return se+1
        else
            return 6
        end
    end
end

function faceNo( nds, v1, v2, v3 )
# Number of tetrahedra facee v1-v2-v3 in element defined by nds
    
    i1 = findfirst( ==(v1), nds )
    i2 = findfirst( ==(v2), nds )
    i3 = findfirst( ==(v3), nds )
    if i1 == nothing || i2 == nothing || i3 == nothing
        return 0
    end
    fi = minimum( [i1, i2, i3 ] )
    th = maximum( [i1, i2, i3 ] )
    se = setdiff( [i1, i2, i3 ], [fi], [th] )
    map=Dict( (1,2,3)=>1, (1,2,4)=>2, (2,3,4)=>3, (1,3,4)=>4 )
    return map[(fi,se,th)]
end

function splitTetra( mesh, e, edge, edgeList, lev="" )  
# Bissects e-th element, cutting edge  (longest if edge==0)
    eev = [ 1 2; 1 3; 1 4; 2 3; 2 4; 3 4 ]  # maps edge -> end nodes
    nds = mesh["v"][e,:]
    
    le = longestEdge3D4N( mesh, e )
    #println( lev, "Split tetra ", e, " nodes ", nds, " on edge ", edge, " longest edges: ", le )
    if ! (edge in le)
        edge = le[1]
        enodes = [ nds[eev[edge,1]], nds[eev[edge,2]] ]
        #println( lev, " -> adding edge ", enodes )
        push!(edgeList,enodes)
        for el in setdiff( intersect( mesh["iv"][enodes[1]], mesh["iv"][enodes[2]] ), Set{Int}(e) )
            eledge = edgeNo(mesh["v"][el,:],enodes[1],enodes[2]) # needs to be checked: why some edges are 0?
            #println( "Elem ", el, "-> edge ", eledge )
            if eledge > 0
                splitTetra( mesh, el, eledge, edgeList, lev * "  " )
            else
                println( lev, "Skipping el $el, no edge ", enodes, " while splitting element ", e, " with nodes ", nds, " size(x)=", size(mesh["x"],1) )
								println( lev, "Element $el = ", mesh["v"][el,:] )
            end
        end
    end
end

function addNewVertex( mesh, newX, enodes )
# Adds new vertex to the mesh     

    mesh["x"] = vcat( mesh["x"],  newX  )
    
    mesh["iv"] = vcat( mesh["iv"], Set{Int}(intersect( mesh["iv"][enodes[1]], mesh["iv"][enodes[2]] )) )

    if "b" in keys(mesh) && mesh["b"] != nothing
        if mesh["b"][enodes[1]] == mesh["b"][enodes[2]]
            mesh["b"] = vcat( mesh["b"], [mesh["b"][enodes[1]]] )
        else
            mesh["b"] = vcat( mesh["b"], [0] )
        end
    end
    
    node = size(mesh["x"],1)

		#println( "Node ", node, " splits edge ", enodes )
    
    return node
end

function splitEdge3D4N( mesh, edge, refinedelements=nothing )
# Bissects all elements adjacent to edge (vector of vertices#)
    if ! ("iv" in keys(mesh)) || mesh["iv"] == nothing
        makeinop(mesh)
    end
    if length( intersect( mesh["iv"][edge[1]], mesh["iv"][edge[2]] ) ) == 0
        return nothing
    end
    newNode= addNewVertex( mesh, sum(mesh["x"][edge,:], dims=1) / 2, edge )
    for el in intersect( mesh["iv"][edge[1]], mesh["iv"][edge[2]] )
        eledge = edgeNo(mesh["v"][el,:],edge[1],edge[2])
        if eledge > 0 
    		#println( "   dividing $el=", mesh["v"][el,:] )
			divideTetra( mesh, el, eledge, newNode, refinedelements )
			#println( "      got $el=", mesh["v"][el,:] )
			nel = size(mesh["v"],1)
			#println( "      and $nel=", mesh["v"][nel,:] )
		end
    end
    return nothing
end

function removeFromInop( mesh, nodes, element )
# Removes element from map node -> elements
    for v in nodes
        delete!( mesh["iv"][v], element )
    end
end

function addToInop( mesh, nodes, element )
# Adds element to map node -> elements
    for v in nodes
        push!( mesh["iv"][v], element )
    end
end

function divideTetra( mesh, e, edge, newnode, refinedelements )
# Split e-th element, cutting edge by inserting newnode
    #println( "Split element ", e, " edge ", edge, " node ", newnode )
    e2v = ( [5 2 3 4; 1 5 3 4], [ 1 2 5 4; 2 3 5 4], [1 2 3 5 ; 5 2 3 4], [1 2 5 4; 1 5 3 4], [1 2 3 5; 1 5 3 4], [ 1 2 5 4; 1 2 3 5] ) # maps edge# -> vertices of bisection#
    efv = [ [1 2 4; 1 3 2 ], [1 4 3; 1 3 2], [1 2 4; 1 4 3], [ 1 3 2; 2 3 4], [1 2 4; 2 3 4], [2 3 4; 1 4 3] ] # maps edge->faces
    eev = [ 1 2; 1 3; 1 4; 2 3; 2 4; 3 4 ]  # maps edge -> end nodes
    
    oldvol = volume3D4N( mesh["x"][mesh["v"][e,:],:] )
    # Split element into two
    nds = [ mesh["v"][e,:]; newnode ]
		removeFromInop( mesh, mesh["v"][e,:], e )
    #println( "+spliting 5 nodes: ", nds)
    mesh["v"][e,:]= nds[e2v[edge][1,:]]'
		addToInop( mesh, mesh["v"][e,:], e )  
    if e <= size(mesh["err"],1)
        newvol =  volume3D4N( mesh["x"][mesh["v"][e,:],:] )
        mesh["err"][e] *= newvol/oldvol
    end
    #println( "redefine elem ",e,": ", nds[e2v[edge][1,:]]')
    mesh["v"]= vcat( mesh["v"], nds[e2v[edge][2,:]]')
	nel = size(mesh["v"],1)
	addToInop( mesh, mesh["v"][nel,:], nel )
    #println( "+elem: ", nds[e2v[edge][2,:]]', " -> ", size(mesh["v"],1))
    mesh["s"]= vcat( mesh["s"], mesh["s"][e] )

	push!( mesh["iv"][newnode], size(mesh["v"],1) )
    
    # Correct "faces" of the mesh -> if the splitted faces were marked - add the new faces to the marking
    if "faces" in keys(mesh) && mesh["faces"] != nothing
        enodes = [ nds[eev[edge,1]], nds[eev[edge,2]] ]
        for s=1:2
            face = nds[efv[edge][s,:]]
            idx = findall(x->length(intersect(x[2:4],face))==3,eachrow(mesh["faces"]))
            for i in idx
                oldface = mesh["faces"][i,:]
				#print( "   Face :", oldface )
                mesh["faces"][i,:] = replace( oldface, enodes[1]=>newnode )
				#print( " -> ", mesh["faces"][i,:] )
                mesh["faces"] = vcat(mesh["faces"], replace( oldface, enodes[2]=>newnode )')
				#println( " + ", mesh["faces"][end,:] )
            end
        end
    end
    
    if refinedelements != nothing
        push!(refinedelements,e)
        push!(refinedelements,size(mesh["v"],1))
    end
    return nothing
end

function refineMesh( mesh, localEE, threshold=0.5 )
# Refine mesh, bissecting elements for which localEE > threshold
    #println( "Refining mesh of ", size(mesh["x"],1), " vertices and ", size(mesh["v"],1), " elements" )
    if size(mesh["v"],1) != size(localEE,1)
        println( "Mesh and Error estimate are not compatible!")
        return
    end
    
    if !("iv" in keys(mesh)) || mesh["iv"] == nothing
        makeinop( mesh )
    end    
    
    oldne = size(localEE,1)
    maxEE= maximum(localEE)
    mesh["err"] = localEE
    for e=sort([1:oldne...],by=x->localEE[x],rev=true)
        if mesh["err"][e]/maxEE > threshold
            #println("-------------------------------------------------------------")
            #println( e, ": ", localEE[e]/maxEE )
            edgeList = Set{Any}()
            splitTetra( mesh, e, 0, edgeList )
            for edg in edgeList
                elemList = Set{Int}()
                splitEdge3D4N( mesh, edg, elemList )
                #makeinop(mesh)   # not necessary for all elements - should be improved 
                #println( "\tEdge: ", edg, " -> elements: ", elemList )
            end
        else
            break
        end
    end
		makeinop( mesh )
    delete!(mesh,"err")
		return nothing
end       

function makeinop( mesh )
# Generate mapping inverse to nop: node->elements
    x= mesh["x"]
    nop= mesh["v"]
    nv = size(x,1)
    ne,nve = size(nop)  
    pon = [ Set{Int}() for i=1:nv ]
    for e=1:ne
        for v=1:nve
            push!(pon[nop[e,v]],e)
        end
    end
    mesh["iv"] = pon
    #println("made inop for $nv vertices")
    return nothing
end

function triangleFaces( tf, nop, e )
# Faces (i.e. edges) of triangle e
    tf[1,1]= nop[e,1];
    tf[1,2]= nop[e,2];
    tf[2,1]= nop[e,2];
    tf[2,2]= nop[e,3];
    tf[3,1]= nop[e,3];
    tf[3,2]= nop[e,1];
    return nothing
end

function tetraFaces( tf, nop, e )
# Faces (i.e. triangles) of tetra e
    tf[1,1]= nop[e,1];
    tf[1,2]= nop[e,2];
    tf[1,3]= nop[e,3];
    tf[2,1]= nop[e,1];
    tf[2,2]= nop[e,2];
    tf[2,3]= nop[e,4];
    tf[3,1]= nop[e,2];
    tf[3,2]= nop[e,4];
    tf[3,3]= nop[e,3];
    tf[4,1]= nop[e,3];
    tf[4,2]= nop[e,4];
    tf[4,3]= nop[e,1];
    return nothing
end

function addBnd( mesh, xconditon::Function, bndNo )
# Adds boundary nbndNo for nodes matching xcondition (mesh["b"] = bndNo for these vertices)
    newBndNodes = findall(xconditon,eachrow(mesh["x"]))
    if ! ("b" in keys(mesh)) || mesh["b"] == nothing
        mesh["b"] = zeros(Int,size(mesh["x"],1))
    end
    mesh["b"][newBndNodes] .= bndNo
    return nothing 
end

function makeFaces( mesh )
# Generate list of mesh["faces"] out of mesh["b"]
    if "b" in keys(mesh) && mesh["b"] != nothing
        nop= mesh["v"]
        x= mesh["x"]
        nv = size(x,1)
        ne,nve = size(nop) 
        if nve == 3
            facesfun = triangleFaces
            faces = Set{NTuple{3,Int}}()
            tf = zeros(Int,3,2)
        else
            facesfun = tetraFaces
            faces = Set{NTuple{4,Int}}()
            tf = zeros(Int,4,3)
        end
        for e=1:ne
            facesfun( tf, nop, e )
            for f=1:size(tf,1)
                nodes= sort(tf[f,:])
                if mesh["b"][tf[f,1]] > 0 && all(x->x==mesh["b"][tf[f,1]],mesh["b"][nodes])
                    push!(faces, tuple( mesh["b"][tf[f,1]], nodes... ))
                end
            end
        end
        mesh["faces"] = (hcat(collect.(faces)...))'
    end
    return nothing
end

function faces2Bnds( mesh )
# Moves bnd markings from mesh["faces"] to mesh["b"]
    if "faces" in keys(mesh) && mesh["faces"] != nothing
        if ! ("b" in keys(mesh)) || mesh["b"] == nothing
            mesh["b"] = zeros(Int,size(mesh["x"],1))
        end
        for i=1:size(mesh["faces"],1)
            bnd = mesh["faces"][i,1]
            for v in mesh["faces"][i,2:end]
                if mesh["b"][v] == 0
                    mesh["b"][v] = bnd
                end
            end
        end
    end
    return nothing
end

function getXSurface( mesh, bndList )
# Generate list of faces out of bndList
# List needs to be further processed - it does not contain surface mark, just node-numbers
    if "b" in keys(mesh) && mesh["b"] != nothing
        nop= mesh["v"]
        x= mesh["x"]
        nv = size(x,1)
        ne,nve = size(nop) 
        if nve == 3
            facesfun = triangleFaces
            faces = Set{NTuple{3,Int}}()
            tf = zeros(Int,4,2)
        else
            facesfun = tetraFaces
            faces = Set{NTuple{3,Int}}()
            tf = zeros(Int,4,3)
        end
        for e=1:ne
            facesfun( tf, nop, e )
            for f in eachrow(tf)
                if ! ( nothing in indexin(mesh["b"][f],bndList) )
                    push!(faces, (f...,) )
                end
            end
        end
        return (hcat(collect.(faces)...))'
    else
        return nothing
    end
end

function markOuterSurface3D( mesh, bndNo )
# Adds outer surface of the whole domain to mesh["faces"]
# Assigns bndNo to this surface
    surf = Set{NTuple{3,Int}}()
    tf = zeros(Int,4,3)
    for e=1:size(mesh["v"],1)
        tetraFaces( tf, mesh["v"], e )
        for f in eachrow(tf)
            sf = Tuple(i for i in sort(f))
            if sf in surf
                delete!(surf,sf)
            else
                push!(surf,sf)
            end
        end
    end
    if "faces" in keys(mesh) && mesh["faces"] != nothing
        mesh["faces"] = vcat( mesh["faces"], [ bndNo*ones(Int, length(surf),1) [  t[k] for t in surf, k=1:3 ]] )
    else
        mesh["faces"] = [ bndNo*ones(Int, length(surf),1) [  t[k] for t in surf, k=1:3 ]]
    end
    if ! ("b" in keys(mesh)) || mesh["b"] == nothing
        mesh["b"] = zeros(Int,size(mesh["x"],1),1)
    end
    for f in eachrow(mesh["faces"])
        if f[1] == bndNo
            mesh["b"][f[2:end]] .= bndNo
        end
    end
    return nothing
end

function markSubvolBnd3D( mesh, sub, bndNo )
# Adds outer boundary of subdomain sub to the mesh["faces"].
# Assignd bndNo to this surface.
    surf = Set{NTuple{3,Int}}()
    tf = zeros(Int,4,3)
    for e=1:size(mesh["v"],1)
        if mesh["s"][e] == sub
            tetraFaces( tf, mesh["v"], e )
            for f in eachrow(tf)
                sf = Tuple(i for i in sort(f))
                if sf in surf
                    delete!(surf,sf)
                else
                    push!(surf,sf)
                end
            end
        end
    end
    if "faces" in keys(mesh) && mesh["faces"] != nothing
        mesh["faces"] = vcat( mesh["faces"], [ bndNo*ones(Int, length(surf),1) [  t[k] for t in surf, k=1:3 ]] )
    else
        mesh["faces"] = [ bndNo*ones(Int, length(surf),1) [  t[k] for t in surf, k=1:3 ]]
    end
    if ! ("b" in keys(mesh)) || mesh["b"] == nothing
        mesh["b"] = zeros(Int,size(mesh["x"],1),1)
    end
    for f in eachrow(mesh["faces"])
        if f[1] == bndNo
            mesh["b"][f[2:end]] .= bndNo
        end
    end
    return nothing
end

function makeSurfaceView( mesh, V, bndList )
# Creates surface mesh for bndList and associated vie of V (field) on this surface
    nop = getXSurface(mesh,bndList)
    if length(nop) == 0
        return nothing,nothing
    end
    nodes=sort(collect(Set(vec(sort(reshape(nop,:,1),dims=1)))))
    bnd = mesh["b"][nodes]
    x = mesh["x"][nodes,:]
    surfV = V[nodes]
    for i=1:size(nop,1)
            for j=1:size(nop,2)
                nop[i,j]= findfirst(x->nop[i,j]==nodes[x],1:size(nodes,1))
            end
    end
    return Dict("v"=>nop,"x"=>x,"b"=>bnd),surfV
end

function addXSurface( mesh, bndList, number )
# Adds surface spanning bndList to mesh["faces"].
# Assigns number to this surface
    facelist = getXSurface(mesh,bndList)
    if "faces" in keys(mesh) && mesh["faces"] != nothing
        mesh["faces"] = vcat( mesh["faces"], [ number*ones(Int,size(facelist,1),1) facelist ] )
    else
        mesh["faces"] = [ number*ones(Int,size(facelist,1),1) facelist ]
    end
    return nothing
end
    
function face2face( f1, f2 )
# Check if f1 and f2 are same
    return length(f1) == length(f2) && length(f1) == length(intersect(f1,f2))
end

function extractSubMesh( mesh; vertexList=nothing, elementList=nothing )
# Extract sub-mesh with vertices and elements listed
    if vertexList == nothing && elementList == nothing
        return nothing
    end
    eList = Set{Int}()
    if elementList != nothing
        for e in elementList
            push!(eList,e)
        end
    end
    if vertexList != nothing
        makeinop( mesh )
        for v in vertexList
            for e in mesh["iv"][v]
                push!(eList,e)
            end
        end
    end
    elements = sort(collect(eList))
    nop= mesh["v"][elements,:]
    mat = mesh["s"][elements]
    nodes=sort(collect(Set(vec(sort(reshape(nop,:,1),dims=1)))))
    bnd = mesh["b"][nodes]
    x = mesh["x"][nodes,:]
    for i=1:size(nop,1)
            for j=1:size(nop,2)
                nop[i,j]= findfirst(x->nop[i,j]==nodes[x],1:size(nodes,1))
            end
    end
    return Dict("v"=>nop,"x"=>x,"b"=>bnd,"s"=>mat)
end

function makeSurfOutOfSubdomain( msh, surfID, sbdID )
# Reorder points of triangular faces of surf, so their normals will point out of the subdomain.
    if "faces" in keys(msh) && msh["faces"] != nothing
        if !("iv" in keys(msh)) || msh["iv"] == nothing
            makeinop( msh )
        end 
        nds=[0,0,0]
        ndsr=[0,0,0,0]
        nflips=0
        for f in eachrow(msh["faces"])
            if f[1] == surfID
                # find element(s) adjacent to the face
                nds = f[2:4]
                elems = [e for e in intersect(msh["iv"][nds]...)]
                subs = msh["s"][elems]
                # select the element from sbdID subdomain 
                e = findall(x->x==sbdID,subs)
                if length(e) > 1
                    println( "Warning: $f is not the outer face of subdomain $sbdID" )
                else
                    # find the node not belonging to the face
                    elem = elems[e[1]]
                    fourth = [v for v in setdiff(msh["v"][elem,:], nds)]
                    if length(fourth) == 1
                        # check if the normal of the face points towards the fourth node of tetrahedra
                        ndsr = [nds...,fourth[1]]
                        test = det( [ msh["x"][ndsr,:] ones(4,1) ] )
                        # swap orientatnion of the face if it does
                        if test < 0
                            f[4],f[3] = f[3],f[4]
                            nflips += 1
                        end
                    end
                end
            end
        end
        return nflips
    end
    return nothing
end

function skinElectrode( msh, line, rad, side, skinID, eleID; adjustShape=true )
# Creates new surface (eleID) taking those patches of skinID, points of which are not 
# further than rad from the line (like cutting skin with cylinder). If adjustShape is set,
# the points of skinID are moved to form the circle.
    faces = Set{NTuple{4,Int}}()
    side = sign(side)
    for i=1:size(msh["faces"],1)
        if msh["faces"][i,1] == skinID
            nin = 0
            sgn = 0
            for j=2:size(msh["faces"][i,:],1)
                r,s,nearest = pointLineDist( line, msh["x"][msh["faces"][i,j],:] )
                if r <= rad*1.1
                    nin += 1
                    sgn += s
                end
            end
            if nin  == 3
                if sign(sgn) == side 
                    push!(faces, tuple( eleID, msh["faces"][i,2:4]... ))
                end
            end
        end
    end
    if adjustShape
        out = Set{NTuple{2,Int}}()
        for f in faces
            nds = [f[2:4]...,f[2]]
            for e=1:3
                edg = Tuple( i for i in sort([nds[e],nds[e+1]]) )
                if edg in out
                    delete!(out,edg)
                else
                    push!(out,edg)
                end
            end        
        end
        verts = Set{Int}()
        for e in out
            push!(verts,e[1])
            push!(verts,e[2])
        end
        for v in verts
            curpos = msh["x"][v,:]
            r,s,nearest = pointLineDist( line, curpos )
            newpos = nearest + rad/r*(curpos-nearest)
            msh["x"][v,:] = newpos
        end
    end
    return (hcat(collect.(faces)...))'
end

function box2D( x_min, y_min, x_max, y_max, nx, ny )
# Makes 2D regular grid in a box
    xt = range( x_min, stop=x_max, length=nx )
    yt = range( y_min, stop=y_max, length=ny )

    x = [ repeat(xt,inner=ny,outer=1) repeat(yt,outer=nx) ]

    nop = repeat( [ 1:(nx-1)*(ny-1) convert.(Int32,floor.(((1:(nx-1)*(ny-1)).-1)./(ny-1))).+1  ((1:(nx-1)*(ny-1)).-1).%(ny-1).+1 ], inner=(2,1), outer=(1,1))

    for i=1:2:size(nop,1)
        r = nop[i,2]
        c = nop[i,3]
        nop[i,1] = (r-1)*ny+c
        nop[i,2] = (r-1)*ny+c+ny
        nop[i,3] = (r-1)*ny+c+1
        nop[i+1,1] = nop[i,2]
        nop[i+1,2] = nop[i,2]+1
        nop[i+1,3] = nop[i,3]
    end
    
    return Dict( "x"=>x, "v"=>nop )
end

function box3D( x_min, y_min, x_max, y_max, z_min, z_max, nx, ny, nz )  # NOT READY YET!!!
# Make 3D regular grid in a box
    xt = range( x_min, stop=x_max, length=nx )
    yt = range( y_min, stop=y_max, length=ny )
    zt = range( z_min, stop=z_max, length=nz )

    x = [ repeat(xt,inner=ny,outer=1) repeat(yt,outer=nx) repeat(zt,inner=nz,outer=nx*ny) ]

    nop = repeat( [ 1:(nx-1)*(ny-1) convert.(Int32,floor.(((1:(nx-1)*(ny-1)).-1)./(ny-1))).+1  ((1:(nx-1)*(ny-1)).-1).%(ny-1).+1 ], inner=(2,1), outer=(1,1))

    for i=1:2:size(nop,1)
        r = nop[i,2]
        c = nop[i,3]
        nop[i,1] = (r-1)*ny+c
        nop[i,2] = (r-1)*ny+c+ny
        nop[i,3] = (r-1)*ny+c+1
        nop[i+1,1] = nop[i,2]
        nop[i+1,2] = nop[i,2]+1
        nop[i+1,3] = nop[i,3]
    end
    
    return Dict( "x"=>x, "v"=>nop )
end

function transf( p, m )
# 2D transformation of point x given transformation matrix M
	return ([p  ones(size(p,1),1)] * m)[:,1:2]
end

function rad2cart( x )
# 2D transform (r,phi) -> (x,y)   
	r2c(ra) = [ ra[1]*cos(ra[2]) ra[1]*sin(ra[2]) ]
	return vcat(r2c.( [ x[i,:] for i=1:size(x,1)] )...) 
end


############## LEGACY !!!! #################################################################################################################
function addVertex( mesh, newVertex, edge, enodes )
# Adds new vertex to the mesh     
    efv = [ [1 2 4; 1 3 2 ], [1 4 3; 1 3 2], [1 2 4; 1 4 3], [ 1 3 2; 2 3 4], [1 2 4; 2 3 4], [2 3 4; 1 4 3] ]
    mesh["x"] = vcat( mesh["x"],  newVertex  )
    node = size(mesh["x"],1)
    mesh["iv"] = vcat( mesh["iv"], Set{Int}(intersect( mesh["iv"][enodes[1]], mesh["iv"][enodes[2]] )) )
    mesh["b"] = vcat( mesh["b"], [0] )
    #println( "+node: ", sum(mesh["x"][enodes,:], dims=1) / 2, " -> ", node  )
    if "faces" in keys(mesh) && mesh["faces"] != nothing
        for s=1:2
            face = efv[edge][s,:]
            idx = findall(x->length(intersect(x[2:4],face))==3,eachrow(mesh["faces"]))
            for i in idx
                oldface = mesh["faces"][i,:]
                mesh["faces"][i,:] = replace( oldface, enodes[1]=>node )
                mesh["faces"] = vcat(mesh["faces"], replace( oldface, enodes[2]=>node )')
            end
        end
    end
    if "b" in keys(mesh) && mesh["b"] != nothing
        if mesh["b"][enodes[1]] == mesh["b"][enodes[2]]
            mesh["b"] = vcat( mesh["b"], [mesh["b"][enodes[1]]] )
        else
            mesh["b"] = vcat( mesh["b"], [0] )
        end
    end
    return node
end

function bissectTetra( mesh, e, edge, newnode, refinedelements=nothing  )   #### !!!! Does not work correct !!!!
 # Bissect e-th element, cutting edge  (longest if edge==0)
    #println( "element ", e, " edge ", edge, " node ", newnode )
    e2v = ( [5 2 3 4; 1 5 3 4], [ 1 2 5 4; 2 3 5 4], [1 2 3 5 ; 5 2 3 4], [1 2 5 4; 1 5 3 4], [1 2 3 5; 1 5 3 4], [ 1 2 5 4; 1 2 3 5] ) # maps edge# -> vertices of bisection#
    eev = [ 1 2; 1 3; 1 4; 2 3; 2 4; 3 4 ]  # maps edge -> end nodes
            
    le = longestEdge( mesh, e )
    if edge in le
        # split e 
        #println( "Edge ", edge, " is one of the longest in element ",e, " - just splitting it")
        nds = [ mesh["v"][e,:]; newnode ]
        #println( "+spliting 5 nodes: ", nds)
        mesh["v"][e,:]= nds[e2v[edge][1,:]]'
        #println( "redefine elem ",e,": ", nds[e2v[edge][1,:]]')
        mesh["v"]= vcat( mesh["v"], nds[e2v[edge][2,:]]')
        #println( "+elem: ", nds[e2v[edge][2,:]]', " -> ", size(mesh["v"],1))
        mesh["s"]= vcat( mesh["s"], mesh["s"][e] )
        if refinedelements != nothing
            push!(refinedelements,e)
            push!(refinedelements,size(mesh["v"],1))
        end
    else
        nds = mesh["v"][e,:]
        savedNodes = nothing
        if edge > 0
            savedNodes = [ nds[eev[edge,1]], nds[eev[edge,2]] ]
        end
        edge = le[1]
        #println(nds)
        #println( " bissecting edge ",eev[edge,:] )
        enodes = [ nds[eev[edge,1]], nds[eev[edge,2]] ]
        #println( " with nodes: ", enodes )
        #println( " (coordinates: ", mesh["x"][enodes,:], " )" )
        node= addVertex( mesh, sum(mesh["x"][enodes,:], dims=1) / 2, edge, enodes )
        for el in setdiff( intersect( mesh["iv"][enodes[1]], mesh["iv"][enodes[2]] ), Set{Int}(e) )
            eledge = edgeNo(mesh["v"][el,:],enodes[1],enodes[2]) # needs to be checked: why some edges are 0?
            #println( "Elem ", el, "-> edge ", eledge )
            if eledge > 0
                bissectTetra( mesh, el, eledge, node, refinedelements )
            else
                #println( "Skipping el $el, no edge ", enodes )
            end
        end
        nds = [ mesh["v"][e,:]; node ]
        #println( "After recursion: spliting 5 nodes: ", nds)
        mesh["v"][e,:]= nds[e2v[edge][1,:]]'
        #println( "redefine elem ",e,": ", nds[e2v[edge][1,:]]')
        mesh["v"]= vcat( mesh["v"], nds[e2v[edge][2,:]]')
        #println( "+elem: ", nds[e2v[edge][2,:]]', " -> ", size(mesh["v"],1))
        mesh["s"]= vcat( mesh["s"], mesh["s"][e] )
        e2 = size(mesh["v"],1)
        
        if savedNodes != nothing
            ee = edgeNo( mesh["v"][e,:], savedNodes[1], savedNodes[2] )
            if ee > 0
                bissectTetra( mesh, e, ee, newnode, refinedelements )
            end
            ee = edgeNo( mesh["v"][e2,:], savedNodes[1], savedNodes[2] )
            if ee > 0
                bissectTetra( mesh, e2, ee, newnode, refinedelements )
            end
        end
        
        if refinedelements != nothing
            push!(refinedelements,e)
            push!(refinedelements,size(mesh["v"],1))
        end
    end
end
