using DelimitedFiles

function saveArr( file, arr )
    open(file,"w") do f
        writedlm(f,arr)
    end
end

function loadArr( file )
    return readdlm(file)
end

function loadArr( file, atype )
    return readdlm(file, atype)
end

function loadNetgenNeutralFormat( file )
    open(file) do f
        nv = parse(Int,readline(f))
        x=zeros(nv,3)
        for i=1:nv
            x[i,:]= parse.(Float64,split(readline(f)))
        end
        ne = parse(Int,readline(f))
        nop = zeros(Int,ne,4)
        mat = zeros(Int,ne,1)
        for i=1:ne
            l = parse.(Int,split(readline(f)))
            mat[i]= l[1]
            nop[i,:]= l[2:5]
        end
        nf = parse(Int,readline(f))
        faces = zeros(Int,nf,4)
        bnd = zeros(Int,nv,1)
        for i=1:nf
            l = parse.(Int,split(readline(f)))
            faces[i,:] = l
            for v in l[2:4]
                bnd[v]= l[1]
            end
        end
        return Dict("x"=>x, "v"=>nop, "s"=>vec(mat), "b"=>vec(bnd), "faces"=>faces )
    end
end

function loadMeditINRIA( file )
    open(file) do f
        if ! occursin( "MeshVersionFormatted 2", readline(f) )
            println( "Wrong headder (Version)" )
            return nothing
        end
        if ! occursin( "Dimension", readline(f) )
            println( "Wrong headder (dimension)" )
            return nothing
        end
        dim = parse(Int,readline(f))
        if dim < 2 || dim > 3
            println( "Wrong dimension" )
            return nothing
        end
        
        if ! occursin( "Vertices", readline(f) )
            println( "No vertices" )
            return nothing
        end
        nv = parse(Int,readline(f))
        x = zeros(nv,dim)
        bnd = zeros(Int,nv,1)
        for i=1:nv
            x[i,:] = parse.(Float64,split(readline(f)))[1:dim]
        end
        
        headline = readline(f)
        if occursin( "Edges", headline )
            nedge = parse(Int,readline(f))
            edges = zeros(Int,nedge,3)
            for i=1:nedge
                l = parse.(Int,split(readline(f)))
                edges[i,:] = [ l[3] l[1:2]' ]
            end
            headline = readline(f)
        else
            edges = nothing     
        end
        
        if occursin( "Triangles", headline )
            nt = parse(Int,readline(f))
            triangles = zeros(Int,nt,4)
            for i=1:nt
                l = parse.(Int,split(readline(f)))
                triangles[i,:] = [l[4] l[1:3]']
            end
            headline = readline(f)
        else
            triangles = nothing
        end
        
        if occursin( "Tetrahedra", headline )
            nt = parse(Int,readline(f))
            tetra = zeros(Int,nt,4)
            mat = zeros(Int,nt,1)
            for i=1:nt
                l = parse.(Int,split(readline(f)))
                tetra[i,:] = l[1:4]
                mat[i]= l[5]
            end
        else
            tetra = nothing
        end
        
        if dim == 2
            for i=1:size(edges,1)
                for v in edges[i,2:3]
                    bnd[v] = edges[i,1]
                end
            end
            mat = triangles[:,1]
            return Dict{String,Any}("x"=>x, "v"=>triangles[:,2:4], "s"=>vec(mat), "b"=>vec(bnd), "faces"=>edges )
        else
            for i=1:size(triangles,1)
                for v in triangles[i,2:4]
                    bnd[v] = triangles[i,1]
                end
            end
            #return Dict{String,Array}("x"=>x, "v"=>tetra, "s"=>vec(mat), "b"=>vec(bnd), "faces"=>triangles, "edges"=>edges )
            return Dict{String,Any}("x"=>x, "v"=>tetra, "s"=>vec(mat), "b"=>vec(bnd), "faces"=>triangles, "edges"=>edges )
        end
    end
end

function saveMeditINRIA( file, mesh )
    open(file,"w") do f
        x = mesh["x"]
        nop = mesh["v"]
        if "b" in keys(mesh) && mesh["b"] != nothing
            bnd = mesh["b"]
        else
            bnd = [ 0 for i=1:size(x,1) ]
        end
        if "s" in keys(mesh) && mesh["s"] != nothing
            sbd = mesh["s"]
        else
            sbd = [ 0 for i=1:size(nop,1) ]
        end
        nv,dim = size(x)
        write(f, " MeshVersionFormatted 2\n" )
        write(f, " Dimension\n" )
        write(f, " $dim\n" )
        
        write(f, " Vertices\n" )
        write(f, " $nv\n" )
        for i=1:nv
            v=x[i,:]
            b = bnd[i]
            for v in x[i,:]
                write( f, " $v")
            end
            write( f, " $b\n")
        end
        
        write(f, " Edges\n" )
        if "edges" in keys(mesh) && mesh["edges"] != nothing
            ne = size(mesh["edges"],1)
            write( f, " $ne\n")
            for i=1:ne
                b,v1,v2 = mesh["edges"][i,:]
                write( f, " $v1 $v2 $b\n")
            end
        else
            write( f, " 0\n" )
        end
        
        write(f, " Triangles\n" )
        if dim == 2
            ne = size(mesh["v"],1)
            write( f, " $ne\n")
            for i=1:ne
                v1,v2,v3 = mesh["v"][i,:]
                m = sbd[i]
                write( f, " $v1 $v2 $v3 $m\n")
            end
        else
            if "faces" in keys(mesh) && mesh["faces"] != nothing
                ne = size(mesh["faces"],1)
                write( f, " $ne\n")
                for i=1:ne
                    b,v1,v2,v3 = mesh["faces"][i,:]
                    write( f, " $v1 $v2 $v3 $b\n")
                end
            else
                write( f, " 0\n")
            end
        end
        
        write(f, " Tetrahedra\n" )
        if dim == 2
            write( f, " 0\n")
        else
            ne = size(mesh["v"],1)
            write( f, " $ne\n")
            for i=1:ne
                v1,v2,v3,v4 = mesh["v"][i,:]
                m = mesh["s"][i]
                write( f, " $v1 $v2 $v3 $v4 $m\n")
            end
        end
        
        write( f, "End\n" )
    end
    return nothing
end

using Pickle
function loadNG( file )
# Load mesh from pickle generated by netgen
    d=Pickle.load(file)
    nop = Int.(reduce(hcat,d["v"])')
    x= Float64.(reduce(hcat,d["x"])')
    mat= Int.(d["s"])
    if "b" in keys(d)
        bnd=zeros(Int,size(x,1),1)
        for v in d["b"]
            bnd[v] = 1;
        end
    else
        bnd = nothing
    end
    return Dict{String,Array}("x"=>x, "v"=>nop, "s"=>vec(mat), "b"=>vec(bnd) )
end

function loadOldNG( file )
# Load mesh from pickle generated by netgen
    d=Pickle.load(file)
    nop = Int.(reduce(hcat,d["nop"])')
    x= Float64.(reduce(hcat,d["x"])')
    mat= Int.(d["mat"])
    if "bnd" in keys(d)
        bnd=zeros(Int,size(x,1),1)
        for v in d["bnd"]
            bnd[v] = 1;
        end
    else
        bnd = nothing
    end
    return Dict{String,Array}("x"=>x, "v"=>nop, "s"=>vec(mat), "b"=>vec(bnd) )
end