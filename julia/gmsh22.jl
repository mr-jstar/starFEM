# ============================================================
# Strumieniowy importer Gmsh MSH 2.2 ASCII dla StarMES
#
# Wynik:
#   msh["x"]                  :: Matrix{Float64}   (nn × 3)
#   msh["v"]                  :: Matrix{Int64}     (ne × 4)
#   msh["faces"]              :: Matrix{Int64}     (nf × 4)
#                                 kolumny: geometrical_tag, n1, n2, n3
#   msh["s"]                  :: Vector{Int64}     (ne)
#                                 geometryczny tag objętości elementu
#   msh["b"]                  :: Vector{Int64}     (nn)
#                                 jeden geometryczny tag brzegu na węzeł
#   msh["physical names"]     :: Dict{Tuple{Int64,Int64},String}
#                                 (dim, physical_tag) => nazwa
#   msh["physical tags"]      :: Dict{String,Tuple{Int64,Int64}}
#                                 nazwa => (dim, physical_tag)
#   msh["physical entities"]  :: Dict{String,Set{Int64}}
#                                 nazwa => zbiór tagów geometrycznych
#
# Obsługiwane elementy:
#   Gmsh type 2 : trójkąt 3-węzłowy
#   Gmsh type 4 : czworościan 4-węzłowy
#
# Plik musi być zapisany jako MSH 2.2 ASCII, np. w Python API:
#   gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
#   gmsh.option.setNumber("Mesh.Binary", 0)
#   gmsh.write("model.msh")
# ============================================================


# ------------------------------------------------------------
# Funkcje pomocnicze wejścia
# ------------------------------------------------------------

function _gmsh_expect_line(io::IO, expected::AbstractString)
    eof(io) && error("Nieoczekiwany koniec pliku; oczekiwano: $expected")
    line = strip(readline(io))
    line == expected || error("Oczekiwano '$expected', otrzymano '$line'")
    return nothing
end


function _gmsh_skip_section!(io::IO, start_marker::AbstractString)
    startswith(start_marker, "\$") ||
        error("Niepoprawny znacznik sekcji: $start_marker")

    section_name = start_marker[2:end]
    end_marker = "\$End" * section_name

    while !eof(io)
        strip(readline(io)) == end_marker && return nothing
    end

    error("Brak znacznika końca sekcji: $end_marker")
end


@inline function _gmsh_node_index(
    tag_to_index::Dict{Int64,Int64},
    node_tag::Int64
)
    idx = get(tag_to_index, node_tag, Int64(0))
    idx != 0 || error("Element odwołuje się do nieznanego węzła Gmsh: $node_tag")
    return idx
end


# ------------------------------------------------------------
# $PhysicalNames
# ------------------------------------------------------------

function _gmsh_read_physical_names22!(
    io::IO,
    names::Dict{Tuple{Int64,Int64},String}
)
    n = parse(Int, strip(readline(io)))
    sizehint!(names, n)

    for _ in 1:n
        line = strip(readline(io))

        m = match(
            r"^\s*(-?\d+)\s+(-?\d+)\s+\"(.*)\"\s*$",
            line
        )

        m === nothing && error("Niepoprawny wpis w \$PhysicalNames: $line")

        dim  = parse(Int64, m.captures[1])
        tag  = parse(Int64, m.captures[2])
        name = m.captures[3]

        names[(dim, tag)] = name
    end

    _gmsh_expect_line(io, "\$EndPhysicalNames")
    return nothing
end


# ------------------------------------------------------------
# $Nodes
# ------------------------------------------------------------

function _gmsh_read_nodes22(io::IO)
    nn = parse(Int, strip(readline(io)))

    node_tags = Vector{Int64}(undef, nn)
    x = Matrix{Float64}(undef, nn, 3)

    for i in 1:nn
        fields = split(strip(readline(io)))
        length(fields) >= 4 || error("Niepoprawny rekord węzła nr $i")

        node_tags[i] = parse(Int64, fields[1])
        x[i,1] = parse(Float64, fields[2])
        x[i,2] = parse(Float64, fields[3])
        x[i,3] = parse(Float64, fields[4])
    end

    _gmsh_expect_line(io, "\$EndNodes")
    return node_tags, x
end


# ------------------------------------------------------------
# $Elements
# ------------------------------------------------------------

function _gmsh_read_elements22(io::IO)
    number_of_elements = parse(Int, strip(readline(io)))

    tet_nodes       = NTuple{4,Int64}[]
    tet_physical    = Int64[]
    tet_geometrical = Int64[]

    face_nodes       = NTuple{3,Int64}[]
    face_physical    = Int64[]
    face_geometrical = Int64[]

    sizehint!(tet_nodes, number_of_elements)
    sizehint!(tet_physical, number_of_elements)
    sizehint!(tet_geometrical, number_of_elements)

    for k in 1:number_of_elements
        fields = split(strip(readline(io)))
        length(fields) >= 3 || error("Niepoprawny rekord elementu nr $k")

        element_type = parse(Int, fields[2])
        number_of_tags = parse(Int, fields[3])

        physical_tag = number_of_tags >= 1 ? parse(Int64, fields[4]) : Int64(0)
        geometrical_tag = number_of_tags >= 2 ? parse(Int64, fields[5]) : Int64(0)

        first_node = 4 + number_of_tags

        if element_type == 2
            # trójkąt 3-węzłowy
            length(fields) >= first_node + 2 ||
                error("Niepoprawny trójkąt liniowy w rekordzie $k")

            push!(face_nodes, (
                parse(Int64, fields[first_node]),
                parse(Int64, fields[first_node + 1]),
                parse(Int64, fields[first_node + 2])
            ))
            push!(face_physical, physical_tag)
            push!(face_geometrical, geometrical_tag)

        elseif element_type == 4
            # czworościan 4-węzłowy
            length(fields) >= first_node + 3 ||
                error("Niepoprawny czworościan liniowy w rekordzie $k")

            push!(tet_nodes, (
                parse(Int64, fields[first_node]),
                parse(Int64, fields[first_node + 1]),
                parse(Int64, fields[first_node + 2]),
                parse(Int64, fields[first_node + 3])
            ))
            push!(tet_physical, physical_tag)
            push!(tet_geometrical, geometrical_tag)
        end
    end

    _gmsh_expect_line(io, "\$EndElements")

    return (
        tet_nodes,
        tet_physical,
        tet_geometrical,
        face_nodes,
        face_physical,
        face_geometrical
    )
end


# ------------------------------------------------------------
# physical name -> geometryczne encje
# ------------------------------------------------------------

function _gmsh_make_physical_entities(
    physical_names::Dict{Tuple{Int64,Int64},String},
    tet_physical::Vector{Int64},
    tet_geometrical::Vector{Int64},
    face_physical::Vector{Int64},
    face_geometrical::Vector{Int64}
)
    result = Dict{String,Set{Int64}}()

    for (_, name) in physical_names
        result[name] = Set{Int64}()
    end

    for i in eachindex(face_physical)
        name = get(physical_names, (Int64(2), face_physical[i]), nothing)
        name === nothing && continue
        push!(result[name], face_geometrical[i])
    end

    for i in eachindex(tet_physical)
        name = get(physical_names, (Int64(3), tet_physical[i]), nothing)
        name === nothing && continue
        push!(result[name], tet_geometrical[i])
    end

    return result
end


# ------------------------------------------------------------
# Składanie wyniku
# ------------------------------------------------------------

function _gmsh_assemble22(
    node_tags::Vector{Int64},
    x::Matrix{Float64},
    tet_nodes::Vector{NTuple{4,Int64}},
    tet_physical::Vector{Int64},
    tet_geometrical::Vector{Int64},
    face_nodes::Vector{NTuple{3,Int64}},
    face_physical::Vector{Int64},
    face_geometrical::Vector{Int64},
    physical_names::Dict{Tuple{Int64,Int64},String}
)
    nn = length(node_tags)
    ne = length(tet_nodes)
    nf = length(face_nodes)

    length(tet_physical) == ne || error("Niezgodna liczba tagów fizycznych czworościanów")
    length(tet_geometrical) == ne || error("Niezgodna liczba tagów geometrycznych czworościanów")
    length(face_physical) == nf || error("Niezgodna liczba tagów fizycznych ścian")
    length(face_geometrical) == nf || error("Niezgodna liczba tagów geometrycznych ścian")

    tag_to_index = Dict{Int64,Int64}()
    sizehint!(tag_to_index, nn)

    for i in eachindex(node_tags)
        tag = node_tags[i]
        haskey(tag_to_index, tag) && error("Powtórzony tag węzła Gmsh: $tag")
        tag_to_index[tag] = Int64(i)
    end

    v = Matrix{Int64}(undef, ne, 4)
    s = Vector{Int64}(undef, ne)

    for e in 1:ne
        nodes = tet_nodes[e]
        v[e,1] = _gmsh_node_index(tag_to_index, nodes[1])
        v[e,2] = _gmsh_node_index(tag_to_index, nodes[2])
        v[e,3] = _gmsh_node_index(tag_to_index, nodes[3])
        v[e,4] = _gmsh_node_index(tag_to_index, nodes[4])
        s[e] = tet_geometrical[e]
    end

    faces = Matrix{Int64}(undef, nf, 4)

    for f in 1:nf
        nodes = face_nodes[f]
        faces[f,1] = face_geometrical[f]
        faces[f,2] = _gmsh_node_index(tag_to_index, nodes[1])
        faces[f,3] = _gmsh_node_index(tag_to_index, nodes[2])
        faces[f,4] = _gmsh_node_index(tag_to_index, nodes[3])
    end

    # b zachowuje tylko jeden tag geometryczny na węzeł.
    # Dla węzłów należących do kilku powierzchni wybierany jest
    # najmniejszy niezerowy tag. Pełna informacja pozostaje w "faces".
    b = zeros(Int64, nn)

    for f in 1:nf
        geometrical_tag = faces[f,1]
        geometrical_tag == 0 && continue

        @inbounds for j in 2:4
            node = faces[f,j]
            if b[node] == 0 || geometrical_tag < b[node]
                b[node] = geometrical_tag
            end
        end
    end

    physical_tags = Dict{String,Tuple{Int64,Int64}}()
    sizehint!(physical_tags, length(physical_names))

    for ((dim, physical_tag), name) in physical_names
        haskey(physical_tags, name) &&
            error("Powtórzona nazwa grupy fizycznej: '$name'")
        physical_tags[name] = (dim, physical_tag)
    end

    physical_entities = _gmsh_make_physical_entities(
        physical_names,
        tet_physical,
        tet_geometrical,
        face_physical,
        face_geometrical
    )

    return Dict{String,Any}(
        "x"                  => x,
        "v"                  => v,
        "faces"              => faces,
        "s"                  => s,
        "b"                  => b,
        "physical names"     => physical_names,
        "physical tags"      => physical_tags,
        "physical entities"  => physical_entities
    )
end


# ------------------------------------------------------------
# Główna funkcja importera
# ------------------------------------------------------------

function loadGmsh22(filename::AbstractString)
    physical_names = Dict{Tuple{Int64,Int64},String}()

    node_tags = Int64[]
    x = Matrix{Float64}(undef, 0, 3)

    tet_nodes       = NTuple{4,Int64}[]
    tet_physical    = Int64[]
    tet_geometrical = Int64[]

    face_nodes       = NTuple{3,Int64}[]
    face_physical    = Int64[]
    face_geometrical = Int64[]

    version_found = false
    nodes_found = false
    elements_found = false

    open(filename, "r") do io
        while !eof(io)
            marker = strip(readline(io))
            isempty(marker) && continue

            if marker == "\$MeshFormat"
                fields = split(strip(readline(io)))
                length(fields) >= 3 || error("Niepoprawna sekcja \$MeshFormat")

                version = fields[1]
                filetype = parse(Int, fields[2])

                startswith(version, "2.") ||
                    error("Oczekiwano MSH 2.x, otrzymano wersję $version")

                filetype == 0 ||
                    error("Obsługiwany jest tylko format MSH 2.x ASCII")

                _gmsh_expect_line(io, "\$EndMeshFormat")
                version_found = true

            elseif marker == "\$PhysicalNames"
                _gmsh_read_physical_names22!(io, physical_names)

            elseif marker == "\$Nodes"
                node_tags, x = _gmsh_read_nodes22(io)
                nodes_found = true

            elseif marker == "\$Elements"
                (
                    tet_nodes,
                    tet_physical,
                    tet_geometrical,
                    face_nodes,
                    face_physical,
                    face_geometrical
                ) = _gmsh_read_elements22(io)
                elements_found = true

            elseif startswith(marker, "\$") && !startswith(marker, "\$End")
                _gmsh_skip_section!(io, marker)

            else
                error("Nieoczekiwana zawartość poza sekcją Gmsh: $marker")
            end
        end
    end

    version_found || error("Brak sekcji \$MeshFormat")
    nodes_found || error("Brak sekcji \$Nodes")
    elements_found || error("Brak sekcji \$Elements")

    return _gmsh_assemble22(
        node_tags,
        x,
        tet_nodes,
        tet_physical,
        tet_geometrical,
        face_nodes,
        face_physical,
        face_geometrical,
        physical_names
    )
end


# Alias zgodny z planowaną nazwą w bibliotece.
loadGmsh(filename::AbstractString) = loadGmsh22(filename)


# ------------------------------------------------------------
# Funkcje pomocnicze do pracy z grupami fizycznymi
# ------------------------------------------------------------

"""
    physicalTag(msh, name)

Zwraca `(dimension, physical_tag)` dla grupy o podanej nazwie.
"""
function physicalTag(msh::AbstractDict, name::AbstractString)
    tags = msh["physical tags"]
    haskey(tags, name) || error("Brak grupy fizycznej o nazwie '$name'")
    return tags[name]
end


"""
    physicalEntities(msh, name)

Zwraca zbiór geometrycznych tagów encji należących do grupy fizycznej.
"""
function physicalEntities(msh::AbstractDict, name::AbstractString)
    entities = msh["physical entities"]
    haskey(entities, name) || error("Brak grupy fizycznej o nazwie '$name'")
    return entities[name]
end


"""
    physicalFaces(msh, name)

Zwraca indeksy wierszy `msh["faces"]`, które należą do wskazanej
fizycznej grupy powierzchniowej.
"""
function physicalFaces(msh::AbstractDict, name::AbstractString)
    dim, _ = physicalTag(msh, name)
    dim == 2 || error("Grupa '$name' ma wymiar $dim, a nie 2")

    entities = physicalEntities(msh, name)
    tags = @view msh["faces"][:,1]
    return findall(tag -> tag in entities, tags)
end


"""
    physicalElements(msh, name)

Zwraca indeksy czworościanów należących do wskazanej fizycznej
grupy objętościowej.
"""
function physicalElements(msh::AbstractDict, name::AbstractString)
    dim, _ = physicalTag(msh, name)
    dim == 3 || error("Grupa '$name' ma wymiar $dim, a nie 3")

    entities = physicalEntities(msh, name)
    return findall(tag -> tag in entities, msh["s"])
end


"""
    physicalFaceNodes(msh, name)

Zwraca posortowany wektor unikalnych numerów węzłów należących do
wskazanej fizycznej grupy powierzchniowej.
"""
function physicalFaceNodes(msh::AbstractDict, name::AbstractString)
    face_ids = physicalFaces(msh, name)
    faces = msh["faces"]

    nodes = Set{Int64}()
    sizehint!(nodes, 3 * length(face_ids))

    for f in face_ids
        push!(nodes, faces[f,2])
        push!(nodes, faces[f,3])
        push!(nodes, faces[f,4])
    end

    result = collect(nodes)
    sort!(result)
    return result
end


"""
    printPhysicalGroups(msh)

Wypisuje nazwy grup fizycznych, ich wymiar, tag fizyczny i tagi
geometrycznych encji.
"""
function printPhysicalGroups(msh::AbstractDict)
    tags = msh["physical tags"]
    entities = msh["physical entities"]

    names = sort!(collect(keys(tags)))

    for name in names
        dim, physical_tag = tags[name]
        geometry_tags = sort!(collect(entities[name]))
        println(
            name,
            ": dim=", dim,
            ", physical tag=", physical_tag,
            ", geometrical tags=", geometry_tags
        )
    end

    return nothing
end
