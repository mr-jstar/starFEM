 ============================================================
# Strumieniowy importer Gmsh MSH ASCII 2.2 i 4.1
#
# Wynik:
#
# msh["x"]              :: Matrix{Float64}   nn × 3
# msh["v"]              :: Matrix{Int64}     ne × 4
# msh["faces"]          :: Matrix{Int64}     nf × 4
#                         kolumny: physical_tag, n1, n2, n3
# msh["s"]              :: Vector{Int64}     ne
# msh["b"]              :: Vector{Int64}     nn
# msh["physical names"] :: Dict{Tuple{Int64,Int64},String}
#
# Obsługiwane elementy:
#   typ 2: trójkąt 3-węzłowy
#   typ 4: czworościan 4-węzłowy
# ============================================================


# ------------------------------------------------------------
# Pomocnicze operacje wejścia
# ------------------------------------------------------------

"""
    next_nonempty_line(io)

Czyta następną niepustą linię. Zwraca `nothing` po osiągnięciu EOF.
"""
function next_nonempty_line(io::IO)
    while !eof(io)
        line = strip(readline(io))
        isempty(line) || return line
    end
    return nothing
end


"""
    expect_line(io, expected)

Czyta następną niepustą linię i sprawdza jej zawartość.
"""
function expect_line(io::IO, expected::AbstractString)
    line = next_nonempty_line(io)

    line === nothing &&
        error("Unexpected end of file; expected $expected")

    line == expected ||
        error("Expected $expected, got: $line")

    return nothing
end


"""
    skip_section!(io, section_start)

Pomija nieobsługiwaną sekcję, np. `$Periodic`.
"""
function skip_section!(io::IO, section_start::AbstractString)
    startswith(section_start, "\$") ||
        error("Invalid section marker: $section_start")

    section_name = section_start[2:end]
    end_marker = "\$End" * section_name

    while !eof(io)
        strip(readline(io)) == end_marker && return nothing
    end

    error("Missing $end_marker")
end


"""
    read_int_tokens(io, n)

Czyta dokładnie `n` liczb całkowitych. Liczby mogą znajdować się
w jednej lub kilku kolejnych liniach.
"""
function read_int_tokens(io::IO, n::Integer)
    result = Vector{Int64}(undef, n)
    k = 1

    while k <= n
        eof(io) && error("Unexpected EOF while reading integer values")

        fields = split(strip(readline(io)))

        for field in fields
            k > n &&
                error("Too many integer values in Gmsh block")

            result[k] = parse(Int64, field)
            k += 1
        end
    end

    return result
end


# ------------------------------------------------------------
# $MeshFormat
# ------------------------------------------------------------

function read_mesh_format!(io::IO)
    fields = split(strip(readline(io)))

    length(fields) >= 3 ||
        error("Invalid \$MeshFormat line")

    version = VersionNumber(fields[1])
    file_type = parse(Int, fields[2])
    data_size = parse(Int, fields[3])

    file_type == 0 ||
        error("Binary MSH files are not supported; save the mesh as ASCII")

    expect_line(io, "\$EndMeshFormat")

    return version, data_size
end


# ------------------------------------------------------------
# $PhysicalNames
# ------------------------------------------------------------

function read_physical_names!(io::IO)
    names = Dict{Tuple{Int64,Int64},String}()

    n = parse(Int, strip(readline(io)))
    sizehint!(names, n)

    for _ in 1:n
        line = strip(readline(io))

        # Nazwa może zawierać spacje.
        m = match(
            r"^\s*(-?\d+)\s+(-?\d+)\s+\"(.*)\"\s*$",
            line
        )

        m === nothing &&
            error("Invalid \$PhysicalNames entry: $line")

        dim = parse(Int64, m.captures[1])
        tag = parse(Int64, m.captures[2])
        name = m.captures[3]

        names[(dim, tag)] = name
    end

    expect_line(io, "\$EndPhysicalNames")

    return names
end


# ------------------------------------------------------------
# $Entities dla MSH 4.x
#
# Wynik:
#   (dimension, entity tag) => lista physical tags
# ------------------------------------------------------------

function read_entities41!(io::IO)
    entity_physical =
        Dict{Tuple{Int64,Int64},Vector{Int64}}()

    header = split(strip(readline(io)))

    length(header) == 4 ||
        error("Invalid \$Entities header")

    npoints   = parse(Int, header[1])
    ncurves   = parse(Int, header[2])
    nsurfaces = parse(Int, header[3])
    nvolumes  = parse(Int, header[4])

    # Punkty:
    # pointTag X Y Z numPhysicalTags physicalTag ...
    for _ in 1:npoints
        fields = split(strip(readline(io)))

        length(fields) >= 5 ||
            error("Invalid point entity")

        entity_tag = parse(Int64, fields[1])
        nphys = parse(Int, fields[5])

        length(fields) >= 5 + nphys ||
            error("Invalid physical tags for point entity $entity_tag")

        if nphys > 0
            tags = Vector{Int64}(undef, nphys)

            for k in 1:nphys
                tags[k] = parse(Int64, fields[5 + k])
            end

            entity_physical[(Int64(0), entity_tag)] = tags
        end
    end

    # Krzywe, powierzchnie i objętości:
    #
    # entityTag minX minY minZ maxX maxY maxZ
    # numPhysicalTags physicalTag ...
    # numBoundingEntities ...
    for (dim, number_of_entities) in (
        (Int64(1), ncurves),
        (Int64(2), nsurfaces),
        (Int64(3), nvolumes)
    )
        for _ in 1:number_of_entities
            fields = split(strip(readline(io)))

            length(fields) >= 8 ||
                error("Invalid entity definition for dimension $dim")

            entity_tag = parse(Int64, fields[1])
            nphys = parse(Int, fields[8])

            length(fields) >= 8 + nphys ||
                error("Invalid physical tags for entity ($dim,$entity_tag)")

            if nphys > 0
                tags = Vector{Int64}(undef, nphys)

                for k in 1:nphys
                    tags[k] = parse(Int64, fields[8 + k])
                end

                entity_physical[(dim, entity_tag)] = tags
            end
        end
    end

    expect_line(io, "\$EndEntities")

    return entity_physical
end


# ------------------------------------------------------------
# Węzły MSH 2.2
# ------------------------------------------------------------

function read_nodes22!(io::IO)
    number_of_nodes = parse(Int, strip(readline(io)))

    node_tags = Vector{Int64}(undef, number_of_nodes)
    coordinates = Matrix{Float64}(undef, number_of_nodes, 3)

    for local_node in 1:number_of_nodes
        fields = split(strip(readline(io)))

        length(fields) >= 4 ||
            error("Invalid node record in MSH 2.x")

        node_tags[local_node] = parse(Int64, fields[1])
        coordinates[local_node,1] = parse(Float64, fields[2])
        coordinates[local_node,2] = parse(Float64, fields[3])
        coordinates[local_node,3] = parse(Float64, fields[4])
    end

    expect_line(io, "\$EndNodes")

    return node_tags, coordinates
end


# ------------------------------------------------------------
# Węzły MSH 4.1
# ------------------------------------------------------------

function read_nodes41!(io::IO)
    header = split(strip(readline(io)))

    length(header) == 4 ||
        error("Invalid \$Nodes header for MSH 4.x")

    number_of_blocks = parse(Int, header[1])
    number_of_nodes = parse(Int, header[2])

    node_tags = Vector{Int64}(undef, number_of_nodes)
    coordinates = Matrix{Float64}(undef, number_of_nodes, 3)

    local_node = 0

    for _ in 1:number_of_blocks
        block_header = split(strip(readline(io)))

        length(block_header) == 4 ||
            error("Invalid node block header")

        entity_dim = parse(Int, block_header[1])
        # entity_tag = parse(Int64, block_header[2])
        parametric = parse(Int, block_header[3])
        nodes_in_block = parse(Int, block_header[4])

        block_tags = read_int_tokens(io, nodes_in_block)

        for k in 1:nodes_in_block
            fields = split(strip(readline(io)))

            # Przy węzłach parametrycznych po xyz występują
            # jeszcze współrzędne parametryczne. Nie są potrzebne.
            minimum_fields = 3 + (parametric == 1 ? entity_dim : 0)

            length(fields) >= minimum_fields ||
                error("Invalid coordinate record in node block")

            local_node += 1

            node_tags[local_node] = block_tags[k]
            coordinates[local_node,1] = parse(Float64, fields[1])
            coordinates[local_node,2] = parse(Float64, fields[2])
            coordinates[local_node,3] = parse(Float64, fields[3])
        end
    end

    local_node == number_of_nodes ||
        error(
            "Expected $number_of_nodes nodes, read $local_node"
        )

    expect_line(io, "\$EndNodes")

    return node_tags, coordinates
end


# ------------------------------------------------------------
# Elementy MSH 2.2
#
# W tym formacie pierwszy tag elementu jest physical tagiem.
# ------------------------------------------------------------

function read_elements22!(io::IO)
    number_of_elements = parse(Int, strip(readline(io)))

    tet_nodes = NTuple{4,Int64}[]
    tet_physical = Int64[]

    face_nodes = NTuple{3,Int64}[]
    face_physical = Int64[]

    for _ in 1:number_of_elements
        fields = split(strip(readline(io)))

        length(fields) >= 3 ||
            error("Invalid element record in MSH 2.x")

        element_type = parse(Int, fields[2])
        number_of_tags = parse(Int, fields[3])

        physical_tag =
            number_of_tags > 0 ? parse(Int64, fields[4]) : Int64(0)

        first_node_field = 4 + number_of_tags

        if element_type == 2
            # 3-węzłowy trójkąt
            length(fields) >= first_node_field + 2 ||
                error("Invalid 3-node triangle")

            push!(
                face_nodes,
                (
                    parse(Int64, fields[first_node_field]),
                    parse(Int64, fields[first_node_field + 1]),
                    parse(Int64, fields[first_node_field + 2])
                )
            )

            push!(face_physical, physical_tag)

        elseif element_type == 4
            # 4-węzłowy czworościan
            length(fields) >= first_node_field + 3 ||
                error("Invalid 4-node tetrahedron")

            push!(
                tet_nodes,
                (
                    parse(Int64, fields[first_node_field]),
                    parse(Int64, fields[first_node_field + 1]),
                    parse(Int64, fields[first_node_field + 2]),
                    parse(Int64, fields[first_node_field + 3])
                )
            )

            push!(tet_physical, physical_tag)
        end
    end

    expect_line(io, "\$EndElements")

    return tet_nodes, tet_physical, face_nodes, face_physical
end


# ------------------------------------------------------------
# Elementy MSH 4.1
#
# Fizyczny tag nie występuje bezpośrednio przy elemencie.
# Zapamiętujemy encję, do której należy blok elementów.
# ------------------------------------------------------------

function read_elements41!(io::IO)
    header = split(strip(readline(io)))

    length(header) == 4 ||
        error("Invalid \$Elements header for MSH 4.x")

    number_of_blocks = parse(Int, header[1])
    # total_elements = parse(Int, header[2])

    tet_nodes = NTuple{4,Int64}[]
    tet_entities = Tuple{Int64,Int64}[]

    face_nodes = NTuple{3,Int64}[]
    face_entities = Tuple{Int64,Int64}[]

    for _ in 1:number_of_blocks
        block_header = split(strip(readline(io)))

        length(block_header) == 4 ||
            error("Invalid element block header")

        entity_dim = parse(Int64, block_header[1])
        entity_tag = parse(Int64, block_header[2])
        element_type = parse(Int, block_header[3])
        elements_in_block = parse(Int, block_header[4])

        if element_type == 2
            # trójkąt liniowy
            for _ in 1:elements_in_block
                fields = split(strip(readline(io)))

                # elementTag n1 n2 n3
                length(fields) >= 4 ||
                    error("Invalid 3-node triangle")

                push!(
                    face_nodes,
                    (
                        parse(Int64, fields[2]),
                        parse(Int64, fields[3]),
                        parse(Int64, fields[4])
                    )
                )

                push!(face_entities, (entity_dim, entity_tag))
            end

        elseif element_type == 4
            # czworościan liniowy
            for _ in 1:elements_in_block
                fields = split(strip(readline(io)))

                # elementTag n1 n2 n3 n4
                length(fields) >= 5 ||
                    error("Invalid 4-node tetrahedron")

                push!(
                    tet_nodes,
                    (
                        parse(Int64, fields[2]),
                        parse(Int64, fields[3]),
                        parse(Int64, fields[4]),
                        parse(Int64, fields[5])
                    )
                )

                push!(tet_entities, (entity_dim, entity_tag))
            end

        else
            # Nieobsługiwany typ elementu: pomijamy blok.
            for _ in 1:elements_in_block
                eof(io) &&
                    error("Unexpected EOF inside element block")
                readline(io)
            end
        end
    end

    expect_line(io, "\$EndElements")

    return tet_nodes, tet_entities, face_nodes, face_entities
end


# ------------------------------------------------------------
# Mapowanie fizycznych tagów w MSH 4.1
# ------------------------------------------------------------

function entity_physical_tag(
    entity_physical::Dict{
        Tuple{Int64,Int64},
        Vector{Int64}
    },
    entity::Tuple{Int64,Int64}
)
    tags = get(entity_physical, entity, nothing)

    tags === nothing && return Int64(0)
    isempty(tags) && return Int64(0)

    # Element może teoretycznie należeć przez encję do kilku
    # grup fizycznych, ale s i faces mają pojedynczy tag.
    # Przyjmujemy pierwszy tag zapisany przez Gmsh.
    return tags[1]
end


# ------------------------------------------------------------
# Budowa mapy tag węzła Gmsh -> indeks w tablicy Julii
# ------------------------------------------------------------

function make_node_index(node_tags::Vector{Int64})
    tag_to_index = Dict{Int64,Int64}()
    sizehint!(tag_to_index, length(node_tags))

    for local_index in eachindex(node_tags)
        tag = node_tags[local_index]

        haskey(tag_to_index, tag) &&
            error("Duplicate Gmsh node tag: $tag")

        tag_to_index[tag] = Int64(local_index)
    end

    return tag_to_index
end


# ------------------------------------------------------------
# Składanie wynikowych macierzy
# ------------------------------------------------------------

function assemble_gmsh_mesh(
    coordinates::Matrix{Float64},
    node_tags::Vector{Int64},
    tet_node_tags::Vector{NTuple{4,Int64}},
    tet_physical::Vector{Int64},
    face_node_tags::Vector{NTuple{3,Int64}},
    face_physical::Vector{Int64},
    physical_names::Dict{Tuple{Int64,Int64},String}
)
    tag_to_index = make_node_index(node_tags)

    number_of_tets = length(tet_node_tags)
    number_of_faces = length(face_node_tags)
    number_of_nodes = size(coordinates, 1)

    length(tet_physical) == number_of_tets ||
        error("Inconsistent number of tetrahedral physical tags")

    length(face_physical) == number_of_faces ||
        error("Inconsistent number of face physical tags")

    v = Matrix{Int64}(undef, number_of_tets, 4)
    s = Vector{Int64}(undef, number_of_tets)

    for e in 1:number_of_tets
        nodes = tet_node_tags[e]

        v[e,1] = get_node_index(tag_to_index, nodes[1])
        v[e,2] = get_node_index(tag_to_index, nodes[2])
        v[e,3] = get_node_index(tag_to_index, nodes[3])
        v[e,4] = get_node_index(tag_to_index, nodes[4])

        s[e] = tet_physical[e]
    end

    faces = Matrix{Int64}(undef, number_of_faces, 4)

    for f in 1:number_of_faces
        nodes = face_node_tags[f]

        faces[f,1] = face_physical[f]
        faces[f,2] = get_node_index(tag_to_index, nodes[1])
        faces[f,3] = get_node_index(tag_to_index, nodes[2])
        faces[f,4] = get_node_index(tag_to_index, nodes[3])
    end

    # Pojedynczy węzeł może należeć do kilku powierzchni.
    # Ponieważ b ma tylko jeden Int na węzeł, zapisujemy
    # najmniejszy niezerowy tag powierzchni.
    b = zeros(Int64, number_of_nodes)

    for f in 1:number_of_faces
        physical_tag = faces[f,1]
        physical_tag == 0 && continue

        for j in 2:4
            node = faces[f,j]

            if b[node] == 0 || physical_tag < b[node]
                b[node] = physical_tag
            end
        end
    end

    return Dict{String,Any}(
        "x"              => coordinates,
        "v"              => v,
        "faces"          => faces,
        "s"              => s,
        "b"              => b,
        "physical names" => physical_names
    )
end


@inline function get_node_index(
    tag_to_index::Dict{Int64,Int64},
    tag::Int64
)
    index = get(tag_to_index, tag, Int64(0))

    index != 0 ||
        error("Element refers to unknown Gmsh node tag $tag")

    return index
end


# ------------------------------------------------------------
# Główna funkcja
# ------------------------------------------------------------

function loadGmsh(filename::AbstractString)
    version = nothing

    physical_names =
        Dict{Tuple{Int64,Int64},String}()

    entity_physical =
        Dict{Tuple{Int64,Int64},Vector{Int64}}()

    node_tags = Int64[]
    coordinates = Matrix{Float64}(undef, 0, 3)

    tet_node_tags = NTuple{4,Int64}[]
    face_node_tags = NTuple{3,Int64}[]

    # Dla MSH 2.x są to od razu tagi fizyczne.
    tet_physical = Int64[]
    face_physical = Int64[]

    # Dla MSH 4.x przechowujemy encje do późniejszego mapowania.
    tet_entities = Tuple{Int64,Int64}[]
    face_entities = Tuple{Int64,Int64}[]

    open(filename, "r") do io
        while !eof(io)
            marker = next_nonempty_line(io)
            marker === nothing && break

            if marker == "\$MeshFormat"
                version, _ = read_mesh_format!(io)

            elseif marker == "\$PhysicalNames"
                physical_names = read_physical_names!(io)

            elseif marker == "\$Entities"
                version === nothing &&
                    error("\$Entities encountered before \$MeshFormat")

                version.major >= 4 ||
                    error("\$Entities section in unsupported MSH version")

                entity_physical = read_entities41!(io)

            elseif marker == "\$Nodes"
                version === nothing &&
                    error("\$Nodes encountered before \$MeshFormat")

                if version.major == 2
                    node_tags, coordinates = read_nodes22!(io)
                elseif version.major == 4
                    node_tags, coordinates = read_nodes41!(io)
                else
                    error("Unsupported MSH version $version")
                end

            elseif marker == "\$Elements"
                version === nothing &&
                    error("\$Elements encountered before \$MeshFormat")

                if version.major == 2
                    tet_node_tags,
                    tet_physical,
                    face_node_tags,
                    face_physical = read_elements22!(io)

                elseif version.major == 4
                    tet_node_tags,
                    tet_entities,
                    face_node_tags,
                    face_entities = read_elements41!(io)

                else
                    error("Unsupported MSH version $version")
                end

            elseif startswith(marker, "\$") &&
                   !startswith(marker, "\$End")
                skip_section!(io, marker)

            else
                error("Unexpected content outside a section: $marker")
            end
        end
    end

    version === nothing &&
        error("Missing \$MeshFormat section")

    isempty(node_tags) &&
        error("No nodes found in $filename")

    if version.major == 4
        tet_physical = Vector{Int64}(undef, length(tet_entities))

        for i in eachindex(tet_entities)
            tet_physical[i] =
                entity_physical_tag(entity_physical, tet_entities[i])
        end

        face_physical =
            Vector{Int64}(undef, length(face_entities))

        for i in eachindex(face_entities)
            face_physical[i] =
                entity_physical_tag(entity_physical, face_entities[i])
        end
    end

    return assemble_gmsh_mesh(
        coordinates,
        node_tags,
        tet_node_tags,
        tet_physical,
        face_node_tags,
        face_physical,
        physical_names
    )
end