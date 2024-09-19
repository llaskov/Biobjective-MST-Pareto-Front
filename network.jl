# Represent a network with adjacency list and adjacency matrix

# input matrices from text files
using DelimitedFiles

# Record to store in the adjacency list of the network
mutable struct NtwrkRec

    # data fields
    vrtx_end::Int64                         # end vertex id of the directed edge
    wght_one::Float64                       # first criteria weight value
    wght_two::Float64                       # second criteria weight value

    # Constructor
    # vrtx_end      end vertex id of the directed edge
    # wght_one      first criteria weight value
    # whgt_two      second criteria weight value
    # return        the constructed record
    function NtwrkRec(vrtx_end::Int64, wght_one::Float64, wght_two::Float64)
        ar = new(vrtx_end, wght_one, wght_two)

        return ar
    end
end

# Network adjacency list
mutable struct AdjLst{T}

    # data fields
    list::Vector{Vector{T}}                 # vector of vectors of records

    # constructor
    # return        the constructed empty adjacency list
    function AdjLst{T}() where T
        al = new{T}()
        al.list = fill(Vector{T}([]), 0)
        return al
    end

    # constructor
    # numb          number of lists
    # return        the constructed adjacency list with empty lists allocated
    function AdjLst{T}(numb::Int64) where T
        al = new{T}()
        al.list = fill(Vector{T}([]), numb)

        return al
    end
end

# Number of vertices in the network
# al            adjacency list of the network
Base.size(al::AdjLst) = size(al.list[:, 1], 1)

# Number of vertices in the network
# al            adjacency list of the network
numbvert(al::AdjLst) = size(al.list[:, 1], 1)

# Number of edges in the network
# al            adjacency list of the network
function numbedge(
        al::AdjLst{T}
    ) where T
    result = 0
    for i in eachindex(al.list)
        result += size(al.list[i], 1) 
    end

    return result
end

# Push an adjacent element into the adjacency list
# al            network adjacency list
# vrtx_sta      starting vertex of the edge
# rec           record to insert in the adjacency list
function pushadj!(
        al::AdjLst{T}, 
        vrtx_sta::Int64, 
        rec::T
    ) where T
    if size(al.list[vrtx_sta], 1) == 0
        # the adjacency list for starting vertex is empty, 
        # create a list of single element   
        al.list[vrtx_sta] = [rec]
    else
         # the adjacency list for starting vertex is not empty, push the record
        push!(al.list[vrtx_sta], rec)
    end
end

# Insert in front into the adjacency list an adjacent element 
# al            network adjacency list
# vrtx_sta      starting vertex of the edge
# rec           record to insert in the adjacency list
function prependadj!(
        al::AdjLst{T}, 
        vrtx_sta::Int64, 
        rec::T
    ) where T
    if size(al.list[vrtx_sta], 1) == 0
        # the adjacency list for starting vertex is empty, 
        # create a list of single element   
        al.list[vrtx_sta] = [rec]
    else
         # the adjacency list for starting vertex is not empty, push the record
        prepend!(al.list[vrtx_sta], rec)
    end
end

# Append new list on the adjacent list
pushlst!(al::AdjLst{T}) where T = push!(al.list, Vector{T}([]))

# New adjacent element that substitutes previous
# al            network adjacency list
# vrtx_sta      starting vertex of the edge
# rec           record to insert in the adjacency list
function newadj!(
        al::AdjLst{T}, 
        vrtx_sta::Int64, 
        rec::T
    ) where T
    al.list[vrtx_sta] = [rec]
end

# Maximal vertex number to find the total number of vertices
# in_adjlst_mat matrix that stores the adjacency list
# return        the maximal vertex number
function maxvrtx(in_adjlst_mat::Array{Int64})
    result = 0
    j = 1
    while j < size(in_adjlst_mat[1, :], 1)
        i = 1
        while i < size(in_adjlst_mat[:, 1], 1)
            if in_adjlst_mat[i, j] > result
                result = in_adjlst_mat[i, j]
            end
            i += 1
        end
        j += 3
    end
    
    return result
end

# Read network adjacency list from text file that contains a list in 
# rectangular matrix, the columns are truncated with 0
function readadj(filename::AbstractString)
    in_adjlst_mat = readdlm(filename, Int64)
    numb = maxvrtx(in_adjlst_mat)               # number of vertices in the graph
    alnet = AdjLst{NtwrkRec}(numb)
    for i in eachindex(in_adjlst_mat[:, 1])
        j = 1
        while j <= size(in_adjlst_mat, 2) && in_adjlst_mat[i, j] != 0
            nr = NtwrkRec(
                in_adjlst_mat[i, j],            # ending edge of the vertex
                1.0 * in_adjlst_mat[i, j + 1],  # weight value of the first criteria
                1.0 * in_adjlst_mat[i, j + 2])  # weight value of the second criteria
            pushadj!(
                alnet,                          # adjacency list 
                i,                              # starting vertex of the edge 
                nr)                             # record to insert in adjacency list
            j += 3
        end
    end

    return alnet
end

# Transform network representation from adjacency list to adjacency matrix
# al            network adjacency list
# return        network adjacency matrix
function toadjmtrx(al::AdjLst{NtwrkRec})
    numb = size(al)                         # number of vertices

    adj_mtrx = Array{Union{Nothing, Pair{Float64, Float64}}}(nothing, numb, numb)

    for i in eachindex(al.list), j in eachindex(al.list[i])
        adj_mtrx[i, al.list[i][j].vrtx_end] = 
            Pair(al.list[i][j].wght_one, al.list[i][j].wght_two)
    end

    return adj_mtrx
end
