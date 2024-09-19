# Minimum spanning tree using Prim's algorithm
# Implementation uses Fibonacci heap priority queue

# input matrices from text files
using DelimitedFiles

# import Fibonacci heap data structure
include("fibonacci-heap.jl")

# Record to store in the Fibonacci heap
mutable struct VertRec

    # data fields
    id::Int64                               # the id of the vertex
    dist::Float64                           # minimum weight of an edge connecting
                                            # the vertex with the MST constructed
                                            # the key value for the Fibonacci heap

    # Constructor
    # id            id of the vertex
    # dist          shortest distance estimate
    # return        the constructed record
    function VertRec(id::Int64, dist::Float64)
        vr = new(id, dist)

        return vr
    end
end

# overload less than operator for the Fibonacci heap key value
function Base.:<(vr_frst::VertRec, vr_scnd::VertRec)
    return vr_frst.dist < vr_scnd.dist
end

# Represent a network using adjacency list
# in_mtrx       input adjacency matrix of weights
# return        the adjacency list of outgoing edges
function adjlist(in_mtrx::Array{Int64})
    result = fill(Vector{Pair{Int64, Int64}}([]), size(in_mtrx, 1))

    for j in eachindex(in_mtrx[1, :]), i in eachindex(in_mtrx[:, 1])
        if in_mtrx[i, j] != 0
            edge = Pair(j, in_mtrx[i, j])
            if size(result[i], 1) == 0
                # the adjacency list for vertex i is empty, create a list of single 
                # element
                result[i] = [edge]
            else
                # the adjacency list for vertex i is not empty, push the edge
                push!(result[i], edge)
            end
        end
    end

    return result
end

# Initialize the shortest path estimates and predecessors vectors
# numb          number of vertices in the graph
# srce          source vertex
# return        dist   vector of min weight of a connecting edge with constructing MST
#               pred   vector of predecessors
#               vstd   Boolean array that stores whether vertex is visited 
function initsingsrc(numb::Int64, srce::Int64)
    dist = Vector{Union{Nothing, Float64}}(nothing, numb)
    pred = Vector{Union{Nothing, Int64}}(nothing, numb)
    vstd = Array{Bool}(undef, numb)

    for i in eachindex(dist)
        dist[i] = Inf
        vstd[i] = false
    end
    dist[srce] = 0

    return dist, pred, vstd
end

# Prim's algorithm
# adj_list      adjacency list of the graph
# srce          root vertex
function prim(
        adj_list::Vector{Vector{Pair{Int64, Int64}}},
        srce::Int64
    )
    # number of vertices in the graph
    numb = size(adj_list[:, 1], 1)

    # shortest distance estimates and predecessors vectors
    dist, pred, vstd = initsingsrc(numb, srce)

    # Fibonacci heap queued by the shortest path estimates of vertices
    fheap = FibonacciHeap{VertRec}()

    # table that stores the nodes of the heap indexed by vertex id
    table = Vector{Union{Nothing, FHeapNode{VertRec}}}(nothing, numb)

    for i = 1:numb
        vr = VertRec(i, dist[i])
        table[i] = insert!(fheap, vr)
    end

    # main loop of the algorithm
    while !isempty(fheap)

        # extract the minimum node from priority queue
        node = extractmin!(fheap)
        sta_vrtx = node.key.id      # the id of the current min vertex

        # mark node as visited
        vstd[sta_vrtx] = true

        # visit each adjacent vertex of the current
        for adj in adj_list[sta_vrtx]
            end_vrtx = adj.first    # vertex id
            wght = adj.second       # edge weight

            if vstd[end_vrtx] == false && wght < dist[end_vrtx]
                
                # update min weight of a connecting edge in the vector
                dist[end_vrtx] = wght

                # reorder the Fibonacci heap by decreasekey function
                new_key = VertRec(end_vrtx, wght * 1.0)
                decreasekey!(fheap, table[end_vrtx], new_key)

                # set predecessor
                pred[end_vrtx] = sta_vrtx
            end
        end
    end

    return dist, pred
end

# input adjacency matrix of the graph
in_adj_mtrx = readdlm("mtrx10.txt", Int64)
println("Graph input adjacency matrix:")
display(in_adj_mtrx)
println()

# represent the graph with adjacency list
adj_list = adjlist(in_adj_mtrx)
println("Graph adjacency list:")
display(adj_list)
println()

# Prims's algorithm
dist, pred = prim(adj_list, 1)

println("Vector of minimum weights:")
display(dist)
println("Total weight of a MST: ", sum(dist))
println("Vector of predecessors:")
display(pred)
