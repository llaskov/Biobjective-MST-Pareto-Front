# Solves the minimum length minimum risk minimum spanning trees problem
# Computes the minimum complete Pareto optimal set

# include
include("network.jl")                       # Representation of network
include("fibonacci-heap.jl")                # Fibonacci heap data structure

# Record to store in the Fibonacci heap
mutable struct VertRec

    # data fields
    id::Int64                               # the id of the vertex
    dist::Array{Int64}                      # array of two components: length, risk 

    # Constructor
    # id            id of the vertex
    # dist          array of two components: length, risk
    # return        the constructed record
    function VertRec(id::Int64, dist::Array{Int64})
        vr = new(id, dist)

        return vr
    end
end

# overload less than operator for the Fibonacci heap key value
function Base.:<(vr_frst::VertRec, vr_scnd::VertRec)
    return vr_frst.dist < vr_scnd.dist
end

# Structure to store a spanning tree with its weight values
mutable struct SpanTree

    # data fields
    pred::Vector{Union{Nothing, Int64}}     # list of predecessors
    lngth::Int64                            # total length of the tree
    risk::Int64                             # total risk of the tree

    # Constructor
    # pred          list of predecessors
    # lngth         total length of the tree
    # risk          total risk of the tree    
    # return        the constructed record
    function SpanTree(
            pred::Vector{Union{Nothing, Int64}}, 
            lngth::Int64, 
            risk::Int64
        )
        st = new(pred, lngth, risk)

        return st
    end
end

# Display a spanning tree
# st        spanning tree
function Base.display(st::SpanTree)
    println("Predecesors list: ", st.pred)
    println("Totoal length: ", st.lngth, "\tTotal risk: ", st.risk)
end

# Initialize the shortest distance estimates and predecessors vectors
# numb          number of vertices in the graph
# srce          source vertex
# return        dist   vector of pairs determine distance to constructing MST
#               pred   vector of predecessors
#               vstd   Boolean array that stores whether vertex is visited 
function initsingsrc(
        numb::Int64,
        srce::Int64
    )
    dist = Vector{Union{Nothing, Array{Int64}}}(nothing, numb)
    pred = Vector{Union{Nothing, Int64}}(nothing, numb)
    vstd = BitArray(undef,  numb)

    for i in eachindex(dist)
        dist[i] = [typemax(Int64), typemax(Int64)]
        vstd[i] = false
    end
    dist[srce] = [0, 0]

    return dist, pred, vstd
end

# Calculate the two weight total values of a biobjective MST
# dist          vector of pairs that determine distance attribute of the MST
# return        length  total length of the MST
#               risk    total risk of the MST
function mstweights(dist::Vector{Union{Nothing, Array{Int64}}})
    lngth::Int64 = 0
    risk::Int64 = 0

    for wght_pair in dist
        lngth += wght_pair[1]
        risk = max(risk, wght_pair[2]) 
    end

    return lngth, risk
end

# Find minimum length spanning tree with minimum risk
# Extension of Prim's algorithm
# alnet         adjacency list of the input network
# root          root vertex
# return        dist    vector of pairs determine distance to constructing MST
#               pred    vector of predecessors
function lengthrisk(
        alnet::AdjLst{NtwrkRec},
        root::Int64
    )
    numb = size(alnet)                      # number of vertices in the network

    # shortest distance estimates and predecessors vectors
    dist, pred, vstd = initsingsrc(numb, root)

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
        for adj in alnet.list[sta_vrtx]
            end_vrtx = adj.vrtx_end # vertex id
            wght = [adj.wght_one, adj.wght_two] # edge weight pair

            if vstd[end_vrtx] == false && wght < dist[end_vrtx]
                
                # update min weight of a connecting edge in the vector
                dist[end_vrtx] = wght

                # reorder the Fibonacci heap by decreasekey function
                new_key = VertRec(end_vrtx, wght)
                decreasekey!(fheap, table[end_vrtx], new_key)

                # set predecessor
                pred[end_vrtx] = sta_vrtx
            end
        end
    end
    
    # construct the MST structure to return
    lngth, risk = mstweights(dist)
    st = SpanTree(pred, lngth, risk)
    
    return st
end

# Visit each reachable vertex of the network using DFS
# srce          source vertex
# alnet         adjacency list of the input network
# vstd          a bit array that stores for each vertex whether it is visited
function dfsvisit!(
        srce::Int64,
        alnet::AdjLst{NtwrkRec},
        vstd::BitArray
    )
    vstd[srce] = true

    # visit each adjacent vertex of the current
    for adj in alnet.list[srce]
        if vstd[adj.vrtx_end] == false
            dfsvisit!(adj.vrtx_end, alnet, vstd)
        end
    end    
end

# Verify weather a network is connected using DFS
# alnet         adjacency list of the input network
# return        true    if the network is connected
#               false   if the network is not connected
function isconnect(alnet::AdjLst{NtwrkRec})

    # bit array to store for each vertex whether it is visited
    vstd = BitArray(undef,  numbvert(alnet))

    for i in eachindex(vstd)
        vstd[i] = false
    end

    dfsvisit!(1, al_in_net, vstd)           # call DFS traverse of the network
    result = true                           # if all vertices are visited, true
    for vrtx in vstd
        if vrtx == false                    # just one is not visited, false
            result = false
            break
        end
    end

    return result
end

# Traverse the edges of the input network and compose a subnetwork with edges
# that have risk less than a given boundary risk
# alnet         adjacency list of the input network, will be modified
# risk          boundary risk
function restrict!(
        alnet::AdjLst{NtwrkRec},
        risk::Int64
    )
    for vrtx_indx = 1:numbvert(alnet)
        filter!(adj -> adj.wght_two < risk, alnet.list[vrtx_indx])
    end
end

# Compose the minimum complete Pareto front of the biobjective MSTs
# alnet         adjacency list of the input network, will be modified
# root          the root vertex
# return        pf      list of the minimum complete Pareto front
function mincomplpf(
        alnet::AdjLst{NtwrkRec},
        root::Int64
    )
    pf = Vector{SpanTree}()                 # vector to contain the Pareto front
    while isconnect(alnet)
        st = lengthrisk(alnet, root)
        push!(pf, st)
        restrict!(alnet, st.risk)
    end

    return pf
end

# input
println("--> Input network:")
al_in_net = readadj("net100_hard.txt")       # input network adjacency list
println("Size (number of vertics): ", numbvert(al_in_net))
println("Size (number of edges): ", numbedge(al_in_net))

# run minimum Pareto optimal front algorithm
println("--> Minimum Pareto optimal front:")
@time pf = mincomplpf(al_in_net, 1)
display(pf)
