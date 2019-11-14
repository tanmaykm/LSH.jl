module LSH

using Primes
using SparseArrays
using Statistics

import Base: match, show

export LSHModel
export match, signature

const DEFAULT_THRESHOLD = 0.75

struct LSHModel{T <: Integer}
    A::Vector{T}
    B::Vector{T}
    c::T
    termhashes::Array{T,2}
    dochashes::Array{T,2}
    bandhash::Dict{Float64,Any}
end

function LSHModel(nhashes::Integer, dtm::SparseMatrixCSC{T,T}) where {T <: Integer}
    ndocs,nterms = size(dtm)
    c = T(nextprime(nterms))
    A = rand(T(1):T(c), nhashes)
    B = rand(T(1):T(c), nhashes)

    termhashes = Array{T,2}(undef,nterms,nhashes)
    dochashes = Array{T,2}(undef,ndocs,nhashes)

    termhash!(termhashes, A, B, c)
    apply_term_weights!(termhashes, dtm)
    dochash!(dochashes, termhashes, dtm)

    LSHModel{T}(A, B, c, termhashes, dochashes, Dict{Float64,Any}())
end

function show(io::IO, model::LSHModel{T}) where T
    nterms,nhashes = size(model.termhashes)
    ndocs,_ = size(model.dochashes)
    print(io, "LSHModel{$T}(nhashes=$nhashes, ndocs=$ndocs, nterms=$nterms)")
end

# Weighted Minhash
function termhash!(termhashes::Array{T,2}, A::Vector{T}, B::Vector{T}, c::T) where {T <: Integer}
    @info("generating term hashes...")
    nterms,nhashes = size(termhashes)
    terms = 1:nterms
    for h in 1:nhashes
        a = A[h]
        b = B[h]
        termhashes[:,h] = (x->T((a*x+b)%c)).(terms)
    end
    nothing
end

function apply_term_weights!(termhashes::Array{T,2}, dtm::SparseMatrixCSC{T,T}) where {T <: Integer}
    @info("applying term weights...")
    ndocs,nterms = size(dtm)
    twf = [log(1+ndocs/sum(dtm[:,t])) for t in 1:nterms]
    mw = minimum(twf)
    termweights = [ceil(T, w/mw) for w in twf]
    termhashes .*= termweights
    nothing
end

function dochash!(dochashes::Array{T,2}, termhashes::Array{T,2}, dtm::SparseMatrixCSC{T,T}) where {T <: Integer}
    @info("generating signatures...")
    ndocs,nhashes = size(dochashes)

    for h in 1:nhashes
        th = termhashes[:,h]
        for d in 1:ndocs
            tm = dtm[d,:]
            nzinds = SparseArrays.nonzeroinds(tm)
            dochashes[d,h] = minimum(th[nzinds] .* tm[nzinds])
        end
    end
    nothing
end

### LSH
similarity(sig1::Vector{T}, sig2::Vector{T}) where {T <: Integer} = mean(.==(sig1, sig2))
threshold(nbands::T, nrows::T) where {T <: Integer} = (1/nbands) ^ (1/nrows)

function factors(n::T) where {T <: Integer}
    allfactors = T[]

    for i in T(1):n
        if n % i == 0
            push!(allfactors, i)
        end
    end

    allfactors
end

# nbands and nrows are across hash algorithms
function best_partition(nhashes::T, threslim::Float64) where {T <: Integer}
    allfactors = factors(nhashes)
    mindiff = Inf
    bestpart = (T(0),T(0))

    for nbands in allfactors
        nrows = div(nhashes, nbands)
        t = threshold(nbands, nrows)
        tdiff = abs(t - threslim)
        if tdiff < mindiff
            mindiff = tdiff
            bestpart = (nbands, nrows)
        end
    end

    bestpart
end

function bandhash(model::LSHModel{T}, threshold::Float64) where {T <: Integer}
    get!(model.bandhash, threshold) do
        dochashes = model.dochashes
        ndocs,nhashes = size(dochashes)
        nbands,nrows = best_partition(T(nhashes), threshold)

        ht = Dict{Vector{T},Vector{T}}()

        for d in T(1):ndocs
            dh = dochashes[d,:]
            for b in 1:nbands
                brange = (1:nrows) .+ ((b-1)*nrows)
                bsig = dh[brange]
                bsigdocs = get!(ht, bsig) do
                    T[]
                end
                append!(bsigdocs, d)
            end
        end
        return ht
    end
end

match(model::LSHModel{T}, sig::Array{T,2}, threshold::Float64=DEFAULT_THRESHOLD) where {T <: Integer} = [match(model, sig[d,:], threshold) for d in 1:size(sig,1)]

function match(model::LSHModel{T}, sig::Vector{T}, threshold::Float64=DEFAULT_THRESHOLD) where {T <: Integer}
    nhashes = size(model.dochashes, 2)
    @assert(length(sig) == nhashes)

    nbands,nrows = best_partition(nhashes, threshold)
    bh = bandhash(model, threshold)
    similars = Set{Int}()

    for h in 1:nhashes
        for b in 1:nbands
            brange = (1:nrows) .+ ((b-1)*nrows)
            bsig = sig[brange]
            (x->push!(similars,x)).(get(bh, bsig, Int[]))
        end
    end
    result = []
    for d in similars
        sim = similarity(sig, model.dochashes[d,:])
        (sim < threshold) && continue
        push!(result, (similarity=sim, id=d))
    end
    sort!(result, by=x->x.similarity; rev=true)
end

function signature(model::LSHModel{T}, tm::SparseVector{T,T}) where {T <: Integer}
    termhashes = model.termhashes
    nterms,nhashes = size(termhashes)
    nzinds = SparseArrays.nonzeroinds(tm)
    T[minimum(termhashes[nzinds,h] .* tm[nzinds]) for h in 1:nhashes]
end

function signature(model::LSHModel{T}, tm::SparseMatrixCSC{T,T}) where {T <: Integer}
    nterms1,nhashes = size(model.termhashes)
    ndocs,nterms2 = size(tm)
    @assert(nterms1 == nterms2)
    dochashes = Array{T,2}(undef, ndocs, nhashes)

    for d in 1:ndocs
        dochashes[d,:] = signature(model, tm[d,:])
    end

    dochashes 
end

function selftest(model::LSHModel{T}) where {T <: Integer}
    @info("running selftest...")
    ndocs,nhashes = size(model.dochashes)
    nfailures = 0
    for d in 1:ndocs
        sig = model.dochashes[d,:]
        matches = match(model, sig)
        (d in [x.id for x in matches]) || (nfailures += 1; @warn("docid did not match", docid=d))
        
        if (length(matches) > 1) && (d != first(matches).id)
            # check if there are multiple matches with same score
            topscore = first(matches).similarity
            docscore = matches[findfirst(x->x.id == d, matches)].similarity
            (docscore == topscore) || (nfailures += 1; @warn("doc not best batch for its own signature", docid=d, matches=matches))
        end
    end
    @info("selftest done", nfailures)
    nfailures
end

end # module LSH
