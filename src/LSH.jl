module LSH

using Primes
using SparseArrays
using Statistics

const DEFAULT_THRESHOLD = 0.75

struct LSHModel
    A::Vector{Int}
    B::Vector{Int}
    c::Int
    termhashes::Array{Int,2}
    dochashes::Array{Int,2}
    bandhash::Dict{Float64,Any}

    function LSHModel(nhashes::Int, dtm::SparseMatrixCSC{Int,Int})
        ndocs,nterms = size(dtm)
        c = nextprime(nterms)
        A = rand(1:c, nhashes)
        B = rand(1:c, nhashes)

        termhashes = Array{Int,2}(undef,nterms,nhashes)
        dochashes = Array{Int,2}(undef,ndocs,nhashes)

        termhash!(termhashes, A, B, c)
        apply_term_weights!(termhashes, dtm)
        dochash!(dochashes, termhashes, dtm)

        new(A, B, c, termhashes, dochashes, Dict{Float64,Any}())
    end
end

# Weighted Minhash
function termhash!(termhashes::Array{Int,2}, A::Vector{Int}, B::Vector{Int}, c::Int)
    @info("generating term hashes...")
    nterms,nhashes = size(termhashes)
    terms = 1:nterms
    for h in 1:nhashes
        a = A[h]
        b = B[h]
        termhashes[:,h] = (x->(a*x+b)%c).(terms)
    end
    nothing
end

function apply_term_weights!(termhashes::Array{Int,2}, dtm::SparseMatrixCSC{Int,Int})
    @info("applying term weights...")
    ndocs,nterms = size(dtm)
    twf = [log(1+ndocs/sum(dtm[:,t])) for t in 1:nterms]
    mw = minimum(twf)
    termweights = [ceil(Int, w/mw) for w in twf]
    termhashes .*= termweights
    nothing
end

function dochash!(dochashes::Array{Int,2}, termhashes::Array{Int,2}, dtm::SparseMatrixCSC{Int,Int})
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
similarity(sig1::Vector{Int}, sig2::Vector{Int}) = mean(.==(sig1, sig2))
threshold(nbands::Int, nrows::Int) = (1/nbands) ^ (1/nrows)

function factors(n::Int)
    allfactors = Int[]

    for i in 1:n
        if n % i == 0
            push!(allfactors, i)
        end
    end

    allfactors
end

# nbands and nrows are across hash algorithms
function best_partition(nhashes::Int, threslim::Float64)
    allfactors = factors(nhashes)
    mindiff = Inf
    bestpart = (0,0)

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

function bandhash(model::LSHModel, threshold::Float64)
    get!(model.bandhash, threshold) do
        dochashes = model.dochashes
        ndocs,nhashes = size(dochashes)
        nbands,nrows = best_partition(nhashes, threshold)

        ht = Dict{Vector{Int},Vector{Int}}()

        for d in 1:ndocs
            dh = dochashes[d,:]
            for b in 1:nbands
                brange = (1:nrows) .+ ((b-1)*nrows)
                bsig = dh[brange]
                bsigdocs = get!(ht, bsig) do
                    Int[]
                end
                append!(bsigdocs, d)
            end
        end
        return ht
    end
end

match(model::LSHModel, sig::Array{Int,2}, threshold::Float64=DEFAULT_THRESHOLD) = [match(model, sig[d,:], threshold) for d in 1:size(sig,1)]

function match(model::LSHModel, sig::Vector{Int}, threshold::Float64=DEFAULT_THRESHOLD)
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

function signature(model::LSHModel, tm::SparseVector{Int,Int})
    termhashes = model.termhashes
    nterms,nhashes = size(termhashes)
    nzinds = SparseArrays.nonzeroinds(tm)
    [minimum(termhashes[nzinds,h] .* tm[nzinds]) for h in 1:nhashes]
end

function signature(model::LSHModel, tm::SparseMatrixCSC{Int,Int})
    nterms1,nhashes = size(model.termhashes)
    ndocs,nterms2 = size(tm)
    @assert(nterms1 == nterms2)
    dochashes = Array{Int,2}(undef, ndocs, nhashes)

    for d in 1:ndocs
        dochashes[d,:] = signature(model, tm[d,:])
    end

    dochashes 
end

function selftest(model::LSHModel)
    ndocs,nhashes = size(model.dochashes)
    for d in 1:ndocs
        sig = model.dochashes[d,:]
        matches = match(model, sig)
        (d in [x.id for x in matches]) || @warn("docid did not match", docid=d)
        
        if (length(matches) > 1) && (d != first(matches).id)
            # check if there are multiple matches with same score
            topscore = first(matches).similarity
            docscore = matches[findfirst(x->x.id == d, matches)].similarity
            (docscore == topscore) || @warn("doc not best batch for its own signature", docid=d, matches=matches)
        end
   end
end

end # module LSH
