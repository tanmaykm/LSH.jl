using Test
using LSH
using DelimitedFiles
using SparseArrays

function test_typed(T)
    @info("testing with type $T")

    # testinput.csv prepared from document-term matrix of few SPDX files:
    # BSD-1-Clause.txt, GPL-1.0-only.txt, LGPL-2.0-only.txt, MIT.txt, Vim.txt
    dtm = convert(SparseMatrixCSC{T,T}, sparse(convert(Array{T,2}, readdlm("testinput.csv", ','))))
    model = LSHModel(20, dtm)

    @info("testing signatures...")
    signatures = signature(model, dtm)
    @test signatures == model.dochashes

    nfailures = LSH.selftest(model)
    @test nfailures == 0

    allmatches = LSH.match(model, model.dochashes)
    for idx in 1:length(allmatches)
        @test first(allmatches[idx]).id == idx
    end
end

test_typed(Int64)
test_typed(Int32)
