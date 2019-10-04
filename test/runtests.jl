using Test
using LSH
using DelimitedFiles
using SparseArrays

# testinput.csv prepared from document-term matrix of few SPDX files:
# BSD-1-Clause.txt, GPL-1.0-only.txt, LGPL-2.0-only.txt, MIT.txt, Vim.txt
dtm = sparse(convert(Array{Int64,2}, readdlm("testinput.csv", ',')))
model = LSHModel(20, dtm)
nfailures = LSH.selftest(model)
@test nfailures == 0
