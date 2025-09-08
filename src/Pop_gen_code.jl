
using Distributions


function WrightFisher(N,p,t,μ;every=1)
    g=1
    freqs=Float64[]
    append!(freqs,p)
    while g<t
        current = freqs[g]
        current = current*(1-μ)+(1-current)*μ
        next = only(rand(Binomial(N,current),1)/N)
        append!(freqs,next)
        g += 1
    end
    ret = union(1:every:t,t)
    return(freqs[ret])
end


function WrightFisher_long(N,p)
    freqs=Float64[]
    append!(freqs,p)
    while last(freqs)*(1-last(freqs)) > 0
        next = only(rand(Binomial(N,last(freqs)),1)/N)
        append!(freqs,next)
    end
    return(freqs)
end

function fake_alignment(n,L;gaps=false)
    if gaps
        geno_matrix = rand([missing,"A","T","C","G"],n,L)
    else
        geno_matrix = rand(["A","T","C","G"],n,L)
    end
    return(geno_matrix)
end

function is_segregating(allele)
    if(length(unique(allele))>1)
        return(true)
    else
        return(false)
    end
end

function fake_alignment_biallelic(ancestral::Vector{String},n::Int64;gaps=0,adj=0)
   alignment = repeat(reshape(ancestral,1,length(ancestral)),n)
   sfs_rates = [1/i for i in 1:n]
   sfs_rates = sfs_rates ./ sum(sfs_rates)
   mutation_rates = (rand(Categorical(sfs_rates),length(ancestral)) .- 1)/(n+adj)
   mutations =  reduce(hcat,[rand(Bernoulli(mutation_rates[i]),n) for i in 1:length(ancestral)])
   gaps = rand(Bernoulli(gaps),size(alignment))
   alt_alleles = [rand(setdiff(["A","T","C","G"],[x]),1) for x in ancestral]
   for i in 1:length(ancestral)
        alignment[findall(mutations[:,i]),i] .= alt_alleles[i]
        alignment[findall(gaps[:,i]),i] .= "-"
   end
   return(alignment)
end


function only_segregating(gm)
    idx = findall([is_segregating(i) for i in eachcol(gm)])
    return(gm[:,idx])
end

function geno_mat_to_Int(gm)
    geno_mat = zeros(size(gm))
    for i in 1:size(gm)[2]
        cases = unique(gm[:,i])
        rep = Dict(cases[x]=>x for x in 1:length(cases))
        geno_mat[:,i] = [rep[x] for x in gm[:,i]] 
    end
    return(geno_mat)
end


function fis(counts::Vector{Int64})
    p = (counts[1]+0.5*counts[2])/sum(counts)
    hexp = 2*p*(1-p)
    fis = 1.0-(counts[2]/sum(counts)/hexp)
    return(fis)
end