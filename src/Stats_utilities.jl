using MultivariateStats
using Loess
using StatsBase
using Statistics

function scaled_loadings(PCA)
    scaled_loadings = loadings(PCA).*transpose(principalvars(PCA).^0.5)
    return(scaled_loadings)
end

function broken_stick_dist(r,n)
    1/n*sum([1/(n-i+1) for i in 1:r])
end

function broken_stick_pca(pca)
    prinvars = principalvars(pca)/sum(principalvars(pca))
    vals = [broken_stick_dist(i,length(prinvars)) for i in 1:length(prinvars)]
    sig_pcs = findfirst((vals .- prinvars) .< 0)-1
    return(reverse(vals), sig_pcs)
end

function scree(pca;broken_stick=false,margins=5mm)
    prinvars = principalvars(pca)/sum(principalvars(pca))
    ret_plot = bar(prinvars,xlabel="PC",ylabel="Variance explained",legend=false,margins=margins)
    if broken_stick
        vals, sig_pcs = broken_stick_pca(pca)
        scatter!(ret_plot,vals,marker_z=reduce(vcat,[ones(sig_pcs),zeros(length(prinvars)-sig_pcs)]))
    end
    return(ret_plot)
end

"""
    loess!(plot,x,y)
    adds a loess curve to a plot
"""

function loess!(p,x,y)
    model = loess(x,y)
    us = range(extrema(x)...;step=(extrema(x)[2]-extrema(x)[1])/100)
    vs = predict(model,us)
    plot!(p,us,vs,legend=false)
end

function nonunique(x::AbstractArray{T}) where T
    xs = sort(x)
    duplicatedvector = T[]
    for i=2:length(xs)
        if (isequal(xs[i],xs[i-1]) && (length(duplicatedvector)==0 || !isequal(duplicatedvector[end], xs[i])))
            push!(duplicatedvector,xs[i])
        end
    end
    duplicatedvector
end
