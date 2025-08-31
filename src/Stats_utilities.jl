
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

function pca_plot(data_frame)
    pca_model = fit(PCA,data_frame)
    

end