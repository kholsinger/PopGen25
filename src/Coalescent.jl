include("Stats_utilities.jl")

function coalesce(N,samples)
    ancestor = [(x,rand(1:N,1)[1]) for x in samples]
    return(ancestor)
end

function coal_sim(N,samples;stop=missing)
    lineages = samples
    coalescent = [[(x,x) for x in samples]]
    while(keep_running(coalescent,lineages,stop))
        ancs = coalesce(N,lineages)
        coalescent = push!(coalescent,ancs)
        lineages = unique([x[2] for x in ancs])
    end
    return(coalescent)
end

function keep_running(coalescent,lineages,stop)
    if !ismissing(stop)
        length(coalescent) < stop
    else
        length(lineages) >1
    end
end

function disentangle(coalescent)
    for g in 2:length(coalescent)
        gen = coalescent[g]
        for sample in gen

        end
    end
    new_coal
end



function col_func_coal(samples,coals,N;def_col=:black,sel_col=:white,coal_col=:red)
    cols = fill(def_col,N)
    cols[samples] .= sel_col
    cols[coals] .= coal_col
    return(cols)
end

function modify_pop_values(coalescent,new_idx_start)
    new_coal = [[(y[1] + new_idx_start, y[2] + new_idx_start) for y in x] for x in coalescent]
    return(new_coal)
end

function fuse_pops(coal_1,coal_2)
    samples = vcat([x[2] for x in coal_1[length(coal_1)]],[x[2] for x in coal_2[length(coal_2)]])
    return(samples)
end

function coal_plot(coalescent,pop;t=length(coalescent))
    samples = [x[1] for x in coalescent[1]]
    N = length(pop)
    p = scatter(pop,fill(1,N),c=col_func_coal(samples,[],N),leg=false,ylim=(0,max(5,length(coalescent)+1));markersize=200/max(N,length(coalescent)))
    for i in 2:t
        samples = [x[2] for x in coalescent[i]]
        coals = nonunique(samples)
        scatter!(p,pop,fill(i,N),c=col_func_coal(samples,coals,N);markersize=200/max(N,length(coalescent)))
        arrow = Plots.arrow(:closed,:head,0.25/max(N,t),0.25/max(N,t))
        for x in coalescent[i]
            plot!(p,[x[1],x[2]],[i-1+0.25,i-0.25],arrow=arrow,c=:white)
        end
    end
    return(p)
end